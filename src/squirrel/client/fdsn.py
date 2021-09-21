# http://pyrocko.org - GPLv3
#
# The Pyrocko Developers, 21st Century
# ---|P------/S----------~Lg----------

from __future__ import absolute_import, print_function

import time
import os
import copy
import logging
import tempfile
from collections import defaultdict
try:
    import cPickle as pickle
except ImportError:
    import pickle
import os.path as op
from .base import Source, Constraint
from ..model import make_waveform_promise_nut, ehash, InvalidWaveform, \
    order_summary, WaveformOrder
from pyrocko.client import fdsn

from pyrocko import util, trace, io
from pyrocko.io_common import FileLoadError
from pyrocko.io import stationxml

from pyrocko.guts import Object, String, Timestamp, List, Tuple, Int, Dict, \
    Duration, Bool

guts_prefix = 'squirrel'

fdsn.g_timeout = 60.

logger = logging.getLogger('pyrocko.squirrel.client.fdsn')

sites_not_supporting = {
    'startbefore': ['geonet'],
    'includerestricted': ['geonet']}


def plural_s(x):
    if not isinstance(x, int):
        x = len(x)

    return 's' if x != 1 else ''


def diff(fn_a, fn_b):
    try:
        if os.stat(fn_a).st_size != os.stat(fn_b).st_size:
            return True

    except OSError:
        return True

    with open(fn_a, 'rb') as fa:
        with open(fn_b, 'rb') as fb:
            while True:
                a = fa.read(1024)
                b = fb.read(1024)
                if a != b:
                    return True

                if len(a) == 0 or len(b) == 0:
                    return False


class Archive(Object):

    def add(self):
        raise NotImplementedError()


class MSeedArchive(Archive):
    template = String.T(default=op.join(
        '%(tmin_year)s',
        '%(tmin_month)s',
        '%(tmin_day)s',
        'trace_%(network)s_%(station)s_%(location)s_%(channel)s'
        + '_%(tmin_us)s_%(tmax_us)s.mseed'))

    def __init__(self, **kwargs):
        Archive.__init__(self, **kwargs)
        self._base_path = None

    def set_base_path(self, path):
        self._base_path = path

    def add(self, trs):
        path = op.join(self._base_path, self.template)
        return io.save(trs, path, overwrite=True)


def combine_selections(selection):
    out = []
    last = None
    for this in selection:
        if last and this[:4] == last[:4] and this[4] == last[5]:
            last = last[:5] + (this[5],)
        else:
            if last:
                out.append(last)

            last = this

    if last:
        out.append(last)

    return out


def orders_sort_key(order):
    return (order.codes, order.tmin)


def orders_to_selection(orders):
    selection = []
    for order in sorted(orders, key=orders_sort_key):
        selection.append(
            order.codes[1:5] + (order.tmin, order.tmax))

    return combine_selections(selection)


class ErrorEntry(Object):
    time = Timestamp.T()
    order = WaveformOrder.T()
    kind = String.T()
    details = String.T(optional=True)


class ErrorAggregate(Object):
    site = String.T()
    kind = String.T()
    details = String.T()
    entries = List.T(ErrorEntry.T())
    codes_list = List.T(Tuple.T(None, String.T()))
    time_spans = List.T(Tuple.T(2, Timestamp.T()))

    def __str__(self):
        codes = ['.'.join(x) for x in self.codes_list]
        scodes = '\n' + util.ewrap(codes, indent='    ') if codes else '<none>'
        tss = self.time_spans
        sspans = '\n' + util.ewrap(('%s - %s' % (
            util.time_to_str(ts[0]), util.time_to_str(ts[1])) for ts in tss),
            indent='   ')

        return ('FDSN "%s": download error summary for "%s" (%i)\n%s  '
                'Codes:%s\n  Time spans:%s') % (
            self.site,
            self.kind,
            len(self.entries),
            '  Details: %s\n' % self.details if self.details else '',
            scodes,
            sspans)


class ErrorLog(Object):
    site = String.T()
    entries = List.T(ErrorEntry.T())
    checkpoints = List.T(Int.T())

    def append_checkpoint(self):
        self.checkpoints.append(len(self.entries))

    def append(self, time, order, kind, details=''):
        entry = ErrorEntry(time=time, order=order, kind=kind, details=details)
        self.entries.append(entry)

    def iter_aggregates(self):
        by_kind_details = defaultdict(list)
        for entry in self.entries:
            by_kind_details[entry.kind, entry.details].append(entry)

        kind_details = sorted(by_kind_details.keys())

        for kind, details in kind_details:
            entries = by_kind_details[kind, details]
            codes_list = sorted(set(entry.order.codes for entry in entries))
            selection = orders_to_selection(entry.order for entry in entries)
            time_spans = sorted(set(row[-2:] for row in selection))
            yield ErrorAggregate(
                site=self.site,
                kind=kind,
                details=details,
                entries=entries,
                codes_list=codes_list,
                time_spans=time_spans)

    def summarize_recent(self):
        ioff = self.checkpoints[-1] if self.checkpoints else 0
        recent = self.entries[ioff:]
        kinds = sorted(set(entry.kind for entry in recent))
        if recent:
            return '%i error%s (%s)' % (
                len(recent), plural_s(recent), '; '.join(kinds))
        else:
            return ''


class FDSNSource(Source):

    '''
    Squirrel data-source to transparently get data from FDSN web services.

    Attaching an :py:class:`FDSNSource` object to a :py:class:`Squirrel` allows
    the latter to download station and waveform data from an FDSN web service
    should the data not already happen to be available locally.
    '''

    site = String.T(
        help='FDSN site url or alias name (see '
             ':py:mod:`pyrocko.client.fdsn`).')

    query_args = Dict.T(
        String.T(), String.T(),
        optional=True,
        help='Common query arguments, which are appended to all queries.')

    expires = Duration.T(
        optional=True,
        help='Expiration time [s]. Information older than this will be '
             'refreshed. This only applies to station-metadata. Waveforms do '
             'not expire. If set to ``None`` neither type of data  expires.')

    cache_path = String.T(
        optional=True,
        help='Directory path where any downloaded waveforms and station '
             'meta-data are to be kept. By default the Squirrel '
             'environment\'s cache directory is used.')

    shared_waveforms = Bool.T(
        default=True,
        help='If ``True``, waveforms are shared with other FDSN sources in '
             'the same Squirrel environment. If ``False``, they are kept '
             'separate.')

    user_credentials = Tuple.T(
        2, String.T(),
        optional=True,
        help='User and password for FDSN servers requiring password '
             'authentication')

    auth_token = String.T(
        optional=True,
        help='Authentication token to be presented to the FDSN server.')

    auth_token_path = String.T(
        optional=True,
        help='Path to file containing the authentication token to be '
             'presented to the FDSN server.')

    def __init__(self, site, query_args=None, **kwargs):
        Source.__init__(self, site=site, query_args=query_args, **kwargs)

        self._constraint = None
        self._hash = self.make_hash()
        self._source_id = 'client:fdsn:%s' % self._hash
        self._error_infos = []

    def make_hash(self):
        s = self.site
        s += 'notoken' \
            if (self.auth_token is None and self.auth_token_path is None) \
            else 'token'

        if self.user_credentials is not None:
            s += self.user_credentials[0]
        else:
            s += 'nocred'

        if self.query_args is not None:
            s += ','.join(
                '%s:%s' % (k, self.query_args[k])
                for k in sorted(self.query_args.keys()))
        else:
            s += 'noqueryargs'

        return ehash(s)

    def get_hash(self):
        return self._hash

    def get_auth_token(self):
        if self.auth_token:
            return self.auth_token

        elif self.auth_token_path is not None:
            try:
                with open(self.auth_token_path, 'rb') as f:
                    return f.read().decode('ascii')

            except OSError as e:
                raise FileLoadError(
                    'Cannot load auth token file (%s): %s'
                    % (str(e), self.auth_token_path))

        else:
            raise Exception(
                'FDSNSource: auth_token and auth_token_path are mutually '
                'exclusive.')

    def setup(self, squirrel, check=True, progress_viewer='terminal'):
        self._cache_path = op.join(
            self.cache_path or squirrel._cache_path, 'fdsn')

        util.ensuredir(self._cache_path)
        self._load_constraint()
        self._archive = MSeedArchive()
        waveforms_path = self._get_waveforms_path()
        util.ensuredir(waveforms_path)
        self._archive.set_base_path(waveforms_path)

        squirrel.add(
            self._get_waveforms_path(),
            check=check,
            progress_viewer=progress_viewer)

    def _get_constraint_path(self):
        return op.join(self._cache_path, self._hash, 'constraint.pickle')

    def _get_channels_path(self):
        return op.join(self._cache_path, self._hash, 'channels.stationxml')

    def _get_waveforms_path(self):
        if self.shared_waveforms:
            return op.join(self._cache_path, 'waveforms')
        else:
            return op.join(self._cache_path, self._hash, 'waveforms')

    def _log_meta(self, message, target=logger.info):
        log_prefix = 'FDSN "%s" metadata:' % self.site
        target(' '.join((log_prefix, message)))

    def _log_info_data(self, *args):
        log_prefix = 'FDSN "%s" waveforms:' % self.site
        logger.info(' '.join((log_prefix,) + args))

    def _str_expires(self, now):
        t = self._get_expiration_time()
        if t is None:
            return 'expires: never'
        else:
            expire = 'expires' if t > now else 'expired'
            return '%s: %s' % (
                expire,
                util.time_to_str(t, format='%Y-%m-%d %H:%M:%S'))

    def update_channel_inventory(self, squirrel, constraint=None):
        if constraint is None:
            constraint = Constraint()

        expiration_time = self._get_expiration_time()
        now = time.time()

        log_target = logger.info
        if self._constraint and self._constraint.contains(constraint) \
                and (expiration_time is None or now < expiration_time):

            s_case = 'using cached'

        else:
            if self._constraint:
                constraint_temp = copy.deepcopy(self._constraint)
                constraint_temp.expand(constraint)
                constraint = constraint_temp

            try:
                channel_sx = self._do_channel_query(constraint)

                channel_sx.created = None  # timestamp would ruin diff

                fn = self._get_channels_path()
                util.ensuredirs(fn)
                fn_temp = fn + '.%i.temp' % os.getpid()
                channel_sx.dump_xml(filename=fn_temp)

                if op.exists(fn):
                    if diff(fn, fn_temp):
                        os.rename(fn_temp, fn)
                        s_case = 'updated'
                    else:
                        os.unlink(fn_temp)
                        squirrel.get_database().silent_touch(fn)
                        s_case = 'upstream unchanged'

                else:
                    os.rename(fn_temp, fn)
                    s_case = 'new'

                self._constraint = constraint
                self._dump_constraint()

            except OSError as e:
                s_case = 'update failed (%s)' % str(e)
                log_target = logger.error

        self._log_meta(
            '%s (%s)' % (s_case, self._str_expires(now)),
            target=log_target)

        fn = self._get_channels_path()
        if os.path.exists(fn):
            squirrel.add(fn)

    def _do_channel_query(self, constraint):
        extra_args = {}

        if self.site in sites_not_supporting['startbefore']:
            if constraint.tmin is not None:
                extra_args['starttime'] = constraint.tmin
            if constraint.tmax is not None:
                extra_args['endtime'] = constraint.tmax

        else:
            if constraint.tmin is not None:
                extra_args['endafter'] = constraint.tmin
            if constraint.tmax is not None:
                extra_args['startbefore'] = constraint.tmax

        if self.site not in sites_not_supporting['includerestricted']:
            extra_args.update(
                includerestricted=(
                    self.user_credentials is not None
                    or self.auth_token is not None))

        if self.query_args is not None:
            extra_args.update(self.query_args)

        self._log_meta('querying...')

        try:
            channel_sx = fdsn.station(
                site=self.site,
                format='text',
                level='channel',
                **extra_args)
            return channel_sx

        except fdsn.EmptyResult:
            return stationxml.FDSNStationXML(source='dummy-emtpy-result')

    def _load_constraint(self):
        fn = self._get_constraint_path()
        if op.exists(fn):
            with open(fn, 'rb') as f:
                self._constraint = pickle.load(f)
        else:
            self._constraint = None

    def _dump_constraint(self):
        with open(self._get_constraint_path(), 'wb') as f:
            pickle.dump(self._constraint, f, protocol=2)

    def _get_expiration_time(self):
        if self.expires is None:
            return None

        try:
            path = self._get_channels_path()
            t = os.stat(path)[8]
            return t + self.expires

        except OSError:
            return 0.0

    def update_waveform_promises(self, squirrel, constraint):
        from pyrocko.squirrel import Squirrel

        # get meta information of stuff available through this source
        sub_squirrel = Squirrel(database=squirrel.get_database())
        fn = self._get_channels_path()
        if os.path.exists(fn):
            sub_squirrel.add([fn], check=False)

        nuts = sub_squirrel.iter_nuts(
            'channel', constraint.tmin, constraint.tmax)

        path = self._source_id
        squirrel.add_virtual(
            (make_waveform_promise_nut(
                file_path=path,
                **nut.waveform_promise_kwargs) for nut in nuts),
            virtual_paths=[path])

    def _get_user_credentials(self):
        d = {}
        if self.user_credentials is not None:
            d['user'], d['passwd'] = self.user_credentials

        if self.auth_token is not None:
            d['token'] = self.get_auth_token()

        return d

    def download_waveforms(
            self, squirrel, orders, success, error_permanent, error_temporary):

        elog = ErrorLog(site=self.site)
        orders.sort(key=orders_sort_key)
        neach = 20
        i = 0
        while i < len(orders):
            orders_now = orders[i:i+neach]
            selection_now = orders_to_selection(orders_now)

            nsuccess = 0
            elog.append_checkpoint()
            self._log_info_data(
                'downloading, %s' % order_summary(orders_now))

            with tempfile.NamedTemporaryFile() as f:
                try:
                    data = fdsn.dataselect(
                        site=self.site, selection=selection_now,
                        **self._get_user_credentials())

                    now = time.time()

                    while True:
                        buf = data.read(1024)
                        if not buf:
                            break
                        f.write(buf)

                    f.flush()

                    trs = io.load(f.name)

                    by_nslc = defaultdict(list)
                    for tr in trs:
                        by_nslc[tr.nslc_id].append(tr)

                    for order in orders_now:
                        trs_order = []
                        err_this = None
                        for tr in by_nslc[order.codes[1:5]]:
                            try:
                                order.validate(tr)
                                trs_order.append(tr.chop(
                                    order.tmin, order.tmax, inplace=False))

                            except trace.NoData:
                                err_this = (
                                    'empty result', 'empty sub-interval')

                            except InvalidWaveform as e:
                                err_this = ('invalid waveform', str(e))

                        if len(trs_order) == 0:
                            if err_this is None:
                                err_this = ('empty result', '')

                            elog.append(now, order, *err_this)
                            error_permanent(order)
                        else:
                            if len(trs_order) != 1:
                                if err_this:
                                    elog.append(
                                        now, order,
                                        'partial result, %s' % err_this[0],
                                        err_this[1])
                                else:
                                    elog.append(now, order, 'partial result')

                            paths = self._archive.add(trs_order)
                            squirrel.add(paths)
                            nsuccess += 1
                            success(order)

                except fdsn.EmptyResult:
                    now = time.time()
                    for order in orders_now:
                        elog.append(now, order, 'empty result')
                        error_permanent(order)

                except util.HTTPError as e:
                    now = time.time()
                    for order in orders_now:
                        elog.append(now, order, 'http error', str(e))
                        error_temporary(order)

            emessage = elog.summarize_recent()
            self._log_info_data(
                '%i download%s successful' % (nsuccess, plural_s(nsuccess))
                + (', %s' % emessage if emessage else ''))

            i += neach

        for agg in elog.iter_aggregates():
            logger.warn(str(agg))


__all__ = [
    'FDSNSource',
]
