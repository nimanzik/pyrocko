# http://pyrocko.org - GPLv3
#
# The Pyrocko Developers, 21st Century
# ---|P------/S----------~Lg----------
from __future__ import division, absolute_import

from struct import unpack
import os
import re
import math
import logging

from pyrocko import trace
from pyrocko.util import reuse, ensuredirs
from .io_common import FileLoadError, FileSaveError

logger = logging.getLogger('pyrocko.io.mseed')

MSEED_HEADER_BYTES = 64
VALID_RECORD_LENGTHS = tuple(2**exp for exp in range(8, 20))


class CodeTooLong(FileSaveError):
    pass


def iload(filename, load_data=True, segment=-1, segment_nrecords=512):
    from pyrocko import mseed_ext

    have_zero_rate_traces = False
    try:
        traces = []
        segments = []
        for tr in mseed_ext.get_traces(
                filename, load_data, segment, segment_nrecords):
            network, station, location, channel = tr[1:5]
            tmin = float(tr[5])/float(mseed_ext.HPTMODULUS)
            tmax = float(tr[6])/float(mseed_ext.HPTMODULUS)
            try:
                deltat = reuse(1.0/float(tr[7]))
            except ZeroDivisionError:
                have_zero_rate_traces = True
                continue

            ydata = tr[8]

            traces.append(trace.Trace(
                network, station, location, channel, tmin, tmax,
                deltat, ydata))

            segments.append(tr[9])

        for tr in traces:
            yield tr

    except (OSError, mseed_ext.MSeedError) as e:
        raise FileLoadError(str(e)+' (file: %s)' % filename)

    if have_zero_rate_traces:
        logger.warning(
            'Ignoring traces with sampling rate of zero in file %s '
            '(maybe LOG traces)' % filename)


def as_tuple(tr, dataquality='D'):
    from pyrocko import mseed_ext
    itmin = int(round(tr.tmin*mseed_ext.HPTMODULUS))
    itmax = int(round(tr.tmax*mseed_ext.HPTMODULUS))
    srate = 1.0/tr.deltat
    return (tr.network, tr.station, tr.location, tr.channel,
            itmin, itmax, srate, dataquality, tr.get_ydata())


def save(traces, filename_template, additional={}, overwrite=True,
         dataquality='D', record_length=4096, append=False):
    from pyrocko import mseed_ext

    assert record_length in VALID_RECORD_LENGTHS
    assert dataquality in ('D', 'E', 'C', 'O', 'T', 'L'), 'invalid dataquality'
    overwrite = True if append else overwrite

    fn_tr = {}
    for tr in traces:
        for code, maxlen, val in zip(
                ['network', 'station', 'location', 'channel'],
                [2, 5, 2, 3],
                tr.nslc_id):

            if len(val) > maxlen:
                raise CodeTooLong(
                    '%s code too long to be stored in MSeed file: %s' %
                    (code, val))

        fn = tr.fill_template(filename_template, **additional)
        if not overwrite and os.path.exists(fn):
            raise FileSaveError('File exists: %s' % fn)

        if fn not in fn_tr:
            fn_tr[fn] = []

        fn_tr[fn].append(tr)

    for fn, traces_thisfile in fn_tr.items():
        trtups = []
        traces_thisfile.sort(key=lambda a: a.full_id)
        for tr in traces_thisfile:
            trtups.append(as_tuple(tr, dataquality))

        ensuredirs(fn)
        try:
            mseed_ext.store_traces(trtups, fn, record_length, append)
        except mseed_ext.MSeedError as e:
            raise FileSaveError(
                str(e) + ' (while storing traces to file \'%s\')' % fn)

    return list(fn_tr.keys())


tcs = {}


def get_bytes(traces, dataquality='D', record_length=4096):
    from pyrocko import mseed_ext

    assert record_length in VALID_RECORD_LENGTHS
    assert dataquality in ('D', 'E', 'C', 'O', 'T', 'L'), 'invalid dataquality'

    nbytes_approx = 0
    rl = record_length
    trtups = []
    for tr in traces:
        for code, maxlen, val in zip(
                ['network', 'station', 'location', 'channel'],
                [2, 5, 2, 3],
                tr.nslc_id):

            if len(val) > maxlen:
                raise CodeTooLong(
                    '%s code too long to be stored in MSeed file: %s' %
                    (code, val))

        nbytes_approx += math.ceil(
            tr.ydata.nbytes / (rl-MSEED_HEADER_BYTES)) * rl
        trtups.append(as_tuple(tr, dataquality))

    return mseed_ext.mseed_bytes(trtups, nbytes_approx, record_length)


def detect(first512):

    if len(first512) < 256:
        return False

    rec = first512

    try:
        sequence_number = int(rec[:6])
    except Exception:
        return False
    if sequence_number < 0:
        return False

    type_code = rec[6:7]
    if type_code in b'DRQM':
        bads = []
        for sex in '<>':
            bad = False
            fmt = sex + '6s1s1s5s2s3s2s10sH2h4Bl2H'
            vals = unpack(fmt, rec[:48])
            fmt_btime = sex + 'HHBBBBH'
            tvals = unpack(fmt_btime, vals[7])
            if tvals[1] < 1 or tvals[1] > 367 or tvals[2] > 23 or \
                    tvals[3] > 59 or tvals[4] > 60 or tvals[6] > 9999:
                bad = True

            bads.append(bad)

        if all(bads):
            return False

    else:
        if type_code not in b'VAST':
            return False

        continuation_code = rec[7:8]
        if continuation_code != b' ':
            return False

        blockette_type = rec[8:8+3].decode()
        if not re.match(r'^\d\d\d$', blockette_type):
            return False

        try:
            blockette_length = int(rec[11:11+4])
        except Exception:
            return False

        if blockette_length < 7:
            return False

    return True
