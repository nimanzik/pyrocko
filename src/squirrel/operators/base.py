# http://pyrocko.org - GPLv3
#
# The Pyrocko Developers, 21st Century
# ---|P------/S----------~Lg----------

from __future__ import absolute_import, print_function

import re
import fnmatch
from collections import defaultdict

from ..model import QuantityType, separator
from .. import error

from pyrocko.guts import Object, String, Duration, Float


def _cglob_translate(creg):
    dd = []
    for c in creg:
        if c == '*':
            d = r'[^%s]*' % separator
        elif c == '?':
            d = r'[^%s]' % separator
        elif c == '.':
            d = separator
        else:
            d = c

        dd.append(d)
    reg = ''.join(dd)
    return reg


_compiled_patterns = {}


def compiled(pattern):
    if pattern not in _compiled_patterns:
        rpattern = re.compile(fnmatch.translate(pattern), re.I)
        _compiled_patterns[pattern] = rpattern
    else:
        rpattern = _compiled_patterns[pattern]

    return rpattern


class Filtering(Object):

    def filter(self, it):
        return list(it)


class RegexFiltering(Object):
    pattern = String.T(default=r'(.*)')

    def __init__(self, **kwargs):
        Filtering.__init__(self, **kwargs)
        self._compiled_pattern = re.compile(self.pattern)

    def filter(self, it):
        return [
            x for x in it if self._compiled_pattern.fullmatch(x)]


class Grouping(Object):

    def key(self, codes):
        return codes

    def group(self, it):
        grouped = defaultdict(list)
        for codes in it:
            k = self.key(codes)
            grouped[k].append(codes)

        for k in grouped:
            grouped[k] = tuple(grouped[k])

        return grouped


class RegexGrouping(Grouping):
    pattern = String.T(default=r'(.*)')

    def __init__(self, **kwargs):
        Grouping.__init__(self, **kwargs)
        self._compiled_pattern = re.compile(self.pattern)

    def key(self, codes):
        return self._compiled_pattern.fullmatch(codes).groups()


class ComponentGrouping(RegexGrouping):
    pattern = String.T(default=_cglob_translate('(*.*.*.*.*)?(.*)'))


class Naming(Object):
    suffix = String.T(default='')

    def get_name(self, base):
        return base + self.suffix


class RegexNaming(Naming):
    pattern = String.T(default=r'(.*)')
    replacement = String.T(default=r'\1')

    def __init__(self, **kwargs):
        Naming.__init__(self, **kwargs)
        self._compiled_pattern = re.compile(self.pattern)

    def get_name(self, base):
        return self._compiled_pattern.sub(
            self.replacement, base) + self.suffix


class ReplaceComponentNaming(RegexNaming):
    pattern = String.T(default=_cglob_translate('(*.*.*.*.*)?(.*)'))
    replacement = String.T(default=r'\1{component}\2')


class Operator(Object):

    input_codes_filtering = Filtering.T(default=Filtering.D())
    input_codes_grouping = Grouping.T(default=Grouping.D())
    output_codes_naming = Naming.T(default=Naming.D())

    def __init__(self, **kwargs):
        Object.__init__(self, **kwargs)
        self._output_names_cache = {}

    @property
    def name(self):
        return self.__class__.__name__

    def iter_mappings(self, available):
        filtered = self.input_codes_filtering.filter(available)
        grouped = self.input_codes_grouping.group(filtered)
        for k, group in grouped.items():
            yield (group, self._get_output_names_cached(group))

    def _get_output_names_cached(self, group):
        if group not in self._output_names_cache:
            self._output_names_cache[group] = self._get_output_names(group)

        return self._output_names_cache[group]

    def _get_output_names(self, group):
        return [self.output_codes_naming.get_name(codes) for codes in group]

    def get_channels(self, *args, **kwargs):
        pass

    def get_waveforms(self, squirrel, in_codes, out_codes, **kwargs):
        assert len(in_codes) == 1 and len(out_codes) == 1
        in_codes_tup = in_codes[0].split(separator)
        trs = squirrel.get_waveforms(codes=in_codes_tup, **kwargs)
        agn, net, sta, loc, cha, ext = out_codes[0].split(separator)
        for tr in trs:
            tr.set_codes(
                agency=agn,
                network=net,
                station=sta,
                location=loc,
                channel=cha,
                extra=ext)

        return trs


class Parameters(Object):
    pass


class RestitutionParameters(Parameters):
    frequency_min = Float.T()
    frequency_max = Float.T()
    frequency_taper_factor = Float.T(default=1.5)
    time_taper_factor = Float.T(default=2.0)


class Restitution(Operator):
    output_codes_naming = Naming(suffix='R{quantity}')
    quantity = QuantityType.T(default='displacement')

    @property
    def name(self):
        return 'Restitution(%s)' % self.quantity[0]

    def _get_output_names(self, group):
        return [
            self.output_codes_naming.get_name(codes).format(
                quantity=self.quantity[0])
            for codes in group]

    def get_tpad(self, params):
        return params.time_taper_factor / params.frequency_min

    def get_waveforms(
            self, squirrel, in_codes, out_codes, params, tmin, tmax, **kwargs):

        assert len(in_codes) == 1 and len(out_codes) == 1
        in_codes_tup = in_codes[0].split(separator)

        tpad = self.get_tpad(params)

        tmin_raw = tmin - tpad
        tmax_raw = tmax + tpad

        trs = squirrel.get_waveforms(
            codes=in_codes_tup, tmin=tmin_raw, tmax=tmax_raw, **kwargs)

        try:
            resp = squirrel.get_response(
                codes=in_codes_tup,
                tmin=tmin_raw,
                tmax=tmax_raw).get_effective()

        except error.NotAvailable:
            return []

        freqlimits = (
            params.frequency_min / params.frequency_taper_factor,
            params.frequency_min,
            params.frequency_max,
            params.frequency_max * params.frequency_taper_factor)

        agn, net, sta, loc, cha, ext = out_codes[0].split(separator)
        trs_rest = []
        for tr in trs:
            tr_rest = tr.transfer(
                tfade=tpad,
                freqlimits=freqlimits,
                transfer_function=resp,
                invert=True)

            tr_rest.set_codes(
                agency=agn,
                network=net,
                station=sta,
                location=loc,
                channel=cha,
                extra=ext)

            trs_rest.append(tr_rest)

        return trs_rest


class Shift(Operator):
    output_codes_naming = Naming(suffix='S')
    delay = Duration.T()


class Transform(Operator):
    input_codes_grouping = Grouping.T(default=ComponentGrouping.D())
    output_codes_naming = ReplaceComponentNaming(suffix='T{system}')

    def _get_output_names(self, group):
        return [
            self.output_codes_naming.get_name(group[0]).format(
                component=c, system=self.components.lower())
            for c in self.components]


class ToENZ(Transform):
    components = 'ENZ'


class ToRTZ(Transform):
    components = 'RTZ'
    backazimuth = Float.T()


class ToLTQ(Transform):
    components = 'LTQ'


class Composition(Operator):
    g = Operator.T()
    f = Operator.T()

    def __init__(self, g, f, **kwargs):
        Operator.__init__(self, g=g, f=f, **kwargs)

    @property
    def name(self):
        return '(%s ○ %s)' % (self.g.name, self.f.name)

    def iter_mappings(self, available):
        g_available = []
        deps = {}
        for f_in_codes, f_out_codes in self.f.iter_mappings(available):
            g_available.extend(f_out_codes)
            for codes in f_out_codes:
                deps[codes] = f_in_codes

        for g_in_codes, g_out_codes in self.g.iter_mappings(g_available):
            gf_in_codes = []
            for codes in g_in_codes:
                gf_in_codes.extend(deps[codes])

            yield gf_in_codes, g_out_codes


__all__ = [
    'Operator',
    'Restitution',
    'Shift',
    'ToENZ',
    'ToRTZ',
    'ToLTQ',
    'Composition']
