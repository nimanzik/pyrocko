# http://pyrocko.org - GPLv3
#
# The Pyrocko Developers, 21st Century
# ---|P------/S----------~Lg----------

from __future__ import absolute_import, print_function

import re
import fnmatch
from collections import defaultdict

from ..model import QuantityType

from pyrocko.guts import Object, String, Duration

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
    pattern = String.T(default=r'(.*)\Z')

    def __init__(self, **kwargs):
        Filtering.__init__(self, **kwargs)
        self._compiled_pattern = re.compile(self.pattern)

    def filter(self, it):
        return [
            x for x in it if self._compiled_pattern.match(x)]


class Grouping(Object):

    def key(self, codes):
        return codes

    def group(self, it):
        grouped = defaultdict(list)
        for codes in it:
            k = self.key(codes)
            grouped[k].append(codes)

        return grouped


class RegexGrouping(Grouping):
    pattern = String.T(default=r'(.*)\Z')

    def __init__(self, **kwargs):
        Grouping.__init__(self, **kwargs)
        self._compiled_pattern = re.compile(self.pattern)

    def key(self, codes):
        return self._compiled_pattern.match(codes).groups()


class ComponentGrouping(RegexGrouping):
    pattern = String.T(default=r'(.*\t.*\t.*\t.*\t.*).(\t.*)\Z')


class Naming(Object):
    suffix = String.T(optional=True)

    def get_name(self, base, **params):
        return base + self.suffix.format(**params) if self.suffix else ''


class RegexNaming(Naming):
    pattern = String.T(default=r'(.*)\Z')
    replacement = String.T(default=r'\1')

    def __init__(self, **kwargs):
        Naming.__init__(self, **kwargs)
        self._compiled_pattern = re.compile(self.pattern)

    def get_name(self, base, **params):
        return self._compiled_pattern.sub(
            self.replacement.format(**params), base) \
                + self.suffix.format(**params) if self.suffix else ''


class ReplaceComponentNaming(RegexNaming):
    pattern = String.T(default=r'(.*\t.*\t.*\t.*\t.*).(\t.*)\Z')
    replacement = String.T(default=r'\1{component}\2')


class Operator(Object):

    input_codes_filtering = Filtering.T(default=Filtering.D())
    input_codes_grouping = Grouping.T(default=Grouping.D())
    output_codes_naming = Naming.T()

    @property
    def name(self):
        return self.__class__.__name__

    def iter_mappings(self, available):
        filtered = self.input_codes_filtering.filter(available)
        grouped = self.input_codes_grouping.group(filtered)
        for k, group in grouped.items():
            yield (group, self._get_output_names(group))

    def _get_output_names(self, group):
        return [self.output_codes_naming.get_name(codes) for codes in group]

    def get_channels(self, *args, **kwargs):
        pass

    def get_waveforms(self, *args, **kwargs):
        pass


class Restitution(Operator):
    quantity = QuantityType.T(default='displacement')
    output_codes_naming = Naming(suffix='R{quantity[0]}')

    def _get_output_names(self, group):
        return [
            self.output_codes_naming.get_name(codes, quantity=self.quantity)
            for codes in group]


class Shift(Operator):
    output_codes_naming = Naming(suffix='S')
    delay = Duration.T()


class Transform(Operator):
    input_codes_grouping = Grouping.T(default=ComponentGrouping.D())
    output_codes_naming = ReplaceComponentNaming(suffix='T{system}')

    def _get_output_names(self, group):
        return [
            self.output_codes_naming.get_name(
                group[0], component=c, system=self.components.lower())
            for c in self.components]


class ToENZ(Transform):
    components = 'ENZ'


class ToRTZ(Transform):
    components = 'RTZ'


class ToLTQ(Transform):
    components = 'LTQ'


__all__ = [
    'Operator',
    'Restitution',
    'Shift',
    'ToENZ',
    'ToRTZ',
    'ToLTQ']
