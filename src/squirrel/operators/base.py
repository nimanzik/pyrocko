# http://pyrocko.org - GPLv3
#
# The Pyrocko Developers, 21st Century
# ---|P------/S----------~Lg----------

from __future__ import absolute_import, print_function

import re
import fnmatch
from collections import defaultdict

from ..model import QuantityType, separator

from pyrocko.guts import Object, Tuple, String, Duration

_compiled_patterns = {}


def compiled(pattern):
    if pattern not in _compiled_patterns:
        rpattern = re.compile(fnmatch.translate(pattern), re.I)
        _compiled_patterns[pattern] = rpattern
    else:
        rpattern = _compiled_patterns[pattern]

    return rpattern


class Grouping(Object):
    pattern = re.compile(r'(.*)\Z')

    def key(self, codes):
        return self.pattern.match(codes).groups()


class ComponentGrouping(Grouping):
    pattern = re.compile(r'(.*\t.*\t.*\t.*\t.*).(\t.*)\Z')


class Operator(Object):

    input_codes_pattern = Tuple.T(
        6, String.T(), default=('*', '*', '*', '*', '*', '*'))
    input_codes_grouping = Grouping.T(default=Grouping.D())

    def iter_mappings(self, squirrel):
        available = [
            separator.join(codes)
            for codes in squirrel.get_codes(kind='waveform')]

        c_in_pat = compiled(self.input_codes_pattern)

        filtered = [codes for codes in available if c_in_pat.match(codes)]

        grouped = defaultdict(list)
        for codes in filtered:
            k = self.input_codes_grouping.key(codes)
            grouped[k].append(codes)

        for k, group in grouped.items():
            yield(group)

    def get_channels(self, *args, **kwargs):
        pass

    def get_waveforms(self, *args, **kwargs):
        pass



class Restitution(Operator):
    target = QuantityType.T()


class Shift(Operator):
    delay = Duration.T()


class ToENZ(Operator):
    input_codes_grouping = Grouping.T(default=ComponentGrouping.D())


class ToRTZ(Operator):
    input_codes_grouping = Grouping.T(default=ComponentGrouping.D())


class ToLTQ(Operator):
    input_codes_grouping = Grouping.T(default=ComponentGrouping.D())


__all__ = [
    'Operator',
    'Restitution',
    'Shift',
    'ToENZ',
    'ToRTZ',
    'ToLTQ']
