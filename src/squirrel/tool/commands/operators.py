# http://pyrocko.org - GPLv3
#
# The Pyrocko Developers, 21st Century
# ---|P------/S----------~Lg----------

from __future__ import absolute_import, print_function

from pyrocko import util
from .. import common
from pyrocko.squirrel.operators import Restitution, ToENZ, ToRTZ, Shift, \
    Operator, Composition
from pyrocko.squirrel.model import separator


guts_prefix = 'squirrel'


def setup(subparsers):
    p = common.add_parser(
        subparsers, 'operators',
        help='Print available operator mappings.')

    common.add_selection_arguments(p)
    common.add_query_arguments(p)
    return p


def call(parser, args):
    d = common.squirrel_query_from_arguments(args)
    d
    squirrel = common.squirrel_from_selection_arguments(args)
    # squirrel.add_operator(
    #      Composition(ToRTZ(), Restitution(quantity='velocity')))
    #squirrel.add_operator(Restitution(quantity='velocity'))

    def scodes(codes):
        css = list(zip(*(c.split(separator) for c in codes)))
        if sum(not all(c == cs[0] for c in cs) for cs in css) == 1:
            return '.'.join(
                cs[0] if all(c == cs[0] for c in cs) else '(%s)' % ','.join(cs)
                for cs in css)
        else:
            return ', '.join(c.replace(separator, '.') for c in codes)

    tmin = util.str_to_time('2020-03-10 00:00:00')
    tmax = tmin + 60.

    for operator, in_codes, out_codes in squirrel.get_operator_mappings():
        print('%s <- %s <- %s' % (
            scodes(out_codes), operator.name, scodes(in_codes)))

        # trs = operator.get_waveforms(
        #     squirrel, in_codes, out_codes, tmin=tmin, tmax=tmax)
        #
        # for tr in trs:
        #     print(tr)
