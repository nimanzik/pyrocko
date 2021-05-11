# http://pyrocko.org - GPLv3
#
# The Pyrocko Developers, 21st Century
# ---|P------/S----------~Lg----------

from __future__ import absolute_import, print_function

import time

from .. import common
from pyrocko.squirrel.operators import Restitution, ToENZ, ToRTZ, Shift
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
    squirrel.add_operator(Restitution(quantity='velocity'))
    squirrel.add_operator(ToENZ())
    squirrel.add_operator(ToRTZ())
    squirrel.add_operator(Shift(delay=100.))

    def scodes(codes):
        return ','.join(c.replace(separator, '.') for c in codes)

    t0 = time.time()
    ops = squirrel.get_operator_mappings()
    t1 = time.time()
    ops2 = squirrel.get_operator_mappings()
    t2 = time.time()

    for operator, in_codes, out_codes in squirrel.get_operator_mappings():
        print('%s => %s => %s' % (
            scodes(in_codes), operator.name, scodes(out_codes)))

    print(t1-t0, t2-t1)
