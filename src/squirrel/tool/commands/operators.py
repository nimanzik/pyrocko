# http://pyrocko.org - GPLv3
#
# The Pyrocko Developers, 21st Century
# ---|P------/S----------~Lg----------

from __future__ import absolute_import, print_function

from .. import common
from pyrocko.squirrel.operators import Restitution, ToZEN


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
    squirrel.add_operator(Restitution(target='displacement'))
    squirrel.add_operator(ToZEN())

    for mapping in squirrel.get_operator_mappings():
        print(mapping)
