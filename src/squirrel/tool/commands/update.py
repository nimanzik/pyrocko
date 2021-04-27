# http://pyrocko.org - GPLv3
#
# The Pyrocko Developers, 21st Century
# ---|P------/S----------~Lg----------

from __future__ import absolute_import, print_function

from .. import common
from pyrocko.squirrel.error import SquirrelError


def setup(subparsers):
    p = common.add_parser(
        subparsers, 'update',
        help='Update remote sources inventories.')

    common.add_selection_arguments(p)
    common.add_query_arguments(p)

    p.add_argument(
        '--promises',
        action='store_true',
        dest='promises',
        default=False,
        help='Update waveform promises.')

    return p


def call(parser, args):
    d = common.squirrel_query_from_arguments(args)
    squirrel = common.squirrel_from_selection_arguments(args)

    if 'tmin' not in d or 'tmax' not in d:
        raise SquirrelError('Time span required.')

    squirrel.update(tmin=d['tmin'], tmax=d['tmax'])
    if args.promises:
        squirrel.update_waveform_promises()

    print(squirrel)
