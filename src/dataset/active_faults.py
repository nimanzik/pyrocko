# http://pyrocko.org - GPLv3
#
# The Pyrocko Developers, 21st Century
# ---|P------/S----------~Lg----------
from __future__ import absolute_import

import csv
import logging
import numpy as num
from os import path as op
from collections import OrderedDict
import json

from pyrocko import config, util


logger = logging.getLogger('ActiveFaults')


class Fault(object):
    __fields__ = OrderedDict()
    __slots__ = list(__fields__.keys())

    def __init__(self, *args):
        for f in args:
            nodes = len(f['geometry']['coordinates'])
            for (attr, attr_type) in self.__fields__.items():
                values = []
                if attr == 'average_dip':
                    try:
                        values = float(f['properties']['average_dip'][1:3])
                    except:
                        values = -999
                if attr == 'average_dip':
                    try:
                        values = float(f['properties']['average_dip'][1:3])
                    except:
                        values = -999
                elif attr == 'average_rake':
                    try:
                        values = float(f['properties']['average_rake'][1:3])
                    except:
                        values = -999
                elif attr == 'lat':
                    for i in range(0, nodes):
                        values.append(f['geometry']['coordinates'][i][1])
                elif attr == 'lon':
                    for i in range(0, nodes):
                        values.append(f['geometry']['coordinates'][i][0])
                else:
                        values = 0.
                try:
                    setattr(self, attr, attr_type(values))
                except ValueError as e:
                    print(list(zip(self.__fields__.keys(), args)))
                    raise e

    def __str__(self):
        d = {attr: getattr(self, attr) for attr in self.__fields__.keys()}
        return '\n'.join(['%s: %s' % (attr, val) for attr, val in d.items()])


class ActiveFault(Fault):
    __fields__ = OrderedDict([

        ('lat', list),
        ('lon', list),
        ('average_dip', float),
        ('average_rake', float),
        ('lower_seis_depth', float),
        ('upper_seis_depth', float),
    ])


class ActiveFaults(object):
    URL_GEM_ACTIVE_FAULTS = 'https://raw.githubusercontent.com/cossatot/gem-global-active-faults/master/geojson/gem_active_faults.geojson'  # noqa

    def __init__(self):
        self.fname_active_faults = op.join(
            config.config().fault_lines_dir, 'gem_active_faults.geojson')

        if not op.exists(self.fname_active_faults):
                        self.download()

        self.active_faults = []
        self._load_faults(self.fname_active_faults, ActiveFault)

    def _load_faults(self, fname, cls):
        with open(fname, 'r') as f:
            gj = json.load(f)
            faults = gj['features']
            for f in faults:
                fault = cls(f)
                self.active_faults.append(fault)
        logger.debug('loaded %d fault', self.nactive_faults)

    def download(self):
        logger.info('Downloading GEM active faults database...')
        util.download_file(self.URL_GEM_ACTIVE_FAULTS,
                           self.fname_active_faults)

    @property
    def nactive_faults(self):
        return len(self.active_faults)

    def nactive_faults_nodes(self):
        return int(sum(len(f.lat) for f in self.active_faults))

    def get_coords(self):
        return num.array([f['coordinates'] for f in self.active_faults])


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    activefaults = ActiveFaults()
