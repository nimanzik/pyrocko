from __future__ import division, print_function, absolute_import

import unittest
import numpy as num
from pyrocko import gf, util

km = 1e3
d2r = num.pi / 180.


class GFSourceTypesTestCase(unittest.TestCase):

    def test_rectangular_source(self):
        # WIP
        nsrc = 5
        rect_sources = []

        for n in range(nsrc):
            src = gf.RectangularSource(
                lat=0., lon=0.,
                anchor='bottom',
                north_shift=5000., east_shift=9000., depth=4.*km,
                width=2.*km, length=8.*km,
                dip=0., rake=0., strike=(180./nsrc + 1) * n,
                slip=1.)
            rect_sources.append(src)

    @staticmethod
    def plot_rectangular_source(src, store):
        from matplotlib import pyplot as plt
        from matplotlib.patches import Polygon
        ax = plt.gca()
        ne = src.outline(cs='xy')
        p = Polygon(num.fliplr(ne), fill=False, color='r', alpha=.7)
        ax.add_artist(p)

        mt = src.discretize_basesource(store)
        ax.scatter(mt.east_shifts, mt.north_shifts, alpha=1)
        ax.scatter(src.east_shift, src.north_shift, color='r')

        plt.axis('equal')
        plt.show()

    def test_double_sf_source(self):
        dsf_source = gf.DoubleSFSource(
            rfn1=0.1, rfe1=-0.3, rfd1=1.,
            rfn2=0., rfe2=.94, rfd2=-.48,
            force=1.2e6,
            mix=0.3,
            time=util.stt('2020-10-30 14:28:11.3'),
            delta_time=10., delta_depth=500.,
            distance=2000., azimuth=35.,
            lat=0., lon=0.,
            north_shift=3000., east_shift=6700., depth=2.*km)

        a1 = 1.0 - dsf_source.mix
        a2 = dsf_source.mix

        delta_north = num.cos(dsf_source.azimuth * d2r) * dsf_source.distance
        delta_east = num.sin(dsf_source.azimuth * d2r) * dsf_source.distance

        f1 = num.array([dsf_source.rfn1, dsf_source.rfe1, dsf_source.rfd1])
        f2 = num.array([dsf_source.rfn2, dsf_source.rfe2, dsf_source.rfd2])

        scale_factor = dsf_source.force / num.linalg.norm(a1 * f1 + a2 * f2)

        sf_source1 = gf.SFSource(
            lat=dsf_source.lat, lon=dsf_source.lon,
            north_shift=dsf_source.north_shift - delta_north * a2,
            east_shift=dsf_source.east_shift - delta_east * a2,
            depth=dsf_source.depth - dsf_source.delta_depth * a2,
            time=dsf_source.time - dsf_source.delta_time * a2,
            fn=f1[0] * scale_factor * a1,
            fe=f1[1] * scale_factor * a1,
            fd=f1[2] * scale_factor * a1)

        sf_source2 = gf.SFSource(
            lat=dsf_source.lat, lon=dsf_source.lon,
            north_shift=dsf_source.north_shift + delta_north * a1,
            east_shift=dsf_source.east_shift + delta_east * a1,
            depth=dsf_source.depth + dsf_source.delta_depth * a1,
            time=dsf_source.time + dsf_source.delta_time * a1,
            fn=f2[0] * scale_factor * a2,
            fe=f2[1] * scale_factor * a2,
            fd=f2[2] * scale_factor * a2)

        assert sf_source1.fn == dsf_source.fn1
        assert sf_source1.fe == dsf_source.fe1
        assert sf_source1.fd == dsf_source.fd1
        assert sf_source2.fn == dsf_source.fn2
        assert sf_source2.fe == dsf_source.fe2
        assert sf_source2.fd == dsf_source.fd2

        assert num.linalg.norm(
            num.add(*dsf_source.forces)) == dsf_source.force

        dsf_source1, dsf_source2 = dsf_source.split()
        for p in sf_source1.T.propnames:
            assert getattr(sf_source1, p) == getattr(dsf_source1, p)
            assert getattr(sf_source2, p) == getattr(dsf_source2, p)


if __name__ == '__main__':
    util.setup_logging('test_gf_source_types', 'warning')
    unittest.main(defaultTest='GFSourceTypesTestCase')
