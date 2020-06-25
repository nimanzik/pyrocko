from pyrocko import moment_tensor as pmt
from pyrocko.guts import Object, Float, load_all, load_all_xml
from .seismosizer import source_classes, Source, Cloneable, Location


class UncertaintySource(Location, Cloneable):
    time = Float.T(
        default=0.,
        help='Error on source origin time')

    def __init___(self, **kwargs):
        Location.__init__(self, **kwargs)


class UncertaintySourceWithMagntiude(UncertaintySource):
    magnitude = Float.T(
        default=6.0,
        help='error on moment magnitude Mw as in [Hanks and Kanamori, 1979]')

    def __init__(self, **kwargs):
        if 'moment' in kwargs:
            mom = kwargs.pop('moment')
            if 'magnitude' not in kwargs:
                kwargs['magnitude'] = float(pmt.moment_to_magnitude(mom))

        UncertaintySource.__init__(self, **kwargs)

    @property
    def moment(self):
        return float(pmt.magnitude_to_moment(self.magnitude))

    @moment.setter
    def moment(self, value):
        self.magnitude = float(pmt.moment_to_magnitude(value))


class UncertaintyDCSource(UncertaintySource):
    strike = Float.T(
        default=0.0,
        help='error on strike direction in [deg], '
             'measured clockwise from north')

    dip = Float.T(
        default=90.0,
        help='error on dip angle in [deg], measured downward from horizontal')

    rake = Float.T(
        default=0.0,
        help='error on rake angle in [deg], '
             'measured counter-clockwise from right-horizontal '
             'in on-plane view')


uncertainty_classes = [
    UncertaintySource,
    UncertaintySourceWithMagntiude,
    UncertaintyDCSource
]


def load_source_container(filename):
    return load_all(filename=filename)


def load_source_container_from_xml(filename):
    return load_all_xml(filename=filename)


class SourceContainer(Object):
    source__ = Source.T(optional=True)
    uncertainty__ = UncertaintySource.T(optional=True)

    @property
    def source(self):
        return self.source__

    @source.setter
    def source(self, source):
        if source is None:
            return

        if not (type(source) in source_classes):
            raise AttributeError('Source is not defined in pyrocko.gf.')

        self.source__ = source

    @property
    def uncertainty(self):
        return self.uncertainty__

    @uncertainty.setter
    def uncertainty(self, uncertainty):
        if not (type(uncertainty) in uncertainty_classes):
            raise AttributeError('No valid uncertainty class given.')

        if self.source__ is not None:
            if uncertainty.T.classname != (
                    'Uncertainty' + self.source.T.classname):
                raise AttributeError(
                    'Uncertainty is not matching used source.')

        self.uncertainty__ = uncertainty

    def load_grond_results(self):
        pass


__all__ = [
    'load_source_container',
    'load_source_container_from_xml',
    'SourceContainer'] + [u.__name__ for u in uncertainty_classes]
