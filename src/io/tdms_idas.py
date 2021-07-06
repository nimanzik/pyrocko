import logging

import numpy as num

from pyrocko import trace
from .io_common import FileLoadError

logger = logging.getLogger(__name__)


META_KEYS = {
    'measure_length': 'MeasureLength[m]',
    'start_position': 'StartPosition[m]',
    'spatial_resolution': 'SpatialResolution[m]',
    'fibre_index': 'FibreIndex',
    'unit_calibration': 'Unit Calibration (nm)',
    'start_distance': 'Start Distance (m)',
    'stop_distance': 'Stop Distance (m)',
    'normalization': 'Normalization',
    'decimation_filter': 'Decimation Filter',
    'gauge_length': 'GaugeLength',
    'norm_offset': 'Norm Offset',
    'source_mode': 'Source Mode',
    'time_decimation': 'Time Decimation',
    'zero_offset': 'Zero Offset (m)',
    'p_parameter': 'P',
    'p_coefficients': 'P Coefficients',
    'idas_version': 'iDASVersion',
    'precice_sampling_freq': 'Precise Sampling Frequency (Hz)',
    'receiver_gain': 'Receiver Gain',
    'continuous_mode': 'Continuous Mode'
}


def get_meta(tdms_properties):
    prop = tdms_properties

    deltat = 1. / prop['SamplingFrequency[Hz]']
    tmin = prop['GPSTimeStamp'].astype(num.float) * 1e-6

    fibre_meta = {key: prop.get(key_map, -1)
                  for key, key_map in META_KEYS.items()}

    coeff = fibre_meta['p_coefficients']
    coeff = tuple(map(float, coeff.split(';')))

    gain = fibre_meta['receiver_gain']
    gain = tuple(map(int, gain.split(';')))
    fibre_meta['receiver_gain'] = coeff

    return deltat, tmin, fibre_meta


def iload(filename, load_data=True):
    try:
        from nptdms import TdmsFile
    except ImportError:
        raise FileLoadError('Could not import Python module nptdms.')

    try:
        if load_data:
            tdms = TdmsFile.read(filename)
        else:
            tdms = TdmsFile.read_metadata(filename)
    except ValueError as e:
        raise FileLoadError('Cannot load %s: %s' % (filename, str(e)))

    deltat, tmin, meta = get_meta(tdms.properties)

    for group in tdms.groups():
        for channel in group.channels():

            assert int(channel.name) < 99999
            station = '%05i' % int(channel.name)
            nsamples = channel._length

            tr = trace.Trace(
                network='DA',
                station=station,
                ydata=None,
                deltat=deltat,
                tmin=tmin,
                tmax=tmin + nsamples*deltat,
                meta=meta)

            if load_data:
                tr.set_ydata(channel[:])

            yield tr


def detect(first512):
    return first512.startswith(b'TDSm.')
