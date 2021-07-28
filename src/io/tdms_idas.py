import logging
import os
import struct
import datetime
import mmap
import numpy as num
import pandas as pd

from pyrocko import trace
from .io_common import FileLoadError

logger = logging.getLogger(__name__)


def write_property_dict(prop_dict, out_file):
    from pprint import pformat

    f = open(out_file, "w")
    f.write("tdms_property_map=")
    f.write(pformat(prop_dict))
    f.close()


def type_not_supported(vargin):
    """Function raises a NotImplementedException."""
    raise NotImplementedError("Reading of this tdsDataType is not implemented")


def parse_time_stamp(fractions, seconds):
    """
    Convert time TDMS time representation to datetime
    fractions   -- fractional seconds (2^-64)
    seconds     -- The number of seconds since 1/1/1904
    @rtype : datetime.datetime
    """
    if (
        fractions is not None
        and seconds is not None
        and fractions + seconds > 0
    ):
        return datetime.timedelta(
            0, fractions * 2 ** -64 + seconds
        ) + datetime.datetime(1904, 1, 1)
    else:
        return None


# Enum mapping TDM data types to description string, numpy type where exists
# See Ref[2] for enum values
TDS_DATA_TYPE = dict(
    {
        0x00: "void",  # tdsTypeVoid
        0x01: "int8",  # tdsTypeI8
        0x02: "int16",  # tdsTypeI16
        0x03: "int32",  # tdsTypeI32
        0x04: "int64",  # tdsTypeI64
        0x05: "uint8",  # tdsTypeU8
        0x06: "uint16",  # tdsTypeU16
        0x07: "uint32",  # tdsTypeU32
        0x08: "uint64",  # tdsTypeU64
        0x09: "float32",  # tdsTypeSingleFloat
        0x0A: "float64",  # tdsTypeDoubleFloat
        0x0B: "float128",  # tdsTypeExtendedFloat
        0x19: "singleFloatWithUnit",  # tdsTypeSingleFloatWithUnit
        0x1A: "doubleFloatWithUnit",  # tdsTypeDoubleFloatWithUnit
        0x1B: "extendedFloatWithUnit",  # tdsTypeExtendedFloatWithUnit
        0x20: "str",  # tdsTypeString
        0x21: "bool",  # tdsTypeBoolean
        0x44: "datetime",  # tdsTypeTimeStamp
        0xFFFFFFFF: "raw",  # tdsTypeDAQmxRawData
    }
)

# Function mapping for reading TDMS data types
TDS_READ_VAL = dict(
    {
        "void": lambda f: None,  # tdsTypeVoid
        "int8": lambda f: struct.unpack("<b", f.read(1))[0],
        "int16": lambda f: struct.unpack("<h", f.read(2))[0],
        "int32": lambda f: struct.unpack("<i", f.read(4))[0],
        "int64": lambda f: struct.unpack("<q", f.read(8))[0],
        "uint8": lambda f: struct.unpack("<B", f.read(1))[0],
        "uint16": lambda f: struct.unpack("<H", f.read(2))[0],
        "uint32": lambda f: struct.unpack("<I", f.read(4))[0],
        "uint64": lambda f: struct.unpack("<Q", f.read(8))[0],
        "float32": lambda f: struct.unpack("<f", f.read(4))[0],
        "float64": lambda f: struct.unpack("<d", f.read(8))[0],
        "float128": type_not_supported,
        "singleFloatWithUnit": type_not_supported,
        "doubleFloatWithUnit": type_not_supported,
        "extendedFloatWithUnit": type_not_supported,
        "str": lambda f: f.read(struct.unpack("<i", f.read(4))[0]),
        "bool": lambda f: struct.unpack("<?", f.read(1))[0],
        "datetime": lambda f: parse_time_stamp(
            struct.unpack("<Q", f.read(8))[0],
            struct.unpack("<q", f.read(8))[0],
        ),
        "raw": type_not_supported,
    }
)

DECIMATE_MASK = 0b00100000
LEAD_IN_LENGTH = 28
FILEINFO_NAMES = (
    "file_tag",
    "toc",
    "version",
    "next_segment_offset",
    "raw_data_offset",
)


class TdmsReader(object):
    """A TDMS file reader object for reading properties and data"""

    def __init__(self, filename):
        self._properties = None
        self._end_of_properties_offset = None
        self._data_type = None
        self._chunk_size = None

        self._raw_data = None
        self._raw_data2 = None  # The mapped data in the 'Next Segment'
        self._raw_last_chunk = None
        self._raw2_last_chunk = None

        self.file_size = os.path.getsize(filename)
        self._channel_length = None
        self._seg1_length = None
        self._seg2_length = None

        # TODO: Error if file not big enough to hold header
        self._tdms_file = open(filename, "rb")
        # Read lead in (28 bytes):
        lead_in = self._tdms_file.read(LEAD_IN_LENGTH)
        # lead_in is 28 bytes:
        # [string of length 4][int32][int32][int64][int64]
        fields = struct.unpack("<4siiQQ", lead_in)

        # TODO: validate file
        if fields[0].decode() not in "TDSm":
            msg = "Not a TDMS file (TDSm tag not found)"
            raise (TypeError, msg)

        self.fileinfo = dict(zip(FILEINFO_NAMES, fields))
        self.fileinfo["decimated"] = not bool(
            self.fileinfo["toc"] & DECIMATE_MASK
        )
        # Make offsets relative to beginning of file:
        self.fileinfo["next_segment_offset"] += LEAD_IN_LENGTH
        self.fileinfo["raw_data_offset"] += LEAD_IN_LENGTH
        self.fileinfo["file_size"] = os.path.getsize(self._tdms_file.name)

        # TODO: Validate lead in:
        if self.fileinfo["next_segment_offset"] > self.file_size:
            self.fileinfo["next_segment_offset"] = self.file_size
            # raise(ValueError, "Next Segment Offset too large in TDMS header")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._tdms_file.close()

    def _get_channel_length(self):
        if not self._channel_length:
            self._initialise_data()

        return self._channel_length

    channel_length = property(_get_channel_length)

    def get_properties(self, mapped=False):
        """
        Return a dictionary of properties. Read from file only if necessary.
        """
        # Check if already hold properties in memory
        if self._properties is None:
            self._properties = self._read_properties()
        else:
            return self._properties.loc[:, "Value"].to_dict()

    def _read_property(self):
        """
        Read a single property from the TDMS file.
        Return the name, type and value of the property as a list.
        """
        # Read length of object path:
        var = struct.unpack("<i", self._tdms_file.read(4))[0]
        # Read property name and type:
        name, data_type = struct.unpack(
            "<{0}si".format(var), self._tdms_file.read(var + 4)
        )
        # Lookup function to read and parse property value based on type:
        value = TDS_READ_VAL[TDS_DATA_TYPE[data_type]](self._tdms_file)
        name = name.decode()
        if data_type == 32:
            value = value.decode()

        return name, data_type, value

    def _read_properties(self):
        """Read the properties from the file"""
        self._tdms_file.seek(LEAD_IN_LENGTH, 0)
        # Number of channels is total objects - file objects - group objects
        self.fileinfo["n_channels"] = (
            struct.unpack("i", self._tdms_file.read(4))[0] - 2
        )
        # Read length of object path:
        var = struct.unpack("<i", self._tdms_file.read(4))[0]
        # skip over object path and raw data index:
        self._tdms_file.seek(var + 4, 1)
        # Read number of properties in this group:
        var = struct.unpack("<i", self._tdms_file.read(4))[0]

        # loop through and read each property
        properties = [self._read_property() for _ in range(var)]
        df = pd.DataFrame(properties)
        df.columns = ["Property", "Type", "Value"]
        df.set_index("Property", inplace=True)

        self._end_of_properties_offset = self._tdms_file.tell()

        self._read_chunk_size()
        # TODO: Add number of channels to properties
        return df

    def _read_chunk_size(self):
        """Read the data chunk size from the TDMS file header."""
        if self._end_of_properties_offset is None:
            self._read_properties()

        self._tdms_file.seek(self._end_of_properties_offset, 0)

        # skip over Group Information:
        var = struct.unpack("<i", self._tdms_file.read(4))[0]
        self._tdms_file.seek(var + 8, 1)

        # skip over first channel path and length of index information:
        var = struct.unpack("<i", self._tdms_file.read(4))[0]
        self._tdms_file.seek(var + 4, 1)

        self._data_type = TDS_DATA_TYPE.get(
            struct.unpack("<i", self._tdms_file.read(4))[0]
        )
        if self._data_type not in ("int16", "float32"):
            raise Exception("Unsupported TDMS data type: " + self._data_type)

        # Read Dimension of the raw data array (has to be 1):
        # dummy = struct.unpack("<i", self._tdms_file.read(4))[0]

        self._chunk_size = struct.unpack("<i", self._tdms_file.read(4))[0]

    def get_data(self, first_ch=0, last_ch=None, first_s=0, last_s=None):
        """
        Get a block of data from the TDMS file.
        first_ch -- The first channel to load
        last_ch  -- The last channel to load
        first_s  -- The first sample to load
        last_s   -- The last sample to load
        """
        if self._raw_data is None:
            self._initialise_data()
        if first_ch is None or first_ch < 0:
            first_ch = 0
        if last_ch is None or last_ch >= self.fileinfo["n_channels"]:
            last_ch = self.fileinfo["n_channels"]
        else:
            last_ch += 1
        if last_s is None or last_s > self._channel_length:
            last_s = self._channel_length
        else:
            last_s += 1
        nch = num.int(max(last_ch - first_ch, 0))
        ns = num.int(max(last_s - first_s, 0))

        # Allocate output container
        data = num.empty((ns, nch), dtype=num.dtype(self._data_type))
        if data.size == 0:
            return data

        # 1. Index first block & reshape?
        first_blk = first_s // self._chunk_size
        last_blk = last_s // self._chunk_size
        last_full_blk = min(last_blk + 1, self._raw_data.shape[1])
        nchunk = min(
            max(last_full_blk - first_blk, 0), self._raw_data.shape[1]
        )
        first_s_1a = max(first_s - first_blk * self._chunk_size, 0)
        last_s_1a = min(
            last_s - first_blk * self._chunk_size, nchunk * self._chunk_size
        )
        ind_s = 0
        ind_e = ind_s + max(last_s_1a - first_s_1a, 0)

        # data_1a = self._raw_data[:, first_blk:last_full_blk, first_ch:last_ch] \
        #   .reshape((self._chunk_size*nchunk, nch), order='F')[first_s_1a:last_s_1a, :]
        d = self._raw_data[:, first_blk:last_full_blk, first_ch:last_ch]
        d.shape = (self._chunk_size * nchunk, nch)
        d.reshape((self._chunk_size * nchunk, nch), order="F")
        data[ind_s:ind_e, :] = d[first_s_1a:last_s_1a, :]

        # 2. Index first additional samples
        first_s_1b = max(
            first_s - self._raw_data.shape[1] * self._chunk_size, 0
        )
        last_s_1b = min(
            last_s - self._raw_data.shape[1] * self._chunk_size,
            self._raw_last_chunk.shape[0],
        )
        ind_s = ind_e
        ind_e = ind_s + max(last_s_1b - first_s_1b, 0)
        # data_1b = self._raw_last_chunk[first_s_1b:last_s_1b,first_ch:last_ch]
        if ind_e > ind_s:
            data[ind_s:ind_e, :] = self._raw_last_chunk[
                first_s_1b:last_s_1b, first_ch:last_ch
            ]

        # 3. Index second block
        first_s_2 = max(first_s - self._seg1_length, 0)
        last_s_2 = last_s - self._seg1_length
        if (first_s_2 > 0 or last_s_2 > 0) and self._raw_data2 is not None:
            first_blk_2 = max(first_s_2 // self._chunk_size, 0)
            last_blk_2 = max(last_s_2 // self._chunk_size, 0)
            last_full_blk_2 = min(last_blk_2 + 1, self._raw_data2.shape[1])
            nchunk_2 = min(
                max(last_full_blk_2 - first_blk_2, 0), self._raw_data2.shape[1]
            )
            first_s_2a = max(first_s_2 - first_blk_2 * self._chunk_size, 0)
            last_s_2a = min(
                last_s_2 - first_blk_2 * self._chunk_size,
                nchunk_2 * self._chunk_size,
            )
            ind_s = ind_e
            ind_e = ind_s + max(last_s_2a - first_s_2a, 0)
            # data_2a = self._raw_data2[:, first_blk_2:last_full_blk_2,
            #             first_ch:last_ch]\
            #   .reshape((self._chunk_size*nchunk_2, nch), order='F')\
            #   [first_s_2a:last_s_2a, :]
            if ind_e > ind_s:
                data[ind_s:ind_e, :] = self._raw_data2[
                    :, first_blk_2:last_full_blk_2, first_ch:last_ch
                ].reshape((self._chunk_size * nchunk_2, nch), order="F")[
                    first_s_2a:last_s_2a, :
                ]
        # 4. Index second additional samples
        if (
            first_s_2 > 0 or last_s_2 > 0
        ) and self._raw2_last_chunk is not None:
            first_s_2b = max(
                first_s_2 - self._raw_data2.shape[1] * self._chunk_size, 0
            )
            last_s_2b = min(
                last_s_2 - self._raw_data2.shape[1] * self._chunk_size,
                self._raw2_last_chunk.shape[0],
            )
            ind_s = ind_e
            ind_e = ind_s + max(last_s_2b - first_s_2b, 0)
            # data_2b = \
            #   self._raw2_last_chunk[first_s_2b:last_s_2b,first_ch:last_ch]
            if ind_e > ind_s:
                data[ind_s:ind_e, :] = self._raw2_last_chunk[
                    first_s_2b:last_s_2b, first_ch:last_ch
                ]
        # 5. Concatenate blocks
        # data = num.concatenate((data_1a, data_1b, data_2a, data_2b))
        if data.size == 0:
            data = data.reshape(0, 0)
        return data

    def _initialise_data(self):
        """Initialise the memory map for the data array."""
        if self._chunk_size is None:
            self._read_chunk_size()

        dmap = mmap.mmap(self._tdms_file.fileno(), 0, access=mmap.ACCESS_READ)
        rdo = num.int(self.fileinfo["raw_data_offset"])
        nch = num.int(self.fileinfo["n_channels"])

        # TODO: Support streaming file type?
        # TODO: Is this a valid calculation for ChannelLength?
        nso = self.fileinfo["next_segment_offset"]
        self._seg1_length = num.int(
            (nso - rdo) / nch / num.dtype(self._data_type).itemsize
        )
        self._channel_length = self._seg1_length

        if self.fileinfo["decimated"]:
            n_complete_blk = num.int(self._seg1_length / self._chunk_size)
            ax_ord = "C"
        else:
            n_complete_blk = 0
            ax_ord = "F"
        self._raw_data = num.ndarray(
            (n_complete_blk, nch, self._chunk_size),
            dtype=self._data_type,
            buffer=dmap,
            offset=rdo,
        )
        # Rotate the axes to [chunk_size, nblk, nch]
        self._raw_data = num.rollaxis(self._raw_data, 2)
        additional_samples = num.int(
            self._seg1_length - n_complete_blk * self._chunk_size
        )
        additional_samples_offset = (
            rdo
            + n_complete_blk
            * nch
            * self._chunk_size
            * num.dtype(self._data_type).itemsize
        )
        self._raw_last_chunk = num.ndarray(
            (nch, additional_samples),
            dtype=self._data_type,
            buffer=dmap,
            offset=additional_samples_offset,
            order=ax_ord,
        )
        # Rotate the axes to [samples, nch]
        self._raw_last_chunk = num.rollaxis(self._raw_last_chunk, 1)

        if self.file_size == nso:
            self._seg2_length = 0
        else:
            self._tdms_file.seek(nso + 12, 0)
            (seg2_nso, seg2_rdo) = struct.unpack(
                "<qq", self._tdms_file.read(2 * 8)
            )
            self._seg2_length = (
                (seg2_nso - seg2_rdo)
                / nch
                / num.dtype(self._data_type).itemsize
            )
            if self.fileinfo["decimated"]:
                n_complete_blk2 = num.int(self._seg2_length / self._chunk_size)
            else:
                n_complete_blk2 = num.int(0)
            self._raw_data2 = num.ndarray(
                (n_complete_blk2, nch, self._chunk_size),
                dtype=self._data_type,
                buffer=dmap,
                offset=(nso + LEAD_IN_LENGTH + seg2_rdo),
            )
            self._raw_data2 = num.rollaxis(self._raw_data2, 2)
            additional_samples = num.int(
                self._seg2_length - n_complete_blk2 * self._chunk_size
            )
            additional_samples_offset = (
                nso
                + LEAD_IN_LENGTH
                + seg2_rdo
                + n_complete_blk2
                * nch
                * self._chunk_size
                * num.dtype(self._data_type).itemsize
            )
            self._raw2_last_chunk = num.ndarray(
                (nch, additional_samples),
                dtype=self._data_type,
                buffer=dmap,
                offset=additional_samples_offset,
                order=ax_ord,
            )
            # Rotate the axes to [samples, nch]
            self._raw2_last_chunk = num.rollaxis(self._raw2_last_chunk, 1)

            if self._raw_data2.size != 0 or self._raw2_last_chunk.size != 0:
                pass
                # raise Exception('Second segment contains some data, \
                #                 not currently supported')
            self._channel_length = self._seg1_length + self._seg2_length
        # else:
        #     print "Not decimated"
        #     raise Exception('Reading file with decimated flag not set is not'
        #                     ' supported yet')


if __name__ == "__main__":
    print("TDMS Reader demo.")

    file_path = "path_to_tdms_file.tdms"

    print("File: {0}".format(file_path))

    tdms = TdmsReader(file_path)
    props = tdms.get_properties()

    zero_offset = props.get("Zero Offset (m)")
    channel_spacing = props.get("SpatialResolution[m]") * props.get(
        "Fibre Length Multiplier"
    )
    n_channels = tdms.fileinfo["n_channels"]
    depth = zero_offset + num.arange(n_channels) * channel_spacing
    fs = props.get("SamplingFrequency[Hz]")

    print("Number of channels in file: {0}".format(n_channels))
    print("Time samples in file: {0}".format(tdms.channel_length))
    print("Sampling frequency (Hz): {0}".format(fs))

    first_channel = 250
    last_channel = 2275
    first_time_sample = 0
    last_time_sample = 3999

    some_data = tdms.get_data(
        first_channel, last_channel, first_time_sample, last_time_sample
    )
    print("Size of data loaded: {0}".format(some_data.shape))

    import matplotlib.pyplot as plt

    fig1 = plt.figure()
    img1 = plt.imshow(
        some_data,
        aspect="auto",
        interpolation="none",
        extent=(
            depth[first_channel],
            depth[last_channel],
            last_time_sample / fs,
            first_time_sample / fs,
        ),
    )
    plt.ylabel("Time (seconds)")
    plt.xlabel("Depth (m)")
    plt.show(block=False)


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
    'continuous_mode': 'Continuous Mode',
    'geo_lat': 'SystemInfomation.GPS.Latitude',
    'geo_lon': 'SystemInfomation.GPS.Longitude',
    'geo_elevation': 'SystemInfomation.GPS.Altitude',

    'channel': None,
    'unit': None
}


def get_meta(tdms_properties):
    prop = tdms_properties

    deltat = 1. / prop['SamplingFrequency[Hz]']
    tmin = prop['GPSTimeStamp'].astype(num.float) * 1e-6

    fibre_meta = {key: prop.get(key_map, -1)
                  for key, key_map in META_KEYS.items()
                  if key_map is not None}

    coeff = fibre_meta['p_coefficients']
    coeff = tuple(map(float, coeff.split('\t')))

    gain = fibre_meta['receiver_gain']
    gain = tuple(map(float, gain.split('\t')))
    fibre_meta['receiver_gain'] = coeff

    fibre_meta['unit'] = 'radians'

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

            meta_cha = meta.copy()

            assert int(channel.name) < 99999
            station = '%05i' % int(channel.name)
            meta_cha['channel'] = int(channel.name)
            nsamples = channel._length

            tr = trace.Trace(
                network='DA',
                station=station,
                ydata=None,
                deltat=deltat,
                tmin=tmin,
                tmax=tmin + nsamples*deltat,
                meta=meta_cha)

            if load_data:
                tr.set_ydata(channel[:])

            yield tr


def detect(first512):
    return first512.startswith(b'TDSm.')
