#!/usr/bin/env python

import sys
import os
import struct
import numpy as np
from astropy.time import Time
import astropy.units as u



def get_time_stamp_from_spectrum(byte_seconds_from_1970, byte_microseconds):
        str_seconds_from_1970 = struct.unpack('Q', byte_seconds_from_1970)[0]
        str_byte_microseconds = struct.unpack('Q', byte_microseconds     )[0]


        int_seconds_from_1970    = int(str_seconds_from_1970)
        int_microseconds         = int(str_byte_microseconds)

        float_seconds_from_1970  = np.float64(int_seconds_from_1970)
        float_milliseconds       = np.float64(int_microseconds) / 1000.
        float_fractional_second  = float_milliseconds / 1000.

        float_seconds_from_1970 = float_seconds_from_1970 + float_fractional_second

        return float_seconds_from_1970

class skarabrawfile:

    def __init__(
    self,
    filename,
    filepath,
    projID,
    telescope,
    source_name,
    point_ra,
    point_dec,
    channel_band_MHz,
    freq_top_MHz,
    nchans,
    bit_depth,
    tsamp_us,
    nspectra_per_bin
    ):

        self.filename = filename
        self.filepath = filepath
        self.projID = projID
        self.source_name = source_name
        self.point_ra  = point_ra
        self.point_dec = point_dec
        self.telescope = telescope
        self.channel_band_MHz = channel_band_MHz
        self.freq_top_MHz = freq_top_MHz
        self.nchans = nchans
        self.bit_depth = bit_depth
        self.tsamp_us = tsamp_us
        self.nspectra_per_bin = nspectra_per_bin

    def get_tstamp(self, start = 0):

        try:
            file_abspath = os.path.abspath(os.path.join(self.filepath, self.filename ) )
            with open(file_abspath, 'rb') as file_skarab_obs:
                file_skarab_obs.seek(start)
                byte_seconds_from_1970_spectrum_0 = file_skarab_obs.read(self.bit_depth)
                byte_microseconds_spectrum_0 = file_skarab_obs.read(self.bit_depth)


        except FileNotFoundError:
            print(f"File '{filename}' not found.")

        float_seconds_from_1970 = get_time_stamp_from_spectrum(byte_seconds_from_1970_spectrum_0, byte_microseconds_spectrum_0)


        tstamp = Time(float_seconds_from_1970, format='unix')

        return tstamp

    def get_all_timestamps(self):

        filename = self.filename
        filepath = self.filepath
        nchans = self.nchans
        npols = self.nspectra_per_bin

        timestamps = []
        total_length = os.path.getsize(os.path.join(filepath, filename))
        start = 0

        chunk_size = npols * nchans

        try:
            while start < total_length:
                tstamp = self.get_tstamp(start)
                if tstamp is None:
                    break  # End of file

                timestamps.append(tstamp)
                start += chunk_size + 2 * self.bit_depth

        except FileNotFoundError:
            print(f"File '{self.filename}' not found.")

        timestamps = np.array(timestamps)

        return timestamps

    def get_spectra_per_bin(self, start = 0, mode = "XX,YY,Re(XY),Im(YX)"):

        spectrum_format = f'<{self.nchans}b'
        spectrum_size   = struct.calcsize(spectrum_format)

        try:
            file_abspath = os.path.abspath(os.path.join(self.filepath, self.filename ) )
            with open(file_abspath, 'rb') as file_skarab_obs:
                file_skarab_obs.seek(start)
                byte_seconds_from_1970_spectrum_0 = file_skarab_obs.read(self.bit_depth)
                byte_microseconds_spectrum_0 = file_skarab_obs.read(self.bit_depth)
                if mode == "XX,YY,Re(XY),Im(YX)":
                    spectrum_bytes_xx = file_skarab_obs.read(spectrum_size)
                    spectrum_bytes_yy = file_skarab_obs.read(spectrum_size)
                    spectrum_bytes_xy = file_skarab_obs.read(spectrum_size)
                    spectrum_bytes_yx = file_skarab_obs.read(spectrum_size)
                    spectrum_xx = struct.unpack(f'<{self.nchans}B', spectrum_bytes_xx) # B takes the absolute value apparently
                    spectrum_yy = struct.unpack(f'<{self.nchans}B', spectrum_bytes_yy)
                    spectrum_xy = struct.unpack(f'<{self.nchans}b', spectrum_bytes_xy) # b no
                    spectrum_yx = struct.unpack(f'<{self.nchans}b', spectrum_bytes_yx)

                    return spectrum_xx, spectrum_yy, spectrum_xy, spectrum_yx

        except FileNotFoundError:
            print(f"File '{self.filename}' not found.")


    def get_intensity_dynspec(self, mode = "XX,YY,Re(XY),Im(YX)"):


        filename = self.filename
        filepath = self.filepath
        nchans = self.nchans
        npols = self.nspectra_per_bin

        dynspec = []
        total_length = os.path.getsize(os.path.join(filepath, filename))
        start = 0

        chunk_size = npols * nchans

        try:
            while start < total_length:
                spectrum_xx, spectrum_yy, spectrum_xy, spectrum_yx = self.get_spectra_per_bin(start = start, mode = mode)

                dynspec.append(spectrum_xx + spectrum_yy)
                start += chunk_size + 2 * self.bit_depth

        except FileNotFoundError:
            print(f"File '{self.filename}' not found.")

        dynspec = np.asarray(dynspec)

        return dynspec
