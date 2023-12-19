#!/usr/bin/env python

import sys
import os
import struct
import numpy as np
from astropy.time import Time
import astropy.units as u
import argparse
import matplotlib.pyplot as plt


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

def get_tstamp(filename, start = 0):

    try:
        file_abspath = os.path.abspath(filename)
        with open(file_abspath, 'rb') as file_skarab_obs:
            file_skarab_obs.seek(start)
            byte_seconds_from_1970_spectrum_0 = file_skarab_obs.read(8)
            byte_microseconds_spectrum_0 = file_skarab_obs.read(8)


    except FileNotFoundError:
        print(f"File '{filename}' not found.")

    #file_skarab_obs.close()

    float_seconds_from_1970 = get_time_stamp_from_spectrum(byte_seconds_from_1970_spectrum_0, byte_microseconds_spectrum_0)


    tstamp = Time(float_seconds_from_1970, format='unix')

    return tstamp

def get_all_timestamps(filename, nchans = 2048, npols = 2):
    timestamps = []
    total_length = os.path.getsize(filename)
    start = 0

    chunk_size = npols * nchans

    try:
        while start < total_length:
            tstamp = get_tstamp(filename, start)
            if tstamp is None:
                break  # End of file

            timestamps.append(tstamp)
            start += chunk_size + 16

    except FileNotFoundError:
        print(f"File '{filename}' not found.")

    timestamps = np.array(timestamps)

    return timestamps



def _get_parser():
    """
    Argument parser.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Read a SKARAB raw data file and convert it into a SIGPROC filterbank file.",
    )


    parser.add_argument(
        "-f",
        "--filename",
        action = "store",
        help = "Raw data file to be converted.",
        required=True,
    )

    parser.add_argument(
        "-nc",
        "--nchans",
        type = int,
        help = "Number of spectral channels of the raw data.",
        default = 2048,
    )


    return parser.parse_args()

if __name__ == "__main__":

    args = _get_parser()

    filename = args.filename
    nchans = args.nchans

    tstamps = get_all_timestamps(filename, nchans = nchans )

    tstamps = tstamps

    dts = np.diff(tstamps)

    dts_us = []
    for dt in dts:
        dt = dt.to_value('s', 'long')
        dt = dt * u.s
        dts_us.append(dt.to(u.us).value)

    dts_us = np.array(dts_us)
    print(dts_us.mean())
    print(dts_us.std())

    plt.figure()
    plt.hist(dts_us, bins = 100)
    plt.xlabel("Time (us)")
    plt.show()
