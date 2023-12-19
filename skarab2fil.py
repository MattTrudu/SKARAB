#!/usr/bin/env python

import sys
import os
import struct
import numpy as np
from astropy.time import Time
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

def get_sampling_time(filename, nchan = 4096, npol = 2):
        file_sardara_obs = open(filename, 'rb')

        N_spectra = 100

        array_byte_seconds = np.zeros(N_spectra)

        byte_seconds_from_1970_spectrum_0    = file_sardara_obs.read(8)
        byte_microseconds_spectrum_0         = file_sardara_obs.read(8)
        array_byte_seconds[0] = get_time_stamp_from_spectrum(byte_seconds_from_1970_spectrum_0, byte_microseconds_spectrum_0)

        for k in range(1, N_spectra):
                file_sardara_obs.seek( k * (npol*nchan + 16))
                byte_seconds_from_1970    = file_sardara_obs.read(8)
                byte_microseconds         = file_sardara_obs.read(8)
                array_byte_seconds[k] = get_time_stamp_from_spectrum(byte_seconds_from_1970, byte_microseconds)

        file_sardara_obs.close()

        for k in range(1, N_spectra):
                print("array_byte_seconds[%d] - array_byte_seconds[%d] = %.4e s" % (k, k-1, array_byte_seconds[k]-array_byte_seconds[k-1] ))



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

    return parser.parse_args()

if __name__ == "__main__":

    args = _get_parser()

    filename = args.filename

    get_sampling_time(filename)


    tstamp = get_tstamp(filename, start = 0)
    print("0",tstamp.mjd)
    tstamp = get_tstamp(filename, start = 8192 - 16)
    print("8192 - 16", tstamp.mjd)
    tstamp = get_tstamp(filename, start = 8192)
    print("8192", tstamp.mjd)
    tstamp = get_tstamp(filename, start = 8192 + 16)
    print("8192 + 16", tstamp.mjd)


    #get_all_timestamps(filename)
