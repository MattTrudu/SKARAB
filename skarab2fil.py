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

def get_tstart(filename):

    """
    file_abspath = os.path.abspath(filename)
    file_sardara_obs = open(file_abspath, 'rb')

    byte_seconds_from_1970_spectrum_0 = file_sardara_obs.read(8)
    byte_microseconds_spectrum_0      = file_sardara_obs.read(8)
    file_sardara_obs.close()
    """

    try:
        file_abspath = os.path.abspath(filename)
        with open(file_abspath, 'rb') as file_sardara_obs:
            byte_seconds_from_1970_spectrum_0 = file_sardara_obs.read(8)
            byte_microseconds_spectrum_0 = file_sardara_obs.read(8)

    except FileNotFoundError:
        print(f"File '{filename}' not found.")

    float_seconds_from_1970 = get_time_stamp_from_spectrum(byte_seconds_from_1970_spectrum_0, byte_microseconds_spectrum_0)


    tstart = Time(float_seconds_from_1970, format='unix')

    print("Tstart:", tstart.mjd)

    return tstart

def check_sampling_time(filename, nchans = 4096, npol = 2):

    bytes_per_tstamp  = 16
    bytes_per_spectra = npol * nchans
    chunk_size = bytes_per_spectra + bytes_per_tstamp

    with open(filename, 'rb') as file:
        while True:
            chunk = file.read(bytes_per_tstamp)

            if not chunk:
                break  # End of file

            # Process the current chunk (16 bytes)
            data_to_process = chunk[:bytes_per_tstamp]
            print("Processing:", data_to_process)
            byte_seconds_from_1970_spectrum_0 = chunk.read(8)
            byte_microseconds_spectrum_0 = chunk.read(8)


            float_seconds_from_1970 = get_time_stamp_from_spectrum(byte_seconds_from_1970_spectrum_0, byte_microseconds_spectrum_0)


            tstart = Time(float_seconds_from_1970, format='unix')

            print("Tstart:", tstart.mjd)
            # Move the file pointer to the next chunk
            file.seek(chunk_size - bytes_to_read, 1)



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

    #tstart = get_tstart(filename)

    check_sampling_time(filename, nchans = 4096, npol = 2)
