#!/usr/bin/env python

import sys
import os
import struct
import numpy as np
from astropy.time import Time
import astropy.units as u
import argparse
import matplotlib.pyplot as plt
from skarab import skarabrawfile
from sigpyproc import Header, Filterbank
from astropy.coordinates import SkyCoord


def convert_to_radians(rawdatafile):

    ra_str  = rawdatafile.point_ra
    dec_str = rawdatafile.point_dec

    coords = SkyCoord(ra=ra_str, dec=dec_str, unit=(u.hourangle, u.deg))

    # Accessing the right ascension and declination in radians
    ra_rad = coords.ra.rad
    dec_rad = coords.dec.rad

    return ra_rad, dec_rad

def convert_to_hoursdeg(rawdatafile):

    ra_str  = rawdatafile.point_ra
    dec_str = rawdatafile.point_dec

    coords = SkyCoord(ra=ra_str, dec=dec_str, unit=(u.hourangle, u.deg))

    # Accessing the right ascension and declination in radians
    ra_hour = coords.ra.hour
    dec_deg = coords.dec.deg

    return ra_hour, dec_deg

def make_sigpyproc_header(rawdatafile):

    dynspec = rawdatafile.get_intensity_dynspec()

    nsamp = dynspec.shape[0]

    mjdstart = rawdatafile.get_tstamp()
    mjdstart = mjdstart.mjd

    ra_hour, dec_deg = convert_to_hoursdeg(rawdatafile)


    header = Header.Header(  { "telescope_id":     10, #Matteo: it always write SRT now, I'll change it later...
                               "machine_id":       0, # DONE
                               "data_type":        1, # DONE
                               "source_name":      rawdatafile.source_name, # DONE
                               "src_raj":          float(rawdatafile.point_ra.replace(":", "")) ,
                               "src_dej":          float(rawdatafile.point_dec.replace(":", "")),
                               "tstart":           mjdstart  ,
                               "tsamp":            rawdatafile.tsamp_us * 1e-6, # DONE
                               "nsamples":         nsamp, # DONE
                               "nifs":             1,
                               "fch1":             rawdatafile.freq_top_MHz, # DONE        ,
                               "foff":             -np.abs(rawdatafile.channel_band_MHz), # DONE
                               "nchans":           rawdatafile.nchans, # DONE            ,
          }
                         )
    return header




def _get_parser():
    """
    Argument parser.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Read a SKARAB raw data file and convert it into a SIGPROC filterbank file.",
    )


#    parser.add_argument(
#        "-f",
#        "--filename",
#        #nargs = "*",
#        action = "store",
#        help = "Raw data files to be converted",
#        required=True,
#    )
    parser.add_argument(
        "-f",
        "--filenames",
        nargs='*',
        action="store",
        help="List of raw data files to be converted either into a single merged filterbank (default one) or into a list of filterbanks",
        required=True,
    )
    parser.add_argument(
        "-m",
        "--mode",
        action = "store",
        help = "Converting mode. Merge the list of raw into a single filterbank or Separate into a respective list of filterbank (Merge or Separate)",
        default = "Merge",
    )
    parser.add_argument(
        "-nc",
        "--nchans",
        type = int,
        help = "Number of spectral channels of the raw data",
        default = 2048,
    )
    parser.add_argument(
        "-np",
        "--npols",
        type = int,
        help = "Number of polarisations",
        default = 1,
    )
    parser.add_argument('-s',
                        '--source_name',
                        action = "store" ,
                        help = "Name of the source to store",
                        default = "Test"
                        )
    parser.add_argument('-r',
                        '--ra_source',
                        action = "store" ,
                        help = "R.A of the source (format \"hh:mm:ss.ss\")",
                        default = "00:00:00.00"
                        )
    parser.add_argument('-d',
                        '--dec_source',
                        action = "store" ,
                        help = "Dec. of the source (format \"dd:mm:ss.ss\")",
                        default = "00:00:00.00"
                        )
    parser.add_argument('-t',
                        '--telescope',
                        action = "store" ,
                        help = "Name of the telescope to store",
                        default = "SRT"
                        )
    parser.add_argument('-p',
                        '--polarisation',
                        action = "store" ,
                        help = "Polarisation to pick (Left, Right, Both)",
                        default = "Both"
                        )
    parser.add_argument('-pr',
                        '--project',
                        action = "store" ,
                        help = "Project ID",
                        default = "Test"
                        )
    parser.add_argument('-dt',
                        '--time_resolution',
                        type = float,
                        default = 64,
                        action = "store" ,
                        help = "Time resolution (in us) to store"
                        )
    parser.add_argument('-df',
                        '--frequency_resolution',
                        type = float,
                        default = 0.1,
                        action = "store" ,
                        help = "Frequency resolution (in MHz) to store"
                        )
    parser.add_argument('-ft',
                        '--frequency_top',
                        type = float,
                        default = 1500,
                        action = "store" ,
                        help = "Highest observational Frequency (in MHz) ot the band to store"
                        )
    parser.add_argument(
                        "-bd",
                        "--bit_depth",
                        type = int,
                        help = "Bit depth to store the data (current values allowed 8,16 and 32)",
                        default = 8,
                        )
    parser.add_argument('-o',
                        '--output_dir',
                        action = "store" ,
                        help = "Output directory (Default: your current path)",
                        default = "%s/"%(os.getcwd())
                        )
    parser.add_argument('-n',
                        '--output_names',
                        action = "store" ,
                        nargs = "+",
                        help = "Output File Names: (Default: filename_original.fil), if you merge the raws it will consider the first file",
                        default = None
                        )

    return parser.parse_args()

if __name__ == "__main__":

    args = _get_parser()

    filenames    = args.filenames
    mode         = args.mode
    nchans       = args.nchans
    npols        = args.npols
    bd           = args.bit_depth
    ft           = args.frequency_top
    df           = args.frequency_resolution
    dt           = args.time_resolution
    telescope    = args.telescope
    project      = args.project
    source       = args.source_name
    ra_source    = args.ra_source
    dec_source   = args.dec_source
    output_names = args.output_names
    output_dir   = args.output_dir
    pol          = args.polarisation

    if mode == "Merge":
        if output_names[0] is None:
            output_names[0] = filenames[0].replace(".raw","") + ".fil"
        datawrite = []
        filepath, filename = os.path.split(filenames[0])
        rawdatafile = skarabrawfile(filename = filenames[0],
                                    filepath = filepath,
                                    projID = project,
                                    source_name = source,
                                    point_ra  = ra_source,
                                    point_dec = dec_source,
                                    telescope = telescope,
                                    channel_band_MHz = df,
                                    freq_top_MHz = ft,
                                    nchans = nchans,
                                    bit_depth = bd,
                                    tsamp_us = dt,
                                    nspectra_per_bin = npols)
        header = make_sigpyproc_header(rawdatafile)
        nbits = rawdatafile.bit_depth

        for i,filename in enumerate(filenames):
            print(f"Progress: {int((i + 1)/len(filenames)*100)} %", end='\r', flush=True)
            filepath, filename = os.path.split(filename)
            rawdatafile = skarabrawfile(filename = filename,
                                        filepath = filepath,
                                        projID = project,
                                        source_name = source,
                                        point_ra  = ra_source,
                                        point_dec = dec_source,
                                        telescope = telescope,
                                        channel_band_MHz = df,
                                        freq_top_MHz = ft,
                                        nchans = nchans,
                                        bit_depth = bd,
                                        tsamp_us = dt,
                                        nspectra_per_bin = npols)

            dynspec = rawdatafile.get_intensity_dynspec(pol = pol)

            datawrite.append(dynspec)

        datawrite = np.concatenate(datawrite, axis = 0)

        outfile = header.prepOutfile(os.path.join(output_dir,output_names[0]), back_compatible = True, nbits = nbits)

        if int(nbits) == int(8):
            datawrite = datawrite.astype("uint8")
        if int(nbits) == int(16):
            datawrite = datawrite.astype("uint16")
        if int(nbits) == int(32):
            datawrite = datawrite.astype("uint32")

        outfile.cwrite(datawrite.ravel())
        outfile.close()
        print("\nDone.")

    if mode == "Separate":

        if output_names is None:
            output_names = filenames.replace(".raw","") + ".fil"

        for filename, output_name in enumerate(zip(filenames, output_names)):
            if output_name is None:
                output_name = filename.replace(".raw","") + ".fil"
            print(f"Progress: {int((i + 1)/len(filenames)*100)} %", end='\r', flush=True)

            filepath, filename = os.path.split(filename)

            rawdatafile = skarabrawfile(filename = filename,
                                        filepath = filepath,
                                        projID = project,
                                        source_name = source,
                                        point_ra  = ra_source,
                                        point_dec = dec_source,
                                        telescope = telescope,
                                        channel_band_MHz = df,
                                        freq_top_MHz = ft,
                                        nchans = nchans,
                                        bit_depth = bd,
                                        tsamp_us = dt,
                                        nspectra_per_bin = npols)

            header = make_sigpyproc_header(rawdatafile)
            dynspec = rawdatafile.get_intensity_dynspec(pol = pol)
            nbits = rawdatafile.bit_depth

            outfile = header.prepOutfile(os.path.join(output_dir,output_name), back_compatible = True, nbits = nbits)

            if int(nbits) == int(8):
                dynspec = dynspec.astype("uint8")
            if int(nbits) == int(16):
                dynspec = dynspec.astype("uint16")
            if int(nbits) == int(32):
                dynspec = dynspec.astype("uint32")

            outfile.cwrite(dynspec.ravel())
            outfile.close()

        print("\nDone.")

"""
    for filename, output_name in filenames, output_names:

        filepath, filename = os.path.split(filename)

        rawdatafile = skarabrawfile(filename = filename,
                                    filepath = filepath,
                                    projID = project,
                                    source_name = source,
                                    point_ra  = ra_source,
                                    point_dec = dec_source,
                                    telescope = telescope,
                                    channel_band_MHz = df,
                                    freq_top_MHz = ft,
                                    nchans = nchans,
                                    bit_depth = bd,
                                    tsamp_us = dt,
                                    nspectra_per_bin = npols)

        dynspec = rawdatafile.get_intensity_dynspec(pol = pol)

        if mode == "Merge":
            header = make_sigpyproc_header(rawdatafile)
            first_file = False

            nbits = rawdatafile.bit_depth

            if output_names[0] is None:
                output_name = filenames[0].replace(".raw","") + ".fil"

            outfile = header.prepOutfile(os.path.join(output_dir,output_name), back_compatible = True, nbits = nbits)

            if int(nbits) == int(8):
                dynspec = dynspec.astype("uint8")
            if int(nbits) == int(16):
                dynspec = dynspec.astype("uint16")
            if int(nbits) == int(32):
                dynspec = dynspec.astype("uint32")

    outfile.cwrite(dynspec.ravel())
    outfile.close()


        if mode == "Merge":
            header = make_sigpyproc_header(rawdatafile)
            first_file = False

            nbits = rawdatafile.bit_depth

        if output_name is None:
            output_name = filename.replace(".raw","") + ".fil"

    outfile = header.prepOutfile(os.path.join(output_dir,output_name), back_compatible = True, nbits = nbits)

    if int(nbits) == int(8):
        dynspec = dynspec.astype("uint8")
    if int(nbits) == int(16):
        dynspec = dynspec.astype("uint16")
    if int(nbits) == int(32):
        dynspec = dynspec.astype("uint32")


    outfile.cwrite(dynspec.ravel())

    outfile.close()
"""



"""
    plt.figure()
    plt.imshow(dynspec[:,0:1000], aspect = "auto")
    plt.xlabel("Time (bins)")
    plt.ylabel("Frequency (channels)")
    plt.show()


    spectrum_xx, spectrum_yy, spectrum_xy, spectrum_yx = rawdatafile.get_spectra_per_bin()

    channels = np.arange(1, 2049)


    plt.subplot(2, 2, 1)
    plt.plot(channels, spectrum_xx)
    plt.title('Pol X')
    plt.xlabel("Channel")
    plt.ylabel("Intensity")

    plt.subplot(2, 2, 2)
    plt.plot(channels, spectrum_yy)
    plt.title('Pol Y')
    plt.xlabel("Channel")
    plt.ylabel("Intensity")

    plt.subplot(2, 2, 3)
    plt.plot(channels, spectrum_xy)
    plt.title('Pol XY')
    plt.xlabel("Channel")
    plt.ylabel("Intensity")

    plt.subplot(2, 2, 4)
    plt.plot(channels, spectrum_yx)
    plt.title('Pol YX')
    plt.xlabel("Channel")
    plt.ylabel("Intensity")

    plt.tight_layout()

    plt.show()
"""

"""
    tstamps = rawdatafile.get_all_timestamps()

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
"""

"""
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

def get_spectrum(filename, nchans = 2048, npols = 2, start = 0):

    try:
        file_abspath = os.path.abspath(filename)
        with open(file_abspath, 'rb') as file_skarab_obs:
            file_skarab_obs.seek(start)
            file_skarab_obs.read(16)
            specx = np.fromfile(file_skarab_obs, dtype="uint8", count = nchans) / 2
            specy = np.fromfile(file_skarab_obs, dtype="uint8", count = nchans) / 2

    except FileNotFoundError:
        print(f"File '{filename}' not found.")

    return specx, specy

def get_bandpasses(filename, nchans = 2048, npols = 2):

    dyspecx = []
    dyspecy = []
    total_length = os.path.getsize(filename)
    start = 0

    chunk_size = npols * nchans

    try:
        while start < total_length:
            specx, specy = get_spectrum(filename, start)

            dyspecx.append(specx)
            dyspecy.append(specy)
            start += chunk_size + 16

    except FileNotFoundError:
        print(f"File '{filename}' not found.")

    dyspecx = np.array(dyspecx)
    dyspecy = np.array(dyspecy)

    specx = np.mean(dyspecx, axis = 0)
    specy = np.mean(dyspecy, axis = 0)

    return specx, specy

"""
