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
        #datawrite = []
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

        outfile = header.prepOutfile(os.path.join(output_dir,output_names[0]), back_compatible = True, nbits = nbits)

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

            #datawrite.append(dynspec)

        #datawrite = np.concatenate(datawrite, axis = 0)



            if int(nbits) == int(8):
                dynspec = dynspec.astype("uint8")
            if int(nbits) == int(16):
                dynspec = dynspec.astype("uint16")
            if int(nbits) == int(32):
                dynspec = dynspec.astype("uint32")

            outfile.cwrite(dynspec.ravel())

        outfile.close()
        print("\nDone.")

    if mode == "Separate":

        if output_names is None:
            output_names = [None] * len(filenames)


        for i,filename in enumerate(filenames):
            if output_names[i] is None:
                output_names[i] = filename.replace(".raw","") + ".fil"
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

            outfile = header.prepOutfile(os.path.join(output_dir,output_names[i]), back_compatible = True, nbits = nbits)

            if int(nbits) == int(8):
                dynspec = dynspec.astype("uint8")
            if int(nbits) == int(16):
                dynspec = dynspec.astype("uint16")
            if int(nbits) == int(32):
                dynspec = dynspec.astype("uint32")

            outfile.cwrite(dynspec.ravel())
            outfile.close()

        print("\nDone.")


