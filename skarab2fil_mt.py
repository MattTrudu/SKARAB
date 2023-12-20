#!/usr/bin/env python

#################### ALESSANDRO RIDOLFI ########################
#                        Version 0.1                           #
#                    Cagliari, Feb 2020                        #
################################################################


'''
Formato dati SARDARA al 20 Febbraio 2020

Ogni tot tempo, uno spettro viene scritto, ed e' composto cosi':

[  16 Bytes  ] [         4096 Bytes         ]   [         4096 Bytes         ]
[ Time stamp ] [        Spettro pol X       ]   [        Spettro pol Y       ]

1) TIME STAMP
   I 16 bytes del time stamp sono divisi in due blocchi da 8 bytes, entrambi  unsigned long long int ('Q' in struck.unpack)

   Primi   8 bytes = Numero di secondi dal 1970 (GMtime)
   Secondi 8 bytes = Numero di microsecondi aggiuntivi


2) SPETTRI
   Ogni spettro e' composto da 4096 canali, ognuno rappresentato da un unsigned int da 1 Byte (= 8 bit = 256 livelli, unsigned char o 'B' in struck.unpack)


   Nei dati del 15 Febbraio 2020, i dati sono stati campionati ad 1200 MHz, SENZA downconversione.
   Tuttavia la banda di osservazione del ricevitore in banda L va da 1300 a 1800 MHz, che quindi cade nella seconda finestra di Nyquist.
   Il che significa che la vera banda sara' visible nella prima finestra di Nyquist, ma flippata.

   Dunque, i 4096 canali avranno bordi a 2400 MHz per il primo canale e 1200 MHz per l'ultimo canale.
   La larghezza dei canali sara' dunque 1200 MHz / 4096 = 0.29296875 MHz

   La frequenza centrale del top della banda sara' dunque 1200 -0.5*0.29296875 MHz = 1199.853515625 MHz


'''

import sys, os, os.path, glob, shutil, time, struct, datetime
import numpy as np
from astropy.time import Time
from sigpyproc import Header, Filterbank
import pylab as plt


def get_start_time_astropy(file_name):
        file_abspath = os.path.abspath(file_name)
        file_sardara_obs = open(file_abspath, 'rb')

        byte_seconds_from_1970_spectrum_0    = file_sardara_obs.read(8)
        file_sardara_obs.close()

        float_seconds_from_1970 = get_time_stamp_from_spectrum(byte_seconds_from_1970_spectrum_0, byte_microseconds_spectrum_0)


        T_start = Time(float_seconds_from_1970, format='unix')   #unix format is number of seconds from 1970-01-01

        return T_start



def get_time_stamp_from_spectrum(byte_seconds_from_1970, byte_microseconds):
        str_seconds_from_1970 = struct.unpack('Q', byte_seconds_from_1970)[0]
        str_byte_microseconds = struct.unpack('Q', byte_microseconds     )[0]


        int_seconds_from_1970    = int(str_seconds_from_1970)
        int_microseconds         = int(str_byte_microseconds)

        float_seconds_from_1970  = float64(int_seconds_from_1970)
        float_milliseconds       = float64(int_microseconds) / 1000.
        float_fractional_second  = float_milliseconds / 1000.

        float_seconds_from_1970 = float_seconds_from_1970 + float_fractional_second

        #print "get_time_stamp_from_spectrum::  int_seconds_from_1970 = ", int_seconds_from_1970
        #print "get_time_stamp_from_spectrum::  int_microseconds      = ", int_microseconds
        #print "get_time_stamp_from_spectrum::  float_seconds_from_1970 = ", float_seconds_from_1970

        return float_seconds_from_1970


class Sardara_Observation(object):
        def __init__(self, file_name, nchan_per_spectrum, npol, freq_low_edge_first_chan_rawfile_MHz, bw_total_MHz, t_samp_s, flag_descending_freqs, source_name, ra_str, dec_str, flag_check_packet_loss=1):
                self.file_abspath = os.path.abspath(file_name)
                self.file_nameonly = self.file_abspath.split("/")[-1]
                self.file_basename, self.file_extension = os.path.splitext(self.file_nameonly)

                self.T_start_astropy_Time     = get_start_time_astropy(self.file_abspath)
                self.MJD_start                = self.T_start_astropy_Time.mjd
                self.MJD_int                  = int(self.MJD_start)
                self.source_name              = source_name
                self.ra_str                   = ra_str
                self.dec_str                  = dec_str


                self.telescope                = "SRT"
                self.backend                  = "SKARAB"
                self.date_obs                 = datetime.datetime.fromtimestamp(self.T_start_astropy_Time.unix).strftime("%d/%m/%Y")      #30/10/2019
                self.nbits                    = 8
                self.npol                     = npol
                self.nchan                    = nchan_per_spectrum

                if flag_descending_freqs == 1:
                        self.band_orientation_factor = -1
                        print("WARNING: You set flag_descending_freqs = 1")
                elif flag_descending_freqs == 0:
                        self.band_orientation_factor = +1


                self.bw_total_MHz_fabs        = bw_total_MHz
                self.bw_total_MHz_with_sign   = self.band_orientation_factor * self.bw_total_MHz_fabs

                self.chanbw_MHz_fabs          = np.fabs(self.bw_total_MHz_fabs/ self.nchan)   #Double fabs, redundant
                self.chanbw_MHz_with_sign     = self.band_orientation_factor * self.chanbw_MHz_fabs


                self.freq_low_edge_first_chan_rawfile_MHz = freq_low_edge_first_chan_rawfile_MHz
                if flag_descending_freqs == 1:
                        self.freq_observing_band_top_edge_MHz     = freq_low_edge_first_chan_rawfile_MHz
                        self.freq_observing_band_low_edge_MHz     = self.freq_observing_band_top_edge_MHz - self.bw_total_MHz_fabs
                elif flag_descending_freqs == 0:
                        self.freq_observing_band_low_edge_MHz     = freq_low_edge_first_chan_rawfile_MHz
                        self.freq_observing_band_top_edge_MHz     = self.freq_observing_band_low_edge_MHz + self.bw_total_MHz_fabs




                self.freq_central_MHz         = self.freq_observing_band_low_edge_MHz + (0.5*self.bw_total_MHz_fabs)
                self.freq_filterbank_fch1_MHz = self.freq_observing_band_top_edge_MHz - (0.5*self.chanbw_MHz_fabs)


                print("This way I interpret:")
                print("Top edge of actual observing band: %.3f MHz" % (self.freq_observing_band_top_edge_MHz))
                print("Central frequency:                 %.3f MHz" % (self.freq_central_MHz))
                print("Low edge of actual observing band: %.3f MHz" % (self.freq_observing_band_low_edge_MHz))

                print("First channel low edge raw_files:  %.3f MHz" % (self.freq_low_edge_first_chan_rawfile_MHz))
                print("First channel filterbank (fch1):   %.3f MHz" % (self.freq_filterbank_fch1_MHz))


                self.numbytes                 = os.path.getsize(self.file_abspath)
                self.N_samples                = self.numbytes / (16 + self.npol*self.nchan)

                if t_samp_s == -1:
                        self.t_samp_s                 = self.get_sampling_time(self.nchan, self.npol)
                        print("NOW self.t_samp_s = ", self.t_samp_s )
                else:
                        self.t_samp_s                 = t_samp_s
                        print("t_samp_s != -1 ---> set from shell ---> self.t_samp_s = %.20f" % (self.t_samp_s))


                        #print "Sampling frequency = %.5f Hz" % (1./self.t_samp_s)

                #self.coordinate_sys           = self.dict_properties['Coordinate Sys']
                #self.drift_rate               = self.dict_properties['Drift Rate']


                self.T_obs_s                  = self.N_samples*self.t_samp_s


                self.header_filterbank        = self.make_filterbank_header()






        def get_sampling_time(self, nchan, npol):
                file_sardara_obs = open(self.file_abspath, 'rb')

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



                #byte_seconds_from_1970_spectrum_0    = file_sardara_obs.read(8)
                #byte_microseconds_spectrum_0         = file_sardara_obs.read(8)

                #file_sardara_obs.seek(npol*nchan + 16)

                #byte_seconds_from_1970_spectrum_1    = file_sardara_obs.read(8)
                #byte_microseconds_spectrum_1         = file_sardara_obs.read(8)


                #print "byte_seconds_from_1970_spectrum_0 = ", byte_seconds_from_1970_spectrum_0
                #time_stamp_spectrum_1_seconds = get_time_stamp_from_spectrum(byte_seconds_from_1970_spectrum_1, byte_microseconds_spectrum_1)

                #print "time_stamp_spectrum_0_seconds = ", time_stamp_spectrum_0_seconds
                #print "time_stamp_spectrum_1_seconds = ", time_stamp_spectrum_1_seconds

                #t_samp_s = time_stamp_spectrum_1_seconds - time_stamp_spectrum_0_seconds

                #date_time_UTC_spectrum_0 = datetime.datetime.fromtimestamp(time_stamp_spectrum_0_seconds)
                #date_time_UTC_spectrum_1 = datetime.datetime.fromtimestamp(time_stamp_spectrum_1_seconds)

                #print "date_time_UTC_spectrum_0 = ", date_time_UTC_spectrum_0
                #print "date_time_UTC_spectrum_1 = ", date_time_UTC_spectrum_1


                return t_samp_s

        def make_filterbank_header(self):
                filterbank_header = Header.Header(  { "telescope_id":     10,
                                           "machine_id":       0,
                                           "data_type":        1,
                                           "source_name":      self.source_name,
                                           "src_raj":          float64(self.ra_str.replace(":", ""))   ,
                                           "src_dej":          float64(self.dec_str.replace(":", ""))   ,
                                           "tstart":           self.MJD_start  ,
                                           "tsamp":            self.t_samp_s,
                                           "nsamples":         self.N_samples,
                                           "nifs":             1         ,
                                           "fch1":             self.freq_filterbank_fch1_MHz       ,
                                           "foff":             -np.fabs(self.chanbw_MHz_fabs),
                                           "fchannel":         self.freq_central_MHz,
                                           "nchans":           self.nchan            ,
                                           "obs_date":         self.date_obs           }
                )




                return filterbank_header





#DEFAULT PARAMETERS
string_version = "0.1 (26Feb2019)"
flag_flipband = 0
verbosity_level = 1
N_spectra_per_block = 100000
outfilename = ""
t_samp_s = -1
nchan   = 4096
npol    = 2
flag_descending_freqs = 0

#SHELL ARGUMENT
if (len(sys.argv) == 1 or ("-h" in sys.argv) or ("-help" in sys.argv) or ("--help" in sys.argv)):

        print("USAGE: %s -rawfiles \"sardara*.raw\" -source_name \"J2145-0750\" -nchan_per_spectrum 4096 -npol {1,2} -freq_low_edge_first_chan_rawfile_MHz [2400] -bw_total_MHz [MHz] [-descending_freqs] [-o outfilename] [-Q]" % (os.path.basename(sys.argv[0])))

        print("%20s  %-33s:  %-50s" % ("-rawfiles", "<files*.raw>", "Original SARDARA observation files"))
        print("%20s  %-33s:  %-50s" % ("-nchan_per_spectrum",         "N (default: 4096)", "Input file is Total Intensity (onlyI) or Full Stokes (FS) data"))
        print("%20s  %-33s:  %-50s" % ("-npol", "{1,2} (default: 2)", "Number of polarization channels"))
        print("%20s  %-33s:  %-50s" % ("-freq_low_edge_first_chan_rawfile_MHz", "[float]", "Frequency of the edge of the first frequency channel of raw files"))
        print("%20s  %-33s:  %-50s" % ("-bw_total_MHz", "[MHz]", "Total observing bandwidth"))
        print("%20s  %-33s:  %-50s" % ("-N_spectra_per_block", "[int] (default: 100000)", "Number of spectra per block to write"))
        print("%20s  %-33s:  %-50s" % ("-descending_freqs", "(default: NO)", "Channels are in descending order (highest frequency is first channel)"))
        print("%20s  %-33s:  %-50s" % ("-source_name", "", "Name of the source"))
        print("%20s  %-33s:  %-50s" % ("-ra_dec", "\"hh:mm:ss.xx,deg:arcm:arcs.xx\"", "Coordinates of the source "))
        print("%20s  %-33s:  %-50s" % ("-o", "(default: 'basename.fil')", "Output filename"))
        print("%20s  %-33s:  %-50s" % ("-Q", "", "Quiet mode: do not print any information"))
        print("%20s  %-33s:  %-50s" % ("-version", "", "Print code version"))

        exit()
elif (("-v" in sys.argv) or ("-version" in sys.argv) or ("--version" in sys.argv)):
        print("Version: %s" % (string_version))
        exit()
else:
        for n in range(1, len(sys.argv)):
                if (sys.argv[n] == "-rawfiles"):
                        string_files = sys.argv[n+1]
                        if ("*" in string_files) or ("?" in string_files):
                                list_files = [ os.path.abspath(x) for x in sorted(glob.glob(string_files.strip("\"")))]
                        else:
                                list_files = [ os.path.abspath(x) for x in string_files.replace(" ", "").split(",")]



                elif (sys.argv[n] == "-nchan_per_spectrum"):
                        nchan = int(sys.argv[n+1])
                elif (sys.argv[n] == "-freq_low_edge_first_chan_rawfile_MHz"):
                        freq_low_edge_first_chan_rawfile_MHz = float(sys.argv[n+1])
                elif (sys.argv[n] == "-bw_total_MHz"):
                        bw_total_MHz = float(sys.argv[n+1])
                elif (sys.argv[n] == "-t_samp_s"):
                        t_samp_s = float(sys.argv[n+1])

                elif (sys.argv[n] == "-npol"):
                        npol = int(sys.argv[n+1])
                        if npol != 1 and npol != 2:
                                print("ERROR! npol must be either 1 or 2!")
                                exit()
                elif (sys.argv[n] == "-descending_freqs"):
                        flag_descending_freqs = 1
                elif (sys.argv[n] == "-N_spectra_per_block"):
                        N_spectra_per_block = int(sys.argv[n+1])
                elif (sys.argv[n] == "-source_name"):
                        source_name = sys.argv[n+1]
                elif (sys.argv[n] == "-ra_dec"):
                        ra_str, dec_str = sys.argv[n+1].split(",")
                elif (sys.argv[n] == "-o"):
                        outfilename = sys.argv[n+1]
                elif (sys.argv[n] == "-Q"):
                        verbosity_level = 0
N_files = len(list_files)
###################################################
# LOOP OVER RAW DATA FILES
###################################################


first_sardara_observation_filename = list_files[0]
first_sardara_observation_abspath = os.path.abspath(first_sardara_observation_filename)

first_sardara_observation          = Sardara_Observation(first_sardara_observation_abspath, nchan, npol, freq_low_edge_first_chan_rawfile_MHz, bw_total_MHz, t_samp_s, flag_descending_freqs, source_name, ra_str, dec_str, flag_check_packet_loss=1)



if verbosity_level >= 1:

        print("#"*62)
        print("#" + " "*26 + "%s" % ("%s") % (os.path.basename(sys.argv[0])) + " "*25 + "#")
        print("#" + " "*23 + "%s" % (string_version) + " "*22 + "#")
        print("#"*62)

        print("======================================================================")
        print("Observation Properties:")

        print("%40s: %s" % ("Source name", first_sardara_observation.source_name))
        print("%40s: %s" % ("RA", first_sardara_observation.ra_str))
        print("%40s: %s" % ("DEC", first_sardara_observation.dec_str))

        print("%40s: %s" % ("Num of channels", first_sardara_observation.nchan))
        print("%40s: %s" % ("Channel width (MHz)", first_sardara_observation.chanbw_MHz_fabs))
        print("%40s: %s" % ("Lower edge of first channel of each spectrum (MHz)", first_sardara_observation.freq_low_edge_first_chan_rawfile_MHz))
        print("%40s: %s" % ("Central freq. of observing band (MHz)", first_sardara_observation.freq_central_MHz))
        print("%40s: %s" % ("Total observing bandwidth (MHz)", first_sardara_observation.bw_total_MHz_fabs))
        print("%40s: %s" % ("Sampling time (s)", first_sardara_observation.t_samp_s))
        print("%40s: %s" % ("Bits per sample", first_sardara_observation.nbits))
        print("%40s: %s" % ("MJD of first time sample", first_sardara_observation.MJD_start))
        print("%40s: %s" % ("Number of samples", first_sardara_observation.N_samples))
        print("%40s: %s" % ("Observation length (s)", first_sardara_observation.T_obs_s))




if outfilename == "":
        sardara_observation_filterbank_filename = os.path.splitext(first_sardara_observation_filename)[0] + ".fil"
else:
        sardara_observation_filterbank_filename = outfilename

sardara_observation_filterbank = first_sardara_observation.header_filterbank.prepOutfile( sardara_observation_filterbank_filename, nbits=8)
sardara_observation_filterbank.cwrite(np.array([]))
sardara_observation_filterbank.flush()
sardara_observation_filterbank.close()

sardara_observation_filterbank_file = open(sardara_observation_filterbank_filename, 'ab')

for k in range(N_files):
        raw_filename = list_files[k]
        raw_filename_abspath = os.path.abspath(raw_filename)

        sardara_observation = Sardara_Observation(raw_filename_abspath, nchan, npol, freq_low_edge_first_chan_rawfile_MHz, bw_total_MHz, t_samp_s, flag_descending_freqs, source_name, ra_str, dec_str, flag_check_packet_loss=1)


        fil_header_size = os.path.getsize(sardara_observation_filterbank_filename)
        N_tot = os.path.getsize(raw_filename_abspath)

        N_spectra = sardara_observation.N_samples
        nchan     = sardara_observation.nchan
        npol      = sardara_observation.npol





        with open(raw_filename_abspath, "r") as f:
                for i in range(N_spectra):
                        f.seek( (16 + npol*nchan)*i )
                        f.read(16)
                        pol_x = np.fromfile(f, dtype="uint8", count=nchan) / 2
                        pol_y = np.fromfile(f, dtype="uint8", count=nchan) / 2

                        total_I = pol_x #+ pol_y
                        #plt.plot(total_I)
                        #plt.show()
                        sardara_observation_filterbank_file.write(total_I)
                        sys.stdout.write('\rFile %d / %d   ---- > Now I am at spectrum # %d / %d --> f.tell() = %d' % (k, N_files, i, N_spectra, f.tell()))
                        sys.stdout.flush()

        sardara_observation_filterbank_file.flush()

sardara_observation_filterbank_file.close()
