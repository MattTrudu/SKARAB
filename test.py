# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import datetime



def read_spectrum(file_path, num_spectra, samples_per_packet=1024):

    bufferp_right_im = []
    bufferp_right_re = []
    bufferp_left_im = []
    bufferp_left_re = []

    count = 0

    try:
        with open(file_path, 'rb') as file:
            while count < num_spectra:
                # Leggi il timestamp (8 byte per `tv_sec` e 8 byte per `tv_usec`)
                timestamp_data = np.fromfile(file, dtype=np.int64, count=2)
                
                if len(timestamp_data) != 2:
                    print("Error: incomplete timestamp data.")
                    break
                
                timestamp_sec = timestamp_data[0]
                timestamp_usec = timestamp_data[1]

                dt = datetime.datetime.utcfromtimestamp(timestamp_sec)
                dt = dt.replace(microsecond=timestamp_usec)
                readable_timestamp = dt.strftime("%Y-%m-%d %H:%M:%S.%f")
                print("Timestamp:", readable_timestamp)

                # Lettura degli spettri come array numpy di int16
                right_im = np.fromfile(file, dtype=np.int16, count=samples_per_packet)
                right_re = np.fromfile(file, dtype=np.int16, count=samples_per_packet)
                left_im = np.fromfile(file, dtype=np.int16, count=samples_per_packet)
                left_re = np.fromfile(file, dtype=np.int16, count=samples_per_packet)

                if (len(right_im) < samples_per_packet or len(right_re) < samples_per_packet or
                    len(left_im) < samples_per_packet or len(left_re) < samples_per_packet):
                    break

                bufferp_right_im.extend(right_im)
                bufferp_right_re.extend(right_re)
                bufferp_left_im.extend(left_im)
                bufferp_left_re.extend(left_re)

                count += 1
        
        if count > 0:
            plot(bufferp_right_im, bufferp_right_re, bufferp_left_im, bufferp_left_re)
        else:
            print("No spectra were read.")

    except IOError:
        print("File not found: {}".format(file_path))
    except Exception as e:
        print("Error reading file: {}".format(e))



def plot(right_im, right_re, left_im, left_re):

    num_points = len(right_im)
    x = range(num_points)

    plt.figure(figsize=(10, 6))

    plt.subplot(2, 2, 1)
    plt.plot(x, right_im)
    plt.title('Right Imaginary')
    plt.xlabel("Samples")
    plt.ylabel("Intensity")

    plt.subplot(2, 2, 2)
    plt.plot(x, right_re)
    plt.title('Right Real')
    plt.xlabel("Samples")
    plt.ylabel("Intensity")

    plt.subplot(2, 2, 3)
    plt.plot(x, left_im)
    plt.title('Left Imaginary')
    plt.xlabel("Samples")
    plt.ylabel("Intensity")

    plt.subplot(2, 2, 4)
    plt.plot(x, left_re)
    plt.title('Left Real')
    plt.xlabel("Samples")
    plt.ylabel("Intensity")

    plt.tight_layout()
    plt.savefig('skarab.png')


if __name__ == "__main__":
    if len(sys.argv) < 1:
        print("Usage: python script.py <file_name>")
        sys.exit(1)

    base_path = r"C:\Users\Alessandro\Downloads"
    file_name = sys.argv[1]
    file_path = os.path.join(base_path, file_name)

    read_spectrum(file_path, 10)
