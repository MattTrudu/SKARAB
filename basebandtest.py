import numpy as np
from baseband import dada
import astropy.units as u
from astropy.time import Time

# Define parameters for the DADA file
sample_rate = 16 * u.MHz  # Sampling rate
samples_per_frame = 16000  # Number of samples per frame
npol = 2  # Number of polarizations
nchan = 1  # Number of frequency channels
bps = 8  # Bits per sample
complex_data = True  # Data is complex
start_time = Time('2023-01-01T00:00:00.000')  # Start time of observation

# Calculate total number of samples
total_samples = samples_per_frame * 10  # 10 frames of data

# Create synthetic data: a simple sine wave
time = np.arange(total_samples) / sample_rate.to(u.Hz).value  # Time array in seconds
frequency = 1 * u.kHz  # Signal frequency
signal = 0.5 * np.exp(2j * np.pi * (frequency.to(u.Hz).value) * time)  # Complex sine wave

# Reshape the data to ensure it matches the expected structure
data = signal.astype(np.complex64).view(np.float32)  # Convert to float32 for writing
data = data.reshape(-1, npol)  # Ensure trailing dimension matches npol

custom_header = {'TELESCOPE': 'My_Custom_Telescope'}

# Define the output file template
output_template = '{utc_start}.{obs_offset:016d}.000000.dada'

# Open a DADA file for writing
with dada.open(output_template, 'ws', sample_rate=sample_rate,
               samples_per_frame=samples_per_frame, npol=npol,
               nchan=nchan, bps=bps, complex_data=complex_data,
               time=start_time, custom = custom_header) as fh:
    fh.write(data)

print("DADA file created successfully.")






