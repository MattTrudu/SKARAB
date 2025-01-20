import numpy as np
import baseband
from baseband.dada import open as dada_open

# Define parameters for the DADA file
header = {
    "OBS_ID": "test_observation",
    "SOURCE": "test_source",
    "TELESCOPE": "test_telescope",
    "INSTRUMENT": "test_instrument",
    "UTC_START": "2023-01-01-00:00:00",
    "FREQ": 1400.0,  # MHz
    "BW": 100.0,     # MHz
    "TSAMP": 1e-6,   # Time sample resolution (1 microsecond)
    "NBIT": 8,       # Number of bits per sample
    "NPOL": 2,       # Number of polarizations
    "NCHAN": 1,      # Number of frequency channels
    "SAMPLES_PER_FRAME": 1024,  # Number of samples per frame
    "FILE_SIZE": 1024 * 1024,  # Approx. file size
}

# Create synthetic data: a sine wave signal
samples_per_frame = header["SAMPLES_PER_FRAME"]
nframes = 100  # Total number of frames
total_samples = samples_per_frame * nframes
time = np.arange(total_samples) * header["TSAMP"]
signal = np.exp(2j * np.pi * 1000 * time)  # Complex exponential at 1 kHz

# Reshape the data into frames
data = signal.astype(np.complex64).reshape(-1, samples_per_frame)

# Write to a DADA file
output_file = "example.dada"
with dada_open(output_file, mode="ws", header=header) as f:
    for frame in data:
        f.write(frame.tobytes())

print(f"Saved DADA file: {output_file}")



