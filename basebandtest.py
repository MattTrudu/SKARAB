import numpy as np
import baseband
from baseband.dada import DADAFile

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
    "RESOLUTION": 16,  # Bytes per sample (complex int8 = 2 bytes per sample * 8)
    "SAMPLES_PER_FRAME": 1024,  # Number of samples per frame
    "FILE_SIZE": 1024 * 1024 * 1,  # 1 MB file size
}

# Create some synthetic data: a simple sine wave
samples_per_frame = header["SAMPLES_PER_FRAME"]
nframes = 100  # Total number of frames
time = np.arange(samples_per_frame * nframes) * header["TSAMP"]
data = (0.5 * np.sin(2 * np.pi * 1000 * time)).astype(np.float32)  # 1 kHz sine wave

# Reshape the data to match the frame structure
data = data.view(np.complex64).reshape(-1, samples_per_frame)  # Create complex samples

# Write to a DADA file
output_file = "example.dada"
with DADAFile(output_file, mode="ws", header=header) as dada:
    for frame in data:
        dada.write(frame.tobytes())  # Write each frame as bytes

print(f"Saved DADA file: {output_file}")


