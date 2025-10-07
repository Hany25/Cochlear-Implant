import numpy as np
import soundfile as sf
import subprocess
from scipy.signal import resample
import os
import matplotlib.pyplot as plt

FS_TARGET = 16000  # target sample rate (Hz)

def ensure_pcm_wav(filename):
    """Ensure the input file is a readable PCM WAV. Convert if needed."""
    try:
        info = sf.info(filename)
        if info.format != "WAV" or "PCM" not in info.subtype:
            raise ValueError("Not PCM WAV")
    except Exception:
        # Use FFmpeg to convert it
        base, _ = os.path.splitext(filename)
        fixed_file = base + "_fixed.wav"
        print(f"Converting {filename} → {fixed_file}")
        subprocess.run([
            "ffmpeg", "-y", "-i", filename,
            "-acodec", "pcm_s16le", "-ar", "44100", "-ac", "1", fixed_file
        ], check=True)
        filename = fixed_file
    return filename

def read_and_prepare(filename, fs_target=FS_TARGET):
    """Read, convert to mono if needed, normalize, and resample."""
    filename = ensure_pcm_wav(filename)
    x, fs = sf.read(filename)

    # If stereo → mono
    if x.ndim > 1:
        print("Stereo detected, converting to mono.")
        x = np.mean(x, axis=1)

    # Normalize
    x = x / np.max(np.abs(x))

    # Resample to 16 kHz
    if fs != fs_target:
        print(f"Resampling from {fs} Hz → {fs_target} Hz")
        num_samples = int(len(x) * fs_target / fs)
        x = resample(x, num_samples)
        fs = fs_target

    return x, fs

def plot_waveform(x, fs, title):
    """Plot waveform vs. sample number."""
    plt.figure()
    plt.plot(np.arange(len(x)), x)
    plt.title(title)
    plt.xlabel("Sample number")
    plt.ylabel("Amplitude")
    plt.show()

def generate_cosine(fs, duration, freq=1000):
    """Generate 1 kHz cosine with same duration as input."""
    t = np.arange(int(fs * duration)) / fs
    signal = np.cos(2 * np.pi * freq * t)
    return t, signal

# ------------------ MAIN ------------------
if __name__ == "__main__":
    filename = "recording-11s.wav"
    print(f"Processing: {filename}")

    x, fs = read_and_prepare(filename, FS_TARGET)
    print(f"Loaded file → {len(x)} samples at {fs} Hz")

    # Play sound if desired (requires sounddevice)
    # import sounddevice as sd; sd.play(x, fs); sd.wait()

    # Write normalized + resampled output
    sf.write("processed.wav", x, fs)
    print("Wrote processed.wav")

    # Plot input waveform
    plot_waveform(x, fs, "Processed Input Waveform")

    # Generate 1 kHz cosine and plot two cycles
    duration = len(x) / fs
    t, cos_signal = generate_cosine(fs, duration)
    plt.figure()
    plt.plot(t[:int(fs / 500)], cos_signal[:int(fs / 500)])  # two cycles of 1 kHz
    plt.title("1 kHz Cosine (two cycles)")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.show()
