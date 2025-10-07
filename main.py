import numpy as np
import soundfile as sf
import scipy.signal as signal
import matplotlib.pyplot as plt
import os

FS_TARGET = 16000

#will be multiple files, rn is just one sample
AUDIO_FILES = ["recording-24s.wav"]


# --- Helper functions ---

def ensure_mono(x):
    """Convert stereo to mono by averaging channels."""
    if x.ndim == 2:
        return np.mean(x, axis=1)
    return x


def read_and_prepare(filename, fs_target=16000):
    """Read, convert to mono, resample if needed, and return (x, fs)."""
    x, fs = sf.read(filename)
    x = ensure_mono(x)

    if fs != fs_target:
        # Use polyphase resampling for best quality
        x = signal.resample_poly(x, fs_target, fs)
        fs = fs_target
    return x, fs


def write_audio(filename, data, fs):
    """Write audio file to disk."""
    outname = os.path.splitext(filename)[0] + f"_mono_{fs}Hz.wav"
    sf.write(outname, data, fs)
    return outname


def plot_waveform(x, title, fs=None):
    """Plot waveform vs sample or time."""
    plt.figure()
    if fs:
        t = np.arange(len(x)) / fs
        plt.plot(t, x)
        plt.xlabel("Time [s]")
    else:
        plt.plot(x)
        plt.xlabel("Sample number")
    plt.title(title)
    plt.ylabel("Amplitude")
    plt.tight_layout()
    plt.show()


def generate_cosine(fs, duration, freq=1000):
    """Generate a cosine wave at freq Hz for the given duration."""
    t = np.arange(0, duration, 1/fs)
    cos_signal = np.cos(2 * np.pi * freq * t)
    return cos_signal, t


# --- Phase 2 helper functions ---
def design_filterbank(N, fs, f_low=100, f_high=8000, kind="fir", order=256):
    """Design a logarithmic-spaced FIR or IIR bandpass filterbank."""
    edges = np.logspace(np.log10(f_low), np.log10(f_high), N + 1)
    filters = []
    for i in range(N):
        low, high = edges[i], edges[i + 1]
        band = [low / (fs / 2), high / (fs / 2)]
        if kind == "fir":
            b = signal.firwin(order + 1, band, pass_zero=False)
            a = [1.0]
        else:
            b, a = signal.butter(4, band, btype='band')
        filters.append((b, a, low, high))
    return filters


def envelope_extraction(filtered, fs, cutoff=400):
    """Rectify and lowpass filter to extract envelope."""
    rectified = np.abs(filtered)
    b, a = signal.butter(4, cutoff / (fs / 2), btype='low')
    envelope = signal.filtfilt(b, a, rectified)
    return envelope


# --- MAIN SCRIPT ---

for filename in AUDIO_FILES:
    print(f"\nProcessing: {filename}")
    x, fs = read_and_prepare(filename, FS_TARGET)
    outname = write_audio(filename, x, fs)
    print(f"Saved mono/resampled audio as {outname}")

    # Plot waveform
    plot_waveform(x, f"Waveform of {filename}", fs)

    # Generate and plot cosine
    cosine, t = generate_cosine(fs, len(x) / fs, freq=1000)
    plot_waveform(cosine[:int(fs / 500)], "Two cycles of 1 kHz cosine", fs)

# Example: run filterbank and envelope extraction on first file
if AUDIO_FILES:
    x, fs = read_and_prepare(AUDIO_FILES[0], FS_TARGET)
    N = 8  # number of bands
    fb = design_filterbank(N, fs, 100, 8000, "fir", order=256)

    filtered_signals = []
    envelopes = []

    for (b, a, f1, f2) in fb:
        y = signal.lfilter(b, a, x)
        env = envelope_extraction(y, fs)
        filtered_signals.append(y)
        envelopes.append(env)

    # Plot lowest and highest band results
    t = np.arange(len(x)) / fs

    plt.figure(figsize=(10, 6))
    plt.subplot(2, 1, 1)
    plt.plot(t, filtered_signals[0])
    plt.title("Filtered Signal - Lowest Band")
    plt.xlabel("Time [s]")

    plt.subplot(2, 1, 2)
    plt.plot(t, filtered_signals[-1])
    plt.title("Filtered Signal - Highest Band")
    plt.xlabel("Time [s]")
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.subplot(2, 1, 1)
    plt.plot(t, envelopes[0])
    plt.title("Envelope - Lowest Band")
    plt.xlabel("Time [s]")

    plt.subplot(2, 1, 2)
    plt.plot(t, envelopes[-1])
    plt.title("Envelope - Highest Band")
    plt.xlabel("Time [s]")
    plt.tight_layout()
    plt.show()
