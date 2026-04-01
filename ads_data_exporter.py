import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

SHOW_PLOT = True

CSV_FILE = "2.csv"
CHANNEL_NAME = "EEF3"          # choose electrode column here
FS_HZ = 500                    # set this to your actual sampling rate
DURATION_SEC = 60 * 5          # 5 minutes
SCALE = 1.0                    # use 1.0 if CSV is already integer-like / ADC-like

# AXI / kernel assumptions
AXI_BITS = 256
SAMPLE_BITS = 16
SAMPLES_PER_WORD = AXI_BITS // SAMPLE_BITS   # 16
BLOCK_SAMPLES = 2048

OUT_DIR = f"ads_cpp_data/{BLOCK_SAMPLES}_samples/"
OUT_FILE = os.path.join(OUT_DIR, "ecg_dram_image.h")

os.makedirs(OUT_DIR, exist_ok=True)

# -----------------------------
# Load CSV
# -----------------------------
df = pd.read_csv(CSV_FILE)

print("Columns found:", list(df.columns))

if CHANNEL_NAME not in df.columns:
    raise ValueError(f"Column '{CHANNEL_NAME}' not found in CSV")

# Extract one electrode column
signal = df[CHANNEL_NAME].dropna().to_numpy()

# Remove NaNs if present
signal = signal[~np.isnan(signal)]

# Truncate to desired duration
num_samples = min(len(signal), DURATION_SEC * FS_HZ)
signal = signal[:num_samples]

# Convert to int16
# If your CSV is already in ADC counts and safely inside int16 range,
# SCALE can stay 1.0
signal_i16 = np.round(signal * SCALE).astype(np.int16)

# Pad so total samples is a multiple of BLOCK_SAMPLES
pad_needed = (BLOCK_SAMPLES - (len(signal_i16) % BLOCK_SAMPLES)) % BLOCK_SAMPLES
if pad_needed > 0:
    signal_i16 = np.pad(signal_i16, (0, pad_needed), mode="constant")

total_samples = len(signal_i16)
num_blocks = total_samples // BLOCK_SAMPLES
num_axi_words = total_samples // SAMPLES_PER_WORD

# -----------------------------
# Pack 16 x int16 into one 256-bit word
# -----------------------------
packed_words = []

for w in range(num_axi_words):
    base = w * SAMPLES_PER_WORD
    word_val = 0

    for i in range(SAMPLES_PER_WORD):
        sample = int(signal_i16[base + i]) & 0xFFFF
        word_val |= sample << (SAMPLE_BITS * i)

    packed_words.append(word_val)


if SHOW_PLOT:
    # window length in seconds
    WINDOW_SEC = 5
    WINDOW_SAMPLES = int(WINDOW_SEC * FS_HZ)

    signal_plot = signal_i16
    total_samples = len(signal_plot)

    # time axis
    time = np.arange(total_samples) / FS_HZ

    # initial window
    start = 0
    end = start + WINDOW_SAMPLES

    fig, ax = plt.subplots(figsize=(12,4))
    plt.subplots_adjust(bottom=0.25)

    line, = ax.plot(time[start:end], signal_plot[start:end])
    ax.set_xlabel("Time (seconds)")
    ax.set_ylabel("Amplitude")
    ax.set_title("Scrollable ECG Viewer")
    ax.grid(True)

    # slider axis
    ax_slider = plt.axes([0.1, 0.1, 0.8, 0.03])
    slider = Slider(
        ax=ax_slider,
        label="Scroll (seconds)",
        valmin=0,
        valmax=(total_samples/FS_HZ) - WINDOW_SEC,
        valinit=0
    )

    def update(val):
        start_sec = slider.val
        start = int(start_sec * FS_HZ)
        end = start + WINDOW_SAMPLES

        line.set_xdata(time[start:end])
        line.set_ydata(signal_plot[start:end])

        ax.set_xlim(time[start], time[end-1])
        ax.relim()
        ax.autoscale_view(True, True, True)

        fig.canvas.draw_idle()

    slider.on_changed(update)

    plt.show()



# -----------------------------
# Create C++ header
# -----------------------------
with open(OUT_FILE, "w") as f:
    f.write("#ifndef ECG_DRAM_IMAGE_H\n")
    f.write("#define ECG_DRAM_IMAGE_H\n\n")
    f.write("#include <ap_int.h>\n")
    f.write("#include <stdint.h>\n\n")

    f.write(f"#define ECG_FS_HZ {FS_HZ}\n")
    f.write(f"#define ECG_DURATION_SEC {DURATION_SEC}\n")
    f.write(f"#define ECG_TOTAL_SAMPLES {total_samples}\n")
    f.write(f"#define ECG_BLOCK_SAMPLES {BLOCK_SAMPLES}\n")
    f.write(f"#define ECG_SAMPLES_PER_WORD {SAMPLES_PER_WORD}\n")
    f.write(f"#define ECG_NUM_BLOCKS {num_blocks}\n")
    f.write(f"#define ECG_NUM_AXI_WORDS {num_axi_words}\n\n")

    f.write(f"typedef ap_uint<{AXI_BITS}> axi_word_t;\n\n")
    f.write("static const axi_word_t ecg_dram[ECG_NUM_AXI_WORDS] = {")

    for idx, val in enumerate(packed_words):
        hex_str = f"{val:064x}"
        sep = "," if idx != len(packed_words) - 1 else ""
        if idx % 2 == 0:
            f.write("\n    ")
        f.write(f'axi_word_t("0x{hex_str}"){sep}')

    f.write("\n};\n\n")
    f.write("#endif\n")

print("Export Complete for:")
print(f"CSV file          : {CSV_FILE}")
print(f"Channel           : {CHANNEL_NAME}")
print(f"Sampling rate     : {FS_HZ} Hz")
print(f"Samples exported  : {total_samples}")
print(f"AXI words         : {num_axi_words}")
print(f"2048-sample blocks: {num_blocks}")
print(f"Output file       : {OUT_FILE}")