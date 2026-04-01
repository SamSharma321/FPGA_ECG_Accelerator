import os
import math
import numpy as np
import wfdb
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

SHOW_PLOT = True

DB_NAME = 'mitdb'
RECORD_NAME = "100"
CHANNEL_INDEX = 0
DURATION_SEC = 60 * 5 # 5 minutes
SCALE = 1000.0 # float ECG to int 16 scaling factor

# AXI / kernel assumptions
AXI_BITS = 256
SAMPLE_BITS = 16
SAMPLES_PER_WORD = AXI_BITS // SAMPLE_BITS # 16
BLOCK_SAMPLES = 2048 # 2048 samples per block

OUT_DIR = f"cpp_data/{BLOCK_SAMPLES}_samples/"
OUT_FILE = os.path.join(OUT_DIR, "ecg_dram_image.h")

os.makedirs(OUT_DIR, exist_ok=True) # Make a directory, if it exists, do nothing


# LOAD RECORDS - using the wfdb library
record = wfdb.rdrecord(RECORD_NAME, pn_dir=DB_NAME)
signal = record.p_signal[:, CHANNEL_INDEX] # Get the specified channel
fs = int(record.fs) # Sampling frequency

num_smaples = min(len(signal), DURATION_SEC * fs)
signal = signal[:num_smaples] # Truncate to desired duration

# convert to int16 - changeable according to the recorded data 
signal_i16 = np.array(signal * SCALE).astype(np.int16)

# Make it mutiple of 16 samples so that it can be taken in blocks of 2048 samples
pad_needed = (BLOCK_SAMPLES - (len(signal_i16) % BLOCK_SAMPLES)) % BLOCK_SAMPLES
if pad_needed > 0:
    signal_i16 = np.pad(signal_i16, (0, pad_needed), mode='constant') # pad sigal with zeros - only at the end

total_samples = len(signal_i16)
num_blocks = total_samples // BLOCK_SAMPLES # 2048 // 16 = 16 in ideal situation
num_axi_words = total_samples // SAMPLES_PER_WORD

packed_words = []

# prepare the data in 
for w in range(num_axi_words):
    base = w * SAMPLES_PER_WORD
    word_val = 0

    for i in range(SAMPLES_PER_WORD):
        sample = int(signal_i16[base + i]) & 0xFFFF
        word_val |= (sample << (16 * i))

    packed_words.append(word_val)

if SHOW_PLOT:
    # window length in seconds
    WINDOW_SEC = 5
    WINDOW_SAMPLES = int(WINDOW_SEC * fs)

    signal_plot = signal_i16
    total_samples = len(signal_plot)

    # time axis
    time = np.arange(total_samples) / fs

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
        valmax=(total_samples/fs) - WINDOW_SEC,
        valinit=0
    )

    def update(val):
        start_sec = slider.val
        start = int(start_sec * fs)
        end = start + WINDOW_SAMPLES

        line.set_xdata(time[start:end])
        line.set_ydata(signal_plot[start:end])

        ax.set_xlim(time[start], time[end-1])
        ax.relim()
        ax.autoscale_view(True, True, True)

        fig.canvas.draw_idle()

    slider.on_changed(update)

    plt.show()


# Create C++ header
with open(OUT_FILE, "w") as f:
    f.write("#ifndef ECG_DRAM_IMAGE_H\n")
    f.write("#define ECG_DRAM_IMAGE_H\n\n")
    f.write("#include <ap_int.h>\n")
    f.write("#include <stdint.h>\n\n")

    f.write(f"#define ECG_FS_HZ {fs}\n")
    f.write(f"#define ECG_DURATION_SEC {DURATION_SEC}\n")
    f.write(f"#define ECG_TOTAL_SAMPLES {total_samples}\n")
    f.write(f"#define ECG_BLOCK_SAMPLES {BLOCK_SAMPLES}\n")
    f.write(f"#define ECG_SAMPLES_PER_WORD {SAMPLES_PER_WORD}\n")
    f.write(f"#define ECG_NUM_BLOCKS {num_blocks}\n")
    f.write(f"#define ECG_NUM_AXI_WORDS {num_axi_words}\n")

    f.write("typedef ap_uint<256> axi_word_t;\n\n")
    f.write("static const axi_word_t ecg_dram[ECG_NUM_AXI_WORDS] = {")

    for idx, val in enumerate(packed_words):
        hex_str = f"{val:064x}"   
        sep = ", " if idx != len(packed_words) - 1 else ""
        if not(idx % 4):
            f.write("\n    ")
        if SAMPLES_PER_WORD > 1:
            f.write(rf'axi_word_t("0x{hex_str}"){sep}')
        else:
            f.write(rf'0x{hex_str}{sep}')

    f.write("};\n\n")
    f.write('#endif\n')

print("Export Complete for:")
print(f"Record            : {RECORD_NAME}")
print(f"Sampling rate     : {fs} Hz")
print(f"Samples exported  : {total_samples}")
print(f"AXI words         : {num_axi_words}")
print(f"2048-sample blocks: {num_blocks}")
print(f"Output file       : {OUT_FILE}")








