### FPGA ECG Accelerator

A high-performance hardware accelerator design for real-time ECG signal processing on FPGA using Xilinx Vitis HLS. This design processes ECG data blocks to detect heartbeat peaks and compute heart rate statistics.

#### Overview

This repository contains a fully optimized FPGA kernel tha- Careful buffer sizing to prevent backpressure in FIFO chains

##### Circular Buffer Patternrocesses ECG (electrocardiogram) data streams to:
- Remove baseline wandering with high-pass filtering
- Compute signal derivatives
- Detect R-wave peaks using energy integration
- Analyze heart rate variability (RR intervals)
- Generate cardiac arrhythmia indicators

The design achieves **strict II=1 (initiation interval) on all pipeline stages**, enabling maximum throughput with minimal resource utilization.

#### Repository Structure

##### Core Design Files

#### **`top.cpp`** – Main FPGA Kernel
The optimized Vitis HLS kernel implementing a 9-stage streaming pipeline architecture.

**Key Components:**
- **Stage 1 (k_load)**: Burst-reads AXI words from DRAM into a BRAM register file
- **Stage 2 (k_unpack)**: Unpacks 32 samples per AXI word in parallel (II=FACTOR)
- **Stage 3 (k_baseline)**: Causal moving-average high-pass filter (128-sample window)
- **Stage 4 (k_deriv)**: First-order backward difference filter
- **Stage 5 (k_abs_vec)**: Vectorized absolute-value computation (16-lane parallel processing)
- **Stage 6 (k_integrate)**: Moving-window integrator (32-sample window)
- **Stage 7 (k_detect_peaks)**: Threshold crossing detector with refractory period gating
- **Stage 8 (k_rr_analysis)**: RR interval statistical analysis (O(peaks) complexity)
- **Stage 9 (k_finalize_summary)**: Generates final output metrics

**Architecture Diagram:**
```
  ┌────────────┐   s_wide    ┌────────────┐   s_ecg (ecg_t, 1/cyc)
  │  k_load    │────────────▶│ k_unpack   │──────────────────────────────┐
  │ (burst AXI)│             │ (FACTOR    │                              │
  └────────────┘             │  parallel) │                              ▼
                             └────────────┘               ┌─────────────────────┐
                                                           │   k_baseline        │
                                                           │ (circular buf II=1) │
                                                           └──────────┬──────────┘
                                                                      │ s_hp
                                                                      ▼
                                                           ┌─────────────────────┐
                                                           │   k_deriv           │
                                                           │   (diff,  II=1)     │
                                                           └──────────┬──────────┘
                                                                      │ s_deriv
                                                         ┌────────────▼───────────┐
                                                         │  k_abs_vec             │
                                                         │  (FACTOR lanes parallel│
                                                         │   packed wide stream)  │
                                                         └────────────┬───────────┘
                                                                      │ s_energy_wide
                                                         ┌────────────▼───────────┐
                                                         │  k_integrate           │
                                                         │  (circular buf II=1)   │
                                                         └────────────┬───────────┘
                                                                      │ s_integrated
                                                         ┌────────────▼───────────┐
                                                         │  k_detect_peaks        │
                                                         │  (refractory gate)     │
                                                         └────────────┬───────────┘
                                                                      │ peak_positions[]
                                                                      │ peak_count
                                                                      ▼
                                              ┌──────────────────────────────────────┐
                                              │  k_rr_analysis  (O(peaks), unrolled) │
                                              └──────────────────┬───────────────────┘
                                                                 │
                                                                 ▼
                                              ┌──────────────────────────────────────┐
                                              │  k_finalize_summary                  │
                                              └──────────────────────────────────────┘
```

###### **`dcl.h`** – Type Definitions & Configuration
Defines all data types, constants, and kernel interface.

**Type Definitions:**
- `ecg_t`: signed 16-bit sample values
- `deriv_t`: signed 16-bit derivative values
- `mag_t`: unsigned 16-bit magnitude values
- `accum_t`: unsigned 32-bit accumulator
- `rr_t`: unsigned 16-bit RR interval
- `bpm_t`: unsigned 16-bit BPM value
- `axi_word_t`: 512-bit AXI data word (32 × 16-bit samples)
- `ECGSummaryOut`: Output struct containing:
  - `num_peaks_detected`: Number of detected heartbeat peaks
  - `avg_rr_samples`: Average RR interval in samples
  - `avg_rr_ms`: Average RR interval in milliseconds
  - `avg_bpm`: Average beats per minute
  - `low_hr`: Flag (low heart rate < 60 BPM)
  - `high_hr`: Flag (high heart rate > 100 BPM)
  - `irregular`: Flag (irregular heartbeat pattern)

**Configuration Parameters:**
- `FS_HZ`: Sampling frequency (500 Hz)
- `BASELINE_WIN`: High-pass filter window (128 samples)
- `INTEG_WIN`: Integration window (32 samples)
- `REFRACTORY_SAMPLES`: Minimum samples between peaks (256)
- `RR_HIST_LEN`: RR interval histogram length (16)
- `BLOCK_SAMPLES`: Samples per processing block (2048)
- `SAMPLES_PER_WORD`: Samples per AXI word (32)
- `LOW_BPM_THRESH` / `HIGH_BPM_THRESH`: 60 / 100 BPM thresholds
- `DETECT_THRESH`: Threshold multiplier for peak detection (6)

###### **`baseline.cpp`** – Reference Implementation
Single-block baseline implementation without HLS optimizations. Used for correctness validation and initial algorithm development.

**Key Functions:**
- `kernel_load_unpack_baseline_remove()`: Loads and processes one ECG block
- Provides algorithmic reference for validating optimized implementations

###### **`host.cpp`** – Testbench & Verification
CPU-side test harness that runs the kernel and validates results against expected outputs.

**Features:**
- Pre-computed expected results for 74 ECG blocks
- Automated verification with ±1% tolerance on numeric metrics
- Detailed error reporting for mismatches
- Test passes when all blocks match expected outputs

###### **`ecg_dram_image.h`** – Test Dataset
Pre-computed ECG data in DRAM format.

**Dataset Specifications:**
- **Total Samples**: 151,552 (300 seconds at 500 Hz)
- **Block Size**: 2,048 samples per block
- **Total Blocks**: 74
- **Data Format**: 16-bit signed samples packed into 512-bit AXI words
- **AXI Words**: 4,736 total words (32 samples × 512-bit per word)

##### Build & Configuration Files

###### **`makefile`** – Build Script
Compiles the C++ design and testbench using Vitis HLS.

**Commands:**
```bash
make        # Compile all sources
make clean  # Remove build artifacts
```

**Configuration:**
- Includes Xilinx Vitis headers and libraries
- Enables simulation support for FFT, FIR, DSP48E1
- C++11 standard with optimization disabled for simulation

###### **`script.tcl`** – Vitis HLS Synthesis Script
Automates HLS synthesis, implementation, and RTL generation.

##### Documentation

###### **`baseline_explanation.pdf`**
Detailed documentation of the baseline algorithm, including:
- Signal processing theory for ECG analysis
- Baseline removal methodology
- Peak detection approach
- Statistical calculations for heart rate metrics

###### **`report.pdf`**
Comprehensive design report containing:
- Architecture overview
- Performance analysis and optimization results
- Resource utilization metrics
- Timing and throughput analysis
- Design trade-offs and justifications

#### Algorithm Overview

##### Signal Processing Pipeline

1. **Baseline Removal (High-Pass Filter)**
   - Uses circular buffer moving average
   - Window: 128 samples
   - Removes low-frequency wander from ECG signal

2. **Derivative Computation**
   - First-order backward difference: `d[n] = hp[n] - hp[n-1]`
   - Emphasizes signal transitions

3. **Magnitude & Energy**
   - Absolute value of derivative
   - Indicates signal strength

4. **Integration (Moving Window)**
   - Window: 32 samples
   - Integrates energy to detect QRS complexes

5. **Peak Detection**
   - Threshold crossing on integrated signal
   - Refractory period: 256 samples minimum between peaks
   - Prevents multiple detections of same heartbeat

6. **RR Analysis**
   - Computes intervals between consecutive peaks
   - Calculates statistics: average, BPM, variability
   - Detects irregular patterns (>25% outliers)

##### Performance Characteristics

- **Throughput**: ~2,048 cycles per block (one sample per cycle)
- **Latency**: ~2,300 cycles (pipeline fill + processing)
- **Block Processing**: 74 blocks × ~2,048 cycles ≈ 151,552 total cycles
- **Initiation Interval**: Strict II=1 on all pipeline stages
- **Resource Efficiency**: Leverages parallelism without excessive area overhead

#### Data Types & Bit-Widths

| Type | Width | Purpose |
|------|-------|---------|
| `ecg_t` | 16-bit signed | Raw ECG samples |
| `deriv_t` | 16-bit signed | Derivative values |
| `mag_t` | 16-bit unsigned | Magnitude/energy values |
| `accum_t` | 32-bit unsigned | Accumulator for filtering |
| `rr_t` | 16-bit unsigned | RR intervals in samples |
| `bpm_t` | 16-bit unsigned | Heart rate in BPM |
| `axi_word_t` | 512-bit | AXI data bus (32 × 16-bit) |

#### Key Optimizations

##### HLS Pragmas & Techniques

1. **PIPELINE II=1**: Achieves one sample per cycle throughput
2. **DATAFLOW**: Enables concurrent execution of all pipeline stages
3. **ARRAY_PARTITION (cyclic factor=2)**: Reduces memory access conflicts
4. **UNROLL**: Fully unrolls inner loops for parallelism (k_abs_vec)
5. **INLINE OFF**: Separates burst-read from unpack to isolate AXI latency
6. **LUTRAM**: Uses distributed RAM for FIFO buffers (not BRAM36)

##### Vectorization Strategy

The `k_abs_vec` stage demonstrates **16-way parallelism**:
- Accumulates 16 serial samples into packed word
- Applies absolute value to all 16 lanes simultaneously
- Emits one packed output per cycle
- No inter-sample dependency enables full unrolling

#### Building the Design

##### Prerequisites
- Xilinx Vitis HLS 2021.1+ (or compatible version)
- C++ compiler with C++11 support
- FPGA board or simulation environment

##### Compilation
```bash
cd /Users/samsharma/Work/git/FPGA/FPGA_ECG_Accelerator
make clean
make
```

##### Running Testbench
```bash
./result
```

Expected output (on success):
```
TEST PASSED!!!
```

#### Output Format

The kernel produces an `ECGSummaryOut` structure for each block:

```cpp
struct ECGSummaryOut {
    ap_uint<16> num_peaks_detected;    // [0..127] peaks per block
    rr_t        avg_rr_samples;        // Average RR interval in samples
    ap_uint<32> avg_rr_ms;             // Average RR interval in milliseconds
    bpm_t       avg_bpm;               // Heart rate in beats per minute
    ap_uint<1>  low_hr;                // 1 if BPM < 60
    ap_uint<1>  high_hr;               // 1 if BPM > 100
    ap_uint<1>  irregular;             // 1 if >25% RR outliers
};
```

#### Verification

The testbench validates:
- Peak count ±1% tolerance
- Average RR interval (samples) ±1% tolerance
- Average RR interval (ms) ±1% tolerance
- Average BPM ±1% tolerance
- Boolean flags (low_hr, high_hr, irregular) exact match

All 74 blocks must pass for successful test execution.

#### Files Summary

| File | Lines | Purpose |
|------|-------|---------|
| `top.cpp` | 543 | Main HLS kernel with 9-stage pipeline |
| `baseline.cpp` | 341 | Reference single-block implementation |
| `host.cpp` | 169 | Testbench with validation logic |
| `dcl.h` | 50 | Type definitions & configuration |
| `ecg_dram_image.h` | 4,757 | Pre-computed ECG dataset (74 blocks) |
| `makefile` | - | Build automation |
| `script.tcl` | - | HLS synthesis automation |
| `baseline_explanation.pdf` | - | Algorithm documentation |
| `report.pdf` | - | Design report |

#### Design Insights

##### Why II=1 on Every Stage?

Each pipeline stage maintains strict II=1 initiation interval by:
1. **No loop-carried dependencies** (except the necessary ones managed by circular buffers)
2. **Single sample per cycle throughput** at input and output
3. **Cyclic array partitioning** to avoid memory access conflicts
4. **Careful buffer sizing** to prevent backpressure in FIFO chains

### Circular Buffer Pattern

Baseline and integration stages use circular buffers with the pattern:
```
new_val = (accum + input - old_val) / window_size
```

This enables:
- No accumulator reset between iterations
- Minimal resource usage (single accumulator, not array)
- Natural II=1 dependency on `accum` carry

##### Vectorization Without Hazards

The `k_abs_vec` stage exploits the **absence of inter-sample dependencies** in absolute-value computation:
- Pack 16 samples into single wide word
- Apply abs to all 16 lanes in parallel (fully unrolled)
- Maintain II=1 by using two sub-loops, each at II=1

#### References

For more details, consult:
- **baseline_explanation.pdf** – Algorithm theory
- **report.pdf** – Design analysis and results
- **top.cpp** – Detailed implementation comments
- Xilinx Vitis HLS User Guide – HLS pragmas and optimization techniques

#### Contact

For questions or contributions, please contact the repository maintainers.

---

**Repository**: FPGA_ECG_Accelerator  
**Owner**: SamSharma321  
**Branch**: main  
**Last Updated**: April 2026
