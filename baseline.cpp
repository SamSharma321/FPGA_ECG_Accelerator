#include <ap_int.h>
#include <hls_stream.h>
#include <stdint.h>
#include "dcl.h"

static ecg_t unpack_sample(axi_word_t word, int idx)
{
    ap_uint<16> bits = word.range((idx * 16) + 15, idx * 16);
    return (ecg_t)bits;
}

static const int MAX_PEAKS =
    (BLOCK_SAMPLES + REFRACTORY_SAMPLES - 1) / REFRACTORY_SAMPLES;

// -----------------------------------------------------------------------------
// Kernel 1: Load, unpack, and baseline remove one ECG block
// -----------------------------------------------------------------------------
static void kernel_load_unpack_baseline_remove(
    const axi_word_t *ecg_dram,
    int blk,
    ecg_t hp_out[BLOCK_SAMPLES])
{
    ecg_t raw_block[BLOCK_SAMPLES];
    ecg_t baseline_buf[BASELINE_WIN];
    accum_t baseline_sum = 0;
    int baseline_idx = 0;
    const int words_per_block =
        (BLOCK_SAMPLES + SAMPLES_PER_WORD - 1) / SAMPLES_PER_WORD;
    int block_word_base = blk * words_per_block;

    for (int i = 0; i < BASELINE_WIN; i++) {
        baseline_buf[i] = 0;
    }

    for (int i = 0; i < BLOCK_SAMPLES; i++) {
        raw_block[i] = 0;
    }

    for (int w = 0; w < words_per_block; w++)
    {
        int dram_index = block_word_base + w;
        axi_word_t word = ecg_dram[dram_index];
        int sample_base = w * SAMPLES_PER_WORD;

        for (int lane = 0; lane < SAMPLES_PER_WORD; lane++)
        {
            int sample_index = sample_base + lane;
            if (sample_index < BLOCK_SAMPLES)
            {
                raw_block[sample_index] = unpack_sample(word, lane);
            } else {
                break;
            }
        }
    }

    for (int n = 0; n < BLOCK_SAMPLES; n++)
    {
        ecg_t sample = raw_block[n];
        ecg_t old_sample = baseline_buf[baseline_idx];

        baseline_buf[baseline_idx] = sample;
        baseline_idx = (baseline_idx + 1) % BASELINE_WIN;

        baseline_sum = baseline_sum + (accum_t)sample;
        baseline_sum = baseline_sum - (accum_t)old_sample;

        accum_t mean_acc = baseline_sum / BASELINE_WIN;
        ecg_t baseline_avg = (ecg_t)mean_acc;
        hp_out[n] = sample - baseline_avg;
    }
}

// -----------------------------------------------------------------------------
// Kernel 2 + 3: Derivative, abs-energy, and moving-window integration
// -----------------------------------------------------------------------------
static void kernel_derivative_abs_energy_integrate(
    ecg_t hp_in[BLOCK_SAMPLES],
    mag_t integrated_out[BLOCK_SAMPLES])
{
    deriv_t deriv_block[BLOCK_SAMPLES];
    mag_t integ_buf[INTEG_WIN];
    accum_t integ_sum = 0;
    int integ_idx = 0;

    for (int i = 0; i < INTEG_WIN; i++) {
        integ_buf[i] = 0;
    }

    if (BLOCK_SAMPLES > 0) {
        deriv_block[0] = (deriv_t)(hp_in[0] - (ecg_t)0);
    }

    for (int n = 1; n < BLOCK_SAMPLES; n++) {
        deriv_block[n] = (deriv_t)(hp_in[n] - hp_in[n - 1]);
    }

    for (int n = 0; n < BLOCK_SAMPLES; n++)
    {
        deriv_t d = deriv_block[n];
        mag_t energy;
        if (d < 0) {
            energy = (mag_t)(-d);
        }
        else {
            energy = (mag_t)d;
        }

        mag_t old_energy = integ_buf[integ_idx];
        integ_buf[integ_idx] = energy;
        integ_idx = (integ_idx + 1) % INTEG_WIN;

        integ_sum = integ_sum + (accum_t)energy;
        integ_sum = integ_sum - (accum_t)old_energy;
        integrated_out[n] = (mag_t)(integ_sum / INTEG_WIN);
    }
}

// -----------------------------------------------------------------------------
// Kernel 4: Peak event detection with refractory control
// Outputs accepted peak indices.
// Kept mostly same because it is already sequential in nature.
// -----------------------------------------------------------------------------
static void kernel_detect_peaks(
    mag_t integrated_in[BLOCK_SAMPLES],
    ap_uint<16> peak_positions[MAX_PEAKS],
    ap_uint<16> &peak_count)
{
    ap_uint<1> prev_above_thresh = 0;
    ap_uint<16> refractory_count = REFRACTORY_SAMPLES;
    ap_uint<16> local_peak_count = 0;

    for (int n = 0; n < BLOCK_SAMPLES; n++)
    {
        ap_uint<1> above_thresh;
        if (integrated_in[n] > DETECT_THRESH) {
            above_thresh = 1;
        }
        else {
            above_thresh = 0;
        }

        ap_uint<1> rising_edge;
        if ((above_thresh == 1) && (prev_above_thresh == 0)) {
            rising_edge = 1;
        }
        else {
            rising_edge = 0;
        }

        if (refractory_count < REFRACTORY_SAMPLES) {
            refractory_count = refractory_count + 1;
        }

        if ((rising_edge == 1) && (refractory_count >= REFRACTORY_SAMPLES))
        {
            if (local_peak_count < MAX_PEAKS) {
                peak_positions[local_peak_count] = (ap_uint<16>)n;
                local_peak_count = local_peak_count + 1;
            }
            refractory_count = 0;
        }

        prev_above_thresh = above_thresh;
    }
    peak_count = local_peak_count;
}

// -----------------------------------------------------------------------------
// Kernel 5: RR interval statistics and irregularity estimation
// -----------------------------------------------------------------------------
static void kernel_rr_analysis(
    ap_uint<16> peak_positions[MAX_PEAKS],
    ap_uint<16> peak_count,
    ap_uint<32> &rr_sum_samples,
    ap_uint<16> &rr_count,
    ap_uint<16> &irregular_count)
{
    rr_t rr_hist[RR_HIST_LEN];
    ap_uint<32> rr_hist_sum = 0;
    ap_uint<16> rr_hist_count = 0;
    int rr_idx = 0;

    rr_sum_samples = 0;
    rr_count = 0;
    irregular_count = 0;

    for (int i = 0; i < RR_HIST_LEN; i++) {
        rr_hist[i] = 0;
    }

    if (peak_count < 2) {
        return;
    }

    for (int p = 1; p < peak_count; p++) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=7 avg=4
        ap_uint<16> curr_peak = peak_positions[p];
        ap_uint<16> prev_peak = peak_positions[p - 1];
        rr_t current_rr = (rr_t)(curr_peak - prev_peak);

        if (current_rr != 0) {
            rr_sum_samples = rr_sum_samples + current_rr;
            rr_count = rr_count + 1;

            rr_t old_rr = rr_hist[rr_idx];
            if (old_rr != 0) {
                rr_hist_sum = rr_hist_sum - old_rr;
            }

            rr_hist[rr_idx] = current_rr;
            rr_hist_sum = rr_hist_sum + current_rr;
            if (rr_hist_count < RR_HIST_LEN) {
                rr_hist_count = rr_hist_count + 1;
            }
            rr_idx = (rr_idx + 1) % RR_HIST_LEN;

            if (rr_hist_count >= 2) {
                rr_t rr_mean = (rr_t)(rr_hist_sum / rr_hist_count);
                rr_t diff;

                if (current_rr > rr_mean) {
                    diff = (rr_t)(current_rr - rr_mean);
                }
                else {
                    diff = (rr_t)(rr_mean - current_rr);
                }

                if (rr_mean != 0) {
                    ap_uint<32> lhs = (ap_uint<32>)diff * 100;
                    ap_uint<32> rhs = (ap_uint<32>)rr_mean * 20;

                    if (lhs > rhs) {
                        irregular_count = irregular_count + 1;
                    }
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Kernel 6: Final summary generation
// -----------------------------------------------------------------------------
static void kernel_finalize_summary(
    ap_uint<16> peak_count,
    ap_uint<32> rr_sum_samples,
    ap_uint<16> rr_count,
    ap_uint<16> irregular_count,
    ECGSummaryOut &out)
{
    rr_t avg_rr_samples = 0;
    ap_uint<32> avg_rr_ms = 0;
    bpm_t avg_bpm = 0;
    ap_uint<1> low_flag = 0;
    ap_uint<1> high_flag = 0;
    ap_uint<1> irregular_flag = 0;

    if (rr_count != 0) {
        avg_rr_samples = (rr_t)(rr_sum_samples / rr_count);
        avg_rr_ms = ((ap_uint<32>)avg_rr_samples * 1000) / FS_HZ;

        if (avg_rr_samples != 0) {
            avg_bpm = (bpm_t)((60 * FS_HZ) / avg_rr_samples);
        }

        if (avg_bpm < LOW_BPM_THRESH) {
            low_flag = 1;
        }

        if (avg_bpm > HIGH_BPM_THRESH) {
            high_flag = 1;
        }

        ap_uint<32> irregular_lhs = (ap_uint<32>)irregular_count * 100;
        ap_uint<32> irregular_rhs = (ap_uint<32>)rr_count * 25;

        if (irregular_lhs > irregular_rhs) {
            irregular_flag = 1;
        }
    }

    out.num_peaks_detected = peak_count;
    out.avg_rr_samples = avg_rr_samples;
    out.avg_rr_ms = avg_rr_ms;
    out.avg_bpm = avg_bpm;
    out.low_hr = low_flag;
    out.high_hr = high_flag;
    out.irregular = irregular_flag;
}

// -----------------------------------------------------------------------------
// Sequential baseline block processor
// -----------------------------------------------------------------------------
static void process_block_baseline(
    const axi_word_t *ecg_dram,
    int blk,
    ECGSummaryOut &out)
{
    ecg_t hp_block[BLOCK_SAMPLES];
    mag_t integrated_block[BLOCK_SAMPLES];
    ap_uint<16> peak_positions[MAX_PEAKS];

    ap_uint<16> peak_count = 0;
    ap_uint<32> rr_sum_samples = 0;
    ap_uint<16> rr_count = 0;
    ap_uint<16> irregular_count = 0;

    kernel_load_unpack_baseline_remove(ecg_dram, blk, hp_block);
    kernel_derivative_abs_energy_integrate(hp_block, integrated_block);
    kernel_detect_peaks(integrated_block, peak_positions, peak_count);
    kernel_rr_analysis(peak_positions, peak_count, rr_sum_samples, rr_count, irregular_count);
    kernel_finalize_summary(peak_count, rr_sum_samples, rr_count, irregular_count, out);
}

// -----------------------------------------------------------------------------
// Top-level baseline kernel
// -----------------------------------------------------------------------------
void top_kernel(
    const axi_word_t *ecg_dram,
    ECGSummaryOut *out_mem,
    int num_blocks)
{
#pragma HLS INTERFACE m_axi port=ecg_dram offset=slave bundle=gmem0 \
    depth=4736 max_read_burst_length=64 num_read_outstanding=4 latency=96
#pragma HLS INTERFACE m_axi port=out_mem offset=slave bundle=gmem1 \
    depth=74 max_write_burst_length=16 num_write_outstanding=2 latency=32
#pragma HLS INTERFACE s_axilite port=ecg_dram
#pragma HLS INTERFACE s_axilite port=out_mem
#pragma HLS INTERFACE s_axilite port=num_blocks
#pragma HLS INTERFACE s_axilite port=return

    for (int blk = 0; blk < num_blocks; blk++)
    {
#pragma HLS LOOP_TRIPCOUNT min=0 max=74 avg=37
        ECGSummaryOut result;
        process_block_baseline(ecg_dram, blk, result);
        out_mem[blk] = result;
    }
}
