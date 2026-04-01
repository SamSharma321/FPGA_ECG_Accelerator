// =============================================================================
// top.cpp  –  Vitis HLS Optimized ECG Kernel
// =============================================================================
//
// ARCHITECTURE
// ─────────────────────────────────────────────────────────────────────────────
//
//  ┌────────────┐   s_wide    ┌────────────┐   s_ecg (ecg_t, 1/cyc)
//  │  k_load    │────────────▶│ k_unpack   │──────────────────────────────┐
//  │ (burst AXI)│             │ (FACTOR    │                              │
//  └────────────┘             │  parallel) │                              ▼
//                             └────────────┘               ┌─────────────────────┐
//                                                           │   k_baseline        │
//                                                           │ (circular buf II=1) │
//                                                           └──────────┬──────────┘
//                                                                      │ s_hp
//                                                                      ▼
//                                                           ┌─────────────────────┐
//                                                           │   k_deriv           │
//                                                           │   (diff,  II=1)     │
//                                                           └──────────┬──────────┘
//                                                                      │ s_deriv
//                                                         ┌────────────▼───────────┐
//                                                         │  k_abs_vec             │
//                                                         │  (FACTOR lanes parallel│
//                                                         │   packed wide stream)  │
//                                                         └────────────┬───────────┘
//                                                                      │ s_energy_wide
//                                                         ┌────────────▼───────────┐
//                                                         │  k_integrate           │
//                                                         │  (circular buf II=1)   │
//                                                         └────────────┬───────────┘
//                                                                      │ s_integrated
//                                                         ┌────────────▼───────────┐
//                                                         │  k_detect_peaks        │
//                                                         │  (refractory gate)     │
//                                                         └────────────┬───────────┘
//                                                                      │ peak_positions[]
//                                                                      │ peak_count
//                                                                      ▼
//                                              ┌──────────────────────────────────────┐
//                                              │  k_rr_analysis  (O(peaks), unrolled) │
//                                              └──────────────────────┬───────────────┘
//                                                                     │
//                                                                     ▼
//                                              ┌──────────────────────────────────────┐
//                                              │  k_finalize_summary                  │
//                                              └──────────────────────────────────────┘
//
// =============================================================================

#include "dcl.h"
#include <ap_int.h>
#include <hls_stream.h>
#include <stdint.h>

// -----------------------------------------------------------------------------
// SAMPLES_PER_WORD = 256 / 16 = 16  (16 × 16-bit samples per axi_word_t)
// FACTOR           – parallelism for vectorised abs stage; must divide
//                    BLOCK_SAMPLES and SAMPLES_PER_WORD.
// We pick FACTOR = SAMPLES_PER_WORD so one AXI word = one packed element.
// -----------------------------------------------------------------------------
static const int FACTOR         = SAMPLES_PER_WORD;           // e.g., 16
static const int MAG_BITS       = sizeof(mag_t) * 8;          // width of mag_t in bits
static const int DERIV_BITS     = sizeof(deriv_t) * 8;        // width of deriv_t in bits
static const int CHUNKS_PER_BLK = BLOCK_SAMPLES / FACTOR;     // packed chunks per block

static_assert((BLOCK_SAMPLES % FACTOR) == 0,
              "BLOCK_SAMPLES must be divisible by FACTOR=SAMPLES_PER_WORD");

// Wide packed stream types: FACTOR lanes of mag_t / deriv_t per token
typedef ap_uint<DERIV_BITS * FACTOR> wide_deriv_t;   // FACTOR deriv_t packed
typedef ap_uint<MAG_BITS   * FACTOR> wide_mag_t;     // FACTOR mag_t  packed

// Convenience bit-field accessor (same macro pattern as reference)
#define LANE_RANGE(stream_var, idx, W) \
    (stream_var).range(((idx) * (W)) + (W) - 1, (idx) * (W))

// =============================================================================
// Stage 1 – k_load
// Burst-reads all AXI words for block `blk` into a local register file,
// then streams them as axi_word_t tokens (one word per cycle).
// Separating burst-read from unpack keeps both inner loops at strict II=1
// and avoids mixing m_axi latency with compute logic.
// =============================================================================
static void k_load(
    const axi_word_t *__restrict__ ecg_dram,
    int                             blk,
    hls::stream<axi_word_t>        &s_wide
) {
#pragma HLS INLINE off
    const int WPB  = (BLOCK_SAMPLES + SAMPLES_PER_WORD - 1) / SAMPLES_PER_WORD;
    const int base = blk * WPB;

    // Burst read: BRAM register file absorbs full burst in one shot.
    axi_word_t wbuf[WPB];
#pragma HLS ARRAY_PARTITION variable=wbuf type=complete dim=1

    burst_rd: for (int w = 0; w < WPB; w++) {
#pragma HLS PIPELINE II=1
        wbuf[w] = ecg_dram[base + w];
    }

    // Stream packed AXI words downstream (one word = SAMPLES_PER_WORD samples)
    stream_out: for (int w = 0; w < WPB; w++) {
#pragma HLS PIPELINE II=1
        s_wide.write(wbuf[w]);
    }
}

// =============================================================================
// Stage 2 – k_unpack
// Unpacks FACTOR=SAMPLES_PER_WORD samples from each AXI word in parallel
// (fully unrolled inner loop → FACTOR concurrent 16-bit extractions per cycle)
// and serialises them as individual ecg_t tokens to the filter chain.
// =============================================================================
static void k_unpack(
    hls::stream<axi_word_t> &s_wide,
    hls::stream<ecg_t>      &s_ecg
) {
#pragma HLS INLINE off
    const int WPB = BLOCK_SAMPLES / FACTOR;   // exact after static_assert

    unpack_loop: for (int w = 0; w < WPB; w++) {
#pragma HLS PIPELINE II=FACTOR   // emits FACTOR tokens; II = FACTOR so s_ecg sees 1/cyc
        axi_word_t word = s_wide.read();

        // Unrolled extraction: FACTOR concurrent 16-bit slices in RTL
        for (int lane = 0; lane < FACTOR; lane++) {
#pragma HLS UNROLL
            ap_uint<16> bits = word.range((lane << 4) + 15, lane << 4);
            s_ecg.write((ecg_t)bits);
        }
    }
}

// =============================================================================
// Stage 3 – k_baseline
// Causal moving-average high-pass filter using a circular buffer.
//
// Analysis:
//   Loop carries:  sum[n] = sum[n-1] + x[n] - cbuf[idx]
//                  cbuf[idx] = x[n]    (write after the dependent read)
//   The same cbuf slot recurs every BASELINE_WIN iterations.
//   Provided BASELINE_WIN > pipeline depth (typically 3-5 stages), there
//   is NO true inter-iteration hazard → DEPENDENCE inter RAW false is safe.
//
//   Division by BASELINE_WIN (compile-time constant) is folded into a
//   multiply-by-reciprocal by HLS (no hardware divider, no extra latency).
// =============================================================================
static void k_baseline(
    hls::stream<ecg_t> &s_in,
    hls::stream<ecg_t> &s_out
) {
#pragma HLS INLINE off

    ecg_t   cbuf[BASELINE_WIN];
#pragma HLS ARRAY_PARTITION variable=cbuf type=cyclic factor=2 dim=1
    accum_t sum = 0;
    int     idx = 0;

    // Initialise before the hot loop so HLS sees no reset overhead in II=1 path
    init_bl: for (int i = 0; i < BASELINE_WIN; i++) {
#pragma HLS PIPELINE II=1
        cbuf[i] = (ecg_t)0;
    }

    bl_main: for (int n = 0; n < BLOCK_SAMPLES; n++) {
#pragma HLS PIPELINE II=1
// #pragma HLS DEPENDENCE variable=cbuf inter RAW false
        ecg_t x   = s_in.read();
        ecg_t old = cbuf[idx];         // read from slot before writing
        cbuf[idx] = x;
        idx       = (idx + 1 == BASELINE_WIN) ? 0 : idx + 1;
        sum      += (accum_t)x - (accum_t)old;
        ecg_t avg = (ecg_t)(sum / (accum_t)BASELINE_WIN);
        s_out.write(x - avg);
    }
}

// =============================================================================
// Stage 4 – k_deriv
// d[n] = hp[n] – hp[n-1]   (first-order backward difference)
// Pure carry on `prev`: natural II=1 dependency with no resource conflict.
// =============================================================================
static void k_deriv(
    hls::stream<ecg_t>   &s_in,
    hls::stream<deriv_t> &s_out
) {
#pragma HLS INLINE off
    ecg_t prev = (ecg_t)0;

    deriv_main: for (int n = 0; n < BLOCK_SAMPLES; n++) {
#pragma HLS PIPELINE II=1
        ecg_t curr = s_in.read();
        s_out.write((deriv_t)(curr - prev));
        prev = curr;
    }
}

// =============================================================================
// Stage 5 – k_abs_vec   (VECTORISED absolute-value stage)
// This is the only stage with NO inter-sample dependency.
// We exploit this by accumulating FACTOR deriv_t tokens into one wide_deriv_t
// packed word (via a shift-register accumulator at II=1), then applying abs()
// to all FACTOR lanes simultaneously in a single fully-unrolled pass,
// and emitting a packed wide_mag_t word for the integrator to consume.
//
// Result: FACTOR-way parallelism at II=1 — the abs computation for an entire
// AXI-word's worth of samples occupies exactly one pipeline cycle.
// =============================================================================
static void k_abs_vec(
    hls::stream<deriv_t>    &s_in,        // serial input, 1 sample/cycle
    hls::stream<wide_mag_t> &s_out        // packed output, FACTOR samples/token
) {
#pragma HLS INLINE off

    const int WPB = CHUNKS_PER_BLK;      // tokens per block = BLOCK_SAMPLES/FACTOR

    // Accumulate FACTOR serial samples into one packed word, then process in bulk.
    // Two-phase design keeps each sub-loop at strict II=1.

    pack_loop: for (int chunk = 0; chunk < WPB; chunk++) {
        // --- Phase A: pack FACTOR deriv_t tokens into one wide register ---
        wide_deriv_t packed = (wide_deriv_t)0;
        pack_inner: for (int lane = 0; lane < FACTOR; lane++) {
#pragma HLS PIPELINE II=1
            deriv_t d = s_in.read();
            LANE_RANGE(packed, lane, DERIV_BITS) = d.range(DERIV_BITS - 1, 0);
        }

        // --- Phase B: apply |abs| to all FACTOR lanes in one cycle ---
        wide_mag_t result = (wide_mag_t)0;
        abs_unroll: for (int lane = 0; lane < FACTOR; lane++) {
#pragma HLS UNROLL
            // Reinterpret stored bits as signed deriv_t
            deriv_t d;
            d.range(DERIV_BITS - 1, 0) = LANE_RANGE(packed, lane, DERIV_BITS);
            mag_t e = (d < (deriv_t)0) ? (mag_t)(-d) : (mag_t)d;
            LANE_RANGE(result, lane, MAG_BITS) = e.range(MAG_BITS - 1, 0);
        }
        s_out.write(result);
    }
}

// =============================================================================
// Stage 6 – k_integrate
// Causal moving-window integrator; unpacks each wide_mag_t token from
// k_abs_vec one lane at a time and feeds the circular-buffer accumulator.
// Same DEPENDENCE / ARRAY_PARTITION treatment as k_baseline.
// =============================================================================
static void k_integrate(
    hls::stream<wide_mag_t> &s_in,
    hls::stream<mag_t>      &s_out
) {
#pragma HLS INLINE off

    mag_t   cbuf[INTEG_WIN];
#pragma HLS ARRAY_PARTITION variable=cbuf type=cyclic factor=2 dim=1
    accum_t sum = 0;
    int     idx = 0;

    init_ig: for (int i = 0; i < INTEG_WIN; i++) {
#pragma HLS PIPELINE II=1
        cbuf[i] = (mag_t)0;
    }

    const int WPB = CHUNKS_PER_BLK;

    ig_outer: for (int chunk = 0; chunk < WPB; chunk++) {
        wide_mag_t packed = s_in.read();   // read one packed token (FACTOR lanes)

        ig_inner: for (int lane = 0; lane < FACTOR; lane++) {
#pragma HLS PIPELINE II=1
// #pragma HLS DEPENDENCE variable=cbuf inter RAW false
            mag_t e;
            e.range(MAG_BITS - 1, 0) = LANE_RANGE(packed, lane, MAG_BITS);

            mag_t old = cbuf[idx];
            cbuf[idx] = e;
            idx       = (idx + 1 == INTEG_WIN) ? 0 : idx + 1;
            sum      += (accum_t)e - (accum_t)old;
            s_out.write((mag_t)(sum / (accum_t)INTEG_WIN));
        }
    }
}

// =============================================================================
// Stage 7 – k_detect_peaks
// Rising-edge threshold crossing gated by a refractory counter.
// Writes accepted peak sample-positions into a flat array (sequential writes
// → HLS DATAFLOW ping-pong buffer compatible).
// peak_count is written exactly once at the end (valid DATAFLOW scalar output).
// =============================================================================
static void k_detect_peaks(
    hls::stream<mag_t> &s_in,
    ap_uint<16>         peak_positions[BLOCK_SAMPLES],
    ap_uint<16>        &peak_count
) {
#pragma HLS INLINE off

    ap_uint<1>  prev_above  = 0;
    ap_uint<16> refract_cnt = (ap_uint<16>)REFRACTORY_SAMPLES;  // starts "cooled"
    ap_uint<16> cnt         = 0;

    detect_main: for (int n = 0; n < BLOCK_SAMPLES; n++) {
#pragma HLS PIPELINE II=1
        mag_t v = s_in.read();

        ap_uint<1> above  = (v > (mag_t)DETECT_THRESH) ? ap_uint<1>(1) : ap_uint<1>(0);
        ap_uint<1> rising = above & (~prev_above & ap_uint<1>(1));

        if (refract_cnt < (ap_uint<16>)REFRACTORY_SAMPLES)
            refract_cnt++;

        if (rising == ap_uint<1>(1) && refract_cnt >= (ap_uint<16>)REFRACTORY_SAMPLES) {
            peak_positions[cnt] = (ap_uint<16>)n;
            cnt++;
            refract_cnt = 0;
        }
        prev_above = above;
    }
    peak_count = cnt;
}

// =============================================================================
// Stage 8 – k_rr_analysis  (O(num_peaks), not O(BLOCK_SAMPLES))
// Runs sequentially AFTER the DATAFLOW pipeline completes.
//
// Inner RR_HIST_LEN loop: fully unrolled → balanced binary adder tree.
// Outer peak loop: PIPELINE II=1 (history array is partitioned completely
// so all RR_HIST_LEN entries are simultaneously readable as registers).
//
// Division by valid_cnt: variable but bounded to RR_HIST_LEN (small, e.g. 8).
// HLS will synthesise this as a small LUT-based divider; the loop II is 1
// because the division result is not carried into the next iteration.
// =============================================================================
static void k_rr_analysis(
    ap_uint<16>  peak_positions[BLOCK_SAMPLES],
    ap_uint<16>  peak_count,
    ap_uint<32> &rr_sum_samples,
    ap_uint<16> &rr_count,
    ap_uint<16> &irregular_count
) {
#pragma HLS INLINE off

    rr_t rr_hist[RR_HIST_LEN];
#pragma HLS ARRAY_PARTITION variable=rr_hist type=complete dim=1

    int rr_idx = 0;
    rr_sum_samples  = 0;
    rr_count        = 0;
    irregular_count = 0;

    // Initialise history array (fully unrolled → single cycle in RTL)
    for (int i = 0; i < RR_HIST_LEN; i++) {
#pragma HLS UNROLL
        rr_hist[i] = (rr_t)0;
    }

    if (peak_count < 2) return;

    rr_main: for (int p = 1; p < (int)peak_count; p++) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=7 avg=4
#pragma HLS PIPELINE
        rr_t cur_rr = (rr_t)((ap_uint<16>)(peak_positions[p] - peak_positions[p - 1]));
        if (cur_rr == (rr_t)0) continue;

        rr_sum_samples += (ap_uint<32>)cur_rr;
        rr_count++;

        rr_hist[rr_idx] = cur_rr;
        rr_idx = (rr_idx + 1 == RR_HIST_LEN) ? 0 : rr_idx + 1;

        // ---- Parallel adder tree over RR_HIST_LEN entries (fully unrolled) ----
        ap_uint<32> local_sum = 0;
        ap_uint<16> valid_cnt = 0;
        for (int k = 0; k < RR_HIST_LEN; k++) {
#pragma HLS UNROLL
            if (rr_hist[k] != (rr_t)0) {
                local_sum += (ap_uint<32>)rr_hist[k];
                valid_cnt++;
            }
        }
#pragma HLS bind_op variable=local_sum op=add impl=fabric latency=1

        if (valid_cnt >= 2) {
            rr_t rr_mean = (rr_t)(local_sum / (ap_uint<32>)valid_cnt);
            rr_t diff    = (cur_rr > rr_mean) ? (rr_t)(cur_rr - rr_mean)
                                               : (rr_t)(rr_mean - cur_rr);
            if (rr_mean != (rr_t)0 &&
                ((ap_uint<32>)diff * 100) > ((ap_uint<32>)rr_mean * 20)) {
                irregular_count++;
            }
        }
    }
}

// =============================================================================
// Stage 9 – k_finalize_summary
// Pure combinational/single-iteration computation.
// =============================================================================
static void k_finalize_summary(
    ap_uint<16>    peak_count,
    ap_uint<32>    rr_sum_samples,
    ap_uint<16>    rr_count,
    ap_uint<16>    irregular_count,
    ECGSummaryOut &out
) {
#pragma HLS INLINE off
#pragma HLS PIPELINE II=1    // single-iteration body → effectively combinational

    rr_t        avg_rr = 0;
    ap_uint<32> rr_ms  = 0;
    bpm_t       bpm    = 0;
    ap_uint<1>  lo = 0, hi = 0, irr = 0;

    if (rr_count != 0) {
        avg_rr = (rr_t)(rr_sum_samples / (ap_uint<32>)rr_count);
        rr_ms  = ((ap_uint<32>)avg_rr * 1000u) / (ap_uint<32>)FS_HZ;

        if (avg_rr != (rr_t)0)
            bpm = (bpm_t)(((ap_uint<32>)60 * (ap_uint<32>)FS_HZ) / (ap_uint<32>)avg_rr);

        lo  = (bpm < (bpm_t)LOW_BPM_THRESH)  ? ap_uint<1>(1) : ap_uint<1>(0);
        hi  = (bpm > (bpm_t)HIGH_BPM_THRESH) ? ap_uint<1>(1) : ap_uint<1>(0);

        if (((ap_uint<32>)irregular_count * 100u) > ((ap_uint<32>)rr_count * 25u))
            irr = 1;
    }

    out.num_peaks_detected = peak_count;
    out.avg_rr_samples     = avg_rr;
    out.avg_rr_ms          = rr_ms;
    out.avg_bpm            = bpm;
    out.low_hr             = lo;
    out.high_hr            = hi;
    out.irregular          = irr;
}

// =============================================================================
// process_block_df – DATAFLOW wrapper (Stages 1-7)
// ─────────────────────────────────────────────────────────────────────────────
// All seven stages are instantiated as concurrent processes by DATAFLOW.
// Each stage communicates exclusively through its output hls::stream (or, for
// the final stage, through the peak_positions array / peak_count scalar).
//
// FIFO depth sizing rationale:
//   • s_wide / s_ecg / s_hp / s_deriv / s_integrated : depth=32 is adequate
//     because every stage runs at exactly II=1 sample/cycle → no skew.
//   • s_energy_wide : depth=4 (one packed token per FACTOR samples; minimal
//     buffering needed between k_abs_vec and k_integrate).
//   • LUTRAM for short FIFOs to avoid consuming BRAM36 blocks.
// =============================================================================
static void process_block_df(
    const axi_word_t *ecg_dram,
    int               blk,
    ap_uint<16>       peak_positions[BLOCK_SAMPLES],
    ap_uint<16>      &peak_count
) {
#pragma HLS DATAFLOW

    // ---- Inter-stage FIFOs ----
    hls::stream<axi_word_t>  s_wide       ("s_wide");
    hls::stream<ecg_t>       s_ecg        ("s_ecg");
    hls::stream<ecg_t>       s_hp         ("s_hp");
    hls::stream<deriv_t>     s_deriv      ("s_deriv");
    hls::stream<wide_mag_t>  s_energy_wide("s_energy_wide");
    hls::stream<mag_t>       s_integrated ("s_integrated");

#pragma HLS stream variable=s_wide        depth=32
#pragma HLS stream variable=s_ecg         depth=32
#pragma HLS stream variable=s_hp          depth=32
#pragma HLS stream variable=s_deriv       depth=32
#pragma HLS stream variable=s_energy_wide depth=4
#pragma HLS stream variable=s_integrated  depth=32

#pragma HLS bind_storage variable=s_wide        type=fifo impl=lutram
#pragma HLS bind_storage variable=s_ecg         type=fifo impl=lutram
#pragma HLS bind_storage variable=s_hp          type=fifo impl=lutram
#pragma HLS bind_storage variable=s_deriv       type=fifo impl=lutram
#pragma HLS bind_storage variable=s_energy_wide type=fifo impl=lutram
#pragma HLS bind_storage variable=s_integrated  type=fifo impl=lutram

    // ---- Seven concurrent pipeline stages ----
    k_load(ecg_dram, blk, s_wide);
    k_unpack(s_wide, s_ecg);
    k_baseline(s_ecg, s_hp);
    k_deriv(s_hp, s_deriv);
    k_abs_vec(s_deriv, s_energy_wide);
    k_integrate(s_energy_wide, s_integrated);
    k_detect_peaks(s_integrated, peak_positions, peak_count);
}

// =============================================================================
// top_kernel – AXI master / AXI-Lite slave top-level
// ─────────────────────────────────────────────────────────────────────────────
// For each block:
//   1. process_block_df  runs the 7-stage streaming pipeline under DATAFLOW.
//      Throughput ≈ BLOCK_SAMPLES cycles (II=1 per stage, all overlapped).
//   2. k_rr_analysis     runs sequentially on the O(peaks) peak list.
//   3. k_finalize_summary computes the ECGSummaryOut struct.
//   4. Result is written to out_mem via a separate AXI master port (gmem1).
//
// =============================================================================
void top_kernel(
    const axi_word_t *ecg_dram,
    ECGSummaryOut    *out_mem,
    int               num_blocks
) {
#pragma HLS INTERFACE m_axi port=ecg_dram  offset=slave bundle=gmem0 \
                                depth=4736 max_widen_bitwidth=256 max_read_burst_length=16
#pragma HLS INTERFACE m_axi port=out_mem   offset=slave bundle=gmem1 depth=74
#pragma HLS INTERFACE s_axilite port=ecg_dram
#pragma HLS INTERFACE s_axilite port=out_mem
#pragma HLS INTERFACE s_axilite port=num_blocks
#pragma HLS INTERFACE s_axilite port=return

    // Local peak buffer – reused across blocks (no reset needed; k_detect_peaks
    // overwrites the valid prefix [0..peak_count-1] each time).
    ap_uint<16> peak_positions[BLOCK_SAMPLES];
#pragma HLS ARRAY_PARTITION variable=peak_positions type=cyclic factor=4 dim=1

    block_loop: for (int blk = 0; blk < num_blocks; blk++) {
#pragma HLS LOOP_TRIPCOUNT min=0 max=74 avg=37
        // ---- Phase A: streaming pipeline (DATAFLOW, ≈BLOCK_SAMPLES cycles) ----
        ap_uint<16> peak_count = 0;
        process_block_df(ecg_dram, blk, peak_positions, peak_count);

        // ---- Phase B: RR statistics (O(peaks), pipelined at II=1) ----
        ap_uint<32> rr_sum  = 0;
        ap_uint<16> rr_cnt  = 0;
        ap_uint<16> irr_cnt = 0;
        k_rr_analysis(peak_positions, peak_count, rr_sum, rr_cnt, irr_cnt);

        // ---- Phase C: summary generation (single iteration) ----
        ECGSummaryOut result;
        k_finalize_summary(peak_count, rr_sum, rr_cnt, irr_cnt, result);

        out_mem[blk] = result;
    }
}
