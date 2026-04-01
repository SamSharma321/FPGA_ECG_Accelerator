
#include <stdio.h>
#include <ap_fixed.h>

typedef ap_int<16>   ecg_t;
typedef ap_int<16>   deriv_t;
typedef ap_uint<16>  mag_t;
typedef ap_uint<32>  accum_t;
typedef ap_uint<16>  rr_t;
typedef ap_uint<16>  bpm_t;
typedef ap_uint<512> axi_word_t;

#define FS_HZ 500
#define BASELINE_WIN 128
#define INTEG_WIN 32
#define REFRACTORY_SAMPLES 256
#define RR_HIST_LEN 16 
#define BLOCK_SAMPLES 2048 
#define SAMPLES_PER_WORD 32 
#define LOW_BPM_THRESH 60 
#define HIGH_BPM_THRESH 100 
#define DETECT_THRESH 6

// 360 Hz
// #define FS_HZ               360
// #define BASELINE_WIN        128
// #define INTEG_WIN           54
// #define REFRACTORY_SAMPLES  72
// #define RR_HIST_LEN         8
// #define BLOCK_SAMPLES       2048
// #define SAMPLES_PER_WORD    16

// #define LOW_BPM_THRESH      60
// #define HIGH_BPM_THRESH     100
// #define DETECT_THRESH       50   // keep only if integrated signal scale supports it

struct ECGSummaryOut {
    ap_uint<16> num_peaks_detected;
    rr_t        avg_rr_samples;
    ap_uint<32> avg_rr_ms;
    bpm_t       avg_bpm;
    ap_uint<1>  low_hr;
    ap_uint<1>  high_hr;
    ap_uint<1>  irregular;
};

void top_kernel(
    const axi_word_t *ecg_dram,
    ECGSummaryOut *out_mem,
    int num_blocks
);
