#include <stdio.h>
#include <stdlib.h>
#include "dcl.h"
#include <iostream>
#include "ecg_dram_image.h"

struct ExpectedSummary {
    unsigned int num_peaks_detected;
    unsigned int avg_rr_samples;
    unsigned int avg_rr_ms;
    unsigned int avg_bpm;
    unsigned int low_hr;
    unsigned int high_hr;
    unsigned int irregular;
};

static const ExpectedSummary kExpected[ECG_NUM_BLOCKS] = {
/* Peaks | Avg RR Samp | Avg RR ms | Avg BPM | Low HR | High HR | Irregular */
    {6,       377,         754,        79,       0,       0,          0},
    {5,       384,         768,        78,       0,       0,          0},
    {6,       358,         716,        83,       0,       0,          0},
    {5,       375,         750,        80,       0,       0,          0},
    {6,       386,         772,        77,       0,       0,          0},
    {5,       413,         826,        72,       0,       0,          0},
    {5,       411,         822,        72,       0,       0,          0},
    {5,       394,         788,        76,       0,       0,          0},
    {5,       420,         840,        71,       0,       0,          0},
    {5,       393,         786,        76,       0,       0,          0},
    {5,       421,         842,        71,       0,       0,          0},
    {5,       397,         794,        75,       0,       0,          0},
    {5,       393,         786,        76,       0,       0,          0},
    {6,       375,         750,        80,       0,       0,          0},
    {5,       403,         806,        74,       0,       0,          0},
    {5,       389,         778,        77,       0,       0,          0},
    {5,       413,         826,        72,       0,       0,          0},
    {5,       421,         842,        71,       0,       0,          0},
    {5,       400,         800,        75,       0,       0,          0},
    {6,       347,         694,        86,       0,       0,          0},
    {4,       429,         858,        69,       0,       0,          0},
    {6,       398,         796,        75,       0,       0,          0},
    {5,       368,         736,        81,       0,       0,          0},
    {5,       419,         838,        71,       0,       0,          0},
    {5,       398,         796,        75,       0,       0,          0},
    {5,       407,         814,        73,       0,       0,          0},
    {5,       424,         848,        70,       0,       0,          0},
    {5,       406,         812,        73,       0,       0,          0},
    {5,       407,         814,        73,       0,       0,          0},
    {5,       406,         812,        73,       0,       0,          0},
    {5,       402,         804,        74,       0,       0,          0},
    {5,       424,         848,        70,       0,       0,          0},
    {5,       427,         854,        70,       0,       0,          0},
    {5,       387,         774,        77,       0,       0,          0},
    {5,       377,         754,        79,       0,       0,          0},
    {5,       372,         744,        80,       0,       0,          0},
    {5,       380,         760,        78,       0,       0,          0},
    {6,       372,         744,        80,       0,       0,          0},
    {5,       403,         806,        74,       0,       0,          0},
    {5,       395,         790,        75,       0,       0,          0},
    {5,       411,         822,        72,       0,       0,          0},
    {6,       386,         772,        77,       0,       0,          0},
    {5,       371,         742,        80,       0,       0,          0},
    {5,       417,         834,        71,       0,       0,          0},
    {5,       420,         840,        71,       0,       0,          0},
    {5,       413,         826,        72,       0,       0,          0},
    {5,       406,         812,        73,       0,       0,          0},
    {5,       366,         732,        81,       0,       0,          0},
    {5,       413,         826,        72,       0,       0,          0},
    {5,       429,         858,        69,       0,       0,          0},
    {5,       424,         848,        70,       0,       0,          0},
    {5,       401,         802,        74,       0,       0,          0},
    {5,       420,         840,        71,       0,       0,          0},
    {6,       363,         726,        82,       0,       0,          0},
    {5,       395,         790,        75,       0,       0,          0},
    {5,       356,         712,        84,       0,       0,          0},
    {5,       413,         826,        72,       0,       0,          0},
    {6,       385,         770,        77,       0,       0,          0},
    {5,       403,         806,        74,       0,       0,          0},
    {5,       403,         806,        74,       0,       0,          0},
    {5,       390,         780,        76,       0,       0,          0},
    {6,       366,         732,        81,       0,       0,          0},
    {5,       398,         796,        75,       0,       0,          0},
    {5,       398,         796,        75,       0,       0,          0},
    {5,       402,         804,        74,       0,       0,          0},
    {6,       378,         756,        79,       0,       0,          0},
    {4,       419,         838,        71,       0,       0,          0},
    {6,       384,         768,        79,       0,       0,          0},
    {5,       386,         772,        77,       0,       0,          0},
    {5,       405,         810,        74,       0,       0,          0},
    {5,       404,         808,        74,       0,       0,          0},
    {5,       406,         812,        73,       0,       0,          0},
    {5,       394,         788,        76,       0,       0,          0},
    {2,       409,         818,        73,       0,       0,          0}
};

static unsigned int uabs_diff(unsigned int a, unsigned int b)
{
    return (a > b) ? (a - b) : (b - a);
}

static unsigned int tolerance_1pct(unsigned int expected)
{
    if (expected == 0) {
        return 0;
    }
    return (expected + 99) / 100;
}

static bool check_metric(
    int block_idx,
    const char *name,
    unsigned int actual,
    unsigned int expected)
{
    unsigned int tolerance = tolerance_1pct(expected);
    unsigned int diff = uabs_diff(actual, expected);
    if (diff <= tolerance) {
        return true;
    }

    std::cout << "Mismatch in block " << block_idx << " for " << name
              << ": expected " << expected
              << ", got " << actual
              << ", tolerance +/-" << tolerance << "\n";
    return false;
}

static bool check_flag(
    int block_idx,
    const char *name,
    unsigned int actual,
    unsigned int expected)
{
    if (actual == expected) {
        return true;
    }

    std::cout << "Mismatch in block " << block_idx << " for " << name
              << ": expected " << expected
              << ", got " << actual << "\n";
    return false;
}

int main()
{
    ECGSummaryOut out[ECG_NUM_BLOCKS];
    bool pass = true;

    top_kernel(ecg_dram, out, ECG_NUM_BLOCKS);

    for (int i = 0; i < ECG_NUM_BLOCKS; ++i) {
        const ExpectedSummary &exp = kExpected[i];
        pass &= check_metric(i, "Peaks", (unsigned int)out[i].num_peaks_detected, exp.num_peaks_detected);
        pass &= check_metric(i, "Avg RR samp", (unsigned int)out[i].avg_rr_samples, exp.avg_rr_samples);
        pass &= check_metric(i, "Avg RR ms", (unsigned int)out[i].avg_rr_ms, exp.avg_rr_ms);
        pass &= check_metric(i, "Avg BPM", (unsigned int)out[i].avg_bpm, exp.avg_bpm);
        pass &= check_flag(i, "Low HR", (unsigned int)out[i].low_hr, exp.low_hr);
        pass &= check_flag(i, "High HR", (unsigned int)out[i].high_hr, exp.high_hr);
        pass &= check_flag(i, "Irregular", (unsigned int)out[i].irregular, exp.irregular);
    }

    if (!pass) {
        std::cout << "TEST FAILED!!!\n";
        return 1;
    }

    std::cout << "TEST PASSED!!!\n";
    return 0;
}
