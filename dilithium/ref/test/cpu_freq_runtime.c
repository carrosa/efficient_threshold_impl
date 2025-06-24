// #define _POSIX_C_SOURCE 200112L
#define _GNU_SOURCE
#include "cpu_freq_runtime.h"
#include <stdint.h>
#include <time.h>
#include <stdio.h>

static inline uint64_t rdtsc(void)
{
    unsigned int lo, hi;
    __asm__ volatile("rdtsc" : "=a"(lo), "=d"(hi));
    return ((uint64_t)hi << 32) | lo;
}

uint64_t calibrate_cpu_freq_hz(void)
{
    struct timespec start_ts, end_ts;
    uint64_t start_cycles, end_cycles;
    // Burn CPU cycles before calibration to "wake up" CPU
    for (volatile int i = 0; i < 100000000; i++)
    {
        __asm__ volatile("");
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &start_ts);
    start_cycles = rdtsc();
    struct timespec delay = {0, 100000000};
    nanosleep(&delay, NULL);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end_ts);
    end_cycles = rdtsc();

    double elapsed_sec = (end_ts.tv_sec - start_ts.tv_sec) + (end_ts.tv_nsec - start_ts.tv_nsec) / 1e9;
    return (uint64_t)((end_cycles - start_cycles) / elapsed_sec);
}