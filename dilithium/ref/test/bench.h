#ifndef BENCH_H
#define BENCH_H

#include <stdint.h>
#include <stdio.h>
#include "cpu_freq_runtime.h"

// Helper for generating unique names
#define CONCAT_(x, y) x##y
#define CONCAT(x, y) CONCAT_(x, y)

#define CPU_FREQ 2600000000 // in Hz (Caroline PC, comment out if measuring runtime freq)

#ifndef CPU_FREQ
#define CPU_FREQ calibrate_cpu_freq_hz()
#endif
// #define CPU_FREQ 2294608000 // in Hz (NTNU Server)

// How many times to run the benchmarks
#define BENCH 10

#define BENCH_ONCE(LABEL, FUNC)             \
    bench_reset();                          \
    printf("BENCH: " LABEL "%*c = ",        \
           (int)(32 - strlen(LABEL)), ' '); \
    bench_before();                         \
    FUNC;                                   \
    bench_after();                          \
    bench_compute(1);                       \
    bench_print();

#define BENCH_ONCE_W_ARGS(LABEL, FUNC, ...) \
    printf("CPU_FREQ not defined, using runtime calibration: %llu Hz\n", CPU_FREQ); \
    BENCH_ONCE(LABEL, FUNC(__VA_ARGS__))

#define BENCH_MANY(LABEL, FUNC)             \
    bench_reset();                          \
    printf("BENCH: " LABEL "%*c = ",        \
           (int)(32 - strlen(LABEL)), ' '); \
    bench_before();                         \
    for (int i = 0; i < BENCH; i++)         \
    {                                       \
        FUNC;                               \
    }                                       \
    bench_after();                          \
    bench_compute(BENCH);                   \
    bench_print();

#define BENCH_ONCE_WRAP(LABEL, FUNC, ...)           \
    static inline void __bench_wrapper_##FUNC(void) \
    {                                               \
        FUNC(__VA_ARGS__);                          \
    }                                               \
    BENCH_ONCE(LABEL, __bench_wrapper_##FUNC())

#define BENCH_MANY_WRAP(LABEL, FUNC, ...)           \
    static inline void __bench_wrapper_##FUNC(void) \
    {                                               \
        FUNC(__VA_ARGS__);                          \
    }                                               \
    BENCH_MANY(LABEL, __bench_wrapper_##FUNC())

#define BENCH_BLOCK_ONCE(LABEL, BLOCK) \
    static inline void CONCAT(__bench_block_wrapper_, __LINE__)(void){BLOCK} BENCH_ONCE(LABEL, CONCAT(__bench_block_wrapper_, __LINE__)())

#define BENCH_BLOCK_MANY(LABEL, BLOCK) \
    static inline void CONCAT(__bench_block_wrapper_, __LINE__)(void){BLOCK} BENCH_MANY(LABEL, CONCAT(__bench_block_wrapper_, __LINE__)())

void bench_reset(void);
void bench_before(void);
void bench_after(void);
void bench_compute(int benches);
void bench_print(void);
unsigned long long bench_get_total(void);

#endif
