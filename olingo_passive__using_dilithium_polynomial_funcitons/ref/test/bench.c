#include "bench.h"
#include "cpucycles.h"

static uint64_t before;
static uint64_t after;
static uint64_t total;
static uint64_t overhead;
static uint64_t grouping_total;

void bench_reset(void)
{
    total = 0;
    overhead = cpucycles_overhead();
}

void bench_before(void)
{
    before = cpucycles();
}

void bench_after(void)
{
    after = cpucycles();
    uint64_t result = after-before;
    total += (result >overhead ? result - overhead : 0);
}

void bench_compute(int benches)
{
    if (benches > 0)
    {
        total /= benches;
    }
}

void bench_print(void)
{
    double time_sec = (double)total / CPU_FREQ;
    printf("%llu cycles (%.6f seconds)\n\n", total, time_sec);
}


unsigned long long bench_get_total(void)
{
    return total;
}