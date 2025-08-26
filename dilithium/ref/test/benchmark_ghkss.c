// #define _POSIX_C_SOURCE 200112L
#define _GNU_SOURCE
#include <time.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "bench.h"
#include "../params.h"
#include "ghkss_inputs.h"
#include "../poly256.h"
#include <assert.h>
#include "cpucycles.h"
#include "../GHKSS256.h"

static inline void keygen_wrapper(keygen_input_t *in)
{
    keygen_e_1d(in->A, in->B, in->S);
}

static inline void bench_keygen(void)
{
    keygen_input_t *in = malloc(sizeof(keygen_input_t));
    if (!in)
    {
        fprintf(stderr, "Memory allocation failed for keygen_input_t\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < 10; i++)
    {
        init_keygen_input(in);
        BENCH_ONCE_W_ARGS("KeyGen", keygen_wrapper, in);
        clear_keygen_input(in);
    }
    free(in);
}

// --- helper: monotonic time in seconds (double) ---
static inline double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

// Forward decl of the function under test (your implementation)
void as_keygen_1(
    poly (*As)[L],    // [K][L]
    poly (*Ae)[LHAT], // [KHAT][LHAT]
    poly (*Be)[M],    // [KHAT][M]
    poly *si,         // [L]
    poly *yi,         // [K]
    uint8_t h_yi[32],
    ctx_t *ctx);

// ---- Benchmark ----
void bench_as_keygen_1(int iters)
{
    if (iters <= 0)
        iters = 1;

    // ---------- Allocate / init matrices and vectors ----------
    poly(*As)[L] = (poly(*)[L])malloc(sizeof(poly[K][L]));
    poly(*Ae)[LHAT] = (poly(*)[LHAT])malloc(sizeof(poly[KHAT][LHAT]));
    poly(*Be)[M] = (poly(*)[M])malloc(sizeof(poly[KHAT][M]));
    poly *si = (poly *)malloc(sizeof(poly[L]));
    poly *yi = (poly *)malloc(sizeof(poly[K]));
    poly(*Se)[M] = (poly(*)[M])malloc(sizeof(poly[LHAT][M])); // temp for keygen_e_1d
    ctx_t ctx = {0};

    ctx.u = (poly *)malloc(sizeof(poly) * LHAT);
    ctx.v = (poly *)malloc(sizeof(poly) * M);

    if (!As || !Ae || !Be || !si || !yi || !Se)
    {
        fprintf(stderr, "Allocation failed\n");
        exit(1);
    }

    // Initialize poly objects
    for (int i = 0; i < K; i++)
        for (int j = 0; j < L; j++)
            poly_init(&As[i][j]);
    for (int i = 0; i < KHAT; i++)
        for (int j = 0; j < LHAT; j++)
            poly_init(&Ae[i][j]);
    for (int i = 0; i < KHAT; i++)
        for (int j = 0; j < M; j++)
            poly_init(&Be[i][j]);
    for (int i = 0; i < L; i++)
        poly_init(&si[i]);
    for (int i = 0; i < K; i++)
        poly_init(&yi[i]);
    for (int i = 0; i < LHAT; i++)
        for (int j = 0; j < M; j++)
            poly_init(&Se[i][j]);

    // // ctx.u, ctx.v expected to be arrays of poly with sizes [LHAT], [M].
    // // If ctx_t already stores them by value, just init them; otherwise allocate if they are pointers.
    for (int i = 0; i < LHAT; i++)
        poly_init(&ctx.u[i]);
    for (int j = 0; j < M; j++)
        poly_init(&ctx.v[j]);

    // ---------- Prepare public matrices ----------
    // As: sample uniformly and push to NTT once (si will be NTT in as_keygen_1)
    {
        uint8_t seed[SEEDBYTES];
        uint16_t nonce = 0;
        gen_randomness(seed);
        for (int i = 0; i < K; i++)
        {
            for (int j = 0; j < L; j++)
            {
                poly_uniform(&As[i][j], seed, nonce++);
                poly_ntt(&As[i][j]); // NTT domain to match pointwise multiply
            }
        }
    }

    // Ae, Be: generate using your earlier keygen (expects Se but we can discard Se afterwards)
    keygen_e_1d(Ae, Be, Se);

    // ---------- Benchmark loop ----------
    unsigned long long *cycles = (unsigned long long *)malloc(sizeof(unsigned long long) * iters);
    double *secs = (double *)malloc(sizeof(double) * iters);
    if (
        !cycles ||
        !secs)
    {
        fprintf(stderr, "Timing array alloc failed\n");
        exit(1);
    }

    // Warm-up (populate caches, JITs, etc.)
    {
        uint8_t h_yi[32] = {0};
        // Zero yi (as_keygen_1 accumulates into yi)
        for (int i = 0; i < K; i++)
            poly_zero(&yi[i]);
        as_keygen_1(As, Ae, Be, si, yi, h_yi, &ctx);
    }

    // Measured runs
    for (int t = 0; t < iters; t++)
    {
        // Reset yi to zero before each run (as_keygen_1 adds into it)
        for (int i = 0; i < K; i++)
            poly_zero(&yi[i]);

        uint8_t h_yi[32];

        unsigned long long before_c = cpucycles();
        double before_s = now_sec();

        as_keygen_1(As, Ae, Be, si, yi, h_yi, &ctx);

        double after_s = now_sec();
        unsigned long long after_c = cpucycles();
        cycles[t] = (after_c - before_c);
        secs[t] = (after_s - before_s);
    }

    // ---------- Stats ----------
    unsigned long long sum_c = 0ULL, min_c = (unsigned long long)-1, max_c = 0ULL;
    double sum_s = 0.0, min_s = 1e300, max_s = 0.0;

    for (int t = 0; t < iters; t++)
    {
        if (cycles[t] < min_c)
            min_c = cycles[t];
        if (cycles[t] > max_c)
            max_c = cycles[t];
        if (secs[t] < min_s)
            min_s = secs[t];
        if (secs[t] > max_s)
            max_s = secs[t];
        sum_c += cycles[t];
        sum_s += secs[t];
    }

    double avg_c = (double)sum_c / (double)iters;
    double avg_s = sum_s / (double)iters;

    // ---------- Print ----------
    for (int t = 0; t < iters; t++)
    {
        printf("RUN %2d: %12llu cycles, %8.6f s\n", t + 1, (unsigned long long)cycles[t], secs[t]);
    }
    printf("--------------------------------------------------\n");
    printf("as_keygen_1 AVG: %12.0f cycles  (min %llu, max %llu)\n", avg_c, (unsigned long long)min_c, (unsigned long long)max_c);
    printf("as_keygen_1 AVG: %8.6f s      (min %.6f, max %.6f)\n", avg_s, min_s, max_s);

    // ---------- Cleanup ----------
    free(cycles);

    for (int i = 0; i < K; i++)
        for (int j = 0; j < L; j++)
            poly_clear(&As[i][j]);
    for (int i = 0; i < KHAT; i++)
        for (int j = 0; j < LHAT; j++)
            poly_clear(&Ae[i][j]);
    for (int i = 0; i < KHAT; i++)
        for (int j = 0; j < M; j++)
            poly_clear(&Be[i][j]);
    for (int i = 0; i < L; i++)
        poly_clear(&si[i]);
    for (int i = 0; i < K; i++)
        poly_clear(&yi[i]);
    for (int i = 0; i < LHAT; i++)
        for (int j = 0; j < M; j++)
            poly_clear(&Se[i][j]);
    for (int i = 0; i < LHAT; i++)
        poly_clear(&ctx.u[i]);
    for (int j = 0; j < M; j++)
        poly_clear(&ctx.v[j]);

    free(As);
    free(Ae);
    free(Be);
    free(si);
    free(yi);
    free(Se);
    free(ctx.u);
    free(ctx.v);
}
void bench_as_keygen_preenc(int iters)
{
    if (iters <= 0)
        iters = 1;

    // Allocate & init As, si, yi
    poly(*As)[L] = (poly(*)[L])malloc(sizeof(poly[K][L]));
    poly *si = (poly *)malloc(sizeof(poly[L]));
    poly *yi = (poly *)malloc(sizeof(poly[K]));
    if (!As || !si || !yi)
    {
        fprintf(stderr, "alloc failed\n");
        exit(1);
    }

    for (int i = 0; i < K; i++)
        for (int j = 0; j < L; j++)
            poly_init(&As[i][j]);
    for (int i = 0; i < L; i++)
        poly_init(&si[i]);
    for (int i = 0; i < K; i++)
        poly_init(&yi[i]);

    // Prepare As in NTT
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    gen_randomness(seed);
    for (int i = 0; i < K; i++)
        for (int j = 0; j < L; j++)
        {
            poly_uniform(&As[i][j], seed, nonce++);
            poly_ntt(&As[i][j]);
        }

    unsigned long long *cycles = (unsigned long long *)malloc(sizeof(*cycles) * iters);
    double *secs = (double *)malloc(sizeof(*secs) * iters);
    if (!cycles || !secs)
    {
        fprintf(stderr, "timing alloc failed\n");
        exit(1);
    }

    // Warm-up
    {
        uint8_t h[32] = {0};
        for (int i = 0; i < K; i++)
            poly_zero(&yi[i]);
        as_keygen_preenc(As, si, yi, h);
    }

    // Measure
    for (int t = 0; t < iters; t++)
    {
        for (int i = 0; i < K; i++)
            poly_zero(&yi[i]);
        uint8_t h[32];
        double t0 = now_sec();
        unsigned long long c0 = cpucycles();
        as_keygen_preenc(As, si, yi, h);
        unsigned long long c1 = cpucycles();
        double t1 = now_sec();
        cycles[t] = c1 - c0;
        secs[t] = t1 - t0;
    }

    // Stats/print
    unsigned long long sumc = 0, minc = ~0ULL, maxc = 0;
    double sums = 0, mins = 1e300, maxs = 0;
    for (int t = 0; t < iters; t++)
    {
        if (cycles[t] < minc)
            minc = cycles[t];
        if (cycles[t] > maxc)
            maxc = cycles[t];
        if (secs[t] < mins)
            mins = secs[t];
        if (secs[t] > maxs)
            maxs = secs[t];
        sumc += cycles[t];
        sums += secs[t];
    }
    double avgc = (double)sumc / iters, avgs = sums / iters;
    for (int t = 0; t < iters; t++)
        printf("[preenc] RUN %2d: %12llu cycles, %8.6f s\n", t + 1, (unsigned long long)cycles[t], secs[t]);
    printf("--------------------------------------------------\n");
    printf("[preenc] AVG: %12.0f cycles (min %llu, max %llu)\n", avgc, (unsigned long long)minc, (unsigned long long)maxc);
    printf("[preenc] AVG: %8.6f s      (min %.6f, max %.6f)\n", avgs, mins, maxs);

    // Cleanup
    free(cycles);
    free(secs);
    for (int i = 0; i < K; i++)
        for (int j = 0; j < L; j++)
            poly_clear(&As[i][j]);
    for (int i = 0; i < L; i++)
        poly_clear(&si[i]);
    for (int i = 0; i < K; i++)
        poly_clear(&yi[i]);
    free(As);
    free(si);
    free(yi);
}

void bench_as_keygen_encrypt_only(int iters)
{
    if (iters <= 0)
        iters = 1;

    // Allocate & init Ae, Be, Se, si, ctx
    poly(*Ae)[LHAT] = (poly(*)[LHAT])malloc(sizeof(poly[KHAT][LHAT]));
    poly(*Be)[M] = (poly(*)[M])malloc(sizeof(poly[KHAT][M]));
    poly(*Se)[M] = (poly(*)[M])malloc(sizeof(poly[LHAT][M]));
    poly *si = (poly *)malloc(sizeof(poly[L]));
    ctx_t ctx = (ctx_t){0};
    ctx.u = (poly *)malloc(sizeof(poly) * LHAT);
    ctx.v = (poly *)malloc(sizeof(poly) * M);
    if (!Ae || !Be || !Se || !si || !ctx.u || !ctx.v)
    {
        fprintf(stderr, "alloc failed\n");
        exit(1);
    }

    for (int i = 0; i < KHAT; i++)
        for (int j = 0; j < LHAT; j++)
            poly_init(&Ae[i][j]);
    for (int i = 0; i < KHAT; i++)
        for (int j = 0; j < M; j++)
            poly_init(&Be[i][j]);
    for (int i = 0; i < LHAT; i++)
        for (int j = 0; j < M; j++)
            poly_init(&Se[i][j]);
    for (int i = 0; i < L; i++)
        poly_init(&si[i]);
    for (int i = 0; i < LHAT; i++)
        poly_init(&ctx.u[i]);
    for (int j = 0; j < M; j++)
        poly_init(&ctx.v[j]);

    // Generate (Ae,Be)
    keygen_e_1d(Ae, Be, Se);

    // Prepare si (NTT) – treat as consumed by encrypt_1d
    for (int i = 0; i < L; i++)
    {
        gaussian_sampler_y(si[i].coeffs, N);
        poly_reduce(&si[i]);
        poly_ntt(&si[i]);
    }

    unsigned long long *cycles = (unsigned long long *)malloc(sizeof(*cycles) * iters);
    double *secs = (double *)malloc(sizeof(*secs) * iters);
    if (!cycles || !secs)
    {
        fprintf(stderr, "timing alloc failed\n");
        exit(1);
    }

    // Warm-up
    as_keygen_encrypt_only(&ctx, Ae, Be, si);

    // IMPORTANT: if encrypt_1d mutates si/ctx.u/ctx.v, re-init them each iteration.
    for (int t = 0; t < iters; t++)
    {
        // refresh si in case encrypt_1d consumed it
        for (int i = 0; i < L; i++)
        {
            poly_zero(&si[i]);
            gaussian_sampler_y(si[i].coeffs, N);
            poly_reduce(&si[i]);
            poly_ntt(&si[i]);
        }
        // refresh ctx.u/ctx.v if your encrypt mutates them (cheap and safe)
        for (int i = 0; i < LHAT; i++)
            poly_zero(&ctx.u[i]);
        for (int j = 0; j < M; j++)
            poly_zero(&ctx.v[j]);

        double t0 = now_sec();
        unsigned long long c0 = cpucycles();
        as_keygen_encrypt_only(&ctx, Ae, Be, si);
        unsigned long long c1 = cpucycles();
        double t1 = now_sec();

        cycles[t] = c1 - c0;
        secs[t] = t1 - t0;
    }

    // Stats/print
    unsigned long long sumc = 0, minc = ~0ULL, maxc = 0;
    double sums = 0, mins = 1e300, maxs = 0;
    for (int t = 0; t < iters; t++)
    {
        if (cycles[t] < minc)
            minc = cycles[t];
        if (cycles[t] > maxc)
            maxc = cycles[t];
        if (secs[t] < mins)
            mins = secs[t];
        if (secs[t] > maxs)
            maxs = secs[t];
        sumc += cycles[t];
        sums += secs[t];
    }
    double avgc = (double)sumc / iters, avgs = sums / iters;
    for (int t = 0; t < iters; t++)
        printf("[encrypt] RUN %2d: %12llu cycles, %8.6f s\n", t + 1, (unsigned long long)cycles[t], secs[t]);
    printf("--------------------------------------------------\n");
    printf("[encrypt] AVG: %12.0f cycles (min %llu, max %llu)\n", avgc, (unsigned long long)minc, (unsigned long long)maxc);
    printf("[encrypt] AVG: %8.6f s      (min %.6f, max %.6f)\n", avgs, mins, maxs);

    // Cleanup
    free(cycles);
    free(secs);
    for (int i = 0; i < KHAT; i++)
        for (int j = 0; j < LHAT; j++)
            poly_clear(&Ae[i][j]);
    for (int i = 0; i < KHAT; i++)
        for (int j = 0; j < M; j++)
            poly_clear(&Be[i][j]);
    for (int i = 0; i < LHAT; i++)
        for (int j = 0; j < M; j++)
            poly_clear(&Se[i][j]);
    for (int i = 0; i < L; i++)
        poly_clear(&si[i]);
    for (int i = 0; i < LHAT; i++)
        poly_clear(&ctx.u[i]);
    for (int j = 0; j < M; j++)
        poly_clear(&ctx.v[j]);
    free(Ae);
    free(Be);
    free(Se);
    free(si);
    free(ctx.u);
    free(ctx.v);
}

extern void as_keygen_2(ctx_t *ctx_s, poly (*yi)[K], uint8_t (*h_yi)[32], ctx_t *ctx, poly *yprime);

// --- helpers to init/clear ctx_t arrays (pointer-owning version) ---
static void ctx_alloc_init(ctx_t *c)
{
    c->u = (poly *)malloc(sizeof(poly) * LHAT);
    c->v = (poly *)malloc(sizeof(poly) * M);
    if (!c->u || !c->v)
    {
        fprintf(stderr, "ctx alloc failed\n");
        exit(1);
    }
    for (int i = 0; i < LHAT; i++)
        poly_init(&c->u[i]);
    for (int j = 0; j < M; j++)
        poly_init(&c->v[j]);
}
static void ctx_clear_free(ctx_t *c)
{
    for (int i = 0; i < LHAT; i++)
        poly_clear(&c->u[i]);
    for (int j = 0; j < M; j++)
        poly_clear(&c->v[j]);
    free(c->u);
    free(c->v);
}

// --- fill helpers (time-domain polys) ---
static void sample_poly_time(poly *p)
{
    gaussian_sampler_y(p->coeffs, N);
    poly_reduce(p);
}
static void zero_ctx(ctx_t *c)
{
    for (int i = 0; i < LHAT; i++)
        poly_zero(&c->u[i]);
    for (int j = 0; j < M; j++)
        poly_zero(&c->v[j]);
}

// Precompute h_yi to MATCH H0(yi[i]) (keeps control-flow stable)
static void fill_yi_and_hash(poly (*yi)[K], uint8_t (*h_yi)[32])
{
    for (int u = 0; u < USERS; u++)
    {
        for (int k = 0; k < K; k++)
        {
            poly_init(&yi[u][k]);
            sample_poly_time(&yi[u][k]);
        }
        Hash(yi[u], h_yi[u]);
    }
}

void bench_as_keygen_2(int iters)
{
    if (iters <= 0)
        iters = 1;

    // --- allocate / init inputs ---
    ctx_t ctx_s = {0};
    ctx_alloc_init(&ctx_s);
    zero_ctx(&ctx_s);

    ctx_t *ctx = (ctx_t *)malloc(sizeof(ctx_t) * USERS);
    if (!ctx)
    {
        fprintf(stderr, "alloc ctx array failed\n");
        exit(1);
    }
    for (int u = 0; u < USERS; u++)
    {
        ctx[u].u = NULL;
        ctx[u].v = NULL;
        ctx_alloc_init(&ctx[u]);
    }

    poly(*yi)[K] = (poly(*)[K])malloc(sizeof(poly[USERS][K]));
    uint8_t (*h_yi)[32] = (uint8_t (*)[32])malloc(sizeof(uint8_t[USERS][32]));
    poly *yprime = (poly *)malloc(sizeof(poly) * K);
    if (!yi || !h_yi || !yprime)
    {
        fprintf(stderr, "alloc yi/h_yi/yprime failed\n");
        exit(1);
    }

    // init yprime polys
    for (int k = 0; k < K; k++)
        poly_init(&yprime[k]);

    // init yi polys + hashes (and also put some data in ctx[u])
    for (int u = 0; u < USERS; u++)
    {
        for (int j = 0; j < LHAT; j++)
            sample_poly_time(&ctx[u].u[j]);
        for (int j = 0; j < M; j++)
            sample_poly_time(&ctx[u].v[j]);
    }
    fill_yi_and_hash(yi, h_yi);

    // --- timing arrays ---
    unsigned long long *cycles = (unsigned long long *)malloc(sizeof(*cycles) * iters);
    double *secs = (double *)malloc(sizeof(*secs) * iters);
    if (!cycles || !secs)
    {
        fprintf(stderr, "timing alloc failed\n");
        exit(1);
    }

    // --- warm-up ---
    {
        // reset ctx_s and yprime (function accumulates into them)
        zero_ctx(&ctx_s);
        for (int k = 0; k < K; k++)
            poly_zero(&yprime[k]);
        as_keygen_2(&ctx_s, yi, h_yi, ctx, yprime);
    }

    // --- measured runs ---
    for (int t = 0; t < iters; t++)
    {
        // refresh accumulators
        zero_ctx(&ctx_s);
        for (int k = 0; k < K; k++)
            poly_zero(&yprime[k]);

        double t0 = now_sec();
        unsigned long long c0 = cpucycles();

        as_keygen_2(&ctx_s, yi, h_yi, ctx, yprime);

        unsigned long long c1 = cpucycles();
        double t1 = now_sec();

        cycles[t] = c1 - c0;
        secs[t] = t1 - t0;
    }

    // --- stats ---
    unsigned long long sumc = 0, minc = ~0ULL, maxc = 0;
    double sums = 0.0, mins = 1e300, maxs = 0.0;
    for (int t = 0; t < iters; t++)
    {
        if (cycles[t] < minc)
            minc = cycles[t];
        if (cycles[t] > maxc)
            maxc = cycles[t];
        if (secs[t] < mins)
            mins = secs[t];
        if (secs[t] > maxs)
            maxs = secs[t];
        sumc += cycles[t];
        sums += secs[t];
    }
    double avgc = (double)sumc / iters;
    double avgs = sums / iters;

    // --- print ---
    for (int t = 0; t < iters; t++)
        printf("[as_keygen_2] RUN %2d: %12llu cycles, %8.6f s\n",
               t + 1, (unsigned long long)cycles[t], secs[t]);
    printf("--------------------------------------------------\n");
    printf("[as_keygen_2] AVG: %12.0f cycles (min %llu, max %llu)\n",
           avgc, (unsigned long long)minc, (unsigned long long)maxc);
    printf("[as_keygen_2] AVG: %8.6f s      (min %.6f, max %.6f)\n",
           avgs, mins, maxs);

    // --- cleanup ---
    free(cycles);
    free(secs);

    for (int u = 0; u < USERS; u++)
    {
        for (int k = 0; k < K; k++)
            poly_clear(&yi[u][k]);
    }
    free(yi);
    free(h_yi);

    for (int k = 0; k < K; k++)
        poly_clear(&yprime[k]);
    free(yprime);

    for (int u = 0; u < USERS; u++)
        ctx_clear_free(&ctx[u]);
    free(ctx);
    ctx_clear_free(&ctx_s);
}

static poly *make_r_ntt(void)
{
    poly *r = (poly *)malloc(sizeof(poly) * M);
    if (!r)
    {
        fprintf(stderr, "alloc r failed\n");
        exit(1);
    }
    poly_1d_init(r, M);
    for (int i = 0; i < L; i++)
    {
        gaussian_sampler_w(r[i].coeffs, N);
        poly_reduce(&r[i]);
        poly_ntt(&r[i]);
    }
    for (int i = L; i < M; i++)
    {
        poly_zero(&r[i]);
        poly_ntt(&r[i]);
    }
    return r;
}

void bench_as_sign_r1(int iters)
{
    if (iters <= 0)
        iters = 1;

    // Allocate & init As
    poly(*As)[L] = (poly(*)[L])malloc(sizeof(poly[K][L]));
    if (!As)
    {
        fprintf(stderr, "alloc As failed\n");
        exit(1);
    }
    for (int i = 0; i < K; i++)
        for (int j = 0; j < L; j++)
            poly_init(&As[i][j]);

    // Prepare As in NTT
    {
        uint8_t seed[SEEDBYTES];
        uint16_t nonce = 0;
        gen_randomness(seed);
        for (int i = 0; i < K; i++)
            for (int j = 0; j < L; j++)
            {
                poly_uniform(&As[i][j], seed, nonce++);
                poly_ntt(&As[i][j]);
            }
    }

    unsigned long long *cycles = (unsigned long long *)malloc(sizeof(*cycles) * iters);
    double *secs = (double *)malloc(sizeof(*secs) * iters);
    if (!cycles || !secs)
    {
        fprintf(stderr, "timing alloc failed\n");
        exit(1);
    }

    // Warm-up
    {
        uint8_t h[32];
        poly *r = NULL;
        as_sign_r1(As, h, &r);
        poly_1d_clear(r, M);
        free(r);
    }

    // Measure
    for (int t = 0; t < iters; t++)
    {
        uint8_t h[32];
        poly *r = NULL;

        double t0 = now_sec();
        unsigned long long c0 = cpucycles();
        as_sign_r1(As, h, &r);
        unsigned long long c1 = cpucycles();
        double t1 = now_sec();

        cycles[t] = c1 - c0;
        secs[t] = t1 - t0;

        // cleanup r from this run
        poly_1d_clear(r, M);
        free(r);
    }

    // Stats/print
    unsigned long long sumc = 0, minc = ~0ULL, maxc = 0;
    double sums = 0, mins = 1e300, maxs = 0;
    for (int t = 0; t < iters; t++)
    {
        if (cycles[t] < minc)
            minc = cycles[t];
        if (cycles[t] > maxc)
            maxc = cycles[t];
        if (secs[t] < mins)
            mins = secs[t];
        if (secs[t] > maxs)
            maxs = secs[t];
        sumc += cycles[t];
        sums += secs[t];
    }
    double avgc = (double)sumc / iters, avgs = sums / iters;

    for (int t = 0; t < iters; t++)
        printf("[as_sign_r1] RUN %2d: %12llu cycles, %8.6f s\n", t + 1, (unsigned long long)cycles[t], secs[t]);
    printf("--------------------------------------------------\n");
    printf("[as_sign_r1] AVG: %12.0f cycles (min %llu, max %llu)\n", avgc, (unsigned long long)minc, (unsigned long long)maxc);
    printf("[as_sign_r1] AVG: %8.6f s      (min %.6f, max %.6f)\n", avgs, mins, maxs);

    // Cleanup
    free(cycles);
    free(secs);
    for (int i = 0; i < K; i++)
        for (int j = 0; j < L; j++)
            poly_clear(&As[i][j]);
    free(As);
}

void bench_as_sign_r2(int iters)
{
    if (iters <= 0)
        iters = 1;

    // Allocate & init Ae, Be, Se (for keygen_e_1d), u, v
    poly(*Ae)[LHAT] = (poly(*)[LHAT])malloc(sizeof(poly[KHAT][LHAT]));
    poly(*Be)[M] = (poly(*)[M])malloc(sizeof(poly[KHAT][M]));
    poly(*Se)[M] = (poly(*)[M])malloc(sizeof(poly[LHAT][M]));
    poly *u = (poly *)malloc(sizeof(poly) * LHAT);
    poly *v = (poly *)malloc(sizeof(poly) * M);
    if (!Ae || !Be || !Se || !u || !v)
    {
        fprintf(stderr, "alloc failed\n");
        exit(1);
    }

    for (int i = 0; i < KHAT; i++)
        for (int j = 0; j < LHAT; j++)
            poly_init(&Ae[i][j]);
    for (int i = 0; i < KHAT; i++)
        for (int j = 0; j < M; j++)
            poly_init(&Be[i][j]);
    for (int i = 0; i < LHAT; i++)
        for (int j = 0; j < M; j++)
            poly_init(&Se[i][j]);
    for (int i = 0; i < LHAT; i++)
        poly_init(&u[i]);
    for (int j = 0; j < M; j++)
        poly_init(&v[j]);

    // Generate (Ae,Be)
    keygen_e_1d(Ae, Be, Se);

    unsigned long long *cycles = (unsigned long long *)malloc(sizeof(*cycles) * iters);
    double *secs = (double *)malloc(sizeof(*secs) * iters);
    if (!cycles || !secs)
    {
        fprintf(stderr, "timing alloc failed\n");
        exit(1);
    }

    // Warm-up
    {
        for (int i = 0; i < LHAT; i++)
            poly_zero(&u[i]);
        for (int j = 0; j < M; j++)
            poly_zero(&v[j]);
        poly *r = make_r_ntt();
        as_sign_r2(Ae, Be, u, v, r, 1);
    }

    // Measure
    for (int t = 0; t < iters; t++)
    {
        // fresh u,v (in case encrypt_1d writes to them) and fresh r
        for (int i = 0; i < LHAT; i++)
            poly_zero(&u[i]);
        for (int j = 0; j < M; j++)
            poly_zero(&v[j]);
        poly *r = make_r_ntt();

        double t0 = now_sec();
        unsigned long long c0 = cpucycles();
        as_sign_r2(Ae, Be, u, v, r, 1); // frees r
        unsigned long long c1 = cpucycles();
        double t1 = now_sec();

        cycles[t] = c1 - c0;
        secs[t] = t1 - t0;
    }

    // Stats/print
    unsigned long long sumc = 0, minc = ~0ULL, maxc = 0;
    double sums = 0, mins = 1e300, maxs = 0;
    for (int t = 0; t < iters; t++)
    {
        if (cycles[t] < minc)
            minc = cycles[t];
        if (cycles[t] > maxc)
            maxc = cycles[t];
        if (secs[t] < mins)
            mins = secs[t];
        if (secs[t] > maxs)
            maxs = secs[t];
        sumc += cycles[t];
        sums += secs[t];
    }
    double avgc = (double)sumc / iters, avgs = sums / iters;

    for (int t = 0; t < iters; t++)
        printf("[as_sign_r2] RUN %2d: %12llu cycles, %8.6f s\n", t + 1, (unsigned long long)cycles[t], secs[t]);
    printf("--------------------------------------------------\n");
    printf("[as_sign_r2] AVG: %12.0f cycles (min %llu, max %llu)\n", avgc, (unsigned long long)minc, (unsigned long long)maxc);
    printf("[as_sign_r2] AVG: %8.6f s      (min %.6f, max %.6f)\n", avgs, mins, maxs);

    // Cleanup
    free(cycles);
    free(secs);
    for (int i = 0; i < KHAT; i++)
        for (int j = 0; j < LHAT; j++)
            poly_clear(&Ae[i][j]);
    for (int i = 0; i < KHAT; i++)
        for (int j = 0; j < M; j++)
            poly_clear(&Be[i][j]);
    for (int i = 0; i < LHAT; i++)
        for (int j = 0; j < M; j++)
            poly_clear(&Se[i][j]);
    for (int i = 0; i < LHAT; i++)
        poly_clear(&u[i]);
    for (int j = 0; j < M; j++)
        poly_clear(&v[j]);
    free(Ae);
    free(Be);
    free(Se);
    free(u);
    free(v);
}

void bench_as_sign_r3(int iters)
{
    if (iters <= 0)
        iters = 1;

    // Allocate
    poly(*A)[L] = (poly(*)[L])malloc(sizeof(poly[K][L]));
    poly *yprime = (poly *)malloc(sizeof(poly) * K);
    poly(*u_r)[LHAT] = (poly(*)[LHAT])malloc(sizeof(poly[THRESHOLD][LHAT]));
    poly(*v_r)[M] = (poly(*)[M])malloc(sizeof(poly[THRESHOLD][M]));
    poly *u_s = (poly *)malloc(sizeof(poly) * LHAT);
    poly *v_s = (poly *)malloc(sizeof(poly) * M);
    poly *u_z = (poly *)malloc(sizeof(poly) * LHAT);
    poly *v_z = (poly *)malloc(sizeof(poly) * M);
    poly(*ski)[M] = (poly(*)[M])malloc(sizeof(poly[LHAT][M]));
    poly *mu = (poly *)malloc(sizeof(poly) * M);  // “L but padded to M” → allocate M
    poly *dsi = (poly *)malloc(sizeof(poly) * M); // output buffer, size M for safety
    poly(*w_i)[K] = (poly(*)[K])malloc(sizeof(poly[THRESHOLD][K]));
    uint8_t **h_wi = (uint8_t **)malloc(sizeof(uint8_t *) * THRESHOLD);
    poly *c = (poly *)malloc(sizeof(poly) * 1);
    poly *z = (poly *)malloc(sizeof(poly) * L);
    poly *hvec = (poly *)malloc(sizeof(poly) * K);
    int32_t *users = (int32_t *)malloc(sizeof(int32_t) * THRESHOLD);
    if (!A || !yprime || !u_r || !v_r || !u_s || !v_s || !u_z || !v_z ||
        !ski || !mu || !dsi || !w_i || !h_wi || !c || !z || !hvec || !users)
    {
        fprintf(stderr, "alloc failed\n");
        exit(1);
    }
    for (int i = 0; i < THRESHOLD; i++)
    {
        h_wi[i] = (uint8_t *)malloc(32);
        if (!h_wi[i])
        {
            fprintf(stderr, "alloc h_wi\n");
            exit(1);
        }
    }

    // Init polys
    for (int i = 0; i < K; i++)
        for (int j = 0; j < L; j++)
            poly_init(&A[i][j]);
    for (int i = 0; i < K; i++)
        poly_init(&yprime[i]);
    for (int t = 0; t < THRESHOLD; t++)
        for (int j = 0; j < LHAT; j++)
            poly_init(&u_r[t][j]);
    for (int t = 0; t < THRESHOLD; t++)
        for (int j = 0; j < M; j++)
            poly_init(&v_r[t][j]);
    for (int j = 0; j < LHAT; j++)
        poly_init(&u_s[j]), poly_init(&u_z[j]);
    for (int j = 0; j < M; j++)
        poly_init(&v_s[j]), poly_init(&v_z[j]);
    for (int i = 0; i < LHAT; i++)
        for (int j = 0; j < M; j++)
            poly_init(&ski[i][j]);
    for (int j = 0; j < M; j++)
        poly_init(&mu[j]), poly_init(&dsi[j]);
    for (int t = 0; t < THRESHOLD; t++)
        for (int k = 0; k < K; k++)
            poly_init(&w_i[t][k]);
    poly_init(&c[0]);
    for (int i = 0; i < L; i++)
        poly_init(&z[i]);
    for (int i = 0; i < K; i++)
        poly_init(&hvec[i]);

    // Fill inputs with time-domain samples; zero accumulators
    for (int i = 0; i < K; i++)
        for (int j = 0; j < L; j++)
            sample_poly_time(&A[i][j]);
    for (int i = 0; i < K; i++)
        sample_poly_time(&yprime[i]);
    for (int t = 0; t < THRESHOLD; t++)
        for (int j = 0; j < LHAT; j++)
            sample_poly_time(&u_r[t][j]);
    for (int t = 0; t < THRESHOLD; t++)
        for (int j = 0; j < M; j++)
            sample_poly_time(&v_r[t][j]);
    for (int j = 0; j < LHAT; j++)
        sample_poly_time(&u_s[j]), poly_zero(&u_z[j]);
    for (int j = 0; j < M; j++)
        sample_poly_time(&v_s[j]), poly_zero(&v_z[j]);
    for (int i = 0; i < LHAT; i++)
        for (int j = 0; j < M; j++)
            sample_poly_time(&ski[i][j]);
    // mu: “L but padded to M” → fill first L, zero rest
    for (int j = 0; j < L; j++)
        sample_poly_time(&mu[j]);
    for (int j = L; j < M; j++)
        poly_zero(&mu[j]);
    // dsi output buffer zeroed
    for (int j = 0; j < M; j++)
        poly_zero(&dsi[j]);
    // w_i inputs
    for (int t = 0; t < THRESHOLD; t++)
        for (int k = 0; k < K; k++)
            sample_poly_time(&w_i[t][k]);

    // h_wi and users
    for (int i = 0; i < THRESHOLD; i++)
    {
        memset(h_wi[i], 0, 32);
        users[i] = i;
    }
    int user = 0; // benchmark with first user

    // timing arrays
    unsigned long long *cycles = (unsigned long long *)malloc(sizeof(*cycles) * iters);
    double *secs = (double *)malloc(sizeof(*secs) * iters);
    if (!cycles || !secs)
    {
        fprintf(stderr, "timing alloc failed\n");
        exit(1);
    }

    // warm-up (reset accumulators)
    for (int j = 0; j < LHAT; j++)
        poly_zero(&u_z[j]);
    for (int j = 0; j < M; j++)
        poly_zero(&v_z[j]);
    for (int j = 0; j < M; j++)
        poly_zero(&dsi[j]);
    as_sign_2(c, z, hvec, A, yprime, u_r, v_r, u_s, v_s, u_z, v_z, ski, mu, dsi, w_i, h_wi, user, users);

    // measured runs
    for (int t = 0; t < iters; t++)
    {
        // refresh accumulators (u_z, v_z, dsi) and outputs (c,z,h)
        for (int j = 0; j < LHAT; j++)
            poly_zero(&u_z[j]);
        for (int j = 0; j < M; j++)
            poly_zero(&v_z[j]);
        for (int j = 0; j < M; j++)
            poly_zero(&dsi[j]);
        poly_zero(&c[0]);
        for (int i = 0; i < L; i++)
            poly_zero(&z[i]);
        for (int i = 0; i < K; i++)
            poly_zero(&hvec[i]);

        double t0 = now_sec();
        unsigned long long c0 = cpucycles();

        as_sign_2(c, z, hvec, A, yprime, u_r, v_r, u_s, v_s, u_z, v_z, ski, mu, dsi, w_i, h_wi, user, users);

        unsigned long long c1 = cpucycles();
        double t1 = now_sec();

        cycles[t] = c1 - c0;
        secs[t] = t1 - t0;
    }

    // stats/print
    unsigned long long sumc = 0, minc = ~0ULL, maxc = 0;
    double sums = 0, mins = 1e300, maxs = 0;
    for (int t = 0; t < iters; t++)
    {
        if (cycles[t] < minc)
            minc = cycles[t];
        if (cycles[t] > maxc)
            maxc = cycles[t];
        if (secs[t] < mins)
            mins = secs[t];
        if (secs[t] > maxs)
            maxs = secs[t];
        sumc += cycles[t];
        sums += secs[t];
    }
    double avgc = (double)sumc / iters, avgs = sums / iters;

    for (int t = 0; t < iters; t++)
        printf("[as_sign_2] RUN %2d: %12llu cycles, %8.6f s\n",
               t + 1, (unsigned long long)cycles[t], secs[t]);
    printf("--------------------------------------------------\n");
    printf("[as_sign_2] AVG: %12.0f cycles (min %llu, max %llu)\n",
           avgc, (unsigned long long)minc, (unsigned long long)maxc);
    printf("[as_sign_2] AVG: %8.6f s      (min %.6f, max %.6f)\n",
           avgs, mins, maxs);

    // cleanup
    free(cycles);
    free(secs);
    for (int i = 0; i < K; i++)
        for (int j = 0; j < L; j++)
            poly_clear(&A[i][j]);
    for (int i = 0; i < K; i++)
        poly_clear(&yprime[i]);
    for (int t = 0; t < THRESHOLD; t++)
        for (int j = 0; j < LHAT; j++)
            poly_clear(&u_r[t][j]);
    for (int t = 0; t < THRESHOLD; t++)
        for (int j = 0; j < M; j++)
            poly_clear(&v_r[t][j]);
    for (int j = 0; j < LHAT; j++)
        poly_clear(&u_s[j]), poly_clear(&u_z[j]);
    for (int j = 0; j < M; j++)
        poly_clear(&v_s[j]), poly_clear(&v_z[j]);
    for (int i = 0; i < LHAT; i++)
        for (int j = 0; j < M; j++)
            poly_clear(&ski[i][j]);
    for (int j = 0; j < M; j++)
        poly_clear(&mu[j]), poly_clear(&dsi[j]);
    for (int t = 0; t < THRESHOLD; t++)
        for (int k = 0; k < K; k++)
            poly_clear(&w_i[t][k]);
    poly_clear(&c[0]);
    for (int i = 0; i < L; i++)
        poly_clear(&z[i]);
    for (int i = 0; i < K; i++)
        poly_clear(&hvec[i]);

    for (int i = 0; i < THRESHOLD; i++)
        free(h_wi[i]);
    free(h_wi);

    free(A);
    free(yprime);
    free(u_r);
    free(v_r);
    free(u_s);
    free(v_s);
    free(u_z);
    free(v_z);
    free(ski);
    free(mu);
    free(dsi);
    free(w_i);
    free(c);
    free(z);
    free(hvec);
    free(users);
}

static inline void dkg1_wrapper(dkg_input_t *in)
{
    dk_gen_1(in->Ae, in->Si, in->Ei, in->Bi, in->h_Bi);
}

static inline void dkg2_wrapper(dkg_input_t *in)
{
    dk_gen_2(in->Si, in->Ei, in->Ae, in->users);
}

static inline void dkg3_wrapper(dkg_input_t *in)
{
    dk_gen_3(in->Bj, in->Sij, in->Be, in->Sei);
}

static inline void bench_dkg(void)
{
    dkg_input_t *in = malloc(sizeof(dkg_input_t));
    init_dkg(in);
    BENCH_MANY_W_ARGS("DKG round 1", dkg1_wrapper, in);
    BENCH_MANY_W_ARGS("DKG round 2", dkg2_wrapper, in);
    BENCH_MANY_W_ARGS("DKG round 3", dkg3_wrapper, in);
    clear_dkg(in);
    free(in);
}

static inline void sig1_wrapper(sig_input_t *in)
{
    as_sign_1(in->A_pks, in->Ae, in->Be, in->u, in->v);
}

static inline void sig2_wrapper(sig_input_t *in)
{
    as_sign_2(in->c_sig, in->z_sig, in->h_sig, in->A_pks, in->yprime, in->u_r, in->v_r, in->u_s, in->v_s, in->u_z, in->v_z, in->ski, in->mu, in->dsi, in->w_i, in->h_wi, in->user, in->users);
}

static inline void sig3_wrapper(sig_input_t *in)
{
    // as_sign_3(in->u_z, in->v_z, in->A_pks, in->c_sig, in->z_sig, in->h_sig, in->yprime, in->wprime, in->dsi, in->users);
}

static inline void verify_wrapper(sig_input_t *in)
{
}

static inline void bench_sign(void)
{
    sig_input_t *in = malloc(sizeof(sig_input_t));
    // init_sig(in);
    // BENCH_ONCE_W_ARGS("Sig round 1", sig1_wrapper, in);
    // BENCH_ONCE_W_ARGS("Sig round 2", sig2_wrapper, in);
    // BENCH_ONCE_W_ARGS("Sig round 3", sig3_wrapper, in);
    // clear_sig(in);
    free(in);
}

int main(void)
{
    init_params();
    init_zetas();

    // bench_dkg();
    // bench_keygen();
    // bench_as_keygen_1(2);
    // bench_as_keygen_preenc(10);
    // bench_as_keygen_encrypt_only(10);
    // bench_as_keygen_2(10);
    // bench_as_sign_r1(1);
    // bench_as_sign_r2(1);
    bench_as_sign_r3(1);

    free_params();
    clear_zetas();
    return 0;
}
