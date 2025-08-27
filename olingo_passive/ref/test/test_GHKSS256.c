#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "../params.h"
#include "../poly256.h"
#include "../GHKSS256.h"
#include "../gaussian_ref.h"
#include <assert.h>

#define CPU_FREQ 2600000000 // in Hz



void print_timing(uint64_t start, uint64_t end, const char *label)
{
    uint64_t cycles = end - start;
    double time_in_seconds = (double)cycles / CPU_FREQ;
    double time_in_ms = time_in_seconds * 1000;
    // printf("%s: %f seconds\n", label, time_in_seconds);
    printf("%s: %f ms\n", label, time_in_ms);
}

static inline uint64_t get_cycles()
{
    unsigned int aux;
    return __rdtscp(&aux);
}

void gen_message_1d(poly *Msg) // Accepts poly Msg[M] as a pointer
{
    uint8_t seed[SEEDBYTES];
    int i, k;

    gen_randomness(seed);

    // Initialize GMP random state
    gmp_randstate_t rand_state;
    gmp_randinit_default(rand_state);
    gmp_randseed_ui(rand_state, *(unsigned long *)seed); // For reproducibility

    for (i = 0; i < M; i++)
    {
        for (k = 0; k < N; k++)
        {
            mpz_urandomm(Msg[i].coeffs[k], rand_state, GMP_q);
        }
        poly_reduce(&Msg[i]);
    }

    gmp_randclear(rand_state);
}

void dummy_poly_w(poly *p)
{
    gaussian_sampler_w(p->coeffs, N); // Fill with Gaussian samples
    poly_reduce(p);
}

void benchmark()
{

    int32_t users[USERS];
    for (int i = 0; i < USERS; i++)
    {
        users[i] = i + 1;
    }
    int user = 1;
    // Gen randomness
    uint8_t seed[SEEDBYTES];
    uint8_t nonce = 0;
    gen_randomness(seed);
    uint64_t start, end;
    int Bz = 100;
    // Allocate a dummy poly
    poly dummy;
    poly_init(&dummy);

    // Allocations
    pke_t *pke = malloc(sizeof(pke_t));
    pke->A = malloc(sizeof(poly[KHAT][LHAT]));
    pke->B = malloc(sizeof(poly[KHAT][M]));
    pke->S = malloc(sizeof(poly[LHAT][M]));
    POLY_2D_INIT(pke->A, KHAT, LHAT);
    POLY_2D_INIT(pke->B, KHAT, M);
    POLY_2D_INIT(pke->S, LHAT, M);

    // Setup encryption keys
    keygen_e_1d(pke->A, pke->B, pke->S);

    poly *mu = malloc(sizeof(poly[M]));
    poly_1d_init(mu, M);

    // Sample mu
    gen_message_1d(mu);

    pks_t *pks = malloc(sizeof(pks_t));
    pks->A = malloc(sizeof(poly[K][L]));
    pks->yprime = malloc(sizeof(poly[K]));
    POLY_2D_INIT(pks->A, K, L);
    poly_1d_init(pks->yprime, K);

    // Generate Ae
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < L; j++)
        {
            poly_uniform(&pks->A[i][j], seed, nonce++);
            poly_ntt(&pks->A[i][j]);
            poly_reduce(&pks->A[i][j]);
        }
    }

    poly(*Be)[M] = malloc(sizeof(poly[KHAT][M]));
    POLY_2D_INIT(Be, KHAT, M);

    poly(*Sei)[M] = malloc(sizeof(poly[LHAT][M]));
    POLY_2D_INIT(Sei, LHAT, M);

    // Init new variables
    poly(*si)[L] = malloc(sizeof(poly[USERS][L]));
    POLY_2D_INIT(si, USERS, L);

    poly(*yi)[K] = malloc(sizeof(poly[USERS][K]));
    POLY_2D_INIT(yi, USERS, K);
    uint8_t (*h_yi)[32] = malloc(sizeof(uint8_t[USERS][32]));

    ctx_t *ctx_si = malloc(sizeof(ctx_t[USERS]));
    for (int i = 0; i < USERS; i++)
    {
        ctx_si[i].u = malloc(sizeof(poly[LHAT]));
        ctx_si[i].v = malloc(sizeof(poly[M]));
        poly_1d_init(ctx_si[i].u, LHAT);
        poly_1d_init(ctx_si[i].v, M);
    }

    ctx_t ctx_s;
    ctx_s.u = malloc(sizeof(poly[LHAT]));
    ctx_s.v = malloc(sizeof(poly[M]));
    poly_1d_init(ctx_s.u, LHAT);
    poly_1d_init(ctx_s.v, M);

    // AS Keygen
    printf("as_keygen_1\n");
    start = get_cycles();
    as_keygen_1(pks->A, pke->A, pke->B, si[0], yi[0], h_yi[0], &ctx_si[0]);
    end = get_cycles();

    printf("as_keygen_2\n");
    start = get_cycles();
    as_keygen_2(&ctx_s, yi, h_yi, ctx_si, pks->yprime);
    end = get_cycles();

    poly_1d_clear((poly *)si, USERS * L);
    free(si);
    poly_1d_clear((poly *)yi, USERS * K);
    free(yi);
    free(h_yi);

    // Clear ctx_si
    for (int i = 0; i < USERS; i++)
    {
        poly_1d_clear(ctx_si[i].u, LHAT); // clear polynomials in u
        poly_1d_clear(ctx_si[i].v, M);    // clear polynomials in v
        free(ctx_si[i].u);                // free u
        free(ctx_si[i].v);                // free v
    }
    free(ctx_si); // free the array of ctx_t structs

    // Init ctx_r
    ctx_t *ctx_r = malloc(sizeof(ctx_t[THRESHOLD]));
    for (int i = 0; i < THRESHOLD; i++)
    {
        ctx_r[i].u = malloc(sizeof(poly[LHAT]));
        ctx_r[i].v = malloc(sizeof(poly[M]));
        poly_1d_init(ctx_r[i].u, LHAT);
        poly_1d_init(ctx_r[i].v, M);
    }

    // Populate ctx_r with dummy values
    for (int i = 1; i < THRESHOLD; i++)
    {
        for (int j = 0; j < LHAT; j++)
        {
            poly_copy(&ctx_r[i].u[j], &dummy);
        }
        for (int j = 0; j < M; j++)
        {
            poly_copy(&ctx_r[i].v[j], &dummy);
        }
    }

    start = get_cycles();
    poly *r = malloc(sizeof(poly[M]));
    poly_1d_init(r, M);
    as_sign_1(pks->A, pke->A, pke->B, ctx_r[0].u, ctx_r[0].v, r);
    end = get_cycles();
    print_timing(start, end, "as_sign_round1");

    start = get_cycles();
    as_sign_round2(pke->A, pke->B, ctx_r[0].u, ctx_r[0].v, r);
    end = get_cycles();
    print_timing(start, end, "as_sign_round2");

    // Init sig, ctx_r, ctx_z, dsi, w_i, h_wi
    sig_t sig;
    poly_init(&sig.c); // Initialize challenge polynomial
    sig.z = malloc(sizeof(poly[L]));
    sig.h = malloc(sizeof(poly[K]));
    poly_1d_init(sig.z, L); // Initialize array of L polys
    poly_1d_init(sig.h, K); // Initialize array of K polys

    ctx_t ctx_z;
    ctx_z.u = malloc(sizeof(poly[LHAT]));
    ctx_z.v = malloc(sizeof(poly[M]));
    poly_1d_init(ctx_z.u, LHAT);
    poly_1d_init(ctx_z.v, M);

    poly(*w_i)[K] = malloc(sizeof(poly[THRESHOLD][K]));
    poly(*dsi)[M] = malloc(sizeof(poly[THRESHOLD][M]));
    for (int i = 0; i < THRESHOLD; i++)
    {
        poly_1d_init(w_i[i], K);
        poly_1d_init(dsi[i], M);
    }
    for (int i = 1; i < THRESHOLD; i++)
    {
        for (int j = 0; j < M; j++)
        {
            poly_copy(dsi[i], &dummy);
        }
    }

    uint8_t *h_wi = malloc(sizeof(uint8_t[THRESHOLD]));

    start = get_cycles();
    as_sign_round3(&sig, *pks, ctx_r, ctx_s, ctx_z, Sei, mu, dsi[0], w_i, h_wi, user, users);
    end = get_cycles();
    print_timing(start, end, "as_sign_round3");

    // Init wprime
    poly *wprime = malloc(sizeof(poly[K]));
    poly_1d_init(wprime, K);

    // Clean all variables not used below
    // ctx_r
    for (int i = 0; i < THRESHOLD; i++)
    {
        poly_1d_clear(ctx_r[i].u, LHAT);
        poly_1d_clear(ctx_r[i].v, M);
        free(ctx_r[i].u);
        free(ctx_r[i].v);
    }
    free(ctx_r);

    // ctx_s
    poly_1d_clear(ctx_s.u, LHAT);
    poly_1d_clear(ctx_s.v, M);
    free(ctx_s.u);
    free(ctx_s.v);

    // Sei
    poly_1d_clear((poly *)Sei, LHAT * M);
    free(Sei);

    // w_i
    for (int i = 0; i < THRESHOLD; i++)
    {
        poly_1d_clear(w_i[i], K);
    }
    free(w_i);

    // h_wi
    free(h_wi);

    start = get_cycles();
    as_sign_comb(&ctx_z, pks, &sig, pks->yprime, wprime, dsi, users);
    end = get_cycles();
    print_timing(start, end, "as_sign_comb");

    // Free variables not used by verify
    // ctx_z
    poly_1d_clear(ctx_z.u, LHAT);
    poly_1d_clear(ctx_z.v, M);
    free(ctx_z.u);
    free(ctx_z.v);

    // wprime
    poly_1d_clear(wprime, K);
    free(wprime);

    start = get_cycles();
    as_verify_sig(&sig, pks, mu, Bz);
    end = get_cycles();
    print_timing(start, end, "as_sign_verify");

    // Free all remaining variables
    // sig
    poly_clear(&sig.c);
    poly_1d_clear(sig.z, L);
    poly_1d_clear(sig.h, K);
    free(sig.z);
    free(sig.h);

    // pks
    poly_1d_clear((poly *)pks->A, K * L);
    free(pks->A);
    poly_1d_clear(pks->yprime, K);
    free(pks->yprime);
    free(pks);

    // mu
    poly_1d_clear(mu, M);
    free(mu);
    poly_clear(&dummy);
}

int main(void)
{
    // Init GMP version of #define variables
    init_params();
    init_zetas();
    benchmark();

    // Clear GMP variables
    clear_zetas();
    free_params();
    // free(timings);
    return 0;
}
