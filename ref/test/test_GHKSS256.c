#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "../params.h"
#include "../poly256.h"
#include "../GHKSS256.h"
#include "../gaussian_ref.h"
#include <assert.h>

#define CPU_FREQ 2600000000 // in Hz

typedef struct
{
    int64_t keygen_e_1d;
    int64_t dk_gen_1;
    int64_t dk_gen_2;
    int64_t dk_gen_3;
    int64_t as_keygen_1;
    int64_t as_keygen_2;
    int64_t as_sign_1;
    int64_t as_sign_2;
    int64_t as_sign_3;
    int64_t as_verify;
} timing_results_t;

void print_avg_timing_struct(timing_results_t *timings, int count)
{
    int64_t total_keygen_e_1d = 0;
    int64_t total_dk_gen_1 = 0;
    int64_t total_dk_gen_2 = 0;
    int64_t total_dk_gen_3 = 0;
    int64_t total_as_keygen_1 = 0;
    int64_t total_as_keygen_2 = 0;
    int64_t total_as_sign_1 = 0;
    int64_t total_as_sign_2 = 0;
    int64_t total_as_sign_3 = 0;
    int64_t total_as_verify = 0;

    for (int i = 0; i < count; i++)
    {
        total_keygen_e_1d += timings[i].keygen_e_1d;
        total_dk_gen_1 += timings[i].dk_gen_1;
        total_dk_gen_2 += timings[i].dk_gen_2;
        total_dk_gen_3 += timings[i].dk_gen_3;
        total_as_keygen_1 += timings[i].as_keygen_1;
        total_as_keygen_2 += timings[i].as_keygen_2;
        total_as_sign_1 += timings[i].as_sign_1;
        total_as_sign_2 += timings[i].as_sign_2;
        total_as_sign_3 += timings[i].as_sign_3;
        total_as_verify += timings[i].as_verify;
    }

#define AVG(label, total)                                                      \
    do                                                                         \
    {                                                                          \
        double avg_cycles = (double)(total) / count;                           \
        double avg_seconds = avg_cycles / CPU_FREQ;                            \
        double avg_ms = avg_seconds * 1000.0;                                  \
        printf("%-16s: %10.3f ms (%10.6f sec)\n", label, avg_ms, avg_seconds); \
    } while (0)

    printf("\n--- Average Benchmark Timings (%d runs) ---\n", count);
    AVG("keygen_e_1d", total_keygen_e_1d);
    AVG("dk_gen_1", total_dk_gen_1);
    AVG("dk_gen_2", total_dk_gen_2);
    AVG("dk_gen_3", total_dk_gen_3);
    AVG("as_keygen_1", total_as_keygen_1);
    AVG("as_keygen_2", total_as_keygen_2);
    AVG("as_sign_1", total_as_sign_1);
    AVG("as_sign_2", total_as_sign_2);
    AVG("as_sign_3", total_as_sign_3);
    AVG("as_verify", total_as_verify);

#undef AVG
}

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

    int num_timings = 100;
    timing_results_t timings[num_timings];
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
    // for (int i = 0; i < num_timings; i++)
    // {
    printf("keygen_e_1d (just to generate an A\n");
    start = get_cycles();
    keygen_e_1d(pke->A, pke->B, pke->S);
    end = get_cycles();
    print_timing(start, end, "keygen_e_1d");
    // timings[i].keygen_e_1d = end - start;
    // }
    // print_timing(start, end, "Keygen_e_1d");

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

    poly(*Si)[M] = malloc(sizeof(poly[LHAT][M]));
    POLY_2D_INIT(Si, LHAT, M);

    poly(*Ei)[M] = malloc(sizeof(poly[KHAT][M]));
    POLY_2D_INIT(Ei, KHAT, M);

    poly(*Bi)[M] = malloc(sizeof(poly[KHAT][M]));
    POLY_2D_INIT(Bi, KHAT, M);

    uint8_t (*h_Bi)[32] = malloc(sizeof(uint8_t[USERS][32]));

    start = get_cycles();
    dk_gen_1(pke->A, Si, Ei, Bi, h_Bi[0]);
    end = get_cycles();
    print_timing(start, end, "dk_gen_1");

    start = get_cycles();
    dk_gen_2(
        Si,
        Ei,
        pke,
        users);
    end = get_cycles();
    print_timing(start, end, "dk_gen_2");

    poly_1d_clear((poly *)Si, LHAT * M);
    free(Si);
    poly_1d_clear((poly *)Ei, KHAT * M);
    free(Ei);
    poly_1d_clear((poly *)Bi, KHAT * M);
    free(Bi);

    // Allocate for 1 user, rest we get from dummy polys
    poly(*Bj)[KHAT][M] = malloc(sizeof(poly[USERS][KHAT][M]));
    POLY_3D_INIT(Bj, USERS, KHAT, M);

    // Allocate for 1 user, rest we get from dummy polys
    poly(*Sij)[LHAT][M] = malloc(sizeof(poly[USERS][LHAT][M]));
    POLY_3D_INIT(Sij, USERS, LHAT, M);

    // Init dummy values
    for (int i = 1; i < USERS; i++)
    {
        for (int j = 0; j < KHAT; j++)
        {
            for (int k = 0; k < M; k++)
            {
                poly_copy(&Bj[i][j][k], &dummy);
                poly_copy(&Sij[i][j][k], &dummy);
            }
        }
    }

    poly(*Be)[M] = malloc(sizeof(poly[KHAT][M]));
    POLY_2D_INIT(Be, KHAT, M);

    poly(*Sei)[M] = malloc(sizeof(poly[LHAT][M]));
    POLY_2D_INIT(Sei, LHAT, M);

    printf("dk_gen_3\n");
    start = get_cycles();
    dk_gen_3(
        Bj,
        Sij,
        Be,
        Sei);
    end = get_cycles();
    print_timing(start, end, "dk_gen_3");

    // Clear variables not used anymore
    poly_1d_clear((poly *)Sij, USERS * LHAT * M);
    free(Sij);
    poly_1d_clear((poly *)Bj, USERS * KHAT * M);
    free(Bj);
    // Immediately after dk_gen_1
    free(h_Bi);

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
    for (int i = 0; i < num_timings; i++)
    {
        printf("as_keygen_1\n");
        start = get_cycles();
        as_keygen_1(pks->A, pke->A, pke->B, si[0], yi[0], h_yi[0], &ctx_si[0]);
        end = get_cycles();
        timings[i].as_keygen_1 = end - start;

        printf("as_keygen_2\n");
        start = get_cycles();
        as_keygen_2(&ctx_s, yi, h_yi, ctx_si, pks->yprime);
        end = get_cycles();
        timings[i].as_keygen_2 = end - start;
    }

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

    for (int i = 0; i < num_timings; i++)
    {
        printf("as_sign_1\n");
        start = get_cycles();
        as_sign_1(pks->A, pke->A, pke->B, ctx_r[0]);
        end = get_cycles();
        timings->as_sign_1 = end - start;
    }

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

    uint8_t *h_wi = malloc(sizeof(uint8_t[THRESHOLD]));

    for (int i = 0; i < num_timings; i++)
    {
        printf("as_sign_2\n");
        start = get_cycles();
        as_sign_2(&sig, *pks, ctx_r, ctx_s, ctx_z, Sei, mu, dsi[0], w_i, h_wi, user, users);
        end = get_cycles();
        timings[i].as_sign_2 = end - start;
    }

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

    for (int i = 0; i < num_timings; i++)
    {
        printf("as_sign_3\n");
        start = get_cycles();
        as_sign_3(&ctx_z, pks, &sig, pks->yprime, wprime, dsi, users);
        end = get_cycles();
        timings[i].as_sign_3 = end - start;
    }

    // Free variables not used by verify
    // ctx_z
    poly_1d_clear(ctx_z.u, LHAT);
    poly_1d_clear(ctx_z.v, M);
    free(ctx_z.u);
    free(ctx_z.v);

    // wprime
    poly_1d_clear(wprime, K);
    free(wprime);

    for (int i = 0; i < num_timings; i++)
    {
        printf("as_verify\n");
        start = get_cycles();
        as_verify(&sig, pks, mu, Bz);
        end = get_cycles();
        timings[i].as_verify = end - start;
    }

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
    print_avg_timing_struct(timings, num_timings);
}

int main(void)
{
    // Init GMP version of #define variables
    init_params();
    init_zetas();
    // int num_timings = 1;
    // timing_results_t *timings = malloc(sizeof(timing_results_t) * num_timings);

    // for (int i = 0; i < num_timings; i++)
    // {
    //     // printf("\n\n ----- \nRunning benchmark %d \n ----- \n\n", i);
    //     benchmark(&timings[i]);
    // }
    // print_avg_timing_struct(timings, num_timings);
    benchmark();

    // Clear GMP variables
    clear_zetas();
    free_params();
    // free(timings);
    return 0;
}
