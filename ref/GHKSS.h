#ifndef GHKSS_H
#define GHKSS_H

#include "GHKSS_params.h"
#include "params.h"
#include "poly.h"
#include "reduce.h"

// #define indcpa_mat_keypair DILITHIUM_NAMESPACE(indcpa_mat_keypair)
// void indcpa_mat_keypair(poly A[MAT_m][MAT_m],
//                         poly B[MAT_m][MAT_l],
//                         poly S[MAT_m][MAT_l],
//                         const uint8_t coins[MAT_COINS]);
// void indcpa_mat_keypair();
// static inline void poly_zeroize(poly *p);
static inline void poly_zeroize(poly *p)
{
    memset(p->coeffs, 0, sizeof(p->coeffs));
}

// Function to copy a polynomial
static inline void poly_copy(poly *dest, const poly *src)
{
    for (int i = 0; i < N; i++)
    {
        dest->coeffs[i] = src->coeffs[i];
    }
}

static inline int32_t to_montgomery(int32_t x)
{
    // R^2 mod Q is often called (1 << 32) % Q, or a precomputed constant
    // Let's call it "R2".
    // For Q=8380417, R2=41978  (this is actually mont^2/256 in some references)
    const int32_t R2 = 41978;
    return montgomery_reduce((int64_t)x * R2);
}

static inline void poly_scale_ntt(poly *r, int32_t scalar_mont)
{
    for (int i = 0; i < N; i++)
    {
        int64_t t = (int64_t)r->coeffs[i] * scalar_mont;
        r->coeffs[i] = montgomery_reduce(t);
    }
}

static inline void poly_add_ntt(poly *r, const poly *a, const poly *b)
{
    for (int i = 0; i < N; i++)
    {
        int64_t tmp = (int64_t)a->coeffs[i] + b->coeffs[i];
        tmp %= Q; // or do "tmp -= Q if tmp >= Q"
        r->coeffs[i] = (int32_t)tmp;
    }
}

static inline void poly_sub_ntt(poly *r, const poly *a, const poly *b)
{
    for (int i = 0; i < N; i++)
    {
        int64_t tmp = (int64_t)a->coeffs[i] - b->coeffs[i];
        tmp %= Q; // or do "tmp -= Q if tmp >= Q"
        r->coeffs[i] = (int32_t)tmp;
    }
}

static inline int32_t fqadd(int32_t a, int32_t b)
{
    int64_t r = a + b;
    if (r >= Q)
        r -= Q;
    return (int32_t)r;
}

static inline int32_t fqsub(int32_t a, int32_t b)
{
    int64_t r = a - b;
    if (r < 0)
        r += Q;
    return (int32_t)r;
}

static inline int32_t fqmul(int32_t a, int32_t b)
{
    int64_t r = (int64_t)a * b;
    r %= Q;
    if (r < 0)
        r += Q;
    return (int32_t)r;
}

static int32_t fqpow(int32_t base, int64_t exp)
{
    int32_t result = 1;
    int32_t cur = base;
    while (exp > 0)
    {
        if (exp & 1)
            result = fqmul(result, cur);
        cur = fqmul(cur, cur);
        exp >>= 1;
    }
    return result;
}

// Compute modular inverse with fermats little theorem
static int32_t fqinv(int32_t a)
{
    return fqpow(a, Q - 2);
}

static void lagrange_coeff(int32_t *lambda_i, int32_t user, const int32_t *users, int t)
{
    int32_t j, lambda = 1, numerator, denominator, denominator_inv, fraction;
    for (int x = 0; x < t; x++)
    {
        j = users[x];
        if (j == user)
            continue;
        numerator = fqsub(0, j);
        denominator = fqsub(user, j);
        denominator_inv = fqinv(denominator);
        fraction = fqmul(numerator, denominator_inv);
        lambda = fqmul(lambda, fraction);
    }
    *lambda_i = lambda;
}

static inline void poly_scale(poly *r, uint32_t c)
{
    for (unsigned j = 0; j < N; j++)
    {
        // r->coeffs[j] = (int32_t)(r->coeffs[j] * c);
        int64_t temp = (int64_t)r->coeffs[j] * c;
        // reduce mod Q
        temp %= Q;
        // if (temp < 0)
        //     temp += Q;
        r->coeffs[j] = (int32_t)temp;
    }
    poly_reduce(r);
}

static void poly_mod(poly *r, uint32_t c)
{
    for (unsigned j = 0; j < N; j++)
    {
        r->coeffs[j] %= c; // Reduce modulo smaller modulus (e.g., p)
        if (r->coeffs[j] < 0)
            r->coeffs[j] += c;
    }
}

void indcpa_mat_keypair(poly A[M][M],
                        poly B[M][L],
                        poly S[M][L]);

void indcpa_mat_enc(poly U[M][K],
                    poly V[L][K],
                    poly A[M][M],
                    poly B[M][L],
                    poly Msg[L][K]);

void indcpa_mat_dec(poly Msg[L][K],
                    poly S[M][L],
                    poly U[M][K],
                    poly V[L][K]);

void indcpa_mat_tdec(poly ds_i[L][K],
                     poly S_i[M][L],
                     poly U[M][K],
                     poly V[L][K],
                     const int32_t user,
                     const int32_t *users,
                     const size_t t);

void indcpa_mat_comb(poly Msg[L][K],
                     poly V[L][K],
                     poly (*ds)[L][K],
                     const size_t usize);

void split_key_S(poly (*S_i)[M][L],
                 poly S[M][L],
                 const size_t usize);

void gen_shr_poly(poly P_coeff[][L][K],
                  poly S[L][K],
                  int t);

void compute_shr(poly Share[L][K],
                 poly P_coeff[][L][K],
                 int u, int t);

void reconstruct_secret(poly S[L][K],
                        int32_t *U, int t,
                        poly Shares[][L][K]);

// #define indcpa_mat_enc DILITHIUM_NAMESPACE(indcpa_mat_enc)
// void indcpa_mat_enc(poly U[MAT_m][MAT_k], poly V[MAT_l][MAT_k], poly A[MAT_m][MAT_m], poly B[MAT_m][MAT_l], poly M[MAT_l][MAT_k]);

// #define indcpa_mat_dec DILITHIUM_NAMESPACE(indcpa_mat_dec)
// void indcpa_mat_dec(poly M[MAT_l][MAT_k], poly S[MAT_m][MAT_l], poly U[MAT_m][MAT_k], poly V[MAT_l][MAT_k]);

// // #define poly_uniform DILITHIUM_NAMESPACE(indcpa_mat_keypair)
// // void poly_uniform(poly *a, const uint8_t seed[SEEDBYTES], uint16_t nonce);

// #define lagrange_coeff DILITHIUM_NAMESPACE(lagrange_coeff)
// void lagrange_coeff(int16_t *lambda_i, int16_t user, const int16_t *users, int t);

// #define gen_sharing_poly DILITHIUM_NAMESPACE(gen_sharing_poly)
// void gen_sharing_poly(poly P_coeff[][MAT_l][MAT_k],
//                       poly S[MAT_l][MAT_k],
//                       uint8_t master_seed[SEEDBYTES],
//                       int t);

// #define compute_share DILITHIUM_NAMESPACE(compute_share)
// void compute_share(poly Share[MAT_l][MAT_k],
//                    poly P_coeff[][MAT_l][MAT_k],
//                    int u, int t);

// #define reconstruct_secret DILITHIUM_NAMESPACE(reconstruct_secret)
// void reconstruct_secret(poly S[MAT_l][MAT_k],
//                         int16_t *U, int t,
//                         poly Shares[][MAT_l][MAT_k]);

// #define GHKSS DILITHIUM_NAMESPACE(GHKSS)
void GHKSS(void);

#endif