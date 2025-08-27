#include "GHKSS256.h"
#include "poly256.h"
#include "string.h"
#include <stdio.h>
#include "gaussian_ref.h"

void keygen_e_1d(
    poly (*A)[LHAT], // size: [KHAT][LHAT]
    poly (*B)[M],    // size: [KHAT][M]
    poly (*S)[M])    // size: [LHAT][M]
{
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    gen_randomness(seed);

    // Sample matrix A
    for (int i = 0; i < KHAT; i++)
    {
        for (int j = 0; j < LHAT; j++)
        {
            poly_uniform(&A[i][j], seed, nonce++);
            poly_ntt(&A[i][j]);
            poly_reduce(&A[i][j]);
        }
    }

    // Sample S
    for (int i = 0; i < LHAT; i++)
    {
        for (int j = 0; j < M; j++)
        {
            gaussian_sampler_y(S[i][j].coeffs, N);
            poly_ntt(&S[i][j]);
            poly_reduce(&S[i][j]);
        }
    }
    poly tmp, t;
    poly_init(&tmp);
    poly_init(&t);

    // Compute B = A * S + q * E
    for (int i = 0; i < KHAT; i++)
    {
        for (int k = 0; k < M; k++)
        {
            for (int j = 0; j < LHAT; j++)
            {
                poly_pointwise_montgomery(&tmp, &A[i][j], &S[j][k]);
                poly_add(&B[i][k], &B[i][k], &tmp);
            }

            poly_reduce(&B[i][k]);
            poly_invntt_tomont(&B[i][k]);

            // Add noise qE
            gaussian_sampler(t.coeffs, N);
            poly_scale(&t, GMP_q);
            poly_add(&B[i][k], &B[i][k], &t);
        }
    }
    poly_clear(&tmp);
    poly_clear(&t);
}

void encrypt_1d(poly *U,         // [LHAT]
                poly *V,         // [M]
                poly (*A)[LHAT], // [KHAT][LHAT]
                poly (*B)[M],    // [KHAT][M]
                poly *Msg)       // [M]
{
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    gen_randomness(seed);

    // Allocate and initialize R[KHAT] on heap
    poly *R = malloc(sizeof(poly[KHAT]));
    poly_1d_init(R, KHAT);

    // Sample R
    for (int i = 0; i < KHAT; i++)
    {
        gaussian_sampler(R[i].coeffs, N);
        poly_ntt(&R[i]);
        poly_reduce(&R[i]);
    }

    // Ensure B is in NTT form
    for (int i = 0; i < KHAT; i++)
    {
        for (int j = 0; j < M; j++)
        {
            poly_ntt(&B[i][j]);
            poly_reduce(&B[i][j]);
        }
    }

    // Compute U = A^T * R + qE_2
    for (int i = 0; i < LHAT; i++)
    {
        for (int j = 0; j < KHAT; j++)
        {
            poly t;
            poly_init(&t);
            poly_pointwise_montgomery(&t, &A[j][i], &R[j]);
            poly_add(&U[i], &U[i], &t);
            poly_clear(&t);
        }
        poly_reduce(&U[i]);
        poly_invntt_tomont(&U[i]);

        poly t;
        poly_init(&t);
        gaussian_sampler(t.coeffs, N);
        poly_scale(&t, GMP_q);
        poly_add(&U[i], &U[i], &t);
        poly_reduce(&U[i]);
        poly_clear(&t);
    }

    // Compute V = B^T * R + qE_3 + Msg
    for (int i = 0; i < M; i++)
    {
        poly t;
        poly_init(&t);

        for (int j = 0; j < KHAT; j++)
        {
            poly tmp;
            poly_init(&tmp);
            poly_pointwise_montgomery(&tmp, &B[j][i], &R[j]);
            poly_add(&t, &t, &tmp);
            poly_clear(&tmp);
        }

        poly_reduce(&t);
        poly_invntt_tomont(&t);

        poly noise;
        poly_init(&noise);
        gaussian_sampler(noise.coeffs, N);
        poly_scale(&noise, GMP_q);
        poly_add(&t, &t, &noise);
        poly_clear(&noise);

        poly_add(&V[i], &t, &Msg[i]);
        poly_reduce(&V[i]);
        poly_clear(&t);
    }

    // Free R
    poly_1d_clear(R, KHAT);
    free(R);
}

void gen_sharing_poly(poly (*P_coeff)[LHAT][M], poly (*S)[M])
{
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    gen_randomness(seed);

    for (int i = 0; i < LHAT; i++)
    {
        fast_polyvec_copy(&P_coeff[0][i], &S[i], M); // Copy S to P_coeff[0]
    }

    for (int k = 1; k < THRESHOLD; k++)
    {
        for (int i = 0; i < LHAT; i++)
        {
            for (int j = 0; j < M; j++)
            {
                poly_uniform(&P_coeff[k][i][j], seed, nonce++);
                poly_reduce(&P_coeff[k][i][j]);
            }
        }
    }
}

// First iteration
void compute_share_horner(
    poly (*Share)[M],         // [LHAT][M]
    poly (*P_coeff)[LHAT][M], // [THRESHOLD][LHAT][M]
    int32_t user)             // users[USERS]
{
    for (int i = 0; i < LHAT; i++)
    {
        for (int j = 0; j < M; j++)
        {
            for (int r = THRESHOLD - 1; r >= 0; r--)
            {
                poly_scale_si(&Share[i][j], user);
                poly_add(&Share[i][j], &Share[i][j], &P_coeff[r][i][j]);
                poly_reduce(&Share[i][j]);
            }
        }
    }
}

void compute_share(poly (*Share)[M], poly (*P_coeff)[LHAT][M], int32_t user)
{
    mpz_t umod, tmp, pwr;
    mpz_inits(umod, tmp, pwr, NULL);
    mpz_set_ui(umod, user);
    mpz_set_ui(pwr, 1);

    for (int r = 0; r < THRESHOLD; r++)
    {
        for (int i = 0; i < LHAT; i++)
        {
            for (int j = 0; j < M; j++)
            {
                poly_scale(&P_coeff[r][i][j], pwr);
                poly_add(&Share[i][j], &Share[i][j], &P_coeff[r][i][j]);
                poly_reduce(&Share[i][j]);
            }
        }

        mpz_mul(tmp, umod, pwr);
        mpz_mod(tmp, tmp, GMP_Q);
        mpz_set(pwr, tmp);
    }

    mpz_clears(umod, tmp, pwr, NULL);
}

void reconstruct_secret(poly (*S)[M], int32_t *users, poly (*Shares)[LHAT][M])
{
    poly tpoly;
    poly_init(&tpoly);
    mpz_t lambda;
    mpz_init(lambda);

    for (int k = 0; k < THRESHOLD; k++)
    {
        int u = users[k];
        lagrange_coeff(lambda, u, users, THRESHOLD);

        for (int i = 0; i < LHAT; i++)
        {
            for (int j = 0; j < M; j++)
            {
                poly_copy(&tpoly, &Shares[u - 1][i][j]);
                poly_scale(&tpoly, lambda);
                poly_add(&S[i][j], &S[i][j], &tpoly);
                poly_reduce(&S[i][j]);
            }
        }
    }

    poly_clear(&tpoly);
    mpz_clear(lambda);
}

void tdecrypt_1d(poly *ds_i, poly (*S_i)[M], poly *U, int32_t user, int32_t *users, int t)
{
    uint8_t seed[SEEDBYTES];
    gen_randomness(seed);
    uint16_t nonce = 0;

    mpz_t lambda;
    mpz_init(lambda);
    lagrange_coeff(lambda, user, users, t);

    for (int i = 0; i < M; i++)
    {
        poly t1, t2;
        poly_inits(&t1, &t2, NULL);
        for (int j = 0; j < LHAT; j++)
        {
            poly tmp;
            poly_init(&tmp);
            poly_pointwise_montgomery(&tmp, &S_i[j][i], &U[j]);
            poly_add(&t2, &t2, &tmp);
            poly_clear(&tmp);
        }

        poly_reduce(&t2);
        poly_invntt_tomont(&t2);
        poly_scale(&t2, lambda);

        gaussian_sampler(t1.coeffs, N);
        poly_scale(&t1, GMP_q);
        poly_add(&ds_i[i], &t1, &t2);
        poly_reduce(&ds_i[i]);

        poly_clears(&t1, &t2, NULL);
    }

    mpz_clear(lambda);
}

void combine_1d(
    poly *Msg,     // Output: poly[L]  — accumulator of the result
    poly *V,       // Input: poly[M]   — original ciphertext vector
    poly (*ds)[M], // Input: [usize][M] — decryption shares
    int usize)     // usize = number of users/shares to combine
{
    // Sum decryption shares: Msg = sum_i ds[i]
    for (int i = 0; i < usize; i++)
    {
        for (int j = 0; j < L; j++)
        {
            poly_add(&Msg[j], &Msg[j], &ds[i][j]);
            poly_reduce(&Msg[j]);
        }
    }

    // Msg = V - sum_i ds[i] mod Q, then mod q
    for (int i = 0; i < L; i++)
    {
        poly_sub(&Msg[i], &V[i], &Msg[i]);
        poly_reduce(&Msg[i]);
        poly_mod(&Msg[i], GMP_q);
    }
}

void gen_randomness(uint8_t seed[SEEDBYTES])
{
    uint8_t seedbuf[2 * SEEDBYTES];
    randombytes(seedbuf, SEEDBYTES);
    seedbuf[SEEDBYTES + 0] = M;
    seedbuf[SEEDBYTES + 1] = L;
    shake256(seedbuf, 2 * SEEDBYTES, seedbuf, SEEDBYTES + 2);
    memcpy(seed, seedbuf, SEEDBYTES);
}

void lagrange_coeff(mpz_t lambda, int32_t user, int32_t *users, int t)
{
    mpz_set_ui(lambda, 1);

    // Define and init
    mpz_t numerator, denominator, denominator_inv, fraction, exp;
    mpz_inits(numerator, denominator, denominator_inv, fraction, exp, NULL);

    mpz_sub_ui(exp, GMP_Q, 2);

    for (int x = 0; x < t; x++)
    {
        if (user == users[x])
            continue;

        mpz_set_si(numerator, -users[x]);
        mpz_mod(numerator, numerator, GMP_Q);

        mpz_set_si(denominator, user - users[x]);
        mpz_mod(denominator, denominator, GMP_Q);

        mpz_powm(denominator_inv, denominator, exp, GMP_Q);

        mpz_mul(fraction, numerator, denominator_inv);
        mpz_mod(fraction, fraction, GMP_Q);

        mpz_mul(lambda, lambda, fraction);
        mpz_mod(lambda, lambda, GMP_Q);
    }
    mpz_clears(numerator, denominator, denominator_inv, fraction, exp, NULL);
}

void poly_round(poly *out, const poly *in, mpz_t qw, mpz_t qq)
{
    mpz_t tmp;
    mpz_init(tmp);

    for (int i = 0; i < N; i++)
    {
        mpz_mul(tmp, in->coeffs[i], qq);
        mpz_fdiv_q(tmp, tmp, qw);
        mpz_set(out->coeffs[i], tmp);
    }
    mpz_clear(tmp);
}

void poly_vec_norm2(mpz_t out, poly *v, int len)
{
    mpz_t sum, tmp;
    mpz_init_set_ui(sum, 0);
    mpz_init(tmp);
    for (int i = 0; i < len; i++)
    {
        for (int j = 0; j < N; j++)
        {
            mpz_mul(tmp, v[i].coeffs[j], v[i].coeffs[j]);
            mpz_add(sum, sum, tmp);
        }
    }
    mpz_set(out, sum);
    mpz_clear(sum);
    mpz_clear(tmp);
}

// Safe, streaming version: no giant buf; fixed-width mpz_export with left padding
void compute_challenge(
    poly *c,
    poly (*A)[L],     // [K][L]
    poly yprime[K],   // [K]
    poly wprime[K],   // [K]
    poly mu[L])       // [L]
{
    // Use the SAME modulus you use for coefficients everywhere
    const mpz_t *QQ = &GMP_q;

    const size_t qbits       = mpz_sizeinbase(*QQ, 2);
    const size_t coeff_bytes = (qbits + 7) / 8;   // ceil(bits(q)/8)

    keccak_state st;
    shake256_init(&st);

    // Fixed-size slot for one coefficient
    uint8_t *slot = (uint8_t*)malloc(coeff_bytes);
    if (!slot) { fprintf(stderr, "OOM in compute_challenge\n"); abort(); }

    // Helper lambda-ish (C version): export coeff in [0,q) into fixed-width slot
    #define EXPORT_COEFF_FIXED(gmpz_) do {                                   \
        /* ensure in [0,q) */                                                \
        if (mpz_sgn((gmpz_)) < 0 || mpz_cmp((gmpz_), *QQ) >= 0)               \
            mpz_mod((gmpz_), (gmpz_), *QQ);                                   \
        size_t written = 0;                                                  \
        mpz_export(NULL, &written, 1, 1, 1, 0, (gmpz_));                     \
        if (written > coeff_bytes) {                                         \
            /* should not happen if < q; clamp and recompute */              \
            mpz_mod((gmpz_), (gmpz_), *QQ);                                   \
            written = 0;                                                     \
            mpz_export(NULL, &written, 1, 1, 1, 0, (gmpz_));                 \
        }                                                                     \
        memset(slot, 0, coeff_bytes);                                        \
        mpz_export(slot + (coeff_bytes - written), NULL, 1, 1, 1, 0, (gmpz_)); \
        shake256_absorb(&st, slot, coeff_bytes);                             \
    } while (0)

    // Serialize A (K x L polynomials)
    for (int i = 0; i < K; i++)
        for (int j = 0; j < L; j++)
            for (int k = 0; k < N; k++)
                EXPORT_COEFF_FIXED(A[i][j].coeffs[k]);

    // Serialize yprime (K polys)
    for (int i = 0; i < K; i++)
        for (int k = 0; k < N; k++)
            EXPORT_COEFF_FIXED(yprime[i].coeffs[k]);

    // Serialize wprime (K polys)
    for (int i = 0; i < K; i++)
        for (int k = 0; k < N; k++)
            EXPORT_COEFF_FIXED(wprime[i].coeffs[k]);

    // Serialize mu (L polys)
    for (int i = 0; i < L; i++)
        for (int k = 0; k < N; k++)
            EXPORT_COEFF_FIXED(mu[i].coeffs[k]);

    #undef EXPORT_COEFF_FIXED

    // Finalize and squeeze 32 bytes (adjust if your poly_challenge expects a different size)
    uint8_t hash[32];
    shake256_finalize(&st);
    shake256_squeeze(hash, sizeof hash, &st);

    poly_challenge(c, hash);

    free(slot);
}


void as_keygen_1(
    poly (*As)[L],    // [K][L]
    poly (*Ae)[LHAT], // [KHAT][LHAT]
    poly (*Be)[M],    // [KHAT][M]
    poly *si,         // [L]
    poly *yi,         // [K]
    uint8_t h_yi[32],
    ctx_t *ctx) // ctx->u: [LHAT], ctx->v: [M]
{
    // Sample si
    for (int i = 0; i < L; i++)
    {
        gaussian_sampler_y(si[i].coeffs, N);
        poly_reduce(&si[i]);
        poly_ntt(&si[i]);
    }

    // Sample ei
    poly *ei = malloc(sizeof(poly[K]));
    poly_1d_init(ei, K);
    for (int i = 0; i < K; i++)
    {
        gaussian_sampler_y(ei[i].coeffs, N);
        poly_reduce(&ei[i]);
        poly_ntt(&ei[i]);
    }

    // Compute yi = As * si + ei
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < L; j++)
        {
            poly tmp;
            poly_init(&tmp);
            poly_pointwise_montgomery(&tmp, &As[i][j], &si[j]);
            poly_add(&yi[i], &yi[i], &tmp);
            poly_clear(&tmp);
        }

        poly_invntt_tomont(&yi[i]);
        poly_reduce(&yi[i]);
        poly_add(&yi[i], &yi[i], &ei[i]);
        poly_reduce(&yi[i]);
    }

    // Hash yi with H0
    Hash(yi, h_yi); // yi is an array of poly[K]

    // // Encrypt si under (Ae, Be)
    encrypt_1d(ctx->u, ctx->v, Ae, Be, si);

    // Free temp ei
    poly_1d_clear(ei, K);
    free(ei);
}

void as_keygen_preenc(
    poly (*As)[L],    // [K][L], NTT domain
    poly *si,         // [L], initialized
    poly *yi,         // [K], initialized
    uint8_t h_yi[32]) // output hash
{
    // Sample si in time domain, reduce, then NTT
    for (int i = 0; i < L; i++) {
        gaussian_sampler_y(si[i].coeffs, N);
        poly_reduce(&si[i]);
        poly_ntt(&si[i]);
    }

    // Sample ei in TIME domain (do NOT NTT; we add after inverse NTT)
    poly *ei = (poly *)malloc(sizeof(poly) * K);
    if (!ei) { fprintf(stderr, "OOM in as_keygen_preenc(ei)\n"); abort(); }
    poly_1d_init(ei, K);
    for (int i = 0; i < K; i++) {
        gaussian_sampler_y(ei[i].coeffs, N);
        poly_reduce(&ei[i]);
    }

    // yi = INTT( sum_j As[i][j] ⊙ si[j] ) + ei
    for (int i = 0; i < K; i++) poly_zero(&yi[i]);

    poly tmp;
    poly_init(&tmp);
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < L; j++) {
            poly_pointwise_montgomery(&tmp, &As[i][j], &si[j]);
            poly_add(&yi[i], &yi[i], &tmp);
        }
        poly_invntt_tomont(&yi[i]);
        poly_reduce(&yi[i]);
        poly_add(&yi[i], &yi[i], &ei[i]);  // add noise in time domain
        poly_reduce(&yi[i]);
    }
    poly_clear(&tmp);

    // Hash yi with H0
    Hash(yi, h_yi);

    // Free temp
    poly_1d_clear(ei, K);
    free(ei);
}

// ===== Encryption-only stage =====
void as_keygen_encrypt_only(
    ctx_t *ctx,       // ctx->u: [LHAT], ctx->v: [M], initialized
    poly (*Ae)[LHAT], // [KHAT][LHAT]
    poly (*Be)[M],    // [KHAT][M]
    poly *si)         // [L], NTT domain (as produced by as_keygen_preenc)
{
    encrypt_1d(ctx->u, ctx->v, Ae, Be, si);
}

void as_keygen_2(
    ctx_t *ctx_s,        // ctx_s->u [LHAT], ctx_s->v [M]
    poly (*yi)[K],       // yi[USERS][K]
    uint8_t (*h_yi)[32], // h_yi[USERS][32]
    ctx_t *ctx,          // ctx[USERS]
    poly *yprime)        // yprime[K]
{
    // Allocate temporary accumulator y
    poly *y = malloc(sizeof(poly[K]));
    poly_1d_init(y, K);

    // Compute y = sum_i yi[i]
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < USERS; j++)
        {
            poly_add(&y[i], &y[i], &yi[j][i]);
            poly_reduce(&y[i]);
        }
    }

    // Check H0 hashes
    for (int i = 0; i < USERS; i++)
    {
        uint8_t h_y[32];
        Hash(yi[i], h_y); // yi[i] is poly[K]
        // Hashes will not match with dummy vars for benchmarking
        if (memcmp(h_y, h_yi[i], 32) != 0)
        {
            continue;
            printf("Hashes do not match\n");
            free(y);
            return;
        }
    }

    // Homomorphic aggregation of ctx_s = sum_i ctx[i]
    for (int i = 0; i < USERS; i++)
    {
        for (int j = 0; j < M; j++)
        {
            poly_add(&ctx_s->v[j], &ctx_s->v[j], &ctx[i].v[j]);
            poly_reduce(&ctx_s->v[j]);
        }
        for (int j = 0; j < LHAT; j++)
        {
            poly_add(&ctx_s->u[j], &ctx_s->u[j], &ctx[i].u[j]);
            poly_reduce(&ctx_s->u[j]);
        }
    }

    // Compute yprime = y >> KAPPA_Y (mod q')
    mpz_t qprime;
    mpz_init(qprime);
    mpz_tdiv_q_2exp(qprime, GMP_q, KAPPA_Y);

    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < N; j++)
        {
            mpz_tdiv_q_2exp(yprime[i].coeffs[j], y[i].coeffs[j], KAPPA_Y);
            mpz_mod(yprime[i].coeffs[j], yprime[i].coeffs[j], qprime);
        }
    }

    mpz_clear(qprime);
    poly_1d_clear(y, K);
    free(y);
}
void as_sign_r1(
    poly (*As)[L],     // [K][L], NTT domain
    uint8_t h_w[32],   // output hash of w
    poly **r_out)      // output: allocated r[M] (NTT domain)
{
    // Allocate r ∈ poly[M] and w ∈ poly[K]
    poly *r = (poly*)malloc(sizeof(poly) * M);
    poly *w = (poly*)malloc(sizeof(poly) * K);
    if (!r || !w) { fprintf(stderr, "alloc failed in as_sign_r1\n"); exit(1); }
    poly_1d_init(r, M);
    poly_1d_init(w, K);

    // Sample r for first L slots (NTT), pad [L..M-1] with zeros then NTT
    for (int i = 0; i < L; i++) {
        gaussian_sampler_w(r[i].coeffs, N);
        poly_reduce(&r[i]);
        poly_ntt(&r[i]);
    }
    for (int i = L; i < M; i++) {
        poly_zero(&r[i]);
        poly_ntt(&r[i]); // keep domain consistent
    }

    // Compute w = A * r + e'
    poly tmp; poly_init(&tmp);
    for (int i = 0; i < K; i++) {
        poly_zero(&w[i]); // accumulator in NTT domain
        for (int j = 0; j < L; j++) {
            poly_zero(&tmp);
            poly_pointwise_montgomery(&tmp, &As[i][j], &r[j]);
            poly_add(&w[i], &w[i], &tmp);
        }
        poly_reduce(&w[i]);
        poly_invntt_tomont(&w[i]);     // to time domain

        // Add e' (time domain)
        poly noise; poly_init(&noise);
        gaussian_sampler_w(noise.coeffs, N);
        poly_reduce(&noise);
        poly_add(&w[i], &w[i], &noise);
        poly_reduce(&w[i]);
        poly_clear(&noise);
    }
    poly_clear(&tmp);

    // Hash w -> h_w
    Hash(w, h_w); // (ensure your Hash/H1 uses the fixed-width mpz_export version)

    // Cleanup w, keep r for the caller
    poly_1d_clear(w, K); free(w);
    *r_out = r; // caller (or as_sign_r2) will clear+free r
}

void as_sign_r2(
    poly (*Ae)[LHAT],  // [KHAT][LHAT]
    poly (*Be)[M],     // [KHAT][M]
    poly *u,           // [LHAT], initialized polys
    poly *v,           // [M],    initialized polys
    poly *r,           // [M], NTT domain (from as_sign_r1)
    int  want_free_r)  // if non-zero, this clears+frees r after encrypt
{
    encrypt_1d(u, v, Ae, Be, r);

    if (want_free_r) {
        poly_1d_clear(r, M);
        free(r);
    }
}

void as_sign_1(
    poly (*As)[L],    // [K][L]
    poly (*Ae)[LHAT], // [KHAT][LHAT]
    poly (*Be)[M],    // [KHAT][M]
    poly *u,          // [LHAT]
    poly *v,           // [M]
    poly *r // M
)
{
    // Allocate r ∈ poly[M] and w ∈ poly[K]
    poly *w = malloc(sizeof(poly[K]));
    poly_1d_init(w, K);

    // Sample r for first L slots, pad with 0s for M-L
    for (int i = 0; i < L; i++)
    {
        gaussian_sampler_w(r[i].coeffs, N);
        poly_reduce(&r[i]);
        poly_ntt(&r[i]);
    }
    for (int i = L; i < M; i++)
    {
        poly_ntt(&r[i]); // Still transform to NTT domain
    }

    // Compute w = A * r + e'
    poly tmp;
    poly_init(&tmp);
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < L; j++)
        {
            poly_zero(&tmp);
            poly_pointwise_montgomery(&tmp, &As[i][j], &r[j]);
            poly_add(&w[i], &w[i], &tmp);
            // poly_reduce(&w[i]);
        }

        poly_reduce(&w[i]);
        poly_invntt_tomont(&w[i]);

        // Sample and add eprime
        poly noise;
        poly_init(&noise);
        gaussian_sampler_w(noise.coeffs, N);
        poly_reduce(&noise);
        poly_add(&w[i], &w[i], &noise);
        poly_reduce(&w[i]);
        poly_clear(&noise);
    }
    poly_clear(&tmp);

    // Compute h_w = H1(w)
    uint8_t h_w[32];
    Hash(w, h_w); // w is an array of poly[K]

    // Encrypt r under Ae, Be
    // encrypt_1d(u, v, Ae, Be, r);

    // Cleanup
    poly_1d_clear(w, K);
    free(w);
}

void as_sign_round2(
    poly (*Ae)[LHAT], // [KHAT][LHAT]
    poly (*Be)[M],    // [KHAT][M]
    poly *u,          // [LHAT]
    poly *v,           // [M]
    poly *r // M
)
{

    // Encrypt r under Ae, Be
    encrypt_1d(u, v, Ae, Be, r);
    // Cleanup
    poly_1d_clear(r, M);
    free(r);
}

void as_sign_1_old(
    poly (*As)[L],    // [K][L]
    poly (*Ae)[LHAT], // [KHAT][LHAT]
    poly (*Be)[M],    // [KHAT][M]
    ctx_t ctx_r)      // ctx_si.u [LHAT], ctx_si.v [M]
{
    // Allocate r ∈ poly[M] and w ∈ poly[K]
    poly *r = malloc(sizeof(poly[M]));
    poly *w = malloc(sizeof(poly[K]));
    poly_1d_init(r, M);
    poly_1d_init(w, K);

    // Sample r for first L slots, pad with 0s for M-L
    for (int i = 0; i < L; i++)
    {
        gaussian_sampler_w(r[i].coeffs, N);
        poly_reduce(&r[i]);
        poly_ntt(&r[i]);
    }
    for (int i = L; i < M; i++)
    {
        poly_ntt(&r[i]); // Still transform to NTT domain
    }

    // Compute w = A * r + e'
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < L; j++)
        {
            poly tmp;
            poly_init(&tmp);
            poly_pointwise_montgomery(&tmp, &As[i][j], &r[j]);
            poly_add(&w[i], &w[i], &tmp);
            poly_clear(&tmp);
        }

        poly_reduce(&w[i]);
        poly_invntt_tomont(&w[i]);

        // Sample and add eprime
        poly noise;
        poly_init(&noise);
        gaussian_sampler_w(noise.coeffs, N);
        poly_reduce(&noise);
        poly_add(&w[i], &w[i], &noise);
        poly_reduce(&w[i]);
        poly_clear(&noise);
    }

    // Compute h_w = H1(w)
    uint8_t h_w[32];
    H1(w, h_w); // w is an array of poly[K]

    // Encrypt r under Ae, Be
    encrypt_1d(ctx_r.u, ctx_r.v, Ae, Be, r);

    // Cleanup
    poly_1d_clear(r, M);
    poly_1d_clear(w, K);
    free(r);
    free(w);
}

void as_sign_2(
    poly *c,           // just one poly
    poly *z,           // L
    poly *h,           // K
    poly (*A)[L],      // KxL
    poly *yprime,      // K
    poly (*u_r)[LHAT], // THRESHOLD x LHAT
    poly (*v_r)[M],    // THRESHOLD x M
    poly *u_s,         // LHAT
    poly *v_s,         // M
    poly *u_z,         // LHAT
    poly *v_z,         // M
    poly (*ski)[M],    // LHAT x M
    poly *mu,          // L, but padded to M
    poly *dsi,         // L, but padded to M
    poly (*w_i)[K],    // THRESHOLD x K
    uint8_t *h_wi[THRESHOLD],
    int user,
    int32_t *users)
{
    // Compute w = sum(w_i) for all i \in U
    // Compute rounded w, round by q_w
    poly w[K];
    poly_1d_init(w, K);
    for (int i = 0; i < THRESHOLD; i++)
    {
        for (int j = 0; j < K; j++)
        {
            poly_add(&w[j], &w[j], &w_i[i][j]);
            poly_reduce(&w[j]);
        }
    }
    // Round w
    mpz_t qprime;
    mpz_init(qprime);
    mpz_tdiv_q_2exp(qprime, GMP_q, KAPPA_W);
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < N; j++)
        {
            mpz_tdiv_q_2exp(w[i].coeffs[j], w[i].coeffs[j], KAPPA_W);
            mpz_mod(w[i].coeffs[j], w[i].coeffs[j], qprime);
        }
    }
    mpz_clear(qprime);

    // Compute c = H2(w, pk_S, mu)
    compute_challenge(c, A, yprime, w, mu);

    // Compute ctx_z = c * ctx_s + sum(ctx_rj) for all j \in U
    for (int i = 0; i < THRESHOLD; i++)
    {
        for (int j = 0; j < LHAT; j++)
        {
            poly_add(&u_z[j], &u_z[j], &u_r[i][j]);
            poly_reduce(&u_z[j]);
        }
        for (int j = 0; j < M; j++)
        {
            poly_add(&v_z[j], &v_z[j], &v_r[i][j]);
            poly_reduce(&v_z[j]);
        }
    }

    poly tmp;
    poly_init(&tmp);
    // Compute c * ctx_s
    // for (int i = 0; i < M; i++) // M = LHAT here
    // {
    //     poly_zero(&tmp);
    //     poly_pointwise_montgomery(&tmp, c, &v_s[i]);
    //     poly_add(&v_z[i], &v_z[i], &tmp);
    //     poly_reduce(&v_z[i]);
    // }
    for (int i = 0; i < LHAT; i++)
    {
        poly_zero(&tmp);
        poly_pointwise_montgomery(&tmp, c, &v_s[i]);
        poly_add(&v_z[i], &v_z[i], &tmp);
        poly_reduce(&v_z[i]);
        poly_zero(&tmp);
        poly_pointwise_montgomery(&tmp, c, &u_s[i]);
        poly_add(&u_z[i], &u_z[i], &tmp);
        poly_reduce(&u_z[i]);
    }

    poly_clear(&tmp);
    // Compute ds_i = tdec(ctx_z, sk_i, U)
    tdecrypt_1d(dsi, ski, u_z, user, users, THRESHOLD);
}

void as_sign_round3(
    sig_t *sig,
    pks_t pks,
    ctx_t ctx_r[THRESHOLD],
    ctx_t ctx_s,
    ctx_t ctx_z,
    poly ski[LHAT][M],
    poly mu[M],  // L, but padded to M
    poly dsi[M], // L, but padded to M
    poly w_i[THRESHOLD][K],
    uint8_t h_wi[THRESHOLD],
    int user,
    int32_t users[USERS])
{
    // Compute w = sum(w_i) for all i \in U
    // Compute rounded w, round by q_w
    poly w[K];
    poly_1d_init(w, K);
    for (int i = 0; i < THRESHOLD; i++)
    {
        for (int j = 0; j < K; j++)
        {
            poly_add(&w[j], &w[j], &w_i[i][j]);
            poly_reduce(&w[j]);
        }
    }
    // Round w
    mpz_t qprime;
    mpz_init(qprime);
    mpz_tdiv_q_2exp(qprime, GMP_q, KAPPA_W);
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < N; j++)
        {
            mpz_tdiv_q_2exp(w[i].coeffs[j], w[i].coeffs[j], KAPPA_W);
            mpz_mod(w[i].coeffs[j], w[i].coeffs[j], qprime);
        }
    }
    mpz_clear(qprime);

    // Compute c = H2(w, pk_S, mu)
    compute_challenge(&sig->c, pks.A, pks.yprime, w, mu);


    for (int i = 0; i < THRESHOLD; i++)
    {
        for (int j = 0; j < LHAT; j++)
        {
            poly_add(&ctx_z.u[j], &ctx_z.u[j], &ctx_r[i].u[j]);
            poly_reduce(&ctx_z.u[j]);
        }
        for (int j = 0; j < M; j++)
        {
            poly_add(&ctx_z.v[j], &ctx_z.v[j], &ctx_r[i].v[j]);
            poly_reduce(&ctx_z.v[j]);
        }
    }

    // Compute c * ctx_s
    for (int i = 0; i < M; i++)
    {
        poly tmp;
        poly_init(&tmp);
        poly_pointwise_montgomery(&tmp, &sig->c, &ctx_s.v[i]);
        poly_add(&ctx_z.v[i], &ctx_z.v[i], &tmp);
        poly_reduce(&ctx_z.v[i]);
        poly_clear(&tmp);
    }
    for (int i = 0; i < LHAT; i++)
    {
        poly tmp;
        poly_init(&tmp);
        poly_pointwise_montgomery(&tmp, &sig->c, &ctx_s.u[i]);
        poly_add(&ctx_z.u[i], &ctx_z.u[i], &tmp);
        poly_reduce(&ctx_z.u[i]);
        poly_clear(&tmp);
    }

    // Compute ds_i = tdec(ctx_z, sk_i, U)
    tdecrypt_1d(dsi, ski, ctx_z.u, user, users, THRESHOLD);
}

void as_sign_comb(
    ctx_t *ctx_z,         // ctx_z->v [M]
    const pks_t *pks,     // pks->A [K][L]
    sig_t *sig,           // sig->c, sig->z [L], sig->h [K]
    const poly *yprime,   // [K]
    poly *wprime,         // [K]
    poly (*dsi)[M],       // [THRESHOLD][M]
    const int32_t *users) // [USERS]
{
    // Compute z = combine(ctx_z->v, dsi)
    combine_1d(sig->z, ctx_z->v, dsi, THRESHOLD);

    // Allocate temporary t vector
    poly *t = malloc(sizeof(poly[K]));
    poly_1d_init(t, K);

    // Compute t = A*z - 2^kappa_y * c * y'
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < L; j++)
        {
            poly tmp;
            poly_init(&tmp);
            poly_pointwise_montgomery(&tmp, &pks->A[i][j], &sig->z[j]);
            poly_add(&t[i], &t[i], &tmp);
            poly_reduce(&t[i]);
            poly_clear(&tmp);
        }
    }

    for (int i = 0; i < K; i++)
    {
        poly tmp;
        poly_init(&tmp);
        poly_pointwise_montgomery(&tmp, &sig->c, &yprime[i]);
        poly_scale(&tmp, GMP_KAPPA_Y);
        poly_reduce(&tmp);
        poly_sub(&t[i], &t[i], &tmp);
        poly_reduce(&t[i]);
        poly_clear(&tmp);
    }

    // Round t[i] = t[i] >> KAPPA_W mod q'
    mpz_t qprime;
    mpz_init(qprime);
    mpz_tdiv_q_2exp(qprime, GMP_q, KAPPA_W);

    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < N; j++)
        {
            mpz_tdiv_q_2exp(t[i].coeffs[j], t[i].coeffs[j], KAPPA_W);
            mpz_mod(t[i].coeffs[j], t[i].coeffs[j], qprime);
        }
    }

    mpz_clear(qprime);

    // Compute h = w' - t
    for (int i = 0; i < K; i++)
    {
        poly_sub(&sig->h[i], &wprime[i], &t[i]);
        poly_reduce(&sig->h[i]);
    }

    poly_1d_clear(t, K);
    free(t);
}

int as_verify_sig(
    sig_t *sig,
    pks_t *pks,
    poly *mu,
    int Bz)
{

    // Compute w* = [As * z -2^{kappa_y} * c * yprime]_qw
    poly wstar[K];
    poly_1d_init(wstar, K);
    // As * z
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < L; j++)
        {
            poly tmp;
            poly_init(&tmp);

            poly_pointwise_montgomery(&tmp, &pks->A[i][j], &sig->z[j]);
            poly_add(&wstar[i], &wstar[i], &tmp);
            poly_reduce(&wstar[i]);

            poly_clear(&tmp);
        }
    }
    // 2^{kappa_y} * c * yprime
    for (int i = 0; i < K; i++)
    {
        poly tmp;
        poly_init(&tmp);

        poly_pointwise_montgomery(&tmp, &sig->c, &pks->yprime[i]);
        poly_scale(&tmp, GMP_KAPPA_Y);
        poly_reduce(&tmp);
        poly_sub(&wstar[i], &wstar[i], &tmp);
        poly_reduce(&wstar[i]);

        poly_clear(&tmp);
    }
    // Round t
    mpz_t qprime;
    mpz_init(qprime);
    mpz_tdiv_q_2exp(qprime, GMP_q, KAPPA_W);
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < N; j++)
        {
            mpz_tdiv_q_2exp(wstar[i].coeffs[j], wstar[i].coeffs[j], KAPPA_W);
            mpz_mod(wstar[i].coeffs[j], wstar[i].coeffs[j], qprime);
        }
    }
    mpz_clear(qprime);
    for (int i = 0; i < K; i++)
    {
        poly_add(&wstar[i], &wstar[i], &sig->h[i]);
        poly_reduce(&wstar[i]);
    }

    poly cstar;
    poly_init(&cstar);

    compute_challenge(&cstar, pks->A, pks->yprime, wstar, mu);

    // Clear
    poly_1d_clear(wstar, K);

    return polycmp(&cstar, &sig->c);
}

void dk_gen_1(
    poly (*Ae)[LHAT], // Ae: KHAT x LHAT
    poly (*Si)[M],    // Si: LHAT x M
    poly (*Ei)[M],    // Ei: KHAT x M
    poly (*Bi)[M],    // Bi: KHAT x M
    uint8_t *h_Bi)    // 32-byte output hash
{

    // Sample Si
    for (int i = 0; i < LHAT; i++)
    {
        for (int j = 0; j < M; j++)
        {
            binary_sampler(Si[i][j].coeffs, N);
            poly_ntt(&Si[i][j]);
            poly_reduce(&Si[i][j]);
        }
    }
    poly tmp;
    poly t;
    poly_init(&tmp);
    poly_init(&t);
    // Compute Bi = pke.A * Si + q* Ei
    for (int i = 0; i < KHAT; i++)
    {
        for (int k = 0; k < M; k++)
        {
            for (int j = 0; j < LHAT; j++)
            {
                poly_zero(&tmp);
                poly_pointwise_montgomery(&tmp, &Ae[i][j], &Si[j][k]);
                poly_add(&Bi[i][k], &Bi[i][k], &tmp);
                poly_reduce(&Bi[i][k]);
            }
            poly_invntt_tomont(&Bi[i][k]);

            // Add noise q*Ei
            binary_sampler(Ei[i][k].coeffs, N);
            poly_reduce(&Ei[i][k]);
            // poly_copy(&Ei[i][k], &t);
            poly_scale(&Ei[i][k], GMP_q);
            poly_add(&Bi[i][k], &Bi[i][k], &Ei[i][k]);
            poly_reduce(&Bi[i][k]);
        }
    }
    poly_clear(&tmp);
    poly_clear(&t);

    // Hash Bi
    H0_matrix(Bi, KHAT, M, h_Bi);
}


void poly_matmul_acc(
    poly (*out)[M],  // [KHAT][M]
    poly (*A)[LHAT], // [KHAT][LHAT]
    poly (*S)[M],    // [LHAT][M]
    poly *tmp        // temporary scratch polynomial
)
{
    for (int i = 0; i < KHAT; i++)
    {
        for (int k = 0; k < M; k++)
        {
            // First term
            poly_pointwise_montgomery(&out[i][k], &A[i][0], &S[0][k]);
            for (int j = 1; j < LHAT; j++)
            {
                poly_pointwise_montgomery(tmp, &A[i][j], &S[j][k]);
                poly_add(&out[i][k], &out[i][k], tmp);
            }
            poly_invntt_tomont(&out[i][k]);
            poly_reduce(&out[i][k]);
        }
    }
}

void dk_gen_2(
    poly (*Si)[M],    // [LHAT][M]
    poly (*Ei)[M],    // [LHAT][M]
    poly (*Ae)[LHAT], // [KHAT][LHAT]
    int32_t *users)   // users[USERS]
{
    // Allocate sharing polynomial coefficients
    poly(*P_coeff_S)[LHAT][M] = malloc(sizeof(poly) * THRESHOLD * LHAT * M);
    poly(*P_coeff_E)[KHAT][M] = malloc(sizeof(poly) * THRESHOLD * KHAT * M);

    // Since benchmarking not needed for other than computation, need to store and send in real application, will be big
    poly(*Shares_S)[M] = malloc(sizeof(poly) * LHAT * M); // [LHAT][M]
    poly(*Shares_E)[M] = malloc(sizeof(poly) * KHAT * M); // [KHAT][M]
    poly(*Bij)[M] = malloc(sizeof(poly) * KHAT * M);      // [USERS][KHAT][M]

    POLY_2D_INIT(Shares_S, LHAT, M);
    POLY_2D_INIT(Shares_E, KHAT, M);
    POLY_2D_INIT(Bij, KHAT, M);

    // Initialize
    POLY_3D_INIT(P_coeff_S, THRESHOLD, LHAT, M);
    POLY_3D_INIT(P_coeff_E, THRESHOLD, KHAT, M);

    // Generate sharing polynomials from Si and Ei
    gen_sharing_poly(P_coeff_S, Si);
    gen_sharing_poly(P_coeff_E, Ei);

    // Generate secret shares for each user
    poly tmp;
    poly_init(&tmp);
    for (int u = 0; u < USERS; u++)
    {
        compute_share_horner(Shares_S, P_coeff_S, users[u]);
        compute_share_horner(Shares_E, P_coeff_E, users[u]);

        // Compute Bij = A * Shares_S[u] + q * Shares_E[u]
        poly_matmul_acc(Bij, Ae, Shares_S, &tmp);

        for (int i = 0; i < KHAT; i++)
        {
            for (int k = 0; k < M; k++)
            {
                poly_scale(&Shares_E[i][k], GMP_q);
                poly_add(&Bij[i][k], &Bij[i][k], &Shares_E[i][k]);
                poly_reduce(&Bij[i][k]);
            }
        }
    }
    poly_clear(&tmp);

    // Clear temporary sharing polynomials
    poly_1d_clear((poly *)P_coeff_S, THRESHOLD * LHAT * M);
    poly_1d_clear((poly *)P_coeff_E, THRESHOLD * KHAT * M);
    poly_1d_clear(Shares_S, LHAT * M);
    poly_1d_clear(Shares_E, KHAT * M);
    poly_1d_clear(Bij, KHAT * M);
    free(P_coeff_S);
    free(P_coeff_E);
    free(Shares_S);
    free(Shares_E);
    free(Bij);
}

void dk_gen_3(
    poly (*Bj)[KHAT][M],  // [USERS][KHAT][M]
    poly (*Sij)[LHAT][M], // [USERS][LHAT][M]
    poly (*Be)[M],        // [KHAT][M]
    poly (*Sei)[M])       // [LHAT][M]
{
    // Compute Be = sum_i Bj[i]
    for (int i = 0; i < USERS; i++)
    {
        for (int j = 0; j < KHAT; j++)
        {
            for (int k = 0; k < M; k++)
            {
                poly_add(&Be[j][k], &Be[j][k], &Bj[0][j][k]);
                poly_reduce(&Be[j][k]);
            }
        }
    }

    // Compute Sei = sum_i Sij[i]
    for (int i = 0; i < USERS; i++)
    {
        for (int j = 0; j < LHAT; j++)
        {
            for (int k = 0; k < M; k++)
            {
                poly_add(&Sei[j][k], &Sei[j][k], &Sij[0][j][k]);
                poly_reduce(&Sei[j][k]);
            }
        }
    }
}
