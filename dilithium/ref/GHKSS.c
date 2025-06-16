#include <stdio.h>
#include "params.h"
#include "GHKSS_params.h"
#include "GHKSS.h"
#include "randombytes.h"
#include "polyvec.h"
#include "poly.h"
#include "symmetric.h"
#include "fips202.h"
#include "reduce.h"

// static inline void poly_mod(poly *r, uint32_t c)
// {
//     for (unsigned j = 0; j < N; j++)
//     {
//         r->coeffs[j] = (int32_t)(r->coeffs[j] % c);
//     }
// }

void indcpa_mat_keypair(poly A[M][M],
                        poly B[M][L],
                        poly S[M][L])
{
    uint8_t seedbuf[2 * SEEDBYTES + CRHBYTES];
    uint8_t tr[TRBYTES];
    const uint8_t *rho, *rhoprime;
    uint16_t nonce = 0;
    /* Get randomness for rho, rhoprime and key */
    randombytes(seedbuf, SEEDBYTES);
    seedbuf[SEEDBYTES + 0] = M;
    seedbuf[SEEDBYTES + 1] = L;
    shake256(seedbuf, 2 * SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES + 2);
    rho = seedbuf;
    rhoprime = rho + SEEDBYTES;

    // For temporary polynomials
    poly t;
    poly_zeroize(&t);

    /* Expand matrix */
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            poly_uniform(&A[i][j], rho, nonce++);
            poly_ntt(&A[i][j]);
            poly_reduce(&A[i][j]);
        }
    }

    // Sample S
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < L; j++)
        {
            poly_uniform_eta(&S[i][j], rhoprime, nonce++);
            poly_ntt(&S[i][j]);
            poly_reduce(&S[i][j]);
        }
    }

    for (int i = 0; i < M; i++)
    {
        for (int k = 0; k < L; k++)
        {
            poly_zeroize(&B[i][k]);
            for (int j = 0; j < M; j++)
            {
                poly_pointwise_montgomery(&t, &A[i][j], &S[j][k]);
                poly_add(&B[i][k], &B[i][k], &t);
            }
            poly_reduce(&B[i][k]);
            poly_invntt_tomont(&B[i][k]);

            // Add noise qE_1
            poly_uniform_eta(&t, rhoprime, nonce++);
            poly_scale(&t, MAT_q);
            poly_add(&B[i][k], &B[i][k], &t);
            // poly_reduce(&B[i][k]);
        }
    }
}

void indcpa_mat_enc(poly U[M][K],
                    poly V[L][K],
                    poly A[M][M],
                    poly B[M][L],
                    poly Msg[L][K])
{
    uint8_t seedbuf[2 * SEEDBYTES + CRHBYTES];
    uint8_t tr[TRBYTES];
    const uint8_t *rho, *rhoprime;
    /* Get randomness for rho, rhoprime and key */
    randombytes(seedbuf, SEEDBYTES);
    seedbuf[SEEDBYTES + 0] = M;
    seedbuf[SEEDBYTES + 1] = L;
    shake256(seedbuf, 2 * SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES + 2);
    rho = seedbuf;
    rhoprime = rho + SEEDBYTES;
    uint16_t nonce = 0;

    int i, j, k;
    poly R[M][K];
    // For temporary polynomials
    poly t;

    // sample R mxk
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < K; j++)
        {
            poly_uniform_eta(&R[i][j], rhoprime, nonce++);
            poly_ntt(&R[i][j]);
        }
    }

    for (i = 0; i < M; i++)
    {
        for (j = 0; j < L; j++)
        {
            poly_ntt(&B[i][j]);
            poly_reduce(&B[i][j]);
        }
    }

    // Compute U = A^\trans R + qE_2
    for (i = 0; i < M; i++)
    {
        for (k = 0; k < K; k++)
        {
            poly_zeroize(&U[i][k]);
            poly_zeroize(&t);
            for (j = 0; j < M; j++)
            {
                poly_pointwise_montgomery(&t, &A[j][i], &R[j][k]);
                poly_add(&U[i][k], &U[i][k], &t);
            }
            poly_reduce(&U[i][k]);
            poly_invntt_tomont(&U[i][k]);

            poly_uniform_eta(&t, rhoprime, nonce++);
            poly_scale(&t, MAT_q);
            poly_add(&U[i][k], &U[i][k], &t);
            // poly_reduce(&U[i][k]);
        }
    }

    // Compute V = B^\trans R +qE_3 + M
    for (i = 0; i < L; i++)
    {
        for (k = 0; k < K; k++)
        {
            poly_zeroize(&V[i][k]);
            poly_zeroize(&t);
            for (j = 0; j < M; j++)
            {
                poly_pointwise_montgomery(&t, &B[j][i], &R[j][k]);
                poly_add(&V[i][k], &V[i][k], &t);
            }
            poly_reduce(&V[i][k]);
            poly_invntt_tomont(&V[i][k]);

            poly_uniform_eta(&t, rhoprime, nonce++);
            poly_scale(&t, MAT_q);
            poly_add(&V[i][k], &V[i][k], &t);
            poly_reduce(&V[i][k]);

            // Add M
            poly_add(&V[i][k], &V[i][k], &Msg[i][k]);
            // poly_reduce(&V[i][k]);
        }
    }
}

void indcpa_mat_dec(poly Msg[L][K],
                    poly S[M][L],
                    poly U[M][K],
                    poly V[L][K])
{
    unsigned int i, j, k;
    poly t, tmp;

    for (i = 0; i < M; i++)
    {
        for (j = 0; j < K; j++)
        {
            poly_ntt(&U[i][j]);
            poly_reduce(&U[i][j]);
        }
    }

    // Compute V - S^\trans U
    for (i = 0; i < L; i++)
    {
        for (k = 0; k < K; k++)
        {
            poly_zeroize(&t);
            poly_zeroize(&tmp);
            poly_zeroize(&Msg[i][k]);
            for (j = 0; j < M; j++)
            {
                poly_pointwise_montgomery(&tmp, &S[j][i], &U[j][k]);
                poly_add(&t, &t, &tmp);
            }
            poly_reduce(&t);
            poly_invntt_tomont(&t);

            // Compute M[i][k] = V[i][k] - t
            poly_sub(&Msg[i][k], &V[i][k], &t);
            poly_reduce(&Msg[i][k]);
            poly_mod(&Msg[i][k], MAT_q);
        }
    }
}

void indcpa_mat_tdec(poly ds_i[L][K],
                     poly S_i[M][L],
                     poly U[M][K],
                     poly V[L][K],
                     const int32_t user,
                     const int32_t *users,
                     const size_t t)
{
    /* Get randomness for rho, rhoprime and key */
    uint8_t seedbuf[2 * SEEDBYTES + CRHBYTES];
    uint8_t tr[TRBYTES];
    const uint8_t *rho, *rhoprime;
    randombytes(seedbuf, SEEDBYTES);
    seedbuf[SEEDBYTES + 0] = M;
    seedbuf[SEEDBYTES + 1] = L;
    shake256(seedbuf, 2 * SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES + 2);
    rho = seedbuf;
    rhoprime = rho + SEEDBYTES;
    uint16_t nonce = 0;

    // Variables
    unsigned int i, j, k;
    poly t1, t2;

    // Compute Lagrange coeff
    int32_t lambda;
    lagrange_coeff(&lambda, user, users, t);

    // Compute ds_i = \lambda_i S_i^\trans U + q E_i
    for (i = 0; i < L; i++)
    {
        for (k = 0; k < K; k++)
        {
            poly_zeroize(&t1);
            poly_zeroize(&t2);
            poly_zeroize(&ds_i[i][k]);

            for (j = 0; j < M; j++)
            {
                poly_pointwise_montgomery(&t1, &S_i[j][i], &U[j][k]);
                poly_add(&t2, &t2, &t1);
            }
            poly_reduce(&t2);
            poly_invntt_tomont(&t2);
            poly_scale(&t2, lambda); // scale includes reduce

            // Sample q \cdot E_ij
            poly_zeroize(&t1);
            poly_uniform_eta(&t1, rhoprime, nonce++);
            poly_scale(&t1, MAT_q);
            poly_add(&ds_i[i][k], &t1, &t2);
        }
    }
}

void indcpa_mat_comb(poly Msg[L][K],
                     poly V[L][K],
                     poly (*ds)[L][K],
                     const size_t usize)
{
    size_t i, j, k;

    for (i = 0; i < L; i++)
    {
        for (j = 0; j < K; j++)
        {
            poly_zeroize(&Msg[i][j]);
        }
    }
    // Sum decryption shares
    for (i = 0; i < usize; i++)
    {
        for (j = 0; j < L; j++)
        {
            for (k = 0; k < K; k++)
            {
                poly_add(&Msg[j][k], &Msg[j][k], &ds[i][j][k]);
                poly_reduce(&Msg[j][k]);
            }
        }
    }

    // (V - \sum{ds}_{i\in\users}) \mod Q \mod q
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < K; j++)
        {
            poly_sub(&Msg[i][j], &V[i][j], &Msg[i][j]);
            poly_reduce(&Msg[i][j]);
            poly_mod(&Msg[i][j], MAT_q);
        }
    }
}

// Really wide gaussian for noise
// Commitments
// Homomorphic properties
// Signature

void split_key_S(poly (*S_i)[M][L],
                 poly S[M][L],
                 const size_t usize)
{
    uint8_t seedbuf[2 * SEEDBYTES + CRHBYTES];
    uint8_t tr[TRBYTES];
    const uint8_t *rho, *rhoprime;
    randombytes(seedbuf, SEEDBYTES);
    seedbuf[SEEDBYTES + 0] = M;
    seedbuf[SEEDBYTES + 1] = L;
    shake256(seedbuf, 2 * SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES + 2);
    rho = seedbuf;
    rhoprime = rho + SEEDBYTES;
    uint16_t nonce = 0;

    for (size_t u = 0; u < usize - 1; u++)
    {
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < L; j++)
            {
                poly_uniform_eta(&S_i[u][i][j], rhoprime, nonce++);
                poly_ntt(&S_i[u][i][j]);
                poly_reduce(&S_i[u][i][j]);
            }
        }
    }

    // Compute last key share
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < L; j++)
        {
            // Init last keyshare to S
            // poly_zeroize(&S_i[usize - 1][i][j]);
            poly_copy(&S_i[usize - 1][i][j], &S[i][j]);
            // Subtract all previous S_i from last keyshare
            for (size_t u = 0; u < usize - 1; u++)
            {
                poly_sub(&S_i[usize - 1][i][j], &S_i[usize - 1][i][j], &S_i[u][i][j]);
            }
        }
    }
}

void gen_shr_poly(poly P_coeff[][M][L],
                  poly S[M][L],
                  int t)
{
    // SEED
    uint8_t seedbuf[2 * SEEDBYTES + CRHBYTES];
    uint8_t tr[TRBYTES];
    const uint8_t *rho, *rhoprime;
    randombytes(seedbuf, SEEDBYTES);
    seedbuf[SEEDBYTES + 0] = M;
    seedbuf[SEEDBYTES + 1] = L;
    shake256(seedbuf, 2 * SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES + 2);
    rho = seedbuf;
    rhoprime = rho + SEEDBYTES;
    uint16_t nonce = 0;

    int i, j, k;
    for (k = 0; k < t; k++)
    {
        for (i = 0; i < M; i++)
        {
            for (j = 0; j < L; j++)
            {
                poly_zeroize(&P_coeff[k][i][j]);
            }
        }
    }
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < L; j++)
        {
            poly_copy(&P_coeff[0][i][j], &S[i][j]);
        }
    }

    for (k = 1; k < t; k++)
    {
        for (i = 0; i < M; i++)
        {
            for (j = 0; j < L; j++)
            {
                poly_uniform(&P_coeff[k][i][j], rhoprime, nonce++);
                poly_ntt(&P_coeff[k][i][j]);
                poly_reduce(&P_coeff[k][i][j]);
            }
        }
    }
}

void compute_shr(poly Share[M][L],
                 poly P_coeff[][M][L],
                 int user, int t)
{
    for (int m = 0; m < M; m++)
    {
        for (int l = 0; l < L; l++)
        {
            poly_zeroize(&Share[m][l]);
        }
    }
    int32_t umod = user % Q;
    if (umod < 0)
        umod += Q;

    // pwr = user^0 = 1 initially
    int32_t pwr = 1;

    for (int r = 0; r < t; r++)
    {
        // Convert pwr to Montgomery domain
        // int32_t pwr_mont = to_montgomery(pwr);

        // Scale P_coeff[r] by pwr^r in NTT domain
        for (int m = 0; m < M; m++)
        {
            for (int l = 0; l < L; l++)
            {
                poly tpoly;
                poly_copy(&tpoly, &P_coeff[r][m][l]);
                poly_scale(&tpoly, pwr);
                poly_add(&Share[m][l], &Share[m][l], &tpoly);
                poly_reduce(&Share[m][l]);
            }
        }

        // Now multiply pwr by user for next iteration
        int64_t tmp = (int64_t)pwr * umod;
        tmp %= Q;
        if (tmp < 0)
            tmp += Q;
        pwr = (int32_t)tmp;
    }
}

void reconstruct_secret(poly S[L][K],
                        int32_t *U, int t,
                        poly Shares[][L][K])
{
    int i, j, k, r;
    int32_t u, lambda;
    int64_t val;
    poly tpoly;

    for (i = 0; i < L; i++)
    {
        for (j = 0; j < K; j++)
        {
            poly_zeroize(&S[i][j]);
        }
    }
    for (k = 0; k < t; k++)
    {
        u = U[k];
        lagrange_coeff(&lambda, u, U, t);
        // lambda = to_montgomery(lambda);
        for (i = 0; i < L; i++)
        {
            for (j = 0; j < K; j++)
            {
                poly_copy(&tpoly, &Shares[u - 1][i][j]);
                // Scale by lambda
                poly_scale(&tpoly, lambda);
                poly_add(&S[i][j], &S[i][j], &tpoly);
                poly_reduce(&S[i][j]);
            }
        }
    }
}
