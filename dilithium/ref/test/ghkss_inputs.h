#ifndef GHKSS_INPUTS_H
#define GHKSS_INPUTS_H

#include "../params.h"
#include "../GHKSS256.h"
#include "../poly256.h"
#include "stdlib.h"
#include "string.h"

// =====================
// KEYGEN
// =====================
typedef struct
{
    poly (*A)[LHAT]; // size: [KHAT][LHAT]
    poly (*B)[M];    // size: [KHAT][M]
    poly (*S)[M];    // size: [LHAT][M]
} keygen_input_t;

void init_keygen_input(keygen_input_t *in)
{
    in->A = malloc(sizeof(poly[KHAT][LHAT]));
    in->B = malloc(sizeof(poly[KHAT][M]));
    in->S = malloc(sizeof(poly[LHAT][M]));
    POLY_2D_INIT(in->A, KHAT, LHAT);
    POLY_2D_INIT(in->B, KHAT, M);
    POLY_2D_INIT(in->S, LHAT, M);
}

void clear_keygen_input(keygen_input_t *in)
{
    POLY_2D_CLEAR(in->A, KHAT, LHAT);
    POLY_2D_CLEAR(in->B, KHAT, M);
    POLY_2D_CLEAR(in->S, LHAT, M);
    free(in->A);
    free(in->B);
    free(in->S);
}

// =====================
// Encryption
// =====================
typedef struct
{
    poly *U;         // [LHAT]
    poly *V;         // [M]
    poly (*A)[LHAT]; // [KHAT][LHAT]
    poly (*B)[M];    // [KHAT][M]
    poly *Msg;       // [M]
} encrypt_input_t;

void init_encrypt_input(encrypt_input_t *in)
{
    in->U = malloc(sizeof(poly[LHAT]));
    in->V = malloc(sizeof(poly[M]));
    in->A = malloc(sizeof(poly[KHAT][LHAT]));
    in->B = malloc(sizeof(poly[KHAT][M]));
    in->Msg = malloc(sizeof(poly[M]));

    if (!in->U || !in->V || !in->A || !in->B || !in->Msg)
    {
        fprintf(stderr, "Memory allocation failed in init_encrypt_input\n");
        exit(1);
    }

    poly_1d_init(in->U, LHAT);
    poly_1d_init(in->V, M);
    POLY_2D_INIT(in->A, KHAT, LHAT);
    POLY_2D_INIT(in->B, KHAT, M);
    poly_1d_init(in->Msg, M);
}

void clear_encrypt_input(encrypt_input_t *in)
{
    poly_1d_clear(in->U, LHAT);
    poly_1d_clear(in->V, M);
    POLY_2D_CLEAR(in->A, KHAT, LHAT);
    POLY_2D_CLEAR(in->B, KHAT, M);
    poly_1d_clear(in->Msg, M);

    free(in->U);
    free(in->V);
    free(in->A);
    free(in->B);
    free(in->Msg);
}

void sample_encrypt_input(encrypt_input_t *in)
{
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    gen_randomness(seed);

    // Sample A
    for (int i = 0; i < KHAT; i++)
    {
        for (int j = 0; j < LHAT; j++)
        {
            poly_uniform(&in->A[i][j], seed, nonce++);
            poly_ntt(&in->A[i][j]);
            poly_reduce(&in->A[i][j]);
        }
    }

    // Sample B
    for (int i = 0; i < KHAT; i++)
    {
        for (int j = 0; j < M; j++)
        {
            gaussian_sampler(in->B[i][j].coeffs, N);
            poly_ntt(&in->B[i][j]);
            poly_reduce(&in->B[i][j]);
        }
    }

    // Sample Msg
    for (int i = 0; i < M; i++)
    {
        gaussian_sampler(in->Msg[i].coeffs, N);
        poly_ntt(&in->Msg[i]);
        poly_reduce(&in->Msg[i]);
    }
}

// =====================
// DKG
// =====================

typedef struct
{
    // DKG 1 inputs
    poly (*Ae)[LHAT]; // size: [KHAT][LHAT]
    poly (*Si)[M];    // size: [LHAT][M]
    poly (*Ei)[M];    // size: [KHAT][M]
    poly (*Bi)[M];    // size: [KHAT][M]

    // DKG 1 outputs
    uint8_t *h_Bi; // 32-byte output hash for dk_gen_1

    // DKG 2 inputs
    // Si, Ei, Ae are reused from 1
    int32_t *users;

    // DKG 3 inputs
    poly (*Bj)[KHAT][M];  // size: [USERS][KHAT][M]
    poly (*Sij)[LHAT][M]; // size: [USERS][LHAT][M]
    poly (*Be)[M];        // size: [KHAT][M]
    poly (*Sei)[M];       // size: [LHAT][M]
} dkg_input_t;

void init_dkg(dkg_input_t *in)
{
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    gen_randomness(seed);

    in->Ae = malloc(sizeof(poly) * KHAT * LHAT);
    in->Si = malloc(sizeof(poly) * LHAT * M);
    in->Ei = malloc(sizeof(poly) * KHAT * M);
    in->Bi = malloc(sizeof(poly) * KHAT * M);
    in->users = malloc(sizeof(int32_t) * USERS);
    in->Bj = malloc(sizeof(poly) * USERS * KHAT * M);
    in->Sij = malloc(sizeof(poly) * USERS * LHAT * M);
    in->Be = malloc(sizeof(poly) * KHAT * M);
    in->Sei = malloc(sizeof(poly) * LHAT * M);

    POLY_2D_INIT(in->Ae, KHAT, LHAT);
    POLY_2D_INIT(in->Si, LHAT, M);
    POLY_2D_INIT(in->Ei, KHAT, M);
    POLY_2D_INIT(in->Bi, KHAT, M);
    POLY_3D_INIT(in->Bj, USERS, KHAT, M);
    POLY_3D_INIT(in->Sij, USERS, LHAT, M);
    POLY_2D_INIT(in->Be, KHAT, M);
    POLY_2D_INIT(in->Sei, LHAT, M);

    // Sample Ae
    for (int i = 0; i < KHAT; i++)
    {
        for (int j = 0; j < LHAT; j++)
        {
            poly_uniform(&in->Ae[i][j], seed, nonce++);
            poly_ntt(&in->Ae[i][j]);
            poly_reduce(&in->Ae[i][j]);
        }
    }
    // Set users
    for (int i = 0; i < USERS; i++)
    {
        in->users[i] = i + 1;
    }

    // Used for DKG3
    for (int i = 0; i < USERS; i++)
    {
        for (int j = 0; j < M; j++)
        {
            for (int k = 0; k < KHAT; k++)
            {
                poly_uniform(&in->Bj[i][k][j], seed, nonce++);
                poly_reduce(&in->Bj[i][k][j]);
            }
            for (int k = 0; k < LHAT; k++)
            {
                poly_uniform(&in->Sij[i][k][j], seed, nonce++);
                poly_reduce(&in->Sij[i][k][j]);
            }
        }
    }
}

// Used inside DKG round 1 benchmark
void init_dkg1(dkg_input_t *in)
{
    in->h_Bi = malloc(sizeof(uint8_t[32]));
}


void clear_dkg(dkg_input_t *in)
{
    POLY_2D_CLEAR(in->Ae, KHAT, LHAT);
    POLY_2D_CLEAR(in->Si, LHAT, M);
    POLY_2D_CLEAR(in->Ei, KHAT, M);
    POLY_2D_CLEAR(in->Bi, KHAT, M);
    POLY_3D_CLEAR(in->Bj, USERS, KHAT, M);
    POLY_3D_CLEAR(in->Sij, USERS, LHAT, M);
    POLY_2D_CLEAR(in->Be, KHAT, M);
    POLY_2D_CLEAR(in->Sei, LHAT, M);

    free(in->Ae);
    free(in->Si);
    free(in->Ei);
    free(in->Bi);
    free(in->h_Bi);
    free(in->users);
    free(in->Bj);
    free(in->Sij);
    free(in->Be);
    free(in->Sei);
}

// ======================
// SIGNATURE
// ======================

typedef struct
{
    poly (*A)[L];    // size: [K][L]
    poly (*Ae)[LHAT]; // size: [KHAT][LHAT]
    poly (*Be)[M];    // size: [KHAT][M]
    poly *u;          // size: [LHAT]
    poly *v;          // size: [M]

    
} sig_input_t;

void init_sig(sig_input_t *in)
{
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    gen_randomness(seed);

    in->A = malloc(sizeof(poly[K][L]));
    in->Ae = malloc(sizeof(poly[KHAT][LHAT]));
    in->Be = malloc(sizeof(poly[KHAT][M]));
    in->u = malloc(sizeof(poly[LHAT]));
    in->v = malloc(sizeof(poly[M]));

    POLY_2D_INIT(in->A, K, L);
    POLY_2D_INIT(in->Ae, KHAT, LHAT);
    POLY_2D_INIT(in->Be, KHAT, M);
    poly_1d_init(in->u, LHAT);
    poly_1d_init(in->v, M);

    // Populate dummy data
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < L; j++)
        {
            poly_uniform(&in->A[i][j], seed, nonce++);
            poly_ntt(&in->A[i][j]);
            poly_reduce(&in->A[i][j]);
        }
    }

    for (int i = 0; i < KHAT; i++)
    {
        for (int j = 0; j < LHAT; j++)
        {
            poly_uniform(&in->Ae[i][j], seed, nonce++);
            poly_ntt(&in->Ae[i][j]);
            poly_reduce(&in->Ae[i][j]);
        }
    }

    for (int i = 0; i < KHAT; i++)
    {
        for (int j = 0; j < M; j++)
        {
            poly_uniform(&in->Be[i][j], seed, nonce++);
            poly_ntt(&in->Be[i][j]);
            poly_reduce(&in->Be[i][j]);
        }
    }
}

void clear_sig(sig_input_t *in)
{
    POLY_2D_CLEAR(in->A, K, L);
    POLY_2D_CLEAR(in->Ae, KHAT, LHAT);
    POLY_2D_CLEAR(in->Be, KHAT, M);
    poly_1d_clear(in->u, LHAT);
    poly_1d_clear(in->v, M);

    free(in->A);
    free(in->Ae);
    free(in->Be);
    free(in->u);
    free(in->v);
}

#endif