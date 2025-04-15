#ifndef POLY256_H
#define POLY256_H

#include <gmp.h>
#include <stdint.h>
#include <stdarg.h>
#include "params.h"
#include "ntt256.h"
#include "reduce256.h"
#include "symmetric.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{
    mpz_t coeffs[N];
} poly;

void H1(
    poly arr[K],
    uint8_t *hash);

void H0(
    poly arr[K],
    uint8_t *hash);

void H0_matrix(poly *arr, size_t rows, size_t cols, uint8_t *hash);

int polycmp(const poly *p1, const poly *p2);
void poly_init(poly *p);
void poly_1d_init(poly *arr, size_t count);
void poly_inits(poly *p, ...);

void poly_clear(poly *p);
void poly_1d_clear(poly *arr, size_t count);
void poly_clears(poly *p, ...);

#define POLY_2D_INIT(arr, d1, d2) poly_1d_init((poly *)(&(arr)[0][0]), (d1) * (d2));
#define POLY_2D_CLEAR(arr, d1, d2) poly_1d_clear((poly *)(&(arr)[0][0]), (d1) * (d2));
#define POLY_3D_INIT(arr, d1, d2, d3) poly_1d_init((poly *)(&(arr)[0][0]), (d1) * (d2) * (d3));
#define POLY_3D_CLEAR(arr, d1, d2, d3) poly_1d_clear((poly *)(&(arr)[0][0]), (d1) * (d2) * (d3));


#define POLY_4D_INIT(ptr, d1, d2, d3, d4) \
    poly_1d_init((poly *)(ptr), (d1) * (d2) * (d3) * (d4))

#define POLY_4D_CLEAR(ptr, d1, d2, d3, d4) \
    poly_1d_clear((poly *)(ptr), (d1) * (d2) * (d3) * (d4))




#define POLY_2D_INIT_SAFE(arr, d1, d2) \
    for (int i = 0; i < (d1); i++)     \
        for (int j = 0; j < (d2); j++) \
    poly_init(&(arr)[i][j])

void poly_reduce(poly *a);
void poly_caddq(poly *a);

void poly_add(poly *c, const poly *a, const poly *b);
void poly_sub(poly *c, const poly *a, const poly *b);
void poly_shiftl(poly *a);

void poly_ntt(poly *a);
void poly_invntt_tomont(poly *a);
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b);

void poly_decompose(poly *a1, poly *a0, const poly *a);
unsigned int poly_make_hint(poly *h, const poly *a0, const poly *a1);
void poly_use_hint(poly *b, const poly *a, const poly *h);

int poly_chknorm(const poly *a, int32_t B);
void poly_uniform(poly *a,
                  const uint8_t seed[SEEDBYTES],
                  uint16_t nonce);
void poly_uniform_eta(poly *a,
                      const uint8_t seed[CRHBYTES],
                      uint16_t nonce);

void poly_copy(poly *dst, const poly *src);

void poly_scale(poly *r, mpz_t c);

void poly_mod(poly *r, mpz_t c);

int poly_equal(poly *a, poly *b);

void poly_challenge(poly *c, const uint8_t seed[CTILDEBYTES]);

void poly_hash_shake256(
    poly *p,
    uint8_t *hash,
    size_t outlen);

void poly_arr_hash_shake256(
    const poly *arr,
    size_t len,
    uint8_t *hash,
    size_t outlen);

void poly_matrix_hash_shake256(
    const poly *mat,
    size_t rows,
    size_t cols,
    uint8_t *hash,
    size_t outlen);

void get_random_dummy_poly(poly *p);

#endif
