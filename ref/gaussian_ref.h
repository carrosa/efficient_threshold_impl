#include <stdint.h>
#include "fastrandombytes.h"
#include <x86intrin.h>
#include "cpucycles.h"
#include <stdio.h>
#include "randombytes.h"
#include <math.h>
#include <gmp.h>

void gaussian_sampler(mpz_t *sample, uint32_t slen);
void gaussian_sampler_scaled(mpz_t *sample, uint32_t slen, uint64_t target_sigma);
// void gaussian_sampler_scaled(mpz_t *sample, uint32_t slen, const mpz_t target_sigma);
void gaussian_sampler_y(mpz_t *s, uint32_t len); // sigma_s
void gaussian_sampler_w(mpz_t *sample, uint32_t len); // sigma_r
void gaussian_sampler_drown(mpz_t *s, uint32_t len); // sigma_tdec
// void gaussian_sampler_bin(mpz_t *s, uint32_t len); //(subscript ctx) sigma sqrt 2/3
void binary_sampler(mpz_t *s, uint32_t len);
