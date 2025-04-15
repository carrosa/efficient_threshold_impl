#ifndef NTT256_H
#define NTT256_H

#include "params.h"
// Definitions, to be moved later
// #define N 256

#include <stdint.h>
#include <gmp.h>

void init_zetas(void);
void clear_zetas(void);

void ntt(mpz_t a[N]);

void invntt_tomont(mpz_t a[N]);

#endif
