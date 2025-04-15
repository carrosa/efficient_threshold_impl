#ifndef REDUCE_H
#define REDUCE_H

#include <gmp.h>
#include "params.h"


void montgomery_reduce(mpz_t r, mpz_t a);

void reduce32(mpz_t r, mpz_t a);

void caddq(mpz_t r, const mpz_t a);

void freeze(mpz_t r, const mpz_t a);

#endif
