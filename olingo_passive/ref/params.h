#ifndef PARAMS_H
#define PARAMS_H

#include <gmp.h>
#include "config.h"

#define SEEDBYTES 64 // 64 // 32
#define CRHBYTES 64
#define CRYPTO_BYTES 32
#define TRBYTES 64
#define RNDBYTES 32
#define D 13
#define SAMPLEBYTES 15

#define THRESHOLD 127
#define USERS 128

#define N 256
#define MONT_BITS 128
#define q "562949953438721"

#define Q "169754086842540358005352209309511681"
#define QINV "-133978373632131801599285206084969574399"
#define DIV "-25978460312424105965743875580309239"
#define MONT "-74577198354954337356572233802708949"
#define DIV_QINV "57740766783911453634293358172553945353"
#define ROOT_OF_UNITY "75574343650629980472360693068435123"

extern mpz_t GMP_Q;
extern mpz_t GMP_q;
extern mpz_t GMP_QINV;
extern mpz_t GMP_DIV;
extern mpz_t GMP_MONT;
extern mpz_t GMP_DIV_QINV;
extern mpz_t GMP_ROOT_OF_UNITY;
extern mpz_t GMP_MONT_BITS;
extern const mp_bitcnt_t MPBITCNT_MONT_BITS;

void init_params(void);
void free_params(void);

#define K 6
#define L 6
#define KHAT 22
#define LHAT 22
#define M 22

#define KAPPA_Y 20
#define KAPPA_Y_STR "1048576"
#define KAPPA_W 38
#define KAPPA_W_STR "274877906944"
extern mpz_t GMP_KAPPA_Y;
extern mpz_t GMP_KAPPA_W;

#define TAIL_BOUND 4.8
#define TAIL_BOUND_2 1.34

#define SIGMA_BINARY sqrt(2/3)
#define SIGMA_SMALL (1 << 20)
#define SIGMA_LARGE (1 << 38)

#define CTILDEBYTES 32
#define TAU 60

#define ETA 2

// Benchmarking
#define BENCHMARK_ITERATION 1
#endif
