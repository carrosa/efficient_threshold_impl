#include "params.h"

mpz_t GMP_Q, GMP_q, GMP_QINV, GMP_DIV, GMP_MONT, GMP_DIV_QINV, GMP_ROOT_OF_UNITY, GMP_MONT_BITS, GMP_Q_Y, GMP_Q_W, GMP_KAPPA_Y, GMP_KAPPA_W;
const mp_bitcnt_t MPBITCNT_MONT_BITS = MONT_BITS;

void init_params(void)
{
    mpz_inits(GMP_Q, GMP_q, GMP_QINV, GMP_DIV, GMP_MONT, GMP_DIV_QINV, GMP_ROOT_OF_UNITY, GMP_MONT_BITS, GMP_Q_Y, GMP_Q_W, NULL);

    mpz_set_ui(GMP_MONT_BITS, MONT_BITS);
    mpz_set_str(GMP_q, q, 10);
    mpz_set_str(GMP_Q, Q, 10);
    mpz_set_str(GMP_QINV, QINV, 10);
    mpz_set_str(GMP_DIV, DIV, 10);
    mpz_set_str(GMP_MONT, MONT, 10);
    mpz_set_str(GMP_DIV_QINV, DIV_QINV, 10);
    mpz_set_str(GMP_ROOT_OF_UNITY, ROOT_OF_UNITY, 10);
    mpz_set_str(GMP_KAPPA_Y, KAPPA_Y_STR, 10);
    mpz_set_str(GMP_KAPPA_W, KAPPA_W_STR, 10);
}
void free_params(void)
{
    mpz_clears(GMP_Q, GMP_q, GMP_QINV, GMP_DIV, GMP_MONT, GMP_DIV_QINV, GMP_ROOT_OF_UNITY, GMP_MONT_BITS, GMP_Q_Y, GMP_Q_W, GMP_KAPPA_Y, GMP_KAPPA_W, NULL);
}