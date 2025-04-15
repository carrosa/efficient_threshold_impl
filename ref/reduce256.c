#include "reduce256.h"

/*************************************************
 * Name:        montgomery_reduce
 *
 * Description: For a finite field element a (stored in an mpz_t)
 *              with 0 <= a < Q*R, compute r = a * R^{-1} mod Q,
 *              using Montgomery reduction with R = 2^(R_BITS).
 *
 *              The algorithm is as follows:
 *                1. t = (a * QINV) mod R, where R = 2^(R_BITS)
 *                2. r = (a + t * Q) / R
 *                3. If r >= Q then subtract Q to get a proper representative.
 *
 * Arguments:
 *   - mpz_t r: output (the reduced result)
 *   - const mpz_t a: input element (should be in [0, Q*R) )
 *
 * Returns:
 *   None (result is stored in r)
 **************************************************/
void montgomery_reduce(mpz_t r, mpz_t a)
{
  mpz_t t, tmp;
  mpz_inits(t, tmp, NULL);

  // t = (a * QINV) mod 2^MONT_BITS
  mpz_mod_2exp(t, a, MONT_BITS);
  mpz_mul(t, t, GMP_QINV);

  // t = t*Q
  mpz_mul(tmp, t, GMP_Q);

  // r = (a - tmp) >> MONT_BITS
  mpz_sub(r, a, tmp);
  mpz_fdiv_q_2exp(r, r, MPBITCNT_MONT_BITS);

  // Ensure result is in [0, Q-1] // It should be in this range now anyway, I don't think this is required.
  // if (mpz_cmp(r, GMP_Q) >= 0)
  //   mpz_sub(r, r, GMP_Q);

  mpz_clears(t, tmp, NULL);
}

/*************************************************
 * Name:        caddq
 *
 * Description: If the input a is negative (stored in an mpz_t),
 *              add Q so that the result is nonnegative.
 *
 * Arguments:
 *   - mpz_t r: output (the result after conditionally adding Q)
 *   - const mpz_t a: input element
 *
 * Returns:
 *   None (result is stored in r)
 **************************************************/
void caddq(mpz_t r, const mpz_t a)
{
  if (mpz_sgn(a) < 0)
  {
    mpz_add(r, a, GMP_Q);
  }
  else
  {
    mpz_set(r, a);
  }
}

/*************************************************
 * Name:        freeze
 *
 * Description: For a finite field element a (stored in an mpz_t),
 *              compute its canonical representative r = a mod Q.
 *              Optionally, this routine can produce a symmetric
 *              representative in the range [ -Q/2, Q/2 ].
 *
 * Arguments:
 *   - mpz_t r: output (the canonical representative)
 *   - const mpz_t a: input element
 *
 * Returns:
 *   None (result is stored in r)
 **************************************************/
void freeze(mpz_t r, const mpz_t a)
{
  mpz_mod(r, a, GMP_Q);
  {
    mpz_t half;
    mpz_init(half);
    mpz_fdiv_q_ui(half, GMP_Q, 2); // half = floor(Q/2)
    if (mpz_cmp(r, half) >= 0)
    {
      mpz_sub(r, r, GMP_Q);
    }
    mpz_clear(half);
  }
}
