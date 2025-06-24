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
int mont_initialized = 0;
mp_limb_t mont_qinv;
mp_limb_t mont_q_limbs[3];
mp_size_t mont_q_size;
mp_bitcnt_t mont_bits = 128;

void poly_init_mont(void)
{
  if (!mont_initialized)
  {
    mont_qinv = mpz_get_ui(GMP_QINV);
    mont_q_size = GMP_Q->_mp_size;
    mpn_copyi(mont_q_limbs, GMP_Q->_mp_d, mont_q_size);
    mont_initialized = 1;
  }
}

void montgomery_reduce_old(mpz_t r, mpz_t a, mpz_t t, mpz_t tmp)
{
  mp_size_t a_size = a->_mp_size;
  mp_size_t abs_a_size = __GMP_ABS(a_size);
  mp_size_t mod_size = mont_bits / GMP_NUMB_BITS;
  if (abs_a_size > mod_size)
  {
    abs_a_size = mod_size;
  }
  mp_ptr t_limbs = t->_mp_d;
}

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

int barrett_initialized = 0;
mpz_t barrett_m, q_half, neg_q_half;
unsigned int barrett_shift;
void init_barrett(void)
{
  if (!barrett_initialized)
  {
    mpz_t barrett_k;
    mpz_inits(
        barrett_m,
        barrett_k,
        q_half,
        neg_q_half,
        NULL);
    unsigned int k_bits = 119; // Hardcoded log2(Q)
    mpz_set_ui(barrett_k, 1);
    mpz_mul_2exp(barrett_k, barrett_k, 2 * k_bits);
    mpz_fdiv_q(barrett_m, barrett_k, GMP_Q); // m = floor(2^(2*k) / Q)
    barrett_shift = 2 * k_bits;              // 238
    mpz_clear(barrett_k);

    mpz_fdiv_q_ui(q_half, GMP_Q, 2); // Q/2
    mpz_neg(neg_q_half, q_half);     // -Q/2

    barrett_initialized = 1;
  }
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
