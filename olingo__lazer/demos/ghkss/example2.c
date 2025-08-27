#include "lazer.h"
#include "params.h"


/*
This example is to show how polynomial multiplication fails when computing c=ab
We have tried using different NTT-friendly primes like our original one from example1 and 12289.
*/

int
main (void)
{

  INT_T (p, 2);
  int_set_i64 (p, 12289);

  // Derive seeds using shake128
  shake128_state_t state;
  uint8_t seed[32] = { 0 };
  uint8_t seed_r[32], seed_e1[32], seed_e2[32], seed_msg[32];
  memcpy (seed_r, seed, 32);
  uint8_t tmp[34];
  memcpy (tmp, seed, 32);
  tmp[32] = 'e';
  tmp[33] = '1';
  shake128_init (state);
  shake128_absorb (state, tmp, 34);
  shake128_squeeze (state, seed_e1, 32);
  shake128_clear (state);

  memcpy (tmp, seed, 32);
  tmp[32] = 'e';
  tmp[33] = '2';
  shake128_init (state);
  shake128_absorb (state, tmp, 34);
  shake128_squeeze (state, seed_e2, 32);
  shake128_clear (state);

  memcpy (tmp, seed, 32);
  tmp[32] = 'm';
  tmp[33] = '\0';
  shake128_init (state);
  shake128_absorb (state, tmp, 34);
  shake128_squeeze (state, seed_msg, 32);
  shake128_clear (state);

  // Define ring Rq and Rp
  POLYRING_T (Rp, p, 256);

  // Attemt polynomial multiplication
  poly_t a, b, c;
  poly_alloc (a, Rp);
  poly_alloc (b, Rp);
  poly_alloc (c, Rp);

  poly_urandom (a, p, 15, seed_r, 0);
  poly_urandom (b, p, 15, seed_e1, 0);
  poly_dump (a);
  poly_dump (b);
  poly_mul (c, a, b); // NTT seems to fail
  poly_dump (c);

  // Free variables
  poly_free(a);
  poly_free(b);
  poly_free(c);

  return 1;
}