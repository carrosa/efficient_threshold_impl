#include "lazer.h"
#include "params.h"

/*
This example is to show how matrix-vector multiplication fails when computing t=As
*/

void
prover (uint8_t *proof, polyvec_t s, polymat_t A, polyvec_t t,
        const uint8_t pp[32])
{
  lin_prover_state_t prover;
  lin_prover_init (prover, pp, params);
  lin_prover_set_statement (prover, A, t);
  lin_prover_set_witness (prover, s);
  lin_prover_prove (prover, proof, NULL, NULL);
  lin_prover_clear (prover);
}

int
verifier (const uint8_t *proof, polymat_t A, polyvec_t t, const uint8_t pp[32])
{
  lin_verifier_state_t verifier;
  lin_verifier_init (verifier, pp, params);
  lin_verifier_set_statement (verifier, A, t);
  int accept = lin_verifier_verify (verifier, proof, NULL);
  lin_verifier_clear (verifier);
  return accept;
}

int
main (void)
{
  // Params
  const unsigned int deg = 256;
  const unsigned int mhat = 22;
  const unsigned int lhat = 22;
  const unsigned int khat = 22;
  const unsigned int m = 44;
  const unsigned int n = 88;

  lazer_init ();
  // Set p to 633544653872304603170707504554316801
  INT_T (p, 4);
  INT_T (q, 4);
  const uint8_t p_bytes[] = { 0x01, 0x14, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                              0xf0, 0xc9, 0x15, 0xbf, 0x29, 0x04, 0x7a };
  const uint8_t q_bytes[] = { 0x01, 0x46, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80 };
  memset (p->limbs, 0, sizeof (limb_t) * p->nlimbs);
  memset (q->limbs, 0, sizeof (limb_t) * q->nlimbs);

  int_import (p, p_bytes, sizeof (p_bytes));
  int_import (q, q_bytes, sizeof (q_bytes));
  printf ("p = ");
  int_dump (p); // Should now print: 633544653872304603170707504554316801
  printf ("\n");
  printf ("q = ");
  int_dump (q);
  printf ("\n");

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
  POLYRING_T (Rp, p, deg);
  POLYRING_T (Rq, q, deg);

  // Allocate matrix A = [A1 | A2 | A3]
  polymat_t A, A1, A2, A3;
  polymat_t I_lhat, qI_lhat, Z_lhat;           // For A2
  polymat_t Z_top, I_mhat, qI_mhat, A3_bottom; // For A3

  polymat_alloc (A, Rp, m, n);         // 44x88
  polymat_alloc (A1, Rp, m, khat);     // 44x22
  polymat_alloc (A2, Rp, m, lhat);     // 44x22
  polymat_alloc (A3, Rp, m, 2 * mhat); // 44x44

  polymat_alloc (I_lhat, Rp, lhat, lhat);  // 22x22
  polymat_alloc (qI_lhat, Rp, lhat, lhat); // 22x22
  polymat_alloc (Z_lhat, Rp, mhat, lhat);  // 22x22

  polymat_alloc (Z_top, Rp, lhat, 2 * mhat);     // 22x44
  polymat_alloc (I_mhat, Rp, mhat, mhat);        // 22x22
  polymat_alloc (qI_mhat, Rp, mhat, mhat);       // 22x22
  polymat_alloc (A3_bottom, Rp, mhat, 2 * mhat); // 22x44

  // Allocate secret s = [r | e1 | e2 | msg]
  polyvec_t s, r, e1, e2, msg;
  polyvec_alloc (s, Rp, n);      // 88
  polyvec_alloc (r, Rp, khat);   // 22
  polyvec_alloc (e1, Rp, lhat);  // 22
  polyvec_alloc (e2, Rp, mhat);  // 22
  polyvec_alloc (msg, Rp, mhat); // 44

  // Allocate t = As
  polyvec_t t;
  polyvec_alloc (t, Rp, m); // 44

  // Set A1
  polymat_urandom (A1, p, 119, seed, 0);

  // Set A2
  polymat_set_one (I_lhat);
  polymat_scale (qI_lhat, q, I_lhat);
  polymat_set_zero (Z_lhat);
  for (unsigned int i = 0; i < lhat; i++)
    {
      polyvec_t row;
      polyvec_alloc (row, Rp, qI_lhat->ncols);
      polymat_get_row (row, qI_lhat, i);
      polymat_set_row (A2, row, i);
      polyvec_free (row);
    }
  for (unsigned int i = 0; i < mhat; i++)
    {
      polyvec_t row;
      polyvec_alloc (row, Rp, Z_lhat->ncols);
      polymat_get_row (row, Z_lhat, i);
      polymat_set_row (A2, row, i + lhat);
      polyvec_free (row);
    }

  // Set A3
  polymat_set_zero (Z_top);
  polymat_set_one (I_mhat);
  polymat_scale (qI_mhat, q, I_mhat);
  for (unsigned int i = 0; i < qI_mhat->ncols; i++)
    {
      polyvec_t col;
      polyvec_alloc (col, Rp, qI_mhat->nrows);
      polymat_get_col (col, qI_mhat, i);
      polymat_set_col (A3_bottom, col, i);
      polyvec_free (col);
    }
  for (unsigned int i = 0; i < I_mhat->ncols; i++)
    {
      polyvec_t col;
      polyvec_alloc (col, Rp, I_mhat->nrows);
      polymat_get_col (col, I_mhat, i);
      polymat_set_col (A3_bottom, col, i + qI_mhat->ncols);
      polyvec_free (col);
    }

  for (unsigned int i = 0; i < lhat; i++)
    {
      polyvec_t row;
      polyvec_alloc (row, Rp, Z_top->ncols);
      polymat_get_row (row, Z_top, i);
      polymat_set_row (A3, row, i);
      polyvec_free (row);
    }
  for (unsigned int i = 0; i < mhat; i++)
    {
      polyvec_t row;
      polyvec_alloc (row, Rp, A3_bottom->ncols);
      polymat_get_row (row, A3_bottom, i);
      polymat_set_row (A3, row, i + lhat);
      polyvec_free (row);
    }

  // Set A = [A1 | A2 | A3]
  for (unsigned int i = 0; i < khat; i++) // From A1
    {
      polyvec_t col;
      polyvec_alloc (col, Rp, A1->nrows);
      polymat_get_col (col, A1, i);
      polymat_set_col (A, col, i);
      polyvec_free (col);
    }
  for (unsigned int i = 0; i < lhat; i++) // From A2
    {
      polyvec_t col;
      polyvec_alloc (col, Rp, A2->nrows);
      polymat_get_col (col, A2, i);
      polymat_set_col (A, col, i + khat);
      polyvec_free (col);
    }
  for (unsigned int i = 0; i < 2 * mhat; i++) // From A3
    {
      polyvec_t col;
      polyvec_alloc (col, Rp, A3->nrows);
      polymat_get_col (col, A3, i);
      polymat_set_col (A, col, i + khat + lhat);
      polyvec_free (col);
    }

  // Set s = [r | e1 | e2 | msg]
  polyvec_brandom (r, 2, seed_r, 0);
  polyvec_brandom (e1, 2, seed_e1, 0);
  polyvec_brandom (e2, 2, seed_e2, 0);
  INT_T (lo, 4);
  int_set_zero (lo);
  polyvec_urandom_bnd (msg, lo, q, seed_msg, 0);

  unsigned int idx = 0;
  for (unsigned int i = 0; i < r->nelems; i++)
    {
      polyvec_set_elem (s, idx++, polyvec_get_elem_src (r, i));
    }
  for (unsigned int i = 0; i < e1->nelems; i++)
    {
      polyvec_set_elem (s, idx++, polyvec_get_elem_src (e1, i));
    }
  for (unsigned int i = 0; i < e2->nelems; i++)
    {
      polyvec_set_elem (s, idx++, polyvec_get_elem_src (e2, i));
    }
  for (unsigned int i = 0; i < msg->nelems; i++)
    {
      polyvec_set_elem (s, idx++, polyvec_get_elem_src (msg, i));
    }
  ASSERT_ERR (idx == s->nelems);

  polyvec_redc (s, s);

  // Compute t = A * s
  polyvec_mul (t, A, s);
  polyvec_dump(t); // Segmentation fault here

  // TODO: Prover and verifier

  // Free variables
  polyvec_free (t);
  polyvec_free (s);
  polyvec_free (r);
  polyvec_free (e1);
  polyvec_free (e2);
  polyvec_free (msg);
  polymat_free (A);
  polymat_free (A1);
  polymat_free (A2);
  polymat_free (A3);
  polymat_free (I_lhat);
  polymat_free (qI_lhat);
  polymat_free (Z_lhat);
  polymat_free (Z_top);
  polymat_free (I_mhat);
  polymat_free (qI_mhat);
  polymat_free (A3_bottom);

  return 1;
}