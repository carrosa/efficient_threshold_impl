#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "params.h"
#include "indcpa.h"
#include "polyvec.h"
#include "poly.h"
#include "ntt.h"
#include "symmetric.h"
#include "randombytes.h"
#include "fips202.h"

/*************************************************
 * Name:        pack_pk
 *
 * Description: Serialize the public key as concatenation of the
 *              serialized vector of polynomials pk
 *              and the public seed used to generate the matrix A.
 *
 * Arguments:   uint8_t *r: pointer to the output serialized public key
 *              polyvec *pk: pointer to the input public-key polyvec
 *              const uint8_t *seed: pointer to the input public seed
 **************************************************/
static void pack_pk(uint8_t r[KYBER_INDCPA_PUBLICKEYBYTES],
                    polyvec *pk,
                    const uint8_t seed[KYBER_SYMBYTES])
{
  polyvec_tobytes(r, pk);
  memcpy(r + KYBER_POLYVECBYTES, seed, KYBER_SYMBYTES);
}

/*************************************************
 * Name:        unpack_pk
 *
 * Description: De-serialize public key from a byte array;
 *              approximate inverse of pack_pk
 *
 * Arguments:   - polyvec *pk: pointer to output public-key polynomial vector
 *              - uint8_t *seed: pointer to output seed to generate matrix A
 *              - const uint8_t *packedpk: pointer to input serialized public key
 **************************************************/
static void unpack_pk(polyvec *pk,
                      uint8_t seed[KYBER_SYMBYTES],
                      const uint8_t packedpk[KYBER_INDCPA_PUBLICKEYBYTES])
{
  polyvec_frombytes(pk, packedpk);
  memcpy(seed, packedpk + KYBER_POLYVECBYTES, KYBER_SYMBYTES);
}

/*************************************************
 * Name:        pack_sk
 *
 * Description: Serialize the secret key
 *
 * Arguments:   - uint8_t *r: pointer to output serialized secret key
 *              - polyvec *sk: pointer to input vector of polynomials (secret key)
 **************************************************/
static void pack_sk(uint8_t r[KYBER_INDCPA_SECRETKEYBYTES], polyvec *sk)
{
  polyvec_tobytes(r, sk);
}

/*************************************************
 * Name:        unpack_sk
 *
 * Description: De-serialize the secret key; inverse of pack_sk
 *
 * Arguments:   - polyvec *sk: pointer to output vector of polynomials (secret key)
 *              - const uint8_t *packedsk: pointer to input serialized secret key
 **************************************************/
static void unpack_sk(polyvec *sk, const uint8_t packedsk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec_frombytes(sk, packedsk);
}

/*************************************************
 * Name:        pack_ciphertext
 *
 * Description: Serialize the ciphertext as concatenation of the
 *              compressed and serialized vector of polynomials b
 *              and the compressed and serialized polynomial v
 *
 * Arguments:   uint8_t *r: pointer to the output serialized ciphertext
 *              poly *pk: pointer to the input vector of polynomials b
 *              poly *v: pointer to the input polynomial v
 **************************************************/
static void pack_ciphertext(uint8_t r[KYBER_INDCPA_BYTES], polyvec *b, poly *v)
{
  polyvec_compress(r, b);
  poly_compress(r + KYBER_POLYVECCOMPRESSEDBYTES, v);
}

/*************************************************
 * Name:        unpack_ciphertext
 *
 * Description: De-serialize and decompress ciphertext from a byte array;
 *              approximate inverse of pack_ciphertext
 *
 * Arguments:   - polyvec *b: pointer to the output vector of polynomials b
 *              - poly *v: pointer to the output polynomial v
 *              - const uint8_t *c: pointer to the input serialized ciphertext
 **************************************************/
static void unpack_ciphertext(polyvec *b, poly *v, const uint8_t c[KYBER_INDCPA_BYTES])
{
  polyvec_decompress(b, c);
  poly_decompress(v, c + KYBER_POLYVECCOMPRESSEDBYTES);
}

/*************************************************
 * Name:        rej_uniform
 *
 * Description: Run rejection sampling on uniform random bytes to generate
 *              uniform random integers mod q
 *
 * Arguments:   - int16_t *r: pointer to output buffer
 *              - unsigned int len: requested number of 16-bit integers (uniform mod q)
 *              - const uint8_t *buf: pointer to input buffer (assumed to be uniformly random bytes)
 *              - unsigned int buflen: length of input buffer in bytes
 *
 * Returns number of sampled 16-bit integers (at most len)
 **************************************************/
static unsigned int rej_uniform(int16_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  uint16_t val0, val1;

  ctr = pos = 0;
  while (ctr < len && pos + 3 <= buflen)
  {
    val0 = ((buf[pos + 0] >> 0) | ((uint16_t)buf[pos + 1] << 8)) & 0xFFF;
    val1 = ((buf[pos + 1] >> 4) | ((uint16_t)buf[pos + 2] << 4)) & 0xFFF;
    pos += 3;

    if (val0 < KYBER_Q)
      r[ctr++] = val0;
    if (ctr < len && val1 < KYBER_Q)
      r[ctr++] = val1;
  }

  return ctr;
}

#define gen_a(A, B) gen_matrix(A, B, 0)
#define gen_at(A, B) gen_matrix(A, B, 1)

/*************************************************
 * Name:        gen_matrix
 *
 * Description: Deterministically generate matrix A (or the transpose of A)
 *              from a seed. Entries of the matrix are polynomials that look
 *              uniformly random. Performs rejection sampling on output of
 *              a XOF
 *
 * Arguments:   - polyvec *a: pointer to ouptput matrix A
 *              - const uint8_t *seed: pointer to input seed
 *              - int transposed: boolean deciding whether A or A^T is generated
 **************************************************/
#if (XOF_BLOCKBYTES % 3)
#error "Implementation of gen_matrix assumes that XOF_BLOCKBYTES is a multiple of 3"
#endif

#define GEN_MATRIX_NBLOCKS ((12 * KYBER_N / 8 * (1 << 12) / KYBER_Q + XOF_BLOCKBYTES) / XOF_BLOCKBYTES)
// Not static for benchmarking
void gen_matrix(polyvec *a, const uint8_t seed[KYBER_SYMBYTES], int transposed)
{
  unsigned int ctr, i, j;
  unsigned int buflen;
  uint8_t buf[GEN_MATRIX_NBLOCKS * XOF_BLOCKBYTES];
  xof_state state;

  for (i = 0; i < MAT_COLS; i++)
  {
    for (j = 0; j < MAT_COLS; j++)
    {
      if (transposed)
        xof_absorb(&state, seed, i, j);
      else
        xof_absorb(&state, seed, j, i);

      xof_squeezeblocks(buf, GEN_MATRIX_NBLOCKS, &state);
      buflen = GEN_MATRIX_NBLOCKS * XOF_BLOCKBYTES;
      ctr = rej_uniform(a[i].vec[j].coeffs, KYBER_N, buf, buflen);

      while (ctr < KYBER_N)
      {
        xof_squeezeblocks(buf, 1, &state);
        buflen = XOF_BLOCKBYTES;
        ctr += rej_uniform(a[i].vec[j].coeffs + ctr, KYBER_N - ctr, buf, buflen);
      }
    }
  }
}
void print_buffer(uint8_t *buf, size_t length)
{
  for (size_t i = 0; i < length; i++)
  {
    printf("%u ", buf[i]); // Print each byte in hexadecimal, two digits with leading zero
    if ((i + 1) % 16 == 0) // Print 16 bytes per line for readability
      printf("\n");
  }
  printf("\n\n\n");
}

void gen_matrix_c(poly A[MAT_m][MAT_m], const uint8_t seed[MAT_COINS], int transposed)
{
  unsigned int ctr, i, j;
  unsigned int buflen;
  uint8_t buf[GEN_MATRIX_NBLOCKS * XOF_BLOCKBYTES];
  xof_state state;

  for (i = 0; i < MAT_m; i++)
  {
    for (j = 0; j < MAT_m; j++)
    {
      if (transposed)
        xof_absorb(&state, seed, i, j);
      else
        xof_absorb(&state, seed, j, i);

      xof_squeezeblocks(buf, GEN_MATRIX_NBLOCKS, &state);
      buflen = GEN_MATRIX_NBLOCKS * XOF_BLOCKBYTES;
      ctr = rej_uniform(A[i][j].coeffs, KYBER_N, buf, buflen);

      while (ctr < KYBER_N)
      {
        xof_squeezeblocks(buf, 1, &state);
        buflen = XOF_BLOCKBYTES;
        ctr += rej_uniform(A[i][j].coeffs + ctr, KYBER_N - ctr, buf, buflen);
      }
    }
  }
}

/*************************************************
 * Name:        indcpa_keypair_derand
 *
 * Description: Generates public and private key for the CPA-secure
 *              public-key encryption scheme underlying Kyber
 *
 * Arguments:   - uint8_t *pk: pointer to output public key
 *                             (of length KYBER_INDCPA_PUBLICKEYBYTES bytes)
 *              - uint8_t *sk: pointer to output private key
 *                             (of length KYBER_INDCPA_SECRETKEYBYTES bytes)
 *              - const uint8_t *coins: pointer to input randomness
 *                             (of length KYBER_SYMBYTES bytes)
 **************************************************/
void indcpa_keypair_derand(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                           uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES],
                           const uint8_t coins[KYBER_SYMBYTES])
{
  unsigned int i;
  uint8_t buf[2 * KYBER_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf + KYBER_SYMBYTES;
  uint8_t nonce = 0;
  polyvec a[KYBER_K], e, pkpv, skpv;

  memcpy(buf, coins, KYBER_SYMBYTES);
  buf[KYBER_SYMBYTES] = KYBER_K;
  hash_g(buf, buf, KYBER_SYMBYTES + 1);

  gen_a(a, publicseed);

  for (i = 0; i < KYBER_K; i++)
    poly_getnoise_eta1(&skpv.vec[i], noiseseed, nonce++);
  for (i = 0; i < KYBER_K; i++)
    poly_getnoise_eta1(&e.vec[i], noiseseed, nonce++);

  polyvec_ntt(&skpv);
  polyvec_ntt(&e);

  // matrix-vector multiplication
  for (i = 0; i < KYBER_K; i++)
  {
    polyvec_basemul_acc_montgomery(&pkpv.vec[i], &a[i], &skpv);
    poly_tomont(&pkpv.vec[i]);
  }

  polyvec_add(&pkpv, &pkpv, &e);
  polyvec_reduce(&pkpv);

  pack_sk(sk, &skpv);
  pack_pk(pk, &pkpv, publicseed);
}

/*************************************************
 * Name:        indcpa_enc
 *
 * Description: Encryption function of the CPA-secure
 *              public-key encryption scheme underlying Kyber.
 *
 * Arguments:   - uint8_t *c: pointer to output ciphertext
 *                    poly_reduce(r);              (of length KYBER_INDCPA_MSGBYTES bytes)
 *              - const uint8_t *pk: pointer to input public key
 *                                   (of length KYBER_INDCPA_PUBLICKEYBYTES)
 *              - const uint8_t *coins: pointer to input random coins used as seed
 *                                      (of length KYBER_SYMBYTES) to deterministically
 *                                      generate all randomness
 **************************************************/
void indcpa_enc(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES])
{
  unsigned int i;
  uint8_t seed[KYBER_SYMBYTES];
  uint8_t nonce = 0;
  polyvec sp, pkpv, ep, at[KYBER_K], b;
  poly v, k, epp;

  unpack_pk(&pkpv, seed, pk);
  poly_frommsg(&k, m);
  gen_at(at, seed);

  for (i = 0; i < KYBER_K; i++)
    poly_getnoise_eta1(sp.vec + i, coins, nonce++);
  for (i = 0; i < KYBER_K; i++)
    poly_getnoise_eta2(ep.vec + i, coins, nonce++);
  poly_getnoise_eta2(&epp, coins, nonce++);

  polyvec_ntt(&sp);

  // matrix-vector multiplication
  for (i = 0; i < KYBER_K; i++)
    polyvec_basemul_acc_montgomery(&b.vec[i], &at[i], &sp);

  polyvec_basemul_acc_montgomery(&v, &pkpv, &sp);

  polyvec_invntt_tomont(&b);
  poly_invntt_tomont(&v);

  polyvec_add(&b, &b, &ep);
  poly_add(&v, &v, &epp);
  poly_add(&v, &v, &k);
  polyvec_reduce(&b);
  poly_reduce(&v);

  pack_ciphertext(c, &b, &v);
}

/*************************************************
 * Name:        indcpa_dec
 *
 * Description: Decryption function of the CPA-secure
 *              public-key encryption scheme underlying Kyber.
 *
 * Arguments:   - uint8_t *m: pointer to output decrypted message
 *                            (of length KYBER_INDCPA_MSGBYTES)
 *              - const uint8_t *c: pointer to input ciphertext
 *                                  (of length KYBER_INDCPA_BYTES)
 *              - const uint8_t *sk: pointer to input secret key
 *                                   (of length KYBER_INDCPA_SECRETKEYBYTES)
 **************************************************/
void indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec b, skpv;
  poly v, mp;

  unpack_ciphertext(&b, &v, c);
  unpack_sk(&skpv, sk);

  polyvec_ntt(&b);
  polyvec_basemul_acc_montgomery(&mp, &skpv, &b);
  poly_invntt_tomont(&mp);

  poly_sub(&mp, &v, &mp);
  poly_reduce(&mp);

  poly_tomsg(m, &mp);
}

// Compute the l2-norm squared of a polynomial
static double poly_l2_norm_sq(const poly *p)
{
  double norm_sq = 0.0;
  for (int i = 0; i < KYBER_N; i++)
  {
    norm_sq += (double)(p->coeffs[i] * p->coeffs[i]);
  }
  return norm_sq;
}

// Compute the l2-norm squared of a matrix of polynomials
static double poly_matrix_l2_norm_sq(poly **matrix, int rows, int cols)
{
  double norm_sq = 0.0;
  for (int i = 0; i < rows; i++)
  {
    for (int j = 0; j < cols; j++)
    {
      norm_sq += poly_l2_norm_sq(&matrix[i][j]);
    }
  }
  return norm_sq;
}

// Compute the infinity norm of a polynomial
static int16_t poly_inf_norm(const poly *p)
{
  int16_t max = 0;
  for (int i = 0; i < KYBER_N; i++)
  {
    int16_t abs_coeff = abs(p->coeffs[i]);
    if (abs_coeff > max)
    {
      max = abs_coeff;
    }
  }
  return max;
}

// Compute the infinity norm of a matrix of polynomials
static int16_t poly_matrix_inf_norm(poly **matrix, int rows, int cols)
{
  int16_t max = 0;
  for (int i = 0; i < rows; i++)
  {
    for (int j = 0; j < cols; j++)
    {
      int16_t norm = poly_inf_norm(&matrix[i][j]);
      if (norm > max)
      {
        max = norm;
      }
    }
  }
  return max;
}

static inline void poly_zeroize(poly *p)
{
  memset(p->coeffs, 0, sizeof(p->coeffs));
}

static inline void poly_scale(poly *r, uint16_t c)
{
  for (unsigned j = 0; j < KYBER_N; j++)
  {
    r->coeffs[j] = (int16_t)(r->coeffs[j] * c);
  }
}

void indcpa_mat_keypair(poly A[MAT_m][MAT_m],
                        poly B[MAT_m][MAT_l],
                        poly S[MAT_m][MAT_l],
                        const uint8_t coins[MAT_COINS])
{
  unsigned int i, j, k;
  uint8_t buf[2 * MAT_COINS];
  const uint8_t *pubseed = buf;
  const uint8_t *noiseseed = buf + MAT_COINS;
  uint8_t nonce = 0;

  poly t;
  poly_zeroize(&t);

  memcpy(buf, coins, MAT_COINS);
  buf[MAT_COINS] = MAT_COLS;
  hash_g(buf, buf, MAT_COINS + 1);

  // Compute mxm A in NTT domain
  gen_matrix_c(A, pubseed, 0);

  //  sample S
  for (i = 0; i < MAT_m; i++)
  {
    for (j = 0; j < MAT_l; j++)
    {
      poly_getnoise_eta1(&S[i][j], noiseseed, nonce++);
      poly_ntt(&S[i][j]);
    }
  }

  // TODO: Check norm condition
  // Compute AS + qE
  for (i = 0; i < MAT_m; i++)
  {
    for (k = 0; k < MAT_l; k++)
    {
      poly_zeroize(&B[i][k]);
      for (j = 0; j < MAT_m; j++)
      {
        poly_basemul_montgomery(&t, &A[i][j], &S[j][k]);
        poly_add(&B[i][k], &B[i][k], &t);
      }
      poly_invntt_tomont(&B[i][k]);
      poly_reduce(&B[i][k]);

      // Add noise qE_1
      poly_getnoise_eta1(&t, noiseseed, nonce++);
      poly_scale(&t, MAT_q);
      poly_add(&B[i][k], &B[i][k], &t);
      // poly_reduce(&B[i][k]);
    }
  }
}

void indcpa_mat_enc(poly U[MAT_m][MAT_k],
                    poly V[MAT_l][MAT_k],
                    poly A[MAT_m][MAT_m],
                    poly B[MAT_m][MAT_l],
                    poly M[MAT_l][MAT_k])
{
  unsigned int i, j, k;
  uint8_t buf[2 * MAT_COINS];
  const uint8_t *pubseed = buf;
  const uint8_t *noiseseed = buf + MAT_COINS;
  uint8_t nonce = 0;
  poly R[MAT_m][MAT_k], t;

  poly_zeroize(&t);

  // sample R mxk
  for (i = 0; i < MAT_m; i++)
  {
    for (j = 0; j < MAT_k; j++)
    {
      poly_getnoise_eta1(&R[i][j], noiseseed, nonce++);
      poly_ntt(&R[i][j]);
    }
  }

  for (i = 0; i < MAT_m; i++)
  {
    for (j = 0; j < MAT_l; j++)
    {
      poly_ntt(&B[i][j]);
    }
  }

  // Compute U = A^\trans R + qE_2
  for (i = 0; i < MAT_m; i++)
  {
    for (k = 0; k < MAT_k; k++)
    {
      poly_zeroize(&U[i][k]);
      poly_zeroize(&t);
      for (j = 0; j < MAT_m; j++)
      {
        poly_basemul_montgomery(&t, &A[j][i], &R[j][k]);
        poly_add(&U[i][k], &U[i][k], &t);
      }
      poly_invntt_tomont(&U[i][k]);
      poly_reduce(&U[i][k]);

      poly_getnoise_eta1(&t, noiseseed, nonce++);
      poly_scale(&t, MAT_q);
      poly_add(&U[i][k], &U[i][k], &t);
      // poly_reduce(&U[i][k]);
    }
  }

  // Compute V = B^\trans R +qE_3 + M
  for (i = 0; i < MAT_l; i++)
  {
    for (k = 0; k < MAT_k; k++)
    {
      poly_zeroize(&V[i][k]);
      poly_zeroize(&t);
      for (j = 0; j < MAT_m; j++)
      {
        poly_basemul_montgomery(&t, &B[j][i], &R[j][k]);
        poly_add(&V[i][k], &V[i][k], &t);
      }
      poly_invntt_tomont(&V[i][k]);
      poly_reduce(&V[i][k]);

      poly_getnoise_eta2(&t, noiseseed, nonce++);
      poly_scale(&t, MAT_q);
      poly_add(&V[i][k], &V[i][k], &t);
      // poly_reduce(&V[i][k]);

      // Add M
      poly_add(&V[i][k], &V[i][k], &M[i][k]);
      // poly_reduce(&V[i][k]);
    }
  }
}

void indcpa_mat_dec(poly M[MAT_l][MAT_k], poly S[MAT_m][MAT_l], poly U[MAT_m][MAT_k], poly V[MAT_l][MAT_k])
{
  unsigned int i, j, k;
  poly t, tmp;

  poly_zeroize(&t);
  poly_zeroize(&tmp);

  for (i = 0; i < MAT_m; i++)
  {
    for (j = 0; j < MAT_k; j++)
    {
      poly_ntt(&U[i][j]);
    }
  }

  // Compute S^\trans U
  for (i = 0; i < MAT_l; i++)
  {
    for (k = 0; k < MAT_k; k++)
    {
      poly_zeroize(&t);
      poly_zeroize(&M[i][k]);
      for (j = 0; j < MAT_m; j++)
      {
        poly_basemul_montgomery(&tmp, &S[j][i], &U[j][k]);
        poly_add(&t, &t, &tmp);
      }
      poly_invntt_tomont(&t);
      poly_reduce(&t);

      // Compute M[i][k] = V[i][k] - t
      poly_sub(&M[i][k], &V[i][k], &t);
      poly_reduce(&M[i][k]);
    }
  }
}

static inline int16_t fqadd(int16_t a, int16_t b)
{
  int32_t r = a + b;
  if (r >= MAT_Q)
    r -= MAT_Q;
  return (int16_t)r;
}

static inline int16_t fqsub(int16_t a, int16_t b)
{
  int32_t r = a - b;
  if (r < 0)
    r += MAT_Q;
  return (int16_t)r;
}

static inline int16_t fqmul(int16_t a, int16_t b)
{
  int32_t r = (int32_t)a * b;
  r %= MAT_Q;
  if (r < 0)
    r += MAT_Q;
  return (int16_t)r;
}

static int16_t fqpow(int16_t base, int32_t exp)
{
  int16_t result = 1;
  int16_t cur = base;
  while (exp > 0)
  {
    if (exp & 1)
      result = fqmul(result, cur);
    cur = fqmul(cur, cur);
    exp >>= 1;
  }
  return result;
}

// Compute modular inverse with fermats little theorem
static int16_t fqinv(int16_t a)
{
  return fqpow(a, MAT_Q - 2);
}

void lagrange_coeff(int16_t *lambda_i, int16_t user, const int16_t *users, int t)
{
  int16_t j, lambda = 1, numerator, denominator, denominator_inv, fraction;
  for (int x = 0; x < t; x++)
  {
    j = users[x];
    if (j == user)
      continue;
    numerator = fqsub(0, j);
    denominator = fqsub(user, j);
    denominator_inv = fqinv(denominator);
    fraction = fqmul(numerator, denominator_inv);
    lambda = fqmul(lambda, fraction);
  }
  *lambda_i = lambda;
}

// Threshold Deryption
static void indcpa_mat_tdec(poly M[MAT_l][MAT_k],
                            poly S[MAT_m][MAT_l],
                            poly U[MAT_m][MAT_k],
                            const int16_t user,
                            const int16_t *users,
                            const size_t t,
                            const uint8_t coins[MAT_COINS])
{
  int i, j, k;
  poly tpoly, tmp;

  uint8_t buf[2 * MAT_COINS];
  const uint8_t *pubseed = buf;
  const uint8_t *noiseseed = buf + MAT_COINS;
  uint8_t nonce = 0;
  int16_t lambda;

  lagrange_coeff(&lambda, user, users, t);

  poly_zeroize(&tpoly);
  poly_zeroize(&tmp);

  for (i = 0; i < MAT_m; i++)
  {
    for (j = 0; j < MAT_k; j++)
    {
      poly_ntt(&U[i][j]);
    }
  }

  // Compute S^\trans U
  for (i = 0; i < MAT_l; i++)
  {
    for (k = 0; k < MAT_k; k++)
    {
      poly_zeroize(&tpoly);
      poly_zeroize(&tmp);
      poly_zeroize(&M[i][k]);
      for (j = 0; j < MAT_m; j++)
      {
        poly_basemul_montgomery(&tmp, &S[j][i], &U[j][k]);
        poly_add(&tpoly, &tpoly, &tmp);
      }
      poly_invntt_tomont(&tpoly);
      poly_reduce(&tpoly);

      poly_scale(&tpoly, lambda);

      poly_zeroize(&tmp);
      poly_getnoise_eta1(&tmp, noiseseed, nonce++);
      poly_scale(&tmp, MAT_q);

      poly_add(&M[i][k], &tpoly, &tmp);
    }
  }
}

static void indcpa_mat_comb() {}

void poly_uniform(poly *a, const uint8_t seed[KYBER_SYMBYTES], uint16_t nonce)
{
  uint8_t extseed[KYBER_SYMBYTES + 2];
  uint8_t buf[SHAKE128_RATE];
  uint32_t ctr = 0;
  keccak_state state;

  shake128_init(&state);

  memcpy(extseed, seed, KYBER_SYMBYTES);
  extseed[KYBER_SYMBYTES] = (uint8_t)(nonce & 0xFF);
  extseed[KYBER_SYMBYTES + 1] = (uint8_t)((nonce >> 8) & 0xFF);
  shake128_absorb(&state, extseed, sizeof(extseed));
  shake128_finalize(&state);

  do
  {
    shake128_squeezeblocks(buf, 1, &state);
    ctr += rej_uniform(a->coeffs + ctr, KYBER_N - ctr, buf, SHAKE128_RATE);
  } while (ctr < KYBER_N);
}

static inline void poly_copy(poly *dst, const poly *src)
{
  for (int i = 0; i < KYBER_N; i++)
  {
    dst->coeffs[i] = src->coeffs[i];
  }
}

void gen_sharing_poly(poly P_coeff[][MAT_l][MAT_k],
                      poly S[MAT_l][MAT_k],
                      uint8_t master_seed[KYBER_SYMBYTES],
                      int t)
{
  // P(0) = S
  int i, j, k;
  for (i = 0; i < MAT_l; i++)
  {
    for (j = 0; j < MAT_k; j++)
    {
      poly_copy(&P_coeff[0][i][j], &S[i][j]);
    }
  }
  uint16_t nonce = 0;
  for (k = 1; k < t; k++)
  {
    for (i = 0; i < MAT_l; i++)
    {
      for (j = 0; j < MAT_k; j++)
      {
        poly_uniform(&P_coeff[k][i][j], master_seed, nonce++);
        poly_reduce(&P_coeff[k][i][j]);
      }
    }
  }
}

void compute_share(poly Share[MAT_l][MAT_k],
                   poly P_coeff[][MAT_l][MAT_k],
                   int u, int t)
{
  int i, j, k, r;
  int32_t tmp, val;
  poly tpoly;

  for (i = 0; i < MAT_l; i++)
  {
    for (j = 0; j < MAT_k; j++)
    {
      poly_copy(&Share[i][j], &P_coeff[0][i][j]);
      poly_reduce(&Share[i][j]);
    }
  }

  // Acc higher degree terms
  int16_t umod = (int16_t)(u % MAT_Q);
  if (umod < 0)
    umod += MAT_Q;
  int16_t pwr = 1;

  for (k = 1; k < t; k++)
  {
    tmp = (int32_t)pwr * umod;
    tmp %= MAT_Q;
    if (tmp < 0)
      tmp += MAT_Q;
    pwr = (int16_t)tmp;

    for (i = 0; i < MAT_l; i++)
    {
      for (j = 0; j < MAT_k; j++)
      {
        poly_copy(&tpoly, &P_coeff[k][i][j]);
        // scale tpoly by pwr
        for (r = 0; r < KYBER_N; r++)
        {
          val = tpoly.coeffs[r];
          val = val * pwr % MAT_Q;
          if (val < 0)
            val += MAT_Q;
          tpoly.coeffs[r] = (int16_t)val;
        }
        poly_add(&Share[i][j], &Share[i][j], &tpoly);
        poly_reduce(&Share[i][j]);
      }
    }
  }
}

void reconstruct_secret(poly S[MAT_l][MAT_k],
                        int16_t *U, int t,
                        poly Shares[][MAT_l][MAT_k])
{
  int i, j, k, r;
  int16_t u, lambda;
  int32_t val;
  poly tpoly;

  for (i = 0; i < MAT_l; i++)
  {
    for (j = 0; j < MAT_k; j++)
    {
      poly_zeroize(&S[i][j]);
    }
  }
  for (k = 0; k < t; k++)
  {
    u = U[k];
    lagrange_coeff(&lambda, u, U, t);
    for (i = 0; i < MAT_l; i++)
    {
      for (j = 0; j < MAT_k; j++)
      {
        poly_copy(&tpoly, &Shares[u - 1][i][j]);
        // Scale by lambda
        for (r = 0; r < KYBER_N; r++)
        {
          val = tpoly.coeffs[r];
          val = val * lambda % MAT_Q;
          if (val < 0)
            val += MAT_Q;
          tpoly.coeffs[r] = (int16_t)val;
        }
        poly_add(&S[i][j], &S[i][j], &tpoly);
      }
    }
  }
  for (i = 0; i < MAT_l; i++)
  {
    for (j = 0; j < MAT_k; j++)
    {
      poly_reduce(&S[i][j]);
    }
  }
}

// 1. ntt
// 2. SS
// 3. DKG

// Signature Scheme

// MAT_p ~ 49 bit
// Max number of signatures issues ~ 60 bit
// plaintext modulus MAT_q ~ 49 bits
// ciphertext modulus MAT_Q ~