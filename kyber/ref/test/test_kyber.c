#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "../kem.h"
#include "../randombytes.h"
#include "../indcpa.h"

// #define NTESTS 1000
#define NTESTS 1

static int test_keys(void)
{
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t ct[CRYPTO_CIPHERTEXTBYTES];
  uint8_t key_a[CRYPTO_BYTES];
  uint8_t key_b[CRYPTO_BYTES];

  // Alice generates a public key
  crypto_kem_keypair(pk, sk);

  // Bob derives a secret key and creates a response
  crypto_kem_enc(ct, key_b, pk);

  // Alice uses Bobs response to get her shared key
  crypto_kem_dec(key_a, ct, sk);

  if (memcmp(key_a, key_b, CRYPTO_BYTES))
  {
    printf("ERROR keys\n");
    return 1;
  }

  return 0;
}

static int test_invalid_sk_a(void)
{
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t ct[CRYPTO_CIPHERTEXTBYTES];
  uint8_t key_a[CRYPTO_BYTES];
  uint8_t key_b[CRYPTO_BYTES];

  // Alice generates a public key
  crypto_kem_keypair(pk, sk);

  // Bob derives a secret key and creates a response
  crypto_kem_enc(ct, key_b, pk);

  // Replace secret key with random values
  randombytes(sk, CRYPTO_SECRETKEYBYTES);

  // Alice uses Bobs response to get her shared key
  crypto_kem_dec(key_a, ct, sk);

  if (!memcmp(key_a, key_b, CRYPTO_BYTES))
  {
    printf("ERROR invalid sk\n");
    return 1;
  }

  return 0;
}

static int test_invalid_ciphertext(void)
{
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t ct[CRYPTO_CIPHERTEXTBYTES];
  uint8_t key_a[CRYPTO_BYTES];
  uint8_t key_b[CRYPTO_BYTES];
  uint8_t b;
  size_t pos;

  do
  {
    randombytes(&b, sizeof(uint8_t));
  } while (!b);
  randombytes((uint8_t *)&pos, sizeof(size_t));

  // Alice generates a public key
  crypto_kem_keypair(pk, sk);

  // Bob derives a secret key and creates a response
  crypto_kem_enc(ct, key_b, pk);

  // Change some byte in the ciphertext (i.e., encapsulated key)
  ct[pos % CRYPTO_CIPHERTEXTBYTES] ^= b;

  // Alice uses Bobs response to get her shared key
  crypto_kem_dec(key_a, ct, sk);

  if (!memcmp(key_a, key_b, CRYPTO_BYTES))
  {
    printf("ERROR invalid ciphertext\n");
    return 1;
  }

  return 0;
}

static void print_matrix_mxm(poly A[MAT_m][MAT_m], unsigned int coeff_count, char var)
{
  unsigned int i, j, k;
  for (i = 0; i < MAT_m; i++)
  {
    for (j = 0; j < MAT_m; j++)
    {
      printf("%c%d = [", var, i * 2 + j + 1);
      for (k = 0; k < coeff_count; k++)
      {
        printf("%d", A[i][j].coeffs[k]); // Adjust %g if coefficients are integers
        if (k < coeff_count - 1)
          printf(", ");
      }
      printf("]\n\n");
    }
  }
}

static void print_matrix_mxl(poly A[MAT_m][MAT_l], unsigned int coeff_count, char var)
{
  unsigned int i, j, k;
  for (i = 0; i < MAT_m; i++)
  {
    for (j = 0; j < MAT_l; j++)
    {
      printf("%c%d = [", var, i * 2 + j + 1);
      for (k = 0; k < coeff_count; k++)
      {
        printf("%d", A[i][j].coeffs[k]); // Adjust %g if coefficients are integers
        if (k < coeff_count - 1)
          printf(", ");
      }
      printf("]\n\n");
    }
  }
}

static void print_matrix_lxk(poly A[MAT_l][MAT_k], unsigned int coeff_count, char var)
{
  unsigned int i, j, k;
  for (i = 0; i < MAT_m; i++)
  {
    for (j = 0; j < MAT_l; j++)
    {
      printf("%c%d = [", var, i * 2 + j + 1);
      for (k = 0; k < coeff_count; k++)
      {
        printf("%d", A[i][j].coeffs[k]); // Adjust %g if coefficients are integers
        if (k < coeff_count - 1)
          printf(", ");
      }
      printf("]\n\n");
    }
  }
}

static void generate_random_01_poly(poly *p)
{
  for (int i = 0; i < KYBER_N; i++)
  {
    uint8_t random_byte;
    randombytes(&random_byte, 1);
    p->coeffs[i] = random_byte & 1;
  }
}

static void generate_binary_lxk(poly A[MAT_l][MAT_k])
{
  unsigned int i, j;
  for (i = 0; i < MAT_l; i++)
  {
    for (j = 0; j < MAT_k; j++)
    {
      generate_random_01_poly(&A[i][j]);
    }
  }
}
static int matcmp(poly M[MAT_l][MAT_k], poly M2[MAT_l][MAT_k])
{
  for (size_t i = 0; i < MAT_l; i++)
  {
    for (size_t j = 0; j < MAT_k; j++)
    {
      // Compare the coefficient arrays of two polynomials
      if (memcmp(M[i][j].coeffs, M2[i][j].coeffs, sizeof(uint16_t) * KYBER_N) != 0)
      {
        return 0; // Matrices are not identical
      }
    }
  }
  return 1; // Matrices are identical
}

static void print_message(const uint8_t msg[KYBER_SYMBYTES])
{
  for (unsigned int i = 0; i < KYBER_SYMBYTES; i++)
  {
    printf("%02x", msg[i]);
  }
  printf("\n");
}

// Caroline's test function
static int test_indcpa(void)
{
  poly A[MAT_m][MAT_m], B[MAT_m][MAT_l], S[MAT_m][MAT_l];
  poly U[MAT_m][MAT_k], V[MAT_l][MAT_k], M_decrypted[MAT_l][MAT_k];
  uint8_t coins[MAT_COINS];
  poly M[MAT_l][MAT_k];

  // Generate random bytes for coins
  randombytes(coins, MAT_COINS);

  // Generate keypair
  indcpa_mat_keypair(A, B, S, coins);

  // Generate random messages
  generate_binary_lxk(M);

  // Encrypt
  indcpa_mat_enc(U, V, A, B, M);

  // Decrypt
  indcpa_mat_dec(M_decrypted, S, U, V);

  if (matcmp(M_decrypted, M))
  {
    printf("EQUAL\n");
  }
  else
  {
    printf("NOT EQUAL\n");
  }

  return 0;
}

// Helper function to generate a random secret S
static void generate_random_secret(poly S[MAT_l][MAT_k])
{
  for (unsigned int i = 0; i < MAT_l; i++)
  {
    for (unsigned int j = 0; j < MAT_k; j++)
    {
      // Generate noise or uniform random polynomial as the secret
      uint8_t seed[KYBER_SYMBYTES];
      randombytes(seed, KYBER_SYMBYTES);
      poly_uniform(&S[i][j], seed, (uint16_t)(i * MAT_k + j));
      poly_reduce(&S[i][j]);
    }
  }
}

static int test_poly_uniform(void)
{
  // Known seed
  uint8_t seed[KYBER_SYMBYTES];
  memset(seed, 0x42, KYBER_SYMBYTES);

  uint16_t nonce = 123;

  poly a;
  poly_uniform(&a, seed, nonce);

  // Check that all coefficients are within [0, Q-1]
  for (int i = 0; i < KYBER_N; i++)
  {
    if (a.coeffs[i] < 0 || a.coeffs[i] >= KYBER_Q)
    {
      printf("Coefficient out of range: a.coeffs[%d] = %d\n", i, a.coeffs[i]);
      return 1;
    }
  }

  printf("Test poly_uniform: first few coeffs: ");
  for (int i = 0; i < 8; i++)
  {
    printf("%d ", a.coeffs[i]);
  }
  printf("\n");

  return 0;
}

static int test_sharing_and_reconstruction(void)
{
  // Parameters
  int t = 2;
  int n = 3;
  int16_t U[2] = {3, 2};

  poly S[MAT_l][MAT_k];
  poly S_rec[MAT_l][MAT_k];
  poly P_coeff[t][MAT_l][MAT_k];
  poly Shares[n][MAT_l][MAT_k];

  uint8_t master_seed[KYBER_SYMBYTES];
  randombytes(master_seed, KYBER_SYMBYTES);

  // Generate a random secret S
  generate_random_secret(S);

  // Generate sharing polynomial P(X)
  // P(0)=S and other coefficients random
  gen_sharing_poly(P_coeff, S, master_seed, t);

  // Compute shares for each user u in [1..n]
  for (int u = 1; u <= n; u++)
  {
    compute_share(Shares[u - 1], P_coeff, u, t);
  }

  // Reconstruct secret from subset U
  reconstruct_secret(S_rec, U, t, Shares);

  for (int i = 0; i < MAT_l; i++)
  {
    for (int j = 0; j < MAT_k; j++)
    {
      // Print first 8 coefficients for illustration
      printf("S[%d][%d].coeffs[0..7]: ", i, j);
      for (int r = 0; r < 8 && r < KYBER_N; r++)
      {
        printf("%d ", S[i][j].coeffs[r]);
      }
      printf("\n");

      printf("S_rec[%d][%d].coeffs[0..7]: ", i, j);
      for (int r = 0; r < 8 && r < KYBER_N; r++)
      {
        printf("%d ", S_rec[i][j].coeffs[r]);
      }
      printf("\n\n");
    }
  }

  // Check if S_rec == S
  if (matcmp(S, S_rec))
  {
    printf("test_sharing_and_reconstruction PASSED\n");
    return 0;
  }
  else
  {
    printf("test_sharing_and_reconstruction FAILED\n");
    return 1;
  }
}

int main(void)
{
  unsigned int i;
  int r;

  // for (i = 0; i < NTESTS; i++)
  // {
  //   r = test_keys();
  //   r |= test_invalid_sk_a();
  //   r |= test_invalid_ciphertext();
  //   if (r)
  //     return 1;
  // }

  for (i = 0; i < NTESTS; i++)
  {
    r = test_indcpa();
    r |= test_poly_uniform();
    r |= test_sharing_and_reconstruction();
    if (r)
      return 1;
  }

  // printf("CRYPTO_SECRETKEYBYTES:  %d\n", CRYPTO_SECRETKEYBYTES);
  // printf("CRYPTO_PUBLICKEYBYTES:  %d\n", CRYPTO_PUBLICKEYBYTES);
  // printf("CRYPTO_CIPHERTEXTBYTES: %d\n", CRYPTO_CIPHERTEXTBYTES);

  return 0;
}
