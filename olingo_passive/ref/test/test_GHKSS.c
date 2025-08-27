#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "../params.h"
#include "../GHKSS.h"
#include "../poly.h"

void generate_random_binary_msg(poly Msg[L][K])
{
    uint8_t random_byte;
    int i, j, k;

    // Iterate over the Msg[L][K] matrix
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < K; j++)
        {
            // Iterate over coefficients of each polynomial
            for (k = 0; k < N; k++)
            {
                // Generate a random byte
                randombytes(&random_byte, 1);

                // Set coefficient to 0 or 1 based on the least significant bit
                Msg[i][j].coeffs[k] = random_byte & 1;
            }
        }
    }
}

static int matcmp(poly Msg1[L][K], poly Msg2[L][K])
{
    for (size_t i = 0; i < L; i++)
    {
        for (size_t j = 0; j < K; j++)
        {
            // Compare the coefficient arrays of two polynomials
            if (memcmp(Msg1[i][j].coeffs, Msg2[i][j].coeffs, sizeof(uint32_t) * N) != 0)
            {
                return 0; // Matrices are not identical
            }
        }
    }
    return 1; // Matrices are identical
}

static void print_matrix(poly A[M][M], unsigned int coeff_count, char var)
{
    unsigned int i, j, k;
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < M; j++)
        {
            printf("%c[%d][%d] = [", var, i, j);
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

int test_split_key(void)
{
    const size_t usize = 3;
    poly A[M][M], B[M][L], S[M][L], S_rec[M][L], S_i[usize][M][L], oS[L][K];

    indcpa_mat_keypair(A, B, S);
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < L; j++)
        {
            poly_copy(&oS[i][j], &S[i][j]);
        }
    }
    split_key_S(S_i, S, usize);

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < L; j++)
        {
            poly_zeroize(&S_rec[i][j]);
        }
    }
    for (size_t u = 0; u < usize; u++)
    {
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < L; j++)
            {
                poly_add(&S_rec[i][j], &S_rec[i][j], &S_i[u][i][j]);
            }
        }
    }
    if (matcmp(oS, S_rec))
    {
        printf("S split/reconstruction SUCCESS\n");
    }
    else
    {
        printf("S split/reconstruction FAILURE\n");
    }
}

static void generate_random_secret(poly S[L][K])
{
    int nonce = 0;
    uint8_t seed[SEEDBYTES];
    randombytes(seed, SEEDBYTES);
    for (unsigned int i = 0; i < L; i++)
    {
        for (unsigned int j = 0; j < K; j++)
        {
            // Generate noise or uniform random polynomial as the secret
            poly_uniform(&S[i][j], seed, nonce++);
            poly_reduce(&S[i][j]);
        }
    }
}

static int test_shamir(void)
{
    // Parameters
    int t = 2;
    int n = 5;
    int32_t users[2] = {5, 2};

    poly A[M][M], B[M][L], S[M][L], U[M][K], V[L][K], Msg[L][K], Msg_dec[L][K];
    poly S_rec[M][L];
    poly P_coeff[t][L][K];
    poly Shares[n][L][K];

    indcpa_mat_keypair(A, B, S);
    generate_random_binary_msg(Msg);
    indcpa_mat_enc(U, V, A, B, Msg);

    // Generate a random secret S

    // Generate sharing polynomial P(X)
    // P(0)=S and other coefficients random
    gen_shr_poly(P_coeff, S, t);

    // Compute shares for each user u in [1..n]
    for (int u = 1; u <= n; u++)
    {
        compute_shr(Shares[u - 1], P_coeff, u, t);
    }

    // Reconstruct secret from subset U
    reconstruct_secret(S_rec, users, t, Shares);

    indcpa_mat_dec(Msg_dec, S_rec, U, V);

    // Check if S_rec == S
    if (matcmp(Msg, Msg_dec))
    {
        printf("test_sharing_and_reconstruction SUCCESS\n");
        return 0;
    }
    else
    {
        printf("test_sharing_and_reconstruction FAILURE\n");
        return 1;
    }
}

int test_enc_dec(void)
{
    poly A[M][M], B[M][L], S[M][L], U[M][K], V[L][K], Msg[L][K], Msg_dec[L][K];

    indcpa_mat_keypair(A, B, S);
    generate_random_binary_msg(Msg);
    indcpa_mat_enc(U, V, A, B, Msg);
    indcpa_mat_dec(Msg_dec, S, U, V);

    if (matcmp(Msg, Msg_dec))
    {
        printf("Enc/Dec SUCCESS\n");
    }
    else
    {
        printf("Enc/Dec FAILURE\n");
    }
}

int test_tdec(void)
{
    const size_t usize = 2;
    const size_t threshold = 2;
    const int32_t users[2] = {1, 2};
    poly A[M][M],
        B[M][L],
        S[M][L],
        U[M][K],
        V[L][K],
        Msg[L][K],
        Msg_dec[L][K],
        ds_i[usize][L][K],
        S_i[usize][M][L],
        P_coeff[threshold][M][L];

    poly Share1[M][L], Share2[M][L];
    poly pDec1[L][K], pDec2[L][K];

    indcpa_mat_keypair(A, B, S);

    gen_shr_poly(P_coeff, S, threshold);
    // print_matrix(P_coeff, 2, 'P');

    // Compute shares for each user u in [1..n]
    for (int u = 0; u < usize; u++)
    {
        compute_shr(S_i[u], P_coeff, users[u], threshold);
    }
    // split_key_S(S_i, S, usize);

    // Compute shares
    // compute_shr(Share1, P_coeff, 1, threshold);
    // compute_shr(Share2, P_coeff, 2, threshold);

    generate_random_binary_msg(Msg);
    // print_matrix(Msg, N, 'N');
    indcpa_mat_enc(U, V, A, B, Msg);

    // Convert U to ntt domain
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < K; j++)
        {
            poly_ntt(&U[i][j]);
            poly_reduce(&U[i][j]);
        }
    }

    for (size_t u = 0; u < threshold; u++)
    {
        indcpa_mat_tdec(ds_i[u], S_i[u], U, V, users[u], users, threshold);
    }
    // Decryption
    // indcpa_mat_tdec(pDec1, Share1, U, V, 1, users, threshold);
    // indcpa_mat_tdec(pDec2, Share2, U, V, 2, users, threshold);
    // static poly ds[2][L][K];

    // // Copy pDec1 into ds[0], pDec2 into ds[1]
    // memcpy(ds[0], pDec1, sizeof(ds[0]));
    // memcpy(ds[1], pDec2, sizeof(ds[1]));
    indcpa_mat_comb(Msg_dec, V, ds_i, threshold);

    // print_matrix(Msg, 8, 'N');
    // print_matrix(Msg_dec, 8, 'M');

    if (matcmp(Msg, Msg_dec))
    {
        printf("TDec SUCCESS\n");
    }
    else
    {
        printf("TDec FAILURE\n");
    }
}

int main(void)
{
    int r;

    r = test_enc_dec();
    r |= test_shamir();
    r |= test_tdec();

    if (r)
        return 1;

    return 0;
}
