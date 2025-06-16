#ifndef INDCPA_H
#define INDCPA_H

#include <stdint.h>
#include "params.h"
#include "polyvec.h"

#define gen_matrix KYBER_NAMESPACE(gen_matrix)
void gen_matrix(polyvec *a, const uint8_t seed[KYBER_SYMBYTES], int transposed);

#define indcpa_keypair_derand KYBER_NAMESPACE(indcpa_keypair_derand)
void indcpa_keypair_derand(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                           uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES],
                           const uint8_t coins[KYBER_SYMBYTES]);

#define indcpa_enc KYBER_NAMESPACE(indcpa_enc)
void indcpa_enc(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES]);

#define indcpa_dec KYBER_NAMESPACE(indcpa_dec)
void indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]);

#define indcpa_mat_keypair KYBER_NAMESPACE(indcpa_mat_keypair)
void indcpa_mat_keypair(poly A[MAT_m][MAT_m],
                        poly B[MAT_m][MAT_l],
                        poly S[MAT_m][MAT_l],
                        const uint8_t coins[MAT_COINS]);

#define indcpa_mat_enc KYBER_NAMESPACE(indcpa_mat_enc)
void indcpa_mat_enc(poly U[MAT_m][MAT_k], poly V[MAT_l][MAT_k], poly A[MAT_m][MAT_m], poly B[MAT_m][MAT_l], poly M[MAT_l][MAT_k]);

#define indcpa_mat_dec KYBER_NAMESPACE(indcpa_mat_dec)
void indcpa_mat_dec(poly M[MAT_l][MAT_k], poly S[MAT_m][MAT_l], poly U[MAT_m][MAT_k], poly V[MAT_l][MAT_k]);

#define poly_uniform KYBER_NAMESPACE(poly_uniform)
void poly_uniform(poly *a, const uint8_t seed[KYBER_SYMBYTES], uint16_t nonce);

#define lagrange_coeff KYBER_NAMESPACE(lagrange_coeff)
void lagrange_coeff(int16_t *lambda_i, int16_t user, const int16_t *users, int t);

#define gen_sharing_poly KYBER_NAMESPACE(gen_sharing_poly)
void gen_sharing_poly(poly P_coeff[][MAT_l][MAT_k],
                      poly S[MAT_l][MAT_k],
                      uint8_t master_seed[KYBER_SYMBYTES],
                      int t);

#define compute_share KYBER_NAMESPACE(compute_share)
void compute_share(poly Share[MAT_l][MAT_k],
                   poly P_coeff[][MAT_l][MAT_k],
                   int u, int t);

#define reconstruct_secret KYBER_NAMESPACE(reconstruct_secret)
void reconstruct_secret(poly S[MAT_l][MAT_k],
                        int16_t *U, int t,
                        poly Shares[][MAT_l][MAT_k]);

#endif
