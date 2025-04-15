#ifndef GHKSS_PARAMS_H
#define GHKSS_PARAMS_H

#define GHKSS_NAMESPACE(s) GHKSS_ref_##s

#define MAT_COINS 32
#define MAT_ROWS 2
#define MAT_COLS 2
#define MAT_POLY_BYTES 384
#define MAT_POLYVEC_BYTES (MAT_ROWS * MAT_POLY_BYTES)
#define MAT_PUBLICKEY_BYTES (MAT_COLS * MAT_POLYVEC_BYTES)
#define MAT_SECRETKEY_BYTES (MAT_COLS * MAT_POLYVEC_BYTES)
#define MAT_Q 3329
#define MAT_q 2 // 3329
#define MAT_N 256

#define MAT_m 2
#define MAT_l 2
#define MAT_k 2

#define SEEDBYTES 32

#define M 4

#endif