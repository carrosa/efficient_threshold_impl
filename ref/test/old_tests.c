
int test_enc_dec(void)
{
    poly A[M][M], B[M][L], S[M][L], U[M][K], V[L][K], Msg[L][K], Msg_dec[L][K];

    // Initialize all polys
    POLY_2D_INIT(A, M, M);
    POLY_2D_INIT(B, M, L);
    POLY_2D_INIT(S, M, L);
    POLY_2D_INIT(U, M, K);
    POLY_2D_INIT(V, L, K);
    POLY_2D_INIT(Msg, L, K);
    POLY_2D_INIT(Msg_dec, L, K);

    gen_message(Msg);

    gen_keys(A, B, S);
    encrypt(U, V, A, B, Msg);
    decrypt(Msg_dec, S, U, V);

    int ret;
    if (matcmp(&Msg[0][0], &Msg_dec[0][0], L, K))
    {
        printf("Enc/Dec SUCCESS\n");
        ret = 1;
    }
    else
    {
        printf("Enc/Dec FAILURE\n");
        ret = 0;
    }
    POLY_2D_CLEAR(A, M, M);
    POLY_2D_CLEAR(B, M, L);
    POLY_2D_CLEAR(S, M, L);
    POLY_2D_CLEAR(U, M, K);
    POLY_2D_CLEAR(V, L, K);
    POLY_2D_CLEAR(Msg, L, K);
    POLY_2D_CLEAR(Msg_dec, L, K);
    return ret;
}

static int test_shamir(void)
{
    // Parameters
    int t = 2;
    int n = 5;
    int32_t users[6] = {1, 2, 3, 4, 5};

    poly A[M][M], B[M][L], S[M][L], S_rec[M][L], P_coeff[2][M][L], Shares[5][M][L];

    // Initialize all polys
    POLY_2D_INIT(A, M, M);
    POLY_2D_INIT(B, M, L);
    POLY_2D_INIT(S, M, L);
    POLY_2D_INIT(S_rec, M, L);
    POLY_3D_INIT(P_coeff, t, M, L);
    POLY_3D_INIT(Shares, n, L, K);

    gen_keys(A, B, S);
    gen_sharing_poly(P_coeff, S, t);

    // Compute shares
    for (int u = 1; u <= n; u++)
    {
        compute_share(Shares[u - 1], P_coeff, u, t);
    }

    // Reconstruct secret from subset U
    reconstruct_secret(S_rec, users, t, Shares);

    // Check if S_rec == S
    int ret;
    if (matcmp(&S[0][0], &S_rec[0][0], M, L))
    {
        printf("test_sharing_and_reconstruction SUCCESS\n");
        ret = 1;
    }
    else
    {
        printf("test_sharing_and_reconstruction FAILURE\n");
        ret = 0;
    }

    POLY_2D_CLEAR(A, M, M);
    POLY_2D_CLEAR(B, M, L);
    POLY_2D_CLEAR(S, M, L);
    POLY_2D_CLEAR(S_rec, M, L);
    POLY_3D_CLEAR(P_coeff, t, L, K);
    POLY_3D_CLEAR(Shares, n, L, K);
    return ret;
}

static int test_tdec(void)
{
    int usize = 5;
    int threshold = 2;
    poly A[M][M],
        B[M][L],
        S[M][L],
        U[M][K],
        V[L][K],
        Msg[L][K],
        Msg_dec[L][K],
        ds_i[5][L][K],
        S_i[5][M][L],
        P_coeff[2][M][L];

    int32_t users[6] = {1, 2, 3, 4, 5};

    // Initialize all polys
    POLY_2D_INIT(A, M, M);
    POLY_2D_INIT(B, M, L);
    POLY_2D_INIT(S, M, L);
    POLY_2D_INIT(U, M, K);
    POLY_2D_INIT(V, L, K);
    POLY_2D_INIT(Msg, L, K);
    POLY_2D_INIT(Msg_dec, L, K);
    POLY_3D_INIT(ds_i, usize, L, K);
    POLY_3D_INIT(S_i, usize, M, L);
    POLY_3D_INIT(P_coeff, threshold, M, L);

    gen_keys(A, B, S);

    gen_sharing_poly(P_coeff, S, threshold);

    // Compute shares for each user u in [1..n]
    for (int u = 0; u < usize; u++)
    {
        compute_share(S_i[u], P_coeff, users[u], threshold);
    }

    gen_message(Msg);
    encrypt(U, V, A, B, Msg);

    // Convert U to ntt domain
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < K; j++)
        {
            poly_ntt(&U[i][j]);
            poly_reduce(&U[i][j]);
        }
    }

    for (int u = 0; u < threshold; u++)
    {
        tdecrypt(ds_i[u], S_i[u], U, users[u], users, threshold);
    }
    combine(Msg_dec, V, ds_i, threshold);

    int ret;
    if (matcmp(&Msg[0][0], &Msg_dec[0][0], L, K))
    {
        printf("TDec SUCCESS\n");
        ret = 1;
    }
    else
    {
        printf("TDec FAILURE\n");
        ret = 0;
    }

    POLY_2D_CLEAR(A, M, M);
    POLY_2D_CLEAR(B, M, L);
    POLY_2D_CLEAR(S, M, L);
    POLY_2D_CLEAR(U, M, K);
    POLY_2D_CLEAR(V, L, K);
    POLY_2D_CLEAR(Msg, L, K);
    POLY_2D_CLEAR(Msg_dec, L, K);
    POLY_3D_CLEAR(ds_i, usize, L, K);
    POLY_3D_CLEAR(S_i, usize, M, L);
    POLY_3D_CLEAR(P_coeff, threshold, M, L);

    return ret;
}

static int test_tdec_not_enough_keys(void)
{
    int usize = 5;
    int threshold = 2;
    int32_t users[6] = {1, 2, 3, 4, 5};
    poly A[M][M],
        B[M][L],
        S[M][L],
        U[M][K],
        V[L][K],
        Msg[L][K],
        Msg_dec[L][K],
        ds_i[5][L][K],
        S_i[5][M][L],
        P_coeff[2][M][L];

    // Initialize all polys
    POLY_2D_INIT(A, M, M);
    POLY_2D_INIT(B, M, L);
    POLY_2D_INIT(S, M, L);
    POLY_2D_INIT(U, M, K);
    POLY_2D_INIT(V, L, K);
    POLY_2D_INIT(Msg, L, K);
    POLY_2D_INIT(Msg_dec, L, K);
    POLY_3D_INIT(ds_i, usize, L, K);
    POLY_3D_INIT(S_i, usize, M, L);
    POLY_3D_INIT(P_coeff, threshold, M, L);

    gen_keys(A, B, S);

    gen_sharing_poly(P_coeff, S, threshold);

    // Compute shares for each user u in [1..n]
    for (int u = 0; u < usize - 1; u++)
    {
        compute_share(S_i[u], P_coeff, users[u], threshold);
    }

    gen_message(Msg);
    encrypt(U, V, A, B, Msg);

    // Convert U to ntt domain
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < K; j++)
        {
            poly_ntt(&U[i][j]);
            poly_reduce(&U[i][j]);
        }
    }

    for (int u = 0; u < threshold - 1; u++)
    {
        tdecrypt(ds_i[u], S_i[u], U, users[u], users, threshold);
    }
    combine(Msg_dec, V, ds_i, threshold);

    int ret;
    if (!matcmp(&Msg[0][0], &Msg_dec[0][0], L, K))
    {
        printf("TDec with not enough users SUCCESS (wrong decryption)\n");
        ret = 1;
    }
    else
    {
        printf("TDec with not enough users FAILURE (correct decryption)\n");
        ret = 0;
    }

    POLY_2D_CLEAR(A, M, M);
    POLY_2D_CLEAR(B, M, L);
    POLY_2D_CLEAR(S, M, L);
    POLY_2D_CLEAR(U, M, K);
    POLY_2D_CLEAR(V, L, K);
    POLY_2D_CLEAR(Msg, L, K);
    POLY_2D_CLEAR(Msg_dec, L, K);
    POLY_3D_CLEAR(ds_i, usize, L, K);
    POLY_3D_CLEAR(S_i, usize, M, L);
    POLY_3D_CLEAR(P_coeff, threshold, M, L);

    return ret;
}

static int test_sign_verify(void)
{
    poly s[L];
    poly A[K][L];
    poly yprime[K];
    poly c;
    poly z[L];
    poly h[K];
    uint8_t seed[SEEDBYTES];
    uint8_t mu[32]; // 32 bytes as an arbitrary message length
    int ret = 0;

    // Initialize polynomials
    POLY_2D_INIT(A, K, L);
    poly_1d_init(s, L);
    poly_1d_init(yprime, K);
    poly_init(&c);
    poly_1d_init(z, L);
    poly_1d_init(h, K);

    // Generate random seed and message
    randombytes(seed, SEEDBYTES);
    randombytes(mu, 32);

    // // Generate keys
    sigkeygen(A, yprime, s);

    // // Sign the message
    sign1p(s, A, yprime, mu, 32, c, z, h);

    // Verify the signature
    if (verify(A, yprime, c, z, h, mu, 32) == 1)
    {
        printf("Sign/Verify SUCCESS\n");
        ret = 1;
    }
    else
    {
        printf("Sign/Verify FAILURE\n");
        ret = 0;
    }

    // Clear polynomials
    poly_1d_clear(s, L);
    POLY_2D_CLEAR(A, K, L);
    poly_1d_clear(yprime, K);
    poly_clear(&c);
    poly_1d_clear(z, L);
    poly_1d_clear(h, K);

    return ret;
}

static int test_sign_verify_tampered_signature(void)
{
    poly s[L];
    poly A[K][L];
    poly yprime[K];
    poly c;
    poly z[L];
    poly h[K];
    uint8_t seed[SEEDBYTES];
    uint8_t mu[32];
    int ret;

    // Initialize polynomials
    poly_1d_init(s, L);
    POLY_2D_INIT(A, K, L);
    poly_1d_init(yprime, K);
    poly_init(&c);
    poly_1d_init(z, L);
    poly_1d_init(h, K);

    // Generate random seed and message
    randombytes(seed, SEEDBYTES);
    randombytes(mu, 32);

    // Generate keys
    sigkeygen(A, yprime, s);

    // Sign the message
    sign1p(s, A, yprime, mu, 32, c, z, h);

    // Tamper with the signature (e.g., increment h[0].coeffs[0])
    mpz_add_ui(h[0].coeffs[0], h[0].coeffs[0], 1);
    mpz_mod(h[0].coeffs[0], h[0].coeffs[0], GMP_q);

    // Verify the tampered signature
    if (verify(A, yprime, c, z, h, mu, 32) == 0)
    {
        printf("Sign/Verify with tampered signature SUCCESS (verification failed as expected)\n");
        ret = 1;
    }
    else
    {
        printf("Sign/Verify with tampered signature FAILURE (verification succeeded unexpectedly)\n");
        ret = 0;
    }

    // Clear polynomials
    poly_1d_clear(s, L);
    POLY_2D_CLEAR(A, K, L);
    poly_1d_clear(yprime, K);
    poly_clear(&c);
    poly_1d_clear(z, L);
    poly_1d_clear(h, K);

    return ret;
}

static int test_sign_verify_tampered_message(void)
{
    poly s[L];
    poly A[K][L];
    poly yprime[K];
    poly c;
    poly z[L];
    poly h[K];
    uint8_t seed[SEEDBYTES];
    uint8_t mu[32];
    int ret;

    // Initialize polynomials
    poly_1d_init(s, L);
    POLY_2D_INIT(A, K, L);
    poly_1d_init(yprime, K);
    poly_init(&c);
    poly_1d_init(z, L);
    poly_1d_init(h, K);

    // Generate random seed and message
    randombytes(seed, SEEDBYTES);
    randombytes(mu, 32);

    // Generate keys
    sigkeygen(A, yprime, s);

    // Sign the message
    sign1p(s, A, yprime, mu, 32, c, z, h);

    // Tamper with the message (e.g., flip a bit)
    mu[0] ^= 1;

    // Verify with tampered message
    if (verify(A, yprime, c, z, h, mu, 32) == 0)
    {
        printf("Sign/Verify with tampered message SUCCESS (verification failed as expected)\n");
        ret = 1;
    }
    else
    {
        printf("Sign/Verify with tampered message FAILURE (verification succeeded unexpectedly)\n");
        ret = 0;
    }

    // Clear polynomials
    poly_1d_clear(s, L);
    POLY_2D_CLEAR(A, K, L);
    poly_1d_clear(yprime, K);
    poly_clear(&c);
    poly_1d_clear(z, L);
    poly_1d_clear(h, K);

    return ret;
}

// static int test_as_scheme(void)
// {
//     printf("\n\n ----- \nRunning test_as_scheme(void) \n ----- \n\n");
//     uint64_t start, end;
//     int Bz = 100; // What should this be?
//     poly As[K][L];
//     pke_t *pke = malloc(sizeof(pke_t));
//     poly si[USERS][L];
//     poly yi[USERS][K];
//     poly yprime[K];
//     uint8_t h_yi[USERS][32];
//     ctx_t ctx_si[USERS];
//     // For signing need these as well
//     ctx_t ctx_s;
//     ctx_t ctx_z;
//     sig_t sig;
//     // pks_t pks;
//     pks_t *pks = malloc(sizeof(pks_t));
//     ctx_t ctx_r[THRESHOLD];
//     poly mu[M];
//     poly dsi[THRESHOLD][M];
//     poly wprime[K];
//     poly w_i[THRESHOLD][K];
//     uint8_t h_wi[THRESHOLD];

//     int user = 1;
//     int32_t users[USERS] = {1};
//     uint8_t h_Bi[USERS][32];

//     poly Si[USERS][LHAT][M];
//     poly Ei[USERS][KHAT][M];
//     poly Bi[USERS][KHAT][M];

//     poly(*Shares_Si)[USERS][LHAT][M] = malloc(sizeof(poly[USERS][USERS][LHAT][M]));
//     poly(*Shares_Ei)[USERS][KHAT][M] = malloc(sizeof(poly[USERS][USERS][KHAT][M]));
//     poly(*Bij_u)[USERS][LHAT][M] = malloc(sizeof(poly[USERS][USERS][LHAT][M]));

//     poly(*Bj)[USERS][KHAT][M] = malloc(sizeof(poly[USERS][USERS][KHAT][M]));
//     poly(*Sij)[USERS][LHAT][M] = malloc(sizeof(poly[USERS][USERS][KHAT][M]));
//     poly(*Be)[KHAT][M] = malloc(sizeof(poly[USERS][KHAT][M]));
//     poly(*Sei)[LHAT][M] = malloc(sizeof(poly[USERS][LHAT][M]));

//     POLY_3D_INIT(Si, USERS, LHAT, KHAT);
//     POLY_3D_INIT(Ei, USERS, KHAT, KHAT);
//     POLY_3D_INIT(Bi, USERS, KHAT, M);

//     poly_1d_init((poly *)Shares_Si, USERS * USERS * LHAT * M);
//     poly_1d_init((poly *)Shares_Ei, USERS * USERS * KHAT * M);
//     poly_1d_init((poly *)Bij_u, USERS * USERS * LHAT * M);
//     poly_1d_init((poly *)Bj, USERS * USERS * KHAT * M);
//     poly_1d_init((poly *)Sij, USERS * USERS * LHAT * M);
//     poly_1d_init((poly *)Be, KHAT * M);
//     poly_1d_init((poly *)Sei, LHAT * M);

//     // Initialize variables
//     POLY_2D_INIT(As, K, L);
//     POLY_2D_INIT(pke->A, KHAT, LHAT);
//     POLY_2D_INIT(pke->B, KHAT, M);
//     POLY_2D_INIT(pke->S, LHAT, M);
//     POLY_2D_INIT(si, USERS, L);
//     POLY_2D_INIT(yi, USERS, K);
//     poly_1d_init(yprime, K);
//     // Init sig
//     poly_init(&sig.c);
//     poly_1d_init(sig.z, L);
//     poly_1d_init(sig.h, K);
//     // Init pks
//     POLY_2D_INIT(pks->A, K, L);
//     poly_1d_init(pks->yprime, K);
//     // Init ctx_si
//     for (int i = 0; i < USERS; i++)
//     {
//         poly_1d_init(ctx_si[i].u, LHAT);
//         poly_1d_init(ctx_si[i].v, M);
//     }
//     // Init ctx_r
//     for (int i = 0; i < THRESHOLD; i++)
//     {
//         poly_1d_init(ctx_r[i].u, LHAT);
//         poly_1d_init(ctx_r[i].v, M);
//     }
//     // Init ctx_s
//     poly_1d_init(ctx_s.u, LHAT);
//     poly_1d_init(ctx_s.v, M);
//     // Init ski
//     poly_1d_init(mu, M);
//     // Init w_i & dsi
//     for (int i = 0; i < THRESHOLD; i++)
//     {
//         poly_1d_init(w_i[i], K);
//         poly_1d_init(dsi[i], M);
//     }
//     poly_1d_init(wprime, K);
//     // Init ctx_z
//     poly_1d_init(ctx_z.u, LHAT);
//     poly_1d_init(ctx_z.v, M);

//     int ret = 0;

//     // Setup encryption keys
//     start = get_cycles();
//     keygen_e_1d(pke->A, pke->B, pke->S);
//     end = get_cycles();
//     print_timing(start, end, "Keygen_e_1d");

//     // Sample mu
//     gen_message_1d(mu);

//     // Generate Ae
//     uint8_t seed[SEEDBYTES];
//     uint8_t nonce = 0;
//     gen_randomness(seed);
//     for (int i = 0; i < K; i++)
//     {
//         for (int j = 0; j < L; j++)
//         {
//             poly_uniform(&pks->A[i][j], seed, nonce++);
//             poly_ntt(&pks->A[i][j]);
//             poly_reduce(&pks->A[i][j]);
//         }
//     }

//     // Run dk_gen_1
//     for (int i = 0; i < USERS; i++)
//     {
//         start = get_cycles();
//         dk_gen_1(pke->A, Si[i], Ei[i], Bi[i], h_Bi[i]);
//         end = get_cycles();
//         print_timing(start, end, "dk_gen_1");
//     }
//     for (int i = 0; i < USERS; i++)
//     {
//         start = get_cycles();
//         dk_gen_2(
//             Si[i],
//             Ei[i],
//             Shares_Si[i],
//             Shares_Ei[i],
//             Bij_u[i],
//             pke,
//             users);
//         end = get_cycles();
//         print_timing(start, end, "dk_gen_2");
//     }
//     for (int i = 0; i < USERS; i++)
//     {
//         start = get_cycles();
//         dk_gen_3(
//             Bj,
//             Sij,
//             Be,
//             Sei);
//         end = get_cycles();
//         print_timing(start, end, "dk_gen_3");
//     }
//     printf("\n\n------ DKGEN SUCCESS ------\n\n");

//     // AS Keygen

//     // Run as_keygen_1
//     for (int i = 0; i < USERS; i++)
//     {
//         start = get_cycles();
//         as_keygen_1(pks->A, pke->A, pke->B, si[i], yi[i], h_yi[i], &ctx_si[i]);
//         end = get_cycles();
//         print_timing(start, end, "as_keygen_1");
//     }
//     // Same done for every user
//     start = get_cycles();
//     as_keygen_2(&ctx_s, yi, h_yi, ctx_si, pks->yprime);
//     end = get_cycles();
//     print_timing(start, end, "as_keygen_2");
//     printf("\n\n------ KEYGEN SUCCESS ------\n\n");

//     // Signing

//     for (int i = 0; i < USERS; i++)
//     {
//         start = get_cycles();
//         as_sign_1(pks->A, pke->A, pke->B, ctx_si[i]);
//         end = get_cycles();
//         print_timing(start, end, "as_sign_1");
//     }
//     printf("\n\n------ SIGN 1 ------\n\n");

//     for (int i = 0; i < THRESHOLD; i++)
//     {
//         start = get_cycles();
//         as_sign_2(&sig, *pks, ctx_r, ctx_s, ctx_z, pke->S, mu, dsi[i], w_i, h_wi, user, users);
//         end = get_cycles();
//     }
//     printf("\n\n------ SIGN 2 ------\n\n");
//     print_timing(start, end, "as_sign_2");

//     // Last sign, done same per user
//     start = get_cycles();
//     as_sign_3(&ctx_z, pks, &sig, pks->yprime, wprime, dsi, users);
//     end = get_cycles();
//     printf("\n\n------ SIGN 3 ------\n\n");
//     print_timing(start, end, "as_sign_3");

//     // print_sig(&sig);

//     // Verify
//     start = get_cycles();
//     as_verify(&sig, pks, mu, Bz);
//     printf("\n\n------ VERIFY ------\n\n");
//     end = get_cycles();
//     print_timing(start, end, "as_verify");

//     // Clear
//     POLY_3D_CLEAR(Si, USERS, LHAT, M);
//     POLY_3D_CLEAR(Ei, USERS, KHAT, M);
//     POLY_3D_CLEAR(Bi, USERS, KHAT, M);

//     poly_1d_clear((poly *)Shares_Si, USERS * USERS * LHAT * M);
//     poly_1d_clear((poly *)Shares_Ei, USERS * USERS * KHAT * M);
//     poly_1d_clear((poly *)Bij_u, USERS * USERS * LHAT * M);
//     poly_1d_clear((poly *)Bj, USERS * USERS * KHAT * M);
//     poly_1d_clear((poly *)Sij, USERS * USERS * LHAT * M);
//     poly_1d_clear((poly *)Be, KHAT * M);
//     poly_1d_clear((poly *)Sei, LHAT * M);

//     free(Shares_Si);
//     free(Shares_Ei);
//     free(Bij_u);
//     free(Bj);
//     free(Sij);
//     free(Be);
//     free(Sei);

//     POLY_2D_CLEAR(As, K, L);
//     POLY_2D_CLEAR(pke->A, KHAT, LHAT);
//     POLY_2D_CLEAR(pke->B, KHAT, M);
//     POLY_2D_CLEAR(pke->S, LHAT, M);
//     POLY_2D_CLEAR(si, USERS, L);
//     POLY_2D_CLEAR(yi, USERS, K);
//     poly_clear(&sig.c);
//     poly_1d_clear(sig.z, L);
//     poly_1d_clear(sig.h, K);
//     POLY_2D_CLEAR(pks->A, K, L);
//     poly_1d_clear(pks->yprime, K);
//     for (int i = 0; i < USERS; i++)
//     {
//         poly_1d_clear(ctx_si[i].u, LHAT);
//         poly_1d_clear(ctx_si[i].v, M);
//     }
//     for (int i = 0; i < THRESHOLD; i++)
//     {
//         poly_1d_clear(ctx_r[i].u, LHAT);
//         poly_1d_clear(ctx_r[i].v, M);
//     }
//     poly_1d_clear(ctx_s.u, LHAT);
//     poly_1d_clear(ctx_s.v, M);
//     poly_1d_clear(mu, M);
//     for (int i = 0; i < THRESHOLD; i++)
//     {
//         poly_1d_clear(w_i[i], K);
//         poly_1d_clear(dsi[i], M);
//     }
//     poly_1d_clear(wprime, K);
//     free(pke);
//     free(pks);
//     // poly_1d_clear(ctx_z.u, LHAT);
//     // poly_1d_clear(ctx_z.v, M);
//     // ret = 1;
//     // return ret;
// }



/// Keep for me

    // poly(*Shares_Si)[LHAT][M] = malloc(sizeof(poly[USERS][LHAT][M]));
    // poly_1d_init((poly *)Shares_Si, USERS * LHAT * M);

    // poly(*Shares_Ei)[KHAT][M] = malloc(sizeof(poly[USERS][KHAT][M]));
    // poly_1d_init((poly *)Shares_Ei, USERS * KHAT * M);

    // poly(*Bij_u)[LHAT][M] = malloc(sizeof(poly[USERS][LHAT][M]));
    // poly_1d_init((poly *)Bij_u, USERS * LHAT * M);

    // poly(*As)[L] = malloc(sizeof(poly[K][L]));
    // poly *yprime = malloc(sizeof(poly[K]));
    // ctx_t ctx_s, ctx_z;
    // ctx_t *ctx_r = malloc(sizeof(ctx_t[THRESHOLD]));

    // sig_t sig;
    // poly_init(&sig.c);
    // sig.z = malloc(sizeof(poly[L]));
    // sig.h = malloc(sizeof(poly[K]));

    // poly(*dsi)[M] = malloc(sizeof(poly[THRESHOLD][M]));
    // poly *wprime = malloc(sizeof(poly[K]));
    // poly(*w_i)[K] = malloc(sizeof(poly[THRESHOLD][K]));
    // uint8_t *h_wi = malloc(sizeof(uint8_t[THRESHOLD]));

    // POLY_2D_INIT(As, K, L);
    // POLY_2D_INIT(si, USERS, L);
    // POLY_2D_INIT(yi, USERS, K);
    // poly_1d_init(yprime, K);
    // poly_1d_init(sig.z, L);
    // poly_1d_init(sig.h, K);

    // for (int i = 0; i < THRESHOLD; i++)
    // {
    //     ctx_r[i].u = malloc(sizeof(poly[LHAT]));
    //     ctx_r[i].v = malloc(sizeof(poly[M]));
    //     poly_1d_init(ctx_r[i].u, LHAT);
    //     poly_1d_init(ctx_r[i].v, M);
    // }
    // ctx_s.u = malloc(sizeof(poly[LHAT]));
    // poly_1d_init(ctx_s.u, LHAT);
    // ctx_s.v = malloc(sizeof(poly[M]));
    // poly_1d_init(ctx_s.v, M);
    // ctx_z.u = malloc(sizeof(poly[LHAT]));
    // poly_1d_init(ctx_z.u, LHAT);
    // ctx_z.v = malloc(sizeof(poly[M]));
    // poly_1d_init(ctx_z.v, M);

    // for (int i = 0; i < THRESHOLD; i++)
    // {
    //     poly_1d_init(w_i[i], K);
    //     poly_1d_init(dsi[i], M);
    // }
    // poly_1d_init(wprime, K);

    // poly_1d_init((poly *)Bj, USERS * KHAT * M);
    // poly_1d_init((poly *)Sij, USERS * KHAT * M);
    // poly_1d_init((poly *)Be, KHAT * M);
    // poly_1d_init((poly *)Sei, LHAT * M);




    // Run as_keygen_1
    // for (int i = 0; i < USERS; i++)
    // {
    //     start = get_cycles();
    //     as_keygen_1(pks->A, pke->A, pke->B, si[i], yi[i], h_yi[i], &ctx_si[i]);
    //     end = get_cycles();
    //     print_timing(start, end, "as_keygen_1");
    // }
    // // Same done for every user
    // start = get_cycles();
    // as_keygen_2(&ctx_s, yi, h_yi, ctx_si, pks->yprime);
    // end = get_cycles();
    // print_timing(start, end, "as_keygen_2");
    // printf("\n\n------ KEYGEN SUCCESS ------\n\n");

    // // Signing

    // for (int i = 0; i < USERS; i++)
    // {
    //     start = get_cycles();
    //     as_sign_1(pks->A, pke->A, pke->B, ctx_si[i]);
    //     end = get_cycles();
    //     print_timing(start, end, "as_sign_1");
    // }
    // printf("\n\n------ SIGN 1 ------\n\n");

    // for (int i = 0; i < THRESHOLD; i++)
    // {
    //     start = get_cycles();
    //     as_sign_2(&sig, *pks, ctx_r, ctx_s, ctx_z, pke->S, mu, dsi[i], w_i, h_wi, user, users);
    //     end = get_cycles();
    // }
    // printf("\n\n------ SIGN 2 ------\n\n");
    // print_timing(start, end, "as_sign_2");

    // // Last sign, done same per user
    // start = get_cycles();
    // as_sign_3(&ctx_z, pks, &sig, pks->yprime, wprime, dsi, users);
    // end = get_cycles();
    // printf("\n\n------ SIGN 3 ------\n\n");
    // print_timing(start, end, "as_sign_3");

    // // print_sig(&sig);

    // // Verify
    // start = get_cycles();
    // as_verify(&sig, pks, mu, Bz);
    // printf("\n\n------ VERIFY ------\n\n");
    // end = get_cycles();
    // print_timing(start, end, "as_verify");


// OLD GHKSS256.c


// void decrypt_1d(poly Msg[L],
//                 poly S[L][L],
//                 poly U[L],
//                 poly V[L])
void decrypt_1d(poly Msg[M],
                poly S[LHAT][M],
                poly U[LHAT],
                poly V[M])
{
    unsigned int i, j, k;
    poly tmp;
    poly_init(&tmp);

    for (i = 0; i < LHAT; i++)
    {
        poly_ntt(&U[i]);
        poly_reduce(&U[i]);
    }

    // Compute V - S^\trans U
    for (i = 0; i < M; i++)
    {
        poly t;
        poly_init(&t);
        for (j = 0; j < LHAT; j++)
        {
            poly_pointwise_montgomery(&tmp, &S[j][i], &U[j]);
            poly_add(&t, &t, &tmp);
        }
        poly_reduce(&t);
        poly_invntt_tomont(&t);

        // Compute M[i][k] = V[i][k] - t
        poly_sub(&Msg[i], &V[i], &t);
        // poly_copy(&Msg[i], &V[i]);
        poly_reduce(&Msg[i]);
        poly_mod(&Msg[i], GMP_q);
        poly_clear(&t);
    }
    poly_clear(&tmp);
}

void gen_keys(
    poly A[M][M],
    poly B[M][L],
    poly S[M][L])
{
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    gen_randomness(seed);

    // Sample matrix A
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            poly_uniform(&A[i][j], seed, nonce++);
            poly_ntt(&A[i][j]);
            poly_reduce(&A[i][j]);
        }
    }

    // Sample S
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < L; j++)
        {
            gaussian_sampler(S[i][j].coeffs, N);
            poly_ntt(&S[i][j]);
            poly_reduce(&S[i][j]);
        }
    }

    for (int i = 0; i < M; i++)
    {
        for (int k = 0; k < L; k++)
        {
            poly t;
            poly_init(&t);

            for (int j = 0; j < M; j++)
            {
                poly_pointwise_montgomery(&t, &A[i][j], &S[j][k]);
                poly_add(&B[i][k], &B[i][k], &t);
            }
            poly_reduce(&B[i][k]);
            poly_invntt_tomont(&B[i][k]);

            // Add noise qE_1
            gaussian_sampler(t.coeffs, N);
            poly_scale(&t, GMP_q);
            poly_add(&B[i][k], &B[i][k], &t);
            poly_clear(&t);
        }
    }
}

void encrypt(poly U[M][K],
             poly V[L][K],
             poly A[M][M],
             poly B[M][L],
             poly Msg[L][K])
{
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    gen_randomness(seed);

    int i, j, k;
    poly R[M][K];
    // For temporary polynomials

    POLY_2D_INIT(R, M, K);

    // sample R mxk
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < K; j++)
        {
            gaussian_sampler(R[i][j].coeffs, N);
            poly_ntt(&R[i][j]);
        }
    }

    for (i = 0; i < M; i++)
    {
        for (j = 0; j < L; j++)
        {
            poly_ntt(&B[i][j]);
            poly_reduce(&B[i][j]);
        }
    }

    // Compute U = A^\trans R + qE_2
    for (i = 0; i < M; i++)
    {
        for (k = 0; k < K; k++)
        {
            poly t;
            poly_init(&t);
            for (j = 0; j < M; j++)
            {
                poly_pointwise_montgomery(&t, &A[j][i], &R[j][k]);
                poly_add(&U[i][k], &U[i][k], &t);
            }
            poly_reduce(&U[i][k]);
            poly_invntt_tomont(&U[i][k]);

            gaussian_sampler(t.coeffs, N);
            poly_scale(&t, GMP_q);
            poly_add(&U[i][k], &U[i][k], &t);
            poly_reduce(&U[i][k]);
            poly_clear(&t);
        }
    }

    // Compute V = B^\trans R +qE_3 + M
    for (i = 0; i < L; i++)
    {
        for (k = 0; k < K; k++)
        {
            poly t;
            poly_init(&t);
            for (j = 0; j < M; j++)
            {
                poly_pointwise_montgomery(&t, &B[j][i], &R[j][k]);
                poly_add(&V[i][k], &V[i][k], &t);
            }
            poly_reduce(&V[i][k]);
            poly_invntt_tomont(&V[i][k]);

            gaussian_sampler(t.coeffs, N);
            poly_scale(&t, GMP_q);
            poly_add(&V[i][k], &V[i][k], &t);
            poly_reduce(&V[i][k]);

            // Add M
            poly_add(&V[i][k], &V[i][k], &Msg[i][k]);
            poly_reduce(&V[i][k]);
            poly_clear(&t);
        }
    }
    POLY_2D_CLEAR(R, M, K);
}

void decrypt(poly Msg[L][K],
             poly S[M][L],
             poly U[M][K],
             poly V[L][K])
{
    unsigned int i, j, k;
    poly tmp;
    poly_init(&tmp);

    for (i = 0; i < M; i++)
    {
        for (j = 0; j < K; j++)
        {
            poly_ntt(&U[i][j]);
            poly_reduce(&U[i][j]);
        }
    }

    // Compute V - S^\trans U
    for (i = 0; i < L; i++)
    {
        for (k = 0; k < K; k++)
        {
            poly t;
            poly_init(&t);
            for (j = 0; j < M; j++)
            {
                poly_pointwise_montgomery(&tmp, &S[j][i], &U[j][k]);
                poly_add(&t, &t, &tmp);
            }
            poly_reduce(&t);
            poly_invntt_tomont(&t);

            // Compute M[i][k] = V[i][k] - t
            poly_sub(&Msg[i][k], &V[i][k], &t);
            // poly_copy(&Msg[i][k], &V[i][k]);
            poly_reduce(&Msg[i][k]);
            poly_mod(&Msg[i][k], GMP_q);
            poly_clear(&t);
        }
    }
    poly_clear(&tmp);
}

void tdecrypt(poly ds_i[L][K],
              poly S_i[M][L],
              poly U[M][K],
              int32_t user,
              int32_t *users,
              int t)
{
    uint8_t seed[SEEDBYTES];
    gen_randomness(seed);
    uint16_t nonce = 0;

    // Compute Lagrange coeff
    mpz_t lambda;
    mpz_init(lambda);
    lagrange_coeff(lambda, user, users, t);

    // Compute ds_i = \lambda_i S_i^\trans U + q E_i
    for (int l = 0; l < L; l++)
    {
        for (int k = 0; k < K; k++)
        {
            poly t1, t2;
            poly_inits(&t1, &t2, NULL);
            for (int m = 0; m < M; m++)
            {
                poly_pointwise_montgomery(&t1, &S_i[m][l], &U[m][k]);
                poly_add(&t2, &t2, &t1);
            }
            poly_reduce(&t2);
            poly_invntt_tomont(&t2);
            poly_scale(&t2, lambda); // scale includes reduce

            // Sample q \cdot E_ij
            gaussian_sampler(t1.coeffs, N);
            poly_scale(&t1, GMP_q);
            poly_add(&ds_i[l][k], &t1, &t2);
            poly_clears(&t1, &t2, NULL);
        }
    }
    mpz_clear(lambda);
}

// void combine_1d(poly Msg[L],
//                 poly V[L],
//                 poly (*ds)[L],
//                 int usize)

// void combine_1d(poly *Msg, // Should be L?
//                 poly V[M],
//                 poly (*ds)[M],
//                 int usize)
// {
//     int i, j, k;

//     // Sum decryption shares
//     for (i = 0; i < usize; i++)
//     {
//         for (j = 0; j < M; j++)
//         {
//             poly_add(&Msg[j], &Msg[j], &ds[i][j]);
//             poly_reduce(&Msg[j]);
//         }
//     }

//     // (V - \sum{ds}_{i\in\users}) \mod Q \mod q
//     for (i = 0; i < M; i++)
//     {
//         poly_sub(&Msg[i], &V[i], &Msg[i]);
//         poly_reduce(&Msg[i]);
//         poly_mod(&Msg[i], GMP_q);
//     }
// }

void combine(poly Msg[L][K],
             poly V[L][K],
             poly (*ds)[L][K],
             int usize)
{
    int i, j, k;

    // Sum decryption shares
    for (i = 0; i < usize; i++)
    {
        for (j = 0; j < L; j++)
        {
            for (k = 0; k < K; k++)
            {
                poly_add(&Msg[j][k], &Msg[j][k], &ds[i][j][k]);
                poly_reduce(&Msg[j][k]);
            }
        }
    }

    // (V - \sum{ds}_{i\in\users}) \mod Q \mod q
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < K; j++)
        {
            poly_sub(&Msg[i][j], &V[i][j], &Msg[i][j]);
            poly_reduce(&Msg[i][j]);
            poly_mod(&Msg[i][j], GMP_q);
        }
    }
}

void sigkeygen(
    poly A[K][L],
    poly yprime[L],
    poly s[L])
{
    poly e[K];

    poly_1d_init(e, K);

    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    gen_randomness(seed);

    for (int i = 0; i < L; i++)
    {
        gaussian_sampler(s[i].coeffs, N);
        poly_reduce(&s[i]);
        poly_ntt(&s[i]);
    }
    for (int i = 0; i < K; i++)
    {
        gaussian_sampler(e[i].coeffs, N);
        poly_reduce(&e[i]);
        poly_ntt(&e[i]);
        for (int j = 0; j < L; j++)
        {
            poly_uniform(&A[i][j], seed, nonce++);
            poly_reduce(&A[i][j]);
            poly_ntt(&A[i][j]);
        }
    }
    for (int i = 0; i < K; i++)
    {
        poly_pointwise_montgomery(&yprime[i], &A[0][i], &s[i]);
        for (int j = 1; j < L; j++)
        {
            poly_pointwise_montgomery(&s[i], &A[j][i], &s[i]);
            poly_add(&yprime[i], &yprime[i], &s[i]);
        }
        poly_invntt_tomont(&yprime[i]);
        poly_reduce(&yprime[i]);
    }

    for (int i = 0; i < K; i++)
    {
        poly_round(&yprime[i], &yprime[i], GMP_Q_Y, GMP_q);
    }

    poly_1d_clear(e, K);
}

void sign1p(
    poly s[L],
    poly A[K][L],
    poly yprime[K],
    uint8_t *mu,
    size_t mu_len,
    poly c,
    poly z[L],
    poly h[K])
{
    poly r[L];
    poly eprime[K];
    poly w[K];
    poly wprime[K];
    poly t[K];

    poly_1d_init(r, L);
    poly_1d_init(eprime, K);
    poly_1d_init(w, K);
    poly_1d_init(wprime, K);
    poly_1d_init(t, K);

    // Sample r
    for (int i = 0; i < L; i++)
    {
        gaussian_sampler(r[i].coeffs, N);
        poly_reduce(&r[i]);
        poly_ntt(&r[i]);
    }

    // w = Ar + eprime
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < L; j++)
        {
            poly tmp;
            poly_init(&tmp);
            poly_pointwise_montgomery(&tmp, &A[i][j], &r[j]);
            poly_add(&w[i], &w[i], &tmp);
            poly_clear(&tmp);
        }
        poly_reduce(&w[i]);
        poly_invntt_tomont(&w[i]);

        poly noise;
        poly_init(&noise);
        gaussian_sampler(noise.coeffs, N);
        poly_reduce(&noise);

        poly_add(&w[i], &w[i], &noise);
        poly_reduce(&w[i]);
        poly_clear(&noise);
    }

    for (int i = 0; i < K; i++)
    {
        poly_round(&wprime[i], &w[i], GMP_Q_W, GMP_q);
    }

    // Compute c = H(pk, wprime, msg)
    // compute_challenge(c, A, yprime, wprime, mu, mu_len);
    poly_ntt(&c);
    poly_reduce(&c);

    // Compute z = cs + r
    for (int i = 0; i < L; i++)
    {
        poly tmp;
        poly_init(&tmp);
        poly_pointwise_montgomery(&tmp, &c, &s[i]);
        poly_add(&z[i], &tmp, &r[i]);
        poly_reduce(&z[i]);
        poly_clear(&tmp);
    }

    // Compute t = round(A*z - 2^{k_y} * c * yprime) mod q_w
    poly Az[K];
    poly_1d_init(Az, K);
    // Convert z to ntt comain
    for (int i = 0; i < L; i++)
    {
        poly_ntt(&z[i]);
    }
    // Compute Az
    for (int i = 0; i < K; i++)
    {
        poly_pointwise_montgomery(&Az[i], &A[i][0], &z[0]);
        for (int j = 1; j < L; j++)
        {
            poly tmp;
            poly_init(&tmp);
            poly_pointwise_montgomery(&tmp, &A[i][j], &z[j]);
            poly_add(&Az[i], &Az[i], &tmp);
            poly_reduce(&Az[i]);
            poly_clear(&tmp);
        }
        poly_invntt_tomont(&Az[i]);
    }
    poly cy[K];
    poly_1d_init(cy, K);
    for (int i = 0; i < K; i++)
    {
        poly tmp;
        poly_init(&tmp);
        poly_ntt(&yprime[i]);
        poly_pointwise_montgomery(&tmp, &c, &yprime[i]);
        poly_invntt_tomont(&yprime[i]);
        poly_add(&cy[i], &tmp, &cy[i]);
        poly_reduce(&cy[i]);
        poly_clear(&tmp);
    }
    for (int i = 0; i < K; i++)
    {
        poly_sub(&t[i], &Az[i], &cy[i]);
        poly_reduce(&t[i]);
        poly_round(&t[i], &t[i], GMP_Q_W, GMP_q);
    }

    // Compute h = w' - t
    for (int i = 0; i < K; i++)
    {
        poly_sub(&h[i], &wprime[i], &t[i]);
        poly_reduce(&h[i]);
    }

    // Clean up
    poly_1d_clear(r, L);
    poly_1d_clear(eprime, K);
    poly_1d_clear(w, K);
    poly_1d_clear(wprime, K);
    poly_1d_clear(t, K);
    poly_1d_clear(Az, K);
    poly_1d_clear(cy, K);
}

int verify(
    poly A[K][L],
    poly yprime[K],
    poly c,
    poly z[L],
    poly h[K],
    uint8_t *mu,
    size_t mu_len)
{
    // Compute t = \lfloor A z - 2^{\kappa_y} c y' \rceil_{q_w}
    poly Az[K];
    poly_1d_init(Az, K);
    for (int i = 0; i < L; i++)
    {
        poly_ntt(&z[i]);
    }
    for (int i = 0; i < K; i++)
    {
        poly_pointwise_montgomery(&Az[i], &A[i][0], &z[0]);
        for (int j = 1; j < L; j++)
        {
            poly tmp;
            poly_init(&tmp);
            poly_pointwise_montgomery(&tmp, &A[i][j], &z[j]);
            poly_add(&Az[i], &Az[i], &tmp);
            poly_clear(&tmp);
        }
        poly_invntt_tomont(&Az[i]);
    }
    poly cy[K];
    poly_1d_init(cy, K);
    for (int i = 0; i < K; i++)
    {
        poly tmp;
        poly_init(&tmp);
        poly_pointwise_montgomery(&tmp, &c, &yprime[i]);
        poly_add(&cy[i], &tmp, &cy[i]);
        poly_scale(&cy[i], GMP_KAPPA_Y);
        poly_reduce(&cy[i]);
        poly_clear(&tmp);
    }
    poly t[K];
    poly_1d_init(t, K);
    for (int i = 0; i < K; i++)
    {
        poly_sub(&t[i], &Az[i], &cy[i]);
        poly_reduce(&t[i]);
        poly_round(&t[i], &t[i], GMP_Q_W, GMP_q);
    }

    // Compute w' = t + h
    poly wprime[K];
    poly_1d_init(wprime, K);
    for (int i = 0; i < K; i++)
    {
        poly_add(&wprime[i], &t[i], &h[i]);
        poly_reduce(&wprime[i]);
    }

    // Compute c' = H(pk, w', mu)
    poly c_prime;
    poly_init(&c_prime);
    // Made changes to below TODO; fix if wanted, not needed
    // compute_challenge(c_prime, A, yprime, wprime, mu, mu_len);

    // Check c == c'
    int c_equal = poly_equal(&c, &c_prime);

    // Clean up
    poly_1d_clear(Az, K);
    poly_1d_clear(cy, K);
    poly_1d_clear(t, K);
    poly_1d_clear(wprime, K);
    poly_clear(&c_prime);

    return c_equal;
}

void pp(const char *name, poly *p, int n)
{
    printf("%s = [", name);
    for (int i = 0; i < n; i++)
    {
        gmp_printf("%Zd", p->coeffs[i]);
        if (i != n - 1)
            printf(", ");
    }
    printf("]\n");
}


// void as_keygen_1(
//     poly As[K][L],
//     poly Ae[KHAT][LHAT],
//     poly Be[KHAT][M],
//     poly *si, // L
//     poly *yi, // K
//     uint8_t h_yi[32],
//     ctx_t *ctx)
// {
//     for (int i = 0; i < L; i++)
//     {
//         gaussian_sampler_y(si[i].coeffs, N);
//         poly_reduce(&si[i]);
//         poly_ntt(&si[i]);
//     }

//     poly ei[K];
//     poly_1d_init(ei, K);
//     for (int i = 0; i < K; i++)
//     {
//         gaussian_sampler_y(ei[i].coeffs, N);
//         poly_reduce(&ei[i]);
//         poly_ntt(&ei[i]);
//     }

//     for (int i = 0; i < K; i++)
//     {
//         for (int j = 0; j < L; j++)
//         {
//             poly tmp;
//             poly_init(&tmp);
//             poly_pointwise_montgomery(&tmp, &As[i][j], &si[j]);
//             poly_add(&yi[i], &yi[i], &tmp);
//             poly_clear(&tmp);
//         }
//         poly_invntt_tomont(&yi[i]);
//         poly_reduce(&yi[i]);
//         poly_add(&yi[i], &yi[i], &ei[i]);
//         poly_reduce(&yi[i]);
//     }

//     // Hash yi with H0
//     // poly_arr_hash_shake256(yi, K, h_yi, SHAKE256_RATE);
//     H0(yi, h_yi);

//     // Encryption of si
//     encrypt_1d(ctx->u, ctx->v, Ae, Be, si);

//     // // Print ctx.u
//     // pp("ctx.u", &ctx->u[0], N);
//     poly_1d_clear(ei, K);

//     // Compute NIZK proof of correct encryption of si
//     // Will not do in this implementation
//     // Do this function n times, then do keygen_final.
// }


// void as_keygen_2(
//     ctx_t *ctx_s,
//     poly yi[USERS][K],
//     uint8_t h_yi[USERS][32],
//     ctx_t ctx[USERS],
//     poly yprime[K])
// {
//     poly y[K];
//     poly_1d_init(y, K);
//     // Compute y
//     for (int i = 0; i < K; i++)
//     {
//         for (int j = 0; j < USERS; j++)
//         {
//             poly_add(&y[i], &y[i], &yi[j][i]);
//             poly_reduce(&y[i]);
//         }
//     }
//     // Check hashes:
//     for (int i = 0; i < USERS; i++)
//     {
//         uint8_t h_y[32];
//         // poly_arr_hash_shake256(&y, K, &h_y, SHAKE256_RATE);
//         H0(yi[i], h_y);
//         if (memcmp(h_y, h_yi[i], 32) != 0)
//         {
//             printf("Hashes do not match\n");
//             return;
//         }
//     }

//     // Compute homomorhpic addition of ctx_s
//     for (int i = 0; i < USERS; i++)
//     {
//         for (int j = 0; j < M; j++)
//         {
//             poly_add(&ctx_s->v[j], &ctx_s->v[j], &ctx[i].v[j]);
//             poly_reduce(&ctx_s->v[j]);
//         }
//         for (int j = 0; j < LHAT; j++)
//         {
//             poly_add(&ctx_s->u[j], &ctx_s->u[j], &ctx[i].u[j]);
//             poly_reduce(&ctx_s->u[j]);
//         }
//     }
//     mpz_t qprime;
//     mpz_init(qprime);
//     mpz_tdiv_q_2exp(qprime, GMP_q, KAPPA_Y);
//     for (int i = 0; i < K; i++)
//     {
//         for (int j = 0; j < N; j++)
//         {
//             mpz_tdiv_q_2exp(yprime[i].coeffs[j], y[i].coeffs[j], KAPPA_Y);
//             mpz_mod(yprime[i].coeffs[j], yprime[i].coeffs[j], qprime);
//         }
//     }
//     mpz_clear(qprime);
// }

// void as_sign_1(
//     poly As[K][L],
//     poly Ae[KHAT][LHAT],
//     poly Be[KHAT][M],
//     ctx_t ctx_si)
// {
//     poly r[M]; // r is padded with 0s from L to M
//     poly w[K];

//     poly_1d_init(r, M);
//     poly_1d_init(w, K);

//     // Sample r + NTT
//     for (int i = 0; i < L; i++)
//     {
//         gaussian_sampler_w(r[i].coeffs, N);
//         poly_reduce(&r[i]);
//         poly_ntt(&r[i]);
//     }

//     // Compute w = Ar + eprime
//     for (int i = 0; i < K; i++)
//     {
//         for (int j = 0; j < L; j++)
//         {
//             poly tmp;
//             poly_init(&tmp);
//             poly_pointwise_montgomery(&tmp, &As[i][j], &r[j]);
//             poly_add(&w[i], &w[i], &tmp);
//             poly_clear(&tmp);
//         }
//         poly_reduce(&w[i]);
//         poly_invntt_tomont(&w[i]);

//         poly noise; // Also known as eprime
//         poly_init(&noise);
//         gaussian_sampler_w(noise.coeffs, N);
//         poly_reduce(&noise);

//         poly_add(&w[i], &w[i], &noise);
//         poly_reduce(&w[i]);
//         poly_clear(&noise);
//     }

//     // Compute h_w = H(pk_e, w)
//     // TODO: Also hash pk_e (only w for now)
//     uint8_t h_w[32];
//     H1(w, h_w);
//     // poly_arr_hash_shake256(&w, K, &h_w, SHAKE256_RATE);

//     // // Compute encryption of randomness r
//     encrypt_1d(ctx_si.u, ctx_si.v, Ae, Be, r);
// }

// void as_sign_3(
//     ctx_t *ctx_z,
//     const pks_t *pks,
//     sig_t *sig,
//     const poly *yprime, // K
//     poly *wprime,       // K
//     poly dsi[THRESHOLD][M],
//     const int32_t *users) // USERS
// {
//     // Compute z = Comb(ctx_z, {ds_j}_{j \in U})
//     combine_1d(sig->z, ctx_z->v, dsi, THRESHOLD);
//     // Compute t = round(As * z - 2^{kappa_y} * c * yprime)_qw
//     poly t[K];
//     poly_1d_init(t, K);
//     // As * z
//     for (int i = 0; i < K; i++)
//     {
//         for (int j = 0; j < L; j++)
//         {
//             poly tmp;
//             poly_init(&tmp);

//             poly_pointwise_montgomery(&tmp, &pks->A[i][j], &sig->z[j]);
//             poly_add(&t[i], &t[i], &tmp);
//             poly_reduce(&t[i]);

//             poly_clear(&tmp);
//         }
//     }
//     // 2^{kappa_y} * c * yprime
//     for (int i = 0; i < K; i++)
//     {
//         poly tmp;
//         poly_init(&tmp);

//         poly_pointwise_montgomery(&tmp, &sig->c, &yprime[i]);
//         poly_scale(&tmp, GMP_KAPPA_Y);
//         poly_reduce(&tmp);
//         poly_sub(&t[i], &t[i], &tmp);
//         poly_reduce(&t[i]);

//         poly_clear(&tmp);
//     }
//     // Round t
//     mpz_t qprime;
//     mpz_init(qprime);
//     mpz_tdiv_q_2exp(qprime, GMP_q, KAPPA_W);
//     for (int i = 0; i < K; i++)
//     {
//         for (int j = 0; j < N; j++)
//         {
//             mpz_tdiv_q_2exp(t[i].coeffs[j], t[i].coeffs[j], KAPPA_W);
//             mpz_mod(t[i].coeffs[j], t[i].coeffs[j], qprime);
//         }
//     }
//     mpz_clear(qprime);

//     // Compute h = w' - t
//     for (int i = 0; i < K; i++)
//     {
//         poly_sub(&sig->h[i], &wprime[i], &t[i]);
//         poly_reduce(&sig->h[i]);
//     }

//     poly_1d_clear(t, K);
//     // sig = (c, z, h) [done, "returned" as sig]
//     // Clean up variables [TODO]
// }

// void dk_gen_3(
//     poly Bj[USERS][KHAT][M],
//     poly Sij[USERS][LHAT][M],
//     poly Be[KHAT][M],
//     poly Sei[LHAT][M])
// {
//     // Compute Be = sum_j (Bj)
//     for (int i = 0; i < USERS; i++)
//     {
//         for (int j = 0; j < KHAT; j++)
//         {
//             for (int k = 0; k < M; k++)
//             {
//                 poly_add(&Be[j][k], &Be[j][k], &Bj[i][j][k]);
//                 poly_reduce(&Be[j][k]);
//             }
//         }
//     }
//     // Compute Sei = sum_j (Sji)
//     for (int i = 0; i < USERS; i++)
//     {
//         for (int j = 0; j < LHAT; j++)
//         {
//             for (int k = 0; k < M; k++)
//             {
//                 poly_add(&Sei[j][k], &Sei[j][k], &Sij[i][j][k]);
//                 poly_reduce(&Sei[j][k]);
//             }
//         }
//     }
// }


int test_encryption(void)
{
    poly A[KHAT][LHAT],
        B[KHAT][M],
        S[LHAT][M],
        U[LHAT],
        V[M],
        Msg[M],
        Msg_dec[M];

    // Initialize all polys
    POLY_2D_INIT(A, KHAT, LHAT);
    POLY_2D_INIT(B, KHAT, M);
    POLY_2D_INIT(S, LHAT, M);
    poly_1d_init(U, LHAT);
    poly_1d_init(V, M);
    poly_1d_init(Msg, M);
    poly_1d_init(Msg_dec, M);

    gen_message_1d(Msg);

    keygen_e_1d(A, B, S);
    encrypt_1d(U, V, A, B, Msg);
    decrypt_1d(Msg_dec, S, U, V);

    int ret = 1;
    for (int i = 0; i < M; i++)
    {
        if (polycmp(&Msg[i], &Msg_dec[i]) != 1)
        {
            ret = 0;
            break;
        }
    }
    if (ret == 1)
    {
        printf("Enc/Dec SUCCESS\n");
    }
    else
    {
        printf("Enc/Dec FAILURE\n");
    }

    POLY_2D_CLEAR(A, KHAT, LHAT);
    POLY_2D_CLEAR(B, KHAT, M);
    POLY_2D_CLEAR(S, LHAT, M);
    poly_1d_clear(U, LHAT);
    poly_1d_clear(V, M);
    poly_1d_clear(Msg, M);
    poly_1d_clear(Msg_dec, M);
    return ret;
}
void pke_init(pke_t *pke)
{
    for (int i = 0; i < KHAT; i++)
    {
        for (int j = 0; j < LHAT; j++)
        {
            poly_init(&pke->A[i][j]);
        }
    }
    for (int i = 0; i < KHAT; i++)
    {
        for (int j = 0; j < M; j++)
        {
            poly_init(&pke->B[i][j]);
        }
    }
    for (int i = 0; i < LHAT; i++)
    {
        for (int j = 0; j < M; j++)
        {
            poly_init(&pke->S[i][j]);
        }
    }
}

void print_poly(const char *name, poly *p, int n)
{
    printf("%s = [", name);
    for (int i = 0; i < n; i++)
    {
        gmp_printf("%Zd", p->coeffs[i]);
        if (i != n - 1)
            printf(", ");
    }
    printf("]\n");
}

void print_sig(sig_t *sig)
{
    printf("Signature:\n");

    // Print c
    printf("c:\n");
    print_poly("c", &sig->c, N);

    // Print z
    for (int i = 0; i < L; i++)
    {
        char buf[32];
        snprintf(buf, sizeof(buf), "z[%d]", i);
        print_poly(buf, &sig->z[i], N);
    }

    // Print h
    for (int i = 0; i < K; i++)
    {
        char buf[32];
        snprintf(buf, sizeof(buf), "h[%d]", i);
        print_poly(buf, &sig->h[i], N);
    }
}

void print_timing(uint64_t start, uint64_t end, const char *label)
{
    uint64_t cycles = end - start;
    double time_in_seconds = (double)cycles / CPU_FREQ;
    double time_in_ms = time_in_seconds * 1000;
    printf("%s: %f seconds\n", label, time_in_seconds);
    printf("%s: %f ms\n", label, time_in_ms);
}

void print_avg_timing(uint64_t *timings, int count, const char *label)
{
    uint64_t total_cycles = 0;
    for (int i = 0; i < count; i++)
    {
        total_cycles += timings[i];
    }
    double avg_cycles = (double)total_cycles / count;
    double avg_time_in_seconds = (double)avg_cycles / CPU_FREQ;
    double avg_time_in_ms = avg_time_in_seconds * 1000;
    printf("Average %s: %f ms\n", label, avg_time_in_ms);
    printf("Average %s: %f seconds\n", label, avg_time_in_seconds);
}

int matcmp(poly *A, poly *B, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if (!polycmp(&A[i * cols + j], &B[i * cols + j]))
                return 0; // matrices differ
        }
    }
    return 1;
}

void print_two_vectors(poly A[L], poly B[L],
                       int coeff_count,
                       char labelA, char labelB)
{
    int i, j;
    for (i = 0; i < L; i++)
    {
        // Print a header for the current pair of polynomials:
        printf("%c[%d] and %c[%d]:\n", labelA, i, labelB, i);
        for (j = 0; j < coeff_count; j++)
        {
            printf("Coefficient %d: ", j);
            // Print the k-th coefficient from the polynomial in A
            mpz_out_str(stdout, 10, A[i].coeffs[j]);
            printf(" | ");
            // Print the j-th coefficient from the polynomial in B
            mpz_out_str(stdout, 10, B[i].coeffs[j]);
            printf("\n");
        }
        printf("\n"); // Separate different polynomial pairs
    }
}

// ignore warnings on the following function
void print_two_matrices(poly A[M][M], poly B[M][M],
                        int coeff_count,
                        char labelA, char labelB)
{
    int i, j, k;
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < M; j++)
        {
            // Print a header for the current pair of polynomials:
            printf("%c[%d][%d] and %c[%d][%d]:\n", labelA, i, j, labelB, i, j);
            for (k = 0; k < coeff_count; k++)
            {
                printf("Coefficient %d: ", k);
                // Print the k-th coefficient from the polynomial in A
                mpz_out_str(stdout, 10, A[i][j].coeffs[k]);
                printf(" | ");
                // Print the k-th coefficient from the polynomial in B
                mpz_out_str(stdout, 10, B[i][j].coeffs[k]);
                printf("\n");
            }
            printf("\n"); // Separate different polynomial pairs
        }
    }
}

void print_matrix(poly A[K][L], int coeff_count, char var)
{
    int i, j, k;
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < L; j++)
        {
            printf("%c[%d][%d] = [", var, i, j);
            for (k = 0; k < coeff_count; k++)
            {
                // Print each mpz_t coefficient in base 10.
                mpz_out_str(stdout, 10, A[i][j].coeffs[k]);
                if (k < coeff_count - 1)
                    printf(", ");
            }
            printf("]\n\n");
        }
    }
}
void encrypt_1d(poly U[LHAT],
                poly V[M],
                poly A[KHAT][LHAT],
                poly B[KHAT][M],
                poly Msg[M])
{
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    gen_randomness(seed);

    int i, j, k;
    poly R[KHAT];
    poly_1d_init(R, KHAT);

    // sample R mxk
    for (i = 0; i < KHAT; i++)
    {
        gaussian_sampler(R[i].coeffs, N);
        poly_ntt(&R[i]);
        poly_reduce(&R[i]);
    }

    for (i = 0; i < KHAT; i++)
    {
        for (j = 0; j < M; j++)
        {
            poly_ntt(&B[i][j]);
            poly_reduce(&B[i][j]);
        }
    }

    // Compute U = A^\trans R + qE_2
    for (i = 0; i < LHAT; i++)
    {
        for (j = 0; j < KHAT; j++)
        {
            poly t;
            poly_init(&t);
            poly_pointwise_montgomery(&t, &A[j][i], &R[j]);
            poly_add(&U[i], &U[i], &t);
            poly_clear(&t);
        }
        poly_reduce(&U[i]);
        poly_invntt_tomont(&U[i]);

        poly t;
        poly_init(&t);
        gaussian_sampler(t.coeffs, N);
        poly_scale(&t, GMP_q);
        poly_add(&U[i], &U[i], &t);
        poly_reduce(&U[i]);
    }

    // Compute V = B^\trans R +qE_3 + M
    for (i = 0; i < M; i++)
    {
        poly t;
        poly_init(&t);
        for (j = 0; j < KHAT; j++)
        {
            poly_pointwise_montgomery(&t, &B[j][i], &R[j]);
            poly_add(&V[i], &V[i], &t);
            poly_reduce(&V[i]);
        }
        poly_reduce(&V[i]);
        poly_invntt_tomont(&V[i]);

        gaussian_sampler(t.coeffs, N);
        poly_scale(&t, GMP_q);
        poly_add(&V[i], &V[i], &t);
        poly_reduce(&V[i]);

        // Add M
        poly_add(&V[i], &V[i], &Msg[i]);
        poly_reduce(&V[i]);
        poly_clear(&t);
    }
    poly_1d_clear(R, KHAT);
}
void gen_sharing_poly(poly P_coeff[][LHAT][M],
                      poly S[LHAT][M])
{
    uint8_t seed[SEEDBYTES];
    uint16_t nonce = 0;
    gen_randomness(seed);

    for (int i = 0; i < LHAT; i++)
    {
        for (int j = 0; j < M; j++)
        {
            poly_copy(&P_coeff[0][i][j], &S[i][j]);
        }
    }

    for (int k = 1; k < THRESHOLD; k++)
    {
        for (int i = 0; i < LHAT; i++)
        {
            for (int j = 0; j < M; j++)
            {
                poly_uniform(&P_coeff[k][i][j], seed, nonce++);
                poly_ntt(&P_coeff[k][i][j]);
                poly_reduce(&P_coeff[k][i][j]);
            }
        }
    }
}

void compute_share(poly Share[LHAT][M],
                   poly P_coeff[][LHAT][M],
                   int32_t user)
{
    mpz_t umod, tmp, pwr;
    mpz_inits(umod, tmp, pwr, NULL);
    mpz_set_ui(umod, user);
    mpz_set_ui(pwr, 1);
    // If wanted we can do umod % Q here, but also just define user in test with 0 < user < Q (which we are doing).

    // pwr = user^0 = 1 initially

    for (int r = 0; r < THRESHOLD; r++)
    {
        for (int i = 0; i < LHAT; i++)
        {
            for (int j = 0; j < M; j++)
            {
                poly tpoly;
                poly_init(&tpoly);
                poly_copy(&tpoly, &P_coeff[r][i][j]);
                poly_scale(&tpoly, pwr);
                poly_add(&Share[i][j], &Share[i][j], &tpoly);
                poly_reduce(&Share[i][j]);
                poly_clear(&tpoly);
            }
        }

        // Now multiply pwr by user for next iteration
        mpz_mul(tmp, umod, pwr);
        mpz_mod(tmp, tmp, GMP_Q);
        mpz_set(pwr, tmp);
    }
    mpz_clears(umod, tmp, NULL);
}

void reconstruct_secret(poly S[LHAT][M],
                        int32_t *users,
                        poly (*Shares)[LHAT][M])
{
    int32_t u;

    poly tpoly;
    mpz_t lambda;
    mpz_init(lambda);
    poly_init(&tpoly);

    for (int k = 0; k < THRESHOLD; k++)
    {
        u = users[k];
        lagrange_coeff(lambda, u, users, THRESHOLD);
        // lambda = to_montgomery(lambda);
        for (int i = 0; i < LHAT; i++)
        {
            for (int j = 0; j < M; j++)
            {
                poly_copy(&tpoly, &Shares[u - 1][i][j]);
                // Scale by lambda
                poly_scale(&tpoly, lambda);
                poly_add(&S[i][j], &S[i][j], &tpoly);
                poly_reduce(&S[i][j]);
            }
        }
    }
    mpz_clear(lambda);
    poly_clear(&tpoly);
}

void tdecrypt_1d(poly ds_i[M],
                 poly S_i[LHAT][M],
                 poly U[LHAT],
                 int32_t user,
                 int32_t *users,
                 int t)
{
    uint8_t seed[SEEDBYTES];
    gen_randomness(seed);
    uint16_t nonce = 0;

    // Compute Lagrange coeff
    mpz_t lambda;
    mpz_init(lambda);
    lagrange_coeff(lambda, user, users, t);

    // Compute ds_i = \lambda_i S_i^\trans U + q E_i
    for (int i = 0; i < M; i++)
    {
        poly t1, t2;
        poly_inits(&t1, &t2, NULL);
        for (int j = 0; j < LHAT; j++)
        {
            poly_pointwise_montgomery(&t1, &S_i[j][i], &U[j]);
            poly_add(&t2, &t2, &t1);
        }
        poly_reduce(&t2);
        poly_invntt_tomont(&t2);
        poly_scale(&t2, lambda); // scale includes reduce

        // Sample q \cdot E_ij
        gaussian_sampler(t1.coeffs, N);
        poly_scale(&t1, GMP_q);
        poly_add(&ds_i[i], &t1, &t2);
        poly_clears(&t1, &t2, NULL);
    }
    mpz_clear(lambda);
}
