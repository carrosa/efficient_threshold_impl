#include <math.h>
#include <stdlib.h>

#include <flint/flint.h>
#include <flint/fmpz_mod_poly.h>

#include "blake3.h"
#include "common.h"
#include "test.h"
#include "bench.h"
#include "assert.h"
#include "sample_z_small.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

#define MSGS 22 // Not using shuffle prover for this implementation.

static void lin_hash(params::poly_q &beta, comkey_t &key, commit_t x,
                     commit_t y[SIZE], params::poly_q alpha[SIZE], params::poly_q u[SIZE],
                     params::poly_q t, params::poly_q _t[SIZE])
{
    uint8_t hash[BLAKE3_OUT_LEN];
    blake3_hasher hasher;

    blake3_hasher_init(&hasher);

    /* Hash public key. */
    for (size_t i = 0; i < HEIGHT; i++)
    {
        for (int j = 0; j < WIDTH - HEIGHT; j++)
        {
            blake3_hasher_update(&hasher, (const uint8_t *)key.A1[i][j].data(),
                                 16 * DEGREE);
        }
    }
    for (size_t j = 0; j < WIDTH; j++)
    {
        blake3_hasher_update(&hasher, (const uint8_t *)key.A2[0][j].data(),
                             16 * DEGREE);
    }

    /* Hash alpha, beta from linear relation. */
    for (size_t i = 0; i < SIZE; i++)
    {
        blake3_hasher_update(&hasher, (const uint8_t *)alpha[i].data(),
                             16 * DEGREE);
    }

    blake3_hasher_update(&hasher, (const uint8_t *)x.c1.data(), 16 * DEGREE);
    for (int i = 0; i < SIZE; i++)
    {
        blake3_hasher_update(&hasher, (const uint8_t *)y[i].c1.data(), 16 * DEGREE);
    }
    for (size_t i = 0; i < x.c2.size(); i++)
    {
        blake3_hasher_update(&hasher, (const uint8_t *)x.c2[i].data(),
                             16 * DEGREE);
        for (int j = 0; j < SIZE; j++)
        {
            blake3_hasher_update(&hasher, (const uint8_t *)y[j].c2[i].data(),
                                 16 * DEGREE);
        }
    }

    for (int i = 0; i < SIZE; i++)
    {
        blake3_hasher_update(&hasher, (const uint8_t *)u[i].data(), 16 * DEGREE);
    }
    blake3_hasher_update(&hasher, (const uint8_t *)t.data(), 16 * DEGREE);
    for (int i = 0; i < SIZE; i++)
    {
        blake3_hasher_update(&hasher, (const uint8_t *)_t[i].data(), 16 * DEGREE);
    }

    blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);

    /* Sample challenge from RNG seeded with hash. */
    nfl::fastrandombytes_seed(hash);
    bdlop_sample_chal(beta);
    nfl::fastrandombytes_reseed();
}

static void lin_hash(params::poly_q &beta, comkey_t &key, commit_t x,
                     commit_t y, params::poly_q alpha[22], params::poly_q &u,
                     params::poly_q t, params::poly_q _t)
{
    uint8_t hash[BLAKE3_OUT_LEN];
    blake3_hasher hasher;

    blake3_hasher_init(&hasher);

    /* Hash public key. */
    for (size_t i = 0; i < HEIGHT; i++)
    {
        for (int j = 0; j < WIDTH - HEIGHT; j++)
        {
            blake3_hasher_update(&hasher, (const uint8_t *)key.A1[i][j].data(),
                                 16 * DEGREE);
        }
    }
    for (size_t j = 0; j < WIDTH; j++)
    {
        blake3_hasher_update(&hasher, (const uint8_t *)key.A2[0][j].data(),
                             16 * DEGREE);
    }

    /* Hash alpha, beta from linear relation. */
    for (size_t i = 0; i < SIZE; i++)
    {
        blake3_hasher_update(&hasher, (const uint8_t *)alpha[i].data(),
                             16 * DEGREE);
    }

    blake3_hasher_update(&hasher, (const uint8_t *)x.c1.data(), 16 * DEGREE);
    blake3_hasher_update(&hasher, (const uint8_t *)y.c1.data(), 16 * DEGREE);
    for (size_t i = 0; i < x.c2.size(); i++)
    {
        blake3_hasher_update(&hasher, (const uint8_t *)x.c2[i].data(),
                             16 * DEGREE);
        blake3_hasher_update(&hasher, (const uint8_t *)y.c2[i].data(),
                             16 * DEGREE);
    }

    blake3_hasher_update(&hasher, (const uint8_t *)u.data(), 16 * DEGREE);
    blake3_hasher_update(&hasher, (const uint8_t *)t.data(), 16 * DEGREE);
    blake3_hasher_update(&hasher, (const uint8_t *)_t.data(), 16 * DEGREE);

    blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);

    /* Sample challenge from RNG seeded with hash. */
    nfl::fastrandombytes_seed(hash);
    bdlop_sample_chal(beta);
    nfl::fastrandombytes_reseed();
}

static void poly_inverse(params::poly_q &inv, params::poly_q p)
{
    std::array<mpz_t, params::poly_q::degree> coeffs;
    fmpz_t q;
    fmpz_mod_poly_t poly, irred;
    fmpz_mod_ctx_t ctx_q;

    fmpz_init(q);
    for (size_t i = 0; i < params::poly_q::degree; i++)
    {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    fmpz_set_mpz(q, params::poly_q::moduli_product());
    fmpz_mod_ctx_init(ctx_q, q);
    fmpz_mod_poly_init(poly, ctx_q);
    fmpz_mod_poly_init(irred, ctx_q);

    p.poly2mpz(coeffs);
    fmpz_mod_poly_set_coeff_ui(irred, params::poly_q::degree, 1, ctx_q);
    fmpz_mod_poly_set_coeff_ui(irred, 0, 1, ctx_q);

    for (size_t i = 0; i < params::poly_q::degree; i++)
    {
        fmpz_mod_poly_set_coeff_mpz(poly, i, coeffs[i], ctx_q);
    }
    fmpz_mod_poly_invmod(poly, poly, irred, ctx_q);

    for (size_t i = 0; i < params::poly_q::degree; i++)
    {
        fmpz_mod_poly_get_coeff_mpz(coeffs[i], poly, i, ctx_q);
    }

    inv.mpz2poly(coeffs);

    fmpz_clear(q);
    for (size_t i = 0; i < params::poly_q::degree; i++)
    {
        mpz_clear(coeffs[i]);
    }
}

static void simul_inverse(params::poly_q inv[MSGS], params::poly_q m[MSGS])
{
    params::poly_q u, t[MSGS];
    inv[0] = m[0];
    t[0] = m[0];

    for (size_t i = 1; i < MSGS; i++)
    {
        t[i] = m[i];
        inv[i] = inv[i - 1] * m[i];
    }

    u = inv[MSGS - 1];
    u.invntt_pow_invphi();
    poly_inverse(u, u);
    u.ntt_pow_phi();

    for (size_t i = MSGS - 1; i > 0; i--)
    {
        inv[i] = u * inv[i - 1];
        u = u * t[i];
    }
    inv[0] = u;
}

static int rej_sampling_ptrs(const vector<const params::poly_q *> &zptrs,
                             const vector<const params::poly_q *> &vptrs,
                             uint64_t s2, int T)
{
    array<mpz_t, params::poly_q::degree> coeffs0, coeffs1;
    params::poly_q t;
    mpz_t dot, norm, qDivBy2, tmp;
    double r, M = exp((T * 1.14 * sqrt(0.667) * 23 * sqrt(WIDTH * N)) / (2.0 * s2));
    mpf_t u;
    uint8_t buf[8];
    int64_t seed;
    gmp_randstate_t state;

    mpf_init(u);
    gmp_randinit_mt(state);
    mpz_inits(dot, norm, qDivBy2, tmp, nullptr);
    for (size_t k = 0; k < params::poly_q::degree; ++k)
    {
        mpz_init2(coeffs0[k], (params::poly_q::bits_in_moduli_product() << 2));
        mpz_init2(coeffs1[k], (params::poly_q::bits_in_moduli_product() << 2));
    }
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);
    mpz_set_ui(norm, 0);
    mpz_set_ui(dot, 0);

    // Accumulate across all pointers, flatten
    for (size_t idx = 0; idx < zptrs.size(); ++idx)
    {
        // z
        t = *zptrs[idx];
        t.invntt_pow_invphi();
        t.poly2mpz(coeffs0);

        // v
        t = *vptrs[idx];
        t.invntt_pow_invphi();
        t.poly2mpz(coeffs1);

        // center + update dot/norm
        for (size_t k = 0; k < params::poly_q::degree; ++k)
        {
            util::center(coeffs0[k], coeffs0[k],
                         params::poly_q::moduli_product(), qDivBy2);
            util::center(coeffs1[k], coeffs1[k],
                         params::poly_q::moduli_product(), qDivBy2);

            mpz_mul(tmp, coeffs0[k], coeffs1[k]);
            mpz_add(dot, dot, tmp);

            mpz_mul(tmp, coeffs1[k], coeffs1[k]);
            mpz_add(norm, norm, tmp);
        }
    }

    getrandom(buf, sizeof buf, 0);
    memcpy(&seed, buf, sizeof seed);
    gmp_randseed_ui(state, (unsigned long)seed);
    mpf_urandomb(u, state, mpf_get_default_prec());

    r = -2.0 * mpz_get_d(dot) + mpz_get_d(norm);
    r = r / (2.0 * s2);
    r = exp(r) / M;
    int result = mpf_get_d(u) > r;

    // Cleanup
    mpf_clear(u);
    mpz_clears(dot, norm, qDivBy2, tmp, nullptr);
    for (size_t k = 0; k < params::poly_q::degree; ++k)
    {
        mpz_clear(coeffs0[k]);
        mpz_clear(coeffs1[k]);
    }
    return result;
}

// Rejection sampler on vectors of polynomials:
// returns 1 = REJECT, 0 = ACCEPT
static int rej_sampling_vec(const std::vector<params::poly_q> &z, // length WIDTH
                            const std::vector<params::poly_q> &v, // length WIDTH (typically chal * r)
                            uint64_t s2)
{
    // Basic shape checks (optional)
    if (z.size() != v.size())
    {
        // Mismatch: reject conservatively
        return 1;
    }

    std::array<mpz_t, params::poly_q::degree> coeffs0, coeffs1;
    params::poly_q t;
    mpz_t dot, norm, qDivBy2, tmp;
    double r, M = 3.18;
    int64_t seed;
    mpf_t u;
    uint8_t buf[8];
    gmp_randstate_t state;

    // Init big ints/floats
    mpf_init(u);
    gmp_randinit_mt(state);
    mpz_inits(dot, norm, qDivBy2, tmp, nullptr);
    for (size_t k = 0; k < params::poly_q::degree; ++k)
    {
        mpz_init2(coeffs0[k], (params::poly_q::bits_in_moduli_product() << 2));
        mpz_init2(coeffs1[k], (params::poly_q::bits_in_moduli_product() << 2));
    }

    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);
    mpz_set_ui(norm, 0);
    mpz_set_ui(dot, 0);

    // Accumulate over coordinates
    for (size_t i = 0; i < z.size(); ++i)
    {
        // Bring to coeff domain and extract big-int coeffs
        t = z[i];
        t.invntt_pow_invphi();
        t.poly2mpz(coeffs0);
        t = v[i];
        t.invntt_pow_invphi();
        t.poly2mpz(coeffs1);

        // Center coefficients and update dot/norm
        for (size_t k = 0; k < params::poly_q::degree; ++k)
        {
            util::center(coeffs0[k], coeffs0[k],
                         params::poly_q::moduli_product(), qDivBy2);
            util::center(coeffs1[k], coeffs1[k],
                         params::poly_q::moduli_product(), qDivBy2);

            // dot += coeffs0[k] * coeffs1[k]
            mpz_mul(tmp, coeffs0[k], coeffs1[k]);
            mpz_add(dot, dot, tmp);

            // norm += coeffs1[k]^2
            mpz_mul(tmp, coeffs1[k], coeffs1[k]);
            mpz_add(norm, norm, tmp);
        }
    }

    // Uniform u in [0,1)
    getrandom(buf, sizeof buf, 0); // (you can replace with a checked fill helper if you prefer)
    memcpy(&seed, buf, sizeof seed);
    gmp_randseed_ui(state, (unsigned long)seed);
    mpf_urandomb(u, state, mpf_get_default_prec());

    // Acceptance prob r = exp( ( -2<z,v> + ||v||^2 ) / (2 s^2) ) / M  (M=1 here)
    r = -2.0 * mpz_get_d(dot) + mpz_get_d(norm);
    r /= (2.0 * (double)s2);
    r = exp(r) / M;

    int reject = (mpf_get_d(u) > r) ? 1 : 0;

    // Cleanup
    mpf_clear(u);
    mpz_clears(dot, norm, qDivBy2, tmp, nullptr);
    for (size_t k = 0; k < params::poly_q::degree; ++k)
    {
        mpz_clear(coeffs0[k]);
        mpz_clear(coeffs1[k]);
    }
    return reject;
}
static int rej_sampling(params::poly_q z[WIDTH], params::poly_q v[WIDTH],
                        uint64_t s2)
{
    array<mpz_t, params::poly_q::degree> coeffs0, coeffs1;
    params::poly_q t;
    mpz_t dot, norm, qDivBy2, tmp;
    double r, M = 3.18;
    int64_t seed;
    mpf_t u;
    uint8_t buf[8];
    gmp_randstate_t state;
    int result;

    /// Constructors
    mpf_init(u);
    gmp_randinit_mt(state);
    mpz_inits(dot, norm, qDivBy2, tmp, nullptr);
    for (size_t i = 0; i < params::poly_q::degree; i++)
    {
        mpz_init2(coeffs0[i], (params::poly_q::bits_in_moduli_product() << 2));
        mpz_init2(coeffs1[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);
    mpz_set_ui(norm, 0);
    mpz_set_ui(dot, 0);
    for (int i = 0; i < WIDTH; i++)
    {
        t = z[i];
        t.invntt_pow_invphi();
        t.poly2mpz(coeffs0);
        t = v[i];
        t.invntt_pow_invphi();
        t.poly2mpz(coeffs1);
        for (size_t i = 0; i < params::poly_q::degree; i++)
        {
            util::center(coeffs0[i], coeffs0[i],
                         params::poly_q::moduli_product(), qDivBy2);
            util::center(coeffs1[i], coeffs1[i],
                         params::poly_q::moduli_product(), qDivBy2);
            mpz_mul(tmp, coeffs0[i], coeffs1[i]);
            mpz_add(dot, dot, tmp);
            mpz_mul(tmp, coeffs1[i], coeffs1[i]);
            mpz_add(norm, norm, tmp);
        }
    }

    getrandom(buf, sizeof(buf), 0);
    memcpy(&seed, buf, sizeof(buf));
    gmp_randseed_ui(state, seed);
    mpf_urandomb(u, state, mpf_get_default_prec());

    r = -2.0 * mpz_get_d(dot) + mpz_get_d(norm);
    r = r / (2.0 * s2);
    r = exp(r) / M;
    result = mpf_get_d(u) > r;

    mpf_clear(u);
    mpz_clears(dot, norm, qDivBy2, tmp, nullptr);
    for (size_t i = 0; i < params::poly_q::degree; i++)
    {
        mpz_clear(coeffs0[i]);
        mpz_clear(coeffs1[i]);
    }
    return result;
}

// New
static void lin_prover(params::poly_q y_E[WIDTH], params::poly_q y_S[SIZE][WIDTH],
                       params::poly_q &t_E, params::poly_q t_S[SIZE], params::poly_q u[SIZE],
                       commit_t com_E, commit_t com_S[SIZE], params::poly_q alpha[SIZE],
                       comkey_t &key, vector<params::poly_q> r_E,
                       vector<vector<params::poly_q>> r_S)
{
    params::poly_q chal, tmp[WIDTH], _tmp[SIZE][WIDTH];
    array<mpz_t, params::poly_q::degree> coeffs;
    mpz_t qDivBy2;
    int rej0, rej1;

    mpz_init(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++)
    {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

    do
    {
        /* Prover samples y,y' from Gaussian. */
        for (int i = 0; i < WIDTH; i++)
        {
            for (size_t k = 0; k < params::poly_q::degree; k++)
            {
                int64_t coeff = sample_z(0.0, SIGMA_C);
                mpz_set_si(coeffs[k], coeff);
            }
            y_E[i].mpz2poly(coeffs);
            y_E[i].ntt_pow_phi();
            for (int j = 0; j < SIZE; j++)
            {
                for (size_t k = 0; k < params::poly_q::degree; k++)
                {
                    int64_t coeff = sample_z(0.0, SIGMA_C);
                    mpz_set_si(coeffs[k], coeff);
                }
                y_S[j][i].mpz2poly(coeffs);
                y_S[j][i].ntt_pow_phi();
            }
        }

        // Compute t_E = A1*y_E and t_S = A1*y_S
        t_E = y_E[0];
        for (int i = 0; i < SIZE; i++)
        {
            t_S[i] = y_S[i][0];
        }
        for (int j = 0; j < WIDTH - HEIGHT; j++)
        {
            t_E = t_E + key.A1[0][j] * y_E[j + HEIGHT];
            for (int k = 0; k < SIZE; k++)
            {
                t_S[k] = t_S[k] + key.A1[0][j] * y_S[k][j + HEIGHT];
            }
        }

        for (int i = 0; i < SIZE; ++i)
            u[i] = 0;
        for (int k = 0; k < SIZE; ++k)
        {
            const params::poly_q alph = alpha[k];
            for (int i = 0; i < SIZE; ++i)
            {
                // row-dot: A2[i,:] · y_S[k,:]
                params::poly_q acc = key.A2[i][0] * y_S[k][0];
                for (int j = 1; j < WIDTH; ++j)
                {
                    acc = acc + key.A2[i][j] * y_S[k][j];
                }
                u[i] = u[i] + alph * acc; // scale by alpha_k and add
            }
        }

        for (int i = 0; i < SIZE; ++i)
        {
            params::poly_q acc = key.A2[i][0] * y_E[0];
            for (int j = 1; j < WIDTH; j++)
            {
                acc = acc + key.A2[i][j] * y_E[j];
            }
            u[i] = u[i] - acc; // subtract A2[i,:] · y_E
        }

        /* Sample challenge. */
        lin_hash(chal, key, com_E, com_S, alpha, u, t_E, t_S); // TODO: send whole u and t

        /* Prover */
        for (int i = 0; i < WIDTH; i++)
        {
            // Compute d*r and d*r'
            tmp[i] = chal * r_E[i];
            for (int j = 0; j < SIZE; j++)
            {
                _tmp[j][i] = chal * r_S[j][i];
            }
            // Compute z = y + d*r and z' = y' + d*r'
            y_E[i] = y_E[i] + tmp[i];
            for (int j = 0; j < SIZE; j++)
            {
                y_S[j][i] = y_S[j][i] + _tmp[j][i];
            }
        }

        std::vector<const params::poly_q *> zptrs;
        std::vector<const params::poly_q *> vptrs;
        zptrs.reserve(WIDTH * (1 + SIZE));
        vptrs.reserve(WIDTH * (1 + SIZE));

        for (int row = 0; row < SIZE; ++row)
        {
            zptrs.push_back(&y_E[row]);
            vptrs.push_back(&tmp[row]);
            for (int col = 0; col < WIDTH; ++col)
            {
                zptrs.push_back(&y_S[row][col]);
                vptrs.push_back(&_tmp[row][col]);
            }
        }

        double sigma_c = 0.675 * (SIZE + 1) * 1.14 * sqrt(0.667) * 23 * sqrt(WIDTH * N);
        rej1 = rej_sampling_ptrs(zptrs, vptrs, (uint64_t)sigma_c * (uint64_t)sigma_c, SIZE + 1);
    } while (rej0 || rej1);

    for (size_t i = 0; i < params::poly_q::degree; i++)
    {
        mpz_clear(coeffs[i]);
    }
    mpz_clear(qDivBy2);
}

// New
static int lin_verifier(params::poly_q z_E[WIDTH], params::poly_q z_S[SIZE][WIDTH],
                        params::poly_q t_E, params::poly_q t_S[SIZE], params::poly_q u[SIZE],
                        commit_t com_E, commit_t com_S[SIZE], params::poly_q alpha[SIZE], comkey_t &key, params::poly_q beta[SIZE])
{
    params::poly_q chal, v, tmp, zero = 0;
    params::poly_q _v[SIZE], _tmp[SIZE];
    int result = 1;

    /* Sample challenge. */
    lin_hash(chal, key, com_E, com_S, alpha, u, t_E, t_S); // TODO: send whole u

    /* Verifier checks norm, reconstruct from NTT representation. */
    for (int i = 0; i < WIDTH; i++)
    {
        v = z_E[i];
        v.invntt_pow_invphi();
        result &= bdlop_test_norm(v, 4 * SIGMA_C * SIGMA_C);
        for (int j = 0; j < SIZE; j++)
        {
            v = z_S[j][i];
            v.invntt_pow_invphi();
            result &= bdlop_test_norm(v, 4 * SIGMA_C * SIGMA_C);
        }
    }

    /* Verifier computes A1z and A1z'. */
    v = z_E[0];
    // _v = _z[0];
    for (int i = 0; i < HEIGHT; i++)
    {
        for (int j = 0; j < WIDTH - HEIGHT; j++)
        {
            v = v + key.A1[i][j] * z_E[j + HEIGHT];
            // _v = _v + key.A1[i][j] * _z[j + HEIGHT];
        }
    }
    for (int b = 0; b < SIZE; b++)
    {
        _v[b] = z_S[b][0];
        for (int j = 0; j < WIDTH - HEIGHT; j++)
        {
            _v[b] = _v[b] + key.A1[0][j] * z_S[b][j + HEIGHT];
        }
    }

    tmp = t_E + chal * com_E.c1 - v;
    tmp.invntt_pow_invphi();
    result &= (tmp == zero);
    for (int i = 0; i < SIZE; i++)
    {
        _tmp[i] = 0;
        _tmp[i] = t_S[i] + chal * com_S[i].c1 - _v[i];
        _tmp[i].invntt_pow_invphi();
        result &= (_tmp[i] == zero);
    }

    // check lin
    for (int i = 0; i < SIZE; ++i)
    {
        // sum_k alpha[k] * (A2[i,:] · z_S[k,:])
        params::poly_q sumAlphaA2zS = 0;
        for (int k = 0; k < SIZE; ++k)
        {
            params::poly_q dot = 0;
            // row-dot: A2[i,:] · z_S[k,:]
            dot = key.A2[i][0] * z_S[k][0];
            for (int j = 1; j < WIDTH; ++j)
                dot = dot + key.A2[i][j] * z_S[k][j];

            sumAlphaA2zS = sumAlphaA2zS + alpha[k] * dot;
        }

        // Az[i] = A2[i,:] · z_E
        params::poly_q Az_i = key.A2[i][0] * z_E[0];
        for (int j = 1; j < WIDTH; ++j)
            Az_i = Az_i + key.A2[i][j] * z_E[j];

        // ac2[i] = Σ_j alpha[j] * com_S[j].c2[i]
        params::poly_q ac2_i = 0;
        for (int j = 0; j < SIZE; ++j)
            ac2_i = ac2_i + alpha[j] * com_S[j].c2[i];

        // diff = (LHS - u) - chal * (ac2 + beta - com_E.c2)
        params::poly_q lhs_minus_u = (sumAlphaA2zS - Az_i) - u[i];
        params::poly_q rhs_known = (ac2_i + beta[i] - com_E.c2[i]) * chal;
        params::poly_q diff = lhs_minus_u - rhs_known;

        // Compare in coeff domain
        diff.invntt_pow_invphi();
        result &= (diff == zero);
    }
    return result;
}


#ifdef MAIN

inline params::poly_q make_scale_poly(uint64_t lambda = 1)
{
    const char *q_str = "633544653872304603170707504554316801";

    const mpz_t &qq = params::poly_q::moduli_product(); // NFLlib modulus product

    mpz_t q, invq, lam, s;
    mpz_init_set_str(q, q_str, 10); // load q
    mpz_init(invq);
    mpz_init(lam);
    mpz_init(s);

    // invq = q^{-1} mod Q
    if (mpz_invert(invq, q, qq) == 0)
    {
        throw std::runtime_error("q is not invertible modulo NFL Q");
    }

    mpz_set_ui(lam, lambda);
    mpz_mul(s, lam, invq);
    mpz_mod(s, s, qq); // s = lambda * q^{-1} mod Q

    // Load s as a constant polynomial
    std::array<mpz_t, params::poly_q::degree> coeffs;
    mpz_init(coeffs[0]);
    mpz_set(coeffs[0], s);
    for (size_t k = 1; k < coeffs.size(); ++k)
    {
        mpz_init(coeffs[k]);
        mpz_set_ui(coeffs[k], 0);
    }

    params::poly_q scale;
    scale.mpz2poly(coeffs); // coeff domain
    scale.ntt_pow_phi();    // convert to NTT

    for (auto &c : coeffs)
        mpz_clear(c);
    mpz_clear(q);
    mpz_clear(invq);
    mpz_clear(lam);
    mpz_clear(s);
    return scale;
}

inline params::poly_q make_q_inv_const()
{
    const char *q_str = "633544653872304603170707504554316801";
    const mpz_t &qq = params::poly_q::moduli_product();

    mpz_t q, invq;
    mpz_init_set_str(q, q_str, 10);
    mpz_init(invq);

    if (mpz_invert(invq, q, qq) == 0)
    {
        mpz_clear(q);
        mpz_clear(invq);
        throw std::runtime_error("q is not invertible modulo NFL Q");
    }

    // Build constant polynomial with value invq
    std::array<mpz_t, params::poly_q::degree> coeffs;
    mpz_init(coeffs[0]);
    mpz_set(coeffs[0], invq);
    for (size_t k = 1; k < coeffs.size(); ++k)
    {
        mpz_init(coeffs[k]);
        mpz_set_ui(coeffs[k], 0);
    }

    params::poly_q qinv;
    qinv.mpz2poly(coeffs); // coeff domain constant
    qinv.ntt_pow_phi();    // NTT constant

    for (auto &c : coeffs)
        mpz_clear(c);
    mpz_clear(q);
    mpz_clear(invq);
    return qinv;
}

static void matvec_mul(vector<params::poly_q> result, vector<vector<params::poly_q>> M, vector<params::poly_q> v, int size)
{
}

inline void mul_alpha_c2(vector<params::poly_q> &result, vector<params::poly_q> &c2, params::poly_q &alpha)
{
    assert(result.size() == c2.size());
    for (int i = 0; i < result.size(); i++)
    {
        result[i] = alpha * c2[i];
    }
}

inline void mul_A2_y(vector<params::poly_q> &result, vector<vector<params::poly_q>> &A2, vector<params::poly_q> &y)
{
    assert(result.size() == A2.size() && A2[0].size() == y.size());
    for (int i = 0; i < A2.size(); i++)
    {
        params::poly_q acc = A2[i][0] * y[0];
        for (int j = 1; j < A2[i].size(); j++)
        {
            acc = acc + A2[i][j] * y[j];
        }
        result[i] = acc;
    }
}

static void bench_lin()
{
    comkey_t key;
    commit_t com_E;
    commit_t com_S[SIZE];
    params::poly_q alpha[SIZE], u[SIZE], t_E, t_S[SIZE], y_E[WIDTH], y_S[SIZE][WIDTH];
    vector<params::poly_q> r_E(WIDTH), s(WIDTH);                             // Randomness vectors.
    vector<vector<params::poly_q>> r_S(SIZE, vector<params::poly_q>(WIDTH)); // Randomness vectors for secret matrix.
    vector<params::poly_q> E(SIZE);                                          // noise.
    vector<params::poly_q> S(SIZE);                                          // secret

    vector<vector<params::poly_q>> S_cols(SIZE, vector<params::poly_q>(SIZE)); // secret matrix
    params::poly_q uvec[SIZE];
    params::poly_q ds[SIZE], beta[SIZE];

    std::array<mpz_t, params::poly_q::degree> coeffs;
    for (size_t i = 0; i < params::poly_q::degree; i++)
    {
        mpz_init(coeffs[i]);
    }
    int64_t coeff;

    int lambda = 1;
    params::poly_q alpha_scale = make_scale_poly(lambda);
    params::poly_q invq = make_q_inv_const();

    params::poly_q lambda_poly = lambda;
    lambda_poly.ntt_pow_phi();
    const char *q_str = "633544653872304603170707504554316801";
    std::array<mpz_t, params::poly_q::degree> q_coeffs;
    for (size_t k = 0; k < q_coeffs.size(); k++)
    {
        mpz_init(q_coeffs[k]);
    }
    mpz_set_str(q_coeffs[0], q_str, 10);
    for (size_t k = 1; k < q_coeffs.size(); k++)
    {
        mpz_set_ui(q_coeffs[k], 0);
    }
    params::poly_q q_poly;
    q_poly.mpz2poly(q_coeffs);
    q_poly.ntt_pow_phi();

    for (auto &c : q_coeffs)
    {
        mpz_clear(c);
    }

    // Sample S_cols
    for (int i = 0; i < SIZE; i++)
    {
        uvec[i] = nfl::uniform();
        uvec[i].ntt_pow_phi();
        for (int j = 0; j < SIZE; j++)
        {
            S_cols[i][j] = nfl::ZO_dist();
            S_cols[i][j].ntt_pow_phi();
        }
    }

    // Sample E with coeffs gaussian 1208
    for (int i = 0; i < SIZE; i++)
    {
        for (int j = 0; j < params::poly_q::degree; j++)
        {
            int64_t coeff = sample_z(0.0, SIGMAE);
            mpz_set_si(coeffs[j], coeff);
        }
        E[i].mpz2poly(coeffs);
        E[i].ntt_pow_phi();
    }
    params::poly_q zero = 0;

    bdlop_keygen(key);
    bdlop_sample_rand(r_E);
    bdlop_commit(com_E, E, key, r_E);

    for (int i = 0; i < SIZE; i++)
    {
        bdlop_sample_rand(r_S[i]);
        bdlop_commit(com_S[i], S_cols[i], key, r_S[i]);
    }

    for (int i = 0; i < SIZE; i++)
    {
        alpha[i] = alpha_scale * uvec[i]; // nfl::uniform();
        // alpha[i].ntt_pow_phi();
    }
    for (int i = 0; i < SIZE; ++i)
    {
        params::poly_q sumAlphaS = 0;
        for (int j = 0; j < SIZE; ++j)
        {
            sumAlphaS = sumAlphaS + alpha[j] * S_cols[j][i];
        }
        beta[i] = E[i] - sumAlphaS;
    }

    // Set ds = S_cols * uvec + E

    // bdlop_sample_chal(beta);

    // BENCH_BEGIN("linear hash")
    // {
    //     BENCH_ADD(lin_hash(beta, key, com_E, com_S, alpha, u, t, _t));
    // }
    // BENCH_END;

    printf("\nSign round 3 linearity proof:\n\n");
    BENCH_BEGIN("linear proof")
    {
        BENCH_ADD(lin_prover(y_E, y_S, t_E, t_S, u, com_E, com_S, alpha, key, r_E, r_S));
    }
    BENCH_END;
    // lin_prover(y_E, y_S, t_E, t_S, u, com_E, com_S, alpha, key, r_E, r_S);
    // int ver = lin_verifier(y_E, y_S, t_E, t_S, u, com_E, com_S, alpha, key, beta);

    BENCH_BEGIN("linear verifier")
    {
        BENCH_ADD(lin_verifier(y_E, y_S, t_E, t_S, u, com_E, com_S, alpha, key, beta));
    }
    BENCH_END;
}

// NFLlib exposes only non-const data(), so we read bytes via a safe const_cast.
// We do NOT mutate, just treat the memory as const bytes.
static inline void hash_poly_bytes(blake3_hasher *h, const params::poly_q &p)
{
    auto *ptr = const_cast<params::poly_q &>(p).data(); // uint64_t*
    constexpr size_t BYTES_PER_POLY = 16 * DEGREE;      // 2 moduli * 8 bytes * DEGREE
    blake3_hasher_update(h, reinterpret_cast<const uint8_t *>(ptr), BYTES_PER_POLY);
}
static void lin_hash_dkg(params::poly_q &chal,
                         comkey_t &key,
                         commit_t &x,
                         vector<commit_t> &y,
                         vector<params::poly_q> &alpha,
                         vector<params::poly_q> &u,
                         params::poly_q &t,
                         vector<params::poly_q> &_t)
{
    uint8_t hash[BLAKE3_OUT_LEN];
    blake3_hasher hasher;
    blake3_hasher_init(&hasher);

    // key.A1 / key.A2
    for (size_t i = 0; i < HEIGHT; ++i)
        for (size_t j = 0; j < WIDTH - HEIGHT; ++j)
            hash_poly_bytes(&hasher, key.A1[i][j]);

    for (size_t i = 0; i < SIZE; ++i)
        for (size_t j = 0; j < WIDTH; ++j)
            hash_poly_bytes(&hasher, key.A2[i][j]);

    // encode sizes to lock transcript
    auto mix_u64 = [&](uint64_t v)
    { blake3_hasher_update(&hasher, (const uint8_t *)&v, sizeof v); };
    mix_u64((uint64_t)alpha.size());
    mix_u64((uint64_t)y.size());
    mix_u64((uint64_t)x.c2.size());
    for (const auto &yy : y)
        mix_u64((uint64_t)yy.c2.size());
    mix_u64((uint64_t)u.size());
    mix_u64((uint64_t)_t.size());

    for (size_t i = 0; i < alpha.size(); ++i)
        hash_poly_bytes(&hasher, alpha[i]);

    hash_poly_bytes(&hasher, x.c1);
    for (size_t i = 0; i < y.size(); ++i)
        hash_poly_bytes(&hasher, y[i].c1);

    for (size_t i = 0; i < x.c2.size(); ++i)
        hash_poly_bytes(&hasher, x.c2[i]);
    for (size_t j = 0; j < y.size(); ++j)
        for (size_t i = 0; i < y[j].c2.size(); ++i)
            hash_poly_bytes(&hasher, y[j].c2[i]);

    for (size_t i = 0; i < u.size(); ++i)
        hash_poly_bytes(&hasher, u[i]);
    hash_poly_bytes(&hasher, t);
    for (size_t i = 0; i < _t.size(); ++i)
        hash_poly_bytes(&hasher, _t[i]);

    blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);
    nfl::fastrandombytes_seed(hash);
    bdlop_sample_chal(chal);
    nfl::fastrandombytes_reseed();
}

static void lin_hash_dkg2(params::poly_q &chal,
                          comkey_t &key,
                          vector<commit_t> &y,
                          vector<params::poly_q> &codewordi,
                          vector<params::poly_q> &u,
                          vector<params::poly_q> &_t)
{
    uint8_t hash[BLAKE3_OUT_LEN];
    blake3_hasher hasher;
    blake3_hasher_init(&hasher);

    // key.A1 / key.A2
    for (size_t i = 0; i < HEIGHT; ++i)
        for (size_t j = 0; j < WIDTH - HEIGHT; ++j)
            hash_poly_bytes(&hasher, key.A1[i][j]);

    for (size_t i = 0; i < SIZE; ++i)
        for (size_t j = 0; j < WIDTH; ++j)
            hash_poly_bytes(&hasher, key.A2[i][j]);

    // encode sizes to lock transcript
    auto mix_u64 = [&](uint64_t v)
    { blake3_hasher_update(&hasher, (const uint8_t *)&v, sizeof v); };
    mix_u64((uint64_t)codewordi.size());
    mix_u64((uint64_t)y.size());
    for (const auto &yy : y)
        mix_u64((uint64_t)yy.c2.size());
    mix_u64((uint64_t)u.size());
    mix_u64((uint64_t)_t.size());

    for (size_t i = 0; i < codewordi.size(); ++i)
        hash_poly_bytes(&hasher, codewordi[i]);

    for (size_t i = 0; i < y.size(); ++i)
        hash_poly_bytes(&hasher, y[i].c1);

    for (size_t j = 0; j < y.size(); ++j)
        for (size_t i = 0; i < y[j].c2.size(); ++i)
            hash_poly_bytes(&hasher, y[j].c2[i]);

    for (size_t i = 0; i < u.size(); ++i)
        hash_poly_bytes(&hasher, u[i]);
    for (size_t i = 0; i < _t.size(); ++i)
        hash_poly_bytes(&hasher, _t[i]);

    blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);
    nfl::fastrandombytes_seed(hash);
    bdlop_sample_chal(chal);
    nfl::fastrandombytes_reseed();
}

static void lin_prover_dkg(
    vector<params::poly_q> &y_S,
    vector<vector<params::poly_q>> &y_Si,
    params::poly_q &t_S,
    vector<params::poly_q> &t_Si,
    vector<params::poly_q> &u,
    commit_t &com_S,
    vector<commit_t> &com_Si,
    vector<params::poly_q> &lambdai,
    comkey_t &key,
    vector<params::poly_q> &r_S,
    vector<vector<params::poly_q>> &r_Si, int THRESHOLD)
{
    params::poly_q chal;
    vector<params::poly_q> tmp(WIDTH);
    vector<vector<params::poly_q>> tmpi(THRESHOLD, vector<params::poly_q>(WIDTH));
    array<mpz_t, params::poly_q::degree> coeffs;
    mpz_t qDivBy2;
    int rej0, rej1;
    mpz_init(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++)
    {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

    do
    {
        // Sample y_S, y_Si
        for (int i = 0; i < WIDTH; i++)
        {
            for (size_t k = 0; k < params::poly_q::degree; k++)
            {
                int64_t coeff = sample_z(0.0, SIGMA_C);
                mpz_set_si(coeffs[k], coeff);
            }
            y_S[i].mpz2poly(coeffs);
            y_S[i].ntt_pow_phi();
            for (int j = 0; j < THRESHOLD; j++)
            {
                for (size_t k = 0; k < params::poly_q::degree; k++)
                {
                    int64_t coeff = sample_z(0.0, SIGMA_C);
                    mpz_set_si(coeffs[k], coeff);
                }
                y_Si[j][i].mpz2poly(coeffs);
                y_Si[j][i].ntt_pow_phi();
            }
        }

        // Compute t_S, t_Si
        t_S = y_S[0];
        for (int i = 0; i < THRESHOLD; i++)
        {
            t_Si[i] = y_Si[i][0];
        }
        for (int i = 0; i < WIDTH - HEIGHT; i++)
        {
            t_S = t_S + key.A1[0][i] * y_S[i + HEIGHT];
            for (int j = 0; j < THRESHOLD; j++)
            {
                t_Si[j] = t_Si[j] + key.A1[0][i] * y_Si[j][i + HEIGHT];
            }
        }

        for (size_t i = 0; i < SIZE; ++i)
        {
            params::poly_q acc_row = 0;

            for (size_t t = 0; t < THRESHOLD; ++t)
            {
                params::poly_q dot = key.A2[i][0] * y_Si[t][0];
                for (size_t j = 1; j < WIDTH; ++j)
                    dot = dot + key.A2[i][j] * y_Si[t][j];
                acc_row = acc_row + lambdai[t] * dot;
            }

            params::poly_q dot_yS = key.A2[i][0] * y_S[0];
            for (size_t j = 1; j < WIDTH; ++j)
                dot_yS = dot_yS + key.A2[i][j] * y_S[j];

            u[i] = acc_row - dot_yS;
        }

        lin_hash_dkg(chal, key, com_S, com_Si, lambdai, u, t_S, t_Si);

        for (int i = 0; i < WIDTH; i++)
        {
            tmp[i] = chal * r_S[i];
            for (int j = 0; j < THRESHOLD; j++)
            {
                tmpi[j][i] = chal * r_Si[j][i];
            }
            y_S[i] = y_S[i] + tmp[i];
            for (int j = 0; j < THRESHOLD; j++)
            {
                y_Si[j][i] = y_Si[j][i] + tmpi[j][i];
            }
        }
        // rej0 = rej_sampling_vec(y_S, tmp, (uint64_t)SIGMA_C * (uint64_t)SIGMA_C);
        vector<const params::poly_q *> zptrs;
        vector<const params::poly_q *> vptrs;
        zptrs.reserve(WIDTH * (1 + THRESHOLD));
        vptrs.reserve(WIDTH * (1 + THRESHOLD));

        for (int i = 0; i < WIDTH; ++i)
        {
            vptrs.push_back(&tmp[i]);
            zptrs.push_back(&y_S[i]);
            for (int j = 0; j < THRESHOLD; ++j)
            {
                vptrs.push_back(&tmpi[j][i]);
                zptrs.push_back(&y_Si[j][i]);
            }
        }

        double sigma_c = 0.675 * (THRESHOLD + 1) * 1.14 * sqrt(0.667) * 23 * sqrt(WIDTH * N);
        rej1 = rej_sampling_ptrs(zptrs, vptrs, (uint64_t)sigma_c * (uint64_t)sigma_c, THRESHOLD + 1);
    } while (rej0 || rej1);
    mpz_clear(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++)
    {
        mpz_clear(coeffs[i]);
    }
}

static int lin_verifier_dkg(
    vector<params::poly_q> &z_S,
    vector<vector<params::poly_q>> &z_Si,
    params::poly_q &t_S,
    vector<params::poly_q> &t_Si,
    vector<params::poly_q> &u,
    commit_t &com_S,
    vector<commit_t> &com_Si,
    vector<params::poly_q> &lambdai,
    comkey_t &key,
    vector<params::poly_q> &beta, int THRESHOLD)
{
    int result = 1;
    params::poly_q chal, v, tmp, zero = 0;
    vector<params::poly_q> _v(THRESHOLD);
    // vector<params::poly_q> lhs(WIDTH), rhs(WIDTH);
    vector<params::poly_q> _tmp(THRESHOLD);

    lin_hash_dkg(chal, key, com_S, com_Si, lambdai, u, t_S, t_Si);

    for (int i = 0; i < WIDTH; i++)
    {
        v = z_S[i];
        v.invntt_pow_invphi();
        result &= bdlop_test_norm(v, 4 * SIGMA_C * SIGMA_C);
        for (int j = 0; j < THRESHOLD; j++)
        {
            v = z_Si[j][i];
            v.invntt_pow_invphi();
            result &= bdlop_test_norm(v, 4 * SIGMA_C * SIGMA_C);
        }
    }

    /* Verifier computes A1z and A1z'. */
    v = z_S[0];
    for (int j = 0; j < WIDTH - HEIGHT; j++)
    {
        v = v + key.A1[0][j] * z_S[j + HEIGHT];
    }
    for (int b = 0; b < THRESHOLD; b++)
    {
        _v[b] = z_Si[b][0];
        for (int j = 0; j < WIDTH - HEIGHT; j++)
        {
            _v[b] = _v[b] + key.A1[0][j] * z_Si[b][j + HEIGHT];
        }
    }

    tmp = t_S + chal * com_S.c1 - v;
    tmp.invntt_pow_invphi();
    result &= (tmp == zero);
    for (int i = 0; i < THRESHOLD; i++)
    {
        _tmp[i] = 0;
        _tmp[i] = t_Si[i] + chal * com_Si[i].c1 - _v[i];
        _tmp[i].invntt_pow_invphi();
        result &= (_tmp[i] == zero);
    }

    for (size_t i = 0; i < SIZE; ++i)
    {
        // LHS: sum_t lambda[t]*(A2*z_Si[t])  -  (A2*z_S)
        params::poly_q lhs = 0;
        for (size_t t = 0; t < THRESHOLD; ++t)
        {
            params::poly_q dot = key.A2[i][0] * z_Si[t][0];
            for (size_t j = 1; j < WIDTH; ++j)
                dot = dot + key.A2[i][j] * z_Si[t][j];
            lhs = lhs + (lambdai[t] * dot);
        }
        params::poly_q dot_zS = key.A2[i][0] * z_S[0];
        for (size_t j = 1; j < WIDTH; ++j)
            dot_zS = dot_zS + key.A2[i][j] * z_S[j];
        lhs = lhs - dot_zS;

        // RHS: (sum_t lambda[t]*com_Si[t].c2[i] + beta[i] - com_S.c2[i]) * chal + u[i]
        params::poly_q c2_combo = 0;
        for (size_t t = 0; t < THRESHOLD; ++t)
            c2_combo = c2_combo + (lambdai[t] * com_Si[t].c2[i]);
        c2_combo = c2_combo - com_S.c2[i];

        params::poly_q rhs = (c2_combo)*chal + u[i];

        params::poly_q diff = lhs - rhs;
        diff.invntt_pow_invphi();
        result &= (diff == zero);
    }
    return result;
}

void bench_lin_dkg3(int THRESHOLD)
{
    cout << "\nDKG lin proof 3 with" << " THRESHOLD = " << THRESHOLD << "\n"
         << endl;
    ;
    comkey_t key;
    commit_t com_S;
    vector<commit_t> com_Si(THRESHOLD);
    vector<vector<params::poly_q>> Si(THRESHOLD, vector<params::poly_q>(SIZE)); // Secret shares (1 columns of shares)
    vector<params::poly_q> S(SIZE);                                             // Secret (1 column)
    vector<params::poly_q> lambdai(THRESHOLD);

    // Set lambdai
    std::array<mpz_t, params::poly_q::degree> coeffs;
    for (size_t k = 0; k < coeffs.size(); k++)
    {
        mpz_init(coeffs[k]);
    }
    for (int i = 1; i < THRESHOLD + 1; i++)
    {
        mpz_set_ui(coeffs[0], i);
        for (size_t k = 1; k < coeffs.size(); k++)
        {
            mpz_set_ui(coeffs[k], 0);
        }
        lambdai[i - 1].mpz2poly(coeffs);
        lambdai[i - 1].ntt_pow_phi();
    }

    // Sample Si's
    for (int i = 0; i < THRESHOLD; i++)
    {
        for (int j = 0; j < SIZE; j++)
        {
            Si[i][j] = nfl::uniform();
            Si[i][j].ntt_pow_phi();
        }
    }
    // Compute S = Si * lambdai
    for (int i = 0; i < SIZE; i++)
    {
        params::poly_q acc = Si[0][i] * lambdai[0];
        for (int j = 1; j < THRESHOLD; j++)
        {
            acc = acc + Si[j][i] * lambdai[j];
        }
        S[i] = acc;
    }

    // Commitment randomness
    vector<params::poly_q> r_S(WIDTH);
    vector<vector<params::poly_q>> r_Si(THRESHOLD, vector<params::poly_q>(WIDTH));
    bdlop_sample_rand(r_S);
    for (int i = 0; i < THRESHOLD; i++)
    {
        bdlop_sample_rand(r_Si[i]);
    }
    // commitment keygen
    bdlop_keygen(key);
    // Commit to S and Si
    bdlop_commit(com_S, S, key, r_S);
    for (int i = 0; i < THRESHOLD; i++)
    {
        bdlop_commit(com_Si[i], Si[i], key, r_Si[i]);
    }
    // Prove
    vector<params::poly_q> y_S(WIDTH);
    vector<vector<params::poly_q>> y_Si(THRESHOLD, vector<params::poly_q>(WIDTH));
    params::poly_q t_S;
    vector<params::poly_q> t_Si(THRESHOLD);
    vector<params::poly_q> u(SIZE);

    BENCH_BEGIN("linear proof"){
        // BENCH_ADD(lin_prover(y_E, y_S, t_E, t_S, u, com_E, com_S, alpha, key, r_E, r_S));
        BENCH_ADD(lin_prover_dkg(y_S, y_Si, t_S, t_Si, u, com_S, com_Si, lambdai, key, r_S, r_Si, THRESHOLD))} BENCH_END;
    BENCH_BEGIN("verify")
    {
        BENCH_ADD(lin_verifier_dkg(y_S, y_Si, t_S, t_Si, u, com_S, com_Si, lambdai, key, S, THRESHOLD));
    }
    BENCH_END;
}

static inline params::poly_q sample_const_poly_from_uniform()
{
    params::poly_q p = nfl::uniform(); // coeff-domain random poly

    // Extract coeffs to mpz_t[], zero out all but coeff 0
    std::array<mpz_t, params::poly_q::degree> coeffs;
    for (auto &c : coeffs)
        mpz_init2(c, (params::poly_q::bits_in_moduli_product() << 2));

    p.poly2mpz(coeffs); // coeff domain
    for (size_t k = 1; k < coeffs.size(); ++k)
        mpz_set_ui(coeffs[k], 0);

    // Rebuild the constant polynomial and ntt
    params::poly_q cpoly;
    cpoly.mpz2poly(coeffs);
    cpoly.ntt_pow_phi();

    for (auto &c : coeffs)
        mpz_clear(c);
    return cpoly;
}

static void lin_prover_dkg2(
    vector<vector<params::poly_q>> &y_Si,
    vector<params::poly_q> &t_Si,
    vector<params::poly_q> &u,
    vector<commit_t> &com_Si,
    comkey_t &key,
    vector<vector<params::poly_q>> &r_Si,
    vector<params::poly_q> &codewordi,
    int FULL)
{

    params::poly_q chal;
    vector<params::poly_q> tmp(WIDTH);
    vector<vector<params::poly_q>> tmpi(FULL, vector<params::poly_q>(WIDTH));
    array<mpz_t, params::poly_q::degree> coeffs;
    mpz_t qDivBy2;
    int rej0, rej1;
    mpz_init(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++)
    {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

    do
    {
        // Sample y_Si
        for (int i = 0; i < WIDTH; i++)
        {
            for (int j = 0; j < FULL; j++)
            {
                for (size_t k = 0; k < params::poly_q::degree; k++)
                {
                    int64_t coeff = sample_z(0.0, SIGMA_C);
                    mpz_set_si(coeffs[k], coeff);
                }
                y_Si[j][i].mpz2poly(coeffs);
                y_Si[j][i].ntt_pow_phi();
            }
        }

        // Compute t_Si
        for (int i = 0; i < FULL; i++)
        {
            t_Si[i] = y_Si[i][0];
        }
        for (int i = 0; i < WIDTH - HEIGHT; i++)
        {
            for (int j = 0; j < FULL; j++)
            {
                t_Si[j] = t_Si[j] + key.A1[0][i] * y_Si[j][i + HEIGHT];
            }
        }

        // Compute u
        for (int r = 0; r < SIZE; r++)
        {
            params::poly_q acc = 0;
            for (int t = 0; t < FULL; t++)
            {
                params::poly_q dot = key.A2[r][0] * y_Si[t][0];
                for (int j = 1; j < WIDTH; j++)
                    dot = dot + key.A2[r][j] * y_Si[t][j];
                acc = acc + codewordi[t] * dot;
            }
            u[r] = acc;
        }

        // Compute z_Si

        lin_hash_dkg2(chal, key, com_Si, codewordi, u, t_Si);

        for (int i = 0; i < WIDTH; i++)
        {
            for (int j = 0; j < FULL; j++)
            {
                tmpi[j][i] = chal * r_Si[j][i];
            }
            for (int j = 0; j < FULL; j++)
            {
                y_Si[j][i] = y_Si[j][i] + tmpi[j][i];
            }
        }

        rej1 = 0;
        vector<const params::poly_q *> zptrs;
        vector<const params::poly_q *> vptrs;
        zptrs.reserve(WIDTH * (FULL));
        vptrs.reserve(WIDTH * (FULL));

        for (int i = 0; i < WIDTH; ++i)
        {
            for (int j = 0; j < FULL; ++j)
            {
                zptrs.push_back(&y_Si[j][i]);
                vptrs.push_back(&tmpi[j][i]);
            }
        }
        double sigma_c = 0.675 * FULL * 1.14 * sqrt(0.667) * 23 * sqrt(WIDTH * N);
        rej1 = rej_sampling_ptrs(zptrs, vptrs, (uint64_t)sigma_c * (uint64_t)sigma_c, FULL);
    } while (rej1);
    mpz_clear(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++)
    {
        mpz_clear(coeffs[i]);
    }
}

static int lin_verifier_dkg2(
    vector<vector<params::poly_q>> &z_Si,
    vector<params::poly_q> &t_Si,
    vector<params::poly_q> &u,
    vector<commit_t> &com_Si,
    vector<params::poly_q> &codewordi,
    vector<params::poly_q> &W,
    comkey_t &key,
    int THRESHOLD)
{
    int result = 1;
    params::poly_q chal, v, tmp, zero = 0;
    vector<params::poly_q> _v(THRESHOLD);
    // vector<params::poly_q> lhs(WIDTH), rhs(WIDTH);
    vector<params::poly_q> _tmp(THRESHOLD);

    lin_hash_dkg2(chal, key, com_Si, codewordi, u, t_Si);

    for (int i = 0; i < WIDTH; i++)
    {
        for (int j = 0; j < THRESHOLD; j++)
        {
            v = z_Si[j][i];
            v.invntt_pow_invphi();
            result &= bdlop_test_norm(v, 4 * SIGMA_C * SIGMA_C);
        }
    }

    /* Verifier computes A1z and A1z'. */
    for (int b = 0; b < THRESHOLD; b++)
    {
        _v[b] = z_Si[b][0];
        for (int j = 0; j < WIDTH - HEIGHT; j++)
        {
            _v[b] = _v[b] + key.A1[0][j] * z_Si[b][j + HEIGHT];
        }
    }

    for (int i = 0; i < THRESHOLD; i++)
    {
        _tmp[i] = 0;
        _tmp[i] = t_Si[i] + chal * com_Si[i].c1 - _v[i];
        _tmp[i].invntt_pow_invphi();
        result &= (_tmp[i] == zero);
    }

    for (size_t r = 0; r < SIZE; r++)
    {
        // LHS := sum_t codewordi[t] * (A2 · z_Si[t])
        params::poly_q LHS = 0;
        for (int t = 0; t < THRESHOLD; t++)
        {
            params::poly_q dot = key.A2[r][0] * z_Si[t][0];
            for (size_t j = 1; j < WIDTH; ++j)
            {
                dot = dot + key.A2[r][j] * z_Si[t][j];
            }
            LHS = LHS + (codewordi[t] * dot);
        }

        params::poly_q c2_sum = 0;
        for (size_t t = 0; t < (size_t)THRESHOLD; t++)
        {
            c2_sum = c2_sum + (codewordi[t] * com_Si[t].c2[r]);
        }
        params::poly_q RHS = (c2_sum - W[r]) * chal;
        RHS = RHS + u[r];

        // Compare LHS_row - RHS_row == 0 (in coeff domain)
        params::poly_q diff = LHS - RHS;
        diff.invntt_pow_invphi();
        result &= (diff == zero);
    }
    return result;
}

void bench_lin_dkg2(int FULL)
{
    cout << "\nDKG lin proof 2 with" << " FULL_AMOUNT = " << FULL << "\n"
         << endl;
    ;
    comkey_t key;
    vector<commit_t> com_Si(FULL);
    vector<vector<params::poly_q>> Si(FULL, vector<params::poly_q>(SIZE)); // Secret shares (1 columns of shares)
    vector<params::poly_q> W(SIZE);                                        // Secret (1 column)
    vector<params::poly_q> codewordi(FULL);

    // Sample Si's
    for (int i = 0; i < FULL; i++)
    {
        for (int j = 0; j < SIZE; j++)
        {
            Si[i][j] = nfl::uniform();
            Si[i][j].ntt_pow_phi();
        }
    }

    // Sample codewordi
    for (int i = 0; i < FULL; i++)
    {
        // Should in actuality sample codewords here, not uniform constant polynomials
        codewordi[i] = sample_const_poly_from_uniform();
    }

    // Compute dummy W used for benchmarking only (not in real protocol, but gives same amount of computations)
    for (int j = 0; j < SIZE; j++)
    {
        params::poly_q acc = 0; // accumulator for this column
        for (int i = 0; i < FULL; i++)
        {
            acc = acc + codewordi[i] * Si[i][j];
        }
        W[j] = acc;
    }

    // Commit to Si:
    vector<vector<params::poly_q>> r_Si(FULL, vector<params::poly_q>(WIDTH));
    vector<vector<params::poly_q>> y_Si(FULL, vector<params::poly_q>(WIDTH));
    vector<params::poly_q> t_Si(FULL);
    vector<params::poly_q> u(SIZE);
    bdlop_keygen(key);
    for (int i = 0; i < FULL; i++)
    {
        bdlop_sample_rand(r_Si[i]);
        bdlop_commit(com_Si[i], Si[i], key, r_Si[i]);
    }
    BENCH_BEGIN("linear proof")
    {
        BENCH_ADD(lin_prover_dkg2(y_Si, t_Si, u, com_Si, key, r_Si, codewordi, FULL));
    }
    BENCH_END;
    BENCH_BEGIN("linear verifier")
    {
        BENCH_ADD(lin_verifier_dkg2(y_Si, t_Si, u, com_Si, codewordi, W, key, FULL));
    }
    BENCH_END;
}

int main(int argc, char *argv[])
{

    printf("\n** Running tests for linearity proofs:\n\n");
    bench_lin();
    bench_lin_dkg2(16);
    bench_lin_dkg2(32);
    bench_lin_dkg2(64);
    bench_lin_dkg3(15);
    bench_lin_dkg3(31);
    bench_lin_dkg3(63);
}
#endif