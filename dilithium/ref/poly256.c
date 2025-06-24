#include "poly256.h"

/* Return 1 if the two polynomials are equal (all coefficients match), or 0 otherwise. */
int polycmp(const poly *p1, const poly *p2)
{
    int i;
    for (i = 0; i < N; i++)
    {
        if (mpz_cmp(p1->coeffs[i], p2->coeffs[i]) != 0)
            return 0; // not equal if any coefficient differs
    }
    return 1;
}

// Init polynomial
void poly_init(poly *p)
{
    for (int i = 0; i < N; i++)
    {
        mpz_init2(p->coeffs[i], 128);
    }
}
void poly_zero(poly *p)
{
    unsigned int i;
    for (i = 0; i < N; i++)
    {
        mpz_set_ui(p->coeffs[i], 0);
    }
}

void poly_1d_init(poly *arr, size_t count)
{
    for (size_t i = 0; i < count; i++)
    {
        poly_init(&arr[i]);
    }
}

void poly_inits(poly *p, ...)
{
    va_list args;
    poly *f;
    va_start(args, p);
    f = p;
    while (f != NULL)
    {
        poly_init(f);
        f = va_arg(args, poly *);
    }
    va_end(args);
}

// Clear polynomial
void poly_clear(poly *p)
{
    unsigned int i;
    for (i = 0; i < N; i++)
    {
        mpz_clear(p->coeffs[i]);
    }
}
void poly_1d_clear(poly *arr, size_t count)
{
    for (size_t i = 0; i < count; i++)
    {
        poly_clear(&arr[i]);
    }
}
void poly_clears(poly *p, ...)
{
    va_list args;
    poly *f;
    va_start(args, p);
    f = p;
    while (f != NULL)
    {
        poly_clear(f);
        f = va_arg(args, poly *);
    }
    va_end(args);
}

void poly_reduce(poly *a)
{
    mpz_t t, tmp;
    mpz_inits(t, tmp, NULL);

    for (int i = 0; i < N; i++)
    {
        // t = floor(a * m / 2^(2*log2(Q)))
        mpz_mul(t, a->coeffs[i], barrett_m);
        mpz_fdiv_q_2exp(t, t, barrett_shift);
        // a = a - t * Q
        mpz_mul(tmp, t, GMP_Q);
        mpz_sub(a->coeffs[i], a->coeffs[i], tmp);
        // Center the result
        if (mpz_cmp(a->coeffs[i], q_half) >= 0)
        {
            mpz_sub(a->coeffs[i], a->coeffs[i], GMP_Q);
        }
        else if (mpz_cmp(a->coeffs[i], neg_q_half) < 0)
        {
            mpz_add(a->coeffs[i], a->coeffs[i], GMP_Q);
        }
    }    
    mpz_clears(t, tmp, NULL);
}

// Reduce polynomial
void poly_reduce_old(poly *a)
{
    unsigned int i;
    // mpz_t tmp;
    // mpz_init(tmp);
    for (i = 0; i < N; i++)
    {
        /* Use the GMP freeze routine (in reduce256.c/h) to reduce modulo Q */
        mpz_mod(a->coeffs[i], a->coeffs[i], GMP_Q);
        // freeze(tmp, a->coeffs[i]);
        // mpz_set(a->coeffs[i], tmp);
    }
    // mpz_clear(tmp);
}

// Add Q to negative coefficients i a polynomial
void poly_caddq(poly *a)
{
    unsigned int i;
    mpz_t tmp;
    mpz_init(tmp);
    for (i = 0; i < N; i++)
    {
        caddq(tmp, a->coeffs[i]);
        mpz_set(a->coeffs[i], tmp);
    }
    mpz_clear(tmp);
}

// Polynomial arithmetic
void poly_add(poly *c, const poly *a, const poly *b)
{
    unsigned int i;
    for (i = 0; i < N; i++)
    {
        mpz_add(c->coeffs[i], a->coeffs[i], b->coeffs[i]);
    }
}

void poly_sub(poly *c, const poly *a, const poly *b)
{
    unsigned int i;
    for (i = 0; i < N; i++)
    {
        mpz_sub(c->coeffs[i], a->coeffs[i], b->coeffs[i]);
    }
}

void poly_shiftl(poly *a)
{
    unsigned int i;
    for (i = 0; i < N; i++)
    {
        mpz_mul_2exp(a->coeffs[i], a->coeffs[i], D);
    }
}

// NTT
void poly_ntt(poly *a)
{
    ntt(a->coeffs);
}

// Inverse NTT
void poly_invntt_tomont(poly *a)
{
    invntt_tomont(a->coeffs);
}

// Pointwise multiplication in NTT Domain
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b)
{
    unsigned int i;
    mpz_t product;
    mpz_init(product);
    for (i = 0; i < N; i++)
    {
        mpz_mul(product, a->coeffs[i], b->coeffs[i]);
        montgomery_reduce(c->coeffs[i], product);
    }
    mpz_clear(product);
}

// Norm check
int poly_chknorm(const poly *a, int32_t B)
{
    unsigned int i;
    for (i = 0; i < N; i++)
    {
        if (mpz_cmpabs_ui(a->coeffs[i], (unsigned long)B) >= 0)
            return 1;
    }
    return 0;
}

/*============================================================
 * Uniform and Noise Sampling
 *
 * For a 256-bit modulus Q, we generate 32 bytes (256 bits) per coefficient.
 * The following helper performs rejection sampling: It reads 32-byte blocks,
 * converts them to an mpz_t candidate, and accepts the candidate if it is less than Q.
 *============================================================*/
static unsigned int rej_uniform_gmp(mpz_t *a, unsigned int len,
                                    const uint8_t *buf, unsigned int buflen)
{
    unsigned int ctr = 0, pos = 0;
    mpz_t candidate;
    mpz_init(candidate);
    while (ctr < len && pos + SAMPLEBYTES <= buflen)
    {
        /* Import 32 bytes (big-endian) as an mpz_t number. */
        mpz_import(candidate, SAMPLEBYTES, 1, 1, 0, 0, buf + pos);
        pos += SAMPLEBYTES;
        if (mpz_cmp(candidate, GMP_Q) < 0)
        {
            mpz_set(a[ctr], candidate);
            ctr++;
        }
    }
    mpz_clear(candidate);
    return ctr;
}
#define POLY_UNIFORM_NBLOCKS ((768 + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES)
void poly_uniform(poly *a, const uint8_t seed[SEEDBYTES], uint16_t nonce)
{
    unsigned int i, ctr, off;
    unsigned int buflen = POLY_UNIFORM_NBLOCKS * STREAM256_BLOCKBYTES;
    /* Allocate a buffer large enough to hold several blocks plus extra bytes */
    uint8_t buf[POLY_UNIFORM_NBLOCKS * STREAM256_BLOCKBYTES + 3];
    stream256_state state;
    stream256_init(&state, seed, nonce);
    stream256_squeezeblocks(buf, POLY_UNIFORM_NBLOCKS, &state);
    ctr = rej_uniform_gmp(a->coeffs, N, buf, buflen);
    while (ctr < N)
    {
        off = buflen % 3;
        for (i = 0; i < off; i++)
            buf[i] = buf[buflen - off + i];
        stream256_squeezeblocks(buf + off, 1, &state);
        buflen = STREAM256_BLOCKBYTES + off;
        ctr += rej_uniform_gmp(a->coeffs + ctr, N - ctr, buf, buflen);
    }
}

#if defined(ETA)
/* Sampling coefficients in [-ETA, ETA] via rejection.
   We first sample into a temporary int32_t array and then set the GMP coefficient.
*/
#if ETA == 2
#define POLY_UNIFORM_ETA_NBLOCKS ((136 + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES)
#elif ETA == 4
#define POLY_UNIFORM_ETA_NBLOCKS ((227 + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES)
#endif

static unsigned int rej_eta(int32_t *a, unsigned int len,
                            const uint8_t *buf, unsigned int buflen)
{
    unsigned int ctr = 0, pos = 0;
    uint32_t t0, t1;
    while (ctr < len && pos < buflen)
    {
        t0 = buf[pos] & 0x0F;
        t1 = buf[pos++] >> 4;
#if ETA == 2
        if (t0 < 15)
        {
            t0 = t0 - ((205 * t0) >> 10) * 5;
            a[ctr++] = 2 - t0;
        }
        if (t1 < 15 && ctr < len)
        {
            t1 = t1 - ((205 * t1) >> 10) * 5;
            a[ctr++] = 2 - t1;
        }
#elif ETA == 4
        if (t0 < 9)
            a[ctr++] = 4 - t0;
        if (t1 < 9 && ctr < len)
            a[ctr++] = 4 - t1;
#endif
    }
    return ctr;
}

void poly_uniform_eta(poly *a, const uint8_t seed[CRHBYTES], uint16_t nonce)
{
    unsigned int ctr;
    unsigned int buflen = POLY_UNIFORM_ETA_NBLOCKS * STREAM256_BLOCKBYTES;
    uint8_t buf[POLY_UNIFORM_ETA_NBLOCKS * STREAM256_BLOCKBYTES];
    stream256_state state;
    stream256_init(&state, seed, nonce);
    stream256_squeezeblocks(buf, POLY_UNIFORM_ETA_NBLOCKS, &state);
    int32_t tmp[N];
    ctr = rej_eta(tmp, N, buf, buflen);
    while (ctr < N)
    {
        stream256_squeezeblocks(buf, 1, &state);
        ctr += rej_eta(tmp + ctr, N - ctr, buf, STREAM256_BLOCKBYTES);
    }
    for (unsigned int i = 0; i < N; i++)
    {
        mpz_set_si(a->coeffs[i], tmp[i]);
    }
}

void poly_copy(poly *dst, const poly *src)
{
    for (int i = 0; i < N; i++)
    {
        mpz_set(dst->coeffs[i], src->coeffs[i]);
    }
}

void fast_polyvec_copy(poly *dst, const poly *src, size_t length)
{
    for (size_t i = 0; i < length; i++)
    {
        for (unsigned int j = 0; j < N; j++)
        {
            mpz_t *dst_coeff = &dst[i].coeffs[j];
            const mpz_t *src_coeff = &src[i].coeffs[j];
            mp_size_t size = (*src_coeff)->_mp_size;
            mp_ptr dst_limbs = (*dst_coeff)->_mp_d;
            mp_srcptr src_limbs = (*src_coeff)->_mp_d;
            mp_size_t abs_size = __GMP_ABS(size);

            if (abs_size > (*dst_coeff)->_mp_alloc)
            {
                dst_limbs = mpz_realloc(*dst_coeff, abs_size);
                (*dst_coeff)->_mp_d = dst_limbs;
            }
            if (abs_size > 0)
            {
                mpn_copyi(dst_limbs, src_limbs, abs_size);
            }
            (*dst_coeff)->_mp_size = size;
        }
    }
}

void poly_scale_si(poly *r, int32_t c)
{
    for (int j = 0; j < N; j++)
    {
        // Multiply the j-th coefficient by c
        mpz_mul_si(r->coeffs[j], r->coeffs[j], c);
    }

    poly_reduce(r);
}

void poly_scale(poly *r, mpz_t c)
{
    for (int j = 0; j < N; j++)
    {
        // Multiply the j-th coefficient by c
        mpz_mul(r->coeffs[j], r->coeffs[j], c);
    }

    poly_reduce(r);
}

void poly_mod(poly *r, mpz_t c)
{
    unsigned int j;
    for (j = 0; j < N; j++)
    {
        mpz_mod(r->coeffs[j], r->coeffs[j], c);
    }
}

int poly_equal(poly *a, poly *b)
{
    for (int i = 0; i < N; i++)
    {
        if (mpz_cmp(a->coeffs[i], b->coeffs[i]) != 0)
        {
            return 0;
        }
    }
    return 1;
}

/*************************************************
 * Name:        challenge
 *
 * Description: Implementation of H. Samples polynomial with TAU nonzero
 *              coefficients in {-1,1} using the output stream of
 *              SHAKE256(seed).
 *
 * Arguments:   - poly *c: pointer to output polynomial
 *              - const uint8_t mu[]: byte array containing seed of length CTILDEBYTES
 **************************************************/
void poly_challenge(poly *c, const uint8_t seed[CTILDEBYTES])
{
    unsigned int i, b, pos;
    uint64_t signs;
    uint8_t buf[SHAKE256_RATE];
    keccak_state state;

    shake256_init(&state);
    shake256_absorb(&state, seed, CTILDEBYTES);
    shake256_finalize(&state);
    shake256_squeezeblocks(buf, 1, &state);

    signs = 0;
    for (i = 0; i < 8; ++i)
        signs |= (uint64_t)buf[i] << 8 * i;
    pos = 8;

    for (i = N - TAU; i < N; ++i)
    {
        do
        {
            if (pos >= SHAKE256_RATE)
            {
                shake256_squeezeblocks(buf, 1, &state);
                pos = 0;
            }

            b = buf[pos++];
        } while (b > i);

        mpz_set(c->coeffs[i], c->coeffs[b]);
        mpz_set_si(c->coeffs[b], 1 - 2 * (signs & 1));
        signs >>= 1;
    }
}

void poly_hash_shake256(
    poly *p,
    uint8_t *hash,
    size_t outlen)
{
    keccak_state state;
    shake256_init(&state);
    // Compute space needed to be allocated
    size_t coeff_bytes = (mpz_sizeinbase(GMP_Q, 2) + 7) / 8; // Ceiling of bits to bytes
    size_t count = 0;

    for (int i = 0; i < N; i++)
    {
        uint8_t tmp[coeff_bytes];
        memset(tmp, 0, coeff_bytes);
        size_t *data = mpz_export(tmp, NULL, 1, coeff_bytes, 1, 0, p->coeffs[i]);

        int sign = mpz_sgn(p->coeffs[i]);
        uint8_t sign_byte = (sign < 0) ? 1 : 0;
        // shake256_absorb(&state, data, 1);
        shake256_absorb(&state, &sign_byte, 1);
        if (count > 0)
        {
            shake256_absorb(&state, (const uint8_t *)data, count * sizeof(mp_limb_t));
            free(data);
        }
    }
    shake256_finalize(&state);
    shake256_squeeze(hash, outlen, &state);
}

void poly_arr_hash_shake256(
    const poly *arr,
    size_t len,
    uint8_t *hash,
    size_t outlen)
{
    keccak_state state;
    shake256_init(&state);
    for (int i = 0; i < len; i++)
    {
        uint8_t tmp_hash[SHAKE256_RATE];
        poly_hash_shake256(&arr[i], tmp_hash, SHAKE256_RATE);
        shake256_absorb(&state, tmp_hash, SHAKE256_RATE);
    }
    shake256_finalize(&state);
    shake256_squeeze(hash, outlen, &state);
}

void poly_matrix_hash_shake256(
    const poly *mat,
    size_t rows,
    size_t cols,
    uint8_t *hash,
    size_t outlen)
{
    keccak_state state;
    shake256_init(&state);

    for (size_t i = 0; i < rows * cols; i++)
    {
        uint8_t temp_hash[32];
        poly_hash_shake256(&mat[i], temp_hash, 32);
        shake256_absorb(&state, temp_hash, 32);
    }

    shake256_finalize(&state);
    shake256_squeeze(hash, outlen, &state);
}

void H1(
    poly arr[K],
    uint8_t *hash)
{
    keccak_state state;
    shake256_init(&state);

    size_t coeff_bytes = (mpz_sizeinbase(GMP_Q, 2) + 7) / 8; // Ceiling of bits to bytes
    size_t poly_bytes = N * coeff_bytes;
    size_t total_len = K * poly_bytes;

    uint8_t *buf = malloc(total_len);
    size_t pos = 0;

    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            mpz_export(buf + pos + j * coeff_bytes, NULL, 1, coeff_bytes, 1, 0, arr[i].coeffs[j]);
        }
        pos += poly_bytes;
    }
    shake256_absorb(&state, buf, total_len);
    shake256_finalize(&state);
    shake256_squeeze(hash, 32, &state);
    free(buf);
}

void H0(
    poly arr[K],
    uint8_t *hash)
{
    keccak_state state;
    shake256_init(&state);

    size_t coeff_bytes = (mpz_sizeinbase(GMP_Q, 2) + 7) / 8; // Ceiling of bits to bytes
    size_t poly_bytes = N * coeff_bytes;
    size_t total_len = K * poly_bytes;

    uint8_t *buf = malloc(total_len);
    size_t pos = 0;

    for (size_t i = 0; i < K; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            mpz_export(buf + pos + j * coeff_bytes, NULL, 1, coeff_bytes, 1, 0, arr[i].coeffs[j]);
        }
        pos += poly_bytes;
    }
    shake256_absorb(&state, buf, total_len);
    shake256_finalize(&state);
    shake256_squeeze(hash, 32, &state);
    free(buf);
}

void H0_matrix(poly *arr, size_t rows, size_t cols, uint8_t *hash)
{
    keccak_state state;
    shake256_init(&state);

    size_t coeff_bytes = (mpz_sizeinbase(GMP_Q, 2) + 7) / 8;
    size_t total_len = rows * cols * N * coeff_bytes;

    uint8_t *buf = malloc(total_len);
    if (!buf)
        abort();

    for (size_t i = 0; i < rows * cols; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            size_t offset = (i * N + j) * coeff_bytes;
            mpz_export(buf + offset, NULL, 1, coeff_bytes, 1, 0, arr[i].coeffs[j]);
        }
    }

    shake256_absorb(&state, buf, total_len);
    shake256_finalize(&state);
    shake256_squeeze(hash, 32, &state);

    free(buf);
}

void get_random_dummy_poly(poly *p)
{
    gaussian_sampler_w(p->coeffs, N); // Fill with Gaussian samples
    poly_reduce(p);
}

#endif
