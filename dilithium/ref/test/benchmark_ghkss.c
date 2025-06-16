#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "bench.h"
#include "../params.h"
#include "ghkss_inputs.h"
#include "../poly256.h"
#include <assert.h>

static inline void keygen_wrapper(void)
{
    keygen_input_t *in = malloc(sizeof(keygen_input_t));
    if (!in)
    {
        fprintf(stderr, "Memory allocation failed for keygen_input_t\n");
        exit(EXIT_FAILURE);
    }
    init_keygen_input(in);
    keygen_e_1d(in->A, in->B, in->S);
    clear_keygen_input(in);
    free(in);
}

static inline void dkg1_wrapper(dkg_input_t *in)
{
    init_dkg1(in);
    dk_gen_1(in->Ae, in->Si, in->Ei, in->Bi, in->h_Bi);
}

static inline void dkg2_wrapper(dkg_input_t *in)
{
    dk_gen_2(in->Si, in->Ei, in->Ae, in->users);
}

static inline void dkg3_wrapper(dkg_input_t *in)
{
    dk_gen_3(in->Bj, in->Sij, in->Be, in->Sei);
}

static inline void bench_dkg(void)
{
    dkg_input_t *in = malloc(sizeof(dkg_input_t));
    init_dkg(in);
    BENCH_ONCE_W_ARGS("DKG round 1", dkg1_wrapper, in);
    BENCH_ONCE_W_ARGS("DKG round 2", dkg2_wrapper, in);
    BENCH_ONCE_W_ARGS("DKG round 3", dkg3_wrapper, in);
    clear_dkg(in);    
    free(in);
}

static inline void sig1_wrapper(sig_input_t *in)
{
    as_sign_1(in->A, in->Ae, in->Be, in->u, in->v);
}

static inline void sig2_wrapper(sig_input_t *in)
{
    as_sign_2(in->Ae, in->Be, in->u, in->v, in->A, in->Ae, in->Be);
}

static inline void sig3_wrapper(sig_input_t *in)
{
}

static inline void bench_sign(void)
{
    sig_input_t *in = malloc(sizeof(sig_input_t));
    init_sig(in);
    BENCH_ONCE_W_ARGS("Sig round 1", sig1_wrapper, in);
    BENCH_ONCE_W_ARGS("Sig round 2", sig2_wrapper, in);
    BENCH_ONCE_W_ARGS("Sig round 3", sig3_wrapper, in);
    clear_sig(in);
    free(in);
}


int main(void)
{
    init_params();
    init_zetas();

    bench_dkg();

    free_params();
    clear_zetas();
    return 0;
}
