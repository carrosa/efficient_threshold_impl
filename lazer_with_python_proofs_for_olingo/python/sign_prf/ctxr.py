import sys
import time
sys.path.append('..')   # path to LaZer module
from lazer import *      # import LaZer Python module
import hashlib           # for SHAKE128
from lazer import lin_prover_state_t, lin_verifier_state_t  # Explicit import
from _ctxr_params_cffi import lib


def decompose_vector(m: polyvec_t, B: int, ring: polyring_t):
    m0 = []
    for pol in m.to_pol_list():
        tmp = []
        for coeff in pol.to_list():
           tmp.append(coeff % B) 
        m0.append(poly_t(ring, tmp))
    m1 = []
    for pol in m.to_pol_list():
        tmp = []
        for coeff in pol.to_list():
           tmp.append(coeff // B) 
        m1.append(poly_t(ring, tmp))
    return m0, m1

def main():
    # === Setup ===
    from ctxr_params import mod, deg, dim, moddimEnc, secdimEnc, mhat, q_small, secdim, moddim, sigma_w, tail_bound_2, beta

    d, p, m, n, khat, lhat, mhat, q = deg, mod, dim[0], dim[1], moddimEnc, secdimEnc, mhat, q_small

    seed = b'\0' * 32  # public randomness (proof seed)
    prover = lin_prover_state_t(seed, lib.get_params("param"))
    verifier = lin_verifier_state_t(seed, lib.get_params("param"))

    Rp = polyring_t(d, p)


    # === Matrix A = [A1 | A2 | A3] ===
    dom_A1 = 0  # domain separator

    # A1 ∈ Rp^{m × khat} : uniformly random
    A1 = polymat_t.urandom_static(Rp, m, khat, p, seed, dom_A1)


    # A2 ∈ Rp^{m × lhat} : stacked [q*I_lhat ; 0]
    I_lhat = polymat_t.identity(Rp, lhat)
    qI_lhat = q * I_lhat
    Z_lhat = polymat_t(Rp, mhat, lhat)

    A2 = polymat_t(Rp, m, lhat)
    for i in range(lhat):
        A2.set_row(i, qI_lhat.get_row(i).copy())
    for i in range(mhat):
        A2.set_row(lhat + i, Z_lhat.get_row(i).copy())

    # A3 ∈ Rp^{m × mhat + 2*secdim} : stacked [0 ; q*I | I] (stacked not correct)
    Z_top = polymat_t(Rp, lhat, mhat + 2*secdim)
    I_mhat = polymat_t.identity(Rp, mhat)
    qI_mhat = q * I_mhat
    m_mat_top = polymat_t.identity(Rp, 2*secdim)
    m_mat_bot = polymat_t(Rp, mhat-2*secdim, 2*secdim)
    m_mat = polymat_t(Rp, mhat, 2*secdim)
    for i in range(2*secdim):
        m_mat.set_row(i, m_mat_top.get_row(i).copy())
    for i in range(mhat - 2*secdim):
        m_mat.set_row(2*secdim + i, m_mat_bot.get_row(i).copy())
    A3_bottom = polymat_t(Rp, mhat, mhat + 2*secdim, [qI_mhat, m_mat])

    A3 = polymat_t(Rp, m, mhat + 2*secdim)
    for i in range(lhat):
        A3.set_row(i, Z_top.get_row(i).copy())
    for i in range(mhat):
        A3.set_row(lhat + i, A3_bottom.get_row(i).copy())

    
    
    # Final matrix A ∈ Rp^{m × n}
    A = polymat_t(Rp, m, n, [A1, A2, A3])

    # === Build s = [r || e1 || e2 || m] ===
    # Derive seeds for reproducible randomness
    seed_r  = seed
    seed_e1 = hashlib.shake_128(seed + b"e1").digest(32)
    seed_e2 = hashlib.shake_128(seed + b"e2").digest(32)
    seed_m  = hashlib.shake_128(seed + b"m").digest(32)

    # r, e1, e2: binomial distribution with η = 2
    r = polyvec_t.brandom_static(Rp, khat, 2, seed_r, 0)
    e1 = polyvec_t.brandom_static(Rp, lhat, 2, seed_e1, 0)
    e2 = polyvec_t.brandom_static(Rp, mhat, 2, seed_e2, 0)

    # m: uniform in [0, q)
    # msg = polyvec_t.urandom_bnd_static(Rp, mhat, 0, q, seed_m, 0)
    log2o = math.log2(sigma_w / 1.55)
    # tail_bound_2*sigma_w*math.sqrt(secdim*deg)
    msg_og = polyvec_t.grandom_static(Rp, secdim, int(log2o), seed_m, 0, 0)
    # msg = polyvec_t.brandom_static(Rp, mhat, 2, seed_m, 0)
    m0_lst, m1_lst = decompose_vector(msg_og, beta, Rp)
    for i, pol in enumerate(msg_og.to_pol_list()):
        for j, coeff in enumerate(pol.to_list()):
            assert(coeff == m0_lst[i][j] + beta*m1_lst[i][j])
    m0 = polyvec_t(Rp, secdim, m0_lst)
    m1 = polyvec_t(Rp, secdim, m1_lst)

    # Full ciphertext vector s ∈ Rp^n (r || e1 || e2 || m)
    # s = polyvec_t(Rp, n, [r, e1, e2, msg])
    s = polyvec_t(Rp, n, [r, e1, e2, m0, m1])
    # A = polymat_t.urandom_static(Rp, m, n, p, seed_e1, 0)
    # s = polyvec_t.urandom_bnd_static(Rp, n, 0, 1, seed, 0)

    t = A*s

    # === Prover setup ===
    print("set_statement...")
    prover.set_statement(A, -t)

    print("set_witness...")
    print(type(s))

    prover.set_witness(w = s)


    print("generate proof ...")
    start_time = time.time()
    proof = prover.prove()
    print("proof done")
    end_time = time.time()

    # Report proof time and size
    elapsed_ms = (end_time - start_time) * 1000
    proof_size_bits = len(proof) * 8
    proof_size_kb = proof_size_bits / 8000  # 1 KB = 8000 bits
    print(f"Prover time: {elapsed_ms:.2f} ms")
    print(f"Proof size: {proof_size_bits} bits ({proof_size_kb:.2f} KB)")

    # === Verifier ===
    verifier.set_statement(A, -t)

    print("verify proof ... ")
    try:
        start_time = time.time()
        verifier.verify(proof)
        end_time = time.time()
        print(f"Verify time: {(end_time - start_time) * 1000:.2f} ms")
    except VerificationError as e:
        print("reject")
        print(e)


if __name__ == "__main__":
    main()
