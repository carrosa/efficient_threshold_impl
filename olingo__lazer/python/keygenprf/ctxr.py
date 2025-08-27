import sys
import time

sys.path.append("..")  # path to LaZer module
from lazer import *  # import LaZer Python module
import hashlib  # for SHAKE128
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


def check_large_coeffs(label, vec, qhat):
    print(f"Checking {label}...")
    for i, pol in enumerate(vec.to_pol_list()):
        coeffs = pol.to_list()
        for j, c in enumerate(coeffs):
            if not (0 <= c < qhat * 2):
                print(f"[{label}][{i}][{j}] = {c} is too large")
                return (i, j, c)
    print(f"{label}: all coefficients < 2·qhat")
    return None


def prove_verify(A: polymat_t, w: polyvec_t, t: polyvec_t):
    seed = b"\0" * 32  # public randomness (proof seed)
    prover = lin_prover_state_t(seed, lib.get_params("param"))
    verifier = lin_verifier_state_t(seed, lib.get_params("param"))
    # === Prover setup ===
    print("set_statement...")
    prover.set_statement(A, -t)

    print("set_witness...")
    prover.set_witness(w)

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
    except VerificationError as e:
        print("reject")
        print(e)
    else:
        print("accept")


def main():
    # === Setup ===
    from ctxr_params import dim, deg, qhat, lhat, khat, q, m, sigma_ctx

    seed_A11 = hashlib.shake_128(b"\0A11").digest(32)
    seed_A21 = hashlib.shake_128(b"\0A21").digest(32)

    seed_r = hashlib.shake_128(b"\0r").digest(32)
    seed_e1 = hashlib.shake_128(b"\0e1").digest(32)
    seed_e2 = hashlib.shake_128(b"\0e2").digest(32)
    seed_msg = hashlib.shake_128(b"\0msg").digest(32)

    log2o = int(math.log2(sigma_ctx / 1.55))
    Rqhat = polyring_t(deg, qhat)

    # === Generate A ===
    A = polymat_t(Rqhat, dim[0], dim[1])

    A11 = polymat_t.urandom_static(Rqhat, lhat, khat, qhat, seed_A11, 0)
    A12 = polymat_t.identity(Rqhat, lhat) * q
    A13 = polymat_t(Rqhat, m, m)
    A14 = polymat_t(Rqhat, m, m)
    Arow1 = polymat_t(Rqhat, lhat, khat + lhat + 2 * m, [A11, A12, A13, A14])
    for i in range(lhat):
        A.set_row(i, Arow1.get_row(i))

    A21 = polymat_t.urandom_static(Rqhat, m, khat, qhat, seed_A21, 0)
    A22 = polymat_t(Rqhat, lhat, lhat)
    A23 = polymat_t.identity(Rqhat, m) * q
    A24 = polymat_t.identity(Rqhat, m)
    Arow2 = polymat_t(Rqhat, m, khat + lhat + 2 * m, [A21, A22, A23, A24])
    for i in range(m):
        A.set_row(lhat + i, Arow2.get_row(i))

    # === Generate w ===
    # r = polyvec_t.grandom_static(Rqhat, khat, log2o, seed_r, 0, 0)
    # e1 = polyvec_t.grandom_static(Rqhat, lhat, log2o, seed_e1, 0, 0)
    # e2 = polyvec_t.grandom_static(Rqhat, m, log2o, seed_e2, 0, 0)
    r = polyvec_t.brandom_static(Rqhat, khat, 2, seed_r, 0)
    e1 = polyvec_t.brandom_static(Rqhat, lhat, 2, seed_e1, 0)
    e2 = polyvec_t.brandom_static(Rqhat, m, 2, seed_e2, 0)
    # msg = polyvec_t.urandom_bnd_static(Rq, m, 0, q, seed_msg, 0).lift(Rqhat)
    # msg = polyvec_t.urandom_bnd_static(Rqhat, m, 0, 10, seed_msg, 0)
    msg = polyvec_t.grandom_static(Rqhat, m, int(log2o), seed_msg, 0, 0)
    w = polyvec_t(Rqhat, khat + lhat + 2 * m, [r, e1, e2, msg])

    # A = polymat_t.urandom_static(Rqhat, dim[0], dim[1], q, seed_A11, 0)
    # w = polyvec_t.urandom_bnd_static(Rqhat, dim[1], 0, 1, seed_r, 0)
    # A = polymat_t(Rqhat, dim[0], dim[1])
    # w = polyvec_t(Rqhat, dim[1])
    t = A * w

    print("dim: ", dim)
    print("dim(w):", w.dim)
    print("dim(A):", f"({A.rows}, {A.cols})")
    print("dim(t):", t.dim)
    
    # prove_verify(A, w, t)

    seed = b"\0" * 32  # public randomness (proof seed)
    prover = lin_prover_state_t(seed, lib.get_params("param"))
    verifier = lin_verifier_state_t(seed, lib.get_params("param"))
    # === Prover setup ===
    print("set_statement...")
    prover.set_statement(A, -t)

    print("set_witness...")
    prover.set_witness(w)

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

    # # === Verifier ===
    # verifier.set_statement(A, -t)

    # print("verify proof ... ")
    # try:
    #     start_time = time.time()
    #     verifier.verify(proof)
    #     end_time = time.time()
    # except VerificationError as e:
    #     print("reject")
    #     print(e)
    # else:
    #     print("accept")


if __name__ == "__main__":
    main()
