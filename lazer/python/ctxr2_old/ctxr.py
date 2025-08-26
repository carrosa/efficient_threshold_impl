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

def main():
    # === Setup ===
    from ctxr_params import deg, m, q, qhat, lhat, beta, sigma_tdec, dim, mod

    seed = b"\0" * 32  # public randomness (proof seed)
    seed_s = b"\0s"
    seed_m = b"\0m"
    prover = lin_prover_state_t(seed, lib.get_params("param"))
    verifier = lin_verifier_state_t(seed, lib.get_params("param"))

    log2o = math.log2(sigma_tdec / 1.55)
    Rqhat = polyring_t(deg, qhat)

    A_right = polymat_t.identity(Rqhat, m) * q
                
    uvec = polyvec_t.urandom_bnd_static(Rqhat, lhat, 0, qhat, seed, 0)
    # uvec_beta = polyvec_t(Rqhat, lhat, [uvec * beta])

    # # Verify uvec and uvec_beta:
    # for i, poly in enumerate(uvec.to_pol_list()):
    #     assert uvec_beta.to_pol_list()[i] == poly * beta
        
    # # uvec_decomposed_1, uvec_decomposed_2 = decompose_vector(uvec, beta, Rqhat)
    # uvec_full = polyvec_t(Rqhat, lhat)
    # for i in range(lhat):
    #     if i % 2 == 0:
    #         uvec_full.set_elem(uvec[i // 2].copy(), i)
    #     else:
    #         uvec_full.set_elem(uvec_beta[i // 2].copy(), i)
    
    # # Verify uvec_full:
    # for i, poly in enumerate(uvec_decomposed_1):
    #     assert uvec_full.to_pol_list()[i * 2] == poly
    # for i, poly in enumerate(uvec_decomposed_2):
    #     assert uvec_full.to_pol_list()[i * 2 + 1] == poly


    # A_left = polymat_t(Rqhat, m, lhat*m)
    # for row in range(m):
    #     tmp_row = polyvec_t(Rqhat, lhat * m)
    #     for col in range(lhat):
    #         tmp_row.set_elem(uvec.get_elem(col), row * lhat + col)
    #     A_left.set_row(row, tmp_row)
    A_left = polymat_t.urandom_static(Rqhat, m, lhat*m, q, seed_s, 0)

    
    
    # A_left = polymat_t(Rp, m, 2 * lhat * m)
    # for coli in range(m):
    #     tmp_row = polyvec_t(Rp, 2 * lhat * m)
    #     for i, elem in enumerate(uvec_full_sub.to_pol_list()):
    #         tmp_row.set_elem(elem, lhat * coli + i)
    #     A_left.set_row(coli, tmp_row)

    A = polymat_t(Rqhat, m, lhat * m + m, [A_left, A_right])


    # s_og = polyvec_t.urandom_bnd_static(Rqhat, lhat * m, 0, q, seed_s, 0)
    # print("s_og.l2sqr =", s_og.l2sqr())
    # print("s_og.linf =", s_og.linf())
    # s_og0, s_og1 = decompose_vector(s_og, beta, Rqhat)
    # s_og0 = polyvec_t(Rqhat, lhat * m, s_og0)
    # s_og1 = polyvec_t(Rqhat, lhat * m, s_og1)
    # s_og0.print()
    # s_og1.print()
    # print(beta)
    # s_top = polyvec_t(Rqhat, lhat*m)
    # for i in range(lhat*m):
    #     if i % 2 == 0:
    #         s_top.set_elem(s_og0[i // 2].copy(), i)
    #     else:
    #         s_top.set_elem(s_og1[i // 2].copy(), i)

    # print("s_top.l2sqr =", s_top.l2sqr())
    # print("s_top.linf =", s_top.linf())
    binary_poly = poly_t(Rqhat, [1] + [0]*(deg-1)) #
    s_top = [binary_poly for _ in range(lhat*m)]
    s_top = polyvec_t(Rqhat, lhat * m, s_top)

    
    s_bot = polyvec_t.grandom_static(Rqhat, m, int(log2o), seed_m, 0, 0) # Ei
    print("s_bot coeffs max abs:", max(abs(c) for p in s_bot.to_pol_list() for c in p.to_list()))

    s = polyvec_t(Rqhat, lhat * m + m, s_top.to_pol_list() + s_bot.to_pol_list())

    
    print("slen", len(s.to_pol_list()))
    t = A * s

    print("s_bot.l2sqr =", s_bot.l2sqr())
    print("s_bot.linf =", s_bot.linf())

    # exit()

    # === Prover setup ===
    print("set_statement...")
    prover.set_statement(A, -t)

    print("set_witness...")
    print(type(s))

    prover.set_witness(w=s)

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


if __name__ == "__main__":
    main()
