from sage.all import *
from dualcode_params import mod, USERS, THRESHOLD
from hashlib import sha256


class DotDict(dict):
    """
    A dictionary that supports dot access:
    d.key is equivalent to d['key']
    """

    def __getattr__(self, name):
        if name in self:
            return self[name]
        raise AttributeError(f"'DotDict' object has no attribute '{name}'")

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError(f"'DotDict' object has no attribute '{name}'")

    def __repr__(self):
        return f"DotDict({super().__repr__()})"

    def recursive(self):
        """
        Recursively convert all nested dicts to DotDicts.
        """
        for key, value in self.items():
            if isinstance(value, dict) and not isinstance(value, DotDict):
                self[key] = DotDict(value).recursive()
        return self


class PublicLedger:
    def __init__(self):
        self.public_keys = DotDict()
        self.shares = DotDict()
        self.commitments = DotDict()
        self.decrypted_shares = DotDict()

    def post_pk(self, i, pki):
        self.public_keys[i] = pki

    def post_share(self, dealer_j, i, shat_i, vi, proof):
        self.shares.setdefault(dealer_j, DotDict())
        self.shares[dealer_j][i] = DotDict(shat=shat_i, v=vi, proof=proof)

    def post_commitment(self, dealer_j, commitment):
        self.commitments[dealer_j] = commitment

    def post_decrypted_share(self, dealer_j, i, stilde_i, proof):
        self.decrypted_shares.setdefault(dealer_j, DotDict())
        self.decrypted_shares[dealer_j][i] = DotDict(stilde=stilde_i, proof=proof)

    def get_shares(self, dealer_j):
        return self.shares.get(dealer_j, DotDict())

    def get_commitment(self, dealer_j):
        return self.commitments.get(dealer_j, None)

    def get_decrypted_shares(self, dealer_j):
        return self.decrypted_shares.get(dealer_j, DotDict())

    def get_public_key(self, i):
        return self.public_keys.get(i, None)

    def all_public_keys(self):
        return dict(self.public_keys)


class ReedSolomon:
    def __init__(self, q, n, t):
        assert n < q, "Reed-Solomon code requires n < q"
        self.q = q
        self.n = n
        self.t = t
        self.field = GF(q)
        self.P = PolynomialRing(self.field, "x")
        self.x = self.P.gen()
        self.eval_points = [self.field(i + 1) for i in range(n)]

        # Precompute dual code coeffs
        self.dual_coeffs = []
        for i in range(n):
            xi = self.eval_points[i]
            others = [self.eval_points[j] for j in range(n) if j != i]
            vi = prod([1 / (xi - xj) for xj in others])
            self.dual_coeffs.append(vi)

    def share(self, msg):
        F = self.field
        coeffs = [F.random_element() for _ in range(self.t - 1)]
        poly = F(msg) + sum([coeffs[i] * self.x ** (i + 1) for i in range(self.t - 1)])
        shares = [poly(xi) for xi in self.eval_points]
        return vector(shares)

    def sample_dual_codeword(self):
        deg = self.n - self.t - 1
        f = sum([self.field.random_element() * self.x**i for i in range(deg + 1)])
        return vector(
            [self.dual_coeffs[i] * f(self.eval_points[i]) for i in range(self.n)]
        )

    def verify(self, v):
        assert len(v) == self.n, "Vector length must match code length n"
        dual = self.sample_dual_codeword()
        return v.inner_product(dual) == 0

    def reconstruct(self, points):
        assert len(points) == self.t, "Reconstruct requires exactly t points"
        F = self.field
        x_vals = [F(x) for (x, _) in points]
        y_vals = [F(y) for (_, y) in points]

        def lagrange_basis(j):
            xj = x_vals[j]
            num = prod([0 - x_vals[m] for m in range(self.t) if m != j])
            den = prod([xj - x_vals[m] for m in range(self.t) if m != j])
            return num / den

        result = sum([y_vals[j] * lagrange_basis(j) for j in range(self.t)])
        return result


class PVSS_DDH_System:
    def __init__(self, n, q, ledger: PublicLedger, rs: ReedSolomon):
        self.n = n
        self.q = q
        self.F = GF(self.q)
        self.g = self.F.multiplicative_generator()
        self.h = self._random_generator(self.g)
        self.ledger = ledger
        self.rs = rs
        self.secret_keys = []

    def _random_generator(self, exclude):
        while True:
            h = self.F.random_element()
            if h == 0 or h == exclude:
                continue
            if h.multiplicative_order() == self.F.order() - 1:
                return h

    def generate_keys(self):
        self.secret_keys = [self.F.random_element() for _ in range(self.n)]
        for i, sk in enumerate(self.secret_keys):
            pk = power_mod(self.h, sk, self.q)
            self.ledger.post_pk(i, pk)

    def get_keypair(self, i):
        pk = self.ledger.get_public_key(i)
        sk = self.secret_keys[i]
        return sk, pk

    def distribute(self):
        s = self.F.random_element()
        S = power_mod(self.h, s, self.q)
        shares = self.rs.share(s)
        commitments = []
        encrypted_shares = []
        a1_list = []
        a2_list = []
        w_list = []
        for i in range(self.n):
            pki = self.ledger.get_public_key(i)
            if pki is None:
                raise ValueError(f"Public key for user {i} not found in ledger.")
            shat_i = power_mod(pki, shares[i], self.q)
            v_i = power_mod(self.g, shares[i], self.q)
            w_i = self.F.random_element()
            a1_i = power_mod(self.g, w_i, self.q)
            a2_i = power_mod(pki, w_i, self.q)
            commitments.append(v_i)
            encrypted_shares.append(shat_i)
            a1_list.append(a1_i)
            a2_list.append(a2_i)
            w_list.append(w_i)

        hash_input = "_".join(
            [
                str(x)
                for x in [self.g, self.h]
                + commitments
                + encrypted_shares
                + a1_list
                + a2_list
            ]
        ).encode()
        e = self.F(int(sha256(hash_input).hexdigest(), 16) % self.q)
        for i in range(self.n):
            z_i = self.F(w_list[i] - shares[i] * e)
            proof = (e, z_i)
            self.ledger.post_share(
                dealer_j="D",
                i=i,
                shat_i=encrypted_shares[i],
                vi=commitments[i],
                proof=proof,
            )
        return s  # Return for testing

    def verify(self, dealer_j="D"):
        shares = self.ledger.get_shares(dealer_j)
        if len(shares) != self.n:
            return False

        commitments = []
        encrypted_shares = []
        a1_list = []
        a2_list = []
        for i in range(self.n):
            share = shares.get(i)
            if share is None:
                print("Share is None")
                return False
            v_i = share.v
            shat_i = share.shat
            e, z_i = share.proof
            pki = self.ledger.get_public_key(i)
            if pki is None:
                print("pki is None")
                return False
            a1_i = (power_mod(self.g, z_i, self.q) * power_mod(v_i, e, self.q)) % self.q
            a2_i = (power_mod(pki, z_i, self.q) * power_mod(shat_i, e, self.q)) % self.q
            commitments.append(v_i)
            encrypted_shares.append(shat_i)
            a1_list.append(a1_i)
            a2_list.append(a2_i)

        hash_input = "_".join(
            str(x)
            for x in [self.g, self.h]
            + commitments
            + encrypted_shares
            + a1_list
            + a2_list
        ).encode()
        e_computed = self.F(int(sha256(hash_input).hexdigest(), 16) % self.q)

        for i in range(self.n):
            e, _ = shares[i].proof
            if e_computed != e:
                print(f"Challenge does not match for proof {i}")
                return False

        v = vector(self.F, [shares[i].v for i in range(self.n)])
        return self.rs.verify(v)

    def reconstruct(self, party_indices, dealer_j="D"):
        if len(party_indices) < self.rs.t:
            raise ValueError(f"Need at least {self.rs.t} shares to reconstruct")

        decrypted_shares = self.ledger.get_decrypted_shares(dealer_j)
        points = []
        for i in party_indices:
            if i >= self.n or decrypted_shares[i] is None:
                raise ValueError(f"Invalid or missing decrypted share for party {i}")
            stilde_i = decrypted_shares[i].stilde
            # TODO: Verify DLEQ proof for DLEQ(h, pk_i, stilde_i, shat_i)
            points.append((i + 1, stilde_i))  # Use 1-based indices for ReedSolomon

        # Reconstruct h^s using Lagrange interpolation
        return self.rs.reconstruct(points)


def main():
    ledger = PublicLedger()
    rs = ReedSolomon(q=mod, n=USERS, t=THRESHOLD)
    pvss = PVSS_DDH_System(USERS, mod, ledger, rs)
    pvss.generate_keys()

    s = pvss.distribute()
    print(f"Secret s: {s}, h^s: {power_mod(pvss.h, s, mod)}")

    assert pvss.verify(), "Verification Failed"

    # for i in range(THRESHOLD):
    #     sk_i, pk_i = pvss.get_keypair(i)
    #     shat_i = pvss.ledger.get_shares("D")[i].shat
    #     stilde_i = power_mod(shat_i, sk_i ** (-1), mod)  # stilde_i = h^s_i
    #     # TODO: Generate DLEQ proof for DLEQ(h, pk_i, stilde_i, shat_i)
    #     proof = None
    #     pvss.ledger.post_decrypted_share("D", i, stilde_i, proof)

    # # Reconstruct secret
    # recovered = pvss.reconstruct(range(THRESHOLD))
    # print(f"Recovered h^s: {recovered}")
    # assert recovered == power_mod(pvss.h, s, mod), "Reconstruction failed"


if __name__ == "__main__":
    main()
