from math import sqrt, exp
# from sage.all import var, solve

vname = "param"

q = 2101003516513793
qhat = 633544653872304603170707504554316801
# q = 12289
# qhat = 786433

lhat = 22
khat = 22
m = 22

deg   = 256  # ring Rp degree d
mod   = qhat    # ring Rp modulus p
dim   = (lhat + m, khat + lhat + 2*m)  # dimensions of A in Rp^(m,n)
#dim = (lhat, khat)

sigma_ctx = 1208

rej_rate = 0.01
tail_bound_2 = 1.01
while ((tail_bound_2**deg) * exp((deg / 2) * (1 - tail_bound_2**2))) > rej_rate:
    tail_bound_2 += 0.01

# Partition of s
# wpart = (
#     [[i] for i in range(khat)] +
#     [[i] for i in range(khat, khat+ lhat)] + 
#     [[i] for i in range(khat + lhat, khat + lhat + m)] +
#     [[i] for i in range(khat + lhat + m, khat + lhat + 2*m)]
# )
wpart = [
    list(range(khat + lhat + 2*m))
]

# B = q * sqrt(deg)
# x = var('x')
# eq = x**2 *sqrt(deg)+sqrt(deg)*x - B == 0
# solutions = solve(eq, x)
# beta = int(solutions[0].rhs().n() if solutions[0].rhs().n() > 0 else solutions[1].rhs().n())

# l2-norm bounds
# wl2 = (
#     [3*tail_bound_2 * sigma_ctx * sqrt(deg) for _ in range(khat)] +
#     [3*tail_bound_2 * sigma_ctx * sqrt(deg) for _ in range(khat, khat+ lhat)] + 
#     [3*tail_bound_2 * sigma_ctx * sqrt(deg)  for _ in range(khat + lhat, khat + lhat + m)] +
#     [3*q * sqrt(deg)  for _ in range(khat + lhat + m, khat + lhat + 2*m)]
# )

wl2 = [
    3*tail_bound_2*sigma_ctx*sqrt(deg)
]


# binary coefficients (0 = no restriction)
# wbin = [0 for _ in range(khat + lhat + 2*m)]
wbin = [0]

# rejection sampling flags (0 = not used)
# wrej = [0 for _ in range(khat + lhat + 2*m)]
wrej = [0]


# Optional: linf-norm bound on s (commented out)
# wlinf = 1

