from math import sqrt, exp
from sage.all import var, solve

vname = "param"

q = 2101003516513793
qhat = 633544653872304603170707504554316801
# q = 12289
# qhat = 786433

lhat = 2
m = 2

deg   = 256  # ring Rp degree d
mod   = qhat    # ring Rp modulus p
# mod = qhat
dim   = (m, lhat*m + m)  # dimensions of A in Rp^(m,n)

sigma_tdec = 1208

rej_rate = 0.01
tail_bound_2 = 1.01
while ((tail_bound_2**deg) * exp((deg / 2) * (1 - tail_bound_2**2))) > rej_rate:
    tail_bound_2 += 0.01

# Partition of s
wpart = [
    list(range(lhat*m)),
    list(range(lhat*m, lhat*m + m))
]

B = q * sqrt(deg)
x = var('x')
eq = x**2 *sqrt(deg)+sqrt(deg)*x - B == 0
solutions = solve(eq, x)
beta = int(solutions[0].rhs().n() if solutions[0].rhs().n() > 0 else solutions[1].rhs().n())

# l2-norm bounds
wl2 = [
    0,
    tail_bound_2 * sigma_tdec * sqrt(deg)
]


# binary coefficients (0 = no restriction)
wbin = [1, 0]

# rejection sampling flags (0 = not used)
wrej = [0, 0]

# Optional: linf-norm bound on s (commented out)
# wlinf = 1

