from math import sqrt, exp
from sage.all import var, solve

vname = "param"

q = 2101003516513793
qhat = 633544653872304603170707504554316801
# q = 12289
# qhat = 786433

lhat = 22
m = 22

deg   = 256  # ring Rp degree d
mod   = qhat    # ring Rp modulus p
dim   = (23, 47)  # dimensions of A in Rp^(m,n)

sigma_tdec = 1208

tail_bound_2 = 1.14

# Partition of s
wpart = [
    list(range(25)),
    list(range(25, 47))
]


# l2-norm bounds
wl2 = [0, 0 
       #    6*tail_bound_2 * sigma_tdec * sqrt(deg)
]


# binary coefficients (0 = no restriction)
wbin = [1, 1]

# rejection sampling flags (0 = not used)
wrej = [0 for _ in range(m+1)]

