from math import sqrt, exp
from sage.all import var, solve

vname = "param"

ringdim = 256

q_small = 2101003516513793

moddimEnc = 22
secdimEnc = 22
mhat = 22
q = 633544653872304603170707504554316801

moddim = 6
secdim = 6
# sigma_w = 2**38 for signing keygen proof, sigma_w = 2**20 for keygen proof
# sigma_w = 2**38 
sigma_w = 2**20 # Actually sigma_y

rej_rate = 0.01
tail_bound_2 = 1.01
while ((tail_bound_2**ringdim) * exp((ringdim / 2) * (1 - tail_bound_2**2))) > rej_rate:
    tail_bound_2 += 0.01


deg   = 256  # ring Rp degree d
mod   = q    # ring Rp modulus p
dim   = (secdimEnc + mhat, moddimEnc + secdimEnc + mhat + 2*secdim)  # dimensions of A in Rp^(m,n)

# Partition of s
wpart = (
    [[i] for i in range(moddimEnc + secdimEnc + mhat)] +
    [[i] for i in range(moddimEnc + secdimEnc + mhat, moddimEnc + secdimEnc + mhat + secdim)] + # first 6 elements of m0
    [[i] for i in range(moddimEnc + secdimEnc + mhat + secdim, moddimEnc + secdimEnc + mhat + 2*secdim)] # first 6 elements of m1
)

B = tail_bound_2 * sigma_w * sqrt(ringdim)
x = var('x')
eq = x**2 *sqrt(ringdim)+sqrt(ringdim)*x - B == 0
solutions = solve(eq, x)
beta = int(solutions[0].rhs().n() if solutions[0].rhs().n() > 0 else solutions[1].rhs().n())

# l2-norm bounds
wl2 = (
    [tail_bound_2 * sqrt(2/3) * sqrt(ringdim)*4 for _ in range(moddimEnc + secdimEnc + mhat)] +
    [beta * sqrt(ringdim)*4 for _ in range(moddimEnc + secdimEnc + mhat, moddimEnc + secdimEnc + mhat + secdim)] +
    [beta * sqrt(ringdim)*4 for _ in range(moddimEnc + secdimEnc + mhat + secdim, moddimEnc + secdimEnc + mhat + 2*secdim)]
)


# binary coefficients (0 = no restriction)
wbin = [0 for _ in range(moddimEnc + secdimEnc + mhat + 2*secdim)]

# rejection sampling flags (0 = not used)
wrej = [0 for _ in range(moddimEnc + secdimEnc + mhat + 2*secdim)]

# Optional: linf-norm bound on s (commented out)
# wlinf = 1

