from math import sqrt, exp
# Unused for the sake of the examples.
vname = "param"

ringdim = 256

q_small = 2101003516513793

moddimEnc = 22
secdimEnc = 22
mhat = 22
q = 633544653872304603170707504554316801

rej_rate = 0.01

tail_bound_2 = 1.01
while ((tail_bound_2**ringdim) * exp((ringdim / 2) * (1 - tail_bound_2**2))) > rej_rate:
    tail_bound_2 += 0.01

deg   = 256  # ring Rp degree d
mod   = q    # ring Rp modulus p
dim   = (secdimEnc + mhat, moddimEnc + secdimEnc + 2*mhat)  # dimensions of A in Rp^(m,n)

# Partition of s
wpart = (
    [[i] for i in range(moddimEnc + secdimEnc + mhat)] +
    [[i] for i in range(moddimEnc + secdimEnc + mhat, moddimEnc + secdimEnc + 2*mhat)]
)

# l2-norm bounds
wl2 = (
    [tail_bound_2 * sqrt(2/3) * sqrt(ringdim) for _ in range(moddimEnc + secdimEnc + mhat)] +
    [q_small * sqrt(ringdim) for _ in range(moddimEnc + secdimEnc + mhat, moddimEnc + secdimEnc + 2*mhat)]
)

# binary coefficients (0 = no restriction)
wbin = [0 for _ in range(moddimEnc + secdimEnc + 2*mhat)]

# rejection sampling flags (0 = not used)
wrej = [0 for _ in range(moddimEnc + secdimEnc + 2*mhat)]

# Optional: linf-norm bound on s (commented out)
# wlinf = 1
