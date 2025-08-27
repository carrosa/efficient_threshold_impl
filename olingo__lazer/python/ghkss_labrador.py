from lazer import *
from labrador import *
import time
import hashlib  # for SHAKE128


m = 22
lhat =22
sigma_tdec = 1208

big_mod = 633544653872304603170707504554316801
small_mod = 2101003516513793
deg = 256 # Degree of the polynomial ring
witness_num_polys = 1 # Number of polynomials per witness
statement_num_polys = m*lhat # Number of polynomials in the statement


shake128 = hashlib.shake_128(bytes.fromhex("00"))
TARGPP = shake128.digest(32)

qhat_ring = polyring_t(deg, big_mod)
q_ring = polyring_t(deg, small_mod)

qhat_primesize = str(math.ceil(math.log2(qhat_ring.mod)))
q_primesize = str(math.ceil(math.log2(q_ring.mod)))

SAME_KEY = 1 # Use the same secret/public key for all signatures to save time on key generation
ID = int_to_poly(1, qhat_ring)

deg_list = [deg] * statement_num_polys

# Number of polynomials per witness
num_pols_list = [witness_num_polys] * statement_num_polys

# Norm bounds for each witness
norms = [q * math.sqrt(deg)] + [t2 * sigma_tdec * math.sqrt(deg)]
norm_list = norms * statement_num_polys

# Number of constraints
num_constraints = 50

# Initialize the proof statement
PS = proof_statement(deg_list, num_pols_list, norm_list, num_constraints, qhat_primesize)
