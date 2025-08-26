# Benchmarking Olingo

## Passive security scheme:

We build upon dilithiums codebase using GMP for larger integers.
Passive signing benchmarks can be generated with `make ghkss`, from `test_GHKSS256.c`.
Passive DKG and keygen benchmarks can be generated with `make bench`, from `benchmark_ghkss.c`.

Information about requirements to run is found inside the folder `olingo_passive__using_dilithium_polynomial_funcitons`.

To benchmark with different thresholds and total parties, one needs to change the THRESHOLD and USERS variables in params.h. Keep in mind benchmarking DKG with large THRESHOLD and USERS is slow and requires quite a bit of memory.

We simulate 1 party only, meaning we use dummy data for the benchmarks where applicable to simulate having multiple parties.

## LIN proofs:

We adapt https://github.com/dfaranha/lattice-verifiable-mixnet for BDLOP commitments and linear proofs. See its readme for more info on how to build.

Information about requirements to run is found inside the folder `lin_proofs_for_olingo__based_on_lattice-verifiable-mixnet_code`

## LNP proofs with lazer:

We use the python interface for Lazer for LNP style proofs.
The folders for our benchmarks:

* python/bnd_E_prf
* python/dkg_prf1
* python/sign_prf

To build each of the benchmarks, run `make` in their respective folders.

More information about requirements to run is found inside the folder `lazer_with_python_proofs_for_olingo`, in the file `lazer_getting_started.html` and in the github repository for lazer: https://github.com/lazer-crypto/lazer


__WARNING__: This is an academic proof of concept/benchmarking, and in particular has not received code review. This implementation is NOT ready for any type of production use.
