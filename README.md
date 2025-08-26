# Benchmarking Olingo

## Passive security scheme:

We build upon dilithiums codebase using GMP for larger integers.
Passive signing benchmarks can be generated with `make ghkss`, from `test_GHKSS256.c`.
Passive DKG and keygen benchmarks can be generated with `make bench`, from `benchmark_ghkss.c`.

We simulate 1 party only, meaning for the benchmarks, we use dummy data where applicable.

## LIN proofs:

We adapt https://github.com/dfaranha/lattice-verifiable-mixnet for BDLOP commitments and linear proofs. See its readme for more info on how to build.

## LNP proofs with lazer:

We use the python interface for Lazer for LNP style proofs.
The folders for our benchmarks:

* python/bnd_E_prf
* python/dkg_prf1
* python/sign_prf



__WARNING__: This is an academic proof of concept/benchmarking, and in particular has not received code review. This implementation is NOT ready for any type of production use.
