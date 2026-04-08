# Architecture

This document describes the structure of the SmithMassagerC implementation.
The algorithm is from Birmpilis, Labahn, and Storjohann (ISSAC 2020, JSC 2023);
the C implementation is by Ziwen Wang.

## High-level flow

```
main(p)
  └─ interestHerm(A)          // generate test matrix A[i,j] = i^j mod p
  └─ smithMassager(U, M, T, S, A)
       └─ smithMassagerHelper(...)     // retry wrapper (Las Vegas)
```

`smithMassagerHelper` does three things:

**Phase 1.** `largestInvariantFactor` — find s_n (or a multiple) via IML's
Dixon lifting solver.  Solve AX = B for random B, extract LCM of denominators.

**Phase 2.** Main loop (`while sumR < n`) — peel off invariant factors in
batches of size r:

```
indexMassager(Sp, Up, Mp, Tp, B, n, sumR, r, s, k, Q)
  ├─ specialIntCert(...)     // compute P = sB^{-1}J mod s
  │    └─ fastIntCert(...)
  │         └─ horBase(...)  // high-order residue via p-adic lifting
  │              └─ initLift / lift  (double-plus-one lifting)
  │                   └─ rnsMatrix_gemm  (BLAS-accelerated)
  │
  └─ computeProjBasis(...)   // extract SNF from the projection
       └─ SNF(...)           // recursive divide-and-conquer
            └─ stab, vectorGCD, extractMatrixGCDHelper
```

Then update B and accumulate into (U, M, T, S).

**Phase 3.** `fmpz_uniCert` — certify |det B| = 1 via p-adic lifting until
the residue vanishes.

## Source tree

### src/core/ — main algorithm and linear algebra infrastructure

| File | Role |
|------|------|
| `main.c` | Entry point: parse CLI arg, generate test matrix, call `smithMassager` |
| `smithMassager.c` | Top-level algorithm, iteration loop, B matrix updates |
| `indexMassager.c` | Compute index-(m,r) Smith massager for one batch |
| `projection.c` | Recursive divide-and-conquer SNF of the projection mod s |
| `basic.c` | `stab` primitive, `vectorGCD`, `extractMatrixGCDHelper` |
| `largestInvariantFactor.c` | Phase 1: compute s_n via IML |
| `specialIntCert.c` | Integrality certification exploiting block structure of B |
| `fastIntCert.c` | Fast integrality certificate via high-order residue |
| `rns_matrix.h` | RNS matrix type (`rnsMatrix_t`) |
| `basis.c` | RNS basis: prime selection, Lagrange interpolants, CRT weights |
| `openblas_residue.c` | Modular matrix stored as doubles, multiplied via `cblas_dgemm` |
| `reconstruct.c` | CRT reconstruction (RNS → integer) |
| `pl_matrix.c` | Partial linearization: bit-slice large-entry matrices |
| `mpz_matrix.c` | GMP-based integer matrix and operations |
| `fmpz_conv.c` | Conversion between FLINT `fmpz_mat_t` and GMP `mpzMatrix_t` |
| `fmpz_mat_helper.c` | Utilities on `fmpz_mat_t`: cmod, imlSolve, uniCert wrapper |
| `arith_utils.c` | GMP arithmetic: `gcdex`, `stab`, `modsL`, `prevprime` |
| `mods.c` | Inline modular reduction for doubles (included into `openblas_residue.c`) |
| `genHerm.c` | Test matrix generator: A[i,j] = i^j mod n |
| `thresholds.h` | Tuning constants for algorithm crossovers |
| `timer.h` | `REAL_TIMER` macro for timing instrumentation |

### src/lift/ — p-adic lifting engine

| File | Role |
|------|------|
| `lift.c` | Double-plus-one lifting: `initLift`, `lift` (each call doubles precision) |
| `rns_conversion.c` | Convert between RNS bases (Lagrange interpolation or reconstruct-reduce) |
| `imlsolve.c` | Wrapper around IML's `nonsingSolvLlhsMM` |

### src/spinv/ — sparse inverse and certification

| File | Role |
|------|------|
| `highorder.c` | High-order residue: iterate lifting until R = 0 or stabilize |
| `unicert.c` | Unimodularity certification (Phase 3) |
| `spinv.c` | Sparse p-adic expansion of A^{-1} (used by HNF, not by Smith massager) |

### src/iherm/ — Hermite Normal Form (independent module)

| File | Role |
|------|------|
| `iherm.c` | Iterative HNF algorithm |
| `pk_matrix.c` | Packed triangular matrices (only non-identity columns stored) |
| `hnfproj.c` | HNF of a single projection |

### src/iml-1.0.3/ — vendored IML library

Dixon lifting for exact linear system solving over Z.

## Key data structures

| Type | Defined in | Purpose |
|------|-----------|---------|
| `fmpz_mat_t` | FLINT (external) | Primary matrix type in the Smith massager |
| `mpzMatrix_t` | `mpz_matrix.h` | GMP-based matrix, used in lifting/RNS subsystems |
| `rnsMatrix_t` | `rns_matrix.h` | Matrix as ℓ residue images mod word-sized primes |
| `basis_t` | `basis.h` | RNS basis: primes, Lagrange interpolants, CRT weights |
| `openblasResidue_t` | `openblas_residue.h` | Single residue matrix (doubles), multiplied via `cblas_dgemm` |
| `lift_info_t` | `lift.h` | State for double-plus-one lifting (two RNS bases, R, M) |
| `pk_t` | `pk_matrix.h` | Packed upper triangular (only non-identity columns) |
| `plMatrix_t` | `pl_matrix.h` | Partial linearization: large entries split into bit slices |

## The fundamental primitive: stab

The projection/SNF computation rests on `stab(c, a, b, N)`: find the
minimal c ≥ 0 such that gcd(a + cb, N) = gcd(a, b, N).  This drives
`vectorGCD` (find u with ⟨u, p⟩ = gcd(p_1, …, p_n) mod s) and
`extractMatrixGCDHelper` (find a column combination achieving the content).
