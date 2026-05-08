# Changes from upstream

This fork is based on
[SmithMassager/SmithMassagerC](https://github.com/SmithMassager/SmithMassagerC)
(commit `f7f92ca`).  The following changes were made to fix bugs, enable
building on macOS without Maple, and clean up the code.

## Bug fixes

- **Fixed buffer overflow in threaded RNS operations** (`src/core/rns_matrix.c`,
  `src/lift/rns_conversion.c`):  The thread and argument arrays were declared
  as fixed-size `[32]` on the stack.  If the RNS basis has more than 32
  moduli (which occurs for large matrices with big entries), this overflows
  and causes memory corruption.  Changed to dynamic allocation via `malloc`.
- **Non-deterministic PRNG seeding** (`src/core/fmpz_rand.c`):  The FLINT
  random state was initialized with a fixed default seed, making Las Vegas
  restarts deterministic.  Now seeded from `time()` and `getpid()`.  Also
  fixed the `initialized` flag which was never set to 1.
- **Use `fmpz_mat_randbits` instead of `fmpz_mat_randtest` for random
  projections** (`src/core/largestInvariantFactor.c`,
  `src/core/indexMassager.c`):  `fmpz_mat_randtest` is FLINT's
  edge-case generator and intentionally produces sparse / structured
  matrices.  Using it for the random projection violated the algorithm's
  randomness assumption and caused the index massager to fail repeatedly
  (15-20% per-attempt success rate, dozens of restarts).  With
  `fmpz_mat_randbits` (uniform random) the algorithm essentially always
  succeeds in one attempt, as predicted by the theory.

## Algorithm correctness improvements

- **Multi-solve LIF** (`src/core/largestInvariantFactor.c`):  The
  `largestInvariantFactor` routine now performs `nsolves=5` independent
  IML solves with fresh random `B`, taking the LCM of the resulting
  denominator estimates.  A single solve only captures a divisor of
  `s_n` (with high probability equal to `s_n` for generic matrices,
  but sometimes a proper divisor); taking the LCM across multiple
  solves makes the estimate of `s_n` essentially exact.
- **Dynamic projection width** (`src/core/indexMassager.c`,
  `src/core/smithMassager.c`):  The number of "extra" projection
  columns `k` was hardcoded to 5.  Theory requires
  `k = ceil(log2 r) + c` for per-call success probability
  >= 1 - 1/2^c.  Now `k = 2*ceil(log2 r) + 20`, which gives a
  per-call success probability essentially 1.

## Build system

- **Removed unconditional Maple dependency** (`defs.mk`):  The base `LDLIBS`
  included `-lmaplec`, causing link failures on systems without Maple.
  Removed it; the Maple flag is already added conditionally in `src/Makefile`
  when `MAPLEDIR` is set.
- **macOS ARM64 defaults** (`defs.mk`):  Set Homebrew paths for GMP, OpenBLAS,
  and FLINT under `/opt/homebrew/opt/`.  Set `CC` to `cc` (Apple Clang).
  Set `FLINT_INCLUDE_DIR` to `.../include/flint/` to match Homebrew's layout.
- **Conditional Maple file copy** (`Makefile`):  Wrapped the `cp
  src/maple/extern.mpl lib/` step in `ifdef MAPLEDIR` so it does not fail
  when Maple is not installed.

## Code quality

- **Guarded debug prints** (`src/core/smithMassager.c`,
  `src/core/projection.c`):  Wrapped diagnostic `printf` calls in
  `#ifdef DEBUG` so they only appear when compiling with `make DEBUG=1`.
- **Removed dead code** (`src/core/projection.c`):  Removed three blocks of
  manually duplicated `gettimeofday` timing code marked "TODO: DELETE THIS
  AFTER".  Replaced with the existing `REAL_TIMER` macro.
- **Fixed typo** (`src/core/projection.c`):  "divisble" → "divisible".
