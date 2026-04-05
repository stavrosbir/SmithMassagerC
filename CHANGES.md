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
