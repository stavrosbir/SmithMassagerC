# SmithMassagerC

A C implementation of the Smith massager algorithm for computing the
Smith normal form of nonsingular integer matrices.  The algorithm runs
in time softly equivalent to matrix multiplication:
O(n^ω · B(log n + log ‖A‖) · (log n)²) bit operations.

Based on the algorithms described in:

- S. Birmpilis, G. Labahn, A. Storjohann.  *A Las Vegas algorithm for
  computing the Smith form of a nonsingular integer matrix.*  ISSAC 2020.
- S. Birmpilis, G. Labahn, A. Storjohann.  *A fast algorithm for computing
  the Smith normal form with multipliers for a nonsingular integer matrix.*
  Journal of Symbolic Computation, 118:1–40, 2023.

## Origin

This is a fork of
[SmithMassager/SmithMassagerC](https://github.com/SmithMassager/SmithMassagerC),
originally developed at the University of Waterloo by Colton Pauderis,
Arne Storjohann, and collaborators.  See [CHANGES.md](CHANGES.md) for
the modifications made in this fork.

## Dependencies

- [GMP](https://gmplib.org/) — arbitrary-precision arithmetic
- [FLINT](https://flintlib.org/) — integer matrix operations (v3.x)
- [OpenBLAS](https://www.openblas.net/) — BLAS-accelerated linear algebra
- [IML](https://cs.uwaterloo.ca/~astorjoh/iml.html) — vendored (no install needed)

## Quick start

### macOS (Apple Silicon)

```bash
brew install gmp openblas flint
make lib
./bin/smith-massager 53
```

### macOS (Intel)

```bash
brew install gmp openblas flint
```

Edit `defs.mk` and change the Homebrew prefix from `/opt/homebrew/opt/`
to `/usr/local/opt/`, then:

```bash
make lib
./bin/smith-massager 53
```

### Linux (Ubuntu/Debian)

```bash
sudo apt install libgmp-dev libopenblas-dev libflint-dev
```

Edit `defs.mk` to point to your system paths (typically
`/usr/lib/` and `/usr/include/`), then:

```bash
make lib
./bin/smith-massager 53
```

## Build targets

| Target | Description |
|--------|-------------|
| `make lib` | Build static library (`lib/libhnfproj.a`) and `bin/smith-massager` |
| `make examples` | Build example binaries (`bin/hor_main`, `bin/simple_iherm`, etc.) |
| `make shared` | Build shared library (`lib/libhnfproj.so`) |
| `make clean` | Remove build artifacts |

## Usage

The `smith-massager` binary takes a prime number `n` as argument,
constructs an n×n test matrix A with entries A[i,j] = i^j mod n,
computes the Smith massager (U, M, T, S), and prints the result.

```
$ ./bin/smith-massager 5
```

Output includes the matrices U, M, T, S and a `success: 1` line
confirming the unimodularity certification passed.

## Configuration

All build configuration is in `defs.mk`:

- **Library paths**: Set `GMP_LIB_DIR`, `OPENBLAS_LIB_DIR`,
  `FLINT_LIB_DIR` and corresponding `_INCLUDE_DIR` variables.
- **Compiler**: `CC` defaults to `cc`.  Change to `gcc` if preferred.
- **Debug mode**: `make DEBUG=1` enables `-O0 -g` and debug print statements.
- **Maple**: Set `MAPLEDIR` to enable optional Maple integration (not required).

## Source layout

```
src/core/     Main algorithm: Smith massager, RNS arithmetic, projections
src/lift/     Double-plus-one lifting, RNS basis conversion
src/spinv/    Sparse inverse, high-order residues, unimodularity certification
src/iherm/    Packed Hermite form computation
src/iml-1.0.3/ Vendored IML library (Dixon lifting)
examples/     Example programs
tests/        Tuning benchmarks
```

## License

BSD 3-clause.  See [LICENSE](LICENSE).
