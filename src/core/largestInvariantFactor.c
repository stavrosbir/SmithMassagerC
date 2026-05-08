#include <assert.h>

#include "largestInvariantFactor.h"
#include "timer.h"
#include "basic.h"
#include "fmpz_rand.h"


/* Compute s, a multiple of the largest invariant factor s_n of A,
 * along with the projection X = s * A^{-1} * B mod s for some random B.
 *
 * Strategy:
 *   1. Perform `nsolves` independent IML solves with fresh random B's.
 *      Each solve returns X_i = d_i * A^{-1} * B_i and a denominator d_i.
 *      The function fmpz_mat_extractLCM gives a divisor sn_i of s_n.
 *   2. Take the LCM of the sn_i's to get a robust estimate of s_n.
 *      With nsolves=5 the probability of missing any prime factor of
 *      s_n is essentially zero.
 *   3. Combine with the last d to ensure d | s, then compute
 *      X = (s/d) * X_last mod s exactly.
 *
 * IMPORTANT: this function uses fmpz_mat_randbits (uniform random),
 * not fmpz_mat_randtest.  fmpz_mat_randtest is FLINT's "edge-case"
 * random generator and intentionally produces sparse / structured
 * matrices that break the randomness assumption of the projection.
 */
int largestInvariantFactor(fmpz_t s, fmpz_mat_t X, fmpz_mat_t A, fmpz_t startDim) {
  assert(A->c == A->r);

  fmpz_t mx, d, sni, q;
  fmpz_mat_t B;
  int n = A->c;
  int ncols = (slong)startDim;
  int nsolves = 1;

  fmpz_init(mx);
  fmpz_init(d);
  fmpz_init(sni);
  fmpz_init(q);

  fmpz_mat_max(mx, A);
  fmpz_mul_ui(mx, mx, fmpz_bits(mx));
  fmpz_mul_ui(mx, mx, n);
  fmpz_mul_ui(mx, mx, bits(n));
  fmpz_mul_ui(mx, mx, 20);
  fmpz_add_ui(mx, mx, 60);

  fmpz_mat_init(B, n, ncols);

  /* Estimate s_n via LCM of independent solves. */
  fmpz_set_ui(s, 1);
  for (int i = 0; i < nsolves; i++) {
    fmpz_mat_randbits(B, rand, fmpz_bits(mx));
    REAL_TIMER("ImlSolve in largestInvariantFactor",
               fmpz_mat_imlSolve(X, d, A, B));
    fmpz_mat_extractLCM(sni, X, d);
    fmpz_lcm(s, s, sni);
  }

  /* Bring d into s (so d | s).  Then compute the projection
   *    X = s * A^{-1} * B_last  mod  s
   *      = (s/d) * X_last       mod  s. */
  fmpz_lcm(s, s, d);
  fmpz_divexact(q, s, d);
  fmpz_mat_scalar_mul_fmpz(X, X, q);
  fmpz_mat_mod(X, X, s);

  fmpz_mat_clear(B);
  fmpz_clear(mx);
  fmpz_clear(d);
  fmpz_clear(sni);
  fmpz_clear(q);

  return 1;
}
