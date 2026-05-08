#include <assert.h>

#include "largestInvariantFactor.h"
#include "timer.h"
#include "basic.h"
#include "fmpz_rand.h"


/* Compute s, a multiple of the largest invariant factor s_n of A,
 * along with the projection X = s * A^{-1} * B mod s for a random B.
 *
 * One IML solve gives X_0 = d * A^{-1} * B, where d divides s_n.
 * fmpz_mat_extractLCM gives sn_est = d / gcd(content(X_0), d), the
 * LCM of the column denominators of A^{-1} * B; for uniform random
 * B this equals s_n with overwhelming probability.
 *
 * We set s = lcm(sn_est, d) so that d | s, then compute the
 * projection X = (s / d) * X_0 mod s exactly.
 *
 * IMPORTANT: this routine uses fmpz_mat_randbits (uniform random),
 * not fmpz_mat_randtest.  fmpz_mat_randtest is FLINT's "edge-case"
 * random generator and intentionally produces sparse / structured
 * matrices that break the randomness assumption of the projection.
 */
int largestInvariantFactor(fmpz_t s, fmpz_mat_t X, fmpz_mat_t A, slong startDim) {
  assert(A->c == A->r);
  assert(startDim > 0);

  fmpz_t mx, d, q;
  fmpz_mat_t B;
  int n = A->c;

  fmpz_init(mx);
  fmpz_init(d);
  fmpz_init(q);

  /* Bit-bound for the random entries of B.  This formula bounds
   * the size needed so that A^{-1} * B captures s_n with high
   * probability for the test matrices in use. */
  fmpz_mat_max(mx, A);
  fmpz_mul_ui(mx, mx, fmpz_bits(mx));
  fmpz_mul_ui(mx, mx, n);
  fmpz_mul_ui(mx, mx, bits(n));
  fmpz_mul_ui(mx, mx, 20);
  fmpz_add_ui(mx, mx, 60);

  fmpz_mat_init(B, n, startDim);
  fmpz_mat_randbits(B, rand, fmpz_bits(mx));

  REAL_TIMER("ImlSolve in largestInvariantFactor",
             fmpz_mat_imlSolve(X, d, A, B));

  /* sn_est = d / gcd(content(X), d). */
  fmpz_mat_extractLCM(s, X, d);

  /* Ensure d | s so the projection X = (s/d) * X_0 mod s is exact. */
  fmpz_lcm(s, s, d);
  fmpz_divexact(q, s, d);
  fmpz_mat_scalar_mul_fmpz(X, X, q);
  fmpz_mat_mod(X, X, s);

  fmpz_mat_clear(B);
  fmpz_clear(mx);
  fmpz_clear(d);
  fmpz_clear(q);

  return 1;
}
