#include "fmpz_rand.h"
#include <time.h>
#include <unistd.h>

int initialized = 0;
flint_rand_t rand;

void fmpz_rand_init() {
  if (!initialized) {
    flint_rand_init(rand);
    ulong seed = (ulong)time(NULL) ^ ((ulong)getpid() << 16);
    flint_rand_set_seed(rand, seed, seed * 1103515245UL + 12345UL);
    initialized = 1;
  }
}

void fmpz_rand_clear() {
  if (initialized) {
    flint_rand_clear(rand);
  }
}
