// Compiles the same NaN-clamp assertions as run_tests but under the WASM release
// math flags (-O3 -ffast-math -fno-finite-math-only, set in tests/CMakeLists.txt)
// so the hs::clamp NaN->hi contract is verified against the codegen that ships,
// not only the default-IEEE-semantics build.
#include <cstdio>

#include "core/engine/engine.h"
#include "tests/test_color.h"

int main() {
  using namespace hs_test::color_tests;

#define HS_RUN_CLAMP_TEST(fn) fn();
  HS_FASTMATH_CLAMP_TESTS(HS_RUN_CLAMP_TEST)
#undef HS_RUN_CLAMP_TEST

  const int failed = hs_test::stats().failed;
  const int total = hs_test::stats().passed + failed;
  std::printf("=== fastmath_clamp: %d passed, %d failed (-ffast-math "
              "-fno-finite-math-only) ===\n",
              hs_test::stats().passed, failed);
  // Floor against silent drift: a gutted test body would otherwise stay green.
  // The case list itself is shared with test_color.h via HS_FASTMATH_CLAMP_TESTS,
  // so a new clamp test extends this pass automatically. Bump when adding cases.
  constexpr int MIN_ASSERTIONS = 26;
  if (total < MIN_ASSERTIONS) {
    std::printf("=== fastmath_clamp: only %d assertions ran, expected >= %d "
                "(a check was dropped) ===\n",
                total, MIN_ASSERTIONS);
    return 1;
  }
  return failed ? 1 : 0;
}
