/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
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
  if (total == 0) {
    std::printf("=== fastmath_clamp: NO ASSERTIONS RAN ===\n");
    return 1;
  }
  return failed ? 1 : 0;
}
