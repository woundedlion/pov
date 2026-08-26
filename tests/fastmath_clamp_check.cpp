/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// Compiles the same NaN-clamp assertions as run_tests but under the WASM release
// math flags (-O3 -ffast-math -fno-finite-math-only, set in tests/CMakeLists.txt)
// so the hs::clamp NaN->hi contract is verified against the codegen that ships,
// not only the default-IEEE-semantics build.
#include "core/engine/engine.h"
#include "tests/test_color.h"
#include "tests/test_harness.h"

int main() {
  using namespace hs_test::color_tests;

  const hs_test::ModuleScope scope = hs_test::begin_module(
      "fastmath_clamp (-ffast-math -fno-finite-math-only)");
#define HS_RUN_CLAMP_TEST(fn) fn();
  HS_FASTMATH_CLAMP_TESTS(HS_RUN_CLAMP_TEST)
#undef HS_RUN_CLAMP_TEST
  return hs_test::end_module(scope) ? 1 : 0;
}
