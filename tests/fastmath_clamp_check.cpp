/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// Runs clamp-before-cast, NaN-saturation and affine-form cases under WASM
// release math flags (-O3 -ffast-math -fno-finite-math-only).
#include "core/engine/engine.h"
#include "effects/DisplacementField.h"
#include "tests/test_color.h"
#include "tests/test_harness.h"
#include "tests/test_interpolate.h"
#include "tests/test_pullback.h"

#include <limits>

namespace hs_test::effects_tests {
struct DisplacementFieldWhiteBox {
  static bool sweep_uses_cache(bool present) {
    using Field = DisplacementField<96, 20>;
    static Field::SweepCache cache;
    Field::Sweep sweep{};
    sweep.cache = present ? &cache : nullptr;
    return sweep.uses_cache();
  }
};
} // namespace hs_test::effects_tests

int main() {
  static_assert(HS_ENABLE_TEST_ORACLES && !HS_ENABLE_TEST_HOOKS);
  using namespace hs_test::color_tests;

  const hs_test::ModuleScope scope = hs_test::begin_module(
      "fastmath_clamp (-ffast-math -fno-finite-math-only)");
  HS_EXPECT_FALSE(
      hs_test::effects_tests::DisplacementFieldWhiteBox::sweep_uses_cache(
          false));
  HS_EXPECT_TRUE(
      hs_test::effects_tests::DisplacementFieldWhiteBox::sweep_uses_cache(
          true));
#define HS_RUN_CLAMP_TEST(fn) fn();
  HS_FASTMATH_CLAMP_TESTS(HS_RUN_CLAMP_TEST)
#undef HS_RUN_CLAMP_TEST
  hs_test::interpolate_tests::test_clamp_progress();
  const volatile float from = 0.9817478656768799f;
  const volatile float to = -2.1598448753356934f;
  const volatile float fraction = 0.0625f;
  HS_EXPECT_EQ(interp::shortest_periodic(from, to, fraction, math::TWO_PI_F,
                                         -math::PI_F),
               0.7853983044624329f);
  HS_EXPECT_EQ(math::fast_acos(std::numeric_limits<float>::quiet_NaN()), 0.0f);
  hs_test::pullback_tests::test_pullback_affine_forms_match_general();
  return hs_test::end_module(scope) ? 1 : 0;
}
