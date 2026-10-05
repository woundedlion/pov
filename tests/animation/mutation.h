/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_animation.h.

// ============================================================================
// Mutation
// ============================================================================

/**
 * @brief Verifies each step Mutation writes f(easing(t/duration)) into its
 * bound float.
 */
inline void test_mutation_applies_function_of_eased_time() {
  float v = 0.0f;
  const int duration = 5;
  Animation::Mutation m(
      v, [](float e) { return e * 10.0f; }, duration,
      [](float t) { return t * t; });
  constexpr float EXPECTED[] = {0.4f, 1.6f, 3.6f, 6.4f, 10.0f};
  for (int i = 0; i < duration; ++i) {
    m.step(fake_canvas());
    HS_EXPECT_NEAR(v, EXPECTED[i], 1e-3f);
  }
  HS_EXPECT_NEAR(v, 10.0f, 1e-3f);
  HS_EXPECT_TRUE(m.done());
}

/**
 * @brief Verifies a zero-duration Mutation coerces to 1 frame and yields a
 * finite result.
 */
inline void test_mutation_duration_zero_finite() {
  float v = 0.0f;
  Animation::Mutation m(v, [](float e) { return e; }, 0, math::ease_linear);
  m.step(fake_canvas());
  HS_EXPECT_TRUE(std::isfinite(v));
  HS_EXPECT_NEAR(v, 1.0f, 1e-3f);
}
