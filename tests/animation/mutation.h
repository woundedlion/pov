/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

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
  HS_EXPECT_FALSE(m.done());
  m.step(fake_canvas());
  HS_EXPECT_TRUE(m.done());
  HS_EXPECT_TRUE(std::isfinite(v));
  HS_EXPECT_NEAR(v, 1.0f, 1e-3f);
}

/** @brief Non-finite mutation results retain the preceding finite value. */
inline void test_mutation_nonfinite_result_holds_prior_value() {
  float value = 0.0f;
  Animation::Mutation mutation(
      value,
      [](float t) {
        return t == 0.5f ? std::numeric_limits<float>::quiet_NaN() : t;
      },
      4, math::ease_linear);
  for (float expected : {0.25f, 0.25f, 0.75f, 1.0f}) {
    mutation.step(fake_canvas());
    HS_EXPECT_EQ(value, expected);
  }
}
