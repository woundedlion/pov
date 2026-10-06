/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Driver
// ============================================================================

/**
 * @brief Verifies a wrapping Driver adds its per-frame increment each step and
 * wraps back into [0,1) when it reaches 1.
 */
inline void test_driver_increments_and_wraps() {
  float v = 0.0f;
  Animation::Driver d(v, 0.25f, /*wrap=*/true);
  d.step(fake_canvas());
  HS_EXPECT_NEAR(v, 0.25f, 1e-5f);
  d.step(fake_canvas());
  d.step(fake_canvas());
  HS_EXPECT_NEAR(v, 0.75f, 1e-5f);
  d.step(fake_canvas()); // 1.0 -> wraps to 0.0
  HS_EXPECT_NEAR(v, 0.0f, 1e-5f);
  HS_EXPECT_NEAR(d.get_mutant(), 0.0f, 1e-5f);
}

/**
 * @brief Verifies a non-wrapping Driver accumulates its increment past 1
 * without bound.
 */
inline void test_driver_no_wrap_accumulates() {
  float v = 0.0f;
  Animation::Driver d(v, 0.5f, /*wrap=*/false);
  for (int i = 0; i < 5; ++i)
    d.step(fake_canvas());
  HS_EXPECT_NEAR(v, 2.5f, 1e-5f);
}

/**
 * @brief Verifies a live-bound Driver ignores a non-finite slider frame instead
 * of permanently poisoning the wrapped mutant via wrap_t(NaN).
 */
inline void test_driver_nan_source_does_not_poison() {
  float v = 0.0f;
  float slider = 0.25f;
  Animation::Driver d(v, &slider, /*scale=*/1.0f, /*wrap=*/true);
  d.step(fake_canvas());
  HS_EXPECT_NEAR(v, 0.25f, 1e-5f);

  slider = std::numeric_limits<float>::quiet_NaN();
  d.step(fake_canvas());
  HS_EXPECT_TRUE(std::isfinite(v)); // NaN frame kept the last good speed
  HS_EXPECT_NEAR(v, 0.5f, 1e-5f);

  slider = 0.25f;
  d.step(fake_canvas());
  HS_EXPECT_NEAR(v, 0.75f, 1e-5f);
}
