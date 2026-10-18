/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Orientation upsampling and collapse.

/**
 * @brief Verifies Orientation::upsample SLERP-interpolates the recorded
 * sub-frames up to a target count (preserving endpoints) and collapse()
 * discards all but the newest.
 */
inline void test_orientation_upsample_then_collapse() {
  math::Orientation<8> o; // identity, 1 frame
  o.push(math::make_rotation(
      math::Z_AXIS, math::PI_F / 2)); // 2 frames: identity, +90 about Z
  HS_EXPECT_EQ(o.length(), 2);

  o.upsample(5);
  HS_EXPECT_EQ(o.length(), 5);

  // Endpoints preserved: frame 0 ~ identity (+X stays +X); frame 4 ~ the +90
  // rotation about Z (+X -> +Y).
  math::Vector f0 = o.orient(math::X_AXIS, 0);
  HS_EXPECT_NEAR(f0.x, 1.0f, 1e-3f);
  math::Vector f4 = o.orient(math::X_AXIS, 4);
  HS_EXPECT_NEAR(f4.y, 1.0f, 1e-3f);

  // SLERP is monotone: +X decreases across the interpolated frames.
  float prevx = 2.0f;
  for (int i = 0; i < 5; ++i) {
    float x = o.orient(math::X_AXIS, i).x;
    HS_EXPECT_LE(x, prevx + 1e-4f);
    prevx = x;
  }

  o.collapse();
  HS_EXPECT_EQ(o.length(), 1);
  math::Vector c = o.orient(math::X_AXIS, 0);
  HS_EXPECT_NEAR(c.y, 1.0f, 1e-3f);
}

/**
 * @brief `Orientation::get(float)` hits every history frame at its fraction,
 *        slerps between them, and returns the lone frame of a still history.
 */
inline void test_orientation_get_fraction_follows_history() {
  math::Orientation<> o(math::make_rotation(math::Y_AXIS, 0.25f));
  o.push(math::make_rotation(math::Y_AXIS, 0.35f));
  o.push(math::make_rotation(math::Y_AXIS, 0.75f));
  for (int i = 0; i < o.length(); ++i)
    HS_EXPECT_NEAR(std::fabs(math::dot(o.get(0.5f * i), o.get(i))), 1.0f,
                   1e-5f);
  const math::Quaternion mid = o.get(0.75f);
  const math::Quaternion want = math::make_rotation(math::Y_AXIS, 0.55f);
  HS_EXPECT_NEAR(std::fabs(math::dot(mid, want)), 1.0f, 1e-5f);
  const math::Quaternion early = o.get(0.25f);
  const math::Quaternion early_want = math::make_rotation(math::Y_AXIS, 0.30f);
  HS_EXPECT_NEAR(std::fabs(math::dot(early, early_want)), 1.0f, 1e-5f);

  const math::Orientation<> still(math::make_rotation(math::X_AXIS, 0.3f));
  HS_EXPECT_NEAR(std::fabs(math::dot(still.get(0.4f), still.get())), 1.0f,
                 1e-6f);
}
