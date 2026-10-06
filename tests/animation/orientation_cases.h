/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_animation.h.

// Orientation upsampling and collapse.

/**
 * @brief Verifies Orientation::upsample SLERP-interpolates the recorded
 * sub-frames up to a target count (preserving endpoints) and collapse()
 * discards all but the newest.
 * @details These are the two primitives multi-animation motion blur is built
 * on.
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
