/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Rotation subdivision and accumulated frame deltas.

/**
 * @brief Verifies rotation_substeps returns a tight ceil with each sub-interval
 * within MAX.
 */
inline void test_rotation_substeps_shared_and_tight() {
  constexpr float MAX = 0.1f;
  // Always at least 1, even for a sub-threshold angle.
  HS_EXPECT_EQ(Animation::rotation_substeps(0.0f, MAX), 1);
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX * 0.5f, MAX), 1);
  // Tight ceil: N*MAX needs exactly N subdivisions.
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX, MAX), 1);
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX * 3.0f, MAX), 3);
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX * 3.2f, MAX), 4);
  // Every sub-interval stays within MAX.
  for (float a = 0.0f; a < 2.0f; a += 0.013f) {
    int n = Animation::rotation_substeps(a, MAX);
    HS_EXPECT_GE(n, 1);
    HS_EXPECT_LE(a / n, MAX + 1e-6f);
  }
}

/**
 * @brief Verifies Rotation::step does not discard sub-MIN_STEP_ANGLE increments.
 * @details A rotation slow enough that each frame's delta is below MIN_STEP_ANGLE
 * (1e-4 rad) must still accumulate those deltas so the orientation actually
 * turns over many frames. ease_linear is linear here, so the per-frame delta is
 * total_angle / duration.
 */
inline void test_rotation_accumulates_subthreshold_deltas() {
  using Ori = math::Orientation<16>;
  Ori o; // identity
  // 0.05 rad over 1000 frames => 5e-5 rad/frame, half of MIN_STEP_ANGLE, so every
  // frame's raw delta is below the early-out threshold.
  Animation::Rotation<288, 16> rot(o, math::Z_AXIS, 0.05f, 1000,
                                   math::ease_linear);
  for (int i = 0; i < 20; ++i)
    rot.step(fake_canvas());
  // A Z rotation sends +X toward +Y.
  math::Vector v = o.orient(math::X_AXIS, o.length() - 1);
  HS_EXPECT_GT(v.y, 5e-4f);
  HS_EXPECT_NEAR(v.x, 1.0f, 1e-3f);
}

/**
 * @brief Verifies a repeating Rotation lands its full sweep every cycle.
 * @details An easing with zero slope at t=1 leaves a sub-MIN_STEP_ANGLE residual on
 * the final frame. That frame has no successor to accumulate into, so dropping
 * it slips the residual (~1e-4 rad here) once per cycle.
 */
inline void test_rotation_applies_final_frame_residual() {
  using Ori = math::Orientation<16>;
  constexpr float ANGLE = 0.2f;
  constexpr int DURATION = 1000;
  constexpr int CYCLES = 10;
  Ori o; // identity
  Animation::Rotation<288, 16> rot(o, math::Z_AXIS, ANGLE, DURATION,
                                   math::ease_in_out_sin);
  for (int c = 0; c < CYCLES; ++c) {
    for (int i = 0; i < DURATION; ++i)
      rot.step(fake_canvas());
    rot.rewind();
  }
  // A Z rotation sends +X to (cos, sin) of the accumulated angle.
  math::Vector v = o.orient(math::X_AXIS, o.length() - 1);
  HS_EXPECT_NEAR(std::atan2(v.y, v.x), ANGLE * CYCLES, 2e-4f);
}
