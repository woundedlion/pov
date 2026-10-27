/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief RingSpin white-box cases. */

/** @brief White-box accessor for RingSpin's ring pool. */
struct RingSpinWhiteBox {
  using RS = RingSpin<DEFAULT_W, DEFAULT_H>;

  /** @brief Sub-frame count of ring @p i's current orientation step. */
  static int substeps(const RS &fx, int i) {
    return fx.rings[i].orientation.length();
  }
  /** @brief World-space axis of ring @p i's great circle at sub-frame @p s. */
  static math::Vector axis(const RS &fx, int i, int s) {
    return fx.rings[i].orientation.orient(math::Y_AXIS, s).normalized();
  }
  /** @brief Ring-pool size. */
  static int num_rings() { return RS::NUM_RINGS; }
  /** @brief Half-width in radians of the trail-head stroke. */
  static float head_half_width(const RS &fx) {
    const float pixel_w = math::coarse_pixel_pitch<DEFAULT_W, DEFAULT_H>();
    return 2.0f * pixel_w * fx.params.thickness;
  }
};

/** @brief Strobe configuration preserves a nonblack RingSpin render. */
inline void test_ringspin_strobe_configuration_preserves_rendering() {
  std::vector<Pixel> expected;
  bool lit = false;
  render_capture<RingSpin, SMALL_W, SMALL_H>(expected, 2, nullptr, &lit);
  HS_EXPECT_TRUE(lit);
  reset_effect_globals();
  pin_frame_clock(0);
  {
    RingSpin<SMALL_W, SMALL_H> roster_effect;
    HS_EXPECT_TRUE(roster_effect.strobe_columns());
  }
  reset_effect_globals();
  pin_frame_clock(0);
  RingSpin<SMALL_W, SMALL_H> firmware_effect(
      RingSpin<SMALL_W, SMALL_H>::Strobe::DISABLED);
  HS_EXPECT_FALSE(firmware_effect.strobe_columns());
  firmware_effect.init();
  for (int frame = 0; frame < 2; ++frame) {
    pin_frame_clock(frame);
    firmware_effect.draw_frame();
    firmware_effect.advance_display();
  }
  HS_EXPECT_EQ(expected.size(), size_t{SMALL_W * SMALL_H});
  for (int y = 0; y < SMALL_H; ++y)
    for (int x = 0; x < SMALL_W; ++x) {
      HS_EXPECT_EQ(firmware_effect.get_pixel(x, y), expected[y * SMALL_W + x]);
    }
}

/**
 * @brief Verifies every lit RingSpin pixel sits on one of its great circles.
 * @details On the first frame, where the trail holds a single orientation
 *          step, every lit pixel must lie within the stroke half-width of some
 *          ring's equator. Alpha 0 must then blank the frame.
 */
inline void test_ringspin_trail_hugs_its_great_circles() {
  reset_effect_globals();
  pin_frame_clock(0);
  using WB = RingSpinWhiteBox;
  WB::RS fx;
  fx.init();
  pin_frame_clock(1);
  fx.draw_frame();
  fx.advance_display();

  // Polar deviation from the equator that the stroke half-width admits.
  const float band = std::sin(WB::head_half_width(fx));
  int lit = 0, off_circle = 0;
  float worst = 0.0f;
  for (int y = 0; y < DEFAULT_H; ++y)
    for (int x = 0; x < DEFAULT_W; ++x) {
      if (is_black(fx.get_pixel(x, y)))
        continue;
      ++lit;
      const math::Vector v = math::pixel_to_vector<DEFAULT_W, DEFAULT_H>(x, y);
      float nearest = 1.0f;
      for (int i = 0; i < WB::num_rings(); ++i)
        for (int s = 0; s < WB::substeps(fx, i); ++s)
          nearest =
              std::min(nearest, std::abs(math::dot(v, WB::axis(fx, i, s))));
      worst = hs_test::fold_worst(worst, nearest);
      if (nearest > band)
        ++off_circle;
    }

  std::printf("  [info] RingSpin: lit=%d off-circle=%d worst |dot|=%.6f band "
              "%.6f\n",
              lit, off_circle, worst, band);
  HS_EXPECT_GT(lit, 0);
  HS_EXPECT_EQ(off_circle, 0);

  HS_EXPECT_EQ(fx.updateParameter("Alpha", 0.0f), ParamSetResult::APPLIED);
  fx.draw_frame();
  fx.advance_display();
  HS_EXPECT_EQ((count_lit_region<DEFAULT_W, DEFAULT_H>(fx)), size_t{0});
}

// PetalFlow spawn-gap, hue-cursor and ring-pool bounds on every frame.
