/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file temporal_sampling.h
 * @brief Column-time interpolation and retained-resource tests.
 */

/** @brief Both arms use absolute column timing through clock rollover. */
inline void test_temporal_clock_sampling() {
  HS_CONTEXT("temporal clock sampling");
  for (int x = 0; x < 144; ++x) {
    HS_EXPECT_NEAR(Pullback::sweep_fraction(x, 288), x / 144.0f, 0.0f);
    HS_EXPECT_NEAR(Pullback::sweep_fraction(x + 144, 288), x / 144.0f, 0.0f);
  }
  HS_EXPECT_NEAR(Pullback::sample_clock(0.99f, 0.01f, 0.5f, 1.0f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(Pullback::sample_clock(0.01f, 0.99f, 0.5f, 1.0f), 0.0f, 1e-6f);
  const float LAST = Pullback::sample_clock(0.98f, 0.0f, 143.0f / 144, 1.0f);
  HS_EXPECT_NEAR(1.0f - LAST, 0.02f / 144, 1e-6f);

  using FX = GridSpace<96, 20>;
  typename FX::FrameState from{}, to{};
  from.params = FX::initial_params();
  to.params = from.params;
  from.resources.template get<"outer_warp">().phase = 0.99f;
  to.resources.template get<"outer_warp">().phase = 0.01f;
  to.outer_conjugate = math::make_rotation(math::Y_AXIS, 0.4f);
  const auto middle = Pullback::sample_frame(from, to, 0.5f);
  HS_EXPECT_NEAR(middle.resources.template get<"outer_warp">().phase, 1.0f,
                 1e-6f);
  const math::Quaternion EXPECTED = math::make_rotation(math::Y_AXIS, 0.2f);
  HS_EXPECT_NEAR(math::dot(middle.outer_conjugate, EXPECTED), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(Pullback::sample_frame(from, to, 0.0f)
                     .resources.template get<"outer_warp">()
                     .phase,
                 0.99f, 0.0f);
  HS_EXPECT_NEAR(Pullback::sample_frame(from, to, 1.0f)
                     .resources.template get<"outer_warp">()
                     .phase,
                 0.01f, 0.0f);
}

/** @brief Retained color tables survive subsequent bakes. */
inline void test_temporal_owned_resources() {
  HS_CONTEXT("temporal owned color resources");
  reset_globals();
  using FX = GridSpace<96, 20>;
  FX effect;
  effect.init();
  effect.draw_frame();
  effect.advance_display();
  const auto current = effect.frame_for_test();
  Pullback::TemporalHistory<FX::FrameState, Pullback::HueMode::NOISE, true>
      history;
  history.capture(current, persistent_arena);
  HS_EXPECT_TRUE(history.previous.palette != current.palette);
  HS_EXPECT_TRUE(history.previous.hue_noise_lut != current.hue_noise_lut);
  HS_EXPECT_TRUE(history.previous.hue_rotation_lut != current.hue_rotation_lut);
  const Pixel PALETTE = history.previous.palette->get(0.4f).color;
  const auto NOISE = history.hue_noise;
  const auto ROTATION = history.hue_rotation;
  for (int frame = 0; frame < 4; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_TRUE(history.previous.palette->get(0.4f).color == PALETTE);
  HS_EXPECT_TRUE(history.hue_noise == NOISE);
  HS_EXPECT_TRUE(history.hue_rotation == ROTATION);
}

/** @brief GridSpace shades each clipped column at its own sweep time. */
inline void test_gridspace_temporal_pixels() {
  HS_CONTEXT("GridSpace column-time pixels");
  reset_globals();
  constexpr int W = 96;
  constexpr int H = 20;
  using FX = GridSpace<W, H>;
  FX effect;
  effect.init();
  effect.draw_frame();
  effect.advance_display();
  Pullback::TemporalHistory<FX::FrameState, Pullback::HueMode::NOISE, true>
      history;
  history.capture(effect.frame_for_test(), persistent_arena);
  effect.set_clip(3, 17, 57, 83);
  effect.draw_frame();
  effect.advance_display();
  const auto current = effect.frame_for_test();
  const auto frozen = FX::RenderPipeline::prepare(current);
  int changed = 0;
  for (int x = 57; x < 83; ++x) {
    const float TIME = static_cast<float>(x - W / 2) / (W / 2);
    const auto frame = FX::RenderPipeline::prepare(
        Pullback::sample_frame(history.previous, current, TIME));
    for (int y = 3; y < 17; ++y) {
      const float SIN_PHI = math::TrigLUT<W, H>::sin_phi[y];
      const math::Vector VIEW(SIN_PHI * math::TrigLUT<W, H>::cos_theta(x),
                              math::TrigLUT<W, H>::cos_phi[y],
                              SIN_PHI * math::TrigLUT<W, H>::sin_theta[x]);
      const Color4 EXPECTED = FX::shade(VIEW, frame);
      const Pixel PIXEL = EXPECTED.color * EXPECTED.alpha;
      HS_EXPECT_TRUE(effect.get_pixel(x, y) == PIXEL);
      const Color4 CURRENT = FX::shade(VIEW, frozen);
      if (PIXEL != CURRENT.color * CURRENT.alpha)
        ++changed;
    }
  }
  HS_EXPECT_TRUE(changed > 0);
}
