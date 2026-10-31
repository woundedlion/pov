/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file ashcloud_shared_color.h
 * @brief AshCloud shared-column dispatch and camera qualification. */

namespace ashcloud_shared_tests {
using FX = AshCloud<288, 144>;
using Pipeline = FX::RenderPipeline;

inline Color4 unfactored_shared(const math::Vector &view,
                                const FX::Frame &frame) {
  const auto ROTATED =
      Pipeline::NodeAt<0>::run(Pullback::SphereSample{view, 0.0f}, frame.ctx,
                               std::get<0>(frame.prepared));
  const auto RIPPLED =
      Pipeline::NodeAt<1>::run(ROTATED, frame.ctx, std::get<1>(frame.prepared));
  const auto PROJECTED =
      Pipeline::NodeAt<2>::run(RIPPLED, frame.ctx, std::get<2>(frame.prepared));
  const auto SAMPLED = Pipeline::NodeAt<3>::run(PROJECTED, frame.ctx,
                                                std::get<3>(frame.prepared));
  const auto COVERED =
      Pipeline::NodeAt<4>::run(SAMPLED, frame.ctx, std::get<4>(frame.prepared));
  Pullback::ColorSample color;
  const bool READY = Pipeline::NodeAt<5>::run_shared_base(
      COVERED, frame.ctx, std::get<5>(frame.prepared), color);
  HS_EXPECT_TRUE(READY);
  Pipeline::NodeAt<6>::run_phase(color, frame.ctx, std::get<6>(frame.prepared));
  Pipeline::NodeAt<7>::run_phase(color, frame.ctx, std::get<7>(frame.prepared));
  return Pullback::finish_color_sample(color);
}

inline math::Vector direction(int x, float sine, float cosine) {
  return {sine * math::TrigLUT<288, 144>::cos_theta(x), cosine,
          sine * math::TrigLUT<288, 144>::sin_theta[x]};
}
} // namespace ashcloud_shared_tests

/** @brief The shared color path is exact before camera factoring. */
inline void test_ash_cloud_shared_color_identity() {
  using namespace ashcloud_shared_tests;
  reset_globals();
  pin_frame_clock(0);
  FX effect;
  effect.init();
  for (int frame_index = 1; frame_index <= 16; ++frame_index) {
    pin_frame_clock(frame_index);
    effect.draw_frame();
    effect.advance_display();
    const auto *frame = ComposedFrameWhiteBox::prepared_frame(effect);
    HS_EXPECT_TRUE(frame != nullptr);
    if (!frame)
      return;
    for (int x = 0; x < 288; x += 18) {
      const auto SHADER = FX::prepare_column_shader(*frame, x);
      if (!SHADER.shared)
        continue;
      for (int y = 0; y < 144; y += 18) {
        const auto VIEW = direction(x, math::TrigLUT<288, 144>::sin_phi[y],
                                    math::TrigLUT<288, 144>::cos_phi[y]);
        const auto EXPECTED = FX::shade(VIEW, *frame);
        const auto ACTUAL = unfactored_shared(VIEW, *frame);
        HS_EXPECT_TRUE(ACTUAL.color == EXPECTED.color);
        HS_EXPECT_EQ(bits(ACTUAL.alpha), bits(EXPECTED.alpha));
      }
    }
  }
}

/** @brief Factored camera columns stay within their independent numeric qualification. */
inline void test_ash_cloud_shared_camera_qualification() {
  using namespace ashcloud_shared_tests;
  reset_globals();
  pin_frame_clock(0);
  FX effect;
  effect.init();
  uint64_t samples = 0, sum = 0, factored_samples = 0,
           unfactored_mismatches = 0;
  int maximum_rgb = 0, maximum_premultiplied = 0;
  float maximum_alpha = 0.0f, maximum_direction = 0.0f;
  for (int frame_index = 1; frame_index <= 600; ++frame_index) {
    pin_frame_clock(frame_index);
    effect.draw_frame();
    effect.advance_display();
    const auto *frame = ComposedFrameWhiteBox::prepared_frame(effect);
    HS_EXPECT_TRUE(frame != nullptr);
    if (!frame)
      return;
    for (int x = 0; x < 288; x += 18) {
      const auto SHADER = FX::prepare_column_shader(*frame, x);
      if (!SHADER.shared)
        continue;
      for (int y = 0; y < 144; ++y) {
        const float SINE = math::TrigLUT<288, 144>::sin_phi[y];
        const float COSINE = math::TrigLUT<288, 144>::cos_phi[y];
        const auto VIEW = direction(x, SINE, COSINE);
        const auto EXPECTED = FX::shade(VIEW, *frame);
        const auto ACTUAL = SHADER(SINE, COSINE);
        factored_samples += SHADER.factored;
        if (!SHADER.factored)
          unfactored_mismatches += ACTUAL.color != EXPECTED.color ||
                                   bits(ACTUAL.alpha) != bits(EXPECTED.alpha);
        const Pixel EXPECTED_PREMULTIPLIED = EXPECTED.color * EXPECTED.alpha;
        const Pixel ACTUAL_PREMULTIPLIED = ACTUAL.color * ACTUAL.alpha;
        const auto check = [&](int expected, int actual,
                               int expected_premultiplied,
                               int actual_premultiplied) {
          const int ERROR = std::abs(expected - actual);
          maximum_rgb = std::max(maximum_rgb, ERROR);
          maximum_premultiplied =
              std::max(maximum_premultiplied,
                       std::abs(expected_premultiplied - actual_premultiplied));
          sum += ERROR;
        };
        check(EXPECTED.color.r, ACTUAL.color.r, EXPECTED_PREMULTIPLIED.r,
              ACTUAL_PREMULTIPLIED.r);
        check(EXPECTED.color.g, ACTUAL.color.g, EXPECTED_PREMULTIPLIED.g,
              ACTUAL_PREMULTIPLIED.g);
        check(EXPECTED.color.b, ACTUAL.color.b, EXPECTED_PREMULTIPLIED.b,
              ACTUAL_PREMULTIPLIED.b);
        maximum_alpha =
            std::max(maximum_alpha, fabsf(EXPECTED.alpha - ACTUAL.alpha));
        const auto ROTATED =
            Pipeline::NodeAt<0>::run(Pullback::SphereSample{VIEW, 0.0f},
                                     frame->ctx, std::get<0>(frame->prepared));
        const auto ACTUAL_ROTATED =
            SHADER.factored ? SINE * SHADER.equator + COSINE * SHADER.up
                            : ROTATED.dir;
        const auto DELTA = ACTUAL_ROTATED - ROTATED.dir;
        maximum_direction = std::max(maximum_direction, DELTA.length());
        ++samples;
      }
    }
  }
  HS_EXPECT_EQ(samples, uint64_t{599 * 16 * 144});
  HS_EXPECT_EQ(factored_samples, FX::camera_column_admitted() ? samples : 0);
  HS_EXPECT_EQ(unfactored_mismatches, uint64_t{0});
  if (!FX::camera_column_admitted()) {
    HS_EXPECT_EQ(maximum_rgb, 0);
    HS_EXPECT_EQ(maximum_premultiplied, 0);
    HS_EXPECT_EQ(maximum_alpha, 0.0f);
    HS_EXPECT_EQ(maximum_direction, 0.0f);
  }
  HS_EXPECT_LE(maximum_rgb, 256);
  HS_EXPECT_LE(maximum_premultiplied, 128);
  HS_EXPECT_LE(maximum_alpha, .002f);
  HS_EXPECT_LE(maximum_direction, 1e-6f);
  std::printf(
      "  AshCloud camera samples=%llu RGB max=%d mean=%.9g premultiplied=%d alpha=%.9g direction=%.9g\n",
      static_cast<unsigned long long>(samples), maximum_rgb,
      double(sum) / (samples * 3), maximum_premultiplied, double(maximum_alpha),
      double(maximum_direction));
}

/** @brief Endpoints and unready color caches select the canonical shader exactly. */
inline void test_ash_cloud_shared_column_fallbacks() {
  using namespace ashcloud_shared_tests;
  reset_globals();
  pin_frame_clock(0);
  FX effect;
  effect.init();
  for (int frame_index = 1; frame_index <= 2; ++frame_index) {
    pin_frame_clock(frame_index);
    effect.draw_frame();
    effect.advance_display();
  }
  const auto *current = ComposedFrameWhiteBox::prepared_frame(effect);
  HS_EXPECT_TRUE(current != nullptr);
  if (!current)
    return;
  HS_EXPECT_TRUE(FX::prepare_column_shader(*current, 0).shared);
  for (int mode = 0; mode < 7; ++mode) {
    FX::Frame frame(current->ctx);
    auto &base = std::get<5>(frame.prepared);
    if (mode == 0)
      frame.ctx.temporal.fraction = base.fraction = 0.0f;
    if (mode == 1)
      frame.ctx.temporal.fraction = base.fraction = 1.0f;
    if (mode == 2)
      base.lab.low = nullptr;
    if (mode == 3)
      base.lab.rgb_low = nullptr;
    if (mode == 4)
      base.lab.rgb_high = nullptr;
    if (mode == 5)
      std::get<6>(frame.prepared).ready = false;
    if (mode == 6)
      std::get<7>(frame.prepared).ready = false;
    for (int x = 0; x < 288; x += 18) {
      const auto SHADER = FX::prepare_column_shader(frame, x);
      HS_EXPECT_FALSE(SHADER.shared);
      for (int y = 0; y < 144; ++y) {
        const float SINE = math::TrigLUT<288, 144>::sin_phi[y];
        const float COSINE = math::TrigLUT<288, 144>::cos_phi[y];
        const auto EXPECTED = FX::shade(direction(x, SINE, COSINE), frame);
        const auto ACTUAL = SHADER(SINE, COSINE);
        HS_EXPECT_TRUE(ACTUAL.color == EXPECTED.color);
        HS_EXPECT_EQ(bits(ACTUAL.alpha), bits(EXPECTED.alpha));
      }
    }
  }
}

/** @brief Camera factoring admits only the qualified physical resolution and endpoints. */
inline void test_ash_cloud_camera_column_admission() {
  using namespace ashcloud_shared_tests;
  HS_EXPECT_FALSE((AshCloud<96, 20>::camera_column_admitted()));
  math::init_geometry_luts<288, 144>();
  const bool PHYSICAL =
      math::DISPLAY_NORTH_PHI == .02f * math::PI_F &&
      math::DISPLAY_SOUTH_PHI == .98f * math::PI_F &&
      math::DisplayGeometry<144>::NORTH_PHI == .02f * math::PI_F &&
      math::DisplayGeometry<144>::SOUTH_PHI == .98f * math::PI_F;
  HS_EXPECT_EQ(FX::camera_column_admitted(), PHYSICAL);
#if HS_RUNTIME_DISPLAY_GEOMETRY && !defined(HS_TEST_H_OFFSET)
  const float SAVED_NORTH = math::DISPLAY_NORTH_PHI;
  const float SAVED_SOUTH = math::DISPLAY_SOUTH_PHI;
  for (const auto ENDPOINTS :
       {std::pair{0.0f, math::PI_F},
        std::pair{.02f * math::PI_F, .98f * math::PI_F},
        std::pair{.01f * math::PI_F, .99f * math::PI_F},
        std::pair{std::nextafter(.02f * math::PI_F, math::PI_F),
                  .98f * math::PI_F}}) {
    HS_EXPECT_TRUE(
        math::set_display_geometry(ENDPOINTS.first, ENDPOINTS.second));
    const bool EXPECTED = ENDPOINTS.first == .02f * math::PI_F &&
                          ENDPOINTS.second == .98f * math::PI_F;
    if (!EXPECTED)
      HS_EXPECT_FALSE(FX::camera_column_admitted());
    math::init_geometry_luts<288, 144>();
    HS_EXPECT_EQ(FX::camera_column_admitted(), EXPECTED);
  }
  HS_EXPECT_TRUE(math::set_display_geometry(SAVED_NORTH, SAVED_SOUTH));
  math::init_geometry_luts<288, 144>();
#endif
}
