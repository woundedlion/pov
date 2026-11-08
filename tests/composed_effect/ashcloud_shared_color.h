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
        const auto COLUMN = SHADER(math::TrigLUT<288, 144>::sin_phi[y],
                                   math::TrigLUT<288, 144>::cos_phi[y]);
#if HS_TEST_FAST_MATH
        HS_EXPECT_LE(std::abs(int(COLUMN.color.r) - int(ACTUAL.color.r)), 1);
        HS_EXPECT_LE(std::abs(int(COLUMN.color.g) - int(ACTUAL.color.g)), 1);
        HS_EXPECT_LE(std::abs(int(COLUMN.color.b) - int(ACTUAL.color.b)), 1);
#else
        HS_EXPECT_TRUE(COLUMN.color == ACTUAL.color);
#endif
        HS_EXPECT_EQ(bits(COLUMN.alpha), bits(ACTUAL.alpha));
      }
    }
  }
}

/** @brief Most columns take the shared path, and its colors match canonical
 * colors within each math tier. */
inline void test_ash_cloud_shared_column_qualification() {
  using namespace ashcloud_shared_tests;
  constexpr int FRAMES = 600, COLUMN_STRIDE = 18;
  reset_globals();
  pin_frame_clock(0);
  FX effect;
  effect.init();
  uint64_t samples = 0, sum = 0, mismatches = 0;
  int maximum_rgb = 0, maximum_premultiplied = 0;
  float maximum_alpha = 0.0f;
  for (int frame_index = 1; frame_index <= FRAMES; ++frame_index) {
    pin_frame_clock(frame_index);
    effect.draw_frame();
    effect.advance_display();
    const auto *frame = ComposedFrameWhiteBox::prepared_frame(effect);
    HS_EXPECT_TRUE(frame != nullptr);
    if (!frame)
      return;
    for (int x = 0; x < 288; x += COLUMN_STRIDE) {
      const auto SHADER = FX::prepare_column_shader(*frame, x);
      if (!SHADER.shared)
        continue;
      for (int y = 0; y < 144; ++y) {
        const float SINE = math::TrigLUT<288, 144>::sin_phi[y];
        const float COSINE = math::TrigLUT<288, 144>::cos_phi[y];
        const auto VIEW = direction(x, SINE, COSINE);
        const auto EXPECTED = FX::shade(VIEW, *frame);
        const auto ACTUAL = SHADER(SINE, COSINE);
        mismatches += ACTUAL.color != EXPECTED.color ||
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
        ++samples;
      }
    }
  }
  constexpr uint64_t SAMPLED = uint64_t{FRAMES} * (288 / COLUMN_STRIDE) * 144;
  HS_EXPECT_GE(samples, SAMPLED * 9 / 10);
#if HS_TEST_FAST_MATH
  HS_EXPECT_LE(maximum_rgb, 16);
  HS_EXPECT_LE(maximum_premultiplied, 8);
  HS_EXPECT_LE(maximum_alpha, 1e-5f);
#else
  HS_EXPECT_EQ(mismatches, uint64_t{0});
  HS_EXPECT_EQ(maximum_rgb, 0);
  HS_EXPECT_EQ(maximum_premultiplied, 0);
  HS_EXPECT_EQ(maximum_alpha, 0.0f);
#endif
  std::printf(
      "  AshCloud shared samples=%llu mismatches=%llu RGB max=%d mean=%.9g premultiplied=%d alpha=%.9g\n",
      static_cast<unsigned long long>(samples),
      static_cast<unsigned long long>(mismatches), maximum_rgb,
      double(sum) / (samples * 3), maximum_premultiplied,
      double(maximum_alpha));
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
