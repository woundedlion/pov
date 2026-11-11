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
  const auto rotated =
      Pipeline::NodeAt<0>::run(Pullback::SphereSample{view, 0.0f}, frame.ctx,
                               std::get<0>(frame.prepared));
  const auto rippled =
      Pipeline::NodeAt<1>::run(rotated, frame.ctx, std::get<1>(frame.prepared));
  const auto projected =
      Pipeline::NodeAt<2>::run(rippled, frame.ctx, std::get<2>(frame.prepared));
  const auto sampled = Pipeline::NodeAt<3>::run(projected, frame.ctx,
                                                std::get<3>(frame.prepared));
  const auto covered =
      Pipeline::NodeAt<4>::run(sampled, frame.ctx, std::get<4>(frame.prepared));
  Pullback::ColorSample color;
  const bool ready = Pipeline::NodeAt<5>::run_shared_base(
      covered, frame.ctx, std::get<5>(frame.prepared), color);
  HS_EXPECT_TRUE(ready);
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
      const auto shader = FX::prepare_column_shader(*frame, x);
      if (!shader.shared)
        continue;
      for (int y = 0; y < 144; y += 18) {
        const auto view = direction(x, math::TrigLUT<288, 144>::sin_phi[y],
                                    math::TrigLUT<288, 144>::cos_phi[y]);
        const auto expected = FX::shade(view, *frame);
        const auto actual = unfactored_shared(view, *frame);
        HS_EXPECT_EQ(actual.color, expected.color);
        HS_EXPECT_EQ(actual.alpha, expected.alpha);
        const auto column = shader(math::TrigLUT<288, 144>::sin_phi[y],
                                   math::TrigLUT<288, 144>::cos_phi[y]);
        HS_EXPECT_EQ(column.color, actual.color);
        HS_EXPECT_EQ(column.alpha, actual.alpha);
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
      const auto shader = FX::prepare_column_shader(*frame, x);
      if (!shader.shared)
        continue;
      for (int y = 0; y < 144; ++y) {
        const float sine = math::TrigLUT<288, 144>::sin_phi[y];
        const float cosine = math::TrigLUT<288, 144>::cos_phi[y];
        const auto view = direction(x, sine, cosine);
        const auto expected_value = FX::shade(view, *frame);
        const auto actual_value = shader(sine, cosine);
        mismatches += actual_value.color != expected_value.color ||
                      bits(actual_value.alpha) != bits(expected_value.alpha);
        const Pixel expected_premultiplied_value =
            expected_value.color * expected_value.alpha;
        const Pixel actual_premultiplied_value =
            actual_value.color * actual_value.alpha;
        const auto check = [&](int expected, int actual,
                               int expected_premultiplied,
                               int actual_premultiplied) {
          const int error = std::abs(expected - actual);
          maximum_rgb = std::max(maximum_rgb, error);
          maximum_premultiplied =
              std::max(maximum_premultiplied,
                       std::abs(expected_premultiplied - actual_premultiplied));
          sum += error;
        };
        check(expected_value.color.r, actual_value.color.r,
              expected_premultiplied_value.r, actual_premultiplied_value.r);
        check(expected_value.color.g, actual_value.color.g,
              expected_premultiplied_value.g, actual_premultiplied_value.g);
        check(expected_value.color.b, actual_value.color.b,
              expected_premultiplied_value.b, actual_premultiplied_value.b);
        maximum_alpha = std::max(
            maximum_alpha, fabsf(expected_value.alpha - actual_value.alpha));
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
      const auto shader = FX::prepare_column_shader(frame, x);
      HS_EXPECT_FALSE(shader.shared);
      for (int y = 0; y < 144; ++y) {
        const float sine = math::TrigLUT<288, 144>::sin_phi[y];
        const float cosine = math::TrigLUT<288, 144>::cos_phi[y];
        const auto expected = FX::shade(direction(x, sine, cosine), frame);
        const auto actual = shader(sine, cosine);
        HS_EXPECT_TRUE(actual.color == expected.color);
        HS_EXPECT_EQ(bits(actual.alpha), bits(expected.alpha));
      }
    }
  }
}
