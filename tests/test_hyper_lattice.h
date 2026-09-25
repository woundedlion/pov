/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "effects/HyperLattice.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace hyper_lattice_tests {

namespace HL = HyperLatticeDetail;

struct HyperLatticeWhiteBox {
  using Effect = HyperLattice<96, 20>;

  static math::Vec4 origin(const Effect &effect) { return effect.origin; }
  static std::array<float, 6> rotation_phase(const Effect &effect) {
    return effect.rotation_phase;
  }
  static HL::Params &params(Effect &effect) { return effect.params; }
  static void step_depth_palette(Effect &effect) {
    effect.depth_palette.step();
  }
  static bool depth_palette_fading(const Effect &effect) {
    return effect.depth_palette.fading();
  }
  static Pixel depth_color(const Effect &effect, float amount) {
    return effect.depth_palette.palette().get(amount).color;
  }
  static Pixel axis_color(const Effect &effect, float amount) {
    return effect.axis_palette.get(amount).color;
  }
  static const BakedPalette *depth_palette(const Effect &effect) {
    return &effect.depth_palette.palette();
  }
  static const BakedPalette *axis_palette(const Effect &effect) {
    return &effect.axis_palette.view();
  }
};

inline void test_periodic_distance() {
  HS_EXPECT_EQ(HL::periodic_distance(0.0f), 0.0f);
  HS_EXPECT_EQ(HL::periodic_distance(1.0f), 0.0f);
  HS_EXPECT_NEAR(HL::periodic_distance(-1.25f), 0.25f, 1e-6f);
  HS_EXPECT_EQ(HL::periodic_distance(-0.5f), 0.5f);
  HS_EXPECT_EQ(HL::periodic_distance(1.5f), 0.5f);
  HS_EXPECT_EQ(HL::periodic_distance(2.5f), 0.5f);
  HS_EXPECT_NEAR(HL::periodic_distance(3.5f), 0.5f, 1e-6f);
}

inline void test_edge_metrics() {
  const math::Vec4 origin{{0.0f, 0.0f, 0.0f, 0.0f}};
  const math::Vec4 direction{{0.0f, 0.02f, 0.4f, 0.1f}};
  const HL::EdgeMetric cubic =
      HL::edge_metric_3d_at(origin, direction, 0, 1.0f);
  HS_EXPECT_NEAR(cubic.distance_sq, 0.0004f, 1e-7f);
  HS_EXPECT_EQ(cubic.free_axis, uint8_t(2));

  HL::EdgeMetric hyper{};
  HS_EXPECT_TRUE(HL::edge_metric_4d_at_bounded(origin, direction, 0, 1.0f, 0.5f,
                                               0.25f, hyper));
  HS_EXPECT_NEAR(hyper.distance_sq, 0.0104f, 1e-6f);
  HS_EXPECT_EQ(hyper.free_axis, uint8_t(2));

  // Neither in-plane component inside the band.
  HS_EXPECT_FALSE(HL::edge_metric_4d_at_bounded(origin, direction, 0, 1.0f,
                                                0.01f, 0.25f, hyper));
  // One in-plane component inside the band, the w component outside it.
  HS_EXPECT_FALSE(HL::edge_metric_4d_at_bounded(origin, direction, 0, 1.0f,
                                                0.05f, 0.25f, hyper));
  // Inside the band, beyond the squared-distance limit.
  HS_EXPECT_FALSE(HL::edge_metric_4d_at_bounded(origin, direction, 0, 1.0f,
                                                0.5f, 1.0e-6f, hyper));

  HL::EdgeMetric no_axis{};
  HS_EXPECT_TRUE(HL::edge_metric_4d_at_bounded<false>(
      origin, direction, 0, 1.0f, 0.5f, 0.25f, no_axis));
  HS_EXPECT_NEAR(no_axis.distance_sq, 0.0104f, 1e-6f);
  HS_EXPECT_EQ(no_axis.free_axis, uint8_t(0));
}

inline void test_so4_rotation() {
  HL::FrameState frame{};
  frame.params.mode = HL::LatticeMode::FOUR_D_SLICE;
  frame.params.far_distance = 8.0f;
  frame.rotation_phase[3] = 0.5f * math::PI_F;
  const HL::PreparedTrace prepared = HL::prepare_trace(frame);
  const math::Vec4 rotated =
      prepared.world_to_lattice.apply({{1.0f, 0.0f, 0.0f, 0.0f}});
  HS_EXPECT_NEAR(rotated[0], 0.0f, 2e-4f);
  HS_EXPECT_NEAR(rotated[3], 1.0f, 2e-4f);
  float norm_sq = 0.0f;
  for (int axis = 0; axis < HL::DIMENSIONS; ++axis)
    norm_sq += rotated[axis] * rotated[axis];
  HS_EXPECT_NEAR(norm_sq, 1.0f, 2e-4f);

  frame.params.mode = HL::LatticeMode::THREE_D;
  const math::Vec4 cubic = HL::prepare_trace(frame).world_to_lattice.apply(
      {{1.0f, 0.0f, 0.0f, 0.0f}});
  HS_EXPECT_EQ(cubic[3], 0.0f);
}

inline void test_dimensional_rotation_wrap_is_continuous() {
  HL::FrameState frame{};
  frame.params.mode = HL::LatticeMode::DIMENSIONAL_RIFT;
  frame.rotation_phase[3] = math::TWO_PI_F - 1.0e-4f;
  const math::Vec4 before = HL::prepare_trace(frame).world_to_lattice.apply(
      {{1.0f, 0.0f, 0.0f, 0.0f}});
  frame.rotation_phase[3] = 0.0f;
  const math::Vec4 after = HL::prepare_trace(frame).world_to_lattice.apply(
      {{1.0f, 0.0f, 0.0f, 0.0f}});
  for (int axis = 0; axis < HL::DIMENSIONS; ++axis)
    HS_EXPECT_NEAR(before[axis], after[axis], 2.0e-4f);
}

inline void test_resolution_aware_wire_coverage() {
  constexpr float LOW_RES = HL::pixel_half_angle<96, 20>();
  constexpr float HIGH_RES = HL::pixel_half_angle<288, 144>();
  static_assert(LOW_RES > HIGH_RES);

  HL::FrameState frame{};
  frame.pixel_half_angle = LOW_RES;
  const HL::PreparedTrace low = HL::prepare_trace(frame);
  frame.pixel_half_angle = HIGH_RES;
  const HL::PreparedTrace high = HL::prepare_trace(frame);
  HS_EXPECT_GT(low.aa_scale, high.aa_scale);
  HS_EXPECT_NEAR(low.aa_scale / high.aa_scale, LOW_RES / HIGH_RES, 1e-5f);

  constexpr float RADIUS = 0.1f;
  constexpr float HALF_WIDTH = 0.02f;
  constexpr float OFFSET = 0.01f;
  const float inside = HL::wire_coverage((RADIUS - OFFSET) * (RADIUS - OFFSET),
                                         RADIUS, HALF_WIDTH);
  const float boundary = HL::wire_coverage(RADIUS * RADIUS, RADIUS, HALF_WIDTH);
  const float outside = HL::wire_coverage((RADIUS + OFFSET) * (RADIUS + OFFSET),
                                          RADIUS, HALF_WIDTH);
  HS_EXPECT_NEAR(boundary, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(inside + outside, 1.0f, 1e-5f);

  constexpr float DISTANCE = 1.0f;
  const float low_half_width = HALF_WIDTH + low.aa_scale * DISTANCE;
  const float high_half_width = HALF_WIDTH + high.aa_scale * DISTANCE;
  HS_EXPECT_GT(low_half_width, high_half_width);
  const float low_outside = HL::wire_coverage(
      (RADIUS + OFFSET) * (RADIUS + OFFSET), RADIUS, low_half_width);
  const float high_outside = HL::wire_coverage(
      (RADIUS + OFFSET) * (RADIUS + OFFSET), RADIUS, high_half_width);
  HS_EXPECT_GT(low_outside, high_outside);
  HS_EXPECT_GT(high_outside, outside);
}

inline void test_near_field_fade() {
  constexpr float RADIUS = 0.1f;
  constexpr float NEAR_START = 1.5f * RADIUS;
  constexpr float NEAR_INV_SPAN = 1.0f / (2.5f * RADIUS);
  HS_EXPECT_EQ(HL::near_field_coverage(0.15f, NEAR_START, NEAR_INV_SPAN), 0.0f);
  HS_EXPECT_GT(HL::near_field_coverage(0.275f, NEAR_START, NEAR_INV_SPAN),
               0.0f);
  HS_EXPECT_LT(HL::near_field_coverage(0.275f, NEAR_START, NEAR_INV_SPAN),
               1.0f);
  HS_EXPECT_EQ(HL::near_field_coverage(0.4f, NEAR_START, NEAR_INV_SPAN), 1.0f);
}

inline void test_far_shell_fade() {
  HS_EXPECT_EQ(HL::shell_horizon_coverage(0, 2, 0.75f, 1.0f), 1.0f);
  HS_EXPECT_EQ(HL::shell_horizon_coverage(1, 2, 1.0f, 1.0f), 1.0f);
  HS_EXPECT_NEAR(HL::shell_horizon_coverage(1, 2, 1.5f, 1.0f), 0.5f, 1e-6f);
  HS_EXPECT_EQ(HL::shell_horizon_coverage(1, 2, 2.0f, 1.0f), 0.0f);
}

inline void test_pause_does_not_stop_motion() {
  reset_globals();
  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  effect.setAnimationsPaused(true);

  const math::Vec4 origin_before = HyperLatticeWhiteBox::origin(effect);
  const auto rotation_before = HyperLatticeWhiteBox::rotation_phase(effect);
  effect.draw_frame();
  effect.advance_display();
  const math::Vec4 origin_after = HyperLatticeWhiteBox::origin(effect);
  const auto rotation_after = HyperLatticeWhiteBox::rotation_phase(effect);
  HS_EXPECT_NE(origin_after[0], origin_before[0]);
  HS_EXPECT_NE(rotation_after[0], rotation_before[0]);

  HyperLatticeWhiteBox::params(effect).speed = 0.0f;
  const math::Vec4 stopped_before = HyperLatticeWhiteBox::origin(effect);
  effect.draw_frame();
  effect.advance_display();
  const math::Vec4 stopped_after = HyperLatticeWhiteBox::origin(effect);
  for (int axis = 0; axis < HL::DIMENSIONS; ++axis)
    HS_EXPECT_EQ(stopped_after[axis], stopped_before[axis]);
}

inline void test_depth_palette_mutates_slowly_while_paused() {
  reset_globals();
  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  effect.setAnimationsPaused(true);
  const Pixel initial = HyperLatticeWhiteBox::depth_color(effect, 0.5f);
  const Pixel axis = HyperLatticeWhiteBox::axis_color(effect, 0.5f);
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_TRUE(HyperLatticeWhiteBox::depth_palette_fading(effect));
  HS_EXPECT_EQ(HyperLatticeWhiteBox::depth_color(effect, 0.5f), initial);

  for (int frame = 0; frame < 480; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_NE(HyperLatticeWhiteBox::depth_color(effect, 0.5f), initial);
  HS_EXPECT_EQ(HyperLatticeWhiteBox::axis_color(effect, 0.5f), axis);
}

inline void test_depth_palette_keeps_cool_character() {
  reset_globals();
  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  const Pixel initial = HyperLatticeWhiteBox::depth_color(effect, 0.5f);
  for (int frame = 0; frame <= 2880; ++frame) {
    if (frame % 120 == 0) {
      float previous_lightness = 0.0f;
      for (int sample = 0; sample <= 16; ++sample) {
        const OKLCH color = pixel_to_oklch(
            HyperLatticeWhiteBox::depth_color(effect, sample / 16.0f));
        const float hue = math::wrap_t(color.h / math::TWO_PI_F) * 360.0f;
        HS_EXPECT_TRUE(hue >= 195.0f && hue <= 290.0f);
        HS_EXPECT_GT(color.L, previous_lightness);
        previous_lightness = color.L;
      }
    }
    if (frame < 2880)
      HyperLatticeWhiteBox::step_depth_palette(effect);
  }
  HS_EXPECT_EQ(HyperLatticeWhiteBox::depth_color(effect, 0.5f), initial);
}

inline void test_next_plane_is_strict() {
  HS_EXPECT_EQ(HL::next_plane_offset(0.0f, true), 1.0f);
  HS_EXPECT_EQ(HL::next_plane_offset(0.0f, false), 1.0f);
  HS_EXPECT_EQ(HL::next_plane_offset(0.25f, true), 0.75f);
  HS_EXPECT_EQ(HL::next_plane_offset(0.25f, false), 0.25f);
}

inline void test_trace_layers_are_front_to_back() {
  HL::FrameState frame{};
  frame.params = HyperLattice<96, 20>::preset_params(0);
  frame.origin = {{0.25f, 0.0f, 0.31f, 0.43f}};
  const HL::PreparedTrace prepared = HL::prepare_trace(frame);
  float previous_distance = 0.0f;
  int layers = 0;
  HL::trace_layers(math::X_AXIS, prepared, [&](const HL::TraceHit &hit) {
    HS_EXPECT_GT(hit.coverage, 0.0f);
    HS_EXPECT_LE(hit.coverage, 1.0f);
    HS_EXPECT_GT(hit.distance, previous_distance);
    previous_distance = hit.distance;
    ++layers;
    return true;
  });
  HS_EXPECT_EQ(layers, 2);
}

inline void test_layer_composite_reveals_background() {
  LayerComposite boundary;
  boundary.remaining = MIN_ENCODABLE_ALPHA;
  HS_EXPECT_FALSE(boundary.saturated());
  LayerComposite empty;
  HS_EXPECT_EQ(empty.finish().alpha, 0.0f);
  HS_EXPECT_FALSE(empty.saturated());
  empty.add(Pixel(60000, 0, 0), 0.5f);
  empty.add(Pixel(0, 60000, 0), 0.5f);
  const Color4 partial = empty.finish();
  HS_EXPECT_EQ(partial.alpha, 0.75f);
  HS_EXPECT_EQ(partial.color.r, 40000);
  HS_EXPECT_EQ(partial.color.g, 20000);
  HS_EXPECT_EQ(partial.color.b, 0);
  HS_EXPECT_FALSE(empty.saturated());
  empty.add(Pixel(0, 0, 60000), 1.0f);
  HS_EXPECT_TRUE(empty.saturated());
  LayerComposite composite;
  composite.add(Pixel(65535, 0, 0), 0.5f);
  composite.add(Pixel(0, 65535, 0), 1.0f);
  const Color4 result = composite.finish();
  HS_EXPECT_EQ(result.alpha, 1.0f);
  HS_EXPECT_NEAR(result.color.r, 32768, 1);
  HS_EXPECT_NEAR(result.color.g, 32768, 1);
  HS_EXPECT_EQ(result.color.b, 0);
}

inline void test_surface_origin_parallax() {
  HL::FrameState frame{};
  frame.params = HyperLattice<96, 20>::preset_params(0);
  frame.params.sphere_radius = 0.0f;
  frame.origin = {{0.25f, 0.0f, 0.31f, 0.43f}};
  const HL::TraceHit centered =
      HL::trace(math::X_AXIS, HL::prepare_trace(frame));
  frame.params.cell_size = 2.0f;
  const HL::PreparedTrace scaled_trace = HL::prepare_trace(frame);
  const HL::TraceHit scaled = HL::trace(math::X_AXIS, scaled_trace);
  frame.params.sphere_radius = 0.4f;
  const HL::TraceHit scaled_surface =
      HL::trace(math::X_AXIS, HL::prepare_trace(frame));
  frame.params.cell_size = 1.0f;
  const HL::TraceHit surfaced =
      HL::trace(math::X_AXIS, HL::prepare_trace(frame));
  HS_EXPECT_NEAR(centered.distance, 0.75f, 1e-6f);
  HS_EXPECT_NEAR(scaled.distance, 1.5f, 1e-6f);
  HS_EXPECT_NEAR(scaled_surface.distance, 0.7f, 1e-6f);
  HS_EXPECT_NEAR(scaled_trace.far_distance, frame.params.far_distance, 1e-6f);
  HS_EXPECT_NEAR(surfaced.distance, 0.35f, 1e-6f);
}

inline void test_hyperplane_event() {
  HL::FrameState frame{};
  frame.params = HyperLattice<96, 20>::preset_params(1);
  frame.params.sphere_radius = 0.4f;
  frame.origin = {{0.0f, 0.0f, 0.31f, 0.25f}};
  frame.rotation_phase[3] = 0.5f * math::PI_F;
  const HL::PreparedTrace prepared = HL::prepare_trace(frame);
  const HL::TraceHit hit = HL::trace(math::X_AXIS, prepared);
  HS_EXPECT_GT(hit.coverage, 0.0f);
  HS_EXPECT_EQ(hit.free_axis, uint8_t(2));
  HS_EXPECT_NEAR(hit.distance, 0.35f, 3e-4f);
}

inline void test_coincident_planes_form_one_layer() {
  HL::FrameState frame{};
  frame.params = HyperLattice<96, 20>::preset_params(0);
  frame.params.sphere_radius = 0.0f;
  frame.origin = {{0.25f, 0.25f, 0.31f, 0.43f}};
  const HL::PreparedTrace prepared = HL::prepare_trace(frame);
  constexpr float INV_SQRT_TWO = 0.707106781f;
  float previous_distance = 0.0f;
  int layers = 0;
  HL::trace_layers(math::Vector(INV_SQRT_TWO, INV_SQRT_TWO, 0.0f), prepared,
                   [&](const HL::TraceHit &hit) {
                     HS_EXPECT_GT(hit.distance, previous_distance);
                     previous_distance = hit.distance;
                     ++layers;
                     return true;
                   });
  HS_EXPECT_EQ(layers, 2);
}

/** @brief One shade() result in signature fold order: RGB then Q16 alpha. */
struct ShadeSample {
  uint16_t r;     /**< Linear red channel. */
  uint16_t g;     /**< Linear green channel. */
  uint16_t b;     /**< Linear blue channel. */
  uint16_t alpha; /**< Coverage as Q16. */
};

/**
 * @brief Per-channel slack, in 16-bit linear units, for a libm difference.
 * @details The palette lookup runs cbrtf/powf through the OKLab gamut search,
 * whose last bits differ between libm builds; sixteen linear units stay below
 * one sRGB code step at the dark end. Sized as in
 * tests/mindsplatter_palette_check.cpp, which holds the same path to the same
 * band.
 */
constexpr uint16_t MAX_SHADE_CHANNEL_DELTA = 16;

/**
 * @brief Folds a sample table into an FNV-1a 64 signature.
 * @param samples Table start.
 * @param count Row count.
 * @return The folded signature.
 */
inline uint64_t shade_signature(const ShadeSample *samples, size_t count) {
  uint64_t signature = hs_test::FNV1A64_BASIS;
  for (size_t row = 0; row < count; ++row) {
    signature = hs_test::fnv1a64_channel(signature, samples[row].r);
    signature = hs_test::fnv1a64_channel(signature, samples[row].g);
    signature = hs_test::fnv1a64_channel(signature, samples[row].b);
    signature = hs_test::fnv1a64_channel(signature, samples[row].alpha);
  }
  return signature;
}

/**
 * @brief Scores rendered samples against a golden table, signature first.
 * @param label Context label printed with every out-of-band channel.
 * @param rendered Samples in fold order, one row per (preset, direction).
 * @param golden Rows the pin was recorded from.
 * @param count Row count of both tables.
 * @param per_preset Rows per preset, so a failing row names preset and sample.
 * @param pin Signature the golden table folds to.
 * @details The signature is one comparison over the whole table and is all a
 * matching run pays. A run whose hash moved falls through to the per-channel
 * band, which names the preset, sample and channel that drifted rather than
 * printing two 64-bit numbers.
 */
inline void expect_shade_samples(const char *label, const ShadeSample *rendered,
                                 const ShadeSample *golden, size_t count,
                                 size_t per_preset, uint64_t pin) {
  const uint64_t golden_signature = shade_signature(golden, count);
  HS_EXPECT_EQ(golden_signature, pin);
  if (shade_signature(rendered, count) == golden_signature)
    return;
  for (size_t row = 0; row < count; ++row) {
    const ShadeSample &got = rendered[row];
    const ShadeSample &want = golden[row];
    HS_CONTEXT(label, static_cast<long long>(row / per_preset),
               static_cast<long long>(row % per_preset));
    HS_EXPECT_NEAR(got.r, want.r, MAX_SHADE_CHANNEL_DELTA);
    HS_EXPECT_NEAR(got.g, want.g, MAX_SHADE_CHANNEL_DELTA);
    HS_EXPECT_NEAR(got.b, want.b, MAX_SHADE_CHANNEL_DELTA);
    HS_EXPECT_NEAR(got.alpha, want.alpha, MAX_SHADE_CHANNEL_DELTA);
  }
}

/**
 * @brief Pins shade() over a fixed direction and preset sample.
 * @details GOLDEN is the oracle; the FNV-1a 64 signature over it is the
 * one-comparison pre-check a matching run pays. Provenance: no generator emits
 * either. Re-derive by printing the RGB and Q16 alpha of every sample from this
 * case built by the native clang test toolchain
 * (cmake/toolchain-native-clang.cmake) and pasting the table and its fold back.
 * Unlike test_specialized_render_signature(), this path takes its coverage ramp
 * through an IEEE division rather than fast_reciprocal()'s Newton step, so it
 * reproduces under the shipping -ffast-math -fno-finite-math-only pair the
 * fast-math CI leg builds this module with, and carries no skip.
 */
inline void test_render_signature() {
  reset_globals();
  static constexpr math::Vector DIRECTIONS[] = {
      {1.0f, 0.0f, 0.0f},
      {-1.0f, 0.0f, 0.0f},
      {0.0f, 1.0f, 0.0f},
      {0.0f, 0.0f, 1.0f},
      {0.577350269f, 0.577350269f, 0.577350269f},
      {-0.707106781f, 0.707106781f, 0.0f},
      {0.301511345f, -0.904534034f, 0.301511345f},
      {0.267261242f, 0.534522484f, 0.801783726f},
      {-0.408248290f, 0.816496581f, -0.408248290f},
      {0.923879533f, 0.382683432f, 0.0f},
      {-0.270598050f, 0.653281482f, 0.707106781f},
      {0.5f, -0.5f, 0.707106781f},
  };
  static constexpr math::Vec4 ORIGINS[] = {
      {{0.17f, 0.31f, 0.43f, 0.59f}},
      {{0.91f, 0.07f, 0.73f, 0.37f}},
  };
  static constexpr std::array<float, 6> ROTATIONS[] = {
      {0.0f, 0.3f, 0.7f, 0.0f, 0.0f, 0.0f},
      {1.1f, 2.3f, 0.4f, 0.0f, 0.0f, 0.0f},
  };
  static constexpr ShadeSample GOLDEN[] = {
      {11755, 29236, 38374, 21052},
      {6289, 20810, 33584, 7906},
      {7018, 22013, 34137, 89},
      {8005, 23674, 35054, 44437},
      {443, 411, 3629, 73},
      {16091, 34459, 41337, 23209},
      {25323, 42646, 46187, 52135},
      {13974, 29680, 38765, 13992},
      {578, 1541, 12313, 3619},
      {2037, 10615, 27205, 18},
      {522, 569, 5266, 27},
      {24216, 41698, 45666, 720},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {6971, 21237, 32816, 3214},
      {534, 605, 5610, 593},
      {7245, 22485, 34485, 361},
      {2024, 10399, 26729, 40},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
  };

  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  ShadeSample rendered[std::size(GOLDEN)];
  size_t row = 0;
  for (size_t preset = 0; preset < std::size(ORIGINS); ++preset) {
    const HL::FrameState frame{
        HyperLatticeWhiteBox::Effect::preset_params(preset),
        ORIGINS[preset],
        ROTATIONS[preset],
        HL::pixel_half_angle<288, 144>(),
        HyperLatticeWhiteBox::depth_palette(effect),
        HyperLatticeWhiteBox::axis_palette(effect),
    };
    const HL::PreparedTrace prepared = HL::prepare_trace(frame);
    for (size_t sample = 0; sample < std::size(DIRECTIONS); ++sample) {
      const Color4 color =
          HL::shade({DIRECTIONS[sample], 0.0f}, frame, prepared);
      rendered[row++] = {color.color.r, color.color.g, color.color.b,
                         frac_to_q16(color.alpha)};
    }
  }
  expect_shade_samples("render_signature", rendered, GOLDEN, std::size(GOLDEN),
                       std::size(DIRECTIONS), 9115850422160364994ull);
}

inline void test_specialized_slice_transition() {
  reset_globals();
  static constexpr float AMOUNTS[] = {0.5f, 0.625f, 0.75f, 0.875f};

  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  const HL::Params start = HyperLatticeWhiteBox::Effect::preset_params(0);
  const HL::Params target = HyperLatticeWhiteBox::Effect::preset_params(1);
  const auto compare_shells = [&]<uint8_t SHELL_COUNT>(HL::ShellCount shells) {
    float max_visible_error = 0.0f;
    float max_alpha_error = 0.0f;
    for (float amount : AMOUNTS) {
      HL::FrameState frame{};
      frame.params.lerp(start, target, amount);
      frame.params.shells = shells;
      frame.origin = {{0.17f, 0.31f, 0.43f, 0.59f}};
      frame.rotation_phase = {0.2f, 1.7f, 2.8f, 0.9f, 1.3f, 2.1f};
      frame.pixel_half_angle = HL::pixel_half_angle<288, 144>();
      frame.depth_palette = HyperLatticeWhiteBox::depth_palette(effect);
      frame.axis_palette = HyperLatticeWhiteBox::axis_palette(effect);
      const HL::PreparedTrace prepared = HL::prepare_trace(frame);
      for (int y = 0; y < 144; ++y)
        for (int x = 0; x < 288; ++x) {
          const math::Vector direction = math::pixel_to_vector<288, 144>(x, y);
          const Color4 exact = HL::shade({direction, 0.0f}, frame, prepared);
          const Color4 specialized = HL::shade_mode<true, SHELL_COUNT>(
              {direction, 0.0f}, frame, prepared);
          max_visible_error =
              std::max(max_visible_error,
                       fabsf(static_cast<float>(specialized.color.r) *
                                 specialized.alpha -
                             static_cast<float>(exact.color.r) * exact.alpha));
          max_visible_error =
              std::max(max_visible_error,
                       fabsf(static_cast<float>(specialized.color.g) *
                                 specialized.alpha -
                             static_cast<float>(exact.color.g) * exact.alpha));
          max_visible_error =
              std::max(max_visible_error,
                       fabsf(static_cast<float>(specialized.color.b) *
                                 specialized.alpha -
                             static_cast<float>(exact.color.b) * exact.alpha));
          max_alpha_error =
              std::max(max_alpha_error, fabsf(specialized.alpha - exact.alpha));
        }
    }
    HS_EXPECT_NEAR(max_visible_error, 0.0f, 1.0f);
    HS_EXPECT_NEAR(max_alpha_error, 0.0f, 5.0e-6f);
  };
  compare_shells.operator()<2>(HL::ShellCount::TWO);
  compare_shells.operator()<3>(HL::ShellCount::THREE);
}

/**
 * @brief Pins the specialized 4D-slice pipeline over the same style of sample.
 * @details Same table and pre-check as test_render_signature(), over
 * SpecializedRenderPipeline<2>'s prepare/evaluate pair at preset 1.
 * Provenance: no generator emits either. Re-derive by printing the RGB and Q16
 * alpha of every sample from an IEEE build of this case and pasting the table
 * and its fold back.
 */
inline void test_specialized_render_signature() {
  reset_globals();
  static constexpr math::Vector DIRECTIONS[] = {
      {1.0f, 0.0f, 0.0f},
      {0.0f, 1.0f, 0.0f},
      {0.0f, 0.0f, 1.0f},
      {0.577350269f, 0.577350269f, 0.577350269f},
      {-0.707106781f, 0.707106781f, 0.0f},
      {0.301511345f, -0.904534034f, 0.301511345f},
      {0.267261242f, 0.534522484f, 0.801783726f},
      {-0.408248290f, 0.816496581f, -0.408248290f},
      {0.923879533f, 0.382683432f, 0.0f},
      {-0.270598050f, 0.653281482f, 0.707106781f},
      {0.5f, -0.5f, 0.707106781f},
      {-1.0f, 0.0f, 0.0f},
  };
  static constexpr math::Vec4 ORIGINS[] = {
      {{0.17f, 0.31f, 0.43f, 0.59f}},
      {{0.91f, 0.07f, 0.73f, 0.37f}},
  };
  static constexpr std::array<float, 6> ROTATIONS[] = {
      {0.0f, 0.3f, 0.7f, 0.11f, 0.23f, 0.41f},
      {1.1f, 2.3f, 0.4f, 0.53f, 0.19f, 0.87f},
  };
  static constexpr ShadeSample GOLDEN[] = {
      {9816, 26581, 36812, 9757},
      {0, 0, 0, 0},
      {3644, 15194, 30453, 1837},
      {0, 0, 0, 0},
      {15610, 33375, 40320, 5667},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {1301, 5768, 18384, 48},
      {595, 969, 8822, 592},
      {1996, 10112, 26026, 182},
      {726, 2921, 14849, 107},
      {6757, 21039, 32748, 4104},
      {608, 1338, 10913, 355},
      {1417, 6578, 20172, 4073},
      {2069, 10950, 28047, 1675},
      {3522, 15066, 30626, 9197},
      {4923, 18222, 32574, 0},
      {0, 0, 0, 0},
      {2011, 10321, 26590, 116},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {758, 3076, 15069, 992},
      {3106, 14001, 30108, 830},
  };

  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  ShadeSample rendered[std::size(GOLDEN)];
  size_t row = 0;
  for (size_t index = 0; index < std::size(ORIGINS); ++index) {
    const HL::FrameState context{
        HyperLatticeWhiteBox::Effect::preset_params(1),
        ORIGINS[index],
        ROTATIONS[index],
        HL::pixel_half_angle<288, 144>(),
        HyperLatticeWhiteBox::depth_palette(effect),
        HyperLatticeWhiteBox::axis_palette(effect),
    };
    const auto frame = HL::SpecializedRenderPipeline<2>::prepare(context);
    for (size_t sample = 0; sample < std::size(DIRECTIONS); ++sample) {
      const Color4 color = HL::SpecializedRenderPipeline<2>::evaluate(
          DIRECTIONS[sample], frame.ctx, frame.prepared);
      rendered[row++] = {color.color.r, color.color.g, color.color.b,
                         frac_to_q16(color.alpha)};
    }
  }
#if defined(HS_TEST_FAST_MATH)
  // fast_wire_coverage's Newton reciprocal reassociates under the shipping
  // flag pair; the IEEE legs carry the table.
  (void)rendered;
  hs_test::skip_case(__func__,
                     "HS_TEST_FAST_MATH: fast_wire_coverage reassociates");
#else
  expect_shade_samples("specialized_render_signature", rendered, GOLDEN,
                       std::size(GOLDEN), std::size(DIRECTIONS),
                       17325614847026740988ull);
#endif
}

inline void test_presets_and_pipeline() {
  using Effect = HyperLattice<96, 20>;
  static_assert(HL::RenderPipeline::STAGE_COUNT == 1);
  static_assert(HL::RenderPipeline::Validation::MONOTONE);
  static_assert(HL::RenderPipeline::Validation::ENTRY);
  static_assert(HL::RenderPipeline::Validation::EXIT);
  for (size_t index = 0; index < Effect::PRESET_IDS.size(); ++index)
    HS_EXPECT_TRUE(Effect::valid_params(Effect::preset_params(index)));
  static_assert(Effect::PRESET_IDS.size() == 3);
  static_assert(Effect::PRESET_IDS[0] == "cubic-flight");
  static_assert(Effect::PRESET_IDS[1] == "hypercube-flight");
  static_assert(Effect::PRESET_IDS[2] == "deep-grid");

  constexpr HL::Params CUBIC_PRESET = Effect::preset_params(0);
  static_assert(CUBIC_PRESET.mode == HL::LatticeMode::THREE_D);
  static_assert(CUBIC_PRESET.color == HL::ColorMode::DEPTH);
  static_assert(CUBIC_PRESET.shells == HL::ShellCount::TWO);

  constexpr HL::Params SLICE_PRESET = Effect::preset_params(1);
  static_assert(SLICE_PRESET.mode == HL::LatticeMode::FOUR_D_SLICE);
  static_assert(SLICE_PRESET.color == HL::ColorMode::DEPTH);
  static_assert(SLICE_PRESET.shells == HL::ShellCount::TWO);

  constexpr HL::Params DEEP_PRESET = Effect::preset_params(2);
  static_assert(DEEP_PRESET.mode == HL::LatticeMode::THREE_D);
  static_assert(DEEP_PRESET.color == HL::ColorMode::DEPTH);
  static_assert(DEEP_PRESET.shells == HL::ShellCount::TWO);
  reset_globals();
  Effect effect;
  effect.init();
  std::vector<std::vector<float>> values;
  for (size_t index = 0; index < Effect::PRESET_IDS.size(); ++index) {
    HS_EXPECT_TRUE(effect.selectPreset(index));
    std::vector<float> current;
    for (const auto &def : effect.getParameters()) {
      HS_EXPECT_GE(def.get(), def.min);
      HS_EXPECT_LE(def.get(), def.max);
      current.push_back(def.get());
    }
    for (const auto &previous : values)
      HS_EXPECT_TRUE(current != previous);
    values.push_back(std::move(current));
  }
}

inline void test_dimension_dropdown_and_mode_lerp() {
  reset_globals();
  using Effect = HyperLatticeWhiteBox::Effect;
  Effect effect;
  effect.init();
  const ParamDef *cell_size = effect.getParameters().find("Cell Size");
  HS_EXPECT_TRUE(cell_size != nullptr);
  HS_EXPECT_EQ(cell_size->max, 10.0f);
  const ParamDef *far_distance = effect.getParameters().find("Far Distance");
  HS_EXPECT_TRUE(far_distance != nullptr);
  HS_EXPECT_TRUE(effect.getParameters().find("Far Cells") == nullptr);
  const ParamDef *dimension = effect.getParameters().find("Dimension");
  HS_EXPECT_TRUE(dimension != nullptr);
  HS_EXPECT_TRUE(dimension->is_enum());
  HS_EXPECT_EQ(dimension->option_count, 3);
  HS_EXPECT_EQ(std::string_view(dimension->options[0]), std::string_view("3D"));
  HS_EXPECT_EQ(std::string_view(dimension->options[1]),
               std::string_view("Dimensional Rift"));
  HS_EXPECT_EQ(std::string_view(dimension->options[2]),
               std::string_view("4D Slice"));
  HS_EXPECT_EQ(std::string_view(dimension->export_options[2]),
               std::string_view("LatticeMode::FOUR_D_SLICE"));
  HS_EXPECT_EQ(effect.updateParameter("Dimension", 2.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(HyperLatticeWhiteBox::params(effect).mode,
               HL::LatticeMode::FOUR_D_SLICE);

  HL::Params start;
  HL::Params target;
  target.mode = HL::LatticeMode::FOUR_D_SLICE;
  HL::Params blended;
  blended.lerp(start, target, 0.49f);
  HS_EXPECT_EQ(blended.mode, HL::LatticeMode::THREE_D);
  blended.lerp(start, target, 0.5f);
  HS_EXPECT_EQ(blended.mode, HL::LatticeMode::FOUR_D_SLICE);
}

/**
 * @brief DIMENSIONAL_RIFT traces the blended edge metric: the pure modes'
 * planes, at coverage between theirs.
 * @details The hyper metric adds the w component the cubic metric ignores, so
 * a small w puts the rift's blend strictly between the two on the same plane.
 */
inline void test_dimensional_rift_layers() {
  HL::FrameState frame{};
  frame.params = HyperLattice<96, 20>::preset_params(0);
  frame.params.sphere_radius = 0.0f;
  frame.origin = {{0.25f, 0.02f, 0.03f, 0.05f}};
  const auto nearest = [&](HL::LatticeMode mode) {
    frame.params.mode = mode;
    return HL::trace(math::X_AXIS, HL::prepare_trace(frame));
  };
  const HL::TraceHit cubic = nearest(HL::LatticeMode::THREE_D);
  const HL::TraceHit hyper = nearest(HL::LatticeMode::FOUR_D_SLICE);
  const HL::TraceHit rift = nearest(HL::LatticeMode::DIMENSIONAL_RIFT);
  HS_EXPECT_EQ(HL::prepare_trace(frame).dimension_mix,
               HL::DIMENSIONAL_RIFT_MIX);
  HS_EXPECT_GT(hyper.coverage, 0.0f);
  HS_EXPECT_NEAR(rift.distance, cubic.distance, 1e-6f);
  HS_EXPECT_NEAR(rift.distance, hyper.distance, 1e-6f);
  HS_EXPECT_GT(cubic.coverage, rift.coverage);
  HS_EXPECT_GT(rift.coverage, hyper.coverage);
  HS_EXPECT_EQ(cubic.free_axis, uint8_t(2));
  HS_EXPECT_EQ(hyper.free_axis, uint8_t(3));
  HS_EXPECT_EQ(rift.free_axis, hyper.free_axis);

  frame.params.sphere_radius = 1.0f;
  frame.origin = {{0.17f, 0.31f, 0.43f, 0.59f}};
  frame.rotation_phase = {0.2f, 1.7f, 2.8f, 0.9f, 1.3f, 2.1f};
  const HL::PreparedTrace prepared = HL::prepare_trace(frame);
  static constexpr math::Vector DIRECTIONS[] = {
      {1.0f, 0.0f, 0.0f},
      {0.0f, 1.0f, 0.0f},
      {0.0f, 0.0f, 1.0f},
      {0.577350269f, 0.577350269f, 0.577350269f},
      {-0.707106781f, 0.707106781f, 0.0f},
      {0.301511345f, -0.904534034f, 0.301511345f},
  };
  int layers = 0;
  for (const math::Vector &direction : DIRECTIONS) {
    float previous_distance = 0.0f;
    HL::trace_layers(direction, prepared, [&](const HL::TraceHit &hit) {
      HS_EXPECT_GT(hit.coverage, 0.0f);
      HS_EXPECT_LE(hit.coverage, 1.0f);
      HS_EXPECT_GT(hit.distance, previous_distance);
      HS_EXPECT_LT(hit.free_axis, uint8_t(HL::DIMENSIONS));
      previous_distance = hit.distance;
      ++layers;
      return true;
    });
  }
  HS_EXPECT_GT(layers, 0);
}

/**
 * @brief ColorMode::AXIS shades through the axis palette and ShellCount::ONE
 * fades its only shell at the horizon; neither takes the specialized slice.
 */
inline void test_axis_color_and_single_shell() {
  reset_globals();
  using Effect = HyperLatticeWhiteBox::Effect;
  Effect effect;
  effect.init();
  HL::FrameState frame{};
  frame.params = Effect::preset_params(1);
  frame.params.sphere_radius = 0.0f;
  frame.origin = {{0.25f, 0.02f, 0.01f, 0.43f}};
  frame.pixel_half_angle = HL::pixel_half_angle<96, 20>();
  frame.depth_palette = HyperLatticeWhiteBox::depth_palette(effect);
  frame.axis_palette = HyperLatticeWhiteBox::axis_palette(effect);
  HS_EXPECT_TRUE(Effect::uses_specialized_slice(frame.params));
  frame.params.color = HL::ColorMode::AXIS;
  HS_EXPECT_FALSE(Effect::uses_specialized_slice(frame.params));
  frame.params.color = HL::ColorMode::DEPTH;
  frame.params.shells = HL::ShellCount::ONE;
  HS_EXPECT_FALSE(Effect::uses_specialized_slice(frame.params));

  const auto layers = [&](HL::ShellCount shells) {
    frame.params.shells = shells;
    std::vector<HL::TraceHit> hits;
    HL::trace_layers(math::X_AXIS, HL::prepare_trace(frame),
                     [&](const HL::TraceHit &hit) {
                       hits.push_back(hit);
                       return true;
                     });
    return hits;
  };
  const std::vector<HL::TraceHit> two = layers(HL::ShellCount::TWO);
  const std::vector<HL::TraceHit> one = layers(HL::ShellCount::ONE);
  HS_EXPECT_EQ(two.size(), size_t{2});
  HS_EXPECT_EQ(one.size(), size_t{1});
  HS_EXPECT_EQ(one.front().distance, two.front().distance);
  HS_EXPECT_GT(one.front().coverage, 0.0f);
  HS_EXPECT_LT(one.front().coverage, two.front().coverage);
  HS_EXPECT_NEAR(one.front().coverage,
                 two.front().coverage *
                     HL::shell_horizon_coverage(0, 1, one.front().distance,
                                                1.0f / frame.params.cell_size),
                 1e-6f);

  frame.params.color = HL::ColorMode::AXIS;
  const HL::PreparedTrace prepared = HL::prepare_trace(frame);
  const HL::TraceHit hit = HL::trace(math::X_AXIS, prepared);
  const Color4 axis = HL::shade({math::X_AXIS, 0.0f}, frame, prepared);
  frame.params.color = HL::ColorMode::DEPTH;
  const Color4 depth =
      HL::shade({math::X_AXIS, 0.0f}, frame, HL::prepare_trace(frame));
  HS_EXPECT_NEAR(hit.coverage, one.front().coverage, 1e-6f);
  HS_EXPECT_NEAR(axis.alpha, hit.coverage, 1e-6f);
  HS_EXPECT_EQ(axis.alpha, depth.alpha);
  HS_EXPECT_TRUE(axis.color.r != depth.color.r ||
                 axis.color.g != depth.color.g ||
                 axis.color.b != depth.color.b);
  const float depth_fraction = hit.distance * prepared.inv_far;
  const Pixel expected =
      frame.axis_palette->get_color_unit(
          (static_cast<float>(hit.free_axis) + 0.75f * depth_fraction) / 4.0f) *
      (0.45f + 0.55f * (1.0f - depth_fraction));
  HS_EXPECT_NEAR(axis.color.r, expected.r, 1);
  HS_EXPECT_NEAR(axis.color.g, expected.g, 1);
  HS_EXPECT_NEAR(axis.color.b, expected.b, 1);
}

inline int run_hyper_lattice_tests() {
  hs_test::ModuleFixture fixture("hyper_lattice");
  test_periodic_distance();
  test_edge_metrics();
  test_so4_rotation();
  test_dimensional_rotation_wrap_is_continuous();
  test_resolution_aware_wire_coverage();
  test_near_field_fade();
  test_far_shell_fade();
  test_pause_does_not_stop_motion();
  test_depth_palette_mutates_slowly_while_paused();
  test_depth_palette_keeps_cool_character();
  test_next_plane_is_strict();
  test_trace_layers_are_front_to_back();
  test_layer_composite_reveals_background();
  test_surface_origin_parallax();
  test_hyperplane_event();
  test_coincident_planes_form_one_layer();
  test_render_signature();
  test_specialized_slice_transition();
  test_specialized_render_signature();
  test_presets_and_pipeline();
  test_dimension_dropdown_and_mode_lerp();
  test_dimensional_rift_layers();
  test_axis_color_and_single_shell();
  return fixture.result();
}

} // namespace hyper_lattice_tests
} // namespace hs_test
