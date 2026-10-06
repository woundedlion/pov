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

template <typename Effect>
inline const ParamDef *required_param(Effect &effect, const char *name) {
  const auto *parameter = effect.getParameters().find(name);
  HS_EXPECT_TRUE(parameter != nullptr);
  return parameter;
}

namespace Trace = HyperLatticeDetail::Trace;

/** @brief Prepares a traced frame over the tests' crossing scratch. */
inline Trace::Prepared prepare_traced(Trace::Settings settings) {
  static SDF::OctetTrace::CrossingStorage crossings;
  settings.crossings = &crossings;
  return Trace::prepare(settings);
}

namespace HL {
using Params = HyperLatticeDetail::Params;
using FrameState = HyperLatticeDetail::FrameState;
using LatticeMode = HyperLatticeDetail::LatticeMode;
using ShellCount = HyperLatticeDetail::ShellCount;
using RenderPipeline = HyperLatticeDetail::RenderPipeline;
template <uint8_t N>
using SpecializedRenderPipeline =
    HyperLatticeDetail::SpecializedRenderPipeline<N>;
using HyperLatticeDetail::pixel_half_angle;
using namespace SDF::Lattice;
struct PreparedTrace : SDF::Lattice::PreparedTrace {
  Raycast::Appearance appearance;
};
/** @brief Crossing scratch the helpers bind to frames that carry none. */
inline HyperLatticeDetail::CrossingList crossings;
inline PreparedTrace prepare_trace(const FrameState &frame) {
  FrameState bound = frame;
  if (!bound.crossings)
    bound.crossings = &crossings;
  const auto p = HyperLatticeDetail::prepare_trace(bound);
  return {p.lattice, p.appearance};
}
template <bool SLICE = false>
auto trace_plane(const math::Vec4 &origin, const math::Vec4 &direction,
                 int axis, float distance, float step,
                 const PreparedTrace &prepared) {
  auto result = SDF::Lattice::trace_plane<SLICE>(origin, direction, axis,
                                                 distance, step, prepared);
  result.coverage *= prepared.appearance.opacity(distance);
  return result;
}
template <typename Consume>
void trace_layers(const math::Vector &normal, const PreparedTrace &prepared,
                  Consume consume) {
  SDF::Lattice::Events<> events(normal, prepared);
  Raycast::trace_events(
      events, {0, prepared.far_distance}, {},
      [&](const Raycast::Contribution &hit) {
        const float coverage =
            hit.coverage * prepared.appearance.opacity(hit.t);
        return coverage <= 0 ||
               consume(TraceHit{coverage, hit.t,
                                static_cast<uint8_t>(hit.feature)});
      });
}
inline TraceHit trace(const math::Vector &normal,
                      const PreparedTrace &prepared) {
  TraceHit result;
  trace_layers(normal, prepared, [&](const TraceHit &hit) {
    result = hit;
    return false;
  });
  return result;
}
template <bool SLICE = false, uint8_t SHELLS = 0>
Color4 shade_mode(const Pullback::SphereSample &input, const FrameState &frame,
                  const PreparedTrace &prepared) {
  return HyperLatticeDetail::Renderer<SLICE, SHELLS>::shade(
      input.dir, frame, {prepared, prepared.appearance, &crossings});
}
inline Color4 shade(const Pullback::SphereSample &input,
                    const FrameState &frame, const PreparedTrace &prepared) {
  return shade_mode(input, frame, prepared);
}
} // namespace HL

struct HyperLatticeWhiteBox {
  using Effect = HyperLattice<96, 20>;

  static math::Vec4 origin(const Effect &effect) { return effect.origin; }
  static std::array<float, 6> rotation_phase(const Effect &effect) {
    return effect.rotation_phase;
  }
  static HL::Params &params(Effect &effect) { return effect.params; }
  /** @brief Steps the timeline and the preset choreography one frame. */
  static void step_choreography(Effect &effect, Canvas &canvas) {
    effect.timeline.step(canvas);
    effect.step_choreography();
  }
  static float preset_gain(const Effect &effect) { return effect.preset_gain; }
  static math::Vec4 trace_center(const Effect &effect) {
    return effect.trace_center;
  }
  static void advance_state(Effect &effect) { effect.advance_state(); }
  static void center_camera(Effect &effect) {
    effect.origin = {};
    effect.trace_center = {};
    effect.rotation_phase = {};
  }
  static HL::FrameState frame(const Effect &effect) {
    return {effect.params,
            effect.origin,
            effect.rotation_phase,
            HyperLatticeDetail::pixel_half_angle<96, 20>(),
            &effect.depth_palette.palette(),
            effect.preset_gain,
            effect.crossing_list};
  }
  static void step_depth_palette(Effect &effect) {
    effect.depth_palette.step();
  }
  static bool depth_palette_fading(const Effect &effect) {
    return effect.depth_palette.fading();
  }
  static Pixel depth_color(const Effect &effect, float amount) {
    return effect.depth_palette.palette().get(amount).color;
  }
  static const BakedPalette *depth_palette(const Effect &effect) {
    return &effect.depth_palette.palette();
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
  frame.params.mode = HL::LatticeMode::FOUR_D_SLICE;
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
  const Raycast::Appearance FADE{
      .inv_far = 0, .near_start = NEAR_START, .near_inv_span = NEAR_INV_SPAN};
  HS_EXPECT_EQ(FADE.opacity(0.15f), 0.0f);
  HS_EXPECT_GT(FADE.opacity(0.275f), 0.0f);
  HS_EXPECT_LT(FADE.opacity(0.275f), 1.0f);
  HS_EXPECT_EQ(FADE.opacity(0.4f), 1.0f);
  for (HL::LatticeMode mode :
       {HL::LatticeMode::THREE_D, HL::LatticeMode::FOUR_D_SLICE}) {
    for (float radius : {0.0f, 1.0f, 2.0f}) {
      for (float cell_size : {0.25f, 1.0f, 10.0f}) {
        HL::FrameState frame{};
        frame.params.mode = mode;
        frame.params.sphere_radius = radius;
        frame.params.cell_size = cell_size;
        frame.params.far_distance = 16.0f;
        frame.params.near_fade = 0.2f;
        const auto prepared = HL::prepare_trace(frame);
        const float span = frame.params.near_fade * cell_size * (1.0f + radius);
        for (float fraction :
             {0.0f, 0.001f, 0.25f, 0.5f, 0.75f, 0.999f, 1.0f}) {
          const float distance =
              prepared.appearance.near_start + fraction * span;
          const math::Vec4 origin{{-distance / cell_size, 0.0f, 0.0f, 0.0f}};
          const math::Vec4 direction{{1.0f / cell_size, 0.0f, 0.0f, 0.0f}};
          const auto hit = HL::trace_plane<false>(
              origin, direction, 0, distance, cell_size, prepared);
          const float fog = 1.0f - distance / frame.params.far_distance;
          HS_EXPECT_NEAR(hit.coverage / (fog * fog),
                         fraction * fraction * (3.0f - 2.0f * fraction), 1e-5f);
          if (mode == HL::LatticeMode::FOUR_D_SLICE) {
            const auto specialized = HL::trace_plane<true>(
                origin, direction, 0, distance, cell_size, prepared);
            HS_EXPECT_NEAR(specialized.coverage, hit.coverage, 1e-5f);
          }
        }
      }
    }
  }
  HL::FrameState frame{};
  frame.params.sphere_radius = 0.0f;
  const auto centered = HL::prepare_trace(frame);
  frame.params.sphere_radius = 1.0f;
  const auto surface = HL::prepare_trace(frame);
  HS_EXPECT_NEAR(surface.appearance.near_start,
                 centered.appearance.near_start * 2.0f, 1e-6f);
  HS_EXPECT_NEAR(surface.appearance.near_inv_span,
                 centered.appearance.near_inv_span * 0.5f, 1e-6f);
  const float distance = centered.appearance.near_start +
                         0.25f / centered.appearance.near_inv_span;
  HS_EXPECT_LT(surface.appearance.opacity(distance),
               centered.appearance.opacity(distance));
  HS_EXPECT_EQ(surface.appearance.opacity(0.0f), 0.0f);
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
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_TRUE(HyperLatticeWhiteBox::depth_palette_fading(effect));
  HS_EXPECT_EQ(HyperLatticeWhiteBox::depth_color(effect, 0.5f), initial);

  for (int frame = 0; frame < 480; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_NE(HyperLatticeWhiteBox::depth_color(effect, 0.5f), initial);
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
  frame.params = HyperLattice<96, 20>::preset(0).params;
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
  frame.params = HyperLattice<96, 20>::preset(0).params;
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
  frame.params =
      HyperLattice<96, 20>::preset(HyperLattice<96, 20>::HYPERCUBE_PRESET_INDEX)
          .params;
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
  frame.params = HyperLattice<96, 20>::preset(0).params;
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
 * whose last bits differ between libm builds.
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
 * @param label Context label for out-of-band channels.
 * @param rendered Samples in fold order, one row per (preset, direction).
 * @param golden Rows the pin was recorded from.
 * @param count Row count of both tables.
 * @param per_preset Rows per preset, so a failing row names preset and sample.
 * @param pin Signature the golden table folds to.
 * @details A run whose signature moved falls through to the per-channel band,
 * which names the preset, sample and channel that drifted.
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
    HS_EXPECT_NEAR(got.alpha, want.alpha, 1);
  }
}

/**
 * @brief Pins shade() over a fixed direction and preset sample.
 * @details GOLDEN is the oracle; its FNV-1a 64 signature is the pre-check.
 * Re-derive by printing the sample table's RGB and Q16 alpha from this case
 * built by the native clang test toolchain (cmake/toolchain-native-clang.cmake)
 * and pasting the table and its fold back. The coverage ramp here is an IEEE
 * division, so the table reproduces under -ffast-math -fno-finite-math-only.
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
      {11664, 29080, 38288, 17981},
      {6289, 20810, 33585, 7906},
      {7018, 22012, 34137, 89},
      {8004, 23671, 35053, 44366},
      {443, 411, 3630, 73},
      {15859, 34189, 41182, 14730},
      {23945, 40773, 44710, 11477},
      {9649, 23679, 35425, 9199},
      {578, 1541, 12313, 3619},
      {2037, 10614, 27205, 18},
      {522, 569, 5266, 27},
      {22397, 39389, 44357, 142},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {2986, 13626, 29849, 4010},
      {335, 272, 2239, 6},
      {2856, 13378, 29840, 1664},
      {1290, 5736, 18677, 7},
      {958, 3979, 15939, 110},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {1621, 7785, 22182, 95},
  };

  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  ShadeSample rendered[std::size(GOLDEN)];
  size_t row = 0;
  static constexpr size_t PRESETS[] = {
      HyperLatticeWhiteBox::Effect::CUBIC_PRESET_INDEX,
      HyperLatticeWhiteBox::Effect::HYPERCUBE_PRESET_INDEX};
  for (size_t preset = 0; preset < std::size(ORIGINS); ++preset) {
    const HL::FrameState frame{
        HyperLatticeWhiteBox::Effect::preset(PRESETS[preset]).params,
        ORIGINS[preset],
        ROTATIONS[preset],
        HL::pixel_half_angle<288, 144>(),
        HyperLatticeWhiteBox::depth_palette(effect),
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
                       std::size(DIRECTIONS), 4451398320181246065ull);
}

inline void test_specialized_slice_transition() {
  reset_globals();
  static constexpr float AMOUNTS[] = {0.5f, 0.625f, 0.75f, 0.875f};

  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  const HL::Params start = HyperLatticeWhiteBox::Effect::preset(0).params;
  const HL::Params target =
      HyperLatticeWhiteBox::Effect::preset(
          HyperLatticeWhiteBox::Effect::HYPERCUBE_PRESET_INDEX)
          .params;
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
      const HL::PreparedTrace prepared = HL::prepare_trace(frame);
      int lit = 0;
      for (int y = 0; y < 144; ++y)
        for (int x = 0; x < 288; ++x) {
          const math::Vector direction = math::pixel_to_vector<288, 144>(x, y);
          const Color4 exact = HL::shade({direction, 0.0f}, frame, prepared);
          lit += exact.alpha > 0.0f ? 1 : 0;
          const Color4 specialized = HL::shade_mode<true, SHELL_COUNT>(
              {direction, 0.0f}, frame, prepared);
          max_visible_error = hs_test::fold_worst(
              max_visible_error,
              fabsf(static_cast<float>(specialized.color.r) *
                        specialized.alpha -
                    static_cast<float>(exact.color.r) * exact.alpha));
          max_visible_error = hs_test::fold_worst(
              max_visible_error,
              fabsf(static_cast<float>(specialized.color.g) *
                        specialized.alpha -
                    static_cast<float>(exact.color.g) * exact.alpha));
          max_visible_error = hs_test::fold_worst(
              max_visible_error,
              fabsf(static_cast<float>(specialized.color.b) *
                        specialized.alpha -
                    static_cast<float>(exact.color.b) * exact.alpha));
          max_alpha_error = hs_test::fold_worst(
              max_alpha_error, fabsf(specialized.alpha - exact.alpha));
        }
      HS_EXPECT_GT(lit, 288 * 144 / 2);
    }
    HS_EXPECT_NEAR(max_visible_error, 0.0f, 1.0f);
    HS_EXPECT_NEAR(max_alpha_error, 0.0f, 5.0e-6f);
  };
  compare_shells.operator()<2>(HL::ShellCount::TWO);
  compare_shells.operator()<0>(HL::ShellCount::THREE);
}

/**
 * @brief Pins the specialized 4D-slice pipeline over the same style of sample.
 * @details Same table layout and signature pre-check as test_render_signature(), over
 * SpecializedRenderPipeline<2>'s prepare/evaluate pair at HYPERCUBE_PRESET_INDEX.
 * Re-derive by printing the RGB and Q16
 * alpha of the sample table from an IEEE build of this case and pasting the table
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
      {3791, 15608, 30853, 10910},
      {22451, 40404, 44803, 4202},
      {2061, 10983, 28120, 2087},
      {0, 0, 0, 0},
      {6428, 20914, 33453, 7343},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {591, 1572, 12363, 73},
      {382, 326, 2757, 67},
      {1528, 7343, 21515, 646},
      {544, 634, 5876, 10},
      {2776, 13147, 29664, 4776},
      {426, 386, 3330, 5},
      {563, 1769, 13611, 1740},
      {1660, 7986, 22501, 1517},
      {2053, 10938, 28005, 7814},
      {22640, 40562, 44900, 57472},
      {0, 0, 0, 0},
      {1277, 5601, 18399, 68},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {2383, 10967, 25672, 1072},
      {2037, 10670, 27344, 1402},
  };

  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  ShadeSample rendered[std::size(GOLDEN)];
  size_t row = 0;
  for (size_t index = 0; index < std::size(ORIGINS); ++index) {
    const HL::FrameState context{
        HyperLatticeWhiteBox::Effect::preset(
            HyperLatticeWhiteBox::Effect::HYPERCUBE_PRESET_INDEX)
            .params,
        ORIGINS[index],
        ROTATIONS[index],
        HL::pixel_half_angle<288, 144>(),
        HyperLatticeWhiteBox::depth_palette(effect),
        1.0f,
        &HL::crossings,
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
                       17793107560919685780ull);
#endif
}

inline void test_presets_and_pipeline() {
  using Effect = HyperLattice<96, 20>;
  static_assert(HL::RenderPipeline::STAGE_COUNT == 1);
  static_assert(HL::RenderPipeline::Validation::MONOTONE);
  static_assert(HL::RenderPipeline::Validation::ENTRY);
  static_assert(HL::RenderPipeline::Validation::EXIT);
  for (size_t index = 0; index < Effect::PRESET_IDS.size(); ++index)
    HS_EXPECT_TRUE(Effect::valid_params(Effect::preset(index).params));
  static_assert(Effect::PRESET_IDS.size() == Effect::SHELL_4D_PRESET_INDEX + 1);
  static_assert(Effect::PRESET_IDS[Effect::CUBIC_PRESET_INDEX] ==
                "cubic-flight");
  static_assert(Effect::PRESET_IDS[Effect::WIDE_PRESET_INDEX] ==
                "cubic-wide-flight");
  static_assert(Effect::PRESET_IDS[Effect::HYPERCUBE_PRESET_INDEX] ==
                "hypercube-flight");

  constexpr HL::Params CUBIC_PRESET = Effect::preset(0).params;
  static_assert(CUBIC_PRESET.mode == HL::LatticeMode::THREE_D);
  static_assert(CUBIC_PRESET.shells == HL::ShellCount::TWO);

  constexpr HL::Params SLICE_PRESET =
      Effect::preset(Effect::HYPERCUBE_PRESET_INDEX).params;
  static_assert(SLICE_PRESET.mode == HL::LatticeMode::FOUR_D_SLICE);
  static_assert(SLICE_PRESET.shells == HL::ShellCount::TWO);

  static_assert(SLICE_PRESET.sphere_radius == 0.0f);
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

inline void test_configuration_adoption_and_snapshots() {
  reset_globals();
  using Effect = HyperLatticeWhiteBox::Effect;
  Effect effect;
  effect.init();
  const auto initial = effect.serialize_parameters();
  auto stale = initial;
  stale.schema_version = Effect::PARAMETER_SCHEMA_VERSION - 1;
  stale.params.cell_size = 5;
  HS_EXPECT_FALSE(effect.restore_parameters(stale));
  HS_EXPECT_EQ(effect.serialize_parameters().params.cell_size,
               initial.params.cell_size);
  auto invalid = initial;
  invalid.params.mode = static_cast<HL::LatticeMode>(2);
  HS_EXPECT_FALSE(effect.restore_parameters(invalid));
  HS_EXPECT_EQ(effect.serialize_parameters().params.mode, initial.params.mode);
  HS_EXPECT_EQ(effect.updateParameter("4D Spin", .01f),
               ParamSetResult::READONLY);
  const auto schema = effect.getParameterSchemaGeneration();
  HS_EXPECT_EQ(effect.updateParameter("Cell Size", 5), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("View", 1), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.serialize_parameters().params.cell_size,
               Effect::preset(Effect::HYPERCUBE_PRESET_INDEX).params.cell_size);
  HS_EXPECT_EQ(
      effect.serialize_parameters().params.wire_radius,
      Effect::preset(Effect::HYPERCUBE_PRESET_INDEX).params.wire_radius);
  HS_EXPECT_GT(effect.getParameterSchemaGeneration(), schema);
  HS_EXPECT_EQ(effect.updateParameter("4D Spin", .01f),
               ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(effect.restore_parameters(initial));
  const auto *spin_4d = required_param(effect, "4D Spin");
  if (!spin_4d)
    return;
  HS_EXPECT_TRUE(spin_4d->readonly);
  HL::Params start = Effect::preset(0).params;
  HL::Params target = Effect::preset(Effect::HYPERCUBE_PRESET_INDEX).params;
  target.cell_size = 4;
  HL::Params mixed;
  target.shells = HL::ShellCount::THREE;
  for (float t : {.0f, .1f, .49f, .5f, .75f, 1.f}) {
    mixed.lerp(start, target, t);
    const auto &expected = t < .5f ? start : target;
    HS_EXPECT_EQ(mixed.mode, expected.mode);
    HS_EXPECT_EQ(mixed.pattern, expected.pattern);
    HS_EXPECT_EQ(mixed.shells, expected.shells);
    HS_EXPECT_EQ(mixed.cell_size,
                 hs::lerp(start.cell_size, target.cell_size, t));
    HS_EXPECT_EQ(mixed.sphere_radius,
                 hs::lerp(start.sphere_radius, target.sphere_radius, t));
    HS_EXPECT_EQ(mixed.wire_radius,
                 hs::lerp(start.wire_radius, target.wire_radius, t));
    HS_EXPECT_EQ(mixed.speed, hs::lerp(start.speed, target.speed, t));
    HS_EXPECT_EQ(mixed.spin_3d, hs::lerp(start.spin_3d, target.spin_3d, t));
    HS_EXPECT_EQ(mixed.spin_4d, hs::lerp(start.spin_4d, target.spin_4d, t));
  }
}

/**
 * @brief Segues within one pattern and view morph at full brightness; any
 *        other segue holds each preset whole and dims through black.
 */
inline void test_family_segues() {
  using Effect = HyperLatticeWhiteBox::Effect;
  constexpr size_t COUNT = Effect::PRESET_IDS.size();
  for (size_t index = 0; index < COUNT; ++index) {
    const auto FROM = Effect::preset(index).params;
    const auto TO = Effect::preset((index + 1) % COUNT).params;
    const bool SAME_FAMILY = FROM.pattern == TO.pattern && FROM.mode == TO.mode;
    const auto DEPARTURE = Effect::preset(index).segue;
    HS_EXPECT_EQ(std::holds_alternative<Segue::Preset::Lerp>(DEPARTURE),
                 SAME_FAMILY);
    HS_EXPECT_EQ(std::holds_alternative<Segue::Preset::Fade>(DEPARTURE),
                 !SAME_FAMILY);
  }

  reset_globals();
  Effect effect;
  effect.init();
  const auto &params = HyperLatticeWhiteBox::params(effect);
  const auto WIDE = Effect::preset(Effect::WIDE_PRESET_INDEX).params;
  const auto HYPERCUBE = Effect::preset(Effect::HYPERCUBE_PRESET_INDEX).params;
  HS_EXPECT_TRUE(effect.selectPreset(Effect::WIDE_PRESET_INDEX));
  effect.setAnimationsPaused(false);
  Canvas canvas(effect);
  for (int frame = 0; frame < Effect::PRESET_DWELL_FRAMES; ++frame)
    HyperLatticeWhiteBox::step_choreography(effect, canvas);
  HS_EXPECT_EQ(HyperLatticeWhiteBox::preset_gain(effect), 1.0f);
  const int FADE =
      Segue::Preset::frames(Effect::preset(Effect::WIDE_PRESET_INDEX).segue);
  float darkest = 1.0f;
  for (int frame = 1; frame <= FADE; ++frame) {
    HyperLatticeWhiteBox::step_choreography(effect, canvas);
    const float GAIN = HyperLatticeWhiteBox::preset_gain(effect);
    darkest = fminf(darkest, GAIN);
    const auto &expected = frame < FADE / 2 ? WIDE : HYPERCUBE;
    HS_EXPECT_EQ(params.mode, expected.mode);
    HS_EXPECT_EQ(params.cell_size, expected.cell_size);
  }
  HS_EXPECT_EQ(darkest, 0.0f);
  HS_EXPECT_EQ(HyperLatticeWhiteBox::preset_gain(effect), 1.0f);
  HS_EXPECT_EQ(effect.getPresetIndex(), Effect::HYPERCUBE_PRESET_INDEX);
}

inline void test_dimension_dropdown_and_mode_lerp() {
  reset_globals();
  using Effect = HyperLatticeWhiteBox::Effect;
  Effect effect;
  effect.init();
  const auto *near_fade = required_param(effect, "Near Fade");
  if (!near_fade)
    return;
  HS_EXPECT_EQ(near_fade->get(), 0.5f);
  HS_EXPECT_EQ(effect.updateParameter("Near Fade", 1.2f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(HyperLatticeWhiteBox::params(effect).near_fade, 1.2f);
  HL::Params fade_start;
  HL::Params fade_target;
  fade_target.near_fade = 1.5f;
  HL::Params fade_blended;
  fade_blended.lerp(fade_start, fade_target, 0.5f);
  HS_EXPECT_EQ(fade_blended.near_fade, 1.0f);
  for (float invalid : {0.0f, 2.01f}) {
    fade_blended.near_fade = invalid;
    HS_EXPECT_FALSE(Effect::valid_params(fade_blended));
  }
  const auto *cell_size = required_param(effect, "Cell Size");
  if (!cell_size)
    return;
  HS_EXPECT_EQ(cell_size->max, 10.0f);
  const auto *far_distance = required_param(effect, "Far Distance");
  if (!far_distance)
    return;
  HS_EXPECT_TRUE(effect.getParameters().find("Far Cells") == nullptr);
  const auto *dimension = required_param(effect, "View");
  if (!dimension)
    return;
  HS_EXPECT_TRUE(dimension->is_enum());
  HS_EXPECT_EQ(dimension->option_count, 2);
  HS_EXPECT_EQ(std::string_view(dimension->options[0]),
               std::string_view("3D perspective"));
  HS_EXPECT_EQ(std::string_view(dimension->options[1]),
               std::string_view("4D slice"));
  HS_EXPECT_EQ(std::string_view(dimension->export_options[1]),
               std::string_view("LatticeMode::FOUR_D_SLICE"));
  HS_EXPECT_EQ(effect.updateParameter("View", 1.0f), ParamSetResult::APPLIED);
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
 * @brief ShellCount::ONE fades its only shell at the horizon and does not take
 * the specialized slice.
 */
inline void test_single_shell() {
  reset_globals();
  using Effect = HyperLatticeWhiteBox::Effect;
  Effect effect;
  effect.init();
  HL::FrameState frame{};
  frame.params = Effect::preset(Effect::HYPERCUBE_PRESET_INDEX).params;
  frame.params.sphere_radius = 0.0f;
  frame.origin = {{0.25f, 0.02f, 0.01f, 0.43f}};
  frame.pixel_half_angle = HL::pixel_half_angle<96, 20>();
  frame.depth_palette = HyperLatticeWhiteBox::depth_palette(effect);
  HS_EXPECT_TRUE(Effect::uses_specialized_slice(frame.params));
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
}

inline void test_octet_prepared_projection() {
  using Effect = HyperLatticeWhiteBox::Effect;
  using Domain = Raycast::SamplingDomain;
  reset_globals();
  Effect effect;
  effect.init();
  Trace::Settings settings;
  settings.palette = HyperLatticeWhiteBox::depth_palette(effect);
  settings.pixel_half_angle = .018f;
  const auto expect_same = [](const Raycast::ShadedTrace &actual,
                              const Raycast::ShadedTrace &expected) {
    HS_EXPECT_EQ(actual.trace.status, expected.trace.status);
    HS_EXPECT_NEAR(actual.color.alpha, expected.color.alpha, 3e-4f);
    HS_EXPECT_NEAR(actual.color.color.r, expected.color.color.r, 2);
    HS_EXPECT_NEAR(actual.color.color.g, expected.color.color.g, 2);
    HS_EXPECT_NEAR(actual.color.color.b, expected.color.color.b, 2);
  };
  const auto expect_premultiplied = [](const SDF::OctetTrace::Sample &actual,
                                       const Raycast::ShadedTrace &expected) {
    HS_EXPECT_EQ(actual.status, expected.trace.status);
    const Pixel PREMULTIPLIED = expected.color.color * expected.color.alpha;
    HS_EXPECT_NEAR(actual.color.r, PREMULTIPLIED.r, 1);
    HS_EXPECT_NEAR(actual.color.g, PREMULTIPLIED.g, 1);
    HS_EXPECT_NEAR(actual.color.b, PREMULTIPLIED.b, 1);
  };
  size_t lit = 0;
  for (Domain domain : {Domain::SPATIAL_3D, Domain::SLICE_4D}) {
    settings.domain = domain;
    settings.center = {
        {.231f, -.437f, .719f, domain == Domain::SLICE_4D ? .383f : 0.0f}};
    for (float angle : {0.0f, .31f, 1.23f}) {
      settings.embedding = math::Mat4::identity();
      math::rotate_plane(settings.embedding, 0, 1, angle);
      math::rotate_plane(settings.embedding, 1, 2, angle * .73f);
      if (domain == Domain::SLICE_4D) {
        math::rotate_plane(settings.embedding, 0, 3, angle * .41f);
        math::rotate_plane(settings.embedding, 2, 3, angle * 1.9f);
      }
      for (float cell : {.4f, 1.7f}) {
        settings.cell_size = cell;
        settings.wire_radius = .055f * cell;
        for (float radial : {0.0f, .7f, 2.3f}) {
          settings.radial_start = radial;
          auto prepared = prepare_traced(settings);
          HS_EXPECT_TRUE(prepared.valid);
          const auto &camera = prepared.camera;
          for (float near : {0.0f, .19f}) {
            prepared.camera.interval.near = near;
            for (int sample = 0; sample < 48; ++sample) {
              const float PHASE = sample * 2.39996323f;
              const float Z = (sample + .5f) / 24.0f - 1.0f;
              const float R = sqrtf(1.0f - Z * Z);
              const math::Vector DIRECTION(R * cosf(PHASE), R * sinf(PHASE), Z);
              const auto RAY = camera.ray(DIRECTION);
              const auto AMBIENT = camera.embedding.apply(
                  {{DIRECTION.x, DIRECTION.y, DIRECTION.z, 0.0f}});
              Raycast::ShadedTrace world;
              Raycast::ShadedTrace projected;
              SDF::OctetTrace::Sample traced;
              if (domain == Domain::SLICE_4D) {
                SDF::OctetEvents4 generic(prepared.octet4,
                                          camera.point4(RAY.origin), AMBIENT,
                                          RAY.interval, prepared.footprint);
                world =
                    Raycast::shade_events(generic, RAY.interval,
                                          prepared.limits, prepared.appearance);
                const SDF::OctetEvents4 events(
                    prepared.octet4, prepared.octet4_projection, DIRECTION,
                    camera.radial_start, near, prepared.footprint);
                auto copy = events;
                projected = Raycast::shade_events(
                    copy, RAY.interval, prepared.limits, prepared.appearance);
                traced = SDF::OctetTrace::trace_events(
                    events, RAY.interval, prepared.limits, prepared.appearance);
              } else {
                const Raycast::Ray WORLD{camera.point3(RAY.origin),
                                         {AMBIENT[0], AMBIENT[1], AMBIENT[2]},
                                         RAY.interval};
                SDF::OctetEvents generic(prepared.octet, WORLD,
                                         prepared.footprint);
                world =
                    Raycast::shade_events(generic, RAY.interval,
                                          prepared.limits, prepared.appearance);
                const SDF::OctetEvents events(prepared.octet_projection,
                                              DIRECTION, camera.radial_start,
                                              near, prepared.footprint);
                auto copy = events;
                projected = Raycast::shade_events(
                    copy, RAY.interval, prepared.limits, prepared.appearance);
                traced = SDF::OctetTrace::trace_events(
                    events, RAY.interval, prepared.limits, prepared.appearance);
              }
              expect_same(projected, world);
              expect_premultiplied(traced, projected);
              const auto SHADED =
                  domain == Domain::SLICE_4D
                      ? Trace::shade<true>(DIRECTION, prepared)
                      : Trace::shade<false>(DIRECTION, prepared);
              HS_EXPECT_EQ(SHADED.status, traced.status);
              // The 4D trace sums plane positions in its canonical frame's
              // order, which moves crossing distances by an ulp.
              if (domain == Domain::SLICE_4D) {
                HS_EXPECT_NEAR(SHADED.color.r, traced.color.r, 1);
                HS_EXPECT_NEAR(SHADED.color.g, traced.color.g, 1);
                HS_EXPECT_NEAR(SHADED.color.b, traced.color.b, 1);
              } else {
                HS_EXPECT_EQ(SHADED.color, traced.color);
              }
              lit += traced.color != Pixel{};
            }
          }
        }
      }
    }
  }
  HS_EXPECT_GT(lit, size_t{0});
}

/**
 * @brief The canonical-frame 4D octet trace against OctetEvents4 over random
 *        flights, including its plane-share crossing skip.
 */
inline void test_octet_4d_canonical_trace() {
  using Effect = HyperLatticeWhiteBox::Effect;
  reset_globals();
  Effect effect;
  effect.init();
  Trace::Settings settings;
  settings.palette = HyperLatticeWhiteBox::depth_palette(effect);
  settings.pixel_half_angle = HL::pixel_half_angle<288, 144>();
  settings.domain = Raycast::SamplingDomain::SLICE_4D;
  uint32_t state = 0x0c7e74du;
  const auto uniform = [&state](float low, float high) {
    state = state * 1664525u + 1013904223u;
    return low + (high - low) * static_cast<float>(state >> 8) * 0x1p-24f;
  };
  size_t compared = 0, lit = 0, differing = 0;
  for (int frame = 0; frame < 40; ++frame) {
    settings.cell_size = uniform(.6f, 3.0f);
    settings.wire_radius = uniform(.015f, .09f) * settings.cell_size;
    settings.far_distance = uniform(3.0f, 9.0f);
    settings.aa_strength = uniform(0.0f, 2.0f);
    settings.radial_start = frame % 3 ? 0.0f : uniform(0.0f, 1.5f);
    settings.center = {
        {uniform(-3, 3), uniform(-3, 3), uniform(-3, 3), uniform(-3, 3)}};
    settings.embedding = math::Mat4::identity();
    for (int a = 0; a < 4; ++a)
      for (int b = a + 1; b < 4; ++b)
        math::rotate_plane(settings.embedding, a, b, uniform(0.0f, 6.3f));
    const auto prepared = prepare_traced(settings);
    HS_EXPECT_TRUE(prepared.valid);
    const auto &camera = prepared.camera;
    for (int ray = 0; ray < 400; ++ray) {
      const math::Vector DIRECTION =
          math::Vector{uniform(-1, 1), uniform(-1, 1), uniform(-1, 1)}
              .normalized();
      const SDF::OctetEvents4 events(
          prepared.octet4, prepared.octet4_projection, DIRECTION,
          camera.radial_start, camera.interval.near, prepared.footprint);
      const auto EXPECTED = SDF::OctetTrace::trace_events(
          events, camera.interval, prepared.limits, prepared.appearance);
      const auto ACTUAL = SDF::OctetTrace::trace_4d(
          DIRECTION, camera, prepared.octet4, prepared.octet4_projection,
          prepared.footprint, prepared.limits, prepared.appearance,
          *prepared.crossings);
      HS_EXPECT_EQ(ACTUAL.status, EXPECTED.status);
      const int DELTA = std::max({abs(ACTUAL.color.r - EXPECTED.color.r),
                                  abs(ACTUAL.color.g - EXPECTED.color.g),
                                  abs(ACTUAL.color.b - EXPECTED.color.b)});
      differing += DELTA > 1;
      HS_EXPECT_LE(DELTA, 16);
      ++compared;
      lit += EXPECTED.color != Pixel{};
    }
  }
  HS_EXPECT_GT(lit, compared / 10);
  // Counted crossings and the fixed-point bound move a far-boundary
  // crossing or a rounding tie in a few rays by a few codes.
  HS_EXPECT_LE(differing, compared / 4000);
}

inline void test_traced_presets() {
  using Effect = HyperLatticeWhiteBox::Effect;
  reset_globals();
  Effect effect;
  effect.init();
  auto invalid = effect.serialize_parameters();
  invalid.params.pattern = static_cast<Effect::Pattern>(255);
  HS_EXPECT_FALSE(effect.restore_parameters(invalid));
  effect.setAnimationsPaused(true);
  const auto initial = effect.serialize_parameters();
  for (size_t i : {Effect::OCTET_PRESET_INDEX, Effect::OCTET_4D_PRESET_INDEX}) {
    const bool SLICE = i == Effect::OCTET_4D_PRESET_INDEX;
    HL::FrameState frame{};
    frame.params = Effect::preset(i).params;
    frame.params.sphere_radius = .7f;
    frame.depth_palette = HyperLatticeWhiteBox::depth_palette(effect);
    const auto before = prepare_traced(
        HyperLatticeDetail::trace_settings(frame, {{.4f, .7f, .2f, .8f}}));
    frame.params.cell_size *= 2;
    const auto after = prepare_traced(
        HyperLatticeDetail::trace_settings(frame, {{.4f, .7f, .2f, .8f}}));
    HS_EXPECT_TRUE(before.valid && after.valid);
    const auto wrong_domain = !SLICE
                                  ? Trace::shade<true>(math::X_AXIS, before)
                                  : Trace::shade<false>(math::X_AXIS, before);
    HS_EXPECT_EQ(wrong_domain.status, Raycast::TraceStatus::INVALID_QUERY);
    HS_EXPECT_EQ(before.camera.domain, !SLICE
                                           ? Raycast::SamplingDomain::SPATIAL_3D
                                           : Raycast::SamplingDomain::SLICE_4D);
    if (SLICE) {
      HS_EXPECT_NE(before.camera.center[3], 0);
      frame.rotation_phase[3] = .7f;
      const auto rotated = prepare_traced(
          HyperLatticeDetail::trace_settings(frame, {{.4f, .7f, .2f, .8f}}));
      HS_EXPECT_TRUE(rotated.valid);
      HS_EXPECT_NE(rotated.camera.point4(math::X_AXIS)[3],
                   after.camera.point4(math::X_AXIS)[3]);
    }
    for (int axis = 0; axis < 4; ++axis)
      HS_EXPECT_EQ(before.camera.center[axis], after.camera.center[axis]);
    HS_EXPECT_EQ(before.camera.radial_start, .7f);
    HS_EXPECT_EQ(before.camera.radial_start, after.camera.radial_start);
    HS_EXPECT_FALSE(Effect::PRESET_IDS[i].starts_with("experimental-"));
    auto selected = initial;
    selected.params = Effect::preset(i).params;
    HS_EXPECT_TRUE(effect.restore_parameters(selected));
    const auto *pattern = required_param(effect, "Pattern");
    if (!pattern)
      return;
    HS_EXPECT_EQ(pattern->get(), 1);
    const auto *view = required_param(effect, "View");
    if (!view)
      return;
    HS_EXPECT_EQ(view->get(), float(SLICE));
    const auto *spin_4d = required_param(effect, "4D Spin");
    if (!spin_4d)
      return;
    HS_EXPECT_EQ(spin_4d->readonly, !SLICE);
    const auto *lattice_planes = required_param(effect, "Lattice Planes");
    if (!lattice_planes)
      return;
    HS_EXPECT_TRUE(lattice_planes->readonly);
    const auto *wire_radius = required_param(effect, "Wire Radius");
    if (!wire_radius)
      return;
    HS_EXPECT_FALSE(wire_radius->readonly);
    const auto *aa_strength = required_param(effect, "AA Strength");
    if (!aa_strength)
      return;
    HS_EXPECT_FALSE(aa_strength->readonly);
    std::vector<Pixel> previous;
    for (int frame = 0; frame < 2; ++frame) {
      effect.draw_frame();
      effect.advance_display();
      int lit = 0, changed = 0;
      for (int y = 0; y < 20; ++y)
        for (int x = 0; x < 96; ++x) {
          const auto PIXEL = effect.get_pixel(x, y);
          lit += PIXEL.r != 0 || PIXEL.g != 0 || PIXEL.b != 0;
          if (frame == 0)
            previous.push_back(PIXEL);
          else
            changed += PIXEL != previous[y * 96 + x];
        }
      HS_EXPECT_GT(lit, 0);
      if (frame)
        HS_EXPECT_GT(changed, 0);
      const auto *unfinished = required_param(effect, "Unfinished Rays");
      if (!unfinished)
        return;
      HS_EXPECT_LE(unfinished->get(), unfinished->max);
    }
    const auto target = selected.params;
    HL::Params blend;
    blend.lerp(initial.params, target, .49f);
    HS_EXPECT_EQ(blend.pattern, initial.params.pattern);
    HS_EXPECT_EQ(blend.cell_size,
                 hs::lerp(initial.params.cell_size, target.cell_size, .49f));
    blend.lerp(initial.params, target, .5f);
    HS_EXPECT_EQ(blend.pattern, target.pattern);
    HS_EXPECT_EQ(blend.cell_size,
                 hs::lerp(initial.params.cell_size, target.cell_size, .5f));
  }
  HS_EXPECT_TRUE(effect.restore_parameters(initial));
  const auto *view = required_param(effect, "View");
  if (!view)
    return;
  HS_EXPECT_EQ(view->get(), 0);
  const auto *lattice_planes = required_param(effect, "Lattice Planes");
  if (!lattice_planes)
    return;
  HS_EXPECT_FALSE(lattice_planes->readonly);
  HS_EXPECT_EQ(effect.updateParameter("Pattern", 1), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.serialize_parameters().params.pattern,
               Effect::Pattern::OCTET);
  HS_EXPECT_FALSE(effect.selectPreset(Effect::PRESET_IDS.size()));
}

inline void test_pattern_view_controls() {
  using Effect = HyperLatticeWhiteBox::Effect;
  reset_globals();
  Effect effect;
  effect.init();
  HS_EXPECT_TRUE(effect.getParameters().find("Configuration") == nullptr);
  const auto *pattern = required_param(effect, "Pattern");
  if (!pattern)
    return;
  const auto *view = required_param(effect, "View");
  if (!view)
    return;
  HS_EXPECT_EQ(pattern->option_count, 3);
  HS_EXPECT_EQ(view->option_count, 2);
  HS_EXPECT_EQ(effect.updateParameter("View", 1), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(pattern->get(), 0);
  HS_EXPECT_EQ(view->get(), 1);
  const auto *spin_4d = required_param(effect, "4D Spin");
  if (!spin_4d)
    return;
  HS_EXPECT_FALSE(spin_4d->readonly);
  HS_EXPECT_EQ(effect.updateParameter("Cell Size", 3), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("Pattern", 0), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.serialize_parameters().params.cell_size, 3);
  HS_EXPECT_EQ(view->get(), 1);
  auto invalid = effect.serialize_parameters();
  invalid.params.pattern = static_cast<Effect::Pattern>(255);
  HS_EXPECT_FALSE(effect.restore_parameters(invalid));
  HS_EXPECT_EQ(pattern->get(), 0);
  HS_EXPECT_EQ(view->get(), 1);
  static_assert(Effect::PRESET_IDS[Effect::OCTET_PRESET_INDEX] ==
                "octet-flight");
  static_assert(Effect::PRESET_IDS[Effect::OCTET_WIDE_PRESET_INDEX] ==
                "octet-wide-flight");
  HS_EXPECT_EQ(std::string_view(pattern->export_options[1]),
               std::string_view("Pattern::OCTET"));
  const auto four_d = effect.serialize_parameters();
  HS_EXPECT_EQ(effect.updateParameter("Near Fade", .37f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("Pattern", 1), ParamSetResult::APPLIED);
  const auto octet4 = effect.serialize_parameters();
  HS_EXPECT_EQ(octet4.params.pattern, Effect::Pattern::OCTET);
  HS_EXPECT_EQ(octet4.params.mode, HL::LatticeMode::FOUR_D_SLICE);
  HS_EXPECT_EQ(octet4.params.cell_size,
               Effect::preset(Effect::OCTET_4D_PRESET_INDEX).params.cell_size);
  HS_EXPECT_EQ(octet4.params.near_fade, .37f);
  const auto *octet_slice_spin_4d = required_param(effect, "4D Spin");
  if (!octet_slice_spin_4d)
    return;
  HS_EXPECT_FALSE(octet_slice_spin_4d->readonly);
  HS_EXPECT_EQ(effect.updateParameter("View", 0), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(pattern->get(), 1);
  HS_EXPECT_EQ(view->get(), 0);
  const auto *octet_shell_spin_4d = required_param(effect, "4D Spin");
  if (!octet_shell_spin_4d)
    return;
  HS_EXPECT_TRUE(octet_shell_spin_4d->readonly);
  HS_EXPECT_EQ(effect.serialize_parameters().params.spin_4d, 0);
  const auto *lattice_planes = required_param(effect, "Lattice Planes");
  if (!lattice_planes)
    return;
  HS_EXPECT_TRUE(lattice_planes->readonly);
  const auto *wire_radius = required_param(effect, "Wire Radius");
  if (!wire_radius)
    return;
  HS_EXPECT_FALSE(wire_radius->readonly);
  const auto *aa_strength = required_param(effect, "AA Strength");
  if (!aa_strength)
    return;
  HS_EXPECT_FALSE(aa_strength->readonly);
  HS_EXPECT_TRUE(effect.restore_parameters(four_d));
  HS_EXPECT_EQ(pattern->get(), 0);
  HS_EXPECT_EQ(view->get(), 1);
  HS_EXPECT_TRUE(effect.restore_parameters(octet4));
  HS_EXPECT_EQ(pattern->get(), 1);
  HS_EXPECT_EQ(view->get(), 1);
  HS_EXPECT_EQ(effect.updateParameter("Pattern", 0), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(view->get(), 1);
  HS_EXPECT_EQ(
      effect.serialize_parameters().params.wire_radius,
      Effect::preset(Effect::HYPERCUBE_PRESET_INDEX).params.wire_radius);
  const auto *four_d_lattice_planes = required_param(effect, "Lattice Planes");
  if (!four_d_lattice_planes)
    return;
  HS_EXPECT_FALSE(four_d_lattice_planes->readonly);
}

/** @brief Pins normal pattern views and rejects unsupported numeric IDs. */
inline void test_regular_patterns() {
  using Effect = HyperLatticeWhiteBox::Effect;
  static_assert(static_cast<uint8_t>(Effect::Pattern::CUBIC_WIRE) == 0);
  static_assert(static_cast<uint8_t>(Effect::Pattern::OCTET) == 1);
  static_assert(static_cast<uint8_t>(Effect::Pattern::SHELLS) == 6);
  static_assert(Effect::PARAMETER_SCHEMA_VERSION == 15);
  reset_globals();
  Effect effect;
  effect.init();
  effect.setAnimationsPaused(true);
  HS_EXPECT_TRUE(effect.getParameters().find("Shear") == nullptr);
  HS_EXPECT_TRUE(effect.getParameters().find("Stretch") == nullptr);
  const auto *pattern = required_param(effect, "Pattern");
  if (!pattern)
    return;
  constexpr int64_t VALUES[] = {0, 1, 6};
  HS_EXPECT_EQ(pattern->option_count, static_cast<int>(std::size(VALUES)));
  for (size_t i = 0; i < std::size(VALUES); ++i)
    HS_EXPECT_EQ(pattern->option_values[i], VALUES[i]);
  for (const auto &configuration : Effect::CONFIGURATIONS) {
    auto snapshot = effect.serialize_parameters();
    snapshot.params =
        Effect::pattern_defaults(configuration.pattern, configuration.domain);
    HS_EXPECT_TRUE(effect.restore_parameters(snapshot));
    HS_EXPECT_EQ(Effect::configuration_id(snapshot.params), configuration.id);
    const auto *view = required_param(effect, "View");
    if (!view)
      return;
    HS_EXPECT_FALSE(view->readonly);
    const bool SHELLS = configuration.pattern == Effect::Pattern::SHELLS;
    const auto *shell_radius = required_param(effect, "Shell Radius");
    if (!shell_radius)
      return;
    HS_EXPECT_EQ(shell_radius->readonly, !SHELLS);
    const auto *wire_radius = required_param(effect, "Wire Radius");
    if (!wire_radius)
      return;
    HS_EXPECT_EQ(wire_radius->readonly, SHELLS);
    effect.draw_frame();
    effect.advance_display();
    size_t lit = 0;
    for (int y = 0; y < 20; ++y)
      for (int x = 0; x < 96; ++x) {
        const auto pixel = effect.get_pixel(x, y);
        lit += pixel.r || pixel.g || pixel.b;
      }
    HS_EXPECT_GT(lit, size_t{0});
    HS_EXPECT_TRUE(effect.restore_parameters(snapshot));
    const auto current = effect.serialize_parameters();
    for (int retired : {2, 3, 4, 5}) {
      auto invalid = current;
      invalid.params.pattern = static_cast<Effect::Pattern>(retired);
      HS_EXPECT_FALSE(Effect::supported_combination(invalid.params));
      HS_EXPECT_FALSE(effect.restore_parameters(invalid));
      HS_EXPECT_EQ(effect.updateParameter("Pattern", retired),
                   ParamSetResult::INADMISSIBLE);
      HS_EXPECT_EQ(effect.serialize_parameters().params.pattern,
                   current.params.pattern);
    }
    auto invalid_radius = current;
    invalid_radius.params.shell_radius = .33f;
    HS_EXPECT_FALSE(effect.restore_parameters(invalid_radius));
    auto old_schema = current;
    old_schema.schema_version = Effect::PARAMETER_SCHEMA_VERSION - 1;
    HS_EXPECT_FALSE(effect.restore_parameters(old_schema));
    HS_EXPECT_TRUE(effect.restore_parameters(current));
    HL::Params blend;
    const auto start = Effect::preset(Effect::CUBIC_PRESET_INDEX).params;
    blend.lerp(start, snapshot.params, .49f);
    HS_EXPECT_EQ(blend.pattern, start.pattern);
    blend.lerp(start, snapshot.params, .5f);
    HS_EXPECT_EQ(blend.pattern, snapshot.params.pattern);
    HS_EXPECT_EQ(
        blend.shell_radius,
        hs::lerp(start.shell_radius, snapshot.params.shell_radius, .5f));
    auto common = current;
    common.params = start;
    common.params.pattern = configuration.pattern;
    common.params.mode = configuration.domain;
    common.params.speed = common.params.spin_3d = common.params.spin_4d = 0;
    HS_EXPECT_TRUE(effect.restore_parameters(common));
    HyperLatticeWhiteBox::center_camera(effect);
    effect.draw_frame();
    effect.advance_display();
    if (configuration.pattern != Effect::Pattern::CUBIC_WIRE) {
      std::vector<Pixel> actual;
      for (int y = 0; y < 20; ++y)
        for (int x = 0; x < 96; ++x)
          actual.push_back(effect.get_pixel(x, y));
      const auto frame = HyperLatticeWhiteBox::frame(effect);
      const auto prepared = HyperLatticeDetail::prepare_trace(frame);
      {
        Canvas canvas(effect);
        const auto draw_cubic = [&]<bool SLICE>() {
          Scan::Shader::draw<96, 20, 1>(canvas, [&](const math::Vector &view) {
            return HyperLatticeDetail::Renderer<SLICE, 2>::shade_premultiplied(
                view, prepared);
          });
        };
        if (Effect::uses_specialized_slice(common.params))
          draw_cubic.template operator()<true>();
        else
          draw_cubic.template operator()<false>();
      }
      effect.advance_display();
      size_t differing = 0;
      for (int y = 0; y < 20; ++y)
        for (int x = 0; x < 96; ++x)
          differing += effect.get_pixel(x, y) != actual[y * 96 + x];
      HS_EXPECT_GT(differing, size_t{0});
    }
  }
  HS_EXPECT_EQ(effect.updateParameter("Pattern", 6), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.serialize_parameters().params.pattern,
               Effect::Pattern::SHELLS);
  HS_EXPECT_EQ(pattern->get(), 6);
  HS_EXPECT_EQ(std::string_view(pattern->export_options[2]),
               std::string_view("Pattern::SHELLS"));
}

inline void test_speed_range() {
  using Effect = HyperLatticeWhiteBox::Effect;
  reset_globals();
  Effect effect;
  effect.init();
  const auto *speed = required_param(effect, "Speed");
  if (!speed)
    return;
  HS_EXPECT_EQ(speed->max, 0.3f);
  HS_EXPECT_EQ(effect.updateParameter("Speed", 0.3f), ParamSetResult::APPLIED);
  auto snapshot = effect.serialize_parameters();
  HS_EXPECT_EQ(snapshot.params.speed, 0.3f);
  HS_EXPECT_TRUE(effect.restore_parameters(snapshot));
  snapshot.params.speed = 0.31f;
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
}

inline void test_octet_continuous_flight() {
  using Effect = HyperLatticeWhiteBox::Effect;
  for (size_t preset :
       {Effect::OCTET_PRESET_INDEX, Effect::OCTET_WIDE_PRESET_INDEX}) {
    reset_globals();
    Effect effect;
    effect.init();
    HS_EXPECT_TRUE(effect.selectPreset(preset));
    auto &params = HyperLatticeWhiteBox::params(effect);
    params.spin_3d = params.spin_4d = 0;
    const float INITIAL_SPEED = params.speed;
    const auto initial = HyperLatticeWhiteBox::trace_center(effect);
    const float period = std::sqrt(2.0f) * params.cell_size;
    math::Vec4 previous = initial;
    math::Vec4 increment{};
    int wraps = 0;
    for (int frame = 0; frame < 2000; ++frame) {
      HyperLatticeWhiteBox::advance_state(effect);
      const auto current = HyperLatticeWhiteBox::trace_center(effect);
      for (int axis = 0; axis < 4; ++axis) {
        float delta = current[axis] - previous[axis];
        if (delta < 0) {
          delta += period;
          ++wraps;
        }
        HS_EXPECT_GT(delta, 0);
        if (frame == 0)
          increment[axis] = delta;
        else
          HS_EXPECT_NEAR(delta, increment[axis], 3e-7f);
      }
      previous = current;
    }
    HS_EXPECT_GT(wraps, 4);
    params.speed = 0;
    HyperLatticeWhiteBox::advance_state(effect);
    for (int axis = 0; axis < 4; ++axis)
      HS_EXPECT_EQ(HyperLatticeWhiteBox::trace_center(effect)[axis],
                   previous[axis]);
    params.speed = 0.3f;
    params.cell_size = .25f;
    const float small_period = std::sqrt(2.0f) * params.cell_size;
    HyperLatticeWhiteBox::advance_state(effect);
    const auto fast = HyperLatticeWhiteBox::trace_center(effect);
    for (int axis = 0; axis < 4; ++axis) {
      const float expected =
          std::fmod(previous[axis] + increment[axis] * (.3f / INITIAL_SPEED),
                    small_period);
      HS_EXPECT_NEAR(fast[axis], expected, 3e-5f);
      HS_EXPECT_GE(fast[axis], 0);
      HS_EXPECT_LT(fast[axis], small_period);
    }
  }
  for (float cell_size : {.25f, 1.5f, 10.0f}) {
    const float period = std::sqrt(2.0f) * cell_size;
    SDF::OctetFramework octet{cell_size, .055f * cell_size};
    SDF::OctetFramework4 octet4{cell_size, .055f * cell_size};
    for (math::Vec4 point : {math::Vec4{{.21f, .39f, .67f, .14f}},
                             math::Vec4{{-.45f, .62f, -.31f, .87f}}}) {
      const auto sample4 = octet4.sample(point);
      const auto sample3 = octet.sample({point[0], point[1], point[2]});
      for (int axis = 0; axis < 4; ++axis) {
        auto translated = point;
        translated[axis] += period;
        HS_EXPECT_NEAR(octet4.sample(translated).field, sample4.field, 2e-6f);
        if (axis < 3)
          HS_EXPECT_NEAR(
              octet.sample({translated[0], translated[1], translated[2]}).field,
              sample3.field, 2e-6f);
      }
    }
  }
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
  test_configuration_adoption_and_snapshots();
  test_dimension_dropdown_and_mode_lerp();
  test_single_shell();
  test_family_segues();
  test_octet_prepared_projection();
  test_octet_4d_canonical_trace();
  test_regular_patterns();
  test_octet_continuous_flight();
  test_traced_presets();
  test_pattern_view_controls();
  test_speed_range();
  return fixture.result();
}

} // namespace hyper_lattice_tests
} // namespace hs_test
