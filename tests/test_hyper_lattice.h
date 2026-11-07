/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief HyperLattice rendering regressions.
 */
#pragma once

#include "effects/HyperLattice.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/pixel_test_util.h"

namespace hs_test {
namespace hyper_lattice_tests {

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
/** @brief Composites lattice crossings along a view ray. */
template <bool SLICE = false, uint8_t SHELLS = 0>
Color4 shade(const math::Vector &view, const PreparedTrace &prepared) {
  return HyperLatticeDetail::composite_crossings<SLICE, SHELLS>(
             view, {prepared, prepared.appearance, &crossings})
      .finish();
}
} // namespace HL

struct HyperLatticeWhiteBox {
  using Effect = HyperLattice<96, 20>;

  static math::Vec4 origin(const Effect &effect) { return effect.origin; }
  static std::array<float, 6> rotation_phase(const Effect &effect) {
    return effect.rotation_phase;
  }
  static HL::Params &params(Effect &effect) { return effect.params; }
  static float iridescent_phase(const Effect &effect) {
    return effect.iridescent_phase;
  }
  static Pixel iridescent_color(const Effect &effect, float amount) {
    return effect.iridescent_palette.view().get(amount).color;
  }
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
  static Trace::Prepared prepare_traced(const Effect &effect,
                                        const HL::FrameState &frame) {
    return effect.prepare_traced(frame, trace_center(effect));
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

  /** @brief Render state and traced camera of the latest drawn frame. */
  struct Drawn {
    HL::FrameState frame; ///< Render inputs.
    math::Vec4 center;    ///< Traced camera position.
  };
  /** @brief The latest drawn frame of an effect. */
  template <int W, int H> static Drawn drawn(const HyperLattice<W, H> &effect) {
    return {effect.history.previous.state, effect.history.previous.center};
  }
  /** @brief Live parameters of an effect at any resolution. */
  template <int W, int H> static HL::Params &live(HyperLattice<W, H> &effect) {
    return effect.params;
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

/**
 * @brief Cubic and bounded 4D metrics report the expected clearance and free
 * axis.
 */
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
  frame.params.mode = HL::LatticeMode::SLICE_4D;
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

  frame.params.mode = HL::LatticeMode::SPATIAL_3D;
  const math::Vec4 cubic = HL::prepare_trace(frame).world_to_lattice.apply(
      {{1.0f, 0.0f, 0.0f, 0.0f}});
  HS_EXPECT_EQ(cubic[3], 0.0f);
}

inline void test_dimensional_rotation_wrap_is_continuous() {
  HL::FrameState frame{};
  frame.params.mode = HL::LatticeMode::SLICE_4D;
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
       {HL::LatticeMode::SPATIAL_3D, HL::LatticeMode::SLICE_4D}) {
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
          if (mode == HL::LatticeMode::SLICE_4D) {
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

/** @brief The final shell layer fades linearly to zero at the horizon. */
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

inline void test_iridescent_controls_drive_palette() {
  reset_globals();
  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  effect.setAnimationsPaused(true);
  HS_EXPECT_EQ(effect.updateParameter("Iridescent Speed", 0.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("Iridescent Frequency", 1.0f),
               ParamSetResult::APPLIED);
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_EQ(HyperLatticeWhiteBox::iridescent_phase(effect), 0.0f);
  const Pixel initial = HyperLatticeWhiteBox::iridescent_color(effect, .25f);
  HS_EXPECT_EQ(effect.updateParameter("Iridescent Frequency", 2.0f),
               ParamSetResult::APPLIED);
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_EQ(HyperLatticeWhiteBox::iridescent_phase(effect), 0.0f);
  HS_EXPECT_NE(HyperLatticeWhiteBox::iridescent_color(effect, .25f), initial);
  HS_EXPECT_EQ(effect.updateParameter("Iridescent Speed", -.1f),
               ParamSetResult::APPLIED);
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_NEAR(HyperLatticeWhiteBox::iridescent_phase(effect),
                 math::TWO_PI_F - .1f, 1e-6f);
  HS_EXPECT_EQ(effect.updateParameter("Pattern", 1.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("View", 1.0f), ParamSetResult::APPLIED);
  const auto saved = effect.serialize_parameters();
  HS_EXPECT_EQ(saved.params.iridescent_frequency, 2.0f);
  HS_EXPECT_EQ(saved.params.iridescent_speed, -.1f);
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

  run_frames(effect, 480);
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

/**
 * @brief Cell size scales hit distance, and the surface origin advances the
 * ray.
 */
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

/**
 * @brief A rotated slice reaches the expected hyperplane distance and free
 * axis.
 */
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
 * @details The depth palette is baked through cbrtf/powf in the OKLab gamut
 * search, whose last bits differ between libm builds.
 */
constexpr uint16_t MAX_SHADE_CHANNEL_DELTA = 16;

/**
 * @brief Scores rendered samples against a golden table.
 * @param label Context label for out-of-band channels.
 * @param rendered Samples in fold order, one row per (preset, direction).
 * @param golden Expected sample rows.
 * @param count Row count of both tables.
 * @param per_preset Rows per preset, so a failing row names preset and sample.
 */
inline void expect_shade_samples(const char *label, const ShadeSample *rendered,
                                 const ShadeSample *golden, size_t count,
                                 size_t per_preset) {
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
      const Color4 color = HL::shade(DIRECTIONS[sample], prepared);
      rendered[row++] = {color.color.r, color.color.g, color.color.b,
                         frac_to_q16(color.alpha)};
    }
  }
  expect_shade_samples("render_signature", rendered, GOLDEN, std::size(GOLDEN),
                       std::size(DIRECTIONS));
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
          const Color4 exact = HL::shade(direction, prepared);
          lit += exact.alpha > 0.0f ? 1 : 0;
          const Color4 specialized =
              HL::shade<true, SHELL_COUNT>(direction, prepared);
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
 * @brief Pins specialized 4D-slice composite_crossings samples.
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
    const auto prepared = HL::prepare_trace(context);
    for (size_t sample = 0; sample < std::size(DIRECTIONS); ++sample) {
      const Color4 color = HL::shade<true, 2>(DIRECTIONS[sample], prepared);
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
                       std::size(GOLDEN), std::size(DIRECTIONS));
#endif
}

/**
 * @brief Renders forty frames of a preset with column-time sampling disabled.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @param preset Preset index.
 * @return Pin folded over the first two frames and the last.
 */
template <int W, int H> PixelGolden instantaneous_render(size_t preset) {
  reset_globals();
  HyperLattice<W, H> effect;
  effect.init();
  if (preset)
    HS_EXPECT_TRUE(effect.selectPreset(preset));
  Spacetime::set_flat_surface(true);
  PixelGolden pin;
  for (int frame = 0; frame < 40; ++frame) {
    effect.draw_frame();
    effect.advance_display();
    if (frame < 2 || frame == 39)
      fold_pixel_golden<W, H>(pin, effect);
  }
  Spacetime::set_flat_surface(false);
  return pin;
}

/** @brief Every preset renders as recorded with column-time sampling off. */
inline void test_instantaneous_render_pins() {
  static constexpr RecordedGolden SMALL[] = {
      {{0xa94c916ce098b710ull,
        {{15879938ull, 33806139ull, 49901045ull}, 929681480ull, 4747878914ull}},
       {0xe83aa7de601edb1bull,
        {{15879941ull, 33806144ull, 49901036ull},
         929681465ull,
         4747879454ull}}},
      {{0x48dacd6716f2e38aull,
        {{13620838ull, 27247961ull, 41856810ull}, 788158262ull, 3993009095ull}},
       {0x54b396c8d6c90c84ull,
        {{13620839ull, 27247967ull, 41856810ull},
         788158425ull,
         3993009814ull}}},
      {{0xe4e4b214301d1fbbull,
        {{14016703ull, 28343629ull, 38782553ull}, 794147276ull, 3840080868ull}},
       {0xb60da53ffc9d96e7ull,
        {{14016703ull, 28343635ull, 38782547ull},
         794147220ull,
         3840080922ull}}},
      {{0xe18b869cecc87b1aull,
        {{6877621ull, 15388262ull, 25339230ull}, 452146061ull, 2245807307ull}},
       {0x216ead7f6a941b29ull,
        {{6877623ull, 15388258ull, 25339219ull}, 452145868ull, 2245806308ull}}},
      {{0x642134b4d9e14f12ull,
        {{15574978ull, 32052190ull, 46757765ull}, 865498826ull, 4441257817ull}},
       {0xcd5c27f243cbf481ull,
        {{15574982ull, 32052282ull, 46757771ull},
         865498963ull,
         4441262711ull}}},
      {{0x20ca3708641fc294ull,
        {{15094188ull, 29775474ull, 49246986ull}, 950669762ull, 4630249745ull}},
       {0xa19f572ef67c31ecull,
        {{15094181ull, 29775480ull, 49246981ull},
         950669704ull,
         4630249183ull}}},
      {{0x832cd9176e473075ull,
        {{6022067ull, 12585294ull, 17867684ull}, 349572479ull, 1658865216ull}},
       {0x75701bceb98520fbull,
        {{6022070ull, 12585302ull, 17867698ull}, 349572695ull, 1658866632ull}}},
      {{0xf44109a19d638e14ull,
        {{13884469ull, 26887004ull, 33744219ull}, 731916001ull, 3728240502ull}},
       {0x9231ccaafe6e82acull,
        {{13884469ull, 26887017ull, 33744232ull},
         731916192ull,
         3728241718ull}}},
      {{0x234d1dcf958db5aeull,
        {{1683359ull, 3121210ull, 4839967ull}, 87818621ull, 461967067ull}},
       {0x42772e7edbe6e402ull,
        {{1683363ull, 3121211ull, 4839966ull}, 87818636ull, 461967398ull}}},
  };
  static constexpr RecordedGolden FULL[] = {
      {{0x6efc9ab887e3687eull,
        {{286420393ull, 604026957ull, 880428212ull},
         125165022407ull,
         255031996507ull}},
       {0x48531a97457d671cull,
        {{286420451ull, 604027075ull, 880428028ull},
         125165022266ull,
         255031987260ull}}},
      {{0x9a6c8e560b9e215eull,
        {{226325234ull, 466226264ull, 694760587ull},
         99255641409ull,
         202137900576ull}},
       {0xa4194985b8074047ull,
        {{226325227ull, 466226402ull, 694760504ull},
         99255645653ull,
         202137915890ull}}},
      {{0xcfd58f213c8d43f5ull,
        {{56724327ull, 114293350ull, 149493309ull},
         23553561504ull,
         45995717139ull}},
       {0x46f34c3e3e99c108ull,
        {{56724362ull, 114293392ull, 149493300ull},
         23553565214ull,
         45995734529ull}}},
      {{0x33127babce0433d0ull,
        {{98796890ull, 216440996ull, 339016827ull},
         47134284500ull,
         92777983011ull}},
       {0x51f59dbbb558b7dfull,
        {{98796887ull, 216441042ull, 339016800ull},
         47134287388ull,
         92777982356ull}}},
      {{0xf471ca401bf88bb0ull,
        {{203014955ull, 412599828ull, 570236991ull},
         82903792178ull,
         169550270676ull}},
       {0x1718f69a430a4cf2ull,
        {{203014970ull, 412599843ull, 570236906ull},
         82903784158ull,
         169550257129ull}}},
      {{0x07f3bda104eefd64ull,
        {{103486630ull, 181217087ull, 244413948ull},
         41465102299ull,
         81990154845ull}},
       {0xa11a4dd81e2755fdull,
        {{103486623ull, 181217097ull, 244413924ull},
         41465101267ull,
         81990151803ull}}},
      {{0xf89b62bf071d6a88ull,
        {{136909400ull, 281503680ull, 371235824ull},
         56122681005ull,
         106198756306ull}},
       {0xa41e404771940ce1ull,
        {{136909454ull, 281503670ull, 371235770ull},
         56122675450ull,
         106198748654ull}}},
      {{0x198d7ceb0de43b30ull,
        {{338516795ull, 661866865ull, 832048094ull},
         133695350754ull,
         278181293068ull}},
       {0x2926cc5f952d7723ull,
        {{338516979ull, 661866986ull, 832048186ull},
         133695373713ull,
         278181343504ull}}},
      {{0x2c8abb9d70d06307ull,
        {{28033812ull, 50416081ull, 63394720ull},
         8643288134ull,
         19672679156ull}},
       {0xa392e46781270164ull,
        {{28033859ull, 50416082ull, 63394693ull},
         8643285665ull,
         19672684866ull}}},
  };
  static_assert(std::size(SMALL) == HyperLattice<96, 20>::PRESET_IDS.size());
  static_assert(std::size(FULL) == HyperLattice<96, 20>::PRESET_IDS.size());
  for (size_t preset = 0; preset < std::size(SMALL); ++preset) {
    HS_CONTEXT("preset", static_cast<long long>(preset));
    const PixelGolden small = instantaneous_render<96, 20>(preset);
    const PixelGolden full = instantaneous_render<288, 144>(preset);
    expect_pixel_golden("HyperLattice 96x20", small, SMALL[preset]);
    expect_pixel_golden("HyperLattice 288x144", full, FULL[preset]);
  }
}

/** @brief Unwrapped displacement between two periodic coordinates.
 *  @param from Earlier coordinate.
 *  @param to Later coordinate.
 *  @param period Coordinate period.
 *  @return The shortest signed displacement. */
inline float unwrapped(float from, float to, float period) {
  const float DELTA = to - from;
  return DELTA - period * std::round(DELTA / period);
}

/**
 * @brief Render state part way between two drawn frames, evaluated
 *        independently of the effect's sampler.
 * @param from Previous drawn frame.
 * @param to Current drawn frame.
 * @param fraction Fraction of the frame.
 * @return State and traced camera at that time.
 */
inline HyperLatticeWhiteBox::Drawn
state_between(const HyperLatticeWhiteBox::Drawn &from,
              const HyperLatticeWhiteBox::Drawn &to, float fraction) {
  HyperLatticeWhiteBox::Drawn out = to;
  out.frame.params.lerp(from.frame.params, to.frame.params, fraction);
  const auto &p = to.frame.params;
  const float PERIOD =
      p.cell_size * (p.pattern == HyperLatticeDetail::Pattern::OCTET
                         ? std::sqrt(2.0f)
                         : 1.0f);
  for (int axis = 0; axis < HL::DIMENSIONS; ++axis) {
    out.frame.origin[axis] =
        from.frame.origin[axis] +
        unwrapped(from.frame.origin[axis], to.frame.origin[axis], 1.0f) *
            fraction;
    out.center[axis] =
        from.center[axis] +
        unwrapped(from.center[axis], to.center[axis], PERIOD) * fraction;
  }
  for (size_t plane = 0; plane < out.frame.rotation_phase.size(); ++plane)
    out.frame.rotation_phase[plane] =
        from.frame.rotation_phase[plane] +
        unwrapped(from.frame.rotation_phase[plane],
                  to.frame.rotation_phase[plane], math::TWO_PI_F) *
            fraction;
  out.frame.gain = hs::lerp(from.frame.gain, to.frame.gain, fraction);
  return out;
}

/**
 * @brief One pixel of a state, traced by the generic renderers.
 * @tparam W,H Canvas size.
 * @param state Render state and traced camera.
 * @param x,y Pixel.
 * @return Premultiplied colour.
 */
template <int W, int H>
inline Pixel traced_pixel(const HyperLatticeWhiteBox::Drawn &state, int x,
                          int y) {
  const math::Vector VIEW = math::pixel_to_vector<W, H>(x, y);
  const auto &params = state.frame.params;
  const bool SLICE = params.mode == HL::LatticeMode::SLICE_4D;
  if (params.pattern == HyperLatticeDetail::Pattern::CUBIC_WIRE) {
    const auto prepared = HyperLatticeDetail::prepare_trace(state.frame);
    return SLICE ? HyperLatticeDetail::Renderer<true>::shade_premultiplied(
                       VIEW, prepared)
                 : HyperLatticeDetail::Renderer<>::shade_premultiplied(
                       VIEW, prepared);
  }
  static SDF::ShellLayerStorage layers;
  auto settings = HyperLatticeDetail::trace_settings(state.frame, state.center);
  settings.shell_layers = &layers;
  const auto prepared = prepare_traced(settings);
  return (SLICE ? Trace::shade<true>(VIEW, prepared)
                : Trace::shade<false>(VIEW, prepared))
      .color;
}

/** @brief Agreement of rendered columns with a time evaluator. */
struct ColumnAgreement {
  int worst = 0;   ///< Largest channel gap, 16-bit codes.
  int changed = 0; ///< Pixels whose column-time state differs from the end.
};

/**
 * @brief Compares a rendered window with each column's state between two
 *        drawn frames.
 * @tparam W,H Canvas size.
 * @param effect Effect after draw_frame() and advance_display().
 * @param from Previous drawn frame.
 * @param to Current drawn frame.
 * @param x0,x1,y0,y1 Window.
 * @return Worst gap and the count of pixels the motion changed.
 */
template <int W, int H>
inline ColumnAgreement column_agreement(const HyperLattice<W, H> &effect,
                                        const HyperLatticeWhiteBox::Drawn &from,
                                        const HyperLatticeWhiteBox::Drawn &to,
                                        int x0, int x1, int y0, int y1) {
  ColumnAgreement out;
  for (int x = x0; x < x1; ++x) {
    const auto STATE =
        state_between(from, to, static_cast<float>(x % (W / 2)) / (W / 2));
    for (int y = y0; y < y1; ++y) {
      const Pixel EXPECTED = traced_pixel<W, H>(STATE, x, y);
      out.worst = std::max<int>(
          out.worst, max_channel_gap(effect.get_pixel(x, y), EXPECTED));
      out.changed += EXPECTED != traced_pixel<W, H>(to, x, y);
    }
  }
  return out;
}

/** @brief Largest per-channel gap between a column-time render and its time
 *  evaluator, in 16-bit codes. */
constexpr int COLUMN_TOLERANCE = 16;

/**
 * @brief Checks a preset's columns against the time evaluator across two
 *        frame seams.
 * @tparam W,H Canvas size.
 * @param preset Preset index.
 * @param x0,x1 Clipped columns; the render margin is zero.
 */
template <int W, int H>
inline void check_column_time(size_t preset, int x0, int x1) {
  using WB = HyperLatticeWhiteBox;
  HS_CONTEXT("preset", static_cast<long long>(preset), x0);
  reset_globals();
  HyperLattice<W, H> effect;
  effect.init();
  effect.set_margin(0);
  effect.set_clip(0, H, x0, x1);
  if (preset)
    HS_EXPECT_TRUE(effect.selectPreset(preset));
  effect.draw_frame();
  effect.advance_display();
  const auto first = WB::drawn(effect);
  const auto primed = column_agreement(effect, first, first, x0, x1, 0, H);
  HS_EXPECT_LE(primed.worst, COLUMN_TOLERANCE);
  run_frames(effect, 40);
  auto from = WB::drawn(effect);
  for (int frame = 0; frame < 2; ++frame) {
    HS_CONTEXT("frame", frame);
    effect.draw_frame();
    effect.advance_display();
    const auto to = WB::drawn(effect);
    const auto agreement = column_agreement(effect, from, to, x0, x1, 0, H);
    HS_EXPECT_LE(agreement.worst, COLUMN_TOLERANCE);
    HS_EXPECT_GT(agreement.changed, 0);
    for (int x : {0, W / 2})
      if (x >= x0 && x < x1)
        for (int y = 0; y < H; ++y)
          HS_EXPECT_LE(max_channel_gap(effect.get_pixel(x, y),
                                       traced_pixel<W, H>(from, x, y)),
                       COLUMN_TOLERANCE);
    from = to;
  }
}

/** @brief Each column traces the state at its sweep time, through the seam
 *  into the next frame, on every trace path. */
inline void test_column_time_matches_time_evaluator() {
  using Effect = HyperLatticeWhiteBox::Effect;
  for (size_t preset :
       {Effect::CUBIC_PRESET_INDEX, Effect::HYPERCUBE_PRESET_INDEX,
        Effect::OCTET_PRESET_INDEX, Effect::OCTET_4D_PRESET_INDEX,
        Effect::SHELL_PRESET_INDEX, Effect::SHELL_4D_PRESET_INDEX})
    check_column_time<96, 20>(preset, 0, 96);
  for (size_t preset :
       {Effect::SHELL_PRESET_INDEX, Effect::SHELL_4D_PRESET_INDEX})
    check_column_time<288, 144>(preset, 140, 150);
}

/**
 * @brief Renders three column-time frames of a preset into a clip band.
 * @param preset Preset index.
 * @param x0,x1,y0,y1 Band.
 * @return The band's pixels.
 */
inline std::vector<Pixel> column_time_band(size_t preset, int x0, int x1,
                                           int y0, int y1) {
  reset_globals();
  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  effect.set_margin(3);
  effect.set_clip(y0, y1, x0, x1);
  if (preset)
    HS_EXPECT_TRUE(effect.selectPreset(preset));
  for (int frame = 0; frame < 3; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  std::vector<Pixel> band;
  for (int y = y0; y < y1; ++y)
    for (int x = x0; x < x1; ++x)
      band.push_back(effect.get_pixel(x, y));
  return band;
}

/** @brief Segment bands crossing both arm boundaries match the single-engine
 *  render. */
inline void test_column_time_segments_match_full_frame() {
  using Effect = HyperLatticeWhiteBox::Effect;
  struct Band {
    int x0, x1, y0, y1;
  };
  const Band bands[] = {{0, 30, 0, 20}, {40, 96, 0, 20}, {45, 51, 3, 17}};
  for (size_t preset :
       {Effect::CUBIC_PRESET_INDEX, Effect::OCTET_PRESET_INDEX}) {
    const auto full = column_time_band(preset, 0, 96, 0, 20);
    for (const Band &b : bands) {
      HS_CONTEXT("band", b.x0, static_cast<long long>(preset));
      const auto band = column_time_band(preset, b.x0, b.x1, b.y0, b.y1);
      size_t i = 0, different = 0;
      for (int y = b.y0; y < b.y1; ++y)
        for (int x = b.x0; x < b.x1; ++x, ++i)
          different += band[i] != full[static_cast<size_t>(y) * 96 + x];
      HS_EXPECT_EQ(different, size_t{0});
    }
  }
}

/**
 * @brief Draws the frame after a change and checks it shows the current
 *        state in every column.
 * @param change Applies the change to the effect.
 */
template <typename Change> inline void check_history_reset(Change &&change) {
  using WB = HyperLatticeWhiteBox;
  reset_globals();
  WB::Effect effect;
  effect.init();
  for (int frame = 0; frame < 2; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  const auto before = WB::drawn(effect);
  change(effect);
  effect.draw_frame();
  effect.advance_display();
  const auto after = WB::drawn(effect);
  const auto held = column_agreement(effect, after, after, 0, 96, 0, 20);
  HS_EXPECT_LE(held.worst, COLUMN_TOLERANCE);
  int moved = 0;
  for (int x = 0; x < 96; ++x)
    for (int y = 0; y < 20; ++y)
      moved += max_channel_gap(effect.get_pixel(x, y),
                               traced_pixel<96, 20>(before, x, y)) >
               COLUMN_TOLERANCE;
  HS_EXPECT_GT(moved, 0);
}

/** @brief Discrete changes, manual edits and manual presets draw their first
 *  frame without motion. */
inline void test_column_time_resets() {
  using WB = HyperLatticeWhiteBox;
  check_history_reset([](WB::Effect &effect) {
    HS_EXPECT_EQ(effect.updateParameter("Pattern", 1), ParamSetResult::APPLIED);
  });
  check_history_reset([](WB::Effect &effect) {
    WB::live(effect).shells = HL::ShellCount::THREE;
  });
  check_history_reset([](WB::Effect &effect) {
    HS_EXPECT_TRUE(effect.selectPreset(WB::Effect::WIDE_PRESET_INDEX));
  });
  check_history_reset([](WB::Effect &effect) {
    HS_EXPECT_EQ(effect.updateParameter("Cell Size", 1.5f),
                 ParamSetResult::APPLIED);
  });
}

/** @brief Disabled sampling draws the current state while history advances,
 *  and resumed sampling starts from the disabled frame. */
inline void test_column_time_toggle() {
  using WB = HyperLatticeWhiteBox;
  reset_globals();
  WB::Effect effect;
  effect.init();
  effect.draw_frame();
  effect.advance_display();
  Spacetime::set_flat_surface(true);
  effect.draw_frame();
  effect.advance_display();
  const auto held = WB::drawn(effect);
  for (int x = 0; x < 96; ++x)
    for (int y = 0; y < 20; ++y)
      HS_EXPECT_LE(max_channel_gap(effect.get_pixel(x, y),
                                   traced_pixel<96, 20>(held, x, y)),
                   COLUMN_TOLERANCE);
  Spacetime::set_flat_surface(false);
  effect.draw_frame();
  effect.advance_display();
  const auto resumed =
      column_agreement(effect, held, WB::drawn(effect), 0, 96, 0, 20);
  HS_EXPECT_LE(resumed.worst, COLUMN_TOLERANCE);
  HS_EXPECT_GT(resumed.changed, 0);
}

/** @brief Columns interpolate continuous parameters, and a shell frame whose
 *  columns straddle a trace threshold traces every column validly. */
inline void test_column_time_parameter_motion() {
  using WB = HyperLatticeWhiteBox;
  {
    reset_globals();
    WB::Effect effect;
    effect.init();
    effect.draw_frame();
    effect.advance_display();
    const auto from = WB::drawn(effect);
    WB::live(effect).cell_size *= 1.05f;
    effect.draw_frame();
    effect.advance_display();
    const auto agreement =
        column_agreement(effect, from, WB::drawn(effect), 0, 96, 0, 20);
    HS_EXPECT_LE(agreement.worst, COLUMN_TOLERANCE);
    HS_EXPECT_GT(agreement.changed, 0);
  }
  constexpr int W = 288, H = 144;
  reset_globals();
  HyperLattice<W, H> effect;
  effect.init();
  effect.set_margin(0);
  effect.set_clip(0, H, 1, 21);
  HS_EXPECT_TRUE(effect.selectPreset(WB::Effect::SHELL_PRESET_INDEX));
  const HL::Params SHELLS = WB::live(effect);
  WB::live(effect).shell_radius = .32f;
  WB::live(effect).far_distance = 16.0f;
  effect.draw_frame();
  effect.advance_display();
  const auto from = WB::drawn(effect);
  WB::live(effect) = SHELLS;
  effect.draw_frame();
  effect.advance_display();
  const auto to = WB::drawn(effect);
  const auto start = prepare_traced(
      HyperLatticeDetail::trace_settings(from.frame, from.center));
  HS_EXPECT_FALSE(start.periodic_shells.single_owner);
  HS_EXPECT_FALSE(start.periodic_shells.march);
  HS_EXPECT_TRUE(
      prepare_traced(HyperLatticeDetail::trace_settings(to.frame, to.center))
          .periodic_shells.single_owner);
  const auto agreement = column_agreement(effect, from, to, 1, 21, 0, H);
  HS_EXPECT_LE(agreement.worst, COLUMN_TOLERANCE);
}

/** @brief Frame equality ignores padding and retains value differences. */
inline void test_frame_equality_ignores_padding() {
  using Endpoint = HyperLatticeDetail::FrameEndpoint;
  using State = HL::FrameState;
  using Params = HL::Params;
  Endpoint from{};
  Endpoint to = from;
  auto *bytes = reinterpret_cast<unsigned char *>(&to.state);
  auto change_padding = [&](size_t begin, size_t end) {
    for (size_t i = begin; i < end; ++i)
      bytes[i] ^= 0xff;
  };
  change_padding(offsetof(Params, pattern) + sizeof(to.state.params.pattern),
                 offsetof(Params, sphere_radius));
  change_padding(offsetof(Params, shells) + sizeof(to.state.params.shells),
                 offsetof(Params, shell_radius));
  change_padding(offsetof(State, pixel_half_angle) + sizeof(float),
                 offsetof(State, depth_palette));
  change_padding(offsetof(State, gain) + sizeof(float),
                 offsetof(State, crossings));
  HS_EXPECT_TRUE(from == to);
  const HyperLatticeDetail::FrameInterval constant(from.state, to.state,
                                                   from.center, {});
  HS_EXPECT_TRUE(constant.params_constant());
  to.state.params.speed += 0.1f;
  HS_EXPECT_FALSE(from == to);
  const HyperLatticeDetail::FrameInterval moving(from.state, to.state,
                                                 from.center, {});
  HS_EXPECT_FALSE(moving.params_constant());
  to = from;
  to.center[3] = 1.0f;
  HS_EXPECT_FALSE(from == to);
}

inline void test_presets_are_valid_and_distinct() {
  using Effect = HyperLattice<96, 20>;
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
  static_assert(CUBIC_PRESET.mode == HL::LatticeMode::SPATIAL_3D);
  static_assert(CUBIC_PRESET.shells == HL::ShellCount::TWO);

  constexpr HL::Params SLICE_PRESET =
      Effect::preset(Effect::HYPERCUBE_PRESET_INDEX).params;
  static_assert(SLICE_PRESET.mode == HL::LatticeMode::SLICE_4D);
  static_assert(SLICE_PRESET.shells == HL::ShellCount::TWO);

  static_assert(SLICE_PRESET.sphere_radius == 0.0f);
  reset_globals();
  Effect effect;
  effect.init();
  std::vector<std::vector<float>> values;
  for (size_t index = 0; index < Effect::PRESET_IDS.size(); ++index) {
    HS_EXPECT_TRUE(effect.selectPreset(index));
    hs_test::expect_params_in_range(effect, "HyperLattice preset");
    std::vector<float> current;
    for (const auto &def : effect.getParameters())
      current.push_back(def.get());
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
  HS_REQUIRED_PARAM(spin_4d, effect, "4D Spin");
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

  HS_EXPECT_TRUE(effect.selectPreset(Effect::CUBIC_PRESET_INDEX));
  effect.setAnimationsPaused(false);
  const auto CUBIC = Effect::preset(Effect::CUBIC_PRESET_INDEX).params;
  for (int frame = 0; frame < Effect::PRESET_DWELL_FRAMES; ++frame)
    HyperLatticeWhiteBox::step_choreography(effect, canvas);
  const int LERP =
      Segue::Preset::frames(Effect::preset(Effect::CUBIC_PRESET_INDEX).segue);
  for (int frame = 1; frame <= LERP; ++frame) {
    HyperLatticeWhiteBox::step_choreography(effect, canvas);
    HS_EXPECT_EQ(HyperLatticeWhiteBox::preset_gain(effect), 1.0f);
    if (frame < LERP) {
      HS_EXPECT_GT(params.cell_size, CUBIC.cell_size);
      HS_EXPECT_LT(params.cell_size, WIDE.cell_size);
    }
  }
  HS_EXPECT_EQ(params.cell_size, WIDE.cell_size);
  HS_EXPECT_EQ(effect.getPresetIndex(), Effect::WIDE_PRESET_INDEX);
}

/** @brief Near Fade exposes live controls and interpolates linearly. */
inline void test_near_fade_control() {
  reset_globals();
  using Effect = HyperLatticeWhiteBox::Effect;
  Effect effect;
  effect.init();
  HS_REQUIRED_PARAM(near_fade, effect, "Near Fade");
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
}

inline void test_dimension_dropdown_and_mode_lerp() {
  reset_globals();
  using Effect = HyperLatticeWhiteBox::Effect;
  Effect effect;
  effect.init();
  HS_REQUIRED_PARAM(cell_size, effect, "Cell Size");
  HS_EXPECT_EQ(cell_size->max, 10.0f);
  HS_REQUIRED_PARAM(far_distance, effect, "Far Distance");

  HS_REQUIRED_PARAM(dimension, effect, "View");
  HS_EXPECT_TRUE(dimension->is_enum());
  HS_EXPECT_EQ(dimension->option_count, 2);
  HS_EXPECT_STREQ(dimension->options[0], "3D perspective");
  HS_EXPECT_STREQ(dimension->options[1], "4D slice");
  HS_EXPECT_STREQ(dimension->export_options[1], "LatticeMode::SLICE_4D");
  HS_EXPECT_EQ(effect.updateParameter("View", 1.0f), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(HyperLatticeWhiteBox::params(effect).mode,
               HL::LatticeMode::SLICE_4D);

  HL::Params start;
  HL::Params target;
  target.mode = HL::LatticeMode::SLICE_4D;
  HL::Params blended;
  blended.lerp(start, target, 0.49f);
  HS_EXPECT_EQ(blended.mode, HL::LatticeMode::SPATIAL_3D);
  blended.lerp(start, target, 0.5f);
  HS_EXPECT_EQ(blended.mode, HL::LatticeMode::SLICE_4D);
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

/** @brief Trace geometry maps each persisted pattern explicitly. */
inline void test_trace_geometry_mapping() {
  static_assert(static_cast<uint8_t>(Trace::Geometry::OCTET) == 0);
  static_assert(static_cast<uint8_t>(Trace::Geometry::SHELLS) == 1);
  HL::FrameState frame{};
  frame.params.pattern = HyperLatticeDetail::Pattern::OCTET;
  HS_EXPECT_EQ(HyperLatticeDetail::trace_settings(frame, {}).geometry,
               Trace::Geometry::OCTET);
  frame.params.pattern = HyperLatticeDetail::Pattern::SHELLS;
  HS_EXPECT_EQ(HyperLatticeDetail::trace_settings(frame, {}).geometry,
               Trace::Geometry::SHELLS);
}

/** @brief Shared trace paths preserve specialization, domain guards and fallback. */
inline void test_shared_trace_path_dispatch() {
  reset_globals();
  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  SDF::ShellLayerStorage layers;
  Trace::Settings settings;
  settings.palette = HyperLatticeWhiteBox::depth_palette(effect);
  settings.shell_layers = &layers;
  settings.far_distance = 3.0f;
  const auto same = [](const Trace::Sample &actual,
                       const Trace::Sample &expected) {
    HS_EXPECT_EQ(actual.status, expected.status);
    HS_EXPECT_EQ(actual.color, expected.color);
  };
  const auto check = [&]<bool SLICE_4D>() {
    settings.domain = SLICE_4D ? Raycast::SamplingDomain::SLICE_4D
                               : Raycast::SamplingDomain::SPATIAL_3D;
    settings.geometry = Trace::Geometry::OCTET;
    auto prepared = prepare_traced(settings);
    HS_EXPECT_TRUE(prepared.valid);
    HS_EXPECT_EQ(Trace::trace_path(prepared), Trace::TracePath::OCTET);
    same(Trace::shade_path<Trace::TracePath::OCTET, SLICE_4D>(math::X_AXIS,
                                                              prepared),
         Trace::shade<SLICE_4D>(math::X_AXIS, prepared));
    settings.geometry = Trace::Geometry::SHELLS;
    for (const float radius : {.1f, .35f, .48f}) {
      settings.shell_radius = radius;
      settings.pixel_half_angle = radius == .48f ? .1f : 0.0f;
      prepared = prepare_traced(settings);
      HS_EXPECT_TRUE(prepared.valid);
      const auto expected_path = radius == .48f ? Trace::TracePath::GENERIC
                                 : !SLICE_4D && radius == .1f
                                     ? Trace::TracePath::SINGLE_SHELLS
                                     : Trace::TracePath::MARCHED_SHELLS;
      HS_EXPECT_EQ(Trace::trace_path(prepared), expected_path);
      if constexpr (SLICE_4D) {
        prepared.periodic_shells.single_owner = true;
        HS_EXPECT_NE(Trace::trace_path(prepared),
                     Trace::TracePath::SINGLE_SHELLS);
        prepared.periodic_shells.single_owner = false;
      }
      for (const auto &view : {math::X_AXIS, math::Y_AXIS,
                               math::Vector(.2f, .3f, 1).normalized()}) {
        const auto actual = Trace::shade<SLICE_4D>(view, prepared);
        same(Trace::shade_path<Trace::TracePath::GENERIC, SLICE_4D>(view,
                                                                    prepared),
             actual);
        if constexpr (!SLICE_4D)
          if (radius == .1f)
            same(Trace::shade_path<Trace::TracePath::SINGLE_SHELLS, false>(
                     view, prepared),
                 actual);
        if (expected_path == Trace::TracePath::MARCHED_SHELLS)
          same(Trace::shade_path<Trace::TracePath::MARCHED_SHELLS, SLICE_4D>(
                   view, prepared),
               actual);
        if (expected_path == Trace::TracePath::GENERIC)
          same(actual,
               (SDF::trace_periodic_shells_march<SLICE_4D ? 4 : 3, true, true>(
                   prepared.periodic_shells, prepared.camera, view,
                   prepared.limits, prepared.appearance, layers)));
      }
      const auto wrong_domain =
          Trace::shade_path<Trace::TracePath::GENERIC, !SLICE_4D>(math::X_AXIS,
                                                                  prepared);
      HS_EXPECT_EQ(wrong_domain.status, Raycast::TraceStatus::INVALID_QUERY);
    }
    prepared.valid = false;
    HS_EXPECT_EQ(Trace::trace_path(prepared), Trace::TracePath::GENERIC);
    const auto invalid = Trace::shade_path<Trace::TracePath::GENERIC, SLICE_4D>(
        math::X_AXIS, prepared);
    HS_EXPECT_EQ(invalid.status, Raycast::TraceStatus::INVALID_QUERY);
  };
  check.template operator()<false>();
  check.template operator()<true>();
}

/** @brief Traced preset cameras preserve origins across scale and reject wrong domains. */
inline void test_traced_preset_camera_setup() {
  using Effect = HyperLatticeWhiteBox::Effect;
  reset_globals();
  Effect effect;
  effect.init();
  effect.setAnimationsPaused(true);
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
  }
}

/** @brief Traced presets expose the matching editable controls and retain fallback state. */
inline void test_traced_preset_control_state() {
  using Effect = HyperLatticeWhiteBox::Effect;
  reset_globals();
  Effect effect;
  effect.init();
  effect.setAnimationsPaused(true);
  const auto initial = effect.serialize_parameters();
  for (size_t i : {Effect::OCTET_PRESET_INDEX, Effect::OCTET_4D_PRESET_INDEX}) {
    const bool SLICE = i == Effect::OCTET_4D_PRESET_INDEX;
    auto selected = initial;
    selected.params = Effect::preset(i).params;
    HS_EXPECT_TRUE(effect.restore_parameters(selected));
    HS_REQUIRED_PARAM(pattern, effect, "Pattern");
    HS_EXPECT_EQ(pattern->get(), 1);
    HS_REQUIRED_PARAM(view, effect, "View");
    HS_EXPECT_EQ(view->get(), float(SLICE));
    HS_REQUIRED_PARAM(spin_4d, effect, "4D Spin");
    HS_EXPECT_EQ(spin_4d->readonly, !SLICE);
    HS_REQUIRED_PARAM(lattice_planes, effect, "Lattice Planes");
    HS_EXPECT_TRUE(lattice_planes->readonly);
    HS_REQUIRED_PARAM(wire_radius, effect, "Wire Radius");
    HS_EXPECT_FALSE(wire_radius->readonly);
    HS_REQUIRED_PARAM(aa_strength, effect, "AA Strength");
    HS_EXPECT_FALSE(aa_strength->readonly);
  }
  HS_EXPECT_TRUE(effect.restore_parameters(initial));
  HS_REQUIRED_PARAM(view, effect, "View");
  HS_EXPECT_EQ(view->get(), 0);
  HS_REQUIRED_PARAM(lattice_planes, effect, "Lattice Planes");
  HS_EXPECT_FALSE(lattice_planes->readonly);
  HS_EXPECT_EQ(effect.updateParameter("Pattern", 1), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.serialize_parameters().params.pattern,
               Effect::Pattern::OCTET);
  HS_EXPECT_FALSE(effect.selectPreset(Effect::PRESET_IDS.size()));
}

/** @brief Traced octet presets render changing, complete frames. */
inline void test_traced_preset_render() {
  using Effect = HyperLatticeWhiteBox::Effect;
  reset_globals();
  Effect effect;
  effect.init();
  effect.setAnimationsPaused(true);
  const auto initial = effect.serialize_parameters();
  for (size_t i : {Effect::OCTET_PRESET_INDEX, Effect::OCTET_4D_PRESET_INDEX}) {
    auto selected = initial;
    selected.params = Effect::preset(i).params;
    HS_EXPECT_TRUE(effect.restore_parameters(selected));
    std::vector<Pixel> previous;
    for (int draw = 0; draw < 2; ++draw) {
      effect.draw_frame();
      effect.advance_display();
      std::vector<Pixel> current;
      capture_frame<96, 20>(effect, current);
      const size_t lit = count_lit_region<96, 20>(effect);
      const size_t changed =
          draw == 0 ? 0 : compare_buffers(previous, current).different_pixels;
      previous = std::move(current);
      HS_EXPECT_GT(lit, size_t{0});
      if (draw)
        HS_EXPECT_GT(changed, size_t{0});
      HS_REQUIRED_PARAM(unfinished, effect, "Unfinished Rays");
      HS_EXPECT_EQ(unfinished->get(), 0.0f);
    }
  }
}

inline void test_pattern_view_controls() {
  using Effect = HyperLatticeWhiteBox::Effect;
  reset_globals();
  Effect effect;
  effect.init();

  HS_REQUIRED_PARAM(pattern, effect, "Pattern");
  HS_REQUIRED_PARAM(view, effect, "View");
  HS_EXPECT_EQ(pattern->option_count, 3);
  HS_EXPECT_EQ(view->option_count, 2);
  HS_EXPECT_EQ(effect.updateParameter("View", 1), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(pattern->get(), 0);
  HS_EXPECT_EQ(view->get(), 1);
  HS_REQUIRED_PARAM(spin_4d, effect, "4D Spin");
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
  HS_EXPECT_STREQ(pattern->export_options[1], "Pattern::OCTET");
  const auto four_d = effect.serialize_parameters();
  HS_EXPECT_EQ(effect.updateParameter("Near Fade", .37f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("Pattern", 1), ParamSetResult::APPLIED);
  const auto octet4 = effect.serialize_parameters();
  HS_EXPECT_EQ(octet4.params.pattern, Effect::Pattern::OCTET);
  HS_EXPECT_EQ(octet4.params.mode, HL::LatticeMode::SLICE_4D);
  HS_EXPECT_EQ(octet4.params.cell_size,
               Effect::preset(Effect::OCTET_4D_PRESET_INDEX).params.cell_size);
  HS_EXPECT_EQ(octet4.params.near_fade, .37f);
  HS_REQUIRED_PARAM(octet_slice_spin_4d, effect, "4D Spin");
  HS_EXPECT_FALSE(octet_slice_spin_4d->readonly);
  HS_EXPECT_EQ(effect.updateParameter("View", 0), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(pattern->get(), 1);
  HS_EXPECT_EQ(view->get(), 0);
  HS_REQUIRED_PARAM(octet_shell_spin_4d, effect, "4D Spin");
  HS_EXPECT_TRUE(octet_shell_spin_4d->readonly);
  HS_EXPECT_EQ(effect.serialize_parameters().params.spin_4d, 0);
  HS_REQUIRED_PARAM(lattice_planes, effect, "Lattice Planes");
  HS_EXPECT_TRUE(lattice_planes->readonly);
  HS_REQUIRED_PARAM(wire_radius, effect, "Wire Radius");
  HS_EXPECT_FALSE(wire_radius->readonly);
  HS_REQUIRED_PARAM(aa_strength, effect, "AA Strength");
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
  HS_REQUIRED_PARAM(four_d_lattice_planes, effect, "Lattice Planes");
  HS_EXPECT_FALSE(four_d_lattice_planes->readonly);
}

/** @brief Shipping traced presets finish every ray at their supported grids. */
inline void test_traced_presets_finish_rays() {
  const auto check = []<int W, int H>(std::initializer_list<size_t> presets,
                                      int x0, int x1) {
    for (size_t preset : presets) {
      HS_CONTEXT("grid", W, H);
      HS_CONTEXT("preset", preset);
      reset_globals();
      HyperLattice<W, H> effect;
      effect.init();
      effect.setAnimationsPaused(true);
      effect.set_clip(0, H, x0, x1);
      HS_EXPECT_TRUE(effect.selectPreset(preset));
      for (int frame = 0; frame < 3; ++frame) {
        effect.draw_frame();
        effect.advance_display();
        HS_REQUIRED_PARAM(unfinished, effect, "Unfinished Rays");
        HS_EXPECT_EQ(unfinished->get(), 0.0f);
      }
    }
  };
  using Effect = HyperLatticeWhiteBox::Effect;
  check.template operator()<96, 20>({Effect::OCTET_PRESET_INDEX,
                                     Effect::OCTET_WIDE_PRESET_INDEX,
                                     Effect::OCTET_4D_PRESET_INDEX},
                                    0, 96);
  check.template operator()<288, 144>({Effect::SHELL_PRESET_INDEX,
                                       Effect::SHELL_CLOSE_PRESET_INDEX,
                                       Effect::SHELL_4D_PRESET_INDEX},
                                      140, 150);
}

/** @brief Restores cubic, octet and shell views and rejects pattern IDs 2–5. */
inline void test_regular_patterns() {
  using Effect = HyperLatticeWhiteBox::Effect;
  static_assert(static_cast<uint8_t>(Effect::Pattern::CUBIC_WIRE) == 0);
  static_assert(static_cast<uint8_t>(Effect::Pattern::OCTET) == 1);
  static_assert(static_cast<uint8_t>(Effect::Pattern::SHELLS) == 6);
  reset_globals();
  Effect effect;
  effect.init();
  effect.setAnimationsPaused(true);

  HS_REQUIRED_PARAM(pattern, effect, "Pattern");
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
    HS_REQUIRED_PARAM(view, effect, "View");
    HS_EXPECT_FALSE(view->readonly);
    const bool SHELLS = configuration.pattern == Effect::Pattern::SHELLS;
    HS_REQUIRED_PARAM(shell_radius, effect, "Shell Radius");
    HS_EXPECT_EQ(shell_radius->readonly, !SHELLS);
    HS_REQUIRED_PARAM(wire_radius, effect, "Wire Radius");
    HS_EXPECT_EQ(wire_radius->readonly, SHELLS);
    effect.draw_frame();
    effect.advance_display();
    const size_t lit = count_lit_region<96, 20>(effect);
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
      capture_frame<96, 20>(effect, actual);
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
      std::vector<Pixel> reference;
      capture_frame<96, 20>(effect, reference);
      const size_t differing =
          compare_buffers(reference, actual).different_pixels;
      HS_EXPECT_GT(differing, size_t{0});
    }
  }
  HS_EXPECT_EQ(effect.updateParameter("Pattern", 6), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.serialize_parameters().params.pattern,
               Effect::Pattern::SHELLS);
  HS_EXPECT_EQ(pattern->get(), 6);
  HS_EXPECT_STREQ(pattern->export_options[2], "Pattern::SHELLS");
}

/**
 * @brief The speed ceiling is accepted and a higher restored speed is
 * rejected.
 */
inline void test_speed_range() {
  using Effect = HyperLatticeWhiteBox::Effect;
  reset_globals();
  Effect effect;
  effect.init();
  HS_REQUIRED_PARAM(speed, effect, "Speed");
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
  HS_RUN(test_periodic_distance);
  HS_RUN(test_edge_metrics);
  HS_RUN(test_so4_rotation);
  HS_RUN(test_dimensional_rotation_wrap_is_continuous);
  HS_RUN(test_resolution_aware_wire_coverage);
  HS_RUN(test_near_field_fade);
  HS_RUN(test_far_shell_fade);
  HS_RUN(test_pause_does_not_stop_motion);
  HS_RUN(test_iridescent_controls_drive_palette);
  HS_RUN(test_depth_palette_mutates_slowly_while_paused);
  HS_RUN(test_depth_palette_keeps_cool_character);
  HS_RUN(test_next_plane_is_strict);
  HS_RUN(test_trace_layers_are_front_to_back);
  HS_RUN(test_layer_composite_reveals_background);
  HS_RUN(test_surface_origin_parallax);
  HS_RUN(test_hyperplane_event);
  HS_RUN(test_coincident_planes_form_one_layer);
  HS_RUN(test_render_signature);
  HS_RUN(test_specialized_slice_transition);
  HS_RUN(test_specialized_render_signature);
  HS_RUN(test_instantaneous_render_pins);
  HS_RUN(test_column_time_matches_time_evaluator);
  HS_RUN(test_column_time_segments_match_full_frame);
  HS_RUN(test_column_time_resets);
  HS_RUN(test_column_time_toggle);
  HS_RUN(test_column_time_parameter_motion);
  HS_RUN(test_frame_equality_ignores_padding);
  HS_RUN(test_presets_are_valid_and_distinct);
  HS_RUN(test_configuration_adoption_and_snapshots);
  HS_RUN(test_near_fade_control);
  HS_RUN(test_dimension_dropdown_and_mode_lerp);
  HS_RUN(test_single_shell);
  HS_RUN(test_family_segues);
  HS_RUN(test_regular_patterns);
  HS_RUN(test_octet_continuous_flight);
  HS_RUN(test_trace_geometry_mapping);
  HS_RUN(test_shared_trace_path_dispatch);
  HS_RUN(test_traced_preset_camera_setup);
  HS_RUN(test_traced_preset_control_state);
  HS_RUN(test_traced_preset_render);
  HS_RUN(test_traced_presets_finish_rays);
  HS_RUN(test_pattern_view_controls);
  HS_RUN(test_speed_range);
  return fixture.result();
}

} // namespace hyper_lattice_tests
} // namespace hs_test
