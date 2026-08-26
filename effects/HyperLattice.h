/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file HyperLattice.h
 * @brief Reflective flight through cubic and four-dimensional hypercubic
 *        lattices.
 */

#include <array>
#include <cmath>
#include <cstdint>
#include <string_view>
#include <tuple>

#include "core/color/effect_palette_recipes.h"
#include "core/color/layer_composite.h"
#include "core/control/choreography.h"
#include "core/engine/engine.h"
#include "core/math/4dmath.h"
#include "core/render/pullback.h"

namespace hs_test {
namespace hyper_lattice_tests {
struct HyperLatticeWhiteBox;
} // namespace hyper_lattice_tests
} // namespace hs_test

namespace HyperLatticeDetail {

constexpr int DIMENSIONS = VEC4_DIMENSIONS;
constexpr int MAX_SHELLS = 3;
constexpr float DIRECTION_EPSILON = 1.0e-4f;

enum class ReflectionMode : uint8_t { CHROME, RADIAL };
enum class ColorMode : uint8_t { DEPTH, AXIS };
enum class ShellCount : uint8_t { ONE, TWO, THREE };
enum class LatticeMode : uint8_t { THREE_D, DIMENSIONAL_RIFT, FOUR_D_SLICE };

constexpr float DIMENSIONAL_RIFT_MIX = 0.55f;

__attribute__((always_inline)) inline float periodic_distance(float value) {
  return fabsf(value - nearbyintf(value));
}

struct EdgeMetric {
  float distance_sq;
  uint8_t free_axis;
};

struct TransitionalMetrics {
  EdgeMetric cubic;
  EdgeMetric hyper;
};

__attribute__((always_inline)) inline float
periodic_distance_at(const Vec4 &ray_origin, const Vec4 &direction, int axis,
                     float distance) {
  return periodic_distance(ray_origin[axis] + distance * direction[axis]);
}

template <int AXIS0, int AXIS1>
__attribute__((always_inline)) EdgeMetric edge_metric_3d_axes(
    const Vec4 &ray_origin, const Vec4 &direction, float distance) {
  const float component0 =
      periodic_distance_at(ray_origin, direction, AXIS0, distance);
  const float component1 =
      periodic_distance_at(ray_origin, direction, AXIS1, distance);
  const float component0_sq = component0 * component0;
  const float component1_sq = component1 * component1;
  return {std::min(component0_sq, component1_sq),
          static_cast<uint8_t>(component1_sq > component0_sq ? AXIS1 : AXIS0)};
}

inline EdgeMetric edge_metric_3d_at(const Vec4 &ray_origin,
                                    const Vec4 &direction, int plane_axis,
                                    float distance) {
  switch (plane_axis) {
  case 0:
    return edge_metric_3d_axes<1, 2>(ray_origin, direction, distance);
  case 1:
    return edge_metric_3d_axes<0, 2>(ray_origin, direction, distance);
  default:
    return edge_metric_3d_axes<0, 1>(ray_origin, direction, distance);
  }
}

template <int AXIS0, int AXIS1, int AXIS2, bool NEED_AXIS = true>
__attribute__((always_inline)) bool
edge_metric_4d_axes_bounded(const Vec4 &ray_origin, const Vec4 &direction,
                            float distance, float limit, float limit_sq,
                            EdgeMetric &result) {
  const float component0 =
      periodic_distance_at(ray_origin, direction, AXIS0, distance);
  const float component1 =
      periodic_distance_at(ray_origin, direction, AXIS1, distance);
  const bool near0 = component0 < limit;
  const bool near1 = component1 < limit;
  if (!near0 && !near1)
    return false;
  const float component2 =
      periodic_distance_at(ray_origin, direction, AXIS2, distance);
  if (component2 >= limit && (!near0 || !near1))
    return false;
  const float component0_sq = component0 * component0;
  const float component1_sq = component1 * component1;
  const float component2_sq = component2 * component2;
  const float sum = component0_sq + component1_sq + component2_sq;
  if constexpr (!NEED_AXIS) {
    result = {
        sum - std::max(component0_sq, std::max(component1_sq, component2_sq)),
        0};
    return result.distance_sq < limit_sq;
  }
  float largest = component0_sq;
  uint8_t free_axis = AXIS0;
  if (component1_sq > largest) {
    largest = component1_sq;
    free_axis = AXIS1;
  }
  if (component2_sq > largest) {
    largest = component2_sq;
    free_axis = AXIS2;
  }
  result = {sum - largest, free_axis};
  return result.distance_sq < limit_sq;
}

template <bool NEED_AXIS = true>
inline bool edge_metric_4d_at_bounded(const Vec4 &ray_origin,
                                      const Vec4 &direction, int plane_axis,
                                      float distance, float limit,
                                      float limit_sq, EdgeMetric &result) {
  switch (plane_axis) {
  case 0:
    return edge_metric_4d_axes_bounded<1, 2, 3, NEED_AXIS>(
        ray_origin, direction, distance, limit, limit_sq, result);
  case 1:
    return edge_metric_4d_axes_bounded<0, 2, 3, NEED_AXIS>(
        ray_origin, direction, distance, limit, limit_sq, result);
  case 2:
    return edge_metric_4d_axes_bounded<0, 1, 3, NEED_AXIS>(
        ray_origin, direction, distance, limit, limit_sq, result);
  default:
    return edge_metric_4d_axes_bounded<0, 1, 2, NEED_AXIS>(
        ray_origin, direction, distance, limit, limit_sq, result);
  }
}

template <int AXIS0, int AXIS1>
__attribute__((always_inline)) TransitionalMetrics transitional_metrics_axes(
    const Vec4 &ray_origin, const Vec4 &direction, float distance) {
  const float component0 =
      periodic_distance_at(ray_origin, direction, AXIS0, distance);
  const float component1 =
      periodic_distance_at(ray_origin, direction, AXIS1, distance);
  const float component2 =
      periodic_distance_at(ray_origin, direction, 3, distance);
  const float component0_sq = component0 * component0;
  const float component1_sq = component1 * component1;
  const float component2_sq = component2 * component2;
  const float sum_4d = component0_sq + component1_sq + component2_sq;
  float farthest_4d = component0_sq;
  uint8_t free_axis_4d = AXIS0;
  if (component1_sq > farthest_4d) {
    farthest_4d = component1_sq;
    free_axis_4d = AXIS1;
  }
  if (component2_sq > farthest_4d) {
    farthest_4d = component2_sq;
    free_axis_4d = 3;
  }
  return {{std::min(component0_sq, component1_sq),
           static_cast<uint8_t>(component1_sq > component0_sq ? AXIS1 : AXIS0)},
          {sum_4d - farthest_4d, free_axis_4d}};
}

inline TransitionalMetrics transitional_metrics_at(const Vec4 &ray_origin,
                                                   const Vec4 &direction,
                                                   int plane_axis,
                                                   float distance) {
  switch (plane_axis) {
  case 0:
    return transitional_metrics_axes<1, 2>(ray_origin, direction, distance);
  case 1:
    return transitional_metrics_axes<0, 2>(ray_origin, direction, distance);
  default:
    return transitional_metrics_axes<0, 1>(ray_origin, direction, distance);
  }
}

struct Params {
  LatticeMode mode = LatticeMode::THREE_D;
  float sphere_radius = 0.4f;
  float wire_radius = 0.055f;
  float softness = 0.012f;
  float far_cells = 7.0f;
  float aa_strength = 1.0f;
  float speed = 0.018f;
  float spin_3d = 0.0024f;
  float spin_4d = 0.0f;
  float chrome_warp = 0.65f;
  ReflectionMode reflection = ReflectionMode::CHROME;
  ColorMode color = ColorMode::DEPTH;
  ShellCount shells = ShellCount::TWO;

  void lerp(const Params &start, const Params &target, float amount) {
    mode = amount < 0.5f ? start.mode : target.mode;
    sphere_radius = hs::lerp(start.sphere_radius, target.sphere_radius, amount);
    wire_radius = hs::lerp(start.wire_radius, target.wire_radius, amount);
    softness = hs::lerp(start.softness, target.softness, amount);
    far_cells = hs::lerp(start.far_cells, target.far_cells, amount);
    aa_strength = hs::lerp(start.aa_strength, target.aa_strength, amount);
    speed = hs::lerp(start.speed, target.speed, amount);
    spin_3d = hs::lerp(start.spin_3d, target.spin_3d, amount);
    spin_4d = hs::lerp(start.spin_4d, target.spin_4d, amount);
    chrome_warp = hs::lerp(start.chrome_warp, target.chrome_warp, amount);
    reflection = amount < 0.5f ? start.reflection : target.reflection;
    color = amount < 0.5f ? start.color : target.color;
    shells = amount < 0.5f ? start.shells : target.shells;
  }
};

// lerp() and valid_params() enumerate every Params field by hand; a field added
// without touching both pins silently at the start value across every
// transition. Every enum has a fixed uint8_t base, so the size holds under ARM
// -fshort-enums.
static_assert(sizeof(Params) == 44,
              "HyperLattice::Params field set changed — update lerp() and "
              "valid_params() to match");

struct FrameState {
  Params params;
  Vec4 origin;
  std::array<float, 6> rotation_phase;
  float pixel_half_angle;
  const BakedPalette *depth_palette;
  const BakedPalette *axis_palette;
};

struct Binding {
  using FrameState = HyperLatticeDetail::FrameState;
  using Instrumentation = Pullback::NoInstrumentation;
};

struct PreparedTrace {
  Params params;
  Vec4 origin;
  Mat4 world_to_lattice;
  float inv_far;
  float aa_scale;
  float outer_radius_base;
  float near_start;
  float near_inv_span;
  float dimension_mix;
  LatticeMode mode;
};

template <int W, int H> constexpr float pixel_half_angle() {
  constexpr float HORIZONTAL = TWO_PI_F / static_cast<float>(W);
  constexpr float VERTICAL = PI_F / static_cast<float>(H + hs::H_OFFSET - 1);
  return 0.5f * (HORIZONTAL > VERTICAL ? HORIZONTAL : VERTICAL);
}

__attribute__((always_inline)) inline float
lattice_ramp(float edge0, float edge1, float value) {
  return cubic_kernel((value - edge0) / (edge1 - edge0));
}

inline float wire_coverage(float metric_sq, float radius, float half_width) {
  const float signed_distance = sqrtf(metric_sq) - radius;
  return 1.0f - lattice_ramp(-half_width, half_width, signed_distance);
}

__attribute__((always_inline)) inline float
fast_wire_coverage(float metric_sq, float radius, float half_width) {
  const float signed_distance = sqrtf(metric_sq) - radius;
  const float ramp_position =
      0.5f + 0.5f * signed_distance * fast_reciprocal(half_width);
  return 1.0f - cubic_kernel(ramp_position);
}

inline float near_field_coverage(float distance, float near_start,
                                 float near_inv_span) {
  return cubic_kernel((distance - near_start) * near_inv_span);
}

inline float shell_horizon_coverage(uint8_t shell, uint8_t shell_count,
                                    float distance, float magnitude) {
  if (shell + 1 < shell_count)
    return 1.0f;
  return 1.0f - lattice_ramp(static_cast<float>(shell_count - 1),
                             static_cast<float>(shell_count),
                             distance * magnitude);
}

inline float next_plane_offset(float origin, bool positive) {
  const float fraction = wrap_t(origin);
  if (fraction == 0.0f)
    return 1.0f;
  return positive ? 1.0f - fraction : fraction;
}

inline float dimension_mix(LatticeMode mode) {
  if (mode == LatticeMode::THREE_D)
    return 0.0f;
  if (mode == LatticeMode::DIMENSIONAL_RIFT)
    return DIMENSIONAL_RIFT_MIX;
  return 1.0f;
}

inline PreparedTrace prepare_trace(const FrameState &frame) {
  PreparedTrace prepared;
  prepared.params = frame.params;
  prepared.origin = frame.origin;
  prepared.world_to_lattice = Mat4::identity();
  rotate_plane(prepared.world_to_lattice, 0, 1, frame.rotation_phase[0]);
  rotate_plane(prepared.world_to_lattice, 0, 2, frame.rotation_phase[1]);
  rotate_plane(prepared.world_to_lattice, 1, 2, frame.rotation_phase[2]);
  prepared.mode = frame.params.mode;
  prepared.dimension_mix = dimension_mix(frame.params.mode);
  rotate_plane(prepared.world_to_lattice, 0, 3,
               prepared.dimension_mix * frame.rotation_phase[3]);
  rotate_plane(prepared.world_to_lattice, 1, 3,
               prepared.dimension_mix * frame.rotation_phase[4]);
  rotate_plane(prepared.world_to_lattice, 2, 3,
               prepared.dimension_mix * frame.rotation_phase[5]);
  prepared.inv_far = 1.0f / frame.params.far_cells;
  prepared.aa_scale = frame.params.aa_strength * frame.pixel_half_angle;
  prepared.outer_radius_base = frame.params.wire_radius + frame.params.softness;
  prepared.near_start = 1.5f * frame.params.wire_radius;
  const float near_end = 4.0f * frame.params.wire_radius;
  prepared.near_inv_span = 1.0f / (near_end - prepared.near_start);
  return prepared;
}

inline Vec4 reflected_direction(const Vector &normal, ReflectionMode mode,
                                float chrome_warp) {
  if (mode == ReflectionMode::RADIAL)
    return {{normal.x, normal.y, normal.z, 0.0f}};
  const float twice_x = 2.0f * normal.x;
  const Vec4 reflected{{twice_x * normal.x - 1.0f, twice_x * normal.y,
                        twice_x * normal.z, 0.0f}};
  const float blend = 0.4f * chrome_warp;
  Vec4 direction{{hs::lerp(normal.x, reflected[0], blend),
                  hs::lerp(normal.y, reflected[1], blend),
                  hs::lerp(normal.z, reflected[2], blend), 0.0f}};
  const float inverse_length =
      fast_rsqrt(direction[0] * direction[0] + direction[1] * direction[1] +
                 direction[2] * direction[2]);
  for (int axis = 0; axis < 3; ++axis)
    direction[axis] *= inverse_length;
  return direction;
}

struct TraceHit {
  float coverage = 0.0f;
  float distance = 0.0f;
  uint8_t free_axis = 0;
};

template <bool SLICE_4D, bool SPECIALIZED_SLICE = false>
__attribute__((always_inline)) inline TraceHit
trace_plane(const Vec4 &ray_origin, const Vec4 &direction, int plane_axis,
            float distance, float plane_step, const PreparedTrace &prepared) {
  HS_PROFILE_DEEP(hl_plane_eval);
  if (distance <= prepared.near_start)
    return {0.0f, distance, 0};

  const float coverage_outer_radius =
      prepared.outer_radius_base + prepared.aa_scale * distance * plane_step;
  const float coverage_half_width =
      coverage_outer_radius - prepared.params.wire_radius;
  const float depth = distance * prepared.inv_far;
  const float fog = std::max(0.0f, 1.0f - depth);
  const float outer_radius = coverage_outer_radius;
  const float outer_radius_sq = outer_radius * outer_radius;
  float metric_sq;
  uint8_t free_axis;
  float dimensional_coverage = 1.0f;
  if constexpr (SLICE_4D) {
    EdgeMetric metric_4d;
    if (!edge_metric_4d_at_bounded<!SPECIALIZED_SLICE>(
            ray_origin, direction, plane_axis, distance, outer_radius,
            outer_radius_sq, metric_4d))
      return {0.0f, distance, 0};
    metric_sq = metric_4d.distance_sq;
    free_axis = metric_4d.free_axis;
  } else if (prepared.mode == LatticeMode::THREE_D) {
    const EdgeMetric metric_3d =
        edge_metric_3d_at(ray_origin, direction, plane_axis, distance);
    metric_sq = metric_3d.distance_sq;
    free_axis = metric_3d.free_axis;
    if (metric_sq >= outer_radius_sq)
      return {0.0f, distance, free_axis};
  } else if (plane_axis < 3 && prepared.mode == LatticeMode::DIMENSIONAL_RIFT) {
    const TransitionalMetrics metrics =
        transitional_metrics_at(ray_origin, direction, plane_axis, distance);
    metric_sq = hs::lerp(metrics.cubic.distance_sq, metrics.hyper.distance_sq,
                         prepared.dimension_mix);
    free_axis = prepared.dimension_mix < 0.5f ? metrics.cubic.free_axis
                                              : metrics.hyper.free_axis;
    if (metric_sq >= outer_radius_sq)
      return {0.0f, distance, free_axis};
  } else {
    EdgeMetric metric_4d;
    if (!edge_metric_4d_at_bounded(ray_origin, direction, plane_axis, distance,
                                   outer_radius, outer_radius_sq, metric_4d))
      return {0.0f, distance, 0};
    metric_sq = metric_4d.distance_sq;
    free_axis = metric_4d.free_axis;
    if (plane_axis == 3)
      dimensional_coverage = prepared.dimension_mix;
  }

  const float edge =
      SPECIALIZED_SLICE
          ? fast_wire_coverage(metric_sq, prepared.params.wire_radius,
                               coverage_half_width)
          : wire_coverage(metric_sq, prepared.params.wire_radius,
                          coverage_half_width);
  const float near = near_field_coverage(distance, prepared.near_start,
                                         prepared.near_inv_span);
  return {edge * fog * fog * near * dimensional_coverage, distance, free_axis};
}

struct TraceCursor {
  float distance;
  float step;
  float magnitude;
  uint8_t shell;
  bool active;
};

template <bool SLICE_4D = false, bool SPECIALIZED_SLICE = false,
          typename ConsumeFn>
__attribute__((always_inline)) inline void
trace_layers_mode(const Vector &normal, const PreparedTrace &prepared,
                  ConsumeFn consume) {
  HS_PROFILE_DEEP(hl_trace_layers);
  constexpr float GROUP_EPSILON = 1.0e-4f;
  Vec4 ray_origin = prepared.origin;
  if (prepared.params.sphere_radius != 0.0f) {
    const Vec4 surface_normal =
        prepared.world_to_lattice.apply({{normal.x, normal.y, normal.z, 0.0f}});
    for (int axis = 0; axis < DIMENSIONS; ++axis)
      ray_origin[axis] += prepared.params.sphere_radius * surface_normal[axis];
  }
  const Vec4 direction = prepared.world_to_lattice.apply(reflected_direction(
      normal, prepared.params.reflection, prepared.params.chrome_warp));
  const uint8_t shell_count =
      SPECIALIZED_SLICE ? 2 : static_cast<uint8_t>(prepared.params.shells) + 1;
  TraceCursor cursors[DIMENSIONS];
  for (int axis = 0; axis < DIMENSIONS; ++axis) {
    TraceCursor &cursor = cursors[axis];
    cursor.active = false;
    if constexpr (!SLICE_4D)
      if (axis == 3 && prepared.mode == LatticeMode::THREE_D)
        continue;
    const float component = direction[axis];
    const float magnitude = fabsf(component);
    if (magnitude < DIRECTION_EPSILON)
      continue;
    cursor.step = 1.0f / magnitude;
    cursor.magnitude = magnitude;
    cursor.distance =
        next_plane_offset(ray_origin[axis], component > 0.0f) * cursor.step;
    cursor.shell = 0;
    cursor.active = cursor.distance < prepared.params.far_cells;
  }

  for (int event = 0; event < DIMENSIONS * MAX_SHELLS; ++event) {
    HS_PROFILE_DEEP(hl_event_step);
    float nearest = prepared.params.far_cells;
    float second_nearest = prepared.params.far_cells;
    int nearest_axis = -1;
    const auto consider_cursor = [&](int axis) __attribute__((always_inline)) {
      const TraceCursor &cursor = cursors[axis];
      if (!cursor.active)
        return;
      if (cursor.distance < nearest) {
        second_nearest = nearest;
        nearest = cursor.distance;
        nearest_axis = axis;
      } else if (cursor.distance < second_nearest) {
        second_nearest = cursor.distance;
      }
    };
    consider_cursor(0);
    consider_cursor(1);
    consider_cursor(2);
    consider_cursor(3);
    if (nearest >= prepared.params.far_cells)
      break;

    const float tolerance = GROUP_EPSILON * std::max(1.0f, nearest);
    uint8_t pending = static_cast<uint8_t>(1u << nearest_axis);
    if (second_nearest - nearest <= tolerance) {
      pending = 0;
      for (int axis = 0; axis < DIMENSIONS; ++axis) {
        const TraceCursor &cursor = cursors[axis];
        if (cursor.active && cursor.distance - nearest <= tolerance)
          pending |= static_cast<uint8_t>(1u << axis);
      }
    }

    TraceHit layer;
    do {
      const int axis = __builtin_ctz(static_cast<unsigned>(pending));
      pending &= static_cast<uint8_t>(pending - 1);
      TraceCursor &cursor = cursors[axis];
      TraceHit candidate = trace_plane<SLICE_4D, SPECIALIZED_SLICE>(
          ray_origin, direction, axis, cursor.distance, cursor.step, prepared);
      if (candidate.coverage > 0.0f)
        candidate.coverage *= shell_horizon_coverage(
            cursor.shell, shell_count, cursor.distance, cursor.magnitude);
      if (candidate.coverage > layer.coverage)
        layer = candidate;
      ++cursor.shell;
      cursor.distance += cursor.step;
      cursor.active = cursor.shell < shell_count &&
                      cursor.distance < prepared.params.far_cells;
    } while (pending != 0);
    if (layer.coverage > 0.0f && !consume(layer))
      return;
  }
}

template <typename ConsumeFn>
inline void trace_layers(const Vector &normal, const PreparedTrace &prepared,
                         ConsumeFn consume) {
  trace_layers_mode(normal, prepared, consume);
}

inline TraceHit trace(const Vector &normal, const PreparedTrace &prepared) {
  TraceHit nearest;
  trace_layers(normal, prepared, [&](const TraceHit &hit) {
    nearest = hit;
    return false;
  });
  return nearest;
}

template <bool SLICE_4D = false, bool SPECIALIZED_SLICE = false>
__attribute__((always_inline)) inline Color4
shade_mode(const Pullback::SphereSample &input, const FrameState &frame,
           const PreparedTrace &prepared) {
  HS_PROFILE_DEEP(hl_shade);
  LayerComposite composite;
  const BakedPalette &palette =
      *(SPECIALIZED_SLICE || prepared.params.color == ColorMode::DEPTH
            ? frame.depth_palette
            : frame.axis_palette);
  trace_layers_mode<SLICE_4D, SPECIALIZED_SLICE>(
      input.dir, prepared,
      [&](const TraceHit &hit) __attribute__((always_inline)) {
        HS_PROFILE_DEEP(hl_layer_composite);
        const float depth = hit.distance * prepared.inv_far;
        const float value =
            SPECIALIZED_SLICE || prepared.params.color == ColorMode::DEPTH
                ? 1.0f - depth
                : (static_cast<float>(hit.free_axis) + 0.75f * depth) / 4.0f;
        Pixel color = palette.get_color_unit(value);
        color = color * (0.45f + 0.55f * (1.0f - depth));
        composite.add(color, hit.coverage);
        return composite.remaining > MIN_ENCODABLE_ALPHA;
      });
  return composite.finish();
}

inline Color4 shade(const Pullback::SphereSample &input,
                    const FrameState &frame, const PreparedTrace &prepared) {
  return shade_mode(input, frame, prepared);
}

template <bool SLICE_4D = false, bool SPECIALIZED_SLICE = false>
struct ModeShadeStage
    : Pullback::Stage::Contract<ModeShadeStage<SLICE_4D, SPECIALIZED_SLICE>,
                                Pullback::SphereSample, Color4> {
  using Policies = std::tuple<>;

  template <typename PipelineBinding>
  static PreparedTrace
  prepare(const typename PipelineBinding::FrameState &frame) {
    return prepare_trace(frame);
  }

  template <typename PipelineBinding>
  __attribute__((always_inline)) static Color4
  run(const Pullback::SphereSample &input,
      const typename PipelineBinding::FrameState &frame,
      const PreparedTrace &prepared) {
    return shade_mode<SLICE_4D, SPECIALIZED_SLICE>(input, frame, prepared);
  }
};

using RenderPipeline = Pullback::Pipeline<Binding, ModeShadeStage<>>;
using SpecializedRenderPipeline =
    Pullback::Pipeline<Binding, ModeShadeStage<true, true>>;

} // namespace HyperLatticeDetail

template <int W, int H>
class HyperLattice : public ChoreographedEffect<HyperLattice<W, H>,
                                                HyperLatticeDetail::Params> {
  using Choreography =
      ChoreographedEffect<HyperLattice<W, H>, HyperLatticeDetail::Params>;
  friend Choreography;

public:
  using Params = HyperLatticeDetail::Params;
  using LatticeMode = HyperLatticeDetail::LatticeMode;
  using ReflectionMode = HyperLatticeDetail::ReflectionMode;
  using ColorMode = HyperLatticeDetail::ColorMode;
  using ShellCount = HyperLatticeDetail::ShellCount;

  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "cubic-flight", "hypercube-flight"};
  static constexpr Segue::Preset::Lerp PRESET_SEGUE{240, ease_in_out_sin,
                                                    /*pausable=*/true};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 320;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 6;

  static constexpr Params preset_params(size_t index) {
    Params value;
    switch (index) {
    case 0:
      value.mode = LatticeMode::THREE_D;
      value.sphere_radius = 0.4f;
      value.wire_radius = 0.055f;
      value.softness = 0.08f;
      value.far_cells = 4.198f;
      value.aa_strength = 1.0f;
      value.speed = 0.05f;
      value.spin_3d = 0.015f;
      value.spin_4d = 0.0f;
      value.chrome_warp = 0.65f;
      value.reflection = ReflectionMode::RADIAL;
      value.color = ColorMode::DEPTH;
      value.shells = ShellCount::TWO;
      break;
    case 1:
      value.mode = LatticeMode::FOUR_D_SLICE;
      value.sphere_radius = 0.0f;
      value.wire_radius = 0.03546f;
      value.softness = 0.029612f;
      value.far_cells = 8.0f;
      value.aa_strength = 1.0f;
      value.speed = 0.03f;
      value.spin_3d = 0.01089f;
      value.spin_4d = 0.015f;
      value.chrome_warp = 0.65f;
      value.reflection = ReflectionMode::CHROME;
      value.color = ColorMode::DEPTH;
      value.shells = ShellCount::TWO;
      break;
    default:
      break;
    }
    return value;
  }

  static constexpr bool valid_params(const Params &value) {
    return static_cast<uint8_t>(value.mode) <=
               static_cast<uint8_t>(LatticeMode::FOUR_D_SLICE) &&
           value.sphere_radius >= SPHERE_RADIUS_MIN &&
           value.sphere_radius <= SPHERE_RADIUS_MAX &&
           value.wire_radius >= WIRE_RADIUS_MIN &&
           value.wire_radius <= WIRE_RADIUS_MAX &&
           value.softness >= SOFTNESS_MIN && value.softness <= SOFTNESS_MAX &&
           value.far_cells >= FAR_CELLS_MIN &&
           value.far_cells <= FAR_CELLS_MAX &&
           value.aa_strength >= AA_STRENGTH_MIN &&
           value.aa_strength <= AA_STRENGTH_MAX && value.speed >= SPEED_MIN &&
           value.speed <= SPEED_MAX && value.spin_3d >= SPIN_3D_MIN &&
           value.spin_3d <= SPIN_3D_MAX && value.spin_4d >= SPIN_4D_MIN &&
           value.spin_4d <= SPIN_4D_MAX &&
           value.chrome_warp >= CHROME_WARP_MIN &&
           value.chrome_warp <= CHROME_WARP_MAX &&
           static_cast<uint8_t>(value.reflection) <=
               static_cast<uint8_t>(ReflectionMode::RADIAL) &&
           static_cast<uint8_t>(value.color) <=
               static_cast<uint8_t>(ColorMode::AXIS) &&
           static_cast<uint8_t>(value.shells) <=
               static_cast<uint8_t>(ShellCount::THREE);
  }

  /**
   * @brief Whether a parameter set is the shape SpecializedRenderPipeline bakes
   *        in.
   * @param value Parameters to test.
   * @return true when the parameters match the specialized trace assumptions.
   * @details The shape is preset 1's; the assert in draw_frame() ties the two,
   *          so retuning that preset cannot leave the gate behind.
   */
  static constexpr bool uses_specialized_slice(const Params &value) {
    return value.mode == LatticeMode::FOUR_D_SLICE &&
           value.reflection == ReflectionMode::CHROME &&
           value.color == ColorMode::DEPTH && value.shells == ShellCount::TWO;
  }

  HS_COLD_MEMBER HyperLattice() : Choreography(W, H, {.strobe = true}) {}

  void init() override {
    begin_choreography();
    register_animated_param("Dimension", &params.mode, MODE_OPTIONS,
                            MODE_EXPORT_OPTIONS, std::size(MODE_OPTIONS));
    register_animated_param("Sphere Radius", &params.sphere_radius,
                            SPHERE_RADIUS_MIN, SPHERE_RADIUS_MAX);
    register_animated_param("Wire Radius", &params.wire_radius, WIRE_RADIUS_MIN,
                            WIRE_RADIUS_MAX);
    register_animated_param("Softness", &params.softness, SOFTNESS_MIN,
                            SOFTNESS_MAX);
    register_animated_param("Far Cells", &params.far_cells, FAR_CELLS_MIN,
                            FAR_CELLS_MAX);
    register_animated_param("AA Strength", &params.aa_strength, AA_STRENGTH_MIN,
                            AA_STRENGTH_MAX);
    register_animated_param("Speed", &params.speed, SPEED_MIN, SPEED_MAX);
    register_animated_param("3D Spin", &params.spin_3d, SPIN_3D_MIN,
                            SPIN_3D_MAX);
    register_animated_param("4D Spin", &params.spin_4d, SPIN_4D_MIN,
                            SPIN_4D_MAX);
    register_animated_param("Chrome Warp", &params.chrome_warp, CHROME_WARP_MIN,
                            CHROME_WARP_MAX);
    register_animated_param("Reflection", &params.reflection,
                            REFLECTION_OPTIONS, REFLECTION_EXPORT_OPTIONS,
                            std::size(REFLECTION_OPTIONS));
    register_animated_param("Color", &params.color, COLOR_OPTIONS,
                            COLOR_EXPORT_OPTIONS, std::size(COLOR_OPTIONS));
    register_animated_param("Shells", &params.shells, SHELL_OPTIONS,
                            SHELL_EXPORT_OPTIONS, std::size(SHELL_OPTIONS));
    depth_palette.init_generated(persistent_arena, next_depth_palette, nullptr,
                                 0, PALETTE_FADE_FRAMES, ease_in_out_sin);
    const GenerativePalette fixed_axis_palette{
        EffectPaletteRecipes::raymarch()};
    axis_palette.bake(persistent_arena, fixed_axis_palette);
  }

  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(hl_timeline_step);
      timeline.step(canvas);
    }
    step_choreography();
    advance_state();
    depth_palette.step();
    const HyperLatticeDetail::FrameState context{
        params,
        origin,
        rotation_phase,
        HyperLatticeDetail::pixel_half_angle<W, H>(),
        &depth_palette.palette(),
        &axis_palette};
    const auto frame = HyperLatticeDetail::RenderPipeline::prepare(context);
    {
      HS_PROFILE(hl_shader_draw);
      static_assert(uses_specialized_slice(preset_params(1)),
                    "preset 1 no longer selects the specialized slice trace");
      if (uses_specialized_slice(frame.ctx.params)) {
        Scan::Shader::draw_cached<W, H, 1>(
            canvas, [&frame](const Vector &view) {
              return HyperLatticeDetail::SpecializedRenderPipeline::evaluate(
                  view, frame.ctx, frame.prepared);
            });
      } else {
        Scan::Shader::draw_cached<W, H, 1>(
            canvas, [&frame](const Vector &view) HS_HOT_FLASH_MEMBER {
              return HyperLatticeDetail::RenderPipeline::evaluate(
                  view, frame.ctx, frame.prepared);
            });
      }
    }
  }

#if HS_ENABLE_EFFECT_CONTROL_API
  void profile_select_preset(size_t index) {
    if (index >= PRESET_IDS.size() || !this->selectPreset(index)) {
      hs::log("Profile preset selection failed: %u/%u",
              static_cast<unsigned>(index),
              static_cast<unsigned>(PRESET_IDS.size()));
      return;
    }
    hs::log("Profile preset: %u/%u", static_cast<unsigned>(index),
            static_cast<unsigned>(PRESET_IDS.size()));
  }
#endif

private:
  using Choreography::begin_choreography;
  using Choreography::params;
  using Choreography::step_choreography;
  using Choreography::register_animated_param;
  using Choreography::timeline;
  using Choreography::transition;

  void blend_params(float progress) {
    params.lerp(transition.from, transition.to, progress);
  }

  void advance_state() {
    static constexpr float VELOCITY[HyperLatticeDetail::DIMENSIONS] = {
        0.4815434f, 0.2993373f, 0.4034555f, 0.7223151f};
    static constexpr float RATE[6] = {1.0f, 0.731f, 0.517f,
                                      1.0f, 0.707f, 0.419f};
    for (int axis = 0; axis < HyperLatticeDetail::DIMENSIONS; ++axis)
      origin[axis] = wrap_t(origin[axis] + params.speed * VELOCITY[axis]);
    for (int plane = 0; plane < 3; ++plane)
      rotation_phase[plane] =
          wrap(rotation_phase[plane] + params.spin_3d * RATE[plane], TWO_PI_F);
    for (int plane = 3; plane < 6; ++plane)
      rotation_phase[plane] =
          wrap(rotation_phase[plane] + params.spin_4d * RATE[plane], TWO_PI_F);
  }

  static void next_depth_palette(void *, uint32_t sequence,
                                 GenerativePalette &out) {
    static constexpr uint32_t BASE_HUE = 219;
    static constexpr uint32_t HUE_STEP = 159;
    out = GenerativePalette{EffectPaletteRecipes::raymarch_at(
        PaletteRecipes::hue_turns(BASE_HUE + sequence * HUE_STEP))};
  }

  static constexpr float SPHERE_RADIUS_MIN = 0.0f, SPHERE_RADIUS_MAX = 1.5f;
  static constexpr float WIRE_RADIUS_MIN = 0.015f, WIRE_RADIUS_MAX = 0.18f;
  static constexpr float SOFTNESS_MIN = 0.002f, SOFTNESS_MAX = 0.08f;
  static constexpr float FAR_CELLS_MIN = 2.0f, FAR_CELLS_MAX = 16.0f;
  static constexpr float AA_STRENGTH_MIN = 0.0f, AA_STRENGTH_MAX = 2.0f;
  static constexpr float SPEED_MIN = 0.0f, SPEED_MAX = 0.05f;
  static constexpr float SPIN_3D_MIN = 0.0f, SPIN_3D_MAX = 0.015f;
  static constexpr float SPIN_4D_MIN = 0.0f, SPIN_4D_MAX = 0.015f;
  static constexpr float CHROME_WARP_MIN = 0.0f, CHROME_WARP_MAX = 1.0f;

  static constexpr int PALETTE_FADE_FRAMES = 960;

  static constexpr const char *MODE_OPTIONS[] = {"3D", "Dimensional Rift",
                                                 "4D Slice"};
  static constexpr const char *MODE_EXPORT_OPTIONS[] = {
      "LatticeMode::THREE_D", "LatticeMode::DIMENSIONAL_RIFT",
      "LatticeMode::FOUR_D_SLICE"};
  static constexpr const char *REFLECTION_OPTIONS[] = {"Mirror Ball",
                                                       "Embedded Sphere"};
  static constexpr const char *REFLECTION_EXPORT_OPTIONS[] = {
      "ReflectionMode::CHROME", "ReflectionMode::RADIAL"};
  static constexpr const char *COLOR_OPTIONS[] = {"Depth", "Axis"};
  static constexpr const char *COLOR_EXPORT_OPTIONS[] = {"ColorMode::DEPTH",
                                                         "ColorMode::AXIS"};
  static constexpr const char *SHELL_OPTIONS[] = {"1", "2", "3"};
  static constexpr const char *SHELL_EXPORT_OPTIONS[] = {
      "ShellCount::ONE", "ShellCount::TWO", "ShellCount::THREE"};

  Vec4 origin{{0.17f, 0.31f, 0.43f, 0.59f}};
  std::array<float, 6> rotation_phase{};
  PaletteCycler depth_palette;
  BakedPalette axis_palette;

  friend struct hs_test::hyper_lattice_tests::HyperLatticeWhiteBox;

  static constexpr size_t FOOTPRINT_BYTES =
      PaletteCycler::generated_arena_bytes() +
      BakedPalette::required_arena_bytes();
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "HyperLattice persistent footprint exceeds the default "
                "partition");
};

static_assert(
    [] {
      using Effect = HyperLattice<1, 1>;
      for (size_t index = 0; index < Effect::PRESET_IDS.size(); ++index)
        if (!Effect::valid_params(Effect::preset_params(index)))
          return false;
      return true;
    }(),
    "HyperLattice preset is outside a registered slider range");

#include "core/control/registry.h"
REGISTER_EFFECT(HyperLattice)
