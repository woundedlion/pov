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
#include "core/control/choreography.h"
#include "core/engine/engine.h"
#include "core/render/pullback.h"

namespace hs_test {
namespace hyper_lattice_tests {
struct HyperLatticeWhiteBox;
} // namespace hyper_lattice_tests
} // namespace hs_test

namespace HyperLatticeDetail {

constexpr int DIMENSIONS = 4;
constexpr int MAX_SHELLS = 3;
constexpr int MAX_PROJECTED_LINES = 32;
constexpr float DIRECTION_EPSILON = 1.0e-4f;

enum class ReflectionMode : uint8_t { CHROME, RADIAL };
enum class ColorMode : uint8_t { DEPTH, AXIS };
enum class ShellCount : uint8_t { ONE, TWO, THREE };
enum class LatticeMode : uint8_t {
  THREE_D,
  DIMENSIONAL_RIFT,
  FOUR_D_SLICE,
  FOUR_D_PROJECTED
};

constexpr float DIMENSIONAL_RIFT_MIX = 0.55f;

struct Vec4 {
  float v[DIMENSIONS];

  constexpr float &operator[](int index) { return v[index]; }
  constexpr float operator[](int index) const { return v[index]; }
};

struct Mat4 {
  float m[DIMENSIONS][DIMENSIONS]{};

  static constexpr Mat4 identity() {
    Mat4 result;
    for (int i = 0; i < DIMENSIONS; ++i)
      result.m[i][i] = 1.0f;
    return result;
  }

  __attribute__((always_inline)) Vec4 apply(const Vec4 &input) const {
    return {{m[0][0] * input[0] + m[0][1] * input[1] + m[0][2] * input[2] +
                 m[0][3] * input[3],
             m[1][0] * input[0] + m[1][1] * input[1] + m[1][2] * input[2] +
                 m[1][3] * input[3],
             m[2][0] * input[0] + m[2][1] * input[1] + m[2][2] * input[2] +
                 m[2][3] * input[3],
             m[3][0] * input[0] + m[3][1] * input[1] + m[3][2] * input[2] +
                 m[3][3] * input[3]}};
  }
};

HS_COLD static void rotate_plane(Mat4 &matrix, int a, int b, float angle) {
  const float c = fast_cosf(angle);
  const float s = fast_sinf(angle);
  for (int column = 0; column < DIMENSIONS; ++column) {
    const float av = matrix.m[a][column];
    const float bv = matrix.m[b][column];
    matrix.m[a][column] = c * av - s * bv;
    matrix.m[b][column] = s * av + c * bv;
  }
}

__attribute__((always_inline)) inline float periodic_distance(float value) {
  return fabsf(value - nearbyintf(value));
}

struct EdgeMetric {
  float distance_sq;
  uint8_t free_axis;
};

inline EdgeMetric edge_metric_3d(const Vec4 &point, int plane_axis) {
  float best_fixed = 1.0f;
  float farthest = -1.0f;
  uint8_t free_axis = 0;
  for (int axis = 0; axis < 3; ++axis) {
    if (axis == plane_axis)
      continue;
    const float distance = periodic_distance(point[axis]);
    const float distance_sq = distance * distance;
    best_fixed = std::min(best_fixed, distance_sq);
    if (distance_sq > farthest) {
      farthest = distance_sq;
      free_axis = static_cast<uint8_t>(axis);
    }
  }
  return {best_fixed, free_axis};
}

inline EdgeMetric edge_metric_4d(const Vec4 &point, int plane_axis) {
  float sum = 0.0f;
  float largest = -1.0f;
  uint8_t free_axis = 0;
  for (int axis = 0; axis < DIMENSIONS; ++axis) {
    if (axis == plane_axis)
      continue;
    const float distance = periodic_distance(point[axis]);
    const float distance_sq = distance * distance;
    sum += distance_sq;
    if (distance_sq > largest) {
      largest = distance_sq;
      free_axis = static_cast<uint8_t>(axis);
    }
  }
  return {sum - largest, free_axis};
}

struct TransitionalMetrics {
  EdgeMetric cubic;
  EdgeMetric hyper;
};

inline TransitionalMetrics transitional_metrics(const Vec4 &point,
                                                int plane_axis) {
  float best_fixed = 1.0f;
  float farthest_3d = -1.0f;
  float sum_4d = 0.0f;
  float farthest_4d = -1.0f;
  uint8_t free_axis_3d = 0;
  uint8_t free_axis_4d = 0;
  for (int axis = 0; axis < DIMENSIONS; ++axis) {
    if (axis == plane_axis)
      continue;
    const float distance = periodic_distance(point[axis]);
    const float distance_sq = distance * distance;
    sum_4d += distance_sq;
    if (distance_sq > farthest_4d) {
      farthest_4d = distance_sq;
      free_axis_4d = static_cast<uint8_t>(axis);
    }
    if (axis < 3) {
      best_fixed = std::min(best_fixed, distance_sq);
      if (distance_sq > farthest_3d) {
        farthest_3d = distance_sq;
        free_axis_3d = static_cast<uint8_t>(axis);
      }
    }
  }
  return {{best_fixed, free_axis_3d}, {sum_4d - farthest_4d, free_axis_4d}};
}

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

template <int AXIS0, int AXIS1, int AXIS2>
__attribute__((always_inline)) bool
edge_metric_4d_axes_bounded(const Vec4 &ray_origin, const Vec4 &direction,
                            float distance, float limit_sq,
                            EdgeMetric &result) {
  const float component0 =
      periodic_distance_at(ray_origin, direction, AXIS0, distance);
  const float component1 =
      periodic_distance_at(ray_origin, direction, AXIS1, distance);
  const float component0_sq = component0 * component0;
  const float component1_sq = component1 * component1;
  if (std::min(component0_sq, component1_sq) >= limit_sq)
    return false;
  const float component2 =
      periodic_distance_at(ray_origin, direction, AXIS2, distance);
  const float component2_sq = component2 * component2;
  const float sum = component0_sq + component1_sq + component2_sq;
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

inline bool edge_metric_4d_at_bounded(const Vec4 &ray_origin,
                                      const Vec4 &direction, int plane_axis,
                                      float distance, float limit_sq,
                                      EdgeMetric &result) {
  switch (plane_axis) {
  case 0:
    return edge_metric_4d_axes_bounded<1, 2, 3>(ray_origin, direction, distance,
                                                limit_sq, result);
  case 1:
    return edge_metric_4d_axes_bounded<0, 2, 3>(ray_origin, direction, distance,
                                                limit_sq, result);
  case 2:
    return edge_metric_4d_axes_bounded<0, 1, 3>(ray_origin, direction, distance,
                                                limit_sq, result);
  default:
    return edge_metric_4d_axes_bounded<0, 1, 2>(ray_origin, direction, distance,
                                                limit_sq, result);
  }
}

template <int AXIS0, int AXIS1, int AXIS2>
__attribute__((always_inline)) EdgeMetric edge_metric_4d_axes(
    const Vec4 &ray_origin, const Vec4 &direction, float distance) {
  const float component0 =
      periodic_distance_at(ray_origin, direction, AXIS0, distance);
  const float component1 =
      periodic_distance_at(ray_origin, direction, AXIS1, distance);
  const float component2 =
      periodic_distance_at(ray_origin, direction, AXIS2, distance);
  const float component0_sq = component0 * component0;
  const float component1_sq = component1 * component1;
  const float component2_sq = component2 * component2;
  const float sum = component0_sq + component1_sq + component2_sq;
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
  return {sum - largest, free_axis};
}

inline EdgeMetric edge_metric_4d_at(const Vec4 &ray_origin,
                                    const Vec4 &direction, int plane_axis,
                                    float distance) {
  switch (plane_axis) {
  case 0:
    return edge_metric_4d_axes<1, 2, 3>(ray_origin, direction, distance);
  case 1:
    return edge_metric_4d_axes<0, 2, 3>(ray_origin, direction, distance);
  case 2:
    return edge_metric_4d_axes<0, 1, 3>(ray_origin, direction, distance);
  default:
    return edge_metric_4d_axes<0, 1, 2>(ray_origin, direction, distance);
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

struct ProjectedLine {
  Vector anchor;
  uint8_t free_axis;
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
  std::array<Vector, DIMENSIONS> projected_axes;
  std::array<ProjectedLine, MAX_PROJECTED_LINES> projected_lines;
  uint8_t projected_line_count;
};

inline void prepare_projected_lines(PreparedTrace &prepared);

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

inline float projected_half_width(float distance, float plane_step,
                                  float aa_scale, float softness) {
  return softness + aa_scale * distance * plane_step;
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
  prepared.projected_line_count = 0;
  if (prepared.mode == LatticeMode::FOUR_D_PROJECTED)
    prepare_projected_lines(prepared);
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

__attribute__((always_inline)) inline TraceHit
trace_plane(const Vec4 &ray_origin, const Vec4 &direction, int plane_axis,
            float distance, float plane_step, const PreparedTrace &prepared) {
  HS_PROFILE_DEEP(hl_plane_eval);
  if (distance <= prepared.near_start)
    return {0.0f, distance, 0};

  const float outer_radius =
      prepared.outer_radius_base + prepared.aa_scale * distance * plane_step;
  const float outer_radius_sq = outer_radius * outer_radius;
  float metric_sq;
  uint8_t free_axis;
  float dimensional_coverage = 1.0f;
  if (prepared.mode == LatticeMode::THREE_D) {
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
                                   outer_radius_sq, metric_4d))
      return {0.0f, distance, 0};
    metric_sq = metric_4d.distance_sq;
    free_axis = metric_4d.free_axis;
    if (plane_axis == 3)
      dimensional_coverage = prepared.dimension_mix;
  }

  const float half_width = outer_radius - prepared.params.wire_radius;
  const float edge =
      wire_coverage(metric_sq, prepared.params.wire_radius, half_width);
  const float fog = std::max(0.0f, 1.0f - distance * prepared.inv_far);
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

struct ProjectedDistance {
  float distance_sq;
  float ray_distance;
};

inline Vector projected_lattice_axis(const PreparedTrace &prepared, int axis) {
  return Vector(prepared.world_to_lattice.m[axis][0],
                prepared.world_to_lattice.m[axis][1],
                prepared.world_to_lattice.m[axis][2]);
}

inline ProjectedDistance projected_line_distance(const Vector &ray_origin,
                                                 const Vector &direction,
                                                 const Vector &anchor,
                                                 const Vector &line_direction) {
  const Vector delta = anchor - ray_origin;
  const float line_length_sq = dot(line_direction, line_direction);
  if (line_length_sq < DIRECTION_EPSILON * DIRECTION_EPSILON) {
    const float ray_distance = dot(delta, direction);
    const Vector separation = delta - direction * ray_distance;
    return {dot(separation, separation), ray_distance};
  }
  const float coupling = dot(direction, line_direction);
  const float ray_projection = dot(direction, delta);
  const float line_projection = dot(line_direction, delta);
  const float denominator = line_length_sq - coupling * coupling;
  float line_position;
  float ray_distance;
  if (denominator > DIRECTION_EPSILON * DIRECTION_EPSILON) {
    line_position = (coupling * ray_projection - line_projection) / denominator;
    ray_distance = ray_projection + coupling * line_position;
  } else {
    ray_distance = ray_projection;
    line_position =
        (coupling * ray_distance - line_projection) / line_length_sq;
  }
  const Vector separation =
      delta + line_direction * line_position - direction * ray_distance;
  return {dot(separation, separation), ray_distance};
}

inline bool projected_lines_coincident(const ProjectedLine &left,
                                       const ProjectedLine &right,
                                       const Vector &line_direction) {
  if (left.free_axis != right.free_axis)
    return false;
  const Vector delta = left.anchor - right.anchor;
  const float line_length_sq = dot(line_direction, line_direction);
  Vector separation = delta;
  if (line_length_sq >= DIRECTION_EPSILON * DIRECTION_EPSILON)
    separation -=
        line_direction * (dot(delta, line_direction) / line_length_sq);
  return dot(separation, separation) < 1.0e-8f;
}

inline void prepare_projected_lines(PreparedTrace &prepared) {
  constexpr float HIDDEN_DEPTH_WEIGHT = 0.0625f;
  constexpr uint8_t QUOTAS[MAX_SHELLS] = {3, 6, 8};
  struct RankedLine {
    ProjectedLine line;
    float score;
    bool used;
  };

  const Vec4 hidden_direction =
      prepared.world_to_lattice.apply({{0.0f, 0.0f, 0.0f, 1.0f}});
  for (int axis = 0; axis < DIMENSIONS; ++axis)
    prepared.projected_axes[axis] = projected_lattice_axis(prepared, axis);
  const uint8_t quota = QUOTAS[static_cast<uint8_t>(prepared.params.shells)];
  for (uint8_t free_axis = 0; free_axis < DIMENSIONS; ++free_axis) {
    RankedLine ranked[27];
    for (int ordinal = 0; ordinal < 27; ++ordinal) {
      int digits = ordinal;
      Vector anchor(0.0f, 0.0f, 0.0f);
      float hidden_anchor = 0.0f;
      for (int axis = 0; axis < DIMENSIONS; ++axis) {
        if (axis == free_axis)
          continue;
        const int offset = digits % 3 - 1;
        digits /= 3;
        const float fixed =
            nearbyintf(prepared.origin[axis]) + static_cast<float>(offset);
        const float delta = fixed - prepared.origin[axis];
        anchor += prepared.projected_axes[axis] * delta;
        hidden_anchor += hidden_direction[axis] * delta;
      }
      const Vector line_direction = prepared.projected_axes[free_axis];
      const float line_length_sq = dot(line_direction, line_direction);
      float line_position = 0.0f;
      float projected_distance_sq = dot(anchor, anchor);
      if (line_length_sq >= DIRECTION_EPSILON * DIRECTION_EPSILON) {
        line_position = -dot(anchor, line_direction) / line_length_sq;
        const Vector separation = anchor + line_direction * line_position;
        projected_distance_sq = dot(separation, separation);
      }
      const float hidden_depth =
          hidden_anchor + hidden_direction[free_axis] * line_position;
      ranked[ordinal] = {{anchor, free_axis},
                         projected_distance_sq +
                             HIDDEN_DEPTH_WEIGHT * hidden_depth * hidden_depth,
                         false};
    }

    uint8_t selected = 0;
    while (selected < quota) {
      int best = -1;
      for (int ordinal = 0; ordinal < 27; ++ordinal)
        if (!ranked[ordinal].used &&
            (best < 0 || ranked[ordinal].score < ranked[best].score))
          best = ordinal;
      if (best < 0)
        break;
      ranked[best].used = true;
      bool duplicate = false;
      for (uint8_t index = 0; index < prepared.projected_line_count; ++index)
        duplicate |= projected_lines_coincident(
            ranked[best].line, prepared.projected_lines[index],
            prepared.projected_axes[free_axis]);
      if (duplicate)
        continue;
      prepared.projected_lines[prepared.projected_line_count++] =
          ranked[best].line;
      ++selected;
    }
  }
}

template <typename ConsumeFn>
__attribute__((always_inline)) inline void
trace_projected_layers(const Vector &normal, const PreparedTrace &prepared,
                       ConsumeFn consume) {
  const Vec4 reflected = reflected_direction(normal, prepared.params.reflection,
                                             prepared.params.chrome_warp);
  const Vector ray_origin = normal * prepared.params.sphere_radius;
  const Vector direction(reflected[0], reflected[1], reflected[2]);
  TraceHit hits[MAX_PROJECTED_LINES];
  uint8_t hit_count = 0;
  for (uint8_t index = 0; index < prepared.projected_line_count; ++index) {
    const ProjectedLine &line = prepared.projected_lines[index];
    const ProjectedDistance projected =
        projected_line_distance(ray_origin, direction, line.anchor,
                                prepared.projected_axes[line.free_axis]);
    if (projected.ray_distance <= prepared.near_start ||
        projected.ray_distance >= prepared.params.far_cells)
      continue;
    const float half_width =
        projected_half_width(projected.ray_distance, 1.0f, prepared.aa_scale,
                             prepared.params.softness);
    const float outer_radius = prepared.params.wire_radius + half_width;
    if (projected.distance_sq >= outer_radius * outer_radius)
      continue;
    const float edge = wire_coverage(projected.distance_sq,
                                     prepared.params.wire_radius, half_width);
    const float fog =
        std::max(0.0f, 1.0f - projected.ray_distance * prepared.inv_far);
    const float near = near_field_coverage(
        projected.ray_distance, prepared.near_start, prepared.near_inv_span);
    const TraceHit hit{edge * fog * fog * near, projected.ray_distance,
                       line.free_axis};
    uint8_t position = hit_count;
    while (position > 0 && hits[position - 1].distance > hit.distance) {
      hits[position] = hits[position - 1];
      --position;
    }
    hits[position] = hit;
    ++hit_count;
  }
  for (uint8_t index = 0; index < hit_count; ++index)
    if (!consume(hits[index]))
      return;
}

template <bool PROJECTED, typename ConsumeFn>
__attribute__((always_inline)) inline void
trace_layers_mode(const Vector &normal, const PreparedTrace &prepared,
                  ConsumeFn consume) {
  HS_PROFILE_DEEP(hl_trace_layers);
  if constexpr (PROJECTED) {
    trace_projected_layers(normal, prepared, consume);
    return;
  }
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
  const uint8_t shell_count = static_cast<uint8_t>(prepared.params.shells) + 1;
  TraceCursor cursors[DIMENSIONS];
  for (int axis = 0; axis < DIMENSIONS; ++axis) {
    TraceCursor &cursor = cursors[axis];
    cursor.active = false;
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
      TraceHit candidate = trace_plane(ray_origin, direction, axis,
                                       cursor.distance, cursor.step, prepared);
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
  if (prepared.mode == LatticeMode::FOUR_D_PROJECTED)
    trace_layers_mode<true>(normal, prepared, consume);
  else
    trace_layers_mode<false>(normal, prepared, consume);
}

inline TraceHit trace(const Vector &normal, const PreparedTrace &prepared) {
  TraceHit nearest;
  trace_layers(normal, prepared, [&](const TraceHit &hit) {
    nearest = hit;
    return false;
  });
  return nearest;
}

struct LayerComposite {
  float red = 0.0f;
  float green = 0.0f;
  float blue = 0.0f;
  float remaining = 1.0f;

  void add(const Pixel &color, float coverage) {
    const float weight = remaining * coverage;
    red += static_cast<float>(color.r) * weight;
    green += static_cast<float>(color.g) * weight;
    blue += static_cast<float>(color.b) * weight;
    remaining *= 1.0f - coverage;
  }

  Color4 finish() const {
    const float alpha = 1.0f - remaining;
    if (alpha <= 0.0f)
      return {};
    const float inverse_alpha = 1.0f / alpha;
    const Pixel color{static_cast<uint16_t>(hs::clamp(
                          red * inverse_alpha + 0.5f, 0.0f, 65535.0f)),
                      static_cast<uint16_t>(hs::clamp(
                          green * inverse_alpha + 0.5f, 0.0f, 65535.0f)),
                      static_cast<uint16_t>(hs::clamp(
                          blue * inverse_alpha + 0.5f, 0.0f, 65535.0f))};
    return {color, alpha};
  }
};

template <bool PROJECTED>
__attribute__((always_inline)) inline Color4
shade_mode(const Pullback::SphereSample &input, const FrameState &frame,
           const PreparedTrace &prepared) {
  HS_PROFILE_DEEP(hl_shade);
  LayerComposite composite;
  const BakedPalette &palette =
      *(prepared.params.color == ColorMode::DEPTH ? frame.depth_palette
                                                  : frame.axis_palette);
  trace_layers_mode<PROJECTED>(
      input.dir, prepared,
      [&](const TraceHit &hit) __attribute__((always_inline)) {
        HS_PROFILE_DEEP(hl_layer_composite);
        const float depth = hit.distance * prepared.inv_far;
        const float value =
            prepared.params.color == ColorMode::DEPTH
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
  if (prepared.mode == LatticeMode::FOUR_D_PROJECTED)
    return shade_mode<true>(input, frame, prepared);
  return shade_mode<false>(input, frame, prepared);
}

template <bool PROJECTED>
struct ModeShadeStage
    : Pullback::Stage::Contract<ModeShadeStage<PROJECTED>,
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
    return shade_mode<PROJECTED>(input, frame, prepared);
  }
};

struct ShadeStage
    : Pullback::Stage::Contract<ShadeStage, Pullback::SphereSample, Color4> {
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
    return shade(input, frame, prepared);
  }
};

using RenderPipeline = Pullback::Pipeline<Binding, ShadeStage>;
using SliceRenderPipeline = Pullback::Pipeline<Binding, ModeShadeStage<false>>;
using ProjectedRenderPipeline =
    Pullback::Pipeline<Binding, ModeShadeStage<true>>;

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
  static constexpr Segue::Lerp PRESET_SEGUE{240, ease_in_out_sin,
                                            /*pausable=*/true};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 320;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 5;

  static constexpr Params preset_params(size_t index) {
    Params value;
    switch (index) {
    case 1:
      value.mode = LatticeMode::FOUR_D_SLICE;
      value.sphere_radius = 0.0f;
      value.wire_radius = 0.03546f;
      value.softness = 0.029612f;
      value.far_cells = 7.264f;
      value.aa_strength = 1.0f;
      value.speed = 0.03f;
      value.spin_3d = 0.01089f;
      value.spin_4d = 0.015f;
      value.chrome_warp = 0.65f;
      value.reflection = ReflectionMode::CHROME;
      value.color = ColorMode::DEPTH;
      value.shells = ShellCount::THREE;
      break;
    default:
      break;
    }
    return value;
  }

  static constexpr bool valid_params(const Params &value) {
    return static_cast<uint8_t>(value.mode) <=
               static_cast<uint8_t>(LatticeMode::FOUR_D_PROJECTED) &&
           value.sphere_radius >= 0.0f && value.sphere_radius <= 1.5f &&
           value.wire_radius >= 0.015f && value.wire_radius <= 0.18f &&
           value.softness >= 0.002f && value.softness <= 0.08f &&
           value.far_cells >= 2.0f && value.far_cells <= 16.0f &&
           value.aa_strength >= 0.0f && value.aa_strength <= 2.0f &&
           value.speed >= 0.0f && value.speed <= 0.05f &&
           value.spin_3d >= 0.0f && value.spin_3d <= 0.015f &&
           value.spin_4d >= 0.0f && value.spin_4d <= 0.015f &&
           value.chrome_warp >= 0.0f && value.chrome_warp <= 1.0f &&
           static_cast<uint8_t>(value.reflection) <=
               static_cast<uint8_t>(ReflectionMode::RADIAL) &&
           static_cast<uint8_t>(value.color) <=
               static_cast<uint8_t>(ColorMode::AXIS) &&
           static_cast<uint8_t>(value.shells) <=
               static_cast<uint8_t>(ShellCount::THREE);
  }

  HS_COLD_MEMBER HyperLattice() : Choreography(W, H, {.strobe = true}) {}

  void init() override {
    begin_choreography();
    register_animated_param("Dimension", &params.mode, MODE_OPTIONS,
                            MODE_EXPORT_OPTIONS, std::size(MODE_OPTIONS));
    register_animated_param("Sphere Radius", &params.sphere_radius, 0.0f, 1.5f);
    register_animated_param("Wire Radius", &params.wire_radius, 0.015f, 0.18f);
    register_animated_param("Softness", &params.softness, 0.002f, 0.08f);
    register_animated_param("Far Cells", &params.far_cells, 2.0f, 16.0f);
    register_animated_param("AA Strength", &params.aa_strength, 0.0f, 2.0f);
    register_animated_param("Speed", &params.speed, 0.0f, 0.05f);
    register_animated_param("3D Spin", &params.spin_3d, 0.0f, 0.015f);
    register_animated_param("4D Spin", &params.spin_4d, 0.0f, 0.015f);
    register_animated_param("Chrome Warp", &params.chrome_warp, 0.0f, 1.0f);
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
      if (frame.ctx.params.mode != LatticeMode::FOUR_D_PROJECTED) {
        Scan::Shader::draw_cached<W, H, 1>(
            canvas, [&frame](const Vector &view) {
              return HyperLatticeDetail::SliceRenderPipeline::evaluate(
                  view, frame.ctx, frame.prepared);
            });
      } else {
        Scan::Shader::draw_cached<W, H, 1>(
            canvas, [&frame](const Vector &view) HS_HOT_FLASH_MEMBER {
              return HyperLatticeDetail::ProjectedRenderPipeline::evaluate(
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

  static constexpr int PALETTE_FADE_FRAMES = 960;

  static constexpr const char *MODE_OPTIONS[] = {"3D", "Dimensional Rift",
                                                 "4D Slice", "4D Projected"};
  static constexpr const char *MODE_EXPORT_OPTIONS[] = {
      "LatticeMode::THREE_D", "LatticeMode::DIMENSIONAL_RIFT",
      "LatticeMode::FOUR_D_SLICE", "LatticeMode::FOUR_D_PROJECTED"};
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

  HyperLatticeDetail::Vec4 origin{{0.17f, 0.31f, 0.43f, 0.59f}};
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
