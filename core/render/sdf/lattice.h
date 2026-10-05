/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file lattice.h
 * @brief Prepared periodic lattice crossings and shading. */

#include <array>
#include "math/4dmath.h"
#include "render/ray/contract.h"
#include "render/ray/shade.h"

namespace SDF::Lattice {
constexpr int DIMENSIONS = math::VEC4_DIMENSIONS;
constexpr int MAX_SHELLS = 3;
constexpr float DIRECTION_EPSILON = 1.0e-4f;
enum class Domain : uint8_t { THREE_D, FOUR_D_SLICE };
enum class ShellCount : uint8_t { ONE, TWO, THREE };
using LatticeMode = Domain;
struct Settings {
  Domain mode = Domain::THREE_D;
  float sphere_radius = 1;
  float cell_size = 1;
  float wire_radius = .055f;
  float softness = .012f;
  float aa_strength = 1;
  ShellCount shells = ShellCount::TWO;
};
struct PreparedTrace {
  Settings params;
  math::Vec4 origin;
  math::Mat4 world_to_lattice;
  float far_distance;
  float aa_scale;
  float outer_radius_base;
  float sphere_radius_world;
  Domain mode;
};
inline PreparedTrace prepare(const Settings &settings, const math::Vec4 &origin,
                             const math::Mat4 &embedding, float far_distance,
                             float pixel_half_angle) {
  HS_CHECK(static_cast<uint8_t>(settings.shells) < MAX_SHELLS,
           "lattice shell count exceeds crossing capacity");
  HS_CHECK(
      settings.softness > 0 && settings.cell_size > 0 &&
          settings.aa_strength >= 0,
      "lattice requires positive softness and cell size and nonnegative AA");
  const float INV_CELL = 1.0f / settings.cell_size;
  PreparedTrace result{settings,
                       origin,
                       embedding,
                       far_distance,
                       settings.aa_strength * pixel_half_angle * INV_CELL *
                           INV_CELL,
                       settings.wire_radius + settings.softness,
                       settings.sphere_radius * settings.cell_size,
                       settings.mode};
  for (int i = 0; i < DIMENSIONS; ++i)
    for (int j = 0; j < DIMENSIONS; ++j)
      result.world_to_lattice.m[i][j] *= INV_CELL;
  return result;
}
/** @brief One ray's covered plane crossings, sorted by distance. */
struct CrossingList {
  static constexpr int CAPACITY = DIMENSIONS * SDF::Lattice::MAX_SHELLS;
  std::array<float, CAPACITY> distances;
  std::array<float, CAPACITY> coverages;
};

struct PreparedShading {
  SDF::Lattice::PreparedTrace lattice;
  Raycast::Appearance appearance;
  CrossingList *crossings = nullptr;
};

struct TraceHit {
  float coverage = 0;
  float distance = 0;
  uint8_t free_axis = 0;
};
__attribute__((always_inline)) inline float periodic_distance(float value) {
  return fabsf(value - nearbyintf(value));
}

struct EdgeMetric {
  float distance_sq;
  uint8_t free_axis;
};

__attribute__((always_inline)) inline float
periodic_distance_at(const math::Vec4 &ray_origin, const math::Vec4 &direction,
                     int axis, float distance) {
  return periodic_distance(ray_origin[axis] + distance * direction[axis]);
}

template <int AXIS0, int AXIS1>
__attribute__((always_inline)) EdgeMetric edge_metric_3d_axes(
    const math::Vec4 &ray_origin, const math::Vec4 &direction, float distance) {
  const float component0 =
      periodic_distance_at(ray_origin, direction, AXIS0, distance);
  const float component1 =
      periodic_distance_at(ray_origin, direction, AXIS1, distance);
  const float component0_sq = component0 * component0;
  const float component1_sq = component1 * component1;
  return {fminf(component0_sq, component1_sq),
          static_cast<uint8_t>(component1_sq > component0_sq ? AXIS1 : AXIS0)};
}

inline EdgeMetric edge_metric_3d_at(const math::Vec4 &ray_origin,
                                    const math::Vec4 &direction, int plane_axis,
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

/**
 * @brief Finds the nearest 4D edge within the supplied coverage radius.
 * @details An edge fixes two of the three remaining coordinates. Two outside
 * the radius reject the hit before evaluating its squared distance.
 */
template <int AXIS0, int AXIS1, int AXIS2, bool NEED_AXIS = true>
__attribute__((always_inline)) bool
edge_metric_4d_axes_bounded(const math::Vec4 &ray_origin,
                            const math::Vec4 &direction, float distance,
                            float limit, float limit_sq, EdgeMetric &result) {
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
    result = {sum - fmaxf(component0_sq, fmaxf(component1_sq, component2_sq)),
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
__attribute__((always_inline)) inline bool edge_metric_4d_at_bounded(
    const math::Vec4 &ray_origin, const math::Vec4 &direction, int plane_axis,
    float distance, float limit, float limit_sq, EdgeMetric &result) {
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

// Unlike flash-resident smooth_ramp, this stays inline on the hot crossing
// path.
__attribute__((always_inline)) inline float
lattice_ramp(float edge0, float edge1, float value) {
  return math::cubic_kernel((value - edge0) / (edge1 - edge0));
}

inline float wire_coverage(float metric_sq, float radius, float half_width) {
  const float signed_distance = sqrtf(metric_sq) - radius;
  return 1.0f - lattice_ramp(-half_width, half_width, signed_distance);
}

__attribute__((always_inline)) inline float
fast_wire_coverage(float metric_sq, float radius, float half_width) {
  const float signed_distance = sqrtf(metric_sq) - radius;
  const float ramp_position =
      0.5f + 0.5f * signed_distance * math::fast_reciprocal(half_width);
  return 1.0f - math::cubic_kernel(ramp_position);
}

inline float shell_horizon_coverage(uint8_t shell, uint8_t shell_count,
                                    float distance, float magnitude) {
  if (shell + 1 < shell_count)
    return 1.0f;
  // The ramp spans one unit, so it needs no normalizing division.
  return 1.0f - math::cubic_kernel(distance * magnitude -
                                   static_cast<float>(shell_count - 1));
}

inline float next_plane_offset(float origin, bool positive) {
  const float fraction = math::wrap_t(origin);
  if (fraction == 0.0f)
    return 1.0f;
  return positive ? 1.0f - fraction : fraction;
}

template <bool SLICE_4D>
__attribute__((always_inline)) inline TraceHit
trace_plane(const math::Vec4 &ray_origin, const math::Vec4 &direction,
            int plane_axis, float distance, float plane_step,
            const PreparedTrace &prepared) {
  HS_PROFILE_DEEP(hl_plane_eval);

  const float coverage_outer_radius =
      prepared.outer_radius_base + prepared.aa_scale * distance * plane_step;
  const float coverage_half_width =
      coverage_outer_radius - prepared.params.wire_radius;
  const float outer_radius_sq = coverage_outer_radius * coverage_outer_radius;
  float metric_sq;
  uint8_t free_axis;
  if constexpr (SLICE_4D) {
    EdgeMetric metric_4d;
    if (!edge_metric_4d_at_bounded<!SLICE_4D>(ray_origin, direction, plane_axis,
                                              distance, coverage_outer_radius,
                                              outer_radius_sq, metric_4d))
      return {0.0f, distance, 0};
    metric_sq = metric_4d.distance_sq;
    free_axis = metric_4d.free_axis;
  } else if (prepared.mode == Domain::THREE_D) {
    const EdgeMetric metric_3d =
        edge_metric_3d_at(ray_origin, direction, plane_axis, distance);
    metric_sq = metric_3d.distance_sq;
    free_axis = metric_3d.free_axis;
    if (metric_sq >= outer_radius_sq)
      return {0.0f, distance, free_axis};
  } else {
    EdgeMetric metric_4d;
    if (!edge_metric_4d_at_bounded(ray_origin, direction, plane_axis, distance,
                                   coverage_outer_radius, outer_radius_sq,
                                   metric_4d))
      return {0.0f, distance, 0};
    metric_sq = metric_4d.distance_sq;
    free_axis = metric_4d.free_axis;
  }

  const float edge =
      SLICE_4D ? fast_wire_coverage(metric_sq, prepared.params.wire_radius,
                                    coverage_half_width)
               : wire_coverage(metric_sq, prepared.params.wire_radius,
                               coverage_half_width);
  return {edge, distance, free_axis};
}

/** @brief Approximate plane-crossing coverage for cubic and hypercubic edges. */
template <bool SLICE_4D = false, uint8_t FIXED_SHELL_COUNT = 0> struct Events {
  static_assert(FIXED_SHELL_COUNT <= MAX_SHELLS);
  static constexpr size_t STREAM_COUNT = DIMENSIONS;
  static constexpr size_t GROUP_CAPACITY = 1;
  /** @brief Stream state; float fields are read only while active. */
  struct Cursor {
    float distance, step, magnitude;
    uint8_t shell = 0;
    bool active = false;
  };
  const PreparedTrace &prepared;
  math::Vec4 direction, origin;
  std::array<Cursor, STREAM_COUNT> cursors;
  uint8_t shell_count;

  __attribute__((always_inline)) Events(const math::Vector &normal,
                                        const PreparedTrace &prepared)
      : prepared(prepared), direction(prepared.world_to_lattice.apply(
                                {{normal.x, normal.y, normal.z, 0}})),
        origin(prepared.origin),
        shell_count(FIXED_SHELL_COUNT
                        ? FIXED_SHELL_COUNT
                        : static_cast<uint8_t>(prepared.params.shells) + 1) {
    for (int axis = 0; axis < DIMENSIONS; ++axis) {
      origin[axis] += prepared.sphere_radius_world * direction[axis];
      if (!SLICE_4D && axis == 3 && prepared.mode == Domain::THREE_D)
        continue;
      auto &cursor = cursors[axis];
      cursor.magnitude = fabsf(direction[axis]);
      if (cursor.magnitude < DIRECTION_EPSILON)
        continue;
      cursor.step = 1.0f / cursor.magnitude;
      cursor.distance =
          next_plane_offset(origin[axis], direction[axis] > 0) * cursor.step;
      cursor.active = cursor.distance < prepared.far_distance;
    }
  }
  __attribute__((always_inline)) bool active(size_t i) const {
    return cursors[i].active;
  }
  __attribute__((always_inline)) float distance(size_t i) const {
    return cursors[i].distance;
  }
  __attribute__((always_inline)) Raycast::Contribution
  candidate(size_t i) const {
    const auto &cursor = cursors[i];
    const auto hit =
        trace_plane<SLICE_4D>(origin, direction, static_cast<int>(i),
                              cursor.distance, cursor.step, prepared);
    Raycast::Contribution result;
    result.t = hit.distance;
    result.coverage = hit.coverage *
                      shell_horizon_coverage(cursor.shell, shell_count,
                                             cursor.distance, cursor.magnitude);
    result.feature = hit.free_axis;
    result.merge_identity = 0;
    return result;
  }
  __attribute__((always_inline)) void advance(size_t i) {
    auto &cursor = cursors[i];
    ++cursor.shell;
    cursor.distance += cursor.step;
    cursor.active =
        cursor.shell < shell_count && cursor.distance < prepared.far_distance;
  }
};
/**
 * @brief Composites one ray's plane crossings front to back.
 * @details Matches Raycast::shade_events over SDF::Lattice::Events. Each axis
 * evaluates its crossings with the axis fixed, and only covered crossings enter
 * the distance-ordered layer list: an uncovered crossing never closes a merge
 * group, so the groups are runs of covered crossings within the relative
 * tolerance of their first distance, each one layer at that distance with the
 * run's largest coverage. The stream capacity bounds every ray below the
 * candidate and layer budgets, so neither is tracked.
 */
template <bool SLICE_4D, uint8_t SHELLS>
__attribute__((always_inline)) inline LayerComposite
composite_crossings(const math::Vector &normal,
                    const PreparedShading &prepared) {
  static_assert(SHELLS <= MAX_SHELLS);
  HS_PROFILE_DEEP(hl_shade);
  using SDF::Lattice::DIRECTION_EPSILON;
  constexpr float RELATIVE_TOLERANCE = 1.0e-4f;
  static_assert(CrossingList::CAPACITY <= Raycast::TraceLimits{}.max_layers);
  const auto &lattice = prepared.lattice;
  const math::Vec4 direction =
      lattice.world_to_lattice.apply({{normal.x, normal.y, normal.z, 0}});
  math::Vec4 origin = lattice.origin;
  std::array<float, DIMENSIONS> magnitude;
  float product = 1.0f;
  for (int axis = 0; axis < DIMENSIONS; ++axis) {
    origin[axis] += lattice.sphere_radius_world * direction[axis];
    magnitude[axis] = !SLICE_4D && axis == 3 && lattice.mode == Domain::THREE_D
                          ? 0.0f
                          : fabsf(direction[axis]);
    if (magnitude[axis] >= DIRECTION_EPSILON)
      product *= magnitude[axis];
  }
  // One division serves every axis: each step is the product of the other
  // magnitudes over the product of all of them.
  const float INVERSE_PRODUCT = 1.0f / product;
  const uint8_t SHELL_COUNT =
      SHELLS ? SHELLS : static_cast<uint8_t>(lattice.params.shells) + 1;
  auto &distances = prepared.crossings->distances;
  auto &coverages = prepared.crossings->coverages;
  int count = 0;
  const auto cross = [&](int axis) __attribute__((always_inline)) {
    if (!(magnitude[axis] >= DIRECTION_EPSILON))
      return;
    float step = INVERSE_PRODUCT;
    for (int other = 0; other < DIMENSIONS; ++other)
      if (other != axis && magnitude[other] >= DIRECTION_EPSILON)
        step *= magnitude[other];
    float distance =
        SDF::Lattice::next_plane_offset(origin[axis], direction[axis] > 0) *
        step;
    for (uint8_t shell = 0;
         shell < SHELL_COUNT && distance < lattice.far_distance;
         ++shell, distance += step) {
      HS_PROFILE_DEEP(hl_event_step);
      const float COVERAGE =
          SDF::Lattice::trace_plane<SLICE_4D>(origin, direction, axis, distance,
                                              step, lattice)
              .coverage *
          SDF::Lattice::shell_horizon_coverage(shell, SHELL_COUNT, distance,
                                               magnitude[axis]);
      if (!(COVERAGE > 0.0f))
        continue;
      int slot = count++;
      for (; slot > 0 && distances[slot - 1] > distance; --slot) {
        distances[slot] = distances[slot - 1];
        coverages[slot] = coverages[slot - 1];
      }
      distances[slot] = distance;
      coverages[slot] = COVERAGE;
    }
  };
  cross(0);
  cross(1);
  cross(2);
  cross(3);
  LayerComposite composite;
  for (int layer = 0; layer < count;) {
    HS_PROFILE_DEEP(hl_layer_composite);
    const float T = distances[layer];
    const float GROUP_END = T + RELATIVE_TOLERANCE * fmaxf(1.0f, T);
    float coverage = coverages[layer];
    for (++layer; layer < count && distances[layer] <= GROUP_END; ++layer)
      coverage = fmaxf(coverage, coverages[layer]);
    prepared.appearance.composite(composite, T, coverage);
    if (composite.saturated())
      break;
  }
  return composite;
}

} // namespace SDF::Lattice
