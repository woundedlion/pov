/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file spacetime.h
 * @brief Frame motion, sampling surfaces, poses and time-aware fields.
 */

#include "core/render/spacetime/surface.h"
#include "core/render/spacetime/worldlines.h"
#include "core/render/spacetime/field.h"
#include "core/render/spacetime/pose_table.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <span>

namespace Spacetime {

/** @brief Shortest quaternion arc with frame-invariant trigonometry. */
struct TemporalRotation : math::SlerpArc {
  /** @brief Prepares the shortest arc between unit orientations.
   * @param a Previous unit orientation.
   * @param b Current unit orientation. */
  HS_FLASH_INLINE TemporalRotation(const math::Quaternion &a,
                                   const math::Quaternion &b)
      : math::SlerpArc(a, b) {}
};

/** @brief Most history nodes capture() accepts. */
inline constexpr int MAX_CAPTURE_NODES = 16;

/** @brief Frame fractions at which a display arm passes a moving point. */
struct Captures {
  static constexpr int CAPACITY = 2;      ///< Most arm passes in one frame.
  std::array<float, CAPACITY> fraction{}; ///< Frame fraction per pass.
  int count = 0;                          ///< Live passes.
};

/** @brief Tagged construction of a pipeline with explicit-time entry points. */
struct Timed {};

/** @brief Work performed by geometry-callback capture and rasterization. */
struct GeometryCounters {
  uint32_t column_path_evaluations =
      0; /**< Columns evaluating their edge's endpoint paths. */
  uint32_t column_shader_calls =
      0; /**< Source shaders for emitted column samples. */
  uint32_t column_boundary_samples =
      0; /**< Additional arm-boundary evaluations. */
  uint32_t column_single_samples = 0; /**< Samples emitted through full AA. */
  uint32_t column_tap_samples = 0; /**< Fast-motion samples emitted per tap. */
  uint32_t column_draw_setups =
      0; /**< Shared column-render contexts prepared. */
  uint32_t column_mesh_culls = 0; /**< Mesh bounding-cap decisions. */
  uint32_t column_mesh_rejected =
      0; /**< Meshes rejected before vertex caching. */
  uint32_t column_edges_prepared =
      0; /**< Visible edges with prepared coverage. */
  uint32_t column_edges_culled =
      0; /**< Source edges rejected before per-edge setup. */
  uint32_t column_projections = 0; /**< Column-loop screen projections. */
  uint32_t column_generated_samples =
      0; /**< Samples advanced before clipping. */
  uint32_t column_samples =
      0; /**< Coverage samples emitted in owned columns. */
  uint32_t columns_visited = 0;  /**< Destination columns examined. */
  uint32_t column_crossings = 0; /**< Meridian intersections found. */
  uint32_t column_span_fills =
      0; /**< Columns receiving multiple coverage samples. */
  uint32_t midpoint_evaluations =
      0; /**< Source vertices evaluated at half-frame. */
  uint32_t endpoint_evaluations =
      0; /**< Source vertices evaluated at the frame's start. */
};

/** @brief Host instrumentation for the geometry-callback route. */
inline GeometryCounters geometry_counters;

/** @brief Increments a geometry counter in instrumented builds.
 * @tparam field Counter to increment.
 * @param amount Work items performed. */
template <uint32_t GeometryCounters::*field>
HS_HOT_INLINE inline void count_geometry(uint32_t amount = 1) {
#if HS_ENABLE_TEST_HOOKS
  geometry_counters.*field += amount;
#else
  (void)amount;
#endif
}

/**
 * @brief Solves `u = (theta(u) / pi) mod 1` for a point sampled at evenly
 *        spaced frame fractions.
 * @details Azimuth is linear between nodes. The arm-relative position
 * `theta / pi - u` must fall strictly from node to node; each integer it
 * crosses is one pass, located by binary search over the nodes.
 * @param nodes The point at fractions `i / (size - 1)`; at least two and at
 *        most `MAX_CAPTURE_NODES`, unit length.
 * @return Zero, one or two passes with u in [0, 1). A rise between nodes,
 * from a swing near a pole or motion with the sweep faster than the arm, or a
 * fall of two or more, returns one pass at u = 1.
 */
HS_FLASH_MEMBER inline Captures capture(std::span<const math::Vector> nodes) {
  assert(nodes.size() >= 2 && nodes.size() <= MAX_CAPTURE_NODES);
  const int last = static_cast<int>(nodes.size()) - 1;
  std::array<float, MAX_CAPTURE_NODES> relative;
  float azimuth = math::half_turns(nodes[0]);
  float unwrapped = azimuth;
  relative[0] = unwrapped;
  Captures out;
  for (int i = 1; i <= last; ++i) {
    const float next = math::half_turns(nodes[i]);
    unwrapped += std::remainder(next - azimuth, 2.0f);
    azimuth = next;
    relative[i] = unwrapped - static_cast<float>(i) / last;
    if (!(relative[i] < relative[i - 1])) {
      out.fraction[0] = 1.0f;
      out.count = 1;
      return out;
    }
  }
  if (!(relative[0] - relative[last] < Captures::CAPACITY)) {
    out.fraction[0] = 1.0f;
    out.count = 1;
    return out;
  }
  const auto begin = relative.begin();
  const auto end = begin + last + 1;
  for (float n = std::floor(relative[0]); n > relative[last]; n -= 1.0f) {
    const int j = static_cast<int>(
        std::partition_point(begin, end, [n](float r) { return r >= n; }) -
        begin);
    const int i = j - 1;
    const float t = (relative[i] - n) / (relative[i] - relative[j]);
    out.fraction[out.count++] = (static_cast<float>(i) + t) / last;
  }
  return out;
}

/**
 * @brief Arm passes of a point carried along an orientation's sub-frame path.
 * @tparam CAP Orientation capacity.
 * @param orientation History of at least two frames spanning one frame.
 * @param v Unrotated point, unit length.
 * @return capture() over @p v rotated by each history frame.
 */
template <int CAP>
HS_FLASH_MEMBER Captures capture(const math::Orientation<CAP> &orientation,
                                 const math::Vector &v) {
  static_assert(CAP <= MAX_CAPTURE_NODES);
  std::array<math::Vector, CAP> nodes;
  for (int i = 0; i < orientation.length(); ++i)
    nodes[i] = math::rotate(v, orientation.get(i));
  return capture(std::span<const math::Vector>(
      nodes.data(), static_cast<size_t>(orientation.length())));
}

} // namespace Spacetime
