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
#include <cstdint>

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

/** @brief Tagged construction of a pipeline with explicit-time entry points. */
struct Timed {};

/** @brief Work performed by geometry-callback capture and rasterization. */
struct GeometryCounters {
  /** @brief Point chords certified without projection curvature. */
  uint32_t point_small_chords = 0;
  /** @brief Closed-form coherent anchor captures. */
  uint32_t chord_capture_solves = 0;
  /** @brief Closed-form point-piece capture calls. */
  uint32_t point_capture_solves = 0;
  /** @brief Coherent moving-chord crossings. */
  uint32_t chord_fast_samples = 0;
  /** @brief Moving chart anchor intervals examined. */
  uint32_t chord_intervals = 0;
  /** @brief Invisible anchor intervals skipped. */
  uint32_t chord_intervals_culled = 0;
  /** @brief Moving screen-chord samples emitted. */
  uint32_t chord_screen_samples = 0;
  /** @brief Supported columns examining polar geometry. */
  uint32_t polar_column_visits = 0;
  /** @brief Shared polar-run segment captures. */
  uint32_t polar_geometry_samples = 0;
  /** @brief Polar chart coverage samples. */
  uint32_t polar_column_samples = 0;
  /** @brief Point trajectory knots projected. */
  uint32_t point_piece_evaluations = 0;
  /** @brief Columns evaluating their edge's endpoint paths. */
  uint32_t column_path_evaluations = 0;
  /** @brief Source shaders for emitted column samples. */
  uint32_t column_shader_calls = 0;
  /** @brief Samples entering the source-history route. */
  uint32_t column_history_visits = 0;
  /** @brief Additional arm-boundary evaluations. */
  uint32_t column_boundary_samples = 0;
  /** @brief Samples emitted through full AA. */
  uint32_t column_single_samples = 0;
  /** @brief Fast-motion samples emitted per tap. */
  uint32_t column_tap_samples = 0;
  /** @brief Shared column-render contexts prepared. */
  uint32_t column_draw_setups = 0;
  /** @brief Mesh bounding-cap decisions. */
  uint32_t column_mesh_culls = 0;
  /** @brief Meshes rejected before vertex caching. */
  uint32_t column_mesh_rejected = 0;
  /** @brief Visible edges with prepared coverage. */
  uint32_t column_edges_prepared = 0;
  /** @brief Source edges rejected before per-edge setup. */
  uint32_t column_edges_culled = 0;
  /** @brief Column-loop screen projections. */
  uint32_t column_projections = 0;
  /** @brief Samples advanced before clipping. */
  uint32_t column_generated_samples = 0;
  /** @brief Coverage samples emitted in owned columns. */
  uint32_t column_samples = 0;
  /** @brief Destination columns examined. */
  uint32_t columns_visited = 0;
  /** @brief Meridian intersections found. */
  uint32_t column_crossings = 0;
  /** @brief Columns receiving multiple coverage samples. */
  uint32_t column_span_fills = 0;
  /** @brief Source vertices evaluated at half-frame. */
  uint32_t midpoint_evaluations = 0;
  /** @brief Source vertices evaluated at the frame's start. */
  uint32_t endpoint_evaluations = 0;
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

} // namespace Spacetime
