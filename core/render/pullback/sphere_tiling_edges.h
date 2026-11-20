/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file sphere_tiling_edges.h @brief Sphere tiling edge queries. */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>

#include "core/math/cube_projection.h"
#include "core/math/scalar_kernels.h"
#include "core/math/vector.h"

namespace Pullback::Source {

/// Finite great-circle arc oriented from the lower vertex index.
struct SphereTilingEdge {
  uint16_t a;          ///< Start vertex.
  uint16_t b;          ///< End vertex.
  math::Vector normal; ///< Arc plane normal.
  math::Vector from_a; ///< Arc tangent at the start.
  math::Vector from_b; ///< Reverse arc tangent at the end.
  float inv_length;    ///< Reciprocal angular length, or zero for a point.
};

/// Edge palette coordinate and coverage.
struct SphereTilingEdgeValue {
  float value;    ///< Arc fraction.
  float coverage; ///< Coverage after applying width and softness.
};

/// Prepared angular limits for an edge and its coverage fringe.
struct SphereTilingEdgeBand {
  float half_width;   ///< Opaque half-width in radians.
  float outer;        ///< Angular distance beyond which coverage is zero.
  float sine_squared; ///< Squared sine of the outer angle.

  /** @brief Prepares the edge width and outward fade.
   * @param line_thickness Full edge width in degrees.
   * @param line_softness Fade width in degrees. */
  SphereTilingEdgeBand(float line_thickness, float line_softness);
};

namespace tiling_edges_detail {

inline float chord_squared(float angle) {
  const float half = sinf(0.5f * angle);
  return 4.0f * half * half;
}

inline float chord_angle(float distance_squared) {
  return 2.0f * asinf(0.5f * sqrtf(std::min(4.0f, distance_squared)));
}

inline int cell_at(const math::Vector &direction, int steps) {
  const auto [face, u, v, major] = math::project_cube_face(direction);
  const float half_steps = 0.5f * steps;
  const float scale = half_steps / major;
  const auto x =
      hs::clamp(static_cast<int>(u * scale + half_steps), 0, steps - 1);
  const auto y =
      hs::clamp(static_cast<int>(v * scale + half_steps), 0, steps - 1);
  return (face * steps + y) * steps + x;
}

inline bool projects_to_arc(const SphereTilingEdge &edge,
                            const math::Vector &direction) {
  return edge.inv_length > 0.0f && math::dot(edge.from_a, direction) >= 0.0f &&
         math::dot(edge.from_b, direction) >= 0.0f;
}

inline float endpoint_chord_squared(const SphereTilingEdge &edge,
                                    const math::Vector &direction,
                                    std::span<const math::Vector> vertices) {
  return std::min(math::distance_squared(direction, vertices[edge.a]),
                  math::distance_squared(direction, vertices[edge.b]));
}

inline float edge_chord_squared(const SphereTilingEdge &edge,
                                const math::Vector &direction,
                                std::span<const math::Vector> vertices) {
  if (projects_to_arc(edge, direction)) {
    const float sine = math::dot(edge.normal, direction);
    const float sine_sq = std::min(1.0f, sine * sine);
    return 2.0f * sine_sq / (1.0f + sqrtf(1.0f - sine_sq));
  }
  return endpoint_chord_squared(edge, direction, vertices);
}

inline float edge_sine_squared(const SphereTilingEdge &edge,
                               const math::Vector &direction,
                               std::span<const math::Vector> vertices) {
  if (projects_to_arc(edge, direction)) {
    const float sine = math::dot(edge.normal, direction);
    return sine * sine;
  }
  const float chord =
      std::min(2.0f, endpoint_chord_squared(edge, direction, vertices));
  return chord * (1.0f - 0.25f * chord);
}

inline float edge_position(const SphereTilingEdge &edge,
                           const math::Vector &direction,
                           std::span<const math::Vector> vertices) {
  if (edge.inv_length == 0.0f)
    return 0.0f;
  const float from_start = math::dot(edge.from_a, direction);
  const float start_dot = math::dot(vertices[edge.a], direction);
  if (from_start > 0.0f && math::dot(edge.from_b, direction) > 0.0f)
    return hs::clamp(atan2f(from_start, start_dot) * edge.inv_length, 0.0f,
                     1.0f);
  return start_dot >= math::dot(vertices[edge.b], direction) ? 0.0f : 1.0f;
}

} // namespace tiling_edges_detail

inline SphereTilingEdgeBand::SphereTilingEdgeBand(float line_thickness,
                                                  float line_softness) {
  constexpr float DEGREES_TO_RADIANS = math::PI_F / 180.0f;
  const float half_width_degrees = 0.5f * line_thickness;
  half_width = half_width_degrees * DEGREES_TO_RADIANS;
  outer = (half_width_degrees + line_softness) * DEGREES_TO_RADIANS;
  const float sine = sinf(outer);
  sine_squared = sine * sine;
}

/// Non-owning cubemap index over projected mesh edges.
struct SphereTilingEdgeView {
  std::span<const math::Vector> vertices;  ///< Unit vertex directions.
  std::span<const SphereTilingEdge> edges; ///< Arcs in vertex index order.
  std::span<const uint32_t> offsets;       ///< Cell offsets into candidate_ids.
  std::span<const uint16_t> candidate_ids; ///< Candidate arc IDs.
  int steps = 0;                           ///< Cells along one cube face axis.

  /** @brief Samples edge color position and coverage.
   * @param unit_direction Direction in the solid frame.
   * @param band Prepared edge width and fringe.
   * @return Edge palette coordinate and coverage. */
  SphereTilingEdgeValue sample(const math::Vector &unit_direction,
                               const SphereTilingEdgeBand &band) const {
    using namespace tiling_edges_detail;
    if (edges.empty())
      return {};
    // Squared sine orders distances within the admitted band below 90 degrees.
    float best = band.sine_squared;
    const SphereTilingEdge *selected = nullptr;
    const int cell = cell_at(unit_direction, steps);
    const auto candidates =
        candidate_ids.subspan(offsets[cell], offsets[cell + 1] - offsets[cell]);
    for (uint16_t id : candidates) {
      const auto &edge = edges[id];
      const float distance = edge_sine_squared(edge, unit_direction, vertices);
      if (distance < best) {
        best = distance;
        selected = &edge;
      }
    }
    if (!selected)
      return {};
    const float distance = std::min(band.outer, asinf(sqrtf(best)));
    return {edge_position(*selected, unit_direction, vertices),
            1.0f - math::smooth_ramp(band.half_width, band.outer, distance)};
  }
};

} // namespace Pullback::Source
