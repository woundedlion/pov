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

inline int cell_at(const math::CubeFaceProjection &point, int steps) {
  const auto [face, u, v, major] = point;
  const float half_steps = 0.5f * steps;
  const float scale = half_steps / major;
  const auto x =
      hs::clamp(static_cast<int>(u * scale + half_steps), 0, steps - 1);
  const auto y =
      hs::clamp(static_cast<int>(v * scale + half_steps), 0, steps - 1);
  return (face * steps + y) * steps + x;
}

inline int cell_at(const math::Vector &direction, int steps) {
  return cell_at(math::project_cube_face(direction), steps);
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

inline float edge_position(const SphereTilingEdge &edge,
                           const math::Vector &direction,
                           std::span<const math::Vector> vertices,
                           bool narrow = false, bool interior = false) {
  if (edge.inv_length == 0.0f)
    return 0.0f;
  const float from_start = math::dot(edge.from_a, direction);
  const float start_dot = math::dot(vertices[edge.a], direction);
  if (from_start > 0.0f &&
      (interior || math::dot(edge.from_b, direction) > 0.0f)) {
    float angle;
    if (narrow && from_start <= 0.25f * start_dot) {
      const float t = from_start / start_dot;
      const float t2 = t * t;
      // The first omitted atan term is at most 0.25^11 / 11 radians.
      angle = t * (1.0f +
                   t2 * (-1.0f / 3.0f +
                         t2 * (1.0f / 5.0f + t2 * (-1.0f / 7.0f + t2 / 9.0f))));
    } else
      angle = atan2f(from_start, start_dot);
    return hs::clamp(angle * edge.inv_length, 0.0f, 1.0f);
  }
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

/// Narrow-band cubemap shared by the denser baked solids.
struct SphereTilingNarrowIndex {
  static constexpr int STEPS = 48; ///< Cells along one cube face axis.
  static constexpr float OUTER =
      1.5f * (math::PI_F / 180.0f); ///< Outer angle limit.
  static constexpr uint16_t INTERIOR =
      0x4000; ///< Arc contains every projected foot.
  static constexpr uint16_t LAST = 0x8000; ///< Last candidate in the cell.
  std::span<const uint16_t> cells; ///< Empty, single candidate, or list index.
  std::span<const uint16_t> candidates; ///< Multi-edge lists with flagged IDs.
};

/// Non-owning cubemap index over projected mesh edges.
struct SphereTilingEdgeView {
  std::span<const math::Vector> vertices;  ///< Unit vertex directions.
  std::span<const SphereTilingEdge> edges; ///< Arcs in vertex index order.
  std::span<const uint32_t> offsets;       ///< Cell offsets into candidate_ids.
  std::span<const uint16_t> candidate_ids; ///< Candidate arc IDs.
  int steps = 0;                           ///< Cells along one cube face axis.
  SphereTilingNarrowIndex narrow{};        ///< Optional narrow-band candidates.

  /** @brief Samples edge color position and coverage.
   * @param unit_direction Direction in the solid frame.
   * @param band Prepared edge width and fringe.
   * @return Edge palette coordinate and coverage. */
  SphereTilingEdgeValue sample(const math::Vector &unit_direction,
                               const SphereTilingEdgeBand &band) const {
    return sample(unit_direction, band,
                  math::project_cube_face(unit_direction));
  }

  /** @brief Samples edges using a shared cube-face projection.
   * @param unit_direction Direction in the solid frame.
   * @param band Prepared edge width and fringe.
   * @param projection Dominant cube face for the unit direction.
   * @return Edge palette coordinate and coverage. */
  SphereTilingEdgeValue
  sample(const math::Vector &unit_direction, const SphereTilingEdgeBand &band,
         const math::CubeFaceProjection &projection) const {
    using namespace tiling_edges_detail;
    if (edges.empty())
      return {};
    std::span<const uint16_t> candidates;
    const bool use_narrow =
        !narrow.cells.empty() && band.outer <= SphereTilingNarrowIndex::OUTER;
    uint16_t direct = 0;
    if (use_narrow) {
      const int cell = cell_at(projection, SphereTilingNarrowIndex::STEPS);
      direct = narrow.cells[cell];
      if (direct == 0)
        return {};
      if (direct & SphereTilingNarrowIndex::LAST)
        candidates = std::span<const uint16_t>(&direct, 1);
      else
        candidates = narrow.candidates.subspan(direct - 1);
    } else {
      const int cell = cell_at(projection, steps);
      const size_t begin = offsets[cell];
      candidates = candidate_ids.subspan(begin, offsets[cell + 1] - begin);
    }

    // Squared sine orders distances within the admitted band below 90 degrees.
    float best = band.sine_squared;
    float plane_limit = best * 1.00001f + 1.0e-9f;
    const SphereTilingEdge *selected = nullptr;
    bool selected_interior = false;
    const uint16_t interior_flag =
        use_narrow ? SphereTilingNarrowIndex::INTERIOR : 0;
    const uint16_t last_flag = use_narrow ? SphereTilingNarrowIndex::LAST : 0;
    const uint16_t index_mask = ~(interior_flag | last_flag);
    for (uint16_t candidate : candidates) {
      const auto &edge = edges[candidate & index_mask];
      const float sine = math::dot(edge.normal, unit_direction);
      const float plane = sine * sine;
      if (plane <= plane_limit) {
        float distance = plane;
        const bool interior = (candidate & interior_flag) != 0;
        if (!interior && !projects_to_arc(edge, unit_direction)) {
          const float chord = std::min(
              2.0f, endpoint_chord_squared(edge, unit_direction, vertices));
          distance = chord * (1.0f - 0.25f * chord);
        }
        if (distance < best) {
          best = distance;
          plane_limit = best * 1.00001f + 1.0e-9f;
          selected = &edge;
          selected_interior = interior;
        }
      }
      if (candidate & last_flag)
        break;
    }
    if (!selected)
      return {};
    const float sine = sqrtf(best);
    const float angle = use_narrow ? sine * (1.0f + best / 6.0f) : asinf(sine);
    const float distance = std::min(band.outer, angle);
    return {edge_position(*selected, unit_direction, vertices, use_narrow,
                          selected_interior),
            1.0f - math::smooth_ramp(band.half_width, band.outer, distance)};
  }
};

} // namespace Pullback::Source
