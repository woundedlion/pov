/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/canvas.h"
#include "core/render/filter/screen_anti_alias.h"
#include "core/render/plot/raster_config.h"
#include "core/render/spacetime/crossing.h"
#include "core/render/spacetime/pose_table.h"
#include "core/render/spacetime/swept_bound.h"
#include "core/render/spacetime/worldlines.h"

/** @file column_geometry.h @brief Geometry callbacks, edge descriptors and meridian crossings. */

namespace Plot {

HS_O3_BEGIN

/** @brief Borrowed geometric control-point evaluator and path bounds. */
struct GeometryCallback {
  /** Evaluates one original control point at a frame fraction. */
  FunctionRef<math::Vector(size_t, float)> evaluate;
  Spacetime::PathBound path{}; /**< Derivative bounds of every control point's
                                  path over the frame. */

  bool source_positions =
      true; /**< Shaders read source-space positions under world stages;
                 false hands them world positions. */
  /** Bounds each vertex's path, in place of `path`, when set. */
  FunctionRef<Spacetime::PathBound(size_t)> vertex_path{};
  bool chained = false; /**< Edges run along polylines, each vertex starting
                             at most one edge and ending at most one other:
                             the column walk plots each edge's start but not
                             its end, so a closed chain plots every vertex once
                             and an open one all but its last. */

  bool rigid_lattice =
      false; ///< Resolve lattice and AA support per column, skipping frame-end pitch measurement.
  bool uniform_color =
      false; ///< Shader writes only color, constant over each edge and frame.
  float capture_error =
      0.025f; ///< Maximum projected vertex capture error, pixels.

  /** @brief Constructs an exact evaluator.
   * @param fn Borrowed control-point callback.
   * @return The exact geometry descriptor. */
  static GeometryCallback exact(FunctionRef<math::Vector(size_t, float)> fn) {
    return {fn};
  }
};

#if HS_ENABLE_TEST_HOOKS
/** @brief Spacing checks using the directional certificate. */
inline uint32_t g_column_lattice_checks = 0;
/** @brief Spacing checks settled by the directional certificate. */
inline uint32_t g_column_lattice_certificates = 0;
#endif

/** @brief Angular bound between a meridian crossing's arc fraction and the
 * column walk sampler's parameter, radians. */
inline constexpr float COLUMN_ARC_PARAMETER_PAD = 1e-5f;

/** @brief Largest projected motion between adjacent column times, pixels, at
 * which a column holding an edge end splats its taps at its own time. */
inline constexpr float COLUMN_SINGLE_TAP_ERROR = 0.02f;

/** @brief Intersects a geodesic basis with a forward meridian.
 * @param start_side,end_side Endpoint distances from the meridian plane.
 * @param tangent_side Start tangent distance from the meridian plane.
 * @param start_along,tangent_along Basis projections along the meridian.
 * @param inverse_length Reciprocal arc length, or zero for a collapsed edge.
 * @param fraction Receives the material arc fraction.
 * @return Whether the minor arc crosses the forward meridian. */
HS_HOT_INLINE inline bool
column_basis_crossing(float start_side, float end_side, float tangent_side,
                      float start_along, float tangent_along,
                      float inverse_length, float &fraction) {
  if (start_side * end_side > 0.0f || fabsf(start_side - end_side) < 1e-9f)
    return false;
  float angle = math::precise_atan2(-start_side, tangent_side);
  const float forward = tangent_side * start_along - start_side * tangent_along;
  if (angle < 0.0f) {
    angle += math::PI_F;
    if (forward > 0.0f)
      return false;
  } else if (forward < 0.0f)
    return false;
  fraction = hs::clamp(angle * inverse_length, 0.0f, 1.0f);
  return true;
}

/** @brief Intersects a minor arc with a destination meridian.
 * @param a,b Arc endpoint directions; a is unit length.
 * @param direction Horizontal meridian direction.
 * @param length Minor arc length.
 * @param tangent Unit tangent at the first endpoint.
 * @param fraction Resulting arc fraction.
 * @return Whether the meridian crosses the minor arc. */
HS_HOT_INLINE inline bool
column_edge_crossing(const math::Vector &a, const math::Vector &b,
                     const math::Vector &direction, float length,
                     const math::Vector &tangent, float &fraction) {
  const math::Vector normal(-direction.z, 0.0f, direction.x);
  return column_basis_crossing(
      math::dot(a, normal), math::dot(b, normal), math::dot(tangent, normal),
      math::dot(a, direction), math::dot(tangent, direction),
      length > 1e-7f ? 1.0f / length : 0.0f, fraction);
}

/** @brief A mesh edge as two vertex indices. */
struct MeshEdge {
  uint16_t u; ///< First endpoint vertex index.
  uint16_t v; ///< Second endpoint vertex index.
};

/** @brief Mesh edges sharing one column walk. */
struct EdgeBatch {
  size_t count; /**< Edge count. */
  void (*prepare)(Fragments &, const math::Vector &, const math::Vector &,
                  size_t,
                  float); /**< Writes an edge's endpoint registers into the
                                walk's fragments from its frame-end
                                vertices, index and frame-end length in radians. */
  const math::Vector
      *vertices; /**< Frame-end positions for motion.size() vertices. */
  const Spacetime::Worldlines<>
      &motion;                /**< Source captures and interpolants of the
                              vertices this frame. */
  const MeshEdge *pairs;      /**< Vertex indices of each edge. */
  bool chained = false;       /**< Edges leave their end vertex unplotted. */
  bool rigid_lattice = false; ///< Resolve lattice and AA support per column.
  bool uniform_color = false; ///< Shader writes only one color per edge.
  float age = 0.0f; ///< Age offset added to the prepared endpoint registers.
};

/** @brief Prepares endpoint fragments for a timed edge.
 * @param points Bound two-fragment output buffer.
 * @param from,to Endpoint vertices.
 * @param edge_index Edge register value.
 * @param length Prepared frame-end edge length, radians. */
FASTRUN HS_O3_FN HS_NOINLINE_NOCLONE inline void
prepare_timed_edge(Fragments &points, const math::Vector &from,
                   const math::Vector &to, size_t edge_index, float length) {
  points.clear();
  Fragment a, b;
  a.pos = from;
  b.pos = length > 0.0f ? to : from;
  a.v2 = b.v2 = static_cast<float>(edge_index);
  b.v0 = length > 0.0f ? 1.0f : 0.0f;
  b.v1 = length;
  points.push_back(a);
  points.push_back(b);
}

HS_O3_END

} // namespace Plot
