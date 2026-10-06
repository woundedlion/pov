/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file wasm_predicates.h
 * @brief Pure (no-Emscripten) boundary predicates for the WASM bridge.
 *
 * Validates and clamps untyped JS scalars before they reach engine code that
 * would otherwise trap or run unbounded.
 */
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace hs_wasm {

/** Maximum near-pole LOD aggressiveness accepted from live WASM controls. */
inline constexpr float MAX_POLE_LOD_AGGRESSIVENESS = 8.0f;

/**
 * @brief Clamps near-pole LOD aggressiveness to the live-tuning range.
 * @param aggressiveness Requested columns-per-footprint multiplier.
 * @return A finite value in [0, MAX_POLE_LOD_AGGRESSIVENESS].
 */
inline float clamp_pole_lod_aggressiveness(float aggressiveness) {
  if (std::isnan(aggressiveness) || aggressiveness < 0.0f)
    return 0.0f;
  return aggressiveness > MAX_POLE_LOD_AGGRESSIVENESS
             ? MAX_POLE_LOD_AGGRESSIVENESS
             : aggressiveness;
}

/**
 * @brief Validates a clip band against the canvas extent.
 * @param x0 Inclusive left column.
 * @param x1 Exclusive right column.
 * @param y0 Inclusive top row.
 * @param y1 Exclusive bottom row.
 * @param width Canvas width; x1 must not exceed it.
 * @param height Canvas height; y1 must not exceed it.
 * @return true iff every bound is integral and the band is non-negative,
 *         ordered, and within the canvas.
 * @details Negatives would feed ClipRegion's modulo arithmetic. Each axis is
 *          checked against its own extent. Doubles because an i32 embind
 *          parameter coerces NaN and multiples of 2^32 to 0 without a range
 *          check in a release build.
 */
inline bool clip_bounds_valid(double x0, double x1, double y0, double y1,
                              int width, int height) {
  const double bounds[4] = {x0, x1, y0, y1};
  for (const double bound : bounds) {
    // NaN fails the ordered comparison, so it needs no separate test.
    if (!(bound >= 0.0) || bound != std::floor(bound))
      return false;
  }
  return x0 <= x1 && x1 <= width && y0 <= y1 && y1 <= height;
}

/**
 * @brief Validates a preset index against the active effect's roster.
 * @param index Requested preset index from the JS boundary, as a double.
 * @param preset_count Number of presets the effect exposes.
 * @return true iff the index is integral and in [0, preset_count).
 * @details A double because a uint32_t embind parameter wraps modulo 2^32
 *          without a range check in a release build (NaN arrives as 0).
 *          Non-integral requests are rejected, not truncated.
 */
inline bool preset_index_valid(double index, size_t preset_count) {
  // NaN fails the ordered comparison, so it needs no separate test.
  if (!(index >= 0.0) || index != std::floor(index))
    return false;
  return index < static_cast<double>(preset_count);
}

/**
 * @brief Clamps a relax iteration count into [0, max_iterations].
 * @param iterations Requested pass count from the JS boundary, as a double.
 * @param max_iterations Inclusive upper bound.
 * @return The clamped count.
 * @details Negative and NaN requests floor at 0. A double because a JS number
 *          past INT32_MAX wraps negative through an i32 parameter.
 */
inline int clamp_relax_iterations(double iterations, int max_iterations) {
  if (!(iterations > 0.0))
    return 0;
  if (iterations > static_cast<double>(max_iterations))
    return max_iterations;
  return static_cast<int>(iterations);
}

/**
 * @brief True when a [0,1] boundary fraction falls outside its domain.
 */
inline bool unit_fraction_out_of_range(double t) {
  return t < 0.0f || t > 1.0f;
}

/**
 * @brief Clamps a [0,1] boundary fraction into range.
 * @param t Requested fraction from the JS boundary.
 * @return t clamped to [0, 1]; a NaN passes through unchanged.
 * @details Callers reject non-finite args first.
 */
inline float clamp_unit_fraction(double t) {
  if (t < 0.0f)
    return 0.0f;
  if (t > 1.0f)
    return 1.0f;
  return static_cast<float>(t);
}

/// Largest float strictly below 1 — the top of a half-open [0,1) domain.
inline constexpr float LARGEST_FRACTION_BELOW_ONE = 0x1.fffffep-1f;
static_assert(LARGEST_FRACTION_BELOW_ONE < 1.0f,
              "LARGEST_FRACTION_BELOW_ONE must satisfy a t < 1 domain check");

/**
 * @brief True when a fraction exceeds the representable [0,1) domain.
 */
inline bool half_open_fraction_out_of_range(double t) {
  return t < 0.0f || t > LARGEST_FRACTION_BELOW_ONE;
}

/**
 * @brief Clamps a boundary fraction into a half-open [0,1) domain.
 * @param t Requested fraction from the JS boundary.
 * @return t clamped to [0, LARGEST_FRACTION_BELOW_ONE]; a NaN passes through
 *         unchanged.
 * @details For operators that assert `t < 1.0f`. Callers reject non-finite
 *          args first.
 */
inline float clamp_half_open_fraction(double t) {
  if (t < 0.0f)
    return 0.0f;
  if (t > LARGEST_FRACTION_BELOW_ONE)
    return LARGEST_FRACTION_BELOW_ONE;
  return static_cast<float>(t);
}

/**
 * @brief Saturates a finite JS scalar to the engine float range.
 * @param value Finite value supplied by the JS caller.
 * @return A finite float, including either saturation endpoint.
 */
inline float clamp_finite_float(double value) {
  constexpr float MAX = std::numeric_limits<float>::max();
  if (value > MAX)
    return MAX;
  if (value < -MAX)
    return -MAX;
  return static_cast<float>(value);
}

/**
 * @brief The largest of a mesh's three element counts.
 * @param verts Mesh vertex count.
 * @param faces Mesh face count.
 * @param indices Mesh flat face-index count.
 */
inline size_t mesh_largest_element_count(size_t verts, size_t faces,
                                         size_t indices) {
  const size_t largest = verts > faces ? verts : faces;
  return largest > indices ? largest : indices;
}

/**
 * @brief True when a mesh operator's expansion would carry some stage past an
 *        element ceiling.
 * @param verts Input mesh vertex count.
 * @param faces Input mesh face count.
 * @param indices Input mesh flat face-index count.
 * @param expansion Largest multiple of the input's biggest element count that
 *        any of this operator's intermediate or output stages reaches; >= 1.
 * @param max_elements Ceiling every stage must stay within.
 * @return true when the operator must be rejected.
 * @details Divides rather than multiplies so the prediction cannot overflow.
 */
inline bool mesh_op_expansion_over_ceiling(size_t verts, size_t faces,
                                           size_t indices, size_t expansion,
                                           size_t max_elements) {
  if (expansion == 0)
    return true;
  return mesh_largest_element_count(verts, faces, indices) >
         max_elements / expansion;
}

/**
 * @brief True when a mesh operator's finalized output will not fit in what is
 *        left of an arena.
 * @param verts Input mesh vertex count.
 * @param faces Input mesh face count.
 * @param indices Input mesh flat face-index count.
 * @param expansion Operator's element expansion (see
 *        mesh_op_expansion_over_ceiling); >= 1.
 * @param bytes_per_element Arena bytes one output element retains.
 * @param used_bytes Bytes already committed in the arena.
 * @param capacity_bytes Arena capacity in bytes.
 * @return true when the operator must be rejected.
 * @details Prices each of the output's counts at the predicted peak, which
 *          bounds the whole mesh.
 */
inline bool mesh_op_output_over_arena(size_t verts, size_t faces,
                                      size_t indices, size_t expansion,
                                      size_t bytes_per_element,
                                      size_t used_bytes,
                                      size_t capacity_bytes) {
  if (expansion == 0 || used_bytes >= capacity_bytes)
    return true;
  const size_t remaining = capacity_bytes - used_bytes;
  return mesh_largest_element_count(verts, faces, indices) >
         remaining / (expansion * bytes_per_element);
}

/**
 * @brief The largest side count in a mesh's per-face count list.
 * @param counts Per-face side counts.
 * @param num_faces Length of @p counts.
 * @return The widest face's side count, or 0 for a mesh with no faces.
 */
inline size_t mesh_max_face_degree(const uint8_t *counts, size_t num_faces) {
  size_t largest = 0;
  for (size_t i = 0; i < num_faces; ++i) {
    if (counts[i] > largest)
      largest = counts[i];
  }
  return largest;
}

/**
 * @brief The largest number of faces meeting at any one vertex.
 * @param faces Flat per-face vertex index list.
 * @param total_indices Length of @p faces.
 * @param incidence Scratch for one counter per vertex, cleared here; must hold
 *        at least @p num_verts entries.
 * @param num_verts Vertex count; an index at or past it is skipped.
 * @return The highest incidence count, or 0 for a mesh with no faces.
 * @details On a closed manifold a vertex's incidence count is its valence.
 */
inline size_t mesh_max_vertex_valence(const uint16_t *faces,
                                      size_t total_indices, uint32_t *incidence,
                                      size_t num_verts) {
  for (size_t i = 0; i < num_verts; ++i)
    incidence[i] = 0;
  size_t largest = 0;
  for (size_t i = 0; i < total_indices; ++i) {
    const size_t v = faces[i];
    if (v >= num_verts)
      continue;
    const size_t n = ++incidence[v];
    if (n > largest)
      largest = n;
  }
  return largest;
}

/**
 * @brief True when a mesh operator would emit a face with more sides than the
 *        mesh's 8-bit per-face count can hold.
 * @param max_face_degree Widest face in the input mesh.
 * @param max_vertex_valence Highest vertex valence in the input mesh.
 * @param face_degree_factor Multiple of @p max_face_degree that the operator's
 *        widest face-derived face reaches; 0 when it emits none.
 * @param valence_factor Multiple of @p max_vertex_valence that the operator's
 *        widest vertex-derived face reaches; 0 when it emits none.
 * @param max_degree Inclusive side-count ceiling (UINT8_MAX).
 * @return true when the operator must be rejected.
 * @details Divides rather than multiplies so the prediction cannot overflow.
 */
inline bool mesh_op_face_degree_overflows(size_t max_face_degree,
                                          size_t max_vertex_valence,
                                          size_t face_degree_factor,
                                          size_t valence_factor,
                                          size_t max_degree) {
  const bool face_over = face_degree_factor != 0 &&
                         max_face_degree > max_degree / face_degree_factor;
  const bool valence_over =
      valence_factor != 0 && max_vertex_valence > max_degree / valence_factor;
  return face_over || valence_over;
}

/**
 * @brief True when a Hankin contact angle falls outside its [0, max] domain.
 * @param radians Contact angle from the JS boundary.
 * @param max_radians Inclusive upper bound of the operator's domain (pi/2).
 * @return true when the angle is out of domain.
 * @details An out-of-domain angle aliases onto an in-domain pattern. Callers
 *          reject non-finite args first.
 */
inline bool hankin_angle_out_of_range(double radians, float max_radians) {
  return radians < 0.0f || radians > max_radians;
}

} // namespace hs_wasm
