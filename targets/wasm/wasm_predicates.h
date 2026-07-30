/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file wasm_predicates.h
 * @brief Pure (no-Emscripten) boundary predicates for the WASM bridge.
 *
 * The JS frontend passes untyped integers across the embind boundary; wasm.cpp
 * validates and clamps them before they reach engine code that would otherwise
 * trap or run unbounded. Those checks are plain scalar arithmetic with no
 * Emscripten dependency, so they live here and are host-unit-testable without an
 * Emscripten toolchain — see tests/test_wasm_predicates.h. wasm.cpp keeps only
 * the logging/embind shell on top.
 */
#pragma once

#include <cstddef>
#include <cstdint>

namespace hs_wasm {

/**
 * @brief Validates a clip band against the canvas extent.
 * @param x0 Inclusive left column.
 * @param x1 Exclusive right column.
 * @param y0 Inclusive top row.
 * @param y1 Exclusive bottom row.
 * @param width Canvas width; x1 must not exceed it.
 * @param height Canvas height; y1 must not exceed it.
 * @return true iff the band is non-negative, ordered, and within the canvas.
 * @details Negatives would feed ClipRegion's modulo arithmetic; a transposed
 *          (y-first) JS call must fail the range check rather than clip the
 *          wrong axis. Rejecting here keeps the untyped boundary from trapping
 *          the whole WASM module.
 */
inline bool clip_bounds_valid(int x0, int x1, int y0, int y1, int width,
                              int height) {
  return x0 >= 0 && y0 >= 0 && x0 <= x1 && x1 <= width && y0 <= y1 &&
         y1 <= height;
}

/**
 * @brief Clamps a relax iteration count into [0, max_iterations].
 * @param iterations Requested pass count from the JS boundary.
 * @param max_iterations Inclusive upper bound.
 * @return The clamped count.
 * @details relax(1e9) would freeze the main thread for billions of passes, so
 *          the unbounded JS count is clamped rather than trusted. Negative
 *          requests floor at 0.
 */
inline int clamp_relax_iterations(int iterations, int max_iterations) {
  if (iterations < 0)
    return 0;
  if (iterations > max_iterations)
    return max_iterations;
  return iterations;
}

/**
 * @brief True when a [0,1] boundary fraction falls outside its domain.
 */
inline bool unit_fraction_out_of_range(float t) { return t < 0.0f || t > 1.0f; }

/**
 * @brief Clamps a [0,1] boundary fraction into range.
 * @param t Requested fraction from the JS boundary.
 * @return t clamped to [0, 1]; a NaN passes through unchanged.
 * @details truncate/bevel/chamfer/snub feed their fraction to an always-on
 *          HS_CHECK, so an out-of-range value from a direct/API caller would
 *          abort the whole module. Callers reject non-finite args first.
 */
inline float clamp_unit_fraction(float t) {
  if (t < 0.0f)
    return 0.0f;
  if (t > 1.0f)
    return 1.0f;
  return t;
}

/**
 * @brief True when a mesh's element counts exceed a tooling ceiling.
 * @param verts Mesh vertex count.
 * @param faces Mesh face count.
 * @param indices Mesh flat face-index count.
 * @param max_elements Ceiling each of the three counts must stay within.
 * @return true when any of the three counts exceeds the ceiling.
 * @details A JS-driven operator chain can grow a mesh past what the tooling
 *          scratch arenas hold or what the engine's 16-bit connectivity can
 *          address. Rejecting it at the boundary keeps the engine's fail-fast
 *          traps out of reach of the JS caller.
 */
inline bool tooling_mesh_over_ceiling(size_t verts, size_t faces,
                                      size_t indices, size_t max_elements) {
  return verts > max_elements || faces > max_elements || indices > max_elements;
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
 * @details A ceiling on the input alone does not bound the intermediates a
 *          composition builds — bevel = truncate of ambo reaches six times the
 *          input index count, and gyro's snub stage five — so the admissible
 *          input scales down by the operator's expansion. Division, not
 *          multiplication, so the prediction itself cannot overflow.
 */
inline bool mesh_op_expansion_over_ceiling(size_t verts, size_t faces,
                                           size_t indices, size_t expansion,
                                           size_t max_elements) {
  return mesh_largest_element_count(verts, faces, indices) >
         max_elements / expansion;
}

/**
 * @brief True when a Hankin contact angle falls outside its [0, max] domain.
 * @param radians Contact angle from the JS boundary.
 * @param max_radians Inclusive upper bound of the operator's domain (pi/2).
 * @return true when the angle is out of domain.
 * @details The Hankin construction is 2*pi-periodic in the angle and mirrors
 *          under negation, so an out-of-domain angle silently aliases onto an
 *          in-domain pattern instead of failing. Callers reject non-finite args
 *          first.
 */
inline bool hankin_angle_out_of_range(float radians, float max_radians) {
  return radians < 0.0f || radians > max_radians;
}

/**
 * @brief True when a gradient-shape index falls outside [lo, hi].
 * @details bakeLut casts the JS int into the GradientShape enum; an out-of-range
 *          value is UB, so the boundary clamps it to the default shape (lo).
 */
inline bool gradient_shape_out_of_range(int shape, int lo, int hi) {
  return shape < lo || shape > hi;
}

/**
 * @brief Clamps an out-of-range gradient-shape index to lo (the default shape).
 */
inline int clamp_gradient_shape(int shape, int lo, int hi) {
  return gradient_shape_out_of_range(shape, lo, hi) ? lo : shape;
}

/**
 * @brief True when an untyped HSV key byte falls outside [0, 255].
 */
inline bool hsv_key_out_of_range(int v) { return v < 0 || v > 255; }

/**
 * @brief Clamps an untyped HSV key into [0, 255] before the uint8_t narrowing.
 * @details Without the clamp the uint8_t cast would wrap mod 256, turning an
 *          out-of-range hue/sat/value into a wrong-but-valid byte.
 */
inline uint8_t clamp_hsv_key(int v) {
  if (v < 0)
    return 0;
  if (v > 255)
    return 255;
  return static_cast<uint8_t>(v);
}

} // namespace hs_wasm
