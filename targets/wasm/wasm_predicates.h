/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file wasm_predicates.h
 * @brief Pure (no-Emscripten) boundary predicates for the WASM bridge.
 *
 * The JS frontend passes untyped integers across the embind boundary; wasm.cpp
 * validates and clamps them before they reach engine code that would otherwise
 * trap or run unbounded. Those checks are plain integer arithmetic with no
 * Emscripten dependency, so they live here and are host-unit-testable without an
 * Emscripten toolchain — see tests/test_wasm_predicates.h. wasm.cpp keeps only
 * the logging/embind shell on top.
 */
#pragma once

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
