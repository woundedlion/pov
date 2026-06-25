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

} // namespace hs_wasm
