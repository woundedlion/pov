/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host unit tests for the WASM boundary predicates
 * (targets/wasm/wasm_predicates.h). These checks gate untyped integers crossing
 * the embind boundary before they reach engine code that would trap (clip
 * bounds) or run unbounded (relax iterations). They compile only under
 * Emscripten inside wasm.cpp, so the pure predicates are extracted and exercised
 * here without the toolchain.
 */
#pragma once

#include "targets/wasm/wasm_predicates.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace wasm_predicates_tests {

/**
 * @brief Exercises clip_bounds_valid across in-range, edge, and malformed bands.
 */
inline void check_clip_bounds() {
  constexpr int W = 96, H = 20;

  // A typical sub-canvas band is accepted.
  HS_EXPECT_TRUE(hs_wasm::clip_bounds_valid(0, 48, 0, 10, W, H));
  // Full-canvas band (exclusive ends equal the extent) is accepted.
  HS_EXPECT_TRUE(hs_wasm::clip_bounds_valid(0, W, 0, H, W, H));
  // Empty but ordered band (x0 == x1) is accepted.
  HS_EXPECT_TRUE(hs_wasm::clip_bounds_valid(5, 5, 3, 3, W, H));

  // Negative origin rejected (would feed ClipRegion modulo arithmetic).
  HS_EXPECT_TRUE(!hs_wasm::clip_bounds_valid(-1, 10, 0, 10, W, H));
  HS_EXPECT_TRUE(!hs_wasm::clip_bounds_valid(0, 10, -1, 10, W, H));
  // Inverted order rejected on each axis.
  HS_EXPECT_TRUE(!hs_wasm::clip_bounds_valid(10, 5, 0, 10, W, H));
  HS_EXPECT_TRUE(!hs_wasm::clip_bounds_valid(0, 10, 10, 5, W, H));
  // Out of canvas rejected on each axis.
  HS_EXPECT_TRUE(!hs_wasm::clip_bounds_valid(0, W + 1, 0, 10, W, H));
  HS_EXPECT_TRUE(!hs_wasm::clip_bounds_valid(0, 10, 0, H + 1, W, H));
  // A transposed (y-first) call that swaps the extents must fail the check
  // rather than clip the wrong axis.
  HS_EXPECT_TRUE(!hs_wasm::clip_bounds_valid(0, H, 0, W, W, H));
}

/**
 * @brief Exercises clamp_relax_iterations across negative, in-range, and over.
 */
inline void check_relax_clamp() {
  constexpr int kMax = 1000;

  // In-range counts pass through unchanged.
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(0, kMax), 0);
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(1, kMax), 1);
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(kMax, kMax), kMax);
  // Negative floors at 0.
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(-1, kMax), 0);
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(-1000000, kMax), 0);
  // Over-large counts clamp to the cap (relax(1e9) would freeze the thread).
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(kMax + 1, kMax), kMax);
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(1000000000, kMax), kMax);
}

/**
 * @brief Module entry point: runs the boundary-predicate checks.
 * @return The module's failure count.
 */
inline int run_wasm_predicates_tests() {
  auto scope = hs_test::begin_module("wasm_predicates");
  check_clip_bounds();
  check_relax_clamp();
  return hs_test::end_module(scope);
}

} // namespace wasm_predicates_tests
} // namespace hs_test
