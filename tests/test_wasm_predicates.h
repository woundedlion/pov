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
#include "tests/test_fixture.h"
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
  constexpr int MAX = 1000;

  // In-range counts pass through unchanged.
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(0, MAX), 0);
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(1, MAX), 1);
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(MAX, MAX), MAX);
  // Negative floors at 0.
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(-1, MAX), 0);
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(-1000000, MAX), 0);
  // Over-large counts clamp to the cap (relax(1e9) would freeze the thread).
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(MAX + 1, MAX), MAX);
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(1000000000, MAX), MAX);
}

/**
 * @brief Exercises the bakeLut gradient-shape range check and clamp.
 */
inline void check_gradient_shape_clamp() {
  constexpr int LO = 0, HI = 3; // STRAIGHT .. FALLOFF

  // In-range shapes pass through untouched and read as in range.
  for (int s = LO; s <= HI; ++s) {
    HS_EXPECT_TRUE(!hs_wasm::gradient_shape_out_of_range(s, LO, HI));
    HS_EXPECT_EQ(hs_wasm::clamp_gradient_shape(s, LO, HI), s);
  }
  // Out-of-range (below and above) is flagged and folds to the default (LO).
  HS_EXPECT_TRUE(hs_wasm::gradient_shape_out_of_range(-1, LO, HI));
  HS_EXPECT_TRUE(hs_wasm::gradient_shape_out_of_range(HI + 1, LO, HI));
  HS_EXPECT_EQ(hs_wasm::clamp_gradient_shape(-1, LO, HI), LO);
  HS_EXPECT_EQ(hs_wasm::clamp_gradient_shape(HI + 1, LO, HI), LO);
  HS_EXPECT_EQ(hs_wasm::clamp_gradient_shape(1000000, LO, HI), LO);
}

/**
 * @brief Exercises the bakeLut HSV-key range check and [0,255] clamp.
 */
inline void check_hsv_key_clamp() {
  // In-range keys, including the boundaries, pass through unchanged.
  HS_EXPECT_TRUE(!hs_wasm::hsv_key_out_of_range(0));
  HS_EXPECT_TRUE(!hs_wasm::hsv_key_out_of_range(255));
  HS_EXPECT_EQ(hs_wasm::clamp_hsv_key(0), 0);
  HS_EXPECT_EQ(hs_wasm::clamp_hsv_key(128), 128);
  HS_EXPECT_EQ(hs_wasm::clamp_hsv_key(255), 255);
  // Out-of-range keys are flagged and saturate rather than wrap mod 256.
  HS_EXPECT_TRUE(hs_wasm::hsv_key_out_of_range(-1));
  HS_EXPECT_TRUE(hs_wasm::hsv_key_out_of_range(256));
  HS_EXPECT_EQ(hs_wasm::clamp_hsv_key(-1), 0);
  HS_EXPECT_EQ(hs_wasm::clamp_hsv_key(256), 255);
  HS_EXPECT_EQ(hs_wasm::clamp_hsv_key(300), 255);
}

/**
 * @brief Module entry point: runs the boundary-predicate checks.
 * @return The module's failure count.
 */
inline int run_wasm_predicates_tests() {
  hs_test::ModuleFixture fixture("wasm_predicates");
  check_clip_bounds();
  check_relax_clamp();
  check_gradient_shape_clamp();
  check_hsv_key_clamp();
  return fixture.result();
}

} // namespace wasm_predicates_tests
} // namespace hs_test
