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
 * @brief Exercises the [0,1] operator-fraction range check and clamp.
 */
inline void check_unit_fraction_clamp() {
  // In-range fractions, including the boundaries, pass through unchanged.
  HS_EXPECT_TRUE(!hs_wasm::unit_fraction_out_of_range(0.0f));
  HS_EXPECT_TRUE(!hs_wasm::unit_fraction_out_of_range(0.5f));
  HS_EXPECT_TRUE(!hs_wasm::unit_fraction_out_of_range(1.0f));
  HS_EXPECT_EQ(hs_wasm::clamp_unit_fraction(0.0f), 0.0f);
  HS_EXPECT_EQ(hs_wasm::clamp_unit_fraction(0.5f), 0.5f);
  HS_EXPECT_EQ(hs_wasm::clamp_unit_fraction(1.0f), 1.0f);

  // Out-of-range fractions are flagged and saturate, so truncate/bevel/chamfer/
  // snub cannot trip their always-on HS_CHECK and abort the module.
  HS_EXPECT_TRUE(hs_wasm::unit_fraction_out_of_range(-0.001f));
  HS_EXPECT_TRUE(hs_wasm::unit_fraction_out_of_range(1.001f));
  HS_EXPECT_EQ(hs_wasm::clamp_unit_fraction(-0.001f), 0.0f);
  HS_EXPECT_EQ(hs_wasm::clamp_unit_fraction(1.001f), 1.0f);
  HS_EXPECT_EQ(hs_wasm::clamp_unit_fraction(-1e30f), 0.0f);
  HS_EXPECT_EQ(hs_wasm::clamp_unit_fraction(1e30f), 1.0f);
}

/**
 * @brief Exercises the per-operator tooling mesh-size ceiling.
 */
inline void check_tooling_mesh_ceiling() {
  constexpr size_t MAX = 65536;

  // A mesh at the ceiling on every count is still accepted.
  HS_EXPECT_TRUE(!hs_wasm::tooling_mesh_over_ceiling(0, 0, 0, MAX));
  HS_EXPECT_TRUE(!hs_wasm::tooling_mesh_over_ceiling(12, 20, 60, MAX));
  HS_EXPECT_TRUE(!hs_wasm::tooling_mesh_over_ceiling(MAX, MAX, MAX, MAX));

  // Each count is checked independently.
  HS_EXPECT_TRUE(hs_wasm::tooling_mesh_over_ceiling(MAX + 1, 20, 60, MAX));
  HS_EXPECT_TRUE(hs_wasm::tooling_mesh_over_ceiling(12, MAX + 1, 60, MAX));
  HS_EXPECT_TRUE(hs_wasm::tooling_mesh_over_ceiling(12, 20, MAX + 1, MAX));
}

/**
 * @brief Exercises the per-operator expansion ceiling.
 */
inline void check_mesh_op_expansion_ceiling() {
  constexpr size_t MAX = 65535; // the engine's 16-bit connectivity range

  HS_EXPECT_EQ(hs_wasm::mesh_largest_element_count(8, 6, 24), 24u);
  HS_EXPECT_EQ(hs_wasm::mesh_largest_element_count(70000, 6, 24), 70000u);
  HS_EXPECT_EQ(hs_wasm::mesh_largest_element_count(8, 70000, 24), 70000u);

  // A cube is far inside every operator's bound.
  for (size_t e :
       {size_t(1), size_t(2), size_t(3), size_t(4), size_t(5), size_t(6)}) {
    HS_EXPECT_TRUE(!hs_wasm::mesh_op_expansion_over_ceiling(8, 6, 24, e, MAX));
  }

  // A non-expanding operator (dual, relax) admits a mesh right at the ceiling.
  HS_EXPECT_TRUE(!hs_wasm::mesh_op_expansion_over_ceiling(1, 1, MAX, 1, MAX));
  HS_EXPECT_TRUE(
      hs_wasm::mesh_op_expansion_over_ceiling(1, 1, MAX + 1, 1, MAX));

  // The admissible input scales down by the expansion: the same mesh passes at
  // 3x (kis) and fails at 6x (meta, bevel).
  HS_EXPECT_TRUE(
      !hs_wasm::mesh_op_expansion_over_ceiling(1, 1, MAX / 3, 3, MAX));
  HS_EXPECT_TRUE(
      hs_wasm::mesh_op_expansion_over_ceiling(1, 1, MAX / 3, 6, MAX));
  HS_EXPECT_TRUE(
      !hs_wasm::mesh_op_expansion_over_ceiling(1, 1, MAX / 6, 6, MAX));

  // The bound is per-stage, so it bites one element past MAX / expansion.
  HS_EXPECT_TRUE(
      !hs_wasm::mesh_op_expansion_over_ceiling(1, 1, MAX / 5, 5, MAX));
  HS_EXPECT_TRUE(
      hs_wasm::mesh_op_expansion_over_ceiling(1, 1, MAX / 5 + 1, 5, MAX));

  // Any of the three counts can be the binding one.
  HS_EXPECT_TRUE(hs_wasm::mesh_op_expansion_over_ceiling(30000, 6, 24, 4, MAX));
  HS_EXPECT_TRUE(hs_wasm::mesh_op_expansion_over_ceiling(8, 30000, 24, 4, MAX));
  HS_EXPECT_TRUE(hs_wasm::mesh_op_expansion_over_ceiling(8, 6, 30000, 4, MAX));
}

/**
 * @brief Exercises the Hankin contact-angle domain check.
 */
inline void check_hankin_angle_domain() {
  constexpr float MAX = 1.5707963f; // pi/2
  constexpr float FIFTY_FOUR_DEG = 0.9424778f;

  // In-domain angles, including both endpoints, are accepted.
  HS_EXPECT_TRUE(!hs_wasm::hankin_angle_out_of_range(0.0f, MAX));
  HS_EXPECT_TRUE(!hs_wasm::hankin_angle_out_of_range(FIFTY_FOUR_DEG, MAX));
  HS_EXPECT_TRUE(!hs_wasm::hankin_angle_out_of_range(MAX, MAX));

  // Negative and past-pi/2 angles are rejected rather than mirroring an
  // in-domain pattern.
  HS_EXPECT_TRUE(hs_wasm::hankin_angle_out_of_range(-0.001f, MAX));
  HS_EXPECT_TRUE(hs_wasm::hankin_angle_out_of_range(MAX + 0.01f, MAX));
  HS_EXPECT_TRUE(hs_wasm::hankin_angle_out_of_range(12.0f, MAX));
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
  check_unit_fraction_clamp();
  check_tooling_mesh_ceiling();
  check_mesh_op_expansion_ceiling();
  check_hankin_angle_domain();
  check_gradient_shape_clamp();
  check_hsv_key_clamp();
  return fixture.result();
}

} // namespace wasm_predicates_tests
} // namespace hs_test
