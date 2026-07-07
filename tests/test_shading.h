/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Direct unit tests for core/render/shading.h — the Fragment register carrier
 * (lerp across pos/v0-v3/age/size/color), the edge-distance and topology-slot
 * helpers shared by the mesh effects, and the null vertex/fragment shaders.
 * These helpers are otherwise exercised only incidentally through the rasterizer
 * tests; this module pins their contract directly.
 *
 * Self-contained header. run_shading_tests() returns the module failure count.
 */
#pragma once

#include "core/render/shading.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace shading_tests {

// --- Fragment::lerp ---------------------------------------------------------

/**
 * @brief Verifies Fragment::lerp reproduces each endpoint at t=0 and t=1.
 */
inline void test_fragment_lerp_endpoints() {
  Fragment a;
  a.pos = Vector(1, 0, 0);
  a.v0 = 1.0f; a.v1 = 2.0f; a.v2 = 3.0f; a.v3 = 4.0f;
  a.size = 5.0f; a.age = 6.0f; a.color = Color4(Pixel(0, 0, 0), 0.2f);

  Fragment b;
  b.pos = Vector(0, 1, 0);
  b.v0 = 10.0f; b.v1 = 20.0f; b.v2 = 30.0f; b.v3 = 40.0f;
  b.size = 50.0f; b.age = 60.0f; b.color = Color4(Pixel(0, 0, 0), 0.8f);

  Fragment lo = Fragment::lerp(a, b, 0.0f);
  HS_EXPECT_NEAR(lo.v0, a.v0, 1e-6f);
  HS_EXPECT_NEAR(lo.size, a.size, 1e-6f);
  HS_EXPECT_NEAR(lo.color.alpha, a.color.alpha, 1e-6f);

  Fragment hi = Fragment::lerp(a, b, 1.0f);
  HS_EXPECT_NEAR(hi.v3, b.v3, 1e-6f);
  HS_EXPECT_NEAR(hi.age, b.age, 1e-6f);
  HS_EXPECT_NEAR(hi.color.alpha, b.color.alpha, 1e-6f);
}

/**
 * @brief Verifies the midpoint interpolates every register, so a sample between
 *        two control points carries pos/v0-v3/age/size/color rather than
 *        resetting size and color to their struct defaults.
 */
inline void test_fragment_lerp_midpoint_carries_registers() {
  Fragment a;
  a.pos = Vector(2, 0, 0);
  a.v0 = 0.0f; a.v1 = 0.0f; a.v2 = 0.0f; a.v3 = 0.0f;
  a.size = 4.0f; a.age = 0.0f; a.color = Color4(Pixel(0, 0, 0), 0.0f);

  Fragment b;
  b.pos = Vector(0, 4, 0);
  b.v0 = 8.0f; b.v1 = 12.0f; b.v2 = 16.0f; b.v3 = 20.0f;
  b.size = 8.0f; b.age = 10.0f; b.color = Color4(Pixel(0, 0, 0), 1.0f);

  Fragment m = Fragment::lerp(a, b, 0.5f);
  HS_EXPECT_NEAR(m.pos.x, 1.0f, 1e-6f);
  HS_EXPECT_NEAR(m.pos.y, 2.0f, 1e-6f);
  HS_EXPECT_NEAR(m.v0, 4.0f, 1e-6f);
  HS_EXPECT_NEAR(m.v1, 6.0f, 1e-6f);
  HS_EXPECT_NEAR(m.v2, 8.0f, 1e-6f);
  HS_EXPECT_NEAR(m.v3, 10.0f, 1e-6f);
  HS_EXPECT_NEAR(m.size, 6.0f, 1e-6f);
  HS_EXPECT_NEAR(m.age, 5.0f, 1e-6f);
  HS_EXPECT_NEAR(m.color.alpha, 0.5f, 1e-6f);
}

// --- fragment_edge_dist -----------------------------------------------------

/**
 * @brief Verifies fragment_edge_dist returns -v1/size (inward depth) for a
 *        well-sized face.
 */
inline void test_fragment_edge_dist_normal() {
  Fragment f;
  f.v1 = -0.5f; // negative = inside the face
  f.size = 2.0f;
  HS_EXPECT_NEAR(fragment_edge_dist(f), 0.25f, 1e-6f);
}

/**
 * @brief Verifies a degenerate (near-zero-size) face yields 0 rather than
 *        dividing by ~0.
 */
inline void test_fragment_edge_dist_degenerate_face() {
  Fragment f;
  f.v1 = -1.0f;
  f.size = 0.0f;
  HS_EXPECT_NEAR(fragment_edge_dist(f), 0.0f, 1e-6f);
}

// --- mesh_topology_slot -----------------------------------------------------

/**
 * @brief Verifies an in-range face index resolves to its topology class wrapped
 *        modulo the palette count.
 */
inline void test_mesh_topology_slot_in_range_wraps() {
  const int topology[] = {0, 5, 2};
  Fragment f;
  f.v2 = 1.0f; // face index 1 -> class 5
  HS_EXPECT_EQ((mesh_topology_slot<4>(f, topology, 3)), 1); // wrap(5, 4) == 1
}

/**
 * @brief Verifies an out-of-range or negative face index falls back to class 0
 *        rather than reading out of bounds.
 */
inline void test_mesh_topology_slot_out_of_range_falls_back() {
  const int topology[] = {3, 3, 3}; // class 3 everywhere
  Fragment over;
  over.v2 = 9.0f; // >= num_faces -> class 0
  HS_EXPECT_EQ((mesh_topology_slot<4>(over, topology, 3)), 0);

  Fragment neg;
  neg.v2 = -1.0f; // negative -> class 0
  HS_EXPECT_EQ((mesh_topology_slot<4>(neg, topology, 3)), 0);
}

// --- null shaders -----------------------------------------------------------

/**
 * @brief Verifies NullVertexShader leaves every fragment register untouched.
 */
inline void test_null_vertex_shader_is_identity() {
  Fragment f;
  f.pos = Vector(0.3f, -0.4f, 0.866f);
  f.v0 = 7.0f; f.size = 9.0f; f.age = 11.0f;
  f.color = Color4(Pixel(1, 2, 3), 0.5f);

  NullVertexShader()(f);
  HS_EXPECT_NEAR(f.pos.x, 0.3f, 1e-6f);
  HS_EXPECT_NEAR(f.v0, 7.0f, 1e-6f);
  HS_EXPECT_NEAR(f.size, 9.0f, 1e-6f);
  HS_EXPECT_NEAR(f.age, 11.0f, 1e-6f);
  HS_EXPECT_NEAR(f.color.alpha, 0.5f, 1e-6f);
}

/**
 * @brief Verifies NullFragmentShader emits a fully transparent color regardless
 *        of input.
 */
inline void test_null_fragment_shader_is_transparent() {
  Fragment f;
  f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
  Color4 out = NullFragmentShader()(Vector(1, 0, 0), f);
  HS_EXPECT_NEAR(out.alpha, 0.0f, 1e-6f);
  HS_EXPECT_EQ(out.color.r, 0);
  HS_EXPECT_EQ(out.color.g, 0);
  HS_EXPECT_EQ(out.color.b, 0);
}

/**
 * @brief Runs every shading test case.
 * @return The module's failure count.
 */
inline int run_shading_tests() {
  hs_test::ModuleFixture fixture("shading");

  test_fragment_lerp_endpoints();
  test_fragment_lerp_midpoint_carries_registers();
  test_fragment_edge_dist_normal();
  test_fragment_edge_dist_degenerate_face();
  test_mesh_topology_slot_in_range_wraps();
  test_mesh_topology_slot_out_of_range_falls_back();
  test_null_vertex_shader_is_identity();
  test_null_fragment_shader_is_transparent();

  return fixture.result();
}

} // namespace shading_tests
} // namespace hs_test
