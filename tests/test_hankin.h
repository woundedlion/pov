/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/mesh/hankin.h.
 *
 * Coverage:
 *   - compile_hankin builds base_vertices, static_vertices, dynamic_vertices,
 *     dynamic_instructions, and face arrays consistently.
 *   - update_hankin populates an output PolyMesh for both flat (angle=0)
 *     and twisted (angle≠0) configurations.
 *   - one-shot MeshOps::hankin convenience wrapper produces a valid mesh.
 *   - the far-star guard keeps star points local at a resonance angle where
 *     contact planes go near-parallel.
 *   - CompiledHankin::clone makes an independent deep copy.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include "core/mesh/hankin.h"
#include "core/mesh/solids.h"
#include "tests/mesh_test_util.h"
#include "tests/test_conway.h" // check_euler_genus0, check_consistent_winding
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace hankin_tests {

inline uint8_t hankin_target_buf[256 * 1024];
inline uint8_t hankin_temp_buf[256 * 1024];

// ---------------------------------------------------------------------------
// compile_hankin
// ---------------------------------------------------------------------------

/**
 * @brief Verifies compile_hankin fills every CompiledHankin array with the
 *        expected sizes and self-consistent contents for a cube input.
 */
inline void test_compile_hankin_populates_arrays() {
  Arena target(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, target, temp);

  HS_EXPECT_EQ(compiled.base_vertices.size(), cube.vertices.size());
  for (size_t i = 0; i < cube.vertices.size(); ++i) {
    HS_EXPECT_NEAR(compiled.base_vertices[i].x, cube.vertices[i].x, 1e-6f);
    HS_EXPECT_NEAR(compiled.base_vertices[i].y, cube.vertices[i].y, 1e-6f);
    HS_EXPECT_NEAR(compiled.base_vertices[i].z, cube.vertices[i].z, 1e-6f);
  }

  // one static vertex per unique edge; cube has 12 edges
  HS_EXPECT_EQ(compiled.static_vertices.size(), (size_t)12);

  // static vertices are edge midpoints normalised onto the unit sphere
  for (size_t i = 0; i < compiled.static_vertices.size(); ++i) {
    HS_EXPECT_NEAR(compiled.static_vertices[i].length(), 1.0f, 1e-3f);
  }

  // one dynamic vertex per half-edge (twice the edge count)
  HS_EXPECT_EQ(compiled.dynamic_vertices.size(), (size_t)24);
  HS_EXPECT_EQ(compiled.dynamic_instructions.size(),
               compiled.dynamic_vertices.size());

  HS_EXPECT_EQ((size_t)compiled.static_offset,
               compiled.static_vertices.size());

  size_t total = 0;
  for (size_t i = 0; i < compiled.face_counts.size(); ++i)
    total += compiled.face_counts[i];
  HS_EXPECT_EQ(total, compiled.faces.size());

  size_t max_v = compiled.static_vertices.size() + compiled.dynamic_vertices.size();
  for (size_t i = 0; i < compiled.faces.size(); ++i)
    HS_EXPECT_TRUE(compiled.faces[i] < max_v);
}

/**
 * @brief Verifies every dynamic instruction references base-vertex and
 *        static-vertex indices within their respective array bounds.
 */
inline void test_compile_hankin_instruction_indices_in_range() {
  Arena target(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, target, temp);

  size_t V = cube.vertices.size();
  size_t S = compiled.static_vertices.size();
  for (size_t i = 0; i < compiled.dynamic_instructions.size(); ++i) {
    const auto &ins = compiled.dynamic_instructions[i];
    HS_EXPECT_TRUE(ins.v_corner < V);
    HS_EXPECT_TRUE(ins.v_prev < V);
    HS_EXPECT_TRUE(ins.v_next < V);
    HS_EXPECT_TRUE(ins.idx_m1 < S);
    HS_EXPECT_TRUE(ins.idx_m2 < S);
  }
}

/**
 * @brief Verifies each static vertex is the normalised midpoint of the input
 *        edge it represents, not merely some unit-length point.
 * @details The size and unit-length checks above accept any twelve points on the
 *          sphere; this pins their geometry. Every dynamic instruction names a
 *          corner and its two flanking edges (idx_m1 = edge v_prev–v_corner,
 *          idx_m2 = edge v_corner–v_next), so recompute each edge's normalised
 *          midpoint straight from base_vertices and require the referenced static
 *          vertex to match. Walking all instructions touches every edge midpoint.
 */
inline void test_compile_hankin_static_vertices_are_edge_midpoints() {
  Arena target(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, target, temp);

  HS_EXPECT_TRUE(compiled.dynamic_instructions.size() > 0);
  for (size_t i = 0; i < compiled.dynamic_instructions.size(); ++i) {
    const auto &ins = compiled.dynamic_instructions[i];

    const Vector p_corner = compiled.base_vertices[ins.v_corner];
    const Vector p_prev = compiled.base_vertices[ins.v_prev];
    const Vector p_next = compiled.base_vertices[ins.v_next];

    const Vector exp_m1 = ((p_prev + p_corner) * 0.5f).normalized();
    const Vector exp_m2 = ((p_corner + p_next) * 0.5f).normalized();

    const Vector got_m1 = compiled.static_vertices[ins.idx_m1];
    const Vector got_m2 = compiled.static_vertices[ins.idx_m2];

    HS_EXPECT_NEAR(got_m1.x, exp_m1.x, 1e-5f);
    HS_EXPECT_NEAR(got_m1.y, exp_m1.y, 1e-5f);
    HS_EXPECT_NEAR(got_m1.z, exp_m1.z, 1e-5f);
    HS_EXPECT_NEAR(got_m2.x, exp_m2.x, 1e-5f);
    HS_EXPECT_NEAR(got_m2.y, exp_m2.y, 1e-5f);
    HS_EXPECT_NEAR(got_m2.z, exp_m2.z, 1e-5f);
  }
}

/**
 * @brief Exercises compile_hankin on the icosahedron (triangular faces,
 *        5-valent vertices) so the star-face and rosette-orbit arithmetic is
 *        covered on non-quad geometry, not just the cube's 4-valent faces.
 * @details Icosahedron: 12 vertices, 30 edges, 20 triangular faces, 60
 *          half-edges. static_vertices == edge count (30), dynamic_vertices ==
 *          half-edge count (60), and the same self-consistency invariants hold.
 */
inline void test_compile_hankin_icosahedron_triangular_faces() {
  Arena target(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh ico;
  build_solid<Solids::Icosahedron>(ico, temp);

  CompiledHankin compiled;
  MeshOps::compile_hankin(ico, compiled, target, temp);

  HS_EXPECT_EQ(compiled.base_vertices.size(), (size_t)12);

  // one static vertex per unique edge; icosahedron has 30 edges
  HS_EXPECT_EQ(compiled.static_vertices.size(), (size_t)30);
  for (size_t i = 0; i < compiled.static_vertices.size(); ++i)
    HS_EXPECT_NEAR(compiled.static_vertices[i].length(), 1.0f, 1e-3f);

  // one dynamic vertex per half-edge (twice the edge count = 60)
  HS_EXPECT_EQ(compiled.dynamic_vertices.size(), (size_t)60);
  HS_EXPECT_EQ(compiled.dynamic_instructions.size(),
               compiled.dynamic_vertices.size());

  HS_EXPECT_EQ((size_t)compiled.static_offset,
               compiled.static_vertices.size());

  size_t total = 0;
  for (size_t i = 0; i < compiled.face_counts.size(); ++i)
    total += compiled.face_counts[i];
  HS_EXPECT_EQ(total, compiled.faces.size());

  size_t max_v =
      compiled.static_vertices.size() + compiled.dynamic_vertices.size();
  for (size_t i = 0; i < compiled.faces.size(); ++i)
    HS_EXPECT_TRUE(compiled.faces[i] < max_v);

  size_t V = ico.vertices.size();
  size_t S = compiled.static_vertices.size();
  for (size_t i = 0; i < compiled.dynamic_instructions.size(); ++i) {
    const auto &ins = compiled.dynamic_instructions[i];
    HS_EXPECT_TRUE(ins.v_corner < V);
    HS_EXPECT_TRUE(ins.v_prev < V);
    HS_EXPECT_TRUE(ins.v_next < V);
    HS_EXPECT_TRUE(ins.idx_m1 < S);
    HS_EXPECT_TRUE(ins.idx_m2 < S);
  }
}

/**
 * @brief Verifies compile_hankin's emission-order contract: exactly one star
 *        face per base face, emitted first, in base-face order, before any
 *        rosette faces.
 * @tparam Solid Seed solid descriptor to compile.
 * @details The bookend hankin↔base palette mapping is the identity because of
 *          this order. Star face fi is pinned to base face fi structurally: it
 *          has 2x the base face's sides, alternates midpoint (< static_offset)
 *          and star-point (>= static_offset) entries starting with a midpoint,
 *          and its star points' instruction corners are exactly base face fi's
 *          corner set. The remaining faces are rosettes (one per base vertex),
 *          distinguished by their leading star-point entry.
 */
template <typename Solid>
inline void check_star_faces_first_in_base_face_order() {
  Arena target(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh base;
  build_solid<Solid>(base, temp);

  CompiledHankin compiled;
  MeshOps::compile_hankin(base, compiled, target, temp);

  const size_t F = base.face_counts.size();
  const size_t V = base.vertices.size();
  const auto static_offset = static_cast<uint16_t>(compiled.static_offset);
  HS_EXPECT_EQ(compiled.face_counts.size(), F + V);

  size_t off = 0;
  size_t base_off = 0;
  for (size_t fi = 0; fi < F; ++fi) {
    const int bc = base.face_counts[fi];
    const int oc = compiled.face_counts[fi];
    HS_EXPECT_EQ(oc, 2 * bc);

    std::vector<uint16_t> corners;
    for (int k = 0; k < oc; ++k) {
      const uint16_t idx = compiled.faces[off + k];
      if (k % 2 == 0) {
        HS_EXPECT_TRUE(idx < static_offset);
      } else {
        HS_EXPECT_TRUE(idx >= static_offset);
        corners.push_back(
            compiled.dynamic_instructions[idx - static_offset].v_corner);
      }
    }

    std::vector<uint16_t> expected;
    for (int k = 0; k < bc; ++k)
      expected.push_back(base.faces[base_off + k]);
    std::sort(corners.begin(), corners.end());
    std::sort(expected.begin(), expected.end());
    HS_EXPECT_EQ(corners.size(), expected.size());
    for (size_t k = 0; k < corners.size() && k < expected.size(); ++k)
      HS_EXPECT_EQ(corners[k], expected[k]);

    off += oc;
    base_off += bc;
  }

  // Everything after the F star faces is a rosette: reversed emission puts a
  // star-point entry first, never a midpoint.
  for (size_t fi = F; fi < compiled.face_counts.size(); ++fi) {
    HS_EXPECT_TRUE(compiled.faces[off] >= static_offset);
    off += compiled.face_counts[fi];
  }
}

/**
 * @brief Pins the star-first, base-face-order emission on a quad seed (cube)
 *        and a triangle seed (icosahedron).
 */
inline void test_compile_hankin_star_faces_first_in_base_face_order() {
  check_star_faces_first_in_base_face_order<Solids::Cube>();
  check_star_faces_first_in_base_face_order<Solids::Icosahedron>();
}

// ---------------------------------------------------------------------------
// update_hankin
// ---------------------------------------------------------------------------

/**
 * @brief Verifies that at angle 0 each dynamic vertex lands on the normalised
 *        corner position.
 */
inline void test_update_hankin_flat_collapses_to_corners() {
  Arena target(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, target, temp);

  PolyMesh out;
  MeshOps::update_hankin(compiled, out, target, /*angle*/ 0.0f);

  for (size_t i = 0; i < compiled.dynamic_instructions.size(); ++i) {
    const auto &ins = compiled.dynamic_instructions[i];
    Vector expected = compiled.base_vertices[ins.v_corner].normalized();
    HS_EXPECT_NEAR(compiled.dynamic_vertices[i].x, expected.x, 1e-4f);
    HS_EXPECT_NEAR(compiled.dynamic_vertices[i].y, expected.y, 1e-4f);
    HS_EXPECT_NEAR(compiled.dynamic_vertices[i].z, expected.z, 1e-4f);
  }
}

/**
 * @brief Verifies update_hankin's degenerate-edge fallback (hankin.h:316-320):
 *        a zero-length corner edge at a NON-zero angle still collapses the
 *        dynamic vertex onto the normalised corner and stays finite/unit-length.
 * @details The clean platonic solids never produce a degenerate edge, so the
 *          `dot(cross,cross) < EPS_CROSS_SQ` branch is otherwise uncovered. Here
 *          we compile a normal cube, then force one instruction's previous corner
 *          to coincide with its corner so cross(p_prev, p_corner) == 0, emulating
 *          a malformed mesh. Without the guard this path would feed a zero vector
 *          into normalized() and emit NaN.
 */
inline void test_update_hankin_degenerate_edge_collapses_to_corner() {
  Arena target(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, target, temp);
  HS_EXPECT_TRUE(compiled.dynamic_instructions.size() > 0);

  // Make instruction 0's previous corner coincide with its corner so
  // cross(p_prev, p_corner) == 0, forcing a zero-length edge.
  const size_t idx = 0;
  const auto ins = compiled.dynamic_instructions[idx];
  HS_EXPECT_TRUE(ins.v_prev != ins.v_corner);
  compiled.base_vertices[ins.v_prev] = compiled.base_vertices[ins.v_corner];

  // Non-zero angle to take the twisted path, not the flat short-circuit.
  PolyMesh out;
  MeshOps::update_hankin(compiled, out, target, /*angle*/ 0.5f);

  const Vector corner = compiled.base_vertices[ins.v_corner].normalized();
  const Vector got = compiled.dynamic_vertices[idx];
  HS_EXPECT_NEAR(got.x, corner.x, 1e-4f);
  HS_EXPECT_NEAR(got.y, corner.y, 1e-4f);
  HS_EXPECT_NEAR(got.z, corner.z, 1e-4f);

  // Stays finite and unit-length — never NaN from normalising a zero.
  HS_EXPECT_TRUE(std::isfinite(got.x) && std::isfinite(got.y) &&
                 std::isfinite(got.z));
  HS_EXPECT_NEAR(got.length(), 1.0f, 1e-4f);
}

/**
 * @brief Verifies update_hankin at a non-zero angle emits an output mesh whose
 *        vertex, face-count, and face arrays match the compiled sizes and stay
 *        consistent.
 */
inline void test_update_hankin_populates_output_mesh() {
  Arena target(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, target, temp);

  PolyMesh out;
  MeshOps::update_hankin(compiled, out, target, /*angle*/ 0.4f);

  HS_EXPECT_EQ(out.vertices.size(),
               compiled.static_vertices.size() +
                   compiled.dynamic_vertices.size());
  HS_EXPECT_EQ(out.face_counts.size(), compiled.face_counts.size());
  HS_EXPECT_EQ(out.faces.size(), compiled.faces.size());

  check_face_counts_consistent(out);
  check_indices_in_range(out);
}

// ---------------------------------------------------------------------------
// hankin() one-shot wrapper
// ---------------------------------------------------------------------------

/**
 * @brief Verifies the one-shot hankin() wrapper returns a non-empty,
 *        self-consistent mesh whose vertices all lie on (or near) the unit
 *        sphere.
 */
inline void test_hankin_one_shot_produces_valid_mesh() {
  Arena target(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  PolyMesh out = MeshOps::hankin(cube, target, temp, /*angle*/ 0.3f);

  HS_EXPECT_TRUE(out.vertices.size() > 0);
  HS_EXPECT_TRUE(out.face_counts.size() > 0);
  HS_EXPECT_TRUE(out.faces.size() > 0);
  check_face_counts_consistent(out);
  check_indices_in_range(out);

  // Looser 1e-2: dynamic vertices are ray/segment intersections that land
  // slightly off the unit sphere by construction.
  check_all_unit_vertices(out, 1e-2f);
}

/**
 * @brief Verifies flat (angle 0) and twisted (angle > 0) meshes share topology
 *        but differ geometrically in at least one vertex.
 */
inline void test_hankin_flat_and_twisted_differ() {
  Arena target_a(hankin_target_buf, sizeof(hankin_target_buf) / 2);
  Arena target_b(hankin_target_buf + sizeof(hankin_target_buf) / 2,
                 sizeof(hankin_target_buf) / 2);
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  // Indices [0, static_count) are angle-independent edge midpoints; the
  // remainder are angle-dependent star points.
  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, target_a, temp);
  const size_t static_count = compiled.static_vertices.size();

  PolyMesh flat;
  MeshOps::update_hankin(compiled, flat, target_a, 0.0f);
  PolyMesh twisted;
  MeshOps::update_hankin(compiled, twisted, target_b, 0.6f);

  HS_EXPECT_EQ(flat.vertices.size(), twisted.vertices.size());
  HS_EXPECT_EQ(flat.face_counts.size(), twisted.face_counts.size());

  // Edge midpoints are angle-independent, so identical across angles.
  for (size_t i = 0; i < static_count; ++i)
    HS_EXPECT_NEAR((flat.vertices[i] - twisted.vertices[i]).length(), 0.0f,
                   1e-6f);

  float max_dynamic_delta = 0.0f;
  for (size_t i = static_count; i < flat.vertices.size(); ++i) {
    const float d = (flat.vertices[i] - twisted.vertices[i]).length();
    if (d > max_dynamic_delta)
      max_dynamic_delta = d;
  }
  HS_EXPECT_TRUE(max_dynamic_delta > 1e-2f);
}

// ---------------------------------------------------------------------------
// Topology of compiled Hankin output
// ---------------------------------------------------------------------------

/**
 * @brief Verifies Hankin output is a closed genus-0 manifold with consistent
 *        winding, for both a quad seed (cube) and a triangle seed (icosahedron).
 * @details The structural smoke above (face-count consistency, index range,
 *          loose unit-length) accepts a non-manifold or backwards-wound mesh.
 *          Hankin is the family most likely to open a seam or invert a face, so
 *          apply the same Euler + manifold-edge-degree + winding oracle Conway
 *          and solids use: every undirected edge bounds exactly two faces,
 *          V - E + F == 2, every face points outward, and no directed edge is
 *          reused in the same direction.
 */
inline void test_hankin_output_is_genus0_manifold() {
  {
    Arena target(hankin_target_buf, sizeof(hankin_target_buf));
    Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));
    PolyMesh cube;
    build_solid<Solids::Cube>(cube, temp);
    PolyMesh out = MeshOps::hankin(cube, target, temp, /*angle*/ 0.4f);
    conway_tests::check_euler_genus0(out);
    conway_tests::check_consistent_winding(out);
  }
  {
    Arena target(hankin_target_buf, sizeof(hankin_target_buf));
    Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));
    PolyMesh ico;
    build_solid<Solids::Icosahedron>(ico, temp);
    PolyMesh out = MeshOps::hankin(ico, target, temp, /*angle*/ 0.4f);
    conway_tests::check_euler_genus0(out);
    conway_tests::check_consistent_winding(out);
  }
}

// ---------------------------------------------------------------------------
// Far-star-point guard
// ---------------------------------------------------------------------------

inline uint8_t hankin_reso_a[512 * 1024];
inline uint8_t hankin_reso_b[512 * 1024];
inline uint8_t hankin_reso_target[512 * 1024];

/**
 * @brief Verifies star points stay local at a resonance angle.
 * @details The dodecahedron hk35/ambo/hk62/ambo/relax prefix at a 43-degree
 *          contact angle puts one corner class's contact planes near-parallel,
 *          so their ray intersections land ~64 degrees from the corner.
 *          Without the far-star guard the output grows sliver faces whose
 *          longest edge is ~24x the median; with it every edge stays within
 *          the healthy ratio.
 */
inline void test_update_hankin_resonance_star_points_stay_local() {
  Arena target(hankin_reso_target, sizeof(hankin_reso_target));
  PolyMesh prefix;
  {
    Arena a(hankin_reso_a, sizeof(hankin_reso_a));
    Arena b(hankin_reso_b, sizeof(hankin_reso_b));
    prefix = Solids::finalize_solid(
        Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
            .hankin(35.0f * Solids::IslamicStarPatterns::D2R)
            .ambo()
            .hankin(62.0f * Solids::IslamicStarPatterns::D2R)
            .ambo()
            .relax(100)
            .build(),
        target);
  }
  Arena out_arena(hankin_reso_a, sizeof(hankin_reso_a));
  Arena temp(hankin_reso_b, sizeof(hankin_reso_b));
  PolyMesh out = MeshOps::hankin(prefix, out_arena, temp, 43.0f * Solids::IslamicStarPatterns::D2R);

  std::vector<float> edges;
  size_t off = 0;
  for (size_t f = 0; f < out.face_counts.size(); ++f) {
    int n = out.face_counts[f];
    for (int i = 0; i < n; ++i) {
      Vector u = out.vertices[out.faces[off + i]].normalized();
      Vector v = out.vertices[out.faces[off + (i + 1) % n]].normalized();
      edges.push_back(std::acos(std::max(-1.0f, std::min(1.0f, dot(u, v)))));
    }
    off += n;
  }
  std::sort(edges.begin(), edges.end());
  float median = edges[edges.size() / 2];
  float max = edges.back();
  HS_EXPECT_TRUE(max <= 6.0f * median);
}

// ---------------------------------------------------------------------------
// CompiledHankin::clone
// ---------------------------------------------------------------------------

/**
 * @brief Verifies clone() produces a structurally equal copy with independent
 *        backing storage (distinct data pointers) and matching element values.
 */
inline void test_compiled_hankin_clone_deep_copies() {
  Arena src_arena(hankin_target_buf, sizeof(hankin_target_buf) / 2);
  Arena dst_arena(hankin_target_buf + sizeof(hankin_target_buf) / 2,
                  sizeof(hankin_target_buf) / 2);
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  CompiledHankin src;
  MeshOps::compile_hankin(cube, src, src_arena, temp);

  CompiledHankin dst;
  CompiledHankin::clone(src, dst, dst_arena);

  HS_EXPECT_EQ(dst.base_vertices.size(), src.base_vertices.size());
  HS_EXPECT_EQ(dst.static_vertices.size(), src.static_vertices.size());
  HS_EXPECT_EQ(dst.dynamic_vertices.size(), src.dynamic_vertices.size());
  HS_EXPECT_EQ(dst.dynamic_instructions.size(), src.dynamic_instructions.size());
  HS_EXPECT_EQ(dst.face_counts.size(), src.face_counts.size());
  HS_EXPECT_EQ(dst.faces.size(), src.faces.size());
  HS_EXPECT_EQ(dst.static_offset, src.static_offset);

  HS_EXPECT_TRUE(dst.base_vertices.data() != src.base_vertices.data());
  HS_EXPECT_TRUE(dst.static_vertices.data() != src.static_vertices.data());

  for (size_t i = 0; i < src.base_vertices.size(); ++i) {
    HS_EXPECT_NEAR(dst.base_vertices[i].x, src.base_vertices[i].x, 1e-6f);
    HS_EXPECT_NEAR(dst.base_vertices[i].y, src.base_vertices[i].y, 1e-6f);
    HS_EXPECT_NEAR(dst.base_vertices[i].z, src.base_vertices[i].z, 1e-6f);
  }
  for (size_t i = 0; i < src.faces.size(); ++i) {
    HS_EXPECT_EQ(dst.faces[i], src.faces[i]);
  }
}

/**
 * @brief Verifies clear() resets every CompiledHankin array back to empty.
 */
inline void test_compiled_hankin_clear() {
  Arena arena(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, arena, temp);
  HS_EXPECT_TRUE(compiled.base_vertices.size() > 0);

  compiled.clear();
  HS_EXPECT_EQ(compiled.base_vertices.size(), (size_t)0);
  HS_EXPECT_EQ(compiled.static_vertices.size(), (size_t)0);
  HS_EXPECT_EQ(compiled.dynamic_vertices.size(), (size_t)0);
  HS_EXPECT_EQ(compiled.dynamic_instructions.size(), (size_t)0);
  HS_EXPECT_EQ(compiled.face_counts.size(), (size_t)0);
  HS_EXPECT_EQ(compiled.faces.size(), (size_t)0);
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

/**
 * @brief Runs all hankin test cases under one module scope.
 * @return Number of test failures recorded by the module.
 */
inline int run_hankin_tests() {
  hs_test::ModuleFixture fixture("hankin");

  test_compile_hankin_populates_arrays();
  test_compile_hankin_instruction_indices_in_range();
  test_compile_hankin_static_vertices_are_edge_midpoints();
  test_compile_hankin_icosahedron_triangular_faces();
  test_compile_hankin_star_faces_first_in_base_face_order();

  test_update_hankin_flat_collapses_to_corners();
  test_update_hankin_degenerate_edge_collapses_to_corner();
  test_update_hankin_populates_output_mesh();

  test_hankin_one_shot_produces_valid_mesh();
  test_hankin_flat_and_twisted_differ();

  test_hankin_output_is_genus0_manifold();

  test_update_hankin_resonance_star_points_stay_local();

  test_compiled_hankin_clone_deep_copies();
  test_compiled_hankin_clear();

  return fixture.result();
}

} // namespace hankin_tests
} // namespace hs_test

