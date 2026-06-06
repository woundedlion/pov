/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/hankin.h.
 *
 * Coverage:
 *   - compile_hankin builds baseVertices, staticVertices, dynamicVertices,
 *     dynamicInstructions, and face arrays consistently.
 *   - update_hankin populates an output PolyMesh for both flat (angle=0)
 *     and twisted (angle≠0) configurations.
 *   - one-shot MeshOps::hankin convenience wrapper produces a valid mesh.
 *   - CompiledHankin::clone makes an independent deep copy.
 */
#pragma once

#include <cstdint>
#include "core/hankin.h"
#include "core/solids.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace hankin_tests {

inline uint8_t hankin_target_buf[256 * 1024];
inline uint8_t hankin_temp_buf[256 * 1024];

// Build a PolyMesh from a Solids descriptor into the given arena.
template <typename Solid>
inline void build_solid(PolyMesh &mesh, Arena &arena) {
  mesh.initialize(arena, Solid::NUM_VERTS, Solid::NUM_FACES,
                  Solid::faces.size());
  for (const auto &v : Solid::vertices) mesh.vertices.push_back(v);
  for (auto fc : Solid::face_counts) mesh.face_counts.push_back(fc);
  for (auto fi : Solid::faces)
    mesh.faces.push_back(static_cast<uint16_t>(fi));
}

inline void check_face_counts_consistent(const PolyMesh &m) {
  size_t total = 0;
  for (size_t i = 0; i < m.face_counts.size(); ++i)
    total += m.face_counts[i];
  HS_EXPECT_EQ(total, m.faces.size());
}

inline void check_indices_in_range(const PolyMesh &m) {
  size_t V = m.vertices.size();
  for (size_t i = 0; i < m.faces.size(); ++i)
    HS_EXPECT_TRUE(m.faces[i] < V);
}

inline void check_all_unit_vertices(const PolyMesh &m) {
  for (size_t i = 0; i < m.vertices.size(); ++i) {
    float len = m.vertices[i].length();
    HS_EXPECT_NEAR(len, 1.0f, 1e-2f);
  }
}

// ---------------------------------------------------------------------------
// compile_hankin
// ---------------------------------------------------------------------------

inline void test_compile_hankin_populates_arrays() {
  Arena target(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, target, temp);

  // baseVertices mirrors the input vertex list exactly
  HS_EXPECT_EQ(compiled.baseVertices.size(), cube.vertices.size());
  for (size_t i = 0; i < cube.vertices.size(); ++i) {
    HS_EXPECT_NEAR(compiled.baseVertices[i].x, cube.vertices[i].x, 1e-6f);
    HS_EXPECT_NEAR(compiled.baseVertices[i].y, cube.vertices[i].y, 1e-6f);
    HS_EXPECT_NEAR(compiled.baseVertices[i].z, cube.vertices[i].z, 1e-6f);
  }

  // staticVertices: one per unique edge of the input. Cube has 12 edges.
  HS_EXPECT_EQ(compiled.staticVertices.size(), (size_t)12);

  // Static vertices are midpoints normalised onto the unit sphere.
  for (size_t i = 0; i < compiled.staticVertices.size(); ++i) {
    HS_EXPECT_NEAR(compiled.staticVertices[i].length(), 1.0f, 1e-3f);
  }

  // dynamicVertices: one per half-edge (= twice the edge count).
  HS_EXPECT_EQ(compiled.dynamicVertices.size(), (size_t)24);
  // dynamicInstructions parallels dynamicVertices.
  HS_EXPECT_EQ(compiled.dynamicInstructions.size(),
               compiled.dynamicVertices.size());

  // staticOffset == |staticVertices| (dynamic indices start at staticOffset).
  HS_EXPECT_EQ((size_t)compiled.staticOffset,
               compiled.staticVertices.size());

  // Face data is internally consistent.
  size_t total = 0;
  for (size_t i = 0; i < compiled.face_counts.size(); ++i)
    total += compiled.face_counts[i];
  HS_EXPECT_EQ(total, compiled.faces.size());

  // All face indices must be < |staticVertices| + |dynamicVertices|.
  size_t max_v = compiled.staticVertices.size() + compiled.dynamicVertices.size();
  for (size_t i = 0; i < compiled.faces.size(); ++i)
    HS_EXPECT_TRUE(compiled.faces[i] < max_v);
}

inline void test_compile_hankin_instruction_indices_in_range() {
  Arena target(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, target, temp);

  size_t V = cube.vertices.size();
  size_t S = compiled.staticVertices.size();
  for (size_t i = 0; i < compiled.dynamicInstructions.size(); ++i) {
    const auto &ins = compiled.dynamicInstructions[i];
    HS_EXPECT_TRUE(ins.vCorner < V);
    HS_EXPECT_TRUE(ins.vPrev < V);
    HS_EXPECT_TRUE(ins.vNext < V);
    HS_EXPECT_TRUE(ins.idxM1 < S);
    HS_EXPECT_TRUE(ins.idxM2 < S);
  }
}

// ---------------------------------------------------------------------------
// update_hankin — flat (angle = 0) collapses dynamic verts to the
// normalised corner positions.
// ---------------------------------------------------------------------------

inline void test_update_hankin_flat_collapses_to_corners() {
  Arena target(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, target, temp);

  // angle = 0 → each dynamic vertex equals the normalised vCorner.
  PolyMesh out;
  MeshOps::update_hankin(compiled, out, target, /*angle*/ 0.0f);

  for (size_t i = 0; i < compiled.dynamicInstructions.size(); ++i) {
    const auto &ins = compiled.dynamicInstructions[i];
    Vector expected = compiled.baseVertices[ins.vCorner].normalized();
    HS_EXPECT_NEAR(compiled.dynamicVertices[i].x, expected.x, 1e-4f);
    HS_EXPECT_NEAR(compiled.dynamicVertices[i].y, expected.y, 1e-4f);
    HS_EXPECT_NEAR(compiled.dynamicVertices[i].z, expected.z, 1e-4f);
  }
}

inline void test_update_hankin_populates_output_mesh() {
  Arena target(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, target, temp);

  PolyMesh out;
  MeshOps::update_hankin(compiled, out, target, /*angle*/ 0.4f);

  // Output vertices = staticVertices + dynamicVertices
  HS_EXPECT_EQ(out.vertices.size(),
               compiled.staticVertices.size() +
                   compiled.dynamicVertices.size());
  HS_EXPECT_EQ(out.face_counts.size(), compiled.face_counts.size());
  HS_EXPECT_EQ(out.faces.size(), compiled.faces.size());

  check_face_counts_consistent(out);
  check_indices_in_range(out);
}

// ---------------------------------------------------------------------------
// hankin() one-shot wrapper
// ---------------------------------------------------------------------------

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

  // Static + dynamic vertices live on (or very near) the unit sphere.
  check_all_unit_vertices(out);
}

inline void test_hankin_flat_and_twisted_differ() {
  Arena target_a(hankin_target_buf, sizeof(hankin_target_buf) / 2);
  Arena target_b(hankin_target_buf + sizeof(hankin_target_buf) / 2,
                 sizeof(hankin_target_buf) / 2);
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube_a;
  build_solid<Solids::Cube>(cube_a, temp);
  PolyMesh flat = MeshOps::hankin(cube_a, target_a, temp, 0.0f);

  PolyMesh cube_b;
  build_solid<Solids::Cube>(cube_b, temp);
  PolyMesh twisted = MeshOps::hankin(cube_b, target_b, temp, 0.6f);

  // Same topology, possibly different geometry.
  HS_EXPECT_EQ(flat.vertices.size(), twisted.vertices.size());
  HS_EXPECT_EQ(flat.face_counts.size(), twisted.face_counts.size());

  // The dynamic vertices should differ — find at least one with a
  // non-trivial delta. Static vertices (midpoints) are angle-independent.
  bool found_diff = false;
  for (size_t i = 0; i < flat.vertices.size(); ++i) {
    Vector d = flat.vertices[i] - twisted.vertices[i];
    if (d.length() > 1e-3f) {
      found_diff = true;
      break;
    }
  }
  HS_EXPECT_TRUE(found_diff);
}

// ---------------------------------------------------------------------------
// CompiledHankin::clone — deep copy into a separate arena
// ---------------------------------------------------------------------------

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

  HS_EXPECT_EQ(dst.baseVertices.size(), src.baseVertices.size());
  HS_EXPECT_EQ(dst.staticVertices.size(), src.staticVertices.size());
  HS_EXPECT_EQ(dst.dynamicVertices.size(), src.dynamicVertices.size());
  HS_EXPECT_EQ(dst.dynamicInstructions.size(), src.dynamicInstructions.size());
  HS_EXPECT_EQ(dst.face_counts.size(), src.face_counts.size());
  HS_EXPECT_EQ(dst.faces.size(), src.faces.size());
  HS_EXPECT_EQ(dst.staticOffset, src.staticOffset);

  // Backing storage is independent (pointers differ — different arenas).
  HS_EXPECT_TRUE(dst.baseVertices.data() != src.baseVertices.data());
  HS_EXPECT_TRUE(dst.staticVertices.data() != src.staticVertices.data());

  // Element values match.
  for (size_t i = 0; i < src.baseVertices.size(); ++i) {
    HS_EXPECT_NEAR(dst.baseVertices[i].x, src.baseVertices[i].x, 1e-6f);
    HS_EXPECT_NEAR(dst.baseVertices[i].y, src.baseVertices[i].y, 1e-6f);
    HS_EXPECT_NEAR(dst.baseVertices[i].z, src.baseVertices[i].z, 1e-6f);
  }
  for (size_t i = 0; i < src.faces.size(); ++i) {
    HS_EXPECT_EQ(dst.faces[i], src.faces[i]);
  }
}

inline void test_compiled_hankin_clear() {
  Arena arena(hankin_target_buf, sizeof(hankin_target_buf));
  Arena temp(hankin_temp_buf, sizeof(hankin_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, arena, temp);
  HS_EXPECT_TRUE(compiled.baseVertices.size() > 0);

  compiled.clear();
  HS_EXPECT_EQ(compiled.baseVertices.size(), (size_t)0);
  HS_EXPECT_EQ(compiled.staticVertices.size(), (size_t)0);
  HS_EXPECT_EQ(compiled.dynamicVertices.size(), (size_t)0);
  HS_EXPECT_EQ(compiled.dynamicInstructions.size(), (size_t)0);
  HS_EXPECT_EQ(compiled.face_counts.size(), (size_t)0);
  HS_EXPECT_EQ(compiled.faces.size(), (size_t)0);
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

inline int run_hankin_tests() {
  auto scope = hs_test::begin_module("hankin");

  test_compile_hankin_populates_arrays();
  test_compile_hankin_instruction_indices_in_range();

  test_update_hankin_flat_collapses_to_corners();
  test_update_hankin_populates_output_mesh();

  test_hankin_one_shot_produces_valid_mesh();
  test_hankin_flat_and_twisted_differ();

  test_compiled_hankin_clone_deep_copies();
  test_compiled_hankin_clear();

  return hs_test::end_module(scope);
}

} // namespace hankin_tests
} // namespace hs_test

