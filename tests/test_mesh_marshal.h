/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host unit tests for the WASM mesh-editor validate-and-build layer
 * (targets/wasm/mesh_marshal.h). fromData is the most input-fragile surface in
 * the bridge — it ingests two untrusted flat arrays and must reject every
 * malformed shape rather than trap and abort the whole module. That logic lives
 * behind emscripten::val in wasm.cpp and was previously untested; these tests
 * drive the extracted pure layer directly and cover each reject branch, the
 * "arena untouched on reject" invariant, the boundary at UINT8_MAX face size,
 * and the subtle face-counting edge cases (empty faces, no trailing delimiter).
 */
#pragma once

#include "core/memory.h"
#include "core/mesh.h"
#include "targets/wasm/mesh_marshal.h"
#include "tests/test_harness.h"

#include <cstdint>
#include <vector>

namespace hs_test {
namespace mesh_marshal_tests {

using hs_wasm::build_mesh_from_flat;
using hs_wasm::MeshBuildStatus;

// A unit tetrahedron as flat (vertex, -1-delimited face) streams. Four verts,
// four triangular faces => 12 face-vertex indices.
inline void valid_tetra(std::vector<float> &v, std::vector<int> &f) {
  v = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
  f = {0, 1, 2, -1, 0, 1, 3, -1, 0, 2, 3, -1, 1, 2, 3, -1};
}

inline void test_valid_mesh_builds() {
  static uint8_t buf[1 << 16];
  Arena arena(buf, sizeof(buf));
  std::vector<float> v;
  std::vector<int> f;
  valid_tetra(v, f);

  PolyMesh m;
  MeshBuildStatus st = build_mesh_from_flat(v, f, arena, m);
  HS_EXPECT(st == MeshBuildStatus::kOk, "valid tetra builds");
  HS_EXPECT_EQ(m.vertices.size(), static_cast<size_t>(4));
  HS_EXPECT_EQ(m.face_counts.size(), static_cast<size_t>(4));
  HS_EXPECT_EQ(m.faces.size(), static_cast<size_t>(12));
  // Each face is a triangle.
  for (size_t i = 0; i < m.face_counts.size(); ++i)
    HS_EXPECT_EQ(m.face_counts[i], static_cast<uint8_t>(3));
}

// A trailing face with no final -1 must still be counted/built, and consecutive
// or leading -1 delimiters (empty faces) must emit nothing — the exact cases the
// counting pass is written to get right.
inline void test_counting_edge_cases() {
  static uint8_t buf[1 << 16];
  Arena arena(buf, sizeof(buf));
  std::vector<float> v = {0, 0, 0, 1, 0, 0, 0, 1, 0};
  // Leading -1 (empty), a real triangle, consecutive -1 (empty), a real edge
  // pair, then a trailing face with NO closing -1.
  std::vector<int> f = {-1, 0, 1, 2, -1, -1, 0, 1, -1, 1, 2};

  PolyMesh m;
  MeshBuildStatus st = build_mesh_from_flat(v, f, arena, m);
  HS_EXPECT(st == MeshBuildStatus::kOk, "counting edge cases build");
  HS_EXPECT_EQ(m.face_counts.size(), static_cast<size_t>(3)); // 3 non-empty faces
  HS_EXPECT_EQ(m.faces.size(), static_cast<size_t>(7));       // 3 + 2 + 2 indices
  HS_EXPECT_EQ(m.face_counts[0], static_cast<uint8_t>(3));
  HS_EXPECT_EQ(m.face_counts[1], static_cast<uint8_t>(2));
  HS_EXPECT_EQ(m.face_counts[2], static_cast<uint8_t>(2));
}

inline void test_reject_vertices_not_multiple_of_3() {
  static uint8_t buf[1 << 16];
  Arena arena(buf, sizeof(buf));
  std::vector<float> v = {0, 0, 0, 1}; // length 4, not a multiple of 3
  std::vector<int> f = {};

  PolyMesh m;
  MeshBuildStatus st = build_mesh_from_flat(v, f, arena, m);
  HS_EXPECT(st == MeshBuildStatus::kVerticesNotMultipleOf3,
            "non-multiple-of-3 vertex stream rejected");
  HS_EXPECT_EQ(arena.get_offset(), static_cast<size_t>(0)); // arena untouched
}

inline void test_reject_face_index_out_of_range() {
  static uint8_t buf[1 << 16];
  Arena arena(buf, sizeof(buf));
  std::vector<float> v = {0, 0, 0, 1, 0, 0}; // 2 verts => valid indices {0,1}

  // Index >= num_verts.
  {
    std::vector<int> f = {0, 1, 2, -1};
    PolyMesh m;
    HS_EXPECT(build_mesh_from_flat(v, f, arena, m) ==
                  MeshBuildStatus::kFaceIndexOutOfRange,
              "face index >= num_verts rejected");
    HS_EXPECT_EQ(arena.get_offset(), static_cast<size_t>(0));
  }
  // Index < -1 (a stray negative that is not the delimiter).
  {
    std::vector<int> f = {0, 1, -2, -1};
    PolyMesh m;
    HS_EXPECT(build_mesh_from_flat(v, f, arena, m) ==
                  MeshBuildStatus::kFaceIndexOutOfRange,
              "face index < -1 rejected");
    HS_EXPECT_EQ(arena.get_offset(), static_cast<size_t>(0));
  }
}

inline void test_face_vertex_count_boundary() {
  static uint8_t buf[1 << 16];
  Arena arena(buf, sizeof(buf));
  std::vector<float> v = {0, 0, 0}; // one vertex; index 0 is the only valid one

  // Exactly UINT8_MAX (255) verts in one face: the largest a uint8_t face_count
  // can hold — must build.
  {
    std::vector<int> f(UINT8_MAX, 0);
    f.push_back(-1);
    PolyMesh m;
    HS_EXPECT(build_mesh_from_flat(v, f, arena, m) == MeshBuildStatus::kOk,
              "255-vertex face builds (boundary)");
    HS_EXPECT_EQ(m.face_counts.size(), static_cast<size_t>(1));
    HS_EXPECT_EQ(m.face_counts[0], static_cast<uint8_t>(UINT8_MAX));
  }
  // One past the boundary (256) would wrap the uint8_t face_count: reject. Use a
  // fresh arena so the "no allocation on reject" check is unambiguous (the 255
  // case above already populated `arena`).
  {
    static uint8_t reject_buf[1 << 16];
    Arena reject_arena(reject_buf, sizeof(reject_buf));
    std::vector<int> f(UINT8_MAX + 1, 0);
    f.push_back(-1);
    PolyMesh m;
    HS_EXPECT(build_mesh_from_flat(v, f, reject_arena, m) ==
                  MeshBuildStatus::kFaceTooManyVerts,
              "256-vertex face rejected");
    HS_EXPECT_EQ(reject_arena.get_offset(), static_cast<size_t>(0));
  }
}

// The half-edge builder indexes everything in uint16_t, so a mesh whose vertex
// count OR half-edge (face-vertex) count exceeds UINT16_MAX must be rejected
// before the (uint16_t) store silently truncates. This bound is checked ahead of
// the arena budget, so a modest arena still reaches it without allocating.
inline void test_reject_index_range_overflow() {
  static uint8_t buf[1 << 16];

  // (a) vertex count > UINT16_MAX: 65536 verts, no faces.
  {
    Arena arena(buf, sizeof(buf));
    std::vector<float> v(static_cast<size_t>(UINT16_MAX + 1) * 3, 0.0f);
    std::vector<int> f; // no faces
    PolyMesh m;
    HS_EXPECT(build_mesh_from_flat(v, f, arena, m) ==
                  MeshBuildStatus::kIndexRangeOverflow,
              ">65535 vertices rejected");
    HS_EXPECT_EQ(arena.get_offset(), static_cast<size_t>(0));
  }

  // (b) half-edge count > UINT16_MAX while every single face stays within the
  //     uint8_t face-size limit: 258 faces of 255 verts (65790 face-verts).
  {
    Arena arena(buf, sizeof(buf));
    std::vector<float> v = {0, 0, 0, 1, 0, 0}; // 2 verts; index 0 is valid
    std::vector<int> f;
    for (int face = 0; face < 258; ++face) {
      for (int k = 0; k < UINT8_MAX; ++k)
        f.push_back(0);
      f.push_back(-1);
    }
    PolyMesh m;
    HS_EXPECT(build_mesh_from_flat(v, f, arena, m) ==
                  MeshBuildStatus::kIndexRangeOverflow,
              ">65535 face-vertex indices rejected");
    HS_EXPECT_EQ(arena.get_offset(), static_cast<size_t>(0));
  }
}

inline void test_reject_arena_overflow() {
  // A tiny arena that cannot fit even a single-vertex mesh's budget. The budget
  // check must reject BEFORE any allocation, so the arena trap never fires.
  static uint8_t buf[16];
  Arena arena(buf, sizeof(buf));
  std::vector<float> v = {0, 0, 0};
  std::vector<int> f = {0, -1};

  PolyMesh m;
  MeshBuildStatus st = build_mesh_from_flat(v, f, arena, m);
  HS_EXPECT(st == MeshBuildStatus::kArenaOverflow,
            "oversize mesh soft-rejected (no trap)");
  HS_EXPECT_EQ(arena.get_offset(), static_cast<size_t>(0)); // never allocated
}

inline int run_mesh_marshal_tests() {
  auto scope = hs_test::begin_module("mesh_marshal");
  test_valid_mesh_builds();
  test_counting_edge_cases();
  test_reject_vertices_not_multiple_of_3();
  test_reject_face_index_out_of_range();
  test_face_vertex_count_boundary();
  test_reject_index_range_overflow();
  test_reject_arena_overflow();
  return hs_test::end_module(scope);
}

} // namespace mesh_marshal_tests
} // namespace hs_test
