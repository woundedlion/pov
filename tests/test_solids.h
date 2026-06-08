/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/solids.h.
 *
 * Coverage:
 *   - Registry integrity: every registered solid (simple + catalan + islamic)
 *     builds into a non-empty mesh with finite vertices, in-range face indices,
 *     and consistent face_counts/faces totals.
 *   - Unit-sphere intent: Platonic/Archimedean/Catalan generators are designed
 *     to live on the unit sphere (seeds are unit vectors and the Conway ops
 *     used here re-normalize), so their vertices are asserted unit-magnitude.
 *     Islamic-pattern seeds may be open / non-spherical (hankin/expand on a
 *     pattern can move points off the sphere), so magnitude is NOT asserted for
 *     those — only finiteness + structural invariants.
 *   - Euler characteristic V - E + F == 2 for the hardcoded closed Platonic
 *     solids, using the half-edge edge count (E = half_edges/2) as in
 *     test_mesh.h.
 *   - Bounds: get_entry() out-of-range and get_by_name() unknown name now TRAP
 *     (fail-fast), so only the valid boundary (last index) is exercised here.
 *   - Determinism: building the same registry entry twice yields identical
 *     vertex counts and positions.
 */
#pragma once

#include <cmath>
#include <cstdint>
#include "core/mesh.h"
#include "core/solids.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace solids_tests {

// Generation budget. The complex Islamic-pattern generators chain many Conway
// operators (relax/ambo/bitruncate/hankin/...) and the WASM tooling path uses
// 4 MB scratch arenas for exactly these, resetting them before each build. We
// mirror that: two large scratch arenas (reset per solid) plus a geometry arena
// that holds the finalized result. A second geometry arena lets us build the
// same solid twice for the determinism check without aliasing storage.
inline uint8_t solids_geom_a[4 * 1024 * 1024];
inline uint8_t solids_geom_b[4 * 1024 * 1024];
inline uint8_t solids_scratch_a[4 * 1024 * 1024];
inline uint8_t solids_scratch_b[4 * 1024 * 1024];

// ---------------------------------------------------------------------------
// Structural invariants (mirrors the helpers in test_conway.h, but operate on
// the finalized PolyMesh returned by the registry).
// ---------------------------------------------------------------------------

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

inline void check_all_finite(const PolyMesh &m) {
  for (size_t i = 0; i < m.vertices.size(); ++i) {
    const Vector &v = m.vertices[i];
    HS_EXPECT_TRUE(std::isfinite(v.x) && std::isfinite(v.y) &&
                   std::isfinite(v.z));
  }
}

inline void check_all_unit(const PolyMesh &m) {
  for (size_t i = 0; i < m.vertices.size(); ++i)
    HS_EXPECT_NEAR(m.vertices[i].length(), 1.0f, 1e-3f);
}

inline void check_nonempty(const PolyMesh &m) {
  HS_EXPECT_TRUE(m.vertices.size() > 0);
  HS_EXPECT_TRUE(m.face_counts.size() > 0);
  HS_EXPECT_TRUE(m.faces.size() > 0);
}

inline void check_basic(const PolyMesh &m) {
  check_nonempty(m);
  check_face_counts_consistent(m);
  check_indices_in_range(m);
  check_all_finite(m);
}

// Build a registry entry by index into fresh scratch arenas (reset each call),
// finalizing into the supplied geometry arena.
inline PolyMesh build_index(size_t index, Arena &geom) {
  Arena a(solids_scratch_a, sizeof(solids_scratch_a));
  Arena b(solids_scratch_b, sizeof(solids_scratch_b));
  return Solids::finalize_solid(Solids::get_entry(index).generate(a, b), geom);
}

// ---------------------------------------------------------------------------
// Registry integrity — spherical families (Platonic, Archimedean, Catalan).
// These are designed to live on the unit sphere, so we additionally assert
// unit magnitude.
// ---------------------------------------------------------------------------

inline void test_simple_registry_solids_are_spherical_and_valid() {
  for (size_t i = 0; i < Solids::Collections::get_simple_solids().size(); ++i) {
    Arena geom(solids_geom_a, sizeof(solids_geom_a));
    PolyMesh m = build_index(i, geom);
    check_basic(m);
    check_all_unit(m); // Platonic/Archimedean seeds + ops stay on unit sphere
  }
}

inline void test_catalan_registry_solids_are_spherical_and_valid() {
  const size_t base = Solids::Collections::get_simple_solids().size();
  for (size_t k = 0; k < Solids::Collections::get_catalan_solids().size();
       ++k) {
    Arena geom(solids_geom_a, sizeof(solids_geom_a));
    PolyMesh m = build_index(base + k, geom);
    check_basic(m);
    check_all_unit(m); // Catalan = dual of an Archimedean; ops re-normalize
  }
}

// ---------------------------------------------------------------------------
// Registry integrity — Islamic patterns. Structural invariants only:
// hankin/expand on a pattern can legitimately move points off the unit sphere
// and may yield open meshes, so magnitude/closure are NOT asserted here.
// ---------------------------------------------------------------------------

inline void test_islamic_registry_solids_are_valid() {
  const size_t base = Solids::Collections::get_simple_solids().size() +
                      Solids::Collections::get_catalan_solids().size();
  for (size_t k = 0; k < Solids::Collections::get_islamic_solids().size();
       ++k) {
    Arena geom(solids_geom_a, sizeof(solids_geom_a));
    PolyMesh m = build_index(base + k, geom);
    check_basic(m);
  }
}

// Sanity: NUM_ENTRIES equals the sum of the three registries, and indexing
// the full range never crashes.
inline void test_registry_count_matches_collections() {
  size_t sum = Solids::Collections::get_simple_solids().size() +
               Solids::Collections::get_catalan_solids().size() +
               Solids::Collections::get_islamic_solids().size();
  HS_EXPECT_EQ(sum, (size_t)Solids::NUM_ENTRIES);
}

// ---------------------------------------------------------------------------
// Euler characteristic V - E + F == 2 for the closed hardcoded Platonic
// solids. We build the half-edge mesh and count edges as half_edges/2, exactly
// as test_mesh.h does.
// ---------------------------------------------------------------------------

inline void check_euler_for_index(size_t index) {
  Arena geom(solids_geom_a, sizeof(solids_geom_a));
  PolyMesh m = build_index(index, geom);

  // Half-edge construction needs its own scratch; reuse geom_b.
  Arena he_arena(solids_geom_b, sizeof(solids_geom_b));
  HalfEdgeMesh he(he_arena, m);

  int V = static_cast<int>(he.vertices.size());
  int E = static_cast<int>(he.half_edges.size()) / 2;
  int F = static_cast<int>(he.faces.size());
  HS_EXPECT_EQ(V - E + F, 2);
}

inline void test_euler_platonic_solids() {
  // simple_registry indices 0-4 are the Platonic solids (closed manifolds).
  for (size_t i = 0; i < Solids::Collections::get_platonic_solids().size(); ++i)
    check_euler_for_index(i);
}

// ---------------------------------------------------------------------------
// Fallbacks (read directly from solids.h).
// ---------------------------------------------------------------------------

// Out-of-range get_entry() and unknown get_by_name() now TRAP (fail-fast) rather
// than substituting a default solid, so those error paths can't be exercised
// here without death-test infrastructure. Instead verify the last valid index
// builds correctly (boundary of the valid range).
inline void test_get_entry_last_valid_index_builds() {
  const Solids::Entry &e = Solids::get_entry(Solids::NUM_ENTRIES - 1);
  HS_EXPECT_TRUE(e.name != nullptr);

  Arena geom(solids_geom_a, sizeof(solids_geom_a));
  Arena a(solids_scratch_a, sizeof(solids_scratch_a));
  Arena b(solids_scratch_b, sizeof(solids_scratch_b));
  PolyMesh m = Solids::finalize_solid(e.generate(a, b), geom);
  check_basic(m);
  check_all_unit(m);
}

inline void test_get_by_name_known_returns_that_solid() {
  Arena geom(solids_geom_a, sizeof(solids_geom_a));
  Arena a(solids_scratch_a, sizeof(solids_scratch_a));
  Arena b(solids_scratch_b, sizeof(solids_scratch_b));
  PolyMesh m = Solids::get_by_name(geom, a, b, "octahedron");
  check_basic(m);
  check_all_unit(m);
  HS_EXPECT_EQ(m.vertices.size(), (size_t)6);
  HS_EXPECT_EQ(m.face_counts.size(), (size_t)8);
}

// ---------------------------------------------------------------------------
// Determinism: building the same entry twice yields identical geometry. We use
// two distinct geometry arenas so the two results coexist for comparison.
// ---------------------------------------------------------------------------

inline void check_determinism_for_index(size_t index) {
  Arena geom1(solids_geom_a, sizeof(solids_geom_a));
  PolyMesh m1 = build_index(index, geom1);

  Arena geom2(solids_geom_b, sizeof(solids_geom_b));
  PolyMesh m2 = build_index(index, geom2);

  HS_EXPECT_EQ(m1.vertices.size(), m2.vertices.size());
  HS_EXPECT_EQ(m1.face_counts.size(), m2.face_counts.size());
  HS_EXPECT_EQ(m1.faces.size(), m2.faces.size());

  if (m1.vertices.size() == m2.vertices.size()) {
    for (size_t i = 0; i < m1.vertices.size(); ++i) {
      HS_EXPECT_NEAR(m1.vertices[i].x, m2.vertices[i].x, 1e-6f);
      HS_EXPECT_NEAR(m1.vertices[i].y, m2.vertices[i].y, 1e-6f);
      HS_EXPECT_NEAR(m1.vertices[i].z, m2.vertices[i].z, 1e-6f);
    }
  }
  if (m1.faces.size() == m2.faces.size()) {
    for (size_t i = 0; i < m1.faces.size(); ++i)
      HS_EXPECT_EQ(m1.faces[i], m2.faces[i]);
  }
}

inline void test_determinism_hardcoded_platonic() {
  // index 1 = cube (pure data, no procedural ops).
  check_determinism_for_index(1);
}

inline void test_determinism_archimedean_with_conway_ops() {
  // index 6 = cuboctahedron (cube -> ambo): exercises a Conway op pipeline.
  check_determinism_for_index(6);
}

inline void test_determinism_complex_islamic() {
  // A long chained Islamic-pattern generator must also be reproducible.
  const size_t base = Solids::Collections::get_simple_solids().size() +
                      Solids::Collections::get_catalan_solids().size();
  check_determinism_for_index(base); // first islamic entry
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

inline int run_solids_tests() {
  auto scope = hs_test::begin_module("solids");

  test_registry_count_matches_collections();

  test_simple_registry_solids_are_spherical_and_valid();
  test_catalan_registry_solids_are_spherical_and_valid();
  test_islamic_registry_solids_are_valid();

  test_euler_platonic_solids();

  test_get_entry_last_valid_index_builds();
  test_get_by_name_known_returns_that_solid();

  test_determinism_hardcoded_platonic();
  test_determinism_archimedean_with_conway_ops();
  test_determinism_complex_islamic();

  return hs_test::end_module(scope);
}

} // namespace solids_tests
} // namespace hs_test

