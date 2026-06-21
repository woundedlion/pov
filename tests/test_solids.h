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
 *   - Bounds: get_entry() out-of-range and get_by_name() unknown name TRAP
 *     (fail-fast), so only the valid boundary (last index) is exercised here.
 *   - Determinism: building the same registry entry twice yields identical
 *     vertex counts and positions.
 */
#pragma once

#include <cmath>
#include <cstdint>
#include "core/mesh.h"
#include "core/solids.h"
#include "tests/mesh_test_util.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace solids_tests {

// Generation budget. The complex Islamic-pattern generators chain many Conway
// operators (relax/ambo/bevel/hankin/...) and the WASM tooling path uses
// 4 MB scratch arenas for exactly these, resetting them before each build. We
// mirror that: two large scratch arenas (reset per solid) plus a geometry arena
// that holds the finalized result. A second geometry arena lets us build the
// same solid twice for the determinism check without aliasing storage.
inline uint8_t solids_geom_a[4 * 1024 * 1024];
inline uint8_t solids_geom_b[4 * 1024 * 1024];
inline uint8_t solids_scratch_a[4 * 1024 * 1024];
inline uint8_t solids_scratch_b[4 * 1024 * 1024];

// ---------------------------------------------------------------------------
// Structural invariants. check_face_counts_consistent(), check_indices_in_range(),
// and the unit-sphere check (check_all_unit_vertices) live in
// tests/mesh_test_util.h; the rest are specific to the registry path here.
// ---------------------------------------------------------------------------

/**
 * @brief Asserts every vertex coordinate is finite (no NaN/Inf from the generators).
 * @param m Mesh whose vertex coordinates are checked.
 */
inline void check_all_finite(const PolyMesh &m) {
  for (size_t i = 0; i < m.vertices.size(); ++i) {
    const Vector &v = m.vertices[i];
    HS_EXPECT_TRUE(std::isfinite(v.x) && std::isfinite(v.y) &&
                   std::isfinite(v.z));
  }
}

/**
 * @brief Asserts the mesh has at least one vertex, face count, and face index.
 * @param m Mesh to check for non-emptiness.
 */
inline void check_nonempty(const PolyMesh &m) {
  HS_EXPECT_TRUE(m.vertices.size() > 0);
  HS_EXPECT_TRUE(m.face_counts.size() > 0);
  HS_EXPECT_TRUE(m.faces.size() > 0);
}

/**
 * @brief Structural validity bundle applied to every registered solid.
 * @param m Mesh to validate.
 * @details Asserts non-empty, consistent face_counts/faces, in-range face
 *          indices, and finite coordinates.
 */
inline void check_basic(const PolyMesh &m) {
  check_nonempty(m);
  check_face_counts_consistent(m);
  check_indices_in_range(m);
  check_all_finite(m);
}

/**
 * @brief Builds a registry entry by index, finalizing into the supplied geometry arena.
 * @param index Registry entry index to build.
 * @param geom Geometry arena that holds the finalized mesh.
 * @return The finalized PolyMesh for the entry.
 * @details Uses fresh scratch arenas (reset each call) for the generation pass.
 */
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

/**
 * @brief Verifies every simple (Platonic/Archimedean) entry builds to a valid
 *        mesh whose vertices sit on the unit sphere within 1e-3.
 */
inline void test_simple_registry_solids_are_spherical_and_valid() {
  for (size_t i = 0; i < Solids::Collections::get_simple_solids().size(); ++i) {
    Arena geom(solids_geom_a, sizeof(solids_geom_a));
    PolyMesh m = build_index(i, geom);
    check_basic(m);
    check_all_unit_vertices(m, 1e-3f); // Platonic/Archimedean seeds + ops stay on unit sphere
  }
}

/**
 * @brief Verifies every Catalan entry builds to a valid mesh whose vertices sit
 *        on the unit sphere within 1e-3.
 * @details Catalan indices follow the simple-solid block in the registry.
 */
inline void test_catalan_registry_solids_are_spherical_and_valid() {
  const size_t base = Solids::Collections::get_simple_solids().size();
  for (size_t k = 0; k < Solids::Collections::get_catalan_solids().size();
       ++k) {
    Arena geom(solids_geom_a, sizeof(solids_geom_a));
    PolyMesh m = build_index(base + k, geom);
    check_basic(m);
    check_all_unit_vertices(m, 1e-3f); // Catalan = dual of an Archimedean; ops re-normalize
  }
}

// ---------------------------------------------------------------------------
// Registry integrity — Islamic patterns. Structural invariants only:
// hankin/expand on a pattern can legitimately move points off the unit sphere
// and may yield open meshes, so magnitude/closure are NOT asserted here.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies every Islamic-pattern entry builds to a structurally valid mesh.
 * @details No sphere assertion (hankin/expand may move points off the sphere).
 *          Islamic indices follow the simple + Catalan blocks in the registry.
 */
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

/**
 * @brief Verifies NUM_ENTRIES equals the sum of the three registries.
 * @details Confirms indexing the full range corresponds to the combined
 *          simple + Catalan + Islamic collection sizes.
 */
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

/**
 * @brief Builds the registry entry at `index`, derives its half-edge mesh, and
 *        asserts the Euler characteristic V - E + F == 2.
 * @param index Registry entry index to check.
 * @details Edges are counted as half_edges/2, matching test_mesh.h.
 */
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

/**
 * @brief Verifies V - E + F == 2 for each Platonic solid (closed manifolds).
 */
inline void test_euler_platonic_solids() {
  // simple_registry indices 0-4 are the Platonic solids (closed manifolds).
  for (size_t i = 0; i < Solids::Collections::get_platonic_solids().size(); ++i)
    check_euler_for_index(i);
}

// ---------------------------------------------------------------------------
// Fallbacks (read directly from solids.h).
// ---------------------------------------------------------------------------

/**
 * @brief Verifies the last valid registry index builds correctly (range boundary).
 * @details Out-of-range get_entry() and unknown get_by_name() TRAP (fail-fast),
 *          so those error paths can't be exercised without death-test
 *          infrastructure; only the valid boundary is checked here.
 */
inline void test_get_entry_last_valid_index_builds() {
  const Solids::Entry &e = Solids::get_entry(Solids::NUM_ENTRIES - 1);
  HS_EXPECT_TRUE(e.name != nullptr);

  Arena geom(solids_geom_a, sizeof(solids_geom_a));
  Arena a(solids_scratch_a, sizeof(solids_scratch_a));
  Arena b(solids_scratch_b, sizeof(solids_scratch_b));
  PolyMesh m = Solids::finalize_solid(e.generate(a, b), geom);
  check_basic(m);
  check_all_unit_vertices(m, 1e-3f);
}

/**
 * @brief Verifies get_by_name("octahedron") returns that specific solid.
 * @details Asserts the result is valid, on the unit sphere, and has the
 *          octahedron's 6 vertices and 8 faces.
 */
inline void test_get_by_name_known_returns_that_solid() {
  Arena geom(solids_geom_a, sizeof(solids_geom_a));
  Arena a(solids_scratch_a, sizeof(solids_scratch_a));
  Arena b(solids_scratch_b, sizeof(solids_scratch_b));
  PolyMesh m = Solids::get_by_name(geom, a, b, "octahedron");
  check_basic(m);
  check_all_unit_vertices(m, 1e-3f);
  HS_EXPECT_EQ(m.vertices.size(), (size_t)6);
  HS_EXPECT_EQ(m.face_counts.size(), (size_t)8);
}

// ---------------------------------------------------------------------------
// Determinism: building the same entry twice yields identical geometry. We use
// two distinct geometry arenas so the two results coexist for comparison.
// ---------------------------------------------------------------------------

/**
 * @brief Builds the entry at `index` twice into separate arenas and asserts the
 *        two meshes match in counts, vertex positions (1e-6), and face indices.
 * @param index Registry entry index to build twice.
 */
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

/**
 * @brief Verifies determinism on a hardcoded solid (the cube).
 * @details The cube is pure data with no procedural ops, so two builds must be
 *          bit-identical.
 */
inline void test_determinism_hardcoded_platonic() {
  // index 1 = cube (pure data, no procedural ops).
  check_determinism_for_index(1);
}

/**
 * @brief Verifies determinism through a Conway op pipeline (cube -> ambo).
 * @details The procedural path must reproduce identical geometry across builds.
 */
inline void test_determinism_archimedean_with_conway_ops() {
  // index 6 = cuboctahedron (cube -> ambo): exercises a Conway op pipeline.
  check_determinism_for_index(6);
}

/**
 * @brief Verifies determinism on the longest path (a deeply chained
 *        Islamic-pattern generator), which must also be reproducible.
 */
inline void test_determinism_complex_islamic() {
  const size_t base = Solids::Collections::get_simple_solids().size() +
                      Solids::Collections::get_catalan_solids().size();
  check_determinism_for_index(base); // first islamic entry
}

// ---------------------------------------------------------------------------
// High-water regression at the real shipping arena configuration.
//
// The IslamicStars effect is the one that ships these recipes: spawn_shape()
// builds get_islamic_solids()[idx].generate(a, b) through scratch_arena_a /
// scratch_arena_b, which init() sizes via configure_arenas(..., 120 KB, 120 KB)
// (effects/IslamicStars.h). SolidBuilder ping-pongs the two arenas WITHOUT
// resetting between ops (see the SCRATCH ARENA CONTRACT in core/conway.h), so a
// whole recipe chain accumulates into that pair and its peak is the recipe's
// high-water mark. The COMPOSITION POLARITY note in core/conway.h calls out
// cube_relax_bevel33_relax_hk68_expand5 as the one recipe whose relax()-after-
// bevel() runs input and output on the same arena; that "measured to fit" was a
// comment, and this test makes it an automated guard. A recipe edit or a Conway
// operator-table change that pushes peak scratch over budget would otherwise
// surface only as a device-only OOM trap, invisible to the host suite.
//
// 64-bit host vs 32-bit device: a recipe's scratch is flat POD arrays (PolyMesh
// is ArenaVector<Vector/uint8_t/uint16_t/int>, the half-edge build uses POD
// records) whose element sizes are identical on both builds, so these arena
// footprints do NOT carry the host/device pointer delta core/memory.h warns
// about for pointer-bearing pooled structs. Where any delta exists at all it can
// only make the 64-bit host figure LARGER (8-byte vs 4-byte pointers), never
// smaller, so a host high-water mark is a conservative UPPER bound on the device
// figure: passing here guarantees the recipe fits the device's 120 KB split.
//
// The backing buffers (solids_scratch_a/b, multi-MB) are deliberately larger
// than the asserted budget so a regression yields a precise high-water assertion
// — the actual byte count vs 120 KB — instead of an opaque mid-build OOM trap;
// the HS_EXPECT_LE against the real 120 KB budget is the guard.
// ---------------------------------------------------------------------------

constexpr size_t kIslamicScratchBudget =
    120 * 1024; /**< IslamicStars' per-arena scratch split (symmetric). */

/**
 * @brief Runs one Islamic recipe through a real-budget arena pair and asserts
 *        each arena's peak usage stays within IslamicStars' 120 KB split.
 * @param entry Registry entry whose generator is exercised.
 * @details A fresh arena pair per recipe isolates each measurement; the recipe
 *          builds through (a, b) exactly as IslamicStars::spawn_shape does.
 */
inline void check_high_water_for_recipe(const Solids::Entry &entry) {
  Arena a(solids_scratch_a, sizeof(solids_scratch_a));
  Arena b(solids_scratch_b, sizeof(solids_scratch_b));
  PolyMesh m = entry.generate(a, b);
  check_nonempty(m); // the build actually produced geometry through this pair

  HS_EXPECT_LE(a.get_high_water_mark(), kIslamicScratchBudget);
  HS_EXPECT_LE(b.get_high_water_mark(), kIslamicScratchBudget);
}

/**
 * @brief Verifies every Islamic-pattern recipe fits IslamicStars' 120 KB scratch
 *        split — the configuration these recipes actually ship through.
 * @details cube_relax_bevel33_relax_hk68_expand5 (the recipe the polarity note
 *          flags as the one running an op on a non-alternating arena) is the
 *          first Islamic entry and so is covered by this sweep.
 */
inline void test_islamic_recipes_fit_islamicstars_budget() {
  for (const Solids::Entry &e : Solids::islamic_registry)
    check_high_water_for_recipe(e);
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

/**
 * @brief Runs all solids tests.
 * @return The module's failure count.
 */
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

  test_islamic_recipes_fit_islamicstars_budget();

  return hs_test::end_module(scope);
}

} // namespace solids_tests
} // namespace hs_test

