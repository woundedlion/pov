/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/mesh/solids.h.
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
 *   - Sliver-face invariant: every Islamic recipe keeps its longest geodesic
 *     edge within 6x the median edge.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include "core/mesh/mesh.h"
#include "core/color/palettes.h"
#include "core/mesh/solids.h"
#include "tests/mesh_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace solids_tests {

// Two scratch arenas (reset per solid) plus two geometry arenas, sized to match
// the WASM tooling path's 4 MB scratch. The second geometry arena lets the
// determinism check build the same solid twice without aliasing storage.
inline uint8_t solids_geom_a[4 * 1024 * 1024];
inline uint8_t solids_geom_b[4 * 1024 * 1024];
inline uint8_t solids_scratch_a[4 * 1024 * 1024];
inline uint8_t solids_scratch_b[4 * 1024 * 1024];

// ---------------------------------------------------------------------------
// Structural invariants. check_face_counts_consistent(),
// check_indices_in_range(), and the unit-sphere check (check_all_unit_vertices)
// live in tests/mesh_test_util.h; the rest are specific to the registry path
// here.
// ---------------------------------------------------------------------------

/**
 * @brief Asserts every vertex coordinate is finite (no NaN/Inf from the
 * generators).
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
 * @brief Builds a registry entry by index, finalizing into the supplied
 * geometry arena.
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
    check_all_unit_vertices(m, 1e-3f);
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
    check_all_unit_vertices(m, 1e-3f);
  }
}

// ---------------------------------------------------------------------------
// Registry integrity — Islamic patterns. Structural invariants only:
// hankin/expand on a pattern can legitimately move points off the unit sphere
// and may yield open meshes, so magnitude/closure are NOT asserted here.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies every Islamic-pattern entry builds to a structurally valid
 * mesh.
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
 * @brief Verifies no Islamic-pattern solid has sliver faces: the longest
 *        geodesic edge stays within 6x the median edge.
 * @details A hankin contact angle near a resonance (contact planes of one
 *          corner class near-parallel) slings star points far from their
 *          corners, producing sliver faces that render as long lines. Healthy
 *          registry recipes measure at most ~3.4x; the broken hk43 recipe
 *          measured 23.8x.
 */
inline void test_islamic_solids_have_no_sliver_edges() {
  const size_t base = Solids::Collections::get_simple_solids().size() +
                      Solids::Collections::get_catalan_solids().size();
  for (size_t k = 0; k < Solids::Collections::get_islamic_solids().size();
       ++k) {
    Arena geom(solids_geom_a, sizeof(solids_geom_a));
    PolyMesh m = build_index(base + k, geom);
    std::vector<float> edges;
    size_t off = 0;
    for (size_t f = 0; f < m.face_counts.size(); ++f) {
      int n = m.face_counts[f];
      for (int i = 0; i < n; ++i) {
        Vector u = m.vertices[m.faces[off + i]].normalized();
        Vector v = m.vertices[m.faces[off + (i + 1) % n]].normalized();
        edges.push_back(std::acos(std::max(-1.0f, std::min(1.0f, dot(u, v)))));
      }
      off += n;
    }
    std::sort(edges.begin(), edges.end());
    float median = edges[edges.size() / 2];
    float max = edges.back();
    HS_EXPECT_TRUE(max <= 6.0f * median);
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
 * @param expected_V,expected_E,expected_F Optional independent count oracles;
 *        pass -1 (the default) to skip a given check. When provided, each is
 *        asserted exactly. V - E + F == 2 alone is blind to a class-wide
 *        omission/duplication that preserves the relation (e.g. every face
 *        split in two), so a fixed expected-count oracle is the real guard.
 * @details Edges are counted as half_edges/2, matching test_mesh.h.
 */
inline void check_euler_for_index(size_t index, int expected_V = -1,
                                  int expected_E = -1, int expected_F = -1) {
  Arena geom(solids_geom_a, sizeof(solids_geom_a));
  PolyMesh m = build_index(index, geom);

  // Half-edge construction needs its own scratch; reuse geom_b.
  Arena he_arena(solids_geom_b, sizeof(solids_geom_b));
  HalfEdgeMesh he(he_arena, m);

  int V = static_cast<int>(m.vertices.size());
  int E = static_cast<int>(he.half_edges.size()) / 2;
  int F = static_cast<int>(he.faces.size());
  HS_EXPECT_EQ(V - E + F, 2);

  if (expected_V >= 0)
    HS_EXPECT_EQ(V, expected_V);
  if (expected_E >= 0)
    HS_EXPECT_EQ(E, expected_E);
  if (expected_F >= 0)
    HS_EXPECT_EQ(F, expected_F);
}

/**
 * @brief Verifies the Euler characteristic and the exact (V, E, F) counts for
 *        each Platonic solid (closed manifolds).
 * @details The per-solid counts are fixed mathematical constants, so they form
 *          an independent oracle: a builder bug that uniformly mis-counts
 *          vertices or faces while still satisfying V - E + F == 2 is caught
 *          here, where the bare relation check could not see it.
 */
inline void test_euler_platonic_solids() {
  // simple_registry indices 0-4 are the Platonic solids, in this fixed order.
  struct Counts {
    int v, e, f;
  };
  static constexpr Counts PLATONIC[] = {
      {4, 6, 4},    // tetrahedron
      {8, 12, 6},   // cube
      {6, 12, 8},   // octahedron
      {20, 30, 12}, // dodecahedron
      {12, 30, 20}, // icosahedron
  };
  const auto platonic = Solids::Collections::get_platonic_solids();
  HS_EXPECT_EQ(platonic.size(), sizeof(PLATONIC) / sizeof(PLATONIC[0]));
  for (size_t i = 0; i < platonic.size(); ++i)
    check_euler_for_index(i, PLATONIC[i].v, PLATONIC[i].e, PLATONIC[i].f);
}

/**
 * @brief Verifies every Archimedean and Catalan entry is a closed 2-manifold
 *        (V-E+F==2).
 * @details Extends the topological oracle over the two spherical families
 *          between the Platonic block and the Islamic block. Archimedean
 * indices follow the Platonic block inside the simple registry; Catalan indices
 *          follow the whole simple block. Exact per-entry counts are not pinned
 *          here — the Euler invariant catches a generator regression that opens
 * a seam, drops a face, or duplicates geometry.
 */
inline void test_euler_archimedean_catalan_solids() {
  const size_t archimedean_base =
      Solids::Collections::get_platonic_solids().size();
  for (size_t k = 0; k < Solids::Collections::get_archimedean_solids().size();
       ++k)
    check_euler_for_index(archimedean_base + k);

  const size_t catalan_base = Solids::Collections::get_simple_solids().size();
  for (size_t k = 0; k < Solids::Collections::get_catalan_solids().size(); ++k)
    check_euler_for_index(catalan_base + k);
}

/**
 * @brief Verifies every Islamic-pattern entry is a closed 2-manifold
 * (V-E+F==2).
 * @details Stronger than test_islamic_registry_solids_are_valid's check_basic
 *          (finite / consistent / in-range), which a wrong-but-self-consistent
 *          generator still passes. Despite the cautious "may yield open meshes"
 *          note on the structural test above, every entry currently in the
 *          registry closes (verified across all of them), so the Euler oracle
 * is enforceable and catches a generator regression that opens a seam, drops a
 * face, or duplicates geometry — the topological equivalent of the exact V/E/F
 * oracle the Platonic solids get. Exact per-entry counts are deliberately NOT
 * pinned: the pattern generators are actively tuned, so a golden count would
 * invert the signal (every intentional retune reds the test) the way
 * test_effects.h rejects golden-frame hashing; the Euler invariant is the
 * stable altitude. If a future entry is intentionally open, exclude it here
 * with a comment rather than weakening the check for all.
 */
inline void test_islamic_registry_solids_are_closed() {
  const size_t base = Solids::Collections::get_simple_solids().size() +
                      Solids::Collections::get_catalan_solids().size();
  for (size_t k = 0; k < Solids::Collections::get_islamic_solids().size(); ++k)
    check_euler_for_index(base + k);
}

// ---------------------------------------------------------------------------
// Fallbacks (read directly from solids.h).
// ---------------------------------------------------------------------------

/**
 * @brief Verifies the last valid registry index builds correctly (range
 * boundary).
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

/**
 * @brief Verifies registry names are globally unique and index/name lookups
 * agree.
 * @details The WASM picker enumerates solids by global index but builds them by
 *          first-name match, so a name duplicated across the three registries
 *          would make those two paths silently diverge. Assert every name is
 *          distinct and that find_entry(get_entry(i).name) resolves back to
 * that same entry.
 */
inline void test_registry_names_unique_and_roundtrip() {
  for (int i = 0; i < Solids::NUM_ENTRIES; ++i) {
    const Solids::Entry &ei = Solids::get_entry(i);
    HS_EXPECT_TRUE(ei.name != nullptr);
    HS_EXPECT_TRUE(Solids::find_entry(ei.name) == &ei);
    for (int j = i + 1; j < Solids::NUM_ENTRIES; ++j) {
      HS_EXPECT_TRUE(std::string_view(ei.name) !=
                     std::string_view(Solids::get_entry(j).name));
    }
  }
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

  for (size_t i = 0; i < m1.vertices.size(); ++i) {
    HS_EXPECT_NEAR(m1.vertices[i].x, m2.vertices[i].x, 1e-6f);
    HS_EXPECT_NEAR(m1.vertices[i].y, m2.vertices[i].y, 1e-6f);
    HS_EXPECT_NEAR(m1.vertices[i].z, m2.vertices[i].z, 1e-6f);
  }
  for (size_t i = 0; i < m1.faces.size(); ++i)
    HS_EXPECT_EQ(m1.faces[i], m2.faces[i]);
}

/**
 * @brief Verifies determinism on a hardcoded solid (the cube).
 * @details The cube is pure data with no procedural ops, so two builds must be
 *          bit-identical.
 */
inline void test_determinism_hardcoded_platonic() {
  // index 1 = cube
  check_determinism_for_index(1);
}

/**
 * @brief Verifies determinism through a Conway op pipeline (cube -> ambo).
 * @details The procedural path must reproduce identical geometry across builds.
 */
inline void test_determinism_archimedean_with_conway_ops() {
  // index 6 = cuboctahedron (cube -> ambo)
  check_determinism_for_index(6);
}

/**
 * @brief Verifies determinism across the whole Islamic-pattern family.
 * @details The Islamic generators are the deepest Conway-op chains — the most
 *          likely to introduce order/RNG-dependent nondeterminism — and a
 *          nondeterministic op might be reached only by a later entry, not the
 *          first. So this re-builds and diffs every islamic index (mirroring
 * the registry-integrity loops above), not just the first.
 */
inline void test_determinism_complex_islamic() {
  const size_t base = Solids::Collections::get_simple_solids().size() +
                      Solids::Collections::get_catalan_solids().size();
  for (size_t k = 0; k < Solids::Collections::get_islamic_solids().size(); ++k)
    check_determinism_for_index(base + k);
}

// ---------------------------------------------------------------------------
// High-water regression at the real shipping arena configuration.
//
// IslamicStars::spawn_shape builds each recipe through a 114 KB / 80 KB scratch
// pair that ping-pongs WITHOUT resetting between ops, so a recipe chain's peak
// is its high-water mark. Over-budget would otherwise surface only as a
// device-only OOM trap. Scratch is flat POD whose only host/device delta
// (64-bit pointers) can only make the host figure larger, so the host
// high-water mark is a conservative upper bound on the device figure.
// ---------------------------------------------------------------------------

constexpr size_t ISLAMIC_SCRATCH_A_BUDGET =
    114 * 1024; /**< IslamicStars' scratch_a split (mirrors init()). */
constexpr size_t ISLAMIC_SCRATCH_B_BUDGET =
    80 * 1024; /**< IslamicStars' scratch_b split (mirrors init()). */

/**
 * @brief Runs one Islamic recipe through a real-budget arena pair and asserts
 *        each arena's peak usage stays within IslamicStars' scratch split.
 * @param entry Registry entry whose generator is exercised.
 * @details A fresh arena pair per recipe isolates each measurement; the recipe
 *          builds through (a, b) exactly as IslamicStars::spawn_shape does.
 */
inline void check_high_water_for_recipe(const Solids::Entry &entry) {
  Arena a(solids_scratch_a, sizeof(solids_scratch_a));
  Arena b(solids_scratch_b, sizeof(solids_scratch_b));
  PolyMesh m = entry.generate(a, b);
  check_nonempty(m);

  HS_EXPECT_LE(a.get_high_water_mark(), ISLAMIC_SCRATCH_A_BUDGET);
  HS_EXPECT_LE(b.get_high_water_mark(), ISLAMIC_SCRATCH_B_BUDGET);
}

/**
 * @brief Verifies every Islamic-pattern recipe fits IslamicStars' scratch
 *        split — the configuration these recipes actually ship through.
 * @details cube_relax_bevel33_relax_hk675_expand5 (the recipe the polarity note
 *          flags as the one running an op on a non-alternating arena) is the
 *          first Islamic entry and so is covered by this sweep.
 */
inline void test_islamic_recipes_fit_islamicstars_budget() {
  for (const Solids::Entry &e : Solids::islamic_registry)
    check_high_water_for_recipe(e);
}

// ---------------------------------------------------------------------------
// Persistent-budget regression for the IslamicStars carousel.
//
// Guards the persistent half of IslamicStars' split (device GLOBAL_ARENA_SIZE
// minus the scratch pools): the baked palette bank plus the double-buffered
// carousel. The native 8 MB GLOBAL_ARENA_SIZE means a device persistent
// overflow can't surface by running the effect here. Peak residents during a
// cross-fade are the palette bank plus the two adjacent carousel slots that
// coexist until the swap (spawn_shape cycles idx % N), so the peak is the
// largest adjacent-pair sum, not twice the largest single slot.
// ---------------------------------------------------------------------------

constexpr size_t ISLAMIC_PERSISTENT_BUDGET =
    DEVICE_GLOBAL_ARENA_SIZE - ISLAMIC_SCRATCH_A_BUDGET -
    ISLAMIC_SCRATCH_B_BUDGET; /**< IslamicStars' persistent split. */

/**
 * @brief Verifies the worst adjacent pair of Islamic shapes, plus the palette
 *        bank, fits IslamicStars' persistent carousel split.
 * @details Compiles + classifies every Islamic solid into a fresh arena exactly
 *          as spawn_shape does, records each slot's footprint, then asserts the
 *          largest registry-adjacent pair (the carousel's worst cross-fade
 *          coexistence) plus one palette bank stays within the device budget,
 *          and that the largest single slot fits through scratch_b (the
 *          compact_keep_front evacuation path).
 */
inline void test_islamic_solids_fit_islamicstars_persistent_budget() {
  size_t palette_bytes;
  {
    Arena pal(solids_geom_b, sizeof(solids_geom_b));
    MeshPaletteBank bank;
    bank.bake_all(pal);
    palette_bytes = pal.get_high_water_mark();
  }

  const size_t base = Solids::Collections::get_simple_solids().size() +
                      Solids::Collections::get_catalan_solids().size();
  const auto islamic = Solids::Collections::get_islamic_solids();
  const size_t N = islamic.size();

  size_t slot_bytes[Solids::NUM_ENTRIES];
  size_t worst_slot = 0, worst_v = 0, worst_f = 0, worst_k = 0;
  for (size_t k = 0; k < N; ++k) {
    Arena geom(solids_geom_a, sizeof(solids_geom_a));
    PolyMesh raw = build_index(base + k, geom);

    // The three arenas are distinct backing buffers (no aliasing with geom,
    // which still holds raw).
    Arena slot_arena(solids_scratch_a, sizeof(solids_scratch_a));
    Arena sa(solids_scratch_b, sizeof(solids_scratch_b));
    Arena sb(solids_geom_b, sizeof(solids_geom_b));
    MeshState slot;
    MeshOps::compile(raw, slot, slot_arena, scratch_arena_a);
    MeshOps::classify_faces_by_topology(slot, sa, sb, slot_arena);

    slot_bytes[k] = slot_arena.get_high_water_mark();
    if (slot_bytes[k] > worst_slot) {
      worst_slot = slot_bytes[k];
      worst_v = slot.vertices.size();
      worst_f = slot.face_counts.size();
      worst_k = k;
    }
  }

  // Largest registry-adjacent pair, including the (N-1, 0) cycle wrap.
  size_t worst_pair = 0, worst_pair_i = 0;
  for (size_t i = 0; i < N; ++i) {
    size_t pair = slot_bytes[i] + slot_bytes[(i + 1) % N];
    if (pair > worst_pair) {
      worst_pair = pair;
      worst_pair_i = i;
    }
  }

  const size_t peak = palette_bytes + worst_pair;
  std::printf(
      "  [islamic persistent] palette=%zu B, worst slot=%zu B (idx %zu, "
      "V=%zu F=%zu), worst adj pair=%zu B (idx %zu+%zu), peak=%zu B / "
      "budget=%zu B\n",
      palette_bytes, worst_slot, worst_k, worst_v, worst_f, worst_pair,
      worst_pair_i, (worst_pair_i + 1) % N, peak,
      (size_t)ISLAMIC_PERSISTENT_BUDGET);
  HS_EXPECT_LE(peak, (size_t)ISLAMIC_PERSISTENT_BUDGET);
  // compact_keep_front evacuates the front slot through scratch_b.
  HS_EXPECT_LE(worst_slot, ISLAMIC_SCRATCH_B_BUDGET);
}

// ---------------------------------------------------------------------------
// High-water regression for HankinSolids at its shipping arena configuration.
//
// HankinSolids::init() splits the device arena as configure_arenas(GLOBAL - 24
// KB
// - 32 KB, 24 KB, 32 KB): a 24 KB scratch_a / 32 KB scratch_b pair and the rest
// persistent. scratch_a hosts two non-overlapping peaks, both measured here:
//   * Load: load_shape() runs the whole generate -> compile_hankin ->
//     update_hankin chain inside one generate() call (scratch ping-pongs
//     without an intervening reset), then classify_faces_by_topology() reuses
//     the scratch after generate() rewinds it.
//   * Render: draw_mesh() transforms the front mesh into scratch_a, then
//     Scan::Mesh::draw() stacks an SDF::FaceScratchBuffer on top. For the
//     heaviest hankin mesh this render peak (transformed vertices + the fixed
//     FaceScratchBuffer) exceeds the load peak — the path that actually decides
//     the scratch_a budget, and the one a load-only check missed.
// The two heaviest Archimedean solids (truncatedIcosidodecahedron,
// snubDodecahedron) only cycle into the simple-solid carousel here. Scratch is
// flat POD whose only host/device delta (64-bit pointers) inflates the host
// figure, so the host high-water mark is a conservative upper bound.
// ---------------------------------------------------------------------------

constexpr size_t HANKIN_SCRATCH_A_BUDGET =
    24 * 1024; /**< HankinSolids scratch_a. */
constexpr size_t HANKIN_SCRATCH_B_BUDGET =
    32 * 1024; /**< HankinSolids scratch_b. */
constexpr float HANKIN_ANGLE =
    PI_F / 4.0f; /**< Mid-sweep; counts are angle-independent. */

/**
 * @brief Runs one simple solid through HankinSolids' full load AND render paths
 *        and asserts each scratch arena's peak stays within the 24 KB / 32 KB
 *        split.
 * @param entry Registry entry whose generator is exercised.
 * @details Mirrors load_shape (generate + compile_hankin + update_hankin
 * sharing the scratch pair without a reset, then classify reusing the rewound
 *          scratch) and then draw_mesh (transform the mesh into scratch_a, then
 *          the SDF::FaceScratchBuffer Scan::Mesh::draw stacks on top).
 */
inline void check_hankin_high_water_for_solid(const Solids::Entry &entry) {
  Arena a(solids_scratch_a, sizeof(solids_scratch_a));
  Arena b(solids_scratch_b, sizeof(solids_scratch_b));
  Arena persist(solids_geom_a, sizeof(solids_geom_a));

  PolyMesh base = Solids::finalize_solid(entry.generate(a, b), a);
  CompiledHankin hankin;
  MeshOps::compile_hankin(base, hankin, persist, a);
  MeshState mesh;
  MeshOps::update_hankin(hankin, mesh, persist, HANKIN_ANGLE);
  size_t a_peak = a.get_high_water_mark();
  size_t b_peak = b.get_high_water_mark();

  // generate() rewinds the scratch pair before classify runs.
  a.reset();
  b.reset();
  MeshOps::classify_faces_by_topology(mesh, a, b, persist);
  if (a.get_high_water_mark() > a_peak)
    a_peak = a.get_high_water_mark();
  if (b.get_high_water_mark() > b_peak)
    b_peak = b.get_high_water_mark();

  // Render path: draw_mesh transforms the front mesh into scratch_a (transform
  // copies only the vertices; topology is borrowed), then Scan::Mesh::draw
  // allocates one SDF::FaceScratchBuffer in the same arena. Their sum is the
  // render-time scratch_a peak.
  a.reset();
  {
    ScratchScope render_scope(a);
    MeshState rotated;
    MeshOps::transform(mesh, rotated, a, [](const Vector &v) { return v; });
    (void)a.allocate(sizeof(SDF::FaceScratchBuffer),
                     alignof(SDF::FaceScratchBuffer));
  }
  if (a.get_high_water_mark() > a_peak)
    a_peak = a.get_high_water_mark();

  HS_EXPECT_LE(a_peak, HANKIN_SCRATCH_A_BUDGET);
  HS_EXPECT_LE(b_peak, HANKIN_SCRATCH_B_BUDGET);
}

/**
 * @brief Verifies every simple solid fits HankinSolids' 24 KB / 32 KB scratch
 *        split — including the two heaviest Archimedean solids the carousel now
 *        cycles through.
 */
inline void test_hankin_solids_fit_hankinsolids_scratch_budget() {
  for (const Solids::Entry &e : Solids::Collections::get_simple_solids())
    check_hankin_high_water_for_solid(e);
}

// ---------------------------------------------------------------------------
// Persistent-budget regression for the HankinSolids carousel.
//
// The persistent half is GLOBAL_ARENA_SIZE - 24 KB - 32 KB (~242 KB on device).
// The native 8 MB GLOBAL_ARENA_SIZE means a device overflow can't surface by
// running the effect here. Peak residents during a morph are the baked palette
// bank plus the two adjacent solids that coexist until compaction — each solid
// contributing its compiled-hankin pattern, the rasterized mesh slot, and its
// per-face topology — so the peak is the largest registry-adjacent pair sum,
// not twice the largest single solid.
// ---------------------------------------------------------------------------

constexpr size_t HANKIN_PERSISTENT_BUDGET =
    DEVICE_GLOBAL_ARENA_SIZE - 24 * 1024 - 32 * 1024; /**< ~242 KB on device. */

/**
 * @brief Verifies the worst adjacent pair of simple solids, plus the palette
 *        bank, fits HankinSolids' persistent carousel split.
 * @details Builds each solid's persistent footprint (compiled hankin + mesh +
 *          topology) exactly as load_shape does, records it, then asserts the
 *          largest registry-adjacent pair (the morph's worst coexistence) plus
 *          one palette bank stays within the device budget.
 */
inline void test_hankin_solids_fit_hankinsolids_persistent_budget() {
  size_t palette_bytes;
  {
    Arena pal(solids_geom_b, sizeof(solids_geom_b));
    MeshPaletteBank bank;
    bank.bake_all(pal);
    palette_bytes = pal.get_high_water_mark();
  }

  const auto simple = Solids::Collections::get_simple_solids();
  const size_t N = simple.size();

  size_t slot_bytes[Solids::NUM_ENTRIES];
  size_t worst_slot = 0, worst_k = 0;
  for (size_t k = 0; k < N; ++k) {
    Arena a(solids_scratch_a, sizeof(solids_scratch_a));
    Arena b(solids_scratch_b, sizeof(solids_scratch_b));
    Arena slot(solids_geom_a, sizeof(solids_geom_a));

    PolyMesh base = Solids::finalize_solid(simple[k].generate(a, b), a);
    CompiledHankin hankin;
    MeshOps::compile_hankin(base, hankin, slot, a);
    MeshState mesh;
    MeshOps::update_hankin(hankin, mesh, slot, HANKIN_ANGLE);
    MeshOps::classify_faces_by_topology(mesh, a, b, slot);

    slot_bytes[k] = slot.get_high_water_mark();
    if (slot_bytes[k] > worst_slot) {
      worst_slot = slot_bytes[k];
      worst_k = k;
    }
  }

  // Largest registry-adjacent pair, including the (N-1, 0) cycle wrap.
  size_t worst_pair = 0, worst_pair_i = 0;
  for (size_t i = 0; i < N; ++i) {
    size_t pair = slot_bytes[i] + slot_bytes[(i + 1) % N];
    if (pair > worst_pair) {
      worst_pair = pair;
      worst_pair_i = i;
    }
  }

  const size_t peak = palette_bytes + worst_pair;
  std::printf(
      "  [hankin persistent] palette=%zu B, worst slot=%zu B (%s), worst "
      "adj pair=%zu B (%s+%s), peak=%zu B / budget=%zu B\n",
      palette_bytes, worst_slot, simple[worst_k].name, worst_pair,
      simple[worst_pair_i].name, simple[(worst_pair_i + 1) % N].name, peak,
      (size_t)HANKIN_PERSISTENT_BUDGET);
  HS_EXPECT_LE(peak, (size_t)HANKIN_PERSISTENT_BUDGET);
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

/**
 * @brief Runs all solids tests.
 * @return The module's failure count.
 */
inline int run_solids_tests() {
  hs_test::ModuleFixture fixture("solids");

  test_registry_count_matches_collections();

  test_simple_registry_solids_are_spherical_and_valid();
  test_catalan_registry_solids_are_spherical_and_valid();
  test_islamic_registry_solids_are_valid();
  test_islamic_solids_have_no_sliver_edges();

  test_euler_platonic_solids();
  test_euler_archimedean_catalan_solids();
  test_islamic_registry_solids_are_closed();

  test_get_entry_last_valid_index_builds();
  test_get_by_name_known_returns_that_solid();
  test_registry_names_unique_and_roundtrip();

  test_determinism_hardcoded_platonic();
  test_determinism_archimedean_with_conway_ops();
  test_determinism_complex_islamic();

  test_islamic_recipes_fit_islamicstars_budget();
  test_islamic_solids_fit_islamicstars_persistent_budget();

  test_hankin_solids_fit_hankinsolids_scratch_budget();
  test_hankin_solids_fit_hankinsolids_persistent_budget();

  return fixture.result();
}

} // namespace solids_tests
} // namespace hs_test
