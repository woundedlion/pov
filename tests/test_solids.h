/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
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
 *     test_mesh.h, plus outward-consistent winding on every registry entry.
 *   - Bounds: get_entry() out-of-range and get_by_name() unknown name TRAP
 *     (fail-fast), so only the valid boundary (last index) is exercised here.
 *   - Determinism: building the same registry entry twice yields identical
 *     vertex counts and positions.
 *   - Sliver-face invariant: every Islamic recipe keeps its longest geodesic
 *     edge within 6x the median edge.
 *   - Topology-hash rounding margin: no registered solid has an interior angle
 *     within TOPOLOGY_ANGLE_MARGIN_DEG of the X.5 boundary the classifier's
 *     whole-degree quantisation rounds on.
 *   - Recipe tables: every entry with a non-null recipe replays bitwise-equal
 *     to its generator, lowered (expand_to_primitives) replay is bitwise-equal
 *     to authored replay, and each composite op's decomposition is pinned.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>
#include "core/mesh/mesh.h"
#include "core/color/palettes.h"
#include "core/mesh/recipe.h"
#include "core/mesh/solids.h"
#include "effects/HankinSolids.h"
#include "tests/mesh_test_util.h"
#include "tests/test_conway.h" // check_consistent_winding
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
    HS_EXPECT_TRUE(!edges.empty());
    if (edges.empty())
      continue;
    std::sort(edges.begin(), edges.end());
    float median = edges[edges.size() / 2];
    float max = edges.back();
    HS_EXPECT_TRUE(max <= 6.0f * median);
  }
}

// ---------------------------------------------------------------------------
// Topology-hash rounding margin. classify_faces_impl quantises every interior
// angle to whole degrees and hashes the sorted result into the topology id that
// drives palette assignment. The region compiles -O3 -ffast-math on the device
// image and unmodified on host, so an angle sitting on an X.5 boundary rounds
// one way on host and can round the other on a board: a different palette class
// for the same shape, on a rig that requires every board to agree.
// ---------------------------------------------------------------------------

// Minimum tolerated distance from an X.5 rounding boundary, in degrees. Derived
// from two measured quantities, and asserted against both below: the roster's
// tightest angle clears a boundary by 1.05e-3 deg
// (truncatedOctahedron_gyro_kis_hk17, 98.501053 deg), and the worst float
// rounding uncertainty of the angle expression itself is 2.47e-4 deg
// (dodecahedron_hk35_ambo_hk62_ambo_relax_hk42). The floor sits ~2x above the
// uncertainty and ~2x below the tightest angle.
inline constexpr float TOPOLOGY_ANGLE_MARGIN_DEG = 5e-4f;

/**
 * @brief Interior angle at corner `k` of a face, in degrees, before rounding.
 * @param m Mesh owning the vertices.
 * @param idx Start of the face's index run.
 * @param count Sides in the face.
 * @param k Corner to measure.
 * @return The value classify_faces_impl rounds to a whole degree.
 * @details Mirrors classify_faces_impl (core/mesh/mesh.h) term for term:
 *          unnormalized edges, the same degenerate-edge cutoff, fast_acos.
 */
inline float classifier_angle_deg(const PolyMesh &m, const uint16_t *idx,
                                  int count, int k) {
  const int prev_k = k == 0 ? count - 1 : k - 1;
  const int next_k = k + 1 == count ? 0 : k + 1;
  const Vector e1 = m.vertices[idx[prev_k]] - m.vertices[idx[k]];
  const Vector e2 = m.vertices[idx[next_k]] - m.vertices[idx[k]];
  const float m1 = dot(e1, e1);
  const float m2 = dot(e2, e2);
  if (!(m1 > math::EPS_LEN_SQ && m2 > math::EPS_LEN_SQ))
    return 0.0f;
  const float d = hs::clamp(dot(e1, e2) / sqrtf(m1 * m2), -1.0f, 1.0f);
  return fast_acos(d) * 180.0f / PI_F;
}

/**
 * @brief classifier_angle_deg evaluated in double from the same float vertices.
 * @param m Mesh owning the vertices.
 * @param idx Start of the face's index run.
 * @param count Sides in the face.
 * @param k Corner to measure.
 * @return The same expression, including fast_acos's polynomial, in double.
 * @details Rounding-uncertainty oracle for the margin floor. fast_acos's
 *          approximation error is a property of the polynomial and cancels in
 *          the difference, so what is left is how far float intermediate
 *          rounding alone can move the angle — the same lever -ffast-math pulls
 *          on device via reassociation, FMA contraction, and reciprocal
 *          division.
 */
inline double classifier_angle_deg_ref(const PolyMesh &m, const uint16_t *idx,
                                       int count, int k) {
  constexpr double PI_D = 3.14159265358979323846;
  const int prev_k = k == 0 ? count - 1 : k - 1;
  const int next_k = k + 1 == count ? 0 : k + 1;
  const Vector &prev = m.vertices[idx[prev_k]];
  const Vector &curr = m.vertices[idx[k]];
  const Vector &next = m.vertices[idx[next_k]];
  const double e1x = (double)prev.x - curr.x;
  const double e1y = (double)prev.y - curr.y;
  const double e1z = (double)prev.z - curr.z;
  const double e2x = (double)next.x - curr.x;
  const double e2y = (double)next.y - curr.y;
  const double e2z = (double)next.z - curr.z;
  const double m1 = e1x * e1x + e1y * e1y + e1z * e1z;
  const double m2 = e2x * e2x + e2y * e2y + e2z * e2z;
  if (!(m1 > math::EPS_LEN_SQ && m2 > math::EPS_LEN_SQ))
    return 0.0;
  double d = (e1x * e2x + e1y * e2y + e1z * e2z) / std::sqrt(m1 * m2);
  d = d < -1.0 ? -1.0 : (d > 1.0 ? 1.0 : d);
  const double ax = d < 0.0 ? -d : d;
  const double r =
      std::sqrt(1.0 - ax) *
      (1.5707963 + ax * (-0.2121144 + ax * (0.0742610 + ax * -0.0187293)));
  return (d < 0.0 ? PI_D - r : r) * 180.0 / PI_D;
}

/**
 * @brief Verifies no registered solid has an interior angle within
 *        TOPOLOGY_ANGLE_MARGIN_DEG of a whole-degree rounding boundary, and
 *        that the floor still exceeds the expression's rounding uncertainty.
 * @details An angle at X.5 rounds on a coin flip between the host build and the
 *          device's -ffast-math classifier region, and the rounded angles are
 *          the topology hash, so the two builds would hand the same solid
 *          different palette classes. Sweeps every entry in all three
 *          registries, naming the offending entry, face, and corner. The second
 *          assertion keeps the floor honest: a generator retune that made the
 *          angle expression less well conditioned than the margin would leave a
 *          passing sweep meaningless.
 */
inline void test_registry_angles_clear_rounding_boundary() {
  size_t measured = 0;
  double worst_uncertainty = 0.0;
  for (size_t i = 0; i < Solids::NUM_ENTRIES; ++i) {
    Arena geom(solids_geom_a, sizeof(solids_geom_a));
    PolyMesh m = build_index(i, geom);
    const uint8_t *fc = m.get_face_counts_data();
    const uint16_t *faces = m.get_faces_data();
    const size_t F = m.get_face_counts_size();
    size_t off = 0;
    for (size_t f = 0; f < F; ++f) {
      const int count = fc[f];
      if (count < 3) {
        off += count;
        continue;
      }
      for (int k = 0; k < count; ++k) {
        const float deg = classifier_angle_deg(m, faces + off, count, k);
        const float margin = std::fabs(deg - std::floor(deg) - 0.5f);
        if (margin < TOPOLOGY_ANGLE_MARGIN_DEG)
          std::printf("  [angle margin] %s: face %zu corner %d angle %.6f deg "
                      "sits %.6f deg from an X.5 rounding boundary\n",
                      Solids::get_entry(i).name, f, k, static_cast<double>(deg),
                      static_cast<double>(margin));
        HS_EXPECT_TRUE(margin >= TOPOLOGY_ANGLE_MARGIN_DEG);

        const double uncertainty =
            std::fabs(static_cast<double>(deg) -
                      classifier_angle_deg_ref(m, faces + off, count, k));
        if (uncertainty > worst_uncertainty)
          worst_uncertainty = uncertainty;
        ++measured;
      }
      off += count;
    }
  }
  if (worst_uncertainty >= TOPOLOGY_ANGLE_MARGIN_DEG)
    std::printf("  [angle margin] rounding uncertainty %.6g deg has reached "
                "the %.6g deg margin floor; the sweep no longer proves "
                "anything\n",
                worst_uncertainty,
                static_cast<double>(TOPOLOGY_ANGLE_MARGIN_DEG));
  HS_EXPECT_TRUE(worst_uncertainty < TOPOLOGY_ANGLE_MARGIN_DEG);
  // Guards against a vacuous pass on an empty or all-degenerate sweep.
  HS_EXPECT_TRUE(measured > 0);
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
 * @details Edges are counted as half_edges/2, matching test_mesh.h. Winding is
 *          asserted too: Euler alone passes on an inside-out solid, which the
 *          renderer would shade as backfaces.
 */
inline void check_euler_for_index(size_t index, int expected_V = -1,
                                  int expected_E = -1, int expected_F = -1) {
  Arena geom(solids_geom_a, sizeof(solids_geom_a));
  PolyMesh m = build_index(index, geom);
  conway_tests::check_consistent_winding(m);

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
 * @brief Asserts two meshes are bitwise identical: same counts and same bytes
 *        for vertices, face_counts, and faces.
 * @param m1,m2 Meshes to compare.
 */
inline void check_bitwise_equal_meshes(const PolyMesh &m1, const PolyMesh &m2) {
  HS_EXPECT_EQ(m1.vertices.size(), m2.vertices.size());
  HS_EXPECT_EQ(m1.face_counts.size(), m2.face_counts.size());
  HS_EXPECT_EQ(m1.faces.size(), m2.faces.size());
  if (m1.vertices.size() != m2.vertices.size() ||
      m1.face_counts.size() != m2.face_counts.size() ||
      m1.faces.size() != m2.faces.size())
    return;
  HS_EXPECT_EQ(std::memcmp(m1.vertices.data(), m2.vertices.data(),
                           m1.vertices.size() * sizeof(Vector)),
               0);
  HS_EXPECT_EQ(std::memcmp(m1.face_counts.data(), m2.face_counts.data(),
                           m1.face_counts.size() * sizeof(uint8_t)),
               0);
  HS_EXPECT_EQ(std::memcmp(m1.faces.data(), m2.faces.data(),
                           m1.faces.size() * sizeof(uint16_t)),
               0);
}

/**
 * @brief Builds the entry at `index` twice into separate arenas and asserts the
 *        two meshes are bitwise identical.
 * @param index Registry entry index to build twice.
 * @details A generator is deterministic or it is not: the same input through
 * the same code path twice must reproduce every float bit. A near-equality
 * window would admit a generator whose output drifts below the tolerance.
 */
inline void check_determinism_for_index(size_t index) {
  Arena geom1(solids_geom_a, sizeof(solids_geom_a));
  PolyMesh m1 = build_index(index, geom1);

  Arena geom2(solids_geom_b, sizeof(solids_geom_b));
  PolyMesh m2 = build_index(index, geom2);

  check_bitwise_equal_meshes(m1, m2);
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
// IslamicStars::spawn_shape builds each recipe through the scratch pair its
// init() splits off, and that pair ping-pongs WITHOUT resetting between ops, so
// a recipe chain's peak is its high-water mark. Over-budget would otherwise
// surface only as a device-only OOM trap. Scratch is flat POD whose only
// host/device delta (64-bit pointers) can only make the host figure larger, so
// the host high-water mark is a conservative upper bound on the device figure.
// ---------------------------------------------------------------------------

/** IslamicStars keeps its split private, so these mirror it by value. */
constexpr size_t ISLAMIC_SCRATCH_A_BUDGET =
    116 * 1024; /**< IslamicStars::SPLIT_SCRATCH_A_DEFAULT. */
constexpr size_t ISLAMIC_SCRATCH_B_BUDGET =
    74 * 1024; /**< IslamicStars::SPLIT_SCRATCH_B_DEFAULT. */

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
 */
inline void test_islamic_recipes_fit_islamicstars_budget() {
  for (const Solids::Entry &e : Solids::islamic_registry)
    check_high_water_for_recipe(e);
}

/**
 * @brief Verifies a chain leaves only its seed and its result resident.
 * @details Every operator returns its output in the arena the builder is
 *          writing, so after the role swap the live mesh is in the other one
 *          and the next step's arena holds nothing but spent intermediates. An
 *          even-length chain lands its result in `b`, so `a` ends back at the
 *          offset it held when the chain started — with only the seed below it.
 */
inline void test_solid_builder_reclaims_intermediates() {
  Arena a(solids_scratch_a, sizeof(solids_scratch_a));
  Arena b(solids_scratch_b, sizeof(solids_scratch_b));
  PolyMesh seed = Solids::Platonic::tetrahedron(a, b);
  const size_t a_mark = a.get_offset();
  const size_t b_mark = b.get_offset();

  PolyMesh out = Solids::SolidBuilder(std::move(seed), a, b)
                     .truncate(0.25f)
                     .ambo()
                     .kis()
                     .dual()
                     .ambo()
                     .kis()
                     .build();

  check_nonempty(out);
  check_face_counts_consistent(out);
  check_indices_in_range(out);
  check_all_unit_vertices(out, 1e-4f);

  HS_EXPECT_EQ(a.get_offset(), a_mark);
  HS_EXPECT_GT(b.get_offset(), b_mark);
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

/** The arena split is canvas-independent; this instantiation names it. */
using HankinFx = HankinSolids<96, 20>;

constexpr size_t HANKIN_SCRATCH_A_BUDGET =
    HankinFx::SCRATCH_A_BYTES; /**< HankinSolids scratch_a. */
constexpr size_t HANKIN_SCRATCH_B_BUDGET =
    HankinFx::SCRATCH_B_BYTES; /**< HankinSolids scratch_b. */
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
// Recipe tables: the declarative chains in solids.h must replay bitwise-equal
// to the generator functions they mirror, and expand_to_primitives' composite
// decompositions must replay bitwise-equal to the authored composites.
// ---------------------------------------------------------------------------

/** Lowered-step buffer capacity (deepest chain is 6 steps, meta emits 3). */
constexpr size_t MAX_LOWERED_STEPS = 32;

/**
 * @brief Verifies build_recipe replays every non-null recipe bitwise-identical
 *        to its entry's generator, across all three registries.
 * @details The anchor gate: it is what stops the recipe tables and the shipping
 *          geometry from silently diverging.
 */
inline void test_recipes_match_generators_bitwise() {
  size_t checked = 0;
  for (std::span<const Solids::Entry> reg : Solids::all_registries()) {
    for (const Solids::Entry &e : reg) {
      if (!e.recipe)
        continue;
      Arena geom1(solids_geom_a, sizeof(solids_geom_a));
      PolyMesh gen;
      {
        Arena a(solids_scratch_a, sizeof(solids_scratch_a));
        Arena b(solids_scratch_b, sizeof(solids_scratch_b));
        gen = Solids::finalize_solid(e.generate(a, b), geom1);
      }
      Arena geom2(solids_geom_b, sizeof(solids_geom_b));
      PolyMesh rec;
      {
        Arena a(solids_scratch_a, sizeof(solids_scratch_a));
        Arena b(solids_scratch_b, sizeof(solids_scratch_b));
        rec = Solids::finalize_solid(Solids::build_recipe(*e.recipe, a, b),
                                     geom2);
      }
      check_bitwise_equal_meshes(gen, rec);
      ++checked;
    }
  }
  // Guards against a vacuous pass if the recipe pointers are dropped. Recipes
  // are carried by the Islamic registry alone.
  HS_EXPECT_EQ(checked, Solids::Collections::get_islamic_solids().size());
}

/**
 * @brief Verifies lowered replay (expand_to_primitives + build_steps) matches
 *        authored replay (build_recipe) bitwise for every non-null recipe.
 */
inline void test_recipe_lowered_replay_matches_authored() {
  for (std::span<const Solids::Entry> reg : Solids::all_registries()) {
    for (const Solids::Entry &e : reg) {
      if (!e.recipe)
        continue;
      Solids::OpStep lowered[MAX_LOWERED_STEPS];
      size_t n =
          Solids::expand_to_primitives(*e.recipe, lowered, MAX_LOWERED_STEPS);
      HS_EXPECT_EQ(n, Solids::lowered_step_count(*e.recipe));
      Arena geom1(solids_geom_a, sizeof(solids_geom_a));
      PolyMesh authored;
      {
        Arena a(solids_scratch_a, sizeof(solids_scratch_a));
        Arena b(solids_scratch_b, sizeof(solids_scratch_b));
        authored = Solids::finalize_solid(Solids::build_recipe(*e.recipe, a, b),
                                          geom1);
      }
      Arena geom2(solids_geom_b, sizeof(solids_geom_b));
      PolyMesh low;
      {
        Arena a(solids_scratch_a, sizeof(solids_scratch_a));
        Arena b(solids_scratch_b, sizeof(solids_scratch_b));
        low = Solids::finalize_solid(
            Solids::build_steps(e.recipe->seed, lowered, n, a, b), geom2);
      }
      check_bitwise_equal_meshes(authored, low);
    }
  }
}

/**
 * @brief Verifies expand_to_primitives is a pass-through on the two shipping
 *        recipes (hankin/ambo chains contain no composites).
 */
inline void test_shipping_recipe_lowering_is_identity() {
  const Solids::Entry *entries[] = {
      Solids::find_entry("dodecahedron_hk62_ambo_hk62"),
      Solids::find_entry("octahedron_hk17_ambo_hk73")};
  for (const Solids::Entry *e : entries) {
    HS_EXPECT_TRUE(e != nullptr && e->recipe != nullptr);
    if (!e || !e->recipe)
      continue;
    Solids::OpStep lowered[MAX_LOWERED_STEPS];
    size_t n =
        Solids::expand_to_primitives(*e->recipe, lowered, MAX_LOWERED_STEPS);
    HS_EXPECT_EQ(n, (size_t)e->recipe->count);
    if (n != e->recipe->count)
      continue;
    for (size_t i = 0; i < n; ++i) {
      HS_EXPECT_TRUE(lowered[i].op == e->recipe->steps[i].op);
      HS_EXPECT_EQ(lowered[i].param, e->recipe->steps[i].param);
      HS_EXPECT_EQ(lowered[i].twist, e->recipe->steps[i].twist);
    }
  }
}

/**
 * @brief Builds a one-step recipe on the tetrahedron seed and asserts lowered
 *        replay equals authored replay bitwise.
 * @param step The (composite) op step to pin.
 */
inline void check_composite_lowering(const Solids::OpStep &step) {
  const Solids::Recipe recipe = {0 /* tetrahedron */, &step, 1};
  Solids::OpStep lowered[MAX_LOWERED_STEPS];
  size_t n = Solids::expand_to_primitives(recipe, lowered, MAX_LOWERED_STEPS);
  Arena geom1(solids_geom_a, sizeof(solids_geom_a));
  PolyMesh authored;
  {
    Arena a(solids_scratch_a, sizeof(solids_scratch_a));
    Arena b(solids_scratch_b, sizeof(solids_scratch_b));
    authored =
        Solids::finalize_solid(Solids::build_recipe(recipe, a, b), geom1);
  }
  Arena geom2(solids_geom_b, sizeof(solids_geom_b));
  PolyMesh low;
  {
    Arena a(solids_scratch_a, sizeof(solids_scratch_a));
    Arena b(solids_scratch_b, sizeof(solids_scratch_b));
    low = Solids::finalize_solid(
        Solids::build_steps(recipe.seed, lowered, n, a, b), geom2);
  }
  check_bitwise_equal_meshes(authored, low);
}

/**
 * @brief Pins every composite decomposition (gyro/meta/needle/zip/bevel),
 *        including bevel(0.5) lowering its truncate half to ambo.
 * @details A transposed decomposition (e.g. needle = kd written as dk) fails
 *          here rather than producing a subtly wrong solid at runtime.
 */
inline void test_composite_lowering_matches_composites() {
  check_composite_lowering({Solids::Op::GYRO});
  check_composite_lowering({Solids::Op::META});
  check_composite_lowering({Solids::Op::NEEDLE});
  check_composite_lowering({Solids::Op::ZIP});
  check_composite_lowering({Solids::Op::BEVEL, 0.2f});
  check_composite_lowering({Solids::Op::BEVEL, 0.5f});

  // bevel(0.5)'s lowered form is exactly ambo, ambo (truncate's ambo
  // short-circuit is not a legal leg target).
  const Solids::OpStep bevel_half = {Solids::Op::BEVEL, 0.5f};
  const Solids::Recipe recipe = {0, &bevel_half, 1};
  Solids::OpStep lowered[MAX_LOWERED_STEPS];
  size_t n = Solids::expand_to_primitives(recipe, lowered, MAX_LOWERED_STEPS);
  HS_EXPECT_EQ(n, (size_t)2);
  if (n == 2) {
    HS_EXPECT_TRUE(lowered[0].op == Solids::Op::AMBO);
    HS_EXPECT_TRUE(lowered[1].op == Solids::Op::AMBO);
  }
}

// ---------------------------------------------------------------------------
// Morph feasibility: every Islamic-pattern entry must ship a recipe whose every
// lowered primitive step a morph leg can cover, so an infeasible shape reds CI
// up front instead of silently cutting to a whole-generate fallback at runtime.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies every islamic_registry entry is morph-feasible: it carries a
 *        non-null recipe whose every lowered primitive step satisfies
 *        is_morphable_step.
 * @details Fail-hard gate for the morphing carousel's roster. A future recipe
 *          that lowers to a step no leg can sweep (EXPAND, or a truncate/chamfer
 *          param outside the characterized range) reds here and names the
 *          offending entry and lowered step, rather than cutting to a
 *          whole-generate fallback unannounced. expand_to_primitives mirrors the
 *          lowering the morph path runs, so this checks the same steps it will.
 */
inline void test_islamic_recipes_are_morph_feasible() {
  size_t checked = 0;
  for (const Solids::Entry &e : Solids::islamic_registry) {
    if (e.recipe == nullptr) {
      std::printf("  [morph feasibility] %s: no recipe\n", e.name);
      HS_EXPECT_TRUE(e.recipe != nullptr);
      continue;
    }
    Solids::OpStep lowered[MAX_LOWERED_STEPS];
    size_t n =
        Solids::expand_to_primitives(*e.recipe, lowered, MAX_LOWERED_STEPS);
    for (size_t i = 0; i < n; ++i) {
      if (!Solids::is_morphable_step(lowered[i]))
        std::printf("  [morph feasibility] %s: lowered step %zu "
                    "(op=%d param=%g) is not morphable\n",
                    e.name, i, static_cast<int>(lowered[i].op),
                    static_cast<double>(lowered[i].param));
      HS_EXPECT_TRUE(Solids::is_morphable_step(lowered[i]));
    }
    ++checked;
  }
  // Guards against a vacuous pass: every Islamic entry must have been checked.
  HS_EXPECT_EQ(checked, Solids::Collections::get_islamic_solids().size());
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
  test_registry_angles_clear_rounding_boundary();

  test_euler_platonic_solids();
  test_euler_archimedean_catalan_solids();
  test_islamic_registry_solids_are_closed();

  test_get_entry_last_valid_index_builds();
  test_get_by_name_known_returns_that_solid();
  test_registry_names_unique_and_roundtrip();

  test_determinism_hardcoded_platonic();
  test_determinism_archimedean_with_conway_ops();
  test_determinism_complex_islamic();

  test_recipes_match_generators_bitwise();
  test_recipe_lowered_replay_matches_authored();
  test_shipping_recipe_lowering_is_identity();
  test_composite_lowering_matches_composites();
  test_islamic_recipes_are_morph_feasible();

  test_islamic_recipes_fit_islamicstars_budget();
  test_solid_builder_reclaims_intermediates();
  test_islamic_solids_fit_islamicstars_persistent_budget();

  test_hankin_solids_fit_hankinsolids_scratch_budget();
  test_hankin_solids_fit_hankinsolids_persistent_budget();

  return fixture.result();
}

} // namespace solids_tests
} // namespace hs_test
