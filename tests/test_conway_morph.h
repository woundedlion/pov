/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Operator-level tests for the ConwayMorph transition design
 * (docs/conway_morph_spec.md §7.3–§7.5).
 *
 * Coverage:
 *   - Settle correspondence: relax output vertex order is the identity over its
 *     input (same counts, byte-identical topology, vertex i stays nearest to
 *     input vertex i), so a relaxed endpoint is per-vertex slerpable.
 *   - Bridge convergence: snub(tetrahedron).relax converges to the regular
 *     icosahedron; ambo(tetrahedron) is the regular octahedron.
 *   - Clean-swap invisibility: truncate(seed, 0.5 - eps) vertices pairwise
 *     merge onto ambo(seed) vertices, and each parameterized op's primary
 *     faces at t = T_EPS geometrically match the seed's faces.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include "core/mesh/conway.h"
#include "core/mesh/solids.h"
#include "tests/mesh_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace conway_morph_tests {

inline uint8_t morph_target_buf[256 * 1024]; /**< Op output arena. */
inline uint8_t morph_temp_buf[256 * 1024];   /**< Op scratch arena. */
inline uint8_t morph_aux_buf[256 * 1024];    /**< Seed / second-result arena. */

/** Sweep clamp epsilon for op-at-epsilon endpoints. */
constexpr float T_EPS = 0.02f;

// ---------------------------------------------------------------------------
// Seeds: the solids the edge table sweeps from, including the two ADOPT seeds
// (cuboctahedron, icosidodecahedron), which are ambo of a platonic seed.
// ---------------------------------------------------------------------------

/** @brief Sweep seeds of the ConwayMorph edge table. */
enum class MorphSeed {
  TETRAHEDRON,
  CUBE,
  OCTAHEDRON,
  DODECAHEDRON,
  ICOSAHEDRON,
  CUBOCTAHEDRON,
  ICOSIDODECAHEDRON,
};

inline constexpr MorphSeed MORPH_SEEDS[] = {
    MorphSeed::TETRAHEDRON,       MorphSeed::CUBE,
    MorphSeed::OCTAHEDRON,        MorphSeed::DODECAHEDRON,
    MorphSeed::ICOSAHEDRON,       MorphSeed::CUBOCTAHEDRON,
    MorphSeed::ICOSIDODECAHEDRON,
};

/**
 * @brief Seed name for failure diagnostics.
 * @param s Seed identifier.
 * @return Static name string.
 */
inline const char *seed_name(MorphSeed s) {
  switch (s) {
  case MorphSeed::TETRAHEDRON:
    return "tetrahedron";
  case MorphSeed::CUBE:
    return "cube";
  case MorphSeed::OCTAHEDRON:
    return "octahedron";
  case MorphSeed::DODECAHEDRON:
    return "dodecahedron";
  case MorphSeed::ICOSAHEDRON:
    return "icosahedron";
  case MorphSeed::CUBOCTAHEDRON:
    return "cuboctahedron";
  case MorphSeed::ICOSIDODECAHEDRON:
    return "icosidodecahedron";
  }
  return "?";
}

/**
 * @brief Builds a sweep seed mesh.
 * @param s Seed identifier.
 * @param target Arena receiving the seed mesh.
 * @param temp Scratch arena (holds the platonic base for the ambo seeds).
 * @return The seed PolyMesh in `target`.
 */
inline PolyMesh build_morph_seed(MorphSeed s, Arena &target, Arena &temp) {
  PolyMesh m;
  switch (s) {
  case MorphSeed::TETRAHEDRON:
    build_solid<Solids::Tetrahedron>(m, target);
    return m;
  case MorphSeed::CUBE:
    build_solid<Solids::Cube>(m, target);
    return m;
  case MorphSeed::OCTAHEDRON:
    build_solid<Solids::Octahedron>(m, target);
    return m;
  case MorphSeed::DODECAHEDRON:
    build_solid<Solids::Dodecahedron>(m, target);
    return m;
  case MorphSeed::ICOSAHEDRON:
    build_solid<Solids::Icosahedron>(m, target);
    return m;
  case MorphSeed::CUBOCTAHEDRON:
    build_solid<Solids::Cube>(m, temp);
    return MeshOps::ambo(m, target, temp);
  case MorphSeed::ICOSIDODECAHEDRON:
    build_solid<Solids::Dodecahedron>(m, temp);
    return MeshOps::ambo(m, target, temp);
  }
  return m;
}

// ---------------------------------------------------------------------------
// Shared oracles
// ---------------------------------------------------------------------------

/**
 * @brief Largest absolute deviation of the mesh's face-loop edge lengths from
 *        their mean.
 * @param m Mesh whose edges are measured (each shared edge counted twice; the
 *          uniform double count leaves mean and max deviation unchanged).
 * @return max_e |len(e) - mean_len|.
 */
inline float max_edge_length_deviation(const PolyMesh &m) {
  double sum = 0.0;
  int n = 0;
  size_t off = 0;
  for (size_t fi = 0; fi < m.face_counts.size(); ++fi) {
    const int c = m.face_counts[fi];
    for (int k = 0; k < c; ++k) {
      sum += distance_between(m.vertices[m.faces[off + k]],
                              m.vertices[m.faces[off + (k + 1) % c]]);
      ++n;
    }
    off += c;
  }
  const float mean = static_cast<float>(sum / n);
  float worst = 0.0f;
  off = 0;
  for (size_t fi = 0; fi < m.face_counts.size(); ++fi) {
    const int c = m.face_counts[fi];
    for (int k = 0; k < c; ++k) {
      const float d = distance_between(m.vertices[m.faces[off + k]],
                                       m.vertices[m.faces[off + (k + 1) % c]]) -
                      mean;
      worst = std::max(worst, std::abs(d));
    }
    off += c;
  }
  return worst;
}

/**
 * @brief Verifies the output's primary faces (emitted first, in source-face
 *        order) geometrically match the seed's faces.
 * @param seed Source mesh the operator ran on.
 * @param out Operator output.
 * @param corners_per_source Output corners expected near each seed corner:
 *        2 for truncate (both edge cuts of a corner), 1 for expand/snub.
 * @param tol Max distance from an output corner to its seed corner.
 * @details Primary face fi must have seed_count(fi) * corners_per_source sides
 *          with exactly corners_per_source of them within tol of each seed
 *          corner — pinning emission order, side counts, and geometry at once.
 */
inline void check_primary_faces_match_seed(const PolyMesh &seed,
                                           const PolyMesh &out,
                                           int corners_per_source, float tol) {
  const size_t F = seed.face_counts.size();
  HS_EXPECT_GE(out.face_counts.size(), F);
  size_t seed_off = 0;
  size_t out_off = 0;
  for (size_t fi = 0; fi < F; ++fi) {
    const int bc = seed.face_counts[fi];
    const int oc = out.face_counts[fi];
    HS_EXPECT_EQ(oc, bc * corners_per_source);
    for (int k = 0; k < bc; ++k) {
      const Vector corner = seed.vertices[seed.faces[seed_off + k]];
      int near_count = 0;
      for (int j = 0; j < oc; ++j) {
        if ((out.vertices[out.faces[out_off + j]] - corner).length() <= tol)
          ++near_count;
      }
      HS_EXPECT_EQ(near_count, corners_per_source);
    }
    seed_off += bc;
    out_off += oc;
  }
}

// ---------------------------------------------------------------------------
// §7.3 Settle correspondence: relax output vertex order is the identity over
// its input, so a relaxed endpoint is per-vertex slerpable.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies relax(50) on a registry form (expand(dodecahedron), the
 *        rhombicosidodecahedron chain) preserves vertex order and topology.
 * @details Counts equal, face_counts/faces byte-identical, and each relaxed
 *          vertex stays strictly nearest to its own input vertex — a relax
 *          rewrite that reorders vertices fails here loudly.
 */
inline void test_relax_is_vertex_order_identity() {
  Arena target(morph_target_buf, sizeof(morph_target_buf));
  Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
  Arena aux(morph_aux_buf, sizeof(morph_aux_buf));

  PolyMesh dodeca;
  build_solid<Solids::Dodecahedron>(dodeca, temp);
  PolyMesh unrelaxed = MeshOps::expand(dodeca, target, temp);
  PolyMesh relaxed = MeshOps::relax(unrelaxed, aux, temp, 50);

  HS_EXPECT_EQ(relaxed.vertices.size(), unrelaxed.vertices.size());
  HS_EXPECT_EQ(relaxed.face_counts.size(), unrelaxed.face_counts.size());
  HS_EXPECT_EQ(relaxed.faces.size(), unrelaxed.faces.size());
  HS_EXPECT_EQ(std::memcmp(relaxed.face_counts.data(),
                           unrelaxed.face_counts.data(),
                           relaxed.face_counts.size() * sizeof(uint8_t)),
               0);
  HS_EXPECT_EQ(std::memcmp(relaxed.faces.data(), unrelaxed.faces.data(),
                           relaxed.faces.size() * sizeof(uint16_t)),
               0);

  for (size_t i = 0; i < relaxed.vertices.size(); ++i) {
    size_t nearest = 0;
    float best = 1e9f;
    for (size_t j = 0; j < unrelaxed.vertices.size(); ++j) {
      const float d = (relaxed.vertices[i] - unrelaxed.vertices[j]).length();
      if (d < best) {
        best = d;
        nearest = j;
      }
    }
    HS_EXPECT_EQ(nearest, i);
  }
}

// ---------------------------------------------------------------------------
// §7.4 Bridge convergence: the tetrahedral edges that cross symmetry families.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies snub(tetrahedron, 0.5, 0).relax(50) is the regular
 *        icosahedron: 12 vertices, 20 triangles, equal edges on the unit
 *        sphere (relax supplies the canonical form, as the registry snub
 *        chains rely on).
 */
inline void test_snub_tetrahedron_relax_converges_to_icosahedron() {
  Arena target(morph_target_buf, sizeof(morph_target_buf));
  Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
  Arena aux(morph_aux_buf, sizeof(morph_aux_buf));

  PolyMesh tetra;
  build_solid<Solids::Tetrahedron>(tetra, temp);
  PolyMesh snubbed = MeshOps::snub(tetra, target, temp, 0.5f, 0.0f);
  PolyMesh relaxed = MeshOps::relax(snubbed, aux, temp, 50);

  HS_EXPECT_EQ(relaxed.vertices.size(), (size_t)12);
  HS_EXPECT_EQ(relaxed.face_counts.size(), (size_t)20);
  for (size_t fi = 0; fi < relaxed.face_counts.size(); ++fi)
    HS_EXPECT_EQ((int)relaxed.face_counts[fi], 3);
  check_face_counts_consistent(relaxed);
  check_indices_in_range(relaxed);
  check_all_unit_vertices(relaxed, 1e-3f);

  // Regular icosahedron: every edge equals the mean (chord ~1.0515).
  HS_EXPECT_LE(max_edge_length_deviation(relaxed), 0.02f);
}

/**
 * @brief Verifies ambo(tetrahedron) is the regular octahedron: 6 vertices,
 *        8 triangles, equal edges, and a vertex-set bijection onto the
 *        Octahedron seed (normalized tetra edge midpoints are the ±axes).
 */
inline void test_ambo_tetrahedron_is_regular_octahedron() {
  Arena target(morph_target_buf, sizeof(morph_target_buf));
  Arena temp(morph_temp_buf, sizeof(morph_temp_buf));

  PolyMesh tetra;
  build_solid<Solids::Tetrahedron>(tetra, temp);
  PolyMesh a = MeshOps::ambo(tetra, target, temp);

  HS_EXPECT_EQ(a.vertices.size(), (size_t)6);
  HS_EXPECT_EQ(a.face_counts.size(), (size_t)8);
  for (size_t fi = 0; fi < a.face_counts.size(); ++fi)
    HS_EXPECT_EQ((int)a.face_counts[fi], 3);
  check_face_counts_consistent(a);
  check_indices_in_range(a);
  check_all_unit_vertices(a, 1e-3f);
  HS_EXPECT_LE(max_edge_length_deviation(a), 1e-4f);

  bool used[Solids::Octahedron::NUM_VERTS] = {};
  for (size_t i = 0; i < a.vertices.size(); ++i) {
    int match = -1;
    for (int j = 0; j < Solids::Octahedron::NUM_VERTS; ++j) {
      if (!used[j] &&
          (a.vertices[i] - Solids::Octahedron::vertices[j]).length() <= 1e-4f) {
        match = j;
        break;
      }
    }
    HS_EXPECT_TRUE(match >= 0);
    if (match >= 0)
      used[match] = true;
  }
}

// ---------------------------------------------------------------------------
// §7.5 Clean-swap invisibility: the boundary swaps exchange geometrically
// matching meshes.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies truncate(seed, 0.5 - eps) vertices pairwise merge onto
 *        ambo(seed) vertices for every sweep seed: each ambo vertex has
 *        exactly 2 truncate vertices within tolerance.
 */
inline void test_truncate_near_half_merges_onto_ambo() {
  constexpr float MERGE_EPS = 1e-3f;
  constexpr float MERGE_TOL = 1e-2f;
  for (MorphSeed s : MORPH_SEEDS) {
    Arena target(morph_target_buf, sizeof(morph_target_buf));
    Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
    Arena aux(morph_aux_buf, sizeof(morph_aux_buf));

    PolyMesh seed = build_morph_seed(s, aux, temp);
    PolyMesh tr = MeshOps::truncate(seed, target, temp, 0.5f - MERGE_EPS);
    PolyMesh am = MeshOps::ambo(seed, temp, target);

    HS_EXPECT_EQ(tr.vertices.size(), 2 * am.vertices.size());
    for (size_t i = 0; i < am.vertices.size(); ++i) {
      int merged = 0;
      for (size_t j = 0; j < tr.vertices.size(); ++j) {
        if ((tr.vertices[j] - am.vertices[i]).length() <= MERGE_TOL)
          ++merged;
      }
      if (merged != 2)
        std::printf("    [swap] %s: ambo vertex %zu has %d truncate vertices "
                    "within tol\n",
                    seed_name(s), i, merged);
      HS_EXPECT_EQ(merged, 2);
    }
  }
}

/**
 * @brief Verifies each parameterized op at t = T_EPS emits primary faces that
 *        geometrically match the seed's faces, for every sweep seed.
 * @details truncate contributes two cut corners per seed corner; expand and
 *          snub (zero twist) contribute one inset corner each. Tolerances
 *          bound the T_EPS displacement plus the unit-sphere renormalization.
 */
inline void test_ops_at_t_eps_primary_faces_match_seed() {
  for (MorphSeed s : MORPH_SEEDS) {
    {
      Arena target(morph_target_buf, sizeof(morph_target_buf));
      Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      PolyMesh seed = build_morph_seed(s, aux, temp);
      PolyMesh out = MeshOps::truncate(seed, target, temp, T_EPS);
      check_primary_faces_match_seed(seed, out, /*corners_per_source*/ 2,
                                     /*tol*/ 0.08f);
    }
    {
      Arena target(morph_target_buf, sizeof(morph_target_buf));
      Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      PolyMesh seed = build_morph_seed(s, aux, temp);
      PolyMesh out = MeshOps::expand(seed, target, temp, T_EPS);
      check_primary_faces_match_seed(seed, out, /*corners_per_source*/ 1,
                                     /*tol*/ 0.06f);
    }
    {
      Arena target(morph_target_buf, sizeof(morph_target_buf));
      Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      PolyMesh seed = build_morph_seed(s, aux, temp);
      PolyMesh out = MeshOps::snub(seed, target, temp, T_EPS, 0.0f);
      check_primary_faces_match_seed(seed, out, /*corners_per_source*/ 1,
                                     /*tol*/ 0.06f);
    }
  }
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

/**
 * @brief Runs all ConwayMorph operator-level tests.
 * @return The module's failure count.
 */
inline int run_conway_morph_tests() {
  hs_test::ModuleFixture fixture("conway_morph");

  test_relax_is_vertex_order_identity();

  test_snub_tetrahedron_relax_converges_to_icosahedron();
  test_ambo_tetrahedron_is_regular_octahedron();

  test_truncate_near_half_merges_onto_ambo();
  test_ops_at_t_eps_primary_faces_match_seed();

  return fixture.result();
}

} // namespace conway_morph_tests
} // namespace hs_test
