/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Operator-level tests for the OpLeg transition design
 * (docs/conway_morph_spec.md §7.1–§7.5).
 *
 * Coverage:
 *   - Endpoint exactness: every ConwayGraph edge endpoint on the registry code
 *     path equals the registry generator output exactly; dual-family-seed and
 *     bridge arrivals match within geometric tolerance instead.
 *   - Topology constancy: per edge, samples across the sweep interval hold
 *     constant V/F/I, closed genus-0 manifold, faces >= 3 sides, unit
 *     vertices; per-edge morph-frame scratch peaks fit the HankinSolids
 *     scratch split.
 *   - Settle correspondence: relax output vertex order is the identity over its
 *     input (same counts, byte-identical topology, vertex i stays nearest to
 *     input vertex i), so a relaxed endpoint is per-vertex slerpable.
 *   - Bridge convergence: snub(tetrahedron).relax converges to the regular
 *     icosahedron; ambo(tetrahedron) is the regular octahedron.
 *   - Jitterbug bridge: snub(tetrahedron) at the tabled icosa point is the
 *     regular icosahedron with no relax; at (0.5, -pi/3) its vertices merge
 *     pairwise onto the octahedron; the clamped leg holds V12/F20/E30.
 *   - Clean-swap invisibility: truncate(seed, 0.5 - eps) vertices pairwise
 *     merge onto ambo(seed) vertices, and each parameterized op's primary
 *     faces at t = T_EPS geometrically match the seed's faces.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <map>
#include <vector>
#include "core/color/palettes.h"
#include "core/mesh/conway.h"
#include "core/mesh/conway_graph.h"
#include "core/mesh/hankin.h"
#include "core/mesh/recipe.h"
#include "core/mesh/solids.h"
#include "core/render/canvas.h"
#include "effects/HankinSolids.h"
#include "effects/IslamicStars.h"
#include "tests/mesh_test_util.h"
#include "tests/pixel_test_util.h"
#include "tests/test_conway.h" // check_euler_genus0, face_type_histogram
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace conway_morph_tests {

inline uint8_t morph_target_buf[256 * 1024]; /**< Op output arena. */
inline uint8_t morph_temp_buf[256 * 1024];   /**< Op scratch arena. */
inline uint8_t morph_aux_buf[256 * 1024];    /**< Seed / second-result arena. */
inline uint8_t morph_persist_buf[64 * 1024]; /**< Persistent-seed arena. */

using ConwayGraph::T_EPS;

// ---------------------------------------------------------------------------
// Seeds: the solids the edge table sweeps from, including the two ADOPT seeds
// (cuboctahedron, icosidodecahedron), which are ambo of a platonic seed.
// ---------------------------------------------------------------------------

/** @brief Sweep seeds of the OpLeg edge table. */
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
 * @brief Newell-sum area of face fi.
 */
inline float poly_face_area(const PolyMesh &m, size_t fi) {
  size_t off = 0;
  for (size_t i = 0; i < fi; ++i)
    off += m.face_counts[i];
  const int n = m.face_counts[fi];
  Vector s(0.0f, 0.0f, 0.0f);
  for (int k = 1; k + 1 < n; ++k) {
    const Vector e1 = m.vertices[m.faces[off + k]] - m.vertices[m.faces[off]];
    const Vector e2 =
        m.vertices[m.faces[off + k + 1]] - m.vertices[m.faces[off]];
    s = s + cross(e1, e2);
  }
  return 0.5f * s.length();
}

/**
 * @brief Verifies got's vertices merge pairwise onto want's: exactly two got
 *        vertices within tol of every want vertex.
 */
inline void check_pairwise_vertex_cover(const PolyMesh &got,
                                        const PolyMesh &want, float tol) {
  HS_EXPECT_EQ(got.vertices.size(), 2 * want.vertices.size());
  for (size_t i = 0; i < want.vertices.size(); ++i) {
    int merged = 0;
    for (size_t j = 0; j < got.vertices.size(); ++j) {
      if ((got.vertices[j] - want.vertices[i]).length() <= tol)
        ++merged;
    }
    HS_EXPECT_EQ(merged, 2);
  }
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
/** Corner-match radius at a truncate's T_EPS end, where each seed corner has
 * split into two cut vertices. Above the widest T_EPS cut on the registry
 * seeds and far below their closest corner spacing, so a match cannot alias
 * onto a neighbouring corner. */
constexpr float PRIMARY_CORNER_TOL_TRUNCATE = 0.08f;
/** Same radius for the expand/snub/chamfer eps ends, whose single corner per
 * source moves less. */
constexpr float PRIMARY_CORNER_TOL_SINGLE = 0.06f;

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
 * @brief Verifies snub(tetrahedron, 0.5, SNUB_BRIDGE_TWIST).relax(50) is the
 *        regular icosahedron: 12 vertices, 20 triangles, equal edges on the
 *        unit sphere (relax supplies the canonical form, as the registry snub
 *        chains rely on), at the bridge's tabled arrival twist.
 */
inline void test_snub_tetrahedron_relax_converges_to_icosahedron() {
  Arena target(morph_target_buf, sizeof(morph_target_buf));
  Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
  Arena aux(morph_aux_buf, sizeof(morph_aux_buf));

  PolyMesh tetra;
  build_solid<Solids::Tetrahedron>(tetra, temp);
  PolyMesh snubbed =
      MeshOps::snub(tetra, target, temp, 0.5f, ConwayGraph::SNUB_BRIDGE_TWIST);
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
// Jitterbug bridge (icosahedron <-> octahedron on the tetra snub family):
// both endpoint parameter pins plus the clamped-leg topology sweep.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies snub(tetrahedron, T_JITTERBUG_ICOSA, TWIST_JITTERBUG_ICOSA)
 *        is the regular icosahedron directly — 12 vertices, 20 triangles, all
 *        30 edges equal on the unit sphere with no relax — pinning the
 *        jitterbug bridge's icosa endpoint parameters.
 */
inline void test_jitterbug_icosa_point_is_regular() {
  Arena target(morph_target_buf, sizeof(morph_target_buf));
  Arena temp(morph_temp_buf, sizeof(morph_temp_buf));

  PolyMesh tetra;
  build_solid<Solids::Tetrahedron>(tetra, temp);
  PolyMesh s =
      MeshOps::snub(tetra, target, temp, ConwayGraph::T_JITTERBUG_ICOSA,
                    ConwayGraph::TWIST_JITTERBUG_ICOSA);

  HS_EXPECT_EQ(s.vertices.size(), (size_t)12);
  HS_EXPECT_EQ(s.face_counts.size(), (size_t)20);
  for (size_t fi = 0; fi < s.face_counts.size(); ++fi)
    HS_EXPECT_EQ((int)s.face_counts[fi], 3);
  check_face_counts_consistent(s);
  check_indices_in_range(s);
  check_all_unit_vertices(s, 1e-3f);
  HS_EXPECT_LE(max_edge_length_deviation(s), 1e-5f);
}

/**
 * @brief Verifies the jitterbug octa endpoint snub(tetrahedron, 0.5, -pi/3):
 *        the 12 vertices merge pairwise onto the registry octahedron's 6 and
 *        exactly the 12 edge-orbit faces are zero-area (the SDF zero-area cull
 *        hides them, so the clean swap to the held octahedron changes no
 *        pixels).
 */
inline void test_jitterbug_octa_end_covers_octahedron() {
  Arena target(morph_target_buf, sizeof(morph_target_buf));
  Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
  Arena aux(morph_aux_buf, sizeof(morph_aux_buf));

  PolyMesh tetra;
  build_solid<Solids::Tetrahedron>(tetra, temp);
  PolyMesh s = MeshOps::snub(tetra, target, temp, 0.5f,
                             ConwayGraph::TWIST_JITTERBUG_OCTA);
  PolyMesh octa;
  build_solid<Solids::Octahedron>(octa, aux);

  check_pairwise_vertex_cover(s, octa, 1e-4f);

  // Emission order: 4 primary + 4 vertex-orbit faces (the octahedron's 8),
  // then the 12 collapsed edge-orbit faces.
  HS_EXPECT_EQ(s.face_counts.size(), (size_t)20);
  int zero_area = 0;
  for (size_t fi = 0; fi < s.face_counts.size(); ++fi) {
    const float a = poly_face_area(s, fi);
    if (a < 1e-6f)
      ++zero_area;
    else
      HS_EXPECT_GT(a, 0.5f); // equilateral sqrt(2)-side triangle: ~0.866
    if (fi < 8)
      HS_EXPECT_GT(a, 0.5f);
  }
  HS_EXPECT_EQ(zero_area, 12);
}

/**
 * @brief Verifies the jitterbug leg exactly as OpLeg runs it — t from
 *        the icosa point to the T_EPS_JITTERBUG clamp with the tabled twist
 *        endpoints — holds constant V12/F20/E30 closed genus-0 topology,
 *        >= 3-side faces, and unit vertices, with the collapsing edge never
 *        shorter than the clamp chord (spec section 7.2 for the new edge).
 */
inline void test_jitterbug_sweep_holds_topology() {
  constexpr int SAMPLES = 17;
  Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
  PolyMesh tetra;
  build_solid<Solids::Tetrahedron>(tetra, aux);

  for (int s = 0; s < SAMPLES; ++s) {
    const float k = static_cast<float>(s) / (SAMPLES - 1);
    const float t =
        ConwayGraph::T_JITTERBUG_ICOSA +
        (ConwayGraph::T_EPS_JITTERBUG - ConwayGraph::T_JITTERBUG_ICOSA) * k;
    const float twist = ConwayGraph::TWIST_JITTERBUG_ICOSA +
                        (ConwayGraph::TWIST_JITTERBUG_OCTA -
                         ConwayGraph::TWIST_JITTERBUG_ICOSA) *
                            k;

    Arena target(morph_target_buf, sizeof(morph_target_buf));
    Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
    PolyMesh out = MeshOps::snub(tetra, target, temp, t, twist);

    HS_EXPECT_EQ(out.vertices.size(), (size_t)12);
    HS_EXPECT_EQ(out.face_counts.size(), (size_t)20);
    HS_EXPECT_EQ(out.faces.size(), (size_t)60); // E = I / 2 = 30
    for (size_t fi = 0; fi < out.face_counts.size(); ++fi)
      HS_EXPECT_TRUE(out.face_counts[fi] >= 3);
    check_face_counts_consistent(out);
    check_indices_in_range(out);
    check_all_unit_vertices(out, 1e-3f);
    conway_tests::check_euler_genus0(out);

    // The collapsing edge shrinks monotonically toward the octa end but the
    // clamp keeps it a positive sliver.
    float min_edge = 1e9f;
    size_t off = 0;
    for (size_t fi = 0; fi < out.face_counts.size(); ++fi) {
      const int c = out.face_counts[fi];
      for (int j = 0; j < c; ++j)
        min_edge = std::min(
            min_edge,
            distance_between(out.vertices[out.faces[off + j]],
                             out.vertices[out.faces[off + (j + 1) % c]]));
      off += c;
    }
    HS_EXPECT_GE(min_edge, 0.019f);
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
                                     PRIMARY_CORNER_TOL_TRUNCATE);
    }
    {
      Arena target(morph_target_buf, sizeof(morph_target_buf));
      Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      PolyMesh seed = build_morph_seed(s, aux, temp);
      PolyMesh out = MeshOps::expand(seed, target, temp, T_EPS);
      check_primary_faces_match_seed(seed, out, /*corners_per_source*/ 1,
                                     PRIMARY_CORNER_TOL_SINGLE);
    }
    {
      Arena target(morph_target_buf, sizeof(morph_target_buf));
      Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      PolyMesh seed = build_morph_seed(s, aux, temp);
      PolyMesh out = MeshOps::snub(seed, target, temp, T_EPS, 0.0f);
      check_primary_faces_match_seed(seed, out, /*corners_per_source*/ 1,
                                     PRIMARY_CORNER_TOL_SINGLE);
    }
  }
}

// ---------------------------------------------------------------------------
// §7.1 Endpoint exactness: sweeping to an edge endpoint arrives at the
// registry generator's output — exactly where the composition is the registry
// chain (same code path, same seed frame), within geometric tolerance where it
// is not (dual-family ambo arrivals, bridge arrivals, flash-baked relax
// arrivals, and t = 0 ends, which emit expanded topology).
// ---------------------------------------------------------------------------

/**
 * @brief Runs an edge's operator on a seed at one parameter point.
 * @param e Edge whose op kind is dispatched.
 * @param seed Seed mesh the op runs on.
 * @param target Arena receiving the output mesh.
 * @param temp Arena for the op's transient scratch.
 * @param t Operator parameter.
 * @param twist Snub twist (snub edges only).
 * @return The swept PolyMesh in `target`.
 */
inline PolyMesh run_edge_op(const ConwayGraph::EdgeSpec &e,
                            const PolyMesh &seed, Arena &target, Arena &temp,
                            float t, float twist) {
  switch (e.op) {
  case ConwayGraph::MorphOp::TRUNCATE:
    return MeshOps::truncate(seed, target, temp, t);
  case ConwayGraph::MorphOp::EXPAND:
    return MeshOps::expand(seed, target, temp, t);
  case ConwayGraph::MorphOp::SNUB:
    return MeshOps::snub(seed, target, temp, t, twist);
  case ConwayGraph::MorphOp::CHAMFER:
    return MeshOps::chamfer(seed, target, temp, t);
  }
  return PolyMesh{};
}

/** How an edge endpoint is compared against its node's registry output. */
enum class EndRegime {
  EXACT,        /**< Same code path: bitwise vertices, identical topology. */
  EPS_PRIMARY,  /**< t = 0 end: op(seed, T_EPS) primaries match the seed. */
  VERTEX_MATCH, /**< Same geometry, different vertex order (dual-family ambo,
                     ambo(tetra) bridge). */
  REGULAR,      /**< Relax-canonical arrival in a walk-dependent orientation
                     (tetra -> icosa bridge). */
  PAIR_COVER,   /**< Jitterbug octa end: vertices merge pairwise onto the node
                     mesh's, and the edge-orbit faces are zero-area. */
  BAKED_RELAX,  /**< Registry node ends in relax_baked: identical topology, and
                     vertices within the relax convergence gate. */
};

/**
 * @brief Whether a simple-registry node's generator ends in relax_baked.
 * @param node Simple-registry index.
 * @return True for the two nodes carrying a flash bake.
 * @details Their registry mesh is a table of float bits captured from a host
 *   IEEE relax(), while the leg runs relax() live under the build's own float
 *   semantics — so the two agree bitwise only on a build whose arithmetic
 *   matches the bake's, not on the -ffast-math shipping targets.
 */
inline bool is_relax_baked_node(uint8_t node) {
  return node == ConwayGraph::SNUB_DODECAHEDRON ||
         node == ConwayGraph::TRUNCATED_ICOSIDODECAHEDRON;
}

/**
 * @brief Comparison regime of an edge's to_node end.
 * @param e Edge to classify.
 * @return EXACT when op(seed, t_to) [+ relax] is the to_node registry chain;
 *         BAKED_RELAX when that chain ends in a flash bake rather than the live
 *         relax the leg runs; VERTEX_MATCH for arrivals off the registry seed
 *         (dual-family ambo, non-settle bridges); REGULAR for the settling
 *         bridge, whose relax orientation tracks the seed frame, not the
 *         registry icosahedron; PAIR_COVER for the jitterbug bridge's collapsed
 *         octa end.
 */
inline EndRegime to_end_regime(const ConwayGraph::EdgeSpec &e) {
  using namespace ConwayGraph;
  if (is_jitterbug_edge(e))
    return EndRegime::PAIR_COVER;
  if (e.to_node == CUBOCTAHEDRON && e.seed_solid == OCTAHEDRON)
    return EndRegime::VERTEX_MATCH;
  if (e.to_node == ICOSIDODECAHEDRON && e.seed_solid == ICOSAHEDRON)
    return EndRegime::VERTEX_MATCH;
  if (e.bridge)
    return e.settle ? EndRegime::REGULAR : EndRegime::VERTEX_MATCH;
  if (e.settle && is_relax_baked_node(e.to_node))
    return EndRegime::BAKED_RELAX;
  return EndRegime::EXACT;
}

/**
 * @brief Asserts two meshes are exactly equal: bitwise vertex floats,
 *        identical face_counts and faces arrays.
 */
inline void check_exactly_equal(const PolyMesh &got, const PolyMesh &want) {
  HS_EXPECT_EQ(got.vertices.size(), want.vertices.size());
  HS_EXPECT_EQ(got.face_counts.size(), want.face_counts.size());
  HS_EXPECT_EQ(got.faces.size(), want.faces.size());
  if (got.vertices.size() != want.vertices.size() ||
      got.face_counts.size() != want.face_counts.size() ||
      got.faces.size() != want.faces.size())
    return;
  for (size_t i = 0; i < got.vertices.size(); ++i) {
    HS_EXPECT_EQ(got.vertices[i].x, want.vertices[i].x);
    HS_EXPECT_EQ(got.vertices[i].y, want.vertices[i].y);
    HS_EXPECT_EQ(got.vertices[i].z, want.vertices[i].z);
  }
  for (size_t i = 0; i < got.face_counts.size(); ++i)
    HS_EXPECT_EQ((int)got.face_counts[i], (int)want.face_counts[i]);
  for (size_t i = 0; i < got.faces.size(); ++i)
    HS_EXPECT_EQ(got.faces[i], want.faces[i]);
}

/**
 * @brief Asserts two meshes share a bitwise-identical face list and agree
 *        vertex-for-vertex within the relax convergence gate.
 * @param got Mesh whose settle end ran a live relax.
 * @param want Registry mesh whose chain ends in a flash bake.
 * @details The bake froze relax()'s converged vertices at the generating host's
 *   IEEE arithmetic; a build with different float semantics reaches the gate at
 *   its own converged point. sqrt(RELAX_CONVERGE_EPS_SQ) bounds the step relax
 *   declines to take, and normalization only shrinks it, so two relaxations
 *   that both reach the gate sit at most that far apart. Everything integral —
 *   vertex order, face_counts, faces — stays exact.
 */
inline void check_equal_within_relax_gate(const PolyMesh &got,
                                          const PolyMesh &want) {
  const float tol = std::sqrt(MeshOps::RELAX_CONVERGE_EPS_SQ);
  HS_EXPECT_EQ(got.vertices.size(), want.vertices.size());
  HS_EXPECT_EQ(got.face_counts.size(), want.face_counts.size());
  HS_EXPECT_EQ(got.faces.size(), want.faces.size());
  if (got.vertices.size() != want.vertices.size() ||
      got.face_counts.size() != want.face_counts.size() ||
      got.faces.size() != want.faces.size())
    return;
  for (size_t i = 0; i < got.vertices.size(); ++i)
    HS_EXPECT_LE(distance_between(got.vertices[i], want.vertices[i]), tol);
  for (size_t i = 0; i < got.face_counts.size(); ++i)
    HS_EXPECT_EQ((int)got.face_counts[i], (int)want.face_counts[i]);
  for (size_t i = 0; i < got.faces.size(); ++i)
    HS_EXPECT_EQ(got.faces[i], want.faces[i]);
}

/**
 * @brief Asserts two meshes carry the same geometry up to vertex order:
 *        equal counts, equal face-type histograms, and a vertex-set bijection
 *        within tol.
 */
inline void check_equal_up_to_vertex_order(const PolyMesh &got,
                                           const PolyMesh &want, float tol) {
  HS_EXPECT_EQ(got.vertices.size(), want.vertices.size());
  HS_EXPECT_EQ(got.face_counts.size(), want.face_counts.size());
  HS_EXPECT_EQ(got.faces.size(), want.faces.size());
  HS_EXPECT_TRUE(conway_tests::face_type_histogram(got) ==
                 conway_tests::face_type_histogram(want));
  if (got.vertices.size() != want.vertices.size())
    return;
  std::vector<bool> used(want.vertices.size(), false);
  for (size_t i = 0; i < got.vertices.size(); ++i) {
    bool matched = false;
    for (size_t j = 0; j < want.vertices.size(); ++j) {
      if (!used[j] && (got.vertices[i] - want.vertices[j]).length() <= tol) {
        used[j] = true;
        matched = true;
        break;
      }
    }
    HS_EXPECT_TRUE(matched);
  }
}

/**
 * @brief Asserts a mesh is the registry solid's regular form in an arbitrary
 *        orientation: equal counts, equal face-type histograms, unit vertices,
 *        and near-equal edge lengths.
 */
inline void check_regular_form(const PolyMesh &got, const PolyMesh &want,
                               float edge_dev_tol) {
  HS_EXPECT_EQ(got.vertices.size(), want.vertices.size());
  HS_EXPECT_EQ(got.face_counts.size(), want.face_counts.size());
  HS_EXPECT_EQ(got.faces.size(), want.faces.size());
  HS_EXPECT_TRUE(conway_tests::face_type_histogram(got) ==
                 conway_tests::face_type_histogram(want));
  check_all_unit_vertices(got, 1e-3f);
  HS_EXPECT_LE(max_edge_length_deviation(got), edge_dev_tol);
}

/**
 * @brief Verifies every ConwayGraph edge endpoint against its node's registry
 *        generator: exact on the registry code path, geometric tolerance for
 *        the t = 0 ends, the off-registry-seed arrivals and the flash-baked
 *        relax arrivals.
 * @details Seeds are built via the registry generators, so the DERIVE_AMBO
 *          rows (cuboctahedron / icosidodecahedron seeds) run the exact bevel
 *          decomposition of their to_node chains. from ends at t = 0 use the
 *          EPS_PRIMARY regime (an op at 0 emits expanded topology with
 *          coincident positions, never the seed mesh itself).
 */
inline void test_edge_endpoints_match_registry() {
  constexpr size_t HALF = sizeof(morph_target_buf) / 2;
  for (int ei = 0; ei < ConwayGraph::NUM_EDGES; ++ei) {
    const ConwayGraph::EdgeSpec &e = ConwayGraph::EDGES[ei];
    const int failed_before = hs_test::stats().failed;

    Arena sa(morph_aux_buf, HALF);
    Arena sb(morph_aux_buf + HALF, HALF);
    PolyMesh seed = Solids::simple_registry[e.seed_solid].generate(sa, sb);

    // from end: t = 0 emits expanded topology, so compare op(seed, T_EPS)
    // primaries against the seed (= the from_node registry mesh); a non-zero
    // t_from is the from_node registry chain itself, except the jitterbug
    // icosa point, which is regular in the tetra frame, not the registry
    // orientation.
    {
      Arena oa(morph_temp_buf, HALF);
      Arena ob(morph_temp_buf + HALF, HALF);
      if (e.t_from == 0.0f) {
        PolyMesh got = run_edge_op(e, seed, oa, ob, T_EPS, e.twist_from);
        const int per_corner = e.op == ConwayGraph::MorphOp::TRUNCATE ? 2 : 1;
        const float tol = e.op == ConwayGraph::MorphOp::TRUNCATE
                              ? PRIMARY_CORNER_TOL_TRUNCATE
                              : PRIMARY_CORNER_TOL_SINGLE;
        check_primary_faces_match_seed(seed, got, per_corner, tol);
      } else {
        Arena ra(morph_target_buf, HALF);
        Arena rb(morph_target_buf + HALF, HALF);
        PolyMesh want = Solids::simple_registry[e.from_node].generate(ra, rb);
        PolyMesh got = run_edge_op(e, seed, oa, ob, e.t_from, e.twist_from);
        if (ConwayGraph::is_jitterbug_edge(e))
          check_regular_form(got, want, 1e-4f);
        else
          check_exactly_equal(got, want);
      }
    }

    // to end.
    {
      Arena ra(morph_target_buf, HALF);
      Arena rb(morph_target_buf + HALF, HALF);
      PolyMesh want = Solids::simple_registry[e.to_node].generate(ra, rb);

      Arena oa(morph_temp_buf, HALF);
      Arena ob(morph_temp_buf + HALF, HALF);
      PolyMesh got = run_edge_op(e, seed, oa, ob, e.t_to, e.twist_to);
      if (e.settle)
        got = MeshOps::relax(got, ob, oa, 50);

      switch (to_end_regime(e)) {
      case EndRegime::EXACT:
        check_exactly_equal(got, want);
        break;
      case EndRegime::BAKED_RELAX:
        check_equal_within_relax_gate(got, want);
        break;
      case EndRegime::VERTEX_MATCH:
        check_equal_up_to_vertex_order(got, want, 1e-4f);
        break;
      case EndRegime::REGULAR:
        check_regular_form(got, want, 0.02f);
        break;
      case EndRegime::PAIR_COVER:
        check_pairwise_vertex_cover(got, want, 1e-4f);
        break;
      case EndRegime::EPS_PRIMARY:
        break; // from-end-only regime
      }
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [endpoint] edge %d: %s -> %s\n", ei,
                  Solids::simple_registry[e.from_node].name,
                  Solids::simple_registry[e.to_node].name);
  }
}

// ---------------------------------------------------------------------------
// §7.2 Topology-constancy sweep: connectivity is fixed on the open interval,
// so classification and palette assignment can hoist to once per leg.
// ---------------------------------------------------------------------------

/** Samples per edge sweep. */
constexpr int SWEEP_SAMPLES = 16;

/**
 * @brief Per-leg accumulator for the invariants every OpLeg draw callback
 *        shares: a compiled face count matching the shading and constant
 *        across the leg, in-range ramp indices, and a vertex list that keeps
 *        its count and order frame to frame.
 */
struct LegDrawProbe {
  size_t drawn = 0;           /**< Frames handed to the callback. */
  size_t faces = 0;           /**< Face count latched on the first frame. */
  float worst_step = 0.0f;    /**< Largest per-vertex motion between frames. */
  std::vector<Vector> prev_v; /**< Previous frame's vertices. */

  /**
   * @brief Folds one drawn frame in.
   * @param m Compiled mesh handed to the callback.
   * @param sh Shading handed alongside it.
   */
  void observe(const MeshState &m, const Animation::OpLeg::Shading &sh) {
    HS_EXPECT_EQ(m.face_counts.size(), sh.faces);
    if (drawn == 0)
      faces = sh.faces;
    else
      HS_EXPECT_EQ(sh.faces, faces);
    for (size_t f = 0; f < sh.faces; ++f)
      HS_EXPECT_LT(static_cast<int>(sh.face_ramp[f]),
                   Animation::OpLeg::MAX_BLEND_PAIRS);
    if (!prev_v.empty()) {
      HS_EXPECT_EQ(m.vertices.size(), prev_v.size());
      for (size_t i = 0; i < m.vertices.size(); ++i)
        worst_step =
            std::max(worst_step, distance_between(m.vertices[i], prev_v[i]));
    }
    prev_v.assign(m.vertices.begin(), m.vertices.end());
    ++drawn;
  }
};

/** Raw mesh counts a Conway operator is required to produce. */
struct OpCounts {
  size_t v; /**< Vertices. */
  size_t f; /**< Faces. */
  size_t i; /**< Flat face indices; on a closed mesh this is 2E. */
};

/**
 * @brief Reads a closed mesh's operand counts.
 * @param m Closed manifold mesh.
 * @return {V, F, I} of @p m.
 */
inline OpCounts mesh_op_counts(const PolyMesh &m) {
  return {m.vertices.size(), m.face_counts.size(), m.faces.size()};
}

/**
 * @brief Counts a closed mesh's degree-2 vertices.
 * @param m Closed manifold mesh; each index is one face incidence.
 * @return Vertices with exactly two incident faces.
 * @details Truncating one produces a 2-gon the op drops, so this is the only
 *          term the Conway forms need beyond (V, E, F). Hankin seeds carry
 *          degree-2 star points; the registry solids carry none.
 */
inline size_t degree_two_vertex_count(const PolyMesh &m) {
  std::vector<uint32_t> degree(m.vertices.size(), 0u);
  for (size_t k = 0; k < m.faces.size(); ++k)
    ++degree[m.faces[k]];
  size_t n = 0;
  for (uint32_t d : degree)
    n += d == 2u ? 1u : 0u;
  return n;
}

/**
 * @brief Closed-form output counts of each Conway operator the graph sweeps.
 * @param op Operator applied.
 * @param seed Seed mesh.
 * @return The counts @p op must produce from @p seed, at every sweep parameter
 *         on the open interval.
 * @details Standard Conway arithmetic in the seed's (V, E, F): truncate
 *          (2E, F+V, 3E), expand (2E, F+V+E, 4E), snub (2E, F+V+2E, 5E),
 *          chamfer (V+2E, F+E, 4E) as (V', F', E'); the returned index count
 *          is 2E'. Truncate additionally drops the 2-gon each degree-2 vertex
 *          would raise, costing one face and one edge apiece; the other three
 *          are the min-degree-3 forms, which is what their seeds are. Deriving
 *          these from the seed is what separates "constant across the sweep"
 *          from "constant and correct": a latched sample cannot tell a
 *          wrong-but-manifold count apart.
 */
inline OpCounts morph_op_counts(ConwayGraph::MorphOp op, const PolyMesh &seed) {
  const size_t v = seed.vertices.size();
  const size_t f = seed.face_counts.size();
  const size_t e = seed.faces.size() / 2;
  switch (op) {
  case ConwayGraph::MorphOp::TRUNCATE: {
    const size_t d2 = degree_two_vertex_count(seed);
    return {2 * e, f + v - d2, 6 * e - 2 * d2};
  }
  case ConwayGraph::MorphOp::EXPAND:
    return {2 * e, f + v + e, 8 * e};
  case ConwayGraph::MorphOp::SNUB:
    return {2 * e, f + v + 2 * e, 10 * e};
  case ConwayGraph::MorphOp::CHAMFER:
    return {v + 2 * e, f + e, 8 * e};
  }
  return {v, f, seed.faces.size()};
}

/**
 * @brief Closed-form output counts of MeshOps::hankin, at any contact angle.
 * @param seed Seed mesh.
 * @return {3E, F+V, 8E}.
 * @details compile_hankin lays down one midpoint per edge and one star point
 *          per half-edge (3E vertices), one star face per seed face and one
 *          rosette per seed vertex (F+V), and each face contributes twice its
 *          degree in indices on both sides (4I = 8E). Degree-2 vertices raise
 *          a quad rosette like any other, so no correction applies.
 */
inline OpCounts hankin_op_counts(const PolyMesh &seed) {
  const size_t e = seed.faces.size() / 2;
  return {3 * e, seed.face_counts.size() + seed.vertices.size(), 8 * e};
}

/**
 * @brief Asserts a sweep sample's raw counts equal a derived expectation.
 * @param v Sample vertex count.
 * @param f Sample face count.
 * @param i Sample flat index count.
 * @param d Expected counts.
 */
inline void expect_op_counts(size_t v, size_t f, size_t i, const OpCounts &d) {
  HS_EXPECT_EQ(v, d.v);
  HS_EXPECT_EQ(f, d.f);
  HS_EXPECT_EQ(i, d.i);
}

/**
 * @brief Sweep interval of an edge, clamped as a leg runs it.
 * @param e Edge to clamp.
 * @param t_lo Out: max(t_from, T_EPS).
 * @param t_hi Out: t_to, additionally capped at 0.5 - T_EPS on truncate legs
 *        (the ambo short-circuit changes emission order and face count) and
 *        held at T_EPS_JITTERBUG on the jitterbug bridge (the t = 0.5 end is
 *        the pairwise-merged octahedron).
 */
inline void edge_sweep_interval(const ConwayGraph::EdgeSpec &e, float &t_lo,
                                float &t_hi) {
  t_lo = std::max(e.t_from, T_EPS);
  t_hi = e.op == ConwayGraph::MorphOp::TRUNCATE ? std::min(e.t_to, 0.5f - T_EPS)
                                                : e.t_to;
  if (ConwayGraph::is_jitterbug_edge(e))
    t_hi = std::max(t_hi, ConwayGraph::T_EPS_JITTERBUG);
}

/**
 * @brief Verifies every edge holds constant topology across its sweep:
 *        fixed V/F/I, closed genus-0 manifold, all faces >= 3 sides, unit
 *        vertices, no traps.
 * @details Snub twist interpolates linearly with t, as a leg sweeps it.
 */
inline void test_edge_sweeps_hold_topology() {
  constexpr size_t HALF = sizeof(morph_aux_buf) / 2;
  for (int ei = 0; ei < ConwayGraph::NUM_EDGES; ++ei) {
    const ConwayGraph::EdgeSpec &e = ConwayGraph::EDGES[ei];
    const int failed_before = hs_test::stats().failed;

    Arena sa(morph_aux_buf, HALF);
    Arena sb(morph_aux_buf + HALF, HALF);
    PolyMesh seed = Solids::simple_registry[e.seed_solid].generate(sa, sb);

    float t_lo, t_hi;
    edge_sweep_interval(e, t_lo, t_hi);

    size_t v0 = 0, f0 = 0, i0 = 0;
    for (int s = 0; s < SWEEP_SAMPLES; ++s) {
      const float u = static_cast<float>(s) / (SWEEP_SAMPLES - 1);
      const float t = t_lo + (t_hi - t_lo) * u;
      const float twist =
          e.twist_from +
          (e.twist_to - e.twist_from) * ((t - e.t_from) / (e.t_to - e.t_from));

      Arena target(morph_target_buf, sizeof(morph_target_buf));
      Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
      PolyMesh out = run_edge_op(e, seed, target, temp, t, twist);

      if (s == 0) {
        v0 = out.vertices.size();
        f0 = out.face_counts.size();
        i0 = out.faces.size();
        expect_op_counts(v0, f0, i0, morph_op_counts(e.op, seed));
      } else {
        HS_EXPECT_EQ(out.vertices.size(), v0);
        HS_EXPECT_EQ(out.face_counts.size(), f0);
        HS_EXPECT_EQ(out.faces.size(), i0);
      }
      for (size_t fi = 0; fi < out.face_counts.size(); ++fi)
        HS_EXPECT_TRUE(out.face_counts[fi] >= 3);
      check_face_counts_consistent(out);
      check_indices_in_range(out);
      check_all_unit_vertices(out, 1e-3f);
      conway_tests::check_euler_genus0(out);
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [sweep] edge %d: %s -> %s\n", ei,
                  Solids::simple_registry[e.from_node].name,
                  Solids::simple_registry[e.to_node].name);
  }
}

// ---------------------------------------------------------------------------
// Morph-frame scratch high-water gate at HankinSolids' shipping split
// (mirrors the test_solids.h HANKIN_SCRATCH_A/B_BUDGET idiom). A morph frame
// runs one op plus MeshOps::compile in the scratch pair under LIFO scopes;
// host high-water marks are a conservative upper bound on the device figure.
// ---------------------------------------------------------------------------

/** The arena split is canvas-independent; this instantiation names it. */
using HankinFx = HankinSolids<96, 20>;

constexpr size_t MORPH_SCRATCH_A_BUDGET =
    HankinFx::SCRATCH_A_BYTES; /**< HankinSolids scratch_a split. */
constexpr size_t MORPH_SCRATCH_B_BUDGET =
    HankinFx::SCRATCH_B_BYTES; /**< HankinSolids scratch_b split. */

/**
 * @brief Verifies every edge's per-frame scratch peak (op + compile into a
 *        fresh arena pair) fits HankinSolids' 24 KB / 32 KB scratch split.
 * @details The seed lives in a persistent arena as the effect holds it;
 *          topology is t-constant, so one mid-sweep sample per edge is the
 *          frame peak. Reports the worst pair across the table.
 */
inline void test_edge_morph_frames_fit_scratch_budget() {
  constexpr size_t HALF = sizeof(morph_aux_buf) / 2;
  size_t worst_a = 0, worst_b = 0;
  int worst_a_edge = 0, worst_b_edge = 0;

  for (int ei = 0; ei < ConwayGraph::NUM_EDGES; ++ei) {
    const ConwayGraph::EdgeSpec &e = ConwayGraph::EDGES[ei];

    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed;
    {
      Arena ga(morph_aux_buf, HALF);
      Arena gb(morph_aux_buf + HALF, HALF);
      seed = Solids::finalize_solid(
          Solids::simple_registry[e.seed_solid].generate(ga, gb), persist);
    }

    float t_lo, t_hi;
    edge_sweep_interval(e, t_lo, t_hi);
    const float t = (t_lo + t_hi) * 0.5f;
    const float twist =
        e.twist_from +
        (e.twist_to - e.twist_from) * ((t - e.t_from) / (e.t_to - e.t_from));

    Arena a(morph_target_buf, MORPH_SCRATCH_A_BUDGET);
    Arena b(morph_temp_buf, MORPH_SCRATCH_B_BUDGET);
    {
      ScratchScope frame_a(a);
      ScratchScope frame_b(b);
      PolyMesh swept = run_edge_op(e, seed, a, b, t, twist);
      MeshState frame;
      MeshOps::compile(swept, frame, a, b);
    }

    const size_t a_peak = a.get_high_water_mark();
    const size_t b_peak = b.get_high_water_mark();
    if (a_peak > worst_a) {
      worst_a = a_peak;
      worst_a_edge = ei;
    }
    if (b_peak > worst_b) {
      worst_b = b_peak;
      worst_b_edge = ei;
    }
    HS_EXPECT_LE(a_peak, MORPH_SCRATCH_A_BUDGET);
    HS_EXPECT_LE(b_peak, MORPH_SCRATCH_B_BUDGET);
  }

  const ConwayGraph::EdgeSpec &wa = ConwayGraph::EDGES[worst_a_edge];
  const ConwayGraph::EdgeSpec &wb = ConwayGraph::EDGES[worst_b_edge];
  std::printf(
      "  [morph scratch] worst a=%zu B (%s -> %s) / budget=%zu B, "
      "worst b=%zu B (%s -> %s) / budget=%zu B\n",
      worst_a, Solids::simple_registry[wa.from_node].name,
      Solids::simple_registry[wa.to_node].name, (size_t)MORPH_SCRATCH_A_BUDGET,
      worst_b, Solids::simple_registry[wb.from_node].name,
      Solids::simple_registry[wb.to_node].name, (size_t)MORPH_SCRATCH_B_BUDGET);
}

// ---------------------------------------------------------------------------
// Walk policy: the recency-weighted random walk visits every node within a
// bounded leg count and keeps long-run visitation balanced (no node above 3x
// the mean share, none below a quarter of it) — the hub-heavy degree-
// proportional bias this replaced gave cuboctahedron ~13x the pendant rate.
// ---------------------------------------------------------------------------

/** Legs within which every node must have been visited (measured worst over
 * 200 seeds: 200; the tested seeds reach it by 166). */
constexpr int WALK_COVERAGE_BOUND = 250;

/**
 * @brief Whether a completed leg reseeds the held seed from its arrival.
 * @param e Edge spec the leg ran on.
 * @param arrived Registry index of the node the leg landed on.
 * @param arrived_at_to True when the leg ran forward (landed on to_node).
 * @return True iff HankinSolids::finish_morph_cycle adopts here.
 * @details Mirror of the production gate; the reverse jitterbug arrival adopts
 * the icosahedron too, holding its canonical relax form rather than the
 * unrelaxed jitterbug mesh.
 */
inline bool leg_adopts_seed(const ConwayGraph::EdgeSpec &e, int arrived,
                            bool arrived_at_to) {
  return e.reseed == ConwayGraph::Reseed::ADOPT &&
         ConwayGraph::is_platonic(arrived) &&
         (arrived_at_to || arrived == ConwayGraph::ICOSAHEDRON);
}

/**
 * @brief Applies one leg's seed reconciliation to a held seed identity.
 * @param edge Index into ConwayGraph::EDGES of the leg about to start.
 * @param node Node the leg departs from.
 * @param held Held seed identity, updated in place by a non-KEEP fix.
 * @details Mirror of HankinSolids::start_morph_cycle's reconciliation, minus
 * the mesh rebuilds; the HS_CHECKs it pins there become expectations here.
 */
inline void reconcile_seed(int edge, int node, int &held) {
  using namespace ConwayGraph;
  const SeedFix fix = seed_fix_at_start(edge, held);
  HS_EXPECT_TRUE(fix != SeedFix::INVALID);
  if (fix == SeedFix::DUAL_SWAP) {
    HS_EXPECT_TRUE(node == CUBOCTAHEDRON || node == ICOSIDODECAHEDRON);
    held = dual_platonic(held);
  } else if (fix == SeedFix::REGEN_TETRA) {
    HS_EXPECT_TRUE(node == OCTAHEDRON || node == ICOSAHEDRON);
    held = TETRAHEDRON;
  }
  HS_EXPECT_TRUE(fix == SeedFix::DERIVE_AMBO || held == EDGES[edge].seed_solid);
}

/**
 * @brief Simulates long walks over several RNG seeds and pins coverage and
 *        per-node share bounds.
 */
inline void test_walk_policy_coverage_and_balance() {
  using namespace ConwayGraph;
  constexpr int LEGS = 10000;
  constexpr int MEAN = LEGS / NUM_NODES;

  for (uint32_t seed : {1u, 2u, 3u, 42u, 1337u}) {
    const int failed_before = hs_test::stats().failed;
    hs::random().seed(seed);
    uint8_t visits[NUM_NODES] = {};
    int counts[NUM_NODES] = {};
    int node = TETRAHEDRON;
    int held = TETRAHEDRON;
    int prev = -1;
    int in_family = 0;
    bool seen[NUM_NODES] = {};
    int seen_count = 1;
    int coverage_leg = -1;
    seen[node] = true;
    record_visit(visits, node);

    for (int leg = 0; leg < LEGS; ++leg) {
      const int e = pick_next_edge(node, prev, in_family, visits,
                                   static_cast<uint32_t>(hs::random()()));
      HS_EXPECT_TRUE(edge_touches(e, node));
      HS_EXPECT_TRUE(e != prev || node_degree(node) == 1);
      reconcile_seed(e, node, held);
      const int next = edge_other_end(e, node);
      in_family = family(next) != family(node) ? 0 : in_family + 1;
      if (leg_adopts_seed(EDGES[e], next, EDGES[e].to_node == next))
        held = next;
      node = next;
      prev = e;
      ++counts[node];
      record_visit(visits, node);
      if (!seen[node]) {
        seen[node] = true;
        if (++seen_count == NUM_NODES)
          coverage_leg = leg + 1;
      }
    }

    HS_EXPECT_GT(coverage_leg, 0);
    HS_EXPECT_LE(coverage_leg, WALK_COVERAGE_BOUND);
    int mn = counts[0], mx = counts[0];
    for (int i = 0; i < NUM_NODES; ++i) {
      HS_EXPECT_LE(counts[i], 3 * MEAN);
      HS_EXPECT_GE(counts[i], MEAN / 4);
      mn = std::min(mn, counts[i]);
      mx = std::max(mx, counts[i]);
    }
    if (hs_test::stats().failed == failed_before)
      std::printf("  [walk] seed %u: coverage@%d legs, share max/min = "
                  "%d/%d, max/mean = %.2f\n",
                  seed, coverage_leg, mx, mn, static_cast<double>(mx) / MEAN);
    else
      std::printf("    [walk] seed %u failed (coverage@%d, max %d, min %d)\n",
                  seed, coverage_leg, mx, mn);
  }
}

// ---------------------------------------------------------------------------
// Ordered profile tour: one ORDERED_TOUR pass covers all 18 nodes with a
// legal seed reconciliation at every leg, and the cycle wraps back to the
// registry start state so repeated passes are identical.
// ---------------------------------------------------------------------------

/**
 * @brief Simulates the walk state machine over two ordered-tour cycles.
 */
inline void test_ordered_tour_full_coverage_and_wrap() {
  using namespace ConwayGraph;
  int node = TETRAHEDRON;
  int held = TETRAHEDRON;
  int prev = -1;
  bool seen[NUM_NODES] = {};
  int seen_count = 1;
  int coverage_leg = -1;
  seen[node] = true;

  for (uint32_t leg = 0; leg < 2u * ORDERED_TOUR_LEN; ++leg) {
    const int e = pick_next_edge_ordered(node, prev, leg);
    HS_EXPECT_TRUE(edge_touches(e, node));

    reconcile_seed(e, node, held);

    const bool reverse = EDGES[e].to_node == node;
    const int arrived = reverse ? EDGES[e].from_node : EDGES[e].to_node;
    if (leg_adopts_seed(EDGES[e], arrived, !reverse))
      held = arrived;
    node = arrived;
    prev = e;
    if (!seen[node]) {
      seen[node] = true;
      if (++seen_count == NUM_NODES)
        coverage_leg = static_cast<int>(leg) + 1;
    }

    // Cycle closure: every pass ends back at the registry start state.
    if ((leg + 1) % ORDERED_TOUR_LEN == 0) {
      HS_EXPECT_EQ(node, (int)TETRAHEDRON);
      HS_EXPECT_EQ(held, (int)TETRAHEDRON);
    }
  }
  HS_EXPECT_GT(coverage_leg, 0);
  HS_EXPECT_LE(coverage_leg, ORDERED_TOUR_LEN);
  std::printf("  [tour] %d legs per cycle, full coverage after %d\n",
              ORDERED_TOUR_LEN, coverage_leg);
}

// ---------------------------------------------------------------------------
// Ambo-on-hankin probe (docs/opchain_morph_spec.md section 10 Phase 0): every
// ambo leg the Islamic recipes run on a hankin mesh is a truncate sweep whose
// compiled face count must not move within the leg and whose swept mesh must
// stay a closed genus-0 manifold at every parameter.
// ---------------------------------------------------------------------------

/** @brief One ambo-on-hankin sweep seed from the Islamic registry chains. */
struct HankinAmboSite {
  const char *name;                     /**< Diagnostic label. */
  PolyMesh (*seed)(Arena &a, Arena &b); /**< Chain prefix up to the ambo. */
};

inline PolyMesh probe_dodeca_hk62(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
      .hankin(62.0f * D2R)
      .build();
}
inline PolyMesh probe_dodeca_hk35(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
      .hankin(35.0f * D2R)
      .build();
}
inline PolyMesh probe_dodeca_hk35_ambo_hk62(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
      .hankin(35.0f * D2R)
      .ambo()
      .hankin(62.0f * D2R)
      .build();
}
inline PolyMesh probe_dodeca_hk54(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
      .hankin(54.0f * D2R)
      .build();
}
inline PolyMesh probe_octa_hk17(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::octahedron(a, b), a, b)
      .hankin(17.0f * D2R)
      .build();
}
inline PolyMesh probe_octa_hk34(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::octahedron(a, b), a, b)
      .hankin(34.0f * D2R)
      .build();
}
inline PolyMesh probe_rhombicubocta_hk63(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Archimedean::rhombicuboctahedron(a, b), a,
                              b)
      .hankin(63.0f * D2R)
      .build();
}
inline PolyMesh probe_ticosa_hk54(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Archimedean::truncatedIcosahedron(a, b),
                              a, b)
      .hankin(54.0f * D2R)
      .build();
}

inline constexpr HankinAmboSite HANKIN_AMBO_SITES[] = {
    {"dodecahedron_hk62", probe_dodeca_hk62},
    {"dodecahedron_hk35", probe_dodeca_hk35},
    {"dodecahedron_hk35_ambo_hk62", probe_dodeca_hk35_ambo_hk62},
    {"dodecahedron_hk54", probe_dodeca_hk54},
    {"octahedron_hk17", probe_octa_hk17},
    {"octahedron_hk34", probe_octa_hk34},
    {"rhombicuboctahedron_hk63", probe_rhombicubocta_hk63},
    {"truncatedIcosahedron_hk54", probe_ticosa_hk54},
};

/**
 * @brief Steps a truncate (ambo-leg) sweep on every hankin seed the Islamic
 *        recipes ambo, asserting constant raw and compiled face counts and a
 *        closed genus-0 manifold at every sampled parameter.
 */
inline void test_ambo_leg_on_hankin_seed_holds_topology() {
  constexpr int SAMPLES = 33;
  constexpr float T_HI = 0.5f - ConwayGraph::T_EPS_AMBO;

  for (const HankinAmboSite &site : HANKIN_AMBO_SITES) {
    const int failed_before = hs_test::stats().failed;

    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed;
    {
      constexpr size_t HALF = sizeof(morph_aux_buf) / 2;
      Arena ga(morph_aux_buf, HALF);
      Arena gb(morph_aux_buf + HALF, HALF);
      seed = Solids::finalize_solid(site.seed(ga, gb), persist);
    }

    size_t v0 = 0, f0 = 0, i0 = 0, compiled0 = 0;
    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    for (int s = 0; s < SAMPLES; ++s) {
      const float t =
          T_EPS + (T_HI - T_EPS) * (static_cast<float>(s) / (SAMPLES - 1));
      ScratchScope frame_a(a);
      ScratchScope frame_b(b);
      PolyMesh swept = MeshOps::truncate(seed, a, b, t);
      MeshState compiled;
      MeshOps::compile(swept, compiled, a, b);
      if (s == 0) {
        v0 = swept.vertices.size();
        f0 = swept.face_counts.size();
        i0 = swept.faces.size();
        compiled0 = compiled.face_counts.size();
        expect_op_counts(v0, f0, i0,
                         morph_op_counts(ConwayGraph::MorphOp::TRUNCATE, seed));
      } else {
        HS_EXPECT_EQ(swept.vertices.size(), v0);
        HS_EXPECT_EQ(swept.face_counts.size(), f0);
        HS_EXPECT_EQ(swept.faces.size(), i0);
        HS_EXPECT_EQ(compiled.face_counts.size(), compiled0);
      }
      check_face_counts_consistent(swept);
      check_indices_in_range(swept);
      check_all_unit_vertices(swept, 1e-3f);
      conway_tests::check_euler_genus0(swept);
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [hankin-ambo] %s failed (raw F=%zu, compiled F=%zu)\n",
                  site.name, f0, compiled0);
    else
      std::printf("  [hankin-ambo] %s: F=%zu compiled=%zu across %d samples\n",
                  site.name, f0, compiled0, SAMPLES);
  }
}

// ---------------------------------------------------------------------------
// Hankin-sweep probe (docs/opchain_morph_spec.md section 5.1): the four
// Phase-1 hankin legs re-run the one-shot MeshOps::hankin (the update-path
// geometry) per sampled angle; V/F/I and the compiled face count must not
// move from THETA_EPS to the recipe's arrival angle, and every sample must
// stay a closed genus-0 manifold.
// ---------------------------------------------------------------------------

inline uint8_t morph_bank_buf[64 * 1024]; /**< Baked palette LUT arena. */

/** @brief One Phase-1 hankin-sweep leg seed and its arrival angle. */
struct HankinSweepSite {
  const char *name;                     /**< Diagnostic label. */
  PolyMesh (*seed)(Arena &a, Arena &b); /**< Chain prefix up to the hankin. */
  float theta_star;                     /**< Arrival contact angle, radians. */
};

inline PolyMesh probe_dodecahedron(Arena &a, Arena &b) {
  return Solids::Platonic::dodecahedron(a, b);
}
inline PolyMesh probe_dodeca_hk62_ambo(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
      .hankin(62.0f * D2R)
      .ambo()
      .build();
}
inline PolyMesh probe_octahedron(Arena &a, Arena &b) {
  return Solids::Platonic::octahedron(a, b);
}
inline PolyMesh probe_octa_hk17_ambo(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::octahedron(a, b), a, b)
      .hankin(17.0f * D2R)
      .ambo()
      .build();
}

inline constexpr HankinSweepSite HANKIN_SWEEP_SITES[] = {
    {"dodecahedron", probe_dodecahedron,
     62.0f * Solids::IslamicStarPatterns::D2R},
    {"dodecahedron_hk62_ambo", probe_dodeca_hk62_ambo,
     62.0f * Solids::IslamicStarPatterns::D2R},
    {"octahedron", probe_octahedron, 17.0f * Solids::IslamicStarPatterns::D2R},
    {"octahedron_hk17_ambo", probe_octa_hk17_ambo,
     73.0f * Solids::IslamicStarPatterns::D2R},
};

/**
 * @brief Steps a hankin sweep on every Phase-1 hankin-leg seed, asserting
 *        constant raw and compiled face counts and a closed genus-0 manifold
 *        at every sampled angle from THETA_EPS to the arrival angle.
 */
inline void test_hankin_sweep_on_islamic_seeds_holds_topology() {
  constexpr int SAMPLES = 17;
  constexpr float THETA_EPS = Animation::OpLeg::THETA_EPS;

  for (const HankinSweepSite &site : HANKIN_SWEEP_SITES) {
    const int failed_before = hs_test::stats().failed;

    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed;
    {
      constexpr size_t HALF = sizeof(morph_aux_buf) / 2;
      Arena ga(morph_aux_buf, HALF);
      Arena gb(morph_aux_buf + HALF, HALF);
      seed = Solids::finalize_solid(site.seed(ga, gb), persist);
    }

    size_t v0 = 0, f0 = 0, i0 = 0, compiled0 = 0;
    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    for (int s = 0; s < SAMPLES; ++s) {
      const float theta =
          THETA_EPS + (site.theta_star - THETA_EPS) *
                          (static_cast<float>(s) / (SAMPLES - 1));
      ScratchScope frame_a(a);
      ScratchScope frame_b(b);
      PolyMesh swept = MeshOps::hankin(seed, a, b, theta);
      MeshState compiled;
      MeshOps::compile(swept, compiled, a, b);
      if (s == 0) {
        v0 = swept.vertices.size();
        f0 = swept.face_counts.size();
        i0 = swept.faces.size();
        compiled0 = compiled.face_counts.size();
        expect_op_counts(v0, f0, i0, hankin_op_counts(seed));
      } else {
        HS_EXPECT_EQ(swept.vertices.size(), v0);
        HS_EXPECT_EQ(swept.face_counts.size(), f0);
        HS_EXPECT_EQ(swept.faces.size(), i0);
        HS_EXPECT_EQ(compiled.face_counts.size(), compiled0);
      }
      check_face_counts_consistent(swept);
      check_indices_in_range(swept);
      check_all_unit_vertices(swept, 1e-3f);
      conway_tests::check_euler_genus0(swept);
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [hankin-sweep] %s failed (raw F=%zu, compiled F=%zu)\n",
                  site.name, f0, compiled0);
    else
      std::printf("  [hankin-sweep] %s: F=%zu compiled=%zu across %d samples\n",
                  site.name, f0, compiled0, SAMPLES);
  }
}

// ---------------------------------------------------------------------------
// Hankin-sweep stability probe: per-step branch, displacement and face-normal
// diagnostics for the four Phase-1 hankin legs. Compares the shipping
// parameterization (update_hankin re-solved per frame) against a
// slerp-from-corner parameterization (solve once at theta_star, then slerp
// each dynamic vertex out of its collapsed corner).
// ---------------------------------------------------------------------------

/** @brief Branch update_hankin took for one dynamic vertex. */
enum class HankinBranch : uint8_t {
  INTERSECT, /**< Contact-plane intersection accepted. */
  COLLAPSED, /**< is_flat or zero-length edge: snapped to the corner. */
  FALLBACK,  /**< Degenerate or far intersection: edge-midpoint mean. */
  BLENDED,   /**< Partial fallback: 0 < fallback_blend < 1. */
};

/** @brief One dynamic vertex's solved position and the branch that made it. */
struct HankinSolve {
  Vector pos;
  HankinBranch branch;
  /** dist^2(star, corner) / max(dist^2(m, corner)); this mirror takes the
   * fallback above STAR_FAR_RATIO_SQ. Zero on the non-intersect branches. */
  float far_ratio = 0;
};

/** @brief Mirror of update_hankin's local blend_above (core/mesh/hankin.h). */
inline float hankin_blend_above(float value, float start, float end) {
  const float t =
      std::max(0.0f, std::min(1.0f, (value - start) / (end - start)));
  return quintic_kernel(t);
}

/**
 * @brief Mirrors MeshOps::update_hankin's per-vertex solve, exposing the
 *        branch each dynamic vertex takes.
 * @param compiled Baked hankin topology.
 * @param angle Contact angle in radians.
 * @param out Filled with one entry per dynamic instruction.
 * @details Reproduces the regularization anchor, the plane_cross_sq parallel
 *   gate and the fallback blend; positions must track update_hankin exactly.
 */
inline void hankin_solve(const CompiledHankin &compiled, float angle,
                         std::vector<HankinSolve> &out) {
  const bool is_flat = std::abs(angle) < math::TOLERANCE;
  const float cos_ha = cosf(angle * 0.5f);
  const float sin_ha = sinf(angle * 0.5f);

  out.assign(compiled.dynamic_instructions.size(), HankinSolve{});
  for (size_t i = 0; i < compiled.dynamic_instructions.size(); ++i) {
    const HankinInstruction &instr = compiled.dynamic_instructions[i];
    const Vector p_corner = compiled.corner(instr.v_corner);
    const Vector cn = normalized_or(p_corner, p_corner);

    if (is_flat) {
      out[i] = {cn, HankinBranch::COLLAPSED};
      continue;
    }

    const Vector m1 = compiled.static_vertices[instr.idx_m1];
    const Vector m2 = compiled.static_vertices[instr.idx_m2];
    const Vector cross1 = cross(compiled.corner(instr.v_prev), p_corner);
    const Vector cross2 = cross(p_corner, compiled.corner(instr.v_next));
    if (dot(cross1, cross1) < math::EPS_CROSS_SQ ||
        dot(cross2, cross2) < math::EPS_CROSS_SQ) {
      out[i] = {cn, HankinBranch::COLLAPSED};
      continue;
    }

    const Quaternion q1(cos_ha, sin_ha * m1.x, sin_ha * m1.y, sin_ha * m1.z);
    const Quaternion q2(cos_ha, -sin_ha * m2.x, -sin_ha * m2.y, -sin_ha * m2.z);
    Vector intersect =
        cross(rotate(cross1.normalized(), q1), rotate(cross2.normalized(), q2));
    const float plane_cross_sq = dot(intersect, intersect);

    Vector fallback = normalized_or(m1 + m2, cn);
    if (dot(fallback, p_corner) < 0.0f)
      fallback = -fallback;
    const Vector oriented_intersect = intersect * dot(intersect, p_corner);
    const Vector raw_star = normalized_or(oriented_intersect, fallback);
    const float local_sq =
        std::max(distance_squared(m1, cn), distance_squared(m2, cn));
    if (!(local_sq > math::EPS_LEN_SQ)) {
      out[i] = {fallback, HankinBranch::FALLBACK, 0.0f};
      continue;
    }

    const float raw_ratio_sq = distance_squared(raw_star, cn) / local_sq;
    const float conditioned = hankin_blend_above(
        raw_ratio_sq, MeshOps::HANKIN_CONDITIONED_NEAR_RATIO_SQ,
        MeshOps::HANKIN_CONDITIONED_FAR_RATIO_SQ);
    const float anchor =
        std::max(0.0f,
                 MeshOps::HANKIN_PARALLEL_REGULARIZATION_SQ - plane_cross_sq) +
        conditioned * std::max(0.0f, MeshOps::HANKIN_CONDITIONED_CLEAR_SQ -
                                         plane_cross_sq);
    intersect = normalized_or(oriented_intersect + fallback * anchor, fallback);

    const float far_ratio = distance_squared(intersect, cn) / local_sq;
    const float parallel_gate = hankin_blend_above(
        -plane_cross_sq, -MeshOps::HANKIN_PARALLEL_GATE_HI_SQ,
        -MeshOps::HANKIN_PARALLEL_GATE_LO_SQ);
    const float fallback_blend =
        hankin_blend_above(far_ratio, MeshOps::STAR_FAR_BLEND_START_RATIO_SQ,
                           MeshOps::STAR_FAR_RATIO_SQ) *
        parallel_gate;
    if (fallback_blend >= 1.0f) {
      out[i] = {fallback, HankinBranch::FALLBACK, far_ratio};
      continue;
    }
    if (fallback_blend > 0.0f) {
      out[i] = {normalized_or(intersect * (1.0f - fallback_blend) +
                                  fallback * fallback_blend,
                              fallback),
                HankinBranch::BLENDED, far_ratio};
      continue;
    }
    out[i] = {intersect, HankinBranch::INTERSECT, far_ratio};
  }
}

/** Chord tolerance between hankin_solve and MeshOps::update_hankin; the two
 * evaluate identical expressions under independently contracted float
 * arithmetic, so bitwise agreement is not guaranteed. */
inline constexpr float HANKIN_MIRROR_TOL = 1e-5f;

/**
 * @brief Pins hankin_solve's star points to MeshOps::update_hankin's.
 * @param compiled Baked hankin topology.
 * @param angle Contact angle in radians.
 * @param arena Scratch arena backing the reference mesh; rewound on return.
 * @return Largest chord between a mirrored and a shipping star point.
 */
inline float hankin_check_mirror(CompiledHankin &compiled, float angle,
                                 Arena &arena) {
  ScratchScope guard(arena);
  std::vector<HankinSolve> mirror;
  hankin_solve(compiled, angle, mirror);
  PolyMesh ref;
  MeshOps::update_hankin(compiled, ref, arena, angle);
  const size_t base = compiled.static_vertices.size();
  HS_EXPECT_EQ(ref.vertices.size(), base + mirror.size());
  float max_chord = 0;
  for (size_t i = 0; i < mirror.size() && base + i < ref.vertices.size(); ++i)
    max_chord = std::max(max_chord,
                         (ref.vertices[base + i] - mirror[i].pos).magnitude());
  HS_EXPECT_LT(max_chord, HANKIN_MIRROR_TOL);
  return max_chord;
}

/**
 * @brief Newell normals of every compiled hankin face for a solved star set.
 * @param compiled Baked hankin topology.
 * @param dyn Solved dynamic vertices.
 * @param out Filled with one unnormalized Newell normal per face.
 */
inline void hankin_face_normals(const CompiledHankin &compiled,
                                const std::vector<HankinSolve> &dyn,
                                std::vector<Vector> &out) {
  auto vertex_at = [&](uint16_t idx) {
    return idx < compiled.static_offset ? compiled.static_vertices[idx]
                                        : dyn[idx - compiled.static_offset].pos;
  };
  out.assign(compiled.face_counts.size(), Vector());
  size_t base = 0;
  for (size_t f = 0; f < compiled.face_counts.size(); ++f) {
    const size_t n = compiled.face_counts[f];
    Vector normal;
    for (size_t k = 0; k < n; ++k) {
      const Vector a = vertex_at(compiled.faces[base + k]);
      const Vector b = vertex_at(compiled.faces[base + (k + 1) % n]);
      normal.x += (a.y - b.y) * (a.z + b.z);
      normal.y += (a.z - b.z) * (a.x + b.x);
      normal.z += (a.x - b.x) * (a.y + b.y);
    }
    out[f] = normal;
    base += n;
  }
}

/** Newell magnitude (twice the face area) below which a face's normal
 * direction is numerical noise and its sign is not evidence of a fold. */
inline constexpr float HANKIN_FLAT_FACE = 1e-4f;

/** @brief Per-step stability metrics of one sweep sample. */
struct HankinStepStats {
  float theta = 0; /**< Contact angle of this sample, radians. */
  int branch_flips =
      0; /**< Vertices whose branch changed vs the previous step. */
  float max_disp = 0;  /**< Largest single-vertex chord vs the previous step. */
  float mean_disp = 0; /**< Mean dynamic-vertex chord vs the previous step. */
  int normal_flips =
      0;              /**< Non-degenerate faces whose Newell normal reversed. */
  int flat_faces = 0; /**< Faces below HANKIN_FLAT_FACE this step. */
  float max_far_ratio = 0;    /**< Largest far_ratio this step. */
  float max_corner_chord = 0; /**< Largest chord(star point, its corner). */
};

/**
 * @brief Fills @p stats from consecutive solved states.
 * @param compiled Baked hankin topology.
 * @param prev Previous step's solve (empty for the first step).
 * @param prev_normals Previous step's face normals.
 * @param curr Current step's solve.
 * @param curr_normals Current step's face normals.
 * @param stats Metrics to fill; theta must already be set.
 */
inline void hankin_step_stats(const CompiledHankin &compiled,
                              const std::vector<HankinSolve> &prev,
                              const std::vector<Vector> &prev_normals,
                              const std::vector<HankinSolve> &curr,
                              const std::vector<Vector> &curr_normals,
                              HankinStepStats &stats) {
  for (size_t i = 0; i < curr.size(); ++i) {
    stats.max_far_ratio = std::max(stats.max_far_ratio, curr[i].far_ratio);
    const Vector cn = normalized_or(
        compiled.base_vertices[compiled.dynamic_instructions[i].v_corner],
        curr[i].pos);
    stats.max_corner_chord =
        std::max(stats.max_corner_chord, (curr[i].pos - cn).magnitude());
  }
  for (const Vector &n : curr_normals)
    if (n.magnitude() < HANKIN_FLAT_FACE)
      ++stats.flat_faces;
  if (prev.empty())
    return;
  double sum = 0;
  for (size_t i = 0; i < curr.size(); ++i) {
    if (curr[i].branch != prev[i].branch)
      ++stats.branch_flips;
    const float d = (curr[i].pos - prev[i].pos).magnitude();
    sum += d;
    stats.max_disp = std::max(stats.max_disp, d);
  }
  stats.mean_disp = static_cast<float>(sum / curr.size());
  for (size_t f = 0; f < curr_normals.size(); ++f) {
    if (prev_normals[f].magnitude() < HANKIN_FLAT_FACE ||
        curr_normals[f].magnitude() < HANKIN_FLAT_FACE)
      continue;
    if (dot(prev_normals[f], curr_normals[f]) < 0)
      ++stats.normal_flips;
  }
}

/** @brief Sweep-wide roll-up of the per-step metrics. */
struct HankinSweepSummary {
  int total_branch_flips = 0;
  int total_normal_flips = 0;
  int steps_with_normal_flips = 0;
  float worst_max_disp = 0;
  float worst_mean_disp = 0;
  int worst_step = 0;        /**< Step index owning worst_max_disp. */
  float spike_ratio = 0;     /**< worst_max_disp / mean_disp at that step. */
  int worst_flat_faces = 0;  /**< Most sub-HANKIN_FLAT_FACE faces in a step. */
  float worst_far_ratio = 0; /**< Largest far_ratio over the sweep. */
  float worst_corner_chord =
      0; /**< Largest chord(star point, corner) over the sweep. */
};

/**
 * @brief Rolls the per-step table into a summary.
 * @param table Per-step metrics, index 0 being the first (baseline) sample.
 * @return The roll-up.
 */
inline HankinSweepSummary
hankin_summarize(const std::vector<HankinStepStats> &table) {
  HankinSweepSummary sum;
  for (const HankinStepStats &r : table) {
    sum.worst_flat_faces = std::max(sum.worst_flat_faces, r.flat_faces);
    sum.worst_far_ratio = std::max(sum.worst_far_ratio, r.max_far_ratio);
    sum.worst_corner_chord =
        std::max(sum.worst_corner_chord, r.max_corner_chord);
  }
  for (size_t s = 1; s < table.size(); ++s) {
    sum.total_branch_flips += table[s].branch_flips;
    sum.total_normal_flips += table[s].normal_flips;
    if (table[s].normal_flips > 0)
      ++sum.steps_with_normal_flips;
    sum.worst_mean_disp = std::max(sum.worst_mean_disp, table[s].mean_disp);
    if (table[s].max_disp > sum.worst_max_disp) {
      sum.worst_max_disp = table[s].max_disp;
      sum.worst_step = static_cast<int>(s);
      sum.spike_ratio = table[s].mean_disp > 0
                            ? table[s].max_disp / table[s].mean_disp
                            : 0.0f;
    }
  }
  return sum;
}

/**
 * @brief Measures per-frame sweep stability of the four Phase-1 hankin legs
 *        under the shipping re-solve (uniform and eased theta) and under a
 *        slerp-from-corner parameterization.
 * @details Reports branch flips, per-step vertex chords and face-normal
 * reversals; asserts the slerp path is normal-flip free and that its endpoints
 * reproduce the collapsed form and the theta_star solve. Slerp carries the
 * arrival branch at every k, so its branch-flip count is zero by construction.
 */
inline void test_hankin_sweep_vertex_stability() {
  constexpr int SAMPLES = 32;
  constexpr float THETA_EPS = Animation::OpLeg::THETA_EPS;

  for (const HankinSweepSite &site : HANKIN_SWEEP_SITES) {
    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed;
    {
      constexpr size_t HALF = sizeof(morph_aux_buf) / 2;
      Arena ga(morph_aux_buf, HALF);
      Arena gb(morph_aux_buf + HALF, HALF);
      seed = Solids::finalize_solid(site.seed(ga, gb), persist);
    }

    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    CompiledHankin compiled;
    MeshOps::compile_hankin(seed, compiled, a, b);

    std::vector<HankinSolve> collapsed, arrival;
    hankin_solve(compiled, 0.0f, collapsed);
    hankin_solve(compiled, site.theta_star, arrival);

    // Every metric below reads hankin_solve, not the shipping solver. Pin the
    // two together over both sample grids first, else the whole suite measures
    // the mirror.
    float mirror_chord =
        std::max(hankin_check_mirror(compiled, 0.0f, b),
                 hankin_check_mirror(compiled, site.theta_star, b));
    for (int s = 0; s < SAMPLES; ++s) {
      const float u = static_cast<float>(s) / (SAMPLES - 1);
      const float span = site.theta_star - THETA_EPS;
      mirror_chord = std::max(
          {mirror_chord, hankin_check_mirror(compiled, THETA_EPS + span * u, b),
           hankin_check_mirror(compiled, THETA_EPS + span * ease_in_out_sin(u),
                               b)});
    }

    // Guard headroom: the far-star guard is scaled by the corner's local edge
    // scale, so on coarse seeds STAR_FAR_RATIO_SQ * local_sq can exceed the
    // 4.0 maximum squared chord, making the guard unreachable for that vertex.
    int unreachable = 0;
    float max_local_sq = 0;
    for (size_t i = 0; i < compiled.dynamic_instructions.size(); ++i) {
      const HankinInstruction &instr = compiled.dynamic_instructions[i];
      const Vector cn = normalized_or(compiled.base_vertices[instr.v_corner],
                                      compiled.base_vertices[instr.v_corner]);
      const float local_sq = std::max(
          distance_squared(compiled.static_vertices[instr.idx_m1], cn),
          distance_squared(compiled.static_vertices[instr.idx_m2], cn));
      max_local_sq = std::max(max_local_sq, local_sq);
      if (MeshOps::STAR_FAR_RATIO_SQ * local_sq >= 4.0f)
        ++unreachable;
    }

    std::printf("  [hankin-stability] %s: theta* = %.1f deg, base F=%zu, "
                "dyn V=%zu, hankin F=%zu, guard-unreachable %d/%zu "
                "(max local_sq %.4f)\n",
                site.name, site.theta_star * 180.0f / PI_F,
                seed.face_counts.size(), collapsed.size(),
                compiled.face_counts.size(), unreachable, collapsed.size(),
                max_local_sq);
    std::printf(
        "      mirror vs update_hankin: max chord %.3e over %d angles\n",
        mirror_chord, 2 * SAMPLES + 2);

    // Three parameterizations over the same sample grid.
    std::vector<HankinStepStats> tables[3];
    std::vector<std::vector<Vector>> resolve_pos(SAMPLES);
    float path_dev = 0;
    const char *mode_name[3] = {"resolve-uniform", "resolve-eased",
                                "slerp-eased    "};
    for (int mode = 0; mode < 3; ++mode) {
      std::vector<HankinSolve> prev, curr;
      std::vector<Vector> prev_normals, curr_normals;
      for (int s = 0; s < SAMPLES; ++s) {
        const float u = static_cast<float>(s) / (SAMPLES - 1);
        const float k = mode == 0 ? u : ease_in_out_sin(u);
        HankinStepStats row;
        if (mode == 2) {
          row.theta = site.theta_star;
          curr.assign(arrival.size(), HankinSolve{});
          for (size_t i = 0; i < arrival.size(); ++i) {
            curr[i] = {slerp(collapsed[i].pos, arrival[i].pos, k),
                       arrival[i].branch, 0.0f};
            path_dev = std::max(path_dev,
                                (curr[i].pos - resolve_pos[s][i]).magnitude());
          }
        } else {
          row.theta = THETA_EPS + (site.theta_star - THETA_EPS) * k;
          hankin_solve(compiled, row.theta, curr);
          if (mode == 1) {
            resolve_pos[s].resize(curr.size());
            for (size_t i = 0; i < curr.size(); ++i)
              resolve_pos[s][i] = curr[i].pos;
          }
        }
        hankin_face_normals(compiled, curr, curr_normals);
        hankin_step_stats(compiled, prev, prev_normals, curr, curr_normals,
                          row);
        tables[mode].push_back(row);
        prev = curr;
        prev_normals = curr_normals;
      }
    }

    for (int mode = 0; mode < 3; ++mode) {
      const HankinSweepSummary sum = hankin_summarize(tables[mode]);
      std::printf("      %s  branch_flips=%d normal_flips=%d (in %d steps) "
                  "worst_max=%.5f @s%d spike=%.1fx worst_mean=%.5f "
                  "flat<=%d far<=%.1f corner<=%.3f\n",
                  mode_name[mode], sum.total_branch_flips,
                  sum.total_normal_flips, sum.steps_with_normal_flips,
                  sum.worst_max_disp, sum.worst_step, sum.spike_ratio,
                  sum.worst_mean_disp, sum.worst_flat_faces,
                  sum.worst_far_ratio, sum.worst_corner_chord);
    }

    std::printf("      path deviation slerp vs resolve at equal k: max=%.5f\n",
                path_dev);

    // Slerp endpoints must reproduce the collapsed form and the theta_star
    // solve; the leg's closing bookend swap is only invisible if k=1 lands on
    // the arrival geometry.
    float end0 = 0, end1 = 0;
    size_t exact0 = 0, exact1 = 0;
    for (size_t i = 0; i < arrival.size(); ++i) {
      const Vector s0 = slerp(collapsed[i].pos, arrival[i].pos, 0.0f);
      const Vector s1 = slerp(collapsed[i].pos, arrival[i].pos, 1.0f);
      end0 = std::max(end0, (s0 - collapsed[i].pos).magnitude());
      end1 = std::max(end1, (s1 - arrival[i].pos).magnitude());
      exact0 += std::memcmp(&s0, &collapsed[i].pos, sizeof(Vector)) == 0;
      exact1 += std::memcmp(&s1, &arrival[i].pos, sizeof(Vector)) == 0;
    }
    std::printf("      endpoints: k=0 max_err=%.3e bitwise=%zu/%zu, "
                "k=1 max_err=%.3e bitwise=%zu/%zu\n",
                end0, exact0, arrival.size(), end1, exact1, arrival.size());
    HS_EXPECT_TRUE(end0 < 1e-5f);
    HS_EXPECT_TRUE(end1 < 1e-5f);

    const HankinSweepSummary slerp_sum = hankin_summarize(tables[2]);
    HS_EXPECT_EQ(slerp_sum.total_normal_flips, 0);

    // Opening bookend: chord between the collapsed form and the leg's first
    // drawn angle, in sphere radii (sub-pixel iff below ~1/display radius).
    float eps_chord = 0;
    std::vector<HankinSolve> at_eps;
    hankin_solve(compiled, THETA_EPS, at_eps);
    for (size_t i = 0; i < at_eps.size(); ++i)
      eps_chord =
          std::max(eps_chord, (at_eps[i].pos - collapsed[i].pos).magnitude());
    std::printf("      theta_eps=%.3f opening chord max=%.6f radii "
                "(%.3f px at r=64)\n",
                THETA_EPS, eps_chord, eps_chord * 64.0f);

    // The slerp path's k = 0 form is the exact collapse, where every rosette
    // face has zero area; a K_EPS floor is its THETA_EPS analog. The smallest
    // k that lifts every face above HANKIN_FLAT_FACE is the first probe on
    // every site, so a leg that opened its first drawn frame on a sub-area
    // face would need a k_eps above one step.
    constexpr int K_STEPS = 200;
    constexpr float K_EPS_BUDGET = 1.0f / K_STEPS;
    float k_eps = 1.0f;
    std::vector<HankinSolve> probe(arrival.size());
    std::vector<Vector> probe_normals;
    for (int q = 1; q <= K_STEPS; ++q) {
      const float k = static_cast<float>(q) / K_STEPS;
      for (size_t i = 0; i < arrival.size(); ++i)
        probe[i] = {slerp(collapsed[i].pos, arrival[i].pos, k),
                    arrival[i].branch, 0.0f};
      hankin_face_normals(compiled, probe, probe_normals);
      bool all_lit = true;
      for (const Vector &n : probe_normals)
        all_lit &= n.magnitude() >= HANKIN_FLAT_FACE;
      if (all_lit) {
        k_eps = k;
        break;
      }
    }
    std::printf("      slerp k_eps (first k with no sub-area face) = %.3f\n",
                k_eps);
    HS_EXPECT_LE(k_eps, K_EPS_BUDGET);
  }
}

/**
 * @brief Smoke-tests an OpLeg HANKIN_SWEEP end to end on the dodecahedron
 *        seed: construction populates the landing (star faces first, in
 *        base-face order), and every step hands the draw callback a compiled
 *        mesh with the constant hankin face count and in-range ramp indices.
 */
inline void test_opleg_hankin_sweep_smoke() {
  reset_globals();
  configure_arenas(GLOBAL_ARENA_SIZE - 24 * 1024 - 32 * 1024, 24 * 1024,
                   32 * 1024);
  hs::random().seed(2026u);

  Arena leg(morph_target_buf, sizeof(morph_target_buf));
  Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));

  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  PolyMesh dodeca;
  build_solid<Solids::Dodecahedron>(dodeca, leg);
  uint8_t pal[16];
  for (size_t f = 0; f < dodeca.face_counts.size(); ++f)
    pal[f] = static_cast<uint8_t>(f % Animation::OpLeg::PALETTES);

  Animation::OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                           .prev_face_palette = pal,
                                           .prev_faces =
                                               dodeca.face_counts.size()};

  // Per-frame motion bound: growing star points out from their corners keeps
  // every step small and unimodal. Re-solving the contact-plane intersection
  // per frame instead sends star points on geodesic excursions (measured to
  // 1.84 chord on ambo-of-hankin seeds), which draws as lines crossing the
  // pattern; the bound is what stops that parameterization coming back.
  constexpr float MAX_STEP_CHORD = 0.15f;
  LegDrawProbe probe;
  auto cb = [&](Canvas &, const MeshState &m,
                const Animation::OpLeg::Shading &sh) { probe.observe(m, sh); };

  using Solids::IslamicStarPatterns::D2R;
  constexpr int SWEEP = 8;
  Animation::OpLeg anim(dodeca,
                        Animation::OpLeg::HankinSweepSpec{
                            .theta_start = Animation::OpLeg::THETA_EPS,
                            .theta_end = 62.0f * D2R,
                            .sweep_frames = SWEEP},
                        leg, cb, handoff);

  // 12 star faces (the base-face-order prefix) + 20 rosette faces.
  const Animation::OpLeg::Landing &landing = anim.landing();
  HS_EXPECT_EQ(landing.primary_faces, dodeca.face_counts.size());
  HS_EXPECT_EQ(landing.faces,
               dodeca.face_counts.size() + dodeca.vertices.size());
  HS_EXPECT_TRUE(landing.topology != nullptr);

  hs_test::StubEffect fx(288, 144);
  for (int f = 0; f < SWEEP; ++f) {
    {
      Canvas c(fx);
      anim.step(c);
    }
    fx.advance_display();
  }
  HS_EXPECT_EQ(probe.drawn, (size_t)SWEEP);
  HS_EXPECT_EQ(probe.faces, landing.faces);
  HS_EXPECT_LT(probe.worst_step, MAX_STEP_CHORD);
  std::printf("  [opleg hankin] worst per-frame vertex step %.4f chord "
              "(bound %.2f)\n",
              (double)probe.worst_step, (double)MAX_STEP_CHORD);
}

// ---------------------------------------------------------------------------
// Recipe-step leg kinds (docs/opchain_morph_spec.md sections 2.1/2.2): the
// truncate, snub and relax legs the pure-inflate recipes need. Topology
// constancy is sampled on the chain prefixes the shipping recipes actually
// reach; the smoke tests drive a whole leg through OpLeg on the same seeds.
// ---------------------------------------------------------------------------

/** @brief One recipe-step leg site: the chain prefix the step sweeps on. */
struct StepLegSite {
  const char *name;                     /**< Diagnostic label. */
  PolyMesh (*seed)(Arena &a, Arena &b); /**< Chain prefix up to the step. */
  float param; /**< Arrival t, or relax iteration count. */
};

inline PolyMesh probe_icosahedron(Arena &a, Arena &b) {
  return Solids::Platonic::icosahedron(a, b);
}
inline PolyMesh probe_icosa_ambo(Arena &a, Arena &b) {
  return Solids::SolidBuilder(Solids::Platonic::icosahedron(a, b), a, b)
      .ambo()
      .build();
}
inline PolyMesh probe_icosa_snub(Arena &a, Arena &b) {
  return Solids::SolidBuilder(Solids::Platonic::icosahedron(a, b), a, b)
      .snub()
      .build();
}
inline PolyMesh probe_icosa_snub_relax(Arena &a, Arena &b) {
  return Solids::SolidBuilder(Solids::Platonic::icosahedron(a, b), a, b)
      .snub()
      .relax()
      .build();
}
inline PolyMesh probe_ticosa_ambo(Arena &a, Arena &b) {
  return Solids::SolidBuilder(Solids::Archimedean::truncatedIcosahedron(a, b),
                              a, b)
      .ambo()
      .build();
}
inline PolyMesh probe_ticosa_ambo_relax217(Arena &a, Arena &b) {
  return Solids::SolidBuilder(Solids::Archimedean::truncatedIcosahedron(a, b),
                              a, b)
      .ambo()
      .relax(217)
      .build();
}
inline PolyMesh probe_dodeca_ambo_bevel33(Arena &a, Arena &b) {
  return Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
      .ambo()
      .bevel(0.33f)
      .build();
}

/** Truncate-leg sites: the three pure-inflate recipes truncating at 0.33. */
inline constexpr StepLegSite TRUNCATE_LEG_SITES[] = {
    {"icosahedron_ambo", probe_icosa_ambo, 0.33f},
    {"truncatedIcosahedron_ambo_relax217", probe_ticosa_ambo_relax217, 0.33f},
    {"icosahedron_snub_relax", probe_icosa_snub_relax, 0.33f},
};

/** Snub-leg sites: the one recipe snubbing, at the .snub() defaults. */
inline constexpr StepLegSite SNUB_LEG_SITES[] = {
    {"icosahedron", probe_icosahedron, 0.5f},
};

/** Standalone relax-leg sites, spanning short to long iteration counts. */
inline constexpr StepLegSite RELAX_LEG_SITES[] = {
    {"icosahedron_snub", probe_icosa_snub, 8.0f},
    {"dodecahedron_ambo_bevel33", probe_dodeca_ambo_bevel33, 100.0f},
    {"truncatedIcosahedron_ambo", probe_ticosa_ambo, 217.0f},
};

/** @brief Topology fingerprint of one sweep sample. */
struct SweepFingerprint {
  size_t v = 0;        /**< Raw vertex count. */
  size_t f = 0;        /**< Raw face count. */
  size_t i = 0;        /**< Raw face-index count. */
  size_t compiled = 0; /**< Compiled face count. */
};

/**
 * @brief Fingerprints a swept mesh and runs the structural checks.
 * @param swept Sweep sample.
 * @param a Arena receiving the compiled mesh.
 * @param b Compile scratch arena.
 * @return The sample's fingerprint.
 */
inline SweepFingerprint check_sweep_sample(const PolyMesh &swept, Arena &a,
                                           Arena &b) {
  MeshState compiled;
  MeshOps::compile(swept, compiled, a, b);
  check_face_counts_consistent(swept);
  check_indices_in_range(swept);
  check_all_unit_vertices(swept, 1e-3f);
  conway_tests::check_euler_genus0(swept);
  return {swept.vertices.size(), swept.face_counts.size(), swept.faces.size(),
          compiled.face_counts.size()};
}

/**
 * @brief Asserts a sweep sample fingerprint matches the opening one.
 * @param s Sample fingerprint.
 * @param first Opening-sample fingerprint.
 */
inline void expect_same_fingerprint(const SweepFingerprint &s,
                                    const SweepFingerprint &first) {
  HS_EXPECT_EQ(s.v, first.v);
  HS_EXPECT_EQ(s.f, first.f);
  HS_EXPECT_EQ(s.i, first.i);
  HS_EXPECT_EQ(s.compiled, first.compiled);
}

/**
 * @brief Builds a leg site's seed into a persistent arena.
 * @param site Site whose chain prefix is built.
 * @param persist Arena receiving the finalized seed.
 * @return The seed mesh in @p persist.
 */
inline PolyMesh build_step_leg_seed(const StepLegSite &site, Arena &persist) {
  constexpr size_t HALF = sizeof(morph_aux_buf) / 2;
  Arena ga(morph_aux_buf, HALF);
  Arena gb(morph_aux_buf + HALF, HALF);
  return Solids::finalize_solid(site.seed(ga, gb), persist);
}

/**
 * @brief Steps a truncate sweep on every seed the recipes truncate, asserting
 *        constant raw and compiled face counts and a closed genus-0 manifold
 *        at every sampled parameter.
 */
inline void test_truncate_leg_on_recipe_seeds_holds_topology() {
  constexpr int SAMPLES = 33;

  for (const StepLegSite &site : TRUNCATE_LEG_SITES) {
    const int failed_before = hs_test::stats().failed;
    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed = build_step_leg_seed(site, persist);

    SweepFingerprint first;
    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    for (int s = 0; s < SAMPLES; ++s) {
      const float t = T_EPS + (site.param - T_EPS) *
                                  (static_cast<float>(s) / (SAMPLES - 1));
      ScratchScope frame_a(a);
      ScratchScope frame_b(b);
      const SweepFingerprint fp =
          check_sweep_sample(MeshOps::truncate(seed, a, b, t), a, b);
      if (s == 0) {
        first = fp;
        expect_op_counts(fp.v, fp.f, fp.i,
                         morph_op_counts(ConwayGraph::MorphOp::TRUNCATE, seed));
      } else {
        expect_same_fingerprint(fp, first);
      }
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [truncate-leg] %s failed (raw F=%zu, compiled F=%zu)\n",
                  site.name, first.f, first.compiled);
    else
      std::printf("  [truncate-leg] %s: t*=%.2f F=%zu compiled=%zu across %d "
                  "samples\n",
                  site.name, (double)site.param, first.f, first.compiled,
                  SAMPLES);
  }
}

/**
 * @brief Steps a snub sweep on every seed the recipes snub, asserting constant
 *        raw and compiled face counts and a closed genus-0 manifold at every
 *        sampled parameter.
 */
inline void test_snub_leg_on_recipe_seeds_holds_topology() {
  constexpr int SAMPLES = 33;

  for (const StepLegSite &site : SNUB_LEG_SITES) {
    const int failed_before = hs_test::stats().failed;
    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed = build_step_leg_seed(site, persist);

    SweepFingerprint first;
    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    for (int s = 0; s < SAMPLES; ++s) {
      const float t = T_EPS + (site.param - T_EPS) *
                                  (static_cast<float>(s) / (SAMPLES - 1));
      ScratchScope frame_a(a);
      ScratchScope frame_b(b);
      const SweepFingerprint fp =
          check_sweep_sample(MeshOps::snub(seed, a, b, t, 0.0f), a, b);
      if (s == 0) {
        first = fp;
        expect_op_counts(fp.v, fp.f, fp.i,
                         morph_op_counts(ConwayGraph::MorphOp::SNUB, seed));
      } else {
        expect_same_fingerprint(fp, first);
      }
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [snub-leg] %s failed (raw F=%zu, compiled F=%zu)\n",
                  site.name, first.f, first.compiled);
    else
      std::printf("  [snub-leg] %s: t*=%.2f F=%zu compiled=%zu across %d "
                  "samples\n",
                  site.name, (double)site.param, first.f, first.compiled,
                  SAMPLES);
  }
}

/**
 * @brief Steps a relax slerp on every seed the recipes relax standalone,
 *        asserting the vertex count/order identity the kind rests on plus
 *        constant raw and compiled face counts across the slerp.
 */
inline void test_relax_leg_on_recipe_seeds_holds_topology() {
  constexpr int SAMPLES = 33;

  for (const StepLegSite &site : RELAX_LEG_SITES) {
    const int failed_before = hs_test::stats().failed;
    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed = build_step_leg_seed(site, persist);

    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    PolyMesh relaxed = MeshOps::relax(seed, a, b, static_cast<int>(site.param));

    // The precondition of a standalone relax leg: same vertex count, same
    // topology bytes, and vertex i still nearest its own seed vertex, so the
    // leg slerps per-vertex with no correspondence pass.
    HS_EXPECT_EQ(relaxed.vertices.size(), seed.vertices.size());
    HS_EXPECT_EQ(relaxed.face_counts.size(), seed.face_counts.size());
    HS_EXPECT_EQ(relaxed.faces.size(), seed.faces.size());
    HS_EXPECT_EQ(std::memcmp(relaxed.face_counts.data(),
                             seed.face_counts.data(),
                             seed.face_counts.size() * sizeof(uint8_t)),
                 0);
    HS_EXPECT_EQ(std::memcmp(relaxed.faces.data(), seed.faces.data(),
                             seed.faces.size() * sizeof(uint16_t)),
                 0);
    for (size_t i = 0; i < relaxed.vertices.size(); ++i) {
      size_t nearest = 0;
      float best = 1e9f;
      for (size_t j = 0; j < seed.vertices.size(); ++j) {
        const float d = distance_between(relaxed.vertices[i], seed.vertices[j]);
        if (d < best) {
          best = d;
          nearest = j;
        }
      }
      HS_EXPECT_EQ(nearest, i);
    }

    SweepFingerprint first;
    for (int s = 0; s < SAMPLES; ++s) {
      const float k = static_cast<float>(s) / (SAMPLES - 1);
      ScratchScope frame_a(a);
      ScratchScope frame_b(b);
      PolyMesh swept;
      MeshOps::clone(seed, swept, a);
      for (size_t i = 0; i < swept.vertices.size(); ++i)
        swept.vertices[i] = slerp(seed.vertices[i], relaxed.vertices[i], k);
      const SweepFingerprint fp = check_sweep_sample(swept, a, b);
      if (s == 0) {
        first = fp;
        // A relax slerp moves vertices only: the seed's own counts are the
        // expectation at every k.
        expect_op_counts(fp.v, fp.f, fp.i, mesh_op_counts(seed));
      } else {
        expect_same_fingerprint(fp, first);
      }
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [relax-leg] %s failed (raw F=%zu, compiled F=%zu)\n",
                  site.name, first.f, first.compiled);
    else
      std::printf("  [relax-leg] %s: iters=%d V=%zu F=%zu compiled=%zu across "
                  "%d samples\n",
                  site.name, (int)site.param, first.v, first.f, first.compiled,
                  SAMPLES);
  }
}

// ---------------------------------------------------------------------------
// Medial (Conway-dual) bridge: MeshOps::medial produces one shared rectified
// connectivity with both endpoint vertex sets a_e (== ambo(P)) and b_e
// (== ambo(dual(P))). The smooth dual replaces the instant partition swap with
// a truncate sweep to ambo(P), a slerp of every medial vertex a_e -> b_e, and a
// truncate sweep down to dual(P). These gate the medial leg (the new slerp
// segment) on the real DUAL-leg seeds: the endpoints match ambo(P)/ambo(dual(P))
// to tolerance (the correspondence proof), and across the slerp there are no
// face inversions, no self-intersection (total solid angle holds at 4pi), no
// degenerate collapse, no antipodal slerp inputs, and the emission order is
// fixed frame to frame.
// ---------------------------------------------------------------------------

inline PolyMesh probe_icosa_kis_snub(Arena &a, Arena &b) {
  return Solids::SolidBuilder(Solids::Platonic::icosahedron(a, b), a, b)
      .kis()
      .snub()
      .build();
}
inline PolyMesh probe_toct_snub(Arena &a, Arena &b) {
  return Solids::SolidBuilder(Solids::Archimedean::truncatedOctahedron(a, b), a,
                              b)
      .snub()
      .build();
}
inline PolyMesh probe_dodeca_hk72_ambo(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
      .hankin(72.0f * D2R)
      .ambo()
      .build();
}
inline PolyMesh probe_icosidodeca_trunc5_ambo(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::TRUNCATE_T_NEAR;
  return Solids::SolidBuilder(Solids::Archimedean::icosidodecahedron(a, b), a,
                              b)
      .truncate(TRUNCATE_T_NEAR)
      .ambo()
      .build();
}
/** DUAL-leg sites: the mesh each recipe applies a smooth dual to (the gyro
 * snub-derived seeds, the ambo-of-hankin and ambo-of-truncate seeds, and the
 * needle's hankin seed). */
inline constexpr StepLegSite DUAL_LEG_SITES[] = {
    {"icosahedron_kis_gyro", probe_icosa_kis_snub, 0.0f},
    {"truncatedOctahedron_gyro", probe_toct_snub, 0.0f},
    {"dodecahedron_hk72_ambo_dual", probe_dodeca_hk72_ambo, 0.0f},
    {"icosidodecahedron_truncate5d_ambo_dual", probe_icosidodeca_trunc5_ambo,
     0.0f},
    {"truncatedIcosahedron_ambo_relax100_hk54_needle",
     build_ticosa_ambo_relax100_hk54, 0.0f},
};

/** @brief Max nearest-vertex distance from every vertex of @p x to @p y. */
inline float medial_vertex_set_dist(const PolyMesh &x, const PolyMesh &y) {
  float worst = 0.0f;
  for (const auto &vx : x.vertices) {
    float best = 1e9f;
    for (const auto &vy : y.vertices)
      best = std::min(best, distance_between(vx, vy));
    worst = std::max(worst, best);
  }
  return worst;
}
inline float medial_vertex_set_dist(const ArenaVector<Vector> &x,
                                    const PolyMesh &y) {
  float worst = 0.0f;
  for (const auto &vx : x) {
    float best = 1e9f;
    for (const auto &vy : y.vertices)
      best = std::min(best, distance_between(vx, vy));
    worst = std::max(worst, best);
  }
  return worst;
}

/** @brief Signed total solid angle of a mesh (per-face fan from its centroid).
 * A simple closed surface sums to 4pi; a self-intersection breaks it. */
inline double medial_total_solid_angle(const PolyMesh &m) {
  double total = 0.0;
  size_t off = 0;
  for (size_t f = 0; f < m.face_counts.size(); ++f) {
    const int n = m.face_counts[f];
    Vector c(0, 0, 0);
    for (int k = 0; k < n; ++k)
      c = c + m.vertices[m.faces[off + k]];
    c = c.normalized();
    for (int k = 0; k < n; ++k) {
      const Vector a = m.vertices[m.faces[off + k]];
      const Vector b = m.vertices[m.faces[off + (k + 1) % n]];
      const double num = dot(c, cross(a, b));
      const double den = 1.0 + dot(c, a) + dot(a, b) + dot(b, c);
      total += 2.0 * std::atan2(num, den);
    }
    off += n;
  }
  return total;
}

/** @brief Smallest spherical-triangle-fan area over all faces. */
inline double medial_min_face_area(const PolyMesh &m) {
  double mn = 1e9;
  size_t off = 0;
  for (size_t f = 0; f < m.face_counts.size(); ++f) {
    const int n = m.face_counts[f];
    Vector c(0, 0, 0);
    for (int k = 0; k < n; ++k)
      c = c + m.vertices[m.faces[off + k]];
    c = c.normalized();
    double area = 0.0;
    for (int k = 0; k < n; ++k) {
      const Vector a = m.vertices[m.faces[off + k]];
      const Vector b = m.vertices[m.faces[off + (k + 1) % n]];
      area += 0.5 * cross(a - c, b - c).length();
    }
    mn = std::min(mn, area);
    off += n;
  }
  return mn;
}

/** @brief Faces whose area-weighted normal points inward (a fold/inversion). */
inline int medial_inverted_faces(const PolyMesh &m) {
  int inv = 0;
  size_t off = 0;
  for (size_t f = 0; f < m.face_counts.size(); ++f) {
    const int n = m.face_counts[f];
    Vector c(0, 0, 0);
    for (int k = 0; k < n; ++k)
      c = c + m.vertices[m.faces[off + k]];
    c = c.normalized();
    Vector nrm(0, 0, 0);
    for (int k = 0; k < n; ++k) {
      const Vector a = m.vertices[m.faces[off + k]];
      const Vector b = m.vertices[m.faces[off + (k + 1) % n]];
      nrm = nrm + cross(a, b);
    }
    if (dot(nrm, c) < 0.0f)
      ++inv;
    off += n;
  }
  return inv;
}

/**
 * @brief Gates the medial dual bridge on every DUAL-leg seed: endpoint match
 *        plus slerp well-formedness (task-validated failure modes).
 */
inline void test_medial_dual_bridge_wellformed() {
  constexpr int SAMPLES = 33;
  // Below any real medial face; snub-derived seeds bottom out ~4e-2, all such
  // minima at the endpoints (never manufactured mid-slerp).
  constexpr double MIN_FACE_AREA = 1e-3;
  // Well clear of an antipodal/coincident slerp singularity.
  constexpr float MIN_ENDPOINT_DOT = 0.9f;
  constexpr float MAX_INTER_SAMPLE_STEP = 0.05f;
  constexpr float ENDPOINT_TOL = 1e-4f;

  for (const StepLegSite &site : DUAL_LEG_SITES) {
    const int failed_before = hs_test::stats().failed;
    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh P = build_step_leg_seed(site, persist);

    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    Arena aux(morph_aux_buf, sizeof(morph_aux_buf));

    PolyMesh med_a;
    ArenaVector<Vector> med_b;
    MeshOps::medial(P, med_a, med_b, a, b);

    // Correspondence proof: s = 0 is ambo(P), s = 1 (b_e positions) is
    // ambo(dual(P)). Both are the same rectified polyhedron (one face per
    // primal face + one per primal vertex), so the FACE count is identical
    // through the bridge and leg 3 sweeps truncate(dual(P)) on that same
    // connectivity. Vertex identity is checked by position, not count: a
    // hankin seed has degree-2 star-point vertices, so its dual is lossy
    // (digon faces drop) and MeshOps::ambo(dual(P)) merges the coincident edge
    // midpoints — leg 3's truncate(dual, 0.5-eps) keeps them apart, matching
    // the medial's 2E vertices exactly, so the bridge stays continuous.
    PolyMesh ambo_p = MeshOps::ambo(P, b, aux);
    PolyMesh dual_p = MeshOps::dual(P, b, aux);
    PolyMesh ambo_dual_p = MeshOps::ambo(dual_p, aux, b);
    HS_EXPECT_EQ(med_a.vertices.size(), ambo_p.vertices.size());
    HS_EXPECT_EQ(med_a.face_counts.size(), ambo_p.face_counts.size());
    HS_EXPECT_EQ(med_a.face_counts.size(), ambo_dual_p.face_counts.size());
    HS_EXPECT_TRUE(std::equal(med_a.face_counts.begin(),
                              med_a.face_counts.end(),
                              ambo_p.face_counts.begin()));
    HS_EXPECT_TRUE(std::equal(med_a.faces.begin(), med_a.faces.end(),
                              ambo_p.faces.begin()));

    // MeshOps::medial documents out_a as bit-identical to ambo(P), and the
    // bridge's ambo(P)->medial seam depends on it. Pin exact float equality
    // index-for-index, before the snorm16 packing below can absorb a drift.
    const size_t pair_n =
        std::min(med_a.vertices.size(), ambo_p.vertices.size());
    size_t bit_diff = 0;
    float worst_bit = 0.0f;
    for (size_t v = 0; v < pair_n; ++v) {
      const Vector &mv = med_a.vertices[v];
      const Vector &av = ambo_p.vertices[v];
      if (mv.x != av.x || mv.y != av.y || mv.z != av.z) {
        ++bit_diff;
        worst_bit = std::max(worst_bit, distance_between(mv, av));
      }
    }
    HS_EXPECT_EQ(bit_diff, size_t(0));
    if (bit_diff != 0)
      std::printf("    [medial] %s out_a != ambo(P) exactly: %zu/%zu vertices, "
                  "worst |d|=%.3e\n",
                  site.name, bit_diff, pair_n, (double)worst_bit);

    // The MEDIAL_SLERP leg stores both endpoint sets snorm16-packed, so gate the
    // quantized-then-decoded positions the leg actually slerps, not the
    // full-precision medial output.
    for (auto &v : med_a.vertices)
      v = Snorm3::encode(v).decode().normalized();
    for (auto &v : med_b)
      v = Snorm3::encode(v).decode().normalized();

    // Every medial vertex sits on an ambo(P) vertex at s=0 and an ambo(dual(P))
    // vertex at s=1 (set containment; a lossy dual makes s=1 many-to-one).
    const float end_a = medial_vertex_set_dist(med_a, ambo_p);
    const float end_b = medial_vertex_set_dist(med_b, ambo_dual_p);
    HS_EXPECT_LT(end_a, ENDPOINT_TOL);
    HS_EXPECT_LT(end_b, ENDPOINT_TOL);

    // No antipodal/coincident slerp inputs.
    float min_dot = 2.0f;
    for (size_t v = 0; v < med_a.vertices.size(); ++v)
      min_dot = std::min(
          min_dot, dot(med_a.vertices[v].normalized(), med_b[v].normalized()));
    HS_EXPECT_GT(min_dot, MIN_ENDPOINT_DOT);

    // Slerp sweep: fixed connectivity, per-vertex slerp; well-formed throughout.
    int total_inv = 0;
    double worst_4pi = 0.0, min_area = 1e9;
    float max_step = 0.0f;
    SweepFingerprint first;
    std::vector<Vector> prev(med_a.vertices.size());
    for (int s = 0; s < SAMPLES; ++s) {
      const float k = static_cast<float>(s) / (SAMPLES - 1);
      ScratchScope frame_a(aux);
      PolyMesh frame;
      frame.vertices.bind(aux, med_a.vertices.size());
      for (size_t v = 0; v < med_a.vertices.size(); ++v)
        frame.vertices.push_back(slerp(med_a.vertices[v], med_b[v], k));
      frame.face_counts.bind(aux, med_a.face_counts.size());
      frame.face_counts.append_bulk(med_a.face_counts.data(),
                                    med_a.face_counts.size());
      frame.faces.bind(aux, med_a.faces.size());
      frame.faces.append_bulk(med_a.faces.data(), med_a.faces.size());

      {
        ScratchScope ca(a);
        ScratchScope cb(b);
        const SweepFingerprint fp = check_sweep_sample(frame, a, b);
        if (s == 0)
          first = fp;
        else
          expect_same_fingerprint(fp, first); // fixed emission order/count
      }
      total_inv += medial_inverted_faces(frame);
      worst_4pi = std::max(worst_4pi, std::abs(medial_total_solid_angle(frame) -
                                               4.0 * 3.14159265358979323846));
      min_area = std::min(min_area, medial_min_face_area(frame));
      if (s > 0)
        for (size_t v = 0; v < frame.vertices.size(); ++v)
          max_step =
              std::max(max_step, distance_between(frame.vertices[v], prev[v]));
      prev.assign(frame.vertices.begin(), frame.vertices.end());
    }

    HS_EXPECT_EQ(total_inv, 0);
    HS_EXPECT_LT(worst_4pi, 1e-5);
    HS_EXPECT_GT(min_area, MIN_FACE_AREA);
    HS_EXPECT_LT(max_step, MAX_INTER_SAMPLE_STEP);

    if (hs_test::stats().failed != failed_before)
      std::printf("    [medial] %s FAILED (endA=%.2e endB=%.2e min a.b=%.4f "
                  "inv=%d 4pi-err=%.2e minA=%.3e step=%.4f)\n",
                  site.name, (double)end_a, (double)end_b, (double)min_dot,
                  total_inv, worst_4pi, min_area, (double)max_step);
    else
      std::printf(
          "  [medial] %s: F=%zu endpoints<%.0e min a.b=%.4f 4pi-err=%.1e "
          "minA=%.2e step=%.4f\n",
          site.name, med_a.face_counts.size(), (double)ENDPOINT_TOL,
          (double)min_dot, worst_4pi, min_area, (double)max_step);
  }
}

/**
 * @brief Drives a MEDIAL_SLERP leg to completion through OpLeg on every DUAL-leg
 *        seed: constant compiled face count, bounded per-frame motion, a total
 *        palette mapping, and a landing describing the whole medial face list.
 * @details The leg departs from ambo(P) (one face per primal face + one per
 * primal vertex), so the handoff is class-keyed on that mesh with its face
 * centroids for the geometric provenance mapping.
 */
inline void test_opleg_medial_leg_smoke() {
  using Animation::OpLeg;
  reset_globals();
  configure_arenas(GLOBAL_ARENA_SIZE - 116 * 1024 - 74 * 1024, 116 * 1024,
                   74 * 1024);
  hs::random().seed(2026u);

  Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  constexpr int SWEEP = 24;
  constexpr float MAX_STEP_CHORD = 0.15f;

  for (const StepLegSite &site : DUAL_LEG_SITES) {
    const int failed_before = hs_test::stats().failed;
    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh P = build_step_leg_seed(site, persist);

    Arena leg(morph_target_buf, sizeof(morph_target_buf));
    Arena temp(morph_temp_buf, sizeof(morph_temp_buf));

    // Departed mesh: ambo(P), class-keyed palettes plus face centroids.
    PolyMesh ambo_p = Solids::finalize_solid(MeshOps::ambo(P, leg, temp), leg);
    {
      ScratchScope ta(scratch_arena_a);
      ScratchScope tb(scratch_arena_b);
      MeshOps::classify_faces_by_topology(ambo_p, scratch_arena_a,
                                          scratch_arena_b, leg);
    }
    const size_t prev_faces = ambo_p.face_counts.size();
    std::vector<uint8_t> pal(prev_faces);
    std::vector<Vector> centroid(prev_faces);
    size_t off = 0;
    for (size_t f = 0; f < prev_faces; ++f) {
      pal[f] = static_cast<uint8_t>(
          wrap(static_cast<int>(ambo_p.topology[f]), OpLeg::PALETTES));
      const int n = ambo_p.face_counts[f];
      centroid[f] = face_centroid_unit(ambo_p, off, n);
      off += n;
    }

    OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                  .prev_face_palette = pal.data(),
                                  .prev_faces = prev_faces,
                                  .prev_face_centroid = centroid.data()};

    LegDrawProbe probe;
    auto cb = [&](Canvas &, const MeshState &m, const OpLeg::Shading &sh) {
      probe.observe(m, sh);
    };

    OpLeg anim(P, OpLeg::MedialSpec{.sweep_frames = SWEEP}, leg, cb, handoff);
    const OpLeg::Landing &landing = anim.landing();
    // The medial connectivity is ambo(P), so the leg's face list is the whole
    // rectified polyhedron and every face lives the whole slerp.
    HS_EXPECT_EQ(landing.faces, prev_faces);
    HS_EXPECT_EQ(landing.primary_faces, prev_faces);

    hs_test::StubEffect fx(288, 144);
    for (int f = 0; f < SWEEP; ++f) {
      {
        Canvas c(fx);
        anim.step(c);
      }
      fx.advance_display();
    }
    HS_EXPECT_EQ(probe.drawn, (size_t)SWEEP);
    HS_EXPECT_EQ(landing.faces, probe.faces);
    HS_EXPECT_LT(probe.worst_step, MAX_STEP_CHORD);
    HS_EXPECT_TRUE(landing.from_palette != nullptr);

    std::printf("  [opleg medial] %s: F=%zu across %d frames, worst step "
                "%.4f%s\n",
                site.name, landing.faces, SWEEP, (double)probe.worst_step,
                hs_test::stats().failed != failed_before ? " FAILED" : "");
  }
}

/**
 * @brief Drives the dual bridge's medial -> closing-leg handoff on every
 *        DUAL-leg seed and pins the face correspondence across the seam.
 * @details The closing leg's face list is block-transposed against the
 * medial's ([D-faces][D-vertex orbits] vs [P-faces][P-vertex orbits]): the
 * k-th emitted P-vertex orbit is the k-th dual face, and each P-face is a
 * dual vertex whose orbit face lands somewhere in the trailing block. Both
 * sides coincide geometrically at the ambo point, so the probe derives the
 * permutation by exact centroid matching there and requires the closing
 * leg's from-palettes to follow it. A rendered A/B (leg-2 last frame vs
 * leg-3 first frame, real ramp LUTs, no camera) then gates the seam
 * at pixel level: the diff must stay near one in-leg step, not a flip. The
 * needle site runs on truncate(X, 1/3) -- the dt-macro bridge seed whose
 * valence-2 star tips make the blocks unequal (362 vs 720) -- on the
 * bridge arena split, so the closing leg's construction peak is covered.
 */
inline void test_opleg_dual_bridge_seam_correspondence() {
  using Animation::OpLeg;
  constexpr int SWEEP = 24;
  constexpr int RW = 288, RH = 144;
  constexpr float SEAM_MATCH_TOL = 0.02f;
  // A continuous seam leaves a core of pixels the palette swap never touches;
  // a mis-keyed one repaints nearly the whole frame. Both that core and the
  // in-leg control's scale together with the shading and the camera, so gate
  // their ratio rather than an absolute share of the frame: measured minimum
  // across sites 0.213 (needle), against 0.133 for a mis-keyed seam.
  constexpr float SEAM_QUIET_RATIO = 0.15f;
  // Measured maximum across sites 6.14M; a mis-keyed needle reaches 8.4M.
  constexpr long long SEAM_SUMABS_MAX = 7300000ll;

  static Pipeline<RW, RH> filters;

  for (const StepLegSite &site : DUAL_LEG_SITES) {
    const int failed_before = hs_test::stats().failed;
    reset_globals();
    configure_arenas(GLOBAL_ARENA_SIZE - 132 * 1024 - 74 * 1024, 132 * 1024,
                     74 * 1024);
    hs::random().seed(2026u);

    Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
    MeshPaletteBank bank;
    bank.bake_all(bank_arena);

    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    Arena leg(morph_target_buf, sizeof(morph_target_buf));
    Arena temp(morph_temp_buf, sizeof(morph_temp_buf));

    // The needle reaches its bridge through the dt macro, so its seam runs on
    // truncate(X, 1/3); every other site duals its recipe mesh directly.
    PolyMesh P = build_step_leg_seed(site, persist);
    if (std::strstr(site.name, "needle")) {
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      P = Solids::finalize_solid(MeshOps::truncate(P, aux, temp, 1.0f / 3.0f),
                                 leg);
    }
    const size_t PF = P.face_counts.size();

    // Departed mesh ambo(P): class-keyed palettes in medial face order.
    PolyMesh ambo_p;
    {
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      ambo_p = Solids::finalize_solid(MeshOps::ambo(P, aux, temp), leg);
    }
    {
      ScratchScope ta(scratch_arena_a);
      ScratchScope tb(scratch_arena_b);
      MeshOps::classify_faces_by_topology(ambo_p, scratch_arena_a,
                                          scratch_arena_b, leg);
    }
    const size_t nf = ambo_p.face_counts.size();
    std::vector<uint8_t> pal2(nf);
    for (size_t f = 0; f < nf; ++f)
      pal2[f] = static_cast<uint8_t>(
          wrap(static_cast<int>(ambo_p.topology[f]), OpLeg::PALETTES));
    hs_test::StubEffect fx(RW, RH);
    std::vector<Pixel> snaps[3]; // leg-2 last, leg-3 first, leg-3 second
    int drawn = 0, rasterize_at = -1;
    auto cb = [&](Canvas &c, const MeshState &m, const OpLeg::Shading &sh) {
      ++drawn;
      if (drawn != rasterize_at)
        return;
      auto shader = [&](const Vector &, Fragment &frag) {
        int fi = static_cast<int>(frag.v2);
        int ramp =
            (fi >= 0 && fi < static_cast<int>(sh.faces)) ? sh.face_ramp[fi] : 0;
        float t = hs::clamp(fragment_edge_dist(frag) * sh.gain, 0.0f, 1.0f);
        frag.color = sh.ramps[ramp].get(t);
        frag.color.alpha = 255;
      };
      Scan::Mesh::draw<RW, RH>(filters, c, m, shader, scratch_arena_b);
    };
    auto snap = [&](std::vector<Pixel> &out) {
      capture_frame<RW, RH>(fx, out);
    };

    // Leg 2: the medial slerp, departed from ambo(P). Crossfading handoffs,
    // as IslamicStars drives the bridge; the seeded RNG keeps the per-leg
    // target shuffles deterministic here.
    OpLeg::PaletteHandoff handoff2{.bank = &bank.bank,
                                   .prev_face_palette = pal2.data(),
                                   .prev_faces = nf,
                                   .correspondence =
                                       OpLeg::FaceCorrespondence::IDENTITY};
    OpLeg leg2(P, OpLeg::MedialSpec{.sweep_frames = SWEEP}, leg, cb, handoff2);
    const OpLeg::Landing &landing2 = leg2.landing();
    HS_EXPECT_EQ(landing2.faces, nf);
    HS_EXPECT_TRUE(landing2.arrival_topology != nullptr);
    HS_EXPECT_TRUE(landing2.arrival_point != nullptr);
    for (size_t f = 0; f < nf; ++f)
      HS_EXPECT_EQ(landing2.from_palette[f], pal2[f]);
    rasterize_at = SWEEP;
    for (int f = 0; f < SWEEP; ++f) {
      {
        Canvas c(fx);
        leg2.step(c);
      }
      fx.advance_display();
    }
    snap(snaps[0]);

    // Leg-3 handoff, exactly as schedule_dual_untruncate builds it: departed
    // centroids and palettes cached by leg 2.
    std::vector<uint8_t> pal3(nf);
    std::vector<Vector> cen3(nf);
    const PolyMesh &arrival = *landing2.arrival_topology;
    size_t off = 0;
    for (size_t f = 0; f < nf; ++f) {
      Vector c(0.0f, 0.0f, 0.0f);
      for (int j = 0; j < arrival.face_counts[f]; ++j) {
        const size_t v = arrival.faces[off + j];
        c = c + landing2.arrival_point[v].decode().normalized();
      }
      cen3[f] = c.normalized();
      pal3[f] = landing2.landed_palette(f);
      off += arrival.face_counts[f];
    }

    PolyMesh D;
    {
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      D = Solids::finalize_solid(MeshOps::dual(P, aux, temp), leg);
    }
    {
      ScratchScope ta(scratch_arena_a);
      ScratchScope tb(scratch_arena_b);
      MeshOps::classify_faces_by_topology(D, scratch_arena_a, scratch_arena_b,
                                          leg);
    }
    // dual drops sub-3 orbits, so D's faces are exactly the medial's P-vertex
    // block: the two blocks partition the shared face count.
    const size_t DF = D.face_counts.size();
    HS_EXPECT_EQ(DF, nf - PF);

    OpLeg::PaletteHandoff handoff3{.bank = &bank.bank,
                                   .prev_face_palette = pal3.data(),
                                   .prev_faces = nf,
                                   .correspondence =
                                       OpLeg::FaceCorrespondence::DUAL_CLOSING};
    OpLeg::BookendClasses bookend3{.topology = D.topology.data(), .faces = DF};
    OpLeg leg3(D,
               OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                     .t_start = 0.5f,
                                     .t_end = 0.0f,
                                     .sweep_frames = SWEEP,
                                     .bridge_provenance = true,
                                     .borrow_seed = true},
               leg, cb, handoff3, bookend3);
    const OpLeg::Landing &landing3 = leg3.landing();
    HS_EXPECT_EQ(landing3.faces, nf);
    HS_EXPECT_TRUE(landing3.from_palette != nullptr);

    // Derive the true seam permutation from the exact geometry: ambo(D) has
    // the closing leg's face order and the seam's exact positions.
    std::vector<int> perm(nf, -1);
    {
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      PolyMesh ambo_d = MeshOps::ambo(D, aux, temp);
      HS_EXPECT_EQ(ambo_d.face_counts.size(), nf);
      std::vector<char> used(nf, 0);
      int bad_match = 0, non_bijective = 0;
      size_t off = 0;
      for (size_t l = 0; l < nf; ++l) {
        Vector c(0.0f, 0.0f, 0.0f);
        for (int j = 0; j < ambo_d.face_counts[l]; ++j)
          c = c + ambo_d.vertices[ambo_d.faces[off + j]];
        c = c.normalized();
        off += ambo_d.face_counts[l];
        size_t best = 0;
        float bd = 1e9f;
        for (size_t m = 0; m < nf; ++m) {
          const float d = distance_between(c, cen3[m]);
          if (d < bd) {
            bd = d;
            best = m;
          }
        }
        if (bd > SEAM_MATCH_TOL)
          ++bad_match;
        if (used[best])
          ++non_bijective;
        used[best] = 1;
        perm[l] = static_cast<int>(best);
      }
      HS_EXPECT_EQ(bad_match, 0);
      HS_EXPECT_EQ(non_bijective, 0);
    }

    // Block structure: D-faces are the medial's P-vertex orbits in emission
    // order; D-vertex orbit faces land in the medial's P-face block.
    int block1_viol = 0, block2_viol = 0, from_mismatch = 0;
    for (size_t l = 0; l < nf; ++l) {
      if (l < DF && perm[l] != static_cast<int>(PF + l))
        ++block1_viol;
      if (l >= DF && perm[l] >= static_cast<int>(PF))
        ++block2_viol;
      if (landing3.from_palette[l] != pal3[perm[l]])
        ++from_mismatch;
    }
    HS_EXPECT_EQ(block1_viol, 0);
    HS_EXPECT_EQ(block2_viol, 0);
    HS_EXPECT_EQ(from_mismatch, 0);

    // Rendered seam A/B plus a one-step in-leg control.
    drawn = 0;
    rasterize_at = 1;
    {
      Canvas c(fx);
      leg3.step(c);
    }
    fx.advance_display();
    snap(snaps[1]);
    rasterize_at = 2;
    {
      Canvas c(fx);
      leg3.step(c);
    }
    fx.advance_display();
    snap(snaps[2]);

    auto diff = [&](const std::vector<Pixel> &a, const std::vector<Pixel> &b,
                    long long &sumabs, int &changed) {
      sumabs = 0;
      changed = 0;
      for (size_t i = 0; i < a.size(); ++i) {
        const int dr = std::abs((a[i].r >> 8) - (b[i].r >> 8));
        const int dg = std::abs((a[i].g >> 8) - (b[i].g >> 8));
        const int db = std::abs((a[i].b >> 8) - (b[i].b >> 8));
        sumabs += dr + dg + db;
        if (dr | dg | db)
          ++changed;
      }
    };
    long long seam_sum = 0, ctrl_sum = 0;
    int seam_px = 0, ctrl_px = 0;
    diff(snaps[0], snaps[1], seam_sum, seam_px);
    diff(snaps[1], snaps[2], ctrl_sum, ctrl_px);

    const int seam_quiet = RW * RH - seam_px;
    const int ctrl_quiet = RW * RH - ctrl_px;
    const float quiet_ratio =
        static_cast<float>(seam_quiet) / static_cast<float>(ctrl_quiet);
    HS_EXPECT_GT(quiet_ratio, SEAM_QUIET_RATIO);
    HS_EXPECT_LT(seam_sum, SEAM_SUMABS_MAX);

    std::printf("  [opleg seam] %s: F=%zu blocks %zu/%zu, seam diff "
                "sumabs=%lld px=%d (control %lld/%d) quiet=%.3f/%.3f%s\n",
                site.name, nf, DF, nf - DF, seam_sum, seam_px, ctrl_sum,
                ctrl_px, static_cast<double>(quiet_ratio),
                static_cast<double>(SEAM_QUIET_RATIO),
                hs_test::stats().failed != failed_before ? " FAILED" : "");
  }
}

/** @brief Recipe-step leg kind driven by the smoke test. */
enum class StepLegKind { TRUNCATE, SNUB, RELAX };

/** Paused redraws check_step_leg_smoke issues at its pause frame. */
inline constexpr int PAUSED_REDRAWS = 2;

/**
 * @brief Drives one recipe-step leg to completion through OpLeg, gating
 *        compiled-face-count constancy, ramp indices and per-frame motion.
 * @param kind Leg kind under test.
 * @param site Seed site and arrival parameter.
 * @param frames Leg length in frames.
 * @param max_step_chord Per-frame max-vertex-motion bound.
 * @param easing Easing driving the sweep parameter.
 * @param frames_out Optional sink receiving each drawn frame's vertices.
 * @param pause_after Leg frame after which two step_paused() redraws run; 0
 *        drives the leg unpaused.
 * @details Mirrors the hankin smoke test's shape: the departed palettes come
 * from the seed's own classification and the bookend from the step's clean
 * endpoint, exactly as IslamicStars derives them.
 */
inline void
check_step_leg_smoke(StepLegKind kind, const StepLegSite &site, int frames,
                     float max_step_chord, EasingFn easing = ease_in_out_sin,
                     std::vector<std::vector<Vector>> *frames_out = nullptr,
                     int pause_after = 0) {
  using Animation::OpLeg;
  const int failed_before = hs_test::stats().failed;

  reset_globals();
  configure_arenas(GLOBAL_ARENA_SIZE - 24 * 1024 - 32 * 1024, 24 * 1024,
                   32 * 1024);
  hs::random().seed(2026u);

  Arena leg_arena(morph_target_buf, sizeof(morph_target_buf));
  Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  PolyMesh seed = build_step_leg_seed(site, leg_arena);
  PolyMesh endpoint;
  {
    constexpr size_t HALF = sizeof(morph_temp_buf) / 2;
    Arena ea(morph_temp_buf, HALF);
    Arena eb(morph_temp_buf + HALF, HALF);
    PolyMesh raw;
    switch (kind) {
    case StepLegKind::TRUNCATE:
      raw = MeshOps::truncate(seed, ea, eb, site.param);
      break;
    case StepLegKind::SNUB:
      raw = MeshOps::snub(seed, ea, eb, site.param, 0.0f);
      break;
    case StepLegKind::RELAX:
      raw = MeshOps::relax(seed, ea, eb, static_cast<int>(site.param));
      break;
    }
    endpoint = Solids::finalize_solid(raw, leg_arena);
  }

  ScratchScope a_guard(scratch_arena_a);
  ScratchScope b_guard(scratch_arena_b);
  MeshOps::classify_faces_by_topology(seed, scratch_arena_a, scratch_arena_b,
                                      leg_arena);
  MeshOps::classify_faces_by_topology(endpoint, scratch_arena_a,
                                      scratch_arena_b, leg_arena);

  std::array<int, OpLeg::PALETTES> slots;
  MeshPaletteBank::shuffle_indices(slots);

  const size_t prev_faces = seed.face_counts.size();
  uint8_t *prev_pal = leg_arena.allocate_n<uint8_t>(prev_faces);
  Vector *prev_centroid = leg_arena.allocate_n<Vector>(prev_faces);
  size_t off = 0;
  for (size_t f = 0; f < prev_faces; ++f) {
    prev_pal[f] = static_cast<uint8_t>(
        slots[wrap(static_cast<int>(seed.topology[f]), OpLeg::PALETTES)]);
    const int n = seed.face_counts[f];
    prev_centroid[f] = face_centroid_unit(seed, off, n);
    off += n;
  }

  OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                .prev_face_palette = prev_pal,
                                .prev_faces = prev_faces,
                                .prev_face_centroid = prev_centroid};
  OpLeg::BookendClasses bookend{.topology = endpoint.topology.data(),
                                .faces = endpoint.face_counts.size()};

  LegDrawProbe probe;
  auto cb = [&](Canvas &, const MeshState &m, const OpLeg::Shading &sh) {
    probe.observe(m, sh);
    if (frames_out)
      frames_out->push_back(probe.prev_v);
  };

  hs_test::StubEffect fx(288, 144);

  auto run = [&](OpLeg &&leg) {
    const OpLeg::Landing &landing = leg.landing();
    HS_EXPECT_EQ(landing.primary_faces, seed.face_counts.size());
    HS_EXPECT_TRUE(landing.topology != nullptr);
    for (int f = 0; f < frames; ++f) {
      {
        Canvas c(fx);
        leg.step(c);
      }
      fx.advance_display();
      if (f + 1 != pause_after)
        continue;
      for (int p = 0; p < PAUSED_REDRAWS; ++p) {
        {
          Canvas c(fx);
          leg.step_paused(c);
        }
        fx.advance_display();
      }
    }
    HS_EXPECT_EQ(probe.drawn,
                 (size_t)(frames + (pause_after > 0 ? PAUSED_REDRAWS : 0)));
    HS_EXPECT_EQ(probe.faces, landing.faces);
  };

  switch (kind) {
  case StepLegKind::TRUNCATE:
    run(OpLeg(seed,
              OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                    .t_start = 0.0f,
                                    .t_end = site.param,
                                    .sweep_frames = frames},
              leg_arena, cb, handoff, bookend, OpLeg::classic_blend, easing));
    break;
  case StepLegKind::SNUB:
    run(OpLeg(seed,
              OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::SNUB,
                                    .t_start = 0.0f,
                                    .t_end = site.param,
                                    .sweep_frames = frames},
              leg_arena, cb, handoff, bookend, OpLeg::classic_blend, easing));
    break;
  case StepLegKind::RELAX:
    run(OpLeg(seed,
              OpLeg::RelaxSpec{.iterations = static_cast<int>(site.param),
                               .sweep_frames = frames},
              leg_arena, cb, handoff, bookend, OpLeg::classic_blend, easing));
    break;
  }

  HS_EXPECT_LT(probe.worst_step, max_step_chord);
  const char *label = kind == StepLegKind::TRUNCATE ? "truncate"
                      : kind == StepLegKind::SNUB   ? "snub"
                                                    : "relax";
  std::printf("  [opleg %s] %s: %zu faces, worst per-frame vertex step %.4f "
              "chord (bound %.2f)%s\n",
              label, site.name, probe.faces, (double)probe.worst_step,
              (double)max_step_chord,
              hs_test::stats().failed != failed_before ? " FAILED" : "");
}

/** Gated-swap sites: the seeds the partition recipes run kis and dual on. */
inline constexpr StepLegSite GATE_LEG_SITES[] = {
    {"icosahedron", probe_icosahedron, 0.0f},
    {"icosahedron_snub", probe_icosa_snub, 0.0f},
};

/**
 * @brief Drives one gated-swap leg to completion through OpLeg, gating the
 *        per-side compiled-face-count constancy, the single swap, the ramp
 *        indices and the gain envelope.
 * @param op Partition operator the leg swaps to.
 * @param site Seed site.
 * @param gate Half-gate length in frames.
 * @details Derives the handoff and bookend exactly as IslamicStars does: the
 * departed palettes from the seed's own classification, the bookend from the
 * op's clean endpoint.
 */
inline void check_gated_leg_smoke(Animation::OpLeg::SwapOp op,
                                  const StepLegSite &site, int gate) {
  using Animation::OpLeg;
  const int failed_before = hs_test::stats().failed;
  const bool is_kis = op == OpLeg::SwapOp::KIS;

  reset_globals();
  configure_arenas(GLOBAL_ARENA_SIZE - 24 * 1024 - 32 * 1024, 24 * 1024,
                   32 * 1024);
  hs::random().seed(2026u);

  Arena leg_arena(morph_target_buf, sizeof(morph_target_buf));
  Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  PolyMesh seed = build_step_leg_seed(site, leg_arena);
  PolyMesh endpoint;
  {
    constexpr size_t HALF = sizeof(morph_temp_buf) / 2;
    Arena ea(morph_temp_buf, HALF);
    Arena eb(morph_temp_buf + HALF, HALF);
    endpoint = Solids::finalize_solid(is_kis ? MeshOps::kis(seed, ea, eb)
                                             : MeshOps::dual(seed, ea, eb),
                                      leg_arena);
  }

  ScratchScope a_guard(scratch_arena_a);
  ScratchScope b_guard(scratch_arena_b);
  MeshOps::classify_faces_by_topology(seed, scratch_arena_a, scratch_arena_b,
                                      leg_arena);
  MeshOps::classify_faces_by_topology(endpoint, scratch_arena_a,
                                      scratch_arena_b, leg_arena);

  std::array<int, OpLeg::PALETTES> slots;
  MeshPaletteBank::shuffle_indices(slots);

  const size_t prev_faces = seed.face_counts.size();
  uint8_t *prev_pal = leg_arena.allocate_n<uint8_t>(prev_faces);
  for (size_t f = 0; f < prev_faces; ++f)
    prev_pal[f] = static_cast<uint8_t>(
        slots[wrap(static_cast<int>(seed.topology[f]), OpLeg::PALETTES)]);

  OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                .prev_face_palette = prev_pal,
                                .prev_faces = prev_faces};
  OpLeg::BookendClasses bookend{.topology = endpoint.topology.data(),
                                .faces = endpoint.face_counts.size()};

  const int frames = 2 * gate + 1;
  int drawn = 0, swaps = 0, swap_frame = -1;
  size_t side_faces = 0;
  auto cb = [&](Canvas &, const MeshState &m, const OpLeg::Shading &sh) {
    HS_EXPECT_EQ(m.face_counts.size(), sh.faces);
    for (size_t f = 0; f < sh.faces; ++f)
      HS_EXPECT_LT(static_cast<int>(sh.face_ramp[f]), OpLeg::MAX_BLEND_PAIRS);

    // De-flash: the gate never dims. Gain is exactly 1 on every frame.
    HS_EXPECT_EQ(sh.gain, 1.0f);

    if (drawn == 0) {
      side_faces = sh.faces;
    } else if (sh.faces != side_faces) {
      ++swaps;
      swap_frame = drawn;
      side_faces = sh.faces;
    }
    ++drawn;
  };

  hs_test::StubEffect fx(288, 144);

  OpLeg leg(seed, OpLeg::GatedSwapSpec{.op = op, .gate_frames = gate},
            leg_arena, cb, handoff, bookend);
  const OpLeg::Landing &landing = leg.landing();
  for (int f = 0; f < frames; ++f) {
    {
      Canvas c(fx);
      leg.step(c);
    }
    fx.advance_display();
  }

  HS_EXPECT_EQ(drawn, frames);
  // Topology changes exactly once, at the swap frame, and the leg lands on the
  // arrival face count.
  HS_EXPECT_EQ(swaps, 1);
  HS_EXPECT_EQ(swap_frame, gate);
  HS_EXPECT_EQ(side_faces, landing.faces);
  HS_EXPECT_EQ(landing.faces, endpoint.face_counts.size());
  HS_EXPECT_EQ(landing.primary_faces, prev_faces);
  HS_EXPECT_TRUE(landing.from_palette != nullptr);

  std::printf("  [opleg %s] %s: F %zu -> %zu, swap at frame %d of %d, gain 1 "
              "throughout%s\n",
              is_kis ? "kis" : "dual", site.name, prev_faces, landing.faces,
              swap_frame, frames,
              hs_test::stats().failed != failed_before ? " FAILED" : "");
}

/**
 * @brief Smoke-tests a gated-swap leg for both partition ops on every gate
 *        site, at IslamicStars' 6-frame half-gate.
 */
inline void test_opleg_gated_swap_smoke() {
  using Animation::OpLeg;
  for (const StepLegSite &site : GATE_LEG_SITES) {
    check_gated_leg_smoke(OpLeg::SwapOp::KIS, site, 6);
    check_gated_leg_smoke(OpLeg::SwapOp::DUAL, site, 6);
  }
}

/**
 * @brief Smoke-tests one leg of each new recipe-step kind end to end.
 */
inline void test_opleg_step_leg_smoke() {
  // Leg lengths mirror IslamicStars' budget (spec section 7): 24 frames for an
  // operator sweep, 16 for a standalone relax.
  check_step_leg_smoke(StepLegKind::TRUNCATE, TRUNCATE_LEG_SITES[0], 24, 0.15f);
  check_step_leg_smoke(StepLegKind::SNUB, SNUB_LEG_SITES[0], 24, 0.15f);
  for (const StepLegSite &site : RELAX_LEG_SITES)
    check_step_leg_smoke(StepLegKind::RELAX, site, 16, 0.15f);
}

/**
 * @brief Pins a paused leg's redraw: the held frame, drawn again, with the
 *        sweep clock untouched.
 * @details A paused timeline drives every event through step_paused(); a leg
 * that dropped the redraw would blank the frame and one that advanced would run
 * the sweep out under the pause. Both paused draws must reproduce the held
 * frame's vertices bitwise, and the resumed frames must match the unpaused run
 * frame for frame.
 */
inline void test_opleg_step_paused_holds_frame() {
  constexpr int FRAMES = 8;
  constexpr int PAUSE_AFTER = 4;
  // Short leg: the chord bound is the smoke test's concern, not this one.
  constexpr float CHORD_MAX = 1.0f;
  std::vector<std::vector<Vector>> unpaused, held;
  check_step_leg_smoke(StepLegKind::TRUNCATE, TRUNCATE_LEG_SITES[0], FRAMES,
                       CHORD_MAX, ease_in_out_sin, &unpaused);
  check_step_leg_smoke(StepLegKind::TRUNCATE, TRUNCATE_LEG_SITES[0], FRAMES,
                       CHORD_MAX, ease_in_out_sin, &held, PAUSE_AFTER);
  HS_EXPECT_EQ(unpaused.size(), (size_t)FRAMES);
  HS_EXPECT_EQ(held.size(), (size_t)(FRAMES + PAUSED_REDRAWS));
  if (held.size() != (size_t)(FRAMES + PAUSED_REDRAWS))
    return;

  auto identical = [](const std::vector<Vector> &a,
                      const std::vector<Vector> &b) {
    if (a.size() != b.size() || a.empty())
      return false;
    for (size_t i = 0; i < a.size(); ++i)
      if (a[i].x != b[i].x || a[i].y != b[i].y || a[i].z != b[i].z)
        return false;
    return true;
  };

  // Every paused draw repeats the frame the leg was holding.
  for (int p = 0; p < PAUSED_REDRAWS; ++p)
    HS_EXPECT_TRUE(identical(held[PAUSE_AFTER - 1], held[PAUSE_AFTER + p]));

  // The clock never moved: the resumed frames are the unpaused run's.
  size_t drifted = 0;
  for (int f = 0; f < FRAMES; ++f) {
    const size_t k = f < PAUSE_AFTER ? (size_t)f : (size_t)(f + PAUSED_REDRAWS);
    if (!identical(unpaused[f], held[k]))
      ++drifted;
  }
  HS_EXPECT_EQ(drifted, (size_t)0);
  std::printf("  [opleg paused] truncate: %d frames, %d paused redraws at "
              "frame %d, %zu resumed frames drifted\n",
              FRAMES, PAUSED_REDRAWS, PAUSE_AFTER, drifted);
}

/**
 * @brief Drives swept legs under an easing whose range leaves [0, 1], gating
 *        the sweep parameter against extrapolation past the arrival.
 * @details ease_out_elastic exceeds 1 over x in [0.075, 0.225], peaking near
 * 1.35. Unclamped, that drives the operator past the endpoint the constructor
 * pinned inside its topology-constant interval, and past the domain the ops
 * themselves guard (a reverse leg reaches negative t). Clamped, frames 2-5 of
 * 24 hold the arrival exactly, which the frame-4-vs-last comparison pins; frame
 * 1 is still mid-sweep, so a leg frozen from the start would not pass. The
 * chord bound is loose because elastic covers half the sweep in one frame by
 * design.
 */
inline void test_opleg_step_leg_overshooting_easing() {
  constexpr StepLegSite NEAR_AMBO{"icosahedron_ambo_truncate049",
                                  probe_icosa_ambo, 0.49f};
  constexpr int FRAMES = 24;
  for (int k = 0; k < 2; ++k) {
    std::vector<std::vector<Vector>> drawn;
    const bool truncate = k == 0;
    check_step_leg_smoke(truncate ? StepLegKind::TRUNCATE : StepLegKind::SNUB,
                         truncate ? NEAR_AMBO : SNUB_LEG_SITES[0], FRAMES, 2.0f,
                         ease_out_elastic, &drawn);
    HS_EXPECT_EQ(drawn.size(), (size_t)FRAMES);
    const std::vector<Vector> &arrival = drawn.back();
    const std::vector<Vector> &peak = drawn[3];
    const std::vector<Vector> &opening = drawn[0];
    HS_EXPECT_EQ(peak.size(), arrival.size());
    float worst_peak = 0.0f, worst_opening = 0.0f;
    for (size_t i = 0; i < arrival.size(); ++i) {
      worst_peak = std::max(worst_peak, distance_between(peak[i], arrival[i]));
      worst_opening =
          std::max(worst_opening, distance_between(opening[i], arrival[i]));
    }
    HS_EXPECT_LT(worst_peak, 1e-5f);
    HS_EXPECT_GT(worst_opening, 1e-3f);
    std::printf("  [opleg elastic] %s: overshoot frame %.2e from arrival, "
                "opening frame %.4f\n",
                truncate ? "truncate" : "snub", (double)worst_peak,
                (double)worst_opening);
  }
}

/**
 * @brief Pins the crossfade colour model on a ConwayGraph edge leg under the
 *        classic_blend default: exact `from` at frame 1, moved off `from` by
 *        mid-leg, exact `to` at the arrival frame.
 */
inline void test_opleg_edge_leg_crossfade() {
  using Animation::OpLeg;
  reset_globals();
  configure_arenas(GLOBAL_ARENA_SIZE - 116 * 1024 - 74 * 1024, 116 * 1024,
                   74 * 1024);
  hs::random().seed(2026u);

  Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  hs_test::StubEffect fx(288, 144);

  // LUT grid-aligned sample coordinates for exact ramp-color comparisons.
  constexpr float SAMPLES[] = {0.0f, 0.5f, 1.0f};
  const OpLeg::Landing *lp = nullptr;
  std::vector<char> all_from, all_to; // per drawn frame: exact from/to palettes
  auto matches = [&](const OpLeg::Shading &sh, size_t f, uint8_t pal) {
    for (float t : SAMPLES) {
      const Color4 got = sh.ramps[sh.face_ramp[f]].get(t);
      const Color4 exp = bank.bank.entries[pal].get(t);
      if (got.color.r != exp.color.r || got.color.g != exp.color.g ||
          got.color.b != exp.color.b)
        return false;
    }
    return true;
  };
  auto cb = [&](Canvas &, const MeshState &, const OpLeg::Shading &sh) {
    bool from_ok = true, to_ok = true;
    for (size_t f = 0; f < sh.faces; ++f) {
      const uint8_t from = lp->from_palette[f];
      const uint8_t to = lp->to_palette[wrap(static_cast<int>(lp->topology[f]),
                                             OpLeg::PALETTES)];
      from_ok = from_ok && matches(sh, f, from);
      to_ok = to_ok && matches(sh, f, to);
    }
    all_from.push_back(from_ok);
    all_to.push_back(to_ok);
  };
  auto divergent_faces = [&]() {
    int n = 0;
    for (size_t f = 0; f < lp->faces; ++f)
      if (lp->from_palette[f] !=
          lp->to_palette[wrap(static_cast<int>(lp->topology[f]),
                              OpLeg::PALETTES)])
        ++n;
    return n;
  };
  auto run_frames = [&](OpLeg &leg, int frames) {
    all_from.clear();
    all_to.clear();
    lp = &leg.landing();
    for (int f = 0; f < frames; ++f) {
      {
        Canvas c(fx);
        leg.step(c);
      }
      fx.advance_display();
    }
    HS_EXPECT_EQ(all_from.size(), (size_t)frames);
  };
  {
    const int failed_before = hs_test::stats().failed;
    Arena leg_arena(morph_target_buf, sizeof(morph_target_buf));
    PolyMesh cube;
    build_solid<Solids::Cube>(cube, leg_arena);
    uint8_t pal[16];
    for (size_t f = 0; f < cube.face_counts.size(); ++f)
      pal[f] = static_cast<uint8_t>(f % OpLeg::PALETTES);
    const int edge = [] {
      for (int e = 0; e < ConwayGraph::NUM_EDGES; ++e)
        if (ConwayGraph::EDGES[e].from_node == ConwayGraph::CUBE &&
            ConwayGraph::EDGES[e].to_node == ConwayGraph::TRUNCATED_CUBE)
          return e;
      return -1;
    }();
    HS_EXPECT_GE(edge, 0);
    OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                  .prev_face_palette = pal,
                                  .prev_faces = cube.face_counts.size()};
    constexpr int EDGE_FRAMES = 24;
    OpLeg leg(cube,
              OpLeg::EdgeSweepSpec{.edge = &ConwayGraph::EDGES[edge],
                                   .sweep_frames = EDGE_FRAMES},
              leg_arena, cb, handoff);
    run_frames(leg, EDGE_FRAMES);
    HS_EXPECT_GT(divergent_faces(), 0);
    HS_EXPECT_TRUE(all_from[0]);
    HS_EXPECT_TRUE(!all_from[EDGE_FRAMES / 2 - 1]);
    HS_EXPECT_TRUE(all_to[EDGE_FRAMES - 1]);
    std::printf("  [crossfade] edge leg: F=%zu mid-leg blend intact%s\n",
                lp->faces,
                hs_test::stats().failed != failed_before ? " FAILED" : "");
  }
}

/**
 * @brief Gates the recipe steps no leg kind covers and reports the clamp
 *        consequence the gate exists to avoid.
 */
inline void test_unsweepable_recipe_steps_are_gated() {
  using Solids::Op;
  using Solids::IslamicStarPatterns::D2R;
  using Solids::IslamicStarPatterns::TRUNCATE_T_FAR;

  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::TRUNCATE, 0.33f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::SNUB, 0.5f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::RELAX, 8.0f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::HANKIN, 62.0f * D2R}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::AMBO}));
  // 0.01 is below T_EPS but sweepable: the leg births at the derived
  // per-arrival floor instead of clamping to a still image
  // (opchain_morph_spec 5.1).
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::TRUNCATE, 0.01f}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::TRUNCATE, 0.001f}));
  // 0.873 is a far-side arrival past the ambo pinch: now a real sweeping leg
  // (opchain_morph_spec 5.1), not a clamped clean-swap. A fully-crossed t == 1
  // stays blocked (cut faces collapse).
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::TRUNCATE, TRUNCATE_T_FAR}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::TRUNCATE, 1.0f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::CHAMFER, 0.63f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::KIS}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::DUAL}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::CHAMFER, 0.001f}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::CHAMFER, 0.9f}));
}

// ---------------------------------------------------------------------------
// Recipe chain build replay (docs/opchain_morph_spec.md sections 9.2/9.3):
// every registry entry with a non-null recipe is lowered and replayed leg by
// leg exactly as IslamicStars builds it — same crossfading handoff and
// bookend — stepping every leg frame by frame with a recording draw callback.
// Gates per-leg compiled-face-count constancy, the crossfade colour model
// (every frame of every leg draws every face's (from, to) ramp bit-exact at
// the leg's blend weight; each leg departs from the palette the previous leg
// landed on), the final per-face sprite handoff, the per-final-class palette
// symmetry of the finished shape, the distinct blend-pair ceiling, and the
// persistent/scratch high-water against IslamicStars' configured split.
// ---------------------------------------------------------------------------

/** The arena split is canvas-independent; this instantiation names it. */
using IslamicFx = IslamicStars<288, 144>;

constexpr size_t ISLAMIC_SCRATCH_A_BUDGET =
    IslamicFx::SPLIT_SCRATCH_A_DEFAULT; /**< IslamicStars scratch_a split. */
constexpr size_t ISLAMIC_SCRATCH_B_BUDGET =
    IslamicFx::SPLIT_SCRATCH_B_BUILD; /**< IslamicStars build scratch_b. */
/** Device persistent budget of IslamicStars' arena split. */
constexpr size_t ISLAMIC_PERSISTENT_BUDGET = DEVICE_GLOBAL_ARENA_SIZE -
                                             ISLAMIC_SCRATCH_A_BUDGET -
                                             ISLAMIC_SCRATCH_B_BUDGET;
/** Scratch capacity the replay runs with, above every budget it gates. */
constexpr size_t REPLAY_SCRATCH_CAPACITY = 512 * 1024;

/**
 * @brief Arena high-waters and shape of one replayed chain.
 */
struct ChainPeaks {
  size_t persistent = 0;  /**< persistent_arena high-water, bytes. */
  size_t scratch_a = 0;   /**< scratch_arena_a high-water, bytes. */
  size_t scratch_b = 0;   /**< scratch_arena_b high-water, bytes. */
  size_t legs = 0;        /**< Lowered primitive step count. */
  size_t faces = 0;       /**< Face count of the finished solid. */
  int palettes = 0;       /**< Distinct palettes on the finished shape. */
  int final_classes = 0;  /**< Distinct newborn classes on the final leg. */
  int blend_pairs = 0;    /**< Max distinct (from, to) pairs on any leg. */
  bool supported = false; /**< Every lowered step has a leg kind. */
};

/**
 * @brief Replays one recipe leg by leg as IslamicStars builds it and gates
 *        continuity, classification, and the arena high-waters.
 * @param name Diagnostic label.
 * @param recipe Recipe replayed.
 * @param gate Whether to assert the budgets and print the per-chain line; false
 * measures a chain the effect would reject, for the arena survey.
 * @return The chain's arena high-waters.
 */
inline ChainPeaks replay_build_chain(const char *name,
                                     const Solids::Recipe &recipe,
                                     bool gate = true) {
  ChainPeaks peaks;
  using Animation::OpLeg;
  constexpr int HANKIN_LEG_FRAMES = 32, SWEEP_LEG_FRAMES = 24,
                RELAX_LEG_FRAMES = 16, GATE_HALF_FRAMES = 6;
  constexpr size_t MAX_FACES = 1152;
  constexpr size_t MAX_STEPS = 8;

  {
    const int failed_before = hs_test::stats().failed;

    reset_globals();
    // Capacities are the host's, not the device split: a chain that overruns a
    // budget must report its high-water, not OOM-trap the replay before the
    // measurement. The budgets below are what the peaks are gated against.
    configure_arenas(GLOBAL_ARENA_SIZE - 2 * REPLAY_SCRATCH_CAPACITY,
                     REPLAY_SCRATCH_CAPACITY, REPLAY_SCRATCH_CAPACITY);
    hs::random().seed(2026u);

    MeshPaletteBank bank;
    bank.bake_all(persistent_arena);

    // Seed solid: persistent PolyMesh plus the compiled/classified carousel
    // slot, exactly as spawn_shape prepares them.
    PolyMesh cur;
    MeshState seed_slot;
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      cur = Solids::finalize_solid(
          Solids::simple_registry[recipe.seed].generate(a, b), target);
      MeshOps::compile(cur, seed_slot, target, a);
    });
    {
      ScratchScope a_guard(scratch_arena_a);
      ScratchScope b_guard(scratch_arena_b);
      MeshOps::classify_faces_by_topology(seed_slot, scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }

    // The spawned seed's shuffled palette order, consumed by class ordinal.
    std::array<int, OpLeg::PALETTES> slots;
    MeshPaletteBank::shuffle_indices(slots);
    std::array<uint8_t, OpLeg::PALETTES> order;
    for (int i = 0; i < OpLeg::PALETTES; ++i)
      order[i] = static_cast<uint8_t>(slots[i]);

    Solids::OpStep steps[MAX_STEPS];
    const size_t count = Solids::expand_to_primitives(recipe, steps, MAX_STEPS);
    HS_EXPECT_GT(count, (size_t)0);
    bool supported = true;
    int leg_frames[MAX_STEPS] = {};
    bool gated[MAX_STEPS] = {};
    for (size_t k = 0; k < count; ++k) {
      if (!Solids::is_morphable_step(steps[k]))
        supported = false;
      gated[k] =
          steps[k].op == Solids::Op::KIS || steps[k].op == Solids::Op::DUAL;
      leg_frames[k] = gated[k] ? 2 * GATE_HALF_FRAMES + 1
                      : steps[k].op == Solids::Op::HANKIN ? HANKIN_LEG_FRAMES
                      : steps[k].op == Solids::Op::RELAX  ? RELAX_LEG_FRAMES
                                                          : SWEEP_LEG_FRAMES;
    }
    peaks.legs = count;
    peaks.supported = supported;
    if (gate) {
      HS_EXPECT_TRUE(supported);
      if (!supported) {
        std::printf("    [chain] %s: unsupported lowered op\n", name);
        return peaks;
      }
    }

    hs_test::StubEffect fx(288, 144);
    uint8_t prev_pal_buf[MAX_FACES];
    Vector prev_centroid[MAX_FACES];
    uint8_t carried_to[MAX_FACES] = {};
    const bool pin_final = supported;
    std::vector<int> full_topo;
    const OpLeg::Landing *prev_landing = nullptr;
    PolyMesh next;
    // A hankin leg's endpoint is rebuilt into scratch_a and stays there until
    // the boundary evacuates it into the fresh persistent arena, as
    // finish_build_leg does. One is live at a time, so each rebuild rewinds to
    // this mark first.
    const size_t endpoint_mark = scratch_arena_a.get_offset();

    for (size_t k = 0; k < count; ++k) {
      const size_t prev_faces = cur.face_counts.size();
      HS_EXPECT_LE(prev_faces, MAX_FACES);
      // The replay threads state leg to leg, so an over-cap leg cannot be
      // skipped: abandon the chain instead of writing past the handoff arrays.
      if (prev_faces > MAX_FACES)
        return peaks;
      {
        size_t off = 0;
        for (size_t f = 0; f < prev_faces; ++f) {
          const int n = cur.face_counts[f];
          prev_centroid[f] = face_centroid_unit(cur, off, n);
          off += n;
        }
      }
      const uint8_t *prev_pal;
      if (k == 0) {
        // Seed keying, as spawn_shape does it: class ids are dense, so class
        // c is the c-th cohort.
        for (size_t f = 0; f < prev_faces; ++f) {
          const int cls = seed_slot.topology[f];
          prev_pal_buf[f] =
              static_cast<uint8_t>(order[wrap(cls, OpLeg::PALETTES)]);
        }
        prev_pal = prev_pal_buf;
      } else {
        for (size_t f = 0; f < prev_faces; ++f)
          prev_pal_buf[f] = carried_to[f];
        prev_pal = prev_pal_buf;
      }

      // Eager clean endpoint, exactly as start_build_leg derives it (for the
      // ambo leg: the AMBO mesh, never the swept truncate form). A hankin step
      // builds none: its leg rebuilds its own arrival once it has run.
      const bool hankin_step = steps[k].op == Solids::Op::HANKIN;
      OpLeg::BookendClasses bookend;
      if (!hankin_step) {
        hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
          PolyMesh nx;
          switch (steps[k].op) {
          case Solids::Op::AMBO:
            nx = MeshOps::ambo(cur, a, b);
            break;
          case Solids::Op::TRUNCATE:
            nx = MeshOps::truncate(cur, a, b, steps[k].param);
            break;
          case Solids::Op::SNUB:
            nx = MeshOps::snub(cur, a, b, steps[k].param, steps[k].twist);
            break;
          case Solids::Op::KIS:
            nx = MeshOps::kis(cur, a, b);
            break;
          case Solids::Op::DUAL:
            nx = MeshOps::dual(cur, a, b);
            break;
          case Solids::Op::CHAMFER:
            nx = MeshOps::chamfer(cur, a, b, steps[k].param);
            break;
          default:
            nx = steps[k].bake
                     ? MeshOps::relax_baked(cur, a, *steps[k].bake)
                     : MeshOps::relax(cur, a, b,
                                      static_cast<int>(steps[k].param));
            break;
          }
          next = Solids::finalize_solid(nx, target);
        });
        {
          ScratchScope a_guard(scratch_arena_a);
          ScratchScope b_guard(scratch_arena_b);
          MeshOps::classify_faces_by_topology(
              next, scratch_arena_a, scratch_arena_b, persistent_arena);
        }
        const size_t bookend_faces = next.face_counts.size();
        HS_EXPECT_LE(bookend_faces, MAX_FACES);
        bookend = {.topology = next.topology.data(), .faces = bookend_faces};
      }

      OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                    .prev_face_palette = prev_pal,
                                    .prev_faces = prev_faces,
                                    .prev_face_centroid = prev_centroid};

      size_t drawn = 0;
      size_t leg_faces = 0;
      int off_palette_frames = 0;
      const OpLeg::Landing *lp = nullptr;
      // LUT grid-aligned sample coordinates for exact ramp-color comparisons.
      constexpr float PROBE_T[] = {0.0f, 0.5f, 1.0f};
      Arena blend_arena(morph_temp_buf, sizeof(morph_temp_buf));
      // A gated leg's face count is constant per side and changes once, at the
      // swap; every other kind holds one count for the whole leg. Every frame
      // must draw every face's (from, to) ramp bit-exact at the frame's blend
      // weight: endpoint weights alias the bank LUTs, interior weights are
      // rebuilt here through the same bake_palette_blend the leg uses.
      auto cb = [&](Canvas &, const MeshState &m, const OpLeg::Shading &sh) {
        HS_EXPECT_EQ(m.face_counts.size(), sh.faces);
        const bool side_start =
            drawn == 0 ||
            (gated[k] && drawn == static_cast<size_t>(GATE_HALF_FRAMES));
        if (side_start)
          leg_faces = sh.faces;
        else
          HS_EXPECT_EQ(sh.faces, leg_faces);
        const bool seed_side =
            gated[k] && drawn < static_cast<size_t>(GATE_HALF_FRAMES);
        const int frame = static_cast<int>(drawn) + 1;
        const float w = gated[k] ? OpLeg::trailing_blend(frame, leg_frames[k])
                                 : OpLeg::classic_blend(frame, leg_frames[k]);
        blend_arena.reset();
        BakedPalette expected[OpLeg::PALETTES][OpLeg::PALETTES];
        bool baked[OpLeg::PALETTES][OpLeg::PALETTES] = {};
        bool frame_ok = true;
        for (size_t f = 0; f < sh.faces && lp; ++f) {
          const uint8_t from =
              seed_side ? prev_pal_buf[f] : lp->from_palette[f];
          const uint8_t to = seed_side ? from : lp->landed_palette(f);
          const BakedPalette *want = &bank.bank.entries[to];
          if (from != to) {
            if (!baked[from][to]) {
              bake_palette_blend(expected[from][to], blend_arena,
                                 bank.bank.entries[from], bank.bank.entries[to],
                                 w);
              baked[from][to] = true;
            }
            want = &expected[from][to];
          }
          for (float t : PROBE_T) {
            const Color4 got = sh.ramps[sh.face_ramp[f]].get(t);
            const Color4 exp = want->get(t);
            if (got.color.r != exp.color.r || got.color.g != exp.color.g ||
                got.color.b != exp.color.b)
              frame_ok = false;
          }
        }
        if (!frame_ok)
          ++off_palette_frames;
        ++drawn;
      };

      auto make_leg = [&]() {
        switch (steps[k].op) {
        case Solids::Op::KIS:
          return OpLeg(cur,
                       OpLeg::GatedSwapSpec{.op = OpLeg::SwapOp::KIS,
                                            .gate_frames = GATE_HALF_FRAMES},
                       persistent_arena, cb, handoff, bookend);
        case Solids::Op::DUAL:
          return OpLeg(cur,
                       OpLeg::GatedSwapSpec{.op = OpLeg::SwapOp::DUAL,
                                            .gate_frames = GATE_HALF_FRAMES},
                       persistent_arena, cb, handoff, bookend);
        case Solids::Op::HANKIN:
          return OpLeg(cur,
                       OpLeg::HankinSweepSpec{.theta_start = 0.0f,
                                              .theta_end = steps[k].param,
                                              .sweep_frames = leg_frames[k]},
                       persistent_arena, cb, handoff, bookend);
        case Solids::Op::AMBO:
          return OpLeg(
              cur,
              OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                    .t_start = 0.0f,
                                    .t_end = 0.5f,
                                    .sweep_frames = leg_frames[k]},
              persistent_arena, cb, handoff, bookend);
        case Solids::Op::TRUNCATE:
          return OpLeg(
              cur,
              OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                    .t_start = 0.0f,
                                    .t_end = steps[k].param,
                                    .sweep_frames = leg_frames[k]},
              persistent_arena, cb, handoff, bookend);
        case Solids::Op::SNUB:
          return OpLeg(cur,
                       OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::SNUB,
                                             .t_start = 0.0f,
                                             .t_end = steps[k].param,
                                             .twist_end = steps[k].twist,
                                             .sweep_frames = leg_frames[k]},
                       persistent_arena, cb, handoff, bookend);
        case Solids::Op::CHAMFER:
          return OpLeg(
              cur,
              OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::CHAMFER,
                                    .t_start = 0.0f,
                                    .t_end = steps[k].param,
                                    .sweep_frames = leg_frames[k]},
              persistent_arena, cb, handoff, bookend);
        default:
          return OpLeg(
              cur,
              OpLeg::RelaxSpec{.iterations = static_cast<int>(steps[k].param),
                               .bake = steps[k].bake,
                               .sweep_frames = leg_frames[k]},
              persistent_arena, cb, handoff, bookend);
        }
      };
      OpLeg leg = make_leg();
      const OpLeg::Landing &landing = leg.landing();
      lp = &landing;

      for (int f = 0; f < leg_frames[k]; ++f) {
        {
          Canvas c(fx);
          leg.step(c);
        }
        fx.advance_display();
      }
      HS_EXPECT_EQ(drawn, (size_t)leg_frames[k]);
      HS_EXPECT_EQ(landing.faces, leg_faces);
      HS_EXPECT_TRUE(landing.from_palette != nullptr);
      // Every frame drew every face's (from, to) ramp bit-exact.
      HS_EXPECT_EQ(off_palette_frames, 0);

      // The landed carry: a surviving face departs the next leg from the
      // palette this leg landed it on. A partition op keeps no emission-order
      // correspondence to its seed, so its from-palettes are the leg's own
      // geometric provenance.
      HS_EXPECT_LE(landing.faces, MAX_FACES);
      if (landing.faces > MAX_FACES)
        return peaks;
      if (!supported) {
        // A clamped leg lands somewhere other than its clean endpoint, so
        // neither the palette correspondence nor the closing handoff below is
        // meaningful for it.
      } else if (k > 0 && !gated[k]) {
        // First-offending face rather than a per-face assertion: the replay
        // sweeps tens of thousands of faces, which would swamp the module's
        // assertion floor.
        int first_broken_carry = -1;
        for (size_t f = 0; f < prev_faces; ++f)
          if ((int)landing.from_palette[f] != (int)carried_to[f] &&
              first_broken_carry < 0)
            first_broken_carry = (int)f;
        HS_EXPECT_EQ(first_broken_carry, -1);
      }
      if (supported) {
        // The leg's fresh target set is a permutation of the bank.
        bool seen_to[OpLeg::PALETTES] = {};
        for (int i = 0; i < OpLeg::PALETTES; ++i) {
          HS_EXPECT_LE((int)landing.to_palette[i], OpLeg::PALETTES - 1);
          seen_to[landing.to_palette[i]] = true;
        }
        for (bool s : seen_to)
          HS_EXPECT_TRUE(s);
        // Final-leg birth grain: distinct newborn classes (the perceptual
        // grouping report).
        if (k + 1 == count) {
          const size_t births_from = gated[k] ? 0 : prev_faces;
          int classes = 0;
          int prev_cls = -1;
          for (;;) {
            int cls = -1;
            for (size_t f = births_from; f < landing.faces; ++f) {
              const int c = landing.topology[f];
              if (c > prev_cls && (cls < 0 || c < cls))
                cls = c;
            }
            if (cls < 0)
              break;
            ++classes;
            prev_cls = cls;
          }
          peaks.final_classes = classes;
        }
      }
      peaks.blend_pairs = std::max(peaks.blend_pairs, landing.blend_pairs);
      HS_EXPECT_LE(landing.blend_pairs, OpLeg::MAX_BLEND_PAIRS);
      for (size_t f = 0; f < landing.faces; ++f)
        carried_to[f] = landing.landed_palette(f);

      // Landing lives in the leg's arena-backed Transients; outlives `leg`.
      prev_landing = &landing;

      if (hankin_step) {
        next = PolyMesh();
        scratch_arena_a.set_offset(endpoint_mark);
        OpLeg::arrival_mesh(landing, next, scratch_arena_a);
      }

      // Full-precision final classification for the symmetry pin, rebuilt as
      // the leg's constructor builds its arrival (before the snorm16 pack).
      // Test-local arenas keep the replay's gated peaks untouched; a
      // non-hankin final leg lands on the clean endpoint, whose full-precision
      // classification the closing compile below provides.
      if (pin_final && k + 1 == count && hankin_step) {
        Arena pa(morph_target_buf, sizeof(morph_target_buf));
        Arena pb(morph_temp_buf, sizeof(morph_temp_buf));
        Arena pc(morph_aux_buf, sizeof(morph_aux_buf));
        CompiledHankin hk;
        MeshOps::compile_hankin(cur, hk, pa, pb);
        PolyMesh full;
        MeshOps::update_hankin(hk, full, pa,
                               std::max(steps[k].param, OpLeg::THETA_EPS));
        MeshOps::classify_faces_by_topology(full, pb, pc, pa);
        full_topo.assign(full.topology.begin(), full.topology.end());
      }

      // Mirror finish_build_leg's boundary compaction: the finished leg's
      // transients are reclaimed and only the endpoint the next leg sweeps
      // from crosses the reset. Without it the replay accumulates every leg
      // and reports peaks the effect never reaches. The last leg's landing is
      // consumed by the final swap below, so it is left standing exactly as
      // the effect leaves it for the next shape's compaction.
      if (k + 1 < count) {
        Persist<PolyMesh> pn(next, scratch_arena_b, persistent_arena);
        cur = PolyMesh();
        seed_slot = MeshState();
        prev_landing = nullptr;
        persistent_arena.reset();
        bank.bake_all(persistent_arena);
      }
      cur = std::move(next);
    }

    // Final sprite handoff, mirroring finish_build: the finished solid's
    // per-face palettes are the last landing's landed palettes, snapshotted
    // before the closing compaction. The compiled slot must match them face
    // for face, so the sprite's first frame draws exactly what the last leg
    // frame drew (the closing w = 1 plateau).
    const size_t landed_faces = cur.face_counts.size();
    HS_EXPECT_LE(landed_faces, prev_landing->faces);
    HS_EXPECT_LE(landed_faces, MAX_FACES);
    if (landed_faces > MAX_FACES)
      return peaks;
    uint8_t sprite_pal[MAX_FACES];
    for (size_t f = 0; f < landed_faces; ++f)
      sprite_pal[f] = prev_landing->landed_palette(f);
    prev_landing = nullptr;
    MeshState final_slot;
    {
      ScratchScope a_guard(scratch_arena_a);
      {
        ScratchScope seed_guard(scratch_arena_b);
        PolyMesh built;
        MeshOps::clone(cur, built, scratch_arena_b);
        cur = PolyMesh();
        seed_slot = MeshState();
        persistent_arena.reset();
        bank.bake_all(persistent_arena);
        MeshOps::compile(built, final_slot, persistent_arena, scratch_arena_a);
      }
      ScratchScope b_guard(scratch_arena_b);
      MeshOps::classify_faces_by_topology(final_slot, scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }
    if (supported) {
      HS_EXPECT_EQ(landed_faces, final_slot.topology.size());
      int first_broken_carry = -1;
      for (size_t f = 0; f < landed_faces; ++f)
        if ((int)sprite_pal[f] != (int)carried_to[f] && first_broken_carry < 0)
          first_broken_carry = (int)f;
      HS_EXPECT_EQ(first_broken_carry, -1);
    }
    if (pin_final && full_topo.empty())
      full_topo.assign(final_slot.topology.begin(), final_slot.topology.end());

    // Variety: distinct palettes on the finished shape.
    peaks.palettes = 0;
    {
      bool seen[OpLeg::PALETTES] = {};
      for (size_t f = 0; f < landed_faces; ++f)
        if (sprite_pal[f] < OpLeg::PALETTES)
          seen[sprite_pal[f]] = true;
      for (bool s : seen)
        peaks.palettes += s;
    }

    // Symmetry pin: on the finished shape any two faces with the same
    // full-precision final class wear the same palette — the landed palette
    // is a function of the final classification alone. A
    // quantized-classification regression shatters classes and breaks this.
    if (pin_final) {
      HS_EXPECT_EQ(full_topo.size(), landed_faces);
      int asymmetric = 0;
      if (full_topo.size() == landed_faces) {
        for (size_t f = 0; f < landed_faces; ++f)
          for (size_t g = f + 1; g < landed_faces; ++g)
            if (full_topo[f] == full_topo[g] && sprite_pal[f] != sprite_pal[g])
              ++asymmetric;
      }
      HS_EXPECT_EQ(asymmetric, 0);
    }

    peaks.persistent = persistent_arena.get_high_water_mark();
    peaks.scratch_a = scratch_arena_a.get_high_water_mark();
    peaks.scratch_b = scratch_arena_b.get_high_water_mark();
    peaks.faces = final_slot.face_counts.size();
    if (gate) {
      std::printf("  [chain] %s: %zu legs, %d/%d palettes, final leg %d "
                  "classes, %d/%d blend pairs, persistent=%zu B / %zu B, "
                  "scratch a=%zu B / %zu B, b=%zu B / %zu B\n",
                  name, count, peaks.palettes, OpLeg::PALETTES,
                  peaks.final_classes, peaks.blend_pairs,
                  OpLeg::MAX_BLEND_PAIRS, peaks.persistent,
                  (size_t)ISLAMIC_PERSISTENT_BUDGET, peaks.scratch_a,
                  (size_t)ISLAMIC_SCRATCH_A_BUDGET, peaks.scratch_b,
                  (size_t)ISLAMIC_SCRATCH_B_BUDGET);
      HS_EXPECT_LE(peaks.persistent, ISLAMIC_PERSISTENT_BUDGET);
      HS_EXPECT_LE(peaks.scratch_a, ISLAMIC_SCRATCH_A_BUDGET);
      HS_EXPECT_LE(peaks.scratch_b, ISLAMIC_SCRATCH_B_BUDGET);

      if (hs_test::stats().failed != failed_before)
        std::printf("    [chain] %s FAILED\n", name);
    }
  }
  return peaks;
}

/** simple_registry seed indices the partition chains build from. */
inline constexpr uint8_t SEED_CUBE = 1;
inline constexpr uint8_t SEED_ICOSAHEDRON = 4;
static_assert(std::string_view(Solids::simple_registry[SEED_CUBE].name) ==
              "cube");
static_assert(
    std::string_view(Solids::simple_registry[SEED_ICOSAHEDRON].name) ==
    "icosahedron");

/** Partition chains: a lone kis, a gated leg departing a gated leg, and a
 * gated leg departing a swept one. No registry entry carries a partition
 * recipe yet, so these stand in for the recipe table. */
inline constexpr Solids::OpStep CHAIN_KIS[] = {{Solids::Op::KIS}};
inline constexpr Solids::OpStep CHAIN_KIS_DUAL[] = {{Solids::Op::KIS},
                                                    {Solids::Op::DUAL}};
inline constexpr Solids::OpStep CHAIN_AMBO_DUAL[] = {{Solids::Op::AMBO},
                                                     {Solids::Op::DUAL}};
inline constexpr Solids::OpStep CHAIN_HK62_DUAL[] = {
    {Solids::Op::HANKIN, 62.0f * Solids::IslamicStarPatterns::D2R},
    {Solids::Op::DUAL}};
inline constexpr Solids::Recipe DODECAHEDRON_KIS_RECIPE = {
    Solids::SEED_DODECAHEDRON, CHAIN_KIS, std::size(CHAIN_KIS)};
inline constexpr Solids::Recipe CUBE_KIS_DUAL_RECIPE = {
    SEED_CUBE, CHAIN_KIS_DUAL, std::size(CHAIN_KIS_DUAL)};
inline constexpr Solids::Recipe ICOSAHEDRON_AMBO_DUAL_RECIPE = {
    SEED_ICOSAHEDRON, CHAIN_AMBO_DUAL, std::size(CHAIN_AMBO_DUAL)};
inline constexpr Solids::Recipe DODECAHEDRON_HK62_DUAL_RECIPE = {
    Solids::SEED_DODECAHEDRON, CHAIN_HK62_DUAL, std::size(CHAIN_HK62_DUAL)};

/**
 * @brief Replays every registry recipe plus the partition chains leg by leg.
 */
inline void test_recipe_chain_build_replay() {
  int chains = 0;
  for (const Solids::Entry &entry : Solids::Collections::get_islamic_solids()) {
    if (!entry.recipe)
      continue;
    ++chains;
    replay_build_chain(entry.name, *entry.recipe);
  }
  // Guards against a vacuous pass if the recipe pointers are dropped.
  HS_EXPECT_GE(chains, 2);

  replay_build_chain("dodecahedron_kis", DODECAHEDRON_KIS_RECIPE);
  replay_build_chain("cube_kis_dual", CUBE_KIS_DUAL_RECIPE);
  replay_build_chain("icosahedron_ambo_dual", ICOSAHEDRON_AMBO_DUAL_RECIPE);
  replay_build_chain("dodecahedron_hk62_dual", DODECAHEDRON_HK62_DUAL_RECIPE);
}

// ---------------------------------------------------------------------------
// Smooth kis/needle reconcile: the identity meshes (dt / dtd) and the authored
// kis/needle meshes are topology-exact (identical V/E/F/I); the residual gap is
// closed by a per-vertex great-circle slerp along the nearest-vertex bijection.
// Gates that the bijection is injective on every affected seed, the residual is
// small, and the identity truncate depth (t = 1/3) is what the effect emits.
// ---------------------------------------------------------------------------

/** One smooth-kis/needle macro site: how X is reached, and which identity vs
 * authored meshes the reconcile leg spans. */
struct ReconcileSite {
  const char *name;
  const Solids::Recipe *recipe;
  bool dtd; /**< true = standalone kis (dtd); false = trailing dual,kis (dt). */
};

inline const ReconcileSite RECONCILE_SITES[] = {
    {"needle", &Solids::TRUNCATED_ICOSAHEDRON_AMBO_RELAX100_HK54_NEEDLE_RECIPE,
     false},
    {"gyro_kis", &Solids::TRUNCATED_OCTAHEDRON_GYRO_KIS_HK17_RECIPE, false},
    {"icosahedron_kis_gyro", &Solids::ICOSAHEDRON_KIS_GYRO_RECIPE, true},
};

/** The identity truncate depth the smooth path sweeps to (the "uniform" Conway
 * depth); IslamicStars::MACRO_TRUNCATE_T mirrors this. */
inline constexpr float RECONCILE_TRUNCATE_T = 1.0f / 3.0f;

/**
 * @brief Verifies the Conway-identity reconcile is well-posed on every smooth
 *        kis/needle seed: the identity mesh (dt/dtd) and the authored kis/needle
 *        mesh share V/E/F, the nearest-vertex map is a bijection, and the
 *        residual chord the reconcile slerp closes is bounded.
 */
inline void test_reconcile_bijection_wellposed() {
  // Residual well under half the vertex spacing on the densest seed; a bijection
  // failure would trip the injectivity check first, this bounds the slerp arc.
  constexpr float MAX_RESIDUAL_CHORD = 0.12f;
  for (const ReconcileSite &site : RECONCILE_SITES) {
    const int failed_before = hs_test::stats().failed;
    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
    Arena keep(morph_persist_buf, sizeof(morph_persist_buf));

    Solids::OpStep lowered[8];
    const size_t count = Solids::expand_to_primitives(*site.recipe, lowered, 8);
    // Locate the kis and its preceding dual (dt pair) or its standalone form.
    size_t kis = count;
    for (size_t k = 0; k < count; ++k)
      if (lowered[k].op == Solids::Op::KIS) {
        kis = k;
        break;
      }
    HS_EXPECT_LT(kis, count);
    const size_t x_prefix = site.dtd ? kis : kis - 1; // steps to reach X

    PolyMesh X =
        Solids::build_steps(site.recipe->seed, lowered, x_prefix, a, b);
    PolyMesh Xc; // stable copy: identity/authored reuse a,b as ping-pong
    MeshOps::clone(X, Xc, aux);

    // identity: dt = dual(truncate(X)); dtd = dual(truncate(dual(X))). Clone the
    // identity out of a/b before building authored, which reuses those arenas.
    PolyMesh identity;
    PolyMesh authored;
    if (site.dtd) {
      PolyMesh dX = MeshOps::dual(Xc, a, b);
      PolyMesh dXc;
      MeshOps::clone(dX, dXc, aux);
      PolyMesh t = MeshOps::truncate(dXc, a, b, RECONCILE_TRUNCATE_T);
      MeshOps::clone(MeshOps::dual(t, b, a), identity, keep);
      authored = MeshOps::kis(Xc, a, b); // kis(X)
    } else {
      PolyMesh t = MeshOps::truncate(Xc, a, b, RECONCILE_TRUNCATE_T);
      MeshOps::clone(MeshOps::dual(t, b, a), identity, keep);
      authored = MeshOps::needle(Xc, a, b); // kis(dual(X))
    }

    const size_t V = identity.vertices.size();
    HS_EXPECT_EQ(V, authored.vertices.size());
    HS_EXPECT_EQ(identity.face_counts.size(), authored.face_counts.size());
    HS_EXPECT_EQ(identity.faces.size(), authored.faces.size());

    // Greedy nearest-vertex map identity -> authored; must be a bijection.
    std::vector<bool> used(V, false);
    std::vector<int> match(V, -1);
    float worst_chord = 0.0f;
    bool injective = true;
    for (size_t i = 0; i < V; ++i) {
      int best = -1;
      float best_dot = -2.0f;
      for (size_t j = 0; j < V; ++j) {
        const float d = dot(identity.vertices[i], authored.vertices[j]);
        if (d > best_dot) {
          best_dot = d;
          best = static_cast<int>(j);
        }
      }
      if (best < 0 || used[best])
        injective = false;
      else
        used[best] = true;
      match[i] = best;
      if (best >= 0)
        worst_chord =
            std::max(worst_chord, distance_between(identity.vertices[i],
                                                   authored.vertices[best]));
    }
    HS_EXPECT_TRUE(injective);
    HS_EXPECT_LT(worst_chord, MAX_RESIDUAL_CHORD);

    if (hs_test::stats().failed != failed_before)
      std::printf("    [reconcile] %s FAILED (V=%zu inj=%d worst_chord=%.4f)\n",
                  site.name, V, injective ? 1 : 0, (double)worst_chord);
    else
      std::printf("  [reconcile] %s: V=%zu F=%zu bijective worst_chord=%.4f\n",
                  site.name, V, identity.face_counts.size(),
                  (double)worst_chord);
  }
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

/**
 * @brief Runs all OpLeg operator-level tests.
 * @return The module's failure count.
 */
inline int run_conway_morph_tests() {
  hs_test::ModuleFixture fixture("conway_morph");

  test_edge_endpoints_match_registry();
  test_edge_sweeps_hold_topology();
  test_edge_morph_frames_fit_scratch_budget();

  test_relax_is_vertex_order_identity();

  test_snub_tetrahedron_relax_converges_to_icosahedron();
  test_ambo_tetrahedron_is_regular_octahedron();

  test_jitterbug_icosa_point_is_regular();
  test_jitterbug_octa_end_covers_octahedron();
  test_jitterbug_sweep_holds_topology();

  test_truncate_near_half_merges_onto_ambo();
  test_ops_at_t_eps_primary_faces_match_seed();

  test_ambo_leg_on_hankin_seed_holds_topology();
  test_hankin_sweep_on_islamic_seeds_holds_topology();
  test_hankin_sweep_vertex_stability();
  test_opleg_hankin_sweep_smoke();

  test_truncate_leg_on_recipe_seeds_holds_topology();
  test_snub_leg_on_recipe_seeds_holds_topology();
  test_relax_leg_on_recipe_seeds_holds_topology();
  test_medial_dual_bridge_wellformed();
  test_opleg_medial_leg_smoke();
  test_opleg_dual_bridge_seam_correspondence();
  test_opleg_step_leg_smoke();
  test_opleg_step_paused_holds_frame();
  test_opleg_step_leg_overshooting_easing();
  test_opleg_gated_swap_smoke();
  test_opleg_edge_leg_crossfade();
  test_unsweepable_recipe_steps_are_gated();

  test_recipe_chain_build_replay();
  test_reconcile_bijection_wellposed();

  test_walk_policy_coverage_and_balance();
  test_ordered_tour_full_coverage_and_wrap();

  return fixture.result();
}

} // namespace conway_morph_tests
} // namespace hs_test
