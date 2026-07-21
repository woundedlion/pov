/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
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
#include "tests/mesh_test_util.h"
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
// §7.1 Endpoint exactness: sweeping to an edge endpoint arrives at the
// registry generator's output — exactly where the composition is the registry
// chain (same code path, same seed frame), within geometric tolerance where it
// is not (dual-family ambo arrivals, bridge arrivals, and t = 0 ends, which
// emit expanded topology).
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
};

/**
 * @brief Comparison regime of an edge's to_node end.
 * @param e Edge to classify.
 * @return EXACT when op(seed, t_to) [+ relax] is the to_node registry chain;
 *         VERTEX_MATCH for arrivals off the registry seed (dual-family ambo,
 *         non-settle bridges); REGULAR for the settling bridge, whose relax
 *         orientation tracks the seed frame, not the registry icosahedron;
 *         PAIR_COVER for the jitterbug bridge's collapsed octa end.
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
 *        the t = 0 ends and the off-registry-seed arrivals.
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
        const float tol =
            e.op == ConwayGraph::MorphOp::TRUNCATE ? 0.08f : 0.06f;
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
        HS_EXPECT_TRUE(v0 > 0 && f0 > 0 && i0 > 0);
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

constexpr size_t MORPH_SCRATCH_A_BUDGET =
    24 * 1024; /**< HankinSolids scratch_a split. */
constexpr size_t MORPH_SCRATCH_B_BUDGET =
    32 * 1024; /**< HankinSolids scratch_b split. */

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

    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
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
      const int next = edge_other_end(e, node);
      in_family = family(next) != family(node) ? 0 : in_family + 1;
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

    // Seed reconciliation must be legal, and the two non-KEEP fixes only
    // fire at the nodes the effect's HS_CHECKs allow them at.
    const SeedFix fix = seed_fix_at_start(e, held);
    HS_EXPECT_TRUE(fix != SeedFix::INVALID);
    if (fix == SeedFix::DUAL_SWAP) {
      HS_EXPECT_TRUE(node == CUBOCTAHEDRON || node == ICOSIDODECAHEDRON);
      held = dual_platonic(held);
    } else if (fix == SeedFix::REGEN_TETRA) {
      HS_EXPECT_TRUE(node == OCTAHEDRON || node == ICOSAHEDRON);
      held = TETRAHEDRON;
    }

    const bool reverse = EDGES[e].to_node == node;
    const int arrived = reverse ? EDGES[e].from_node : EDGES[e].to_node;
    // ADOPT stores every platonic arrival as the held seed, reverse legs
    // included (the jitterbug's octa -> icosa leg adopts the canonical
    // icosahedron; reverse tetra arrivals re-adopt the already-held tetra).
    if (EDGES[e].reseed == Reseed::ADOPT && is_platonic(arrived))
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
        HS_EXPECT_TRUE(v0 > 0 && f0 > 0 && i0 > 0);
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
      const float theta = THETA_EPS + (site.theta_star - THETA_EPS) *
                                          (static_cast<float>(s) /
                                           (SAMPLES - 1));
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
        HS_EXPECT_TRUE(v0 > 0 && f0 > 0 && i0 > 0);
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
      std::printf(
          "  [hankin-sweep] %s: F=%zu compiled=%zu across %d samples\n",
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
};

/** @brief One dynamic vertex's solved position and the branch that made it. */
struct HankinSolve {
  Vector pos;
  HankinBranch branch;
  /** dist^2(star, corner) / max(dist^2(m, corner)); the STAR_FAR_RATIO_SQ
   * guard fires above 36. Zero on the non-intersect branches. */
  float far_ratio = 0;
};

/**
 * @brief Mirrors MeshOps::update_hankin's per-vertex solve, exposing the
 *        branch each dynamic vertex takes.
 * @param compiled Baked hankin topology.
 * @param angle Contact angle in radians.
 * @param out Filled with one entry per dynamic instruction.
 */
inline void hankin_solve(const CompiledHankin &compiled, float angle,
                         std::vector<HankinSolve> &out) {
  const bool is_flat = std::abs(angle) < math::TOLERANCE;
  const float cos_ha = cosf(angle * 0.5f);
  const float sin_ha = sinf(angle * 0.5f);

  out.assign(compiled.dynamic_instructions.size(), HankinSolve{});
  for (size_t i = 0; i < compiled.dynamic_instructions.size(); ++i) {
    const HankinInstruction &instr = compiled.dynamic_instructions[i];
    const Vector p_corner = compiled.base_vertices[instr.v_corner];
    const Vector cn = normalized_or(p_corner, p_corner);

    if (is_flat) {
      out[i] = {cn, HankinBranch::COLLAPSED};
      continue;
    }

    const Vector m1 = compiled.static_vertices[instr.idx_m1];
    const Vector m2 = compiled.static_vertices[instr.idx_m2];
    const Vector cross1 = cross(compiled.base_vertices[instr.v_prev], p_corner);
    const Vector cross2 = cross(p_corner, compiled.base_vertices[instr.v_next]);
    if (dot(cross1, cross1) < math::EPS_CROSS_SQ ||
        dot(cross2, cross2) < math::EPS_CROSS_SQ) {
      out[i] = {cn, HankinBranch::COLLAPSED};
      continue;
    }

    const Quaternion q1(cos_ha, sin_ha * m1.x, sin_ha * m1.y, sin_ha * m1.z);
    const Quaternion q2(cos_ha, -sin_ha * m2.x, -sin_ha * m2.y, -sin_ha * m2.z);
    Vector intersect = cross(rotate(cross1.normalized(), q1),
                             rotate(cross2.normalized(), q2));

    bool degenerate = dot(intersect, intersect) < math::EPS_LEN_SQ;
    float far_ratio = 0;
    if (!degenerate) {
      if (dot(intersect, p_corner) < 0)
        intersect = -intersect;
      intersect = intersect.normalized();
      const float local_sq =
          std::max(distance_squared(m1, cn), distance_squared(m2, cn));
      far_ratio = distance_squared(intersect, cn) / local_sq;
      degenerate = far_ratio > MeshOps::STAR_FAR_RATIO_SQ;
    }
    if (degenerate) {
      Vector fallback = normalized_or(m1 + m2, cn);
      if (dot(fallback, p_corner) < 0)
        fallback = -fallback;
      out[i] = {fallback, HankinBranch::FALLBACK, far_ratio};
      continue;
    }
    out[i] = {intersect, HankinBranch::INTERSECT, far_ratio};
  }
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
    return idx < compiled.static_offset
               ? compiled.static_vertices[idx]
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
  float theta = 0;      /**< Contact angle of this sample, radians. */
  int fallback = 0;     /**< Dynamic vertices on the FALLBACK branch. */
  int branch_flips = 0; /**< Vertices whose branch changed vs the previous step. */
  float max_disp = 0;   /**< Largest single-vertex chord vs the previous step. */
  float mean_disp = 0;  /**< Mean dynamic-vertex chord vs the previous step. */
  int normal_flips = 0; /**< Non-degenerate faces whose Newell normal reversed. */
  int flat_faces = 0;   /**< Faces below HANKIN_FLAT_FACE this step. */
  float max_far_ratio = 0;    /**< Largest far_ratio this step (guard fires at 36). */
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
    if (curr[i].branch == HankinBranch::FALLBACK)
      ++stats.fallback;
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
  int worst_step = 0;    /**< Step index owning worst_max_disp. */
  float spike_ratio = 0; /**< worst_max_disp / mean_disp at that step. */
  int worst_flat_faces = 0;    /**< Most sub-HANKIN_FLAT_FACE faces in a step. */
  float worst_far_ratio = 0;   /**< Largest far_ratio over the sweep. */
  float worst_corner_chord = 0;/**< Largest chord(star point, corner) over the sweep. */
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

/** @brief Prints one per-step table row set. */
inline void hankin_print_table(const char *label,
                               const std::vector<HankinStepStats> &table) {
  std::printf("      %s: s theta_deg  fb flip  max_disp mean_disp nflip flat "
              "far_ratio  corner\n",
              label);
  for (size_t s = 0; s < table.size(); ++s) {
    const HankinStepStats &r = table[s];
    std::printf("        %16s%2zu %9.3f %3d %4d %9.5f %9.5f %5d %4d %9.2f "
                "%7.4f\n",
                "", s, r.theta * 180.0f / PI_F, r.fallback, r.branch_flips,
                r.max_disp, r.mean_disp, r.normal_flips, r.flat_faces,
                r.max_far_ratio, r.max_corner_chord);
  }
}

/**
 * @brief Measures per-frame sweep stability of the four Phase-1 hankin legs
 *        under the shipping re-solve (uniform and eased theta) and under a
 *        slerp-from-corner parameterization.
 * @details Reports fallback-branch population, branch flips, per-step vertex
 * chords and face-normal reversals; asserts the slerp path is branch-flip and
 * normal-flip free and that its endpoints reproduce the collapsed form and the
 * theta_star solve.
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

    // Guard headroom: the far-star guard is scaled by the corner's local edge
    // scale, so on coarse seeds 36 * local_sq can exceed the 4.0 maximum
    // squared chord, making the guard unreachable for that vertex.
    int unreachable = 0;
    float max_local_sq = 0;
    for (size_t i = 0; i < compiled.dynamic_instructions.size(); ++i) {
      const HankinInstruction &instr = compiled.dynamic_instructions[i];
      const Vector cn = normalized_or(compiled.base_vertices[instr.v_corner],
                                      compiled.base_vertices[instr.v_corner]);
      const float local_sq =
          std::max(distance_squared(compiled.static_vertices[instr.idx_m1], cn),
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

    hankin_print_table("resolve-eased", tables[1]);
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
    HS_EXPECT_EQ(slerp_sum.total_branch_flips, 0);
    HS_EXPECT_EQ(slerp_sum.total_normal_flips, 0);

    // Opening bookend: chord between the collapsed form and the leg's first
    // drawn angle, in sphere radii (sub-pixel iff below ~1/display radius).
    float eps_chord = 0;
    std::vector<HankinSolve> at_eps;
    hankin_solve(compiled, THETA_EPS, at_eps);
    for (size_t i = 0; i < at_eps.size(); ++i)
      eps_chord = std::max(eps_chord, (at_eps[i].pos - collapsed[i].pos).magnitude());
    std::printf("      theta_eps=%.3f opening chord max=%.6f radii "
                "(%.3f px at r=64)\n",
                THETA_EPS, eps_chord, eps_chord * 64.0f);

    // The slerp path's k = 0 form is the exact collapse, where every rosette
    // face has zero area; a K_EPS floor is its THETA_EPS analog. Report the
    // smallest k that lifts every face above HANKIN_FLAT_FACE.
    float k_eps = 1.0f;
    std::vector<HankinSolve> probe(arrival.size());
    std::vector<Vector> probe_normals;
    for (int q = 1; q <= 200; ++q) {
      const float k = static_cast<float>(q) / 200.0f;
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
  uint8_t pal[16], sides[16];
  for (size_t f = 0; f < dodeca.face_counts.size(); ++f) {
    pal[f] = static_cast<uint8_t>(f % Animation::OpLeg::PALETTES);
    sides[f] = dodeca.face_counts[f];
  }

  Animation::OpLeg::PaletteHandoff handoff{
      &bank.bank, pal, sides, dodeca.face_counts.size(), false};

  // Per-frame motion bound: growing star points out from their corners keeps
  // every step small and unimodal. Re-solving the contact-plane intersection
  // per frame instead sends star points on geodesic excursions (measured to
  // 1.84 chord on ambo-of-hankin seeds), which draws as lines crossing the
  // pattern; the bound is what stops that parameterization coming back.
  constexpr float MAX_STEP_CHORD = 0.15f;
  size_t drawn_frames = 0;
  size_t drawn_faces = 0;
  float worst_step = 0.0f;
  std::vector<Vector> prev_v;
  auto cb = [&](Canvas &, const MeshState &m,
                const Animation::OpLeg::Shading &sh) {
    HS_EXPECT_EQ(m.face_counts.size(), sh.faces);
    if (drawn_frames == 0)
      drawn_faces = sh.faces;
    else
      HS_EXPECT_EQ(sh.faces, drawn_faces);
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
    ++drawn_frames;
  };

  using Solids::IslamicStarPatterns::D2R;
  constexpr int SWEEP = 8;
  Animation::OpLeg anim(dodeca, Animation::OpLeg::THETA_EPS, 62.0f * D2R, leg,
                        cb, handoff, SWEEP);

  // 12 star faces (the base-face-order prefix) + 20 rosette faces.
  const Animation::OpLeg::Landing &landing = anim.landing();
  HS_EXPECT_EQ(landing.primary_faces, dodeca.face_counts.size());
  HS_EXPECT_EQ(landing.faces,
               dodeca.face_counts.size() + dodeca.vertices.size());
  HS_EXPECT_TRUE(landing.topology != nullptr);

  struct SmokeFx : public Effect {
    SmokeFx() : Effect(288, 144) {}
    void draw_frame() override {}
  };
  SmokeFx fx;
  for (int f = 0; f < SWEEP; ++f) {
    {
      Canvas c(fx);
      anim.step(c);
    }
    fx.advance_display();
  }
  HS_EXPECT_EQ(drawn_frames, (size_t)SWEEP);
  HS_EXPECT_EQ(drawn_faces, landing.faces);
  HS_EXPECT_LT(worst_step, MAX_STEP_CHORD);
  std::printf("  [opleg hankin] worst per-frame vertex step %.4f chord "
              "(bound %.2f)\n",
              (double)worst_step, (double)MAX_STEP_CHORD);
}

// ---------------------------------------------------------------------------
// Recipe chain build replay (docs/opchain_morph_spec.md sections 9.2/9.3):
// every registry entry with a non-null recipe is lowered and replayed leg by
// leg exactly as IslamicStars builds it — same handoff, bookend, pinned
// targets, and blend-range dance — stepping every leg frame by frame with a
// recording draw callback. Gates per-leg compiled-face-count constancy,
// from-palette continuity across leg boundaries, the final clean-swap
// classification, and the persistent/scratch high-water against IslamicStars'
// configured split.
// ---------------------------------------------------------------------------

constexpr size_t ISLAMIC_SCRATCH_A_BUDGET =
    114 * 1024; /**< IslamicStars scratch_a split. */
constexpr size_t ISLAMIC_SCRATCH_B_BUDGET =
    80 * 1024; /**< IslamicStars scratch_b split. */
/** Device persistent budget of IslamicStars' arena split. */
constexpr size_t ISLAMIC_PERSISTENT_BUDGET =
    DEVICE_GLOBAL_ARENA_SIZE - ISLAMIC_SCRATCH_A_BUDGET -
    ISLAMIC_SCRATCH_B_BUDGET;

/**
 * @brief Replays every non-null recipe leg by leg as IslamicStars builds it
 *        and gates continuity, classification, and the arena high-waters.
 */
inline void test_recipe_chain_build_replay() {
  using Animation::OpLeg;
  constexpr int HANKIN_LEG_FRAMES = 32, AMBO_LEG_FRAMES = 24;
  constexpr size_t MAX_FACES = 256;
  constexpr size_t MAX_STEPS = 8;

  struct ChainFx : public Effect {
    ChainFx() : Effect(288, 144) {}
    void draw_frame() override {}
  };

  int chains = 0;
  for (const Solids::Entry &entry : Solids::Collections::get_islamic_solids()) {
    if (!entry.recipe)
      continue;
    ++chains;
    const Solids::Recipe &recipe = *entry.recipe;
    const int failed_before = hs_test::stats().failed;

    reset_globals();
    configure_arenas(GLOBAL_ARENA_SIZE - ISLAMIC_SCRATCH_A_BUDGET -
                         ISLAMIC_SCRATCH_B_BUDGET,
                     ISLAMIC_SCRATCH_A_BUDGET, ISLAMIC_SCRATCH_B_BUDGET);
    hs::random().seed(2026u);

    MeshPaletteBank bank;
    bank.bake_all(persistent_arena);

    // Seed solid: persistent PolyMesh plus the compiled/classified carousel
    // slot, exactly as spawn_shape prepares them.
    PolyMesh cur;
    MeshState seed_slot;
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
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

    std::array<int, OpLeg::PALETTES> slots;
    MeshPaletteBank::shuffle_indices(slots);
    std::array<uint8_t, OpLeg::PALETTES> targets;
    for (int i = 0; i < OpLeg::PALETTES; ++i)
      targets[i] = static_cast<uint8_t>(i);
    hs::shuffle(targets.begin(), targets.end());

    Solids::OpStep steps[MAX_STEPS];
    const size_t count = Solids::expand_to_primitives(recipe, steps, MAX_STEPS);
    HS_EXPECT_GT(count, (size_t)0);
    bool supported = true;
    int leg_frames[MAX_STEPS] = {};
    int total_frames = 0;
    for (size_t k = 0; k < count; ++k) {
      if (steps[k].op != Solids::Op::HANKIN && steps[k].op != Solids::Op::AMBO)
        supported = false;
      leg_frames[k] = steps[k].op == Solids::Op::HANKIN ? HANKIN_LEG_FRAMES
                                                        : AMBO_LEG_FRAMES;
      total_frames += leg_frames[k];
    }
    HS_EXPECT_TRUE(supported);
    if (!supported) {
      std::printf("    [chain] %s: unsupported lowered op\n", entry.name);
      continue;
    }

    // Finished solid's classification: every leg keys its targets on a prefix
    // of it, so a surviving face's target is fixed for the whole chain.
    int final_topo_buf[MAX_FACES] = {};
    size_t final_faces = 0;
    {
      ScratchScope a_guard(scratch_arena_a);
      ScratchScope b_guard(scratch_arena_b);
      PolyMesh fin =
          Solids::build_recipe(recipe, scratch_arena_a, scratch_arena_b);
      MeshState fin_state;
      MeshOps::compile(fin, fin_state, scratch_arena_a, scratch_arena_b);
      MeshOps::classify_faces_by_topology(fin_state, scratch_arena_a,
                                          scratch_arena_b, scratch_arena_a);
      final_faces = fin_state.topology.size();
      HS_EXPECT_LE(final_faces, MAX_FACES);
      for (size_t f = 0; f < final_faces; ++f)
        final_topo_buf[f] = fin_state.topology[f];
    }
    const int *final_topo = final_topo_buf;

    ChainFx fx;
    uint8_t prev_pal_buf[MAX_FACES];
    Vector prev_centroid[MAX_FACES];
    uint8_t carried_from[MAX_FACES] = {};
    uint8_t cur_target[MAX_FACES] = {};
    uint8_t prev_target[MAX_FACES] = {};
    const OpLeg::Landing *prev_landing = nullptr;
    PolyMesh next;
    float c_done = 0.0f;

    for (size_t k = 0; k < count; ++k) {
      const size_t prev_faces = cur.face_counts.size();
      HS_EXPECT_LE(prev_faces, MAX_FACES);
      {
        size_t off = 0;
        for (size_t f = 0; f < prev_faces; ++f) {
          Vector c(0.0f, 0.0f, 0.0f);
          const int n = cur.face_counts[f];
          for (int j = 0; j < n; ++j)
            c = c + cur.vertices[cur.faces[off + j]];
          prev_centroid[f] = c.normalized();
          off += n;
        }
      }
      const uint8_t *prev_pal;
      if (k == 0) {
        for (size_t f = 0; f < prev_faces; ++f)
          prev_pal_buf[f] = static_cast<uint8_t>(
              slots[wrap(seed_slot.topology[f], OpLeg::PALETTES)]);
        prev_pal = prev_pal_buf;
      } else {
        HS_EXPECT_EQ(prev_landing->faces, prev_faces);
        prev_pal = prev_landing->from_palette;
      }

      // Eager clean endpoint, exactly as start_build_leg derives it (for the
      // ambo leg: the AMBO mesh, never the swept truncate form).
      generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
        PolyMesh nx = steps[k].op == Solids::Op::AMBO
                          ? MeshOps::ambo(cur, a, b)
                          : MeshOps::hankin(cur, a, b, steps[k].param);
        next = Solids::finalize_solid(nx, target);
      });
      const size_t bookend_faces = next.face_counts.size();
      HS_EXPECT_LE(bookend_faces, MAX_FACES);
      HS_EXPECT_LE(bookend_faces, final_faces);

      OpLeg::PaletteHandoff handoff{&bank.bank, prev_pal, nullptr,
                                    prev_faces, false,    prev_centroid,
                                    &targets};
      OpLeg::BookendClasses bookend{final_topo, bookend_faces};

      size_t drawn = 0;
      size_t leg_faces = 0;
      auto cb = [&](Canvas &, const MeshState &m, const OpLeg::Shading &sh) {
        HS_EXPECT_EQ(m.face_counts.size(), sh.faces);
        if (drawn == 0)
          leg_faces = sh.faces;
        else
          HS_EXPECT_EQ(sh.faces, leg_faces);
        ++drawn;
      };

      const float c_hi = k + 1 == count
                             ? 1.0f
                             : c_done + static_cast<float>(leg_frames[k]) /
                                            static_cast<float>(total_frames);
      OpLeg leg =
          steps[k].op == Solids::Op::HANKIN
              ? OpLeg(cur, 0.0f, steps[k].param, persistent_arena, cb, handoff,
                      leg_frames[k], bookend)
              : OpLeg(cur, ConwayGraph::MorphOp::TRUNCATE, 0.0f, 0.5f, 0.0f,
                      0.0f, persistent_arena, cb, handoff, leg_frames[k],
                      bookend);
      leg.set_blend_range(c_done, c_hi);
      const OpLeg::Landing &landing = leg.landing();

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

      // From-palette continuity: the emission prefix keeps the ORIGINAL from
      // ids, so every leg's pair set blends (original from -> pinned target)
      // with only the weight range advancing.
      for (size_t f = 0; f < prev_faces; ++f)
        HS_EXPECT_EQ((int)landing.from_palette[f],
                     (int)(k == 0 ? prev_pal_buf[f] : carried_from[f]));
      HS_EXPECT_LE(landing.faces, MAX_FACES);
      for (size_t f = 0; f < landing.faces; ++f)
        carried_from[f] = landing.from_palette[f];

      // Target continuity: a surviving face's target palette must not move at
      // a leg boundary. The blend weight is mid-range there, so a moved target
      // jumps that face's color — the defect keying targets on each leg's own
      // arrival produced (32/32 and 62/62 faces on the dodecahedron chain).
      for (size_t f = 0; f < landing.faces; ++f)
        cur_target[f] =
            targets[wrap(landing.topology[f], OpLeg::PALETTES)];
      if (k > 0) {
        size_t moved = 0;
        for (size_t f = 0; f < prev_faces; ++f)
          if (cur_target[f] != prev_target[f])
            ++moved;
        if (moved)
          std::printf("    [chain target] %s leg %zu: %zu/%zu faces change "
                      "target at w=%.3f\n",
                      entry.name, k, moved, prev_faces, (double)c_done);
        HS_EXPECT_EQ(moved, (size_t)0);
      }
      std::memcpy(prev_target, cur_target, landing.faces);

      // Landing lives in the leg's arena-backed Transients; outlives `leg`.
      prev_landing = &landing;
      cur = std::move(next);
      c_done = c_hi;
    }

    // Final clean swap: the landing's bookend classification must match a
    // fresh classification of the compiled endpoint, and the landed targets
    // are the per-shape pinned set verbatim.
    MeshState final_slot;
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &) {
      MeshOps::compile(cur, final_slot, target, a);
    });
    {
      ScratchScope a_guard(scratch_arena_a);
      ScratchScope b_guard(scratch_arena_b);
      MeshOps::classify_faces_by_topology(final_slot, scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }
    HS_EXPECT_EQ(prev_landing->faces, final_slot.topology.size());
    for (size_t f = 0; f < prev_landing->faces; ++f)
      HS_EXPECT_EQ(prev_landing->topology[f], final_slot.topology[f]);
    for (int i = 0; i < OpLeg::PALETTES; ++i)
      HS_EXPECT_EQ((int)prev_landing->to_palette[i], (int)targets[i]);

    const size_t p_peak = persistent_arena.get_high_water_mark();
    const size_t a_peak = scratch_arena_a.get_high_water_mark();
    const size_t b_peak = scratch_arena_b.get_high_water_mark();
    std::printf("  [chain] %s: %zu legs, persistent=%zu B / %zu B, "
                "scratch a=%zu B / %zu B, b=%zu B / %zu B\n",
                entry.name, count, p_peak, (size_t)ISLAMIC_PERSISTENT_BUDGET,
                a_peak, (size_t)ISLAMIC_SCRATCH_A_BUDGET, b_peak,
                (size_t)ISLAMIC_SCRATCH_B_BUDGET);
    HS_EXPECT_LE(p_peak, ISLAMIC_PERSISTENT_BUDGET);
    HS_EXPECT_LE(a_peak, ISLAMIC_SCRATCH_A_BUDGET);
    HS_EXPECT_LE(b_peak, ISLAMIC_SCRATCH_B_BUDGET);

    if (hs_test::stats().failed != failed_before)
      std::printf("    [chain] %s FAILED\n", entry.name);
  }
  // Guards against a vacuous pass if the recipe pointers are dropped.
  HS_EXPECT_GE(chains, 2);
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
  test_recipe_chain_build_replay();

  test_walk_policy_coverage_and_balance();
  test_ordered_tour_full_coverage_and_wrap();

  return fixture.result();
}

} // namespace conway_morph_tests
} // namespace hs_test
