/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Operator-level tests for the ConwayMorph transition design
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
#include "core/mesh/conway.h"
#include "core/mesh/conway_graph.h"
#include "core/mesh/solids.h"
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
 * @brief Verifies the jitterbug leg exactly as ConwayMorph runs it — t from
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
  return Solids::SolidBuilder(Solids::Archimedean::rhombicuboctahedron(a, b),
                              a, b)
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
// Runner
// ---------------------------------------------------------------------------

/**
 * @brief Runs all ConwayMorph operator-level tests.
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

  test_walk_policy_coverage_and_balance();
  test_ordered_tour_full_coverage_and_wrap();

  return fixture.result();
}

} // namespace conway_morph_tests
} // namespace hs_test
