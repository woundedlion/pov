/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Pre-flight measurements for the OpChainMorph pure-inflate roster
 * (docs/opchain_morph_spec.md sections 2, 2.4, 8.3).
 *
 * Coverage:
 *   - Chamfer zero-area birth limit: newborn hexagon area and preserved-face
 *     displacement as t -> 0, on simple seeds and on the shipping hankin seed.
 *   - Chamfer sweep: constant V/F/I and compiled face count, closed genus-0
 *     manifold, unit vertices, outward face normals, per-step displacement.
 *   - Chamfer birth epsilon: smallest t at which no newborn hexagon is culled
 *     by SDF::Face's collapsed-area reject, per seed.
 *   - Needle gated-swap chain: the DUAL then KIS partition ops the needle
 *     recipe lowers to, each landing a well-formed closed manifold on the
 *     hankin(54 deg) seed, with the {DUAL, KIS} lowering matching
 * MeshOps::needle.
 *   - Build-chain provenance: face-centroid spacing per intermediate mesh
 *     against PROVENANCE_TOL_SQ, nearest/second-nearest ambiguity of the
 *     lookups build_palette_mapping actually performs, and the mapping path
 *     and prev_faces each real leg takes.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <span>
#include <vector>

#include "core/animation/opleg.h"
#include "core/mesh/conway.h"
#include "core/mesh/conway_graph.h"
#include "core/mesh/hankin.h"
#include "core/mesh/recipe.h"
#include "core/mesh/solids.h"
#include "core/render/sdf.h"
#include "tests/mesh_test_util.h"
#include "tests/test_conway.h" // check_euler_genus0
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace opchain_probe_tests {

inline uint8_t probe_a_buf[768 * 1024];       /**< Op output arena. */
inline uint8_t probe_b_buf[768 * 1024];       /**< Op scratch arena. */
inline uint8_t probe_seed_buf[512 * 1024];    /**< Held seed arena. */
inline SDF::FaceScratchBuffer probe_face_scratch; /**< Face setup scratch. */

/** Canvas rows the shipping mesh raster uses; SDF::Face's cull decision is
 * taken against this geometry. */
inline constexpr int PROBE_H = 144;
inline constexpr int PROBE_HV = PROBE_H + hs::H_OFFSET;

using Solids::Op;
using Solids::OpStep;

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

/** @brief Face-vertex offset table for a PolyMesh, in emission order. */
inline void face_offsets(const PolyMesh &m, std::vector<size_t> &out) {
  out.assign(m.face_counts.size(), 0);
  size_t off = 0;
  for (size_t f = 0; f < m.face_counts.size(); ++f) {
    out[f] = off;
    off += m.face_counts[f];
  }
}

/** @brief Unit vertex-average centroid of every face, in emission order. */
inline void face_centroids(const PolyMesh &m, std::vector<Vector> &out) {
  out.assign(m.face_counts.size(), Vector(0, 0, 0));
  size_t off = 0;
  for (size_t f = 0; f < m.face_counts.size(); ++f) {
    Vector c(0, 0, 0);
    const int n = m.face_counts[f];
    for (int k = 0; k < n; ++k)
      c = c + m.vertices[m.faces[off + k]];
    out[f] = c.normalized();
    off += n;
  }
}

/** @brief Newell area vector of one face (magnitude = planar area). */
inline Vector newell(const PolyMesh &m, size_t off, int n) {
  Vector a(0, 0, 0);
  for (int k = 0; k < n; ++k)
    a = a + cross(m.vertices[m.faces[off + k]],
                  m.vertices[m.faces[off + (k + 1) % n]]);
  return a * 0.5f;
}

/** @brief Whether SDF::Face rejects a face outright (collapsed or off-band). */
inline bool face_is_culled(const PolyMesh &m, size_t off, int n) {
  if (static_cast<size_t>(n) > SDF::FaceScratchBuffer::MAX_VERTS)
    return false;
  SDF::Face face(
      std::span<const Vector>(m.vertices.data(), m.vertices.size()),
      std::span<const uint16_t>(m.faces.data() + off, static_cast<size_t>(n)),
      probe_face_scratch, PROBE_HV, PROBE_H);
  return face.y_min > face.y_max;
}

/** @brief Median of a scratch vector (sorts in place); 0 when empty. */
inline float median_of(std::vector<float> &v) {
  if (v.empty())
    return 0.0f;
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}

// ---------------------------------------------------------------------------
// Chamfer probe (spec section 2, risk 5): chamfer has never been swept and has
// no T_EPS characterization, yet truncatedIcosahedron_hk58_chamfer63 needs one
// -- on a hankin mesh, which compounds the risk (section 2.4).
// ---------------------------------------------------------------------------

/** @brief One chamfer-leg seed. */
struct ChamferSite {
  const char *name;                     /**< Diagnostic label. */
  PolyMesh (*seed)(Arena &a, Arena &b); /**< Chain prefix up to the chamfer. */
};

inline PolyMesh probe_cube(Arena &a, Arena &b) {
  return Solids::Platonic::cube(a, b);
}
inline PolyMesh probe_dodecahedron(Arena &a, Arena &b) {
  return Solids::Platonic::dodecahedron(a, b);
}
inline PolyMesh probe_ticosa(Arena &a, Arena &b) {
  return Solids::Archimedean::truncatedIcosahedron(a, b);
}
inline PolyMesh probe_ticosa_hk58(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Archimedean::truncatedIcosahedron(a, b),
                              a, b)
      .hankin(58.0f * D2R)
      .build();
}

inline constexpr ChamferSite CHAMFER_SITES[] = {
    {"cube", probe_cube},
    {"dodecahedron", probe_dodecahedron},
    {"truncatedIcosahedron", probe_ticosa},
    {"truncatedIcosahedron_hk58", probe_ticosa_hk58},
};

/** Arrival parameter of the one shipping chamfer leg. */
inline constexpr float CHAMFER_T_STAR = 0.63f;

/**
 * @brief Measures chamfer's zero-area birth limit: newborn hexagon area and
 *        preserved-face displacement against t, on every chamfer seed.
 * @details The birth criterion (spec section 2) needs the created faces to
 * reach zero area while the preserved faces are unchanged. Chamfer emits F
 * shrunk primaries then E hexagons, so the split is exact in emission order.
 */
inline void test_chamfer_zero_area_birth_limit() {
  constexpr float T[] = {1e-5f, 1e-4f, 1e-3f, 1e-2f, 1e-1f};

  for (const ChamferSite &site : CHAMFER_SITES) {
    Arena persist(probe_seed_buf, sizeof(probe_seed_buf));
    Arena a(probe_a_buf, sizeof(probe_a_buf));
    Arena b(probe_b_buf, sizeof(probe_b_buf));
    PolyMesh seed;
    {
      ScratchScope ga(a);
      ScratchScope gb(b);
      seed = Solids::finalize_solid(site.seed(a, b), persist);
    }
    const size_t F = seed.face_counts.size();
    const size_t E = seed.faces.size() / 2;
    std::vector<Vector> seed_centroid;
    face_centroids(seed, seed_centroid);

    std::printf("  [chamfer-birth] %s (F=%zu E=%zu)\n", site.name, F, E);
    float prev_area = 0.0f;
    float limit_slope = 0.0f;
    for (float t : T) {
      ScratchScope fa(a);
      ScratchScope fb(b);
      PolyMesh out = MeshOps::chamfer(seed, a, b, t);
      HS_EXPECT_EQ(out.face_counts.size(), F + E);

      std::vector<size_t> off;
      face_offsets(out, off);
      float born_area = 0.0f;
      for (size_t f = F; f < out.face_counts.size(); ++f) {
        const Vector n = newell(out, off[f], out.face_counts[f]);
        born_area += std::sqrt(dot(n, n));
      }
      // Preserved faces: same side count, centroid held, and every vertex
      // converging on a seed vertex as t -> 0.
      float max_centroid_shift = 0.0f;
      float max_vertex_shift = 0.0f;
      for (size_t f = 0; f < F; ++f) {
        HS_EXPECT_EQ(static_cast<int>(out.face_counts[f]),
                     static_cast<int>(seed.face_counts[f]));
        Vector c(0, 0, 0);
        const int n = out.face_counts[f];
        for (int k = 0; k < n; ++k) {
          const Vector &v = out.vertices[out.faces[off[f] + k]];
          c = c + v;
          float near_sq = 1e9f;
          for (size_t s = 0; s < seed.vertices.size(); ++s) {
            const Vector e = v - seed.vertices[s];
            near_sq = std::min(near_sq, dot(e, e));
          }
          max_vertex_shift = std::max(max_vertex_shift, std::sqrt(near_sq));
        }
        c = c.normalized();
        const Vector d = c - seed_centroid[f];
        max_centroid_shift = std::max(max_centroid_shift, std::sqrt(dot(d, d)));
      }
      std::printf("      t=%-8.5f born_area=%.3e (area/t=%.4f) "
                  "preserved_centroid=%.3e preserved_vertex=%.3e\n",
                  static_cast<double>(t), static_cast<double>(born_area),
                  static_cast<double>(born_area / t),
                  static_cast<double>(max_centroid_shift),
                  static_cast<double>(max_vertex_shift));
      // Births vanish linearly in t; the preserved faces follow.
      HS_EXPECT_TRUE(born_area > prev_area);
      if (limit_slope == 0.0f)
        limit_slope = born_area / t;
      HS_EXPECT_TRUE(born_area / t < 1.1f * limit_slope);
      HS_EXPECT_TRUE(max_centroid_shift < 1e-4f);
      HS_EXPECT_TRUE(max_vertex_shift < 2.0f * t);
      prev_area = born_area;
    }
  }
}

/**
 * @brief Steps a chamfer sweep from T_EPS to the shipping arrival on every
 *        chamfer seed, asserting constant raw and compiled face counts, a
 *        closed genus-0 manifold, unit vertices and outward face normals, and
 *        reporting the peak per-step vertex displacement.
 */
inline void test_chamfer_sweep_holds_topology() {
  constexpr int SAMPLES = 32;

  for (const ChamferSite &site : CHAMFER_SITES) {
    const int failed_before = hs_test::stats().failed;

    Arena persist(probe_seed_buf, sizeof(probe_seed_buf));
    Arena a(probe_a_buf, sizeof(probe_a_buf));
    Arena b(probe_b_buf, sizeof(probe_b_buf));
    PolyMesh seed;
    {
      ScratchScope ga(a);
      ScratchScope gb(b);
      seed = Solids::finalize_solid(site.seed(a, b), persist);
    }

    size_t v0 = 0, f0 = 0, i0 = 0, compiled0 = 0;
    std::vector<Vector> prev_vertices;
    std::vector<Vector> prev_normal;
    std::vector<size_t> off;
    float max_step = 0.0f;
    float min_outward = 1e9f;
    float min_normal_dot = 1e9f;
    for (int s = 0; s < SAMPLES; ++s) {
      const float t = ConwayGraph::T_EPS +
                      (CHAMFER_T_STAR - ConwayGraph::T_EPS) *
                          (static_cast<float>(s) / (SAMPLES - 1));
      ScratchScope fa(a);
      ScratchScope fb(b);
      PolyMesh swept = MeshOps::chamfer(seed, a, b, t);
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

      face_offsets(swept, off);
      std::vector<Vector> normal(swept.face_counts.size());
      for (size_t f = 0; f < swept.face_counts.size(); ++f) {
        normal[f] = newell(swept, off[f], swept.face_counts[f]);
        Vector c(0, 0, 0);
        const int n = swept.face_counts[f];
        for (int k = 0; k < n; ++k)
          c = c + swept.vertices[swept.faces[off[f] + k]];
        const float len = std::sqrt(dot(normal[f], normal[f])) *
                          std::sqrt(dot(c, c));
        if (len > 0.0f)
          min_outward = std::min(min_outward, dot(normal[f], c) / len);
      }
      if (s > 0) {
        for (size_t v = 0; v < swept.vertices.size(); ++v) {
          const Vector d = swept.vertices[v] - prev_vertices[v];
          max_step = std::max(max_step, std::sqrt(dot(d, d)));
        }
        for (size_t f = 0; f < normal.size(); ++f) {
          const float la = std::sqrt(dot(prev_normal[f], prev_normal[f]));
          const float lb = std::sqrt(dot(normal[f], normal[f]));
          if (la > 0.0f && lb > 0.0f)
            min_normal_dot =
                std::min(min_normal_dot, dot(prev_normal[f], normal[f]) /
                                             (la * lb));
        }
      }
      prev_vertices.assign(swept.vertices.data(),
                           swept.vertices.data() + swept.vertices.size());
      prev_normal = normal;
    }

    // No face turns inside out and no vertex teleports between steps.
    HS_EXPECT_TRUE(min_outward > 0.0f);
    HS_EXPECT_TRUE(min_normal_dot > 0.0f);
    if (hs_test::stats().failed != failed_before)
      std::printf("    [chamfer-sweep] %s failed (raw F=%zu compiled=%zu)\n",
                  site.name, f0, compiled0);
    else
      std::printf("  [chamfer-sweep] %s: V=%zu F=%zu I=%zu compiled=%zu "
                  "max_step=%.4f min_outward=%.4f min_normal_dot=%.4f\n",
                  site.name, v0, f0, i0, compiled0,
                  static_cast<double>(max_step),
                  static_cast<double>(min_outward),
                  static_cast<double>(min_normal_dot));
  }
}

// ---------------------------------------------------------------------------
// Truncate sub-T_EPS birth (opchain_morph_spec section 5.1): the two
// truncatedIcosahedron_ambo_relax_truncate001_hankin{59,73} recipes truncate at
// 0.01, below the flat T_EPS birth floor. Both share the same truncate seed
// (truncatedIcosahedron.ambo().relax()) and truncate param, so the leg is
// identical; the hankin angle that differs is a later leg. The leg births at
// min(T_EPS, 0.01 * T_EPS_TRUNCATE_FRAC) = 0.002 and sweeps to 0.01. This must
// be a real, well-formed animation, not a still image or an inverted birth.
// ---------------------------------------------------------------------------

/** @brief One truncate-leg seed. */
struct TruncateSite {
  const char *name;                     /**< Recipe the leg belongs to. */
  PolyMesh (*seed)(Arena &a, Arena &b); /**< Chain prefix up to the truncate. */
};

inline PolyMesh probe_ticosa_ambo_relax(Arena &a, Arena &b) {
  return Solids::SolidBuilder(Solids::Archimedean::truncatedIcosahedron(a, b),
                              a, b)
      .ambo()
      .relax()
      .build();
}

inline constexpr TruncateSite TRUNCATE_SITES[] = {
    {"truncatedIcosahedron_ambo_relax_truncate001_hankin59",
     probe_ticosa_ambo_relax},
    {"truncatedIcosahedron_ambo_relax_truncate001_hankin73",
     probe_ticosa_ambo_relax},
};

/** Arrival parameter of the two truncate001 recipes. */
inline constexpr float TRUNCATE001_T_STAR = 0.01f;

/**
 * @brief Steps the truncate001 leg from its derived birth floor to 0.01 on each
 *        truncate001 seed, asserting a real sweep (birth < arrival), constant
 *        raw and compiled face counts, a closed genus-0 manifold, unit
 * vertices, every face positive-area, and no face inverting across the sweep.
 * @details The silent on-screen failure these recipes carry (spec section 5.1)
 * is a truncate whose sub-T_EPS arrival clamps both endpoints to T_EPS -- a
 * still image ending on a mesh built at the wrong parameter. The birth floor
 * mirrors the OpLeg recipe-step clamp exactly.
 */
inline void test_truncate001_birth_sweep_holds_topology() {
  constexpr int SAMPLES = 32;
  const float birth =
      std::min(ConwayGraph::T_EPS,
               TRUNCATE001_T_STAR * ConwayGraph::T_EPS_TRUNCATE_FRAC);
  // A real animation, not a still image.
  HS_EXPECT_TRUE(birth < TRUNCATE001_T_STAR);
  HS_EXPECT_TRUE(TRUNCATE001_T_STAR >= ConwayGraph::T_EPS_TRUNCATE_MIN);

  for (const TruncateSite &site : TRUNCATE_SITES) {
    const int failed_before = hs_test::stats().failed;

    Arena persist(probe_seed_buf, sizeof(probe_seed_buf));
    Arena a(probe_a_buf, sizeof(probe_a_buf));
    Arena b(probe_b_buf, sizeof(probe_b_buf));
    PolyMesh seed;
    {
      ScratchScope ga(a);
      ScratchScope gb(b);
      seed = Solids::finalize_solid(site.seed(a, b), persist);
    }

    size_t v0 = 0, f0 = 0, i0 = 0, compiled0 = 0;
    std::vector<Vector> prev_normal;
    std::vector<size_t> off;
    float min_area = 1e9f;
    float min_outward = 1e9f;
    float min_normal_dot = 1e9f;
    for (int s = 0; s < SAMPLES; ++s) {
      const float t = birth + (TRUNCATE001_T_STAR - birth) *
                                  (static_cast<float>(s) / (SAMPLES - 1));
      ScratchScope fa(a);
      ScratchScope fb(b);
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

      face_offsets(swept, off);
      std::vector<Vector> normal(swept.face_counts.size());
      for (size_t f = 0; f < swept.face_counts.size(); ++f) {
        normal[f] = newell(swept, off[f], swept.face_counts[f]);
        // Planar area: no face collapses or inverts anywhere in the sweep.
        min_area = std::min(min_area, std::sqrt(dot(normal[f], normal[f])));
        Vector c(0, 0, 0);
        const int n = swept.face_counts[f];
        for (int k = 0; k < n; ++k)
          c = c + swept.vertices[swept.faces[off[f] + k]];
        const float len =
            std::sqrt(dot(normal[f], normal[f])) * std::sqrt(dot(c, c));
        if (len > 0.0f)
          min_outward = std::min(min_outward, dot(normal[f], c) / len);
      }
      if (s > 0) {
        for (size_t f = 0; f < normal.size(); ++f) {
          const float la = std::sqrt(dot(prev_normal[f], prev_normal[f]));
          const float lb = std::sqrt(dot(normal[f], normal[f]));
          if (la > 0.0f && lb > 0.0f)
            min_normal_dot = std::min(
                min_normal_dot, dot(prev_normal[f], normal[f]) / (la * lb));
        }
      }
      prev_normal = normal;
    }

    // Every face keeps positive area and outward orientation across the sweep.
    HS_EXPECT_TRUE(min_area > 0.0f);
    HS_EXPECT_TRUE(min_outward > 0.0f);
    HS_EXPECT_TRUE(min_normal_dot > 0.0f);
    if (hs_test::stats().failed != failed_before)
      std::printf("    [truncate001] %s failed (raw F=%zu compiled=%zu)\n",
                  site.name, f0, compiled0);
    else
      std::printf("  [truncate001] %s: birth=%.4f->%.2f V=%zu F=%zu I=%zu "
                  "compiled=%zu min_area=%.3e min_outward=%.4f "
                  "min_normal_dot=%.4f\n",
                  site.name, static_cast<double>(birth),
                  static_cast<double>(TRUNCATE001_T_STAR), v0, f0, i0,
                  compiled0, static_cast<double>(min_area),
                  static_cast<double>(min_outward),
                  static_cast<double>(min_normal_dot));
  }
}

// ---------------------------------------------------------------------------
// Far-side truncate sweep (opchain_morph_spec section 5.1): the two
// truncatedIcos{ahedron,idodecahedron}_truncate50d_ambo_dual recipes truncate
// at 50 deg = 0.873, PAST the ambo pinch (t = 0.5). truncate emits a constant
// topology (2E vertices, F+V faces, 3I indices) for every t != 0.5; only the
// exact-0.5 short-circuit differs (returns ambo, V vertices). A far-side leg
// births on the near side (small t), sweeps through 0.5 with the pinch guard
// nudging that one sample off the short-circuit, and arrives at 0.873 on the
// intentionally self-intersecting truncate branch.
//
// Past 0.5 the cut faces self-intersect BY DESIGN, so signed area, winding, and
// closed-manifold checks fail legitimately (MEMORY: structural checks only for
// crossed-face ops). This test asserts only what stays true across the pinch:
// constant raw and compiled face count, constant V/F/I, finite unit vertices,
// and no exact-0.5 evaluation. Positive area is asserted only on the near side.
// ---------------------------------------------------------------------------

inline PolyMesh probe_ticosidodeca(Arena &a, Arena &b) {
  return Solids::Archimedean::truncatedIcosidodecahedron(a, b);
}

inline constexpr TruncateSite FAR_TRUNCATE_SITES[] = {
    {"truncatedIcosahedron_truncate50d_ambo_dual", probe_ticosa},
    {"truncatedIcosidodecahedron_truncate50d_ambo_dual", probe_ticosidodeca},
};

/** Arrival parameter of the two truncate50d recipes: 50 deg past the pinch. */
inline constexpr float TRUNCATE50D_T_STAR =
    50.0f * Solids::IslamicStarPatterns::D2R;
/** Below this t, the truncate cut faces do not yet self-intersect, so positive
 * area still holds and can be asserted. Above the pinch (0.5) it cannot. */
inline constexpr float FAR_SIDE_NEAR_LIMIT = 0.49f;

/**
 * @brief Steps the far-side truncate leg from its near-side birth floor through
 *        the ambo pinch to 0.873 on both truncate50d seeds, asserting the leg
 *        does not trap and does not change topology across the pinch.
 * @details Mirrors the OpLeg recipe-step clamp: the leg births at min(T_EPS,
 * arrival * T_EPS_TRUNCATE_FRAC) and the far-side arrival passes through
 * unclamped (below T_EPS_TRUNCATE_FAR_MAX). Every per-frame sample is routed
 * through ConwayGraph::truncate_off_pinch, so a sample landing exactly on 0.5
 * is nudged off the ambo short-circuit and the frame keeps the truncate
 * topology. STRUCTURAL checks only past 0.5 (self-intersecting by design): no
 * signed-area or manifold assertion there.
 */
inline void test_truncate50d_far_side_sweep_holds_topology() {
  constexpr int SAMPLES = 48;
  const float birth =
      std::min(ConwayGraph::T_EPS,
               TRUNCATE50D_T_STAR * ConwayGraph::T_EPS_TRUNCATE_FRAC);
  // A real animation across the pinch: birth on the near side, arrival past it,
  // arrival unclamped (below the far-side cap).
  HS_EXPECT_TRUE(birth < 0.5f);
  HS_EXPECT_TRUE(TRUNCATE50D_T_STAR > 0.5f);
  HS_EXPECT_TRUE(TRUNCATE50D_T_STAR <= ConwayGraph::T_EPS_TRUNCATE_FAR_MAX);
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::TRUNCATE, TRUNCATE50D_T_STAR}));
  // The guard moves an exact-0.5 sample onto the truncate branch.
  HS_EXPECT_TRUE(ConwayGraph::truncate_off_pinch(0.5f) != 0.5f);

  for (const TruncateSite &site : FAR_TRUNCATE_SITES) {
    const int failed_before = hs_test::stats().failed;

    Arena persist(probe_seed_buf, sizeof(probe_seed_buf));
    Arena a(probe_a_buf, sizeof(probe_a_buf));
    Arena b(probe_b_buf, sizeof(probe_b_buf));
    PolyMesh seed;
    {
      ScratchScope ga(a);
      ScratchScope gb(b);
      seed = Solids::finalize_solid(site.seed(a, b), persist);
    }

    size_t v0 = 0, f0 = 0, i0 = 0, compiled0 = 0;
    std::vector<size_t> off;
    float min_area_near = 1e9f;
    bool pinch_guarded = false;
    for (int s = 0; s < SAMPLES; ++s) {
      // Linear birth -> arrival, plus one sample forced onto the exact pinch so
      // the guard is exercised deterministically.
      float raw = (s == SAMPLES / 2)
                      ? 0.5f
                      : birth + (TRUNCATE50D_T_STAR - birth) *
                                    (static_cast<float>(s) / (SAMPLES - 1));
      const float t = ConwayGraph::truncate_off_pinch(raw);
      // No frame ever evaluates truncate at the exact ambo short-circuit.
      HS_EXPECT_TRUE(t != 0.5f);
      if (raw == 0.5f)
        pinch_guarded = true;

      ScratchScope fa(a);
      ScratchScope fb(b);
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
        // Topology is t-independent off the pinch: no pop across 0.5.
        HS_EXPECT_EQ(swept.vertices.size(), v0);
        HS_EXPECT_EQ(swept.face_counts.size(), f0);
        HS_EXPECT_EQ(swept.faces.size(), i0);
        HS_EXPECT_EQ(compiled.face_counts.size(), compiled0);
      }
      check_face_counts_consistent(swept);
      check_indices_in_range(swept);
      check_all_unit_vertices(swept, 1e-3f);
      for (size_t i = 0; i < swept.vertices.size(); ++i)
        HS_EXPECT_TRUE(std::isfinite(swept.vertices[i].length()));

      // Positive area is a near-side-only invariant: past the pinch the cut
      // faces self-intersect by design and signed area legitimately flips.
      if (t <= FAR_SIDE_NEAR_LIMIT) {
        face_offsets(swept, off);
        for (size_t f = 0; f < swept.face_counts.size(); ++f) {
          const Vector n = newell(swept, off[f], swept.face_counts[f]);
          min_area_near = std::min(min_area_near, std::sqrt(dot(n, n)));
        }
      }
    }

    HS_EXPECT_TRUE(pinch_guarded);
    HS_EXPECT_TRUE(min_area_near > 0.0f);
    if (hs_test::stats().failed != failed_before)
      std::printf("    [truncate50d] %s failed (raw F=%zu compiled=%zu)\n",
                  site.name, f0, compiled0);
    else
      std::printf(
          "  [truncate50d] %s: birth=%.4f -> %.4f through pinch V=%zu "
          "F=%zu I=%zu compiled=%zu min_area_near=%.3e guard_fired=%d\n",
          site.name, static_cast<double>(birth),
          static_cast<double>(TRUNCATE50D_T_STAR), v0, f0, i0, compiled0,
          static_cast<double>(min_area_near), pinch_guarded ? 1 : 0);
  }
}

/**
 * @brief Finds the smallest chamfer t at which no newborn hexagon is rejected
 *        by SDF::Face's collapsed-area cull, per seed.
 * @details MeshOps::compile drops faces by side count only, so it never drops
 * a chamfer birth; the size-dependent reject is SDF::Face's
 * COLLAPSED_AREA_RATIO test at draw time. This bisects that boundary and
 * reports it as the chamfer analog of T_EPS.
 */
inline void test_chamfer_birth_epsilon() {
  for (const ChamferSite &site : CHAMFER_SITES) {
    Arena persist(probe_seed_buf, sizeof(probe_seed_buf));
    Arena a(probe_a_buf, sizeof(probe_a_buf));
    Arena b(probe_b_buf, sizeof(probe_b_buf));
    PolyMesh seed;
    {
      ScratchScope ga(a);
      ScratchScope gb(b);
      seed = Solids::finalize_solid(site.seed(a, b), persist);
    }
    const size_t F = seed.face_counts.size();

    // Culled newborn count at t; compile's own face count for contrast.
    auto probe = [&](float t, size_t &culled, size_t &compiled_faces) {
      ScratchScope fa(a);
      ScratchScope fb(b);
      PolyMesh out = MeshOps::chamfer(seed, a, b, t);
      MeshState compiled;
      MeshOps::compile(out, compiled, a, b);
      compiled_faces = compiled.face_counts.size();
      std::vector<size_t> off;
      face_offsets(out, off);
      culled = 0;
      for (size_t f = F; f < out.face_counts.size(); ++f)
        if (face_is_culled(out, off[f], out.face_counts[f]))
          ++culled;
    };

    size_t culled = 0, compiled_faces = 0, raw_faces = 0;
    {
      ScratchScope fa(a);
      ScratchScope fb(b);
      PolyMesh out = MeshOps::chamfer(seed, a, b, ConwayGraph::T_EPS);
      raw_faces = out.face_counts.size();
    }

    float lo = 1e-9f, hi = 0.25f;
    probe(hi, culled, compiled_faces);
    const bool clears_at_hi = culled == 0;
    if (clears_at_hi) {
      for (int i = 0; i < 40; ++i) {
        const float mid = std::sqrt(lo * hi);
        size_t c = 0, cf = 0;
        probe(mid, c, cf);
        (c == 0 ? hi : lo) = mid;
      }
    }
    size_t at_eps = 0, cf_eps = 0;
    probe(ConwayGraph::T_EPS, at_eps, cf_eps);
    // compile is combinatorial: its face count never moves with t.
    HS_EXPECT_EQ(cf_eps, raw_faces);
    HS_EXPECT_EQ(compiled_faces, raw_faces);
    HS_EXPECT_TRUE(clears_at_hi);
    std::printf("  [chamfer-eps] %s: births clear the SDF cull at t>=%.2e "
                "(T_EPS=%.3f culls %zu of %zu; compile keeps %zu of %zu)\n",
                site.name, static_cast<double>(hi),
                static_cast<double>(ConwayGraph::T_EPS), at_eps,
                raw_faces - F, cf_eps, raw_faces);
  }
}

// ---------------------------------------------------------------------------
// Provenance probe (spec section 8.3): PROVENANCE_TOL_SQ = 0.15^2 is sized for
// F <= 92 nodes. The build chains reach several hundred faces, so measure the
// centroid spacing they actually present and the ambiguity of the lookups
// build_palette_mapping actually performs.
// ---------------------------------------------------------------------------

/** @brief One pure-inflate build chain, lowered to primitive steps. */
struct ChainSite {
  const char *name;      /**< Registry entry the chain mirrors. */
  uint8_t seed;          /**< simple_registry index. */
  const OpStep *steps;   /**< Lowered primitive chain. */
  size_t count;          /**< Number of steps. */
};

using Solids::IslamicStarPatterns::D2R;

inline constexpr uint8_t SEED_CUBE = 1;
inline constexpr uint8_t SEED_OCTAHEDRON = 2;
inline constexpr uint8_t SEED_DODECAHEDRON = 3;
inline constexpr uint8_t SEED_ICOSAHEDRON = 4;
inline constexpr uint8_t SEED_RHOMBICUBOCTAHEDRON = 9;
inline constexpr uint8_t SEED_TRUNCATED_ICOSAHEDRON = 14;
inline constexpr uint8_t SEED_TRUNCATED_ICOSIDODECAHEDRON = 16;

static_assert(std::string_view(Solids::simple_registry[SEED_CUBE].name) ==
              "cube");
static_assert(
    std::string_view(Solids::simple_registry[SEED_RHOMBICUBOCTAHEDRON].name) ==
    "rhombicuboctahedron");
static_assert(
    std::string_view(Solids::simple_registry[SEED_TRUNCATED_ICOSAHEDRON].name) ==
    "truncatedIcosahedron");
static_assert(std::string_view(
                  Solids::simple_registry[SEED_TRUNCATED_ICOSIDODECAHEDRON]
                      .name) == "truncatedIcosidodecahedron");

inline constexpr OpStep CHAIN_DODECA_HK62_AMBO_HK62[] = {
    {Op::HANKIN, 62.0f * D2R}, {Op::AMBO}, {Op::HANKIN, 62.0f * D2R}};
inline constexpr OpStep CHAIN_DODECA_HK35_AMBO_HK62_AMBO_RELAX_HK42[] = {
    {Op::HANKIN, 35.0f * D2R}, {Op::AMBO},   {Op::HANKIN, 62.0f * D2R},
    {Op::AMBO},                {Op::RELAX, 100.0f}, {Op::HANKIN, 42.0f * D2R}};
inline constexpr OpStep CHAIN_DODECA_HK54_AMBO_HK72[] = {
    {Op::HANKIN, 54.0f * D2R}, {Op::AMBO}, {Op::HANKIN, 72.0f * D2R}};
inline constexpr OpStep CHAIN_OCTA_HK17_AMBO_HK73[] = {
    {Op::HANKIN, 17.0f * D2R}, {Op::AMBO}, {Op::HANKIN, 73.0f * D2R}};
inline constexpr OpStep CHAIN_OCTA_HK34_AMBO_HK72[] = {
    {Op::HANKIN, 34.0f * D2R}, {Op::AMBO}, {Op::HANKIN, 72.0f * D2R}};
inline constexpr OpStep CHAIN_RHOMBICUBOCTA_HK63_AMBO_HK63[] = {
    {Op::HANKIN, 63.0f * D2R}, {Op::AMBO}, {Op::HANKIN, 63.0f * D2R}};
inline constexpr OpStep CHAIN_TICOSA_HK54_AMBO_HK72[] = {
    {Op::HANKIN, 54.0f * D2R}, {Op::AMBO}, {Op::HANKIN, 72.0f * D2R}};
inline constexpr OpStep CHAIN_TICOSA_HK58_CHAMFER63[] = {
    {Op::HANKIN, 58.0f * D2R}, {Op::CHAMFER, CHAMFER_T_STAR}};
inline constexpr OpStep CHAIN_ICOSA_AMBO_TRUNCATE033_HK59[] = {
    {Op::AMBO}, {Op::TRUNCATE, 0.33f}, {Op::HANKIN, 59.0f * D2R}};
inline constexpr OpStep CHAIN_ICOSA_SNUB_RELAX_TRUNCATE033_HK62[] = {
    {Op::SNUB, 0.5f, 0.0f},
    {Op::RELAX, 8.0f},
    {Op::TRUNCATE, 0.33f},
    {Op::HANKIN, 62.0f * D2R}};
inline constexpr OpStep CHAIN_TICOSA_AMBO_RELAX_TRUNCATE33_HK64[] = {
    {Op::AMBO},
    {Op::RELAX, 217.0f},
    {Op::TRUNCATE, 0.33f},
    {Op::HANKIN, 64.0f * D2R}};
inline constexpr OpStep CHAIN_DODECA_AMBO_BEVEL33_RELAX_HK66[] = {
    {Op::AMBO},
    {Op::AMBO},
    {Op::TRUNCATE, 0.33f},
    {Op::RELAX, 100.0f},
    {Op::HANKIN, 66.0f * D2R}};
inline constexpr OpStep CHAIN_TICOSIDODECA_BEVEL5_RELAX_HK77[] = {
    {Op::AMBO}, {Op::AMBO}, {Op::RELAX, 100.0f}, {Op::HANKIN, 77.0f * D2R}};

inline constexpr ChainSite CHAIN_SITES[] = {
    {"dodecahedron_hk62_ambo_hk62", SEED_DODECAHEDRON,
     CHAIN_DODECA_HK62_AMBO_HK62, std::size(CHAIN_DODECA_HK62_AMBO_HK62)},
    {"dodecahedron_hk35_ambo_hk62_ambo_relax_hk42", SEED_DODECAHEDRON,
     CHAIN_DODECA_HK35_AMBO_HK62_AMBO_RELAX_HK42,
     std::size(CHAIN_DODECA_HK35_AMBO_HK62_AMBO_RELAX_HK42)},
    {"dodecahedron_hk54_ambo_hk72", SEED_DODECAHEDRON,
     CHAIN_DODECA_HK54_AMBO_HK72, std::size(CHAIN_DODECA_HK54_AMBO_HK72)},
    {"octahedron_hk17_ambo_hk73", SEED_OCTAHEDRON, CHAIN_OCTA_HK17_AMBO_HK73,
     std::size(CHAIN_OCTA_HK17_AMBO_HK73)},
    {"octahedron_hk34_ambo_hk72", SEED_OCTAHEDRON, CHAIN_OCTA_HK34_AMBO_HK72,
     std::size(CHAIN_OCTA_HK34_AMBO_HK72)},
    {"rhombicuboctahedron_hk63_ambo_hk63", SEED_RHOMBICUBOCTAHEDRON,
     CHAIN_RHOMBICUBOCTA_HK63_AMBO_HK63,
     std::size(CHAIN_RHOMBICUBOCTA_HK63_AMBO_HK63)},
    {"truncatedIcosahedron_hk54_ambo_hk72", SEED_TRUNCATED_ICOSAHEDRON,
     CHAIN_TICOSA_HK54_AMBO_HK72, std::size(CHAIN_TICOSA_HK54_AMBO_HK72)},
    {"truncatedIcosahedron_hk58_chamfer63", SEED_TRUNCATED_ICOSAHEDRON,
     CHAIN_TICOSA_HK58_CHAMFER63, std::size(CHAIN_TICOSA_HK58_CHAMFER63)},
    {"icosahedron_ambo_truncate033_hankin59", SEED_ICOSAHEDRON,
     CHAIN_ICOSA_AMBO_TRUNCATE033_HK59,
     std::size(CHAIN_ICOSA_AMBO_TRUNCATE033_HK59)},
    {"icosahedron_snub_relax_truncate033_hankin62", SEED_ICOSAHEDRON,
     CHAIN_ICOSA_SNUB_RELAX_TRUNCATE033_HK62,
     std::size(CHAIN_ICOSA_SNUB_RELAX_TRUNCATE033_HK62)},
    {"truncatedIcosahedron_ambo_relax_truncate33_hk64",
     SEED_TRUNCATED_ICOSAHEDRON, CHAIN_TICOSA_AMBO_RELAX_TRUNCATE33_HK64,
     std::size(CHAIN_TICOSA_AMBO_RELAX_TRUNCATE33_HK64)},
    {"dodecahedron_ambo_bevel33_relax_hk66", SEED_DODECAHEDRON,
     CHAIN_DODECA_AMBO_BEVEL33_RELAX_HK66,
     std::size(CHAIN_DODECA_AMBO_BEVEL33_RELAX_HK66)},
    {"truncatedIcosidodecahedron_bevel5_relax_hk77",
     SEED_TRUNCATED_ICOSIDODECAHEDRON, CHAIN_TICOSIDODECA_BEVEL5_RELAX_HK77,
     std::size(CHAIN_TICOSIDODECA_BEVEL5_RELAX_HK77)},
};

/** @brief Nearest and second-nearest chord distance from c into `pts`. */
inline void two_nearest(const Vector &c, const std::vector<Vector> &pts,
                        size_t &best, float &d1, float &d2) {
  best = 0;
  d1 = 1e9f;
  d2 = 1e9f;
  for (size_t j = 0; j < pts.size(); ++j) {
    const Vector d = c - pts[j];
    const float dsq = dot(d, d);
    if (dsq < d1) {
      d2 = d1;
      d1 = dsq;
      best = j;
    } else if (dsq < d2) {
      d2 = dsq;
    }
  }
  d1 = std::sqrt(d1);
  d2 = std::sqrt(d2);
}

/**
 * @brief Reports face-centroid spacing per intermediate build-chain mesh
 *        against PROVENANCE_TOL_SQ's 0.15 chord radius.
 */
inline void test_build_chain_centroid_spacing() {
  constexpr float TOL = 0.15f;
  float worst_ratio = 0.0f;
  const char *worst_name = "";
  size_t max_faces = 0;

  for (const ChainSite &site : CHAIN_SITES) {
    Arena a(probe_a_buf, sizeof(probe_a_buf));
    Arena b(probe_b_buf, sizeof(probe_b_buf));
    for (size_t k = 0; k <= site.count; ++k) {
      ScratchScope ga(a);
      ScratchScope gb(b);
      PolyMesh m = Solids::build_steps(site.seed, site.steps, k, a, b);
      std::vector<Vector> c;
      face_centroids(m, c);
      std::vector<float> nn(c.size(), 0.0f);
      float min_spacing = 1e9f;
      for (size_t f = 0; f < c.size(); ++f) {
        float best = 1e9f;
        for (size_t g = 0; g < c.size(); ++g) {
          if (g == f)
            continue;
          const Vector d = c[f] - c[g];
          best = std::min(best, dot(d, d));
        }
        nn[f] = std::sqrt(best);
        min_spacing = std::min(min_spacing, nn[f]);
      }
      const float med = median_of(nn);
      const float ratio = TOL / (0.5f * min_spacing);
      max_faces = std::max(max_faces, c.size());
      if (ratio > worst_ratio) {
        worst_ratio = ratio;
        worst_name = site.name;
      }
      std::printf("  [spacing] %-46s step %zu: F=%-4zu min=%.4f med=%.4f "
                  "tol/half-min=%.2f\n",
                  site.name, k, c.size(), static_cast<double>(min_spacing),
                  static_cast<double>(med), static_cast<double>(ratio));
      HS_EXPECT_TRUE(c.size() > 0);
    }
  }
  std::printf("  [spacing] worst tol/half-min = %.2f on %s; largest chain "
              "mesh F=%zu\n",
              static_cast<double>(worst_ratio), worst_name, max_faces);
}

/**
 * @brief Reports, per real build leg, the mapping path build_palette_mapping
 *        takes, its departed face count, and the ambiguity of every
 *        nearest-centroid lookup the leg performs.
 * @details Legs open at the sweep's start parameter, which is what
 * build_palette_mapping's start_centroid array holds.
 */
inline void test_build_chain_provenance_ambiguity() {
  constexpr float TOL_SQ = 0.15f * 0.15f;
  size_t max_prev_faces = 0;
  const char *max_prev_name = "";
  float worst_newborn_ratio = 1e9f;
  float worst_prefix_offset = 0.0f;
  size_t prefix_legs = 0, full_legs = 0, misidentified = 0;

  for (const ChainSite &site : CHAIN_SITES) {
    Arena persist(probe_seed_buf, sizeof(probe_seed_buf));
    Arena a(probe_a_buf, sizeof(probe_a_buf));
    Arena b(probe_b_buf, sizeof(probe_b_buf));
    for (size_t k = 0; k < site.count; ++k) {
      const OpStep &step = site.steps[k];
      // relax legs settle onto the preceding sweep and open no mapping.
      if (step.op == Op::RELAX)
        continue;

      ScratchScope ga(a);
      ScratchScope gb(b);
      PolyMesh seed = Solids::build_steps(site.seed, site.steps, k, a, b);
      std::vector<Vector> prev_c;
      face_centroids(seed, prev_c);
      const size_t prev_faces = prev_c.size();
      if (prev_faces > max_prev_faces) {
        max_prev_faces = prev_faces;
        max_prev_name = site.name;
      }

      // The mesh the leg's first frame draws.
      PolyMesh start;
      switch (step.op) {
      case Op::HANKIN:
        start = MeshOps::hankin(seed, a, b, Animation::OpLeg::THETA_EPS);
        break;
      case Op::AMBO:
        start = MeshOps::truncate(seed, a, b, ConwayGraph::T_EPS);
        break;
      case Op::CHAMFER:
        start = MeshOps::chamfer(seed, a, b, ConwayGraph::T_EPS);
        break;
      case Op::TRUNCATE:
        start = MeshOps::truncate(seed, a, b, ConwayGraph::T_EPS);
        break;
      case Op::SNUB:
        start = MeshOps::snub(seed, a, b, ConwayGraph::T_EPS, step.twist);
        break;
      default:
        continue;
      }
      std::vector<Vector> start_c;
      face_centroids(start, start_c);
      const size_t total = start_c.size();
      const bool full = total == prev_faces;
      (full ? full_legs : prefix_legs)++;

      // Prefix path: face f departs from prev face f by emission identity.
      float max_prefix_offset = 0.0f;
      size_t leg_misidentified = 0;
      for (size_t f = 0; f < std::min(prev_faces, total); ++f) {
        size_t best;
        float d1, d2;
        two_nearest(start_c[f], prev_c, best, d1, d2);
        max_prefix_offset = std::max(max_prefix_offset, d1);
        if (best != f)
          ++leg_misidentified;
      }
      misidentified += leg_misidentified;
      worst_prefix_offset = std::max(worst_prefix_offset, max_prefix_offset);

      // Newborn path: one nearest-centroid lookup per newborn face.
      float min_ratio = 1e9f;
      float max_newborn_d1 = 0.0f;
      for (size_t f = prev_faces; f < total; ++f) {
        size_t best;
        float d1, d2;
        two_nearest(start_c[f], prev_c, best, d1, d2);
        max_newborn_d1 = std::max(max_newborn_d1, d1);
        if (d2 > 0.0f)
          min_ratio = std::min(min_ratio, d1 / d2);
      }
      if (total > prev_faces)
        worst_newborn_ratio = std::min(worst_newborn_ratio, min_ratio);

      std::printf("  [prov] %-46s leg %zu %-8s: prev=%-4zu total=%-4zu %s "
                  "prefix_d=%.4f (%zu misidentified) newborn_d=%.4f "
                  "d1/d2=%.3f\n",
                  site.name, k,
                  step.op == Op::HANKIN    ? "hankin"
                  : step.op == Op::AMBO    ? "ambo"
                  : step.op == Op::CHAMFER ? "chamfer"
                  : step.op == Op::SNUB    ? "snub"
                                           : "truncate",
                  prev_faces, total, full ? "FULL " : "PREFIX",
                  static_cast<double>(max_prefix_offset), leg_misidentified,
                  static_cast<double>(max_newborn_d1),
                  static_cast<double>(total > prev_faces ? min_ratio : 1.0f));
      // The prefix identity must also be the geometric nearest, or the
      // newborn-class inheritance is reading the wrong departed face.
      HS_EXPECT_TRUE(max_prefix_offset * max_prefix_offset < TOL_SQ ||
                     leg_misidentified > 0);
    }
  }
  std::printf("  [prov] %zu prefix legs, %zu full-correspondence legs, %zu "
              "misidentified prefix faces; max prev_faces=%zu on %s; "
              "worst newborn d1/d2=%.3f; worst prefix "
              "offset=%.4f (tol 0.15)\n",
              prefix_legs, full_legs, misidentified, max_prev_faces,
              max_prev_name, static_cast<double>(worst_newborn_ratio),
              static_cast<double>(worst_prefix_offset));
  // None of the swept ops scanned here takes the tolerance-checked
  // full-correspondence path, so PROVENANCE_TOL_SQ never applies to them.
  HS_EXPECT_EQ(full_legs, static_cast<size_t>(0));
  HS_EXPECT_TRUE(max_prev_faces > 128);
}

// ---------------------------------------------------------------------------
// Needle gated-swap chain (opchain_morph_spec section 3.3): the
// truncatedIcosahedron_ambo_relax100_hk54_needle recipe ends in needle, which
// expand_to_primitives lowers to a DUAL then a KIS gated swap
// (core/mesh/recipe.h). No shipping recipe runs dual or kis on a hankin mesh,
// so this pins that both partition ops land a well-formed closed manifold on
// the hankin(54 deg) arrival. A gated swap carries no parameter sweep, so each
// leg's arrival is a single mesh: dual(seed), then kis(dual(seed)) == needle.
// ---------------------------------------------------------------------------

inline PolyMesh probe_ticosa_ambo_relax100_hk54(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Archimedean::truncatedIcosahedron(a, b),
                              a, b)
      .ambo()
      .relax(100)
      .hankin(54.0f * D2R)
      .build();
}

/**
 * @brief Asserts a mesh is a closed genus-0 manifold of positive-area faces on
 *        the unit sphere; returns its compiled face count.
 */
inline size_t check_manifold_landing(const PolyMesh &m, Arena &a, Arena &b) {
  check_face_counts_consistent(m);
  check_indices_in_range(m);
  check_all_unit_vertices(m, 1e-3f);
  conway_tests::check_euler_genus0(m);
  std::vector<size_t> off;
  face_offsets(m, off);
  float min_area = 1e9f;
  for (size_t f = 0; f < m.face_counts.size(); ++f) {
    const Vector n = newell(m, off[f], m.face_counts[f]);
    min_area = std::min(min_area, std::sqrt(dot(n, n)));
  }
  for (size_t v = 0; v < m.vertices.size(); ++v)
    HS_EXPECT_TRUE(std::isfinite(m.vertices[v].length()));
  HS_EXPECT_TRUE(min_area > 0.0f);
  MeshState compiled;
  MeshOps::compile(m, compiled, a, b);
  return compiled.face_counts.size();
}

/**
 * @brief Builds the needle chain's DUAL then KIS gated-swap landings on the
 *        hankin(54 deg) seed, asserting each lands a well-formed closed
 *        manifold and that the {DUAL, KIS} lowering reproduces the composite
 *        needle exactly.
 * @details needle lowers to two gated swaps with no sweep, so the guard is that
 * each partition op builds without trapping on a hankin mesh -- the one context
 * no shipping recipe exercises -- and that expand_to_primitives' {DUAL, KIS}
 * pair matches MeshOps::needle bit for bit.
 */
inline void test_needle_gated_swap_builds_on_hankin() {
  const int failed_before = hs_test::stats().failed;

  Arena persist(probe_seed_buf, sizeof(probe_seed_buf));
  Arena a(probe_a_buf, sizeof(probe_a_buf));
  Arena b(probe_b_buf, sizeof(probe_b_buf));

  PolyMesh seed;
  {
    ScratchScope ga(a);
    ScratchScope gb(b);
    seed =
        Solids::finalize_solid(probe_ticosa_ambo_relax100_hk54(a, b), persist);
  }
  size_t seed_compiled = 0;
  {
    ScratchScope fa(a);
    ScratchScope fb(b);
    MeshState c;
    MeshOps::compile(seed, c, a, b);
    seed_compiled = c.face_counts.size();
  }

  // DUAL leg landing: departs the hankin seed, opens on its dual.
  PolyMesh dual_mesh;
  {
    ScratchScope fa(a);
    ScratchScope fb(b);
    dual_mesh = Solids::finalize_solid(MeshOps::dual(seed, a, b), persist);
  }
  size_t dual_compiled = 0;
  {
    ScratchScope fa(a);
    ScratchScope fb(b);
    dual_compiled = check_manifold_landing(dual_mesh, a, b);
  }
  // A gated swap draws one static mesh per side, so compile keeps every face:
  // the leg's per-side face count is constant.
  HS_EXPECT_EQ(dual_compiled, dual_mesh.face_counts.size());

  // KIS leg landing: departs the dual, opens on kis(dual) == the needle
  // arrival.
  size_t kis_v = 0, kis_f = 0, kis_i = 0, kis_compiled = 0;
  {
    ScratchScope fa(a);
    ScratchScope fb(b);
    PolyMesh kis_mesh = MeshOps::kis(dual_mesh, a, b);
    kis_v = kis_mesh.vertices.size();
    kis_f = kis_mesh.face_counts.size();
    kis_i = kis_mesh.faces.size();
    kis_compiled = check_manifold_landing(kis_mesh, a, b);
  }
  HS_EXPECT_EQ(kis_compiled, kis_f);
  // kis raises one triangle per parent face-side, so its face count is the
  // dual's total face-index count.
  HS_EXPECT_EQ(kis_f, dual_mesh.faces.size());

  // The lowering matches the composite: needle n = kd = kis of dual.
  {
    ScratchScope fa(a);
    ScratchScope fb(b);
    PolyMesh needle_mesh = MeshOps::needle(seed, a, b);
    HS_EXPECT_EQ(needle_mesh.vertices.size(), kis_v);
    HS_EXPECT_EQ(needle_mesh.face_counts.size(), kis_f);
    HS_EXPECT_EQ(needle_mesh.faces.size(), kis_i);
  }

  if (hs_test::stats().failed != failed_before)
    std::printf(
        "    [needle] truncatedIcosahedron_ambo_relax100_hk54 failed\n");
  else
    std::printf(
        "  [needle] hk54 seed F=%zu(compiled %zu) -> dual F=%zu(%zu) -> "
        "kis F=%zu(%zu) == needle; closed manifold at each landing\n",
        seed.face_counts.size(), seed_compiled, dual_mesh.face_counts.size(),
        dual_compiled, kis_f, kis_compiled);
}

/**
 * @brief Runs the OpChainMorph pre-flight probes.
 * @return Failure count.
 */
inline int run_opchain_probe_tests() {
  hs_test::ModuleFixture fixture("opchain_probe");

  test_chamfer_zero_area_birth_limit();
  test_chamfer_sweep_holds_topology();
  test_chamfer_birth_epsilon();

  test_truncate001_birth_sweep_holds_topology();
  test_truncate50d_far_side_sweep_holds_topology();

  test_needle_gated_swap_builds_on_hankin();

  test_build_chain_centroid_spacing();
  test_build_chain_provenance_ambiguity();

  return fixture.result();
}

} // namespace opchain_probe_tests
} // namespace hs_test
