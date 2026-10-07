/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Geometry and provenance checks for the OpChainMorph build recipes
 * (docs/specs/opchain_morph_spec.md).
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <span>
#include <vector>

#include "core/animation/animation.h"
#include "core/mesh/conway.h"
#include "core/mesh/conway_graph.h"
#include "core/mesh/hankin.h"
#include "core/mesh/recipe.h"
#include "core/mesh/solids.h"
#include "core/render/sdf.h"
#include "tests/conway_test_util.h"
#include "tests/mesh_test_util.h"
#include "tests/test_conway.h" // check_euler_characteristic_two
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace opchain_probe_tests {

inline uint8_t probe_a_buf[768 * 1024];           /**< Op output arena. */
inline uint8_t probe_b_buf[768 * 1024];           /**< Op scratch arena. */
inline uint8_t probe_seed_buf[512 * 1024];        /**< Held seed arena. */
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
inline void face_centroids(const PolyMesh &m, std::vector<math::Vector> &out) {
  out.assign(m.face_counts.size(), math::Vector(0, 0, 0));
  size_t off = 0;
  for (size_t f = 0; f < m.face_counts.size(); ++f) {
    const int n = m.face_counts[f];
    out[f] = face_centroid_unit(m, off, n);
    off += n;
  }
}

/** @brief Whether SDF::Face rejects a face outright (collapsed or off-band). */
inline bool face_is_culled(const PolyMesh &m, size_t off, int n) {
  if (static_cast<size_t>(n) > SDF::FaceScratchBuffer::MAX_VERTS)
    return false;
  SDF::Face face(
      std::span<const math::Vector>(m.vertices.data(), m.vertices.size()),
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
// Chamfer sweep over Hankin meshes to characterize the T_EPS boundary.
// ---------------------------------------------------------------------------

/** @brief One chamfer-leg seed. */
struct ChamferSite {
  const char *name;                     /**< Diagnostic label. */
  PolyMesh (*seed)(Arena &a, Arena &b); /**< Chain prefix up to the chamfer. */
};

inline PolyMesh probe_cube(Arena &a, Arena &b) {
  return Solids::Platonic::cube(a, b);
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

/** Arrival parameter of the shipping chamfer leg. */
inline constexpr float CHAMFER_T_STAR = Solids::CHAMFER_T_MAX;

/**
 * @brief Measures chamfer's zero-area birth limit: newborn hexagon area and
 *        preserved-face displacement against t, on every chamfer seed.
 * @details Chamfer emits F shrunk primaries then E hexagons, so the split is
 * exact in emission order.
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
    std::vector<math::Vector> seed_centroid;
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
        const math::Vector n =
            face_area_vector(out, off[f], out.face_counts[f]);
        born_area += std::sqrt(math::dot(n, n));
      }
      // Preserved faces: same side count, centroid held, and every vertex
      // converging on a seed vertex as t -> 0.
      float max_centroid_shift = 0.0f;
      float max_vertex_shift = 0.0f;
      for (size_t f = 0; f < F; ++f) {
        HS_EXPECT_EQ(static_cast<int>(out.face_counts[f]),
                     static_cast<int>(seed.face_counts[f]));
        math::Vector c(0, 0, 0);
        const int n = out.face_counts[f];
        for (int k = 0; k < n; ++k) {
          const math::Vector &v = out.vertices[out.faces[off[f] + k]];
          c = c + v;
          float near_sq = 1e9f;
          for (size_t s = 0; s < seed.vertices.size(); ++s) {
            const math::Vector e = v - seed.vertices[s];
            near_sq = std::min(near_sq, math::dot(e, e));
          }
          max_vertex_shift = std::max(max_vertex_shift, std::sqrt(near_sq));
        }
        c = c.normalized();
        const math::Vector d = c - seed_centroid[f];
        max_centroid_shift =
            std::max(max_centroid_shift, std::sqrt(math::dot(d, d)));
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
      else
        HS_EXPECT_TRUE(born_area / t < 1.1f * limit_slope);
      HS_EXPECT_TRUE(max_centroid_shift < 1e-4f);
      HS_EXPECT_TRUE(max_vertex_shift < 2.0f * t);
      prev_area = born_area;
    }
  }
}

/**
 * @brief Steps a chamfer sweep from T_EPS to the shipping arrival on every
 *        chamfer seed, asserting topology and outward orientation hold.
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
    std::vector<math::Vector> prev_vertices;
    std::vector<math::Vector> prev_normal;
    std::vector<size_t> off;
    float max_step = 0.0f;
    float min_outward = 1e9f;
    float min_normal_dot = 1e9f;
    for (int s = 0; s < SAMPLES; ++s) {
      const float t =
          ConwayGraph::T_EPS + (CHAMFER_T_STAR - ConwayGraph::T_EPS) *
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
        HS_EXPECT_SIZE_OR_RETURN(swept.vertices, v0);
        HS_EXPECT_SIZE_OR_RETURN(swept.face_counts, f0);
        HS_EXPECT_SIZE_OR_RETURN(swept.faces, i0);
        HS_EXPECT_SIZE_OR_RETURN(compiled.face_counts, compiled0);
      }
      check_face_counts_consistent(swept);
      check_indices_in_range(swept);
      check_all_unit_vertices(swept, 1e-3f);
      conway_tests::check_euler_characteristic_two(swept);

      face_offsets(swept, off);
      std::vector<math::Vector> normal(swept.face_counts.size());
      for (size_t f = 0; f < swept.face_counts.size(); ++f) {
        normal[f] = face_area_vector(swept, off[f], swept.face_counts[f]);
        math::Vector c(0, 0, 0);
        const int n = swept.face_counts[f];
        for (int k = 0; k < n; ++k)
          c = c + swept.vertices[swept.faces[off[f] + k]];
        const float len = std::sqrt(math::dot(normal[f], normal[f])) *
                          std::sqrt(math::dot(c, c));
        if (len > 0.0f)
          min_outward = std::min(min_outward, math::dot(normal[f], c) / len);
      }
      if (s > 0) {
        for (size_t v = 0; v < swept.vertices.size(); ++v) {
          const math::Vector d = swept.vertices[v] - prev_vertices[v];
          max_step = std::max(max_step, std::sqrt(math::dot(d, d)));
        }
        for (size_t f = 0; f < normal.size(); ++f) {
          const float la = std::sqrt(math::dot(prev_normal[f], prev_normal[f]));
          const float lb = std::sqrt(math::dot(normal[f], normal[f]));
          if (la > 0.0f && lb > 0.0f)
            min_normal_dot =
                std::min(min_normal_dot,
                         math::dot(prev_normal[f], normal[f]) / (la * lb));
        }
      }
      prev_vertices.assign(swept.vertices.data(),
                           swept.vertices.data() + swept.vertices.size());
      prev_normal = normal;
    }

    // No face turns inside out and no vertex teleports between steps.
    HS_EXPECT_TRUE(min_outward > 0.0f);
    HS_EXPECT_TRUE(min_normal_dot > 0.0f);
    HS_EXPECT_LT(max_step, 0.1f);
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
// Truncate sub-T_EPS birth: the truncate001 recipes arrive below the T_EPS
// birth floor.
// ---------------------------------------------------------------------------

/** @brief One truncate-leg seed. */
struct TruncateSite {
  const char *name;                     /**< Recipe the leg belongs to. */
  PolyMesh (*seed)(Arena &a, Arena &b); /**< Chain prefix up to the truncate. */
};

inline PolyMesh probe_ticosa_ambo_relax(Arena &a, Arena &b) {
  return recipe_step_seed<
      Solids::TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN59_RECIPE,
      Solids::Op::TRUNCATE>(a, b);
}

inline constexpr TruncateSite TRUNCATE_SITES[] = {
    {"truncatedIcosahedron_ambo_relax_truncate001", probe_ticosa_ambo_relax},
};

/** Arrival parameter of the truncate001 recipes. */
inline const float TRUNCATE001_T_STAR = [] {
  constexpr auto &RECIPE =
      Solids::TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE001_HANKIN59_RECIPE;
  constexpr size_t CAPACITY = Solids::lowered_step_count(RECIPE);
  Solids::OpStep lowered[CAPACITY];
  const size_t count = Solids::expand_to_primitives(RECIPE, lowered, CAPACITY);
  for (size_t i = 0; i < count; ++i)
    if (lowered[i].op == Solids::Op::TRUNCATE)
      return lowered[i].param;
  return 0.0f;
}();

/**
 * @brief Steps the truncate001 leg from its derived birth floor to its
 *        arrival, asserting topology holds and no face collapses or inverts.
 * @details The birth floor is OpLeg's recipe-step truncate_birth_floor.
 */
inline void test_truncate001_birth_sweep_holds_topology() {
  constexpr int SAMPLES = 32;
  const float birth = ConwayGraph::truncate_birth_floor(TRUNCATE001_T_STAR);
  // A real animation, not a still image.
  HS_EXPECT_TRUE(birth < TRUNCATE001_T_STAR);
  HS_EXPECT_TRUE(TRUNCATE001_T_STAR >= ConwayGraph::T_TRUNCATE_ARRIVAL_MIN);

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
    std::vector<math::Vector> prev_normal;
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
        HS_EXPECT_SIZE_OR_RETURN(swept.vertices, v0);
        HS_EXPECT_SIZE_OR_RETURN(swept.face_counts, f0);
        HS_EXPECT_SIZE_OR_RETURN(swept.faces, i0);
        HS_EXPECT_SIZE_OR_RETURN(compiled.face_counts, compiled0);
      }
      check_face_counts_consistent(swept);
      check_indices_in_range(swept);
      check_all_unit_vertices(swept, 1e-3f);
      conway_tests::check_euler_characteristic_two(swept);

      face_offsets(swept, off);
      std::vector<math::Vector> normal(swept.face_counts.size());
      for (size_t f = 0; f < swept.face_counts.size(); ++f) {
        normal[f] = face_area_vector(swept, off[f], swept.face_counts[f]);
        min_area =
            std::min(min_area, std::sqrt(math::dot(normal[f], normal[f])));
        math::Vector c(0, 0, 0);
        const int n = swept.face_counts[f];
        for (int k = 0; k < n; ++k)
          c = c + swept.vertices[swept.faces[off[f] + k]];
        const float len = std::sqrt(math::dot(normal[f], normal[f])) *
                          std::sqrt(math::dot(c, c));
        if (len > 0.0f)
          min_outward = std::min(min_outward, math::dot(normal[f], c) / len);
      }
      if (s > 0) {
        for (size_t f = 0; f < normal.size(); ++f) {
          const float la = std::sqrt(math::dot(prev_normal[f], prev_normal[f]));
          const float lb = std::sqrt(math::dot(normal[f], normal[f]));
          if (la > 0.0f && lb > 0.0f)
            min_normal_dot =
                std::min(min_normal_dot,
                         math::dot(prev_normal[f], normal[f]) / (la * lb));
        }
      }
      prev_normal = normal;
    }

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
// Far-side truncate sweep: the truncate50d recipes arrive past the ambo pinch
// (t = 0.5), where truncate short-circuits to ambo. Past 0.5 the cut faces
// self-intersect, so only structural invariants hold.
// ---------------------------------------------------------------------------

inline PolyMesh probe_ticosidodeca(Arena &a, Arena &b) {
  return Solids::Archimedean::truncatedIcosidodecahedron(a, b);
}

inline constexpr TruncateSite FAR_TRUNCATE_SITES[] = {
    {"truncatedIcosahedron_truncate50d_ambo_dual", probe_ticosa},
    {"truncatedIcosidodecahedron_truncate50d_ambo_dual", probe_ticosidodeca},
};

/** Arrival parameter of the truncate50d recipes, past the pinch. */
inline constexpr float TRUNCATE50D_T_STAR =
    Solids::IslamicStarPatterns::TRUNCATE_T_FAR;
/** Below this t, the truncate cut faces do not yet self-intersect, so positive
 * area still holds and can be asserted. Above the pinch (0.5) it cannot. */
inline constexpr float FAR_SIDE_NEAR_LIMIT = 0.49f;

/**
 * @brief Steps the far-side truncate leg from its near-side birth floor through
 *        the ambo pinch on both truncate50d seeds, asserting the leg does not
 *        trap and does not change topology across the pinch.
 * @details Samples go through ConwayGraph::truncate_off_pinch.
 */
inline void test_truncate50d_far_side_sweep_holds_topology() {
  constexpr int SAMPLES = 48;
  const float birth = ConwayGraph::truncate_birth_floor(TRUNCATE50D_T_STAR);
  HS_EXPECT_TRUE(birth < 0.5f);
  HS_EXPECT_TRUE(TRUNCATE50D_T_STAR > 0.5f);
  HS_EXPECT_TRUE(TRUNCATE50D_T_STAR <= ConwayGraph::T_TRUNCATE_FAR_MAX);
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
    for (int s = 0; s < SAMPLES; ++s) {
      // Linear birth -> arrival, plus one sample forced onto the exact pinch so
      // the guard is exercised deterministically.
      float raw = (s == SAMPLES / 2)
                      ? 0.5f
                      : birth + (TRUNCATE50D_T_STAR - birth) *
                                    (static_cast<float>(s) / (SAMPLES - 1));
      const float t = ConwayGraph::truncate_off_pinch(raw);
      HS_EXPECT_TRUE(t != 0.5f);

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
        HS_EXPECT_SIZE_OR_RETURN(swept.vertices, v0);
        HS_EXPECT_SIZE_OR_RETURN(swept.face_counts, f0);
        HS_EXPECT_SIZE_OR_RETURN(swept.faces, i0);
        HS_EXPECT_SIZE_OR_RETURN(compiled.face_counts, compiled0);
      }
      check_face_counts_consistent(swept);
      check_indices_in_range(swept);
      check_all_unit_vertices(swept, 1e-3f);
      int nonfinite = 0;
      for (size_t i = 0; i < swept.vertices.size(); ++i)
        nonfinite += !std::isfinite(swept.vertices[i].length());
      HS_EXPECT_EQ(nonfinite, 0);

      // Positive area holds only on the near side of the pinch.
      if (t <= FAR_SIDE_NEAR_LIMIT) {
        face_offsets(swept, off);
        for (size_t f = 0; f < swept.face_counts.size(); ++f) {
          const math::Vector n =
              face_area_vector(swept, off[f], swept.face_counts[f]);
          const math::Vector c =
              face_centroid_unit(swept, off[f], swept.face_counts[f]);
          min_area_near = std::min(min_area_near, math::dot(n, c));
        }
      }
    }

    HS_EXPECT_TRUE(min_area_near > 0.0f);
    if (hs_test::stats().failed != failed_before)
      std::printf("    [truncate50d] %s failed (raw F=%zu compiled=%zu)\n",
                  site.name, f0, compiled0);
    else
      std::printf("  [truncate50d] %s: birth=%.4f -> %.4f through pinch V=%zu "
                  "F=%zu I=%zu compiled=%zu min_area_near=%.3e\n",
                  site.name, static_cast<double>(birth),
                  static_cast<double>(TRUNCATE50D_T_STAR), v0, f0, i0,
                  compiled0, static_cast<double>(min_area_near));
  }
}

/**
 * @brief Finds the smallest chamfer t at which no newborn hexagon is rejected
 *        by SDF::Face's collapsed-area cull, per seed.
 * @details Bisects the SDF::Face COLLAPSED_AREA_RATIO boundary and bounds it
 * against T_EPS, which must clear the cull outright.
 */
inline void test_chamfer_birth_epsilon() {
  constexpr float MAX_BIRTH_EPSILON = 1e-4f;
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
    // The sweeps open at T_EPS, so every birth must clear the cull there.
    HS_EXPECT_EQ(at_eps, static_cast<size_t>(0));
    HS_EXPECT_LT(hi, MAX_BIRTH_EPSILON);
    std::printf("  [chamfer-eps] %s: births clear the SDF cull at t>=%.2e "
                "(T_EPS=%.3f culls %zu of %zu; compile keeps %zu of %zu)\n",
                site.name, static_cast<double>(hi),
                static_cast<double>(ConwayGraph::T_EPS), at_eps, raw_faces - F,
                cf_eps, raw_faces);
  }
}

// ---------------------------------------------------------------------------
// OpLeg provenance tolerance: core/animation/opleg.h.
// ---------------------------------------------------------------------------

/** @brief One shipping build chain, lowered to primitive steps. */
struct ChainSite {
  const char *name;    /**< Registry entry the chain mirrors. */
  uint8_t seed;        /**< simple_registry index. */
  const OpStep *steps; /**< Lowered primitive chain. */
  size_t count;        /**< Number of steps. */
};

template <typename Fn> inline void for_each_shipping_chain(Fn &&fn) {
  constexpr size_t MAX_STEPS =
      Solids::max_lowered_step_count(Solids::islamic_registry);
  for (const Solids::Entry &entry : Solids::Collections::get_islamic_solids()) {
    OpStep steps[MAX_STEPS];
    const size_t count =
        Solids::expand_to_primitives(*entry.recipe, steps, MAX_STEPS);
    fn(ChainSite{entry.name, entry.recipe->seed, steps, count});
  }
}

/** @brief Nearest and second-nearest chord distance from c into `pts`. */
inline void two_nearest(const math::Vector &c,
                        const std::vector<math::Vector> &pts, size_t &best,
                        float &d1, float &d2) {
  best = 0;
  d1 = 1e9f;
  d2 = 1e9f;
  for (size_t j = 0; j < pts.size(); ++j) {
    const math::Vector d = c - pts[j];
    const float dsq = math::dot(d, d);
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
 * @brief Pins face-centroid spacing per intermediate build-chain mesh against
 *        OpLeg's PROVENANCE_TOL_SQ chord radius.
 * @details The ratio is TOL over half the tightest centroid spacing, so a
 * denser chain mesh raises it. The cap is the measured worst plus margin.
 */
inline void test_build_chain_centroid_spacing() {
  const float TOL = std::sqrt(Animation::OpLeg::PROVENANCE_TOL_SQ);
  constexpr float MAX_MEASURED_TOL_RATIO = 5.17f;
  constexpr float TOL_RATIO_MARGIN = 0.25f;
  constexpr float MAX_TOL_RATIO = MAX_MEASURED_TOL_RATIO + TOL_RATIO_MARGIN;
  float worst_ratio = 0.0f;
  const char *worst_name = "";
  size_t max_faces = 0;

  for_each_shipping_chain([&](const ChainSite &site) {
    Arena a(probe_a_buf, sizeof(probe_a_buf));
    Arena b(probe_b_buf, sizeof(probe_b_buf));
    for (size_t k = 0; k <= site.count; ++k) {
      ScratchScope ga(a);
      ScratchScope gb(b);
      PolyMesh m = Solids::build_steps(site.seed, site.steps, k, a, b);
      std::vector<math::Vector> c;
      face_centroids(m, c);
      std::vector<float> nn(c.size(), 0.0f);
      float min_spacing = 1e9f;
      for (size_t f = 0; f < c.size(); ++f) {
        float best = 1e9f;
        for (size_t g = 0; g < c.size(); ++g) {
          if (g == f)
            continue;
          const math::Vector d = c[f] - c[g];
          best = std::min(best, math::dot(d, d));
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
  });
  std::printf("  [spacing] worst tol/half-min = %.2f on %s; largest chain "
              "mesh F=%zu\n",
              static_cast<double>(worst_ratio), worst_name, max_faces);
  HS_EXPECT_LT(worst_ratio, MAX_TOL_RATIO);
}

/** @brief Pins per-leg face counts and maximum nearest-centroid ambiguity. */
inline void test_build_chain_provenance_ambiguity() {
  constexpr float TOL_SQ = Animation::OpLeg::PROVENANCE_TOL_SQ;
  struct ExpectedLeg {
    const char *name;
    size_t leg, previous, total;
    float ratio;
  };
  constexpr ExpectedLeg EXPECTED[] = {
      {"dodecahedron_hk62_ambo_hk62", 0, 12, 32, 1.000f},
      {"dodecahedron_hk62_ambo_hk62", 1, 32, 62, 1.000f},
      {"dodecahedron_hk62_ambo_hk62", 2, 62, 182, 0.566f},
      {"truncatedIcosahedron_hk58_chamfer63", 0, 32, 92, 0.831f},
      {"truncatedIcosahedron_hk58_chamfer63", 1, 92, 452, 0.784f},
      {"dodecahedron_ambo_bevel33_relax_hk66", 0, 12, 32, 1.000f},
      {"dodecahedron_ambo_bevel33_relax_hk66", 1, 32, 62, 1.000f},
      {"dodecahedron_ambo_bevel33_relax_hk66", 2, 62, 122, 0.602f},
      {"dodecahedron_ambo_bevel33_relax_hk66", 4, 122, 362, 0.679f},
      {"truncatedIcosahedron_ambo_relax_truncate33_hk64", 0, 32, 92, 0.844f},
      {"truncatedIcosahedron_ambo_relax_truncate33_hk64", 2, 92, 182, 1.000f},
      {"truncatedIcosahedron_ambo_relax_truncate33_hk64", 3, 182, 542, 0.879f},
      {"dodecahedron_bevel2_relax_gyro", 0, 12, 32, 1.000f},
      {"dodecahedron_bevel2_relax_gyro", 1, 32, 62, 1.000f},
      {"dodecahedron_bevel2_relax_gyro", 3, 62, 542, 0.868f},
      {"truncatedIcosidodecahedron_bevel5_relax_hk77", 0, 62, 182, 0.703f},
      {"truncatedIcosidodecahedron_bevel5_relax_hk77", 1, 182, 362, 1.000f},
      {"truncatedIcosidodecahedron_bevel5_relax_hk77", 3, 362, 722, 0.798f},
      {"truncatedOctahedron_gyro_kis_hk17", 0, 14, 110, 0.990f},
      {"truncatedOctahedron_gyro_kis_hk17", 3, 360, 542, 1.000f},
      {"truncatedIcosahedron_ambo_relax_truncate001_hankin59", 0, 32, 92,
       0.844f},
      {"truncatedIcosahedron_ambo_relax_truncate001_hankin59", 2, 92, 182,
       1.000f},
      {"truncatedIcosahedron_ambo_relax_truncate001_hankin59", 3, 182, 542,
       0.181f},
      {"truncatedIcosahedron_ambo_relax_truncate001_hankin73", 0, 32, 92,
       0.844f},
      {"truncatedIcosahedron_ambo_relax_truncate001_hankin73", 2, 92, 182,
       1.000f},
      {"truncatedIcosahedron_ambo_relax_truncate001_hankin73", 3, 182, 542,
       0.181f},
      {"icosahedron_ambo_truncate033_hankin59", 0, 20, 32, 1.000f},
      {"icosahedron_ambo_truncate033_hankin59", 1, 32, 62, 1.000f},
      {"icosahedron_ambo_truncate033_hankin59", 2, 62, 182, 0.835f},
      {"dodecahedron_hk35_ambo_hk62_ambo_relax_hk42", 0, 12, 32, 1.000f},
      {"dodecahedron_hk35_ambo_hk62_ambo_relax_hk42", 1, 32, 62, 1.000f},
      {"dodecahedron_hk35_ambo_hk62_ambo_relax_hk42", 2, 62, 182, 0.629f},
      {"dodecahedron_hk35_ambo_hk62_ambo_relax_hk42", 3, 182, 362, 1.000f},
      {"dodecahedron_hk35_ambo_hk62_ambo_relax_hk42", 5, 362, 1082, 0.755f},
      {"octahedron_hk17_ambo_hk73", 0, 8, 14, 1.000f},
      {"octahedron_hk17_ambo_hk73", 1, 14, 26, 1.000f},
      {"octahedron_hk17_ambo_hk73", 2, 26, 74, 0.589f},
      {"icosahedron_kis_gyro", 1, 60, 272, 1.000f},
      {"truncatedIcosidodecahedron_truncate50d_ambo_dual", 0, 62, 182, 0.703f},
      {"truncatedIcosidodecahedron_truncate50d_ambo_dual", 1, 182, 542, 0.320f},
      {"icosidodecahedron_truncate5d_ambo_dual", 0, 32, 62, 1.000f},
      {"icosidodecahedron_truncate5d_ambo_dual", 1, 62, 182, 0.172f},
      {"snubDodecahedron_truncate5d_ambo_dual", 0, 92, 152, 0.997f},
      {"snubDodecahedron_truncate5d_ambo_dual", 1, 152, 452, 0.195f},
      {"octahedron_hk34_ambo_hk72", 0, 8, 14, 1.000f},
      {"octahedron_hk34_ambo_hk72", 1, 14, 26, 1.000f},
      {"octahedron_hk34_ambo_hk72", 2, 26, 74, 0.556f},
      {"rhombicuboctahedron_hk63_ambo_hk63", 0, 26, 50, 0.776f},
      {"rhombicuboctahedron_hk63_ambo_hk63", 1, 50, 98, 1.000f},
      {"rhombicuboctahedron_hk63_ambo_hk63", 2, 98, 290, 0.766f},
      {"truncatedIcosahedron_hk54_ambo_hk72", 0, 32, 92, 0.831f},
      {"truncatedIcosahedron_hk54_ambo_hk72", 1, 92, 182, 1.000f},
      {"truncatedIcosahedron_hk54_ambo_hk72", 2, 182, 542, 0.590f},
      {"dodecahedron_hk54_ambo_hk72", 0, 12, 32, 1.000f},
      {"dodecahedron_hk54_ambo_hk72", 1, 32, 62, 1.000f},
      {"dodecahedron_hk54_ambo_hk72", 2, 62, 182, 0.576f},
      {"dodecahedron_hk72_ambo_dual_hk20", 0, 12, 32, 1.000f},
      {"dodecahedron_hk72_ambo_dual_hk20", 1, 32, 62, 1.000f},
      {"dodecahedron_hk72_ambo_dual_hk20", 3, 120, 182, 1.000f},
      {"truncatedIcosahedron_truncate50d_ambo_dual", 0, 32, 92, 0.844f},
      {"truncatedIcosahedron_truncate50d_ambo_dual", 1, 92, 272, 0.209f},
      {"icosahedron_snub_relax_truncate033_hankin62", 0, 20, 92, 1.000f},
      {"icosahedron_snub_relax_truncate033_hankin62", 2, 92, 152, 0.997f},
      {"icosahedron_snub_relax_truncate033_hankin62", 3, 152, 452, 0.955f},
  };
  size_t checked = 0;
  size_t max_prev_faces = 0;
  const char *max_prev_name = "";
  float worst_newborn_ratio = 0.0f;
  float worst_prefix_offset = 0.0f;
  size_t prefix_legs = 0, full_legs = 0, misidentified = 0;

  for_each_shipping_chain([&](const ChainSite &site) {
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
      std::vector<math::Vector> prev_c;
      face_centroids(seed, prev_c);
      const size_t prev_faces = prev_c.size();
      if (prev_faces > max_prev_faces) {
        max_prev_faces = prev_faces;
        max_prev_name = site.name;
      }

      // The mesh the leg's first frame draws.
      PolyMesh start;
      switch (step.op) {
      case Op::HANKIN: {
        CompiledHankin hk;
        MeshOps::compile_hankin(seed, hk, b, a, true);
        MeshOps::update_hankin(
            hk, start, a, std::max(step.param, Animation::OpLeg::THETA_EPS));
        const size_t statics = hk.static_vertices.size();
        for (size_t i = statics; i < start.vertices.size(); ++i) {
          const math::Vector corner =
              hk.corner(hk.dynamic_instructions[i - statics].v_corner);
          const math::Vector arrival =
              math::Snorm3::encode(start.vertices[i]).decode();
          start.vertices[i] = math::slerp(corner.normalized(), arrival,
                                          Animation::OpLeg::K_EPS);
        }
        break;
      }
      case Op::AMBO:
        start = MeshOps::truncate(seed, a, b, ConwayGraph::T_EPS);
        break;
      case Op::CHAMFER:
        start = MeshOps::chamfer(seed, a, b, ConwayGraph::T_EPS);
        break;
      case Op::TRUNCATE:
        start = MeshOps::truncate(
            seed, a, b, ConwayGraph::truncate_birth_floor(step.param));
        break;
      case Op::SNUB:
        start = MeshOps::snub(seed, a, b, ConwayGraph::T_EPS, 0.0f);
        break;
      default:
        continue;
      }
      std::vector<math::Vector> start_c;
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

      // Every newborn bounds production's first-newborn-per-class ambiguity.
      float max_ratio = 0.0f;
      float max_newborn_d1 = 0.0f;
      for (size_t f = prev_faces; f < total; ++f) {
        size_t best;
        float d1, d2;
        two_nearest(start_c[f], prev_c, best, d1, d2);
        max_newborn_d1 = std::max(max_newborn_d1, d1);
        if (d2 > 0.0f)
          max_ratio = std::max(max_ratio, d1 / d2);
      }
      if (total > prev_faces)
        worst_newborn_ratio = std::max(worst_newborn_ratio, max_ratio);

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
                  static_cast<double>(total > prev_faces ? max_ratio : 1.0f));
      const auto expected = std::find_if(
          std::begin(EXPECTED), std::end(EXPECTED),
          [&](const ExpectedLeg &entry) {
            return std::strcmp(entry.name, site.name) == 0 && entry.leg == k;
          });
      HS_EXPECT_TRUE(expected != std::end(EXPECTED));
      if (expected != std::end(EXPECTED)) {
        ++checked;
        HS_EXPECT_EQ(prev_faces, expected->previous);
        HS_EXPECT_EQ(total, expected->total);
        HS_EXPECT_LE(max_ratio, expected->ratio + 0.01f);
      }
      // The prefix identity must also be the geometric nearest, or the
      // newborn-class inheritance is reading the wrong departed face.
      HS_EXPECT_EQ(leg_misidentified, static_cast<size_t>(0));
      HS_EXPECT_LT(max_prefix_offset * max_prefix_offset, TOL_SQ);
    }
  });
  std::printf("  [prov] %zu prefix legs, %zu full-correspondence legs, %zu "
              "misidentified prefix faces; max prev_faces=%zu on %s; "
              "worst newborn d1/d2=%.3f; worst prefix "
              "offset=%.4f (tol %.2f)\n",
              prefix_legs, full_legs, misidentified, max_prev_faces,
              max_prev_name, static_cast<double>(worst_newborn_ratio),
              static_cast<double>(worst_prefix_offset),
              static_cast<double>(std::sqrt(TOL_SQ)));
  // None of the swept ops scanned here takes the tolerance-checked
  // full-correspondence path, so PROVENANCE_TOL_SQ never applies to them.
  HS_EXPECT_EQ(full_legs, static_cast<size_t>(0));
  HS_EXPECT_EQ(checked, std::size(EXPECTED));
  HS_EXPECT_TRUE(max_prev_faces > 128);
}

// ---------------------------------------------------------------------------
// Needle primitive lowering on the hankin(54 deg) test seed.
// ---------------------------------------------------------------------------

/**
 * @brief Checks two-face edge incidence, Euler characteristic 2, positive-area
 *        faces, and near-unit vertices; returns the compiled face count.
 */
inline size_t check_manifold_landing(const PolyMesh &m, Arena &a, Arena &b) {
  check_face_counts_consistent(m);
  check_indices_in_range(m);
  check_all_unit_vertices(m, 1e-3f);
  conway_tests::check_euler_characteristic_two(m);
  std::vector<size_t> off;
  face_offsets(m, off);
  float min_area = 1e9f;
  for (size_t f = 0; f < m.face_counts.size(); ++f) {
    const math::Vector n = face_area_vector(m, off[f], m.face_counts[f]);
    const math::Vector c = face_centroid_unit(m, off[f], m.face_counts[f]);
    min_area = std::min(min_area, math::dot(n, c));
  }
  for (size_t v = 0; v < m.vertices.size(); ++v)
    HS_EXPECT_TRUE(std::isfinite(m.vertices[v].length()));
  HS_EXPECT_TRUE(min_area > 0.0f);
  MeshState compiled;
  MeshOps::compile(m, compiled, a, b);
  return compiled.face_counts.size();
}

/**
 * @brief Builds needle's DUAL then KIS primitive results on the hankin(54 deg)
 *        test seed, asserting two-face edge incidence, Euler characteristic 2,
 *        and exact element-wise composite parity.
 */
inline void test_needle_partition_lowering_builds_on_hankin() {
  const int failed_before = hs_test::stats().failed;

  Arena persist(probe_seed_buf, sizeof(probe_seed_buf));
  Arena a(probe_a_buf, sizeof(probe_a_buf));
  Arena b(probe_b_buf, sizeof(probe_b_buf));

  PolyMesh seed;
  {
    ScratchScope ga(a);
    ScratchScope gb(b);
    seed =
        Solids::finalize_solid(build_ticosa_ambo_relax100_hk54(a, b), persist);
  }
  size_t seed_compiled = 0;
  {
    ScratchScope fa(a);
    ScratchScope fb(b);
    MeshState c;
    MeshOps::compile(seed, c, a, b);
    seed_compiled = c.face_counts.size();
  }

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
  HS_EXPECT_EQ(dual_compiled, dual_mesh.face_counts.size());

  // KIS on the dual produces the needle result.
  PolyMesh kis_mesh;
  size_t kis_compiled = 0;
  {
    ScratchScope fa(a);
    ScratchScope fb(b);
    kis_mesh = Solids::finalize_solid(MeshOps::kis(dual_mesh, a, b), persist);
    kis_compiled = check_manifold_landing(kis_mesh, a, b);
  }
  HS_EXPECT_EQ(kis_compiled, kis_mesh.face_counts.size());
  // kis raises one triangle per parent face-side, so its face count is the
  // dual's total face-index count.
  HS_EXPECT_EQ(kis_mesh.face_counts.size(), dual_mesh.faces.size());

  // The lowering matches the composite: needle n = kd = kis of dual.
  {
    ScratchScope fa(a);
    ScratchScope fb(b);
    PolyMesh needle_mesh = MeshOps::needle(seed, a, b);
    check_meshes_identical(needle_mesh, kis_mesh);
  }

  if (hs_test::stats().failed != failed_before)
    std::printf(
        "    [needle] truncatedIcosahedron_ambo_relax100_hk54 failed\n");
  else
    std::printf(
        "  [needle] hk54 seed F=%zu(compiled %zu) -> dual F=%zu(%zu) -> "
        "kis F=%zu(%zu) == needle; closed manifold at each landing\n",
        seed.face_counts.size(), seed_compiled, dual_mesh.face_counts.size(),
        dual_compiled, kis_mesh.face_counts.size(), kis_compiled);
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

  test_needle_partition_lowering_builds_on_hankin();

  test_build_chain_centroid_spacing();
  test_build_chain_provenance_ambiguity();

  return fixture.result();
}

} // namespace opchain_probe_tests
} // namespace hs_test
