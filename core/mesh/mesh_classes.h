/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <cfloat>
#include <cmath>
#include <cstdint>

#include "math/3dmath.h"
#include "engine/memory.h"
#include "engine/platform.h"
#include "render/sdf.h"
#include "mesh/spatial.h"

/**
 * @file mesh_classes.h
 * @brief Congruence-class clustering + canonical distance-LUT bake.
 *
 * Clusters a mesh's faces into congruence classes (geometric clustering seeded
 * per topology class) and bakes one signed-distance LUT per class
 * (SDF::build_canonical_distance_lut); Face::distance then serves sign-pure
 * probes from a bilinear lookup instead of the exact edge walk.
 *
 * ONLY FOR EFFECTS WHOSE MESHES HOLD STILL between spawns (rigid orientation is
 * fine — congruence is frame-invariant). Per-frame deformation breaks the
 * canonical premise, so SDF::ALIGN_MAX_DEV_DIAGS drops deviating faces.
 *
 * Clustering is never depended on: an unassigned face (NO_CLASS), a class
 * without a LUT, or a degenerate alignment all degrade to the exact path.
 */
namespace MeshOps {

/** Sentinel class id: face keeps the per-face exact path. */
inline constexpr uint8_t NO_CLASS = 0xFF;
/** Congruence-class capacity per mesh (census max over the registry is 24;
 *  overflow degrades the excess faces to NO_CLASS, it never traps). */
inline constexpr int MAX_CONGRUENCE_CLASSES = 32;
/** Procrustes RMS residual (pixels) under which two faces are congruent. */
inline constexpr float CONGRUENCE_EPS_PX = 0.25f;
/** Target LUT cell diagonal in pixels: the sign-unsafe fallback band is one
 *  cell diagonal, so this pins it to a small fraction of a pixel. */
inline constexpr float LUT_TARGET_DIAG_PX = 0.35f;
/** Per-class LUT grid resolution bounds. */
inline constexpr int CLASS_LUT_MIN_N = 32;
inline constexpr int CLASS_LUT_MAX_N = 64;
/** Per-mesh LUT byte budget. Classes are allocated by descending face count
 *  until it is spent; the remainder run NO_CLASS. Identical on every target
 *  (host/WASM/device) so sim and device output cannot fork. Sized so a
 *  double-buffered pair of bakes plus a palette bank fits a 136 KB device
 *  persistent partition; a consumer effect must add its bakes to its
 *  test_solids-style persistent-budget sweep. */
inline constexpr size_t CLASS_LUT_BUDGET = 18 * 1024;
/** Minimum bake-predicted hit share for a class LUT to be kept: below this
 *  the probes mostly land in the fallback band and pay the guard for
 *  nothing (small faces relative to the cell diagonal). */
inline constexpr float MIN_CLASS_HIT_SHARE = 0.4f;

/**
 * @brief Per-face congruence record, baked once per spawned mesh.
 */
struct FaceClassRec {
  uint8_t class_id;    /**< Congruence class, or NO_CLASS. */
  uint8_t vert_offset; /**< Cyclic offset aligning mesh vertex order to canonical. */
  uint8_t reflected;   /**< Non-zero for the mirror family. */
};

/**
 * @brief One congruence class: canonical shape + optional distance LUT.
 */
struct CongruenceClass {
  SDF::ClassLut lut;              /**< Distance LUT; data == nullptr when the class has none. */
  const float *canon_xy = nullptr; /**< Canonical centered 2D polygon, x/y pairs. */
  float canon_sq = 0.0f;          /**< Sum of squared canonical vertex magnitudes. */
  int n_verts = 0;                /**< Vertex count. */
  int topo_id = -1;               /**< Seeding topology class. */
  uint16_t members = 0;           /**< Faces assigned to this class. */
  bool concave = false;           /**< Canonical shape is concave (LUT-eligible). */
};

/**
 * @brief Congruence-class bake for one spawned mesh (persistent-arena owned).
 * @details Lives in the effect's persistent arena and dies at every arena
 * compaction — rebake unconditionally after compact_keep_front or the fading
 * mesh silently degrades to NO_CLASS (correct but slower).
 */
struct MeshClassBake {
  ArenaVector<CongruenceClass> classes; /**< Congruence classes, dense ids. */
  ArenaVector<FaceClassRec> face_recs;  /**< Per-face records, mesh face order. */
  float worst_residual_px = 0.0f;   /**< Worst accepted RMS residual (pixels). */
  float predicted_hit_share = 0.0f; /**< Safe-cell share on LUT-bound faces (quality; gate >= MIN_CLASS_HIT_SHARE). */
  float lut_face_share = 0.0f;      /**< Faces bound to a LUT / all faces (coverage). */
  uint16_t shared_faces = 0;        /**< Faces in classes with >= 2 members. */
  uint16_t concave_faces = 0;       /**< Faces in concave (LUT-eligible) classes. */
  uint16_t lut_faces = 0;           /**< Faces whose class received a LUT. */
  uint16_t luts_built = 0;          /**< Classes that received a LUT. */
};

/**
 * @brief Tests a centered 2D polygon for concavity.
 * @param xy Polygon vertices, x/y pairs.
 * @param count Vertex count.
 * @return True when successive-edge turns have mixed signs (same relative-turn
 *         epsilon as Face::build_half_planes, so a class is LUT-eligible
 *         exactly when its faces would miss the convex fast path).
 */
inline bool polygon_is_concave(const float *xy, int count) {
  bool pos = false, neg = false;
  for (int i = 0; i < count; ++i) {
    int i1 = (i + 1 == count) ? 0 : i + 1;
    int i2 = (i1 + 1 == count) ? 0 : i1 + 1;
    float e1x = xy[2 * i1] - xy[2 * i], e1y = xy[2 * i1 + 1] - xy[2 * i + 1];
    float e2x = xy[2 * i2] - xy[2 * i1], e2y = xy[2 * i2 + 1] - xy[2 * i1 + 1];
    float cr = e1x * e2y - e1y * e2x;
    float scale = (e1x * e1x + e1y * e1y) * (e2x * e2x + e2y * e2y);
    if (cr * cr > 1e-12f * scale) {
      if (cr > 0)
        pos = true;
      else
        neg = true;
    }
  }
  return pos && neg;
}

/**
 * @brief Clusters a mesh's faces into congruence classes and bakes the
 *        per-class canonical distance LUTs.
 * @param mesh Compiled mesh; topology must already be populated
 *        (classify_faces_by_topology).
 * @param scratch Scratch arena (rewound before return).
 * @param persistent Arena receiving the bake (records, canonical polygons,
 *        LUTs); same lifetime as the mesh's slot.
 * @param pixel_width One pixel in gnomonic plane units (2*pi/W) — sets the
 *        congruence epsilon and the LUT resolution target.
 * @param out Freshly default-constructed bake to populate.
 * @details Greedy clustering seeded per topology class: a face joins the
 * first class whose canonical polygon aligns (over cyclic offset x reflection
 * x optimal rotation) within CONGRUENCE_EPS_PX RMS, else founds a new class
 * from its own centered projection. Gnomonic projection about each face's own
 * centroid is position-covariant, so the clustering is valid for any mesh
 * orientation and is baked once per spawn.
 *
 * LUTs are built for concave classes with >= 2 members, largest first, until
 * CLASS_LUT_BUDGET is spent. Logs the census telemetry (classes, coverage,
 * worst residual, predicted hit share) at the end.
 */
[[maybe_unused]] HS_COLD static void
build_mesh_class_bake(const MeshState &mesh, Arena &scratch, Arena &persistent,
                      float pixel_width, MeshClassBake &out) {
  ScratchScope scratch_guard(scratch);

  const size_t F = mesh.get_face_counts_size();
  const uint8_t *fc = mesh.get_face_counts_data();
  const uint16_t *fi = mesh.get_faces_data();
  const uint16_t *fo = mesh.get_face_offsets_data();
  HS_CHECK(mesh.topology.size() == F,
           "build_mesh_class_bake requires classify_faces_by_topology first");

  out.classes.bind(persistent, MAX_CONGRUENCE_CLASSES);
  out.face_recs.bind(persistent, F);

  constexpr int MAX_VERTS = SDF::FaceScratchBuffer::MAX_VERTS;
  const float eps_plane = CONGRUENCE_EPS_PX * pixel_width;

  // Arena-hosted (not stack): this runs inside the effect spawn path, whose
  // stack high-water is budget-gated (tests/stack_measure.cpp).
  float *zx = scratch.allocate_n<float>(MAX_VERTS);
  float *zy = scratch.allocate_n<float>(MAX_VERTS);
  float worst_res_px = 0.0f;

  for (size_t f = 0; f < F; ++f) {
    const int count = fc[f];
    out.face_recs.push_back({NO_CLASS, 0, 0});
    FaceClassRec &rec = out.face_recs[f];
    if (count < 3 || count > MAX_VERTS)
      continue;

    // Gnomonic projection about the face's own centroid — the same projection
    // Face::setup_frame_and_polygon builds per frame, so the canonical shape
    // and the per-frame polygon differ only by an in-plane rotation
    // (+ reflection), which the alignment correlation absorbs.
    const uint16_t *idx = fi + fo[f];
    Vector center(0, 0, 0);
    for (int k = 0; k < count; ++k)
      center = center + mesh.vertices[idx[k]];
    center.normalize();
    Vector u = cross(center, least_parallel_axis(center)).normalized();
    Vector w = cross(center, u).normalized();
    float mx = 0.0f, my = 0.0f;
    for (int k = 0; k < count; ++k) {
      const Vector &v = mesh.vertices[idx[k]];
      float d = dot(v, center);
      if (fabsf(d) < math::TOLERANCE)
        d = copysignf(math::TOLERANCE, d);
      zx[k] = dot(v, u) / d;
      zy[k] = dot(v, w) / d;
      mx += zx[k];
      my += zy[k];
    }
    mx /= count;
    my /= count;
    float zz = 0.0f;
    for (int k = 0; k < count; ++k) {
      zx[k] -= mx;
      zy[k] -= my;
      zz += zx[k] * zx[k] + zy[k] * zy[k];
    }

    // Best alignment against every same-topology class rep.
    const int topo = mesh.topology[f];
    int best_class = -1, best_off = 0;
    bool best_refl = false;
    float best_res_sq = FLT_MAX;
    for (size_t c = 0; c < out.classes.size(); ++c) {
      const CongruenceClass &cls = out.classes[c];
      if (cls.topo_id != topo || cls.n_verts != count)
        continue;
      for (int refl = 0; refl < 2; ++refl) {
        for (int off = 0; off < count; ++off) {
          SDF::AlignCorr a = SDF::align_correlate(
              cls.canon_xy, count, off, refl != 0,
              [&](int j, float &gx, float &gy) {
                gx = zx[j];
                gy = zy[j];
              });
          float res_sq =
              std::max(cls.canon_sq + zz -
                           2.0f * sqrtf(a.rr * a.rr + a.ri * a.ri),
                       0.0f);
          if (res_sq < best_res_sq) {
            best_res_sq = res_sq;
            best_class = static_cast<int>(c);
            best_off = off;
            best_refl = refl != 0;
          }
        }
      }
    }

    if (best_class >= 0 && best_res_sq <= eps_plane * eps_plane * count) {
      rec = {static_cast<uint8_t>(best_class),
             static_cast<uint8_t>(best_off),
             static_cast<uint8_t>(best_refl ? 1 : 0)};
      out.classes[best_class].members++;
      float res_px = sqrtf(best_res_sq / count) / pixel_width;
      if (res_px > worst_res_px)
        worst_res_px = res_px;
      continue;
    }

    if (out.classes.size() >= MAX_CONGRUENCE_CLASSES)
      continue; // degrade to NO_CLASS; the exact path is always correct

    // Found a new class from this face's own centered projection.
    float *canon = persistent.allocate_n<float>(2 * count);
    for (int k = 0; k < count; ++k) {
      canon[2 * k] = zx[k];
      canon[2 * k + 1] = zy[k];
    }
    CongruenceClass cls;
    cls.canon_xy = canon;
    cls.canon_sq = zz;
    cls.n_verts = count;
    cls.topo_id = topo;
    cls.members = 1;
    cls.concave = polygon_is_concave(canon, count);
    rec = {static_cast<uint8_t>(out.classes.size()), 0, 0};
    out.classes.push_back(cls);
  }
  out.worst_residual_px = worst_res_px;

  for (size_t c = 0; c < out.classes.size(); ++c)
    if (out.classes[c].members >= 2)
      out.shared_faces += out.classes[c].members;

  // LUT build: concave shared classes only (convex faces already have the
  // ~lookup-cheap half-plane path), largest first until the budget is spent.
  int order[MAX_CONGRUENCE_CLASSES];
  int n_elig = 0;
  for (size_t c = 0; c < out.classes.size(); ++c) {
    const CongruenceClass &cls = out.classes[c];
    if (cls.members < 2 || !cls.concave)
      continue;
    int pos = n_elig++;
    while (pos > 0 &&
           out.classes[order[pos - 1]].members < cls.members) {
      order[pos] = order[pos - 1];
      --pos;
    }
    order[pos] = static_cast<int>(c);
  }

  size_t budget = CLASS_LUT_BUDGET;
  const float target_diag = LUT_TARGET_DIAG_PX * pixel_width;
  float hit_share_acc = 0.0f;
  int dropped_classes = 0, dropped_faces = 0, lowq_classes = 0;
  // Staging buffer: LUTs are built here first and promoted to the persistent
  // arena only if their predicted hit share clears the quality bar, so a
  // low-value LUT never spends persistent budget.
  int16_t *staging = scratch.allocate_n<int16_t>(CLASS_LUT_MAX_N * CLASS_LUT_MAX_N);
  for (int e = 0; e < n_elig; ++e) {
    CongruenceClass &cls = out.classes[order[e]];
    out.concave_faces += cls.members;
    // Size the grid so the cell diagonal lands on target_diag (the box gains
    // BOUNDS_MARGIN_WIDE per side inside the builder, mirrored here).
    float bmin_x = FLT_MAX, bmax_x = -FLT_MAX, bmin_y = FLT_MAX,
          bmax_y = -FLT_MAX;
    for (int k = 0; k < cls.n_verts; ++k) {
      bmin_x = std::min(bmin_x, cls.canon_xy[2 * k]);
      bmax_x = std::max(bmax_x, cls.canon_xy[2 * k]);
      bmin_y = std::min(bmin_y, cls.canon_xy[2 * k + 1]);
      bmax_y = std::max(bmax_y, cls.canon_xy[2 * k + 1]);
    }
    float Rx = (bmax_x - bmin_x) * 0.5f + SDF::BOUNDS_MARGIN_WIDE;
    float Ry = (bmax_y - bmin_y) * 0.5f + SDF::BOUNDS_MARGIN_WIDE;
    int n = static_cast<int>(
        hs::clamp(ceilf(2.0f * sqrtf(Rx * Rx + Ry * Ry) / target_diag) + 1.0f,
                  static_cast<float>(CLASS_LUT_MIN_N),
                  static_cast<float>(CLASS_LUT_MAX_N)));
    // Degrade resolution before dropping a class: a coarser grid on more
    // classes buys more served probes than a fine grid on fewer (the fallback
    // band widens with the cell diagonal, but only near the boundary).
    size_t bytes = static_cast<size_t>(n) * n * sizeof(int16_t);
    if (bytes > budget) {
      n = static_cast<int>(sqrtf(static_cast<float>(budget) / sizeof(int16_t)));
      n = std::min(n, CLASS_LUT_MAX_N); // staging buffer is kClassLutMaxN²
      bytes = static_cast<size_t>(n) * n * sizeof(int16_t);
      if (n < CLASS_LUT_MIN_N) {
        ++dropped_classes;
        dropped_faces += cls.members;
        continue;
      }
    }
    SDF::build_canonical_distance_lut(cls.canon_xy, cls.n_verts, n, staging,
                                      cls.lut);

    // Predicted hit share: fraction of cells inside the cull disk whose four
    // corners are sign-pure and beyond the interpolation guard, weighted by
    // the class's face count. A bake-time proxy for the runtime lut_hits
    // ratio on LUT-bound faces (gate: >= MIN_CLASS_HIT_SHARE).
    float circ = 0.0f;
    for (int k = 0; k < cls.n_verts; ++k) {
      float r2 = cls.canon_xy[2 * k] * cls.canon_xy[2 * k] +
                 cls.canon_xy[2 * k + 1] * cls.canon_xy[2 * k + 1];
      circ = std::max(circ, r2);
    }
    float max_dist = sqrtf(circ) + SDF::BOUNDS_MARGIN_WIDE;
    float max_dist_sq = max_dist * max_dist;
    float step_x = (2.0f * cls.lut.Rx) / (n - 1);
    float step_y = (2.0f * cls.lut.Ry) / (n - 1);
    int in_disk = 0, safe = 0;
    for (int gy = 0; gy + 1 < n; ++gy) {
      float ccy = (cls.lut.cy - cls.lut.Ry) + (gy + 0.5f) * step_y;
      for (int gx = 0; gx + 1 < n; ++gx) {
        float ccx = (cls.lut.cx - cls.lut.Rx) + (gx + 0.5f) * step_x;
        if (ccx * ccx + ccy * ccy > max_dist_sq)
          continue;
        ++in_disk;
        int q00 = staging[gy * n + gx], q10 = staging[gy * n + gx + 1];
        int q01 = staging[(gy + 1) * n + gx],
            q11 = staging[(gy + 1) * n + gx + 1];
        bool same_sign = (q00 > 0) == (q10 > 0) && (q00 > 0) == (q01 > 0) &&
                         (q00 > 0) == (q11 > 0);
        int min_q = std::min({std::abs(q00), std::abs(q10), std::abs(q01),
                              std::abs(q11)});
        if (same_sign && min_q * cls.lut.dequant > cls.lut.safe_dist)
          ++safe;
      }
    }
    float safe_frac = in_disk > 0 ? static_cast<float>(safe) / in_disk : 0.0f;

    // A class serving too few probes costs the guard on every probe and then
    // walks anyway (small faces are mostly fallback band) — keep the exact
    // path instead of paying for a LUT that rarely fires.
    if (safe_frac < MIN_CLASS_HIT_SHARE) {
      ++lowq_classes;
      cls.lut = SDF::ClassLut();
      continue;
    }

    budget -= bytes;
    int16_t *data = static_cast<int16_t *>(
        persistent.allocate(bytes, alignof(int16_t)));
    std::copy(staging, staging + static_cast<size_t>(n) * n, data);
    cls.lut.data = data;
    out.luts_built++;
    out.lut_faces += cls.members;
    hit_share_acc += cls.members * safe_frac;
  }
  out.predicted_hit_share =
      out.lut_faces > 0 ? hit_share_acc / out.lut_faces : 0.0f;
  out.lut_face_share = F > 0 ? static_cast<float>(out.lut_faces) / F : 0.0f;

  hs::log("mesh class bake: F=%d classes=%d shared=%.1f%% worst=%.3fpx "
          "concave=%d luts=%d/%d lut_faces=%.1f%% (dropped %d cls/%d faces, "
          "%d low-quality, %dB left) pred_hit=%.1f%%",
          (int)F, (int)out.classes.size(),
          F > 0 ? 100.0f * out.shared_faces / F : 0.0f, worst_res_px,
          (int)out.concave_faces, (int)out.luts_built, n_elig,
          100.0f * out.lut_face_share, dropped_classes, dropped_faces,
          lowq_classes, (int)budget, 100.0f * out.predicted_hit_share);
}

} // namespace MeshOps
