/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include "math/geometry.h"
#include "platform/constants.h"
#include "render/clip.h"
#include "render/sdf/common.h"

/**
 * @file face.h
 * @brief SDF::Face, the mesh-face leaf, and the congruence-class distance LUT
 * its probes are served from.
 */

namespace SDF {

// --- Congruence-class canonical distance LUTs --------------------------------
// Every islamic mesh's faces are near-exact copies of a handful of canonical 2D
// shapes (gnomonic projection about each face's centroid is
// position-covariant). MeshOps bakes one signed-distance LUT per congruence
// class at spawn (core/mesh/mesh_classes.h); Scan::Mesh binds it per frame via
// bind_class_lut, and Face::distance serves sign-pure probes from a bilinear
// lookup instead of the exact per-edge walk.

/** Minimum squared normalized correlation for a valid class-LUT alignment;
 *  below this the face is too deformed and keeps the exact path. */
inline constexpr float ALIGN_MIN_CORR_SQ = 0.25f;
/** Maximum per-vertex deviation from the aligned canonical shape (as a multiple
 *  of the LUT cell diagonal) before a face keeps the exact path. The facility
 *  fits only meshes that hold still per spawn (mesh_classes.h). */
inline constexpr float ALIGN_MAX_DEV_DIAGS = 0.25f;

/**
 * @brief Canonical congruence-class signed-distance LUT, baked once per
 *        spawned mesh.
 * @details Distances are in canonical gnomonic plane units, quantized to int16
 * over the LUT box diameter (step ~1e-5 plane units). The domain is the
 * canonical polygon's bounding box + BOUNDS_MARGIN_WIDE, a different region
 * from the max_dist_sq cull disk (circumradius + the same margin): a probe can
 * survive the cull and still land outside the domain, so Face::distance's grid
 * clamp is required to keep the fetch in bounds.
 */
struct ClassLut {
  const int16_t *data =
      nullptr;          /**< n*n quantized signed distances (row-major). */
  int n = 0;            /**< Grid resolution per axis. */
  float cx = 0, cy = 0; /**< Canonical bounding-box center. */
  float Rx = 0, Ry = 0; /**< Half-extents (+ margin). */
  float inv_step_x = 0; /**< Reciprocal cell width. */
  float inv_step_y = 0; /**< Reciprocal cell height. */
  float safe_dist = 0;  /**< Cell diagonal (sign-pure interpolation bound). */
  float dequant = 0;    /**< int16 -> plane-unit scale. */
};

/**
 * @brief Bakes the signed point-to-polygon distance field of a canonical
 *        centered 2D polygon into an int16 grid.
 * @param poly_xy Centered polygon vertices, x/y pairs.
 * @param count Vertex count (>= 3).
 * @param n Grid resolution per axis (>= 2).
 * @param out Storage for n*n quantized samples.
 * @param lut Receives the domain/quantization parameters, with data = out.
 * @details Exact per-edge walk with crossing-test sign, over the bounding box
 * + BOUNDS_MARGIN_WIDE. That box is not a superset of the runtime cull disk
 * (circumradius + the same margin), so probes landing outside the domain rely
 * on Face::distance's grid clamp. Quantization scale is the box diameter (an
 * upper bound on any in-box distance: the polygon meets its own bounding box),
 * giving a step of ~1e-5 plane units — far below the interpolation bound.
 */
inline void build_canonical_distance_lut(const float *poly_xy, int count, int n,
                                         int16_t *out, ClassLut &lut) {
  HS_CHECK(count >= 3,
           "build_canonical_distance_lut requires at least 3 polygon vertices");
  HS_CHECK(n >= 2, "build_canonical_distance_lut requires a grid resolution "
                   "of at least 2");
  float bb_min_x = FLT_MAX, bb_max_x = -FLT_MAX;
  float bb_min_y = FLT_MAX, bb_max_y = -FLT_MAX;
  for (int i = 0; i < count; ++i) {
    float vx = poly_xy[2 * i], vy = poly_xy[2 * i + 1];
    bb_min_x = std::min(bb_min_x, vx);
    bb_max_x = std::max(bb_max_x, vx);
    bb_min_y = std::min(bb_min_y, vy);
    bb_max_y = std::max(bb_max_y, vy);
  }
  lut.cx = (bb_min_x + bb_max_x) * 0.5f;
  lut.cy = (bb_min_y + bb_max_y) * 0.5f;
  lut.Rx = std::max((bb_max_x - bb_min_x) * 0.5f + BOUNDS_MARGIN_WIDE, 0.01f);
  lut.Ry = std::max((bb_max_y - bb_min_y) * 0.5f + BOUNDS_MARGIN_WIDE, 0.01f);
  lut.n = n;
  lut.inv_step_x = (n - 1) / (2.0f * lut.Rx);
  lut.inv_step_y = (n - 1) / (2.0f * lut.Ry);
  float step_x = (2.0f * lut.Rx) / (n - 1);
  float step_y = (2.0f * lut.Ry) / (n - 1);
  // The plane SDF is 1-Lipschitz, so a zero anywhere in a cell puts every
  // corner within one cell diagonal of it; a min corner magnitude above the
  // diagonal guarantees a sign-pure cell (safe to interpolate).
  lut.safe_dist = sqrtf(step_x * step_x + step_y * step_y);
  float dmax = 2.0f * sqrtf(lut.Rx * lut.Rx + lut.Ry * lut.Ry);
  lut.dequant = dmax / 32767.0f;
  float quant = 32767.0f / dmax;

  for (int gy = 0; gy < n; ++gy) {
    float qy = (lut.cy - lut.Ry) + gy * step_y;
    for (int gx = 0; gx < n; ++gx) {
      float qx = (lut.cx - lut.Rx) + gx * step_x;
      float d_sq = FLT_MAX;
      bool inside = false;
      for (int i = 0; i < count; ++i) {
        float vx = poly_xy[2 * i], vy = poly_xy[2 * i + 1];
        int i2 = (i + 1 == count) ? 0 : i + 1;
        float ex = poly_xy[2 * i2] - vx, ey = poly_xy[2 * i2 + 1] - vy;
        float len_sq = ex * ex + ey * ey;
        float wx = qx - vx, wy = qy - vy;
        float t = len_sq > 1e-12f
                      ? hs::clamp((wx * ex + wy * ey) / len_sq, 0.0f, 1.0f)
                      : 0.0f;
        float bx = wx - ex * t, by = wy - ey * t;
        float dsq = bx * bx + by * by;
        if (dsq < d_sq)
          d_sq = dsq;
        if ((vy > qy) != (poly_xy[2 * i2 + 1] > qy)) {
          float ix = vx + (qy - vy) * ex / ey;
          if (qx < ix)
            inside = !inside;
        }
      }
      float d = (inside ? -1.0f : 1.0f) * sqrtf(d_sq);
      out[gy * n + gx] =
          static_cast<int16_t>(hs::clamp(d * quant, -32767.0f, 32767.0f));
    }
  }
  lut.data = out;
}

/**
 * @brief Accumulated complex correlation between a canonical polygon and a
 *        centered projection (see align_correlate).
 */
struct AlignCorr {
  float rr, ri; /**< Sum of canon_k * conj(z'_k) (real, imaginary). */
  float cc, zz; /**< Power terms: sum |canon_k|^2 and sum |z'_k|^2. */
};

/**
 * @brief Pairs canonical vertices with a centered 2D projection under a cyclic
 *        vertex offset + optional reflection.
 * @tparam GetZ Accessor invoked as get_z(j, zx, zy), returning the centered
 *         projection vertex j.
 * @tparam Visit Invoked as visit(k, zx, zy) for canonical vertex k and its
 *         paired projection vertex, already conjugated for the mirror family.
 * @param count Vertex count.
 * @param vert_offset Projection index corresponding to canonical vertex 0.
 * @param reflected Mirror family: conjugate each vertex and walk the
 *        projection indices in reverse (a mirrored face winds the opposite way
 *        in a consistently-wound mesh).
 * @param get_z Centered-projection vertex accessor.
 * @param visit Per-correspondence sink.
 * @details Single source for the correspondence convention — the correlation
 * below, bake-time clustering (mesh_classes.h) and the per-frame
 * Face::bind_class_lut all route through it, so the (offset, reflected)
 * encoding cannot drift.
 */
template <typename GetZ, typename Visit>
inline void align_walk(int count, int vert_offset, bool reflected, GetZ get_z,
                       Visit visit) {
  int j = vert_offset;
  for (int k = 0; k < count; ++k) {
    float zx, zy;
    get_z(j, zx, zy);
    if (reflected) {
      zy = -zy;
      if (--j < 0)
        j = count - 1;
    } else {
      if (++j == count)
        j = 0;
    }
    visit(k, zx, zy);
  }
}

/**
 * @brief Correlates a canonical polygon against a centered 2D projection under
 *        a cyclic vertex offset + optional reflection.
 * @tparam GetZ Accessor invoked as get_z(j, zx, zy), returning the centered
 *         projection vertex j.
 * @param canon_xy Canonical centered polygon, x/y pairs, canonical order.
 * @param count Vertex count.
 * @param vert_offset Projection index corresponding to canonical vertex 0.
 * @param reflected Mirror family (see align_walk).
 * @param get_z Centered-projection vertex accessor.
 * @return The correlation sums; the least-squares residual is
 *         cc + zz - 2*|r|, and the optimal rotation is r / |r|.
 */
template <typename GetZ>
inline AlignCorr align_correlate(const float *canon_xy, int count,
                                 int vert_offset, bool reflected, GetZ get_z) {
  AlignCorr a{0.0f, 0.0f, 0.0f, 0.0f};
  align_walk(count, vert_offset, reflected, get_z,
             [&](int k, float zx, float zy) {
               float cx = canon_xy[2 * k], cy = canon_xy[2 * k + 1];
               // canon * conj(z)
               a.rr += cx * zx + cy * zy;
               a.ri += cy * zx - cx * zy;
               a.cc += cx * cx + cy * cy;
               a.zz += zx * zx + zy * zy;
             });
  return a;
}

/**
 * @brief Scratch buffer for Face computations to avoid allocations.
 */
struct FaceScratchBuffer {
  static constexpr int MAX_VERTS = 64; /**< Maximum vertices per face. */
  static constexpr size_t MAX_INTERVALS =
      2; /**< Capacity of the azimuth coverage array: one span, or two when the
              covered arc straddles theta=0. */
  std::array<Vector, MAX_VERTS + 1>
      poly_2d; /**< Projected 2D polygon (+1 entry to avoid modulo). */
  std::array<Vector, MAX_VERTS> edge_vectors; /**< Per-edge 2D vectors. */
  std::array<float, MAX_VERTS>
      edge_lengths_sq;                  /**< Per-edge squared lengths. */
  std::array<Vector, MAX_VERTS> planes; /**< Per-edge great-circle normals. */
  std::array<Interval, MAX_INTERVALS>
      intervals;                       /**< Azimuth coverage intervals. */
  std::array<float, MAX_VERTS> thetas; /**< Per-vertex azimuth angles. */
  std::array<float, MAX_VERTS>
      inv_edge_lengths_sq; /**< Reciprocal squared edge lengths. */
  std::array<float, MAX_VERTS>
      inv_edge_j; /**< Reciprocal of each edge's y-component. */
  std::array<Vector, MAX_VERTS + 1>
      verts_3d; /**< 3D vertices (+1 wrap entry). */

  /**
   * @brief Packed per-edge data for the cache-friendly distance() fallback.
   */
  struct EdgePacked {
    float vx, vy, ex, ey, inv_len_sq,
        inv_ej; /**< Edge origin, vector, reciprocals. */
    uint32_t key_vy,
        key_next_vy; /**< angle_key of this and the next vertex's y; equal when
                        the edge is degenerate in y. */
  };
  std::array<EdgePacked, MAX_VERTS> packed_edges; /**< Packed per-edge data. */

  /**
   * @brief Outward unit edge normal and line offset for the convex fast path.
   */
  struct HalfPlane {
    float nx, ny, off, pad; /**< Unit normal, offset (dist = nx*px + ny*py +
                               off), padding to a 16-byte stride. */
  };
  std::array<HalfPlane, MAX_VERTS>
      half_planes; /**< Convex-face edge half-planes. */
  std::array<float, MAX_VERTS + 1>
      pseudo_angles; /**< Unwrapped vertex pseudo-angles for the sector walk. */
  std::array<uint32_t, MAX_VERTS + 1>
      sector_keys; /**< pseudo_angles as order-preserving integer keys. */
  /** Bumped by every Face that finishes building over this buffer, so a Face
   *  holding an older value has had its geometry retargeted. */
  uint32_t claim_seq = 0;
};

static_assert(sdf_max_spans<Face>::value >= FaceScratchBuffer::MAX_INTERVALS,
              "Face's span bound must cover the scratch interval array it "
              "replays");

/**
 * @brief Order-preserving unsigned key for a non-NaN float: key(a) <= key(b)
 *        exactly when a <= b, with -0.0 and +0.0 mapping to the same key.
 * @details Lets the sector search compare in the core registers instead of
 * paying a vcmpe + vmrs FPU-to-core transfer per iteration.
 */
__attribute__((always_inline)) inline uint32_t angle_key(float x) {
  uint32_t u = std::bit_cast<uint32_t>(x);
  return (u & 0x80000000u) ? (0u - u) : (u + 0x80000000u);
}

/**
 * @brief Diamond pseudo-angle of (x, y) in [0, 4), strictly monotonic with
 *        atan2 but trig-free. Used to bin a query point into a face's angular
 *        sector without an atan2 in the per-pixel path.
 */
__attribute__((always_inline)) inline float pseudo_angle(float y, float x) {
  float d = fabsf(x) + fabsf(y);
  if (d < 1e-20f)
    return 0.0f;
  float r = y / d;
  if (y >= 0.0f)
    return (x >= 0.0f) ? r : (2.0f - r);
  return (x < 0.0f) ? (2.0f - r) : (4.0f + r);
}

/**
 * @brief Represents a planar face for SDF rendering.
 * @details Computes the 2D projection and vertical/horizontal bounds used to
 * accelerate rasterization. Every span member below views the
 * FaceScratchBuffer handed to the constructor and owns none of it, so the
 * buffer must outlive the Face AND back no other live Face: building a second
 * Face over the same buffer silently retargets the first one's geometry. Two
 * Faces in one CSG composition (SDF::Union<Face, Face>) therefore need two
 * buffers. Every build stamps the claim, and get_vertical_bounds() traps on a
 * retargeted Face once per draw, ahead of any probe.
 */
struct Face {
  Vector center; /**< Normalized face centroid (projection axis). */
  Vector basis_v, basis_u, basis_w; /**< Local tangent frame (v = center). */
  int count;                        /**< Vertex/edge count; 0 if culled. */
  float size = 0.0f;        /**< Inradius metric for AA normalization; radians
                                 unless linear_dist. */
  float radius = 0.0f;      /**< Circumradius in the 2D projection. */
  float max_dist = 0.0f;    /**< Cull radius (circumradius plus margin). */
  float max_dist_sq = 0.0f; /**< Squared cull radius. */

  std::span<Vector> poly_2d;      /**< Projected 2D polygon (+1 wrap entry). */
  std::span<Vector> edge_vectors; /**< Per-edge 2D vectors. */
  std::span<float> edge_lengths_sq;     /**< Per-edge squared lengths. */
  std::span<float> inv_edge_lengths_sq; /**< Reciprocal squared edge lengths. */
  std::span<float> inv_edge_j; /**< Reciprocal of each edge's y-component. */

  int y_min, y_max; /**< Inclusive vertical row bounds. */
  int build_height; /**< Canvas height the bounds were computed for. */
  int build_width; /**< Clip width the azimuth cull ran against; 0 if unclipped. */
  std::span<Interval> intervals; /**< Azimuth coverage intervals (radians). */
  bool full_width;               /**< True when the face spans all columns. */
  static constexpr bool is_solid =
      true; /**< Face renders as a filled region. */

  using EdgePacked =
      FaceScratchBuffer::EdgePacked;  /**< Packed per-edge record type. */
  std::span<EdgePacked> packed_edges; /**< Packed per-edge data. */
  using HalfPlane =
      FaceScratchBuffer::HalfPlane; /**< Convex edge half-plane record type. */
  std::span<const HalfPlane>
      half_planes;          /**< Outward unit edge normals (convex faces). */
  bool convex = false;      /**< 2D projection is convex; distance() takes the
                               half-plane path. */
  bool linear_dist = false; /**< Face is small enough to report plane distance
                               without the atan. */

  // Sector-walk state: a concave star-shaped-about-centroid face bins each
  // query point into its angular sector (by pseudo-angle) and walks only the
  // sector's edge and its nearest neighbors (sector_kmax on each side) instead
  // of all `count` edges.
  static constexpr int SECTOR_MIN_COUNT =
      10; /**< Below this the full walk is already cheap; skip the sector path.
           */
  static constexpr float SECTOR_MONO_TOL =
      0.05f; /**< Max vertex pseudo-angle backtrack (of 4 per turn) still binned
                on the K2 sector path; beyond it the fan is not star-shaped. */
  static constexpr int SECTOR_KMAX_MAX =
      2; /**< Widest neighbor walk build_sectors assigns. */
  static_assert(SECTOR_KMAX_MAX < SECTOR_MIN_COUNT,
                "plane_dist_sector's ring walk applies one wrap correction");
  std::span<const uint32_t>
      sector_keys; /**< angle_key of each unwrapped vertex pseudo-angle times
                      sector_sgn, count+1; weakly increasing (K2 faces dip <=
                      SECTOR_MONO_TOL). */
  float sector_base = 0.0f; /**< First unwrapped pseudo-angle, sgn-folded. */
  float sector_span = 0.0f; /**< Unsigned span (~4). */
  float sector_sgn = 1.0f;  /**< Winding: +1 CCW, -1 CW. Folded into the
                               table, base and span so the sector search
                               compares one direction. */
  int sector_kmax =
      1; /**< Neighbors walked each side: 1 (strict, bit-exact) or 2 (mildly
            bent, absorbs off-by-one binning). */
  bool sector_ok =
      false; /**< Star-shaped about centroid; sector walk usable. */

  // Congruence-class LUT binding (bind_class_lut), flattened for the probe
  // loop: one affine map takes gnomonic (px, py) straight to LUT grid
  // coordinates (centroid shift, canonical rotation/reflection, and grid
  // scale folded together), and the sign-purity guard compares raw int16
  // magnitudes against a pre-divided quantized threshold.
  const int16_t *lut_data =
      nullptr;                  /**< Class LUT samples; null = exact path. */
  int lut_n = 0;                /**< LUT grid resolution per axis. */
  int32_t lut_q_safe = 0;       /**< safe_dist in quantized units. */
  float lut_ax, lut_bx, lut_cx; /**< Grid-x affine coefficients. */
  float lut_ay, lut_by, lut_cy; /**< Grid-y affine coefficients. */
  float lut_clamp;              /**< Grid clamp bound (n - 2). */
  float lut_dequant;            /**< int16 -> plane-unit scale. */

  const FaceScratchBuffer *scratch_owner =
      nullptr; /**< Buffer the spans view; null when the face culled before
                    claiming one. */
  uint32_t scratch_claim = 0; /**< Claim this face took on that buffer. */

  /**
   * @brief Builds a face's projection, bounds, and edge data.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices from the pool.
   * @param scratch Scratch storage the spans alias; exclusive to this Face for
   *        its whole lifetime, and reusable only once the Face is dead.
   * @param h_virt Virtual row count (height plus pole offset).
   * @param height Canvas height in rows.
   * @param clip Optional render clip used to tighten the face bounds.
   */
  HS_O3_FN Face(std::span<const Vector> vertices,
                std::span<const uint16_t> indices, FaceScratchBuffer &scratch,
                int h_virt, int height, const ClipRegion *clip = nullptr)
      : build_height(height), build_width(clip ? clip->w : 0),
        full_width(true) {

    count = indices.size();
    HS_CHECK(count > 0 && count <= FaceScratchBuffer::MAX_VERTS,
             "Face: vertex count must be in (0, MAX_VERTS]");

    // Early vertical exit: a face whose latitude band (plus AA margin) maps to
    // an empty row range can never be rasterized.
    const bool phi_culled = [&] {
      HS_PROFILE_DEEP(face_phi_extent);
      return compute_phi_extent(vertices, indices, h_virt, height);
    }();
    if (phi_culled) {
      count = 0;
      y_min = 1;
      y_max = 0;
      return;
    }

    {
      HS_PROFILE_DEEP(face_project);
      setup_frame_and_polygon(vertices, indices, scratch);
    }

    // A fully collapsed polygon (signed area ~ 0, e.g. a hankin rosette at
    // angle 0 spiking out to each edge midpoint and back through its corner)
    // encloses no region, but its boundary would still rasterize as a ~1 px
    // AA line; cull it like the phi-extent reject. The residual of an exactly
    // collapsed face is float noise (< ~1e-6 of radius^2), orders of
    // magnitude under the thinnest real sliver a sweep draws, so the
    // threshold decision is identical sim/device. The compare is inclusive so
    // that coincident vertices, which zero both sides, are culled too.
    float area2 = 0.0f;
    for (int i = 0; i < count; ++i)
      area2 +=
          poly_2d[i].x * poly_2d[i + 1].y - poly_2d[i + 1].x * poly_2d[i].y;
    if (fabsf(area2) <= COLLAPSED_AREA_RATIO * radius * radius) {
      // Scratch already holds this face's geometry, so retire any earlier
      // face's claim on the way out.
      ++scratch.claim_seq;
      count = 0;
      y_min = 1;
      y_max = 0;
      return;
    }

    {
      HS_PROFILE_DEEP(face_thetas);
      compute_thetas(scratch);
    }
    {
      HS_PROFILE_DEEP(face_azimuth);
      compute_azimuth_intervals(scratch);
    }

    // Azimuth half of the cull, ahead of the bounds pass. The clip row band is
    // conservative until the exact face bounds are available below.
    if (clip && clip->render_y_start() < clip->render_y_end() &&
        !pole_within_circumcircle() &&
        clip_rejects_azimuth(*clip, clip->render_y_start(),
                             clip->render_y_end() - 1)) {
      ++scratch.claim_seq;
      count = 0;
      y_min = 1;
      y_max = 0;
      return;
    }

    {
      HS_PROFILE_DEEP(face_bounds);

      // Vertical bounds via full arc-extrema + pole analysis. A vertex-only phi
      // span misses the great-circle edge bulge toward a pole, leaving
      // near-pole faces with unscanned rows; the arc-extrema path covers them.
      compute_full_bounds(scratch, count, center, h_virt, height, y_min, y_max);
      compute_inradius(scratch);
    }

    edge_vectors = std::span<Vector>(scratch.edge_vectors.data(), count);
    edge_lengths_sq = std::span<float>(scratch.edge_lengths_sq.data(), count);
    inv_edge_lengths_sq =
        std::span<float>(scratch.inv_edge_lengths_sq.data(), count);
    inv_edge_j = std::span<float>(scratch.inv_edge_j.data(), count);

    {
      HS_PROFILE_DEEP(face_pole);
      apply_pole_containment(height);
    }

    // Whole-face clip cull: y_min/y_max and the azimuth coverage now match what
    // the scan would draw, so a face disjoint from the clip band yields no
    // in-band pixel.
    if (clip && clip_rejects(*clip)) {
      ++scratch.claim_seq;
      count = 0;
      y_min = 1;
      y_max = 0;
      return;
    }

    {
      HS_PROFILE_DEEP(face_edges);
      pack_edges(scratch);
      build_half_planes(scratch, area2);
    }
    {
      HS_PROFILE_DEEP(face_sectors);
      build_sectors(scratch);
    }

    scratch_owner = &scratch;
    scratch_claim = ++scratch.claim_seq;
  }

  /**
   * @brief Tests whether the clip band excludes the whole face.
   * @param cr Clip region (display bounds plus render margin).
   * @return True when no in-band pixel can be produced — the face's vertical
   *         band lies outside the render rows, or its azimuth coverage lies
   *         outside the render columns. A full-width face or an inactive x-clip
   *         never rejects on the horizontal axis.
   * @details Mirrors the scan's own culls (Scan::rasterize's vertical clamp and
   *          per-fragment XClip), so a reject here drops only faces the scan
   *          would have shaded nothing for.
   */
  bool clip_rejects(const ClipRegion &cr) const {
    if (y_max < cr.render_y_start() || y_min > cr.render_y_end() - 1)
      return true;
    return clip_rejects_azimuth(cr, std::max(y_min, cr.render_y_start()),
                                std::min(y_max, cr.render_y_end() - 1));
  }

  /**
   * @brief Horizontal half of the clip cull.
   * @param cr Clip region (display bounds plus render margin).
   * @param band_y_min First row whose azimuth coverage is relevant.
   * @param band_y_max Last row whose azimuth coverage is relevant.
   * @return True when the face's azimuth coverage lies outside the render
   *         columns. A full-width face or an inactive x-clip never rejects.
   */
  bool clip_rejects_azimuth(const ClipRegion &cr, int band_y_min,
                            int band_y_max) const {
    if (full_width)
      return false;
    const ClipRegion::XClip xc = cr.x_clip();
    if (!xc.active)
      return false;
    const int Wd = cr.w;
    const int h_virt = cr.h + hs::H_OFFSET;
    const float phi_scale = PI_F / static_cast<float>(h_virt - 1);
    const float sin_phi =
        std::min(sinf(band_y_min * phi_scale), sinf(band_y_max * phi_scale));
    const float pw = face_azimuth_pad(Wd, sin_phi);
    const int band_len = xc.length(Wd);
    for (const auto &iv : intervals) {
      // Mirrors get_horizontal_intervals' radians->column mapping, so the cull
      // matches the emitted columns exactly.
      int a = static_cast<int>(floorf((iv.first - pw) * Wd / TWO_PI_F));
      int b = static_cast<int>(ceilf((iv.second + pw) * Wd / TWO_PI_F));
      int len = b - a;
      if (len <= 0)
        continue;
      if (len > Wd)
        len = Wd;
      const int s = ((a % Wd) + Wd) % Wd;
      if (ClipRegion::arcs_overlap(xc.rs, band_len, s, len, Wd))
        return false;
    }
    return true;
  }

  // ---------------------------------------------------------------------------
  // Construction phases: single-call-site helpers factored out of the ctor,
  // force-inlined.
  // ---------------------------------------------------------------------------

  /**
   * @brief Latitude-band reject for the face.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices.
   * @param h_virt Virtual row count (height plus pole offset).
   * @param height Canvas height in rows.
   * @return True when the phi extent plus AA margin maps to an empty
   *         canvas-row range.
   */
  __attribute__((always_inline)) bool
  compute_phi_extent(std::span<const Vector> vertices,
                     std::span<const uint16_t> indices, int h_virt,
                     int height) const {
    float min_y_val = 2.0f;
    float max_y_val = -2.0f;

    for (int idx : indices) {
      float y = vertices[idx].y;
      min_y_val = __builtin_fminf(y, min_y_val);
      max_y_val = __builtin_fmaxf(y, max_y_val);
    }

    float min_phi_check = fast_acos(hs::clamp(max_y_val, -1.0f, 1.0f));
    float max_phi_check = fast_acos(hs::clamp(min_y_val, -1.0f, 1.0f));
    Bounds rows =
        phi_bounds_to_rows(min_phi_check - BOUNDS_MARGIN,
                           max_phi_check + BOUNDS_MARGIN, h_virt, height);

    return rows.y_min > rows.y_max;
  }

  /**
   * @brief Builds the local tangent frame, gnomonic 2D projection, and 3D
   * arrays.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices.
   * @param scratch Scratch storage receiving poly_2d and verts_3d.
   * @details Sets basis_u/v/w, poly_2d (with circumradius), and the local 3D
   * vertex array.
   */
  __attribute__((always_inline)) void
  setup_frame_and_polygon(std::span<const Vector> vertices,
                          std::span<const uint16_t> indices,
                          FaceScratchBuffer &scratch) {
    // Single gather over the shared pool: every later pass reads the local
    // copy instead of chasing indices back into `vertices`.
    center = Vector(0, 0, 0);
    for (int i = 0; i < count; ++i) {
      const Vector &v = vertices[indices[i]];
      scratch.verts_3d[i] = v;
      center = center + v;
    }
    center.normalize();

    basis_v = center;
    Vector ref = least_parallel_axis(center);
    basis_u = cross(center, ref).normalized();
    basis_w = cross(center, basis_u).normalized();

    float max_r2 = 0.0f;
    for (int i = 0; i < count; ++i) {
      const Vector &v = scratch.verts_3d[i];
      // Gnomonic projection divides by d = cos(angle from face center),
      // singular near the center's antipode; clamp d away from zero,
      // sign-preserving.
      float d = dot(v, basis_v);
      if (fabsf(d) < math::TOLERANCE)
        d = copysignf(math::TOLERANCE, d);
      float px = dot(v, basis_u) / d;
      float py = dot(v, basis_w) / d;

      scratch.poly_2d[i] = Vector(px, py, 0);

      float r2 = px * px + py * py;
      max_r2 = __builtin_fmaxf(r2, max_r2);
    }
    radius = sqrtf(max_r2);
    max_dist = radius + BOUNDS_MARGIN_WIDE;
    max_dist_sq = max_dist * max_dist;

    scratch.poly_2d[count] = scratch.poly_2d[0];
    poly_2d = std::span<Vector>(scratch.poly_2d.data(), count + 1);

    scratch.verts_3d[count] = scratch.verts_3d[0];
  }

  /**
   * @brief Computes the face "size" (inradius) from the projected polygon.
   * @param scratch Scratch storage holding poly_2d and the per-edge vectors and
   * squared lengths compute_full_bounds stored, which must run first.
   * @details size = minimum distance from the projected centroid to any edge,
   * floored to a fraction of the circumradius for degenerate slivers. A large
   * face converts it to radians, matching the metric distance() reports.
   */
  __attribute__((always_inline)) void
  compute_inradius(const FaceScratchBuffer &scratch) {
    float min_edge_dist = 1e9f;
    for (int i = 0; i < count; ++i) {
      const Vector &v1 = scratch.poly_2d[i];
      const Vector &edge = scratch.edge_vectors[i];
      float t = 0.0f;
      const float inv_edge_len_sq = scratch.inv_edge_lengths_sq[i];
      if (inv_edge_len_sq > 0.0f) {
        t = dot(-v1, edge) * inv_edge_len_sq;
        t = __builtin_fmaxf(0.0f, __builtin_fminf(1.0f, t));
      }
      Vector closest = v1 + edge * t;
      float d_line = closest.magnitude();

      min_edge_dist = __builtin_fminf(d_line, min_edge_dist);
    }
    size = __builtin_fmaxf(min_edge_dist, radius * MIN_SIZE_RADIUS_RATIO);
    linear_dist = size < 0.2f;
    // distance() reports radians for a large face; the same fast_atan2 keeps
    // dist/size exactly 1 at the inradius.
    if (!linear_dist)
      size = fast_atan2(size, 1.0f);
  }

  /**
   * @brief Packs per-edge data contiguously for the distance() fallback.
   * @param scratch Scratch storage receiving packed_edges.
   */
  __attribute__((always_inline)) void pack_edges(FaceScratchBuffer &scratch) {
    for (int i = 0; i < count; ++i) {
      auto &ep = scratch.packed_edges[i];
      ep.vx = poly_2d[i].x;
      ep.vy = poly_2d[i].y;
      ep.ex = edge_vectors[i].x;
      ep.ey = edge_vectors[i].y;
      ep.inv_len_sq = inv_edge_lengths_sq[i];
      ep.inv_ej = inv_edge_j[i];
      ep.key_vy = angle_key(poly_2d[i].y);
      // A y-degenerate edge (inv_ej == 0) has no usable crossing x; equal keys
      // drop it from distance()'s parity test.
      ep.key_next_vy =
          (inv_edge_j[i] != 0.0f) ? angle_key(poly_2d[i + 1].y) : ep.key_vy;
    }
    packed_edges = std::span<EdgePacked>(scratch.packed_edges.data(), count);
  }

  /**
   * @brief Detects a convex 2D projection and builds its edge half-planes.
   * @param scratch Scratch storage receiving half_planes.
   * @param area2 Twice the polygon's signed area, from the collapsed-face cull.
   * @details For a convex polygon the signed distance is max over edges of the
   * half-plane distance: exact everywhere inside and in the edge slabs outside;
   * outside a vertex's normal cone it underestimates (line distance, not vertex
   * distance), which only softens the AA corner within its ~1-pixel band. A
   * concave, degenerate-edged, or wrongly-oriented polygon leaves convex false
   * and distance() on the exact walk.
   */
  __attribute__((always_inline)) void
  build_half_planes(FaceScratchBuffer &scratch, float area2) {
    bool pos = false, neg = false;
    // The turn test carries the previous edge in registers, so the ring closes
    // without indexing edge_vectors through a modulo.
    const Vector *e1 = &edge_vectors[count - 1];
    float l1 = edge_lengths_sq[count - 1];
    for (int i = 0; i < count; ++i) {
      const Vector &e2 = edge_vectors[i];
      float cr = e1->x * e2.y - e1->y * e2.x;
      float scale = l1 * edge_lengths_sq[i];
      e1 = &e2;
      l1 = edge_lengths_sq[i];
      if (cr * cr > TURN_EPS_SQ * scale) {
        if (cr > 0)
          pos = true;
        else
          neg = true;
      }
    }
    if (pos && neg)
      return;

    float sign = area2 >= 0 ? 1.0f : -1.0f;
    // d0 is the origin's worst half-plane distance: the projected face center
    // must be strictly interior or the winding/orientation is untrustworthy.
    float d0 = -FLT_MAX;
    for (int i = 0; i < count; ++i) {
      float len_sq = edge_lengths_sq[i];
      if (len_sq < 1e-12f)
        return;
      float inv = sign / sqrtf(len_sq);
      float nx = edge_vectors[i].y * inv;
      float ny = -edge_vectors[i].x * inv;
      float off = -(nx * poly_2d[i].x + ny * poly_2d[i].y);
      auto &hp = scratch.half_planes[i];
      hp.nx = nx;
      hp.ny = ny;
      hp.off = off;
      d0 = __builtin_fmaxf(off, d0);
    }
    if (d0 >= 0.0f)
      return;
    half_planes = std::span<const HalfPlane>(scratch.half_planes.data(), count);
    convex = true;
  }

  /**
   * @brief Builds the angular sector table for the concave sector walk.
   * @param scratch Scratch storage receiving the unwrapped vertex
   * pseudo-angles.
   * @details Only concave faces with at least SECTOR_MIN_COUNT vertices
   * qualify. A face that is star-shaped about its projected centroid (the
   * gnomonic origin) has monotonic vertex pseudo-angles spanning a full turn;
   * that monotonicity is what lets plane_dist_sector bin a query point into one
   * fan sector by angle alone. Strictly monotonic faces bin exactly (K1);
   * mildly-bent faces whose worst vertex backtracks by no more than
   * SECTOR_MONO_TOL still bin to within a neighbor and take the wider K2 walk.
   * A larger inversion or wrong total turn (not star-shaped, e.g. heavily
   * deformed) leaves sector_ok false and keeps the exact walk.
   */
  __attribute__((always_inline)) void
  build_sectors(FaceScratchBuffer &scratch) {
    sector_ok = false;
    if (convex || count < SECTOR_MIN_COUNT)
      return;
    float prev = pseudo_angle(poly_2d[0].y, poly_2d[0].x);
    float acc = prev, total = 0.0f;
    scratch.pseudo_angles[0] = prev;
    for (int i = 1; i <= count; ++i) {
      float a = pseudo_angle(poly_2d[i].y, poly_2d[i].x);
      float d = a - prev;
      if (d > 2.0f)
        d -= 4.0f;
      if (d < -2.0f)
        d += 4.0f;
      acc += d;
      total += d;
      scratch.pseudo_angles[i] = acc;
      prev = a;
    }
    // Star-shaped about the centroid <=> exactly one full turn with no vertex
    // backtracking. Reject a wrong total turn: the sectors would overlap or
    // leave a gap.
    if (fabsf(fabsf(total) - 4.0f) > 1e-3f)
      return;
    // Fold the winding into the table: pseudo_angles becomes weakly increasing
    // whichever way the polygon is wound, so the probe search compares one
    // direction. Scaling by +/-1 is exact, so every step and bin is unchanged.
    float sgn = (total >= 0.0f) ? 1.0f : -1.0f;
    float min_step = FLT_MAX;
    float prev_s = scratch.pseudo_angles[0] * sgn;
    scratch.pseudo_angles[0] = prev_s;
    for (int i = 1; i <= count; ++i) {
      float cur = scratch.pseudo_angles[i] * sgn;
      scratch.pseudo_angles[i] = cur;
      min_step = __builtin_fminf(cur - prev_s, min_step);
      prev_s = cur;
    }
    // min_step > 0: strictly monotonic, sectors don't overlap, K1 bins exactly.
    // A backtrack up to SECTOR_MONO_TOL overlaps sectors slightly, so the bin
    // can land one neighbor off -> K2's wider walk still reaches the true edge.
    // A larger inversion is not star-shaped; keep the exact walk.
    if (min_step <= -SECTOR_MONO_TOL)
      return;
    sector_kmax = (min_step > 0.0f) ? 1 : SECTOR_KMAX_MAX;
    for (int i = 0; i <= count; ++i)
      scratch.sector_keys[i] = angle_key(scratch.pseudo_angles[i]);
    sector_keys =
        std::span<const uint32_t>(scratch.sector_keys.data(), count + 1);
    sector_base = scratch.pseudo_angles[0];
    sector_span = total * sgn;
    sector_sgn = sgn;
    sector_ok = true;
  }

  /**
   * @brief Aligns the current projection to a canonical class shape and binds
   *        its distance LUT.
   * @param lut Canonical-frame LUT for the face's congruence class.
   * @param canon_xy Canonical centered 2D polygon, x/y pairs.
   * @param vert_offset Cyclic offset aligning mesh vertex order to canonical.
   * @param reflected True for the mirror family.
   * @return False when the correlation is degenerate (badly deformed face);
   *         the face then keeps the exact path.
   * @details One complex correlation over the vertices recovers the in-plane
   * rotation placing the canonical shape at the face's least-squares pose, so
   * the interior gradient follows the face's rigid motion through a ripple
   * while edges stay exact. Rotational-symmetry ambiguity is harmless (the LUT
   * is invariant under the shape's symmetry group).
   */
  HS_COLD_MEMBER bool bind_class_lut(const ClassLut *lut, const float *canon_xy,
                                     int vert_offset, bool reflected) {
    HS_CHECK(vert_offset >= 0 && vert_offset < count,
             "bind_class_lut: vertex offset outside the face");
    float mx = 0.0f, my = 0.0f;
    for (int i = 0; i < count; ++i) {
      mx += poly_2d[i].x;
      my += poly_2d[i].y;
    }
    float inv_n = 1.0f / count;
    mx *= inv_n;
    my *= inv_n;

    AlignCorr a = align_correlate(canon_xy, count, vert_offset, reflected,
                                  [&](int j, float &zx, float &zy) {
                                    zx = poly_2d[j].x - mx;
                                    zy = poly_2d[j].y - my;
                                  });
    float r2 = a.rr * a.rr + a.ri * a.ri;
    if (r2 <= ALIGN_MIN_CORR_SQ * a.cc * a.zz)
      return false;
    float inv_r = 1.0f / sqrtf(r2);
    float c = a.rr * inv_r, s = a.ri * inv_r;

    // |d_true - d_canon| is bounded by the worst aligned vertex deviation, so
    // widen the sign-purity guard by that bound: a bent face falls back to the
    // exact walk near its true edges instead of serving wrong-signed canonical
    // distances (faces separating under ripple). Faces bent beyond range keep
    // the exact path.
    float max_dev_sq = 0.0f;
    align_walk(
        count, vert_offset, reflected,
        [&](int j, float &zx, float &zy) {
          zx = poly_2d[j].x - mx;
          zy = poly_2d[j].y - my;
        },
        [&](int k, float zx, float zy) {
          float ex = canon_xy[2 * k] - (c * zx - s * zy);
          float ey = canon_xy[2 * k + 1] - (s * zx + c * zy);
          float dev_sq = ex * ex + ey * ey;
          if (dev_sq > max_dev_sq)
            max_dev_sq = dev_sq;
        });
    float max_dev = sqrtf(max_dev_sq);
    if (max_dev > ALIGN_MAX_DEV_DIAGS * lut->safe_dist)
      return false;

    // q = rot * (p - m), with the mirror family's conjugation folded into the
    // matrix; then the LUT grid transform (q - box_min) * inv_step folded on
    // top, so the probe loop runs a single affine map on raw (px, py).
    float m00 = c, m01 = reflected ? s : -s;
    float m10 = s, m11 = reflected ? -c : c;
    lut_ax = m00 * lut->inv_step_x;
    lut_bx = m01 * lut->inv_step_x;
    lut_cx = (-(m00 * mx + m01 * my) - lut->cx + lut->Rx) * lut->inv_step_x;
    lut_ay = m10 * lut->inv_step_y;
    lut_by = m11 * lut->inv_step_y;
    lut_cy = (-(m10 * mx + m11 * my) - lut->cy + lut->Ry) * lut->inv_step_y;
    lut_n = lut->n;
    lut_clamp = static_cast<float>(lut->n - 2);
    lut_dequant = lut->dequant;
    lut_q_safe =
        static_cast<int32_t>((lut->safe_dist + max_dev) / lut->dequant);
    lut_data = lut->data;
    return true;
  }

  /**
   * @brief Fills the scratch vertex azimuths the interval pass consumes.
   * @param scratch Scratch storage holding verts_3d and receiving thetas.
   */
  __attribute__((always_inline)) void
  compute_thetas(FaceScratchBuffer &scratch) const {
    for (int i = 0; i < count; ++i) {
      const Vector &v = scratch.verts_3d[i];
      float theta = fast_atan2(v.z, v.x);
      if (theta < 0)
        theta += TWO_PI_F;
      scratch.thetas[i] = theta;
    }
  }

  /**
   * @brief Computes the face's azimuth coverage intervals.
   * @param scratch Scratch storage holding thetas and receiving intervals.
   * @details Finds the largest angular gap between vertices; if it exceeds pi
   * the face does not wrap, so the complementary horizontal interval(s) are
   * emitted, else it spans full width. Coarse: only the single largest gap is
   * excised.
   */
  __attribute__((always_inline)) void
  compute_azimuth_intervals(FaceScratchBuffer &scratch) {
    // Insertion sort: faces carry a few dozen vertices at most, and std::sort
    // stays out of line here (an __introsort_loop plus an __insertion_sort
    // call per face).
    float *th = scratch.thetas.data();
    for (int i = 1; i < count; ++i) {
      float t = th[i];
      float *p = th + i;
      while (p != th && p[-1] > t) {
        p[0] = p[-1];
        --p;
      }
      *p = t;
    }
    float max_gap = 0;
    float gap_start = 0;
    for (int i = 0; i < count; ++i) {
      float next = (i + 1 < count) ? scratch.thetas[i + 1]
                                   : (scratch.thetas[0] + TWO_PI_F);
      float diff = next - scratch.thetas[i];
      if (diff > max_gap) {
        max_gap = diff;
        gap_start = scratch.thetas[i];
      }
    }

    int interval_count = 0;
    if (max_gap > PI_F) {
      full_width = false;
      float start_t = fmodf(gap_start + max_gap, TWO_PI_F);
      // fmodf can leave start_t at ~2*PI instead of ~0, producing a degenerate
      // [~2*PI, 2*PI] sliver below; snap to 0.
      if (start_t > TWO_PI_F - 1e-4f)
        start_t = 0.0f;
      float end_t = gap_start;

      if (start_t <= end_t) {
        scratch.intervals[interval_count++] = {start_t, end_t};
      } else {
        scratch.intervals[interval_count++] = {0.0f, end_t};
        scratch.intervals[interval_count++] = {start_t, TWO_PI_F};
      }
    } else {
      full_width = true;
    }
    intervals = std::span<Interval>(scratch.intervals.data(), interval_count);
  }

  /**
   * @brief Necessary condition for apply_pole_containment to fire.
   * @return True when a pole falls within the gnomonic circumcircle of the
   *         face's vertices. Both poles project to the same radius, so one
   *         test covers each. pole_inside_polygon can only report inside for
   *         points within the vertex convex hull, so a false here rules out
   *         pole containment.
   */
  __attribute__((always_inline)) bool pole_within_circumcircle() const {
    return center.y * center.y * (1.0f + radius * radius) >= 1.0f;
  }

  /**
   * @brief Ray-crossing test of a projected pole against the 2D face polygon.
   * @param ppx Projected pole x in the face's 2D basis.
   * @param ppy Projected pole y in the face's 2D basis.
   * @return true if (ppx, ppy) lies inside the polygon.
   */
  __attribute__((always_inline)) bool pole_inside_polygon(float ppx,
                                                          float ppy) const {
    bool inside = false;
    for (int i = 0; i < count; ++i) {
      if (inv_edge_j[i] != 0.0f &&
          (poly_2d[i].y > ppy) != (poly_2d[i + 1].y > ppy)) {
        float ix = poly_2d[i].x +
                   (ppy - poly_2d[i].y) * edge_vectors[i].x * inv_edge_j[i];
        if (ppx < ix)
          inside = !inside;
      }
    }
    return inside;
  }

  /** How a projected pole sits against the 2D face polygon. */
  enum class PoleHit {
    NONE,     /**< Outside the closed region. */
    BOUNDARY, /**< Within POLE_BOUNDARY_TOL of a vertex or an edge. */
    INTERIOR  /**< Strictly inside, off the boundary band. */
  };

  /**
   * @brief Classifies a projected pole against the 2D face polygon.
   * @param ppx Projected pole x in the face's 2D basis.
   * @param ppy Projected pole y in the face's 2D basis.
   * @return Where the pole falls relative to the closed polygon.
   * @details A pole sitting exactly on a vertex or an edge - where a partition
   * op plants an apex on the projection axis - leaves pole_inside_polygon's
   * crossing parity decided by rounding, so congruent faces disagree on their
   * azimuth coverage and the pole row loses part of its overlap. Testing the
   * boundary band first makes the answer independent of that rounding.
   */
  PoleHit pole_hit(float ppx, float ppy) const {
    if (!pole_within_circumcircle())
      return PoleHit::NONE;
    const float tol = radius * POLE_BOUNDARY_TOL;
    const float tol_sq = tol * tol;
    for (int i = 0; i < count; ++i) {
      const float ax = ppx - poly_2d[i].x;
      const float ay = ppy - poly_2d[i].y;
      if (ax * ax + ay * ay <= tol_sq)
        return PoleHit::BOUNDARY;
      const float t =
          hs::clamp((ax * edge_vectors[i].x + ay * edge_vectors[i].y) *
                        inv_edge_lengths_sq[i],
                    0.0f, 1.0f);
      const float dx = ax - t * edge_vectors[i].x;
      const float dy = ay - t * edge_vectors[i].y;
      if (dx * dx + dy * dy <= tol_sq)
        return PoleHit::BOUNDARY;
    }
    return pole_inside_polygon(ppx, ppy) ? PoleHit::INTERIOR : PoleHit::NONE;
  }

  /**
   * @brief Extends the vertical bounds when the face reaches a pole.
   * @param height Canvas height in rows.
   * @details A face enclosing a pole wraps every azimuth, so it takes full
   * width outright. A face that only meets the pole on its boundary keeps its
   * azimuth wedge, which get_horizontal_intervals widens per row.
   */
  __attribute__((always_inline)) void apply_pole_containment(int height) {
    if (center.y > 0.01f) {
      float inv_c = 1.0f / center.y;
      const PoleHit hit = pole_hit(basis_u.y * inv_c, basis_w.y * inv_c);
      if (hit != PoleHit::NONE) {
        y_min = 0;
        if (hit == PoleHit::INTERIOR)
          full_width = true;
      }
    }
    // South pole (0, -1, 0)
    if (center.y < -0.01f) {
      float inv_c = 1.0f / -center.y;
      const PoleHit hit = pole_hit(-basis_u.y * inv_c, -basis_w.y * inv_c);
      if (hit != PoleHit::NONE) {
        y_max = height - 1;
        if (hit == PoleHit::INTERIOR)
          full_width = true;
      }
    }
  }

  /**
   * @brief Refines phi bounds with an edge's great-circle arc extremum.
   * @param n Normalized great-circle plane normal of the edge.
   * @param v1 First edge endpoint (on the unit sphere).
   * @param v2 Second edge endpoint (on the unit sphere).
   * @param min_phi In/out running minimum phi.
   * @param max_phi In/out running maximum phi.
   * @details The extremum of an edge's arc may lie between its endpoints;
   * project the pole-tangent onto the plane and, if it falls inside the arc,
   * fold its phi into the bounds.
   */
  static __attribute__((always_inline)) void
  refine_phi_from_arc_extremum(const Vector &n, const Vector &v1,
                               const Vector &v2, float &min_phi,
                               float &max_phi) {
    float ny = n.y;
    if (std::abs(ny) < 0.99999f) {
      float nx = n.x;
      float nz = n.z;
      float tx = -nx * ny;
      float ty = 1.0f - ny * ny;
      float tz = -nz * ny;
      float t_len_sq = tx * tx + ty * ty + tz * tz;
      if (t_len_sq > 1e-12f) {
        float inv_len = 1.0f / sqrtf(t_len_sq);
        float ptx = tx * inv_len;
        float pty = ty * inv_len;
        float ptz = tz * inv_len;
        float cx1 = (v1.y * ptz - v1.z * pty) * nx +
                    (v1.z * ptx - v1.x * ptz) * ny +
                    (v1.x * pty - v1.y * ptx) * nz;
        float cx2 = (pty * v2.z - ptz * v2.y) * nx +
                    (ptz * v2.x - ptx * v2.z) * ny +
                    (ptx * v2.y - pty * v2.x) * nz;
        if (cx1 > 0 && cx2 > 0)
          min_phi =
              __builtin_fminf(fast_acos(hs::clamp(pty, -1.0f, 1.0f)), min_phi);
        if (cx1 < 0 && cx2 < 0)
          max_phi =
              __builtin_fmaxf(fast_acos(hs::clamp(-pty, -1.0f, 1.0f)), max_phi);
      }
    }
  }

  /**
   * @brief Snaps phi bounds to a pole when the face's planes enclose it.
   * @param scratch Scratch storage holding the compacted great-circle planes.
   * @param planes_count Number of valid planes.
   * @param center Normalized face centroid.
   * @param min_phi In/out phi minimum, set to 0 if the north pole is enclosed.
   * @param max_phi In/out phi maximum, set to PI if the south pole is enclosed.
   */
  static __attribute__((always_inline)) void
  snap_phi_for_pole_planes(const FaceScratchBuffer &scratch, int planes_count,
                           const Vector &center, float &min_phi,
                           float &max_phi) {
    bool np_inside = (planes_count > 0);
    bool sp_inside = (planes_count > 0);
    // Both flags only ever clear, so the scan is done once neither survives.
    for (int pi = 0; pi < planes_count && (np_inside || sp_inside); ++pi) {
      float py = scratch.planes[pi].y;
      bool center_pos = dot(center, scratch.planes[pi]) > 0;
      if ((py > 0) != center_pos)
        np_inside = false;
      if ((py < 0) != center_pos)
        sp_inside = false;
    }
    if (np_inside)
      min_phi = 0.0f;
    if (sp_inside)
      max_phi = PI_F;
  }

  /**
   * @brief Full-path vertical bounds (arc extrema + pole containment) for large
   * faces.
   * @param scratch Scratch storage holding poly_2d/verts_3d, receiving edge
   * data and planes.
   * @param count Vertex/edge count.
   * @param center Normalized face centroid.
   * @param h_virt Virtual row count (height plus pole offset).
   * @param height Canvas height in rows.
   * @param y_min_out Output: first covered row.
   * @param y_max_out Output: last covered row.
   */
  HS_O3_FN static void compute_full_bounds(FaceScratchBuffer &scratch,
                                           int count, const Vector &center,
                                           int h_virt, int height,
                                           int &y_min_out, int &y_max_out) {
    float min_phi = 100.0f;
    float max_phi = -100.0f;
    int planes_count = 0;
    for (int i = 0; i < count; ++i) {
      const Vector &v1 = scratch.verts_3d[i];
      const Vector &v2 = scratch.verts_3d[i + 1];
      Vector edge = scratch.poly_2d[i + 1] - scratch.poly_2d[i];
      scratch.edge_vectors[i] = edge;
      float edge_len_sq = dot(edge, edge);
      scratch.edge_lengths_sq[i] = edge_len_sq;
      scratch.inv_edge_lengths_sq[i] =
          (edge_len_sq > 1e-12f) ? (1.0f / edge_len_sq) : 0.0f;
      scratch.inv_edge_j[i] =
          (std::abs(edge.y) > 1e-12f) ? (1.0f / edge.y) : 0.0f;
      Vector normal = cross(v1, v2);
      float len_sq = dot(normal, normal);
      // planes[] is COMPACTED: a degenerate edge pushes no plane, so planes[k]
      // does NOT correspond to edge k (unlike the per-edge arrays indexed by
      // i). Downstream consumers treat planes[] as a standalone set, never by
      // edge.
      if (len_sq > 1e-12f)
        scratch.planes[planes_count++] = normal.normalized();
      float phi_val = fast_acos(hs::clamp(v1.y, -1.0f, 1.0f));
      min_phi = __builtin_fminf(phi_val, min_phi);
      max_phi = __builtin_fmaxf(phi_val, max_phi);
      // Arc Extrema Logic: only when this edge pushed its own plane, else
      // planes[planes_count - 1] is a prior edge's normal against these
      // endpoints.
      if (len_sq > 1e-12f)
        refine_phi_from_arc_extremum(scratch.planes[planes_count - 1], v1, v2,
                                     min_phi, max_phi);
    }
    snap_phi_for_pole_planes(scratch, planes_count, center, min_phi, max_phi);
    Bounds rows = phi_bounds_to_rows(min_phi - BOUNDS_MARGIN,
                                     max_phi + BOUNDS_MARGIN, h_virt, height);
    y_min_out = rows.y_min;
    y_max_out = rows.y_max;
  }

  /**
   * @brief Returns the face's precomputed inclusive row bounds.
   * @tparam H Canvas height in rows; must match the construction height the
   * bounds were computed for.
   * @return The stored {y_min, y_max} bounds.
   */
  template <int H> Bounds get_vertical_bounds() const {
    HS_CHECK(H == build_height,
             "Face::get_vertical_bounds: H differs from construction height");
    HS_CHECK(!scratch_owner || scratch_owner->claim_seq == scratch_claim,
             "SDF::Face scanned after a later Face claimed its scratch buffer");
    return {y_min, y_max};
  }

  /**
   * @brief Reports whether rounded azimuth intervals change across a row band.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @param y_lo First row in the band.
   * @param y_hi Last row in the band.
   * @return True when the band requires per-row interval construction.
   */
  template <int W, int H>
  bool horizontal_intervals_vary_by_row(int y_lo, int y_hi) const {
    if (full_width || y_lo >= y_hi)
      return false;
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();

    const float sin_lo = TrigLUT<W, H>::sin_phi[y_lo];
    const float sin_hi = TrigLUT<W, H>::sin_phi[y_hi];
    const float sin_min = std::min(sin_lo, sin_hi);
    float sin_max = std::max(sin_lo, sin_hi);
    constexpr int H_VIRT = H + hs::H_OFFSET;
    constexpr int EQUATOR_LO = (H_VIRT - 1) / 2;
    constexpr int EQUATOR_HI = H_VIRT / 2;
    if (y_lo <= EQUATOR_LO && EQUATOR_LO <= y_hi)
      sin_max = std::max(sin_max, TrigLUT<W, H>::sin_phi[EQUATOR_LO]);
    if (y_lo <= EQUATOR_HI && EQUATOR_HI <= y_hi)
      sin_max = std::max(sin_max, TrigLUT<W, H>::sin_phi[EQUATOR_HI]);

    const float narrow_pad = face_azimuth_pad(W, sin_max);
    const float wide_pad = face_azimuth_pad(W, sin_min);
    const float column_scale = W / TWO_PI_F;
    for (const auto &iv : intervals) {
      if (floorf((iv.first - narrow_pad) * column_scale) !=
              floorf((iv.first - wide_pad) * column_scale) ||
          ceilf((iv.second + narrow_pad) * column_scale) !=
              ceilf((iv.second + wide_pad) * column_scale))
        return true;
    }
    return false;
  }

  /**
   * @brief Emits the face's azimuth-coverage intervals for a row.
   * @tparam W Canvas width in columns; must match the clip width the
   * construction-time azimuth cull ran against.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y Row index, which sets the pole widening below.
   * @param out Sink accepting (float start, float end).
   * @return True if intervals were emitted; false when the row requires a full
   * scan.
   * @details The pad is an azimuth angle, so it holds a whole pixel of AA reach
   * only at the equator; at colatitude phi one pad p of great-circle reach
   * subtends asin(p / sin phi), reaching the whole row once sin phi <= p.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    HS_CHECK(
        build_width == 0 || W == build_width,
        "Face::get_horizontal_intervals: W differs from the clip width the "
        "azimuth cull ran against");
    if (full_width)
      return false;
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    const float sin_phi = TrigLUT<W, H>::sin_phi[y];
    const float base_pad = face_azimuth_pad(W);
    if (sin_phi <= base_pad)
      return false;
    const float pad = asinf(base_pad / sin_phi);
    for (const auto &iv : intervals) {
      float f_x1 = (iv.first - pad) * W / TWO_PI_F;
      float f_x2 = (iv.second + pad) * W / TWO_PI_F;
      out(floorf(f_x1), ceilf(f_x2));
    }
    return true;
  }

  /**
   * @brief Signed planar distance via the convex half-plane max.
   * @param px Gnomonic x of the query point.
   * @param py Gnomonic y of the query point.
   * @return Signed distance in the tangent plane (negative inside).
   */
  HS_O3_FN float plane_dist_convex(float px, float py) const {
    HS_SCAN_METRIC(hs::g_scan_metrics.convex_hits++);
    float d = -FLT_MAX;
    for (int i = 0; i < count; ++i) {
      const auto &hp = half_planes[i];
      float di = hp.nx * px + hp.ny * py + hp.off;
      d = __builtin_fmaxf(di, d);
    }
    return d;
  }

  /**
   * @brief Squared planar distance via the exact per-edge walk.
   * @param px Gnomonic x of the query point.
   * @param py Gnomonic y of the query point.
   * @param inside_out Set true when the query lies inside the polygon; carries
   * the sign the squared return cannot.
   * @return Squared distance to the nearest edge, in the tangent plane.
   */
  HS_O3_FN float plane_dsq_exact(float px, float py, bool &inside_out) const {
    float d = FLT_MAX;
    bool inside = false;
    const uint32_t qk = angle_key(py);
    for (int i = 0; i < count; ++i) {
      const auto &ep = packed_edges[i];
      float wx = px - ep.vx, wy = py - ep.vy;
      float t = (wx * ep.ex + wy * ep.ey) * ep.inv_len_sq;
      float cv = hs::clamp(t, 0.0f, 1.0f);
      float bx = wx - ep.ex * cv, by = wy - ep.ey * cv;
      float dsq = bx * bx + by * by;
      d = __builtin_fminf(dsq, d);
      if ((ep.key_vy > qk) != (ep.key_next_vy > qk)) {
        float isx = ep.vx + (py - ep.vy) * ep.ex * ep.inv_ej;
        if (px < isx)
          inside = !inside;
      }
    }
    inside_out = inside;
    return d;
  }

  /**
   * @brief Squared planar distance via the concave sector walk.
   * @param px Gnomonic x of the query point.
   * @param py Gnomonic y of the query point.
   * @param inside_out Set true when the query lies inside the polygon; carries
   * the sign the squared return cannot.
   * @return Squared distance to the nearest edge, in the tangent plane.
   * @details Bins the query into its fan sector by pseudo-angle (a binary
   * search over the monotonic vertex angle_keys), then takes the exact min
   * segment distance over only that sector's edge and its sector_kmax neighbors
   * each side (K1 = 1 for strict faces, K2 = 2 for mildly-bent faces whose bin
   * can land a neighbor off). The sign comes free from the sector's boundary
   * edge: the polygon is consistently wound, so a query on the interior side of
   * edge s (matched to the winding) is inside. Near-exact for star faces because
   * the true nearest edge is almost always the sector's own edge or an immediate
   * neighbor. Only enabled when build_sectors set sector_ok.
   */
  HS_O3_FN float plane_dsq_sector(float px, float py, bool &inside_out) const {
    float p = pseudo_angle(py, px) * sector_sgn;
    float rel = (p - sector_base) / sector_span; // -> [0, 1) after the fold
    rel -= floorf(rel);
    uint32_t qk = angle_key(sector_base + rel * sector_span);
    int lo = 0, hi = count;
    while (lo + 1 < hi) {
      int mid = (lo + hi) >> 1;
      if (sector_keys[mid] <= qk)
        lo = mid;
      else
        hi = mid;
    }
    int s = lo;

    float d = FLT_MAX;
    // kmax < SECTOR_MIN_COUNT <= count, so one wrap correction suffices.
    int idx = s - sector_kmax;
    if (idx < 0)
      idx += count;
    for (int k = -sector_kmax; k <= sector_kmax; ++k) {
      const auto &ep = packed_edges[idx];
      if (++idx == count)
        idx = 0;
      float wx = px - ep.vx, wy = py - ep.vy;
      float t =
          hs::clamp((wx * ep.ex + wy * ep.ey) * ep.inv_len_sq, 0.0f, 1.0f);
      float bx = wx - ep.ex * t, by = wy - ep.ey * t;
      float dsq = bx * bx + by * by;
      d = __builtin_fminf(dsq, d);
    }
    const auto &e0 = packed_edges[s];
    float cr = e0.ex * (py - e0.vy) - e0.ey * (px - e0.vx);
    inside_out = cr * sector_sgn >= 0.0f;
    return d;
  }

  static constexpr uint32_t PROBE_HAS_LUT = 1u << 0;
  static constexpr uint32_t PROBE_CONVEX = 1u << 1;
  static constexpr uint32_t PROBE_SECTOR = 1u << 2;
  static constexpr uint32_t PROBE_LINEAR = 1u << 3;

  /** @return Packed distance-path flags for repeated probes of this face. */
  uint32_t probe_flags() const {
    return (lut_data ? PROBE_HAS_LUT : 0u) | (convex ? PROBE_CONVEX : 0u) |
           (sector_ok ? PROBE_SECTOR : 0u) | (linear_dist ? PROBE_LINEAR : 0u);
  }

  /**
   * @brief Computes signed distance to the face, writing into res.
   * @tparam ComputeUVs Accepted for interface parity; the face stores no UVs.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = raw_dist = the signed edge distance,
   *        size = inradius, both in gnomonic plane units on a linear_dist
   *        face and in radians otherwise.
   * @param reject_dsq Squared plane distance at/above which an outside probe
   *        is reported as the far sentinel without taking the square root.
   *        Must be conservative: only probes the caller rejects on dist may
   *        cross it. FLT_MAX disables the cull. Consulted only on a
   *        linear_dist face's edge-walk path; large faces and the convex
   *        half-plane path ignore it.
   * @note Distances live in the face's gnomonic tangent plane. Small faces
   *       (linear_dist) report the plane distance directly; large faces convert
   *       via fast_atan2(plane, 1). Do not treat raw_dist as a metric geodesic
   *       angle.
   */
  template <bool ComputeUVs = true>
  HS_O3_FN void distance(const Vector &p, DistanceResult &res,
                         float reject_dsq = FLT_MAX) const {
    distance_with_flags<ComputeUVs>(p, res, reject_dsq, probe_flags());
  }

  /**
   * @brief Computes distance using flags captured by probe_flags().
   * @tparam ComputeUVs Accepted for interface parity; the face stores no UVs.
   * @param p Point on sphere (normalized).
   * @param res Output distance result.
   * @param reject_dsq Conservative squared distance rejection threshold,
   *        honored only when PROBE_LINEAR is set and the probe lands outside;
   *        the PROBE_CONVEX path ignores it.
   * @param probe_flags Distance-path flags captured after the face's last LUT
   *        binding or geometry update.
   */
  template <bool ComputeUVs = true>
  HS_O3_FN void distance_with_flags(const Vector &p, DistanceResult &res,
                                    float reject_dsq,
                                    uint32_t probe_flags) const {
    HS_SCAN_METRIC(hs::g_scan_metrics.pixels_tested++);
    HS_PROBE_TICK();
    HS_PROBE_COUNT(n_probe);
    HS_PROBE_MARK(hs_t);

    float cos_angle = dot(p, center);
    if (cos_angle <= 0.01f) {
      HS_SCAN_METRIC(hs::g_scan_metrics.pixels_culled++);
      HS_PROBE_SPAN(point, hs_t);
      HS_PROBE_COUNT(n_cull_cos);
      res = DistanceResult(FAR_SENTINEL, 0.0f, FAR_SENTINEL, 0.0f, size);
      return;
    }
    HS_PROBE_SPAN(point, hs_t);

    float inv_cos = 1.0f / cos_angle;
    float px = dot(p, basis_u) * inv_cos;
    float py = dot(p, basis_w) * inv_cos;

    float p_r2 = px * px + py * py;
    if (p_r2 > max_dist_sq) {
      HS_SCAN_METRIC(hs::g_scan_metrics.pixels_culled++);
      HS_PROBE_SPAN(project, hs_t);
      HS_PROBE_COUNT(n_cull_r);
      res = DistanceResult(FAR_SENTINEL, 0.0f, FAR_SENTINEL, 0.0f, size);
      return;
    }
    HS_PROBE_SPAN(project, hs_t);

    float plane_dist;
    bool lut_served = false;
    if (probe_flags & PROBE_HAS_LUT) {
      // Affine map into the canonical LUT grid, then a 4-tap bilinear fetch.
      // Only sign-pure cells at least one cell diagonal from the boundary are
      // served; the AA fringe and sign-unsafe cells fall back to the edge walk
      // on the TRUE per-frame edges.
      float fx = lut_ax * px + lut_bx * py + lut_cx;
      float fy = lut_ay * px + lut_by * py + lut_cy;
      // Clamp is memory safety: the cull disk is not contained in the LUT box,
      // so a surviving probe can map outside the grid.
      fx = hs::clamp(fx, 0.0f, lut_clamp);
      fy = hs::clamp(fy, 0.0f, lut_clamp);
      int ix = (int)fx;
      int iy = (int)fy;
      const int16_t *cell = lut_data + iy * lut_n + ix;
      int32_t q00 = cell[0], q10 = cell[1];
      int32_t q01 = cell[lut_n], q11 = cell[lut_n + 1];
      // Sign-purity + magnitude guard in quantized integer units.
      int32_t all_or = q00 | q10 | q01 | q11;
      int32_t all_and = q00 & q10 & q01 & q11;
      int32_t min_q = std::min(
          {std::abs(q00), std::abs(q10), std::abs(q01), std::abs(q11)});
      if ((all_or >= 0 || all_and < 0) && min_q > lut_q_safe) {
        HS_SCAN_METRIC(hs::g_scan_metrics.lut_hits++);
        float tx = fx - ix;
        float ty = fy - iy;
        float d0 = q00 + (q10 - q00) * tx;
        float d1 = q01 + (q11 - q01) * tx;
        plane_dist = lut_dequant * (d0 + (d1 - d0) * ty);
        lut_served = true;
        HS_PROBE_SPAN(edge_lut, hs_t);
        HS_PROBE_COUNT(n_lut);
      }
    }
    if (!lut_served) {
      HS_SCAN_METRIC(hs::g_scan_metrics.exact_hits++);
      if (probe_flags & PROBE_CONVEX) {
        plane_dist = plane_dist_convex(px, py);
        HS_PROBE_SPAN(edge_convex, hs_t);
        HS_PROBE_COUNT(n_convex);
      } else {
        bool inside;
        float dsq;
        if (probe_flags & PROBE_SECTOR) {
          HS_SCAN_METRIC(hs::g_scan_metrics.sector_hits++);
          dsq = plane_dsq_sector(px, py, inside);
          HS_PROBE_SPAN(edge_sector, hs_t);
          HS_PROBE_COUNT(n_sector);
        } else {
          dsq = plane_dsq_exact(px, py, inside);
          HS_PROBE_SPAN(edge_exact, hs_t);
          HS_PROBE_COUNT(n_exact);
        }
        // Outside probes past the (margin-carrying) reject bound skip the
        // sqrt; the caller rejects them on dist alone.
        if ((probe_flags & PROBE_LINEAR) && !inside && dsq >= reject_dsq) {
          res = DistanceResult(FAR_SENTINEL, 0.0f, FAR_SENTINEL, 0.0f, size);
          HS_PROBE_SPAN(pack, hs_t);
          return;
        }
        plane_dist = (inside ? -1.0f : 1.0f) * sqrtf(dsq);
      }
    }

    // Small faces skip the plane->angle conversion: tan(angle) ~ angle to
    // within size^2/3 of the shading gradient (< 1.5% at the 0.2 threshold).
    float raw = (probe_flags & PROBE_LINEAR) ? plane_dist
                                             : fast_atan2(plane_dist, 1.0f);
    res = DistanceResult(raw, 0.0f, raw, 0.0f, size);
    HS_PROBE_SPAN(pack, hs_t);
  }
};

// Leaf roster for the CSG composition contract: a leaf that stops satisfying
// SDFShape fails here rather than at whichever composition happens to use it.
static_assert(SDFShape<Face>);

} // namespace SDF
