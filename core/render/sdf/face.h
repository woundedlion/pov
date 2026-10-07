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
#include <new>
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

#include "render/sdf/face_lut.h"
/**
 * @brief Scratch buffer for Face computations to avoid allocations.
 */
struct FaceScratchBuffer {
  static constexpr int MAX_VERTS = 64; /**< Maximum vertices per face. */
  static constexpr size_t MAX_INTERVALS =
      2; /**< Capacity of the azimuth coverage array: one span, or two when the
              covered arc straddles theta=0. */
  std::array<math::Vector, MAX_VERTS + 1>
      poly_2d; /**< Projected 2D polygon (+1 entry to avoid modulo). */
  std::array<math::Vector, MAX_VERTS> edge_vectors; /**< Per-edge 2D vectors. */
  std::array<float, MAX_VERTS>
      edge_lengths_sq; /**< Per-edge squared lengths. */
  std::array<math::Vector, MAX_VERTS>
      planes; /**< Compacted great-circle normals during bounds; then admitted
                  sector rays (minimum-radius x/y, vertex radius squared z). */
  std::array<Interval, MAX_INTERVALS>
      intervals;                       /**< Azimuth coverage intervals. */
  std::array<float, MAX_VERTS> thetas; /**< Per-vertex azimuth angles. */
  std::array<float, MAX_VERTS>
      inv_edge_lengths_sq; /**< Reciprocal squared edge lengths. */
  std::array<float, MAX_VERTS>
      inv_edge_j; /**< Reciprocal of each edge's y-component. */
  std::array<math::Vector, MAX_VERTS + 1>
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
  /** @brief Sorted vertex rows and their crossing-edge masks. */
  struct YWalkCache {
    std::array<uint64_t, MAX_VERTS>
        masks; /**< Edges crossing each row interval. */
    std::array<uint8_t, MAX_VERTS>
        indices; /**< Original vertex index at each row. */
    std::array<uint8_t, MAX_VERTS>
        edge_flags; /**< Incident edges owned by lower/upper row traversal. */
  };
  static_assert(sizeof(YWalkCache) <= sizeof(HalfPlane) * MAX_VERTS);
  union {
    std::array<HalfPlane, MAX_VERTS>
        half_planes;   /**< Convex edge half-planes. */
    YWalkCache y_walk; /**< Exact non-convex row traversal. */
  };
  std::array<float, MAX_VERTS + 1>
      pseudo_angles; /**< Unwrapped sector angles, or sorted exact-walk rows. */
  std::array<uint32_t, MAX_VERTS + 1>
      sector_keys; /**< pseudo_angles as order-preserving integer keys. */
  /** Bumped by every Face that writes geometry here, including post-projection
   * culls; an older value identifies retargeted geometry. */
  uint32_t claim_seq = 0;
};

static_assert(sdf_max_spans<Face>::value >= FaceScratchBuffer::MAX_INTERVALS,
              "Face's span bound must cover the scratch interval array it "
              "replays");

/**
 * @brief Order-preserving unsigned key for a non-NaN float: key(a) <= key(b)
 *        exactly when a <= b, with -0.0 and +0.0 mapping to the same key.
 */
__attribute__((always_inline)) inline uint32_t angle_key(float x) {
  uint32_t u = std::bit_cast<uint32_t>(x);
  return (u & 0x80000000u) ? (0u - u) : (u + 0x80000000u);
}

/**
 * @brief Diamond pseudo-angle of (x, y) in [0, 4), strictly monotonic with
 *        atan2 but trig-free.
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
 * @details The span members view the FaceScratchBuffer handed to the
 * constructor and own none of it: the buffer must outlive the Face and back no
 * other live Face, since building a second Face over it retargets the first.
 * get_vertical_bounds() traps on a retargeted Face.
 */
struct Face {
  math::Vector center; /**< Normalized face centroid (projection axis). */
  math::Vector basis_v, basis_u,
      basis_w;              /**< Local tangent frame (v = center). */
  int count;                /**< Vertex/edge count; 0 if culled. */
  float size = 0.0f;        /**< Inradius metric for AA normalization; radians
                                 unless linear_dist. */
  float radius = 0.0f;      /**< Circumradius in the 2D projection. */
  float max_dist = 0.0f;    /**< Cull radius (circumradius plus margin). */
  float max_dist_sq = 0.0f; /**< Squared cull radius. */

  std::span<math::Vector> poly_2d; /**< Projected 2D polygon (+1 wrap entry). */
  std::span<math::Vector> edge_vectors; /**< Per-edge 2D vectors. */
  std::span<float> edge_lengths_sq;     /**< Per-edge squared lengths. */
  std::span<float> inv_edge_lengths_sq; /**< Reciprocal squared edge lengths. */
  std::span<float> inv_edge_j; /**< Reciprocal of each edge's y-component. */

  int y_min, y_max; /**< Inclusive vertical row bounds. */
  int build_height; /**< Canvas height the bounds were computed for. */
  math::LatitudeGeometry build_geometry;
  int build_width; /**< Clip width the azimuth cull ran against; 0 if unclipped. */
  const float *build_azimuth_pads; /**< Optional row padding table. */
  std::span<Interval> intervals;   /**< Azimuth coverage intervals (radians). */
  bool full_width;                 /**< True when the face spans all columns. */
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

  static constexpr int SECTOR_MIN_COUNT =
      10; /**< Minimum vertex count for the sector path. */
  std::span<const uint32_t>
      sector_keys; /**< Strictly increasing unwrapped angle keys, count+1. */
  std::span<const math::Vector>
      sector_rays; /**< Minimum-radius ray points (x/y), vertex radius squared (z). */
  float sector_min_radius_sq =
      0.0f; /**< Conservative squared origin-to-boundary distance. */
  float sector_base = 0.0f; /**< First unwrapped pseudo-angle, sgn-folded. */
  float sector_sgn = 1.0f;  /**< Winding: +1 CCW, -1 CW. Folded into the
                               table and base so the sector search
                               compares one direction. */
  std::span<const float> y_coordinates;  /**< Sorted exact-walk vertex rows. */
  std::span<const uint64_t> y_masks;     /**< Row-interval crossing edges. */
  std::span<const uint8_t> y_indices;    /**< Sorted original vertex indices. */
  std::span<const uint8_t> y_edge_flags; /**< Lower bits 0/1, upper bits 2/3:
                                            outgoing/previous incident edges. */
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

  /** @brief Retires the face: no vertices and bounds naming no row. */
  void mark_culled() {
    count = 0;
    y_min = BOUNDS_CULLED.y_min;
    y_max = BOUNDS_CULLED.y_max;
  }

  /** @brief Builds a face on an explicit pole-to-pole virtual grid. */
  Face(std::span<const math::Vector> vertices,
       std::span<const uint16_t> indices, FaceScratchBuffer &scratch,
       int virtual_height, int height, const ClipRegion *clip = nullptr,
       const float *azimuth_pads = nullptr, float bounds_margin = BOUNDS_MARGIN)
      : Face(vertices, indices, scratch,
             math::LatitudeGeometry(height, 0.0f,
                                    math::PI_F * (height - 1) /
                                        (virtual_height - 1)),
             height, clip, azimuth_pads, bounds_margin) {
    const math::LatitudeGeometry DISPLAY_GEOMETRY(height);
    HS_CHECK(build_geometry.row_to_phi(0) == DISPLAY_GEOMETRY.row_to_phi(0) &&
                 build_geometry.row_to_phi(height - 1) ==
                     DISPLAY_GEOMETRY.row_to_phi(height - 1),
             "Face: virtual height must match the display geometry");
  }

  /**
   * @brief Builds a face's projection, bounds, and edge data.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices from the pool.
   * @param scratch Scratch storage the spans alias; exclusive to this Face for
   *        its whole lifetime, and reusable only once the Face is dead.
   * @param geometry Display latitude mapping.
   * @param height Canvas height in rows.
   * @param clip Optional render clip used to tighten the face bounds.
   * @param azimuth_pads Optional latitude-adjusted padding table.
   * @param bounds_margin Angular padding around the vertical bounds. Callers
   *        rendering antialiasing fringes must pass at least `2π / width`;
   *        the default is suitable only when that reach is not required.
   */
  HS_O3_FN Face(std::span<const math::Vector> vertices,
                std::span<const uint16_t> indices, FaceScratchBuffer &scratch,
                const math::LatitudeGeometry &geometry, int height,
                const ClipRegion *clip = nullptr,
                const float *azimuth_pads = nullptr,
                float bounds_margin = BOUNDS_MARGIN)
      : build_height(height), build_geometry(geometry),
        build_width(clip ? clip->w : 0), build_azimuth_pads(azimuth_pads),
        full_width(true) {

    count = indices.size();
    HS_CHECK(count > 0 && count <= FaceScratchBuffer::MAX_VERTS,
             "Face: vertex count must be in (0, MAX_VERTS]");

    // Early vertical exit: a face whose latitude band (plus AA margin) maps to
    // an empty row range can never be rasterized.
    const bool phi_culled = [&] {
      HS_PROFILE_DEEP(face_phi_extent);
      return compute_phi_extent(vertices, indices, geometry, height,
                                bounds_margin);
    }();
    if (phi_culled) {
      mark_culled();
      return;
    }

    {
      HS_PROFILE_DEEP(face_project);
      setup_frame_and_polygon(vertices, indices, scratch);
    }

    // Cull collapsed faces by signed area; the inclusive test also rejects
    // coincident vertices when both sides are zero.
    float area2 = 0.0f;
    for (int i = 0; i < count; ++i)
      area2 +=
          poly_2d[i].x * poly_2d[i + 1].y - poly_2d[i + 1].x * poly_2d[i].y;
    if (fabsf(area2) <= COLLAPSED_AREA_RATIO * radius * radius) {
      // Scratch already holds this face's geometry, so retire any earlier
      // face's claim on the way out.
      ++scratch.claim_seq;
      mark_culled();
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

    // Azimuth half of the cull, against the conservative clip row band.
    if (clip && clip->render_y_start() < clip->render_y_end() &&
        !pole_within_circumcircle() &&
        clip_rejects_azimuth(*clip, clip->render_y_start(),
                             clip->render_y_end() - 1)) {
      ++scratch.claim_seq;
      mark_culled();
      return;
    }

    {
      HS_PROFILE_DEEP(face_bounds);

      // Vertical bounds via arc extrema + pole analysis: great-circle edges
      // bulge poleward past their vertices.
      compute_full_bounds(scratch, count, center, geometry, height, y_min,
                          y_max, bounds_margin);
      compute_inradius(scratch);
    }

    edge_vectors = std::span<math::Vector>(scratch.edge_vectors.data(), count);
    edge_lengths_sq = std::span<float>(scratch.edge_lengths_sq.data(), count);
    inv_edge_lengths_sq =
        std::span<float>(scratch.inv_edge_lengths_sq.data(), count);
    inv_edge_j = std::span<float>(scratch.inv_edge_j.data(), count);

    {
      HS_PROFILE_DEEP(face_pole);
      apply_pole_containment(height);
    }

    // Whole-face clip cull: y_min/y_max and the azimuth coverage match what the
    // scan draws.
    if (clip && clip_rejects(*clip)) {
      ++scratch.claim_seq;
      mark_culled();
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
      build_y_walk(scratch);
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
   * @details Mirrors Scan::rasterize's vertical clamp and per-fragment XClip.
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
    float pw;
    if (build_azimuth_pads) {
      pw = std::max(build_azimuth_pads[band_y_min],
                    build_azimuth_pads[band_y_max]);
    } else {
      const float sin_phi =
          std::min(sinf(build_geometry.row_to_phi(band_y_min)),
                   sinf(build_geometry.row_to_phi(band_y_max)));
      pw = face_azimuth_pad(Wd, sin_phi);
    }
    const int band_len = xc.length(Wd);
    const float COLUMN_SCALE = Wd / math::TWO_PI_F;
    for (const auto &iv : intervals) {
      // Match get_horizontal_intervals' radians-to-column rounding.
      int a = static_cast<int>(floorf((iv.start - pw) * COLUMN_SCALE));
      int b = static_cast<int>(ceilf((iv.end + pw) * COLUMN_SCALE));
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

#include "render/sdf/face_geometry.h"
};

// Leaf roster for the CSG composition contract.
static_assert(SDFShape<Face>);

} // namespace SDF
