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
#include <memory>
#include <type_traits>
#include <utility>
#include "core/math/geometry.h"
#include "core/platform/constants.h"
#include "core/render/clip.h"
#include "core/render/sdf/common.h"
#include "core/render/sdf/face_lut.h"

/**
 * @file face.h
 * @brief SDF::Face, the mesh-face leaf, and the congruence-class distance LUT
 * its probes are served from.
 */

namespace hs_test {
struct PreparedFaceTestAccess;
} // namespace hs_test

namespace SDF {

/** @brief Gnomonic coordinates with a sign-preserving horizon clamp.
 * @param v Direction to project.
 * @param center Projection axis.
 * @param u First tangent axis.
 * @param w Second tangent axis.
 * @param px Receives the first plane coordinate.
 * @param py Receives the second plane coordinate. */
HS_HOT_INLINE inline void gnomonic_project(const math::Vector &v,
                                           const math::Vector &center,
                                           const math::Vector &u,
                                           const math::Vector &w, float &px,
                                           float &py) {
  float depth = math::dot(v, center);
  if (fabsf(depth) < math::TOLERANCE)
    depth = copysignf(math::TOLERANCE, depth);
  px = math::dot(v, u) / depth;
  py = math::dot(v, w) / depth;
}

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
                  sector rays (unit direction x/y, minimum radius z). */
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
    float vx;         ///< Edge origin x.
    float vy;         ///< Edge origin y.
    float ex;         ///< Edge vector x.
    float ey;         ///< Edge vector y.
    float inv_len_sq; ///< Reciprocal squared edge length.
    float inv_ej;     ///< Reciprocal of the edge's y-component.
    uint32_t key_vy;  ///< angle_key of this vertex's y.
    uint32_t
        key_next_vy; ///< angle_key of the next vertex's y; equals key_vy for a y-degenerate edge.
  };
  std::array<EdgePacked, MAX_VERTS> packed_edges; /**< Packed per-edge data. */

  /**
   * @brief Outward unit edge normal and line offset for the convex fast path.
   */
  struct HalfPlane {
    float nx;       ///< Axis or normal x-component.
    float ny;       ///< Axis or normal y-component.
    float off;      ///< Line offset: signed distance is nx*px + ny*py + off.
    uint16_t first; ///< First edge in the collinear run, wrapping at count.
    uint16_t last;  ///< Last edge in the collinear run, wrapping at count.
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
  /** Bumped by every PreparedFace that writes geometry here, including post-projection
   * culls; an older value identifies retargeted geometry. */
  uint32_t claim_seq = 0;
};

static_assert(sdf_max_spans<Face>::value >= FaceScratchBuffer::MAX_INTERVALS,
              "Face's span bound must cover the scratch interval array it "
              "replays");

/**
 * @brief Order-preserving unsigned key for a non-NaN float: key(a) <= key(b)
 *        exactly when a <= b, with -0.0 and +0.0 mapping to the same key.
 * @param x Non-NaN value.
 * @return Order-preserving key.
 */
__attribute__((always_inline)) inline uint32_t angle_key(float x) {
  uint32_t u = std::bit_cast<uint32_t>(x);
  return (u & 0x80000000u) ? (0u - u) : (u + 0x80000000u);
}

/** @brief Rotation and angular padding of a compact face's bounds. */
struct FaceBoundsPose {
  math::Vector x;    ///< Bounds-frame x axis in endpoint space.
  math::Vector y;    ///< Bounds-frame north in endpoint space.
  math::Vector z;    ///< Bounds-frame z axis in endpoint space.
  float residual;    ///< Rotation over the face's seeing columns.
  const float *pads; ///< Shared conservative row azimuth pads.

  /**
   * @brief Expresses an endpoint direction in the bounds frame.
   * @param p Endpoint direction.
   * @return Bounds-frame direction.
   */
  HS_HOT_INLINE math::Vector project(const math::Vector &p) const {
    return {math::dot(p, x), math::dot(p, y), math::dot(p, z)};
  }
};

/**
 * @brief Prepared projection, bounds and distance state for one planar face.
 * @details The span members view the FaceScratchBuffer handed to the
 * constructor and own none of it: the buffer must outlive the PreparedFace and
 * back no other live PreparedFace. get_vertical_bounds() traps on retargeted
 * scratch storage.
 */
struct PreparedFace {
  friend struct Face;
  friend struct ::hs_test::PreparedFaceTestAccess;

  static constexpr float AZIMUTH_ERROR =
      2e-5f; ///< Vertex azimuth error, radians.
  static constexpr bool is_solid =
      true; /**< PreparedFace renders as a filled region. */
  using EdgePacked =
      FaceScratchBuffer::EdgePacked; /**< Packed per-edge record type. */
  using HalfPlane =
      FaceScratchBuffer::HalfPlane; /**< Convex edge half-plane record type. */
  static constexpr int SECTOR_MIN_COUNT =
      10; /**< Minimum vertex count for the sector path. */

private:
  math::Vector center;      /**< Normalized face centroid (projection axis). */
  math::Vector basis_u;     ///< Gnomonic x axis, perpendicular to center.
  math::Vector basis_w;     ///< Gnomonic y axis, cross(center, basis_u).
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

  int y_min;        ///< First row of the inclusive vertical bounds.
  int y_max;        ///< Last row of the inclusive vertical bounds.
  int build_height; /**< Canvas height the bounds were computed for. */
  /// Latitude mapping the bounds were computed for.
  math::LatitudeGeometry build_geometry;
  int build_width; /**< Clip width the azimuth cull ran against; 0 if unclipped. */
  const float *build_azimuth_pads; /**< Optional row padding table. */
  std::span<Interval> intervals;   /**< Azimuth coverage intervals (radians). */
  bool full_width;                 /**< True when the face spans all columns. */

  std::span<EdgePacked> packed_edges; /**< Packed per-edge data. */
  std::span<const HalfPlane>
      half_planes;     /**< Outward unit edge normals (convex faces). */
  bool convex = false; /**< 2D projection is convex; distance() takes the
                               half-plane path. */
  bool linear_dist =
      false; /**< PreparedFace is small enough to report plane distance
                               without the atan. */

  std::span<const uint32_t>
      sector_keys; /**< Strictly increasing unwrapped angle keys, count+1. */
  std::span<const math::Vector>
      sector_rays; /**< Unit ray directions (x/y), minimum boundary radius (z). */
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
  float plane_min_x = 0.0f;              ///< Left projected extent.
  float plane_max_x = 0.0f;              ///< Right projected extent.
  float plane_min_y = 0.0f;              ///< Lower projected extent.
  float plane_max_y = 0.0f;              ///< Upper projected extent.
  bool sector_ok =
      false; /**< Star-shaped about centroid; sector walk usable. */

  // Congruence-class LUT binding (bind_class_lut), flattened for the probe
  // loop: one affine map takes gnomonic (px, py) straight to LUT grid
  // coordinates (centroid shift, canonical rotation/reflection, and grid
  // scale folded together), and the sign-purity guard compares raw int16
  // magnitudes against a pre-divided quantized threshold.
  const int16_t *lut_data =
      nullptr;            /**< Class LUT samples; null = exact path. */
  int lut_n = 0;          /**< LUT grid resolution per axis. */
  int32_t lut_q_safe = 0; /**< safe_dist in quantized units. */
  float lut_ax;           ///< Grid-x coefficient on px.
  float lut_bx;           ///< Grid-x coefficient on py.
  float lut_cx;           ///< Grid-x constant term.
  float lut_ay;           ///< Grid-y coefficient on px.
  float lut_by;           ///< Grid-y coefficient on py.
  float lut_cy;           ///< Grid-y constant term.
  float lut_clamp;        /**< Grid clamp bound (n - 2). */
  float lut_dequant;      /**< int16 -> plane-unit scale. */

  const FaceScratchBuffer *scratch_owner =
      nullptr; /**< Buffer the spans view; null when the face culled before
                    claiming one. */
  uint32_t scratch_claim = 0; /**< Claim this face took on that buffer. */

public:
  /** @return Normalized face centroid, the gnomonic projection axis. */
  const math::Vector &centroid() const { return center; }
  /** @return Gnomonic x axis, perpendicular to centroid(). */
  const math::Vector &tangent_u() const { return basis_u; }
  /** @return Gnomonic y axis, cross(centroid(), tangent_u()). */
  const math::Vector &tangent_w() const { return basis_w; }
  /** @return Vertex/edge count; 0 if culled. */
  int vertex_count() const { return count; }
  /** @return Inradius metric for AA normalization; radians unless the face
   *          reports plane distance. */
  float inradius() const { return size; }
  /** @return Circumradius in the 2D projection. */
  float circumradius() const { return radius; }
  /** @return Squared radial cull radius in the 2D projection. */
  float cull_radius_sq() const { return max_dist_sq; }
  /** @return Plane-unit reach of the bound class LUT past the polygon; read
   *          only when probe_flags() has PROBE_HAS_LUT. */
  float lut_pad() const { return (lut_q_safe + 1) * lut_dequant; }
  /** @return Whether distance() takes the convex half-plane path. */
  bool is_convex() const { return convex; }
  /** @return Whether the bounds name no row. */
  bool culled() const { return y_min > y_max; }
  /** @return Claim this face took on its scratch buffer. */
  uint32_t claim() const { return scratch_claim; }

  /**
   * @brief Builds a face on an explicit pole-to-pole virtual grid.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices from the pool.
   * @param scratch Scratch storage the spans alias; exclusive to this PreparedFace for
   *        its whole lifetime.
   * @param virtual_height Rows spanning pole to pole; must reproduce the
   *        display geometry's first and last rows.
   * @param height Canvas height in rows.
   * @param clip Optional render clip used to tighten the face bounds.
   * @param azimuth_pads Optional latitude-adjusted padding table.
   * @param bounds_margin Angular padding around the vertical bounds.
   */
  HS_HOT_FLASH_MEMBER
  PreparedFace(std::span<const math::Vector> vertices,
               std::span<const uint16_t> indices, FaceScratchBuffer &scratch,
               int virtual_height, int height, const ClipRegion *clip,
               const float *azimuth_pads, float bounds_margin)
      : PreparedFace(vertices, indices, scratch,
                     virtual_grid_geometry(virtual_height, height), height,
                     clip, azimuth_pads, bounds_margin) {
    const math::LatitudeGeometry display_geometry(height);
    HS_CHECK(build_geometry.row_to_phi(0) == display_geometry.row_to_phi(0) &&
                 build_geometry.row_to_phi(height - 1) ==
                     display_geometry.row_to_phi(height - 1),
             "PreparedFace: virtual height must match the display geometry");
  }

  /**
   * @brief Builds a face's projection, bounds, and edge data.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices from the pool.
   * @param scratch Scratch storage the spans alias; exclusive to this PreparedFace for
   *        its whole lifetime, and reusable only once the PreparedFace is dead.
   * @param geometry Display latitude mapping.
   * @param height Canvas height in rows.
   * @param clip Optional render clip used to tighten the face bounds.
   * @param azimuth_pads Optional latitude-adjusted padding table.
   * @param bounds_pose Optional seeing-time frame for exact-distance face bounds;
   *        requires @p clip and replaces @p bounds_margin.
   * @param latitude_checked Whether this vertex/clip latitude band already passed.
   * @param bounds_margin Angular padding around the vertical bounds. Callers
   *        rendering antialiasing fringes must pass at least `2π / width`.
   */
  HS_HOT_FLASH_MEMBER
  PreparedFace(std::span<const math::Vector> vertices,
               std::span<const uint16_t> indices, FaceScratchBuffer &scratch,
               const math::LatitudeGeometry &geometry, int height,
               const ClipRegion *clip, const float *azimuth_pads,
               float bounds_margin, const FaceBoundsPose *bounds_pose = nullptr,
               bool latitude_checked = false)
      : build_height(height), build_geometry(geometry),
        build_width(clip ? clip->w : 0), build_azimuth_pads(azimuth_pads),
        full_width(true) {
    build_geometry_tables(vertices, indices, scratch, geometry, height, clip,
                          bounds_margin, bounds_pose, latitude_checked);
  }

  /**
   * @brief Tests whether the clip band excludes the whole face.
   * @param cr Clip region (display bounds plus render margin).
   * @return True when no in-band pixel can be produced — the face's vertical
   *         band lies outside the render rows, or its azimuth coverage lies
   *         outside the render columns. A full-width face or an inactive x-clip
   *         never rejects on the horizontal axis.
   * @details Mirrors the scan's vertical clamp (`Scan::clamp_rows_to_clip`)
   * and its run-level column-arc clip (`Scan::clip_run`).
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
    const int wd = cr.w;
    float pw;
    if (build_azimuth_pads) {
      pw = std::max(build_azimuth_pads[band_y_min],
                    build_azimuth_pads[band_y_max]);
    } else {
      const float sin_phi =
          std::min(sinf(build_geometry.row_to_phi(band_y_min)),
                   sinf(build_geometry.row_to_phi(band_y_max)));
      pw = face_azimuth_pad(wd, sin_phi);
    }
    const int band_len = xc.length(wd);
    const float column_scale = wd / math::TWO_PI_F;
    for (const auto &iv : intervals) {
      // Match get_horizontal_intervals' radians-to-column rounding.
      int a = static_cast<int>(floorf((iv.start - pw) * column_scale));
      int b = static_cast<int>(ceilf((iv.end + pw) * column_scale));
      int len = b - a;
      if (len <= 0)
        continue;
      if (len > wd)
        len = wd;
      const int s = ((a % wd) + wd) % wd;
      if (ClipRegion::arcs_overlap(xc.rs, band_len, s, len, wd))
        return false;
    }
    return true;
  }

  /**
   * @brief Latitude-band reject for the face.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices.
   * @param geometry Display latitude mapping.
   * @param height Canvas height in rows.
   * @param bounds_margin Angular padding around the vertical bounds.
   * @param clip Optional render-row clip.
   * @return True when the padded face cannot overlap the canvas or render clip.
   */
  __attribute__((always_inline)) static bool
  phi_band_rejects(std::span<const math::Vector> vertices,
                   std::span<const uint16_t> indices,
                   const math::LatitudeGeometry &geometry, int height,
                   float bounds_margin, const ClipRegion *clip);

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
   * @details The boundary band is tested first: a pole exactly on a vertex or
   * edge would otherwise leave pole_inside_polygon's crossing parity to
   * rounding.
   */
  PoleHit pole_hit(float ppx, float ppy) const;

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
  refine_phi_from_arc_extremum(const math::Vector &n, const math::Vector &v1,
                               const math::Vector &v2, float &min_phi,
                               float &max_phi);

  /**
   * @brief Returns the face's precomputed inclusive row bounds.
   * @tparam H Canvas height in rows; must match the construction height the
   * bounds were computed for.
   * @return The stored {y_min, y_max} bounds.
   */
  template <int H> Bounds get_vertical_bounds() const;

  /**
   * @brief Reports whether rounded azimuth intervals change across a row band.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @param y_lo First row in the band.
   * @param y_hi Last row in the band.
   * @return True when the band requires per-row interval construction.
   */
  template <int W, int H>
  bool horizontal_intervals_vary_by_row(int y_lo, int y_hi) const;

  /**
   * @brief Whether two rows emit the same rounded azimuth intervals.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param first_y First raster row.
   * @param second_y Second raster row.
   * @return True when both rows emit identical column intervals.
   */
  template <int W, int H>
  bool horizontal_intervals_equal_rows(int first_y, int second_y) const;

  /**
   * @brief Emits the face's azimuth-coverage intervals for a row.
   * @tparam W Canvas width in columns; must match the clip width the
   * construction-time azimuth cull ran against.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y Row index, which sets the pole widening.
   * @param out Sink accepting (float start, float end).
   * @return True when the row was handled, possibly with no intervals outside
   *   the face's rows; false requests a full scan.
   * @details The pad is an azimuth angle, so it holds a whole pixel of AA reach
   * only at the equator; at colatitude phi one pad p of great-circle reach
   * subtends asin(p / sin phi), reaching the whole row once sin phi <= p.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const;

  /**
   * @brief Returns the latitude-adjusted AA padding for one raster row.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param y Raster row.
   * @return Azimuth padding in radians.
   */
  template <int W, int H> float azimuth_pad_at_row(int y) const;

  /** @brief Squared distance to a clamped packed segment.
   * @param edge Packed segment; a zero inverse length selects its start.
   * @param px,py Query point in the tangent plane.
   * @return Squared point-to-segment distance. */
  HS_HOT_INLINE static float segment_dsq(const EdgePacked &edge, float px,
                                         float py);

  /** @brief Squared distance to a clamped planar segment.
   * @param vx,vy Segment start.
   * @param ex,ey Segment displacement.
   * @param inv_length2 Inverse squared length; zero selects the start.
   * @param px,py Query point.
   * @return Squared point-to-segment distance. */
  HS_HOT_INLINE static float segment_dsq(float vx, float vy, float ex, float ey,
                                         float inv_length2, float px, float py);

  /**
   * @brief Signed planar distance through convex edge slabs and vertex cones.
   * @param px Gnomonic x of the query point.
   * @param py Gnomonic y of the query point.
   * @param reject_distance Distances at or above this cutoff may be capped.
   * @return Signed distance in the tangent plane (negative inside).
   */
  HS_O3_FN float plane_dist_convex(float px, float py,
                                   float reject_distance = FLT_MAX) const;

  /**
   * @brief Squared planar distance via the exact per-edge walk.
   * @param px Gnomonic x of the query point.
   * @param py Gnomonic y of the query point.
   * @param inside_out Set true when the query lies inside the polygon; carries
   * the sign the squared return cannot.
   * @return Squared distance to the nearest edge, in the tangent plane.
   */
  HS_O3_FN float plane_dsq_exact(float px, float py, bool &inside_out) const;

  /**
   * @brief Exact squared distance and parity by sorted vertex-row traversal.
   * @details Checks all edges crossing the query row, then incident edges at
   * neighboring vertex rows. Unvisited segments lie outside the visited strip;
   * its vertical gaps bound their distance. Requires a populated y-walk cache.
   * @param px Gnomonic x of the query point.
   * @param py Gnomonic y of the query point.
   * @param inside_out Set true when the query lies inside the polygon.
   * @return Squared distance to the nearest edge, in the tangent plane.
   */
  HS_O3_FN HS_NOINLINE_NOCLONE float plane_dsq_y_walk(float px, float py,
                                                      bool &inside_out) const;

  /**
   * @brief Squared planar distance via the concave sector walk.
   * @param px Gnomonic x of the query point.
   * @param py Gnomonic y of the query point.
   * @param inside_out Set true when the query lies inside the polygon; carries
   * the sign the squared return cannot.
   * @param reject_dsq Outside squared distances at or above this threshold
   * may return the threshold; FLT_MAX keeps the distance exact.
   * @param column_sector Optional sector carried between probes of one column.
   * @return Squared distance to the nearest edge, capped outside when requested.
   * @details Searches the fan sector and its neighbors, certifying the minimum
   * against the omitted edges' angular wedge and radial bound, expanding until
   * certified or every segment has been evaluated. Requires sector_ok.
   */
  HS_O3_FN float plane_dsq_sector(float px, float py, bool &inside_out,
                                  float reject_dsq = FLT_MAX,
                                  int *column_sector = nullptr) const;

  /** @brief Face has a congruence-class distance LUT (bind_class_lut). */
  static constexpr uint32_t PROBE_HAS_LUT = 1u << 0;

  /** @brief Use convex half-plane distance. */
  static constexpr uint32_t PROBE_CONVEX = 1u << 1;

  /** @brief Use the sector-indexed edge walk. */
  static constexpr uint32_t PROBE_SECTOR = 1u << 2;

  /** @brief Reject against a compact projected rectangle. */
  static constexpr uint32_t PROBE_PLANE_BOX = 1u << 4;

  /** @brief Return distances in the gnomonic plane. */
  static constexpr uint32_t PROBE_LINEAR = 1u << 3;

  /** @return Packed distance-path flags for repeated probes of this face. */
  uint32_t probe_flags() const;

  /**
   * @brief Computes signed distance to the face, writing into res.
   * @tparam ComputeUVs Accepted for interface parity; the face stores no UVs.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = raw_dist = the signed edge distance,
   *        size = max(inradius, 0.25 * circumradius), both in gnomonic plane units on a linear_dist
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
  HS_O3_FN void distance(const math::Vector &p, DistanceResult &res,
                         float reject_dsq = FLT_MAX) const;

  /**
   * @brief Evaluates distance with cached flags, recomputing the radial cull cosine.
   * @tparam ComputeUVs Accepted for interface parity; the face stores no UVs.
   * @param p Unit query direction.
   * @param res Output distance result.
   * @param reject_dsq Conservative squared plane-distance rejection threshold.
   * @param probe_flags Flags captured by probe_flags() after geometry/LUT updates.
   * @details The five-argument overload accepts a cached cosine for repeated
   * probes.
   */
  template <bool ComputeUVs = true>
  HS_O3_FN void distance_with_flags(const math::Vector &p, DistanceResult &res,
                                    float reject_dsq,
                                    uint32_t probe_flags) const;

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
   * @param min_cos Cosine of the radial cull angle.
   * @param reject_plane Cutoff at which the caller rejects a planar distance.
   *        Convex distances at or above it may be capped; PROBE_PLANE_BOX
   *        returns FAR_SENTINEL when its box distance reaches it.
   *        FLT_MAX disables both shortcuts.
   */
  template <bool ComputeUVs = true>
  HS_O3_FN void distance_with_flags(const math::Vector &p, DistanceResult &res,
                                    float reject_dsq, uint32_t probe_flags,
                                    float min_cos,
                                    float reject_plane = FLT_MAX) const;

  /**
   * @brief Evaluates a direction already expressed in the face's frame.
   * @tparam ComputeUVs Accepted for interface parity.
   * @param cos_angle Direction dotted with the face center.
   * @param u Direction dotted with basis_u.
   * @param w Direction dotted with basis_w.
   * @param res Output distance result.
   * @param reject_dsq Conservative squared planar rejection threshold.
   * @param probe_flags Cached distance-path flags.
   * @param min_cos Cosine of the radial cull angle.
   * @param reject_plane Cutoff the caller rejects on; convex distances may be
   * capped and PROBE_PLANE_BOX may return FAR_SENTINEL at or above it.
   * FLT_MAX disables both shortcuts.
   * @param column_sector Optional sector carried down a projected column.
   */
  template <bool ComputeUVs = true>
  HS_O3_FN HS_HOT_INLINE void
  distance_projected(float cos_angle, float u, float w, DistanceResult &res,
                     float reject_dsq, uint32_t probe_flags, float min_cos,
                     float reject_plane = FLT_MAX,
                     int *column_sector = nullptr) const;

private:
  /** @brief Retires the face: no vertices and bounds naming no row. */
  void mark_culled() {
    count = 0;
    y_min = BOUNDS_CULLED.y_min;
    y_max = BOUNDS_CULLED.y_max;
  }

  /**
   * @brief Builds projected geometry and raster tables.
   * @param vertices Shared vertex pool.
   * @param indices PreparedFace vertex indices.
   * @param scratch Exclusive face scratch storage.
   * @param geometry Display latitude mapping.
   * @param height Canvas height in rows.
   * @param clip Optional render clip.
   * @param bounds_margin Angular padding around the bounds.
   * @param bounds_pose Optional seeing-time bounds frame; requires @p clip and
   *        replaces @p bounds_margin.
   * @param latitude_checked Whether this vertex/clip latitude band already passed.
   */
  HS_HOT_FLASH_MEMBER void build_geometry_tables(
      std::span<const math::Vector> vertices, std::span<const uint16_t> indices,
      FaceScratchBuffer &scratch, const math::LatitudeGeometry &geometry,
      int height, const ClipRegion *clip, float bounds_margin,
      const FaceBoundsPose *bounds_pose, bool latitude_checked) {
    HS_CHECK(!bounds_pose || clip, "PreparedFace: bounds pose requires a clip");
    count = indices.size();
    HS_CHECK(count > 0 && count <= FaceScratchBuffer::MAX_VERTS,
             "PreparedFace: vertex count must be in (0, MAX_VERTS]");

    // Early vertical exit: a face whose latitude band (plus AA margin) maps to
    // an empty row range can never be rasterized.
    const bool phi_culled = !latitude_checked && [&]() HS_HOT_INLINE {
      HS_PROFILE_DEEP(face_phi_extent);
      return phi_band_rejects(vertices, indices, geometry, height,
                              bounds_margin, clip);
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

    math::Vector bounds_center = center;
    if (bounds_pose) {
      bounds_center = bounds_pose->project(center);
      for (int i = 0; i <= count; ++i)
        scratch.verts_3d[i] = bounds_pose->project(scratch.verts_3d[i]);
      bounds_margin = math::TWO_PI_F / build_width + bounds_pose->residual;
      build_azimuth_pads = bounds_pose->pads;
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

    compute_edges(scratch);
    {
      HS_PROFILE_DEEP(face_bounds);

      // Vertical bounds via arc extrema + pole analysis: great-circle edges
      // bulge poleward past their vertices.
      compute_full_bounds(scratch, count, bounds_center, geometry, height,
                          y_min, y_max, bounds_margin);
    }

    {
      HS_PROFILE_DEEP(face_pole);
      if (!bounds_pose)
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
      const float convex_radius = build_half_planes(scratch, area2);
      compute_inradius(scratch, convex_radius);
      pack_edges(scratch);
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
   * @brief Builds the local tangent frame, gnomonic 2D projection, and 3D
   * arrays.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices.
   * @param scratch Scratch storage receiving poly_2d and verts_3d.
   * @details Sets basis_u/w, poly_2d (with circumradius), and the local 3D
   * vertex array.
   */
  __attribute__((always_inline)) void
  setup_frame_and_polygon(std::span<const math::Vector> vertices,
                          std::span<const uint16_t> indices,
                          FaceScratchBuffer &scratch);

  /**
   * @brief Builds projected edge vectors and their reciprocals.
   * @param scratch Face geometry receiving the edge tables.
   */
  __attribute__((always_inline)) void compute_edges(FaceScratchBuffer &scratch);

  /**
   * @brief Computes the face "size" (inradius) from the projected polygon.
   * @param scratch Scratch storage holding poly_2d, per-edge vectors and
   * reciprocal squared lengths (inv_edge_lengths_sq) from compute_edges,
   * which must run first.
   * @param convex_radius Origin-to-boundary distance from accepted half-planes,
   * or a negative value for the segment walk.
   * @details size = minimum distance from the projected centroid to any edge,
   * floored to a fraction of the circumradius for degenerate slivers. A large
   * face converts it to radians, matching the metric distance() reports.
   */
  __attribute__((always_inline)) void
  compute_inradius(const FaceScratchBuffer &scratch,
                   float convex_radius = -1.0f);

  /**
   * @brief Packs per-edge data contiguously for the distance() fallback.
   * @param scratch Scratch storage receiving packed_edges.
   */
  __attribute__((always_inline)) void pack_edges(FaceScratchBuffer &scratch);

  /**
   * @brief Detects a convex 2D projection and builds its edge half-planes.
   * @param scratch Scratch storage receiving half_planes.
   * @param area2 Twice the polygon's signed area, from the collapsed-face cull.
   * @return Origin-to-boundary distance for an admitted convex polygon, else -1.
   * @details The maximum edge half-plane distance is exact inside a convex
   * polygon and in exterior edge slabs; vertex cones use endpoint distance.
   * Concave, self-intersecting, degenerate-edged, or wrongly-oriented polygons
   * leave convex false and distance() on the exact walk.
   */
  __attribute__((always_inline)) float
  build_half_planes(FaceScratchBuffer &scratch, float area2);

  /**
   * @brief Builds the angular sector table for the concave sector walk.
   * @param scratch Scratch storage receiving the unwrapped vertex
   * pseudo-angles.
   * @details Qualifies non-convex faces with at least SECTOR_MIN_COUNT
   * vertices that are star-shaped about the projected centroid (vertex
   * pseudo-angles strictly monotonic over one full turn, with every edge
   * facing the origin). Otherwise sector_ok stays false.
   */
  __attribute__((always_inline)) void build_sectors(FaceScratchBuffer &scratch);

  /**
   * @brief Builds vertex-row crossing masks for exact non-convex probes.
   * @param scratch Buffer whose y-walk cache is populated.
   */
  __attribute__((always_inline)) void build_y_walk(FaceScratchBuffer &scratch);

  /**
   * @brief Aligns the current projection to a canonical class shape and binds
   *        its distance LUT.
   * @param lut Canonical-frame LUT for the face's congruence class.
   * @param canon_xy Canonical centered 2D polygon, x/y pairs.
   * @param vert_offset Cyclic offset aligning mesh vertex order to canonical.
   * @param reflected True for the mirror family.
   * @return False when the correlation is degenerate or a vertex deviates past
   *         ALIGN_MAX_DEV_DIAGS; the face then keeps the exact path.
   * @details One complex correlation over the vertices recovers the in-plane
   * rotation placing the canonical shape at the face's least-squares pose.
   * Rotational-symmetry ambiguity is harmless (the LUT is invariant under the
   * shape's symmetry group).
   */
  HS_FLASH_MEMBER bool bind_class_lut(const ClassLut *lut,
                                      const float *canon_xy, int vert_offset,
                                      bool reflected);

  /**
   * @brief Fills the scratch vertex azimuths the interval pass consumes.
   * @param scratch Scratch storage holding verts_3d and receiving thetas.
   */
  __attribute__((always_inline)) void
  compute_thetas(FaceScratchBuffer &scratch) const;

  /**
   * @brief Computes the face's azimuth coverage intervals.
   * @param scratch Scratch storage holding thetas and receiving intervals.
   * @details Finds the largest angular gap between vertices; if it exceeds pi
   * the face does not wrap, so the complementary horizontal interval(s) are
   * emitted, else it spans full width. Coarse: only the single largest gap is
   * excised.
   */
  __attribute__((always_inline)) void
  compute_azimuth_intervals(FaceScratchBuffer &scratch);

  /**
   * @brief Necessary condition for apply_pole_containment to fire.
   * @return True when a pole falls within the gnomonic circumcircle of the
   *         face's vertices (one test covers both poles); false rules out
   *         pole containment.
   */
  __attribute__((always_inline)) bool pole_within_circumcircle() const;

  /**
   * @brief Ray-crossing test of a projected pole against the 2D face polygon.
   * @param ppx Projected pole x in the face's 2D basis.
   * @param ppy Projected pole y in the face's 2D basis.
   * @return true if (ppx, ppy) lies inside the polygon.
   */
  __attribute__((always_inline)) bool pole_inside_polygon(float ppx,
                                                          float ppy) const;

  /**
   * @brief Extends the vertical bounds when the face reaches a pole.
   * @param height Canvas height in rows.
   * @details A face enclosing a pole wraps every azimuth, so it takes full
   * width outright. A face that only meets the pole on its boundary keeps its
   * azimuth wedge, which get_horizontal_intervals widens per row.
   */
  __attribute__((always_inline)) void apply_pole_containment(int height);

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
                           const math::Vector &center, float &min_phi,
                           float &max_phi);

  /**
   * @brief Vertical bounds (arc extrema + pole-plane snap) for every face.
   * @param scratch Scratch storage holding verts_3d; receives the compacted
   * great-circle planes.
   * @param count Vertex/edge count.
   * @param center Normalized face centroid.
   * @param geometry Display latitude mapping.
   * @param height Canvas height in rows.
   * @param y_min_out Output: first covered row.
   * @param y_max_out Output: last covered row.
   * @param bounds_margin Angular padding around the vertical bounds.
   */
  HS_HOT_FLASH_MEMBER static void
  compute_full_bounds(FaceScratchBuffer &scratch, int count,
                      const math::Vector &center,
                      const math::LatitudeGeometry &geometry, int height,
                      int &y_min_out, int &y_max_out, float bounds_margin);
};

} // namespace SDF

#include "core/render/sdf/face_geometry.h"

namespace SDF {

/** @brief Draw-scoped face tables and prepared geometry. */
struct FaceStorage {
  FaceScratchBuffer tables; ///< Reusable vertex and edge tables.
  union {
    PreparedFace geometry; ///< Derived projection, bounds and distance state.
  };
  uint32_t generation = 0; ///< Advances before any geometry replacement.

  /** @brief Leaves geometry uninitialized until a face is prepared. */
  FaceStorage() {}
  FaceStorage(const FaceStorage &) = delete;
  FaceStorage &operator=(const FaceStorage &) = delete;
};

static_assert(std::is_trivially_destructible_v<PreparedFace>);

/** @brief Generation-checked handle to one prepared face. */
struct Face {
private:
  FaceStorage *storage; ///< Arena storage containing the current preparation.
  uint32_t claim;       ///< Generation captured by this handle.

  /** @brief Validates the handle's captured storage generation. */
  void validate() const {
    HS_CHECK(claim == storage->generation, "SDF::Face storage was retargeted");
  }

public:
  static constexpr bool is_solid = true; ///< Filled-region SDF leaf.

  /** @brief Prepares geometry in reusable face storage.
   * @tparam Args Display geometry and bounds argument types.
   * @param vertices Shared vertex pool.
   * @param indices Face vertex indices.
   * @param storage Exclusive draw-scoped face storage.
   * @param args Display geometry and bounds arguments. */
  template <typename... Args>
  Face(std::span<const math::Vector> vertices,
       std::span<const uint16_t> indices, FaceStorage &storage, Args &&...args)
      : storage(&storage), claim(++storage.generation) {
    std::construct_at(&storage.geometry, vertices, indices, storage.tables,
                      std::forward<Args>(args)...);
  }

  /** @brief Binds immutable geometry for this storage generation.
   * @return Current prepared geometry; stale handles trap. */
  const PreparedFace &prepared() const {
    validate();
    return storage->geometry;
  }

  /** @brief Binds a congruence-class distance LUT to this generation's
   *         geometry.
   * @param lut Canonical-frame LUT for the face's congruence class.
   * @param canon_xy Canonical centered 2D polygon, x/y pairs.
   * @param vert_offset Cyclic offset aligning mesh vertex order to canonical.
   * @param reflected True for the mirror family.
   * @return False when the face keeps the exact path; stale handles trap. */
  bool bind_class_lut(const ClassLut *lut, const float *canon_xy,
                      int vert_offset, bool reflected) {
    validate();
    return storage->geometry.bind_class_lut(lut, canon_xy, vert_offset,
                                            reflected);
  }

  /** @brief Returns the prepared inclusive row bounds.
   * @tparam H Canvas height in rows.
   * @return Inclusive row bounds. */
  template <int H> Bounds get_vertical_bounds() const {
    return prepared().template get_vertical_bounds<H>();
  }

  /** @brief Emits prepared azimuth intervals for a row.
   * @tparam W,H Canvas dimensions.
   * @tparam OutputIt Interval sink type.
   * @param y Raster row.
   * @param out Interval sink.
   * @return Whether the row was handled. */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    return prepared().template get_horizontal_intervals<W, H>(y, out);
  }

  /** @brief Reports whether intervals vary across a latitude band.
   * @tparam W,H Canvas dimensions.
   * @param first,last Inclusive row band.
   * @return Whether the scan requires per-row interval construction. */
  template <int W, int H>
  bool horizontal_intervals_vary_by_row(int first, int last) const {
    return prepared().template horizontal_intervals_vary_by_row<W, H>(first,
                                                                      last);
  }

  /** @brief Compares the prepared intervals at two rows.
   * @tparam W,H Canvas dimensions.
   * @param first,second Raster rows.
   * @return Whether both rows emit identical rounded intervals. */
  template <int W, int H>
  bool horizontal_intervals_equal_rows(int first, int second) const {
    return prepared().template horizontal_intervals_equal_rows<W, H>(first,
                                                                     second);
  }

  /** @brief Evaluates the prepared signed distance.
   * @tparam ComputeUVs Whether shader registers are requested.
   * @param p Unit query direction.
   * @param res Output distance and shader registers.
   * @param reject_dsq Outside-distance rejection threshold squared. */
  template <bool ComputeUVs = true>
  void distance(const math::Vector &p, DistanceResult &res,
                float reject_dsq = FLT_MAX) const {
    prepared().template distance<ComputeUVs>(p, res, reject_dsq);
  }
};

static_assert(sizeof(Face) == 2 * sizeof(void *));
static_assert(SDFShape<PreparedFace>);

// Leaf roster for the CSG composition contract.
static_assert(SDFShape<Face>);

} // namespace SDF
