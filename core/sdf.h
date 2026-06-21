/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <cfloat>
#include <span>
#include "geometry.h"
#include "constants.h"
#include "concepts.h"
#include "static_circular_buffer.h"
#include "util.h"

namespace SDF {

// --- Rasterization constants ------------------------------------------------
/** Margin added to bounding boxes for AA kernel width (small faces). */
static constexpr float BOUNDS_MARGIN = 0.05f;
/** Expanded margin for shapes with variable thickness or large faces. */
static constexpr float BOUNDS_MARGIN_WIDE = 0.1f;
/** Minimum horizontal projection length for interval culling. */
static constexpr float MIN_HORIZONTAL_PROJ = 0.01f;
/** Epsilon for near-zero denominators in interval math. */
static constexpr float INTERVAL_DENOM_EPS = 1e-6f;
/** Threshold for near-pole ring approximation safety. */
static constexpr float POLE_SAFE_MARGIN = 0.05f;
/** Inner/outer radius ratio for star shapes (1/φ² ≈ 0.382). */
static constexpr float STAR_INNER_RATIO = 0.382f;
/** Vertex phi-span threshold below which Face uses the fast vertex-only
 *  bounding path instead of the full arc-extrema computation. */
static constexpr float FAST_BOUNDS_PHI_THRESHOLD = 0.3f;
/** Minimum inradius-to-circumradius ratio used to floor Face::size,
 *  preventing degenerate near-zero inradii from collapsing AA. */
static constexpr float MIN_SIZE_RADIUS_RATIO = 0.25f;

/** Maximum disjoint scanline spans budgeted per row. scan_region's seam-split
 *  `norm` buffer is bound to exactly 2x this (one span can split in two at the
 *  x=0 seam); see the static_assert there. */
inline constexpr size_t kIntervalSpanCap = 32;

/** Per-row scanline interval buffer. Fixed capacity, accumulate-only. */
using IntervalBuffer =
    StaticCircularBuffer<std::pair<float, float>, kIntervalSpanCap>;

/**
 * @brief Append a scanline interval, trapping on overflow.
 * @param buf Per-row interval buffer to append to.
 * @param start Interval start column (float).
 * @param end Interval end column (float).
 * @details The CSG ops and scan_region accumulate per-row intervals into a
 * fixed-capacity buffer. StaticCircularBuffer::push_back evicts the OLDEST entry
 * when full (correct for trails, wrong here) — so an overflow would silently
 * drop geometry. A row exceeding the capacity is a sizing bug (e.g. a deeply
 * nested CSG producing more disjoint spans than budgeted), so trap at the
 * violation site instead of dropping coverage (fail-fast).
 */
inline void push_interval(IntervalBuffer &buf, float start, float end) {
  HS_CHECK(!buf.is_full(),
           "SDF scanline interval buffer overflow (>32 spans in one row)");
  buf.push_back({start, end});
}

/**
 * @brief Insertion-sort an interval buffer in place by start coordinate.
 * @param buf Per-row interval buffer to sort in place.
 * @details Uses raw-pointer indexing (the buffer is freshly built, head == 0, so
 * it is contiguous from index 0), avoiding the per-access modulo. Shared by
 * merge_intervals and the Subtract set-difference.
 */
inline void sort_intervals_by_start(IntervalBuffer &buf) {
  auto *data = &buf[0];
  size_t n = buf.size();
  for (size_t i = 1; i < n; ++i) {
    auto key = data[i];
    int j = static_cast<int>(i) - 1;
    while (j >= 0 && data[j].first > key.first) {
      data[j + 1] = data[j];
      --j;
    }
    data[j + 1] = key;
  }
}

/**
 * @brief Sort an interval buffer by start, then emit the union of overlapping
 * intervals via out(start, end). Shared by Union/SmoothUnion.
 * @tparam OutputIt Sink type invoked as out(float start, float end).
 * @param merged Per-row interval buffer (sorted and merged in place).
 * @param out Sink receiving each merged interval.
 * @details Precondition: `merged` is non-empty (callers guard with is_empty()).
 * Templated on the output sink so it inlines at -O3.
 */
template <typename OutputIt>
inline void merge_intervals(IntervalBuffer &merged, OutputIt out) {
  sort_intervals_by_start(merged);
  auto *data = &merged[0];
  size_t n = merged.size();
  float cur_start = data[0].first;
  float cur_end = data[0].second;
  for (size_t i = 1; i < n; ++i) {
    if (data[i].first <= cur_end) {
      cur_end = std::max(cur_end, data[i].second);
    } else {
      out(cur_start, cur_end);
      cur_start = data[i].first;
      cur_end = data[i].second;
    }
  }
  out(cur_start, cur_end);
}

/**
 * @brief Fold an angle into [0, π], equivalent to acosf(cosf(x)) without trig.
 * @param x Angle in radians.
 * @return Folded angle in [0, π].
 */
inline float clamp_phi(float x) {
  if (x < 0.0f)
    return -x;
  if (x > PI_F)
    return 2.0f * PI_F - x;
  return x;
}
/**
 * @brief Vertical scanline bounds (inclusive min/max row index).
 */
struct Bounds {
  int y_min, y_max; /**< Inclusive first/last row covered. */
};
/**
 * @brief Horizontal scanline interval (start/end column index).
 */
struct Interval {
  int start, end; /**< Half-open column range [start, end). */
};
/**
 * @brief Result of a signed distance query.
 *
 * `dist` and `size` have fixed meanings, but `t`, `raw_dist` and `aux` are
 * deliberately overloaded per shape — each shape packs whatever supplementary
 * value its fragment shader needs into them, so the authoritative meaning for a
 * given shape is that shape's own "Returns:" docblock (e.g. Ring, Flower,
 * SphericalPolygon). The scan rasterizer copies them straight into the Fragment
 * register file with no reinterpretation (see Scan::process_pixel):
 *   t        -> Fragment::v0
 *   raw_dist -> Fragment::v1
 *   aux      -> Fragment::v3
 *   size     -> Fragment::size
 * (Fragment::v2 is reserved and always 0.)
 */
struct DistanceResult {
  float dist;        /**< Signed distance (negative inside); always this meaning. */
  float t;           /**< Per-shape: normalized parameter (0-1) or angle. */
  float raw_dist;    /**< Per-shape: unsigned or supplementary distance. */
  float aux;         /**< Per-shape: auxiliary value (e.g. barycentric coordinate). */
  float size = 1.0f; /**< Size metric for AA-falloff normalization. */

  /**
   * @brief Default-constructs an uninitialized result.
   */
  DistanceResult() = default;
  /**
   * @brief Packs the five per-pixel SDF outputs into a result.
   * @param d Signed distance (negative inside).
   * @param t_val Per-shape normalized parameter (0-1) or angle.
   * @param rd Per-shape unsigned or supplementary distance.
   * @param ax Per-shape auxiliary value.
   * @param sz Size metric for AA-falloff normalization.
   */
  DistanceResult(float d, float t_val, float rd, float ax, float sz)
      : dist(d), t(t_val), raw_dist(rd), aux(ax), size(sz) {}
};

/**
 * @brief Emit the single horizontal interval where a row crosses a great-circle
 * "cap" of half-angle `acos(cos_cap)` centred on an axis whose projection onto
 * the scan plane is (ny, R_val, alpha_angle). Shared by PlanarPolygon /
 * SphericalPolygon / Star, whose scanline math is otherwise identical.
 *
 * @tparam W Canvas width in columns.
 * @tparam OutputIt Sink type invoked as out(float start, float end).
 * @param cos_cap cos of the cap's angular radius (cosf(thickness) etc.).
 * @param ny y-component of the cap axis.
 * @param R_val Horizontal projection length of the cap axis.
 * @param alpha_angle Azimuth of the cap axis (radians).
 * @param cos_phi Cosine of the row's polar angle.
 * @param sin_phi Sine of the row's polar angle.
 * @param reject_full_width Star-only: drop the row to a full scan when the
 *        interval would span the whole width (the others never set this).
 * @param out Sink accepting (float start, float end).
 * @return false to request a full-width fallback scan, true if the (possibly
 *         empty) interval was handled.
 */
template <int W, typename OutputIt>
inline bool emit_cap_interval(float cos_cap, float ny, float R_val,
                              float alpha_angle, float cos_phi, float sin_phi,
                              bool reject_full_width, OutputIt out) {
  if (R_val < MIN_HORIZONTAL_PROJ)
    return false;

  float denom = R_val * sin_phi;
  if (std::abs(denom) < INTERVAL_DENOM_EPS)
    return false;

  float C_min = (cos_cap - ny * cos_phi) / denom;
  if (C_min > 1.0f)
    return true; // Completely outside the cap
  if (C_min < -1.0f)
    return false; // Full scan fallback (simplification)

  // fast_acos: ~1.3e-4 rad peak error ≈ 0.006 px at W=288, far under the
  // floor/ceil pad below. Matches the Ring/DistortedRing scanline path.
  float d_alpha = fast_acos(C_min);
  float scale = W / (2.0f * PI_F);
  float x1 = floorf((alpha_angle - d_alpha) * scale);
  float x2 = ceilf((alpha_angle + d_alpha) * scale);
  if (reject_full_width && x2 - x1 >= W)
    return false;
  out(x1, x2);
  return true;
}

/**
 * @brief Map an angular band [phi_min, phi_max] (radians, polar angle) to the
 * inclusive scanline row range it covers, clamped to [0, height-1].
 * @param phi_min Lower polar-angle edge of the band (radians).
 * @param phi_max Upper polar-angle edge of the band (radians).
 * @param h_virt Virtual row count (height plus pole offset).
 * @param height Canvas height in rows.
 * @return Inclusive row bounds covering the band.
 * @details phi spans [0,π] over (h_virt-1) virtual rows; the lower edge floors
 * and the upper edge ceils so a partially-covered row is never dropped.
 * Out-of-range phi needs no pre-clamp: the row clamps fold a sub-0 lower edge to
 * 0 and a past-π upper edge to height-1. This runtime form is the single source
 * for the floor/ceil conversion — the templated `phi_bounds_to_rows<H>` and the
 * Face's construction-time bounds (which only know `h_virt`/`height` at runtime)
 * both route through it.
 */
inline Bounds phi_bounds_to_rows(float phi_min, float phi_max, int h_virt,
                                 int height) {
  int y_min =
      std::max(0, static_cast<int>(floorf((phi_min * (h_virt - 1)) / PI_F)));
  int y_max = std::min(
      height - 1, static_cast<int>(ceilf((phi_max * (h_virt - 1)) / PI_F)));
  return {y_min, y_max};
}

/**
 * @brief Compile-time-H wrapper over the runtime phi-to-row conversion.
 * @tparam H Canvas height in rows.
 * @param phi_min Lower polar-angle edge of the band (radians).
 * @param phi_max Upper polar-angle edge of the band (radians).
 * @return Inclusive row bounds covering the band.
 * @details The leaf shapes' get_vertical_bounds all close with this; phi spans
 * [0,π] over (H+H_OFFSET-1) virtual rows.
 */
template <int H> inline Bounds phi_bounds_to_rows(float phi_min, float phi_max) {
  return phi_bounds_to_rows(phi_min, phi_max, H + hs::H_OFFSET, H);
}

/**
 * @brief Clamp an annular band's edge cosines to the visible scan range and
 * return its angular half-extents.
 * @param cos_outer Great-circle cosine of the band's far edge.
 * @param cos_inner Great-circle cosine of the band's near edge.
 * @param ny y-component of the band axis.
 * @param cos_phi Cosine of the row's polar angle.
 * @param denom Row scale factor R·sinφ.
 * @param angle_min Output: smaller angular half-extent (radians from azimuth).
 * @param angle_max Output: larger angular half-extent (radians from azimuth).
 * @return False if the band misses this row; true with angles written otherwise.
 * @details cos decreases with angle, so the larger cosine (cos_inner) yields the
 * smaller angle. Shared by the annular scanline emitters.
 */
inline bool annular_band_angles(float cos_outer, float cos_inner, float ny,
                                float cos_phi, float denom, float &angle_min,
                                float &angle_max) {
  float C_min = (cos_outer - ny * cos_phi) / denom;
  float C_max = (cos_inner - ny * cos_phi) / denom;
  float min_cos = std::max(-1.0f, C_min);
  float max_cos = std::min(1.0f, C_max);
  if (min_cos > max_cos)
    return false; // Empty row
  angle_min = fast_acos(max_cos);
  angle_max = fast_acos(min_cos);
  return true;
}

/**
 * @brief Emit the horizontal scanline interval(s) where a row crosses an
 * annular band, handling the two pole-wraparound degeneracies.
 * @tparam W Canvas width in columns.
 * @tparam OutputIt Sink type invoked as out(float start, float end).
 * @param cos_outer Great-circle cosine of the band's far edge.
 * @param cos_inner Great-circle cosine of the band's near edge.
 * @param ny y-component of the band axis.
 * @param cos_phi Cosine of the row's polar angle.
 * @param denom Row scale factor R·sinφ.
 * @param alpha_angle Azimuth of the band axis (radians).
 * @param out Sink accepting (float start, float end).
 * @details A band touching the near/far pole collapses its two arcs into a
 * single span. Shared by Ring and DistortedRing, whose annular scanline math is
 * otherwise byte-identical. Emits nothing for a missed row; the caller reports
 * the row handled either way.
 */
template <int W, typename OutputIt>
inline void emit_annular_band(float cos_outer, float cos_inner, float ny,
                              float cos_phi, float denom, float alpha_angle,
                              OutputIt out) {
  float angle_min, angle_max;
  if (!annular_band_angles(cos_outer, cos_inner, ny, cos_phi, denom, angle_min,
                           angle_max))
    return;

  float scale = W / (2.0f * PI_F);
  float safe_threshold = 2.0f * PI_F / W;

  if (angle_min <= safe_threshold) {
    out(floorf((alpha_angle - angle_max) * scale),
        ceilf((alpha_angle + angle_max) * scale));
  } else if (angle_max >= PI_F - safe_threshold) {
    out(floorf((alpha_angle + angle_min) * scale),
        ceilf((alpha_angle + 2 * PI_F - angle_min) * scale));
  } else {
    out(floorf((alpha_angle - angle_max) * scale),
        ceilf((alpha_angle - angle_min) * scale));
    out(floorf((alpha_angle + angle_min) * scale),
        ceilf((alpha_angle + angle_max) * scale));
  }
}

/**
 * @brief Calculates signed distance to a ring.
 * @details DistanceResult fields: dist = signed distance (negative inside);
 * t = normalized parameter (0-1) corresponding to angle/2PI; raw_dist = unsigned
 * distance to centerline.
 */
struct Ring {
  const Basis &basis; /**< Orientation frame (v = ring axis). */
  float radius;       /**< Ring radius as a fraction of the hemisphere. */
  float thickness;    /**< Half-width of the stroke (radians). */
  float phase;        /**< Azimuth phase offset (radians). */

  Vector normal, u, w; /**< Ring axis and the two in-plane basis vectors. */
  float nx, ny, nz;    /**< Components of the ring axis. */
  float target_angle, center_phi; /**< Centerline polar angle and axis colatitude. */
  float cos_max, cos_min, cos_target, inv_sin_target, sin_target; /**< Precomputed band trig. */

  float r_val;       /**< Horizontal projection length of the axis (for full-row check). */
  float alpha_angle; /**< Azimuth angle of the normal vector in the XZ plane. */
  static constexpr bool is_solid = false; /**< Ring renders as a stroke. */

  /**
   * @brief Builds a ring from its basis, radius, thickness, and phase.
   * @param b Orientation frame (v = ring axis).
   * @param r Ring radius as a fraction of the hemisphere.
   * @param th Half-width of the stroke (radians).
   * @param ph Azimuth phase offset (radians).
   */
  Ring(const Basis &b, float r, float th, float ph = 0)
      : basis(b), radius(r), thickness(th), phase(ph) {
    normal = basis.v;
    u = basis.u;
    w = basis.w;
    nx = normal.x;
    ny = normal.y;
    nz = normal.z;

    target_angle = radius * (PI_F / 2.0f);
    // Cold path (ctor, once per shape): full-precision trig, matching the
    // polygon family. The precomputed cos/sin constants below feed the
    // per-pixel reject and centerline distance, so exact values here cost
    // nothing in the hot loop while keeping the SDF accurate.
    center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));

    float ang_min = std::max(0.0f, target_angle - thickness);
    float ang_max = std::min(PI_F, target_angle + thickness);
    cos_max = cosf(ang_min);
    cos_min = cosf(ang_max);
    cos_target = cosf(target_angle);

    sin_target = sinf(target_angle);
    bool safe_approx = (target_angle > POLE_SAFE_MARGIN &&
                        target_angle < PI_F - POLE_SAFE_MARGIN);
    inv_sin_target = safe_approx ? (1.0f / sin_target) : 0.0f;

    // For get_horizontal_bounds
    r_val = sqrtf(nx * nx + nz * nz);
    alpha_angle = atan2f(nz, nx);
  }

  /**
   * @brief Maps the ring's latitude band to its inclusive scanline row range.
   * @tparam H Canvas height in rows.
   * @return Row bounds covering the ring plus its AA falloff.
   */
  template <int H> Bounds get_vertical_bounds() const {
    float a1 = center_phi - target_angle;
    float a2 = center_phi + target_angle;
    float phi_min = 0, phi_max = PI_F;

    if (a1 > 0) {
      float p1 = clamp_phi(a1);
      float p2 = clamp_phi(a2);
      phi_min = std::min(p1, p2);
    }
    if (a2 < PI_F) {
      float p1 = clamp_phi(a1);
      float p2 = clamp_phi(a2);
      phi_max = std::max(p1, p2);
    }

    float eff_th = 0.95f * thickness; // quintic_kernel(0.05) ≈ 0.001
    float f_phi_min = std::max(0.0f, phi_min - eff_th);
    float f_phi_max = std::min(PI_F, phi_max + eff_th);

    return phi_bounds_to_rows<H>(f_phi_min, f_phi_max);
  }

  /**
   * @brief Computes the horizontal scanline intervals for this shape at a given
   * y-coordinate.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The vertical pixel coordinate (row index).
   * @param out Output iterator or callback accepting (float start, float end).
   * @return True if intervals were found and reported; false requests a full scan.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    float cos_phi = TrigLUT<W, H>::cos_phi[y];
    float sin_phi = TrigLUT<W, H>::sin_phi[y];

    if (r_val < MIN_HORIZONTAL_PROJ)
      return false;

    float denom = r_val * sin_phi;
    if (std::abs(denom) < INTERVAL_DENOM_EPS)
      return false;

    float C_target = (cos_target - ny * cos_phi) / denom;
    float scale = W / (2.0f * PI_F);

    // Centerline-crossing: fast path when the annular band edges won't be
    // clamped to [-1,1]. Past the clamp boundary the first-order formula
    // over-estimates width (spikes), so we fall back to the exact annular
    // band arcs there.
    float clamp_bound = 1.0f - sin_target * thickness / denom;
    if (clamp_bound < 1.0f && C_target > -clamp_bound && C_target < clamp_bound &&
        sin_target > 1e-3f) {
      // Floor the radicand: at the upper edge of the guard band C_target nears 1
      // and 1-C^2 nears 0 (and can go slightly negative under -ffast-math, which
      // would yield sqrtf(NaN)). The floor both kills the NaN and caps half_width
      // so a grazing row cannot emit a runaway interval; genuine near-tangent
      // rows are handled by the annular fallback below.
      float sin_cross = sqrtf(std::max(1.0f - C_target * C_target, 1e-6f));
      float acos_C = fast_acos(C_target);
      float eff_th = 0.95f * thickness;
      float half_width = eff_th * sin_target / (denom * sin_cross);

      float hw_px = half_width * scale;
      float t1 = (alpha_angle - acos_C) * scale;
      float t2 = (alpha_angle + acos_C) * scale;

      out(floorf(t1 - hw_px), ceilf(t1 + hw_px));
      out(floorf(t2 - hw_px), ceilf(t2 + hw_px));
      return true;
    }

    // Annular band: exact intervals for near-tangent / out-of-range rows
    emit_annular_band<W>(cos_min, cos_max, ny, cos_phi, denom, alpha_angle, out);
    return true;
  }

  /**
   * @brief Computes signed distance to the ring.
   * @param p Point on sphere (normalized).
   * @return DistanceResult with dist, t, and raw_dist populated.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  /**
   * @brief Computes signed distance to the ring, writing into res.
   * @tparam ComputeUVs When true, also computes the azimuthal t parameter.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed distance, raw_dist = unsigned
   *        centerline distance, t = azimuth in [0,1) when ComputeUVs.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float d = dot(p, normal);
    if (d < cos_min || d > cos_max) {
      res = DistanceResult(100.0f, 0.0f, 100.0f, 0.0f, thickness);
      return;
    }

    float dist = 0;
    if (inv_sin_target != 0) {
      dist = std::abs(d - cos_target) * inv_sin_target;
    } else {
      float polar = fast_acos(hs::clamp(d, -1.0f, 1.0f));
      dist = std::abs(polar - target_angle);
    }

    float t = 0.0f;
    if constexpr (ComputeUVs) {
      float dot_u = dot(p, u);
      float dot_w = dot(p, w);
      float azimuth = fast_atan2(dot_w, dot_u);
      if (azimuth < 0)
        azimuth += 2 * PI_F;
      azimuth += phase;
      t = wrap_t(azimuth / (2 * PI_F));
    }

    res = DistanceResult(dist - thickness, t, dist, 0.0f, thickness);
  }
};

/**
 * @brief Calculates signed distance to a distorted ring.
 * @details DistanceResult fields: dist = signed distance minus thickness;
 * t = normalized parameter (0-1) corresponding to angle/2PI; raw_dist = unsigned
 * distance to centerline.
 */
struct DistortedRing {
  const Basis &basis;     /**< Orientation frame (v = ring axis). */
  float radius;           /**< Ring radius as a fraction of the hemisphere. */
  float thickness;        /**< Half-width of the stroke (radians). */
  ScalarFn shift_fn;      /**< Per-azimuth centerline shift, t in [0,1) -> radians. */
  float max_distortion;   /**< Maximum magnitude of shift_fn (radians). */
  float phase;            /**< Azimuth phase offset (radians). */

  Vector normal, u, w;            /**< Ring axis and the two in-plane basis vectors. */
  float nx, ny, nz;               /**< Components of the ring axis. */
  float target_angle, center_phi; /**< Centerline polar angle and axis colatitude. */
  float max_thickness;            /**< thickness + max_distortion (radians). */

  float r_val;                    /**< Horizontal projection length of the axis. */
  float alpha_angle;              /**< Azimuth of the normal in the XZ plane. */
  float cos_max_limit, cos_min_limit; /**< Cosines of the widened band edges. */
  static constexpr bool is_solid = false; /**< Distorted ring renders as a stroke. */

  /**
   * @brief Builds a distorted ring with a per-azimuth centerline shift.
   * @param b Orientation frame (v = ring axis).
   * @param r Ring radius as a fraction of the hemisphere.
   * @param th Half-width of the stroke (radians).
   * @param sf Per-azimuth centerline shift function, t in [0,1) -> radians.
   * @param md Maximum magnitude of sf (radians), used to widen the band.
   * @param ph Azimuth phase offset (radians).
   */
  DistortedRing(const Basis &b, float r, float th, ScalarFn sf, float md,
                float ph)
      : basis(b), radius(r), thickness(th), shift_fn(sf), max_distortion(md),
        phase(ph) {
    normal = basis.v;
    u = basis.u;
    w = basis.w;
    nx = normal.x;
    ny = normal.y;
    nz = normal.z;
    target_angle = radius * (PI_F / 2.0f);
    // Cold path (ctor): full-precision trig, matching Ring and the polygon
    // family. Runs once per shape, never per pixel.
    center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));
    max_thickness = thickness + max_distortion;

    r_val = sqrtf(nx * nx + nz * nz);
    alpha_angle = atan2f(nz, nx);

    float ang_min = std::max(0.0f, target_angle - max_thickness);
    float ang_max = std::min(PI_F, target_angle + max_thickness);
    cos_max_limit = cosf(ang_min);
    cos_min_limit = cosf(ang_max);
  }

  /**
   * @brief Maps the distorted ring's widened latitude band to its row range.
   * @tparam H Canvas height in rows.
   * @return Inclusive row bounds covering the ring plus distortion margin.
   */
  template <int H> Bounds get_vertical_bounds() const {
    float a1 = center_phi - target_angle;
    float a2 = center_phi + target_angle;
    float phi_min = 0, phi_max = PI_F;

    if (a1 > 0) {
      float p1 = clamp_phi(a1);
      float p2 = clamp_phi(a2);
      phi_min = std::min(p1, p2);
    }
    if (a2 < PI_F) {
      float p1 = clamp_phi(a1);
      float p2 = clamp_phi(a2);
      phi_max = std::max(p1, p2);
    }

    float margin = max_thickness + BOUNDS_MARGIN_WIDE;
    float f_phi_min = std::max(0.0f, phi_min - margin);
    float f_phi_max = std::min(PI_F, phi_max + margin);

    return phi_bounds_to_rows<H>(f_phi_min, f_phi_max);
  }

  /**
   * @brief Emits the widened annular-band intervals for one scanline row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    float cos_phi = TrigLUT<W, H>::cos_phi[y];
    float sin_phi = TrigLUT<W, H>::sin_phi[y];

    if (r_val < MIN_HORIZONTAL_PROJ)
      return false;

    float denom = r_val * sin_phi;
    if (std::abs(denom) < INTERVAL_DENOM_EPS)
      return false;

    emit_annular_band<W>(cos_min_limit, cos_max_limit, ny, cos_phi, denom,
                         alpha_angle, out);
    return true;
  }

  /**
   * @brief Computes signed distance to the distorted ring.
   * @param p Point on sphere (normalized).
   * @return DistanceResult with dist, t, and raw_dist populated.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  /**
   * @brief Computes signed distance to the distorted ring, writing into res.
   * @tparam ComputeUVs When true, also computes the azimuthal t parameter.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed distance minus thickness, raw_dist =
   *        unsigned centerline distance, t = azimuth in [0,1) when ComputeUVs.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float d = dot(p, normal);
    // Early reject: outside bounding annulus
    if (d < cos_min_limit || d > cos_max_limit) {
      res = DistanceResult(100.0f, 0.0f, 100.0f, 0.0f, thickness);
      return;
    }
    float polar = fast_acos(hs::clamp(d, -1.0f, 1.0f));

    float dot_u = dot(p, u);
    float dot_w = dot(p, w);
    float azimuth = fast_atan2(dot_w, dot_u);
    if (azimuth < 0)
      azimuth += 2 * PI_F;

    float t_norm = wrap_t((azimuth + phase) / (2 * PI_F));
    float shift = shift_fn(t_norm);

    if constexpr (!ComputeUVs)
      t_norm = 0.0f;

    float local_target = target_angle + shift;
    float dist = std::abs(polar - local_target);

    res = DistanceResult(dist - thickness, t_norm, dist, 0.0f, thickness);
  }
};

/**
 * @brief CSG Union operation (A + B), taking the minimum distance of two shapes.
 * @tparam A First child shape type.
 * @tparam B Second child shape type.
 */
template <typename A, typename B> struct Union {
  const A &a;         /**< First child shape. */
  const B &b;         /**< Second child shape. */
  float thickness;    /**< Max child thickness (drives AA falloff). */
  // A composite renders solid (1px silhouette AA) only if every child is solid.
  // If any child is a stroke, the result is a stroke so the rasterizer takes
  // the thickness-falloff AA path and each child keeps its own soft edge —
  // matching how that child looks drawn on its own.
  static constexpr bool is_solid = A::is_solid && B::is_solid; /**< Solid iff both children are. */

  /**
   * @brief Builds a union of two child shapes.
   * @param shape_a First child shape.
   * @param shape_b Second child shape.
   */
  Union(const A &shape_a, const B &shape_b)
      : a(shape_a), b(shape_b),
        thickness(std::max(shape_a.thickness, shape_b.thickness)) {}

  /**
   * @brief Row bounds spanning the union of the children's bands.
   * @tparam H Canvas height in rows.
   * @return Inclusive row bounds covering either child.
   */
  template <int H> Bounds get_vertical_bounds() const {
    auto b1 = a.template get_vertical_bounds<H>();
    auto b2 = b.template get_vertical_bounds<H>();
    return {std::min(b1.y_min, b2.y_min), std::max(b1.y_max, b2.y_max)};
  }

  /**
   * @brief Emits the merged union of both children's intervals for one row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan (whenever
   *         either child falls back).
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    StaticCircularBuffer<std::pair<float, float>, 32> merged;

    bool has_a = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { push_interval(merged, start, end); });

    bool has_b = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { push_interval(merged, start, end); });

    // If neither shape provides intervals, fall back to full scan
    if (!has_a && !has_b)
      return false;

    // If only one shape returned intervals and the other fell back,
    // we must do full scan (the fallback shape needs all pixels)
    if (!has_a || !has_b)
      return false;

    if (merged.is_empty())
      return true; // Both reported intervals but none produced

    merge_intervals(merged, out);
    return true;
  }

  /**
   * @brief Signed distance to the union (minimum of the two children).
   * @param p Point on sphere (normalized).
   * @return DistanceResult of the nearer child.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  /**
   * @brief Signed distance to the union, writing the nearer child into res.
   * @tparam ComputeUVs Forwarded to each child's distance().
   * @param p Point on sphere (normalized).
   * @param res Output result; the nearer child's full result is kept.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    a.template distance<ComputeUVs>(p, res);
    DistanceResult res_b;
    b.template distance<ComputeUVs>(p, res_b);
    if (res.dist < res_b.dist)
      return; // A is already closer, keep it
    res = res_b;
  }
};

/**
 * @brief Smooth CSG Union using a polynomial smooth minimum (Inigo Quilez smin).
 * @tparam A First child shape type.
 * @tparam B Second child shape type.
 * @details Shapes organically blend together within radius k (radians).
 */
template <typename A, typename B> struct SmoothUnion {
  const A &a;         /**< First child shape. */
  const B &b;         /**< Second child shape. */
  float k;            /**< Smoothing radius in radians (e.g. 0.1). */
  float thickness;    /**< Max child thickness (drives AA falloff). */
  // Solid only if every child is solid (see Union); a stroke child routes the
  // composite through the thickness-falloff AA path, so smin-blended strokes
  // keep their own soft edge instead of collapsing to a hard 1-px silhouette.
  static constexpr bool is_solid = A::is_solid && B::is_solid; /**< Solid iff both children are. */

  /**
   * @brief Builds a smooth union of two child shapes.
   * @param shape_a First child shape.
   * @param shape_b Second child shape.
   * @param smoothness Blend radius k in radians.
   */
  SmoothUnion(const A &shape_a, const B &shape_b, float smoothness)
      : a(shape_a), b(shape_b), k(smoothness),
        thickness(std::max(shape_a.thickness, shape_b.thickness)) {}

  /**
   * @brief Row bounds spanning both children's bands, padded by the blend radius.
   * @tparam H Canvas height in rows.
   * @return Inclusive row bounds expanded by k (converted to rows).
   */
  template <int H> Bounds get_vertical_bounds() const {
    constexpr int H_VIRT = H + hs::H_OFFSET;
    auto b1 = a.template get_vertical_bounds<H>();
    auto b2 = b.template get_vertical_bounds<H>();
    // Expand bounds by the blending radius (k, in radians) converted to rows:
    // phi spans [0,π] over (H_VIRT-1) rows, mirroring the horizontal pad (k→px).
    // The pad scales with k so large blends are not clipped top/bottom.
    int pad = std::max(1, static_cast<int>(ceilf(k * (H_VIRT - 1) / PI_F)));
    return {std::max(0, std::min(b1.y_min, b2.y_min) - pad),
            std::min(H - 1, std::max(b1.y_max, b2.y_max) + pad)};
  }

  /**
   * @brief Conservative union of the children's intervals, padded by k.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false (full scan) if either child
   *         falls back to full width.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    StaticCircularBuffer<std::pair<float, float>, 32> merged;
    float pad_px = k * W / (2 * PI_F);

    bool has_a = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          push_interval(merged, start - pad_px, end + pad_px);
        });

    bool has_b = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          push_interval(merged, start - pad_px, end + pad_px);
        });

    // If either child falls back to full-width, so must the blend
    if (!has_a || !has_b)
      return false;

    if (merged.is_empty())
      return true;

    merge_intervals(merged, out);
    return true;
  }

  /**
   * @brief Signed distance to the smooth union.
   * @param p Point on sphere (normalized).
   * @return DistanceResult of the nearer child with the smin blend applied.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  /**
   * @brief Signed distance to the smooth union, writing into res.
   * @tparam ComputeUVs Forwarded to each child's distance().
   * @param p Point on sphere (normalized).
   * @param res Output result; the nearer child's result with its dist reduced by
   *        the cubic smin blend term.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    a.template distance<ComputeUVs>(p, res);
    DistanceResult res_b;
    b.template distance<ComputeUVs>(p, res_b);

    // Polynomial smooth min (cubic)
    float h = std::max(k - std::abs(res.dist - res_b.dist), 0.0f) / k;
    float m = h * h * h * k * (1.0f / 6.0f);

    if (res.dist < res_b.dist) {
      res.dist -= m;
    } else {
      res_b.dist -= m;
      res = res_b;
    }
  }
};

/**
 * @brief CSG Subtraction operation (A - B), computing max(A, -B).
 * @tparam A Minuend shape type (the shape carved into).
 * @tparam B Subtrahend shape type (the shape removed).
 */
template <typename A, typename B> struct Subtract {
  const A &a;      /**< Minuend shape. */
  const B &b;      /**< Subtrahend shape (removed from A). */
  float thickness; /**< Inherited from the minuend A. */
  // Solid only if every child is solid (see Union); a stroke child routes the
  // composite through the thickness-falloff AA path.
  static constexpr bool is_solid = A::is_solid && B::is_solid; /**< Solid iff both children are. */

  /**
   * @brief Builds a subtraction (shape_a minus shape_b).
   * @param shape_a Minuend shape.
   * @param shape_b Subtrahend shape.
   */
  Subtract(const A &shape_a, const B &shape_b)
      : a(shape_a), b(shape_b), thickness(shape_a.thickness) {}

  /**
   * @brief Row bounds of the minuend (subtraction never grows the band).
   * @tparam H Canvas height in rows.
   * @return The minuend's inclusive row bounds.
   */
  template <int H> Bounds get_vertical_bounds() const {
    return a.template get_vertical_bounds<H>();
  }

  /**
   * @brief Emits A's intervals with B's overlapping intervals removed.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan whenever A
   *         or B falls back (B's fallback cannot be treated as full coverage).
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    StaticCircularBuffer<std::pair<float, float>, 32> intervals_a;
    StaticCircularBuffer<std::pair<float, float>, 32> intervals_b;

    // The set-difference loop below and scan_region's coalescer both assume
    // intervals arrive in non-decreasing start order; a multi-interval child
    // (a nested CSG) can emit them out of order, which without sorting yields a
    // wrong set difference AND unsorted output the coalescer silently drops.
    // sort_intervals_by_start() (above) is shared with Union/SmoothUnion.
    bool has_a = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { push_interval(intervals_a, start, end); });

    // If A falls back to full-width, we can't restrict
    if (!has_a)
      return false;

    // If A produced no intervals on this row, nothing to subtract from
    if (intervals_a.is_empty())
      return true;

    // A's pieces are emitted in A-interval order, so A must be start-sorted for
    // the output (and the pass-through below) to stay start-ordered.
    sort_intervals_by_start(intervals_a);

    // A stroke subtrahend (e.g. Ring, Line) emits intervals describing only its
    // thin drawn band, and adjacent edge-bands can coalesce into a single chord
    // interval that spans the stroke's hollow interior. Treating those as
    // removable would carve A's interior between the stroke edges, and the
    // skipped columns never get re-filled (scan_region only refines columns it
    // scans). So for a non-solid B, do NOT do interval subtraction: pass A's
    // intervals through unchanged and let scan_region's per-pixel max(A, -B)
    // carve exactly the stroke band and leave A's interior intact.
    if constexpr (!B::is_solid) {
      for (size_t i = 0; i < intervals_a.size(); ++i)
        out(intervals_a[i].first, intervals_a[i].second);
      return true;
    }

    bool has_b = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { push_interval(intervals_b, start, end); });

    // B fell back to a full-row scan: it could not produce intervals, NOT that
    // it covers the row. Without B's intervals we cannot compute the set
    // difference, so request a full-row scan (return false) exactly as
    // Union/SmoothUnion do — scan_region then evaluates distance() = max(A, -B)
    // per pixel, which is correct. Returning true + emitting nothing would skip
    // the row and silently erase all of A (e.g. when B is a near-vertical Ring
    // or any AngularRepeat, both of which fall back unconditionally).
    if (!has_b)
      return false;

    // B produced no intervals on this row: it removes nothing, so pass A's
    // intervals through unchanged.
    if (intervals_b.is_empty()) {
      for (size_t i = 0; i < intervals_a.size(); ++i)
        out(intervals_a[i].first, intervals_a[i].second);
      return true;
    }

    // The per-A shrink/split loop walks B left-to-right, so B must be sorted.
    sort_intervals_by_start(intervals_b);

    // Set difference: A minus B
    // For each A interval, subtract all overlapping B intervals
    for (size_t ai = 0; ai < intervals_a.size(); ++ai) {
      float cur_start = intervals_a[ai].first;
      float cur_end = intervals_a[ai].second;

      for (size_t bi = 0; bi < intervals_b.size(); ++bi) {
        float bs = intervals_b[bi].first;
        float be = intervals_b[bi].second;

        // No overlap
        if (be <= cur_start || bs >= cur_end)
          continue;

        // B covers the left portion: emit nothing, shrink from left
        if (bs <= cur_start) {
          cur_start = be;
          if (cur_start >= cur_end)
            break;
          continue;
        }

        // B covers the right portion: shrink from right
        if (be >= cur_end) {
          cur_end = bs;
          if (cur_start >= cur_end)
            break;
          continue;
        }

        // B splits A: emit left piece, continue with right piece
        out(cur_start, bs);
        cur_start = be;
      }

      if (cur_start < cur_end)
        out(cur_start, cur_end);
    }
    return true;
  }

  /**
   * @brief Signed distance to the subtraction (A - B).
   * @param p Point on sphere (normalized).
   * @return DistanceResult of max(A, -B).
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  /**
   * @brief Signed distance to the subtraction, writing into res.
   * @tparam ComputeUVs Forwarded to each child's distance().
   * @param p Point on sphere (normalized).
   * @param res Output result; max(A, -B) with B's distance negated when it wins.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    a.template distance<ComputeUVs>(p, res);
    DistanceResult res_b;
    b.template distance<ComputeUVs>(p, res_b);
    // Max(A, -B)
    if (-res_b.dist > res.dist) {
      res_b.dist = -res_b.dist;
      res = res_b;
    }
  }
};

/**
 * @brief CSG Intersection operation (A & B), taking the maximum distance.
 * @tparam A First child shape type.
 * @tparam B Second child shape type.
 */
template <typename A, typename B> struct Intersection {
  const A &a;      /**< First child shape. */
  const B &b;      /**< Second child shape. */
  float thickness; /**< Min child thickness (drives AA falloff). */
  // Solid only if every child is solid (see Union); a stroke child routes the
  // composite through the thickness-falloff AA path.
  static constexpr bool is_solid = A::is_solid && B::is_solid; /**< Solid iff both children are. */

  /**
   * @brief Builds an intersection of two child shapes.
   * @param shape_a First child shape.
   * @param shape_b Second child shape.
   */
  Intersection(const A &shape_a, const B &shape_b)
      : a(shape_a), b(shape_b),
        thickness(std::min(shape_a.thickness, shape_b.thickness)) {}

  /**
   * @brief Row bounds of the overlap of the children's bands.
   * @tparam H Canvas height in rows.
   * @return Inclusive row bounds covered by both children.
   */
  template <int H> Bounds get_vertical_bounds() const {
    auto b1 = a.template get_vertical_bounds<H>();
    auto b2 = b.template get_vertical_bounds<H>();
    return {std::max(b1.y_min, b2.y_min), std::min(b1.y_max, b2.y_max)};
  }

  /**
   * @brief Emits the per-row intersection of both children's intervals.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false (full scan) only when both
   *         children fall back to full width.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    StaticCircularBuffer<std::pair<float, float>, 32> intervals_a;
    StaticCircularBuffer<std::pair<float, float>, 32> intervals_b;

    bool has_a = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { push_interval(intervals_a, start, end); });

    bool has_b = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { push_interval(intervals_b, start, end); });

    // A child that returns false wants a full-width scan; intersecting that
    // full width with the other child is just the other child's intervals,
    // which we already collected above. Replay the buffer instead of
    // re-invoking the child — avoids doubling a nested CSG child's per-row
    // interval work, and is robust to a child whose emission is not identical
    // across calls.
    if (!has_a) {
      for (size_t i = 0; i < intervals_b.size(); ++i)
        out(intervals_b[i].first, intervals_b[i].second);
      return has_b; // both fell back (has_b == false) -> full scan
    }
    if (!has_b) {
      for (size_t i = 0; i < intervals_a.size(); ++i)
        out(intervals_a[i].first, intervals_a[i].second);
      return true; // has_a is true here
    }

    if (intervals_a.is_empty() || intervals_b.is_empty())
      return true;

    // The merge-sweep below assumes both lists are start-sorted (it advances by
    // comparing interval ends); a multi-interval child (a nested CSG) can emit
    // them out of order, so sort first — same as Subtract. Single-interval
    // children are already trivially sorted.
    sort_intervals_by_start(intervals_a);
    sort_intervals_by_start(intervals_b);

    // Intersect sorted intervals
    size_t idx_a = 0;
    size_t idx_b = 0;

    while (idx_a < intervals_a.size() && idx_b < intervals_b.size()) {
      auto iv_a = intervals_a[idx_a];
      auto iv_b = intervals_b[idx_b];

      float start = std::max(iv_a.first, iv_b.first);
      float end = std::min(iv_a.second, iv_b.second);

      if (start < end) {
        out(start, end);
      }

      if (iv_a.second < iv_b.second) {
        idx_a++;
      } else {
        idx_b++;
      }
    }

    return true;
  }

  /**
   * @brief Signed distance to the intersection (maximum of the two children).
   * @param p Point on sphere (normalized).
   * @return DistanceResult of the farther child.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  /**
   * @brief Signed distance to the intersection, writing the farther child.
   * @tparam ComputeUVs Forwarded to each child's distance().
   * @param p Point on sphere (normalized).
   * @param res Output result; the farther child's full result is kept.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    a.template distance<ComputeUVs>(p, res);
    DistanceResult res_b;
    b.template distance<ComputeUVs>(p, res_b);
    // Max(A, B)
    if (res.dist > res_b.dist)
      return;
    res = res_b;
  }
};

/**
 * @brief Angular domain repetition modifier.
 * @tparam Shape The child shape type being repeated.
 * @details Folds the azimuthal angle to create N copies of a shape around an
 * arbitrary axis for constant cost (a single distance evaluation).
 */
template <typename Shape> struct AngularRepeat {
  const Shape &shape; /**< Child shape being repeated. */
  Vector axis, u, w;  /**< Rotation axis and the derived perpendicular plane (u, w). */
  int repetitions;    /**< Number of copies around the axis. */
  float thickness;    /**< Inherited from the child shape. */
  static constexpr bool is_solid = Shape::is_solid; /**< Matches the child's solidity. */

  /**
   * @brief Repeats the shape around an arbitrary axis.
   * @param s Child shape to repeat.
   * @param reps Number of copies (must be > 0).
   * @param ax Rotation axis (unit length).
   */
  AngularRepeat(const Shape &s, int reps, const Vector &ax)
      : shape(s), axis(ax), repetitions(reps), thickness(s.thickness) {
    HS_CHECK(reps > 0); // sector folding divides by repetitions
    // Derive perpendicular plane via Gram-Schmidt
    Vector ref = (fabsf(ax.y) < 0.9f) ? Vector(0, 1, 0) : Vector(1, 0, 0);
    float d = ref.x * ax.x + ref.y * ax.y + ref.z * ax.z;
    u = Vector(ref.x - d * ax.x, ref.y - d * ax.y, ref.z - d * ax.z);
    float len = u.length();
    if (len > TOLERANCE)
      u = u * (1.0f / len);
    w = Vector(ax.y * u.z - ax.z * u.y, ax.z * u.x - ax.x * u.z,
               ax.x * u.y - ax.y * u.x);
  }

  /**
   * @brief Repeats the shape around the Y-axis.
   * @param s Child shape to repeat.
   * @param reps Number of copies (must be > 0).
   */
  AngularRepeat(const Shape &s, int reps)
      : AngularRepeat(s, reps, Vector(0, 1, 0)) {}

  /**
   * @brief Row bounds for the repeated shape.
   * @tparam H Canvas height in rows.
   * @return The child's band for a Y-axis fold; the full canvas otherwise.
   */
  template <int H> Bounds get_vertical_bounds() const {
    // Only a fold about the Y axis preserves latitude (the axis component, which
    // becomes a copy's invariant, is then the row axis), so the un-repeated
    // child's latitude band still bounds every copy. For any other axis the
    // copies sweep through latitudes the child never occupies; forwarding the
    // child's narrow band would row-clip those copies out of the scan and drop
    // them silently. Fall back to the full canvas there — conservative, and the
    // common Y-axis path keeps its tight bound. (axis is unit length, as the
    // distance() projection basis requires, so axis.y near ±1 is the Y fold.)
    if (fabsf(axis.y) < 1.0f - TOLERANCE)
      return {0, H - 1};
    return shape.template get_vertical_bounds<H>();
  }

  /**
   * @brief Always requests a full-width scan (copies cover the full azimuth).
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type (unused).
   * @return Always false (full scan).
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int, OutputIt) const {
    // Full scan: repetitions cover the full azimuth
    return false;
  }

  /**
   * @brief Signed distance to the nearest repeated copy.
   * @param p Point on sphere (normalized).
   * @return DistanceResult of the child evaluated in the folded sector.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  /**
   * @brief Folds p into one sector and evaluates the child, writing into res.
   * @tparam ComputeUVs Forwarded to the child's distance().
   * @param p Point on sphere (normalized).
   * @param res Output result of the child at the folded point.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    // Project p into local coordinate system
    float local_u = p.x * u.x + p.y * u.y + p.z * u.z;
    float local_v = p.x * axis.x + p.y * axis.y + p.z * axis.z;
    float local_w = p.x * w.x + p.y * w.y + p.z * w.z;

    // Fold angle in the u-w plane
    float theta = fast_atan2(local_w, local_u);
    if (theta < 0)
      theta += 2 * PI_F;

    float sector = (2 * PI_F) / repetitions;
    float folded_theta = fmodf(theta + sector / 2.0f, sector) - sector / 2.0f;

    // Reconstruct folded local coordinates (preserving axis component)
    float r = sqrtf(local_u * local_u + local_w * local_w);
    float fu = r * fast_cosf(folded_theta);
    float fw = r * fast_sinf(folded_theta);

    // Project back to world space
    Vector folded_p(fu * u.x + local_v * axis.x + fw * w.x,
                    fu * u.y + local_v * axis.y + fw * w.y,
                    fu * u.z + local_v * axis.z + fw * w.z);

    shape.template distance<ComputeUVs>(folded_p, res);
  }
};

/**
 * @brief Scratch buffer for Face computations to avoid allocations.
 */
struct FaceScratchBuffer {
  static constexpr int MAX_VERTS = 64; /**< Maximum vertices per face. */
  static constexpr int LUT_N = 32;     /**< Distance-LUT grid resolution per axis. */
  std::array<Vector, MAX_VERTS + 1> poly_2d; /**< Projected 2D polygon (+1 entry to avoid modulo). */
  std::array<Vector, MAX_VERTS> edge_vectors;       /**< Per-edge 2D vectors. */
  std::array<float, MAX_VERTS> edge_lengths_sq;     /**< Per-edge squared lengths. */
  std::array<Vector, MAX_VERTS> planes;             /**< Per-edge great-circle normals. */
  std::array<std::pair<float, float>, 4> intervals; /**< Azimuth coverage intervals. */
  std::array<float, MAX_VERTS> thetas;              /**< Per-vertex azimuth angles. */
  std::array<float, MAX_VERTS> inv_edge_lengths_sq; /**< Reciprocal squared edge lengths. */
  std::array<float, MAX_VERTS> inv_edge_j;          /**< Reciprocal of each edge's y-component. */
  std::array<Vector, MAX_VERTS + 1> verts_3d;       /**< 3D vertices (+1 wrap entry). */
  std::array<Vector, MAX_VERTS> edge_normals;       /**< Per-edge normalized 3D normals. */
  std::array<float, LUT_N * LUT_N> dist_lut;        /**< Precomputed signed-distance LUT. */

  /**
   * @brief Packed per-edge data for the cache-friendly distance() fallback.
   */
  struct EdgePacked {
    float vx, vy, ex, ey, inv_len_sq, inv_ej, next_vy, _pad; /**< Edge origin, vector, reciprocals, next-vertex y, padding. */
  };
  std::array<EdgePacked, MAX_VERTS> packed_edges; /**< Packed per-edge data. */
};

/**
 * @brief Represents a planar face for SDF rendering.
 * @details Computes the 2D projection and vertical/horizontal bounds used to
 * accelerate rasterization.
 */
struct Face {
  Vector center;                    /**< Normalized face centroid (projection axis). */
  Vector basis_v, basis_u, basis_w; /**< Local tangent frame (v = center). */
  int count;                        /**< Vertex/edge count. */
  float thickness;                  /**< Edge half-width (radians). */
  float size;                       /**< Inradius metric for AA normalization. */
  float max_r2 = 0.0f;              /**< Max squared 2D vertex radius. */
  float radius = 0.0f;              /**< Circumradius in the 2D projection. */
  float max_dist = 0.0f;            /**< Cull radius (circumradius plus margin). */
  float max_dist_sq = 0.0f;         /**< Squared cull radius. */

  std::span<Vector> poly_2d;            /**< Projected 2D polygon (+1 wrap entry). */
  std::span<Vector> edge_vectors;       /**< Per-edge 2D vectors. */
  std::span<float> edge_lengths_sq;     /**< Per-edge squared lengths. */
  std::span<Vector> planes;             /**< Per-edge great-circle normals. */
  std::span<float> inv_edge_lengths_sq; /**< Reciprocal squared edge lengths. */
  std::span<float> inv_edge_j;          /**< Reciprocal of each edge's y-component. */
  std::span<Vector> verts_3d;           /**< 3D vertices (+1 wrap entry). */
  std::span<Vector> edge_normals;       /**< Per-edge normalized 3D normals. */

  int y_min, y_max;                            /**< Inclusive vertical row bounds. */
  std::span<std::pair<float, float>> intervals; /**< Azimuth coverage intervals (radians). */
  bool full_width;                             /**< True when the face spans all columns. */
  static constexpr bool is_solid = true;       /**< Face renders as a filled region. */

  // Packed edge data + LUT for distance computation
  using EdgePacked = FaceScratchBuffer::EdgePacked; /**< Packed per-edge record type. */
  std::span<EdgePacked> packed_edges;               /**< Packed per-edge data. */
  static constexpr int LUT_N = FaceScratchBuffer::LUT_N; /**< Distance-LUT resolution per axis. */
  const float *dist_lut = nullptr;     /**< Pointer into the scratch distance LUT. */
  float lut_cx = 0.0f, lut_cy = 0.0f;  /**< LUT bounding-box center. */
  float lut_Rx = 0.0f, lut_Ry = 0.0f;  /**< LUT bounding-box half-extents. */
  float lut_inv_step_x = 0.0f;         /**< Reciprocal LUT cell width. */
  float lut_inv_step_y = 0.0f;         /**< Reciprocal LUT cell height. */
  float lut_safe_dist = 0.0f;          /**< Per-face cell diagonal (sign-pure bound). */

  /**
   * @brief Builds a face's projection, bounds, edge data, and distance LUT.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices from the pool.
   * @param th Edge half-width (radians).
   * @param scratch Reusable scratch storage backing the spans.
   * @param h_virt Virtual row count (height plus pole offset).
   * @param height Canvas height in rows.
   */
  Face(std::span<const Vector> vertices, std::span<const uint16_t> indices,
       float th, FaceScratchBuffer &scratch, int h_virt, int height)
      : thickness(th), full_width(true) {

    // Early vertical exit: a face whose latitude band (plus thickness margin)
    // maps to an empty row range can never be rasterized.
    float min_phi_check, max_phi_check;
    if (compute_phi_extent(vertices, indices, h_virt, height, min_phi_check,
                           max_phi_check)) {
      y_min = 1;
      y_max = 0;
      return;
    }

    count = indices.size();
    // A face with more vertices than the scratch budget would build wrong
    // geometry from a truncated index list; an empty index list (count == 0)
    // would build a degenerate zero-edge polygon. Trap both bounds instead of
    // silently masking.
    HS_CHECK(count > 0 && count <= FaceScratchBuffer::MAX_VERTS,
             "Face: vertex count must be in (0, MAX_VERTS]");

    setup_frame_and_polygon(vertices, indices, scratch);
    compute_inradius(scratch);

    // Vertical bounds: large faces need the full arc-extrema + pole analysis;
    // small faces use the cheaper vertex-only phi span.
    int planes_count = 0;
    if (max_phi_check - min_phi_check < FAST_BOUNDS_PHI_THRESHOLD) {
      compute_fast_bounds(scratch, vertices, indices, count, thickness,
                          min_phi_check, max_phi_check, h_virt, height, y_min,
                          y_max);
    } else {
      planes_count =
          compute_full_bounds(scratch, vertices, indices, count, center,
                              thickness, h_virt, height, y_min, y_max);
    }

    edge_vectors = std::span<Vector>(scratch.edge_vectors.data(), count);
    edge_lengths_sq = std::span<float>(scratch.edge_lengths_sq.data(), count);
    inv_edge_lengths_sq = std::span<float>(scratch.inv_edge_lengths_sq.data(), count);
    inv_edge_j = std::span<float>(scratch.inv_edge_j.data(), count);
    planes = std::span<Vector>(scratch.planes.data(), planes_count);

    pack_edges(scratch);
    build_distance_lut(scratch);
    compute_azimuth_intervals(scratch);
    apply_pole_containment(height);
  }

  // ---------------------------------------------------------------------------
  // Construction phases. Each is a single-call-site helper factored out of the
  // constructor for readability; all are force-inlined so the per-frame
  // Mesh::draw face-setup cost stays free of call overhead.
  // ---------------------------------------------------------------------------

  /**
   * @brief Latitude-band reject for the face.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices.
   * @param h_virt Virtual row count (height plus pole offset).
   * @param height Canvas height in rows.
   * @param min_phi_check Output: minimum polar angle of the face (radians).
   * @param max_phi_check Output: maximum polar angle of the face (radians).
   * @return True when the phi extent plus thickness margin maps to an empty
   *         canvas-row range; the raw phi extent is always written out for reuse.
   */
  __attribute__((always_inline)) bool
  compute_phi_extent(std::span<const Vector> vertices,
                     std::span<const uint16_t> indices, int h_virt, int height,
                     float &min_phi_check, float &max_phi_check) const {
    float min_y_val = 2.0f;
    float max_y_val = -2.0f;

    for (int idx : indices) {
      float y = vertices[idx].y;
      if (y < min_y_val)
        min_y_val = y;
      if (y > max_y_val)
        max_y_val = y;
    }

    min_phi_check = fast_acos(hs::clamp(max_y_val, -1.0f, 1.0f));
    max_phi_check = fast_acos(hs::clamp(min_y_val, -1.0f, 1.0f));
    float margin_check = thickness + BOUNDS_MARGIN;

    Bounds rows = phi_bounds_to_rows(min_phi_check - margin_check,
                                     max_phi_check + margin_check, h_virt,
                                     height);

    return rows.y_min > rows.y_max;
  }

  /**
   * @brief Builds the local tangent frame, gnomonic 2D projection, and 3D arrays.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices.
   * @param scratch Scratch storage receiving poly_2d, verts_3d, edge_normals.
   * @details Sets basis_u/v/w, poly_2d (with circumradius), and the 3D
   * vertex/edge-normal arrays.
   */
  __attribute__((always_inline)) void
  setup_frame_and_polygon(std::span<const Vector> vertices,
                          std::span<const uint16_t> indices,
                          FaceScratchBuffer &scratch) {
    center = Vector(0, 0, 0);
    for (int i = 0; i < count; ++i)
      center = center + vertices[indices[i]];
    center.normalize();

    basis_v = center;
    Vector ref =
        (std::abs(dot(center, X_AXIS)) > math::COS_AXIS_PARALLEL) ? Y_AXIS
                                                                  : X_AXIS;
    basis_u = cross(center, ref).normalized();
    basis_w = cross(center, basis_u).normalized();

    max_r2 = 0.0f;
    for (int i = 0; i < count; ++i) {
      const Vector &v = vertices[indices[i]];
      // Gnomonic projection divides by d = cos(angle from the face center); it
      // is singular as a vertex approaches the center's antipode (d -> 0),
      // yielding ±inf/NaN coordinates. Registered solids never produce an
      // in-face vertex that far from its own face center, so this cannot fire on
      // shipped geometry — but guard it for parity with the plot-path planar
      // projection (plot.h COS_PLANAR_ANTIPODE) and the fast_atan2 origin guard.
      // Clamp d away from zero, preserving its sign so the projection stays on
      // the correct side.
      float d = dot(v, basis_v);
      if (fabsf(d) < math::TOLERANCE)
        d = copysignf(math::TOLERANCE, d);
      float px = dot(v, basis_u) / d;
      float py = dot(v, basis_w) / d;

      scratch.poly_2d[i] = Vector(px, py, 0);

      float r2 = px * px + py * py;
      if (r2 > max_r2)
        max_r2 = r2;
    }
    radius = sqrtf(max_r2);
    max_dist = radius + 0.1f;
    max_dist_sq = max_dist * max_dist;

    scratch.poly_2d[count] = scratch.poly_2d[0];
    poly_2d = std::span<Vector>(scratch.poly_2d.data(), count + 1);

    // Store 3D vertices and edge normals for angular distance computation
    for (int i = 0; i < count; ++i)
      scratch.verts_3d[i] = vertices[indices[i]];
    scratch.verts_3d[count] = scratch.verts_3d[0];
    verts_3d = std::span<Vector>(scratch.verts_3d.data(), count + 1);

    for (int i = 0; i < count; ++i) {
      Vector n = cross(scratch.verts_3d[i], scratch.verts_3d[i + 1]);
      float len = n.magnitude();
      scratch.edge_normals[i] =
          (len > 1e-9f) ? n * (1.0f / len) : Vector(0, 0, 0);
    }
    edge_normals = std::span<Vector>(scratch.edge_normals.data(), count);
  }

  /**
   * @brief Computes the face "size" (inradius) from the projected polygon.
   * @param scratch Scratch storage holding poly_2d.
   * @details size = minimum distance from the projected centroid to any edge,
   * floored to a fraction of the circumradius for degenerate slivers.
   */
  __attribute__((always_inline)) void
  compute_inradius(FaceScratchBuffer &scratch) {
    float min_edge_dist = 1e9f;
    for (int i = 0; i < count; ++i) {
      Vector v1 = scratch.poly_2d[i];
      Vector v2 = scratch.poly_2d[(i + 1) % count];

      Vector edge = v2 - v1;
      float edge_len_sq = dot(edge, edge);
      float t = 0.0f;
      if (edge_len_sq > 1e-9f) {
        t = dot(-v1, edge) / edge_len_sq;
        t = std::max(0.0f, std::min(1.0f, t));
      }
      Vector closest = v1 + edge * t;
      float d_line = closest.magnitude();

      if (d_line < min_edge_dist)
        min_edge_dist = d_line;
    }
    size = (min_edge_dist > 1e8f) ? 1.0f : min_edge_dist;

    if (size < radius * MIN_SIZE_RADIUS_RATIO)
      size = radius * MIN_SIZE_RADIUS_RATIO;
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
      ep.next_vy = poly_2d[i + 1].y;
    }
    packed_edges = std::span<EdgePacked>(scratch.packed_edges.data(), count);
  }

  /**
   * @brief Precomputes the anisotropic signed-distance LUT over the bounding box.
   * @param scratch Scratch storage receiving dist_lut.
   */
  __attribute__((always_inline)) void
  build_distance_lut(FaceScratchBuffer &scratch) {
    float bb_min_x = FLT_MAX, bb_max_x = -FLT_MAX;
    float bb_min_y = FLT_MAX, bb_max_y = -FLT_MAX;
    for (int i = 0; i < count; ++i) {
      float vx = poly_2d[i].x, vy = poly_2d[i].y;
      if (vx < bb_min_x)
        bb_min_x = vx;
      if (vx > bb_max_x)
        bb_max_x = vx;
      if (vy < bb_min_y)
        bb_min_y = vy;
      if (vy > bb_max_y)
        bb_max_y = vy;
    }
    float margin = 0.1f;
    lut_cx = (bb_min_x + bb_max_x) * 0.5f;
    lut_cy = (bb_min_y + bb_max_y) * 0.5f;
    lut_Rx = (bb_max_x - bb_min_x) * 0.5f + margin;
    lut_Ry = (bb_max_y - bb_min_y) * 0.5f + margin;
    if (lut_Rx < 0.01f)
      lut_Rx = 0.01f;
    if (lut_Ry < 0.01f)
      lut_Ry = 0.01f;
    lut_inv_step_x = (LUT_N - 1) / (2.0f * lut_Rx);
    lut_inv_step_y = (LUT_N - 1) / (2.0f * lut_Ry);
    dist_lut = scratch.dist_lut.data();
    float step_x = (2.0f * lut_Rx) / (LUT_N - 1);
    float step_y = (2.0f * lut_Ry) / (LUT_N - 1);
    // Bilinear interpolation is trusted only when no zero-isocontour can cross
    // the cell. The plane SDF is 1-Lipschitz in (px,py), so a zero anywhere in a
    // cell forces every corner to lie within one cell diagonal of it; requiring
    // the smallest corner magnitude to exceed the diagonal therefore guarantees
    // a sign-pure cell. This is a true per-face bound: it scales with the face's
    // own LUT cell size, so it stays correct (sign-pure edges) for large faces.
    lut_safe_dist = sqrtf(step_x * step_x + step_y * step_y);
    for (int gy = 0; gy < LUT_N; ++gy) {
      float qy = (lut_cy - lut_Ry) + gy * step_y;
      for (int gx = 0; gx < LUT_N; ++gx) {
        float qx = (lut_cx - lut_Rx) + gx * step_x;
        float d_sq = FLT_MAX;
        bool is_inside = false;
        for (int i = 0; i < count; ++i) {
          const auto &ep = packed_edges[i];
          float wx = qx - ep.vx, wy = qy - ep.vy;
          float t = (wx * ep.ex + wy * ep.ey) * ep.inv_len_sq;
          float cv = hs::clamp(t, 0.0f, 1.0f);
          float bx = wx - ep.ex * cv, by = wy - ep.ey * cv;
          float dsq = bx * bx + by * by;
          if (dsq < d_sq)
            d_sq = dsq;
          if ((ep.vy > qy) != (ep.next_vy > qy)) {
            float ix = ep.vx + (qy - ep.vy) * ep.ex * ep.inv_ej;
            if (qx < ix)
              is_inside = !is_inside;
          }
        }
        scratch.dist_lut[gy * LUT_N + gx] =
            (is_inside ? -1.0f : 1.0f) * sqrtf(d_sq);
      }
    }
  }

  /**
   * @brief Computes the face's azimuth coverage intervals.
   * @param scratch Scratch storage holding thetas and receiving intervals.
   * @details Finds the largest angular gap between vertices; if it exceeds pi the
   * face does not wrap, so the complementary horizontal interval(s) are emitted.
   * Otherwise the face spans the full width.
   */
  __attribute__((always_inline)) void
  compute_azimuth_intervals(FaceScratchBuffer &scratch) {
    std::sort(scratch.thetas.begin(), scratch.thetas.begin() + count);
    float max_gap = 0;
    float gap_start = 0;
    for (int i = 0; i < count; ++i) {
      float next = (i + 1 < count) ? scratch.thetas[i + 1]
                                   : (scratch.thetas[0] + 2 * PI_F);
      float diff = next - scratch.thetas[i];
      if (diff > max_gap) {
        max_gap = diff;
        gap_start = scratch.thetas[i];
      }
    }

    int interval_count = 0;
    if (max_gap > PI_F) {
      full_width = false;
      float start_t = fmodf(gap_start + max_gap, 2 * PI_F);
      float end_t = gap_start;

      if (start_t <= end_t) {
        scratch.intervals[interval_count++] = {start_t, end_t};
      } else {
        scratch.intervals[interval_count++] = {0.0f, end_t};
        scratch.intervals[interval_count++] = {start_t, 2 * PI_F};
      }
    } else {
      full_width = true;
    }
    intervals = std::span<std::pair<float, float>>(scratch.intervals.data(),
                                                   interval_count);
  }

  /**
   * @brief Extends the vertical bounds when the face encircles a pole.
   * @param height Canvas height in rows.
   * @details If the face encircles the north or south pole, the vertical bounds
   * are extended to it and full-width azimuth coverage is forced.
   */
  __attribute__((always_inline)) void apply_pole_containment(int height) {
    if (center.y > 0.01f) {
      float inv_c = 1.0f / center.y;
      float ppx = basis_u.y * inv_c;
      float ppy = basis_w.y * inv_c;
      bool pole_inside = false;
      for (int i = 0; i < count; ++i) {
        if ((poly_2d[i].y > ppy) != (poly_2d[i + 1].y > ppy)) {
          float ix = poly_2d[i].x +
                     (ppy - poly_2d[i].y) * edge_vectors[i].x * inv_edge_j[i];
          if (ppx < ix)
            pole_inside = !pole_inside;
        }
      }
      if (pole_inside) {
        y_min = 0;
        full_width = true;
      }
    }
    // South pole (0, -1, 0)
    if (center.y < -0.01f) {
      float cos_sp = -center.y;
      float inv_c = 1.0f / cos_sp;
      float ppx = -basis_u.y * inv_c;
      float ppy = -basis_w.y * inv_c;
      bool pole_inside = false;
      for (int i = 0; i < count; ++i) {
        if ((poly_2d[i].y > ppy) != (poly_2d[i + 1].y > ppy)) {
          float ix = poly_2d[i].x +
                     (ppy - poly_2d[i].y) * edge_vectors[i].x * inv_edge_j[i];
          if (ppx < ix)
            pole_inside = !pole_inside;
        }
      }
      if (pole_inside) {
        y_max = height - 1;
        full_width = true;
      }
    }
  }

  /**
   * @brief Fast-path vertical bounds using vertex-based phi for small faces.
   * @param scratch Scratch storage receiving edge data and thetas.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices.
   * @param count Vertex/edge count.
   * @param thickness Edge half-width (radians).
   * @param min_phi_check Minimum face polar angle (radians).
   * @param max_phi_check Maximum face polar angle (radians).
   * @param h_virt Virtual row count (height plus pole offset).
   * @param height Canvas height in rows.
   * @param y_min_out Output: first covered row.
   * @param y_max_out Output: last covered row.
   */
  static void compute_fast_bounds(FaceScratchBuffer &scratch,
                                  std::span<const Vector> vertices,
                                  std::span<const uint16_t> indices, int count,
                                  float thickness, float min_phi_check,
                                  float max_phi_check, int h_virt, int height,
                                  int &y_min_out, int &y_max_out) {
    for (int i = 0; i < count; ++i) {
      Vector edge = scratch.poly_2d[(i + 1) % count] - scratch.poly_2d[i];
      scratch.edge_vectors[i] = edge;
      float edge_len_sq = dot(edge, edge);
      scratch.edge_lengths_sq[i] = edge_len_sq;
      scratch.inv_edge_lengths_sq[i] =
          (edge_len_sq > 1e-12f) ? (1.0f / edge_len_sq) : 0.0f;
      scratch.inv_edge_j[i] =
          (std::abs(edge.y) > 1e-12f) ? (1.0f / edge.y) : 0.0f;
      float theta = fast_atan2(vertices[indices[i]].z, vertices[indices[i]].x);
      if (theta < 0)
        theta += 2 * PI_F;
      scratch.thetas[i] = theta;
    }
    float margin = thickness + BOUNDS_MARGIN;
    Bounds rows = phi_bounds_to_rows(min_phi_check - margin,
                                     max_phi_check + margin, h_virt, height);
    y_min_out = rows.y_min;
    y_max_out = rows.y_max;
  }

  /**
   * @brief Full-path vertical bounds (arc extrema + pole containment) for large faces.
   * @param scratch Scratch storage receiving edge data, planes, and thetas.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices.
   * @param count Vertex/edge count.
   * @param center Normalized face centroid.
   * @param thickness Edge half-width (radians).
   * @param h_virt Virtual row count (height plus pole offset).
   * @param height Canvas height in rows.
   * @param y_min_out Output: first covered row.
   * @param y_max_out Output: last covered row.
   * @return The number of great-circle planes computed.
   */
  static int compute_full_bounds(FaceScratchBuffer &scratch,
                                 std::span<const Vector> vertices,
                                 std::span<const uint16_t> indices, int count,
                                 const Vector &center, float thickness,
                                 int h_virt, int height, int &y_min_out,
                                 int &y_max_out) {
    float min_phi = 100.0f;
    float max_phi = -100.0f;
    int planes_count = 0;
    for (int i = 0; i < count; ++i) {
      int idx1 = indices[i];
      int idx2 = indices[(i + 1) % count];
      const Vector &v1 = vertices[idx1];
      const Vector &v2 = vertices[idx2];
      Vector edge = scratch.poly_2d[(i + 1) % count] - scratch.poly_2d[i];
      scratch.edge_vectors[i] = edge;
      float edge_len_sq = dot(edge, edge);
      scratch.edge_lengths_sq[i] = edge_len_sq;
      scratch.inv_edge_lengths_sq[i] =
          (edge_len_sq > 1e-12f) ? (1.0f / edge_len_sq) : 0.0f;
      scratch.inv_edge_j[i] =
          (std::abs(edge.y) > 1e-12f) ? (1.0f / edge.y) : 0.0f;
      Vector normal = cross(v1, v2);
      float len_sq = dot(normal, normal);
      if (len_sq > 1e-12f)
        scratch.planes[planes_count++] = normal.normalized();
      float phi_val = fast_acos(hs::clamp(v1.y, -1.0f, 1.0f));
      if (phi_val < min_phi)
        min_phi = phi_val;
      if (phi_val > max_phi)
        max_phi = phi_val;
      // Arc Extrema Logic
      if (planes_count > 0) {
        const Vector &n = scratch.planes[planes_count - 1];
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
            if (cx1 > 0 && cx2 > 0) {
              float phi_top = fast_acos(hs::clamp(pty, -1.0f, 1.0f));
              if (phi_top < min_phi)
                min_phi = phi_top;
            }
            if (cx1 < 0 && cx2 < 0) {
              float phi_bot = fast_acos(hs::clamp(-pty, -1.0f, 1.0f));
              if (phi_bot > max_phi)
                max_phi = phi_bot;
            }
          }
        }
      }
      float theta = fast_atan2(v1.z, v1.x);
      if (theta < 0)
        theta += 2 * PI_F;
      scratch.thetas[i] = theta;
    }
    bool np_inside = (planes_count > 0);
    bool sp_inside = (planes_count > 0);
    for (int pi = 0; pi < planes_count; ++pi) {
      float center_side = dot(center, scratch.planes[pi]);
      if ((scratch.planes[pi].y > 0) != (center_side > 0))
        np_inside = false;
      if ((-scratch.planes[pi].y > 0) != (center_side > 0))
        sp_inside = false;
    }
    if (np_inside)
      min_phi = 0.0f;
    if (sp_inside)
      max_phi = PI_F;
    float margin = thickness + BOUNDS_MARGIN;
    Bounds rows =
        phi_bounds_to_rows(min_phi - margin, max_phi + margin, h_virt, height);
    y_min_out = rows.y_min;
    y_max_out = rows.y_max;
    return planes_count;
  }

  /**
   * @brief Returns the face's precomputed inclusive row bounds.
   * @tparam H Canvas height in rows.
   * @return The stored {y_min, y_max} bounds.
   */
  template <int H> Bounds get_vertical_bounds() const { return {y_min, y_max}; }

  /**
   * @brief Emits the face's azimuth-coverage intervals for a row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param out Sink accepting (float start, float end).
   * @return True if intervals were emitted; false (full scan) for full-width faces.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int, OutputIt out) const {
    if (full_width)
      return false;
    for (const auto &iv : intervals) {
      float f_x1 = iv.first * W / (2 * PI_F);
      float f_x2 = iv.second * W / (2 * PI_F);
      out(floorf(f_x1), ceilf(f_x2));
    }
    return true;
  }

  /**
   * @brief Computes signed distance to the face.
   * @param p Point on sphere (normalized).
   * @return DistanceResult with the signed angular distance and size metric.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  /**
   * @brief Computes signed distance to the face, writing into res.
   * @tparam ComputeUVs Accepted for interface parity; the face stores no UVs.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed angular distance to the edge minus
   *        thickness, raw_dist = unsigned angular distance, size = inradius.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    HS_SCAN_METRIC(hs::g_scan_metrics.pixels_tested++);

    float cos_angle = dot(p, center);
    if (cos_angle <= 0.01f) {
      HS_SCAN_METRIC(hs::g_scan_metrics.pixels_culled++);
      res = DistanceResult(100.0f, 0.0f, 100.0f, 0.0f, size);
      return;
    }

    float inv_cos = 1.0f / cos_angle;
    float px = dot(p, basis_u) * inv_cos;
    float py = dot(p, basis_w) * inv_cos;

    float p_r2 = px * px + py * py;
    if (p_r2 > max_dist_sq) {
      HS_SCAN_METRIC(hs::g_scan_metrics.pixels_culled++);
      res = DistanceResult(100.0f, 0.0f, 100.0f, 0.0f, size);
      return;
    }

    // LUT lookup + same-sign hybrid
    float fx = (px - lut_cx + lut_Rx) * lut_inv_step_x;
    float fy = (py - lut_cy + lut_Ry) * lut_inv_step_y;
    fx = hs::clamp(fx, 0.0f, (float)(LUT_N - 2));
    fy = hs::clamp(fy, 0.0f, (float)(LUT_N - 2));
    int ix = (int)fx;
    int iy = (int)fy;
    float tx = fx - ix;
    float ty = fy - iy;
    float d00 = dist_lut[iy * LUT_N + ix];
    float d10 = dist_lut[iy * LUT_N + ix + 1];
    float d01 = dist_lut[(iy + 1) * LUT_N + ix];
    float d11 = dist_lut[(iy + 1) * LUT_N + ix + 1];

    bool same_sign = (d00 > 0) == (d10 > 0) && (d00 > 0) == (d01 > 0) &&
                     (d00 > 0) == (d11 > 0);
    // Per-face threshold: cell diagonal (tightest correct bound)
    float min_abs =
        std::min({std::abs(d00), std::abs(d10), std::abs(d01), std::abs(d11)});
    float plane_dist;
    if (same_sign && min_abs > lut_safe_dist) {
      HS_SCAN_METRIC(hs::g_scan_metrics.lut_hits++);
      plane_dist = d00 * (1.0f - tx) * (1.0f - ty) + d10 * tx * (1.0f - ty) +
                   d01 * (1.0f - tx) * ty + d11 * tx * ty;
    } else {
      HS_SCAN_METRIC(hs::g_scan_metrics.exact_hits++);
      float d = FLT_MAX;
      bool inside = false;
      for (int i = 0; i < count; ++i) {
        const auto &ep = packed_edges[i];
        float wx = px - ep.vx, wy = py - ep.vy;
        float t = (wx * ep.ex + wy * ep.ey) * ep.inv_len_sq;
        float cv = hs::clamp(t, 0.0f, 1.0f);
        float bx = wx - ep.ex * cv, by = wy - ep.ey * cv;
        float dsq = bx * bx + by * by;
        if (dsq < d)
          d = dsq;
        if ((ep.vy > py) != (ep.next_vy > py)) {
          float isx = ep.vx + (py - ep.vy) * ep.ex * ep.inv_ej;
          if (px < isx)
            inside = !inside;
        }
      }
      plane_dist = (inside ? -1.0f : 1.0f) * sqrtf(d);
    }

    float angular_dist_raw = fast_atan2(plane_dist, 1.0f);
    float angular_dist = angular_dist_raw - thickness;
    res = DistanceResult(angular_dist, 0.0f, angular_dist_raw, 0.0f, size);
  }
};

/**
 * @brief Calculates signed distance to a planar polygon.
 * @details DistanceResult fields: dist = signed distance from edge (negative
 * inside); t = normalized polar distance (polar / thickness); raw_dist = polar
 * distance from center.
 */
struct PlanarPolygon {
  const Basis &basis;                  /**< Orientation frame (v = polygon axis). */
  float thickness;                     /**< Polygon radius / apothem scale (radians). */
  int sides;                           /**< Number of polygon sides. */
  float phase;                         /**< Azimuth phase offset (radians). */
  float apothem;                       /**< Precomputed inradius (radians). */
  float nx, ny, nz, R_val, alpha_angle; /**< Axis components, XZ projection length and azimuth. */
  float phi_min, phi_max;              /**< Vertical bounds as an angular band (radians). */
  static constexpr bool is_solid = true; /**< Polygon renders as a filled region. */

  /**
   * @brief Builds a planar polygon from its basis, thickness, side count, phase.
   * @param b Orientation frame (v = polygon axis).
   * @param th Polygon radius / apothem scale (radians).
   * @param s Number of polygon sides (must be > 0).
   * @param ph Azimuth phase offset (radians).
   */
  PlanarPolygon(const Basis &b, float th, int s, float ph)
      : basis(b), thickness(th), sides(s), phase(ph) {
    HS_CHECK(sides > 0); // sector folding divides by sides
    apothem = thickness * cosf(PI_F / sides);
    nx = basis.v.x;
    ny = basis.v.y;
    nz = basis.v.z;
    R_val = sqrtf(nx * nx + nz * nz);
    alpha_angle = atan2f(nz, nx);

    float center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));
    float margin = thickness + BOUNDS_MARGIN_WIDE;
    phi_min = std::max(0.0f, center_phi - margin);
    phi_max = std::min(PI_F, center_phi + margin);
  }

  /**
   * @brief Maps the polygon's latitude band to its inclusive row range.
   * @tparam H Canvas height in rows.
   * @return Row bounds covering the polygon plus margin.
   */
  template <int H> Bounds get_vertical_bounds() const {
    return phi_bounds_to_rows<H>(phi_min, phi_max);
  }

  /**
   * @brief Emits the bounding-cap interval for one scanline row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    // Pad the cap by one pixel so the outer AA fringe at the polygon tips is
    // scanned, matching Star (process_pixel draws up to d < +pixel_width).
    float pixel_width = 2.0f * PI_F / W;
    return emit_cap_interval<W>(cosf(thickness + pixel_width), ny, R_val,
                                alpha_angle, TrigLUT<W, H>::cos_phi[y],
                                TrigLUT<W, H>::sin_phi[y], false, out);
  }

  /**
   * @brief Signed distance to the planar polygon edge.
   * @param p Point on sphere (normalized).
   * @return DistanceResult with dist, t, and raw_dist populated.
   * @details Treats the polar angle as a flat radius (tangent-plane
   * approximation), so it is accurate for small caps and diverges from the true
   * geodesic distance as the radius grows. Use SphericalPolygon for geodesic
   * great-circle edges.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  /**
   * @brief Signed distance to the planar polygon edge, writing into res.
   * @tparam ComputeUVs When true, also computes the normalized radial t.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = polar*cos(local) - apothem, raw_dist =
   *        polar angle from center, t = polar/thickness when ComputeUVs.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float polar = fast_acos(hs::clamp(dot(p, basis.v), -1.0f, 1.0f));

    float dot_u = dot(p, basis.u);
    float dot_w = dot(p, basis.w);
    float azimuth = fast_atan2(dot_w, dot_u);
    if (azimuth < 0)
      azimuth += 2 * PI_F;
    azimuth += phase;

    float sector = 2 * PI_F / sides;
    float local = wrap(azimuth + sector / 2.0f, sector) - sector / 2.0f;

    float dist_edge = polar * fast_cosf(local) - apothem;
    float t_val = 0.0f;
    if constexpr (ComputeUVs)
      t_val = polar / thickness;

    res = DistanceResult(dist_edge, t_val, polar, 0.0f, apothem);
  }
};

/**
 * @brief Calculates signed distance to a spherical polygon (great-circle edges).
 * @details Uses sector folding plus a precomputed great-circle plane normal for
 * O(1) per-pixel distance, with exact angular distances for smooth AA.
 * DistanceResult fields: dist = signed angular distance to the nearest
 * great-circle edge (negative inside); t = normalized radial position (polar /
 * circumradius); raw_dist = polar angle from the polygon center.
 */
struct SphericalPolygon {
  const Basis &basis;     /**< Orientation frame (v = polygon axis). */
  int sides;              /**< Number of polygon sides. */
  float phase;            /**< Azimuth phase offset (radians). */
  float circumradius;     /**< Angular distance from center to vertex (radians). */
  float edge_nv;          /**< Edge normal dotted with the center axis. */
  float edge_nu;          /**< Edge normal dotted with the u-axis. */
  float phi_min, phi_max; /**< Vertical bounds as an angular band (radians). */
  float nx, ny, nz, R_val, alpha_angle; /**< Axis components, XZ projection length and azimuth. */
  static constexpr bool is_solid = true; /**< Polygon renders as a filled region. */

  /**
   * @brief Builds a spherical polygon from its basis, radius, side count, phase.
   * @param b Orientation frame (v = polygon axis).
   * @param radius Polygon radius as a fraction of the hemisphere.
   * @param s Number of polygon sides (must be > 0).
   * @param ph Azimuth phase offset (radians).
   */
  SphericalPolygon(const Basis &b, float radius, int s, float ph)
      : basis(b), sides(s), phase(ph) {
    HS_CHECK(sides > 0); // sector folding divides by sides
    circumradius = radius * (PI_F / 2.0f);

    // Build canonical edge: between vertices at azimuth ±π/n from
    // the sector bisector (u-axis), at angular distance circumradius
    float half_step = PI_F / sides;
    float sin_r = sinf(circumradius);
    float cos_r = cosf(circumradius);
    float cos_hs = cosf(half_step);
    float sin_hs = sinf(half_step);

    Vector v1 = basis.v * cos_r + (basis.u * cos_hs + basis.w * sin_hs) * sin_r;
    Vector v2 = basis.v * cos_r + (basis.u * cos_hs - basis.w * sin_hs) * sin_r;

    // Normal pointing outward (away from polygon interior)
    Vector en = cross(v2, v1);
    float len = en.magnitude();
    if (len > 1e-9f) {
      en = en * (1.0f / len);
    } else {
      // Degenerate canonical edge: circumradius near 0 or near PI (both
      // reachable — ShapeShifter's radius slider drives circumradius up to PI,
      // where v1==v2==-v so cross()==0), or a near-collinear low-sides polygon.
      // Both degenerate limits drive en toward ±u, so substitute the sector
      // bisector axis as a defined unit normal (edge_nv=0, edge_nu=1) instead of
      // an indeterminate near-zero vector — the polygon then collapses
      // deterministically to its center rather than producing garbage edges.
      en = basis.u;
    }
    // Ensure outward: dot(center, n) should be negative
    if (dot(en, basis.v) > 0)
      en = -en;

    edge_nv = dot(en, basis.v);
    edge_nu = dot(en, basis.u);

    // Vertical/horizontal bounds (same as SDF::PlanarPolygon)
    nx = basis.v.x;
    ny = basis.v.y;
    nz = basis.v.z;
    R_val = sqrtf(nx * nx + nz * nz);
    alpha_angle = atan2f(nz, nx);

    float center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));
    float margin = circumradius + BOUNDS_MARGIN_WIDE;
    phi_min = std::max(0.0f, center_phi - margin);
    phi_max = std::min(PI_F, center_phi + margin);
  }

  /**
   * @brief Maps the polygon's latitude band to its inclusive row range.
   * @tparam H Canvas height in rows.
   * @return Row bounds covering the polygon plus margin.
   */
  template <int H> Bounds get_vertical_bounds() const {
    return phi_bounds_to_rows<H>(phi_min, phi_max);
  }

  /**
   * @brief Emits the bounding-cap interval for one scanline row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    // Pad the cap by one pixel so the outer AA fringe at the polygon tips is
    // scanned, matching Star (process_pixel draws up to d < +pixel_width).
    float pixel_width = 2.0f * PI_F / W;
    return emit_cap_interval<W>(cosf(circumradius + pixel_width), ny, R_val,
                                alpha_angle, TrigLUT<W, H>::cos_phi[y],
                                TrigLUT<W, H>::sin_phi[y], false, out);
  }

  /**
   * @brief Signed distance to the spherical polygon.
   * @param p Point on sphere (normalized).
   * @return DistanceResult with dist, t, and raw_dist populated.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  /**
   * @brief Signed distance to the spherical polygon, writing into res.
   * @tparam ComputeUVs When true, also computes the normalized radial t.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed angular distance to the nearest
   *        great-circle edge, raw_dist = polar angle, t = polar/circumradius
   *        when ComputeUVs.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float polar = angle_between(p, basis.v);

    float dot_u = dot(p, basis.u);
    float dot_w = dot(p, basis.w);
    float azimuth = fast_atan2(dot_w, dot_u);
    if (azimuth < 0)
      azimuth += 2 * PI_F;
    azimuth += phase;

    float sector = 2 * PI_F / sides;
    float local = wrap(azimuth + sector / 2.0f, sector) - sector / 2.0f;

    // Angular distance to the nearest great circle edge via precomputed normal
    // cos(local) is even, so sector folding works automatically
    float sin_p = fast_sinf(polar);
    float cos_p = fast_cosf(polar);
    float dp = edge_nv * cos_p + edge_nu * fast_cosf(local) * sin_p;
    float dist_edge = asinf(hs::clamp(dp, -1.0f, 1.0f));

    float t_val = 0.0f;
    if constexpr (ComputeUVs)
      t_val = polar / circumradius;
    res = DistanceResult(dist_edge, t_val, polar, 0.0f, circumradius);
  }
};

/**
 * @brief Calculates signed distance to a star shape.
 * @details DistanceResult fields: dist = signed distance from edge (negative
 * inside); t = normalized azimuth (azimuth / 2PI); raw_dist = polar distance
 * from center.
 */
struct Star {
  const Basis &basis;                    /**< Orientation frame (v = star axis). */
  int sides;                             /**< Number of star points. */
  float phase;                           /**< Azimuth phase offset (radians). */
  static constexpr bool is_solid = true; /**< Star renders as a filled region. */

  float nx, ny, plane_d; /**< 2D edge plane (normal and offset) for one point. */
  float thickness;       /**< Outer radius / AA scale (radians). */

  float scan_ny, scan_nx, scan_nz, scan_r, scan_alpha; /**< Axis components, XZ projection length and azimuth. */
  float phi_min, phi_max; /**< Vertical bounds as an angular band (radians). */

  /**
   * @brief Builds a star from its basis, radius, point count, and phase.
   * @param b Orientation frame (v = star axis).
   * @param radius Outer radius as a fraction of the hemisphere.
   * @param s Number of star points (must be > 0).
   * @param ph Azimuth phase offset (radians).
   */
  Star(const Basis &b, float radius, int s, float ph)
      : basis(b), sides(s), phase(ph) {
    HS_CHECK(sides > 0); // sector folding divides by sides
    // A zero radius collapses the inner/outer points onto the origin, so the
    // edge vector (dx,dy) below is zero-length and nx/ny = -dy/len, dx/len
    // become NaN, poisoning plane_d and every distance() result. Trap the
    // sizing bug here (cold path, like the sides>0 trap) rather than shipping
    // NaN to the LEDs. Both production consumers keep radius structurally > 0.
    HS_CHECK(radius > 0.0f); // zero radius -> zero-length edge normal (NaN)
    float outer_radius = radius * (PI_F / 2.0f);
    float inner_radius = outer_radius * STAR_INNER_RATIO;
    float angle_step = PI_F / sides;

    float v_t = outer_radius;
    float v_vx = inner_radius * cosf(angle_step);
    float v_vy = inner_radius * sinf(angle_step);

    float dx = v_vx - v_t;
    float dy = v_vy;
    float len = sqrtf(dx * dx + dy * dy);
    nx = -dy / len;
    ny = dx / len;
    plane_d = -(nx * v_t);
    thickness = outer_radius;

    scan_nx = basis.v.x;
    scan_ny = basis.v.y;
    scan_nz = basis.v.z;
    scan_r = sqrtf(scan_nx * scan_nx + scan_nz * scan_nz);
    scan_alpha = atan2f(scan_nz, scan_nx);

    float center_phi = acosf(std::max(-1.0f, std::min(1.0f, basis.v.y)));
    float margin = outer_radius + BOUNDS_MARGIN_WIDE;
    phi_min = std::max(0.0f, center_phi - margin);
    phi_max = std::min(PI_F, center_phi + margin);
  }

  /**
   * @brief Maps the star's latitude band to its inclusive row range.
   * @tparam H Canvas height in rows.
   * @return Row bounds covering the star plus margin.
   */
  template <int H> Bounds get_vertical_bounds() const {
    return phi_bounds_to_rows<H>(phi_min, phi_max);
  }

  /**
   * @brief Emits the bounding-circle interval for one scanline row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan (also when
   *         the cap would span the whole width).
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    // Bounding circle; pad the cap by one pixel and reject full-width rows.
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    float pixel_width = 2.0f * PI_F / W;
    return emit_cap_interval<W>(cosf(thickness + pixel_width), scan_ny, scan_r,
                                scan_alpha, TrigLUT<W, H>::cos_phi[y],
                                TrigLUT<W, H>::sin_phi[y], true, out);
  }

  /**
   * @brief Signed distance to the star.
   * @param p Point on sphere (normalized).
   * @return DistanceResult with dist, t, and raw_dist populated.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  /**
   * @brief Signed distance to the star, writing into res.
   * @tparam ComputeUVs When true, also stores the normalized azimuth t. The
   *        azimuth itself is always computed (the petal-sector geometry needs
   *        it); only the final t store is gated, matching Ring/Flower.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed distance to the nearest point edge,
   *        raw_dist = polar distance from center, t = normalized azimuth when
   *        ComputeUVs (0 otherwise).
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float scan_dist = angle_between(p, basis.v);
    float dot_u = dot(p, basis.u);
    float dot_w = dot(p, basis.w);
    float azimuth = fast_atan2(dot_w, dot_u);
    if (azimuth < 0)
      azimuth += 2 * PI_F;

    azimuth += phase;

    float sector_angle = 2 * PI_F / sides;
    float local_azimuth =
        wrap(azimuth + sector_angle / 2.0f, sector_angle) - sector_angle / 2.0f;
    local_azimuth = std::abs(local_azimuth);

    float px = scan_dist * fast_cosf(local_azimuth);
    float py = scan_dist * fast_sinf(local_azimuth);

    float dist_to_edge = px * nx + py * ny + plane_d;

    float t = 0.0f;
    if constexpr (ComputeUVs)
      t = wrap_t(azimuth / (2 * PI_F));
    res = DistanceResult(-dist_to_edge, t, scan_dist, 0.0f, thickness);
  }
};

/**
 * @brief Calculates signed distance to a flower shape.
 * @details DistanceResult fields: dist = signed distance from the flower edge;
 * t = normalized scan distance (scan_dist / thickness); raw_dist = scan distance
 * from the antipode.
 */
struct Flower {
  const Basis &basis; /**< Orientation frame (v = flower axis). */
  int sides;          /**< Number of petals. */
  float phase;        /**< Azimuth phase offset (radians). */
  float thickness;    /**< Outer radius / AA scale (radians). */
  float apothem;      /**< Petal inradius offset (PI - outer radius). */
  Vector antipode;    /**< Antipode of the flower axis (scan origin). */
  float scan_nx, scan_ny, scan_nz, scan_R, scan_alpha; /**< Antipode components, XZ projection length and azimuth. */
  float phi_min, phi_max; /**< Vertical bounds as an angular band (radians). */
  static constexpr bool is_solid = true; /**< Flower renders as a filled region. */

  /**
   * @brief Builds a flower from its basis, radius, petal count, and phase.
   * @param b Orientation frame (v = flower axis).
   * @param radius Outer radius as a fraction of the hemisphere.
   * @param s Number of petals (must be > 0).
   * @param ph Azimuth phase offset (radians).
   */
  Flower(const Basis &b, float radius, int s, float ph)
      : basis(b), sides(s), phase(ph) {
    HS_CHECK(sides > 0); // sector folding divides by sides
    float outer = radius * (PI_F / 2.0f);
    apothem = PI_F - outer;
    thickness = outer;
    antipode = -basis.v;

    scan_nx = antipode.x;
    scan_ny = antipode.y;
    scan_nz = antipode.z;
    scan_R = sqrtf(scan_nx * scan_nx + scan_nz * scan_nz);
    scan_alpha = atan2f(scan_nz, scan_nx);

    float center_phi = acosf(std::max(-1.0f, std::min(1.0f, antipode.y)));
    float margin = thickness + BOUNDS_MARGIN_WIDE;
    phi_min = std::max(0.0f, center_phi - margin);
    phi_max = std::min(PI_F, center_phi + margin);
  }

  /**
   * @brief Maps the flower's latitude band to its inclusive row range.
   * @tparam H Canvas height in rows.
   * @return Row bounds covering the flower plus margin.
   */
  template <int H> Bounds get_vertical_bounds() const {
    return phi_bounds_to_rows<H>(phi_min, phi_max);
  }

  /**
   * @brief Emits the annular-band interval(s) for one scanline row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    // Phi trig from the static LUT (bit-identical to cosf(y_to_phi(y))) and
    // fast_acos for the band edges — same as Ring's annular fallback path.
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    float cos_phi = TrigLUT<W, H>::cos_phi[y];
    float sin_phi = TrigLUT<W, H>::sin_phi[y];

    if (scan_R < MIN_HORIZONTAL_PROJ)
      return false;

    float denom = scan_R * sin_phi;
    if (std::abs(denom) < INTERVAL_DENOM_EPS)
      return false;

    float ang_max = thickness;
    float cos_limit = cosf(ang_max);

    // Route through the shared annular-band emitter (like Ring/DistortedRing) so
    // the two pole-wraparound degeneracies collapse into a single span instead
    // of emitting holey two-arc intervals near the antipode for large radii.
    emit_annular_band<W>(cos_limit, 1.0f, scan_ny, cos_phi, denom, scan_alpha,
                         out);
    return true;
  }

  /**
   * @brief Signed distance to the flower.
   * @param p Point on sphere (normalized).
   * @return DistanceResult with dist, t, and raw_dist populated.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  /**
   * @brief Signed distance to the flower, writing into res.
   * @tparam ComputeUVs When true, also computes the normalized scan-distance t.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed distance from the petal edge,
   *        raw_dist = scan distance from the antipode, t = scan_dist/thickness
   *        when ComputeUVs.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float scan_dist = angle_between(p, antipode);
    float polar = PI_F - scan_dist;

    float dot_u = dot(p, basis.u);
    float dot_w = dot(p, basis.w);
    float azimuth = fast_atan2(dot_w, dot_u);
    if (azimuth < 0)
      azimuth += 2 * PI_F;
    azimuth += phase;

    float sector = 2 * PI_F / sides;
    float local = wrap(azimuth + sector / 2.0f, sector) - sector / 2.0f;

    float dist_edge = polar * fast_cosf(local) - apothem;
    float t_val = 0.0f;
    if constexpr (ComputeUVs)
      t_val = scan_dist / thickness;

    res = DistanceResult(-dist_edge, t_val, scan_dist, 0.0f, thickness);
  }
};

/**
 * @brief Signed distance to a great-circle arc segment of given thickness.
 * @details DistanceResult fields: dist = signed distance to the segment minus
 * thickness (negative inside); raw_dist = unsigned angular distance to the
 * segment.
 */
struct Line {
  Vector a, b;     /**< Arc endpoints (unit vectors). */
  float thickness; /**< Half-width of the stroke (radians). */

  Vector n;               /**< Great-circle plane normal of the arc. */
  float len;              /**< Arc length (radians). */
  float phi_min, phi_max; /**< Precomputed vertical bounds (radians). */

  // Bounding-cap geometry, loop-invariant across scanlines — hoisted out of
  // get_horizontal_intervals (every other shape precomputes its cap in the
  // constructor and reads the TrigLUT for per-row phi). cap_horiz_valid is
  // false when the cap has no horizontal projection (its axis is ~vertical).
  Vector mid; /**< Arc midpoint axis (bounding-cap center). */
  float mid_ny = 0.0f, mid_r = 0.0f, mid_alpha = 0.0f; /**< Midpoint y, XZ projection length, azimuth. */
  float cap_D_min = 0.0f; /**< Cosine of the bounding-cap radius. */
  bool cap_horiz_valid = false; /**< False when the cap axis is ~vertical (no horizontal projection). */
  static constexpr bool is_solid = false; /**< Line renders as a stroke. */

  /**
   * @brief Builds a great-circle arc segment from two endpoints.
   * @param start First endpoint (unit vector).
   * @param end Second endpoint (unit vector).
   * @param th Half-width of the stroke (radians).
   */
  Line(const Vector &start, const Vector &end, float th)
      : a(start), b(end), thickness(th) {
    len = angle_between(a, b);
    if (len < 1e-6f) {
      n = Vector(0, 0, 0);
    } else {
      // Antipodal endpoints (len ~ π) make cross() degenerate to ~0 with no
      // defined great circle; guard the normalize rather than trap.
      n = normalized_or(cross(a, b), Vector(1, 0, 0));
    }

    // Compute vertical bounds from endpoint phi values + thickness. The
    // endpoints alone are not enough: a great-circle arc bulges toward a pole
    // between them, so its peak latitude can lie well outside [phi_a, phi_b]
    // (e.g. two equal-latitude points whose arc passes through a pole). Extend
    // by the arc's interior latitude extremum when present — the same
    // tangent-sign-flip test plot.h's edge_row_span uses.
    float phi_a = acosf(std::max(-1.0f, std::min(1.0f, a.y)));
    float phi_b = acosf(std::max(-1.0f, std::min(1.0f, b.y)));
    float phi_lo = std::min(phi_a, phi_b);
    float phi_hi = std::max(phi_a, phi_b);
    // y(t) = a.y·cos + v_perp.y·sin turns inside the arc iff the forward
    // tangent's y-component (cross(n, p).y) flips sign between the endpoints.
    // The extremal y is the great circle's peak latitude ±sqrt(1 - n.y²).
    float t0 = cross(n, a).y;
    float t1 = cross(n, b).y;
    if ((t0 > 0.0f) != (t1 > 0.0f)) {
      float peak = sqrtf(std::max(0.0f, 1.0f - n.y * n.y));
      float phi_ext = acosf(t0 > 0.0f ? peak : -peak);
      phi_lo = std::min(phi_lo, phi_ext);
      phi_hi = std::max(phi_hi, phi_ext);
    }
    float margin = thickness + BOUNDS_MARGIN;
    phi_min = std::max(0.0f, phi_lo - margin);
    phi_max = std::min(PI_F, phi_hi + margin);

    // Bounding cap centered on the segment midpoint — all loop-invariant, so
    // precompute here instead of per scanline. Antipodal endpoints sum to ~0
    // (no defined midpoint); guard the normalize.
    mid = normalized_or(a + b, Vector(1, 0, 0));
    mid_ny = mid.y;
    mid_r = sqrtf(mid.x * mid.x + mid.z * mid.z);
    mid_alpha = atan2f(mid.z, mid.x);
    float cap_radius = len * 0.5f + thickness;
    cap_D_min = cosf(cap_radius);
    cap_horiz_valid = mid_r >= MIN_HORIZONTAL_PROJ;
  }

  /**
   * @brief Maps the arc's latitude band to its inclusive row range.
   * @tparam H Canvas height in rows.
   * @return Row bounds covering the arc plus thickness margin.
   */
  template <int H> Bounds get_vertical_bounds() const {
    return phi_bounds_to_rows<H>(phi_min, phi_max);
  }

  /**
   * @brief Signed distance to the arc segment.
   * @param p Point on sphere (normalized).
   * @return DistanceResult with dist and raw_dist populated.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  /**
   * @brief Signed distance to the arc segment, writing into res.
   * @tparam ComputeUVs Accepted for interface parity; the arc stores no UVs.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = unsigned segment distance minus thickness,
   *        raw_dist = unsigned angular distance to the segment.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    if (len < 1e-6f) {
      float dist = angle_between(p, a);
      res = DistanceResult(dist - thickness, 0.0f, dist, 0.0f, thickness);
      return;
    }

    float d_plane = dot(p, n);
    Vector p_proj_plane = p - n * d_plane;
    float proj_mag = p_proj_plane.magnitude();

    if (proj_mag < 1e-6f) {
      float d_a = angle_between(p, a);
      float d_b = angle_between(p, b);
      float dist = std::min(d_a, d_b);
      res = DistanceResult(dist - thickness, 0.0f, dist, 0.0f, thickness);
      return;
    }

    Vector p_proj = p_proj_plane / proj_mag;
    float ang_a = angle_between(a, p_proj);
    float ang_b = angle_between(b, p_proj);

    float dist_seg = 0.0f;
    if (ang_a + ang_b <= len + 1e-4f) {
      dist_seg = asinf(std::abs(d_plane));
    } else {
      float d_a = angle_between(p, a);
      float d_b = angle_between(p, b);
      dist_seg = std::min(d_a, d_b);
    }

    res = DistanceResult(dist_seg - thickness, 0.0f, dist_seg, 0.0f, thickness);
  }

  /**
   * @brief Emits the bounding-cap interval for one scanline row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan (also when
   *         the cap has no horizontal projection).
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    if (!cap_horiz_valid)
      return false;

    // Per-row phi trig from the LUT (indexed by virtual row, so H_OFFSET is
    // baked in), matching every other shape; the cap geometry below was
    // precomputed in the constructor.
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    float cos_phi = TrigLUT<W, H>::cos_phi[y];
    float sin_phi = TrigLUT<W, H>::sin_phi[y];

    float denom = mid_r * sin_phi;
    if (std::abs(denom) < INTERVAL_DENOM_EPS)
      return false;

    float C_min = (cap_D_min - mid_ny * cos_phi) / denom;
    if (C_min > 1.0f)
      return true; // Row is outside the bounding cap
    if (C_min < -1.0f)
      return false; // Cap covers full width

    // fast_acos here matches every other scanline path (Ring/DistortedRing/
    // Flower); C_min is already clamped to [-1,1] by the guards above, and the
    // ~1.3e-4 rad approximation error is far under the floor/ceil pad below.
    float d_alpha = fast_acos(C_min);
    float f_x1 = (mid_alpha - d_alpha) * W / (2 * PI_F);
    float f_x2 = (mid_alpha + d_alpha) * W / (2 * PI_F);
    out(floorf(f_x1), ceilf(f_x2));
    return true;
  }
};

// ============================================================================
// 3D Volumetric SDF Shapes (for Scan::Volume raymarching)
// ============================================================================

/**
 * @brief 3D Torus signed distance field.
 *
 * The torus ring lies in the XZ plane with symmetry axis along Y.
 * Major radius R = ring centerline distance from origin.
 * Minor radius r = tube cross-section radius.
 *
 * Unlike the 2D spherical SDF shapes above (which use DistanceResult),
 * 3D volumetric shapes return plain float distances and operate in
 * Cartesian ray-space.
 */
struct Torus {
  float R; /**< Major radius (ring centerline). */
  float r; /**< Minor radius (tube cross-section). */

  /**
   * @brief Signed distance from a point to the torus surface.
   * @param p Query point in Cartesian ray-space.
   * @return Signed distance (negative inside, positive outside).
   */
  float distance(const Vector &p) const {
    float q = sqrtf(p.x * p.x + p.z * p.z) - R;
    return sqrtf(q * q + p.y * p.y) - r;
  }

  /**
   * @brief Surface normal at a point near the torus surface.
   * @param p Query point in Cartesian ray-space.
   * @return Unit outward normal.
   */
  Vector normal(const Vector &p) const {
    float xz_len = sqrtf(p.x * p.x + p.z * p.z);
    float scale = (xz_len > TOLERANCE) ? R / xz_len : 0.0f;
    return (p - Vector(p.x * scale, 0.0f, p.z * scale)).normalized();
  }

  /**
   * @brief Populates a Fragment's registers for shading.
   * @param p Query point in Cartesian ray-space.
   * @param frag Output fragment; v0 = ring angle (0-1, for palette lookup),
   *        v1/v2/v3 = surface normal (x, y, z).
   */
  void populate(const Vector &p, Fragment &frag) const {
    Vector n = normal(p);
    frag.v0 = (fast_atan2(p.z, p.x) + PI_F) / (2.0f * PI_F);
    frag.v1 = n.x;
    frag.v2 = n.y;
    frag.v3 = n.z;
  }
};

/**
 * @brief Domain warp functions for composing with WarpedSDF.
 */
namespace Warp {

/**
 * @brief Oscillates the Y coordinate sinusoidally around the azimuthal
 * angle θ = atan2(z, x).
 *
 * Produces twisted/undulating geometry when composed with a torus.
 * Provides an analytical Lipschitz bound for safe sphere tracing.
 */
struct Twist {
  int twist;       /**< Number of oscillations around the ring. */
  float amplitude; /**< Vertical displacement magnitude. */
  float R;         /**< Major radius (needed for the Lipschitz bound). */

  /** @brief Precomputed context: s = sqrtf(x² + z²), shared across apply/lipschitz. */
  using Ctx = float;

  /**
   * @brief Precomputes the shared per-point context s = sqrtf(x² + z²).
   * @param p Query point.
   * @return The radial distance s in the XZ plane.
   */
  Ctx make_ctx(const Vector &p) const { return sqrtf(p.x * p.x + p.z * p.z); }

  /**
   * @brief Warps the domain by displacing Y by amplitude * sin(twist * θ).
   * @param p Query point.
   * @return The warped point.
   */
  Vector apply(const Vector &p, Ctx /*s*/) const {
    float theta = fast_atan2(p.z, p.x);
    return Vector(
        p.x, p.y - amplitude * fast_sinf(static_cast<float>(twist) * theta), p.z);
  }

  /**
   * @brief Analytical Lipschitz constant of the warp at a point.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @return The exact operator norm of the warp Jacobian (>= 1).
   */
  float lipschitz(const Vector & /*p*/, Ctx s) const {
    if (twist == 0)
      return 1.0f;
    // The warp Jacobian is the shear I - e_y·gᵀ with e_y ⊥ g and |g| = γ; its
    // exact operator norm (largest singular value) is γ/2 + √(1 + γ²/4). This
    // exact bound prevents Raymarch's sphere tracing from over-estimating the
    // safe step and marching through thin surfaces at high Twist. γ uses
    // |twist·amplitude| so the bound stays conservative regardless of sign.
    const float gamma =
        fabsf(static_cast<float>(twist) * amplitude) / std::max(s, R * 0.5f);
    return 0.5f * gamma + sqrtf(1.0f + 0.25f * gamma * gamma);
  }

  /**
   * @brief Maximum possible inflation of the bounding volume.
   * @return The displacement amplitude (radians of XYZ space).
   */
  float bounding_inflation() const { return amplitude; }

  /**
   * @brief Analytical normal correction via the chain rule (once per hit).
   * @param p Query point.
   * @param base_n Unwarped surface normal.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @return The corrected unit normal accounting for the warp.
   */
  Vector correct_normal(const Vector &p, const Vector &base_n, Ctx s) const {
    if (twist == 0 || amplitude < TOLERANCE)
      return base_n;
    float inv_s = (s > TOLERANCE) ? 1.0f / s : 0.0f;
    float theta = fast_atan2(p.z, p.x);
    float n_theta = static_cast<float>(twist) * theta;
    float dh_dtheta = -amplitude * static_cast<float>(twist) * fast_cosf(n_theta);
    float inv_s2 = inv_s * inv_s;

    float dh_dx = dh_dtheta * (-p.z) * inv_s2;
    float dh_dz = dh_dtheta * p.x * inv_s2;

    return Vector(base_n.x - base_n.y * dh_dx, base_n.y,
                  base_n.z - base_n.y * dh_dz)
        .normalized();
  }
};

} // namespace Warp

/**
 * @brief Composable domain-warped volume SDF.
 *
 * Applies a Warp to the input domain of any 3D volume SDF, with:
 * - Bounding fast-path (skips warp trig when far from surface)
 * - Analytical Lipschitz correction for safe sphere tracing
 * - Normal correction via the warp's chain rule
 *
 * The Warp must provide:
 *   Ctx make_ctx(const Vector &p) const;
 *   Vector apply(const Vector &p, const Ctx &ctx) const;
 *   float lipschitz(const Vector &p, const Ctx &ctx) const;
 *   float bounding_inflation() const;
 * Optionally:
 *   Vector correct_normal(const Vector &p, const Vector &n, const Ctx&) const;
 */
template <typename SDF, typename Warp> struct WarpedVolume {
  SDF base;  /**< The underlying volume SDF being warped. */
  Warp warp; /**< The domain warp applied before the base SDF. */

  /**
   * @brief Cheap lower bound: base distance minus the warp's max displacement.
   * @param p Query point in Cartesian ray-space.
   * @return A conservative lower bound on the warped distance.
   */
  float bounding_distance(const Vector &p) const {
    return base.distance(p) - warp.bounding_inflation();
  }

  /**
   * @brief Raw warped SDF distance with no Lipschitz correction.
   * @param p Query point in Cartesian ray-space.
   * @return The base SDF evaluated at the warped point (use for surface projection).
   */
  float raw_distance(const Vector &p) const {
    auto ctx = warp.make_ctx(p);
    return base.distance(warp.apply(p, ctx));
  }

  /**
   * @brief March-safe distance with bounding fast-path and Lipschitz correction.
   * @param p Query point in Cartesian ray-space.
   * @return A sphere-tracing-safe (under-estimated) distance to the surface.
   */
  float distance(const Vector &p) const {
    float bd = bounding_distance(p);
    if (bd > warp.bounding_inflation())
      return bd;

    auto ctx = warp.make_ctx(p);
    float d = base.distance(warp.apply(p, ctx));

    if (d > 0.0f) {
      float lip = warp.lipschitz(p, ctx);
      if (lip > 1.0f)
        d /= lip;
    }
    return d;
  }

  /**
   * @brief Surface normal with the warp's chain-rule correction.
   * @param p Query point in Cartesian ray-space.
   * @return The corrected unit surface normal.
   */
  Vector normal(const Vector &p) const {
    auto ctx = warp.make_ctx(p);
    Vector base_n = base.normal(warp.apply(p, ctx));
    if constexpr (requires { warp.correct_normal(p, base_n, ctx); }) {
      return warp.correct_normal(p, base_n, ctx);
    }
    return base_n;
  }

  /**
   * @brief Populates a Fragment's registers for shading.
   * @param p Query point in Cartesian ray-space.
   * @param frag Output fragment; v0 = ring angle (0-1), v1/v2/v3 = surface normal.
   */
  void populate(const Vector &p, Fragment &frag) const {
    Vector n = normal(p);
    frag.v0 = (fast_atan2(p.z, p.x) + PI_F) / (2.0f * PI_F);
    frag.v1 = n.x;
    frag.v2 = n.y;
    frag.v3 = n.z;
  }
};

} // namespace SDF
