/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <cfloat>
#include <span>
#include <type_traits>
#include "math/geometry.h"
#include "render/shading.h"
#include "engine/constants.h"
#include "engine/concepts.h"
#include "engine/static_circular_buffer.h"
#include "engine/util.h"

#ifdef HS_AA_AUDIT
#include "tests/aa_audit.h"
#endif

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
/** Inner/outer radius ratio for star shapes. */
static constexpr float STAR_INNER_RATIO = ::STAR_INNER_RATIO;
/** Minimum inradius-to-circumradius ratio used to floor Face::size,
 *  preventing degenerate near-zero inradii from collapsing AA. */
static constexpr float MIN_SIZE_RADIUS_RATIO = 0.25f;

/** Signed-area-to-circumradius-squared ratio below which a Face is culled as
 *  fully collapsed (no enclosed region). Sits orders of magnitude above the
 *  float noise of an exactly collapsed polygon (~1e-7) and below the thinnest
 *  real sliver a mesh sweep draws (~1e-3), so the sim/device decision is
 *  identical under fast-math. */
static constexpr float COLLAPSED_AREA_RATIO = 1e-5f;

/** Maximum disjoint scanline spans a single shape (leaf) emits per row.
 *  scan_region's `intervals` buffer holds the widest top-level CSG emission
 *  (Subtract/Intersection: 2x this + 2), and its seam-split `norm` buffer is
 *  2x intervals (one span can split in two at the x=0 seam); see the
 *  static_asserts there. */
inline constexpr size_t INTERVAL_SPAN_CAP = 32;

/** Per-row scanline interval buffer for a single shape. Fixed capacity,
 *  accumulate-only. */
using IntervalBuffer =
    StaticCircularBuffer<std::pair<float, float>, INTERVAL_SPAN_CAP>;

/** Per-row accumulator for a binary CSG op that merges BOTH children's spans
 *  into one buffer (Union/SmoothUnion). Each child can contribute up to
 *  INTERVAL_SPAN_CAP spans, so the union accumulator is sized to hold both. */
using MergedIntervalBuffer =
    StaticCircularBuffer<std::pair<float, float>, 2 * INTERVAL_SPAN_CAP>;

// Forward-declared so the span-count trait below can pattern-match the binary
// CSG ops (defined later in this header) before their full definitions.
template <typename A, typename B> struct Union;
template <typename A, typename B> struct SmoothUnion;
template <typename A, typename B> struct Subtract;
template <typename A, typename B> struct Intersection;

/** Compile-time upper bound on the scanline spans a shape may emit to its
 * parent in one row. A leaf is capped at INTERVAL_SPAN_CAP; Union/SmoothUnion
 * merge both children into one MergedIntervalBuffer, so their bound is the SUM
 * of the children's, static_asserted against the buffer capacity to reject an
 *  overflowing nesting at compile time. */
template <typename T> struct sdf_max_spans {
  static constexpr size_t value = INTERVAL_SPAN_CAP;
};
template <typename A, typename B> struct sdf_max_spans<Union<A, B>> {
  static constexpr size_t value =
      sdf_max_spans<A>::value + sdf_max_spans<B>::value;
};
template <typename A, typename B> struct sdf_max_spans<SmoothUnion<A, B>> {
  static constexpr size_t value =
      sdf_max_spans<A>::value + sdf_max_spans<B>::value;
};
// Intersection seam-splits each child into a [0, W) frame, then merge-sweeps
// the two start-sorted lists (one span per advance). Each child grows by at
// most one and the sweep advances |norm_a| + |norm_b| times: bound |A| + |B|
// + 2.
template <typename A, typename B> struct sdf_max_spans<Intersection<A, B>> {
  static constexpr size_t value =
      sdf_max_spans<A>::value + sdf_max_spans<B>::value + 2;
};
// Subtract seam-splits each child into a [0, W) frame before differencing; each
// child grows by at most one (norm_a <= |A|+1, norm_b <= |B|+1). The set
// difference splits an A span once per enclosed (disjoint) B span, so the
// output is bounded by |norm_a| + |norm_b| = |A| + |B| + 2.
template <typename A, typename B> struct sdf_max_spans<Subtract<A, B>> {
  static constexpr size_t value =
      sdf_max_spans<A>::value + sdf_max_spans<B>::value + 2;
};

/**
 * @brief Append a scanline interval, trapping on overflow.
 * @tparam N Buffer capacity (deduced); supports both the per-shape and the
 * two-child union accumulators.
 * @param buf Per-row interval buffer to append to.
 * @param start Interval start column (float).
 * @param end Interval end column (float).
 * @details StaticCircularBuffer::push_back evicts the OLDEST entry when full
 * (correct for trails, wrong here), so an overflow would silently drop
 * geometry. A row exceeding capacity is a sizing bug, so trap at the violation
 * site.
 */
template <size_t N>
inline void push_interval(StaticCircularBuffer<std::pair<float, float>, N> &buf,
                          float start, float end) {
  HS_CHECK(!buf.is_full(), "SDF scanline interval buffer overflow in one row");
  buf.push_back({start, end});
}

/**
 * @brief Insertion-sort an interval buffer in place by start coordinate.
 * @param buf Per-row interval buffer to sort in place.
 * @details Raw-pointer indexing (buffer freshly built, head == 0, contiguous)
 * avoids the per-access modulo. Shared by merge_intervals and Subtract.
 */
template <size_t N>
inline void
sort_intervals_by_start(StaticCircularBuffer<std::pair<float, float>, N> &buf) {
  HS_CHECK(buf.is_linear(),
           "sort_intervals_by_start: raw linear indexing requires head==0");
  auto *data = &buf[0];
  size_t n = buf.size();
  for (size_t i = 1; i < n; ++i) {
    auto key = data[i];
    size_t j = i;
    while (j > 0 && data[j - 1].first > key.first) {
      data[j] = data[j - 1];
      --j;
    }
    data[j] = key;
  }
}

/**
 * @brief Wrap each interval start into [0, W) and split any span that crosses
 * the x=0 seam, appending the result to @p dst.
 * @tparam W Canvas width in columns.
 * @param src Source intervals in unwrapped column space (may straddle θ=0).
 * @param dst Destination buffer; must hold up to 2x the source span count (one
 *        span splits into at most two at the seam).
 * @details Mirrors scan_region's normalization so seam math is bit-identical. A
 * span of length >= W is emitted as a single full-row [0, W) span. For the
 * common in-[0,W) case this copies through unchanged.
 */
template <int W, size_t N, size_t M>
inline void normalize_intervals_to_range(
    const StaticCircularBuffer<std::pair<float, float>, N> &src,
    StaticCircularBuffer<std::pair<float, float>, M> &dst) {
  constexpr float Wf = static_cast<float>(W);
  for (size_t i = 0; i < src.size(); ++i) {
    float len = src[i].second - src[i].first;
    HS_CHECK(len >= 0.0f,
             "normalize_intervals_to_range: interval end precedes start");
    if (len >= Wf) {
      push_interval(dst, 0.0f, Wf);
      continue;
    }
    float s = fmodf(src[i].first, Wf);
    if (s < 0.0f)
      s += Wf;
    float e = s + len;
    if (e <= Wf) {
      push_interval(dst, s, e);
    } else {
      push_interval(dst, s, Wf);
      push_interval(dst, 0.0f, e - Wf);
    }
  }
}

/**
 * @brief Sort an interval buffer by start, then emit the union of overlapping
 * intervals via out(start, end). Shared by Union/SmoothUnion.
 * @tparam N Buffer capacity (deduced).
 * @tparam OutputIt Sink type invoked as out(float start, float end).
 * @param merged Per-row interval buffer (sorted and merged in place).
 * @param out Sink receiving each merged interval.
 * @details Precondition: `merged` is non-empty (callers guard with is_empty()).
 * Templated on the output sink so it inlines at -O3.
 */
template <size_t N, typename OutputIt>
inline void
merge_intervals(StaticCircularBuffer<std::pair<float, float>, N> &merged,
                OutputIt out) {
  sort_intervals_by_start(merged);
  HS_CHECK(merged.is_linear(),
           "merge_intervals: raw linear indexing requires head==0");
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
 * @brief Fold any real angle into [0, π], equivalent to acosf(cosf(x)) without
 *        trig.
 * @param x Angle in radians (any finite value).
 * @return Folded angle in [0, π].
 * @details cos is even and 2π-periodic, so fold the sign, reduce mod 2π, then
 *          reflect the upper half-period across the south pole. The full fold
 *          (not the [-π, 2π]-only short form) holds the equivalence for any
 * input (e.g. a Ring radius > 2 driving center_phi ± target_angle past range).
 */
inline float clamp_phi(float x) {
  x = fabsf(x);              // cos(-x) = cos(x): fold negatives
  x = fmodf(x, 2.0f * PI_F); // 2π-periodic -> [0, 2π)
  if (x > PI_F)
    x = 2.0f * PI_F - x; // reflect (π, 2π) across the south pole -> (0, π)
  return x;
}

/** @brief A colatitude band as inclusive [phi_min, phi_max] bounds in [0, π].
 */
struct PhiBand {
  float phi_min, phi_max;
};

/**
 * @brief Folds a `center ± half-angle` colatitude band into clamped [0, π].
 * @param center_phi Band-center colatitude (radians).
 * @param target_angle Half-width of the band (radians).
 * @return {phi_min, phi_max}: the band's folded extent. An edge that runs past
 * a pole pins that side to the pole (0 or π); a fully in-range band returns its
 * folded endpoints.
 * @details Single source for the Ring/DistortedRing get_vertical_bounds
 * latitude fold so the two cannot drift apart.
 */
inline PhiBand clamp_phi_band(float center_phi, float target_angle) {
  float a1 = center_phi - target_angle;
  float a2 = center_phi + target_angle;
  float p1 = clamp_phi(a1);
  float p2 = clamp_phi(a2);
  float phi_min = 0, phi_max = PI_F;

  if (a1 > 0)
    phi_min = std::min(p1, p2);
  if (a2 < PI_F)
    phi_max = std::max(p1, p2);
  return {phi_min, phi_max};
}
/**
 * @brief Vertical scanline bounds (inclusive min/max row index).
 */
struct Bounds {
  int y_min, y_max; /**< Inclusive first/last row covered. */
};
/**
 * @brief Result of a signed distance query.
 *
 * `dist` and `size` have fixed meanings, but `t`, `raw_dist` and `aux` are
 * overloaded per shape — the authoritative meaning is that shape's own
 * "Returns:" docblock. The scan rasterizer copies them into the Fragment
 * register file with no reinterpretation (see Scan::process_pixel): t        ->
 * Fragment::v0 raw_dist -> Fragment::v1 aux      -> Fragment::v3 size     ->
 * Fragment::size Fragment::v2 is generated downstream: Scan writes stroke AA
 * coverage or 0 for solid shapes, and Scan::Mesh replaces it with a face index.
 */
struct DistanceResult {
  float dist; /**< Signed distance (negative inside); always this meaning. */
  float t;    /**< Per-shape: normalized parameter (0-1) or angle. */
  float raw_dist; /**< Per-shape: unsigned or supplementary distance. */
  float aux; /**< Per-shape: auxiliary value (e.g. barycentric coordinate). */
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
 * @brief Structural fingerprint shared by every SDF leaf and CSG combinator:
 * a static is_solid flag and an AA-falloff thickness.
 * @tparam T Candidate shape type.
 * @details The CSG combinators assert this on their children so a wrong-type
 * argument fails at the boundary. distance() and the scanline members vary by
 * render path and are not part of the shared contract.
 */
template <typename T>
concept SDFShape = requires(const T &t) {
  { T::is_solid } -> std::convertible_to<bool>;
  { t.thickness } -> std::convertible_to<float>;
};

/**
 * @brief Axis components plus its scan-plane projection: XZ-projection length
 * R_val and azimuth alpha_angle.
 */
struct AxisProjection {
  float nx, ny, nz, R_val, alpha_angle;
};

/**
 * @brief Decomposes an axis into components and its scan-plane projection.
 * @param axis The axis vector.
 * @return Components nx/ny/nz, XZ-projection length R_val, azimuth alpha_angle.
 */
inline AxisProjection project_axis(const Vector &axis) {
  return {axis.x, axis.y, axis.z, sqrtf(axis.x * axis.x + axis.z * axis.z),
          atan2f(axis.z, axis.x)};
}

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
    return false; // Full scan fallback

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
 * and the upper edge ceils so a partially-covered row is never dropped, and the
 * row clamps fold out-of-range phi. Single source for the floor/ceil
 * conversion, routed through by `phi_bounds_to_rows<H>` and the Face's
 * construction bounds.
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
template <int H>
inline Bounds phi_bounds_to_rows(float phi_min, float phi_max) {
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
 * @return False if the band misses this row; true with angles written
 * otherwise.
 * @details cos decreases with angle, so the larger cosine (cos_inner) yields
 * the smaller angle. Shared by the annular scanline emitters.
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
 * t = normalized parameter (0-1) corresponding to angle/2PI; raw_dist =
 * unsigned distance to centerline.
 */
struct Ring {
  const Basis &basis; /**< Orientation frame (v = ring axis). */
  float radius;       /**< Ring radius as a fraction of the hemisphere. */
  float thickness;    /**< Half-width of the stroke (radians). */
  float phase;        /**< Azimuth phase offset (radians). */

  Vector normal, u, w; /**< Ring axis and the two in-plane basis vectors. */
  float ny;            /**< y-component of the ring axis. */
  float target_angle,
      center_phi; /**< Centerline polar angle and axis colatitude. */
  float cos_max, cos_min, cos_target, inv_sin_target,
      sin_target; /**< Precomputed band trig. */

  float r_val;       /**< Horizontal projection length of the axis (for full-row
                        check). */
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
    AxisProjection ap = project_axis(normal);
    ny = ap.ny;

    target_angle = radius * (PI_F / 2.0f);
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

    r_val = ap.R_val;
    alpha_angle = ap.alpha_angle;
  }

  /**
   * @brief Maps the ring's latitude band to its inclusive scanline row range.
   * @tparam H Canvas height in rows.
   * @return Row bounds covering the ring plus its AA falloff.
   */
  template <int H> Bounds get_vertical_bounds() const {
    PhiBand band = clamp_phi_band(center_phi, target_angle);

    float eff_th = 0.95f * thickness; // quintic_kernel(0.05) ≈ 0.001
    float f_phi_min = std::max(0.0f, band.phi_min - eff_th);
    float f_phi_max = std::min(PI_F, band.phi_max + eff_th);

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
   * @return True if intervals were found and reported; false requests a full
   * scan.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    float cos_phi = TrigLUT<W, H>::cos_phi[y];
    float sin_phi = TrigLUT<W, H>::sin_phi[y];

    if (needs_full_row_scan(sin_phi))
      return false;

    float denom = r_val * sin_phi;
    float ny_cos_phi = ny * cos_phi;
    float C_target = (cos_target - ny_cos_phi) / denom;
    float scale = W / (2.0f * PI_F);

    // Pole-wrap: when a band edge arc collapses within one column of a pole,
    // emit_annular_band merges its two arcs but the centerline fast path below
    // emits them unmerged, leaving a seam gap. cos decreases with angle, so the
    // near edge cos_max gives the smaller angle.
    float C_band_near = (cos_max - ny_cos_phi) / denom;
    float C_band_far = (cos_min - ny_cos_phi) / denom;
    float cos_pole = cosf(2.0f * PI_F / W);
    bool pole_wrap = C_band_near >= cos_pole || C_band_far <= -cos_pole;

    // Centerline fast path, valid only while the band edges stay inside [-1,1];
    // past the clamp boundary the first-order width formula spikes, so fall
    // back to the exact annular arcs.
    float clamp_bound = 1.0f - sin_target * thickness / denom;
    if (!pole_wrap && clamp_bound < 1.0f && C_target > -clamp_bound &&
        C_target < clamp_bound && sin_target > 1e-3f) {
      // Floor the radicand: 1-C^2 nears 0 and can go slightly negative under
      // -ffast-math (sqrtf(NaN)); the floor also caps half_width on grazing
      // rows.
      float sin_cross = sqrtf(std::max(1.0f - C_target * C_target, 1e-6f));
      float acos_C = fast_acos(C_target);
      // Full thickness: the scan interval must span the stroke-AA footprint
      // process_pixel derives from full `size`, else the outermost AA column
      // clips on grazing rows.
      float half_width = thickness * sin_target / (denom * sin_cross);

      float hw_px = half_width * scale;
      float t1 = (alpha_angle - acos_C) * scale;
      float t2 = (alpha_angle + acos_C) * scale;

      // Both padded twin arcs can straddle θ=0; the seam merge is left to
      // scan_region's [0,W) wrap/coalesce pass.
      out(floorf(t1 - hw_px), ceilf(t1 + hw_px));
      out(floorf(t2 - hw_px), ceilf(t2 + hw_px));
      return true;
    }

    // Annular band: exact intervals for near-tangent / out-of-range rows
    emit_annular_band<W>(cos_min, cos_max, ny, cos_phi, denom, alpha_angle,
                         out);
    return true;
  }

  /**
   * @brief Whether row interval math is degenerate and the row must be
   *        full-row scanned.
   * @param sin_phi sin of the row's colatitude.
   * @return True when get_horizontal_intervals would return false at this row.
   */
  bool needs_full_row_scan(float sin_phi) const {
    return r_val < MIN_HORIZONTAL_PROJ ||
           std::abs(r_val * sin_phi) < INTERVAL_DENOM_EPS;
  }

  /**
   * @brief Stroke coverage from a precomputed axis dot, skipping the
   *        DistanceResult round trip.
   * @param d dot(p, normal) for the pixel's unit vector.
   * @return quintic stroke alpha in [0, 1]; 0 outside the band. Same distance
   *         branches and float ops as distance<false>() + process_pixel's
   *         stroke epilogue, so coverage is bit-identical to that path.
   */
  __attribute__((always_inline)) float stroke_alpha(float d) const {
    if (d < cos_min || d > cos_max)
      return 0.0f;
    float dist;
    if (inv_sin_target != 0) {
      dist = std::abs(d - cos_target) * inv_sin_target;
    } else {
      float polar = fast_acos(hs::clamp(d, -1.0f, 1.0f));
      dist = std::abs(polar - target_angle);
    }
    float sd = dist - thickness;
    if (sd >= 0.0f || thickness <= 0.0f)
      return 0.0f;
    return quintic_kernel(-sd / thickness);
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
 * t = normalized parameter (0-1) corresponding to angle/2PI; raw_dist =
 * unsigned distance to centerline.
 */
struct DistortedRing {
  const Basis &basis; /**< Orientation frame (v = ring axis). */
  float radius;       /**< Ring radius as a fraction of the hemisphere. */
  float thickness;    /**< Half-width of the stroke (radians). */
  ScalarFn shift_fn;  /**< Per-azimuth centerline shift, t in [0,1) -> radians;
                         empty in knot mode. */
  const float *knots =
      nullptr;   /**< Optional lut_n + 1 shift knots (entry lut_n repeats entry
                    0); selects exact polyline distance. */
  int lut_n = 0; /**< Knot cell count when knots is set. */
  float max_distortion; /**< Maximum magnitude of the shift (radians). */
  float phase;          /**< Azimuth phase offset (radians). */

  Vector normal, u, w; /**< Ring axis and the two in-plane basis vectors. */
  float ny;            /**< y-component of the ring axis. */
  float target_angle,
      center_phi;      /**< Centerline polar angle and axis colatitude. */
  float max_thickness; /**< thickness + max_distortion (radians). */

  float r_val;       /**< Horizontal projection length of the axis. */
  float alpha_angle; /**< Azimuth of the normal in the XZ plane. */
  float cos_max_limit, cos_min_limit; /**< Cosines of the widened band edges. */
  bool suppress_pole_fill =
      false; /**< Drop the degenerate exact-pole row rather than full-row
                filling it (see get_horizontal_intervals). */
  static constexpr bool is_solid =
      false; /**< Distorted ring renders as a stroke. */

  /**
   * @brief Builds a distorted ring with a per-azimuth centerline shift.
   * @param b Orientation frame (v = ring axis).
   * @param r Ring radius as a fraction of the hemisphere.
   * @param th Half-width of the stroke (radians).
   * @param sf Per-azimuth centerline shift function, t in [0,1) -> radians.
   * @param md Maximum magnitude of sf over t in [0,1) (radians). PRECONDITION:
   *           md must be a true upper bound on |sf|. It widens the reject
   * bands, so an underestimate silently culls genuine arcs. Pinned by the cull
   *           tests, not checked here.
   * @param ph Azimuth phase offset (radians).
   */
  DistortedRing(const Basis &b, float r, float th, ScalarFn sf, float md,
                float ph)
      : basis(b), radius(r), thickness(th), shift_fn(sf), max_distortion(md),
        phase(ph) {
    normal = basis.v;
    u = basis.u;
    w = basis.w;
    AxisProjection ap = project_axis(normal);
    ny = ap.ny;
    target_angle = radius * (PI_F / 2.0f);
    center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));
    max_thickness = thickness + max_distortion;

    r_val = ap.R_val;
    alpha_angle = ap.alpha_angle;

    float ang_min = std::max(0.0f, target_angle - max_thickness);
    float ang_max = std::min(PI_F, target_angle + max_thickness);
    cos_max_limit = cosf(ang_min);
    cos_min_limit = cosf(ang_max);
  }

  /**
   * @brief Builds a distorted ring whose centerline is a shift-knot polyline.
   * @param b Orientation frame (v = ring axis).
   * @param r Ring radius as a fraction of the hemisphere.
   * @param th Half-width of the stroke (radians).
   * @param kn n + 1 centerline shifts (radians), one per equal azimuth cell,
   *           entry n repeating entry 0; must outlive the shape. distance()
   *           returns the exact distance to this polyline (within the local
   *           tangent chart), so steep or sharply curved segments render at
   *           full stroke width with no slope approximation.
   * @param n Number of knot cells; at least 1.
   * @param md Maximum magnitude of the knots (radians); same true-upper-bound
   *           precondition as the ScalarFn overload.
   * @param ph Azimuth phase offset (radians).
   */
  DistortedRing(const Basis &b, float r, float th, const float *kn, int n,
                float md, float ph)
      : DistortedRing(b, r, th, ScalarFn{}, md, ph) {
    knots = kn;
    lut_n = n;
    float min_shift = kn[0];
    float max_shift = kn[0];
    // Per-chunk knot ranges for the per-pixel prefilter. A segment registers
    // its endpoints in every chunk its azimuth extent touches (a straddling
    // segment spans two), so the polyline inside a chunk never leaves
    // [chunk_lo, chunk_hi].
    for (int c = 0; c < PREFILTER_CHUNKS; ++c) {
      chunk_lo[c] = 1e9f;
      chunk_hi[c] = -1e9f;
    }
    for (int k = 0; k < n; ++k) {
      float lo = std::min(kn[k], kn[k + 1]);
      float hi = std::max(kn[k], kn[k + 1]);
      min_shift = std::min(min_shift, lo);
      max_shift = std::max(max_shift, hi);
      int c1 = k * PREFILTER_CHUNKS / n;
      int c2 = std::min((k + 1) * PREFILTER_CHUNKS / n, PREFILTER_CHUNKS - 1);
      for (int c = c1; c <= c2; ++c) {
        chunk_lo[c] = std::min(chunk_lo[c], lo);
        chunk_hi[c] = std::max(chunk_hi[c], hi);
      }
    }

    max_distortion = std::max(std::abs(min_shift), std::abs(max_shift));
    max_thickness = thickness + max_distortion;
    float ang_min = hs::clamp(target_angle + min_shift - thickness, 0.0f, PI_F);
    float ang_max = hs::clamp(target_angle + max_shift + thickness, 0.0f, PI_F);
    cos_max_limit = cosf(ang_min);
    cos_min_limit = cosf(ang_max);
  }

  /**
   * @brief Maps the distorted ring's widened latitude band to its row range.
   * @tparam H Canvas height in rows.
   * @return Inclusive row bounds covering the ring plus distortion margin.
   */
  template <int H> Bounds get_vertical_bounds() const {
    PhiBand band = clamp_phi_band(center_phi, target_angle);

    float margin = max_thickness + BOUNDS_MARGIN_WIDE;
    float f_phi_min = std::max(0.0f, band.phi_min - margin);
    float f_phi_max = std::min(PI_F, band.phi_max + margin);

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
      // r_val cleared MIN_HORIZONTAL_PROJ above, so a vanishing denom is the
      // exact pole row: every column aliases to the one pole point. A displaced
      // stroke reaches the pole at a single azimuth, but the degenerate row
      // math can't recover which column, so the default (return false) full-row
      // scans and fills the whole aliased row -- fine for a lone ring, but a
      // dense ring stack renders it as a solid pole cap. suppress_pole_fill
      // drops the row.
      return suppress_pole_fill;

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
   * @param res Output result; dist = signed distance minus thickness, raw_dist
   * = unsigned centerline distance, t = azimuth in [0,1) when ComputeUVs.
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

    float dist;
    if (knots)
      dist = polyline_distance(t_norm, polar,
                               sqrtf(std::max(1.0f - d * d, POLE_SIN2_FLOOR)));
    else
      dist = std::abs(polar - (target_angle + shift_fn(t_norm)));

    if constexpr (!ComputeUVs)
      t_norm = 0.0f;

    res = DistanceResult(dist - thickness, t_norm, dist, 0.0f, thickness);
  }

  /**
   * @brief distance<true>() from a precomputed pixel frame.
   * @param d Pixel dot ring axis (= dot(p, normal)).
   * @param polar fast_acos(clamp(d, -1, 1)).
   * @param sin_polar sqrtf(max(1 - d * d, POLE_SIN2_FLOOR)).
   * @param t_norm Pixel azimuth in [0, 1), phase applied.
   * @param res Output result, identical to distance<true>() except for
   *        undisplaced knot rings, which take the exact polar distance (the
   *        zero-knot polyline agrees only to within an ulp).
   * @details Same-axis ring stacks share d/polar/sin_polar/t_norm across every
   * ring at a pixel; hoisting them there drops the per-ring dot/acos/atan2
   * recompute.
   */
  HS_O3_FN void distance_from_frame(float d, float polar, float sin_polar,
                                    float t_norm, DistanceResult &res) const {
    if (d < cos_min_limit || d > cos_max_limit) {
      res = DistanceResult(100.0f, 0.0f, 100.0f, 0.0f, thickness);
      return;
    }
    float dist;
    if (knots)
      dist = max_distortion > 0.0f ? polyline_distance(t_norm, polar, sin_polar)
                                   : std::abs(polar - target_angle);
    else
      dist = std::abs(polar - (target_angle + shift_fn(t_norm)));
    res = DistanceResult(dist - thickness, t_norm, dist, 0.0f, thickness);
  }

  static constexpr float POLE_SIN2_FLOOR =
      1e-6f; /**< sin^2 floor keeping the chart's azimuth scale positive at the
                poles. */

private:
  static constexpr int PREFILTER_CHUNKS =
      32; /**< Azimuth chunks in the knot-range prefilter. */
  static constexpr int MAX_SEARCH_CELLS =
      64; /**< Outward search budget per side; only near-pole chart compression
             approaches it. */

  float chunk_lo[PREFILTER_CHUNKS]; /**< Min knot per azimuth chunk (knot mode).
                                     */
  float chunk_hi[PREFILTER_CHUNKS]; /**< Max knot per azimuth chunk (knot mode).
                                     */

  /**
   * @brief Distance from a pixel to the knot polyline.
   * @param t_norm Pixel azimuth in [0, 1) (phase applied).
   * @param polar Pixel polar angle (radians).
   * @param sin_polar sqrtf(max(1 - d * d, POLE_SIN2_FLOOR)) for the pixel;
   *        hoisted to the caller so a ring stack pays it once per pixel.
   * @return Geodesic distance to the nearest polyline point (radians), exact
   *         within the local tangent chart for distances up to `thickness`;
   *         beyond that reach only the > thickness ordering is preserved.
   * @details Works in the chart (azimuth * sin(polar), polar) centered on the
   * pixel: exact point-to-segment distances, searched outward from the
   * pixel's own cell. A segment o cells away is at least (o - 1) * cell_u
   * away, so the search stops once that gap exceeds the best distance found
   * (or the stroke reach, past which alpha is zero regardless). True distance
   * is 1-Lipschitz in screen position, so the stroke's alpha cannot ripple
   * along the centerline the way slope-corrected vertical-distance estimates
   * do at curvature extrema.
   */
  HS_O3_FN float polyline_distance(float t_norm, float polar,
                                   float sin_polar) const {
    const float base =
        target_angle - polar; // knot m sits at v = base + knots[m]

    // Prefilter: when a chunk's arc exceeds the stroke reach, only the pixel's
    // chunk and its neighbours can hold a within-reach curve point; a pixel
    // whose polar offset clears all three knot ranges by more than thickness
    // skips the segment search (most band pixels, in a displaced ring).
    const float chunk_u = (2.0f * PI_F / PREFILTER_CHUNKS) * sin_polar;
    if (chunk_u >= thickness) {
      int c = static_cast<int>(t_norm * PREFILTER_CHUNKS);
      if (c >= PREFILTER_CHUNKS)
        c = PREFILTER_CHUNKS - 1;
      int cl = c == 0 ? PREFILTER_CHUNKS - 1 : c - 1;
      int cr = c == PREFILTER_CHUNKS - 1 ? 0 : c + 1;
      float lo = std::min(chunk_lo[cl], std::min(chunk_lo[c], chunk_lo[cr]));
      float hi = std::max(chunk_hi[cl], std::max(chunk_hi[c], chunk_hi[cr]));
      float gap = std::max(base + lo, -(base + hi));
      if (gap > thickness)
        return gap;
    }

    float x = t_norm * lut_n;
    int j = static_cast<int>(x);
    if (j >= lut_n) // t_norm * lut_n can round up to lut_n at the seam
      j = lut_n - 1;
    float f = x - j;
    const float cell_u = (2.0f * PI_F / lut_n) * sin_polar;
    const float cell_u2 = cell_u * cell_u;

    // Chart v of knot m relative to the pixel (pixel at the origin).
    auto knot_v = [&](int m) {
      if (m >= lut_n)
        m -= lut_n;
      else if (m < 0)
        m += lut_n;
      return base + knots[m];
    };
    // Perpendicular foot on the segment rising cell_u wide from (u0, v0) to
    // v1; contributes only when the foot lands inside the segment, so the
    // division is paid roughly once per pixel.
    auto interior_d2 = [&](float u0, float v0, float v1, float &best2) {
      float dv = v1 - v0;
      float numer = -(u0 * cell_u + v0 * dv);
      float len2 = cell_u2 + dv * dv;
      if (numer > 0.0f && numer < len2) {
        float cross = u0 * dv - v0 * cell_u;
        best2 = std::min(best2, cross * cross / len2);
      }
    };

    // Straddling cell first, then knot-by-knot outward on each arm; endpoint
    // distances are shared between adjacent segments so each step loads one
    // new knot.
    float ul = -f * cell_u; // arm frontiers: knots j (left), j + 1 (right)
    float ur = (1.0f - f) * cell_u;
    float vl = knot_v(j);
    float vr = knot_v(j + 1);
    float best2 = std::min(ul * ul + vl * vl, ur * ur + vr * vr);
    interior_d2(ul, vl, vr, best2);

    const float th2 = thickness * thickness;
    float bound2 = std::min(best2, th2);
    const int max_o = std::min(MAX_SEARCH_CELLS, lut_n / 2 + 1);
    for (int o = 1; o <= max_o; ++o) {
      // A segment past a frontier knot can't beat that knot's |u|.
      bool right = ur * ur < bound2;
      bool left = ul * ul < bound2;
      if (!right && !left)
        break;
      if (right) {
        float un = ur + cell_u;
        float vn = knot_v(j + o + 1);
        best2 = std::min(best2, un * un + vn * vn);
        interior_d2(ur, vr, vn, best2);
        ur = un;
        vr = vn;
      }
      if (left) {
        float un = ul - cell_u;
        float vn = knot_v(j - o);
        best2 = std::min(best2, un * un + vn * vn);
        interior_d2(un, vn, vl, best2);
        ul = un;
        vl = vn;
      }
      bound2 = std::min(bound2, best2);
    }
    // Beyond the stroke reach only the > thickness ordering matters; skip
    // the sqrt for those pixels.
    return best2 >= th2 ? thickness : sqrtf(best2);
  }
};

/**
 * @brief Undisplaced DistortedRing geometry with exact polar distance.
 */
struct FlatDistortedRing : private DistortedRing {
  using DistortedRing::get_horizontal_intervals;
  using DistortedRing::get_vertical_bounds;
  using DistortedRing::is_solid;
  using DistortedRing::suppress_pole_fill;
  using DistortedRing::thickness;

  /**
   * @brief Builds an undisplaced ring using exact polar centerline distance.
   * @param b Orientation frame (v = ring axis).
   * @param r Ring radius as a fraction of the hemisphere.
   * @param th Half-width of the stroke (radians).
   * @param ph Azimuth phase offset (radians).
   */
  FlatDistortedRing(const Basis &b, float r, float th, float ph = 0.0f)
      : DistortedRing(b, r, th, ScalarFn{}, 0.0f, ph) {}

  /**
   * @brief Computes signed distance to the undisplaced ring.
   * @param p Point on sphere (normalized).
   * @return DistanceResult with dist, t, and raw_dist populated.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  /**
   * @brief Computes signed distance to the undisplaced ring, writing into res.
   * @tparam ComputeUVs When true, also computes the azimuthal t parameter.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed distance, raw_dist = unsigned
   *        centerline distance, t = azimuth in [0,1) when ComputeUVs.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float d = dot(p, normal);
    if (d < cos_min_limit || d > cos_max_limit) {
      res = DistanceResult(100.0f, 0.0f, 100.0f, 0.0f, thickness);
      return;
    }

    float polar = fast_acos(hs::clamp(d, -1.0f, 1.0f));
    float t_norm = 0.0f;
    if constexpr (ComputeUVs) {
      float azimuth = fast_atan2(dot(p, w), dot(p, u));
      if (azimuth < 0)
        azimuth += 2 * PI_F;
      t_norm = wrap_t((azimuth + phase) / (2 * PI_F));
    }

    float dist = std::abs(polar - target_angle);
    res = DistanceResult(dist - thickness, t_norm, dist, 0.0f, thickness);
  }
};

/**
 * @brief CSG Union operation (A + B), taking the minimum distance of two
 * shapes.
 * @tparam A First child shape type.
 * @tparam B Second child shape type.
 */
template <typename A, typename B> struct Union {
  const A &a;      /**< First child shape. */
  const B &b;      /**< Second child shape. */
  float thickness; /**< Max child thickness (drives AA falloff). */
  static constexpr bool is_solid =
      A::is_solid || B::is_solid; /**< Solid if either child is; the union keeps
                                     both interiors. */

  static_assert(SDFShape<A> && SDFShape<B>,
                "CSG Union children must be SDF shapes "
                "(is_solid/thickness)");
  static_assert(A::is_solid == B::is_solid,
                "CSG Union children must share solidity; a solid+stroke mix "
                "renders the stroke winner through the solid AA branch");
  static_assert(
      sdf_max_spans<A>::value + sdf_max_spans<B>::value <=
          2 * INTERVAL_SPAN_CAP,
      "nested CSG union exceeds MergedIntervalBuffer capacity; flatten "
      "the union or raise INTERVAL_SPAN_CAP");

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
    // A culled child returns {1,0}; folding its sentinel via min/max would
    // annex rows neither child occupies.
    bool a_culled = b1.y_min > b1.y_max;
    bool b_culled = b2.y_min > b2.y_max;
    if (a_culled && b_culled)
      return {1, 0};
    int lo = a_culled   ? b2.y_min
             : b_culled ? b1.y_min
                        : std::min(b1.y_min, b2.y_min);
    int hi = a_culled   ? b2.y_max
             : b_culled ? b1.y_max
                        : std::max(b1.y_max, b2.y_max);
    return {lo, hi};
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
    MergedIntervalBuffer merged;

    bool has_a = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { push_interval(merged, start, end); });

    bool has_b = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { push_interval(merged, start, end); });

    // One child fell back to full width: the whole row needs the full scan.
    if (!has_a || !has_b)
      return false;

    if (merged.is_empty())
      return true;

    // Emitted spans may straddle θ=0 and are not seam-normalized to [0,W): the
    // union merge is frame-tolerant, so scan_region's wrap+coalesce is the seam
    // authority (Subtract/Intersection normalize first because their pairwise
    // span comparison is frame-sensitive).
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
      return;
    res = res_b;
  }
};

/**
 * @brief Smooth CSG Union using a polynomial smooth minimum (Inigo Quilez
 * smin).
 * @tparam A First child shape type.
 * @tparam B Second child shape type.
 * @details Shapes organically blend together within radius k (radians).
 */
template <typename A, typename B> struct SmoothUnion {
  const A &a;      /**< First child shape. */
  const B &b;      /**< Second child shape. */
  float k;         /**< Smoothing radius in radians (e.g. 0.1). */
  float thickness; /**< Max child thickness (drives AA falloff). */
  static constexpr bool is_solid =
      A::is_solid || B::is_solid; /**< Solid if either child is; the smooth
                                     union keeps both interiors. */

  static_assert(SDFShape<A> && SDFShape<B>,
                "CSG SmoothUnion children must be SDF shapes "
                "(is_solid/thickness)");
  static_assert(A::is_solid == B::is_solid,
                "CSG SmoothUnion children must share solidity; a solid+stroke "
                "mix renders the stroke winner through the solid AA branch");
  static_assert(
      sdf_max_spans<A>::value + sdf_max_spans<B>::value <=
          2 * INTERVAL_SPAN_CAP,
      "nested CSG smooth-union exceeds MergedIntervalBuffer capacity; "
      "flatten the union or raise INTERVAL_SPAN_CAP");

  /**
   * @brief Builds a smooth union of two child shapes.
   * @param shape_a First child shape.
   * @param shape_b Second child shape.
   * @param smoothness Blend radius k in radians.
   */
  SmoothUnion(const A &shape_a, const B &shape_b, float smoothness)
      : a(shape_a), b(shape_b), k(smoothness),
        thickness(std::max(shape_a.thickness, shape_b.thickness)) {
    HS_CHECK(k > 0.0f);
  }

  /**
   * @brief Row bounds spanning both children's bands, padded by the blend
   * radius.
   * @tparam H Canvas height in rows.
   * @return Inclusive row bounds expanded by k (converted to rows).
   */
  template <int H> Bounds get_vertical_bounds() const {
    constexpr int H_VIRT = H + hs::H_OFFSET;
    auto b1 = a.template get_vertical_bounds<H>();
    auto b2 = b.template get_vertical_bounds<H>();
    // A culled child returns {1,0}; folding its sentinel via min/max would
    // annex rows neither child occupies.
    bool a_culled = b1.y_min > b1.y_max;
    bool b_culled = b2.y_min > b2.y_max;
    if (a_culled && b_culled)
      return {1, 0};
    int lo = a_culled   ? b2.y_min
             : b_culled ? b1.y_min
                        : std::min(b1.y_min, b2.y_min);
    int hi = a_culled   ? b2.y_max
             : b_culled ? b1.y_max
                        : std::max(b1.y_max, b2.y_max);
    // Expand by the blend radius k (radians) converted to rows: phi spans [0,π]
    // over (H_VIRT-1) rows.
    int pad = std::max(1, static_cast<int>(ceilf(k * (H_VIRT - 1) / PI_F)));
    return {std::max(0, lo - pad), std::min(H - 1, hi + pad)};
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
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    MergedIntervalBuffer merged;
    // Great-circle weld radius k spans k/sin(phi) columns of azimuth; the
    // equatorial conversion under-covers toward the poles. Clamp to full width
    // where the latitude factor diverges.
    float sin_phi = TrigLUT<W, H>::sin_phi[y];
    float pad_px =
        sin_phi > INTERVAL_DENOM_EPS
            ? std::min(k * W / (2 * PI_F) / sin_phi, static_cast<float>(W))
            : static_cast<float>(W);

    bool has_a = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          push_interval(merged, start - pad_px, end + pad_px);
        });

    bool has_b = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          push_interval(merged, start - pad_px, end + pad_px);
        });

    // If either child falls back to full-width, so must the blend.
    if (!has_a || !has_b)
      return false;

    if (merged.is_empty())
      return true;

    // Emitted spans may straddle θ=0 and are not seam-normalized to [0,W): the
    // union merge is frame-tolerant, so scan_region's wrap+coalesce is the seam
    // authority (Subtract/Intersection normalize first because their pairwise
    // span comparison is frame-sensitive).
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
   * @param res Output result; the nearer child's result with its dist reduced
   * by the cubic smin blend term.
   * @note Only `dist` is blended across the weld; the auxiliary registers
   *       (`t`/`raw_dist`/`size`/UVs) snap to the nearer child, so a shader
   * keying off them sees a hard edge through the weld. Intentional, not a bug.
   * @warning The cubic smin pulls `dist` below the true distance near the weld,
   *          so this SDF is not sphere-tracing-safe (unlike WarpedVolume's
   *          Lipschitz-corrected distance) — scanline rasterization only.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    a.template distance<ComputeUVs>(p, res);
    DistanceResult res_b;
    b.template distance<ComputeUVs>(p, res_b);

    // Polynomial smooth min (cubic).
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
  static constexpr bool is_solid =
      A::is_solid; /**< Tracks the minuend; carving B never changes A's
                      solidity. */

  static_assert(SDFShape<A> && SDFShape<B>,
                "CSG Subtract children must be SDF shapes "
                "(is_solid/thickness)");
  // Each child is collected into an IntervalBuffer (cap INTERVAL_SPAN_CAP)
  // before differencing, so a child that could emit more spans must be rejected
  // at compile time rather than trapping in push_interval at runtime.
  static_assert(sdf_max_spans<A>::value <= INTERVAL_SPAN_CAP &&
                    sdf_max_spans<B>::value <= INTERVAL_SPAN_CAP,
                "nested CSG Subtract child exceeds IntervalBuffer capacity; "
                "flatten the nesting or raise INTERVAL_SPAN_CAP");

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
    IntervalBuffer intervals_a;
    IntervalBuffer intervals_b;

    bool has_a = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          push_interval(intervals_a, start, end);
        });

    if (!has_a)
      return false;

    if (intervals_a.is_empty())
      return true;

    // The set-difference loop and scan_region's coalescer both require
    // start-sorted intervals; a multi-interval child can emit them out of
    // order.
    sort_intervals_by_start(intervals_a);

    // A stroke subtrahend's edge-bands can coalesce into one chord spanning the
    // stroke's hollow interior; subtracting that would carve A's interior. So
    // for a non-solid B, emit A's spans (seam-split into [0, W)) and let
    // per-pixel max(A, -B) carve exactly the stroke band — with no horizontal
    // culling, so Subtract<solid, stroke> pays full A-coverage shading.
    if constexpr (!B::is_solid) {
      constexpr size_t SEAM_SPLIT_CAP = 2 * INTERVAL_SPAN_CAP;
      static_assert(2 * sdf_max_spans<A>::value <= SEAM_SPLIT_CAP,
                    "post-seam-split span count exceeds norm buffer capacity");
      StaticCircularBuffer<std::pair<float, float>, SEAM_SPLIT_CAP> norm_a;
      normalize_intervals_to_range<W>(intervals_a, norm_a);
      sort_intervals_by_start(norm_a);
      for (size_t i = 0; i < norm_a.size(); ++i)
        out(norm_a[i].first, norm_a[i].second);
      return true;
    }

    bool has_b = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          push_interval(intervals_b, start, end);
        });

    // B's fallback means "no intervals produced", not "covers the row": without
    // B's intervals the set difference is undefined, so request a full-row scan
    // and let per-pixel max(A, -B) handle it. Emitting nothing would erase A.
    if (!has_b)
      return false;

    // B produced no intervals: it removes nothing, so pass A through raw.
    if (intervals_b.is_empty()) {
      for (size_t i = 0; i < intervals_a.size(); ++i)
        out(intervals_a[i].first, intervals_a[i].second);
      return true;
    }

    // Normalize both children into [0, W) (seam-split) before differencing: a
    // band straddling θ=0 can be emitted by A and B in different wrap frames,
    // and the raw-coordinate disjointness test below would then miss the
    // overlap and under-carve at the seam. Seam-splitting at most doubles each
    // child's span count, so the buffers are sized 2x.
    //
    // Output bound: the set difference splits an A span once per enclosed B
    // span (B spans disjoint), so at most |norm_a| + |norm_b| spans. Only one
    // span per child straddles θ=0 in a row, so the real maximum is
    // sdf_max_spans<A> + sdf_max_spans<B> + 2; reaching even
    // 2*INTERVAL_SPAN_CAP would require both children to emit INTERVAL_SPAN_CAP
    // disjoint arcs in a single row, which no SDF shape does. push_interval
    // traps (fail-fast) if that ever holds.
    constexpr size_t SEAM_SPLIT_CAP = 2 * INTERVAL_SPAN_CAP;
    static_assert(2 * sdf_max_spans<A>::value <= SEAM_SPLIT_CAP &&
                      2 * sdf_max_spans<B>::value <= SEAM_SPLIT_CAP,
                  "post-seam-split span count exceeds norm buffer capacity");
    StaticCircularBuffer<std::pair<float, float>, SEAM_SPLIT_CAP> norm_a;
    StaticCircularBuffer<std::pair<float, float>, SEAM_SPLIT_CAP> norm_b;
    normalize_intervals_to_range<W>(intervals_a, norm_a);
    normalize_intervals_to_range<W>(intervals_b, norm_b);

    // The set-difference loop and scan_region's coalescer both require
    // start-sorted intervals; the seam split above can reorder them.
    sort_intervals_by_start(norm_a);
    sort_intervals_by_start(norm_b);

    // Set difference: for each A interval, subtract all overlapping B
    // intervals.
    for (size_t ai = 0; ai < norm_a.size(); ++ai) {
      float cur_start = norm_a[ai].first;
      float cur_end = norm_a[ai].second;

      for (size_t bi = 0; bi < norm_b.size(); ++bi) {
        float bs = norm_b[bi].first;
        float be = norm_b[bi].second;

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
   * @param res Output result; max(A, -B) with B's distance negated when it
   * wins.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    a.template distance<ComputeUVs>(p, res);
    DistanceResult res_b;
    b.template distance<ComputeUVs>(p, res_b);
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
  static constexpr bool is_solid =
      A::is_solid && B::is_solid; /**< Solid iff both children are. */

  static_assert(SDFShape<A> && SDFShape<B>,
                "CSG Intersection children must be SDF shapes "
                "(is_solid/thickness)");
  // Each child is collected into an IntervalBuffer (cap INTERVAL_SPAN_CAP)
  // before the merge-sweep, so a child that could emit more spans must be
  // rejected at compile time rather than trapping in push_interval at runtime.
  static_assert(
      sdf_max_spans<A>::value <= INTERVAL_SPAN_CAP &&
          sdf_max_spans<B>::value <= INTERVAL_SPAN_CAP,
      "nested CSG Intersection child exceeds IntervalBuffer capacity; "
      "flatten the nesting or raise INTERVAL_SPAN_CAP");

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
    IntervalBuffer intervals_a;
    IntervalBuffer intervals_b;

    bool has_a = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          push_interval(intervals_a, start, end);
        });

    bool has_b = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          push_interval(intervals_b, start, end);
        });

    // A full-width child intersected with the other child is just the other
    // child's intervals, already collected above; replay the buffer.
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

    // Normalize both children into [0, W) (seam-split) before the merge sweep:
    // a band straddling θ=0 can be emitted by A and B in different wrap frames
    // (A as [-5, 5], B as [W-5, W+5]), and the raw-coordinate overlap below
    // would then miss the shared coverage and under-report at the seam.
    // Splitting into a common [0, W) frame makes the comparison correct.
    // Seam-splitting at most doubles each child's span count, so the buffers
    // are sized 2x.
    constexpr size_t SEAM_SPLIT_CAP = 2 * INTERVAL_SPAN_CAP;
    static_assert(2 * sdf_max_spans<A>::value <= SEAM_SPLIT_CAP &&
                      2 * sdf_max_spans<B>::value <= SEAM_SPLIT_CAP,
                  "post-seam-split span count exceeds norm buffer capacity");
    StaticCircularBuffer<std::pair<float, float>, SEAM_SPLIT_CAP> norm_a;
    StaticCircularBuffer<std::pair<float, float>, SEAM_SPLIT_CAP> norm_b;
    normalize_intervals_to_range<W>(intervals_a, norm_a);
    normalize_intervals_to_range<W>(intervals_b, norm_b);

    // The merge sweep requires both lists start-sorted; the seam split above
    // can reorder them.
    sort_intervals_by_start(norm_a);
    sort_intervals_by_start(norm_b);

    size_t idx_a = 0;
    size_t idx_b = 0;

    while (idx_a < norm_a.size() && idx_b < norm_b.size()) {
      auto iv_a = norm_a[idx_a];
      auto iv_b = norm_b[idx_b];

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
 *
 * UV semantics: distance() evaluates the child at the *folded* point, so the
 * child's UV registers (DistanceResult.t / Fragment::v0) are sector-local — t
 * measures azimuth within one sector and resets at every sector boundary (each
 * copy reuses the child's full UV range).
 */
template <typename Shape> struct AngularRepeat {
  const Shape &shape; /**< Child shape being repeated. */
  Vector axis, u,
      w; /**< Rotation axis and the derived perpendicular plane (u, w). */
  int repetitions; /**< Number of copies around the axis. */
  float thickness; /**< Inherited from the child shape. */
  static constexpr bool is_solid =
      Shape::is_solid; /**< Matches the child's solidity. */

  static_assert(SDFShape<Shape>, "AngularRepeat child must be an SDF shape "
                                 "(is_solid/thickness)");

  /**
   * @brief Repeats the shape around an arbitrary axis.
   * @param s Child shape to repeat.
   * @param reps Number of copies (must be > 0).
   * @param ax Rotation axis (unit length).
   */
  AngularRepeat(const Shape &s, int reps, const Vector &ax)
      : shape(s), axis(ax), repetitions(reps), thickness(s.thickness) {
    HS_CHECK(reps > 0);
    HS_CHECK(fabsf(ax.length() - 1.0f) < 1e-3f);
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
    // Only a Y-axis fold (axis.y near ±1) preserves latitude; any other axis
    // sweeps latitudes the child never occupies, so its band cannot bound the
    // copies and every row must be scanned.
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
   * @param res Output result of the child at the folded point. Note res's UV
   *        registers (t / azimuth) are sector-local: the child sees the folded
   *        point, so t spans one sector and resets at each sector boundary.
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

// --- Congruence-class canonical distance LUTs --------------------------------
// Every islamic mesh's faces are near-exact copies of a handful of canonical 2D
// shapes (gnomonic projection about each face's centroid is
// position-covariant). MeshOps bakes one signed-distance LUT per congruence
// class at spawn (core/mesh/mesh_classes.h); Scan::Mesh binds it per frame via
// bind_class_lut, and Face::distance serves sign-pure probes from a bilinear
// lookup instead of the exact per-edge walk.

/** Minimum squared normalized correlation for a valid class-LUT alignment;
 *  below this the face is too deformed and keeps the exact path. */
static constexpr float ALIGN_MIN_CORR_SQ = 0.25f;
/** Maximum per-vertex deviation from the aligned canonical shape (as a multiple
 *  of the LUT cell diagonal) before a face keeps the exact path. The facility
 *  fits only meshes that hold still per spawn (mesh_classes.h). */
static constexpr float ALIGN_MAX_DEV_DIAGS = 0.25f;

/**
 * @brief Canonical congruence-class signed-distance LUT, baked once per
 *        spawned mesh.
 * @details Distances are in canonical gnomonic plane units, quantized to int16
 * over the LUT box diameter (step ~1e-5 plane units). The domain covers the
 * canonical polygon's bounding box + BOUNDS_MARGIN_WIDE, matching the
 * max_dist_sq cull ring, so any probe surviving the cull lands in-domain.
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
 * @brief Accumulated complex correlation between a canonical polygon and a
 *        centered projection (see align_correlate).
 */
struct AlignCorr {
  float rr, ri; /**< Sum of canon_k * conj(z'_k) (real, imaginary). */
  float cc, zz; /**< Power terms: sum |canon_k|^2 and sum |z'_k|^2. */
};

/**
 * @brief Correlates a canonical polygon against a centered 2D projection under
 *        a cyclic vertex offset + optional reflection.
 * @tparam GetZ Accessor invoked as get_z(j, zx, zy), returning the centered
 *         projection vertex j.
 * @param canon_xy Canonical centered polygon, x/y pairs, canonical order.
 * @param count Vertex count.
 * @param vert_offset Projection index corresponding to canonical vertex 0.
 * @param reflected Mirror family: conjugate each vertex and walk the
 *        projection indices in reverse (a mirrored face winds the opposite way
 *        in a consistently-wound mesh).
 * @param get_z Centered-projection vertex accessor.
 * @return The correlation sums; the least-squares residual is
 *         cc + zz - 2*|r|, and the optimal rotation is r / |r|.
 * @details Single source for the correspondence convention — bake-time
 * clustering (mesh_classes.h) and the per-frame Face::bind_class_lut both
 * route through it, so the (offset, reflected) encoding cannot drift.
 */
template <typename GetZ>
inline AlignCorr align_correlate(const float *canon_xy, int count,
                                 int vert_offset, bool reflected, GetZ get_z) {
  AlignCorr a{0.0f, 0.0f, 0.0f, 0.0f};
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
    float cx = canon_xy[2 * k], cy = canon_xy[2 * k + 1];
    // canon * conj(z)
    a.rr += cx * zx + cy * zy;
    a.ri += cy * zx - cx * zy;
    a.cc += cx * cx + cy * cy;
    a.zz += zx * zx + zy * zy;
  }
  return a;
}

/**
 * @brief Bakes the signed point-to-polygon distance field of a canonical
 *        centered 2D polygon into an int16 grid.
 * @param poly_xy Centered polygon vertices, x/y pairs.
 * @param count Vertex count (>= 3).
 * @param n Grid resolution per axis (>= 2).
 * @param out Storage for n*n quantized samples.
 * @param lut Receives the domain/quantization parameters, with data = out.
 * @details Exact per-edge walk with crossing-test sign, over the bounding box
 * + BOUNDS_MARGIN_WIDE. Quantization scale is the box diameter (an upper bound
 * on any in-box distance: the polygon meets its own bounding box), giving a
 * step of ~1e-5 plane units — far below the interpolation bound.
 */
inline void build_canonical_distance_lut(const float *poly_xy, int count, int n,
                                         int16_t *out, ClassLut &lut) {
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
 * @brief Scratch buffer for Face computations to avoid allocations.
 */
struct FaceScratchBuffer {
  static constexpr int MAX_VERTS = 64; /**< Maximum vertices per face. */
  std::array<Vector, MAX_VERTS + 1>
      poly_2d; /**< Projected 2D polygon (+1 entry to avoid modulo). */
  std::array<Vector, MAX_VERTS> edge_vectors; /**< Per-edge 2D vectors. */
  std::array<float, MAX_VERTS>
      edge_lengths_sq;                  /**< Per-edge squared lengths. */
  std::array<Vector, MAX_VERTS> planes; /**< Per-edge great-circle normals. */
  std::array<std::pair<float, float>, 4>
      intervals;                       /**< Azimuth coverage intervals. */
  std::array<float, MAX_VERTS> thetas; /**< Per-vertex azimuth angles. */
  std::array<float, MAX_VERTS>
      inv_edge_lengths_sq; /**< Reciprocal squared edge lengths. */
  std::array<float, MAX_VERTS>
      inv_edge_j; /**< Reciprocal of each edge's y-component. */
  std::array<Vector, MAX_VERTS + 1>
      verts_3d; /**< 3D vertices (+1 wrap entry). */
  std::array<Vector, MAX_VERTS>
      edge_normals; /**< Per-edge normalized 3D normals. */

  /**
   * @brief Packed per-edge data for the cache-friendly distance() fallback.
   */
  struct EdgePacked {
    float vx, vy, ex, ey, inv_len_sq,
        inv_ej; /**< Edge origin, vector, reciprocals. */
    uint32_t key_vy,
        key_next_vy; /**< angle_key of this and the next vertex's y. */
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
};

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
 * accelerate rasterization.
 */
struct Face {
  Vector center; /**< Normalized face centroid (projection axis). */
  Vector basis_v, basis_u, basis_w; /**< Local tangent frame (v = center). */
  int count;                        /**< Vertex/edge count. */
  float thickness;                  /**< Edge half-width (radians). */
  float size;               /**< Inradius metric for AA normalization. */
  float radius = 0.0f;      /**< Circumradius in the 2D projection. */
  float max_dist = 0.0f;    /**< Cull radius (circumradius plus margin). */
  float max_dist_sq = 0.0f; /**< Squared cull radius. */

  std::span<Vector> poly_2d;      /**< Projected 2D polygon (+1 wrap entry). */
  std::span<Vector> edge_vectors; /**< Per-edge 2D vectors. */
  std::span<float> edge_lengths_sq;     /**< Per-edge squared lengths. */
  std::span<Vector> planes;             /**< Per-edge great-circle normals. */
  std::span<float> inv_edge_lengths_sq; /**< Reciprocal squared edge lengths. */
  std::span<float> inv_edge_j;    /**< Reciprocal of each edge's y-component. */
  std::span<Vector> verts_3d;     /**< 3D vertices (+1 wrap entry). */
  std::span<Vector> edge_normals; /**< Per-edge normalized 3D normals. */

  int y_min, y_max; /**< Inclusive vertical row bounds. */
  int build_height; /**< Canvas height the bounds were computed for. */
  std::span<std::pair<float, float>>
      intervals;   /**< Azimuth coverage intervals (radians). */
  bool full_width; /**< True when the face spans all columns. */
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

  /**
   * @brief Builds a face's projection, bounds, and edge data.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices from the pool.
   * @param th Edge half-width (radians).
   * @param scratch Reusable scratch storage backing the spans.
   * @param h_virt Virtual row count (height plus pole offset).
   * @param height Canvas height in rows.
   * @param clip Optional render clip used to tighten the face bounds.
   */
  HS_O3_FN Face(std::span<const Vector> vertices,
                std::span<const uint16_t> indices, float th,
                FaceScratchBuffer &scratch, int h_virt, int height,
                const ClipRegion *clip = nullptr)
      : thickness(th), build_height(height), full_width(true) {

    // Early vertical exit: a face whose latitude band (plus thickness margin)
    // maps to an empty row range can never be rasterized.
    const bool phi_culled = [&] {
      HS_PROFILE(face_phi_extent);
      return compute_phi_extent(vertices, indices, h_virt, height);
    }();
    if (phi_culled) {
      y_min = 1;
      y_max = 0;
      return;
    }

    count = indices.size();
    HS_CHECK(count > 0 && count <= FaceScratchBuffer::MAX_VERTS,
             "Face: vertex count must be in (0, MAX_VERTS]");

    { HS_PROFILE(face_project);
      setup_frame_and_polygon(vertices, indices, scratch);
    }

    // A fully collapsed polygon (signed area ~ 0, e.g. a hankin rosette at
    // angle 0 spiking out to each edge midpoint and back through its corner)
    // encloses no region, but its boundary would still rasterize as a ~1 px
    // AA line; cull it like the phi-extent reject. The residual of an exactly
    // collapsed face is float noise (< ~1e-6 of radius^2), orders of
    // magnitude under the thinnest real sliver a sweep draws, so the
    // threshold decision is identical sim/device. Stroked shapes keep
    // drawing their outline.
    if (thickness <= 0.0f) {
      float area2 = 0.0f;
      for (int i = 0; i < count; ++i)
        area2 +=
            poly_2d[i].x * poly_2d[i + 1].y - poly_2d[i + 1].x * poly_2d[i].y;
      if (fabsf(area2) < COLLAPSED_AREA_RATIO * radius * radius) {
        y_min = 1;
        y_max = 0;
        return;
      }
    }

    { HS_PROFILE(face_thetas);
      compute_thetas(scratch);
    }
    { HS_PROFILE(face_azimuth);
      compute_azimuth_intervals(scratch);
    }

    // Azimuth half of the cull, ahead of the bounds pass. Only
    // apply_pole_containment can still widen the coverage these intervals
    // describe, so the circumcircle guard must clear first. The row half stays
    // below, where the bounds are exact.
    if (clip && !pole_within_circumcircle() && clip_rejects_azimuth(*clip)) {
      y_min = 1;
      y_max = 0;
      return;
    }

    int planes_count;
    { HS_PROFILE(face_bounds);
      compute_inradius(scratch);

      // Vertical bounds via full arc-extrema + pole analysis. A vertex-only phi
      // span misses the great-circle edge bulge toward a pole, leaving
      // near-pole faces with unscanned rows; the arc-extrema path covers them.
      planes_count =
          compute_full_bounds(scratch, count, center, thickness, h_virt, height,
                              y_min, y_max);
    }

    edge_vectors = std::span<Vector>(scratch.edge_vectors.data(), count);
    edge_lengths_sq = std::span<float>(scratch.edge_lengths_sq.data(), count);
    inv_edge_lengths_sq =
        std::span<float>(scratch.inv_edge_lengths_sq.data(), count);
    inv_edge_j = std::span<float>(scratch.inv_edge_j.data(), count);
    planes = std::span<Vector>(scratch.planes.data(), planes_count);

    { HS_PROFILE(face_pole);
      apply_pole_containment(height);
    }

    // Whole-face clip cull: y_min/y_max and the azimuth coverage now match what
    // the scan would draw, so a face disjoint from the clip band yields no
    // in-band pixel.
    if (clip && clip_rejects(*clip)) {
      y_min = 1;
      y_max = 0;
      return;
    }

    { HS_PROFILE(face_edges);
      pack_edges(scratch);
      build_half_planes(scratch);
    }
    { HS_PROFILE(face_sectors);
      build_sectors(scratch);
    }
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
    return clip_rejects_azimuth(cr);
  }

  /**
   * @brief Horizontal half of the clip cull.
   * @param cr Clip region (display bounds plus render margin).
   * @return True when the face's azimuth coverage lies outside the render
   *         columns. A full-width face or an inactive x-clip never rejects.
   */
  bool clip_rejects_azimuth(const ClipRegion &cr) const {
    if (full_width)
      return false;
    const ClipRegion::XClip xc = cr.x_clip();
    if (!xc.active)
      return false;
    const int Wd = cr.w;
    const float pw = 1.25f * (2.0f * PI_F / static_cast<float>(Wd));
    const int band_len = ((xc.re - xc.rs) % Wd + Wd) % Wd;
    for (const auto &iv : intervals) {
      // Same radians->column mapping and AA pad get_horizontal_intervals scans
      // with, so the cull matches the emitted columns exactly.
      int a = static_cast<int>(floorf((iv.first - pw) * Wd / (2.0f * PI_F)));
      int b = static_cast<int>(ceilf((iv.second + pw) * Wd / (2.0f * PI_F)));
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
   * @return True when the phi extent plus thickness margin maps to an empty
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
    float margin_check = thickness + BOUNDS_MARGIN;

    Bounds rows =
        phi_bounds_to_rows(min_phi_check - margin_check,
                           max_phi_check + margin_check, h_virt, height);

    return rows.y_min > rows.y_max;
  }

  /**
   * @brief Builds the local tangent frame, gnomonic 2D projection, and 3D
   * arrays.
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
      Vector v2 = scratch.poly_2d[i + 1];

      Vector edge = v2 - v1;
      float edge_len_sq = dot(edge, edge);
      float t = 0.0f;
      if (edge_len_sq > 1e-9f) {
        t = dot(-v1, edge) / edge_len_sq;
        t = __builtin_fmaxf(0.0f, __builtin_fminf(1.0f, t));
      }
      Vector closest = v1 + edge * t;
      float d_line = closest.magnitude();

      min_edge_dist = __builtin_fminf(d_line, min_edge_dist);
    }
    size = __builtin_fmaxf(min_edge_dist, radius * MIN_SIZE_RADIUS_RATIO);
    linear_dist = size < 0.2f;
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
      ep.key_next_vy = angle_key(poly_2d[i + 1].y);
    }
    packed_edges = std::span<EdgePacked>(scratch.packed_edges.data(), count);
  }

  /**
   * @brief Detects a convex 2D projection and builds its edge half-planes.
   * @param scratch Scratch storage receiving half_planes.
   * @details For a convex polygon the signed distance is max over edges of the
   * half-plane distance: exact everywhere inside and in the edge slabs outside;
   * outside a vertex's normal cone it underestimates (line distance, not vertex
   * distance), which only softens the AA corner within its ~1-pixel band. A
   * concave, degenerate-edged, or wrongly-oriented polygon leaves convex false
   * and distance() on the exact walk.
   */
  __attribute__((always_inline)) void
  build_half_planes(FaceScratchBuffer &scratch) {
    float area2 = 0.0f;
    bool pos = false, neg = false;
    // The turn test carries the previous edge in registers, so the ring closes
    // without indexing edge_vectors through a modulo.
    const Vector *e1 = &edge_vectors[count - 1];
    float l1 = edge_lengths_sq[count - 1];
    for (int i = 0; i < count; ++i) {
      const Vector &a = poly_2d[i];
      const Vector &b = poly_2d[i + 1];
      area2 += a.x * b.y - b.x * a.y;
      const Vector &e2 = edge_vectors[i];
      float cr = e1->x * e2.y - e1->y * e2.x;
      // Relative turn test: |sin| > 1e-6 between edge directions.
      float scale = l1 * edge_lengths_sq[i];
      e1 = &e2;
      l1 = edge_lengths_sq[i];
      if (cr * cr > 1e-12f * scale) {
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
  bool bind_class_lut(const ClassLut *lut, const float *canon_xy,
                      int vert_offset, bool reflected) {
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
    {
      int j = vert_offset;
      for (int k = 0; k < count; ++k) {
        float zx = poly_2d[j].x - mx;
        float zy = poly_2d[j].y - my;
        if (reflected) {
          zy = -zy;
          if (--j < 0)
            j = count - 1;
        } else {
          if (++j == count)
            j = 0;
        }
        float ex = canon_xy[2 * k] - (c * zx - s * zy);
        float ey = canon_xy[2 * k + 1] - (s * zx + c * zy);
        float dev_sq = ex * ex + ey * ey;
        if (dev_sq > max_dev_sq)
          max_dev_sq = dev_sq;
      }
    }
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
        theta += 2 * PI_F;
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
      // fmodf can leave start_t at ~2*PI instead of ~0, producing a degenerate
      // [~2*PI, 2*PI] sliver below; snap to 0.
      if (start_t > 2 * PI_F - 1e-4f)
        start_t = 0.0f;
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

  /**
   * @brief Extends the vertical bounds when the face encircles a pole.
   * @param height Canvas height in rows.
   * @details If the face encircles the north or south pole, the vertical bounds
   * are extended to it and full-width azimuth coverage is forced.
   */
  __attribute__((always_inline)) void apply_pole_containment(int height) {
    if (center.y > 0.01f) {
      float inv_c = 1.0f / center.y;
      if (pole_inside_polygon(basis_u.y * inv_c, basis_w.y * inv_c)) {
        y_min = 0;
        full_width = true;
      }
    }
    // South pole (0, -1, 0)
    if (center.y < -0.01f) {
      float inv_c = 1.0f / -center.y;
      if (pole_inside_polygon(-basis_u.y * inv_c, -basis_w.y * inv_c)) {
        y_max = height - 1;
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
          min_phi = __builtin_fminf(fast_acos(hs::clamp(pty, -1.0f, 1.0f)),
                                    min_phi);
        if (cx1 < 0 && cx2 < 0)
          max_phi = __builtin_fmaxf(fast_acos(hs::clamp(-pty, -1.0f, 1.0f)),
                                    max_phi);
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
   * @param thickness Edge half-width (radians).
   * @param h_virt Virtual row count (height plus pole offset).
   * @param height Canvas height in rows.
   * @param y_min_out Output: first covered row.
   * @param y_max_out Output: last covered row.
   * @return The number of great-circle planes computed.
   */
  HS_O3_FN static int compute_full_bounds(FaceScratchBuffer &scratch, int count,
                                          const Vector &center, float thickness,
                                          int h_virt, int height, int &y_min_out,
                                          int &y_max_out) {
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
    float margin = thickness + BOUNDS_MARGIN;
    Bounds rows =
        phi_bounds_to_rows(min_phi - margin, max_phi + margin, h_virt, height);
    y_min_out = rows.y_min;
    y_max_out = rows.y_max;
    return planes_count;
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
    return {y_min, y_max};
  }

  /**
   * @brief Emits the face's azimuth-coverage intervals for a row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param out Sink accepting (float start, float end).
   * @return True if intervals were emitted; false (full scan) for full-width
   * faces.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int, OutputIt out) const {
    if (full_width)
      return false;
#ifdef HS_AA_AUDIT
    if (hs_aa::g_audit.full_scan)
      return false;
#endif
    // Pad by the AA fringe reach (one pixel, in the plane units distance()
    // reports, slightly wider than angular) plus the fast_atan2 slop in the
    // vertex thetas the intervals derive from.
    float pad = 1.25f * (2.0f * PI_F / W);
    for (const auto &iv : intervals) {
      float f_x1 = (iv.first - pad) * W / (2 * PI_F);
      float f_x2 = (iv.second + pad) * W / (2 * PI_F);
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
   * @brief Signed planar distance via the exact per-edge walk.
   * @param px Gnomonic x of the query point.
   * @param py Gnomonic y of the query point.
   * @return Signed distance in the tangent plane (negative inside).
   */
  HS_O3_FN float plane_dist_exact(float px, float py) const {
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
    return (inside ? -1.0f : 1.0f) * sqrtf(d);
  }

  /**
   * @brief Signed planar distance via the concave sector walk.
   * @param px Gnomonic x of the query point.
   * @param py Gnomonic y of the query point.
   * @return Signed distance in the tangent plane (negative inside).
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
  HS_O3_FN float plane_dist_sector(float px, float py) const {
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
    bool inside = cr * sector_sgn >= 0.0f;
    return (inside ? -1.0f : 1.0f) * sqrtf(d);
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
   * @param res Output result; dist = signed edge distance minus thickness,
   *        raw_dist = the signed edge distance, size = inradius (gnomonic
   *        plane units).
   * @note Distances live in the face's gnomonic tangent plane. Small faces
   *       (linear_dist) report the plane distance directly; large faces convert
   *       via fast_atan2(plane, 1). Do not treat raw_dist as a metric geodesic
   *       angle.
   */
  template <bool ComputeUVs = true>
  HS_O3_FN void distance(const Vector &p, DistanceResult &res) const {
    HS_SCAN_METRIC(hs::g_scan_metrics.pixels_tested++);
    HS_PROBE_TICK();
    HS_PROBE_COUNT(n_probe);
    HS_PROBE_MARK(hs_t);

    float cos_angle = dot(p, center);
    if (cos_angle <= 0.01f) {
      HS_SCAN_METRIC(hs::g_scan_metrics.pixels_culled++);
      HS_PROBE_SPAN(point, hs_t);
      HS_PROBE_COUNT(n_cull_cos);
      res = DistanceResult(100.0f, 0.0f, 100.0f, 0.0f, size);
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
      res = DistanceResult(100.0f, 0.0f, 100.0f, 0.0f, size);
      return;
    }
    HS_PROBE_SPAN(project, hs_t);

    float plane_dist;
    if (lut_data) {
      // Affine map into the canonical LUT grid, then a 4-tap bilinear fetch.
      // Only sign-pure cells at least one cell diagonal from the boundary are
      // served; the AA fringe and sign-unsafe cells fall back to the exact walk
      // on the TRUE per-frame edges.
      float fx = lut_ax * px + lut_bx * py + lut_cx;
      float fy = lut_ay * px + lut_by * py + lut_cy;
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
        HS_PROBE_SPAN(edge_lut, hs_t);
        HS_PROBE_COUNT(n_lut);
      } else {
        HS_SCAN_METRIC(hs::g_scan_metrics.exact_hits++);
        if (convex) {
          plane_dist = plane_dist_convex(px, py);
          HS_PROBE_SPAN(edge_convex, hs_t);
          HS_PROBE_COUNT(n_convex);
        } else {
          plane_dist = plane_dist_exact(px, py);
          HS_PROBE_SPAN(edge_exact, hs_t);
          HS_PROBE_COUNT(n_exact);
        }
      }
    } else {
      HS_SCAN_METRIC(hs::g_scan_metrics.exact_hits++);
      if (convex) {
        plane_dist = plane_dist_convex(px, py);
        HS_PROBE_SPAN(edge_convex, hs_t);
        HS_PROBE_COUNT(n_convex);
      } else if (sector_ok) {
        HS_SCAN_METRIC(hs::g_scan_metrics.sector_hits++);
        plane_dist = plane_dist_sector(px, py);
        HS_PROBE_SPAN(edge_sector, hs_t);
        HS_PROBE_COUNT(n_sector);
      } else {
        plane_dist = plane_dist_exact(px, py);
        HS_PROBE_SPAN(edge_exact, hs_t);
        HS_PROBE_COUNT(n_exact);
      }
    }

    // Small faces skip the plane->angle conversion: tan(angle) ~ angle to
    // within size^2/3 of the shading gradient (< 1.5% at the 0.2 threshold).
    float raw = linear_dist ? plane_dist : fast_atan2(plane_dist, 1.0f);
    res = DistanceResult(raw - thickness, 0.0f, raw, 0.0f, size);
    HS_PROBE_SPAN(pack, hs_t);
  }
};

/**
 * @brief Calculates signed distance to a planar polygon.
 * @details DistanceResult fields: dist = signed distance from edge (negative
 * inside); t = normalized polar distance (polar / thickness); raw_dist = polar
 * distance from center.
 */
struct PlanarPolygon {
  const Basis &basis; /**< Orientation frame (v = polygon axis). */
  float thickness;    /**< Polygon radius / apothem scale (radians). */
  int sides;          /**< Number of polygon sides. */
  float phase;        /**< Azimuth phase offset (radians). */
  float apothem;      /**< Precomputed inradius (radians). */
  float ny, R_val,
      alpha_angle; /**< Axis y-component, XZ projection length and azimuth. */
  float phi_min, phi_max; /**< Vertical bounds as an angular band (radians). */
  float sign;             /**< +1 fills the polygon, -1 its complement. */
  static constexpr bool is_solid =
      true; /**< Polygon renders as a filled region. */

  /**
   * @brief Builds a planar polygon from its basis, thickness, side count,
   * phase.
   * @param b Orientation frame (v = polygon axis).
   * @param th Polygon radius / apothem scale (radians).
   * @param s Number of polygon sides (must be > 0).
   * @param ph Azimuth phase offset (radians).
   * @param invert When true, fill the complement (a shape spanning more than a
   *        hemisphere, rendered via its antipodal fold).
   */
  PlanarPolygon(const Basis &b, float th, int s, float ph, bool invert = false)
      : basis(b), thickness(th), sides(s), phase(ph),
        sign(invert ? -1.0f : 1.0f) {
    HS_CHECK(sides > 0);
    apothem = thickness * cosf(PI_F / sides);
    AxisProjection ap = project_axis(basis.v);
    ny = ap.ny;
    R_val = ap.R_val;
    alpha_angle = ap.alpha_angle;

    float center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));
    float margin = thickness + BOUNDS_MARGIN_WIDE;
    phi_min = std::max(0.0f, center_phi - margin);
    phi_max = std::min(PI_F, center_phi + margin);
    if (invert) {
      phi_min = 0.0f;
      phi_max = PI_F;
    }
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
    if (sign < 0.0f) // the complement wraps every row
      return false;
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    // Pad the cap by one pixel so the outer AA fringe is scanned.
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
   * @note The `polar*cos(local)` form under-estimates the true distance near
   *       the sector corners (gradient < 1 there), like the tangent-plane
   *       caveat on the public distance() above; for scanline shading, not a
   *       march-safe metric.
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

    res = DistanceResult(sign * dist_edge, t_val, polar, 0.0f, apothem);
  }
};

/**
 * @brief Calculates signed distance to a spherical polygon (great-circle
 * edges).
 * @details Uses sector folding plus a precomputed great-circle plane normal for
 * O(1) per-pixel distance, with exact angular distances for smooth AA.
 * DistanceResult fields: dist = signed angular distance to the nearest
 * great-circle edge (negative inside); t = normalized radial position (polar /
 * circumradius); raw_dist = polar angle from the polygon center.
 */
struct SphericalPolygon {
  const Basis &basis; /**< Orientation frame (v = polygon axis). */
  int sides;          /**< Number of polygon sides. */
  float phase;        /**< Azimuth phase offset (radians). */
  float circumradius; /**< Angular distance from center to vertex (radians). */
  float edge_nv;      /**< Edge normal dotted with the center axis. */
  float edge_nu;      /**< Edge normal dotted with the u-axis. */
  float phi_min, phi_max; /**< Vertical bounds as an angular band (radians). */
  float ny, R_val,
      alpha_angle; /**< Axis y-component, XZ projection length and azimuth. */
  float sign;      /**< +1 fills the polygon, -1 its complement. */
  static constexpr bool is_solid =
      true; /**< Polygon renders as a filled region. */

  /**
   * @brief Builds a spherical polygon from its basis, radius, side count,
   * phase.
   * @param b Orientation frame (v = polygon axis).
   * @param radius Polygon radius as a fraction of the hemisphere.
   * @param s Number of polygon sides (must be > 0).
   * @param ph Azimuth phase offset (radians).
   * @param invert When true, fill the complement (a shape spanning more than a
   *        hemisphere, rendered via its antipodal fold).
   */
  SphericalPolygon(const Basis &b, float radius, int s, float ph,
                   bool invert = false)
      : basis(b), sides(s), phase(ph), sign(invert ? -1.0f : 1.0f) {
    HS_CHECK(sides > 0);
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
      // Degenerate canonical edge (circumradius near 0 or PI): both limits
      // drive en toward ±u, so substitute the sector bisector as a defined unit
      // normal.
      en = basis.u;
    }
    // Ensure outward: dot(center, n) should be negative
    if (dot(en, basis.v) > 0)
      en = -en;

    edge_nv = dot(en, basis.v);
    edge_nu = dot(en, basis.u);

    // Vertical/horizontal bounds (same as SDF::PlanarPolygon)
    AxisProjection ap = project_axis(basis.v);
    ny = ap.ny;
    R_val = ap.R_val;
    alpha_angle = ap.alpha_angle;

    float center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));
    float margin = circumradius + BOUNDS_MARGIN_WIDE;
    phi_min = std::max(0.0f, center_phi - margin);
    phi_max = std::min(PI_F, center_phi + margin);
    if (invert) {
      phi_min = 0.0f;
      phi_max = PI_F;
    }
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
    if (sign < 0.0f) // the complement wraps every row
      return false;
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    // Pad the cap by one pixel so the outer AA fringe is scanned.
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
    res = DistanceResult(sign * dist_edge, t_val, polar, 0.0f, circumradius);
  }
};

/**
 * @brief Calculates signed distance to a star shape.
 * @details DistanceResult fields: dist = signed distance from edge (negative
 * inside); t = normalized azimuth (azimuth / 2PI); raw_dist = polar distance
 * from center.
 */
struct Star {
  const Basis &basis; /**< Orientation frame (v = star axis). */
  int sides;          /**< Number of star points. */
  float phase;        /**< Azimuth phase offset (radians). */
  static constexpr bool is_solid =
      true; /**< Star renders as a filled region. */

  float nx, ny,
      plane_d;     /**< 2D edge plane (normal and offset) for one point. */
  float thickness; /**< Outer radius / AA scale (radians). */

  float scan_ny, scan_r,
      scan_alpha; /**< Axis y-component, XZ projection length and azimuth. */
  float phi_min, phi_max; /**< Vertical bounds as an angular band (radians). */
  float sign;             /**< +1 fills the star, -1 its complement. */

  /**
   * @brief Builds a star from its basis, radius, point count, and phase.
   * @param b Orientation frame (v = star axis).
   * @param radius Outer radius as a fraction of the hemisphere.
   * @param s Number of star points (must be > 0).
   * @param ph Azimuth phase offset (radians).
   * @param invert When true, fill the complement (a shape spanning more than a
   *        hemisphere, rendered via its antipodal fold).
   */
  Star(const Basis &b, float radius, int s, float ph, bool invert = false)
      : basis(b), sides(s), phase(ph), sign(invert ? -1.0f : 1.0f) {
    HS_CHECK(sides > 0);
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

    AxisProjection ap = project_axis(basis.v);
    scan_ny = ap.ny;
    scan_r = ap.R_val;
    scan_alpha = ap.alpha_angle;

    float center_phi = acosf(std::max(-1.0f, std::min(1.0f, basis.v.y)));
    float margin = outer_radius + BOUNDS_MARGIN_WIDE;
    phi_min = std::max(0.0f, center_phi - margin);
    phi_max = std::min(PI_F, center_phi + margin);
    if (invert) {
      phi_min = 0.0f;
      phi_max = PI_F;
    }
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
    if (sign < 0.0f) // the complement wraps every row
      return false;
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
    res = DistanceResult(sign * -dist_to_edge, t, scan_dist, 0.0f, thickness);
  }
};

/**
 * @brief Calculates signed distance to a flower shape.
 * @details DistanceResult fields: dist = signed distance from the flower edge;
 * t = normalized scan distance (scan_dist / thickness); raw_dist = scan
 * distance from the antipode.
 */
struct Flower {
  const Basis &basis; /**< Orientation frame (v = flower axis). */
  int sides;          /**< Number of petals. */
  float phase;        /**< Azimuth phase offset (radians). */
  float thickness;    /**< Outer radius / AA scale (radians). */
  float apothem;      /**< Petal inradius offset (PI - outer radius). */
  Vector antipode;    /**< Antipode of the flower axis (scan origin). */
  float scan_ny, scan_R, scan_alpha; /**< Antipode y-component, XZ projection
                                        length and azimuth. */
  float phi_min, phi_max; /**< Vertical bounds as an angular band (radians). */
  float sign;             /**< +1 fills the flower, -1 its complement. */
  static constexpr bool is_solid =
      true; /**< Flower renders as a filled region. */

  /**
   * @brief Builds a flower from its basis, radius, petal count, and phase.
   * @param b Orientation frame (v = flower axis).
   * @param radius Outer radius as a fraction of the hemisphere.
   * @param s Number of petals (must be > 0).
   * @param ph Azimuth phase offset (radians).
   * @param invert When true, fill the complement (a shape spanning more than a
   *        hemisphere, rendered via its antipodal fold).
   */
  Flower(const Basis &b, float radius, int s, float ph, bool invert = false)
      : basis(b), sides(s), phase(ph), sign(invert ? -1.0f : 1.0f) {
    HS_CHECK(sides > 0);
    float outer = radius * (PI_F / 2.0f);
    apothem = PI_F - outer;
    thickness = outer;
    antipode = -basis.v;

    AxisProjection ap = project_axis(antipode);
    scan_ny = ap.ny;
    scan_R = ap.R_val;
    scan_alpha = ap.alpha_angle;

    float center_phi = acosf(std::max(-1.0f, std::min(1.0f, antipode.y)));
    float margin = thickness + BOUNDS_MARGIN_WIDE;
    phi_min = std::max(0.0f, center_phi - margin);
    phi_max = std::min(PI_F, center_phi + margin);
    if (invert) {
      phi_min = 0.0f;
      phi_max = PI_F;
    }
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
    if (sign < 0.0f) // the complement wraps every row
      return false;
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

    res = DistanceResult(sign * -dist_edge, t_val, scan_dist, 0.0f, thickness);
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

  // Bounding-cap geometry, loop-invariant across scanlines (precomputed in
  // ctor).
  Vector mid; /**< Arc midpoint axis (bounding-cap center). */
  float mid_ny = 0.0f, mid_r = 0.0f,
        mid_alpha = 0.0f; /**< Midpoint y, XZ projection length, azimuth. */
  float cap_D_min = 0.0f; /**< Cosine of the bounding-cap radius. */
  bool cap_horiz_valid = false; /**< False when the cap axis is ~vertical (no
                                   horizontal projection). */
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
      // Antipodal endpoints (len ~ π) make cross() degenerate; guard the
      // normalize.
      n = normalized_or(cross(a, b), Vector(1, 0, 0));
    }

    // A great-circle arc bulges toward a pole between its endpoints, so its
    // peak latitude can lie outside [phi_a, phi_b]; extend by the interior
    // extremum.
    float phi_a = acosf(std::max(-1.0f, std::min(1.0f, a.y)));
    float phi_b = acosf(std::max(-1.0f, std::min(1.0f, b.y)));
    float phi_lo = std::min(phi_a, phi_b);
    float phi_hi = std::max(phi_a, phi_b);
    // The latitude turns inside the arc iff the forward tangent's y-component
    // (cross(n, p).y) flips sign between endpoints; the extremum is
    // ±sqrt(1-n.y²).
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

    // Bounding cap centered on the segment midpoint. Antipodal endpoints sum to
    // ~0 (no defined midpoint); guard the normalize.
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
    // Two geodesic metrics, C0-continuous at the join: in-segment the foot lies
    // on the arc so asinf(|d_plane|) is exact; past an endpoint the closest
    // point is that cap (asinf(|d_plane|) alone would measure the whole great
    // circle).
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

    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    float cos_phi = TrigLUT<W, H>::cos_phi[y];
    float sin_phi = TrigLUT<W, H>::sin_phi[y];

    return emit_cap_interval<W>(cap_D_min, mid_ny, mid_r, mid_alpha, cos_phi,
                                sin_phi, /*reject_full_width=*/false, out);
  }
};

// ============================================================================
// 3D Volumetric SDF Shapes (for Scan::Volume raymarching)
// ============================================================================
// These shapes are reached only from Scan::Volume::draw's march loop (Raymarch
// is their sole instantiation site), so they share that region's -O3 options —
// an -Os distance() here would be a wall inside the promoted loop.

HS_O3_BEGIN
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
   * @brief Outward normal direction at a point, from a precomputed XZ radius.
   * @param p Query point in Cartesian ray-space.
   * @param inv_xz_len 1 / sqrtf(p.x² + p.z²) for `p`, or 0 on the Y axis.
   * @return The outward normal direction, NOT unit length; the caller
   *         normalizes, or folds `p` through a linear map that it then
   *         normalizes.
   */
  Vector normal_raw(const Vector &p, float inv_xz_len) const {
    float scale = R * inv_xz_len;
    return p - Vector(p.x * scale, 0.0f, p.z * scale);
  }

  /**
   * @brief Surface normal at a point near the torus surface.
   * @param p Query point in Cartesian ray-space.
   * @return Unit outward normal.
   */
  Vector normal(const Vector &p) const {
    float xz_len = sqrtf(p.x * p.x + p.z * p.z);
    float inv_xz_len = (xz_len > TOLERANCE) ? 1.0f / xz_len : 0.0f;
    return normal_raw(p, inv_xz_len).normalized();
  }

  /**
   * @brief Populates a Fragment's registers for shading.
   * @param p Query point in Cartesian ray-space.
   * @param frag Output fragment; v0 = ring angle (0-1, for palette lookup),
   *        v1/v2/v3 = surface normal (x, y, z).
   * @note Volumetric register convention (README "Volumetric Path"), distinct
   *        from Scan's v2 stroke-coverage and mesh face-index conventions.
   */
  void populate(const Vector &p, Fragment &frag) const {
    Vector n = normal(p);
    frag.v0 = (fast_atan2(p.z, p.x) + PI_F) / (2.0f * PI_F);
    frag.v1 = n.x;
    frag.v2 = n.y;
    frag.v3 = n.z;
  }
};
HS_O3_END

/**
 * @brief Domain warp functions for composing with WarpedSDF.
 */
namespace Warp {

HS_O3_BEGIN
/**
 * @brief Oscillates the Y coordinate sinusoidally around the azimuthal
 * angle θ = atan2(z, x).
 *
 * Produces twisted/undulating geometry when composed with a torus.
 * Provides an analytical Lipschitz bound for safe sphere tracing.
 */
struct Twist {
  int twist;           /**< Number of oscillations around the ring. */
  float amplitude;     /**< Vertical displacement magnitude. */
  float R;             /**< Major radius (needed for the Lipschitz bound). */
  float twist_amp;     /**< twist * amplitude, the warp's angular gradient. */
  float twist_amp_abs; /**< |twist * amplitude|, the Lipschitz numerator. */
  float two_over_r;    /**< 2/R, the Lipschitz reciprocal clamp. */

  /**
   * @brief Constructs a twist warp around a torus of major radius R.
   * @param twist_ Number of oscillations around the ring.
   * @param amplitude_ Vertical displacement magnitude.
   * @param R_ Major radius; must be > 0. The Lipschitz bound scales by 2/R, so
   *        R == 0 yields a non-finite bound on the XZ axis. Guarded at the cold
   *        construction site, not per-call.
   */
  Twist(int twist_, float amplitude_, float R_)
      : twist(twist_), amplitude(amplitude_), R(R_),
        twist_amp(static_cast<float>(twist_) * amplitude_),
        twist_amp_abs(fabsf(twist_amp)), two_over_r(2.0f / R_) {
    HS_CHECK(R > 0.0f);
  }

  /** @brief Precomputed context: s = sqrtf(x² + z²), shared across
   * apply/lipschitz. */
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
   * @param s Shared context from make_ctx(p).
   * @return The warped point.
   */
  Vector apply(const Vector &p, Ctx s) const {
    return Vector(p.x, p.y - amplitude * sin_ntheta(p, s), p.z);
  }

  /**
   * @brief sin(n*theta) at theta = atan2(z, x), without evaluating either.
   * @param p Query point.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @return sin(twist * theta).
   * @details Chebyshev recurrence sin((k+1)t) = 2cos(t)sin(kt) - sin((k-1)t),
   * seeded from (cos t, sin t) = (x/s, z/s), so the cost is one reciprocal
   * rather than an atan2 and a sine. Exact to float rounding, where
   * fast_atan2/fast_sinf each carry approximation error.
   */
  float sin_ntheta(const Vector &p, Ctx s) const {
    return sin_ntheta_inv(p, s).sin_n;
  }

  /** @brief sin(twist*theta) with the reciprocal radius that seeded it. */
  struct SinInv {
    float sin_n; /**< sin(twist * theta). */
    float inv_s; /**< 1/s; 2/R where the recurrence degenerates. */
  };

  /**
   * @brief sin(n*theta) and the reciprocal 1/s from the same recurrence seed.
   * @param p Query point.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @return sin(twist * theta) and 1/s, the latter feeding lipschitz().
   * @details On the degenerate axis inv_s carries 2/R, the value the Lipschitz
   * clamp would select there anyway, so no caller needs a second branch.
   */
  SinInv sin_ntheta_inv(const Vector &p, Ctx s) const {
    if (twist == 0 || s < TOLERANCE)
      return {0.0f, two_over_r};
    const float inv_s = 1.0f / s;
    const float two_cos = 2.0f * p.x * inv_s;
    float prev = 0.0f, cur = p.z * inv_s;
    for (int k = 1; k < twist; ++k) {
      const float next = two_cos * cur - prev;
      prev = cur;
      cur = next;
    }
    return {cur, inv_s};
  }

  /**
   * @brief cos(n*theta) at theta = atan2(z, x), by the same recurrence.
   * @param p Query point.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @return cos(twist * theta).
   */
  float cos_ntheta(const Vector &p, Ctx s) const {
    if (twist == 0 || s < TOLERANCE)
      return 1.0f;
    const float inv_s = 1.0f / s;
    const float two_cos = 2.0f * p.x * inv_s;
    float prev = 1.0f, cur = p.x * inv_s;
    for (int k = 1; k < twist; ++k) {
      const float next = two_cos * cur - prev;
      prev = cur;
      cur = next;
    }
    return cur;
  }

  /** @brief Both harmonics of twist*theta from one recurrence pass. */
  struct SinCos {
    float sin_n; /**< sin(twist * theta). */
    float cos_n; /**< cos(twist * theta). */
  };

  /**
   * @brief sin(n*theta) and cos(n*theta) together, sharing one recurrence
   *        setup and loop.
   * @param p Query point.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @return Both harmonics, matching sin_ntheta/cos_ntheta bit for bit.
   * @details For the hit path, which needs both: the march path needs only the
   * sine and calls sin_ntheta so it does not pay for the cosine sequence.
   */
  SinCos sincos_ntheta(const Vector &p, Ctx s) const {
    return sincos_ntheta_inv(p, (s > TOLERANCE) ? 1.0f / s : 0.0f);
  }

  /**
   * @brief sin(n*theta) and cos(n*theta) from an already-computed 1/s.
   * @param p Query point.
   * @param inv_s Reciprocal of the XZ radius; 0 marks the degenerate axis.
   * @return Both harmonics, matching sincos_ntheta bit for bit.
   */
  SinCos sincos_ntheta_inv(const Vector &p, float inv_s) const {
    if (twist == 0 || inv_s == 0.0f)
      return {0.0f, 1.0f};
    const float two_cos = 2.0f * p.x * inv_s;
    float sin_prev = 0.0f, sin_cur = p.z * inv_s;
    float cos_prev = 1.0f, cos_cur = p.x * inv_s;
    for (int k = 1; k < twist; ++k) {
      const float sin_next = two_cos * sin_cur - sin_prev;
      sin_prev = sin_cur;
      sin_cur = sin_next;
      const float cos_next = two_cos * cos_cur - cos_prev;
      cos_prev = cos_cur;
      cos_cur = cos_next;
    }
    return {sin_cur, cos_cur};
  }

  /**
   * @brief Analytical Lipschitz constant of the warp at a point.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @return The exact operator norm of the warp Jacobian (>= 1).
   */
  float lipschitz(const Vector & /*p*/, Ctx s) const {
    return lipschitz(1.0f / std::max(s, R * 0.5f));
  }

  /**
   * @brief Analytical Lipschitz constant from an already-computed 1/s.
   * @param inv_s Reciprocal of the XZ radius, from sin_ntheta_inv().
   * @return The exact operator norm of the warp Jacobian (>= 1).
   * @details The warp Jacobian is the shear I - e_y·gᵀ with e_y ⊥ g and
   * |g| = γ; its operator norm (largest singular value) is γ/2 + √(1 + γ²/4).
   * γ uses |twist·amplitude| so the bound stays conservative regardless of
   * sign, and min(1/s, 2/R) is the clamp 1/max(s, R/2).
   */
  float lipschitz(float inv_s) const {
    if (twist == 0)
      return 1.0f;
    const float gamma = twist_amp_abs * std::min(inv_s, two_over_r);
    return 0.5f * gamma + sqrtf(1.0f + 0.25f * gamma * gamma);
  }

  /**
   * @brief Reciprocal of the Lipschitz constant, from an already-computed 1/s.
   * @param inv_s Reciprocal of the XZ radius, from sin_ntheta_inv().
   * @return 1 / lipschitz(inv_s), in (0, 1].
   * @details (√(1+γ²/4) - γ/2)(√(1+γ²/4) + γ/2) = 1 exactly, so the reciprocal
   * needs no divide. The result never exceeds 1, so callers scale by it
   * unconditionally.
   */
  float lipschitz_inv(float inv_s) const {
    if (twist == 0)
      return 1.0f;
    const float gamma = twist_amp_abs * std::min(inv_s, two_over_r);
    return sqrtf(1.0f + 0.25f * gamma * gamma) - 0.5f * gamma;
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
    return correct_normal(p, base_n, s, cos_ntheta(p, s));
  }

  /**
   * @brief Normal correction from an already-computed cos(twist*theta).
   * @param p Query point.
   * @param base_n Unwarped surface normal; need not be unit length.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @param cos_n cos(twist * theta) at `p`.
   * @return The corrected unit normal.
   * @details The correction is a linear map of base_n, so scaling base_n scales
   * the result and the final normalize cancels it — an unnormalized base normal
   * gives the identical unit result and saves a normalize.
   */
  Vector correct_normal(const Vector &p, const Vector &base_n, Ctx s,
                        float cos_n) const {
    return correct_normal_inv(p, base_n, (s > TOLERANCE) ? 1.0f / s : 0.0f,
                              cos_n);
  }

  /**
   * @brief Normal correction from an already-computed 1/s and cos(twist*theta).
   * @param p Query point.
   * @param base_n Unwarped surface normal; need not be unit length.
   * @param inv_s Reciprocal of the XZ radius; 0 marks the degenerate axis.
   * @param cos_n cos(twist * theta) at `p`.
   * @return The corrected unit normal.
   */
  Vector correct_normal_inv(const Vector &p, const Vector &base_n, float inv_s,
                            float cos_n) const {
    float dh_dtheta = -twist_amp * cos_n;
    float inv_s2 = inv_s * inv_s;

    float dh_dx = dh_dtheta * (-p.z) * inv_s2;
    float dh_dz = dh_dtheta * p.x * inv_s2;

    return Vector(base_n.x - base_n.y * dh_dx, base_n.y,
                  base_n.z - base_n.y * dh_dz)
        .normalized();
  }
};
HS_O3_END

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
 *   Vector apply(const Vector &p, Ctx ctx) const;
 *   float lipschitz(const Vector &p, const Ctx &ctx) const;
 *   float bounding_inflation() const;
 * Optionally:
 *   Vector correct_normal(const Vector &p, const Vector &n, const Ctx&) const;
 */
HS_O3_BEGIN
template <typename SDF, typename Warp> struct WarpedVolume {
  SDF base;  /**< The underlying volume SDF being warped. */
  Warp warp; /**< The domain warp applied before the base SDF. */

  /**
   * @brief Smallest distance the caller needs an accurate value for.
   * @details The cheap bound is returned only strictly above this, so a coarse
   * value can never reach a caller's hit test or antialiasing band. Zero means
   * unstated and selects the warp's maximum displacement, which is safe for any
   * caller but forfeits most of the fast path.
   */
  float precision = 0.0f;

  /** True when the base/warp pair admits the tight per-axis bound below. */
  static constexpr bool TORUS_TWIST = std::is_same_v<SDF, ::SDF::Torus> &&
                                      std::is_same_v<Warp, ::SDF::Warp::Twist>;

  /**
   * @brief Cheap lower bound on the warped distance.
   * @param p Query point in Cartesian ray-space.
   * @return A lower bound on the distance to the warped surface.
   */
  float bounding_distance(const Vector &p) const {
    if constexpr (TORUS_TWIST) {
      // Twist moves only y, by at most bounding_inflation(), so the warped
      // surface lies inside the torus swept +-A along y; this is that solid's
      // exact distance, hence a lower bound. Relaxing the y term alone keeps it
      // far tighter than subtracting A from the whole distance.
      const float q = sqrtf(p.x * p.x + p.z * p.z) - base.R;
      const float dy = std::max(fabsf(p.y) - warp.bounding_inflation(), 0.0f);
      return sqrtf(q * q + dy * dy) - base.r;
    } else {
      return base.distance(p) - warp.bounding_inflation();
    }
  }

  /**
   * @brief Raw warped SDF distance with no Lipschitz correction.
   * @param p Query point in Cartesian ray-space.
   * @return The base SDF evaluated at the warped point (use for surface
   * projection).
   */
  float raw_distance(const Vector &p) const {
    auto ctx = warp.make_ctx(p);
    return base.distance(warp.apply(p, ctx));
  }

  /**
   * @brief March-safe distance with bounding fast-path and Lipschitz
   * correction.
   * @param p Query point in Cartesian ray-space.
   * @return A sphere-tracing-safe (under-estimated) distance to the surface.
   */
  float distance(const Vector &p) const {
    const float gate =
        (precision > 0.0f) ? precision : warp.bounding_inflation();
    if constexpr (TORUS_TWIST) {
      // gate + base.r > 0, so bounding_distance(p) > gate is exactly
      // qq > (gate + base.r)²; the sqrt is then only the fast path's result.
      const float s = sqrtf(p.x * p.x + p.z * p.z);
      const float q = s - base.R;
      const float dy = std::max(fabsf(p.y) - warp.bounding_inflation(), 0.0f);
      const float qq = q * q + dy * dy;
      const float t = gate + base.r;
      if (qq > t * t)
        return sqrtf(qq) - base.r;

      const auto h = warp.sin_ntheta_inv(p, s);
      const Vector warped(p.x, p.y - warp.amplitude * h.sin_n, p.z);
      float d = base.distance(warped);
      if (d > 0.0f)
        d *= warp.lipschitz_inv(h.inv_s);
      return d;
    } else {
      const float bd = base.distance(p) - warp.bounding_inflation();
      if (bd > gate)
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
  }

  /**
   * @brief Surface normal with the warp's chain-rule correction.
   * @param p Query point in Cartesian ray-space.
   * @return The corrected unit surface normal.
   */
  Vector normal(const Vector &p) const {
    auto ctx = warp.make_ctx(p);
    if constexpr (TORUS_TWIST) {
      // Twist displaces only y, so the warped point keeps p's XZ radius `ctx`
      // and the torus normal needs no second sqrt; one recurrence yields both
      // the sine the warp needs and the cosine the correction needs; and the
      // correction normalizes, so the base normal can stay unnormalized.
      // twist == 0 needs no special case: it gives cos_n = 1, hence a zero
      // gradient and an unchanged normal.
      const float inv_s = (ctx > TOLERANCE) ? 1.0f / ctx : 0.0f;
      if (warp.amplitude < TOLERANCE)
        return base.normal_raw(p, inv_s).normalized();
      auto h = warp.sincos_ntheta_inv(p, inv_s);
      Vector warped(p.x, p.y - warp.amplitude * h.sin_n, p.z);
      return warp.correct_normal_inv(p, base.normal_raw(warped, inv_s), inv_s,
                                     h.cos_n);
    } else {
      Vector base_n = base.normal(warp.apply(p, ctx));
      if constexpr (requires { warp.correct_normal(p, base_n, ctx); }) {
        return warp.correct_normal(p, base_n, ctx);
      }
      return base_n;
    }
  }

  /**
   * @brief Populates a Fragment's registers for shading.
   * @param p Query point in Cartesian ray-space.
   * @param frag Output fragment; v0 = ring angle (0-1), v1/v2/v3 = surface
   * normal.
   * @note Volumetric register convention (README "Volumetric Path"), distinct
   *        from Scan's v2 stroke-coverage and mesh face-index conventions.
   */
  void populate(const Vector &p, Fragment &frag) const {
    Vector n = normal(p);
    frag.v0 = (fast_atan2(p.z, p.x) + PI_F) / (2.0f * PI_F);
    frag.v1 = n.x;
    frag.v2 = n.y;
    frag.v3 = n.z;
  }
};
HS_O3_END

} // namespace SDF
