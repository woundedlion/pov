/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include "math/geometry.h"
#include "platform/constants.h"
#include "engine/memory.h"
#include "containers/static_circular_buffer.h"

/**
 * @file common.h
 * @brief The scanline contract every SDF shape shares: azimuth intervals and
 * their span bounds, row bounds, the DistanceResult register table and the
 * cap/annular-band emission helpers the leaves stroke rows with.
 */

namespace SDF {

// --- Rasterization constants ------------------------------------------------
/** Margin added to bounding boxes for AA kernel width (small faces). */
inline constexpr float BOUNDS_MARGIN = 0.05f;
/** Expanded margin for shapes with variable thickness or large faces. */
inline constexpr float BOUNDS_MARGIN_WIDE = 0.1f;
/** Minimum horizontal projection length for interval culling. */
inline constexpr float MIN_HORIZONTAL_PROJ = 0.01f;
/** Epsilon for near-zero denominators in interval math. */
inline constexpr float INTERVAL_DENOM_EPS = 1e-6f;
/** Threshold for near-pole ring approximation safety. */
inline constexpr float POLE_SAFE_MARGIN = 0.05f;
/** Ring thickness ceiling as a fraction of tan(target_angle) for the
 *  linearized centerline distance; bounds its relative error at half this. */
inline constexpr float RING_LINEARIZE_TAN_FRAC = 0.1f;
/** Pole-on-boundary tolerance, as a fraction of a face's 2D circumradius. */
inline constexpr float POLE_BOUNDARY_TOL = 1e-3f;
/** Inner/outer radius ratio for star shapes. */
inline constexpr float STAR_INNER_RATIO = ::STAR_INNER_RATIO;
/** Minimum inradius-to-circumradius ratio used to floor Face::size,
 *  preventing degenerate near-zero inradii from collapsing AA. */
inline constexpr float MIN_SIZE_RADIUS_RATIO = 0.25f;
/** Distance reported in place of a true one past a shape's reject band, far
 *  enough that no AA reach or CSG blend reads it as near-surface. The shapes
 *  that clamp to it are exactly those with blends_smoothly == false. */
inline constexpr float FAR_SENTINEL = 100.0f;

/** @brief Folds an angle into the centered interval for one sector. */
inline float centered_sector_angle(float angle, float sector,
                                   float reciprocal_sector) {
  float shifted = angle + sector * 0.5f;
  return shifted - floorf(shifted * reciprocal_sector) * sector - sector * 0.5f;
}

/**
 * @brief Azimuth of a point in a basis' u-w plane, plus a phase offset.
 * @param p Point on sphere (normalized).
 * @param u Basis u axis.
 * @param w Basis w axis.
 * @param phase Azimuth phase offset (radians), added after the fold.
 * @return The azimuth folded into [0, 2*PI), offset by phase.
 */
__attribute__((always_inline)) inline float
basis_azimuth(const Vector &p, const Vector &u, const Vector &w, float phase) {
  float azimuth = fast_atan2(dot(p, w), dot(p, u));
  if (azimuth < 0)
    azimuth += TWO_PI_F;
  return azimuth + phase;
}

/** Signed-area-to-circumradius-squared ratio below which a Face is culled as
 *  fully collapsed (no enclosed region). Sits orders of magnitude above the
 *  float noise of an exactly collapsed polygon (~1e-7) and below the thinnest
 *  real sliver a mesh sweep draws (~1e-3), so the sim/device decision is
 *  identical under fast-math. */
inline constexpr float COLLAPSED_AREA_RATIO = 1e-5f;

/** Squared relative-turn epsilon for the convexity test: a turn registers only
 *  when |sin| between successive edge directions exceeds 1e-6. Shared with
 *  MeshOps::polygon_is_concave so a class is LUT-eligible exactly when its
 *  faces miss Face::build_half_planes' convex fast path. */
inline constexpr float TURN_EPS_SQ = 1e-12f;

/** AA fringe pad in radians applied to a face's azimuth intervals: one pixel
 *  of falloff reach (in the plane units distance() reports, slightly wider
 *  than angular) plus the fast_atan2 slop in the vertex thetas the intervals
 *  derive from.
 *  @param w Canvas width in columns.
 *  @return Half-width pad added to each end of an azimuth interval. */
constexpr float face_azimuth_pad(int w) {
  return 1.25f * (TWO_PI_F / static_cast<float>(w));
}

/**
 * @brief Converts the face AA reach to azimuth at one latitude.
 * @param w Canvas width in columns.
 * @param sin_phi Sine of the row's colatitude.
 * @return Azimuth half-width, or pi when the reach spans the whole row.
 */
inline float face_azimuth_pad(int w, float sin_phi) {
  const float pad = face_azimuth_pad(w);
  return sin_phi <= pad ? PI_F : asinf(pad / sin_phi);
}

// Scanline interval protocol. get_horizontal_intervals returns true when the
// spans it emitted describe the row, and false to request a full-row scan. A
// false return MUST emit nothing: the caller walks every column instead, so a
// span emitted before the fallback is either dropped or shaded twice.
//
// The protocol has no per-scan entry point, so composite scratch setup and the
// LUT initialization guards repeat per row.

/** Maximum disjoint scanline spans a single shape (leaf) emits per row.
 *  scan_region's `intervals` buffer holds a top-level CSG emission, and its
 *  seam-split `norm` buffer is 2x intervals (one span can split in two at the
 *  x=0 seam); render/scan/raster.h statically checks both capacities. */
inline constexpr size_t INTERVAL_SPAN_CAP = 32;

/** A scanline span [first, second) in fractional column units. An aggregate, so
 *  a span buffer's slots are left uninitialized when the buffer is constructed;
 *  std::pair's default constructor value-initializes every slot, which the
 *  per-draw and per-row buffers below would pay on every construction. */
struct Interval {
  float first;  /**< Span start column. */
  float second; /**< Span end column. */
};

/** Per-row scanline interval buffer for a single shape. Fixed capacity,
 *  accumulate-only. */
using IntervalBuffer = StaticCircularBuffer<Interval, INTERVAL_SPAN_CAP>;

/** Spans AngularRepeat may emit in one row before it falls back to a full
 *  scan. Fixed independently of the child so a repeat's compile-time span bound
 *  stays inside the CSG budget whatever it wraps. */
inline constexpr size_t ANGULAR_REPEAT_SPAN_CAP = 8;

/** Azimuth slop, in radians, between the point AngularRepeat's fold hands the
 *  child and the exact sector-shifted one: fast_atan2's ~0.0038 rad plus the
 *  fast_sinf/fast_cosf reconstruction of the folded direction. */
inline constexpr float ANGULAR_REPEAT_FOLD_SLOP = 0.01f;

/** Largest squared off-axis magnitude (|axis x Y|^2 = x^2 + z^2) that still
 *  counts as a Y-axis fold for AngularRepeat's cull: a tilt of 1e-4 rad, whose
 *  induced azimuth displacement stays inside ANGULAR_REPEAT_FOLD_SLOP even on
 *  the near-pole rows where the 1 / sin(phi) amplification peaks. */
inline constexpr float ANGULAR_REPEAT_Y_AXIS_TOL_SQ = 1e-8f;

/** Per-row accumulator for a binary CSG op that merges BOTH children's spans
 *  into one buffer (Union/SmoothUnion). Each child can contribute up to
 *  INTERVAL_SPAN_CAP spans, so the union accumulator is sized to hold both. */
using MergedIntervalBuffer =
    StaticCircularBuffer<Interval, 2 * INTERVAL_SPAN_CAP>;

/**
 * @brief Places a per-row CSG span buffer in a guarded scratch arena.
 * @tparam Buf Buffer type to construct.
 * @param scratch Open scope over the arena; its rewind reclaims the buffer.
 * @return Reference to the fresh buffer, dead once `scratch` exits.
 * @details The combinators run under scan_region on the deepest render chain,
 * where the device DTCM stack is tight, so their spans go to the arena for the
 * same reason scan_region's do. Arena scopes are LIFO, so a nested combinator's
 * buffers rewind above its parent's.
 */
template <typename Buf> inline Buf &scratch_spans(ScratchScope &scratch) {
  Arena &arena = scratch.get_arena();
  return *arena.make<Buf>();
}

// Forward-declared so the span-count trait below can pattern-match the binary
// CSG ops and the leaves (all defined later in this header) before their full
// definitions.
template <typename A, typename B> struct Union;
template <typename A, typename B> struct SmoothUnion;
template <typename A, typename B> struct Subtract;
template <typename A, typename B> struct Intersection;
struct Ring;
struct DistortedRing;
struct FlatDistortedRing;
struct Face;
struct PlanarPolygon;
struct SphericalPolygon;
struct Star;
struct Flower;
struct Line;
template <typename Shape> struct AngularRepeat;

/** Compile-time upper bound on the scanline spans a shape may emit to its
 * parent in one row. An unrecognized shape falls back to INTERVAL_SPAN_CAP, the
 * runtime buffer capacity; the leaves below are pinned to what their
 * get_horizontal_intervals can actually emit. Union/SmoothUnion merge both
 * children into one MergedIntervalBuffer, so their bound is the SUM of the
 * children's, static_asserted against the buffer capacity to reject an
 *  overflowing nesting at compile time. */
template <typename T> struct sdf_max_spans {
  static constexpr size_t value = INTERVAL_SPAN_CAP;
};

// Per-leaf emissions: the number of out() calls one row can make. These bound
// only how deeply the CSG combinators may nest -- every runtime buffer stays at
// INTERVAL_SPAN_CAP (or its 2x / 2x+2 derivatives), so a leaf is free to be
// re-sized here without touching storage.
//
// emit_annular_band: two arcs on a general row, collapsing to one when the band
// touches the near or far pole.
template <> struct sdf_max_spans<Ring> {
  static constexpr size_t value = 2;
};
template <> struct sdf_max_spans<DistortedRing> {
  static constexpr size_t value = 2;
};
template <> struct sdf_max_spans<FlatDistortedRing> {
  static constexpr size_t value = 2;
};
// emit_cap_interval: a single bounding-cap arc, or a full-scan request.
template <> struct sdf_max_spans<PlanarPolygon> {
  static constexpr size_t value = 1;
};
template <> struct sdf_max_spans<Flower> {
  static constexpr size_t value = 1;
};
template <> struct sdf_max_spans<SphericalPolygon> {
  static constexpr size_t value = 1;
};
template <> struct sdf_max_spans<Star> {
  static constexpr size_t value = 1;
};
template <> struct sdf_max_spans<Line> {
  static constexpr size_t value = 1;
};
// Face replays its azimuth-coverage span, which always views
// FaceScratchBuffer::intervals; tied to that array's size where it is defined.
template <> struct sdf_max_spans<Face> {
  static constexpr size_t value = 2;
};
// AngularRepeat replays the child's spans once per copy, capped independently
// of the child (past the cap it requests a full scan and emits nothing).
template <typename Shape> struct sdf_max_spans<AngularRepeat<Shape>> {
  static constexpr size_t value = ANGULAR_REPEAT_SPAN_CAP;
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
// the two start-sorted lists (one span per advance, |norm_a| + |norm_b| - 1
// advances). A child's spans carry no common wrap frame -- a Union can merge
// two that both cover θ=0 -- so every span may split: bound 2·|A| + 2·|B|.
template <typename A, typename B> struct sdf_max_spans<Intersection<A, B>> {
  static constexpr size_t value =
      2 * sdf_max_spans<A>::value + 2 * sdf_max_spans<B>::value;
};
// Subtract emits the minuend seam-split into a [0, W) frame. A child's spans
// carry no common wrap frame -- a Union can merge two that both cover θ=0 --
// so every span may split: bound 2·|A|.
template <typename A, typename B> struct sdf_max_spans<Subtract<A, B>> {
  static constexpr size_t value = 2 * sdf_max_spans<A>::value;
};

/** True when a shape's distance() reports a usable signed distance outside its
 * surface, which is what SmoothUnion's weld term needs. Ring, DistortedRing,
 * FlatDistortedRing and Face instead clamp to FAR_SENTINEL past their reject
 * band. A combinator blends only if every child does; an unrecognized shape
 * is rejected. */
template <typename T>
inline constexpr bool blends_smoothly =
    requires { T::BLENDS_SMOOTHLY; } && T::BLENDS_SMOOTHLY;
template <> inline constexpr bool blends_smoothly<PlanarPolygon> = true;
template <> inline constexpr bool blends_smoothly<SphericalPolygon> = true;
template <> inline constexpr bool blends_smoothly<Star> = true;
template <> inline constexpr bool blends_smoothly<Flower> = true;
template <> inline constexpr bool blends_smoothly<Line> = true;
template <> inline constexpr bool blends_smoothly<Ring> = false;
template <> inline constexpr bool blends_smoothly<DistortedRing> = false;
template <> inline constexpr bool blends_smoothly<FlatDistortedRing> = false;
template <> inline constexpr bool blends_smoothly<Face> = false;
template <typename Shape>
inline constexpr bool blends_smoothly<AngularRepeat<Shape>> =
    blends_smoothly<Shape>;
template <typename A, typename B>
inline constexpr bool blends_smoothly<Union<A, B>> =
    blends_smoothly<A> && blends_smoothly<B>;
template <typename A, typename B>
inline constexpr bool blends_smoothly<SmoothUnion<A, B>> =
    blends_smoothly<A> && blends_smoothly<B>;
template <typename A, typename B>
inline constexpr bool blends_smoothly<Intersection<A, B>> =
    blends_smoothly<A> && blends_smoothly<B>;
template <typename A, typename B>
inline constexpr bool blends_smoothly<Subtract<A, B>> =
    blends_smoothly<A> && blends_smoothly<B>;

/** Signed distance below which a shape's distance() always reports a true
 * value. At or above it the shape may substitute FAR_SENTINEL, an unbounded
 * jump that says nothing about a neighbouring probe, so a walk that vouches for
 * a run of columns from one probe may do so only while the clearance it tests
 * stays under this margin. A shape that never substitutes reports FLT_MAX via
 * REJECT_MARGIN; a combinator takes the tightest of its children. */
template <typename T>
inline constexpr float reject_margin = [] {
  if constexpr (requires { T::REJECT_MARGIN; })
    return T::REJECT_MARGIN;
  return 0.0f;
}();
template <> inline constexpr float reject_margin<PlanarPolygon> = FLT_MAX;
template <> inline constexpr float reject_margin<SphericalPolygon> = FLT_MAX;
template <> inline constexpr float reject_margin<Star> = FLT_MAX;
template <> inline constexpr float reject_margin<Flower> = FLT_MAX;
template <> inline constexpr float reject_margin<Line> = FLT_MAX;
// The bounding annulus is the stroke band itself, so every probe outside the
// stroke is a candidate for the sentinel.
template <> inline constexpr float reject_margin<Ring> = 0.0f;
template <> inline constexpr float reject_margin<DistortedRing> = 0.0f;
template <> inline constexpr float reject_margin<FlatDistortedRing> = 0.0f;
// The cull disk clears the polygon by BOUNDS_MARGIN_WIDE of gnomonic-plane
// distance; a large face reports that as atan(BOUNDS_MARGIN_WIDE) radians.
template <> inline constexpr float reject_margin<Face> = 0.0996f;
template <typename Shape>
inline constexpr float reject_margin<AngularRepeat<Shape>> =
    reject_margin<Shape>;
template <typename A, typename B>
inline constexpr float reject_margin<Union<A, B>> =
    std::min(reject_margin<A>, reject_margin<B>);
template <typename A, typename B>
inline constexpr float reject_margin<SmoothUnion<A, B>> =
    std::min(reject_margin<A>, reject_margin<B>);
template <typename A, typename B>
inline constexpr float reject_margin<Intersection<A, B>> =
    std::min(reject_margin<A>, reject_margin<B>);
template <typename A, typename B>
inline constexpr float reject_margin<Subtract<A, B>> =
    std::min(reject_margin<A>, reject_margin<B>);

/** Stands for a distance() with no finite change-per-arc factor. */
inline constexpr float ARC_STRETCH_UNBOUNDED = FLT_MAX;

/** Most a shape's distance() can change per unit of great-circle arc, over the
 * band within a pixel or two of its surface -- the only band a walk that
 * vouches for a run of columns from one probe has to cross. Such a walk scales
 * the run's arc by this; against ARC_STRETCH_UNBOUNDED no slack suffices and
 * the run must be walked per column. A combinator takes the loosest child. The
 * default covers reporting in plane units, which runs slightly wider than
 * angular. */
template <typename T> inline constexpr float arc_stretch = 1.25f;
// Sector fold in the azimuthal-equidistant chart: the azimuth term carries
// polar/sin(polar). In the band the circumscribed-disc clamp holds polar within
// a few columns of the circumradius, itself <= PI/2 once a radius past a
// hemisphere folds to the antipode, and polar/sin(polar) stays under 2 out to
// 1.89 rad.
template <> inline constexpr float arc_stretch<PlanarPolygon> = 2.0f;
template <> inline constexpr float arc_stretch<Star> = 2.0f;
// Flower measures polar from the antipode of the point it folds about and has
// no disc clamp, so its azimuth term carries (PI - scan_dist)/sin(scan_dist)
// and the fold axis itself is on the surface: every petal meets there.
template <> inline constexpr float arc_stretch<Flower> = ARC_STRETCH_UNBOUNDED;
template <typename Shape>
inline constexpr float arc_stretch<AngularRepeat<Shape>> = arc_stretch<Shape>;
template <typename A, typename B>
inline constexpr float arc_stretch<Union<A, B>> =
    std::max(arc_stretch<A>, arc_stretch<B>);
template <typename A, typename B>
inline constexpr float arc_stretch<SmoothUnion<A, B>> =
    std::max(arc_stretch<A>, arc_stretch<B>);
template <typename A, typename B>
inline constexpr float arc_stretch<Intersection<A, B>> =
    std::max(arc_stretch<A>, arc_stretch<B>);
template <typename A, typename B>
inline constexpr float arc_stretch<Subtract<A, B>> =
    std::max(arc_stretch<A>, arc_stretch<B>);

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
inline void push_interval(StaticCircularBuffer<Interval, N> &buf, float start,
                          float end) {
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
inline void sort_intervals_by_start(StaticCircularBuffer<Interval, N> &buf) {
  const size_t n = buf.size();
  if (n < 2)
    return;
  HS_CHECK(buf.is_linear(),
           "sort_intervals_by_start: raw linear indexing requires head==0");
  auto *data = &buf[0];
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
inline void
normalize_intervals_to_range(const StaticCircularBuffer<Interval, N> &src,
                             StaticCircularBuffer<Interval, M> &dst) {
  constexpr float Wf = static_cast<float>(W);
  for (size_t i = 0; i < src.size(); ++i) {
    float len = src[i].second - src[i].first;
    HS_CHECK(len >= 0.0f,
             "normalize_intervals_to_range: interval end precedes start");
    if (len >= Wf) {
      push_interval(dst, 0.0f, Wf);
      continue;
    }
    const float s = wrap(src[i].first, Wf);
    const float e = s + len;
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
inline void merge_intervals(StaticCircularBuffer<Interval, N> &merged,
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
  x = fabsf(x);           // cos(-x) = cos(x): fold negatives
  x = fmodf(x, TWO_PI_F); // 2π-periodic -> [0, 2π)
  if (x > PI_F)
    x = TWO_PI_F - x; // reflect (π, 2π) across the south pole -> (0, π)
  return x;
}

/** @brief A colatitude band as inclusive [phi_min, phi_max] bounds in [0, π].
 */
struct PhiBand {
  float phi_min, phi_max;
};

/**
 * @brief Exact colatitude extent of the circle of angular radius
 * `target_angle` about an axis at colatitude `center_phi`.
 * @param center_phi Axis colatitude in [0, π] (radians).
 * @param target_angle Angular radius of the circle (radians, any sign or
 * magnitude).
 * @return {phi_min, phi_max}: the tight band the circle occupies.
 * @details On the circle cos φ = cos(center_phi)cos(target_angle) −
 * sin(center_phi)sin(target_angle)cos ψ, so cos φ sweeps exactly
 * [cos(center_phi + target_angle), cos(center_phi − target_angle)] and the two
 * folded endpoints are its extremes; min/max orders them for a target_angle
 * outside [0, π]. Single source for the Ring/DistortedRing
 * get_vertical_bounds latitude fold so the two cannot drift apart.
 */
inline PhiBand clamp_phi_band(float center_phi, float target_angle) {
  float p1 = clamp_phi(center_phi - target_angle);
  float p2 = clamp_phi(center_phi + target_angle);
  return {std::min(p1, p2), std::max(p1, p2)};
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
 * `dist` and `size` have fixed meanings; `t`, `raw_dist` and `aux` are
 * per-shape registers whose authoritative meanings are the table below. The
 * scan rasterizer copies them into the Fragment register file with no
 * reinterpretation (see Scan::process_pixel): t -> Fragment::v0, raw_dist ->
 * Fragment::v1, aux -> Fragment::v3, size -> Fragment::size. Fragment::v2 is
 * generated downstream: Scan writes stroke AA coverage or 0 for solid shapes,
 * and Scan::Mesh replaces it with a face index.
 *
 * Per-producer register semantics (a leaf built with ComputeUVs = false
 * reports t = 0; a bounds/cull miss reports dist = raw_dist = FAR_SENTINEL
 * with t = 0):
 *
 * | Producer          | dist (negative inside)     | t                | raw_dist         | size |
 * |-------------------|----------------------------|------------------|------------------|------|
 * | Ring              | centerline dist - thickness | azimuth turns in [0,1), phase applied | unsigned centerline distance | thickness |
 * | DistortedRing / FlatDistortedRing | displaced-centerline dist - thickness | azimuth turns in [0,1), phase applied | unsigned distance to the displaced centerline | thickness |
 * | Line              | segment dist - thickness   | 0                | unsigned angular distance to the arc segment | thickness |
 * | Face              | signed edge distance; gnomonic plane units on a small face, radians on a large one (see Face::distance) | 0 | = dist | inradius, in dist's metric |
 * | PlanarPolygon     | sign * max(polar*cos(local) - apothem, polar - circumradius) | polar / circumradius | polar angle from center | apothem |
 * | SphericalPolygon  | sign * max(angular distance to the nearest great-circle edge, polar - circumradius) | polar / circumradius | polar angle from center | circumradius |
 * | Star              | sign * max(folded edge half-plane distance, scan_dist - circumradius) | azimuth turns in [0,1), phase applied | polar angle from center | circumradius |
 * | Flower            | sign * -(polar*cos(local) - apothem) | scan_dist / circumradius | scan distance from the antipode | circumradius |
 *
 * Combinators forward a child's registers: Union/Intersection keep the
 * nearer/farther child's result verbatim; SmoothUnion keeps the nearer
 * child's and blends only dist; Subtract keeps the winner's, negating dist
 * and holding size to the minuend's when B wins; AngularRepeat reports the
 * child at the folded point (t is sector-local).
 *
 * Every current producer writes aux = 0; it is a pass-through register
 * reserved for shader-visible per-shape data, riding to Fragment::v3.
 */
struct DistanceResult {
  float dist; /**< Signed distance (negative inside); always this meaning. */
  float t;    /**< Per-shape register; see the table above. */
  float raw_dist;    /**< Per-shape register; see the table above. */
  float aux;         /**< Per-shape register; see the table above. */
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
 * @brief Evaluates a shape's two-argument distance() into a fresh result.
 * @tparam S Shape type exposing distance<ComputeUVs>(const Vector&,
 *         DistanceResult&).
 * @param shape Shape to sample.
 * @param p Point on sphere (normalized).
 * @return The shape's DistanceResult at p, with UVs computed.
 */
template <typename S>
inline DistanceResult distance_of(const S &shape, const Vector &p) {
  DistanceResult res;
  shape.template distance<true>(p, res);
  return res;
}

/**
 * @brief Structural fingerprint of a CSG-composable SDF shape: a static
 * is_solid flag.
 * @tparam T Candidate shape type.
 * @details is_solid selects the rasterizer's AA path, so a composite must
 * expose one; the CSG combinators assert this concept on their children so a
 * wrong-type argument fails at the boundary. distance() and the scanline
 * members vary by render path and are not part of the shared contract.
 */
template <typename T>
concept SDFShape = requires {
  { T::is_solid } -> std::convertible_to<bool>;
};

/**
 * @brief The scan rasterizer's leaf contract: SDFShape plus the three methods
 * Scan::rasterize drives — vertical bounds, per-row horizontal intervals and a
 * distance evaluator.
 * @tparam T Candidate shape type.
 * @tparam W Canvas width in columns.
 * @tparam H Canvas height in rows.
 * @details Asserted by Scan::rasterize so a shape missing one of the three
 * fails at the rasterizer boundary rather than inside scan_region. The CSG
 * combinators stay on SDFShape: a child can be driven through another render
 * path (or, in tests, exercise the interval sweep alone) and carry only
 * is_solid.
 */
template <typename T, int W, int H>
concept ScanShape =
    SDFShape<T> && requires(const T &s, const Vector &p, DistanceResult &r,
                            void (*out)(float, float)) {
      { s.template get_vertical_bounds<H>() } -> std::convertible_to<Bounds>;
      {
        s.template get_horizontal_intervals<W, H>(0, out)
      } -> std::convertible_to<bool>;
      s.template distance<true>(p, r);
    };

/**
 * @brief Axis components plus its scan-plane projection: XZ-projection length
 * r_val and azimuth alpha_angle.
 */
struct AxisProjection {
  float nx, ny, nz, r_val, alpha_angle;
};

/**
 * @brief Decomposes an axis into components and its scan-plane projection.
 * @param axis The axis vector.
 * @return Components nx/ny/nz, XZ-projection length r_val, azimuth alpha_angle.
 */
inline AxisProjection project_axis(const Vector &axis) {
  return {axis.x, axis.y, axis.z, sqrtf(axis.x * axis.x + axis.z * axis.z),
          atan2f(axis.z, axis.x)};
}

/**
 * @brief Scan-plane projection of a cap axis plus the latitude band the cap
 * spans.
 */
struct CapBounds {
  float ny, r_val,
      alpha_angle; /**< Axis y-component, XZ projection length and azimuth. */
  float phi_min, phi_max; /**< Vertical bounds as an angular band (radians). */
  float cos_radius, sin_radius; /**< Cap radius trig, for the scanline pad. */
};

/**
 * @brief Projects a cap axis and derives its margin-widened latitude band.
 * Shared by the leaf shapes bounded by a cap around a single axis.
 * @param axis The cap axis (normalized).
 * @param radius Angular radius of the cap (radians).
 * @param invert When true the shape fills the complement, which touches every
 *        row, so the band opens to the whole sphere.
 * @return Axis projection, the bounding band widened by BOUNDS_MARGIN_WIDE, and
 *         the cap radius' cosine and sine.
 */
inline CapBounds cap_bounds(const Vector &axis, float radius, bool invert) {
  AxisProjection ap = project_axis(axis);
  float center_phi = acosf(std::max(-1.0f, std::min(1.0f, ap.ny)));
  float margin = radius + BOUNDS_MARGIN_WIDE;
  return {ap.ny,
          ap.r_val,
          ap.alpha_angle,
          invert ? 0.0f : std::max(0.0f, center_phi - margin),
          invert ? PI_F : std::min(PI_F, center_phi + margin),
          cosf(radius),
          sinf(radius)};
}

/**
 * @brief Emit the single horizontal interval where a row crosses a great-circle
 * "cap" of half-angle `acos(cos_cap)` centred on an axis whose projection onto
 * the scan plane is (ny, r_val, alpha_angle). Shared by PlanarPolygon /
 * SphericalPolygon / Star, whose scanline math is otherwise identical.
 *
 * @tparam W Canvas width in columns.
 * @tparam OutputIt Sink type invoked as out(float start, float end).
 * @param cos_cap cos of the cap's angular radius (cosf(circumradius) etc.).
 * @param ny y-component of the cap axis.
 * @param r_val Horizontal projection length of the cap axis.
 * @param alpha_angle Azimuth of the cap axis (radians).
 * @param cos_phi Cosine of the row's polar angle.
 * @param sin_phi Sine of the row's polar angle.
 * @param out Sink accepting (float start, float end).
 * @return false to request a full-width fallback scan, true if the (possibly
 *         empty) interval was handled.
 */
template <int W, typename OutputIt>
inline bool emit_cap_interval(float cos_cap, float ny, float r_val,
                              float alpha_angle, float cos_phi, float sin_phi,
                              OutputIt out) {
  if (r_val < MIN_HORIZONTAL_PROJ)
    return false;

  float denom = r_val * sin_phi;
  if (std::abs(denom) < INTERVAL_DENOM_EPS)
    return false;

  float C_min = (cos_cap - ny * cos_phi) / denom;
  if (C_min > 1.0f)
    return true; // Completely outside the cap
  if (C_min < -1.0f)
    return false; // Full scan fallback

  // fast_acos: ~5e-5 rad peak error ≈ 0.002 px at W=288, far under the
  // floor/ceil pad below. Matches the Ring/DistortedRing scanline path.
  float d_alpha = fast_acos(C_min);
  float scale = W / TWO_PI_F;
  float x1 = floorf((alpha_angle - d_alpha) * scale);
  float x2 = ceilf((alpha_angle + d_alpha) * scale);
  out(x1, x2);
  return true;
}

/**
 * @brief Emit one scanline row of a cap-bounded leaf with a one-pixel cap pad.
 * @tparam W Canvas width in columns.
 * @tparam H Canvas height in rows.
 * @tparam OutputIt Sink type invoked as out(float start, float end).
 * @param sign +1 fills the shape, -1 its complement.
 * @param cos_cap Cosine of the bounding cap's angular radius.
 * @param sin_cap Sine of the bounding cap's angular radius.
 * @param ny y-component of the cap axis.
 * @param r_val Horizontal projection length of the cap axis.
 * @param alpha_angle Azimuth of the cap axis (radians).
 * @param y The row index.
 * @param out Sink accepting (float start, float end).
 * @return False to request a full-width fallback scan, true if the (possibly
 *         empty) interval was handled.
 * @details The complement wraps every row, so sign < 0 always requests the full
 * scan. Shared by PlanarPolygon / SphericalPolygon / Star.
 */
template <int W, int H, typename OutputIt>
inline bool emit_padded_cap_row(float sign, float cos_cap, float sin_cap,
                                float ny, float r_val, float alpha_angle, int y,
                                OutputIt out) {
  if (sign < 0.0f)
    return false;
  if (!TrigLUT<W, H>::initialized)
    TrigLUT<W, H>::init();
  // Column 1 of the theta LUT is one pixel of azimuth, so the cap pad is an
  // angle addition rather than a per-row cosf.
  const float cos_pad = TrigLUT<W, H>::cos_theta(1);
  // A cap padded past pi covers the sphere; cos turns back up there, so the
  // addition would report a cap tighter than the shape.
  const float cos_padded =
      cos_cap <= -cos_pad
          ? -1.0f
          : cos_cap * cos_pad - sin_cap * TrigLUT<W, H>::sin_theta[1];
  return emit_cap_interval<W>(cos_padded, ny, r_val, alpha_angle,
                              TrigLUT<W, H>::cos_phi[y],
                              TrigLUT<W, H>::sin_phi[y], out);
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

  float scale = W / TWO_PI_F;
  float safe_threshold = TWO_PI_F / W;

  if (angle_min <= safe_threshold) {
    out(floorf((alpha_angle - angle_max) * scale),
        ceilf((alpha_angle + angle_max) * scale));
  } else if (angle_max >= PI_F - safe_threshold) {
    out(floorf((alpha_angle + angle_min) * scale),
        ceilf((alpha_angle + TWO_PI_F - angle_min) * scale));
  } else {
    out(floorf((alpha_angle - angle_max) * scale),
        ceilf((alpha_angle - angle_min) * scale));
    out(floorf((alpha_angle + angle_min) * scale),
        ceilf((alpha_angle + angle_max) * scale));
  }
}

} // namespace SDF
