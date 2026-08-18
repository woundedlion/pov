/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include "math/geometry.h"
#include "engine/constants.h"
#include "engine/memory.h"
#include "engine/static_circular_buffer.h"
#include "render/sdf/common.h"

/**
 * @file csg.h
 * @brief The CSG operators that combine shapes: Union, SmoothUnion, Subtract,
 * Intersection and the AngularRepeat domain fold.
 */

namespace SDF {

/**
 * @brief Row bounds covering either child's band, padded and clamped.
 * @param a First child's band; culled when y_min > y_max.
 * @param b Second child's band, same convention.
 * @param pad Rows added on each side of the union.
 * @param max_y Highest row the result may name (inclusive).
 * @return Inclusive bounds over the padded union, clamped to [0, max_y], or the
 *         culled sentinel {1, 0} when both children are culled.
 * @details A culled child contributes nothing rather than dragging the union to
 * its sentinel rows.
 */
inline Bounds union_vertical_bounds(Bounds a, Bounds b, int pad, int max_y) {
  const bool a_culled = a.y_min > a.y_max;
  const bool b_culled = b.y_min > b.y_max;
  if (a_culled && b_culled)
    return {1, 0};
  const int lo = a_culled   ? b.y_min
                 : b_culled ? a.y_min
                            : std::min(a.y_min, b.y_min);
  const int hi = a_culled   ? b.y_max
                 : b_culled ? a.y_max
                            : std::max(a.y_max, b.y_max);
  return {std::max(0, lo - pad), std::min(max_y, hi + pad)};
}

/**
 * @brief Tests whether a row lies inside a band widened by @p pad rows.
 * @param b Band to test against; a culled band (y_min > y_max) holds no row.
 * @param y Row index.
 * @param pad Rows of slack allowed on each side.
 * @return True when y is within pad rows of the band.
 */
inline bool row_within_padded_band(Bounds b, int y, int pad) {
  return b.y_min <= b.y_max && y >= b.y_min - pad && y <= b.y_max + pad;
}

/**
 * @brief CSG Union operation (A + B), taking the minimum distance of two
 * shapes.
 * @tparam A First child shape type.
 * @tparam B Second child shape type.
 */
template <typename A, typename B> struct Union {
  const A &a; /**< First child shape. */
  const B &b; /**< Second child shape. */
  static constexpr bool is_solid =
      A::is_solid; /**< Both children share solidity, pinned by the
                        static_assert below. */

  static_assert(SDFShape<A> && SDFShape<B>,
                "CSG Union children must be SDF shapes (is_solid)");
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
  Union(const A &shape_a, const B &shape_b) : a(shape_a), b(shape_b) {}

  /**
   * @brief Deleted constructors from a temporary child.
   * @details Both children are retained by reference and read on every pixel,
   * so binding a temporary would dangle from the first distance() call.
   */
  Union(const A &&, const B &) = delete;
  Union(const A &, const B &&) = delete;
  Union(const A &&, const B &&) = delete;

  /**
   * @brief Row bounds spanning the union of the children's bands.
   * @tparam H Canvas height in rows.
   * @return Inclusive row bounds covering either child.
   */
  template <int H> Bounds get_vertical_bounds() const {
    auto b1 = a.template get_vertical_bounds<H>();
    auto b2 = b.template get_vertical_bounds<H>();
    return union_vertical_bounds(b1, b2, 0, H - 1);
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
    ScratchScope scratch(scratch_arena_b);
    MergedIntervalBuffer &merged = scratch_spans<MergedIntervalBuffer>(scratch);

    // One child fell back to full width: the whole row needs the full scan, so
    // the merged buffer is discarded and B need not be evaluated.
    bool has_a = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { push_interval(merged, start, end); });
    if (!has_a)
      return false;

    bool has_b = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { push_interval(merged, start, end); });
    if (!has_b)
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
  const A &a; /**< First child shape. */
  const B &b; /**< Second child shape. */
  float k;    /**< Smoothing radius in radians (e.g. 0.1). */
  static constexpr bool is_solid =
      A::is_solid; /**< Both children share solidity, pinned by the
                        static_assert below. */

  static_assert(SDFShape<A> && SDFShape<B>,
                "CSG SmoothUnion children must be SDF shapes (is_solid)");
  static_assert(A::is_solid == B::is_solid,
                "CSG SmoothUnion children must share solidity; a solid+stroke "
                "mix renders the stroke winner through the solid AA branch");
  static_assert(
      sdf_max_spans<A>::value + sdf_max_spans<B>::value <=
          2 * INTERVAL_SPAN_CAP,
      "nested CSG smooth-union exceeds MergedIntervalBuffer capacity; "
      "flatten the union or raise INTERVAL_SPAN_CAP");
  static_assert(blends_smoothly<A> && blends_smoothly<B>,
                "CSG SmoothUnion children must report a signed distance "
                "outside their surface; a far-sentinel clamper welds to "
                "nothing, making this identical to Union");

  /**
   * @brief Builds a smooth union of two child shapes.
   * @param shape_a First child shape.
   * @param shape_b Second child shape.
   * @param smoothness Blend radius k in radians.
   */
  SmoothUnion(const A &shape_a, const B &shape_b, float smoothness)
      : a(shape_a), b(shape_b), k(smoothness) {
    HS_CHECK(k > 0.0f);
  }

  /**
   * @brief Deleted constructors from a temporary child.
   * @details Both children are retained by reference and read on every pixel,
   * so binding a temporary would dangle from the first distance() call.
   */
  SmoothUnion(const A &&, const B &, float) = delete;
  SmoothUnion(const A &, const B &&, float) = delete;
  SmoothUnion(const A &&, const B &&, float) = delete;

  /**
   * @brief Row bounds spanning both children's bands, padded by the blend
   * radius.
   * @tparam H Canvas height in rows.
   * @return Inclusive row bounds expanded by k (converted to rows).
   */
  template <int H> Bounds get_vertical_bounds() const {
    auto b1 = a.template get_vertical_bounds<H>();
    auto b2 = b.template get_vertical_bounds<H>();
    return union_vertical_bounds(b1, b2, pad_rows<H>(), H - 1);
  }

  /**
   * @brief The blend radius as a row count.
   * @tparam H Canvas height in rows.
   * @return k (radians) converted to rows, at least 1.
   */
  template <int H> int pad_rows() const {
    constexpr int H_VIRT = H + hs::H_OFFSET;
    // phi spans [0,π] over (H_VIRT-1) rows.
    return std::max(1, static_cast<int>(ceilf(k * (H_VIRT - 1) / PI_F)));
  }

  /**
   * @brief Conservative union of the children's intervals, padded by k.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false (full scan) if either child
   *         falls back to full width, or if neither child covers a row still
   *         within the weld's reach.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    ScratchScope scratch(scratch_arena_b);
    MergedIntervalBuffer &merged = scratch_spans<MergedIntervalBuffer>(scratch);
    // Great-circle weld radius k spans k/sin(phi) columns of azimuth; the
    // equatorial conversion under-covers toward the poles. Clamp to full width
    // where the latitude factor diverges.
    float sin_phi = TrigLUT<W, H>::sin_phi[y];
    float pad_px =
        sin_phi > INTERVAL_DENOM_EPS
            ? std::min(k * W / TWO_PI_F / sin_phi, static_cast<float>(W))
            : static_cast<float>(W);

    // One child fell back to full width: the whole row needs the full scan, so
    // the merged buffer is discarded and B need not be evaluated. The weld
    // still blends both children through distance() on every scanned pixel.
    bool has_a = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          push_interval(merged, start - pad_px, end + pad_px);
        });
    if (!has_a)
      return false;

    bool has_b = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          push_interval(merged, start - pad_px, end + pad_px);
        });
    if (!has_b)
      return false;

    if (merged.is_empty()) {
      // Neither child covers this row, so there is no span to pad, but the weld
      // bulges outside both children's bands: a row within the blend reach must
      // be scanned in full or its fringe never renders. Rows beyond the reach of
      // either band hold no surface.
      const int pad = pad_rows<H>();
      return !row_within_padded_band(a.template get_vertical_bounds<H>(), y,
                                     pad) &&
             !row_within_padded_band(b.template get_vertical_bounds<H>(), y,
                                     pad);
    }

    // Emitted spans may straddle θ=0 and are not seam-normalized to [0,W): the
    // union merge is frame-tolerant, so scan_region's wrap+coalesce is the seam
    // authority (Subtract/Intersection normalize first because their pairwise
    // span comparison is frame-sensitive).
    merge_intervals(merged, out);
    return true;
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
  const A &a; /**< Minuend shape. */
  const B &b; /**< Subtrahend shape (removed from A). */
  static constexpr bool is_solid =
      A::is_solid; /**< Tracks the minuend; carving B never changes A's
                      solidity. */

  static_assert(SDFShape<A> && SDFShape<B>,
                "CSG Subtract children must be SDF shapes (is_solid)");
  // The minuend is collected into an IntervalBuffer (cap INTERVAL_SPAN_CAP)
  // before the seam split, so a minuend that could emit more spans must be
  // rejected at compile time rather than trapping in push_interval at runtime.
  static_assert(sdf_max_spans<A>::value <= INTERVAL_SPAN_CAP,
                "nested CSG Subtract minuend exceeds IntervalBuffer capacity; "
                "flatten the nesting or raise INTERVAL_SPAN_CAP");

  /**
   * @brief Builds a subtraction (shape_a minus shape_b).
   * @param shape_a Minuend shape.
   * @param shape_b Subtrahend shape.
   */
  Subtract(const A &shape_a, const B &shape_b) : a(shape_a), b(shape_b) {}

  /**
   * @brief Deleted constructors from a temporary child.
   * @details Both children are retained by reference and read on every pixel,
   * so binding a temporary would dangle from the first distance() call.
   */
  Subtract(const A &&, const B &) = delete;
  Subtract(const A &, const B &&) = delete;
  Subtract(const A &&, const B &&) = delete;

  /**
   * @brief Row bounds of the minuend (subtraction never grows the band).
   * @tparam H Canvas height in rows.
   * @return The minuend's inclusive row bounds.
   */
  template <int H> Bounds get_vertical_bounds() const {
    return a.template get_vertical_bounds<H>();
  }

  /**
   * @brief Emits the minuend's intervals, seam-split into [0, W).
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan when A
   *         falls back.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    ScratchScope scratch(scratch_arena_b);
    IntervalBuffer &intervals_a = scratch_spans<IntervalBuffer>(scratch);

    bool has_a = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          push_interval(intervals_a, start, end);
        });

    if (!has_a)
      return false;

    if (intervals_a.is_empty())
      return true;

    // A child's spans bound its coverage rather than matching it: a polygon or
    // star emits its circumscribed cap, and every cap row is padded a pixel
    // wide for AA. Differencing those bounds would cull columns B does not
    // cover — a star's notches, the carve-edge AA fringe — and a stroke B's
    // edge-bands can coalesce into one chord spanning its hollow interior. So
    // emit A's spans (seam-split into [0, W)) and let per-pixel max(A, -B)
    // carve exactly B, with no horizontal culling: Subtract pays full
    // A-coverage shading.
    constexpr size_t SEAM_SPLIT_CAP = 2 * INTERVAL_SPAN_CAP;
    static_assert(2 * sdf_max_spans<A>::value <= SEAM_SPLIT_CAP,
                  "post-seam-split span count exceeds norm buffer capacity");
    using NormBuffer = StaticCircularBuffer<Interval, SEAM_SPLIT_CAP>;
    NormBuffer &norm_a = scratch_spans<NormBuffer>(scratch);
    normalize_intervals_to_range<W>(intervals_a, norm_a);
    for (size_t i = 0; i < norm_a.size(); ++i)
      out(norm_a[i].first, norm_a[i].second);
    return true;
  }

  /**
   * @brief Signed distance to the subtraction, writing into res.
   * @tparam ComputeUVs Forwarded to each child's distance().
   * @param p Point on sphere (normalized).
   * @param res Output result; max(A, -B) with B's distance negated when it
   * wins, keeping A's size metric.
   * @note Across the carve edge only `dist` is negated and `size` is held to
   *       A's; `t`/`raw_dist`/`aux` stay B's own parametrization, so a shader
   *       keying off them sees a hard edge there. `raw_dist` is a per-shape
   *       unsigned or supplementary quantity, not a signed metric, so it is
   *       not negated with `dist`.
   * @note The carve edge is anti-aliased on one side only when B clamps to
   *       FAR_SENTINEL past its reject band (blends_smoothly == false). The
   *       sentinel loses the max, so at the band edge the composite jumps from
   *       B's ramp to A's own distance instead of completing the outer half of
   *       the fringe. Ring's band is exactly its stroke, putting that step on
   *       the carve edge itself; the other clampers carry a margin of true
   *       distance past their surface, which absorbs the step wherever that
   *       margin outruns the AA reach.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    a.template distance<ComputeUVs>(p, res);
    DistanceResult res_b;
    b.template distance<ComputeUVs>(p, res_b);
    if (-res_b.dist > res.dist) {
      // is_solid tracks the minuend, so the AA branch normalizes over A's
      // metric; the carve edge must not widen to B's.
      const float size_a = res.size;
      res_b.dist = -res_b.dist;
      res = res_b;
      res.size = size_a;
    }
  }
};

/**
 * @brief CSG Intersection operation (A & B), taking the maximum distance.
 * @tparam A First child shape type.
 * @tparam B Second child shape type.
 */
template <typename A, typename B> struct Intersection {
  const A &a; /**< First child shape. */
  const B &b; /**< Second child shape. */
  static constexpr bool is_solid =
      A::is_solid; /**< Both children share solidity, pinned by the
                        static_assert below. */

  static_assert(SDFShape<A> && SDFShape<B>,
                "CSG Intersection children must be SDF shapes (is_solid)");
  static_assert(A::is_solid == B::is_solid,
                "CSG Intersection children must share solidity; a solid+stroke "
                "mix renders the solid winner through the stroke AA branch");
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
  Intersection(const A &shape_a, const B &shape_b) : a(shape_a), b(shape_b) {}

  /**
   * @brief Deleted constructors from a temporary child.
   * @details Both children are retained by reference and read on every pixel,
   * so binding a temporary would dangle from the first distance() call.
   */
  Intersection(const A &&, const B &) = delete;
  Intersection(const A &, const B &&) = delete;
  Intersection(const A &&, const B &&) = delete;

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
    ScratchScope scratch(scratch_arena_b);
    IntervalBuffer &intervals_a = scratch_spans<IntervalBuffer>(scratch);
    IntervalBuffer &intervals_b = scratch_spans<IntervalBuffer>(scratch);

    bool has_a = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          push_interval(intervals_a, start, end);
        });

    bool has_b = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          push_interval(intervals_b, start, end);
        });

    // Both children fell back: the row is a full scan, and the protocol forbids
    // emitting anything alongside it.
    if (!has_a && !has_b)
      return false;

    // A fallback child bounds nothing, so the other child's spans still bound
    // the intersection; replay the buffer.
    if (!has_a) {
      for (size_t i = 0; i < intervals_b.size(); ++i)
        out(intervals_b[i].first, intervals_b[i].second);
      return true;
    }
    if (!has_b) {
      for (size_t i = 0; i < intervals_a.size(); ++i)
        out(intervals_a[i].first, intervals_a[i].second);
      return true;
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
    using NormBuffer = StaticCircularBuffer<Interval, SEAM_SPLIT_CAP>;
    NormBuffer &norm_a = scratch_spans<NormBuffer>(scratch);
    NormBuffer &norm_b = scratch_spans<NormBuffer>(scratch);
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
  float sector;    /**< Angular width of one sector, 2*PI / repetitions. */
  float reciprocal_sector; /**< 1 / sector, for the per-pixel fold. */
  static constexpr bool is_solid =
      Shape::is_solid; /**< Matches the child's solidity. */

  static_assert(SDFShape<Shape>,
                "AngularRepeat child must be an SDF shape (is_solid)");

  /**
   * @brief Repeats the shape around an arbitrary axis.
   * @param s Child shape to repeat.
   * @param reps Number of copies (must be > 0).
   * @param ax Rotation axis (unit length).
   */
  AngularRepeat(const Shape &s, int reps, const Vector &ax)
      : shape(s), axis(ax), repetitions(reps),
        sector(TWO_PI_F / static_cast<float>(reps)),
        reciprocal_sector(static_cast<float>(reps) / TWO_PI_F) {
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
   * @brief Deleted constructors from a temporary child.
   * @details The child is retained by reference and read on every pixel, so
   * binding a temporary would dangle from the first distance() call. The axis
   * is copied, so a temporary Vector stays legal.
   */
  AngularRepeat(const Shape &&, int, const Vector &) = delete;
  AngularRepeat(const Shape &&, int) = delete;

  /**
   * @brief Row bounds for the repeated shape.
   * @tparam H Canvas height in rows.
   * @return The child's band for a Y-axis fold, or when the child culls; the
   *         full canvas otherwise.
   */
  template <int H> Bounds get_vertical_bounds() const {
    Bounds child = shape.template get_vertical_bounds<H>();
    if (child.y_min > child.y_max)
      return child;
    // Only a Y-axis fold (axis.y near ±1) preserves latitude; any other axis
    // sweeps latitudes the child never occupies, so its band cannot bound the
    // copies and every row must be scanned.
    if (fabsf(axis.y) < 1.0f - TOLERANCE)
      return {0, H - 1};
    return child;
  }

  /**
   * @brief Replays the child's spans once per copy for a Y-axis fold.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True when the emitted spans describe the row; false requests a
   *         full scan.
   * @details A Y-axis fold shifts a pixel's azimuth by a whole sector and holds
   * its latitude, so every column a copy covers is a child column shifted by a
   * multiple of W / repetitions. Replaying every child span at every shift is a
   * superset of the covered columns — a child span outside the sector the fold
   * actually emits has no pre-image — which is the direction the cull contract
   * allows. Any other axis rotates the copies off the row, so the child's spans
   * cannot bound them and the row falls back to a full scan.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    if (fabsf(axis.y) < 1.0f - TOLERANCE)
      return false;

    ScratchScope scratch(scratch_arena_b);
    IntervalBuffer &child = scratch_spans<IntervalBuffer>(scratch);
    bool overflow = false;
    bool handled = shape.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          if (child.is_full())
            overflow = true;
          else
            child.push_back({start, end});
        });
    // Nothing is emitted before every fallback below: the protocol pairs a
    // false return with an empty emission.
    if (!handled || overflow)
      return false;
    if (child.is_empty())
      return true;
    if (child.size() * static_cast<size_t>(repetitions) >
        ANGULAR_REPEAT_SPAN_CAP)
      return false;

    const float step = static_cast<float>(W) / static_cast<float>(repetitions);
    const float pad =
        1.0f + ANGULAR_REPEAT_FOLD_SLOP * static_cast<float>(W) / TWO_PI_F;
    for (size_t i = 0; i < child.size(); ++i)
      // Copies of a span this wide abut, covering every column anyway, and a
      // padded span at least a sector long would break the len <= W contract.
      if (child[i].second - child[i].first + 2.0f * pad >= step)
        return false;

    for (size_t i = 0; i < child.size(); ++i)
      for (int k = 0; k < repetitions; ++k) {
        const float shift = static_cast<float>(k) * step;
        out(child[i].first + shift - pad, child[i].second + shift + pad);
      }
    return true;
  }

  /**
   * @brief Folds p into one sector and evaluates the child, writing into res.
   * @tparam ComputeUVs Forwarded to the child's distance().
   * @param p Point on sphere (normalized).
   * @param res Output result of the child at the folded point. Note res's UV
   *        registers (t / azimuth) are sector-local: the child sees the folded
   *        point, so t spans one sector and resets at each sector boundary.
   * @note The folded distance is exact only while the child stays inside one
   *       sector. A child crossing a sector boundary reports the distance to
   *       the folded copy rather than to the nearest copy, over-estimating near
   *       the boundary, so AA bands and any distance-driven cull must keep the
   *       shape within its sector.
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
      theta += TWO_PI_F;

    float folded_theta =
        centered_sector_angle(theta, sector, reciprocal_sector);

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

} // namespace SDF
