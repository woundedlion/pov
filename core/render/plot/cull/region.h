/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <utility>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <array>
#include "core/math/geometry.h"
#include "core/render/shading.h"
#include "core/render/clip.h"
#include "core/render/filter/splat.h"
#include "core/memory.h"

/**
 * @file region.h
 * @brief Cull tolerances, the arc-length helper, and CullRegion, the
 * motion-padded clip the edge samplers and span tests cull against.
 */

namespace Plot {

/** @brief Arc length between unit sphere points, preserving short chords.
 *  @param a First unit vector.
 *  @param b Second unit vector.
 *  @return Angle between @p a and @p b, in radians. */
__attribute__((always_inline)) inline float
unit_arc_length(const math::Vector &a, const math::Vector &b) {
  const math::Vector chord = a - b;
  const float length_sq = math::dot(chord, chord);
  // Below 1e-3 radians the chord equals the arc to float precision.
  return length_sq < 1e-6f ? sqrtf(length_sq) : math::angle_between(a, b);
}

/**
 * @brief Geodesic segment shorter than this (radians) collapses to a point.
 * @details A slerp-axis stability bound that picks the interpolation strategy,
 * independent of math::EPS_GEOMETRIC.
 */
inline constexpr float EPS_GEODESIC_SEGMENT = 0.001f;

/**
 * @brief Minimum |cross(a, b)|² for which the arc pole of a geodesic edge is
 *        taken from the cross product rather than a stable perpendicular.
 * @details |cross| = sin(angle), so 1e-8 is the angular 1e-4 band without
 * angle_between's ULP amplification near a normalized dot of ±1.
 */
inline constexpr float EPS_ARC_POLE_SQ = 1e-8f;

/**
 * @brief Minimum |axis.y| for which a geodesic edge's endpoint columns bound
 *        its azimuth span.
 * @details Below this the great circle runs near the poles, where longitude is
 * ill-conditioned and the interior leaves the endpoint columns.
 */
inline constexpr float AXIS_Y_EPS = 1e-4f;

/**
 * @brief Minimum worst-case sin(φ) over a curve for which its plotted columns
 *        are trusted enough to cull by.
 * @details The azimuth Lipschitz bound scales as 1/sin(φ), and nearer the poles
 * the plotted column is float noise; below this the column cull bails.
 */
inline constexpr float MIN_SIN_PHI = 0.05f;

/**
 * @brief Floor on the adaptive sub-step length, as a fraction of base_step.
 * @details Lower clamp on the screen-velocity step, which tends to zero where
 * the azimuthal velocity diverges at the poles.
 */
inline constexpr float MIN_POLE_SCALE = 0.05f;

/**
 * @brief Target screen-space spacing (pixels) between adaptive sub-samples.
 * @details Slightly sub-pixel so neighbouring AntiAlias splats overlap without
 * holes.
 */
inline constexpr float SCREEN_STEP_PX = 0.9f;

/**
 * @brief Columns of slack added on each side of a culled column span.
 * @details Absorbs plot rounding and the AntiAlias tap spread.
 */
inline constexpr int COL_PAD = 2;

/**
 * @brief Columns a padded span reaches past its fractional end.
 * @details The pad plus the boundary column ceil() adds.
 */
inline constexpr int COL_FOOTPRINT = COL_PAD + 1;

/**
 * @brief Rows of slack added to the high end of a geodesic row span.
 * @details AntiAlias emits a sample into floor(row) and floor(row) + 1.
 */
inline constexpr float GEODESIC_ROW_AA_PAD = 1.0f;

/**
 * @brief Columns outside the render band a clip cut is placed at.
 * @details A piece ending on the band edge still overlaps it after the
 * COL_FOOTPRINT widening; the extra column absorbs fast-trig error in the cut.
 */
inline constexpr int CLIP_CUT_COL_PAD = COL_FOOTPRINT + 1;

/**
 * @brief Rows outside the render band a clip cut is placed at.
 * @details could_intersect_y takes the band edge itself as an intersection, and
 * the row map's fast-trig round trip moves a cut by a fraction of a row.
 */
inline constexpr int CLIP_CUT_ROW_PAD = 1;

/** @brief A column interval, unwrapped past the width. */
struct ColumnSpan {
  float lo; ///< First column.
  float hi; ///< Last column.

  /** @return The interval's centre column. */
  HS_HOT_INLINE float centre() const { return 0.5f * (lo + hi); }
  /** @return The interval's half-width in columns. */
  HS_HOT_INLINE float half() const { return 0.5f * (hi - lo); }
};

/** @brief Centre and half-width of a padded longitude band. */
struct CentreHalf {
  float centre; ///< Band centre longitude, radians, unwrapped past a turn.
  float half;   ///< Band half-width, radians.
};

/** @brief How a longitude band's two bounding half-planes combine. */
enum class MeridianTest : uint8_t {
  NONE,   ///< The band spans the turn; every longitude passes.
  BOTH,   ///< The band spans at most half a turn; it lies in both.
  EITHER, ///< The band spans more than half a turn; it lies in either.
};

/** @brief A longitude band's bounding meridian half-planes. */
struct HalfPlanes {
  float start; ///< Band start longitude, radians.
  float end;   ///< Band end longitude, radians.
  float width; ///< Band width, radians.

  /** @return How the half-planes combine: none for a band of a full turn. */
  HS_HOT_INLINE MeridianTest test() const {
    if (width >= 2.0f * math::PI_F)
      return MeridianTest::NONE;
    return width <= math::PI_F ? MeridianTest::BOTH : MeridianTest::EITHER;
  }

  /** @return Normal of the start meridian's half-plane holding the band. */
  HS_HOT_INLINE math::Vector first() const {
    return {-sinf(start), 0.0f, cosf(start)};
  }
  /** @return Normal of the end meridian's half-plane holding the band. */
  HS_HOT_INLINE math::Vector second() const {
    return {sinf(end), 0.0f, -cosf(end)};
  }
};

/**
 * @brief A view of a clip's render columns [first(), end<W>()), crossing the
 *        seam when the end passes the width.
 * @details Every site uses first() + (count(w) - 1) / 2 as its centre.
 */
struct ColumnBand {
  const ClipRegion::XClip &columns; ///< Render-band column predicate.

  /** @brief The render band of a column predicate.
   *  @param xc Render-band column predicate; must outlive the band.
   *  @return The band. */
  HS_HOT_INLINE static ColumnBand of(const ClipRegion::XClip &xc) {
    return {xc};
  }

  /** @return The first render column, in [0, w). */
  HS_HOT_INLINE int first() const { return columns.rs; }

  /** @brief Columns in the band.
   *  @param w Canvas width in columns.
   *  @return The band length, counting past the seam; @p w when inactive. */
  HS_HOT_INLINE int count(int w) const { return columns.length(w); }

  /** @brief The exclusive end column, unwrapped past the width.
   *  @tparam W Canvas width in columns.
   *  @return The wrapped end, plus W when the band crosses the seam. */
  template <int W> HS_HOT_INLINE int end() const {
    return columns.wrap ? columns.re + W : columns.re;
  }

  /** @brief The columns whose samples' taps reach the band.
   *  @tparam W Canvas width in columns.
   *  @return (first() - 1, end<W>()). */
  template <int W> HS_HOT_INLINE ColumnSpan tap_reach() const {
    return {static_cast<float>(first() - 1), static_cast<float>(end<W>())};
  }

  /** @brief The columns a cut is placed at, @p pad outside the band.
   *  @tparam W Canvas width in columns.
   *  @param pad Columns between the band and each cut.
   *  @return The low and high cut columns, wrapped into [0, W). */
  template <int W> HS_HOT_INLINE std::array<int, 2> cut_columns(int pad) const {
    return {math::wrap(first() - pad, W), math::wrap(columns.re + pad, W)};
  }

  /** @brief Centre shared by the column culls.
   *  @param w Canvas width in columns.
   *  @return The midpoint of the first and last render columns. */
  HS_HOT_INLINE float centre(int w) const {
    return first() + static_cast<float>(count(w)) * 0.5f - 0.5f;
  }

  /** @param w Canvas width. @return Point-route half-width in columns. */
  HS_HOT_INLINE float point_half(int w) const { return count(w) * 0.5f + 0.5f; }

  /** @param w Canvas width. @return Swept-cell half-width in columns. */
  HS_HOT_INLINE float swept_half(int w) const { return count(w) * 0.5f + 2.0f; }

  /** @param w Canvas width. @return Cap half-width in columns. */
  HS_HOT_INLINE float cap_half(int w) const { return count(w) * 0.5f + 1.5f; }

  /** @param w Canvas width. @return Clip half-width including chart padding. */
  HS_HOT_INLINE float clip_half(int w) const { return count(w) * 0.5f + 2.5f; }

  /** @brief Columns at a site's half-width, unwrapped past the seam.
   *  @param w Canvas width. @param half Half-width in columns.
   *  @return Symmetric interval about centre(w). */
  HS_HOT_INLINE ColumnSpan span(int w, float half) const {
    const float mid = centre(w);
    return {mid - half, mid + half};
  }

  /** @brief Rounding allowance for the swept-cell endpoint evaluation. */
  static constexpr float SWEPT_ROUNDING = 0x1p-19f;

  /** @brief Bounding meridians at a site's half-width, with angular slack.
   *  @tparam W Canvas width in columns.
   *  @param half Half-width in columns.
   *  @param slack Radians added on each side.
   *  @return The widened span's half-planes. */
  template <int W>
  HS_HOT_INLINE HalfPlanes half_planes(float half, float slack = 0.0f) const {
    const ColumnSpan bounds = span(W, half);
    const float start = bounds.lo * math::RADIANS_PER_COLUMN<W> - slack;
    const float end = bounds.hi * math::RADIANS_PER_COLUMN<W> + slack;
    return {start, end, end - start};
  }

  /** @brief Longitude centre and half-width at a site's half-width.
   *  @param w Canvas width in columns.
   *  @param half Half-width in columns.
   *  @return The span's centre and half-width, radians. */
  HS_HOT_INLINE CentreHalf centre_half(int w, float half) const {
    return {centre(w) * math::TWO_PI_F / w, half * math::TWO_PI_F / w};
  }
};

/**
 * @brief The clip a draw culls against: the exact clip widened by how far
 *        the draw's samples stray from their end positions over the frame.
 * @details Rows widen by the bound in rows. A column span widens by the
 * bound over the span's worst sin(phi) less the bound; a span whose
 * remainder falls below MIN_SIN_PHI, or whose widening passes a quarter of
 * the width, culls by rows only. The exact clip is unchanged for the
 * splat taps and for scans sharing the canvas.
 */
class CullRegion {
public:
  /**
   * @brief A region around @p clip.
   * @param clip Exact clip; must outlive the region.
   * @param xc Column predicate of @p clip.
   * @param bound Radians any sample of the draw strays from its end
   *        position; 0 for a still draw.
   */
  CullRegion(const ClipRegion &clip, const ClipRegion::XClip &xc,
             float bound = 0.0f)
      : exact(clip), columns(xc), motion(bound) {}

  /**
   * @brief A region around @p clip, deriving its column predicate.
   * @param clip Exact clip; must outlive the region.
   * @param bound Radians any sample of the draw strays from its end
   *        position.
   */
  explicit CullRegion(const ClipRegion &clip, float bound = 0.0f)
      : CullRegion(clip, clip.x_clip(), bound) {}

  /** @return The exact clip. */
  const ClipRegion &clip() const { return exact; }
  /** @return The exact clip's column predicate. */
  const ClipRegion::XClip &x_clip() const { return columns; }
  /** @return The motion bound, radians. */
  float bound() const { return motion; }

  /**
   * @brief Whether rows [@p lo, @p hi] widened by the bound reach the clip.
   * @tparam H Canvas height in rows.
   * @param lo First row.
   * @param hi Last row.
   * @return The clip's could_intersect_y over the widened rows.
   */
  template <int H> bool could_intersect_rows(float lo, float hi) const {
    const float pad = motion * math::ROWS_PER_RADIAN<H>;
    return exact.could_intersect_y(lo - pad, hi + pad);
  }

  /**
   * @brief Columns a span widens by.
   * @tparam W Canvas width in columns.
   * @param sin_worst Smallest sin(phi) over the span's end positions.
   * @return The widening, or a negative value when the span must cull by
   *         rows only.
   */
  template <int W> float column_pad(float sin_worst) const {
    const float sin_moved = sin_worst - motion;
    if (sin_moved < MIN_SIN_PHI)
      return -1.0f;
    const float pad =
        motion * (static_cast<float>(W) / (2.0f * math::PI_F)) / sin_moved;
    return pad > W * 0.25f ? -1.0f : pad;
  }

  /**
   * @brief Columns a span between two rows widens by.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @param row_lo First row of the span.
   * @param row_hi Last row of the span.
   * @return column_pad() of the rows' smaller sin(phi).
   */
  template <int W, int H>
  float column_pad_rows(float row_lo, float row_hi) const {
    auto sin_at = [](float row) HS_HOT_INLINE {
      const float phi = math::DisplayGeometry<H>::row_to_phi(row);
      float sine, cosine;
      math::poly_sincosf_0_pi(fminf(math::PI_F, fmaxf(0.0f, phi)), sine,
                              cosine);
      return sine;
    };
    return column_pad<W>(fminf(sin_at(row_lo), sin_at(row_hi)));
  }

  /**
   * @brief Whether a padded column arc reaches the clip's columns.
   * @tparam W Canvas width in columns.
   * @param col_s Arc start column, in [0, W).
   * @param col_len Arc length in columns.
   * @param pad Widening from column_pad(); negative keeps the arc visible.
   * @return True when the widened arc overlaps the clip's columns.
   */
  template <int W>
  bool arc_reaches_columns(int col_s, int col_len, float pad) const {
    if (pad < 0.0f)
      return true;
    const int widen = static_cast<int>(ceilf(pad));
    const int len = col_len + 2 * widen;
    if (len >= W)
      return true;
    const int start = ((col_s - widen) % W + W) % W;
    return ClipRegion::arcs_overlap(columns.rs, columns.length(W), start, len,
                                    W);
  }

private:
  const ClipRegion &exact;   /**< Unmodified splat clip. */
  ClipRegion::XClip columns; /**< Unmodified column predicate. */
  float motion;              /**< Angular padding in radians. */
};

} // namespace Plot
