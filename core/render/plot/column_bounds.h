/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/plot/column_geometry.h"

/** @file column_bounds.h @brief Swept edge culling and column walk frame bounds. */

namespace Plot {

HS_O3_BEGIN

/** @brief Rejects edges whose swept arcs miss a draw's rows and columns.
 * @tparam W Display width.
 * @details Rows use the column walk's own sample cut; columns use the render
 * band widened by two columns on each side. */
template <int W> struct SweptEdgeCull {
  float top;           /**< Highest unit-sphere y the walk emits. */
  float bottom;        /**< Lowest unit-sphere y the walk emits. */
  math::Vector first;  /**< First longitude half-plane normal. */
  math::Vector second; /**< Second longitude half-plane normal. */
  int mode = 0; /**< 0: rows only; 1: hidden outside either half-plane; 2:
                   hidden inside both. */

  /** @brief Prepares the bounds for one draw.
   * @param columns Render-band column predicate.
   * @param top_y,bottom_y Unit-sphere y limits of emitted samples. */
  HS_HOT_FLASH_MEMBER SweptEdgeCull(const ClipRegion::XClip &columns,
                                    float top_y, float bottom_y)
      : top(top_y), bottom(bottom_y) {
    if (!columns.active)
      return;
    constexpr float MARGIN = 2.0f;
    const float start =
        (float(columns.rs) - 0.5f - MARGIN) * math::RADIANS_PER_COLUMN<W>;
    const float width = (float(columns.length(W)) + 2.0f * MARGIN) *
                        math::RADIANS_PER_COLUMN<W>;
    if (width >= 2.0f * math::PI_F)
      return;
    const float end = start + width;
    if (width <= math::PI_F) {
      mode = 1;
      first = {-sinf(start), 0.0f, cosf(start)};
      second = {sinf(end), 0.0f, -cosf(end)};
    } else {
      mode = 2;
      first = {-sinf(end), 0.0f, cosf(end)};
      second = {sinf(start), 0.0f, -cosf(start)};
    }
  }

  /** @brief Highest unit-sphere y a canvas's column walk emits.
   * @tparam H Display height.
   * @param canvas Frame whose clip rows bound the walk.
   * @return The height one row above the clip's render margin. */
  template <int H> HS_HOT_INLINE static float top_of(const Canvas &canvas) {
    return cosf(
        fmaxf(0.0f, math::DisplayGeometry<H>::row_to_phi(float(
                        canvas.clip().y_start - canvas.clip().margin - 1))));
  }

  /** @brief Highest direction that can light a physical row without tap motion.
   * @param row_phi Polar angle of the row immediately above the display.
   * @return Unit-sphere height including projection error. */
  HS_HOT_FLASH_MEMBER static float physical_top(float row_phi) {
    return cosf(fmaxf(
        0.0f, row_phi - 2.0f * static_cast<float>(math::FAST_ACOS_MAX_ERROR)));
  }

  /** @brief Lowest unit-sphere y a canvas's column walk emits.
   * @tparam H Display height.
   * @param canvas Frame whose clip rows bound the walk.
   * @return The height one row below the clip's render margin. */
  template <int H> HS_HOT_INLINE static float bottom_of(const Canvas &canvas) {
    return cosf(std::clamp(
        math::DisplayGeometry<H>::row_to_phi(
            float(std::min(H, canvas.clip().y_end + canvas.clip().margin + 1))),
        0.0f, math::PI_F));
  }

  /** @brief The bounds of a canvas's column walk.
   * @tparam H Display height.
   * @param canvas Frame whose clip rows and render band bound the walk.
   * @return The cull. */
  template <int H> static SweptEdgeCull of(const Canvas &canvas) {
    return SweptEdgeCull(canvas.render_x_clip(), top_of<H>(canvas),
                         bottom_of<H>(canvas));
  }

  /** @brief Tests whether every geodesic between points a frame sweeps
   * stays outside the bounds.
   * @param points Positions whose normalized convex hull, padded, holds every
   * drawn sample.
   * @param count Number of points.
   * @param pad Bound on a drawn sample's distance from that hull.
   * @return True when no drawn sample can be emitted. */
  HS_HOT_INLINE bool rejects(const math::Vector *points, int count,
                             float pad) const {
    if (count == 2)
      return rejects_hull(points, 2, pad);
    return rejects_variable(points, count, pad);
  }

  /** @brief Tests a trail chain using its midpoint's enclosing chord bound.
   * @param world Unit positions in chain order.
   * @param count Number of positions.
   * @param pad Bound beyond the chain's captured positions.
   * @return True when the padded chain misses the clip. */
  HS_HOT_INLINE bool rejects_chain(const math::Vector *world, int count,
                                   float pad) const {
    const math::Vector &centre = world[count / 2];
    float reach_cosine = 1.0f;
    for (int j = 0; j < count; ++j)
      reach_cosine = fminf(reach_cosine, math::dot(centre, world[j]));
    return rejects(&centre, 1,
                   pad + sqrtf(fmaxf(0.0f, 2.0f - 2.0f * reach_cosine)));
  }

  /** @brief Tests whether every arc an edge sweeps over the frame stays
   * outside the bounds.
   * @param points Both endpoints at the frame's start and end.
   * @param pad Bound on a drawn endpoint's distance from its chord.
   * @return True when no drawn sample of the edge can be emitted. */
  HS_HOT_FLASH_MEMBER bool rejects(const math::Vector (&points)[4],
                                   float pad) const {
    return rejects_hull(points, 4, pad);
  }

private:
  /** @brief Tests a variable-size padded hull against the swept clip.
   * @param points Unit captures.
   * @param count Number of captures.
   * @param pad Distance from the normalized hull.
   * @return Whether the padded hull misses the clip. */
  HS_HOT_FLASH_MEMBER bool rejects_variable(const math::Vector *points,
                                            int count, float pad) const {
    return rejects_hull(points, count, pad);
  }

  /** @brief rejects() for any count, inlined so a fixed count unrolls.
   * @param points Positions whose padded hull holds every drawn sample.
   * @param count Number of points.
   * @param pad Bound on a drawn sample's distance from that hull.
   * @return True when no drawn sample can be emitted. */
  HS_HOT_INLINE bool rejects_hull(const math::Vector *points, int count,
                                  float pad) const {
    if (mode == 1 && (Spacetime::swept_beyond(points, count, pad, -first) ||
                      Spacetime::swept_beyond(points, count, pad, -second)))
      return true;
    if (mode == 2 && Spacetime::swept_beyond(points, count, pad, first) &&
        Spacetime::swept_beyond(points, count, pad, second))
      return true;
    float highest = points[0].y, lowest = points[0].y;
    for (int i = 0; i < count; ++i) {
      if (points[i].y >= bottom && points[i].y <= top)
        return false;
      highest = fmaxf(highest, points[i].y);
      lowest = fminf(lowest, points[i].y);
    }
    highest += pad;
    lowest -= pad;
    if ((highest <= 0.0f && highest < bottom) ||
        (lowest >= 0.0f && lowest > top)) {
      float least = 1.0f;
      for (int i = 0; i < count; ++i)
        for (int j = i + 1; j < count; ++j)
          least = fminf(least, math::dot(points[i], points[j]));
      return least > 2.0f * pad;
    }
    float low, high;
    return Spacetime::swept_heights(points, count, pad, low, high) &&
           (high < bottom || low > top);
  }
};

/** @brief A draw's column-walk constants: column window trigonometry and the
 *         canvas clip's row and column bounds.
 * @tparam W,H Display dimensions. */
template <int W, int H> struct ColumnWalkFrame {
  float half_sine;       /**< Sine of half a column. */
  float half_cosine;     /**< Cosine of half a column. */
  float full_sine;       /**< Sine of one column. */
  float full_cosine;     /**< Cosine of one column. */
  float top_y;           /**< Highest unit-sphere y the walk emits. */
  float bottom_y;        /**< Lowest unit-sphere y the walk emits. */
  float physical_top_y;  /**< Highest y lighting a physical row without tap
                             motion. */
  SweptEdgeCull<W> cull; /**< The clip's swept-edge bounds. */

  /** @brief Computes the constants of draws into a canvas.
   * @param canvas Frame whose clip bounds the walk.
   * @return The constants. */
  HS_HOT_FLASH_MEMBER static ColumnWalkFrame of(const Canvas &canvas) {
    const float top = SweptEdgeCull<W>::template top_of<H>(canvas);
    const float bottom = SweptEdgeCull<W>::template bottom_of<H>(canvas);
    return {sinf(0.5f * math::RADIANS_PER_COLUMN<W>),
            cosf(0.5f * math::RADIANS_PER_COLUMN<W>),
            sinf(math::RADIANS_PER_COLUMN<W>),
            cosf(math::RADIANS_PER_COLUMN<W>),
            top,
            bottom,
            SweptEdgeCull<W>::physical_top(
                math::DisplayGeometry<H>::row_to_phi(-1.0f)),
            SweptEdgeCull<W>(canvas.render_x_clip(), top, bottom)};
  }
};

namespace detail {
/** @brief One arm pass's walk over the columns an edge can light. */
struct ColumnSegment {
  int start = 0;    /**< First column in walk order. */
  int count = 0;    /**< Columns to walk. */
  int deadline = 0; /**< Columns walked unlit before the walk gives up. */
};

/** @brief A column's position relative to a modelled longitude arc. */
struct ColumnArcPosition {
  /** Wrapped distance from the arc's start in its direction, columns. */
  float offset;
  /** Absolute endpoint separation, columns. */
  float extent;
};

} // namespace detail

HS_O3_END

} // namespace Plot
