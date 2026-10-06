/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/constants.h" // MAX_W/MAX_H
#include "platform/platform.h"  // HS_AUDIT_CHECK

/**
 * @file clip.h
 * @brief ClipRegion: the segment clip rectangle plus its cylindrical
 * render-band topology (wrap, containment, arc overlap).
 */

/**
 * @brief Clip region for segment-based rendering.
 * @details Display bounds define the ISR's pixel range (exact segment).
 *          Render bounds expand by `margin` for filters that spread taps.
 *          `w`/`h` are the active canvas size.
 */
struct ClipRegion {
  int y_start = 0;   /**< Display top row (inclusive), in pixels. */
  int y_end = MAX_H; /**< Display bottom row (exclusive), in pixels. */
  int x_start = 0;   /**< Display left column (inclusive), in pixels. */
  int x_end = MAX_W; /**< Display right column (exclusive), in pixels. */
  int margin =
      1; /**< Render-bound expansion past the display edges, in pixels. */
  int w = MAX_W; /**< Canvas width, in pixels. */
  int h = MAX_H; /**< Canvas height, in pixels. */

  /** @brief Field-wise equality, for cached clip-state stamps. */
  bool operator==(const ClipRegion &) const = default;

  /**
   * @brief Render-region top edge: display top expanded up by `margin`, floored at 0.
   * @return First render row (inclusive), in pixels.
   * @details Only the low side is clamped; `y_start <= h` by the
   *          display-bounds invariant.
   */
  int render_y_start() const {
    return y_start - margin > 0 ? y_start - margin : 0;
  }
  /**
   * @brief Render-region bottom edge: display bottom expanded down by `margin`, capped at h.
   * @return One-past-last render row (exclusive), in pixels.
   * @details Only the high side is clamped; `y_end >= 0` by invariant.
   */
  int render_y_end() const { return y_end + margin < h ? y_end + margin : h; }
  /**
   * @brief Render-region left edge: display left expanded by `margin`, wrapped mod w (cylindrical).
   * @return First render column, in pixels, in [0, w).
   * @pre 0 <= x_start <= w and 0 <= margin < w, so `x_start - margin` is
   *      within one period of [0, w).
   */
  int render_x_start() const {
    const int v = x_start - margin;
    const int lo = v < 0 ? v + w : v;
    return lo >= w ? lo - w : lo;
  }
  /**
   * @brief Render-region right edge: display right expanded by `margin`, wrapped mod w (cylindrical).
   * @return One-past-last render column, in pixels, in [0, w).
   * @pre 0 <= x_end <= w and 0 <= margin < w.
   */
  int render_x_end() const {
    const int v = x_end + margin;
    return v >= w ? v - w : v;
  }

  /**
   * @brief Reports whether this region covers the entire canvas (no clipping).
   * @return True when display bounds equal the full [0,w) x [0,h) canvas.
   */
  bool is_full() const {
    return y_start == 0 && y_end == h && x_start == 0 && x_end == w;
  }

  /**
   * @brief Pixel-level vertical containment against the render (margin-expanded) bounds.
   * @param y Row index, in pixels.
   * @return True when y lies within [render_y_start(), render_y_end()).
   */
  bool contains_y(int y) const {
    return y >= render_y_start() && y < render_y_end();
  }

  /**
   * @brief Pixel-level horizontal containment against the render (margin-expanded) bounds.
   * @param x Column index, in pixels.
   * @return True when x lies within the cylindrical render band.
   * @details Once full coverage is excluded, coinciding ends (rs == re) mean a
   *          zero-width band.
   */
  bool contains_x(int x) const {
    if (covers_all_columns())
      return true;
    int rs = render_x_start();
    int re = render_x_end();
    if (rs == re)
      return false; // empty band
    return (rs < re) ? (x >= rs && x < re) : (x >= rs || x < re);
  }

  /**
   * @brief Precomputed cylindrical x-clip predicate built once per draw, then queried per fragment.
   * @details The full-coverage case (display width + both margins >= w) folds
   *          into `active == false`, so hot loops skip the test entirely. A
   *          zero-width band leaves `active` set with rs == re and wrap false,
   *          which clips every column.
   */
  struct XClip {
    int rs = 0; /**< Render band start column, in pixels, in [0, w). */
    int re =
        0; /**< Render band end column (exclusive), in pixels, in [0, w). */
    bool active = false; /**< False => no x clipping (full coverage). */
    bool wrap = false;   /**< Band crosses the seam (rs > re). */

    /**
     * @brief Tests whether a fragment column falls outside the render band.
     * @param x Column index, in pixels.
     * @return True when x lies outside the render band and must be skipped.
     */
    bool clipped(int x) const {
      if (!active)
        return false;
      return wrap ? (x < rs && x >= re) : (x < rs || x >= re);
    }

    /**
     * @brief Band length in columns, seam-unwrapped.
     * @param w Cylinder width in columns.
     * @return Column count spanned by [rs, re), counting past the seam when the
     *         band wraps, and the full w when no x clipping applies.
     */
    constexpr int length(int w) const {
      if (!active)
        return w;
      return wrap ? re - rs + w : re - rs;
    }
  };

  /**
   * @brief Builds the precomputed cylindrical x-clip predicate for this region.
   * @return An XClip whose `active` flag is false for full-coverage bands.
   */
  XClip x_clip() const {
    XClip c;
    c.rs = render_x_start();
    c.re = render_x_end();
    c.active = !covers_all_columns();
    c.wrap = c.rs > c.re;
    return c;
  }

  /**
   * @brief Conservative AABB test for whether a screen-space segment touches the render region.
   * @param y1 First endpoint's row coordinate, in pixels.
   * @param y2 Second endpoint's row coordinate, in pixels.
   * @return True if the segment's row range overlaps the render region.
   * @details Does not include splat reach; callers must pad the range for
   * downstream filters before testing overlap.
   */
  bool could_intersect_y(float y1, float y2) const {
    float lo = fminf(y1, y2);
    float hi = fmaxf(y1, y2);
    return hi >= render_y_start() && lo < render_y_end();
  }

  /**
   * @brief Tests whether two cylindrical column arcs overlap.
   * @param s1 First arc start column, in [0, w).
   * @param len1 First arc length in columns.
   * @param s2 Second arc start column, in [0, w).
   * @param len2 Second arc length in columns.
   * @param w Cylinder width in columns.
   * @return True if the arcs share at least one column.
   * @pre w > 0 and both starts in [0, w).
   */
  static bool arcs_overlap(int s1, int len1, int s2, int len2, int w) {
    if (len1 <= 0 || len2 <= 0)
      return false;
    if (len1 >= w || len2 >= w)
      return true;
    HS_AUDIT_CHECK(s1 >= 0 && s1 < w && s2 >= 0 && s2 < w);
    auto covers = [w](int s, int len, int p) {
      int d = p - s;
      if (d < 0)
        d += w;
      return d < len;
    };
    return covers(s1, len1, s2) || covers(s2, len2, s1);
  }

private:
  /**
   * @brief Full-coverage predicate.
   * @return True when the render band (display width plus both margins) spans
   *         the whole cylinder, so no x clipping applies.
   */
  bool covers_all_columns() const {
    return (x_end - x_start) + 2 * margin >= w;
  }
};
