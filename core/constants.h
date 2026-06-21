/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

/**
 * @brief Maximum horizontal resolution (width) for effects.
 */
static constexpr int MAX_W = 288;

/**
 * @brief Maximum vertical resolution (height) for effects.
 */
static constexpr int MAX_H = 144;

/**
 * @brief Clip region for segment-based rendering.
 * @details Display bounds define the ISR's pixel range (exact segment).
 *          Render bounds expand by `margin` to accommodate stateful filters
 *          (e.g. AntiAlias ±1, MeshFeedback noise warp ±8).
 */
struct ClipRegion {
  int y_start = 0;      /**< Display top row (inclusive), in pixels. */
  int y_end   = MAX_H;  /**< Display bottom row (exclusive), in pixels. */
  int x_start = 0;      /**< Display left column (inclusive), in pixels. */
  int x_end   = MAX_W;  /**< Display right column (exclusive), in pixels. */
  int margin  = 1;      /**< Render-bound expansion past the display edges, in pixels. */
  int w       = MAX_W;  /**< Canvas width, in pixels. */
  int h       = MAX_H;  /**< Canvas height, in pixels. */

  /**
   * @brief Render-region top edge: display top expanded up by `margin`, clamped to [0,h].
   * @return First render row (inclusive), in pixels.
   */
  int render_y_start() const { return y_start - margin > 0 ? y_start - margin : 0; }
  /**
   * @brief Render-region bottom edge: display bottom expanded down by `margin`, clamped to [0,h].
   * @return One-past-last render row (exclusive), in pixels.
   */
  int render_y_end()   const { return y_end + margin < h ? y_end + margin : h; }
  /**
   * @brief Render-region left edge: display left expanded by `margin`, wrapped mod w (cylindrical).
   * @return First render column, in pixels, in [0, w).
   * @pre margin < w. The single `+ w` corrects only one period of underflow, so
   *      the [0, w) result holds only while margin <= x_start + w;
   *      Canvas::set_margin traps a margin >= w at configuration time. Kept
   *      single-period (not a double-mod) to avoid extra work on the
   *      per-fragment contains_x() path.
   */
  int render_x_start() const { return (x_start - margin + w) % w; }
  /**
   * @brief Render-region right edge: display right expanded by `margin`, wrapped mod w (cylindrical).
   * @return One-past-last render column, in pixels, in [0, w).
   * @pre margin < w (see render_x_start).
   */
  int render_x_end()   const { return (x_end + margin) % w; }

  /**
   * @brief Reports whether this region covers the entire canvas (no clipping).
   * @return True when display bounds equal the full [0,w) x [0,h) canvas.
   */
  bool is_full() const { return y_start == 0 && y_end == h && x_start == 0 && x_end == w; }

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
   * @details Bands wider than w cover everything; a band whose ends coincide
   *          (rs == re) is a full-width wrap; otherwise the band may cross the
   *          seam, so test as a wrapped interval.
   */
  bool contains_x(int x) const {
    if (x_end - x_start >= w) return true;
    int rs = render_x_start();
    int re = render_x_end();
    if (rs == re) return true; // margin expansion wrapped to full width
    return (rs < re) ? (x >= rs && x < re) : (x >= rs || x < re);
  }

  /**
   * @brief Precomputed cylindrical x-clip predicate built once per draw, then queried per fragment.
   * @details The full-width case and the wrap-to-full case (a sub-arc whose
   *          margin expansion sums to exactly w, so rs == re) both fold into
   *          `active == false`, so hot loops can't accidentally blank a full
   *          segment by treating rs == re as an empty band.
   */
  struct XClip {
    int  rs     = 0;     /**< Render band start column, in pixels, in [0, w). */
    int  re     = 0;     /**< Render band end column (exclusive), in pixels, in [0, w). */
    bool active = false; /**< False => no x clipping (full width or wrapped to full). */
    bool wrap   = false; /**< Band crosses the seam (rs > re). */

    /**
     * @brief Tests whether a fragment column falls outside the render band.
     * @param x Column index, in pixels.
     * @return True when x lies outside the render band and must be skipped.
     */
    bool clipped(int x) const {
      if (!active) return false;
      return wrap ? (x < rs && x >= re) : (x < rs || x >= re);
    }
  };

  /**
   * @brief Builds the precomputed cylindrical x-clip predicate for this region.
   * @return An XClip whose `active` flag is false for full-width or wrap-to-full bands.
   */
  XClip x_clip() const {
    XClip c;
    c.rs     = render_x_start();
    c.re     = render_x_end();
    c.active = (x_end - x_start) < w && c.rs != c.re;
    c.wrap   = c.rs > c.re;
    return c;
  }

  /**
   * @brief Conservative AABB test for whether a screen-space segment touches the render region.
   * @param y1 First endpoint's row coordinate, in pixels.
   * @param y2 Second endpoint's row coordinate, in pixels.
   * @return True if a segment between the two points could produce pixels inside the render region.
   */
  bool could_intersect_y(float y1, float y2) const {
    float lo = y1 < y2 ? y1 : y2;
    float hi = y1 > y2 ? y1 : y2;
    return hi >= render_y_start() && lo < render_y_end();
  }
};

