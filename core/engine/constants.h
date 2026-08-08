/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file constants.h
 * @brief Engine-wide compile-time constants: canvas bounds, star geometry,
 *        pole-LOD tuning, and ClipRegion.
 */

#include "engine/platform.h" // CANVAS_W/CANVAS_H

/**
 * @brief Inner/outer radius ratio for star shapes (1/φ² ≈ 0.382).
 */
inline constexpr float STAR_INNER_RATIO = 0.382f;

/**
 * @brief Maximum horizontal resolution (width) for effects, from CANVAS_W.
 * @details Sizes the shared framebuffers (Effect::buffer_a/buffer_b) and bounds
 *          the Effect constructor, so a target that renders one small canvas
 *          reserves only what it draws (-DCANVAS_W/-DCANVAS_H per env). Builds
 *          that instantiate several resolutions — host tests, WASM, Phantasm —
 *          keep the 288x144 default, which bounds every one of them.
 */
inline constexpr int MAX_W = CANVAS_W;

/**
 * @brief Maximum vertical resolution (height) for effects, from CANVAS_H.
 */
inline constexpr int MAX_H = CANVAS_H;

/**
 * @brief Longest column run a single shade may be splatted across.
 * @details Bounds the near-pole run so one shade never covers a visually
 *          significant arc, and keeps a run short relative to the narrowest
 *          clip segment.
 */
inline constexpr int POLE_LOD_MAX_RUN = 32;

/**
 * @brief Aggressiveness of near-pole azimuthal shading decimation; 0 disables.
 * @details A row at colatitude phi has horizontal pixel pitch sin(phi) times
 *          the vertical, so 1/sin(phi) columns share one physical LED
 *          footprint and need only one shade between them. The column run is
 *          `aggressiveness / sin(phi)`, so 1.0 tracks that footprint exactly
 *          and smaller values stay inside it. At 0 every run is one column and
 *          the scan is bit-identical to an undecimated walk.
 *
 *          The true masking width depends on the LED's angular size and the
 *          per-column exposure, so this is a hardware-calibrated knob rather
 *          than a derived constant. Firmware has no setter, so the starting
 *          value comes from HS_POLE_LOD_DEFAULT.
 */
#ifndef HS_POLE_LOD_DEFAULT
#define HS_POLE_LOD_DEFAULT 0.0f
#endif
#ifdef ARDUINO
inline constexpr float pole_lod_aggressiveness = HS_POLE_LOD_DEFAULT;
inline constexpr bool POLE_LOD_ENABLED = HS_POLE_LOD_DEFAULT > 0.0f;
#else
inline float pole_lod_aggressiveness = HS_POLE_LOD_DEFAULT;
inline constexpr bool POLE_LOD_ENABLED = true;
#endif

/**
 * @brief Clip region for segment-based rendering.
 * @details Display bounds define the ISR's pixel range (exact segment).
 *          Render bounds expand by `margin` to accommodate stateful filters
 *          (e.g. AntiAlias ±1, MeshFeedback noise warp ±8).
 *          `w`/`h` are the active canvas size, set from W/H by the Effect
 *          constructor (core/render/canvas.h); the `MAX_W`/`MAX_H` defaults
 *          apply only to a standalone ClipRegion.
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
   * @details Only the low side is clamped: `y_start - margin` can underflow past
   *          0, but no high-side clamp is needed because `y_start <= h` holds by
   *          the driver's display-bounds invariant (with `margin >= 0`). That same
   *          invariant gives `y_start <= y_end`, so paired with render_y_end the
   *          range never inverts; only the two single-sided clamps are needed.
   */
  int render_y_start() const {
    return y_start - margin > 0 ? y_start - margin : 0;
  }
  /**
   * @brief Render-region bottom edge: display bottom expanded down by `margin`, capped at h.
   * @return One-past-last render row (exclusive), in pixels.
   * @details Mirror of render_y_start: only the high side is clamped (to h);
   *          `y_end >= 0` makes a low-side clamp unnecessary.
   */
  int render_y_end() const { return y_end + margin < h ? y_end + margin : h; }
  /**
   * @brief Render-region left edge: display left expanded by `margin`, wrapped mod w (cylindrical).
   * @return First render column, in pixels, in [0, w).
   * @pre 0 <= x_start <= w and 0 <= margin < w, as Canvas::set_clip /
   *      Canvas::set_margin enforce. That puts `x_start - margin` in
   *      [-(w-1), w], one period either side, so the wrap is a conditional add
   *      plus a conditional subtract instead of a `%` (the high branch fires
   *      only at x_start == w with margin == 0). contains_x() runs per plotted
   *      pixel, where a runtime modulo costs a hardware divide.
   */
  int render_x_start() const {
    const int v = x_start - margin;
    const int lo = v < 0 ? v + w : v;
    return lo >= w ? lo - w : lo;
  }
  /**
   * @brief Render-region right edge: display right expanded by `margin`, wrapped mod w (cylindrical).
   * @return One-past-last render column, in pixels, in [0, w).
   * @pre 0 <= x_end <= w and 0 <= margin < w (see render_x_start); `x_end +
   *      margin` is then in [0, 2w-1], so one conditional subtract wraps it.
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
   * @details A render band spanning >= w columns (display width + both margins)
   *          covers everything; past that test the band is a proper sub-arc, so
   *          ends that coincide (rs == re) mean zero width and cover nothing;
   *          otherwise the band may cross the seam, so test as a wrapped
   *          interval.
   */
  bool contains_x(int x) const {
    if ((x_end - x_start) + 2 * margin >= w)
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
     *         band wraps. Feeds arcs_overlap as the clip arc's length.
     */
    constexpr int length(int w) const { return wrap ? re - rs + w : re - rs; }
  };

  /**
   * @brief Builds the precomputed cylindrical x-clip predicate for this region.
   * @return An XClip whose `active` flag is false for full-coverage bands.
   */
  XClip x_clip() const {
    XClip c;
    c.rs = render_x_start();
    c.re = render_x_end();
    c.active = (x_end - x_start) + 2 * margin < w;
    c.wrap = c.rs > c.re;
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

  /**
   * @brief Tests whether two cylindrical column arcs overlap.
   * @param s1 First arc start column, in [0, w).
   * @param len1 First arc length in columns.
   * @param s2 Second arc start column, in [0, w).
   * @param len2 Second arc length in columns.
   * @param w Cylinder width in columns.
   * @return True if the arcs share at least one column.
   * @pre w > 0 and both starts in [0, w), as every caller's column mapping
   *      guarantees. The seam-relative offset is then in (-w, w), so one
   *      conditional add wraps it instead of a `%`; this runs per geodesic edge
   *      and per face azimuth interval, where a runtime modulo costs a hardware
   *      divide.
   */
  static bool arcs_overlap(int s1, int len1, int s2, int len2, int w) {
    if (len1 >= w || len2 >= w)
      return true;
    if (len1 <= 0 || len2 <= 0)
      return false;
    HS_TEST_CHECK(s1 >= 0 && s1 < w && s2 >= 0 && s2 < w);
    auto covers = [w](int s, int len, int p) {
      int d = p - s;
      if (d < 0)
        d += w;
      return d < len;
    };
    return covers(s1, len1, s2) || covers(s2, len2, s1);
  }
};
