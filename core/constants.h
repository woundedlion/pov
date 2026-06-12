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
 *
 * Display bounds define the ISR's pixel range (exact segment).
 * Render bounds expand by `margin` to accommodate stateful filters
 * (e.g. AntiAlias ±1, MeshFeedback noise warp ±8).
 */
struct ClipRegion {
  int y_start = 0;
  int y_end   = MAX_H;
  int x_start = 0;
  int x_end   = MAX_W;
  int margin  = 1;
  int w       = MAX_W;
  int h       = MAX_H;

  // Render bounds (expanded by margin): y clamped to [0,h], x wrapped mod w (cylindrical)
  int render_y_start() const { return y_start - margin > 0 ? y_start - margin : 0; }
  int render_y_end()   const { return y_end + margin < h ? y_end + margin : h; }
  int render_x_start() const { return (x_start - margin + w) % w; }
  int render_x_end()   const { return (x_end + margin) % w; }

  // Full canvas (no clipping)
  bool is_full() const { return y_start == 0 && y_end == h && x_start == 0 && x_end == w; }

  // Pixel-level containment (render bounds)
  bool contains_y(int y) const {
    return y >= render_y_start() && y < render_y_end();
  }

  bool contains_x(int x) const {
    if (x_end - x_start >= w) return true;
    int rs = render_x_start();
    int re = render_x_end();
    if (rs == re) return true; // margin expansion wrapped to full width
    return (rs < re) ? (x >= rs && x < re) : (x >= rs || x < re);
  }

  /// Precomputed cylindrical x-clip predicate. Construct once per draw, then
  /// call clipped(x) per fragment. The full-width case and the wrap-to-full
  /// case (a sub-arc whose margin expansion sums to exactly w, so rs == re)
  /// both fold into `active == false`, so hot loops can't accidentally blank a
  /// full segment by treating rs == re as an empty band.
  struct XClip {
    int  rs     = 0;
    int  re     = 0;
    bool active = false; // false => no x clipping (full width or wrapped to full)
    bool wrap   = false; // band crosses the seam (rs > re)

    /// True when x lies outside the render band and must be skipped.
    bool clipped(int x) const {
      if (!active) return false;
      return wrap ? (x < rs && x >= re) : (x < rs || x >= re);
    }
  };

  XClip x_clip() const {
    XClip c;
    c.rs     = render_x_start();
    c.re     = render_x_end();
    c.active = (x_end - x_start) < w && c.rs != c.re;
    c.wrap   = c.rs > c.re;
    return c;
  }

  /// Conservative AABB test: could a segment between two screen points
  /// produce pixels inside the render region?
  bool could_intersect_y(float y1, float y2) const {
    float lo = y1 < y2 ? y1 : y2;
    float hi = y1 > y2 ? y1 : y2;
    return hi >= render_y_start() && lo < render_y_end();
  }
};

