/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#ifndef HOLOSPHERE_CORE_CONSTANTS_H_
#define HOLOSPHERE_CORE_CONSTANTS_H_

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

  // Render bounds (expanded by margin, clamped to canvas)
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
    int rs = render_x_start();
    int re = render_x_end();
    return (rs <= re) ? (x >= rs && x < re) : (x >= rs || x < re);
  }

  /// Conservative AABB test: could a segment between two screen points
  /// produce pixels inside the render region?
  bool could_intersect_y(float y1, float y2) const {
    float lo = y1 < y2 ? y1 : y2;
    float hi = y1 > y2 ? y1 : y2;
    return hi >= render_y_start() && lo < render_y_end();
  }
};

#endif // HOLOSPHERE_CORE_CONSTANTS_H_
