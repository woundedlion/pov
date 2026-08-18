/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>
#include "render/filter/pipeline.h"
#include "math/geometry.h"
#include "engine/util.h"

/**
 * @file splat.h
 * @brief SplatTaps and splat_taps(): resolves a sub-pixel sample into its four
 * quintic-eased nearest-neighbor taps, shared by the anti-alias stages.
 */

namespace Filter {
namespace Screen {

HS_O3_BEGIN

/** @brief One sub-pixel sample resolved into its four nearest-neighbor taps. */
struct SplatTaps {
  int x0, x1;       /**< Wrapped left and right columns. */
  int y0, y1;       /**< Top and bottom rows, either may lie outside [0, H). */
  bool y0_physical; /**< y0 lies in [0, H). */
  bool y1_physical; /**< y1 lies in [0, H). */
  float v00, v10, v01,
      v11; /**< Bilinear coverage per tap, row-major from (x0, y0). */
};

/**
 * @brief Splat weight below which a tap contributes nothing worth emitting.
 * @details Looser in Blur (1e-5): these are raw bilinear coverage products,
 * Blur's are normalized 3x3 kernel taps.
 */
static constexpr float SPLAT_TAP_CUTOFF = 1e-8f;

/**
 * @brief Resolves a sub-pixel sample into its four nearest-neighbor taps.
 * @tparam W Framebuffer width in pixels.
 * @tparam H Framebuffer height in pixels.
 * @param x Sub-pixel column coordinate.
 * @param y Sub-pixel row coordinate.
 * @return Wrapped columns, rows, row-in-frame flags, and the four tap weights.
 * @details Both axes are eased with a quintic kernel; the splat is uniform in
 * framebuffer space at every latitude (no sin(phi) density compensation). A row
 * off the top or bottom edge donates its whole weight to the surviving row.
 */
template <int W, int H>
__attribute__((always_inline)) inline SplatTaps splat_taps(float x, float y) {
  // Non-finite coords make the int casts below UB and bypass the wrap.
  assert(std::isfinite(x) && std::isfinite(y));
  // fast_wrap below corrects only a single ±W offset on floorf(x).
  assert(x >= -W && x < 2 * W);
  // y never wraps; bounded only so the cast below stays in range.
  assert(y >= -H && y < 2 * H);

  const float y_floor = floorf(y);
  const float x_floor = floorf(x);
  const float xs = quintic_kernel(x - x_floor);
  const float ys = quintic_kernel(y - y_floor);

  SplatTaps t;
  t.y0 = static_cast<int>(y_floor);
  t.y1 = t.y0 + 1;
  t.x0 = fast_wrap(static_cast<int>(x_floor), W);
  t.x1 = fast_wrap(t.x0 + 1, W);
  t.y0_physical = t.y0 >= 0 && t.y0 < H;
  t.y1_physical = t.y1 >= 0 && t.y1 < H;

  float wy0 = 1.0f - ys;
  float wy1 = ys;
  if (t.y0_physical && !t.y1_physical) {
    wy0 = 1.0f;
    wy1 = 0.0f;
  } else if (!t.y0_physical && t.y1_physical) {
    wy0 = 0.0f;
    wy1 = 1.0f;
  }

  t.v00 = (1.0f - xs) * wy0;
  t.v10 = xs * wy0;
  t.v01 = (1.0f - xs) * wy1;
  t.v11 = xs * wy1;
  return t;
}
HS_O3_END

} // namespace Screen
} // namespace Filter
