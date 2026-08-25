/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>
#include "render/filter/pipeline.h"
#include "color/color.h"
#include "math/geometry.h"

/**
 * @file pixel_chromatic_shift.h
 * @brief Filter::Pixel::ChromaticShift: splits RGB into per-channel copies
 * offset by 1/2/3 spreads of columns, producing a chromatic-aberration fringe.
 */

namespace Filter {
namespace Pixel {

/**
 * @brief Splits RGB into per-channel copies offset by 1/2/3 spreads of columns,
 * producing a chromatic-aberration fringe.
 * @tparam W Canvas width in columns.
 * @tparam Spread Columns per fringe step. The fringe subtends
 *         3 * Spread / W of a turn, so scaling it with W holds the aberration
 *         at a fixed angular width across resolutions.
 */
template <int W, int Spread = 1> class ChromaticShift : public IsPixel {
  static_assert(Spread >= 1, "ChromaticShift requires a positive Spread");
  // fast_wrap corrects only one ±W step, so the three fringe offsets stay in a
  // single wrap of [0,W) only while 3 * Spread < W.
  static_assert(W >= 3 * Spread + 1,
                "ChromaticShift requires W > 3 * Spread for fast_wrap offsets");

public:
  // emits_pixel_centers is withheld: only the three fringe taps round, and only
  // their column. The source tap keeps its sub-pixel x and every tap keeps the
  // caller's y, so a downstream sub-pixel stage is not reduced to an identity.
  /**
   * @brief The three fringe taps land outside the plotted position, so a
   *        segment worker needs 3 * Spread columns of render margin to write
   *        them.
   */
  static constexpr int segment_margin = 3 * Spread;
  /** @brief Constructs the chromatic-shift filter (stateless). */
  ChromaticShift() {}

  /**
   * @brief Emits the source pixel plus R/G/B copies offset by 1/2/3 spreads.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param c Source color; split into single-channel copies.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 2D callback.
   */
  template <typename PassFnT>
  void plot(float x, float y, const ::Pixel &c, float age, float alpha,
            PassFnT &&pass) {
    assert(age >= 0.0f && alpha >= 0.0f);
    // Non-finite x makes the int cast below UB and bypasses the wrap.
    assert(std::isfinite(x));
    assert(x > -W - 0.5f && x < 2 * W - 0.5f);
    pass(x, y, c, age, alpha);

    ::Pixel r_col = c;
    r_col.g = 0;
    r_col.b = 0;
    ::Pixel g_col = c;
    g_col.r = 0;
    g_col.b = 0;
    ::Pixel b_col = c;
    b_col.r = 0;
    b_col.g = 0;

    const float xr = std::round(x);
    // fast_wrap corrects only a single ±W offset, so xr must land in [-W, 2W).
    assert(xr >= -W && xr < 2 * W);
    int xi = fast_wrap(static_cast<int>(xr), W);
    pass(static_cast<float>(fast_wrap(xi + Spread, W)), y, r_col, age, alpha);
    pass(static_cast<float>(fast_wrap(xi + 2 * Spread, W)), y, g_col, age,
         alpha);
    pass(static_cast<float>(fast_wrap(xi + 3 * Spread, W)), y, b_col, age,
         alpha);
  }
};

} // namespace Pixel
} // namespace Filter
