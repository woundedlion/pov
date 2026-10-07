/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

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
 * @tparam Spread Columns per fringe step; scale it with W to hold the
 *         fringe's angular width across resolutions.
 */
template <int W, int Spread = 1> class ChromaticShift : public IsPixel {
  static_assert(Spread >= 1, "ChromaticShift requires a positive Spread");
  // fast_wrap corrects only one ±W step, so the three fringe offsets stay in a
  // single wrap of [0,W) only while 3 * Spread < W.
  static_assert(W >= 3 * Spread + 1,
                "ChromaticShift requires W > 3 * Spread for fast_wrap offsets");

public:
  /** @brief Fringe taps land up to 3 * Spread columns from the plotted one. */
  static constexpr int segment_margin = 3 * Spread;
  /** @brief Constructs the chromatic-shift filter (stateless). */
  ChromaticShift() {}

  /**
   * @brief Emits the source pixel plus R/G/B copies offset by 1/2/3 spreads.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param c Source color; split into single-channel copies.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Source blend alpha in [0, 1]; fringe taps use one quarter.
   * @tparam PassFnT Downstream callback type.
   * @param pass Downstream 2D callback.
   */
  template <typename PassFnT>
  void plot(float x, float y, const ::Pixel &c, float age, float alpha,
            PassFnT &&pass) {
    assert(age >= 0.0f && alpha >= 0.0f);
    const int xi = round_wrap_column<W>(x);
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

    const float fringe_alpha = alpha * 0.25f;
    pass(static_cast<float>(math::fast_wrap(xi + Spread, W)), y, r_col, age,
         fringe_alpha);
    pass(static_cast<float>(math::fast_wrap(xi + 2 * Spread, W)), y, g_col, age,
         fringe_alpha);
    pass(static_cast<float>(math::fast_wrap(xi + 3 * Spread, W)), y, b_col, age,
         fringe_alpha);
  }
};

} // namespace Pixel
} // namespace Filter
