/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/filter/splat.h"
#include "color/color.h"

/**
 * @file screen_anti_alias.h
 * @brief Filter::Screen::AntiAlias: splats a sub-pixel sample across its four
 * nearest pixel neighbors.
 */

namespace Filter {
namespace Screen {

HS_O3_BEGIN

/**
 * @brief Applies 2D anti-aliasing to sub-pixel coordinates.
 * @details Distributes intensity to the 4 nearest neighbors using a quintic kernel.
 */
template <int W, int H> class AntiAlias : public Is2D {
public:
  static constexpr int segment_margin = 1;

  /**
   * @brief Splats a sub-pixel sample across its four nearest pixel neighbors.
   * @tparam PassFnT Downstream 2D callback type.
   * @param x Sub-pixel column coordinate.
   * @param y Sub-pixel row coordinate.
   * @param c Source color, forwarded to each tap.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1]; scaled per tap by its quintic-eased splat weight.
   * @param pass Downstream 2D callback receiving each weighted tap.
   * @details @p pass is a forwarding-reference template so the densest fan-out
   * in the family inlines its taps.
   */
  template <typename PassFnT>
  void plot(float x, float y, const ::Pixel &c, float age, float alpha,
            PassFnT &&pass) {
    assert(age >= 0.0f && alpha >= 0.0f);
    const SplatTaps t = splat_taps<W, H>(x, y);

    if (t.y0_physical && t.v00 > SPLAT_TAP_CUTOFF)
      pass(static_cast<float>(t.x0), static_cast<float>(t.y0), c, age,
           alpha * t.v00);
    if (t.y0_physical && t.v10 > SPLAT_TAP_CUTOFF)
      pass(static_cast<float>(t.x1), static_cast<float>(t.y0), c, age,
           alpha * t.v10);
    if (t.y1_physical && t.v01 > SPLAT_TAP_CUTOFF)
      pass(static_cast<float>(t.x0), static_cast<float>(t.y1), c, age,
           alpha * t.v01);
    if (t.y1_physical && t.v11 > SPLAT_TAP_CUTOFF)
      pass(static_cast<float>(t.x1), static_cast<float>(t.y1), c, age,
           alpha * t.v11);
  }
};
HS_O3_END

} // namespace Screen
} // namespace Filter
