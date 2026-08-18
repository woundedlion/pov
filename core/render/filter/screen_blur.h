/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <cmath>
#include "render/filter/pipeline.h"
#include "color/color.h"
#include "engine/util.h"

/**
 * @file screen_blur.h
 * @brief Filter::Screen::Blur: splats a sample across its 3x3 neighborhood
 * through a variable Gaussian kernel.
 */

namespace Filter {
namespace Screen {

/**
 * @brief Applies a variable 3x3 Gaussian blur.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> class Blur : public Is2D {
public:
  static constexpr int segment_margin = 1;

  /**
   * @brief Constructs a blur with the given initial strength.
   * @param factor Blur strength in [0, 1] (0 = identity, 1 = full Gaussian).
   */
  Blur(float factor = 1.0f) { update(factor); }

  /**
   * @brief Rebuilds the 3x3 kernel for a new blur strength.
   * @param factor Blur strength; clamped to [0, 1].
   */
  void update(float factor) {
    float f = hs::clamp(factor, 0.0f, 1.0f);
    float c = 1.0f - (0.75f * f);
    float e = 0.125f * f;
    float d = 0.0625f * f;

    kernel = {d, e, d, e, c, e, d, e, d};
  }

  /**
   * @brief Splats the sample across its 3x3 neighborhood weighted by the kernel.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param color Source color, forwarded to each tap.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1]; scaled per tap by its kernel weight.
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 2D callback.
   */
  template <typename PassFnT>
  void plot(float x, float y, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    // Non-finite coords make the int casts below UB and bypass the wrap.
    assert(std::isfinite(x) && std::isfinite(y));
    // y never wraps; bounded only so the cast below stays in range.
    assert(y >= -H && y < 2 * H);
    const float xr = std::round(x);
    // fast_wrap corrects only a single ±W offset, so xr must land in [-W, 2W).
    assert(xr >= -W && xr < 2 * W);
    int cx = fast_wrap(static_cast<int>(xr), W);
    int cy = static_cast<int>(std::round(y));

    float inv = 1.0f;
    if (cy - 1 < 0 || cy + 1 >= H) {
      float wsum = 0.0f;
      for (int dy = -1; dy <= 1; dy++) {
        int ny = cy + dy;
        if (ny >= 0 && ny < H) {
          int r = (dy + 1) * 3;
          wsum += kernel[r] + kernel[r + 1] + kernel[r + 2];
        }
      }
      if (wsum > 1e-5f)
        inv = 1.0f / wsum;
    }

    int k = 0;
    for (int dy = -1; dy <= 1; dy++) {
      int ny = cy + dy;

      if (ny >= 0 && ny < H) {
        for (int dx = -1; dx <= 1; dx++) {
          float weight = kernel[k++] * inv;
          // Normalized kernel tap; cutoff is tighter in AntiAlias (1e-8) whose
          // weights are raw bilinear coverage products.
          if (weight > 1e-5f) {
            pass(static_cast<float>(fast_wrap(cx + dx, W)),
                 static_cast<float>(ny), color, age, alpha * weight);
          }
        }
      } else {
        k += 3;
      }
    }
  }

private:
  std::array<float, 9> kernel; /**< Row-major 3x3 blur weights, summing to 1. */
};

} // namespace Screen
} // namespace Filter
