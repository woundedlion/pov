/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "color/color.h"
#include <cassert>

/**
 * @file layer_composite.h
 * @brief Front-to-back "over" accumulator for layered coverage.
 */

/**
 * @brief Accumulates partially covering layers front to back into one Color4.
 * @details Each add() weights its color by the transmittance still left in
 * front of it, so callers may feed layers in depth order without tracking the
 * running alpha themselves. `remaining` is that transmittance; saturated() is
 * the early-out test over it. finish() un-premultiplies the sum, so the
 * returned color is the blended hue and the returned alpha its coverage.
 */
struct LayerComposite {
  float red = 0.0f;
  float green = 0.0f;
  float blue = 0.0f;
  float remaining = 1.0f; /**< Transmittance left in front of the next layer. */

  /**
   * @brief Adds one layer behind everything added so far.
   * @param color Layer color.
   * @param coverage Fraction of the remaining transmittance the layer covers,
   *        in [0, 1]; above 1 drives `remaining` negative and every later add()
   *        subtracts.
   */
  void add(const Pixel &color, float coverage) {
    assert(coverage >= 0.0f && coverage <= 1.0f);
    const float weight = remaining * coverage;
    red += static_cast<float>(color.r) * weight;
    green += static_cast<float>(color.g) * weight;
    blue += static_cast<float>(color.b) * weight;
    remaining *= 1.0f - coverage;
  }

  /**
   * @brief Whether the accumulated coverage leaves anything for a later layer.
   * @return True once the remaining transmittance is at or below the smallest
   *   encodable alpha, where no further add() can change finish().
   */
  bool saturated() const { return remaining <= MIN_ENCODABLE_ALPHA; }

  /** @brief Un-premultiplied color and its accumulated alpha. */
  Color4 finish() const {
    const float alpha = 1.0f - remaining;
    if (alpha <= 0.0f)
      return {};
    const float inverse_alpha = 1.0f / alpha;
    const Pixel color{static_cast<uint16_t>(hs::clamp(
                          red * inverse_alpha + 0.5f, 0.0f, 65535.0f)),
                      static_cast<uint16_t>(hs::clamp(
                          green * inverse_alpha + 0.5f, 0.0f, 65535.0f)),
                      static_cast<uint16_t>(hs::clamp(
                          blue * inverse_alpha + 0.5f, 0.0f, 65535.0f))};
    return {color, alpha};
  }
};
