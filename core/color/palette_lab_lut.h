/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/color/baked_palette.h"
#include "core/color/color_space.h"
#include <algorithm>
#include <array>
#include <cassert>

/** @file palette_lab_lut.h
 * @brief Per-interval RGB-first palette blends cached in OKLab. */

/** @brief Borrowed neighboring time rows for opaque palettes. */
struct PaletteLabTimeView {
  const OKLab *low = nullptr;  ///< Lower time row.
  const OKLab *high = nullptr; ///< Upper time row.
  float weight = 0.0f;         ///< Upper time weight.
};

/** @brief Owned palette-coordinate by frame-time base-color cache. */
struct PaletteLabLut {
  static constexpr int VALUES = 64;          ///< Palette-coordinate nodes.
  static constexpr int TIMES = 3;            ///< Frame-time nodes.
  std::array<OKLab, VALUES * TIMES> nodes{}; ///< RGB-first blended colors.
  bool active = false; ///< Whether both palettes are opaque.

  /** @brief Bakes RGB16 endpoint blends before converting each node to OKLab.
   * @param from Previous immutable palette.
   * @param to Current immutable palette. */
  HS_COLD_MEMBER void bake(const BakedPalette &from, const BakedPalette &to) {
    const auto A = from.constant_alpha();
    const auto B = to.constant_alpha();
    active = A.has_value() && B.has_value() && *A == 1.0f && *B == 1.0f;
    if (!active)
      return;
    for (int time = 0; time < TIMES; ++time) {
      const uint16_t WEIGHT = frac_to_q16(float(time) / (TIMES - 1));
      for (int value = 0; value < VALUES; ++value) {
        const float COORDINATE = float(value) / (VALUES - 1);
        const Pixel RGB =
            from.get_color(COORDINATE).lerp16(to.get_color(COORDINATE), WEIGHT);
        const LinRGB LINEAR = pixel_to_linrgb(RGB);
        nodes[time * VALUES + value] =
            linear_rgb_to_oklab_fast(LINEAR.r, LINEAR.g, LINEAR.b);
      }
    }
  }

  /** @brief Prepares the two time rows for one column.
   * @param fraction Frame progress in [0, 1].
   * @return Borrowed rows, or an inactive view. */
  HS_FLASH_INLINE PaletteLabTimeView view(float fraction) const {
    if (!active)
      return {};
    const float POSITION = hs::clamp(fraction, 0.0f, 1.0f) * (TIMES - 1);
    const int LOW_INDEX = std::min(int(POSITION), TIMES - 2);
    return {nodes.data() + LOW_INDEX * VALUES,
            nodes.data() + (LOW_INDEX + 1) * VALUES, POSITION - LOW_INDEX};
  }

  /** @brief Samples neighboring palette nodes from both prepared time rows.
   * @param time Prepared column rows.
   * @param value Palette coordinate in [0, 1].
   * @return Bilinearly interpolated base OKLab. */
  __attribute__((always_inline)) static OKLab
  sample(const PaletteLabTimeView &time, float value) {
    assert(time.low && time.high);
    const float POSITION = hs::clamp(value, 0.0f, 1.0f) * (VALUES - 1);
    const int LOW_INDEX = std::min(int(POSITION), VALUES - 2);
    const float WEIGHT = POSITION - LOW_INDEX;
    const auto lerp = [](const OKLab &a, const OKLab &b, float fraction) {
      return OKLab{hs::lerp(a.L, b.L, fraction), hs::lerp(a.a, b.a, fraction),
                   hs::lerp(a.b, b.b, fraction)};
    };
    return lerp(lerp(time.low[LOW_INDEX], time.low[LOW_INDEX + 1], WEIGHT),
                lerp(time.high[LOW_INDEX], time.high[LOW_INDEX + 1], WEIGHT),
                time.weight);
  }
};
