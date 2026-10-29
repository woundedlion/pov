/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/color/color_space.h"
#include <span>

/** @file hue_rotation_lut.h
 * @brief Palette-coordinate by hue-rotation acceleration. */

/** @brief Largest float below 1: the top of the half-open palette coordinate
 *  the hue-rotation LUT is baked over and sampled at. */
inline constexpr float UNIT_OPEN_MAX = 0x1.fffffep-1f;

/** @brief View over palette-coordinate by hue-rotation colors. */
struct HueRotationLutView {
  static constexpr int VALUE_STEPS = 64; ///< Palette-coordinate samples.
  static constexpr int HUE_STEPS = 16;   ///< Hue-rotation samples per turn.
  static_assert(HUE_STEPS > 0 && (HUE_STEPS & (HUE_STEPS - 1)) == 0,
                "hue mask wrap requires a power of two");
  static constexpr size_t SIZE = VALUE_STEPS * HUE_STEPS; ///< Entries.

  /// Borrowed table, indexed `value_index * HUE_STEPS + hue_index`.
  const Pixel *data;
  bool active; ///< Whether `data` holds a baked table to sample.
  const math::Vector *linear_data = nullptr; ///< Optional float copy of RGB16.
};

/**
 * @brief Bakes all palette-coordinate and hue-rotation combinations.
 * @tparam Source Palette source exposing Color4 get(float) const.
 * @param output Destination LUT.
 * @param palette Palette source to hue-rotate.
 */
template <typename Source>
HS_FLASH_INLINE __attribute__((aligned(32))) inline void
prepare_hue_rotation_lut(std::span<Pixel, HueRotationLutView::SIZE> output,
                         const Source &palette) {
  for (int value_index = 0; value_index < HueRotationLutView::VALUE_STEPS;
       ++value_index) {
    const float value =
        UNIT_OPEN_MAX * value_index / (HueRotationLutView::VALUE_STEPS - 1);
    const HueRotateBase base = make_hue_rotate_base(palette.get(value));
    for (int hue_index = 0; hue_index < HueRotationLutView::HUE_STEPS;
         ++hue_index) {
      const float amount =
          static_cast<float>(hue_index) / HueRotationLutView::HUE_STEPS;
      output[value_index * HueRotationLutView::HUE_STEPS + hue_index] =
          hue_rotate_lut_gamut(base, amount).color;
    }
  }
}

/** @brief Palette-axis coordinates shared by hue-rotation LUT samples. */
struct HueRotationValueSample {
  int low;         ///< Lower palette row.
  int high;        ///< Upper palette row.
  uint16_t weight; ///< Upper row weight in Q16.

  /** @brief Prepares a palette-axis sample.
   * @param value Palette coordinate in [0, 1]. */
  explicit HueRotationValueSample(float value) {
    const float position =
        hs::clamp(value * (HueRotationLutView::VALUE_STEPS - 1), 0.0f,
                  static_cast<float>(HueRotationLutView::VALUE_STEPS - 1));
    low = lut_index_lo(position);
    high = std::min(low + 1, HueRotationLutView::VALUE_STEPS - 1);
    weight = lut_index_weight(position, low);
  }
};

/** @brief Hue-axis coordinates shared by palette endpoint samples. */
struct HueRotationHueSample {
  int low;             ///< Lower hue column, wrapped into the table.
  int high;            ///< Next hue column, wrapped into the table.
  uint16_t weight;     ///< Upper hue column weight in Q16.
  float linear_weight; ///< Q16 weight converted to a float fraction.

  /** @brief Prepares a periodic hue-axis sample.
   * @param amount Hue rotation in turns. */
  explicit HueRotationHueSample(float amount) {
    const float POSITION =
        hs::clamp((amount - floorf(amount)) * HueRotationLutView::HUE_STEPS,
                  0.0f, static_cast<float>(HueRotationLutView::HUE_STEPS));
    const int STEP = lut_index_lo(POSITION);
    low = STEP & (HueRotationLutView::HUE_STEPS - 1);
    high = (low + 1) & (HueRotationLutView::HUE_STEPS - 1);
    weight = lut_index_weight(POSITION, STEP);
    linear_weight = weight * Q16_TO_UNIT;
  }
};

/** @brief Samples hue rotation at prepared palette-axis coordinates.
 * @param view Prepared hue-rotation LUT.
 * @param value Palette-axis coordinates.
 * @param hue Prepared hue-axis coordinates.
 * @return Hue-rotated linear RGB color. */
__attribute__((always_inline)) inline Pixel
sample_hue_rotation_lut(const HueRotationLutView &view,
                        const HueRotationValueSample &value,
                        const HueRotationHueSample &hue) {
  const auto sample_row = [&](int row) {
    const int offset = row * HueRotationLutView::HUE_STEPS;
    return view.data[offset + hue.low].lerp16(view.data[offset + hue.high],
                                              hue.weight);
  };
  return sample_row(value.low).lerp16(sample_row(value.high), value.weight);
}

/** @brief Samples hue rotation at a prepared palette coordinate.
 * @param view Prepared hue-rotation LUT.
 * @param value Palette-axis coordinates.
 * @param amount Hue rotation in turns.
 * @return Hue-rotated linear RGB color. */
__attribute__((always_inline)) inline Pixel
sample_hue_rotation_lut(const HueRotationLutView &view,
                        const HueRotationValueSample &value, float amount) {
  return sample_hue_rotation_lut(view, value, HueRotationHueSample(amount));
}

/** @brief Samples a hue-rotation LUT with interpolation over both axes.
 * @param view Prepared hue-rotation LUT.
 * @param value Palette coordinate in [0, 1].
 * @param amount Hue rotation in turns.
 * @return Hue-rotated linear RGB color. */
__attribute__((always_inline)) inline Pixel
sample_hue_rotation_lut(const HueRotationLutView &view, float value,
                        float amount) {
  return sample_hue_rotation_lut(view, HueRotationValueSample(value), amount);
}

/** @brief Samples hue rotation without rounding intermediate RGB blends.
 * @param view Prepared hue-rotation LUT.
 * @param value Palette-axis coordinates.
 * @param hue Prepared hue-axis coordinates.
 * @return Linear RGB channel values before quantization. */
__attribute__((always_inline)) inline math::Vector
sample_hue_rotation_linear(const HueRotationLutView &view,
                           const HueRotationValueSample &value,
                           const HueRotationHueSample &hue) {
  const auto row = [&](int index) {
    if (view.linear_data) {
      const math::Vector &a =
          view.linear_data[index * HueRotationLutView::HUE_STEPS + hue.low];
      const math::Vector &b =
          view.linear_data[index * HueRotationLutView::HUE_STEPS + hue.high];
      return math::Vector(hs::lerp(a.x, b.x, hue.linear_weight),
                          hs::lerp(a.y, b.y, hue.linear_weight),
                          hs::lerp(a.z, b.z, hue.linear_weight));
    }
    const Pixel &a = view.data[index * HueRotationLutView::HUE_STEPS + hue.low];
    const Pixel &b =
        view.data[index * HueRotationLutView::HUE_STEPS + hue.high];
    return math::Vector(hs::lerp(float(a.r), float(b.r), hue.linear_weight),
                        hs::lerp(float(a.g), float(b.g), hue.linear_weight),
                        hs::lerp(float(a.b), float(b.b), hue.linear_weight));
  };
  const math::Vector a = row(value.low);
  if (value.weight == 0)
    return a;
  const math::Vector b = row(value.high);
  return a + (value.weight * Q16_TO_UNIT) * (b - a);
}

/** @brief Samples hue rotation without intermediate RGB quantization.
 * @param view Prepared hue-rotation LUT.
 * @param value Palette-axis coordinates.
 * @param amount Hue rotation in turns.
 * @return Linear RGB channel values before quantization. */
__attribute__((always_inline)) inline math::Vector
sample_hue_rotation_linear(const HueRotationLutView &view,
                           const HueRotationValueSample &value, float amount) {
  return sample_hue_rotation_linear(view, value, HueRotationHueSample(amount));
}
