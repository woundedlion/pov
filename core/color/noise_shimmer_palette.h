/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "color/noise_hue_palette.h"

/**
 * @file noise_shimmer_palette.h
 * @brief Spatial noise-driven lightness lift for arbitrary palette sources.
 */

/**
 * @brief Raises palette lightness over the positive lobes of a shared noise field.
 * @tparam Source Palette source exposing Color4 get(float) const.
 */
template <typename Source> class NoiseShimmerPalette {
public:
  NoiseShimmerPalette() = default;

  /**
   * @brief Binds a palette source and prepared cube-map noise field.
   * @param source Non-null palette; must outlive this wrapper.
   * @param noise_lut Non-null prepared noise LUT; must outlive this wrapper.
   */
  NoiseShimmerPalette(const Source *source, const int8_t *noise_lut) {
    bind(source, noise_lut);
  }

  /**
   * @brief Binds a palette source and prepared cube-map noise field.
   * @param source Non-null palette; must outlive this wrapper.
   * @param noise_lut Non-null prepared noise LUT; must outlive this wrapper.
   */
  void bind(const Source *source, const int8_t *noise_lut) {
    HS_CHECK(source != nullptr, "NoiseShimmerPalette bound to null source");
    HS_CHECK(noise_lut != nullptr,
             "NoiseShimmerPalette bound to null noise LUT");
    this->source = source;
    noise_field = {noise_lut, true};
  }

  /**
   * @brief Samples the shared noise field at a finite, nonzero direction.
   * @param direction Sample direction; normalization is not required.
   * @return Noise value in [-1, 1].
   */
  float noise(const math::Vector &direction) const {
    assert(source != nullptr && "NoiseShimmerPalette used before bind()!");
    return sample_hue_noise_lut(noise_field, direction);
  }

  /**
   * @brief Resolves a lightness lift from positive noise and an amount.
   * @param direction Sample direction; normalization is not required.
   * @param amount Lift scale in [0, 1].
   * @return Lift in [0, amount].
   */
  float lightness_shift(const math::Vector &direction, float amount) const {
    return fmaxf(0.0f, noise(direction)) * amount;
  }

  /**
   * @brief Moves OKLab lightness toward white by a resolved lift in [0, 1].
   * @details Preserves alpha and hue; reduces chroma only to stay in gamut.
   * @param value Palette coordinate.
   * @param lift Lightness lift; clamped to [0, 1].
   * @return Lifted palette color.
   */
  Color4 get(float value, float lift) const {
    assert(source != nullptr && "NoiseShimmerPalette used before bind()!");
    return lift_color(source->get(value), lift);
  }

  /**
   * @brief Applies a resolved OKLab lightness lift to a color.
   * @param color Source color; alpha is preserved.
   * @param lift Lightness lift; clamped to [0, 1].
   * @return Lifted color.
   */
  static Color4 lift_color(const Color4 &color, float lift) {
    lift = hs::clamp(lift, 0.0f, 1.0f);
    if (lift == 0.0f)
      return color;
    const LinRGB input = pixel_to_linrgb(color.color);
    OKLab lab = linear_rgb_to_oklab_fast(input.r, input.g, input.b);
    lab.L += (1.0f - lab.L) * lift;
    const LinRGB output = oklab_to_linear_rgb_lut_gamut(lab);
    return Color4(linrgb_to_pixel(output), color.alpha);
  }

  /**
   * @brief Samples the palette with a sphere-domain noise lightness lift.
   * @param value Palette coordinate.
   * @param direction Noise sample direction; normalization is not required.
   * @param amount Lift scale in [0, 1].
   * @return Lifted palette color.
   */
  Color4 get(float value, const math::Vector &direction, float amount) const {
    return get(value, lightness_shift(direction, amount));
  }

private:
  const Source *source = nullptr;
  HueNoiseLutView noise_field{nullptr, false};
};
