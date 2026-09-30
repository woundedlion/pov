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

  /** @brief Binds a palette source and prepared cube-map noise field. */
  NoiseShimmerPalette(const Source *source, const int8_t *noise_lut) {
    bind(source, noise_lut);
  }

  /** @brief Binds a palette source and prepared cube-map noise field. */
  void bind(const Source *source, const int8_t *noise_lut) {
    HS_CHECK(source != nullptr, "NoiseShimmerPalette bound to null source");
    HS_CHECK(noise_lut != nullptr,
             "NoiseShimmerPalette bound to null noise LUT");
    this->source = source;
    noise_field = {noise_lut, true};
  }

  /** @brief Samples the shared noise field at a finite, nonzero direction. */
  float noise(const math::Vector &direction) const {
    assert(source != nullptr && "NoiseShimmerPalette used before bind()!");
    return sample_hue_noise_lut(noise_field, direction);
  }

  /** @brief Resolves a lightness lift from positive noise and an amount in [0, 1]. */
  float lightness_shift(const math::Vector &direction, float amount) const {
    return std::max(0.0f, noise(direction)) * amount;
  }

  /**
   * @brief Moves OKLab lightness toward white by a resolved lift in [0, 1].
   * @details Preserves alpha and hue; reduces chroma only to stay in gamut.
   */
  Color4 get(float value, float lift) const {
    assert(source != nullptr && "NoiseShimmerPalette used before bind()!");
    const Color4 color = source->get(value);
    lift = hs::clamp(lift, 0.0f, 1.0f);
    if (lift == 0.0f)
      return color;
    const LinRGB input = pixel_to_linrgb(color.color);
    OKLab lab = linear_rgb_to_oklab_fast(input.r, input.g, input.b);
    lab.L += (1.0f - lab.L) * lift;
    LinRGB output = oklab_to_linear_rgb(lab);
    if (!linear_rgb_in_gamut(output.r, output.g, output.b))
      output = oklab_to_linear_rgb(gamut_scale_to_boundary_lut(lab));
    return Color4(linrgb_to_pixel(output), color.alpha);
  }

  /** @brief Samples the palette with a sphere-domain noise lightness lift. */
  Color4 get(float value, const math::Vector &direction, float amount) const {
    return get(value, lightness_shift(direction, amount));
  }

private:
  const Source *source = nullptr;
  HueNoiseLutView noise_field{nullptr, false};
};
