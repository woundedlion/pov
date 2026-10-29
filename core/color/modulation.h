/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/color/color_space.h"
#include "core/color/modulation_types.h"
#include "core/color/hue_rotation_lut.h"
#include <cmath>

/** @file modulation.h
 * @brief Domain-independent color mechanics and generic palette composition. */

namespace ColorModulation {
/** @brief Strength of a statically selected mechanic. */
struct AmountParams {
  float amount = 0.0f; ///< Rotation turns or unit-range strength.
};
/** @brief Tests finite unit-range strength.
 * @param amount Mechanic strength.
 * @return Whether strength is finite and in [0, 1]. */
inline bool unit_amount(float amount) {
  return Pattern::finite(amount) && amount >= 0.0f && amount <= 1.0f;
}
/** @brief Tests whether a mechanic requires unit-range strength.
 * @param mechanic Selected mechanic.
 * @return Whether strength is restricted to [0, 1]. */
inline constexpr bool unit_mechanic(Mechanic mechanic) {
  return mechanic == Mechanic::OPACITY ||
         mechanic == Mechanic::POSITIVE_LIGHTNESS_LIFT;
}
/** @brief Validates a finite signed rotation amount.
 * @param params Rotation strength.
 * @return Validated controls. */
inline AmountParams prepare_hue(const AmountParams &params) {
  HS_CHECK(Pattern::finite(params.amount), "HueRotation amount must be finite");
  return params;
}
/** @brief Validates a unit-range mechanic amount.
 * @param params Mechanic strength.
 * @return Validated controls. */
inline AmountParams prepare_unit(const AmountParams &params) {
  HS_CHECK(unit_amount(params.amount),
           "modulation amount must be finite and in [0, 1]");
  return params;
}
/** @brief Signed hue rotation preserving straight alpha. */
struct HueRotation {
  using Params = AmountParams;   ///< User controls.
  using Prepared = AmountParams; ///< Validated controls.
  /** @brief Tests whether evaluation can return the exact source color.
   * @param prepared Validated mechanic controls.
   * @return Whether the mechanic is identity for every scalar. */
  static bool is_identity(const Prepared &prepared) {
    return prepared.amount == 0.0f;
  }
  /** @brief Prepares rotation strength.
   * @param params User controls.
   * @return Validated controls. */
  static Prepared prepare(const Params &params) { return prepare_hue(params); }
  /** @brief Applies a periodic hue rotation.
   * @param color Source color.
   * @param scalar Signed pattern sample.
   * @param prepared Validated strength.
   * @return Rotated color with original alpha. */
  static Color4 apply(Color4 color, float scalar, const Prepared &prepared) {
    const float rotation = prepared.amount * scalar;
    if (rotation == 0.0f || rotation - floorf(rotation) == 0.0f)
      return color;
    return hue_rotate_lut_gamut(color, rotation);
  }
};
/** @brief Positive-lobe OKLab lift preserving straight alpha. */
struct PositiveLightnessLift {
  using Params = AmountParams;   ///< User controls.
  using Prepared = AmountParams; ///< Validated controls.
  /** @brief Tests whether evaluation can return the exact source color.
   * @param prepared Validated mechanic controls.
   * @return Whether the mechanic is identity for every scalar. */
  static bool is_identity(const Prepared &prepared) {
    return prepared.amount == 0.0f;
  }
  /** @brief Prepares lift strength.
   * @param params User controls.
   * @return Validated controls. */
  static Prepared prepare(const Params &params) { return prepare_unit(params); }
  /** @brief Raises lightness toward white.
   * @param color Source color.
   * @param scalar Signed pattern sample.
   * @param prepared Validated strength.
   * @return Lifted color with original alpha. */
  static Color4 apply(Color4 color, float scalar, const Prepared &prepared) {
    const float lift = prepared.amount * fmaxf(0.0f, scalar);
    if (lift == 0.0f)
      return color;
    const LinRGB input = pixel_to_linrgb(color.color);
    OKLab lab = linear_rgb_to_oklab_fast(input.r, input.g, input.b);
    lab.L += (1.0f - lab.L) * lift;
    return Color4(linrgb_to_pixel(oklab_to_linear_rgb_lut_gamut(lab)),
                  color.alpha);
  }
};
/** @brief Pattern mask applied only to alpha. */
struct Opacity {
  using Params = AmountParams;   ///< User controls.
  using Prepared = AmountParams; ///< Validated controls.
  /** @brief Tests whether evaluation can return the exact source color.
   * @param prepared Validated mechanic controls.
   * @return Whether the mechanic is identity for every scalar. */
  static bool is_identity(const Prepared &prepared) {
    return prepared.amount == 0.0f;
  }
  /** @brief Prepares mask strength.
   * @param params User controls.
   * @return Validated controls. */
  static Prepared prepare(const Params &params) { return prepare_unit(params); }
  /** @brief Multiplies straight alpha by a signed pattern mask.
   * @param color Source color.
   * @param scalar Signed pattern sample.
   * @param prepared Validated strength.
   * @return Unchanged RGB and masked alpha. */
  static Color4 apply(Color4 color, float scalar, const Prepared &prepared) {
    if (prepared.amount == 0.0f)
      return color;
    color.alpha *= hs::lerp(1.0f, 0.5f * (scalar + 1.0f), prepared.amount);
    return color;
  }
};
/** @brief Tests selected mechanic parameter admission.
 * @param params Mechanic controls.
 * @return Whether the kind and amount are admitted. */
inline bool valid(const Params &params) {
  return unit_mechanic(params.kind) ? unit_amount(params.amount)
                                    : params.kind == Mechanic::HUE_ROTATION &&
                                          Pattern::finite(params.amount);
}
/** @brief Validates a selected mechanic.
 * @param params Selected mechanic controls. */
inline void validate(const Params &params) {
  switch (params.kind) {
  case Mechanic::HUE_ROTATION:
    (void)HueRotation::prepare({params.amount});
    return;
  case Mechanic::POSITIVE_LIGHTNESS_LIFT:
  case Mechanic::OPACITY:
    (void)prepare_unit({params.amount});
    return;
  }
  HS_CHECK(false, "unknown modulation preparation mechanic");
}
/** @brief Applies a prepared selected mechanic.
 * @param color Source color.
 * @param scalar Signed pattern sample.
 * @param params Validated mechanic controls.
 * @return Modulated straight-alpha color. */
inline Color4 apply(Color4 color, float scalar, const Params &params) {
  switch (params.kind) {
  case Mechanic::HUE_ROTATION:
    return HueRotation::apply(color, scalar, {params.amount});
  case Mechanic::POSITIVE_LIGHTNESS_LIFT:
    return PositiveLightnessLift::apply(color, scalar, {params.amount});
  case Mechanic::OPACITY:
    return Opacity::apply(color, scalar, {params.amount});
  }
  HS_CHECK(false, "unknown modulation application mechanic");
  return color;
}
/** @brief Composes any palette, typed scalar field, and color mechanic.
 * @tparam PaletteSource Palette exposing Color4 get(float) const.
 * @tparam Field Immutable typed scalar field.
 * @tparam Mechanic Domain-independent color operation. */
template <class PaletteSource, class Field, class Mechanic>
struct PatternModulatedPalette {
  const PaletteSource *palette = nullptr; ///< Borrowed base palette.
  typename Field::Prepared field;         ///< Prepared scalar field.
  typename Mechanic::Prepared mechanic;   ///< Prepared color operation.
  /** @brief Samples a modulated palette.
   * @param value Palette lookup coordinate.
   * @param coordinate Typed field coordinate.
   * @return Modulated straight-alpha color. */
  Color4 get(float value, const typename Field::Coordinate &coordinate) const {
    assert(palette != nullptr);
    const Color4 color = palette->get(value);
    if constexpr (requires { Mechanic::is_identity(mechanic); }) {
      if (Mechanic::is_identity(mechanic))
        return color;
    }
    return Mechanic::apply(color, Field::sample(coordinate, field), mechanic);
  }
};

} // namespace ColorModulation
