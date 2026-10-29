/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/pattern/types.h"
#include <cstdint>

/** @file modulation_types.h
 * @brief Passive pattern modulation recipes and mechanic controls. */

namespace ColorModulation {
/** @brief Operation applied to a palette color. */
enum class Mechanic : uint8_t {
  HUE_ROTATION,            ///< Signed hue rotation.
  POSITIVE_LIGHTNESS_LIFT, ///< Positive-lobe OKLab lightness lift.
  OPACITY                  ///< Straight-alpha mask.
};
/** @brief Controls of one selected mechanic. */
struct Params {
  Mechanic kind = Mechanic::HUE_ROTATION; ///< Selected operation.
  float amount = 0.0f; ///< Rotation turns or unit-range strength.
  constexpr bool operator==(const Params &) const = default;
};
/** @brief Serializable sphere modulation recipe. */
struct Config {
  bool enabled = true; ///< Whether modulation is applied.
  Pattern::SourceKind source = Pattern::SourceKind::NOISE; ///< Source recipe.
  uint16_t version = 1; ///< Source recipe version.
  Pattern::Evaluation evaluation =
      Pattern::Evaluation::BAKED_CUBE;      ///< Evaluation strategy.
  Pattern::NoiseParams noise;               ///< Noise source controls.
  Pattern::SineParams sine;                 ///< Directional sine controls.
  Pattern::InterferenceParams interference; ///< Interference controls.
  Pattern::GridParams grid;                 ///< Smooth grid controls.
  int32_t seed = 1337;                      ///< Noise generator seed.
  Params mechanic;                          ///< Color operation.
  constexpr bool operator==(const Config &) const = default;
};
} // namespace ColorModulation
