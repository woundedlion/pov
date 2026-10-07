/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file preset_departure.h
 * @brief Automatic preset departure policies. */

#include <cstdint>
#include <variant>
#include "engine/concepts.h"
#include "math/easing.h"

namespace Segue {
/**
 * @brief Preset-transition policies: how ChoreographedEffect carries an
 * AUTOMATIC preset change onto its target parameter set.
 * @details Each preset names the policy it departs by (PresetEntry::segue),
 * whichever preset comes next. MANUAL and SYNCHRONIZED changes always snap.
 */
namespace Preset {

/** @brief Departure: adopt the next preset immediately; no transition state. */
struct Snap {};

/**
 * @brief Departure: parameter-space crossfade. ChoreographedEffect drives
 * Derived::blend_params(progress) from the departing to the next set.
 */
struct Lerp {
  uint16_t frames = 0; /**< Frames the parameter crossfade spans. */
  EasingFn easing =
      math::ease_linear; /**< Easing applied to the blend progress. */
  bool pausable =
      false; /**< Whether anims_paused freezes an in-flight blend. */
};

/**
 * @brief Departure: dims the departing preset, adopts the next at the first
 * progress sample >= 0.5, then brightens it back to full.
 * @details Even positive frame counts sample the zero-opacity midpoint; odd
 * counts skip it. ChoreographedEffect feeds opacity to
 * Derived::set_preset_opacity; the fade holds full opacity while paused.
 */
struct Fade {
  uint16_t frames =
      0; /**< Dim/brighten duration; zero selects an immediate snap. */
};

/** @brief How a preset departs: the policy of the automatic transition that
 * leaves it. */
using Departure = std::variant<Snap, Lerp, Fade>;

/** @brief Frames a departure spans; zero for Snap. */
constexpr uint16_t frames(const Departure &departure) {
  if (const auto *lerp = std::get_if<Lerp>(&departure))
    return lerp->frames;
  if (const auto *fade = std::get_if<Fade>(&departure))
    return fade->frames;
  return 0;
}

} // namespace Preset
} // namespace Segue
