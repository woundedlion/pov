/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Shared authored feedback presets.
 */
#pragma once

#include "core/render/filter/feedback_style.h"

namespace hs_test::styles_tests {

/** @brief Named authored feedback style. */
struct PresetCase {
  const char *name;      /**< Preset name. */
  Feedback::Style style; /**< Preset parameters. */
};
/** @brief Authored styles exercised by preset and warp tests. */
inline const PresetCase PRESETS[] = {
    {"ArcingLightning", Feedback::Style::ArcingLightning()},
    {"SlowFire", Feedback::Style::SlowFire()},
    {"EnergeticFire", Feedback::Style::EnergeticFire()},
    {"SlowDust", Feedback::Style::SlowDust()},
    {"WavyTrails", Feedback::Style::WavyTrails()},
    {"MeltingHi", Feedback::Style::MeltingHi()},
    {"MeltingLo", Feedback::Style::MeltingLo()},
    {"Miasma", Feedback::Style::Miasma()},
    {"LooseWormhole", Feedback::Style::LooseWormhole()},
    {"TightWormhole", Feedback::Style::TightWormhole()},
    {"WigglingWormhole", Feedback::Style::WigglingWormhole()},
    {"Smoke", Feedback::Style::Smoke()},
};

} // namespace hs_test::styles_tests
