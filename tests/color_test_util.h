/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Shared hue-arc predicate for the suites that assert on OKLCH hue — color and
 * palettes. Kept independent of core/color/color.h's wrap_angle_pi so a hue
 * assertion never compares the module under test against itself.
 */
#pragma once

#include "core/math/3dmath.h"

namespace hs_test {

/**
 * @brief Wraps a hue difference into (-PI, PI] for circular comparison.
 * @param dh Raw hue difference in radians; magnitude below a few turns.
 * @return The equivalent difference in (-PI, PI].
 */
inline float wrap_hue_delta(float dh) {
  while (dh > PI_F)
    dh -= 2.0f * PI_F;
  while (dh < -PI_F)
    dh += 2.0f * PI_F;
  return dh;
}

} // namespace hs_test
