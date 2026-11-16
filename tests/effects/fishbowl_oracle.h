/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Fishbowl temporal oracle. */
#pragma once

#include "effects/Fishbowl.h"
#include "tests/effects/fixture.h"
#include "tests/effects/fishbowl_cases.h"

namespace hs_test::effects_tests {

/** @brief Warped Fishbowl endpoints, births and palette match frozen columns. */
inline void test_fishbowl_temporal_trail_oracle() {
  const ScopedTemporalSampling temporal_sampling(true);
  for (bool full : {false, true}) {
    HS_CONTEXT("full", full);
    auto render = [full](bool timed, float fraction, float displacement) {
      return FishbowlWhiteBox::render_temporal_trail(timed, fraction,
                                                     displacement, full);
    };
    hs_test::seam_compare_render<DEFAULT_W, DEFAULT_H>(render,
                                                       "Fishbowl trail");
  }
}

} // namespace hs_test::effects_tests
