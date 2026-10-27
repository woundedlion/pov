/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Gnomonic stars test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/GnomonicStars.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/gnomonic_stars.h"

/** @brief Runs a seeded effect case.
 * @param test Case entry point. */
static void run_case(auto test) {
  hs::random().seed(1337u);
  test();
}

void run_effect_case_gnomonicstars_radius_px_covers_both_axes() {
  HS_RUN(run_case, test_gnomonicstars_radius_px_covers_both_axes);
}

void run_effect_case_gnomonicstars_spiral_cache_invalidation() {
  HS_RUN(run_case, test_gnomonicstars_spiral_cache_invalidation);
}

void run_effect_case_gnomonicstars_stars_drift_from_last_draw() {
  HS_RUN(run_case, test_gnomonicstars_stars_drift_from_last_draw);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the gnomonic stars effect cases.
 * @return Module failure count. */
int run_gnomonic_stars_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_gnomonicstars_radius_px_covers_both_axes();
  run_effect_case_gnomonicstars_spiral_cache_invalidation();
  run_effect_case_gnomonicstars_stars_drift_from_last_draw();
  return fixture.result();
}

} // namespace hs_test::module_runners
