/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Ring spin test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/RingSpin.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/ring_spin.h"

void run_effect_case_ringspin_strobe_configuration_preserves_rendering() {
  HS_RUN(run_case, test_ringspin_strobe_configuration_preserves_rendering);
}

void run_effect_case_ringspin_trail_hugs_its_great_circles() {
  HS_RUN(run_case, test_ringspin_trail_hugs_its_great_circles);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the ring spin effect cases.
 * @return Module failure count. */
int run_ring_spin_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_ringspin_strobe_configuration_preserves_rendering();
  run_effect_case_ringspin_trail_hugs_its_great_circles();
  return fixture.result();
}

} // namespace hs_test::module_runners
