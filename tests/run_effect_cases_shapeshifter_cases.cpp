/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Shapeshifter cases test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/ShapeShifter.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/shape_shifter.h"

void run_effect_case_shapeshifter_preset_defaults() {
  HS_RUN(run_case, test_shapeshifter_preset_defaults);
}

void run_effect_case_shapeshifter_slider_selections_render() {
  HS_RUN(run_case, test_shapeshifter_slider_selections_render);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the shapeshifter cases effect cases.
 * @return Module failure count. */
int run_shapeshifter_cases_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_shapeshifter_preset_defaults();
  run_effect_case_shapeshifter_slider_selections_render();
  return fixture.result();
}

} // namespace hs_test::module_runners
