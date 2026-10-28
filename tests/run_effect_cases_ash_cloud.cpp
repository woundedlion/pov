/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Ash cloud test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/AshCloud.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/ash_cloud.h"

void run_effect_case_ash_cloud_value_cutout_gates_the_frame() {
  HS_RUN(run_case, test_ash_cloud_value_cutout_gates_the_frame);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the ash cloud effect cases.
 * @return Module failure count. */
int run_ash_cloud_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  if (effects_full_suite()) {
    run_effect_case_ash_cloud_value_cutout_gates_the_frame();
  }
  return fixture.result();
}

} // namespace hs_test::module_runners
