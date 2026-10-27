/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Mobiusrings test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/MobiusRings.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/mobius_rings.h"

/** @brief Runs a seeded effect case.
 * @param test Case entry point. */
static void run_case(auto test) {
  hs::random().seed(1337u);
  test();
}

void run_effect_case_mobius_rings_white_box_check_column_time_oracle() {
  HS_RUN(MobiusRingsWhiteBox::check_column_time_oracle);
}

void run_effect_case_mobius_rings_conformal_and_counter_rotation() {
  HS_RUN(run_case, test_mobius_rings_conformal_and_counter_rotation);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the mobiusrings effect cases.
 * @return Module failure count. */
int run_mobiusrings_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_mobius_rings_conformal_and_counter_rotation();
  run_effect_case_mobius_rings_white_box_check_column_time_oracle();
  return fixture.result();
}

} // namespace hs_test::module_runners
