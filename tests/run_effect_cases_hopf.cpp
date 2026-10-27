/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Hopf test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/HopfFibration.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/hopf.h"

/** @brief Runs a seeded effect case.
 * @param test Case entry point. */
static void run_case(auto test) {
  hs::random().seed(1337u);
  test();
}

void run_effect_case_hopf_temporal_trail_oracle() {
  HS_RUN(run_case, test_hopf_temporal_trail_oracle);
}

void run_effect_case_hopf_projection_math() {
  HS_RUN(run_case, test_hopf_projection_math);
}

void run_effect_case_hopf_trail_trim_keeps_a_segment() {
  HS_RUN(run_case, test_hopf_trail_trim_keeps_a_segment);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the hopf effect cases.
 * @return Module failure count. */
int run_hopf_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_hopf_temporal_trail_oracle();
  run_effect_case_hopf_projection_math();
  run_effect_case_hopf_trail_trim_keeps_a_segment();
  return fixture.result();
}

} // namespace hs_test::module_runners
