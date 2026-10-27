/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Petal flow test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/PetalFlow.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/petal_flow.h"
#include "tests/effects/petal_flow_oracle.h"

/** @brief Runs a seeded effect case.
 * @param test Case entry point. */
static void run_case(auto test) {
  hs::random().seed(1337u);
  test();
}

void run_effect_case_petalflow_temporal_source_oracle() {
  HS_RUN(run_case, test_petalflow_temporal_source_oracle);
}

void run_effect_case_petalflow_spawn_gap_bounded() {
  HS_RUN(run_case, test_petalflow_spawn_gap_bounded);
}

void run_effect_case_petalflow_last_free_slot() {
  HS_RUN(run_case, test_petalflow_last_free_slot);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the petal flow effect cases.
 * @return Module failure count. */
int run_petal_flow_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_petalflow_temporal_source_oracle();
  if (effects_full_suite()) {
    run_effect_case_petalflow_spawn_gap_bounded();
    run_effect_case_petalflow_last_free_slot();
  }
  return fixture.result();
}

} // namespace hs_test::module_runners
