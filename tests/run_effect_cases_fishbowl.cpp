/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Fishbowl test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/Fishbowl.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/fishbowl_cases.h"
#include "tests/effects/fishbowl_oracle.h"

void run_effect_case_fishbowl_snapshot_admission() {
  HS_RUN(run_case, test_fishbowl_snapshot_admission);
}

void run_effect_case_fishbowl_preset_restores_params() {
  HS_RUN(run_case, test_fishbowl_preset_restores_params);
}

void run_effect_case_fishbowl_fire_duty_cycle() {
  HS_RUN(run_case, test_fishbowl_fire_duty_cycle);
}

void run_effect_case_fishbowl_palette_phase_compensates_fill() {
  HS_RUN(run_case, test_fishbowl_palette_phase_compensates_fill);
}

void run_effect_case_fishbowl_adaptive_midpoint() {
  HS_RUN(run_case, test_fishbowl_adaptive_midpoint);
}
void run_effect_case_fishbowl_temporal_trail_oracle() {
  HS_RUN(run_case, test_fishbowl_temporal_trail_oracle);
}

void run_effect_case_fishbowl_scratch_estimate_covers_peak() {
  HS_RUN(run_case, test_fishbowl_scratch_estimate_covers_peak);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the fishbowl effect cases.
 * @return Module failure count. */
int run_fishbowl_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_fishbowl_snapshot_admission();
  run_effect_case_fishbowl_preset_restores_params();
  run_effect_case_fishbowl_fire_duty_cycle();
  run_effect_case_fishbowl_palette_phase_compensates_fill();
  run_effect_case_fishbowl_adaptive_midpoint();
  run_effect_case_fishbowl_temporal_trail_oracle();
  run_effect_case_fishbowl_scratch_estimate_covers_peak();
  return fixture.result();
}

} // namespace hs_test::module_runners
