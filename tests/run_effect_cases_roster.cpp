/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Roster test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "targets/effects.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/roster_cases.h"
#include "tests/effects/column_time_parity.h"
#include "tests/effects/lens_cases.h"

void run_effect_case_parameter_probe_targets() {
  HS_RUN(run_case, test_parameter_probe_targets);
}

void run_effect_case_float_dropdown_lint() {
  HS_RUN(run_case, test_float_dropdown_lint);
}

void run_effect_case_manual_preset_navigation() {
  HS_RUN(run_case, test_manual_preset_navigation);
}

void run_effect_case_needs_full_frame_gate() {
  HS_RUN(run_case, test_needs_full_frame_gate);
}

void run_effect_case_column_time_scans_keep_fidelity_off_frames() {
  HS_RUN(run_case, test_column_time_scans_keep_fidelity_off_frames);
}

void run_effect_case_column_time_scans_engage_with_fidelity_on() {
  HS_RUN(run_case, test_column_time_scans_engage_with_fidelity_on);
}

void run_effect_case_glitch_lens_unit_norm() {
  HS_RUN(run_case, test_glitch_lens_unit_norm);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the roster effect cases.
 * @return Module failure count. */
int run_roster_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_parameter_probe_targets();
  run_effect_case_float_dropdown_lint();
  run_effect_case_manual_preset_navigation();
  run_effect_case_needs_full_frame_gate();
#if defined(_WIN32)
  run_effect_case_column_time_scans_keep_fidelity_off_frames();
#endif
  run_effect_case_column_time_scans_engage_with_fidelity_on();
  run_effect_case_glitch_lens_unit_norm();
  return fixture.result();
}

} // namespace hs_test::module_runners
