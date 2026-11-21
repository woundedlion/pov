/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Roster test entry points. */
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

#include "tests/effects/roster_cases.h"
#include "tests/effects/four_stroke_stars.h"
#include "tests/effects/column_time_fidelity.h"
#include "tests/effects/lens_cases.h"

namespace hs_test::effects_tests {

void run_effect_case_mesh_seam_dark_sample_guards() {
  HS_RUN(run_case, test_mesh_seam_dark_sample_guards);
}

void run_effect_case_parameter_lint_failure_guards() {
  HS_RUN(run_case, test_parameter_lint_failure_guards);
}

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

void run_effect_case_column_time_scans_engage_with_fidelity_on() {
  HS_RUN(run_case, test_column_time_scans_engage_with_fidelity_on);
}

void run_effect_case_glitch_lens_doubles_latitude_and_stays_smooth_at_fold() {
  HS_RUN(run_case, test_glitch_lens_doubles_latitude_and_stays_smooth_at_fold);
}

void run_effect_case_configured_glitch_general_path() {
  HS_RUN(run_case, test_configured_glitch_general_path);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the roster effect cases.
 * @return Module failure count. */
int run_roster_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  HS_RUN(run_case, test_four_stroke_stars_fixed_geometry);
  HS_RUN(run_case, test_four_stroke_stars_skipped_window);
  run_effect_case_mesh_seam_dark_sample_guards();
  run_effect_case_parameter_lint_failure_guards();
  run_effect_case_parameter_probe_targets();
  run_effect_case_float_dropdown_lint();
  run_effect_case_manual_preset_navigation();
  run_effect_case_needs_full_frame_gate();
  run_effect_case_column_time_scans_engage_with_fidelity_on();
  run_effect_case_glitch_lens_doubles_latitude_and_stays_smooth_at_fold();
  run_effect_case_configured_glitch_general_path();
  return fixture.result();
}

} // namespace hs_test::module_runners
