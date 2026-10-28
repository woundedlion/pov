/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Comets test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/Comets.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/comets_cases.h"

void run_effect_case_comets_white_box_check_paths_close() {
  HS_RUN(CometsWhiteBox::check_paths_close);
}

void run_effect_case_comets_white_box_check_births_skip_earlier_columns() {
  HS_RUN(CometsWhiteBox::check_births_skip_earlier_columns);
}

void run_effect_case_comets_rollover_skipped_mid_wipe() {
  HS_RUN(run_case, test_comets_rollover_skipped_mid_wipe);
}

void run_effect_case_comets_wipe_shades_at_column_time() {
  HS_RUN(run_case, test_comets_wipe_shades_at_column_time);
}

void run_effect_case_comets_wipe_history() {
  HS_RUN(run_case, test_comets_wipe_history);
}

void run_effect_case_comets_manual_preset_restarts_path() {
  HS_RUN(run_case, test_comets_manual_preset_restarts_path);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the comets effect cases.
 * @return Module failure count. */
int run_comets_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_comets_white_box_check_paths_close();
  run_effect_case_comets_white_box_check_births_skip_earlier_columns();
  run_effect_case_comets_rollover_skipped_mid_wipe();
  run_effect_case_comets_wipe_shades_at_column_time();
  run_effect_case_comets_wipe_history();
  if (effects_full_suite()) {
    run_effect_case_comets_manual_preset_restarts_path();
  }
  return fixture.result();
}

} // namespace hs_test::module_runners
