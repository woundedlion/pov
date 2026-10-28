/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Dreamballs test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/DreamBalls.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/dreamballs_cases.h"

void run_effect_case_dreamballs_temporal_source_oracle() {
  HS_RUN(run_case, test_dreamballs_temporal_source_oracle);
}

void run_effect_case_dreamballs_column_work_is_bounded() {
  HS_RUN(run_case, test_dreamballs_column_work_is_bounded);
}

void run_effect_case_dream_balls_white_box_check_endpoint_error_bound() {
  HS_RUN(run_case, DreamBallsWhiteBox::check_endpoint_error_bound);
}

void run_effect_case_dreamballs_max_edge_solid_render() {
  HS_RUN(run_case, test_dreamballs_max_edge_solid_render);
}

void run_effect_case_dreamballs_retained_history_resets() {
  HS_RUN(run_case, test_dreamballs_retained_history_resets);
}

void run_effect_case_dreamballs_preset_cycle_bookkeeping() {
  HS_RUN(run_case, test_dreamballs_preset_cycle_bookkeeping);
}

void run_effect_case_dreamballs_base_mesh_selector() {
  HS_RUN(run_case, test_dreamballs_base_mesh_selector);
}

void run_effect_case_dreamballs_weave_topology() {
  HS_RUN(run_case, test_dreamballs_weave_topology);
}

void run_effect_case_dreamballs_defect_weave_renders() {
  HS_RUN(run_case, test_dreamballs_defect_weave_renders);
}

void run_effect_case_dreamballs_respawn_fires_and_honors_pause() {
  HS_RUN(run_case, test_dreamballs_respawn_fires_and_honors_pause);
}

void run_effect_case_dreamballs_select_mid_fade_in_finishes_fade() {
  HS_RUN(run_case, test_dreamballs_select_mid_fade_in_finishes_fade);
}

void run_effect_case_dreamballs_select_mid_fade_out_respawns() {
  HS_RUN(run_case, test_dreamballs_select_mid_fade_out_respawns);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the dreamballs effect cases.
 * @return Module failure count. */
int run_dreamballs_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_dreamballs_temporal_source_oracle();
  run_effect_case_dreamballs_column_work_is_bounded();
  run_effect_case_dream_balls_white_box_check_endpoint_error_bound();
  run_effect_case_dreamballs_max_edge_solid_render();
  run_effect_case_dreamballs_retained_history_resets();
  run_effect_case_dreamballs_preset_cycle_bookkeeping();
  if (effects_full_suite()) {
    run_effect_case_dreamballs_base_mesh_selector();
    run_effect_case_dreamballs_weave_topology();
    run_effect_case_dreamballs_defect_weave_renders();
    run_effect_case_dreamballs_respawn_fires_and_honors_pause();
    run_effect_case_dreamballs_select_mid_fade_in_finishes_fade();
    run_effect_case_dreamballs_select_mid_fade_out_respawns();
  }
  return fixture.result();
}

} // namespace hs_test::module_runners
