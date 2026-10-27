/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Mesh feedback test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/MeshFeedback.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/mesh_feedback.h"

/** @brief Runs a seeded effect case.
 * @param test Case entry point. */
static void run_case(auto test) {
  hs::random().seed(1337u);
  test();
}

void run_effect_case_meshfeedback_temporal_source_oracle() {
  HS_RUN(run_case, test_meshfeedback_temporal_source_oracle);
}

void run_effect_case_meshfeedback_stationary_timed_source() {
  HS_RUN(run_case, test_meshfeedback_stationary_timed_source);
}

void run_effect_case_meshfeedback_trails_arm_cut_seam() {
  HS_RUN(run_case, test_meshfeedback_trails_arm_cut_seam);
}

void run_effect_case_meshfeedback_base_mesh_selector() {
  HS_RUN(run_case, test_meshfeedback_base_mesh_selector);
}

void run_effect_case_meshfeedback_preset_export_arity() {
  HS_RUN(run_case, test_meshfeedback_preset_export_arity);
}

void run_effect_case_meshfeedback_mesh_rebuild_reuses_storage() {
  HS_RUN(run_case, test_meshfeedback_mesh_rebuild_reuses_storage);
}

void run_effect_case_meshfeedback_flush_precedes_mesh_draw() {
  HS_RUN(run_case, test_meshfeedback_flush_precedes_mesh_draw);
}

void run_effect_case_meshfeedback_preset_rotation_syncs_noise() {
  HS_RUN(run_case, test_meshfeedback_preset_rotation_syncs_noise);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the mesh feedback effect cases.
 * @return Module failure count. */
int run_mesh_feedback_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_meshfeedback_temporal_source_oracle();
  run_effect_case_meshfeedback_stationary_timed_source();
  run_effect_case_meshfeedback_trails_arm_cut_seam();
  run_effect_case_meshfeedback_base_mesh_selector();
  run_effect_case_meshfeedback_preset_export_arity();
  run_effect_case_meshfeedback_mesh_rebuild_reuses_storage();
  if (effects_full_suite()) {
    run_effect_case_meshfeedback_flush_precedes_mesh_draw();
    run_effect_case_meshfeedback_preset_rotation_syncs_noise();
  }
  return fixture.result();
}

} // namespace hs_test::module_runners
