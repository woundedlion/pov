/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Displacement field test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/DisplacementField.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/displacement_field.h"

/** @brief Runs a seeded effect case.
 * @param test Case entry point. */
static void run_case(auto test) {
  hs::random().seed(1337u);
  test();
}

void run_effect_case_displacement_field_retained_history_resets() {
  HS_RUN(run_case, test_displacement_field_retained_history_resets);
}

void run_effect_case_displacement_field_lazy_hue_table_matches_eager() {
  HS_RUN(run_case, test_displacement_field_lazy_hue_table_matches_eager);
}

void run_effect_case_displacement_field_hue_table_bind_clears_knots() {
  HS_RUN(run_case, test_displacement_field_hue_table_bind_clears_knots);
}

void run_effect_case_displacement_field_ring_pool_abandoned_slot() {
  HS_RUN(run_case, test_displacement_field_ring_pool_abandoned_slot);
}

void run_effect_case_displacement_field_ring_pool_capacity_release() {
  HS_RUN(run_case, test_displacement_field_ring_pool_capacity_release);
}

void run_effect_case_displacement_field_footprint_covers_init() {
  HS_RUN(run_case, test_displacement_field_footprint_covers_init);
}

void run_effect_case_displacement_field_pause_holds_phase_only() {
  HS_RUN(run_case, test_displacement_field_pause_holds_phase_only);
}

void run_effect_case_displacement_field_default_cycle_frames() {
  HS_RUN(run_case, test_displacement_field_default_cycle_frames);
}

void run_effect_case_displacement_field_presets() {
  HS_RUN(run_case, test_displacement_field_presets);
}

void run_effect_case_displacement_field_edit_cancels_departure() {
  HS_RUN(run_case, test_displacement_field_edit_cancels_departure);
}

void run_effect_case_displacement_field_startup_writes_preserve_entry() {
  HS_RUN(run_case, test_displacement_field_startup_writes_preserve_entry);
}

void run_effect_case_displacement_field_paused_balls_deform_render() {
  HS_RUN(run_case, test_displacement_field_paused_balls_deform_render);
}

void run_effect_case_displacement_field_paused_balls_keep_moving() {
  HS_RUN(run_case, test_displacement_field_paused_balls_keep_moving);
}

void run_effect_case_displacement_field_knot_walk_unit_near_ball() {
  HS_RUN(run_case, test_displacement_field_knot_walk_unit_near_ball);
}

void run_effect_case_displacement_field_zero_hue_scale_is_exact() {
  HS_RUN(run_case, test_displacement_field_zero_hue_scale_is_exact);
}

void run_effect_case_displacement_field_octave_bake_tracks_noise() {
  HS_RUN(run_case, test_displacement_field_octave_bake_tracks_noise);
}

void run_effect_case_displacement_field_arm_boundaries_meet_endpoints() {
  HS_RUN(run_case, test_displacement_field_arm_boundaries_meet_endpoints);
}

void run_effect_case_displacement_field_single_pass_reference() {
  HS_RUN(run_case, test_displacement_field_single_pass_reference);
}

void run_effect_case_displacement_field_boundary_stencil() {
  HS_RUN(run_case, test_displacement_field_boundary_stencil);
}

void run_effect_case_displacement_field_turn_follows_large_steps() {
  HS_RUN(run_case, test_displacement_field_turn_follows_large_steps);
}

void run_effect_case_displacement_field_alternate_groups_match_windows() {
  HS_RUN(run_case, test_displacement_field_alternate_groups_match_windows);
}

void run_effect_case_displacement_field_hue_table_fidelity() {
  HS_RUN(run_case, test_displacement_field_hue_table_fidelity);
}

void run_effect_case_displacement_field_hue_table_frame_fidelity() {
  HS_RUN(run_case, test_displacement_field_hue_table_frame_fidelity);
}

void run_effect_case_displacement_field_clip_tiles_full() {
  HS_RUN(run_case, test_displacement_field_clip_tiles_full);
}

void run_effect_case_displacement_field_ball_spans_and_lifecycle() {
  HS_RUN(run_case, test_displacement_field_ball_spans_and_lifecycle);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the displacement field effect cases.
 * @return Module failure count. */
int run_displacement_field_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_displacement_field_retained_history_resets();
  run_effect_case_displacement_field_lazy_hue_table_matches_eager();
  run_effect_case_displacement_field_hue_table_bind_clears_knots();
  run_effect_case_displacement_field_ring_pool_abandoned_slot();
  run_effect_case_displacement_field_ring_pool_capacity_release();
  run_effect_case_displacement_field_footprint_covers_init();
  run_effect_case_displacement_field_pause_holds_phase_only();
  run_effect_case_displacement_field_default_cycle_frames();
  run_effect_case_displacement_field_presets();
  run_effect_case_displacement_field_edit_cancels_departure();
  run_effect_case_displacement_field_startup_writes_preserve_entry();
  run_effect_case_displacement_field_paused_balls_deform_render();
  run_effect_case_displacement_field_paused_balls_keep_moving();
  run_effect_case_displacement_field_knot_walk_unit_near_ball();
  run_effect_case_displacement_field_zero_hue_scale_is_exact();
  run_effect_case_displacement_field_octave_bake_tracks_noise();
  run_effect_case_displacement_field_arm_boundaries_meet_endpoints();
  run_effect_case_displacement_field_single_pass_reference();
  run_effect_case_displacement_field_boundary_stencil();
  run_effect_case_displacement_field_turn_follows_large_steps();
  run_effect_case_displacement_field_alternate_groups_match_windows();
  if (effects_full_suite()) {
    run_effect_case_displacement_field_hue_table_fidelity();
    run_effect_case_displacement_field_hue_table_frame_fidelity();
    run_effect_case_displacement_field_clip_tiles_full();
    run_effect_case_displacement_field_ball_spans_and_lifecycle();
  }
  return fixture.result();
}

} // namespace hs_test::module_runners
