/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Mesh effect test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/HankinSolids.h"
#include "effects/IslamicStars.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/hankin.h"
#include "tests/effects/islamic_stars.h"
#include "tests/effects/mesh_capture_cases.h"

void run_effect_case_hankinsolids_handoff_history() {
  HS_RUN(run_case, test_hankinsolids_handoff_history);
}

void run_effect_case_hankinsolids_interlace_motion_is_seam_free() {
  HS_RUN(run_case, test_hankinsolids_interlace_motion_is_seam_free);
}

void run_effect_case_hankinsolids_seam_copies_are_one_polygon() {
  HS_RUN(run_case, test_hankinsolids_seam_copies_are_one_polygon);
}

void run_effect_case_hankinsolids_morph_motion_is_seam_free() {
  HS_RUN(run_case, test_hankinsolids_morph_motion_is_seam_free);
}

void run_effect_case_islamicstars_ripple_motion_is_seam_free() {
  HS_RUN(run_case, test_islamicstars_ripple_motion_is_seam_free);
}

void run_effect_case_islamicstars_seam_copies_are_one_polygon() {
  HS_RUN(run_case, test_islamicstars_seam_copies_are_one_polygon);
}

void run_effect_case_islamicstars_playback_motion_is_seam_free() {
  HS_RUN(run_case, test_islamicstars_playback_motion_is_seam_free);
}

void run_effect_case_islamicstars_cases() { HS_RUN(run_islamicstars_cases); }

void run_effect_case_hankinsolids_manual_pause_holds_morph() {
  HS_RUN(run_case, test_hankinsolids_manual_pause_holds_morph);
}

void run_effect_case_hankinsolids_handoffs_shade_at_column_time() {
  HS_RUN(run_case, test_hankinsolids_handoffs_shade_at_column_time);
}

void run_effect_case_hankinsolids_conway_crossfade_shades_at_column_time() {
  HS_RUN(run_case, test_hankinsolids_conway_crossfade_shades_at_column_time);
}

void run_effect_case_hankinsolids_morph_draws_carry_previous_vertices() {
  HS_RUN(run_case, test_hankinsolids_morph_draws_carry_previous_vertices);
}

void run_effect_case_hankinsolids_arena_budget_covers_every_solid() {
  HS_RUN(run_case, test_hankinsolids_arena_budget_covers_every_solid);
}

void run_effect_case_islamicstars_recipe_build_smoke() {
  HS_RUN(run_case, test_islamicstars_recipe_build_smoke);
}

void run_effect_case_islamicstars_roster_cycle_fits_budget() {
  HS_RUN(run_case, test_islamicstars_roster_cycle_fits_budget);
}

void run_effect_case_islamicstars_dual_bridge_fits_budget() {
  HS_RUN(run_case, test_islamicstars_dual_bridge_fits_budget);
}

void run_effect_case_islamicstars_motion_fits_scratch() {
  HS_RUN(run_case, test_islamicstars_motion_fits_scratch);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the mesh effect effect cases.
 * @return Module failure count. */
int run_mesh_effect_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_hankinsolids_manual_pause_holds_morph();
  run_effect_case_hankinsolids_handoffs_shade_at_column_time();
  run_effect_case_hankinsolids_handoff_history();
  run_effect_case_hankinsolids_conway_crossfade_shades_at_column_time();
  run_effect_case_hankinsolids_interlace_motion_is_seam_free();
  run_effect_case_hankinsolids_seam_copies_are_one_polygon();
  run_effect_case_hankinsolids_morph_draws_carry_previous_vertices();
  run_effect_case_hankinsolids_morph_motion_is_seam_free();
  run_effect_case_hankinsolids_arena_budget_covers_every_solid();
  run_effect_case_islamicstars_cases();
  if (effects_full_suite()) {
    run_effect_case_islamicstars_recipe_build_smoke();
    run_effect_case_islamicstars_roster_cycle_fits_budget();
    run_effect_case_islamicstars_dual_bridge_fits_budget();
    run_effect_case_islamicstars_motion_fits_scratch();
  }
  return fixture.result();
}

} // namespace hs_test::module_runners
