/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Reaction diffusion bz test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/BZReactionDiffusion.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/reaction_diffusion_bz.h"

void run_effect_case_bz_q16_roundtrip() {
  HS_RUN(run_case, test_bz_q16_roundtrip);
}

void run_effect_case_bz_ordered_neighbors() {
  HS_RUN(run_case, test_bz_ordered_neighbors);
}

void run_effect_case_bz_advance_species_signs_and_clamp() {
  HS_RUN(run_case, test_bz_advance_species_signs_and_clamp);
}

void run_effect_case_bz_perturb_state_draw_count_pinned() {
  HS_RUN(run_case, test_bz_perturb_state_draw_count_pinned);
}

void run_effect_case_bz_authored_palette() {
  HS_RUN(run_case, test_bz_authored_palette);
}

void run_effect_case_bz_min_diffusion_step_survives_quantization() {
  HS_RUN(run_case, test_bz_min_diffusion_step_survives_quantization);
}

void run_effect_case_bz_perturb_state_saturates_and_nudges() {
  HS_RUN(run_case, test_bz_perturb_state_saturates_and_nudges);
}

void run_effect_case_bz_perturb_scales_with_timestep() {
  HS_RUN(run_case, test_bz_perturb_scales_with_timestep);
}

void run_effect_case_bz_substep_diffuses() {
  HS_RUN(run_case, test_bz_substep_diffuses);
}

void run_effect_case_bz_substep_predator_wiring() {
  HS_RUN(run_case, test_bz_substep_predator_wiring);
}

void run_effect_case_bz_raster_matches_reference() {
  HS_RUN(run_case, test_bz_raster_matches_reference);
}

void run_effect_case_bz_render_center_matches_reference() {
  HS_RUN(run_case, test_bz_render_center_matches_reference);
}

void run_effect_case_bz_untimed_frames_are_pinned() {
  HS_RUN(run_case, test_bz_untimed_frames_are_pinned);
}

void run_effect_case_bz_slabs_shade_the_field_after_their_substep() {
  HS_RUN(run_case, test_bz_slabs_shade_the_field_after_their_substep);
}

void run_effect_case_bz_segmented_slabs_match_full_frame() {
  HS_RUN(run_case, test_bz_segmented_slabs_match_full_frame);
}

void run_effect_case_bz_lattice_basis_matches_rotated_lattice() {
  HS_RUN(run_case, test_bz_lattice_basis_matches_rotated_lattice);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the reaction diffusion bz effect cases.
 * @return Module failure count. */
int run_reaction_diffusion_bz_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_bz_q16_roundtrip();
  run_effect_case_bz_ordered_neighbors();
  run_effect_case_bz_advance_species_signs_and_clamp();
  run_effect_case_bz_perturb_state_draw_count_pinned();
  run_effect_case_bz_authored_palette();
  run_effect_case_bz_min_diffusion_step_survives_quantization();
  run_effect_case_bz_perturb_state_saturates_and_nudges();
  run_effect_case_bz_perturb_scales_with_timestep();
  run_effect_case_bz_substep_diffuses();
  run_effect_case_bz_substep_predator_wiring();
  run_effect_case_bz_raster_matches_reference();
  run_effect_case_bz_render_center_matches_reference();
  run_effect_case_bz_untimed_frames_are_pinned();
  run_effect_case_bz_slabs_shade_the_field_after_their_substep();
  run_effect_case_bz_segmented_slabs_match_full_frame();
  run_effect_case_bz_lattice_basis_matches_rotated_lattice();
  return fixture.result();
}

} // namespace hs_test::module_runners
