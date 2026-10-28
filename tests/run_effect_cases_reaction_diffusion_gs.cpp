/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Reaction diffusion gs test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/GSReactionDiffusion.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/reaction_diffusion_gs.h"

void run_effect_case_gs_palette_is_opaque() {
  HS_RUN(run_case, test_gs_palette_is_opaque);
}

void run_effect_case_gs_seed_palettes_and_pigment_transport() {
  HS_RUN(run_case, test_gs_seed_palettes_and_pigment_transport);
}

void run_effect_case_gs_sparse_pigment_matches_dense() {
  HS_RUN(run_case, test_gs_sparse_pigment_matches_dense);
}

void run_effect_case_gs_direct_draw_matches_grid() {
  HS_RUN(run_case, test_gs_direct_draw_matches_grid);
}

void run_effect_case_gs_dot_kernel_matches_squared_distance() {
  HS_RUN(run_case, test_gs_dot_kernel_matches_squared_distance);
}

void run_effect_case_gs_shared_cube_projection() {
  HS_RUN(run_case, test_gs_shared_cube_projection);
}

void run_effect_case_gs_shared_noise_palette_modifiers() {
  HS_RUN(run_case, test_gs_shared_noise_palette_modifiers);
}

void run_effect_case_gs_noise_bake_matches_reference() {
  HS_RUN(run_case, test_gs_noise_bake_matches_reference);
}

void run_effect_case_gs_reseed_generates_palette() {
  HS_RUN(run_case, test_gs_reseed_generates_palette);
}

void run_effect_case_gs_render_certificates_bound_lattice() {
  HS_RUN(run_case, test_gs_render_certificates_bound_lattice);
}

void run_effect_case_gs_shared_stencil_error_is_bounded() {
  HS_RUN(run_case, test_gs_shared_stencil_error_is_bounded);
}

void run_effect_case_gs_nearest_pigment_shader_matches_scalar_reference() {
  HS_RUN(run_case, test_gs_nearest_pigment_shader_matches_scalar_reference);
}

void run_effect_case_gs_signed_coverage_and_concentration() {
  HS_RUN(run_case, test_gs_signed_coverage_and_concentration);
}

void run_effect_case_gs_support_certificate_matches_display_geometry() {
  HS_RUN(run_case, test_gs_support_certificate_matches_display_geometry);
}

void run_effect_case_gs_nearest_pigment_shader_fidelity() {
  HS_RUN(run_case, test_gs_nearest_pigment_shader_fidelity);
}

void run_effect_case_gs_partial_color_palette_rows() {
  HS_RUN(run_case, test_gs_partial_color_palette_rows);
}

void run_effect_case_gs_reseed_invalidates_shimmer_cache() {
  HS_RUN(run_case, test_gs_reseed_invalidates_shimmer_cache);
}

void run_effect_case_gs_dissolve_frontier_fades_before_clear() {
  HS_RUN(run_case, test_gs_dissolve_frontier_fades_before_clear);
}

void run_effect_case_gs_substep_matches_scalar_reference() {
  HS_RUN(run_case, test_gs_substep_matches_scalar_reference);
}

void run_effect_case_gs_inplace_frame_matches_jacobi() {
  HS_RUN(run_case, test_gs_inplace_frame_matches_jacobi);
}

void run_effect_case_gs_hot_flags_match_directed_graph() {
  HS_RUN(run_case, test_gs_hot_flags_match_directed_graph);
}

void run_effect_case_gs_pigment_codec_layout() {
  HS_RUN(run_case, test_gs_pigment_codec_layout);
}

void run_effect_case_gs_rest_state_is_fixed_point() {
  HS_RUN(run_case, test_gs_rest_state_is_fixed_point);
}

void run_effect_case_gs_substep_signs_and_clamp() {
  HS_RUN(run_case, test_gs_substep_signs_and_clamp);
}

void run_effect_case_gs_reaction_edit_starts_dissolve() {
  HS_RUN(run_case, test_gs_reaction_edit_starts_dissolve);
}

void run_effect_case_gs_stall_hold_starts_dissolve() {
  HS_RUN(run_case, test_gs_stall_hold_starts_dissolve);
}

void run_effect_case_gs_preset_pause_does_not_hold_lifecycle() {
  HS_RUN(run_case, test_gs_preset_pause_does_not_hold_lifecycle);
}

void run_effect_case_gs_dissolve_absorbs_edits_and_stages_reseed() {
  HS_RUN(run_case, test_gs_dissolve_absorbs_edits_and_stages_reseed);
}

void run_effect_case_gs_untimed_frames_are_pinned() {
  HS_RUN(run_case, test_gs_untimed_frames_are_pinned);
}

void run_effect_case_gs_lattice_basis_matches_rotated_lattice() {
  HS_RUN(run_case, test_gs_lattice_basis_matches_rotated_lattice);
}

void run_effect_case_gs_slabs_shade_the_field_after_their_substep() {
  HS_RUN(run_case, test_gs_slabs_shade_the_field_after_their_substep);
}

void run_effect_case_gs_slab_noise_keeps_the_previous_bake() {
  HS_RUN(run_case, test_gs_slab_noise_keeps_the_previous_bake);
}

void run_effect_case_gs_segmented_slabs_match_full_frame() {
  HS_RUN(run_case, test_gs_segmented_slabs_match_full_frame);
}

void run_effect_case_gs_evolution_stays_bounded() {
  HS_RUN(run_case, test_gs_evolution_stays_bounded);
}

void run_effect_case_gs_reaction_corner_stays_bounded() {
  HS_RUN(run_case, test_gs_reaction_corner_stays_bounded);
}

void run_effect_case_gs_dissolve_clears_and_reseeds() {
  HS_RUN(run_case, test_gs_dissolve_clears_and_reseeds);
}

void run_effect_case_gs_staged_reseed_matches_synchronous() {
  HS_RUN(run_case, test_gs_staged_reseed_matches_synchronous);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the reaction diffusion gs effect cases.
 * @return Module failure count. */
int run_reaction_diffusion_gs_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_gs_palette_is_opaque();
  run_effect_case_gs_seed_palettes_and_pigment_transport();
  run_effect_case_gs_sparse_pigment_matches_dense();
  run_effect_case_gs_direct_draw_matches_grid();
  run_effect_case_gs_dot_kernel_matches_squared_distance();
  run_effect_case_gs_shared_cube_projection();
  run_effect_case_gs_shared_noise_palette_modifiers();
  run_effect_case_gs_noise_bake_matches_reference();
  run_effect_case_gs_reseed_generates_palette();
  run_effect_case_gs_render_certificates_bound_lattice();
  run_effect_case_gs_shared_stencil_error_is_bounded();
  run_effect_case_gs_nearest_pigment_shader_matches_scalar_reference();
  run_effect_case_gs_signed_coverage_and_concentration();
  run_effect_case_gs_support_certificate_matches_display_geometry();
  run_effect_case_gs_nearest_pigment_shader_fidelity();
  run_effect_case_gs_partial_color_palette_rows();
  run_effect_case_gs_reseed_invalidates_shimmer_cache();
  run_effect_case_gs_dissolve_frontier_fades_before_clear();
  run_effect_case_gs_substep_matches_scalar_reference();
  run_effect_case_gs_inplace_frame_matches_jacobi();
  run_effect_case_gs_hot_flags_match_directed_graph();
  run_effect_case_gs_pigment_codec_layout();
  run_effect_case_gs_rest_state_is_fixed_point();
  run_effect_case_gs_substep_signs_and_clamp();
  run_effect_case_gs_reaction_edit_starts_dissolve();
  run_effect_case_gs_stall_hold_starts_dissolve();
  run_effect_case_gs_preset_pause_does_not_hold_lifecycle();
  run_effect_case_gs_dissolve_absorbs_edits_and_stages_reseed();
  run_effect_case_gs_untimed_frames_are_pinned();
  run_effect_case_gs_lattice_basis_matches_rotated_lattice();
  run_effect_case_gs_slabs_shade_the_field_after_their_substep();
  run_effect_case_gs_slab_noise_keeps_the_previous_bake();
  run_effect_case_gs_segmented_slabs_match_full_frame();
  if (effects_full_suite()) {
    run_effect_case_gs_evolution_stays_bounded();
    run_effect_case_gs_reaction_corner_stays_bounded();
    run_effect_case_gs_dissolve_clears_and_reseeds();
    run_effect_case_gs_staged_reseed_matches_synchronous();
  }
  return fixture.result();
}

} // namespace hs_test::module_runners
