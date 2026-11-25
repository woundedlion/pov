/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Shapeshifter oracle test-module entry. */

#include "tests/test_shapeshifter_surface.h"
#include "tests/shapeshifter_oracle/base_cases.h"

namespace hs_test::shapeshifter_oracle_tests {

/** @brief Runs shapeshifter oracle tests.
 * @return Failure count. */
inline int run_shapeshifter_oracle_tests() {
  ModuleFixture fixture("shapeshifter_oracle");
  if (hs_test::env_var("HS_SURFACE_DARK_DIAGNOSTICS")) {
    HS_RUN(test_surface_footprint_darkening);
    return fixture.result();
  }
  if (hs_test::env_var("HS_SURFACE_LOBE_DIAGNOSTICS")) {
    HS_RUN(test_surface_lobe_outliers);
    return fixture.result();
  }
  HS_RUN(radial_surface_tests::run_radial_surface_checks);
  HS_RUN(radial_surface_mask_tests::run_radial_surface_mask_checks);
  HS_RUN(test_surface_seams);
  HS_RUN(test_surface_clips);
  HS_RUN(test_surface_mask_parity);
  HS_RUN(test_surface_preset_dispatch);
  HS_RUN(test_surface_max_amplitude);
  HS_RUN(test_chart_edge_rotation_radius);
  HS_RUN(test_polar_walk_matches_rotated_chart);
  HS_RUN(test_timed_contours_match_frozen_columns);
  HS_RUN(test_timed_preset_five_preserves_polar_anchor);
  HS_RUN(test_timed_presets_match_frozen_columns);
  HS_RUN(test_timed_flower_clip_culls_anchor_intervals);
  HS_RUN(test_timed_flower_mid_arm_clip_preserves_motion);
  HS_RUN(test_timed_dense_star_clip_samples_polar_columns);
  HS_RUN(test_count_slider_matches_draw_limit);
  HS_RUN(test_contour_layout_storage_alignment);
  HS_RUN(test_buffer_comparator_statistics);
  HS_RUN(test_reference_matrix_is_deterministic_and_nonblack);
  HS_RUN(test_candidate_matrix_stays_within_visual_budget);
  HS_RUN(test_star_projection_policies_render_different_edges);
  HS_RUN(test_screen_balanced_spacing_is_opt_in_for_every_shape);
  HS_RUN(test_star_options_and_shipping_presets_are_planar);
  HS_RUN(test_high_count_uniform_star_stays_within_visual_budget);
  HS_RUN(test_screen_balanced_star_stroke_stays_within_visual_budget);
  HS_RUN(test_timed_dense_star_stroke_stays_within_visual_budget);
  HS_RUN(test_high_count_uniform_star_covers_north_pole);
  HS_RUN(test_high_count_planar_star_caps_cover_chart_centers);
  HS_RUN(test_high_count_spherical_star_contours_reach_display_north);
  HS_RUN(test_amplitude_preserves_sweep_velocity);
  HS_RUN(test_phase_step_stops_at_the_nyquist_limit);
  HS_RUN(test_alpha_falloff_fades_to_equator);
  HS_RUN(test_opposite_halves_direction);
  HS_RUN(test_folded_shapes_draw_from_equator_to_both_poles);
  HS_RUN(test_screen_balanced_spacing_follows_sampling_envelope);
  HS_RUN(test_preset_transition_adopts_next_params);
  HS_RUN(test_preset_transition_fades_through_black);
  return fixture.result();
}
} // namespace hs_test::shapeshifter_oracle_tests
