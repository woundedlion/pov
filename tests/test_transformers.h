/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Transformer test-module entry point. */

#include "tests/transformers/orientation.h"
#include "tests/transformers/ripple.h"
#include "tests/transformers/noise.h"
#include "tests/transformers/pool.h"
#include "tests/transformers/fields.h"

namespace hs_test::transformers_tests {

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every transformers test case.
 * @return The module's failure count.
 */
inline int run_transformers_tests() {
  hs_test::ModuleFixture fixture("transformers");

  HS_RUN(test_orient_transformer_identity);
  HS_RUN(test_orient_transformer_known_rotation);
  HS_RUN(test_mobius_identity_roundtrip);
  HS_RUN(test_mobius_known_rotation);
  HS_RUN(test_mobius_matches_double_precision_oracle);
  HS_RUN(test_mobius_poles_map_to_coefficient_ratios);
  HS_RUN(test_gnomonic_mobius_identity_roundtrip);
  HS_RUN(test_gnomonic_mobius_known_rotation);
  HS_RUN(test_gnomonic_mobius_equator_rotation);
  HS_RUN(test_gnomonic_mobius_signed_zero_equator);
  HS_RUN(test_ripple_zero_amplitude_is_identity);
  HS_RUN(test_ripple_center_point_is_identity);
  HS_RUN(test_ripple_active_rotates_on_sphere);
  HS_RUN(test_ripple_threshold_reject_path);
  HS_RUN(test_ripple_decay_attenuates);
  HS_RUN(test_ripple_threshold_boundary);
  HS_RUN(test_ripple_kernel_endpoint_order);
  HS_RUN(test_ripple_kernel_support_bands);
  HS_RUN(test_ripple_kernel_angle_at_distance);
  HS_RUN(test_ripple_kernel_sample_transfer);
  HS_RUN(test_ripple_small_angle_series_matches_exact);
  HS_RUN(test_ripple_threshold_collapses_past_pi);
  HS_RUN(test_transforms_nonfinite_passes_through_identity);
  HS_RUN(test_noise_zero_amplitude_is_identity);
  HS_RUN(test_noise_active_stays_on_sphere);
  HS_RUN(test_feedback_warp_inputs_preserve_fold_contract);
  HS_RUN(test_noise_cross_hemisphere_cap);
  HS_RUN(test_transformer_no_entities_is_identity);
  HS_RUN(test_transformer_spawn_applies_and_composes);
  HS_RUN(test_noise_transformer_spawn_synchronizes_frequency);
  HS_RUN(test_transformer_live_edit_preserves_instance_phase_and_seed);
  HS_RUN(test_transformer_pinned_repeat_retains_slot);
  HS_RUN(test_mobius_circular_alias_kernel);
  HS_RUN(test_mobius_gnomonic_alias_kernel);
  HS_RUN(test_transformer_spawn_pausable_freezes_start_delay);
  HS_RUN(test_transformer_nonpinned_slot_reclaimed_after_compaction);
  HS_RUN(test_transformer_full_timeline_releases_slot);
  HS_RUN(test_transformer_callback_after_pool_destroyed);
  HS_RUN(test_transformed_angle_rejects_nonfinite_directions);
  HS_RUN(test_ripple_step_keeps_its_previous_endpoint);
  HS_RUN(test_transformer_substitute_params_compose_in_spawn_order);
  HS_RUN(test_transformer_slots_released_by_timeline_clear);
  HS_RUN(test_transformer_recycled_slot_composes_in_spawn_order);
  HS_RUN(test_field_transformer_no_entities_is_zero);
  HS_RUN(test_field_transformer_sums_and_bounds);
  HS_RUN(test_field_transformer_active_params_spawn_order);
  HS_RUN(test_field_transformer_slot_reclaimed);
  HS_RUN(test_field_transformer_storage_survives_arena_rewind);
  HS_RUN(test_dominant_field_identities);
  HS_RUN(test_dominant_field_strongest_wins_without_stacking);
  HS_RUN(test_dominant_field_sign_crossing_is_continuous);
  HS_RUN(test_bump_field_threshold_sync);
  HS_RUN(test_bump_field_drapes_over_ball);
  HS_RUN(test_bump_field_round_bulge_along_ring);
  HS_RUN(test_bump_field_precomputed_y_parity);
  HS_RUN(test_bump_field_with_y_near_center_stays_bounded);
  HS_RUN(test_bump_field_envelope_gates);
  HS_RUN(test_bump_field_bound_is_conservative);
  HS_RUN(test_noise_params_seed_mirrors_generator);
  HS_RUN(test_noise_product_field_is_linear_and_separable);
  HS_RUN(test_noise_product_field_bound_is_conservative);
  HS_RUN(test_ball_drop_traverses_and_reclaims);
  HS_RUN(test_noise_product_integrates_time);

  return fixture.result();
}

} // namespace hs_test::transformers_tests
