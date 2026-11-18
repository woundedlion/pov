/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Composed effect test-module entry. */

#include "tests/composed_effect/generated_catalog_assignments.h"
#include "tests/composed_effect/derivation_reach.h"
#include "tests/composed_effect/document_values.h"
#include "tests/composed_effect/roster_contracts.h"
#include "tests/composed_effect/stage_probes.h"
#include "tests/composed_effect/authored_palette.h"
#include "tests/composed_effect/parameter_schema.h"
#include "tests/composed_effect/temporal_sampling.h"
#include "tests/composed_effect/retained_hue.h"
#include "tests/composed_effect/chamber_noise_cache.h"
#include "tests/composed_effect/dodecahedral_matrix_cache.h"
#include "tests/composed_effect/ashcloud_shared_color.h"
#include "tests/composed_effect/palette_rgb_cache.h"
#include "tests/composed_effect/grid_worms.h"
#include "tests/composed_effect/authored_policies.h"
#include "tests/composed_effect/base_cases.h"

namespace hs_test::composed_effect_tests {

/** @brief Module entry point for the composed-effect base contract.
 * @return Failure count. */
inline int run_composed_effect_tests() {
  ModuleFixture fixture("composed_effect");
  HS_RUN(test_authored_tangent_bindings_match_operators);
  HS_RUN(test_authored_tangent_fields_after_twist);
  HS_RUN(test_authored_tangent_interval_uses_radians);
  HS_RUN(test_authored_potential_interpolation_and_admission);
  HS_RUN(test_family_refusal_warnings);
  HS_RUN(test_authored_interpolation_destination_aliases);
  HS_RUN(test_authored_surface_noise_topology);
  HS_RUN(test_authored_dynamic_modulation_storage);
  HS_RUN(test_authored_projection_bindings_match_operators);
  HS_RUN(test_authored_dynamic_modulation_controls);
  HS_RUN(test_authored_potential_nested_controls);
  HS_RUN(test_authored_noise_seed_identity);
  HS_RUN(test_millefiori_radial_table_interval);
  HS_RUN(test_millefiori_document_pipeline);
  HS_RUN(test_promoted_document_pipelines);
  HS_RUN(test_grid_worms_twist_clocks);
  HS_RUN(test_grid_worms_document_pipeline);
  HS_RUN(test_ash_cloud_document_pipeline);
  HS_RUN(test_ash_cloud_shared_color_identity);
  HS_RUN(test_ash_cloud_shared_column_qualification);
  HS_RUN(test_ash_cloud_shared_column_fallbacks);
  HS_RUN(test_composed_default_gradient_preparation);
  HS_RUN(test_temporal_constant_families_are_exact);
  HS_RUN(test_composed_spec_default_color_paths);
  HS_RUN(test_composed_persistent_prepared_frame);
  HS_RUN(test_composed_document_modulation_members);
  HS_RUN(test_millefiori_column_color_phase);
  HS_RUN(test_grid_worms_sine_column_phase);
  HS_RUN(test_retained_equal_palette_hue_view);
  HS_RUN(test_retained_palette_identity_rejects_changes);
  HS_RUN(test_ash_cloud_noise_column_phase);
  HS_RUN(test_temporal_direct_noise_endpoint_reuse);
  HS_RUN(test_river_stones_color_tail);
  HS_RUN(test_baked_cube_column_phase);
  HS_RUN(test_chamber_noise_cache_bake_lifetime);
  HS_RUN(test_chamber_noise_cache_admission);
  HS_RUN(test_chamber_noise_cache_phase_and_fallback);
  HS_RUN(test_chamber_noise_cache_endpoint_admission);
  HS_RUN(test_chamber_noise_cache_equal_clock_fraction);
  HS_RUN(test_chamber_noise_cache_far_endpoint_falls_back);
  HS_RUN(test_chamber_noise_cache_rotated_endpoint_falls_back);
  HS_RUN(test_chamber_noise_cache_source_phase_invalidates_patch);
  HS_RUN(test_chamber_noise_cache_runtime_history);
  HS_RUN(test_dodecahedral_matrix_cache_lifetime);
  HS_RUN(test_dodecahedral_matrix_cache_fold_parity);
  HS_RUN(test_dodecahedral_matrix_cache_lens_binding);
  HS_RUN(test_dodecahedral_face_cache_lifetime);
  HS_RUN(test_dodecahedral_face_cache_fold_parity);
  HS_RUN(test_dodecahedral_face_cache_lens_binding);
  HS_RUN(test_palette_rgb_cache_axes_and_wrap);
  HS_RUN(test_palette_rgb_cache_guard_mask_aliasing);
  HS_RUN(test_palette_lab_reference_endpoint_lifetime);
  HS_RUN(test_palette_rgb_cache_matches_hue_rotation_table);
  HS_RUN(test_palette_rgb_cache_lazy_guards);
  HS_RUN(test_palette_rgb_cache_runtime_lifetime);
  HS_RUN(test_palette_rgb_cache_bake_reuse);
  HS_RUN(test_palette_rgb_cache_shared_base_endpoints);
  HS_RUN(test_palette_rgb_cache_midpoint_reuse);
  HS_RUN(test_palette_rgb_cache_any_chroma);
  HS_RUN(test_palette_rgb_cache_q16_compatibility);
  HS_RUN(test_palette_rgb_midpoint_only_storage);
  HS_RUN(test_palette_rgb_midpoint_only_fallback);
  HS_RUN(test_palette_lab_cache_color_padding);
  HS_RUN(test_palette_lab_cache_column_controls);
  HS_RUN(test_palette_lab_cache_through_lerp);
  HS_RUN(test_palette_lab_cache_palette_control_step);
  HS_RUN(test_millefiori_palette_lab_lifetime);
  HS_RUN(test_millefiori_ripple_endpoints);
  HS_RUN(test_composed_strobe_roster);
  HS_RUN(test_composed_analytic_modulation_sliders);
  HS_RUN(test_composed_authored_palette_sources);
  HS_RUN(test_composed_authored_palette_harmony_clock);
  HS_RUN(test_composed_authored_palette_admission);
  HS_RUN(test_composed_analytic_modulation_admission);
  HS_RUN(test_composed_path_modulation_admission);
  HS_RUN(test_temporal_clock_sampling);
  HS_RUN(test_temporal_affine_cycle);
  HS_RUN(test_temporal_custom_shade_is_authoritative);
  HS_RUN(test_scoped_temporal_sampling_restores_entering_mode);
  HS_RUN(test_gridspace_column_basis_geometry);
  HS_RUN(test_gridspace_canonical_columns);
  HS_RUN(test_temporal_equal_palettes_color_as_current);
  HS_RUN(test_temporal_moving_families_interpolate);
  HS_RUN(test_composed_modulation_blends_axes);
  HS_RUN(test_temporal_owned_resources);
  HS_RUN(test_temporal_topology_interval_selects_current_endpoint);
  HS_RUN(test_modulation_interpolate_retains_target_topology);
  HS_RUN(test_composed_factored_columns_match_pixels);
  HS_RUN(test_temporal_modulation_admission);
  HS_RUN(test_temporal_modulation_bake_generation);
  HS_RUN(test_temporal_stacked_modulation_resources);
  HS_RUN(test_gridspace_temporal_pixels);
  HS_RUN(test_composed_temporal_roster);
  HS_RUN(test_mobius_grid_instantaneous_clips_match_full_render);
  HS_RUN(test_mobius_grid_lens_follows_phase);
  HS_RUN(test_mobius_interval_exponential);
  HS_RUN(test_mobius_grid_segmented_matches_single_engine);
  HS_RUN(test_alien_brain_preset_dwell);
  HS_RUN(test_mobius_grid_circular_animation);
  HS_RUN(test_json_number_grammar);
  HS_RUN(test_json_string_grammar);
  HS_RUN(test_catalog_semantic_export);
  HS_RUN(test_modulation_schema_defaults_match_controls);
  HS_RUN(test_composed_hand_registered_families);
  HS_RUN(test_composed_slider_registration);
  HS_RUN(test_composed_modulation_is_device_qualified);
  HS_RUN(test_composed_noise_modulation_scale_limits);
  HS_RUN(check_slider_registration<ProjectedNoiseProbe>, "ProjectedNoiseProbe");
  HS_RUN(test_mobius_frame_admission);
  HS_RUN(test_mobius_generator_interpolation);
  HS_RUN(test_mobius_captured_parameter_restore);
  HS_RUN(test_mobius_automatic_departures_interpolate_generator);
  HS_RUN(test_composed_snapshot_contract);
  HS_RUN(test_parameter_layout_reorder);
  HS_RUN(test_composed_parameter_schema_pins);
  HS_RUN(test_composed_preset_choreography);
  HS_RUN(test_composed_log_positive_curve);
  HS_RUN(test_composed_preset_interpolation);
  HS_RUN(test_palette_mapping_crossfade);
  HS_RUN(test_composed_document_values);
  HS_RUN(test_document_numeric_resource_selection);
  HS_RUN(test_flowers_longitude_seam);
  HS_RUN(test_composed_keyed_parameters);
  HS_RUN(test_composed_keyed_value_providers);
  HS_RUN(test_composed_affine_named_source);
  HS_RUN(test_composed_affine_form_sliders);
  HS_RUN(test_composed_repeated_instances);
  HS_RUN(test_composed_projection_frame_identity);
  HS_RUN(test_composed_projection_walk_storage);
  HS_RUN(test_composed_periodic_ripple_surface);
  HS_RUN(test_composed_noise_source);
  HS_RUN(hs_test::catalog_policy_tests::test_palette_selection_schema_parity);
  HS_RUN(hs_test::catalog_policy_tests::test_authored_warp_catalog_parity);
  HS_RUN(hs_test::catalog_policy_tests::test_authored_source_catalog_parity);
  HS_RUN(hs_test::catalog_policy_tests::test_authored_transfer_catalog_parity);
  HS_RUN(hs_test::catalog_policy_tests::test_authored_sample_crossing);
  HS_RUN(
      hs_test::catalog_policy_tests::test_authored_twin_wave_sample_admission);
  HS_RUN(hs_test::catalog_policy_tests::test_edge_fade_squared_curve);
  HS_RUN(hs_test::keyed_rotation_tests::test_keyed_rotation_runtime_parity);
  HS_RUN(test_choreography_fade_departure);
  HS_RUN(test_choreography_fade_resumes_after_pause);
  HS_RUN(test_choreography_lerp_transition_hooks);
  HS_RUN(test_choreography_lerp_pause_policy);
  return fixture.result();
}
} // namespace hs_test::composed_effect_tests
