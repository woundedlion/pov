/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Effect case entry points shared by suite and focused targets. */
#pragma once

namespace hs_test::effects_tests {

/** @brief Runs RingShowerWhiteBox check_temporal_rings with its case context. */
void run_effect_case_ring_shower_white_box_check_temporal_rings();

/** @brief Runs ThrustersWhiteBox check_temporal_rings with its case context. */
void run_effect_case_thrusters_white_box_check_temporal_rings();

/** @brief Runs MobiusRingsWhiteBox check_column_time_oracle with its case context. */
void run_effect_case_mobius_rings_white_box_check_column_time_oracle();

/** @brief Runs test_dynamo_column_time_oracle with its case context. */
void run_effect_case_dynamo_column_time_oracle();

/** @brief Runs test_dynamo_column_time_first_frame_and_reset with its case context. */
void run_effect_case_dynamo_column_time_first_frame_and_reset();

/** @brief Runs test_dynamo_emitted_points_counts_ring_seeds with its case context. */
void run_effect_case_dynamo_emitted_points_counts_ring_seeds();

/** @brief Runs test_dynamo_trail_ceiling_bounds_the_ring with its case context. */
void run_effect_case_dynamo_trail_ceiling_bounds_the_ring();

/** @brief Runs test_hankinsolids_handoff_history with its case context. */
void run_effect_case_hankinsolids_handoff_history();

/** @brief Runs test_hankinsolids_interlace_motion_is_seam_free with its case context. */
void run_effect_case_hankinsolids_interlace_motion_is_seam_free();

/** @brief Runs test_hankinsolids_seam_copies_are_one_polygon with its case context. */
void run_effect_case_hankinsolids_seam_copies_are_one_polygon();

/** @brief Runs test_hankinsolids_morph_motion_is_seam_free with its case context. */
void run_effect_case_hankinsolids_morph_motion_is_seam_free();

/** @brief Runs test_islamicstars_ripple_motion_is_seam_free with its case context. */
void run_effect_case_islamicstars_ripple_motion_is_seam_free();

/** @brief Runs test_islamicstars_seam_copies_are_one_polygon with its case context. */
void run_effect_case_islamicstars_seam_copies_are_one_polygon();

/** @brief Runs test_islamicstars_playback_motion_is_seam_free with its case context. */
void run_effect_case_islamicstars_playback_motion_is_seam_free();

/** @brief Runs run_islamicstars_cases with its case context. */
void run_effect_case_islamicstars_cases();

/** @brief Runs test_hopf_temporal_trail_oracle with its case context. */
void run_effect_case_hopf_temporal_trail_oracle();

/** @brief Runs test_dreamballs_temporal_source_oracle with its case context. */
void run_effect_case_dreamballs_temporal_source_oracle();

/** @brief Runs test_dreamballs_column_work_is_bounded with its case context. */
void run_effect_case_dreamballs_column_work_is_bounded();

/** @brief Runs DreamBallsWhiteBox check_endpoint_error_bound with its case context. */
void run_effect_case_dream_balls_white_box_check_endpoint_error_bound();

/** @brief Runs test_petalflow_temporal_source_oracle with its case context. */
void run_effect_case_petalflow_temporal_source_oracle();

/** @brief Runs test_meshfeedback_temporal_source_oracle with its case context. */
void run_effect_case_meshfeedback_temporal_source_oracle();

/** @brief Runs test_meshfeedback_stationary_timed_source with its case context. */
void run_effect_case_meshfeedback_stationary_timed_source();

/** @brief Runs test_meshfeedback_trails_arm_cut_seam with its case context. */
void run_effect_case_meshfeedback_trails_arm_cut_seam();

/** @brief Runs test_parameter_probe_targets with its case context. */
void run_effect_case_parameter_probe_targets();

/** @brief Runs test_float_dropdown_lint with its case context. */
void run_effect_case_float_dropdown_lint();

/** @brief Runs test_ringspin_strobe_configuration_preserves_rendering with its case context. */
void run_effect_case_ringspin_strobe_configuration_preserves_rendering();

/** @brief Runs test_meshfeedback_base_mesh_selector with its case context. */
void run_effect_case_meshfeedback_base_mesh_selector();

/** @brief Runs test_meshfeedback_preset_export_arity with its case context. */
void run_effect_case_meshfeedback_preset_export_arity();

/** @brief Runs test_meshfeedback_mesh_rebuild_reuses_storage with its case context. */
void run_effect_case_meshfeedback_mesh_rebuild_reuses_storage();

/** @brief Runs test_gs_palette_is_opaque with its case context. */
void run_effect_case_gs_palette_is_opaque();

/** @brief Runs test_gs_seed_palettes_and_pigment_transport with its case context. */
void run_effect_case_gs_seed_palettes_and_pigment_transport();

/** @brief Runs test_gs_sparse_pigment_matches_dense with its case context. */
void run_effect_case_gs_sparse_pigment_matches_dense();

/** @brief Runs test_gs_direct_draw_matches_grid with its case context. */
void run_effect_case_gs_direct_draw_matches_grid();

/** @brief Runs test_gs_dot_kernel_matches_squared_distance with its case context. */
void run_effect_case_gs_dot_kernel_matches_squared_distance();

/** @brief Runs test_gs_shared_cube_projection with its case context. */
void run_effect_case_gs_shared_cube_projection();

/** @brief Runs test_gs_shared_noise_palette_modifiers with its case context. */
void run_effect_case_gs_shared_noise_palette_modifiers();

/** @brief Runs test_gs_noise_bake_matches_reference with its case context. */
void run_effect_case_gs_noise_bake_matches_reference();

/** @brief Runs test_gs_reseed_generates_palette with its case context. */
void run_effect_case_gs_reseed_generates_palette();

/** @brief Runs test_gs_render_certificates_bound_lattice with its case context. */
void run_effect_case_gs_render_certificates_bound_lattice();

/** @brief Runs test_gs_shared_stencil_error_is_bounded with its case context. */
void run_effect_case_gs_shared_stencil_error_is_bounded();

/** @brief Runs test_gs_nearest_pigment_shader_matches_scalar_reference with its case context. */
void run_effect_case_gs_nearest_pigment_shader_matches_scalar_reference();

/** @brief Runs test_gs_signed_coverage_and_concentration with its case context. */
void run_effect_case_gs_signed_coverage_and_concentration();

/** @brief Runs test_gs_support_certificate_matches_display_geometry with its case context. */
void run_effect_case_gs_support_certificate_matches_display_geometry();

/** @brief Runs test_gs_nearest_pigment_shader_fidelity with its case context. */
void run_effect_case_gs_nearest_pigment_shader_fidelity();

/** @brief Runs test_gs_partial_color_palette_rows with its case context. */
void run_effect_case_gs_partial_color_palette_rows();

/** @brief Runs test_gs_reseed_invalidates_shimmer_cache with its case context. */
void run_effect_case_gs_reseed_invalidates_shimmer_cache();

/** @brief Runs test_gs_dissolve_frontier_fades_before_clear with its case context. */
void run_effect_case_gs_dissolve_frontier_fades_before_clear();

/** @brief Runs test_gs_substep_matches_scalar_reference with its case context. */
void run_effect_case_gs_substep_matches_scalar_reference();

/** @brief Runs test_gs_inplace_frame_matches_jacobi with its case context. */
void run_effect_case_gs_inplace_frame_matches_jacobi();

/** @brief Runs test_fishbowl_preset_and_fire_duty_cycle with its case context. */
void run_effect_case_fishbowl_preset_and_fire_duty_cycle();

/** @brief Runs test_fishbowl_temporal_trail_oracle with its case context. */
void run_effect_case_fishbowl_temporal_trail_oracle();

/** @brief Runs test_sh_preset_mode_mapping with its case context. */
void run_effect_case_sh_preset_mode_mapping();

/** @brief Runs test_sh_manual_preset_replaces_inflight_morph with its case context. */
void run_effect_case_sh_manual_preset_replaces_inflight_morph();

/** @brief Runs test_shapeshifter_preset_defaults with its case context. */
void run_effect_case_shapeshifter_preset_defaults();

/** @brief Runs test_shapeshifter_slider_selections_render with its case context. */
void run_effect_case_shapeshifter_slider_selections_render();

/** @brief Runs test_manual_preset_navigation with its case context. */
void run_effect_case_manual_preset_navigation();

/** @brief Runs test_hankinsolids_manual_pause_holds_morph with its case context. */
void run_effect_case_hankinsolids_manual_pause_holds_morph();

/** @brief Runs test_hankinsolids_handoffs_shade_at_column_time with its case context. */
void run_effect_case_hankinsolids_handoffs_shade_at_column_time();

/** @brief Runs test_hankinsolids_conway_crossfade_shades_at_column_time with its case context. */
void run_effect_case_hankinsolids_conway_crossfade_shades_at_column_time();

/** @brief Runs test_hankinsolids_morph_draws_carry_previous_vertices with its case context. */
void run_effect_case_hankinsolids_morph_draws_carry_previous_vertices();

/** @brief Runs test_fishbowl_scratch_estimate_covers_peak with its case context. */
void run_effect_case_fishbowl_scratch_estimate_covers_peak();

/** @brief Runs test_hankinsolids_arena_budget_covers_every_solid with its case context. */
void run_effect_case_hankinsolids_arena_budget_covers_every_solid();

/** @brief Runs test_dreamballs_max_edge_solid_render with its case context. */
void run_effect_case_dreamballs_max_edge_solid_render();

/** @brief Runs test_raymarch_volume_random_walks_are_independent with its case context. */
void run_effect_case_raymarch_volume_random_walks_are_independent();

/** @brief Runs test_raymarch_preset_and_placement_solids with its case context. */
void run_effect_case_raymarch_preset_and_placement_solids();

/** @brief Runs test_raymarch_surface_frame_uv with its case context. */
void run_effect_case_raymarch_surface_frame_uv();

/** @brief Runs test_raymarch_lighting_and_palette_at_column_time with its case context. */
void run_effect_case_raymarch_lighting_and_palette_at_column_time();

/** @brief Runs test_raymarch_hue_noise_at_column_time with its case context. */
void run_effect_case_raymarch_hue_noise_at_column_time();

/** @brief Runs test_gs_hot_flags_match_directed_graph with its case context. */
void run_effect_case_gs_hot_flags_match_directed_graph();

/** @brief Runs test_gs_pigment_codec_layout with its case context. */
void run_effect_case_gs_pigment_codec_layout();

/** @brief Runs test_gs_rest_state_is_fixed_point with its case context. */
void run_effect_case_gs_rest_state_is_fixed_point();

/** @brief Runs test_gs_substep_signs_and_clamp with its case context. */
void run_effect_case_gs_substep_signs_and_clamp();

/** @brief Runs test_bz_q16_roundtrip with its case context. */
void run_effect_case_bz_q16_roundtrip();

/** @brief Runs test_bz_ordered_neighbors with its case context. */
void run_effect_case_bz_ordered_neighbors();

/** @brief Runs test_bz_advance_species_signs_and_clamp with its case context. */
void run_effect_case_bz_advance_species_signs_and_clamp();

/** @brief Runs test_bz_perturb_state_draw_count_pinned with its case context. */
void run_effect_case_bz_perturb_state_draw_count_pinned();

/** @brief Runs test_hopf_projection_math with its case context. */
void run_effect_case_hopf_projection_math();

/** @brief Runs test_raymarch_constexpr_sqrt_converges with its case context. */
void run_effect_case_raymarch_constexpr_sqrt_converges();

/** @brief Runs test_displacement_field_retained_history_resets with its case context. */
void run_effect_case_displacement_field_retained_history_resets();

/** @brief Runs test_dreamballs_retained_history_resets with its case context. */
void run_effect_case_dreamballs_retained_history_resets();

/** @brief Runs test_needs_full_frame_gate with its case context. */
void run_effect_case_needs_full_frame_gate();

/** @brief Runs test_voronoi_axes_use_uniform_sampler with its case context. */
void run_effect_case_voronoi_axes_use_uniform_sampler();

/** @brief Runs test_voronoi_segment_render_matches_full_frame with its case context. */
void run_effect_case_voronoi_segment_render_matches_full_frame();

/** @brief Runs test_voronoi_instantaneous_render_pins with its case context. */
void run_effect_case_voronoi_instantaneous_render_pins();

/** @brief Runs test_voronoi_endpoint_pairs_match_brute_force with its case context. */
void run_effect_case_voronoi_endpoint_pairs_match_brute_force();

/** @brief Runs test_voronoi_temporal_toggle with its case context. */
void run_effect_case_voronoi_temporal_toggle();

/** @brief Runs test_voronoi_reseed_resets_history with its case context. */
void run_effect_case_voronoi_reseed_resets_history();

/** @brief Runs test_voronoi_retained_corners_match_search with its case context. */
void run_effect_case_voronoi_retained_corners_match_search();

/** @brief Runs test_sh_field_blend_endpoints_and_frame with its case context. */
void run_effect_case_sh_field_blend_endpoints_and_frame();

/** @brief Runs test_sh_field_stays_inside_unit_range with its case context. */
void run_effect_case_sh_field_stays_inside_unit_range();

/** @brief Runs test_sh_polarity_split_and_ao_shaping with its case context. */
void run_effect_case_sh_polarity_split_and_ao_shaping();

/** @brief Runs test_sh_cycle_covers_every_mode with its case context. */
void run_effect_case_sh_cycle_covers_every_mode();

/** @brief Runs test_sh_instantaneous_pixel_golden with its case context. */
void run_effect_case_sh_instantaneous_pixel_golden();

/** @brief Runs test_sh_columns_follow_sweep_time with its case context. */
void run_effect_case_sh_columns_follow_sweep_time();

/** @brief Runs test_sh_segmented_columns_match_single_engine with its case context. */
void run_effect_case_sh_segmented_columns_match_single_engine();

/** @brief Runs test_sh_manual_preset_resets_morph_history with its case context. */
void run_effect_case_sh_manual_preset_resets_morph_history();

/** @brief Runs test_sh_fidelity_toggle_keeps_history with its case context. */
void run_effect_case_sh_fidelity_toggle_keeps_history();

/** @brief Runs test_sh_frame_sampler_rollover_and_sign with its case context. */
void run_effect_case_sh_frame_sampler_rollover_and_sign();

/** @brief Runs test_gs_reaction_edit_starts_dissolve with its case context. */
void run_effect_case_gs_reaction_edit_starts_dissolve();

/** @brief Runs test_gs_stall_hold_starts_dissolve with its case context. */
void run_effect_case_gs_stall_hold_starts_dissolve();

/** @brief Runs test_gs_preset_pause_does_not_hold_lifecycle with its case context. */
void run_effect_case_gs_preset_pause_does_not_hold_lifecycle();

/** @brief Runs test_gs_dissolve_absorbs_edits_and_stages_reseed with its case context. */
void run_effect_case_gs_dissolve_absorbs_edits_and_stages_reseed();

/** @brief Runs test_bz_authored_palette with its case context. */
void run_effect_case_bz_authored_palette();

/** @brief Runs test_bz_min_diffusion_step_survives_quantization with its case context. */
void run_effect_case_bz_min_diffusion_step_survives_quantization();

/** @brief Runs test_bz_perturb_state_saturates_and_nudges with its case context. */
void run_effect_case_bz_perturb_state_saturates_and_nudges();

/** @brief Runs test_bz_perturb_scales_with_timestep with its case context. */
void run_effect_case_bz_perturb_scales_with_timestep();

/** @brief Runs test_bz_substep_diffuses with its case context. */
void run_effect_case_bz_substep_diffuses();

/** @brief Runs test_bz_substep_predator_wiring with its case context. */
void run_effect_case_bz_substep_predator_wiring();

/** @brief Runs test_bz_raster_matches_reference with its case context. */
void run_effect_case_bz_raster_matches_reference();

/** @brief Runs test_bz_render_center_matches_reference with its case context. */
void run_effect_case_bz_render_center_matches_reference();

/** @brief Runs test_bz_untimed_frames_are_pinned with its case context. */
void run_effect_case_bz_untimed_frames_are_pinned();

/** @brief Runs test_bz_slabs_shade_the_field_after_their_substep with its case context. */
void run_effect_case_bz_slabs_shade_the_field_after_their_substep();

/** @brief Runs test_bz_segmented_slabs_match_full_frame with its case context. */
void run_effect_case_bz_segmented_slabs_match_full_frame();

/** @brief Runs test_bz_lattice_basis_matches_rotated_lattice with its case context. */
void run_effect_case_bz_lattice_basis_matches_rotated_lattice();

/** @brief Runs test_gs_untimed_frames_are_pinned with its case context. */
void run_effect_case_gs_untimed_frames_are_pinned();

/** @brief Runs test_gs_lattice_basis_matches_rotated_lattice with its case context. */
void run_effect_case_gs_lattice_basis_matches_rotated_lattice();

/** @brief Runs test_gs_slabs_shade_the_field_after_their_substep with its case context. */
void run_effect_case_gs_slabs_shade_the_field_after_their_substep();

/** @brief Runs test_gs_slab_noise_keeps_the_previous_bake with its case context. */
void run_effect_case_gs_slab_noise_keeps_the_previous_bake();

/** @brief Runs test_gs_segmented_slabs_match_full_frame with its case context. */
void run_effect_case_gs_segmented_slabs_match_full_frame();

/** @brief Runs test_dreamballs_preset_cycle_bookkeeping with its case context. */
void run_effect_case_dreamballs_preset_cycle_bookkeeping();

/** @brief Runs CometsWhiteBox check_paths_close with its case context. */
void run_effect_case_comets_white_box_check_paths_close();

/** @brief Runs CometsWhiteBox check_births_skip_earlier_columns with its case context. */
void run_effect_case_comets_white_box_check_births_skip_earlier_columns();

/** @brief Runs test_comets_rollover_skipped_mid_wipe with its case context. */
void run_effect_case_comets_rollover_skipped_mid_wipe();

/** @brief Runs test_comets_wipe_shades_at_column_time with its case context. */
void run_effect_case_comets_wipe_shades_at_column_time();

/** @brief Runs test_comets_wipe_history with its case context. */
void run_effect_case_comets_wipe_history();

/** @brief Runs ThrustersWhiteBox check_warp_endpoints with its case context. */
void run_effect_case_thrusters_white_box_check_warp_endpoints();

/** @brief Runs ThrustersWhiteBox check_fire_spawns_opposed_pair with its case context. */
void run_effect_case_thrusters_white_box_check_fire_spawns_opposed_pair();

/** @brief Runs ThrustersWhiteBox check_collapsed_ring_falls_back_to_a_spin_axis with its case context. */
void run_effect_case_thrusters_white_box_check_collapsed_ring_falls_back_to_a_spin_axis();

/** @brief Runs ThrustersWhiteBox check_fifo_evicts_the_oldest_pair with its case context. */
void run_effect_case_thrusters_white_box_check_fifo_evicts_the_oldest_pair();

/** @brief Runs ThrustersWhiteBox check_expired_slots_retire_by_pair with its case context. */
void run_effect_case_thrusters_white_box_check_expired_slots_retire_by_pair();

/** @brief Runs RingShowerWhiteBox check_radius_endpoints with its case context. */
void run_effect_case_ring_shower_white_box_check_radius_endpoints();

/** @brief Runs DynamoWhiteBox check_gap_reduction_uses_shortest_path with its case context. */
void run_effect_case_dynamo_white_box_check_gap_reduction_uses_shortest_path();

/** @brief Runs DynamoWhiteBox check_overlapping_wipes_stay_in_range with its case context. */
void run_effect_case_dynamo_white_box_check_overlapping_wipes_stay_in_range();

/** @brief Runs DynamoWhiteBox check_wipe_band_colors with its case context. */
void run_effect_case_dynamo_white_box_check_wipe_band_colors();

/** @brief Runs test_ringspin_trail_hugs_its_great_circles with its case context. */
void run_effect_case_ringspin_trail_hugs_its_great_circles();

/** @brief Runs test_hopf_trail_trim_keeps_a_segment with its case context. */
void run_effect_case_hopf_trail_trim_keeps_a_segment();

/** @brief Runs test_gnomonicstars_radius_px_covers_both_axes with its case context. */
void run_effect_case_gnomonicstars_radius_px_covers_both_axes();

/** @brief Runs test_gnomonicstars_spiral_cache_invalidation with its case context. */
void run_effect_case_gnomonicstars_spiral_cache_invalidation();

/** @brief Runs test_gnomonicstars_stars_drift_from_last_draw with its case context. */
void run_effect_case_gnomonicstars_stars_drift_from_last_draw();

/** @brief Runs test_galaxies_spawn_along_spiral_with_orbital_velocity with its case context. */
void run_effect_case_galaxies_spawn_along_spiral_with_orbital_velocity();

/** @brief Runs test_galaxies_alternate_orbits_and_arm_spin with its case context. */
void run_effect_case_galaxies_alternate_orbits_and_arm_spin();

/** @brief Runs test_galaxies_arm_contrast_and_age_fade with its case context. */
void run_effect_case_galaxies_arm_contrast_and_age_fade();

/** @brief Runs test_galaxies_arm_sharpness with its case context. */
void run_effect_case_galaxies_arm_sharpness();

/** @brief Runs test_galaxies_star_color_mix with its case context. */
void run_effect_case_galaxies_star_color_mix();

/** @brief Runs test_galaxies_core_fade with its case context. */
void run_effect_case_galaxies_core_fade();

/** @brief Runs test_galaxies_white_star_stays_white with its case context. */
void run_effect_case_galaxies_white_star_stays_white();

/** @brief Runs test_galaxies_emission_rate_and_capacity with its case context. */
void run_effect_case_galaxies_emission_rate_and_capacity();

/** @brief Runs test_galaxies_live_particle_count with its case context. */
void run_effect_case_galaxies_live_particle_count();

/** @brief Runs test_column_time_scans_keep_fidelity_off_frames with its case context. */
void run_effect_case_column_time_scans_keep_fidelity_off_frames();

/** @brief Runs test_column_time_scans_engage_with_fidelity_on with its case context. */
void run_effect_case_column_time_scans_engage_with_fidelity_on();

/** @brief Runs test_galaxies_black_hole_core_toggle with its case context. */
void run_effect_case_galaxies_black_hole_core_toggle();

/** @brief Runs test_galaxies_packed_orbits_match_float with its case context. */
void run_effect_case_galaxies_packed_orbits_match_float();

/** @brief Runs test_displacement_field_lazy_hue_table_matches_eager with its case context. */
void run_effect_case_displacement_field_lazy_hue_table_matches_eager();

/** @brief Runs test_displacement_field_hue_table_bind_clears_knots with its case context. */
void run_effect_case_displacement_field_hue_table_bind_clears_knots();

/** @brief Runs test_displacement_field_ring_pool_abandoned_slot with its case context. */
void run_effect_case_displacement_field_ring_pool_abandoned_slot();

/** @brief Runs test_displacement_field_ring_pool_capacity_release with its case context. */
void run_effect_case_displacement_field_ring_pool_capacity_release();

/** @brief Runs test_displacement_field_footprint_covers_init with its case context. */
void run_effect_case_displacement_field_footprint_covers_init();

/** @brief Runs test_displacement_field_pause_holds_phase_only with its case context. */
void run_effect_case_displacement_field_pause_holds_phase_only();

/** @brief Runs test_displacement_field_default_cycle_frames with its case context. */
void run_effect_case_displacement_field_default_cycle_frames();

/** @brief Runs test_displacement_field_presets with its case context. */
void run_effect_case_displacement_field_presets();

/** @brief Runs test_displacement_field_edit_cancels_departure with its case context. */
void run_effect_case_displacement_field_edit_cancels_departure();

/** @brief Runs test_displacement_field_startup_writes_preserve_entry with its case context. */
void run_effect_case_displacement_field_startup_writes_preserve_entry();

/** @brief Runs test_displacement_field_paused_balls_deform_render with its case context. */
void run_effect_case_displacement_field_paused_balls_deform_render();

/** @brief Runs test_displacement_field_paused_balls_keep_moving with its case context. */
void run_effect_case_displacement_field_paused_balls_keep_moving();

/** @brief Runs test_displacement_field_knot_walk_unit_near_ball with its case context. */
void run_effect_case_displacement_field_knot_walk_unit_near_ball();

/** @brief Runs test_displacement_field_zero_hue_scale_is_exact with its case context. */
void run_effect_case_displacement_field_zero_hue_scale_is_exact();

/** @brief Runs test_displacement_field_octave_bake_tracks_noise with its case context. */
void run_effect_case_displacement_field_octave_bake_tracks_noise();

/** @brief Runs test_displacement_field_arm_boundaries_meet_endpoints with its case context. */
void run_effect_case_displacement_field_arm_boundaries_meet_endpoints();

/** @brief Runs test_displacement_field_single_pass_reference with its case context. */
void run_effect_case_displacement_field_single_pass_reference();

/** @brief Runs test_displacement_field_boundary_stencil with its case context. */
void run_effect_case_displacement_field_boundary_stencil();

/** @brief Runs test_displacement_field_turn_follows_large_steps with its case context. */
void run_effect_case_displacement_field_turn_follows_large_steps();

/** @brief Runs test_displacement_field_alternate_groups_match_windows with its case context. */
void run_effect_case_displacement_field_alternate_groups_match_windows();

/** @brief Runs test_glitch_lens_unit_norm with its case context. */
void run_effect_case_glitch_lens_unit_norm();

/** @brief Runs test_mobius_rings_conformal_and_counter_rotation with its case context. */
void run_effect_case_mobius_rings_conformal_and_counter_rotation();

/** @brief Runs test_voronoi_union_candidates_cover_nearest with its case context. */
void run_effect_case_voronoi_union_candidates_cover_nearest();

/** @brief Runs test_voronoi_column_time_matches_brute_force with its case context. */
void run_effect_case_voronoi_column_time_matches_brute_force();

/** @brief Runs test_sh_pullback_matches_direct_colorize with its case context. */
void run_effect_case_sh_pullback_matches_direct_colorize();

/** @brief Runs test_sh_morph_chain_rearms with its case context. */
void run_effect_case_sh_morph_chain_rearms();

/** @brief Runs test_gs_evolution_stays_bounded with its case context. */
void run_effect_case_gs_evolution_stays_bounded();

/** @brief Runs test_gs_reaction_corner_stays_bounded with its case context. */
void run_effect_case_gs_reaction_corner_stays_bounded();

/** @brief Runs test_gs_dissolve_clears_and_reseeds with its case context. */
void run_effect_case_gs_dissolve_clears_and_reseeds();

/** @brief Runs test_gs_staged_reseed_matches_synchronous with its case context. */
void run_effect_case_gs_staged_reseed_matches_synchronous();

/** @brief Runs test_dreamballs_base_mesh_selector with its case context. */
void run_effect_case_dreamballs_base_mesh_selector();

/** @brief Runs test_dreamballs_weave_topology with its case context. */
void run_effect_case_dreamballs_weave_topology();

/** @brief Runs test_dreamballs_defect_weave_renders with its case context. */
void run_effect_case_dreamballs_defect_weave_renders();

/** @brief Runs test_dreamballs_respawn_fires_and_honors_pause with its case context. */
void run_effect_case_dreamballs_respawn_fires_and_honors_pause();

/** @brief Runs test_dreamballs_select_mid_fade_in_finishes_fade with its case context. */
void run_effect_case_dreamballs_select_mid_fade_in_finishes_fade();

/** @brief Runs test_dreamballs_select_mid_fade_out_respawns with its case context. */
void run_effect_case_dreamballs_select_mid_fade_out_respawns();

/** @brief Runs test_meshfeedback_flush_precedes_mesh_draw with its case context. */
void run_effect_case_meshfeedback_flush_precedes_mesh_draw();

/** @brief Runs test_meshfeedback_preset_rotation_syncs_noise with its case context. */
void run_effect_case_meshfeedback_preset_rotation_syncs_noise();

/** @brief Runs test_comets_manual_preset_restarts_path with its case context. */
void run_effect_case_comets_manual_preset_restarts_path();

/** @brief Runs test_ash_cloud_value_cutout_gates_the_frame with its case context. */
void run_effect_case_ash_cloud_value_cutout_gates_the_frame();

/** @brief Runs test_raymarch_unit_bounds_contains_twisted_tube with its case context. */
void run_effect_case_raymarch_unit_bounds_contains_twisted_tube();

/** @brief Runs test_petalflow_spawn_gap_bounded with its case context. */
void run_effect_case_petalflow_spawn_gap_bounded();

/** @brief Runs test_petalflow_last_free_slot with its case context. */
void run_effect_case_petalflow_last_free_slot();

/** @brief Runs test_displacement_field_hue_table_fidelity with its case context. */
void run_effect_case_displacement_field_hue_table_fidelity();

/** @brief Runs test_displacement_field_hue_table_frame_fidelity with its case context. */
void run_effect_case_displacement_field_hue_table_frame_fidelity();

/** @brief Runs test_displacement_field_clip_tiles_full with its case context. */
void run_effect_case_displacement_field_clip_tiles_full();

/** @brief Runs test_displacement_field_ball_spans_and_lifecycle with its case context. */
void run_effect_case_displacement_field_ball_spans_and_lifecycle();

/** @brief Runs test_islamicstars_recipe_build_smoke with its case context. */
void run_effect_case_islamicstars_recipe_build_smoke();

/** @brief Runs test_islamicstars_roster_cycle_fits_budget with its case context. */
void run_effect_case_islamicstars_roster_cycle_fits_budget();

/** @brief Runs test_islamicstars_dual_bridge_fits_budget with its case context. */
void run_effect_case_islamicstars_dual_bridge_fits_budget();

/** @brief Runs test_islamicstars_motion_fits_scratch with its case context. */
void run_effect_case_islamicstars_motion_fits_scratch();

} // namespace hs_test::effects_tests
