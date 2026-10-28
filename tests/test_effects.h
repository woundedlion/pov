/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Effects white-box suite and compatibility aggregate. */
#pragma once

#include "targets/effects.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test {
namespace effects_tests {

#include "tests/effects/parameter_probe.h"
#include "tests/effects/spherical_harmonics.h"
#include "tests/effects/reaction_diffusion_gs.h"
#include "tests/effects/reaction_diffusion_bz.h"
#include "tests/effects/hankin.h"
#include "tests/effects/dreamballs_cases.h"
#include "tests/effects/mesh_feedback.h"
#include "tests/effects/comets_cases.h"
#include "tests/effects/ash_cloud.h"
#include "tests/effects/ring_whitebox.h"
#include "tests/effects/dynamo_cases.h"
#include "tests/effects/hopf.h"
#include "tests/effects/rings_temporal.h"
#include "tests/effects/dynamo_temporal.h"
#include "tests/effects/raymarch_cases.h"
#include "tests/effects/gnomonic_stars.h"
#include "tests/effects/fishbowl_cases.h"
#include "tests/effects/ring_spin.h"
#include "tests/effects/petal_flow.h"
#include "tests/effects/lens_cases.h"
#include "tests/effects/mobius_rings.h"
#include "tests/effects/roster_cases.h"
#include "tests/effects/petal_flow_oracle.h"
#include "tests/effects/fishbowl_oracle.h"
#include "tests/effects/galaxies_tests.h"
#include "tests/effects/displacement_field.h"
#include "tests/effects/shape_shifter.h"
#include "tests/effects/voronoi_cases.h"
#include "tests/effects/islamic_stars.h"
#include "tests/effects/column_time_parity.h"
#include "tests/effects/mesh_capture_cases.h"

/**
 * @brief Module entry point for the effects white-box suite.
 * @return Module result code from hs_test::end_module (0 on success).
 * @details Per-effect invariants with an explicit oracle; excluded from the
 * fast-math axis.
 */
inline int run_effects_tests() {
  hs_test::ModuleFixture fixture("effects");

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  const char *mobius_only = std::getenv("HS_MOBIUSRINGS_ONLY");
  const char *islamic_only = std::getenv("HS_ISLAMICSTARS_ONLY");
  const char *mesh_only = std::getenv("HS_MESH_CAPTURE_ONLY");
  const char *rings_only = std::getenv("HS_RING_TEMPORAL_ONLY");
  const char *dynamo_only = std::getenv("HS_DYNAMO_ONLY");
#pragma clang diagnostic pop
  if (rings_only && std::atoi(rings_only) != 0) {
    run_effect_case_ring_shower_white_box_check_temporal_rings();
    run_effect_case_thrusters_white_box_check_temporal_rings();
    return fixture.result();
  }
  if (mobius_only && std::atoi(mobius_only) != 0) {
    run_effect_case_mobius_rings_white_box_check_column_time_oracle();
    return fixture.result();
  }
  if (dynamo_only && std::atoi(dynamo_only) != 0) {
    run_effect_case_dynamo_column_time_oracle();
    run_effect_case_dynamo_column_time_first_frame_and_reset();
    if (std::atoi(dynamo_only) == 1) {
      run_effect_case_dynamo_emitted_points_counts_ring_seeds();
      run_effect_case_dynamo_trail_ceiling_bounds_the_ring();
    }
    return fixture.result();
  }
  if (mesh_only && std::atoi(mesh_only) == 3) {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    const char *preset = std::getenv("HS_REPLAY_PRESET");
    const char *speed = std::getenv("HS_REPLAY_TRANS_SPEED");
    const char *frames = std::getenv("HS_REPLAY_FRAMES");
#pragma clang diagnostic pop
    replay_islamicstars_profile(
        preset ? std::atoi(preset) : 6, speed ? std::atof(speed) : 4.0f,
        frames ? std::atoi(frames) : PROFILE_REPLAY_FRAMES);
    return fixture.result();
  }
  if (mesh_only && std::atoi(mesh_only) == 2) {
    run_effect_case_hankinsolids_handoff_history();
    measure_mesh_capture_work();
    return fixture.result();
  }
  if (mesh_only && std::atoi(mesh_only) != 0) {
    run_effect_case_hankinsolids_interlace_motion_is_seam_free();
    run_effect_case_hankinsolids_seam_copies_are_one_polygon();
    run_effect_case_hankinsolids_morph_motion_is_seam_free();
    run_effect_case_islamicstars_ripple_motion_is_seam_free();
    run_effect_case_islamicstars_seam_copies_are_one_polygon();
    run_effect_case_islamicstars_playback_motion_is_seam_free();
    return fixture.result();
  }
  if (islamic_only && std::atoi(islamic_only) != 0) {
    run_effect_case_islamicstars_cases();
    return fixture.result();
  }

  run_effect_case_hopf_temporal_trail_oracle();
  run_effect_case_ring_shower_white_box_check_temporal_rings();
  run_effect_case_thrusters_white_box_check_temporal_rings();
  run_effect_case_dreamballs_temporal_source_oracle();
  run_effect_case_dreamballs_column_work_is_bounded();
  run_effect_case_dream_balls_white_box_check_endpoint_error_bound();
  run_effect_case_petalflow_temporal_source_oracle();
  run_effect_case_meshfeedback_temporal_source_oracle();
  run_effect_case_meshfeedback_stationary_timed_source();
  run_effect_case_meshfeedback_trails_arm_cut_seam();
  run_effect_case_parameter_probe_targets();
  run_effect_case_float_dropdown_lint();
  run_effect_case_ringspin_strobe_configuration_preserves_rendering();
  run_effect_case_meshfeedback_base_mesh_selector();
  run_effect_case_meshfeedback_preset_export_arity();
  run_effect_case_meshfeedback_mesh_rebuild_reuses_storage();
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
  run_effect_case_fishbowl_preset_and_fire_duty_cycle();
  run_effect_case_fishbowl_temporal_trail_oracle();
  run_effect_case_sh_preset_mode_mapping();
  run_effect_case_sh_manual_preset_replaces_inflight_morph();
  run_effect_case_shapeshifter_preset_defaults();
  run_effect_case_shapeshifter_slider_selections_render();
  run_effect_case_manual_preset_navigation();
  run_effect_case_hankinsolids_manual_pause_holds_morph();
  run_effect_case_hankinsolids_handoffs_shade_at_column_time();
  run_effect_case_hankinsolids_handoff_history();
  run_effect_case_hankinsolids_conway_crossfade_shades_at_column_time();
  run_effect_case_hankinsolids_interlace_motion_is_seam_free();
  run_effect_case_hankinsolids_seam_copies_are_one_polygon();
  run_effect_case_hankinsolids_morph_draws_carry_previous_vertices();
  run_effect_case_hankinsolids_morph_motion_is_seam_free();
  // Both tiers: arena budgets.
  run_effect_case_fishbowl_scratch_estimate_covers_peak();
  run_effect_case_hankinsolids_arena_budget_covers_every_solid();
  run_effect_case_dreamballs_max_edge_solid_render();
  // Both tiers: Raymarch geometry and presets.
  run_effect_case_raymarch_volume_random_walks_are_independent();
  run_effect_case_raymarch_preset_and_placement_solids();
  run_effect_case_raymarch_surface_frame_uv();
  run_effect_case_raymarch_lighting_and_palette_at_column_time();
  run_effect_case_raymarch_hue_noise_at_column_time();
  // Resolution-independent white-box math.
  run_effect_case_gs_hot_flags_match_directed_graph();
  run_effect_case_gs_pigment_codec_layout();
  run_effect_case_gs_rest_state_is_fixed_point();
  run_effect_case_gs_substep_signs_and_clamp();
  run_effect_case_bz_q16_roundtrip();
  run_effect_case_bz_advance_species_signs_and_clamp();
  run_effect_case_bz_perturb_state_draw_count_pinned();
  run_effect_case_hopf_projection_math();
  run_effect_case_raymarch_constexpr_sqrt_converges();
  run_effect_case_displacement_field_retained_history_resets();
  run_effect_case_dreamballs_retained_history_resets();
  // Both tiers: additional white-box cases.
  run_effect_case_needs_full_frame_gate();
  run_effect_case_voronoi_axes_use_uniform_sampler();
  run_effect_case_voronoi_segment_render_matches_full_frame();
  run_effect_case_voronoi_instantaneous_render_pins();
  run_effect_case_voronoi_endpoint_pairs_match_brute_force();
  run_effect_case_voronoi_temporal_toggle();
  run_effect_case_voronoi_reseed_resets_history();
  run_effect_case_voronoi_retained_corners_match_search();
  run_effect_case_sh_field_blend_endpoints_and_frame();
  run_effect_case_sh_field_stays_inside_unit_range();
  run_effect_case_sh_polarity_split_and_ao_shaping();
  run_effect_case_sh_cycle_covers_every_mode();
  run_effect_case_sh_instantaneous_pixel_golden();
  run_effect_case_sh_columns_follow_sweep_time();
  run_effect_case_sh_segmented_columns_match_single_engine();
  run_effect_case_sh_manual_preset_resets_morph_history();
  run_effect_case_sh_fidelity_toggle_keeps_history();
  run_effect_case_sh_frame_sampler_rollover_and_sign();
  run_effect_case_gs_reaction_edit_starts_dissolve();
  run_effect_case_gs_stall_hold_starts_dissolve();
  run_effect_case_gs_preset_pause_does_not_hold_lifecycle();
  run_effect_case_gs_dissolve_absorbs_edits_and_stages_reseed();
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
  run_effect_case_gs_untimed_frames_are_pinned();
  run_effect_case_gs_lattice_basis_matches_rotated_lattice();
  run_effect_case_gs_slabs_shade_the_field_after_their_substep();
  run_effect_case_gs_slab_noise_keeps_the_previous_bake();
  run_effect_case_gs_segmented_slabs_match_full_frame();
  run_effect_case_dreamballs_preset_cycle_bookkeeping();
  run_effect_case_comets_white_box_check_paths_close();
  run_effect_case_comets_white_box_check_births_skip_earlier_columns();
  run_effect_case_comets_rollover_skipped_mid_wipe();
  run_effect_case_comets_wipe_shades_at_column_time();
  run_effect_case_comets_wipe_history();
  run_effect_case_thrusters_white_box_check_warp_endpoints();
  run_effect_case_thrusters_white_box_check_fire_spawns_opposed_pair();
  run_effect_case_thrusters_white_box_check_collapsed_ring_falls_back_to_a_spin_axis();
  run_effect_case_thrusters_white_box_check_fifo_evicts_the_oldest_pair();
  run_effect_case_thrusters_white_box_check_expired_slots_retire_by_pair();
  run_effect_case_ring_shower_white_box_check_radius_endpoints();
  run_effect_case_dynamo_column_time_oracle();
  run_effect_case_dynamo_column_time_first_frame_and_reset();
  run_effect_case_dynamo_white_box_check_gap_reduction_uses_shortest_path();
  run_effect_case_dynamo_white_box_check_overlapping_wipes_stay_in_range();
  run_effect_case_dynamo_white_box_check_wipe_band_colors();
  run_effect_case_dynamo_emitted_points_counts_ring_seeds();
  run_effect_case_ringspin_trail_hugs_its_great_circles();
  run_effect_case_hopf_trail_trim_keeps_a_segment();
  run_effect_case_gnomonicstars_radius_px_covers_both_axes();
  run_effect_case_gnomonicstars_spiral_cache_invalidation();
  run_effect_case_gnomonicstars_stars_drift_from_last_draw();
  run_effect_case_galaxies_spawn_along_spiral_with_orbital_velocity();
  run_effect_case_galaxies_alternate_orbits_and_arm_spin();
  run_effect_case_galaxies_arm_contrast_and_age_fade();
  run_effect_case_galaxies_arm_sharpness();
  run_effect_case_galaxies_star_color_mix();
  run_effect_case_galaxies_core_fade();
  run_effect_case_galaxies_white_star_stays_white();
  run_effect_case_galaxies_emission_rate_and_capacity();
  run_effect_case_galaxies_live_particle_count();
#if defined(_WIN32)
  run_effect_case_column_time_scans_keep_fidelity_off_frames();
#endif
  run_effect_case_column_time_scans_engage_with_fidelity_on();
  run_effect_case_galaxies_black_hole_core_toggle();
  run_effect_case_galaxies_packed_orbits_match_float();
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
  run_effect_case_glitch_lens_unit_norm();
  run_effect_case_mobius_rings_conformal_and_counter_rotation();
  run_effect_case_mobius_rings_white_box_check_column_time_oracle();
  run_effect_case_islamicstars_cases();
  run_effect_case_ash_cloud_value_cutout_gates_the_frame();

  // FULL tier only (HS_EFFECTS_FULL=1).
  if (effects_full_suite()) {
    run_effect_case_voronoi_union_candidates_cover_nearest();
    run_effect_case_voronoi_column_time_matches_brute_force();
    run_effect_case_sh_pullback_matches_direct_colorize();
    run_effect_case_sh_morph_chain_rearms();
    run_effect_case_gs_evolution_stays_bounded();
    run_effect_case_gs_reaction_corner_stays_bounded();
    run_effect_case_gs_dissolve_clears_and_reseeds();
    run_effect_case_gs_staged_reseed_matches_synchronous();
    run_effect_case_dreamballs_base_mesh_selector();
    run_effect_case_dreamballs_weave_topology();
    run_effect_case_dreamballs_defect_weave_renders();
    run_effect_case_dreamballs_respawn_fires_and_honors_pause();
    run_effect_case_dreamballs_select_mid_fade_in_finishes_fade();
    run_effect_case_dreamballs_select_mid_fade_out_respawns();
    run_effect_case_meshfeedback_flush_precedes_mesh_draw();
    run_effect_case_meshfeedback_preset_rotation_syncs_noise();
    run_effect_case_comets_manual_preset_restarts_path();
    run_effect_case_dynamo_trail_ceiling_bounds_the_ring();
    run_effect_case_raymarch_unit_bounds_contains_twisted_tube();
    run_effect_case_petalflow_spawn_gap_bounded();
    run_effect_case_petalflow_last_free_slot();
    run_effect_case_displacement_field_hue_table_fidelity();
    run_effect_case_displacement_field_hue_table_frame_fidelity();
    run_effect_case_displacement_field_clip_tiles_full();
    run_effect_case_displacement_field_ball_spans_and_lifecycle();
    run_effect_case_islamicstars_recipe_build_smoke();
    run_effect_case_islamicstars_roster_cycle_fits_budget();
    run_effect_case_islamicstars_dual_bridge_fits_budget();
    run_effect_case_islamicstars_motion_fits_scratch();
  } else {
    std::printf("  [TIER] full effects cases omitted; set HS_EFFECTS_FULL=1\n");
  }

  return fixture.result();
}

} // namespace effects_tests
} // namespace hs_test
