/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/animation/animation.h.
 */
#pragma once

#include "core/animation/orientation.h"
#include "core/math/mobius.h"
#include <array>
#include <cstring>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>
#include "core/animation/animation.h"
#include "core/render/canvas.h"
#include "core/math/easing.h"
#include "core/mesh/mesh.h" // PolyMesh, MeshOps::compile (mesh test fixtures)
#include "tests/fd_capture_util.h"
#include "tests/mesh_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/vec_test_util.h"

namespace hs_test {
namespace animation_tests {

#include "tests/animation/canvas_fixture.h"
#include "tests/animation/paths.h"
#include "tests/animation/transitions.h"
#include "tests/animation/mutation.h"
#include "tests/animation/driver.h"
#include "tests/animation/lerp.h"
#include "tests/animation/orientation_trail.h"
#include "tests/animation/borrow_contracts.h"
#include "tests/animation/rotation_cases.h"
#include "tests/animation/timeline_cases.h"
#include "tests/animation/orientation_cases.h"
#include "tests/animation/motion_cases.h"
#include "tests/animation/particles.h"
#include "tests/animation/sprite_lifecycle.h"
#include "tests/animation/mesh_carousel.h"
#include "tests/animation/deep_tween.h"
#include "tests/animation/arena_compaction.h"
#include "tests/animation/color_wipe.h"
#include "tests/animation/mobius.h"
#include "tests/animation/ripple_noise.h"
#include "tests/animation/public_api.h"
/**
 * @brief Runs every animation test case in this module.
 * @return The module's failure count.
 */
inline int run_animation_tests() {
  hs_test::ModuleFixture fixture("animation");

  // Module-scoped fake-canvas fixture.
  hs_test::StubEffect fake_fx(8, 8);
  Canvas fake_cv(fake_fx);
  fake_canvas_ptr() = &fake_cv;

  test_progress_pause_and_eased_bounds();
  test_trail_body_records_independent_orientation_history();

  test_path_adjacent_segments_fill_exact_capacity();
  test_path_empty_returns_origin();
  test_path_endpoints_and_clamp();
  test_path_collapse_keeps_last();

  test_transition_reaches_target_linear();
  test_transition_monotonic_and_starts_from_current();
  test_transition_duration_zero_no_divide_by_zero();
  test_transition_quantized_floors_result();
  test_transition_repeat_retraverses_each_cycle();
  test_transition_paused_holds_value();

  test_mutation_applies_function_of_eased_time();
  test_mutation_duration_zero_finite();

  test_driver_increments_and_wraps();
  test_driver_no_wrap_accumulates();
  test_driver_set_speed_ignores_non_finite();
  test_driver_nan_source_does_not_poison();

  test_lerp_drives_subject_to_target();
  test_lerp_midpoint();

  test_orientation_trail_index_zero_is_oldest();
  test_orientation_trail_expire_drops_oldest();
  test_orientation_trail_clear();

  test_rotation_substeps_shared_and_tight();
  test_rotation_accumulates_subthreshold_deltas();
  test_rotation_applies_final_frame_residual();
  test_timeline_shared_orientation_composes_motion_blur();
  test_timeline_collapse_past_id_cache();
  test_timeline_sequences_events_by_start_frame();
  test_timeline_pausable_event_uses_active_time();
  test_timeline_accepts_maximum_start_frame();
  test_timeline_rollover_preserves_eligibility();
  test_timeline_paused_rollover_preserves_active_time();
  test_timers_rollover_preserves_intervals();
  test_timeline_repeating_animation_rewinds_each_cycle();
  test_timeline_cancel_removes_repeating_animation();
  test_timeline_cancel_suppresses_step_side_effects();
  test_timeline_paused_event_redraws_from_start_frame();
  test_timeline_repeating_canceled_in_callback_fires_then_once();
  test_timeline_cancel_fires_post_callback();
  test_timeline_cancel_while_paused_removes_event();
  test_timeline_cancel_before_start();
  test_timeline_compaction_preserves_later_events();
  test_timeline_then_chains_follow_up_event();
  test_repeating_timer_fires_then_each_cycle();
  test_repeating_timer_canceled_in_callback_fires_then_once();
  test_timer_then_self_cancellation_completes_once();
  test_timeline_clear_destroys_events_keeping_frame();
  test_timeline_instance_boundary_reclaims_pinned_event();
  test_timeline_full_guard_rejects_overflow();
  test_timeline_remove_clear_hook_unregisters_by_ctx();
  test_orientation_upsample_then_collapse();
  test_motion_repeating_does_not_drift();
  test_motion_reanchor_after_path_swap();
  test_motion_codriven_survives_repeat_seam();

  test_particle_system_spawn_and_capacity_guard();
  test_particle_system_lifetime_boundaries();
  test_point_particle_storage_and_lifetime();
  test_particle_system_spawn_initializes_and_steps();
  test_particle_system_sparse_trail_sampling();
  test_particle_system_reclaims_at_life_expiry();
  test_particle_system_attractor_kills_within_radius();
  test_particle_system_attractor_kill_radius_boundary();
  test_particle_system_signed_axis_one_step_equivalence();
  test_particle_system_signed_axis_boundaries();
  test_particle_system_signed_axis_trajectory();

  test_sprite_fade_in_plateau_fade_out_envelope();
  test_sprite_clamps_overshooting_fade_in();
  test_sprite_overlapping_fades_stay_continuous();
  test_sprite_paused_holds_frame();
  test_timeline_pause_redraws_held_sprite();
  test_sprite_paused_before_first_step_holds_first_opacity();

  test_crossfade_segue_schedules_overlapping_sprite();
  test_crossfade_segue_clamps_fade_to_half_duration();
  test_crossfade_segue_overlap_is_configurable();
  test_sequential_segue_never_overlaps_sprites();
  test_dissolve_segue_masks_partition_keys();
  test_dissolve_segue_reseeds_per_frame_and_transition();
  test_dissolve_segue_overlaps_the_full_fade_window();
  test_segue_base_hooks_are_identity();
  test_segue_visible_gate_culls_only_dark_phases();
  test_iris_bloom_fill_contracts_to_face_centers();
  test_lace_fill_keeps_edge_band();
  test_sweep_phase_front_ordering();
  test_meshcarousel_face_phases_use_sweep_frame_and_slots();
  test_terminator_sweep_orders_by_axis();
  test_terminator_sweep_fades_faces_over_fixed_frames();
  test_terminator_sweep_fade_sliders_apply_without_reschedule();
  test_terminator_sweep_per_face_fade_random_in_range();
  test_shockwave_orders_by_distance_from_origin();
  test_per_face_segues_satisfy_draw_contract();
  test_segue_policies_forward_pause_gate();
  test_breakdown_fades_classes_sequentially();
  test_breakdown_guards_degenerate_class_inputs();
  test_spin_flip_warp_is_rigid();
  test_gold_convergence_grades_to_gold();

  test_tweenable_rejects_bare_orientation();
  test_deep_tween_global_t_spans_unit_interval();
  test_tween_orientation_skips_shared_boundary();
  test_deep_tween_collapsed_newest_frame_reaches_one();
  test_deep_tween_all_collapsed_reaches_one();
  test_deep_tween_frames_groups_flat_emission();
  test_deep_tween_interior_motionless_frame_no_gap();
  test_deep_tween_oldest_motionless_frame_no_gap();
  test_tween_vectortrail_single_sample_reaches_one();
  test_quantized_vector_trail_roundtrip_and_ring();

  test_meshcarousel_compact_keep_front_drops_back();
  test_meshcarousel_compact_drop_all_frees_both_slots();

  test_colorwipe_reaches_target_keys();
  test_colorwipe_uses_owned_start_snapshot();
  test_colorwipe_slow_fade_resolves_every_frame();
  test_colorwipe_paused_holds_keys();

  test_mobiuswarp_closes_at_completion();
  test_mobiuswarp_bind_scale_reads_live();
  test_mobiuswarp_retains_last_finite_scale();
  test_mobiuswarp_circular_traces_radius();
  test_mobiuswarp_circular_bind_scale_reads_live();
  test_mobiuswarp_evolving_bounded_and_perpetual();
  test_mobiuswarp_evolving_wrapped_live_channels();
  test_mobiuswarp_evolving_value_semantics();
  test_mobiuswarp_evolving_long_uptime();

  test_ripple_envelope_and_done_boundary();
  test_noise_publishes_time_and_is_perpetual();
  test_noise_loop_coordinates_and_samples<Animation::NoiseParams,
                                          Animation::Noise>();
  test_noise_loop_coordinates_and_samples<Animation::NoiseProductParams,
                                          Animation::NoiseProduct>();
  test_noise_loop_live_controls_and_value_semantics<Animation::NoiseParams,
                                                    Animation::Noise>();
  test_noise_loop_live_controls_and_value_semantics<
      Animation::NoiseProductParams, Animation::NoiseProduct>();
  test_noise_loop_long_uptime<Animation::NoiseParams, Animation::Noise>();
  test_noise_loop_long_uptime<Animation::NoiseProductParams,
                              Animation::NoiseProduct>();
  test_noise_loop_preserves_finite_duration();

  test_random_walk_stays_unit_and_travels();
  test_random_walk_stable_rotation_matches_same_state();

  test_random_timer_fires_within_range();
  test_one_shot_timer_ends_by_completion_not_cancel();
  test_finish_terminates_a_repeating_animation();
  test_finished_param_animation_progress_is_finite();
  test_periodic_timer_set_period_reschedules_from_now();
  test_periodic_timer_set_period_unchanged_does_not_defer();
  test_mobiusflow_degenerate_inputs_remain_finite();
  test_mobiusflow_step_preserves_unit_product();
  test_particle_system_emitter_dispatch();
  test_motion_set_duration_below_position_rescales();

  const int result = fixture.result();
  // Unpublish before fake_cv/fake_fx destruct.
  fake_canvas_ptr() = nullptr;
  return result;
}

} // namespace animation_tests
} // namespace hs_test
