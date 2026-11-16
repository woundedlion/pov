/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Canvas test-module entry. */

#include "tests/canvas/effect_transition.h"
#include "tests/canvas/sparse_parameter_specs.h"
#include "tests/canvas/base_cases.h"

namespace hs_test::canvas_tests {

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every canvas test in order.
 * @return The module's failure count reported by end_module().
 */
inline int run_canvas_tests() {
  hs_test::ModuleFixture fixture("canvas");

  HS_RUN(test_construction_dims_and_clear);
  HS_RUN(test_effect_config_strobe_override);
  HS_RUN(test_construction_dimension_boundaries);
  HS_RUN(test_output_envelope_endpoints);
  HS_RUN(test_effect_transition_fenced_commit);
  HS_RUN(test_effect_transition_fences_hold);
  HS_RUN(test_effect_transition_replacement_pause_and_restore);
  HS_RUN(test_effect_transition_restore_none_enters_failsafe);
  HS_RUN(test_effect_transition_shorter_replacement_preserves_progress);
  HS_RUN(test_effect_transition_failsafe_retry);
  HS_RUN(test_effect_transition_refusal_branches);
  HS_RUN(test_preset_state_machine);
  HS_RUN(test_frame_visible_only_after_advance_display);
  HS_RUN(test_non_persist_frame_starts_cleared);
  HS_RUN(test_clip_clear_exact_rectangle);
  HS_RUN(test_clip_clear_alternates_clips_and_buffers);
  HS_RUN(test_clip_clear_full_clip_matches_full_clear);
  HS_RUN(test_clip_clear_full_width_band);
  HS_RUN(test_clip_clear_constructor_zeroes_both_buffers);
  HS_RUN(test_clip_clear_does_not_change_persistence);
  HS_RUN(test_clip_clear_full_frame_without_reads_keeps_margin);
  HS_RUN(test_clip_clear_outside_band_reader_clears_full_buffer);
  HS_RUN(test_persist_pixels_copies_previous_frame);
  HS_RUN(test_double_buffer_handoff_no_aliasing);
  HS_RUN(test_double_buffer_handoff_concurrent);
  HS_RUN(test_ctor_spin_waits_for_buffer_free);
  HS_RUN(test_canvas_2d_and_1d_access_and_prev);
  HS_RUN(test_parameter_group_registered_companion);
  HS_RUN(test_register_float_and_bool_params);
  HS_RUN(test_update_parameter_by_name);
#if HS_ENABLE_PARAM_GUI_BRIDGE
  HS_RUN(test_parameter_display_mirror);
#endif
  HS_RUN(test_update_parameter_rejects_readonly);
  HS_RUN(test_register_and_update_enum_param);
  HS_RUN(test_typed_enum_and_global_param_metadata);
  HS_RUN(test_typed_enum_storage_widths);
  HS_RUN(test_typed_parameter_specs);
  HS_RUN(test_sparse_parameter_fields);
  HS_RUN(test_sparse_parameter_schema_rejections);
  HS_RUN(test_sparse_parameter_factory_and_float_fields);
  HS_RUN(test_sparse_parameter_updates);
  HS_RUN(test_sparse_integer_registration_metadata);
#if HS_ENABLE_PARAM_GUI_BRIDGE
  HS_RUN(test_parameter_spec_preserves_requested_float);
#endif
  HS_RUN(test_paramlist_fills_to_capacity);
  HS_RUN(test_paramlist_schema_generation);
  HS_RUN(test_clip_setters);
  HS_RUN(test_clip_content_key_is_exact);
  HS_RUN(test_canvas_clip_key_tracks_mutators);
  HS_RUN(test_canvas_cached_render_clip_matches_region);
  HS_RUN(test_effect_config_margin);
  HS_RUN(test_pipeline_ref_routes_screen_coordinate_overloads);
  HS_RUN(test_pipeline_ref_requires_plot_overloads);
  HS_RUN(test_pipeline_ref_copy_is_independent_of_source);

  return fixture.result();
}
} // namespace hs_test::canvas_tests
