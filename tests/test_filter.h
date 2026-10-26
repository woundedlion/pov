/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/render/filter.h.
 */
#pragma once

#include "tests/vec_test_util.h"

#include "core/animation/orientation.h"
#include "core/math/mobius.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <span>
#include <type_traits>

#include "core/render/filter.h"
#include "core/render/filter/feedback_warp_cache.h"
#include "core/render/filter/pixel_feedback.h"
#include "core/render/canvas.h"
#include "tests/pixel_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include <vector>

namespace hs_test {
struct FeedbackWhiteBox {
  template <int W, int H>
  static bool holds(const Filter::Pixel::Feedback<W, H> &feedback,
                    const Canvas &canvas) {
    const auto grid = feedback.make_coarse_grid(canvas);
    const auto band = feedback.make_render_band(canvas.clip(), grid);
    return feedback.warp_cache.holds(feedback.warp_cache_key(band));
  }
  template <int W, int H>
  static size_t fills(const Filter::Pixel::Feedback<W, H> &feedback) {
    return feedback.warp_cache.cached_fills;
  }
};

namespace filter_tests {

static_assert(hs::H_OFFSET == 0,
              "test_filter.h assumes row H-1 is the south pole; the legacy "
              "offset-3 mapping belongs in tests/test_h_offset_renorm.h");

template <typename P>
concept RawFramePlotter = requires(P &pipeline, Canvas &canvas) {
  pipeline.plot(canvas, 0, 0, Pixel{}, 0.0f, 1.0f);
};

template <typename P>
concept TerminalFlusher =
    requires(P &pipeline, Canvas &canvas) { pipeline.flush(canvas, 1.0f); };

template <typename P>
concept ReplacementFrameStarter = requires(P &pipeline, Canvas &canvas) {
  pipeline.begin_frame(canvas, 1.0f);
};

#include "tests/filter/traits.h"
#include "tests/filter/trait_inheritance.h"
#include "tests/filter/pipeline_sink.h"
#include "tests/filter/antialias.h"
#include "tests/filter/blur.h"
#include "tests/filter/chromatic_shift.h"
#include "tests/filter/feedback.h"
#include "tests/filter/world_filters.h"
#include "tests/filter/temporal_capture.h"
#include "tests/filter/temporal_copies.h"
#include "tests/filter/canvas_routing.h"
#include "tests/filter/feedback_flush.h"
#include "tests/filter/feedback_column_time.h"
#include "tests/filter/trails.h"
#include "tests/filter/segmented_bounds.h"
/**
 * @brief Runs every filter test case under the "filter" module scope.
 * @return The module's failure count.
 */
inline int run_filter_tests() {
  hs_test::ModuleFixture fixture("filter");

  HS_RUN(test_trait_member_values);
  HS_RUN(test_filter_trait_inheritance);
  HS_RUN(test_crosses_segments_trait_and_fold);
  HS_RUN(test_history_domain_folds);

  HS_RUN(test_replacing_terminal_without_history_flushes);
  HS_RUN(test_pipeline_sink_and_terminal_concept_surface);
  HS_RUN(test_pipeline_get_returns_correct_filter);

  HS_RUN(test_antialias_weights_partition);
  HS_RUN(test_antialias_integer_coord_single_tap);
  HS_RUN(test_antialias_seam_wraps_left_column);
  HS_RUN(test_antialias_far_seam_wraps_both_taps);
  HS_RUN(test_antialias_clips_virtual_subpole_row);
  HS_RUN(test_antialias_folds_north_boundary_row);

  HS_RUN(test_blur_factor_zero_is_identity);
  HS_RUN(test_blur_folds_boundary_and_drops_virtual_rows);
  HS_RUN(test_blur_full_kernel_sums_to_alpha);
  HS_RUN(test_blur_kernel_weights_by_offset);
  HS_RUN(test_blur_update_changes_kernel);
  HS_RUN(test_blur_clamps_factor);
  HS_RUN(test_blur_wraps_column_taps);
  HS_RUN(test_blur_pole_row_renormalizes);

  HS_RUN(test_chromatic_shift_fanout);

  HS_RUN(test_feedback_style_binding);
  HS_RUN(test_feedback_plot_is_passthrough);
  HS_RUN(test_feedback_warp_cache_valid_only_after_fill);

  HS_RUN(test_world_hole_masks_cap);
  HS_RUN(test_world_hole_setters);
  HS_RUN(test_world_hole_radius_limits);
  HS_RUN(test_world_hole_angular_reference);
  HS_RUN(test_world_orient_rotates_and_keeps_static_age);
  HS_RUN(test_world_orient_emits_at_capture_time);
  HS_RUN(test_world_orient_explicit_time);
  HS_RUN(test_world_orient_follows_history_at_column_times);
  HS_RUN(test_world_orient_sample_cache);
  HS_RUN(test_timed_probe_bypasses_history);
  HS_RUN(test_temporal_column_arithmetic_matches_remainder);
  HS_RUN(test_temporal_orientation_bound_is_conservative);
  HS_RUN(test_timed_boundary_tap_ownership);
  HS_RUN(test_timed_cull_region);
  HS_RUN(test_timed_replicate_selected_copy);
  HS_RUN(test_timed_vertex_replicate_selected_copy);
  HS_RUN(test_temporal_capture_matches_exhaustive_oracle);
  HS_RUN(test_temporal_capture_column_matches_composed_time);
  HS_RUN(test_temporal_capture_tilted_axis_within_half_column);
  HS_RUN(test_temporal_capture_unsolvable_motion_falls_back_to_endpoint);
  HS_RUN(test_sweep_slabs_tile_the_half_turn);
  HS_RUN(test_world_orient_cull_edge_mirrors_plot);
  HS_RUN(test_world_orient_slice_selects_by_projection);
  HS_RUN(test_world_orient_slice_cull_edge_bounds_all_slices);
  HS_RUN(test_world_vertex_replicate_fanout_and_age);
  HS_RUN(test_world_vertex_replicate_cull_edge_mirrors_plot);
  HS_RUN(test_replicate_three_planar_basis_copies);
  HS_RUN(test_pipeline_could_intersect_clip_forwards_through_stages);
  HS_RUN(test_world_mobius_identity_and_transform);

  HS_RUN(test_pipeline_sink_2d_plot_writes_wraps_clips);
  HS_RUN(test_pipeline_composition_alpha_and_draw_order);
  HS_RUN(test_pipeline_sink_3d_plot_routes_to_canvas);
  HS_RUN(test_pipeline_world_replicate_fans_out);
  HS_RUN(test_pipeline_2d_into_3d_head_roundtrips);
  HS_RUN(test_pipeline_screen_antialias_routes_to_sink);
  HS_RUN(test_direct_antialias_sink_framebuffer_parity);
  HS_RUN(test_direct_antialias_sink_stale_clip);
  HS_RUN(test_feedback_flush_blends_prev_frame);
  HS_RUN(test_feedback_north_pole_uses_single_physical_sample);
  HS_RUN(test_feedback_south_pole_uses_single_physical_sample);
  HS_RUN(test_feedback_polar_rows_use_spherical_footprint);
  HS_RUN(test_feedback_flush_respects_clip);
  HS_RUN(test_feedback_flush_melt_warp_displaces_south);
  HS_RUN(test_feedback_north_cap_uses_exact_control_rows);
  HS_RUN(test_feedback_animated_cap_controls_match_compositor_lattice);
  HS_RUN(test_feedback_poles_resolve_one_source_longitude);
  HS_RUN(test_feedback_half_res_warp_uses_pair_midpoint);
  HS_RUN(test_feedback_spherical_ring_control_rows);
  HS_RUN(test_feedback_init_storage_matches_storage_bytes);
  HS_RUN(test_spherical_field_layout_angular_error);
  HS_RUN(test_feedback_seam_warp_keeps_its_latitude_row);
  HS_RUN(test_feedback_cached_north_cap_clips_share_control_rows);
  HS_RUN(test_feedback_warp_cache_matches_uncached);
  HS_RUN(test_feedback_polar_rows_hit_their_targets);
  HS_RUN(test_feedback_flush_straddled_taps_stay_on_branch);
  HS_RUN(test_feedback_column_time_reads);
  HS_RUN(test_feedback_fade_table_entries);
  HS_RUN(test_feedback_warp_cache_keys_column_time);
  HS_RUN(test_feedback_column_time_oracle);

  HS_RUN(test_world_trails_int16_quantization_roundtrip);
  HS_RUN(test_world_trails_clamps_out_of_range);
  HS_RUN(test_world_trails_capacity_evicts_one_slot);
  HS_RUN(test_world_trails_ttl_expiry);
  HS_RUN(test_world_trails_set_lifetime_caps_ttl);
  HS_RUN(test_world_trails_midbuffer_expiry_reclaims_slot);
  HS_RUN(test_screen_trails_store_emit_decay);
  HS_RUN(test_screen_trails_negative_age_clamps_t);
  HS_RUN(test_screen_trails_forwards_aged_emission);
  HS_RUN(test_screen_trails_at_capacity_replaces_last_slot);
  HS_RUN(test_screen_trails_set_lifetime_caps_ttl);
  HS_RUN(test_trails_alpha_gates);
  HS_RUN(test_mixed_domain_flush_drains_both_buffers);

  HS_RUN(test_effect_needs_full_frame_default_false);
  HS_RUN(test_screen_trails_banded_matches_full);
  HS_RUN(test_screen_trails_alternating_clip_matches_full);
  HS_RUN(test_feedback_banded_diverges_from_full);

  return fixture.result();
}

} // namespace filter_tests
} // namespace hs_test
