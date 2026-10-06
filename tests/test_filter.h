/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/render/filter.h.
 */
#pragma once

#include "core/animation/orientation.h"
#include "core/math/mobius.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <span>

#include "core/render/filter.h"
#include "core/render/filter/pixel_feedback.h"
#include "core/render/canvas.h"
#include "tests/pixel_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include <vector>

namespace hs_test {
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
#include "tests/filter/canvas_routing.h"
#include "tests/filter/trails.h"
#include "tests/filter/segmented_bounds.h"
/**
 * @brief Runs every filter test case under the "filter" module scope.
 * @return The module's failure count.
 */
inline int run_filter_tests() {
  hs_test::ModuleFixture fixture("filter");

  test_trait_member_values();
  test_filter_trait_inheritance();
  test_crosses_segments_trait_and_fold();
  test_history_domain_folds();

  test_replacing_terminal_without_history_flushes();
  test_pipeline_sink_is_2d();
  test_pipeline_get_returns_correct_filter();

  test_antialias_weights_partition();
  test_antialias_integer_coord_single_tap();
  test_antialias_seam_wraps_left_column();
  test_antialias_far_seam_wraps_both_taps();
  test_antialias_clips_virtual_subpole_row();

  test_blur_factor_zero_is_identity();
  test_blur_folds_boundary_and_drops_virtual_rows();
  test_blur_full_kernel_sums_to_alpha();
  test_blur_kernel_weights_by_offset();
  test_blur_update_changes_kernel();
  test_blur_clamps_factor();
  test_blur_wraps_column_taps();
  test_blur_pole_row_renormalizes();

  test_chromatic_shift_fanout();

  test_feedback_style_binding();
  test_feedback_plot_is_passthrough();

  test_world_hole_masks_cap();
  test_world_hole_setters();
  test_world_orient_rotates_and_keeps_static_age();
  test_world_orient_motion_blur_sweep_ages();
  test_world_orient_cull_edge_mirrors_plot();
  test_world_orient_slice_selects_by_projection();
  test_world_orient_slice_cull_edge_bounds_all_slices();
  test_world_vertex_replicate_fanout_and_age();
  test_world_vertex_replicate_cull_edge_mirrors_plot();
  test_pipeline_could_intersect_clip_forwards_through_stages();
  test_world_mobius_identity_and_transform();

  test_pipeline_sink_2d_plot_blends_wraps_clips();
  test_pipeline_composition_alpha_and_draw_order();
  test_pipeline_sink_3d_plot_routes_to_canvas();
  test_pipeline_world_replicate_fans_out();
  test_pipeline_2d_into_3d_head_roundtrips();
  test_pipeline_screen_antialias_routes_to_sink();
  test_direct_antialias_sink_framebuffer_parity();
  test_direct_antialias_sink_stale_clip();
  test_feedback_flush_blends_prev_frame();
  test_feedback_north_pole_uses_single_physical_sample();
  test_feedback_south_pole_uses_single_physical_sample();
  test_feedback_polar_rows_use_spherical_footprint();
  test_feedback_flush_respects_clip();
  test_feedback_flush_melt_warp_displaces_south();
  test_feedback_north_cap_uses_exact_control_rows();
  test_feedback_animated_cap_controls_match_compositor_lattice();
  test_feedback_poles_resolve_one_source_longitude();
  test_feedback_half_res_warp_uses_pair_midpoint();
  test_feedback_spherical_ring_control_rows();
  test_feedback_spherical_field_angular_error();
  test_feedback_seam_warp_keeps_its_latitude_row();
  test_feedback_cached_north_cap_clips_share_control_rows();
  test_feedback_warp_cache_matches_uncached();
  test_feedback_polar_rows_hit_their_targets();
  test_feedback_flush_straddled_taps_stay_on_branch();

  test_world_trails_int16_quantization_roundtrip();
  test_world_trails_clamps_out_of_range();
  test_world_trails_capacity_evicts_one_slot();
  test_world_trails_ttl_expiry();
  test_world_trails_set_lifetime_caps_ttl();
  test_world_trails_midbuffer_expiry_reclaims_slot();
  test_screen_trails_store_emit_decay();
  test_screen_trails_negative_age_clamps_t();
  test_screen_trails_forwards_aged_emission();
  test_screen_trails_at_capacity_replaces_last_slot();
  test_mixed_domain_flush_drains_both_buffers();

  test_effect_needs_full_frame_default_false();
  test_screen_trails_banded_matches_full();
  test_screen_trails_alternating_clip_matches_full();
  test_feedback_banded_diverges_from_full();

  return fixture.result();
}

} // namespace filter_tests
} // namespace hs_test
