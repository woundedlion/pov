/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/render/plot.h and the ClipRegion helpers in
 * core/render/clip.h.
 */
#pragma once
#include "tests/pixel_test_util.h"

#include "core/animation/orientation.h"
#include "core/render/plot.h"
#include "core/render/scan.h"
#include "core/animation/animation.h" // Segue::Dissolve
#include "core/render/filter.h"
#include "core/math/geometry.h"
#include "core/render/canvas.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/mesh_test_util.h"
#include "tests/vec_test_util.h"

#include <algorithm>
#include <bit>
#include <limits>
#include <vector>

namespace hs_test {
namespace plot_scan_tests {

template <typename Chart>
concept BorrowableRasterBasis = requires(Chart &&basis) {
  Plot::RasterProjection::planar(std::forward<Chart>(basis));
};

static_assert(BorrowableRasterBasis<math::Basis &>);
static_assert(BorrowableRasterBasis<const math::Basis &>);
static_assert(!BorrowableRasterBasis<math::Basis>);
static_assert(!BorrowableRasterBasis<const math::Basis>);

static_assert(std::is_aggregate_v<Plot::RasterOptions>);
static_assert(std::is_trivially_copyable_v<Plot::RasterOptions>);
constexpr size_t RASTER_OPTIONS_BYTE_BUDGET = sizeof(void *) == 8 ? 88 : 48;
static_assert(sizeof(Plot::RasterOptions) <= RASTER_OPTIONS_BYTE_BUDGET);
static_assert(
    !std::is_constructible_v<Plot::RasterLoop, bool, const Fragment *>);
static_assert(
    !std::is_constructible_v<Plot::RasterProjection, const math::Basis *,
                             std::span<const uint8_t>>);
static_assert(std::is_constructible_v<Plot::PointProjections,
                                      const float (&)[3], const float (&)[3]>);
static_assert(!std::is_constructible_v<Plot::PointProjections,
                                       const float (&)[3], const float (&)[2]>);
static_assert(!std::is_constructible_v<Plot::PointProjections, const float *,
                                       const float *, size_t>);
static_assert(!Plot::RasterLoop{}.is_closed());
static_assert(Plot::RasterLoop::closed().is_closed());

#include "tests/plot_scan/sampling_fixture.h"
#include "tests/plot_scan/geodesic_trig.h"
#include "tests/plot_scan/raster_fixture.h"
#include "tests/plot_scan/line.h"
#include "tests/plot_scan/row_clip.h"
#include "tests/plot_scan/row_spans.h"
#include "tests/plot_scan/column_spans.h"
#include "tests/plot_scan/mesh_edges.h"
#include "tests/plot_scan/screen_step.h"
#include "tests/plot_scan/ring.h"
#include "tests/plot_scan/distorted_ring.h"
#include "tests/plot_scan/multiline_sample.h"
#include "tests/plot_scan/multiline_draw.h"
#include "tests/plot_scan/star_flower.h"
#include "tests/plot_scan/rasterize.h"
#include "tests/plot_scan/particles.h"
#include "tests/plot_scan/planar_metric.h"
#include "tests/plot_scan/single_pass.h"
#include "tests/plot_scan/planar_chords.h"
/**
 * @brief Runs every plot/scan sampling test in this module.
 * @return Number of failed assertions reported across the module's tests.
 */
inline int run_plot_scan_tests() {
  hs_test::ModuleFixture fixture("plot_scan");

  test_four_regular_and_medial_edge_extraction();

  test_geodesic_sincos_bit_parity();
  test_line_sample_endpoints_and_unit_length();
  test_line_sample_interior_between_endpoints();
  test_line_sample_degenerate_segment();
  test_line_sample_antipodal_stable_axis();
  test_line_sample_near_antipodal_ulp_stable_axis();

  test_clip_could_intersect_y();
  test_clip_x_band_topologies();
  test_clip_x_wrap_matches_modulo();
  test_row_span_covers_arc_bulge();
  test_cap_may_touch_clip_is_conservative();
  test_clip_arcs_overlap();
  test_col_span_rejects_ill_conditioned_pole();
  test_col_span_covers_arc();
  test_edge_visible_in_clip_is_conservative();
  test_rasterize_column_cull_pixel_parity();
  test_mesh_edge_gate_pixel_parity();
  test_rasterize_window_preserves_terminal_sample();
  test_short_geodesic_arc_lengths();
  test_rasterize_short_edge_windows();
  test_mesh_dissolve_masks_partition_edges();
  test_gate_trail_column_cull_honors_unbounded_edge();
  test_raw_geodesic_edge_gate_parity();
  test_finish_col_span_one_period();
  test_wrap_one_period_matches_modulo();
  test_cartesian_quadrant_gate_classification();
  test_cartesian_quadrant_gate_is_conservative();
  test_gate_trail_edges_matches_edge_visible();
  test_mesh_clip_cut_separates_band();
  test_rasterize_gate_bits_pixel_parity();
  test_screen_step_matches_analytic_unclamped();
  test_edge_fits_one_dot_is_conservative();
  test_antialiased_dot_clip_footprint();
  test_geodesic_edge_gate_keeps_upper_antialias_tap();
  test_antialiased_dot_gate_matches_antialias_taps();

  test_ring_sample_unit_length_and_progress();
  test_ring_sample_lut_matches_direct();
  test_ring_draw_stride_tracks_full_grid();
  test_ring_draw_accepts_direct_sink();
  test_planar_chords_match_rasterize_brightness();
  test_planar_chords_pole_split_matches_whole_star();
  test_planar_band_split_matches_whole_polyline();

  test_distorted_ring_sample_angle_addition_identity();
  test_distorted_ring_shift_matches_fn_point();

  test_multiline_sample_arclength_param();

  test_star_sample_unit_length_closed();
  test_star_sample_radius_trig_parity();
  test_star_continuous_matches_standard_near_side();
  test_star_continuous_crosses_equator();
  test_star_continuous_collapses_at_antipode();
  test_flower_sample_unit_length_closed();

  // rasterize() control-flow coverage.
  test_rasterize_subpixel_open_segment_plots_both_endpoints();
  test_rasterize_open_segment_gap_free();
  test_rasterize_closed_loop_gap_free_no_dup();
  test_rasterize_antipodal_seam_planar_falls_back_geodesic();
  test_rasterize_planar_segment_gap_free_arclength();
  test_rasterize_planar_arc_registers_track_drawn_arc();
  test_planar_sampler_from_cull_parity();
  test_rasterize_planar_policy_parity();
  test_rasterize_sampling_follows_world_transforms();
  test_screen_step_axes_match_stage_walk();
  test_rasterize_cull_follows_filter_orientation();

  test_particle_system_draws_active_trails_with_registers();
  test_particle_system_empty_zero_lifetime_is_noop();
  test_particle_system_skips_unrenderable_trails();
  test_particle_system_direct_trail_materialization_registers();
  test_particle_system_sparse_history_live_tip();
  test_particle_system_v0_zero_at_oldest_sample();
  test_particle_system_custom_v2_mapper();
  test_particle_system_direct_trail_materialization_output_parity();
  test_particle_system_deferred_shader_parity_and_skip();
  test_particle_system_fused_optional_deferred_shader();
  test_particle_system_gate_pixel_parity_random_trails();
  test_particle_system_subpixel_trail_dot_parity();

  test_azimuthal_project_radius_is_geodesic_angle();
  test_azimuthal_roundtrip_identity();
  test_azimuthal_unproject_hits_great_circle_point();
  test_planar_arc_length_matches_fine_quadrature();
  test_dual_metric_radial_vs_azimuthal();
  test_planar_arc_cumul_monotone_and_endpoints();

  test_multiline_draw_covers_only_its_geodesic_edges();
  test_arc_angular_distance_clamps_to_minor_arc();
  test_multiline_draw_closed_adds_the_seam_edge();
  test_plot_line_over_pole_reaches_row0();
  test_plot_line_antipodal_replay_parameter();

  test_planar_one_pass_matches_forward_difference();
  test_planar_one_pass_tangent_is_forward_and_orthogonal();
  test_rasterize_single_pass_planar_matches_two_pass();
  test_rasterize_single_pass_closed_loop_matches_two_pass();
  test_rasterize_single_pass_balances_terminal_interval();
  test_rasterize_step_budget_backstop_finishes_segment();
  test_rasterize_default_sampling_policy_parity();
  test_rasterize_balanced_sampling_scope();
  test_rasterize_balanced_sampling_density_and_alpha();
  test_rasterize_balanced_pole_guard();
  test_rasterize_balanced_geodesic_density_and_alpha();
  test_rasterize_balanced_high_alpha_saturates();
  test_rasterize_balanced_star_visual_budget();
  test_rasterize_single_pass_geodesic_endpoints_and_omit_end();
  test_rasterize_single_pass_geodesic_stress_arcs_are_gap_free();
  test_rasterize_single_pass_geodesic_quadrant_clip_parity();

  return fixture.result();
}

} // namespace plot_scan_tests
} // namespace hs_test
