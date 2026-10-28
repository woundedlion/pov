/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/render/plot.h and the ClipRegion helpers in
 * core/render/clip.h.
 */
#pragma once
#include "tests/pixel_test_util.h"
#include "tests/temporal_oracle.h"

#include "core/animation/orientation.h"
#include "core/render/plot.h"
#include "core/render/scan.h"
#include "core/animation/animation.h" // Segue::Dissolve
#include "core/render/filter.h"
#include "core/math/geometry.h"
#include "core/render/canvas.h"
#include "core/render/spacetime/crossing.h"
#include "core/render/spacetime/pose_table.h"
#include "core/render/spacetime/swept_bound.h"
#include "core/render/spacetime/worldlines.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/mesh_test_util.h"
#include "tests/vec_test_util.h"

#include <algorithm>
#include <bit>
#include <cstring>
#include <functional>
#include <memory>
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
static_assert(sizeof(void *) != 8 || sizeof(Plot::RasterOptions) <= 88);
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
#include "tests/plot_scan/temporal_seam.h"
#include "tests/plot_scan/timed_mesh.h"
#include "tests/plot_scan/temporal_poses.h"
#include "tests/plot_scan/temporal_motion.h"
#include "tests/plot_scan/temporal_crossings.h"
#include "tests/plot_scan/temporal_swept.h"
#include "tests/plot_scan/point_route.h"
/**
 * @brief Runs every plot/scan sampling test in this module.
 * @return Number of failed assertions reported across the module's tests.
 */
inline int run_plot_scan_tests() {
  hs_test::ModuleFixture fixture("plot_scan");

  HS_RUN(hs_test::test_reference_row_displacement_is_analytic);
  HS_RUN(hs_test::test_temporal_oracle_uses_whole_frame_endpoint_budget);
  HS_RUN(test_retained_endpoint_lifecycle);
  HS_RUN(test_timed_projected_entry_bypasses_world);
  HS_RUN(test_timed_mesh_endpoint_cache);
  HS_RUN(test_column_edge_quality);
  HS_RUN(test_column_polar_arc_coverage);
  HS_RUN(test_column_crossing_parameter_bound);
  HS_RUN(test_quadratic_screen_capture_meets_surface);
  HS_RUN(test_column_basis_crossing_orientation);
  HS_RUN(test_column_walk_arm_boundary_coverage);
  HS_RUN(test_column_walk_terminal_on_shrinking_edge);
  HS_RUN(test_column_edge_cull_is_conservative);
  HS_RUN(test_column_edge_cull_keeps_padded_point_sets);
  HS_RUN(test_column_physical_north_cap);
  HS_RUN(test_column_directional_spacing_certificate);
  HS_RUN(test_column_steady_lattice);
  HS_RUN(test_column_source_basis_transport);
  HS_RUN(test_column_uniform_color_preserves_age);
  HS_RUN(test_column_uniform_color_seam_registers);
  HS_RUN(test_column_capture_reads_vertex_paths);
  HS_RUN(test_column_plan_covers_lit_columns);
  HS_RUN(test_column_mesh_batch_cull);
  HS_RUN(test_column_single_splat_motion_gate);
  HS_RUN(test_column_basis_normalization);
  HS_RUN(test_stationary_column_mesh_skips_geometry_cache);
  HS_RUN(test_column_walk_chained_ring_plots_vertices_once);
  HS_RUN(test_column_walk_seam_taps_shade_at_their_times);
  HS_RUN(test_trail_gate_uses_motion_region);
  HS_RUN(test_four_regular_and_medial_edge_extraction);
  HS_RUN(test_point_route_polar_sweep_support);
  HS_RUN(test_geometry_callback_aggregate_contract);
  HS_RUN(test_point_route_piece_tap_policy);
  HS_RUN(test_point_route_physical_speed_covers_equator_cell);
  HS_RUN(test_point_route_small_chord_certificate);
  HS_RUN(test_point_route_physical_curvature_bound);
  HS_RUN(test_point_route_curved_cap_enters_clip);
  HS_RUN(test_point_route_matches_column_oracle);
  HS_RUN(test_point_route_gate_keeps_visible_taps);
  HS_RUN(test_point_route_shades_at_pass_time);
  HS_RUN(test_steady_edge_plots_lattice_points);

  HS_RUN(test_pose_table_follows_history_path);
  HS_RUN(test_pose_table_holds_a_still_history);
  HS_RUN(test_pose_table_rotates_between_column_times);
  HS_RUN(test_column_poses_match_endpoint_slerp);
  HS_RUN(test_interpolants_bound_rotation);
  HS_RUN(test_interpolants_bound_meridian);
  HS_RUN(test_interpolants_bound_orbit);
  HS_RUN(test_interpolants_bound_curved_path);
  HS_RUN(test_interpolants_bound_loxodrome);
  HS_RUN(test_interpolants_bound_warped_great_circle);
  HS_RUN(test_vertex_motion_carries_frame_ends);
  HS_RUN(test_rotation_chord_bound_covers_unit_vectors);
  HS_RUN(test_arm_passes_match_brute_force);
  HS_RUN(test_arm_passes_extend_to_frame_ends);
  HS_RUN(test_piecewise_arm_passes_match_brute_force);
  HS_RUN(test_model_track_guard_bounds_true_passes);
  HS_RUN(test_model_track_splits_or_refuses_near_poles);
  HS_RUN(test_swept_bounds_cover_moving_edges);
  HS_RUN(test_swept_heights_refuse_wide_spreads);

  HS_RUN(test_geodesic_sincos_bit_parity);
  HS_RUN(test_line_sample_endpoints_and_unit_length);
  HS_RUN(test_line_sample_interior_between_endpoints);
  HS_RUN(test_line_sample_degenerate_segment);
  HS_RUN(test_line_sample_antipodal_stable_axis);
  HS_RUN(test_line_sample_near_antipodal_ulp_stable_axis);

  HS_RUN(test_clip_could_intersect_y);
  HS_RUN(test_clip_x_band_topologies);
  HS_RUN(test_clip_x_wrap_matches_modulo);
  HS_RUN(test_row_span_covers_arc_bulge);
  HS_RUN(test_cap_may_touch_clip_is_conservative);
  HS_RUN(test_visible_chunk_mask_covers_sampled_ring);
  HS_RUN(test_for_each_chunk_span_tiles_columns);
  HS_RUN(test_clip_arcs_overlap);
  HS_RUN(test_col_span_rejects_ill_conditioned_pole);
  HS_RUN(test_col_span_covers_arc);
  HS_RUN(test_edge_visible_in_clip_is_conservative);
  HS_RUN(test_rasterize_column_cull_pixel_parity);
  HS_RUN(test_mesh_edge_gate_pixel_parity);
  HS_RUN(test_rasterize_window_preserves_terminal_sample);
  HS_RUN(test_short_geodesic_arc_lengths);
  HS_RUN(test_rasterize_short_edge_windows);
  HS_RUN(test_mesh_dissolve_masks_partition_edges);
  HS_RUN(test_gate_trail_column_cull_honors_unbounded_edge);
  HS_RUN(test_raw_geodesic_edge_gate_parity);
  HS_RUN(test_finish_col_span_one_period);
  HS_RUN(test_wrap_one_period_matches_modulo);
  HS_RUN(test_cartesian_quadrant_gate_classification);
  HS_RUN(test_cartesian_latitude_gate_covers_row_map_error);
  HS_RUN(test_cartesian_quadrant_gate_is_conservative);
  HS_RUN(test_gate_trail_edges_matches_edge_visible);
  HS_RUN(test_mesh_clip_cut_separates_band);
  HS_RUN(test_rasterize_gate_bits_pixel_parity);
  HS_RUN(test_screen_step_matches_analytic_unclamped);
  HS_RUN(test_edge_fits_one_dot_is_conservative);
  HS_RUN(test_antialiased_dot_clip_footprint);
  HS_RUN(test_geodesic_edge_gate_keeps_upper_antialias_tap);
  HS_RUN(test_antialiased_dot_gate_matches_antialias_taps);

  HS_RUN(test_ring_sample_unit_length_and_progress);
  HS_RUN(test_ring_sample_lut_matches_direct);
  HS_RUN(test_ring_draw_stride_tracks_full_grid);
  HS_RUN(test_ring_draw_accepts_direct_sink);
  HS_RUN(test_planar_chords_match_rasterize_brightness);
  HS_RUN(test_planar_chords_pole_split_matches_whole_star);
  HS_RUN(test_planar_band_split_matches_whole_polyline);

  HS_RUN(test_distorted_ring_sample_angle_addition_identity);
  HS_RUN(test_distorted_ring_shift_matches_fn_point);

  HS_RUN(test_multiline_sample_arclength_param);

  HS_RUN(test_star_sample_unit_length_closed);
  HS_RUN(test_star_sample_radius_trig_parity);
  HS_RUN(test_star_continuous_matches_standard_near_side);
  HS_RUN(test_star_continuous_crosses_equator);
  HS_RUN(test_star_continuous_collapses_at_antipode);
  HS_RUN(test_flower_sample_unit_length_closed);

  // rasterize() control-flow coverage.
  HS_RUN(test_rasterize_subpixel_open_segment_plots_both_endpoints);
  HS_RUN(test_rasterize_open_segment_gap_free);
  HS_RUN(test_rasterize_closed_loop_gap_free_no_dup);
  HS_RUN(test_rasterize_antipodal_seam_planar_falls_back_geodesic);
  HS_RUN(test_rasterize_planar_segment_gap_free_arclength);
  HS_RUN(test_rasterize_planar_arc_registers_track_drawn_arc);
  HS_RUN(test_planar_sampler_from_cull_parity);
  HS_RUN(test_rasterize_planar_policy_parity);
  HS_RUN(test_rasterize_register_copy_follows_interpolate_registers);
  HS_RUN(test_rasterize_planar_arc_registers_skip_culled_segments);
  HS_RUN(test_rasterize_sampling_follows_world_transforms);
  HS_RUN(test_screen_step_axes_match_stage_walk);
  HS_RUN(test_rasterize_cull_follows_filter_orientation);

  HS_RUN(test_particle_system_draws_active_trails_with_registers);
  HS_RUN(test_particle_system_empty_zero_lifetime_is_noop);
  HS_RUN(test_particle_system_skips_unrenderable_trails);
  HS_RUN(test_particle_system_direct_trail_materialization_registers);
  HS_RUN(test_particle_system_sparse_history_live_tip);
  HS_RUN(test_particle_system_v0_zero_at_oldest_sample);
  HS_RUN(test_particle_system_custom_v2_mapper);
  HS_RUN(test_particle_system_direct_trail_materialization_output_parity);
  HS_RUN(test_particle_system_deferred_shader_parity_and_skip);
  HS_RUN(test_particle_system_fused_optional_deferred_shader);
  HS_RUN(test_particle_system_gate_pixel_parity_random_trails);
  HS_RUN(test_particle_system_subpixel_trail_dot_parity);

  HS_RUN(test_azimuthal_project_radius_is_geodesic_angle);
  HS_RUN(test_azimuthal_roundtrip_identity);
  HS_RUN(test_azimuthal_unproject_hits_great_circle_point);
  HS_RUN(test_planar_arc_length_matches_fine_quadrature);
  HS_RUN(test_dual_metric_radial_vs_azimuthal);
  HS_RUN(test_planar_arc_cumul_monotone_and_endpoints);

  HS_RUN(test_multiline_draw_covers_only_its_geodesic_edges);
  HS_RUN(test_arc_angular_distance_clamps_to_minor_arc);
  HS_RUN(test_multiline_draw_closed_adds_the_seam_edge);
  HS_RUN(test_plot_line_over_pole_reaches_row0);
  HS_RUN(test_plot_line_antipodal_replay_parameter);

  HS_RUN(test_planar_one_pass_matches_forward_difference);
  HS_RUN(test_planar_one_pass_tangent_is_forward_and_orthogonal);
  HS_RUN(test_rasterize_single_pass_planar_matches_two_pass);
  HS_RUN(test_rasterize_single_pass_closed_loop_matches_two_pass);
  HS_RUN(test_rasterize_single_pass_balances_terminal_interval);
  HS_RUN(test_rasterize_step_budget_backstop_finishes_segment);
  HS_RUN(test_rasterize_default_sampling_policy_parity);
  HS_RUN(test_rasterize_balanced_sampling_scope);
  HS_RUN(test_polar_stroke_alpha_is_bounded);
  HS_RUN(test_rasterize_balanced_sampling_density_and_alpha);
  HS_RUN(test_rasterize_balanced_pole_guard);
  HS_RUN(test_rasterize_balanced_geodesic_density_and_alpha);
  HS_RUN(test_rasterize_balanced_high_alpha_saturates);
  HS_RUN(test_rasterize_balanced_star_visual_budget);
  HS_RUN(test_rasterize_single_pass_geodesic_endpoints_and_omit_end);
  HS_RUN(test_rasterize_single_pass_geodesic_stress_arcs_are_gap_free);
  HS_RUN(test_rasterize_single_pass_geodesic_quadrant_clip_parity);

  return fixture.result();
}

} // namespace plot_scan_tests
} // namespace hs_test
