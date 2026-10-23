/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/render/scan.h using framebuffer readback, recording
 * callbacks, and pixel sinks.
 */
#pragma once

#include "core/render/scan.h"
#include "core/render/sdf/volume.h"
#include "core/render/plot.h"
#include "core/render/filter.h"
#include "core/render/canvas.h"
#include "core/render/spacetime.h"
#include "core/math/geometry.h"
#include "tests/mesh_test_util.h"
#include "tests/pixel_test_util.h"
#include "tests/volume_reference.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cfloat>
#include <cstring>
#include <memory>
#include <vector>

namespace hs_test {
namespace scan_tests {

template <int N, int LUT_N>
inline void build_ring_knots(float (&knots)[N][LUT_N + 1]) {
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < LUT_N; ++k) {
      const float t = 2.0f * math::PI_F * k / LUT_N;
      knots[i][k] = 0.06f * sinf((i + 2) * t) + 0.03f * cosf(3.0f * t + i);
    }
    knots[i][LUT_N] = knots[i][0];
  }
}

/**
 * @brief Pins Render::pole_lod_aggressiveness for a scope and restores it on exit.
 */
struct ScopedPoleLod {
  float saved; /**< Value in force at construction. */
  /**
   * @brief Pins the knob.
   * @param value Aggressiveness held for the scope.
   */
  explicit ScopedPoleLod(float value) : saved(Render::pole_lod_aggressiveness) {
    Render::pole_lod_aggressiveness = value;
  }
  ~ScopedPoleLod() { Render::pole_lod_aggressiveness = saved; }
  ScopedPoleLod(const ScopedPoleLod &) = delete;
  ScopedPoleLod &operator=(const ScopedPoleLod &) = delete;
};

#include "tests/scan/shader_draw.h"
#include "tests/scan/ring_raster.h"
#include "tests/scan/region_pole_lod.h"
#include "tests/scan/stroke_contract.h"
#include "tests/scan/filled_shapes.h"
#include "tests/scan/volume_march.h"
#include "tests/scan/circle_point.h"
#include "tests/scan/column_basis.h"

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs the scan tests under the "scan" module scope.
 * @return The number of test failures recorded by the module.
 */
inline int run_scan_tests() {
  hs_test::ModuleFixture fixture("scan");

  test_replicated_clip_matches_full_frame();
  test_bounding_sphere_initializes_trig();
  test_min_alpha_boundary();
  test_scan_epilogue_contract();
  test_shader_constant_fills_canvas();
  test_shader_ssaa_premultiplies_partial_coverage();
  test_shader_split_ssaa_averages_subsamples();
  test_shader_positional_maps_latitude();
  test_shader_respects_clip_band();
  test_ssaa_grid_sample_positions();
  test_shader_clip_arc_matches_predicate();
  test_shader_render_x_window_matches_predicate();
  test_simulate_in_slabs_orders_substeps_and_columns();
  test_ring_rasterize_produces_bounded_output();
  test_ring_long_radius_azimuth_unflipped();
  test_ring_rasterize_lights_expected_row();
  test_stroke_aa_is_monotone_ramp();
  test_ring_rasterize_empty_clip_draws_nothing();
  test_distorted_ring_flat_matches_zero_knot_raster();
  test_ring_group_matches_sequential();
  test_distorted_ring_candidates_outside_poles();
  test_distorted_ring_stack_matches_sequential();
  test_distorted_ring_stack_empty_clip_skips_table();
  test_fused_walks_ignore_pole_lod();
  test_face_rasterize_matches_scan_region();
  test_scan_shader_v2_contract();
  test_scan_region_seam_no_double_plot();
  test_scan_region_fractional_boundary_no_double_plot();
  test_scan_region_clip_arc_matches_predicate();
  test_pole_lod_runs_are_canvas_anchored();
  test_pole_lod_run_clamps_to_max_run();
  test_pole_lod_shading_matches_undecimated();
  test_pole_lod_concave_face_matches_undecimated();
  test_report_stretch_forwards_through_csg();
  test_prepared_face_reports_match_unprepared();
  test_convex_face_pole_fringe_matches_undecimated();
  test_csg_stroke_aa_uses_winning_child_thickness();

  test_star_pixel_placement();
  test_planar_polygon_pixel_placement();
  test_flower_pixel_placement();
  test_solid_color_path_matches_generic();
  test_spherical_sine_distance_framebuffer_error();
  test_overlapping_fills_composite_blend();

  test_point_draws_the_analytic_cap();
  test_pole_centred_cap_takes_the_full_row_scan();
  test_circle_and_point_keep_exact_pixel_centers();
  test_circle_and_point_match_their_rings();
  test_circle_extent_follows_its_radius();

  test_column_basis_matches_rotated_shapes();
  test_column_basis_follows_history_path();
  test_column_basis_segments_and_endpoints();
  test_column_basis_blocks_cover_their_rotation();
  test_point_birth_skips_earlier_columns();
  test_lit_caps_keep_swept_frames();
  test_lit_caps_follow_history_rate();
  test_swept_proxies_keep_lit_pixels();
  test_rotate_by_matches_exact_rotation();
  test_swept_cap_clip_keeps_lit_clips();
  test_volume_columns_match_rotated_placements();
  test_ring_group_previous_axes_match_column_rings();
  test_normalized_chord_error_bounds_the_turn();
  test_ring_group_column_rows_cover_band();
  test_ring_group_column_rows_window_matches_full();
  test_object_column_basis_matches_rotated_meridian();
  test_ring_stack_honours_render_window();
  test_cap_may_touch_window_covers_window_pixels();
  test_timed_shader_reads_column_sweep_time();
  test_face_palette_shader_moves_at_column_time();
  test_volume_timed_shader_reads_column_sweep_time();
  test_sweep_column_inverts_sweep_fraction();
  test_column_basis_prepares_for_moving_state();
  test_drift_matches_shape_at_column_centre();
  test_settled_drift_matches_column_drift();
  test_mesh_capture_unwrap();
  test_capture_window_runs();
  test_capture_window_faces();
  test_mesh_capture_vertex();
  test_vertex_arc_and_reach_flags();
  test_mesh_capture_reuse();
  test_mesh_seam_face_copies();
  test_mesh_wide_face_copies();
  test_mesh_deformation_shared_edges();
  test_mesh_polar_face_sweep();
  test_mesh_polar_face_windows();
  test_mesh_motion_activation();
  test_recorded_face_caps_match_measured_columns();
  test_moving_faces_outside_clip_rows_hold();

  test_transformed_volume_world_local_roundtrip();
  test_volume_raymarch_silhouette_and_registers();
  test_volume_draw_occluded_edge_blends_over_background();
  test_volume_scalar_state_differential();
  test_volume_dense_reference();
  test_volume_trace_nearly_tied_minimum();
  test_volume_trace_closest_stops_at_first_graze();
  test_volume_probe_occluder_reports_background_graze_point();
  test_volume_trace_closest_overrelax_never_skips_surface();

  return fixture.result();
}

} // namespace scan_tests
} // namespace hs_test
