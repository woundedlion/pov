/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for the SDF shapes, CSG combinators and scanline cull.
 */
#pragma once

#include "core/render/sdf.h"
#include "core/render/sdf/volume.h"
#include "core/render/sdf/face_class_bake.h"
#include "core/render/scan.h"
#include "core/math/geometry.h"
#include "tests/vec_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cfloat>
#include <cmath>
#include <cstring>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

namespace hs_test {
namespace sdf_tests {

/**
 * @brief Builds the canonical equator-facing basis: v = +Y, u = +X, w = +Z.
 * @return A Basis oriented so its pole points along +Y.
 */
inline math::Basis equator_basis() {
  return math::Basis{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                     math::Vector(0, 0, 1)};
}

#include "tests/sdf/primitives.h"
#include "tests/sdf/warped_volume.h"
#include "tests/sdf/csg_combinators.h"
#include "tests/sdf/interval_cull.h"
#include "tests/sdf/face_oracles.h"

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every sdf test case.
 * @return The module's failure count.
 */
inline int run_sdf_tests() {
  hs_test::ModuleFixture fixture("sdf");

  test_clamp_phi_in_range();
  test_clamp_phi_negative_reflects();
  test_clamp_phi_above_pi_reflects();
  test_clamp_phi_full_range();
  test_clamp_phi_band_matches_circle_extent();
  test_clamp_phi_band_pole_crossing_poses();
  test_centered_sector_angle_matches_wrap();

  test_ring_roundoff_at_exact_center();
  test_ring_on_centerline();
  test_ring_inside_band();
  test_ring_outside_band_returns_sentinel();
  test_ring_just_outside_band();
  test_ring_small_radius_distance_symmetric();

  test_distorted_ring_constant_shift_moves_centerline();
  test_distorted_ring_sin_shift_varies_by_azimuth();
  test_distorted_ring_flat_matches_zero_knots();
  test_distorted_ring_polyline_distance_matches_bruteforce();
  test_distorted_ring_closes_without_a_sentinel();
  test_distorted_ring_knot_extrema_tighten_band();
  test_distorted_ring_past_reach_reports_far_sentinel();
  test_distorted_ring_frame_distance_matches_distance();

  test_polygon_at_center_inside();
  test_polygon_far_point_outside();

  test_spherical_polygon_center_inside();
  test_spherical_polygon_far_outside();
  test_spherical_polygon_center_and_edge_magnitude();
  test_spherical_polygon_sine_distance_aa_error();
  test_spherical_polygon_sine_full_interior();
  test_spherical_polygon_composes_under_csg();

  test_star_center_inside();
  test_star_far_outside();
  test_star_tip_on_boundary();

  test_flower_interior_along_petal();
  test_flower_petal_tip_on_boundary();
  test_solid_shape_unit_angle_and_no_uv_paths();

  test_inverted_fill_stays_centered();
  test_inverted_fill_scans_full_sphere();

  test_line_on_arc_is_inside();
  test_line_endpoint_is_on_line();
  test_line_perpendicular_off();
  test_line_just_above_cross_threshold();
  test_line_degenerate_zero_length();
  test_line_near_coincident_endpoints_stay_point_like();

  test_torus_on_centerline_is_inside();
  test_torus_on_surface();
  test_torus_origin_is_outside_hole();
  test_torus_normal_points_outward_on_outer_rim();
  test_torus_normal_points_outward_on_top();

  test_twist_apply_displaces_y();
  test_twist_lipschitz_identity_and_closed_form();
  test_twist_bounding_inflation();
  test_twisted_torus_matches_recurrence();
  test_twist_axis_threshold_siblings_agree();
  test_warped_volume_distance_is_sphere_trace_safe();
  test_warped_volume_distance_matches_lipschitz_correction();
  test_twist_correct_normal_unit_length();

  test_union_picks_closest_shape();

  test_subtract_inside_a_outside_b_remains_inside();
  test_subtract_inside_both_becomes_outside();
  test_subtract_keeps_minuend_size_when_b_wins();
  test_subtract_solid_b_leaves_the_minuend_uncarved();
  test_subtract_star_notch_columns_survive_the_carve();
  test_subtract_empty_b_passes_a_through_verbatim();
  test_subtract_full_width_b_still_emits_the_minuend();
  test_subtract_seam_straddle_forwards_minuend();
  test_subtract_many_arc_preserves_minuend();

  test_intersection_requires_both_inside();
  test_intersection_unsorted_child_yields_sorted_result();
  test_intersection_full_width_child_replays_other();
  test_intersection_full_scan_emits_no_spans();
  test_intersection_seam_straddle_overlaps_across_wrap_frames();

  test_smooth_union_matches_union_far_from_boundary();
  test_smooth_union_blends_inside_band();
  test_smooth_union_solidity_follows_children();
  test_sentinel_clampers_are_not_blendable();
  test_csg_combinators_reject_temporary_children();

  test_union_merges_overlapping_intervals();
  test_union_seam_straddle_merges_overlapping_intervals();
  test_nested_union_emits_every_child_arc();
  test_smooth_union_seam_straddle_merges_padded_intervals();
  test_smooth_union_pad_widens_toward_pole();

  test_angular_repeat_matches_base_at_zero_angle();
  test_angular_repeat_creates_copies();
  test_angular_repeat_t_is_sector_local();

  test_warped_volume_bounding_distance_never_over_estimates();

  test_annular_angles_bound_reference();
  test_cull_covers_interior_over_orientation_grid();
  test_pole_axis_ring_bounds_skip_pole_rows();
  test_linearized_ring_bounds_cover_visible_rows();
  test_intersection_cull_covers_interior_over_polygon_pairs();
  test_subtract_cull_covers_interior_over_leaf_pairs();
  test_smooth_union_cull_covers_fringe_over_leaf_pairs();
  test_smooth_union_scans_rows_past_both_children();
  test_angular_repeat_non_y_axis_cull_covers_copies();
  test_angular_repeat_y_axis_cull_narrows_rows();
  test_angular_repeat_tilted_axis_forfeits_cull();
  test_line_arc_bulge_cull_covers_interior();
  test_line_antipodal_cull_covers_interior();
  test_line_thick_cap_past_pi_cull_covers_interior();
  test_ring_pole_wrap_cull_covers_interior();
  test_distorted_ring_cull_covers_interior_high_freq();
  test_face_cull_covers_aa_fringe();
  test_face_azimuth_cull_matches_boundary_column();
  test_face_vertical_margin_tracks_pixel_width();
  test_face_latitude_pad_reduces_fringe_drops();
  test_star_polygon_cull_covers_aa_fringe();
  test_face_pole_vertex_matches_full_scan();
  test_face_sector_backtrack_sign();
  test_face_distance_matches_exact_oracle();
  test_face_class_lut_matches_oracle();

  return fixture.result();
}

} // namespace sdf_tests
} // namespace hs_test
