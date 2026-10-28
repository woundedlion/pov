/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for the SDF shapes, CSG combinators and scanline cull.
 */
#pragma once

#include "core/render/sdf.h"
#include "core/mesh/solids.h"
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
#include <limits>
#include <numbers>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace hs_test {
namespace sdf_tests {

struct RejectMarginShape {
  static constexpr double REJECT_MARGIN = 0.125;
};
static_assert(SDF::reject_margin<RejectMarginShape> == 0.125f);
static_assert(SDF::reject_margin<int> == 0.0f);

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

  HS_RUN(test_fold_phi_in_range);
  HS_RUN(test_fold_phi_negative_reflects);
  HS_RUN(test_fold_phi_above_pi_reflects);
  HS_RUN(test_fold_phi_full_range);
  HS_RUN(test_fold_phi_band_matches_circle_extent);
  HS_RUN(test_fold_phi_band_pole_crossing_poses);
  HS_RUN(test_centered_sector_angle_matches_wrap);

  HS_RUN(test_ring_roundoff_at_exact_center);
  HS_RUN(test_ring_on_centerline);
  HS_RUN(test_ring_inside_band);
  HS_RUN(test_ring_outside_band_returns_sentinel);
  HS_RUN(test_ring_just_outside_band);
  HS_RUN(test_ring_small_radius_distance_symmetric);

  HS_RUN(test_distorted_ring_copies_basis);
  HS_RUN(test_distorted_ring_constant_shift_moves_centerline);
  HS_RUN(test_distorted_ring_sin_shift_varies_by_azimuth);
  HS_RUN(test_distorted_ring_flat_matches_zero_knots);
  HS_RUN(test_distorted_ring_polyline_distance_matches_bruteforce);
  HS_RUN(test_distorted_ring_capped_search_is_lower_bound);
  HS_RUN(test_distorted_ring_closes_without_a_sentinel);
  HS_RUN(test_distorted_ring_knot_extrema_tighten_band);
  HS_RUN(test_distorted_ring_past_reach_reports_far_sentinel);
  HS_RUN(test_distorted_ring_frame_distance_matches_distance);

  HS_RUN(test_polygon_at_center_inside);
  HS_RUN(test_polygon_far_point_outside);

  HS_RUN(test_spherical_polygon_center_inside);
  HS_RUN(test_spherical_polygon_far_outside);
  HS_RUN(test_spherical_polygon_center_and_edge_magnitude);
  HS_RUN(test_spherical_polygon_sine_distance_aa_error);
  HS_RUN(test_spherical_polygon_sine_full_interior);
  HS_RUN(test_spherical_polygon_composes_under_csg);

  HS_RUN(test_star_center_inside);
  HS_RUN(test_star_far_outside);
  HS_RUN(test_star_tip_on_boundary);

  HS_RUN(test_flower_interior_along_petal);
  HS_RUN(test_flower_petal_tip_on_boundary);
  HS_RUN(test_solid_shape_unit_angle_and_no_uv_paths);

  HS_RUN(test_inverted_fill_stays_centered);
  HS_RUN(test_inverted_fill_scans_full_sphere);

  HS_RUN(test_line_on_arc_is_inside);
  HS_RUN(test_line_endpoint_is_on_line);
  HS_RUN(test_line_perpendicular_off);
  HS_RUN(test_line_just_above_cross_threshold);
  HS_RUN(test_line_degenerate_zero_length);
  HS_RUN(test_line_near_coincident_endpoints_stay_point_like);

  HS_RUN(test_torus_on_centerline_is_inside);
  HS_RUN(test_torus_on_surface);
  HS_RUN(test_torus_origin_is_outside_hole);
  HS_RUN(test_torus_normal_points_outward_on_outer_rim);
  HS_RUN(test_torus_normal_points_outward_on_top);

  HS_RUN(test_twist_apply_displaces_y);
  HS_RUN(test_twist_lipschitz_identity_and_closed_form);
  HS_RUN(test_twist_bounding_inflation);
  HS_RUN(test_twisted_torus_matches_recurrence);
  HS_RUN(test_twist_axis_threshold_siblings_agree);
  HS_RUN(test_warped_volume_distance_is_sphere_trace_safe);
  HS_RUN(test_warped_volume_distance_matches_lipschitz_correction);
  HS_RUN(test_twist_correct_normal_unit_length);

  HS_RUN(test_union_picks_closest_shape);

  HS_RUN(test_subtract_inside_a_outside_b_remains_inside);
  HS_RUN(test_subtract_inside_both_becomes_outside);
  HS_RUN(test_subtract_keeps_minuend_size_when_b_wins);
  HS_RUN(test_subtract_solid_b_leaves_the_minuend_uncarved);
  HS_RUN(test_subtract_star_notch_columns_survive_the_carve);
  HS_RUN(test_subtract_empty_b_passes_a_through_verbatim);
  HS_RUN(test_subtract_full_width_b_still_emits_the_minuend);
  HS_RUN(test_subtract_seam_straddle_forwards_minuend);
  HS_RUN(test_subtract_many_arc_preserves_minuend);

  HS_RUN(test_intersection_requires_both_inside);
  HS_RUN(test_intersection_unsorted_child_yields_sorted_result);
  HS_RUN(test_intersection_full_width_child_replays_other);
  HS_RUN(test_intersection_full_scan_emits_no_spans);
  HS_RUN(test_intersection_seam_straddle_overlaps_across_wrap_frames);

  HS_RUN(test_smooth_union_matches_union_far_from_boundary);
  HS_RUN(test_smooth_union_blends_inside_band);
  HS_RUN(test_smooth_union_solidity_follows_children);
  HS_RUN(test_sentinel_clampers_are_not_blendable);
  HS_RUN(test_csg_combinators_reject_temporary_children);

  HS_RUN(test_union_merges_overlapping_intervals);
  HS_RUN(test_union_seam_straddle_merges_overlapping_intervals);
  HS_RUN(test_nested_union_emits_every_child_arc);
  HS_RUN(test_smooth_union_seam_straddle_merges_padded_intervals);
  HS_RUN(test_smooth_union_pad_widens_toward_pole);

  HS_RUN(test_angular_repeat_matches_base_at_zero_angle);
  HS_RUN(test_angular_repeat_creates_copies);
  HS_RUN(test_angular_repeat_t_is_sector_local);

  HS_RUN(test_warped_volume_bounding_distance_never_over_estimates);

  HS_RUN(test_annular_angles_bound_reference);
  HS_RUN(test_cull_covers_interior_over_orientation_grid);
  HS_RUN(test_pole_axis_ring_bounds_skip_pole_rows);
  HS_RUN(test_linearized_ring_bounds_cover_visible_rows);
  HS_RUN(test_intersection_cull_covers_interior_over_polygon_pairs);
  HS_RUN(test_subtract_cull_covers_interior_over_leaf_pairs);
  HS_RUN(test_smooth_union_cull_covers_fringe_over_leaf_pairs);
  HS_RUN(test_smooth_union_scans_rows_past_both_children);
  HS_RUN(test_angular_repeat_non_y_axis_cull_covers_copies);
  HS_RUN(test_angular_repeat_y_axis_cull_keeps_rows_narrows_columns);
  HS_RUN(test_angular_repeat_tilted_axis_forfeits_cull);
  HS_RUN(test_line_arc_bulge_cull_covers_interior);
  HS_RUN(test_face_arc_extrema_bounds);
  HS_RUN(test_line_antipodal_cull_covers_interior);
  HS_RUN(test_line_thick_cap_past_pi_cull_covers_interior);
  HS_RUN(test_ring_pole_wrap_cull_covers_interior);
  HS_RUN(test_distorted_ring_cull_covers_interior_high_freq);
  HS_RUN(test_face_cull_covers_aa_fringe);
  HS_RUN(test_face_hemisphere_clip_covers_fringe);
  HS_RUN(test_face_azimuth_cull_matches_boundary_column);
  HS_RUN(test_face_vertical_margin_tracks_pixel_width);
  HS_RUN(test_face_latitude_pad_reduces_fringe_drops);
  HS_RUN(test_star_polygon_cull_covers_aa_fringe);
  HS_RUN(test_face_pole_vertex_matches_full_scan);
  HS_RUN(test_face_sector_backtrack_sign);
  HS_RUN(test_face_oracle_rejects_invalid_distances);
  HS_RUN(test_face_asymmetric_sector_matches_oracle);
  HS_RUN(test_face_randomized_sector_matches_oracle);
  HS_RUN(test_face_randomized_backtracking_matches_oracle);
  HS_RUN(test_face_tied_rows_match_oracle);
  HS_RUN(test_face_shipping_recipes_match_oracle);
  HS_RUN(test_face_handle_reuses_arena_storage);
  HS_RUN(test_face_distance_matches_exact_oracle);
  HS_RUN(test_face_convex_rejection_cutoff);
  HS_RUN(test_face_plane_box_rejection);
  HS_RUN(test_face_class_lut_matches_oracle);

  return fixture.result();
}

} // namespace sdf_tests
} // namespace hs_test
