/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/color/color.h and the palette layer built on it.
 */
#pragma once

#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>

#include "core/color/color.h"
#include "core/color/composition.h"
#include "effects/common/palette_recipes.h"
#include "core/color/generative_palette.h"
#include "core/color/modulation.h"
#include "core/pattern/noise.h"
#include "core/pattern/sphere_lut.h"
#include "core/color/palette_cycler.h"
#include "core/color/srgb_decode.h"
#include "core/memory.h"
#include "tests/color_test_util.h"
#include "tests/pixel_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace color_tests {

#include "tests/color/pixel_blend.h"
#include "tests/color/oklab_gamut.h"
#include "tests/color/srgb_transfer.h"
#include "tests/color/gradient_bake.h"
#include "tests/color/generative_palette_cases.h"
#include "tests/color/palette_cycler_cases.h"
#include "tests/color/palette_layer.h"
#include "tests/color/modulation.h"
#include "tests/pattern/sphere_lut.h"

inline void test_wrap_angle_pi_large_arguments() {
  for (float angle : {25700.0f, -25700.0f, 1000.0f * math::TWO_PI_F + 0.25f}) {
    const float wrapped = wrap_angle_pi(angle);
    float expected = fmodf(angle, math::TWO_PI_F);
    if (expected > math::PI_F)
      expected -= math::TWO_PI_F;
    if (expected < -math::PI_F)
      expected += math::TWO_PI_F;
    HS_EXPECT_NEAR(wrapped, expected, 1e-6f);
    HS_EXPECT_TRUE(wrapped >= -math::PI_F && wrapped <= math::PI_F);
  }
}

/**
 * @brief Verifies an exact half turn keeps its sign and the first value past
 *        it lands on the opposite side of the seam.
 */
inline void test_wrap_angle_pi_half_turn_keeps_sign() {
  HS_EXPECT_EQ(wrap_angle_pi(math::PI_F), math::PI_F);
  HS_EXPECT_EQ(wrap_angle_pi(-math::PI_F), -math::PI_F);
  const float past_pi = wrap_angle_pi(std::nextafter(math::PI_F, 4.0f));
  HS_EXPECT_TRUE(past_pi < 0.0f && past_pi >= -math::PI_F);
  const float past_minus_pi = wrap_angle_pi(std::nextafter(-math::PI_F, -4.0f));
  HS_EXPECT_TRUE(past_minus_pi > 0.0f && past_minus_pi <= math::PI_F);
}

inline void test_clamp_finite_bounds_backend_parity() {
  for (float lo : {-2.0f, 0.0f, 0.5f})
    for (float hi : {1.0f, 2.0f})
      for (float value : {-INFINITY, -3.0f, 0.5f, 3.0f, INFINITY, NAN}) {
        volatile float input = value;
        HS_EXPECT_EQ(hs::clamp(input, lo, hi),
                     __builtin_fmaxf(lo, __builtin_fminf(input, hi)));
      }
}

// Clamp-before-cast / NaN-saturation checks that must also hold under
// -ffast-math -fno-finite-math-only; fastmath_clamp_check.cpp iterates this list.
#define HS_FASTMATH_CLAMP_TESTS(X)                                             \
  X(test_clamp_finite_bounds_backend_parity)                                   \
  X(test_blend_alpha_clamps_before_cast)                                       \
  X(test_pixel_scale_clamps_before_cast)                                       \
  X(test_gradient_get_clamps_out_of_range)                                     \
  X(test_generative_palette_get_clamps_out_of_range)                           \
  X(test_hue_rotation_lut_clamps_out_of_range_value)                           \
  X(test_generative_palette_get_nan_saturates_to_endpoint)

inline void test_lms_transform_pair_matches_scalar() {
  const float matrices[][9] = {
      {1, 0, 0, 0, 1, 0, 0, 0, 1},
      {0.8f, -0.3f, 0.5f, 0.2f, 1.2f, -0.4f, -0.1f, 0.3f, 0.8f}};
  const float samples[][3] = {{0, 0, 0},
                              {1, 1, 1},
                              {0.2f, 0.7f, 0.4f},
                              {-0.2f, 1.2f, 0.5f},
                              {0.9f, 0.1f, 0.7f}};
  for (const auto &matrix : matrices) {
    for (const auto &a : samples) {
      for (const auto &b : samples) {
        float scalar[6], paired[6];
        lms_cbrt_transform_rgb(matrix, a[0], a[1], a[2], scalar[0], scalar[1],
                               scalar[2]);
        lms_cbrt_transform_rgb(matrix, b[0], b[1], b[2], scalar[3], scalar[4],
                               scalar[5]);
        lms_cbrt_transform_rgb2(matrix, a[0], a[1], a[2], b[0], b[1], b[2],
                                paired[0], paired[1], paired[2], paired[3],
                                paired[4], paired[5]);
        HS_EXPECT_EQ(std::memcmp(scalar, paired, sizeof(scalar)), 0);
        lms_cbrt_transform_rgb_lut(matrix, a[0], a[1], a[2], scalar[0],
                                   scalar[1], scalar[2]);
        lms_cbrt_transform_rgb_lut(matrix, b[0], b[1], b[2], scalar[3],
                                   scalar[4], scalar[5]);
        lms_cbrt_transform_rgb2_lut(matrix, a[0], a[1], a[2], b[0], b[1], b[2],
                                    paired[0], paired[1], paired[2], paired[3],
                                    paired[4], paired[5]);
#if defined(HS_TEST_FAST_MATH)
        for (int channel = 0; channel < 6; ++channel)
          HS_EXPECT_NEAR(scalar[channel], paired[channel], 2e-6f);
#else
        HS_EXPECT_EQ(std::memcmp(scalar, paired, sizeof(scalar)), 0);
#endif
      }
    }
  }
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every color-module test and reports the aggregate result.
 * @return The module's failure count from end_module() (0 on success).
 */
inline int run_color_tests() {
  hs_test::ModuleFixture fixture("color");
  HS_RUN(test_pixel_quarter_accumulation_rounds_per_sample);
  HS_RUN(test_lms_transform_pair_matches_scalar);
  HS_RUN(test_baked_palette_storage_and_views);
  HS_RUN(test_baked_palette_constant_alpha);
  HS_RUN(test_lerp16_endpoints);
  HS_RUN(test_lerp16_midpoint);
  HS_RUN(test_lerp16_rounds_to_nearest);
  HS_RUN(test_color4_lerp_straight_alpha);
  HS_RUN(test_blend_outputs_denormal_alpha);
  HS_RUN(test_blend_outputs_tiny_normal_alpha);
  HS_RUN(test_blend_outputs_endpoints_verbatim);
  HS_RUN(test_wrap_angle_pi_large_arguments);
  HS_RUN(test_wrap_angle_pi_half_turn_keeps_sign);
  HS_RUN(test_lerp16_full_range_correct);
  HS_RUN(test_lerp16_bounded);

  HS_RUN(test_blend_add_packed_lane_layout);

  HS_RUN(test_oklab_roundtrip);
  HS_RUN(test_oklch_roundtrip);
  HS_RUN(test_oklab_reference_triples);
  HS_RUN(test_oklch_gray_is_achromatic);
  HS_RUN(test_lerp_oklch_achromatic_hue);
  HS_RUN(test_lerp_oklch_shortest_arc_midpoint);
  HS_RUN(test_lerp_oklch_endpoints);
  HS_RUN(test_lerp_oklch_extrapolation_clamped);
  HS_RUN(test_oklch_to_pixel_saturates_and_preserves_in_gamut);
  HS_RUN(test_gamut_clip_preserves_hue);
  HS_RUN(test_gamut_bracket_refine_out_of_gamut_lower_bound);
  HS_RUN(test_gamut_direction_lookup_matches_angle);
  HS_RUN(test_gamut_refine_matrices_match_the_conversions);
  HS_RUN(test_gamut_master_clip_lands_on_first_exit);
  HS_RUN(test_gamut_continuous_chroma_is_smooth_and_in_gamut);
  HS_RUN(test_gamut_lut_clip_lands_on_first_exit);
  HS_RUN(test_gamut_lut_downsample_preserves_bracket);
  HS_RUN(test_gamut_lut_release_and_passthrough);
  HS_RUN(test_gamut_cell_nonfinite_coordinates);
  HS_RUN(test_gamut_lut_boundary_scale_rounding);
  HS_RUN(test_configure_arenas_releases_gamut_lut);
  HS_RUN(test_oklch_to_pixel_holds_hue_out_of_gamut);

  HS_RUN(test_hue_sincos_matches_libm);
  HS_RUN(test_hue_rotate_preserves_gray);
  HS_RUN(test_hue_rotate_full_turn_identity);
  HS_RUN(test_hue_rotate_full_turn_in_steps_holds_hue_and_chroma);
  HS_RUN(test_hue_rotate_base_matches_direct);

  HS_RUN(test_srgb_to_linear_endpoints);
  HS_RUN(test_linear_to_srgb8_decode_matches_lut);
  HS_RUN(test_min_encodable_alpha_is_the_encode_floor);
  HS_RUN(test_srgb_linear_lut_vs_float_reference);
  HS_RUN(test_srgb_linear_roundtrip_lut);
  HS_RUN(test_srgb_to_linear_interp_recovers_subpixel_precision);
  HS_RUN(test_srgb_linear_roundtrip_float);

  HS_RUN(test_gradient_endpoints);
  HS_RUN(test_gradient_in_range_valid_and_monotone);
  HS_RUN(test_gradient_solid_color);
  HS_RUN(test_gradient_interpolates_between_entries);
  HS_RUN(test_gradient_first_stop_offset_flat_fills_prefix);
  HS_RUN(test_gradient_three_stops_interior_and_flanks);
  HS_RUN(test_gradient_hard_stop_is_abrupt);

  HS_RUN(test_baked_palette_matches_source_endpoints);
  HS_RUN(test_baked_palette_in_range);
  HS_RUN(test_baked_palette_rebake_samples_closed_interval);
  HS_RUN(test_baked_palette_samplers_match_get);
  HS_RUN(test_baked_palette_clone_from_matches_source);
  HS_RUN(test_dot_key_inverts_dot_keyed_coordinate);
  HS_RUN(test_dot_keyed_bake_round_trips_through_dot_key);
  HS_RUN(test_bake_palette_blend_nan_weight_stays_finite);
  HS_RUN(test_baked_palette_rebake_crossfade);
  HS_RUN(test_step_wipe_rebake_skips_arming_then_decrements);
  HS_RUN(test_palette_wipe_arm_step_cadence);

  HS_RUN(test_procedural_palette_cosine);
  HS_RUN(test_mutating_palette_blends_endpoints);
  HS_RUN(test_generative_palette_deterministic);
  HS_RUN(test_effect_palette_recipe_roster);
  HS_RUN(test_generative_palette_recipe_validation);
  HS_RUN(test_generative_palette_canonical_ignores_inactive_fields);
  HS_RUN(test_generative_palette_input_window);
  HS_RUN(test_generative_palette_resolves_axes_and_harmony);
  HS_RUN(test_generative_palette_hue_torsion);
  HS_RUN(test_generative_palette_blue_cusp_is_continuous);
  HS_RUN(test_generative_palette_local_gamut_stays_in_gamut);
  HS_RUN(test_generative_palette_domain_invariants);
  HS_RUN(test_generative_palette_morph_policy_contracts);
  HS_RUN(test_generative_palette_snapshot_lerp);
  HS_RUN(test_generative_palette_lerp_accumulates_segment_deltas);
  HS_RUN(test_generative_palette_snapshot_keeps_faint_chroma_chromatic);
  HS_RUN(test_generative_palette_snapshot_keeps_absolute_gray_achromatic);
  HS_RUN(test_generative_palette_lerp_gray_key_keeps_neighbor_hue);
  HS_RUN(test_generative_palette_lerp_target_aliases_this);
  HS_RUN(test_generative_palette_snapshot_lerp_closes_loop);
  HS_RUN(test_generative_palette_cartesian_path_neutralizes_midpoint);
  HS_RUN(test_generative_palette_rejects_unavailable_path_minimum);
  HS_RUN(test_generative_palette_absolute_basis_canonicalizes_headroom);
  HS_RUN(test_generative_palette_morph_compatible);
  HS_RUN(test_generative_palette_lerp_mixed_curves_continuous);
  HS_RUN(test_generative_palette_lerp_interpolates_loop_seam);
  HS_RUN(test_palette_cycler_arena_components);
  HS_RUN(test_palette_cycler_key_morph_cycle);
  HS_RUN(test_palette_cycler_heterogeneous_crossfade);
  HS_RUN(test_palette_cycler_pause_and_static);
  HS_RUN(test_palette_cycler_zero_dwell_chains_fades);
  HS_RUN(test_palette_cycler_generated_cycle);
  HS_RUN(test_palette_cycler_generated_display_cadence);
  HS_RUN(test_palette_cycler_sampled_restore_and_catchup);
  HS_RUN(test_generated_palette_bank_routes_and_rechromas);
  HS_RUN(test_palette_cycler_generated_chroma_keeps_morph);
  HS_RUN(test_palette_cycler_hidden_advance_catches_up);
  HS_RUN(test_palette_cycler_paused_hidden_advance_catches_up);
  HS_RUN(test_palette_cycler_roster_hidden_advance_catches_up);
  HS_RUN(test_palette_cycler_bake_generation);
  HS_RUN(test_standalone_palette_rotations_morph_compatible);
  HS_RUN(test_palette_modifiers);
  HS_RUN(test_noise_warp_modifier);
  HS_RUN(test_drift_modifier);
  HS_RUN(test_hue_spin_shade);
  HS_RUN(test_hue_wobble_shade);
  HS_RUN(test_sparkle_shade);
  HS_RUN(test_chroma_pulse_shade);
  HS_RUN(test_lightness_grain_shade);
  HS_RUN(test_iridescent_shade);
  HS_RUN(test_static_palette_composition);
  HS_RUN(test_palette_wrappers);
  HS_RUN(test_palette_shade_coord_policy);
  HS_RUN(test_color_modulation_mechanics);
  HS_RUN(test_pattern_modulated_palette);
  HS_RUN(test_sphere_lut_seamless_across_faces);
  HS_RUN(test_sphere_lut_all_edges);
  HS_RUN(test_sphere_noise_bake_qualification);
  HS_RUN(test_sphere_face_linear_table);
  HS_RUN(test_sphere_bake_cache);
  HS_RUN(test_sphere_bakes_match_reference);

  // Clamp-before-cast / NaN-saturation checks, shared with the fast-math pass.
  const int before_clamp = hs_test::stats().passed + hs_test::stats().failed;
#define HS_RUN_CLAMP_TEST(fn) fn();
  HS_FASTMATH_CLAMP_TESTS(HS_RUN_CLAMP_TEST)
#undef HS_RUN_CLAMP_TEST
  HS_EXPECT_GT(hs_test::stats().passed + hs_test::stats().failed - before_clamp,
               0);

  return fixture.result();
}

} // namespace color_tests
} // namespace hs_test
