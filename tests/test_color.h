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
#include "core/color/effect_palette_recipes.h"
#include "core/color/generative_palette.h"
#include "core/color/noise_hue_palette.h"
#include "core/color/noise_shimmer_palette.h"
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
 * @return Process exit code from hs_test::end_module: 0 on success, non-zero on
 *         any failure.
 */
inline int run_color_tests() {
  hs_test::ModuleFixture fixture("color");
  test_pixel_quarter_accumulation_rounds_per_sample();
  test_lms_transform_pair_matches_scalar();
  test_baked_palette_storage_and_views();
  test_lerp16_endpoints();
  test_lerp16_midpoint();
  test_lerp16_rounds_to_nearest();
  test_color4_lerp_straight_alpha();
  test_blend_outputs_denormal_alpha();
  test_blend_outputs_tiny_normal_alpha();
  test_blend_outputs_endpoints_verbatim();
  test_wrap_angle_pi_large_arguments();
  test_wrap_angle_pi_half_turn_keeps_sign();
  test_lerp16_full_range_correct();
  test_lerp16_bounded();

  test_blend_add_packed_lane_layout();

  test_oklab_roundtrip();
  test_oklch_roundtrip();
  test_oklab_reference_triples();
  test_oklch_gray_is_achromatic();
  test_lerp_oklch_achromatic_hue();
  test_lerp_oklch_shortest_arc_midpoint();
  test_lerp_oklch_endpoints();
  test_lerp_oklch_extrapolation_clamped();
  test_oklch_to_pixel_saturates_and_preserves_in_gamut();
  test_gamut_clip_preserves_hue();
  test_gamut_bracket_refine_out_of_gamut_lower_bound();
  test_gamut_direction_lookup_matches_angle();
  test_gamut_refine_matrices_match_the_conversions();
  test_gamut_master_clip_lands_on_first_exit();
  test_gamut_continuous_chroma_is_smooth_and_in_gamut();
  test_gamut_lut_clip_lands_on_first_exit();
  test_gamut_lut_downsample_preserves_bracket();
  test_gamut_lut_release_and_passthrough();
  test_gamut_cell_nonfinite_coordinates();
  test_gamut_lut_boundary_scale_rounding();
  test_configure_arenas_releases_gamut_lut();
  test_oklch_to_pixel_holds_hue_out_of_gamut();

  test_hue_sincos_matches_libm();
  test_hue_rotate_preserves_gray();
  test_hue_rotate_full_turn_identity();
  test_hue_rotate_full_turn_in_steps_holds_hue_and_chroma();
  test_hue_rotate_base_matches_direct();

  test_srgb_to_linear_endpoints();
  test_linear_to_srgb_endpoints();
  test_linear_to_srgb8_decode_matches_lut();
  test_min_encodable_alpha_is_the_encode_floor();
  test_srgb_linear_lut_vs_float_reference();
  test_srgb_linear_roundtrip_lut();
  test_srgb_to_linear_interp_recovers_subpixel_precision();
  test_srgb_linear_roundtrip_float();

  test_gradient_endpoints();
  test_gradient_in_range_valid_and_monotone();
  test_gradient_solid_color();
  test_gradient_interpolates_between_entries();
  test_gradient_first_stop_offset_flat_fills_prefix();
  test_gradient_three_stops_interior_and_flanks();
  test_gradient_hard_stop_is_abrupt();

  test_baked_palette_matches_source_endpoints();
  test_baked_palette_in_range();
  test_baked_palette_rebake_samples_closed_interval();
  test_baked_palette_color_sampler_matches_get();
  test_baked_palette_clone_from_matches_source();
  test_dot_key_inverts_dot_keyed_coordinate();
  test_dot_keyed_bake_round_trips_through_dot_key();
  test_bake_palette_blend_nan_weight_stays_finite();
  test_baked_palette_rebake_crossfade();
  test_step_wipe_rebake_skips_arming_then_decrements();
  test_palette_wipe_arm_step_cadence();

  test_procedural_palette_cosine();
  test_mutating_palette_blends_endpoints();
  test_generative_palette_deterministic();
  test_effect_palette_recipe_roster();
  test_generative_palette_recipe_validation();
  test_generative_palette_canonical_ignores_inactive_fields();
  test_generative_palette_input_window();
  test_generative_palette_resolves_axes_and_harmony();
  test_generative_palette_hue_torsion();
  test_generative_palette_blue_cusp_is_continuous();
  test_generative_palette_local_gamut_stays_in_gamut();
  test_generative_palette_domain_invariants();
  test_generative_palette_morph_policy_contracts();
  test_generative_palette_snapshot_lerp();
  test_generative_palette_lerp_accumulates_segment_deltas();
  test_generative_palette_snapshot_keeps_faint_chroma_chromatic();
  test_generative_palette_snapshot_keeps_absolute_gray_achromatic();
  test_generative_palette_lerp_target_aliases_this();
  test_generative_palette_snapshot_lerp_closes_loop();
  test_generative_palette_cartesian_path_neutralizes_midpoint();
  test_generative_palette_rejects_unavailable_path_minimum();
  test_generative_palette_absolute_basis_canonicalizes_headroom();
  test_generative_palette_morph_compatible();
  test_generative_palette_lerp_mixed_curves_continuous();
  test_generative_palette_lerp_interpolates_loop_seam();
  test_palette_cycler_arena_components();
  test_palette_cycler_key_morph_cycle();
  test_palette_cycler_heterogeneous_crossfade();
  test_palette_cycler_pause_and_static();
  test_palette_cycler_zero_dwell_chains_fades();
  test_palette_cycler_generated_cycle();
  test_generated_palette_bank_routes_and_rechromas();
  test_palette_cycler_generated_chroma_keeps_morph();
  test_palette_cycler_hidden_advance_catches_up();
  test_palette_cycler_roster_hidden_advance_catches_up();
  test_palette_cycler_bake_generation();
  test_standalone_palette_rotations_morph_compatible();
  test_palette_modifiers();
  test_noise_warp_modifier();
  test_drift_modifier();
  test_hue_spin_shade();
  test_hue_wobble_shade();
  test_sparkle_shade();
  test_chroma_pulse_shade();
  test_lightness_grain_shade();
  test_iridescent_shade();
  test_static_palette_composition();
  test_palette_wrappers();
  test_palette_shade_coord_policy();
  test_noise_hue_palette();
  test_noise_hue_palette_direct();
  test_noise_shimmer_palette();
  test_hue_noise_lut_seamless_across_faces();
  test_hue_noise_bake_cache();
  test_hue_noise_paired_bakes_match_reference();

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
