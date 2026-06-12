/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Direct unit tests for core/styles.h — the Feedback::Style POD: named presets,
 * scalar lerp with function-pointer/discrete snapping, the transform functions
 * (identity/noise/melt warp, plain/hue fade), and sync_noise().
 *
 * Self-contained header. run_styles_tests() returns the module failure count.
 */
#pragma once

#include "core/styles.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace styles_tests {

// --- Named presets ----------------------------------------------------------

// Spot-check that the constexpr preset factories carry their documented scalar
// values and wire up the expected space/color transforms.
inline void test_named_presets() {
  Feedback::Style smoke = Feedback::Style::Smoke();
  HS_EXPECT_NEAR(smoke.fade, 0.9f, 1e-6f);
  HS_EXPECT_NEAR(smoke.frequency, 0.42f, 1e-6f);
  HS_EXPECT_TRUE(smoke.space_fn == &Feedback::noise_warp);
  HS_EXPECT_TRUE(smoke.color_fn == &Feedback::hue_fade);

  // The two melt-based presets use melt_warp; Swirling drops the hue shift.
  HS_EXPECT_TRUE(Feedback::Style::Melting().space_fn == &Feedback::melt_warp);
  HS_EXPECT_TRUE(Feedback::Style::Swirling().space_fn == &Feedback::melt_warp);
  HS_EXPECT_TRUE(Feedback::Style::Swirling().color_fn == &Feedback::plain_fade);
}

// --- lerp -------------------------------------------------------------------

// Style::lerp interpolates scalar fields linearly but snaps discrete fields
// (transform pointers, downsample) to b at t >= 0.5, and never overwrites the
// subject's own bound noise state.
inline void test_lerp_scalars_and_snapping() {
  NoiseParams na;
  Feedback::Style a{};
  a.fade = 0.0f;
  a.hue_shift = 0.0f;
  a.amplitude = 0.0f;
  a.frequency = 0.0f;
  a.speed = 0.0f;
  a.scale = 0.0f;
  a.space_fn = &Feedback::identity_warp;
  a.color_fn = &Feedback::plain_fade;
  a.downsample = 2;
  a.noise = &na;

  Feedback::Style b{};
  b.fade = 1.0f;
  b.hue_shift = 1.0f;
  b.amplitude = 1.0f;
  b.frequency = 1.0f;
  b.speed = 1.0f;
  b.scale = 1.0f;
  b.space_fn = &Feedback::noise_warp;
  b.color_fn = &Feedback::hue_fade;
  b.downsample = 8;
  b.noise = nullptr;

  NoiseParams subj;
  Feedback::Style mid{};
  mid.noise = &subj; // the subject's own bound state, distinct from a/b
  mid.lerp(a, b, 0.5f);
  // Scalars interpolate linearly.
  HS_EXPECT_NEAR(mid.fade, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(mid.hue_shift, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(mid.scale, 0.5f, 1e-6f);
  // Discrete fields snap to b at t >= 0.5.
  HS_EXPECT_TRUE(mid.space_fn == &Feedback::noise_warp);
  HS_EXPECT_TRUE(mid.color_fn == &Feedback::hue_fade);
  HS_EXPECT_EQ(mid.downsample, 8);
  // noise is the subject's bound state: lerping presets must leave it alone,
  // never pulling a.noise (&na) or b.noise (nullptr) over it.
  HS_EXPECT_TRUE(mid.noise == &subj);

  // Just below the midpoint the discrete fields stay on a.
  Feedback::Style lo{};
  lo.lerp(a, b, 0.4f);
  HS_EXPECT_TRUE(lo.space_fn == &Feedback::identity_warp);
  HS_EXPECT_TRUE(lo.color_fn == &Feedback::plain_fade);
  HS_EXPECT_EQ(lo.downsample, 2);
}

// --- Transform functions ----------------------------------------------------

// identity_warp returns its input direction unchanged.
inline void test_identity_warp() {
  Feedback::Style s{};
  Vector v(0.3f, -0.4f, 0.866f);
  Vector out = Feedback::identity_warp(v, s);
  HS_EXPECT_NEAR(out.x, v.x, 1e-6f);
  HS_EXPECT_NEAR(out.y, v.y, 1e-6f);
  HS_EXPECT_NEAR(out.z, v.z, 1e-6f);
}

// With no NoiseParams bound, noise_warp has nothing to sample and must pass the
// direction through unchanged.
inline void test_noise_warp_null_is_identity() {
  Feedback::Style s{};
  s.noise = nullptr; // unbound → no distortion
  Vector v = Vector(0.6f, 0.4f, 0.69f).normalized();
  Vector out = Feedback::noise_warp(v, s);
  HS_EXPECT_NEAR(out.x, v.x, 1e-6f);
  HS_EXPECT_NEAR(out.y, v.y, 1e-6f);
  HS_EXPECT_NEAR(out.z, v.z, 1e-6f);
}

// melt_warp slerps a direction toward the north pole at a rate set by speed,
// preserving unit length. With noise disabled, an equator point should rise.
inline void test_melt_warp_drifts_toward_north() {
  Feedback::Style s{};
  s.speed = 1.0f;     // positive drip rate
  s.amplitude = 0.0f; // no noise wobble
  s.noise = nullptr;  // skip the noise branch
  Vector v(1.0f, 0.0f, 0.0f); // on the equator (y = 0)
  Vector out = Feedback::melt_warp(v, s);
  // Slerp toward the north pole pulls y upward and shrinks x.
  HS_EXPECT_TRUE(out.y > 0.0f);
  HS_EXPECT_TRUE(out.x < 1.0f);
  HS_EXPECT_NEAR(out.length(), 1.0f, 1e-4f); // slerp preserves unit length
}

// plain_fade scales each channel by the fade factor with no hue change.
inline void test_plain_fade_scales_linearly() {
  Pixel p(10000, 20000, 30000);
  Feedback::Style s{};
  Pixel out = Feedback::plain_fade(p, 0.5f, s);
  HS_EXPECT_EQ(out.r, 5000);
  HS_EXPECT_EQ(out.g, 10000);
  HS_EXPECT_EQ(out.b, 15000);
}

// With a zero hue shift, hue_fade dims a gray pixel while keeping it gray. Gray
// has no chroma to rotate, so this avoids depending on OKLCH round-trip
// precision.
inline void test_hue_fade_zero_shift_preserves_gray() {
  Pixel gray(20000, 20000, 20000);
  Feedback::Style s{};
  s.hue_shift = 0.0f;
  Pixel out = Feedback::hue_fade(gray, 0.5f, s);
  HS_EXPECT_TRUE(out.r < gray.r); // dimmed
  HS_EXPECT_TRUE(out.r > 0);
  // Stays gray within a small round-trip tolerance.
  HS_EXPECT_TRUE(std::abs((int)out.r - (int)out.g) < 64);
  HS_EXPECT_TRUE(std::abs((int)out.g - (int)out.b) < 64);
}

// --- sync_noise -------------------------------------------------------------

// sync_noise copies the Style's noise-related scalars into its bound
// NoiseParams, and is a safe no-op when no NoiseParams is bound.
inline void test_sync_noise_pushes_scalars() {
  NoiseParams np;
  Feedback::Style s{};
  s.amplitude = 7.0f;
  s.frequency = 0.33f;
  s.speed = 2.5f;
  s.scale = 9.0f;
  s.noise = &np;
  s.sync_noise();
  HS_EXPECT_NEAR(np.amplitude, 7.0f, 1e-6f);
  HS_EXPECT_NEAR(np.frequency, 0.33f, 1e-6f);
  HS_EXPECT_NEAR(np.speed, 2.5f, 1e-6f);
  HS_EXPECT_NEAR(np.scale, 9.0f, 1e-6f);

  // Unbound sync is a safe no-op (must not dereference null).
  Feedback::Style unbound{};
  unbound.noise = nullptr;
  unbound.sync_noise();
  HS_EXPECT_TRUE(true);
}

// Run every styles test case; returns the module's failure count.
inline int run_styles_tests() {
  auto scope = hs_test::begin_module("styles");

  test_named_presets();
  test_lerp_scalars_and_snapping();
  test_identity_warp();
  test_noise_warp_null_is_identity();
  test_melt_warp_drifts_toward_north();
  test_plain_fade_scales_linearly();
  test_hue_fade_zero_shift_preserves_gray();
  test_sync_noise_pushes_scalars();

  return hs_test::end_module(scope);
}

} // namespace styles_tests
} // namespace hs_test
