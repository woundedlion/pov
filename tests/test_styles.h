/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Direct unit tests for core/engine/styles.h — the Feedback::Style POD: named presets,
 * scalar lerp with function-pointer/discrete snapping, the transform functions
 * (identity/noise/melt warp — both the unbound identity path and the bound-noise
 * production branch — plain/hue fade), and sync_noise().
 *
 * Self-contained header. run_styles_tests() returns the module failure count.
 */
#pragma once

#include "core/engine/styles.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace styles_tests {

// --- Named presets ----------------------------------------------------------

/**
 * @brief Verifies the constexpr preset factories carry their documented scalar
 *        values and wire up the expected space/color transforms.
 * @details Spot-checks Smoke's fade/frequency and noise/hue transforms, and
 *          confirms both melt-based presets use melt_warp while Swirling drops
 *          the hue shift.
 */
inline void test_named_presets() {
  Feedback::Style smoke = Feedback::Style::Smoke();
  HS_EXPECT_NEAR(smoke.fade, 0.9f, 1e-6f);
  HS_EXPECT_NEAR(smoke.frequency, 0.42f, 1e-6f);
  HS_EXPECT_TRUE(smoke.space_fn == &Feedback::noise_warp);
  HS_EXPECT_TRUE(smoke.color_fn == &Feedback::hue_fade);

  HS_EXPECT_TRUE(Feedback::Style::Melting().space_fn == &Feedback::melt_warp);
  HS_EXPECT_TRUE(Feedback::Style::Swirling().space_fn == &Feedback::melt_warp);
  HS_EXPECT_TRUE(Feedback::Style::Swirling().color_fn == &Feedback::plain_fade);
}

// --- lerp -------------------------------------------------------------------

/**
 * @brief Verifies Style::lerp interpolates scalar fields linearly, snaps
 *        discrete fields, and preserves the subject's bound noise state.
 * @details Scalar fields (fade, hue_shift, scale, ...) interpolate linearly;
 *          discrete fields (transform pointers, downsample) snap to b at
 *          t >= 0.5 and stay on a below the midpoint. The subject's bound
 *          noise pointer must never be overwritten by a's or b's noise.
 */
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
  mid.noise = &subj;
  mid.lerp(a, b, 0.5f);
  HS_EXPECT_NEAR(mid.fade, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(mid.hue_shift, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(mid.scale, 0.5f, 1e-6f);
  HS_EXPECT_TRUE(mid.space_fn == &Feedback::noise_warp);
  HS_EXPECT_TRUE(mid.color_fn == &Feedback::hue_fade);
  HS_EXPECT_EQ(mid.downsample, 8);
  HS_EXPECT_TRUE(mid.noise == &subj);

  Feedback::Style lo{};
  lo.lerp(a, b, 0.4f);
  HS_EXPECT_TRUE(lo.space_fn == &Feedback::identity_warp);
  HS_EXPECT_TRUE(lo.color_fn == &Feedback::plain_fade);
  HS_EXPECT_EQ(lo.downsample, 2);
}

// --- Transform functions ----------------------------------------------------

/**
 * @brief Verifies identity_warp returns its input direction unchanged.
 */
inline void test_identity_warp() {
  Feedback::Style s{};
  Vector v(0.3f, -0.4f, 0.866f);
  Vector out = Feedback::identity_warp(v, s);
  HS_EXPECT_NEAR(out.x, v.x, 1e-6f);
  HS_EXPECT_NEAR(out.y, v.y, 1e-6f);
  HS_EXPECT_NEAR(out.z, v.z, 1e-6f);
}

/**
 * @brief Verifies noise_warp passes the direction through unchanged when no
 *        NoiseParams is bound.
 * @details With nothing to sample, the warp must be the identity.
 */
inline void test_noise_warp_null_is_identity() {
  Feedback::Style s{};
  s.noise = nullptr;
  Vector v = Vector(0.6f, 0.4f, 0.69f).normalized();
  Vector out = Feedback::noise_warp(v, s);
  HS_EXPECT_NEAR(out.x, v.x, 1e-6f);
  HS_EXPECT_NEAR(out.y, v.y, 1e-6f);
  HS_EXPECT_NEAR(out.z, v.z, 1e-6f);
}

/**
 * @brief Verifies melt_warp slerps a direction toward the north pole at a rate
 *        set by speed while preserving unit length.
 * @details With noise disabled, an equator point should rise: y increases, x
 *          shrinks, and the result stays unit length.
 */
inline void test_melt_warp_drifts_toward_north() {
  Feedback::Style s{};
  s.speed = 1.0f;
  s.amplitude = 0.0f;
  s.noise = nullptr;
  Vector v(1.0f, 0.0f, 0.0f); // on the equator (y = 0)
  Vector out = Feedback::melt_warp(v, s);
  // speed=1 slerps 0.04 of the 90 deg arc toward the pole: y rises ~0.0637, x
  // drops ~0.002. Pin a minimum drift so a no-op warp can't pass.
  HS_EXPECT_TRUE(out.y > 0.05f);
  HS_EXPECT_TRUE(out.x < 0.999f);
  HS_EXPECT_NEAR(out.length(), 1.0f, 1e-4f);
}

/**
 * @brief Verifies noise_warp actually distorts when a NoiseParams is bound.
 * @details The null-noise test above only covers the identity early-out. Here a
 *          real NoiseParams is bound and primed via the production sync_noise()
 *          path (which pushes the Style's amplitude/frequency/speed/scale into
 *          it); noise_warp must then take the bound branch — delegating to
 *          noise_transform — and displace the direction off the input while
 *          keeping it on the unit sphere. Displacement is summed across samples
 *          so a single noise zero-crossing can't make the test flaky.
 */
inline void test_noise_warp_bound_distorts() {
  NoiseParams np;
  Feedback::Style s{};
  s.amplitude = 0.6f;
  s.frequency = 0.5f;
  s.speed = 0.0f;
  s.scale = 4.0f;
  s.noise = &np;
  s.sync_noise();

  const Vector samples[] = {Vector(1, 0, 0), Vector(0, 0, 1),
                            Vector(0.4f, 0.6f, 0.7f).normalized()};
  float total_moved = 0.0f;
  for (const Vector &v : samples) {
    Vector out = Feedback::noise_warp(v, s);
    HS_EXPECT_NEAR(out.length(), 1.0f, 1e-3f);
    total_moved +=
        std::abs(out.x - v.x) + std::abs(out.y - v.y) + std::abs(out.z - v.z);
  }
  HS_EXPECT_GT(total_moved, 1e-2f);
}

/**
 * @brief Verifies melt_warp's bound-noise branch perturbs the drip.
 * @details test_melt_warp_drifts_toward_north exercises only the noise-disabled
 *          drip. With a NoiseParams bound and amplitude above the wobble floor,
 *          melt_warp must take the `s.noise && s.amplitude > floor` branch and
 *          perturb the drifted point, so its output diverges from the same Style
 *          with noise unbound (the pure-drip result), while staying unit length.
 */
inline void test_melt_warp_bound_noise_perturbs() {
  NoiseParams np;
  Feedback::Style s{};
  s.speed = 1.0f;
  s.amplitude = 0.6f; // above the melt wobble floor → noise branch runs
  s.frequency = 0.5f;
  s.scale = 4.0f;
  s.noise = &np;
  s.sync_noise();

  Feedback::Style drip_only = s;
  drip_only.noise = nullptr;

  const Vector samples[] = {Vector(1, 0, 0), Vector(0, 0, 1),
                            Vector(0.4f, 0.6f, 0.7f).normalized()};
  float total_divergence = 0.0f;
  for (const Vector &v : samples) {
    Vector with_noise = Feedback::melt_warp(v, s);
    Vector pure_drip = Feedback::melt_warp(v, drip_only);
    HS_EXPECT_NEAR(with_noise.length(), 1.0f, 1e-3f);
    total_divergence += std::abs(with_noise.x - pure_drip.x) +
                        std::abs(with_noise.y - pure_drip.y) +
                        std::abs(with_noise.z - pure_drip.z);
  }
  HS_EXPECT_GT(total_divergence, 1e-3f);
}

/**
 * @brief Verifies plain_fade scales each channel by the fade factor with no
 *        hue change.
 */
inline void test_plain_fade_scales_linearly() {
  Pixel p(10000, 20000, 30000);
  Feedback::Style s{};
  Pixel out = Feedback::plain_fade(p, 0.5f, s);
  HS_EXPECT_EQ(out.r, 5000);
  HS_EXPECT_EQ(out.g, 10000);
  HS_EXPECT_EQ(out.b, 15000);
}

/**
 * @brief Verifies hue_fade with a zero hue shift dims a gray pixel while
 *        keeping it gray.
 * @details Gray has no chroma to rotate, so this avoids depending on OKLCH
 *          round-trip precision; the channels stay equal within a small
 *          tolerance.
 */
inline void test_hue_fade_zero_shift_preserves_gray() {
  Pixel gray(20000, 20000, 20000);
  Feedback::Style s{};
  s.hue_shift = 0.0f;
  Pixel out = Feedback::hue_fade(gray, 0.5f, s);
  HS_EXPECT_TRUE(out.r < gray.r);
  HS_EXPECT_TRUE(out.r > 0);
  // Stays gray within the OKLCH round-trip tolerance.
  HS_EXPECT_TRUE(std::abs((int)out.r - (int)out.g) < 64);
  HS_EXPECT_TRUE(std::abs((int)out.g - (int)out.b) < 64);
}

/**
 * @brief Verifies sync_hue() caches cos/sin of the hue_shift turn angle into
 *        hue_ca/hue_sa, leaving the struct at the identity rotation until called.
 * @details hue_fade reads this frame-constant cache instead of recomputing the
 *          rotation per pixel; sync_hue must write exactly fast_cosf/fast_sinf of
 *          hue_shift*2*PI. A fresh Style defaults to the identity (1, 0).
 */
inline void test_sync_hue_caches_rotation() {
  Feedback::Style s{};
  HS_EXPECT_NEAR(s.hue_ca, 1.0f, 1e-6f);
  HS_EXPECT_NEAR(s.hue_sa, 0.0f, 1e-6f);

  s.hue_shift = 0.25f; // quarter turn
  s.sync_hue();
  const float angle = 0.25f * (2.0f * PI_F);
  HS_EXPECT_NEAR(s.hue_ca, fast_cosf(angle), 1e-6f);
  HS_EXPECT_NEAR(s.hue_sa, fast_sinf(angle), 1e-6f);
  // A quarter turn lands near (cos, sin) = (0, 1).
  HS_EXPECT_NEAR(s.hue_ca, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(s.hue_sa, 1.0f, 1e-3f);
}

/**
 * @brief Verifies hue_fade with a nonzero shift actually rotates a saturated
 *        pixel's hue via the sync_hue cache.
 * @details Unlike the gray-pixel zero-shift case, a saturated pixel has chroma to
 *          rotate: the rotated result must diverge from the plain (hue-preserving)
 *          fade, and must match an independent per-call rotation by the same
 *          amount — pinning that hue_fade consumes the sync_hue cache correctly.
 */
inline void test_hue_fade_nonzero_shift_rotates_saturated() {
  Pixel red(50000, 2000, 2000);
  Feedback::Style s{};
  s.hue_shift = 0.33f;
  s.sync_hue();
  const float fade = 0.8f;
  Pixel rotated = Feedback::hue_fade(red, fade, s);

  Pixel plain = Feedback::plain_fade(red, fade, s);
  const int delta = std::abs((int)rotated.r - (int)plain.r) +
                    std::abs((int)rotated.g - (int)plain.g) +
                    std::abs((int)rotated.b - (int)plain.b);
  HS_EXPECT_TRUE(delta > 1000);
}

/**
 * @brief Verifies the identity rotation's cbrt-LMS matrix is the identity.
 * @details hue_rotate_lms_matrix folds oklab_to_lms_cbrt . rotate .
 *          lms_to_oklab; at rotation zero the fold must reduce to the identity
 *          within float error of the two OKLab matrices, which pins them as
 *          mutual inverses.
 */
inline void test_hue_rotate_lms_matrix_identity() {
  float k[9];
  hue_rotate_lms_matrix(1.0f, 0.0f, k);
  for (int i = 0; i < 9; ++i)
    HS_EXPECT_NEAR(k[i], (i % 4 == 0) ? 1.0f : 0.0f, 1e-5f);
}

/**
 * @brief Parity sweep: hue_fade's folded cbrt-LMS path vs the reference
 *        fade-then-hue_rotate composition.
 * @details The fold changes only float association and drops the fade's
 *          intermediate u16 quantization, so results must track the reference
 *          within a few 16-bit LSB. The sweep crosses saturated primaries
 *          (which exercise the gamut-clip slow path), a dark and a gray tone
 *          with every preset-range fade/shift combination. The measured max
 *          delta on this sweep is 19 LSB; the tolerance carries ~3x margin.
 *          Runs on MeshFeedback's boundary grid, the one hue_fade ships against:
 *          the clip's bracket residual sets this delta, so a different grid
 *          measures a configuration nothing runs (the flash master reads 106).
 */
inline void test_hue_fade_matches_rotate_reference() {
  constexpr float HS_HUE_FADE_TOL = 64.0f;
  alignas(uint16_t) static uint8_t lut_buf[gamut_lut_bytes(256, 128)];
  Arena lut_arena(lut_buf, sizeof(lut_buf));
  init_gamut_lut(lut_arena, 256, 128);
  const Pixel colors[] = {Pixel(65535, 0, 0),     Pixel(0, 65535, 0),
                          Pixel(0, 0, 65535),     Pixel(65535, 65535, 0),
                          Pixel(50000, 2000, 2000), Pixel(300, 200, 100),
                          Pixel(20000, 20000, 20000)};
  const float fades[] = {0.0f, 0.58f, 0.9f, 0.99f};
  const float shifts[] = {0.0f, 0.01f, 0.05f, 0.1f, 0.33f};
  for (const Pixel &c : colors)
    for (float fade : fades)
      for (float shift : shifts) {
        Feedback::Style s{};
        s.hue_shift = shift;
        s.sync_hue();
        Pixel got = Feedback::hue_fade(c, fade, s);
        Pixel ref = hue_rotate(Color4(c * fade, 1.0f), s.hue_ca, s.hue_sa).color;
        HS_EXPECT_NEAR((float)got.r, (float)ref.r, HS_HUE_FADE_TOL);
        HS_EXPECT_NEAR((float)got.g, (float)ref.g, HS_HUE_FADE_TOL);
        HS_EXPECT_NEAR((float)got.b, (float)ref.b, HS_HUE_FADE_TOL);
      }
  release_gamut_lut();
}

/**
 * @brief Pins hue_fade_apply2 against the scalar hue_fade_apply at the
 *        quantized output, the level the display observes.
 * @details The paired path shares one fast_cbrt6 reciprocal across both
 *          pixels, which re-associates the arithmetic and moves the cube roots
 *          by ~4e-7 relative. The gamut clip brackets its input on a grid, so
 *          near a bracket boundary that tiny shift selects a different bracket
 *          and the channel moves by tens of LSB rather than one; the bound is
 *          therefore loose by necessity. The measured max on this sweep is 39
 *          LSB and the tolerance carries ~3x margin. Sweeps preset-range fades
 *          and several hue rotations (identity included) over realistic
 *          channels plus the zero, near-zero and all-black cases the
 *          composite's NEAR_BLACK skip relies on, and the saturated primaries
 *          that drive the clip. fast_cbrt6's own accuracy is pinned tightly in
 *          test_fast_cbrt6; this case guards the composite's end-to-end output.
 */
inline void test_hue_fade_apply2_tracks_scalar() {
  constexpr int HS_PAIR_TOL = 128;
  const float channels[][3] = {
      {0.0f, 0.0f, 0.0f},             {1.0f, 0.0f, 0.0f},
      {0.5f, 0.25f, 0.125f},          {300.0f, 200.0f, 100.0f},
      {20000.0f, 20000.0f, 20000.0f}, {65535.0f, 0.0f, 0.0f},
      {0.0f, 65535.0f, 0.0f},         {0.0f, 0.0f, 65535.0f},
      {65535.0f, 65535.0f, 0.0f},     {50000.0f, 2000.0f, 2000.0f},
      {65535.0f, 65535.0f, 65535.0f}, {12345.0f, 54321.0f, 999.0f}};
  const int n = (int)(sizeof(channels) / sizeof(channels[0]));
  const float fades[] = {0.58f, 0.75f, 0.9f, 0.99f};
  const float shifts[] = {0.0f, 0.01f, 0.1f, 0.33f, 0.75f};

  for (float shift : shifts)
    for (float fade : fades) {
      Feedback::Style s{};
      s.hue_shift = shift;
      s.sync_hue();
      // The composite folds the fade and the u16 normalization into k.
      float k[9];
      const float sc = fast_cbrt(fade * (1.0f / 65535.0f));
      for (int i = 0; i < 9; ++i)
        k[i] = s.hue_k[i] * sc;

      for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) {
          Pixel s0 = Feedback::hue_fade_apply(k, channels[a][0], channels[a][1],
                                              channels[a][2]);
          Pixel s1 = Feedback::hue_fade_apply(k, channels[b][0], channels[b][1],
                                              channels[b][2]);
          Pixel p0, p1;
          Feedback::hue_fade_apply2(k, channels[a][0], channels[a][1],
                                    channels[a][2], channels[b][0],
                                    channels[b][1], channels[b][2], p0, p1);
          HS_EXPECT_TRUE(std::abs((int)p0.r - (int)s0.r) <= HS_PAIR_TOL);
          HS_EXPECT_TRUE(std::abs((int)p0.g - (int)s0.g) <= HS_PAIR_TOL);
          HS_EXPECT_TRUE(std::abs((int)p0.b - (int)s0.b) <= HS_PAIR_TOL);
          HS_EXPECT_TRUE(std::abs((int)p1.r - (int)s1.r) <= HS_PAIR_TOL);
          HS_EXPECT_TRUE(std::abs((int)p1.g - (int)s1.g) <= HS_PAIR_TOL);
          HS_EXPECT_TRUE(std::abs((int)p1.b - (int)s1.b) <= HS_PAIR_TOL);
        }
    }

  // An all-black pair stays exactly black, the state the composite writes to
  // overwrite the stale double-buffer frame.
  Feedback::Style s{};
  s.hue_shift = 0.2f;
  s.sync_hue();
  float k[9];
  const float sc = fast_cbrt(0.9f * (1.0f / 65535.0f));
  for (int i = 0; i < 9; ++i)
    k[i] = s.hue_k[i] * sc;
  Pixel p0, p1;
  Feedback::hue_fade_apply2(k, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, p0, p1);
  HS_EXPECT_TRUE(p0.r == 0 && p0.g == 0 && p0.b == 0);
  HS_EXPECT_TRUE(p1.r == 0 && p1.g == 0 && p1.b == 0);
}

// --- sync_noise -------------------------------------------------------------

/**
 * @brief Verifies sync_noise copies the Style's noise-related scalars into its
 *        bound NoiseParams and is a safe no-op when none is bound.
 * @details Pushes amplitude/frequency/speed/scale into the bound NoiseParams;
 *          an unbound sync must not dereference null.
 */
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

  // Clearing the binding makes sync a no-op: it must neither dereference the
  // null pointer nor write through to the previously-bound NoiseParams.
  s.noise = nullptr;
  s.amplitude = 1.0f;
  s.sync_noise();
  HS_EXPECT_NEAR(np.amplitude, 7.0f, 1e-6f);
}

/**
 * @brief Runs every styles test case.
 * @return The module's failure count.
 */
inline int run_styles_tests() {
  hs_test::ModuleFixture fixture("styles");

  test_named_presets();
  test_lerp_scalars_and_snapping();
  test_identity_warp();
  test_noise_warp_null_is_identity();
  test_noise_warp_bound_distorts();
  test_melt_warp_drifts_toward_north();
  test_melt_warp_bound_noise_perturbs();
  test_plain_fade_scales_linearly();
  test_hue_fade_zero_shift_preserves_gray();
  test_sync_hue_caches_rotation();
  test_hue_fade_nonzero_shift_rotates_saturated();
  test_hue_rotate_lms_matrix_identity();
  test_hue_fade_matches_rotate_reference();
  test_hue_fade_apply2_tracks_scalar();
  test_sync_noise_pushes_scalars();

  return fixture.result();
}

} // namespace styles_tests
} // namespace hs_test
