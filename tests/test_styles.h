/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
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
 * @details Spot-checks Smoke's fade/frequency and noise/hue transforms and
 *          verifies the MeshFeedback melt presets exactly.
 */
inline void test_named_presets() {
  const auto expect_noise_hue_style =
      [](const Feedback::Style &style, float fade, float hue_shift,
         float amplitude, float frequency, float speed, float scale) {
        HS_EXPECT_NEAR(style.fade, fade, 1e-6f);
        HS_EXPECT_NEAR(style.hue_shift, hue_shift, 1e-6f);
        HS_EXPECT_NEAR(style.amplitude, amplitude, 1e-6f);
        HS_EXPECT_NEAR(style.frequency, frequency, 1e-6f);
        HS_EXPECT_NEAR(style.speed, speed, 1e-6f);
        HS_EXPECT_NEAR(style.scale, scale, 1e-6f);
        HS_EXPECT_TRUE(style.space_fn == &Feedback::noise_warp);
        HS_EXPECT_TRUE(style.color_fn == &Feedback::hue_fade);
      };
  const auto expect_melt_hue_style =
      [](const Feedback::Style &style, float fade, float hue_shift,
         float amplitude, float frequency, float speed, float scale) {
        HS_EXPECT_NEAR(style.fade, fade, 1e-6f);
        HS_EXPECT_NEAR(style.hue_shift, hue_shift, 1e-6f);
        HS_EXPECT_NEAR(style.amplitude, amplitude, 1e-6f);
        HS_EXPECT_NEAR(style.frequency, frequency, 1e-6f);
        HS_EXPECT_NEAR(style.speed, speed, 1e-6f);
        HS_EXPECT_NEAR(style.scale, scale, 1e-6f);
        HS_EXPECT_TRUE(style.space_fn == &Feedback::melt_warp);
        HS_EXPECT_TRUE(style.color_fn == &Feedback::hue_fade);
      };

  expect_noise_hue_style(Feedback::Style::ArcingLightning(), 0.5f, 0.14426951f,
                         3.27f, 0.09f, 1.5f, 50.0f);
  expect_noise_hue_style(Feedback::Style::SlowFire(), 0.8732f, 0.12316483f,
                         1.56f, 0.5297f, 0.1f, 50.0f);
  expect_noise_hue_style(Feedback::Style::EnergeticFire(), 0.8732f, 0.12316483f,
                         1.56f, 0.22087f, 0.9f, 50.0f);
  expect_noise_hue_style(Feedback::Style::SlowDust(), 0.83952f, 0.09546948f,
                         1.56f, 0.07237f, 0.6f, 50.0f);
  expect_noise_hue_style(Feedback::Style::WavyTrails(), 0.7257f, 0.22518973f,
                         1.95f, 0.01f, 5.0f, 50.0f);
  expect_melt_hue_style(Feedback::Style::MeltingHi(), 0.59004f, 0.18955015f,
                        4.38f, 0.06346f, 0.2f, 22.3554f);
  expect_melt_hue_style(Feedback::Style::MeltingLo(), 0.59004f, 0.18955015f,
                        1.95f, 0.06346f, 0.2f, 22.3554f);
  expect_noise_hue_style(Feedback::Style::Miasma(), 0.80586f, 0.234f, 2.61f,
                         0.05059f, 0.725f, 26.297501f);
  expect_noise_hue_style(Feedback::Style::LooseWormhole(), 0.7257f, 0.22519f,
                         11.25f, 0.01f, 0.0f, 10.1798f);
  expect_noise_hue_style(Feedback::Style::TightWormhole(), 0.7257f, 0.22519f,
                         6.42f, 0.01f, 0.0f, 7.8844f);
  expect_noise_hue_style(Feedback::Style::WigglingWormhole(), 0.7257f, 0.22519f,
                         7.11f, 0.01f, 0.0f, 29.1917f);

  Feedback::Style smoke = Feedback::Style::Smoke();
  HS_EXPECT_NEAR(smoke.fade, 0.9f, 1e-6f);
  HS_EXPECT_NEAR(smoke.hue_shift, 0.09491219f, 1e-6f);
  HS_EXPECT_NEAR(smoke.frequency, 0.42f, 1e-6f);
  HS_EXPECT_TRUE(smoke.space_fn == &Feedback::noise_warp);
  HS_EXPECT_TRUE(smoke.color_fn == &Feedback::hue_fade);

  HS_EXPECT_TRUE(Feedback::Style::MeltingHi().space_fn == &Feedback::melt_warp);
}

/**
 * @brief Verifies named presets retain their original per-frame hue rotations.
 * @details Drives the pinned shift through sync_hue() and recovers the angle
 *          from the cached cos/sin, so the preset table is compared against what
 *          the production path actually rotates by rather than against a
 *          restatement of its formula.
 */
inline void test_named_presets_preserve_frame_hue() {
  struct Case {
    Feedback::Style style;
    float frame_shift;
  };
  const Case cases[] = {
      {Feedback::Style::ArcingLightning(), 0.1f},
      {Feedback::Style::SlowFire(), 0.0167f},
      {Feedback::Style::EnergeticFire(), 0.0167f},
      {Feedback::Style::Smoke(), 0.01f},
      {Feedback::Style::SlowDust(), 0.0167f},
      {Feedback::Style::WavyTrails(), 0.0722f},
      {Feedback::Style::MeltingHi(), 0.1f},
      {Feedback::Style::MeltingLo(), 0.1f},
      {Feedback::Style::Miasma(), 0.05050779f},
      {Feedback::Style::LooseWormhole(), 0.07220009f},
      {Feedback::Style::TightWormhole(), 0.07220009f},
      {Feedback::Style::WigglingWormhole(), 0.07220009f},
  };
  for (const Case &c : cases) {
    Feedback::Style s = c.style;
    s.sync_hue();
    // The cached angle comes from fast trig, which lands within 2.3e-4 turns
    // over this table; the bound carries ~4x margin.
    HS_EXPECT_NEAR(std::atan2(s.hue_sa, s.hue_ca) / (2.0f * PI_F),
                   c.frame_shift, 1e-3f);
  }
}

/**
 * @brief Pins sync_hue's per-frame rotation to hue_shift turns per e-fold of
 *        feedback brightness decay.
 * @details hue_shift is defined as turns of hue rotation per e-fold decrease in
 *          brightness, so a fade of exp(-n) — which loses exactly n e-folds per
 *          frame — must rotate by n * hue_shift turns. The expected angle comes
 *          from that rate and std cos/sin, independent of the production
 *          expression, so it pins the sign of the log and the hue_shift factor
 *          that a ratio-only invariant cannot see. Rates avoid half-turn
 *          multiples, where a sign flip is unobservable.
 */
inline void test_sync_hue_rotates_per_efold() {
  constexpr float HS_TURN_TOL = 2e-3f;
  const float efolds[] = {0.5f, 1.0f, 2.0f, 3.0f};
  const float shifts[] = {0.05f, 0.2f, 0.33f};
  for (float efold : efolds)
    for (float shift : shifts) {
      Feedback::Style s{};
      s.fade = std::exp(-efold);
      s.hue_shift = shift;
      s.sync_hue();
      const float radians = 2.0f * PI_F * efold * shift;
      HS_EXPECT_NEAR(s.hue_ca, std::cos(radians), HS_TURN_TOL);
      HS_EXPECT_NEAR(s.hue_sa, std::sin(radians), HS_TURN_TOL);
    }
}

// --- lerp -------------------------------------------------------------------

/**
 * @brief Second ColorFn used only to observe lerp's function-pointer snapping.
 * @param p Source pixel color.
 * @param fade Per-frame scalar fade multiplier in [0, 1].
 * @return Pixel scaled by fade.
 */
inline Pixel lerp_probe_fade(const Pixel &p, float fade,
                             const Feedback::Style &) {
  return p * fade;
}

/**
 * @brief Verifies Style::lerp interpolates scalar fields linearly, snaps
 *        discrete fields, and preserves the subject's bound noise state.
 * @details Scalar fields (fade, hue_shift, scale, ...) interpolate linearly;
 *          discrete fields (transform pointers, downsample) snap to b at
 *          t >= 0.5 and stay on a below the midpoint. The subject's bound
 *          noise pointer must never be overwritten by a's or b's noise.
 */
inline void test_lerp_scalars_and_snapping() {
  Animation::NoiseParams na;
  Feedback::Style a{};
  a.fade = 0.0f;
  a.hue_shift = 0.0f;
  a.amplitude = 0.0f;
  a.frequency = 0.0f;
  a.speed = 0.0f;
  a.scale = 0.0f;
  a.space_fn = &Feedback::melt_warp;
  a.color_fn = &lerp_probe_fade;
  a.downsample = 2;
  a.noise = &na;

  Feedback::Style b{};
  b.fade = 1.0f;
  b.hue_shift = 1.0f;
  // amplitude/frequency/speed carry endpoints whose midpoint differs from the
  // Style default, so dropping their interpolation cannot land on the value the
  // assertion expects.
  b.amplitude = 3.0f;
  b.frequency = 0.7f;
  b.speed = 5.0f;
  b.scale = 1.0f;
  b.space_fn = &Feedback::noise_warp;
  b.color_fn = &Feedback::hue_fade;
  b.downsample = 8;
  b.noise = nullptr;

  Animation::NoiseParams subj;
  Feedback::Style mid{};
  mid.noise = &subj;
  mid.lerp(a, b, 0.5f);
  HS_EXPECT_NEAR(mid.fade, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(mid.hue_shift, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(mid.amplitude, 1.5f, 1e-6f);
  HS_EXPECT_NEAR(mid.frequency, 0.35f, 1e-6f);
  HS_EXPECT_NEAR(mid.speed, 2.5f, 1e-6f);
  HS_EXPECT_NEAR(mid.scale, 0.5f, 1e-6f);
  HS_EXPECT_TRUE(mid.space_fn == &Feedback::noise_warp);
  HS_EXPECT_TRUE(mid.color_fn == &Feedback::hue_fade);
  HS_EXPECT_EQ(mid.downsample, 8);
  HS_EXPECT_TRUE(mid.noise == &subj);

  Feedback::Style lo{};
  lo.lerp(a, b, 0.4f);
  HS_EXPECT_TRUE(lo.space_fn == &Feedback::melt_warp);
  HS_EXPECT_TRUE(lo.color_fn == &lerp_probe_fade);
  HS_EXPECT_EQ(lo.downsample, 2);
}

// --- Transform functions ----------------------------------------------------

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
  Animation::NoiseParams np;
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
  Animation::NoiseParams np;
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
 * @brief Verifies equal tail brightness produces equal accumulated hue.
 */
inline void test_sync_hue_matches_hue_at_equal_brightness() {
  struct TailPoint {
    float brightness;
    float ca;
    float sa;
  };
  const auto sample_tail = [](Feedback::Style style, int frames) {
    style.sync_hue();
    TailPoint point{1.0f, 1.0f, 0.0f};
    for (int i = 0; i < frames; ++i) {
      point.brightness *= style.fade;
      float ca = point.ca * style.hue_ca - point.sa * style.hue_sa;
      point.sa = point.sa * style.hue_ca + point.ca * style.hue_sa;
      point.ca = ca;
    }
    return point;
  };

  Feedback::Style long_tail{};
  long_tail.fade = 0.5f;
  long_tail.hue_shift = 0.01f;
  Feedback::Style short_tail = long_tail;
  short_tail.fade = 0.25f;

  for (int short_frames = 1; short_frames <= 3; ++short_frames) {
    TailPoint short_point = sample_tail(short_tail, short_frames);
    TailPoint long_point = sample_tail(long_tail, short_frames * 2);
    HS_EXPECT_NEAR(short_point.brightness, long_point.brightness, 1e-6f);
    HS_EXPECT_NEAR(short_point.ca, long_point.ca, 2e-3f);
    HS_EXPECT_NEAR(short_point.sa, long_point.sa, 2e-3f);
  }

  Feedback::Style no_history{};
  no_history.fade = 0.0f;
  no_history.hue_shift = 0.25f;
  no_history.sync_hue();
  HS_EXPECT_NEAR(no_history.hue_ca, 1.0f, 1e-6f);
  HS_EXPECT_NEAR(no_history.hue_sa, 0.0f, 1e-6f);
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

  Pixel plain = red * fade;
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
  const Pixel colors[] = {Pixel(65535, 0, 0),        Pixel(0, 65535, 0),
                          Pixel(0, 0, 65535),        Pixel(65535, 65535, 0),
                          Pixel(50000, 2000, 2000),  Pixel(300, 200, 100),
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
        Pixel ref =
            hue_rotate(Color4(c * fade, 1.0f), s.hue_ca, s.hue_sa).color;
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
 *          The 2880 pairs report through a worst-delta accumulator rather than
 *          per-pair assertions: per-pair HS_EXPECTs would put ~17k assertions
 *          into the module's floor and leave the floor gate blind to every
 *          other case in the module being deleted.
 */
inline void test_hue_fade_apply2_tracks_scalar() {
  constexpr int HS_PAIR_TOL = 128;
  const float channels[][3] = {{0.0f, 0.0f, 0.0f},
                               {1.0f, 0.0f, 0.0f},
                               {0.5f, 0.25f, 0.125f},
                               {300.0f, 200.0f, 100.0f},
                               {20000.0f, 20000.0f, 20000.0f},
                               {65535.0f, 0.0f, 0.0f},
                               {0.0f, 65535.0f, 0.0f},
                               {0.0f, 0.0f, 65535.0f},
                               {65535.0f, 65535.0f, 0.0f},
                               {50000.0f, 2000.0f, 2000.0f},
                               {65535.0f, 65535.0f, 65535.0f},
                               {12345.0f, 54321.0f, 999.0f}};
  const int n = (int)(sizeof(channels) / sizeof(channels[0]));
  const float fades[] = {0.58f, 0.75f, 0.9f, 0.99f};
  const float shifts[] = {0.0f, 0.01f, 0.1f, 0.33f, 0.75f};

  int worst_delta = 0;
  auto track = [&worst_delta](const Pixel &got, const Pixel &want) {
    worst_delta = std::max(worst_delta, std::abs((int)got.r - (int)want.r));
    worst_delta = std::max(worst_delta, std::abs((int)got.g - (int)want.g));
    worst_delta = std::max(worst_delta, std::abs((int)got.b - (int)want.b));
  };

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
          track(p0, s0);
          track(p1, s1);
        }
    }
  HS_EXPECT_LE(worst_delta, HS_PAIR_TOL);

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
  Animation::NoiseParams np;
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
  test_named_presets_preserve_frame_hue();
  test_sync_hue_rotates_per_efold();
  test_lerp_scalars_and_snapping();
  test_noise_warp_null_is_identity();
  test_noise_warp_bound_distorts();
  test_melt_warp_drifts_toward_north();
  test_melt_warp_bound_noise_perturbs();
  test_hue_fade_zero_shift_preserves_gray();
  test_sync_hue_matches_hue_at_equal_brightness();
  test_hue_fade_nonzero_shift_rotates_saturated();
  test_hue_rotate_lms_matrix_identity();
  test_hue_fade_matches_rotate_reference();
  test_hue_fade_apply2_tracks_scalar();
  test_sync_noise_pushes_scalars();

  return fixture.result();
}

} // namespace styles_tests
} // namespace hs_test
