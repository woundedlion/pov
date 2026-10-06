/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for the Feedback::Style POD.
 */
#pragma once

#include <iterator>

#include "core/render/filter/feedback_style.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace styles_tests {

// --- Named presets ----------------------------------------------------------

/**
 * @brief Verifies preset parameter domains and their space/color transforms.
 */
inline void test_named_presets() {
  struct PresetCase {
    Feedback::Style style;
    Feedback::SpaceFn space;
  };
  const PresetCase presets[] = {
      {Feedback::Style::ArcingLightning(), &Feedback::noise_warp},
      {Feedback::Style::SlowFire(), &Feedback::noise_warp},
      {Feedback::Style::EnergeticFire(), &Feedback::noise_warp},
      {Feedback::Style::SlowDust(), &Feedback::noise_warp},
      {Feedback::Style::WavyTrails(), &Feedback::noise_warp},
      {Feedback::Style::MeltingHi(), &Feedback::melt_warp},
      {Feedback::Style::MeltingLo(), &Feedback::melt_warp},
      {Feedback::Style::Miasma(), &Feedback::noise_warp},
      {Feedback::Style::LooseWormhole(), &Feedback::noise_warp},
      {Feedback::Style::TightWormhole(), &Feedback::noise_warp},
      {Feedback::Style::WigglingWormhole(), &Feedback::noise_warp},
      {Feedback::Style::Smoke(), &Feedback::noise_warp},
  };
  for (size_t index = 0; index < std::size(presets); ++index) {
    const auto &style = presets[index].style;
    HS_CONTEXT("style", static_cast<long long>(index));
    HS_EXPECT_TRUE(std::isfinite(style.fade) && style.fade >= 0.0f &&
                   style.fade <= 1.0f);
    HS_EXPECT_TRUE(std::isfinite(style.hue_shift));
    HS_EXPECT_TRUE(std::isfinite(style.amplitude) && style.amplitude >= 0.0f);
    HS_EXPECT_TRUE(std::isfinite(style.frequency) && style.frequency >= 0.0f);
    HS_EXPECT_TRUE(std::isfinite(style.speed));
    HS_EXPECT_TRUE(std::isfinite(style.scale) && style.scale > 0.0f);
    HS_EXPECT_TRUE(style.color_fn == &Feedback::hue_fade);
    HS_EXPECT_TRUE(style.space_fn == presets[index].space);
  }
}

/**
 * @brief Verifies named presets retain their original per-frame hue rotations.
 * @details Recovers the angle from sync_hue()'s cached cos/sin.
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
    // The cached angle comes from fast trig.
    HS_EXPECT_NEAR(std::atan2(s.hue_sa, s.hue_ca) / (2.0f * math::PI_F),
                   c.frame_shift, 1e-3f);
  }
}

/**
 * @brief Pins sync_hue's per-frame rotation to hue_shift turns per e-fold of
 *        feedback brightness decay.
 * @details A fade of exp(-n) must rotate by n * hue_shift turns. Rates avoid
 *          half-turn multiples, where a sign flip is unobservable.
 */
inline void test_sync_hue_rotates_per_efold() {
  constexpr float TURN_TOL = 2e-3f;
  const float efolds[] = {0.5f, 1.0f, 2.0f, 3.0f};
  const float shifts[] = {0.05f, 0.2f, 0.33f};
  for (float efold : efolds)
    for (float shift : shifts) {
      Feedback::Style s{};
      s.fade = std::exp(-efold);
      s.hue_shift = shift;
      s.sync_hue();
      const float radians = 2.0f * math::PI_F * efold * shift;
      HS_EXPECT_NEAR(s.hue_ca, std::cos(radians), TURN_TOL);
      HS_EXPECT_NEAR(s.hue_sa, std::sin(radians), TURN_TOL);
    }
}

// --- lerp -------------------------------------------------------------------

/**
 * @brief Second ColorFn for observing lerp's function-pointer snapping.
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
 *        discrete fields, and pushes the blend into the subject's bound noise.
 * @details Discrete fields snap to b at t >= 0.5. The subject keeps its own
 *          noise pointer, and that NoiseParams carries the blended scalars.
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
  a.pole_half_res = 0;
  a.noise = &na;

  Feedback::Style b{};
  b.fade = 1.0f;
  b.hue_shift = 1.0f;
  // amplitude/frequency/speed midpoints differ from the Style default.
  b.amplitude = 3.0f;
  b.frequency = 0.7f;
  b.speed = 5.0f;
  b.scale = 1.0f;
  b.space_fn = &Feedback::noise_warp;
  b.color_fn = &Feedback::hue_fade;
  b.downsample = 8;
  b.pole_half_res = 1;
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
  HS_EXPECT_EQ(mid.pole_half_res, 1);
  HS_EXPECT_TRUE(mid.noise == &subj);
  HS_EXPECT_NEAR(subj.amplitude, 1.5f, 1e-6f);
  HS_EXPECT_NEAR(subj.frequency, 0.35f, 1e-6f);
  HS_EXPECT_NEAR(subj.speed, 2.5f, 1e-6f);
  HS_EXPECT_NEAR(subj.scale, 0.5f, 1e-6f);

  Feedback::Style lo{};
  lo.lerp(a, b, 0.4f);
  HS_EXPECT_TRUE(lo.space_fn == &Feedback::melt_warp);
  HS_EXPECT_TRUE(lo.color_fn == &lerp_probe_fade);
  HS_EXPECT_EQ(lo.downsample, 2);
  HS_EXPECT_EQ(lo.pole_half_res, 0);
}

// --- Transform functions ----------------------------------------------------

/**
 * @brief Verifies noise_warp passes the direction through unchanged when no
 *        NoiseParams is bound.
 */
inline void test_noise_warp_null_is_identity() {
  Feedback::Style s{};
  s.noise = nullptr;
  math::Vector v = math::Vector(0.6f, 0.4f, 0.69f).normalized();
  math::Vector out = Feedback::noise_warp(v, s);
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
  math::Vector v(1.0f, 0.0f, 0.0f); // on the equator (y = 0)
  math::Vector out = Feedback::melt_warp(v, s);
  // speed=1 slerps 0.04 of the 90 deg arc toward the pole: y rises ~0.0628, x
  // drops ~0.002.
  HS_EXPECT_TRUE(out.y > 0.05f);
  HS_EXPECT_TRUE(out.x < 0.999f);
  HS_EXPECT_NEAR(out.length(), 1.0f, 1e-4f);
}

/**
 * @brief Verifies noise_warp actually distorts when a NoiseParams is bound.
 * @details NoiseParams is primed via sync_noise(); the output must leave the
 *          input while staying unit length. Displacement is summed across
 *          samples so a single noise zero-crossing cannot pass as identity.
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

  const math::Vector samples[] = {math::Vector(1, 0, 0), math::Vector(0, 0, 1),
                                  math::Vector(0.4f, 0.6f, 0.7f).normalized()};
  float total_moved = 0.0f;
  for (const math::Vector &v : samples) {
    math::Vector out = Feedback::noise_warp(v, s);
    HS_EXPECT_NEAR(out.length(), 1.0f, 1e-3f);
    total_moved +=
        std::abs(out.x - v.x) + std::abs(out.y - v.y) + std::abs(out.z - v.z);
  }
  HS_EXPECT_GT(total_moved, 1e-2f);
}

/**
 * @brief Verifies melt_warp's bound-noise branch perturbs the drip.
 * @details With amplitude above the wobble floor, the output diverges from the
 *          same Style with noise unbound while staying unit length.
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

  const math::Vector samples[] = {math::Vector(1, 0, 0), math::Vector(0, 0, 1),
                                  math::Vector(0.4f, 0.6f, 0.7f).normalized()};
  float total_divergence = 0.0f;
  for (const math::Vector &v : samples) {
    math::Vector with_noise = Feedback::melt_warp(v, s);
    math::Vector pure_drip = Feedback::melt_warp(v, drip_only);
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
 * @details The rotated result must diverge from the hue-preserving fade.
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
 *          lms_to_oklab, so this pins the two OKLab matrices as mutual inverses.
 */
inline void test_hue_rotate_lms_matrix_identity() {
  float k[9];
  hue_rotate_lms_matrix(1.0f, 0.0f, k);
  for (int i = 0; i < 9; ++i)
    HS_EXPECT_NEAR(k[i], (i % 4 == 0) ? 1.0f : 0.0f, 1e-5f);
}

/**
 * @brief Parity sweep: hue_fade's folded cbrt-LMS path vs the reference
 *        fade-then-rotate composition through the tabulated gamut clip.
 * @details Allows 64 u16-channel LSBs. The reference rotates in OKLab with the
 * same cached trig pair and, out of gamut, rescales chroma onto the flash
 * grid's cell minimum.
 */
inline void test_hue_fade_matches_rotate_reference() {
  constexpr float HUE_FADE_TOL = 64.0f;
  const auto rotate_reference = [](const Pixel &faded, float ca, float sa) {
    const LinRGB in = pixel_to_linrgb(faded);
    const OKLab lab = linear_rgb_to_oklab(in.r, in.g, in.b);
    const OKLab rotated{lab.L, lab.a * ca - lab.b * sa,
                        lab.a * sa + lab.b * ca};
    LinRGB out = oklab_to_linear_rgb(rotated);
    if (!linear_rgb_in_gamut(out.r, out.g, out.b))
      out = oklab_to_linear_rgb(gamut_scale_to_boundary_lut(rotated));
    return linrgb_to_pixel(out);
  };
  alignas(uint16_t) static uint8_t
      lut_buf[gamut_lut_bytes(GAMUT_LUT_ANGLE_STEPS, GAMUT_LUT_L_STEPS)];
  Arena lut_arena(lut_buf, sizeof(lut_buf));
  init_gamut_lut(lut_arena, GAMUT_LUT_ANGLE_STEPS, GAMUT_LUT_L_STEPS);
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
        s.fade = fade;
        s.hue_shift = shift;
        s.sync_hue();
        Pixel got = Feedback::hue_fade(c, fade, s);
        Pixel ref = rotate_reference(c * fade, s.hue_ca, s.hue_sa);
        HS_EXPECT_NEAR((float)got.r, (float)ref.r, HUE_FADE_TOL);
        HS_EXPECT_NEAR((float)got.g, (float)ref.g, HUE_FADE_TOL);
        HS_EXPECT_NEAR((float)got.b, (float)ref.b, HUE_FADE_TOL);
      }
  release_gamut_lut();
}

/**
 * @brief Pins hue_fade_apply2 against the scalar hue_fade_apply at the
 *        quantized output, the level the display observes.
 * @details Allows 128 u16-channel LSBs.
 */
inline void test_hue_fade_apply2_tracks_scalar() {
  constexpr int PAIR_TOL = 128;
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
  const int n = (int)std::size(channels);
  const float fades[] = {0.58f, 0.75f, 0.9f, 0.99f};
  const float shifts[] = {0.0f, 0.01f, 0.1f, 0.33f, 0.75f};

  for (int si = 0; si < (int)std::size(shifts); ++si)
    for (int fi = 0; fi < (int)std::size(fades); ++fi) {
      HS_CONTEXT("shift/fade", si, fi);
      const float shift = shifts[si];
      const float fade = fades[fi];
      Feedback::Style s{};
      s.fade = fade;
      s.hue_shift = shift;
      s.sync_hue();
      // The composite folds the fade and the u16 normalization into k.
      float k[9];
      const float sc = math::fast_cbrt(fade * (1.0f / 65535.0f));
      for (int i = 0; i < 9; ++i)
        k[i] = s.hue_k[i] * sc;

      for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) {
          HS_CONTEXT("channel pair", a, b);
          Pixel s0 = Feedback::hue_fade_apply(k, channels[a][0], channels[a][1],
                                              channels[a][2]);
          Pixel s1 = Feedback::hue_fade_apply(k, channels[b][0], channels[b][1],
                                              channels[b][2]);
          Pixel p0, p1;
          Feedback::hue_fade_apply2(k, channels[a][0], channels[a][1],
                                    channels[a][2], channels[b][0],
                                    channels[b][1], channels[b][2], p0, p1);
          HS_EXPECT_LE(std::abs((int)p0.r - (int)s0.r), PAIR_TOL);
          HS_EXPECT_LE(std::abs((int)p0.g - (int)s0.g), PAIR_TOL);
          HS_EXPECT_LE(std::abs((int)p0.b - (int)s0.b), PAIR_TOL);
          HS_EXPECT_LE(std::abs((int)p1.r - (int)s1.r), PAIR_TOL);
          HS_EXPECT_LE(std::abs((int)p1.g - (int)s1.g), PAIR_TOL);
          HS_EXPECT_LE(std::abs((int)p1.b - (int)s1.b), PAIR_TOL);
        }
    }

  // An all-black pair stays exactly black.
  Feedback::Style s{};
  s.hue_shift = 0.2f;
  s.sync_hue();
  float k[9];
  const float sc = math::fast_cbrt(0.9f * (1.0f / 65535.0f));
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
  Animation::NoiseParams reference;
  reference.frequency = 0.33f;
  reference.sync();
  for (const auto &point :
       {math::Vector(1.25f, -2.75f, 4.5f), math::Vector(-7.0f, 0.125f, 3.0f)})
    HS_EXPECT_EQ(np.noise.GetNoise(point.x, point.y, point.z),
                 reference.noise.GetNoise(point.x, point.y, point.z));

  s.noise = nullptr;
  s.amplitude = 1.0f;
  s.sync_noise();
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
