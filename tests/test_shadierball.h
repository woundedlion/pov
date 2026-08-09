/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * ShadierBall white-box invariants: wrapped phase accumulators, the twin-wave
 * closed form, lens-slot dispatch through the projection, and the triadic
 * palette walk's morph-compatibility chain.
 */
#pragma once

#include <cstdint>

#include "tests/test_effects.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace shadierball_tests {

using effects_tests::FRAME_MS;
using effects_tests::FRAME_US;
using effects_tests::reset_effect_globals;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;

/**
 * @brief White-box accessor for ShadierBall's accumulators, pipeline stages,
 *        and palette provider (befriended in effects/ShadierBall.h).
 */
struct ShadierBallWhiteBox {
  using SDB = ShadierBall<SMALL_W, SMALL_H>;
  using Function = SDB::Function;
  using Projection = SDB::Projection;
  using Lens = SDB::Lens;
  using WaveState = SDB::WaveState;
  using Sample = SDB::Sample;
  static constexpr int FADE_FRAMES = SDB::PALETTE_FADE_FRAMES;
  static constexpr uint32_t HUE_STEP = SDB::HUE_STEP;
  static constexpr float AXIS_EPS = SDB::GNOMONIC_AXIS_EPS;
  static float phase(const SDB &sdb) { return sdb.phase; }
  static float wave_angle(const SDB &sdb) { return sdb.wave_angle; }
  static void seed_accumulators(SDB &sdb, float v) {
    sdb.phase = v;
    sdb.wave_angle = v;
  }
  static void set_walk(SDB &sdb, const Quaternion &q) { sdb.walk.set(q); }
  static void update_camera(SDB &sdb, float wander) {
    sdb.update_camera(wander);
  }
  static Quaternion camera(const SDB &sdb) { return sdb.camera; }
  static Complex project(const Vector &v, Projection projection, Lens lens,
                         float mix) {
    return SDB::project(v, projection, lens, mix);
  }
  static Complex project_point(const Vector &v, Projection projection) {
    return SDB::project_point(v, projection);
  }
  static Vector glitch_lens(const Vector &v) { return SDB::glitch_lens(v); }
  static Vector apply_lens(const Vector &v, Lens lens) {
    return SDB::apply_lens(v, lens);
  }
  static float sample_function(Function function, const Complex &p,
                               const WaveState &wave) {
    return SDB::sample_function(function, p, wave);
  }
  static Sample sample_pipeline(const Vector &v, Function function,
                                Projection projection, Lens lens,
                                const WaveState &wave, float pattern_freq,
                                float pole_fade, float lens_mix) {
    return SDB::sample_pipeline(v, {function, projection, lens}, wave,
                                pattern_freq, pole_fade, lens_mix);
  }
  static void make_palette(uint32_t &hue, uint32_t sequence,
                           GenerativePalette &out) {
    SDB::next_triadic_palette(&hue, sequence, out);
  }
  static Pixel palette_color(const SDB &sdb, float t) {
    return sdb.palette_cycler.palette().get(t).color;
  }
};

/**
 * @brief Verifies both phase accumulators stay in [0, 2pi) once wrapped.
 * @details draw_frame wraps phase and wave_angle by fmodf(., 2pi); a dropped
 *          wrap bands fast_sinf range reduction. Seeding both past the ceiling
 *          makes a removed fmodf fail on frame 0.
 */
inline void test_shadierball_phase_wrapped() {
  using WB = ShadierBallWhiteBox;
  reset_effect_globals();
  hs::set_mock_time(0, 0);
  WB::SDB sdb;
  sdb.init();

  const float two_pi = 2.0f * PI_F;
  WB::seed_accumulators(sdb, two_pi * 4.0f);

  for (int f = 0; f < 64; ++f) {
    hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                      static_cast<unsigned long>(f) * FRAME_US);
    sdb.draw_frame();
    sdb.advance_display();
    const float ph = WB::phase(sdb);
    HS_EXPECT_GE(ph, 0.0f);
    HS_EXPECT_LT(ph, two_pi);
    const float wa = WB::wave_angle(sdb);
    HS_EXPECT_GE(wa, 0.0f);
    HS_EXPECT_LT(wa, two_pi);
  }
  hs::clear_mock_time();
}

/**
 * @brief Pins the twin-wave function against its closed form.
 * @details sample_function at TWIN_WAVE must equal the mean of the fixed and
 *          rotated traveling waves bit for bit, and reduce to a single wave
 *          when the rotation angle is zero.
 */
inline void test_shadierball_twin_wave_formula() {
  using WB = ShadierBallWhiteBox;
  const float args[] = {-6.0f, -2.5f, -0.7f, 0.0f, 0.9f, 3.1f, 5.8f};
  const float phases[] = {0.0f, 0.8f, 2.4f, 4.9f};
  const float angles[] = {0.0f, 0.35f, 1.2f, 2.7f};
  for (float re : args) {
    for (float im : args) {
      const Complex p(re, im);
      for (float ph : phases) {
        for (float angle : angles) {
          const WB::WaveState wave{ph, angle, fast_cosf(angle),
                                   fast_sinf(angle)};
          const float rotated = re * wave.angle_cos + im * wave.angle_sin;
          const float expected =
              0.5f * (fast_sinf(re + ph) + fast_sinf(rotated + ph));
          HS_EXPECT_EQ(WB::sample_function(WB::Function::TWIN_WAVE, p, wave),
                       expected);
        }
      }
    }
  }

  // Zero rotation collapses the pair onto one wave exactly.
  const WB::WaveState aligned{1.1f, 0.0f, 1.0f, 0.0f};
  HS_EXPECT_EQ(WB::sample_function(WB::Function::TWIN_WAVE,
                                   Complex(0.7f, -1.1f), aligned),
               fast_sinf(0.7f + 1.1f));
}

/**
 * @brief Verifies lens-slot dispatch through the stereographic projection.
 * @details Every lens maps unit directions to unit directions; the twist
 *          preserves latitude and fixes the equator; the kaleidoscope folds
 *          all azimuths into one half-sector wedge and fixes the azimuth-free
 *          poles; Lens::NONE and the
 *          lens_mix endpoints reduce to the plain projections; a fractional
 *          mix stays linear through its midpoint; and the full pipeline lands
 *          every value in [0, 1].
 */
inline void test_shadierball_lens_projection() {
  using WB = ShadierBallWhiteBox;
  const Vector dirs[] = {Vector(1, 0, 0),
                         Vector(0, 1, 0),
                         Vector(0, -1, 0),
                         Vector(0, 0, 1),
                         Vector(-1, 0, 0),
                         Vector(0, 0, -1),
                         Vector(1, 1, 1).normalized(),
                         Vector(-1, 2, -3).normalized(),
                         Vector(3, -1, 2).normalized(),
                         Vector(0.2f, -0.9f, 0.4f).normalized(),
                         Vector(-0.7f, 0.1f, 0.7f).normalized()};
  for (const Vector &v : dirs) {
    HS_EXPECT_NEAR(WB::glitch_lens(v).length(), 1.0f, 1e-3f);
  }

  // The twist lens preserves latitude and turns the poles most; the
  // kaleidoscope folds every azimuth into one half-sector wedge. Both keep
  // unit norm.
  constexpr float SECTOR_HALF = PI_F / 6.0f;
  for (const Vector &d : dirs) {
    const Vector twisted = WB::apply_lens(d, WB::Lens::TWIST);
    // fast_sinf/fast_cosf are not exactly unit-amplitude, so the trig-built
    // lenses carry ~2e-3 of norm error by construction.
    HS_EXPECT_NEAR(twisted.length(), 1.0f, 4e-3f);
    HS_EXPECT_EQ(twisted.y, d.y);
    const Vector folded = WB::apply_lens(d, WB::Lens::KALEIDOSCOPE);
    HS_EXPECT_NEAR(folded.length(), 1.0f, 4e-3f);
    HS_EXPECT_EQ(folded.y, d.y);
    if (d.x == 0.0f && d.z == 0.0f) {
      // Poles have no azimuth to fold; the kaleidoscope fixes them exactly.
      HS_EXPECT_EQ(folded.x, 0.0f);
      HS_EXPECT_EQ(folded.z, 0.0f);
    } else {
      const float folded_azimuth = fast_atan2(folded.z, folded.x);
      // fast_atan2 and the fold's fast trig each carry ~4e-3 rad of error.
      HS_EXPECT_GE(folded_azimuth, -1e-2f);
      HS_EXPECT_LE(folded_azimuth, SECTOR_HALF + 1e-2f);
    }
  }
  const Vector equator(1, 0, 0);
  const Vector untwisted = WB::apply_lens(equator, WB::Lens::TWIST);
  HS_EXPECT_NEAR(untwisted.x, 1.0f, 1e-6f);
  HS_EXPECT_NEAR(untwisted.z, 0.0f, 1e-6f);

  const Vector v(0.6f, 0.48f, 0.64f);
  constexpr WB::Projection STEREO = WB::Projection::STEREOGRAPHIC;
  const Complex direct = stereo(v);
  const Complex lensed = stereo(WB::glitch_lens(v));
  const Complex none = WB::project(v, STEREO, WB::Lens::NONE, 0.7f);
  HS_EXPECT_EQ(none.re, direct.re);
  HS_EXPECT_EQ(none.im, direct.im);
  const Complex mix_zero = WB::project(v, STEREO, WB::Lens::GLITCH, 0.0f);
  HS_EXPECT_EQ(mix_zero.re, direct.re);
  HS_EXPECT_EQ(mix_zero.im, direct.im);
  const Complex mix_full = WB::project(v, STEREO, WB::Lens::GLITCH, 1.0f);
  HS_EXPECT_EQ(mix_full.re, lensed.re);
  HS_EXPECT_EQ(mix_full.im, lensed.im);

  const Complex before = WB::project(v, STEREO, WB::Lens::GLITCH, 0.499f);
  const Complex middle = WB::project(v, STEREO, WB::Lens::GLITCH, 0.5f);
  const Complex after = WB::project(v, STEREO, WB::Lens::GLITCH, 0.501f);
  HS_EXPECT_NEAR(middle.re - before.re, after.re - middle.re, 1e-5f);
  HS_EXPECT_NEAR(middle.im - before.im, after.im - middle.im, 1e-5f);

  const WB::WaveState wave{0.8f, 0.5f, fast_cosf(0.5f), fast_sinf(0.5f)};
  for (const Vector &d : dirs) {
    const float value =
        WB::sample_pipeline(d, WB::Function::TWIN_WAVE, STEREO,
                            WB::Lens::GLITCH, wave, 6.0f, 2.0f, 1.0f)
            .value;
    HS_EXPECT_GE(value, 0.0f);
    HS_EXPECT_LE(value, 1.0f);
  }
}

/**
 * @brief Pins the projection-slot dispatch and each map's closed form.
 * @details STEREOGRAPHIC must equal stereo() bit for bit; EQUIRECTANGULAR must
 *          equal (|azimuth| * cos(latitude), latitude), stay continuous across
 *          the longitude wrap, and collapse toward zero at the poles; GNOMONIC
 *          must equal the central division, map antipodal directions to the
 *          same point, and stay finite on the clamped equator ring; the full
 *          pipeline lands every value in [0, 1] under every projection.
 */
inline void test_shadierball_projection_maps() {
  using WB = ShadierBallWhiteBox;
  const Vector dirs[] = {Vector(1, 0, 0),
                         Vector(0, 1, 0),
                         Vector(0, -1, 0),
                         Vector(0, 0, 1),
                         Vector(1, 1, 1).normalized(),
                         Vector(-1, 2, -3).normalized(),
                         Vector(3, -1, 2).normalized(),
                         Vector(0.2f, -0.9f, 0.4f).normalized(),
                         Vector(-0.7f, 0.1f, 0.7f).normalized()};

  constexpr float HALF_PI = 0.5f * PI_F;
  for (const Vector &v : dirs) {
    const Complex st = WB::project_point(v, WB::Projection::STEREOGRAPHIC);
    const Complex st_expected = stereo(v);
    HS_EXPECT_EQ(st.re, st_expected.re);
    HS_EXPECT_EQ(st.im, st_expected.im);

    const Complex eq = WB::project_point(v, WB::Projection::EQUIRECTANGULAR);
    const float cyl_radius = sqrtf(v.x * v.x + v.z * v.z);
    HS_EXPECT_EQ(eq.re, std::fabs(fast_atan2(v.z, v.x)) * cyl_radius);
    HS_EXPECT_EQ(eq.im, HALF_PI - fast_acos(v.y));
    HS_EXPECT_GE(eq.re, 0.0f);
    HS_EXPECT_LE(eq.re, PI_F + 1e-2f);
    HS_EXPECT_GE(eq.im, -HALF_PI - 1e-2f);
    HS_EXPECT_LE(eq.im, HALF_PI + 1e-2f);

    const Complex gn = WB::project_point(v, WB::Projection::GNOMONIC);
    if (std::fabs(v.y) >= 0.1f) {
      HS_EXPECT_EQ(gn.re, v.x / v.y);
      HS_EXPECT_EQ(gn.im, v.z / v.y);
      const Complex gn_anti =
          WB::project_point(Vector(-v.x, -v.y, -v.z), WB::Projection::GNOMONIC);
      HS_EXPECT_EQ(gn.re, gn_anti.re);
      HS_EXPECT_EQ(gn.im, gn_anti.im);
    } else {
      // Equator ring: the |y| floor keeps the division finite.
      HS_EXPECT_LE(std::fabs(gn.re), 1.0f / WB::AXIS_EPS);
      HS_EXPECT_LE(std::fabs(gn.im), 1.0f / WB::AXIS_EPS);
    }
  }

  // The azimuth fold keeps the map continuous across the longitude wrap (the
  // z = 0, x < 0 half-plane) ...
  constexpr float WRAP_EPS = 1e-4f;
  const Complex wrap_above =
      WB::project_point(Vector(-0.8f, 0.6f, WRAP_EPS).normalized(),
                        WB::Projection::EQUIRECTANGULAR);
  const Complex wrap_below =
      WB::project_point(Vector(-0.8f, 0.6f, -WRAP_EPS).normalized(),
                        WB::Projection::EQUIRECTANGULAR);
  HS_EXPECT_NEAR(wrap_above.re, wrap_below.re, 1e-2f);
  HS_EXPECT_NEAR(wrap_above.im, wrap_below.im, 1e-6f);
  // ... and the cos(latitude) taper collapses the azimuth singularity at the
  // poles.
  const Complex near_pole =
      WB::project_point(Vector(0.03f, 0.999f, -0.03f).normalized(),
                        WB::Projection::EQUIRECTANGULAR);
  HS_EXPECT_LE(near_pole.re, PI_F * 0.05f);

  const WB::WaveState wave{0.8f, 0.5f, fast_cosf(0.5f), fast_sinf(0.5f)};
  const WB::Projection projections[] = {WB::Projection::EQUIRECTANGULAR,
                                        WB::Projection::STEREOGRAPHIC,
                                        WB::Projection::GNOMONIC};
  for (WB::Projection projection : projections) {
    for (const Vector &d : dirs) {
      const float value =
          WB::sample_pipeline(d, WB::Function::TWIN_WAVE, projection,
                              WB::Lens::GLITCH, wave, 6.0f, 2.0f, 1.0f)
              .value;
      HS_EXPECT_GE(value, 0.0f);
      HS_EXPECT_LE(value, 1.0f);
    }
  }
}

/**
 * @brief Pins the rings, spiral, and grid function bodies to their closed
 *        forms.
 * @details Each body must match its formula bit for bit; rings must ignore the
 *          rotation state; the spiral's integer arm count must keep the
 *          azimuth seam continuous; and the pipeline lands every value in
 *          [0, 1] for every function and lens.
 */
inline void test_shadierball_field_functions() {
  using WB = ShadierBallWhiteBox;
  const float args[] = {-6.0f, -2.5f, -0.7f, 0.0f, 0.9f, 3.1f, 5.8f};
  const WB::WaveState wave{0.9f, 0.7f, fast_cosf(0.7f), fast_sinf(0.7f)};
  const WB::WaveState spun{0.9f, 2.3f, fast_cosf(2.3f), fast_sinf(2.3f)};
  for (float re : args) {
    for (float im : args) {
      const Complex p(re, im);
      const float radius = sqrtf(re * re + im * im);

      HS_EXPECT_EQ(WB::sample_function(WB::Function::RINGS, p, wave),
                   fast_sinf(radius - wave.phase));
      HS_EXPECT_EQ(WB::sample_function(WB::Function::RINGS, p, wave),
                   WB::sample_function(WB::Function::RINGS, p, spun));

      const float azimuth = fast_atan2(im, re);
      HS_EXPECT_EQ(
          WB::sample_function(WB::Function::SPIRAL, p, wave),
          fast_sinf(radius - 3.0f * (azimuth + wave.angle) - wave.phase));

      const float a = re * wave.angle_cos + im * wave.angle_sin;
      const float b = -re * wave.angle_sin + im * wave.angle_cos;
      HS_EXPECT_EQ(WB::sample_function(WB::Function::GRID, p, wave),
                   fast_sinf(a + wave.phase) * fast_cosf(b - wave.phase));
    }
  }

  // Spiral azimuth seam: the +-pi wrap differs by 3 whole turns, so the value
  // is continuous across it.
  constexpr float EPSILON = 1e-4f;
  const float radius = 2.0f;
  const Complex above(radius * fast_cosf(PI_F - EPSILON),
                      radius * fast_sinf(PI_F - EPSILON));
  const Complex below(radius * fast_cosf(PI_F + EPSILON),
                      radius * fast_sinf(PI_F + EPSILON));
  HS_EXPECT_NEAR(WB::sample_function(WB::Function::SPIRAL, above, wave),
                 WB::sample_function(WB::Function::SPIRAL, below, wave), 5e-2f);

  const Vector dirs[] = {Vector(1, 0, 0), Vector(1, 1, 1).normalized(),
                         Vector(0.2f, -0.9f, 0.4f).normalized(),
                         Vector(-0.7f, 0.1f, 0.7f).normalized()};
  const WB::Function functions[] = {WB::Function::TWIN_WAVE,
                                    WB::Function::RINGS, WB::Function::SPIRAL,
                                    WB::Function::GRID};
  const WB::Lens lenses[] = {WB::Lens::NONE, WB::Lens::GLITCH, WB::Lens::TWIST,
                             WB::Lens::KALEIDOSCOPE};
  for (WB::Function function : functions) {
    for (WB::Lens lens : lenses) {
      for (const Vector &d : dirs) {
        const float value =
            WB::sample_pipeline(d, function, WB::Projection::STEREOGRAPHIC,
                                lens, wave, 6.0f, 2.0f, 1.0f)
                .value;
        HS_EXPECT_GE(value, 0.0f);
        HS_EXPECT_LE(value, 1.0f);
      }
    }
  }
}

/**
 * @brief Bounds the wander camera's motion by the walk step times the gain.
 * @details Integrating a partial gain must stay continuous through full turns:
 *          every per-frame camera step is at most the walk step scaled by the
 *          gain, a zero gain freezes the camera under a moving walk, and a
 *          stationary walk moves the camera at no gain.
 */
inline void test_shadierball_wander_camera() {
  using WB = ShadierBallWhiteBox;
  constexpr float WALK_STEP = 0.02f;
  constexpr float WANDER = 0.25f;
  constexpr float MAX_CAMERA_STEP = WALK_STEP * WANDER + 1e-4f;
  WB::SDB sdb;
  Quaternion prev = WB::camera(sdb);

  for (int frame = 0; frame < 400; ++frame) {
    const float angle = static_cast<float>(frame) * WALK_STEP;
    WB::set_walk(sdb, make_rotation(Y_AXIS, angle));
    WB::update_camera(sdb, WANDER);
    const Quaternion current = WB::camera(sdb);
    if (frame > 0) {
      const float step =
          2.0f * acosf(hs::clamp(std::fabs(dot(prev, current)), 0.0f, 1.0f));
      HS_EXPECT_LE(step, MAX_CAMERA_STEP);
    }
    prev = current;
  }

  // A stationary walk moves the camera at no gain.
  for (float wander : {0.0f, 0.25f, 0.7f, 1.0f}) {
    WB::update_camera(sdb, wander);
    HS_EXPECT_NEAR(std::fabs(dot(prev, WB::camera(sdb))), 1.0f, 1e-6f);
  }

  // A zero gain freezes the camera under a moving walk.
  const Quaternion frozen = WB::camera(sdb);
  for (int frame = 0; frame < 32; ++frame) {
    WB::set_walk(sdb, make_rotation(X_AXIS, 0.1f * frame));
    WB::update_camera(sdb, 0.0f);
    HS_EXPECT_NEAR(std::fabs(dot(frozen, WB::camera(sdb))), 1.0f, 1e-6f);
  }
}

/**
 * @brief Verifies the triadic palette walk morphs and advances.
 * @details Every provider successor must be morph-compatible with its
 *          predecessor (the cycler fail-fast checks each retarget), the hue
 *          index must advance by HUE_STEP per palette, a full in-effect fade
 *          must land on a different palette, and "Pause Animation" must never
 *          freeze the cycle.
 */
inline void test_shadierball_palette_walk() {
  using WB = ShadierBallWhiteBox;
  uint32_t hue = 0;
  GenerativePalette previous;
  WB::make_palette(hue, 0, previous);
  HS_EXPECT_EQ(hue, uint32_t(0));
  for (uint32_t sequence = 1; sequence <= 8; ++sequence) {
    GenerativePalette next;
    WB::make_palette(hue, sequence, next);
    HS_EXPECT_TRUE(previous.morph_compatible(next));
    previous = next;
  }
  HS_EXPECT_EQ(hue, uint32_t(8) * WB::HUE_STEP);

  reset_effect_globals();
  hs::set_mock_time(0, 0);
  WB::SDB sdb;
  sdb.init();
  const Pixel before = WB::palette_color(sdb, 0.25f);
  for (int f = 0; f < WB::FADE_FRAMES + 8; ++f) {
    hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                      static_cast<unsigned long>(f) * FRAME_US);
    sdb.draw_frame();
    sdb.advance_display();
  }
  const Pixel after = WB::palette_color(sdb, 0.25f);
  HS_EXPECT_TRUE(before.r != after.r || before.g != after.g ||
                 before.b != after.b);

  // The pause stops preset advancement only; the palette keeps cycling.
  sdb.setAnimationsPaused(true);
  for (int f = 0; f < WB::FADE_FRAMES + 8; ++f) {
    const unsigned long tick =
        static_cast<unsigned long>(f) + WB::FADE_FRAMES + 8;
    hs::set_mock_time(tick * FRAME_MS, tick * FRAME_US);
    sdb.draw_frame();
    sdb.advance_display();
  }
  const Pixel paused = WB::palette_color(sdb, 0.25f);
  HS_EXPECT_TRUE(after.r != paused.r || after.g != paused.g ||
                 after.b != paused.b);
  hs::clear_mock_time();
}

/**
 * @brief Module entry point for the ShadierBall white-box invariants.
 * @return Module result code from hs_test::end_module (0 on success).
 */
inline int run_shadierball_tests() {
  ModuleFixture fixture("shadierball");
  test_shadierball_phase_wrapped();
  test_shadierball_twin_wave_formula();
  test_shadierball_lens_projection();
  test_shadierball_projection_maps();
  test_shadierball_field_functions();
  test_shadierball_wander_camera();
  test_shadierball_palette_walk();
  return fixture.result();
}

} // namespace shadierball_tests
} // namespace hs_test
