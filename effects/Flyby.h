/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

// Unit-test accessor for the noise-time and trig-phase wrap invariants.
namespace hs_test {
namespace effects_tests {
struct FlybyWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Stereographic fly-through effect.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Projects a noise-warped Cartesian grid onto a sphere, rotating
 * around Y, blending continuously between camera/warp presets.
 * @note `Liquid2D` is the sibling stereographic effect, sharing the core
 *       primitives (`stereo()`, `stereo_noise_warp()`, `pole_attenuation()`,
 *       `pole_normalize_pattern()`) but with its own `project()`, `sample()`
 *       pattern, and `Params::lerp`; the two are independent effects, so a
 *       change here need not propagate there.
 */
template <int W, int H> class Flyby : public Effect {
public:
  /**
   * @brief Constructs the effect at W x H and disables pixel persistence.
   */
  FLASHMEM Flyby() : Effect(W, H, {.strobe = true}) {}

  /**
   * @brief Initializes animated params, orientation, palette, and preset loop.
   * @details Registers animated params, seeds orientation and rotation, bakes
   * the palette LUT, and kicks off the looping preset blend.
   */
  void init() override {

    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

    register_param("Warp Scale", &params.warp_scale, 0.1f, 100.0f);
    register_param("Warp Strength", &params.warp_strength, 0.0f, 30.0f);
    register_param("Pattern Freq", &params.pattern_freq, 1.0f, 20.0f);
    register_param("Speed", &params.speed, 0.0f, 2.0f);
    register_param("Pole Fade", &params.pole_fade, 1.0f, 20.0f);
    register_param("Drift", &drift, 0.0f, 2.0f);
    register_param("Hue Shift", &params.hue_shift, 0.0f, 1.0f);
    // Flag every preset-driven param so "Pause Animation" lets the user take a
    // slider over. Drift is a standalone live control, not preset-driven, so it
    // is omitted and edits apply during normal playback.
    for (const char *n :
         {"Warp Scale", "Warp Strength", "Pattern Freq", "Speed", "Pole Fade",
          "Hue Shift"})
      mark_animated(n);

    orientation.set(make_rotation(Vector(0, 0, -1), Vector(0, -1, 0)));
    timeline.add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 300,
                                           ease_linear, true));

    palette.bake(persistent_arena,
                 GenerativePalette{
                     GradientShape::STRAIGHT, HarmonyType::SPLIT_COMPLEMENTARY,
                     BrightnessProfile::FLAT, SaturationProfile::MID, 42});

    next_preset();
  }

  /**
   * @brief Advances to the next preset and schedules a blend into it.
   * @details Schedules a LERP_FRAMES blend into the next preset, re-arming
   * itself on completion to loop forever.
   */
  void next_preset() {
    constexpr int LERP_FRAMES = 480;
    presets.next();
    timeline.add(0, Animation::Lerp(params, presets.prev_get(), presets.get(),
                                    LERP_FRAMES, ease_in_out_sin, &anims_paused_)
                        .then([this]() { next_preset(); }));
  }

  /**
   * @brief Advances animation phases and shades every pixel for one frame.
   * @details Projects to stereographic space, warps by noise, samples the grid
   * pattern, attenuates near the pole, and colors via the palette with hue
   * rotated by the local warp displacement.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    // Wrap the noise-time accumulator so the float ULP never swallows the
    // increment and freezes the warp; OpenSimplex2 is aperiodic so the wrap pops
    // the field once per period (far apart at TIME_PERIOD).
    noise_time = fmodf(noise_time + params.speed, TIME_PERIOD);
    // Wrap the trig phases to 2pi so fast_sinf/fast_cosf keep precise range
    // reduction; each tracks its own coefficient (sin +t, cos -drift*t).
    constexpr float TWO_PI_F = 2.0f * PI_F;
    sin_phase = fmodf(sin_phase + params.speed, TWO_PI_F);
    drift_phase = fmodf(drift_phase + params.speed * drift, TWO_PI_F);

    auto shader = [&](const Vector &v) -> Color4 {
      Complex z = project(v);
      float r_sq = z.re * z.re + z.im * z.im;
      auto [w, displacement] = warp(z, r_sq, noise_time);
      float pattern = sample(w, sin_phase, drift_phase);
      float value = pole_normalize_pattern(pattern, r_sq, params.pole_fade);
      Color4 c = palette.get(value);
      c.alpha *= (1.0f - value);
      c = hue_rotate(c, -displacement * params.hue_shift);
      return c;
    };

    Scan::Shader::draw<W, H, 1>(canvas, shader);
  }

private:
  // Test seam: reaches the noise-time and trig-phase wrap invariants.
  friend struct ::hs_test::effects_tests::FlybyWhiteBox;

  /**
   * @brief Projects a sphere point to stereographic space.
   * @param v World-space unit vector on the sphere.
   * @return Stereographic image of @p v under the animated orientation.
   * @details Applies the animated orientation before stereographic projection.
   */
  Complex project(const Vector &v) const {
    return stereo(orientation.unorient(v));
  }

  /**
   * @brief Applies a noise-based warp in stereographic space.
   * @param z Stereographic coordinate to warp.
   * @param r_sq Squared radius z.re^2 + z.im^2, used for pole attenuation.
   * @param t Noise time axis; scaled by 0.3 before sampling.
   * @return Warped coordinate and warp displacement.
   * @details The warp is attenuated near the pole.
   */
  StereoWarpResult warp(const Complex &z, float r_sq, float t) const {
    return stereo_noise_warp(z, r_sq, noise, params.warp_scale,
                             params.warp_strength, params.pole_fade, t * 0.3f);
  }

  /**
   * @brief Samples the Cartesian grid pattern from warped coordinates.
   * @param w Warped stereographic coordinate.
   * @param sin_phase Wrapped +t term in [0, 2pi) (see draw_frame).
   * @param drift_phase Wrapped drift*t term in [0, 2pi) (see draw_frame).
   * @return Product of two sinusoids in [-1, 1] forming the grid pattern.
   */
  float sample(const Complex &w, float sin_phase, float drift_phase) const {
    // Near the pole |w| -> STEREO_INF, so w*pattern_freq can reach ~2e5 where
    // fast_sinf range reduction bands; clamp the (pole-attenuated) argument.
    float pu = hs::clamp(w.re * params.pattern_freq, -STEREO_PATTERN_ARG_LIMIT,
                         STEREO_PATTERN_ARG_LIMIT);
    float pv = hs::clamp(w.im * params.pattern_freq, -STEREO_PATTERN_ARG_LIMIT,
                         STEREO_PATTERN_ARG_LIMIT);
    return fast_sinf(pu + sin_phase) * fast_cosf(pv - drift_phase);
  }

  /**
   * @brief Wrap period for the noise-time accumulator (see draw_frame).
   * @details Large enough that warp pops are far apart, small enough that the
   *          float ULP never swallows the per-frame Speed increment.
   */
  static constexpr float TIME_PERIOD = 65536.0f;

  Timeline timeline;
  Orientation<> orientation;
  FastNoiseLite noise;
  float noise_time = 0.0f;   /**< Noise-time axis, wrapped to TIME_PERIOD (see draw_frame). */
  float sin_phase = 0.0f;    /**< Wrapped to [0, 2pi): the pattern's +t term. */
  float drift_phase = 0.0f;  /**< Wrapped to [0, 2pi): pattern's drift*t term. */
  float drift = 0.7f;        /**< Live drift-rate control; scales the cos phase advance. */

  BakedPalette palette;

  /**
   * @brief Tunable warp/pattern/color state, one snapshot per preset.
   */
  struct Params {
    float warp_scale = 1.5f;
    float warp_strength = 0.5f;
    float pattern_freq = 8.0f;
    float speed = 0.30f;
    float pole_fade = 2.0f;
    float hue_shift = 0.15f;

    /**
     * @brief Interpolates every field in parallel from a to b.
     * @param a Source params (t = 0).
     * @param b Destination params (t = 1).
     * @param t Blend factor in [0, 1].
     * @details Every field morphs simultaneously over t. A straight
     * simultaneous blend keeps the fly-through fluid, because each preset reads
     * as a single coherent camera/warp pose; sequencing the fields one at a
     * time would make the motion stutter field by field.
     */
    void lerp(const Params &a, const Params &b, float t) {
      constexpr int N = 6;
      // Trips if the field set changes, so the per-field lerp below and the
      // bare-float presets can't silently fall out of sync with Params.
      static_assert(sizeof(Params) == N * sizeof(float),
                    "Flyby::Params field set changed — update lerp and the "
                    "preset float lists to match");
      warp_scale = hs::lerp(a.warp_scale, b.warp_scale, t);
      warp_strength = hs::lerp(a.warp_strength, b.warp_strength, t);
      pattern_freq = hs::lerp(a.pattern_freq, b.pattern_freq, t);
      speed = hs::lerp(a.speed, b.speed, t);
      pole_fade = hs::lerp(a.pole_fade, b.pole_fade, t);
      hue_shift = hs::lerp(a.hue_shift, b.hue_shift, t);
    }
  };
  Params params;

  Presets<Params, 5> presets = {{{
      {{47.752f, 11.55f, 2.7f, 0.586f, 1.55f, 0.097f}},
      {{0.1f, 0.87f, 14.262f, 0.586f, 3.527f, 0.097f}},
      {{1.5f, 0.5f, 8.0f, 0.30f, 2.0f, 0.15f}},
      {{47.752f, 2.55f, 7.878f, 0.562f, 2.843f, 0.0f}},
      {{100.0f, 8.67f, 1.0f, 0.586f, 3.432f, 0.636f}},
  }}};
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(Flyby)
