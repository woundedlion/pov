/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"

/**
 * @brief Stereographic fly-through effect.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Projects a noise-warped Cartesian grid onto a sphere, rotating
 * around Y, blending continuously between camera/warp presets.
 * @note Sibling stereographic shader — `Liquid2D` — shares this pipeline
 *       (stereo project → noise warp → sin/cos pattern → pole attenuate) and the
 *       per-pixel `sample()` (kept in sync, finding 399). They are not unified
 *       into a `StereoShaderBase`, so propagate shader fixes across both. Known
 *       divergences (finding 407): the two `Params::lerp` use different
 *       interpolation strategies and the warp time-scales differ by undocumented
 *       factors.
 */
template <int W, int H> class Flyby : public Effect {
public:
  /**
   * @brief Constructs the effect at W x H and disables pixel persistence.
   */
  FLASHMEM Flyby() : Effect(W, H) { persist_pixels = false; }

  /**
   * @brief Initializes animated params, orientation, palette, and preset loop.
   * @details Registers animated params, seeds orientation and rotation, bakes
   * the palette LUT, and kicks off the looping preset blend.
   */
  void init() override {

    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

    registerParam("Warp Scale", &params.warp_scale, 0.1f, 100.0f);
    registerParam("Warp Strength", &params.warp_strength, 0.0f, 30.0f);
    registerParam("Pattern Freq", &params.pattern_freq, 1.0f, 20.0f);
    registerParam("Speed", &params.speed, 0.0f, 2.0f);
    registerParam("Pole Fade", &params.pole_fade, 1.0f, 20.0f);
    registerParam("Drift", &params.drift, 0.0f, 2.0f);
    registerParam("Hue Shift", &params.hue_shift, 0.0f, 1.0f);
    // Every param is driven by the continuous preset lerp; flag them so the
    // standard "Pause Animation" toggle lets the user take a slider over.
    for (const char *n :
         {"Warp Scale", "Warp Strength", "Pattern Freq", "Speed", "Pole Fade",
          "Drift", "Hue Shift"})
      markAnimated(n);

    // Start with tangent plane at -Z, then rotate around Y
    orientation.set(make_rotation(Vector(0, 0, -1), Vector(0, -1, 0)));
    timeline.add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 300,
                                           ease_mid, true));

    // Bake the generative palette into a fast 16-bit LUT
    palette.bake(persistent_arena,
                 GenerativePalette{
                     GradientShape::STRAIGHT, HarmonyType::SPLIT_COMPLEMENTARY,
                     BrightnessProfile::FLAT, SaturationProfile::MID, 42});

    params = presets.get();
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
   * @brief Requests a trailing black display frame (a hardware POV-background
   * frame), not a software alpha-composite over a background buffer.
   * @return true; the driver emits a black frame after this effect's image so
   * the bright stereographic pattern stands against a dark POV background
   * instead of smearing over the previous revolution. The shader itself
   * OVERWRITES every pixel — `Scan::Shader::draw` writes `color * alpha`, and
   * `c.alpha *= (1 - value)` only darkens toward black; no persistent background
   * buffer is blended in or revealed.
   */
  bool show_bg() const override { return true; }

  /**
   * @brief Advances animation phases and shades every pixel for one frame.
   * @details Projects to stereographic space, warps by noise, samples the grid
   * pattern, attenuates near the pole, and colors via the palette with hue
   * rotated by the local warp displacement.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    // Noise time axis: grows unbounded. OpenSimplex2 is not periodic, so this
    // cannot be wrapped without a visible jump in the warp; it feeds GetNoise
    // (which degrades gracefully with large coordinates), not fast_sinf, so it
    // is not the range-reduction hazard the trig phases are.
    noise_time += params.speed;
    // Trig phases ARE wrapped to 2pi so fast_sinf/fast_cosf keep precise range
    // reduction over multi-hour runs. Each tracks its own time coefficient so
    // the wrap stays invisible (the pattern uses +t in sin and -drift*t in cos);
    // for the shipped presets drift == 0, so drift_phase stays 0.
    constexpr float TWO_PI = 2.0f * PI_F;
    sin_phase = fmodf(sin_phase + params.speed, TWO_PI);
    drift_phase = fmodf(drift_phase + params.speed * params.drift, TWO_PI);

    auto shader = [&](const Vector &v) -> Color4 {
      Complex z = project(v);
      float r_sq = z.re * z.re + z.im * z.im;
      auto [w, displacement] = warp(z, r_sq, noise_time);
      float pattern = sample(w, sin_phase, drift_phase);
      float value = attenuate(pattern, r_sq);
      Color4 c = palette.get(value);
      c.alpha *= (1.0f - value); // Linear alpha falloff toward bright values
      c = hue_rotate(c, -displacement * params.hue_shift);
      return c;
    };

    Scan::Shader::draw<W, H, 1>(canvas, shader);
  }

private:
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
    // Soft-limit the trig argument: near the pole |w| -> STEREO_INF, so
    // w*pattern_freq can reach ~2e5 where fast_sinf range reduction bands. The
    // pole cap is pole-attenuated anyway, so clamp rather than feed the trig a
    // coordinate it cannot resolve (see STEREO_PATTERN_ARG_LIMIT). Kept in sync
    // with Liquid2D::sample (finding 407).
    float pu = hs::clamp(w.re * params.pattern_freq, -STEREO_PATTERN_ARG_LIMIT,
                         STEREO_PATTERN_ARG_LIMIT);
    float pv = hs::clamp(w.im * params.pattern_freq, -STEREO_PATTERN_ARG_LIMIT,
                         STEREO_PATTERN_ARG_LIMIT);
    return fast_sinf(pu + sin_phase) * fast_cosf(pv - drift_phase);
  }

  /**
   * @brief Applies pole attenuation and normalizes the pattern to [0, 1].
   * @param pattern Raw grid pattern value in [-1, 1].
   * @param r_sq Squared stereographic radius driving pole attenuation.
   * @return Attenuated, normalized value in [0, 1].
   */
  float attenuate(float pattern, float r_sq) const {
    float fade = pole_attenuation(r_sq, params.pole_fade);
    return (pattern * fade + 1.0f) * 0.5f;
  }

  Timeline timeline;
  Orientation<> orientation;
  FastNoiseLite noise;
  float noise_time = 0.0f;   /**< Unbounded noise-time axis (see draw_frame). */
  float sin_phase = 0.0f;    /**< Wrapped to [0, 2pi): the pattern's +t term. */
  float drift_phase = 0.0f;  /**< Wrapped to [0, 2pi): pattern's drift*t term. */

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
    float drift = 0.7f;
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
      warp_scale = hs::lerp(a.warp_scale, b.warp_scale, t);
      warp_strength = hs::lerp(a.warp_strength, b.warp_strength, t);
      pattern_freq = hs::lerp(a.pattern_freq, b.pattern_freq, t);
      speed = hs::lerp(a.speed, b.speed, t);
      pole_fade = hs::lerp(a.pole_fade, b.pole_fade, t);
      drift = hs::lerp(a.drift, b.drift, t);
      hue_shift = hs::lerp(a.hue_shift, b.hue_shift, t);
    }
  };
  Params params;

  Presets<Params, 4> presets = {{{
      {{47.752f, 11.55f, 2.7f, 0.586f, 1.55f, 0.0f, 0.097f}},
      {{0.1f, 0.87f, 14.262f, 0.586f, 3.527f, 0.0f, 0.097f}},
      {{1.5f, 0.5f, 8.0f, 0.30f, 2.0f, 0.0f, 0.15f}},
      {{47.752f, 2.55f, 7.878f, 0.562f, 2.843f, 0.0f, 0.0f}},
  }}};
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Flyby)
