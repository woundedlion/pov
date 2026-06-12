/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

/// Stereographic fly-through: noise-warped Cartesian grid on a sphere,
/// rotating around Y, blending continuously between camera/warp presets.
template <int W, int H> class Flyby : public Effect {
public:
  FLASHMEM Flyby() : Effect(W, H) { persist_pixels = false; }

  /// Register animated params, seed orientation and rotation, bake the
  /// palette LUT, and kick off the looping preset blend.
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

  /// Advance to the next preset and schedule a LERP_FRAMES blend into it,
  /// re-arming itself on completion to loop forever.
  void next_preset() {
    constexpr int LERP_FRAMES = 480;
    presets.next();
    timeline.add(0, Animation::Lerp(params, presets.prev_get(), presets.get(),
                                    LERP_FRAMES, ease_in_out_sin, &anims_paused_)
                        .then([this]() { next_preset(); }));
  }

  /// Draw atop the persistent background buffer (this effect is translucent).
  bool show_bg() const override { return true; }

  /// Advance phases, then shade each pixel: project to stereographic space,
  /// warp by noise, sample the grid pattern, attenuate near the pole, and
  /// color via the palette with hue rotated by the local warp displacement.
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
  /// Stereographic projection with animated orientation.
  Complex project(const Vector &v) const {
    return stereo(orientation.unorient(v));
  }

  /// Noise-based warp in stereographic space, attenuated near pole.
  StereoWarpResult warp(const Complex &z, float r_sq, float t) const {
    return stereo_noise_warp(z, r_sq, noise, params.warp_scale,
                             params.warp_strength, params.pole_fade, t * 0.3f);
  }

  /// Cartesian grid pattern from warped coordinates. `sin_phase` is the wrapped
  /// `+t` term; `drift_phase` is the wrapped `drift*t` term (see draw_frame).
  float sample(const Complex &w, float sin_phase, float drift_phase) const {
    return fast_sinf(w.re * params.pattern_freq + sin_phase) *
           fast_cosf(w.im * params.pattern_freq - drift_phase);
  }

  /// Pole attenuation applied to pattern, normalized to [0,1].
  float attenuate(float pattern, float r_sq) const {
    float fade = pole_attenuation(r_sq, params.pole_fade);
    return (pattern * fade + 1.0f) * 0.5f;
  }

  Timeline timeline;
  Orientation<> orientation;
  FastNoiseLite noise;
  float noise_time = 0.0f;   // unbounded noise-time axis (see draw_frame)
  float sin_phase = 0.0f;    // wrapped to [0, 2pi): the pattern's +t term
  float drift_phase = 0.0f;  // wrapped to [0, 2pi): the pattern's drift*t term

  BakedPalette palette;

  /// Tunable warp/pattern/color state, one snapshot per preset.
  struct Params {
    float warp_scale = 1.5f;
    float warp_strength = 0.5f;
    float pattern_freq = 8.0f;
    float speed = 0.30f;
    float pole_fade = 2.0f;
    float drift = 0.7f;
    float hue_shift = 0.15f;

    /// Parallel interpolation: every field morphs simultaneously over t.
    /// Deliberately unlike Liquid2D's staggered, one-field-at-a-time slicing —
    /// Flyby's presets read as a single coherent camera/warp pose, so a
    /// straight simultaneous blend keeps the fly-through fluid, whereas
    /// sequencing the fields would make the motion stutter field by field.
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
