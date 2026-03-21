/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include "../scan.h"

template <int W, int H> class Flyby : public Effect {
public:
  FLASHMEM Flyby() : Effect(W, H) { persist_pixels = false; }

  void init() override {

    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

    registerParam("Warp Scale", &params.warp_scale, 0.1f, 100.0f);
    registerParam("Warp Strength", &params.warp_strength, 0.0f, 30.0f);
    registerParam("Pattern Freq", &params.pattern_freq, 1.0f, 20.0f);
    registerParam("Speed", &params.speed, 0.0f, 2.0f);
    registerParam("Pole Fade", &params.pole_fade, 1.0f, 20.0f);
    registerParam("Falloff", &params.falloff, 0.0f, 10.0f);
    registerParam("Drift", &params.drift, 0.0f, 2.0f);
    registerParam("Hue Shift", &params.hue_shift, 0.0f, 1.0f);

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
  }

  bool show_bg() const override { return true; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    float t = static_cast<float>(timeline.t) * params.speed;

    auto shader = [&](const Vector &v) -> Color4 {
      Complex z = project(v);
      float r_sq = z.re * z.re + z.im * z.im;
      auto [w, displacement] = warp(z, r_sq, t);
      float pattern = sample(w, t);
      float value = attenuate(pattern, r_sq);
      Color4 c = palette.get(value);
      c.alpha *= powf(1.0f - value, params.falloff);
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

  /// Cartesian grid pattern from warped coordinates.
  float sample(const Complex &w, float t) const {
    return sinf(w.re * params.pattern_freq + t) *
           cosf(w.im * params.pattern_freq - t * params.drift);
  }

  /// Pole attenuation applied to pattern, normalized to [0,1].
  float attenuate(float pattern, float r_sq) const {
    float fade = 1.0f / (1.0f + (r_sq / (params.pole_fade * params.pole_fade)));
    return (pattern * fade + 1.0f) * 0.5f;
  }

  Timeline<W> timeline;
  Orientation<W> orientation;
  FastNoiseLite noise;

  BakedPalette palette;

  struct Params {
    float warp_scale = 1.5f;
    float warp_strength = 0.5f;
    float pattern_freq = 8.0f;
    float speed = 0.30f;
    float pole_fade = 2.0f;
    float falloff = 0.0f;
    float drift = 0.7f;
    float hue_shift = 0.15f;

    void lerp(const Params &a, const Params &b, float t) {
      constexpr int N = 8;
      float *dst[N] = {&warp_scale, &warp_strength, &pattern_freq, &speed,
                       &pole_fade,  &falloff,       &drift,        &hue_shift};
      const float src[N] = {a.warp_scale, a.warp_strength, a.pattern_freq,
                            a.speed,      a.pole_fade,     a.falloff,
                            a.drift,      a.hue_shift};
      const float tgt[N] = {b.warp_scale, b.warp_strength, b.pattern_freq,
                            b.speed,      b.pole_fade,     b.falloff,
                            b.drift,      b.hue_shift};
      int active = 0;
      for (int i = 0; i < N; ++i)
        if (src[i] != tgt[i])
          ++active;
      if (active == 0)
        return;
      float slice = 1.0f / active;
      int slot = 0;
      for (int i = 0; i < N; ++i) {
        if (src[i] == tgt[i]) {
          *dst[i] = tgt[i];
          continue;
        }
        float tl = hs::clamp((t - slot * slice) / slice, 0.0f, 1.0f);
        *dst[i] = ::lerp(src[i], tgt[i], tl);
        ++slot;
      }
    }
  };
  Params params;

  Presets<Params, 4> presets = {{{
      {{47.752f, 11.55f, 1.0f, 0.586f, 2.843f, 1.2f, 0.0f, 0.097f}},
      {{1.5f, 0.5f, 8.0f, 0.30f, 2.0f, 1.2f, 0.7f, 0.15f}},
      {{47.752f, 2.55f, 7.878f, 0.562f, 2.843f, 1.2f, 0.0f, 0.0f}},
  }}};
};

#include "../effect_registry.h"
REGISTER_EFFECT(Flyby)
