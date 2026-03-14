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

    registerParam("Warp Scale", &params.warp_scale, 0.1f, 10.0f);
    registerParam("Warp Strength", &params.warp_strength, 0.0f, 3.0f);
    registerParam("Pattern Freq", &params.pattern_freq, 1.0f, 20.0f);
    registerParam("Time Speed", &params.time_speed, 0.01f, 2.0f);
    registerParam("Complexity", &params.complexity, 0.5f, 3.0f);
    registerParam("Pole Fade", &params.pole_fade, 1.0f, 20.0f);

    timeline.add(0, Animation::RandomWalk<W>(orientation, UP, noise));
    timeline.add(0, Animation::RandomWalk<W>(global_orientation, UP, noise));

    params = presets.get();
  }

  bool show_bg() const override { return true; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    float t = static_cast<float>(timeline.t);

    auto shader = [&](const Vector &v) -> Color4 {
      // Transform
      Vector rotated_v = global_orientation.unorient(v);
      Complex z = stereo(orientation.orient(rotated_v));

      // Warp
      float noise_time = t * 0.5f;
      float warp_x =
          noise.GetNoise(z.re * params.warp_scale, z.im * params.warp_scale,
                         noise_time) *
          params.warp_strength;
      float warp_y =
          noise.GetNoise(z.re * params.warp_scale + 100.0f,
                         z.im * params.warp_scale + 100.0f, noise_time) *
          params.warp_strength;

      // Sample
      float u = z.re + warp_x;
      float v_coord = z.im + warp_y;
      float pattern = sinf(u * params.pattern_freq + t) *
                      cosf(v_coord * params.pattern_freq - t * 0.7f);
      float r_sq = u * u + v_coord * v_coord;
      float attenuation =
          1.0f / (1.0f + (r_sq / (params.pole_fade * params.pole_fade)));

      float normalized = (pattern * attenuation + 1.0f) * 0.5f;
      return palette.get(normalized);
    };

    Scan::Shader::draw<W, H>(canvas, shader);
  }

private:
  Timeline<W> timeline;
  Orientation<W> orientation;
  Orientation<W> global_orientation;
  FastNoiseLite noise;
  GenerativePalette palette{GradientShape::CIRCULAR, HarmonyType::ANALOGOUS,
                            BrightnessProfile::CUP, SaturationProfile::VIBRANT,
                            77};

  struct Params {
    float warp_scale = 1.5f;
    float warp_strength = 0.5f;
    float pattern_freq = 5.0f;
    float time_speed = 0.05f;
    float complexity = 1.0f;
    float pole_fade = 2.0f;

    void lerp(const Params &a, const Params &b, float t) {
      constexpr int N = 6;
      float *dst[N] = {&warp_scale, &warp_strength, &pattern_freq,
                       &time_speed, &complexity,    &pole_fade};
      const float src[N] = {a.warp_scale, a.warp_strength, a.pattern_freq,
                            a.time_speed, a.complexity,    a.pole_fade};
      const float tgt[N] = {b.warp_scale, b.warp_strength, b.pattern_freq,
                            b.time_speed, b.complexity,    b.pole_fade};
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
      {"Default", {1.5f, 0.5f, 5.0f, 0.05f, 1.0f, 2.0f}},
  }}};
};

#include "../effect_registry.h"
REGISTER_EFFECT(Flyby)
