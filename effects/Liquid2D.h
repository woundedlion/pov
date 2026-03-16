/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include "../scan.h"

template <int W, int H> class Liquid2D : public Effect {
public:
  FLASHMEM Liquid2D() : Effect(W, H) { persist_pixels = false; }

  void init() override {
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

    registerParam("Warp Scale", &params.warp_scale, 0.1f, 10.0f);
    registerParam("Warp Strength", &params.warp_strength, 0.0f, 3.0f);
    registerParam("Pattern Freq", &params.pattern_freq, 1.0f, 20.0f);
    registerParam("Time Speed", &params.time_speed, 0.1f, 5.0f);
    registerParam("Complexity", &params.complexity, 0.5f, 3.0f);
    registerParam("Pole Fade", &params.pole_fade, 1.0f, 20.0f);

    timeline.add(0, Animation::RandomWalk<W>(orientation, UP, noise));
    timeline.add(0, Animation::RandomWalk<W>(global_orientation, UP, noise));
    time_driver_ = timeline.add_get(
        0, Animation::Driver(accumulated_time, params.time_speed, false));

    // Cycle presets every 3-5 seconds via a 2 second lerp
    timeline.add(0, Animation::RandomTimer(
                        90, 150,
                        [this](Canvas &) {
                          presets.next();
                          timeline.add(0, Animation::Lerp(params,
                                                          presets.prev_get(),
                                                          presets.get(), 60,
                                                          ease_in_out_sin));
                        },
                        true));

    params = presets.get();
  }

  bool show_bg() const override { return true; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    time_driver_->set_speed(params.time_speed);
    float t = accumulated_time;
    float noise_time = t * 0.5f;

    struct WarpCtx {
      float warp_x, warp_y;
    };

    auto setup = [&](const Vector &v) -> WarpCtx {
      Vector rotated_v = global_orientation.unorient(v);
      Vector sample_v = apply_glitch_lens(rotated_v);
      Complex z = stereo(orientation.orient(sample_v));
      return {noise.GetNoise(z.re * params.warp_scale,
                             z.im * params.warp_scale, noise_time) *
                  params.warp_strength,
              noise.GetNoise(z.re * params.warp_scale + 100.0f,
                             z.im * params.warp_scale + 100.0f, noise_time) *
                  params.warp_strength};
    };

    auto sample = [&](const Vector &v, const WarpCtx &ctx) -> float {
      Vector rotated_v = global_orientation.unorient(v);
      Vector sample_v = apply_glitch_lens(rotated_v);
      Complex z = stereo(orientation.orient(sample_v));

      float u = z.re + ctx.warp_x;
      float v_coord = z.im + ctx.warp_y;
      float pu = u * params.pattern_freq;
      float pv = v_coord * params.pattern_freq;
      float pattern = sinf(pu + params.complexity * sinf(pv + t)) *
                      cosf(pv + params.complexity * cosf(pu - t * 0.8f));
      float r_sq = u * u + v_coord * v_coord;
      float attenuation =
          1.0f / (1.0f + (r_sq / (params.pole_fade * params.pole_fade)));
      return (pattern * attenuation + 1.0f) * 0.5f;
    };

    Scan::Shader::draw<W, H>(canvas, palette, setup, sample);
  }

private:
  Timeline<W> timeline;
  Orientation<W> orientation;
  Orientation<W> global_orientation;
  FastNoiseLite noise;
  GenerativePalette palette{
      GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
      BrightnessProfile::CUP, SaturationProfile::VIBRANT, 141};

  Vector apply_glitch_lens(Vector v) const {
    // 1. Mirror Southern Hemisphere
    if (v.y < 0.0f) {
      v.y = -v.y;
      v.z = -v.z; // X-axis reflection
    }

    float R2 = v.x * v.x + v.z * v.z;

    // North pole singularity protection
    if (R2 < 1e-6f) {
      return Vector(0.0f, 1.0f, 0.0f);
    }

    // 2. Trig-less Squish (phi * 2) and Warp (theta * 3)
    float inv_R2 = 1.0f / R2;

    return Vector(2.0f * v.y * v.x * (4.0f * v.x * v.x * inv_R2 - 3.0f), // i
                  2.0f * v.y * v.y - 1.0f,                               // j
                  2.0f * v.y * v.z * (3.0f - 4.0f * v.z * v.z * inv_R2)  // k
    );
  }

  struct Params {
    float warp_scale = 1.5f;
    float warp_strength = 0.5f;
    float pattern_freq = 5.0f;
    float time_speed = 0.1f;
    float complexity = 0.5f;
    float pole_fade = 1.8f;

    void lerp(const Params &a, const Params &b, float t) {
      constexpr int N = 6;
      float *dst[N] = {&warp_scale, &warp_strength, &pattern_freq,
                       &time_speed, &complexity,    &pole_fade};
      const float src[N] = {a.warp_scale, a.warp_strength, a.pattern_freq,
                            a.time_speed, a.complexity,    a.pole_fade};
      const float tgt[N] = {b.warp_scale, b.warp_strength, b.pattern_freq,
                            b.time_speed, b.complexity,    b.pole_fade};
      // Count parameters that actually change
      int active = 0;
      for (int i = 0; i < N; ++i)
        if (src[i] != tgt[i])
          ++active;
      if (active == 0)
        return;
      // Divide time equally among active params, lerp sequentially
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
  float accumulated_time = 0.0f;
  Animation::Driver *time_driver_ = nullptr;

  Presets<Params, 4> presets = {{{
      {"Geometric", {1.5f, 0.5f, 5.0f, 0.1f, 0.5f, 1.8f}},
      {"SlowDrip", {1.5f, 0.5f, 1.2f, 0.05f, 3.0f, 2.0f}},
  }}};
};

#include "../effect_registry.h"
REGISTER_EFFECT(Liquid2D)
