#pragma once

#include "../ShaderEffect.h"

template <int W, int H> class Liquid2D : public ShaderEffect<W, H> {
  using Base = ShaderEffect<W, H>;

public:
  Liquid2D()
      : Base(GenerativePalette(
            GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
            BrightnessProfile::CUP, SaturationProfile::VIBRANT, 141)) {}

  void init() override {
    this->noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

    this->registerParam("Warp Scale", &params.warp_scale, 0.1f, 10.0f);
    this->registerParam("Warp Strength", &params.warp_strength, 0.0f, 3.0f);
    this->registerParam("Pattern Freq", &params.pattern_freq, 1.0f, 20.0f);
    this->registerParam("Time Speed", &params.time_speed, 0.1f, 5.0f);
    this->registerParam("Complexity", &params.complexity, 0.5f, 3.0f);
    this->registerParam("Pole Fade", &params.pole_fade, 1.0f, 20.0f);

    // Cycle presets every 3-5 seconds via a 2 second lerp
    this->timeline.add(0, Animation::RandomTimer(
                              90, 150,
                              [this](Canvas &) {
                                presets.next();
                                this->timeline.add(
                                    0, Animation::Lerp(
                                           params, presets.prev_get(),
                                           presets.get(), 60, ease_in_out_sin));
                              },
                              true));

    params = presets.get();
  }

  // --- Virtual hooks ---

  Complex transform(const Vector &v) const override {
    Vector rotated_v = this->global_orientation.unorient(v);
    Vector sample_v = apply_glitch_lens(rotated_v);
    return stereo(this->orientation.orient(sample_v));
  }

  float sample(const Complex &z, float warp_x, float warp_y,
               float t) const override {
    float u = z.re + warp_x;
    float v_coord = z.im + warp_y;

    float pu = u * params.pattern_freq;
    float pv = v_coord * params.pattern_freq;

    float pattern = sinf(pu + params.complexity * sinf(pv + t)) *
                    cosf(pv + params.complexity * cosf(pu - t * 0.8f));

    float r_sq = u * u + v_coord * v_coord;
    float attenuation =
        1.0f / (1.0f + (r_sq / (params.pole_fade * params.pole_fade)));

    return pattern * attenuation;
  }

  void calc_warp(const Complex &z, float t, float &warp_x,
                 float &warp_y) const override {
    float noise_time = t * 0.5f;
    warp_x = this->noise.GetNoise(z.re * params.warp_scale,
                                  z.im * params.warp_scale, noise_time) *
             params.warp_strength;
    warp_y =
        this->noise.GetNoise(z.re * params.warp_scale + 100.0f,
                             z.im * params.warp_scale + 100.0f, noise_time) *
        params.warp_strength;
  }

  float get_time() const override {
    return static_cast<float>(this->timeline.t) * params.time_speed;
  }

private:
  // Pure algebraic 3D Glitch Lens - Zero Trigonometry!
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
      warp_scale = ::lerp(a.warp_scale, b.warp_scale, t);
      warp_strength = ::lerp(a.warp_strength, b.warp_strength, t);
      pattern_freq = ::lerp(a.pattern_freq, b.pattern_freq, t);
      time_speed = ::lerp(a.time_speed, b.time_speed, t);
      complexity = ::lerp(a.complexity, b.complexity, t);
      pole_fade = ::lerp(a.pole_fade, b.pole_fade, t);
    }
  };
  Params params;

  Presets<Params, 4> presets = {{{
      {"Geometric", {1.5f, 0.5f, 5.0f, 0.1f, 0.5f, 1.8f}},
      {"SlowDrip", {1.5f, 0.5f, 1.2f, 0.05f, 3.0f, 2.0f}},
  }}};
};

#include "../effect_registry.h"
REGISTER_EFFECT(Liquid2D)
