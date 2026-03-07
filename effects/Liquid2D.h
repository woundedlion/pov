#pragma once

#include "../effects_engine.h"
#include <cmath>

template <int W, int H> class Liquid2D : public Effect {
public:
  Liquid2D()
      : Effect(W, H),
        myPalette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
                  BrightnessProfile::CUP, SaturationProfile::VIBRANT, 141) {
    persist_pixels = false;

    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

    registerParam("Warp Scale", &params.warp_scale, 0.1f, 10.0f);
    registerParam("Warp Strength", &params.warp_strength, 0.0f, 3.0f);
    registerParam("Pattern Freq", &params.pattern_freq, 1.0f, 20.0f);
    registerParam("Time Speed", &params.time_speed, 0.1f, 5.0f);
    registerParam("Complexity", &params.complexity, 0.5f, 3.0f);
    registerParam("Pole Fade", &params.pole_fade, 1.0f, 20.0f);

    timeline.add(0, Animation::RandomWalk<W>(orientation, UP));
    timeline.add(0, Animation::RandomWalk<W>(global_orientation, UP));
  }

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

  Complex transform_point(const Vector &true_v) const {
    Vector rotated_v = global_orientation.unorient(true_v);
    Vector sample_v = apply_glitch_lens(rotated_v);
    return stereo(orientation.orient(sample_v));
  }

  float sample_pattern(const Complex &z, float warp_x, float warp_y,
                       float t) const {
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

  void calc_warp_noise(const Complex &z, float t, float &warp_x,
                       float &warp_y) const {
    float noise_time = t * 0.5f;
    warp_x = noise.GetNoise(z.re * params.warp_scale, z.im * params.warp_scale,
                            noise_time) *
             params.warp_strength;
    warp_y = noise.GetNoise(z.re * params.warp_scale + 100.0f,
                            z.im * params.warp_scale + 100.0f, noise_time) *
             params.warp_strength;
  }

  virtual void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    float t = timeline.t * params.time_speed;

    float eps = 0.5f;
    float offsets_x[4] = {eps, -eps, eps, -eps};
    float offsets_y[4] = {eps, eps, -eps, -eps};

    constexpr float h_virt_minus_1 = static_cast<float>(H + hs::H_OFFSET - 1);
    constexpr float w_float = static_cast<float>(W);

    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {

        // --- CALC THE NOISE ONCE PER LED ---
        Vector true_center_v = pixel_to_vector<W, H>(x, y);
        Complex center_z = transform_point(true_center_v);

        float warp_x, warp_y;
        calc_warp_noise(center_z, t, warp_x, warp_y);

        float total_pattern = 0.0f;

        // --- 4x SSAA LOOP ---
        for (int i = 0; i < 4; ++i) {
          float px = static_cast<float>(x) + offsets_x[i];
          float py = static_cast<float>(y) + offsets_y[i];

          // Direct calculation to bypass Spherical struct and clamping overhead
          float theta = (px * 2.0f * PI_F) / w_float;
          float phi = (py * PI_F) / h_virt_minus_1;
          float sin_phi = sinf(phi);
          Vector true_v(sin_phi * cosf(theta), cosf(phi),
                        sin_phi * sinf(theta));

          Complex z = transform_point(true_v);
          total_pattern += sample_pattern(z, warp_x, warp_y, t);
        }

        // We always take 4 samples, so * 0.25f is faster than division
        float avg_pattern = total_pattern * 0.25f;
        float normalized_pattern = (avg_pattern + 1.0f) * 0.5f;
        canvas(x, y) = myPalette.get(normalized_pattern).color;
      }
    }
  }

  // Preserved for physical POV inter-pixel blanking
  virtual bool show_bg() const override { return true; }

private:
  Timeline<W> timeline;
  Orientation<W> orientation;
  Orientation<W> global_orientation;
  FastNoiseLite noise;
  GenerativePalette myPalette;

  struct Params {
    float warp_scale = 1.5f;
    float warp_strength = 0.5f;
    float pattern_freq = 5.0f;
    float time_speed = 0.1f;
    float complexity = 0.5f;
    float pole_fade = 1.8f;
  };
  Params params;
};