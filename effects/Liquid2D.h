#pragma once

#include "../effects_engine.h"
#include <cmath>

template <int W, int H> class Liquid2D : public Effect {
public:
  Liquid2D()
      : Effect(W, H),
        myPalette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
                  BrightnessProfile::CUP) {
    persist_pixels = false;

    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

    registerParam("Warp Scale", &params.warp_scale, 0.1f, 10.0f);
    registerParam("Warp Strength", &params.warp_strength, 0.0f, 3.0f);
    registerParam("Pattern Freq", &params.pattern_freq, 1.0f, 20.0f);
    registerParam("Time Speed", &params.time_speed, 0.1f, 5.0f);
    registerParam("Complexity", &params.complexity, 0.5f, 3.0f);
    registerParam("Pole Fade", &params.pole_fade, 1.0f, 20.0f);

    timeline.add(0, Animation::RandomWalk<W>(orientation, UP));
  }

  // Emulates the JS/C++ memory interlacing glitch as a 2D Post-Process Lens
  void get_virtual_coords(float px, float py, float &vx, float &vy) {
    float h_half = static_cast<float>(H) / 2.0f;

    // 1. Mirror Southern Hemisphere to Northern Hemisphere
    float mirrored_py = py;
    if (py >= h_half) {
      mirrored_py = static_cast<float>(H - 1) - py;
    }

    // 2. The Horizontal Warp: In the glitch, C++ rows were half the width of JS
    // rows, causing the pattern to wrap around the physical sphere twice.
    vx = fmodf(px * 2.0f, static_cast<float>(W));

    // 3. The Vertical Squish: In the glitch, one JS row consumed two C++ rows.
    // This forces the math from North Pole to South Pole in half the space!
    vy = mirrored_py * 2.0f;

    // Note: The original glitch actually had a 1-pixel vertical seam because
    // the right half of the screen pulled from an odd C++ row while the left
    // pulled from an even. Leaving vy as `mirrored_py * 2.0f` fixes that seam
    // permanently!

    vy = hs::clamp(vy, 0.0f, static_cast<float>(H - 1));
  }

  virtual void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    float t = timeline.t * params.time_speed;

    float eps = 0.5f;
    float offsets_x[4] = {eps, -eps, eps, -eps};
    float offsets_y[4] = {eps, eps, -eps, -eps};

    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {

        // --- CALC THE NOISE ONCE PER LED USING THE VIRTUAL LENS ---
        float center_vx, center_vy;
        get_virtual_coords(static_cast<float>(x), static_cast<float>(y),
                           center_vx, center_vy);

        Vector center_v = pixel_to_vector<W, H>(center_vx, center_vy);
        Vector center_rot = orientation.orient(center_v);
        Complex center_z = stereo(center_rot);

        float noise_time = t * 0.5f;
        float warp_x =
            noise.GetNoise(center_z.re * params.warp_scale,
                           center_z.im * params.warp_scale, noise_time) *
            params.warp_strength;
        float warp_y = noise.GetNoise(center_z.re * params.warp_scale + 100.0f,
                                      center_z.im * params.warp_scale + 100.0f,
                                      noise_time) *
                       params.warp_strength;

        float total_pattern = 0.0f;
        float valid_samples = 0.0f;

        // --- 4x SSAA LOOP ---
        for (int i = 0; i < 4; ++i) {
          float px = static_cast<float>(x) + offsets_x[i];
          float py = static_cast<float>(y) + offsets_y[i];

          float vx, vy;
          get_virtual_coords(px, py, vx, vy);

          Vector sample_v = pixel_to_vector<W, H>(vx, vy);
          Vector sample_rot = orientation.orient(sample_v);
          Complex z = stereo(sample_rot);

          float u = z.re + warp_x;
          float v_coord = z.im + warp_y;

          float pu = u * params.pattern_freq;
          float pv = v_coord * params.pattern_freq;

          float pattern = sinf(pu + params.complexity * sinf(pv + t)) *
                          cosf(pv + params.complexity * cosf(pu - t * 0.8f));

          // Smoothly fade out both mathematical poles to prevent chaotic
          // strobing
          float r_sq = u * u + v_coord * v_coord;
          float attenuation =
              1.0f / (1.0f + (r_sq / (params.pole_fade * params.pole_fade)));

          total_pattern += pattern * attenuation;
          valid_samples += 1.0f;
        }

        if (valid_samples > 0.0f) {
          float avg_pattern = total_pattern / valid_samples;
          float normalized_pattern = (avg_pattern + 1.0f) * 0.5f;
          canvas(x, y) = myPalette.get(normalized_pattern).color;
        } else {
          canvas(x, y) = Pixel(0, 0, 0);
        }
      }
    }
  }

  virtual bool show_bg() const override { return true; }

private:
  Timeline<W> timeline;
  Orientation<W> orientation;
  FastNoiseLite noise;
  GenerativePalette myPalette;

  struct Params {
    float warp_scale = 1.5f;
    float warp_strength = 0.5f;
    float pattern_freq = 5.0f;
    float time_speed = 0.1f;
    float complexity = 1.8f;
    float pole_fade = 2.0f;

    Params lerp(const Params &p, float t) const {
      return {::lerp(warp_scale, p.warp_scale, t),
              ::lerp(warp_strength, p.warp_strength, t),
              ::lerp(pattern_freq, p.pattern_freq, t),
              ::lerp(time_speed, p.time_speed, t),
              ::lerp(complexity, p.complexity, t),
              ::lerp(pole_fade, p.pole_fade, t)};
    }
  };
  Params params;
};