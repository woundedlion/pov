/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include "../effects_engine.h"

template <int W, int H> class Metaballs : public Effect {
public:
  struct Ball {
    Vector p;
    Vector v;
    float r;
  };

  FLASHMEM Metaballs() : Effect(W, H), palette(Palettes::richSunset) {}

  void init() override {
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise.SetSeed(hs::rand_int(0, 10000));
    init_balls();

    registerParam("Max Infl", &params.max_influence, 1.0f, 50.0f);
    registerParam("Gravity", &params.gravity, 0.0f, 0.02f);
    registerParam("Noise Str", &params.noise_strength, 0.0f, 0.05f);
    registerParam("Noise Spd", &params.noise_speed, 0.0f, 10.0f);
  }

  bool show_bg() const override { return false; } // We overwrite all pixels

  void draw_frame() override {
    Canvas canvas(*this);
    t += 0.01f;

    // 1. Physics
    for (size_t i = 0; i < balls.size(); i++) {
      auto &b = balls[i];

      // Gravity to center
      Vector center_force = b.p * (-params.gravity);

      // Noise force for wandering
      // Use ball index to offset noise so they don't all move the same way
      float nx = noise.GetNoise(t * params.noise_speed, (float)i * 10.0f, 0.0f);
      float ny =
          noise.GetNoise(t * params.noise_speed, (float)i * 10.0f, 100.0f);
      float nz =
          noise.GetNoise(t * params.noise_speed, (float)i * 10.0f, 200.0f);
      Vector noise_force(nx, ny, nz);
      noise_force = noise_force * params.noise_strength;

      b.v += center_force + noise_force;
      b.p += b.v;
    }

    // 2. Render
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        Vector v = pixel_to_vector<W, H>(x, y);

        float sum = 0.0f;
        for (const auto &b : balls) {
          float dist_sq = distance_squared(v, b.p);
          // Avoid div by zero
          if (dist_sq < 0.0001f)
            dist_sq = 0.0001f;
          sum += (b.r * b.r) / dist_sq;
        }

        float t_val = sum / params.max_influence;
        if (t_val > 1.0f)
          t_val = 1.0f;

        Color4 c = get_color(palette, t_val);
        canvas(x, y) = c.color;
      }
    }
  }

  // Params
  struct Params {
    float max_influence = 10.0f;
    float gravity = 0.004f;
    float num_balls = 25.0f;
    float radius_scale = 1.0f;
    float velocity_scale = 0.7f;
    float noise_strength = 0.0077f;
    float noise_speed = 4.0f;
  } params;

private:
  StaticCircularBuffer<Ball, 50> balls;
  ProceduralPalette palette = Palettes::richSunset;
  FastNoiseLite noise;
  float t = 0.0f;

  void init_balls() {
    balls.clear();
    for (int i = 0; i < (int)params.num_balls && i < 50; ++i) {
      Vector p = random_vector() * 0.5f; // Start inside
      float r = hs::rand_f(0.5f, 0.8f) * params.radius_scale;
      Vector v = random_vector() * (0.05f * params.velocity_scale);
      balls.push_back({p, v, r});
    }
  }

  float distance_squared(const Vector &a, const Vector &b) {
    float dx = a.x - b.x;
    float dy = a.y - b.y;
    float dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
  }
};

#include "../effect_registry.h"
REGISTER_EFFECT(Metaballs)
