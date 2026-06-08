/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include "core/effects_engine.h"

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
    registerParam("Num Balls", &params.num_balls, 1.0f,
                  static_cast<float>(BALL_CAPACITY));
    registerParam("Radius", &params.radius_scale, 0.1f, 3.0f);
    registerParam("Velocity", &params.velocity_scale, 0.0f, 3.0f);
    registerParam("Noise Str", &params.noise_strength, 0.0f, 0.05f);
    registerParam("Noise Spd", &params.noise_speed, 0.0f, 10.0f);
  }

  bool show_bg() const override { return false; } // We overwrite all pixels

  void draw_frame() override {
    Canvas canvas(*this);

    // Re-seed when the GUI changes the ball count or initial-velocity scale
    // (both are set at spawn time). Radius is applied live below, so it doesn't
    // require a re-seed. Integer count compare avoids thrashing on sub-unit drag.
    if (active_ball_count() != seeded_num_balls ||
        params.velocity_scale != seeded_velocity_scale)
      init_balls();

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

    // 2. Render via the shared full-screen shader path (clip-aware, instrumented).
    // SAMPLES=1: one sample at pixel center, identical to the former raw loop —
    // the metaball field is O(balls) per sample, too heavy to supersample.
    auto shader = [&](const Vector &v) -> Color4 {
      float sum = 0.0f;
      for (const auto &b : balls) {
        float dist_sq = distance_squared(v, b.p);
        // Avoid div by zero
        if (dist_sq < 0.0001f)
          dist_sq = 0.0001f;
        // Radius slider applied live (b.r is the base radius); a size tweak
        // grows the field without re-seeding ball positions.
        float br = b.r * params.radius_scale;
        sum += (br * br) / dist_sq;
      }

      float t_val = sum / params.max_influence;
      if (t_val > 1.0f)
        t_val = 1.0f;

      return palette.get(t_val);
    };
    Scan::Shader::draw<W, H, 1>(canvas, shader);
  }

  // Registered (live-tunable) params.
  struct Params {
    float max_influence = 10.0f;
    float gravity = 0.004f;
    float num_balls = 25.0f;      // GUI slider; re-seeds on integer change
    float radius_scale = 1.0f;    // GUI slider; applied live at render
    float velocity_scale = 0.7f;  // GUI slider; initial velocity, re-seeds
    float noise_strength = 0.0077f;
    float noise_speed = 4.0f;
  } params;

private:
  static constexpr int BALL_CAPACITY = 50;

  StaticCircularBuffer<Ball, BALL_CAPACITY> balls;
  ProceduralPalette palette; // initialized in the ctor (Palettes::richSunset)
  FastNoiseLite noise;
  float t = 0.0f;

  // Snapshot of the seed-time config, to detect GUI edits that need a re-seed.
  int seeded_num_balls = 0;
  float seeded_velocity_scale = 0.0f;

  int active_ball_count() const {
    return std::clamp(static_cast<int>(params.num_balls), 1, BALL_CAPACITY);
  }

  void init_balls() {
    const int n = active_ball_count();
    balls.clear();
    for (int i = 0; i < n; ++i) {
      Vector p = random_vector() * 0.5f; // Start inside
      // Base radius (unscaled); the Radius slider is applied live at render.
      float r = hs::rand_f(0.5f, 0.8f);
      Vector v = random_vector() * (0.05f * params.velocity_scale);
      balls.push_back({p, v, r});
    }
    seeded_num_balls = n;
    seeded_velocity_scale = params.velocity_scale;
  }
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Metaballs)
