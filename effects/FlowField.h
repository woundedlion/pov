/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <functional>
#include <memory>
#include <vector>
#include <map>
#include "../effects_engine.h"
#include "../FastNoiseLite.h"

template <int W>
class FlowField : public Effect {
public:
  FlowField() :
    Effect(W),
    palette(GradientShape::STRAIGHT, HarmonyType::ANALOGOUS, BrightnessProfile::ASCENDING),
    trails(k_trail_length),
    filters(
      FilterOrient<W>(orientation),
      FilterAntiAlias<W>()
    )
  {
    persist_pixels = false;
    noise_generator.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise_generator.SetSeed(hs::rand_int(0, 65535));
    reset_particles();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    // timeline.step(canvas); // Timeline logic omitted for simple port
    t += k_time_scale;
    dots.clear();

    for (Particle& p : particles) {
      // 1. Calculate Noise Force (Flow Field)
      // 4D noise: x, y, z, t
      float fx = noise_generator.GetNoise(p.pos.i * k_noise_scale, p.pos.j * k_noise_scale, p.pos.k * k_noise_scale + t) * k_force_scale;
      float fy = noise_generator.GetNoise(p.pos.i * k_noise_scale + 100, p.pos.j * k_noise_scale, p.pos.k * k_noise_scale + t) * k_force_scale;
      float fz = noise_generator.GetNoise(p.pos.i * k_noise_scale + 200, p.pos.j * k_noise_scale, p.pos.k * k_noise_scale + t) * k_force_scale;
      Vector force(fx, fy, fz);

      // 2. Update Velocity with Damping (Friction)
      p.vel = p.vel + force;
      p.vel = p.vel * 0.96f; // Friction prevents runaway speed

      float speed = p.vel.length();
      if (speed > k_max_speed) {
        p.vel = (p.vel / speed) * k_max_speed;
      }
      else if (speed < 0) {
        p.vel = Vector(0, 0, 0);
      }

      // 3. Update Position
      p.pos = p.pos + p.vel;
      p.pos.normalize(); // Snap back to sphere

      // 4. Respawn Logic (Prevent sinks/clumping)
      if (hs::rand_f() < 0.005f) {
        p.pos = random_vector();
        p.vel = Vector(0, 0, 0);
      }

      // 5. Create Dot
      // Map Y (-1 to 1) to (0 to 1) for palette
      float palette_t = (p.pos.j + 1.0f) / 2.0f;
      dots.emplace_back(Dot(p.pos, palette.get(palette_t)));
    }

    // 6. Render with Trails
    trails.record(dots, 0, k_alpha);

    trails.render(canvas, filters,
      [this](const Vector& v, float t) {
        return palette.get(t);
      }
    );
  }

private:

  struct Particle {
    Vector pos;
    Vector vel;

    Particle() : pos(random_vector()), vel(0, 0, 0) {}
  };

  static constexpr int k_num_particles = 600;
  static constexpr int k_trail_length = 14;
  static constexpr float k_noise_scale = 2.0f;
  static constexpr float k_force_scale = 0.005f;
  static constexpr float k_max_speed = 0.03f;
  static constexpr float k_alpha = 0.8f;
  static constexpr float k_time_scale = 0.005f;

  static constexpr int k_max_trail_dots = k_num_particles * k_trail_length + k_num_particles; // Capacity buffer

  float t = 0;
  FastNoiseLite noise_generator;
  GenerativePalette palette;
  std::array<Particle, k_num_particles> particles;

  Orientation orientation;
  DecayBuffer<W, k_max_trail_dots> trails;

  Pipeline<W,
    FilterOrient<W>,
    FilterAntiAlias<W>
  > filters;
  Dots dots;

  void reset_particles() {
    for (int i = 0; i < k_num_particles; ++i) {
      particles[i] = Particle();
    }
  }
};