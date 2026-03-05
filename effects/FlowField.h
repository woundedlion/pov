/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include <functional>
#include <memory>

#include <map>
#include "../effects_engine.h"
#include "../FastNoiseLite.h"

template <int W, int H> class FlowField : public Effect {
public:
  FLASHMEM FlowField()
      : Effect(W, H), palette(GradientShape::STRAIGHT, HarmonyType::ANALOGOUS,
                              BrightnessProfile::ASCENDING),
        filters(Filter::World::Trails<W, k_max_trail_dots>(k_trail_length),
                Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>()) {
    registerParam("Scale", &params.noise_scale, 0.1f, 10.0f);
    registerParam("Force", &params.force_scale, 0.001f, 0.05f);
    registerParam("Max Spd", &params.max_speed, 0.01f, 0.1f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Time Spd", &params.time_scale, 0.001f, 0.05f);

    noise_generator.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise_generator.SetSeed(hs::rand_int(0, 65535));
    reset_particles();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    t += params.time_scale;

    for (Particle &p : particles) {
      // Calculate Noise Force (Flow Field)
      // 4D noise: x, y, z, t
      float fx = noise_generator.GetNoise(p.pos.i * params.noise_scale,
                                          p.pos.j * params.noise_scale,
                                          p.pos.k * params.noise_scale + t) *
                 params.force_scale;
      float fy = noise_generator.GetNoise(p.pos.i * params.noise_scale + 100,
                                          p.pos.j * params.noise_scale,
                                          p.pos.k * params.noise_scale + t) *
                 params.force_scale;
      float fz = noise_generator.GetNoise(p.pos.i * params.noise_scale + 200,
                                          p.pos.j * params.noise_scale,
                                          p.pos.k * params.noise_scale + t) *
                 params.force_scale;
      Vector force(fx, fy, fz);

      // Update Velocity with Damping (Friction)
      p.vel = p.vel + force;
      p.vel = p.vel * 0.96f;

      float speed = p.vel.length();
      if (speed > params.max_speed) {
        p.vel = (p.vel / speed) * params.max_speed;
      } else if (speed < 0) {
        p.vel = Vector(0, 0, 0);
      }

      p.pos = p.pos + p.vel;
      p.pos.normalize();

      // Respawn Logic (Prevent sinks/clumping)
      if (hs::rand_f() < 0.005f) {
        p.pos = random_vector();
        p.vel = Vector(0, 0, 0);
      }

      // Map Y (-1 to 1) to (0 to 1) for palette
      float palette_t = (p.pos.j + 1.0f) / 2.0f;
      Color4 c = palette.get(palette_t);
      filters.plot(canvas, p.pos, c.color, 0, c.alpha);
    }

    // Render with Trails
    filters.flush(
        canvas,
        [this](const Vector &v, float t_trail) -> Color4 {
          return palette.get(t_trail);
        },
        params.alpha);
  }

private:
  struct Particle {
    Vector pos;
    Vector vel;

    Particle() : pos(random_vector()), vel(0, 0, 0) {}
  };

  static constexpr int k_num_particles = 600;
  static constexpr int k_trail_length = 14;

  struct Params {
    float noise_scale = 2.0f;
    float force_scale = 0.005f;
    float max_speed = 0.03f;
    float alpha = 0.8f;
    float time_scale = 0.005f;
  } params;

  static constexpr int k_max_trail_dots =
      k_num_particles * k_trail_length + k_num_particles; // Capacity buffer

  float t = 0;
  FastNoiseLite noise_generator;
  GenerativePalette palette;
  std::array<Particle, k_num_particles> particles;

  Orientation<W> orientation;

  Pipeline<W, H, Filter::World::Trails<W, k_max_trail_dots>,
           Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;

  void reset_particles() {
    for (int i = 0; i < k_num_particles; ++i) {
      particles[i] = Particle();
    }
  }
};