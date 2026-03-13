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
        filters(Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>()) {}

  void init() override {
    registerParam("Scale", &params.noise_scale, 0.1f, 10.0f);
    registerParam("Force", &params.force_scale, 0.001f, 0.05f);
    registerParam("Max Spd", &params.max_speed, 0.01f, 0.1f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Time Spd", &params.time_scale, 0.001f, 0.05f);

    noise_generator.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise_generator.SetSeed(hs::rand_int(0, 65535));

    particle_system.init(persistent_arena, 0.96f, 0.0f, 300.0f);

    particle_system.add_emitter([this](ParticleSystem &sys) mutable {
      while (sys.active_count < sys.pool.capacity()) {
        sys.spawn(random_vector(), Vector(0, 0, 0), 0);
      }

      for (size_t i = 0; i < sys.active_count; ++i) {
        auto &p = sys.pool[i];

        if (hs::rand_f() < 0.005f) {
          p.position = random_vector();
          p.velocity = Vector(0, 0, 0);
          p.history.clear();
        }

        float fx =
            noise_generator.GetNoise(p.position.x * params.noise_scale,
                                     p.position.y * params.noise_scale,
                                     p.position.z * params.noise_scale + t) *
            params.force_scale;
        float fy =
            noise_generator.GetNoise(p.position.x * params.noise_scale + 100.0f,
                                     p.position.y * params.noise_scale,
                                     p.position.z * params.noise_scale + t) *
            params.force_scale;
        float fz =
            noise_generator.GetNoise(p.position.x * params.noise_scale + 200.0f,
                                     p.position.y * params.noise_scale,
                                     p.position.z * params.noise_scale + t) *
            params.force_scale;
        Vector force(fx, fy, fz);

        p.velocity = p.velocity + force;

        float speed = p.velocity.magnitude();
        if (speed > params.max_speed) {
          p.velocity = (p.velocity / speed) * params.max_speed;
        }
      }
    });
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    t += params.time_scale;

    particle_system.step(canvas);

    auto vertex_shader = [&](Fragment &f) {
      f.pos = orientation.orient(f.pos);
    };

    auto fragment_shader = [&](const Vector &v, Fragment &f) {
      float alpha = f.v0;
      float palette_t = (v.y + 1.0f) / 2.0f;
      Color4 c = palette.get(palette_t);
      c.alpha = c.alpha * alpha * params.alpha;
      f.color = c;
    };

    Plot::ParticleSystem::draw<W, H>(filters, canvas, particle_system,
                                     fragment_shader, vertex_shader);
  }

private:
  static constexpr int k_num_particles = 600;
  static constexpr int k_trail_length = 14;

  typedef Animation::ParticleSystem<W, k_num_particles, k_trail_length, 1, 0>
      ParticleSystem;

  struct Params {
    float noise_scale = 2.0f;
    float force_scale = 0.005f;
    float max_speed = 0.03f;
    float alpha = 0.8f;
    float time_scale = 0.005f;
  } params;

  float t = 0;
  FastNoiseLite noise_generator;
  GenerativePalette palette;
  ParticleSystem particle_system;
  Orientation<W> orientation;

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "../effect_registry.h"
REGISTER_EFFECT(FlowField)
