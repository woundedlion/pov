/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

// Particle flow field over the sphere: particles advect along a 3D OpenSimplex2
// noise field (time-evolving via `t`), leaving SLERP-blurred trails. The whole
// field slowly drifts via a RandomWalk on `orientation`, applied by the
// Filter::World::Orient pipeline stage. W/H are the equirectangular canvas dims.
template <int W, int H> class FlowField : public Effect {
public:
  // Build the analogous-harmony palette and the Orient + AntiAlias filter
  // pipeline; noise seeding and emitter setup happen in init().
  FLASHMEM FlowField()
      : Effect(W, H), palette(GradientShape::STRAIGHT, HarmonyType::ANALOGOUS,
                              BrightnessProfile::ASCENDING),
        filters(Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>()) {}

  // Register UI params, seed the noise generators, start the orientation drift,
  // and install the emitter that keeps the pool full and advances particles
  // along the noise force field each step.
  void init() override {
    registerParam("Scale", &params.noise_scale, 0.1f, 10.0f);
    registerParam("Force", &params.force_scale, 0.001f, 0.05f);
    registerParam("Max Spd", &params.max_speed, 0.01f, 0.1f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Time Spd", &params.time_scale, 0.001f, 0.05f);

    noise_generator.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise_generator.SetSeed(hs::rand_int(0, 65535));

    // Drive the whole-sphere drift consumed by the Filter::World::Orient stage.
    // Without an animation the orientation stays at identity and the filter
    // (plus its SLERP motion-blur sweep) is a no-op. `orient_noise` is a
    // dedicated generator so RandomWalk never reconfigures the flow-field noise.
    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, random_vector(), orient_noise,
                        Animation::RandomWalk<W>::Options::Languid()));

    particle_system.init(persistent_arena, 0.96f, 0.0f, 300.0f);

    // Per-step emitter: keep the pool saturated, then push each particle by the
    // noise force field clamped to max_speed.
    particle_system.add_emitter([this](ParticleSystem &sys) mutable {
      while (sys.active_count < sys.pool.capacity()) {
        sys.spawn(random_vector(), Vector(0, 0, 0), 0);
      }

      for (size_t i = 0; i < sys.active_count; ++i) {
        auto &p = sys.pool[i];

        // Rare random respawn so particles caught in field sinks redistribute.
        if (hs::rand_f() < 0.005f) {
          p.position = random_vector();
          p.velocity = Vector(0, 0, 0);
          p.history.clear();
        }

        // Sample three decorrelated noise channels (offset along x by 100/200)
        // for the force vector's components.
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

  // Trails fully render the frame; no background fill is needed behind them.
  bool show_bg() const override { return false; }

  // Advance noise time, step the orientation drift and particles, then plot the
  // particle trails with a latitude-mapped palette through the orient/AA filters.
  void draw_frame() override {
    Canvas canvas(*this);
    // Wrap the noise-time accumulator so it never grows large enough for the
    // float ULP to swallow the increment and freeze the field. OpenSimplex2 is
    // not periodic, so the wrap produces a brief flow-field pop once per period;
    // at 2048 the ULP (~2.4e-4) stays well under the slowest Time Spd increment
    // (0.001) so t always advances, and the pop is hours apart at the default
    // rate (and power/show cycles reset it far sooner).
    t = fmodf(t + params.time_scale, TIME_PERIOD);

    timeline.step(canvas);
    particle_system.step(canvas);

    // Rotation is owned by the Filter::World::Orient stage in `filters`, which
    // also sweeps the orientation history for SLERP motion blur. Applying it
    // again in a vertex shader would rotate every point twice, so none is used.
    // Color by latitude (v.y in [-1,1] -> palette t in [0,1]); the trail's
    // per-vertex alpha (f.v0) and the global Alpha param fade the tail.
    auto fragment_shader = [&](const Vector &v, Fragment &f) {
      float alpha = f.v0;
      float palette_t = (v.y + 1.0f) / 2.0f;
      Color4 c = palette.get(palette_t);
      c.alpha = c.alpha * alpha * params.alpha;
      f.color = c;
    };

    Plot::ParticleSystem::draw<W, H>(filters, canvas, particle_system,
                                     fragment_shader);
  }

private:
  static constexpr int k_num_particles = 600;
  static constexpr int k_trail_length = 14;
  // Wrap period for the noise-time accumulator (see draw_frame): large enough
  // that pops are far apart, small enough that the ULP never swallows the
  // increment.
  static constexpr float TIME_PERIOD = 2048.0f;

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
  FastNoiseLite orient_noise;
  GenerativePalette palette;
  ParticleSystem particle_system;
  Orientation<> orientation;
  Timeline timeline;

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(FlowField)
