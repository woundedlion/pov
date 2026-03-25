/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

template <int W, int H> class MindSplatter : public Effect {
public:
  FLASHMEM MindSplatter()
      : Effect(W, H), presets{{{{{0.85f, 1.0f, 0.025f, 0.2f}},
                                {{0.85f, 1.0f, 0.025f, 0.52f}},
                                {{0.85f, 1.0f, 0.025f, 0.92f}},
                                {{0.85f, 2.0f, 0.094f, 0.2f}},
                                {{0.85f, 3.0f, 0.035f, 1.0f}}}}},
        filters(Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>()),
        particle_system() {}

  void init() override {
    configure_arenas(GLOBAL_ARENA_SIZE - 16384, 16384, 0);

    registerParam("Friction", &params.friction, 0.5f, 1.0f);
    registerParam("Well Str", &params.well_strength, 0.0f, 20.0f);
    registerParam("Init Spd", &params.initial_speed, 0.0f, 0.1f);
    registerParam("Ang Spd", &params.angular_speed, 0.0f, 1.0f);
    registerParam("Particles", &params.active_count, 0.0f,
                  (float)NUM_PARTICLES);
    registerParam("Run Presets", &run_presets, true);

    timeline.add(0, Animation::RandomWalk<W>(orientation, Y_AXIS, noise));

    auto preset_timer = Animation::PeriodicTimer(
        160,
        [this](Canvas &) {
          if (!run_presets)
            return;
          presets.next();
          timeline.add(0, Animation::Lerp(params, presets.prev_get(),
                                          presets.get(), 48, ease_mid));
        },
        true);
    timeline.add(0, preset_timer);

    rebuild();
    start_warp();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    particle_system.friction = params.friction;
    particle_system.step(canvas);
    params.active_count = (float)particle_system.active_count;

    draw_particles(canvas, 1.0f);
  }

private:
  static const int NUM_PARTICLES = 1024;

  typedef Solids::Cube EmitSolid;
  typedef Solids::Octahedron AttractSolid;

  typedef Animation::ParticleSystem<W, NUM_PARTICLES, 23, EmitSolid::NUM_VERTS,
                                    AttractSolid::NUM_VERTS>
      ParticleSystem;

  // Params
  struct Params {
    float friction = 0.85f;
    float well_strength = 1.0f;
    float initial_speed = 0.025f;
    float angular_speed = 0.2f;
    float active_count = 0.0f;

    void lerp(const Params &start, const Params &target, float t) {
      friction = start.friction + (target.friction - start.friction) * t;
      well_strength = start.well_strength +
                      (target.well_strength - start.well_strength) * t;
      initial_speed = start.initial_speed +
                      (target.initial_speed - start.initial_speed) * t;
      angular_speed = start.angular_speed +
                      (target.angular_speed - start.angular_speed) * t;
    }
  } params;

  Presets<Params, 5> presets;

  Orientation<W> orientation;
  FastNoiseLite noise;
  Timeline<W> timeline;
  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;
  ParticleSystem particle_system;
  GenerativePalette base_palette{
      GradientShape::CIRCULAR, HarmonyType::COMPLEMENTARY,
      BrightnessProfile::FLAT, SaturationProfile::MID};
  std::array<float, EmitSolid::NUM_VERTS> emitter_hues;
  std::array<int, EmitSolid::NUM_VERTS> emit_counters;

  // Warp params
  MobiusParams mobius;
  float warp_scale = 0.6f;
  bool run_presets = true;

  FLASHMEM void rebuild() {
    particle_system.init(persistent_arena, params.friction, 0.001f, 160.0f);

    // Add Attractors
    for (const auto &v : AttractSolid::vertices) {
      particle_system.add_attractor(v, params.well_strength, 0.003f, 0.2f);
    }

    // Emitter Hues
    for (size_t i = 0; i < EmitSolid::NUM_VERTS; ++i) {
      emitter_hues[i] = hs::rand_f();
    }

    emit_counters.fill(0);

    // Add Emitters
    for (size_t i = 0; i < EmitSolid::NUM_VERTS; ++i) {
      Vector axis = EmitSolid::vertices[i];

      particle_system.add_emitter([this, i, axis](ParticleSystem &sys) mutable {
        float angle = (emit_counters[i]++) * params.angular_speed;

        // Basis
        auto basis = make_basis(Quaternion(), axis);
        Vector vel = (basis.u * cosf(angle) + basis.w * sinf(angle)) *
                     params.initial_speed;

        // Update Hue
        emitter_hues[i] = fmodf(emitter_hues[i] + G * 0.1f, 1.0f);

        // Spawn with hue seed
        if (particle_system.active_count < particle_system.pool.capacity()) {
          uint16_t seed_u16 = static_cast<uint16_t>(emitter_hues[i] * 65535.0f);
          particle_system.spawn(axis, vel, seed_u16);
        }
      });
    }
  }

  void draw_particles(Canvas &canvas, float opacity = 1.0f) {
    auto vertex_shader = [&](Fragment &f) {
      Vector original_pos = f.pos;
      float holeAlpha = 1.0f;
      for (const auto &attr : particle_system.attractors) {
        float d = angle_between(original_pos, attr.position);
        if (d < attr.event_horizon) {
          float t = d / attr.event_horizon;
          holeAlpha *= quintic_kernel(t);
        }
      }

      f.pos = mobius_transform(f.pos, mobius);
      f.pos = orientation.orient(f.pos);
      f.v3 *= holeAlpha;
    };

    auto fragment_shader = [&](const Vector &v, Fragment &f) {
      float alpha = std::min(f.v0, f.v3);
      size_t p_idx = static_cast<size_t>(f.v2 + 0.5f);

      if (p_idx < 0 || p_idx >= particle_system.active_count) {
#ifdef DEBUG
        assert(false);
#endif
        f.color = Color4(CRGB(0, 0, 0), 0.0f);
        return;
      }

      const auto &p = particle_system.pool[p_idx];
      float seed_f = static_cast<float>(p.color_seed) / 65535.0f;
      float t_shifted = wrap_t(f.v0 + seed_f);
      Color4 c = hue_rotate(base_palette.get(t_shifted), seed_f);
      c.alpha = c.alpha * alpha * alpha * opacity;
      f.color = c;
    };

    Plot::ParticleSystem::draw<W, H>(filters, canvas, particle_system,
                                     fragment_shader, vertex_shader);
  }

  void start_warp() { schedule_warp(); }

  void schedule_warp() {
    auto timer =
        Animation::RandomTimer(180, 300, [this](Canvas &) { perform_warp(); });
    timeline.add(0, timer);
  }

  void perform_warp() {
    auto warp = Animation::MobiusWarp(mobius, warp_scale, 160, false);
    warp.then([this]() { schedule_warp(); });
    timeline.add(0, warp);
  }
};

#include "core/effect_registry.h"
REGISTER_EFFECT(MindSplatter)
