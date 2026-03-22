/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <cmath>
#include <algorithm>
#include "effects_engine.h"

template <int W, int H> class PetalFlow : public Effect {
public:
  struct Params {
    float twist_factor = 0.35f;
    float speed = 2.5f;
    float alpha = 0.2f;
  } params;

  FLASHMEM PetalFlow()
      : Effect(W, H),
        palette({0.029f, 0.029f, 0.029f}, {0.500f, 0.500f, 0.500f},
                {0.461f, 0.461f, 0.461f}, {0.539f, 0.701f, 0.809f}),
        orientation(), filters(Filter::World::Orient<W>(orientation),
                               Filter::Screen::AntiAlias<W, H>()),
        // Spawner: Manage ring creation based on gap accumulation
        spawner(1, [this](Canvas &) { this->check_spawn(); }, true) {}

  void init() override {
    // Register Params
    registerParam("Twist", &params.twist_factor, 0.0f, 5.0f);
    registerParam("Speed", &params.speed, 0.0f, 20.0f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);

    // Initialize Rings
    for (int i = 0; i < MAX_RINGS; ++i) {
      rings[i].active = false;
    }

    // Initialize Timeline
    init_timeline();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    update_and_draw_rings(canvas);
  }

private:
  static constexpr int MAX_RINGS = 64;
  static constexpr float START_RHO = -3.75f;
  static constexpr float END_RHO = 3.75f;
  static constexpr float SPACING = 0.3f; // Spacing in rho units

  struct Ring {
    float rho;
    float hue;
    bool active;
  };

  Ring rings[MAX_RINGS];
  float gap_accumulator = 0.0f;

  ProceduralPalette palette;
  Orientation<W> orientation;

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;

  Timeline<W> timeline;
  Animation::PeriodicTimer spawner;

  FLASHMEM void init_timeline() {
    timeline.add(0, Animation::Rotation<W>(orientation, UP, PI_F / 4.0f, 160,
                                           ease_mid, true));
    gap_accumulator = 0.0f;
    timeline.add(0, spawner);

    // Pre-warm rings along the path
    for (float r = END_RHO - 0.01f; r > START_RHO; r -= SPACING) {
      spawn_ring_at_pos(r);
    }
  }

  void check_spawn() {
    float move_dist = params.speed * 0.009375f;
    gap_accumulator += move_dist;

    while (gap_accumulator >= SPACING) {
      gap_accumulator -= SPACING;
      spawn_ring_at_pos(START_RHO + gap_accumulator);
    }
  }

  static float next_hue;

  void spawn_ring_at_pos(float initial_rho) {
    for (int i = 0; i < MAX_RINGS; ++i) {
      if (!rings[i].active) {
        rings[i].active = true;
        rings[i].rho = initial_rho;
        rings[i].hue = next_hue;
        next_hue = wrap(next_hue + 0.13f, 1.0f);
        return;
      }
    }
  }

  void update_and_draw_rings(Canvas &canvas) {
    float move_dist = params.speed * 0.009375f;
    for (int i = 0; i < MAX_RINGS; ++i) {
      if (!rings[i].active) continue;
      rings[i].rho += move_dist;
      if (rings[i].rho > END_RHO) {
        rings[i].active = false;
        continue;
      }
      draw_ring(canvas, rings[i]);
    }
  }

  void draw_ring(Canvas &canvas, const Ring &ring) {

    // Calculate Opacity (Distance Based)
    float dist = std::abs(ring.rho);
    float opacity = 1.0f;
    if (dist > 2.5f) {
      opacity = std::max(0.0f, 1.0f - (dist - 2.5f) / 1.0f);
    }

    // Combined Opacity
    float effective_opacity = opacity * params.alpha;
    if (effective_opacity <= 0.01f)
      return;

    const int num_samples = W / 2;
    const float step = 2.0f * PI_F / num_samples;

    Color4 base_col = palette.get(ring.hue).color;
    base_col.alpha = effective_opacity;

    float twist_angle = (ring.rho / SPACING) * params.twist_factor;

    auto get_shift = [](float t) -> float {
      return 0.6f * std::abs(sinf(3.0f * PI_F * t));
    };

    ScratchScope _(scratch_arena_a);
    Fragments fragments;
    fragments.bind(scratch_arena_a, num_samples + 1);
    for (int i = 0; i < num_samples; ++i) {
      float t_norm = static_cast<float>(i) / num_samples;
      float theta = i * step;

      float rho_val = ring.rho + get_shift(t_norm);
      float final_theta = theta + twist_angle;

      float R = expf(rho_val);

      float r2 = R * R;
      float denom = 1.0f + r2;
      float x = 2.0f * R * cosf(final_theta) / denom;
      float y = 2.0f * R * sinf(final_theta) / denom;
      float z = (r2 - 1.0f) / denom;

      Fragment f;
      f.pos = Vector(x, y, z);
      f.v0 = t_norm;
      f.age = 0;

      fragments.push_back(f);
    }

    if (!fragments.is_empty()) {
      Fragment f = fragments[0];
      f.v0 = 1.0f;
      fragments.push_back(f);
    }

    auto fragment_shader = [&](const Vector &, Fragment &f) {
      f.color = base_col;
    };

    Plot::rasterize<W, H>(filters, canvas, fragments, fragment_shader, true);
  }
};

template <int W, int H>
float PetalFlow<W, H>::next_hue = 0.0f;

#include "effect_registry.h"
REGISTER_EFFECT(PetalFlow)
