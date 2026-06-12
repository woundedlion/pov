/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <cmath>
#include <algorithm>
#include "core/effects_engine.h"

// Concentric rings advance along a log-radial path and are projected onto the
// sphere via inverse stereographic projection. A per-ring radial wobble shapes
// each ring into petals, and a rho-dependent twist swirls the whole flow.
template <int W, int H> class PetalFlow : public Effect {
public:
  // User-tunable controls (registered as sliders in init()).
  //   twist_factor: swirl strength applied per unit of rho spacing.
  //   speed:        unitless flow rate; scaled by RHO_PER_SPEED into per-frame motion.
  //   alpha:        global opacity multiplier.
  struct Params {
    float twist_factor = 0.35f;
    float speed = 2.5f;
    float alpha = 0.2f;
  } params;

  // Wire up the palette, orientation filters, and the per-frame ring spawner.
  FLASHMEM PetalFlow()
      : Effect(W, H),
        palette({0.029f, 0.029f, 0.029f}, {0.500f, 0.500f, 0.500f},
                {0.461f, 0.461f, 0.461f}, {0.539f, 0.701f, 0.809f}),
        orientation(), filters(Filter::World::Orient<W>(orientation),
                               Filter::Screen::AntiAlias<W, H>()),
        // Fires every frame to spawn new rings as the flow accumulates gap.
        spawner(1, [this](Canvas &) { this->check_spawn(); }, true) {}

  // Register sliders, clear all rings, and seed the timeline. next_hue resets
  // here so hue assignment is deterministic on every (re)init.
  void init() override {
    registerParam("Twist", &params.twist_factor, 0.0f, 5.0f);
    registerParam("Speed", &params.speed, 0.0f, 20.0f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);

    for (int i = 0; i < MAX_RINGS; ++i) {
      rings[i].active = false;
    }
    next_hue = 0.0f;

    init_timeline();
  }

  bool show_bg() const override { return false; }

  // Step the timeline (rotation + spawner) then advance and render all rings.
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
  // Rho advanced per frame per unit of the Speed slider — the conversion from
  // the (unitless) Speed control to per-frame motion along the path.
  static constexpr float RHO_PER_SPEED = 0.009375f;
  // Petal shaping for the per-ring radial wobble: lobe count per revolution and
  // wobble depth in rho units.
  static constexpr float PETAL_LOBES = 3.0f;
  static constexpr float PETAL_DEPTH = 0.6f;

  // One ring on the flow. rho is its position along the log-radial path;
  // hue selects its palette color; inactive slots are free for reuse.
  struct Ring {
    float rho;
    float hue;
    bool active;
  };

  Ring rings[MAX_RINGS];
  float gap_accumulator = 0.0f;

  ProceduralPalette palette;
  Orientation<> orientation;

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;

  Timeline timeline;
  Animation::PeriodicTimer spawner;

  // Seed the timeline with the looping orientation rotation and the spawner,
  // then pre-fill the entire path with evenly spaced rings so the flow is full
  // from frame zero rather than filling in over time.
  FLASHMEM void init_timeline() {
    timeline.add(0, Animation::Rotation<W>(orientation, UP, PI_F / 4.0f, 160,
                                           ease_mid, true));
    gap_accumulator = 0.0f;
    timeline.add(0, spawner);

    for (float r = END_RHO - 0.01f; r > START_RHO; r -= SPACING) {
      spawn_ring_at_pos(r);
    }
  }

  // Accumulate this frame's travel and emit a ring each time a full SPACING of
  // gap opens up at the start of the path, keeping ring density constant.
  void check_spawn() {
    float move_dist = params.speed * RHO_PER_SPEED;
    gap_accumulator += move_dist;

    while (gap_accumulator >= SPACING) {
      gap_accumulator -= SPACING;
      spawn_ring_at_pos(START_RHO + gap_accumulator);
    }
  }

  // Per-instance hue cursor, reset in init() so hue assignment stays
  // deterministic for the fixed-seed segmented driver. Advanced per spawn.
  float next_hue = 0.0f;

  // Claim the first inactive ring slot, place it at initial_rho, assign the
  // next hue, and advance the hue cursor. No-op if all slots are in use.
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

  // Advance every active ring along rho; retire rings that pass END_RHO,
  // otherwise draw them.
  void update_and_draw_rings(Canvas &canvas) {
    float move_dist = params.speed * RHO_PER_SPEED;
    for (int i = 0; i < MAX_RINGS; ++i) {
      if (!rings[i].active)
        continue;
      rings[i].rho += move_dist;
      if (rings[i].rho > END_RHO) {
        rings[i].active = false;
        continue;
      }
      draw_ring(canvas, rings[i]);
    }
  }

  // Build and rasterize one ring: fade by distance from the equator (rho=0),
  // shape it into petals via a radial wobble, twist it by rho, and project each
  // sample onto the sphere. Skipped entirely when its opacity is negligible.
  void draw_ring(Canvas &canvas, const Ring &ring) {

    // Fade rings out past |rho|=2.5 over a window of 1.0 rho unit (near the poles).
    float dist = std::abs(ring.rho);
    float opacity = 1.0f;
    if (dist > 2.5f) {
      opacity = std::max(0.0f, 1.0f - (dist - 2.5f) / 1.0f);
    }

    float effective_opacity = opacity * params.alpha;
    if (effective_opacity <= 0.01f)
      return;

    const int num_samples = W / 2;
    const float step = 2.0f * PI_F / num_samples;

    Color4 base_col = palette.get(ring.hue).color;
    base_col.alpha = effective_opacity;

    float twist_angle = (ring.rho / SPACING) * params.twist_factor;

    // Radial wobble around the ring: rectified sine gives PETAL_LOBES lobes per
    // revolution, the petals. t is the normalized angle [0,1).
    auto get_shift = [](float t) -> float {
      return PETAL_DEPTH * std::abs(fast_sinf(PETAL_LOBES * PI_F * t));
    };

    ScratchScope _(scratch_arena_a);
    Fragments fragments;
    fragments.bind(scratch_arena_a, num_samples + 1);
    for (int i = 0; i < num_samples; ++i) {
      float t_norm = static_cast<float>(i) / num_samples;
      float theta = i * step;

      float rho_val = ring.rho + get_shift(t_norm);
      float final_theta = theta + twist_angle;

      // Inverse stereographic projection of the planar point at
      // (radius R=exp(rho_val), angle final_theta) onto the unit sphere.
      float R = expf(rho_val);

      float r2 = R * R;
      float denom = 1.0f + r2;
      float x = 2.0f * R * fast_cosf(final_theta) / denom;
      float y = 2.0f * R * fast_sinf(final_theta) / denom;
      float z = (r2 - 1.0f) / denom;

      Fragment f;
      f.pos = Vector(x, y, z);
      f.v0 = t_norm;
      f.age = 0;

      fragments.push_back(f);
    }

    // Close the loop: duplicate the first sample with v0=1.0 so the strip wraps
    // seamlessly back to the start.
    if (!fragments.is_empty()) {
      Fragment f = fragments[0];
      f.v0 = 1.0f;
      fragments.push_back(f);
    }

    // Flat-shade every fragment with the ring's distance-faded palette color.
    auto fragment_shader = [&](const Vector &, Fragment &f) {
      f.color = base_col;
    };

    Plot::rasterize<W, H>(filters, canvas, fragments, fragment_shader, true);
  }
};

#include "core/effect_registry.h"
REGISTER_EFFECT(PetalFlow)
