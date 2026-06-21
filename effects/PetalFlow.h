/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <cmath>
#include <algorithm>
#include "core/engine.h"

/**
 * @brief Flowing petal effect built from concentric rings on the sphere.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Concentric rings advance along a log-radial path and are projected
 * onto the sphere via inverse stereographic projection. A per-ring radial
 * wobble shapes each ring into petals, and a rho-dependent twist swirls the
 * whole flow.
 */
template <int W, int H> class PetalFlow : public Effect {
public:
  /**
   * @brief User-tunable controls, registered as sliders in init().
   */
  struct Params {
    float twist_factor = 0.35f; /**< Swirl strength applied per unit of rho spacing. */
    float speed = 2.5f;         /**< Unitless flow rate; scaled by RHO_PER_SPEED into per-frame motion. */
    float alpha = 0.2f;         /**< Global opacity multiplier. */
  } params;

  /**
   * @brief Constructs the effect, wiring up palette, orientation filters, and
   * the per-frame ring spawner.
   * @details The spawner fires every frame to emit new rings as the flow
   * accumulates gap.
   */
  FLASHMEM PetalFlow()
      : Effect(W, H),
        palette({0.029f, 0.029f, 0.029f}, {0.500f, 0.500f, 0.500f},
                {0.461f, 0.461f, 0.461f}, {0.539f, 0.701f, 0.809f}),
        orientation(), filters(Filter::World::Orient<W>(orientation),
                               Filter::Screen::AntiAlias<W, H>()),
        // Fires every frame to spawn new rings as the flow accumulates gap.
        spawner(1, [this](Canvas &) { this->check_spawn(); }, true) {}

  /**
   * @brief Registers sliders, clears all rings, and seeds the timeline.
   * @details next_hue resets here so hue assignment is deterministic on every
   * (re)init.
   */
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

  /**
   * @brief Reports whether the engine should draw its background first.
   * @return false; this effect paints over a cleared frame.
   */
  bool show_bg() const override { return false; }

  /**
   * @brief Steps the timeline (rotation + spawner), then advances and renders
   * all rings.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    update_and_draw_rings(canvas);
  }

private:
  static constexpr int MAX_RINGS = 64;        /**< Maximum number of concurrent rings. */
  static constexpr float START_RHO = -3.75f;  /**< Rho at which rings are spawned (path start). */
  static constexpr float END_RHO = 3.75f;     /**< Rho past which rings are retired (path end). */
  static constexpr float SPACING = 0.3f;      /**< Inter-ring spacing in rho units. */
  /**
   * @brief Rho advanced per frame per unit of the Speed slider.
   * @details Converts the unitless Speed control to per-frame motion along the
   * path.
   */
  static constexpr float RHO_PER_SPEED = 0.009375f;
  static constexpr float PETAL_LOBES = 3.0f;  /**< Petal lobe count per revolution for the radial wobble. */
  static constexpr float PETAL_DEPTH = 0.6f;  /**< Petal wobble depth in rho units. */

  /**
   * @brief One ring on the flow.
   * @details Inactive slots are free for reuse.
   */
  struct Ring {
    float rho;    /**< Position along the log-radial path, in rho units. */
    float hue;    /**< Palette hue selector in [0, 1). */
    bool active;  /**< Whether this slot holds a live ring. */
  };

  Ring rings[MAX_RINGS];          /**< Fixed pool of ring slots, active or free. */
  float gap_accumulator = 0.0f;   /**< Accumulated path travel awaiting the next spawn, in rho units. */

  ProceduralPalette palette;      /**< Color palette sampled by ring hue. */
  Orientation<> orientation;      /**< Shared orientation driven by the timeline rotation. */

  /** @brief Render pipeline: world-space orientation then screen-space anti-aliasing. */
  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;

  Timeline timeline;                  /**< Schedules the orientation rotation and the spawner. */
  Animation::PeriodicTimer spawner;   /**< Per-frame timer that triggers check_spawn(). */

  /**
   * @brief Seeds the timeline and pre-fills the path with rings.
   * @details Adds the looping orientation rotation and the spawner, then
   * pre-fills the entire path with evenly spaced rings so the flow is full from
   * frame zero rather than filling in over time.
   */
  FLASHMEM void init_timeline() {
    timeline.add(0, Animation::Rotation<W>(orientation, UP, PI_F / 4.0f, 160,
                                           ease_mid, true));
    gap_accumulator = 0.0f;
    timeline.add(0, spawner);

    for (float r = END_RHO - 0.01f; r > START_RHO; r -= SPACING) {
      spawn_ring_at_pos(r);
    }
  }

  /**
   * @brief Accumulates this frame's travel and spawns rings as gap opens up.
   * @details Emits a ring each time a full SPACING of gap opens up at the start
   * of the path, keeping ring density constant.
   */
  void check_spawn() {
    float move_dist = params.speed * RHO_PER_SPEED;
    gap_accumulator += move_dist;

    while (gap_accumulator >= SPACING) {
      gap_accumulator -= SPACING;
      spawn_ring_at_pos(START_RHO + gap_accumulator);
    }
  }

  /**
   * @brief Per-instance hue cursor, advanced per spawn.
   * @details Reset in init() so hue assignment stays deterministic for the
   * fixed-seed segmented driver.
   */
  float next_hue = 0.0f;

  /**
   * @brief Claims the first inactive ring slot and initializes it.
   * @param initial_rho Position along the path at which to place the ring, in
   * rho units.
   * @details Assigns the next hue and advances the hue cursor. No-op if all
   * slots are in use.
   */
  void spawn_ring_at_pos(float initial_rho) {
    for (int i = 0; i < MAX_RINGS; ++i) {
      if (!rings[i].active) {
        rings[i].active = true;
        rings[i].rho = initial_rho;
        rings[i].hue = next_hue;
        // Advance the hue cursor by a fixed step that is not a simple fraction
        // of 1, so successive rings spread around the wheel instead of repeating
        // a short cycle (0.13 -> ~7.7 rings per loop, no early alias).
        constexpr float kHueStep = 0.13f;
        next_hue = wrap(next_hue + kHueStep, 1.0f);
        return;
      }
    }
  }

  /**
   * @brief Advances every active ring along rho and draws it.
   * @param canvas Target canvas to render into.
   * @details Retires rings that pass END_RHO; otherwise draws them.
   */
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

  /**
   * @brief Builds and rasterizes one ring.
   * @param canvas Target canvas to render into.
   * @param ring Ring to draw, supplying its rho position and hue.
   * @details Fades the ring by distance from the equator (rho=0), shapes it
   * into petals via a radial wobble, twists it by rho, and projects each sample
   * onto the sphere. Skipped entirely when its opacity is negligible.
   */
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
