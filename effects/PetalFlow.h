/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file PetalFlow.h
 * @brief Concentric rings advanced along a log-radial path and shaped into
 *        flowing petals.
 */

#include <cmath>
#include <algorithm>
#include "core/engine/engine.h"

// Unit-test accessor for the spawn-gap accumulator and hue-cursor wrap invariants.
namespace hs_test {
namespace effects_tests {
struct PetalFlowWhiteBox;
} // namespace effects_tests
} // namespace hs_test

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
   * @brief Constructs the effect, wiring up palette and orientation filters.
   */
  HS_COLD_MEMBER PetalFlow()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})),
        palette({0.029f, 0.029f, 0.029f}, {0.500f, 0.500f, 0.500f},
                {0.461f, 0.461f, 0.461f}, {0.539f, 0.701f, 0.809f}),
        orientation(), filters(Filter::World::Orient(orientation),
                               Filter::Screen::AntiAlias<W, H>()) {}

  /**
   * @brief Registers sliders, clears all rings, and seeds the timeline.
   * @details next_hue resets here so hue assignment is deterministic on every
   * (re)init.
   */
  void init() override {
    register_param("Twist", &params.twist_factor, 0.0f, 5.0f);
    register_param("Speed", &params.speed, 0.0f, SPEED_MAX);
    register_param("Alpha", &params.alpha, 0.0f, 1.0f);
    register_param("Density", &params.density, 0.5f, DENSITY_MAX);

    for (int i = 0; i < MAX_RINGS; ++i) {
      rings[i].active = false;
    }
    next_hue = 0.0f;

    build_shift_table();
    init_timeline();
  }

  /**
   * @brief Steps the timeline (rotation + spawner), then advances and renders
   * all rings.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(pf_timeline_step);
      timeline.step(canvas);
    }
    {
      HS_PROFILE(pf_draw_rings);
      update_and_draw_rings(canvas);
    }
  }

private:
  // Test seam: reaches the spawn-gap accumulator and hue-cursor wrap invariants.
  friend struct ::hs_test::effects_tests::PetalFlowWhiteBox;

  static constexpr int MAX_RINGS =
      65; /**< Maximum number of concurrent rings. */
  static constexpr float START_RHO =
      -3.75f; /**< Rho at which rings are spawned (path start). */
  static constexpr float END_RHO =
      3.75f; /**< Rho past which rings are retired (path end). */
  static constexpr float SPACING =
      0.3f; /**< Base inter-ring spacing in rho units, before the Density scaling. */
  static constexpr float DENSITY_MAX =
      2.5f; /**< Density slider maximum; sets the tightest live spacing, SPACING / DENSITY_MAX. */
  static constexpr float SPEED_MAX =
      20.0f; /**< Speed slider maximum; sets the longest per-frame travel. */
  /**
   * @brief Rho advanced per frame per unit of the Speed slider.
   * @details Converts the unitless Speed control to per-frame motion along the
   * path.
   */
  static constexpr float RHO_PER_SPEED = SPACING / 32.0f;
  /**
   * @brief Peak concurrent rings, at the tightest spacing and the fastest flow.
   * @details check_spawn() runs before retirement and emits one ring per gap the
   * frame's travel opens, so the newest ring can sit a whole frame of travel
   * below START_RHO while the oldest still sits on END_RHO. Rings hold exactly
   * one gap apart, so the peak is one per whole gap across that span plus the
   * ring closing it.
   */
  static constexpr int RINGS_ON_PATH =
      static_cast<int>((END_RHO - START_RHO + SPEED_MAX * RHO_PER_SPEED) *
                       DENSITY_MAX / SPACING) +
      1;
  static_assert(MAX_RINGS >= RINGS_ON_PATH,
                "PetalFlow ring pool is smaller than the path carries at "
                "maximum Speed x Density; raise MAX_RINGS or retune SPACING/"
                "END_RHO/START_RHO/DENSITY_MAX/SPEED_MAX");
  static constexpr float FADE_START_RHO =
      2.5f; /**< |rho| at which the pole fade begins. */
  /**
   * @brief Rho span over which a ring fades from full opacity to zero.
   * @details The fade completes short of END_RHO, so rings travel the remaining
   * |rho| fully transparent and are never drawn there.
   */
  static constexpr float FADE_WIDTH = 1.0f;
  static_assert(FADE_START_RHO + FADE_WIDTH <= END_RHO &&
                    -(FADE_START_RHO + FADE_WIDTH) >= START_RHO,
                "PetalFlow pole fade must complete inside the ring path");
  static constexpr float PETAL_LOBES =
      3.0f; /**< Petal lobe count per revolution for the radial wobble. */
  static constexpr float PETAL_DEPTH =
      0.6f; /**< Petal wobble depth in rho units. */
  static constexpr int NUM_SAMPLES =
      W / 2; /**< Angular samples drawn per ring. */
  // draw_ring stages one ring's NUM_SAMPLES fragments in scratch_a at a time,
  // alongside rasterize's own sub-step cache.
  static_assert(NUM_SAMPLES * sizeof(Fragment) +
                        Plot::rasterize_scratch_a_bytes<W>() <=
                    DEFAULT_SCRATCH_A_SIZE,
                "PetalFlow ring staging exceeds the default scratch_a budget; "
                "retune NUM_SAMPLES or carve a larger scratch arena");

  /**
   * @brief One ring on the flow.
   * @details Inactive slots are free for reuse.
   */
  struct Ring {
    float rho;   /**< Position along the log-radial path, in rho units. */
    float hue;   /**< Palette hue selector in [0, 1). */
    bool active; /**< Whether this slot holds a live ring. */
  };

  Ring rings[MAX_RINGS]; /**< Fixed pool of ring slots, active or free. */

  /**
   * @brief Per-sample exp(radial wobble), cached once in init().
   * @details The petal shift depends only on the normalized angle, so
   * exp(shift) is identical for every ring and frame; caching it leaves a single
   * exp(rho) per ring instead of one exp per sample.
   */
  float exp_shift[NUM_SAMPLES];
  float gap_accumulator =
      0.0f; /**< Accumulated path travel awaiting the next spawn, in rho units. */
  float next_hue =
      0.0f; /**< Per-instance hue cursor, advanced per spawn; reset in init() so hue assignment stays deterministic for the fixed-seed segmented driver. */

#ifdef HS_TEST_BUILD
  int spawns = 0;         /**< Rings placed into a free slot, cumulative. */
  int dropped_spawns = 0; /**< Spawn requests that found no free slot. */
#endif

  ProceduralPalette palette; /**< Color palette sampled by ring hue. */
  Orientation<>
      orientation; /**< Shared orientation driven by the timeline rotation. */

  /** @brief Render pipeline: world-space orientation then screen-space anti-aliasing. */
  Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>>
      filters;

  Timeline timeline; /**< Schedules the orientation rotation and the spawner. */

  /**
   * @brief Precomputes the geometry-static exp(petal wobble) per sample.
   * @details The radial wobble is a function of the normalized angle alone, so
   * its exponential is constant across rings and frames.
   */
  HS_COLD_MEMBER void build_shift_table() {
    for (int i = 0; i < NUM_SAMPLES; ++i) {
      float t_norm = static_cast<float>(i) / NUM_SAMPLES;
      float shift =
          PETAL_DEPTH * std::abs(fast_sinf(PETAL_LOBES * PI_F * t_norm));
      exp_shift[i] = expf(shift);
    }
  }

  /**
   * @brief Seeds the timeline and pre-fills the path with rings.
   * @details Adds the looping orientation rotation and the spawner, then
   * pre-fills the entire path with evenly spaced rings so the flow is full from
   * frame zero rather than filling in over time.
   */
  HS_COLD_MEMBER void init_timeline() {
    timeline.add(0, Animation::Rotation<W>(orientation, UP, PI_F / 4.0f, 160,
                                           ease_linear, true));
    gap_accumulator = 0.0f;
    timeline.add(0, Animation::PeriodicTimer(
                        1, [this](Canvas &) { this->check_spawn(); }, true));

    // Pre-fill the path with evenly spaced rings so frame zero looks like a
    // running flow. Start one epsilon inside END_RHO so the oldest ring sits just
    // short of the retire boundary (draw_frame() culls rho > END_RHO).
    for (float r = END_RHO - 0.01f; r > START_RHO; r -= spacing()) {
      spawn_ring_at_pos(r);
    }
  }

  /**
   * @brief Rho advanced this frame: the Speed slider scaled into per-frame
   * motion. Spawn and render must step by the identical amount, so both read it
   * here.
   */
  float move_dist() const { return params.speed * RHO_PER_SPEED; }

  /** @brief Live inter-ring spacing: the base SPACING scaled by the Density slider. */
  float spacing() const { return SPACING / params.density; }

  /**
   * @brief Accumulates this frame's travel and spawns rings as gap opens up.
   * @details Emits a ring each time a full SPACING of gap opens up at the start
   * of the path, keeping ring density constant.
   */
  void check_spawn() {
    const float move = move_dist();
    const float gap = spacing();
    gap_accumulator += move;

    // update_and_draw_rings() advances every ring this same frame, so spawn one
    // move short to land the first drawn position at START_RHO + the residual.
    while (gap_accumulator >= gap) {
      gap_accumulator -= gap;
      spawn_ring_at_pos(START_RHO + gap_accumulator - move);
    }
  }

  /**
   * @brief Claims the first inactive ring slot and initializes it.
   * @param initial_rho Position along the path at which to place the ring, in
   * rho units.
   * @details Assigns the next hue and advances the hue cursor. No-op if all
   * slots are in use, which RINGS_ON_PATH sizes the pool to prevent.
   */
  HS_COLD_MEMBER void spawn_ring_at_pos(float initial_rho) {
    for (int i = 0; i < MAX_RINGS; ++i) {
      if (!rings[i].active) {
        rings[i].active = true;
        rings[i].rho = initial_rho;
        rings[i].hue = next_hue;
        constexpr float HUE_STEP = 0.13f;
        next_hue = wrap(next_hue + HUE_STEP, 1.0f);
#ifdef HS_TEST_BUILD
        ++spawns;
#endif
        return;
      }
    }
#ifdef HS_TEST_BUILD
    ++dropped_spawns;
#endif
  }

  /**
   * @brief Advances every active ring along rho and draws it.
   * @param canvas Target canvas to render into.
   * @details Retires rings that pass END_RHO; otherwise draws them. A sub-LSB
   * Alpha paints nothing, so rings still advance and retire but are not
   * rasterized, and the flow resumes in place when Alpha rises.
   */
  void update_and_draw_rings(Canvas &canvas) {
    const float move = move_dist();
    const bool visible = params.alpha >= MIN_VISIBLE_ALPHA;
    for (int i = 0; i < MAX_RINGS; ++i) {
      if (!rings[i].active)
        continue;
      rings[i].rho += move;
      if (rings[i].rho > END_RHO) {
        rings[i].active = false;
        continue;
      }
      if (visible)
        draw_ring(canvas, rings[i]);
    }
  }

  /**
   * @brief Builds and rasterizes one ring.
   * @param canvas Target canvas to render into.
   * @param ring Ring to draw, supplying its rho position and hue.
   * @details Fades the ring by distance from the equator (rho=0), shapes it
   * into petals via a radial wobble, twists it by rho, and projects each sample
   * onto the sphere. Skipped entirely once the pole fade has taken its opacity
   * to ~0; the Alpha slider is gated by the caller.
   */
  void draw_ring(Canvas &canvas, const Ring &ring) {
    float dist = std::abs(ring.rho);
    float opacity = 1.0f;
    if (dist > FADE_START_RHO) {
      opacity = std::max(0.0f, 1.0f - (dist - FADE_START_RHO) / FADE_WIDTH);
    }

    if (opacity <= 0.01f)
      return;
    float effective_opacity = opacity * params.alpha;

    constexpr int num_samples = NUM_SAMPLES;
    const float step = 2.0f * PI_F / num_samples;

    Color4 base_col = palette.get(ring.hue).color;
    base_col.alpha = effective_opacity;

    float twist_angle = (ring.rho / spacing()) * params.twist_factor;
    const float exp_rho = expf(ring.rho);

    ScratchScope scratch_a_guard(scratch_arena_a);
    Fragments fragments;
    fragments.bind(scratch_arena_a, num_samples);
    {
      HS_PROFILE(pf_ring_build);
      for (int i = 0; i < num_samples; ++i) {
        float theta = i * step;

        float final_theta = theta + twist_angle;

        // Inverse stereographic projection of (radius exp(rho + wobble), angle
        // final_theta) onto the unit sphere, poles on +/-Z; core's inv_stereo()
        // is +/-Y-poled. exp_shift caches the geometry-static exp(wobble).
        float R = exp_rho * exp_shift[i];

        float r2 = R * R;
        float denom = 1.0f + r2;
        float x = 2.0f * R * fast_cosf(final_theta) / denom;
        float y = 2.0f * R * fast_sinf(final_theta) / denom;
        float z = (r2 - 1.0f) / denom;

        Fragment f;
        f.pos = Vector(x, y, z);
        fragments.push_back(f);
      }
    }

    auto fragment_shader = [&](const Vector &, Fragment &f) {
      f.color = base_col;
    };

    {
      HS_PROFILE(pf_ring_scan);
      Plot::rasterize<W, H>(filters, canvas, fragments, fragment_shader,
                            {.close_loop = true});
    }
  }

  /**
   * @brief User-tunable controls, registered as sliders in init().
   */
  struct Params {
    float twist_factor =
        0.35f; /**< Swirl strength applied per unit of rho spacing. */
    float speed =
        2.5f; /**< Unitless flow rate; scaled by RHO_PER_SPEED into per-frame motion. */
    float alpha = 0.2f; /**< Global opacity multiplier. */
    float density =
        1.0f; /**< Ring packing; live spacing is SPACING / density. */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(PetalFlow)
