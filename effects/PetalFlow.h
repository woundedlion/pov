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

#include "core/platform/effect_identity.h"
#include "core/animation/orientation.h"
#include <cmath>
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
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY =
      hs::effect_identities::PetalFlow::VALUE;
  /// Canonical machine ID.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;

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
   * @details Runs once per effect instance: init_timeline() adds its events
   * unconditionally.
   */
  HS_COLD_MEMBER void init() override {
    register_param("Twist", &params.twist_factor, 0.0f, 5.0f);
    register_param("Speed", &params.speed, 0.0f, SPEED_MAX);
    register_param("Alpha", &params.alpha, 0.0f, 1.0f);
    register_param("Density", &params.density, 0.5f, DENSITY_MAX);

    for (int i = 0; i < MAX_RINGS; ++i) {
      rings[i].active = false;
    }
    next_hue = 0.0f;
    history_ready = false;

    build_shift_table();
    init_timeline();
  }

  /**
   * @brief Steps the timeline (rotation + spawner), then advances and renders
   * all rings.
   */
  HS_HOT_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(pf_timeline_step);
      timeline.step(canvas);
      if (!history_ready)
        orientation.collapse();
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
  /** @brief Rho advanced per frame per unit of the Speed slider. */
  static constexpr float RHO_PER_SPEED = SPACING / 32.0f;
  /**
   * @brief Peak concurrent rings, at the tightest spacing and the fastest flow.
   * @details Spawning precedes retirement, so the span includes one frame of
   * travel below START_RHO.
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
  /**
   * @brief Pole-fade envelope at or below which a ring is not drawn.
   * @details A cut on the fade envelope, before the Alpha slider scales it, so
   * it is not comparable to the per-sample MIN_ENCODABLE_ALPHA floor.
   */
  static constexpr float MIN_FADE_OPACITY = 0.01f;
  static constexpr float PETAL_LOBES =
      3.0f; /**< Petal lobe count per revolution for the radial wobble. */
  static constexpr float PETAL_DEPTH =
      0.6f; /**< Petal wobble depth in rho units. */
  static constexpr int NUM_SAMPLES =
      W / 2; /**< Angular samples drawn per ring. */
  // scratch_a: one ring's fragments plus rasterize's sub-step cache.
  static_assert(NUM_SAMPLES * sizeof(Fragment) +
                        Plot::rasterize_scratch_a_bytes<W>() <=
                    DEFAULT_SCRATCH_A_SIZE,
                "PetalFlow ring staging exceeds the default scratch_a budget; "
                "retune NUM_SAMPLES or carve a larger scratch arena");

  /** @brief Ring staging, temporal controls and replay steps. */
  static constexpr size_t TIMED_SCRATCH_A_BYTES =
      NUM_SAMPLES * sizeof(Fragment) +
      std::max(NUM_SAMPLES * sizeof(math::PixelCoords),
               2 * sizeof(math::Vector) +
                   Plot::MAX_RIGID_CONTROLS * sizeof(Plot::TimedRigidNode)) +
      Plot::rasterize_step_budget<W>() * sizeof(float);
  static_assert(TIMED_SCRATCH_A_BYTES <= DEFAULT_SCRATCH_A_SIZE,
                "PetalFlow temporal lattice exceeds scratch_a");
  static_assert(NUM_SAMPLES * sizeof(Plot::DenseControl) <=
                    DEFAULT_SCRATCH_B_SIZE,
                "PetalFlow dense controls exceed scratch_b");

  /**
   * @brief One ring on the flow.
   * @details Inactive slots are free for reuse.
   */
  struct Ring {
    float rho;   /**< Position along the log-radial path, in rho units. */
    float hue;   /**< Palette hue selector in [0, 1). */
    bool active; /**< Whether this slot holds a live ring. */
    bool born;   /**< Whether this slot has not yet been drawn. */
  };

  Ring rings[MAX_RINGS]{}; /**< Fixed pool of ring slots, active or free. */

  /** @brief Per-sample exp(radial wobble), cached once in init(). */
  float exp_shift[NUM_SAMPLES];
  float gap_accumulator =
      0.0f; /**< Accumulated path travel awaiting the next spawn, in rho units. */
  float next_hue =
      0.0f; /**< Per-instance hue cursor, advanced per spawn; reset in init(). */

#if HS_ENABLE_TEST_HOOKS
  int spawns = 0; /**< Rings placed into a free slot, cumulative. */
#endif

  ProceduralPalette palette; /**< Color palette sampled by ring hue. */
  math::Orientation<>
      orientation; /**< Shared orientation driven by the timeline rotation. */

  /** @brief Render pipeline: world-space orientation then screen-space anti-aliasing. */
  Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>>
      filters;

  Timeline timeline; /**< Schedules the orientation rotation and the spawner. */

  /** @brief Precomputes the geometry-static exp(petal wobble) per sample. */
  HS_COLD_MEMBER void build_shift_table() {
    for (int i = 0; i < NUM_SAMPLES; ++i) {
      float t_norm = static_cast<float>(i) / NUM_SAMPLES;
      float shift =
          PETAL_DEPTH *
          std::abs(math::fast_sinf(PETAL_LOBES * math::PI_F * t_norm));
      exp_shift[i] = expf(shift);
    }
  }

  /**
   * @brief Seeds the timeline and pre-fills the path with evenly spaced rings.
   */
  HS_COLD_MEMBER void init_timeline() {
    timeline.add(
        0, Animation::Rotation<W>(orientation, math::UP, math::PI_F / 4.0f, 160,
                                  math::ease_linear, {.repeat = true}));
    gap_accumulator = 0.0f;
    timeline.add(
        0, Animation::PeriodicTimer(
               1, [this](Canvas &) { this->check_spawn(); }, {.repeat = true}));

    float last_r = END_RHO - 0.01f;
    for (float r = last_r; r > START_RHO; r -= spacing()) {
      spawn_ring_at_pos(r);
      last_r = r;
    }
    gap_accumulator = last_r - START_RHO;
  }

  /** @brief Rho advanced this frame; spawn and render both step by it. */
  float move_dist() const { return params.speed * RHO_PER_SPEED; }

  /** @brief Live inter-ring spacing: the base SPACING scaled by the Density slider. */
  float spacing() const { return SPACING / params.density; }

  /**
   * @brief Accumulates this frame's travel and spawns rings as gap opens up.
   * @details Emits a ring each time a full live spacing() of gap opens up at the start
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
   * @details Assigns the next hue and advances the hue cursor; traps if full.
   */
  HS_COLD_MEMBER void spawn_ring_at_pos(float initial_rho) {
    for (int i = 0; i < MAX_RINGS; ++i) {
      if (!rings[i].active) {
        rings[i].active = true;
        rings[i].born = true;
        rings[i].rho = initial_rho;
        rings[i].hue = next_hue;
        constexpr float HUE_STEP = 0.13f;
        next_hue = math::wrap(next_hue + HUE_STEP, 1.0f);
#if HS_ENABLE_TEST_HOOKS
        ++spawns;
#endif
        return;
      }
    }
    HS_CHECK(false, "PetalFlow: ring pool full");
  }

  /**
   * @brief Advances every active ring along rho and draws it.
   * @param canvas Target canvas to render into.
   * @details Retires rings that pass END_RHO; otherwise draws them. Below one
   * slider LSB of master alpha, rings advance and retire without rasterizing.
   */
  HS_HOT_FLASH_MEMBER void update_and_draw_rings(Canvas &canvas) {
    const bool continuous =
        history_ready && history_params.twist_factor == params.twist_factor &&
        history_params.speed == params.speed &&
        history_params.alpha == params.alpha &&
        history_params.density == params.density;
    const float move = move_dist();
    const bool visible = params.alpha >= MIN_VISIBLE_ALPHA;
    for (int i = 0; i < MAX_RINGS; ++i) {
      if (!rings[i].active)
        continue;
      const float previous_rho = rings[i].rho;
      rings[i].rho += move;
      if (rings[i].rho > END_RHO) {
        rings[i].active = false;
        continue;
      }
      if (visible)
        draw_ring(canvas, rings[i],
                  continuous && !rings[i].born ? previous_rho : rings[i].rho);
      rings[i].born = false;
    }
    history_params = params;
    history_ready = true;
  }

  /** @brief Pole-fade envelope at a log-radial position.
   * @param rho Log-radial position.
   * @return Opacity before the master alpha. */
  HS_HOT_INLINE static float pole_opacity(float rho) {
    const float dist = std::abs(rho);
    if (dist > FADE_START_RHO)
      return fmaxf(0.0f, 1.0f - (dist - FADE_START_RHO) / FADE_WIDTH);
    return 1.0f;
  }

  /** @brief Evaluates a ring's prepared stereographic geometry.
   * @param index Angular control point.
   * @param exp_rho Exponential of the log-radial position.
   * @param twist_angle Ring twist in radians.
   * @return Unit-sphere control point. */
  HS_HOT_FLASH_MEMBER math::Vector ring_vertex(size_t index, float exp_rho,
                                               float twist_angle) const {
    constexpr float STEP = 2.0f * math::PI_F / NUM_SAMPLES;
    const float theta = index * STEP;
    const float final_theta = theta + twist_angle;
    const float radius = exp_rho * exp_shift[index];
    const float r2 = radius * radius;
    const float denom = 1.0f + r2;
    const float x = 2.0f * radius * math::fast_cosf(final_theta) / denom;
    const float y = 2.0f * radius * math::fast_sinf(final_theta) / denom;
    const float z = (r2 - 1.0f) / denom;
    return math::Vector(x, y, z);
  }

  /** @brief Builds and rasterizes a ring's two temporal endpoints.
   * @param canvas Target canvas.
   * @param ring Current ring position and hue.
   * @param previous_rho Position drawn in the preceding frame. */
  HS_HOT_FLASH_MEMBER void draw_ring(Canvas &canvas, const Ring &ring,
                                     float previous_rho) {
    const bool timed = Temporal::temporal_sampling_enabled();
    if (!timed)
      previous_rho = ring.rho;
    const float opacity = pole_opacity(ring.rho);
    if (opacity <= MIN_FADE_OPACITY &&
        (!timed || pole_opacity(previous_rho) <= MIN_FADE_OPACITY))
      return;
    Color4 base_col = palette.get(ring.hue).color;
    base_col.alpha = opacity * params.alpha;
    // One spacing() division per branch keeps the still twist's rounding.
    float twist_angle;
    float previous_twist;
    if (timed) {
      twist_angle = (ring.rho / spacing()) * params.twist_factor;
      previous_twist = (previous_rho / spacing()) * params.twist_factor;
    } else {
      twist_angle = (ring.rho / spacing()) * params.twist_factor;
      previous_twist = twist_angle;
    }
    const float exp_rho = expf(ring.rho);
    const float previous_exp =
        previous_rho == ring.rho ? exp_rho : expf(previous_rho);
    const float callback_bound =
        fabsf(ring.rho - previous_rho) + fabsf(twist_angle - previous_twist);

    ScratchScope scratch_a_guard(scratch_arena_a);
    Fragments fragments;
    fragments.bind(scratch_arena_a, NUM_SAMPLES);
    {
      HS_PROFILE(pf_ring_build);
      for (int i = 0; i < NUM_SAMPLES; ++i) {
        Fragment f;
        f.pos = ring_vertex(i, exp_rho, twist_angle);
        fragments.push_back(f);
      }
    }
    auto geometry = [&](size_t index, float t) HS_HOT_INLINE {
      return t == 1.0f ? fragments[index].pos
                       : ring_vertex(index, previous_exp, previous_twist);
    };
    auto fragment_shader = [&](const math::Vector &, Fragment &f,
                               float t) HS_HOT_INLINE {
      f.color = base_col;
      if (t != 1.0f)
        f.color.alpha =
            pole_opacity(previous_rho + (ring.rho - previous_rho) * t) *
            params.alpha;
    };
    {
      HS_PROFILE(pf_ring_scan);
      Plot::rasterize<W, H, {.timed = true}>(
          filters, canvas, fragments, fragment_shader,
          {.loop = Plot::RasterLoop::closed()},
          Plot::GeometryCallback::interpolating(geometry), callback_bound);
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

  Params history_params;      /**< Controls at the preceding frame endpoint. */
  bool history_ready = false; /**< Whether a preceding endpoint is available. */
};
