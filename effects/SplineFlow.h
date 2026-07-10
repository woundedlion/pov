/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include "core/engine/engine.h"

/**
 * @brief Closed geodesic spline whose control points drift over the sphere,
 *        drawn as a fading trail that sweeps around a moving draw head.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Control points follow independent random walks across the sphere.
 *          The spline is colored from a baked palette by spline progress.
 */
template <int W, int H> class SplineFlow : public Effect {
public:
  static constexpr int MAX_TRAILS = 30000;
  static constexpr int MAX_POINTS = 12;
  static constexpr int lifetime = 60;

  // The trail pool dominates the persistent partition (~234 KiB at 30000×8 B).
  // Effect keeps the default arena split, so it must fit the device budget.
  static constexpr size_t TRAIL_BYTES =
      MAX_TRAILS * sizeof(typename Filter::World::Trails<W, MAX_TRAILS>::Item);
  static constexpr size_t PERSISTENT_BUDGET =
      DEVICE_GLOBAL_ARENA_SIZE - DEFAULT_SCRATCH_A_SIZE - DEFAULT_SCRATCH_B_SIZE;
  static_assert(TRAIL_BYTES <= PERSISTENT_BUDGET,
                "SplineFlow trail pool exceeds the default persistent "
                "partition; retune MAX_TRAILS or carve arenas");

  /**
   * @brief Tunable parameters exposed to the UI sliders.
   */
  struct Params {
    float tension = 0.3f;   /**< Catmull-Rom spline tension in [0, 1]. */
    float speed = 0.05f;    /**< Draw-head phase rate per step. */
    float drift = 0.5f;     /**< Control-point random-walk speed scale in [0, 1]. */
    float num_points = 6.0f; /**< Number of active control points in [4, 12]. */
    float alpha = 0.6f;     /**< Peak trail opacity in [0.1, 1]. */
  } params;

  /**
   * @brief Constructs the effect and wires up the trail/orient/anti-alias filters.
   */
  FLASHMEM SplineFlow()
      : Effect(W, H, {.strobe = true, .full_frame = decltype(filters)::any_crosses_segments}),
        filters(Filter::World::Trails<W, MAX_TRAILS>(lifetime),
                              Filter::World::Orient<W>(orientation),
                              Filter::Screen::AntiAlias<W, H>()) {}

  /**
   * @brief Allocates trail storage, bakes the palette, registers sliders, and
   *        seeds the control-point/orientation/draw-head animations.
   * @details Runs once before the first frame; sets up persistent state on the
   *          arena and the timeline.
   */
  void init() override {
    filters.template get<Filter::World::Trails<W, MAX_TRAILS>>().init_storage(
        persistent_arena);

    baked_palette_.bake(persistent_arena, Palettes::LAVENDER_LAKE);

    register_param("Tension", &params.tension, 0.0f, 1.0f);
    register_param("Speed", &params.speed, 0.01f, 0.2f);
    register_param("Drift", &params.drift, 0.0f, 1.0f);
    register_param("Num Pts", &params.num_points, 4.0f, 12.0f);
    register_param("Alpha", &params.alpha, 0.1f, 1.0f);

    // Animated control points driven by independent random walks. All MAX_POINTS
    // walks step every frame regardless of the live Num Pts, so raising Num Pts
    // mid-run reveals already-drifted points rather than ones frozen at their seed.
    for (int i = 0; i < MAX_POINTS; ++i) {
      point_orientations[i].set(
          make_rotation(random_vector(), hs::rand_f() * 2 * PI_F));
      point_walks_[i] = timeline.add_get(
          0, Animation::RandomWalk<W>(point_orientations[i], random_vector(),
                                      point_noises_[i],
                                      {params.drift * 0.02f, 0.15f, 0.03f}));
    }

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS, orient_noise_,
                        Animation::RandomWalk<W>::Options::Languid()));

    timeline.add(0, Animation::Driver(draw_head, &params.speed, 1.0f, true));
  }

  /**
   * @brief Advances animations, rebuilds the closed spline from the drifting
   *        control points, draws it, and flushes faded trails.
   * @details Called once per frame; also live-applies the Drift slider to the
   *          control-point walk speeds.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    apply_if_changed(params.drift, last_drift_, [&](float drift) {
      float walk_speed = drift * 0.02f;
      for (auto *w : point_walks_)
        if (w)
          w->set_speed(walk_speed);
    });

    int n = hs::clamp(static_cast<int>(params.num_points), 4, MAX_POINTS);

    ScratchScope frag_guard(scratch_arena_a);
    Fragments control_points;
    control_points.bind(scratch_arena_a, n);
    for (int i = 0; i < n; ++i) {
      Fragment f;
      f.pos = point_orientations[i].orient(Y_AXIS);
      control_points.push_back(f);
    }

    auto shader = [&](const Vector &, Fragment &f) {
      // Color by progress, fade alpha by distance from draw head
      float dist = shortest_distance(f.v0, draw_head, 1.0f);
      float brightness = quintic_kernel(1.0f - dist * 4.0f);
      f.color = baked_palette_.get(f.v0);
      f.color.alpha = brightness * params.alpha;
    };

    Plot::SplineChain::draw<W, H>(filters, canvas, control_points,
                                  params.tension, shader, {},
                                  true, // closed loop
                                  24,   // samples per segment
                                  SplineMode::Geodesic);

    filters.flush(
        canvas,
        [this](const Vector &, float t) {
          Color4 c = baked_palette_.get(t);
          c.alpha *= (1.0f - t);
          return c;
        },
        1.0f);
  }

private:
  float draw_head = 0.0f; /**< Spline progress phase of the draw head in [0, 1). */
  Orientation<> orientation; /**< Global orientation animated by a languid walk. */
  Orientation<> point_orientations[MAX_POINTS]; /**< Per-control-point orientations. */
  Animation::RandomWalk<W> *point_walks_[MAX_POINTS] = {}; /**< Handles for live Drift updates. */
  float last_drift_ = -1.0f; /**< Last applied Drift value for change detection. */
  /**
   * @brief One noise generator per control-point walk.
   * @details The RandomWalk ctor stamps its frequency and seed onto the
   *          generator it is handed, so each walk needs its own to keep its
   *          intended frequency and an independent noise field.
   */
  std::array<FastNoiseLite, MAX_POINTS> point_noises_;
  FastNoiseLite orient_noise_; /**< Noise field for the global orientation walk. */
  Timeline timeline; /**< Drives all animations each frame. */
  BakedPalette baked_palette_; /**< 256-entry palette LUT for coloring. */

  /** @brief Render pipeline chaining the trail, orient, and anti-alias filters. */
  Pipeline<W, H, Filter::World::Trails<W, MAX_TRAILS>, Filter::World::Orient<W>,
           Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(SplineFlow)
