/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include "core/effects_engine.h"

// Closed geodesic spline whose control points drift over the sphere via
// independent random walks, drawn as a fading trail that sweeps around a moving
// draw head. Colored from a baked palette by spline progress.
template <int W, int H> class SplineFlow : public Effect {
public:
  static constexpr int MAX_TRAILS = 30000;
  static constexpr int MAX_POINTS = 12;
  static constexpr int lifetime = 60;

  struct Params {
    float tension = 0.3f;
    float speed = 0.05f;
    float drift = 0.5f;
    float num_points = 6.0f;
    float alpha = 0.6f;
  } params;

  FLASHMEM SplineFlow()
      : Effect(W, H), filters(Filter::World::Trails<W, MAX_TRAILS>(lifetime),
                              Filter::World::Orient<W>(orientation),
                              Filter::Screen::AntiAlias<W, H>()) {}

  // Allocate trail storage, bake the palette, register sliders, and seed the
  // control-point/orientation/draw-head animations on the timeline.
  void init() override {
    filters.template get<Filter::World::Trails<W, MAX_TRAILS>>().init_storage(
        persistent_arena);

    // Bake the immutable palette once into a 256-entry LUT. The trail flush
    // evaluates the palette per live trail item (~25k/frame); a ProceduralPalette
    // get() costs 3 cosf + sRGB-LUT interp + a virtual call each, whereas
    // BakedPalette::get is a table lookup + lerp. Reused by the spline shader too.
    baked_palette_.bake(persistent_arena, Palettes::lavenderLake);

    registerParam("Tension", &params.tension, 0.0f, 1.0f);
    registerParam("Speed", &params.speed, 0.01f, 0.2f);
    registerParam("Drift", &params.drift, 0.0f, 1.0f);
    registerParam("Num Pts", &params.num_points, 4.0f, 12.0f);
    registerParam("Alpha", &params.alpha, 0.1f, 1.0f);

    // Animated control points driven by independent random walks. All walks
    // are infinite (never removed), so the timeline never compacts and these
    // add_get() handles stay valid for live Drift updates.
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

    // Phase accumulator for draw-head position. A bound Driver pulls the Speed
    // slider each step, so it stays live with no retained handle or re-sync.
    timeline.add(0, Animation::Driver(draw_head, &params.speed, 1.0f, true));
  }

  bool show_bg() const override { return false; }

  // Advance animations, rebuild the closed spline from the drifting control
  // points, draw it through the trail pipeline, and flush faded trails.
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    // Live-apply the Drift slider to the control-point walk speeds.
    apply_if_changed(params.drift, last_drift_, [&](float drift) {
      float walk_speed = drift * 0.02f;
      for (auto *w : point_walks_)
        if (w)
          w->set_speed(walk_speed);
    });

    int n = hs::clamp(static_cast<int>(params.num_points), 4, MAX_POINTS);

    // Build control points from animated orientations
    ScratchScope _frag(scratch_arena_a);
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

    // Flush trails
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
  float draw_head = 0.0f;
  Orientation<> orientation;
  Orientation<> point_orientations[MAX_POINTS];
  Animation::RandomWalk<W> *point_walks_[MAX_POINTS] = {};
  float last_drift_ = -1.0f;
  // One noise generator per control-point walk. The RandomWalk ctor stamps its
  // frequency and seed onto the generator it is handed, so each walk needs its
  // own to keep its intended frequency and an independent noise field.
  std::array<FastNoiseLite, MAX_POINTS> point_noises_;
  FastNoiseLite orient_noise_;
  Timeline timeline;
  BakedPalette baked_palette_;

  Pipeline<W, H, Filter::World::Trails<W, MAX_TRAILS>, Filter::World::Orient<W>,
           Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(SplineFlow)
