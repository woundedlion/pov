/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include "../effects_engine.h"

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

  void init() override {
    static_cast<Filter::World::Trails<W, MAX_TRAILS> &>(filters).init_storage(
        persistent_arena);

    registerParam("Tension", &params.tension, 0.0f, 1.0f);
    registerParam("Speed", &params.speed, 0.01f, 0.2f);
    registerParam("Drift", &params.drift, 0.0f, 1.0f);
    registerParam("Num Pts", &params.num_points, 4.0f, 12.0f);
    registerParam("Alpha", &params.alpha, 0.1f, 1.0f);

    // Animated control points driven by independent random walks
    for (int i = 0; i < MAX_POINTS; ++i) {
      point_orientations[i].set(
          make_rotation(random_vector(), hs::rand_f() * 2 * PI_F));
      timeline.add(0, Animation::RandomWalk<W>(
                          point_orientations[i], random_vector(), noise,
                          {params.drift * 0.02f, 0.15f, 0.03f}));
    }

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));

    // Phase accumulator for draw-head position
    timeline.add(0, Animation::Driver(draw_head, params.speed, true));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

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

    auto shader = [&](const Vector &v, Fragment &f) {
      // Color by progress, fade alpha by distance from draw head
      float dist = shortest_distance(f.v0, draw_head, 1.0f);
      float brightness = quintic_kernel(1.0f - dist * 4.0f);
      f.color = palette.get(f.v0);
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
        [](const Vector &, float t) {
          Color4 c = Palettes::lavenderLake.get(t);
          c.alpha *= (1.0f - t);
          return c;
        },
        1.0f);
  }

private:
  float draw_head = 0.0f;
  Orientation<W> orientation;
  Orientation<W> point_orientations[MAX_POINTS];
  FastNoiseLite noise;
  Timeline<W> timeline;
  ProceduralPalette palette = Palettes::lavenderLake;

  Pipeline<W, H, Filter::World::Trails<W, MAX_TRAILS>, Filter::World::Orient<W>,
           Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "../effect_registry.h"
REGISTER_EFFECT(SplineFlow)
