/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

// Scatters polygon "stars" over a Fibonacci spiral on the sphere, then pushes
// every point through an evolving Möbius warp so the field drifts and inflates.
// A Languid RandomWalk slowly reorients the whole field. W/H are the cubemap
// face resolution.
template <int W, int H> class GnomonicStars : public Effect {
public:
  // Live-tunable controls: star count, per-star radius, polygon side count, and
  // a toggle that draws each star's bounding box for debugging.
  struct Params {
    float points = 600.0f;
    float star_radius = 0.02f;
    float star_sides = 4.0f;
    bool debug_bb = false;
  } params;

  FLASHMEM GnomonicStars()
      : Effect(W, H), orientation(), timeline(), transformer(timeline) {}

  // Exposes the user params and arms the timeline: an infinite Möbius warp whose
  // speed is bound to a slider, plus a Languid RandomWalk that reorients the
  // field.
  void init() override {
    registerParam("Points", &params.points, 100.0f, 2000.0f);
    registerParam("Radius", &params.star_radius, 0.01f, 0.1f);
    registerParam("Sides", &params.star_sides, 3.0f, 8.0f);
    registerParam("Debug BB", &params.debug_bb);

    // Spawn the evolving warp and retain its animation reference for the live
    // "Warp Speed" param. The warp is infinite and added before any other
    // timeline event, so pinning it is the safe retained-handle path: if a
    // future edit ever inserts a finite event ahead of it, compaction traps
    // loudly rather than dangling this pointer the GUI reads every frame.
    auto *anim = transformer.spawn_pinned(0, 0.5f, 0.035f);
    // A pinned spawn into a fresh transformer always has a free slot, so a null
    // here is a structural bug, not a runtime condition. Trap rather than
    // silently drop the "Warp Speed" slider.
    HS_CHECK(anim, "GnomonicStars: pinned warp spawn must succeed");
    registerParam("Warp Speed", &anim->speed, 0.0f, 1.0f);

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, UP, noise,
                        Animation::RandomWalk<W>::Options::Languid()));
  }

  bool show_bg() const override { return true; }

  // Advances the timeline, then draws each spiral point as a star whose color is
  // a Y-based gradient and whose basis carries the current orientation and warp.
  void draw_frame() override {
    Canvas canvas(*this);

    timeline.step(canvas);

    // Tints each fragment along a Mango Peel gradient keyed to its world-space Y.
    auto fragment_shader = [](const Vector &p, Fragment &frag) {
      float t = (p.y + 1.0f) * 0.5f; // map Y in [-1, 1] to gradient t in [0, 1]
      Color4 c = Palettes::mangoPeel.get(t);
      frag.color = c;
    };

    const int points = (int)params.points;
    const float radius = params.star_radius;
    const int sides = (int)params.star_sides;

    for (int i = 0; i < points; i++) {
      Vector v = fib_spiral(points, 0.0f, i);

      v = transformer.transform(v);

      // Build the basis at the star position. make_basis() rotates its normal
      // by the orientation, so passing the raw warp output applies the latest
      // orientation exactly once — orienting v separately first would rotate it
      // twice (q*q*v) and drift the field at double the RandomWalk rate.
      Basis basis = make_basis(orientation.get(), v);

      Scan::Star::draw<W, H>(filters, canvas, basis, radius, sides,
                             fragment_shader, 0.0f, params.debug_bb);
    }
  }

private:
  Orientation<> orientation;
  FastNoiseLite noise;
  Timeline timeline;
  Pipeline<W, H> filters;

  MobiusWarpGnomonicTransformer<1> transformer;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(GnomonicStars)
