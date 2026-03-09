/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"

template <int W, int H> class GnomonicStars : public Effect {
public:
  struct Params {
    float points = 600.0f;
    float star_radius = 0.02f;
    float star_sides = 4.0f;
    bool debug_bb = false;
  } params;

  FLASHMEM GnomonicStars()
      : Effect(W, H), orientation(), timeline(), transformer(timeline) {
    registerParam("Points", &params.points, 100.0f, 2000.0f);
    registerParam("Radius", &params.star_radius, 0.01f, 0.1f);
    registerParam("Sides", &params.star_sides, 3.0f, 8.0f);
    registerParam("Debug BB", &params.debug_bb);

    // Spawn the evolving warp and get its animation reference
    auto *anim = transformer.spawn(0, 0.5f, 0.035f);
    if (anim) {
      registerParam("Warp Speed", &anim->speed, 0.0f, 1.0f);
    }

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, UP, noise,
                        Animation::RandomWalk<W>::Options::Languid()));
  }

  bool show_bg() const override { return true; }

  void draw_frame() override {
    Canvas canvas(*this);

    timeline.step(canvas);

    // Fragment Shader for Stars (Mango Peel Gradient based on Y)
    auto fragment_shader = [](const Vector &p, Fragment &frag) {
      float t = (p.y + 1.0f) * 0.5f; // Map [-1, 1] to [0, 1]
      Color4 c = Palettes::mangoPeel.get(t);
      frag.color = c;
      frag.blend = BLEND_OVER;
    };

    const int points = (int)params.points;
    const float radius = params.star_radius;
    const int sides = (int)params.star_sides;

    for (int i = 0; i < points; i++) {
      Vector v = fib_spiral(points, 0.0f, i);

      // Apply Transformer
      v = transformer.transform(v);

      // Orient the star position (using latest orientation)
      v = orientation.orient(v);

      // Create Basis at the star position (aligned with normal v)
      Basis basis = make_basis(orientation.get(), v);

      Scan::Star::draw<W, H>(filters, canvas, basis, radius, sides,
                             fragment_shader, 0.0f, params.debug_bb);
    }
  }

private:
  Orientation<W> orientation;
  FastNoiseLite noise;
  Timeline<W> timeline;
  Pipeline<W, H> filters;

  MobiusWarpGnomonicTransformer<W, 1> transformer;
};
