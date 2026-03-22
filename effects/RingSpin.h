/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include "effects_engine.h"

template <int W, int H> class RingSpin : public Effect {
public:
  static constexpr int TRAIL_LENGTH = 19;
  static constexpr int NUM_RINGS = 4;

  struct Ring {
    Vector normal;
    TransparentVignette palette;
    Orientation<W> orientation;
    Animation::OrientationTrail<Orientation<W>, TRAIL_LENGTH> trail;
    FastNoiseLite noise;
    Ring() : normal(X_AXIS), palette(nullptr) {}

    Ring(const Vector &n, const Palette *p) : normal(n), palette(p) {}
  };

  FLASHMEM RingSpin() : Effect(W, H) {}

  void init() override {
    rings = static_cast<Ring *>(
        persistent_arena.allocate(NUM_RINGS * sizeof(Ring), alignof(Ring)));

    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Thickness", &params.thickness, 0.1f, 10.0f);
    registerParam("Show Bounding", &params.show_bounding_box);

    for (int i = 0; i < NUM_RINGS; ++i) {
      int p_idx = i % source_palettes.size();
      spawn_ring(X_AXIS, &source_palettes[p_idx]);
    }
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    for (int i = 0; i < num_rings; ++i) {
      Ring &ring = rings[i];
      ring.trail.record(ring.orientation);
      deep_tween(ring.trail, [&](const Quaternion &q, float t) {
        Color4 c = ring.palette.get(1.0f - t);
        c.alpha = c.alpha * params.alpha;
        Basis basis = make_basis(q, ring.normal);
        auto fragment_shader = [&](const Vector &, Fragment &f) {
          f.color = c;
        };
        Scan::Ring::draw<W, H, false>(filters, canvas, basis, 1.0f,
                                      params.thickness, fragment_shader, 0.0f,
                                      params.show_bounding_box);
      });
    }
  }

private:
  void spawn_ring(const Vector &normal, const Palette *palette) {
    if (num_rings >= NUM_RINGS)
      return;
    Ring &r = rings[num_rings];
    new (&r) Ring(normal, palette);
    num_rings++;
    timeline.add(0, Animation::RandomWalk<W>(
                        r.orientation, r.normal, r.noise,
                        Animation::RandomWalk<W>::Options::Energetic()));
  }

  Timeline<W> timeline;
  Pipeline<W, H> filters;
  Ring *rings = nullptr;
  int num_rings = 0;

  std::array<ProceduralPalette, 4> source_palettes = {
      Palettes::iceMelt, Palettes::undersea, Palettes::mangoPeel,
      Palettes::richSunset};

  struct Params {
    float alpha = 0.5f;
    float thickness = 2.0f * (2 * PI_F / W);
    bool show_bounding_box = false;
  } params;
};

#include "effect_registry.h"
REGISTER_EFFECT(RingSpin)
