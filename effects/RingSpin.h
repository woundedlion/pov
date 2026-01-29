/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <vector>
#include <array>
#include "../effects_engine.h"

template <int W>
class RingSpin : public Effect {
public:
  static constexpr int TRAIL_LENGTH = 19;
  static constexpr int NUM_RINGS = 4;

  struct Ring {
    Vector normal;
    TransparentVignette palette;
    Orientation orientation;
    OrientationTrail<TRAIL_LENGTH> trail;
    Ring() : normal(X_AXIS), palette(iceMelt) {}

    Ring(const Vector& n, const Palette& p) :
      normal(n), palette(p) {
    }
  };

  RingSpin() :
    Effect(W),
    alpha(0.5f),
    thickness(2.0f * (2 * PI_F / W))
  {
    persist_pixels = false;

    std::vector<const Palette*> source_palettes = {
        &iceMelt,
        &undersea,
        &mangoPeel,
        &richSunset
    };

    for (int i = 0; i < NUM_RINGS; ++i) {
      int p_idx = i % source_palettes.size();
      spawn_ring(X_AXIS, *source_palettes[p_idx]);
    }
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    for (auto& ring : rings) {
      ring.trail.record(ring.orientation);
      deep_tween(ring.trail, [&](const Quaternion& q, float t) {
          Color4 c = ring.palette.get(t);
          c.alpha = c.alpha * alpha;
          Basis basis = make_basis(q, ring.normal);
          Scan::Ring::draw<W>(filters, canvas, basis, 1.0f, thickness, [&](const Vector&, float) { return c; });
        });
    }
  }

private:
  void spawn_ring(const Vector& normal, const Palette& palette) {
    if (rings.is_full()) return;
    rings.push_back(Ring(normal, palette));
    Ring& r = rings.back();
    timeline.add(0, RandomWalk<W>(r.orientation, r.normal));
  }

  Timeline timeline;
  Pipeline<W, FilterAntiAlias<W>> filters;
  StaticCircularBuffer<Ring, NUM_RINGS> rings;

  float alpha;
  float thickness;
};