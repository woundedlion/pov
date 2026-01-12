/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "../effects_engine.h"

template <int W>
class TestPlotPolygon : public Effect {
public:
  struct Ring {
    Vector normal;
    TransparentVignette palette;
    Orientation orientation;

    Ring(const Vector& n, const Palette& p) : normal(n), palette(p) {}
  };

  TestPlotPolygon() :
    Effect(W),
    alpha(0.5f),
    radius(1.0f),
    sides(6),
    phase(0.0f)
  {
    const Palette* palettes[] = { &iceMelt, &undersea, &mangoPeel, &richSunset };
    int num_palettes = 4;
    int num_rings = 1;

    for (int i = 0; i < num_rings; ++i) {
      spawnRing(X_AXIS, *palettes[i % num_palettes]);
    }
  }

  void spawnRing(const Vector& normal, const Palette& palette) {
    auto ring = std::make_unique<Ring>(normal, palette);
    timeline.add(0, RandomWalk<W>(ring->orientation, ring->normal));
    rings.push_back(std::move(ring));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    for (const auto& ring : rings) {
      auto color_fn = [&](const Vector& p, float t) {
        Color4 c = ring->palette.get(t); // Use t for gradient
        c.alpha *= alpha;
        return c;
      };

      Plot<W>::Polygon::draw(filters, canvas, ring->orientation.get(), ring->normal,
        radius, sides, color_fn, phase);
    }
  }

private:
  std::vector<std::unique_ptr<Ring>> rings;
  Timeline timeline;
  Pipeline<W, FilterAntiAlias<W>> filters;
  float alpha;
  float radius;
  int sides;
  float phase;
};
