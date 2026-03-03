/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include <vector>

template <int W, int H> class SpinShapes : public Effect {
public:
  FLASHMEM SpinShapes()
      : Effect(W, H), filters(Filter::World::Orient<W>(camera),
                              Filter::Screen::AntiAlias<W, H>()) {
    registerParam("Sides", &params.sides, 3.0f, 12.0f);
    registerParam("Radius", &params.radius, 0.05f, 2.0f);
    registerParam("Count", &params.count, 1.0f, 100.0f);

    rebuild();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);

    // Step animations (manual management to avoid Timeline limits)
    for (auto &rot : rotations) {
      rot.step(canvas);
    }

    // Draw shapes
    for (auto &shape : shapes) {
      draw_shape(canvas, shape);
    }
  }

private:
  struct Shape {
    Vector normal;
    Orientation<W> orientation;
    int layer;
  };

  struct Params {
    float sides = 3.0f;
    float radius = 0.2f;
    float count = 40.0f;
  } params;

  std::vector<Shape> shapes;
  std::vector<Animation::Rotation<W>> rotations;

  Orientation<W> camera;
  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;

  void rebuild() {
    shapes.clear();
    rotations.clear();

    // Ensure stability for references
    shapes.reserve((int)params.count);
    rotations.reserve((int)params.count);

    for (int i = 0; i < (int)params.count; i++) {
      Vector normal = fib_spiral((int)params.count, 0, i);

      shapes.push_back({normal, Orientation<W>(), 0});
      Shape &s = shapes.back();
      // Rotation(orientation, axis, angle, duration, easing, repeat, space)
      rotations.emplace_back(s.orientation, s.normal, 2 * PI_F, 300 + i * 2,
                             ease_mid, true, Animation::Space::Local);
    }
  }

  void draw_shape(Canvas &canvas, Shape &shape) {
    float t = (shape.normal.j + 1.0f) / 2.0f;
    auto c = Palettes::brightSunrise.get(t);

    auto fragment_shader = [&](const Vector &p, Fragment &f) {
      f.color = Color4(c, 0.6f); // alpha fixed
    };

    Basis basis = make_basis(shape.orientation.get(), shape.normal);
    float phase = 0.0f;

    Scan::Star::draw<W, H>(filters, canvas, basis, params.radius,
                           (int)params.sides, fragment_shader, phase);
  }
};
