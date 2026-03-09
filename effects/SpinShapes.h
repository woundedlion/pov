/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"

template <int W, int H> class SpinShapes : public Effect {
public:
  FLASHMEM SpinShapes()
      : Effect(W, H), filters(Filter::World::Orient<W>(camera),
                              Filter::Screen::AntiAlias<W, H>()) {
    registerParam("Sides", &params.sides, 3.0f, 12.0f);
    registerParam("Radius", &params.radius, 0.05f, 2.0f);
    registerParam("Count", &params.count, 1.0f, 100.0f);
    registerParam("Debug BB", &params.debug_bb);

    rebuild();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
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
    bool debug_bb = false;
  } params;

  StaticCircularBuffer<Shape, 128> shapes;
  Timeline<W, 128> timeline;

  Orientation<W> camera;
  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;

  FLASHMEM void rebuild() {
    shapes.clear();
    timeline = Timeline<W, 128>();

    for (int i = 0; i < (int)params.count; i++) {
      Vector normal = fib_spiral((int)params.count, 0, i);

      shapes.push_back({normal, Orientation<W>(), 0});
      Shape &s = shapes.back();
      timeline.add(0, Animation::Rotation<W>(s.orientation, s.normal, 2 * PI_F,
                                             300 + i * 2, ease_mid, true,
                                             Animation::Space::Local));

      timeline.add(0, Animation::Sprite(
                          [this, i](Canvas &c, float opacity) {
                            this->draw_shape(c, this->shapes[i], opacity);
                          },
                          -1));
    }
  }

  void draw_shape(Canvas &canvas, Shape &shape, float opacity) {
    float t = (shape.normal.y + 1.0f) / 2.0f;
    Color4 c = Palettes::brightSunrise.get(t);
    c.alpha *= opacity * 0.6f;

    auto fragment_shader = [&](const Vector &p, Fragment &f) { f.color = c; };

    Basis basis = make_basis(shape.orientation.get(), shape.normal);
    float phase = 0.0f;

    Scan::Star::draw<W, H>(filters, canvas, basis, params.radius,
                           (int)params.sides, fragment_shader, phase,
                           params.debug_bb);
  }
};
