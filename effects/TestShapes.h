/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"

template <int W, int H> class TestShapes : public Effect {
public:
  enum class ShapeType { PlanarPolygon, SphericalPolygon, Flower, Star };
  enum class RenderMode { Plot, Scan };

  struct Ring {
    Vector normal;
    float scale;
    Color4 color;
    RenderMode mode;
    int layer_index;
    Orientation<W> orientation;

    Ring()
        : scale(1.0f), color(Color4(0, 0, 0, 0)), mode(RenderMode::Plot),
          layer_index(0) {}

    Ring(const Vector &n, float s, const Color4 &c, RenderMode m, int l)
        : normal(n), scale(s), color(c), mode(m), layer_index(l) {}
  };

  TestShapes()
      : Effect(W, H), current_shape(ShapeType::PlanarPolygon), num_shapes(25),
        debug_bb(false) {
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Radius", &params.radius, 0.1f, 5.0f);
    registerParam("Sides", &params.sides, 3.0f, 12.0f);
    registerParam("Twist", &params.twist, -5.0f, 5.0f);

    rebuild();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);

    // Cycle shapes every n frames
    static int frame_count = 0;
    frame_count++;
    if (frame_count % 500 == 0) {
      int next = (static_cast<int>(current_shape) + 1) % 4;
      current_shape = static_cast<ShapeType>(next);
    }

    timeline.step(canvas);
  }

  void rebuild() {
    rings.clear();
    // Re-create timeline to clear old animations
    timeline = Timeline<W, 256>();

    // Twist Mutation: sin wave
    timeline.add(0, Animation::Mutation(
                        params.twist,
                        [](float t) { return (PI_F / 4.0f) * sinf(t * PI_F); },
                        480, ease_mid, true));

    int total_shapes = num_shapes;
    int seed1 = hs::rand_int(0, 65535);

    for (int i = total_shapes - 1; i >= 0; --i) {
      float t =
          static_cast<float>(i) / (total_shapes > 1 ? total_shapes - 1 : 1);
      Color4 color = Palettes::richSunset.get(t);
      spawnRing(X_AXIS, t, color, seed1, RenderMode::Plot, i);
      spawnRing(-X_AXIS, t, color, seed1, RenderMode::Scan, i);
    }
  }

  void spawnRing(const Vector &normal, float scale, const Color4 &color,
                 int seed, RenderMode mode, int layer_index) {
    if (rings.is_full())
      return;
    rings.emplace_back(normal, scale, color, mode, layer_index);
    Ring &ring = rings.back();

    Vector antipode = (normal.i < -0.5f) ? -normal : normal;

    // Animations
    typename Animation::RandomWalk<W>::Options rw_opts;
    rw_opts.seed = seed;
    rw_opts.space = Animation::Space::World;
    timeline.add(0,
                 Animation::RandomWalk<W>(ring.orientation, antipode, rw_opts));
    timeline.add(0, Animation::Rotation<W>(ring.orientation, ring.normal,
                                           2 * PI_F, 160, ease_mid, true,
                                           Animation::Space::Local));
    timeline.add(0, Animation::Sprite(
                        [this, &ring](Canvas &canvas, float opacity) {
                          this->drawShape(canvas, ring, opacity);
                        },
                        -1, 0, ease_mid, 0, ease_mid));
  }

  void drawShape(Canvas &canvas, const Ring &ring, float sprite_alpha) {
    auto fragment_shader = [&](const Vector &p, Fragment &f) {
      Color4 c = ring.color;
      c.alpha = c.alpha * this->params.alpha * sprite_alpha;
      f.color = c;
    };

    float phase = ring.layer_index * this->params.twist;
    float r = this->params.radius * ring.scale;

    Basis basis = make_basis(ring.orientation.get(), ring.normal);
    int sides_int = (int)params.sides;
    if (ring.mode == RenderMode::Plot) {
      switch (current_shape) {
      case ShapeType::Flower:
        Plot::Flower::draw<W, H>(plot_filters, canvas, basis, r, sides_int,
                                 fragment_shader, NullVertexShader{}, phase);
        break;
      case ShapeType::Star:
        Plot::Star::draw<W, H>(plot_filters, canvas, basis, r, sides_int,
                               fragment_shader, phase);
        break;
      case ShapeType::PlanarPolygon:
        Plot::PlanarPolygon::draw<W, H>(plot_filters, canvas, basis, r,
                                        sides_int, fragment_shader, phase);
        break;
      default: // SphericalPolygon
        Plot::SphericalPolygon::draw<W, H>(plot_filters, canvas, basis, r,
                                           sides_int, fragment_shader, phase);
        break;
      }
    } else {
      switch (current_shape) {
      case ShapeType::Flower:
        debug_bb = true;
        Scan::Flower::draw<W, H>(scan_filters, canvas, basis, r, sides_int,
                                 fragment_shader, phase, debug_bb);
        break;
      case ShapeType::Star:
        Scan::Star::draw<W, H>(scan_filters, canvas, basis, r, sides_int,
                               fragment_shader, phase, debug_bb);
        break;
      case ShapeType::PlanarPolygon:
        Scan::PlanarPolygon::draw<W, H>(scan_filters, canvas, basis, r,
                                        sides_int, fragment_shader, phase,
                                        debug_bb);
        break;
      default: // SphericalPolygon
        Scan::SphericalPolygon::draw<W, H>(scan_filters, canvas, basis, r,
                                           sides_int, fragment_shader, phase,
                                           debug_bb);
        break;
      }
    }
  }

private:
  Timeline<W, 256> timeline;
  StaticCircularBuffer<Ring, 128> rings;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> plot_filters;
  Pipeline<W, H> scan_filters;

  struct Params {
    float alpha = 0.5f;
    float radius = 1.0f;
    float sides = 5.0f;
    float twist = 0.0f;
  } params;

  ShapeType current_shape;
  int num_shapes;
  bool debug_bb;
};
