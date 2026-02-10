/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include "../effects_engine.h"

#include "../palettes.h"

template <int W>
class TestShapes : public Effect {
public:
  enum class ShapeType { Polygon, Star, Flower };
  enum class RenderMode { Plot, Scan };

  struct Ring {
    Vector normal;
    float scale;
    Color4 color;
    RenderMode mode;
    int layer_index;
    Orientation orientation;

    Ring() : scale(1.0f), color(Color4(0,0,0,0)), mode(RenderMode::Plot), layer_index(0) {}

    Ring(const Vector& n, float s, const Color4& c, RenderMode m, int l)
      : normal(n), scale(s), color(c), mode(m), layer_index(l) {}
  };

  TestShapes() :
    Effect(W),
    alpha(0.5f),
    radius(1.0f),
    sides(5),
    twist(0.0f),
    current_shape(ShapeType::Polygon),
    num_shapes(25),
    debug_bb(false)
  {
    rebuild();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    
    // Cycle shapes every 500 frames (~5-10s depending on FPS)
    static int frame_count = 0;
    frame_count++;
    if (frame_count % 500 == 0) {
      int next = (static_cast<int>(current_shape) + 1) % 3;
      current_shape = static_cast<ShapeType>(next);
    }

    timeline.step(canvas);
  }

  void rebuild() {
    rings.clear();
    // Re-create timeline to clear old animations
    timeline = Timeline(); 

    // Twist Mutation: sin wave
    // In JS: (Math.PI / 4) * Math.sin(t * Math.PI), duration 480
    // C++ Mutation takes a reference.
    timeline.add(0, Animation::Mutation(twist, [](float t) {
      return (PI_F / 4.0f) * sinf(t * PI_F); 
    }, 480, ease_mid, true));

    int total_shapes = num_shapes;
    int seed1 = hs::rand_int(0, 65535);

    for (int i = total_shapes - 1; i >= 0; --i) {
      float t = static_cast<float>(i) / (total_shapes > 1 ? total_shapes - 1 : 1);
      Color4 color = Palettes::richSunset.get(t);
      spawnRing(X_AXIS, t, color, seed1, RenderMode::Plot, i);
      spawnRing(-X_AXIS, t, color, seed1, RenderMode::Scan, i);
    }
  }

  void spawnRing(const Vector& normal, float scale, const Color4& color, int seed, RenderMode mode, int layer_index) {
    if (rings.is_full()) return;
    rings.emplace_back(normal, scale, color, mode, layer_index);
    Ring& ring = rings.back();

    Vector antipode = (normal.i < -0.5f) ? -normal : normal;
    
    // Animations
    typename Animation::RandomWalk<W>::Options rw_opts;
    rw_opts.seed = seed;
    rw_opts.space = Animation::Space::World;
    timeline.add(0, Animation::RandomWalk<W>(ring.orientation, antipode, rw_opts));
    timeline.add(0, Animation::Rotation<W>(ring.orientation, ring.normal, 2 * PI_F, 160, ease_mid, true, Animation::Space::Local));
    timeline.add(0, Animation::Sprite(
      [this, &ring](Canvas& canvas, float opacity) {
        this->drawShape(canvas, ring, opacity);
      },
      -1, 0, ease_mid, 0, ease_mid));
  }

  void drawShape(Canvas& canvas, const Ring& ring, float sprite_alpha) {
    auto fragment_shader = [&](const Vector& p, const Fragment& f) -> Fragment {
      Color4 c = ring.color;
      c.alpha = c.alpha * this->alpha * sprite_alpha;
      Fragment f_out = f;
      f_out.color = c;
      return f_out;
    };

    float phase = ring.layer_index * this->twist;
    float r = this->radius * ring.scale;
    
    Basis basis = make_basis(ring.orientation.get(), ring.normal);

    if (ring.mode == RenderMode::Plot) {
      switch (current_shape) {
        case ShapeType::Flower:
           Plot::Flower::draw<W>(plot_filters, canvas, basis, r, this->sides, fragment_shader, phase);
           break;
        case ShapeType::Star:
           Plot::Star::draw<W>(plot_filters, canvas, basis, r, this->sides, fragment_shader, phase);
           break;
        default: // Polygon
           Plot::SphericalPolygon::draw<W>(plot_filters, canvas, basis, r, this->sides, fragment_shader, phase);
           break;
      }
    } else {
      switch (current_shape) {
        case ShapeType::Flower:
           Scan::Flower::draw<W>(scan_filters, canvas, basis, r, this->sides, fragment_shader, phase, debug_bb);
           break;
        case ShapeType::Star:
           Scan::Star::draw<W>(scan_filters, canvas, basis, r, this->sides, fragment_shader, phase, debug_bb);
           break;
        default:
           Scan::SphericalPolygon::draw<W>(scan_filters, canvas, basis, r, this->sides, fragment_shader, phase, debug_bb);
           break;
      }
    }
  }

private:
  Timeline timeline;
  StaticCircularBuffer<Ring, 128> rings;
  Pipeline<W, Filter::Screen::AntiAlias<W>> plot_filters;
  Pipeline<W> scan_filters;

  float alpha;
  float radius;
  int sides;
  float twist;
  ShapeType current_shape;
  int num_shapes;
  bool debug_bb;
};
