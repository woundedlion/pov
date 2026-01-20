/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <functional>
#include <memory>
#include <vector>
#include <map>
#include "../effects_engine.h"

template <int W>
class RingShower : public Effect {
public:

  RingShower() :
    Effect(W),
    palette(GradientShape::STRAIGHT, HarmonyType::ANALOGOUS, BrightnessProfile::ASCENDING)
  {
    persist_pixels = false;

    timeline.add(0,
      RandomTimer(4, 48,
        [this](auto&) { this->spawn_ring(); },
        true)
    );

  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

private:

  struct Ring {
    Vector normal;
    float radius;
    float duration;

    Ring() :
      normal(random_vector()),
      radius(0.0),
      duration(0.0)
    {
    }
  };

  void spawn_ring() {
    for (size_t i = 0; i < MAX_RINGS; ++i) {
      if (rings[i].duration <= 0) {
        Ring& ring = rings[i];
        ring.normal = random_vector();
        ring.duration = hs::rand_f() * 72.0f + 8.0f;
        ring.radius = 0;

        timeline.add(0,
          Sprite(
            [this, i](Canvas& canvas, float opacity) { this->draw_ring(canvas, opacity, i); },
            static_cast<int>(ring.duration),
            4, ease_mid,
            0, ease_mid)
        );

        timeline.add(0,
          Transition(ring.radius, 2.0f, static_cast<int>(ring.duration), ease_mid, false, false)
          .then([&ring]() { ring.duration = 0; })
        );

        return;
      }
    }
  }

  void draw_ring(Canvas& canvas, float opacity, size_t index) {
    Ring& ring = rings[index];
    Basis basis = make_basis(orientation.get(), ring.normal);
    Plot<W>::Ring::draw(filters, canvas, basis, ring.radius,
      [&](const Vector& v, float t) {
        Vector z = orientation.orient(X_AXIS);
        Color4 c = palette.get(angle_between(z, v) / PI_F);
        c.alpha *= opacity * alpha;
        return c;
      },
      0);
  }


  static constexpr size_t MAX_RINGS = 16;
  Ring rings[MAX_RINGS];
  Pipeline<W, FilterAntiAlias<W>> filters;
  Orientation orientation;
  Timeline timeline;

  GenerativePalette palette;
  static constexpr float alpha = 0.2f;
};