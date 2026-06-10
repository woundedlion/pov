/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include <functional>
#include "core/effects_engine.h"

template <int W, int H> class RingShower : public Effect {
public:
  FLASHMEM RingShower() : Effect(W, H) {}

  void init() override {
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);

    timeline.add(0, Animation::RandomTimer(
                        4, 48, [this](auto &) { this->spawn_ring(); }, true));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas); // drives the spawn timer only

    // Advance and draw each live ring directly from its slot. Radius growth,
    // fade-in and lifetime are pure functions of `age` computed here rather
    // than by per-ring Sprite/Transition animations capturing the slot: a slot
    // is recyclable, so an animation outliving its reuse would draw whatever
    // ring later lands in it. Driving everything from `age` makes the slot's
    // lifetime explicit and self-contained (mirrors Thrusters).
    for (size_t i = 0; i < MAX_RINGS; ++i) {
      Ring &ring = rings[i];
      if (ring.age >= ring.life)
        continue; // free slot
      draw_ring(canvas, ring.opacity_at(), i);
      ++ring.age;
    }
  }

private:
  struct Ring {
    // Frames spent fading in from 0 to full opacity (was the Sprite's
    // fade_in=4; there was no fade-out).
    static constexpr int FADE_IN_FRAMES = 4;
    // Ring radius grows from 0 to RADIUS_MAX over its whole life.
    static constexpr float RADIUS_MAX = 2.0f;

    Vector normal;
    GenerativePalette palette;
    int life = 0; // total visible frames
    int age = 0;  // frames elapsed; age >= life marks the slot free

    Ring() : normal(random_vector()) {}

    // Eased radius for the current age, mirroring the prior Transition. The
    // Transition stepped once before the first draw, so age + 1 reproduces its
    // value sequence (the whole ring is shifted one frame earlier overall).
    float radius_at() const {
      float t = hs::clamp(static_cast<float>(age + 1) / life, 0.0f, 1.0f);
      return RADIUS_MAX * ease_mid(t);
    }

    // Fade-in over the first FADE_IN_FRAMES frames, then hold full; no fade-out.
    float opacity_at() const {
      if (age + 1 < FADE_IN_FRAMES)
        return ease_mid(static_cast<float>(age + 1) / FADE_IN_FRAMES);
      return 1.0f;
    }
  };

  void spawn_ring() {
    for (size_t i = 0; i < MAX_RINGS; ++i) {
      if (rings[i].age >= rings[i].life) { // free slot
        Ring &ring = rings[i];
        ring.normal = random_vector();
        ring.life = static_cast<int>(hs::rand_f() * 72.0f + 8.0f);
        ring.age = 0;
        ring.palette =
            GenerativePalette(GradientShape::CIRCULAR, HarmonyType::ANALOGOUS,
                              BrightnessProfile::FLAT);
        return;
      }
    }
  }

  void draw_ring(Canvas &canvas, float opacity, size_t index) {
    Ring &ring = rings[index];
    Basis basis = make_basis(orientation.get(), ring.normal);
    auto fragment_shader = [&](const Vector &v, Fragment &f) {
      Vector z = orientation.orient(X_AXIS);
      f.color = ring.palette.get(angle_between(z, v) / PI_F);
      f.color.alpha *= opacity * params.alpha;
    };
    Plot::Ring::draw<W, H>(filters, canvas, basis, ring.radius_at(),
                           fragment_shader);
  }

  static constexpr size_t MAX_RINGS = 16;
  Ring rings[MAX_RINGS];
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  Orientation<W> orientation;
  Timeline<W> timeline;

  struct Params {
    float alpha = 0.2f;
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(RingShower)
