/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include "core/effects_engine.h"

template <int W, int H> class RingSpin : public Effect {
public:
  static constexpr int TRAIL_LENGTH = 19;
  static constexpr int NUM_RINGS = 4;

  struct Ring {
    Vector normal;
    BakedPalette *palette;
    Orientation<W> orientation;
    Animation::OrientationTrail<Orientation<W>, TRAIL_LENGTH> trail;
    FastNoiseLite noise;
    Ring() : normal(X_AXIS), palette(nullptr) {}

    Ring(const Vector &n, BakedPalette *p) : normal(n), palette(p) {}
  };

  FLASHMEM RingSpin() : Effect(W, H) {}

  void init() override {
    rings = static_cast<Ring *>(
        persistent_arena.allocate(NUM_RINGS * sizeof(Ring), alignof(Ring)));

    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Thickness", &params.thickness, 0.1f, 10.0f);
    registerParam("Show Bounding", &params.show_bounding_box);

    // Bake vignette-wrapped palettes into fast LUTs
    for (int i = 0; i < (int)source_palettes.size(); ++i) {
      TransparentVignette v(&source_palettes[i]);
      baked_palettes[i].bake(persistent_arena, v);
    }

    for (int i = 0; i < NUM_RINGS; ++i) {
      int p_idx = i % source_palettes.size();
      spawn_ring(Y_AXIS, &baked_palettes[p_idx]);
    }
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    int draw_calls = 0;
    int total_pixels = 0;
    unsigned long basis_us = 0;
    unsigned long raster_us = 0;

    for (int i = 0; i < num_rings; ++i) {
      Ring &ring = rings[i];
      ring.trail.record(ring.orientation);
      deep_tween(ring.trail, [&](const Quaternion &q, float t) {
        draw_calls++;
        Color4 c = ring.palette->get(1.0f - t);
        c.alpha = c.alpha * params.alpha;

        unsigned long t0 = micros();
        Basis basis = make_basis(q, ring.normal);
        basis_us += micros() - t0;

        auto fragment_shader = [&](const Vector &, Fragment &f) {
          total_pixels++;
          f.color = c;
        };

        // Adaptive thickness: head/tail of trail = 2px, intermediate = 1px
        constexpr float pixel_w = 2.0f * PI_F / W;
        float th = (t < 0.01f || t > 0.95f) ? 2.0f * pixel_w : 1.0f * pixel_w;

        unsigned long t1 = micros();
        Scan::Ring::draw<W, H, false>(filters, canvas, basis, 1.0f,
                                      th, fragment_shader, 0.0f,
                                      params.show_bounding_box);
        raster_us += micros() - t1;
      });
    }
    Serial.print("draws:");
    Serial.print(draw_calls);
    Serial.print(" px:");
    Serial.print(total_pixels);
    Serial.print(" basis:");
    Serial.print(basis_us);
    Serial.print("us rast:");
    Serial.print(raster_us);
    Serial.println("us");
  }

private:
  void spawn_ring(const Vector &normal, BakedPalette *palette) {
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
  std::array<BakedPalette, 4> baked_palettes;

  struct Params {
    float alpha = 0.5f;
    float thickness = 1.0f * (2 * PI_F / W);
    bool show_bounding_box = false;
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(RingSpin)
