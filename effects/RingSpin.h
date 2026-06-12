/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include "core/effects_engine.h"

// Spinning great-circle rings that wander the sphere via random-walk
// orientation, each leaving a motion-blur trail that fades along its length.
template <int W, int H> class RingSpin : public Effect {
public:
  static constexpr int TRAIL_LENGTH = 19; // trail samples per ring
  static constexpr int NUM_RINGS = 4;

  // One ring: its great-circle plane, palette, current orientation, and the
  // history trail used to render the fading motion blur.
  struct Ring {
    Vector normal;
    BakedPalette *palette;
    Orientation<> orientation;
    Animation::OrientationTrail<Orientation<>, TRAIL_LENGTH> trail;
    FastNoiseLite noise;
    Ring() : normal(X_AXIS), palette(nullptr) {}

    Ring(const Vector &n, BakedPalette *p) : normal(n), palette(p) {}
  };

  FLASHMEM RingSpin() : Effect(W, H) {}

  // Allocate rings, register params, bake the vignette palettes, and spawn the
  // initial set of rings.
  void init() override {
    rings = static_cast<Ring *>(
        persistent_arena.allocate(NUM_RINGS * sizeof(Ring), alignof(Ring)));

    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Thickness", &params.thickness, 0.01f, 10.0f);
    registerParam("Show Bounding", &params.show_bounding_box);

    // Bake transparent-vignette palettes into fast LUTs: inset the source into
    // the middle band, then fade alpha at the edges. Wrap=false so the top edge
    // resolves to the source's last stop (wrap_t(1)==0 would fold it to black).
    InsetModifier inset;
    EdgeAlphaShade edge_fade;
    for (int i = 0; i < (int)source_palettes.size(); ++i) {
      StaticPalette<ProceduralPalette, Coords<InsetModifier>,
                    Colors<EdgeAlphaShade>, /*Wrap=*/false>
          v;
      v.bind(&source_palettes[i], &inset, &edge_fade);
      baked_palettes[i].bake(persistent_arena, v);
    }

    for (int i = 0; i < NUM_RINGS; ++i) {
      int p_idx = i % source_palettes.size();
      spawn_ring(Y_AXIS, &baked_palettes[p_idx]);
    }
  }

  bool show_bg() const override { return false; }

  // Advance the timeline and draw each ring's trail back-to-front, fading color
  // and alpha along the trail.
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    for (int i = 0; i < num_rings; ++i) {
      Ring &ring = rings[i];
      ring.trail.record(ring.orientation);
      deep_tween(ring.trail, [&](const Quaternion &q, float t) {
        Color4 c = ring.palette->get(1.0f - t);
        c.alpha = c.alpha * params.alpha;
        if (c.alpha <= 0.001f) return;

        Basis basis = make_basis(q, ring.normal);

        auto fragment_shader = [&](const Vector &, Fragment &f) {
          f.color = c;
        };

        // Adaptive thickness: head/tail of trail = 2px, intermediate = 1px
        constexpr float pixel_w = 2.0f * PI_F / W;
        float th = ((t < 0.01f || t > 0.95f) ? 2.0f * pixel_w : 1.0f * pixel_w) * params.thickness;

        Scan::Ring::draw<W, H, false>(filters, canvas, basis, 1.0f, th,
                                      fragment_shader, 0.0f,
                                      params.show_bounding_box);
      });
    }
  }

private:
  // Construct one more ring and start its energetic random-walk on the timeline;
  // no-op once NUM_RINGS are live.
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

  Timeline timeline;
  Pipeline<W, H> filters;
  Ring *rings = nullptr;
  int num_rings = 0;

  std::array<ProceduralPalette, 4> source_palettes = {
      Palettes::iceMelt, Palettes::undersea, Palettes::mangoPeel,
      Palettes::richSunset};
  std::array<BakedPalette, 4> baked_palettes;

  struct Params {
    float alpha = 0.5f;
    float thickness = 0.8f;
    bool show_bounding_box = false;
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(RingSpin)
