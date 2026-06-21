/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include "core/effects_engine.h"

/**
 * @brief Spinning great-circle rings that wander the sphere.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Each ring's orientation follows a random-walk over the sphere and
 * leaves a motion-blur trail that fades in color and alpha along its length.
 * @note Sibling trail effects — `Comets` and `ChaoticStrings` — share the same
 *       scaffolding (orientation random-walk → motion-blur/fading trail) and are
 *       not unified into a shared base, so propagate trail-rendering fixes across
 *       all three. Known divergences from the two siblings:
 *         - This effect does not apply the `Screen::AntiAlias` filter the
 *           siblings use.
 *         - This effect uses `Orientation<>` (CAP 4) for the ring orientation and
 *           trail, where the siblings use `Orientation<16>` — i.e. up to 4
 *           motion-blur sub-frames per recorded frame versus 16. Intentional: a
 *           full great-circle ring is a spatially huge shape whose successive
 *           trail frames overlap almost completely, so 4 interpolated sub-frames
 *           read identically to 16; the siblings smear a single fast-moving
 *           point/string where the sub-frame count sets the blur smoothness.
 */
template <int W, int H> class RingSpin : public Effect {
public:
  static constexpr int TRAIL_LENGTH = 19; // trail samples per ring
  static constexpr int NUM_RINGS = 4;

  /**
   * @brief One ring: great-circle plane, palette, orientation, and trail.
   * @details Bundles the ring's great-circle plane normal, palette, current
   * orientation, and the history trail used to render the fading motion blur.
   */
  struct Ring {
    Vector normal;
    BakedPalette *palette;
    Orientation<> orientation;
    Animation::OrientationTrail<Orientation<>, TRAIL_LENGTH> trail;
    FastNoiseLite noise;
    /**
     * @brief Constructs a default ring on the X axis with no palette.
     */
    Ring() : normal(X_AXIS), palette(nullptr) {}

    /**
     * @brief Constructs a ring with the given plane normal and palette.
     * @param n Unit normal of the ring's great-circle plane.
     * @param p Baked palette used to color the ring's trail.
     */
    Ring(const Vector &n, BakedPalette *p) : normal(n), palette(p) {}
  };

  /**
   * @brief Constructs the effect at the W x H canvas resolution.
   */
  FLASHMEM RingSpin() : Effect(W, H) {}

  /**
   * @brief Allocates rings, registers params, bakes palettes, and spawns rings.
   * @details Allocates the ring storage from the persistent arena, registers the
   * tunable parameters, bakes the vignette palettes into fast LUTs, and spawns
   * the initial set of rings.
   */
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

  /**
   * @brief Reports whether the effect draws over a background.
   * @return false; the effect renders on a clear background.
   */
  bool show_bg() const override { return false; }

  /**
   * @brief Advances the timeline and draws each ring's trail.
   * @details Steps the timeline, then draws each ring's trail back-to-front,
   * fading color and alpha along the trail.
   */
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
  /**
   * @brief Constructs one more ring and starts its energetic random-walk.
   * @param normal Unit normal of the new ring's great-circle plane.
   * @param palette Baked palette used to color the new ring's trail.
   * @details Adds the ring to the timeline with an energetic random-walk
   * orientation; no-op once NUM_RINGS rings are already live.
   */
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

  /**
   * @brief Tunable rendering parameters for the effect.
   */
  struct Params {
    float alpha = 0.5f;          /**< Global trail opacity multiplier in [0, 1]. */
    float thickness = 0.8f;      /**< Ring line thickness multiplier (unitless). */
    bool show_bounding_box = false; /**< Whether to draw each ring's bounding box. */
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(RingSpin)
