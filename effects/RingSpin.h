/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include "core/engine.h"

// Unit-test accessor for the live ring-count pool-bound invariant.
namespace hs_test {
namespace effects_tests {
struct RingSpinWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Spinning great-circle rings that wander the sphere.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Each ring's orientation follows a random-walk over the sphere and
 * leaves a motion-blur trail that fades in color and alpha along its length.
 * @note Sibling trail effects `Comets` and `ChaoticStrings` share the same
 *       scaffolding (orientation random-walk → motion-blur/fading trail), but
 *       only the common `OrientationTrail::record` + `deep_tween` skeleton lives
 *       in the engine; the per-effect bodies diverge in draw primitive,
 *       transform, color/fade, and accumulate-vs-draw model, so each renders
 *       independently — propagate trail-rendering fixes by hand across all three.
 *       Two intentional differences from the siblings:
 *         - No `Screen::AntiAlias` filter.
 *         - `Orientation<>` (CAP 4) instead of `Orientation<16>`: a full
 *           great-circle ring is spatially huge and its successive trail frames
 *           overlap almost completely, so 4 sub-frames read identically to 16,
 *           where the siblings smear a fast-moving point/string.
 */
template <int W, int H> class RingSpin : public Effect {
public:
  static constexpr int TRAIL_LENGTH = 19; // trail samples per ring
  static constexpr int NUM_RINGS = 4;
  static constexpr int NUM_PALETTES = 4;

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

    // Inset the source into the middle band, then fade alpha at the edges.
    // Wrap=false so the top edge resolves to the source's last stop
    // (wrap_t(1)==0 would fold it to black).
    InsetModifier inset;
    EdgeAlphaShade edge_fade;
    const ProceduralPalette sources[NUM_PALETTES] = {
        Palettes::iceMelt, Palettes::undersea, Palettes::mangoPeel,
        Palettes::richSunset};
    for (int i = 0; i < NUM_PALETTES; ++i) {
      StaticPalette<ProceduralPalette, Coords<InsetModifier>,
                    Colors<EdgeAlphaShade>, /*Wrap=*/false>
          v;
      v.bind(&sources[i], &inset, &edge_fade);
      baked_palettes[i].bake(persistent_arena, v);
    }

    for (int i = 0; i < NUM_RINGS; ++i) {
      int p_idx = i % NUM_PALETTES;
      spawn_ring(&baked_palettes[p_idx]);
    }
  }

  /// POV column-strobe flag; strobes (see Effect::strobe_columns).
  bool strobe_columns() const override { return true; }

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
        // The trail's length-fade comes entirely from the palette: sampling at
        // 1-t walks a transparent-vignette palette whose alpha tapers to ~0
        // toward the tail, so there is deliberately no explicit t-fade here.
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
  // Test seam: reaches the live ring-count pool-bound invariant.
  friend struct ::hs_test::effects_tests::RingSpinWhiteBox;

  /**
   * @brief Constructs one more ring and starts its energetic random-walk.
   * @param palette Baked palette used to color the new ring's trail.
   * @details Adds the ring to the timeline with an energetic random-walk
   * orientation; no-op once NUM_RINGS rings are already live. Append-one
   * primitive: `num_rings` tracks the live count and the `>= NUM_RINGS` guard
   * bounds it against the pool. init() spawns all NUM_RINGS up front so the
   * guard never fires today — it exists to support runtime spawning (a future
   * slider/event) without overrunning `rings`.
   */
  void spawn_ring(BakedPalette *palette) {
    if (num_rings >= NUM_RINGS)
      return;
    Ring &r = rings[num_rings];
    new (&r) Ring(Y_AXIS, palette);
    num_rings++;
    timeline.add(0, Animation::RandomWalk<W>(
                        r.orientation, r.normal, r.noise,
                        Animation::RandomWalk<W>::Options::Energetic()));
  }

  Timeline timeline;
  Pipeline<W, H> filters;
  Ring *rings = nullptr;
  int num_rings = 0;

  std::array<BakedPalette, NUM_PALETTES> baked_palettes;

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
