/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file RingSpin.h
 * @brief Spinning great-circle rings that wander the sphere behind motion-blur
 *        trails.
 */

#include <array>
#include <new> // std::launder
#include "core/engine/engine.h"

// Unit-test accessor for the ring pool's orientations and stroke geometry.
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
 * @note Sibling trail effects `Comets` and `Fishbowl` share the
 *       record + deep_tween skeleton; draw primitive, transform chain and
 *       colour/fade are hand-propagated. Ring carries a plane normal, palette
 *       and noise alongside the orientation + trail, so it does not use their
 *       `Animation::TrailBody`. Differences here: no `Screen::AntiAlias`, and
 *       `Orientation<>` (CAP 4) not `Orientation<16>` — a great-circle ring's
 *       successive trail frames overlap almost completely, so 4 sub-frames read
 *       identically to 16.
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
     * @brief Constructs a ring on the Y-axis great-circle plane.
     * @param p Baked palette used to color the ring's trail.
     */
    Ring(BakedPalette *p) : normal(Y_AXIS), palette(p) {}
  };

  /**
   * @brief Constructs the effect at the W x H canvas resolution.
   */
  HS_COLD_MEMBER RingSpin()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})) {}

  /**
   * @brief Allocates rings, registers params, bakes palettes, and spawns rings.
   * @details Allocates the ring storage from the persistent arena, registers the
   * tunable parameters, bakes the vignette palettes into fast LUTs, and starts
   * each ring's energetic random-walk.
   */
  void init() override {
    register_param("Alpha", &params.alpha, 0.0f, 1.0f);
    register_param("Thickness", &params.thickness, 0.01f, 10.0f);
    register_param("Debug BB", &params.debug_bb);

    // Inset the source into the middle band, then fade alpha at the edges.
    // Wrap=false so the top edge resolves to the source's last stop
    // (wrap_t(1)==0 would fold it to black).
    InsetModifier inset;
    EdgeAlphaShade edge_fade;
    const ProceduralPalette sources[NUM_PALETTES] = {
        Palettes::ICE_MELT, Palettes::UNDERSEA, Palettes::MANGO_PEEL,
        Palettes::RICH_SUNSET};
    for (int i = 0; i < NUM_PALETTES; ++i) {
      StaticPalette<ProceduralPalette, Coords<InsetModifier>,
                    Colors<EdgeAlphaShade>, /*Wrap=*/false>
          v;
      v.bind(&sources[i], &inset, &edge_fade);
      baked_palettes[i].bake(persistent_arena, v);
    }

    rings = persistent_arena.make_n_indexed<Ring>(NUM_RINGS, [&](size_t i) {
      return Ring(&baked_palettes[i % NUM_PALETTES]);
    });
    for (int i = 0; i < NUM_RINGS; ++i) {
      Ring &r = rings[i];
      timeline.add(0, Animation::RandomWalk<W>(
                          r.orientation, r.normal, r.noise,
                          Animation::RandomWalk<W>::Options::Energetic()));
    }
  }

  /**
   * @brief Advances the timeline and draws each ring's trail.
   * @details Steps the timeline, then draws each ring's trail back-to-front,
   * fading color and alpha along the trail.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(rs_timeline_step);
      timeline.step(canvas);
    }

    HS_PROFILE(rs_draw_rings);
    for (int i = 0; i < NUM_RINGS; ++i) {
      Ring &ring = rings[i];
      ring.trail.record(ring.orientation);
      // One fused scan per trail frame (see RingGroup for the blend-order and
      // AA-tail contract).
      deep_tween_frames(ring.trail, [&](const Quaternion *qs, const float *ts,
                                        int count) {
        constexpr int SUB_CAP = decltype(ring.orientation)::CAPACITY;
        Basis bases[SUB_CAP];
        Color4 colors[SUB_CAP];
        alignas(SDF::Ring) unsigned char shape_mem[SUB_CAP * sizeof(SDF::Ring)];
        int slots = 0;
        constexpr float pixel_w = 2.0f * PI_F / W;
        constexpr float MIN_SLOT_ALPHA = 0.001f;
        for (int j = 0; j < count; ++j) {
          float t = ts[j];
          // Length-fade comes from the palette's alpha vignette, not a t term.
          Color4 c = ring.palette->get(1.0f - t);
          c.alpha = c.alpha * params.alpha;
          if (c.alpha <= MIN_SLOT_ALPHA)
            continue;

          // Adaptive thickness: SDF::Ring takes a half-width, so the drawn
          // band is 4px at the trail head/tail and 2px between, before
          // params.thickness scales it.
          float th =
              ((t < 0.01f || t > 0.95f) ? 2.0f * pixel_w : 1.0f * pixel_w) *
              params.thickness;
          bases[slots] = make_basis(qs[j], ring.normal);
          ::new (shape_mem + slots * sizeof(SDF::Ring))
              SDF::Ring(bases[slots], 1.0f, th);
          colors[slots] = c;
          ++slots;
        }
        if (slots == 0)
          return;
        auto *shapes = std::launder(reinterpret_cast<SDF::Ring *>(shape_mem));

        HS_PROFILE(rs_ring_scan);
        Scan::RingGroup::draw<W, H>(
            filters, canvas, shapes, slots,
            [&](int s, const Vector &, Fragment &f) { f.color = colors[s]; },
            params.debug_bb);
      });
    }
  }

private:
  friend struct ::hs_test::effects_tests::RingSpinWhiteBox;

  Timeline timeline;
  Pipeline<W, H> filters;
  Ring *rings = nullptr;

  std::array<BakedPalette, NUM_PALETTES> baked_palettes;

  /**
   * @brief Tunable rendering parameters for the effect.
   */
  struct Params {
    float alpha = 0.5f;     /**< Global trail opacity multiplier in [0, 1]. */
    float thickness = 0.8f; /**< Ring line thickness multiplier (unitless). */
    bool debug_bb = false;  /**< Whether to draw each ring's bounding box. */
  } params;

  // init() allocates the ring pool (each ring carries its TRAIL_LENGTH trail)
  // and one vignette palette LUT per palette, from the persistent arena.
  static constexpr size_t FOOTPRINT_BYTES =
      NUM_RINGS * sizeof(Ring) +
      NUM_PALETTES * BakedPalette::required_arena_bytes();
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "RingSpin persistent footprint exceeds the default partition; "
                "retune TRAIL_LENGTH/NUM_RINGS or carve arenas");
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(RingSpin)
