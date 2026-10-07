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

#include "core/animation/orientation.h"
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
 */
template <int W, int H> class RingSpin : public Effect {
public:
  static constexpr const char *EFFECT_ID = "RingSpin";

  /**
   * @brief Constructs the effect at the W x H canvas resolution.
   * @param strobe Whether the POV driver blanks after each column.
   */
  HS_COLD_MEMBER explicit RingSpin(bool strobe = true)
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = strobe})) {}

  /**
   * @brief Allocates rings, registers params, bakes palettes, and spawns rings.
   */
  HS_COLD_MEMBER void init() override {
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
      return Ring(&baked_palettes[i % NUM_PALETTES].view());
    });
    for (int i = 0; i < NUM_RINGS; ++i) {
      Ring &r = rings[i];
      timeline.add(0, Animation::RandomWalk<W>(
                          r.orientation, math::Y_AXIS, r.noise,
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
      deep_tween_frames(ring.trail, [&](const math::Quaternion *qs,
                                        const float *ts, int count) {
        constexpr int SUB_CAP = decltype(ring.orientation)::CAPACITY;
        math::Basis bases[SUB_CAP];
        Color4 colors[SUB_CAP];
        // SDF::Ring binds its Basis by reference, so it is neither default-
        // constructible nor assignable: slots are placement-new'd into raw
        // storage referencing a bases[] entry that must outlive the slot.
        static_assert(std::is_trivially_destructible_v<SDF::Ring>);
        alignas(SDF::Ring) unsigned char shape_mem[SUB_CAP * sizeof(SDF::Ring)];
        int slots = 0;
        const float pixel_w = math::coarse_pixel_pitch<W, H>();
        // Whole-slot rasterization cut, deliberately above the per-sample
        // MIN_ENCODABLE_ALPHA floor to skip faint ring raster passes.
        constexpr float MIN_SLOT_ALPHA = 0.001f;
        for (int j = 0; j < count; ++j) {
          float t = ts[j];
          // Length-fade comes from the palette's alpha vignette, not a t term.
          Color4 c = ring.palette->get(1.0f - t);
          c.alpha = c.alpha * params.alpha;
          if (c.alpha <= MIN_SLOT_ALPHA)
            continue;

          // SDF::Ring takes a half-width; the band thickens at the trail
          // head and tail.
          float th =
              ((t < 0.01f || t > 0.95f) ? 2.0f * pixel_w : 1.0f * pixel_w) *
              params.thickness;
          bases[slots] = math::make_basis(qs[j], math::Y_AXIS);
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
            [&](int s, const math::Vector &, Fragment &f) {
              f.color = colors[s];
            },
            params.debug_bb);
      });
    }
  }

private:
  friend struct ::hs_test::effects_tests::RingSpinWhiteBox;

  static constexpr int TRAIL_LENGTH = 19; // trail samples per ring
  static constexpr int NUM_RINGS = 4;
  static constexpr int NUM_PALETTES = 4;

  /**
   * @brief One ring: palette, orientation, trail, and random-walk noise.
   * @details The great circle is the Y_AXIS plane under each trail
   * orientation.
   */
  struct Ring {
    const BakedPalette *palette;
    math::Orientation<> orientation;
    Animation::OrientationTrail<math::Orientation<>, TRAIL_LENGTH> trail;
    FastNoiseLite noise;
    /**
     * @brief Constructs a ring drawing from palette @p p.
     * @param p Baked palette used to color the ring's trail.
     */
    Ring(const BakedPalette *p) : palette(p) {}
  };

  Timeline timeline;
  Pipeline<W, H> filters;
  Ring *rings = nullptr;

  std::array<BakedPaletteStorage, NUM_PALETTES> baked_palettes;

  /**
   * @brief Tunable rendering parameters for the effect.
   */
  struct Params {
    float alpha = 0.5f;     /**< Global trail opacity multiplier in [0, 1]. */
    float thickness = 0.8f; /**< Ring line thickness multiplier (unitless). */
    bool debug_bb = false;  /**< Whether to draw each ring's bounding box. */
  } params;

  // Persistent: the ring pool and one vignette palette LUT per palette.
  static constexpr size_t FOOTPRINT_BYTES =
      NUM_RINGS * sizeof(Ring) +
      NUM_PALETTES * BakedPalette::required_arena_bytes();
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "RingSpin persistent footprint exceeds the default partition; "
                "retune TRAIL_LENGTH/NUM_RINGS or carve arenas");
};
