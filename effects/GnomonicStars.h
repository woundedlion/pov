/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file GnomonicStars.h
 * @brief Polygon stars scattered over a Fibonacci spiral and warped by an
 *        evolving Mobius transform.
 */

#include "core/animation/orientation.h"
#include "core/engine/engine.h"

namespace hs_test {
namespace effects_tests {
struct GnomonicStarsWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Scatters polygon "stars" over a Fibonacci spiral on the sphere and
 *        warps the field with an evolving Möbius transform.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details A Languid RandomWalk slowly reorients the whole field.
 */
template <int W, int H> class GnomonicStars : public Effect {
public:
  static constexpr const char *EFFECT_ID = "GnomonicStars";

  /**
   * @brief Constructs the effect at face resolution W x H.
   */
  HS_COLD_MEMBER GnomonicStars()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})),
        orientation(), timeline(), transformer(timeline) {}

  /**
   * @brief Registers the user params and arms the timeline.
   * @details Arms an infinite Möbius warp and a Languid RandomWalk.
   */
  HS_COLD_MEMBER void init() override {
    transformer.init_storage(persistent_arena);
    // Sized to MAX_POINTS so a live "Points" change never reallocates.
    spiral_cache = persistent_arena.allocate_n<math::Vector>(MAX_POINTS);

    register_int_param("Points", &params.points, 100, MAX_POINTS);
    register_param("Radius", &params.star_radius, 0.7f * radius_px(),
                   7.0f * radius_px());
    register_int_param("Sides", &params.star_sides, 3, 8);

    // Args are (scale, speed); speed is overwritten from params each frame.
    warp = transformer.spawn_pinned(0, 0.5f, 0.0f);
    HS_CHECK(warp, "GnomonicStars: pinned warp spawn must succeed");
    register_param("Warp Speed", &params.warp_speed, 0.0f, 1.0f);
    register_param("Debug BB", &params.debug_bb);

    baked_palette.bake(persistent_arena, Palettes::MANGO_PEEL);

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, math::Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));
  }

  /**
   * @brief Advances the timeline and renders the warped star field.
   * @details Draws each spiral point as a star whose color is a Y-based gradient
   *          and whose basis carries the current orientation and warp.
   */
  void draw_frame() override {
    Canvas canvas(*this);

    // Mirror the slider into the warp before the timeline advances it.
    warp->set_speed(params.warp_speed);

    {
      HS_PROFILE(gn_timeline_step);
      timeline.step(canvas);
    }

    auto fragment_shader = [this](const math::Vector &p, Fragment &frag) {
      float t = (p.y + 1.0f) * 0.5f; // Y in [-1, 1] -> gradient t in [0, 1]
      Color4 c = baked_palette.get(t);
      frag.color = c;
    };

    const int points = params.points;
    const float radius = params.star_radius;
    const int sides = params.star_sides;

    // The base spiral depends only on (points, i).
    if (points != cached_points) {
      HS_PROFILE(gn_spiral_build);
      // eps 0.5 puts both endpoints half a step off their pole.
      for (int i = 0; i < points; i++) {
        spiral_cache[i] = math::fib_spiral(points, /*eps=*/0.5f, i);
      }
      cached_points = points;
    }

    {
      HS_PROFILE(gn_draw_stars);
      for (int i = 0; i < points; i++) {
        math::Vector v = transformer.transform(spiral_cache[i]);

        // make_basis() applies the orientation; pass the unrotated point.
        math::Basis basis = math::make_basis(orientation.get(), v);

        {
          HS_PROFILE(gn_star_scan);
          Scan::Star::draw<W, H>(filters, canvas, basis, radius, sides,
                                 fragment_shader, 0.0f, params.debug_bb);
        }
      }
    }
  }

private:
  friend struct ::hs_test::effects_tests::GnomonicStarsWhiteBox;

  /** @brief Spiral-cache capacity; equals the "Points" slider's upper bound. */
  static constexpr int MAX_POINTS = 2000;

  /**
   * @brief The coarser pixel pitch expressed in Scan::Star radius units.
   * @details A radius of 1 spans pi/2 radians. Multiples of this pitch cover
   *          at least that many rows and columns at the equator.
   */
  static constexpr float radius_px() {
    return math::coarse_pixel_pitch<W, H>() * 2.0f / math::PI_F;
  }

  // Persistent allocations: warp pool and phases, spiral lattice, palette LUT.
  using MobiusEntity = typename MobiusWarpGnomonicTransformer<1>::Entity;
  static constexpr size_t FOOTPRINT_BYTES =
      sizeof(MobiusEntity) + alignof(MobiusEntity) + sizeof(int) +
      alignof(int) + MAX_POINTS * sizeof(math::Vector) + alignof(math::Vector) +
      BakedPalette::required_arena_bytes() + 8 * sizeof(double) +
      alignof(double) - 1;
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "GnomonicStars persistent footprint exceeds the default "
                "partition; retune MAX_POINTS or carve arenas");

  math::Orientation<> orientation; /**< Current field orientation quaternion. */
  FastNoiseLite noise;             /**< Noise source driving the RandomWalk. */
  Timeline timeline;               /**< Animation timeline for warp and walk. */
  Pipeline<W, H> filters;          /**< Render filter pipeline for star scan. */

  math::Vector *spiral_cache =
      nullptr;           /**< Persistent base lattice, MAX_POINTS slots. */
  int cached_points = 0; /**< Point count the cache holds (0 = unbuilt). */
  BakedPaletteStorage
      baked_palette; /**< LUT-baked MANGO_PEEL sampled by the shader. */

  MobiusWarpGnomonicTransformer<1>
      transformer; /**< Evolving Möbius warp applied per point. */
  Animation::MobiusWarpEvolving *warp =
      nullptr; /**< Pinned warp handle; mirrors params.warp_speed each frame. */

  /**
   * @brief Live-tunable controls for the star field.
   */
  struct Params {
    int points = 600; /**< Number of stars scattered on the spiral. */
    float star_radius =
        1.4f * radius_px();    /**< Per-star circumradius, ~1.4 px at any W. */
    int star_sides = 4;        /**< Polygon side count per star. */
    float warp_speed = 0.035f; /**< Möbius warp evolution speed. */
    bool debug_bb = false; /**< When true, draws each star's bounding box. */
  } params;
};
