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

#include "core/platform/effect_identity.h"
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
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY =
      hs::effect_identities::GnomonicStars::VALUE;
  /// Canonical machine ID.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;

  /**
   * @brief Constructs the effect at the W x H canvas resolution.
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
    drawn_centres = persistent_arena.allocate_n<math::Vector>(MAX_POINTS);

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
   */
  void draw_frame() override {
    Canvas canvas(*this);

    // Mirror the slider into the warp before the timeline advances it.
    warp->set_speed(params.warp_speed);

    {
      HS_PROFILE(gn_timeline_step);
      timeline.step(canvas);
    }

    refresh_spiral();
    draw_stars(canvas);
  }

private:
  /**
   * @brief Rebuilds the base spiral after a Points change, which starts the
   *        stars' drift history over.
   */
  HS_FLASH_MEMBER void refresh_spiral() {
    const int points = params.points;
    // The base spiral depends only on (points, i).
    if (points == cached_points)
      return;
    HS_PROFILE(gn_spiral_build);
    // eps 0.5 puts both endpoints half a step off their pole.
    for (int i = 0; i < points; i++)
      spiral_cache[i] = math::fib_spiral(points, /*eps=*/0.5f, i);
    cached_points = points;
    warp_drawn = false;
  }

  /**
   * @brief Scans every star.
   * @param canvas Target canvas.
   * @details Draws each spiral point as a star whose color is a Y-based
   * gradient and whose basis carries the current orientation and warp. Each
   * column sees the orientation at its sweep time and each star's warped
   * centre moved linearly from where the last frame drew it; a star whose
   * motion barely changes across the columns that see it is drawn where
   * their middle column saw it.
   */
  HS_HOT_FLASH_MEMBER void draw_stars(Canvas &canvas) {
    auto fragment_shader = [this](const math::Vector &p, Fragment &frag) {
      float t = (p.y + 1.0f) * 0.5f; // Y in [-1, 1] -> gradient t in [0, 1]
      Color4 c = baked_palette.get(t);
      frag.color = c;
    };

    const int points = params.points;
    const float radius = params.star_radius;
    const int sides = params.star_sides;
    const bool drifting = warp_drawn && !Spacetime::flat_surface();

    HS_PROFILE(gn_draw_stars);
    ScratchScope column_guard(scratch_arena_a);
    Scan::ColumnBasis column_basis = Scan::prepare_column_basis<W, H>(
        scratch_arena_a, orientation, drifting);
    Scan::ColumnBasisScope column_scope(canvas, column_basis);
    // A star that cannot reach the clip over the frame, drifting included,
    // draws nothing.
    const bool culls = radius <= 1.0f && !params.debug_bb && !canvas.debug();
    const float reach = radius * (math::PI_F / 2.0f);
    const float sweep = column_basis.azimuth ? column_basis.sweep_arc : 0.0f;
    const ClipRegion &clip = Scan::source_clip<W, H>(filters, canvas);
    const auto cap =
        Plot::make_swept_cap_clip<W, H>(clip, reach, sweep + CULL_DRIFT);
    const auto warp = [this](const math::Vector &point) HS_O3_FN HS_HOT_INLINE
        __attribute__((flatten)) { return transformer.transform(point); };
    for (int i = 0; i < points; i++) {
      math::Vector v = warp(spiral_cache[i]);
      const math::Vector drawn = drawn_centres[i];
      drawn_centres[i] = v;

      // make_basis() applies the orientation; pass the unrotated point.
      math::Basis basis = math::make_basis(orientation.get(), v);
      if (culls) {
        // A chord c subtends at most pi / 2 * c of arc.
        const math::Vector chord = v - drawn;
        const float chord_sq = drifting ? math::dot(chord, chord) : 0.0f;
        if (!Plot::cap_may_touch_clip<H>(
                chord_sq <= CULL_CHORD_SQ
                    ? cap
                    : Plot::make_swept_cap_clip<W, H>(
                          clip, reach,
                          sweep + 0.5f * math::PI_F * sqrtf(chord_sq)),
                basis.v))
          continue;
      }
      if (drifting) {
        const math::Vector drift = star_drift(drawn, v);
        column_basis.drift =
            Scan::settle_drift<W, H>(column_basis, drift, radius, basis)
                ? math::Vector(0.0f, 0.0f, 0.0f)
                : drift;
      }

      {
        HS_PROFILE(gn_star_scan);
        Scan::Star::draw<W, H>(filters, canvas, basis, radius, sides,
                               fragment_shader, 0.0f, params.debug_bb);
      }
    }
    warp_drawn = true;
  }

  /**
   * @brief A star's displacement over the frame at the current orientation.
   * @param previous The star's warped centre where the last frame drew it.
   * @param warped The star's warped centre under the current warp.
   * @return Rotation vector carrying @p previous to @p warped along their
   *         great circle, at the current orientation; zero without motion or
   *         for antipodal centres, which name no circle.
   */
  HS_FLASH_MEMBER math::Vector star_drift(const math::Vector &previous,
                                          const math::Vector &warped) const {
    const math::Vector axis = math::cross(previous, warped);
    const float sine_sq = math::dot(axis, axis);
    if (!(sine_sq > 0.0f))
      return math::Vector(0.0f, 0.0f, 0.0f);
    const float cosine = math::dot(previous, warped);
    // asin(s) / s by its series below CULL_DRIFT; the angle over s past it.
    const float scale =
        cosine > 0.0f && sine_sq <= CULL_SINE_SQ
            ? 1.0f + sine_sq *
                         (1.0f / 6.0f +
                          sine_sq * (3.0f / 40.0f + sine_sq * (5.0f / 112.0f)))
            : atan2f(sqrtf(sine_sq), cosine) / sqrtf(sine_sq);
    return math::rotate(axis * scale, orientation.get());
  }

  /** @brief Drift the draw's shared cull pad covers: eight columns. A star
   *  moving further is culled against its own swept cap. */
  static constexpr float CULL_DRIFT = 8.0f * math::RADIANS_PER_COLUMN<W>;
  /** @brief Squared sine of CULL_DRIFT, by its series. */
  static constexpr float CULL_SINE_SQ =
      (CULL_DRIFT * (1.0f - CULL_DRIFT * CULL_DRIFT / 6.0f)) *
      (CULL_DRIFT * (1.0f - CULL_DRIFT * CULL_DRIFT / 6.0f));
  /** @brief Squared chord whose arc the shared cull pad covers. */
  static constexpr float CULL_CHORD_SQ =
      (2.0f * CULL_DRIFT / math::PI_F) * (2.0f * CULL_DRIFT / math::PI_F);

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

  // Persistent allocations: warp pool and phases, spiral lattice, drawn
  // centres, palette LUT.
  using MobiusEntity = typename MobiusWarpGnomonicTransformer<1>::Entity;
  static constexpr size_t FOOTPRINT_BYTES =
      sizeof(MobiusEntity) + alignof(MobiusEntity) + sizeof(int) +
      alignof(int) +
      2 * (MAX_POINTS * sizeof(math::Vector) + alignof(math::Vector)) +
      BakedPalette::required_arena_bytes() +
      Animation::MobiusWarpEvolving::CHANNEL_COUNT * sizeof(double) +
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
  /// Persistent warped centres of the last star draw, MAX_POINTS slots.
  math::Vector *drawn_centres = nullptr;
  bool warp_drawn = false; /**< Whether drawn_centres holds this spiral's last
                              draw. */
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
