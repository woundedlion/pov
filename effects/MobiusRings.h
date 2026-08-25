/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file MobiusRings.h
 * @brief Rotating sphere of latitude rings and longitude lines warped through
 *        a Mobius transform.
 */

#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"

// Unit-test accessor for the conformal-radius pole branch and the
// counter-rotation singularity guard.
namespace hs_test {
namespace effects_tests {
struct MobiusRingsWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Renders a rotating sphere of latitude rings and longitude lines warped
 *        through a Möbius transform.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Every sampled point is pushed through a Möbius transform so the grid
 *          continuously inflates/deflates between two "holes" on the sphere.
 */
template <int W, int H> class MobiusRings : public Effect {
  // Distinct type keys for Pipeline::get<T>; no per-subclass behaviour.
  class NorthHole : public Filter::World::Hole {
  public:
    using Filter::World::Hole::Hole;
  };

  class SouthHole : public Filter::World::Hole {
  public:
    using Filter::World::Hole::Hole;
  };

public:
  /**
   * @brief Builds the palettes, Möbius generator, and the render filter
   *        pipeline.
   * @details The two Hole filters fade geometry near the rotated north/south
   *          holes; Orient applies the spinning orientation; AntiAlias smooths
   *          the rasterized lines.
   */
  HS_COLD_MEMBER MobiusRings()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})),
        palette(EffectPaletteRecipes::mobius_rings(
            EffectPaletteRecipes::random_base_turns())),
        mobius_gen(timeline),
        filters(NorthHole(Y_AXIS, 1.2f), SouthHole(-Y_AXIS, 1.2f),
                Filter::World::Orient(orientation),
                Filter::Screen::AntiAlias<W, H>()) {}

  // Scratch A holds one curve's fragment buffer (W/4 + 2 samples) and, during
  // the rasterize call it stays live across, rasterize's own sub-step cache.
  static constexpr size_t SCRATCH_A_BYTES = 8 * 1024;
  static_assert(SCRATCH_A_BYTES >= (W / 4 + 2) * sizeof(Fragment) +
                                       Plot::rasterize_scratch_a_bytes<W>(),
                "scratch arena A must fit a curve's fragment buffer and "
                "rasterize's sub-step cache at once");
  using MobiusEntity = typename MobiusWarpCircularTransformer<1>::Entity;
  static constexpr size_t FOOTPRINT_BYTES =
      sizeof(MobiusEntity) + alignof(MobiusEntity) + sizeof(int) +
      alignof(int) + BakedPalette::required_arena_bytes();
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_GLOBAL_ARENA_SIZE - SCRATCH_A_BYTES,
      "MobiusRings persistent footprint exceeds its device partition");

  /**
   * @brief Sizes the arenas, exposes the user params, and arms the timeline
   *        animations.
   * @details Arms a steady Y-axis spin, a periodic palette wipe, and
   *          out-of-phase sine mutations driving the ring/line counts.
   */
  HS_COLD_MEMBER void init() override {
    configure_arenas(GLOBAL_ARENA_SIZE - SCRATCH_A_BYTES, SCRATCH_A_BYTES, 0);

    mobius_gen.init_storage(persistent_arena);
    baked_palette.bake(persistent_arena, palette);

    register_animated_param("Rings", &params.num_rings, 0.0f, 20.0f);
    register_animated_param("Lines", &params.num_lines, 0.0f, 20.0f);
    register_param("Alpha", &params.alpha, 0.0f, 1.0f);

    auto *warp = mobius_gen.spawn_pinned(0, 1.0f, 160, true);
    HS_CHECK(warp, "MobiusRings: pinned warp spawn must succeed");

    timeline
        .add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 400,
                                       ease_linear, true))
        .add(0, Animation::PeriodicTimer(
                    WIPE_PERIOD, [this](Canvas &) { wipe_palette(); }, true))
        .add_pausable(0,
                      Animation::Mutation(params.num_rings,
                                          sin_wave(1.0f, 12.0f, 1.0f, 0.5f),
                                          320, ease_linear, true),
                      &anims_paused)
        .add_pausable(0,
                      Animation::Mutation(params.num_lines,
                                          sin_wave(1.0f, 12.0f, 1.0f, 0.0f),
                                          320, ease_linear, true),
                      &anims_paused);
  }

  /**
   * @brief Advances the phase and draws one frame of rings and longitudes.
   * @details Recomputes the counter-rotation that keeps the two holes
   *          diametrically opposed under the Möbius warp before drawing.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(mg_timeline_step);
      timeline.step(canvas);
    }

    {
      HS_PROFILE(mg_wipe_rebake);
      wipe.step(baked_palette, palette);
    }

    // Integer modulo before the float cast: a float frame counter would lose
    // integer precision past 2^24.
    float phase = static_cast<float>(timeline.frame() % SCROLL_PERIOD) /
                  static_cast<float>(SCROLL_PERIOD);

    Vector n_trans = mobius_gen.transform(Y_AXIS);
    Vector s_trans = mobius_gen.transform(-Y_AXIS);
    Quaternion q = counter_rotation(n_trans + s_trans);

    filters.template get<NorthHole>().set_origin(
        normalized_or(rotate(n_trans, q), Vector(1, 0, 0)));
    filters.template get<SouthHole>().set_origin(
        normalized_or(rotate(s_trans, q), Vector(1, 0, 0)));

    {
      HS_PROFILE(mg_draw_grid);
      {
        HS_PROFILE(mg_rings_draw);
        draw_axis_rings(canvas, Y_AXIS, params.num_rings, phase, q);
      }
      {
        HS_PROFILE(mg_lines_draw);
        draw_longitudes(canvas, params.num_lines, phase, q);
      }
    }
  }

private:
  // Test seam: reaches the conformal-radius pole branch and the
  // counter-rotation singularity guard.
  friend struct ::hs_test::effects_tests::MobiusRingsWhiteBox;

  static constexpr float CONFORMAL_LOG_MIN = -2.5f;
  static constexpr float CONFORMAL_LOG_MAX = 2.5f;
  static constexpr int SCROLL_PERIOD =
      120; /**< Frames per full grid-scroll cycle. */

  /**
   * @brief Counter-rotation swinging the transformed-pole midpoint back onto +Z.
   * @param mid Sum of the two Möbius-transformed poles.
   * @return Rotation mapping the normalized midpoint to +Z, or identity when the
   *         two poles cancel (mid ~ 0 at the strip singularity) and the direction
   *         is undefined — feeding a zero vector to make_rotation is avoided.
   */
  static Quaternion counter_rotation(Vector mid) {
    if (mid.length() > 0.001f) {
      mid.normalize();
      return make_rotation(mid, Z_AXIS);
    }
    return Quaternion();
  }

  /**
   * @brief Maps a longitude height to its scrolled conformal-radius palette
   *        coordinate.
   * @param y Longitude height along the pole axis, in [-1, 1].
   * @param phase Scroll offset subtracted from the coordinate.
   * @return Palette coordinate in [0, 1]; the poles y = +/-1 saturate to 1.0
   *         before the singular conformal radius R = sqrt((1+y)/(1-y)) is
   *         formed, so no non-finite intermediate is produced.
   */
  static float conformal_coord(float y, float phase) {
    constexpr float POLE_EPS = 1e-6f;
    if (1.0f - fabsf(y) < POLE_EPS)
      return 1.0f;
    float log_r = 0.5f * logf((1.0f + y) / (1.0f - y));
    float t =
        (log_r - CONFORMAL_LOG_MIN) / (CONFORMAL_LOG_MAX - CONFORMAL_LOG_MIN);
    return wrap(t - phase, 1.0f);
  }

  /**
   * @brief Generates a fresh palette and schedules a 60-frame cross-fade into
   *        it.
   */
  HS_COLD_MEMBER void wipe_palette() {
    wipe.arm(palette,
             GenerativePalette{EffectPaletteRecipes::mobius_rings(
                                   EffectPaletteRecipes::random_base_turns())}
                 .snapshot(),
             WIPE_FRAMES);
    timeline.add(0, Animation::ColorWipe(palette, wipe.start, wipe.target,
                                         WIPE_FRAMES, ease_linear));
  }

  /**
   * @brief A single Möbius-warped curve to draw.
   * @details Holds the spherical-polygon basis and the sampling radius for a
   *          given ring/line index.
   */
  struct Curve {
    Basis basis;  /**< Orthonormal basis of the spherical polygon. */
    float radius; /**< Sampling radius of the curve, in unit-sphere coords. */
  };

  /**
   * @brief Shared body for both grid passes: samples, warps, and rasterizes
   *        ceil(num) curves.
   * @tparam CurveFn Callable (int i) -> Curve yielding the {basis, radius}.
   * @tparam ShaderFn Callable (int i, float opacity) -> fragment shader
   *         (const Vector&, Fragment&).
   * @param canvas Render target for this frame.
   * @param num Fractional curve count; ceil(num) curves are drawn.
   * @param q Counter-rotation applied after the Möbius warp.
   * @param curve_fn Supplies the {basis, radius} for curve i.
   * @param make_shader Builds curve i's fragment shader, once per curve.
   * @details For each curve it asks curve_fn(i) for the {basis, radius}, samples
   *          a spherical polygon, warps every point through the Möbius transform
   *          plus counter-rotation q, then rasterizes through the shader
   *          make_shader(i, opacity) returns. The two callbacks are all that
   *          differs between the ring and longitude passes.
   */
  template <typename CurveFn, typename ShaderFn>
  void draw_curves(Canvas &canvas, float num, const Quaternion &q,
                   CurveFn curve_fn, ShaderFn make_shader) {
    int count = static_cast<int>(std::ceil(num));
    for (int i = 0; i < count; ++i) {
      ScratchScope frag_guard(scratch_arena_a);
      Curve curve = curve_fn(i);

      Fragments m_fragments;
      // Polygon::sample emits W/4 + 1 points (one closing overlap).
      m_fragments.bind(scratch_arena_a, W / 4 + 2);
      Plot::Polygon<Plot::GeodesicProjection>::sample(m_fragments, curve.basis,
                                                      curve.radius, W / 4);

      // Warp the sampled points in place; v0 carries the polyline parameter and
      // the sampler's arc-length/index registers are unused here.
      size_t n = m_fragments.size();
      for (size_t k = 0; k < n; ++k) {
        Vector transformed = mobius_gen.transform(m_fragments[k].pos);
        Fragment &f = m_fragments[k];
        f.pos = normalized_or(rotate(transformed, q), Vector(1, 0, 0));
        f.v1 = 0.0f;
        f.v2 = 0.0f;
      }

      float opacity = hs::clamp(num - static_cast<float>(i), 0.0f, 1.0f);

      auto fragment_shader = make_shader(i, opacity);

      // Polygon::sample closes the ring with an overlap vertex, so the seam
      // pixel belongs to the opening segment.
      Plot::rasterize<W, H>(filters, canvas, m_fragments, fragment_shader,
                            {.omit_end = true});
    }
  }

  /**
   * @brief Draws num latitude rings around normal.
   * @param canvas Render target for this frame.
   * @param normal Axis the rings encircle, a unit vector.
   * @param num Fractional ring count; ceil(num) rings are drawn.
   * @param phase Scroll offset in [0, 1) advancing the ring radii.
   * @param q Counter-rotation applied after the Möbius warp.
   * @details Each ring's radius spans the shared conformal log range and is
   *          scrolled by phase, mapped through an atan so spacing matches the
   *          stereographic projection.
   */
  void draw_axis_rings(Canvas &canvas, const Vector &normal, float num,
                       float phase, const Quaternion &q) {
    const float range = CONFORMAL_LOG_MAX - CONFORMAL_LOG_MIN;

    // normal is loop-invariant, so the basis is identical for every ring.
    const Basis ring_basis = make_basis(Quaternion(), normal);

    draw_curves(
        canvas, num, q,
        [&](int i) -> Curve {
          float t = wrap((static_cast<float>(i) / num) + phase, 1.0f);
          float r_val = expf(CONFORMAL_LOG_MIN + t * range);
          float radius = (4.0f / PI_F) * atanf(1.0f / r_val);
          return {ring_basis, radius};
        },
        [&](int i, float opacity) {
          Color4 c = baked_palette.get(static_cast<float>(i) / num);
          c.alpha *= opacity * params.alpha;
          return [c](const Vector &, Fragment &f_val) { f_val.color = c; };
        });
  }

  /**
   * @brief Draws num longitude lines (great circles through the poles).
   * @param canvas Render target for this frame.
   * @param num Fractional line count; ceil(num) lines are drawn.
   * @param phase Scroll offset in [0, 1) advancing the hue gradient.
   * @param q Counter-rotation applied after the Möbius warp.
   * @details The color comes from the conformal radius of the source great
   *          circle's parameter, not of the warped fragment position, so the
   *          hue gradient is anchored to the pre-warp curve and scrolls with
   *          phase.
   */
  void draw_longitudes(Canvas &canvas, float num, float phase,
                       const Quaternion &q) {
    draw_curves(
        canvas, num, q,
        [&](int i) -> Curve {
          float theta = (static_cast<float>(i) / num) * PI_F;
          Vector normal(cosf(theta), 0.0f, -sinf(theta));
          // Explicit basis construction to match JS texture alignment.
          Vector v = normal;
          Vector w = Y_AXIS;
          Vector u = cross(v, w);
          return {Basis{u, v, w}, 1.0f};
        },
        [&](int, float opacity) {
          const float alpha = opacity * params.alpha;
          return [this, phase, alpha](const Vector &, Fragment &f_val) {
            float y = fast_sinf(f_val.v0 * 2.0f * PI_F);
            Color4 c = baked_palette.get(conformal_coord(y, phase));
            c.alpha *= alpha;
            f_val.color = c;
          };
        });
  }

  static constexpr int WIPE_FRAMES =
      60; /**< Palette cross-fade duration, in frames. */
  static constexpr int WIPE_PERIOD =
      120; /**< Frames between palette rollovers. */
  static_assert(WIPE_PERIOD > WIPE_FRAMES,
                "MobiusRings: a rollover firing mid-wipe would clobber the "
                "snapshots the live ColorWipe still references");

  GenerativePalette palette; /**< Currently displayed palette. */
  PaletteWipe wipe;          /**< Cross-fade state of the palette rollover. */
  BakedPalette
      baked_palette; /**< LUT-baked copy of `palette` the shaders sample. */
  /**
   * @brief User-tunable parameters.
   * @details num_rings/num_lines are animated counts driven by the Mutation
   *          timelines; alpha is the overall opacity multiplier. Declared before
   *          `timeline` so the counts outlive the Mutations that point here.
   */
  struct Params {
    /** Animated latitude-ring count; the default is the ring Mutation's opening
     * sample, so a pause held from before init() still draws a grid. */
    float num_rings = 12.0f;
    /** Animated longitude-line count; the default is the line Mutation's
     * opening sample, so the lines are drawn before that Mutation starts and
     * under a pause held from before init(). */
    float num_lines = 1.0f;
    float alpha = 0.2f; /**< Overall opacity multiplier in [0, 1]. */
  } params;

  /**
   * @brief Spinning render orientation.
   * @details Declared before `timeline` so it outlives the Rotation that points
   * here, which ~Timeline clears on teardown.
   */
  Orientation<> orientation;
  Timeline timeline; /**< Drives spin, palette wipe, and mutations. */
  MobiusWarpCircularTransformer<1> mobius_gen; /**< Möbius warp generator. */

  /** @brief Render filter pipeline: two hole fades, orientation, anti-alias. */
  Pipeline<W, H, NorthHole, SouthHole, Filter::World::Orient,
           Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "core/control/registry.h"
REGISTER_EFFECT(MobiusRings)
