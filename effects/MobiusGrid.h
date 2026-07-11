/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include "core/engine/engine.h"

/**
 * @brief Renders a rotating sphere of latitude rings and longitude lines warped
 *        through a Möbius transform.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Every sampled point is pushed through a Möbius transform so the grid
 *          continuously inflates/deflates between two "holes" on the sphere.
 */
template <int W, int H> class MobiusGrid : public Effect {

public:
  /**
   * @brief Builds the palettes, Möbius generator, and the render filter
   *        pipeline.
   * @details The two HoleRef filters fade geometry near the rotated north/south
   *          holes; Orient applies the spinning orientation; AntiAlias smooths
   *          the rasterized lines.
   */
  HS_COLD_MEMBER MobiusGrid()
      : Effect(W, H, {.strobe = true, .full_frame = decltype(filters)::any_crosses_segments}),
        palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
                BrightnessProfile::FLAT),
        next_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
                     BrightnessProfile::FLAT),
        mobius_gen(timeline), hole_n(Y_AXIS), hole_s(-Y_AXIS),
        filters(Filter::World::HoleRef<W>(hole_n, 1.2f),
                Filter::World::HoleRef<W>(hole_s, 1.2f),
                Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>()) {}

  /**
   * @brief Sizes the arenas, exposes the user params, and arms the timeline
   *        animations.
   * @details Arms a steady Y-axis spin, a periodic palette wipe, and
   *          out-of-phase sine mutations driving the ring/line counts.
   */
  void init() override {
    configure_arenas(GLOBAL_ARENA_SIZE - 8 * 1024, 8 * 1024, 0);

    baked_palette.bake(persistent_arena, palette);

    register_animated_param("Rings", &params.num_rings, 0.0f, 20.0f);
    register_animated_param("Lines", &params.num_lines, 0.0f, 20.0f);
    register_param("Alpha", &params.alpha, 0.0f, 1.0f);

    auto *warp = mobius_gen.spawn_pinned(0, 1.0f, 160, true);
    HS_CHECK(warp, "MobiusGrid: pinned warp spawn must succeed");

    timeline
        .add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 400,
                                       ease_linear, true))
        .add(0, Animation::PeriodicTimer(
                    120, [this](Canvas &) { wipe_palette(); }, true))
        .add(0, Animation::Mutation(params.num_rings,
                                    sin_wave(12.0f, 1.0f, 1.0f, 0.0f), 320,
                                    ease_linear, true, &anims_paused_))
        .add(160, Animation::Mutation(params.num_lines,
                                      sin_wave(12.0f, 1.0f, 1.0f, 0.5f),
                                      320, ease_linear, true, &anims_paused_));
  }

  /**
   * @brief Advances the phase and draws one frame of rings and longitudes.
   * @details Recomputes the counter-rotation that keeps the two holes
   *          diametrically opposed under the Möbius warp before drawing.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    step_wipe_rebake(wipe_pending_, wipe_frames_remaining_, baked_palette, palette);

    // int % 120 before the float cast: a float frame counter would lose integer
    // precision past 2^24.
    float phase = static_cast<float>(timeline.frame() % 120) / 120.0f;

    Vector n_in = Y_AXIS;
    Vector n_trans = mobius_gen.transform(n_in);
    Vector s_in = -Y_AXIS;
    Vector s_trans = mobius_gen.transform(s_in);
    Vector mid = (n_trans + s_trans);
    // Counter-rotation swinging the pole midpoint back onto +Z. When the two
    // transformed poles cancel (mid ~ 0 at the strip singularity) the direction
    // is undefined; leave q at identity rather than feed a zero vector to
    // make_rotation.
    Quaternion q;
    if (mid.length() > 0.001f) {
      mid.normalize();
      q = make_rotation(mid, Z_AXIS);
    }

    hole_n = normalized_or(rotate(n_trans, q), Vector(1, 0, 0));
    hole_s = normalized_or(rotate(s_trans, q), Vector(1, 0, 0));

    draw_axis_rings(canvas, Y_AXIS, params.num_rings, phase, q);
    draw_longitudes(canvas, params.num_lines, phase, q);
  }

private:
  /**
   * @brief Generates a fresh palette and schedules a 60-frame cross-fade into
   *        it.
   */
  void wipe_palette() {
    next_palette = GenerativePalette(GradientShape::CIRCULAR,
                                     HarmonyType::SPLIT_COMPLEMENTARY,
                                     BrightnessProfile::FLAT);
    timeline.add(
        0, Animation::ColorWipe(palette, next_palette, WIPE_FRAMES, ease_linear));
    wipe_frames_remaining_ = WIPE_FRAMES;
    wipe_pending_ = true;
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
   * @tparam ShadeFn Callable (int i, float opacity, Fragment&) coloring a
   *         fragment.
   * @param canvas Render target for this frame.
   * @param num Fractional curve count; ceil(num) curves are drawn.
   * @param q Counter-rotation applied after the Möbius warp.
   * @param curve_fn Supplies the {basis, radius} for curve i.
   * @param shade Colors each fragment for curve i.
   * @details For each curve it asks curve_fn(i) for the {basis, radius}, samples
   *          a spherical polygon, warps every point through the Möbius transform
   *          plus counter-rotation q, then rasterizes via shade. The two
   *          callbacks are all that differs between the ring and longitude passes.
   */
  template <typename CurveFn, typename ShadeFn>
  void draw_curves(Canvas &canvas, float num, const Quaternion &q,
                   CurveFn curve_fn, ShadeFn shade) {
    int count = static_cast<int>(std::ceil(num));
    for (int i = 0; i < count; ++i) {
      ScratchScope frag_guard(scratch_arena_a);
      Curve curve = curve_fn(i);

      Fragments m_fragments;
      // SphericalPolygon::sample emits W/4 + 1 points (one closing overlap).
      m_fragments.bind(scratch_arena_a, W / 4 + 2);
      Plot::SphericalPolygon::sample(m_fragments, curve.basis, curve.radius,
                                     W / 4);

      // Warp the sampled points in place; v0 carries the polyline parameter and
      // the sampler's arc-length/index registers are unused here.
      size_t n = m_fragments.size();
      for (size_t k = 0; k < n; ++k) {
        Vector transformed = mobius_gen.transform(m_fragments[k].pos);
        Fragment &f = m_fragments[k];
        f.pos = normalized_or(rotate(transformed, q), Vector(1, 0, 0));
        f.v0 = (n > 1) ? (float)k / (n - 1) : 0.0f;
        f.v1 = 0.0f;
        f.v2 = 0.0f;
      }

      float opacity = hs::clamp(num - static_cast<float>(i), 0.0f, 1.0f);

      auto fragment_shader = [&](const Vector &, Fragment &f_val) {
        shade(i, opacity, f_val);
      };

      Plot::rasterize<W, H>(filters, canvas, m_fragments, fragment_shader,
                            true);
    }
  }

  /**
   * @brief Draws num latitude rings around normal.
   * @param canvas Render target for this frame.
   * @param normal Axis the rings encircle, a unit vector.
   * @param num Fractional ring count; ceil(num) rings are drawn.
   * @param phase Scroll offset in [0, 1) advancing the ring radii.
   * @param q Counter-rotation applied after the Möbius warp.
   * @details Each ring's radius is spaced logarithmically (log_min..log_max) and
   *          scrolled by phase, mapped through an atan so spacing matches the
   *          stereographic projection.
   */
  void draw_axis_rings(Canvas &canvas, const Vector &normal, float num,
                       float phase, const Quaternion &q) {
    const float log_min = -2.5f;
    const float log_max = 2.5f;
    const float range = log_max - log_min;

    // normal is loop-invariant, so the basis is identical for every ring.
    const Basis ring_basis = make_basis(Quaternion(), normal);

    int cached_i = -1;
    Color4 cached_c;
    draw_curves(
        canvas, num, q,
        [&](int i) -> Curve {
          float t = wrap((static_cast<float>(i) / num) + phase, 1.0f);
          float r_val = expf(log_min + t * range);
          float radius = (4.0f / PI_F) * atanf(1.0f / r_val);
          return {ring_basis, radius};
        },
        [&](int i, float opacity, Fragment &f_val) {
          if (i != cached_i) {
            cached_i = i;
            cached_c = palette.get(static_cast<float>(i) / num);
            cached_c.alpha *= opacity * params.alpha;
          }
          f_val.color = cached_c;
        });
  }

  /**
   * @brief Draws num longitude lines (great circles through the poles).
   * @param canvas Render target for this frame.
   * @param num Fractional line count; ceil(num) lines are drawn.
   * @param phase Scroll offset in [0, 1) advancing the hue gradient.
   * @param q Counter-rotation applied after the Möbius warp.
   * @details The per-line color comes from each fragment's stereographic
   *          conformal radius, so the hue gradient runs pole-to-pole and
   *          scrolls with phase.
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
          // Fall back to X when the normal aligns with Y (cross collapses to 0),
          // matching the normalized_or guards on the sibling bases above.
          Vector u = normalized_or(cross(v, w), Vector(1, 0, 0));
          return {Basis{u, v, w}, 1.0f};
        },
        [&](int, float opacity, Fragment &f_val) {
          float t_line = f_val.v0;
          float z = fast_sinf(t_line * 2.0f * PI_F);

          // Conformal radius R = sqrt((1+z)/(1-z)) is singular at the poles
          // z = +/-1. Branch on the pole before the division so no non-finite
          // intermediate is produced; both poles saturate to coord = 1.0.
          constexpr float POLE_EPS = 1e-6f;
          float coord;
          if (1.0f - fabsf(z) < POLE_EPS) {
            coord = 1.0f;
          } else {
            float R = sqrtf((1.0f + z) / (1.0f - z));
            float log_r = logf(R);
            const float log_min = -2.5f;
            const float log_max = 2.5f;
            float t = (log_r - log_min) / (log_max - log_min);
            coord = wrap(t - phase, 1.0f);
          }
          Color4 c = baked_palette.get(coord);
          c.alpha *= opacity * params.alpha;
          f_val.color = c;
        });
  }

  static constexpr int WIPE_FRAMES = 60; /**< Palette cross-fade duration, in frames. */

  GenerativePalette palette;      /**< Currently displayed palette. */
  GenerativePalette next_palette; /**< Palette being cross-faded toward. */
  BakedPalette baked_palette;     /**< LUT-baked copy of `palette` the shaders sample. */
  int wipe_frames_remaining_ = 0; /**< Frames left to rebake `palette` for an in-flight wipe. */
  bool wipe_pending_ = false;     /**< Wipe armed this frame; it first steps next frame. */
  Timeline timeline;         /**< Drives spin, palette wipe, and mutations. */
  MobiusWarpCircularTransformer<1> mobius_gen; /**< Möbius warp generator. */

  /**
   * @brief User-tunable parameters.
   * @details num_rings/num_lines are animated counts driven by the Mutation
   *          timelines; alpha is the overall opacity multiplier.
   */
  struct Params {
    float num_rings = 0.0f; /**< Animated latitude-ring count. */
    float num_lines = 0.0f; /**< Animated longitude-line count. */
    float alpha = 0.2f;     /**< Overall opacity multiplier in [0, 1]. */
  } params;

  Orientation<> orientation; /**< Spinning render orientation. */

  Vector hole_n; /**< North hole origin, tracking the rotated geometry. */
  Vector hole_s; /**< South hole origin, tracking the rotated geometry. */

  /** @brief Render filter pipeline: two hole fades, orientation, anti-alias. */
  Pipeline<W, H, Filter::World::HoleRef<W>, Filter::World::HoleRef<W>,
           Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(MobiusGrid)
