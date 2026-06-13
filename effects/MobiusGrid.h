/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include "core/effects_engine.h"

/**
 * @brief Renders a rotating sphere of latitude rings and longitude lines warped
 *        through a Möbius transform.
 * @tparam W Cubemap face width in pixels.
 * @tparam H Cubemap face height in pixels.
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
  FLASHMEM MobiusGrid()
      : Effect(W, H),
        palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
                BrightnessProfile::FLAT),
        next_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
                     BrightnessProfile::FLAT),
        mobius_gen(timeline), holeN(Z_AXIS), holeS(-Z_AXIS),
        filters(Filter::World::HoleRef<W>(holeN, 1.2f),
                Filter::World::HoleRef<W>(holeS, 1.2f),
                Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>()) {}

  /**
   * @brief Reports whether the effect wants a background pass.
   * @return Always false; this effect draws no background.
   */
  bool show_bg() const override { return false; }

  /**
   * @brief Sizes the arenas, exposes the user params, and arms the timeline
   *        animations.
   * @details Arms a steady Y-axis spin, a periodic palette wipe, and
   *          out-of-phase sine mutations driving the ring/line counts.
   */
  void init() override {
    // MobiusGrid requires very little persistent memory.
    // Give it 64KB for Scratch A and 64KB for Scratch B to comfortably handle
    // rasterization arrays.
    configure_arenas(GLOBAL_ARENA_SIZE - 128 * 1024, 64 * 1024, 64 * 1024);

    registerParam("Rings", &params.num_rings, 0.0f, 20.0f);
    registerParam("Lines", &params.num_lines, 0.0f, 20.0f);
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    // Rings/Lines are driven by the Mutations below; flag them so the GUI
    // auto-pauses the animation when the user grabs the slider.
    markAnimated("Rings");
    markAnimated("Lines");

    mobius_gen.spawn(0, 1.0f, 160, true);

    timeline
        .add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 400,
                                       ease_mid, true))
        .add(0, Animation::PeriodicTimer(
                    120, [this](auto &) { wipe_palette(); }, true))
        .add(0, Animation::Mutation(params.num_rings,
                                    sin_wave(12.0f, 1.0f, 1.0f, 0.0f), 320,
                                    ease_mid, true, &anims_paused_))
        .add(160, Animation::Mutation(params.num_lines,
                                      sin_wave(12.0f, 1.0f, 1.0f, 0.5f),
                                      320, ease_mid, true, &anims_paused_));
  }

  /**
   * @brief Advances the phase and draws one frame of rings and longitudes.
   * @details Recomputes the counter-rotation that keeps the two holes
   *          diametrically opposed under the Möbius warp before drawing.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    // Effect-local frame counter wrapped to the 120-frame phase period: an int
    // accumulator stays exact and bounded forever, where a float derived from
    // the global frame count would lose integer precision past 2^24 (~77 h).
    // timeline.step() above advanced one frame, so advance in lockstep here.
    phase_frame_ = (phase_frame_ + 1) % 120;
    float phase = static_cast<float>(phase_frame_) / 120.0f;

    // Calculate Stabilizing Counter-Rotation
    Vector n_in = Y_AXIS;
    Vector n_trans = mobius_gen.transform(n_in);
    Vector s_in = -Y_AXIS;
    Vector s_trans = mobius_gen.transform(s_in);
    Vector mid = (n_trans + s_trans);
    Quaternion q;
    if (mid.length() > 0.001f) {
      mid.normalize();
      q = make_rotation(mid, Z_AXIS);
    }

    // Update hole origins to match the rotated geometry. The Möbius transform
    // can collapse a point to ~0 at the strip singularity; guard the normalize.
    holeN = normalized_or(rotate(n_trans, q), Vector(1, 0, 0));
    holeS = normalized_or(rotate(s_trans, q), Vector(1, 0, 0));

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
    timeline.add(0, Animation::ColorWipe(palette, next_palette, 60, ease_mid));
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
   *          callbacks are all that differs between the ring and longitude
   *          passes; templating on them keeps the helper fully inlined.
   */
  template <typename CurveFn, typename ShadeFn>
  void draw_curves(Canvas &canvas, float num, const Quaternion &q,
                   CurveFn curve_fn, ShadeFn shade) {
    int count = static_cast<int>(std::ceil(num));
    for (int i = 0; i < count; ++i) {
      ScratchScope _frag(scratch_arena_a);
      Curve curve = curve_fn(i);

      Fragments m_points;
      // SphericalPolygon::sample emits W/4 + 1 points (one closing overlap);
      // size to that, not a fixed constant, so the bind scales with resolution.
      m_points.bind(scratch_arena_a, W / 4 + 2);
      Plot::SphericalPolygon::sample(m_points, curve.basis, curve.radius,
                                     W / 4);

      Fragments m_fragments;
      m_fragments.bind(scratch_arena_a, W / 4 + 2);
      for (size_t k = 0; k < m_points.size(); ++k) {
        Vector transformed = mobius_gen.transform(m_points[k].pos);
        Fragment f;
        f.pos = normalized_or(rotate(transformed, q), Vector(1, 0, 0));
        f.v0 = (m_points.size() > 1) ? (float)k / (m_points.size() - 1) : 0.0f;
        m_fragments.push_back(f);
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

    // A ring's color depends only on its index i, but the shader runs ~W/4 times
    // per ring, so resolve palette.get once and reuse it. Rings draw
    // sequentially, so a single cache entry keyed on i suffices.
    int cached_i = -1;
    Color4 cached_c;
    draw_curves(
        canvas, num, q,
        [&](int i) -> Curve {
          float t = wrap((static_cast<float>(i) / num) + phase, 1.0f);
          float r_val = expf(log_min + t * range);
          float radius = (4.0f / PI_F) * atanf(1.0f / r_val);
          return {make_basis(Quaternion(), normal), radius};
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
          Vector u = cross(v, w).normalized();
          return {Basis{u, v, w}, 1.0f};
        },
        [&](int, float opacity, Fragment &f_val) {
          float t_line = f_val.v0;
          float z = sinf(t_line * 2.0f * PI_F);

          // Conformal radius of the stereographic longitude. Singular where
          // z = ±1 (t_line = 0.25 and 0.75): the division by (1 - z) blows R to
          // ±inf (or 0), log_r follows, and t goes non-finite. The intended
          // visual is for the line to saturate to its terminal palette color at
          // the pole. Do NOT "fix" by nudging (1 - z) off zero — that changes the
          // endpoint color and splits the two poles to different ends.
          float R = sqrtf((1.0f + z) / (1.0f - z));
          float log_r = logf(R);
          const float log_min = -2.5f;
          const float log_max = 2.5f;
          float t = (log_r - log_min) / (log_max - log_min);

          // Make the pole saturation explicit instead of relying on palette.get
          // clamping a NaN coordinate (fragile if the -fno-finite-math-only
          // guard of finding 369 ever regressed): a non-finite t collapses to the
          // hi endpoint, exactly matching the prior NaN->hi behavior (both poles
          // reach the same terminal color). palette.get now always sees a finite
          // coordinate.
          float coord = wrap(t - phase, 1.0f);
          if (!std::isfinite(t))
            coord = 1.0f;
          Color4 c = palette.get(coord);
          c.alpha *= opacity * params.alpha;
          f_val.color = c;
        });
  }

  GenerativePalette palette;      /**< Currently displayed palette. */
  GenerativePalette next_palette; /**< Palette being cross-faded toward. */
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
  Timeline timeline;         /**< Drives spin, palette wipe, and mutations. */

  Vector holeN; /**< North hole origin, tracking the rotated geometry. */
  Vector holeS; /**< South hole origin, tracking the rotated geometry. */

  /** @brief Effect-local frame counter wrapped to the ring/line phase period. */
  int phase_frame_ = 0;

  /** @brief Render filter pipeline: two hole fades, orientation, anti-alias. */
  Pipeline<W, H, Filter::World::HoleRef<W>, Filter::World::HoleRef<W>,
           Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(MobiusGrid)
