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

#include "core/platform/effect_identity.h"
#include "core/animation/orientation.h"
#include "effects/common/palette_recipes.h"
#include "core/engine/engine.h"
#include <memory>

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
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY =
      hs::effect_identities::MobiusRings::VALUE;
  /// Canonical machine ID.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;

  /**
   * @brief Builds the palettes, Möbius generator, and the render filter
   *        pipeline.
   */
  HS_COLD_MEMBER MobiusRings()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})),
        palette(EffectPaletteRecipes::mobius_rings(
            EffectPaletteRecipes::random_base_turns())),
        mobius_gen(timeline),
        filters(NorthHole(math::Y_AXIS, 1.2f), SouthHole(-math::Y_AXIS, 1.2f),
                Filter::World::Orient(orientation),
                Filter::Screen::AntiAlias<W, H>()) {}

  /**
   * @brief Sizes the arenas, exposes the user params, and arms the timeline
   *        animations.
   * @details Arms a steady Y-axis spin, a periodic palette wipe, and
   *          out-of-phase sine mutations driving the ring/line counts.
   */
  HS_COLD_MEMBER void init() override {
    ArenaSplit{SCRATCH_A_BYTES, SCRATCH_B_BYTES}.configure();
    history.reset();

    mobius_gen.init_storage(persistent_arena);
    filters.init_storage(persistent_arena);
    baked_palette.bake(persistent_arena, palette);

    register_animated_param("Rings", &params.num_rings, 0.0f, 20.0f);
    register_animated_param("Lines", &params.num_lines, 0.0f, 20.0f);
    register_param("Alpha", &params.alpha, 0.0f, 1.0f);

    auto *warp = mobius_gen.spawn_pinned(0, 1.0f, 160);
    HS_CHECK(warp, "MobiusRings: pinned warp spawn must succeed");

    timeline
        .add(0,
             Animation::Rotation<W>(orientation, math::Y_AXIS, 2 * math::PI_F,
                                    400, math::ease_linear, {.repeat = true}))
        .add(0, Animation::PeriodicTimer(WIPE_PERIOD,
                                         [this](Canvas &) { wipe_palette(); },
                                         {.repeat = true}))
        .add_pausable(
            0,
            Animation::Mutation(params.num_rings,
                                math::sin_wave(1.0f, 12.0f, 1.0f, 0.5f), 320,
                                math::ease_linear, {.repeat = true}),
            &anims_paused)
        .add_pausable(
            0,
            Animation::Mutation(params.num_lines,
                                math::sin_wave(1.0f, 12.0f, 1.0f, 0.0f), 320,
                                math::ease_linear, {.repeat = true}),
            &anims_paused);
  }

  /**
   * @brief Advances the phase and draws one frame of rings and longitudes.
   * @details Recomputes the holes' counter-rotation before drawing.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(mr_timeline_step);
      timeline.step(canvas);
    }

    {
      HS_PROFILE(mr_wipe_rebake);
      wipe.step(baked_palette, palette);
    }

    // Integer modulo before the float cast: a float frame counter would lose
    // integer precision past 2^24.
    float phase = static_cast<float>(timeline.frame() % SCROLL_PERIOD) /
                  static_cast<float>(SCROLL_PERIOD);

    math::Vector n_trans = mobius_gen.transform(math::Y_AXIS);
    math::Vector s_trans = mobius_gen.transform(-math::Y_AXIS);
    math::Quaternion q = counter_rotation(n_trans + s_trans);

    filters.template get<NorthHole>().set_origin(
        math::normalized_or(math::rotate(n_trans, q), math::Vector(1, 0, 0)));
    filters.template get<SouthHole>().set_origin(
        math::normalized_or(math::rotate(s_trans, q), math::Vector(1, 0, 0)));

    const Endpoint end{mobius_gen.active_params(0),
                       q,
                       phase,
                       params.alpha,
                       params.num_rings,
                       params.num_lines};
    const Endpoint start = history.start(end);
    history.capture(end);
    if (!Spacetime::flat_surface()) {
      draw_timed_grid(canvas, start, end);
      return;
    }

    {
      HS_PROFILE(mr_draw_grid);
      {
        HS_PROFILE(mr_rings_draw);
        draw_axis_rings(canvas, math::Y_AXIS, params.num_rings, phase, q);
      }
      {
        HS_PROFILE(mr_lines_draw);
        draw_longitudes(canvas, params.num_lines, phase, q);
      }
    }
  }

private:
  // Test seam: reaches the conformal-radius pole branch and the
  // counter-rotation singularity guard.
  friend struct ::hs_test::effects_tests::MobiusRingsWhiteBox;

  // Scratch A: one curve's fragment buffer plus rasterize's sub-step cache.
  static constexpr size_t SCRATCH_A_BYTES = 8 * 1024;
  static_assert(SCRATCH_A_BYTES >= (W / 4 + 2) * sizeof(Fragment) +
                                       Plot::rasterize_scratch_a_bytes<W>(),
                "scratch arena A must fit a curve's fragment buffer and "
                "rasterize's sub-step cache at once");
  /// Vertex captures and column-time pole origins.
  static constexpr size_t SCRATCH_B_BYTES = 24 * 1024;
  using MobiusEntity = typename MobiusWarpCircularTransformer<1>::Entity;
  static constexpr size_t FOOTPRINT_BYTES =
      sizeof(MobiusEntity) + alignof(MobiusEntity) + sizeof(int) +
      alignof(int) + Filter::World::Orient::storage_bytes(W) +
      BakedPalette::required_arena_bytes();
  static_assert(
      FOOTPRINT_BYTES <=
          ArenaSplit{SCRATCH_A_BYTES, SCRATCH_B_BYTES}.device_persistent(),
      "MobiusRings persistent footprint exceeds its device partition");

  static constexpr float CONFORMAL_LOG_MIN = -2.5f;
  static constexpr float CONFORMAL_LOG_MAX = 2.5f;
  static constexpr int SCROLL_PERIOD =
      120; /**< Frames per full grid-scroll cycle. */

  /**
   * @brief Counter-rotation swinging the transformed-pole midpoint back onto +Z.
   * @param mid Sum of the two Möbius-transformed poles.
   * @return Rotation mapping the normalized midpoint to +Z, or identity when the
   *         two poles cancel (mid ~ 0 at the strip singularity).
   */
  static math::Quaternion counter_rotation(math::Vector mid) {
    if (mid.length() > 0.001f) {
      mid.normalize();
      return math::make_rotation(mid, math::Z_AXIS);
    }
    return math::Quaternion();
  }

  /**
   * @brief Maps a longitude height to its scrolled conformal-radius palette
   *        coordinate.
   * @param y Longitude height along the pole axis, in [-1, 1].
   * @param phase Scroll offset subtracted from the coordinate.
   * @return Palette coordinate in [0, 1]; the poles y = +/-1 saturate to 1.0
   *         before the singular conformal radius R = sqrt((1+y)/(1-y)) is
   *         formed.
   */
  static float conformal_coord(float y, float phase) {
    constexpr float POLE_EPS = 1e-6f;
    if (1.0f - fabsf(y) < POLE_EPS)
      return 1.0f;
    float log_r = 0.5f * logf((1.0f + y) / (1.0f - y));
    float t =
        (log_r - CONFORMAL_LOG_MIN) / (CONFORMAL_LOG_MAX - CONFORMAL_LOG_MIN);
    return math::wrap(t - phase, 1.0f);
  }

  /**
   * @brief Generates a fresh palette and schedules a WIPE_FRAMES cross-fade
   *        into it.
   */
  HS_COLD_MEMBER void wipe_palette() {
    wipe.arm(palette,
             GenerativePalette{EffectPaletteRecipes::mobius_rings(
                                   EffectPaletteRecipes::random_base_turns())}
                 .snapshot(),
             WIPE_FRAMES);
    timeline.add(0, Animation::ColorWipe(palette, wipe.start, wipe.target,
                                         WIPE_FRAMES, math::ease_linear));
  }

  static constexpr int CURVE_SIDES = W / 4;

  /**
   * @brief A single Möbius-warped curve to draw.
   * @details Holds the spherical-polygon basis and the sampling radius for a
   *          given ring/line index.
   */
  struct Curve {
    math::Basis basis; /**< Orthonormal basis of the spherical polygon. */
    float radius; /**< Angular radius as a fraction of a hemisphere, in [0, 2];
                     1 is a great circle. */
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
   * @details Each curve is sampled as a spherical polygon, warped through the
   *          Möbius transform plus counter-rotation q, then rasterized.
   */
  template <typename CurveFn, typename ShaderFn>
  void draw_curves(Canvas &canvas, float num, const math::Quaternion &q,
                   CurveFn curve_fn, ShaderFn make_shader) {
    int count = static_cast<int>(std::ceil(num));
    for (int i = 0; i < count; ++i) {
      ScratchScope frag_guard(scratch_arena_a);
      Curve curve = curve_fn(i);

      Fragments m_fragments;
      // Polygon::sample emits CURVE_SIDES + 1 points (one closing overlap).
      m_fragments.bind(scratch_arena_a, CURVE_SIDES + 2);
      Plot::Polygon<Plot::GeodesicProjection>::sample(
          m_fragments, curve.basis, curve.radius, CURVE_SIDES);

      // Warp the sampled points in place; v0 carries the polyline parameter and
      // the sampler's arc-length/index registers are unused here.
      size_t n = m_fragments.size();
      for (size_t k = 0; k < n; ++k) {
        math::Vector transformed = mobius_gen.transform(m_fragments[k].pos);
        Fragment &f = m_fragments[k];
        f.pos = math::normalized_or(math::rotate(transformed, q),
                                    math::Vector(1, 0, 0));
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
  void draw_axis_rings(Canvas &canvas, const math::Vector &normal, float num,
                       float phase, const math::Quaternion &q) {
    const float range = CONFORMAL_LOG_MAX - CONFORMAL_LOG_MIN;

    const math::Basis ring_basis = math::make_basis(math::Quaternion(), normal);

    draw_curves(
        canvas, num, q,
        [&](int i) -> Curve {
          float t = math::wrap((static_cast<float>(i) / num) + phase, 1.0f);
          float r_val = expf(CONFORMAL_LOG_MIN + t * range);
          float radius = (4.0f / math::PI_F) * atanf(1.0f / r_val);
          return {ring_basis, radius};
        },
        [&](int i, float opacity) {
          Color4 c = baked_palette.get(static_cast<float>(i) / num);
          c.alpha *= opacity * params.alpha;
          return
              [c](const math::Vector &, Fragment &f_val) { f_val.color = c; };
        });
  }

  /**
   * @brief Draws num longitude lines (great circles through the poles).
   * @param canvas Render target for this frame.
   * @param num Fractional line count; ceil(num) lines are drawn.
   * @param phase Scroll offset in [0, 1) advancing the hue gradient.
   * @param q Counter-rotation applied after the Möbius warp.
   * @details The color comes from the conformal radius of the source great
   *          circle's parameter, so the hue gradient is anchored to the pre-warp
   *          curve and scrolls with phase.
   */
  void draw_longitudes(Canvas &canvas, float num, float phase,
                       const math::Quaternion &q) {
    draw_curves(
        canvas, num, q,
        [&](int i) -> Curve {
          float theta = (static_cast<float>(i) / num) * math::PI_F;
          math::Vector normal(cosf(theta), 0.0f, -sinf(theta));
          math::Vector v = normal;
          math::Vector w = math::Y_AXIS;
          math::Vector u = math::cross(v, w);
          return {math::Basis{u, v, w}, 1.0f};
        },
        [&](int, float opacity) {
          const float alpha = opacity * params.alpha;
          return [this, phase, alpha](const math::Vector &, Fragment &f_val) {
            float y = math::fast_sinf((f_val.v0 + 0.5f / CURVE_SIDES) * 2.0f *
                                      math::PI_F);
            Color4 c = baked_palette.get(conformal_coord(y, phase));
            c.alpha *= alpha;
            f_val.color = c;
          };
        });
  }

  /** @brief Retained warp, counter-rotation and colour endpoint. */
  struct Endpoint {
    math::MobiusParams warp;   ///< Plane transformation coefficients.
    math::Quaternion rotation; ///< Counter-rotation.
    float phase = 0.0f;        ///< Palette scroll in turns.
    float alpha = 0.0f;        ///< Overall opacity.
    float rings = 0.0f;        ///< Latitude count.
    float lines = 0.0f;        ///< Longitude count.
  };
  Spacetime::Retained<Endpoint> history; ///< Previous animation endpoint.

  /** @brief Blends the coefficients of two plane transformations.
   * @param a,b Frame endpoints.
   * @param t Frame fraction.
   * @return Coefficients at t. */
  static math::MobiusParams warp_at(const Endpoint &a, const Endpoint &b,
                                    float t) {
    return {a.warp.a + (b.warp.a - a.warp.a) * math::Complex(t, 0),
            a.warp.b + (b.warp.b - a.warp.b) * math::Complex(t, 0),
            a.warp.c + (b.warp.c - a.warp.c) * math::Complex(t, 0),
            a.warp.d + (b.warp.d - a.warp.d) * math::Complex(t, 0)};
  }

  /** @brief Draws the warped grid through column-time vertex captures.
   * @param canvas Destination frame.
   * @param start,end Retained animation endpoints. */
  HS_HOT_FLASH_MEMBER void
  draw_timed_grid(Canvas &canvas, const Endpoint &start, const Endpoint &end) {
    ScratchScope frame_scope(scratch_arena_b);
    auto *north = scratch_arena_b.allocate_n<NorthHole>(W / 2 + 1);
    auto *south = scratch_arena_b.allocate_n<SouthHole>(W / 2 + 1);
    const Spacetime::TemporalRotation rotation(start.rotation, end.rotation);
    for (int c = 0; c <= W / 2; ++c) {
      const float t = static_cast<float>(c) / (W / 2);
      const auto warp = warp_at(start, end, t);
      const auto q = rotation.sample(t);
      std::construct_at(&north[c], filters.template get<NorthHole>());
      std::construct_at(&south[c], filters.template get<SouthHole>());
      north[c].set_origin(math::normalized_or(
          math::rotate(math::mobius_transform(math::Y_AXIS, warp), q),
          math::X_AXIS));
      south[c].set_origin(math::normalized_or(
          math::rotate(math::mobius_transform(-math::Y_AXIS, warp), q),
          math::X_AXIS));
    }
    const auto difference = end.warp.b - start.warp.b;
    const float shift =
        sqrtf(difference.re * difference.re + difference.im * difference.im);
    const float reach = fmaxf(
        sqrtf(start.warp.b.re * start.warp.b.re +
              start.warp.b.im * start.warp.b.im),
        sqrtf(end.warp.b.re * end.warp.b.re + end.warp.b.im * end.warp.b.im));
    const float turn =
        2.0f * acosf(hs::clamp(fabsf(math::dot(start.rotation, end.rotation)),
                               0.0f, 1.0f));
    const auto warp_path =
        Spacetime::PathBound::warped_great_circle(0, shift, reach);
    const float speed = turn + 2.0f * shift;
    const Spacetime::PathBound bound{warp_path.second + turn * turn +
                                         4.0f * turn * shift,
                                     warp_path.third + speed * speed * speed +
                                         3.0f * turn * warp_path.second};
    const float displacement = turn + 2.0f * shift;
    const float phase_delta = std::remainder(end.phase - start.phase, 1.0f);
    auto &sink = filters.next.next;
    for (int axis = 0; axis < 2; ++axis) {
      const float first_num = axis == 0 ? start.rings : start.lines;
      const float num = axis == 0 ? end.rings : end.lines;
      const int count = static_cast<int>(ceilf(fmaxf(first_num, num)));
      for (int i = 0; i < count; ++i) {
        const bool had_curve = i < ceilf(first_num);
        const bool has_curve = i < ceilf(num);
        float first_index = static_cast<float>(i) / fmaxf(first_num, 1e-6f);
        float last_index = static_cast<float>(i) / fmaxf(num, 1e-6f);
        float first = axis == 0 ? math::wrap(first_index + start.phase, 1.0f)
                                : first_index;
        float last =
            axis == 0 ? math::wrap(last_index + end.phase, 1.0f) : last_index;
        if (!had_curve) {
          first = last;
          first_index = last_index;
        }
        if (!has_curve) {
          last = first;
          last_index = first_index;
        }
        const bool wraps = axis == 0 && fabsf(last - first) > 0.5f;
        if (wraps)
          last = first + std::remainder(last - first, 1.0f);
        const float initial = first, final = last;
        for (int part = 0; part < (wraps ? 2 : 1); ++part) {
          ScratchScope curve_scope(scratch_arena_a);
          const float offset = part == 0         ? 0.0f
                               : final > initial ? -1.0f
                                                 : 1.0f;
          first = initial + offset;
          last = final + offset;
          const float source_rate =
              (last - first) *
              (axis == 0 ? CONFORMAL_LOG_MAX - CONFORMAL_LOG_MIN : math::PI_F);
          const auto source_path = Spacetime::PathBound::warped_great_circle(
              source_rate, shift, reach);
          const float rate =
              turn + fabsf(source_rate) * (2.0f + reach * reach) + 2.0f * shift;
          const Spacetime::PathBound curve_bound{
              bound.second + 2.0f * source_path.second,
              bound.third + 4.0f * source_path.third + rate * rate * rate};
          const float curve_displacement =
              displacement + fabsf(source_rate) * (2.0f + reach * reach);
          auto ring_basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
          if (axis == 0 && CONFORMAL_LOG_MIN + last * (CONFORMAL_LOG_MAX -
                                                       CONFORMAL_LOG_MIN) <
                               0.0f)
            ring_basis.u = -ring_basis.u;
          auto evaluate = [&](size_t vertex, float t) HS_HOT_INLINE {
            const float coordinate = first + (last - first) * t;
            const float theta = (static_cast<float>(vertex) * 2.0f + 1.0f) *
                                math::PI_F / CURVE_SIDES;
            math::Vector source;
            if (axis == 0) {
              const float radius =
                  2.0f * atanf(expf(-CONFORMAL_LOG_MIN -
                                    coordinate * (CONFORMAL_LOG_MAX -
                                                  CONFORMAL_LOG_MIN)));
              source =
                  ring_basis.v * cosf(radius) +
                  (ring_basis.u * cosf(theta) + ring_basis.w * sinf(theta)) *
                      sinf(radius);
            } else {
              const float angle = coordinate * math::PI_F;
              const math::Vector normal(cosf(angle), 0.0f, -sinf(angle));
              source = math::cross(normal, math::Y_AXIS) * cosf(theta) +
                       math::Y_AXIS * sinf(theta);
            }
            return math::normalized_or(
                math::rotate(math::mobius_transform(source.normalized(),
                                                    warp_at(start, end, t)),
                             rotation.sample(t)),
                math::X_AXIS);
          };
          struct View {
            std::span<const math::Vector> vertices;
          };
          auto *vertices =
              scratch_arena_a.allocate_n<math::Vector>(CURVE_SIDES);
          ArenaVector<Plot::Mesh::Edge> edges;
          edges.bind(scratch_arena_a, CURVE_SIDES);
          for (int k = 0; k < CURVE_SIDES; ++k) {
            vertices[k] = evaluate(k, 1.0f);
            edges.push_back({static_cast<uint16_t>(k),
                             static_cast<uint16_t>((k + 1) % CURVE_SIDES)});
          }
          const float opacity =
              hs::clamp(num - static_cast<float>(i), 0.0f, 1.0f);
          const float first_opacity =
              hs::clamp(first_num - static_cast<float>(i), 0.0f, 1.0f);
          auto shader = [&](const math::Vector &v, Fragment &f,
                            float t) HS_HOT_INLINE {
            const float source_coordinate = first + (last - first) * t;
            if (axis == 0 &&
                (source_coordinate < 0.0f || source_coordinate >= 1.0f)) {
              f.color.alpha = 0.0f;
              return;
            }
            const float phase = start.phase + phase_delta * t;
            const float coordinate =
                axis == 0 ? hs::lerp(first_index, last_index, t)
                          : conformal_coord(
                                math::fast_sinf(((f.v2 + f.v0) / CURVE_SIDES +
                                                 0.5f / CURVE_SIDES) *
                                                2.0f * math::PI_F),
                                phase);
            f.color = baked_palette.get(coordinate);
            f.color.alpha *= hs::lerp(first_opacity, opacity, t) *
                             hs::lerp(start.alpha, end.alpha, t);
            const int column =
                std::clamp(static_cast<int>(t * (W / 2) + 0.5f), 0, W / 2);
            auto &north_hole = north[column];
            auto &south_hole = south[column];
            float alpha = 0.0f;
            north_hole.plot(
                v, f.color.color, f.age, f.color.alpha,
                [&](const math::Vector &p, const Pixel &color, float age,
                    float a) HS_HOT_INLINE {
                  south_hole.plot(
                      p, color, age, a,
                      [&](const math::Vector &, const Pixel &, float, float a2)
                          HS_HOT_INLINE { alpha = a2; });
                });
            f.color.alpha = alpha;
          };
          auto geometry = Plot::GeometryCallback::exact(evaluate);
          geometry.path = curve_bound;
          geometry.chained = true;
          geometry.rigid_lattice = true;
          Plot::Mesh::draw_column_edges<W, H>(
              sink, canvas, View{{vertices, CURVE_SIDES}}, edges, shader,
              geometry, curve_displacement);
        }
      }
    }
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
  BakedPaletteStorage
      baked_palette; /**< LUT-baked copy of `palette` the shaders sample. */
  /**
   * @brief User-tunable parameters.
   * @details num_rings/num_lines are animated counts driven by the Mutation
   *          timelines; alpha is the overall opacity multiplier. Declared before
   *          `timeline` so the counts outlive the Mutations that point here.
   */
  struct Params {
    /** Animated latitude-ring count; defaults to the ring Mutation's opening
     * sample. */
    float num_rings = 12.0f;
    /** Animated longitude-line count; defaults to the line Mutation's opening
     * sample. */
    float num_lines = 1.0f;
    float alpha = 0.2f; /**< Overall opacity multiplier in [0, 1]. */
  } params;

  /**
   * @brief Spinning render orientation.
   * @details Declared before `timeline` so it outlives the Rotation that points
   * here, which ~Timeline clears on teardown.
   */
  math::Orientation<> orientation;
  Timeline timeline; /**< Drives spin, palette wipe, and mutations. */
  MobiusWarpCircularTransformer<1> mobius_gen; /**< Möbius warp generator. */

  /** @brief Render filter pipeline: two hole fades, orientation, anti-alias. */
  Pipeline<W, H, NorthHole, SouthHole, Filter::World::Orient,
           Filter::Screen::AntiAlias<W, H>>
      filters;
};
