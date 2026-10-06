/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file Comets.h
 * @brief A comet head tracing spherical Lissajous curves, dragging a fading
 *        trail.
 */

#include "core/animation/orientation.h"
#include <array>
#include "core/color/effect_palette_recipes.h"
#include "core/control/choreography.h"
#include "core/engine/engine.h"

namespace hs_test {
namespace effects_tests {
struct CometsWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/** @brief Comets' live parameter set: the preset-driven Lissajous function
 *  plus the slider-bound render values. */
struct CometsParams {
  math::LissajousParams function; /**< Path function the comet head traces. */
  float alpha = 1.0f;     /**< Overall trail opacity multiplier in [0, 1]. */
  float thickness = 0.0f; /**< Comet body half-width; initial_params() seeds the
                               resolution-derived default. */
  float cycle_duration = 80.0f; /**< Duration of one motion cycle, in frames. */
  bool debug_bb = false; /**< When true, draws the fragment bounding box. */
};

/**
 * @brief Comet whose head traces a spherical Lissajous curve, dragging a
 *        fading trail behind it.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Automatic preset changes snap to the next Lissajous entry and
 *          ColorWipe to a freshly generated palette. Manual selection restarts
 *          the path and keeps the palette.
 */
template <int W, int H>
class Comets : public ChoreographedEffect<Comets<W, H>, CometsParams> {
  using Choreography = ChoreographedEffect<Comets<W, H>, CometsParams>;
  friend Choreography;

public:
  static constexpr const char *EFFECT_ID = "Comets";

  using Params = CometsParams;

  static constexpr int TRAIL_LENGTH = Animation::
      TRAIL_HISTORY_LENGTH; /**< Number of past orientations retained in the comet trail. */
  static constexpr int ORIENTATION_SUBSTEPS = Animation::
      TRAIL_ORIENTATION_SUBSTEPS; /**< Interpolation slots per Orientation, shared by the recorded trail
               and Motion. */

  /** Snaps the path function; the palette rolls over via a ColorWipe. */
  static constexpr Segue::Preset::Snap DEPARTURE{};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  /** Preset cadence: two default-duration motion cycles. */
  static constexpr uint16_t PRESET_DWELL_FRAMES = 160;

  /** @brief Angular width of one canvas column, the comet thickness unit. */
  static constexpr float THICKNESS_PX = math::RADIANS_PER_COLUMN<W>;

  static constexpr float ALPHA_MIN = 0.0f, ALPHA_MAX = 1.0f;
  static constexpr float THICKNESS_MIN = 0.0f,
                         THICKNESS_MAX = 7.0f * THICKNESS_PX;
  static constexpr float CYCLE_DURATION_MIN = 10.0f,
                         CYCLE_DURATION_MAX = 200.0f;

  static Params initial_params() {
    return {.function = PRESETS[0].params,
            .alpha = 1.0f,
            .thickness = 2.1f * THICKNESS_PX,
            .cycle_duration = 80.0f,
            .debug_bb = false};
  }

  /** @brief Shared registration, validation and interpolation descriptions. */
  static constexpr auto parameter_fields() {
    return std::tuple{
        Control::Field<Params, float>{
            .id = "alpha",
            .member = &Params::alpha,
            .name = "Alpha",
            .spec = {.min = ALPHA_MIN, .max = ALPHA_MAX}},
        Control::Field<Params, float>{
            .id = "thickness",
            .member = &Params::thickness,
            .name = "Thickness",
            .spec = {.min = THICKNESS_MIN, .max = THICKNESS_MAX}},
        Control::Field<Params, float>{
            .id = "cycle_duration",
            .member = &Params::cycle_duration,
            .name = "Cycle Dur",
            .spec = {.min = CYCLE_DURATION_MIN, .max = CYCLE_DURATION_MAX}},
        Control::Field<Params, bool>{.id = "debug_bb",
                                     .member = &Params::debug_bb,
                                     .name = "Debug BB",
                                     .spec = {.min = 0, .max = 1}}};
  }

  static bool valid_params(const Params &p) {
    return std::isfinite(p.function.m1) && std::isfinite(p.function.m2) &&
           std::isfinite(p.function.a) && std::isfinite(p.function.domain) &&
           p.function.m2 > 0.0f && Control::valid_fields(p, parameter_fields());
  }

  /** @brief Comet head state: world orientation, recorded trail, body axis. */
  using Node = Animation::TrailBody<TRAIL_LENGTH, ORIENTATION_SUBSTEPS>;

  /**
   * @brief Constructs the effect at the templated canvas resolution.
   */
  HS_COLD_MEMBER Comets()
      : Choreography(W, H,
                     pipeline_config<decltype(filters)>({.strobe = true})),
        palette(EffectPaletteRecipes::comets(
            EffectPaletteRecipes::random_base_turns())) {}

  /**
   * @brief Allocates state and wires up the animation timeline.
   */
  HS_COLD_MEMBER void init() override {
    begin_choreography();
    node = persistent_arena.make<Node>();

    baked_palette.bake(persistent_arena, palette);

    this->register_described_params();

    // motion is null here; Motion captures `path` by reference.
    update_path();
    timeline.add(
        0, Animation::RandomWalk<W>(orientation, math::random_vector(), noise));
    // Infinite and pinned, so the retained handle stays valid.
    motion = timeline.add_get(
        0,
        Animation::Motion<W, ORIENTATION_SUBSTEPS>(
            node->orientation, path, (int)params.cycle_duration, true),
        Timeline::Pin::PINNED);
  }

  /**
   * @brief Advances and renders one frame of the comet.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(cm_timeline_step);
      timeline.step(canvas);
    }
    step_choreography();

    apply_if_changed((int)params.cycle_duration, last_cycle_dur, [&](int cd) {
      if (motion)
        motion->set_duration(cd);
    });

    {
      HS_PROFILE(cm_wipe_rebake);
      wipe.step(baked_palette, palette);
    }

    node->trail.record(node->orientation);

    // Alpha below one slider LSB: skip rasterizing; the trail still records.
    if (params.alpha < MIN_VISIBLE_ALPHA)
      return;

    HS_PROFILE(cm_draw_trail);
    deep_tween(node->trail, [&](const math::Quaternion &q, float t) {
      auto fragment_shader = [&](const math::Vector &, Fragment &f) {
        f.color = baked_palette.get(t);
        f.color.alpha *= math::quintic_kernel(t) * params.alpha;
      };

      math::Vector v_local = math::rotate(node->v, q);
      math::Vector v_final = orientation.orient(v_local);
      HS_PROFILE_DEEP(cm_point_scan);
      Scan::Point::draw<W, H>(filters, canvas, v_final, params.thickness,
                              fragment_shader, params.debug_bb);
    });
  }

private:
  using Choreography::begin_choreography;
  using Choreography::params;
  using Choreography::register_param;
  using Choreography::step_choreography;
  using Choreography::timeline;

  /** @brief Params for preset @p index: the live slider values with the
   *  table's function patched in, so a rollover swaps only the path. */
  Params preset_params(size_t index) {
    Params p = params;
    p.function = PRESETS[index].params;
    return p;
  }

  /** @brief Adopts a snap or snapshot target and rebuilds the path from it. */
  void adopt_params(const Params &target) {
    params = target;
    update_path();
  }

  /** @brief A manual selection restarts the authored path deterministically. */
  HS_FLASH_MEMBER bool
  apply_preset(const Effect::PresetChange &change) override {
    if (!Choreography::apply_preset(change))
      return false;
    if (change.origin == Effect::PresetChangeOrigin::MANUAL) {
      node->orientation.set(math::Quaternion());
      node->trail.clear();
      if (motion) {
        motion->rewind();
        motion->reanchor();
      }
    }
    return true;
  }

  /** @brief The automatic cadence also rolls the palette; manual selection and
   *  snapshot restores keep the live palette. */
  HS_FLASH_MEMBER void
  preset_changed(const Effect::PresetChange &change) override {
    if (change.origin == Effect::PresetChangeOrigin::AUTOMATIC)
      update_palette();
  }

  friend struct ::hs_test::effects_tests::CometsWhiteBox;

  /**
   * @brief Snaps an authored domain to the nearest length that closes the curve.
   * @param config The Lissajous parameters whose domain is being snapped.
   * @return The closing domain: lissajous(m1, m2, a, closed_domain) equals the
   *         t=0 start (0,1,0) up to float error.
   * @details The curve returns to (0,1,0) only when m2*domain is a multiple
   *          of 2*PI. The cycle count floors at 1 so m2*domain < PI does not
   *          freeze the head at path_fn(0).
   */
  static float closing_domain(const math::LissajousParams &config) {
    HS_CHECK(config.m2 > 0,
             "Comets Lissajous entry needs m2 > 0; m2 divides the domain");
    float closing_cycles =
        std::round(config.m2 * config.domain / (2 * math::PI_F));
    if (closing_cycles < 1.0f)
      closing_cycles = 1.0f;
    return 2 * math::PI_F * closing_cycles / config.m2;
  }

  /**
   * @brief Rebuilds the path function from the live function params.
   * @details Snaps the traversal length so the spherical Lissajous curve
   *          closes exactly, keeping the trace continuous across loops and
   *          function switches.
   */
  void update_path() {
    math::LissajousParams config = params.function;
    float closed_domain = closing_domain(config);
    // Four scalars fill PlotFn's 16 B inline capacity (no heap fallback).
    const float m1 = config.m1, m2 = config.m2, a = config.a;
    path.f = [m1, m2, a, closed_domain](float t) {
      return math::lissajous(m1, m2, a, t * closed_domain);
    };
    // Without a re-anchor the head teleports for one frame at the path swap.
    if (motion)
      motion->reanchor();
  }

  /**
   * @brief Rolls the palette over to a freshly generated one via a ColorWipe.
   * @details Arms the rebake gate for the wipe's duration and skips the
   *          rollover while a previous wipe is still in flight.
   */
  void update_palette() {
    // A second wipe would clobber the snapshots the live one references.
    if (wipe.in_flight())
      return;
    wipe.arm(palette,
             GenerativePalette{EffectPaletteRecipes::comets(
                                   EffectPaletteRecipes::random_base_turns())}
                 .snapshot(),
             WIPE_FRAMES);
    timeline.add(0, Animation::ColorWipe(palette, wipe.start, wipe.target,
                                         WIPE_FRAMES, math::ease_linear));
  }

  static constexpr int WIPE_FRAMES =
      48; /**< Duration of a palette cross-fade ColorWipe, in frames. */

  // Persistent allocations: the comet Node and the baked palette LUT.
  static constexpr size_t FOOTPRINT_BYTES =
      BakedPalette::required_arena_bytes() + sizeof(Node);
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "Comets persistent footprint exceeds the default partition; "
                "retune TRAIL_LENGTH or carve arenas");

  FastNoiseLite noise; /**< Noise source driving the head's RandomWalk. */
  Pipeline<W, H>
      filters; /**< Render filter pipeline applied to drawn fragments. */
  ProceduralPath path; /**< Current path function the comet head traces. */
  math::Orientation<>
      orientation; /**< World orientation walked by the RandomWalk. */
  GenerativePalette
      palette; /**< Active color palette (mutated by an in-flight ColorWipe). */
  BakedPaletteStorage
      baked_palette; /**< LUT-baked copy of `palette` sampled by the shader. */
  /** @brief Authored Lissajous preset table; preset_params() patches each
   *  entry into the live parameter set.
   *  @details Each row is a LissajousParams {m1, m2, a, domain}: m1 axial (X/Z)
   *           frequency, m2 orbital (Y) frequency, a phase shift in radians,
   *           domain the traversal length t (closing_domain() snaps it so the
   *           curve closes). */
  static constexpr std::array<PresetEntry<math::LissajousParams>, 12> PRESETS =
      {{// {m1, m2, a, domain}
        {{1.06f, 1.06f, 0, 5.909f}, DEPARTURE},
        {{6.06f, 1.0f, 0, 2 * math::PI_F}, DEPARTURE},
        {{6.02f, 4.01f, 0, 3.132f}, DEPARTURE},
        {{46.62f, 62.16f, 0, 0.404f}, DEPARTURE},
        {{46.26f, 69.39f, 0, 0.272f}, DEPARTURE},
        {{19.44f, 9.72f, 0, 0.646f}, DEPARTURE},
        {{8.51f, 17.01f, 0, 0.739f}, DEPARTURE},
        {{7.66f, 6.38f, 0, 4.924f}, DEPARTURE},
        {{8.75f, 5.0f, 0, 5.027f}, DEPARTURE},
        {{11.67f, 14.58f, 0, 2.154f}, DEPARTURE},
        {{11.67f, 8.75f, 0, 2.154f}, DEPARTURE},
        {{10.94f, 8.75f, 0, 2.872f}, DEPARTURE}}};
  Node *node = nullptr; /**< Arena-allocated comet head state. */
  PaletteWipe wipe;     /**< Cross-fade state of the palette rollover. */
  Animation::Motion<W, ORIENTATION_SUBSTEPS> *motion =
      nullptr; /**< Handle to the infinite Motion driving the head along `path`. */
  int last_cycle_dur =
      -1; /**< Last applied Cycle Dur, in frames; -1 forces a first apply. */
};
