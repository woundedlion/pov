/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file Comets.h
 * @brief Comet heads tracing spherical Lissajous curves, dragging fading
 *        trails.
 */

#include <array>
#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"

// Unit-test accessor verifying each authored Lissajous entry closes
// (path_fn(domain) == path_fn(0)).
namespace hs_test {
namespace effects_tests {
struct CometsWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Comet whose head traces a spherical Lissajous curve, dragging a
 *        fading trail behind it.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details The path function and color palette periodically roll over to the
 *          next entry in the function table, cross-fading via a ColorWipe.
 * @note Sibling trail effects `Fishbowl` and `RingSpin` share the
 *       record + deep_tween skeleton and, with Fishbowl, the
 *       `Animation::TrailBody` aggregate at the same capacity and substep
 *       count. Draw primitive, transform chain and colour/fade are
 *       hand-propagated. Comets uses an empty pipeline (no Screen::AntiAlias)
 *       since Scan::Point glows carry their own softness.
 */
template <int W, int H> class Comets : public Effect {
public:
  static constexpr int TRAIL_LENGTH =
      115; /**< Number of past orientations retained in the comet trail. */
  static constexpr int ORIENTATION_SUBSTEPS =
      16; /**< Interpolation slots per Orientation, shared by the recorded trail
               and Motion. */

  /** @brief Comet head state: world orientation, recorded trail, body axis. */
  using Node = Animation::TrailBody<TRAIL_LENGTH, ORIENTATION_SUBSTEPS>;

  /**
   * @brief Constructs the effect at the templated canvas resolution.
   * @details Initializes the base Effect with the W x H dimensions and selects
   *          the first path/palette function table entry.
   */
  HS_COLD_MEMBER Comets()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})),
        palette(EffectPaletteRecipes::comets(
            EffectPaletteRecipes::random_base_turns())) {}

  /**
   * @brief Allocates state and wires up the animation timeline.
   * @details Allocates the node, bakes the palette LUT, registers params, and
   *          builds the timeline: an infinite RandomWalk + Motion + cycle
   *          timer, plus the periodic palette/path rollover.
   */
  void init() override {
    configure_presets(FUNCTIONS.size());
    node = persistent_arena.make<Node>();

    baked_palette.bake(persistent_arena, palette);

    register_param("Alpha", &params.alpha, 0.0f, 1.0f);
    register_param("Thickness", &params.thickness, 0.0f, 7.0f * THICKNESS_PX);
    register_param("Cycle Dur", &params.cycle_duration, 10.0f, 200.0f);
    register_param("Debug BB", &params.debug_bb);

    // Runs before motion exists, so its reanchor() is a no-op here; the path it
    // sets is still live because Motion below captures `path` by reference.
    update_path();
    timeline.add(0,
                 Animation::RandomWalk<W>(orientation, random_vector(), noise));
    // Motion + cycle timer are infinite and added before any finite animation,
    // so the timeline never relocates them and the retained handles stay valid.
    motion = timeline.add_get(
        0,
        Animation::Motion<W, ORIENTATION_SUBSTEPS>(
            node->orientation, path, (int)params.cycle_duration, true),
        Timeline::Pin::PINNED);
    cycle_timer = timeline.add_get(0,
                                   Animation::PeriodicTimer(
                                       2 * (int)params.cycle_duration,
                                       [this](Canvas &) {
                                         const bool advanced = advancePreset();
                                         HS_CHECK(advanced);
                                         update_palette();
                                       },
                                       true),
                                   Timeline::Pin::PINNED, &anims_paused);
  }

  /**
   * @brief Advances and renders one frame of the comet.
   * @details Steps the timeline, live-applies Cycle Dur, rebakes the palette
   *          while a wipe is in flight, records the trail, and draws the comet
   *          body along the trail.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(cm_timeline_step);
      timeline.step(canvas);
    }

    apply_if_changed((int)params.cycle_duration, last_cycle_dur, [&](int cd) {
      if (motion)
        motion->set_duration(cd);
      if (cycle_timer)
        cycle_timer->set_period(2 * cd);
    });

    {
      HS_PROFILE(cm_wipe_rebake);
      wipe.step(baked_palette, palette);
    }

    node->trail.record(node->orientation);

    // A sub-LSB alpha paints nothing; skip rasterizing. The trail is still
    // recorded above, so motion resumes when alpha rises.
    if (params.alpha < MIN_VISIBLE_ALPHA)
      return;

    HS_PROFILE(cm_draw_trail);
    deep_tween(node->trail, [&](const Quaternion &q, float t) {
      auto fragment_shader = [&](const Vector &, Fragment &f) {
        f.color = baked_palette.get(t);
        f.color.alpha *= quintic_kernel(t) * params.alpha;
      };

      Vector v_local = rotate(node->v, q);
      Vector v_final = orientation.orient(v_local);
      HS_PROFILE_DEEP(cm_point_scan);
      Scan::Point::draw<W, H>(filters, canvas, v_final, params.thickness,
                              fragment_shader, params.debug_bb);
    });
  }

private:
  HS_FLASH_MEMBER bool apply_preset(const PresetChange &change) override {
    if (!functions.select(change.to))
      return false;
    update_path();
    if (change.origin == PresetChangeOrigin::MANUAL) {
      node->orientation.set(Quaternion());
      node->trail.clear();
      motion->rewind();
      motion->reanchor();
    }
    return true;
  }

  // Test seam: asserts the closing-loop invariant the smoke harness cannot
  // observe.
  friend struct ::hs_test::effects_tests::CometsWhiteBox;

  static constexpr float THICKNESS_PX = 2.0f * PI_F / W;

  /**
   * @brief Snaps an authored domain to the nearest length that closes the curve.
   * @param config The Lissajous parameters whose domain is being snapped.
   * @return The closing domain: lissajous(m1, m2, a, closed_domain) equals the
   *         t=0 start (0,1,0) up to float error.
   * @details A spherical Lissajous point returns to the t=0 start (0,1,0) only
   *          when m2*domain is an exact multiple of 2*PI. Authored domains miss
   *          that by up to ~1.4 deg; snapping is a <0.3% nudge. Floor the cycle
   *          count at 1 so m2*domain < PI does not round to 0 and freeze the head
   *          at path_fn(0) (the table is authored data that gets extended).
   */
  static float closing_domain(const LissajousParams &config) {
    HS_CHECK(config.m2 > 0,
             "Comets Lissajous entry needs m2 > 0; m2 divides the domain");
    float closing_cycles = std::round(config.m2 * config.domain / (2 * PI_F));
    if (closing_cycles < 1.0f)
      closing_cycles = 1.0f;
    return 2 * PI_F * closing_cycles / config.m2;
  }

  /**
   * @brief Rebuilds the path function from the current function table entry.
   * @details Snaps the traversal length so the spherical Lissajous curve
   *          closes exactly, keeping the trace continuous across loops and
   *          function switches.
   */
  void update_path() {
    LissajousParams config = functions.get();
    // Snap so path_fn(domain) == path_fn(0); an unclosed endpoint pinches the
    // curve to a stray point each cycle.
    float closed_domain = closing_domain(config);
    // Capture only the three scalars + closed_domain (16 B): the whole
    // LissajousParams (20 B) overflows PlotFn's Fn<Vector(float), 16> inline
    // capacity (no heap fallback on Arduino).
    const float m1 = config.m1, m2 = config.m2, a = config.a;
    path.f = [m1, m2, a, closed_domain](float t) {
      return lissajous(m1, m2, a, t * closed_domain);
    };
    // Re-anchor Motion's baseline to the freshly-swapped path: the two curves'
    // travel-tangent frames differ at the seam, so a missing re-anchor teleports
    // the head for one frame.
    if (motion)
      motion->reanchor();
  }

  /**
   * @brief Rolls the palette over to a freshly generated one via a ColorWipe.
   * @details Arms the rebake gate for the wipe's duration and skips the
   *          rollover while a previous wipe is still in flight.
   */
  void update_palette() {
    // Skip while a wipe is in flight: at the Cycle Dur floor the cycle period
    // (20) is shorter than WIPE_FRAMES (48), so the timer can fire mid-wipe and a
    // second wipe would clobber the snapshots the live one still references.
    if (wipe.in_flight())
      return;
    wipe.arm(palette,
             GenerativePalette{EffectPaletteRecipes::comets(
                                   EffectPaletteRecipes::random_base_turns())}
                 .snapshot(),
             WIPE_FRAMES);
    timeline.add(0, Animation::ColorWipe(palette, wipe.start, wipe.target,
                                         WIPE_FRAMES, ease_linear));
  }

  static constexpr int WIPE_FRAMES =
      48; /**< Duration of a palette cross-fade ColorWipe, in frames. */

  // init() allocates the comet Node (holds the OrientationTrail) and one baked
  // palette LUT from the persistent arena.
  static constexpr size_t FOOTPRINT_BYTES =
      BakedPalette::required_arena_bytes() + sizeof(Node);
  // Effect keeps the default arena split, so the footprint must fit the device
  // persistent partition. Guards a TRAIL_LENGTH retune.
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "Comets persistent footprint exceeds the default partition; "
                "retune TRAIL_LENGTH or carve arenas");

  FastNoiseLite noise; /**< Noise source driving the head's RandomWalk. */
  Timeline timeline; /**< Animation timeline owning all scheduled animations. */
  Pipeline<W, H>
      filters; /**< Render filter pipeline applied to drawn fragments. */
  ProceduralPath path; /**< Current path function the comet head traces. */
  Orientation<> orientation; /**< World orientation walked by the RandomWalk. */
  GenerativePalette
      palette; /**< Active color palette (mutated by an in-flight ColorWipe). */
  BakedPalette
      baked_palette; /**< LUT-baked copy of `palette` sampled by the shader. */
  /** @brief Authored Lissajous parameter table backing `functions`.
   *  @details Each row is a LissajousParams {m1, m2, a, domain}: m1 axial (X/Z)
   *           frequency, m2 orbital (Y) frequency, a phase shift in radians,
   *           domain the traversal length t (closing_domain() snaps it so the
   *           curve closes). */
  static constexpr std::array<PresetEntry<LissajousParams>, 12> FUNCTIONS = {
      {// {m1, m2, a, domain}
       {{1.06f, 1.06f, 0, 5.909f}},
       {{6.06f, 1.0f, 0, 2 * PI_F}},
       {{6.02f, 4.01f, 0, 3.132f}},
       {{46.62f, 62.16f, 0, 0.404f}},
       {{46.26f, 69.39f, 0, 0.272f}},
       {{19.44f, 9.72f, 0, 0.646f}},
       {{8.51f, 17.01f, 0, 0.739f}},
       {{7.66f, 6.38f, 0, 4.924f}},
       {{8.75f, 5.0f, 0, 5.027f}},
       {{11.67f, 14.58f, 0, 2.154f}},
       {{11.67f, 8.75f, 0, 2.154f}},
       {{10.94f, 8.75f, 0, 2.872f}}}};
  /** @brief Cyclic selector over FUNCTIONS for the active path/palette entry. */
  Presets<LissajousParams, 12> functions{FUNCTIONS};
  Node *node = nullptr; /**< Arena-allocated comet head state. */
  PaletteWipe wipe;     /**< Cross-fade state of the palette rollover. */
  Animation::Motion<W, ORIENTATION_SUBSTEPS> *motion =
      nullptr; /**< Handle to the infinite Motion driving the head along `path`. */
  Animation::PeriodicTimer *cycle_timer =
      nullptr; /**< Handle to the timer that rolls path/palette over. */
  int last_cycle_dur =
      -1; /**< Last applied Cycle Dur, in frames; -1 forces a first apply. */

  /**
   * @brief User-tunable parameters exposed as effect sliders.
   */
  struct Params {
    float alpha = 1.0f; /**< Overall trail opacity multiplier in [0, 1]. */
    float thickness = 2.1f * THICKNESS_PX; /**< Comet body half-width. */
    float cycle_duration =
        80.0f; /**< Duration of one motion cycle, in frames. */
    bool debug_bb =
        false; /**< When true, draws the fragment bounding box for debugging. */
  } params;    /**< Live parameter block bound to the registered sliders. */
};

#include "core/control/registry.h"
REGISTER_EFFECT(Comets)
