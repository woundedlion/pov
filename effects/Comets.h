/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include "core/engine.h"

// Forward declaration of the unit-test accessor (tests/test_effects.h) that
// reaches closing_domain and the function table to verify each authored
// Lissajous entry closes (path_fn(domain) == path_fn(0)).
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
 * @note Sibling trail effects — `ChaoticStrings` and `RingSpin` — share the same
 *       scaffolding (orientation random-walk → motion-blur/fading trail →
 *       position-colored multiline). Only the record + deep_tween skeleton is
 *       genuinely common, and it already lives in the engine
 *       (`OrientationTrail::record` + the `deep_tween` free function); the
 *       per-effect bodies diverge in draw primitive, transform chain, color/fade,
 *       and accumulate-vs-draw model, so a unifying base would have to
 *       parameterize all of those — a worse abstraction than the two-line idiom
 *       and a fresh coupling across three independently-tuned effects. They are
 *       therefore deliberately not unified; a trail-rendering fix must be
 *       propagated by hand across all three. Known divergences:
 *         - `RingSpin` omits the `Screen::AntiAlias` filter this effect and
 *           `ChaoticStrings` apply.
 *         - `RingSpin` uses `Orientation<>` (CAP 4) where this effect and
 *           `ChaoticStrings` use `Orientation<16>` (up to 4 motion-blur
 *           sub-frames versus 16). Intentional: its full great-circle rings
 *           overlap heavily frame-to-frame, so 4 sub-frames read identically to
 *           16, while a single fast-moving point/string needs the finer fidelity.
 */
template <int W, int H> class Comets : public Effect {
public:
  static constexpr int TRAIL_LENGTH = 115; /**< Number of past orientations retained in the comet trail. */

  /**
   * @brief Comet head state: world orientation, recorded trail, and body axis.
   * @details Holds the head's world orientation, the recorded trail of past
   *          orientations, and the local direction vector being drawn (the
   *          body axis).
   */
  struct Node {
    Orientation<16> orientation; /**< Current world orientation of the comet head. */
    Animation::OrientationTrail<Orientation<16>, TRAIL_LENGTH> trail; /**< Recorded trail of past orientations. */
    Vector v; /**< Local direction vector drawn as the comet body axis. */

    /**
     * @brief Constructs a node with its body axis aligned to the Y axis.
     */
    Node() : v(Y_AXIS) {}
  };

  /**
   * @brief Constructs the effect at the templated canvas resolution.
   * @details Initializes the base Effect with the W x H dimensions and selects
   *          the first path/palette function table entry.
   */
  FLASHMEM Comets() : Effect(W, H), cur_function_idx(0) {}

  /**
   * @brief Allocates state and wires up the animation timeline.
   * @details Allocates the node, bakes the palette LUT, registers params, and
   *          builds the timeline: an infinite RandomWalk + Motion + cycle
   *          timer, plus the periodic palette/path rollover.
   */
  void init() override {
    node = static_cast<Node *>(
        persistent_arena.allocate(sizeof(Node), alignof(Node)));
    new (node) Node();

    // Allocate baked palette LUT (rebaked only while a ColorWipe is animating
    // `palette`; see draw_frame()).
    baked_palette.bake(persistent_arena, palette);

    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    // Thickness is an absolute angular half-width in radians, so this [0, 0.5]
    // rad range is intentionally resolution-independent: the slider maximum is a
    // fixed fraction (~0.5/2π ≈ 8%) of the full 2π ring at any W. (The *default*
    // below is instead tuned to a fixed pixel size for legibility; the two need
    // not share a W-dependence — one fixes the body's pixel look, the other fixes
    // its maximum angular share of the circle.)
    registerParam("Thickness", &params.thickness, 0.0f, 0.5f);
    registerParam("Cycle Dur", &params.cycle_duration, 10.0f, 200.0f);
    registerParam("Debug BB", &params.debug_bb);

    update_path();
    timeline.add(0,
                 Animation::RandomWalk<W>(orientation, random_vector(), noise));
    // Motion + cycle timer are infinite and added before the only finite
    // animation (the periodic ColorWipe, always appended after), so the
    // timeline never moves them — these add_get() handles stay valid for live
    // Cycle Dur updates.
    motion_ = timeline.add_get(
        0, Animation::Motion<W, 16>(node->orientation, path,
                                    (int)params.cycle_duration, true));
    cycle_timer_ = timeline.add_get(
        0, Animation::PeriodicTimer(
               2 * (int)params.cycle_duration,
               [this](Canvas &) {
                 cur_function_idx = (cur_function_idx + 1) % functions.size();
                 update_path();
                 update_palette();
               },
               true));
  }

  /**
   * @brief Reports whether the engine should clear to the background each frame.
   * @return Always false; this effect manages its own framebuffer contents.
   */
  bool show_bg() const override { return false; }

  /**
   * @brief Advances and renders one frame of the comet.
   * @details Steps the timeline, live-applies Cycle Dur, rebakes the palette
   *          while a wipe is in flight, records the trail, and draws the comet
   *          body along the trail.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    // Live-apply the Cycle Dur slider to the motion duration + cycle timer.
    apply_if_changed((int)params.cycle_duration, last_cycle_dur_, [&](int cd) {
      if (motion_)
        motion_->set_duration(cd);
      if (cycle_timer_)
        cycle_timer_->set_period(2 * cd);
    });

    // `palette` only changes while a ColorWipe is mutating it (armed by
    // update_palette() when the cycle timer rolls the palette over); the rest of
    // the time it's static, so skip the 256-entry rebake instead of redoing it
    // every frame.
    if (wipe_frames_remaining_ > 0) {
      baked_palette.rebake(palette);
      --wipe_frames_remaining_;
    }

    node->trail.record(node->orientation);

    deep_tween(node->trail, [&](const Quaternion &q, float t) {
      auto fragment_shader = [&](const Vector &, Fragment &f) {
        f.color = baked_palette.get(t);
        f.color.alpha *= quintic_kernel(t) * params.alpha;
      };

      Vector v_local = rotate(node->v, q);
      Vector v_final = orientation.orient(v_local);
      Scan::Point::draw<W, H>(filters, canvas, v_final, params.thickness,
                              fragment_shader, params.debug_bb);
    });
  }

private:
  // Test seam: reaches the closing_domain snap and the function table so the
  // unit test can assert path_fn(domain) == path_fn(0) for every entry — the
  // closing-loop invariant the smoke harness cannot observe.
  friend struct ::hs_test::effects_tests::CometsWhiteBox;

  /**
   * @brief Snaps an authored domain to the nearest length that closes the curve.
   * @param config The Lissajous parameters whose domain is being snapped.
   * @return The closing domain: lissajous(m1, m2, a, closed_domain) equals the
   *         t=0 start (0,1,0) up to float error.
   * @details A spherical Lissajous point is (sin(m2 t)*cos(...), cos(m2 t),
   *          sin(m2 t)*sin(...)), so it returns to the t=0 start (0,1,0) only when
   *          sin(m2 t)=0 and cos(m2 t)=1, i.e. when m2*domain is an exact multiple
   *          of 2*PI (the m1/phase terms drop out once sin(m2 t)=0). The authored
   *          domains miss that by up to ~1.4 deg. Snapping to the nearest closing
   *          length is a <0.3% domain nudge, invisible to the shape. Floor the
   *          cycle count at 1: when m2*domain < PI the round() collapses to 0,
   *          which would zero the result and freeze the head at path_fn(0). All 12
   *          current entries clear the threshold, but the table is authored data
   *          that gets extended.
   */
  static float closing_domain(const LissajousParams &config) {
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
    LissajousParams config = functions[cur_function_idx];
    // Snap the traversal length so the curve closes exactly (closing_domain
    // carries the geometry). Motion advances the head from path_fn(0) and then
    // hard-resets to the start anchor every cycle (drift correction); when the
    // endpoint doesn't coincide with the start, that reset teleports the head by
    // the gap — a small jump that reads as a path discontinuity, most visibly
    // right when the function switches. The snap makes path_fn(domain) ==
    // path_fn(0) so the reset is a no-op and the trace stays continuous across
    // loops and switches.
    float closed_domain = closing_domain(config);
    // Capture only the three scalars lissajous() needs plus closed_domain (16 B
    // total). Capturing the whole LissajousParams (16 B) + closed_domain (4 B)
    // is 20 B, overflowing PlotFn's Fn<Vector(float), 16> inline capacity — a
    // hard static_assert on the ARDUINO build, where Fn is teensy's
    // inplace_function with no heap fallback (host/WASM std::function masks it).
    const float m1 = config.m1, m2 = config.m2, a = config.a;
    path.f = [m1, m2, a, closed_domain](float t) {
      return lissajous(m1, m2, a, t * closed_domain);
    };
  }

  /**
   * @brief Rolls the palette over to a freshly generated one via a ColorWipe.
   * @details Arms the rebake gate for the wipe's duration and skips the
   *          rollover while a previous wipe is still in flight.
   */
  void update_palette() {
    // Skip the rollover while a wipe is still in flight. The ColorWipe reads
    // `next_palette_` as its target and mutates `palette` for WIPE_FRAMES frames;
    // at the Cycle Dur floor the cycle period (2*cycle_dur = 20) is shorter than
    // WIPE_FRAMES (48), so the timer can fire again mid-wipe. Starting a second
    // wipe would overwrite the next_palette_ the live wipe still references and
    // queue a second mutator fighting over `palette`. The next timer tick picks
    // up the rollover once the in-flight wipe drains.
    if (wipe_frames_remaining_ > 0)
      return;
    next_palette_ =
        GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                          BrightnessProfile::ASCENDING);
    timeline.add(0, Animation::ColorWipe(palette, next_palette_, WIPE_FRAMES,
                                         ease_mid));
    // Arm the rebake gate for the life of the wipe. The wipe is scheduled inside
    // timeline.step() and first steps next frame, mutating `palette` for
    // WIPE_FRAMES frames; +1 covers the arming frame so the count safely spans
    // every mutated frame (a redundant rebake of an unchanged palette is
    // harmless).
    wipe_frames_remaining_ = WIPE_FRAMES + 1;
  }

  static constexpr int WIPE_FRAMES = 48; /**< Duration of a palette cross-fade ColorWipe, in frames. */

  FastNoiseLite noise; /**< Noise source driving the head's RandomWalk. */
  Timeline timeline; /**< Animation timeline owning all scheduled animations. */
  // Empty pipeline by design: Comets deliberately omits Screen::AntiAlias. The
  // comet points are drawn via Scan::Point with their own thickness/softness and
  // overlapping trail alpha, so a screen-space AA pass would only blur the glow
  // without improving the look. Add AntiAlias here if crisp edges are ever wanted.
  Pipeline<W, H> filters; /**< Render filter pipeline applied to drawn fragments. */
  ProceduralPath path; /**< Current path function the comet head traces. */
  Orientation<> orientation; /**< World orientation walked by the RandomWalk. */
  GenerativePalette palette; /**< Active color palette (mutated by an in-flight ColorWipe). */
  BakedPalette baked_palette; /**< LUT-baked copy of `palette` sampled by the shader. */
  /** @brief Function table of Lissajous parameters cycled through over time.
   *  @details Immutable authored data, so it is `static constexpr` — shared
   *           across instances (no per-effect copy) and readable without an
   *           instance (the closing-loop unit test reads it directly). */
  static constexpr std::array<LissajousParams, 12> functions = {{{1.06f, 1.06f, 0, 5.909f},
                                                {6.06f, 1.0f, 0, 2 * PI_F},
                                                {6.02f, 4.01f, 0, 3.132f},
                                                {46.62f, 62.16f, 0, 0.404f},
                                                {46.26f, 69.39f, 0, 0.272f},
                                                {19.44f, 9.72f, 0, 0.646f},
                                                {8.51f, 17.01f, 0, 0.739f},
                                                {7.66f, 6.38f, 0, 4.924f},
                                                {8.75f, 5.0f, 0, 5.027f},
                                                {11.67f, 14.58f, 0, 2.154f},
                                                {11.67f, 8.75f, 0, 2.154f},
                                                {10.94f, 8.75f, 0, 2.872f}}};
  int cur_function_idx; /**< Index into `functions` of the active path/palette entry. */
  Node *node = nullptr; /**< Arena-allocated comet head state. */
  GenerativePalette next_palette_; /**< Target palette a ColorWipe fades toward. */
  Animation::Motion<W, 16> *motion_ = nullptr; /**< Handle to the infinite Motion driving the head along `path`. */
  Animation::PeriodicTimer *cycle_timer_ = nullptr; /**< Handle to the timer that rolls path/palette over. */
  int last_cycle_dur_ = -1; /**< Last applied Cycle Dur, in frames; -1 forces a first apply. */
  int wipe_frames_remaining_ = 0; /**< Frames left to rebake `palette` for the in-flight wipe. */

  /**
   * @brief User-tunable parameters exposed as effect sliders.
   */
  struct Params {
    float alpha = 1.0f; /**< Overall trail opacity multiplier in [0, 1]. */
    float thickness = 2.1f * 2 * PI_F / W; /**< Comet body half-width, in radians; default scaled to ≈2.1 px at this build's W for legibility (the Thickness slider's [0,0.5] rad range is, by contrast, an absolute angular span — see registerParam). */
    float cycle_duration = 80.0f; /**< Duration of one motion cycle, in frames. */
    bool debug_bb = false; /**< When true, draws the fragment bounding box for debugging. */
  } params; /**< Live parameter block bound to the registered sliders. */
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Comets)
