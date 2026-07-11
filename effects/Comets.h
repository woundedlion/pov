/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
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
 * @note Sibling trail effects `ChaoticStrings` and `RingSpin` share only the
 *       record + deep_tween skeleton; their draw/transform/fade diverge, so
 *       trail fixes must be propagated by hand. Comets uses an empty pipeline
 *       (no Screen::AntiAlias) since Scan::Point glows carry their own softness.
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
  FLASHMEM Comets() : Effect(W, H, {.strobe = true}), cur_function_idx(0) {}

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

    baked_palette.bake(persistent_arena, palette);

    register_param("Alpha", &params.alpha, 0.0f, 1.0f);
    register_param("Thickness", &params.thickness, 0.0f, 0.5f);
    register_param("Cycle Dur", &params.cycle_duration, 10.0f, 200.0f);
    register_param("Debug BB", &params.debug_bb);

    // Runs before motion_ exists, so its reanchor() is a no-op here; the path it
    // sets is still live because Motion below captures `path` by reference.
    update_path();
    timeline.add(0,
                 Animation::RandomWalk<W>(orientation, random_vector(), noise));
    // Motion + cycle timer are infinite and added before any finite animation,
    // so the timeline never relocates them and the retained handles stay valid.
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
   * @brief Advances and renders one frame of the comet.
   * @details Steps the timeline, live-applies Cycle Dur, rebakes the palette
   *          while a wipe is in flight, records the trail, and draws the comet
   *          body along the trail.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    apply_if_changed((int)params.cycle_duration, last_cycle_dur_, [&](int cd) {
      if (motion_)
        motion_->set_duration(cd);
      if (cycle_timer_)
        cycle_timer_->set_period(2 * cd);
    });

    // `palette` changes only while a ColorWipe steps. The wipe is armed mid-step
    // and first steps next frame, so skip the redundant rebake on the arming
    // frame.
    if (wipe_pending_) {
      wipe_pending_ = false;
    } else if (wipe_frames_remaining_ > 0) {
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
  // Test seam: asserts the closing-loop invariant the smoke harness cannot
  // observe.
  friend struct ::hs_test::effects_tests::CometsWhiteBox;

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
    LissajousParams config = functions[cur_function_idx];
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
    if (motion_)
      motion_->reanchor();
  }

  /**
   * @brief Rolls the palette over to a freshly generated one via a ColorWipe.
   * @details Arms the rebake gate for the wipe's duration and skips the
   *          rollover while a previous wipe is still in flight.
   */
  void update_palette() {
    // Skip while a wipe is in flight: at the Cycle Dur floor the cycle period
    // (20) is shorter than WIPE_FRAMES (48), so the timer can fire mid-wipe and a
    // second wipe would clobber the next_palette_ the live one still references.
    if (wipe_frames_remaining_ > 0)
      return;
    next_palette_ =
        GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                          BrightnessProfile::ASCENDING);
    timeline.add(0, Animation::ColorWipe(palette, next_palette_, WIPE_FRAMES,
                                         ease_linear));
    wipe_frames_remaining_ = WIPE_FRAMES;
    wipe_pending_ = true;
  }

  static constexpr int WIPE_FRAMES = 48; /**< Duration of a palette cross-fade ColorWipe, in frames. */

  // init() allocates the comet Node (holds the OrientationTrail) and one baked
  // palette LUT from the persistent arena.
  static constexpr size_t FOOTPRINT_BYTES =
      BakedPalette::LUT_SIZE * sizeof(Color4) + sizeof(Node);
  // Effect keeps the default arena split, so the footprint must fit the device
  // persistent partition. Guards a TRAIL_LENGTH retune.
  static constexpr size_t PERSISTENT_BUDGET =
      DEVICE_GLOBAL_ARENA_SIZE - DEFAULT_SCRATCH_A_SIZE - DEFAULT_SCRATCH_B_SIZE;
  static_assert(FOOTPRINT_BYTES <= PERSISTENT_BUDGET,
                "Comets persistent footprint exceeds the default partition; "
                "retune TRAIL_LENGTH or carve arenas");

  FastNoiseLite noise; /**< Noise source driving the head's RandomWalk. */
  Timeline timeline; /**< Animation timeline owning all scheduled animations. */
  Pipeline<W, H> filters; /**< Render filter pipeline applied to drawn fragments. */
  ProceduralPath path; /**< Current path function the comet head traces. */
  Orientation<> orientation; /**< World orientation walked by the RandomWalk. */
  GenerativePalette palette; /**< Active color palette (mutated by an in-flight ColorWipe). */
  BakedPalette baked_palette; /**< LUT-baked copy of `palette` sampled by the shader. */
  /** @brief Function table of Lissajous parameters cycled through over time.
   *  @details Each row is a LissajousParams {m1, m2, a, domain}: m1 axial (X/Z)
   *           frequency, m2 orbital (Y) frequency, a phase shift in radians,
   *           domain the traversal length t (closing_domain() snaps it so the
   *           curve closes). */
  static constexpr std::array<LissajousParams, 12> functions = {{// {m1, m2, a, domain}
                                                {1.06f, 1.06f, 0, 5.909f},
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
  bool wipe_pending_ = false; /**< Wipe armed this frame; it first steps next frame. */

  /**
   * @brief User-tunable parameters exposed as effect sliders.
   */
  struct Params {
    float alpha = 1.0f; /**< Overall trail opacity multiplier in [0, 1]. */
    float thickness = 2.1f * 2 * PI_F / W; /**< Comet body half-width, in radians; default scaled to ≈2.1 px at this build's W for legibility. */
    float cycle_duration = 80.0f; /**< Duration of one motion cycle, in frames. */
    bool debug_bb = false; /**< When true, draws the fragment bounding box for debugging. */
  } params; /**< Live parameter block bound to the registered sliders. */
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(Comets)
