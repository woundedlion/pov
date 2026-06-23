/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"

// Forward declaration of the unit-test accessor (tests/test_effects.h) that
// reaches palette_boundaries to stage the documented overlapping-wipe band
// inversion and assert color() stays memory-safe and in-range.
namespace hs_test {
namespace effects_tests {
struct DynamoWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief A snaking strand of nodes pulled across the sphere, leaving fading
 *        trails, with palettes that sweep in via angular color wipes.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Speed, gap, trail length, and wipe duration are live sliders;
 *          direction, rotation, and wipes are driven by random timers.
 */
template <int W, int H> class Dynamo : public Effect {
public:
  /**
   * @brief One point on the strand: grid position (x,y) and per-step velocity.
   */
  struct Node {
    /**
     * @brief Constructs a node at the origin with zero velocity.
     */
    Node() : x(0), y(0), v(0) {}

    int x; /**< Grid column. */
    int y; /**< Grid row. */
    int v; /**< Per-step velocity along x. */
  };

  /**
   * @brief The effect's canonical generative-palette recipe.
   * @details The initial palette and every color-wipe palette share the same
   *          VIGNETTE / ANALOGOUS / ASCENDING recipe; centralized here so the
   *          three construction arguments are not duplicated across the
   *          constructor and color_wipe().
   */
  static GenerativePalette make_palette() {
    return GenerativePalette(GradientShape::VIGNETTE, HarmonyType::ANALOGOUS,
                             BrightnessProfile::ASCENDING);
  }

  /**
   * @brief Constructs the effect, seeding the initial palette, palette normal,
   *        and filter pipeline.
   */
  FLASHMEM Dynamo()
      : Effect(W, H),
        palettes({make_palette()}),
        palette_normal(Z_AXIS),
        filters(Filter::World::Trails<W, TRAIL_CAPACITY>(
                    (uint32_t)params.trail_length),
                Filter::World::Replicate<W>(3),
                Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>()) {}

  /**
   * @brief POV column-strobe flag — see Effect::strobe_columns.
   * @return false; each lit column persists in the POV sweep until the next
   *         column overwrites it, with no black strobe between columns.
   */
  bool strobe_columns() const override { return false; }

  /**
   * @brief Forces full-canvas rendering per simulator worker.
   * @return decltype(filters)::any_crosses_segments — true, because the
   *         World::Trails stage reprojects trail samples under rotation and so
   *         moves state across segment bands; a band-clipped worker would drop
   *         cross-band trails. Trait-derived, so adding/removing a filter keeps
   *         the gate correct. See docs/segmented_stateful_effects_spec.md.
   */
  bool needs_full_frame() const override {
    return decltype(filters)::any_crosses_segments;
  }

  /**
   * @brief Registers sliders, seeds node rows, primes the baked-palette LUT
   *        pool, and schedules the random reverse/wipe/rotate timers.
   */
  void init() override {
    filters.template get<Filter::World::Trails<W, TRAIL_CAPACITY>>()
        .init_storage(persistent_arena);

    registerParam("Speed", &params.speed, -10.0f, 10.0f);
    registerParam("Gap", &params.gap, 1.0f, 20.0f);
    registerParam("Trail Len", &params.trail_length, 1.0f, 100.0f);
    registerParam("Wipe Dur", &params.wipe_duration, 1.0f, 100.0f);

    for (size_t i = 0; i < NUM_NODES; ++i) {
      nodes[i].y = i;
    }

    // Allocate the baked-palette LUT pool once (one 256-entry Color4 LUT per
    // possible live palette). bake() allocates; subsequent rebake() refills in
    // place with no allocation, so push/pop churn never grows the arena. Seed
    // every slot from the initial palette, then fill the active set.
    for (auto &bp : baked_palettes_)
      bp.bake(persistent_arena, palettes[0]);
    rebake_active_palettes();

    timeline
        .add(0, Animation::RandomTimer(
                    4, 64, [this](Canvas &) { reverse(); }, true))
        .add(0, Animation::RandomTimer(
                    20, 64, [this](Canvas &) { color_wipe(); }, true))
        .add(0, Animation::RandomTimer(
                    48, 160, [this](Canvas &) { rotate(); }, true));
  }

  /**
   * @brief Flips travel direction via a private sign.
   * @details Toggles speed_direction_ rather than mutating the registered
   *          "Speed" slider, so animation never overwrites the user's value.
   *          Effective speed is params.speed * speed_direction_ (applied in
   *          draw_frame).
   */
  void reverse() { speed_direction_ *= -1; }

  /**
   * @brief Schedules a half-turn rotation about a random axis, eased in/out.
   */
  void rotate() {
    timeline.add(0, Animation::Rotation<W>(orientation, random_vector(), PI_F,
                                           40, ease_in_out_sin, false));
  }

  /**
   * @brief Pushes a fresh palette at the front and animates its boundary angle
   *        from 0 up to PI, sweeping the new colors across the sphere.
   * @details Drops the wipe if the palette buffer is full.
   */
  void color_wipe() {
    if (palettes.is_full()) {
      hs::log("Palettes full, dropping color wipe!");
      return;
    }

    palettes.push_front(make_palette());
    palette_boundaries.push_front(0);
    // Re-bake the active set into the LUT pool: the per-pixel color() path reads
    // baked_palettes_, never the GenerativePalettes directly. A push_front
    // shifts every logical index, so rebake the whole active range.
    rebake_active_palettes();

    // On completion, stamp THIS wipe's own boundary slot with the WIPE_COMPLETE
    // sentinel rather than pop_back here: overlapping wipes can have different
    // durations and finish out of order, so a pop_back would evict the oldest
    // still-animating boundary. reap_completed_wipes() then collapses finished
    // wipes from the back in FIFO order without touching a live one. The slot
    // address stays valid because the palettes.is_full() guard above caps
    // palette_boundaries at its capacity, so the ring never wraps to overwrite a
    // live front slot the captured pointer aliases — loosen that guard and this
    // capture becomes a dangling write.
    float *boundary_slot = &palette_boundaries.front();
    timeline.add(0, Animation::Transition(palette_boundaries.front(), PI_F,
                                          (int)params.wipe_duration, ease_linear)
                        .then([boundary_slot]() { *boundary_slot = WIPE_COMPLETE; }));
  }

  /**
   * @brief Collapses color wipes whose Transition has finished (boundary
   *        stamped with WIPE_COMPLETE).
   * @details Pops only from the back (oldest first), so a wipe that finished
   *          early while an older one is still animating waits its FIFO turn
   *          and no live boundary is ever evicted. pop_back removes the highest
   *          logical index, leaving the front-indexed palettes and their baked
   *          LUTs in place, so no rebake is needed.
   */
  void reap_completed_wipes() {
    while (!palette_boundaries.is_empty() &&
           palette_boundaries.back() >= WIPE_COMPLETE) {
      palette_boundaries.pop_back();
      palettes.pop_back();
    }
  }

  /**
   * @brief Refills the baked LUT pool from the live GenerativePalettes so the
   *        hot per-pixel color() path is a table lookup instead of a per-call
   *        OKLCH lerp.
   * @details Called only when the palette set changes (a wipe push/pop) — never
   *          per frame and never per pixel.
   */
  void rebake_active_palettes() {
    for (size_t i = 0; i < palettes.size(); ++i)
      baked_palettes_[i].rebake(palettes[i]);
  }

  /**
   * @brief Picks the color for a direction at a palette parameter.
   * @param v Sample direction; the angle between it and palette_normal selects
   *          a palette band.
   * @param t Palette parameter in [0, 1] indexing into the baked LUT.
   * @return The blended Color4 for the sampled band, with a blend_width-wide
   *         crossfade across each boundary.
   * @details Walks the active palette boundaries from the newest (front) to the
   *          oldest to find the band containing the angle.
   */
  Color4 color(const Vector &v, float t) {
    // Half-width of the cross-fade region on each side of a palette boundary, in
    // radians. PI_F/4 (45 deg) gives a broad, soft seam between adjacent palette
    // bands relative to the [0, PI] angular span — wide enough that boundaries
    // read as gradients, not hard edges.
    constexpr float blend_width = PI_F / 4;
    // +inf sentinel for "no next boundary": `a` is an angle in [0, PI], so any
    // value well above PI makes the `a < next_boundary_lower_edge` test pass.
    constexpr float NO_NEXT_BOUNDARY = 100.0f;
    float a = angle_between(v, palette_normal);

    // palettes[i+1] stays in range via the invariant palettes.size() ==
    // palette_boundaries.size()+1 (color_wipe and the reap push/pop one of each
    // together). As a second guard, palette_boundaries' capacity is
    // MAX_PALETTES-1, so i+1 < palettes capacity even if that invariant broke;
    // it must stay below MAX_PALETTES.
    //
    // This scan assumes palette_boundaries is monotonically non-decreasing in
    // index (newest wipe at front=0 with the smallest angle, animating up to PI);
    // both the short-circuit and the gap test rely on it. A live Wipe-Dur change
    // can transiently invert the order when a newer, shorter wipe overtakes an
    // older one — purely cosmetic, self-heals as wipes drain, and stays in bounds
    // regardless.
    for (size_t i = 0; i < palette_boundaries.size(); ++i) {
      float boundary = palette_boundaries[i];
      auto lower_edge = boundary - blend_width;
      auto upper_edge = boundary + blend_width;

      if (a < lower_edge) {
        return baked_palettes_[i].get(t);
      }

      if (a >= lower_edge && a <= upper_edge) {
        auto blend_factor = (a - lower_edge) / (2 * blend_width);
        auto clamped_blend_factor = hs::clamp(blend_factor, 0.0f, 1.0f);

        Color4 c1 = baked_palettes_[i].get(t);
        Color4 c2 = baked_palettes_[i + 1].get(t);

        uint16_t fract = float_to_pixel16(clamped_blend_factor);
        return Color4(c1.color.lerp16(c2.color, fract),
                      hs::lerp(c1.alpha, c2.alpha, clamped_blend_factor));
      }

      auto next_boundary_lower_edge =
          (i + 1 < palette_boundaries.size()
               ? palette_boundaries[i + 1] - blend_width
               : NO_NEXT_BOUNDARY);

      if (a > upper_edge && a < next_boundary_lower_edge) {
        return baked_palettes_[i + 1].get(t);
      }
    }

    return baked_palettes_[0].get(t);
  }

  /**
   * @brief Renders one frame.
   * @details Syncs the live trail length, steps the timeline, advances the
   *          strand by the accumulated whole steps, plots it, then flushes the
   *          filter pipeline through color().
   */
  void draw_frame() override {
    Canvas canvas(*this);

    // Live "Trail Len" slider: push the (clamped) param into the Trails filter
    // each frame so dragging it retunes the trail length. Clamp into the
    // filter's [1,255] domain (the slider tops out at 100) so set_lifetime's
    // trap only ever fires on a genuine authoring error.
    filters.template get<Filter::World::Trails<W, TRAIL_CAPACITY>>()
        .set_lifetime(hs::clamp((int)params.trail_length, 1, 255));

    timeline.step(canvas);

    // Collapse any color wipes that finished this step, in FIFO order, before
    // their boundaries are read below.
    reap_completed_wipes();

    // Carry the fractional part of |speed| across frames so slow speeds still
    // advance the strand instead of dead-zoning to zero whenever |speed| < 1.
    // The integer part is consumed as whole steps this frame; the remainder
    // rolls forward.
    const float effective_speed = params.speed * speed_direction_;
    speed_accumulator_ += std::abs(effective_speed);
    const int steps = static_cast<int>(speed_accumulator_);
    speed_accumulator_ -= static_cast<float>(steps);
    // The i/steps division is provably guarded: the loop body only runs when
    // steps >= 1, so a sub-1 accumulator produces zero iterations (no
    // divide-by-zero).
    if (steps == 0) {
      // Zero-step frame (|speed| < 1, accumulator still sub-unit): re-emit the
      // strand at its current position without advancing it. Otherwise no plots
      // are produced this frame, and a Trail Len of 1 (ttl popped before
      // drawing) blanks the strand, so sub-unit speeds flicker.
      draw_nodes(canvas, 0.0f);
    } else {
      for (int i = steps - 1; i >= 0; --i) {
        pull(0, effective_speed);
        draw_nodes(canvas, static_cast<float>(i) / steps);
      }
    }

    // Trail re-draw: the Trails filter replays every buffered point and calls
    // this trailFn with t = the point's age fraction (1 - ttl/lifetime). We feed
    // that age straight in as color()'s palette parameter, so a trail fades
    // ALONG the palette with age (newest at t=0, oldest at t=1) rather than just
    // dimming — the trail's color sweep is its age axis.
    filters.flush(
        canvas, [this](const Vector &v, float t) { return color(v, t); }, 1.0f);
  }

private:
  // White-box test reaches palette_boundaries to stage the documented
  // overlapping-wipe band inversion and assert color() stays in-range.
  friend struct ::hs_test::effects_tests::DynamoWhiteBox;

  /**
   * @brief Plots the strand for this sub-step.
   * @param canvas Target canvas to plot into.
   * @param age Trail age fed to the Trails filter (0 = newest).
   * @details The head node is a single half-alpha point; each following node is
   *          a half-alpha line back to its predecessor.
   */
  void draw_nodes(Canvas &canvas, float age) {
    for (size_t i = 0; i < nodes.size(); ++i) {
      if (i == 0) {
        auto from = pixel_to_vector<W, H>(nodes[i].x, nodes[i].y);
        Color4 c = color(from, 0);
        c.alpha *= 0.5f;
        filters.plot(canvas, from, c.color, age, c.alpha);
      } else {
        auto from = pixel_to_vector<W, H>(nodes[i - 1].x, nodes[i - 1].y);
        auto to = pixel_to_vector<W, H>(nodes[i].x, nodes[i].y);
        auto fragment_shader = [this](const Vector &v, Fragment &f) {
          f.color = color(v, 0);
          f.color.alpha *= 0.5f;
        };
        Fragment f_from; f_from.pos = from;
        Fragment f_to;   f_to.pos = to;
        Plot::Line::draw<W, H>(filters, canvas, f_from, f_to,
                               fragment_shader);
      }
    }
  }

  /**
   * @brief Advances the strand one whole step.
   * @param leader Index of the leader node to move first.
   * @param effective_speed Signed speed whose direction the leader moves in.
   * @details Moves the leader node in the signed direction of effective_speed,
   *          then drags every other node toward it from both sides so the chain
   *          follows, keeping each link within `gap`.
   */
  void pull(int leader, float effective_speed) {
    nodes[leader].v = dir(effective_speed);
    move(nodes[leader]);
    for (int i = leader - 1; i >= 0; --i) {
      drag(nodes[i + 1], nodes[i]);
    }
    for (size_t i = leader + 1; i < nodes.size(); ++i) {
      drag(nodes[i - 1], nodes[i]);
    }
  }

  /**
   * @brief Pulls `follower` toward `leader`.
   * @param leader Node the follower chases.
   * @param follower Node moved this step.
   * @details If moving one step would leave the gap too wide, the follower
   *          adopts the leader's velocity and closes the slack until within
   *          `gap`; otherwise it just steps once.
   */
  void drag(Node &leader, Node &follower) {
    int dest = wrap(follower.x + follower.v, W);
    if (shortest_distance(dest, leader.x, W) > (int)params.gap) {
      follower.v = leader.v;
      while (shortest_distance(follower.x, leader.x, W) > (int)params.gap) {
        move(follower);
      }
    } else {
      move(follower);
    }
  }

  /**
   * @brief Advances a node by its velocity, wrapping x into [0, W).
   * @param node Node to move in place.
   */
  void move(Node &node) { node.x = wrap(node.x + node.v, W); }

  /**
   * @brief Computes the unit travel direction for a signed speed.
   * @param speed Signed speed value.
   * @return -1 for negative speed, otherwise +1.
   * @note `speed == 0` maps to +1, not a "stand still" direction, but this is
   *       unobservable: the sole caller pull() only runs dir() on frames with a
   *       whole step to take (|effective_speed| >= 1). The +1 tie-break is a
   *       harmless default; a future caller passing exactly 0 should handle its
   *       own stationary case rather than rely on it.
   */
  int dir(float speed) const { return speed < 0 ? -1 : 1; }

  Timeline timeline; /**< Drives reverse/wipe/rotate animations and timers. */

  static constexpr size_t MAX_PALETTES = 16; /**< Max live palettes. */
  /**
   * @brief Sentinel a completed wipe writes into its boundary slot so
   *        reap_completed_wipes() can collapse it.
   * @details Set well above the live boundary range [0, PI] so it can never
   *          collide with an in-flight value.
   */
  static constexpr float WIPE_COMPLETE = 100.0f;
  static constexpr int H_VIRT = H + hs::H_OFFSET; /**< Virtual row count. */
  static constexpr size_t NUM_NODES = H_VIRT; /**< Strand node count. */
  /**
   * @brief Compile-time Trails storage capacity (max buffered trail points).
   * @details DELIBERATE MEMORY-BUDGET CAP, not sized to the worst case. The
   *          active trail length is the runtime "Trail Len" slider; emission is
   *          slider-driven (up to floor(|Speed|) full-strand plots per frame,
   *          each point stored for up to 100 frames), so the true worst case
   *          reaches hundreds of thousands of points and has no honest small
   *          compile-time bound — unlike ChaoticStrings (TRAIL_LENGTH *
   *          ORIENTATION_SUBSTEPS). Over-budget is graceful and memory-safe:
   *          World::Trails' ring (filter.h push_back) drops the OLDEST point
   *          when full, only shortening the visible tail — it never traps,
   *          overruns the arena, or corrupts a frame. 10000 holds the trail at
   *          typical settings; a static_assert is intentionally omitted as no
   *          truthful worst-case bound exists.
   */
  static constexpr int TRAIL_CAPACITY = 10000;
  StaticCircularBuffer<GenerativePalette, MAX_PALETTES> palettes; /**< Live palettes. */
  StaticCircularBuffer<float, MAX_PALETTES - 1> palette_boundaries; /**< Wipe boundary angles. */
  /**
   * @brief Baked 256-entry LUTs mirroring palettes[] in logical order.
   * @details The per-pixel color() path reads these (fast lerp16 table lookup)
   *          instead of evaluating GenerativePalette's per-call OKLCH lerp
   *          thousands of times per frame. Kept in sync by
   *          rebake_active_palettes() on every wipe push/pop; one slot per
   *          possible live palette so churn never reallocates.
   */
  std::array<BakedPalette, MAX_PALETTES> baked_palettes_;

  Vector palette_normal; /**< Reference axis for band angle selection. */
  std::array<Node, NUM_NODES> nodes; /**< The strand nodes. */

  /**
   * @brief Travel direction toggled by reverse().
   * @details Kept separate from the registered "Speed" slider so animation
   *          never clobbers the user's value.
   */
  int speed_direction_ = 1;
  /**
   * @brief Fractional-step carry so |speed| < 1 still advances the strand over
   *        multiple frames instead of truncating to zero.
   */
  float speed_accumulator_ = 0.0f;

  /**
   * @brief Live slider-backed parameters for the effect.
   */
  struct Params {
    float speed = 2.0f; /**< Strand travel speed. */
    float gap = 5.0f; /**< Target spacing between adjacent nodes. */
    float trail_length = 8.0f; /**< Active trail length. */
    float wipe_duration = 20.0f; /**< Color-wipe transition duration. */
  } params;

  Orientation<> orientation; /**< Current sphere orientation. */

  /**
   * @brief Filter pipeline applied to plotted points before color resolution.
   */
  Pipeline<W, H, Filter::World::Trails<W, TRAIL_CAPACITY>,
           Filter::World::Replicate<W>,
           Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Dynamo)
