/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file Dynamo.h
 * @brief A snaking strand of trailing nodes under palettes that sweep in via
 *        angular color wipes.
 */

#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"

// Unit-test accessor reaching palette_boundaries to stage the overlapping-wipe
// band inversion and assert color() stays memory-safe and in-range.
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
   * @brief The effect's canonical generative-palette recipe, shared by the
   *        initial palette and every color-wipe palette.
   */
  static GenerativePalette make_palette() {
    return GenerativePalette{EffectPaletteRecipes::dynamo(
        EffectPaletteRecipes::random_base_turns())};
  }

  /**
   * @brief Constructs the effect, seeding the initial palette, palette normal,
   *        and filter pipeline.
   */
  HS_COLD_MEMBER Dynamo()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})),
        palettes{make_palette()}, palette_normal(Z_AXIS),
        // draw_frame() pushes the live "Trail Len" slider before the first
        // flush(), so this seed lifetime never reaches the output.
        filters(Filter::World::Trails<TRAIL_CAPACITY>(1),
                Filter::World::Replicate<W>(STRAND_COPIES),
                Filter::World::Orient(orientation),
                Filter::Screen::AntiAlias<W, H>()) {}

  /**
   * @brief Registers sliders, seeds node rows, primes the baked-palette LUT
   *        pool, and schedules the random reverse/wipe/rotate timers.
   */
  void init() override {
    filters.template get<Filter::World::Trails<TRAIL_CAPACITY>>().init_storage(
        persistent_arena);

    nodes = persistent_arena.make_n<Node>(NUM_NODES);

    register_param("Speed", &params.speed, -10.0f, 10.0f);
    register_param("Gap", &params.gap, 1.0f, GAP_MAX);
    register_param("Trail Len", &params.trail_length, 1.0f,
                   static_cast<float>(TRAIL_LEN_MAX));
    register_readonly_param("Trail Cap", &params.trail_ceiling, 1.0f,
                            static_cast<float>(TRAIL_LEN_MAX));
    register_param("Wipe Dur", &params.wipe_duration, 1.0f, 100.0f);

    for (size_t i = 0; i < NUM_NODES; ++i) {
      nodes[i].y = i;
    }

    // Allocate the LUT pool once: bake() allocates, rebake() refills in place,
    // so push/pop churn never grows the arena. The bake also seeds slot 0, the
    // only live palette at this point.
    for (auto &bp : baked_palettes)
      bp.bake(persistent_arena, palettes[0]);

    timeline
        .add(0, Animation::RandomTimer(
                    4, 64, [this](Canvas &) { reverse(); }, true))
        .add(0, Animation::RandomTimer(
                    20, 64, [this](Canvas &) { color_wipe(); }, true))
        .add(0, Animation::RandomTimer(
                    48, 160, [this](Canvas &) { rotate(); }, true));
  }

  /**
   * @brief Renders one frame.
   * @details Syncs the live trail length, steps the timeline, advances the
   *          strand by the accumulated whole steps, plots it, then flushes the
   *          filter pipeline through color().
   */
  void draw_frame() override {
    Canvas canvas(*this);

    // Push the live "Trail Len" slider into the Trails filter, capped to what
    // the ring can hold for the current emission rate. The cap is published as
    // read-only telemetry so the GUI can show why a long trail renders short.
    const int ceiling = trail_length_ceiling();
    params.trail_ceiling = static_cast<float>(ceiling);
    filters.template get<Filter::World::Trails<TRAIL_CAPACITY>>().set_lifetime(
        hs::clamp((int)params.trail_length, 1, ceiling));

    {
      HS_PROFILE(dy_timeline_step);
      timeline.step(canvas);
    }

    // Collapse finished wipes (FIFO) before their boundaries are read below.
    reap_completed_wipes();

    // color() walks palette_boundaries and reads baked_palettes[i] and [i+1];
    // paired push (color_wipe) and pop (reap) keep palettes one ahead of the
    // boundaries. Guard that pairing once here before the per-pixel band walk.
    HS_CHECK(palettes.size() == palette_boundaries.size() + 1);

    // Carry the fractional part of |speed| across frames so |speed| < 1 still
    // advances the strand instead of truncating to zero.
    const float effective_speed = params.speed * speed_direction;
    speed_accumulator += std::abs(effective_speed);
    const int steps = static_cast<int>(speed_accumulator);
    speed_accumulator -= static_cast<float>(steps);
    emitted_points = 0;
    {
      HS_PROFILE(dy_draw_nodes);
      if (steps == 0) {
        // Re-emit the strand in place; otherwise a Trail Len of 1 blanks it on
        // zero-step frames and sub-unit speeds flicker.
        draw_nodes(canvas, 0.0f);
      } else {
        for (int i = steps - 1; i >= 0; --i) {
          pull(effective_speed);
          draw_nodes(canvas, static_cast<float>(i) / steps);
        }
      }
    }
    emission_points = emitted_points / (steps > 0 ? steps : 1);

    // The Trails filter replays each buffered point with t = its age fraction;
    // feeding that as color()'s palette parameter fades the trail along the
    // palette with age (newest t=0, oldest t=1) rather than just dimming.
    {
      HS_PROFILE(dy_filter_flush);
      filters.flush(
          canvas, [this](const Vector &v, float t) { return color(v, t); },
          1.0f);
    }
  }

private:
  friend struct ::hs_test::effects_tests::DynamoWhiteBox;

  /**
   * @brief Flips travel direction via a private sign.
   * @details Toggles speed_direction so animation never overwrites the "Speed"
   *          slider; effective speed is params.speed * speed_direction.
   */
  void reverse() { speed_direction *= -1; }

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
   * @details Drops the wipe if the boundary buffer is full, logging once until
   *          a wipe lands again. That buffer is the
   *          binding capacity (MAX_PALETTES - 1) and is what keeps boundary_slot
   *          from aliasing a reissued ring slot; palettes stays one ahead of it,
   *          so its own push is covered by the same guard.
   */
  void color_wipe() {
    if (palette_boundaries.is_full()) {
      if (!logged_palettes_full) {
        logged_palettes_full = true;
        hs::log("Dynamo: palettes full, dropping color wipe");
      }
      return;
    }
    logged_palettes_full = false;

    palettes.push_front(make_palette());
    palette_boundaries.push_front(0);
    rotate_and_rebake_front();

    // Stamp WIPE_COMPLETE on completion (not pop_back): overlapping wipes can
    // finish out of order, so pop_back would evict a still-animating boundary.
    float *boundary_slot = &palette_boundaries.front();
    timeline.add(
        0, Animation::Transition(palette_boundaries.front(), PI_F,
                                 (int)params.wipe_duration, ease_linear)
               .then([boundary_slot]() { *boundary_slot = WIPE_COMPLETE; }));
  }

  /**
   * @brief Collapses color wipes whose Transition has finished (boundary
   *        stamped with WIPE_COMPLETE).
   * @details Pops only from the back (FIFO, oldest first), so a wipe that
   *          finished early waits its turn and no live boundary is evicted;
   *          front-indexed palettes and LUTs stay in place, needing no rebake.
   */
  void reap_completed_wipes() {
    while (!palette_boundaries.is_empty() &&
           palette_boundaries.back() >= WIPE_COMPLETE) {
      palette_boundaries.pop_back();
      palettes.pop_back();
    }
  }

  /**
   * @brief Realigns the LUT pool after a palettes push_front and bakes the new
   *        front.
   * @details baked_palettes mirrors palettes[] by logical index, so a push_front
   *          shifts every entry by one. Rotating the pointer-sized LUT handles
   *          carries each existing bake to its new index, leaving one
   *          GenerativePalette evaluation pass per wipe instead of one per live
   *          palette. Slot MAX_PALETTES-1 is always dead here (color_wipe()
   *          drops the wipe when the boundary buffer is full, so palettes holds
   *          at most MAX_PALETTES-1 entries before its push_front), so recycling
   *          it to the front strands no LUT storage.
   */
  void rotate_and_rebake_front() {
    BakedPalette recycled = baked_palettes[MAX_PALETTES - 1];
    for (size_t i = MAX_PALETTES - 1; i > 0; --i)
      baked_palettes[i] = baked_palettes[i - 1];
    baked_palettes[0] = recycled;
    baked_palettes[0].rebake(palettes[0]);
  }

  /**
   * @brief Picks the color for a direction at a palette parameter.
   * @param v Unit sample direction; the angle between it and palette_normal
   *          selects a palette band.
   * @param t Palette parameter in [0, 1] indexing into the baked LUT.
   * @return The blended Color4 for the sampled band, with a blend_width-wide
   *         crossfade across each boundary.
   * @details Walks the active palette boundaries from the newest (front) to the
   *          oldest to find the band containing the angle.
   */
  Color4 color(const Vector &v, float t) {
    if (palette_boundaries.size() == 0)
      return baked_palettes[0].get(t);

    // Cross-fade half-width per boundary side, in radians.
    constexpr float blend_width = PI_F / 4;
    // Sentinel for "no next boundary": `a` is in [0, PI], so any value above PI
    // makes the `a < next_boundary_lower_edge` test pass.
    constexpr float NO_NEXT_BOUNDARY = 100.0f;
    float a = fast_acos(hs::clamp(dot(v, palette_normal), -1.0f, 1.0f));

    // The scan assumes palette_boundaries is monotonically non-decreasing; a live
    // Wipe-Dur change can transiently invert it, picking a stale palette for a few
    // frames. Stays in bounds and self-heals as wipes drain.
    for (size_t i = 0; i < palette_boundaries.size(); ++i) {
      float boundary = palette_boundaries[i];
      auto lower_edge = boundary - blend_width;
      auto upper_edge = boundary + blend_width;

      if (a < lower_edge) {
        return baked_palettes[i].get(t);
      }

      if (a >= lower_edge && a <= upper_edge) {
        auto blend_factor = (a - lower_edge) / (2 * blend_width);
        auto clamped_blend_factor = hs::clamp(blend_factor, 0.0f, 1.0f);

        Color4 c1 = baked_palettes[i].get(t);
        Color4 c2 = baked_palettes[i + 1].get(t);

        uint16_t fract = float_to_pixel16(clamped_blend_factor);
        return Color4(c1.color.lerp16(c2.color, fract),
                      hs::lerp(c1.alpha, c2.alpha, clamped_blend_factor));
      }

      auto next_boundary_lower_edge =
          (i + 1 < palette_boundaries.size()
               ? palette_boundaries[i + 1] - blend_width
               : NO_NEXT_BOUNDARY);

      if (a > upper_edge && a < next_boundary_lower_edge) {
        return baked_palettes[i + 1].get(t);
      }
    }

    return baked_palettes[0].get(t);
  }

  /**
   * @brief Plots the strand for this sub-step.
   * @param canvas Target canvas to plot into.
   * @param age Trail age fed to the Trails filter (0 = newest).
   * @details The head node is a single half-alpha point; each following node is
   *          a half-alpha line back to its predecessor. Trails buffers these
   *          points and the same-frame flush() re-emits them at t = 0, so the
   *          live strand composites twice: once here at half alpha, once at
   *          color()'s full alpha.
   */
  void draw_nodes(Canvas &canvas, float age) {
    for (size_t i = 0; i < NUM_NODES; ++i) {
      if (i == 0) {
        auto from = pixel_to_vector<W, H>(nodes[i].x, nodes[i].y);
        Color4 c = color(from, 0);
        c.alpha *= 0.5f;
        ++emitted_points;
        filters.plot(canvas, from, c.color, age, c.alpha);
      } else {
        auto from = pixel_to_vector<W, H>(nodes[i - 1].x, nodes[i - 1].y);
        auto to = pixel_to_vector<W, H>(nodes[i].x, nodes[i].y);
        auto fragment_shader = [this](const Vector &v, Fragment &f) {
          f.color = color(v, 0);
          f.color.alpha *= 0.5f;
          ++emitted_points;
        };
        Fragment f_from;
        f_from.pos = from;
        f_from.age = age;
        Fragment f_to;
        f_to.pos = to;
        f_to.age = age;
        Plot::Line::draw<W, H>(filters, canvas, f_from, f_to, fragment_shader);
      }
    }
  }

  /**
   * @brief Longest trail, in frames, the ring can hold at the current rate.
   * @return Frame count in [1, TRAIL_LEN_MAX]; TRAIL_LEN_MAX before the first
   *         frame is measured.
   * @details Steady-state occupancy is points-per-frame x lifetime, so a longer
   *          trail than this would overrun the ring and evict live points of
   *          arbitrary age (flush()'s compaction leaves the ring unordered),
   *          punching holes in the tail rather than shortening it. The step
   *          count is bounded by |speed| + 1: speed_accumulator carries at most
   *          one extra whole step into any frame.
   */
  int trail_length_ceiling() const {
    if (emission_points == 0)
      return TRAIL_LEN_MAX;
    const uint32_t emissions =
        static_cast<uint32_t>(std::abs(params.speed)) + 1;
    return hs::clamp(
        static_cast<int>(TRAIL_CAPACITY / (emission_points * emissions)), 1,
        TRAIL_LEN_MAX);
  }

  /**
   * @brief Advances the strand one whole step.
   * @param effective_speed Signed speed whose direction the head node moves in.
   * @details Moves node 0 in the signed direction of effective_speed, then drags
   *          each following node toward its predecessor so the chain follows,
   *          keeping every link within `gap`.
   */
  void pull(float effective_speed) {
    nodes[0].v = dir(effective_speed);
    move(nodes[0]);
    for (size_t i = 1; i < NUM_NODES; ++i) {
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
   * @note The slack loop terminates only while `leader.v != 0` — move() is a
   *       no-op at zero velocity. It holds by construction: a node still at
   *       v == 0 has never moved, so it shares its leader's column and the gap
   *       test is false.
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
   * @note `speed == 0` maps to +1, but this is unobservable: a zero speed never
   *       advances the step accumulator, so pull() never runs.
   */
  int dir(float speed) const { return speed < 0 ? -1 : 1; }

  /**
   * @brief Current sphere orientation.
   * @details Declared before `timeline` so it outlives the Rotations that point
   * here, which ~Timeline clears on teardown.
   */
  Orientation<> orientation;
  Timeline timeline; /**< Drives reverse/wipe/rotate animations and timers. */

  static constexpr size_t MAX_PALETTES = 16; /**< Max live palettes. */
  static constexpr int TRAIL_LEN_MAX =
      100; /**< "Trail Len" slider max, and the ceiling "Trail Cap" reports. */
  /**
   * @brief Upper bound of the "Gap" slider (target node spacing).
   * @details Held below W/2 so drag()'s slack-closing loop always terminates:
   *          shortest_distance() saturates at W/2, so a gap that can actually be
   *          reached keeps the loop making progress instead of circling forever.
   */
  static constexpr float GAP_MAX = 20.0f;
  static_assert(2.0f * GAP_MAX < static_cast<float>(W),
                "Gap max must stay below W/2 so drag() terminates");
  /**
   * @brief Sentinel a completed wipe writes into its boundary slot so
   *        reap_completed_wipes() can collapse it.
   * @details Set well above the live boundary range [0, PI] so it can never
   *          collide with an in-flight value.
   */
  static constexpr float WIPE_COMPLETE = 100.0f;
  static constexpr int H_VIRT = H + hs::H_OFFSET; /**< Virtual row count. */
  static constexpr size_t NUM_NODES = H_VIRT;     /**< Strand node count. */
  /** @brief Evenly spaced Y-axis copies of the strand the pipeline emits. */
  static constexpr int STRAND_COPIES = 3;
  /**
   * @brief Compile-time Trails storage capacity (max buffered trail points).
   * @details Sized to the persistent partition left by the nodes and the baked
   *          palette LUTs; trail_length_ceiling() keeps the live trail inside it
   *          so the ring never evicts.
   */
  static constexpr int TRAIL_CAPACITY = 29000;
  StaticCircularBuffer<GenerativePalette, MAX_PALETTES>
      palettes; /**< Live palettes. */
  StaticCircularBuffer<float, MAX_PALETTES - 1>
      palette_boundaries; /**< Wipe boundary angles. */
  /**
   * @brief Baked 256-entry LUTs mirroring palettes[] in logical order.
   * @details Read by the per-pixel color() path (lerp16 lookup, not OKLCH lerp);
   *          a wipe push rotates them to match (rotate_and_rebake_front) and a
   *          reap pops from the back, which shifts nothing. One slot per possible
   *          live palette, so churn never reallocates.
   */
  std::array<BakedPalette, MAX_PALETTES> baked_palettes;

  // init() allocates the nodes, Trails ring buffer, and baked palette LUTs from
  // the persistent arena.
  static constexpr size_t FOOTPRINT_BYTES =
      NUM_NODES * sizeof(Node) +
      TRAIL_CAPACITY *
          sizeof(typename Filter::World::Trails<TRAIL_CAPACITY>::Item) +
      MAX_PALETTES * BakedPalette::required_arena_bytes();
  // Effect keeps the default arena split, so the footprint must fit the device
  // persistent partition. Guards a TRAIL_CAPACITY/MAX_PALETTES retune.
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "Dynamo persistent footprint exceeds the default partition; "
                "retune TRAIL_CAPACITY/MAX_PALETTES or carve arenas");

  Vector palette_normal; /**< Reference axis for band angle selection. */
  Node *nodes = nullptr; /**< Arena-backed strand nodes. */

  /**
   * @brief Travel direction toggled by reverse(); kept separate from the "Speed"
   *        slider so animation never clobbers the user's value.
   */
  int speed_direction = 1;
  /**
   * @brief Fractional-step carry so |speed| < 1 still advances the strand over
   *        multiple frames instead of truncating to zero.
   */
  float speed_accumulator = 0.0f;

  uint32_t emitted_points = 0;  /**< Points plotted this frame. */
  uint32_t emission_points = 0; /**< Points per strand emission, last frame. */

  /** @brief Palettes-full log latch; cleared when a wipe lands. */
  bool logged_palettes_full = false;

  /**
   * @brief Live slider-backed parameters for the effect.
   * @details trail_ceiling is engine-written (read-only); the rest are sliders.
   */
  struct Params {
    float speed = 2.0f;          /**< Strand travel speed. */
    float gap = 5.0f;            /**< Target spacing between adjacent nodes. */
    float trail_length = 8.0f;   /**< Active trail length. */
    float wipe_duration = 20.0f; /**< Color-wipe transition duration. */
    /** Longest trail the ring can hold, in frames (engine-written). */
    float trail_ceiling = static_cast<float>(TRAIL_LEN_MAX);
  } params;

  /**
   * @brief Filter pipeline applied to plotted points before color resolution.
   */
  Pipeline<W, H, Filter::World::Trails<TRAIL_CAPACITY>,
           Filter::World::Replicate<W>, Filter::World::Orient,
           Filter::Screen::AntiAlias<W, H>>
      filters;
};

#include "core/control/registry.h"
REGISTER_EFFECT(Dynamo)
