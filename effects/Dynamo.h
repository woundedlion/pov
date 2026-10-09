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

#include "core/animation/orientation.h"
#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"

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
  static constexpr const char *EFFECT_ID = "Dynamo";

  /**
   * @brief Constructs the effect, seeding the initial palette and filter
   *        pipeline.
   */
  HS_COLD_MEMBER Dynamo()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})),
        palettes{make_palette()},
        filters(Filter::World::Trails<TRAIL_CAPACITY>(1),
                Filter::World::Replicate<W>(STRAND_COPIES),
                Filter::World::Orient(orientation),
                Filter::Screen::AntiAlias<W, H>()) {}

  /**
   * @brief Registers sliders, seeds node rows, primes the baked-palette LUT
   *        pool, and schedules the random reverse/wipe/rotate timers.
   */
  void init() override {
    filters.init_storage(persistent_arena);

    nodes = persistent_arena.make_n<Node>(NUM_NODES);

    register_param("Speed", &params.speed, -SPEED_MAX, SPEED_MAX);
    register_param("Gap", &params.gap, 1.0f, GAP_MAX);
    register_int_param("Trail Len", &params.trail_length, 1, TRAIL_LEN_MAX);
    register_readonly_param("Trail Cap", &params.trail_ceiling, 1.0f,
                            static_cast<float>(TRAIL_LEN_MAX));
    register_param("Wipe Dur", &params.wipe_duration, 1.0f, 100.0f);

    for (size_t i = 0; i < NUM_NODES; ++i) {
      nodes[i].y = math::phi_to_y<H>(static_cast<float>(i) * math::PI_F /
                                     (NUM_NODES - 1));
    }

    // Allocate the LUT pool once; rebake() refills in place.
    for (auto &bp : baked_palettes)
      bp.bake(persistent_arena, palettes[0]);

    timeline
        .add(0, Animation::RandomTimer({.min = 4, .max = 64, .repeat = true},
                                       [this](Canvas &) { reverse(); }))
        .add(0, Animation::RandomTimer({.min = 20, .max = 64, .repeat = true},
                                       [this](Canvas &) { color_wipe(); }))
        .add(0, Animation::RandomTimer({.min = 48, .max = 160, .repeat = true},
                                       [this](Canvas &) { rotate(); }));
  }

  /**
   * @brief Renders one frame.
   */
  void draw_frame() override {
    Canvas canvas(*this);

    // Cap "Trail Len" to what the ring holds at the current emission rate;
    // the cap is published as read-only telemetry.
    const int ceiling = trail_length_ceiling();
    params.trail_ceiling = static_cast<float>(ceiling);
    filters.template get<Filter::World::Trails<TRAIL_CAPACITY>>().set_lifetime(
        hs::clamp((int)params.trail_length, 1, ceiling));

    {
      HS_PROFILE(dy_timeline_step);
      timeline.step(canvas);
    }

    // Collapse finished wipes (FIFO).
    reap_completed_wipes();

    HS_CHECK(palettes.size() == palette_boundaries.size() + 1,
             "Dynamo: palettes must stay one ahead of palette_boundaries "
             "(color_wipe push / reap pop); color() reads [i] and [i+1]");

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
    points_per_emission = emitted_points / (steps > 0 ? steps : 1);

    // Trails replays each point with t = its age fraction (newest 0, oldest
    // 1), so the trail fades along the palette.
    {
      HS_PROFILE(dy_filter_flush);
      filters.flush(
          canvas,
          [this](const math::Vector &v, float t) { return color(v, t); }, 1.0f);
    }
  }

private:
  friend struct ::hs_test::effects_tests::DynamoWhiteBox;

  /**
   * @brief One point on the strand: grid position (x,y) and per-step velocity.
   */
  struct Node {
    /**
     * @brief Constructs a node at the origin with zero velocity.
     */
    Node() : x(0), y(0), v(0) {}

    int x;   /**< Grid column. */
    float y; /**< Fractional display row. */
    int v;   /**< Per-step velocity along x. */
  };

  /**
   * @brief The effect's canonical generative-palette recipe.
   */
  static GenerativePalette make_palette() {
    return GenerativePalette{EffectPaletteRecipes::dynamo(
        EffectPaletteRecipes::random_base_turns())};
  }

  /**
   * @brief Flips travel direction via a private sign.
   * @details Effective speed is params.speed * speed_direction.
   */
  void reverse() { speed_direction *= -1; }

  /**
   * @brief Schedules a half-turn rotation about a random axis, eased in/out.
   */
  void rotate() {
    timeline.add(0, Animation::Rotation<W>(orientation, math::random_vector(),
                                           math::PI_F, 40,
                                           math::ease_in_out_sin, false));
  }

  /**
   * @brief Pushes a fresh palette at the front and animates its boundary angle
   *        from -WIPE_BLEND_WIDTH up to PI + WIPE_BLEND_WIDTH, sweeping the
   *        new colors and their blend band across the whole sphere.
   * @details Drops the wipe, logging once, while the boundary buffer is full;
   *          that bound keeps boundary_slot from aliasing a reissued ring slot.
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
    palette_boundaries.push_front(-WIPE_BLEND_WIDTH);
    rotate_and_rebake_front();

    // Stamp WIPE_COMPLETE on completion (not pop_back): overlapping wipes can
    // finish out of order, so pop_back would evict a still-animating boundary.
    float *boundary_slot = &palette_boundaries.front();
    timeline.add(
        0, Animation::Transition(palette_boundaries.front(),
                                 math::PI_F + WIPE_BLEND_WIDTH,
                                 (int)params.wipe_duration, math::ease_linear)
               .then([boundary_slot]() { *boundary_slot = WIPE_COMPLETE; }));
  }

  /**
   * @brief Collapses color wipes whose Transition has finished (boundary
   *        stamped with WIPE_COMPLETE).
   * @details Pops only from the back (FIFO), so a wipe that finished early
   *          waits its turn.
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
   * @details Slot MAX_PALETTES-1 is free to recycle: palettes holds at most
   *          MAX_PALETTES-1 entries before its push_front.
   */
  void rotate_and_rebake_front() {
    BakedPaletteStorage recycled = std::move(baked_palettes[MAX_PALETTES - 1]);
    for (size_t i = MAX_PALETTES - 1; i > 0; --i)
      baked_palettes[i] = std::move(baked_palettes[i - 1]);
    baked_palettes[0] = std::move(recycled);
    baked_palettes[0].rebake(palettes[0]);
  }

  /**
   * @brief Picks the color for a direction at a palette parameter.
   * @param v Unit sample direction; the angle between it and PALETTE_NORMAL
   *          selects a palette band.
   * @param t Palette parameter in [0, 1] indexing into the baked LUT.
   * @return The blended Color4 for the sampled band, with a WIPE_BLEND_WIDTH-wide
   *         crossfade across each boundary.
   */
  Color4 color(const math::Vector &v, float t) {
    if (palette_boundaries.size() == 0)
      return baked_palettes[0].get(t);

    // Sentinel for "no next boundary": `a` is in [0, PI], so any value above PI
    // makes the `a < next_boundary_lower_edge` test pass.
    constexpr float NO_NEXT_BOUNDARY = 100.0f;
    float a =
        math::fast_acos(hs::clamp(math::dot(v, PALETTE_NORMAL), -1.0f, 1.0f));

    // Assumes non-decreasing palette_boundaries; a live Wipe-Dur change can
    // transiently invert them, picking a stale palette until wipes drain.
    for (size_t i = 0; i < palette_boundaries.size(); ++i) {
      float boundary = palette_boundaries[i];
      auto lower_edge = boundary - WIPE_BLEND_WIDTH;
      auto upper_edge = boundary + WIPE_BLEND_WIDTH;

      if (a < lower_edge) {
        return baked_palettes[i].get(t);
      }

      if (a >= lower_edge && a <= upper_edge) {
        auto blend_factor = (a - lower_edge) / (2 * WIPE_BLEND_WIDTH);
        auto clamped_blend_factor = hs::clamp(blend_factor, 0.0f, 1.0f);

        Color4 c1 = baked_palettes[i].get(t);
        Color4 c2 = baked_palettes[i + 1].get(t);

        uint16_t fract = float_to_pixel16(clamped_blend_factor);
        return Color4(c1.color.lerp16(c2.color, fract),
                      hs::lerp(c1.alpha, c2.alpha, clamped_blend_factor));
      }

      auto next_boundary_lower_edge =
          (i + 1 < palette_boundaries.size()
               ? palette_boundaries[i + 1] - WIPE_BLEND_WIDTH
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
   * @details The head node is a half-alpha point; each following node is a
   *          half-alpha line back to its predecessor.
   */
  void draw_nodes(Canvas &canvas, float age) {
    for (size_t i = 0; i < NUM_NODES; ++i) {
      if (i == 0) {
        auto from = math::pixel_to_vector<W, H>(static_cast<float>(nodes[i].x),
                                                nodes[i].y);
        Color4 c = color(from, 0);
        c.alpha *= 0.5f;
        ++emitted_points;
        filters.plot(canvas, from, c.color, age, c.alpha);
      } else {
        auto from = math::pixel_to_vector<W, H>(
            static_cast<float>(nodes[i - 1].x), nodes[i - 1].y);
        auto to = math::pixel_to_vector<W, H>(static_cast<float>(nodes[i].x),
                                              nodes[i].y);
        auto fragment_shader = [this](const math::Vector &v, Fragment &f) {
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
   * @details Steady-state occupancy is points-per-frame x lifetime; a longer
   *          trail overruns the unordered ring and evicts points of arbitrary
   *          age. A frame takes at most |speed| + 1 steps.
   */
  int trail_length_ceiling() const {
    if (points_per_emission == 0)
      return TRAIL_LEN_MAX;
    const uint32_t emissions =
        static_cast<uint32_t>(std::abs(params.speed)) + 1;
    return hs::clamp(
        static_cast<int>(TRAIL_CAPACITY / (points_per_emission * emissions)), 1,
        TRAIL_LEN_MAX);
  }

  /**
   * @brief Advances the strand one whole step.
   * @param effective_speed Signed speed whose direction the head node moves in.
   * @details Drags each following node toward its predecessor, keeping every
   *          link within `gap`.
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
   * @note The slack loop needs `leader.v != 0`; a node still at v == 0 shares
   *       its leader's column, so the gap test is false.
   */
  void drag(Node &leader, Node &follower) {
    int dest = math::wrap(follower.x + follower.v, W);
    if (math::shortest_distance(dest, leader.x, W) > (int)params.gap) {
      follower.v = leader.v;
      while (math::shortest_distance(follower.x, leader.x, W) >
             (int)params.gap) {
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
  void move(Node &node) { node.x = math::wrap(node.x + node.v, W); }

  /**
   * @brief Computes the unit travel direction for a signed speed.
   * @param speed Signed speed value.
   * @return -1 for negative speed, otherwise +1.
   */
  int dir(float speed) const { return speed < 0 ? -1 : 1; }

  /**
   * @brief Current sphere orientation.
   * @details Declared before `timeline` so it outlives the Rotations that point
   * here, which ~Timeline clears on teardown.
   */
  math::Orientation<> orientation;
  Timeline timeline; /**< Drives reverse/wipe/rotate animations and timers. */

  static constexpr size_t MAX_PALETTES = 16; /**< Max live palettes. */
  static_assert(MAX_PALETTES + 3 <= Timeline::MAX_EVENTS,
                "Dynamo needs three timers, one rotation, and all live wipes");
  static constexpr int TRAIL_LEN_MAX =
      100; /**< "Trail Len" slider max, and the ceiling "Trail Cap" reports. */
  static constexpr float SPEED_MAX = 10.0f;
  /** @brief Upper bound of the "Gap" slider; below W/2 so the gap constraint can bind. */
  static constexpr float GAP_MAX = 20.0f;
  static_assert(2.0f * GAP_MAX < static_cast<float>(W),
                "Gap max must stay below W/2 so the gap constraint can bind");
  static constexpr float WIPE_BLEND_WIDTH = math::PI_F / 4;
  /**
   * @brief Sentinel a completed wipe writes into its boundary slot so
   *        reap_completed_wipes() can collapse it.
   * @details The live range is [-WIPE_BLEND_WIDTH, PI + WIPE_BLEND_WIDTH].
   */
  static constexpr float WIPE_COMPLETE = 100.0f;
#if HS_RUNTIME_DISPLAY_GEOMETRY
  static constexpr size_t NODE_CAPACITY = 2 * (H - 1) + 1;
  const size_t NUM_NODES =
      static_cast<size_t>(math::PI_F * math::ROWS_PER_RADIAN<H> + 0.999f) + 1;
#else
  static constexpr size_t NUM_NODES =
      static_cast<size_t>(math::PI_F * math::ROWS_PER_RADIAN<H> + 0.999f) + 1;
  static constexpr size_t NODE_CAPACITY = NUM_NODES;
#endif
  /** @brief Evenly spaced Y-axis copies of the strand the pipeline emits. */
  static constexpr int STRAND_COPIES = 3;
  /** @brief Reference axis for band angle selection. */
  static constexpr math::Vector PALETTE_NORMAL = math::Z_AXIS;
  /**
   * @brief Compile-time Trails storage capacity (max buffered trail points).
   * @details trail_length_ceiling() keeps the live trail inside it.
   */
  static constexpr int TRAIL_CAPACITY = 29000;
  StaticCircularBuffer<GenerativePalette, MAX_PALETTES>
      palettes; /**< Live palettes. */
  StaticCircularBuffer<float, MAX_PALETTES - 1>
      palette_boundaries; /**< Wipe boundary angles. */
  /**
   * @brief Baked 256-entry LUTs mirroring palettes[] in logical order.
   * @details A wipe push rotates them to match (rotate_and_rebake_front); a
   *          reap pops from the back, which shifts nothing.
   */
  std::array<BakedPaletteStorage, MAX_PALETTES> baked_palettes;

  // Persistent allocations: nodes, the Trails ring and the baked palette LUTs.
  static constexpr size_t FOOTPRINT_BYTES =
      NODE_CAPACITY * sizeof(Node) +
      TRAIL_CAPACITY *
          sizeof(typename Filter::World::Trails<TRAIL_CAPACITY>::Item) +
      MAX_PALETTES * BakedPalette::required_arena_bytes();
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "Dynamo persistent footprint exceeds the default partition; "
                "retune TRAIL_CAPACITY/MAX_PALETTES or carve arenas");

  Node *nodes = nullptr; /**< Arena-backed strand nodes. */

  /**
   * @brief Travel direction toggled by reverse(); separate from the "Speed"
   *        slider.
   */
  int speed_direction = 1;
  /**
   * @brief Fractional-step carry so |speed| < 1 still advances the strand.
   */
  float speed_accumulator = 0.0f;

  uint32_t emitted_points = 0; /**< Points plotted this frame. */
  /** @brief Points one strand emission plotted, measured last frame. */
  uint32_t points_per_emission = 0;

  /** @brief Palettes-full log latch; cleared when a wipe lands. */
  bool logged_palettes_full = false;

  /**
   * @brief Filter pipeline applied to plotted points before color resolution.
   */
  Pipeline<W, H, Filter::World::Trails<TRAIL_CAPACITY>,
           Filter::World::Replicate<W>, Filter::World::Orient,
           Filter::Screen::AntiAlias<W, H>>
      filters;

  /**
   * @brief Live slider-backed parameters for the effect.
   * @details trail_ceiling is engine-written (read-only); the rest are sliders.
   */
  struct Params {
    float speed = 2.0f;          /**< Strand travel speed. */
    float gap = 5.0f;            /**< Target spacing between adjacent nodes. */
    int trail_length = 8;        /**< Active trail length. */
    float wipe_duration = 20.0f; /**< Color-wipe transition duration. */
    /** Longest trail the ring can hold, in frames (engine-written). */
    float trail_ceiling = static_cast<float>(TRAIL_LEN_MAX);
  } params;
};
