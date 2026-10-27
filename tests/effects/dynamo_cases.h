/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Dynamo white-box cases. */

/**
 * @brief White-box accessor for Dynamo's leader drag, palette wipes and trail
 * ring.
 */
struct DynamoWhiteBox {
  /** @brief Lowering the gap closes slack toward the nearest leader column. */
  static void check_gap_reduction_uses_shortest_path() {
    reset_effect_globals();
    Dynamo<DEFAULT_W, DEFAULT_H> effect;
    effect.params.gap = 2;
    for (int direction : {-1, 1}) {
      typename Dynamo<DEFAULT_W, DEFAULT_H>::Node leader, follower;
      leader.x = 0;
      leader.v = direction;
      follower.x = math::wrap(direction * 6, DEFAULT_W);
      follower.v = direction;
      effect.drag(leader, follower);
      HS_EXPECT_EQ(follower.x, math::wrap(direction * 2, DEFAULT_W));
      HS_EXPECT_EQ(follower.v, -direction);
    }
  }

  /**
   * @brief Verifies color() is memory-safe and in-range under inverted bands.
   * @details Stages inverted overlapping-wipe boundaries and pins finite alpha
   * in [0,1] across the angular span. Palette indexing is bounded by the boundary
   * container's capacity.
   */
  static void check_overlapping_wipes_stay_in_range() {
    // init() bakes from persistent_arena and schedules on the global timeline.
    reset_effect_globals();

    Dynamo<DEFAULT_W, DEFAULT_H> effect;
    effect.init();

    const Color4 opening = effect.color(math::Z_AXIS, 0.5f);
    effect.color_wipe();
    HS_EXPECT_EQ(effect.color(math::Z_AXIS, 0.5f).color, opening.color);
    effect.palette_boundaries.front() = math::PI_F + effect.WIPE_BLEND_WIDTH;
    HS_EXPECT_EQ(effect.color(math::Z_AXIS * -1.0f, 0.5f).color,
                 effect.baked_palettes[0].get(0.5f).color);
    effect.color_wipe();
    HS_EXPECT_SIZE_OR_RETURN(effect.palette_boundaries, 2);

    // The newer band (index 0) overtakes the older (index 1): non-monotonic
    // order. These magnitudes also reach the second boundary and the
    // baked_palettes[i+1] access for part of the sweep.
    effect.palette_boundaries[0] = 1.5f; // newer, overtaken ahead of the older
    effect.palette_boundaries[1] = 0.5f; // older, left behind

    // PALETTE_NORMAL is Z_AXIS, so v = (sin theta, 0, cos theta) sweeps
    // angle_between(v, normal) across the full [0, PI] band span.
    constexpr int STEPS = 256;
    for (int i = 0; i <= STEPS; ++i) {
      float theta =
          math::PI_F * static_cast<float>(i) / static_cast<float>(STEPS);
      math::Vector v(std::sin(theta), 0.0f, std::cos(theta));
      Color4 c = effect.color(v, 0.5f);
      HS_EXPECT_TRUE(std::isfinite(c.alpha));
      HS_EXPECT_GE(c.alpha, 0.0f);
      HS_EXPECT_LE(c.alpha, 1.0f);
    }
  }

  /**
   * @brief Verifies color() picks each band's palette under two ordered
   *        boundaries and crossfades to the lerp16 midpoint on each boundary.
   */
  static void check_wipe_band_colors() {
    reset_effect_globals();

    Dynamo<DEFAULT_W, DEFAULT_H> effect;
    effect.init();
    effect.color_wipe();
    effect.color_wipe();
    HS_EXPECT_SIZE_OR_RETURN(effect.palette_boundaries, 2);

    const auto direction = [](float theta) {
      return math::Vector(std::sin(theta), 0.0f, std::cos(theta));
    };
    const auto angle = [&](float theta) {
      return math::fast_acos(hs::clamp(
          math::dot(direction(theta), effect.PALETTE_NORMAL), -1.0f, 1.0f));
    };
    constexpr float FIRST = 0.9f, SECOND = 2.6f;
    effect.palette_boundaries[0] = angle(FIRST);
    effect.palette_boundaries[1] = angle(SECOND);
    HS_EXPECT_LT(effect.palette_boundaries[0] + effect.WIPE_BLEND_WIDTH,
                 effect.palette_boundaries[1] - effect.WIPE_BLEND_WIDTH);

    constexpr float T = 0.5f;
    Pixel bands[3];
    for (int i = 0; i < 3; ++i)
      bands[i] = effect.baked_palettes[i].get(T).color;
    HS_EXPECT_NE(bands[0], bands[1]);
    HS_EXPECT_NE(bands[1], bands[2]);

    HS_EXPECT_EQ(effect.color(math::Z_AXIS, T).color, bands[0]);
    HS_EXPECT_EQ(effect.color(direction(1.75f), T).color, bands[1]);
    const float boundaries[2] = {FIRST, SECOND};
    for (int i = 0; i < 2; ++i) {
      const Pixel mid = effect.color(direction(boundaries[i]), T).color;
      const Pixel expected =
          bands[i].lerp16(bands[i + 1], float_to_pixel16(0.5f));
      HS_EXPECT_LE(max_channel_gap(mid, expected), 1);
    }
  }

  using D = Dynamo<DEFAULT_W, DEFAULT_H>;
  using Ring = Filter::World::Trails<D::TRAIL_CAPACITY>;
  static constexpr int trail_capacity() { return D::TRAIL_CAPACITY; }
  static constexpr int trail_len_max() { return D::TRAIL_LEN_MAX; }
  static constexpr float speed_max() { return D::SPEED_MAX; }
  static void set_speed(D &d, float v) { d.params.speed = v; }
  static void set_trail_length(D &d, float v) { d.params.trail_length = v; }
  static float trail_ceiling(const D &d) { return d.params.trail_ceiling; }
  static uint32_t emitted_points(const D &d) { return d.emitted_points; }
  static uint32_t points_per_emission(const D &d) {
    return d.points_per_emission;
  }
  static size_t trail_points(D &d) { return d.filters.get<Ring>().size(); }
};

/**
 * @brief Pins Dynamo's trail-ring ceiling as the thing that keeps the ring from
 *        evicting.
 * @details Past trail_length_ceiling() the ring would evict points of
 *          arbitrary age. At both sliders' maxima the requested trail would
 *          overrun while steady-state occupancy stays inside the ring.
 */
inline void test_dynamo_trail_ceiling_bounds_the_ring() {
  using WB = DynamoWhiteBox;
  reset_effect_globals();

  WB::D effect;
  effect.init();

  constexpr float MAX_SPEED = WB::speed_max();
  constexpr float MAX_TRAIL = WB::trail_len_max();
  WB::set_speed(effect, MAX_SPEED);
  WB::set_trail_length(effect, MAX_TRAIL);

  const auto capacity = static_cast<uint32_t>(WB::trail_capacity());
  for (int f = 0; f < static_cast<int>(MAX_TRAIL) + 2; ++f) {
    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_LT(WB::trail_points(effect), static_cast<size_t>(capacity));
  }

  // Points buffered per frame at this speed: points_per_emission per whole step,
  // and speed_accumulator carries at most one extra step into a frame.
  const uint32_t per_frame =
      WB::points_per_emission(effect) * (static_cast<uint32_t>(MAX_SPEED) + 1);
  HS_EXPECT_GT(per_frame, 0u);
  HS_EXPECT_GT(per_frame * static_cast<uint32_t>(MAX_TRAIL), capacity);
  HS_EXPECT_GT(WB::trail_points(effect), size_t{0});
  HS_EXPECT_LT(WB::trail_points(effect), static_cast<size_t>(capacity));
}

/**
 * @brief Pins what Dynamo's emitted_points counter actually measures.
 * @details draw_nodes() bumps it from the strand's fragment shader, so the
 *          trail-ring ceiling holds only while the rasterizer shades exactly
 *          once per point it hands the pipeline. With nothing aging out, the
 *          ring's growth across a frame must equal emitted_points.
 */
inline void test_dynamo_emitted_points_counts_ring_seeds() {
  using WB = DynamoWhiteBox;
  reset_effect_globals();

  WB::D effect;
  effect.init();
  WB::set_trail_length(effect, 100.0f);
  WB::set_speed(effect, 1.0f);
  constexpr int FRAMES = 3;

  for (int f = 0; f < FRAMES; ++f) {
    const size_t before = WB::trail_points(effect);
    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_GT(WB::trail_ceiling(effect), static_cast<float>(FRAMES));
    const size_t seeded = WB::trail_points(effect) - before;
    HS_EXPECT_GT(WB::emitted_points(effect), 0u);
    HS_EXPECT_EQ(seeded, static_cast<size_t>(WB::emitted_points(effect)));
  }
}
