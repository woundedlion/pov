/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Mobius warps
// ----------------------------------------------------------------------------
// MobiusWarp and MobiusWarpCircular drive b along an eased angle closing at 2π.
// MobiusWarpEvolving modulates all eight coefficients perpetually.
// ============================================================================

/**
 * @brief Verifies MobiusWarp leaves the origin, closes near b == 0 at
 * completion, and reports done() only on the final frame.
 */
inline void test_mobiuswarp_closes_at_completion() {
  math::MobiusParams params;
  const float scale = 0.4f;
  const int duration = 8;
  Animation::MobiusWarp warp(params, scale, duration, /*repeat=*/false,
                             math::ease_linear);
  HS_EXPECT_FALSE(warp.done());

  warp.step(fake_canvas()); // t=1: b lifts off the origin
  HS_EXPECT_GT(std::abs(params.b.re) + std::abs(params.b.im), 1e-4f);

  for (int i = 1; i < duration - 1; ++i)
    warp.step(fake_canvas());
  HS_EXPECT_FALSE(warp.done());

  warp.step(fake_canvas()); // t=duration: angle 2π -> b back to 0
  HS_EXPECT_TRUE(warp.done());
  HS_EXPECT_NEAR(params.b.re, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(params.b.im, 0.0f, 1e-4f);
}

/** @brief Verifies a NaN live scale leaves MobiusWarp using its last finite value. */
inline void test_mobiuswarp_retains_last_finite_scale() {
  math::MobiusParams params, reference;
  float live = 0.4f;
  Animation::MobiusWarp warp(params, 0.1f, 8, false, math::ease_linear);
  Animation::MobiusWarp expected(reference, 0.4f, 8, false, math::ease_linear);
  warp.bind_scale(live);
  warp.step(fake_canvas());
  expected.step(fake_canvas());
  live = std::numeric_limits<float>::quiet_NaN();
  warp.step(fake_canvas());
  expected.step(fake_canvas());
  HS_EXPECT_EQ(params.b.re, reference.b.re);
  HS_EXPECT_EQ(params.b.im, reference.b.im);
}

/**
 * @brief Verifies bind_scale makes step() read the live referent instead of the
 * captured construction-time scale.
 */
inline void test_mobiuswarp_bind_scale_reads_live() {
  math::MobiusParams params;
  float live = 1.0f;
  const int duration = 4;
  Animation::MobiusWarp warp(params, /*scale=*/0.0f, duration, /*repeat=*/false,
                             math::ease_linear);
  warp.bind_scale(live);
  warp.step(
      fake_canvas()); // captured scale is 0, so any motion comes from live
  HS_EXPECT_GT(std::abs(params.b.re) + std::abs(params.b.im), 1e-4f);
}

/**
 * @brief Verifies MobiusWarpCircular traces param b on the |b| == scale circle,
 * landing at (scale, 0) at completion and reporting the done() boundary.
 */
inline void test_mobiuswarp_circular_traces_radius() {
  math::MobiusParams params;
  const float scale = 0.3f;
  const int duration = 8;
  Animation::MobiusWarpCircular warp(params, scale, duration, /*repeat=*/false,
                                     math::ease_linear);
  warp.step(fake_canvas()); // |b| sits on the scale-radius circle every frame
  HS_EXPECT_NEAR(
      std::sqrt(params.b.re * params.b.re + params.b.im * params.b.im), scale,
      1e-4f);

  for (int i = 1; i < duration; ++i) {
    warp.step(fake_canvas());
    if (i == 1) {
      HS_EXPECT_NEAR(params.b.re, 0.0f, 1e-4f);
      HS_EXPECT_NEAR(params.b.im, -scale, 1e-4f);
    } else if (i == 3) {
      HS_EXPECT_NEAR(params.b.re, -scale, 1e-4f);
      HS_EXPECT_NEAR(params.b.im, 0.0f, 1e-4f);
    } else if (i == 5) {
      HS_EXPECT_NEAR(params.b.re, 0.0f, 1e-4f);
      HS_EXPECT_NEAR(params.b.im, scale, 1e-4f);
    }
  }
  HS_EXPECT_TRUE(warp.done());
  // angle 2π: b.re == scale, b.im == 0.
  HS_EXPECT_NEAR(params.b.re, scale, 1e-4f);
  HS_EXPECT_NEAR(params.b.im, 0.0f, 1e-4f);
}

/**
 * @brief Verifies MobiusWarpCircular's bind_scale makes step() read the live
 * referent instead of the captured construction-time scale.
 */
inline void test_mobiuswarp_circular_bind_scale_reads_live() {
  math::MobiusParams params;
  float live = 0.5f;
  const int duration = 4;
  Animation::MobiusWarpCircular warp(params, /*scale=*/0.0f, duration,
                                     /*repeat=*/false, math::ease_linear);
  warp.bind_scale(live);
  warp.step(fake_canvas()); // captured scale is 0: any radius comes from live
  HS_EXPECT_NEAR(
      std::sqrt(params.b.re * params.b.re + params.b.im * params.b.im), live,
      1e-4f);
}

/**
 * @brief Verifies MobiusWarpEvolving modulates all eight coefficients within
 * ±scale of their captured baseline, drives each of them, and, being perpetual,
 * never reports done().
 */
inline void test_mobiuswarp_evolving_bounded_and_perpetual() {
  const auto saved_rng = hs::random();
  hs::random().seed(1337);
  math::MobiusParams params(2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 3.0f, 0.0f);
  const float scale = 0.5f;
  Animation::MobiusWarpEvolving warp(params, scale, /*speed=*/0.05f);
  HS_EXPECT_FALSE(warp.done());

  const math::MobiusParams base = params;
  const float baseline[8] = {base.a.re, base.a.im, base.b.re, base.b.im,
                             base.c.re, base.c.im, base.d.re, base.d.im};
  float peak[8] = {};
  for (int i = 0; i < 50; ++i) {
    warp.step(fake_canvas());
    HS_EXPECT_FALSE(warp.done()); // perpetual: duration -1
    const float live[8] = {params.a.re, params.a.im, params.b.re, params.b.im,
                           params.c.re, params.c.im, params.d.re, params.d.im};
    for (int c = 0; c < 8; ++c) {
      float delta = std::abs(live[c] - baseline[c]);
      HS_EXPECT_LE(delta, scale + 1e-4f);
      peak[c] = std::max(peak[c], delta);
    }
  }
  // 50 frames at speed 0.05 sweeps even the slowest channel past half amplitude.
  for (int c = 0; c < 8; ++c)
    HS_EXPECT_GT(peak[c], 0.5f * scale);
  hs::random() = saved_rng;
}

/** @brief Wrapped channels follow live speeds after large temporal advances. */
inline void test_mobiuswarp_evolving_wrapped_live_channels() {
  const auto saved_rng = hs::random();
  hs::random().seed(1337);
  math::MobiusParams params;
  const math::MobiusParams base = params;
  const size_t allocated_before = persistent_arena.get_offset();
  Animation::MobiusWarpEvolving warp(params, 0.5f, 0.05f);
  const size_t allocated = persistent_arena.get_offset() - allocated_before;
  HS_EXPECT_GE(allocated, 8 * sizeof(double));
  HS_EXPECT_LE(allocated, 8 * sizeof(double) + alignof(double) - 1);
  const float baseline[] = {base.a.re, base.a.im, base.b.re, base.b.im,
                            base.c.re, base.c.im, base.d.re, base.d.im};
  constexpr float FREQUENCIES[] = {1.0f,  1.13f, 1.27f, 1.39f,
                                   0.71f, 0.83f, 0.97f, 1.09f};
  constexpr double PERIOD = 6.283185307179586476925286766559;
  double phases[8] = {};
  float effective_speed = 0.05f;
  for (float speed :
       {0.05f, 0.05f, 1.0f, 1000000000.0f, 0.05f, -0.03f, 0.0f,
        std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(),
        0.05f, std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::infinity()}) {
    warp.set_speed(speed);
    if (std::isfinite(speed))
      effective_speed = speed;
    warp.step(fake_canvas());
    const float observed[] = {params.a.re, params.a.im, params.b.re,
                              params.b.im, params.c.re, params.c.im,
                              params.d.re, params.d.im};
    for (int i = 0; i < 8; ++i) {
      const double increment = static_cast<double>(effective_speed) *
                               static_cast<double>(FREQUENCIES[i]);
      phases[i] = std::fmod(phases[i] + std::fmod(increment, PERIOD), PERIOD);
      if (phases[i] < 0.0)
        phases[i] += PERIOD;
      const float angle = static_cast<float>(phases[i]) + warp.phase(i);
      const float expected =
          baseline[i] + (i % 2 ? cosf(angle) : sinf(angle)) * 0.5f;
      HS_EXPECT_TRUE(std::isfinite(observed[i]));
      HS_EXPECT_NEAR(observed[i], expected, 1e-5f);
    }
  }
  hs::random() = saved_rng;
}

/** @brief Copies keep independent phases; moves transfer their phase blocks. */
inline void test_mobiuswarp_evolving_value_semantics() {
  const auto saved_rng = hs::random();
  hs::random().seed(1337);
  math::MobiusParams params;
  int completions = 0;
  Animation::MobiusWarpEvolving original(params, 0.5f, 0.05f);
  original.then([&]() { ++completions; });
  original.step(fake_canvas());
  const size_t before_copy = persistent_arena.get_offset();
  Animation::MobiusWarpEvolving copied(original);
  HS_EXPECT_GE(persistent_arena.get_offset() - before_copy, 8 * sizeof(double));
  HS_EXPECT_LE(persistent_arena.get_offset() - before_copy,
               8 * sizeof(double) + alignof(double) - 1);
  original.step(fake_canvas());
  const math::MobiusParams expected_copy = params;
  copied.step(fake_canvas());
  HS_EXPECT_TRUE(std::memcmp(&params, &expected_copy, sizeof(params)) == 0);

  Animation::MobiusWarpEvolving assigned(params, 0.1f, 0.2f);
  const size_t before_assignment = persistent_arena.get_offset();
  assigned = original;
  HS_EXPECT_EQ(persistent_arena.get_offset(), before_assignment);
  original.step(fake_canvas());
  const math::MobiusParams expected_assignment = params;
  assigned.step(fake_canvas());
  HS_EXPECT_TRUE(std::memcmp(&params, &expected_assignment, sizeof(params)) ==
                 0);
  const auto *self = &assigned;
  assigned = *self;
  HS_EXPECT_EQ(persistent_arena.get_offset(), before_assignment);
  assigned.rewind();
  original.step(fake_canvas());
  const math::MobiusParams expected_rewind = params;
  assigned.step(fake_canvas());
  HS_EXPECT_TRUE(std::memcmp(&params, &expected_rewind, sizeof(params)) == 0);

  Animation::MobiusWarpEvolving donor(params, 0.5f, 0.05f);
  const size_t before_move = persistent_arena.get_offset();
  Animation::MobiusWarpEvolving moved(std::move(copied));
  HS_EXPECT_EQ(persistent_arena.get_offset(), before_move);
  copied = donor;
  HS_EXPECT_GE(persistent_arena.get_offset() - before_move, 8 * sizeof(double));
  const size_t before_move_assignment = persistent_arena.get_offset();
  assigned = std::move(moved);
  HS_EXPECT_EQ(persistent_arena.get_offset(), before_move_assignment);
  moved = donor;
  HS_EXPECT_GE(persistent_arena.get_offset() - before_move_assignment,
               8 * sizeof(double));
  assigned.step(fake_canvas());
  HS_EXPECT_TRUE(std::memcmp(&params, &expected_assignment, sizeof(params)) ==
                 0);
  original.cancel();
  original.post_callback();
  HS_EXPECT_EQ(completions, 1);
  assigned.cancel();
  assigned.post_callback();
  HS_EXPECT_EQ(completions, 2);
  copied.post_callback();
  moved.post_callback();
  HS_EXPECT_EQ(completions, 2);
  hs::random() = saved_rng;
}

/** @brief Default modulation still moves after the float clock's uptime limit. */
inline void test_mobiuswarp_evolving_long_uptime() {
  const auto saved_rng = hs::random();
  hs::random().seed(1337);
  math::MobiusParams params;
  const math::MobiusParams base = params;
  Animation::MobiusWarpEvolving warp(params);
  constexpr uint32_t FRAMES = (1u << 25) + 3;
  for (uint32_t frame = 0; frame < FRAMES; ++frame)
    warp.step(fake_canvas());
  const math::MobiusParams aged = params;
  warp.step(fake_canvas());
  HS_EXPECT_TRUE(std::memcmp(&params, &aged, sizeof(params)) != 0);
  constexpr float FREQUENCIES[] = {1.0f,  1.13f, 1.27f, 1.39f,
                                   0.71f, 0.83f, 0.97f, 1.09f};
  const float baseline[] = {base.a.re, base.a.im, base.b.re, base.b.im,
                            base.c.re, base.c.im, base.d.re, base.d.im};
  const float observed[] = {params.a.re, params.a.im, params.b.re, params.b.im,
                            params.c.re, params.c.im, params.d.re, params.d.im};
  constexpr double PERIOD = 6.283185307179586476925286766559;
  for (int i = 0; i < 8; ++i) {
    const double accumulated = static_cast<double>(FRAMES + 1) *
                               static_cast<double>(0.01f) *
                               static_cast<double>(FREQUENCIES[i]);
    const float angle =
        static_cast<float>(std::fmod(accumulated, PERIOD)) + warp.phase(i);
    const float expected =
        baseline[i] + (i % 2 ? cosf(angle) : sinf(angle)) * 0.5f;
    HS_EXPECT_NEAR(observed[i], expected, 1e-5f);
  }
  hs::random() = saved_rng;
}
