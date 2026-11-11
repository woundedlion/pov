/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Ripple / Noise (the Animation:: classes, distinct from the pure transforms)
// ============================================================================

/**
 * @brief Verifies a Ripple seeds its center, ramps amplitude up then down across
 * its life, pins amplitude to 0 at the duration boundary, and reports done()
 * only on the final frame.
 */
inline void test_ripple_envelope_and_done_boundary() {
  Animation::RippleParams params;
  params.amplitude = 2.0f; // peak captured at construction
  const math::Vector CENTER(0.0f, 1.0f, 0.0f);
  const int DURATION = 20;
  Animation::Ripple ripple(params, CENTER, /*speed=*/0.2f, DURATION);

  // Construction zeroes the live amplitude and seats the center.
  HS_EXPECT_NEAR(params.amplitude, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(params.center.y, 1.0f, 1e-6f);

  std::vector<float> amps;
  float prev_phase = params.phase;
  bool phase_advances = true;
  for (int i = 0; i < DURATION; ++i) {
    ripple.step(fake_canvas());
    amps.push_back(params.amplitude);
    if (i < DURATION - 1 && params.phase <= prev_phase)
      phase_advances = false;
    prev_phase = params.phase;
    if (i < DURATION - 1)
      HS_EXPECT_FALSE(ripple.done());
  }
  HS_EXPECT_TRUE(ripple.done()); // t == duration

  HS_EXPECT_TRUE(phase_advances);

  size_t argmax = 0;
  for (size_t i = 0; i < amps.size(); ++i) {
    HS_CONTEXT("frame", static_cast<long long>(i + 1));
    HS_EXPECT_GE(amps[i], 0.0f);
    HS_EXPECT_LE(amps[i], 2.0f);
    if (amps[i] > amps[argmax])
      argmax = i;
  }
  // Unimodal: strictly up to the peak, strictly down after it.
  for (size_t i = 1; i <= argmax; ++i)
    HS_EXPECT_GT(amps[i], amps[i - 1]);
  for (size_t i = argmax + 1; i < amps.size(); ++i)
    HS_EXPECT_LT(amps[i], amps[i - 1]);
  HS_EXPECT_GT(argmax, (size_t)0);
  HS_EXPECT_LT(argmax, amps.size() - 1);
  HS_EXPECT_NEAR(amps.back(), 0.0f, 1e-6f);
  // Phase stops advancing at the boundary, with the envelope.
  HS_EXPECT_NEAR(params.phase, 0.2f * (DURATION - 1), 1e-5f);
}

/**
 * @brief Noise publishes its initial linear trajectory and remains perpetual.
 */
inline void test_noise_publishes_time_and_is_perpetual() {
  Animation::NoiseParams params;
  Animation::Noise noise(params); // default duration -1
  HS_EXPECT_FALSE(noise.done());
  HS_EXPECT_NEAR(params.time, 0.0f, 1e-6f);

  for (int i = 1; i <= 5; ++i) {
    noise.step(fake_canvas());
    HS_EXPECT_NEAR(params.time, static_cast<float>(i), 1e-6f);
    HS_EXPECT_FALSE(noise.done());
  }

  params.speed = 0.25f;
  noise.step(fake_canvas());
  HS_EXPECT_NEAR(params.time, 5.25f, 1e-6f);
}

template <typename Params> float noise_loop_sample(const Params &params) {
  const math::Vector direction = math::Vector(0.3f, 0.4f, 0.5f).normalized();
  if constexpr (std::is_same_v<Params, Animation::NoiseParams>)
    return noise_transform(direction, params).x;
  else
    return noise_product_field(direction, params);
}

/** @brief Noise samplers follow the loop's linear, turning, and seam coordinates. */
template <typename Params, typename NoiseAnimation>
void test_noise_loop_coordinates_and_samples() {
  struct Fixture {
    float phase;
    float coordinate;
  };
  constexpr Fixture FIXTURES[] = {
      {0.0f, 0.0f},      {960.0f, 960.0f},    {1024.0f, 1000.0f},
      {1088.0f, 960.0f}, {2047.0f, 1.0f},     {2048.0f, 0.0f},
      {2049.0f, -1.0f},  {3072.0f, -1000.0f}, {4095.0f, -1.0f},
      {4096.0f, 0.0f},   {8193.0f, 1.0f}};
  for (const auto &fixture : FIXTURES) {
    for (float sign : {-1.0f, 1.0f}) {
      Params params;
      params.amplitude = 0.5f;
      if constexpr (std::is_same_v<Params, Animation::NoiseParams>)
        params.sync();
      NoiseAnimation animation(params);
      params.speed = sign * fixture.phase;
      animation.step(fake_canvas());
      HS_EXPECT_NEAR(params.time, sign * fixture.coordinate, 1e-5f);
      Params expected;
      expected.amplitude = 0.5f;
      if constexpr (std::is_same_v<Params, Animation::NoiseParams>)
        expected.sync();
      expected.time = sign * fixture.coordinate;
      HS_EXPECT_NEAR(noise_loop_sample(params), noise_loop_sample(expected),
                     1e-6f);
    }
  }

  auto coordinate = [](float phase) {
    Params params;
    NoiseAnimation animation(params);
    params.speed = phase;
    animation.step(fake_canvas());
    return params.time;
  };
  // Inside the turn the coordinate stays between the linear segment's end and
  // the peak, symmetric about the peak and odd in the phase.
  for (float offset : {16.0f, 32.0f, 48.0f}) {
    const float rising = coordinate(1024.0f - offset);
    HS_EXPECT_GT(rising, 960.0f);
    HS_EXPECT_LT(rising, 1000.0f);
    HS_EXPECT_NEAR(coordinate(1024.0f + offset), rising, 1e-5f);
    HS_EXPECT_NEAR(coordinate(-(1024.0f - offset)), -rising, 1e-5f);
    HS_EXPECT_LT(coordinate(1024.0f - offset - 16.0f), rising);
  }
  auto sample = [](float phase) {
    Params params;
    params.amplitude = 0.5f;
    if constexpr (std::is_same_v<Params, Animation::NoiseParams>)
      params.sync();
    NoiseAnimation animation(params);
    params.speed = phase;
    animation.step(fake_canvas());
    return noise_loop_sample(params);
  };
  struct SlopeFixture {
    float phase;
    float slope;
  };
  for (float sign : {-1.0f, 1.0f}) {
    for (const auto &fixture :
         {SlopeFixture{960.0f, 1.0f}, SlopeFixture{1024.0f, 0.0f},
          SlopeFixture{1088.0f, -1.0f}, SlopeFixture{2048.0f, -1.0f}}) {
      const float phase = sign * fixture.phase;
      const float before = coordinate(phase - 1.0f);
      const float current = coordinate(phase);
      const float after = coordinate(phase + 1.0f);
      HS_EXPECT_NEAR((after - before) * 0.5f, fixture.slope, 0.001f);
      if (fixture.phase != 1024.0f)
        HS_EXPECT_NEAR(after - 2.0f * current + before, 0.0f, 0.001f);
      else
        HS_EXPECT_NEAR(after - 2.0f * current + before, sign * (-3.0f / 128.0f),
                       0.0002f);
      constexpr float SAMPLE_STEP = 1.0f / 256.0f;
      const float sample_before = sample(phase - SAMPLE_STEP);
      const float sample_current = sample(phase);
      const float sample_after = sample(phase + SAMPLE_STEP);
      constexpr float SAMPLE_BOUND =
          std::is_same_v<Params, Animation::NoiseParams> ? 0.002f : 0.02f;
      HS_EXPECT_LE(std::fabs(sample_after - sample_before), SAMPLE_BOUND);
      HS_EXPECT_LE(
          std::fabs(sample_after - 2.0f * sample_current + sample_before),
          SAMPLE_BOUND);
    }
  }
}

/** @brief Live speed edits advance the retained noise clock. */
template <typename Params, typename NoiseAnimation>
void test_noise_loop_live_speed() {
  Params params;
  params.time = 7.25f;
  params.speed = 0.25f;
  NoiseAnimation animation(params);
  animation.step(fake_canvas());
  HS_EXPECT_EQ(params.time, 7.5f);
  params.speed = 0.0f;
  animation.step(fake_canvas());
  HS_EXPECT_EQ(params.time, 7.5f);
  params.speed = -0.125f;
  animation.step(fake_canvas());
  HS_EXPECT_EQ(params.time, 7.375f);
}

/** @brief Non-finite speeds hold the last noise publication. */
template <typename Params, typename NoiseAnimation>
void test_noise_loop_nonfinite_speed_holds() {
  Params params;
  params.time = 7.25f;
  params.speed = 0.25f;
  NoiseAnimation animation(params);
  animation.step(fake_canvas());
  for (float speed : {std::numeric_limits<float>::quiet_NaN(),
                      std::numeric_limits<float>::infinity(),
                      -std::numeric_limits<float>::infinity()}) {
    params.speed = speed;
    animation.step(fake_canvas());
    HS_EXPECT_EQ(params.time, 7.5f);
  }
}

/** @brief Publication writes and refreshed parameters preserve the noise clock. */
template <typename Params, typename NoiseAnimation>
void test_noise_loop_published_time_and_refresh() {
  Params params;
  params.time = 7.25f;
  params.speed = 0.25f;
  NoiseAnimation animation(params);
  animation.step(fake_canvas());
  params.speed = -0.125f;
  animation.step(fake_canvas());
  params.time = 123.0f;
  params.speed = 0.0f;
  animation.step(fake_canvas());
  HS_EXPECT_EQ(params.time, 123.0f);
  params.speed = std::numeric_limits<float>::quiet_NaN();
  animation.step(fake_canvas());
  HS_EXPECT_EQ(params.time, 123.0f);
  params.speed = 0.25f;
  animation.step(fake_canvas());
  HS_EXPECT_EQ(params.time, 7.625f);
  Params config;
  config.time = 888.0f;
  config.speed = 0.5f;
  params.refresh_from(config);
  animation.step(fake_canvas());
  HS_EXPECT_EQ(params.time, 8.125f);
  animation.rewind();
  animation.step(fake_canvas());
  HS_EXPECT_EQ(params.time, 8.625f);
}

/** @brief Noise copies and assignments preserve the clock without allocation. */
template <typename Params, typename NoiseAnimation>
void test_noise_loop_copy_and_assignment() {
  Params params;
  params.time = 7.25f;
  params.speed = 0.25f;
  NoiseAnimation animation(params);
  animation.step(fake_canvas());
  const size_t before_copy = persistent_arena.get_offset();
  NoiseAnimation copied(animation);
  animation.step(fake_canvas());
  const float expected_copy = params.time;
  copied.step(fake_canvas());
  HS_EXPECT_EQ(params.time, expected_copy);
  NoiseAnimation assigned(params);
  assigned = animation;
  animation.step(fake_canvas());
  const float expected_assignment = params.time;
  assigned.step(fake_canvas());
  HS_EXPECT_EQ(params.time, expected_assignment);
  HS_EXPECT_EQ(persistent_arena.get_offset(), before_copy);
}

/** @brief Non-finite and extreme seeds recover to a finite noise clock. */
template <typename Params, typename NoiseAnimation>
void test_noise_loop_seed_recovery() {
  for (float seed : {std::numeric_limits<float>::quiet_NaN(),
                     std::numeric_limits<float>::infinity(),
                     -std::numeric_limits<float>::infinity()}) {
    Params seeded;
    seeded.time = seed;
    seeded.speed = 0.0f;
    NoiseAnimation seeded_animation(seeded);
    seeded_animation.step(fake_canvas());
    if (std::isnan(seed))
      HS_EXPECT_TRUE(std::isnan(seeded.time));
    else
      HS_EXPECT_EQ(seeded.time, seed);
    seeded.speed = 0.125f;
    seeded_animation.step(fake_canvas());
    HS_EXPECT_EQ(seeded.time, 0.125f);
  }
  for (float seed : {std::numeric_limits<float>::max(),
                     -std::numeric_limits<float>::max()}) {
    Params seeded;
    seeded.time = seed;
    seeded.speed = 0.0f;
    NoiseAnimation seeded_animation(seeded);
    seeded_animation.step(fake_canvas());
    HS_EXPECT_EQ(seeded.time, seed);
    seeded.speed = 0.125f;
    seeded_animation.step(fake_canvas());
    HS_EXPECT_EQ(seeded.time, 0.125f);
    seeded.speed = seed;
    seeded_animation.step(fake_canvas());
    HS_EXPECT_EQ(seeded.time, 0.125f);
  }
}

/** @brief Default noise clocks continue moving after long continuous uptime. */
template <typename Params, typename NoiseAnimation>
void test_noise_loop_long_uptime() {
  Params params;
  const float speed = params.speed;
  NoiseAnimation animation(params);
  constexpr uint32_t FRAMES = (1u << 24) + 2;
  for (uint32_t frame = 0; frame < FRAMES; ++frame)
    animation.step(fake_canvas());
  double expected = std::fmod(static_cast<double>(FRAMES) * speed, 4096.0);
  if (expected >= 2048.0)
    expected -= 4096.0;
  HS_EXPECT_LE(std::fabs(expected), 960.0);
  HS_EXPECT_NEAR(params.time, static_cast<float>(expected), 1e-4f);
  const float aged = params.time;
  animation.step(fake_canvas());
  HS_EXPECT_NEAR(params.time, static_cast<float>(expected + speed), 1e-4f);
  HS_EXPECT_TRUE(params.time != aged);
}

/** @brief Noise's finite duration still completes on its configured frame. */
inline void test_noise_loop_preserves_finite_duration() {
  Animation::NoiseParams params;
  Animation::Noise animation(params, 2);
  animation.step(fake_canvas());
  HS_EXPECT_FALSE(animation.done());
  HS_EXPECT_EQ(params.time, 1.0f);
  animation.step(fake_canvas());
  HS_EXPECT_TRUE(animation.done());
  HS_EXPECT_EQ(params.time, 2.0f);
}

/** Frames the RandomWalk oracles drive. */
constexpr int RANDOM_WALK_FRAMES = 50;

/**
 * @brief A seeded RandomWalk moves the orientation as a rigid rotation, one
 * bounded step per frame, and follows its seed.
 * @details The orientation is a rotation, so the angle between two probes it
 * carries is invariant and no frame turns a probe farther than the configured
 * speed. The seed pair is the negative control: the same seed replays the
 * trajectory, a different one leaves it.
 */
inline void test_random_walk_stays_unit_and_travels() {
  using Walk = Animation::RandomWalk<288, 4>;
  const Walk::Options options = Walk::Options::Energetic();

  // Records one seeded walk's trace of the probe pair, asserting the rigid-body
  // invariants as it goes.
  auto run = [&](int seed, std::vector<math::Vector> &trace) {
    math::Orientation<4> o; // identity
    FastNoiseLite noise;
    Walk walk(o, math::Y_AXIS, noise, options, seed);
    const math::Vector PROBE = math::X_AXIS, COMPANION = math::Z_AXIS;
    const float rest =
        static_cast<float>(chord_angle_between(PROBE, COMPANION));
    math::Vector prev = o.orient(PROBE);
    float travel = 0.0f;
    for (int fr = 0; fr < RANDOM_WALK_FRAMES; ++fr) {
      HS_CONTEXT("frame", fr);
      walk.step(fake_canvas());
      const math::Vector cur = o.orient(PROBE);
      const math::Vector other = o.orient(COMPANION);
      trace.push_back(cur);
      HS_EXPECT_NEAR(cur.length(), 1.0f, 1e-4f);
      HS_EXPECT_NEAR(other.length(), 1.0f, 1e-4f);
      // Rigid: the probe pair keeps its opening angle.
      HS_EXPECT_NEAR(chord_angle_between(cur, other), rest, 1e-4f);
      const float step = static_cast<float>(chord_angle_between(prev, cur));
      HS_EXPECT_LE(step, options.speed + 1e-4f);
      travel += step;
      prev = cur;
    }
    return travel;
  };

  std::vector<math::Vector> trace, replay, divergent;
  const float travel = run(1234, trace);
  HS_EXPECT_GT(travel, 0.01f);
  HS_EXPECT_LE(travel, options.speed * RANDOM_WALK_FRAMES);

  HS_EXPECT_NEAR(run(1234, replay), travel, 0.0f);
  run(4321, divergent);
  HS_EXPECT_SIZE_OR_RETURN(replay, trace.size());
  HS_EXPECT_SIZE_OR_RETURN(divergent, trace.size());
  double parted = 0.0;
  for (size_t i = 0; i < trace.size(); ++i) {
    HS_CONTEXT("frame", static_cast<long long>(i));
    HS_EXPECT_EQ(replay[i].x, trace[i].x);
    HS_EXPECT_EQ(replay[i].y, trace[i].y);
    HS_EXPECT_EQ(replay[i].z, trace[i].z);
    parted = std::max(parted, chord_angle_between(divergent[i], trace[i]));
  }
  HS_EXPECT_GT(parted, 0.01);
}
