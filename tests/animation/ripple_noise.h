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
  const math::Vector center(0.0f, 1.0f, 0.0f);
  const int duration = 20;
  Animation::Ripple ripple(params, center, /*speed=*/0.2f, duration);

  // Construction zeroes the live amplitude and seats the center.
  HS_EXPECT_NEAR(params.amplitude, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(params.center.y, 1.0f, 1e-6f);

  std::vector<float> amps;
  float prev_phase = params.phase;
  bool phase_advances = true;
  for (int i = 0; i < duration; ++i) {
    ripple.step(fake_canvas());
    amps.push_back(params.amplitude);
    if (i < duration - 1 && params.phase <= prev_phase)
      phase_advances = false;
    prev_phase = params.phase;
    if (i < duration - 1)
      HS_EXPECT_FALSE(ripple.done());
  }
  HS_EXPECT_TRUE(ripple.done()); // t == duration

  HS_EXPECT_TRUE(phase_advances);

  // Reference envelope: a parabolic ease-in over the first tenth of the life,
  // scaled by a linear decay to zero at the duration boundary.
  auto reference = [&](int frame) {
    if (frame >= duration)
      return 0.0f;
    const float progress = static_cast<float>(frame) / duration;
    const float attack = std::min(progress / 0.1f, 1.0f);
    return 2.0f * attack * attack * (1.0f - progress);
  };
  size_t argmax = 0;
  for (size_t i = 0; i < amps.size(); ++i) {
    HS_CONTEXT("frame", static_cast<long long>(i + 1));
    HS_EXPECT_NEAR(amps[i], reference(static_cast<int>(i) + 1), 1e-6f);
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
  HS_EXPECT_NEAR(params.phase, 0.2f * (duration - 1), 1e-5f);
}

/**
 * @brief Verifies the Noise animation integrates speed into params.time each
 * step and, being perpetual, never reports done().
 * @details A mid-run speed edit carries the phase on from where it was.
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
    const math::Vector probe = math::X_AXIS, companion = math::Z_AXIS;
    const float rest =
        static_cast<float>(small_angle_between(probe, companion));
    math::Vector prev = o.orient(probe);
    float travel = 0.0f;
    for (int fr = 0; fr < RANDOM_WALK_FRAMES; ++fr) {
      HS_CONTEXT("frame", fr);
      walk.step(fake_canvas());
      const math::Vector cur = o.orient(probe);
      const math::Vector other = o.orient(companion);
      trace.push_back(cur);
      HS_EXPECT_NEAR(cur.length(), 1.0f, 1e-4f);
      HS_EXPECT_NEAR(other.length(), 1.0f, 1e-4f);
      // Rigid: the probe pair keeps its opening angle.
      HS_EXPECT_NEAR(small_angle_between(cur, other), rest, 1e-4f);
      const float step = static_cast<float>(small_angle_between(prev, cur));
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
    parted = std::max(parted, small_angle_between(divergent[i], trace[i]));
  }
  HS_EXPECT_GT(parted, 0.01);
}

/**
 * @brief Compares rotation updates and their composition under identical inputs.
 */
inline void test_random_walk_stable_rotation_matches_same_state() {
  constexpr int FRAMES = 50;
  constexpr double ANGULAR_ERROR_RADIANS = 1e-6;
  constexpr double COMPOSED_DRIFT_RADIANS = 1e-4;
  const Animation::RandomWalkOptions options =
      Animation::RandomWalkOptions::Energetic();

  FastNoiseLite default_noise;
  FastNoiseLite stable_noise;
  for (FastNoiseLite *noise : {&default_noise, &stable_noise}) {
    noise->SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise->SetFrequency(options.noise_scale);
    noise->SetSeed(1234);
  }
  math::Vector default_position = math::Y_AXIS;
  math::Vector default_direction = math::perpendicular_axis(math::Y_AXIS);
  float default_velocity = 0.0f;
  math::Vector default_probe = math::X_AXIS;
  math::Vector stable_probe = default_probe;

  for (uint32_t frame = 1; frame <= FRAMES; ++frame) {
    HS_CONTEXT("frame", static_cast<long long>(frame));
    math::Vector stable_position = default_position;
    math::Vector stable_direction = default_direction;
    float stable_velocity = default_velocity;
    const Animation::RandomWalkDelta expected =
        Animation::step_random_walk<false>(default_position, default_direction,
                                           default_velocity, default_noise,
                                           options, frame);
    const Animation::RandomWalkDelta actual = Animation::step_random_walk<true>(
        stable_position, stable_direction, stable_velocity, stable_noise,
        options, frame);
    HS_EXPECT_LE(small_angle_between(stable_position, default_position),
                 ANGULAR_ERROR_RADIANS);
    HS_EXPECT_LE(small_angle_between(stable_direction, default_direction),
                 ANGULAR_ERROR_RADIANS);
    HS_EXPECT_LE(small_angle_between(actual.axis, expected.axis),
                 ANGULAR_ERROR_RADIANS);
    default_probe = math::rotate(default_probe, expected.rotation).normalized();
    stable_probe = math::rotate(stable_probe, actual.rotation).normalized();
    HS_EXPECT_LE(small_angle_between(stable_probe, default_probe),
                 COMPOSED_DRIFT_RADIANS);
    HS_EXPECT_NEAR(stable_position.length(), 1.0f, 1e-4f);
  }
}
