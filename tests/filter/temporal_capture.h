/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Capture-time plotting — arm-pass solving against an exact column-time oracle
// ============================================================================

/**
 * @brief Orientation history turning about +Y by @p half_turns over one frame.
 * @param start Starting angle, radians.
 * @param half_turns Frame rotation in half-turns.
 * @param steps Substeps after the previous end pose.
 */
inline math::Orientation<> spin_history(float start, float half_turns,
                                        int steps) {
  math::Orientation<> history(math::make_rotation(math::Y_AXIS, start));
  for (int i = 1; i <= steps; ++i)
    history.push(math::make_rotation(
        math::Y_AXIS, start + half_turns * math::PI_F * i / steps));
  return history;
}

/**
 * @brief Arm passes of a moving point found by exhaustive search with exact
 *        `atan2`, independent of `Temporal::capture`.
 * @param position The point at frame fraction u.
 * @return Fractions u in [0, 1) at which `(theta(u) / pi - u)` crosses an
 *         integer.
 */
template <typename PositionFn>
std::vector<float> oracle_captures(PositionFn &&position) {
  constexpr int N = 20000;
  auto offset = [&](float u) {
    const math::Vector p = position(u);
    return std::atan2(p.z, p.x) / math::PI_F - u;
  };
  std::vector<float> found;
  float prev = offset(0.0f);
  if (std::fabs(prev - std::round(prev)) < 1e-6f)
    found.push_back(0.0f);
  for (int i = 1; i <= N; ++i) {
    const float u = static_cast<float>(i) / N;
    const float cur = offset(u);
    const float jump = cur - prev;
    const float unwrapped = prev + std::remainder(jump, 2.0f);
    if (std::floor(prev) != std::floor(unwrapped) && u < 1.0f) {
      const float lo = static_cast<float>(i - 1) / N;
      found.push_back(0.5f * (lo + u));
    }
    prev = cur;
  }
  return found;
}

/**
 * @brief Capture fractions agree with the exhaustive oracle across motion with
 *        and against the sweep, both arm boundaries, and zero and two passes.
 */
inline void test_temporal_capture_matches_exhaustive_oracle() {
  int zero = 0, two = 0, cases = 0;
  for (float drift : {-0.3f, -0.05f, -0.004f, 0.004f, 0.05f, 0.3f}) {
    for (int k = 0; k < 64; ++k) {
      const float start = 2.0f * math::PI_F * k / 64 + 0.013f;
      const math::Orientation<> history = spin_history(start, drift, 3);
      const math::Vector v(0.8f, 0.6f, 0.0f);
      auto position = [&](float u) { return math::rotate(v, history.get(u)); };
      Temporal::Captures got = Temporal::capture(position);
      if (got.count == 2 && got.fraction[1] < got.fraction[0]) {
        std::swap(got.fraction[0], got.fraction[1]);
        std::swap(got.position[0], got.position[1]);
      }
      const std::vector<float> want = oracle_captures(position);
      ++cases;
      zero += want.empty();
      two += want.size() == 2;
      HS_EXPECT_EQ(got.count, static_cast<int>(want.size()));
      for (int i = 0; i < std::min<int>(got.count, want.size()); ++i) {
        HS_EXPECT_NEAR(got.fraction[i], want[i], 2e-3f);
        const math::Vector p = position(got.fraction[i]);
        HS_EXPECT_NEAR(got.position[i].x, p.x, 1e-4f);
        HS_EXPECT_NEAR(got.position[i].z, p.z, 1e-4f);
      }
    }
  }
  HS_EXPECT_GT(zero, 0);
  HS_EXPECT_GT(two, 0);
  HS_EXPECT_GT(cases - zero - two, 0);
}

/**
 * @brief A captured point lands on a column whose composed-path sweep time
 *        shows the point at that same column.
 */
inline void test_temporal_capture_column_matches_composed_time() {
  constexpr int W = 288;
  int checked = 0;
  for (float drift : {-0.2f, 0.02f, 0.2f}) {
    for (int k = 0; k < 48; ++k) {
      const math::Orientation<> history =
          spin_history(2.0f * math::PI_F * k / 48, drift, 3);
      const math::Vector v(0.6f, -0.8f, 0.0f);
      auto position = [&](float u) { return math::rotate(v, history.get(u)); };
      const Temporal::Captures got = Temporal::capture(position);
      for (int i = 0; i < got.count; ++i) {
        const int x =
            math::fast_wrap(static_cast<int>(std::round(
                                math::vector_to_theta<W>(got.position[i]))),
                            W);
        const math::Vector field = position(Temporal::sweep_fraction(x, W));
        const float field_x = math::vector_to_theta<W>(field);
        float gap = std::fabs(field_x - x);
        gap = std::min(gap, W - gap);
        HS_EXPECT_LE(gap, 1.0f);
        ++checked;
      }
    }
  }
  HS_EXPECT_GT(checked, 100);
}

/**
 * @brief Motion faster than `Temporal::MAX_CAPTURE_DRIFT` is emitted once at
 *        the current endpoint.
 */
inline void test_temporal_capture_fast_motion_falls_back_to_endpoint() {
  const math::Orientation<> history = spin_history(0.4f, 0.7f, 3);
  const math::Vector v = math::X_AXIS;
  const Temporal::Captures got = Temporal::capture(
      [&](float u) { return math::rotate(v, history.get(u)); });
  HS_EXPECT_EQ(got.count, 1);
  HS_EXPECT_NEAR(got.fraction[0], 1.0f, 0.0f);
  const math::Vector end = math::rotate(v, history.get());
  HS_EXPECT_NEAR(got.position[0].x, end.x, 1e-6f);
  HS_EXPECT_NEAR(got.position[0].z, end.z, 1e-6f);
}

/**
 * @brief Orient emits each arm pass with age offset by the remaining frame,
 *        and only the current pose with temporal sampling disabled.
 */
inline void test_world_orient_emits_at_capture_time() {
  const math::Orientation<> history = spin_history(-0.5f, -0.08f, 3);
  const math::Vector v(0.0f, 0.0f, 1.0f);
  const Temporal::Captures want = Temporal::capture(
      [&](float u) { return math::rotate(v, history.get(u)); });
  HS_EXPECT_EQ(want.count, 1);
  Filter::World::Orient orient(history);
  std::vector<Tap3D> taps;
  auto record = [&](const math::Vector &p, const Pixel &c, float age,
                    float alpha) { taps.push_back({p, c, age, alpha}); };
  orient.plot(v, Pixel(1, 2, 3), 4.0f, 0.5f, record);
  HS_EXPECT_EQ(taps.size(), 1u);
  HS_EXPECT_LT(want.fraction[0], 0.99f);
  HS_EXPECT_NEAR(taps[0].age, 4.0f + 1.0f - want.fraction[0], 1e-5f);
  HS_EXPECT_NEAR(taps[0].v.x, want.position[0].x, 1e-6f);
  HS_EXPECT_NEAR(taps[0].v.z, want.position[0].z, 1e-6f);
  HS_EXPECT_NEAR(taps[0].alpha, 0.5f, 0.0f);

  Temporal::set_temporal_sampling_enabled(false);
  taps.clear();
  orient.plot(v, Pixel(1, 2, 3), 4.0f, 0.5f, record);
  Temporal::set_temporal_sampling_enabled(true);
  HS_EXPECT_EQ(taps.size(), 1u);
  const math::Vector end = math::rotate(v, history.get());
  HS_EXPECT_NEAR(taps[0].age, 4.0f, 0.0f);
  HS_EXPECT_NEAR(taps[0].v.x, end.x, 1e-6f);
  HS_EXPECT_NEAR(taps[0].v.z, end.z, 1e-6f);
}
