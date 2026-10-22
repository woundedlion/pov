/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief An arm pass found by brute force. */
struct TruePass {
  double time; ///< Frame fraction.
  int arm;     ///< Branch whose column `(arm + time) W / 2` meets the track.
};

/**
 * @brief Arm passes of a column track by dense sampling and bisection.
 * @tparam W Canvas width.
 * @tparam Column Callable `double(double t)` giving the unwrapped column.
 * @param column The track.
 * @param t0,t1 Span to search.
 * @return Every root of `column(t) - (n + t) W / 2` over the span.
 */
template <int W, typename Column>
inline std::vector<TruePass> brute_force_passes(Column &&column, double t0,
                                                double t1) {
  constexpr double HALF = W / 2.0;
  constexpr int STEPS = 4096;
  std::vector<TruePass> passes;
  const auto relative = [&](double t) { return column(t) / HALF - t; };
  double lo_t = t0, lo = relative(t0);
  for (int i = 1; i <= STEPS; ++i) {
    const double hi_t = t0 + (t1 - t0) * i / STEPS;
    const double hi = relative(hi_t);
    for (int n = static_cast<int>(std::floor(std::min(lo, hi)));
         n <= static_cast<int>(std::floor(std::max(lo, hi))); ++n) {
      if (!((lo - n) * (hi - n) <= 0.0) || hi == n)
        continue;
      double a = lo_t, b = hi_t;
      for (int k = 0; k < 60; ++k) {
        const double mid = 0.5 * (a + b);
        if ((relative(a) - n) * (relative(mid) - n) <= 0.0)
          b = mid;
        else
          a = mid;
      }
      passes.push_back({0.5 * (a + b), n});
    }
    lo_t = hi_t;
    lo = hi;
  }
  return passes;
}

/**
 * @brief Collects a solver's passes.
 * @return The passes in emission order.
 */
struct PassList {
  std::vector<Temporal::ArmPass> passes; ///< Emitted passes.
  /** @brief Records a pass.
   * @param pass The pass. */
  void operator()(const Temporal::ArmPass &pass) { passes.push_back(pass); }
};

/**
 * @brief Checks closed-form passes of a linear track against brute force.
 * @param x0,x1 Unwrapped columns at @p t0 and @p t1.
 * @param t0,t1 Span.
 * @param expected Number of passes.
 */
inline void check_linear_passes(float x0, float x1, float t0, float t1,
                                int expected) {
  constexpr int W = 288;
  HS_CONTEXT("track columns x100", std::llround(x0 * 100.0f),
             std::llround(x1 * 100.0f));
  HS_CONTEXT("span x100", std::llround(t0 * 100.0f), std::llround(t1 * 100.0f));
  PassList solved;
  Temporal::arm_passes<W>(x0, x1, t0, t1, false, false, solved);
  const auto truth = brute_force_passes<W>(
      [&](double t) {
        return x0 + (static_cast<double>(x1) - x0) * (t - t0) / (t1 - t0);
      },
      t0, t1);
  HS_EXPECT_EQ(static_cast<int>(solved.passes.size()), expected);
  HS_EXPECT_EQ(static_cast<int>(truth.size()), expected);
  for (const TruePass &want : truth) {
    bool found = false;
    for (const Temporal::ArmPass &got : solved.passes)
      found =
          found || (got.arm == want.arm && fabs(got.time - want.time) < 2e-6 &&
                    fabs(t0 + (t1 - t0) * got.along - got.time) < 1e-6);
    HS_EXPECT_TRUE(found);
  }
}

/** @brief Linear tracks: zero, one and two passes, both arm cuts and partial
 *         spans match brute force. */
inline void test_arm_passes_match_brute_force() {
  check_linear_passes(40.0f, 42.0f, 0.0f, 1.0f, 1);
  check_linear_passes(200.0f, 190.0f, 0.0f, 1.0f, 1);
  // Moving with the sweep across the 143|144 cut: no arm passes it.
  check_linear_passes(143.5f, 145.5f, 0.0f, 1.0f, 0);
  // Against the sweep across 143|144: each arm passes it once.
  check_linear_passes(144.3f, 142.3f, 0.0f, 1.0f, 2);
  // The same across the 287|0 cut, unwrapped past the width.
  check_linear_passes(287.6f, 289.6f, 0.0f, 1.0f, 0);
  check_linear_passes(288.4f, 286.4f, 0.0f, 1.0f, 2);
  // Partial spans of x(t) = 100 + 4 t.
  check_linear_passes(101.0f, 103.0f, 0.25f, 0.75f, 1);
  check_linear_passes(100.0f, 102.0f, 0.0f, 0.5f, 0);
  check_linear_passes(102.0f, 104.0f, 0.5f, 1.0f, 1);
}

/** @brief Splat extensions emit a pass at the frame's start or end only when
 *         its column lies within a column of the frame. */
inline void test_arm_passes_extend_to_frame_ends() {
  constexpr int W = 288;
  PassList without, with;
  Temporal::arm_passes<W>(-0.5f, -0.3f, 0.0f, 1.0f, false, true, without);
  Temporal::arm_passes<W>(-0.5f, -0.3f, 0.0f, 1.0f, true, true, with);
  HS_EXPECT_EQ(static_cast<int>(without.passes.size()), 1);
  HS_EXPECT_EQ(static_cast<int>(with.passes.size()), 2);
  bool extended = false;
  for (const auto &pass : with.passes)
    extended = extended || (pass.arm == 0 && pass.time == 0.0f);
  HS_EXPECT_TRUE(extended);
  PassList far;
  Temporal::arm_passes<W>(-1.5f, -1.3f, 0.0f, 1.0f, true, true, far);
  for (const auto &pass : far.passes)
    HS_EXPECT_NE(pass.arm, 0);
}

/** @brief Piecewise tracks emit each pass once, inner knots included. */
inline void test_piecewise_arm_passes_match_brute_force() {
  constexpr int W = 288;
  const float knots[] = {10.0f, 12.0f, 9.0f};
  PassList solved;
  Temporal::piecewise_arm_passes<W>(knots, 3, false, false, solved);
  const auto truth = brute_force_passes<W>(
      [&](double t) {
        return t < 0.5 ? 10.0 + 4.0 * t : 12.0 - 6.0 * (t - 0.5);
      },
      0.0, 1.0);
  HS_EXPECT_EQ(solved.passes.size(), truth.size());
  for (size_t i = 0; i < truth.size() && i < solved.passes.size(); ++i)
    HS_EXPECT_NEAR(solved.passes[i].time, static_cast<float>(truth[i].time),
                   2e-6f);
  // A pass exactly at the inner knot: x = 72 at t = 1/2 meets arm 0.
  const float through[] = {70.0f, 72.0f, 74.0f};
  PassList knot;
  Temporal::piecewise_arm_passes<W>(through, 3, false, false, knot);
  HS_EXPECT_EQ(static_cast<int>(knot.passes.size()), 1);
}

/**
 * @brief Checks a modelled sphere track's passes against brute force.
 * @param exact Unit path over the span.
 * @param t0,t1 Span.
 * @param deviation Bound on the path's distance from its chord.
 * @param track Receives the model.
 * @return Whether the model holds.
 */
inline bool check_track_passes(const std::function<math::Vector(double)> &exact,
                               double t0, double t1, float deviation,
                               Temporal::LongitudeTrack &track) {
  constexpr int W = 288;
  constexpr double HALF = W / 2.0;
  const math::Vector from = exact(t0), to = exact(t1);
  if (!Temporal::model_track<W>(from, to, deviation,
                                static_cast<float>(t1 - t0), track))
    return false;
  // Unwrap each sampled column against the one before, from the origin.
  std::vector<double> samples;
  double last = track.origin;
  for (int i = 0; i <= 8192; ++i) {
    const math::Vector p = exact(t0 + (t1 - t0) * i / 8192.0);
    double x = std::atan2(static_cast<double>(p.z), static_cast<double>(p.x)) *
               (W / (2.0 * 3.14159265358979323846));
    x += W * std::round((last - x) / W);
    samples.push_back(x);
    last = x;
  }
  const auto track_column = [&](double t) {
    const double s = (t - t0) / (t1 - t0) * 8192.0;
    const int i = std::min(static_cast<int>(s), 8191);
    return samples[i] + (samples[i + 1] - samples[i]) * (s - i);
  };
  const auto truth = brute_force_passes<W>(track_column, t0, t1);
  PassList solved;
  Temporal::track_passes<W>(track, static_cast<float>(t0),
                            static_cast<float>(t1), true, true, solved);
  const double window = track.guard / HALF;
  for (const TruePass &want : truth) {
    bool found = false;
    for (const auto &got : solved.passes)
      found = found ||
              (got.arm == want.arm && fabs(got.time - want.time) <= window);
    HS_CONTEXT("true pass arm, time x1e6", want.arm,
               std::llround(want.time * 1e6));
    HS_EXPECT_TRUE(found);
  }
  for (const auto &got : solved.passes) {
    bool found = false;
    for (const TruePass &want : truth)
      found = found ||
              (got.arm == want.arm && fabs(got.time - want.time) <= window);
    // A pass with no true counterpart is a splat extension at a span end.
    found = found || got.time - t0 <= window || t1 - got.time <= window;
    HS_CONTEXT("solved pass arm, time x1e6", got.arm,
               std::llround(got.time * 1e6));
    HS_EXPECT_TRUE(found);
  }
  return true;
}

/** @brief The linear longitude model's guard bounds every true pass of
 *         rotating points, across arm cuts and partial spans. */
inline void test_model_track_guard_bounds_true_passes() {
  int modelled = 0;
  for (const float start_column : {0.3f, 40.0f, 143.6f, 144.4f, 287.7f}) {
    for (const float height : {-0.6f, 0.0f, 0.4f, 0.85f}) {
      for (const float angle : {-0.06f, -0.01f, 0.004f, 0.03f, 0.08f}) {
        for (const float tilt : {0.0f, 0.7f}) {
          HS_CONTEXT("case", modelled);
          const float theta = start_column * math::TWO_PI_F / 288.0f;
          const float r = sqrtf(1.0f - height * height);
          const math::Vector p0(r * cosf(theta), height, r * sinf(theta));
          const math::Vector axis =
              math::Vector(sinf(tilt), cosf(tilt), 0.0f).normalized();
          const std::function<math::Vector(double)> exact = [&](double t) {
            return math::rotate(
                p0, math::make_rotation(axis, angle * static_cast<float>(t)));
          };
          for (const auto &span :
               {std::pair{0.0, 1.0}, std::pair{0.25, 0.75}}) {
            const float deviation =
                Temporal::PathBound::rotation(
                    angle * static_cast<float>(span.second - span.first))
                    .chord_error();
            Temporal::LongitudeTrack track;
            HS_EXPECT_TRUE(check_track_passes(exact, span.first, span.second,
                                              deviation, track));
            ++modelled;
          }
        }
      }
    }
  }
  HS_EXPECT_EQ(modelled, 400);
}

/** @brief Near a pole the model refuses a track it cannot bound, and a fast
 *         high-latitude track it refuses whole holds in two closed-form
 *         pieces. */
inline void test_model_track_splits_or_refuses_near_poles() {
  constexpr int W = 288;
  // A point carried over the pole: its longitude swings through half a turn.
  const math::Vector near_pole =
      math::Vector(0.02f, 0.999f, 0.03f).normalized();
  const float over = 0.08f;
  const math::Vector over_end =
      math::rotate(near_pole, math::make_rotation(math::X_AXIS, over));
  Temporal::LongitudeTrack track;
  HS_EXPECT_FALSE(Temporal::model_track<W>(
      near_pole, over_end, Temporal::PathBound::rotation(over).chord_error(),
      1.0f, track));
  // A fast spin at high latitude: refused whole, modelled in halves.
  const float height = 0.97f, angle = 1.3f;
  const float r = sqrtf(1.0f - height * height);
  const math::Vector p0(r, height, 0.0f);
  const std::function<math::Vector(double)> spin = [&](double t) {
    return math::rotate(
        p0, math::make_rotation(math::Y_AXIS, angle * static_cast<float>(t)));
  };
  HS_EXPECT_FALSE(Temporal::model_track<W>(
      spin(0.0), spin(1.0),
      Temporal::PathBound::rotation(angle, r).chord_error(), 1.0f, track));
  const float half_deviation =
      Temporal::PathBound::rotation(0.5f * angle, r).chord_error();
  HS_EXPECT_TRUE(check_track_passes(spin, 0.0, 0.5, half_deviation, track));
  HS_EXPECT_TRUE(check_track_passes(spin, 0.5, 1.0, half_deviation, track));
}
