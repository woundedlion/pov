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
 * @brief A point rotated by each frame of an orientation history.
 * @param history Orientation history.
 * @param v Unrotated point.
 */
inline std::vector<math::Vector>
history_nodes(const math::Orientation<> &history, const math::Vector &v) {
  std::vector<math::Vector> nodes;
  for (int i = 0; i < history.length(); ++i)
    nodes.push_back(math::rotate(v, history.get(i)));
  return nodes;
}

/**
 * @brief Captures sorted by fraction.
 * @param got Captures to sort.
 */
inline Temporal::Captures sorted(Temporal::Captures got) {
  if (got.count == 2 && got.fraction[1] < got.fraction[0])
    std::swap(got.fraction[0], got.fraction[1]);
  return got;
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
  for (float drift :
       {-0.6f, -0.3f, -0.05f, -0.004f, 0.004f, 0.05f, 0.3f, 0.6f}) {
    for (int k = 0; k < 64; ++k) {
      const float start = 2.0f * math::PI_F * k / 64 + 0.013f;
      const math::Orientation<> history = spin_history(start, drift, 3);
      const math::Vector v(0.8f, 0.6f, 0.0f);
      auto position = [&](float u) { return math::rotate(v, history.get(u)); };
      const Temporal::Captures got =
          sorted(Temporal::capture(history_nodes(history, v)));
      const std::vector<float> want = oracle_captures(position);
      ++cases;
      zero += want.empty();
      two += want.size() == 2;
      HS_EXPECT_EQ(got.count, static_cast<int>(want.size()));
      for (int i = 0; i < std::min<int>(got.count, want.size()); ++i)
        HS_EXPECT_NEAR(got.fraction[i], want[i], 0.5f / 144);
    }
  }
  HS_EXPECT_GT(zero, 0);
  HS_EXPECT_GT(two, 0);
  HS_EXPECT_GT(cases - zero - two, 0);
}

/**
 * @brief Under tilted-axis rotation, linear azimuth between history frames
 *        places each pass within half a column of the exhaustive oracle.
 */
inline void test_temporal_capture_tilted_axis_within_half_column() {
  constexpr int W = 288;
  hs::random().seed(0x7E3B);
  int checked = 0, boundary = 0;
  float worst = 0.0f;
  for (int trial = 0; trial < 400; ++trial) {
    const math::Vector axis = rand_unit();
    const float frame_angle = hs::rand_f(-12.0f, 12.0f) * 2.0f * math::PI_F / W;
    const math::Quaternion start =
        math::make_rotation(rand_unit(), hs::rand_f(0.0f, 6.0f));
    math::Orientation<> history(start);
    for (int i = 1; i <= 3; ++i)
      history.push(math::make_rotation(axis, frame_angle * i / 3) * start);
    const math::Vector v = rand_unit();
    if (std::fabs(math::rotate(v, start).y) > 0.9f)
      continue;
    auto position = [&](float u) { return math::rotate(v, history.get(u)); };
    const Temporal::Captures got =
        sorted(Temporal::capture(history_nodes(history, v)));
    const std::vector<float> want = oracle_captures(position);
    if (got.count != static_cast<int>(want.size())) {
      // A pass this close to a frame boundary may land on either side of it.
      float edge = 1.0f;
      for (float u : want)
        edge = std::min({edge, u, 1.0f - u});
      for (int i = 0; i < got.count; ++i)
        edge = std::min({edge, got.fraction[i], 1.0f - got.fraction[i]});
      HS_EXPECT_LT(edge * (W / 2), 0.5f);
      ++boundary;
      continue;
    }
    for (int i = 0; i < got.count; ++i) {
      worst = std::max(worst, std::fabs(got.fraction[i] - want[i]) * (W / 2));
      ++checked;
    }
  }
  HS_EXPECT_GT(checked, 200);
  HS_EXPECT_LT(boundary, 10);
  HS_EXPECT_LT(worst, 0.5f);
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
      const Temporal::Captures got =
          Temporal::capture(history_nodes(history, v));
      for (int i = 0; i < got.count; ++i) {
        const int x = math::fast_wrap(
            static_cast<int>(std::round(
                math::vector_to_theta<W>(position(got.fraction[i])))),
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
 * @brief Point on the equator at an azimuth given in half-turns.
 * @param s Azimuth in half-turns.
 */
inline math::Vector at_half_turns(float s) {
  return {std::cos(s * math::PI_F), 0.0f, std::sin(s * math::PI_F)};
}

/**
 * @brief A rise in arm-relative position between history frames, or a fall of
 *        two or more, is emitted once at the current endpoint.
 */
inline void test_temporal_capture_unsolvable_motion_falls_back_to_endpoint() {
  const std::array<math::Vector, 3> rising{
      at_half_turns(0.1f), at_half_turns(0.2f), at_half_turns(0.9f)};
  const std::array<math::Vector, 5> falling{
      at_half_turns(0.1f), at_half_turns(-0.4f), at_half_turns(-0.9f),
      at_half_turns(0.6f), at_half_turns(0.1f)};
  const std::array<math::Vector, 3> solvable{
      at_half_turns(0.4f), at_half_turns(0.35f), at_half_turns(0.3f)};
  for (const Temporal::Captures got :
       {Temporal::capture(rising), Temporal::capture(falling)}) {
    HS_EXPECT_EQ(got.count, 1);
    HS_EXPECT_NEAR(got.fraction[0], 1.0f, 0.0f);
  }
  const Temporal::Captures got = Temporal::capture(solvable);
  HS_EXPECT_EQ(got.count, 1);
  HS_EXPECT_NEAR(got.fraction[0], 0.4f / 0.55f * 0.5f, 2e-3f);
}

/**
 * @brief Orient emits each arm pass at the history pose of its fraction, with
 *        age offset by the remaining frame, and only the current pose with
 *        temporal sampling disabled.
 */
inline void test_world_orient_emits_at_capture_time() {
  const math::Orientation<> history = spin_history(-0.5f, -0.08f, 3);
  const math::Vector v(0.0f, 0.0f, 1.0f);
  const Temporal::Captures want = Temporal::capture(history_nodes(history, v));
  HS_EXPECT_EQ(want.count, 1);
  Filter::World::Orient orient(history);
  std::vector<Tap3D> taps;
  auto record = [&](const math::Vector &p, const Pixel &c, float age,
                    float alpha) { taps.push_back({p, c, age, alpha}); };
  orient.plot(v, Pixel(1, 2, 3), 4.0f, 0.5f, record);
  HS_EXPECT_EQ(taps.size(), 1u);
  HS_EXPECT_LT(want.fraction[0], 0.99f);
  const math::Vector at = math::rotate(v, history.get(want.fraction[0]));
  HS_EXPECT_NEAR(taps[0].age, 4.0f + 1.0f - want.fraction[0], 1e-5f);
  HS_EXPECT_NEAR(taps[0].v.x, at.x, 1e-6f);
  HS_EXPECT_NEAR(taps[0].v.z, at.z, 1e-6f);
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

/**
 * @brief Sweep slabs tile the half-turn in sweep order, and a slab of an even
 *        division holds exactly the columns whose sweep fractions fall in its
 *        share of the frame.
 */
inline void test_sweep_slabs_tile_the_half_turn() {
  for (int width : {20, 96, 288}) {
    const int half = width / 2;
    for (int count : {1, 2, 3, 5, 6, 7, half}) {
      HS_CONTEXT("width, slabs", width, count);
      int next = 0;
      for (int slab = 0; slab < count; ++slab) {
        const Temporal::SweepSlab columns =
            Temporal::sweep_slab(slab, count, width);
        HS_EXPECT_EQ(columns.begin, next);
        HS_EXPECT_LT(columns.begin, columns.end);
        next = columns.end;
        if (half % count != 0)
          continue;
        for (int x = 0; x < width; ++x) {
          const float u = Temporal::sweep_fraction(x, width);
          const bool inside =
              x % half >= columns.begin && x % half < columns.end;
          HS_EXPECT_EQ(inside, u >= static_cast<float>(slab) / count &&
                                   u < static_cast<float>(slab + 1) / count);
        }
      }
      HS_EXPECT_EQ(next, half);
    }
  }
}

/** @brief Explicit fractions follow the pose lattice and refresh next frame. */
inline void test_world_orient_explicit_time() {
  math::Orientation<> history = spin_history(-0.5f, -0.08f, 3);
  Filter::World::Orient orient(history);
  const math::Vector v(0.6f, -0.8f, 0.0f);
  StubEffect effect(288, 144);
  Pipeline<288, 144, Filter::World::Orient, Filter::Screen::AntiAlias<288, 144>>
      timed{Filter::World::Orient(history),
            Filter::Screen::AntiAlias<288, 144>()};
  Pipeline<288, 144, Filter::Screen::AntiAlias<288, 144>> plain;
  for (int frame = 0; frame < 2; ++frame) {
    {
      Canvas canvas(effect);
      HS_EXPECT_EQ(canvas.frame_number(), static_cast<uint32_t>(frame + 1));
      for (float t : {0.0f, 0.13f, 0.5f, 0.73f, 1.0f}) {
        const math::Vector want = math::rotate(v, history.get(t));
        int calls = 0;
        orient.plot_at(
            canvas, v, Pixel(1, 2, 3), 7.0f, 0.25f, t,
            [&](const math::Vector &p, const Pixel &c, float age, float alpha) {
              ++calls;
              HS_EXPECT_NEAR((p - want).magnitude(), 0.0f, 2e-5f);
              HS_EXPECT_EQ(c, Pixel(1, 2, 3));
              HS_EXPECT_EQ(age, 7.0f);
              HS_EXPECT_EQ(alpha, 0.25f);
            });
        HS_EXPECT_EQ(calls, 1);
        std::fill_n(canvas.data(), 288 * 144, Pixel{});
        plain.plot(canvas, orient.at_time(canvas, v, t), Pixel(65535, 0, 0),
                   7.0f, 1.0f);
        std::vector<Pixel> expected(canvas.data(), canvas.data() + 288 * 144);
        std::fill_n(canvas.data(), 288 * 144, Pixel{});
        timed.plot_at(canvas, v, Pixel(65535, 0, 0), 7.0f, 1.0f, t);
        for (size_t i = 0; i < expected.size(); ++i)
          HS_EXPECT_EQ(canvas.data()[i], expected[i]);
      }
      const math::Vector end = math::rotate(v, history.get());
      const math::Vector at_end = orient.at_time(canvas, v, 1.0f);
      HS_EXPECT_EQ(at_end.x, end.x);
      HS_EXPECT_EQ(at_end.y, end.y);
      HS_EXPECT_EQ(at_end.z, end.z);
      Temporal::set_temporal_sampling_enabled(false);
      const math::Vector disabled = orient.at_time(canvas, v, 0.13f);
      Temporal::set_temporal_sampling_enabled(true);
      HS_EXPECT_EQ(disabled.x, end.x);
      HS_EXPECT_EQ(disabled.y, end.y);
      HS_EXPECT_EQ(disabled.z, end.z);
    }
    effect.advance_display();
    history = spin_history(0.7f, 0.12f, 3);
  }
}

/** @brief Callback capture agrees with exhaustive search on each arm branch. */
inline void test_geometry_capture_branches() {
  int zero = 0, one = 0, two = 0;
  for (float delta : {-30.0f, -12.0f, -3.0f, 0.0f, 3.0f, 12.0f, 30.0f}) {
    for (int i = 0; i < 96; ++i) {
      const float start = i * 3.0f + 0.23f;
      for (bool curved : {false, true}) {
        auto position = [&](float t) {
          return at_half_turns((start + delta * (curved ? t * t : t)) / 144.0f);
        };
        auto probe = [&](float t) {
          return math::vector_to_pixel<288, 144>(position(t));
        };
        const auto got = sorted(Temporal::solve_capture<288>(
            probe, probe(1.0f), fmaxf(4.1f, fabsf(delta)), false));
        const auto want = oracle_captures(position);
        zero += want.empty();
        one += want.size() == 1;
        two += want.size() == 2;
        HS_EXPECT_EQ(got.count, static_cast<int>(want.size()));
        for (int pass = 0; pass < std::min<int>(got.count, want.size());
             ++pass) {
          const auto endpoint = probe(1.0f);
          const float unwrapped =
              endpoint.x +
              std::remainder(probe(want[pass]).x - endpoint.x, 288.0f);
          const int branch =
              static_cast<int>(lroundf(unwrapped / 144.0f - want[pass]));
          const float branch_time = Temporal::solve_capture_branch<288>(
              probe, endpoint, fmaxf(4.1f, fabsf(delta)), false, branch);
          HS_EXPECT_NEAR(branch_time, got.fraction[pass], 1e-6f);
          const float error = std::remainder(probe(got.fraction[pass]).x -
                                                 got.fraction[pass] * 144.0f,
                                             144.0f);
          HS_EXPECT_LE(fabsf(error), 0.501f);
          HS_EXPECT_LE(
              fabsf(std::remainder(
                  probe(got.fraction[pass]).x - probe(want[pass]).x, 288.0f)),
              0.5f);
        }
      }
    }
  }
  HS_EXPECT_GT(zero, 0);
  HS_EXPECT_GT(one, 0);
  HS_EXPECT_GT(two, 0);
}

/** @brief Shortcuts limit probes and polar anchors are held and counted. */
inline void test_geometry_capture_shortcuts_and_poles() {
  int probes = 0;
  auto probe = [&](float t) {
    ++probes;
    return math::PixelCoords{143.0f + 3.0f * t, 72.0f};
  };
  const math::PixelCoords branch_end{146.0f, 72.0f};
  HS_EXPECT_EQ(
      Temporal::solve_capture_branch<288>(probe, branch_end, 3.0f, true, 1),
      0.0f);
  HS_EXPECT_EQ(
      Temporal::solve_capture_branch<288>(probe, branch_end, 3.0f, true, 0),
      1.0f);
  probes = 0;
  Temporal::geometry_counters = {};
  const auto no_pass = Temporal::solve_capture<288>(
      probe, math::PixelCoords{146.0f, 72.0f}, 3.0f, true);
  HS_EXPECT_EQ(no_pass.count, 0);
  HS_EXPECT_EQ(probes, 1);
  probes = 0;
  const auto shortcut = Temporal::solve_capture<288>(
      probe, math::PixelCoords{80.0f, 72.0f}, 3.0f, false);
  HS_EXPECT_EQ(shortcut.count, 1);
  HS_EXPECT_EQ(probes, 0);
  HS_EXPECT_EQ(shortcut.fraction[0], 80.0f / 144.0f);
  auto reverse = [&](float t) {
    ++probes;
    return math::PixelCoords{145.0f - 3.0f * t, 72.0f};
  };
  const auto two = sorted(Temporal::solve_capture<288>(
      reverse, math::PixelCoords{142.0f, 72.0f}, 3.0f, true));
  HS_EXPECT_EQ(two.count, 2);
  HS_EXPECT_EQ(probes, 1);
  HS_EXPECT_NEAR(two.fraction[0], 1.0f / 147.0f, 1e-6f);
  HS_EXPECT_NEAR(two.fraction[1], 145.0f / 147.0f, 1e-6f);
  probes = 0;
  const Temporal::MotionBound bound{0.1f, false};
  for (float row : {0.0f, 0.25f, 142.75f, 143.0f}) {
    const float columns = bound.columns<288, 144>(row);
    HS_EXPECT_GE(columns, Temporal::MAX_CAPTURE_COLUMNS<288>);
    const auto held = Temporal::solve_capture<288>(
        probe, math::PixelCoords{143.5f, row}, columns, false);
    HS_EXPECT_EQ(held.count, 1);
    HS_EXPECT_EQ(held.fraction[0], 143.5f / 144.0f);
  }
  HS_EXPECT_EQ(probes, 0);
  HS_EXPECT_EQ(Temporal::geometry_counters.held_anchors, 4u);
  HS_EXPECT_LT((bound.columns<288, 144>(71.5f)), 5.0f);
}

/** @brief Probes transform at their fraction without recording trail samples. */
inline void test_timed_probe_bypasses_history() {
  const math::Orientation<> history = spin_history(-0.5f, -0.08f, 3);
  using Trails = Filter::World::Trails<16>;
  Pipeline<288, 144, Trails, Filter::World::Orient,
           Filter::Screen::AntiAlias<288, 144>>
      pipeline{Trails(10), Filter::World::Orient(history),
               Filter::Screen::AntiAlias<288, 144>()};
  pipeline.init_storage(persistent_arena);
  PipelineRef ref(pipeline, Temporal::Timed{});
  StubEffect effect(288, 144);
  Canvas canvas(effect);
  const math::Vector v(0.6f, -0.8f, 0.0f);
  for (float t : {0.0f, 0.3f, 0.6f, 1.0f}) {
    const auto got = ref.probe(canvas, v, t);
    const auto want =
        math::vector_to_pixel<288, 144>(math::rotate(v, history.get(t)));
    HS_EXPECT_NEAR(got.x, want.x, 0.002f);
    HS_EXPECT_NEAR(got.y, want.y, 0.002f);
  }
  int records = 0;
  auto shader = [&](const math::Vector &, float) {
    ++records;
    return Color4(1, 1, 1, 1);
  };
  pipeline.flush(canvas, WorldTrailFn(shader), 1.0f);
  HS_EXPECT_EQ(records, 0);
  for (int i = 0; i < 288 * 144; ++i)
    HS_EXPECT_EQ(canvas.data()[i], Pixel{});
  ref.plot_at(canvas, v, Pixel(65535, 0, 0), 0.0f, 1.0f, 0.3f);
  pipeline.flush(canvas, WorldTrailFn(shader), 1.0f);
  HS_EXPECT_EQ(records, 1);
}

/** @brief Predicted landings reject writes and history before one corrected plot. */
inline void test_timed_predicted_landing() {
  StubEffect effect(288, 144);
  Canvas canvas(effect);
  const Pixel color(65535, 0, 0);
  const auto point = math::pixel_to_vector<288, 144>(287.5f, 70.0f);
  Filter::Screen::DirectAntiAliasSink<288, 144> direct;
  direct.prepare(canvas);
  Pipeline<288, 144, Filter::Screen::AntiAlias<288, 144>> screen;
  PipelineRef direct_ref(direct, canvas, Temporal::Timed{});
  PipelineRef screen_ref(screen, Temporal::Timed{});
  for (PipelineRef *ref : {&direct_ref, &screen_ref}) {
    std::fill_n(canvas.data(), 288 * 144, Pixel{});
    const Temporal::ColumnPrediction wrong{2.0f, 0.1f};
    const float landed =
        ref->plot_at(canvas, point, color, 0.0f, 0.4f, 1.0f, &wrong);
    HS_EXPECT_NEAR(landed, ref->probe(canvas, point, 1.0f).x - 288.0f, 0.0f);
    HS_EXPECT_TRUE(wrong.needs_correction(landed));
    int lit = 0;
    for (int i = 0; i < 288 * 144; ++i)
      lit += canvas.data()[i].r != 0;
    HS_EXPECT_EQ(lit, 0);
    const Temporal::ColumnPrediction outside{landed, 0.1f, 0.0f, 144.0f};
    ref->plot_at(canvas, point, color, 0.0f, 0.4f, 1.0f, &outside);
    HS_EXPECT_EQ(canvas(287, 70), Pixel{});
    const Temporal::ColumnPrediction corrected{landed, 0.1f};
    HS_EXPECT_NEAR(
        ref->plot_at(canvas, point, color, 0.0f, 0.4f, 1.0f, &corrected),
        landed, 0.0f);
    HS_EXPECT_GT(canvas(287, 70).r, 0);
  }
  using Trails = Filter::World::Trails<16>;
  const auto orientation = spin_history(0.0f, 0.02f, 3);
  Pipeline<288, 144, Trails, Filter::World::Orient,
           Filter::Screen::AntiAlias<288, 144>>
      history{Trails(10), Filter::World::Orient(orientation),
              Filter::Screen::AntiAlias<288, 144>()};
  history.init_storage(persistent_arena);
  PipelineRef ref(history, Temporal::Timed{});
  const auto projected = ref.probe(canvas, point, 0.3f);
  const Temporal::ColumnPrediction wrong{projected.x + 4.0f, 0.1f};
  HS_EXPECT_NEAR(ref.plot_at(canvas, point, color, 0.0f, 0.4f, 0.3f, &wrong),
                 projected.x, 0.0f);
  int records = 0;
  auto shader = [&](const math::Vector &, float) {
    ++records;
    return Color4(1, 1, 1, 1);
  };
  history.flush(canvas, WorldTrailFn(shader), 1.0f);
  HS_EXPECT_EQ(records, 0);
  const Temporal::ColumnPrediction corrected{projected.x, 0.1f};
  ref.plot_at(canvas, point, color, 0.0f, 0.4f, 0.3f, &corrected);
  history.flush(canvas, WorldTrailFn(shader), 1.0f);
  HS_EXPECT_EQ(records, 1);
}

/** @brief Unborn predicted samples never enter world history. */
inline void test_timed_unborn_samples_skip_history() {
  constexpr int W = 288, H = 144;
  StubEffect effect(W, H);
  Canvas canvas(effect);
  using Trails = Filter::World::Trails<16>;
  Pipeline<W, H, Trails, Filter::Screen::AntiAlias<W, H>> pipeline{
      Trails(10), Filter::Screen::AntiAlias<W, H>()};
  pipeline.init_storage(persistent_arena);
  ScratchScope scratch(scratch_arena_a);
  Fragments points;
  points.bind(scratch_arena_a, 2);
  auto geometry = [](size_t index, float time) {
    return math::pixel_to_vector<W, H>(40.0f + 30.0f * index + 3.0f * time,
                                       65.0f);
  };
  for (size_t index = 0; index < 2; ++index) {
    Fragment point;
    point.pos = geometry(index, 1.0f);
    points.push_back(point);
  }
  const float births[] = {1.0f, 1.0f};
  auto callback = Plot::GeometryCallback::exact(geometry);
  callback.births = births;
  auto shade = [](const math::Vector &, Fragment &fragment) {
    fragment.color = Color4(1, 1, 1, 0.4f);
  };
  Plot::rasterize<W, H, {.timed = true}>(pipeline, canvas, points, shade, {},
                                         callback,
                                         3.0f * math::RADIANS_PER_COLUMN<W>);
  int records = 0;
  auto history_shader = [&](const math::Vector &, float) {
    ++records;
    return Color4(1, 1, 1, 1);
  };
  pipeline.flush(canvas, WorldTrailFn(history_shader), 1.0f);
  HS_EXPECT_EQ(records, 0);
}

/** @brief Both timed entries emit no-pass points zero times and two-pass twice. */
inline void test_timed_point_entries() {
  for (bool direct : {false, true}) {
    StubEffect effect(288, 144);
    Canvas canvas(effect);
    Filter::Screen::DirectAntiAliasSink<288, 144> sink;
    sink.prepare(canvas);
    Pipeline<288, 144, Filter::Screen::AntiAlias<288, 144>> pipeline;
    PipelineRef ref(pipeline, Temporal::Timed{});
    for (float delta : {3.0f, -3.0f}) {
      auto geometry = [&](float t) {
        return math::pixel_to_vector<288, 144>(
            (delta > 0.0f ? 143.0f : 145.0f) + delta * t, 72.0f);
      };
      const Temporal::MotionBound bound{3.0f * 2.0f * math::PI_F / 288.0f,
                                        true};
      Temporal::geometry_counters = {};
      std::fill_n(canvas.data(), 288 * 144, Pixel{});
      if (direct)
        sink.plot_timed(canvas, geometry, bound, Pixel(65535, 0, 0), 0.0f, 1.0f,
                        true);
      else
        ref.plot_timed<288, 144>(canvas, geometry, bound, Pixel(65535, 0, 0),
                                 0.0f, 1.0f, true);
      int lit = 0;
      for (int i = 0; i < 288 * 144; ++i)
        lit += canvas.data()[i].r != 0;
      if (delta > 0.0f)
        HS_EXPECT_EQ(lit, 0);
      else {
        HS_EXPECT_GT(canvas(145, 72).r, 0);
        HS_EXPECT_GT(canvas(142, 72).r, 0);
      }
      HS_EXPECT_EQ(Temporal::geometry_counters.end_probes, 1u);
      HS_EXPECT_EQ(Temporal::geometry_counters.solve_probes, 1u);
      HS_EXPECT_EQ(Temporal::geometry_counters.callback_evaluations,
                   delta > 0.0f ? 2u : 4u);
    }
  }
}

/** @brief Complementary boundary parts preserve the full still splat weights. */
inline void test_timed_boundary_tap_ownership() {
  StubEffect effect(288, 144);
  Canvas canvas(effect);
  Pipeline<288, 144, Filter::Screen::AntiAlias<288, 144>> pipeline;
  Filter::Screen::DirectAntiAliasSink<288, 144> sink;
  sink.prepare(canvas);
  const Pixel color(65535, 12000, 30000);
  for (float column : {143.5f, 287.5f}) {
    const auto v = math::pixel_to_vector<288, 144>(column, 72.5f);
    std::fill_n(canvas.data(), 288 * 144, Pixel{});
    pipeline.plot(canvas, v, color, 0.0f, 0.4f);
    std::vector<Pixel> expected(canvas.data(), canvas.data() + 288 * 144);
    for (bool direct : {false, true}) {
      std::fill_n(canvas.data(), 288 * 144, Pixel{});
      for (float t : {0.0f, 1.0f}) {
        if (direct)
          sink.plot_at(canvas, v, color, 0.0f, 0.4f, t);
        else
          pipeline.plot_at(canvas, v, color, 0.0f, 0.4f, t);
      }
      for (size_t i = 0; i < expected.size(); ++i)
        HS_EXPECT_EQ(canvas.data()[i], expected[i]);
    }
  }
}

/** @brief Motion padding admits an entering edge while preserving the clip. */
inline void test_timed_cull_region() {
  const auto history = spin_history(0.0f, 0.10f, 3);
  Pipeline<288, 144, Filter::World::Orient, Filter::Screen::AntiAlias<288, 144>>
      pipeline{Filter::World::Orient(history),
               Filter::Screen::AntiAlias<288, 144>()};
  PipelineRef ref(pipeline, Temporal::Timed{});
  constexpr float ANGLE = 0.10f * math::PI_F;
  HS_EXPECT_GE(ref.cull_bound(), ANGLE);
  HS_EXPECT_LE(ref.cull_bound(), ANGLE + (history.length() - 1) * 0.000203f);
  StubEffect effect(288, 144);
  effect.set_margin(0);
  effect.set_clip(36, 72, 0, 144);
  Canvas canvas(effect);
  const auto exact = canvas.clip();
  const Plot::CullRegion still(canvas.clip());
  const Plot::CullRegion swept(canvas.clip(), ref.cull_bound());
  const auto a = math::pixel_to_vector<288, 144>(5.0f, 52.0f);
  const auto b = math::pixel_to_vector<288, 144>(6.0f, 54.0f);
  HS_EXPECT_FALSE((Plot::edge_visible_in_region<288, 144>(ref, still, a, b)));
  HS_EXPECT_TRUE((Plot::edge_visible_in_region<288, 144>(ref, swept, a, b)));
  const auto below_a = math::pixel_to_vector<288, 144>(60.0f, 75.0f);
  const auto below_b = math::pixel_to_vector<288, 144>(62.0f, 76.0f);
  HS_EXPECT_FALSE(
      (Plot::edge_visible_in_region<288, 144>(ref, still, below_a, below_b)));
  HS_EXPECT_TRUE(
      (Plot::edge_visible_in_region<288, 144>(ref, swept, below_a, below_b)));
  HS_EXPECT_TRUE(canvas.clip() == exact);
  HS_EXPECT_FALSE(still.could_intersect_rows<144>(74.0f, 75.0f));
  HS_EXPECT_TRUE(swept.could_intersect_rows<144>(74.0f, 75.0f));
  HS_EXPECT_LT(swept.column_pad<288>(0.05f), 0.0f);
  const float expected = ref.cull_bound() * (288.0f / (2.0f * math::PI_F)) /
                         (1.0f - ref.cull_bound());
  HS_EXPECT_NEAR(swept.column_pad<288>(1.0f), expected, 1e-6f);
}

/** @brief Cached explicit poses invalidate on fraction, input, frame and collapse. */
inline void test_world_orient_sample_cache() {
  math::Orientation<> history = spin_history(-0.5f, -0.08f, 3);
  Filter::World::Orient cached(history);
  StubEffect effect(288, 144);
  const math::Vector a(0.6f, -0.8f, 0.0f);
  const math::Vector b(-0.2f, 0.3f, 0.7f);
  for (int frame = 0; frame < 2; ++frame) {
    {
      Canvas canvas(effect);
      auto compare = [&](const math::Vector &v, float t) {
        Filter::World::Orient fresh(history);
        const auto expected = fresh.at_time(canvas, v, t);
        const auto actual = cached.at_time(canvas, v, t);
        HS_EXPECT_EQ(actual.x, expected.x);
        HS_EXPECT_EQ(actual.y, expected.y);
        HS_EXPECT_EQ(actual.z, expected.z);
      };
      compare(a, 0.13f);
      compare(a, 0.13f);
      compare(b, 0.13f);
      compare(b, 0.73f);
      compare(a, 0.13f);
      compare(a, 1.0f);
      Temporal::set_temporal_sampling_enabled(false);
      compare(a, 0.13f);
      Temporal::set_temporal_sampling_enabled(true);
      compare(a, 0.13f);
      if (frame == 0) {
        history.collapse();
        compare(b, 0.13f);
        history = spin_history(0.4f, 0.2f, 3);
        compare(a, 0.13f);
      }
    }
    effect.advance_display();
    history = spin_history(0.7f, 0.12f, 3);
  }
}
/** @brief Column arithmetic matches centered IEEE remainders across arm branches. */
inline void test_temporal_column_arithmetic_matches_remainder() {
  for (float period : {12.0f, 288.0f, 1024.0f}) {
    for (int quarter = -64; quarter <= 64; ++quarter) {
      const float value = period * quarter * 0.25f;
      HS_EXPECT_EQ(Temporal::centered_remainder(value, period),
                   std::remainder(value, period));
    }
    for (int branch = -4; branch <= 4; ++branch)
      for (float offset : {-0.25f, 0.0f, 0.25f}) {
        const float value = (branch + 0.5f) * period + offset;
        HS_EXPECT_EQ(Temporal::centered_remainder(value, period),
                     std::remainder(value, period));
      }
  }
  for (int quarter = -2304; quarter <= 2304; ++quarter) {
    const float column = quarter * 0.25f;
    const float turns = column / 144.0f;
    const float expected = turns - floorf(turns);
    HS_EXPECT_NEAR(Temporal::column_fraction<288>(column), expected, 0.000001f);
    HS_EXPECT_GE(Temporal::column_fraction<288>(column), 0.0f);
    HS_EXPECT_LT(Temporal::column_fraction<288>(column), 1.0f);
  }
}

/** @brief Approximate orientation motion bounds remain conservative and preserve exact stillness. */
inline void test_temporal_orientation_bound_is_conservative() {
  const bool enabled = Temporal::temporal_sampling_enabled();
  Temporal::set_temporal_sampling_enabled(true);
  for (int index = 0; index <= 1024; ++index)
    for (float sign : {-1.0f, 1.0f}) {
      const float cosine = float(index) / 1024;
      math::Orientation<> history;
      history.set(math::Quaternion{});
      history.push(math::Quaternion(
          sign * cosine, sqrtf(fmaxf(0.0f, 1.0f - cosine * cosine)), 0, 0));
      Filter::World::Orient orient(history);
      const float expected =
          2.0f *
          acosf(fminf(fabsf(history.get(0).r * history.get(1).r +
                            math::dot(history.get(0).v, history.get(1).v)),
                      1.0f));
      const float actual = orient.cull_bound();
      HS_EXPECT_GE(actual, expected);
      HS_EXPECT_LE(actual - expected, 0.000203f);
      if (index == 1024)
        HS_EXPECT_NEAR(actual, 0.0f, 0.0f);
    }
  math::Orientation<> history;
  history.set(math::Quaternion{});
  history.push(math::make_rotation(math::Y_AXIS, .5f));
  Filter::World::Orient orient(history);
  Temporal::set_temporal_sampling_enabled(false);
  HS_EXPECT_NEAR(orient.cull_bound(), 0.0f, 0.0f);
  Temporal::set_temporal_sampling_enabled(enabled);
}
