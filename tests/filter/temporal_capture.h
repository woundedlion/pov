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
inline Spacetime::Captures sorted(Spacetime::Captures got) {
  if (got.count == 2 && got.fraction[1] < got.fraction[0])
    std::swap(got.fraction[0], got.fraction[1]);
  return got;
}

/**
 * @brief Arm passes of a moving point found by exhaustive search with exact
 *        `atan2`, independent of `Spacetime::capture`.
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
      const Spacetime::Captures got =
          sorted(Spacetime::capture(history_nodes(history, v)));
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
    const Spacetime::Captures got =
        sorted(Spacetime::capture(history_nodes(history, v)));
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
      const Spacetime::Captures got =
          Spacetime::capture(history_nodes(history, v));
      for (int i = 0; i < got.count; ++i) {
        const int x = math::fast_wrap(
            static_cast<int>(std::round(
                math::vector_to_theta<W>(position(got.fraction[i])))),
            W);
        const math::Vector field =
            position(Spacetime::SamplingSurface<W>::time(x));
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
  for (const Spacetime::Captures got :
       {Spacetime::capture(rising), Spacetime::capture(falling)}) {
    HS_EXPECT_EQ(got.count, 1);
    HS_EXPECT_NEAR(got.fraction[0], 1.0f, 0.0f);
  }
  const Spacetime::Captures got = Spacetime::capture(solvable);
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
  const Spacetime::Captures want =
      Spacetime::capture(history_nodes(history, v));
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

  Spacetime::set_flat_surface(true);
  taps.clear();
  orient.plot(v, Pixel(1, 2, 3), 4.0f, 0.5f, record);
  Spacetime::set_flat_surface(false);
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
    const double U = std::numeric_limits<float>::epsilon() / 2.0;
    const double GAMMA = 2 * U / (1 - 2 * U);
    for (int x = 0; x < width; ++x) {
      const double FRACTION = static_cast<double>(x % half) / half;
      const float actual = Spacetime::surface_time(x, width);
      HS_EXPECT_NEAR(actual, FRACTION, GAMMA * FRACTION);
      HS_EXPECT_EQ(Spacetime::surface_column(actual, width), x % half);
    }
    for (int count : {1, 2, 3, 5, 6, 7, half}) {
      HS_CONTEXT("width, slabs", width, count);
      int next = 0;
      for (int slab = 0; slab < count; ++slab) {
        const Spacetime::TimeSlab columns =
            Spacetime::time_slab(slab, count, width);
        HS_EXPECT_EQ(columns.begin, next);
        HS_EXPECT_LT(columns.begin, columns.end);
        next = columns.end;
        if (half % count != 0)
          continue;
        for (int x = 0; x < width; ++x) {
          const bool inside =
              x % half >= columns.begin && x % half < columns.end;
          const int sweep = (x % half) * count;
          HS_EXPECT_EQ(inside,
                       sweep >= slab * half && sweep < (slab + 1) * half);
        }
      }
      HS_EXPECT_EQ(next, half);
    }
  }
}

/**
 * @brief Pushes three sub-steps about changing axes onto a collapsed
 *        orientation, so its frame path leaves its endpoint slerp.
 * @param history Orientation to advance.
 * @param frame Frame number, varying the axes.
 */
inline void push_turning_steps(math::Orientation<> &history, int frame) {
  history.collapse();
  for (int step = 1; step <= 3; ++step) {
    const float phase = 1.3f * static_cast<float>(frame) + 2.1f * step;
    const math::Vector axis =
        math::Vector(cosf(phase), 0.5f, sinf(phase)).normalized();
    history.push(math::make_rotation(axis, 0.1f) * history.get());
  }
}

/** @brief Explicit-time rotation follows the sub-frame history at every
 *         column time, frame after frame. */
inline void test_world_orient_follows_history_at_column_times() {
  math::Orientation<> history(math::make_rotation(math::X_AXIS, 0.3f));
  Filter::World::Orient orient(history);
  orient.init_storage(persistent_arena, 288);
  StubEffect effect(288, 144);
  const math::Vector v = math::Vector(0.4f, -0.3f, 0.86f).normalized();
  float worst = 0.0f, endpoint = 0.0f;
  for (int frame = 0; frame < 4; ++frame) {
    push_turning_steps(history, frame);
    {
      Canvas canvas(effect);
      for (int c = 0; c < 144; ++c) {
        const float t = Spacetime::SamplingSurface<288>::time(c);
        const math::Vector want = math::rotate(v, history.get(t));
        worst = fmaxf(worst, (orient.at_time(canvas, v, t) - want).magnitude());
        const math::Quaternion slerped =
            math::slerp(history.get(0), history.get(), t);
        endpoint =
            fmaxf(endpoint, (math::rotate(v, slerped) - want).magnitude());
      }
    }
    effect.advance_display();
  }
  HS_EXPECT_LT(worst, 2e-6f);
  HS_EXPECT_GT(endpoint, 100.0f * 2e-6f);
}

/** @brief Explicit fractions follow the pose table and refresh next frame. */
inline void test_world_orient_explicit_time() {
  math::Orientation<> history = spin_history(-0.5f, -0.08f, 3);
  Filter::World::Orient orient(history);
  orient.init_storage(persistent_arena, 288);
  const math::Vector v(0.6f, -0.8f, 0.0f);
  StubEffect effect(288, 144);
  Pipeline<288, 144, Filter::World::Orient, Filter::Screen::AntiAlias<288, 144>>
      timed{Filter::World::Orient(history),
            Filter::Screen::AntiAlias<288, 144>()};
  timed.init_storage(persistent_arena);
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
      Spacetime::set_flat_surface(true);
      const math::Vector disabled = orient.at_time(canvas, v, 0.13f);
      Spacetime::set_flat_surface(false);
      HS_EXPECT_EQ(disabled.x, end.x);
      HS_EXPECT_EQ(disabled.y, end.y);
      HS_EXPECT_EQ(disabled.z, end.z);
    }
    effect.advance_display();
    history = spin_history(0.7f, 0.12f, 3);
  }
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
  PipelineRef ref(pipeline, Spacetime::Timed{});
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
  pipeline.init_storage(persistent_arena);
  PipelineRef ref(pipeline, Spacetime::Timed{});
  constexpr float ANGLE = 0.10f * math::PI_F;
  HS_EXPECT_GE(ref.sweep_arc(), ANGLE);
  HS_EXPECT_LE(ref.sweep_arc(), ANGLE + (history.length() - 1) * 0.000203f);
  StubEffect effect(288, 144);
  effect.set_margin(0);
  effect.set_clip(36, 72, 0, 144);
  Canvas canvas(effect);
  const auto exact = canvas.clip();
  const Plot::CullRegion still(canvas.clip());
  const Plot::CullRegion swept(canvas.clip(), ref.sweep_arc());
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
  const float expected = ref.sweep_arc() * (288.0f / (2.0f * math::PI_F)) /
                         (1.0f - ref.sweep_arc());
  HS_EXPECT_NEAR(swept.column_pad<288>(1.0f), expected, 1e-6f);
}

/** @brief Cached explicit poses invalidate on fraction, input, frame and collapse. */
inline void test_world_orient_sample_cache() {
  math::Orientation<> history = spin_history(-0.5f, -0.08f, 3);
  Filter::World::Orient cached(history);
  cached.init_storage(persistent_arena, 288);
  StubEffect effect(288, 144);
  const math::Vector a(0.6f, -0.8f, 0.0f);
  const math::Vector b(-0.2f, 0.3f, 0.7f);
  for (int frame = 0; frame < 2; ++frame) {
    {
      Canvas canvas(effect);
      auto compare = [&](const math::Vector &v, float t) {
        Filter::World::Orient fresh(history);
        fresh.init_storage(persistent_arena, 288);
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
      Spacetime::set_flat_surface(true);
      compare(a, 0.13f);
      Spacetime::set_flat_surface(false);
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
      HS_EXPECT_EQ(math::centered_remainder(value, period),
                   std::remainder(value, period));
    }
    for (int branch = -4; branch <= 4; ++branch)
      for (float offset : {-0.25f, 0.0f, 0.25f}) {
        const float value = (branch + 0.5f) * period + offset;
        HS_EXPECT_EQ(math::centered_remainder(value, period),
                     std::remainder(value, period));
      }
  }
  for (int quarter = -2304; quarter <= 2304; ++quarter) {
    const float column = quarter * 0.25f;
    const float turns = column / 144.0f;
    const float expected = turns - floorf(turns);
    HS_EXPECT_NEAR(Spacetime::SamplingSurface<288>::time(column), expected,
                   0.000001f);
    HS_EXPECT_GE(Spacetime::SamplingSurface<288>::time(column), 0.0f);
    HS_EXPECT_LT(Spacetime::SamplingSurface<288>::time(column), 1.0f);
  }
}

/** @brief Approximate orientation motion bounds remain conservative and preserve exact stillness. */
inline void test_temporal_orientation_bound_is_conservative() {
  const bool enabled = !Spacetime::flat_surface();
  Spacetime::set_flat_surface(false);
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
      const float actual = orient.sweep_arc();
      HS_EXPECT_GE(actual, expected);
      HS_EXPECT_LE(actual - expected, 0.000203f);
      if (index == 1024)
        HS_EXPECT_NEAR(actual, 0.0f, 0.0f);
    }
  math::Orientation<> history;
  history.set(math::Quaternion{});
  history.push(math::make_rotation(math::Y_AXIS, .5f));
  Filter::World::Orient orient(history);
  Spacetime::set_flat_surface(true);
  HS_EXPECT_NEAR(orient.sweep_arc(), 0.0f, 0.0f);
  Spacetime::set_flat_surface(!enabled);
}
