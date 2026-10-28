/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Explicit-time plotting and conservative orientation culling
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

/** @brief Whether a sampled history edge lights the exact clip.
 * @param clip Exact raster clip.
 * @param history Frame rotation history.
 * @param a,b Source edge endpoints.
 * @return True when any independently projected sample has a visible tap. */
inline bool history_edge_lights_clip(const ClipRegion &clip,
                                     const math::Orientation<> &history,
                                     const math::Vector &a,
                                     const math::Vector &b) {
  constexpr double PI_VALUE = 3.14159265358979323846;
  for (int time = 0; time <= 144; ++time) {
    const auto pose = history.get(static_cast<float>(time) / 144.0f);
    for (int sample = 0; sample <= 32; ++sample) {
      const float s = static_cast<float>(sample) / 32.0f;
      const auto v = math::rotate((a * (1.0f - s) + b * s).normalized(), pose);
      const double norm =
          std::sqrt(double(v.x) * v.x + double(v.y) * v.y + double(v.z) * v.z);
      double azimuth = std::atan2(double(v.z), double(v.x));
      if (azimuth < 0.0)
        azimuth += 2.0 * PI_VALUE;
      const double column = azimuth * 288.0 / (2.0 * PI_VALUE);
      const float row = math::phi_to_y<144>(
          static_cast<float>(std::acos(std::clamp(v.y / norm, -1.0, 1.0))));
      const int x = static_cast<int>(std::floor(column));
      const int y = static_cast<int>(std::floor(row));
      for (int dx : {0, 1})
        for (int dy : {0, 1})
          if (clip.contains_y(y + dy) && !clip.x_clip().clipped((x + dx) % 288))
            return true;
    }
  }
  return false;
}

/** @brief End-pose padding keeps an edge visible earlier in its history. */
inline void test_oriented_cull_keeps_visible_history_edges() {
  StubEffect effect(288, 144);
  effect.set_margin(0);
  effect.set_clip(64, 72, 64, 80);
  Canvas canvas(effect);
  const auto entering = spin_history(0.0f, 0.2f, 3);
  Pipeline<288, 144, Filter::World::Orient> explicit_case{
      Filter::World::Orient(entering)};
  const auto a = math::pixel_to_vector<288, 144>(67.0f, 69.0f);
  const auto b = math::pixel_to_vector<288, 144>(69.0f, 73.0f);
  HS_EXPECT_TRUE(history_edge_lights_clip(canvas.clip(), entering, a, b));
  const Plot::CullRegion endpoint(canvas.clip());
  HS_EXPECT_FALSE((Plot::geodesic_reaches_region<288, 144>(
      endpoint, math::rotate(a, entering.get()),
      math::rotate(b, entering.get()))));
  HS_EXPECT_TRUE((Plot::edge_visible_in_clip<288, 144>(
      explicit_case, canvas.clip(), canvas.clip().x_clip(), a, b, nullptr)));
  PipelineRef erased(explicit_case);
  HS_EXPECT_TRUE((Plot::edge_visible_in_clip<288, 144>(
      erased, canvas.clip(), canvas.clip().x_clip(), a, b, nullptr)));
  const std::array<math::Orientation<>, 2> candidates{
      math::Orientation<>(math::make_rotation(math::Y_AXIS, 2.0f)), entering};
  Pipeline<288, 144, Filter::World::OrientSlice> explicit_slices{
      Filter::World::OrientSlice(candidates, math::Y_AXIS)};
  HS_EXPECT_NEAR(explicit_slices.sweep_arc(), explicit_case.sweep_arc(), 1e-6f);
  HS_EXPECT_TRUE((Plot::edge_visible_in_clip<288, 144>(
      explicit_slices, canvas.clip(), canvas.clip().x_clip(), a, b, nullptr)));

  hs::random().seed(0xC011B0u);
  int visible = 0, rejected = 0;
  for (int trial = 0; trial < 160; ++trial) {
    math::Orientation<> history;
    for (int step = 1; step <= 3; ++step)
      history.push(math::make_rotation(rand_unit(), hs::rand_f(-0.3f, 0.3f)) *
                   history.get());
    const auto from = rand_unit();
    const auto to = (from + rand_unit() * 0.3f).normalized();
    Pipeline<288, 144, Filter::World::Orient> pipeline{
        Filter::World::Orient(history)};
    Pipeline<288, 144, Filter::World::OrientSlice> slices{
        Filter::World::OrientSlice(
            std::span<const math::Orientation<>>(&history, 1), math::Y_AXIS)};
    const bool oracle =
        history_edge_lights_clip(canvas.clip(), history, from, to);
    const bool got = Plot::edge_visible_in_clip<288, 144>(
        pipeline, canvas.clip(), canvas.clip().x_clip(), from, to, nullptr);
    const bool sliced = Plot::edge_visible_in_clip<288, 144>(
        slices, canvas.clip(), canvas.clip().x_clip(), from, to, nullptr);
    PipelineRef erased_history(pipeline);
    const bool ERASED_VISIBLE = Plot::edge_visible_in_clip<288, 144>(
        erased_history, canvas.clip(), canvas.clip().x_clip(), from, to,
        nullptr);
    if (oracle) {
      ++visible;
      HS_EXPECT_TRUE(got);
      HS_EXPECT_TRUE(sliced);
      HS_EXPECT_TRUE(ERASED_VISIBLE);
    }
    rejected += !got;
  }
  HS_EXPECT_GT(visible, 0);
  HS_EXPECT_GT(rejected, 0);
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
  const ScopedTemporalSampling temporal_sampling(true);
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
}
