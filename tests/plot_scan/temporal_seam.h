/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Boundary taps use their own geometry, shade and birth fractions. */
inline void test_boundary_taps_use_column_time() {
  constexpr int W = 288, H = 144;
  const bool enabled = Temporal::temporal_sampling_enabled();
  Temporal::set_temporal_sampling_enabled(true);
  for (float boundary : {0.0f, W / 2.0f})
    for (float part_time : {0.0f, 1.0f})
      for (float birth : {0.0f, 0.5f}) {
        int evaluations = 0, taps = 0;
        const int column = part_time == 0.0f
                               ? static_cast<int>(boundary)
                               : (static_cast<int>(boundary) + W - 1) % W;
        const float at = float(column % (W / 2)) / (W / 2);
        const double x = boundary - 0.7f + 0.3f * at;
        const double y = 52.2f + 7.0f * at;
        const double fx = x - std::floor(x), fy = y - std::floor(y);
        const auto ease = [](double value) {
          return value * value * value * (value * (value * 6.0 - 15.0) + 10.0);
        };
        const double horizontal = part_time == 0.0f ? ease(fx) : 1.0 - ease(fx);
        const double opacity = 0.2f + 0.3f * at;
        auto evaluate = [&](float time, math::PixelCoords &position,
                            Pixel &color, float &age, float &alpha) {
          ++evaluations;
          position = {boundary - 0.7f + 0.3f * time, 52.2f + 7.0f * time};
          color = Pixel(static_cast<uint16_t>(10000.0f + 20000.0f * time), 2000,
                        3000);
          age = 2.0f;
          alpha = 0.2f + 0.3f * time;
          return time >= birth;
        };
        Temporal::BoundaryEvaluator evaluator(evaluate);
        Filter::Screen::AntiAlias<W, H> sink;
        sink.plot_at(
            boundary - 0.4f, 45.1f, Pixel(111, 222, 333), 9.0f, 0.7f, part_time,
            [&](float px, float py, const Pixel &color, float age,
                float alpha) {
              HS_EXPECT_EQ(static_cast<int>(px), column);
              HS_EXPECT_GE(static_cast<int>(py),
                           static_cast<int>(std::floor(y)));
              HS_EXPECT_LE(static_cast<int>(py),
                           static_cast<int>(std::floor(y)) + 1);
              const double vertical =
                  py == std::floor(y) ? 1.0 - ease(fy) : ease(fy);
              HS_EXPECT_NEAR(alpha, opacity * horizontal * vertical, 2e-6);
              HS_EXPECT_EQ(color.r,
                           static_cast<uint16_t>(10000.0f + 20000.0f * at));
              HS_EXPECT_NEAR(age, 2.0f, 0.0f);
              ++taps;
            },
            &evaluator);
        HS_EXPECT_EQ(evaluations, 1);
        HS_EXPECT_EQ(taps, at >= birth ? 2 : 0);
      }
  Temporal::set_temporal_sampling_enabled(enabled);
}

/** @brief Work counts from a diagnostic temporal world stage. */
struct RigidWorkCounts {
  int world = 0; /**< World evaluations. */
  int culls = 0; /**< Edge visibility tests. */
};

/** @brief Instrumented temporal stage with a selectable rigidity contract.
 * @tparam Rigid Whether control-point evaluation preserves its geometry. */
template <bool Rigid> struct CountedTemporalWorld : Filter::Is3D {
  static constexpr bool rigid_world_transform =
      Rigid; /**< Control-safe transform. */
  static constexpr bool has_temporal_motion =
      true;                /**< Varies during the frame. */
  RigidWorkCounts *counts; /**< Borrowed operation counters. */
  /** @brief Binds operation counters. */
  explicit CountedTemporalWorld(RigidWorkCounts &counts) : counts(&counts) {}
  /** @brief Bounds the frame's angular motion. */
  float cull_bound() const { return 0.01f; }
  /** @brief Evaluates the world transform at a frame fraction. */
  math::Vector at_time(const Canvas &, const math::Vector &p, float t) const {
    ++counts->world;
    if constexpr (Rigid)
      return math::rotate(p, math::make_rotation(math::Y_AXIS, 0.01f * t));
    else
      return math::Vector(p.x, p.y + 0.005f * t * p.x * p.x, p.z).normalized();
  }
  /** @brief Emits the endpoint transform. */
  template <typename Pass>
  void plot(const math::Vector &p, const Pixel &color, float age, float alpha,
            Pass &&pass) const {
    ++counts->world;
    if constexpr (Rigid)
      pass(math::rotate(p, math::make_rotation(math::Y_AXIS, 0.01f)), color,
           age, alpha);
    else
      pass(math::Vector(p.x, p.y + 0.005f * p.x * p.x, p.z).normalized(), color,
           age, alpha);
  }
  /** @brief Emits a transform at the supplied time. */
  template <typename Pass>
  void plot_at(Canvas &cv, const math::Vector &p, const Pixel &color, float age,
               float alpha, float t, Pass &&pass) const {
    pass(at_time(cv, p, t), color, age, alpha);
  }
  /** @brief Counts a source-edge visibility test. */
  template <typename Forward>
  bool cull_edge(const math::Vector &a, const math::Vector &b,
                 const math::Basis *basis, Forward &&forward) const {
    ++counts->culls;
    return forward(a, b, basis);
  }
  /** @brief Counts a timed source-edge visibility test. */
  template <typename Forward>
  bool cull_edge_timed(const math::Vector &a, const math::Vector &b,
                       const math::Basis *basis, Forward &&forward) const {
    return cull_edge(a, b, basis, forward);
  }
};

/** @brief One post-AA sample retained by the projected-entry fixture. */
struct ProjectedTap {
  float x, y;       /**< Pixel coordinates. */
  Pixel color;      /**< Sample color. */
  float age, alpha; /**< Sample age and opacity. */
};

/** @brief Records samples after screen filtering. */
struct ProjectedTapRecorder : Filter::Is2D {
  std::vector<ProjectedTap> *taps; /**< Borrowed output. */
  /** @brief Binds the destination sample list. */
  explicit ProjectedTapRecorder(std::vector<ProjectedTap> &taps)
      : taps(&taps) {}
  /** @brief Records a filtered screen sample. */
  template <typename Pass>
  void plot(float x, float y, const Pixel &color, float age, float alpha,
            Pass &&) {
    taps->push_back({x, y, color, age, alpha});
  }
};

/** @brief Projected typed and erased entries bypass world stages and retain tap timing. */
inline void test_timed_projected_entry_bypasses_world() {
  constexpr int W = 288, H = 144;
  for (bool erased : {false, true})
    for (float column : {40.3f, 143.6f, 287.6f})
      for (float time : {0.0f, 1.0f}) {
        StubEffect effect(W, H);
        Canvas canvas(effect);
        RigidWorkCounts counts;
        std::vector<ProjectedTap> actual, expected;
        Pipeline<W, H, CountedTemporalWorld<true>,
                 Filter::Screen::AntiAlias<W, H>, ProjectedTapRecorder>
            pipeline(CountedTemporalWorld<true>(counts),
                     Filter::Screen::AntiAlias<W, H>{},
                     ProjectedTapRecorder(actual));
        PipelineRef ref(pipeline, Temporal::Timed{});
        auto evaluate = [&](float at, math::PixelCoords &position, Pixel &color,
                            float &age, float &alpha) {
          position = {column - 0.1f + 0.2f * at, 50.2f + 5.0f * at};
          color = Pixel(static_cast<uint16_t>(10000 + 20000 * at), 2000, 3000);
          age = 2.0f;
          alpha = 0.2f + 0.3f * at;
          return at >= 0.5f;
        };
        Temporal::BoundaryEvaluator evaluator(evaluate);
        Temporal::ColumnPrediction prediction{
            column, std::numeric_limits<float>::infinity()};
        prediction.boundary = &evaluator;
        const Pixel color(111, 222, 333);
        Filter::Screen::AntiAlias<W, H> aa;
        aa.plot_at(
            column, 45.1f, color, 9.0f, 0.7f, time,
            [&](float x, float y, const Pixel &c, float age, float alpha) {
              expected.push_back({x, y, c, age, alpha});
            },
            &evaluator);
        counts.world = 0;
        const math::PixelCoords point{column, 45.1f};
        const float landed =
            erased ? ref.plot_projected_at(canvas, point, color, 9.0f, 0.7f,
                                           time, &prediction)
                   : pipeline.plot_projected_at(canvas, point, color, 9.0f,
                                                0.7f, time, &prediction);
        HS_EXPECT_NEAR(landed, column, 0.0f);
        HS_EXPECT_EQ(counts.world, 0);
        HS_EXPECT_EQ(actual.size(), expected.size());
        for (size_t i = 0; i < std::min(actual.size(), expected.size()); ++i) {
          HS_EXPECT_NEAR(actual[i].x, expected[i].x, 0.0f);
          HS_EXPECT_NEAR(actual[i].y, expected[i].y, 0.0f);
          HS_EXPECT_EQ(actual[i].color.r, expected[i].color.r);
          HS_EXPECT_EQ(actual[i].color.g, expected[i].color.g);
          HS_EXPECT_EQ(actual[i].color.b, expected[i].color.b);
          HS_EXPECT_NEAR(actual[i].age, expected[i].age, 0.0f);
          HS_EXPECT_NEAR(actual[i].alpha, expected[i].alpha, 0.000002f);
        }
      }
}

/** @brief Forced tap timing matches independently frozen interior columns. */
inline void test_interior_taps_use_column_time() {
  constexpr int W = 288, H = 144;
  const bool enabled = Temporal::temporal_sampling_enabled();
  Temporal::set_temporal_sampling_enabled(true);
  for (float birth : {0.0f, 40.5f / 144.0f}) {
    const auto render = [&](int route) {
      StubEffect effect(W, H);
      {
        Canvas canvas(effect);
        auto evaluate = [&](float at, math::PixelCoords &position, Pixel &color,
                            float &age, float &alpha) {
          position = {40.3f + 0.2f * at, 50.2f + 5.0f * at};
          color = Pixel(static_cast<uint16_t>(10000 + 20000 * at), 2000, 3000);
          age = 2.0f;
          alpha = 0.2f + 0.3f * at;
          return at >= birth;
        };
        Temporal::BoundaryEvaluator evaluator(evaluate);
        evaluator.force_timing = true;
        Temporal::ColumnPrediction prediction{
            40.3f, std::numeric_limits<float>::infinity()};
        prediction.boundary = &evaluator;
        if (route == 0) {
          Filter::Screen::AntiAlias<W, H> aa;
          Pipeline<W, H> terminal;
          for (int column : {40, 41}) {
            math::PixelCoords position;
            Pixel color;
            float age, alpha;
            if (evaluate(float(column) / (W / 2), position, color, age, alpha))
              aa.plot(position.x, position.y, color, age, alpha,
                      [&](float x, float y, const Pixel &c, float a,
                          float opacity) {
                        if (static_cast<int>(x) == column)
                          terminal.plot(canvas, x, y, c, a, opacity);
                      });
          }
        } else if (route == 1) {
          Filter::Screen::DirectAntiAliasSink<W, H> sink;
          sink.prepare(canvas);
          sink.plot_projected_at(canvas, {40.3f, 45.1f}, Pixel(111, 222, 333),
                                 9.0f, 0.7f, 40.3f / 144.0f, &prediction);
        } else {
          Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
          pipeline.plot_projected_at(canvas, {40.3f, 45.1f},
                                     Pixel(111, 222, 333), 9.0f, 0.7f,
                                     40.3f / 144.0f, &prediction);
        }
      }
      effect.advance_display();
      std::vector<uint16_t> pixels;
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
          const Pixel pixel = effect.get_pixel(x, y);
          pixels.insert(pixels.end(), {pixel.r, pixel.g, pixel.b});
        }
      return pixels;
    };
    const auto expected = render(0);
    HS_EXPECT_TRUE(std::any_of(expected.begin(), expected.end(),
                               [](uint16_t value) { return value != 0; }));
    HS_EXPECT_TRUE(render(1) == expected);
    HS_EXPECT_TRUE(render(2) == expected);
  }
  Temporal::set_temporal_sampling_enabled(enabled);
}

/** @brief Owned taps preserve support, wrapping, row coverage and tiny-weight rejection. */
inline void test_owned_tap_kernel_matches_ordinary_splat() {
  constexpr int W = 288, H = 144;
  const bool enabled = Temporal::temporal_sampling_enabled();
  Temporal::set_temporal_sampling_enabled(true);
  for (float x : {40.3f, -0.3f, 287.3f, 40.0001f})
    for (float shift : {-2.0f, 0.0f, 2.0f})
      for (float y : {-0.2f, 50.0001f, 143.8f})
        for (float time : {0.0f, 1.0f}) {
          std::vector<ProjectedTap> expected, actual;
          auto evaluate = [&](float at, math::PixelCoords &position,
                              Pixel &color, float &age, float &alpha) {
            position = {x + shift, y + 0.1f * at};
            color =
                Pixel(static_cast<uint16_t>(10000 + 20000 * at), 2000, 3000);
            age = 2.0f;
            alpha = 0.2f + 0.3f * at;
            return true;
          };
          Temporal::BoundaryEvaluator evaluator(evaluate);
          evaluator.force_timing = true;
          Filter::Screen::AntiAlias<W, H> aa;
          const int left = (static_cast<int>(std::floor(x)) + W) % W;
          for (int column : {left, (left + 1) % W}) {
            const float at = float(column % (W / 2)) / (W / 2);
            if (fabsf(at - time) > 0.5f)
              continue;
            math::PixelCoords position;
            Pixel color;
            float age, alpha;
            evaluate(at, position, color, age, alpha);
            aa.plot(position.x, position.y, color, age, alpha,
                    [&](float px, float py, const Pixel &c, float a,
                        float opacity) {
                      if (static_cast<int>(px) == column)
                        expected.push_back({px, py, c, a, opacity});
                    });
          }
          aa.plot_at(
              x, 45.1f, Pixel(111, 222, 333), 9.0f, 0.7f, time,
              [&](float px, float py, const Pixel &c, float age, float alpha) {
                actual.push_back({px, py, c, age, alpha});
              },
              &evaluator);
          HS_EXPECT_EQ(actual.size(), expected.size());
          for (size_t i = 0; i < std::min(actual.size(), expected.size());
               ++i) {
            HS_EXPECT_NEAR(actual[i].x, expected[i].x, 0.0f);
            HS_EXPECT_NEAR(actual[i].y, expected[i].y, 0.0f);
            HS_EXPECT_EQ(actual[i].color.r, expected[i].color.r);
            HS_EXPECT_EQ(actual[i].color.g, expected[i].color.g);
            HS_EXPECT_EQ(actual[i].color.b, expected[i].color.b);
            HS_EXPECT_NEAR(actual[i].age, expected[i].age, 0.0f);
            HS_EXPECT_NEAR(actual[i].alpha, expected[i].alpha, 0.0f);
          }
        }
  {
    Filter::Screen::AntiAlias<1, 2> aa;
    std::vector<float> expected, actual;
    auto evaluate = [](float, math::PixelCoords &position, Pixel &color,
                       float &age, float &alpha) {
      position = {0.3f, 0.4f};
      color = Pixel(10000, 2000, 3000);
      age = 2.0f;
      alpha = 0.4f;
      return true;
    };
    Temporal::BoundaryEvaluator evaluator(evaluate);
    evaluator.force_timing = true;
    for (int duplicate = 0; duplicate < 2; ++duplicate)
      aa.plot(0.3f, 0.4f, Pixel(10000, 2000, 3000), 2.0f, 0.4f,
              [&](float, float, const Pixel &, float, float alpha) {
                expected.push_back(alpha);
              });
    aa.plot_at(
        0.3f, 0.4f, Pixel(10000, 2000, 3000), 2.0f, 0.4f, 0.0f,
        [&](float, float, const Pixel &, float, float alpha) {
          actual.push_back(alpha);
        },
        &evaluator);
    HS_EXPECT_EQ(expected.size(), size_t{8});
    HS_EXPECT_TRUE(actual == expected);
  }
  Temporal::set_temporal_sampling_enabled(enabled);
}
