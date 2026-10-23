/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

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
        const Pixel color(111, 222, 333);
        Filter::Screen::AntiAlias<W, H> aa;
        aa.plot_at(
            column, 45.1f, color, 9.0f, 0.7f, time,
            [&](float x, float y, const Pixel &c, float age, float alpha) {
              expected.push_back({x, y, c, age, alpha});
            });
        counts.world = 0;
        const math::PixelCoords point{column, 45.1f};
        const float landed =
            erased
                ? ref.plot_projected_at(canvas, point, color, 9.0f, 0.7f, time)
                : pipeline.plot_projected_at(canvas, point, color, 9.0f, 0.7f,
                                             time);
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
