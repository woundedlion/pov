/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Direct sink selecting rigid controls or the per-sample fallback.
 * @tparam W,H Canvas dimensions.
 * @tparam Rigid Enables transformed control-point interpolation. */
template <int W, int H, bool Rigid>
struct SeamRenderSink : Filter::Screen::DirectAntiAliasSink<W, H> {
  static constexpr bool rigid_world_stages =
      Rigid; /**< Control interpolation policy. */
};

/** @brief Renders a callback's closed contour at a fixed or capture time.
 * @param geometry Control-point positions at each frame fraction.
 * @param count Number of contour controls.
 * @param bound Angular motion bound.
 * @param timed Selects the timed walk.
 * @param fraction Fraction used by the frozen renderer.
 * @param displacement Frozen geometry's vertical displacement in pixels.
 * @param closed Connects the final control point to the first.
 * @return Row-major linear luminance. */
template <int W, int H, bool Rigid = true, typename Geometry>
std::vector<float> seam_render(Geometry &geometry, size_t count, float bound,
                               bool timed, float fraction,
                               float displacement = 0.0f, bool closed = true) {
  hs_test::StubEffect effect(W, H);
  SeamRenderSink<W, H, Rigid> sink;
  ScratchScope scratch(plot_arena());
  Fragments points;
  points.bind(plot_arena(), count);
  for (size_t i = 0; i < count; ++i) {
    Fragment point;
    point.pos = geometry(i, timed ? 1.0f : fraction);
    if (displacement != 0.0f)
      point.pos = hs_test::displace_reference_rows<H>(point.pos, displacement);
    points.push_back(point);
  }
  {
    Canvas canvas(effect);
    sink.prepare(canvas);
    if (timed)
      Plot::rasterize<W, H, {.timed = true}>(
          sink, canvas, points, white_shader(0.4f),
          {.loop = closed ? Plot::RasterLoop::closed() : Plot::RasterLoop{}},
          Plot::GeometryCallback::interpolating(geometry), bound);
    else
      Plot::rasterize<W, H>(
          sink, canvas, points, white_shader(0.4f),
          {.loop = closed ? Plot::RasterLoop::closed() : Plot::RasterLoop{}});
  }
  effect.advance_display();
  std::vector<float> luminance(W * H);
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const Pixel pixel = effect.get_pixel(x, y);
      luminance[y * W + x] = (pixel.r + pixel.g + pixel.b) / (3.0f * 65535.0f);
    }
  return luminance;
}

/** @brief Compares a source contour with its frozen column renders.
 * @param geometry Continuous source geometry.
 * @param count Number of contour controls.
 * @param bound Angular motion bound.
 * @param name Diagnostic label. */
template <int W, int H, typename Geometry>
void seam_compare(Geometry &geometry, size_t count, float bound,
                  const char *name) {
  const auto render = [&](bool timed, float fraction, float displacement) {
    return seam_render<W, H>(geometry, count, bound, timed, fraction,
                             displacement);
  };
  seam_compare_render<W, H>(render, name);
}

/** @brief A closed rotating ring agrees with its per-column frozen reference. */
inline void test_timed_ring_matches_column_time_reference() {
  constexpr int W = 288, H = 144, COUNT = 144;
  const math::Vector axis = math::Vector(0.25f, 1.0f, 0.15f).normalized();
  const math::Vector normal = math::Vector(0.5f, 0.7f, 0.3f).normalized();
  const math::Vector a = math::cross(normal, math::Y_AXIS).normalized();
  const math::Vector b = math::cross(normal, a);
  for (float direction : {-1.0f, 1.0f}) {
    const float angle = direction * 3.0f * math::RADIANS_PER_COLUMN<W>;
    auto geometry = [&](size_t index, float t) {
      const float theta = 2.0f * math::PI_F * index / COUNT;
      const math::Vector point = a * cosf(theta) + b * sinf(theta);
      const math::Vector previous =
          math::rotate(point, math::make_rotation(axis, -angle));
      return (previous * (1.0f - t) + point * t).normalized();
    };
    seam_compare<W, H>(geometry, COUNT, fabsf(angle),
                       direction < 0.0f ? "ring against" : "ring with");
  }
}

/** @brief A still boundary straddler retains each column's ordinary coverage. */
inline void test_timed_still_straddler_preserves_coverage() {
  constexpr int W = 288, H = 144;
  for (float boundary : {0.0f, W / 2.0f}) {
    auto geometry = [&](size_t index, float) {
      return math::pixel_to_vector<W, H>(boundary + (index ? 4.0f : -4.0f),
                                         58.25f);
    };
    const auto actual =
        seam_render<W, H>(geometry, 2, 0.0f, true, 1.0f, 0.0f, false);
    const auto frozen =
        seam_render<W, H>(geometry, 2, 0.0f, false, 1.0f, 0.0f, false);
    for (int delta : {-1, 0}) {
      const int column = (static_cast<int>(boundary) + delta + W) % W;
      float timed_sum = 0.0f, still_sum = 0.0f;
      for (int y = 0; y < H; ++y) {
        timed_sum += actual[y * W + column];
        still_sum += frozen[y * W + column];
      }
      std::printf("  straddler column=%d timed=%.6f still=%.6f\n", column,
                  timed_sum, still_sum);
      HS_EXPECT_GT(still_sum, 0.0f);
      HS_EXPECT_NEAR(timed_sum, still_sum, 1.0f / 65535.0f);
    }
  }
}

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

/** @brief Timed sink retaining the time and position of emitted samples. */
struct SeamTimedSink : Filter::Screen::DirectAntiAliasSink<288, 144> {
  std::vector<float> times;            /**< Emitted frame fractions. */
  std::vector<math::Vector> positions; /**< Emitted world points. */
  /** @brief Records a timed emission without splatting it.
   * @param position World point.
   * @param time Sample fraction. */
  void plot_at(Canvas &, const math::Vector &position, const Pixel &, float,
               float, float time) {
    times.push_back(time);
    positions.push_back(position);
  }
  /** @brief Records an accepted transformed sample and returns its column. */
  float plot_transformed_at(Canvas &canvas, const math::Vector &position,
                            const Pixel &color, float age, float alpha,
                            float time,
                            const Temporal::ColumnPrediction *prediction) {
    float column = math::vector_to_pixel<288, 144>(position).x;
    if (prediction)
      column = prediction->unwrap(column, 288.0f);
    if (!prediction ||
        (prediction->record_sample && prediction->accepts(column)))
      plot_at(canvas, position, color, age, alpha, time);
    return column;
  }
};

/** @brief Traces the per-sample fallback for a nonrigid world path. */
struct SeamPredictionSink : SeamTimedSink {
  static constexpr bool rigid_world_stages =
      false; /**< Requires sample-time evaluation. */
};

/** @brief Timed samples preserve birth and capture time through control refinement. */
inline void test_timed_registers_preserve_birth_and_time() {
  constexpr int W = 288, H = 144;
  for (bool interpolating : {false, true}) {
    hs_test::StubEffect effect(W, H);
    SeamTimedSink sink;
    Canvas canvas(effect);
    sink.prepare(canvas);
    ScratchScope scratch(plot_arena());
    Fragments points;
    points.bind(plot_arena(), 2);
    auto geometry = [&](size_t index, float t) {
      HS_EXPECT_LT(index, size_t{2});
      return math::pixel_to_vector<W, H>(
          (index ? 80.0f : 40.0f) - 6.0f * (1.0f - t), 65.0f);
    };
    for (size_t i = 0; i < 2; ++i) {
      Fragment point;
      point.pos = geometry(i, 1.0f);
      points.push_back(point);
    }
    const float births[] = {0.4f, 0.4f};
    auto callback = interpolating
                        ? Plot::GeometryCallback::interpolating(geometry)
                        : Plot::GeometryCallback::exact(geometry);
    callback.births = births;
    std::vector<float> shaded;
    auto shader = [&](const math::Vector &, Fragment &fragment, float t) {
      shaded.push_back(t);
      fragment.color = Color4(Pixel(65535, 65535, 65535), 0.4f);
    };
    Temporal::geometry_counters = {};
    Plot::rasterize<W, H, {.interpolate_registers = false, .timed = true}>(
        sink, canvas, points, shader, {}, callback,
        6.0f * math::RADIANS_PER_COLUMN<W>);
    HS_EXPECT_GT(sink.times.size(), size_t{5});
    HS_EXPECT_GE(shaded.size(), sink.times.size());
    for (size_t i = 0; i < sink.times.size(); ++i) {
      HS_EXPECT_GE(sink.times[i], 0.4f);
      HS_EXPECT_LT(sink.times[i], 0.6f);
      HS_EXPECT_TRUE(std::find(shaded.begin(), shaded.end(), sink.times[i]) !=
                     shaded.end());
      const auto pixel = math::vector_to_pixel<W, H>(sink.positions[i]);
      HS_EXPECT_NEAR(sink.times[i], Temporal::column_fraction<W>(pixel.x),
                     0.004f);
    }
  }
}

/** @brief A closed interpolating contour evaluates each previous lattice vertex once. */
inline void test_timed_lattice_vertices_evaluated_once() {
  constexpr int W = 288, H = 144, COUNT = 16;
  std::array<int, COUNT> previous_calls{};
  auto geometry = [&](size_t index, float t) {
    HS_EXPECT_LT(index, size_t{COUNT});
    if (t == 0.0f)
      ++previous_calls[index];
    return math::pixel_to_vector<W, H>(
        12.0f + 18.0f * index - 6.0f * (1.0f - t), 65.0f);
  };
  const auto output = seam_render<W, H>(
      geometry, COUNT, 6.0f * math::RADIANS_PER_COLUMN<W>, true, 1.0f);
  float energy = 0.0f;
  for (float value : output)
    energy += value;
  HS_EXPECT_GT(energy, 1.0f);
  for (int calls : previous_calls)
    HS_EXPECT_EQ(calls, 1);
}

/** @brief A closed contour starts callback evaluation at its first boundary cut. */
inline void test_timed_closed_ring_starts_at_cut() {
  constexpr int W = 288, H = 144, COUNT = 16;
  size_t first_previous = COUNT;
  auto geometry = [&](size_t index, float t) {
    if (t == 0.0f && first_previous == COUNT)
      first_previous = index;
    return math::pixel_to_vector<W, H>(
        12.0f + 18.0f * index - 6.0f * (1.0f - t), 65.0f);
  };
  seam_render<W, H>(geometry, COUNT, 6.0f * math::RADIANS_PER_COLUMN<W>, true,
                    1.0f);
  HS_EXPECT_EQ(first_previous, size_t{6});
}

/** @brief Moving cuts retain distinct previous-start and current-end crossings. */
inline void test_timed_cut_uses_separate_endpoint_crossings() {
  constexpr int W = 288, H = 144;
  for (float boundary : {0.0f, W / 2.0f}) {
    for (float motion : {-6.0f, 6.0f}) {
      auto geometry = [&](size_t index, float t) {
        return math::pixel_to_vector<W, H>(
            boundary + (index ? 4.0f : -4.0f) - motion * (1.0f - t), 65.0f);
      };
      const auto output =
          seam_render<W, H>(geometry, 2, 6.0f * math::RADIANS_PER_COLUMN<W>,
                            true, 1.0f, 0.0f, false);
      const int early_column = static_cast<int>(boundary) + 3;
      const int late_column = (static_cast<int>(boundary) + W - 2) % W;
      float early = 0.0f, late = 0.0f;
      for (int row = 0; row < H; ++row) {
        early += output[row * W + early_column];
        late += output[row * W + late_column];
      }
      if (motion < 0.0f)
        HS_EXPECT_GT(early, 0.2f);
      else
        HS_EXPECT_LT(early, 1.0f / 65535.0f);
      HS_EXPECT_GT(late, 0.2f);
    }
  }
}

/** @brief Polar held anchors still draw both sides of an arm-boundary cut. */
inline void test_timed_held_anchor_keeps_boundary_parts() {
  constexpr int W = 288, H = 144;
  hs_test::StubEffect effect(W, H);
  Canvas canvas(effect);
  Filter::Screen::DirectAntiAliasSink<W, H> sink;
  sink.prepare(canvas);
  ScratchScope guard(scratch_arena_a);
  Fragments points;
  points.bind(scratch_arena_a, 2);
  for (float column : {-4.0f, 4.0f}) {
    Fragment point;
    point.pos = math::pixel_to_vector<W, H>(column, 0.5f);
    points.push_back(point);
  }
  auto geometry = [&](size_t index, float) { return points[index].pos; };
  Temporal::geometry_counters = {};
  Plot::rasterize<W, H, {.timed = true}>(
      sink, canvas, points, white_shader(0.4f), {},
      Plot::GeometryCallback::exact(geometry), 0.1f);
  const auto counters = Temporal::geometry_counters;
  HS_EXPECT_GT(counters.held_anchors, 0u);
  HS_EXPECT_EQ(counters.parts, 2u);
  HS_EXPECT_GT(counters.solve_probes, 0u);
}

/** @brief One correction retains dense coverage and bounded capture-time error. */
inline void test_timed_prediction_correction_preserves_coverage() {
  constexpr int W = 288, H = 144;
  constexpr float MOTION = 30.0f;
  hs_test::StubEffect effect(W, H);
  SeamPredictionSink sink;
  Canvas canvas(effect);
  sink.prepare(canvas);
  ScratchScope scratch(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 2);
  auto geometry = [](size_t index, float t) {
    return math::pixel_to_vector<W, H>(
        (index ? 120.0f : 40.0f) - MOTION * (1.0f - t), 10.0f);
  };
  for (size_t i = 0; i < 2; ++i) {
    Fragment point;
    point.pos = geometry(i, 1.0f);
    points.push_back(point);
  }
  auto shader = [](const math::Vector &, Fragment &fragment, float) {
    fragment.color = Color4(Pixel(65535, 65535, 65535), 0.4f);
  };
  Temporal::geometry_counters = {};
  Plot::rasterize<W, H, {.timed = true}>(
      sink, canvas, points, shader, {}, Plot::GeometryCallback::exact(geometry),
      MOTION * math::RADIANS_PER_COLUMN<W> * sinf(math::y_to_phi<H>(10.0f)));
  HS_EXPECT_GT(Temporal::geometry_counters.corrections, uint32_t{0});
  HS_EXPECT_EQ(Temporal::geometry_counters.held_anchors, uint32_t{0});
  HS_EXPECT_GT(sink.positions.size(), size_t{2});
  for (size_t i = 0; i < sink.positions.size(); ++i) {
    const auto pixel = math::vector_to_pixel<W, H>(sink.positions[i]);
    const float delta =
        fabsf(sink.times[i] - Temporal::column_fraction<W>(pixel.x));
    HS_EXPECT_LE(MOTION * fminf(delta, 1.0f - delta), 0.1f);
    if (i == 0)
      continue;
    const auto previous = math::vector_to_pixel<W, H>(sink.positions[i - 1]);
    const float dx = std::remainder(pixel.x - previous.x, float(W));
    const float dy = pixel.y - previous.y;
    HS_EXPECT_LT(dx * dx + dy * dy, 4.0f);
  }
}

/** @brief A coherent interior contour probes only its first source point. */
inline void test_timed_prediction_reuses_landed_columns() {
  constexpr int W = 288, H = 144, COUNT = 16;
  auto geometry = [](size_t index, float t) {
    return math::pixel_to_vector<W, H>(30.0f + 5.0f * index - 6.0f * (1.0f - t),
                                       65.0f);
  };
  Temporal::geometry_counters = {};
  const auto output = seam_render<W, H, false>(
      geometry, COUNT, 6.0f * math::RADIANS_PER_COLUMN<W>, true, 1.0f, 0.0f,
      false);
  float energy = 0.0f;
  for (float value : output)
    energy += value;
  HS_EXPECT_GT(energy, 1.0f);
  HS_EXPECT_EQ(Temporal::geometry_counters.end_probes, uint32_t{1});
  HS_EXPECT_EQ(Temporal::geometry_counters.refined_controls, uint32_t{0});
}
/** @brief Cut parts retain the original edge's sampling phase. */
inline void test_timed_cut_preserves_sampling_phase() {
  constexpr int W = 288, H = 144;
  for (float boundary : {0.0f, W / 2.0f}) {
    auto geometry = [&](size_t index, float) {
      return math::pixel_to_vector<W, H>(boundary + (index ? 4.0f : -4.0f),
                                         58.25f);
    };
    const auto actual =
        seam_render<W, H>(geometry, 2, 6.0f * math::RADIANS_PER_COLUMN<W>, true,
                          1.0f, 0.0f, false);
    const auto frozen =
        seam_render<W, H>(geometry, 2, 0.0f, false, 1.0f, 0.0f, false);
    for (int delta : {-1, 0, 1}) {
      const int column = (static_cast<int>(boundary) + delta + W) % W;
      float expected = 0.0f, observed = 0.0f;
      for (int row = 0; row < H; ++row) {
        expected += frozen[row * W + column];
        observed += actual[row * W + column];
      }
      HS_EXPECT_NEAR(observed, expected, 1.0f / 65535.0f);
    }
  }
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

/** @brief Rigid work scales with controls while nonrigid work remains per sample. */
inline void test_timed_rigid_work_scales_with_controls() {
  constexpr int W = 288, H = 144;
  auto check = [&]<bool Rigid>() {
    StubEffect effect(W, H);
    effect.set_clip(36, 110, 0, H);
    Canvas canvas(effect);
    RigidWorkCounts counts;
    Pipeline<W, H, CountedTemporalWorld<Rigid>, Filter::Screen::AntiAlias<W, H>>
        pipeline(CountedTemporalWorld<Rigid>(counts),
                 Filter::Screen::AntiAlias<W, H>{});
    ScratchScope scratch(plot_arena());
    Fragments points;
    points.bind(plot_arena(), 2);
    int calls = 0, samples = 0;
    auto geometry = [&](size_t i, float t) {
      ++calls;
      return math::pixel_to_vector<W, H>(40.0f + 60.0f * i + 0.5f * t, 60.0f);
    };
    for (size_t i = 0; i < 2; ++i) {
      Fragment point;
      point.pos = geometry(i, 1.0f);
      points.push_back(point);
    }
    calls = 0;
    auto shader = [&](const math::Vector &, Fragment &f, float) {
      ++samples;
      f.color = Color4(1, 1, 1, 0.4f);
    };
    Temporal::geometry_counters = {};
    Plot::rasterize<W, H, {.timed = true}>(
        pipeline, canvas, points, shader, {},
        Plot::GeometryCallback::exact(geometry),
        0.5f * math::RADIANS_PER_COLUMN<W>);
    HS_EXPECT_GT(samples, 40);
    if constexpr (Rigid) {
      const int taps = int(Temporal::geometry_counters.boundary_evaluations);
      HS_EXPECT_LT(calls - 2 * taps, samples - taps);
      HS_EXPECT_LT(counts.world - taps, samples - taps);
      HS_EXPECT_GT(counts.culls, 0);
      HS_EXPECT_LE(counts.culls, 2);
      HS_EXPECT_GT(Temporal::geometry_counters.rigid_samples, 40u);
    } else {
      HS_EXPECT_GT(calls, 40);
      HS_EXPECT_GT(counts.world, 40);
      HS_EXPECT_EQ(Temporal::geometry_counters.control_evaluations, 0u);
      HS_EXPECT_EQ(Temporal::geometry_counters.rigid_samples, 0u);
    }
  };
  check.template operator()<true>();
  check.template operator()<false>();
}
/** @brief Timed rigid replication preserves independent spatial copies. */
inline void test_timed_rigid_replicate_keeps_copies() {
  constexpr int W = 288, H = 144;
  StubEffect effect(W, H);
  Canvas canvas(effect);
  Pipeline<W, H, Filter::World::Replicate<W>, Filter::Screen::AntiAlias<W, H>>
      pipeline(Filter::World::Replicate<W>(2),
               Filter::Screen::AntiAlias<W, H>{});
  ScratchScope scratch(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 2);
  auto geometry = [](size_t i, float t) {
    return math::pixel_to_vector<W, H>(40.0f + 60.0f * i + 0.5f * t, 60.0f);
  };
  for (size_t i = 0; i < 2; ++i) {
    Fragment point;
    point.pos = geometry(i, 1.0f);
    points.push_back(point);
  }
  auto shader = [](const math::Vector &, Fragment &f, float) {
    f.color = Color4(1, 1, 1, 0.4f);
  };
  Temporal::geometry_counters = {};
  Plot::rasterize<W, H, {.timed = true}>(
      pipeline, canvas, points, shader, {},
      Plot::GeometryCallback::exact(geometry),
      0.5f * math::RADIANS_PER_COLUMN<W>);
  for (int start : {40, 184}) {
    int lit_columns = 0;
    for (int x = start; x < start + 61; ++x) {
      bool lit = false;
      for (int y = 0; y < H; ++y)
        lit |= canvas(x, y).r != 0;
      lit_columns += lit;
    }
    HS_EXPECT_GT(lit_columns, 55);
  }
  HS_EXPECT_GT(Temporal::geometry_counters.rigid_samples, 100u);
}
/** @brief Rigid endpoint transforms preserve the shader's source-space position. */
inline void test_timed_rigid_shader_reads_source_position() {
  constexpr int W = 288, H = 144;
  StubEffect effect(W, H);
  Canvas canvas(effect);
  math::Orientation<> orientation(math::make_rotation(math::X_AXIS, 0.5f));
  orientation.push(math::make_rotation(math::X_AXIS, 0.51f));
  Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>>
      pipeline(Filter::World::Orient(orientation),
               Filter::Screen::AntiAlias<W, H>{});
  pipeline.init_storage(persistent_arena);
  ScratchScope scratch(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 2);
  auto geometry = [](size_t i, float t) {
    const float angle =
        (40.0f + 60.0f * i + 0.5f * t) * math::RADIANS_PER_COLUMN<W>;
    return math::Vector(sinf(angle), 0.0f, cosf(angle));
  };
  for (size_t i = 0; i < 2; ++i) {
    Fragment point;
    point.pos = geometry(i, 1.0f);
    point.v0 = static_cast<float>(i);
    points.push_back(point);
  }
  int samples = 0;
  auto shader = [&](const math::Vector &position, Fragment &fragment, float t) {
    ++samples;
    (void)t;
    HS_EXPECT_NEAR(position.y, 0.0f, 0.000001f);
    HS_EXPECT_NEAR((fragment.pos - position).magnitude(), 0.0f, 0.0f);
    fragment.color = Color4(1, 1, 1, 0.4f);
  };
  Plot::rasterize<W, H, {.timed = true}>(
      pipeline, canvas, points, shader, {},
      Plot::GeometryCallback::exact(geometry),
      0.5f * math::RADIANS_PER_COLUMN<W>);
  HS_EXPECT_GT(samples, 40);
}
/** @brief Refined rigid spans keep capture time bounded on a curved interior edge. */
inline void test_timed_rigid_refinement_bounds_capture_time() {
  constexpr int W = 288, H = 144;
  constexpr float MOTION = 3.0f;
  StubEffect effect(W, H);
  Canvas canvas(effect);
  SeamTimedSink sink;
  sink.prepare(canvas);
  ScratchScope scratch(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 2);
  auto geometry = [](size_t index, float t) {
    return math::pixel_to_vector<W, H>(
        30.0f + 85.0f * index - MOTION * (1.0f - t), 20.0f);
  };
  for (size_t i = 0; i < 2; ++i) {
    Fragment point;
    point.pos = geometry(i, 1.0f);
    points.push_back(point);
  }
  Temporal::geometry_counters = {};
  Plot::rasterize<W, H, {.timed = true}>(
      sink, canvas, points, white_shader(0.4f), {},
      Plot::GeometryCallback::exact(geometry),
      MOTION * math::RADIANS_PER_COLUMN<W> * sinf(math::y_to_phi<H>(20.0f)));
  HS_EXPECT_GT(Temporal::geometry_counters.refined_controls, 0u);
  HS_EXPECT_EQ(Temporal::geometry_counters.rigid_fallback_edges, 0u);
  HS_EXPECT_GT(sink.times.size(), size_t{30});
  for (size_t i = 0; i < sink.times.size(); ++i) {
    const auto projected = math::vector_to_pixel<W, H>(sink.positions[i]);
    const float time_error =
        fabsf(sink.times[i] - Temporal::column_fraction<W>(projected.x));
    HS_EXPECT_LE(time_error * MOTION, 0.1f);
  }
}
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

/** @brief Dense selection rejects sparse controls and fills visible sample gaps. */
inline void test_timed_dense_selection_and_gap_coverage() {
  constexpr int W = 288, H = 144;
  for (float spacing : {0.7f, 1.8f, 5.0f}) {
    auto geometry = [spacing](size_t index, float t) {
      return math::pixel_to_vector<W, H>(40.2f + spacing * index + 0.1f * t,
                                         60.2f);
    };
    Temporal::geometry_counters = {};
    const auto pixels =
        seam_render<W, H>(geometry, 5, 0.1f * math::RADIANS_PER_COLUMN<W>, true,
                          1.0f, 0.0f, false);
    const auto counts = Temporal::geometry_counters;
    HS_EXPECT_EQ(counts.dense_runs, spacing < 2.0f ? 1u : 0u);
    if (spacing < 2.0f) {
      HS_EXPECT_EQ(counts.dense_controls, 5u);
      HS_EXPECT_EQ(counts.parts, 0u);
      HS_EXPECT_EQ(counts.solve_probes, 0u);
      HS_EXPECT_EQ(counts.dense_culls, 1u);
    }
    if (spacing > 1.0f && spacing < 2.0f) {
      HS_EXPECT_GT(counts.dense_gap_samples, 0u);
      for (int x = 41; x <= 47; ++x) {
        float brightness = 0;
        for (int y = 58; y <= 63; ++y)
          brightness += pixels[y * W + x];
        HS_EXPECT_GT(brightness, 0.16f);
      }
    }
  }
}

/** @brief Dense boundary taps evaluate births and shading at their owning column time. */
inline void test_timed_dense_boundary_matches_frozen_taps() {
  constexpr int W = 288, H = 144;
  for (float boundary : {0.0f, 144.0f}) {
    auto geometry = [boundary](size_t i, float t) {
      return math::pixel_to_vector<W, H>(boundary - 0.6f + 0.2f * t,
                                         50.2f + 0.3f * i + 5.0f * t);
    };
    auto render = [&](bool reference) {
      StubEffect effect(W, H);
      {
        Canvas canvas(effect);
        Filter::Screen::DirectAntiAliasSink<W, H> sink;
        sink.prepare(canvas);
        ScratchScope scratch(plot_arena());
        Fragments points;
        points.bind(plot_arena(), 3);
        for (size_t i = 0; i < 3; ++i) {
          Fragment f;
          f.pos = geometry(i, 1.0f);
          points.push_back(f);
        }
        auto shader = [](const math::Vector &, Fragment &f, float t) {
          f.color = Color4(static_cast<uint8_t>(50 + 120 * t), 100, 30,
                           0.2f + 0.3f * t);
        };
        if (!reference) {
          const std::array<float, 3> births{0.5f, 0.5f, 0.5f};
          auto callback = Plot::GeometryCallback::interpolating(geometry);
          callback.births = births;
          Plot::rasterize<W, H, {.timed = true}>(sink, canvas, points, shader,
                                                 {}, callback, 0.12f);
        } else {
          Filter::Screen::AntiAlias<W, H> aa;
          for (size_t i = 0; i < 3; ++i) {
            auto a = math::vector_to_pixel<W, H>(geometry(i, 0.0f));
            auto b = math::vector_to_pixel<W, H>(geometry(i, 1.0f));
            b.x = a.x + std::remainder(b.x - a.x, float(W));
            for (int branch = -1; branch <= 2; ++branch) {
              const float raw =
                  (a.x - branch * (W / 2.0f)) / (W / 2.0f - (b.x - a.x));
              if (raw < 0 && a.x < branch * (W / 2.0f) - 1)
                continue;
              if (raw > 1 && b.x > (branch + 1) * (W / 2.0f) + 1)
                continue;
              const float capture = hs::clamp(raw, 0.0f, 1.0f);
              const int center =
                  static_cast<int>(floorf(a.x + (b.x - a.x) * capture));
              for (int offset = 0; offset < 2; ++offset) {
                const int column = (center + offset + W) % W;
                const float t = float(column % (W / 2)) / (W / 2);
                if (fabsf(t - capture) > 0.5f || t < 0.5f)
                  continue;
                Fragment f;
                shader({}, f, t);
                aa.plot(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t,
                        f.color.color, 0.0f, f.color.alpha,
                        [&](float x, float y, const Pixel &color, float age,
                            float alpha) {
                          if (static_cast<int>(x) == column)
                            sink.plot(canvas, static_cast<int>(x),
                                      static_cast<int>(y), color, age, alpha);
                        });
              }
            }
          }
        }
      }
      effect.advance_display();
      std::vector<Pixel> pixels(W * H);
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          pixels[y * W + x] = effect.get_pixel(x, y);
      return pixels;
    };
    Temporal::geometry_counters = {};
    const auto actual = render(false);
    HS_EXPECT_EQ(Temporal::geometry_counters.dense_runs, 1u);
    HS_EXPECT_EQ(Temporal::geometry_counters.dense_gap_samples, 0u);
    const auto expected = render(true);
    int different = 0, lit = 0;
    for (size_t i = 0; i < actual.size(); ++i) {
      different += actual[i].r != expected[i].r ||
                   actual[i].g != expected[i].g || actual[i].b != expected[i].b;
      lit += expected[i].r > 0;
    }
    HS_EXPECT_GT(lit, 0);
    HS_EXPECT_EQ(different, 0);
  }
}

/** @brief Mixed dense and sparse runs own each open or closed junction once. */
inline void test_timed_dense_mixed_run_junctions() {
  constexpr int W = 288, H = 144;
  constexpr std::array<float, 6> COLUMNS{40.2f, 41.2f, 42.2f,
                                         52.2f, 53.2f, 54.2f};
  for (bool closed : {false, true}) {
    StubEffect effect(W, H);
    Canvas canvas(effect);
    Filter::Screen::DirectAntiAliasSink<W, H> sink;
    sink.prepare(canvas);
    ScratchScope scratch(plot_arena());
    Fragments points;
    points.bind(plot_arena(), COLUMNS.size());
    auto geometry = [&](size_t i, float t) {
      return math::pixel_to_vector<W, H>(COLUMNS[i] + 0.1f * t, 60.2f);
    };
    for (size_t i = 0; i < COLUMNS.size(); ++i) {
      Fragment f;
      f.pos = geometry(i, 1);
      f.v0 = float(i);
      points.push_back(f);
    }
    std::array<int, 6> endpoints{};
    auto shader = [&](const math::Vector &, Fragment &f, float time) {
      int i = int(roundf(f.v0));
      if (i >= 0 && i < 6 && fabsf(f.v0 - i) < 0.000001f &&
          fabsf(time * (W / 2) - roundf(time * (W / 2))) > 0.00001f)
        ++endpoints[i];
      f.color = Color4(255, 255, 255, 0.4f);
    };
    Temporal::geometry_counters = {};
    Plot::rasterize<W, H, {.timed = true}>(
        sink, canvas, points, shader,
        {.loop = closed ? Plot::RasterLoop::closed() : Plot::RasterLoop{}},
        Plot::GeometryCallback::interpolating(geometry),
        0.1f * math::RADIANS_PER_COLUMN<W>);
    HS_EXPECT_EQ(Temporal::geometry_counters.dense_runs, 2u);
    HS_EXPECT_EQ(Temporal::geometry_counters.dense_run_edges, 4u);
    HS_EXPECT_EQ(Temporal::geometry_counters.dense_sparse_edges,
                 closed ? 2u : 1u);
    for (int count : endpoints)
      HS_EXPECT_EQ(count, 1);
  }
}

/** @brief Insufficient dense cache space falls back to identical edge rendering. */
inline void test_timed_dense_cache_falls_back_to_edges() {
  constexpr int W = 288, H = 144;
  auto geometry = [](size_t i, float t) {
    return math::pixel_to_vector<W, H>(40.2f + 1.5f * i + 0.1f * t, 60.2f);
  };
  Plot::g_force_timed_edge_walk = true;
  const auto expected = seam_render<W, H>(
      geometry, 5, 0.1f * math::RADIANS_PER_COLUMN<W>, true, 1, 0, false);
  Plot::g_force_timed_edge_walk = false;
  std::vector<float> actual;
  {
    ScratchScope scratch(scratch_arena_b);
    const size_t available =
        scratch_arena_b.get_capacity() - scratch_arena_b.get_offset();
    scratch_arena_b.allocate_n<uint8_t>(available - 64);
    Temporal::geometry_counters = {};
    actual = seam_render<W, H>(geometry, 5, 0.1f * math::RADIANS_PER_COLUMN<W>,
                               true, 1, 0, false);
    HS_EXPECT_EQ(Temporal::geometry_counters.dense_reject_scratch, 1u);
    HS_EXPECT_EQ(Temporal::geometry_counters.dense_runs, 0u);
  }
  HS_EXPECT(actual == expected,
            "dense cache fallback preserves the edge frame");
}

/** @brief Dense gaps crossing an arm boundary retain the first column's coverage. */
inline void test_timed_dense_crossing_gap_keeps_boundary_coverage() {
  constexpr int W = 288, H = 144;
  for (float boundary : {0.0f, 144.0f}) {
    auto geometry = [boundary](size_t i, float t) {
      return math::pixel_to_vector<W, H>(boundary - 1.2f + 1.8f * i + 0.01f * t,
                                         60.2f);
    };
    const auto actual = seam_render<W, H>(
        geometry, 3, 0.01f * math::RADIANS_PER_COLUMN<W>, true, 1, 0, false);
    const auto frozen = seam_render<W, H>(
        geometry, 3, 0.01f * math::RADIANS_PER_COLUMN<W>, false, 0, 0, false);
    float observed = 0, reference = 0;
    const int column = int(boundary);
    for (int y = 58; y < 64; ++y) {
      observed += actual[y * W + column];
      reference += frozen[y * W + column];
    }
    HS_EXPECT_GT(reference, 0.25f);
    HS_EXPECT_GE(observed, reference * 0.8f);
  }
}

/** @brief Balanced density policies retain their edge sampler for dense controls. */
inline void test_timed_dense_preserves_balanced_policy() {
  constexpr int W = 288, H = 144;
  auto check = []<Plot::RasterSamplingPolicy Policy>() {
    StubEffect effect(W, H);
    Canvas canvas(effect);
    Filter::Screen::DirectAntiAliasSink<W, H> sink;
    sink.prepare(canvas);
    ScratchScope scratch(plot_arena());
    Fragments points;
    points.bind(plot_arena(), 3);
    auto geometry = [](size_t i, float t) {
      return math::pixel_to_vector<W, H>(40.2f + i + 0.1f * t, 60.2f);
    };
    for (size_t i = 0; i < 3; ++i) {
      Fragment f;
      f.pos = geometry(i, 1);
      points.push_back(f);
    }
    Temporal::geometry_counters = {};
    Plot::rasterize<
        W, H, {.single_pass = true, .sampling_policy = Policy, .timed = true}>(
        sink, canvas, points, white_shader(.4f), {.balanced_sampling = true},
        Plot::GeometryCallback::interpolating(geometry),
        0.1f * math::RADIANS_PER_COLUMN<W>);
    HS_EXPECT_EQ(Temporal::geometry_counters.dense_runs, 0u);
    HS_EXPECT_GT(Temporal::geometry_counters.rigid_samples, 0u);
  };
  check.template operator()<Plot::RasterSamplingPolicy::BALANCED>();
  check.template operator()<Plot::RasterSamplingPolicy::SELECTABLE>();
}

/** @brief Dense multi-turn paths wrap emitted centers and boundary taps. */
inline void test_timed_dense_wraps_multiple_windings() {
  constexpr int W = 96, H = 20;
  for (float direction : {-1.0f, 1.0f}) {
    auto geometry = [direction](size_t i, float t) {
      return math::pixel_to_vector<W, H>(
          40.2f + direction * 0.9f * i + 0.01f * t, 10.2f);
    };
    Temporal::geometry_counters = {};
    const auto pixels = seam_render<W, H>(
        geometry, 300, 0.01f * math::RADIANS_PER_COLUMN<W>, true, 1, 0, false);
    HS_EXPECT_EQ(Temporal::geometry_counters.dense_controls, 300u);
    HS_EXPECT_GT(Temporal::geometry_counters.boundary_evaluations, 0u);
    float energy = 0;
    for (float value : pixels)
      energy += value;
    HS_EXPECT_GT(energy, 1.0f);
  }
}

/** @brief Wider dense edges require a flat projected midpoint at both frame endpoints. */
inline void test_timed_wider_dense_edges_reject_curvature() {
  constexpr int W = 288, H = 144;
  for (bool curved : {false, true}) {
    auto geometry = [curved](size_t i, float t) {
      return math::pixel_to_vector<W, H>(40.2f + (curved ? 2.1f : 2.3f) * i +
                                             0.01f * t,
                                         curved ? .3f + 1.0f * i : 60.2f);
    };
    Temporal::geometry_counters = {};
    const auto pixels = seam_render<W, H>(
        geometry, 3, 0.01f * math::RADIANS_PER_COLUMN<W>, true, 1, 0, false);
    (void)pixels;
    HS_EXPECT_GT(Temporal::geometry_counters.dense_midpoint_probes, 0u);
    HS_EXPECT_EQ(Temporal::geometry_counters.dense_runs, curved ? 0u : 1u);
    HS_EXPECT_EQ(Temporal::geometry_counters.dense_wide_edges,
                 curved ? 0u : 2u);
    HS_EXPECT_EQ(Temporal::geometry_counters.dense_curved_edges,
                 curved ? 2u : 0u);
  }
}

/** @brief Wider edges retain the edge route when endpoint sample counts differ. */
inline void test_timed_wider_dense_edges_preserve_sample_count() {
  constexpr int W = 288, H = 144;
  auto geometry = [](size_t i, float t) {
    return math::pixel_to_vector<W, H>(40.2f + (1.7f + .6f * t) * i, 60.2f);
  };
  Plot::g_force_timed_edge_walk = true;
  const auto expected = seam_render<W, H>(
      geometry, 3, 1.2f * math::RADIANS_PER_COLUMN<W>, true, 1, 0, false);
  Plot::g_force_timed_edge_walk = false;
  Temporal::geometry_counters = {};
  const auto pixels = seam_render<W, H>(
      geometry, 3, 1.2f * math::RADIANS_PER_COLUMN<W>, true, 1, 0, false);
  HS_EXPECT(pixels == expected,
            "changing endpoint sample counts preserves the edge frame");
  HS_EXPECT_EQ(Temporal::geometry_counters.dense_runs, 0u);
  HS_EXPECT_EQ(Temporal::geometry_counters.dense_wide_edges, 0u);
  HS_EXPECT_GT(Temporal::geometry_counters.rigid_samples, 0u);
}
