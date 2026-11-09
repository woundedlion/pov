/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Feedback at column time: elapsed reads, the fade table, the analytic
// rotation oracle and arm-cut seams
// ============================================================================

/**
 * @brief Verifies FeedbackColumnTime's reads span their lit times.
 * @details A read within the arm spans one frame plus the column distance; a
 * read between an arm's last column and the next arm's first blends the two
 * taps' lit times; past the forward cut the source is displaced twice as far;
 * behind the backward cut the read slides from the previous arm's last column.
 */
inline void test_feedback_column_time_reads() {
  using Time = Filter::Pixel::FeedbackColumnTime<288>;
  constexpr float HALF = Time::HALF;
  HS_EXPECT_EQ(Time::arm_start(143), 0);
  HS_EXPECT_EQ(Time::arm_start(144), 144);
  HS_EXPECT_NEAR(Time::elapsed_scale(0.0f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(Time::elapsed_scale(16.0f), HALF / (HALF + 16.0f), 1e-6f);
  HS_EXPECT_NEAR(Time::elapsed_scale(-HALF), 2.0f, 1e-6f);

  // Column offsets relative to the arm wrap below its first column.
  HS_EXPECT_NEAR(Time::arm_column(150.0f, 144), 6.0f, 1e-6f);
  HS_EXPECT_NEAR(Time::arm_column(140.0f, 144), 284.0f, 1e-6f);
  HS_EXPECT_NEAR(Time::arm_column(-0.25f, 0), -0.25f, 1e-6f);

  HS_EXPECT_NEAR(Time::lit_time(10.0f), 10.0f / HALF, 1e-6f);
  HS_EXPECT_NEAR(Time::lit_time(HALF - 1.0f), (HALF - 1.0f) / HALF, 1e-6f);
  HS_EXPECT_NEAR(Time::lit_time(HALF + 3.0f), 3.0f / HALF, 1e-6f);

  auto check_tap_times = []<int W>() {
    using Model = Filter::Pixel::FeedbackColumnTime<W>;
    constexpr int H = W / 2;
    constexpr double U = std::numeric_limits<float>::epsilon() / 2.0;
    constexpr double GAMMA = 5 * U / (1 - 5 * U);
    for (const float column :
         {-1.0f, -0.75f, -0.5f, -0.25f, 0.0f, 0.25f, H - 1.0f, H - 0.75f,
          H - 0.5f, H - 0.25f, float(H), H + 0.25f, 2 * H - 1.0f,
          2 * H - 0.25f}) {
      const int lower = static_cast<int>(std::floor(double(column)));
      const int first = (lower % H + H) % H;
      const int second = (first + 1) % H;
      const double WEIGHT = double(column) - lower;
      const double EXPECTED = ((1 - WEIGHT) * first + WEIGHT * second) / H;
      const float actual = Model::lit_time(column);
      if constexpr (W == 256) {
        HS_EXPECT_EQ(actual, EXPECTED);
      } else {
        const double OFFSET = column < 0 ? H : column >= H ? -H : 0;
        const double IN_ARM = double(column) + OFFSET;
        const double TERMS = std::abs(double(column)) + std::abs(OFFSET);
        // Reciprocal, coefficient, products, two additions: gamma(5).
        const double SCALE =
            IN_ARM >= H - 1 ? (H + TERMS) * (H - 1) / H : TERMS / H;
        HS_EXPECT_NEAR(actual, EXPECTED, GAMMA * SCALE);
      }
    }
  };
  check_tap_times.template operator()<256>();
  check_tap_times.template operator()<288>();

  // Target column 20 of the second arm reading 7 columns ahead.
  const auto within = Time::read_within(164.0f, 144, 27.0f, 9.0f);
  HS_EXPECT_NEAR(within.x, 171.0f, 1e-5f);
  HS_EXPECT_NEAR(within.y, 9.0f, 1e-6f);
  HS_EXPECT_NEAR(within.elapsed, 1.0f - 7.0f / HALF, 1e-6f);

  // Past the forward cut: the twice-displaced source, lit a sweep earlier.
  int beyond_calls = 0;
  const auto forward =
      Time::read_past_cut(140.0f, 30.0f, 0, {147.0f, 31.0f}, [&]() {
        ++beyond_calls;
        return Time::Source{154.0f, 32.0f};
      });
  HS_EXPECT_EQ(beyond_calls, 1);
  HS_EXPECT_NEAR(forward.x, 154.0f, 1e-5f);
  HS_EXPECT_NEAR(forward.y, 32.0f, 1e-6f);
  HS_EXPECT_NEAR(forward.elapsed, 1.0f + 140.0f / HALF - 10.0f / HALF, 1e-5f);

  // Behind the backward cut: partway across the slide reads between the
  // previous arm's last column and the arm's first, their lit times blended.
  beyond_calls = 0;
  const auto gap =
      Time::read_past_cut(145.0f, 30.0f, 144, {141.5f, 31.0f}, [&]() {
        ++beyond_calls;
        return Time::Source{0.0f, 0.0f};
      });
  HS_EXPECT_EQ(beyond_calls, 0);
  const float slide = 1.5f / 3.5f;
  HS_EXPECT_NEAR(gap.x, 143.0f + slide, 1e-5f);
  HS_EXPECT_NEAR(gap.y, 30.0f + slide, 1e-5f);
  HS_EXPECT_NEAR(gap.elapsed,
                 1.0f + 1.0f / HALF - (1.0f - slide) * (HALF - 1.0f) / HALF,
                 1e-5f);
}

/** @brief Persistent fade tables rebuild for fade, hue and storage changes. */
inline void test_feedback_fade_table_invalidates_stored_entries() {
  using Table = Filter::Pixel::FeedbackFadeTable;
  alignas(std::max_align_t)
      std::array<uint8_t, 2 * Table::BYTES + 4 * alignof(float)>
          storage{};
  Arena arena(storage.data(), storage.size());
  ScratchScope scratch(scratch_arena_a);
  Table table;
  auto style = ::Feedback::Style::ArcingLightning();
  style.sync_hue();
  table.init_storage(arena);
  const auto first = table.acquire(style, scratch.get_arena());
  style.fade = 0.8f;
  style.sync_hue();
  const auto faded = table.acquire(style, scratch.get_arena());
  HS_EXPECT_EQ(faded.fades, first.fades);
  HS_EXPECT_EQ(faded.hue_k, first.hue_k);
  HS_EXPECT_EQ(faded.fades[Table::UNIT_LEVEL], style.fade);
  style.hue_shift += 0.25f;
  style.sync_hue();
  const auto turned = table.acquire(style, scratch.get_arena());
  HS_EXPECT_EQ(turned.hue_k, first.hue_k);
  const float SCALE = math::fast_cbrt(style.fade * Q16_TO_UNIT);
  for (int i = 0; i < Table::MATRIX; ++i)
    HS_EXPECT_EQ(turned.hue_k[Table::UNIT_LEVEL * Table::MATRIX + i],
                 style.hue_k[i] * SCALE);
  table.init_storage(arena);
  const auto rebound = table.acquire(style, scratch.get_arena());
  HS_EXPECT_NE(rebound.fades, first.fades);
  HS_EXPECT_EQ(rebound.fades[Table::UNIT_LEVEL], style.fade);
  for (int i = 0; i < Table::MATRIX; ++i)
    HS_EXPECT_EQ(rebound.hue_k[Table::UNIT_LEVEL * Table::MATRIX + i],
                 style.hue_k[i] * SCALE);
}

/**
 * @brief Verifies the fade table's entries apply their elapsed frames.
 * @details The one-frame entry is the per-frame transform bit for bit; the
 * two-frame entry fades twice and turns the hue twice as far; the zero-frame
 * entry neither fades nor turns.
 */
inline void test_feedback_fade_table_entries() {
  using Table = Filter::Pixel::FeedbackFadeTable;
  ::Feedback::Style style = ::Feedback::Style::ArcingLightning();
  style.sync_hue();
  ScratchScope scope(scratch_arena_a);
  Table table;
  const Table::View view = table.acquire(style, scope.get_arena());

  HS_EXPECT_EQ(Table::level(0.0f), 0);
  HS_EXPECT_EQ(Table::level(1.0f), Table::UNIT_LEVEL);
  HS_EXPECT_EQ(Table::level(2.0f), Table::LEVELS - 1);
  HS_EXPECT_EQ(Table::level(1.0f + 0.4f / Table::STEPS_PER_FRAME),
               Table::UNIT_LEVEL);

  HS_EXPECT_EQ(view.fades[Table::UNIT_LEVEL], style.fade);
  const float scale = math::fast_cbrt(style.fade * Q16_TO_UNIT);
  for (int i = 0; i < Table::MATRIX; ++i)
    HS_EXPECT_EQ(view.hue_k[Table::UNIT_LEVEL * Table::MATRIX + i],
                 style.hue_k[i] * scale);
  HS_EXPECT_NEAR(view.fades[0], 1.0f, 1e-6f);
  HS_EXPECT_NEAR(view.fades[Table::LEVELS - 1], style.fade * style.fade, 1e-5f);

  auto apply = [&](int level, const Pixel &p) {
    return ::Feedback::hue_fade_apply(view.hue_k + level * Table::MATRIX, p.r,
                                      p.g, p.b);
  };
  const Pixel source(40000, 20000, 10000);
  const Pixel none = apply(0, source);
  HS_EXPECT_NEAR(none.r, 40000.0f, 80.0f);
  HS_EXPECT_NEAR(none.g, 20000.0f, 80.0f);
  HS_EXPECT_NEAR(none.b, 10000.0f, 80.0f);
  const Pixel twice = apply(Table::LEVELS - 1, source);
  const Pixel two_steps =
      apply(Table::UNIT_LEVEL, apply(Table::UNIT_LEVEL, source));
  HS_EXPECT_NEAR(twice.r, two_steps.r, 120.0f);
  HS_EXPECT_NEAR(twice.g, two_steps.g, 120.0f);
  HS_EXPECT_NEAR(twice.b, two_steps.b, 120.0f);
}

/**
 * @brief Verifies the cached warp field is keyed by column time.
 * @details Offsets carry their elapsed scale only with column time, so a
 * field filled with it on must not serve a frame with it off.
 */
inline void test_feedback_warp_cache_keys_column_time() {
  constexpr int W = 64, H = 64;
  hs_test::StubEffect fx(W, H);
  ScratchScope persistent_scope(persistent_arena);
  Animation::NoiseParams noise;
  noise.noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  noise.set_seed(7);
  ::Feedback::Style style = ::Feedback::Style::LooseWormhole();
  style.noise = &noise;
  style.sync_noise();
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};
  auto &feedback = pipe.get<Filter::Pixel::Feedback<W, H>>();
  feedback.init_storage(persistent_arena);
  auto flush = [&]() {
    {
      Canvas c(fx);
      c(W / 2, H / 2) = Pixel(40000, 40000, 40000);
    }
    fx.advance_display();
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
    return FeedbackWhiteBox::holds(feedback, c);
  };
  {
    const ScopedTemporalSampling column_basis(true);
    HS_EXPECT_TRUE(flush());
  }
  fx.advance_display();
  const size_t fills = FeedbackWhiteBox::fills(feedback);
  {
    const ScopedTemporalSampling one_frame(false);
    HS_EXPECT_TRUE(flush());
  }
  fx.advance_display();
  HS_EXPECT_EQ(FeedbackWhiteBox::fills(feedback), fills + 1);
}

/** @brief Per-frame source rotation of rotation_feedback_warp(). */
inline math::Quaternion feedback_frame_rotation{};

/** @brief Space transform reading each direction's source one frame back
 *  under a rigid rotation.
 *  @param v Target direction.
 *  @return Source direction. */
inline math::Vector rotation_feedback_warp(const math::Vector &v,
                                           const ::Feedback::Style &) {
  return math::rotate(v, feedback_frame_rotation);
}

/**
 * @brief Content transported by a rigid rotation and a fade, at any time.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @details Content at direction v and time T came from rotate(v, q) at time
 * T - 1, q the per-frame rotation by `omega` about `axis`, faded by `fade`,
 * so F(T, v) = fade^T g(rotate(v, q^T)) for a smooth pattern g.
 */
template <int W, int H> struct RotatingFeedbackScene {
  math::Vector axis; ///< Rotation axis, unit length.
  float omega;       ///< Per-frame rotation angle, radians.
  float fade;        ///< Per-frame fade.

  /** @brief Content at a direction and time on the [0, 65535] scale.
   *  @param time Frames since the pattern was unrotated.
   *  @param v Unit direction.
   *  @return Channels r, g, b. */
  std::array<float, 3> at(float time, const math::Vector &v) const {
    const math::Vector p =
        math::rotate(v, math::make_rotation(axis, omega * time));
    const float scale = 65535.0f * powf(fade, time);
    return {scale * (0.5f + 0.35f * sinf(2.1f * p.x + 0.7f * p.z + 0.3f)),
            scale * (0.5f + 0.35f * sinf(1.7f * p.y - 1.3f * p.x + 1.1f)),
            scale * (0.5f + 0.35f * sinf(1.9f * p.z + 0.9f * p.y - 0.4f))};
  }

  /** @brief A frame as the display holds it: every column at its own lit
   *  time.
   *  @param frame Frame index.
   *  @param x Column.
   *  @param y Row.
   *  @return Channels r, g, b. */
  std::array<float, 3> column_basis(int frame, int x, int y) const {
    return at(static_cast<float>(frame) +
                  Spacetime::SamplingSurface<W>::time(x),
              math::pixel_to_vector<W, H>(x, y));
  }
};

/** @brief Error of a feedback run against the column-time oracle. */
struct ColumnTimeOracleError {
  float bulk; ///< Mean channel error away from the arm cuts.
  float cut;  ///< Mean channel error of the columns beside the cuts.
  float seam; ///< Perceived step across the cuts over the steps beside them.
};

/**
 * @brief Runs the feedback filter on an analytic rotation and compares every
 * column with the content at its lit time.
 * @param scene Rotation and fade.
 * @param frames Flushes to run after the column-time initial frame.
 * @param first_row First row compared.
 * @param end_row One past the last row compared.
 * @return Errors relative to full scale, and the perceived seam ratio: the
 *         step from a frame's last column to the next frame's first over the
 *         steps between the columns beside them.
 */
template <int W, int H>
ColumnTimeOracleError
run_column_time_oracle(const RotatingFeedbackScene<W, H> &scene, int frames,
                       int first_row, int end_row) {
  // A custom space transform composites uncached, from scratch.
  hs_test::ScopedArenaSplit split(1u << 20, 1u << 18, DEFAULT_SCRATCH_B_SIZE);
  hs_test::StubEffect fx(
      W, H, pipeline_config<Pipeline<W, H, Filter::Pixel::Feedback<W, H>>>());
  ::Feedback::Style style{};
  style.space_fn = &rotation_feedback_warp;
  style.fade = scene.fade;
  style.hue_shift = 0.0f;
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};
  feedback_frame_rotation = math::make_rotation(scene.axis, scene.omega);

  auto quantize = [](float c) {
    return static_cast<uint16_t>(hs::clamp(c, 0.0f, 65535.0f) + 0.5f);
  };
  // Reads behind a backward cut reach two frames back.
  for (int seed = -1; seed <= 0; ++seed) {
    {
      Canvas c(fx);
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
          const auto rgb = scene.column_basis(seed, x, y);
          c(x, y) = Pixel(quantize(rgb[0]), quantize(rgb[1]), quantize(rgb[2]));
        }
    }
    fx.advance_display();
  }
  std::vector<Pixel> previous(W * H);
  for (int frame = 1; frame <= frames; ++frame) {
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        previous[y * W + x] = fx.get_pixel(x, y);
    {
      Canvas c(fx);
      (void)pipe.begin_split_frame(c, 1.0f);
      pipe.finish_frame(c);
    }
    fx.advance_display();
  }

  constexpr int HALF = W / 2;
  constexpr int NEAR_CUT = 6;
  double bulk = 0.0, cut = 0.0;
  int bulk_count = 0, cut_count = 0;
  for (int x = 0; x < W; ++x) {
    double error = 0.0;
    for (int y = first_row; y < end_row; ++y) {
      const Pixel p = fx.get_pixel(x, y);
      const auto want = scene.column_basis(frames, x, y);
      error +=
          (fabs(p.r - want[0]) + fabs(p.g - want[1]) + fabs(p.b - want[2])) /
          (3.0 * 65535.0);
    }
    error /= end_row - first_row;
    if (std::min(x % HALF, HALF - x % HALF) < NEAR_CUT) {
      cut += error;
      ++cut_count;
    } else {
      bulk += error;
      ++bulk_count;
    }
  }

  auto luminance = [](const Pixel &p) {
    return (static_cast<double>(p.r) + p.g + p.b) / (3.0 * 65535.0);
  };
  double across = 0.0, beside = 0.0;
  for (int y = first_row; y < end_row; ++y)
    accumulate_arm_cut_seam_steps<W>(
        [&](int x) { return luminance(previous[y * W + x]); },
        [&](int x) { return luminance(fx.get_pixel(x, y)); }, NEAR_CUT - 1,
        across, beside);
  return {static_cast<float>(bulk / bulk_count),
          static_cast<float>(cut / cut_count),
          static_cast<float>(across / beside)};
}

/**
 * @brief A rotation about the vertical axis whose frames each plot a fixed
 *        pattern over the feedback at opacity `alpha`.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @details Every read within an arm spans the same `a` frames, so content k
 * plots back is `k a` frames old:
 * F(v) = sum_k alpha (1 - alpha)^k fade^(k a) g(rotate(v, q^(k a))).
 */
template <int W, int H> struct SourcedFeedbackScene {
  float omega; ///< Per-frame rotation about the vertical axis, radians.
  float fade;  ///< Per-frame fade.
  float alpha; ///< Opacity of each frame's plot.

  /** @brief The plotted pattern on the [0, 65535] scale.
   *  @param p Unit direction.
   *  @return Channels r, g, b. */
  static std::array<float, 3> pattern(const math::Vector &p) {
    return {65535.0f * (0.5f + 0.35f * sinf(2.1f * p.x + 0.7f * p.z + 0.3f)),
            65535.0f * (0.5f + 0.35f * sinf(1.7f * p.y - 1.3f * p.x + 1.1f)),
            65535.0f * (0.5f + 0.35f * sinf(1.9f * p.z + 0.9f * p.y - 0.4f))};
  }

  /** @brief Frames one read within an arm spans.
   *  @return HALF / (HALF + the source's column offset). */
  float step_frames() const {
    constexpr float HALF = W / 2.0f;
    const math::Vector v = math::pixel_to_vector<W, H>(W / 4, H / 2);
    const math::Vector s =
        math::rotate(v, math::make_rotation(math::Y_AXIS, omega));
    float offset = math::Spherical(s).theta * W / (2.0f * math::PI_F) - W / 4;
    offset -= W * std::round(offset / W);
    return HALF / (HALF + offset);
  }

  /** @brief Steady-state content at a pixel.
   *  @param x Column.
   *  @param y Row.
   *  @return Channels r, g, b. */
  std::array<float, 3> truth(int x, int y) const {
    const math::Vector v = math::pixel_to_vector<W, H>(x, y);
    const float a = step_frames();
    std::array<float, 3> sum{};
    float weight = alpha;
    for (int k = 0; k < 400 && weight > 1e-7f; ++k) {
      const float age = k * a;
      const auto g = pattern(
          math::rotate(v, math::make_rotation(math::Y_AXIS, omega * age)));
      const float w = weight * powf(fade, age);
      for (int c = 0; c < 3; ++c)
        sum[c] += w * g[c];
      weight *= 1.0f - alpha;
    }
    return sum;
  }
};

/** @brief Errors of a sourced feedback run against its steady state. */
struct SourcedOracleError {
  float bulk; ///< Mean channel error away from the arm cuts.
  float cut;  ///< Mean channel error of the columns beside the cuts.
  float seam; ///< Step across the cuts over the steps beside them.
  float ramp; ///< Brightness of each arm's first twelve columns over its last
              ///< twelve, each relative to the steady state.
};

/**
 * @brief Runs the feedback filter with a per-frame pattern plot until steady
 * and compares every column with the steady state.
 * @param scene Rotation, fade and plot opacity.
 * @param frames Frames to run.
 * @param first_row First row compared.
 * @param end_row One past the last row compared.
 * @return Errors relative to full scale, the seam ratio and the ramp.
 */
template <int W, int H>
SourcedOracleError run_sourced_oracle(const SourcedFeedbackScene<W, H> &scene,
                                      int frames, int first_row, int end_row) {
  hs_test::ScopedArenaSplit split(1u << 20, 1u << 18, DEFAULT_SCRATCH_B_SIZE);
  hs_test::StubEffect fx(
      W, H, pipeline_config<Pipeline<W, H, Filter::Pixel::Feedback<W, H>>>());
  ::Feedback::Style style{};
  style.space_fn = &rotation_feedback_warp;
  style.fade = scene.fade;
  style.hue_shift = 0.0f;
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};
  feedback_frame_rotation = math::make_rotation(math::Y_AXIS, scene.omega);

  auto quantize = [](float c) {
    return static_cast<uint16_t>(hs::clamp(c, 0.0f, 65535.0f) + 0.5f);
  };
  std::vector<std::array<float, 3>> plot(W * H);
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      plot[y * W + x] = SourcedFeedbackScene<W, H>::pattern(
          math::pixel_to_vector<W, H>(x, y));
  std::vector<Pixel> previous(W * H);
  for (int frame = 0; frame < frames; ++frame) {
    for (int i = 0; i < W * H; ++i)
      previous[i] = fx.get_pixel(i % W, i / W);
    {
      Canvas c(fx);
      auto &plotter = pipe.begin_split_frame(c, 1.0f);
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
          const auto &g = plot[y * W + x];
          plotter.plot(c, x, y,
                       Pixel(quantize(g[0]), quantize(g[1]), quantize(g[2])),
                       0.0f, scene.alpha);
        }
      pipe.finish_frame(c);
    }
    fx.advance_display();
  }

  constexpr int HALF = W / 2;
  constexpr int NEAR_CUT = 6;
  constexpr int RAMP_COLUMNS = 12;
  double bulk = 0.0, cut = 0.0, start = 0.0, end = 0.0;
  int bulk_count = 0, cut_count = 0;
  for (int x = 0; x < W; ++x) {
    double error = 0.0, shown = 0.0, wanted = 0.0;
    for (int y = first_row; y < end_row; ++y) {
      const Pixel p = fx.get_pixel(x, y);
      const auto want = scene.truth(x, y);
      error +=
          (fabs(p.r - want[0]) + fabs(p.g - want[1]) + fabs(p.b - want[2])) /
          (3.0 * 65535.0);
      shown += double(p.r) + p.g + p.b;
      wanted += double(want[0]) + want[1] + want[2];
    }
    error /= end_row - first_row;
    const int column = x % HALF;
    if (column < RAMP_COLUMNS)
      start += shown / wanted;
    else if (column >= HALF - RAMP_COLUMNS)
      end += shown / wanted;
    if (std::min(column, HALF - column) < NEAR_CUT) {
      cut += error;
      ++cut_count;
    } else {
      bulk += error;
      ++bulk_count;
    }
  }
  auto luminance = [](const Pixel &p) {
    return (static_cast<double>(p.r) + p.g + p.b) / (3.0 * 65535.0);
  };
  double across = 0.0, beside = 0.0;
  for (int y = first_row; y < end_row; ++y)
    accumulate_arm_cut_seam_steps<W>(
        [&](int x) { return luminance(previous[y * W + x]); },
        [&](int x) { return luminance(fx.get_pixel(x, y)); }, NEAR_CUT - 1,
        across, beside);
  return {static_cast<float>(bulk / bulk_count),
          static_cast<float>(cut / cut_count),
          static_cast<float>(across / beside), static_cast<float>(start / end)};
}

/**
 * @brief Split-frame feedback with a per-frame plot under azimuthal rotations
 * matches its steady state, with no brightness ramp across the arms and no
 * step at the cuts.
 */
inline void test_feedback_sourced_column_time_oracle() {
  constexpr int W = 288, H = 144;
  const ScopedTemporalSampling column_basis(true);
  struct Case {
    const char *name;
    float omega;
  };
  const Case cases[] = {{"forward", 0.1f},
                        {"backward", -0.1f},
                        {"forward fast", 0.3f},
                        {"backward fast", -0.3f}};
  for (const Case &c : cases) {
    const SourcedFeedbackScene<W, H> scene{c.omega, 0.85f, 0.3f};
    const auto error = run_sourced_oracle<W, H>(scene, 40, 24, H - 24);
    std::printf("  sourced oracle %s: bulk=%.5f cut=%.5f seam=%.3f ramp=%.4f\n",
                c.name, error.bulk, error.cut, error.seam, error.ramp);
    HS_EXPECT_LT(error.seam, 1.2f);
    HS_EXPECT_LT(fabsf(error.ramp - 1.0f), 0.01f);
    HS_EXPECT_LE(error.cut, 1.5f * error.bulk);
  }
}

/**
 * @brief Feedback under rigid rotations matches the content at each column's
 * lit time, and steps no further across the arm cuts than beside them.
 * @details Azimuthal rotations carry content across both cuts in both
 * directions, through the reads past the forward cut and the slide behind
 * the backward one; a rotation about a horizontal axis carries it through
 * the cap-plane rows; a tilted one through every row. Polar rows composite
 * column pairs, and the pairs beside a cut must not step. A one-frame warp holds
 * content at its previous frame's time, a sweep's worth of error that steps
 * at the cuts.
 */
inline void test_feedback_column_time_oracle() {
  constexpr int W = 288, H = 144;
  const ScopedTemporalSampling column_basis(true);
  struct Case {
    const char *name;
    math::Vector axis;
    float omega;
    int first_row;
    int end_row;
    float seam_bound; ///< Largest seam ratio; 0 leaves it unchecked.
  };
  const Case cases[] = {
      {"azimuthal forward", math::Y_AXIS, 0.1f, 24, H - 24, 1.5f},
      {"azimuthal backward", math::Y_AXIS, -0.1f, 24, H - 24, 1.5f},
      {"cap-plane rows", math::X_AXIS, 0.1f, 4, 24, 1.2f},
      {"north infill rows", math::Y_AXIS, 0.1f, 1, 6, 1.5f},
      {"south infill rows", math::Y_AXIS, -0.1f, H - 6, H - 1, 1.5f},
      {"half-resolution polar rows", math::Y_AXIS, 0.1f, 6, 24, 1.2f},
      {"tilted", math::Vector(0.6f, 0.0f, 0.8f), 0.08f, 4, H - 4, 0.0f},
  };
  for (const Case &c : cases) {
    const RotatingFeedbackScene<W, H> scene{c.axis, c.omega, 0.85f};
    const auto error =
        run_column_time_oracle<W, H>(scene, 6, c.first_row, c.end_row);
    std::printf("  column-time oracle %s: bulk=%.5f cut=%.5f seam=%.3f\n",
                c.name, error.bulk, error.cut, error.seam);
    HS_EXPECT_LT(error.bulk, 0.002f);
    HS_EXPECT_LT(error.cut, 0.0015f);
    if (c.seam_bound > 0.0f)
      HS_EXPECT_LT(error.seam, c.seam_bound);
  }
}

/**
 * @brief A split frame's forward-band coverage composites the band beneath its
 * taps as plotting them over the band does, within a few 16-bit rounding
 * steps.
 * @details Each draw runs once through the feedback pipeline over a split
 * frame, and once through the same stages without the terminal over the same
 * frame. Pixels outside the held columns match bit for bit. Draws: sub-pixel
 * points over every column, and each solid's wireframe through the column
 * walk under a moving orientation.
 */
inline void test_feedback_split_frame_covers_forward_bands() {
  constexpr int W = 288, H = 144;
  const ScopedTemporalSampling column_basis(true);
  hs_test::ScopedArenaSplit split(1u << 21, 1u << 19, DEFAULT_SCRATCH_B_SIZE);
  using Terminal = Filter::Pixel::Feedback<W, H>;
  using AntiAlias = Filter::Screen::AntiAlias<W, H>;
  hs_test::StubEffect fx(
      W, H, pipeline_config<Pipeline<W, H, AntiAlias, Terminal>>());
  ::Feedback::Style style{};
  style.space_fn = &rotation_feedback_warp;
  style.fade = 0.85f;
  style.hue_shift = 0.0f;
  feedback_frame_rotation = math::make_rotation(math::Y_AXIS, 0.15f);
  const RotatingFeedbackScene<W, H> scene{math::Y_AXIS, 0.15f, 0.85f};
  auto quantize = [](float c) {
    return static_cast<uint16_t>(hs::clamp(c, 0.0f, 65535.0f) + 0.5f);
  };
  auto seed_history = [&]() {
    for (int frame = -1; frame <= 0; ++frame) {
      {
        Canvas c(fx);
        for (int y = 0; y < H; ++y)
          for (int x = 0; x < W; ++x) {
            const auto rgb = scene.column_basis(frame, x, y);
            c(x, y) =
                Pixel(quantize(rgb[0]), quantize(rgb[1]), quantize(rgb[2]));
          }
      }
      fx.advance_display();
    }
  };
  auto snapshot = [&](Canvas &c) {
    std::vector<Pixel> pixels(W * H);
    for (int i = 0; i < W * H; ++i)
      pixels[i] = c(i % W, i / W);
    return pixels;
  };

  // Runs one split frame: `held_draw` plots through the frame, `plain_draw`
  // plots the same through the stages before the terminal.
  auto check = [&](const char *name, Terminal &terminal, auto &&begin,
                   auto &&held_draw, auto &&plain_draw) {
    seed_history();
    {
      Canvas c(fx);
      begin(c);
      std::vector<Pixel> under = snapshot(c);
      std::vector<bool> held(W * H);
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
          held[y * W + x] =
              x % (W / 2) >= FeedbackWhiteBox::held_from(terminal, y);
          if (held[y * W + x])
            under[y * W + x] = Pixel(static_cast<uint16_t>(9000 + 197 * x),
                                     static_cast<uint16_t>(60000 - 331 * y),
                                     static_cast<uint16_t>(30000 + 89 * x));
        }
      held_draw(c);
      const int covered = FeedbackWhiteBox::covered(terminal);
      FeedbackWhiteBox::release(terminal, under.data());
      const std::vector<Pixel> composited = snapshot(c);
      for (int i = 0; i < W * H; ++i)
        c(i % W, i / W) = under[i];
      plain_draw(c);
      int outside = 0, worst = 0;
      for (int i = 0; i < W * H; ++i) {
        const Pixel a = composited[i], b = c(i % W, i / W);
        if (!held[i]) {
          outside += !(a == b);
          continue;
        }
        worst = std::max({worst, std::abs(int(a.r) - int(b.r)),
                          std::abs(int(a.g) - int(b.g)),
                          std::abs(int(a.b) - int(b.b))});
      }
      std::printf("  split frame %s: covered=%d outside=%d worst=%d\n", name,
                  covered, outside, worst);
      HS_EXPECT_GT(covered, 0);
      HS_EXPECT_EQ(outside, 0);
      HS_EXPECT_LE(worst, 4);
    }
    fx.advance_display();
  };

  {
    Pipeline<W, H, AntiAlias, Terminal> pipe{AntiAlias{}, Terminal(style)};
    Pipeline<W, H, AntiAlias> plain{AntiAlias{}};
    auto points = [](auto &target, Canvas &c) {
      for (int i = 0; i < 6000; ++i) {
        const float x = fmodf(i * 0.6180339f * 7.0f, float(W));
        const float y = 1.0f + fmodf(i * 0.4142136f * 5.0f, float(H - 2));
        const Pixel color(static_cast<uint16_t>(i * 977 % 65536),
                          static_cast<uint16_t>(i * 331 % 65536), 40000);
        target.plot(c, x, y, color, 0.0f, 0.25f + 0.5f * ((i % 7) / 6.0f));
      }
    };
    auto &terminal = pipe.get<Terminal>();
    PreparedTerminalFrame<decltype(pipe)> *frame = nullptr;
    check(
        "points", terminal,
        [&](Canvas &c) { frame = &pipe.begin_split_frame(c, 1.0f); },
        [&](Canvas &c) { points(*frame, c); },
        [&](Canvas &c) { points(plain, c); });
  }

  math::Orientation<> orientation;
  const math::Vector axis = math::Vector(0.2f, 1.0f, 0.1f).normalized();
  const auto start_pose = math::make_rotation(math::X_AXIS, 0.27f);
  orientation.set(start_pose);
  orientation.push(math::make_rotation(axis, 0.05f) * start_pose);
  using Orient = Filter::World::Orient;
  Pipeline<W, H, Orient, AntiAlias, Terminal> pipe{
      Orient(orientation), AntiAlias{}, Terminal(style)};
  Pipeline<W, H, Orient, AntiAlias> plain{Orient(orientation), AntiAlias{}};
  ScratchScope persistent_scope(persistent_arena);
  pipe.init_storage(persistent_arena);
  plain.init_storage(persistent_arena);
  Plot::GeometryCallback geometry;
  geometry.source_positions = false;
  const Pixel shade(52000, 30000, 12000);
  auto shader = [&](const math::Vector &, Fragment &f) {
    f.color = Color4(shade, 1.0f);
  };
  auto &terminal = pipe.get<Terminal>();
  for (size_t solid = 0; solid < Solids::BASE_MESH_COUNT; ++solid) {
    ScratchScope mesh_scope(persistent_arena);
    MeshState mesh;
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      MeshOps::compile(Solids::get_entry(solid).generate(a, b), mesh, target,
                       a);
    });
    ArenaVector<Plot::Mesh::Edge> edges;
    edges.bind(persistent_arena, Solids::MAX_SOLID_EDGES);
    Plot::Mesh::extract_edges(mesh, edges);
    PreparedTerminalFrame<decltype(pipe)> *frame = nullptr;
    check(
        Solids::get_entry(solid).name, terminal,
        [&](Canvas &c) { frame = &pipe.begin_split_frame(c, 1.0f); },
        [&](Canvas &c) {
          Plot::Mesh::draw_timed<W, H>(*frame, c, mesh, edges, shader,
                                       geometry);
        },
        [&](Canvas &c) {
          Plot::Mesh::draw_timed<W, H>(plain, c, mesh, edges, shader, geometry);
        });
  }
}
