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
        HS_EXPECT_NEAR(actual, EXPECTED, 1e-6);
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
  hs_test::StubEffect fx(W, H);
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
  {
    Canvas c(fx);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        const auto rgb = scene.column_basis(0, x, y);
        c(x, y) = Pixel(quantize(rgb[0]), quantize(rgb[1]), quantize(rgb[2]));
      }
  }
  fx.advance_display();
  std::vector<Pixel> previous(W * H);
  for (int frame = 1; frame <= frames; ++frame) {
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        previous[y * W + x] = fx.get_pixel(x, y);
    {
      Canvas c(fx);
      (void)pipe.begin_frame(c, 1.0f);
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
    for (int start : {0, HALF}) {
      const int last = (start + W - 1) % W;
      across += fabs(luminance(fx.get_pixel(start, y)) -
                     luminance(previous[y * W + last]));
      for (int k = 1; k < NEAR_CUT; ++k) {
        beside += fabs(luminance(fx.get_pixel(start + k, y)) -
                       luminance(fx.get_pixel(start + k - 1, y))) /
                  (2 * (NEAR_CUT - 1));
        beside += fabs(luminance(previous[y * W + last - k]) -
                       luminance(previous[y * W + last - k + 1])) /
                  (2 * (NEAR_CUT - 1));
      }
    }
  return {static_cast<float>(bulk / bulk_count),
          static_cast<float>(cut / cut_count),
          static_cast<float>(across / beside)};
}

/**
 * @brief Feedback under rigid rotations matches the content at each column's
 * lit time, and steps no further across the arm cuts than beside them.
 * @details Azimuthal rotations carry content across both cuts in both
 * directions, through the reads past the forward cut and the slide behind
 * the backward one; a rotation about a horizontal axis carries it through
 * the cap-plane rows; a tilted one through every row. A one-frame warp holds
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
    bool seam;
  };
  const Case cases[] = {
      {"azimuthal forward", math::Y_AXIS, 0.1f, 24, H - 24, true},
      {"azimuthal backward", math::Y_AXIS, -0.1f, 24, H - 24, true},
      {"cap-plane rows", math::X_AXIS, 0.1f, 4, 24, false},
      {"tilted", math::Vector(0.6f, 0.0f, 0.8f), 0.08f, 4, H - 4, false},
  };
  for (const Case &c : cases) {
    const RotatingFeedbackScene<W, H> scene{c.axis, c.omega, 0.85f};
    const auto error =
        run_column_time_oracle<W, H>(scene, 6, c.first_row, c.end_row);
    std::printf("  column-time oracle %s: bulk=%.5f cut=%.5f seam=%.3f\n",
                c.name, error.bulk, error.cut, error.seam);
    HS_EXPECT_LT(error.bulk, 0.002f);
    HS_EXPECT_LT(error.cut, 0.0015f);
    if (c.seam)
      HS_EXPECT_LT(error.seam, 1.5f);
  }
}
