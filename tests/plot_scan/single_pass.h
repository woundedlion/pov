/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// SinglePass rasterizer / PlanarEdgeSampler::one_pass
// ============================================================================

/**
 * @brief Builds the planar sampler the rasterizer would build for edge a->b.
 * @param a Edge start (unit vector).
 * @param b Edge end (unit vector).
 * @param basis Azimuthal-equidistant projection basis.
 * @return The sampler, with its cumulative-arc table already filled.
 */
inline Plot::PlanarEdgeSampler planar_sampler(const math::Vector &a,
                                              const math::Vector &b,
                                              const math::Basis &basis) {
  Fragment fa, fb;
  fa.pos = a;
  fb.pos = b;
  Plot::PlanarEdgeSampler out;
  Plot::rasterize_planar_strategy(fa, fb, basis, /*is_last_segment=*/true,
                                  [&](const Plot::PlanarEdgeSampler &s,
                                      const Fragment &, const Fragment &, float,
                                      bool) { out = s; });
  return out;
}

/**
 * @brief one_pass's analytic tangent agrees with a forward difference of pos(),
 *        at the same position.
 * @details projection_fraction feeds both sides: a reversed or plateaued
 *          mapping shows here, a monotone one that is not arc-uniform does not.
 */
inline void test_planar_one_pass_matches_forward_difference() {
  hs::random().seed(0x51F1);
  int checked = 0;
  float worst_len = 0.0f, worst_pos = 0.0f, worst_tan_len = 0.0f;
  float worst_tan_dot = 1.0f;
  for (int trial = 0; trial < 400; ++trial) {
    math::Basis basis = basis_from_normal(rand_unit());
    math::Vector a = rand_unit();
    math::Vector b = rand_unit();
    // Antipodal-seam edges do not take the planar strategy.
    if (math::dot(a, basis.v) < -Plot::COS_PLANAR_ANTIPODE ||
        math::dot(b, basis.v) < -Plot::COS_PLANAR_ANTIPODE)
      continue;
    Plot::PlanarEdgeSampler s = planar_sampler(a, b, basis);
    if (s.dist < 0.05f)
      continue;

    for (int k = 0; k <= 8; ++k) {
      const float t = static_cast<float>(k) / 8.0f;
      Plot::SamplePT one = s.one_pass(t);
      // Short enough to track the arc, long enough to clear float cancellation.
      constexpr float DT = 1.0f / 256.0f;
      const bool fwd = (t + DT <= 1.0f);
      const math::Vector fd_pos = s.pos(t);
      const math::Vector step = s.pos(fwd ? t + DT : t - DT);
      const math::Vector fd_tan =
          (fwd ? (step - fd_pos) : (fd_pos - step)).normalized();
      worst_len = fold_worst(worst_len, std::abs(one.pos.length() - 1.0f));
      worst_pos = fold_worst(worst_pos, math::angle_between(one.pos, fd_pos));
      worst_tan_len =
          fold_worst(worst_tan_len, std::abs(one.tan.length() - 1.0f));
      worst_tan_dot = -fold_worst(-worst_tan_dot, -math::dot(one.tan, fd_tan));
      ++checked;
    }
  }
  HS_EXPECT_GT(checked, 1000);
  // Bounds are the fast_sinf/fast_cosf budget the two derivations each pay.
  HS_EXPECT_LE(worst_len, 5e-3f);
  HS_EXPECT_LE(worst_pos, 2e-2f);
  HS_EXPECT_LE(worst_tan_len, 5e-3f);
  // Same direction, not merely the same line.
  HS_EXPECT_GT(worst_tan_dot, 0.9998f);
}

/**
 * @brief one_pass's tangent is tangent to the sphere and points along the
 *        rendered chart line.
 * @details Independent of pos(): the tangent must be orthogonal to the
 *          position and must advance toward the edge's far endpoint.
 */
inline void test_planar_one_pass_tangent_is_forward_and_orthogonal() {
  hs::random().seed(0x51F2);
  int checked = 0;
  float worst_orth = 0.0f;
  for (int trial = 0; trial < 400; ++trial) {
    math::Basis basis = basis_from_normal(rand_unit());
    math::Vector a = rand_unit();
    math::Vector b = rand_unit();
    if (math::dot(a, basis.v) < -Plot::COS_PLANAR_ANTIPODE ||
        math::dot(b, basis.v) < -Plot::COS_PLANAR_ANTIPODE)
      continue;
    Plot::PlanarEdgeSampler s = planar_sampler(a, b, basis);
    if (s.dist < 0.2f)
      continue;

    for (int k = 0; k <= 4; ++k) {
      const float t = static_cast<float>(k) / 4.0f;
      Plot::SamplePT one = s.one_pass(t);
      worst_orth =
          fold_worst(worst_orth, std::abs(math::dot(one.pos, one.tan)));
      const float step = 1.0f / 64.0f;
      const math::Vector ahead = s.pos(std::min(1.0f, t + step));
      const math::Vector behind = s.pos(std::max(0.0f, t - step));
      HS_EXPECT_GT(math::dot(one.tan, ahead - behind), 0.0f);
      ++checked;
    }
  }
  HS_EXPECT_GT(checked, 500);
  HS_EXPECT_LE(worst_orth, 2e-2f);
}

/**
 * @brief The SinglePass rasterizer draws the same gap-free planar edge as the
 *        cached two-pass path.
 */
inline void test_rasterize_single_pass_planar_matches_two_pass() {
  constexpr int W = 128, H = 64;
  constexpr float base_step = (2.0f * math::PI_F) / W;
  hs_test::StubEffect fx(W, H);
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 4);

  math::Basis basis = basis_from_normal(math::Vector(0, 1, 0));
  Fragment a, b;
  a.pos = disk_point(basis, 0.3f, 0.0f);
  b.pos = disk_point(basis, 1.3f, 1.0f);
  points.push_back(a);
  points.push_back(b);

  CapturePipeline single, cached;
  {
    Canvas c(fx);
    Plot::rasterize<W, H, Plot::RasterConfig{.single_pass = true}>(
        single, c, points, noop_shader,
        {.projection = Plot::RasterProjection::planar(basis)});
  }
  fx.advance_display();
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(
        cached, c, points, noop_shader,
        {.projection = Plot::RasterProjection::planar(basis)});
  }
  fx.advance_display();

  HS_EXPECT_GT(single.plotted.size(), (size_t)10);
  HS_EXPECT_LE(max_consecutive_gap(single.plotted, /*wrap=*/false),
               1.5f * base_step);
  for (const math::Vector &p : single.plotted)
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
  if (single.plotted.empty())
    return;
  HS_EXPECT_NEAR(math::angle_between(single.plotted.front(), a.pos), 0.0f,
                 1e-2f);
  HS_EXPECT_NEAR(math::angle_between(single.plotted.back(), b.pos), 0.0f,
                 1e-2f);

  // Same curve: every sample lies within a dot of the other path's.
  HS_EXPECT_LE(single.plotted.size(), cached.plotted.size() + 2);
  HS_EXPECT_GE(single.plotted.size() + 2, cached.plotted.size());
  for (const math::Vector &p : single.plotted) {
    float nearest = math::PI_F;
    for (const math::Vector &q : cached.plotted)
      nearest = std::min(nearest, math::angle_between(p, q));
    HS_EXPECT_LE(nearest, base_step);
  }
}

/**
 * @brief SinglePass closes a planar loop's seam the way the cached path does.
 */
inline void test_rasterize_single_pass_closed_loop_matches_two_pass() {
  constexpr int W = 128, H = 64;
  constexpr float base_step = (2.0f * math::PI_F) / W;
  hs_test::StubEffect fx(W, H);
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 8);

  math::Basis basis = basis_from_normal(math::Vector(0, 1, 0));
  for (int i = 0; i < 5; ++i) {
    float az = (2.0f * math::PI_F * i) / 5.0f;
    math::Vector dir = basis.u * cosf(az) + basis.w * sinf(az);
    Fragment f;
    f.pos = (basis.v * cosf(0.7f) + dir * sinf(0.7f)).normalized();
    points.push_back(f);
  }

  const Plot::RasterOptions opts = {.loop = Plot::RasterLoop::closed(),
                                    .projection =
                                        Plot::RasterProjection::planar(basis)};
  CapturePipeline single, cached;
  {
    Canvas c(fx);
    Plot::rasterize<W, H, Plot::RasterConfig{.single_pass = true}>(
        single, c, points, noop_shader, opts);
  }
  fx.advance_display();
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(cached, c, points, noop_shader, opts);
  }
  fx.advance_display();

  HS_EXPECT_GT(single.plotted.size(), (size_t)20);
  HS_EXPECT_LE(max_consecutive_gap(single.plotted, /*wrap=*/true),
               1.5f * base_step);
  for (const math::Vector &p : single.plotted)
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
  HS_EXPECT_LE(single.plotted.size(), cached.plotted.size() + points.size());
  HS_EXPECT_GE(single.plotted.size() + points.size(), cached.plotted.size());
  if (single.plotted.empty())
    return;
  HS_EXPECT_LE(
      math::angle_between(single.plotted.front(), single.plotted.back()),
      1.5f * base_step);
  const math::Vector seam_delta = single.plotted.back() - points[0].pos;
  HS_EXPECT_GT(math::dot(seam_delta, seam_delta), 1e-10f);
  for (const math::Vector &p : single.plotted) {
    float nearest = math::PI_F;
    for (const math::Vector &q : cached.plotted)
      nearest = std::min(nearest, math::angle_between(p, q));
    HS_EXPECT_LE(nearest, base_step);
  }
}

/** @brief Endpoint-aware single-pass steps match constant-speed replay. */
inline void test_rasterize_single_pass_balances_terminal_interval() {
  constexpr int W = 128, H = 64;
  constexpr float BASE_STEP = (2.0f * math::PI_F) / W;
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 2);

  math::Vector start(1.0f, 0.0f, 0.0f);
  math::Vector tangent(0.0f, 0.0f, 1.0f);
  float desired_step = Plot::screen_step<W, H>(start, tangent, BASE_STEP);
  float arc = desired_step * 1.01f;
  Fragment a, b;
  a.pos = start;
  a.v0 = 0.0f;
  b.pos = math::Vector(cosf(arc), 0.0f, sinf(arc));
  b.v0 = 1.0f;
  points.push_back(a);
  points.push_back(b);

  auto draw = [&](bool single_pass) {
    hs_test::StubEffect fx(W, H);
    CapturePipeline pipeline;
    Canvas canvas(fx);
    std::vector<float> samples;
    auto shader = [&](const math::Vector &, Fragment &f) {
      samples.push_back(f.v0);
    };
    if (single_pass)
      Plot::rasterize<W, H, Plot::RasterConfig{.single_pass = true}>(
          pipeline, canvas, points, shader);
    else
      Plot::rasterize<W, H>(pipeline, canvas, points, shader);
    return samples;
  };
  std::vector<float> single = draw(true);
  std::vector<float> cached = draw(false);

  HS_EXPECT_SIZE_OR_RETURN(single, size_t{3});
  HS_EXPECT_SIZE_OR_RETURN(single, cached.size());
  HS_EXPECT_NEAR(single[1], 0.5f, 1e-4f);
  for (size_t i = 0; i < single.size(); ++i)
    HS_EXPECT_NEAR(single[i], cached[i], 1e-4f);
}

/**
 * @brief Exhausting the sub-step budget coarsens a segment in both rasterizer
 *        paths, and truncates it in neither.
 */
inline void test_rasterize_step_budget_backstop_finishes_segment() {
  constexpr int W = 128, H = 64;
  constexpr float BASE_STEP = (2.0f * math::PI_F) / W;
  constexpr size_t BUDGET = 16;
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 2);

  // Equatorial quarter turn: well past the lowered budget.
  Fragment a, b;
  a.pos = math::Vector(1.0f, 0.0f, 0.0f);
  b.pos = math::Vector(0.0f, 0.0f, 1.0f);
  points.push_back(a);
  points.push_back(b);

  struct Capture {
    std::vector<math::Vector> plotted;
    uint32_t backstops = 0;
  };
  auto draw = [&](bool single_pass) {
    hs_test::StubEffect fx(W, H);
    CapturePipeline pipeline;
    hs::g_scan_metrics.reset();
    {
      Canvas canvas(fx);
      if (single_pass)
        Plot::rasterize<W, H, Plot::RasterConfig{.single_pass = true}>(
            pipeline, canvas, points, noop_shader);
      else
        Plot::rasterize<W, H>(pipeline, canvas, points, noop_shader);
    }
    fx.advance_display();
    return Capture{pipeline.plotted, hs::g_scan_metrics.plot_backstop_hits};
  };

  struct RestoreBudget {
    decltype(Plot::g_step_budget_override) saved = Plot::g_step_budget_override;
    ~RestoreBudget() { Plot::g_step_budget_override = saved; }
  } restore_budget;
  Plot::g_step_budget_override = BUDGET;
  const Capture single = draw(true);
  const Capture cached = draw(false);

  HS_EXPECT_EQ(single.backstops, uint32_t{1});
  HS_EXPECT_EQ(cached.backstops, uint32_t{1});

  // Both endpoints drawn, no hole in between, and the emitted count still bound
  // by the budget rather than by the unstretched cadence.
  HS_EXPECT_GT(single.plotted.size(), size_t{0});
  if (single.plotted.empty())
    return;
  HS_EXPECT_NEAR(math::angle_between(single.plotted.front(), a.pos), 0.0f,
                 1e-3f);
  HS_EXPECT_NEAR(math::angle_between(single.plotted.back(), b.pos), 0.0f,
                 1e-3f);
  const float single_gap = max_consecutive_gap(single.plotted, /*wrap=*/false);
  const float cached_gap = max_consecutive_gap(cached.plotted, /*wrap=*/false);
  HS_EXPECT_LE(single_gap, 2.5f * BASE_STEP);
  HS_EXPECT_LE(single_gap, 1.25f * cached_gap);
  HS_EXPECT_LE(single.plotted.size(), 2 * BUDGET);
  HS_EXPECT_GT(single.plotted.size(), BUDGET);
}

/** @brief SELECTABLE with balanced sampling off matches DEFAULT within floating-point tolerance. */
inline void test_rasterize_default_sampling_policy_parity() {
  static_assert(
      std::is_same_v<
          std::integral_constant<Plot::RasterConfig,
                                 Plot::RasterConfig{.single_pass = true}>,
          std::integral_constant<
              Plot::RasterConfig,
              Plot::RasterConfig{.single_pass = true,
                                 .sampling_policy =
                                     Plot::RasterSamplingPolicy::DEFAULT}>>);
  constexpr int W = 128, H = 64;
  constexpr float POLICY_TOL = 1e-4f;
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 18);
  const math::Basis shape_basis = math::make_basis(
      math::Quaternion(0.91f, -0.17f, 0.31f, 0.21f).normalized(), math::X_AXIS);
  constexpr float RADIUS = 0.83f;
  Plot::Star<Plot::PlanarProjection>::sample(points, shape_basis, RADIUS, 8,
                                             0.73f);
  const math::Basis planar_basis =
      Plot::planar_chart_basis(math::get_antipode(shape_basis, RADIUS).first.v);

  auto capture = [&]() {
    hs_test::StubEffect fx(W, H);
    AlphaCapturePipeline pipeline;
    Canvas canvas(fx);
    auto shader = [](const math::Vector &, Fragment &f) {
      f.color = Color4(Pixel(10000, 20000, 30000), 0.37f);
    };
    Plot::rasterize<W, H, Plot::RasterConfig{.single_pass = true}>(
        pipeline, canvas, points, shader,
        {.projection = Plot::RasterProjection::planar(planar_basis),
         .omit_end = true});
    return pipeline;
  };

  const AlphaCapturePipeline implicit_default = capture();
  hs_test::StubEffect selectable_fx(W, H);
  AlphaCapturePipeline selectable_default;
  {
    Canvas canvas(selectable_fx);
    auto shader = [](const math::Vector &, Fragment &f) {
      f.color = Color4(Pixel(10000, 20000, 30000), 0.37f);
    };
    Plot::rasterize<W, H,
                    Plot::RasterConfig{
                        .single_pass = true,
                        .sampling_policy =
                            Plot::RasterSamplingPolicy::SELECTABLE}>(
        selectable_default, canvas, points, shader,
        {.projection = Plot::RasterProjection::planar(planar_basis),
         .omit_end = true,
         .balanced_sampling = false});
  }
  HS_EXPECT_SIZE_OR_RETURN(implicit_default.plotted,
                           selectable_default.plotted.size());
  HS_EXPECT_SIZE_OR_RETURN(implicit_default.alphas,
                           selectable_default.alphas.size());
  HS_EXPECT_GT(implicit_default.plotted.size(), size_t{0});
  const size_t compared = implicit_default.plotted.size();
  for (size_t i = 0; i < compared; ++i) {
    HS_EXPECT_NEAR(implicit_default.plotted[i].x,
                   selectable_default.plotted[i].x, POLICY_TOL);
    HS_EXPECT_NEAR(implicit_default.plotted[i].y,
                   selectable_default.plotted[i].y, POLICY_TOL);
    HS_EXPECT_NEAR(implicit_default.plotted[i].z,
                   selectable_default.plotted[i].z, POLICY_TOL);
    HS_EXPECT_NEAR(implicit_default.alphas[i], selectable_default.alphas[i],
                   POLICY_TOL);
  }
}

/** @brief Balanced sampling leaves one-dot edges unchanged. */
inline void test_rasterize_balanced_sampling_scope() {
  constexpr int W = 128, H = 64;
  const math::Basis basis = basis_from_normal(math::Vector(0, 1, 0));

  auto compare = [&](const math::Vector &start, const math::Vector &end) {
    ScratchScope sc(plot_arena());
    Fragments points;
    points.bind(plot_arena(), 2);
    Fragment a, b;
    a.pos = start;
    b.pos = end;
    points.push_back(a);
    points.push_back(b);
    auto capture = [&]<Plot::RasterSamplingPolicy Policy>() {
      hs_test::StubEffect fx(W, H);
      AlphaCapturePipeline pipeline;
      Canvas canvas(fx);
      auto shader = [](const math::Vector &, Fragment &f) {
        f.color = Color4(Pixel(65535, 65535, 65535), 0.4f);
      };
      Plot::rasterize<W, H,
                      Plot::RasterConfig{.single_pass = true,
                                         .sampling_policy = Policy}>(
          pipeline, canvas, points, shader,
          {.projection = Plot::RasterProjection::planar(basis),
           .balanced_sampling =
               Policy == Plot::RasterSamplingPolicy::SELECTABLE});
      return pipeline;
    };
    const AlphaCapturePipeline standard =
        capture.template operator()<Plot::RasterSamplingPolicy::DEFAULT>();
    const AlphaCapturePipeline balanced =
        capture.template operator()<Plot::RasterSamplingPolicy::SELECTABLE>();
    HS_EXPECT_SIZE_OR_RETURN(standard.plotted, balanced.plotted.size());
    HS_EXPECT_SIZE_OR_RETURN(standard.alphas, balanced.alphas.size());
    HS_EXPECT_GT(standard.plotted.size(), size_t{0});
    const size_t compared = standard.plotted.size();
    for (size_t i = 0; i < compared; ++i) {
      HS_EXPECT_EQ(std::bit_cast<uint32_t>(standard.plotted[i].x),
                   std::bit_cast<uint32_t>(balanced.plotted[i].x));
      HS_EXPECT_EQ(std::bit_cast<uint32_t>(standard.plotted[i].y),
                   std::bit_cast<uint32_t>(balanced.plotted[i].y));
      HS_EXPECT_EQ(std::bit_cast<uint32_t>(standard.plotted[i].z),
                   std::bit_cast<uint32_t>(balanced.plotted[i].z));
      HS_EXPECT_EQ(std::bit_cast<uint32_t>(standard.alphas[i]),
                   std::bit_cast<uint32_t>(balanced.alphas[i]));
    }
  };

  compare(disk_point(basis, 0.7f, 0.0f), disk_point(basis, 0.702f, 0.001f));
}

/** @brief Balanced long edges trade sample density for alpha-weighted coverage. */
inline void test_rasterize_balanced_sampling_density_and_alpha() {
  constexpr int W = 128, H = 64;
  constexpr float BASE_STEP = 2.0f * math::PI_F / W;
  const math::Basis basis = basis_from_normal(math::Vector(0, 1, 0));
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 2);
  Fragment a, b;
  a.pos = disk_point(basis, 0.35f, -0.2f);
  b.pos = disk_point(basis, 1.4f, 1.1f);
  points.push_back(a);
  points.push_back(b);

  auto capture = [&]<Plot::RasterSamplingPolicy Policy>() {
    hs_test::StubEffect fx(W, H);
    AlphaCapturePipeline pipeline;
    Canvas canvas(fx);
    auto shader = [](const math::Vector &, Fragment &f) {
      f.color = Color4(Pixel(65535, 65535, 65535), 0.4f);
    };
    Plot::rasterize<W, H,
                    Plot::RasterConfig{.single_pass = true,
                                       .sampling_policy = Policy}>(
        pipeline, canvas, points, shader,
        {.projection = Plot::RasterProjection::planar(basis),
         .balanced_sampling =
             Policy == Plot::RasterSamplingPolicy::SELECTABLE});
    return pipeline;
  };
  const AlphaCapturePipeline standard =
      capture.template operator()<Plot::RasterSamplingPolicy::DEFAULT>();
  const AlphaCapturePipeline balanced =
      capture.template operator()<Plot::RasterSamplingPolicy::SELECTABLE>();

  HS_EXPECT_GT(standard.plotted.size(), balanced.plotted.size());
  HS_EXPECT_GE(balanced.plotted.size() * 5, standard.plotted.size() * 3);
  HS_EXPECT_LE((max_projected_gap<W, H>(balanced.plotted)), 1.3f);
  const Plot::PlanarEdgeSampler sampler = planar_sampler(a.pos, b.pos, basis);
  const Plot::SamplePT first = sampler.one_pass(0.0f);
  const float default_step =
      Plot::screen_step<W, H>(first.pos, first.tan, BASE_STEP);
  const float candidate_step =
      std::min(BASE_STEP, default_step * (Plot::BALANCED_SCREEN_STEP_PX /
                                          Plot::SCREEN_STEP_PX));
  HS_EXPECT_TRUE(!standard.alphas.empty() && !balanced.alphas.empty());
  if (standard.alphas.empty() || balanced.alphas.empty())
    return;
  HS_EXPECT_NEAR(standard.alphas.front(), 0.4f, 1e-6f);
  HS_EXPECT_NEAR(
      balanced.alphas.front(),
      Plot::balanced_sample_alpha(0.4f, candidate_step / default_step), 1e-6f);
  // Newton correction re-normalizes the off-unit Bhaskara sin/cos positions.
  for (const math::Vector &point : balanced.plotted)
    HS_EXPECT_NEAR(point.length(), 1.0f, 2e-5f);
}

/** @brief Balanced sampling keeps default cadence in the polar clamp region. */
inline void test_rasterize_balanced_pole_guard() {
  constexpr int W = 144, H = 72;
  const math::Basis basis = basis_from_normal(math::Y_AXIS);
  auto near_pole = [&](float azimuth) {
    const math::Vector radial =
        basis.u * cosf(azimuth) + basis.w * sinf(azimuth);
    return (basis.v * cosf(0.02f) + radial * sinf(0.02f)).normalized();
  };
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 2);
  Fragment a, b;
  a.pos = near_pole(-2.4f);
  b.pos = near_pole(2.4f);
  points.push_back(a);
  points.push_back(b);

  constexpr float POLICY_TOL = 1e-4f;
  auto shader = [](const math::Vector &, Fragment &f) {
    f.color = Color4(Pixel(65535, 65535, 65535), 0.4f);
  };
  AlphaCapturePipeline pipeline;
  {
    hs_test::StubEffect fx(W, H);
    Canvas canvas(fx);
    Plot::rasterize<W, H,
                    Plot::RasterConfig{
                        .single_pass = true,
                        .derive_planar_arc_registers = false,
                        .interpolate_registers = false,
                        .sampling_policy =
                            Plot::RasterSamplingPolicy::SELECTABLE}>(
        pipeline, canvas, points, shader,
        {.projection = Plot::RasterProjection::planar(basis),
         .balanced_sampling = true});
  }
  AlphaCapturePipeline standard;
  {
    hs_test::StubEffect fx(W, H);
    Canvas canvas(fx);
    Plot::rasterize<W, H,
                    Plot::RasterConfig{
                        .single_pass = true,
                        .derive_planar_arc_registers = false,
                        .interpolate_registers = false,
                        .sampling_policy =
                            Plot::RasterSamplingPolicy::DEFAULT}>(
        standard, canvas, points, shader,
        {.projection = Plot::RasterProjection::planar(basis)});
  }
  HS_EXPECT_GT(standard.plotted.size(), size_t{2});
  HS_EXPECT_SIZE_OR_RETURN(pipeline.plotted, standard.plotted.size());
  for (size_t i = 0; i < standard.plotted.size(); ++i) {
    HS_EXPECT_NEAR(pipeline.plotted[i].x, standard.plotted[i].x, POLICY_TOL);
    HS_EXPECT_NEAR(pipeline.plotted[i].y, standard.plotted[i].y, POLICY_TOL);
    HS_EXPECT_NEAR(pipeline.plotted[i].z, standard.plotted[i].z, POLICY_TOL);
  }
  for (float alpha : pipeline.alphas)
    HS_EXPECT_NEAR(alpha, 0.4f, 1e-6f);
}

/**
 * @brief Balanced geodesic edges take the sparser steps and the alpha gain
 * without the planar step reuse.
 * @details BALANCED and SELECTABLE-on are distinct instantiations of the same
 * expressions, so the pair is held to POLICY_TOL.
 */
inline void test_rasterize_balanced_geodesic_density_and_alpha() {
  constexpr int W = 128, H = 64;
  constexpr float BASE_STEP = 2.0f * math::PI_F / W;
  constexpr float POLICY_TOL = 1e-4f;
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 2);
  Fragment a, b;
  a.pos = math::Vector(0.8f, 0.3f, 0.5196152f).normalized();
  b.pos = math::Vector(-0.2f, -0.7f, 0.6855655f).normalized();
  points.push_back(a);
  points.push_back(b);

  auto capture = [&]<Plot::RasterSamplingPolicy Policy>() {
    hs_test::StubEffect fx(W, H);
    AlphaCapturePipeline pipeline;
    Canvas canvas(fx);
    Plot::g_planar_full_samples = 0;
    Plot::g_planar_position_samples = 0;
    auto shader = [](const math::Vector &, Fragment &f) {
      f.color = Color4(Pixel(65535, 65535, 65535), 0.4f);
    };
    Plot::rasterize<W, H,
                    Plot::RasterConfig{.single_pass = true,
                                       .sampling_policy = Policy}>(
        pipeline, canvas, points, shader,
        {.balanced_sampling =
             Policy == Plot::RasterSamplingPolicy::SELECTABLE});
    HS_EXPECT_EQ(Plot::g_planar_full_samples, uint32_t{0});
    HS_EXPECT_EQ(Plot::g_planar_position_samples, uint32_t{0});
    return pipeline;
  };
  const AlphaCapturePipeline standard =
      capture.template operator()<Plot::RasterSamplingPolicy::DEFAULT>();
  const AlphaCapturePipeline balanced =
      capture.template operator()<Plot::RasterSamplingPolicy::SELECTABLE>();
  const AlphaCapturePipeline always_balanced =
      capture.template operator()<Plot::RasterSamplingPolicy::BALANCED>();

  HS_EXPECT_GT(standard.plotted.size(), balanced.plotted.size());
  HS_EXPECT_GE(balanced.plotted.size() * 5, standard.plotted.size() * 3);
  HS_EXPECT_LE((max_projected_gap<W, H>(balanced.plotted)), 1.3f);
  const Plot::GeodesicEdgeSpan es = Plot::make_geodesic_edge_span(a.pos, b.pos);
  const Plot::GeodesicEdgeSampler sampler{a.pos, math::cross(es.axis, a.pos),
                                          es.total};
  const Plot::SamplePT first = sampler(0.0f);
  const float default_step =
      Plot::screen_step<W, H>(first.pos, first.tan, BASE_STEP);
  HS_EXPECT_GT(default_step, BASE_STEP * Plot::MIN_POLE_SCALE *
                                 Plot::BALANCED_POLE_GUARD_SCALE);
  const float candidate_step =
      std::min(BASE_STEP, default_step * (Plot::BALANCED_SCREEN_STEP_PX /
                                          Plot::SCREEN_STEP_PX));
  HS_EXPECT_GT(candidate_step, default_step);
  HS_EXPECT_TRUE(!standard.alphas.empty() && !balanced.alphas.empty());
  if (standard.alphas.empty() || balanced.alphas.empty())
    return;
  HS_EXPECT_NEAR(standard.alphas.front(), 0.4f, 1e-6f);
  HS_EXPECT_NEAR(
      balanced.alphas.front(),
      Plot::balanced_sample_alpha(0.4f, candidate_step / default_step), 1e-6f);
  for (const math::Vector &point : balanced.plotted)
    HS_EXPECT_NEAR(point.length(), 1.0f, 2e-5f);

  HS_EXPECT_SIZE_OR_RETURN(always_balanced.plotted, balanced.plotted.size());
  HS_EXPECT_SIZE_OR_RETURN(always_balanced.alphas, balanced.alphas.size());
  const size_t compared = always_balanced.plotted.size();
  for (size_t i = 0; i < compared; ++i) {
    HS_EXPECT_NEAR(always_balanced.plotted[i].x, balanced.plotted[i].x,
                   POLICY_TOL);
    HS_EXPECT_NEAR(always_balanced.plotted[i].y, balanced.plotted[i].y,
                   POLICY_TOL);
    HS_EXPECT_NEAR(always_balanced.plotted[i].z, balanced.plotted[i].z,
                   POLICY_TOL);
    HS_EXPECT_NEAR(always_balanced.alphas[i], balanced.alphas[i], POLICY_TOL);
  }
}

/**
 * @brief A near-opaque balanced stroke saturates every sample's alpha at 1.
 * @details A high-latitude geodesic keeps the default step under 0.8 base
 * steps, so the balanced ratio is the full 1.25 on every sample.
 */
inline void test_rasterize_balanced_high_alpha_saturates() {
  constexpr int W = 128, H = 64;
  constexpr float BASE_STEP = 2.0f * math::PI_F / W;
  constexpr float ALPHA = 0.95f;
  auto sphere_point = [](float colatitude, float longitude) {
    const float radial = sinf(colatitude);
    return math::Vector(radial * cosf(longitude), cosf(colatitude),
                        radial * sinf(longitude))
        .normalized();
  };
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 2);
  Fragment a, b;
  a.pos = sphere_point(0.35f, -0.5f);
  b.pos = sphere_point(0.35f, 0.5f);
  points.push_back(a);
  points.push_back(b);

  auto capture = [&]<Plot::RasterSamplingPolicy Policy>() {
    hs_test::StubEffect fx(W, H);
    AlphaCapturePipeline pipeline;
    Canvas canvas(fx);
    auto shader = [](const math::Vector &, Fragment &f) {
      f.color = Color4(Pixel(65535, 65535, 65535), ALPHA);
    };
    Plot::rasterize<W, H,
                    Plot::RasterConfig{.single_pass = true,
                                       .sampling_policy = Policy}>(
        pipeline, canvas, points, shader,
        {.omit_end = true,
         .balanced_sampling =
             Policy == Plot::RasterSamplingPolicy::SELECTABLE});
    return pipeline;
  };
  const AlphaCapturePipeline standard =
      capture.template operator()<Plot::RasterSamplingPolicy::DEFAULT>();
  const AlphaCapturePipeline balanced =
      capture.template operator()<Plot::RasterSamplingPolicy::SELECTABLE>();

  const Plot::GeodesicEdgeSpan es = Plot::make_geodesic_edge_span(a.pos, b.pos);
  const Plot::GeodesicEdgeSampler sampler{a.pos, math::cross(es.axis, a.pos),
                                          es.total};
  const Plot::SamplePT first = sampler(0.0f);
  const float default_step =
      Plot::screen_step<W, H>(first.pos, first.tan, BASE_STEP);
  HS_EXPECT_GT(default_step, BASE_STEP * Plot::MIN_POLE_SCALE *
                                 Plot::BALANCED_POLE_GUARD_SCALE);
  HS_EXPECT_LE(default_step *
                   (Plot::BALANCED_SCREEN_STEP_PX / Plot::SCREEN_STEP_PX),
               BASE_STEP);
  HS_EXPECT_GT(standard.plotted.size(), balanced.plotted.size());
  HS_EXPECT_GT(balanced.alphas.size(), size_t{4});
  for (float alpha : standard.alphas)
    HS_EXPECT_EQ(alpha, ALPHA);
  for (float alpha : balanced.alphas)
    HS_EXPECT_EQ(alpha, 1.0f);
}

/**
 * @brief Balanced planar stars retain coverage and energy within the clipped budget.
 * @details The clip bounds the visited column span. Render-band pixels match
 * exactly under IEEE math; -ffast-math reassociates accumulated coverage, so
 * the tile allows a 16-level channel difference there.
 */
inline void test_rasterize_balanced_star_visual_budget() {
  constexpr int W = 144, H = 72;
#if defined(HS_TEST_FAST_MATH)
  constexpr int CLIP_CHANNEL_TOL = 16;
#else
  constexpr int CLIP_CHANNEL_TOL = 0;
#endif
  struct StarState {
    math::Quaternion orientation;
    float radius;
    int sides;
    float phase;
  };
  const std::array<StarState, 4> states = {{
      {math::Quaternion(), 0.45f, 4, 0.0f},
      {math::Quaternion(0.93f, -0.11f, 0.24f, 0.25f).normalized(), 0.98f, 7,
       0.37f},
      {math::Quaternion(0.81f, 0.32f, -0.29f, 0.39f).normalized(), 1.02f, 7,
       1.2f},
      {math::Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized(), 1.72f, 16,
       2.1f},
  }};
  struct Frame {
    std::vector<Pixel> pixels;
    uint32_t backstops = 0;
    uint32_t full_samples = 0;
    uint32_t position_samples = 0;
  };

  auto render = [&]<Plot::RasterSamplingPolicy Policy>(const StarState &state,
                                                       const ClipRegion *clip) {
    hs_test::StubEffect fx(W, H);
    if (clip != nullptr)
      fx.set_clip(clip->y_start, clip->y_end, clip->x_start, clip->x_end);
    hs::g_scan_metrics.reset();
    Plot::g_planar_full_samples = 0;
    Plot::g_planar_position_samples = 0;
    {
      ScratchScope sc(plot_arena());
      Fragments points;
      points.bind(plot_arena(), static_cast<size_t>(state.sides * 2 + 2));
      const math::Basis basis =
          math::make_basis(state.orientation, math::X_AXIS);
      Plot::Star<Plot::PlanarProjection>::sample_positions(
          points, basis, state.radius, state.sides, state.phase);
      const math::Basis planar_basis = Plot::planar_chart_basis(
          math::get_antipode(basis, state.radius).first.v);
      Filter::Screen::DirectAntiAliasSink<W, H> sink;
      Canvas canvas(fx);
      initialize_parity_frame<W, H>(canvas);
      sink.prepare(canvas);
      auto shader = [](const math::Vector &, Fragment &f) {
        f.color = Color4(Pixel(65535, 65535, 65535), 0.32f);
      };
      Plot::rasterize<W, H,
                      Plot::RasterConfig{.single_pass = true,
                                         .derive_planar_arc_registers = false,
                                         .interpolate_registers = false,
                                         .sampling_policy = Policy}>(
          sink, canvas, points, shader,
          {.projection = Plot::RasterProjection::planar(planar_basis),
           .omit_end = true,
           .balanced_sampling =
               Policy == Plot::RasterSamplingPolicy::SELECTABLE});
    }
    fx.advance_display();
    Frame frame;
    frame.backstops = hs::g_scan_metrics.plot_backstop_hits;
    frame.full_samples = Plot::g_planar_full_samples;
    frame.position_samples = Plot::g_planar_position_samples;
    frame.pixels.resize(static_cast<size_t>(W) * H);
    hs_test::capture_frame<W, H>(fx, frame.pixels);
    return frame;
  };

  auto energy = [](const Frame &frame) {
    uint64_t total = 0;
    for (const Pixel &pixel : frame.pixels)
      total += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
    return total;
  };
  auto covered = [](const Pixel &pixel) {
    return static_cast<uint32_t>(pixel.r) + pixel.g + pixel.b > 512;
  };

  uint32_t balanced_full_samples = 0;
  uint32_t balanced_position_samples = 0;
  size_t margin_reference_lit = 0;
  for (const StarState &state : states) {
    const Frame standard =
        render.template operator()<Plot::RasterSamplingPolicy::DEFAULT>(
            state, nullptr);
    const Frame balanced =
        render.template operator()<Plot::RasterSamplingPolicy::SELECTABLE>(
            state, nullptr);
    HS_EXPECT_EQ(standard.position_samples, uint32_t{0});
    HS_EXPECT_GT(balanced.position_samples, uint32_t{0});
    balanced_full_samples += balanced.full_samples;
    balanced_position_samples += balanced.position_samples;
    HS_EXPECT_EQ(standard.backstops, uint32_t{0});
    HS_EXPECT_EQ(balanced.backstops, uint32_t{0});
    const double energy_ratio =
        static_cast<double>(energy(balanced)) / energy(standard);
    HS_EXPECT_GT(energy_ratio, 0.98);
    HS_EXPECT_LT(energy_ratio, 1.06);

    size_t standard_coverage = 0;
    size_t balanced_coverage = 0;
    size_t uncovered = 0;
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        const size_t index = static_cast<size_t>(y) * W + x;
        standard_coverage += covered(standard.pixels[index]);
        balanced_coverage += covered(balanced.pixels[index]);
        if (!covered(standard.pixels[index]))
          continue;
        bool neighborhood = false;
        for (int dy = -1; dy <= 1 && !neighborhood; ++dy) {
          const int sy = y + dy;
          if (sy < 0 || sy >= H)
            continue;
          for (int dx = -1; dx <= 1; ++dx) {
            const int sx = (x + dx + W) % W;
            if (covered(balanced.pixels[static_cast<size_t>(sy) * W + sx])) {
              neighborhood = true;
              break;
            }
          }
        }
        uncovered += !neighborhood;
      }
    }
    HS_EXPECT_GT(standard_coverage, size_t{20});
    HS_EXPECT_GT(balanced_coverage * 20, standard_coverage * 19);
    HS_EXPECT_LT(balanced_coverage * 25, standard_coverage * 27);
    HS_EXPECT_EQ(uncovered, size_t{0});

    const std::array<ClipRegion, 4> clips = {{
        {0, H / 2, 0, W / 2},
        {0, H / 2, W / 2, W},
        {H / 2, H, 0, W / 2},
        {H / 2, H, W / 2, W},
    }};
    for (const ClipRegion &clip : clips) {
      const Frame tile =
          render.template operator()<Plot::RasterSamplingPolicy::SELECTABLE>(
              state, &clip);
      HS_EXPECT_EQ(tile.backstops, uint32_t{0});
      ClipRegion render_clip = clip;
      render_clip.w = W;
      render_clip.h = H;
      for (int y = render_clip.render_y_start(); y < render_clip.render_y_end();
           ++y)
        for (int x = 0; x < W; ++x) {
          if (!render_clip.contains_x(x))
            continue;
          const size_t index = static_cast<size_t>(y) * W + x;
          const bool in_display = y >= clip.y_start && y < clip.y_end &&
                                  x >= clip.x_start && x < clip.x_end;
          margin_reference_lit +=
              !in_display && covered(balanced.pixels[index]);
          HS_CONTEXT("balanced star render band", x, y);
          HS_EXPECT_NEAR(tile.pixels[index].r, balanced.pixels[index].r,
                         CLIP_CHANNEL_TOL);
          HS_EXPECT_NEAR(tile.pixels[index].g, balanced.pixels[index].g,
                         CLIP_CHANNEL_TOL);
          HS_EXPECT_NEAR(tile.pixels[index].b, balanced.pixels[index].b,
                         CLIP_CHANNEL_TOL);
        }
    }
  }
  HS_EXPECT_GT(margin_reference_lit, size_t{20});
  HS_EXPECT_GT(balanced_position_samples * 4,
               balanced_full_samples + balanced_position_samples);
}

/** @brief Single-pass geodesics preserve the open-line endpoint contract. */
inline void test_rasterize_single_pass_geodesic_endpoints_and_omit_end() {
  constexpr int W = 128, H = 64;
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 2);

  Fragment a, b;
  a.pos = math::Vector(0.8f, 0.3f, 0.5196152f).normalized();
  b.pos = math::Vector(-0.2f, -0.7f, 0.6855655f).normalized();
  points.push_back(a);
  points.push_back(b);

  auto draw = [&](bool omit_end) {
    hs_test::StubEffect fx(W, H);
    CapturePipeline pipeline;
    Canvas canvas(fx);
    Plot::rasterize<W, H, Plot::RasterConfig{.single_pass = true}>(
        pipeline, canvas, points, noop_shader, {.omit_end = omit_end});
    return pipeline.plotted;
  };
  const std::vector<math::Vector> complete = draw(false);
  const std::vector<math::Vector> omitted = draw(true);

  HS_EXPECT_GT(omitted.size(), size_t{2});
  HS_EXPECT_SIZE_OR_RETURN(complete, omitted.size() + 1);
  HS_EXPECT_NEAR(math::angle_between(complete.front(), a.pos), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(math::angle_between(complete.back(), b.pos), 0.0f, 1e-3f);
  for (size_t i = 0; i < omitted.size(); ++i)
    HS_EXPECT_NEAR(math::angle_between(complete[i], omitted[i]), 0.0f, 1e-3f);
}

/** @brief Single-pass geodesics remain gap-free at poles, seams, and long arcs. */
inline void test_rasterize_single_pass_geodesic_stress_arcs_are_gap_free() {
  constexpr int W = 128, H = 64;
  auto sphere_point = [](float colatitude, float longitude) {
    const float radial = sinf(colatitude);
    return math::Vector(radial * cosf(longitude), cosf(colatitude),
                        radial * sinf(longitude))
        .normalized();
  };
  const std::array<std::pair<math::Vector, math::Vector>, 6> arcs = {{
      {math::Vector(0, 1, 0), sphere_point(0.9f, 1.1f)},
      {sphere_point(math::PI_F - 0.8f, -0.7f), math::Vector(0, -1, 0)},
      {sphere_point(0.45f, 0.0f), sphere_point(0.45f, math::PI_F)},
      {sphere_point(1.2f, math::PI_F - 0.12f),
       sphere_point(1.3f, -math::PI_F + 0.14f)},
      {sphere_point(0.35f, -1.0f), sphere_point(math::PI_F - 0.4f, 2.05f)},
      {sphere_point(1.0f, 0.25f), sphere_point(2.1f, math::PI_F - 0.25f)},
  }};

  for (const auto &[start, end] : arcs) {
    hs_test::StubEffect fx(W, H);
    ScratchScope sc(plot_arena());
    Fragments points;
    points.bind(plot_arena(), 2);
    Fragment a, b;
    a.pos = start;
    b.pos = end;
    points.push_back(a);
    points.push_back(b);

    CapturePipeline pipeline;
    Canvas canvas(fx);
    Plot::rasterize<W, H, Plot::RasterConfig{.single_pass = true}>(
        pipeline, canvas, points, noop_shader);

    HS_EXPECT_GT(pipeline.plotted.size(), size_t{2});
    HS_EXPECT_LE((max_projected_gap<W, H>(pipeline.plotted)), 1.5f);
    if (pipeline.plotted.empty())
      continue;
    HS_EXPECT_NEAR(math::angle_between(pipeline.plotted.front(), start), 0.0f,
                   1e-3f);
    HS_EXPECT_NEAR(math::angle_between(pipeline.plotted.back(), end), 0.0f,
                   1e-3f);
    for (const math::Vector &p : pipeline.plotted)
      HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
  }
}

/**
 * @brief Single-pass direct-AA quadrant render bands match the full render.
 */
inline void test_rasterize_single_pass_geodesic_quadrant_clip_parity() {
  constexpr int W = 96, H = 48;
  const std::array<std::pair<int, int>, 10> control_pixels = {{
      {4, 4},
      {W / 2 - 5, H / 2 - 5},
      {W / 2 + 5, 5},
      {W - 5, H / 2 - 5},
      {W - 5, H / 2 + 5},
      {W / 2 + 5, H - 5},
      {W / 2 - 5, H / 2 + 5},
      {5, H - 5},
      {2, H / 2},
      {W - 2, H / 2},
  }};
  auto shade = [](const math::Vector &, Fragment &f) {
    f.color = Color4(Pixel(65535, 65535, 65535), 0.8f);
  };
  auto render = [&](hs_test::StubEffect &fx) {
    Filter::Screen::DirectAntiAliasSink<W, H> sink;
    ScratchScope sc(plot_arena());
    Fragments points;
    points.bind(plot_arena(), control_pixels.size());
    for (const auto &[x, y] : control_pixels) {
      Fragment f;
      f.pos = math::pixel_to_vector<W, H>(x, y);
      points.push_back(f);
    }

    Canvas canvas(fx);
    initialize_parity_frame<W, H>(canvas);
    sink.prepare(canvas);
    const ClipRegion &clip = canvas.clip();
    const ClipRegion::XClip x_clip = clip.x_clip();
    std::array<uint8_t, control_pixels.size() - 1> bits{};
    std::span<const uint8_t> visible;
    if (!clip.is_full()) {
      if (!Plot::gate_trail_edges<W, H>(sink, clip, x_clip, points,
                                        bits.data()))
        return;
      visible = bits;
    }
    Plot::rasterize<W, H, Plot::RasterConfig{.single_pass = true}>(
        sink, canvas, points, shade,
        {.projection = Plot::RasterProjection::geodesic(visible)});
  };

  std::vector<Pixel> reference(static_cast<size_t>(W) * H);
  {
    hs_test::StubEffect fx(W, H);
    render(fx);
    fx.advance_display();
    hs_test::capture_frame<W, H>(fx, reference);
  }

  const int quadrants[4][4] = {
      {0, H / 2, 0, W / 2},
      {0, H / 2, W / 2, W},
      {H / 2, H, 0, W / 2},
      {H / 2, H, W / 2, W},
  };
  int margin_lit = 0;
  for (const auto &quadrant : quadrants) {
    hs_test::StubEffect fx(W, H);
    fx.set_clip(quadrant[0], quadrant[1], quadrant[2], quadrant[3]);
    render(fx);
    fx.advance_display();
    const RenderBandDiff diff = render_band_diff<W>(fx, reference);
    margin_lit += diff.margin_lit;
    HS_EXPECT_GT(diff.lit, 10);
    expect_render_band_parity("single-pass direct AA", diff);
  }
  HS_EXPECT_GT(margin_lit, 10);
}
