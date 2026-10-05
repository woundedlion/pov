/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_plot_scan.h.

// ============================================================================
// Plot::rasterize — control-flow coverage through a capturing pipeline
// ============================================================================

/**
 * @brief A sub-base_step open segment takes the fast path and plots BOTH
 *        endpoints (start, then the final vertex since the loop is open).
 */
inline void test_rasterize_subpixel_open_segment_plots_both_endpoints() {
  constexpr int W = 128, H = 64;
  hs_test::StubEffect fx(W, H);
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 4);

  Fragment a, b;
  a.pos = math::Vector(1, 0, 0);
  // ~0.02 rad apart, well under base_step (2*pi/W = ~0.049 rad).
  b.pos = math::Vector(1.0f, 0.02f, 0.0f).normalized();
  points.push_back(a);
  points.push_back(b);

  CapturePipeline pipe;
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(pipe, c, points, noop_shader);
  }
  fx.advance_display();

  // Fast path on an open last segment plots curr and next.
  HS_EXPECT_SIZE_OR_RETURN(pipe.plotted, (size_t)2);
  // Chord, not angle_between: acos' derivative diverges at |dot| = 1, so the
  // angle a unit pair reports quantizes in ~3.5e-4 steps and cannot resolve a
  // tolerance this tight. The chord tracks the angle to within angle^3/24.
  HS_EXPECT_NEAR((pipe.plotted.front() - a.pos).length(), 0.0f, 1e-4f);
  HS_EXPECT_NEAR((pipe.plotted.back() - b.pos).length(), 0.0f, 1e-4f);
  for (const math::Vector &p : pipe.plotted)
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
}

/**
 * @brief A normal-length open segment is sampled densely enough to be gap-free
 *        (no inter-sample gap exceeds one pixel column) and lands on both
 *        endpoints.
 */
inline void test_rasterize_open_segment_gap_free() {
  constexpr int W = 128, H = 64;
  constexpr float base_step = (2.0f * math::PI_F) / W;
  hs_test::StubEffect fx(W, H);
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 4);

  Fragment a, b;
  a.pos = math::Vector(1, 0, 0);
  b.pos = math::Vector(0, 1, 0);
  points.push_back(a);
  points.push_back(b);

  CapturePipeline pipe;
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(pipe, c, points, noop_shader);
  }
  fx.advance_display();

  HS_EXPECT_GT(pipe.plotted.size(), (size_t)10);
  // Adaptive sub-stepping caps each advance at ~base_step (slack for quantization).
  HS_EXPECT_LE(max_consecutive_gap(pipe.plotted, /*wrap=*/false),
               1.5f * base_step);
  HS_EXPECT_TRUE(!pipe.plotted.empty());
  if (pipe.plotted.empty())
    return;
  HS_EXPECT_NEAR(math::angle_between(pipe.plotted.front(), a.pos), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(math::angle_between(pipe.plotted.back(), b.pos), 0.0f, 1e-3f);
  for (const math::Vector &p : pipe.plotted)
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
}

/**
 * @brief A closed loop omits each segment's terminal vertex, so the plotted
 *        sequence wraps continuously (gap-free across the seam) with no shared
 *        vertex plotted twice.
 */
inline void test_rasterize_closed_loop_gap_free_no_dup() {
  constexpr int W = 128, H = 64;
  constexpr float base_step = (2.0f * math::PI_F) / W;
  hs_test::StubEffect fx(W, H);
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 4);

  // A spherical triangle with well-separated vertices.
  Fragment v0, v1, v2;
  v0.pos = math::Vector(1, 0, 0);
  v1.pos = math::Vector(0, 1, 0);
  v2.pos = math::Vector(0, 0, 1);
  points.push_back(v0);
  points.push_back(v1);
  points.push_back(v2);

  CapturePipeline pipe;
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(pipe, c, points, noop_shader,
                          {.loop = Plot::RasterLoop::closed()});
  }
  fx.advance_display();

  HS_EXPECT_GT(pipe.plotted.size(), (size_t)20);
  // Continuity all the way around, including the last->first seam.
  HS_EXPECT_LE(max_consecutive_gap(pipe.plotted, /*wrap=*/true),
               1.5f * base_step);
  // No vertex plotted twice: consecutive samples stay distinct.
  for (size_t i = 0; i < pipe.plotted.size(); ++i) {
    const math::Vector delta =
        pipe.plotted[i] - pipe.plotted[(i + 1) % pipe.plotted.size()];
    HS_EXPECT_GT(math::dot(delta, delta), 1e-10f);
  }
}

/**
 * @brief A planar segment whose endpoint sits at the basis antipode falls back
 *        to the geodesic strategy: rasterizing with that planar_basis must
 *        produce exactly the same samples as rasterizing geodesically.
 */
inline void test_rasterize_antipodal_seam_planar_falls_back_geodesic() {
  constexpr int W = 128, H = 64;
  hs_test::StubEffect fx(W, H);
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 4);

  // planar basis centered on +Z; the second endpoint sits within the antipodal
  // margin of -Z (dot < -COS_PLANAR_ANTIPODE), tripping the seam guard.
  math::Basis basis = basis_from_normal(math::Vector(0, 0, 1));
  Fragment a, b;
  a.pos = math::Vector(1, 0, 0);
  b.pos = math::Vector(0.02f, 0.0f, -1.0f).normalized();
  HS_EXPECT_LT(math::dot(b.pos, basis.v), -Plot::COS_PLANAR_ANTIPODE);
  points.push_back(a);
  points.push_back(b);

  CapturePipeline planar_pipe, geo_pipe;
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(
        planar_pipe, c, points, noop_shader,
        {.projection = Plot::RasterProjection::planar(basis)});
  }
  fx.advance_display();
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(geo_pipe, c, points, noop_shader);
  }
  fx.advance_display();

  HS_EXPECT_SIZE_OR_RETURN(planar_pipe.plotted, geo_pipe.plotted.size());
  size_t n = std::min(planar_pipe.plotted.size(), geo_pipe.plotted.size());
  for (size_t i = 0; i < n; ++i)
    HS_EXPECT_NEAR((planar_pipe.plotted[i] - geo_pipe.plotted[i]).length(),
                   0.0f, 1e-5f);
}

/**
 * @brief A non-seam planar segment renders gap-free in ARC length: the
 *        arc-length parameterization (PlanarEdgeSampler's cumulative-arc inversion)
 *        keeps every plotted step near one pixel column and lands on both
 *        endpoints, with no clustering or gaps the projection-linear chord would
 *        otherwise leave.
 * @details Exercises the planar strategy path (rasterize_planar_strategy +
 *          PlanarEdgeSampler); the antipodal-seam case falls back to geodesic. Pins the
 *          end-to-end arc-uniform sampling the PLANAR_LEN_SAMPLES table provides; it
 *          does not isolate the table's contribution from the rasterizer's
 *          adaptive (sin-phi) sub-stepping, which also shapes local density.
 */
inline void test_rasterize_planar_segment_gap_free_arclength() {
  constexpr int W = 128, H = 64;
  constexpr float base_step = (2.0f * math::PI_F) / W;
  hs_test::StubEffect fx(W, H);
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 4);

  // Planar disk about +Y; endpoints sweep colatitude 0.3 -> 1.3 across azimuths
  // so the chord crosses regions of differing azimuthal stretch (r / sin r).
  math::Basis basis = basis_from_normal(math::Vector(0, 1, 0));
  Fragment a, b;
  a.pos = disk_point(basis, 0.3f, 0.0f);
  b.pos = disk_point(basis, 1.3f, 1.0f);
  // Stays out of the antipodal-seam fallback so the planar strategy is used.
  HS_EXPECT_GT(math::dot(a.pos, basis.v), -Plot::COS_PLANAR_ANTIPODE);
  HS_EXPECT_GT(math::dot(b.pos, basis.v), -Plot::COS_PLANAR_ANTIPODE);
  points.push_back(a);
  points.push_back(b);

  CapturePipeline pipe;
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(
        pipe, c, points, noop_shader,
        {.projection = Plot::RasterProjection::planar(basis)});
  }
  fx.advance_display();

  HS_EXPECT_GT(pipe.plotted.size(), (size_t)10);
  if (pipe.plotted.empty())
    return;
  HS_EXPECT_LE(max_consecutive_gap(pipe.plotted, /*wrap=*/false),
               1.5f * base_step);
  // Endpoints land within PlanarEdgeSampler's project/unproject round-trip error.
  HS_EXPECT_NEAR(math::angle_between(pipe.plotted.front(), a.pos), 0.0f, 1e-2f);
  HS_EXPECT_NEAR(math::angle_between(pipe.plotted.back(), b.pos), 0.0f, 1e-2f);
  for (const math::Vector &p : pipe.plotted)
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
}

/**
 * @brief A planar segment's v0/v1 registers measure the RENDERED azimuthal arc,
 *        not the geodesic chord: v1 rises monotonically to the planar arc length
 *        (strictly longer than the endpoints' great-circle distance) and v0 spans
 *        0..1 over that same arc. Locks in the rasterizer's arc-register
 *        re-derivation for planar-policy polygons and stars plus Flower, which
 *        the sample()-level tests cannot observe.
 */
inline void test_rasterize_planar_arc_registers_track_drawn_arc() {
  constexpr int W = 128, H = 64;
  hs_test::StubEffect fx(W, H);
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 4);

  // Same non-seam planar disk edge as the gap-free test: colatitude 0.3 -> 1.3
  // across azimuths, so the rendered edge bows well clear of its chord.
  math::Basis basis = basis_from_normal(math::Vector(0, 1, 0));
  Fragment a, b;
  a.pos = disk_point(basis, 0.3f, 0.0f);
  b.pos = disk_point(basis, 1.3f, 1.0f);
  HS_EXPECT_GT(math::dot(a.pos, basis.v), -Plot::COS_PLANAR_ANTIPODE);
  HS_EXPECT_GT(math::dot(b.pos, basis.v), -Plot::COS_PLANAR_ANTIPODE);
  // Bare control points default v0/v1 to 0, so any nonzero arc below comes
  // solely from the rasterizer's rendered-arc override.
  points.push_back(a);
  points.push_back(b);

  CapturePipeline pipe;
  std::vector<float> v0s, v1s;
  auto capture = [&](const math::Vector &, Fragment &f) {
    v0s.push_back(f.v0);
    v1s.push_back(f.v1);
  };
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(
        pipe, c, points, capture,
        {.projection = Plot::RasterProjection::planar(basis)});
  }
  fx.advance_display();

  HS_EXPECT_GT(v1s.size(), (size_t)10);

  for (size_t i = 1; i < v1s.size(); ++i) {
    HS_EXPECT_GE(v1s[i], v1s[i - 1] - 1e-6f);
    HS_EXPECT_GE(v0s[i], v0s[i - 1] - 1e-6f);
  }

  const float planar = Plot::planar_arc_length(a.pos, b.pos, basis);
  const Plot::PlanarEdgeSampler sampler =
      Plot::make_planar_edge_sampler(a.pos, b.pos, basis);
  float rendered = 0.0f;
  double independent_rendered = 0.0;
  const auto independent_angle = [](const math::Vector &u,
                                    const math::Vector &v) {
    const double uv = static_cast<double>(u.x) * v.x +
                      static_cast<double>(u.y) * v.y +
                      static_cast<double>(u.z) * v.z;
    const double uu = static_cast<double>(u.x) * u.x +
                      static_cast<double>(u.y) * u.y +
                      static_cast<double>(u.z) * u.z;
    const double vv = static_cast<double>(v.x) * v.x +
                      static_cast<double>(v.y) * v.y +
                      static_cast<double>(v.z) * v.z;
    return std::acos(std::clamp(uv / std::sqrt(uu * vv), -1.0, 1.0));
  };
  math::Vector previous = sampler.unproject(0.0f);
  const math::Vector mapped_start = previous;
  for (int i = 1; i <= 32; ++i) {
    const math::Vector current =
        sampler.unproject(static_cast<float>(i) / 32.0f);
    rendered += math::angle_between(previous, current);
    independent_rendered += independent_angle(previous, current);
    previous = current;
  }
  const double bow =
      independent_rendered - independent_angle(mapped_start, previous);
  HS_EXPECT_GT(bow, 1e-3);
  HS_EXPECT_LT(bow, 5e-3);
  HS_EXPECT_NEAR(planar, rendered, 2e-3f);
  HS_EXPECT_TRUE(!v1s.empty() && !v0s.empty());
  if (v1s.empty() || v0s.empty())
    return;
  HS_EXPECT_NEAR(v1s.front(), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(v1s.back(), planar, 1e-3f);
  HS_EXPECT_GT(v1s.back(), planar - bow * 0.5);

  // v0 is v1 normalized by the single-segment total arc: 0 at the start, ~1 end.
  HS_EXPECT_NEAR(v0s.front(), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(v0s.back(), 1.0f, 2e-2f);
}

/**
 * @brief Cull-sample reuse rebuilds the same planar sampler as a fresh build.
 * @details Both overloads evaluate identical expressions on identical inputs
 * (the reuse path's stride-2 read of the cull span's eighth points IS the
 * rebuild path's quarter points), so they agree bit for bit under IEEE. The
 * shipping -ffast-math pair reassociates each overload in its own inlining
 * context, which moves the shared unproject/angle chains, so PARITY_TOL is the
 * asserted contract; a wrong sample index or a dropped endpoint moves these
 * registers by ~0.1 rad.
 */
inline void test_planar_sampler_from_cull_parity() {
  constexpr int W = 128;
  constexpr float PARITY_TOL = 1e-4f;
  hs::random().seed(0xC0115A9u);
  for (int trial = 0; trial < 400; ++trial) {
    const float normal_x = hs::rand_f(-1.0f, 1.0f);
    const float normal_y = hs::rand_f(-1.0f, 1.0f);
    const float normal_z = hs::rand_f(-1.0f, 1.0f);
    math::Vector normal(normal_x, normal_y, normal_z);
    if (normal.length() < 0.01f)
      normal = math::Y_AXIS;
    const math::Basis basis = basis_from_normal(normal.normalized());
    auto point = [&](float radius, float angle) {
      const math::Vector radial = basis.u * cosf(angle) + basis.w * sinf(angle);
      return (basis.v * cosf(radius) + radial * sinf(radius)).normalized();
    };
    const float a_radius = hs::rand_f(0.01f, 2.8f);
    const float a_angle = hs::rand_f(-math::PI_F, math::PI_F);
    const math::Vector a = point(a_radius, a_angle);
    const float b_radius = hs::rand_f(0.01f, 2.8f);
    const float b_angle = hs::rand_f(-math::PI_F, math::PI_F);
    const math::Vector b = point(b_radius, b_angle);
    const Plot::PlanarEdgeSpan span = Plot::make_planar_edge_span(a, b, basis);
    int col_s, col_len;
    math::Vector end;
    Plot::planar_col_span<W>(a, basis, span, col_s, col_len, &end);
    const Plot::PlanarEdgeSampler rebuilt =
        Plot::make_planar_edge_sampler(a, b, basis);
    const Plot::PlanarEdgeSampler reused =
        Plot::make_planar_edge_sampler(span, end, basis);

    auto expect_parity = [](float lhs, float rhs) {
      HS_EXPECT_NEAR(lhs, rhs, PARITY_TOL);
    };
    expect_parity(rebuilt.proj1.first, reused.proj1.first);
    expect_parity(rebuilt.proj1.second, reused.proj1.second);
    expect_parity(rebuilt.dx, reused.dx);
    expect_parity(rebuilt.dy, reused.dy);
    expect_parity(rebuilt.chart_tangent.x, reused.chart_tangent.x);
    expect_parity(rebuilt.chart_tangent.y, reused.chart_tangent.y);
    expect_parity(rebuilt.chart_tangent.z, reused.chart_tangent.z);
    for (size_t i = 0; i < rebuilt.arc_cumul.size(); ++i)
      expect_parity(rebuilt.arc_cumul[i], reused.arc_cumul[i]);
    expect_parity(rebuilt.dist, reused.dist);
    int interval = 0;
    int position_interval = 0;
    for (float t : {0.0f, 0.17f, 0.51f, 0.88f}) {
      const Plot::SamplePT original = rebuilt.one_pass(t);
      const Plot::SamplePT optimized = reused.one_pass_monotonic(t, interval);
      const math::Vector position =
          reused.position_monotonic(t, position_interval);
      expect_parity(original.pos.x, optimized.pos.x);
      expect_parity(original.pos.y, optimized.pos.y);
      expect_parity(original.pos.z, optimized.pos.z);
      expect_parity(optimized.pos.x, position.x);
      expect_parity(optimized.pos.y, position.y);
      expect_parity(optimized.pos.z, position.z);
      expect_parity(original.tan.x, optimized.tan.x);
      expect_parity(original.tan.y, optimized.tan.y);
      expect_parity(original.tan.z, optimized.tan.z);
    }
  }
}

/**
 * @brief Planar compile-time policies preserve positions and default registers.
 * @details The register policies select what is derived, never how a position
 * is computed, so every stream walks the same expressions and the positions
 * agree bit for bit under IEEE. Each instantiation is reassociated in its own
 * inlining context under -ffast-math, so POLICY_TOL is the asserted contract.
 */
inline void test_rasterize_planar_policy_parity() {
  constexpr int W = 128, H = 64;
  constexpr float POLICY_TOL = 1e-4f;
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 16);
  const math::Basis shape_basis = math::make_basis(
      math::Quaternion(0.93f, -0.11f, 0.24f, 0.25f).normalized(), math::X_AXIS);
  constexpr float RADIUS = 0.74f;
  Plot::Star<Plot::PlanarProjection>::sample(points, shape_basis, RADIUS, 7,
                                             1.37f);
  const math::Basis planar_basis =
      Plot::planar_chart_basis(math::get_antipode(shape_basis, RADIUS).first.v);

  struct Stream {
    std::vector<math::Vector> positions;
    std::vector<std::array<uint32_t, 3>> registers;
  };
  auto capture = [&]<bool DerivePlanarArcRegisters, bool InterpolateRegisters>(
                     bool rebuild_sampler = false) {
    hs_test::StubEffect fx(W, H);
    fx.set_clip(0, H, 0, W / 2);
    DirectCapturePipeline pipeline;
    Stream stream;
    auto shader = [&](const math::Vector &, Fragment &f) {
      stream.registers.push_back({std::bit_cast<uint32_t>(f.v0),
                                  std::bit_cast<uint32_t>(f.v1),
                                  std::bit_cast<uint32_t>(f.v2)});
    };
    {
      Canvas canvas(fx);
      Plot::rasterize<W, H,
                      Plot::RasterConfig{.single_pass = true,
                                         .derive_planar_arc_registers =
                                             DerivePlanarArcRegisters,
                                         .interpolate_registers =
                                             InterpolateRegisters}>(
          pipeline, canvas, points, shader,
          {.projection = Plot::RasterProjection::planar(planar_basis),
           .omit_end = true,
           .rebuild_planar_sampler = rebuild_sampler});
    }
    fx.advance_display();
    stream.positions = std::move(pipeline.plotted);
    return stream;
  };

  const Stream derived = capture.template operator()<true, true>();
  const Stream source_registers = capture.template operator()<false, true>();
  const Stream positions_only = capture.template operator()<false, false>();
  const Stream rebuilt_positions_only =
      capture.template operator()<false, false>(true);
  HS_EXPECT_SIZE_OR_RETURN(derived.positions,
                           source_registers.positions.size());
  HS_EXPECT_SIZE_OR_RETURN(derived.positions, positions_only.positions.size());
  HS_EXPECT_SIZE_OR_RETURN(positions_only.positions,
                           rebuilt_positions_only.positions.size());
  HS_EXPECT_SIZE_OR_RETURN(derived.registers,
                           source_registers.registers.size());
  HS_EXPECT_SIZE_OR_RETURN(derived.registers, positions_only.registers.size());
  HS_EXPECT_SIZE_OR_RETURN(positions_only.registers,
                           rebuilt_positions_only.registers.size());
  size_t derived_differences = 0;
  size_t source_differences = 0;
  for (size_t i = 0; i < derived.positions.size(); ++i) {
    for (const Stream *stream : {&source_registers, &positions_only}) {
      HS_EXPECT_NEAR(derived.positions[i].x, stream->positions[i].x,
                     POLICY_TOL);
      HS_EXPECT_NEAR(derived.positions[i].y, stream->positions[i].y,
                     POLICY_TOL);
      HS_EXPECT_NEAR(derived.positions[i].z, stream->positions[i].z,
                     POLICY_TOL);
    }
    HS_EXPECT_NEAR(positions_only.positions[i].x,
                   rebuilt_positions_only.positions[i].x, POLICY_TOL);
    HS_EXPECT_NEAR(positions_only.positions[i].y,
                   rebuilt_positions_only.positions[i].y, POLICY_TOL);
    HS_EXPECT_NEAR(positions_only.positions[i].z,
                   rebuilt_positions_only.positions[i].z, POLICY_TOL);
    derived_differences +=
        derived.registers[i][0] != source_registers.registers[i][0] ||
        derived.registers[i][1] != source_registers.registers[i][1];
    HS_EXPECT_EQ(derived.registers[i][2], source_registers.registers[i][2]);
    source_differences += source_registers.registers[i][0] != 0 ||
                          source_registers.registers[i][1] != 0 ||
                          source_registers.registers[i][2] != 0;
    HS_EXPECT_EQ(positions_only.registers[i][0], 0u);
    HS_EXPECT_EQ(positions_only.registers[i][1], 0u);
    HS_EXPECT_EQ(positions_only.registers[i][2], 0u);
  }
  HS_EXPECT_GT(derived_differences, (size_t)0);
  HS_EXPECT_GT(source_differences, (size_t)0);
}
