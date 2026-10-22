/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Plot::PointRoute: closed-form column-time points against a per-column oracle.
// ---------------------------------------------------------------------------

/** @brief A point stepped along a great circle under a slerped rotation, the
 *         motion of a particle drawn through an orientation. */
struct PointRoutePath {
  math::Vector from;         ///< Object-space position at t = 0.
  math::Vector axis;         ///< Unit axis of the object-space step.
  float angle;               ///< Object-space step angle, radians.
  math::Quaternion previous; ///< Orientation at t = 0.
  math::Quaternion current;  ///< Orientation at t = 1.

  /** @brief World position at a frame fraction.
   * @param t Frame fraction in [0, 1].
   * @return Unit world position. */
  math::Vector at(float t) const {
    const math::Vector p =
        math::rotate(from, math::make_rotation(axis, angle * t));
    return math::rotate(p, math::slerp(previous, current, t));
  }
};

/** @brief Builds a path ending at a pixel.
 * @param x,y Pixel the path reaches at t = 1.
 * @param heading Direction of the object-space step, radians from east.
 * @param angle Object-space step angle, radians.
 * @param spin Orientation turn over the frame, radians.
 * @param spin_axis Unit axis of the orientation turn.
 * @return The path. */
inline PointRoutePath point_route_path(float x, float y, float heading,
                                       float angle, float spin,
                                       const math::Vector &spin_axis) {
  constexpr int W = 288, H = 144;
  PointRoutePath path;
  path.previous = math::make_rotation(spin_axis, 0.3f);
  path.current = math::make_rotation(spin_axis, 0.3f + spin);
  const math::Vector world = math::pixel_to_vector<W, H>(x, y);
  // Increasing column and increasing row at the world position.
  const math::Vector east = math::Vector(-world.z, 0.0f, world.x).normalized();
  const math::Vector south = math::cross(east, world);
  const math::Vector direction = math::rotate(
      east * cosf(heading) + south * sinf(heading), path.current.conjugate());
  const math::Vector end = math::rotate(world, path.current.conjugate());
  path.axis = math::cross(end, direction).normalized();
  path.angle = angle;
  path.from = math::rotate(end, math::make_rotation(path.axis, -angle));
  return path;
}

/** @brief Per-column oracle: each column splats the point frozen at the
 *         fraction that column is lit.
 * @param path Point motion.
 * @param alpha Point opacity.
 * @return Weight per pixel, row-major. */
inline std::vector<float> point_route_oracle(const PointRoutePath &path,
                                             float alpha) {
  constexpr int W = 288, H = 144;
  std::vector<float> weights(W * H, 0.0f);
  for (int column = 0; column < W; ++column) {
    const math::PixelCoords p = math::vector_to_pixel<W, H>(
        path.at(Temporal::sweep_fraction(column, W)));
    const Filter::Screen::SplatTaps t =
        Filter::Screen::splat_taps<W, H>(p.x, p.y);
    const bool left = t.x0 == column;
    if (!left && t.x1 != column)
      continue;
    if (t.y0_physical)
      weights[t.y0 * W + column] += alpha * (left ? t.v00 : t.v10);
    if (t.y1_physical)
      weights[t.y1 * W + column] += alpha * (left ? t.v01 : t.v11);
  }
  return weights;
}

/** @brief Arm passes of a path, by exhaustive search of its true column.
 * @param path Point motion.
 * @return Integer crossings of `column / (W / 2) - t` over the frame. */
inline int point_route_true_passes(const PointRoutePath &path) {
  constexpr int W = 288, N = 4000;
  auto relative = [&](float t, float &unwrapped, float &previous) {
    const float column = math::vector_to_pixel<W, 144>(path.at(t)).x;
    unwrapped += std::remainder(column - previous, float(W));
    previous = column;
    return unwrapped / (W / 2.0f) - t;
  };
  float previous = math::vector_to_pixel<W, 144>(path.at(0.0f)).x;
  float unwrapped = previous;
  const float r0 = relative(0.0f, unwrapped, previous);
  float r = r0;
  for (int i = 1; i <= N; ++i)
    r = relative(float(i) / N, unwrapped, previous);
  return static_cast<int>(floorf(r0) - floorf(r));
}

/** @brief Plots one path through the route.
 * @param canvas Destination canvas; cleared first.
 * @param path Point motion.
 * @param alpha Point opacity.
 * @param motion Motion bound handed to the route.
 * @return Frame pixels. */
inline std::vector<Pixel> point_route_render(Canvas &canvas,
                                             const PointRoutePath &path,
                                             float alpha, float motion) {
  constexpr int W = 288, H = 144;
  std::fill_n(canvas.data(), W * H, Pixel{});
  Filter::Screen::DirectAntiAliasSink<W, H> sink;
  sink.prepare(canvas);
  Plot::PointRoute<W, H> route;
  route.prepare(canvas, motion);
  const math::Vector end = path.at(1.0f);
  route.plot(
      sink, canvas, end, [&] { return path.at(0.0f); },
      [&] { return Color4(Pixel(65535, 65535, 65535), alpha); });
  return std::vector<Pixel>(canvas.data(), canvas.data() + W * H);
}

/** @brief Discrepancy between a route render and the oracle. */
struct PointRouteError {
  float weight = 0.0f; ///< Largest per-pixel weight difference.
  float arc = 0.0f;    ///< Largest centroid separation per arm, in columns.
  float mass = 0.0f;   ///< Largest total-weight difference per arm.
};

/** @brief Compares a route render with the oracle, pixel by pixel and, per
 *         arm, by total weight and by the arc between weighted centroids.
 * @param got Route output, white on black.
 * @param want Oracle weights.
 * @return The discrepancy. */
inline PointRouteError point_route_error(const std::vector<Pixel> &got,
                                         const std::vector<float> &want) {
  constexpr int W = 288, H = 144;
  PointRouteError error;
  std::array<math::Vector, 4> centroid{};
  std::array<float, 4> mass{};
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const size_t i = size_t(y) * W + x;
      const float g = got[i].r / 65535.0f;
      error.weight = fmaxf(error.weight, fabsf(g - want[i]));
      const math::Vector v = math::pixel_to_vector<W, H>(x, y);
      const int arm = x < W / 2 ? 0 : 2;
      centroid[arm] = centroid[arm] + v * g;
      mass[arm] += g;
      centroid[arm + 1] = centroid[arm + 1] + v * want[i];
      mass[arm + 1] += want[i];
    }
  for (int arm : {0, 2}) {
    error.mass = fmaxf(error.mass, fabsf(mass[arm] - mass[arm + 1]));
    if (mass[arm] > 0.05f && mass[arm + 1] > 0.05f)
      error.arc =
          fmaxf(error.arc, math::angle_between(centroid[arm].normalized(),
                                               centroid[arm + 1].normalized()) /
                               math::RADIANS_PER_COLUMN<W>);
  }
  return error;
}

/** @brief Points land in the column lit when they pass it, as the per-column
 *         oracle shows them, across the frame, the arm cut and the poles.
 * @details Away from the poles the render matches the oracle pixel by pixel.
 * Near a pole, columns are a small fraction of a pixel of arc, so the
 * comparison is of each arm's weighted centroid, in arc, and total weight;
 * a point moving with or against the sweep at a speed comparable to the arm's
 * is lit by more or fewer than two columns of the oracle there. */
inline void test_point_route_matches_column_oracle() {
  constexpr int W = 288;
  hs_test::StubEffect effect(W, 144);
  Canvas canvas(effect);
  const math::Vector tilted = math::Vector(0.3f, 1.0f, -0.2f).normalized();
  std::array<int, 3> passes{};
  PointRouteError interior, rapid, cut, polar;
  auto check = [&](const PointRoutePath &path, PointRouteError &worst) {
    const auto want = point_route_oracle(path, 1.0f);
    const auto got = point_route_render(
        canvas, path, 1.0f, (path.at(1.0f) - path.at(0.0f)).magnitude());
    const PointRouteError e = point_route_error(got, want);
    worst.weight = fmaxf(worst.weight, e.weight);
    worst.arc = fmaxf(worst.arc, e.arc);
    worst.mass = fmaxf(worst.mass, e.mass);
  };
  uint32_t state = 17u;
  auto uniform = [&](float lo, float hi) {
    state = state * 1664525u + 1013904223u;
    return lo + (hi - lo) * float(state >> 8) * (1.0f / 16777216.0f);
  };
  for (int i = 0; i < 400; ++i) {
    const float x = uniform(0.0f, float(W));
    const float y = uniform(16.0f, 128.0f);
    const float heading = uniform(0.0f, math::TWO_PI_F);
    for (float spin : {0.0f, 0.02f}) {
      for (float angle : {0.0134f, 0.04f})
        check(point_route_path(x, y, heading, angle, spin, tilted), interior);
      check(point_route_path(x, y, heading, 0.18f, spin, tilted), rapid);
    }
  }
  // Motion along the azimuth across both arm cuts, with and against the
  // sweep: zero, one and two passes.
  for (float base : {143.0f, 287.0f})
    for (float offset = -2.5f; offset <= 2.5f; offset += 0.125f)
      for (float columns : {-2.5f, -1.0f, -0.3f, 0.3f, 1.0f, 2.5f}) {
        const float x = std::fmod(base + offset + float(W), float(W));
        const PointRoutePath path = point_route_path(
            x, 71.7f, 0.0f, columns * math::RADIANS_PER_COLUMN<W>, 0.0f,
            tilted);
        ++passes[std::min(point_route_true_passes(path), 2)];
        check(path, cut);
      }
  for (int i = 0; i < 400; ++i) {
    const float x = uniform(0.0f, float(W));
    const float row = uniform(0.5f, 8.0f);
    const float y = (i & 1) ? row : 143.0f - row;
    const float heading = uniform(0.0f, math::TWO_PI_F);
    for (float angle : {0.0134f, 0.04f})
      check(point_route_path(x, y, heading, angle, 0.01f, tilted), polar);
  }
  std::printf("  point route oracle (weight/arc/mass): interior %.4f/%.4f/%.4f "
              "rapid %.4f/%.4f/%.4f cut %.4f/%.4f/%.4f polar %.4f/%.4f/%.4f "
              "cut passes 0/1/2=%d/%d/%d\n",
              interior.weight, interior.arc, interior.mass, rapid.weight,
              rapid.arc, rapid.mass, cut.weight, cut.arc, cut.mass,
              polar.weight, polar.arc, polar.mass, passes[0], passes[1],
              passes[2]);
  HS_EXPECT_GT(passes[0], 0);
  HS_EXPECT_GT(passes[1], 0);
  HS_EXPECT_GT(passes[2], 0);
  HS_EXPECT_LE(interior.weight, 0.1f);
  HS_EXPECT_LE(interior.arc, 0.06f);
  HS_EXPECT_LE(interior.mass, 0.08f);
  HS_EXPECT_LE(rapid.weight, 0.2f);
  HS_EXPECT_LE(rapid.arc, 0.12f);
  HS_EXPECT_LE(rapid.mass, 0.08f);
  HS_EXPECT_LE(cut.weight, 0.04f);
  HS_EXPECT_LE(cut.mass, 0.04f);
  HS_EXPECT_LE(polar.arc, 0.2f);
  HS_EXPECT_LE(polar.mass, 1.0f);
}

/** @brief The clip gate drops only points that write no visible tap: a
 *         clipped render equals the unclipped one inside the clip, timed and
 *         still, and the still route equals a plain projected splat. */
inline void test_point_route_gate_keeps_visible_taps() {
  constexpr int W = 288, H = 144;
  const math::Vector tilted = math::Vector(0.3f, 1.0f, -0.2f).normalized();
  uint32_t state = 5u;
  auto uniform = [&](float lo, float hi) {
    state = state * 1664525u + 1013904223u;
    return lo + (hi - lo) * float(state >> 8) * (1.0f / 16777216.0f);
  };
  std::vector<PointRoutePath> paths;
  for (int i = 0; i < 3000; ++i)
    paths.push_back(
        point_route_path(uniform(0.0f, float(W)), uniform(0.0f, float(H - 1)),
                         uniform(0.0f, math::TWO_PI_F), uniform(0.0f, 0.06f),
                         uniform(0.0f, 0.03f), tilted));
  // Every path into one frame through the route, or through a plain
  // projected splat at t = 1.
  std::vector<uint8_t> visible(size_t(W) * H);
  auto render = [&](const std::array<int, 4> *clip, bool plain) {
    hs_test::StubEffect effect(W, H);
    if (clip)
      effect.set_clip((*clip)[0], (*clip)[1], (*clip)[2], (*clip)[3]);
    const ClipRegion &region = effect.clip();
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        visible[size_t(y) * W + x] =
            region.contains_y(y) && !region.x_clip().clipped(x);
    std::vector<Pixel> frame;
    {
      Canvas canvas(effect);
      std::fill_n(canvas.data(), W * H, Pixel{});
      Filter::Screen::DirectAntiAliasSink<W, H> sink;
      sink.prepare(canvas);
      Plot::PointRoute<W, H> route;
      for (const PointRoutePath &path : paths) {
        const math::Vector end = path.at(1.0f);
        const Pixel white(65535, 65535, 65535);
        if (plain) {
          sink.plot(canvas, end, white, 0.0f, 0.7f);
          continue;
        }
        route.prepare(canvas, (end - path.at(0.0f)).magnitude());
        route.plot(
            sink, canvas, end, [&] { return path.at(0.0f); },
            [&] { return Color4(white, 0.7f); });
      }
      frame.assign(canvas.data(), canvas.data() + W * H);
    }
    return frame;
  };
  const std::array<std::array<int, 4>, 6> clips = {{{36, 72, 0, 144},
                                                    {72, 108, 144, 288},
                                                    {0, 36, 144, 288},
                                                    {108, 144, 0, 144},
                                                    {50, 51, 100, 101},
                                                    {0, 144, 30, 250}}};
  int drawn = 0;
  for (bool timed : {true, false}) {
    Temporal::set_temporal_sampling_enabled(timed);
    const std::vector<Pixel> full = render(nullptr, false);
    for (const auto &clip : clips) {
      const std::vector<Pixel> got = render(&clip, false);
      for (size_t i = 0; i < got.size(); ++i) {
        drawn += visible[i] && full[i].r > 0;
        HS_EXPECT_EQ(got[i], visible[i] ? full[i] : Pixel{});
      }
      if (!timed) {
        const std::vector<Pixel> plain = render(&clip, true);
        for (size_t i = 0; i < got.size(); ++i)
          HS_EXPECT_EQ(got[i], plain[i]);
      }
    }
  }
  Temporal::set_temporal_sampling_enabled(true);
  HS_EXPECT_GT(drawn, 1000);
}

/** @brief A point whose color changes over the frame is shaded at each
 *         pass's own time, the time its column is lit. */
inline void test_point_route_shades_at_pass_time() {
  constexpr int W = 288, H = 144;
  const math::Vector tilted = math::Vector(0.3f, 1.0f, -0.2f).normalized();
  hs_test::StubEffect effect(W, H);
  Canvas canvas(effect);
  Filter::Screen::DirectAntiAliasSink<W, H> sink;
  sink.prepare(canvas);
  Plot::PointRoute<W, H> route;
  uint32_t state = 11u;
  auto uniform = [&](float lo, float hi) {
    state = state * 1664525u + 1013904223u;
    return lo + (hi - lo) * float(state >> 8) * (1.0f / 16777216.0f);
  };
  int passes = 0;
  for (int i = 0; i < 500; ++i) {
    const PointRoutePath path =
        point_route_path(uniform(0.0f, float(W)), uniform(24.0f, 120.0f),
                         uniform(0.0f, math::TWO_PI_F), 0.02f, 0.01f, tilted);
    const math::Vector end = path.at(1.0f);
    route.prepare(canvas, (end - path.at(0.0f)).magnitude());
    route.plot(
        sink, canvas, end, [&] { return path.at(0.0f); },
        [&](float t) {
          ++passes;
          const float lit = Temporal::column_fraction<W>(
              math::vector_to_pixel<W, H>(path.at(t)).x);
          // Passes reaching into the frame from an arm cut take its ends.
          const float gap = fabsf(lit - t);
          HS_EXPECT_LE(fminf(gap, 1.0f - gap), 0.02f);
          return Color4(Pixel(65535, 65535, 65535), 0.5f);
        });
  }
  HS_EXPECT_GT(passes, 400);
}

/** @brief A chained edge keeping a short lattice over the frame plots its
 *         lattice points before its end vertex; a longer one is left to the
 *         column walk. */
inline void test_steady_edge_plots_lattice_points() {
  constexpr int W = 288, H = 144;
  hs_test::StubEffect effect(W, H);
  Canvas canvas(effect);
  Filter::Screen::DirectAntiAliasSink<W, H> sink;
  sink.prepare(canvas);
  Plot::PointRoute<W, H> route;
  route.prepare(canvas, 0.02f);
  const auto turn = math::make_rotation(math::Y_AXIS, 0.01f);
  for (const float columns : {0.6f, 2.0f, 3.1f, 12.0f}) {
    HS_CONTEXT("columns x10", static_cast<int>(columns * 10.0f));
    const math::Vector ends[2] = {
        math::pixel_to_vector<W, H>(100.0f, 70.0f),
        math::pixel_to_vector<W, H>(100.0f + columns, 70.0f)};
    const math::Vector starts[2] = {math::rotate(ends[0], turn),
                                    math::rotate(ends[1], turn)};
    std::vector<float> fractions;
    const bool plotted = Plot::Mesh::plot_steady_edge<W, H>(
        route, sink, canvas, starts, ends, 4, [&](float s, float) {
          if (fractions.empty() || fractions.back() != s)
            fractions.push_back(s);
          return Color4(Pixel(65535, 65535, 65535), 0.5f);
        });
    Plot::GeodesicEdgeSampler arcs[2];
    const int lattice =
        Plot::ColumnEdgeFrame::steady_lattice<W, H>(starts, ends, 4, arcs);
    HS_EXPECT_EQ(plotted, lattice != 0);
    HS_EXPECT_EQ(fractions.size(), static_cast<size_t>(lattice));
    for (size_t k = 0; k < fractions.size(); ++k)
      HS_EXPECT_NEAR(fractions[k], float(k) / float(lattice), 1e-6f);
    if (columns > 10.0f)
      HS_EXPECT_FALSE(plotted);
    if (columns < 1.0f)
      HS_EXPECT_TRUE(plotted);
  }
}
