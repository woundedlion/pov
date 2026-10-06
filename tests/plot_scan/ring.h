/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_plot_scan.h.

// ============================================================================
// Plot::Ring::sample
// ============================================================================

/**
 * @brief Verifies Ring::sample emits N unit-length fragments plus one closing
 *        overlap fragment, with v0 progress running 0..1 and the closing
 *        fragment coinciding with sample 0.
 */
inline void test_ring_sample_unit_length_and_progress() {
  ScratchScope sc(plot_arena());
  Fragments points;
  const int N = 32;
  points.bind(plot_arena(), N + 2);

  math::Basis b =
      math::make_basis(math::Quaternion(1, 0, 0, 0), math::Vector(0, 1, 0));
  Plot::Ring::sample(points, b, 0.5f, N, 0.0f);

  // N samples + 1 manual-close overlap fragment.
  HS_EXPECT_SIZE_OR_RETURN(points, (size_t)(N + 1));

  for (size_t i = 0; i < points.size(); ++i) {
    HS_EXPECT_NEAR(points[i].pos.length(), 1.0f, 1e-3f);
  }
  HS_EXPECT_NEAR(points[0].v0, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(points[points.size() - 1].v0, 1.0f, 1e-6f);

  HS_EXPECT_NEAR(points.back().pos.x, points[0].pos.x, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.y, points[0].pos.y, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.z, points[0].pos.z, 1e-3f);
}

/**
 * @brief Reconstructs a ring's W control vertices with libm cos/sin, bypassing
 *        the TrigLUT angle-addition identity Plot::Ring builds them from.
 * @param b Ring basis, as handed to Plot::Ring::sample.
 * @param radius Ring radius in [0,2], in hemisphere units as sample().
 * @param phase Angular offset added to every step.
 * @param W Number of control vertices (the close vertex is not emitted).
 * @return The W expected unit positions, in sample order.
 */
inline std::vector<math::Vector>
ring_vertices_direct(const math::Basis &b, float radius, float phase, int W) {
  auto res = math::get_antipode(b, radius);
  const math::Basis &wb = res.first;
  const float theta_eq = res.second * (math::PI_F / 2.0f);
  const float r_val = sinf(theta_eq);
  const float d_val = cosf(theta_eq);
  const float step = 2.0f * math::PI_F / W;

  std::vector<math::Vector> expected;
  expected.reserve(static_cast<size_t>(W));
  for (int i = 0; i < W; ++i) {
    const float t = i * step + phase;
    const math::Vector u_temp = (wb.u * cosf(t)) + (wb.w * sinf(t));
    expected.push_back(((wb.v * d_val) + (u_temp * r_val)).normalized());
  }
  return expected;
}

/** @brief sin of a ring's polar radius: the arc-length scale of its v1. */
inline float ring_arc_scale(const math::Basis &b, float radius) {
  return sinf(math::get_antipode(b, radius).second * (math::PI_F / 2.0f));
}

/**
 * @brief Verifies Ring::sample<W,H> built from the TrigLUT angle-addition
 *        identity matches a direct cos/sin(theta+phase) construction of the same
 *        ring, in both position and the analytic arc-length register.
 */
inline void test_ring_sample_lut_matches_direct() {
  constexpr int W = 64;
  constexpr int H = 64;

  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), W + 2);

  math::Basis b =
      math::make_basis(math::Quaternion(1, 0, 0, 0), math::Vector(0, 1, 0));
  const float radius = 0.5f;
  const float phase = 0.7f;
  Plot::Ring::sample<W, H>(points, b, radius, phase);

  HS_EXPECT_SIZE_OR_RETURN(points, (size_t)(W + 1));

  const std::vector<math::Vector> expected =
      ring_vertices_direct(b, radius, phase, W);
  const float r_val = ring_arc_scale(b, radius);
  const float step = 2.0f * math::PI_F / W;

  for (int i = 0; i < W; ++i) {
    HS_EXPECT_NEAR(points[i].pos.x, expected[i].x, 2e-3f);
    HS_EXPECT_NEAR(points[i].pos.y, expected[i].y, 2e-3f);
    HS_EXPECT_NEAR(points[i].pos.z, expected[i].z, 2e-3f);
    HS_EXPECT_NEAR(points[i].pos.length(), 1.0f, 1e-3f);
    HS_EXPECT_NEAR(points[i].v1, i * step * r_val, 2e-3f);
  }

  HS_EXPECT_NEAR(points.back().pos.x, points[0].pos.x, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.y, points[0].pos.y, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.z, points[0].pos.z, 1e-3f);
  HS_EXPECT_NEAR(points.back().v0, 1.0f, 1e-6f);

  // Close vertex carries the full-perimeter arc length: the ring's true
  // circumference 2*pi*sin(theta_eq).
  HS_EXPECT_NEAR(points.back().v1, 2.0f * math::PI_F * r_val, 2e-3f);
}

/**
 * @brief Verifies a ring drawn on the strided LUT grid stays an unbroken curve
 *        tracking the full-W control grid, for radii where the stride thins.
 * @details The lit set stays one 8-connected component (columns wrap) and
 *          every pixel the full-W runtime sampler lights has a lit neighbour
 *          within one pixel. Coverage is not compared pixel-for-pixel: the full-W
 *          grid's accumulated splat tails widen its stroke by a pixel.
 */
inline void test_ring_draw_stride_tracks_full_grid() {
  constexpr int W = 96, H = 48;
  auto shade = [](const math::Vector &, Fragment &f) {
    f.color = Color4(Pixel(65535, 65535, 65535), 1.0f);
  };

  const math::Vector centers[] = {math::Vector(0, 1, 0),
                                  math::Vector(0.3f, 0.8f, 0.5f),
                                  math::Vector(1, 0, 0)};
  for (const math::Vector &center : centers) {
    math::Basis b = math::make_basis(math::Quaternion(1, 0, 0, 0), center);
    for (float radius : {0.06f, 0.15f, 0.3f}) {
      HS_EXPECT_GT((Plot::ring_lut_stride<W>(sinf(radius * math::PI_F / 2.0f))),
                   1);

      auto render = [&](auto &&draw_one) {
        std::vector<uint8_t> mask(static_cast<size_t>(W) * H, 0);
        hs_test::StubEffect fx(W, H);
        Pipeline<W, H> filters;
        {
          Canvas c(fx);
          draw_one(filters, c);
        }
        fx.advance_display();
        for (int y = 0; y < H; ++y)
          for (int x = 0; x < W; ++x) {
            Pixel p = fx.get_pixel(x, y);
            mask[static_cast<size_t>(y) * W + x] = (p.r | p.g | p.b) ? 1 : 0;
          }
        return mask;
      };

      std::vector<uint8_t> ref =
          render([&](Pipeline<W, H> &filters, Canvas &c) {
            Plot::draw_fragments<W, H>(
                filters, c, {}, shade, {.capacity = W + 2, .omit_end = true},
                [&](Fragments &pts) {
                  Plot::Ring::sample(pts, b, radius, W, 0.0f);
                });
          });
      std::vector<uint8_t> cur =
          render([&](Pipeline<W, H> &filters, Canvas &c) {
            Plot::Ring::draw<W, H>(filters, c, b, radius, shade);
          });

      int lit = 0, seed = -1;
      for (int i = 0; i < W * H; ++i)
        if (cur[static_cast<size_t>(i)]) {
          ++lit;
          if (seed < 0)
            seed = i;
        }
      HS_EXPECT_GT(lit, 0);
      if (seed < 0)
        continue;

      // One 8-connected component: a beaded ring splits into many.
      std::vector<uint8_t> seen(static_cast<size_t>(W) * H, 0);
      std::vector<int> stack{seed};
      seen[static_cast<size_t>(seed)] = 1;
      int reached = 0;
      while (!stack.empty()) {
        int i = stack.back();
        stack.pop_back();
        ++reached;
        int y = i / W, x = i % W;
        for (int dy = -1; dy <= 1; ++dy)
          for (int dx = -1; dx <= 1; ++dx) {
            int ny = y + dy;
            if (ny < 0 || ny >= H)
              continue;
            int nx = ((x + dx) % W + W) % W;
            size_t j = static_cast<size_t>(ny) * W + nx;
            if (cur[j] && !seen[j]) {
              seen[j] = 1;
              stack.push_back(static_cast<int>(j));
            }
          }
      }
      HS_EXPECT_EQ(reached, lit);

      // The drawn curve stays on the reference ring, to within one pixel.
      int drifted = 0;
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
          if (!ref[static_cast<size_t>(y) * W + x])
            continue;
          bool near = false;
          for (int dy = -1; dy <= 1 && !near; ++dy) {
            int ny = y + dy;
            if (ny < 0 || ny >= H)
              continue;
            for (int dx = -1; dx <= 1 && !near; ++dx) {
              int nx = ((x + dx) % W + W) % W;
              near = cur[static_cast<size_t>(ny) * W + nx] != 0;
            }
          }
          if (!near)
            ++drifted;
        }
      HS_EXPECT_EQ(drifted, 0);
    }
  }
}

/** @brief Verifies a stroke primitive draws through a direct anti-alias sink. */
inline void test_ring_draw_accepts_direct_sink() {
  constexpr int W = 64, H = 32;
  hs_test::StubEffect fx(W, H);
  Filter::Screen::DirectAntiAliasSink<W, H> sink;
  {
    Canvas canvas(fx);
    sink.prepare(canvas);
    const math::Basis basis =
        math::make_basis(math::Quaternion(), math::Vector(0, 1, 0));
    auto shade = [](const math::Vector &, Fragment &f) {
      f.color = Color4(Pixel(65535, 65535, 65535), 1.0f);
    };
    Plot::Ring::draw<W, H>(sink, canvas, basis, 0.5f, shade);
  }
  fx.advance_display();

  size_t lit = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      if (!is_black(fx.get_pixel(x, y)))
        ++lit;
  HS_EXPECT_GT(lit, size_t{0});
}
