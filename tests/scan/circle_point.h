/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Circle/Point cap drawing, raster epilogues, and replicated clip parity
//
// A radius-0 SDF::Ring draws a spherical cap of angular radius `thickness`
// centred on the basis axis, with quintic coverage from 1 at the centre to 0 at
// the rim.
// ============================================================================

/**
 * @brief Analytic coverage of a radius-0 ring at one direction.
 * @param v Unit direction of the pixel centre.
 * @param axis Cap centre (the ring's basis axis).
 * @param thickness Cap angular radius in radians.
 * @return Stroke coverage in [0, 1]; 0 at or beyond the rim.
 * @details Uses double-precision trigonometry and the analytic quintic.
 */
inline float cap_coverage(const math::Vector &v, const math::Vector &axis,
                          float thickness) {
  if (thickness <= 0.0f)
    return 0.0f;
  const double d = static_cast<double>(v.x) * axis.x +
                   static_cast<double>(v.y) * axis.y +
                   static_cast<double>(v.z) * axis.z;
  const double angle = std::acos(std::clamp(d, -1.0, 1.0));
  const double t = std::clamp(1.0 - angle / thickness, 0.0, 1.0);
  return static_cast<float>(t * t * t * (10.0 + t * (6.0 * t - 15.0)));
}

/**
 * @brief Verifies Point::draw paints exactly the analytic spherical cap.
 * @details Pole-centred, equatorial and oblique axes each reproduce cap_coverage
 * per pixel: every direction with coverage above 0.02 lit, every uncovered one
 * black, and the plotted channel proportional to the coverage.
 */
inline void test_point_draws_the_analytic_cap() {
  constexpr int W = 96, H = 64;
  constexpr uint16_t LEVEL = 60000;
  constexpr float THICKNESS = 0.35f;
  math::TrigLUT<W, H>::init();

  const math::Vector axes[] = {math::Y_AXIS, math::Vector(1.0f, 0.0f, 0.0f),
                               math::Vector(0.3f, -0.6f, 0.74f).normalized()};
  for (const math::Vector &axis : axes) {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas c(fx);
      Scan::Point::draw<W, H>(
          pipe, c, axis, THICKNESS, [](const math::Vector &, Fragment &f) {
            f.color = Color4(Pixel(LEVEL, LEVEL, LEVEL), 1.0f);
          });
    }
    fx.advance_display();

    // make_basis can reorient, so the cap centre is read back from the basis.
    const math::Basis basis = math::make_basis(math::Quaternion(), axis);
    size_t lit = 0, covered = 0;
    float worst_value_error = 0.0f;
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        const math::Vector v = math::pixel_to_vector<W, H>(x, y);
        const float coverage = cap_coverage(v, basis.v, THICKNESS);
        const bool is_lit = !is_black(fx.get_pixel(x, y));
        if (is_lit)
          ++lit;
        if (coverage > 0.0f)
          ++covered;
        // Coverage 0 is a hard cut in the shape, so an outside pixel is black.
        if (coverage == 0.0f)
          HS_EXPECT_FALSE(is_lit);
        // Inside the cut, the plotted channel is the level scaled by coverage.
        if (coverage > 0.02f) {
          HS_EXPECT_TRUE(is_lit);
          const float want = static_cast<float>(LEVEL) * coverage;
          worst_value_error = std::max(
              worst_value_error,
              std::fabs(static_cast<float>(fx.get_pixel(x, y).r) - want));
        }
      }
    HS_EXPECT_GT(covered, (size_t)0);
    HS_EXPECT_GT(lit, (size_t)0);
    // 60000 * max(quintic derivative) * acos error / thickness, plus rounding.
    HS_EXPECT_LT(worst_value_error, 18.0f);
  }
}

/**
 * @brief Verifies a pole-centred cap really exercises the full-row-scan regime.
 * @details Point::draw at +Y builds a ring whose axis has no horizontal
 * projection, so get_horizontal_intervals refuses every row.
 */
inline void test_pole_centred_cap_takes_the_full_row_scan() {
  constexpr int W = 96, H = 64;
  math::TrigLUT<W, H>::init();
  const math::Basis pole = math::make_basis(math::Quaternion(), math::Y_AXIS);
  const SDF::Ring degenerate(pole, 0.0f, 0.35f);

  const auto rows = degenerate.get_vertical_bounds<H>();
  HS_EXPECT_EQ(rows.y_min, 0);
  HS_EXPECT_GT(rows.y_max, 0);
  int full_scan_rows = 0;
  for (int y = rows.y_min; y <= rows.y_max; ++y)
    if (degenerate.needs_full_row_scan(math::TrigLUT<W, H>::sin_phi[y]))
      ++full_scan_rows;
  HS_EXPECT_EQ(full_scan_rows, rows.y_max - rows.y_min + 1);

  // An equatorial cap of the same size answers with intervals on most rows.
  const math::Basis equator =
      math::make_basis(math::Quaternion(), math::Vector(1.0f, 0.0f, 0.0f));
  const SDF::Ring ordinary(equator, 0.0f, 0.35f);
  const auto eq_rows = ordinary.get_vertical_bounds<H>();
  int interval_rows = 0;
  for (int y = eq_rows.y_min; y <= eq_rows.y_max; ++y)
    if (!ordinary.needs_full_row_scan(math::TrigLUT<W, H>::sin_phi[y]))
      ++interval_rows;
  HS_EXPECT_GT(interval_rows, 0);
}

/** @brief Circle and point scans retain exact pixel centers. */
inline void test_circle_and_point_keep_exact_pixel_centers() {
  constexpr int W = 96, H = 64;
  math::TrigLUT<W, H>::init();
  for (int y : {1, 17, 31, 48, 62})
    for (int x : {1, 7, 13, 23, 47, 71, 95})
      for (bool circle : {false, true}) {
        const auto center = math::pixel_to_vector<W, H>(x, y);
        const auto basis = math::make_basis(math::Quaternion(), center);
        hs_test::StubEffect fx(W, H);
        Pipeline<W, H> pipe;
        {
          Canvas canvas(fx);
          auto shader = [](const math::Vector &, Fragment &f) {
            f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
          };
          if (circle)
            Scan::Circle::draw<W, H>(pipe, canvas, basis, 0.2f, shader);
          else
            Scan::Point::draw<W, H>(pipe, canvas, center, 0.3f, shader);
        }
        fx.advance_display();
        HS_EXPECT_GT(fx.get_pixel(x, y).r, 59900);
      }
}

/**
 * @brief Verifies the Circle and Point wrappers are their documented rings.
 * @details Circle is a radius-0 ring whose stroke half-width is radius * pi/2;
 * Point is a radius-0 ring of the given thickness. Both rasterize
 * bit-identically to that ring.
 */
inline void test_circle_and_point_match_their_rings() {
  constexpr int W = 96, H = 64;
  Pipeline<W, H> pipe;
  auto shader = [](const math::Vector &p, Fragment &f) {
    // Position-dependent so a mismatched basis or phase shows up as color, not
    // only as coverage.
    f.color =
        Color4(Pixel(static_cast<uint16_t>(30000.0f + 20000.0f * p.x),
                     static_cast<uint16_t>(30000.0f + 20000.0f * p.y), 50000),
               1.0f);
  };

  const math::Basis basis =
      math::make_basis(math::Quaternion(), math::Vector(0.2f, 0.8f, -0.5f));
  std::vector<Pixel> circle_frame, ring_frame;
  {
    hs_test::StubEffect fx(W, H);
    {
      Canvas c(fx);
      Scan::Circle::draw<W, H>(pipe, c, basis, 0.4f, shader);
    }
    fx.advance_display();
    capture_frame<W, H>(fx, circle_frame);
  }
  {
    hs_test::StubEffect fx(W, H);
    {
      Canvas c(fx);
      Scan::Ring::draw<W, H>(pipe, c, basis, 0.0f, 0.4f * (math::PI_F / 2.0f),
                             shader);
    }
    fx.advance_display();
    capture_frame<W, H>(fx, ring_frame);
  }
  size_t circle_lit = 0, circle_diff = 0;
  for (size_t i = 0; i < circle_frame.size(); ++i) {
    if (!is_black(circle_frame[i]))
      ++circle_lit;
    if (!(circle_frame[i] == ring_frame[i]))
      ++circle_diff;
  }
  HS_EXPECT_GT(circle_lit, (size_t)0);
  HS_EXPECT_EQ(circle_diff, (size_t)0);

  const math::Vector center(-0.4f, 0.5f, 0.766f);
  std::vector<Pixel> point_frame, point_ring_frame;
  {
    hs_test::StubEffect fx(W, H);
    {
      Canvas c(fx);
      Scan::Point::draw<W, H>(pipe, c, center, 0.3f, shader);
    }
    fx.advance_display();
    capture_frame<W, H>(fx, point_frame);
  }
  {
    hs_test::StubEffect fx(W, H);
    const math::Basis point_basis =
        math::make_basis(math::Quaternion(), center);
    {
      Canvas c(fx);
      Scan::Ring::draw<W, H>(pipe, c, point_basis, 0.0f, 0.3f, shader);
    }
    fx.advance_display();
    capture_frame<W, H>(fx, point_ring_frame);
  }
  size_t point_lit = 0, point_diff = 0;
  for (size_t i = 0; i < point_frame.size(); ++i) {
    if (!is_black(point_frame[i]))
      ++point_lit;
    if (!(point_frame[i] == point_ring_frame[i]))
      ++point_diff;
  }
  HS_EXPECT_GT(point_lit, (size_t)0);
  HS_EXPECT_EQ(point_diff, (size_t)0);
}

/**
 * @brief Verifies a Circle's painted extent tracks its radius argument.
 * @details radius maps to a cap of angular half-width radius * pi/2, so the
 * lit area grows with radius and the farthest lit direction sits at that angle.
 */
inline void test_circle_extent_follows_its_radius() {
  constexpr int W = 96, H = 64;
  const math::Basis basis =
      math::make_basis(math::Quaternion(), math::Vector(0.0f, 0.0f, 1.0f));
  size_t previous_lit = 0;
  for (float radius : {0.15f, 0.3f, 0.6f}) {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas c(fx);
      Scan::Circle::draw<W, H>(
          pipe, c, basis, radius, [](const math::Vector &, Fragment &f) {
            f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
          });
    }
    fx.advance_display();

    const float rim = radius * (math::PI_F / 2.0f);
    size_t lit = 0;
    float widest = 0.0f;
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        if (!is_black(fx.get_pixel(x, y))) {
          ++lit;
          const math::Vector v = math::pixel_to_vector<W, H>(x, y);
          widest = std::max(widest, math::fast_acos(hs::clamp(
                                        math::dot(v, basis.v), -1.0f, 1.0f)));
        }
    HS_EXPECT_GT(lit, previous_lit);
    // No lit pixel past the rim, and the cap is sampled close enough to it
    // that the widest lit direction is within two rows of the rim.
    HS_EXPECT_LT(widest, rim);
    HS_EXPECT_GT(widest, rim - 2.0f * (math::PI_F / (H - 1)));
    previous_lit = lit;
  }
}

/** @brief Captures the coordinates, color, age, and alpha emitted by scans. */
struct EpilogueCapture {
  struct Point {
    int x, y;
    Pixel color;
    float age, alpha;
  };
  std::vector<Point> points;
  void plot(Canvas &, int x, int y, const Pixel &color, float age,
            float alpha) {
    points.push_back({x, y, color, age, alpha});
  }
};

/** @brief Every scan tail preserves shader output and multiplies coverage once. */
inline void test_scan_epilogue_contract() {
  constexpr int W = 96, H = 48;
  const ScopedPoleLod lod(0.0f);
  const Pixel color(45000, 17000, 32000);
  constexpr float OPACITY = 0.37f;
  const auto basis = math::make_basis(math::Quaternion(), math::X_AXIS);
  SDF::Ring ring(basis, 1.0f, 0.18f);
  float knots[9]{};
  SDF::DistortedRing distorted(basis, 1.0f, 0.18f, knots, 8, 0.0f, nullptr);
  const int8_t slots[1] = {0};
  Scan::DistortedRingStack::CandidateTable<W, H> stack_table;
  math::Vector vertices[4];
  const uint16_t indices[4] = {0, 1, 2, 3};
  for (int i = 0; i < 4; ++i) {
    const float angle = i * math::PI_F * 0.5f;
    vertices[i] = (basis.v * cosf(0.7f) +
                   (basis.u * cosf(angle) + basis.w * sinf(angle)) * sinf(0.7f))
                      .normalized();
  }
  for (int path = 0; path < 5; ++path) {
    HS_CONTEXT("epilogue path", path);
    StubEffect effect(W, H);
    EpilogueCapture capture;
    std::vector<EpilogueCapture::Point> expected;
    size_t shaded = 0;
    {
      Canvas canvas(effect);
      SDF::FaceScratchBuffer scratch;
      SDF::Face face(vertices, indices, scratch, H + hs::H_OFFSET, H,
                     &canvas.clip());
      const auto face_coverage = [&](const math::Vector &point) {
        const float d = SDF::distance_of(face, point).dist;
        return Scan::solid_coverage(d, math::TWO_PI_F / W);
      };
      auto shader = [&](const math::Vector &point, Fragment &fragment) {
        HS_EXPECT_PIXEL(fragment.color.color, 0, 0, 0);
        HS_EXPECT_EQ(fragment.color.alpha, 0.0f);
        HS_EXPECT_EQ(fragment.age, 0.0f);
        const size_t variant = shaded++ % 3;
        if (variant == 0)
          return;
        fragment.color =
            Color4(color, variant == 1 ? Scan::MIN_ALPHA : OPACITY);
        fragment.age = 7.0f;
        if (variant == 2) {
          const float coverage = path == 2 ? face_coverage(point) : fragment.v2;
          expected.push_back({0, 0, color, 7.0f, OPACITY * coverage});
        }
      };
      if (path == 0)
        Scan::rasterize<W, H, false>(capture, canvas, ring, shader);
      else if (path == 1)
        Scan::rasterize_solid<W, H>(capture, canvas, face,
                                    Color4(color, OPACITY));
      else if (path == 2)
        Scan::rasterize_face<W, H>(capture, canvas, face, shader);
      else if (path == 3)
        Scan::DistortedRingStack::draw<W, H>(
            capture, canvas, 1, &distorted, slots, 1, stack_table,
            [&](int, const math::Vector &point, Fragment &fragment) {
              shader(point, fragment);
            });
      else
        Scan::RingGroup::draw<W, H>(
            capture, canvas, &ring, 1,
            [&](int, const math::Vector &point, Fragment &fragment) {
              shader(point, fragment);
            });
      if (path == 1) {
        EpilogueCapture reference;
        Scan::rasterize<W, H, false>(
            reference, canvas, face,
            [&](const math::Vector &, Fragment &fragment) {
              fragment.color = Color4(color, OPACITY);
            });
        expected = std::move(reference.points);
      }
    }
    effect.advance_display();
    HS_EXPECT_GT(capture.points.size(), size_t{20});
    HS_EXPECT_EQ(capture.points.size(), expected.size());
    for (size_t index = 0;
         index < std::min(capture.points.size(), expected.size()); ++index) {
      const auto &actual = capture.points[index];
      const auto &want = expected[index];
      if (path == 1) {
        HS_EXPECT_EQ(actual.x, want.x);
        HS_EXPECT_EQ(actual.y, want.y);
      }
      HS_EXPECT_PIXEL(actual.color, want.color.r, want.color.g, want.color.b);
      HS_EXPECT_EQ(actual.age, want.age);
      HS_EXPECT_NEAR(actual.alpha, want.alpha, 1e-5f);
    }
  }
}

/** @brief Replicated source pixels survive destination-band clipping. */
inline void test_replicated_clip_matches_full_frame() {
  constexpr int W = 96, H = 48;
  const Color4 COLOR(Pixel(50000, 30000, 10000), 0.8f);
  const math::Vector NORMAL(0, 0, -1);
  const math::Basis BASIS = math::make_basis(math::Quaternion(), NORMAL);
  for (int count : {2, 3}) {
    for (bool solid : {false, true}) {
      std::vector<Pixel> expected(W * H);
      for (int band = -1; band < 5; ++band) {
        StubEffect effect(W, H);
        effect.set_margin(0);
        const int X0 = band < 0   ? 0
                       : band < 2 ? band * W / 2
                                  : (band - 2) * W / 3;
        const int X1 = band < 0   ? W
                       : band < 2 ? (band + 1) * W / 2
                                  : (band - 1) * W / 3;
        effect.set_clip(0, H, X0, X1);
        Pipeline<W, H, Filter::World::Replicate<W>> pipeline(count);
        {
          Canvas canvas(effect);
          auto shader = [&](const math::Vector &, Fragment &fragment) {
            fragment.color = COLOR;
          };
          if (solid) {
            SDF::PlanarPolygon shape(BASIS, 0.2f, 5, 0.0f);
            Scan::rasterize_solid<W, H>(pipeline, canvas, shape, COLOR);
          } else {
            Scan::Circle::draw<W, H>(pipeline, canvas, NORMAL, 0.1f, shader);
          }
        }
        effect.advance_display();
        int lit = 0;
        for (int y = 0; y < H; ++y) {
          for (int x = 0; x < W; ++x) {
            const Pixel actual = effect.get_pixel(x, y);
            if (band < 0) {
              expected[y * W + x] = actual;
              if (x < W / 2 && (actual.r || actual.g || actual.b))
                ++lit;
            } else {
              const Pixel want =
                  x >= X0 && x < X1 ? expected[y * W + x] : Pixel(0, 0, 0);
              HS_EXPECT_PIXEL(actual, want.r, want.g, want.b);
            }
          }
        }
        if (band < 0)
          HS_EXPECT_GT(lit, 0);
      }
    }
  }
}
