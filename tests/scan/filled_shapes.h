/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Scan::Star / PlanarPolygon / Flower — filled-shape placement oracle
// ============================================================================

/**
 * @brief Reports whether any lit pixel lies in the given row.
 * @tparam W Canvas width in pixels.
 * @param fx Effect whose canvas is sampled.
 * @param y Row to scan.
 * @return True if a non-black pixel is found in row y.
 */
template <int W> inline bool row_has_lit(const hs_test::StubEffect &fx, int y) {
  for (int x = 0; x < W; ++x)
    if (!is_black(fx.get_pixel(x, y)))
      return true;
  return false;
}

/**
 * @brief Placement oracle for a filled cap: the pole the shape is centred on is
 *        lit, the opposite pole stays dark, and the lit area matches the
 *        shape's angular extent.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param fx Effect whose canvas was drawn into.
 * @param cap_north True if the shape caps the +Y pole (row 0); false if it caps
 *        the -Y pole (row H-1). Star/PlanarPolygon centre on basis.v (north);
 *        Flower centres on its antipode -basis.v (south).
 * @param r_min Angular radius of the largest cap inscribed in the shape (the
 *        shape's smallest boundary radius, radians).
 * @param r_max Angular radius of the cap circumscribing the shape (radians).
 * @details Rows are uniform in polar angle, so a pole cap of angular radius r
 *          lights W*H*r/PI pixels and the shape's own count falls between its
 *          inscribed and circumscribed caps, with 3 rows of slack either way.
 */
template <int W, int H>
inline void expect_filled_cap(const hs_test::StubEffect &fx, bool cap_north,
                              float r_min, float r_max) {
  const int near_row = cap_north ? 0 : H - 1;
  const int far_row = cap_north ? H - 1 : 0;
  HS_EXPECT_TRUE((row_has_lit<W>(fx, near_row)));
  HS_EXPECT_FALSE((row_has_lit<W>(fx, far_row)));

  const float rows_per_radian = static_cast<float>(H) / math::PI_F;
  const float slack_rows = 3.0f;
  const int lit = static_cast<int>(count_lit_region<W, H>(fx));
  HS_EXPECT_GE(lit,
               static_cast<int>(W * (r_min * rows_per_radian - slack_rows)));
  HS_EXPECT_LE(lit,
               static_cast<int>(W * (r_max * rows_per_radian + slack_rows)));
}

/** @brief Verifies a filled Star caps the basis.v (+Y) pole and not the other. */
inline void test_star_pixel_placement() {
  constexpr int W = 96, H = 64;
  constexpr float RADIUS = 0.6f;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    math::Basis basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
    Scan::Star::draw<W, H, false>(pipe, c, basis, RADIUS, /*sides=*/5,
                                  [](const math::Vector &, Fragment &f) {
                                    f.color = Color4(Pixel(60000, 60000, 60000),
                                                     1.0f);
                                  });
  }
  fx.advance_display();
  // Tips at the circumradius, notches at Render::STAR_INNER_RATIO of it.
  const float r_max = RADIUS * (math::PI_F / 2.0f);
  expect_filled_cap<W, H>(fx, /*cap_north=*/true,
                          r_max * Render::STAR_INNER_RATIO, r_max);
  const math::Basis basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
  const float probe_radius = r_max * 0.7f;
  const auto probe = [&](float azimuth) {
    const math::Vector direction =
        basis.v * cosf(probe_radius) +
        (basis.u * cosf(azimuth) + basis.w * sinf(azimuth)) *
            sinf(probe_radius);
    const auto pixel = math::vector_to_pixel<W, H>(direction);
    return fx.get_pixel(static_cast<int>(std::round(pixel.x)) % W,
                        static_cast<int>(std::round(pixel.y)));
  };
  HS_EXPECT_FALSE(is_black(probe(0.0f)));
  HS_EXPECT_TRUE(is_black(probe(math::PI_F / 5.0f)));
}

/** @brief Verifies a filled PlanarPolygon caps the basis.v (+Y) pole, not the other. */
inline void test_planar_polygon_pixel_placement() {
  constexpr int W = 96, H = 64;
  constexpr float RADIUS = 0.6f;
  constexpr int SIDES = 6;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    math::Basis basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
    Scan::PlanarPolygon::draw<W, H, false>(
        pipe, c, basis, RADIUS, SIDES, [](const math::Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        });
  }
  fx.advance_display();
  // Vertices at the circumradius, edge midpoints at the apothem.
  const float r_max = RADIUS * (math::PI_F / 2.0f);
  expect_filled_cap<W, H>(fx, /*cap_north=*/true,
                          r_max * cosf(math::PI_F / SIDES), r_max);
}

/**
 * @brief Verifies a filled Flower caps the antipode (-Y) pole, not basis.v.
 * @details Flower's SDF scans from antipode = -basis.v.
 */
inline void test_flower_pixel_placement() {
  constexpr int W = 96, H = 64;
  constexpr float RADIUS = 0.6f;
  constexpr int SIDES = 6;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    math::Basis basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
    Scan::Flower::draw<W, H, false>(
        pipe, c, basis, RADIUS, SIDES, [](const math::Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        });
  }
  fx.advance_display();
  // Petal boundary (PI - r) * cos(local) == PI - r_max: the tip at local 0,
  // narrowest at the sector edge.
  const float r_max = RADIUS * (math::PI_F / 2.0f);
  const float r_min =
      math::PI_F - (math::PI_F - r_max) / cosf(math::PI_F / SIDES);
  expect_filled_cap<W, H>(fx, /*cap_north=*/false, r_min, r_max);
}

enum class SolidShape { PLANAR_POLYGON, SPHERICAL_POLYGON, FLOWER, STAR };

/**
 * @brief Verifies the typed solid-color path matches generic Scan output.
 */
inline void test_solid_color_path_matches_generic() {
  constexpr int W = 64, H = 48;
  constexpr float RADIUS = 0.72f;
  constexpr int SIDES = 5;
  constexpr float PHASE = 0.31f;
  const Color4 color(Pixel(51000, 23000, 9000), 0.37f);
  const SolidShape shapes[] = {SolidShape::PLANAR_POLYGON,
                               SolidShape::SPHERICAL_POLYGON,
                               SolidShape::FLOWER, SolidShape::STAR};

  for (SolidShape shape : shapes) {
    for (bool debug_bb : {false, true}) {
      std::vector<Pixel> generic_pixels;
      generic_pixels.reserve(W * H);

      {
        hs_test::StubEffect generic_fx(W, H);
        Pipeline<W, H> generic_pipeline;
        {
          Canvas canvas(generic_fx);
          math::Basis basis = math::make_basis(
              math::make_rotation(math::Vector(0.3f, 0.7f, 0.2f).normalized(),
                                  0.4f),
              math::Y_AXIS);
          auto shader = [&](const math::Vector &, Fragment &f) {
            f.color = color;
          };
          switch (shape) {
          case SolidShape::PLANAR_POLYGON:
            Scan::PlanarPolygon::draw<W, H, false>(generic_pipeline, canvas,
                                                   basis, RADIUS, SIDES, shader,
                                                   PHASE, debug_bb);
            break;
          case SolidShape::SPHERICAL_POLYGON:
            Scan::SphericalPolygon::draw<W, H, false>(generic_pipeline, canvas,
                                                      basis, RADIUS, SIDES,
                                                      shader, PHASE, debug_bb);
            break;
          case SolidShape::FLOWER:
            Scan::Flower::draw<W, H, false>(generic_pipeline, canvas, basis,
                                            RADIUS, SIDES, shader, PHASE,
                                            debug_bb);
            break;
          case SolidShape::STAR:
            Scan::Star::draw<W, H, false>(generic_pipeline, canvas, basis,
                                          RADIUS, SIDES, shader, PHASE,
                                          debug_bb);
            break;
          }
        }
        generic_fx.advance_display();
        // Guard against both paths drawing nothing.
        const size_t generic_lit = count_lit_region<W, H>(generic_fx);
        HS_EXPECT_GT(generic_lit, (size_t)0);
        capture_frame<W, H>(generic_fx, generic_pixels);
      }

      {
        hs_test::StubEffect solid_fx(W, H);
        Pipeline<W, H> solid_pipeline;
        {
          Canvas canvas(solid_fx);
          math::Basis basis = math::make_basis(
              math::make_rotation(math::Vector(0.3f, 0.7f, 0.2f).normalized(),
                                  0.4f),
              math::Y_AXIS);
          switch (shape) {
          case SolidShape::PLANAR_POLYGON:
            Scan::PlanarPolygon::draw_solid<W, H>(solid_pipeline, canvas, basis,
                                                  RADIUS, SIDES, color, PHASE,
                                                  debug_bb);
            break;
          case SolidShape::SPHERICAL_POLYGON:
            Scan::SphericalPolygon::draw_solid<W, H>(solid_pipeline, canvas,
                                                     basis, RADIUS, SIDES,
                                                     color, PHASE, debug_bb);
            break;
          case SolidShape::FLOWER:
            Scan::Flower::draw_solid<W, H>(solid_pipeline, canvas, basis,
                                           RADIUS, SIDES, color, PHASE,
                                           debug_bb);
            break;
          case SolidShape::STAR:
            Scan::Star::draw_solid<W, H>(solid_pipeline, canvas, basis, RADIUS,
                                         SIDES, color, PHASE, debug_bb);
            break;
          }
        }
        solid_fx.advance_display();
        for (int y = 0; y < H; ++y)
          for (int x = 0; x < W; ++x) {
            const Pixel &generic = generic_pixels[y * W + x];
            const Pixel &solid = solid_fx.get_pixel(x, y);
            HS_EXPECT_EQ(generic.r, solid.r);
            HS_EXPECT_EQ(generic.g, solid.g);
            HS_EXPECT_EQ(generic.b, solid.b);
          }
      }
    }
  }
}

/**
 * @brief Bounds spherical sine-distance framebuffer error at device resolution.
 * @details The two paths' distance gap is fast_acos' ~5e-5 rad wherever the
 *   circumscribed-disc clamp wins; the coverage ramp scales that by the
 *   quintic kernel's slope over 2*pixel_width, so a channel may swing up to
 *   ~2e-3 of full scale (about 100 codes at this color) near a vertex.
 */
inline void test_spherical_sine_distance_framebuffer_error() {
  constexpr int W = 288;
  constexpr int H = 144;
  const Color4 color(Pixel(61000, 43000, 17000), 0.73f);
  math::Basis basis = math::make_basis(
      math::make_rotation(math::Vector(-0.2f, 0.9f, 0.4f).normalized(), 0.63f),
      math::Y_AXIS);
  struct Case {
    float radius;
    int sides;
    float phase;
  };
  const Case cases[] = {{0.24f, 3, -2.2f},
                        {0.74f, 5, 0.37f},
                        {0.98f, 12, 4.8f},
                        {1.38f, 7, -0.61f}};

  auto render = [&]<bool SineDistance>(const Case &c) {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(fx);
      Scan::SphericalPolygon::draw_solid<W, H, SineDistance>(
          pipeline, canvas, basis, c.radius, c.sides, color, c.phase);
    }
    fx.advance_display();
    // Guard against both paths drawing nothing.
    const size_t lit = count_lit_region<W, H>(fx);
    HS_EXPECT_GT(lit, (size_t)0);
    std::vector<Pixel> pixels;
    capture_frame<W, H>(fx, pixels);
    return pixels;
  };

  size_t different_pixels = 0;
  int max_channel_error = 0;
  for (const Case &c : cases) {
    const std::vector<Pixel> exact = render.template operator()<false>(c);
    const std::vector<Pixel> sine = render.template operator()<true>(c);
    for (size_t i = 0; i < exact.size(); ++i) {
      int dr = std::abs(static_cast<int>(exact[i].r) - sine[i].r);
      int dg = std::abs(static_cast<int>(exact[i].g) - sine[i].g);
      int db = std::abs(static_cast<int>(exact[i].b) - sine[i].b);
      int pixel_error = std::max(dr, std::max(dg, db));
      if (pixel_error != 0)
        ++different_pixels;
      max_channel_error = std::max(max_channel_error, pixel_error);
    }
  }

  std::printf("spherical sine framebuffer samples=%d different=%zu max=%d\n",
              W * H * static_cast<int>(std::size(cases)), different_pixels,
              max_channel_error);
  HS_EXPECT_GT(different_pixels, static_cast<size_t>(0));
  HS_EXPECT_LE(different_pixels, static_cast<size_t>(512));
  HS_EXPECT_LE(max_channel_error, 128);
}

/**
 * @brief Verifies overlapping fills composite via the over operator at the
 *        shared pixel.
 * @details The 2D sink blends dst*(1-a) + src*a.
 */
inline void test_overlapping_fills_composite_blend() {
  constexpr int W = 96, H = 64;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe;

  constexpr uint16_t red = 60000, green = 50000;
  {
    Canvas c(fx);
    math::Basis basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
    // First fill: red at half alpha over black -> red*0.5.
    Scan::PlanarPolygon::draw<W, H, false>(
        pipe, c, basis, /*radius=*/0.6f, /*sides=*/6,
        [](const math::Vector &, Fragment &f) {
          f.color = Color4(Pixel(red, 0, 0), 0.5f);
        });
    // Second fill: green at half alpha over the red -> red*0.25 + green*0.5.
    Scan::PlanarPolygon::draw<W, H, false>(
        pipe, c, basis, /*radius=*/0.6f, /*sides=*/6,
        [](const math::Vector &, Fragment &f) {
          f.color = Color4(Pixel(0, green, 0), 0.5f);
        });
  }
  fx.advance_display();

  // Pole row: both fills fully covered, so plot alpha is frag.alpha.
  const Pixel &p = fx.get_pixel(W / 2, 0);
  HS_EXPECT_NEAR((int)p.r, (int)(red * 0.25f), 2);
  HS_EXPECT_NEAR((int)p.g, (int)(green * 0.5f), 2);
  HS_EXPECT_EQ((int)p.b, 0);
}
