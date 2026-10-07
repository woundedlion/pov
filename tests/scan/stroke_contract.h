/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Scan shader v2 contract, CSG stretch and stroke AA
// ============================================================================

/**
 * @brief Capturing plot sink that records the AA alpha process_pixel forwards.
 * @details Bypasses the canvas-blend round trip. The recorded alpha is
 * frag.alpha (=1) times the AA alpha; count tracks whether the pixel was drawn
 * at all.
 */
struct AlphaSink {
  float last_alpha =
      -1.0f;     /**< AA alpha from the most recent plot, -1 if none. */
  int count = 0; /**< Number of times plot() was invoked. */
  /**
   * @brief Records the forwarded AA alpha and increments the plot count.
   * @param a Anti-aliasing alpha forwarded by process_pixel.
   */
  void plot(Canvas &, int, int, const Pixel &, float, float a) {
    last_alpha = a;
    ++count;
  }
};

/**
 * @brief Runs process_pixel for a shape at a surface point and returns its AA alpha.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param shape SDF shape to rasterize at the sample point.
 * @param p Surface point on the unit sphere to evaluate.
 * @param c Canvas providing clip and projection context.
 * @param count Optional out-param; receives how many times the pixel was plotted (0 or 1).
 * @return The forwarded AA alpha, or -1 if the pixel was not drawn.
 */
template <int W, int H>
inline float scan_alpha_at(const auto &shape, const math::Vector &p, Canvas &c,
                           int *count = nullptr) {
  AlphaSink sink;
  SDF::DistanceResult res;
  Fragment frag;
  Scan::process_pixel<W, H, false>(
      0, 0, p, sink, c, shape,
      [](const math::Vector &, Fragment &f) {
        f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
      },
      /*debug_bb=*/false, res, frag);
  if (count)
    *count = sink.count;
  return sink.last_alpha;
}

/**
 * @brief Verifies Scan's v2 contract for stroke and solid shaders.
 * @details A stroke shader that writes v2 as alpha produces coverage squared
 * at the sink; solid shaders receive zero.
 */
inline void test_scan_shader_v2_contract() {
  constexpr int W = 288, H = 144;
  hs_test::StubEffect fx(W, H);
  Canvas canvas(fx);
  math::Basis basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
  constexpr float RADIUS = 0.5f;
  constexpr float THICKNESS = 0.1f;
  const float phi = RADIUS * (math::PI_F / 2.0f) + THICKNESS * 0.5f;
  const math::Vector p(sinf(phi), cosf(phi), 0.0f);
  SDF::Ring ring(basis, RADIUS, THICKNESS);

  SDF::DistanceResult expected_result;
  ring.distance<false>(p, expected_result);
  const float expected_coverage =
      math::quintic_kernel(-expected_result.dist / expected_result.size);
  const float RAW_DIST_FACTOR = math::quintic_kernel(
      1.0f -
      hs::clamp(expected_result.raw_dist / expected_result.size, 0.0f, 1.0f));

  AlphaSink sink;
  SDF::DistanceResult result;
  Fragment frag;
  float shader_coverage = -1.0f;
  Scan::process_pixel<W, H, false>(
      0, 0, p, sink, canvas, ring,
      [&](const math::Vector &, Fragment &f) {
        shader_coverage = f.v2;
        f.color = Color4(Pixel(60000, 60000, 60000), f.v2);
      },
      false, result, frag);

  HS_EXPECT_EQ(sink.count, 1);
  HS_EXPECT_GT(expected_coverage, 0.1f);
  HS_EXPECT_LT(expected_coverage, 0.9f);
  HS_EXPECT_NEAR(shader_coverage, expected_coverage, 1e-6f);
  HS_EXPECT_NEAR(sink.last_alpha, expected_coverage * expected_coverage, 1e-6f);
  HS_EXPECT_NEAR(sink.last_alpha, expected_coverage * RAW_DIST_FACTOR, 2e-6f);

  SDF::PlanarPolygon solid(basis, 0.5f / (math::PI_F / 2.0f), 6, 0.0f);
  AlphaSink solid_sink;
  float solid_v2 = -1.0f;
  Scan::process_pixel<W, H, false>(
      0, 0, basis.v, solid_sink, canvas, solid,
      [&](const math::Vector &, Fragment &f) {
        solid_v2 = f.v2;
        f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
      },
      false, result, frag);
  HS_EXPECT_EQ(solid_sink.count, 1);
  HS_EXPECT_EQ(solid_v2, 0.0f);
}

/** @brief Metadata fixture for CSG distance-report stretch forwarding. */
struct CsgStretchFixture {
  static constexpr bool is_solid = true;
  static constexpr bool BLENDS_SMOOTHLY = true;
  float stretch;

  friend float report_stretch(const CsgStretchFixture &shape) {
    return shape.stretch;
  }
};

/**
 * @brief A composite's report_stretch bounds every child, not only the first.
 * @details SDF::Face reports gnomonic-plane distance, whose stretch over an
 * angular step is 1 + max_dist_sq. A Union of two differently sized faces
 * reports the larger factor either way round.
 */
inline void test_report_stretch_forwards_through_csg() {
  constexpr int H = 144, HV = H + hs::H_OFFSET;
  const math::Basis basis =
      math::make_basis(math::Quaternion(), math::Vector(0, 1, 0));
  math::Vector verts[6];
  uint16_t small_idx[3], large_idx[3];
  for (int i = 0; i < 3; ++i) {
    const float a = (2.0f * math::PI_F * static_cast<float>(i)) / 3.0f;
    const auto ring = [&](float polar) {
      return (basis.v * cosf(polar) +
              (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(polar))
          .normalized();
    };
    verts[i] = ring(0.15f);
    verts[i + 3] = ring(0.90f);
    small_idx[i] = static_cast<uint16_t>(i);
    large_idx[i] = static_cast<uint16_t>(i + 3);
  }
  static SDF::FaceScratchBuffer small_scratch, large_scratch;
  const SDF::Face small(std::span<const math::Vector>(verts, 6),
                        std::span<const uint16_t>(small_idx, 3), small_scratch,
                        HV, H);
  const SDF::Face large(std::span<const math::Vector>(verts, 6),
                        std::span<const uint16_t>(large_idx, 3), large_scratch,
                        HV, H);

  const float small_stretch = Scan::report_stretch(small);
  const float large_stretch = Scan::report_stretch(large);
  HS_EXPECT_GT(large_stretch, small_stretch);

  HS_EXPECT_EQ(
      Scan::report_stretch(SDF::Union<SDF::Face, SDF::Face>{small, large}),
      large_stretch);
  HS_EXPECT_EQ(
      Scan::report_stretch(SDF::Union<SDF::Face, SDF::Face>{large, small}),
      large_stretch);
  HS_EXPECT_EQ(Scan::report_stretch(
                   SDF::Intersection<SDF::Face, SDF::Face>{small, large}),
               large_stretch);
  HS_EXPECT_EQ(
      Scan::report_stretch(SDF::Subtract<SDF::Face, SDF::Face>{small, large}),
      large_stretch);
  HS_EXPECT_EQ(Scan::report_stretch(
                   SDF::Intersection<SDF::Face, SDF::Face>{large, small}),
               large_stretch);
  HS_EXPECT_EQ(
      Scan::report_stretch(SDF::Subtract<SDF::Face, SDF::Face>{large, small}),
      large_stretch);
  const CsgStretchFixture smooth_small{small_stretch},
      smooth_large{large_stretch};
  HS_EXPECT_EQ(Scan::report_stretch(
                   SDF::SmoothUnion<CsgStretchFixture, CsgStretchFixture>{
                       smooth_small, smooth_large, 0.1f}),
               large_stretch);
  HS_EXPECT_EQ(Scan::report_stretch(
                   SDF::SmoothUnion<CsgStretchFixture, CsgStretchFixture>{
                       smooth_large, smooth_small, 0.1f}),
               large_stretch);
  HS_EXPECT_EQ(
      Scan::report_stretch(SDF::AngularRepeat<SDF::Face>(large, 3, basis.v)),
      large_stretch);
}

/**
 * @brief Verifies a CSG composite renders each stroke with its own thickness.
 * @details A Union<thin, thick> evaluated at a point inside the thin stroke
 * yields the same AA alpha as the bare thin Line and a different alpha from the
 * bare thick Line.
 */
inline void test_csg_stroke_aa_uses_winning_child_thickness() {
  constexpr int W = 288, H = 144;
  hs_test::StubEffect fx(W, H);
  Canvas c(fx);

  const float thin = 0.05f, thick = 0.30f;
  // Thin line: equatorial arc +X -> +Z (great circle in the y=0 plane).
  SDF::Line thin_line(math::Vector(1, 0, 0), math::Vector(0, 0, 1), thin);
  // The distant thick sibling makes sibling-max AA width observable while
  // the Union selects the thin line.
  SDF::Line thick_line(math::Vector(-1, 0, 0), math::Vector(0, 0, -1), thick);
  SDF::Union<SDF::Line, SDF::Line> u(thin_line, thick_line);
  // Same geometry as `thin_line` but standalone, for the contrast check.
  SDF::Line thick_solo(math::Vector(1, 0, 0), math::Vector(0, 0, 1), thick);

  static_assert(!SDF::Union<SDF::Line, SDF::Line>::is_solid,
                "a Union of strokes is not solid");

  // Point at geodesic distance 0.025 (half the thin thickness) north of the
  // thin arc, projecting to azimuth 45 deg (well inside the arc).
  const float dist = 0.025f;
  math::Vector p(cosf(dist) * cosf(math::PI_F / 4), sinf(dist),
                 cosf(dist) * sinf(math::PI_F / 4));

  int n_union = 0, n_bare = 0;
  float a_union = scan_alpha_at<W, H>(u, p, c, &n_union);
  float a_thin = scan_alpha_at<W, H>(thin_line, p, c, &n_bare);
  float a_thick = scan_alpha_at<W, H>(thick_solo, p, c);

  HS_EXPECT_EQ(n_union, 1);
  HS_EXPECT_EQ(n_bare, 1);
  // The composite reproduces the winning child's own AA.
  HS_EXPECT_NEAR(a_union, a_thin, 1e-4f);
  HS_EXPECT_GT(fabsf(a_thin - a_thick), 0.1f);
}

/**
 * @brief Verifies the rasterized ring lights the analytically predicted row.
 * @details A ring of normalized radius r is centered on the basis axis at polar
 * angle target = r*(PI/2), lighting a single latitude band whose center row is
 * phi_to_y(target). Rows well away from the band stay dark.
 */
inline void test_ring_rasterize_lights_expected_row() {
  constexpr int W = 96, H = 48;

  auto centroid_and_band = [](float radius) {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipe;
    {
      Canvas c(fx);
      math::Basis basis = math::make_basis(
          math::Quaternion(), math::Y_AXIS); // axis = north pole (+Y)
      Scan::Ring::draw<W, H, false>(pipe, c, basis, radius, /*thickness=*/0.05f,
                                    [](const math::Vector &, Fragment &f) {
                                      f.color = Color4(
                                          Pixel(60000, 60000, 60000), 1.0f);
                                    });
    }
    fx.advance_display();

    int lit[H] = {0};
    long total = 0, weighted = 0;
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        if (!is_black(fx.get_pixel(x, y))) {
          lit[y]++;
          total++;
          weighted += y;
        }
    /** @brief Per-radius result: lit-pixel centroid row, total lit count, and per-row lit counts. */
    struct R {
      float centroid;
      int total;
      int lit[H];
    };
    R r;
    r.centroid = total ? static_cast<float>(weighted) / total : -1.0f;
    r.total = static_cast<int>(total);
    for (int y = 0; y < H; ++y)
      r.lit[y] = lit[y];
    return r;
  };

  for (float radius : {0.5f, 1.0f}) {
    float target = radius * (math::PI_F / 2.0f);
    float expected_y = math::phi_to_y<H>(target);
    auto r = centroid_and_band(radius);

    // The ring is a full circle of latitude: most columns light up its row.
    HS_EXPECT_GT(r.total, W / 2);
    // Lit-pixel centroid lands on the analytically predicted row.
    HS_EXPECT_NEAR(r.centroid, expected_y, 1.0f);
    // Rows far from the band (> 4 px away) are dark.
    int ey = static_cast<int>(expected_y + 0.5f);
    for (int y = 0; y < H; ++y)
      if (std::abs(y - ey) > 4)
        HS_EXPECT_EQ(r.lit[y], 0);
  }
}

/**
 * @brief Verifies the stroke anti-aliasing alpha is a monotone ramp.
 * @details The alpha falls from ~1 at the ring centerline to 0 at its outer
 * surface through intermediate values, sampled radially outward across the
 * band.
 */
inline void test_stroke_aa_is_monotone_ramp() {
  constexpr int W = 288, H = 144;
  hs_test::StubEffect fx(W, H);
  Canvas c(fx);

  const float radius = 0.5f;     // centerline at polar PI/4
  const float thickness = 0.10f; // band half-width in radians
  const float target = radius * (math::PI_F / 2.0f);
  math::Basis basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
  SDF::Ring ring(basis, radius, thickness);

  // March outward from the centerline along the az=0 meridian.
  const int N = 12;
  float prev = 1.0f;
  bool saw_one = false, saw_zero = false, saw_mid = false;
  for (int i = 0; i < N; ++i) {
    float delta = (thickness * 1.4f) * i / (N - 1); // 0 .. 1.4*thickness
    float ph = target + delta;
    math::Vector p(sinf(ph), cosf(ph), 0.0f); // az=0, polar angle ph
    int count = 0;
    float a = scan_alpha_at<W, H>(ring, p, c, &count);
    if (count == 0)
      a = 0.0f; // outside the stroke -> not drawn -> alpha 0

    HS_EXPECT_LE(a, prev + 1e-4f); // monotone non-increasing outward
    prev = a;

    if (a > 0.95f)
      saw_one = true;
    if (a < 0.05f)
      saw_zero = true;
    if (a > 0.1f && a < 0.9f)
      saw_mid = true;
  }

  // Centerline ~opaque, far edge ~transparent, with a genuine ramp between.
  HS_EXPECT_TRUE(saw_one);
  HS_EXPECT_TRUE(saw_zero);
  HS_EXPECT_TRUE(saw_mid);
}
