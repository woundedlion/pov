/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Plot::Star<Plot::PlanarProjection>::sample / Plot::Flower::sample
// ============================================================================

/**
 * @brief Verifies Star::sample emits 2*sides unit-length vertices plus one
 *        closing fragment matching the first vertex (closed loop, v0 == 1 at
 *        close), and that the vertices form a genuine star — not a flower.
 * @details Colatitude alternates outer/inner with inner/outer ==
 *          STAR_INNER_RATIO.
 */
inline void test_star_sample_unit_length_closed() {
  ScratchScope sc(plot_arena());
  Fragments points;
  const int sides = 5;
  points.bind(plot_arena(), sides * 2 + 2);

  math::Basis b = math::make_basis(
      math::Quaternion(0.91f, 0.13f, -0.27f, 0.28f).normalized(), math::X_AXIS);
  Plot::Star<Plot::PlanarProjection>::sample(points, b, 0.5f, sides, 0.0f);

  // 2*sides vertices + 1 close fragment.
  HS_EXPECT_SIZE_OR_RETURN(points, (size_t)(sides * 2 + 1));
  for (size_t i = 0; i < points.size(); ++i) {
    HS_EXPECT_NEAR(points[i].pos.length(), 1.0f, 1e-3f);
  }
  HS_EXPECT_NEAR(points.back().pos.x, points[0].pos.x, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.y, points[0].pos.y, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.z, points[0].pos.z, 1e-3f);
  HS_EXPECT_NEAR(points.back().v0, 1.0f, 1e-6f);

  // Alternating outer/inner colatitude about the center axis.
  const math::Vector axis = math::get_antipode(b, 0.5f).first.v;
  const float outer =
      math::angle_between(points[0].pos, axis); // even index -> outer
  const float inner =
      math::angle_between(points[1].pos, axis); // odd index  -> inner
  HS_EXPECT_GT(outer, inner + 1e-3f);
  HS_EXPECT_NEAR(inner / outer, Plot::STAR_INNER_RATIO, 1e-2f);
  for (int i = 0; i < sides * 2; ++i) {
    const float colat = math::angle_between(points[i].pos, axis);
    HS_EXPECT_NEAR(colat, (i % 2 == 0) ? outer : inner, 1e-3f);
  }
}

// Absolute drift of sample_positions()' unnormalized angle-addition recurrence.
constexpr float STAR_RECURRENCE_DRIFT =
    32.0f * std::numeric_limits<float>::epsilon();

/**
 * @brief Star radius-trig reuse reproduces per-vertex evaluation, the two
 *        sampling entry points place the same vertices, and the ring registers
 *        carry vertex fraction, cumulative arc length and vertex index.
 * @details Cached and uncached sample_positions() are compared bit for bit.
 */
inline void test_star_sample_radius_trig_parity() {
#if defined(HS_TEST_FAST_MATH)
  constexpr float ARC_STEP_TOL = 4e-5f;
  constexpr double ARC_TOTAL_REL_TOL = 6e-4;
#else
  constexpr float ARC_STEP_TOL = 1e-5f;
  constexpr double ARC_TOTAL_REL_TOL = 4e-6;
#endif
  ScratchScope sc(plot_arena());
  const math::Basis b = math::make_basis(
      math::Quaternion(0.91f, 0.13f, -0.27f, 0.28f).normalized(), math::X_AXIS);

  for (int sides : {3, 7, 16}) {
    for (float radius : {0.01f, 0.57f, 1.0f, 1.73f, 1.99f}) {
      for (float phase : {-7.1f, 0.0f, 2.37f}) {
        ScratchScope iteration(plot_arena());
        Fragments actual;
        Fragments positions;
        Fragments cached_positions;
        actual.bind(plot_arena(), sides * 2 + 2);
        positions.bind(plot_arena(), sides * 2 + 2);
        cached_positions.bind(plot_arena(), sides * 2 + 2);
        Plot::Star<Plot::PlanarProjection>::sample(actual, b, radius, sides,
                                                   phase);
        Plot::Star<Plot::PlanarProjection>::sample_positions(
            positions, b, radius, sides, phase);
        Plot::Star<Plot::PlanarProjection>::sample_positions(
            cached_positions, b, radius, sides, phase,
            Plot::Star<Plot::PlanarProjection>::radius_trig(radius),
            Plot::Star<Plot::PlanarProjection>::step_trig(sides));

        auto res = math::get_antipode(b, radius);
        const math::Basis &work_basis = res.first;
        const float outer_radius = res.second * (math::PI_F / 2.0f);
        const float inner_radius = outer_radius * Plot::STAR_INNER_RATIO;
        const float angle_step = math::PI_F / sides;
        const int verts = sides * 2;
        std::vector<math::Vector> reference(verts);
        for (int i = 0; i < verts; ++i) {
          const float theta = phase + i * angle_step;
          const float r = (i % 2 == 0) ? outer_radius : inner_radius;
          const float sin_r = sinf(r);
          const float cos_r = cosf(r);
          const float cos_t = cosf(theta);
          const float sin_t = sinf(theta);
          math::Vector p = (work_basis.v * cos_r) +
                           (work_basis.u * (cos_t * sin_r)) +
                           (work_basis.w * (sin_t * sin_r));
          p.normalize();
          reference[i] = p;
        }

        HS_EXPECT_SIZE_OR_RETURN(actual, static_cast<size_t>(verts + 1));
        HS_EXPECT_SIZE_OR_RETURN(actual, positions.size());
        HS_EXPECT_SIZE_OR_RETURN(positions, cached_positions.size());
        double arc_total = 0.0;
        for (size_t i = 0; i < actual.size(); ++i) {
          HS_CONTEXT("star vertex", static_cast<long long>(i));
          const math::Vector &ref = reference[i % verts];
          // Separately-written expression tree: tolerance, not bit equality.
          HS_EXPECT_NEAR(actual[i].pos.x, ref.x, 1e-6f);
          HS_EXPECT_NEAR(actual[i].pos.y, ref.y, 1e-6f);
          HS_EXPECT_NEAR(actual[i].pos.z, ref.z, 1e-6f);
          HS_EXPECT_NEAR(actual[i].v0, static_cast<float>(i) / verts, 1e-6f);
          HS_EXPECT_EQ(actual[i].v2, static_cast<float>(i));
          if (i > 0) {
            const float arc =
                math::angle_between(actual[i - 1].pos, actual[i].pos);
            arc_total += arc;
            HS_EXPECT_NEAR(actual[i].v1 - actual[i - 1].v1, arc, ARC_STEP_TOL);
          } else {
            HS_EXPECT_EQ(actual[i].v1, 0.0f);
          }
          // Recurrence positions omit per-vertex normalization.
          HS_EXPECT_NEAR(actual[i].pos.x, positions[i].pos.x,
                         STAR_RECURRENCE_DRIFT);
          HS_EXPECT_NEAR(actual[i].pos.y, positions[i].pos.y,
                         STAR_RECURRENCE_DRIFT);
          HS_EXPECT_NEAR(actual[i].pos.z, positions[i].pos.z,
                         STAR_RECURRENCE_DRIFT);
          HS_EXPECT_EQ(std::bit_cast<uint32_t>(positions[i].pos.x),
                       std::bit_cast<uint32_t>(cached_positions[i].pos.x));
          HS_EXPECT_EQ(std::bit_cast<uint32_t>(positions[i].pos.y),
                       std::bit_cast<uint32_t>(cached_positions[i].pos.y));
          HS_EXPECT_EQ(std::bit_cast<uint32_t>(positions[i].pos.z),
                       std::bit_cast<uint32_t>(cached_positions[i].pos.z));
          HS_EXPECT_EQ(std::bit_cast<uint32_t>(positions[i].v0), 0u);
          HS_EXPECT_EQ(std::bit_cast<uint32_t>(positions[i].v1), 0u);
          HS_EXPECT_EQ(std::bit_cast<uint32_t>(positions[i].v2), 0u);
        }
        HS_EXPECT_NEAR_REL(actual.back().v1, arc_total, ARC_TOTAL_REL_TOL);
      }
    }
  }
}

/**
 * @brief Continuous Star levels preserve the standard near-side geometry.
 * @details Compared componentwise: acos-derived angles between near-parallel
 * unit vectors bottom out far above this tolerance.
 */
inline void test_star_continuous_matches_standard_near_side() {
  constexpr float NEAR_SIDE_TOL = 1e-5f;
  ScratchScope sc(plot_arena());
  const math::Basis basis = math::make_basis(
      math::Quaternion(0.91f, 0.13f, -0.27f, 0.28f).normalized(), math::X_AXIS);
  for (float radius : {0.0f, 0.31f, 0.87f, 1.0f}) {
    ScratchScope iteration(plot_arena());
    Fragments standard;
    Fragments continuous;
    standard.bind(plot_arena(), 16);
    continuous.bind(plot_arena(), 16);
    Plot::Star<Plot::PlanarProjection>::sample(standard, basis, radius, 7,
                                               0.37f);
    Plot::Star<Plot::PlanarProjection>::sample_continuous(continuous, basis,
                                                          radius, 7, 0.37f);
    HS_EXPECT_SIZE_OR_RETURN(standard, continuous.size());
    for (size_t i = 0; i < standard.size(); ++i) {
      HS_EXPECT_NEAR(standard[i].pos.x, continuous[i].pos.x, NEAR_SIDE_TOL);
      HS_EXPECT_NEAR(standard[i].pos.y, continuous[i].pos.y, NEAR_SIDE_TOL);
      HS_EXPECT_NEAR(standard[i].pos.z, continuous[i].pos.z, NEAR_SIDE_TOL);
    }
  }
}

/** @brief Continuous Star vertices do not jump at the equatorial level. */
inline void test_star_continuous_crosses_equator() {
  ScratchScope sc(plot_arena());
  const math::Basis basis = math::make_basis(math::Quaternion(), math::X_AXIS);
  Fragments below;
  Fragments seam;
  Fragments above;
  below.bind(plot_arena(), 16);
  seam.bind(plot_arena(), 16);
  above.bind(plot_arena(), 16);
  Plot::Star<Plot::PlanarProjection>::sample_continuous_positions(
      below, basis, 0.9999f, 7, 0.37f);
  Plot::Star<Plot::PlanarProjection>::sample_continuous_positions(
      seam, basis, 1.0f, 7, 0.37f);
  Plot::Star<Plot::PlanarProjection>::sample_continuous_positions(
      above, basis, 1.0001f, 7, 0.37f);
  HS_EXPECT_SIZE_OR_RETURN(below, seam.size());
  HS_EXPECT_SIZE_OR_RETURN(above, seam.size());
  for (size_t i = 0; i < seam.size(); ++i) {
    HS_EXPECT_LT(math::angle_between(below[i].pos, seam[i].pos), 0.001f);
    HS_EXPECT_LT(math::angle_between(seam[i].pos, above[i].pos), 0.001f);
  }
}

/** @brief The final continuous Star level collapses to the opposite pole. */
inline void test_star_continuous_collapses_at_antipode() {
  ScratchScope sc(plot_arena());
  const math::Basis basis = math::make_basis(
      math::Quaternion(0.91f, 0.13f, -0.27f, 0.28f).normalized(), math::X_AXIS);
  Fragments points;
  points.bind(plot_arena(), 16);
  const int sides = 7;
  Plot::Star<Plot::PlanarProjection>::sample_continuous_positions(
      points, basis, 2.0f, sides, 0.37f);
  HS_EXPECT_SIZE_OR_RETURN(points, (size_t)(sides * 2 + 1));
  for (const Fragment &point : points)
    HS_EXPECT_LT(math::angle_between(point.pos, -basis.v), 0.001f);
}

/**
 * @brief Verifies Flower::sample emits 2*sides unit-length vertices plus one
 *        closing fragment matching the first vertex (closed loop), and that the
 *        vertices form a genuine flower — constant radius, not a star.
 */
inline void test_flower_sample_unit_length_closed() {
  ScratchScope sc(plot_arena());
  Fragments points;
  const int sides = 6;
  points.bind(plot_arena(), sides * 2 + 2);

  math::Basis b = math::make_basis(
      math::Quaternion(0.91f, 0.13f, -0.27f, 0.28f).normalized(), math::X_AXIS);
  Plot::Flower::sample(points, b, 0.5f, sides, 0.0f);

  HS_EXPECT_SIZE_OR_RETURN(points, (size_t)(sides * 2 + 1));
  for (size_t i = 0; i < points.size(); ++i) {
    HS_EXPECT_NEAR(points[i].pos.length(), 1.0f, 1e-3f);
  }
  HS_EXPECT_NEAR(points.back().pos.x, points[0].pos.x, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.y, points[0].pos.y, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.z, points[0].pos.z, 1e-3f);

  // Constant colatitude about the center axis.
  const math::Vector axis = math::get_antipode(b, 0.5f).first.v;
  const float colat0 = math::angle_between(points[0].pos, axis);
  for (int i = 0; i < sides * 2; ++i) {
    const float colat = math::angle_between(points[i].pos, axis);
    HS_EXPECT_NEAR(colat, colat0, 1e-3f);
  }
}
