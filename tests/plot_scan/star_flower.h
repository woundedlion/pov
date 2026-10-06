/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_plot_scan.h.

// ============================================================================
// Plot::Star<Plot::PlanarProjection>::sample / Plot::Flower::sample
// ============================================================================

/**
 * @brief Verifies Star::sample emits 2*sides unit-length vertices plus one
 *        closing fragment matching the first vertex (closed loop, v0 == 1 at
 *        close), and that the vertices form a genuine star — not a flower.
 * @details Shape-discriminating check: the colatitude (angle from the shape's
 *          center axis) alternates between an outer and an inner radius, with
 *          inner/outer == STAR_INNER_RATIO.
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

// sample_positions() uses an angle-addition recurrence without normalization
// (core/render/plot/shapes.h: ~4e-6 off-unit drift); the bound is absolute.
constexpr float STAR_RECURRENCE_DRIFT =
    32.0f * std::numeric_limits<float>::epsilon();

/**
 * @brief Star radius-trig reuse reproduces per-vertex evaluation, and the two
 *        sampling entry points place the same vertices.
 * @details Positions against the per-vertex reference and sample() use float
 * tolerance. Cached and uncached sample_positions() share sample_positions_impl
 * and are compared bit for bit, as are the untouched zero registers.
 */
inline void test_star_sample_radius_trig_parity() {
  ScratchScope sc(plot_arena());
  const math::Basis b = math::make_basis(
      math::Quaternion(0.91f, 0.13f, -0.27f, 0.28f).normalized(), math::X_AXIS);

  for (int sides : {3, 7, 16}) {
    for (float radius : {0.01f, 0.57f, 1.0f, 1.73f, 1.99f}) {
      for (float phase : {-7.1f, 0.0f, 2.37f}) {
        ScratchScope iteration(plot_arena());
        Fragments actual;
        Fragments reference;
        Fragments positions;
        Fragments cached_positions;
        actual.bind(plot_arena(), sides * 2 + 2);
        reference.bind(plot_arena(), sides * 2 + 2);
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
        Plot::sample_closed_ring(reference, sides * 2, [&](int i) {
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
          return p;
        });

        HS_EXPECT_SIZE_OR_RETURN(actual, reference.size());
        HS_EXPECT_SIZE_OR_RETURN(actual, positions.size());
        HS_EXPECT_SIZE_OR_RETURN(positions, cached_positions.size());
        for (size_t i = 0; i < actual.size(); ++i) {
          HS_CONTEXT("star vertex", static_cast<long long>(i));
          // Separately-written expression tree: tolerance, not bit equality.
          HS_EXPECT_NEAR(actual[i].pos.x, reference[i].pos.x, 1e-6f);
          HS_EXPECT_NEAR(actual[i].pos.y, reference[i].pos.y, 1e-6f);
          HS_EXPECT_NEAR(actual[i].pos.z, reference[i].pos.z, 1e-6f);
          HS_EXPECT_NEAR(actual[i].v0, reference[i].v0, 1e-6f);
          HS_EXPECT_NEAR(actual[i].v1, reference[i].v1, 1e-6f);
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
      }
    }
  }
}

/**
 * @brief Continuous Star levels preserve the standard near-side geometry.
 * @details Compared componentwise: an acos-derived angle between two
 * near-parallel unit vectors bottoms out around 3.5e-4 rad.
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
 * @details Shape-discriminating check: every vertex sits at the SAME colatitude
 *          about the center axis (a constant polar radius).
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
