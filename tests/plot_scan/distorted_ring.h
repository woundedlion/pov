/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_plot_scan.h.

// ============================================================================
// Plot::DistortedRing::sample  — angle-addition identity (LUT) vs direct
// ============================================================================

/**
 * @brief Verifies DistortedRing::sample with a zero shift function reduces to a
 *        plain ring built via the TrigLUT angle-addition identity, matching a
 *        direct cos/sin construction of the same ring per fragment.
 */
inline void test_distorted_ring_sample_angle_addition_identity() {
  constexpr int W = 64;
  constexpr int H = 64;

  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), W + 2);

  math::Basis b =
      math::make_basis(math::Quaternion(1, 0, 0, 0), math::Vector(0, 1, 0));
  const float radius = 0.5f;
  const float phase = 0.7f;
  ScalarFn zero_shift = [](float) { return 0.0f; };

  Plot::DistortedRing::sample<W, H>(points, b, radius, zero_shift, phase);

  HS_EXPECT_SIZE_OR_RETURN(points, (size_t)(W + 1));

  const std::vector<math::Vector> expected =
      ring_vertices_direct(b, radius, phase, W);
  for (int i = 0; i < W; ++i) {
    HS_EXPECT_NEAR(points[i].pos.x, expected[i].x, 2e-3f);
    HS_EXPECT_NEAR(points[i].pos.y, expected[i].y, 2e-3f);
    HS_EXPECT_NEAR(points[i].pos.z, expected[i].z, 2e-3f);
    HS_EXPECT_NEAR(points[i].pos.length(), 1.0f, 1e-3f);
  }
}

/**
 * @brief Verifies a non-zero shift function moves each sampled vertex to the
 *        colatitude it asks for, and that fn_point lands on the same ring at
 *        phase 0.
 */
inline void test_distorted_ring_shift_matches_fn_point() {
  constexpr int W = 64;
  constexpr int H = 64;

  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), W + 2);

  const math::Basis b =
      math::make_basis(math::Quaternion(0.92387953f, 0.0f, 0.38268343f, 0.0f),
                       math::Vector(0, 1, 0));
  const float radius = 0.5f;
  auto shift_shape = [](float t) {
    return 0.3f * sinf(2.0f * math::PI_F * t) + 0.1f;
  };
  ScalarFn shift_fn = shift_shape;

  Plot::DistortedRing::sample<W, H>(points, b, radius, shift_fn, 0.0f);
  HS_EXPECT_SIZE_OR_RETURN(points, (size_t)(W + 1));

  const Plot::RingFrame frame = Plot::ring_frame(b, radius);
  const float step = 2.0f * math::PI_F / W;
  for (int i = 0; i < W; ++i) {
    HS_CONTEXT("vertex", static_cast<long long>(i));
    const float angle = i * step;
    const float polar =
        frame.theta_eq + shift_shape(angle / (2.0f * math::PI_F));
    HS_EXPECT_NEAR(math::dot(points[i].pos, frame.basis.v), cosf(polar), 2e-3f);

    const math::Vector direct =
        Plot::DistortedRing::fn_point(shift_fn, b, radius, angle);
    HS_EXPECT_NEAR(direct.length(), 1.0f, 1e-3f);
    HS_EXPECT_NEAR(points[i].pos.x, direct.x, 2e-3f);
    HS_EXPECT_NEAR(points[i].pos.y, direct.y, 2e-3f);
    HS_EXPECT_NEAR(points[i].pos.z, direct.z, 2e-3f);
  }
}
