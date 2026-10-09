/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Local arena for sampling.
// ---------------------------------------------------------------------------

/** @brief Backing storage for the module-local sampling arena. */
inline uint8_t plot_scan_arena_buf[256 * 1024];

/**
 * @brief Returns the module-local arena used to bind Fragments for sampling.
 * @return Reference to a function-static Arena over plot_scan_arena_buf.
 * @details The arena is not reset between tests: every caller must wrap its
 *          allocations in a ScratchScope.
 */
inline Arena &plot_arena() {
  static Arena a(plot_scan_arena_buf, sizeof(plot_scan_arena_buf));
  return a;
}

/**
 * @brief Builds an orthonormal basis from a unit normal without a quaternion.
 * @param n Direction used as the basis normal; need not be pre-normalized.
 * @return Basis whose v is the normalized normal and whose u, w span the
 *         tangent plane.
 * @details Mirrors make_basis's construction.
 */
inline math::Basis basis_from_normal(const math::Vector &n) {
  math::Vector v = n.normalized();
  math::Vector ref =
      std::abs(math::dot(v, math::X_AXIS)) > math::COS_AXIS_PARALLEL
          ? math::Y_AXIS
          : math::X_AXIS;
  math::Vector u = math::cross(v, ref).normalized();
  math::Vector w = math::cross(v, u).normalized();
  return {u, v, w};
}

/** @brief Unit-sphere point on a basis-centered angular disk. */
inline math::Vector disk_point(const math::Basis &basis, float colat,
                               float az) {
  const math::Vector dir = basis.u * cosf(az) + basis.w * sinf(az);
  return (basis.v * cosf(colat) + dir * sinf(colat)).normalized();
}

/** @brief Random planar disk edge with the cull sweep's angular limits. */
inline void random_disk_edge(const math::Basis &basis, math::Vector &a,
                             math::Vector &b) {
  const float radius = hs::rand_f(0.2f, 1.4f);
  const float az = hs::rand_f(0, 2 * math::PI_F);
  a = disk_point(basis, radius, az);
  b = disk_point(basis, radius, az + hs::rand_f(0.3f, 2.3f));
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
