/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Math spatial death fixtures and guard cases.

/**
 * @brief Death case: normalizing a degenerate (zero-length) vector must trap.
 * @details Math-core surface — length below epsilon fires the normalize guard.
 */
inline void case_normalize_zero() {
  math::Vector z{opaque(0.0f), opaque(0.0f), opaque(0.0f)};
  math::Vector n = z.normalized(); // length < eps -> HS_CHECK
  if (n.x == 42.0f)
    std::printf("x");
}

/** @brief In-place normalization rejects a zero-length vector. */
inline void case_vector_normalize_in_place_zero() {
  math::Vector v{opaque(0.0f), opaque(0.0f), opaque(0.0f)};
  v.normalize();
  if (v.x == 42.0f)
    std::printf("x");
}

/** @brief Death case: rotating in a degenerate coordinate plane must trap. */
inline void case_rotate_plane_degenerate() {
  math::Mat4 m = math::Mat4::identity();
  math::rotate_plane(m, opaque(1), opaque(1), 0.5f); // a == b -> HS_CHECK
  if (m.m[0][0] == 42.0f)
    std::printf("x");
}

/** @brief Death case: measuring an angle from a degenerate vector must trap. */
inline void case_angle_between_zero() {
  math::Vector zero{opaque(0.0f), opaque(0.0f), opaque(0.0f)};
  math::Vector x{opaque(1.0f), opaque(0.0f), opaque(0.0f)};
  float angle = math::angle_between(zero, x);
  if (angle == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: normalizing a NaN vector must trap.
 * @details Math-core surface — a NaN coordinate poisons the length to NaN, and
 *          `NaN >= epsilon` is false, so the normalize guard fires. The suite's
 *          NaN/Inf fault case: proves a non-finite producer is trapped at the
 *          math seam rather than silently propagating NaN into geometry.
 */
inline void case_normalize_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  math::Vector bad{nan, opaque(0.0f), opaque(0.0f)};
  math::Vector n = bad.normalized(); // length is NaN -> HS_CHECK fails
  if (n.x == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: requesting more KDTree neighbors than MAX_K must trap.
 * @details Spatial surface — k beyond the MAX_K-sized result/heap buffers makes
 *          nearest() trap rather than silently capping the result and masking
 *          the caller's sizing mistake.
 */
inline void case_spatial_knn_over_max() {
  static uint8_t buf[512];
  Arena a(buf, sizeof(buf));
  math::Vector pts[2] = {math::Vector(1.0f, 0.0f, 0.0f),
                         math::Vector(0.0f, 1.0f, 0.0f)};
  KDTree tree(a, std::span<const math::Vector>(pts, 2));
  // Tree is non-empty and k > 0, so the k <= MAX_K guard is reached.
  auto r =
      tree.nearest(math::Vector(1.0f, 0.0f, 0.0f),
                   opaque<size_t>(KDTree::MAX_K + 1)); // k > MAX_K -> HS_CHECK
  if (r.size() == static_cast<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: a lattice index outside [0, RD_N) must trap.
 * @details Spatial surface — node()'s index maps affinely onto the sphere, so an
 *          out-of-range one would silently return a direction off the lattice
 *          instead of naming the caller's mistake.
 */
inline void case_reaction_graph_node_index_out_of_range() {
  math::Vector v = ReactionGraph::node(opaque(ReactionGraph::RD_N));
  if (v.x == 0x7fff)
    std::printf("x");
}

/**
 * @brief Death case: a neighbor-table slot outside the lattice must trap.
 * @details Spatial surface — CubemapLUT's hill-climb and the reaction-diffusion
 *          Laplacian subscript neighbors[] rows unguarded, so validate_neighbors()
 *          traps on a slot that is not a node index before the first such read.
 */
inline void case_reaction_graph_slot_out_of_range() {
  static int16_t table[ReactionGraph::RD_N][ReactionGraph::RD_K] = {};
  table[0][0] = opaque<int16_t>(-1);
  ReactionGraph::validate_neighbors(table);
}

/**
 * @brief Death case: a Flywheel period of zero must trap at construction.
 * @details POV-sync surface — position() divides the int32 elapsed window by the
 *          period, so a zero divides by zero and an over-large one voids the
 *          signed-safe coast window; the constructor rejects both before the
 *          driver ever schedules a column.
 */
inline void case_flywheel_period_zero() {
  pov::sync::Config cfg;
  cfg.cycles_per_half_rev = opaque<uint32_t>(0);
  pov::sync::Flywheel fw(cfg); // period 0 -> HS_CHECK
  (void)fw;
}

/**
 * @brief Death case: a virtual height of one row must trap in the phi mapping.
 * @details Geometry surface — the row-to-angle scale divides by (h_virt - 1),
 *          so a single-row canvas would map every row to a non-finite phi.
 */
inline void case_y_to_phi_degenerate_height() {
  float phi =
      math::y_to_phi_virtual(opaque(0.0f), opaque(1)); // divisor 0 -> trap
  if (phi == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: make_basis with a non-unit quaternion must trap.
 * @details Geometry surface — the rotation assumes a unit quaternion, so a
 *          finite but over-long one would scale and shear the frame rather than
 *          rotate it; the guard fires before the axes are built.
 */
inline void case_make_basis_nonunit_quaternion() {
  math::Quaternion q(opaque(2.0f), opaque(0.0f), opaque(0.0f), opaque(0.0f));
  math::Vector normal{opaque(0.0f), opaque(1.0f), opaque(0.0f)};
  math::Basis b = math::make_basis(q, normal); // |q| = 2 -> HS_CHECK
  if (b.u.x == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: parallel transport between antipodal endpoints must trap.
 * @details Geometry surface — the great circle through antipodes is
 *          ill-determined and the transport divides by 1 + dot, so the guard
 *          fires before the tangent is amplified.
 */
inline void case_parallel_transport_antipodal() {
  math::Vector from{opaque(1.0f), opaque(0.0f), opaque(0.0f)};
  math::Vector to{opaque(-1.0f), opaque(0.0f), opaque(0.0f)};
  math::Vector tangent{opaque(0.0f), opaque(1.0f), opaque(0.0f)};
  math::Vector t =
      math::parallel_transport(from, to, tangent); // dot = -1 -> HS_CHECK
  if (t.x == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a polyhedral fold that never converges must trap.
 * @details Lens surface — two opposed mirrors are not a chamber: each pass
 *          reflects the direction back across the other, so the bounded
 *          reflection loop exhausts its passes and fires the guard.
 */
inline void case_polyhedral_kaleidoscope_no_converge() {
  const std::array<math::Vector, 3> mirrors = {
      math::Vector(opaque(1.0f), 0.0f, 0.0f),
      math::Vector(opaque(-1.0f), 0.0f, 0.0f),
      math::Vector(0.0f, opaque(1.0f), 0.0f)};
  math::Vector v = lenses::polyhedral_kaleidoscope_lens(
      math::Vector(opaque(0.5f), opaque(0.5f), 0.0f), mirrors);
  if (v.x == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a negative equator sample count must trap.
 * @details Spherical-field surface — the count sizes every ring's longitude
 *          walk, so a negative one underflows the per-ring sample allocation.
 */
inline void case_spherical_field_negative_equator_samples() {
  hs::SphericalFieldLayout<32, 16, 0> layout(4, 0, 0, opaque(-1));
  if (layout.sample_count() == 42)
    std::printf("x");
}

/**
 * @brief Death case: a NaN endpoint fed to slerp must trap.
 * @details Math-core surface — the NaN poisons interpolation through both
 *          branches into the final strict normalized(), which traps rather than
 *          emitting a NaN direction into geometry. Proves the non-finite input
 *          is caught at the slerp seam, not just at bare normalize().
 */
inline void case_slerp_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  math::Vector bad{nan, opaque(0.0f), opaque(0.0f)};
  math::Vector dst{opaque(0.0f), opaque(0.0f), opaque(1.0f)};
  math::Vector v =
      math::slerp(bad, dst, opaque(0.5f)); // NaN -> normalized() -> HS_CHECK
  if (v.x == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: make_rotation(from, to) with a NaN source must trap.
 * @details A NaN component fails the unit-vector precondition before rotation
 *          arithmetic, complementing the finite non-unit input case.
 */
inline void case_make_rotation_vectors_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  math::Vector from{nan, opaque(0.0f), opaque(0.0f)};
  math::Vector to{opaque(0.0f), opaque(0.0f), opaque(1.0f)};
  math::Quaternion q = math::make_rotation(from, to);
  if (q.r == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: make_rotation(axis, theta) with a NaN angle must trap.
 * @details Math-core surface — cos/sin of a NaN poison the quaternion, and its
 *          normalized() traps on the NaN magnitude.
 */
inline void case_make_rotation_angle_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  math::Vector axis{opaque(0.0f), opaque(1.0f), opaque(0.0f)};
  math::Quaternion q =
      math::make_rotation(axis, nan); // NaN quat -> normalized() -> HS_CHECK
  if (q.r == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: make_basis with a NaN normal must trap.
 * @details Geometry surface — rotate(normal,.).normalized() is the first strict
 *          normalize in the basis construction and traps on the NaN-poisoned
 *          vector rather than returning a garbage frame.
 */
inline void case_make_basis_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  math::Vector normal{nan, opaque(0.0f), opaque(0.0f)};
  math::Basis b = math::make_basis(math::Quaternion(),
                                   normal); // NaN -> normalized() -> HS_CHECK
  if (b.u.x == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: an active noise_transform fed a non-finite direction must
 *        trap, not propagate NaN/Inf into the rendered geometry.
 * @details The structural audit rejects non-finite directions before noise sampling.
 */
inline void case_noise_transform_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  Animation::NoiseParams p;
  p.amplitude = opaque(0.5f); // active path (skips the zero-amplitude no-op)
  p.scale = opaque(4.0f);
  p.time = opaque(1.0f);
  math::Vector v{nan, opaque(0.0f), opaque(0.0f)};
  math::Vector r = noise_transform(v, p);
  if (r.x == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: make_rotation(from, to) with a non-unit source must trap.
 * @details Math-core surface — the d-based parallel/antiparallel branches assume
 *          |from| = |to| = 1, so a finite but non-unit input must trap at the
 *          unit-vector guard rather than silently skewing the rotation angle.
 */
inline void case_make_rotation_nonunit() {
  math::Vector from{opaque(2.0f), opaque(0.0f), opaque(0.0f)}; // |from| = 2
  math::Vector to{opaque(0.0f), opaque(1.0f), opaque(0.0f)};
  math::Quaternion q = math::make_rotation(from, to); // |from| != 1 -> HS_CHECK
  if (q.r == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: a ring index past the last ring must trap.
 * @details SphericalFieldLayout surface — the chain walk saturates at row H-1
 *          while the offset keeps accumulating, so an out-of-range index would
 *          otherwise hand back a Ring pointing past the sample array.
 */
inline void case_spherical_field_ring_index_oob() {
  hs::SphericalFieldLayout<32, 16, 0> layout(4);
  const auto ring = layout.ring(opaque(layout.ring_count()));
  if (ring.offset == 42)
    std::printf("x");
}

/**
 * @brief Death case: populating past the last ring must trap.
 * @details next_ring() saturates at the last ring, so an overrunning band would
 *          re-populate it and leave the caller believing it wrote fresh rings.
 */
inline void case_spherical_field_populate_ring_end_oob() {
  constexpr hs::SphericalFieldLayout<32, 16, 0> layout(4);
  static float values[layout.sample_count()];
  hs::SphericalField<float, 32, 16, 0> field(values, layout);
  field.populate(0, opaque(layout.ring_count()),
                 [](const math::Vector &v, const auto &) { return v.y; });
  if (values[0] == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: an order past the degree must trap.
 * @details reduced_legendre() has no term to recur on for |m| > l and returns
 *          0, so every sample of the mode comes back black.
 */
inline void case_spherical_harmonic_order_over_degree() {
  const float n = SHMath::normalization(opaque(2), 3);
  if (n == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: a negative flat harmonic index must trap.
 * @details sqrtf of a negative argument is NaN, and the cast of a NaN to int
 *          is undefined, so the decoded level would be arbitrary.
 */
inline void case_spherical_harmonic_decode_negative_index() {
  auto [l, m] = SHMath::decode_lm(opaque(-1));
  if (l == 42 && m == 42)
    std::printf("x");
}

/**
 * @brief Death case: an infill band past the rendered domain must trap.
 * @details A south_infill wider than H puts every row at full longitude
 *          resolution, multiplying sample_count() by the spacing; the arena
 *          would then overflow at an unrelated call site.
 */
inline void case_spherical_field_infill_over_domain() {
  hs::SphericalFieldLayout<32, 16, 0> layout(4, 0, opaque(17));
  if (layout.sample_count() == 42)
    std::printf("x");
}

inline void case_latitude_geometry_degenerate_height() {
  math::LatitudeGeometry geometry(opaque(1), 0.1f, 3.0f);
  if (geometry.row_to_phi(0) == 42.0f)
    std::printf("x");
}

inline void case_latitude_geometry_reversed_span() {
  math::LatitudeGeometry geometry(16, opaque(2.0f), 1.0f);
  if (geometry.row_to_phi(0) == 42.0f)
    std::printf("x");
}

inline void case_peirce_invalid_layout() {
  (void)projections::peirce_projection(
      math::Vector(0, 1, 0), 0,
      static_cast<projections::PeirceLayout>(opaque<uint8_t>(255)), 0);
}

/** @brief Rejects a square-wave duty cycle outside its unit interval. */
inline void case_square_wave_invalid_duty() {
  (void)math::square_wave(0.f, 1.f, 1.f, opaque(-.1f), 0.f);
}

/** @brief Rejects a Fibonacci spiral with no points. */
inline void case_fib_spiral_zero_points() {
  (void)math::fib_spiral(opaque(0), .5f, 0);
}
