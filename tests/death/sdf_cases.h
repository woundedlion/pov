/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Sdf death fixtures and guard cases.

/**
 * @brief Death case: a polygon with fewer than three sides must trap.
 * @details SDF surface — the sector fold divides a full turn by the side count,
 *          so a 2-gon has no interior for the distance to be measured against.
 */
inline void case_sdf_polygon_side_count() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::PlanarPolygon poly(b, opaque(0.5f), opaque(2), opaque(0.0f));
  if (poly.apothem == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: an angular repeat around a non-unit axis must trap.
 * @details SDF surface — the sector fold rotates the query point about the
 *          axis, so a non-unit one scales every folded copy off the sphere.
 */
inline void case_sdf_angular_repeat_nonunit_axis() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::Ring ring(b, opaque(1.0f), opaque(0.1f));
  math::Vector axis{opaque(0.0f), opaque(2.0f), opaque(0.0f)};
  SDF::AngularRepeat<SDF::Ring> rep(ring, opaque(4), axis); // non-unit -> trap
  if (rep.sector == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a knot ring with no cells must trap.
 * @details SDF surface — the per-pixel cell index divides the azimuth by
 *          2π/n, so n == 0 wraps to knots[-1] on every probe.
 */
inline void case_sdf_distorted_ring_zero_knots() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  const float knots[1] = {0.0f};
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(b, opaque(0.5f), opaque(0.05f), knots, opaque(0),
                          opaque(0.0f), pf); // no knot cells -> trap
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

inline void case_sdf_distorted_ring_one_knots() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  const float knots[1] = {0.0f};
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(b, opaque(0.5f), opaque(0.05f), knots, opaque(1),
                          opaque(0.0f), pf);
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

inline void case_sdf_distorted_ring_two_knots() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  const float knots[2] = {0.0f};
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(b, opaque(0.5f), opaque(0.05f), knots, opaque(2),
                          opaque(0.0f), pf);
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a twist warp around a zero-radius torus must trap.
 * @details SDF warp surface — the Lipschitz bound scales by 2/R, so a zero
 *          major radius hands the rasterizer a non-finite step bound.
 */
inline void case_sdf_twist_zero_major_radius() {
  SDF::Warp::Twist tw(opaque(2), opaque(0.1f), opaque(0.0f)); // R = 0 -> trap
  if (tw.two_over_r == opaque(42.0f))
    std::printf("x");
}

inline void case_transformed_torus_invalid_minor_radius() {
  SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{
      {opaque(1.0f), opaque(0.6f)}, {2, 0.1f, 1.0f}};
  Scan::TransformedVolume volume(torus, math::Vector(), math::Quaternion());
  volume.check_trace_preconditions();
}

/**
 * @brief Death case: a spherical polygon wider than a hemisphere must trap.
 * @details SDF surface — beyond the hemisphere the cap fold changes sign, so
 *          the shape must be built inverted about its antipode instead.
 */
inline void case_sdf_spherical_polygon_radius_over_hemisphere() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::SphericalPolygon poly(b, opaque(1.5f), opaque(5), opaque(0.0f));
  if (poly.circumradius == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a flower wider than a hemisphere must trap.
 * @details SDF surface — the petal cap bound is taken about the antipode, so a
 *          radius past the hemisphere inverts the band it derives.
 */
inline void case_sdf_flower_radius_over_hemisphere() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::Flower flower(b, opaque(1.5f), opaque(5), opaque(0.0f));
  if (flower.circumradius == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a zero-radius flower must trap.
 * @details SDF surface — the petal parameter divides by the circumradius, so a
 *          zero radius hands every probe a non-finite distance.
 */
inline void case_sdf_flower_zero_radius() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::Flower flower(b, opaque(0.0f), opaque(5), opaque(0.0f));
  if (flower.circumradius == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: baking a class LUT for a degenerate polygon must trap.
 * @details SDF class-LUT surface — fewer than three vertices leaves no closed
 *          boundary for the crossing test, so every sample would read as
 *          outside.
 */
inline void case_sdf_class_lut_too_few_vertices() {
  static const float poly_xy[4] = {-0.5f, -0.5f, 0.5f, -0.5f};
  static int16_t grid[16];
  SDF::ClassLut lut;
  SDF::build_canonical_distance_lut(poly_xy, opaque(2), opaque(4), grid, lut);
  if (lut.n == opaque(42))
    std::printf("x");
}

/**
 * @brief Death case: baking a class LUT on a single-cell grid must trap.
 * @details SDF class-LUT surface — the cell step divides by (n - 1), so a
 *          resolution below 2 makes the whole domain non-finite.
 */
inline void case_sdf_class_lut_grid_too_small() {
  static const float poly_xy[6] = {-0.5f, -0.5f, 0.5f, -0.5f, 0.0f, 0.5f};
  static int16_t grid[16];
  SDF::ClassLut lut;
  SDF::build_canonical_distance_lut(poly_xy, opaque(3), opaque(1), grid, lut);
  if (lut.n == opaque(42))
    std::printf("x");
}

/**
 * @brief Death case: binding a class LUT at a vertex offset outside the face
 *        must trap.
 * @details SDF class-LUT surface — the offset indexes the canonical polygon
 *          cyclically, so an out-of-range one correlates the face against
 *          storage past the shape.
 */
inline void case_sdf_bind_class_lut_offset_out_of_range() {
  constexpr int H = 16, HV = H + hs::H_OFFSET;
  math::Basis basis =
      math::make_basis(math::Quaternion(), math::Vector(0, 1, 0));
  math::Vector verts[3];
  uint16_t idx[3];
  for (int i = 0; i < 3; ++i) {
    float a = (2.0f * math::PI_F * i) / 3.0f;
    verts[i] = (basis.v * cosf(0.6f) +
                (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(0.6f))
                   .normalized();
    idx[i] = static_cast<uint16_t>(i);
  }
  static SDF::FaceScratchBuffer scratch;
  SDF::Face face(std::span<const math::Vector>(verts, 3),
                 std::span<const uint16_t>(idx, 3), scratch, HV, H);
  static const float canon_xy[6] = {-0.5f, -0.5f, 0.5f, -0.5f, 0.0f, 0.5f};
  SDF::ClassLut lut;
  if (face.bind_class_lut(&lut, canon_xy, opaque(7), false))
    std::printf("x");
}

/**
 * @brief Death case: a ring wider than the antipode must trap.
 * @details SDF ring surface — target_angle is radius * PI/2, so past 2 the
 *          band's cosine limits wrap and the stroke lands at the wrong
 *          latitude.
 */
inline void case_sdf_ring_radius_past_antipode() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::Ring ring(b, opaque(2.5f), opaque(0.05f));
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a ring with a negative stroke half-width must trap.
 * @details SDF ring surface — a negative thickness inverts the angular band,
 *          so every probe returns the far sentinel and the ring renders
 *          nothing.
 */
inline void case_sdf_ring_negative_thickness() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::Ring ring(b, opaque(1.0f), opaque(-0.05f));
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a distorted ring wider than the antipode must trap.
 * @details SDF ring surface — the shared ring geometry derives its band from
 *          radius * PI/2, which past 2 wraps its cosine limits.
 */
inline void case_sdf_distorted_ring_radius_past_antipode() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::FlatDistortedRing ring(b, opaque(2.5f), opaque(0.05f));
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a distorted ring with a negative half-width must trap.
 * @details SDF ring surface — a negative thickness inverts the angular band,
 *          so every probe returns the far sentinel and the ring renders
 *          nothing.
 */
inline void case_sdf_distorted_ring_negative_thickness() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::FlatDistortedRing ring(b, opaque(0.5f), opaque(-0.05f));
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a line with a negative stroke half-width must trap.
 * @details SDF line surface — a negative thickness inverts the angular band
 *          and shrinks the bounding cap below the arc's own half-length, so
 *          the cull drops rows the arc covers.
 */
inline void case_sdf_line_negative_thickness() {
  SDF::Line line(math::Vector(1, 0, 0), math::Vector(0, 0, 1), opaque(-0.05f));
  if (line.thickness == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a distorted ring built with a null shift callback must
 *        trap.
 * @details SDF ring surface — the callback is invoked per azimuth on every
 *          probe, so a null one faults deep inside the rasterizer instead.
 */
inline void case_sdf_distorted_ring_null_shift() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  ScalarFn shift; // default-constructed -> empty
  SDF::DistortedRing ring(b, opaque(0.5f), opaque(0.05f), shift, opaque(0.1f),
                          opaque(0.0f));
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

inline void case_sdf_distorted_ring_negative_distortion() {
  const math::Basis basis{math::X_AXIS, math::Y_AXIS, math::Z_AXIS};
  ScalarFn shift = [](float) { return 0.0f; };
  SDF::DistortedRing ring(basis, opaque(0.5f), opaque(0.05f), shift,
                          opaque(-0.1f), 0.0f);
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}
