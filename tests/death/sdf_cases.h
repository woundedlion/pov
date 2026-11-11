/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// SDF death cases.

/**
 * @brief Death case: a polygon with fewer than three sides must trap.
 */
inline void case_sdf_polygon_side_count() {
  SDF::PlanarPolygon poly(DEATH_IDENTITY_BASIS, opaque(0.5f), opaque(2),
                          opaque(0.0f));
  (void)poly;
}

/**
 * @brief Death case: an angular repeat around a non-unit axis must trap.
 */
inline void case_sdf_angular_repeat_nonunit_axis() {
  SDF::Ring ring(DEATH_IDENTITY_BASIS, opaque(1.0f), opaque(0.1f));
  math::Vector axis{opaque(0.0f), opaque(2.0f), opaque(0.0f)};
  SDF::AngularRepeat<SDF::Ring> rep(ring, opaque(4), axis); // non-unit -> trap
  (void)rep;
}

/**
 * @brief Death case: a knot ring with no cells must trap.
 */
inline void case_sdf_distorted_ring_zero_knots() {
  const float KNOTS[1] = {0.0f};
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(DEATH_IDENTITY_BASIS, opaque(0.5f), opaque(0.05f),
                          KNOTS, opaque(0), opaque(0.0f),
                          pf); // no knot cells -> trap
  (void)ring;
}

inline void case_sdf_distorted_ring_one_knots() {
  const float KNOTS[1] = {0.0f};
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(DEATH_IDENTITY_BASIS, opaque(0.5f), opaque(0.05f),
                          KNOTS, opaque(1), opaque(0.0f), pf);
  (void)ring;
}

inline void case_sdf_distorted_ring_two_knots() {
  const float KNOTS[2] = {0.0f};
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(DEATH_IDENTITY_BASIS, opaque(0.5f), opaque(0.05f),
                          KNOTS, opaque(2), opaque(0.0f), pf);
  (void)ring;
}

/**
 * @brief Death case: a twist warp around a zero-radius torus must trap.
 */
inline void case_sdf_twist_zero_major_radius() {
  SDF::Warp::Twist tw(opaque(2), opaque(0.1f), opaque(0.0f)); // R = 0 -> trap
  (void)tw;
}

inline void case_transformed_torus_invalid_minor_radius() {
  SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{
      {opaque(1.0f), opaque(0.6f)}, {2, 0.1f, 1.0f}};
  Scan::TransformedVolume volume(torus, math::Vector(), math::Quaternion());
  volume.check_trace_preconditions();
}

/**
 * @brief Death case: a spherical polygon wider than a hemisphere must trap.
 * @details Build such a shape inverted about its antipode instead.
 */
inline void case_sdf_spherical_polygon_radius_over_hemisphere() {
  SDF::SphericalPolygon poly(DEATH_IDENTITY_BASIS, opaque(1.5f), opaque(5),
                             opaque(0.0f));
  (void)poly;
}

/**
 * @brief Death case: a flower wider than a hemisphere must trap.
 */
inline void case_sdf_flower_radius_over_hemisphere() {
  SDF::Flower flower(DEATH_IDENTITY_BASIS, opaque(1.5f), opaque(5),
                     opaque(0.0f));
  (void)flower;
}

/**
 * @brief Death case: a zero-radius flower must trap.
 */
inline void case_sdf_flower_zero_radius() {
  SDF::Flower flower(DEATH_IDENTITY_BASIS, opaque(0.0f), opaque(5),
                     opaque(0.0f));
  (void)flower;
}

/**
 * @brief Death case: baking a class LUT for a degenerate polygon must trap.
 */
inline void case_sdf_class_lut_too_few_vertices() {
  static const float POLY_XY[4] = {-0.5f, -0.5f, 0.5f, -0.5f};
  static int16_t grid[16];
  SDF::ClassLut lut;
  SDF::build_canonical_distance_lut(POLY_XY, opaque(2), opaque(4), grid, lut);
  (void)lut;
}

/**
 * @brief Death case: baking a class LUT on a single-cell grid must trap.
 */
inline void case_sdf_class_lut_grid_too_small() {
  static const float POLY_XY[6] = {-0.5f, -0.5f, 0.5f, -0.5f, 0.0f, 0.5f};
  static int16_t grid[16];
  SDF::ClassLut lut;
  SDF::build_canonical_distance_lut(POLY_XY, opaque(3), opaque(1), grid, lut);
  (void)lut;
}

/**
 * @brief Death case: binding a class LUT at a vertex offset outside the face
 *        must trap.
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
  SDF::PreparedFace face(std::span<const math::Vector>(verts, 3),
                         std::span<const uint16_t>(idx, 3), scratch, HV, H,
                         nullptr, nullptr, SDF::BOUNDS_MARGIN);
  static const float CANON_XY[6] = {-0.5f, -0.5f, 0.5f, -0.5f, 0.0f, 0.5f};
  SDF::ClassLut lut;
  (void)face.bind_class_lut(&lut, CANON_XY, opaque(7), false);
}

/**
 * @brief Death case: a ring wider than the antipode must trap.
 */
inline void case_sdf_ring_radius_past_antipode() {
  SDF::Ring ring(DEATH_IDENTITY_BASIS, opaque(2.5f), opaque(0.05f));
  (void)ring;
}

/**
 * @brief Death case: a ring with a negative stroke half-width must trap.
 */
inline void case_sdf_ring_negative_thickness() {
  SDF::Ring ring(DEATH_IDENTITY_BASIS, opaque(1.0f), opaque(-0.05f));
  (void)ring;
}

/**
 * @brief Death case: a distorted ring wider than the antipode must trap.
 */
inline void case_sdf_distorted_ring_radius_past_antipode() {
  SDF::FlatDistortedRing ring(DEATH_IDENTITY_BASIS, opaque(2.5f),
                              opaque(0.05f));
  (void)ring;
}

/**
 * @brief Death case: a distorted ring with a negative half-width must trap.
 */
inline void case_sdf_distorted_ring_negative_thickness() {
  SDF::FlatDistortedRing ring(DEATH_IDENTITY_BASIS, opaque(0.5f),
                              opaque(-0.05f));
  (void)ring;
}

/**
 * @brief Death case: a line with a negative stroke half-width must trap.
 */
inline void case_sdf_line_negative_thickness() {
  SDF::Line line(math::Vector(1, 0, 0), math::Vector(0, 0, 1), opaque(-0.05f));
  (void)line;
}

/**
 * @brief Death case: a distorted ring built with a null shift callback must
 *        trap.
 * @details The callback is invoked per azimuth on every probe.
 */
inline void case_sdf_distorted_ring_null_shift() {
  ScalarFn shift; // default-constructed -> empty
  SDF::DistortedRing ring(DEATH_IDENTITY_BASIS, opaque(0.5f), opaque(0.05f),
                          shift, opaque(0.1f), opaque(0.0f));
  (void)ring;
}

inline void case_sdf_distorted_ring_negative_distortion() {
  ScalarFn shift = [](float) { return 0.0f; };
  SDF::DistortedRing ring(DEATH_IDENTITY_BASIS, opaque(0.5f), opaque(0.05f),
                          shift, opaque(-0.1f), 0.0f);
  (void)ring;
}

/** @brief Rejects a non-positive smooth-union blend width. */
inline void case_sdf_smooth_union_zero_radius() {
  const SDF::PlanarPolygon polygon(DEATH_IDENTITY_BASIS, .5f, 3, 0.f);
  SDF::SmoothUnion<SDF::PlanarPolygon, SDF::PlanarPolygon> shape(
      polygon, polygon, opaque(0.f));
  (void)shape;
}

/** @brief Rejects an angular repeat with no copies. */
inline void case_sdf_angular_repeat_zero_copies() {
  const SDF::Ring ring(DEATH_IDENTITY_BASIS, 1.f, .1f);
  SDF::AngularRepeat<SDF::Ring> shape(ring, opaque(0), math::Y_AXIS);
  (void)shape;
}

inline void case_framework_invalid_geometry() {
  SDF::TriangularFramework geometry;
  geometry.cell_size = opaque(0.0f);
  SDF::FrameworkEvents events(geometry, {{0, 0, 0}, {1, 0, 0}, {0, 1}});
  (void)events.active(0);
}

inline void case_framework_invalid_ray() {
  SDF::FrameworkEvents events({}, {{0, 0, 0}, {opaque(0.0f), 0, 0}, {0, 1}});
  (void)events.active(0);
}

inline void case_octet_invalid_geometry() {
  SDF::OctetFramework geometry;
  geometry.cell_size = opaque(0.0f);
  SDF::OctetEvents events(geometry, {{0, 0, 0}, {1, 0, 0}, {0, 1}});
  (void)events.active(0);
}

inline void case_octet_invalid_ray() {
  SDF::OctetEvents events(SDF::OctetFramework{},
                          {{0, 0, 0}, {opaque(0.0f), 0, 0}, {0, 1}});
  (void)events.active(0);
}

inline void case_octet4_invalid_geometry() {
  SDF::OctetFramework4 geometry;
  geometry.wire_radius = opaque(0.0f);
  SDF::OctetEvents4 events(geometry, {}, {{1, 0, 0, 0}}, {0, 1});
  (void)events.active(0);
}

inline void case_octet4_invalid_interval() {
  SDF::OctetEvents4 events({}, {}, {{1, 0, 0, 0}}, {opaque(2.0f), 1});
  (void)events.active(0);
}

inline void case_octet4_nonfinite_ray() {
  SDF::OctetEvents4 events(
      {}, {{opaque(std::numeric_limits<float>::infinity()), 0, 0, 0}},
      {{1, 0, 0, 0}}, {0, 1});
  (void)events.active(0);
}

inline void case_octet4_nonunit_direction() {
  SDF::OctetEvents4 events({}, {}, {{opaque(2.0f), 0, 0, 0}}, {0, 1});
  (void)events.active(0);
}

/** @brief Early latitude rejection invalidates old prepared face handles. */
inline void case_face_handle_retargeted_by_early_cull() {
  constexpr int W = 96, H = 20;
  const math::Vector verts[] = {math::pixel_to_vector<W, H>(20.0f, 2.0f),
                                math::pixel_to_vector<W, H>(25.0f, 3.0f),
                                math::pixel_to_vector<W, H>(18.0f, 4.0f)};
  const uint16_t INDICES[] = {0, 1, 2};
  static SDF::FaceStorage storage;
  SDF::Face first(verts, INDICES, storage, math::LatitudeGeometry(H), H,
                  nullptr, nullptr, SDF::BOUNDS_MARGIN);
  const SDF::Face copy = first;
  const uint32_t table_claim = storage.tables.claim_seq;
  ClipRegion clip{.y_start = H - 3,
                  .y_end = H,
                  .x_start = 0,
                  .x_end = W,
                  .margin = 0,
                  .w = W,
                  .h = H};
  SDF::Face culled(verts, INDICES, storage, math::LatitudeGeometry(H), H, &clip,
                   nullptr, SDF::BOUNDS_MARGIN);
  const auto bounds = culled.get_vertical_bounds<H>();
  if (!(bounds.y_min > bounds.y_max) || storage.tables.claim_seq != table_claim)
    std::abort();
  (void)copy.get_vertical_bounds<H>();
}

/** @brief Mutable geometry binding rejects a retargeted face handle. */
inline void case_face_handle_mutable_retargeted() {
  constexpr int H = 20;
  const math::Vector VERTS[] = {math::X_AXIS, math::Y_AXIS, math::Z_AXIS};
  const uint16_t INDICES[] = {0, 1, 2};
  static SDF::FaceStorage storage;
  SDF::Face first(VERTS, INDICES, storage, math::LatitudeGeometry(H), H,
                  nullptr, nullptr, SDF::BOUNDS_MARGIN);
  SDF::Face second(VERTS, INDICES, storage, math::LatitudeGeometry(H), H,
                   nullptr, nullptr, SDF::BOUNDS_MARGIN);
  (void)second;
  (void)first.prepared();
}

/**
 * @brief Death case: scanning a Face whose scratch buffer a later Face claimed.
 * @details The second build retargets the first face's spans, so the first no
 *          longer describes its own geometry.
 */
inline void case_face_scratch_retargeted() {
  constexpr int H = 16, HV = H + hs::H_OFFSET;
  math::Basis basis =
      math::make_basis(math::Quaternion(), math::Vector(0, 1, 0));
  math::Vector verts[6];
  uint16_t idx_a[3], idx_b[3];
  for (int i = 0; i < 6; ++i) {
    float a = (2.0f * math::PI_F * i) / 6.0f;
    verts[i] = (basis.v * cosf(0.6f) +
                (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(0.6f))
                   .normalized();
  }
  for (int i = 0; i < 3; ++i) {
    idx_a[i] = static_cast<uint16_t>(i);
    idx_b[i] = static_cast<uint16_t>(i + 3);
  }
  static SDF::FaceScratchBuffer scratch;
  SDF::PreparedFace first(std::span<const math::Vector>(verts, 6),
                          std::span<const uint16_t>(idx_a, 3), scratch, HV, H,
                          nullptr, nullptr, SDF::BOUNDS_MARGIN);
  SDF::PreparedFace second(std::span<const math::Vector>(verts, 6),
                           std::span<const uint16_t>(idx_b, 3), scratch, HV, H,
                           nullptr, nullptr, SDF::BOUNDS_MARGIN);
  (void)second;
  (void)first.get_vertical_bounds<H>();
}

/** @brief Death case: a virtual-height Face must use the active display grid. */
inline void case_face_virtual_height_mismatched_geometry() {
  constexpr int H = 16;
  const math::Vector VERTS[] = {math::X_AXIS, math::Y_AXIS, math::Z_AXIS};
  const uint16_t INDICES[] = {0, 1, 2};
  static SDF::FaceScratchBuffer scratch;
  SDF::PreparedFace face(VERTS, INDICES, scratch, opaque(H + hs::H_OFFSET + 1),
                         H, nullptr, nullptr, SDF::BOUNDS_MARGIN);
  (void)face;
}

/** @brief A seeing-time face bounds pose requires a render clip. */
inline void case_face_bounds_pose_without_clip() {
  constexpr int H = 16;
  const math::Vector VERTS[] = {math::X_AXIS, math::Y_AXIS, math::Z_AXIS};
  const uint16_t INDICES[] = {0, 1, 2};
  static SDF::FaceScratchBuffer scratch;
  const SDF::FaceBoundsPose POSE{};
  SDF::PreparedFace face(VERTS, INDICES, scratch, math::LatitudeGeometry(H), H,
                         nullptr, nullptr, SDF::BOUNDS_MARGIN, &POSE);
  (void)face;
}

/**
 * @brief Death case: scanning a Face whose scratch buffer a later, culled Face
 *        claimed.
 * @details The second build fills the scratch buffer and only then culls on
 *          collapsed area, so the first face's spans are retargeted even though
 *          nothing will be drawn for the second.
 */
inline void case_face_scratch_retargeted_by_culled_face() {
  constexpr int H = 16, HV = H + hs::H_OFFSET;
  math::Basis basis =
      math::make_basis(math::Quaternion(), math::Vector(0, 1, 0));
  math::Vector verts[4];
  uint16_t idx_a[3], idx_b[3];
  for (int i = 0; i < 3; ++i) {
    float a = (2.0f * math::PI_F * i) / 3.0f;
    verts[i] = (basis.v * cosf(0.6f) +
                (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(0.6f))
                   .normalized();
    idx_a[i] = static_cast<uint16_t>(i);
    // Coincident vertices enclose no area: culled after the scratch fill.
    idx_b[i] = 3;
  }
  verts[3] = basis.v;
  static SDF::FaceScratchBuffer scratch;
  SDF::PreparedFace first(std::span<const math::Vector>(verts, 4),
                          std::span<const uint16_t>(idx_a, 3), scratch, HV, H,
                          nullptr, nullptr, SDF::BOUNDS_MARGIN);
  SDF::PreparedFace second(std::span<const math::Vector>(verts, 4),
                           std::span<const uint16_t>(idx_b, 3), scratch, HV, H,
                           nullptr, nullptr, SDF::BOUNDS_MARGIN);
  (void)second;
  (void)first.get_vertical_bounds<H>();
}
