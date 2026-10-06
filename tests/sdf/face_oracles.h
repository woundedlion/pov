/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Face distance vs an independent exact oracle
//
// Face distances match an independent oracle within 1e-4 inside and on
// concave faces; convex exteriors stay within [0, oracle] with that tolerance.
// ============================================================================

/**
 * @brief Exact signed point-to-polygon distance in a face's tangent plane.
 * @param verts Face vertices on the unit sphere, in ring order.
 * @param p Probe on the unit sphere.
 * @param linear_dist Whether the face reports plane units rather than radians.
 * @return The value Face::distance must reproduce: negative inside, mapped
 *         through fast_atan2 unless the face carries linear distance.
 * @details Rebuilds the gnomonic frame from the vertices alone (normalized
 *          centroid as the projection axis), independent of the face's own
 *          projection, edge packing and basis. Per-edge segment distance plus a
 *          crossing-parity inside test.
 */
inline float exact_plane_distance(std::span<const math::Vector> verts,
                                  const math::Vector &p, bool linear_dist) {
  math::Vector center(0, 0, 0);
  for (const math::Vector &v : verts)
    center = center + v;
  center.normalize();
  const math::Vector u =
      (verts[0] - center * math::dot(verts[0], center)).normalized();
  const math::Vector w = math::cross(center, u).normalized();
  auto project = [&](const math::Vector &v, float &x, float &y) {
    const float d = math::dot(v, center);
    x = math::dot(v, u) / d;
    y = math::dot(v, w) / d;
  };
  float px, py;
  project(p, px, py);
  const size_t n = verts.size();
  float dmin = FLT_MAX;
  bool inside = false;
  for (size_t i = 0; i < n; ++i) {
    float ax, ay, bx, by;
    project(verts[i], ax, ay);
    project(verts[(i + 1) % n], bx, by);
    const float ex = bx - ax, ey = by - ay;
    const float wx = px - ax, wy = py - ay;
    const float len_sq = ex * ex + ey * ey;
    const float t = len_sq > 0.0f
                        ? hs::clamp((wx * ex + wy * ey) / len_sq, 0.0f, 1.0f)
                        : 0.0f;
    const float cx = wx - ex * t, cy = wy - ey * t;
    dmin = std::min(dmin, cx * cx + cy * cy);
    if ((ay > py) != (by > py) && px < ax + (py - ay) * ex / ey)
      inside = !inside;
  }
  const float plane_exact = (inside ? -1.0f : 1.0f) * sqrtf(dmin);
  return linear_dist ? plane_exact : math::fast_atan2(plane_exact, 1.0f);
}

/**
 * @brief Builds one tilted N-gon face and scans its gnomonic box against an exact oracle.
 * @param sample_total Accumulates the number of non-culled samples evaluated.
 * @param sides Number of polygon points (must be <= 8).
 * @param rho Angular circumradius of the face, in radians.
 * @param axis Pole direction the face's basis is built around.
 * @param rho_inner Inner-vertex radius; > 0 interleaves star points (concave).
 * @details A convex face must match the oracle exactly inside and stay within
 * [0, oracle] outside (the half-plane path underestimates in vertex cones); a
 * concave face must match the oracle everywhere via the exact walk.
 */
inline void check_face_distance_oracle(int &sample_total, int sides, float rho,
                                       const math::Vector &axis,
                                       float rho_inner = 0.0f) {
  constexpr int H = 144;
  constexpr int HV = H + hs::H_OFFSET;
  HS_EXPECT_TRUE(sides >= 3 && sides <= 8);
  if (sides < 3 || sides > 8)
    return;

  math::Basis basis = math::make_basis(math::Quaternion(), axis);
  math::Vector verts3d[16];
  uint16_t idx[16];
  const int n_verts = rho_inner > 0.0f ? 2 * sides : sides;
  for (int i = 0; i < n_verts; ++i) {
    float a = (2.0f * math::PI_F * i) / n_verts + 0.37f;
    float r = (rho_inner > 0.0f && (i & 1)) ? rho_inner : rho;
    verts3d[i] =
        (basis.v * cosf(r) + (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(r))
            .normalized();
    idx[i] = static_cast<uint16_t>(i);
  }

  SDF::FaceScratchBuffer scratch;
  SDF::Face face(std::span<const math::Vector>(verts3d, n_verts),
                 std::span<const uint16_t>(idx, n_verts), scratch, HV, H);
  HS_EXPECT_EQ(face.convex, rho_inner <= 0.0f);
  const uint32_t probe_flags = face.probe_flags();

  int samples = 0;
  // Gnomonic point normalize(center + u*px + w*py): distance() recovers (px,py)
  // exactly after dividing by dot(p,center).
  const float reach = face.max_dist * 0.98f;
  constexpr int G = 64;
  for (int gi = 0; gi <= G; ++gi) {
    for (int gj = 0; gj <= G; ++gj) {
      float px = -reach + (2.0f * reach) * gi / G;
      float py = -reach + (2.0f * reach) * gj / G;
      math::Vector p =
          (face.basis_v + face.basis_u * px + face.basis_w * py).normalized();

      hs::g_scan_metrics.exact_hits = 0;
      SDF::DistanceResult res = SDF::distance_of(face, p);
      SDF::DistanceResult cached_res;
      face.distance_with_flags(p, cached_res, FLT_MAX, probe_flags);
      HS_EXPECT_EQ(std::memcmp(&res, &cached_res, sizeof(res)), 0);
      if (hs::g_scan_metrics.exact_hits == 0)
        continue; // culled (outside max_dist / behind the center)

      const float expected = exact_plane_distance(
          std::span<const math::Vector>(verts3d, n_verts), p, face.linear_dist);

      // fast_atan2 preserves sign, so the oracle's sign is the plane sign.
      if (face.convex && expected > 0.0f) {
        // Outside a convex face the half-plane max is a lower bound (line
        // distance, not vertex distance) that never crosses zero.
        HS_EXPECT_TRUE(res.raw_dist >= -1e-4f);
        HS_EXPECT_TRUE(res.raw_dist <= expected + 1e-4f);
      } else {
        HS_EXPECT_NEAR(res.raw_dist, expected, 1e-4f);
      }
      ++samples;
    }
  }
  sample_total += samples;
  HS_EXPECT_GT(samples, 100);
}

/** @brief Checks inside/outside signs across a backtracking sector. */
inline void test_face_sector_backtrack_sign() {
  int checked = 0;
  for (float bend : {-0.08f, -0.04f, 0.0f, 0.04f, 0.08f}) {
    math::Vector vertices[12];
    uint16_t indices[12];
    for (int i = 0; i < 12; ++i) {
      const float angle = (i == 1 ? bend : i * math::TWO_PI_F / 12.0f);
      const float radius = (i & 1) ? 0.24f : 0.5f;
      vertices[i] =
          math::Vector(radius * cosf(angle), radius * sinf(angle), 1.0f)
              .normalized();
      indices[i] = static_cast<uint16_t>(i);
    }
    SDF::FaceScratchBuffer scratch;
    SDF::Face face(std::span<const math::Vector>(vertices, 12),
                   std::span<const uint16_t>(indices, 12), scratch,
                   144 + hs::H_OFFSET, 144);
    if (!face.sector_ok || face.sector_kmax != 2)
      continue;
    for (int x = -100; x <= 100; ++x)
      for (int y = -100; y <= 100; ++y) {
        const float px = x * 0.005f;
        const float py = y * 0.005f;
        bool exact_inside, sector_inside;
        const float exact = face.plane_dsq_exact(px, py, exact_inside);
        const float sector = face.plane_dsq_sector(px, py, sector_inside);
        if (exact < 0.0004f)
          HS_EXPECT_NEAR(sector, exact, 1e-6f);
        if (exact > 1e-10f)
          HS_EXPECT_EQ(sector_inside, exact_inside);
      }
    ++checked;
  }
  HS_EXPECT_GT(checked, 0);
}

/**
 * @brief Verifies Face::distance agrees with the polygon oracle inside and bounds it outside.
 * @details Drives check_face_distance_oracle across a spread of polygons
 *   (triangle, pentagon, hexagon and concave stars) and tilts.
 */
inline void test_face_distance_matches_exact_oracle() {
  int samples = 0;
  check_face_distance_oracle(samples, /*sides=*/3, 0.45f,
                             math::Vector(0.4f, 0.3f, 1.0f));
  check_face_distance_oracle(samples, /*sides=*/5, 0.50f,
                             math::Vector(0.4f, 0.3f, 1.0f));
  check_face_distance_oracle(samples, /*sides=*/6, 0.40f,
                             math::Vector(-0.6f, 0.5f, 0.7f));
  check_face_distance_oracle(samples, /*sides=*/3, 0.12f,
                             math::Vector(0.4f, 0.3f, 1.0f));
  check_face_distance_oracle(samples, /*sides=*/4, 0.50f,
                             math::Vector(0.4f, 0.3f, 1.0f),
                             /*rho_inner=*/0.25f);
  // Concave star: convexity detection must reject it and the exact walk must
  // reproduce the oracle everywhere.
  check_face_distance_oracle(samples, /*sides=*/6, 0.50f,
                             math::Vector(0.4f, 0.3f, 1.0f),
                             /*rho_inner=*/0.25f);
  // The grid actually exercised the distance path.
  HS_EXPECT_GT(samples, 1000);
}

// ============================================================================
// Face congruence-class LUT vs the exact oracle
//
// Face::distance with a bound ClassLut serves sign-pure probes >= one cell
// diagonal from the boundary via a bilinear lookup in the canonical class
// frame; everything else falls back to the exact walk on the true edges.
// ============================================================================

/**
 * @brief Rotates a vector about a unit axis by an angle (Rodrigues).
 * @param v Vector to rotate.
 * @param k Unit rotation axis.
 * @param theta Rotation angle (radians).
 * @return The rotated vector.
 */
inline math::Vector rotate_about(const math::Vector &v, const math::Vector &k,
                                 float theta) {
  float c = cosf(theta), s = sinf(theta);
  return v * c + math::cross(k, v) * s + k * (math::dot(k, v) * (1.0f - c));
}

/**
 * @brief Bakes a canonical LUT from one concave star face and sweeps a
 *        transformed congruent copy against the exact oracle.
 * @param lut_total Accumulates the number of samples served by the LUT path.
 * @param cyc Cyclic shift applied to the copy's vertex order.
 * @param reflected Mirror the copy (and reverse its order, preserving winding).
 * @param rot_angle 3D rotation applied to the copy's vertices.
 */
inline void check_face_class_lut(int &lut_total, int cyc, bool reflected,
                                 float rot_angle) {
  constexpr int H = 144;
  constexpr int HV = H + hs::H_OFFSET;
  constexpr int sides = 6, n_verts = 2 * sides;
  constexpr float rho = 0.40f, rho_inner = 0.20f;
  const math::Vector axis(0.4f, 0.3f, 1.0f);

  math::Basis basis = math::make_basis(math::Quaternion(), axis);
  math::Vector orig[n_verts];
  for (int i = 0; i < n_verts; ++i) {
    float a = (2.0f * math::PI_F * i) / n_verts + 0.37f;
    float r = (i & 1) ? rho_inner : rho;
    orig[i] =
        (basis.v * cosf(r) + (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(r))
            .normalized();
  }

  // Canonical polygon: the untransformed face's own centered 2D projection.
  SDF::FaceScratchBuffer canon_scratch;
  uint16_t canon_idx[n_verts];
  for (int i = 0; i < n_verts; ++i)
    canon_idx[i] = static_cast<uint16_t>(i);
  SDF::Face canon_face(std::span<const math::Vector>(orig, n_verts),
                       std::span<const uint16_t>(canon_idx, n_verts),
                       canon_scratch, HV, H);
  HS_EXPECT_TRUE(!canon_face.convex);
  float canon[2 * n_verts];
  float mx = 0.0f, my = 0.0f;
  for (int i = 0; i < n_verts; ++i) {
    mx += canon_face.poly_2d[i].x;
    my += canon_face.poly_2d[i].y;
  }
  mx /= n_verts;
  my /= n_verts;
  for (int i = 0; i < n_verts; ++i) {
    canon[2 * i] = canon_face.poly_2d[i].x - mx;
    canon[2 * i + 1] = canon_face.poly_2d[i].y - my;
  }

  static int16_t lut_data[64 * 64];
  SDF::ClassLut lut;
  SDF::build_canonical_distance_lut(canon, n_verts, 64, lut_data, lut);
  // Cell diagonal of the 64x64 grid over this star's box, ~0.023; the sweep
  // tolerance is a multiple of it.
  HS_EXPECT_GT(lut.safe_dist, 0.0f);
  HS_EXPECT_LT(lut.safe_dist, 0.03f);

  // Congruent copy: rotate the sphere, then reindex. A cyclic shift by cyc
  // maps canonical vertex k to copy index (n - cyc + k) % n; the mirror family
  // reverses the order (keeping winding consistent) and binds with offset 0.
  math::Vector verts[n_verts];
  const math::Vector rot_axis = math::Vector(0.3f, -0.8f, 0.52f).normalized();
  for (int j = 0; j < n_verts; ++j) {
    int src = reflected ? (n_verts - j) % n_verts : (j + cyc) % n_verts;
    math::Vector v = rotate_about(orig[src], rot_axis, rot_angle);
    if (reflected)
      v.x = -v.x;
    verts[j] = v;
  }
  int off = reflected ? 0 : (n_verts - cyc) % n_verts;

  SDF::FaceScratchBuffer scratch;
  SDF::Face face(std::span<const math::Vector>(verts, n_verts),
                 std::span<const uint16_t>(canon_idx, n_verts), scratch, HV, H);
  HS_EXPECT_TRUE(face.bind_class_lut(&lut, canon, off, reflected));
  const uint32_t probe_flags = face.probe_flags();

  // A degenerate canonical shape must be rejected by the correlation guard.
  {
    SDF::FaceScratchBuffer reject_scratch;
    static const float zeros[2 * n_verts] = {};
    SDF::Face reject_face(std::span<const math::Vector>(verts, n_verts),
                          std::span<const uint16_t>(canon_idx, n_verts),
                          reject_scratch, HV, H);
    HS_EXPECT_TRUE(!reject_face.bind_class_lut(&lut, zeros, 0, false));
  }

  int lut_samples = 0, sign_mismatches = 0;
  float min_lut_mag = FLT_MAX;
  const float reach = face.max_dist * 0.98f;
  constexpr int G = 96;
  for (int gi = 0; gi <= G; ++gi) {
    for (int gj = 0; gj <= G; ++gj) {
      float px = -reach + (2.0f * reach) * gi / G;
      float py = -reach + (2.0f * reach) * gj / G;
      math::Vector p =
          (face.basis_v + face.basis_u * px + face.basis_w * py).normalized();

      hs::g_scan_metrics.lut_hits = 0;
      hs::g_scan_metrics.exact_hits = 0;
      SDF::DistanceResult res = SDF::distance_of(face, p);
      SDF::DistanceResult cached_res;
      face.distance_with_flags(p, cached_res, FLT_MAX, probe_flags);
      HS_EXPECT_EQ(std::memcmp(&res, &cached_res, sizeof(res)), 0);
      bool took_lut = hs::g_scan_metrics.lut_hits > 0;
      if (!took_lut && hs::g_scan_metrics.exact_hits == 0)
        continue; // culled

      const float expected = exact_plane_distance(
          std::span<const math::Vector>(verts, n_verts), p, face.linear_dist);

      if (took_lut) {
        ++lut_samples;
        if ((res.raw_dist < 0.0f) != (expected < 0.0f))
          ++sign_mismatches;
        float mag = std::abs(res.raw_dist);
        if (mag < min_lut_mag)
          min_lut_mag = mag;
        // Bilinear over a sign-pure cell of exact corner samples stays within
        // ~2 cell diagonals of the 1-Lipschitz field; slack covers the int16
        // quantization and the fast_atan2 conversion.
        HS_EXPECT_NEAR(res.raw_dist, expected, 3.0f * lut.safe_dist);
      } else {
        HS_EXPECT_NEAR(res.raw_dist, expected, 1e-4f);
      }
    }
  }
  HS_EXPECT_GT(lut_samples, 100);
  HS_EXPECT_EQ(sign_mismatches, 0);
  // The sign-purity guard keeps every served magnitude a cell diagonal from
  // zero — outside the AA ramp.
  if (lut_samples > 0) {
    const float plane_floor = std::max(
        lut.safe_dist, static_cast<float>(face.lut_q_safe) * face.lut_dequant);
    const float floor_mag =
        face.linear_dist ? plane_floor : math::fast_atan2(plane_floor, 1.0f);
    HS_EXPECT_GE(min_lut_mag, floor_mag - 1e-5f);
  }
  lut_total += lut_samples;
}

/**
 * @brief Verifies the class-LUT path across identity, cyclic-offset + 3D
 *        rotation, and mirror-family bindings.
 */
inline void test_face_class_lut_matches_oracle() {
  int lut_samples = 0;
  check_face_class_lut(lut_samples, /*cyc=*/0, /*reflected=*/false, 0.0f);
  check_face_class_lut(lut_samples, /*cyc=*/5, /*reflected=*/false, 1.1f);
  check_face_class_lut(lut_samples, /*cyc=*/0, /*reflected=*/true, 0.7f);
  // The sweep actually fired the LUT path.
  HS_EXPECT_GT(lut_samples, 300);
}
