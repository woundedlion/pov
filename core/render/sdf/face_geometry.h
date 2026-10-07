/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/sdf/face.h.

// ---------------------------------------------------------------------------
// Face construction, per-frame binding, bounds and distance queries.
// ---------------------------------------------------------------------------

/**
   * @brief Latitude-band reject for the face.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices.
   * @param geometry Display latitude mapping.
   * @param height Canvas height in rows.
   * @param bounds_margin Angular padding around the vertical bounds.
   * @return True when the phi extent plus AA margin maps to an empty
   *         canvas-row range.
   */
__attribute__((always_inline)) bool
compute_phi_extent(std::span<const math::Vector> vertices,
                   std::span<const uint16_t> indices,
                   const math::LatitudeGeometry &geometry, int height,
                   float bounds_margin) const {
  float min_y_val = 2.0f;
  float max_y_val = -2.0f;

  for (int idx : indices) {
    float y = vertices[idx].y;
    min_y_val = __builtin_fminf(y, min_y_val);
    max_y_val = __builtin_fmaxf(y, max_y_val);
  }

  float min_phi_check = math::fast_acos(hs::clamp(max_y_val, -1.0f, 1.0f));
  float max_phi_check = math::fast_acos(hs::clamp(min_y_val, -1.0f, 1.0f));
  Bounds rows =
      phi_bounds_to_rows(min_phi_check - bounds_margin,
                         max_phi_check + bounds_margin, geometry, height);

  return rows.y_min > rows.y_max;
}

/**
   * @brief Builds the local tangent frame, gnomonic 2D projection, and 3D
   * arrays.
   * @param vertices Shared vertex pool.
   * @param indices Indices selecting this face's vertices.
   * @param scratch Scratch storage receiving poly_2d and verts_3d.
   * @details Sets basis_u/v/w, poly_2d (with circumradius), and the local 3D
   * vertex array.
   */
__attribute__((always_inline)) void
setup_frame_and_polygon(std::span<const math::Vector> vertices,
                        std::span<const uint16_t> indices,
                        FaceScratchBuffer &scratch) {
  center = math::Vector(0, 0, 0);
  for (int i = 0; i < count; ++i) {
    const math::Vector &v = vertices[indices[i]];
    scratch.verts_3d[i] = v;
    center = center + v;
  }
  center.normalize();

  basis_v = center;
  basis_u = math::perpendicular_axis(center);
  basis_w = math::cross(center, basis_u).normalized();

  float max_r2 = 0.0f;
  for (int i = 0; i < count; ++i) {
    const math::Vector &v = scratch.verts_3d[i];
    // Gnomonic projection divides by d = cos(angle from face center),
    // singular at 90 degrees from the center; clamp d away from zero,
    // sign-preserving.
    float d = math::dot(v, basis_v);
    if (fabsf(d) < math::TOLERANCE)
      d = copysignf(math::TOLERANCE, d);
    float px = math::dot(v, basis_u) / d;
    float py = math::dot(v, basis_w) / d;

    scratch.poly_2d[i] = math::Vector(px, py, 0);

    float r2 = px * px + py * py;
    max_r2 = __builtin_fmaxf(r2, max_r2);
  }
  radius = sqrtf(max_r2);
  max_dist = radius + BOUNDS_MARGIN_WIDE;
  max_dist_sq = max_dist * max_dist;

  scratch.poly_2d[count] = scratch.poly_2d[0];
  poly_2d = std::span<math::Vector>(scratch.poly_2d.data(), count + 1);

  scratch.verts_3d[count] = scratch.verts_3d[0];
}

/**
   * @brief Computes the face "size" (inradius) from the projected polygon.
   * @param scratch Scratch storage holding poly_2d, per-edge vectors and
   * reciprocal squared lengths (inv_edge_lengths_sq) from compute_full_bounds,
   * which must run first.
   * @details size = minimum distance from the projected centroid to any edge,
   * floored to a fraction of the circumradius for degenerate slivers. A large
   * face converts it to radians, matching the metric distance() reports.
   */
__attribute__((always_inline)) void
compute_inradius(const FaceScratchBuffer &scratch) {
  float min_edge_dist = 1e9f;
  for (int i = 0; i < count; ++i) {
    const math::Vector &v1 = scratch.poly_2d[i];
    const math::Vector &edge = scratch.edge_vectors[i];
    float t = 0.0f;
    const float inv_edge_len_sq = scratch.inv_edge_lengths_sq[i];
    if (inv_edge_len_sq > 0.0f) {
      t = math::dot(-v1, edge) * inv_edge_len_sq;
      t = __builtin_fmaxf(0.0f, __builtin_fminf(1.0f, t));
    }
    math::Vector closest = v1 + edge * t;
    float d_line = closest.magnitude();

    min_edge_dist = __builtin_fminf(d_line, min_edge_dist);
  }
  const float min_radius =
      __builtin_fmaxf(0.0f, min_edge_dist - radius * 1e-5f);
  sector_min_radius_sq = min_radius * min_radius;
  size = __builtin_fmaxf(min_edge_dist, radius * MIN_SIZE_RADIUS_RATIO);
  linear_dist = size < 0.2f;
  // distance() reports radians for a large face; the same fast_atan2 keeps
  // dist/size exactly 1 at the inradius.
  if (!linear_dist)
    size = math::fast_atan2(size, 1.0f);
}

/**
   * @brief Packs per-edge data contiguously for the distance() fallback.
   * @param scratch Scratch storage receiving packed_edges.
   */
__attribute__((always_inline)) void pack_edges(FaceScratchBuffer &scratch) {
  for (int i = 0; i < count; ++i) {
    auto &ep = scratch.packed_edges[i];
    ep.vx = poly_2d[i].x;
    ep.vy = poly_2d[i].y;
    ep.ex = edge_vectors[i].x;
    ep.ey = edge_vectors[i].y;
    ep.inv_len_sq = inv_edge_lengths_sq[i];
    ep.inv_ej = inv_edge_j[i];
    ep.key_vy = angle_key(poly_2d[i].y);
    // A y-degenerate edge (inv_ej == 0) has no usable crossing x; equal keys
    // drop it from distance()'s parity test.
    ep.key_next_vy =
        (inv_edge_j[i] != 0.0f) ? angle_key(poly_2d[i + 1].y) : ep.key_vy;
  }
  packed_edges = std::span<EdgePacked>(scratch.packed_edges.data(), count);
}

/**
   * @brief Detects a convex 2D projection and builds its edge half-planes.
   * @param scratch Scratch storage receiving half_planes.
   * @param area2 Twice the polygon's signed area, from the collapsed-face cull.
   * @details The maximum edge half-plane distance is exact inside a convex
   * polygon and in exterior edge slabs, and underestimates in a vertex's
   * exterior normal cone. Concave, degenerate-edged, or wrongly-oriented
   * polygons leave convex false and distance() on the exact walk.
   */
__attribute__((always_inline)) void
build_half_planes(FaceScratchBuffer &scratch, float area2) {
  ::new (static_cast<void *>(&scratch.half_planes))
      std::array<HalfPlane, FaceScratchBuffer::MAX_VERTS>;
  bool pos = false, neg = false;
  const math::Vector *e1 = &edge_vectors[count - 1];
  float l1 = edge_lengths_sq[count - 1];
  for (int i = 0; i < count; ++i) {
    const math::Vector &e2 = edge_vectors[i];
    float cr = e1->x * e2.y - e1->y * e2.x;
    float scale = l1 * edge_lengths_sq[i];
    e1 = &e2;
    l1 = edge_lengths_sq[i];
    if (cr * cr > TURN_EPS_SQ * scale) {
      if (cr > 0)
        pos = true;
      else
        neg = true;
    }
  }
  if (pos && neg)
    return;

  float sign = area2 >= 0 ? 1.0f : -1.0f;
  // d0 is the origin's worst half-plane distance: the projected face center
  // must be strictly interior or the winding/orientation is untrustworthy.
  float d0 = -FLT_MAX;
  for (int i = 0; i < count; ++i) {
    float len_sq = edge_lengths_sq[i];
    if (len_sq < 1e-12f)
      return;
    float inv = sign / sqrtf(len_sq);
    float nx = edge_vectors[i].y * inv;
    float ny = -edge_vectors[i].x * inv;
    float off = -(nx * poly_2d[i].x + ny * poly_2d[i].y);
    auto &hp = scratch.half_planes[i];
    hp.nx = nx;
    hp.ny = ny;
    hp.off = off;
    d0 = __builtin_fmaxf(off, d0);
  }
  if (d0 >= 0.0f)
    return;
  half_planes = std::span<const HalfPlane>(scratch.half_planes.data(), count);
  convex = true;
}

/**
   * @brief Builds the angular sector table for the concave sector walk.
   * @param scratch Scratch storage receiving the unwrapped vertex
   * pseudo-angles.
   * @details Qualifies non-convex faces with at least SECTOR_MIN_COUNT
   * vertices that are star-shaped about the projected centroid (vertex
   * pseudo-angles strictly monotonic over one full turn, with every edge
   * facing the origin). Otherwise sector_ok stays false.
   */
__attribute__((always_inline)) void build_sectors(FaceScratchBuffer &scratch) {
  sector_ok = false;
  if (convex || count < SECTOR_MIN_COUNT)
    return;
  float prev = pseudo_angle(poly_2d[0].y, poly_2d[0].x);
  float acc = prev, total = 0.0f;
  scratch.pseudo_angles[0] = prev;
  for (int i = 1; i <= count; ++i) {
    float a = pseudo_angle(poly_2d[i].y, poly_2d[i].x);
    float d = a - prev;
    if (d > 2.0f)
      d -= 4.0f;
    if (d < -2.0f)
      d += 4.0f;
    acc += d;
    total += d;
    scratch.pseudo_angles[i] = acc;
    prev = a;
  }
  if (fabsf(fabsf(total) - 4.0f) > 1e-3f)
    return;
  float sgn = (total >= 0.0f) ? 1.0f : -1.0f;
  float min_step = FLT_MAX;
  float prev_s = scratch.pseudo_angles[0] * sgn;
  scratch.pseudo_angles[0] = prev_s;
  for (int i = 1; i <= count; ++i) {
    float cur = scratch.pseudo_angles[i] * sgn;
    scratch.pseudo_angles[i] = cur;
    min_step = __builtin_fminf(cur - prev_s, min_step);
    prev_s = cur;
  }
  if (min_step <= 0.0f)
    return;
  for (int i = 0; i < count; ++i) {
    const auto &a = poly_2d[i];
    const auto &b = poly_2d[i + 1];
    if ((a.x * b.y - a.y * b.x) * sgn <= 0.0f)
      return;
  }
  for (int i = 0; i <= count; ++i)
    scratch.sector_keys[i] = angle_key(scratch.pseudo_angles[i]);
  sector_keys =
      std::span<const uint32_t>(scratch.sector_keys.data(), count + 1);
  sector_base = scratch.pseudo_angles[0];
  sector_sgn = sgn;
  for (int i = 0; i < count; ++i)
    if (inv_edge_lengths_sq[i] == 0.0f) {
      sector_min_radius_sq = 0.0f;
      break;
    }
  for (int i = 0; i < count; ++i) {
    const auto &v = poly_2d[i];
    const float len_sq = v.x * v.x + v.y * v.y;
    const float scale = sqrtf(sector_min_radius_sq / len_sq);
    scratch.planes[i] = math::Vector(v.x * scale, v.y * scale, len_sq);
  }
  sector_rays = std::span<const math::Vector>(scratch.planes.data(), count);
  sector_ok = true;
}

/** @brief Builds vertex-row crossing masks for exact non-convex probes. */
__attribute__((always_inline)) void build_y_walk(FaceScratchBuffer &scratch) {
  if (convex || sector_ok || count < SECTOR_MIN_COUNT)
    return;
  ::new (static_cast<void *>(&scratch.y_walk)) FaceScratchBuffer::YWalkCache;
  auto &cache = scratch.y_walk;
  for (int i = 0; i < count; ++i)
    cache.indices[i] = static_cast<uint8_t>(i);
  std::sort(cache.indices.begin(), cache.indices.begin() + count,
            [&](uint8_t a, uint8_t b) { return poly_2d[a].y < poly_2d[b].y; });
  uint64_t active = 0;
  for (int first = 0; first < count;) {
    int last = first + 1;
    const float y = poly_2d[cache.indices[first]].y;
    while (last < count && poly_2d[cache.indices[last]].y == y)
      ++last;
    for (int j = first; j < last; ++j) {
      const int vertex = cache.indices[j];
      const int previous = vertex == 0 ? count - 1 : vertex - 1;
      if (poly_2d[vertex].y != poly_2d[vertex + 1].y)
        active ^= uint64_t{1} << vertex;
      if (poly_2d[previous].y != poly_2d[previous + 1].y)
        active ^= uint64_t{1} << previous;
    }
    for (int j = first; j < last; ++j) {
      scratch.pseudo_angles[j] = y;
      cache.masks[j] = active;
    }
    first = last;
  }
  y_coordinates = std::span<const float>(scratch.pseudo_angles.data(), count);
  y_masks = std::span<const uint64_t>(cache.masks.data(), count);
  y_indices = std::span<const uint8_t>(cache.indices.data(), count);
}

/**
   * @brief Aligns the current projection to a canonical class shape and binds
   *        its distance LUT.
   * @param lut Canonical-frame LUT for the face's congruence class.
   * @param canon_xy Canonical centered 2D polygon, x/y pairs.
   * @param vert_offset Cyclic offset aligning mesh vertex order to canonical.
   * @param reflected True for the mirror family.
   * @return False when the correlation is degenerate (badly deformed face);
   *         the face then keeps the exact path.
   * @details One complex correlation over the vertices recovers the in-plane
   * rotation placing the canonical shape at the face's least-squares pose.
   * Rotational-symmetry ambiguity is harmless (the LUT is invariant under the
   * shape's symmetry group).
   */
HS_COLD_MEMBER bool bind_class_lut(const ClassLut *lut, const float *canon_xy,
                                   int vert_offset, bool reflected) {
  HS_CHECK(vert_offset >= 0 && vert_offset < count,
           "bind_class_lut: vertex offset outside the face");
  float mx = 0.0f, my = 0.0f;
  for (int i = 0; i < count; ++i) {
    mx += poly_2d[i].x;
    my += poly_2d[i].y;
  }
  float inv_n = 1.0f / count;
  mx *= inv_n;
  my *= inv_n;

  AlignCorr a = align_correlate(canon_xy, count, vert_offset, reflected,
                                [&](int j, float &zx, float &zy) {
                                  zx = poly_2d[j].x - mx;
                                  zy = poly_2d[j].y - my;
                                });
  float r2 = a.rr * a.rr + a.ri * a.ri;
  if (r2 <= ALIGN_MIN_CORR_SQ * a.cc * a.zz)
    return false;
  float inv_r = 1.0f / sqrtf(r2);
  float c = a.rr * inv_r, s = a.ri * inv_r;

  // |d_true - d_canon| is bounded by the worst aligned vertex deviation, so
  // the sign-purity guard widens by that bound.
  float max_dev_sq = 0.0f;
  align_walk(
      count, vert_offset, reflected,
      [&](int j, float &zx, float &zy) {
        zx = poly_2d[j].x - mx;
        zy = poly_2d[j].y - my;
      },
      [&](int k, float zx, float zy) {
        float ex = canon_xy[2 * k] - (c * zx - s * zy);
        float ey = canon_xy[2 * k + 1] - (s * zx + c * zy);
        float dev_sq = ex * ex + ey * ey;
        if (dev_sq > max_dev_sq)
          max_dev_sq = dev_sq;
      });
  float max_dev = sqrtf(max_dev_sq);
  if (max_dev > ALIGN_MAX_DEV_DIAGS * lut->safe_dist)
    return false;

  // q = rot * (p - m), with the mirror family's conjugation folded into the
  // matrix; then the LUT grid transform (q - box_min) * inv_step folded on
  // top, so the probe loop runs a single affine map on raw (px, py).
  float m00 = c, m01 = reflected ? s : -s;
  float m10 = s, m11 = reflected ? -c : c;
  lut_ax = m00 * lut->inv_step_x;
  lut_bx = m01 * lut->inv_step_x;
  lut_cx = (-(m00 * mx + m01 * my) - lut->cx + lut->Rx) * lut->inv_step_x;
  lut_ay = m10 * lut->inv_step_y;
  lut_by = m11 * lut->inv_step_y;
  lut_cy = (-(m10 * mx + m11 * my) - lut->cy + lut->Ry) * lut->inv_step_y;
  lut_n = lut->n;
  lut_clamp = static_cast<float>(lut->n - 2);
  lut_dequant = lut->dequant;
  lut_q_safe = static_cast<int32_t>((lut->safe_dist + max_dev) / lut->dequant);
  lut_data = lut->data;
  return true;
}

/**
   * @brief Fills the scratch vertex azimuths the interval pass consumes.
   * @param scratch Scratch storage holding verts_3d and receiving thetas.
   */
__attribute__((always_inline)) void
compute_thetas(FaceScratchBuffer &scratch) const {
  for (int i = 0; i < count; ++i) {
    const math::Vector &v = scratch.verts_3d[i];
    float theta = math::fast_atan2(v.z, v.x);
    if (theta < 0)
      theta += math::TWO_PI_F;
    scratch.thetas[i] = theta;
  }
}

/**
   * @brief Computes the face's azimuth coverage intervals.
   * @param scratch Scratch storage holding thetas and receiving intervals.
   * @details Finds the largest angular gap between vertices; if it exceeds pi
   * the face does not wrap, so the complementary horizontal interval(s) are
   * emitted, else it spans full width. Coarse: only the single largest gap is
   * excised.
   */
__attribute__((always_inline)) void
compute_azimuth_intervals(FaceScratchBuffer &scratch) {
  float *th = scratch.thetas.data();
  for (int i = 1; i < count; ++i) {
    float t = th[i];
    float *p = th + i;
    while (p != th && p[-1] > t) {
      p[0] = p[-1];
      --p;
    }
    *p = t;
  }
  float max_gap = 0;
  float gap_start = 0;
  for (int i = 0; i < count; ++i) {
    float next = (i + 1 < count) ? scratch.thetas[i + 1]
                                 : (scratch.thetas[0] + math::TWO_PI_F);
    float diff = next - scratch.thetas[i];
    if (diff > max_gap) {
      max_gap = diff;
      gap_start = scratch.thetas[i];
    }
  }

  int interval_count = 0;
  if (max_gap > math::PI_F) {
    full_width = false;
    float start_t = fmodf(gap_start + max_gap, math::TWO_PI_F);
    // fmodf can leave start_t at ~2*PI instead of ~0, producing a degenerate
    // [~2*PI, 2*PI] sliver; snap to 0.
    if (start_t > math::TWO_PI_F - 1e-4f)
      start_t = 0.0f;
    float end_t = gap_start;

    if (start_t <= end_t) {
      scratch.intervals[interval_count++] = {start_t, end_t};
    } else {
      scratch.intervals[interval_count++] = {0.0f, end_t};
      scratch.intervals[interval_count++] = {start_t, math::TWO_PI_F};
    }
  } else {
    full_width = true;
  }
  intervals = std::span<Interval>(scratch.intervals.data(), interval_count);
}

/**
   * @brief Necessary condition for apply_pole_containment to fire.
   * @return True when a pole falls within the gnomonic circumcircle of the
   *         face's vertices (one test covers both poles); false rules out
   *         pole containment.
   */
__attribute__((always_inline)) bool pole_within_circumcircle() const {
  return center.y * center.y * (1.0f + radius * radius) >= 1.0f;
}

/**
   * @brief Ray-crossing test of a projected pole against the 2D face polygon.
   * @param ppx Projected pole x in the face's 2D basis.
   * @param ppy Projected pole y in the face's 2D basis.
   * @return true if (ppx, ppy) lies inside the polygon.
   */
__attribute__((always_inline)) bool pole_inside_polygon(float ppx,
                                                        float ppy) const {
  bool inside = false;
  for (int i = 0; i < count; ++i) {
    if (inv_edge_j[i] != 0.0f &&
        (poly_2d[i].y > ppy) != (poly_2d[i + 1].y > ppy)) {
      float ix = poly_2d[i].x +
                 (ppy - poly_2d[i].y) * edge_vectors[i].x * inv_edge_j[i];
      if (ppx < ix)
        inside = !inside;
    }
  }
  return inside;
}

/** How a projected pole sits against the 2D face polygon. */
enum class PoleHit {
  NONE,     /**< Outside the closed region. */
  BOUNDARY, /**< Within POLE_BOUNDARY_TOL of a vertex or an edge. */
  INTERIOR  /**< Strictly inside, off the boundary band. */
};

/**
   * @brief Classifies a projected pole against the 2D face polygon.
   * @param ppx Projected pole x in the face's 2D basis.
   * @param ppy Projected pole y in the face's 2D basis.
   * @return Where the pole falls relative to the closed polygon.
   * @details The boundary band is tested first: a pole exactly on a vertex or
   * edge would otherwise leave pole_inside_polygon's crossing parity to
   * rounding.
   */
PoleHit pole_hit(float ppx, float ppy) const {
  if (!pole_within_circumcircle())
    return PoleHit::NONE;
  const float tol = radius * POLE_BOUNDARY_TOL;
  const float tol_sq = tol * tol;
  for (int i = 0; i < count; ++i) {
    const float ax = ppx - poly_2d[i].x;
    const float ay = ppy - poly_2d[i].y;
    if (ax * ax + ay * ay <= tol_sq)
      return PoleHit::BOUNDARY;
    const float t =
        hs::clamp((ax * edge_vectors[i].x + ay * edge_vectors[i].y) *
                      inv_edge_lengths_sq[i],
                  0.0f, 1.0f);
    const float dx = ax - t * edge_vectors[i].x;
    const float dy = ay - t * edge_vectors[i].y;
    if (dx * dx + dy * dy <= tol_sq)
      return PoleHit::BOUNDARY;
  }
  return pole_inside_polygon(ppx, ppy) ? PoleHit::INTERIOR : PoleHit::NONE;
}

/**
   * @brief Extends the vertical bounds when the face reaches a pole.
   * @param height Canvas height in rows.
   * @details A face enclosing a pole wraps every azimuth, so it takes full
   * width outright. A face that only meets the pole on its boundary keeps its
   * azimuth wedge, which get_horizontal_intervals widens per row.
   */
__attribute__((always_inline)) void apply_pole_containment(int height) {
  if (center.y > 0.01f) {
    float inv_c = 1.0f / center.y;
    const PoleHit hit = pole_hit(basis_u.y * inv_c, basis_w.y * inv_c);
    if (hit != PoleHit::NONE) {
      y_min = 0;
      if (hit == PoleHit::INTERIOR)
        full_width = true;
    }
  }
  // South pole (0, -1, 0)
  if (center.y < -0.01f) {
    float inv_c = 1.0f / -center.y;
    const PoleHit hit = pole_hit(-basis_u.y * inv_c, -basis_w.y * inv_c);
    if (hit != PoleHit::NONE) {
      y_max = height - 1;
      if (hit == PoleHit::INTERIOR)
        full_width = true;
    }
  }
}

/**
   * @brief Refines phi bounds with an edge's great-circle arc extremum.
   * @param n Normalized great-circle plane normal of the edge.
   * @param v1 First edge endpoint (on the unit sphere).
   * @param v2 Second edge endpoint (on the unit sphere).
   * @param min_phi In/out running minimum phi.
   * @param max_phi In/out running maximum phi.
   * @details The extremum of an edge's arc may lie between its endpoints;
   * project the pole-tangent onto the plane and, if it falls inside the arc,
   * fold its phi into the bounds.
   */
static __attribute__((always_inline)) void
refine_phi_from_arc_extremum(const math::Vector &n, const math::Vector &v1,
                             const math::Vector &v2, float &min_phi,
                             float &max_phi) {
  float ny = n.y;
  if (std::abs(ny) < 0.99999f) {
    float nx = n.x;
    float nz = n.z;
    float tx = -nx * ny;
    float ty = 1.0f - ny * ny;
    float tz = -nz * ny;
    float t_len_sq = tx * tx + ty * ty + tz * tz;
    if (t_len_sq > 1e-12f) {
      float inv_len = 1.0f / sqrtf(t_len_sq);
      float ptx = tx * inv_len;
      float pty = ty * inv_len;
      float ptz = tz * inv_len;
      float cx1 = (v1.y * ptz - v1.z * pty) * nx +
                  (v1.z * ptx - v1.x * ptz) * ny +
                  (v1.x * pty - v1.y * ptx) * nz;
      float cx2 = (pty * v2.z - ptz * v2.y) * nx +
                  (ptz * v2.x - ptx * v2.z) * ny +
                  (ptx * v2.y - pty * v2.x) * nz;
      if (cx1 > 0 && cx2 > 0)
        min_phi = __builtin_fminf(math::fast_acos(hs::clamp(pty, -1.0f, 1.0f)),
                                  min_phi);
      if (cx1 < 0 && cx2 < 0)
        max_phi = __builtin_fmaxf(math::fast_acos(hs::clamp(-pty, -1.0f, 1.0f)),
                                  max_phi);
    }
  }
}

/**
   * @brief Snaps phi bounds to a pole when the face's planes enclose it.
   * @param scratch Scratch storage holding the compacted great-circle planes.
   * @param planes_count Number of valid planes.
   * @param center Normalized face centroid.
   * @param min_phi In/out phi minimum, set to 0 if the north pole is enclosed.
   * @param max_phi In/out phi maximum, set to PI if the south pole is enclosed.
   */
static __attribute__((always_inline)) void
snap_phi_for_pole_planes(const FaceScratchBuffer &scratch, int planes_count,
                         const math::Vector &center, float &min_phi,
                         float &max_phi) {
  bool np_inside = (planes_count > 0);
  bool sp_inside = (planes_count > 0);
  // Both flags only ever clear, so the scan is done once neither survives.
  for (int pi = 0; pi < planes_count && (np_inside || sp_inside); ++pi) {
    float py = scratch.planes[pi].y;
    bool center_pos = math::dot(center, scratch.planes[pi]) > 0;
    if ((py > 0) != center_pos)
      np_inside = false;
    if ((py < 0) != center_pos)
      sp_inside = false;
  }
  if (np_inside)
    min_phi = 0.0f;
  if (sp_inside)
    max_phi = math::PI_F;
}

/**
   * @brief Vertical bounds (arc extrema + pole-plane snap) for every face.
   * @details Fills the per-edge 2D vectors, squared lengths and reciprocals
   * consumed by later passes.
   * @param scratch Scratch storage holding poly_2d/verts_3d, receiving edge
   * data and planes.
   * @param count Vertex/edge count.
   * @param center Normalized face centroid.
   * @param geometry Display latitude mapping.
   * @param height Canvas height in rows.
   * @param y_min_out Output: first covered row.
   * @param y_max_out Output: last covered row.
   * @param bounds_margin Angular padding around the vertical bounds.
   */
HS_O3_FN static void compute_full_bounds(FaceScratchBuffer &scratch, int count,
                                         const math::Vector &center,
                                         const math::LatitudeGeometry &geometry,
                                         int height, int &y_min_out,
                                         int &y_max_out, float bounds_margin) {
  float min_phi = 100.0f;
  float max_phi = -100.0f;
  int planes_count = 0;
  for (int i = 0; i < count; ++i) {
    const math::Vector &v1 = scratch.verts_3d[i];
    const math::Vector &v2 = scratch.verts_3d[i + 1];
    math::Vector edge = scratch.poly_2d[i + 1] - scratch.poly_2d[i];
    scratch.edge_vectors[i] = edge;
    float edge_len_sq = math::dot(edge, edge);
    scratch.edge_lengths_sq[i] = edge_len_sq;
    scratch.inv_edge_lengths_sq[i] =
        (edge_len_sq > 1e-12f) ? (1.0f / edge_len_sq) : 0.0f;
    scratch.inv_edge_j[i] =
        (std::abs(edge.y) > 1e-12f) ? (1.0f / edge.y) : 0.0f;
    math::Vector normal = math::cross(v1, v2);
    float len_sq = math::dot(normal, normal);
    // planes[] is compacted: a degenerate edge pushes no plane, so planes[k]
    // is not edge k.
    if (len_sq > 1e-12f)
      scratch.planes[planes_count++] = normal.normalized();
    float phi_val = math::fast_acos(hs::clamp(v1.y, -1.0f, 1.0f));
    min_phi = __builtin_fminf(phi_val, min_phi);
    max_phi = __builtin_fmaxf(phi_val, max_phi);
    // Arc Extrema Logic: only when this edge pushed its own plane, else
    // planes[planes_count - 1] is a prior edge's normal against these
    // endpoints.
    if (len_sq > 1e-12f)
      refine_phi_from_arc_extremum(scratch.planes[planes_count - 1], v1, v2,
                                   min_phi, max_phi);
  }
  snap_phi_for_pole_planes(scratch, planes_count, center, min_phi, max_phi);
  Bounds rows = phi_bounds_to_rows(min_phi - bounds_margin,
                                   max_phi + bounds_margin, geometry, height);
  y_min_out = rows.y_min;
  y_max_out = rows.y_max;
}

/**
   * @brief Returns the face's precomputed inclusive row bounds.
   * @tparam H Canvas height in rows; must match the construction height the
   * bounds were computed for.
   * @return The stored {y_min, y_max} bounds.
   */
template <int H> Bounds get_vertical_bounds() const {
  HS_CHECK(H == build_height,
           "Face::get_vertical_bounds: H differs from construction height");
  HS_CHECK(!scratch_owner || scratch_owner->claim_seq == scratch_claim,
           "SDF::Face scanned after a later Face claimed its scratch buffer");
  return {y_min, y_max};
}

/**
   * @brief Reports whether rounded azimuth intervals change across a row band.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @param y_lo First row in the band.
   * @param y_hi Last row in the band.
   * @return True when the band requires per-row interval construction.
   */
template <int W, int H>
bool horizontal_intervals_vary_by_row(int y_lo, int y_hi) const {
  if (full_width || y_lo >= y_hi)
    return false;
  const float pad_lo = azimuth_pad_at_row<W, H>(y_lo);
  const float pad_hi = azimuth_pad_at_row<W, H>(y_hi);
  float narrow_pad = std::min(pad_lo, pad_hi);
  const float wide_pad = std::max(pad_lo, pad_hi);
  const float EQUATOR_ROW =
      math::DisplayGeometry<H>::phi_to_row(math::PI_F * 0.5f);
  const int EQUATOR_LO = static_cast<int>(floorf(EQUATOR_ROW));
  const int EQUATOR_HI = static_cast<int>(ceilf(EQUATOR_ROW));
  if (y_lo <= EQUATOR_LO && EQUATOR_LO <= y_hi)
    narrow_pad = std::min(narrow_pad, azimuth_pad_at_row<W, H>(EQUATOR_LO));
  if (y_lo <= EQUATOR_HI && EQUATOR_HI <= y_hi)
    narrow_pad = std::min(narrow_pad, azimuth_pad_at_row<W, H>(EQUATOR_HI));
  const float column_scale = W / math::TWO_PI_F;
  for (const auto &iv : intervals) {
    if (floorf((iv.start - narrow_pad) * column_scale) !=
            floorf((iv.start - wide_pad) * column_scale) ||
        ceilf((iv.end + narrow_pad) * column_scale) !=
            ceilf((iv.end + wide_pad) * column_scale))
      return true;
  }
  return false;
}

/** @brief Whether two rows emit the same rounded azimuth intervals. */
template <int W, int H>
bool horizontal_intervals_equal_rows(int first_y, int second_y) const {
  if (full_width || first_y == second_y)
    return true;
  const float first_pad = azimuth_pad_at_row<W, H>(first_y);
  const float second_pad = azimuth_pad_at_row<W, H>(second_y);
  if (first_pad == math::PI_F || second_pad == math::PI_F)
    return first_pad == second_pad;
  const float column_scale = W / math::TWO_PI_F;
  for (const auto &iv : intervals) {
    if (floorf((iv.start - first_pad) * column_scale) !=
            floorf((iv.start - second_pad) * column_scale) ||
        ceilf((iv.end + first_pad) * column_scale) !=
            ceilf((iv.end + second_pad) * column_scale))
      return false;
  }
  return true;
}

/**
   * @brief Emits the face's azimuth-coverage intervals for a row.
   * @tparam W Canvas width in columns; must match the clip width the
   * construction-time azimuth cull ran against.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y Row index, which sets the pole widening.
   * @param out Sink accepting (float start, float end).
   * @return True when the row was handled, possibly with no intervals outside
   *   the face's rows; false requests a full scan.
   * @details The pad is an azimuth angle, so it holds a whole pixel of AA reach
   * only at the equator; at colatitude phi one pad p of great-circle reach
   * subtends asin(p / sin phi), reaching the whole row once sin phi <= p.
   */
template <int W, int H, typename OutputIt>
bool get_horizontal_intervals(int y, OutputIt out) const {
  if (y < y_min || y > y_max)
    return true;
  HS_CHECK(build_width == 0 || W == build_width,
           "Face::get_horizontal_intervals: W differs from the clip width the "
           "azimuth cull ran against");
  if (full_width)
    return false;
  const float pad = azimuth_pad_at_row<W, H>(y);
  if (pad == math::PI_F)
    return false;
  const float column_scale = W / math::TWO_PI_F;
  for (const auto &iv : intervals) {
    float f_x1 = (iv.start - pad) * column_scale;
    float f_x2 = (iv.end + pad) * column_scale;
    out(floorf(f_x1), ceilf(f_x2));
  }
  return true;
}

/** @brief Returns the latitude-adjusted AA padding for one raster row. */
template <int W, int H> float azimuth_pad_at_row(int y) const {
  if (build_azimuth_pads)
    return build_azimuth_pads[y];
  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  return face_azimuth_pad(W, math::TrigLUT<W, H>::sin_phi[y]);
}

/**
   * @brief Signed planar distance via the convex half-plane max.
   * @param px Gnomonic x of the query point.
   * @param py Gnomonic y of the query point.
   * @return Signed distance in the tangent plane (negative inside).
   */
HS_O3_FN float plane_dist_convex(float px, float py) const {
  HS_SCAN_METRIC(hs::g_scan_metrics.convex_hits++);
  float d = -FLT_MAX;
  for (int i = 0; i < count; ++i) {
    const auto &hp = half_planes[i];
    float di = hp.nx * px + hp.ny * py + hp.off;
    d = __builtin_fmaxf(di, d);
  }
  return d;
}

/**
   * @brief Squared planar distance via the exact per-edge walk.
   * @param px Gnomonic x of the query point.
   * @param py Gnomonic y of the query point.
   * @param inside_out Set true when the query lies inside the polygon; carries
   * the sign the squared return cannot.
   * @return Squared distance to the nearest edge, in the tangent plane.
   */
HS_O3_FN float plane_dsq_exact(float px, float py, bool &inside_out) const {
  if (!y_coordinates.empty())
    return plane_dsq_y_walk(px, py, inside_out);
  float d = FLT_MAX;
  bool inside = false;
  const uint32_t qk = angle_key(py);
  for (int i = 0; i < count; ++i) {
    const auto &ep = packed_edges[i];
    float wx = px - ep.vx, wy = py - ep.vy;
    float t = (wx * ep.ex + wy * ep.ey) * ep.inv_len_sq;
    float cv = hs::clamp(t, 0.0f, 1.0f);
    float bx = wx - ep.ex * cv, by = wy - ep.ey * cv;
    float dsq = bx * bx + by * by;
    d = __builtin_fminf(dsq, d);
    if ((ep.key_vy > qk) != (ep.key_next_vy > qk)) {
      float isx = ep.vx + (py - ep.vy) * ep.ex * ep.inv_ej;
      if (px < isx)
        inside = !inside;
    }
  }
  inside_out = inside;
  return d;
}

/**
   * @brief Exact squared distance and parity by sorted vertex-row traversal.
   * @details Checks all edges crossing the query row, then incident edges at
   * neighboring vertex rows. Unvisited segments lie outside the visited strip;
   * its vertical gaps bound their distance. Requires a populated y-walk cache.
   */
HS_O3_FN float plane_dsq_y_walk(float px, float py, bool &inside_out) const {
  float d = FLT_MAX;
  uint64_t visited = 0;
  auto edge_dsq = [&](int i) {
    const uint64_t bit = uint64_t{1} << i;
    if (visited & bit)
      return;
    visited |= bit;
    const auto &ep = packed_edges[i];
    const float wx = px - ep.vx, wy = py - ep.vy;
    const float t =
        hs::clamp((wx * ep.ex + wy * ep.ey) * ep.inv_len_sq, 0.0f, 1.0f);
    const float bx = wx - ep.ex * t, by = wy - ep.ey * t;
    d = __builtin_fminf(d, bx * bx + by * by);
  };
  int lo = -1, hi = count;
  while (lo + 1 < hi) {
    const int mid = (lo + hi) >> 1;
    if (y_coordinates[mid] <= py)
      lo = mid;
    else
      hi = mid;
  }
  bool inside = false;
  uint64_t active = lo < 0 ? 0 : y_masks[lo];
  while (active) {
    const int i = std::countr_zero(active);
    active &= active - 1;
    edge_dsq(i);
    const auto &ep = packed_edges[i];
    if (ep.key_vy != ep.key_next_vy) {
      const float isx = ep.vx + (py - ep.vy) * ep.ex * ep.inv_ej;
      if (px < isx)
        inside = !inside;
    }
  }
  inside_out = inside;
  int left = lo, right = lo + 1;
  while (left >= 0 || right < count) {
    const float left_gap = left >= 0 ? py - y_coordinates[left] : FLT_MAX;
    const float right_gap = right < count ? y_coordinates[right] - py : FLT_MAX;
    const bool left_far = left < 0 || left_gap * left_gap > d * (1.0f + 1e-5f);
    const bool right_far =
        right >= count || right_gap * right_gap > d * (1.0f + 1e-5f);
    if (left_far && right_far)
      break;
    if (!left_far) {
      const int vertex = y_indices[left--];
      edge_dsq(vertex);
      edge_dsq(vertex == 0 ? count - 1 : vertex - 1);
    }
    if (!right_far) {
      const int vertex = y_indices[right++];
      edge_dsq(vertex);
      edge_dsq(vertex == 0 ? count - 1 : vertex - 1);
    }
  }
  return d;
}

/**
   * @brief Squared planar distance via the concave sector walk.
   * @param px Gnomonic x of the query point.
   * @param py Gnomonic y of the query point.
   * @param inside_out Set true when the query lies inside the polygon; carries
   * the sign the squared return cannot.
   * @return Squared distance to the nearest edge, in the tangent plane.
   * @details Searches the fan sector and its neighbors, certifying the minimum
   * against the omitted edges' angular wedge and radial bound, expanding until
   * certified or every segment has been evaluated. Requires sector_ok.
   */
HS_O3_FN float plane_dsq_sector(float px, float py, bool &inside_out) const {
  float p = pseudo_angle(py, px) * sector_sgn;
  if (p < sector_base)
    p += 4.0f;
  uint32_t qk = angle_key(p);
  int lo = 0, hi = count;
  while (lo + 1 < hi) {
    int mid = (lo + hi) >> 1;
    if (sector_keys[mid] <= qk)
      lo = mid;
    else
      hi = mid;
  }
  int s = lo;
  const auto &sector_start = poly_2d[s];
  const auto &sector_end = poly_2d[s + 1];
  const float ax = sector_start.x * py, ay = sector_start.y * px;
  const float bx = px * sector_end.y, by = py * sector_end.x;
  if ((ax - ay) * sector_sgn <= (fabsf(ax) + fabsf(ay)) * 1e-6f ||
      (bx - by) * sector_sgn <= (fabsf(bx) + fabsf(by)) * 1e-6f)
    return plane_dsq_exact(px, py, inside_out);

  float d = FLT_MAX;
  auto edge_dsq = [&](int idx) {
    const auto &ep = packed_edges[idx];
    float wx = px - ep.vx, wy = py - ep.vy;
    float t = hs::clamp((wx * ep.ex + wy * ep.ey) * ep.inv_len_sq, 0.0f, 1.0f);
    float bx = wx - ep.ex * t, by = wy - ep.ey * t;
    d = __builtin_fminf(bx * bx + by * by, d);
  };
  int first = s - 1;
  if (first < 0)
    first += count;
  int last = s + 2;
  if (last >= count)
    last -= count;
  edge_dsq(first);
  edge_dsq(s);
  edge_dsq(s + 1 == count ? 0 : s + 1);
  int remaining = count - 3;
  // Unvisited segments lie beyond the minimum-radius disk in this arc.
  auto ray_excludes = [&](int vertex) {
    const auto &v = poly_2d[vertex];
    const auto &ray = sector_rays[vertex];
    const float rx = ray.x, ry = ray.y;
    if (px * rx + py * ry < sector_min_radius_sq) {
      const float dx = px - rx, dy = py - ry;
      return d < (dx * dx + dy * dy) * (1.0f - 1e-5f);
    }
    if (px * v.x + py * v.y <= 0.0f)
      return d < (px * px + py * py) * (1.0f - 1e-5f);
    const float cross = px * v.y - py * v.x;
    return d * ray.z < cross * cross * (1.0f - 1e-5f);
  };
  bool left_excluded = false, right_excluded = false;
  while (remaining > 0) {
    left_excluded = left_excluded || ray_excludes(first);
    right_excluded = right_excluded || ray_excludes(last);
    if (left_excluded && right_excluded)
      break;
    if (!left_excluded) {
      if (--first < 0)
        first += count;
      edge_dsq(first);
      if (--remaining == 0)
        break;
    }
    if (!right_excluded) {
      edge_dsq(last);
      if (++last == count)
        last = 0;
      --remaining;
    }
  }
  const auto &edge = packed_edges[s];
  inside_out =
      (edge.ex * (py - edge.vy) - edge.ey * (px - edge.vx)) * sector_sgn >=
      0.0f;
  return d;
}

/** @brief Face has a congruence-class distance LUT (bind_class_lut). */
static constexpr uint32_t PROBE_HAS_LUT = 1u << 0;
/** @brief Use convex half-plane distance. */
static constexpr uint32_t PROBE_CONVEX = 1u << 1;
/** @brief Use the sector-indexed edge walk. */
static constexpr uint32_t PROBE_SECTOR = 1u << 2;
/** @brief Return distances in the gnomonic plane. */
static constexpr uint32_t PROBE_LINEAR = 1u << 3;

/** @return Packed distance-path flags for repeated probes of this face. */
uint32_t probe_flags() const {
  return (lut_data ? PROBE_HAS_LUT : 0u) | (convex ? PROBE_CONVEX : 0u) |
         (sector_ok ? PROBE_SECTOR : 0u) | (linear_dist ? PROBE_LINEAR : 0u);
}

/**
   * @brief Computes signed distance to the face, writing into res.
   * @tparam ComputeUVs Accepted for interface parity; the face stores no UVs.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = raw_dist = the signed edge distance,
   *        size = inradius, both in gnomonic plane units on a linear_dist
   *        face and in radians otherwise.
   * @param reject_dsq Squared plane distance at/above which an outside probe
   *        is reported as the far sentinel without taking the square root.
   *        Must be conservative: only probes the caller rejects on dist may
   *        cross it. FLT_MAX disables the cull. Consulted only on a
   *        linear_dist face's edge-walk path; large faces and the convex
   *        half-plane path ignore it.
   * @note Distances live in the face's gnomonic tangent plane. Small faces
   *       (linear_dist) report the plane distance directly; large faces convert
   *       via fast_atan2(plane, 1). Do not treat raw_dist as a metric geodesic
   *       angle.
   */
template <bool ComputeUVs = true>
HS_O3_FN void distance(const math::Vector &p, DistanceResult &res,
                       float reject_dsq = FLT_MAX) const {
  distance_with_flags<ComputeUVs>(p, res, reject_dsq, probe_flags());
}

/**
   * @brief Evaluates distance with cached flags, recomputing the radial cull cosine.
   * @tparam ComputeUVs Accepted for interface parity; the face stores no UVs.
   * @param p Unit query direction.
   * @param res Output distance result.
   * @param reject_dsq Conservative squared plane-distance rejection threshold.
   * @param probe_flags Flags captured by probe_flags() after geometry/LUT updates.
   * @details The five-argument overload accepts a cached cosine for repeated
   * probes.
   */
template <bool ComputeUVs = true>
HS_O3_FN void distance_with_flags(const math::Vector &p, DistanceResult &res,
                                  float reject_dsq,
                                  uint32_t probe_flags) const {
  const float min_cos = 1.0f / sqrtf(1.0f + max_dist_sq);
  distance_with_flags<ComputeUVs>(p, res, reject_dsq, probe_flags, min_cos);
}

/**
   * @brief Computes distance using flags captured by probe_flags().
   * @tparam ComputeUVs Accepted for interface parity; the face stores no UVs.
   * @param p Point on sphere (normalized).
   * @param res Output distance result.
   * @param reject_dsq Conservative squared distance rejection threshold,
   *        honored only when PROBE_LINEAR is set and the probe lands outside;
   *        the PROBE_CONVEX path ignores it.
   * @param probe_flags Distance-path flags captured after the face's last LUT
   *        binding or geometry update.
   * @param min_cos Cosine of the radial cull angle.
   */
template <bool ComputeUVs = true>
HS_O3_FN void distance_with_flags(const math::Vector &p, DistanceResult &res,
                                  float reject_dsq, uint32_t probe_flags,
                                  float min_cos) const {
  HS_SCAN_METRIC(hs::g_scan_metrics.pixels_tested++);
  HS_PROBE_TICK();
  HS_PROBE_COUNT(n_probe);
  HS_PROBE_MARK(hs_t);

  float cos_angle = math::dot(p, center);
  if (cos_angle < min_cos) {
    HS_SCAN_METRIC(hs::g_scan_metrics.pixels_culled++);
    HS_PROBE_SPAN(point, hs_t);
    HS_PROBE_COUNT(n_cull_r);
    res = DistanceResult(FAR_SENTINEL, 0.0f, FAR_SENTINEL, 0.0f, size);
    return;
  }
  HS_PROBE_SPAN(point, hs_t);

  float inv_cos = 1.0f / cos_angle;
  float px = math::dot(p, basis_u) * inv_cos;
  float py = math::dot(p, basis_w) * inv_cos;
  HS_PROBE_SPAN(project, hs_t);

  float plane_dist;
  bool lut_served = false;
  if (probe_flags & PROBE_HAS_LUT) {
    // Affine map into the canonical LUT grid, then a 4-tap bilinear fetch.
    // Only sign-pure cells at least one cell diagonal from the boundary are
    // served; the AA fringe and sign-unsafe cells fall back to the edge walk
    // on the TRUE per-frame edges.
    float fx = lut_ax * px + lut_bx * py + lut_cx;
    float fy = lut_ay * px + lut_by * py + lut_cy;
    // Clamp is memory safety: the cull disk is not contained in the LUT box,
    // so a surviving probe can map outside the grid.
    fx = hs::clamp(fx, 0.0f, lut_clamp);
    fy = hs::clamp(fy, 0.0f, lut_clamp);
    int ix = (int)fx;
    int iy = (int)fy;
    const int16_t *cell = lut_data + iy * lut_n + ix;
    int32_t q00 = cell[0], q10 = cell[1];
    int32_t q01 = cell[lut_n], q11 = cell[lut_n + 1];
    // Sign-purity + magnitude guard in quantized integer units.
    int32_t all_or = q00 | q10 | q01 | q11;
    int32_t all_and = q00 & q10 & q01 & q11;
    int32_t min_q =
        std::min({std::abs(q00), std::abs(q10), std::abs(q01), std::abs(q11)});
    if ((all_or >= 0 || all_and < 0) && min_q > lut_q_safe) {
      HS_SCAN_METRIC(hs::g_scan_metrics.lut_hits++);
      float tx = fx - ix;
      float ty = fy - iy;
      float d0 = q00 + (q10 - q00) * tx;
      float d1 = q01 + (q11 - q01) * tx;
      plane_dist = lut_dequant * (d0 + (d1 - d0) * ty);
      lut_served = true;
      HS_PROBE_SPAN(edge_lut, hs_t);
      HS_PROBE_COUNT(n_lut);
    }
  }
  if (!lut_served) {
    HS_SCAN_METRIC(hs::g_scan_metrics.exact_hits++);
    if (probe_flags & PROBE_CONVEX) {
      plane_dist = plane_dist_convex(px, py);
      HS_PROBE_SPAN(edge_convex, hs_t);
      HS_PROBE_COUNT(n_convex);
    } else {
      bool inside;
      float dsq;
      if (probe_flags & PROBE_SECTOR) {
        HS_SCAN_METRIC(hs::g_scan_metrics.sector_hits++);
        dsq = plane_dsq_sector(px, py, inside);
        HS_PROBE_SPAN(edge_sector, hs_t);
        HS_PROBE_COUNT(n_sector);
      } else {
        dsq = plane_dsq_exact(px, py, inside);
        HS_PROBE_SPAN(edge_exact, hs_t);
        HS_PROBE_COUNT(n_exact);
      }
      // Outside probes past the (margin-carrying) reject bound skip the
      // sqrt; the caller rejects them on dist alone.
      if ((probe_flags & PROBE_LINEAR) && !inside && dsq >= reject_dsq) {
        res = DistanceResult(FAR_SENTINEL, 0.0f, FAR_SENTINEL, 0.0f, size);
        HS_PROBE_SPAN(pack, hs_t);
        return;
      }
      plane_dist = (inside ? -1.0f : 1.0f) * sqrtf(dsq);
    }
  }

  // Small faces skip the plane->angle conversion: tan(angle) ~ angle.
  float raw = (probe_flags & PROBE_LINEAR) ? plane_dist
                                           : math::fast_atan2(plane_dist, 1.0f);
  res = DistanceResult(raw, 0.0f, raw, 0.0f, size);
  HS_PROBE_SPAN(pack, hs_t);
}
