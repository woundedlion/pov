/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/sdf/face.h.

// --- Congruence-class canonical distance LUTs --------------------------------
// Optional static-mesh LUTs; see face_class_bake.h.

/** Minimum squared normalized correlation for a valid class-LUT alignment;
 *  below this the face is too deformed and keeps the exact path. */
inline constexpr float ALIGN_MIN_CORR_SQ = 0.25f;
/** Maximum per-vertex deviation from the aligned canonical shape (as a multiple
 *  of the LUT cell diagonal) before a face keeps the exact path. The facility
 *  fits only meshes that hold still per spawn (face_class_bake.h). */
inline constexpr float ALIGN_MAX_DEV_DIAGS = 0.25f;

/**
 * @brief Canonical congruence-class signed-distance LUT, baked once per
 *        spawned mesh.
 * @details Distances are in canonical gnomonic plane units, quantized to int16
 * over the LUT box diameter (step ~1e-5 plane units). The domain is the
 * canonical polygon's bounding box + BOUNDS_MARGIN_WIDE, a different region
 * from the max_dist_sq cull disk (circumradius + the same margin): a probe can
 * survive the cull and still land outside the domain, so Face::distance's grid
 * clamp is required to keep the fetch in bounds.
 */
struct ClassLut {
  const int16_t *data =
      nullptr;          /**< n*n quantized signed distances (row-major). */
  int n = 0;            /**< Grid resolution per axis. */
  float cx = 0, cy = 0; /**< Canonical bounding-box center. */
  float Rx = 0, Ry = 0; /**< Half-extents (+ margin). */
  float inv_step_x = 0; /**< Reciprocal cell width. */
  float inv_step_y = 0; /**< Reciprocal cell height. */
  float safe_dist = 0;  /**< Cell diagonal (sign-pure interpolation bound). */
  float dequant = 0;    /**< int16 -> plane-unit scale. */
};

/**
 * @brief Bakes the signed point-to-polygon distance field of a canonical
 *        centered 2D polygon into an int16 grid.
 * @param poly_xy Centered polygon vertices, x/y pairs.
 * @param count Vertex count (>= 3).
 * @param n Grid resolution per axis (>= 2).
 * @param out Storage for n*n quantized samples.
 * @param lut Receives the domain/quantization parameters, with data = out.
 * @details Exact per-edge walk with crossing-test sign, over the bounding box
 * + BOUNDS_MARGIN_WIDE. That box is not a superset of the runtime cull disk
 * (circumradius + the same margin), so probes landing outside the domain rely
 * on Face::distance's grid clamp. Quantization scale is the box diameter (an
 * upper bound on any in-box distance: the polygon meets its own bounding box),
 * giving a step of ~1e-5 plane units — far below the interpolation bound.
 */
inline void build_canonical_distance_lut(const float *poly_xy, int count, int n,
                                         int16_t *out, ClassLut &lut) {
  HS_CHECK(count >= 3,
           "build_canonical_distance_lut requires at least 3 polygon vertices");
  HS_CHECK(n >= 2, "build_canonical_distance_lut requires a grid resolution "
                   "of at least 2");
  float bb_min_x = FLT_MAX, bb_max_x = -FLT_MAX;
  float bb_min_y = FLT_MAX, bb_max_y = -FLT_MAX;
  for (int i = 0; i < count; ++i) {
    float vx = poly_xy[2 * i], vy = poly_xy[2 * i + 1];
    bb_min_x = std::min(bb_min_x, vx);
    bb_max_x = std::max(bb_max_x, vx);
    bb_min_y = std::min(bb_min_y, vy);
    bb_max_y = std::max(bb_max_y, vy);
  }
  lut.cx = (bb_min_x + bb_max_x) * 0.5f;
  lut.cy = (bb_min_y + bb_max_y) * 0.5f;
  lut.Rx = std::max((bb_max_x - bb_min_x) * 0.5f + BOUNDS_MARGIN_WIDE, 0.01f);
  lut.Ry = std::max((bb_max_y - bb_min_y) * 0.5f + BOUNDS_MARGIN_WIDE, 0.01f);
  lut.n = n;
  lut.inv_step_x = (n - 1) / (2.0f * lut.Rx);
  lut.inv_step_y = (n - 1) / (2.0f * lut.Ry);
  float step_x = (2.0f * lut.Rx) / (n - 1);
  float step_y = (2.0f * lut.Ry) / (n - 1);
  // The plane SDF is 1-Lipschitz, so a zero anywhere in a cell puts every
  // corner within one cell diagonal of it; a min corner magnitude above the
  // diagonal guarantees a sign-pure cell (safe to interpolate).
  lut.safe_dist = sqrtf(step_x * step_x + step_y * step_y);
  float dmax = 2.0f * sqrtf(lut.Rx * lut.Rx + lut.Ry * lut.Ry);
  lut.dequant = dmax / 32767.0f;
  float quant = 32767.0f / dmax;

  for (int gy = 0; gy < n; ++gy) {
    float qy = (lut.cy - lut.Ry) + gy * step_y;
    for (int gx = 0; gx < n; ++gx) {
      float qx = (lut.cx - lut.Rx) + gx * step_x;
      float d_sq = FLT_MAX;
      bool inside = false;
      for (int i = 0; i < count; ++i) {
        float vx = poly_xy[2 * i], vy = poly_xy[2 * i + 1];
        int i2 = (i + 1 == count) ? 0 : i + 1;
        float ex = poly_xy[2 * i2] - vx, ey = poly_xy[2 * i2 + 1] - vy;
        float len_sq = ex * ex + ey * ey;
        float wx = qx - vx, wy = qy - vy;
        float t = len_sq > 1e-12f
                      ? hs::clamp((wx * ex + wy * ey) / len_sq, 0.0f, 1.0f)
                      : 0.0f;
        float bx = wx - ex * t, by = wy - ey * t;
        float dsq = bx * bx + by * by;
        if (dsq < d_sq)
          d_sq = dsq;
        if ((vy > qy) != (poly_xy[2 * i2 + 1] > qy)) {
          float ix = vx + (qy - vy) * ex / ey;
          if (qx < ix)
            inside = !inside;
        }
      }
      float d = (inside ? -1.0f : 1.0f) * sqrtf(d_sq);
      out[gy * n + gx] =
          static_cast<int16_t>(hs::clamp(d * quant, -32767.0f, 32767.0f));
    }
  }
  lut.data = out;
}

/**
 * @brief Accumulated complex correlation between a canonical polygon and a
 *        centered projection (see align_correlate).
 */
struct AlignCorr {
  float rr, ri; /**< Sum of canon_k * conj(z'_k) (real, imaginary). */
  float cc, zz; /**< Power terms: sum |canon_k|^2 and sum |z'_k|^2. */
};

/**
 * @brief Pairs canonical vertices with a centered 2D projection under a cyclic
 *        vertex offset + optional reflection.
 * @tparam GetZ Accessor invoked as get_z(j, zx, zy), returning the centered
 *         projection vertex j.
 * @tparam Visit Invoked as visit(k, zx, zy) for canonical vertex k and its
 *         paired projection vertex, already conjugated for the mirror family.
 * @param count Vertex count.
 * @param vert_offset Projection index corresponding to canonical vertex 0.
 * @param reflected Mirror family: conjugate each vertex and walk the
 *        projection indices in reverse (a mirrored face winds the opposite way
 *        in a consistently-wound mesh).
 * @param get_z Centered-projection vertex accessor.
 * @param visit Per-correspondence sink.
 * @details Single source for the correspondence convention — the correlation
 * below, bake-time clustering (face_class_bake.h) and the per-frame
 * Face::bind_class_lut all route through it, so the (offset, reflected)
 * encoding cannot drift.
 */
template <typename GetZ, typename Visit>
inline void align_walk(int count, int vert_offset, bool reflected, GetZ get_z,
                       Visit visit) {
  int j = vert_offset;
  for (int k = 0; k < count; ++k) {
    float zx, zy;
    get_z(j, zx, zy);
    if (reflected) {
      zy = -zy;
      if (--j < 0)
        j = count - 1;
    } else {
      if (++j == count)
        j = 0;
    }
    visit(k, zx, zy);
  }
}

/**
 * @brief Correlates a canonical polygon against a centered 2D projection under
 *        a cyclic vertex offset + optional reflection.
 * @tparam GetZ Accessor invoked as get_z(j, zx, zy), returning the centered
 *         projection vertex j.
 * @param canon_xy Canonical centered polygon, x/y pairs, canonical order.
 * @param count Vertex count.
 * @param vert_offset Projection index corresponding to canonical vertex 0.
 * @param reflected Mirror family (see align_walk).
 * @param get_z Centered-projection vertex accessor.
 * @return The correlation sums; the least-squares residual is
 *         cc + zz - 2*|r|, and the optimal rotation is r / |r|.
 */
template <typename GetZ>
inline AlignCorr align_correlate(const float *canon_xy, int count,
                                 int vert_offset, bool reflected, GetZ get_z) {
  AlignCorr a{0.0f, 0.0f, 0.0f, 0.0f};
  align_walk(count, vert_offset, reflected, get_z,
             [&](int k, float zx, float zy) {
               float cx = canon_xy[2 * k], cy = canon_xy[2 * k + 1];
               // canon * conj(z)
               a.rr += cx * zx + cy * zy;
               a.ri += cy * zx - cx * zy;
               a.cc += cx * cx + cy * cy;
               a.zz += zx * zx + zy * zy;
             });
  return a;
}
