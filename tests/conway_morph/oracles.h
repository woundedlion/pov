/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Shared oracles
// ---------------------------------------------------------------------------

/**
 * @brief Largest absolute deviation of the mesh's face-loop edge lengths from
 *        their mean.
 * @param m Mesh whose edges are measured (each shared edge counted twice; the
 *          uniform double count leaves mean and max deviation unchanged).
 * @return max_e |len(e) - mean_len|.
 */
inline float max_edge_length_deviation(const PolyMesh &m) {
  double sum = 0.0;
  int n = 0;
  size_t off = 0;
  for (size_t fi = 0; fi < m.face_counts.size(); ++fi) {
    const int c = m.face_counts[fi];
    for (int k = 0; k < c; ++k) {
      sum += math::distance_between(m.vertices[m.faces[off + k]],
                                    m.vertices[m.faces[off + (k + 1) % c]]);
      ++n;
    }
    off += c;
  }
  const float mean = static_cast<float>(sum / n);
  float worst = 0.0f;
  off = 0;
  for (size_t fi = 0; fi < m.face_counts.size(); ++fi) {
    const int c = m.face_counts[fi];
    for (int k = 0; k < c; ++k) {
      const float d =
          math::distance_between(m.vertices[m.faces[off + k]],
                                 m.vertices[m.faces[off + (k + 1) % c]]) -
          mean;
      worst = fold_worst(worst, std::abs(d));
    }
    off += c;
  }
  return worst;
}

/**
 * @brief Verifies got's vertices merge pairwise onto want's: exactly two got
 *        vertices within tol of every want vertex.
 */
inline void check_pairwise_vertex_cover(const PolyMesh &got,
                                        const PolyMesh &want, float tol) {
  HS_EXPECT_EQ(got.vertices.size(), 2 * want.vertices.size());
  for (size_t i = 0; i < want.vertices.size(); ++i) {
    HS_CONTEXT("want vertex", static_cast<long long>(i));
    int merged = 0;
    for (size_t j = 0; j < got.vertices.size(); ++j) {
      if ((got.vertices[j] - want.vertices[i]).length() <= tol)
        ++merged;
    }
    HS_EXPECT_EQ(merged, 2);
  }
}

/**
 * @brief Checks identical topology and nearest-vertex index correspondence.
 * @return False when mesh sizes differ; true when all comparisons ran.
 */
inline bool check_vertex_order_identity(const PolyMesh &got,
                                        const PolyMesh &src) {
  HS_EXPECT_EQ(got.vertices.size(), src.vertices.size());
  HS_EXPECT_EQ(got.face_counts.size(), src.face_counts.size());
  HS_EXPECT_EQ(got.faces.size(), src.faces.size());
  if (got.vertices.size() != src.vertices.size() ||
      got.face_counts.size() != src.face_counts.size() ||
      got.faces.size() != src.faces.size())
    return false;
  HS_EXPECT_EQ(std::memcmp(got.face_counts.data(), src.face_counts.data(),
                           got.face_counts.size() * sizeof(uint8_t)),
               0);
  HS_EXPECT_EQ(std::memcmp(got.faces.data(), src.faces.data(),
                           got.faces.size() * sizeof(uint16_t)),
               0);

  for (size_t i = 0; i < got.vertices.size(); ++i) {
    size_t nearest = 0;
    float best = 1e9f;
    for (size_t j = 0; j < src.vertices.size(); ++j) {
      const float d = math::distance_between(got.vertices[i], src.vertices[j]);
      if (d < best) {
        best = d;
        nearest = j;
      }
    }
    HS_CONTEXT("got vertex", static_cast<long long>(i));
    HS_EXPECT_EQ(nearest, i);
  }
  return true;
}

/** Corner-match radius at a truncate's T_EPS end, where each seed corner has
 * split into two cut vertices. Above the widest T_EPS cut on the registry
 * seeds and far below their closest corner spacing. */
constexpr float PRIMARY_CORNER_TOL_TRUNCATE = 0.08f;
/** Same radius for the expand/snub/chamfer eps ends, whose single corner per
 * source moves less. */
constexpr float PRIMARY_CORNER_TOL_SINGLE = 0.06f;

/**
 * @brief Verifies the output's primary faces (emitted first, in source-face
 *        order) geometrically match the seed's faces.
 * @param seed Source mesh the operator ran on.
 * @param out Operator output.
 * @param corners_per_source Output corners expected near each seed corner:
 *        2 for truncate (both edge cuts of a corner), 1 for expand/snub/chamfer.
 * @param tol Max distance from an output corner to its seed corner.
 * @details Primary face fi must have seed_count(fi) * corners_per_source sides
 *          with exactly corners_per_source of them within tol of each seed
 *          corner — pinning emission order, side counts, and geometry at once.
 */
inline void check_primary_faces_match_seed(const PolyMesh &seed,
                                           const PolyMesh &out,
                                           int corners_per_source, float tol) {
  const size_t F = seed.face_counts.size();
  HS_EXPECT_GE(out.face_counts.size(), F);
  if (out.face_counts.size() < F)
    return;
  size_t seed_off = 0;
  size_t out_off = 0;
  for (size_t fi = 0; fi < F; ++fi) {
    const int bc = seed.face_counts[fi];
    const int oc = out.face_counts[fi];
    HS_EXPECT_EQ(oc, bc * corners_per_source);
    for (int k = 0; k < bc; ++k) {
      const math::Vector corner = seed.vertices[seed.faces[seed_off + k]];
      int near_count = 0;
      for (int j = 0; j < oc; ++j) {
        if ((out.vertices[out.faces[out_off + j]] - corner).length() <= tol)
          ++near_count;
      }
      HS_EXPECT_EQ(near_count, corners_per_source);
    }
    seed_off += bc;
    out_off += oc;
  }
}
