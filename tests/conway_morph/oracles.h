/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "tests/test_harness.h"
#include "core/mesh/mesh.h"
#include "tests/mesh_test_util.h"
#include <algorithm>
#include <numeric>
#include <optional>
#include <span>
#include <vector>
#include <cstring>
#include <cmath>
#include <limits>

namespace hs_test::conway_morph_tests {

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
  const float mean = edge_length_mean(m);
  float worst = 0.0f;
  for_each_face_edge(m, [&](const math::Vector &a, const math::Vector &b) {
    worst = fold_worst(worst, std::abs(math::distance_between(a, b) - mean));
  });
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
  if (!expect_same_mesh_sizes(got, src))
    return false;
  HS_EXPECT_EQ(std::memcmp(got.face_counts.data(), src.face_counts.data(),
                           got.face_counts.size() * sizeof(uint8_t)),
               0);
  HS_EXPECT_EQ(std::memcmp(got.faces.data(), src.faces.data(),
                           got.faces.size() * sizeof(uint16_t)),
               0);

  for (size_t i = 0; i < got.vertices.size(); ++i) {
    size_t nearest = src.vertices.size();
    float best = std::numeric_limits<float>::infinity();
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
  const size_t f = seed.face_counts.size();
  HS_EXPECT_GE(out.face_counts.size(), f);
  if (out.face_counts.size() < f)
    return;
  size_t seed_off = 0;
  size_t out_off = 0;
  for (size_t fi = 0; fi < f; ++fi) {
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

/** @brief Sorted face cycles after remapping each vertex index. */
inline std::optional<std::vector<std::vector<size_t>>>
face_cycles(const PolyMesh &mesh, std::span<const size_t> indices) {
  std::vector<std::vector<size_t>> faces;
  size_t offset = 0;
  for (const size_t count : mesh.face_counts) {
    if (count < 3 || count > mesh.faces.size() - offset)
      return std::nullopt;
    std::vector<size_t> cycle;
    for (size_t corner = 0; corner < count; ++corner) {
      const size_t index = mesh.faces[offset++];
      if (index >= indices.size())
        return std::nullopt;
      cycle.push_back(indices[index]);
    }
    auto canonical = cycle;
    for (size_t corner = 0; corner < count; ++corner) {
      canonical = std::min(canonical, cycle);
      const std::vector<size_t> reversed(cycle.rbegin(), cycle.rend());
      canonical = std::min(canonical, reversed);
      std::rotate(cycle.begin(), cycle.begin() + 1, cycle.end());
    }
    faces.push_back(std::move(canonical));
  }
  if (offset != mesh.faces.size())
    return std::nullopt;
  std::sort(faces.begin(), faces.end());
  return faces;
}

/** @brief Checks equivalent vertices and face cycles in either traversal direction. */
inline void check_equal_up_to_vertex_order(const PolyMesh &got,
                                           const PolyMesh &want, float tol) {
  if (!expect_same_mesh_sizes(got, want))
    return;
  std::vector<size_t> correspondence(got.vertices.size());
  std::vector<bool> used(want.vertices.size(), false);
  for (size_t i = 0; i < got.vertices.size(); ++i) {
    size_t j = 0;
    for (; j < want.vertices.size(); ++j)
      if (!used[j] && (got.vertices[i] - want.vertices[j]).length() <= tol)
        break;
    HS_CONTEXT("vertex", i);
    HS_EXPECT_LT(j, want.vertices.size());
    if (j == want.vertices.size())
      return;
    used[j] = true;
    correspondence[i] = j;
  }
  std::vector<size_t> identity(want.vertices.size());
  std::iota(identity.begin(), identity.end(), size_t{0});
  const auto got_faces = face_cycles(got, correspondence);
  const auto want_faces = face_cycles(want, identity);
  HS_EXPECT_TRUE(got_faces && want_faces && got_faces == want_faces);
}

} // namespace hs_test::conway_morph_tests
