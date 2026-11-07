/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file polyhedral_atlas.h
 * @brief Central projection of the sphere onto a convex solid whose faces are
 *        laid out as separate tiles in a plane.
 */

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "core/math/3dmath.h"
#include "core/math/projections.h"

namespace projections {

/** @brief One side of an atlas tile as a signed-distance line. */
struct AtlasSide {
  float x = 0.0f;      ///< Inward unit normal, x.
  float y = 0.0f;      ///< Inward unit normal, y.
  float offset = 0.0f; ///< Distance from the side at the tile center.
};

/**
 * @brief A convex solid's faces as separate tiles, as
 *        polyhedral_atlas_projection() reads them.
 * @tparam FACES Face count.
 * @tparam MAX_SIDES Side count of the largest face.
 * @tparam PAIRS Antipodal face pairs: `FACES / 2` for a centrally symmetric
 *         solid, otherwise 0.
 * @details Coordinates are scaled so the solid's first edge has unit length.
 */
template <size_t FACES, size_t MAX_SIDES, size_t PAIRS = 0>
struct PolyhedralAtlas {
  static_assert(FACES * MAX_SIDES <= 256, "edge classes are one byte");
  static_assert(PAIRS == 0 || 2 * PAIRS == FACES);
  /** Face normal over the face plane's distance from the origin. */
  std::array<math::Vector, FACES> selectors{};
  /** Atlas x coordinate, per face, of a point on the face plane. */
  std::array<math::Vector, FACES> x_axes{};
  /** Atlas y coordinate, per face, of a point on the face plane. */
  std::array<math::Vector, FACES> y_axes{};
  /** Atlas position of each face plane's closest point to the origin. */
  std::array<math::Complex, FACES> centers{};
  /** Each face plane's distance from the origin, with unit vertices. */
  std::array<float, FACES> heights{};
  /** Each tile's sides; side `k` runs from vertex `k` to vertex `k + 1`.
   *  Sides past a face's count never bind. */
  std::array<std::array<AtlasSide, MAX_SIDES>, FACES> sides{};
  /** Side count of each face. */
  std::array<uint8_t, FACES> side_counts{};
  /** Each antipodal pair as {face, opposite face}. */
  std::array<std::array<uint8_t, 2>, PAIRS> antipodes{};
  /** Largest distance of a tile vertex from the origin. */
  float radius = 0.0f;
};

namespace polyhedral_atlas_detail {

/** @brief Unit vector at compile time. @param v Nonzero vector.
 *  @return v / |v|. */
constexpr math::Vector unit(const math::Vector &v) {
  return v / math::constexpr_sqrt(math::dot(v, v));
}

/** @brief Side count of a solid's largest face. @tparam Solid Solid table.
 *  @return The largest `face_counts` entry. */
template <typename Solid> constexpr size_t max_sides() {
  size_t sides = 0;
  for (uint8_t count : Solid::face_counts)
    sides = count > sides ? count : sides;
  return sides;
}

/** @brief Vertex index of a face corner. @tparam Solid Solid table.
 *  @param f Face. @param k Corner, taken modulo the face's side count.
 *  @return Index into `Solid::vertices`. */
template <typename Solid> constexpr size_t corner(size_t f, size_t k) {
  size_t first = 0;
  for (size_t g = 0; g < f; ++g)
    first += Solid::face_counts[g];
  return static_cast<size_t>(Solid::faces[first + k % Solid::face_counts[f]]);
}

/** @brief Outward unit normal of a face. @tparam Solid Solid table.
 *  @param vertices The solid's vertices in any orientation. @param f Face.
 *  @return The normal of counter-clockwise face `f`. */
template <typename Solid, typename Vertices>
constexpr math::Vector face_normal(const Vertices &vertices, size_t f) {
  math::Vector sum(0.0f, 0.0f, 0.0f);
  for (size_t k = 0; k < Solid::face_counts[f]; ++k)
    sum = sum + math::cross(vertices[corner<Solid>(f, k)],
                            vertices[corner<Solid>(f, k + 1)]);
  return unit(sum);
}

/** @brief Equal-height antipodal face pairs of a solid. @tparam Solid Solid table.
 *  @return Half the face count when each face has an equal-height opposite, else 0. */
template <typename Solid> constexpr size_t antipodal_pairs() {
  constexpr size_t F = Solid::face_counts.size();
  for (size_t f = 0; f < F; ++f) {
    bool opposed = false;
    const auto normal = face_normal<Solid>(Solid::vertices, f);
    const float height =
        math::dot(normal, Solid::vertices[corner<Solid>(f, 0)]);
    for (size_t g = 0; g < F; ++g) {
      const auto opposite = face_normal<Solid>(Solid::vertices, g);
      const float difference =
          math::dot(opposite, Solid::vertices[corner<Solid>(g, 0)]) - height;
      opposed = opposed || (math::dot(normal, opposite) < -0.9999f &&
                            difference > -1e-4f && difference < 1e-4f);
    }
    if (!opposed)
      return 0;
  }
  return F / 2;
}

} // namespace polyhedral_atlas_detail

/**
 * @brief Lays out a convex solid's faces as tiles at compile time.
 * @tparam Solid Table with unit `vertices`, `face_counts` and counter-clockwise
 *         `faces`, enclosing the origin.
 * @return The atlas polyhedral_atlas_projection() reads.
 * @details The face whose normal points most nearly up is turned to the north
 * pole with its first vertex on the prime meridian. Tiles fill two rows, upper
 * first, in order of descending latitude then ascending longitude, spaced by
 * twice the largest face circumradius around a centered grid. Each tile turns
 * north up; a pole tile turns its prime-meridian side down in the north and
 * up in the south.
 */
template <typename Solid> constexpr auto make_polyhedral_atlas() {
  using namespace polyhedral_atlas_detail;
  constexpr size_t F = Solid::face_counts.size();
  constexpr size_t K = max_sides<Solid>();
  constexpr size_t V = Solid::vertices.size();
  constexpr size_t PAIRS = antipodal_pairs<Solid>();
  PolyhedralAtlas<F, K, PAIRS> atlas{};

  size_t pole = 0;
  for (size_t f = 1; f < F; ++f)
    if (face_normal<Solid>(Solid::vertices, f).y >
        face_normal<Solid>(Solid::vertices, pole).y)
      pole = f;
  const math::Vector up = face_normal<Solid>(Solid::vertices, pole);
  // Rodrigues rotation of `up` onto +Y; `up.y > 0` for a solid around the origin.
  const math::Vector axis = math::cross(up, math::Vector(0.0f, 1.0f, 0.0f));
  const float tilt = 1.0f / (1.0f + up.y);
  const auto tilted = [&](const math::Vector &v) {
    const math::Vector k_v = math::cross(axis, v);
    return v + k_v + math::cross(axis, k_v) * tilt;
  };
  const math::Vector anchor = tilted(Solid::vertices[corner<Solid>(pole, 0)]);
  const float anchor_radius =
      math::constexpr_sqrt(anchor.x * anchor.x + anchor.z * anchor.z);
  const float spin_cos = anchor.x / anchor_radius;
  const float spin_sin = anchor.z / anchor_radius;
  std::array<math::Vector, V> vertices{};
  for (size_t i = 0; i < V; ++i) {
    const math::Vector t = tilted(Solid::vertices[i]);
    vertices[i] = math::Vector(t.x * spin_cos + t.z * spin_sin, t.y,
                               t.z * spin_cos - t.x * spin_sin);
  }

  std::array<math::Vector, F> normals{};
  std::array<math::Vector, F> e1{};
  std::array<math::Vector, F> e2{};
  std::array<size_t, F> order{};
  float circumradius = 0.0f;
  for (size_t f = 0; f < F; ++f) {
    normals[f] = face_normal<Solid>(vertices, f);
    const float distance = math::dot(normals[f], vertices[corner<Solid>(f, 0)]);
    atlas.heights[f] = distance;
    atlas.selectors[f] = normals[f] / distance;
    e1[f] = unit(vertices[corner<Solid>(f, 0)] - normals[f] * distance);
    e2[f] = math::cross(normals[f], e1[f]);
    atlas.side_counts[f] = Solid::face_counts[f];
    order[f] = f;
    for (size_t k = 0; k < Solid::face_counts[f]; ++k) {
      const math::Vector &v = vertices[corner<Solid>(f, k)];
      const float x = math::dot(v, e1[f]), y = math::dot(v, e2[f]);
      const float r = math::constexpr_sqrt(x * x + y * y);
      circumradius = r > circumradius ? r : circumradius;
    }
  }
  const auto before = [&](size_t a, size_t b) {
    const float dy = normals[a].y - normals[b].y;
    if (dy > 1e-3f || dy < -1e-3f)
      return dy > 0.0f;
    return math::diamond_angle<false>(normals[a].z, normals[a].x) <
           math::diamond_angle<false>(normals[b].z, normals[b].x);
  };
  size_t pair = 0;
  for (size_t f = 0; f < F && PAIRS != 0; ++f)
    for (size_t g = f + 1; g < F; ++g)
      if (math::dot(normals[f], normals[g]) < -0.9999f)
        atlas.antipodes[pair++] = {static_cast<uint8_t>(f),
                                   static_cast<uint8_t>(g)};
  for (size_t i = 1; i < F; ++i)
    for (size_t j = i; j > 0 && before(order[j], order[j - 1]); --j) {
      const size_t swap = order[j];
      order[j] = order[j - 1];
      order[j - 1] = swap;
    }

  const math::Vector edge =
      vertices[corner<Solid>(0, 1)] - vertices[corner<Solid>(0, 0)];
  const float scale = 1.0f / math::constexpr_sqrt(math::dot(edge, edge));
  constexpr size_t COLUMNS = (F + 1) / 2;
  for (size_t i = 0; i < F; ++i) {
    const size_t f = order[i];
    math::Vector north =
        math::Vector(0.0f, 1.0f, 0.0f) - normals[f] * normals[f].y;
    if (math::dot(north, north) < 1e-6f)
      north = math::Vector(normals[f].y > 0.0f ? -1.0f : 1.0f, 0.0f, 0.0f);
    const float nx = math::dot(north, e1[f]), ny = math::dot(north, e2[f]);
    const float length = math::constexpr_sqrt(nx * nx + ny * ny);
    const float c = ny / length, s = nx / length;
    const float cell = 2.0f * circumradius * scale;
    const math::Complex center(cell * (static_cast<float>(i % COLUMNS) -
                                       0.5f * static_cast<float>(COLUMNS - 1)),
                               cell * (0.5f - static_cast<float>(i / COLUMNS)));
    atlas.x_axes[f] = (e1[f] * c - e2[f] * s) * scale;
    atlas.y_axes[f] = (e1[f] * s + e2[f] * c) * scale;
    atlas.centers[f] = center;
    for (size_t k = 0; k < K; ++k) {
      if (k >= Solid::face_counts[f]) {
        atlas.sides[f][k] = {0.0f, 0.0f, NO_EDGE_DISTANCE};
        continue;
      }
      const math::Vector &a = vertices[corner<Solid>(f, k)];
      const math::Vector &b = vertices[corner<Solid>(f, k + 1)];
      const float ax = math::dot(a, atlas.x_axes[f]);
      const float ay = math::dot(a, atlas.y_axes[f]);
      const float dx = math::dot(b, atlas.x_axes[f]) - ax;
      const float dy = math::dot(b, atlas.y_axes[f]) - ay;
      const float side_length = math::constexpr_sqrt(dx * dx + dy * dy);
      const float mx = -dy / side_length, my = dx / side_length;
      atlas.sides[f][k] = {mx, my, -(mx * ax + my * ay)};
      const float rx = ax + center.re, ry = ay + center.im;
      const float r = math::constexpr_sqrt(rx * rx + ry * ry);
      atlas.radius = r > atlas.radius ? r : atlas.radius;
    }
  }
  return atlas;
}

/** @brief The compile-time atlas of a solid. @tparam Solid Solid table. */
template <typename Solid>
inline constexpr auto POLYHEDRAL_ATLAS = make_polyhedral_atlas<Solid>();

/**
 * @brief Projects a direction onto a solid's atlas tiles.
 * @tparam FACES Face count.
 * @tparam MAX_SIDES Side count of the largest face.
 * @tparam PAIRS Antipodal face pairs.
 * @param atlas Atlas from make_polyhedral_atlas().
 * @param v Unit direction on the sphere.
 * @param meridian_cos Cosine of the central meridian.
 * @param meridian_sin Sine of the central meridian.
 * @param area_correction Blend in [0, 1] from central projection toward the
 *        equal-area azimuthal projection about each face center.
 * @return Atlas coordinates, with `region_id` set to the face, the distance to
 *         its nearest side, and `edge_class` `face * MAX_SIDES + side`.
 * @details Each direction is centrally projected onto the plane of the face
 * its ray crosses, so great circles are straight within an uncorrected tile.
 * Correction pulls points toward the tile center, matching the central scale
 * there; at 1 the tile is the face's Lambert azimuthal equal-area image, with
 * curved sides. The side distance shrinks with its point. Every side is a cut.
 */
template <size_t FACES, size_t MAX_SIDES, size_t PAIRS>
HS_O3_FN inline ProjectionKernelResult polyhedral_atlas_projection(
    const PolyhedralAtlas<FACES, MAX_SIDES, PAIRS> &atlas,
    const math::Vector &v, float meridian_cos, float meridian_sin,
    float area_correction) {
  const math::Vector p(v.x * meridian_cos + v.z * meridian_sin, v.y,
                       v.z * meridian_cos - v.x * meridian_sin);
  uint8_t face = 0;
  float nearest;
  if constexpr (PAIRS != 0) {
    const float first = math::dot(p, atlas.selectors[atlas.antipodes[0][0]]);
    nearest = fabsf(first);
    face = atlas.antipodes[0][first < 0.0f];
    for (size_t pair = 1; pair < PAIRS; ++pair) {
      const float d = math::dot(p, atlas.selectors[atlas.antipodes[pair][0]]);
      if (fabsf(d) > nearest) {
        nearest = fabsf(d);
        face = atlas.antipodes[pair][d < 0.0f];
      }
    }
  } else {
    nearest = math::dot(p, atlas.selectors[0]);
    for (uint8_t candidate = 1; candidate < FACES; ++candidate) {
      const float d = math::dot(p, atlas.selectors[candidate]);
      if (d > nearest) {
        nearest = d;
        face = candidate;
      }
    }
  }
  // Equal-area radius over central radius is cosine * sqrt(2 / (1 + cosine)).
  const float cosine = nearest * atlas.heights[face];
  const float rise = 1.0f + cosine;
  float blend = rise;
  if (area_correction != 0.0f)
    blend += area_correction * (cosine * sqrtf(2.0f * rise) - rise);
  const float reciprocal = 1.0f / (nearest * rise);
  const float inverse = rise * reciprocal;
  const float u = math::dot(p, atlas.x_axes[face]);
  const float w = math::dot(p, atlas.y_axes[face]);
  const float x = u * inverse, y = w * inverse;
  uint8_t side = 0;
  float distance = NO_EDGE_DISTANCE;
  for (size_t k = 0; k < MAX_SIDES; ++k) {
    const AtlasSide &line = atlas.sides[face][k];
    const float d = line.x * x + line.y * y + line.offset;
    if (d < distance) {
      distance = d;
      side = static_cast<uint8_t>(k);
    }
  }
  const float scale = blend * reciprocal;
  const float shrink = scale * nearest;
  const math::Complex &center = atlas.centers[face];
  return {.coords = math::Complex(u * scale + center.re, w * scale + center.im),
          .region_id = face,
          .component_id = 0,
          .boundary_flags = projection_boundary(ProjectionBoundary::CUT),
          .fade_edge_distance = fmaxf(0.0f, distance) * shrink,
          .flags = 0,
          .traits = projection_traits(ProjectionTrait::CUT),
          .edge_class = static_cast<uint8_t>(face * MAX_SIDES + side)};
}

} // namespace projections
