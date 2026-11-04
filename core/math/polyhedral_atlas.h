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
  float offset = 0.0f; ///< Distance from the side at the origin.
};

/**
 * @brief A convex solid's faces as separate tiles, as
 *        polyhedral_atlas_projection() reads them.
 * @tparam FACES Face count.
 * @tparam MAX_SIDES Side count of the largest face.
 * @details Coordinates are scaled so the solid's first edge has unit length.
 */
template <size_t FACES, size_t MAX_SIDES> struct PolyhedralAtlas {
  static_assert(FACES * MAX_SIDES <= 256, "edge classes are one byte");
  /** Face normal over the face plane's distance from the origin. */
  std::array<math::Vector, FACES> selectors{};
  /** Atlas x coordinate, per face, of a point on the face plane. */
  std::array<math::Vector, FACES> x_axes{};
  /** Atlas y coordinate, per face, of a point on the face plane. */
  std::array<math::Vector, FACES> y_axes{};
  /** Atlas position of each face plane's closest point to the origin. */
  std::array<math::Complex, FACES> centers{};
  /** Reciprocal of each face plane's distance from the origin, in atlas
   *  units. */
  std::array<float, FACES> inverse_heights{};
  /** Each tile's sides; side `k` runs from vertex `k` to vertex `k + 1`. */
  std::array<std::array<AtlasSide, MAX_SIDES>, FACES> sides{};
  /** Side count of each face. */
  std::array<uint8_t, FACES> side_counts{};
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

/** @brief Monotonic stand-in for an azimuth, in [0, 4).
 *  @param x Cosine-side component; not both zero with @p z.
 *  @param z Sine-side component.
 *  @return A value increasing with atan2(z, x) taken in [0, 2 pi). */
constexpr float pseudo_azimuth(float x, float z) {
  const float t = x / ((x < 0.0f ? -x : x) + (z < 0.0f ? -z : z));
  return z >= 0.0f ? 1.0f - t : 3.0f + t;
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
  PolyhedralAtlas<F, K> atlas{};

  std::array<size_t, F> first{};
  for (size_t f = 1; f < F; ++f)
    first[f] = first[f - 1] + Solid::face_counts[f - 1];
  const auto index = [&](size_t f, size_t k) {
    return static_cast<size_t>(
        Solid::faces[first[f] + k % Solid::face_counts[f]]);
  };
  const auto normal_of = [&](const auto &vertices, size_t f) {
    math::Vector sum(0.0f, 0.0f, 0.0f);
    for (size_t k = 0; k < Solid::face_counts[f]; ++k)
      sum = sum + math::cross(vertices[index(f, k)], vertices[index(f, k + 1)]);
    return unit(sum);
  };

  size_t pole = 0;
  for (size_t f = 1; f < F; ++f)
    if (normal_of(Solid::vertices, f).y > normal_of(Solid::vertices, pole).y)
      pole = f;
  const math::Vector up = normal_of(Solid::vertices, pole);
  // Rodrigues rotation of `up` onto +Y; `up.y > 0` for a solid around the origin.
  const math::Vector axis = math::cross(up, math::Vector(0.0f, 1.0f, 0.0f));
  const float tilt = 1.0f / (1.0f + up.y);
  const auto tilted = [&](const math::Vector &v) {
    const math::Vector k_v = math::cross(axis, v);
    return v + k_v + math::cross(axis, k_v) * tilt;
  };
  const math::Vector anchor = tilted(Solid::vertices[index(pole, 0)]);
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
  std::array<float, F> heights{};
  float circumradius = 0.0f;
  for (size_t f = 0; f < F; ++f) {
    normals[f] = normal_of(vertices, f);
    const float distance = math::dot(normals[f], vertices[index(f, 0)]);
    heights[f] = distance;
    atlas.selectors[f] = normals[f] / distance;
    e1[f] = unit(vertices[index(f, 0)] - normals[f] * distance);
    e2[f] = math::cross(normals[f], e1[f]);
    atlas.side_counts[f] = Solid::face_counts[f];
    order[f] = f;
    for (size_t k = 0; k < Solid::face_counts[f]; ++k) {
      const math::Vector &v = vertices[index(f, k)];
      const float x = math::dot(v, e1[f]), y = math::dot(v, e2[f]);
      const float r = math::constexpr_sqrt(x * x + y * y);
      circumradius = r > circumradius ? r : circumradius;
    }
  }
  const auto before = [&](size_t a, size_t b) {
    const float dy = normals[a].y - normals[b].y;
    if (dy > 1e-3f || dy < -1e-3f)
      return dy > 0.0f;
    return pseudo_azimuth(normals[a].x, normals[a].z) <
           pseudo_azimuth(normals[b].x, normals[b].z);
  };
  for (size_t i = 1; i < F; ++i)
    for (size_t j = i; j > 0 && before(order[j], order[j - 1]); --j) {
      const size_t swap = order[j];
      order[j] = order[j - 1];
      order[j - 1] = swap;
    }

  const math::Vector edge = vertices[index(0, 1)] - vertices[index(0, 0)];
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
    atlas.inverse_heights[f] = 1.0f / (heights[f] * scale);
    for (size_t k = 0; k < Solid::face_counts[f]; ++k) {
      const math::Vector &a = vertices[index(f, k)];
      const math::Vector &b = vertices[index(f, k + 1)];
      const float ax = math::dot(a, atlas.x_axes[f]) + center.re;
      const float ay = math::dot(a, atlas.y_axes[f]) + center.im;
      const float dx = math::dot(b, atlas.x_axes[f]) + center.re - ax;
      const float dy = math::dot(b, atlas.y_axes[f]) + center.im - ay;
      const float side_length = math::constexpr_sqrt(dx * dx + dy * dy);
      const float mx = -dy / side_length, my = dx / side_length;
      atlas.sides[f][k] = {mx, my, -(mx * ax + my * ay)};
      const float r = math::constexpr_sqrt(ax * ax + ay * ay);
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
template <size_t FACES, size_t MAX_SIDES>
HS_O3_FN inline ProjectionKernelResult
polyhedral_atlas_projection(const PolyhedralAtlas<FACES, MAX_SIDES> &atlas,
                            const math::Vector &v, float meridian_cos,
                            float meridian_sin, float area_correction) {
  const math::Vector p(v.x * meridian_cos + v.z * meridian_sin, v.y,
                       v.z * meridian_cos - v.x * meridian_sin);
  uint8_t face = 0;
  float nearest = math::dot(p, atlas.selectors[0]);
  for (uint8_t candidate = 1; candidate < FACES; ++candidate) {
    const float d = math::dot(p, atlas.selectors[candidate]);
    if (d > nearest) {
      nearest = d;
      face = candidate;
    }
  }
  const float inverse = 1.0f / nearest;
  const math::Complex &center = atlas.centers[face];
  const float u = math::dot(p, atlas.x_axes[face]) * inverse;
  const float w = math::dot(p, atlas.y_axes[face]) * inverse;
  const float x = u + center.re, y = w + center.im;
  uint8_t side = 0;
  float distance = NO_EDGE_DISTANCE;
  for (uint8_t k = 0; k < atlas.side_counts[face]; ++k) {
    const AtlasSide &line = atlas.sides[face][k];
    const float d = line.x * x + line.y * y + line.offset;
    if (d < distance) {
      distance = d;
      side = k;
    }
  }
  float shrink = 1.0f;
  if (area_correction != 0.0f) {
    const float reciprocal_height = atlas.inverse_heights[face];
    const float secant =
        sqrtf(1.0f + (u * u + w * w) * reciprocal_height * reciprocal_height);
    shrink +=
        area_correction * (sqrtf(2.0f / (secant * (secant + 1.0f))) - 1.0f);
  }
  return {.coords =
              math::Complex(u * shrink + center.re, w * shrink + center.im),
          .region_id = face,
          .component_id = 0,
          .boundary_flags = projection_boundary(ProjectionBoundary::CUT),
          .fade_edge_distance = fmaxf(0.0f, distance) * shrink,
          .flags = 0,
          .traits = projection_traits(ProjectionTrait::CUT),
          .edge_class = static_cast<uint8_t>(face * MAX_SIDES + side)};
}

} // namespace projections
