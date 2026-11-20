/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file sphere_tiling.h @brief Radially projected solid faces and edges. */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <utility>
#include <vector>
#include "core/mesh/solids.h"
#include "core/pattern/sphere_lut.h"
#include "core/render/pullback/sphere_tiling_edges.h"
#include "core/render/sdf/face.h"

namespace Pullback::Source {

/** @brief Immutable cubemap index over a solid's spherical faces and edges. */
class SphereTilingGeometry {
public:
  static constexpr uint16_t NO_FACE = UINT16_MAX; ///< Uncovered direction.

  /// Cubemap candidate storage for the projected edge geometry.
  struct EdgeIndex {
    int steps;                     ///< Cells along each cube face axis.
    std::vector<uint32_t> offsets; ///< Candidate start offset per cell.
    std::vector<uint16_t> ids;     ///< Candidate edge IDs.
  };

  /** @brief Inward angular distances from a face boundary. */
  struct FaceDepth {
    float sample;      ///< Distance at the sampled direction, in radians.
    float center;      ///< Distance at the face center, in radians.
    float across_seam; ///< Center depth on the other side of a crossed seam.
  };

  /** @brief Builds the geometry and index for a registry solid.
   * @param registry_index Index into the solids registry. */
  explicit SphereTilingGeometry(size_t registry_index) {
    build_geometry(Solids::get_entry(registry_index));
    build_face_index();
    edge_index = build_edge_index(CELL_STEPS);
  }

  /** @brief Finds the nearest projected face containing a unit direction.
   * @param unit_direction Unit direction in the solid's frame.
   * @return Face ID or NO_FACE. */
  uint16_t face_at(const math::Vector &unit_direction) const {
    const auto candidates = face_index.candidates(
        tiling_edges_detail::cell_at(unit_direction, CELL_STEPS));
    uint16_t owner = NO_FACE;
    float nearest = 0.0f;
    for (uint16_t id : candidates) {
      const Face &face = faces[id];
      const float score = math::dot(face.selector, unit_direction);
      if (score > nearest && contains(face, unit_direction)) {
        owner = id;
        nearest = score;
      }
    }
    return owner;
  }

  /** @brief Measures inward depth within a projected face.
   * @param face Face ID returned by face_at(), excluding NO_FACE.
   * @param unit_direction Unit direction in that face.
   * @param filter_radius Angular radius across the nearest seam, in radians.
   * @return Sample depth and each covered face's center depth. */
  FaceDepth face_depth(uint16_t face, const math::Vector &unit_direction,
                       float filter_radius) const {
    const Face &polygon = faces[face];
    const auto [distance, seam] = boundary_at(polygon, unit_direction);
    float across = polygon.center_depth;
    if (distance < filter_radius) {
      const auto [first, second] = edge_faces[seam];
      const uint16_t neighbor = first == face ? second : first;
      if (neighbor != NO_FACE)
        across = faces[neighbor].center_depth;
    }
    return {distance, polygon.center_depth, across};
  }

  /** @brief Views the edge geometry with its default candidate index.
   * @return Indexed projected edges owned by this geometry. */
  SphereTilingEdgeView edge_view() const { return edge_view(edge_index); }

  /** @brief Views the edge geometry with an external candidate index.
   * @param index Index built by build_edge_index().
   * @return Indexed projected edges. */
  SphereTilingEdgeView edge_view(const EdgeIndex &index) const {
    return {vertices, edges, index.offsets, index.ids, index.steps};
  }

  /** @brief Builds the candidate index at the specified resolution.
   * @param steps Cells along one axis of each cube face.
   * @return Indexed edge candidates. */
  EdgeIndex build_edge_index(int steps) const {
    using namespace tiling_edges_detail;
    EdgeIndex index{steps, {}, {}};
    index.offsets.push_back(0);
    std::vector<float> edge_chords(edges.size());
    for_each_cell(steps, [&](int, const math::Vector &center, float radius) {
      float nearest = 4.0f;
      for (size_t id = 0; id < edges.size(); ++id) {
        edge_chords[id] = edge_chord_squared(edges[id], center, vertices);
        nearest = std::min(nearest, edge_chords[id]);
      }
      // Cell radius bounds the change in the distance to any edge.
      const float limit = chord_squared(
          std::min(math::PI_F, chord_angle(nearest) + 2.0f * radius));
      for (size_t id = 0; id < edges.size(); ++id)
        if (edge_chords[id] <= limit)
          index.ids.push_back(static_cast<uint16_t>(id));
      index.offsets.push_back(static_cast<uint32_t>(index.ids.size()));
    });
    return index;
  }

private:
  static constexpr int CELL_STEPS = 12;
  static constexpr int CELL_COUNT = 6 * CELL_STEPS * CELL_STEPS;
  static constexpr float CHART_EPSILON = 1.0e-6f;
  static constexpr float CAP_ANGLE_PAD = 1.0e-4f;

  struct Index {
    std::array<uint32_t, CELL_COUNT + 1> offsets{};
    std::vector<uint16_t> ids;

    std::span<const uint16_t> candidates(int cell) const {
      return std::span<const uint16_t>(ids).subspan(
          offsets[cell], offsets[cell + 1] - offsets[cell]);
    }
  };

  struct Point {
    float x;
    float y;
    uint16_t edge;
  };

  struct Face {
    math::Vector axis;
    math::Vector u;
    math::Vector v;
    math::Vector selector;
    float radius;
    float center_depth;
    size_t first;
    uint16_t count;
  };

  using EdgeIds = std::pair<uint16_t, uint16_t>;

  std::vector<math::Vector> vertices;
  std::vector<Face> faces;
  std::vector<Point> corners;
  std::vector<SphereTilingEdge> edges;
  std::vector<EdgeIds> edge_faces;
  Index face_index;
  EdgeIndex edge_index{};

  bool contains(const Face &face, const math::Vector &direction) const {
    if (face.radius < 0.0f || math::dot(face.axis, direction) <= CHART_EPSILON)
      return false;

    float px;
    float py;
    SDF::gnomonic_project(direction, face.axis, face.u, face.v, px, py);
    bool inside = false;
    Point a = corners[face.first + face.count - 1];
    for (size_t index = face.first; index < face.first + face.count; ++index) {
      const Point b = corners[index];
      const float dx = b.x - a.x;
      const float dy = b.y - a.y;
      const float cross = dx * (py - a.y) - dy * (px - a.x);
      if (fabsf(cross) <= CHART_EPSILON * (fabsf(dx) + fabsf(dy)) &&
          px >= std::min(a.x, b.x) - CHART_EPSILON &&
          px <= std::max(a.x, b.x) + CHART_EPSILON &&
          py >= std::min(a.y, b.y) - CHART_EPSILON &&
          py <= std::max(a.y, b.y) + CHART_EPSILON)
        return true;
      if ((a.y > py) != (b.y > py) && (cross > 0.0f) == (dy > 0.0f))
        inside = !inside;
      a = b;
    }
    return inside;
  }

  std::pair<float, uint16_t> boundary_at(const Face &face,
                                         const math::Vector &direction) const {
    using namespace tiling_edges_detail;
    uint16_t nearest = corners[face.first].edge;
    float best = edge_chord_squared(edges[nearest], direction, vertices);
    for (size_t i = face.first + 1; i < face.first + face.count; ++i) {
      const uint16_t edge = corners[i].edge;
      const float distance =
          edge_chord_squared(edges[edge], direction, vertices);
      if (distance < best) {
        best = distance;
        nearest = edge;
      }
    }
    return {chord_angle(best), nearest};
  }

  Face build_face(const PolyMesh &mesh, std::span<const uint16_t> ids,
                  std::span<const EdgeIds> edge_ids) {
    math::Vector center{};
    math::Vector axis{};
    math::Vector normal{};
    const math::Vector &origin = mesh.vertices[ids.front()];
    math::Vector previous = mesh.vertices[ids.back()] - origin;
    for (uint16_t id : ids) {
      const auto &v = mesh.vertices[id];
      center += v;
      axis += vertices[id];
      const math::Vector local = v - origin;
      normal += math::cross(previous, local);
      previous = local;
    }
    center = center / static_cast<float>(ids.size());
    axis = math::normalized_or(axis, vertices[ids.front()]);
    const float normal_length = normal.length();
    normal = normal_length > 0.0f ? normal / normal_length : axis;
    const auto min_depth = [&](const math::Vector &candidate) {
      float depth = 1.0f;
      for (uint16_t id : ids)
        depth = std::min(depth, math::dot(candidate, vertices[id]));
      return depth;
    };
    float depth = min_depth(axis);
    if (depth <= CHART_EPSILON) {
      const math::Vector aligned =
          math::dot(normal, center) < 0.0f ? -normal : normal;
      axis = aligned;
      depth = min_depth(axis);
    }

    const auto count = static_cast<uint16_t>(ids.size());
    Face face{axis, {}, {}, {}, -1.0f, 0.0f, corners.size(), count};
    if (depth <= CHART_EPSILON)
      return face;

    const float plane_offset = math::dot(normal, center);
    face.selector = fabsf(plane_offset) > CHART_EPSILON
                        ? normal / plane_offset
                        : axis / math::dot(axis, center);
    face.u = math::tangent_axis(axis);
    face.v = math::cross(axis, face.u);
    face.radius = acosf(std::min(1.0f, depth));
    const auto face_id = static_cast<uint16_t>(faces.size());
    uint16_t previous_id = ids.back();
    for (uint16_t id : ids) {
      Point point;
      SDF::gnomonic_project(vertices[id], axis, face.u, face.v, point.x,
                            point.y);
      const EdgeIds key{std::min(previous_id, id), std::max(previous_id, id)};
      const auto match =
          std::lower_bound(edge_ids.begin(), edge_ids.end(), key);
      point.edge = static_cast<uint16_t>(match - edge_ids.begin());
      auto &[first, second] = edge_faces[point.edge];
      if (first == NO_FACE)
        first = face_id;
      else if (first != face_id && second == NO_FACE)
        second = face_id;
      corners.push_back(point);
      previous_id = id;
    }
    const auto center_direction = math::normalized_or(center, axis);
    face.center_depth = boundary_at(face, center_direction).first;
    return face;
  }

  void build_geometry(const Solids::Entry &entry) {
    constexpr size_t ARENA_BYTES = 4 * 1024 * 1024;
    auto storage = std::make_unique_for_overwrite<uint8_t[]>(2 * ARENA_BYTES);
    Arena a(storage.get(), ARENA_BYTES);
    Arena b(storage.get() + ARENA_BYTES, ARENA_BYTES);
    const PolyMesh mesh = entry.generate(a, b);
    // no death case: the tiling test builds every solid in the fixed registry.
    HS_CHECK(mesh.num_faces() < NO_FACE, "sphere tiling: too many faces");
    vertices.reserve(mesh.vertices.size());
    faces.reserve(mesh.num_faces());
    corners.reserve(mesh.faces.size());
    for (const auto &vertex : mesh.vertices)
      vertices.push_back(vertex.normalized());

    std::vector<EdgeIds> edge_ids;
    edge_ids.reserve(mesh.faces.size());
    const std::span<const uint16_t> indices(mesh.faces.data(),
                                            mesh.faces.size());
    size_t offset = 0;
    for (uint8_t count : mesh.face_counts) {
      const auto ids = indices.subspan(offset, count);
      uint16_t previous = ids.back();
      for (uint16_t id : ids) {
        edge_ids.emplace_back(std::min(previous, id), std::max(previous, id));
        previous = id;
      }
      offset += count;
    }
    std::sort(edge_ids.begin(), edge_ids.end());
    edge_ids.erase(std::unique(edge_ids.begin(), edge_ids.end()),
                   edge_ids.end());
    // no death case: the tiling test indexes edges for every fixed registry solid.
    HS_CHECK(edge_ids.size() <= UINT16_MAX, "sphere tiling: too many edges");
    edges.reserve(edge_ids.size());
    edge_faces.resize(edge_ids.size(), {NO_FACE, NO_FACE});
    for (auto [first, last] : edge_ids) {
      const auto &from = vertices[first];
      const auto &to = vertices[last];
      const auto cross = math::cross(from, to);
      const float length = cross.length();
      if (length > 0.0f) {
        const auto normal = cross / length;
        const float arc = atan2f(length, math::dot(from, to));
        edges.push_back({first, last, normal, math::cross(normal, from),
                         math::cross(to, normal), 1.0f / arc});
      } else {
        edges.push_back({first, last, {}, {}, {}, 0.0f});
      }
    }

    offset = 0;
    for (uint8_t count : mesh.face_counts) {
      faces.push_back(
          build_face(mesh, indices.subspan(offset, count), edge_ids));
      offset += count;
    }
  }

  template <class Visitor> static void for_each_cell(int steps, Visitor visit) {
    const float step = 2.0f / steps;
    for (int cell = 0; cell < 6 * steps * steps; ++cell) {
      const int side = cell / (steps * steps);
      const float u = -1.0f + step * (cell % steps);
      const float v = -1.0f + step * ((cell / steps) % steps);
      const auto center = Pattern::sphere_face_direction(side, u + 0.5f * step,
                                                         v + 0.5f * step);
      float largest_chord = 0.0f;
      for (float x : {u, u + step})
        for (float y : {v, v + step})
          largest_chord =
              std::max(largest_chord,
                       math::distance_squared(
                           center, Pattern::sphere_face_direction(side, x, y)));
      const float radius =
          tiling_edges_detail::chord_angle(largest_chord) + CAP_ANGLE_PAD;
      visit(cell, center, radius);
    }
  }

  void build_face_index() {
    for_each_cell(CELL_STEPS,
                  [&](int cell, const math::Vector &center, float radius) {
                    for (size_t id = 0; id < faces.size(); ++id) {
                      const Face &face = faces[id];
                      if (face.radius >= 0.0f && math::dot(face.axis, center) >=
                                                     cosf(face.radius + radius))
                        face_index.ids.push_back(static_cast<uint16_t>(id));
                    }
                    face_index.offsets[cell + 1] =
                        static_cast<uint32_t>(face_index.ids.size());
                  });
  }
};

} // namespace Pullback::Source
