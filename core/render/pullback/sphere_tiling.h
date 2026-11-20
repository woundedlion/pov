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
#include "core/render/sdf/face.h"

namespace Pullback::Source {

/** @brief Immutable cubemap index over a solid's spherical faces and edges. */
class SphereTilingGeometry {
public:
  static constexpr uint16_t NO_FACE = UINT16_MAX; ///< Uncovered direction.

  /** @brief Nearest seam sample. */
  struct EdgeSample {
    float distance; ///< Angular distance in radians.
    float position; ///< Arc fraction from the lower vertex index.
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
    build_index();
  }

  /** @brief Finds the nearest projected face containing a unit direction.
   * @param unit_direction Unit direction in the solid's frame.
   * @return Face ID or NO_FACE. */
  uint16_t face_at(const math::Vector &unit_direction) const {
    const auto candidates = face_index.candidates(cell_at(unit_direction));
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
      const auto &edge = edges[seam];
      const uint16_t neighbor =
          edge.first_face == face ? edge.second_face : edge.first_face;
      if (neighbor != NO_FACE)
        across = faces[neighbor].center_depth;
    }
    return {distance, polygon.center_depth, across};
  }

  /** @brief Samples the nearest finite seam arc within the distance limit.
   * @param unit_direction Unit direction in the solid's frame.
   * @param max_distance_radians Maximum reported distance in [0, pi].
   * @return Distance and arc fraction, or the distance limit and zero. */
  EdgeSample edge_at(const math::Vector &unit_direction,
                     float max_distance_radians) const {
    float best = chord_squared(max_distance_radians);
    const Edge *nearest = nullptr;
    for (uint16_t id : edge_index.candidates(cell_at(unit_direction))) {
      const float distance = edge_chord_squared(edges[id], unit_direction);
      if (distance < best) {
        best = distance;
        nearest = &edges[id];
      }
    }
    if (!nearest)
      return {max_distance_radians, 0.0f};
    const float distance = std::min(max_distance_radians, chord_angle(best));
    return {distance, edge_position(*nearest, unit_direction)};
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

  struct Edge {
    uint16_t a;
    uint16_t b;
    math::Vector normal;
    math::Vector from_a;
    math::Vector from_b;
    float inv_length;
    uint16_t first_face = NO_FACE;
    uint16_t second_face = NO_FACE;
  };

  using EdgeIds = std::pair<uint16_t, uint16_t>;

  std::vector<math::Vector> vertices;
  std::vector<Face> faces;
  std::vector<Point> corners;
  std::vector<Edge> edges;
  Index face_index;
  Index edge_index;

  static float chord_squared(float angle) {
    const float half = sinf(0.5f * angle);
    return 4.0f * half * half;
  }

  static float chord_angle(float distance_squared) {
    return 2.0f * asinf(0.5f * sqrtf(std::min(4.0f, distance_squared)));
  }

  static int cell_at(const math::Vector &direction) {
    const auto [face, u, v, major] = math::project_cube_face(direction);
    const float scale = (0.5f * CELL_STEPS) / major;
    const auto x = hs::clamp(static_cast<int>(u * scale + 0.5f * CELL_STEPS), 0,
                             CELL_STEPS - 1);
    const auto y = hs::clamp(static_cast<int>(v * scale + 0.5f * CELL_STEPS), 0,
                             CELL_STEPS - 1);
    return (face * CELL_STEPS + y) * CELL_STEPS + x;
  }

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

  float edge_chord_squared(const Edge &edge,
                           const math::Vector &direction) const {
    if (edge.inv_length > 0.0f && math::dot(edge.from_a, direction) >= 0.0f &&
        math::dot(edge.from_b, direction) >= 0.0f) {
      const float sine = math::dot(edge.normal, direction);
      const float sine_sq = std::min(1.0f, sine * sine);
      return 2.0f * sine_sq / (1.0f + sqrtf(1.0f - sine_sq));
    }
    return std::min(math::distance_squared(direction, vertices[edge.a]),
                    math::distance_squared(direction, vertices[edge.b]));
  }

  float edge_position(const Edge &edge, const math::Vector &direction) const {
    if (edge.inv_length == 0.0f)
      return 0.0f;
    const float from_start = math::dot(edge.from_a, direction);
    if (from_start > 0.0f && math::dot(edge.from_b, direction) > 0.0f) {
      const float angle =
          atan2f(from_start, math::dot(vertices[edge.a], direction));
      return hs::clamp(angle * edge.inv_length, 0.0f, 1.0f);
    }
    return math::dot(vertices[edge.a], direction) >=
                   math::dot(vertices[edge.b], direction)
               ? 0.0f
               : 1.0f;
  }

  std::pair<float, uint16_t> boundary_at(const Face &face,
                                         const math::Vector &direction) const {
    uint16_t nearest = corners[face.first].edge;
    float best = edge_chord_squared(edges[nearest], direction);
    for (size_t i = face.first + 1; i < face.first + face.count; ++i) {
      const uint16_t edge = corners[i].edge;
      const float distance = edge_chord_squared(edges[edge], direction);
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
      auto &edge = edges[point.edge];
      if (edge.first_face == NO_FACE)
        edge.first_face = face_id;
      else if (edge.first_face != face_id && edge.second_face == NO_FACE)
        edge.second_face = face_id;
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

  void build_index() {
    constexpr float STEP = 2.0f / CELL_STEPS;
    std::vector<float> edge_chords(edges.size());
    for (int cell = 0; cell < CELL_COUNT; ++cell) {
      const int side = cell / (CELL_STEPS * CELL_STEPS);
      const float u = -1.0f + STEP * (cell % CELL_STEPS);
      const float v = -1.0f + STEP * ((cell / CELL_STEPS) % CELL_STEPS);
      const auto center = Pattern::sphere_face_direction(side, u + 0.5f * STEP,
                                                         v + 0.5f * STEP);
      float largest_chord = 0.0f;
      for (float x : {u, u + STEP})
        for (float y : {v, v + STEP})
          largest_chord =
              std::max(largest_chord,
                       math::distance_squared(
                           center, Pattern::sphere_face_direction(side, x, y)));
      const float radius = chord_angle(largest_chord) + CAP_ANGLE_PAD;

      for (size_t id = 0; id < faces.size(); ++id) {
        const Face &face = faces[id];
        if (face.radius >= 0.0f &&
            math::dot(face.axis, center) >= cosf(face.radius + radius))
          face_index.ids.push_back(static_cast<uint16_t>(id));
      }
      face_index.offsets[cell + 1] =
          static_cast<uint32_t>(face_index.ids.size());

      float nearest = 4.0f;
      for (size_t id = 0; id < edges.size(); ++id) {
        edge_chords[id] = edge_chord_squared(edges[id], center);
        nearest = std::min(nearest, edge_chords[id]);
      }
      // Cell radius bounds the change in the distance to any edge.
      const float limit = chord_squared(
          std::min(math::PI_F, chord_angle(nearest) + 2.0f * radius));
      for (size_t id = 0; id < edges.size(); ++id)
        if (edge_chords[id] <= limit)
          edge_index.ids.push_back(static_cast<uint16_t>(id));
      edge_index.offsets[cell + 1] =
          static_cast<uint32_t>(edge_index.ids.size());
    }
  }
};

} // namespace Pullback::Source
