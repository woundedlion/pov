/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include "math/3dmath.h"
#include "render/ray/contract.h"

namespace SDF {

/** @brief Triangular-prism edges with equilateral XY cells and Z layers. */
struct TriangularFramework {
  struct PlaneFamily {
    math::Vector normal;
    float spacing;
  };
  static constexpr size_t STREAM_COUNT = 4;
  static constexpr float TRIANGLE_HEIGHT = 0.8660254037844386f;
  float cell_size = 1.0f;
  float layer_height = 1.0f;
  float wire_radius = 0.05f;
  math::Vector origin{};

  bool valid() const {
    return Raycast::finite(cell_size) && cell_size > 0.0f &&
           Raycast::finite(layer_height) && layer_height > 0.0f &&
           Raycast::finite(wire_radius) && wire_radius > 0.0f &&
           Raycast::finite(origin.x) && Raycast::finite(origin.y) &&
           Raycast::finite(origin.z);
  }

  std::array<PlaneFamily, STREAM_COUNT> plane_families() const {
    const float HEIGHT = TRIANGLE_HEIGHT * cell_size;
    return {{{{0.0f, 1.0f, 0.0f}, HEIGHT},
             {{TRIANGLE_HEIGHT, 0.5f, 0.0f}, HEIGHT},
             {{-TRIANGLE_HEIGHT, 0.5f, 0.0f}, HEIGHT},
             {{0.0f, 0.0f, 1.0f}, layer_height}}};
  }

  math::Vector edge_offset(const math::Vector &p, uint32_t &feature) const {
    const math::Vector Q = p - origin;
    const float HEIGHT = TRIANGLE_HEIGHT * cell_size;
    const float B = roundf(Q.y / HEIGHT);
    const float A = roundf(Q.x / cell_size - 0.5f * B);
    math::Vector result{cell_size, cell_size, 0.0f};
    float best = math::dot(result, result);
    feature = 3;
    for (int b = -1; b <= 1; ++b) {
      for (int a = -1; a <= 1; ++a) {
        const math::Vector OFFSET{Q.x - cell_size * (A + a + 0.5f * (B + b)),
                                  Q.y - HEIGHT * (B + b), 0.0f};
        const float SQUARED = math::dot(OFFSET, OFFSET);
        if (SQUARED < best) {
          best = SQUARED;
          result = OFFSET;
        }
      }
    }
    const float Z = Q.z - layer_height * roundf(Q.z / layer_height);
    const auto FAMILIES = plane_families();
    for (uint32_t i = 0; i < 3; ++i) {
      const float U = math::dot(Q, FAMILIES[i].normal);
      const float D = U - HEIGHT * roundf(U / HEIGHT);
      const float SQUARED = D * D + Z * Z;
      if (SQUARED < best) {
        best = SQUARED;
        result = FAMILIES[i].normal * D + math::Vector{0.0f, 0.0f, Z};
        feature = i;
      }
    }
    return result;
  }

  float distance(const math::Vector &p) const {
    uint32_t feature;
    return edge_offset(p, feature).magnitude() - wire_radius;
  }

  math::Vector normal(const math::Vector &p) const {
    uint32_t feature;
    const math::Vector OFFSET = edge_offset(p, feature);
    const float LENGTH = OFFSET.magnitude();
    return LENGTH > 0.0f ? OFFSET * (1.0f / LENGTH) : math::Vector{};
  }

  Raycast::QueryCapabilities capabilities() const {
    return {true, true, true, 0.0f};
  }

  Raycast::QuerySample sample(const math::Vector &p) const {
    uint32_t feature;
    const float VALUE = edge_offset(p, feature).magnitude() - wire_radius;
    return {VALUE, fabsf(VALUE), VALUE == 0.0f, 0, feature};
  }
};

/** @brief Four plane streams producing approximate framework coverage layers. */
struct FrameworkEvents {
  static constexpr size_t STREAM_COUNT = TriangularFramework::STREAM_COUNT;
  struct Cursor {
    float next = 0.0f;
    float step = 0.0f;
    bool active = false;
  };
  TriangularFramework geometry;
  Raycast::Ray ray;
  Raycast::Footprint footprint;
  std::array<Cursor, STREAM_COUNT> cursors{};

  FrameworkEvents(const TriangularFramework &geometry, const Raycast::Ray &ray,
                  Raycast::Footprint footprint = {})
      : geometry(geometry), ray(ray), footprint(footprint) {
    if (!geometry.valid() || !ray.valid())
      return;
    const auto FAMILIES = geometry.plane_families();
    const math::Vector START = ray.at(ray.interval.near) - geometry.origin;
    for (size_t i = 0; i < STREAM_COUNT; ++i) {
      const float SPEED = math::dot(ray.direction, FAMILIES[i].normal);
      if (SPEED == 0.0f)
        continue;
      const float POSITION = math::dot(START, FAMILIES[i].normal);
      const float CELL = POSITION / FAMILIES[i].spacing;
      const float PLANE = SPEED > 0.0f ? ceilf(CELL) : floorf(CELL);
      Cursor &cursor = cursors[i];
      cursor.next =
          ray.interval.near + (PLANE * FAMILIES[i].spacing - POSITION) / SPEED;
      cursor.step = FAMILIES[i].spacing / fabsf(SPEED);
      cursor.active = Raycast::finite(cursor.next) &&
                      Raycast::finite(cursor.step) && cursor.step > 0.0f;
    }
  }

  bool active(size_t index) const { return cursors[index].active; }
  float distance(size_t index) const { return cursors[index].next; }

  Raycast::Contribution candidate(size_t index) const {
    Raycast::Contribution result;
    result.t = distance(index);
    const math::Vector P = ray.at(result.t);
    const auto SAMPLE = geometry.sample(P);
    const float WIDTH = footprint.at(result.t);
    result.coverage = WIDTH > 0.0f
                          ? std::clamp(0.5f - SAMPLE.field / WIDTH, 0.0f, 1.0f)
                          : (SAMPLE.field <= 0.0f ? 1.0f : 0.0f);
    result.feature = SAMPLE.feature;
    // Coincident plane reports define one approximate junction layer.
    result.merge_identity = 0;
    return result;
  }

  void advance(size_t index) {
    Cursor &cursor = cursors[index];
    const float NEXT = cursor.next + cursor.step;
    cursor.active = Raycast::finite(NEXT) && NEXT > cursor.next;
    cursor.next = NEXT;
  }
};

} // namespace SDF
