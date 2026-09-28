/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include "math/3dmath.h"
#include "math/4dmath.h"
#include "render/ray/contract.h"

namespace SDF {

struct FrameworkPlane {
  math::Vector normal;
  float spacing;
};

/** @brief Triangular-prism edges with equilateral XY cells and Z layers. */
struct TriangularFramework {
  using PlaneFamily = FrameworkPlane;
  static constexpr size_t STREAM_COUNT = 4;
  static constexpr float TRIANGLE_HEIGHT = 0.8660254037844386f;
  float cell_size = 1.0f;
  float layer_height = 1.0f;
  float wire_radius = 0.05f;
  math::Vector origin{};

  HS_HOT_FLASH_MEMBER bool valid() const {
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

/** @brief FCC nearest-neighbor edges; cell_size is the strut length. */
struct OctetFramework {
  using PlaneFamily = FrameworkPlane;
  static constexpr size_t STREAM_COUNT = 4;
  static constexpr float NORMAL_COMPONENT = 0.5773502691896258f;
  static constexpr float PLANE_SPACING = 0.8164965809277260f;
  static constexpr float INVERSE_PLANE_SINE = 1.0606601717798213f;
  float cell_size = 1.0f;
  float wire_radius = 0.05f;
  math::Vector origin{};

  bool valid() const {
    return Raycast::finite(cell_size) && cell_size > 0.0f &&
           Raycast::finite(wire_radius) && wire_radius > 0.0f &&
           Raycast::finite(origin);
  }

  std::array<PlaneFamily, STREAM_COUNT> plane_families() const {
    const float N = NORMAL_COMPONENT;
    const float SPACING = PLANE_SPACING * cell_size;
    return {{{{N, N, N}, SPACING},
             {{N, -N, -N}, SPACING},
             {{-N, N, -N}, SPACING},
             {{-N, -N, N}, SPACING}}};
  }

  math::Vector edge_offset(const math::Vector &p, uint32_t &feature) const {
    const auto FAMILIES = plane_families();
    const math::Vector Q = p - origin;
    const float SPACING = FAMILIES[0].spacing;
    std::array<float, STREAM_COUNT> residual;
    for (size_t i = 0; i < STREAM_COUNT; ++i) {
      const float U = math::dot(Q, FAMILIES[i].normal);
      residual[i] = U - SPACING * roundf(U / SPACING);
    }
    math::Vector result{};
    float best = INFINITY;
    uint32_t pair = 0;
    for (size_t i = 0; i < STREAM_COUNT; ++i) {
      for (size_t j = i + 1; j < STREAM_COUNT; ++j, ++pair) {
        float a = residual[i];
        float b = residual[j];
        if (a * b > 0.0f) {
          float &larger = fabsf(a) > fabsf(b) ? a : b;
          const float SMALLER = std::min(fabsf(a), fabsf(b));
          if (fabsf(larger) + SMALLER / 3.0f > 0.5f * SPACING)
            larger -= copysignf(SPACING, larger);
        }
        const float SQUARED = 1.125f * (a * a + b * b + (2.0f / 3.0f) * a * b);
        if (SQUARED < best) {
          best = SQUARED;
          result = FAMILIES[i].normal * (1.125f * a + 0.375f * b) +
                   FAMILIES[j].normal * (1.125f * b + 0.375f * a);
          feature = pair;
        }
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

  /** @brief Distance to struts contained in the selected plane through p. */
  Raycast::QuerySample plane_sample(size_t index, const math::Vector &p) const {
    const auto FAMILIES = plane_families();
    const math::Vector Q = p - origin;
    const float SPACING = FAMILIES[0].spacing;
    float best = INFINITY;
    uint32_t feature = 0;
    uint32_t pair = 0;
    for (size_t i = 0; i < STREAM_COUNT; ++i) {
      for (size_t j = i + 1; j < STREAM_COUNT; ++j, ++pair) {
        if (i != index && j != index)
          continue;
        const size_t OTHER = i == index ? j : i;
        const float U = math::dot(Q, FAMILIES[OTHER].normal);
        const float D = fabsf(U - SPACING * roundf(U / SPACING));
        if (D < best) {
          best = D;
          feature = pair;
        }
      }
    }
    const float VALUE = best * INVERSE_PLANE_SINE - wire_radius;
    return {VALUE, 0.0f, false, 0, feature};
  }
};

/** @brief D4 nearest-neighbor struts joining integer vertices of even sum. */
struct OctetFramework4 {
  struct PlaneFamily {
    math::Vec4 normal;
    float spacing;
  };
  static constexpr size_t STREAM_COUNT = 8;
  static constexpr float HALF_CUBE = 0.7071067811865475f;
  float cell_size = 1.0f;
  float wire_radius = 0.05f;
  math::Vec4 origin{};

  bool valid() const {
    if (!Raycast::finite(cell_size) || cell_size <= 0.0f ||
        !Raycast::finite(wire_radius) || wire_radius <= 0.0f)
      return false;
    for (int i = 0; i < 4; ++i)
      if (!Raycast::finite(origin[i]))
        return false;
    return true;
  }

  std::array<PlaneFamily, STREAM_COUNT> plane_families() const {
    std::array<PlaneFamily, STREAM_COUNT> result;
    for (size_t i = 0; i < STREAM_COUNT; ++i)
      result[i] = {{{0.5f, (i & 1) ? -0.5f : 0.5f, (i & 2) ? -0.5f : 0.5f,
                     (i & 4) ? -0.5f : 0.5f}},
                   HALF_CUBE * cell_size};
    return result;
  }

  HS_HOT_FLASH_MEMBER math::Vec4 edge_offset(const math::Vec4 &p,
                                             uint32_t &feature) const {
    const float SCALE = HALF_CUBE * cell_size;
    const float INVERSE_SCALE = 1.0f / SCALE;
    math::Vec4 q, rounded, residual, result;
    for (int i = 0; i < 4; ++i) {
      q[i] = (p[i] - origin[i]) * INVERSE_SCALE;
      rounded[i] = roundf(q[i]);
      residual[i] = q[i] - rounded[i];
    }
    float best = INFINITY;
    uint32_t pair = 0;
    feature = 0;
    for (int i = 0; i < 4; ++i)
      for (int j = i + 1; j < 4; ++j)
        for (int sign : {-1, 1}) {
          const float U = q[i] - sign * q[j];
          const float PLANE = roundf(U);
          float along = U - PLANE;
          float parity = PLANE;
          float squared = 0.0f;
          float flip_cost = 0.5f - fabsf(along);
          int flip_axis = -1;
          math::Vec4 offset = residual;
          for (int k = 0; k < 4; ++k) {
            if (k == i || k == j)
              continue;
            parity += rounded[k];
            squared += residual[k] * residual[k];
            const float COST = 1.0f - 2.0f * fabsf(residual[k]);
            if (COST < flip_cost) {
              flip_cost = COST;
              flip_axis = k;
            }
          }
          squared += 0.5f * along * along;
          // D4 requires an even sum of the line invariant and fixed coordinates.
          if (fabsf(parity - 2.0f * roundf(parity * 0.5f)) > 0.5f) {
            squared += flip_cost;
            if (flip_axis < 0)
              along -= copysignf(1.0f, along);
            else
              offset[flip_axis] -= copysignf(1.0f, offset[flip_axis]);
          }
          if (squared < best) {
            best = squared;
            offset[i] = 0.5f * along;
            offset[j] = -0.5f * sign * along;
            for (int k = 0; k < 4; ++k)
              result[k] = offset[k] * SCALE;
            feature = pair;
          }
          ++pair;
        }
    return result;
  }

  static float magnitude(const math::Vec4 &p) {
    float squared = 0.0f;
    for (int i = 0; i < 4; ++i)
      squared += p[i] * p[i];
    return sqrtf(squared);
  }

  float distance(const math::Vec4 &p) const {
    uint32_t feature;
    return magnitude(edge_offset(p, feature)) - wire_radius;
  }

  math::Vec4 normal(const math::Vec4 &p) const {
    uint32_t feature;
    math::Vec4 offset = edge_offset(p, feature);
    const float LENGTH = magnitude(offset);
    for (int i = 0; i < 4; ++i)
      offset[i] = LENGTH > 0.0f ? offset[i] / LENGTH : 0.0f;
    return offset;
  }

  Raycast::QueryCapabilities capabilities() const {
    return {true, true, true, 0.0f};
  }

  Raycast::QuerySample sample(const math::Vec4 &p) const {
    uint32_t feature;
    const float VALUE = magnitude(edge_offset(p, feature)) - wire_radius;
    return {VALUE, fabsf(VALUE), VALUE == 0.0f, 0, feature};
  }
};

struct FrameworkPlaneCursor {
  float next = 0.0f;
  float step = 0.0f;
  bool active = false;

  struct Projection {
    float position;
    float speed;
    float spacing;
  };

  HS_HOT_FLASH_MEMBER static void initialize(FrameworkPlaneCursor *cursors,
                                             const Projection *projections,
                                             size_t count, float near) {
    for (size_t i = 0; i < count; ++i) {
      const float SPEED = projections[i].speed;
      if (SPEED == 0.0f)
        continue;
      const float POSITION = projections[i].position;
      const float SPACING = projections[i].spacing;
      const float CELL = POSITION / SPACING;
      const float PLANE = SPEED > 0.0f ? ceilf(CELL) : floorf(CELL);
      auto &cursor = cursors[i];
      cursor.next = near + (PLANE * SPACING - POSITION) / SPEED;
      cursor.step = SPACING / fabsf(SPEED);
      cursor.active = Raycast::finite(cursor.next) &&
                      Raycast::finite(cursor.step) && cursor.step > 0.0f;
    }
  }
};

/** @brief Shared monotone plane streams and approximate coverage layers. */
template <size_t Count> struct FrameworkPlaneStreams {
  static constexpr size_t STREAM_COUNT = Count;
  Raycast::Footprint footprint;
  std::array<FrameworkPlaneCursor, STREAM_COUNT> cursors{};

  explicit FrameworkPlaneStreams(Raycast::Footprint footprint)
      : footprint(footprint) {}

  bool active(size_t index) const { return cursors[index].active; }
  float distance(size_t index) const { return cursors[index].next; }

  Raycast::Contribution contribution(size_t index,
                                     const Raycast::QuerySample &sample) const {
    Raycast::Contribution result;
    result.t = distance(index);
    const float WIDTH = footprint.at(result.t);
    result.coverage = WIDTH > 0.0f
                          ? std::clamp(0.5f - sample.field / WIDTH, 0.0f, 1.0f)
                          : (sample.field <= 0.0f ? 1.0f : 0.0f);
    result.feature = sample.feature;
    // Coincident plane reports define one approximate junction layer.
    result.merge_identity = 0;
    return result;
  }

  void advance(size_t index) {
    auto &cursor = cursors[index];
    const float NEXT = cursor.next + cursor.step;
    cursor.active = Raycast::finite(NEXT) && NEXT > cursor.next;
    cursor.next = NEXT;
  }
};

struct FrameworkPlaneEvents : FrameworkPlaneStreams<4> {
  Raycast::Ray ray;

  FrameworkPlaneEvents(const Raycast::Ray &ray, Raycast::Footprint footprint)
      : FrameworkPlaneStreams(footprint), ray(ray) {}

  void initialize(const std::array<FrameworkPlane, STREAM_COUNT> &families,
                  const math::Vector &origin, bool valid) {
    if (!valid || !ray.valid())
      return;
    const math::Vector START = ray.at(ray.interval.near) - origin;
    std::array<FrameworkPlaneCursor::Projection, STREAM_COUNT> projections;
    for (size_t i = 0; i < STREAM_COUNT; ++i)
      projections[i] = {math::dot(START, families[i].normal),
                        math::dot(ray.direction, families[i].normal),
                        families[i].spacing};
    FrameworkPlaneCursor::initialize(cursors.data(), projections.data(),
                                     STREAM_COUNT, ray.interval.near);
  }
};

/** @brief Four plane streams producing approximate framework coverage layers. */
struct FrameworkEvents : FrameworkPlaneEvents {
  TriangularFramework geometry;

  FrameworkEvents(const TriangularFramework &geometry, const Raycast::Ray &ray,
                  Raycast::Footprint footprint = {})
      : FrameworkPlaneEvents(ray, footprint), geometry(geometry) {
    initialize(geometry.plane_families(), geometry.origin, geometry.valid());
  }

  Raycast::Contribution candidate(size_t index) const {
    return contribution(index, geometry.sample(ray.at(distance(index))));
  }
};

/** @brief Approximate octet coverage from struts in each crossed plane. */
struct OctetEvents : FrameworkPlaneEvents {
  OctetFramework geometry;

  OctetEvents(const OctetFramework &geometry, const Raycast::Ray &ray,
              Raycast::Footprint footprint = {})
      : FrameworkPlaneEvents(ray, footprint), geometry(geometry) {
    initialize(geometry.plane_families(), geometry.origin, geometry.valid());
  }

  Raycast::Contribution candidate(size_t index) const {
    return contribution(index,
                        geometry.plane_sample(index, ray.at(distance(index))));
  }
};

/** @brief Approximate D4 strut coverage evaluated along an ambient 4D ray. */
struct OctetEvents4 : FrameworkPlaneStreams<8> {
  const OctetFramework4 &geometry;
  math::Vec4 origin;
  math::Vec4 direction;

  OctetEvents4(const OctetFramework4 &geometry, const math::Vec4 &origin,
               const math::Vec4 &direction, Raycast::Interval interval,
               Raycast::Footprint footprint = {})
      : FrameworkPlaneStreams(footprint), geometry(geometry), origin(origin),
        direction(direction) {
    if (!geometry.valid() || !interval.valid())
      return;
    float length2 = 0.0f;
    math::Vec4 start;
    for (int i = 0; i < 4; ++i) {
      if (!Raycast::finite(origin[i]) || !Raycast::finite(direction[i]))
        return;
      length2 += direction[i] * direction[i];
      start[i] = origin[i] + direction[i] * interval.near - geometry.origin[i];
    }
    if (fabsf(length2 - 1.0f) >= 1e-4f)
      return;
    const auto FAMILIES = geometry.plane_families();
    std::array<FrameworkPlaneCursor::Projection, STREAM_COUNT> projections{};
    for (size_t i = 0; i < STREAM_COUNT; ++i) {
      projections[i].spacing = FAMILIES[i].spacing;
      for (int j = 0; j < 4; ++j) {
        projections[i].position += start[j] * FAMILIES[i].normal[j];
        projections[i].speed += direction[j] * FAMILIES[i].normal[j];
      }
    }
    FrameworkPlaneCursor::initialize(cursors.data(), projections.data(),
                                     STREAM_COUNT, interval.near);
  }

  Raycast::Contribution candidate(size_t index) const {
    math::Vec4 p;
    const float T = distance(index);
    for (int i = 0; i < 4; ++i)
      p[i] = origin[i] + direction[i] * T;
    return contribution(index, geometry.sample(p));
  }
};

} // namespace SDF
