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

  template <bool WithOffset, bool Normalized = false>
  HS_HOT_FLASH_MEMBER auto edge_query(const math::Vec4 &p,
                                      uint32_t &feature) const {
    const float SCALE = Normalized ? 1.0f : HALF_CUBE * cell_size;
    const float INVERSE_SCALE = 1.0f / SCALE;
    math::Vec4 residual, absolute;
    float parity = 0.0f;
    int halves = 0;
    for (int k = 0; k < 4; ++k) {
      const float Q = Normalized ? p[k] : (p[k] - origin[k]) * INVERSE_SCALE;
      const float ROUNDED = roundf(Q);
      parity += ROUNDED;
      residual[k] = Q - ROUNDED;
      absolute[k] = fabsf(residual[k]);
      halves += absolute[k] == 0.5f;
    }
    const bool ODD = fabsf(parity - 2.0f * roundf(parity * 0.5f)) > 0.5f;
    // The nearest strut uses the two largest fractional coordinates.
    int i = 0;
    int j = 1;
    if (absolute[j] > absolute[i])
      std::swap(i, j);
    float transverse_squared = 0.0f;
    for (int k = 2; k < 4; ++k) {
      int transverse = k;
      if (absolute[k] > absolute[i]) {
        transverse = j;
        j = i;
        i = k;
      } else if (absolute[k] > absolute[j]) {
        transverse = j;
        j = k;
      }
      if constexpr (!WithOffset)
        transverse_squared += residual[transverse] * residual[transverse];
    }
    if (i > j)
      std::swap(i, j);
    feature = 2 * (i * (7 - i) / 2 + j - i - 1);
    float along = 0.0f;
    if (halves >= 3) {
      // Three half-grid coordinates tie both strut orientations.
      if constexpr (WithOffset) {
        if (ODD != (residual[i] + residual[j] != 0.0f))
          for (int k = 0; k < 4; ++k)
            if (k != i && k != j && absolute[k] == 0.5f) {
              residual[k] = -residual[k];
              break;
            }
        residual[i] = 0.0f;
        residual[j] = 0.0f;
      }
    } else {
      const bool SAME_SIGN = (residual[i] < 0.0f) == (residual[j] < 0.0f);
      const float SIGN =
          absolute[i] != 0.0f && absolute[j] != 0.0f && SAME_SIGN != ODD
              ? 1.0f
              : -1.0f;
      along = residual[i] - SIGN * residual[j];
      if (ODD)
        along -= copysignf(1.0f, along);
      if constexpr (WithOffset) {
        residual[i] = 0.5f * along;
        residual[j] = -0.5f * SIGN * along;
      }
      feature += SIGN > 0.0f;
    }
    if constexpr (WithOffset) {
      for (int k = 0; k < 4; ++k)
        residual[k] *= SCALE;
      return residual;
    } else {
      return sqrtf(transverse_squared + 0.5f * along * along) * SCALE;
    }
  }

  math::Vec4 edge_offset(const math::Vec4 &p, uint32_t &feature) const {
    return edge_query<true>(p, feature);
  }

  static float magnitude(const math::Vec4 &p) {
    float squared = 0.0f;
    for (int i = 0; i < 4; ++i)
      squared += p[i] * p[i];
    return sqrtf(squared);
  }

  float distance(const math::Vec4 &p) const {
    uint32_t feature;
    return edge_query<false>(p, feature) - wire_radius;
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
    const float VALUE = edge_query<false>(p, feature) - wire_radius;
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
  static constexpr size_t GROUP_CAPACITY = 1;
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

/** @brief Single-owner strut coverage at monotone octet plane crossings. */
struct OctetEvents : FrameworkPlaneStreams<4> {
  /** @brief Validated frame projections into the four octet plane families. */
  struct PreparedProjection {
    std::array<math::Vector, STREAM_COUNT> normals;
    std::array<float, STREAM_COUNT> offsets;
    float spacing;
    float inverse_spacing;
    float wire_radius;
  };
  struct Pair {
    float scale = 0.0f;
    uint8_t owner = 0;
    uint8_t other = 0;
  };
  std::array<Pair, 6> pairs{};
  std::array<float, STREAM_COUNT> positions{};
  std::array<float, STREAM_COUNT> speeds{};
  float spacing = 0.0f;
  float inverse_spacing = 0.0f;
  float wire_radius = 0.0f;

  HS_HOT_FLASH_MEMBER OctetEvents(const OctetFramework &geometry,
                                  const Raycast::Ray &ray,
                                  Raycast::Footprint footprint = {})
      : FrameworkPlaneStreams(footprint) {
    if (!geometry.valid() || !ray.valid())
      return;
    const auto FAMILIES = geometry.plane_families();
    spacing = FAMILIES[0].spacing;
    inverse_spacing = 1.0f / spacing;
    wire_radius = geometry.wire_radius;
    const math::Vector ORIGIN = ray.origin - geometry.origin;
    for (size_t i = 0; i < STREAM_COUNT; ++i) {
      positions[i] = math::dot(ORIGIN, FAMILIES[i].normal);
      speeds[i] = math::dot(ray.direction, FAMILIES[i].normal);
    }
    initialize(ray.interval.near);
  }

  /** @brief Initializes a validated unit ray from validated frame projections. */
  HS_HOT_FLASH_MEMBER
  OctetEvents(const PreparedProjection &prepared, const math::Vector &direction,
              float radial_start, Raycast::Interval interval,
              Raycast::Footprint footprint = {})
      : FrameworkPlaneStreams(footprint), spacing(prepared.spacing),
        inverse_spacing(prepared.inverse_spacing),
        wire_radius(prepared.wire_radius) {
    for (size_t i = 0; i < STREAM_COUNT; ++i) {
      speeds[i] = math::dot(direction, prepared.normals[i]);
      positions[i] = prepared.offsets[i] + radial_start * speeds[i];
    }
    initialize(interval.near);
  }

  __attribute__((always_inline)) inline void initialize(float near) {
    std::array<FrameworkPlaneCursor::Projection, STREAM_COUNT> projections;
    for (size_t i = 0; i < STREAM_COUNT; ++i)
      projections[i] = {positions[i] + near * speeds[i], speeds[i], spacing};
    size_t pair = 0;
    std::array<bool, STREAM_COUNT> owned{};
    for (uint8_t i = 0; i < STREAM_COUNT; ++i)
      for (uint8_t j = i + 1; j < STREAM_COUNT; ++j, ++pair) {
        const float A = speeds[i];
        const float B = speeds[j];
        const uint8_t OWNER = fabsf(A) >= fabsf(B) ? i : j;
        if (speeds[OWNER] == 0.0f)
          continue;
        // Ray-to-line distance for tetrahedral plane normals (dot = -1/3).
        pairs[pair] = {speeds[OWNER] * speeds[OWNER] /
                           (A * A + B * B + (2.0f / 3.0f) * A * B),
                       OWNER, OWNER == i ? j : i};
        owned[OWNER] = true;
      }
    FrameworkPlaneCursor::initialize(cursors.data(), projections.data(),
                                     STREAM_COUNT, near);
    for (size_t i = 0; i < STREAM_COUNT; ++i) {
      cursors[i].active = cursors[i].active && owned[i];
      positions[i] *= inverse_spacing;
      speeds[i] *= inverse_spacing;
    }
    for (auto &pair : pairs)
      pair.scale *= spacing * spacing;
  }

  Raycast::Contribution candidate(size_t index) const {
    const float T = distance(index);
    float best = INFINITY;
    uint32_t feature = 0;
    for (uint32_t i = 0; i < pairs.size(); ++i) {
      const auto &pair = pairs[i];
      if (pair.owner != index || pair.scale == 0.0f)
        continue;
      const float U = positions[pair.other] + T * speeds[pair.other];
      const float RESIDUAL = U - roundf(U);
      const float D = RESIDUAL * RESIDUAL * pair.scale;
      if (D < best) {
        best = D;
        feature = i;
      }
    }
    const float WIDTH = footprint.at(T);
    const float SUPPORT = wire_radius + .5f * WIDTH;
    if (best > SUPPORT * SUPPORT) {
      Raycast::Contribution miss;
      miss.t = T;
      miss.coverage = 0;
      return miss;
    }
    return contribution(index,
                        {sqrtf(best) - wire_radius, 0.0f, false, 0, feature});
  }
};

/** @brief Approximate D4 strut coverage evaluated along an ambient 4D ray. */
struct OctetEvents4 : FrameworkPlaneStreams<8> {
  const OctetFramework4 &geometry;
  math::Vec4 local_origin;
  math::Vec4 local_direction;
  float scale = 0.0f;

  OctetEvents4(const OctetFramework4 &geometry, const math::Vec4 &origin,
               const math::Vec4 &direction, Raycast::Interval interval,
               Raycast::Footprint footprint = {})
      : FrameworkPlaneStreams(footprint), geometry(geometry) {
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
    scale = OctetFramework4::HALF_CUBE * geometry.cell_size;
    const float INVERSE_SCALE = 1.0f / scale;
    for (int i = 0; i < 4; ++i) {
      local_origin[i] = (origin[i] - geometry.origin[i]) * INVERSE_SCALE;
      local_direction[i] = direction[i] * INVERSE_SCALE;
    }
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
      p[i] = local_origin[i] + local_direction[i] * T;
    uint32_t feature;
    const float VALUE = scale * geometry.edge_query<false, true>(p, feature) -
                        geometry.wire_radius;
    return contribution(index, {VALUE, 0.0f, false, 0, feature});
  }
};

} // namespace SDF
