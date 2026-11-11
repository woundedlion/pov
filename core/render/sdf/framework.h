/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file framework.h
 * @brief Lattice framework geometries and monotone plane-crossing event streams. */

#include <algorithm>
#include <array>
#include <cmath>
#include "core/math/3dmath.h"
#include "core/math/4dmath.h"
#include "core/render/ray/contract.h"

namespace SDF {

/**
 * @brief Footprint coverage of a signed field value or interior depth.
 * @tparam POSITIVE_INSIDE Whether positive values lie inside the surface.
 * @param width Full ramp span in world units. Nonpositive values select a hard step.
 * @param field Signed value in world units; negative inside by default.
 * @return Coverage in [0, 1].
 */
template <bool POSITIVE_INSIDE = false>
__attribute__((always_inline)) inline float footprint_coverage(float width,
                                                               float field) {
  constexpr float MIDPOINT = 0.5f;
  if (width > 0.0f) {
    const float ratio = field / width;
    return hs::clamp(POSITIVE_INSIDE ? MIDPOINT + ratio : MIDPOINT - ratio,
                     0.0f, 1.0f);
  }
  return (POSITIVE_INSIDE ? field >= 0.0f : field <= 0.0f) ? 1.0f : 0.0f;
}

/** @brief One family of parallel planes. */
struct FrameworkPlane {
  math::Vector normal; ///< Unit plane normal.
  float spacing;       ///< Distance between adjacent planes.
};

namespace detail {
/** @brief Exact interior and exterior clearance of framework strut fields. */
inline constexpr Raycast::QueryCapabilities FRAMEWORK_CAPABILITIES{true, true,
                                                                   true, 0.0f};

/** @brief Field sample of a framework's nearest strut.
 * @tparam Framework Geometry exposing edge_offset and wire_radius.
 * @param framework Framework to query.
 * @param p World-space point.
 * @return Signed distance, clearance and nearest-edge identity. */
template <typename Framework>
inline Raycast::QuerySample framework_sample(const Framework &framework,
                                             const math::Vector &p) {
  uint32_t feature;
  const float value =
      framework.edge_offset(p, feature).magnitude() - framework.wire_radius;
  return {value, fabsf(value), value == 0.0f, 0, feature};
}

/** @brief Outward normal of a framework's nearest strut.
 * @tparam Framework Geometry exposing edge_offset.
 * @param framework Framework to query.
 * @param p World-space point.
 * @return Unit normal, or zero on an edge axis. */
template <typename Framework>
inline math::Vector framework_normal(const Framework &framework,
                                     const math::Vector &p) {
  uint32_t feature;
  const math::Vector offset = framework.edge_offset(p, feature);
  const float length = offset.magnitude();
  return length > 0.0f ? offset * (1.0f / length) : math::Vector{};
}
} // namespace detail

/** @brief Triangular-prism edges with equilateral XY cells and Z layers. */
struct TriangularFramework {
  using PlaneFamily = FrameworkPlane;       ///< Plane family type.
  static constexpr size_t STREAM_COUNT = 4; ///< Number of plane families.
  /// Row height per unit cell size (sqrt(3)/2).
  static constexpr float TRIANGLE_HEIGHT = 0.8660254037844386f;
  float cell_size = 1.0f;    ///< Triangle edge length.
  float layer_height = 1.0f; ///< Spacing of the Z layers.
  float wire_radius = 0.05f; ///< Strut radius.
  math::Vector origin{};     ///< Lattice origin.

  /** @brief Whether all parameters are finite and sizes positive.
   * @return True when the geometry is usable. */
  HS_HOT_FLASH_MEMBER bool valid() const {
    return Raycast::finite(cell_size) && cell_size > 0.0f &&
           Raycast::finite(layer_height) && layer_height > 0.0f &&
           Raycast::finite(wire_radius) && wire_radius > 0.0f &&
           Raycast::finite(origin);
  }

  /** @brief Three XY edge-line families and the Z layer planes.
   * @return Unit normals and spacings of the four families. */
  std::array<PlaneFamily, STREAM_COUNT> plane_families() const {
    const float height = TRIANGLE_HEIGHT * cell_size;
    return {{{{0.0f, 1.0f, 0.0f}, height},
             {{TRIANGLE_HEIGHT, 0.5f, 0.0f}, height},
             {{-TRIANGLE_HEIGHT, 0.5f, 0.0f}, height},
             {{0.0f, 0.0f, 1.0f}, layer_height}}};
  }

  /**
   * @brief Offset from the nearest edge to a point.
   * @param p Query point.
   * @param feature Set to the edge id: 0-2 horizontal family, 3 vertical.
   * @return Vector from the nearest edge point to p.
   */
  math::Vector edge_offset(const math::Vector &p, uint32_t &feature) const {
    const math::Vector q = p - origin;
    const float height = TRIANGLE_HEIGHT * cell_size;
    const float b_value = roundf(q.y / height);
    const float a_value = roundf(q.x / cell_size - 0.5f * b_value);
    math::Vector result{cell_size, cell_size, 0.0f};
    float best = math::dot(result, result);
    feature = 3;
    for (int b = -1; b <= 1; ++b) {
      for (int a = -1; a <= 1; ++a) {
        const math::Vector offset{
            q.x - cell_size * (a_value + a + 0.5f * (b_value + b)),
            q.y - height * (b_value + b), 0.0f};
        const float squared = math::dot(offset, offset);
        if (squared < best) {
          best = squared;
          result = offset;
        }
      }
    }
    const float z = q.z - layer_height * roundf(q.z / layer_height);
    const auto families = plane_families();
    for (uint32_t i = 0; i < 3; ++i) {
      const float u = math::dot(q, families[i].normal);
      const float d = u - height * roundf(u / height);
      const float squared = d * d + z * z;
      if (squared < best) {
        best = squared;
        result = families[i].normal * d + math::Vector{0.0f, 0.0f, z};
        feature = i;
      }
    }
    return result;
  }

  /** @brief Signed distance to the strut surface.
   * @param p Query point.
   * @return Distance to the nearest edge minus the wire radius. */
  float distance(const math::Vector &p) const {
    return detail::framework_sample(*this, p).field;
  }

  /** @brief Outward surface normal.
   * @param p Query point.
   * @return Unit normal, or zero on an edge axis. */
  math::Vector normal(const math::Vector &p) const {
    return detail::framework_normal(*this, p);
  }

  /** @brief Exact field with interior and exterior clearance.
   * @return The query capabilities. */
  Raycast::QueryCapabilities capabilities() const {
    return detail::FRAMEWORK_CAPABILITIES;
  }

  /** @brief Field sample with the nearest edge as feature.
   * @param p Query point.
   * @return Signed distance, clearance and edge id. */
  Raycast::QuerySample sample(const math::Vector &p) const {
    return detail::framework_sample(*this, p);
  }
};

/** @brief FCC nearest-neighbor edges; cell_size is the strut length. */
struct OctetFramework {
  using PlaneFamily = FrameworkPlane;       ///< Plane family type.
  static constexpr size_t STREAM_COUNT = 4; ///< Number of plane families.
  /// Normal component magnitude (1/sqrt(3)).
  static constexpr float NORMAL_COMPONENT = 0.5773502691896258f;
  /// Plane spacing per unit cell size.
  static constexpr float PLANE_SPACING = 0.8164965809277260f;
  /// Reciprocal sine of the dihedral angle between plane families.
  static constexpr float INVERSE_PLANE_SINE = 1.0606601717798213f;
  float cell_size = 1.0f;    ///< Strut length.
  float wire_radius = 0.05f; ///< Strut radius.
  math::Vector origin{};     ///< Lattice origin.

  /** @brief Whether all parameters are finite and sizes positive.
   * @return True when the geometry is usable. */
  bool valid() const {
    return Raycast::finite(cell_size) && cell_size > 0.0f &&
           Raycast::finite(wire_radius) && wire_radius > 0.0f &&
           Raycast::finite(origin);
  }

  /** @brief The four tetrahedral plane families containing the struts.
   * @return Unit normals and spacings of the four families. */
  std::array<PlaneFamily, STREAM_COUNT> plane_families() const {
    const float N = NORMAL_COMPONENT;
    const float spacing = PLANE_SPACING * cell_size;
    return {{{{N, N, N}, spacing},
             {{N, -N, -N}, spacing},
             {{-N, N, -N}, spacing},
             {{-N, -N, N}, spacing}}};
  }

  /**
   * @brief Offset from the nearest strut to a point.
   * @param p Query point.
   * @param feature Set to the strut id: the plane-family pair index 0-5.
   * @return Vector from the nearest strut point to p.
   */
  math::Vector edge_offset(const math::Vector &p, uint32_t &feature) const {
    const auto families = plane_families();
    const math::Vector q = p - origin;
    const float spacing = families[0].spacing;
    std::array<float, STREAM_COUNT> residual;
    for (size_t i = 0; i < STREAM_COUNT; ++i) {
      const float u = math::dot(q, families[i].normal);
      residual[i] = u - spacing * roundf(u / spacing);
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
          const float smaller = fminf(fabsf(a), fabsf(b));
          if (fabsf(larger) + smaller / 3.0f > 0.5f * spacing)
            larger -= copysignf(spacing, larger);
        }
        const float squared = 1.125f * (a * a + b * b + (2.0f / 3.0f) * a * b);
        if (squared < best) {
          best = squared;
          result = families[i].normal * (1.125f * a + 0.375f * b) +
                   families[j].normal * (1.125f * b + 0.375f * a);
          feature = pair;
        }
      }
    }
    return result;
  }

  /** @brief Signed distance to the strut surface.
   * @param p Query point.
   * @return Distance to the nearest strut minus the wire radius. */
  float distance(const math::Vector &p) const {
    return detail::framework_sample(*this, p).field;
  }

  /** @brief Outward surface normal.
   * @param p Query point.
   * @return Unit normal, or zero on a strut axis. */
  math::Vector normal(const math::Vector &p) const {
    return detail::framework_normal(*this, p);
  }

  /** @brief Exact field with interior and exterior clearance.
   * @return The query capabilities. */
  Raycast::QueryCapabilities capabilities() const {
    return detail::FRAMEWORK_CAPABILITIES;
  }

  /** @brief Field sample with the nearest strut as feature.
   * @param p Query point.
   * @return Signed distance, clearance and strut id. */
  Raycast::QuerySample sample(const math::Vector &p) const {
    return detail::framework_sample(*this, p);
  }
};

namespace detail::reference {
/**
   * @brief Distance to struts contained in the selected plane through p.
   * @param framework Octet framework to sample.
   * @param index Plane family index, below OctetFramework::STREAM_COUNT.
   * @param p Query point.
   * @return Signed distance and strut id; clearance is not provided.
   */
inline Raycast::QuerySample plane_sample(const OctetFramework &framework,
                                         size_t index, const math::Vector &p) {
  const auto families = framework.plane_families();
  const math::Vector q = p - framework.origin;
  const float spacing = families[0].spacing;
  float best = INFINITY;
  uint32_t feature = 0;
  uint32_t pair = 0;
  for (size_t i = 0; i < OctetFramework::STREAM_COUNT; ++i) {
    for (size_t j = i + 1; j < OctetFramework::STREAM_COUNT; ++j, ++pair) {
      if (i != index && j != index)
        continue;
      const size_t other = i == index ? j : i;
      const float u = math::dot(q, families[other].normal);
      const float d = fabsf(u - spacing * roundf(u / spacing));
      if (d < best) {
        best = d;
        feature = pair;
      }
    }
  }
  const float value =
      best * OctetFramework::INVERSE_PLANE_SINE - framework.wire_radius;
  return {value, 0.0f, false, 0, feature};
}
} // namespace detail::reference

/** @brief D4 nearest-neighbor struts joining integer vertices of even sum. */
struct OctetFramework4 {
  /** @brief One family of parallel 3-planes. */
  struct PlaneFamily {
    math::Vec4 normal; ///< Unit plane normal.
    float spacing;     ///< Distance between adjacent planes.
  };
  static constexpr size_t STREAM_COUNT = 8; ///< Number of plane families.
  /// Lattice unit per unit cell size (1/sqrt(2)).
  static constexpr float HALF_CUBE = 0.7071067811865475f;
  float cell_size = 1.0f;    ///< Strut length.
  float wire_radius = 0.05f; ///< Strut radius.
  math::Vec4 origin{};       ///< Lattice origin.

  /** @brief Whether all parameters are finite and sizes positive.
   * @return True when the geometry is usable. */
  bool valid() const {
    if (!Raycast::finite(cell_size) || cell_size <= 0.0f ||
        !Raycast::finite(wire_radius) || wire_radius <= 0.0f)
      return false;
    for (int i = 0; i < 4; ++i)
      if (!Raycast::finite(origin[i]))
        return false;
    return true;
  }

  /** @brief The eight (1, +-1, +-1, +-1)/2 plane families.
   * @return Unit normals and spacings of the eight families. */
  std::array<PlaneFamily, STREAM_COUNT> plane_families() const {
    std::array<PlaneFamily, STREAM_COUNT> result;
    for (size_t i = 0; i < STREAM_COUNT; ++i)
      result[i] = {{{0.5f, (i & 1) ? -0.5f : 0.5f, (i & 2) ? -0.5f : 0.5f,
                     (i & 4) ? -0.5f : 0.5f}},
                   HALF_CUBE * cell_size};
    return result;
  }

  /**
   * @brief Nearest strut query.
   * @tparam WithOffset Return the offset vector instead of the distance.
   * @tparam Normalized p is already in lattice units relative to origin;
   *         results are then in lattice units.
   * @param p Query point.
   * @param feature Set to the strut id: 2 * axis-pair index + orientation.
   * @return Vec4 from the nearest strut point to p, or its length.
   */
  template <bool WithOffset, bool Normalized = false>
  HS_HOT_FLASH_MEMBER auto edge_query(const math::Vec4 &p,
                                      uint32_t &feature) const {
    const float scale = Normalized ? 1.0f : HALF_CUBE * cell_size;
    const float inverse_scale = 1.0f / scale;
    math::Vec4 residual, absolute;
    float parity = 0.0f;
    int halves = 0;
    for (int k = 0; k < 4; ++k) {
      const float q = Normalized ? p[k] : (p[k] - origin[k]) * inverse_scale;
      const float rounded = roundf(q);
      parity += rounded;
      residual[k] = q - rounded;
      absolute[k] = fabsf(residual[k]);
      halves += absolute[k] == 0.5f;
    }
    const bool odd = fabsf(parity - 2.0f * roundf(parity * 0.5f)) > 0.5f;
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
    float across = 0.0f;
    if (halves >= 3) {
      // Three half-grid coordinates tie both strut orientations.
      if constexpr (WithOffset) {
        if (odd != (residual[i] + residual[j] != 0.0f))
          for (int k = 0; k < 4; ++k)
            if (k != i && k != j && absolute[k] == 0.5f) {
              residual[k] = -residual[k];
              break;
            }
        residual[i] = 0.0f;
        residual[j] = 0.0f;
      }
    } else {
      const bool same_sign = (residual[i] < 0.0f) == (residual[j] < 0.0f);
      const float sign =
          absolute[i] != 0.0f && absolute[j] != 0.0f && same_sign != odd
              ? 1.0f
              : -1.0f;
      across = residual[i] - sign * residual[j];
      if (odd)
        across -= copysignf(1.0f, across);
      if constexpr (WithOffset) {
        residual[i] = 0.5f * across;
        residual[j] = -0.5f * sign * across;
      }
      feature += sign > 0.0f;
    }
    if constexpr (WithOffset) {
      for (int k = 0; k < 4; ++k)
        residual[k] *= scale;
      return residual;
    } else {
      return sqrtf(transverse_squared + 0.5f * across * across) * scale;
    }
  }

  /** @brief Offset from the nearest strut to a point.
   * @param p Query point.
   * @param feature Set to the strut id.
   * @return Vector from the nearest strut point to p. */
  math::Vec4 edge_offset(const math::Vec4 &p, uint32_t &feature) const {
    return edge_query<true>(p, feature);
  }

  /** @brief Euclidean length.
   * @param p Vector.
   * @return Length of p. */
  static float magnitude(const math::Vec4 &p) {
    float squared = 0.0f;
    for (int i = 0; i < 4; ++i)
      squared += p[i] * p[i];
    return sqrtf(squared);
  }

  /** @brief Signed distance to the strut surface.
   * @param p Query point.
   * @return Distance to the nearest strut minus the wire radius. */
  float distance(const math::Vec4 &p) const {
    uint32_t feature;
    return edge_query<false>(p, feature) - wire_radius;
  }

  /** @brief Outward surface normal.
   * @param p Query point.
   * @return Unit normal, or zero on a strut axis. */
  math::Vec4 normal(const math::Vec4 &p) const {
    uint32_t feature;
    math::Vec4 offset = edge_offset(p, feature);
    const float length = magnitude(offset);
    for (int i = 0; i < 4; ++i)
      offset[i] = length > 0.0f ? offset[i] / length : 0.0f;
    return offset;
  }

  /** @brief Exact field with interior and exterior clearance.
   * @return The query capabilities. */
  Raycast::QueryCapabilities capabilities() const {
    return {true, true, true, 0.0f};
  }

  /** @brief Field sample with the nearest strut as feature.
   * @param p Query point.
   * @return Signed distance, clearance and strut id. */
  Raycast::QuerySample sample(const math::Vec4 &p) const {
    uint32_t feature;
    const float value = edge_query<false>(p, feature) - wire_radius;
    return {value, fabsf(value), value == 0.0f, 0, feature};
  }
};

/** @brief Monotone plane streams for periodic framework adapters. */
template <size_t Count> struct PlaneStreams {
  static constexpr size_t STREAM_COUNT = Count; ///< Number of plane streams.
  static constexpr size_t GROUP_CAPACITY = 1; ///< One merge identity per group.
  Raycast::Footprint footprint;               ///< Pixel footprint for coverage.
  std::array<float, STREAM_COUNT> next; ///< Ray distance of each next crossing.
  std::array<float, STREAM_COUNT> step; ///< Ray distance between crossings.
  std::array<bool, STREAM_COUNT> live;  ///< Whether each stream is active.

  /** @brief Whether a stream has a pending crossing.
   * @param index Stream index.
   * @return True when active. */
  bool active(size_t index) const { return live[index]; }
  /** @brief Ray distance of a stream's next crossing.
   * @param index Stream index.
   * @return Distance of the pending crossing. */
  float distance(size_t index) const { return next[index]; }

  /** @brief Moves a stream to its next crossing.
   * @param index Stream index. */
  void advance(size_t index) {
    const float next_value = next[index] + step[index];
    live[index] = Raycast::finite(next_value) && next_value > next[index];
    next[index] = next_value;
  }

  /**
   * @brief Starts a stream at the first plane at or beyond near.
   * @param index Stream index.
   * @param position Plane coordinate at near, in plane spacings.
   * @param speed Plane coordinate rate per unit distance, in plane spacings.
   * @param near Distance of the interval start.
   */
  __attribute__((always_inline)) void start(size_t index, float position,
                                            float speed, float near) {
    if (speed == 0.0f) {
      live[index] = false;
      return;
    }
    const float plane = speed > 0.0f ? ceilf(position) : floorf(position);
    const float inverse = 1.0f / speed;
    next[index] = near + (plane - position) * inverse;
    step[index] = fabsf(inverse);
    live[index] = Raycast::finite(next[index]) &&
                  Raycast::finite(step[index]) && step[index] > 0.0f;
  }

  /** @brief Footprint coverage of a field value at a ray distance.
   * @param t Ray distance.
   * @param field Signed field value in world units.
   * @return Coverage in [0, 1]. */
  __attribute__((always_inline)) float coverage_of(float t, float field) const {
    const float width = footprint.at(t);
    return footprint_coverage(width, field);
  }

  /**
   * @brief Coverage layer at a stream's next crossing.
   * @param index Stream index.
   * @param sample Geometry sample at the crossing.
   * @return Contribution with footprint coverage of the sample field.
   */
  Raycast::Contribution contribution(size_t index,
                                     const Raycast::QuerySample &sample) const {
    Raycast::Contribution result;
    result.t = next[index];
    result.coverage = coverage_of(result.t, sample.field);
    result.feature = sample.feature;
    return result;
  }
};

/** @brief Four plane streams started from a world-space ray. */
struct FrameworkPlaneEvents : PlaneStreams<4> {
  Raycast::Ray ray; ///< Traced ray.

  /** @brief Constructs inactive plane streams.
   * @param ray Traced ray.
   * @param footprint Pixel footprint for coverage. */
  FrameworkPlaneEvents(const Raycast::Ray &ray, Raycast::Footprint footprint)
      : ray(ray) {
    this->footprint = footprint;
    live.fill(false);
  }

  /** @brief Starts the streams at the ray interval start.
   * @param families Plane families of the geometry.
   * @param origin Lattice origin. */
  void initialize(const std::array<FrameworkPlane, STREAM_COUNT> &families,
                  const math::Vector &origin) {
    HS_CHECK(ray.valid(), "framework event ray must be valid");
    const math::Vector start_value = ray.at(ray.interval.near) - origin;
    for (size_t i = 0; i < STREAM_COUNT; ++i) {
      const auto &family = families[i];
      start(i, math::dot(start_value, family.normal) / family.spacing,
            math::dot(ray.direction, family.normal) / family.spacing,
            ray.interval.near);
    }
  }
};

/** @brief Four plane streams producing approximate framework coverage layers. */
struct FrameworkEvents : FrameworkPlaneEvents {
  TriangularFramework geometry; ///< Framework being traced.

  /** @brief Starts the plane streams for a ray.
   * @param geometry Framework to trace; must be valid.
   * @param ray Traced ray; must be valid.
   * @param footprint Pixel footprint for coverage. */
  FrameworkEvents(const TriangularFramework &geometry, const Raycast::Ray &ray,
                  Raycast::Footprint footprint = {})
      : FrameworkPlaneEvents(ray, footprint), geometry(geometry) {
    HS_CHECK(geometry.valid(), "framework event geometry must be valid");
    initialize(geometry.plane_families(), geometry.origin);
  }

  /** @brief Coverage layer at a stream's next crossing.
   * @param index Stream index.
   * @return Contribution sampled from the geometry. */
  Raycast::Contribution candidate(size_t index) const {
    return contribution(index, geometry.sample(ray.at(distance(index))));
  }
};

/** @brief Single-owner strut coverage at monotone octet plane crossings. */
struct OctetEvents : PlaneStreams<4> {
  /** @brief Greatest number of streams that own a strut pair. */
  static constexpr size_t OWNER_CAPACITY = STREAM_COUNT - 1;
  /** @brief Maximum strut pairs owned by one stream. */
  static constexpr size_t PAIR_CAPACITY = STREAM_COUNT - 1;
  /** @brief Validated frame projections into the four octet plane families. */
  struct PreparedProjection {
    std::array<math::Vector, STREAM_COUNT> normals; /**< Per unit spacing. */
    std::array<float, STREAM_COUNT> offsets;        /**< In plane spacings. */
    float spacing2;                                 ///< Squared plane spacing.
    float wire_radius;                              ///< Strut radius.
  };
  /** @brief A strut pair evaluated at its owner's crossings. */
  struct OwnedPair {
    float numerator;   ///< Owner speed^2 * spacing^2.
    float denominator; ///< a^2 + b^2 + (2/3)ab of the pair speeds.
    uint8_t other;     ///< Non-owner stream index.
    uint8_t pair;      ///< Strut feature id.
  };
  /// Owned pairs per stream.
  std::array<std::array<OwnedPair, PAIR_CAPACITY>, STREAM_COUNT> owned;
  /// Owned pair count per stream.
  std::array<uint8_t, STREAM_COUNT> owned_count;
  /// Plane coordinate at t = 0, in spacings.
  std::array<float, STREAM_COUNT> positions;
  /// Plane coordinate per unit distance.
  std::array<float, STREAM_COUNT> speeds;
  float wire_radius; ///< Strut radius.

  /** @brief Starts the owning streams for a ray.
   * @param geometry Framework to trace; must be valid.
   * @param ray Traced ray; must be valid.
   * @param footprint Pixel footprint for coverage. */
  OctetEvents(const OctetFramework &geometry, const Raycast::Ray &ray,
              Raycast::Footprint footprint = {}) {
    this->footprint = footprint;
    live.fill(false);
    HS_CHECK(geometry.valid(), "octet event geometry must be valid");
    HS_CHECK(ray.valid(), "octet event ray must be valid");
    const auto families = geometry.plane_families();
    const float spacing = families[0].spacing;
    wire_radius = geometry.wire_radius;
    const math::Vector origin = ray.origin - geometry.origin;
    for (size_t i = 0; i < STREAM_COUNT; ++i) {
      positions[i] = math::dot(origin, families[i].normal) / spacing;
      speeds[i] = math::dot(ray.direction, families[i].normal) / spacing;
    }
    initialize(ray.interval.near, spacing * spacing);
  }

  /**
   * @brief Initializes a unit ray from validated frame projections.
   * @param prepared Per-frame plane projections.
   * @param direction Unit view direction.
   * @param radial_start Distance from the projection origin to t = 0.
   * @param near Parameter where the streams start.
   * @param footprint Pixel footprint for coverage.
   */
  __attribute__((always_inline))
  OctetEvents(const PreparedProjection &prepared, const math::Vector &direction,
              float radial_start, float near,
              Raycast::Footprint footprint = {}) {
    this->footprint = footprint;
    wire_radius = prepared.wire_radius;
    for (size_t i = 0; i < STREAM_COUNT; ++i) {
      speeds[i] = math::dot(direction, prepared.normals[i]);
      positions[i] = prepared.offsets[i] + radial_start * speeds[i];
    }
    initialize(near, prepared.spacing2);
  }

  /**
   * @brief Assigns strut pairs to owners and starts the owning streams.
   * @param near Parameter where the streams start.
   * @param spacing2 Squared plane spacing.
   */
  __attribute__((always_inline)) void initialize(float near, float spacing2) {
    owned_count.fill(0);
    uint8_t pair = 0;
    for (uint8_t i = 0; i < STREAM_COUNT; ++i)
      for (uint8_t j = i + 1; j < STREAM_COUNT; ++j, ++pair) {
        const float a = speeds[i];
        const float b = speeds[j];
        const uint8_t owner = fabsf(a) >= fabsf(b) ? i : j;
        const float speed = speeds[owner];
        if (speed == 0.0f)
          continue;
        // Ray-to-line distance for tetrahedral plane normals (dot = -1/3).
        owned[owner][owned_count[owner]++] = {
            speed * speed * spacing2, a * a + b * b + (2.0f / 3.0f) * a * b,
            static_cast<uint8_t>(owner == i ? j : i), pair};
      }
    for (size_t i = 0; i < STREAM_COUNT; ++i) {
      if (owned_count[i])
        start(i, positions[i] + near * speeds[i], speeds[i], near);
      else
        live[i] = false;
    }
  }

  /**
   * @brief Coverage of the nearest owned strut at distance t on a stream.
   * @param index Stream index.
   * @param t Ray parameter.
   * @param feature Receives the nearest strut's feature id.
   * @return Coverage in [0, 1]; 0 when no strut reaches the footprint.
   */
  __attribute__((always_inline)) float coverage(size_t index, float t,
                                                uint32_t &feature) const {
    const auto &pairs = owned[index];
    float numerator = INFINITY;
    float denominator = 1.0f;
    for (size_t k = 0; k < owned_count[index]; ++k) {
      const auto &pair = pairs[k];
      const float u = positions[pair.other] + t * speeds[pair.other];
      const float residual = u - roundf(u);
      const float n = residual * residual * pair.numerator;
      if (n * denominator < numerator * pair.denominator) {
        numerator = n;
        denominator = pair.denominator;
        feature = pair.pair;
      }
    }
    const float support = wire_radius + .5f * footprint.at(t);
    if (numerator > support * support * denominator)
      return 0.0f;
    return coverage_of(t, sqrtf(numerator / denominator) - wire_radius);
  }

  /** @brief Coverage layer at a stream's next crossing.
   * @param index Stream index.
   * @return Contribution with the nearest owned strut's coverage. */
  Raycast::Contribution candidate(size_t index) const {
    Raycast::Contribution result;
    result.t = next[index];
    result.coverage = coverage(index, result.t, result.feature);
    return result;
  }
};

/**
 * @brief Single-owner D4 strut coverage at monotone plane crossings of an
 *        ambient 4D ray.
 * @details Each strut class (e_i + s e_j)/sqrt(2) lies in four plane families
 * and is evaluated only at crossings of the family with the largest
 * |normal . direction|. Coverage uses the ray-to-line closest approach of the
 * class's nearest strut in the crossed plane.
 */
struct OctetEvents4 : PlaneStreams<8> {
  /** @brief Greatest number of families that own a strut class. */
  static constexpr size_t OWNER_CAPACITY = 4;
  /** @brief Greatest number of strut classes one family owns. */
  static constexpr size_t CLASS_CAPACITY = 6;
  /** @brief Validated frame projection of view directions into the lattice. */
  struct PreparedProjection {
    /// Per ambient axis k, the 3D vector dotted with a view direction to
    /// give its k-th ambient component.
    std::array<math::Vector, 4> embedding;
    math::Vec4 origin;   /**< Camera center in lattice units. */
    float inverse_scale; ///< 1 / lattice unit length.
  };
  /** @brief Per-ray constants of a strut class at its owner's crossings. */
  struct OwnedClass {
    float sign;        /**< s in (e_i + s e_j). */
    float transverse;  /**< (d_i - s d_j) / 2. */
    float denominator; /**< 1 - (d . u)^2. */
    float dk;          ///< Direction component on axis k.
    float dl;          ///< Direction component on axis l.
    uint8_t i;         ///< First strut axis.
    uint8_t j;         ///< Second strut axis.
    uint8_t k;         ///< First transverse axis.
    uint8_t l;         ///< Second transverse axis.
    uint8_t feature;   ///< Strut feature id.
  };
  /** @brief Strut classes owned by one plane family. */
  struct Owner {
    uint8_t count; ///< Number of owned classes.
    std::array<OwnedClass, CLASS_CAPACITY> classes; ///< Owned classes.
  };
  /// owner_of value of an unowned family.
  static constexpr uint8_t UNOWNED = 0xff;
  std::array<Owner, OWNER_CAPACITY> owners; ///< Owned classes per owner.
  /// Owner index per family, or UNOWNED.
  std::array<uint8_t, STREAM_COUNT> owner_of;
  math::Vec4 local_origin;    /**< Ray origin in lattice units. */
  math::Vec4 local_direction; /**< Direction in lattice units per distance. */
  float scale;                ///< Lattice unit length in world units.
  float inverse_scale;        ///< 1 / scale.
  float wire_radius;          ///< Strut radius.

  /**
   * @brief Starts the owning streams for an ambient 4D ray.
   * @param geometry Framework to trace; must be valid.
   * @param origin Ray origin.
   * @param direction Unit ray direction.
   * @param interval Traced distance interval; must be valid.
   * @param footprint Pixel footprint for coverage.
   */
  OctetEvents4(const OctetFramework4 &geometry, const math::Vec4 &origin,
               const math::Vec4 &direction, Raycast::Interval interval,
               Raycast::Footprint footprint = {}) {
    this->footprint = footprint;
    live.fill(false);
    HS_CHECK(geometry.valid(), "octet4 event geometry must be valid");
    HS_CHECK(interval.valid(), "octet4 event interval must be valid");
    float length2 = 0.0f;
    for (int i = 0; i < 4; ++i) {
      HS_CHECK(Raycast::finite(origin[i]) && Raycast::finite(direction[i]),
               "octet4 event ray components must be finite");
      length2 += direction[i] * direction[i];
    }
    HS_CHECK(fabsf(length2 - 1.0f) < 1e-4f,
             "octet4 event direction must be unit length");
    scale = OctetFramework4::HALF_CUBE * geometry.cell_size;
    inverse_scale = 1.0f / scale;
    wire_radius = geometry.wire_radius;
    for (int i = 0; i < 4; ++i) {
      local_origin[i] = (origin[i] - geometry.origin[i]) * inverse_scale;
      local_direction[i] = direction[i] * inverse_scale;
    }
    initialize(direction, interval.near);
  }

  /**
   * @brief Initializes a unit view ray from validated frame projections.
   * @param geometry Framework to trace.
   * @param prepared Per-frame 4D embedding.
   * @param direction Unit view direction.
   * @param radial_start Distance from the projection origin to t = 0.
   * @param near Parameter where the streams start.
   * @param footprint Pixel footprint for coverage.
   */
  __attribute__((always_inline))
  OctetEvents4(const OctetFramework4 &geometry,
               const PreparedProjection &prepared,
               const math::Vector &direction, float radial_start, float near,
               Raycast::Footprint footprint = {})
      : scale(OctetFramework4::HALF_CUBE * geometry.cell_size),
        inverse_scale(prepared.inverse_scale),
        wire_radius(geometry.wire_radius) {
    this->footprint = footprint;
    math::Vec4 ambient;
    for (int i = 0; i < 4; ++i) {
      ambient[i] = math::dot(direction, prepared.embedding[i]);
      local_direction[i] = ambient[i] * inverse_scale;
      local_origin[i] = prepared.origin[i] + radial_start * local_direction[i];
    }
    initialize(ambient, near);
  }

  /**
   * @brief Assigns strut classes to owners and starts the owning streams.
   * @param direction Unit ambient 4D ray direction.
   * @param near Parameter where the streams start.
   */
  __attribute__((always_inline)) void initialize(const math::Vec4 &direction,
                                                 float near) {
    // The owners are the family matching the direction's sign pattern and
    // that family with one coordinate flipped.
    const bool flip = direction[0] < 0.0f;
    std::array<float, 4> sign{1.0f, 1.0f, 1.0f, 1.0f};
    uint8_t star = 0;
    for (int k = 1; k < 4; ++k)
      if ((direction[k] < 0.0f) != flip) {
        sign[k] = -1.0f;
        star |= static_cast<uint8_t>(1u << (k - 1));
      }
    owner_of.fill(UNOWNED);
    uint8_t count = 0;
    const auto add = [&](uint8_t family, uint8_t i, uint8_t j, uint8_t k,
                         uint8_t l, float s,
                         uint8_t feature) __attribute__((always_inline)) {
      const float along = direction[i] + s * direction[j];
      const float denominator = 1.0f - 0.5f * along * along;
      if (!(denominator > 0.0f))
        return;
      if (owner_of[family] == UNOWNED) {
        owner_of[family] = count;
        owners[count++].count = 0;
      }
      auto &owner = owners[owner_of[family]];
      auto &strut = owner.classes[owner.count++];
      strut.sign = s;
      strut.transverse = 0.5f * (direction[i] - s * direction[j]);
      strut.denominator = denominator;
      strut.dk = direction[k];
      strut.dl = direction[l];
      strut.i = i;
      strut.j = j;
      strut.k = k;
      strut.l = l;
      strut.feature = feature;
    };
    uint8_t pair = 0;
    for (uint8_t i = 0; i < 4; ++i)
      for (uint8_t j = i + 1; j < 4; ++j, ++pair) {
        uint8_t k = 0;
        while (k == i || k == j)
          ++k;
        uint8_t l = k + 1;
        while (l == i || l == j)
          ++l;
        const float same = sign[i] * sign[j];
        const uint8_t flipped =
            fabsf(direction[j]) < fabsf(direction[i]) ? j : i;
        const uint8_t neighbor = static_cast<uint8_t>(
            star ^ (flipped == 0 ? 7u : 1u << (flipped - 1)));
        add(star, i, j, k, l, -same, 2 * pair + (same < 0.0f));
        add(neighbor, i, j, k, l, same, 2 * pair + (same > 0.0f));
      }
    for (uint8_t family = 0; family < STREAM_COUNT; ++family) {
      if (owner_of[family] == UNOWNED) {
        live[family] = false;
        continue;
      }
      float position = local_origin[0] + near * local_direction[0];
      float speed = local_direction[0];
      for (int k = 1; k < 4; ++k) {
        const float p = local_origin[k] + near * local_direction[k];
        const bool negative = family & (1u << (k - 1));
        position += negative ? -p : p;
        speed += negative ? -local_direction[k] : local_direction[k];
      }
      start(family, 0.5f * position, 0.5f * speed, near);
    }
  }

  /**
   * @brief Coverage of the nearest owned strut at distance t on a stream.
   * @param index Stream index.
   * @param t Ray parameter.
   * @param feature Receives the nearest strut's feature id.
   * @return Coverage in [0, 1]; 0 when no strut reaches the footprint.
   */
  __attribute__((always_inline)) float coverage(size_t index, float t,
                                                uint32_t &feature) const {
    const auto &owner = owners[owner_of[index]];
    math::Vec4 residual;
    int parity = 0;
    for (int k = 0; k < 4; ++k) {
      const float q = local_origin[k] + local_direction[k] * t;
      const float rounded = roundf(q);
      residual[k] = q - rounded;
      parity += static_cast<int>(rounded);
    }
    const bool odd = parity & 1;
    float numerator = INFINITY;
    float denominator = 1.0f;
    for (size_t c = 0; c < owner.count; ++c) {
      const auto &strut = owner.classes[c];
      float across = residual[strut.i] - strut.sign * residual[strut.j];
      float rk = residual[strut.k];
      float rl = residual[strut.l];
      const float step = across > 0.5f ? 1.0f : across < -0.5f ? -1.0f : 0.0f;
      across -= step;
      // The nearest vertex on this class's coset flips the cheapest coordinate.
      if ((step != 0.0f) != odd) {
        const float cost = 0.5f - fabsf(across);
        const float cost_k = 1.0f - 2.0f * fabsf(rk);
        const float cost_l = 1.0f - 2.0f * fabsf(rl);
        if (cost <= cost_k && cost <= cost_l)
          across -= copysignf(1.0f, across);
        else if (cost_k <= cost_l)
          rk -= copysignf(1.0f, rk);
        else
          rl -= copysignf(1.0f, rl);
      }
      const float offset2 = 0.5f * across * across + rk * rk + rl * rl;
      const float dot =
          across * strut.transverse + rk * strut.dk + rl * strut.dl;
      const float n = fmaxf(0.0f, offset2 * strut.denominator - dot * dot);
      if (n * denominator < numerator * strut.denominator) {
        numerator = n;
        denominator = strut.denominator;
        feature = strut.feature;
      }
    }
    const float support = (wire_radius + .5f * footprint.at(t)) * inverse_scale;
    if (numerator > support * support * denominator)
      return 0.0f;
    return coverage_of(t, scale * sqrtf(numerator / denominator) - wire_radius);
  }

  /** @brief Coverage layer at a stream's next crossing.
   * @param index Stream index.
   * @return Contribution with the nearest owned strut's coverage. */
  Raycast::Contribution candidate(size_t index) const {
    Raycast::Contribution result;
    result.t = next[index];
    result.coverage = coverage(index, result.t, result.feature);
    return result;
  }
};

} // namespace SDF
