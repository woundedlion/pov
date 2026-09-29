/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/ray/camera.h"
#include "render/ray/shade.h"

namespace SDF {

HS_O3_BEGIN

/** @brief Cubic basis with XY shear and reciprocal X/Z stretch. */
struct AffineLattice {
  float cell_size = 1;
  float shear = .55f;
  float stretch = 1.4f;

  bool valid() const {
    return Raycast::finite(cell_size) && cell_size > 0 &&
           Raycast::finite(shear) && fabsf(shear) <= 1 &&
           Raycast::finite(stretch) && stretch >= .5f && stretch <= 2;
  }
  math::Vec4 point(const math::Vec4 &p) const {
    return {{cell_size * (stretch * p[0] + shear * p[1]), cell_size * p[1],
             cell_size * p[2] / stretch, cell_size * p[3]}};
  }
  math::Vec4 inverse(const math::Vec4 &p) const {
    return {{(p[0] - shear * p[1]) / (cell_size * stretch), p[1] / cell_size,
             p[2] * stretch / cell_size, p[3] / cell_size}};
  }
  /** @brief Wraps a camera center by exact lattice translation vectors. */
  math::Vec4 wrap(const math::Vec4 &p) const {
    auto q = inverse(p);
    for (int i = 0; i < 4; ++i)
      q[i] -= floorf(q[i]);
    return point(q);
  }
};

struct AffineLatticeEvents {
  static constexpr size_t STREAM_COUNT = 4;
  static constexpr size_t GROUP_CAPACITY = 1;
  AffineLattice geometry;
  math::Vec4 origin, direction, ambient_direction;
  std::array<float, 4> next{}, step{};
  std::array<int, 4> owner{};
  struct Metric {
    int first = 0;
    int second = 0;
    float aa = 0;
    float ab = 0;
    float bb = 0;
    bool active = false;
  };
  std::array<Metric, 4> metrics{};
  int dimensions;
  float radius;
  Raycast::Footprint footprint;

  AffineLatticeEvents(const Raycast::PreparedCamera &camera,
                      const math::Vector &view, AffineLattice geometry,
                      float radius, Raycast::Footprint footprint)
      : geometry(geometry),
        dimensions(camera.domain == Raycast::SamplingDomain::SLICE_4D ? 4 : 3),
        radius(radius), footprint(footprint) {
    ambient_direction = camera.embedding.apply({{view.x, view.y, view.z, 0}});
    direction = geometry.inverse(ambient_direction);
    origin = geometry.inverse(camera.point4(view * camera.radial_start));
    for (int axis = 0; axis < 4; ++axis) {
      next[axis] = INFINITY;
      owner[axis] = -1;
      if (axis >= dimensions)
        continue;
      float speed = 0;
      for (int other = 0; other < dimensions; ++other)
        if (other != axis && fabsf(direction[other]) > speed) {
          speed = fabsf(direction[other]);
          owner[axis] = other;
        }
      if (direction[axis] == 0)
        continue;
      const float P = origin[axis] + camera.interval.near * direction[axis];
      const float PLANE = direction[axis] > 0 ? floorf(P) + 1 : ceilf(P) - 1;
      next[axis] = camera.interval.near + (PLANE - P) / direction[axis];
      step[axis] = 1 / fabsf(direction[axis]);
    }
    for (int free = 0; free < dimensions; ++free) {
      if (owner[free] < 0)
        continue;
      auto &metric = metrics[free];
      math::Vec4 axis{};
      axis[free] = 1;
      const auto U = geometry.point(axis);
      float uu = 0, ud = 0;
      for (int k = 0; k < dimensions; ++k) {
        uu += U[k] * U[k];
        ud += U[k] * ambient_direction[k];
      }
      if (uu - ud * ud <= 1e-12f)
        continue;
      math::Vec4 transverse{};
      float length2 = 0;
      for (int k = 0; k < dimensions; ++k) {
        transverse[k] = ambient_direction[k] - U[k] * (ud / uu);
        length2 += transverse[k] * transverse[k];
      }
      std::array<math::Vec4, 2> projected{};
      int count = 0;
      for (int other = 0; other < dimensions; ++other) {
        if (other == free || other == owner[free])
          continue;
        (count == 0 ? metric.first : metric.second) = other;
        math::Vec4 basis{};
        basis[other] = 1;
        const auto B = geometry.point(basis);
        float bu = 0, bt = 0;
        for (int k = 0; k < dimensions; ++k) {
          bu += B[k] * U[k];
          bt += B[k] * transverse[k];
        }
        for (int k = 0; k < dimensions; ++k)
          projected[count][k] =
              B[k] - U[k] * (bu / uu) - transverse[k] * (bt / length2);
        ++count;
      }
      for (int k = 0; k < dimensions; ++k) {
        metric.aa += projected[0][k] * projected[0][k];
        metric.ab += projected[0][k] * projected[1][k];
        metric.bb += projected[1][k] * projected[1][k];
      }
      metric.active = true;
    }
  }
  bool active(size_t axis) const { return Raycast::finite(next[axis]); }
  float distance(size_t axis) const { return next[axis]; }
  void advance(size_t axis) { next[axis] += step[axis]; }

  Raycast::Contribution candidate(size_t plane) const {
    const float T = next[plane];
    float best = INFINITY;
    uint32_t feature = 0;
    for (int free = 0; free < dimensions; ++free) {
      if (owner[free] != static_cast<int>(plane))
        continue;
      const auto &metric = metrics[free];
      if (!metric.active)
        continue;
      const float P = origin[metric.first] + T * direction[metric.first];
      const float A = P - roundf(P);
      if (dimensions == 3) {
        const float SQUARED = metric.aa * A * A;
        if (SQUARED < best) {
          best = SQUARED;
          feature = free;
        }
      } else {
        const float Q = origin[metric.second] + T * direction[metric.second];
        const float B = Q - roundf(Q);
        for (int j = -1; j <= 1; ++j) {
          const float Y = B + j;
          for (int i = -1; i <= 1; ++i) {
            const float X = A + i;
            const float SQUARED =
                std::max(0.0f, metric.aa * X * X + 2 * metric.ab * X * Y +
                                   metric.bb * Y * Y);
            if (SQUARED < best) {
              best = SQUARED;
              feature = free;
            }
          }
        }
      }
    }
    const float WIDTH = footprint.at(T);
    const float DISTANCE = sqrtf(best) - radius;
    Raycast::Contribution hit;
    hit.t = T;
    hit.coverage = WIDTH > 0 ? std::clamp(.5f - DISTANCE / WIDTH, 0.0f, 1.0f)
                             : (DISTANCE <= 0 ? 1.0f : 0.0f);
    hit.feature = feature;
    return hit;
  }
};

struct AffineLatticeRenderer {
  HS_HOT_FLASH_MEMBER static Raycast::ShadedTrace
  shade(const Raycast::PreparedCamera &camera, const math::Vector &direction,
        float cell_size, float wire_radius, Raycast::Footprint footprint,
        const Raycast::TraceLimits &limits,
        const Raycast::Appearance &appearance, float shear = .55f,
        float stretch = 1.4f) {
    const AffineLattice GEOMETRY{cell_size, shear, stretch};
    if (!camera.valid() || !camera.ray(direction).valid() ||
        !GEOMETRY.valid() || !Raycast::finite(wire_radius) ||
        wire_radius <= 0 || !appearance.palette ||
        !Raycast::finite(footprint.angular_radius) ||
        footprint.angular_radius < 0 ||
        !Raycast::finite(footprint.radial_start) ||
        footprint.radial_start < 0) {
      Raycast::ShadedTrace result{};
      result.trace.status = Raycast::TraceStatus::INVALID_QUERY;
      return result;
    }
    AffineLatticeEvents events(camera, direction, GEOMETRY, wire_radius,
                               footprint);
    return Raycast::shade_events(events, camera.interval, limits, appearance);
  }
};

HS_O3_END

inline Raycast::ShadedTrace shade_affine_lattice(
    const Raycast::PreparedCamera &camera, const math::Vector &direction,
    float cell_size, float wire_radius, Raycast::Footprint footprint,
    const Raycast::TraceLimits &limits, const Raycast::Appearance &appearance,
    float shear = .55f, float stretch = 1.4f) {
  return AffineLatticeRenderer::shade(camera, direction, cell_size, wire_radius,
                                      footprint, limits, appearance, shear,
                                      stretch);
}

} // namespace SDF
