/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/ray/camera.h"
#include "render/ray/shade.h"

namespace SDF {

/** @brief Analytic intersections with disjoint periodic spherical sheets. */
struct PeriodicShells {
  float cell_size = 1;
  float shell_radius = .3f;

  bool valid() const {
    return Raycast::finite(cell_size) && cell_size > 0 &&
           Raycast::finite(shell_radius) && shell_radius > 0 &&
           shell_radius <= .480001f;
  }
  struct Intersections {
    float near = 0;
    float far = 0;
    bool hit = false;
  };
  /** @brief Roots in world-distance units for a cell-centered ambient ray. */
  Intersections intersect(const math::Vec4 &origin, const math::Vec4 &direction,
                          int dimensions) const {
    float a = 0, b = 0,
          c = -shell_radius * shell_radius * cell_size * cell_size;
    for (int k = 0; k < dimensions; ++k) {
      const float O = origin[k], D = direction[k];
      a += D * D;
      b += O * D;
      c += O * O;
    }
    const float DISCRIMINANT = b * b - a * c;
    if (DISCRIMINANT < 0 || a <= 0)
      return {};
    const float ROOT = sqrtf(DISCRIMINANT);
    const float Q = -b - copysignf(ROOT, b);
    const float FIRST = Q == 0 ? 0 : Q / a;
    const float SECOND = Q == 0 ? 0 : c / Q;
    return {std::min(FIRST, SECOND), std::max(FIRST, SECOND), true};
  }
};

/** @brief Validated frame constants; shade with the same camera used to prepare. */
struct PreparedPeriodicShells {
  PeriodicShells geometry;
  Raycast::Footprint footprint;
  float radius_squared = 0;
  bool valid = false;
};

inline PreparedPeriodicShells
prepare_periodic_shells(const Raycast::PreparedCamera &camera, float cell_size,
                        float shell_radius, Raycast::Footprint footprint) {
  PreparedPeriodicShells prepared;
  prepared.geometry = {cell_size, shell_radius};
  prepared.footprint = footprint;
  prepared.valid = camera.valid() && prepared.geometry.valid() &&
                   Raycast::finite(footprint.angular_radius) &&
                   footprint.angular_radius >= 0 &&
                   Raycast::finite(footprint.radial_start) &&
                   footprint.radial_start >= 0;
  if (prepared.valid) {
    prepared.radius_squared =
        shell_radius * shell_radius * cell_size * cell_size;
  }
  return prepared;
}

template <int DIMENSIONS>
HS_HOT_FLASH_MEMBER Raycast::ShadedTrace shade_periodic_shells_dimension(
    const PreparedPeriodicShells &prepared,
    const Raycast::PreparedCamera &camera, const math::Vector &direction,
    const Raycast::TraceLimits &limits, const Raycast::Appearance &appearance) {
  static_assert(DIMENSIONS == 3 || DIMENSIONS == 4);
  Raycast::ShadedTrace result{};
  const float cell_size = prepared.geometry.cell_size;
  const float shell_radius = prepared.geometry.shell_radius;
  const auto footprint = prepared.footprint;
  if (!prepared.valid || !Raycast::finite(direction) ||
      fabsf(math::dot(direction, direction) - 1.0f) >= 1e-4f ||
      !appearance.palette) {
    result.trace.status = Raycast::TraceStatus::INVALID_QUERY;
    return result;
  }
  const auto DIRECTION =
      camera.embedding.apply({{direction.x, direction.y, direction.z, 0}});
  const auto ORIGIN = camera.point4(direction * camera.radial_start);
  float a = 0;
  for (int k = 0; k < DIMENSIONS; ++k)
    a += DIRECTION[k] * DIRECTION[k];
  const float INVERSE_A = 1 / a;
  math::Vec4 cell{}, boundary{}, step{};
  for (int k = 0; k < DIMENSIONS; ++k) {
    const float P = ORIGIN[k] + camera.interval.near * DIRECTION[k];
    cell[k] = floorf(P / cell_size + .5f);
    if (DIRECTION[k] < 0 && P == (cell[k] - .5f) * cell_size)
      --cell[k];
    boundary[k] = DIRECTION[k] == 0
                      ? INFINITY
                      : ((cell[k] + copysignf(.5f, DIRECTION[k])) * cell_size -
                         ORIGIN[k]) /
                            DIRECTION[k];
    step[k] = DIRECTION[k] == 0 ? INFINITY : cell_size / fabsf(DIRECTION[k]);
  }
  LayerComposite composite;
  float start = camera.interval.near;
  while (start < camera.interval.far) {
    if (result.trace.counters.steps >= limits.max_steps) {
      result.trace.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
      break;
    }
    ++result.trace.counters.steps;
    float end = camera.interval.far;
    math::Vec4 local{};
    for (int k = 0; k < DIMENSIONS; ++k) {
      end = std::min(end, boundary[k]);
      local[k] = ORIGIN[k] - cell[k] * cell_size;
    }
    float b = 0, c = -prepared.radius_squared;
    for (int k = 0; k < DIMENSIONS; ++k) {
      const float O = local[k];
      b += O * DIRECTION[k];
      c += O * O;
    }
    const float DISCRIMINANT = b * b - a * c;
    PeriodicShells::Intersections roots;
    if (DISCRIMINANT >= 0 && a > 0) {
      const float Q = -b - copysignf(sqrtf(DISCRIMINANT), b);
      const float FIRST = Q == 0 ? 0 : Q * INVERSE_A;
      const float SECOND = Q == 0 ? 0 : c / Q;
      roots = {std::min(FIRST, SECOND), std::max(FIRST, SECOND), true};
    }
    const bool VERIFIED = roots.hit;
    float filtered_coverage = 0;
    if (!roots.hit && footprint.angular_radius > 0) {
      const float T = -b * INVERSE_A;
      if (T >= start && T <= end) {
        float squared = 0;
        math::Vec4 gradient{};
        for (int k = 0; k < DIMENSIONS; ++k) {
          const float P = local[k] + T * DIRECTION[k];
          squared += P * P;
          gradient[k] = P;
        }
        const float WIDTH = footprint.at(T);
        if (squared > 0 && WIDTH > 0) {
          const float LENGTH = sqrtf(squared);
          float distance = LENGTH - shell_radius * cell_size;
          if constexpr (DIMENSIONS == 4) {
            float projected_squared = 0;
            for (int j = 0; j < 3; ++j) {
              float component = 0;
              for (int k = 0; k < DIMENSIONS; ++k)
                component += camera.embedding.m[k][j] * gradient[k];
              projected_squared += component * component;
            }
            distance = projected_squared > 0
                           ? distance * LENGTH / sqrtf(projected_squared)
                           : INFINITY;
          }
          filtered_coverage = std::clamp(.5f - distance / WIDTH, 0.0f, 1.0f);
          if (filtered_coverage > 0)
            roots = {T, T, true};
        }
      }
    }
    if (roots.hit) {
      const float TIMES[2] = {roots.near, roots.far};
      for (int root = 0; root < 2; ++root) {
        const float T = TIMES[root];
        if (T < start || T > end || (root == 1 && T == TIMES[0]))
          continue;
        if (result.trace.counters.candidates >= limits.max_candidates ||
            result.trace.counters.layers >= limits.max_layers) {
          result.trace.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
          result.color = composite.finish();
          return result;
        }
        ++result.trace.counters.candidates;
        ++result.trace.counters.layers;
        Raycast::Contribution hit;
        hit.t = T;
        hit.verified = VERIFIED;
        math::Vec4 gradient{};
        for (int k = 0; k < DIMENSIONS; ++k)
          gradient[k] = local[k] + T * DIRECTION[k];
        float components[3]{};
        for (int j = 0; j < 3; ++j)
          for (int k = 0; k < DIMENSIONS; ++k)
            components[j] += camera.embedding.m[k][j] * gradient[k];
        const math::Vector PROJECTED{components[0], components[1],
                                     components[2]};
        const float LENGTH2 = math::dot(PROJECTED, PROJECTED);
        const float INVERSE_LENGTH =
            LENGTH2 > 1e-20f && Raycast::finite(LENGTH2) ? 1 / sqrtf(LENGTH2)
                                                         : 0;
        hit.has_normal = INVERSE_LENGTH > 0;
        hit.normal = PROJECTED * INVERSE_LENGTH;
        const float WIDTH = footprint.at(T);
        if (!VERIFIED) {
          hit.coverage = filtered_coverage;
        } else if (WIDTH > 0) {
          const float INCIDENCE = b + T * a;
          const float DEPTH = .25f * (roots.far - roots.near) *
                              fabsf(INCIDENCE) * INVERSE_LENGTH;
          hit.coverage = std::clamp(.5f + DEPTH / WIDTH, 0.0f, 1.0f);
        }
        result.trace.has_surface = true;
        result.trace.contribution = hit;
        composite.add(appearance.color(T),
                      hit.coverage * appearance.opacity(T));
        if (composite.saturated()) {
          result.trace.status = Raycast::TraceStatus::SATURATED;
          result.color = composite.finish();
          return result;
        }
      }
    }
    if (end >= camera.interval.far)
      break;
    if (!(end > start)) {
      result.trace.status = Raycast::TraceStatus::UNRESOLVED;
      break;
    }
    for (int k = 0; k < DIMENSIONS; ++k)
      if (boundary[k] <= end) {
        cell[k] += copysignf(1.0f, DIRECTION[k]);
        boundary[k] += step[k];
      }
    start = end;
  }
  result.color = composite.finish();
  return result;
}

inline Raycast::ShadedTrace shade_periodic_shells(
    const PreparedPeriodicShells &prepared,
    const Raycast::PreparedCamera &camera, const math::Vector &direction,
    const Raycast::TraceLimits &limits, const Raycast::Appearance &appearance) {
  return camera.domain == Raycast::SamplingDomain::SLICE_4D
             ? shade_periodic_shells_dimension<4>(prepared, camera, direction,
                                                  limits, appearance)
             : shade_periodic_shells_dimension<3>(prepared, camera, direction,
                                                  limits, appearance);
}

inline Raycast::ShadedTrace shade_periodic_shells(
    const Raycast::PreparedCamera &camera, const math::Vector &direction,
    float cell_size, float shell_radius, Raycast::Footprint footprint,
    const Raycast::TraceLimits &limits, const Raycast::Appearance &appearance) {
  return shade_periodic_shells(
      prepare_periodic_shells(camera, cell_size, shell_radius, footprint),
      camera, direction, limits, appearance);
}

} // namespace SDF
