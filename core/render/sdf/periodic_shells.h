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
  float inverse_cell = 0;
  /**
   * @brief Whether trace_periodic_shells_3d() may serve this frame: a 3D
   *        camera whose widest filtered sphere stays inside one lattice layer's
   *        rounding cell for every ray direction.
   */
  bool single_owner = false;
  /**
   * @brief Whether trace_periodic_shells_march() may serve this frame: its
   *        widest filtered sphere stays within half a cell.
   */
  bool march = false;
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
    prepared.inverse_cell = 1 / cell_size;
    // A contributing sphere lies within shell_radius plus half the footprint
    // at the far distance of the ray. Seen from the ray's crossing of the
    // sphere's layer across the dominant axis, that clearance grows by at most
    // sqrt(3); below half a cell the crossing rounds to the sphere's center.
    const float WIDEST = shell_radius + .5f *
                                            footprint.at(camera.interval.far) *
                                            prepared.inverse_cell;
    prepared.single_owner =
        camera.domain == Raycast::SamplingDomain::SPATIAL_3D &&
        1.7320508f * WIDEST < .499f;
    prepared.march = WIDEST < .499f;
  }
  return prepared;
}

/** @brief Premultiplied shell color and how its traversal ended. */
struct ShellSample {
  Pixel color;
  Raycast::TraceStatus status = Raycast::TraceStatus::RANGE_COMPLETE;
};

/**
 * @brief 3D shell composite by a march over the lattice layers across the
 *        ray's dominant axis.
 * @details Requires PreparedPeriodicShells::single_owner. Each layer crossing
 * then has one candidate sphere, centered at the crossing's rounded lattice
 * point, and a sphere's roots and filtered closest approach both lie inside its
 * own cell, so a crossing whose ray-to-center clearance exceeds the widest
 * contribution is rejected before any root solve. Matches
 * shade_periodic_shells_dimension<3>, with the step budget counting layers.
 * @return The premultiplied composite.
 */
__attribute__((always_inline)) inline ShellSample trace_periodic_shells_3d(
    const PreparedPeriodicShells &prepared,
    const Raycast::PreparedCamera &camera, const math::Vector &direction,
    const Raycast::TraceLimits &limits, const Raycast::Appearance &appearance) {
  const auto &E = camera.embedding.m;
  const float D[3] = {
      E[0][0] * direction.x + E[0][1] * direction.y + E[0][2] * direction.z,
      E[1][0] * direction.x + E[1][1] * direction.y + E[1][2] * direction.z,
      E[2][0] * direction.x + E[2][1] * direction.y + E[2][2] * direction.z};
  const float CELL = prepared.geometry.cell_size;
  const float INVERSE_CELL = prepared.inverse_cell;
  const float RADIAL_START = camera.radial_start;
  int k = fabsf(D[1]) > fabsf(D[0]) ? 1 : 0;
  if (fabsf(D[2]) > fabsf(D[k]))
    k = 2;
  const int I = k == 0 ? 1 : 0;
  const int J = k == 2 ? 1 : 2;
  const float A = D[0] * D[0] + D[1] * D[1] + D[2] * D[2];
  // A is 1 to rounding, so one Newton step from 1 is its reciprocal.
  const float INVERSE_A = 2 - A;
  const float ORIGIN_K =
      (camera.center[k] + D[k] * RADIAL_START) * INVERSE_CELL;
  const float INVERSE_DK = 1 / D[k];
  const float STEP = fabsf(INVERSE_DK);
  const float NEAR = camera.interval.near;
  const float FAR = camera.interval.far;
  const float LAST = FAR * INVERSE_CELL + .5f * STEP;
  float tau =
      (rintf(ORIGIN_K + NEAR * INVERSE_CELL * D[k]) - ORIGIN_K) * INVERSE_DK;
  float y =
      (camera.center[I] + D[I] * RADIAL_START) * INVERSE_CELL + tau * D[I];
  float z =
      (camera.center[J] + D[J] * RADIAL_START) * INVERSE_CELL + tau * D[J];
  const float DY = D[I] * STEP;
  const float DZ = D[J] * STEP;
  // Clearance bound in cells at the layer, padded by the farthest a closest
  // approach can sit from its layer crossing.
  const float HALF_RATE = .5f * prepared.footprint.angular_radius;
  float reach = prepared.geometry.shell_radius + 1e-4f +
                HALF_RATE * (RADIAL_START * INVERSE_CELL + tau + .75f);
  const float REACH_STEP = HALF_RATE * STEP;
  const float RADIUS = prepared.geometry.shell_radius * CELL;
  const float HALF_INVERSE_RADIUS = .5f / RADIUS;
  const bool FILTERED = prepared.footprint.angular_radius > 0;

  ShellSample result;
  LayerComposite composite;
  int layers = 0;
  const auto emit = [&](float t,
                        float coverage) __attribute__((always_inline)) {
    if (layers >= limits.max_candidates || layers >= limits.max_layers) {
      result.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
      return false;
    }
    ++layers;
    appearance.composite(composite, t, coverage);
    if (composite.saturated()) {
      result.status = Raycast::TraceStatus::SATURATED;
      return false;
    }
    return true;
  };
  // Counting the layers once keeps float compares out of the march; a layer
  // within rounding of LAST holds nothing nearer than FAR.
  const int LAYERS =
      tau <= LAST ? static_cast<int>((LAST - tau) * fabsf(D[k])) + 1 : 0;
  for (int step = 0; step < std::min(LAYERS, limits.max_steps); ++step) {
    const float EY = y - rintf(y);
    const float EZ = z - rintf(z);
    const float OFFSET2 = EY * EY + EZ * EZ;
    const float ALONG = EY * D[I] + EZ * D[J];
    if (OFFSET2 - ALONG * ALONG < reach * reach) {
      const float CROSSING = tau * CELL;
      const float B = ALONG * CELL;
      const float C = OFFSET2 * CELL * CELL - prepared.radius_squared;
      const float DISCRIMINANT = B * B - A * C;
      if (DISCRIMINANT >= 0) {
        const float ROOT = sqrtf(DISCRIMINANT);
        const float TIMES[2] = {CROSSING + (-B - ROOT) * INVERSE_A,
                                CROSSING + (-B + ROOT) * INVERSE_A};
        const float DEPTH = DISCRIMINANT * HALF_INVERSE_RADIUS;
        for (int root = 0; root < 2; ++root) {
          const float T = TIMES[root];
          if (T < NEAR || T > FAR || (root == 1 && T == TIMES[0]))
            continue;
          const float WIDTH = prepared.footprint.at(T);
          const float COVERAGE =
              WIDTH > 0 ? hs::clamp(.5f + DEPTH / WIDTH, 0.0f, 1.0f) : 1.0f;
          if (!emit(T, COVERAGE)) {
            result.color = composite.premultiplied();
            return result;
          }
        }
      } else if (FILTERED) {
        const float T = CROSSING - B * INVERSE_A;
        const float SQUARED = C + prepared.radius_squared - B * B * INVERSE_A;
        const float WIDTH = prepared.footprint.at(T);
        if (T >= NEAR && T <= FAR && SQUARED > 0 && WIDTH > 0) {
          const float COVERAGE =
              hs::clamp(.5f - (sqrtf(SQUARED) - RADIUS) / WIDTH, 0.0f, 1.0f);
          if (COVERAGE > 0 && !emit(T, COVERAGE)) {
            result.color = composite.premultiplied();
            return result;
          }
        }
      }
    }
    y += DY;
    z += DZ;
    tau += STEP;
    reach += REACH_STEP;
  }
  if (LAYERS > limits.max_steps)
    result.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
  result.color = composite.premultiplied();
  return result;
}

/**
 * @brief Shell composite by a march over the lattice layers across the
 *        ambient ray's dominant axis, one step to either side included.
 * @tparam DIMENSIONS 3 for a spatial camera, 4 for a slice.
 * @details Requires PreparedPeriodicShells::march for the camera's domain. A
 * contributing sphere lies within half a cell of the ambient ray, so its
 * closest approach and roots stay inside its own cell, and its center lies
 * within twice that clearance of the ray's crossing of its layer (sqrt(3)
 * times in 3D): the crossing's rounded lattice point, or a neighbor one step
 * across each coordinate whose rounding offset leaves it within reach. A
 * layer's candidates composite in distance order. Matches
 * shade_periodic_shells_dimension<DIMENSIONS>, with the step budget counting
 * layers.
 * @return The premultiplied composite.
 */
template <int DIMENSIONS>
__attribute__((always_inline)) inline ShellSample trace_periodic_shells_march(
    const PreparedPeriodicShells &prepared,
    const Raycast::PreparedCamera &camera, const math::Vector &direction,
    const Raycast::TraceLimits &limits, const Raycast::Appearance &appearance) {
  static_assert(DIMENSIONS == 3 || DIMENSIONS == 4);
  constexpr int OTHERS = DIMENSIONS - 1;
  const auto &E = camera.embedding.m;
  float D[DIMENSIONS];
  for (int axis = 0; axis < DIMENSIONS; ++axis)
    D[axis] = E[axis][0] * direction.x + E[axis][1] * direction.y +
              E[axis][2] * direction.z;
  const float CELL = prepared.geometry.cell_size;
  const float INVERSE_CELL = prepared.inverse_cell;
  const float RADIAL_START = camera.radial_start;
  int k = 0;
  for (int axis = 1; axis < DIMENSIONS; ++axis)
    if (fabsf(D[axis]) > fabsf(D[k]))
      k = axis;
  int other[OTHERS];
  for (int axis = 0, n = 0; axis < DIMENSIONS; ++axis)
    if (axis != k)
      other[n++] = axis;
  float A = 0;
  for (int axis = 0; axis < DIMENSIONS; ++axis)
    A += D[axis] * D[axis];
  // A is 1 to rounding, so one Newton step from 1 is its reciprocal.
  const float INVERSE_A = 2 - A;
  const float INVERSE_DK = 1 / D[k];
  const float STEP = fabsf(INVERSE_DK);
  const float NEAR = camera.interval.near;
  const float FAR = camera.interval.far;
  const float LAST = FAR * INVERSE_CELL + .5f * STEP;
  const float ORIGIN_K =
      (camera.center[k] + D[k] * RADIAL_START) * INVERSE_CELL;
  float tau =
      (rintf(ORIGIN_K + NEAR * INVERSE_CELL * D[k]) - ORIGIN_K) * INVERSE_DK;
  float position[OTHERS], rate[OTHERS], d[OTHERS];
  for (int n = 0; n < OTHERS; ++n) {
    const int AXIS = other[n];
    d[n] = D[AXIS];
    rate[n] = D[AXIS] * STEP;
    position[n] =
        (camera.center[AXIS] + D[AXIS] * RADIAL_START) * INVERSE_CELL +
        tau * D[AXIS];
  }
  // Clearance bound in cells at the layer, padded by the farthest a closest
  // approach can sit from its layer crossing; over |D_k| it bounds the
  // in-layer offset of a contributing center.
  const float HALF_RATE = .5f * prepared.footprint.angular_radius;
  float reach = prepared.geometry.shell_radius + 1e-4f +
                HALF_RATE * (RADIAL_START * INVERSE_CELL + tau + .9f);
  const float REACH_STEP = HALF_RATE * STEP;
  const float RADIUS = prepared.geometry.shell_radius * CELL;
  const bool FILTERED = prepared.footprint.angular_radius > 0;

  ShellSample result;
  LayerComposite composite;
  int layers = 0;
  struct Layer {
    float t, coverage;
  };
  std::array<Layer, 2 << OTHERS> pending;
  const auto projected_length2 =
      [&](const float (&gradient)[DIMENSIONS]) __attribute__((always_inline)) {
        float length2 = 0;
        for (int j = 0; j < 3; ++j) {
          float component = 0;
          for (int axis = 0; axis < DIMENSIONS; ++axis)
            component += E[axis][j] * gradient[axis];
          length2 += component * component;
        }
        return length2;
      };
  // Counting the layers once keeps float compares out of the march; a layer
  // within rounding of LAST holds nothing nearer than FAR.
  const int LAYERS =
      tau <= LAST ? static_cast<int>((LAST - tau) * fabsf(D[k])) + 1 : 0;
  for (int step = 0; step < std::min(LAYERS, limits.max_steps); ++step) {
    std::array<float, OTHERS> offset;
    float farthest = 0;
    for (int n = 0; n < OTHERS; ++n) {
      offset[n] = position[n] - rintf(position[n]);
      farthest = fmaxf(farthest, fabsf(offset[n]));
    }
    const float REACH2 = reach * reach;
    const float SPREAD = reach * STEP;
    int count = 0;
    const auto candidate =
        [&](const std::array<float, OTHERS> &e) __attribute__((always_inline)) {
          float offset2 = 0, along = 0;
          for (int n = 0; n < OTHERS; ++n) {
            offset2 += e[n] * e[n];
            along += e[n] * d[n];
          }
          if (!(offset2 - along * along < REACH2))
            return;
          const float CROSSING = tau * CELL;
          const float B = along * CELL;
          const float C = offset2 * CELL * CELL - prepared.radius_squared;
          const float DISCRIMINANT = B * B - A * C;
          float gradient[DIMENSIONS];
          gradient[k] = 0;
          for (int n = 0; n < OTHERS; ++n)
            gradient[other[n]] = e[n] * CELL;
          const auto at = [&](float u, float (&point)[DIMENSIONS])
                              __attribute__((always_inline)) {
                                for (int axis = 0; axis < DIMENSIONS; ++axis)
                                  point[axis] = gradient[axis] + u * D[axis];
                              };
          if (DISCRIMINANT >= 0) {
            const float ROOT = sqrtf(DISCRIMINANT);
            const float TIMES[2] = {CROSSING + (-B - ROOT) * INVERSE_A,
                                    CROSSING + (-B + ROOT) * INVERSE_A};
            for (int root = 0; root < 2; ++root) {
              const float T = TIMES[root];
              if (T < NEAR || T > FAR || (root == 1 && T == TIMES[0]))
                continue;
              const float WIDTH = prepared.footprint.at(T);
              float coverage = 1.0f;
              if (WIDTH > 0) {
                float point[DIMENSIONS];
                at(T - CROSSING, point);
                const float LENGTH2 = projected_length2(point);
                const float INVERSE_LENGTH =
                    LENGTH2 > 1e-20f && Raycast::finite(LENGTH2)
                        ? 1 / sqrtf(LENGTH2)
                        : 0;
                const float DEPTH =
                    .5f * DISCRIMINANT * INVERSE_A * INVERSE_LENGTH;
                coverage = hs::clamp(.5f + DEPTH / WIDTH, 0.0f, 1.0f);
              }
              pending[count++] = {T, coverage};
            }
          } else if (FILTERED) {
            const float U = -B * INVERSE_A;
            const float T = CROSSING + U;
            const float SQUARED =
                C + prepared.radius_squared - B * B * INVERSE_A;
            const float WIDTH = prepared.footprint.at(T);
            if (T >= NEAR && T <= FAR && SQUARED > 0 && WIDTH > 0) {
              const float LENGTH = sqrtf(SQUARED);
              float distance = LENGTH - RADIUS;
              if constexpr (DIMENSIONS == 4) {
                float point[DIMENSIONS];
                at(U, point);
                const float PROJECTED2 = projected_length2(point);
                distance = PROJECTED2 > 0
                               ? distance * LENGTH / sqrtf(PROJECTED2)
                               : INFINITY;
              }
              const float COVERAGE =
                  hs::clamp(.5f - distance / WIDTH, 0.0f, 1.0f);
              if (COVERAGE > 0)
                pending[count++] = {T, COVERAGE};
            }
          }
        };
    candidate(offset);
    if (farthest > 1 - SPREAD) {
      // Neighbors one step across each coordinate the reach can cross.
      std::array<float, OTHERS> across;
      unsigned crossing = 0;
      for (int n = 0; n < OTHERS; ++n) {
        across[n] = offset[n] - copysignf(1.0f, offset[n]);
        crossing |= static_cast<unsigned>(fabsf(across[n]) < SPREAD) << n;
      }
      for (unsigned mask = 1; mask < 1u << OTHERS; ++mask)
        if ((mask & crossing) == mask) {
          std::array<float, OTHERS> neighbor;
          for (int n = 0; n < OTHERS; ++n)
            neighbor[n] = mask & (1u << n) ? across[n] : offset[n];
          candidate(neighbor);
        }
    }
    for (int layer = 1; layer < count; ++layer)
      for (int slot = layer; slot > 0 && pending[slot - 1].t > pending[slot].t;
           --slot)
        std::swap(pending[slot - 1], pending[slot]);
    for (int layer = 0; layer < count; ++layer) {
      if (layers >= limits.max_candidates || layers >= limits.max_layers) {
        result.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
        result.color = composite.premultiplied();
        return result;
      }
      ++layers;
      appearance.composite(composite, pending[layer].t,
                           pending[layer].coverage);
      if (composite.saturated()) {
        result.status = Raycast::TraceStatus::SATURATED;
        result.color = composite.premultiplied();
        return result;
      }
    }
    for (int n = 0; n < OTHERS; ++n)
      position[n] += rate[n];
    tau += STEP;
    reach += REACH_STEP;
  }
  if (LAYERS > limits.max_steps)
    result.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
  result.color = composite.premultiplied();
  return result;
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
