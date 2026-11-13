/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file periodic_shells.h
 * @brief Periodic sphere shells and ray shading. */

#include "core/render/ray/camera.h"
#include "core/render/ray/shade.h"
#include "core/render/sdf/framework.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

namespace SDF {

/** @brief Analytic intersections with disjoint periodic spherical sheets. */
struct PeriodicShells {
  float cell_size = 1;      ///< World edge length of one lattice cell; > 0.
  float shell_radius = .3f; ///< Sphere radius in cells; in (0, 0.48].

  /** @brief True when both parameters are finite and in range.
   *  @return Whether the shells can be traced. */
  bool valid() const {
    return Raycast::finite(cell_size) && cell_size > 0 &&
           Raycast::finite(shell_radius) && shell_radius > 0 &&
           shell_radius <= .480001f;
  }
  /** @brief Ordered ray parameters of a ray's two sphere roots. */
  struct Intersections {
    float near = 0;   ///< Smaller root.
    float far = 0;    ///< Larger root.
    bool hit = false; ///< Whether the ray meets the sphere.
  };
  /**
   * @brief Roots in world-distance units for a cell-centered ambient ray.
   * @param origin Ray origin relative to the cell center, world units.
   * @param direction Ray direction.
   * @param dimensions Coordinates used, 3 or 4.
   * @return Ordered roots and whether the ray meets the sphere.
   */
  Intersections intersect(const math::Vec4 &origin, const math::Vec4 &direction,
                          int dimensions) const {
    float a = 0, b = 0,
          c = -shell_radius * shell_radius * cell_size * cell_size;
    for (int k = 0; k < dimensions; ++k) {
      const float o = origin[k], d = direction[k];
      a += d * d;
      b += o * d;
      c += o * o;
    }
    const float discriminant = b * b - a * c;
    if (discriminant < 0 || a <= 0)
      return {};
    const float root = sqrtf(discriminant);
    const float q = -b - copysignf(root, b);
    const float first = q == 0 ? 0 : q / a;
    const float second = q == 0 ? 0 : c / q;
    return {std::min(first, second), std::max(first, second), true};
  }
};

/** @brief Validated frame constants; shade with the same camera used to prepare. */
struct PreparedPeriodicShells {
  PeriodicShells geometry;      ///< Shell lattice parameters.
  Raycast::Footprint footprint; ///< Ray filter footprint.
  float radius_squared = 0;     ///< Squared sphere radius in world units.
  float inverse_cell = 0;       ///< Reciprocal of the cell size.
  /**
   * @brief Whether trace_periodic_shells_3d() may serve this frame: a 3D
   *        camera whose widest filtered sphere stays inside one lattice layer's
   *        rounding cell for every ray direction.
   */
  bool single_owner = false;
  /**
   * @brief Whether the narrow trace_periodic_shells_march() and
   *        shade_periodic_shells_dimension() may serve this frame: its widest
   *        filtered sphere stays within half a cell.
   */
  bool march = false;
  bool valid = false; ///< Whether camera, geometry, and footprint are valid.
};

/**
 * @brief Validates the shell parameters and derives per-frame constants.
 * @param camera Camera the frame is shaded with.
 * @param cell_size World edge length of one lattice cell.
 * @param shell_radius Sphere radius in cells.
 * @param footprint Ray filter footprint.
 * @return Prepared constants; derived fields are set only when `valid`.
 */
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
    const float widest = shell_radius + .5f *
                                            footprint.at(camera.interval.far) *
                                            prepared.inverse_cell;
    prepared.single_owner =
        camera.domain == Raycast::SamplingDomain::SPATIAL_3D &&
        1.7320508f * widest < .499f;
    prepared.march = widest < .499f;
  }
  return prepared;
}

/** @brief Premultiplied shell color and how its traversal ended. */
using ShellSample = Raycast::TraceSample;

/** @brief Layers trace_periodic_shells_march() sorts per ray; the caller owns it. */
struct ShellLayerStorage {
  static constexpr int CAPACITY = 16; ///< Maximum pending layers.
  /** @brief One shell layer awaiting compositing. */
  struct Layer {
    float t;                 ///< Ray distance of the layer.
    float coverage;          ///< Coverage of the layer.
    float brightness = 1.0f; ///< Surface color multiplier.
  };
  std::array<Layer, CAPACITY> pending; ///< One ray's layers to composite.

  /**
   * @brief Adds a layer nearer than @p cutoff; when full, drops the farthest
   *        of it and the held layers and lowers @p cutoff to that distance.
   * @param count Layers held; grows by one unless full.
   * @param cutoff Distance from which layers are dropped.
   * @param t Ray distance of the layer.
   * @param coverage Coverage of the layer.
   * @param brightness Surface color multiplier.
   */
  void retain(int &count, float &cutoff, float t, float coverage,
              float brightness = 1.0f) {
    if (t >= cutoff)
      return;
    if (count < CAPACITY) {
      pending[count++] = {t, coverage, brightness};
      return;
    }
    int last = 0;
    for (int slot = 1; slot < count; ++slot)
      if (pending[slot].t > pending[last].t)
        last = slot;
    cutoff = fmaxf(t, pending[last].t);
    if (t < pending[last].t)
      pending[last] = {t, coverage, brightness};
  }
};

/** @brief Composites one shell layer within the candidate and layer budgets.
 * @param t Ray distance.
 * @param coverage Shell coverage.
 * @param brightness Surface color multiplier.
 * @param limits Candidate and layer budgets.
 * @param appearance Layer appearance.
 * @param layers Composited layer count.
 * @param result Receives termination status.
 * @param composite Accumulated color.
 * @return Whether tracing may continue. */
__attribute__((always_inline)) inline bool
composite_shell_layer(float t, float coverage, float brightness,
                      const Raycast::TraceLimits &limits,
                      const Raycast::Appearance &appearance, int &layers,
                      ShellSample &result, LayerComposite &composite) {
  if (layers >= limits.max_candidates || layers >= limits.max_layers) {
    result.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
    return false;
  }
  ++layers;
  appearance.composite(composite, t, coverage, brightness);
  if (composite.saturated()) {
    result.status = Raycast::TraceStatus::SATURATED;
    return false;
  }
  return true;
}

/** @brief Emits the roots or filtered tangent of one single-owner shell.
 * @param prepared Prepared shell constants.
 * @param tau Layer crossing in cells.
 * @param cell World cell size.
 * @param along Along-ray cell offset.
 * @param offset2 Squared cell offset.
 * @param a Squared ray length.
 * @param inverse_a Reciprocal ray-length estimate.
 * @param half_inverse_radius Reciprocal diameter.
 * @param radius Sphere radius in world units.
 * @param filtered Whether the footprint is filtered.
 * @param near Near ray bound.
 * @param far Far ray bound.
 * @param emit Layer sink.
 * @return Whether tracing may continue. */
template <bool SHADED, typename Emit>
__attribute__((always_inline)) inline bool emit_single_owner_shell(
    const PreparedPeriodicShells &prepared, float tau, float cell, float along,
    float offset2, float a, float inverse_a, float half_inverse_radius,
    float radius, bool filtered, float near, float far, Emit &&emit) {
  const float crossing = tau * cell;
  const float b = along * cell;
  const float c = offset2 * cell * cell - prepared.radius_squared;
  const float discriminant = b * b - a * c;
  if (discriminant >= 0) {
    const float root_value = sqrtf(discriminant);
    const float incidence = fminf(root_value * (2 * half_inverse_radius), 1.0f);
    const float shade = SHADED ? .1f + .9f * incidence * incidence : 1.0f;
    const float times[2] = {crossing + (-b - root_value) * inverse_a,
                            crossing + (-b + root_value) * inverse_a};
    const float depth = discriminant * half_inverse_radius;
    for (int root = 0; root < 2; ++root) {
      const float t = times[root];
      if (t < near || t > far || (root == 1 && t == times[0]))
        continue;
      const float width = prepared.footprint.at(t);
      const float coverage =
          width > 0 ? footprint_coverage<true>(width, depth) : 1.0f;
      if (!emit(t, coverage, shade)) {
        return false;
      }
    }
  } else if (filtered) {
    const float t = crossing - b * inverse_a;
    const float squared = c + prepared.radius_squared - b * b * inverse_a;
    const float width = prepared.footprint.at(t);
    if (t >= near && t <= far && squared > 0 && width > 0) {
      const float coverage = footprint_coverage(width, sqrtf(squared) - radius);
      if (coverage > 0 && !emit(t, coverage, SHADED ? .1f : 1.0f)) {
        return false;
      }
    }
  }
  return true;
}

/** @brief Projects an ambient shell gradient into the view subspace.
 * @param camera Prepared embedding.
 * @param gradient Ambient gradient.
 * @return Squared projected length. */
template <int DIMENSIONS>
__attribute__((always_inline)) inline float
shell_projected_length2(const Raycast::PreparedCamera &camera,
                        const float (&gradient)[DIMENSIONS]) {
  float length2 = 0;
  for (int j = 0; j < 3; ++j) {
    float component = 0;
    for (int axis = 0; axis < DIMENSIONS; ++axis)
      component += camera.embedding.m[axis][j] * gradient[axis];
    length2 += component * component;
  }
  return length2;
}

/** @brief Enumerates lattice centers whose layer can reach the ray.
 * @param position Layer crossing coordinates.
 * @param offset Nearest-center residual.
 * @param spread Coordinate reach.
 * @param farthest Largest nearest-center residual.
 * @param candidate Candidate sink. */
template <bool WIDE, size_t AXES, typename Candidate>
__attribute__((always_inline)) inline void
for_each_shell_neighbor(const float (&position)[AXES],
                        const std::array<float, AXES> &offset, float spread,
                        float farthest, Candidate &&candidate) {
  constexpr int OTHERS = static_cast<int>(AXES);
  if constexpr (WIDE) {
    // Every lattice point of the layer within spread on each coordinate.
    std::array<float, OTHERS> e, first;
    int span[OTHERS], index[OTHERS];
    bool empty = false;
    for (int n = 0; n < OTHERS; ++n) {
      const float upper = floorf(position[n] + spread);
      first[n] = e[n] = position[n] - upper;
      span[n] = static_cast<int>(upper - ceilf(position[n] - spread)) + 1;
      index[n] = 0;
      empty |= span[n] <= 0;
    }
    for (int n = 0; !empty && n < OTHERS;) {
      candidate(e);
      for (n = 0; n < OTHERS; ++n) {
        if (++index[n] < span[n]) {
          e[n] = first[n] + static_cast<float>(index[n]);
          break;
        }
        index[n] = 0;
        e[n] = first[n];
      }
    }
  } else {
    candidate(offset);
    if (farthest > 1 - spread) {
      // Neighbors one step across each coordinate the reach can cross.
      std::array<float, OTHERS> across;
      unsigned crossing = 0;
      for (int n = 0; n < OTHERS; ++n) {
        across[n] = offset[n] - copysignf(1.0f, offset[n]);
        crossing |= static_cast<unsigned>(fabsf(across[n]) < spread) << n;
      }
      for (unsigned mask = 1; mask < 1u << OTHERS; ++mask)
        if ((mask & crossing) == mask) {
          std::array<float, OTHERS> neighbor;
          for (int n = 0; n < OTHERS; ++n)
            neighbor[n] = mask & (1u << n) ? across[n] : offset[n];
          candidate(neighbor);
        }
    }
  }
}

/** @brief Sorts retained shell contributions by ray distance.
 * @param storage Retained shell layers.
 * @param count Number of retained layers. */
__attribute__((always_inline)) inline void
sort_shell_layers(ShellLayerStorage &storage, int count) {
  for (int layer = 1; layer < count; ++layer)
    for (int slot = layer;
         slot > 0 && storage.pending[slot - 1].t > storage.pending[slot].t;
         --slot)
      std::swap(storage.pending[slot - 1], storage.pending[slot]);
}

/** @brief Ray constants for one periodic-shell layer march. */
template <int DIMENSIONS> struct ShellMarchRay {
  const Raycast::PreparedCamera &camera; ///< Prepared view embedding.
  const float (&d_value)[DIMENSIONS];    ///< Ambient ray direction.
  const int (&other)[DIMENSIONS - 1];    ///< Transverse axis indices.
  const float (&d)[DIMENSIONS - 1];      ///< Transverse ray direction.
  int k;                                 ///< Dominant axis.
  float cell;                            ///< World cell size.
  float a;                               ///< Squared ray length.
  float inverse_a;                       ///< Reciprocal ray-length estimate.
  float near;                            ///< Near ray bound.
  float far;                             ///< Far ray bound.
  float inverse_radius; ///< Reciprocal sphere radius for incidence.
  float radius;         ///< Sphere radius in world units.
  bool filtered;        ///< Whether the footprint is filtered.
};

/** @brief Retains roots and filtered tangencies of one shell candidate.
 * @param prepared Prepared shell constants.
 * @param ray Layer-march ray constants.
 * @param e Candidate center residual.
 * @param tau Dominant-axis crossing.
 * @param reach2 Squared clearance bound.
 * @param storage Retained contributions.
 * @param count Number of retained contributions.
 * @param cutoff Overflow cutoff distance. */
template <int DIMENSIONS, bool WIDE, bool SHADED>
__attribute__((always_inline)) inline void
retain_shell_candidate(const PreparedPeriodicShells &prepared,
                       const ShellMarchRay<DIMENSIONS> &ray,
                       const std::array<float, DIMENSIONS - 1> &e, float tau,
                       float reach2, ShellLayerStorage &storage, int &count,
                       float &cutoff) {
  constexpr int OTHERS = DIMENSIONS - 1;
  auto &pending = storage.pending;
  float offset2 = 0, along = 0;
  for (int n = 0; n < OTHERS; ++n) {
    offset2 += e[n] * e[n];
    along += e[n] * ray.d[n];
  }
  if (!(offset2 - along * along < reach2))
    return;
  const float crossing = tau * ray.cell;
  const float b = along * ray.cell;
  const float c = offset2 * ray.cell * ray.cell - prepared.radius_squared;
  const float discriminant = b * b - ray.a * c;
  float gradient[DIMENSIONS];
  gradient[ray.k] = 0;
  for (int n = 0; n < OTHERS; ++n)
    gradient[ray.other[n]] = e[n] * ray.cell;
  const auto at = [&](float u, float (&point)[DIMENSIONS])
                      __attribute__((always_inline)) {
                        for (int axis = 0; axis < DIMENSIONS; ++axis)
                          point[axis] = gradient[axis] + u * ray.d_value[axis];
                      };
  if (discriminant >= 0) {
    const float root_value = sqrtf(discriminant);
    const float incidence = fminf(root_value * ray.inverse_radius, 1.0f);
    const float shade = SHADED ? .1f + .9f * incidence * incidence : 1.0f;
    const float times[2] = {crossing + (-b - root_value) * ray.inverse_a,
                            crossing + (-b + root_value) * ray.inverse_a};
    for (int root = 0; root < 2; ++root) {
      const float t = times[root];
      if (t < ray.near || t > ray.far || (root == 1 && t == times[0]))
        continue;
      const float width = prepared.footprint.at(t);
      float coverage = 1.0f;
      if (width > 0) {
        float point[DIMENSIONS];
        at(t - crossing, point);
        const float length2 =
            shell_projected_length2<DIMENSIONS>(ray.camera, point);
        const float inverse_length =
            length2 > 1e-20f && Raycast::finite(length2) ? 1 / sqrtf(length2)
                                                         : 0;
        const float depth = .5f * discriminant * ray.inverse_a * inverse_length;
        coverage = footprint_coverage<true>(width, depth);
      }
      if constexpr (WIDE)
        storage.retain(count, cutoff, t, coverage, shade);
      else
        pending[count++] = {t, coverage, shade};
    }
  } else if (ray.filtered) {
    const float u = -b * ray.inverse_a;
    const float t = crossing + u;
    const float squared = c + prepared.radius_squared - b * b * ray.inverse_a;
    const float width = prepared.footprint.at(t);
    if (t >= ray.near && t <= ray.far && squared > 0 && width > 0) {
      const float length = sqrtf(squared);
      float distance = length - ray.radius;
      if constexpr (DIMENSIONS == 4) {
        float point[DIMENSIONS];
        at(u, point);
        const float projected2 =
            shell_projected_length2<DIMENSIONS>(ray.camera, point);
        distance =
            projected2 > 0 ? distance * length / sqrtf(projected2) : INFINITY;
      }
      const float coverage = footprint_coverage(width, distance);
      if (coverage > 0) {
        if constexpr (WIDE)
          storage.retain(count, cutoff, t, coverage, SHADED ? .1f : 1.0f);
        else
          pending[count++] = {t, coverage, SHADED ? .1f : 1.0f};
      }
    }
  }
}

/** @brief Admits a filtered tangent when the cell has no exact roots.
 * @param prepared Prepared shell constants.
 * @param camera Prepared view embedding.
 * @param local Cell-relative ray origin.
 * @param direction_value Ambient ray direction.
 * @param b Ray-origin dot direction.
 * @param inverse_a Reciprocal squared direction.
 * @param start Cell interval start.
 * @param end Cell interval end.
 * @param roots Receives an admitted tangent.
 * @return Filtered tangent coverage. */
template <int DIMENSIONS>
HS_O3_FN __attribute__((always_inline)) inline float
admit_shell_tangent(const PreparedPeriodicShells &prepared,
                    const Raycast::PreparedCamera &camera,
                    const math::Vec4 &local, const math::Vec4 &direction_value,
                    float b, float inverse_a, float start, float end,
                    PeriodicShells::Intersections &roots) {
  float filtered_coverage = 0;
  if (!roots.hit && prepared.footprint.angular_radius > 0) {
    const float t = -b * inverse_a;
    if (t >= start && t <= end) {
      float squared = 0;
      math::Vec4 gradient{};
      for (int k = 0; k < DIMENSIONS; ++k) {
        const float p = local[k] + t * direction_value[k];
        squared += p * p;
        gradient[k] = p;
      }
      const float width = prepared.footprint.at(t);
      if (squared > 0 && width > 0) {
        const float length = sqrtf(squared);
        float distance = length - prepared.geometry.shell_radius *
                                      prepared.geometry.cell_size;
        if constexpr (DIMENSIONS == 4) {
          float projected_squared = 0;
          for (int j = 0; j < 3; ++j) {
            float component = 0;
            for (int k = 0; k < DIMENSIONS; ++k)
              component += camera.embedding.m[k][j] * gradient[k];
            projected_squared += component * component;
          }
          distance = projected_squared > 0
                         ? distance * length / sqrtf(projected_squared)
                         : INFINITY;
        }
        filtered_coverage = footprint_coverage(width, distance);
        if (filtered_coverage > 0)
          roots = {t, t, true};
      }
    }
  }
  return filtered_coverage;
}

/** @brief Composites ordered roots of one traversed shell cell.
 * @param prepared Prepared shell constants.
 * @param camera Prepared view embedding.
 * @param limits Candidate and layer budgets.
 * @param appearance Layer appearance.
 * @param roots Cell intersections.
 * @param verified Whether roots are exact.
 * @param filtered_coverage Tangent coverage.
 * @param local Cell-relative ray origin.
 * @param direction_value Ambient ray direction.
 * @param a Squared direction length.
 * @param b Ray-origin dot direction.
 * @param start Cell interval start.
 * @param end Cell interval end.
 * @param result Receives trace details and termination.
 * @param composite Accumulated color.
 * @return Whether tracing may continue. */
template <int DIMENSIONS>
HS_O3_FN __attribute__((always_inline)) inline bool composite_shell_cell(
    const PreparedPeriodicShells &prepared,
    const Raycast::PreparedCamera &camera, const Raycast::TraceLimits &limits,
    const Raycast::Appearance &appearance,
    const PeriodicShells::Intersections &roots, bool verified,
    float filtered_coverage, const math::Vec4 &local,
    const math::Vec4 &direction_value, float a, float b, float start, float end,
    Raycast::ShadedTrace &result, LayerComposite &composite) {

  const float times[2] = {roots.near, roots.far};
  for (int root = 0; root < 2; ++root) {
    const float t = times[root];
    if (t < start || t > end || (root == 1 && t == times[0]))
      continue;
    if (result.trace.counters.candidates >= limits.max_candidates ||
        result.trace.counters.layers >= limits.max_layers) {
      result.trace.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
      result.color = composite.finish();
      return false;
    }
    ++result.trace.counters.candidates;
    ++result.trace.counters.layers;
    Raycast::Contribution hit;
    hit.t = t;
    hit.verified = verified;
    math::Vec4 gradient{};
    for (int k = 0; k < DIMENSIONS; ++k)
      gradient[k] = local[k] + t * direction_value[k];
    float inverse_length;
    hit.has_normal = camera.project_normal_dimension<DIMENSIONS>(
        gradient, hit.normal, inverse_length);
    const float width = prepared.footprint.at(t);
    if (!verified) {
      hit.coverage = filtered_coverage;
    } else if (width > 0) {
      const float incidence = b + t * a;
      const float depth =
          .25f * (roots.far - roots.near) * fabsf(incidence) * inverse_length;
      hit.coverage = footprint_coverage<true>(width, depth);
    }
    result.trace.contribution = hit;
    appearance.composite(composite, t, hit.coverage);
    if (composite.saturated()) {
      result.trace.status = Raycast::TraceStatus::SATURATED;
      result.color = composite.finish();
      return false;
    }
  }

  return true;
}

/**
 * @brief 3D shell composite by a march over the lattice layers across the
 *        ray's dominant axis.
 * @details Both max_candidates and max_layers cap composited layers.
 * Requires PreparedPeriodicShells::single_owner. Each layer crossing
 * then has one candidate sphere, centered at the crossing's rounded lattice
 * point, and a sphere's roots and filtered closest approach both lie inside its
 * own cell, so a crossing whose ray-to-center clearance exceeds the widest
 * contribution is rejected before any root solve. Matches
 * shade_periodic_shells_dimension<3>, with the step budget counting layers.
 * @tparam SHADED Darkens shell color by the squared sphere-chord incidence.
 * @param prepared Frame shell constants.
 * @param camera Prepared camera embedding.
 * @param direction Unit view direction.
 * @param limits Step and layer budgets.
 * @param appearance Depth appearance applied to each layer.
 * @return The premultiplied composite.
 */
template <bool SHADED = false>
__attribute__((always_inline)) inline ShellSample trace_periodic_shells_3d(
    const PreparedPeriodicShells &prepared,
    const Raycast::PreparedCamera &camera, const math::Vector &direction,
    const Raycast::TraceLimits &limits, const Raycast::Appearance &appearance) {
  const auto &e = camera.embedding.m;
  const float d[3] = {
      e[0][0] * direction.x + e[0][1] * direction.y + e[0][2] * direction.z,
      e[1][0] * direction.x + e[1][1] * direction.y + e[1][2] * direction.z,
      e[2][0] * direction.x + e[2][1] * direction.y + e[2][2] * direction.z};
  const float cell = prepared.geometry.cell_size;
  const float inverse_cell = prepared.inverse_cell;
  const float radial_start = camera.radial_start;
  int k = fabsf(d[1]) > fabsf(d[0]) ? 1 : 0;
  if (fabsf(d[2]) > fabsf(d[k]))
    k = 2;
  const int i = k == 0 ? 1 : 0;
  const int j = k == 2 ? 1 : 2;
  const float a = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  // a is 1 to rounding, so one Newton step from 1 is its reciprocal.
  const float inverse_a = 2 - a;
  const float origin_k =
      (camera.center[k] + d[k] * radial_start) * inverse_cell;
  const float inverse_dk = 1 / d[k];
  const float step_value = fabsf(inverse_dk);
  const float near = camera.interval.near;
  const float far = camera.interval.far;
  const float last = far * inverse_cell + .5f * step_value;
  float tau =
      (rintf(origin_k + near * inverse_cell * d[k]) - origin_k) * inverse_dk;
  float y =
      (camera.center[i] + d[i] * radial_start) * inverse_cell + tau * d[i];
  float z =
      (camera.center[j] + d[j] * radial_start) * inverse_cell + tau * d[j];
  const float dy = d[i] * step_value;
  const float dz = d[j] * step_value;
  // Clearance bound in cells at the layer, padded by the farthest a closest
  // approach can sit from its layer crossing.
  const float half_rate = .5f * prepared.footprint.angular_radius;
  float reach =
      prepared.geometry.shell_radius + 1e-4f +
      half_rate * (prepared.footprint.radial_start * inverse_cell + tau + .75f);
  const float reach_step = half_rate * step_value;
  const float radius = prepared.geometry.shell_radius * cell;
  const float half_inverse_radius = .5f / radius;
  const bool filtered = prepared.footprint.angular_radius > 0;

  ShellSample result;
  LayerComposite composite;
  int layers = 0;
  const auto emit = [&](float t, float coverage,
                        float brightness) __attribute__((always_inline)) {
    return composite_shell_layer(t, coverage, brightness, limits, appearance,
                                 layers, result, composite);
  };
  // a layer within rounding of last holds nothing nearer than far.
  const int layer_limit =
      tau <= last ? static_cast<int>((last - tau) * fabsf(d[k])) + 1 : 0;
  for (int step = 0; step < std::min(layer_limit, limits.max_steps); ++step) {
    const float ey = y - rintf(y);
    const float ez = z - rintf(z);
    const float offset2 = ey * ey + ez * ez;
    const float along = ey * d[i] + ez * d[j];
    if (offset2 - along * along < reach * reach) {
      if (!emit_single_owner_shell<SHADED>(prepared, tau, cell, along, offset2,
                                           a, inverse_a, half_inverse_radius,
                                           radius, filtered, near, far, emit)) {
        result.color = composite.premultiplied();
        return result;
      }
    }
    y += dy;
    z += dz;
    tau += step_value;
    reach += reach_step;
  }
  if (layer_limit > limits.max_steps)
    result.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
  result.color = composite.premultiplied();
  return result;
}

/**
 * @brief Shell composite by a march over the lattice layers across the
 *        ambient ray's dominant axis.
 * @tparam DIMENSIONS 3 for a spatial camera, 4 for a slice.
 * @tparam WIDE Whether filtered spheres may reach past half a cell.
 * @tparam SHADED Darkens shell color by the squared sphere-chord incidence.
 * @details Both max_candidates and max_layers cap composited layers.
 * Without WIDE, requires PreparedPeriodicShells::march: candidates
 * are the ray's layer crossing's rounded lattice point and its neighbors one
 * step to either side within reach, composited in distance order per layer,
 * matching shade_periodic_shells_dimension<DIMENSIONS>. WIDE serves any
 * footprint: each layer's candidates are every lattice point whose filtered
 * sphere can reach the ray, the march extends past both interval ends by the
 * widest reach, and retained contributions composite in distance order once
 * no later layer can precede them. When more than ShellLayerStorage::CAPACITY
 * contributions await compositing, the farthest is dropped, the ray is cut
 * at its distance, and the status is BUDGET_EXHAUSTED. The step budget
 * counts layers.
 * @param prepared Frame shell constants.
 * @param camera Prepared camera embedding.
 * @param direction Unit view direction.
 * @param limits Step and layer budgets.
 * @param appearance Depth appearance applied to each layer.
 * @param storage Per-ray layer sort scratch.
 * @return The premultiplied composite.
 */
template <int DIMENSIONS, bool WIDE = false, bool SHADED = false>
__attribute__((always_inline)) inline ShellSample trace_periodic_shells_march(
    const PreparedPeriodicShells &prepared,
    const Raycast::PreparedCamera &camera, const math::Vector &direction,
    const Raycast::TraceLimits &limits, const Raycast::Appearance &appearance,
    ShellLayerStorage &storage) {
  static_assert(DIMENSIONS == 3 || DIMENSIONS == 4);
  constexpr int OTHERS = DIMENSIONS - 1;
  static_assert((2 << OTHERS) <= ShellLayerStorage::CAPACITY);
  const auto &e_value = camera.embedding.m;
  float d_value[DIMENSIONS];
  for (int axis = 0; axis < DIMENSIONS; ++axis)
    d_value[axis] = e_value[axis][0] * direction.x +
                    e_value[axis][1] * direction.y +
                    e_value[axis][2] * direction.z;
  const float cell = prepared.geometry.cell_size;
  const float inverse_cell = prepared.inverse_cell;
  const float radial_start = camera.radial_start;
  int k = 0;
  for (int axis = 1; axis < DIMENSIONS; ++axis)
    if (fabsf(d_value[axis]) > fabsf(d_value[k]))
      k = axis;
  int other[OTHERS];
  for (int axis = 0, n = 0; axis < DIMENSIONS; ++axis)
    if (axis != k)
      other[n++] = axis;
  float a = 0;
  for (int axis = 0; axis < DIMENSIONS; ++axis)
    a += d_value[axis] * d_value[axis];
  // a is 1 to rounding, so one Newton step from 1 is its reciprocal.
  const float inverse_a = 2 - a;
  const float inverse_dk = 1 / d_value[k];
  const float step_value = fabsf(inverse_dk);
  const float near = camera.interval.near;
  const float far = camera.interval.far;
  // Bound in cells on a contributing center's distance from the ray.
  const float widest = WIDE ? prepared.geometry.shell_radius + 1e-4f +
                                  .5f * prepared.footprint.angular_radius *
                                      (prepared.footprint.radial_start + far) *
                                      inverse_cell
                            : 0.0f;
  const float last =
      far * inverse_cell + (WIDE ? fmaxf(.5f, widest) : .5f) * step_value;
  const float origin_k =
      (camera.center[k] + d_value[k] * radial_start) * inverse_cell;
  float tau = (rintf(origin_k + near * inverse_cell * d_value[k]) - origin_k) *
              inverse_dk;
  if constexpr (WIDE)
    tau -= fmaxf(0.0f, ceilf(widest - .5f)) * step_value;
  float position[OTHERS], rate[OTHERS], d[OTHERS];
  for (int n = 0; n < OTHERS; ++n) {
    const int axis = other[n];
    d[n] = d_value[axis];
    rate[n] = d_value[axis] * step_value;
    position[n] =
        (camera.center[axis] + d_value[axis] * radial_start) * inverse_cell +
        tau * d_value[axis];
  }
  // Clearance bound in cells at the layer, padded by the farthest a closest
  // approach can sit from its layer crossing; over |D_k| it bounds the
  // in-layer offset of a contributing center.
  const float half_rate = .5f * prepared.footprint.angular_radius;
  float reach = prepared.geometry.shell_radius + 1e-4f +
                half_rate * (prepared.footprint.radial_start * inverse_cell +
                             tau + (WIDE ? widest * step_value : .9f));
  const float reach_step = half_rate * step_value;
  const float radius = prepared.geometry.shell_radius * cell;
  const bool filtered = prepared.footprint.angular_radius > 0;
  const float inverse_radius = SHADED ? 1 / radius : 0.0f;

  const ShellMarchRay<DIMENSIONS> ray{
      camera, d_value, other,          d,      k,       cell, a, inverse_a,
      near,   far,     inverse_radius, radius, filtered};

  ShellSample result;
  LayerComposite composite;
  int layers = 0;
  auto &pending = storage.pending;
  int held = 0;
  float cutoff = INFINITY;
  // a layer within rounding of last holds nothing nearer than far.
  const int layer_limit =
      tau <= last ? static_cast<int>((last - tau) * fabsf(d_value[k])) + 1 : 0;
  for (int step = 0; step < std::min(layer_limit, limits.max_steps); ++step) {
    std::array<float, OTHERS> offset;
    float farthest = 0;
    for (int n = 0; n < OTHERS; ++n) {
      offset[n] = position[n] - rintf(position[n]);
      farthest = fmaxf(farthest, fabsf(offset[n]));
    }
    const float reach_value = WIDE ? fminf(reach, widest) : reach;
    const float reach2 = reach_value * reach_value;
    const float spread = reach_value * step_value;
    int count = held;
    const auto candidate =
        [&](const std::array<float, OTHERS> &e) __attribute__((always_inline)) {
          retain_shell_candidate<DIMENSIONS, WIDE, SHADED>(
              prepared, ray, e, tau, reach2, storage, count, cutoff);
        };
    for_each_shell_neighbor<WIDE>(position, offset, spread, farthest,
                                  candidate);
    sort_shell_layers(storage, count);
    // Later layers contribute no nearer than HORIZON.
    const float horizon =
        WIDE && step + 1 < std::min(layer_limit, limits.max_steps)
            ? (tau + (1 - widest) * step_value) * cell
            : INFINITY;
    int ready = count;
    if constexpr (WIDE)
      for (ready = 0; ready < count && pending[ready].t < horizon;)
        ++ready;
    for (int layer = 0; layer < ready; ++layer) {
      if (!composite_shell_layer(pending[layer].t, pending[layer].coverage,
                                 pending[layer].brightness, limits, appearance,
                                 layers, result, composite)) {
        result.color = composite.premultiplied();
        return result;
      }
    }
    if constexpr (WIDE) {
      held = count - ready;
      for (int slot = 0; slot < held; ++slot)
        pending[slot] = pending[ready + slot];
      if (cutoff <= horizon)
        break;
    }
    for (int n = 0; n < OTHERS; ++n)
      position[n] += rate[n];
    tau += step_value;
    reach += reach_step;
  }
  if (layer_limit > limits.max_steps || (WIDE && cutoff < INFINITY))
    result.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
  result.color = composite.premultiplied();
  return result;
}

/**
 * @brief Shades one ray by a cell walk through the 3D or 4D shell lattice.
 * @tparam DIMENSIONS Lattice dimension, 3 or 4.
 * @details Each traversed cell contributes only its own sphere, so the walk
 * requires PreparedPeriodicShells::march.
 * @param prepared Frame constants from prepare_periodic_shells().
 * @param camera Camera `prepared` was built with.
 * @param direction Unit view direction.
 * @param limits Step and layer budgets.
 * @param appearance Layer shading and compositing; needs a palette.
 * @return Shaded color and trace record; INVALID_QUERY on bad input or a
 *         frame without PreparedPeriodicShells::march.
 */
template <int DIMENSIONS>
HS_HOT_FLASH_MEMBER Raycast::ShadedTrace shade_periodic_shells_dimension(
    const PreparedPeriodicShells &prepared,
    const Raycast::PreparedCamera &camera, const math::Vector &direction,
    const Raycast::TraceLimits &limits, const Raycast::Appearance &appearance) {
  static_assert(DIMENSIONS == 3 || DIMENSIONS == 4);
  Raycast::ShadedTrace result{};
  const float cell_size = prepared.geometry.cell_size;
  if (!prepared.valid || !prepared.march || !Raycast::finite(direction) ||
      fabsf(math::dot(direction, direction) - 1.0f) >= 1e-4f ||
      !appearance.palette) {
    result.trace.status = Raycast::TraceStatus::INVALID_QUERY;
    return result;
  }
  const auto direction_value =
      camera.embedding.apply({{direction.x, direction.y, direction.z, 0}});
  const auto origin = camera.point4(direction * camera.radial_start);
  float a = 0;
  for (int k = 0; k < DIMENSIONS; ++k)
    a += direction_value[k] * direction_value[k];
  const float inverse_a = 1 / a;
  math::Vec4 cell{}, boundary{}, step{};
  for (int k = 0; k < DIMENSIONS; ++k) {
    const float p = origin[k] + camera.interval.near * direction_value[k];
    cell[k] = floorf(p / cell_size + .5f);
    if (direction_value[k] < 0 && p == (cell[k] - .5f) * cell_size)
      --cell[k];
    boundary[k] =
        direction_value[k] == 0
            ? INFINITY
            : ((cell[k] + copysignf(.5f, direction_value[k])) * cell_size -
               origin[k]) /
                  direction_value[k];
    step[k] = direction_value[k] == 0 ? INFINITY
                                      : cell_size / fabsf(direction_value[k]);
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
      end = fminf(end, boundary[k]);
      local[k] = origin[k] - cell[k] * cell_size;
    }
    float b = 0;
    for (int k = 0; k < DIMENSIONS; ++k)
      b += local[k] * direction_value[k];
    PeriodicShells::Intersections roots =
        prepared.geometry.intersect(local, direction_value, DIMENSIONS);
    const bool verified = roots.hit;
    const float filtered_coverage = admit_shell_tangent<DIMENSIONS>(
        prepared, camera, local, direction_value, b, inverse_a, start, end,
        roots);
    if (roots.hit && !composite_shell_cell<DIMENSIONS>(
                         prepared, camera, limits, appearance, roots, verified,
                         filtered_coverage, local, direction_value, a, b, start,
                         end, result, composite))
      return result;
    if (end >= camera.interval.far)
      break;
    if (!(end > start)) {
      result.trace.status = Raycast::TraceStatus::UNRESOLVED;
      break;
    }
    for (int k = 0; k < DIMENSIONS; ++k)
      if (boundary[k] <= end) {
        cell[k] += copysignf(1.0f, direction_value[k]);
        boundary[k] += step[k];
      }
    start = end;
  }
  result.color = composite.finish();
  return result;
}

namespace detail::reference {

/**
 * @brief shade_periodic_shells_dimension() at the camera's domain dimension.
 * @param prepared Frame constants from prepare_periodic_shells().
 * @param camera Camera `prepared` was built with.
 * @param direction Unit view direction.
 * @param limits Step and layer budgets.
 * @param appearance Layer shading and compositing.
 * @return Shaded color and trace record.
 */
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

/**
 * @brief Prepares and shades one ray in a single call.
 * @param camera Camera to shade with.
 * @param direction Unit view direction.
 * @param cell_size World edge length of one lattice cell.
 * @param shell_radius Sphere radius in cells.
 * @param footprint Ray filter footprint.
 * @param limits Step and layer budgets.
 * @param appearance Layer shading and compositing.
 * @return Shaded color and trace record.
 */
inline Raycast::ShadedTrace shade_periodic_shells(
    const Raycast::PreparedCamera &camera, const math::Vector &direction,
    float cell_size, float shell_radius, Raycast::Footprint footprint,
    const Raycast::TraceLimits &limits, const Raycast::Appearance &appearance) {
  return shade_periodic_shells(
      prepare_periodic_shells(camera, cell_size, shell_radius, footprint),
      camera, direction, limits, appearance);
}

} // namespace detail::reference

} // namespace SDF
