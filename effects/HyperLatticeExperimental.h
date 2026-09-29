/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/render/ray/camera.h"
#include "core/render/ray/shade.h"
#include "core/render/sdf/framework.h"
#include "core/render/sdf/cellular_wire.h"
#include "core/render/sdf/affine_lattice.h"
#include "core/render/sdf/periodic_shells.h"

namespace HyperLatticeExperimental {

enum class Geometry : uint8_t {
  OCTET,
  DIAMOND,
  HEXAGONAL,
  RHOMBIC,
  AFFINE_CUBIC,
  SHELLS
};

/** @brief Frame settings; camera distances and near fading use world units. */
struct Settings {
  Raycast::SamplingDomain domain = Raycast::SamplingDomain::SPATIAL_3D;
  float cell_size = 1.0f;
  float wire_radius = 0.055f;
  float radial_start = 0.0f;
  float far_distance = 7.0f;
  float near_fade = 0.5f;
  float aa_strength = 1.0f;
  math::Vec4 center{};
  math::Mat4 embedding = math::Mat4::identity();
  float pixel_half_angle = 0.0f;
  const BakedPalette *palette = nullptr;
  Geometry geometry = Geometry::OCTET;
  float shear = .55f;
  float stretch = 1.4f;
  float shell_radius = .30f;
};

struct Prepared {
  Raycast::PreparedCamera camera;
  Raycast::Footprint footprint;
  Raycast::Appearance appearance;
  Raycast::TraceLimits limits;
  SDF::OctetFramework octet;
  SDF::OctetEvents::PreparedProjection octet_projection{};
  SDF::OctetFramework4 octet4;
  SDF::OctetEvents4::PreparedProjection octet4_projection{};
  const SDF::CellularWire::Geometry *cellular = nullptr;
  SDF::PreparedPeriodicShells periodic_shells;
  bool valid = false;
  Geometry geometry = Geometry::OCTET;
  float shear = .55f;
  float stretch = 1.4f;
  float shell_radius = .30f;
};

/** @brief Premultiplied ray color and how its traversal ended. */
struct Sample {
  Pixel color;
  Raycast::TraceStatus status = Raycast::TraceStatus::RANGE_COMPLETE;
};

HS_FLASH_INLINE inline Prepared prepare(const Settings &settings) {
  Prepared result;
  if (!Raycast::finite(settings.cell_size) || settings.cell_size <= 0.0f ||
      !Raycast::finite(settings.near_fade) || settings.near_fade <= 0.0f ||
      !Raycast::finite(settings.aa_strength) || settings.aa_strength < 0.0f ||
      !Raycast::finite(settings.pixel_half_angle) ||
      settings.pixel_half_angle < 0.0f || !settings.palette)
    return result;

  result.camera.domain = settings.domain;
  result.camera.center = settings.center;
  result.camera.embedding = settings.embedding;
  result.camera.radial_start = settings.radial_start;
  result.camera.interval = {0.0f, settings.far_distance};
  result.footprint = {settings.pixel_half_angle * settings.aa_strength,
                      settings.radial_start};
  result.appearance = {1.0f / settings.far_distance, 0.0f,
                       1.0f / settings.near_fade, settings.palette};
  result.limits.max_candidates = 64;
  result.limits.max_layers = 32;
  result.octet.cell_size = settings.cell_size;
  result.octet.wire_radius = settings.wire_radius;
  result.octet4.cell_size = settings.cell_size;
  result.octet4.wire_radius = settings.wire_radius;
  result.geometry = settings.geometry;
  result.shear = settings.shear;
  result.stretch = settings.stretch;
  result.shell_radius = settings.shell_radius;
  result.valid = result.camera.valid() &&
                 Raycast::finite(result.footprint.angular_radius) &&
                 Raycast::finite(result.appearance.inv_far) &&
                 Raycast::finite(result.appearance.near_inv_span) &&
                 (settings.domain == Raycast::SamplingDomain::SLICE_4D
                      ? result.octet4.valid()
                      : result.octet.valid());
  if (!result.valid)
    return result;
  if (settings.geometry != Geometry::OCTET) {
    const bool CELLULAR = settings.geometry == Geometry::DIAMOND ||
                          settings.geometry == Geometry::HEXAGONAL ||
                          settings.geometry == Geometry::RHOMBIC;
    result.valid =
        (CELLULAR && settings.domain == Raycast::SamplingDomain::SPATIAL_3D) ||
        settings.geometry == Geometry::AFFINE_CUBIC ||
        settings.geometry == Geometry::SHELLS;
    if (CELLULAR)
      result.cellular = &SDF::CellularWire::geometry(
          settings.geometry == Geometry::DIAMOND
              ? SDF::CellularWire::Kind::DIAMOND
          : settings.geometry == Geometry::HEXAGONAL
              ? SDF::CellularWire::Kind::HEXAGONAL
              : SDF::CellularWire::Kind::RHOMBIC);
    if (settings.geometry == Geometry::SHELLS) {
      result.periodic_shells =
          SDF::prepare_periodic_shells(result.camera, settings.cell_size,
                                       settings.shell_radius, result.footprint);
      result.valid = result.periodic_shells.valid;
    }
    return result;
  }
  const auto &E = settings.embedding.m;
  if (settings.domain == Raycast::SamplingDomain::SPATIAL_3D) {
    const auto FAMILIES = result.octet.plane_families();
    auto &projection = result.octet_projection;
    const float SPACING = FAMILIES[0].spacing;
    projection.spacing2 = SPACING * SPACING;
    projection.wire_radius = result.octet.wire_radius;
    const math::Vector ORIGIN(settings.center[0], settings.center[1],
                              settings.center[2]);
    for (size_t i = 0; i < FAMILIES.size(); ++i) {
      const math::Vector NORMAL = FAMILIES[i].normal / SPACING;
      projection.offsets[i] = math::dot(ORIGIN - result.octet.origin, NORMAL);
      projection.normals[i] = {
          NORMAL.x * E[0][0] + NORMAL.y * E[1][0] + NORMAL.z * E[2][0],
          NORMAL.x * E[0][1] + NORMAL.y * E[1][1] + NORMAL.z * E[2][1],
          NORMAL.x * E[0][2] + NORMAL.y * E[1][2] + NORMAL.z * E[2][2]};
    }
  } else {
    auto &projection = result.octet4_projection;
    projection.inverse_scale =
        1.0f / (SDF::OctetFramework4::HALF_CUBE * settings.cell_size);
    for (int i = 0; i < 4; ++i) {
      projection.embedding[i] = math::Vector(E[i][0], E[i][1], E[i][2]);
      projection.origin[i] = (settings.center[i] - result.octet4.origin[i]) *
                             projection.inverse_scale;
    }
  }
  return result;
}

/**
 * @brief Composites an octet adapter's plane crossings front to back.
 * @details Matches Raycast::trace_events over the adapter, whose single merge
 * group keeps the first distance and largest coverage of coincident crossings.
 */
template <typename Events>
__attribute__((always_inline)) inline Sample
trace(const Events &events, Raycast::Interval interval,
      const Raycast::TraceLimits &limits,
      const Raycast::Appearance &appearance) {
  constexpr size_t SLOTS = Events::OWNER_CAPACITY;
  constexpr float RELATIVE_TOLERANCE = 1.0e-4f;
  std::array<float, SLOTS> next;
  std::array<float, SLOTS> step;
  std::array<uint8_t, SLOTS> stream;
  size_t count = 0;
  for (uint8_t i = 0; i < Events::STREAM_COUNT; ++i)
    if (events.active(i)) {
      next[count] = events.next[i];
      step[count] = events.step[i];
      stream[count++] = i;
    }
  for (; count < SLOTS; ++count)
    next[count] = INFINITY;
  Sample result;
  LayerComposite composite;
  int candidates = 0;
  int layers = 0;
  bool pending = false;
  float pending_t = 0.0f;
  float pending_coverage = 0.0f;
  float group_end = 0.0f;
  const auto flush = [&]() __attribute__((always_inline)) {
    if (layers >= limits.max_layers) {
      result.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
      return false;
    }
    ++layers;
    HS_PROFILE_DEEP(hl_layer_composite);
    composite.add(appearance.color(pending_t),
                  pending_coverage * appearance.opacity(pending_t));
    pending = false;
    if (composite.saturated()) {
      result.status = Raycast::TraceStatus::SATURATED;
      return false;
    }
    return true;
  };
  while (true) {
    HS_PROFILE_DEEP(hl_event_step);
    size_t first = 0;
    for (size_t k = 1; k < SLOTS; ++k)
      if (next[k] < next[first])
        first = k;
    const float T = next[first];
    if (!(T <= interval.far))
      break;
    if (pending && T > group_end && !flush()) {
      result.color = composite.premultiplied();
      return result;
    }
    if (candidates >= limits.max_candidates) {
      result.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
      break;
    }
    ++candidates;
    if (T >= interval.near) {
      uint32_t feature;
      const float COVERAGE = events.coverage(stream[first], T, feature);
      if (COVERAGE <= 0.0f) {
        HS_PROFILE_DEEP(hl_event_miss);
      } else if (!pending) {
        pending = true;
        pending_t = T;
        pending_coverage = COVERAGE;
        group_end = T + RELATIVE_TOLERANCE * std::max(1.0f, T);
      } else if (COVERAGE > pending_coverage) {
        pending_coverage = COVERAGE;
      }
    }
    next[first] = T + step[first];
  }
  if (pending)
    flush();
  result.color = composite.premultiplied();
  return result;
}

/**
 * @brief trace() by per-stream walks and a distance-ordered covered list.
 * @details Each stream walks its own crossings with the same accumulated
 * distances trace() pops, and only covered crossings are insertion-sorted.
 * An uncovered crossing never closes a merge group, so the groups are runs of
 * covered crossings within the relative tolerance of their first distance. A
 * ray with more crossings than the candidate budget defers to trace(), which
 * truncates them in distance order.
 */
template <typename Events>
__attribute__((always_inline)) inline Sample
trace_sorted(const Events &events, Raycast::Interval interval,
             const Raycast::TraceLimits &limits,
             const Raycast::Appearance &appearance) {
  constexpr float RELATIVE_TOLERANCE = 1.0e-4f;
  constexpr int CAPACITY = 64;
  std::array<float, CAPACITY> distances;
  std::array<float, CAPACITY> coverages;
  int covered = 0;
  int crossings = 0;
  const int BUDGET = std::min(limits.max_candidates, CAPACITY);
  for (uint8_t stream = 0; stream < Events::STREAM_COUNT; ++stream) {
    if (!events.active(stream))
      continue;
    const float STEP = events.step[stream];
    for (float t = events.next[stream]; t <= interval.far; t += STEP) {
      if (++crossings > BUDGET)
        return trace(events, interval, limits, appearance);
      if (t < interval.near)
        continue;
      uint32_t feature;
      const float COVERAGE = events.coverage(stream, t, feature);
      if (!(COVERAGE > 0.0f))
        continue;
      int slot = covered++;
      for (; slot > 0 && distances[slot - 1] > t; --slot) {
        distances[slot] = distances[slot - 1];
        coverages[slot] = coverages[slot - 1];
      }
      distances[slot] = t;
      coverages[slot] = COVERAGE;
    }
  }
  Sample result;
  LayerComposite composite;
  int layers = 0;
  for (int index = 0; index < covered;) {
    if (layers >= limits.max_layers) {
      result.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
      break;
    }
    ++layers;
    const float T = distances[index];
    const float GROUP_END = T + RELATIVE_TOLERANCE * std::max(1.0f, T);
    float coverage = coverages[index];
    for (++index; index < covered && distances[index] <= GROUP_END; ++index)
      coverage = std::max(coverage, coverages[index]);
    composite.add(appearance.color(T), coverage * appearance.opacity(T));
    if (composite.saturated()) {
      result.status = Raycast::TraceStatus::SATURATED;
      break;
    }
  }
  result.color = composite.premultiplied();
  return result;
}

template <bool SLICE_4D>
HS_HOT_FLASH_MEMBER Sample shade(const math::Vector &direction,
                                 const Prepared &prepared) {
  const auto &camera = prepared.camera;
  if (!prepared.valid ||
      (camera.domain == Raycast::SamplingDomain::SLICE_4D) != SLICE_4D)
    return {{}, Raycast::TraceStatus::INVALID_QUERY};
  if (prepared.geometry != Geometry::OCTET) {
    if constexpr (!SLICE_4D)
      if (prepared.geometry == Geometry::SHELLS &&
          prepared.periodic_shells.single_owner) {
        const auto SHELLS = SDF::trace_periodic_shells_3d(
            prepared.periodic_shells, camera, direction, prepared.limits,
            prepared.appearance);
        return {SHELLS.color, SHELLS.status};
      }
    Raycast::ShadedTrace sample;
    if (prepared.geometry == Geometry::AFFINE_CUBIC)
      sample = SDF::shade_affine_lattice(
          camera, direction, prepared.octet.cell_size,
          prepared.octet.wire_radius, prepared.footprint, prepared.limits,
          prepared.appearance, prepared.shear, prepared.stretch);
    else if (prepared.geometry == Geometry::SHELLS)
      sample = SDF::shade_periodic_shells(prepared.periodic_shells, camera,
                                          direction, prepared.limits,
                                          prepared.appearance);
    else {
      sample = SDF::CellularWire::shade(
          *prepared.cellular, prepared.octet.cell_size,
          prepared.octet.wire_radius, camera, prepared.footprint,
          prepared.limits, prepared.appearance, direction);
    }
    return {sample.color.color * sample.color.alpha, sample.trace.status};
  }
  if constexpr (SLICE_4D) {
    const SDF::OctetEvents4 events(prepared.octet4, prepared.octet4_projection,
                                   direction, camera.radial_start,
                                   camera.interval.near, prepared.footprint);
    return trace_sorted(events, camera.interval, prepared.limits,
                        prepared.appearance);
  } else {
    const SDF::OctetEvents events(prepared.octet_projection, direction,
                                  camera.radial_start, camera.interval.near,
                                  prepared.footprint);
    return trace_sorted(events, camera.interval, prepared.limits,
                        prepared.appearance);
  }
}

} // namespace HyperLatticeExperimental
