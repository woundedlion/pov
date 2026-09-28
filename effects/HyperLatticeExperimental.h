/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/render/ray/camera.h"
#include "core/render/ray/shade.h"
#include "core/render/sdf/framework.h"

namespace HyperLatticeExperimental {

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
  bool valid = false;
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
  result.valid = result.camera.valid() &&
                 Raycast::finite(result.footprint.angular_radius) &&
                 Raycast::finite(result.appearance.inv_far) &&
                 Raycast::finite(result.appearance.near_inv_span) &&
                 (settings.domain == Raycast::SamplingDomain::SLICE_4D
                      ? result.octet4.valid()
                      : result.octet.valid());
  if (!result.valid)
    return result;
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
    const auto FAMILIES = result.octet4.plane_families();
    auto &projection = result.octet4_projection;
    const float INVERSE_SCALE =
        1.0f / (SDF::OctetFramework4::HALF_CUBE * settings.cell_size);
    const float INVERSE_SPACING = 1.0f / FAMILIES[0].spacing;
    for (int i = 0; i < 4; ++i) {
      projection.embedding[i] = math::Vector(E[i][0], E[i][1], E[i][2]);
      projection.origin[i] =
          (settings.center[i] - result.octet4.origin[i]) * INVERSE_SCALE;
    }
    for (size_t f = 0; f < FAMILIES.size(); ++f) {
      math::Vector normal{};
      float offset = 0.0f;
      for (int i = 0; i < 4; ++i) {
        const float N = FAMILIES[f].normal[i] * INVERSE_SPACING;
        normal = normal + projection.embedding[i] * N;
        offset += (settings.center[i] - result.octet4.origin[i]) * N;
      }
      projection.normals[f] = normal;
      projection.offsets[f] = offset;
    }
    for (auto &row : projection.embedding)
      row = row * INVERSE_SCALE;
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

template <bool SLICE_4D>
HS_HOT_FLASH_MEMBER Sample shade(const math::Vector &direction,
                                 const Prepared &prepared) {
  const auto &camera = prepared.camera;
  if (!prepared.valid ||
      (camera.domain == Raycast::SamplingDomain::SLICE_4D) != SLICE_4D)
    return {{}, Raycast::TraceStatus::INVALID_QUERY};
  if constexpr (SLICE_4D) {
    const SDF::OctetEvents4 events(prepared.octet4, prepared.octet4_projection,
                                   direction, camera.radial_start,
                                   camera.interval.near, prepared.footprint);
    return trace(events, camera.interval, prepared.limits, prepared.appearance);
  } else {
    const SDF::OctetEvents events(prepared.octet_projection, direction,
                                  camera.radial_start, camera.interval.near,
                                  prepared.footprint);
    return trace(events, camera.interval, prepared.limits, prepared.appearance);
  }
}

} // namespace HyperLatticeExperimental
