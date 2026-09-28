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
  const BakedPalette *feature_palette = nullptr;
  Raycast::ColorMode color = Raycast::ColorMode::DEPTH;
};

struct Prepared {
  Raycast::PreparedCamera camera;
  Raycast::Footprint footprint;
  Raycast::Appearance appearance;
  Raycast::TraceLimits limits;
  SDF::OctetFramework octet;
  SDF::OctetFramework4 octet4;
  bool valid = false;
};

HS_FLASH_INLINE inline Prepared prepare(const Settings &settings) {
  Prepared result;
  if (!Raycast::finite(settings.cell_size) || settings.cell_size <= 0.0f ||
      !Raycast::finite(settings.near_fade) || settings.near_fade <= 0.0f ||
      !Raycast::finite(settings.aa_strength) || settings.aa_strength < 0.0f ||
      !Raycast::finite(settings.pixel_half_angle) ||
      settings.pixel_half_angle < 0.0f || !settings.palette ||
      (settings.color != Raycast::ColorMode::DEPTH &&
       settings.color != Raycast::ColorMode::AXIS) ||
      (settings.color == Raycast::ColorMode::AXIS && !settings.feature_palette))
    return result;

  result.camera.domain = settings.domain;
  result.camera.center = settings.center;
  result.camera.embedding = settings.embedding;
  result.camera.radial_start = settings.radial_start;
  result.camera.interval = {0.0f, settings.far_distance};
  result.footprint = {settings.pixel_half_angle * settings.aa_strength,
                      settings.radial_start};
  result.appearance = {
      1.0f / settings.far_distance,
      0.0f,
      1.0f / settings.near_fade,
      settings.color,
      settings.palette,
      settings.feature_palette,
      settings.domain == Raycast::SamplingDomain::SLICE_4D ? 12.0f : 6.0f};
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
  return result;
}

template <bool SLICE_4D>
HS_HOT_FLASH_MEMBER Raycast::ShadedTrace shade(const math::Vector &direction,
                                               const Prepared &prepared) {
  const auto RAY = prepared.camera.ray(direction);
  if (!prepared.valid || !RAY.valid() ||
      (prepared.camera.domain == Raycast::SamplingDomain::SLICE_4D) !=
          SLICE_4D) {
    Raycast::ShadedTrace result;
    result.trace.status = Raycast::TraceStatus::INVALID_QUERY;
    return result;
  }
  const auto AMBIENT_DIRECTION = prepared.camera.embedding.apply(
      {{direction.x, direction.y, direction.z, 0.0f}});
  if constexpr (SLICE_4D) {
    SDF::OctetEvents4 events(
        prepared.octet4, prepared.camera.point4(RAY.origin), AMBIENT_DIRECTION,
        RAY.interval, prepared.footprint);
    return Raycast::shade_events(events, RAY.interval, prepared.limits,
                                 prepared.appearance);
  } else {
    const Raycast::Ray WORLD_RAY{
        prepared.camera.point3(RAY.origin),
        {AMBIENT_DIRECTION[0], AMBIENT_DIRECTION[1], AMBIENT_DIRECTION[2]},
        RAY.interval};
    SDF::OctetEvents events(prepared.octet, WORLD_RAY, prepared.footprint);
    return Raycast::shade_events(events, RAY.interval, prepared.limits,
                                 prepared.appearance);
  }
}

} // namespace HyperLatticeExperimental
