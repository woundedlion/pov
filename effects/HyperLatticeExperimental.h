/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/render/ray/query.h"
#include "core/render/ray/shade.h"
#include "core/render/sdf/framework.h"
#include "core/render/sdf/periodic_surface.h"

namespace HyperLatticeExperimental {

enum class Pattern : uint8_t { TRIANGULAR, COSINE, GYROID };

/** @brief Frame settings; camera distances and near fading use world units. */
struct Settings {
  Pattern pattern = Pattern::TRIANGULAR;
  float cell_size = 1.0f;
  float wire_radius = 0.055f;
  float radial_start = 0.0f;
  float far_distance = 7.0f;
  float near_fade = 0.5f;
  float aa_strength = 1.0f;
  math::Vector center{};
  math::Mat4 embedding = math::Mat4::identity();
  float pixel_half_angle = 0.0f;
  const BakedPalette *palette = nullptr;
  const BakedPalette *feature_palette = nullptr;
  Raycast::ColorMode color = Raycast::ColorMode::DEPTH;
};

struct Prepared {
  Pattern pattern = Pattern::TRIANGULAR;
  Raycast::PreparedCamera camera;
  Raycast::Footprint footprint;
  Raycast::Appearance appearance;
  Raycast::TraceLimits limits;
  SDF::TriangularFramework framework;
  SDF::CosineSurface cosine;
  SDF::GyroidSurface gyroid;
  bool valid = false;
};

HS_FLASH_INLINE inline Prepared prepare(const Settings &settings) {
  Prepared result;
  result.pattern = settings.pattern;
  if ((settings.pattern != Pattern::TRIANGULAR &&
       settings.pattern != Pattern::COSINE &&
       settings.pattern != Pattern::GYROID) ||
      !Raycast::finite(settings.cell_size) || settings.cell_size <= 0.0f ||
      !Raycast::finite(settings.near_fade) || settings.near_fade <= 0.0f ||
      !Raycast::finite(settings.aa_strength) || settings.aa_strength < 0.0f ||
      !Raycast::finite(settings.pixel_half_angle) ||
      settings.pixel_half_angle < 0.0f || !settings.palette ||
      (settings.color != Raycast::ColorMode::DEPTH &&
       settings.color != Raycast::ColorMode::AXIS) ||
      (settings.color == Raycast::ColorMode::AXIS && !settings.feature_palette))
    return result;

  result.camera.center = {
      {settings.center.x, settings.center.y, settings.center.z, 0}};
  result.camera.embedding = settings.embedding;
  result.camera.radial_start = settings.radial_start;
  result.camera.interval = {0.0f, settings.far_distance};
  result.footprint = {settings.pixel_half_angle * settings.aa_strength,
                      settings.radial_start};
  result.appearance = {1.0f / settings.far_distance,
                       0.0f,
                       1.0f / settings.near_fade,
                       settings.color,
                       settings.palette,
                       settings.feature_palette,
                       4.0f};
  result.limits.max_candidates = 64;
  result.limits.max_layers = 32;
  result.limits.max_steps = settings.pattern == Pattern::GYROID ? 512 : 256;
  result.limits.max_queries = result.limits.max_steps;
  result.limits.max_refinements = result.limits.max_queries;
  result.limits.position_tolerance = settings.cell_size / 10000.0f;
  result.framework.cell_size = settings.cell_size;
  result.framework.layer_height = settings.cell_size;
  result.framework.wire_radius = settings.wire_radius;
  result.cosine.period = settings.cell_size;
  result.gyroid.period = settings.cell_size;
  result.valid =
      result.camera.valid() &&
      Raycast::finite(result.footprint.angular_radius) &&
      Raycast::finite(result.appearance.inv_far) &&
      Raycast::finite(result.appearance.near_inv_span) &&
      result.limits.position_tolerance > 0.0f &&
      (settings.pattern == Pattern::TRIANGULAR ? result.framework.valid()
       : settings.pattern == Pattern::COSINE   ? result.cosine.valid()
                                               : result.gyroid.valid());
  return result;
}

template <Pattern PATTERN>
HS_HOT_FLASH_MEMBER Raycast::ShadedTrace shade(const math::Vector &direction,
                                               const Prepared &prepared) {
  static_assert(PATTERN == Pattern::TRIANGULAR || PATTERN == Pattern::COSINE ||
                PATTERN == Pattern::GYROID);
  if (!prepared.valid || prepared.pattern != PATTERN) {
    Raycast::ShadedTrace result;
    result.trace.status = Raycast::TraceStatus::INVALID_QUERY;
    return result;
  }
  const auto RAY = prepared.camera.ray(direction);
  if constexpr (PATTERN == Pattern::TRIANGULAR) {
    const auto AMBIENT_DIRECTION = prepared.camera.embedding.apply(
        {{direction.x, direction.y, direction.z, 0.0f}});
    const Raycast::Ray WORLD_RAY{
        prepared.camera.point3(RAY.origin),
        {AMBIENT_DIRECTION[0], AMBIENT_DIRECTION[1], AMBIENT_DIRECTION[2]},
        RAY.interval};
    if (!WORLD_RAY.valid()) {
      Raycast::ShadedTrace result;
      result.trace.status = Raycast::TraceStatus::INVALID_QUERY;
      return result;
    }
    SDF::FrameworkEvents events(prepared.framework, WORLD_RAY,
                                prepared.footprint);
    return Raycast::shade_events(events, WORLD_RAY.interval, prepared.limits,
                                 prepared.appearance);
  } else if constexpr (PATTERN == Pattern::COSINE) {
    const Raycast::DomainQuery3<SDF::CosineSurface> QUERY{prepared.cosine,
                                                          prepared.camera};
    return Raycast::shade_surface(QUERY, RAY, prepared.footprint,
                                  prepared.limits, prepared.appearance);
  } else {
    const Raycast::DomainQuery3<SDF::GyroidSurface> QUERY{prepared.gyroid,
                                                          prepared.camera};
    return Raycast::shade_surface(QUERY, RAY, prepared.footprint,
                                  prepared.limits, prepared.appearance);
  }
}

} // namespace HyperLatticeExperimental
