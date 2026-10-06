/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file lattice_trace.h
 * @brief Lattice event tracing adapters. */

#include "render/sdf/octet_trace.h"
#include "render/sdf/periodic_shells.h"

/** @brief Event adapters for octet trusses and periodic shells. */
namespace SDF::LatticeTrace {

using SDF::OctetTrace::CrossingStorage;
using SDF::OctetTrace::Sample;

/** @brief Geometry selected by a lattice trace. */
enum class Geometry : uint8_t { OCTET = 0, SHELLS = 5 };

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
  float shell_radius = .30f;
  float gain = 1.0f; /**< Brightness scale of the whole frame. */
  /** Scratch the octet traces sort crossings in; required for OCTET. */
  CrossingStorage *crossings = nullptr;
  /** Scratch the shell march sorts layers in; required for SHELLS. */
  SDF::ShellLayerStorage *shell_layers = nullptr;
};

/** @brief Frame state the traces read, built by prepare(). */
struct Prepared {
  Raycast::PreparedCamera camera;
  Raycast::Footprint footprint;
  Raycast::Appearance appearance;
  Raycast::TraceLimits limits;
  SDF::OctetFramework octet;
  SDF::OctetEvents::PreparedProjection octet_projection{};
  SDF::OctetFramework4 octet4;
  SDF::OctetEvents4::PreparedProjection octet4_projection{};
  SDF::PreparedPeriodicShells periodic_shells;
  bool valid = false;
  Geometry geometry = Geometry::OCTET;
  float shell_radius = .30f;
  CrossingStorage *crossings = nullptr;
  SDF::ShellLayerStorage *shell_layers = nullptr;
};

/** @brief Validates settings and precomputes the frame's trace state. */
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
                       1.0f / settings.near_fade, settings.palette,
                       settings.gain};
  result.limits.max_candidates = 64;
  result.limits.max_layers = 32;
  result.octet.cell_size = settings.cell_size;
  result.octet.wire_radius = settings.wire_radius;
  result.octet4.cell_size = settings.cell_size;
  result.octet4.wire_radius = settings.wire_radius;
  result.geometry = settings.geometry;
  result.shell_radius = settings.shell_radius;
  result.crossings = settings.crossings;
  result.shell_layers = settings.shell_layers;
  result.valid = result.camera.valid() &&
                 (settings.geometry != Geometry::OCTET || settings.crossings) &&
                 Raycast::finite(result.footprint.angular_radius) &&
                 Raycast::finite(result.appearance.inv_far) &&
                 Raycast::finite(result.appearance.near_inv_span) &&
                 (settings.domain == Raycast::SamplingDomain::SLICE_4D
                      ? result.octet4.valid()
                      : result.octet.valid());
  if (!result.valid)
    return result;
  if (settings.geometry != Geometry::OCTET) {
    result.valid = settings.geometry == Geometry::SHELLS;
    if (settings.geometry == Geometry::SHELLS) {
      result.periodic_shells =
          SDF::prepare_periodic_shells(result.camera, settings.cell_size,
                                       settings.shell_radius, result.footprint);
      result.valid = result.periodic_shells.valid && settings.shell_layers;
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

/** @brief shade() for a valid octet frame of the matching domain. */
template <bool SLICE_4D>
HS_HOT_FLASH_MEMBER Sample shade_octet(const math::Vector &direction,
                                       const Prepared &prepared) {
  if constexpr (SLICE_4D) {
    return SDF::OctetTrace::trace_4d(
        direction, prepared.camera, prepared.octet4, prepared.octet4_projection,
        prepared.footprint, prepared.limits, prepared.appearance,
        *prepared.crossings);
  } else {
    return SDF::OctetTrace::trace_3d(direction, prepared.camera,
                                     prepared.octet_projection,
                                     prepared.footprint, prepared.limits,
                                     prepared.appearance, *prepared.crossings);
  }
}

/** @brief shade() for a valid shell frame one of the shell traces serves. */
template <bool SLICE_4D>
HS_HOT_FLASH_MEMBER Sample shade_shells(const math::Vector &direction,
                                        const Prepared &prepared) {
  const SDF::ShellSample SHELLS =
      !SLICE_4D && prepared.periodic_shells.single_owner
          ? SDF::trace_periodic_shells_3d(prepared.periodic_shells,
                                          prepared.camera, direction,
                                          prepared.limits, prepared.appearance)
          : SDF::trace_periodic_shells_march<SLICE_4D ? 4 : 3>(
                prepared.periodic_shells, prepared.camera, direction,
                prepared.limits, prepared.appearance, *prepared.shell_layers);
  return {SHELLS.color, SHELLS.status};
}

/** @brief One ray's premultiplied color for the frame's geometry. */
template <bool SLICE_4D>
HS_HOT_FLASH_MEMBER Sample shade(const math::Vector &direction,
                                 const Prepared &prepared) {
  const auto &camera = prepared.camera;
  if (!prepared.valid ||
      (camera.domain == Raycast::SamplingDomain::SLICE_4D) != SLICE_4D)
    return {{}, Raycast::TraceStatus::INVALID_QUERY};
  if (prepared.geometry != Geometry::OCTET) {
    if (prepared.geometry == Geometry::SHELLS &&
        ((!SLICE_4D && prepared.periodic_shells.single_owner) ||
         prepared.periodic_shells.march))
      return shade_shells<SLICE_4D>(direction, prepared);
    if (prepared.geometry != Geometry::SHELLS)
      return {{}, Raycast::TraceStatus::INVALID_QUERY};
    const auto sample =
        SDF::shade_periodic_shells(prepared.periodic_shells, camera, direction,
                                   prepared.limits, prepared.appearance);
    return {sample.color.color * sample.color.alpha, sample.trace.status};
  }
  return shade_octet<SLICE_4D>(direction, prepared);
}

} // namespace SDF::LatticeTrace
