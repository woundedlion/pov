/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file lattice_trace.h
 * @brief Lattice event tracing adapters. */

#include "core/render/sdf/octet_trace.h"
#include "core/render/sdf/periodic_shells.h"

/** @brief Event adapters for octet trusses and periodic shells. */
namespace SDF::LatticeTrace {

using SDF::OctetTrace::CrossingStorage;
using SDF::OctetTrace::Sample;

/** @brief Geometry selected by a lattice trace. */
enum class Geometry : uint8_t { OCTET, SHELLS };

/** @brief Execution path for a prepared lattice trace. */
enum class TracePath : uint8_t {
  SINGLE_SHELLS,  ///< 3D shells, one candidate sphere per layer.
  MARCHED_SHELLS, ///< Shells by a layer march.
  OCTET,          ///< Octet truss events.
  GENERIC         ///< Per-ray dispatch on the prepared frame.
};

/** @brief Frame settings; camera distances and near fading use world units. */
struct Settings {
  /// 3D space or a 3D slice of 4D space.
  Raycast::SamplingDomain domain = Raycast::SamplingDomain::SPATIAL_3D;
  float cell_size = 1.0f;     ///< World edge length of one lattice cell.
  float wire_radius = 0.055f; /**< Octet wire radius in world units. */
  float radial_start = 0.0f;  ///< Ray start distance from the center.
  float far_distance = 7.0f;  ///< Ray length beyond the start.
  float near_fade = 0.5f;     ///< Near fade-in span; must be positive.
  float aa_strength = 1.0f;   ///< Pixel-footprint multiplier; nonnegative.
  math::Vec4 center{};        ///< Camera position.
  /// Rotation embedding view directions in the domain.
  math::Mat4 embedding = math::Mat4::identity();
  float pixel_half_angle = 0.0f;         ///< Angular half-width of one pixel.
  const BakedPalette *palette = nullptr; ///< Required depth palette.
  Geometry geometry = Geometry::OCTET;   ///< Geometry to trace.
  float shell_radius = .30f;             ///< Shell sphere radius in cells.
  float gain = 1.0f; /**< Brightness scale of the whole frame. */
  /** Scratch the octet traces sort crossings in; required for OCTET. */
  CrossingStorage *crossings = nullptr;
  /** Scratch the shell march sorts layers in; required for SHELLS. */
  SDF::ShellLayerStorage *shell_layers = nullptr;
};

/** @brief Frame state the traces read, built by prepare(). */
struct Prepared {
  Raycast::PreparedCamera camera; ///< Ray camera.
  Raycast::Footprint footprint;   ///< Pixel filter footprint.
  Raycast::Appearance appearance; ///< Depth fade, palette and gain.
  Raycast::TraceLimits limits;    ///< Candidate and layer budgets.
  SDF::OctetFramework octet;      ///< 3D octet truss.
  /// 3D octet plane families projected into view space.
  SDF::OctetEvents::PreparedProjection octet_projection{};
  SDF::OctetFramework4 octet4; ///< 4D octet truss.
  /// 4D octet embedding and origin in half-cube units.
  SDF::OctetEvents4::PreparedProjection octet4_projection{};
  SDF::PreparedPeriodicShells periodic_shells; ///< Shell frame constants.
  bool valid = false;                   ///< Whether the frame can be traced.
  Geometry geometry = Geometry::OCTET;  ///< Geometry to trace.
  CrossingStorage *crossings = nullptr; ///< Octet crossing scratch.
  SDF::ShellLayerStorage *shell_layers = nullptr; ///< Shell layer scratch.
};

/**
 * @brief Selects a prepared frame's specialized trace.
 * @param prepared Trace state.
 * @return Specialized path, or GENERIC for invalid or analytic-fallback frames.
 */
HS_FLASH_INLINE inline TracePath trace_path(const Prepared &prepared) {
  if (prepared.valid && prepared.geometry == Geometry::SHELLS &&
      prepared.camera.domain != Raycast::SamplingDomain::SLICE_4D &&
      prepared.periodic_shells.single_owner)
    return TracePath::SINGLE_SHELLS;
  if (prepared.valid && prepared.geometry == Geometry::SHELLS &&
      prepared.periodic_shells.march)
    return TracePath::MARCHED_SHELLS;
  if (prepared.valid && prepared.geometry == Geometry::OCTET)
    return TracePath::OCTET;
  return TracePath::GENERIC;
}

/**
 * @brief Traces a ray using a frame-selected path.
 * @tparam PATH Specialized path, or GENERIC for guarded per-ray dispatch.
 * @tparam SLICE_4D Whether the camera samples a 4D slice.
 * @param direction Unit view direction.
 * @param prepared Trace state; specialized paths require a matching valid frame.
 * @return Premultiplied color and trace status.
 */
template <TracePath PATH, bool SLICE_4D>
HS_HOT_INLINE Sample shade_path(const math::Vector &direction,
                                const Prepared &prepared);

/**
 * @brief Validates settings and precomputes the frame's trace state.
 * @param settings Frame settings.
 * @return Trace state; `valid` is false for unusable settings.
 */
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
  result.crossings = settings.crossings;
  result.shell_layers = settings.shell_layers;
  result.valid = result.camera.valid() &&
                 (settings.geometry != Geometry::OCTET || settings.crossings) &&
                 Raycast::finite(result.footprint.angular_radius) &&
                 result.appearance.valid_for(result.camera.interval) &&
                 Raycast::finite(result.appearance.near_inv_span) &&
                 (settings.geometry != Geometry::OCTET ||
                  (settings.domain == Raycast::SamplingDomain::SLICE_4D
                       ? result.octet4.valid()
                       : result.octet.valid()));
  if (!result.valid)
    return result;
  if (settings.geometry != Geometry::OCTET) {
    result.periodic_shells =
        SDF::prepare_periodic_shells(result.camera, settings.cell_size,
                                     settings.shell_radius, result.footprint);
    result.valid = result.periodic_shells.valid && settings.shell_layers;
    return result;
  }
  const auto &e = settings.embedding.m;
  if (settings.domain == Raycast::SamplingDomain::SPATIAL_3D) {
    const auto families = result.octet.plane_families();
    auto &projection = result.octet_projection;
    const float spacing = families[0].spacing;
    projection.spacing2 = spacing * spacing;
    projection.wire_radius = result.octet.wire_radius;
    const math::Vector origin(settings.center[0], settings.center[1],
                              settings.center[2]);
    for (size_t i = 0; i < families.size(); ++i) {
      const math::Vector normal = families[i].normal / spacing;
      projection.offsets[i] = math::dot(origin - result.octet.origin, normal);
      projection.normals[i] = {
          normal.x * e[0][0] + normal.y * e[1][0] + normal.z * e[2][0],
          normal.x * e[0][1] + normal.y * e[1][1] + normal.z * e[2][1],
          normal.x * e[0][2] + normal.y * e[1][2] + normal.z * e[2][2]};
    }
  } else {
    auto &projection = result.octet4_projection;
    projection.inverse_scale =
        1.0f / (SDF::OctetFramework4::HALF_CUBE * settings.cell_size);
    for (int i = 0; i < 4; ++i) {
      projection.embedding[i] = math::Vector(e[i][0], e[i][1], e[i][2]);
      projection.origin[i] =
          fmodf((settings.center[i] - result.octet4.origin[i]) *
                    projection.inverse_scale,
                2.0f);
    }
  }
  return result;
}

/**
 * @brief shade() for a valid octet frame of the matching domain.
 * @tparam SLICE_4D Whether the camera samples a 4D slice.
 * @param direction Unit view direction.
 * @param prepared Frame trace state.
 * @return Premultiplied color and trace status.
 */
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

/**
 * @brief shade() for a valid shell frame of the matching domain.
 * @tparam SLICE_4D Whether the camera samples a 4D slice.
 * @param direction Unit view direction.
 * @param prepared Frame trace state.
 * @return Premultiplied color and trace status.
 */
template <bool SLICE_4D>
HS_HOT_FLASH_MEMBER Sample shade_shells(const math::Vector &direction,
                                        const Prepared &prepared) {
  constexpr int DIMENSIONS = SLICE_4D ? 4 : 3;
  const auto path = trace_path(prepared);
  if constexpr (!SLICE_4D)
    if (path == TracePath::SINGLE_SHELLS)
      return shade_path<TracePath::SINGLE_SHELLS, false>(direction, prepared);
  if (path == TracePath::MARCHED_SHELLS)
    return shade_path<TracePath::MARCHED_SHELLS, SLICE_4D>(direction, prepared);
  return SDF::trace_periodic_shells_march<DIMENSIONS, true, true>(
      prepared.periodic_shells, prepared.camera, direction, prepared.limits,
      prepared.appearance, *prepared.shell_layers);
}

/**
 * @brief One ray's premultiplied color for the frame's geometry.
 * @tparam SLICE_4D Whether the camera samples a 4D slice.
 * @param direction Unit view direction.
 * @param prepared Frame trace state.
 * @return Premultiplied color and trace status; INVALID_QUERY
 * for an invalid frame or mismatched domain.
 */
template <bool SLICE_4D>
HS_HOT_FLASH_MEMBER Sample shade(const math::Vector &direction,
                                 const Prepared &prepared) {
  const auto &camera = prepared.camera;
  if (!prepared.valid ||
      (camera.domain == Raycast::SamplingDomain::SLICE_4D) != SLICE_4D)
    return {{}, Raycast::TraceStatus::INVALID_QUERY};
  if (prepared.geometry != Geometry::OCTET)
    return shade_shells<SLICE_4D>(direction, prepared);
  return shade_path<TracePath::OCTET, SLICE_4D>(direction, prepared);
}

template <TracePath PATH, bool SLICE_4D>
HS_HOT_INLINE Sample shade_path(const math::Vector &direction,
                                const Prepared &prepared) {
  if constexpr (PATH == TracePath::SINGLE_SHELLS) {
    static_assert(!SLICE_4D);
    return SDF::trace_periodic_shells_3d<true>(
        prepared.periodic_shells, prepared.camera, direction, prepared.limits,
        prepared.appearance);
  } else if constexpr (PATH == TracePath::MARCHED_SHELLS) {
    return SDF::trace_periodic_shells_march<SLICE_4D ? 4 : 3, false, true>(
        prepared.periodic_shells, prepared.camera, direction, prepared.limits,
        prepared.appearance, *prepared.shell_layers);
  } else if constexpr (PATH == TracePath::OCTET) {
    return shade_octet<SLICE_4D>(direction, prepared);
  } else {
    return shade<SLICE_4D>(direction, prepared);
  }
}

} // namespace SDF::LatticeTrace
