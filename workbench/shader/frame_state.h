/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_SHADER_WORKBENCH

/**
 * @file frame_state.h
 * @brief Per-frame runtime records: the prepared stage payloads, the resource bindings, and the immutable FrameState a shading pass reads.
 */

#include "workbench/shader/config.h"

namespace Workbench {

struct Blend {
  Params params;
  PaletteMappingWeights palette_mapping;
};

struct Choreo {
  uint16_t dwell_min;
  uint16_t dwell_max;
  uint16_t blend_frames;
  bool staggered;
};

struct SourceState {
  float primary;
  float secondary;
  float angle;
  float angle_cos;
  float angle_sin;
};

using ProjectedLookup = Pullback::PlaneSample;

struct SourceTraits {
  bool y_periodic;
  bool polar_angle_compatible;
};

/** Warp-program output: source-side coordinates plus accumulated path. */
struct PlanarWarpResult {
  Complex coords;
  float path_length;
};
using SurfaceNoiseResult = Pullback::SurfaceResult;
using PlanarWarpStageResult = Pullback::WarpStepResult;
using FieldSample = Pullback::FieldSample;

struct ClockState {
  float source_primary = 0.0f;
  float source_secondary = 0.0f;
  float source_angle = 0.0f;
  float projection_spin = 0.0f;
  float hue_noise_phase = 0.0f;
  float source_noise_time = 0.0f;
  float surface_noise_time = 0.0f;
  float warp_outer_phase = 0.0f;
  float warp_inner_phase = 0.0f;
  float warp_outer_rotation = 0.0f;
  float warp_inner_rotation = 0.0f;
  float palette_oscillation_phase = 0.0f;

  HS_COLD_MEMBER constexpr ClockState() = default;
  constexpr ClockState(float source_primary, float source_secondary,
                       float source_angle, float projection_spin,
                       float hue_noise_phase)
      : source_primary(source_primary), source_secondary(source_secondary),
        source_angle(source_angle), projection_spin(projection_spin),
        hue_noise_phase(hue_noise_phase) {}
};

struct PreparedTransforms {
  Quaternion projection_conj;
  Quaternion outer_conj;
};

struct PreparedAffineFrame {
  float translation_x;
  float translation_y;
  float scale_x;
  float scale_y;
  float shear;
};

struct PreparedMirrorTile {
  float offset_x;
  float offset_y;
};

struct PreparedVortex {
  float center_x;
  float center_y;
  float radius_sq;
  float angle_numerator;
};

struct PreparedNoiseLoop {
  Vector offset;
};

union PreparedWarpTransform {
  // Vector's user-provided default ctor deletes the union's; an initialized
  // variant member restores it.
  PreparedAffineFrame affine{};
  PreparedMirrorTile mirror;
  PreparedVortex vortex;
  PreparedNoiseLoop noise_loop;
};

struct PreparedWarpStage {
  float rotation_cos;
  float rotation_sin;
  PreparedWarpTransform transform;
};

struct PreparedWarpProgram {
  PreparedWarpStage outer;
  PreparedWarpStage inner;
};

struct PreparedSurfaceNoise {
  Vector loop_offset;
  float direction_cos;
  float direction_sin;
};

struct PreparedHueRotation {
  static constexpr size_t LUT_SIZE = Pullback::Color::HueRotationLutView::SIZE;

  Pixel *lut;
  bool active;
};

struct PreparedHueNoise {
  static constexpr size_t LUT_SIZE = Pullback::Color::HueNoiseLutView::SIZE;

  int8_t *lut;
  bool active;
};

struct ResourceBindings {
  const FastNoiseLite *outer_warp_noise;
  const FastNoiseLite *inner_warp_noise;
  const FastNoiseLite *source_noise;
  const FastNoiseLite *surface_noise;
  const FastNoiseLite *color_noise;
  const BakedPalette *generated_palette;
};

inline constexpr size_t MAX_NOISE_RESOURCES = 9;

/** @brief The interpretive backend's and stage kernels' per-frame
    scratch; the compiled pipelines carry theirs in the program's prepared
    blob instead. */
struct DynamicPrepared {
  SourceState source;
  Pullback::Source::PreparedSphericalRings spherical_rings;
  PreparedWarpProgram warp;
  PreparedSurfaceNoise surface_noise;
};

struct FrameState {
  Slots slots;
  Params params;
  PaletteMappingWeights palette_mapping;
  ClockState clocks;
  PreparedTransforms transforms;
  /** cos and sin of params.projection.central_meridian; the projections that
      rotate by it are per-pixel, the meridian is not. */
  float meridian_cos = 1.0f;
  float meridian_sin = 0.0f;
  PreparedHueRotation prepared_hue_rotation;
  PreparedHueNoise prepared_hue_noise;
  ResourceBindings resources;
  DynamicPrepared dynamic;

  /** @brief Writes the central meridian and its trig pair together.
      @details The pair is what the per-pixel projections rotate by;
               assigning params.projection.central_meridian alone leaves
               them stale and the shade diverges from the interpreter. */
  HS_COLD_MEMBER void set_central_meridian(float radians) {
    params.projection.central_meridian = radians;
    meridian_cos = cosf(radians);
    meridian_sin = sinf(radians);
  }
};

/** @brief Per-sample functor a rasterizer walks over the canvas. */
using ShadeFunction = Color4 (*)(const Vector &, const FrameState &,
                                 const void *);
struct FrameShader {
  const FrameState *frame;
  float alpha;
  ShadeFunction shade_function;
  const void *prepared;

  HS_FLASH_MEMBER Color4 operator()(const Vector &view) const {
    Color4 color = shade_function(view, *frame, prepared);
    color.alpha *= alpha;
    return color;
  }
};

struct EndpointRuntime {
  ClockState clocks{};
  Quaternion projection_wander;
  Quaternion outer_wander;
  Quaternion source_wander;
  PreparedTransforms transforms;

  HS_COLD_MEMBER EndpointRuntime() = default;
};

HS_FLASH_MEMBER inline SourceState
prepare_source_state(const ClockState &clocks) {
  return {clocks.source_primary, clocks.source_secondary, clocks.source_angle,
          fast_cosf(clocks.source_angle), fast_sinf(clocks.source_angle)};
}

HS_FLASH_MEMBER inline Pullback::Source::PreparedSphericalRings
prepare_spherical_rings(const EndpointRuntime &endpoint) {
  const Quaternion orientation =
      make_rotation(X_AXIS, endpoint.clocks.source_angle) *
      endpoint.source_wander;
  return {rotate(Y_AXIS, orientation), endpoint.clocks.source_primary};
}

HS_FLASH_MEMBER inline PreparedSurfaceNoise
prepare_surface_noise(const ClockState &clocks, const Params &params) {
  const float surface_direction = TWO_PI_F * params.surface_noise.direction;
  return {noise_sphere_loop_offset(clocks.surface_noise_time),
          cosf(surface_direction), sinf(surface_direction)};
}

HS_FLASH_MEMBER inline PreparedWarpStage
prepare_warp_stage(const WarpStageSpec &spec, const WarpStageParams &params,
                   float stage_phase, const Complex &source_period = Complex(),
                   float affine_rotation = 0.0f) {
  PreparedWarpStage prepared{};
  float rotation = params.rotation;
  if (spec.kind == WarpStageKind::VECTOR_NOISE)
    rotation = params.vector_angle;
  else if (spec.kind == WarpStageKind::WAVE_SHEAR)
    rotation = params.field_angle;
  if (spec.kind == WarpStageKind::AFFINE_FRAME) {
    const float phase = TWO_PI_F * wrap_t(stage_phase);
    const float phase_cos = cosf(phase);
    rotation = affine_rotation;
    prepared.transform.affine = {
        wrap_t(stage_phase) * params.translation_x * source_period.re,
        wrap_t(stage_phase) * params.translation_y * source_period.im,
        powf(params.scale_x, phase_cos),
        powf(params.scale_y, phase_cos),
        params.shear * phase_cos,
    };
  } else if (spec.kind == WarpStageKind::MIRROR_TILE) {
    prepared.transform.mirror = {
        wrap_t(params.offset_x / params.cell_x + stage_phase) * params.cell_x,
        wrap_t(params.offset_y / params.cell_y) * params.cell_y,
    };
  } else if (spec.kind == WarpStageKind::VORTEX) {
    const float orbit_phase = TWO_PI_F * stage_phase;
    prepared.transform.vortex = {
        params.center_x + params.center_orbit_radius * cosf(orbit_phase),
        params.center_y + params.center_orbit_radius * sinf(orbit_phase),
        params.radius * params.radius,
        TWO_PI_F * params.turns,
    };
  } else if (spec.kind == WarpStageKind::VECTOR_NOISE ||
             spec.kind == WarpStageKind::CURL_FLOW) {
    prepared.transform.noise_loop = {noise_projected_loop_offset(stage_phase)};
  }
  prepared.rotation_cos = cosf(rotation);
  prepared.rotation_sin = sinf(rotation);
  return prepared;
}

struct PreparedEndpoint {
  FrameState *frame;
  ShadeFunction shade;
  const void *prepared;
  InversePipelineId pipeline;
  float alpha;
#if defined(HS_PROFILE_ENABLE)
  size_t preset;
#endif
};

} // namespace Workbench

#endif // HS_ENABLE_SHADER_WORKBENCH
