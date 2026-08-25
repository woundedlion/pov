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
  float diagonal;
  float z;
};

union PreparedWarpTransform {
  PreparedAffineFrame affine;
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

} // namespace Workbench

#endif // HS_ENABLE_SHADER_WORKBENCH
