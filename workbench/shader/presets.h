/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_SHADER_WORKBENCH

/**
 * @file presets.h
 * @brief The authored presets, their shared generators, and the static assertions holding every preset and every preset edge to the admission rules.
 */

#include "workbench/shader/admission.h"
#include "workbench/shader/config.h"
#include "workbench/shader/limits.h"
#include "workbench/shader/options.h"
#include "workbench/shader/resources.h"

namespace Workbench {

inline constexpr Slots GENERATED_SURFACE_NOISE_SLOTS{
    Function::GRID,
    Projection::STEREOGRAPHIC,
    ProjectionFramePolicy::SPIN_WANDER,
    SurfaceLens::GLITCH,
    {{WarpStageKind::NONE}, {WarpStageKind::NONE}},
    SignalWeight::PROJECTION,
    ValueTransfer::NONE,
    CoveragePolicy::OPAQUE,
    PaletteMode::TRIADIC,
    PeirceLayout::SQUARE,
    AiroceanLayout::VERTICAL,
    BonneHemisphere::NORTH,
    GnomonicHemispherePolicy::FOLDED,
    SurfaceNoise::DIRECT,
    SurfaceNoisePlacement::AFTER_LENS};
/** @brief Index of a warp parameter name in the per-position name tables. */
enum WarpParamName : uint8_t {
  WARP_NAME_TRANSLATION_X,
  WARP_NAME_TRANSLATION_Y,
  WARP_NAME_ROTATION,
  WARP_NAME_SCALE_X,
  WARP_NAME_SCALE_Y,
  WARP_NAME_SHEAR,
  WARP_NAME_FREQUENCY,
  WARP_NAME_FIELD_ANGLE,
  WARP_NAME_CENTER_X,
  WARP_NAME_CENTER_Y,
  WARP_NAME_RADIUS,
  WARP_NAME_TURNS,
  WARP_NAME_VECTOR_ANGLE,
  WARP_NAME_CELL_X,
  WARP_NAME_CELL_Y,
  WARP_NAME_OFFSET_X,
  WARP_NAME_OFFSET_Y,
  WARP_NAME_RADIAL_SCALE,
  WARP_NAME_RADIAL_PHASE,
  WARP_NAME_ANGULAR_PHASE,
  WARP_NAME_EDGE_WIDTH,
  WARP_NAME_CENTER_ORBIT,
  WARP_NAME_COUNT,
};

inline constexpr const char *OUTER_WARP_PARAM_NAMES[] = {
    "Planar Warp 1 Translation X", "Planar Warp 1 Translation Y",
    "Planar Warp 1 Rotation",      "Planar Warp 1 Scale X",
    "Planar Warp 1 Scale Y",       "Planar Warp 1 Shear",
    "Planar Warp 1 Frequency",     "Planar Warp 1 Field Angle",
    "Planar Warp 1 Center X",      "Planar Warp 1 Center Y",
    "Planar Warp 1 Radius",        "Planar Warp 1 Turns",
    "Planar Warp 1 Vector Angle",  "Planar Warp 1 Cell X",
    "Planar Warp 1 Cell Y",        "Planar Warp 1 Offset X",
    "Planar Warp 1 Offset Y",      "Planar Warp 1 Radial Scale",
    "Planar Warp 1 Radial Phase",  "Planar Warp 1 Angular Phase",
    "Planar Warp 1 Edge Width",    "Planar Warp 1 Center Orbit"};
inline constexpr const char *INNER_WARP_PARAM_NAMES[] = {
    "Planar Warp 2 Translation X", "Planar Warp 2 Translation Y",
    "Planar Warp 2 Rotation",      "Planar Warp 2 Scale X",
    "Planar Warp 2 Scale Y",       "Planar Warp 2 Shear",
    "Planar Warp 2 Frequency",     "Planar Warp 2 Field Angle",
    "Planar Warp 2 Center X",      "Planar Warp 2 Center Y",
    "Planar Warp 2 Radius",        "Planar Warp 2 Turns",
    "Planar Warp 2 Vector Angle",  "Planar Warp 2 Cell X",
    "Planar Warp 2 Cell Y",        "Planar Warp 2 Offset X",
    "Planar Warp 2 Offset Y",      "Planar Warp 2 Radial Scale",
    "Planar Warp 2 Radial Phase",  "Planar Warp 2 Angular Phase",
    "Planar Warp 2 Edge Width",    "Planar Warp 2 Center Orbit"};
static_assert(sizeof(OUTER_WARP_PARAM_NAMES) / sizeof(const char *) ==
                  WARP_NAME_COUNT,
              "outer warp name table must match WarpParamName");
static_assert(sizeof(INNER_WARP_PARAM_NAMES) / sizeof(const char *) ==
                  WARP_NAME_COUNT,
              "inner warp name table must match WarpParamName");
inline constexpr Params
authored_params(SourceParams source, WarpStageParams outer_warp,
                ProjectionParams projection, SurfaceLensParams surface_lens,
                ColorParams color, OuterCameraParams outer_camera) {
  const WarpStageParams inner_warp{0.1f, 0.0f, 0.0f};
  projection.wander = outer_camera.wander;
  color.hue_noise_speed = hs::clamp(color.hue_noise_speed, -HUE_NOISE_SPEED_MAX,
                                    HUE_NOISE_SPEED_MAX);
  return {source,       {outer_warp, inner_warp},
          projection,   surface_lens,
          {},           color,
          outer_camera, {}};
}

inline constexpr float snap_affine_winding(float value) {
  if (value <= -4.0f)
    return -4.0f;
  if (value >= 4.0f)
    return 4.0f;
  if (value != value)
    return value;
  return value < 0.0f ? static_cast<float>(static_cast<int>(value - 0.5f))
                      : static_cast<float>(static_cast<int>(value + 0.5f));
}

inline constexpr void normalize_config_ranges(Config &config) {
  config.params.color.hue_noise_speed =
      hs::clamp(config.params.color.hue_noise_speed, -HUE_NOISE_SPEED_MAX,
                HUE_NOISE_SPEED_MAX);
  auto snap_affine = [](const WarpStageSpec &spec, WarpStageParams &params) {
    if (spec.kind != WarpStageKind::AFFINE_FRAME)
      return;
    params.translation_x = snap_affine_winding(params.translation_x);
    params.translation_y = snap_affine_winding(params.translation_y);
  };
  snap_affine(config.slots.warp_program.outer, config.params.warp.outer);
  snap_affine(config.slots.warp_program.inner, config.params.warp.inner);
}

inline constexpr Config
wave_shear_generated_preset(float pattern_freq = 4.439f,
                            float complexity = 0.5f, float warp_strength = 0.5f,
                            float warp_speed = 0.015625f) {
  Slots slots = GENERATED_SURFACE_NOISE_SLOTS;
  slots.warp_program.outer.kind = WarpStageKind::WAVE_SHEAR;
  slots.surface_noise = SurfaceNoise::NONE;
  slots.coverage = CoveragePolicy::PROJECTION_WEIGHT_SQUARED;
  slots.palette_mapping = PaletteMapping::CUP;
  Params params =
      authored_params({pattern_freq, 0.245f, complexity, 0.0f, 0.0f, 0.0f},
                      {1.0f, warp_strength, warp_speed}, {1.0f, 0.0f, 0.0f}, {},
                      {0.292f, 0.6304219f, 0.0f}, {0.8f});
  params.projection.wander = 0.0f;
  params.color.palette_chroma = 0.788f;
  params.color.mapping_phase = -0.0f;
  return {slots, params};
}

inline constexpr Config kaleidoscope_mirror_preset() {
  const Slots slots{Function::TWIN_WAVE,
                    Projection::STEREOGRAPHIC,
                    ProjectionFramePolicy::SPIN_WANDER,
                    SurfaceLens::KALEIDOSCOPE,
                    {{WarpStageKind::NONE}, {WarpStageKind::MIRROR_TILE}},
                    SignalWeight::PROJECTION,
                    ValueTransfer::NONE,
                    CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
                    PaletteMode::TRIADIC};
  Params params = authored_params({4.9755f, 0.125f, 0.513f, 0.0f, 0.8f, 0.05f},
                                  {0.1f, 0.0f, 0.5f}, {4.971f, 0.0f, 1.0f}, {},
                                  {0.27f, 2.2033439f, -0.00040800002f}, {1.0f});
  params.color.palette_chroma = 0.361f;
  return {slots, params};
}

inline constexpr Config gnomonic_grid_mirror_preset(SurfaceLens lens) {
  Slots slots{Function::GRID,
              Projection::GNOMONIC,
              ProjectionFramePolicy::IDENTITY,
              lens,
              {{WarpStageKind::MIRROR_TILE}, {WarpStageKind::NONE}},
              SignalWeight::PROJECTION,
              ValueTransfer::NONE,
              CoveragePolicy::EDGE_FADE,
              PaletteMode::TRIADIC};
  slots.gnomonic_hemisphere = GnomonicHemispherePolicy::FOLDED;
  WarpStageParams outer_warp;
  outer_warp.rotation = 0.29530972f;
  outer_warp.cell_x = 5.381125f;
  outer_warp.cell_y = 1.0f;
  outer_warp.offset_x = 1.344f;
  outer_warp.offset_y = -1.456f;
  Params params = authored_params({3.565f, 0.235f, 0.0f, 1.0f, 1.0f, 0.0f},
                                  outer_warp, {1.4f, 0.0f}, {}, {}, {1.0f});
  params.value.edge_width = 0.5f;
  return {slots, params};
}

inline constexpr Config gnomonic_kaleidoscope_grid_mirror_preset() {
  Config config = gnomonic_grid_mirror_preset(SurfaceLens::KALEIDOSCOPE);
  config.params.color.palette_chroma = 0.4f;
  config.params.color.hue_shift_amount = 0.424f;
  config.params.color.hue_noise_scale = 2.2033439f;
  return config;
}

inline constexpr Config peirce_dodecahedral_generated_preset() {
  Slots slots{Function::GRID,
              Projection::PEIRCE_QUINCUNCIAL,
              ProjectionFramePolicy::SPIN_WANDER,
              SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
              {{WarpStageKind::NONE}, {WarpStageKind::NONE}},
              SignalWeight::PROJECTION,
              ValueTransfer::NONE,
              CoveragePolicy::EDGE_FADE,
              PaletteMode::TRIADIC};
  slots.peirce_layout = PeirceLayout::SQUARE;
  Params params =
      authored_params({5.0f, 0.1f, 0.5f, 0.0f, 0.8f, 0.0f}, {}, {1.0f, 0.0f},
                      {}, {0.319f, 1.0f, 0.05f / TWO_PI_F}, {1.0f});
  params.projection.central_meridian = 0.0f;
  params.projection.coordinate_scale = 1.0f;
  params.value.edge_width = 0.1f;
  return {slots, params};
}

inline constexpr Config gnomonic_wave_shear_grid_preset() {
  Slots slots{Function::GRID,
              Projection::GNOMONIC,
              ProjectionFramePolicy::IDENTITY,
              SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
              {{WarpStageKind::WAVE_SHEAR}, {WarpStageKind::MIRROR_TILE}},
              SignalWeight::PROJECTION,
              ValueTransfer::NONE,
              CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
              PaletteMode::TRIADIC};
  slots.gnomonic_hemisphere = GnomonicHemispherePolicy::FOLDED;
  WarpStageParams outer_warp;
  outer_warp.strength = -0.176f;
  outer_warp.speed = -0.00325f;
  outer_warp.frequency = 1.408f;
  outer_warp.field_angle = 2.2305307f;
  WarpStageParams inner_warp;
  inner_warp.rotation = 0.0f;
  inner_warp.cell_x = 1.0f;
  inner_warp.cell_y = 1.0f;
  inner_warp.offset_x = 0.0f;
  inner_warp.offset_y = 0.0f;
  Params params =
      authored_params({6.3287f, 0.04f, 1.704f, 0.0f, 0.8f, 0.027f}, outer_warp,
                      {2.311f, 0.0f}, {}, {0.721f, 1.0f, 0.0f}, {1.0f});
  params.warp.inner = inner_warp;
  return {slots, params};
}

inline constexpr Config gnomonic_affine_lattice_contour_preset() {
  const Slots slots{Function::PRIMITIVE_LATTICE,
                    Projection::GNOMONIC,
                    ProjectionFramePolicy::SPIN_WANDER,
                    SurfaceLens::NONE,
                    {{WarpStageKind::AFFINE_FRAME}, {WarpStageKind::NONE}},
                    SignalWeight::PROJECTION,
                    ValueTransfer::ISO_CONTOUR,
                    CoveragePolicy::PROJECTION_WEIGHT,
                    PaletteMode::TRIADIC};
  Params params;
  params.source.pattern_freq = 8.0f;
  params.source.speed = 0.075f;
  params.source.complexity = 0.009122372f;
  params.source.pattern_mix = 1.0f;
  params.source.secondary_rate = 1.146f;
  params.source.lattice_cell_scale = 1.22925f;
  params.source.lattice_shape_blend = 1.0f;
  params.source.lattice_softness = 0.1608203f;
  params.source.lattice_radius = 0.332981884f;
  params.warp.outer.scale = 50.7493f;
  params.warp.outer.strength = 30.0f;
  params.warp.outer.speed = 0.015625f;
  params.warp.outer.translation_x = 4.0f;
  params.warp.outer.translation_y = 4.0f;
  params.warp.outer.shear = -0.0f;
  params.warp.inner.scale = 0.1f;
  params.projection.spin_rate = 0.0208791979f;
  params.projection.wander = 0.00309175253f;
  params.value.iso_level = 0.138f;
  params.value.iso_width = 0.227034181f;
  params.value.band_count = 19;
  params.value.band_phase = 6.10725641f;
  params.value.cutout_threshold = 0.5f;
  params.value.cutout_softness = 0.05f;
  params.value.edge_width = 0.327f;
  params.color.hue_shift_amount = 0.398f;
  params.color.hue_noise_scale = 0.8300313f;
  params.color.hue_noise_speed = 0.000212000014f;
  params.outer_camera.wander = 1.0f;
  params.surface_noise.scale = 0.507492959f;
  params.surface_noise.strength = 0.5f;
  params.surface_noise.rate = 5.377579e-7f;
  return {slots, params};
}

inline constexpr Config sinusoidal_lattice_curl_preset(float noise_scale) {
  Slots slots{Function::PRIMITIVE_LATTICE,
              Projection::SINUSOIDAL,
              ProjectionFramePolicy::SPIN_WANDER,
              SurfaceLens::NONE,
              {{WarpStageKind::NONE}, {WarpStageKind::NONE}},
              SignalWeight::PROJECTION,
              ValueTransfer::NONE,
              CoveragePolicy::PROJECTION_WEIGHT,
              PaletteMode::TRIADIC};
  slots.palette_mapping = PaletteMapping::CUP;
  slots.brightness_envelope = BrightnessEnvelope::CUP;
  slots.surface_noise = SurfaceNoise::CURL;
  slots.surface_noise_placement = SurfaceNoisePlacement::BEFORE_LENS;
  Params params;
  params.source.pattern_freq = 3.52279997f;
  params.source.speed = 0.1f;
  params.source.complexity = 0.9f;
  params.source.pattern_mix = 1.0f;
  params.source.secondary_rate = 0.8f;
  params.source.lattice_cell_scale = 0.710265636f;
  params.source.lattice_shape_blend = 1.0f;
  params.source.lattice_softness = 0.455532223f;
  params.source.lattice_radius = 0.290762514f;
  params.warp.outer.strength = 1.0f;
  params.warp.outer.speed = 0.000343749998f;
  params.projection.singularity_fade = 20.0f;
  params.projection.wander = 1.0f;
  params.color.hue_shift_amount = 0.268000007f;
  params.color.hue_noise_scale = 2.0f;
  params.color.palette_chroma = 1.0f;
  params.color.mapping_phase = -0.165999994f;
  params.outer_camera.wander = 1.0f;
  params.surface_noise.scale = noise_scale;
  params.surface_noise.strength = 0.0759999976f;
  return {slots, params};
}

inline constexpr Config stereographic_prism_polar_wave_lattice_preset() {
  Slots slots{Function::PRIMITIVE_LATTICE,
              Projection::STEREOGRAPHIC,
              ProjectionFramePolicy::SPIN_WANDER,
              SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM,
              {{WarpStageKind::POLAR_CHART}, {WarpStageKind::WAVE_SHEAR}},
              SignalWeight::PROJECTION,
              ValueTransfer::NONE,
              CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
              PaletteMode::ANALOGOUS};
  slots.palette_mapping = PaletteMapping::CUP;
  slots.surface_noise_placement = SurfaceNoisePlacement::BEFORE_LENS;
  Params params;
  params.source.pattern_freq = 3.52279997f;
  params.source.speed = 0.1f;
  params.source.complexity = 0.9f;
  params.source.pattern_mix = 1.0f;
  params.source.secondary_rate = 0.8f;
  params.source.lattice_cell_scale = 0.774140596f;
  params.source.lattice_shape_blend = 1.0f;
  params.source.lattice_softness = 0.377608389f;
  params.source.lattice_radius = 0.290762514f;
  params.warp.outer.strength = 1.0f;
  params.warp.outer.speed = 0.000343749998f;
  params.warp.outer.translation_x = 4.0f;
  params.warp.inner.speed = 0.000999999931f;
  params.warp.inner.translation_x = -0.0f;
  params.warp.inner.shear = 0.75f;
  params.projection.singularity_fade = 2.27300000f;
  params.projection.wander = 1.0f;
  params.color.hue_shift_amount = 0.268000007f;
  params.color.hue_noise_scale = 2.0f;
  params.color.palette_chroma = 1.0f;
  params.color.mapping_phase = -0.165999994f;
  params.outer_camera.wander = 1.0f;
  params.surface_noise.scale = 3.73634386f;
  params.surface_noise.strength = 0.0759999976f;
  return {slots, params};
}

inline constexpr Config gnomonic_dodecahedral_vector_mirror_grid_preset() {
  Slots slots{Function::GRID,
              Projection::GNOMONIC,
              ProjectionFramePolicy::IDENTITY,
              SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
              {{WarpStageKind::VECTOR_NOISE}, {WarpStageKind::MIRROR_TILE}},
              SignalWeight::PROJECTION,
              ValueTransfer::NONE,
              CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
              PaletteMode::TRIADIC};
  slots.palette_mapping = PaletteMapping::CUP;
  slots.brightness_envelope = BrightnessEnvelope::CUP;
  Params params;
  params.source.pattern_freq = 4.9755f;
  params.source.speed = 0.04f;
  params.source.complexity = 1.704f;
  params.source.secondary_rate = 0.8f;
  params.source.angle_rate = 0.027f;
  params.warp.outer.strength = 0.138f;
  params.warp.outer.speed = -0.00005f;
  params.warp.outer.frequency = 1.408f;
  params.warp.outer.field_angle = 2.23053074f;
  params.warp.inner.speed = 0.00327999983f;
  params.projection.singularity_fade = 2.311f;
  params.projection.wander = 1.0f;
  params.color.hue_shift_amount = 0.721f;
  params.color.palette_chroma = 1.0f;
  params.color.brightness_depth = 0.655f;
  params.outer_camera.wander = 1.0f;
  return {slots, params};
}

inline constexpr Config stereographic_dodecahedral_grid_inner_mirror_preset() {
  Slots slots{Function::GRID,
              Projection::STEREOGRAPHIC,
              ProjectionFramePolicy::SPIN_WANDER,
              SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
              {{WarpStageKind::NONE}, {WarpStageKind::MIRROR_TILE}},
              SignalWeight::PROJECTION,
              ValueTransfer::NONE,
              CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
              PaletteMode::ANALOGOUS};
  slots.palette_mapping = PaletteMapping::CUP;
  Params params;
  params.source.pattern_freq = 2.82629991f;
  params.source.complexity = 0.513f;
  params.source.secondary_rate = 0.8f;
  params.source.angle_rate = 0.0269999988f;
  params.warp.outer.scale = 0.1f;
  params.warp.outer.speed = 0.5f;
  params.warp.inner.scale = 0.1f;
  params.warp.inner.speed = 0.00013f;
  params.warp.inner.cell_y = 0.997703135f;
  params.projection.singularity_fade = 3.432f;
  params.projection.wander = 1.0f;
  params.color.hue_shift_amount = 0.366f;
  params.color.hue_noise_scale = 1.47215629f;
  params.color.palette_chroma = 1.0f;
  params.outer_camera.wander = 1.0f;
  return {slots, params};
}

inline constexpr Config
stereographic_dodecahedral_complex_grid_inner_mirror_preset() {
  Config config = stereographic_dodecahedral_grid_inner_mirror_preset();
  config.params.source.complexity = 3.0f;
  config.params.source.pattern_mix = 1.0f;
  return config;
}

inline constexpr Config
stereographic_dodecahedral_double_mapping_grid_inner_mirror_preset() {
  Config config = stereographic_dodecahedral_complex_grid_inner_mirror_preset();
  config.params.source.pattern_freq = 3.9407f;
  config.params.projection.wander = 0.165f;
  config.params.color.mapping_frequency = 2.0f;
  return config;
}

inline constexpr Config
equirectangular_dodecahedral_double_mapping_grid_inner_mirror_preset() {
  Config config =
      stereographic_dodecahedral_double_mapping_grid_inner_mirror_preset();
  config.slots.projection = Projection::EQUIRECTANGULAR;
  config.params.projection.singularity_fade = 2.14f;
  return config;
}

inline constexpr Config
equirectangular_dodecahedral_grid_inner_mirror_preset() {
  Config config =
      equirectangular_dodecahedral_double_mapping_grid_inner_mirror_preset();
  config.params.color.mapping_frequency = 1.0f;
  return config;
}

inline constexpr Config
equirectangular_dodecahedral_fine_grid_inner_mirror_preset() {
  Config config = equirectangular_dodecahedral_grid_inner_mirror_preset();
  config.params.source.pattern_freq = 0.3985f;
  config.params.warp.inner.speed = 0.00058f;
  config.params.warp.inner.cell_y = 0.901890635f;
  config.params.color.mapping_frequency = 21.212f;
  return config;
}

inline constexpr Config
stereographic_hexagonal_prism_twin_wave_mirror_preset() {
  Slots slots{Function::TWIN_WAVE,
              Projection::STEREOGRAPHIC,
              ProjectionFramePolicy::SPIN_WANDER,
              SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM,
              {{WarpStageKind::NONE}, {WarpStageKind::MIRROR_TILE}},
              SignalWeight::PROJECTION,
              ValueTransfer::NONE,
              CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
              PaletteMode::ANALOGOUS};
  slots.palette_mapping = PaletteMapping::BELL;
  Params params = authored_params(
      {3.881f, 0.128598228f, 0.513f, 0.0f, 0.8f, 0.027f}, {0.1f, 0.0f, 0.5f},
      {4.971f, 0.0f, 1.0f}, {}, {0.226f, 1.47215629f, 0.000138f}, {1.0f});
  params.color.palette_chroma = 1.0f;
  params.color.mapping_frequency = 1.341f;
  params.color.mapping_phase = -1.0f;
  return {slots, params};
}

inline constexpr Config stereographic_alien_core_mirror_preset() {
  Slots slots{Function::GRID,
              Projection::STEREOGRAPHIC,
              ProjectionFramePolicy::IDENTITY,
              SurfaceLens::GLITCH,
              {{WarpStageKind::MIRROR_TILE}, {WarpStageKind::NONE}},
              SignalWeight::PROJECTION,
              ValueTransfer::NONE,
              CoveragePolicy::EDGE_FADE,
              PaletteMode::TRIADIC};
  slots.hue_shift = HueShiftMode::WARP_DISPLACEMENT;
  WarpStageParams outer_warp;
  outer_warp.rotation = 0.295309722f;
  outer_warp.cell_x = 5.381125f;
  outer_warp.offset_x = 1.344f;
  outer_warp.offset_y = -1.456f;
  Params params =
      authored_params({2.5477f, 0.235f, 1.854f, 0.0f, 1.0f, 0.0f}, outer_warp,
                      {1.4f, 0.0f}, {}, {2.048f, 1.0f, 0.0f}, {1.0f});
  params.value.edge_width = 0.5f;
  params.color.palette_chroma = 0.292f;
  return {slots, params};
}

inline constexpr Config stereographic_mobius_twin_wave_inner_mirror_preset() {
  Slots slots{Function::TWIN_WAVE,
              Projection::STEREOGRAPHIC,
              ProjectionFramePolicy::SPIN_WANDER,
              SurfaceLens::MOBIUS,
              {{WarpStageKind::NONE}, {WarpStageKind::MIRROR_TILE}},
              SignalWeight::PROJECTION,
              ValueTransfer::NONE,
              CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
              PaletteMode::COMPLEMENTARY};
  slots.brightness_envelope = BrightnessEnvelope::CUP;
  slots.hue_shift = HueShiftMode::WARP_DISPLACEMENT;
  Params params = authored_params({10.158f, 0.245f, 0.513f, 0.0f, 0.8f, 0.027f},
                                  {0.1f, 0.0f, 0.5f}, {2.102f, 0.0f, 1.0f}, {},
                                  {0.312f, 1.0f, 0.0f}, {1.0f});
  params.surface_lens.mobius = {-1.072f, 0.304f, 0.416f,      0.0f,
                                0.0f,    0.0f,   0.70710677f, 0.0f};
  params.color.palette_chroma = 0.398f;
  return {slots, params};
}

inline constexpr Config stereographic_mobius_animated_inner_mirror_preset() {
  Config config = stereographic_mobius_twin_wave_inner_mirror_preset();
  config.params.warp.inner.speed = 0.005875f;
  config.params.warp.inner.cell_x = 0.2791094f;
  config.params.warp.inner.cell_y = 6.810328f;
  return config;
}

inline constexpr std::array<Preset, 24> PRESETS = {{
    {wave_shear_generated_preset(),
     InversePipelineId::GLITCH_NOISE_GRID_WAVE_SHEAR},
    {kaleidoscope_mirror_preset(),
     InversePipelineId::KALEIDOSCOPE_TWIN_WAVE_INNER_MIRROR},
    {gnomonic_kaleidoscope_grid_mirror_preset(),
     InversePipelineId::GNOMONIC_KALEIDOSCOPE_GRID_MIRROR},
    {gnomonic_grid_mirror_preset(SurfaceLens::GLITCH),
     InversePipelineId::GNOMONIC_ALIEN_CORE_MIRROR},
    {peirce_dodecahedral_generated_preset(),
     InversePipelineId::PEIRCE_DODECAHEDRAL_GRID},
    {gnomonic_wave_shear_grid_preset(),
     InversePipelineId::GNOMONIC_DODECAHEDRAL_GRID_WAVE_MIRROR},
    {gnomonic_affine_lattice_contour_preset(),
     InversePipelineId::GNOMONIC_AFFINE_LATTICE_CONTOUR},
    {sinusoidal_lattice_curl_preset(1.78815627f),
     InversePipelineId::SINUSOIDAL_LATTICE_MELT},
    {sinusoidal_lattice_curl_preset(3.29720306f),
     InversePipelineId::SINUSOIDAL_LATTICE_MELT},
    {stereographic_prism_polar_wave_lattice_preset(),
     InversePipelineId::STEREOGRAPHIC_PRISM_POLAR_WAVE_LATTICE},
    {gnomonic_dodecahedral_vector_mirror_grid_preset(),
     InversePipelineId::GNOMONIC_DODECAHEDRAL_GRID_VECTOR_MIRROR},
    {stereographic_dodecahedral_grid_inner_mirror_preset(),
     InversePipelineId::STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR},
    {stereographic_hexagonal_prism_twin_wave_mirror_preset(),
     InversePipelineId::STEREOGRAPHIC_HEXAGONAL_PRISM_TWIN_WAVE_INNER_MIRROR},
    {stereographic_dodecahedral_complex_grid_inner_mirror_preset(),
     InversePipelineId::STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR},
    {stereographic_dodecahedral_double_mapping_grid_inner_mirror_preset(),
     InversePipelineId::STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR},
    {equirectangular_dodecahedral_double_mapping_grid_inner_mirror_preset(),
     InversePipelineId::EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR},
    {equirectangular_dodecahedral_grid_inner_mirror_preset(),
     InversePipelineId::EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR},
    {equirectangular_dodecahedral_fine_grid_inner_mirror_preset(),
     InversePipelineId::EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR},
    {stereographic_alien_core_mirror_preset(),
     InversePipelineId::STEREOGRAPHIC_ALIEN_CORE_MIRROR},
    {stereographic_mobius_twin_wave_inner_mirror_preset(),
     InversePipelineId::STEREOGRAPHIC_MOBIUS_TWIN_WAVE_INNER_MIRROR},
    {stereographic_mobius_animated_inner_mirror_preset(),
     InversePipelineId::STEREOGRAPHIC_MOBIUS_TWIN_WAVE_INNER_MIRROR},
    {wave_shear_generated_preset(3.1447f, 0.5f, 2.72f, 0.00690625f),
     InversePipelineId::GLITCH_NOISE_GRID_WAVE_SHEAR},
    {wave_shear_generated_preset(7.5227f, 1.698f, 0.0f, 0.00690625f),
     InversePipelineId::GLITCH_NOISE_GRID_WAVE_SHEAR},
    {wave_shear_generated_preset(8.8162f, 1.698f, 1.376f, 0.00559375f),
     InversePipelineId::GLITCH_NOISE_GRID_WAVE_SHEAR},
}};
static_assert(
    [] {
      for (const Preset &preset : PRESETS)
        if (!valid_config(preset.config) ||
            preset.pipeline == InversePipelineId::NONE)
          return false;
      return true;
    }(),
    "a ShaderWorkbench preset lies outside its registered range");
static_assert(
    [] {
      for (size_t index = 0; index < PRESETS.size(); ++index)
        if (!transition_admitted(PRESETS[index].config,
                                 PRESETS[(index + 1) % PRESETS.size()].config))
          return false;
      return true;
    }(),
    "a ShaderWorkbench preset edge lacks continuous transition admission");

} // namespace Workbench

#endif // HS_ENABLE_SHADER_WORKBENCH
