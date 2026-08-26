/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_SHADER_WORKBENCH

/**
 * @file config.h
 * @brief The workbench's authored vocabulary: slot enums, the per-stage parameter families, and the Config they compose into.
 */

#include "core/math/noise_field.h"
#include "core/render/pullback.h"

namespace Workbench {

enum class Function : uint8_t {
  TWIN_WAVE,
  RINGS,
  SPIRAL,
  GRID,
  NOISE_CONTOUR,
  PRIMITIVE_LATTICE,
  NOISE_CONTOUR_SPHERE,
  SPHERICAL_RINGS,
  FRACTAL,
  TESSELLATION
};
enum class Projection : uint8_t {
  SINUSOIDAL,
  STEREOGRAPHIC,
  GNOMONIC,
  BONNE,
  PEIRCE_QUINCUNCIAL,
  AIROCEAN,
  EQUIRECTANGULAR
};
enum class PeirceLayout : uint8_t { DIAMOND, SQUARE, HORIZONTAL, VERTICAL };
enum class AiroceanLayout : uint8_t { VERTICAL, HORIZONTAL };
enum class BonneHemisphere : uint8_t { NORTH, SOUTH };
enum class GnomonicHemispherePolicy : uint8_t {
  FOLDED,
  FRONT_HEMISPHERE,
  BACK_HEMISPHERE
};
enum class SurfaceLens : uint8_t {
  NONE,
  GLITCH,
  TWIST,
  KALEIDOSCOPE,
  MOBIUS,
  KALEIDOSCOPE_TETRAHEDRAL,
  KALEIDOSCOPE_OCTAHEDRAL,
  KALEIDOSCOPE_DODECAHEDRAL,
  KALEIDOSCOPE_TRIANGULAR_PRISM,
  KALEIDOSCOPE_SQUARE_PRISM,
  KALEIDOSCOPE_PENTAGONAL_PRISM,
  KALEIDOSCOPE_HEXAGONAL_PRISM,
  KALEIDOSCOPE_OCTAGONAL_PRISM,
  TANGENT_NOISE = 255
};
enum class WarpEnvelope : uint8_t { FLAT, PROJECTION_WEIGHT, EDGE_FADE };
enum class PolarMode : uint8_t { LINEAR, LOGARITHMIC };
enum class CurlIntegrator : uint8_t { EULER_1, MIDPOINT_2, MIDPOINT_4 };
enum class SurfaceCurlIntegrator : uint8_t { EULER, MIDPOINT, MIDPOINT_2X };
enum class SurfaceNoise : uint8_t { NONE, DIRECT, CURL };
enum class SurfaceNoisePlacement : uint8_t { BEFORE_LENS, AFTER_LENS };
enum class WarpStageKind : uint8_t {
  NONE,
  AFFINE_FRAME,
  WAVE_SHEAR,
  VORTEX,
  VECTOR_NOISE,
  CURL_FLOW,
  MIRROR_TILE,
  POLAR_CHART,
  LEGACY_STEREO_NOISE = 255
};
struct WarpStageSpec {
  WarpStageKind kind;
  NoiseBasis basis = NoiseBasis::SIMPLEX;
  WarpEnvelope envelope = WarpEnvelope::FLAT;
  PolarMode polar_mode = PolarMode::LINEAR;
  CurlIntegrator curl_integrator = CurlIntegrator::EULER_1;
  uint8_t polar_harmonic = 1;
  int32_t seed = 1337;

  constexpr bool operator==(const WarpStageSpec &) const = default;
};
struct WarpProgram {
  WarpStageSpec outer;
  WarpStageSpec inner;

  constexpr bool operator==(const WarpProgram &) const = default;
};
enum class ProjectionFramePolicy : uint8_t { IDENTITY, SPIN_WANDER };
enum class SignalWeight : uint8_t { NONE, PROJECTION };
enum class ValueTransfer : uint8_t { NONE, RIDGE, ISO_CONTOUR, SMOOTH_BANDS };
enum class CoveragePolicy : uint8_t {
  OPAQUE,
  PROJECTION_WEIGHT_SQUARED,
  VALUE_CUTOUT,
  EDGE_FADE,
  PROJECTION_WEIGHT
};
enum class PaletteMode : uint8_t { TRIADIC, COMPLEMENTARY, ANALOGOUS };
enum class PaletteMapping : uint8_t { CUP, BELL, LINEAR, REVERSE };
using PaletteMappingWeights = Pullback::Color::PaletteMappingWeights;
enum class BrightnessEnvelope : uint8_t {
  NONE,
  CUP,
  BELL,
  ASCENDING,
  DESCENDING
};
using HueShiftMode = Pullback::Color::HueMode;

struct Slots {
  Function function;
  Projection projection;
  ProjectionFramePolicy projection_frame;
  SurfaceLens surface_lens;
  WarpProgram warp_program;
  SignalWeight signal_weight;
  ValueTransfer value_transfer;
  CoveragePolicy coverage;
  PaletteMode palette;
  PeirceLayout peirce_layout = PeirceLayout::SQUARE;
  AiroceanLayout airocean_layout = AiroceanLayout::VERTICAL;
  BonneHemisphere bonne_hemisphere = BonneHemisphere::NORTH;
  GnomonicHemispherePolicy gnomonic_hemisphere =
      GnomonicHemispherePolicy::FOLDED;
  SurfaceNoise surface_noise = SurfaceNoise::NONE;
  SurfaceNoisePlacement surface_noise_placement =
      SurfaceNoisePlacement::AFTER_LENS;
  HueShiftMode hue_shift = HueShiftMode::NOISE;
  PaletteMapping palette_mapping = PaletteMapping::LINEAR;
  BrightnessEnvelope brightness_envelope = BrightnessEnvelope::NONE;

  constexpr bool operator==(const Slots &) const = default;
};

struct SourceParams {
  float pattern_freq = 1.0f;
  float speed = 0.0f;
  float complexity = 0.0f;
  float pattern_mix = 0.0f;
  float secondary_rate = 0.0f;
  float angle_rate = 0.0f;
  float noise_scale = 1.0f;
  float noise_contrast = 0.0f;
  float noise_time_rate = 0.0f;
  float lattice_cell_scale = 1.0f;
  float lattice_shape_blend = 0.0f;
  float lattice_softness = 0.05f;
  float lattice_radius = 0.25f;
  NoiseBasis noise_basis = NoiseBasis::SIMPLEX;
  int32_t noise_seed = 2927;
  uint8_t ring_count = 6;
  float ring_thickness = 0.08f;
  float ring_softness = 0.02f;
  float ring_wander = 0.0f;
  float fractal_scale = 0.5f;
  uint8_t fractal_iterations = 8;
  float julia_mix = 0.0f;
  float julia_real = -0.8f;
  float julia_imaginary = 0.156f;
  float fractal_contours = 4.0f;
  float tessellation_cell_scale = 1.0f;
  float tessellation_line_thickness = 0.04f;
  float tessellation_line_softness = 0.02f;
  Pullback::Source::TessellationKind tessellation_kind =
      Pullback::Source::TessellationKind::TRIANGULAR;

  HS_COLD_MEMBER constexpr SourceParams() = default;

  constexpr SourceParams(float pattern_freq, float speed, float complexity,
                         float pattern_mix, float secondary_rate,
                         float angle_rate = 0.0f)
      : pattern_freq(pattern_freq), speed(speed), complexity(complexity),
        pattern_mix(pattern_mix), secondary_rate(secondary_rate),
        angle_rate(angle_rate) {}

  HS_COLD_MEMBER bool operator==(const SourceParams &) const = default;

  HS_COLD_MEMBER void lerp(const SourceParams &a, const SourceParams &b,
                           float t) {
    // Trips if the field set changes, so a new field cannot silently go
    // uninterpolated and unsnapped.
    static_assert(sizeof(SourceParams) == 116,
                  "SourceParams field set changed - update lerp");
    pattern_freq = hs::lerp(a.pattern_freq, b.pattern_freq, t);
    speed = hs::lerp(a.speed, b.speed, t);
    complexity = hs::lerp(a.complexity, b.complexity, t);
    pattern_mix = hs::lerp(a.pattern_mix, b.pattern_mix, t);
    secondary_rate = hs::lerp(a.secondary_rate, b.secondary_rate, t);
    angle_rate = hs::lerp(a.angle_rate, b.angle_rate, t);
    noise_scale = hs::lerp(a.noise_scale, b.noise_scale, t);
    noise_contrast = hs::lerp(a.noise_contrast, b.noise_contrast, t);
    noise_time_rate = hs::lerp(a.noise_time_rate, b.noise_time_rate, t);
    lattice_cell_scale =
        hs::lerp(a.lattice_cell_scale, b.lattice_cell_scale, t);
    lattice_shape_blend =
        hs::lerp(a.lattice_shape_blend, b.lattice_shape_blend, t);
    lattice_softness = hs::lerp(a.lattice_softness, b.lattice_softness, t);
    lattice_radius = hs::lerp(a.lattice_radius, b.lattice_radius, t);
    noise_basis = t < 1.0f ? a.noise_basis : b.noise_basis;
    noise_seed = t < 1.0f ? a.noise_seed : b.noise_seed;
    ring_count = t < 1.0f ? a.ring_count : b.ring_count;
    ring_thickness = hs::lerp(a.ring_thickness, b.ring_thickness, t);
    ring_softness = hs::lerp(a.ring_softness, b.ring_softness, t);
    ring_wander = hs::lerp(a.ring_wander, b.ring_wander, t);
    fractal_scale = hs::lerp(a.fractal_scale, b.fractal_scale, t);
    fractal_iterations = t < 1.0f ? a.fractal_iterations : b.fractal_iterations;
    julia_mix = hs::lerp(a.julia_mix, b.julia_mix, t);
    julia_real = hs::lerp(a.julia_real, b.julia_real, t);
    julia_imaginary = hs::lerp(a.julia_imaginary, b.julia_imaginary, t);
    fractal_contours = hs::lerp(a.fractal_contours, b.fractal_contours, t);
    tessellation_cell_scale =
        hs::lerp(a.tessellation_cell_scale, b.tessellation_cell_scale, t);
    tessellation_line_thickness = hs::lerp(a.tessellation_line_thickness,
                                           b.tessellation_line_thickness, t);
    tessellation_line_softness =
        hs::lerp(a.tessellation_line_softness, b.tessellation_line_softness, t);
    tessellation_kind = t < 1.0f ? a.tessellation_kind : b.tessellation_kind;
  }
};

struct WarpStageParams {
  float scale = 1.0f;
  float strength = 0.0f;
  float speed = 0.0f;
  float translation_x = 0.0f;
  float translation_y = 0.0f;
  float rotation = 0.0f;
  float scale_x = 1.0f;
  float scale_y = 1.0f;
  float shear = 0.0f;
  float frequency = 1.0f;
  float field_angle = 0.0f;
  float center_x = 0.0f;
  float center_y = 0.0f;
  float radius = 1.0f;
  float turns = 0.0f;
  float center_orbit_radius = 0.0f;
  float vector_angle = 0.0f;
  float cell_x = 1.0f;
  float cell_y = 1.0f;
  float offset_x = 0.0f;
  float offset_y = 0.0f;
  float radial_scale = 1.0f;
  float radial_phase = 0.0f;
  float angular_phase = 0.0f;
  float edge_width = 0.1f;

  HS_COLD_MEMBER constexpr WarpStageParams() = default;

  constexpr WarpStageParams(float scale, float strength, float speed)
      : scale(scale), strength(strength), speed(speed) {}

  HS_COLD_MEMBER bool operator==(const WarpStageParams &) const = default;

  HS_COLD_MEMBER void lerp(const WarpStageParams &a, const WarpStageParams &b,
                           float t, bool rotation_is_rate = false) {
    static_assert(sizeof(WarpStageParams) == 100,
                  "WarpStageParams field set changed - update lerp");
    scale = hs::lerp(a.scale, b.scale, t);
    strength = hs::lerp(a.strength, b.strength, t);
    speed = hs::lerp(a.speed, b.speed, t);
    translation_x = hs::lerp(a.translation_x, b.translation_x, t);
    translation_y = hs::lerp(a.translation_y, b.translation_y, t);
    rotation =
        rotation_is_rate
            ? hs::lerp(a.rotation, b.rotation, t)
            : interp::shortest_periodic(a.rotation, b.rotation, t, TWO_PI_F);
    scale_x = expf(hs::lerp(logf(a.scale_x), logf(b.scale_x), t));
    scale_y = expf(hs::lerp(logf(a.scale_y), logf(b.scale_y), t));
    shear = hs::lerp(a.shear, b.shear, t);
    frequency = hs::lerp(a.frequency, b.frequency, t);
    field_angle =
        interp::shortest_periodic(a.field_angle, b.field_angle, t, TWO_PI_F);
    center_x = hs::lerp(a.center_x, b.center_x, t);
    center_y = hs::lerp(a.center_y, b.center_y, t);
    radius = hs::lerp(a.radius, b.radius, t);
    turns = hs::lerp(a.turns, b.turns, t);
    center_orbit_radius =
        hs::lerp(a.center_orbit_radius, b.center_orbit_radius, t);
    vector_angle =
        interp::shortest_periodic(a.vector_angle, b.vector_angle, t, TWO_PI_F);
    cell_x = hs::lerp(a.cell_x, b.cell_x, t);
    cell_y = hs::lerp(a.cell_y, b.cell_y, t);
    offset_x = hs::lerp(a.offset_x, b.offset_x, t);
    offset_y = hs::lerp(a.offset_y, b.offset_y, t);
    radial_scale = hs::lerp(a.radial_scale, b.radial_scale, t);
    radial_phase = hs::lerp(a.radial_phase, b.radial_phase, t);
    angular_phase = hs::lerp(a.angular_phase, b.angular_phase, t);
    edge_width = hs::lerp(a.edge_width, b.edge_width, t);
  }
};

struct WarpParams {
  WarpStageParams outer;
  WarpStageParams inner;

  HS_COLD_MEMBER bool operator==(const WarpParams &) const = default;

  HS_COLD_MEMBER void lerp(const WarpParams &a, const WarpParams &b, float t,
                           const WarpProgram &program = WarpProgram{}) {
    outer.lerp(a.outer, b.outer, t,
               program.outer.kind == WarpStageKind::AFFINE_FRAME);
    inner.lerp(a.inner, b.inner, t,
               program.inner.kind == WarpStageKind::AFFINE_FRAME);
  }
};

struct ProjectionParams {
  float singularity_fade = 1.0f;
  float spin_rate = 0.0f;
  float wander = 0.0f;
  float central_meridian = 0.0f;
  float coordinate_scale = 1.0f;
  float bonne_standard_parallel = PI_F * 0.25f;
  float layout_scroll = 0.0f;

  HS_COLD_MEMBER constexpr ProjectionParams() = default;

  constexpr ProjectionParams(float singularity_fade, float spin_rate)
      : singularity_fade(singularity_fade), spin_rate(spin_rate) {}
  constexpr ProjectionParams(float singularity_fade, float spin_rate,
                             float wander)
      : singularity_fade(singularity_fade), spin_rate(spin_rate),
        wander(wander) {}

  HS_COLD_MEMBER bool operator==(const ProjectionParams &) const = default;

  HS_COLD_MEMBER void lerp(const ProjectionParams &a, const ProjectionParams &b,
                           float t) {
    static_assert(sizeof(ProjectionParams) == 28,
                  "ProjectionParams field set changed - update lerp");
    singularity_fade = hs::lerp(a.singularity_fade, b.singularity_fade, t);
    spin_rate = hs::lerp(a.spin_rate, b.spin_rate, t);
    wander = hs::lerp(a.wander, b.wander, t);
    central_meridian = interp::shortest_periodic(
        a.central_meridian, b.central_meridian, t, TWO_PI_F);
    coordinate_scale = hs::lerp(a.coordinate_scale, b.coordinate_scale, t);
    bonne_standard_parallel =
        hs::lerp(a.bonne_standard_parallel, b.bonne_standard_parallel, t);
    layout_scroll =
        interp::shortest_periodic(a.layout_scroll, b.layout_scroll, t, 1.0f);
  }
};

struct SurfaceLensParams {
  MobiusParams mobius{0.7071067811865475f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                      0.7071067811865475f, 0.0f};

  HS_COLD_MEMBER constexpr SurfaceLensParams() = default;

  HS_COLD_MEMBER bool operator==(const SurfaceLensParams &other) const {
    return mobius.a.re == other.mobius.a.re &&
           mobius.a.im == other.mobius.a.im &&
           mobius.b.re == other.mobius.b.re &&
           mobius.b.im == other.mobius.b.im &&
           mobius.c.re == other.mobius.c.re &&
           mobius.c.im == other.mobius.c.im &&
           mobius.d.re == other.mobius.d.re && mobius.d.im == other.mobius.d.im;
  }

  HS_COLD_MEMBER void lerp(const SurfaceLensParams &a,
                           const SurfaceLensParams &b, float t) {
    // Trips if the field set changes, so a new field cannot silently go
    // uninterpolated and unsnapped.
    static_assert(sizeof(SurfaceLensParams) == 32,
                  "SurfaceLensParams field set changed - update lerp");
    mobius = t < 1.0f ? a.mobius : b.mobius;
  }
};

struct SurfaceNoiseParams {
  NoiseBasis basis = NoiseBasis::SIMPLEX;
  SurfaceCurlIntegrator integrator = SurfaceCurlIntegrator::EULER;
  int32_t seed = 1337;
  float scale = 1.0f;
  float strength = 0.0f;
  float rate = 0.0f;
  float direction = 0.0f;

  HS_COLD_MEMBER bool operator==(const SurfaceNoiseParams &) const = default;

  HS_COLD_MEMBER void lerp(const SurfaceNoiseParams &a,
                           const SurfaceNoiseParams &b, float t) {
    static_assert(sizeof(SurfaceNoiseParams) == 24,
                  "SurfaceNoiseParams field set changed - update lerp");
    basis = t < 1.0f ? a.basis : b.basis;
    integrator = t < 1.0f ? a.integrator : b.integrator;
    seed = t < 1.0f ? a.seed : b.seed;
    scale = hs::lerp(a.scale, b.scale, t);
    strength = hs::lerp(a.strength, b.strength, t);
    rate = hs::lerp(a.rate, b.rate, t);
    direction = interp::shortest_periodic(a.direction, b.direction, t, 1.0f);
  }
};

struct ValueParams {
  float iso_level = 0.5f;
  float iso_width = 0.05f;
  uint8_t band_count = 4;
  float band_phase = 0.0f;
  float cutout_threshold = 0.5f;
  float cutout_softness = 0.05f;
  float edge_width = 0.1f;

  HS_COLD_MEMBER bool operator==(const ValueParams &) const = default;
  HS_COLD_MEMBER void lerp(const ValueParams &a, const ValueParams &b,
                           float t) {
    static_assert(sizeof(ValueParams) == 28,
                  "ValueParams field set changed - update lerp");
    iso_level = hs::lerp(a.iso_level, b.iso_level, t);
    iso_width = hs::lerp(a.iso_width, b.iso_width, t);
    band_count = t < 1.0f ? a.band_count : b.band_count;
    band_phase =
        interp::shortest_periodic(a.band_phase, b.band_phase, t, TWO_PI_F);
    cutout_threshold = hs::lerp(a.cutout_threshold, b.cutout_threshold, t);
    cutout_softness = hs::lerp(a.cutout_softness, b.cutout_softness, t);
    edge_width = hs::lerp(a.edge_width, b.edge_width, t);
  }
};

using ColorParams = Pullback::Color::ColorParams;

struct OuterCameraParams {
  float wander = 0.0f;

  HS_COLD_MEMBER bool operator==(const OuterCameraParams &) const = default;

  HS_COLD_MEMBER void lerp(const OuterCameraParams &a,
                           const OuterCameraParams &b, float t) {
    static_assert(sizeof(OuterCameraParams) == 4,
                  "OuterCameraParams field set changed - update lerp");
    wander = hs::lerp(a.wander, b.wander, t);
  }
};

struct Params {
  SourceParams source;
  WarpParams warp;
  ProjectionParams projection;
  SurfaceLensParams surface_lens;
  ValueParams value;
  ColorParams color;
  OuterCameraParams outer_camera;
  SurfaceNoiseParams surface_noise;

  HS_COLD_MEMBER Params() = default;

  constexpr Params(SourceParams source, WarpParams warp,
                   ProjectionParams projection, SurfaceLensParams surface_lens,
                   ValueParams value, ColorParams color,
                   OuterCameraParams outer_camera,
                   SurfaceNoiseParams surface_noise)
      : source(source), warp(warp), projection(projection),
        surface_lens(surface_lens), value(value), color(color),
        outer_camera(outer_camera), surface_noise(surface_noise) {}

  HS_COLD_MEMBER bool operator==(const Params &) const = default;

  HS_COLD_MEMBER void lerp(const Params &a, const Params &b, float t,
                           const Slots &slots = Slots{}) {
    source.lerp(a.source, b.source, t);
    warp.lerp(a.warp, b.warp, t, slots.warp_program);
    projection.lerp(a.projection, b.projection, t);
    surface_lens.lerp(a.surface_lens, b.surface_lens, t);
    surface_noise.lerp(a.surface_noise, b.surface_noise, t);
    value.lerp(a.value, b.value, t);
    color = Pullback::Fields::interpolate(a.color, b.color, t);
    color.palette_mapping =
        t < 1.0f ? a.color.palette_mapping : b.color.palette_mapping;
    outer_camera.lerp(a.outer_camera, b.outer_camera, t);
  }

  HS_COLD_MEMBER void lerp_staggered(const Params &a, const Params &b, float t,
                                     const Slots &slots = Slots{}) {
    const int phase_count = (a.source != b.source) + (a.warp != b.warp) +
                            (a.projection != b.projection) +
                            (a.surface_lens != b.surface_lens) +
                            (a.value != b.value) + (a.color != b.color) +
                            (a.outer_camera != b.outer_camera) +
                            (a.surface_noise != b.surface_noise);
    int phase = 0;
    source = a.source;
    warp = a.warp;
    projection = a.projection;
    surface_lens = a.surface_lens;
    value = a.value;
    color = a.color;
    outer_camera = a.outer_camera;
    surface_noise = a.surface_noise;
    if (a.source != b.source)
      source.lerp(a.source, b.source, phase_t(t, phase++, phase_count));
    if (a.warp != b.warp)
      warp.lerp(a.warp, b.warp, phase_t(t, phase++, phase_count),
                slots.warp_program);
    if (a.projection != b.projection)
      projection.lerp(a.projection, b.projection,
                      phase_t(t, phase++, phase_count));
    if (a.surface_lens != b.surface_lens)
      surface_lens.lerp(a.surface_lens, b.surface_lens,
                        phase_t(t, phase++, phase_count));
    if (a.value != b.value)
      value.lerp(a.value, b.value, phase_t(t, phase++, phase_count));
    if (a.color != b.color) {
      const float color_t = phase_t(t, phase++, phase_count);
      color = Pullback::Fields::interpolate(a.color, b.color, color_t);
      color.palette_mapping =
          color_t < 1.0f ? a.color.palette_mapping : b.color.palette_mapping;
    }
    if (a.outer_camera != b.outer_camera)
      outer_camera.lerp(a.outer_camera, b.outer_camera,
                        phase_t(t, phase++, phase_count));
    if (a.surface_noise != b.surface_noise)
      surface_noise.lerp(a.surface_noise, b.surface_noise,
                         phase_t(t, phase, phase_count));
  }

  HS_COLD_MEMBER static float phase_t(float t, int phase, int phase_count) {
    return ease_in_out_sin(hs::clamp(t * phase_count - phase, 0.0f, 1.0f));
  }
};

struct Config {
  Slots slots;
  Params params;

  HS_COLD_MEMBER constexpr Config() = default;
  constexpr Config(const Slots &slots, const Params &params)
      : slots(slots), params(params) {}

  HS_COLD_MEMBER bool operator==(const Config &) const = default;
};
using RequestedConfig = Config;

/** @brief Whether the stage scales its amplitude by the warp envelope. */
inline constexpr bool warp_uses_envelope(WarpStageKind kind) {
  return kind == WarpStageKind::WAVE_SHEAR ||
         kind == WarpStageKind::VECTOR_NOISE ||
         kind == WarpStageKind::CURL_FLOW;
}

HS_COLD_MEMBER inline constexpr bool warp_uses_noise(WarpStageKind kind) {
  return kind == WarpStageKind::VECTOR_NOISE ||
         kind == WarpStageKind::CURL_FLOW;
}

HS_COLD_MEMBER inline constexpr bool seam_sensitive_warp(WarpStageKind kind) {
  return kind == WarpStageKind::VECTOR_NOISE ||
         kind == WarpStageKind::CURL_FLOW;
}

inline constexpr bool strict_projection(Projection projection) {
  return projection == Projection::BONNE ||
         projection == Projection::PEIRCE_QUINCUNCIAL ||
         projection == Projection::AIROCEAN;
}

enum class InversePipelineId : uint8_t {
  GLITCH_NOISE_GRID_WAVE_SHEAR,
  KALEIDOSCOPE_TWIN_WAVE_INNER_MIRROR,
  GNOMONIC_KALEIDOSCOPE_GRID_MIRROR,
  GNOMONIC_ALIEN_CORE_MIRROR,
  PEIRCE_DODECAHEDRAL_GRID,
  GNOMONIC_DODECAHEDRAL_GRID_WAVE_MIRROR,
  GNOMONIC_AFFINE_LATTICE_CONTOUR,
  SINUSOIDAL_LATTICE_MELT,
  STEREOGRAPHIC_PRISM_POLAR_WAVE_LATTICE,
  GNOMONIC_DODECAHEDRAL_GRID_VECTOR_MIRROR,
  STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR,
  STEREOGRAPHIC_HEXAGONAL_PRISM_TWIN_WAVE_INNER_MIRROR,
  EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR,
  STEREOGRAPHIC_ALIEN_CORE_MIRROR,
  STEREOGRAPHIC_MOBIUS_TWIN_WAVE_INNER_MIRROR,
  COUNT,
  NONE = 0xff
};

// The manifest carries one row per enumerator.
inline constexpr size_t INVERSE_PROGRAM_COUNT =
    static_cast<size_t>(InversePipelineId::COUNT);

struct SelectedConfig {
  Config config;
  InversePipelineId pipeline;
};
using Preset = SelectedConfig;

HS_COLD_MEMBER inline constexpr bool is_noise_contour(Function function) {
  return function == Function::NOISE_CONTOUR ||
         function == Function::NOISE_CONTOUR_SPHERE;
}

HS_COLD_MEMBER inline constexpr bool is_sphere_source(Function function) {
  return function == Function::NOISE_CONTOUR_SPHERE ||
         function == Function::SPHERICAL_RINGS;
}

} // namespace Workbench

#endif // HS_ENABLE_SHADER_WORKBENCH
