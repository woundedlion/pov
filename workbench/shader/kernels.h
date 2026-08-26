/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_SHADER_WORKBENCH

/**
 * @file kernels.h
 * @brief The interpretive shading path: pulling a sphere sample back through camera, lens, projection and warp, evaluating the source function, and colorizing the result.
 */

#include "workbench/shader/frame_state.h"
#include "workbench/shader/limits.h"
#include "workbench/shader/options.h"

namespace Workbench {

// Declared ahead of first use: at namespace scope an unqualified
// call would otherwise bind to a same-named function at global scope.
__attribute__((always_inline)) inline Vector
apply_frame_free_lens(const Vector &v, SurfaceLens lens);
HS_FLASH_MEMBER inline SurfaceNoiseResult
apply_surface_noise_result(const Vector &v, const FrameState &frame);
HS_FLASH_MEMBER inline Complex gnomonic(const Vector &v);
HS_O3_FN inline float grid(const Complex &p, const SourceParams &params,
                           const SourceState &source);
__attribute__((always_inline)) inline Complex
mirror_tile(const Complex &input, const WarpStageParams &params,
            const PreparedWarpStage &prepared);
HS_FLASH_MEMBER inline Vector mobius_lens(const Vector &v,
                                          const MobiusParams &params);
HS_FLASH_MEMBER inline float primitive_lattice(const Complex &p,
                                               const SourceParams &params);
__attribute__((always_inline)) inline Vector
profiled_apply_lens(const Vector &v, const FrameState &frame);
__attribute__((always_inline)) inline ProjectedLookup
profiled_project_branch(const Vector &v, const FrameState &frame);
HS_FLASH_MEMBER inline float rings(const Complex &p, const SourceState &source);
HS_FLASH_MEMBER inline float
sample_function(Function function, const Complex &p, const SourceState &source);
HS_FLASH_MEMBER inline Color4 shade_projected(const ProjectedLookup &projected,
                                              const FrameState &frame);
HS_FLASH_MEMBER inline float spiral(const Complex &p,
                                    const SourceState &source);
HS_FLASH_MEMBER inline float twin_wave(const Complex &p,
                                       const SourceState &source);
inline float warp_envelope(const Pullback::ProjectionProvenance &provenance,
                           WarpEnvelope envelope, float edge_width);
HS_FLASH_MEMBER inline PlanarWarpStageResult
warp_stage_lookup(const Complex &input,
                  const Pullback::ProjectionProvenance &provenance,
                  const WarpStageSpec &spec, const WarpStageParams &params,
                  float stage_phase, const FastNoiseLite *stage_noise,
                  const PreparedWarpStage &prepared, bool path_length_required);

struct ShaderWorkbenchBinding;

__attribute__((always_inline)) inline Pullback::ProjectionResult
stereographic_lookup(const Vector &local, const FrameState &frame) {
  return Pullback::Projection::stereographic(
      local, frame.params.projection.singularity_fade);
}

/** @brief Whether this frame's colorizer reads displacement metadata. */
__attribute__((always_inline)) inline bool
tracks_displacement(const FrameState &frame) {
  return frame.prepared_hue_rotation.active &&
         frame.slots.hue_shift == HueShiftMode::WARP_DISPLACEMENT;
}

inline bool projection_edge_distance_required(const FrameState &frame) {
  const WarpProgram &program = frame.slots.warp_program;
  return frame.slots.coverage == CoveragePolicy::EDGE_FADE ||
         (program.outer.kind != WarpStageKind::NONE &&
          program.outer.envelope == WarpEnvelope::EDGE_FADE) ||
         (program.inner.kind != WarpStageKind::NONE &&
          program.inner.envelope == WarpEnvelope::EDGE_FADE);
}

__attribute__((always_inline)) inline Vector
outer_camera_lookup(const Vector &view, const FrameState &frame) {
  return rotate(view, frame.transforms.outer_conj);
}

#if HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND ||                              \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
inline ProjectedLookup surface_lens_project_lookup(const Vector &v,
                                                   const FrameState &frame) {
  const Slots &slots = frame.slots;
  Vector pre_lens = v;
  float surface_path_length = 0.0f;
  if (slots.surface_noise != SurfaceNoise::NONE &&
      slots.surface_noise_placement == SurfaceNoisePlacement::BEFORE_LENS) {
    HS_SB_STAGE_MARK(surface_start);
    const SurfaceNoiseResult displaced = apply_surface_noise_result(v, frame);
    pre_lens = displaced.sphere;
    surface_path_length = displaced.path_length;
    HS_SB_STAGE_SPAN(surface_noise, surface_start);
  }
  const Vector lensed = slots.surface_lens == SurfaceLens::NONE
                            ? pre_lens
                            : profiled_apply_lens(pre_lens, frame);
  Vector post_lens = lensed;
  if (slots.surface_noise != SurfaceNoise::NONE &&
      slots.surface_noise_placement == SurfaceNoisePlacement::AFTER_LENS) {
    HS_SB_STAGE_MARK(surface_start);
    const SurfaceNoiseResult displaced =
        apply_surface_noise_result(lensed, frame);
    post_lens = displaced.sphere;
    surface_path_length = displaced.path_length;
    HS_SB_STAGE_SPAN(surface_noise, surface_start);
  }
  ProjectedLookup projected = profiled_project_branch(post_lens, frame);
  projected.path_length = surface_path_length;
  return projected;
}
#endif

HS_FLASH_MEMBER inline Pullback::ProjectionResult
project_bonne(const Vector &local, const FrameState &frame) {
  return Pullback::Projection::bonne(
      local, frame.params.projection.central_meridian,
      (frame.slots.bonne_hemisphere == BonneHemisphere::NORTH ? 1.0f : -1.0f) *
          frame.params.projection.bonne_standard_parallel,
      frame.params.projection.coordinate_scale);
}

HS_FLASH_MEMBER inline Pullback::ProjectionResult
project_peirce(const Vector &local, const FrameState &frame) {
  if (frame.slots.peirce_layout == PeirceLayout::SQUARE &&
      frame.params.projection.central_meridian == 0.0f &&
      projection_edge_distance_required(frame))
    return Pullback::Projection::peirce_fast_square(
        local, frame.params.projection.coordinate_scale,
        frame.params.projection.singularity_fade);
  return Pullback::Projection::peirce(
      local, frame.params.projection.central_meridian,
      static_cast<uint8_t>(frame.slots.peirce_layout),
      frame.params.projection.layout_scroll,
      projection_edge_distance_required(frame),
      frame.params.projection.coordinate_scale,
      frame.params.projection.singularity_fade, frame.meridian_cos,
      frame.meridian_sin);
}

HS_FLASH_MEMBER inline Pullback::ProjectionResult
project_airocean(const Vector &local, const FrameState &frame) {
  return Pullback::Projection::airocean(
      local, frame.params.projection.central_meridian,
      frame.slots.airocean_layout == AiroceanLayout::HORIZONTAL,
      projection_edge_distance_required(frame),
      frame.params.projection.coordinate_scale);
}

HS_FLASH_MEMBER inline Pullback::ProjectionResult
project_sinusoidal(const Vector &local, const FrameState &frame) {
  return Pullback::Projection::folded_sinusoidal(
      local, frame.params.projection.central_meridian);
}

HS_FLASH_MEMBER inline Pullback::ProjectionResult
project_equirectangular(const Vector &local, const FrameState &frame) {
  return Pullback::Projection::equirectangular(
      local, frame.params.projection.central_meridian,
      frame.params.projection.singularity_fade);
}

inline constexpr Pullback::Projection::GnomonicHemisphere
pullback_gnomonic_hemisphere(GnomonicHemispherePolicy policy) {
  switch (policy) {
  case GnomonicHemispherePolicy::FOLDED:
    return Pullback::Projection::GnomonicHemisphere::FOLDED;
  case GnomonicHemispherePolicy::FRONT_HEMISPHERE:
    return Pullback::Projection::GnomonicHemisphere::FRONT;
  case GnomonicHemispherePolicy::BACK_HEMISPHERE:
    return Pullback::Projection::GnomonicHemisphere::BACK;
  }
  __builtin_unreachable();
}

HS_FLASH_MEMBER inline Pullback::ProjectionResult
project_gnomonic(const Vector &local, const FrameState &frame) {
  return Pullback::Projection::gnomonic(
      local, frame.params.projection.singularity_fade,
      pullback_gnomonic_hemisphere(frame.slots.gnomonic_hemisphere));
}

HS_FLASH_MEMBER inline Pullback::ProjectionResult
project_nonstereographic(const Vector &local, const FrameState &frame) {
  if (frame.slots.projection == Projection::BONNE)
    return project_bonne(local, frame);
  if (frame.slots.projection == Projection::PEIRCE_QUINCUNCIAL)
    return project_peirce(local, frame);
  if (frame.slots.projection == Projection::AIROCEAN)
    return project_airocean(local, frame);
  if (frame.slots.projection == Projection::SINUSOIDAL)
    return project_sinusoidal(local, frame);
  if (frame.slots.projection == Projection::EQUIRECTANGULAR)
    return project_equirectangular(local, frame);
  if (frame.slots.projection == Projection::GNOMONIC)
    return project_gnomonic(local, frame);

  __builtin_unreachable();
}

HS_FLASH_MEMBER inline ProjectedLookup project_branch(const Vector &v,
                                                      const FrameState &frame) {
  const Vector local = rotate(v, frame.transforms.projection_conj);
  const Pullback::ProjectionResult result =
      frame.slots.projection != Projection::STEREOGRAPHIC
          ? project_nonstereographic(local, frame)
          : stereographic_lookup(local, frame);
  return {result.coords, result.provenance, local, 0.0f};
}

__attribute__((always_inline)) inline ProjectedLookup
profiled_project_branch(const Vector &v, const FrameState &frame) {
  HS_SB_STAGE_MARK(stage_start);
  const ProjectedLookup projected = project_branch(v, frame);
  HS_SB_STAGE_SPAN(projection, stage_start);
  return projected;
}

HS_FLASH_MEMBER inline Pullback::ProjectionResult
finalize_projection(const Vector &local, const Complex &coords,
                    Projection projection, float singularity_fade,
                    GnomonicHemispherePolicy gnomonic_hemisphere =
                        GnomonicHemispherePolicy::FOLDED) {
  switch (projection) {
  case Projection::SINUSOIDAL:
    return {coords,
            {static_cast<uint8_t>(local.z < 0.0f), 0, 0, 1.0f, 1.0f,
             PROJECTION_FLAG_FOLDED}};
  case Projection::EQUIRECTANGULAR:
    return {
        coords,
        {0, 0, BOUNDARY_CUT, PI_F - std::fabs(coords.re),
         Pullback::Projection::equirectangular_weight(local, singularity_fade),
         0}};
  case Projection::GNOMONIC: {
    const bool in_domain =
        gnomonic_hemisphere == GnomonicHemispherePolicy::FOLDED ||
        (gnomonic_hemisphere == GnomonicHemispherePolicy::FRONT_HEMISPHERE
             ? local.y >= 0.0f
             : local.y < 0.0f);
    return {coords,
            {static_cast<uint8_t>(local.y < 0.0f),
             static_cast<uint8_t>(local.y < 0.0f),
             static_cast<uint8_t>(BOUNDARY_CUT | BOUNDARY_SINGULAR),
             std::fabs(local.y),
             Pullback::Projection::singularity_attenuation(
                 local.y * local.y, local.x * local.x + local.z * local.z,
                 singularity_fade),
             0, 0, 0, in_domain ? 1.0f : 0.0f}};
  }
  case Projection::STEREOGRAPHIC:
  case Projection::BONNE:
  case Projection::PEIRCE_QUINCUNCIAL:
  case Projection::AIROCEAN:
    break;
  }
  __builtin_unreachable();
}

/**
 * @brief Pulls plane coordinates back through both warp stages.
 * @param projected Projection output; supplies the stage input coordinates
 *        and the weight and edge distance the stage envelopes read.
 * @param frame Frame snapshot.
 * @return Source-side coordinates plus the accumulated path length,
 *         carried in from the projected sample and advanced per stage.
 * @details Pullback order is Planar Warp 1 then Planar Warp 2, the reverse
 * of the authored order `source -> Warp 2 -> Warp 1 -> projection`.
 */
#if HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND ||                              \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
HS_FLASH_MEMBER inline PlanarWarpResult
planar_warp_lookup(const ProjectedLookup &projected, const FrameState &frame) {
  const bool path_length_required = tracks_displacement(frame);
  // Accumulate like the typed chain: ((carrier + outer) + inner).
  float path_length = projected.path_length;
  const PlanarWarpStageResult outer = warp_stage_lookup(
      projected.coords, projected.provenance, frame.slots.warp_program.outer,
      frame.params.warp.outer, frame.clocks.warp_outer_phase,
      frame.resources.outer_warp_noise, frame.dynamic.warp.outer,
      path_length_required);
  path_length += outer.path_length;
  const PlanarWarpStageResult inner = warp_stage_lookup(
      outer.coords, projected.provenance, frame.slots.warp_program.inner,
      frame.params.warp.inner, frame.clocks.warp_inner_phase,
      frame.resources.inner_warp_noise, frame.dynamic.warp.inner,
      path_length_required);
  path_length += inner.path_length;
  return {inner.coords, path_length};
}
#endif

HS_FLASH_MEMBER inline PlanarWarpStageResult
finish_closed_form_warp(const Complex &input, const Complex &output,
                        bool path_length_required) {
  return Pullback::Warp::finish_closed_form(input, output,
                                            path_length_required);
}

HS_FLASH_MEMBER inline PlanarWarpStageResult
warp_affine_frame(const Complex &input, const WarpStageParams &,
                  const PreparedWarpStage &prepared,
                  bool path_length_required) {
  return Pullback::Warp::affine_frame(input, prepared, path_length_required);
}

HS_FLASH_MEMBER inline PlanarWarpStageResult
warp_wave_shear(const Complex &input, const WarpStageParams &params,
                float stage_phase, float amplitude,
                const PreparedWarpStage &prepared, bool path_length_required) {
  return Pullback::Warp::wave_shear(input, params, stage_phase, amplitude,
                                    prepared, path_length_required);
}

HS_FLASH_MEMBER inline PlanarWarpStageResult
warp_vortex(const Complex &input, const PreparedWarpStage &prepared,
            bool path_length_required) {
  return Pullback::Warp::vortex(input, prepared, path_length_required);
}

HS_FLASH_MEMBER inline PlanarWarpStageResult
warp_vector_noise(const Complex &input, const WarpStageSpec &spec,
                  const WarpStageParams &params, float amplitude,
                  const FastNoiseLite &noise, const PreparedWarpStage &prepared,
                  bool path_length_required) {
  return Pullback::Warp::vector_noise(input, params, amplitude, noise,
                                      spec.basis, prepared,
                                      path_length_required);
}

HS_FLASH_MEMBER inline PlanarWarpStageResult
warp_curl_flow(const Complex &input, const WarpStageSpec &spec,
               const WarpStageParams &params, float amplitude,
               const FastNoiseLite &noise, float stage_phase,
               bool path_length_required) {
  const uint8_t intervals = spec.curl_integrator == CurlIntegrator::EULER_1 ? 1
                            : spec.curl_integrator == CurlIntegrator::MIDPOINT_2
                                ? 2
                                : 4;
  return Pullback::Warp::curl_flow(input, noise, spec.basis, intervals,
                                   params.scale, amplitude, stage_phase,
                                   path_length_required);
}

HS_FLASH_MEMBER inline PlanarWarpStageResult
warp_polar_chart(const Complex &input, const WarpStageSpec &spec,
                 const WarpStageParams &params, float stage_phase) {
  return Pullback::Warp::polar_chart(input, params, stage_phase,
                                     spec.polar_mode == PolarMode::LOGARITHMIC,
                                     spec.polar_harmonic);
}

/**
 * @brief Pulls plane coordinates back through one warp stage.
 * @param input Coordinates entering the stage.
 * @param projected Projection output, read only for the envelope weight and
 *        edge distance.
 * @param spec Stage kind and its discrete options.
 * @param params Stage parameters, already canonicalized.
 * @param stage_phase Wrapped noise phase for this stage's clock.
 * @param stage_noise Noise resource bound to this stage; may be null for
 *        kinds that sample no noise.
 * @param prepared Per-frame precomputation for this stage.
 * @param path_length_required Whether the frame's colorizer reads the
 *        displacement scalar.
 * @return Stage output coordinates, the delta it applied, and the path length
 *         travelled.
 * @details Path length is the direct displacement for the closed-form kinds
 * and the integrated arc length for curl flow. It is zero when
 * @p path_length_required is false.
 */
HS_FLASH_MEMBER inline PlanarWarpStageResult warp_stage_lookup(
    const Complex &input, const Pullback::ProjectionProvenance &provenance,
    const WarpStageSpec &spec, const WarpStageParams &params, float stage_phase,
    const FastNoiseLite *stage_noise, const PreparedWarpStage &prepared,
    bool path_length_required) {
  if (spec.kind == WarpStageKind::NONE)
    return {input, Complex(), 0.0f};
  const float envelope =
      warp_envelope(provenance, spec.envelope, params.edge_width);
  const float amplitude = params.strength * envelope;
  switch (spec.kind) {
  case WarpStageKind::NONE:
  case WarpStageKind::LEGACY_STEREO_NOISE:
    break;
  case WarpStageKind::AFFINE_FRAME:
    return warp_affine_frame(input, params, prepared, path_length_required);
  case WarpStageKind::WAVE_SHEAR:
    return warp_wave_shear(input, params, stage_phase, amplitude, prepared,
                           path_length_required);
  case WarpStageKind::VORTEX:
    return warp_vortex(input, prepared, path_length_required);
  case WarpStageKind::VECTOR_NOISE:
    if (amplitude == 0.0f)
      return {input, Complex(), 0.0f};
    HS_CHECK(stage_noise != nullptr,
             "ShaderWorkbench vector warp has no noise resource");
    return warp_vector_noise(input, spec, params, amplitude, *stage_noise,
                             prepared, path_length_required);
  case WarpStageKind::CURL_FLOW:
    if (amplitude == 0.0f)
      return {input, Complex(), 0.0f};
    HS_CHECK(stage_noise != nullptr,
             "ShaderWorkbench curl warp has no noise resource");
    return warp_curl_flow(input, spec, params, amplitude, *stage_noise,
                          stage_phase, path_length_required);
  case WarpStageKind::MIRROR_TILE: {
    HS_SB_STAGE_MARK(mirror_start);
    const PlanarWarpStageResult result = finish_closed_form_warp(
        input, mirror_tile(input, params, prepared), path_length_required);
    HS_SB_STAGE_SPAN(mirror_tile, mirror_start);
    return result;
  }
  case WarpStageKind::POLAR_CHART:
    return warp_polar_chart(input, spec, params, stage_phase);
  }
  __builtin_unreachable();
}

inline float warp_envelope(const Pullback::ProjectionProvenance &provenance,
                           WarpEnvelope envelope, float edge_width) {
  return Pullback::Warp::envelope(provenance, edge_width,
                                  envelope == WarpEnvelope::PROJECTION_WEIGHT,
                                  envelope == WarpEnvelope::EDGE_FADE);
}

HS_FLASH_MEMBER inline Complex curl_vector(const Complex &p,
                                           const FastNoiseLite &noise,
                                           NoiseBasis basis, float scale,
                                           float phase) {
  return Pullback::Warp::curl_vector(p, noise, basis, scale, phase);
}

__attribute__((always_inline)) inline Complex
mirror_tile(const Complex &input, const WarpStageParams &params,
            const PreparedWarpStage &prepared) {
  return Pullback::Warp::mirror_tile_coords(input, params, prepared);
}

#if HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND ||                              \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
inline Complex condition_source_coords(const Complex &coords,
                                       const FrameState &frame) {
  if (is_noise_contour(frame.slots.function) ||
      frame.slots.function == Function::PRIMITIVE_LATTICE ||
      frame.slots.function == Function::FRACTAL ||
      frame.slots.function == Function::TESSELLATION)
    return coords;
  return stereo_pattern_args(coords, frame.params.source.pattern_freq);
}

__attribute__((always_inline)) inline float
sample_coverage(const ProjectedLookup &projected, const FrameState &frame) {
  float coverage = 1.0f;
  switch (frame.slots.coverage) {
  case CoveragePolicy::OPAQUE:
  case CoveragePolicy::VALUE_CUTOUT:
    break;
  case CoveragePolicy::PROJECTION_WEIGHT_SQUARED:
    coverage =
        projected.provenance.value_weight * projected.provenance.value_weight;
    break;
  case CoveragePolicy::EDGE_FADE:
    coverage = Pullback::ProjectionCoverage::edge_fade(
        projected.provenance, frame.params.value.edge_width);
    break;
  case CoveragePolicy::PROJECTION_WEIGHT:
    coverage = projected.provenance.value_weight;
    break;
  }
  return coverage;
}

HS_FLASH_MEMBER inline FieldSample
shape_nontrivial_material(FieldSample sampled, const FrameState &frame) {
  switch (frame.slots.value_transfer) {
  case ValueTransfer::NONE:
    break;
  case ValueTransfer::RIDGE:
    sampled.value = unit_bell(sampled.value);
    break;
  case ValueTransfer::ISO_CONTOUR:
    sampled.value = Pullback::Transfer::iso_contour(
        sampled.value, frame.params.value.iso_level,
        frame.params.value.iso_width);
    break;
  case ValueTransfer::SMOOTH_BANDS:
    sampled.value = Pullback::Transfer::smooth_bands(
        sampled.value, static_cast<float>(frame.params.value.band_count),
        frame.params.value.band_phase);
    break;
  }
  if (frame.slots.coverage != CoveragePolicy::VALUE_CUTOUT)
    return sampled;
  return Pullback::Kernel::coverage(
      sampled, Pullback::ValueCoverage::value_cutout(
                   sampled.value, frame.params.value.cutout_threshold,
                   frame.params.value.cutout_softness));
}

inline FieldSample shape_material(float field, const ProjectedLookup &projected,
                                  const PlanarWarpResult &warped,
                                  const FrameState &frame) {
  const float weight = frame.slots.signal_weight == SignalWeight::PROJECTION
                           ? projected.provenance.value_weight
                           : 1.0f;
  const ProjectedLookup sample_input{warped.coords, projected.provenance,
                                     projected.sphere, warped.path_length};
  const FieldSample sampled = Pullback::Kernel::sample(
      sample_input, field * weight, sample_coverage(projected, frame));
  if (frame.slots.value_transfer == ValueTransfer::NONE &&
      frame.slots.coverage == CoveragePolicy::OPAQUE)
    return sampled;
  return shape_nontrivial_material(sampled, frame);
}
#endif

/**
 * @brief Samples the selected noise-contour source coordinate.
 * @param q Projected or sphere noise coordinate.
 * @param frame Frame snapshot.
 * @return Signed field value in [-1, 1].
 */
HS_FLASH_MEMBER inline float sample_noise_contour(const Vector &q,
                                                  const FrameState &frame) {
  return Pullback::Source::noise_contour(*frame.resources.source_noise,
                                         frame.params.source.noise_basis, q,
                                         frame.params.source.noise_contrast);
}

#if HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND ||                              \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
inline Pullback::Source::SphericalRingsSourceParams
spherical_rings_params(const SourceParams &params) {
  return {static_cast<float>(params.ring_count),
          params.ring_thickness,
          params.ring_softness,
          params.speed,
          params.angle_rate,
          params.ring_wander};
}

inline Pullback::Source::FractalSourceParams
fractal_params(const SourceParams &params) {
  return {params.fractal_scale,   static_cast<float>(params.fractal_iterations),
          params.julia_mix,       params.julia_real,
          params.julia_imaginary, params.fractal_contours,
          params.speed,           params.angle_rate};
}

inline Pullback::Source::TessellationSourceParams
tessellation_params(const SourceParams &params) {
  return {params.tessellation_cell_scale, params.tessellation_line_thickness,
          params.tessellation_line_softness, params.angle_rate};
}

HS_FLASH_MEMBER inline float sample_source(const Complex &p,
                                           const ProjectedLookup &projected,
                                           const FrameState &frame) {
  if (frame.slots.function == Function::SPHERICAL_RINGS)
    return Pullback::Source::spherical_rings(
        projected.sphere, spherical_rings_params(frame.params.source),
        frame.dynamic.spherical_rings);
  if (frame.slots.function == Function::FRACTAL)
    return Pullback::Source::escape_fractal(
        p, fractal_params(frame.params.source), frame.dynamic.source);
  if (frame.slots.function == Function::TESSELLATION)
    return Pullback::Source::tessellation(
        p, tessellation_params(frame.params.source),
        frame.params.source.tessellation_kind, frame.dynamic.source);
  if (frame.slots.function == Function::GRID)
    return grid(p, frame.params.source, frame.dynamic.source);
  if (frame.slots.function == Function::NOISE_CONTOUR)
    return sample_noise_contour(
        noise_projected_coordinate(p, frame.params.source.noise_scale,
                                   frame.clocks.source_noise_time),
        frame);
  if (frame.slots.function == Function::NOISE_CONTOUR_SPHERE)
    return sample_noise_contour(
        noise_sphere_coordinate(projected.sphere,
                                frame.params.source.noise_scale,
                                frame.clocks.source_noise_time),
        frame);
  if (frame.slots.function == Function::PRIMITIVE_LATTICE)
    return primitive_lattice(p, frame.params.source);
  return sample_function(frame.slots.function, p, frame.dynamic.source);
}
#endif

HS_FLASH_MEMBER inline float primitive_lattice(const Complex &p,
                                               const SourceParams &params) {
  return Pullback::Source::primitive_lattice(p, params);
}

struct ColorStateProvider {
  using Binding = ShaderWorkbenchBinding;
  using FrameState = Workbench::FrameState;

  static Pullback::Color::PaletteMappingWeights
  mapping_weights(const FrameState &frame) {
    return frame.palette_mapping;
  }
  static float mapping_frequency(const FrameState &frame) {
    return frame.params.color.mapping_frequency;
  }
  static float mapping_phase(const FrameState &frame) {
    return frame.params.color.mapping_phase;
  }
  static float oscillation_depth(const FrameState &frame) {
    return frame.params.color.phase_oscillation_depth;
  }
  static float oscillation_phase(const FrameState &frame) {
    return frame.clocks.palette_oscillation_phase;
  }
  static const BakedPalette &palette(const FrameState &frame) {
    return *frame.resources.generated_palette;
  }
  static Pullback::Color::HueMode hue_mode(const FrameState &frame) {
    return static_cast<Pullback::Color::HueMode>(frame.slots.hue_shift);
  }
  static float hue_shift_amount(const FrameState &frame) {
    return frame.params.color.hue_shift_amount;
  }
  static Pullback::Color::HueRotationLutView
  hue_rotation(const FrameState &frame) {
    return {frame.prepared_hue_rotation.lut,
            frame.prepared_hue_rotation.active};
  }
  static Pullback::Color::HueNoiseLutView hue_noise(const FrameState &frame) {
    return {frame.prepared_hue_noise.lut, frame.prepared_hue_noise.active};
  }
  static Pullback::Color::BrightnessEnvelope
  brightness_envelope(const FrameState &frame) {
    return static_cast<Pullback::Color::BrightnessEnvelope>(
        frame.slots.brightness_envelope);
  }
  static float brightness_depth(const FrameState &frame) {
    return frame.params.color.brightness_depth;
  }
  static float opacity_low(const FrameState &frame) {
    return frame.params.color.opacity_low;
  }
  static float opacity_high(const FrameState &frame) {
    return frame.params.color.opacity_high;
  }
};

/**
 * @brief Maps a shaped material sample to a palette colour.
 * @param sample Shaped value, coverage, and surface coordinate.
 * @param frame Frame snapshot.
 * @return Colour whose alpha is the palette alpha scaled by the coverage.
 */
HS_FLASH_MEMBER inline Color4 colorize_generated(const FieldSample &sample,
                                                 const FrameState &frame) {
  using Policy = Pullback::Color::GeneratedPalette<ColorStateProvider>;
  return Policy::apply(sample, frame, Policy::prepare(frame));
}

__attribute__((always_inline)) inline float
palette_mapping_coordinate(float value, PaletteMapping mapping, float frequency,
                           float offset) {
  return Pullback::Color::palette_mapping_coordinate(
      value, static_cast<Pullback::Color::PaletteMapping>(mapping), frequency,
      offset);
}

__attribute__((always_inline)) inline float
brightness_envelope_gain(float value, BrightnessEnvelope envelope,
                         float depth) {
  return Pullback::Color::brightness_envelope_gain(
      value, static_cast<Pullback::Color::BrightnessEnvelope>(envelope), depth);
}

#if HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND ||                              \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
HS_O3_FN inline Color4 colorize(const FieldSample &sample,
                                const FrameState &frame) {
  return colorize_generated(sample, frame);
}
#endif

HS_FLASH_MEMBER inline void
prepare_hue_rotation_lut(PreparedHueRotation &prepared,
                         const BakedPalette &palette) {
  Pullback::Color::prepare_hue_rotation_lut(
      std::span<Pixel, Pullback::Color::HueRotationLutView::SIZE>(
          prepared.lut, PreparedHueRotation::LUT_SIZE),
      palette);
}

HS_FLASH_MEMBER inline Vector hue_noise_face_direction(int face, float u,
                                                       float v) {
  return Pullback::Color::hue_noise_face_direction(face, u, v);
}

HS_FLASH_MEMBER inline float
sample_hue_noise_lut(const PreparedHueNoise &prepared, const Vector &v) {
  return Pullback::Color::sample_hue_noise_lut({prepared.lut, prepared.active},
                                               v);
}

__attribute__((always_inline)) inline Pixel
sample_hue_rotation_lut(const PreparedHueRotation &prepared, float value,
                        float amount) {
  return Pullback::Color::sample_hue_rotation_lut(
      {prepared.lut, prepared.active}, value, amount);
}

HS_FLASH_MEMBER inline Vector apply_lens(const Vector &v,
                                         const FrameState &frame) {
  switch (frame.slots.surface_lens) {
  case SurfaceLens::NONE:
  case SurfaceLens::GLITCH:
  case SurfaceLens::TWIST:
  case SurfaceLens::KALEIDOSCOPE:
  case SurfaceLens::KALEIDOSCOPE_TETRAHEDRAL:
  case SurfaceLens::KALEIDOSCOPE_OCTAHEDRAL:
  case SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL:
  case SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM:
  case SurfaceLens::KALEIDOSCOPE_SQUARE_PRISM:
  case SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM:
  case SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM:
  case SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM:
    return apply_frame_free_lens(v, frame.slots.surface_lens);
  case SurfaceLens::MOBIUS:
    return mobius_lens(v, frame.params.surface_lens.mobius);
  case SurfaceLens::TANGENT_NOISE:
    break;
  }
  __builtin_unreachable();
}

__attribute__((always_inline)) inline Vector
profiled_apply_lens(const Vector &v, const FrameState &frame) {
  HS_SB_STAGE_MARK(stage_start);
  const Vector lensed = apply_lens(v, frame);
  HS_SB_STAGE_SPAN(lens, stage_start);
  return lensed;
}

HS_FLASH_MEMBER inline Vector surface_curl_field(const Vector &v,
                                                 const FrameState &frame) {
  return Pullback::Surface::curl_field(v, *frame.resources.surface_noise,
                                       frame.params.surface_noise.basis,
                                       frame.params.surface_noise.scale,
                                       frame.dynamic.surface_noise.loop_offset);
}

HS_FLASH_MEMBER inline SurfaceNoiseResult
apply_surface_noise_result(const Vector &v, const FrameState &frame) {
  const SurfaceNoiseParams &params = frame.params.surface_noise;
  const bool path_length_required = tracks_displacement(frame);
  if (frame.slots.surface_noise == SurfaceNoise::DIRECT) {
    return Pullback::Surface::direct_noise(
        v, *frame.resources.surface_noise, params.basis, params.scale,
        frame.dynamic.surface_noise.loop_offset, params.strength,
        frame.dynamic.surface_noise.direction_cos,
        frame.dynamic.surface_noise.direction_sin, path_length_required);
  }
  const Pullback::Surface::Integrator integrator =
      params.integrator == SurfaceCurlIntegrator::EULER
          ? Pullback::Surface::Integrator::EULER
      : params.integrator == SurfaceCurlIntegrator::MIDPOINT
          ? Pullback::Surface::Integrator::MIDPOINT
          : Pullback::Surface::Integrator::MIDPOINT_2X;
  return Pullback::Surface::curl_noise(v, *frame.resources.surface_noise,
                                       params.basis, integrator, params.scale,
                                       frame.dynamic.surface_noise.loop_offset,
                                       params.strength, path_length_required);
}

HS_FLASH_MEMBER inline Vector apply_surface_noise(const Vector &v,
                                                  const FrameState &frame) {
  return apply_surface_noise_result(v, frame).sphere;
}

/**
 * @brief Applies a lens whose image depends on the direction alone.
 * @param v Unit direction on the sphere.
 * @param lens Lens to apply; MOBIUS reads FrameState parameters and is
 *        rejected here.
 * @return The lensed direction.
 */
__attribute__((always_inline)) inline Vector
apply_frame_free_lens(const Vector &v, SurfaceLens lens) {
  switch (lens) {
  case SurfaceLens::NONE:
    return v;
  case SurfaceLens::GLITCH:
    return lenses::glitch_lens(v);
  case SurfaceLens::TWIST:
    return lenses::twist_lens(v);
  case SurfaceLens::KALEIDOSCOPE:
    return lenses::kaleidoscope_lens(v);
  case SurfaceLens::KALEIDOSCOPE_TETRAHEDRAL:
    return lenses::polyhedral_kaleidoscope_lens(v, lenses::TETRAHEDRAL_MIRRORS);
  case SurfaceLens::KALEIDOSCOPE_OCTAHEDRAL:
    return lenses::polyhedral_kaleidoscope_lens(v, lenses::OCTAHEDRAL_MIRRORS);
  case SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL:
    return lenses::dodecahedral_kaleidoscope_lens(v);
  case SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM:
    return lenses::polyhedral_kaleidoscope_lens(
        v, lenses::TRIANGULAR_PRISM_MIRRORS);
  case SurfaceLens::KALEIDOSCOPE_SQUARE_PRISM:
    return lenses::polyhedral_kaleidoscope_lens(v,
                                                lenses::SQUARE_PRISM_MIRRORS);
  case SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM:
    return lenses::polyhedral_kaleidoscope_lens(
        v, lenses::PENTAGONAL_PRISM_MIRRORS);
  case SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM:
    return lenses::polyhedral_kaleidoscope_lens(
        v, lenses::HEXAGONAL_PRISM_MIRRORS);
  case SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM:
    return lenses::polyhedral_kaleidoscope_lens(
        v, lenses::OCTAGONAL_PRISM_MIRRORS);
  case SurfaceLens::MOBIUS:
  case SurfaceLens::TANGENT_NOISE:
    HS_CHECK(false, "frame-parameterized lens needs the FrameState overload");
    __builtin_unreachable();
  }
  __builtin_unreachable();
}

HS_FLASH_MEMBER inline Vector mobius_lens(const Vector &v,
                                          const MobiusParams &params) {
  return Pullback::Lens::mobius(v, params);
}

/**
 * @brief Projects a sphere direction with the projections that take no
 *        `ProjectionParams`.
 * @param v Direction in the projection frame.
 * @param projection Must be equirectangular, stereographic, or gnomonic; the
 *        cartographic kernels read live parameters and are reached through
 *        `project_branch` instead.
 * @return Plane coordinates in the projection's native units.
 */
HS_FLASH_MEMBER inline Complex project_point(const Vector &v,
                                             Projection projection) {
  switch (projection) {
  case Projection::SINUSOIDAL:
    return projections::folded_sinusoidal(v);
  case Projection::EQUIRECTANGULAR:
    return projections::equirectangular(v);
  case Projection::STEREOGRAPHIC:
    return stereo(v);
  case Projection::GNOMONIC:
    return Workbench::gnomonic(v);
  case Projection::BONNE:
  case Projection::PEIRCE_QUINCUNCIAL:
  case Projection::AIROCEAN:
    break;
  }
  __builtin_unreachable();
}

HS_FLASH_MEMBER inline Complex gnomonic(const Vector &v) {
  return Pullback::Projection::gnomonic(
             v, 0.0f, Pullback::Projection::GnomonicHemisphere::FOLDED)
      .coords;
}

HS_FLASH_MEMBER inline float sample_function(Function function,
                                             const Complex &p,
                                             const SourceState &source) {
  switch (function) {
  case Function::TWIN_WAVE:
    return twin_wave(p, source);
  case Function::RINGS:
    return rings(p, source);
  case Function::SPIRAL:
    return spiral(p, source);
  case Function::GRID:
  case Function::NOISE_CONTOUR:
  case Function::NOISE_CONTOUR_SPHERE:
  case Function::PRIMITIVE_LATTICE:
  case Function::SPHERICAL_RINGS:
  case Function::FRACTAL:
  case Function::TESSELLATION:
    break;
  }
  __builtin_unreachable();
}

HS_FLASH_MEMBER inline float twin_wave(const Complex &p,
                                       const SourceState &source) {
  return Pullback::Source::twin_wave(p, source);
}

HS_FLASH_MEMBER inline float rings(const Complex &p,
                                   const SourceState &source) {
  return Pullback::Source::rings(p, source);
}

HS_FLASH_MEMBER inline float spiral(const Complex &p,
                                    const SourceState &source) {
  return Pullback::Source::spiral(p, source);
}

HS_O3_FN inline float grid(const Complex &p, const SourceParams &params,
                           const SourceState &source) {
  return Pullback::Source::grid(p, params, source);
}

#if HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND ||                              \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
/**
 * @brief Shades one sphere sample by pulling it back to a source coordinate.
 * @param view Unit direction of the visible sphere point.
 * @param frame Immutable snapshot of slots, parameters, clocks, transforms,
 *        and palette resources for this frame.
 * @return Straight-alpha colour for the sample.
 * @details Walks outer camera, surface lens, and projection backward. A
 * strict projection whose two lens branches land in different regions cannot
 * be joined in the plane, so the branches are shaded separately and their
 * outputs blended instead.
 */
inline Color4 shade_dynamic(const Vector &view, const FrameState &frame,
                            const void *) {
  const Vector outer_local = outer_camera_lookup(view, frame);
  const ProjectedLookup projected =
      surface_lens_project_lookup(outer_local, frame);
  return shade_projected(projected, frame);
}

/**
 * @brief Runs the planar half of the pullback: warps, samples, shapes, and
 *        colorizes.
 * @param projected Plane coordinates and seam metadata for the sample.
 * @param frame Frame snapshot.
 * @return Straight-alpha colour for the sample.
 */
HS_FLASH_MEMBER inline Color4 shade_projected(const ProjectedLookup &projected,
                                              const FrameState &frame) {
  HS_SB_STAGE_MARK(stage_start);
  const PlanarWarpResult warped = planar_warp_lookup(projected, frame);
  HS_SB_STAGE_SPAN(planar_warp, stage_start);
  const Complex source_coords = condition_source_coords(warped.coords, frame);
  const float field = sample_source(source_coords, projected, frame);
  HS_SB_STAGE_SPAN(source, stage_start);
  const FieldSample material = shape_material(field, projected, warped, frame);
  HS_SB_STAGE_SPAN(material, stage_start);
  const Color4 color = colorize(material, frame);
  HS_SB_STAGE_SPAN(color, stage_start);
  return color;
}
#endif

} // namespace Workbench

#endif // HS_ENABLE_SHADER_WORKBENCH
