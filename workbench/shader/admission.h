/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_SHADER_WORKBENCH

/**
 * @file admission.h
 * @brief Structural legality: which authored configurations are well formed, which coordinate bounds they must respect, and which transitions between two of them are admitted.
 */

#include "workbench/shader/config.h"
#include "workbench/shader/resources.h"
#include "workbench/shader/frame_state.h"
#include "workbench/shader/limits.h"

namespace Workbench {

// Declared ahead of first use: at namespace scope an unqualified
// call would otherwise bind to a same-named function at global scope.
HS_COLD_MEMBER inline constexpr bool safe_program_bounds(const Config &config);
HS_COLD_MEMBER inline constexpr bool valid_mobius(const MobiusParams &params);
inline constexpr bool valid_slot_enums(const Slots &slots);

// valid_config is the entry point the rest of the workbench calls, so it is
// written first; these are defined below it, or in resources.h.
constexpr bool valid_slot_enums(const Slots &slots);
constexpr bool valid_stage_tuple(const WarpStageSpec &spec,
                                 const WarpStageParams &params);
constexpr bool safe_program_bounds(const Config &config);
constexpr bool valid_mobius(const MobiusParams &params);
constexpr bool coefficient_in_range(const Complex &coefficient);
constexpr bool curl_pair_stable(const WarpStageSpec &spec,
                                const WarpStageParams &a,
                                const WarpStageParams &b);
constexpr bool resource_union_fits(const Config &from, const Config &to);

HS_COLD_MEMBER inline constexpr SourceTraits source_traits(Function function) {
  switch (function) {
  case Function::GRID:
    return {true, true};
  case Function::PRIMITIVE_LATTICE:
    return {true, true};
  case Function::TWIN_WAVE:
  case Function::RINGS:
  case Function::SPIRAL:
  case Function::NOISE_CONTOUR:
  case Function::NOISE_CONTOUR_SPHERE:
  case Function::SPHERICAL_RINGS:
  case Function::FRACTAL:
  case Function::TESSELLATION:
    return {false, false};
  }
  return {false, false};
}

HS_COLD_MEMBER inline constexpr Complex
source_cartesian_period(Function function, float lattice_cell_scale) {
  if (function != Function::PRIMITIVE_LATTICE || !(lattice_cell_scale > 0.0f))
    return {};
  const float period = 1.0f / lattice_cell_scale;
  return {period, period};
}

HS_COLD_MEMBER inline constexpr Complex
source_cartesian_period(const Config &config) {
  return source_cartesian_period(config.slots.function,
                                 config.params.source.lattice_cell_scale);
}

HS_COLD_MEMBER inline constexpr bool
affine_has_translation(const WarpStageSpec &spec,
                       const WarpStageParams &params) {
  return spec.kind == WarpStageKind::AFFINE_FRAME &&
         (params.translation_x != 0.0f || params.translation_y != 0.0f);
}

HS_COLD_MEMBER inline constexpr bool whole_affine_winding(float value) {
  if (!(value >= -4.0f && value <= 4.0f))
    return false;
  return value == static_cast<float>(static_cast<int>(value));
}

HS_COLD_MEMBER inline constexpr float hue_shift_amount_max(HueShiftMode mode) {
  return mode == HueShiftMode::NOISE ? HUE_NOISE_AMOUNT_MAX
                                     : HUE_SHIFT_AMOUNT_MAX;
}

HS_COLD_MEMBER inline constexpr bool
hue_shift_amount_in_range(const Config &config) {
  return config.params.color.hue_shift_amount >=
             -hue_shift_amount_max(config.slots.hue_shift) &&
         config.params.color.hue_shift_amount <=
             hue_shift_amount_max(config.slots.hue_shift);
}

HS_COLD_MEMBER inline constexpr bool
affine_translation_compatible(const Config &config) {
  const bool outer = affine_has_translation(config.slots.warp_program.outer,
                                            config.params.warp.outer);
  const bool inner = affine_has_translation(config.slots.warp_program.inner,
                                            config.params.warp.inner);
  if (!outer && !inner)
    return true;
  if (config.slots.function != Function::PRIMITIVE_LATTICE ||
      (outer && config.slots.warp_program.inner.kind != WarpStageKind::NONE) ||
      (config.slots.hue_shift == HueShiftMode::WARP_DISPLACEMENT &&
       config.params.color.hue_shift_amount != 0.0f))
    return false;
  return (!outer ||
          (whole_affine_winding(config.params.warp.outer.translation_x) &&
           whole_affine_winding(config.params.warp.outer.translation_y))) &&
         (!inner ||
          (whole_affine_winding(config.params.warp.inner.translation_x) &&
           whole_affine_winding(config.params.warp.inner.translation_y)));
}

/// Source periods spanned by one angular seam jump of a polar chart.
/// Primitive Lattice is periodic in its own cell scale and ignores the
/// pattern frequency, so its seam spans `2*pi*harmonic*cell_scale` cells.
HS_COLD_MEMBER inline constexpr float
polar_seam_periods(const RequestedConfig &config, const WarpStageSpec &polar) {
  const float harmonic = static_cast<float>(polar.polar_harmonic);
  if (config.slots.function == Function::PRIMITIVE_LATTICE)
    return TWO_PI_F * harmonic * config.params.source.lattice_cell_scale;
  return config.params.source.pattern_freq * harmonic;
}

HS_COLD_MEMBER inline constexpr bool
polar_source_compatible(const RequestedConfig &config,
                        const WarpStageSpec &polar) {
  const SourceTraits traits = source_traits(config.slots.function);
  if (!traits.y_periodic || !traits.polar_angle_compatible)
    return false;
  const float periods = polar_seam_periods(config, polar);
  return periods == static_cast<float>(static_cast<int>(periods));
}

HS_COLD_MEMBER inline constexpr bool
strict_seam_compatible(const Config &config) {
  if (!strict_projection(config.slots.projection))
    return true;
  return config.slots.function != Function::NOISE_CONTOUR &&
         !seam_sensitive_warp(config.slots.warp_program.outer.kind) &&
         !seam_sensitive_warp(config.slots.warp_program.inner.kind);
}

HS_COLD_MEMBER inline constexpr bool
valid_config(const RequestedConfig &candidate) {
  const Slots &slots = candidate.slots;
  if (!valid_slot_enums(slots))
    return false;
  if (is_sphere_source(slots.function) &&
      (slots.warp_program.outer.kind != WarpStageKind::NONE ||
       slots.warp_program.inner.kind != WarpStageKind::NONE))
    return false;
  if (slots.surface_lens == SurfaceLens::TANGENT_NOISE ||
      slots.warp_program.outer.kind == WarpStageKind::LEGACY_STEREO_NOISE ||
      slots.warp_program.inner.kind == WarpStageKind::LEGACY_STEREO_NOISE)
    return false;
  const bool outer_polar =
      slots.warp_program.outer.kind == WarpStageKind::POLAR_CHART;
  const bool inner_polar =
      slots.warp_program.inner.kind == WarpStageKind::POLAR_CHART;
  if ((outer_polar && slots.warp_program.inner.kind != WarpStageKind::NONE &&
       slots.warp_program.inner.kind != WarpStageKind::WAVE_SHEAR) ||
      (inner_polar && slots.warp_program.outer.kind != WarpStageKind::NONE))
    return false;
  if (inner_polar &&
      !polar_source_compatible(candidate, slots.warp_program.inner))
    return false;
  if (outer_polar &&
      slots.warp_program.inner.kind == WarpStageKind::WAVE_SHEAR) {
    const SourceTraits traits = source_traits(slots.function);
    if (!traits.y_periodic || !traits.polar_angle_compatible)
      return false;
  }
  if (outer_polar && slots.warp_program.inner.kind == WarpStageKind::NONE &&
      !polar_source_compatible(candidate, slots.warp_program.outer))
    return false;
  if (!affine_translation_compatible(candidate) ||
      !strict_seam_compatible(candidate) || !preset_in_ranges(candidate) ||
      !hue_shift_amount_in_range(candidate))
    return false;
  if (!valid_stage_tuple(slots.warp_program.outer,
                         candidate.params.warp.outer) ||
      !valid_stage_tuple(slots.warp_program.inner,
                         candidate.params.warp.inner) ||
      !safe_program_bounds(candidate))
    return false;
  if (slots.surface_lens == SurfaceLens::MOBIUS &&
      !valid_mobius(candidate.params.surface_lens.mobius))
    return false;
  const SurfaceNoiseParams &surface_noise = candidate.params.surface_noise;
  if (!enum_at_most(surface_noise.basis, NoiseBasis::RIDGED3) ||
      !enum_at_most(surface_noise.integrator,
                    SurfaceCurlIntegrator::MIDPOINT_2X) ||
      surface_noise.scale < LENS_NOISE_SCALE_MIN ||
      surface_noise.scale > LENS_NOISE_SCALE_MAX ||
      surface_noise.strength <
          (slots.surface_noise == SurfaceNoise::CURL ? -0.5f : 0.0f) ||
      surface_noise.strength > 0.5f || surface_noise.rate < NOISE_RATE_MIN ||
      surface_noise.rate > NOISE_RATE_MAX || surface_noise.direction < 0.0f ||
      surface_noise.direction > 1.0f)
    return false;
  return resource_union_fits(candidate, candidate);
}

HS_COLD_MEMBER inline constexpr bool
valid_warp_spec(const WarpStageSpec &spec) {
  return enum_at_most(spec.basis, NoiseBasis::RIDGED3) &&
         enum_at_most(spec.envelope, WarpEnvelope::EDGE_FADE) &&
         enum_at_most(spec.polar_mode, PolarMode::LOGARITHMIC) &&
         enum_at_most(spec.curl_integrator, CurlIntegrator::MIDPOINT_4) &&
         spec.polar_harmonic >= 1 && spec.polar_harmonic <= POLAR_HARMONIC_MAX;
}

HS_COLD_MEMBER inline constexpr float abs_value(float value) {
  return value < 0.0f ? -value : value;
}

HS_COLD_MEMBER inline constexpr int curl_intervals(CurlIntegrator integrator) {
  return integrator == CurlIntegrator::EULER_1      ? 1
         : integrator == CurlIntegrator::MIDPOINT_2 ? 2
                                                    : 4;
}

/** @brief Maximum component emitted by the bounded curl vector field. */
HS_COLD_MEMBER inline constexpr float curl_vector_component_bound(NoiseBasis) {
  return CURL_VECTOR_COMPONENT_MAX;
}

/**
 * @brief Largest curl-flow strength the stage stability inequality admits at
 *        a stage's live scale, basis, and integrator.
 * @details Solves `scale * |strength| * G / n <= 1/2` — the same inequality
 * `valid_stage_tuple` enforces — for `|strength|`, so the registered slider
 * spans exactly the admissible set instead of a range whose bulk is rejected.
 */
HS_COLD_MEMBER inline constexpr float
curl_strength_limit(const WarpStageSpec &spec, const WarpStageParams &params) {
  const float scale =
      params.scale > WARP_SCALE_MIN ? params.scale : WARP_SCALE_MIN;
  const float stable_limit =
      0.5f * static_cast<float>(curl_intervals(spec.curl_integrator)) /
      (scale * curl_vector_component_bound(spec.basis));
  return stable_limit < CURL_WARP_STRENGTH_MAX ? stable_limit
                                               : CURL_WARP_STRENGTH_MAX;
}

HS_COLD_MEMBER inline constexpr bool
valid_stage_tuple(const WarpStageSpec &spec, const WarpStageParams &params) {
  switch (spec.kind) {
  case WarpStageKind::NONE:
    return true;
  case WarpStageKind::LEGACY_STEREO_NOISE:
    return false;
  case WarpStageKind::AFFINE_FRAME:
    return params.translation_x >= -4.0f && params.translation_x <= 4.0f &&
           params.translation_y >= -4.0f && params.translation_y <= 4.0f &&
           params.rotation >= -TWO_PI_F && params.rotation <= TWO_PI_F &&
           params.scale_x >= 0.25f && params.scale_x <= 4.0f &&
           params.scale_y >= 0.25f && params.scale_y <= 4.0f &&
           params.shear >= -0.75f && params.shear <= 0.75f &&
           params.speed >= NOISE_SPEED_MIN && params.speed <= NOISE_SPEED_MAX;
  case WarpStageKind::WAVE_SHEAR:
    return params.strength >= -4.0f && params.strength <= 4.0f &&
           params.frequency >= 0.0f && params.frequency <= 64.0f &&
           params.speed >= NOISE_SPEED_MIN && params.speed <= NOISE_SPEED_MAX;
  case WarpStageKind::VORTEX:
    return params.radius >= 1.0f / 64.0f && params.radius <= 8.0f &&
           params.turns >= -4.0f && params.turns <= 4.0f &&
           params.center_orbit_radius >= 0.0f &&
           params.center_orbit_radius <= 4.0f &&
           params.speed >= NOISE_SPEED_MIN && params.speed <= NOISE_SPEED_MAX;
  case WarpStageKind::VECTOR_NOISE:
    return params.strength >= 0.0f &&
           params.strength <= VECTOR_WARP_STRENGTH_MAX &&
           params.scale >= 1.0f / 64.0f &&
           params.scale <= VECTOR_WARP_SCALE_MAX &&
           params.speed >= NOISE_SPEED_MIN && params.speed <= NOISE_SPEED_MAX;
  case WarpStageKind::CURL_FLOW:
    return params.strength >= -CURL_WARP_STRENGTH_MAX &&
           params.strength <= CURL_WARP_STRENGTH_MAX &&
           params.scale >= 1.0f / 64.0f &&
           params.scale <= CURL_WARP_SCALE_MAX &&
           params.speed >= NOISE_SPEED_MIN && params.speed <= NOISE_SPEED_MAX &&
           params.scale * abs_value(params.strength) *
                   curl_vector_component_bound(spec.basis) /
                   curl_intervals(spec.curl_integrator) <=
               0.5f;
  case WarpStageKind::MIRROR_TILE:
    return params.rotation >= 0.0f && params.rotation <= TWO_PI_F &&
           params.cell_x >= CELL_MIN && params.cell_x <= CELL_MAX &&
           params.cell_y >= CELL_MIN && params.cell_y <= CELL_MAX &&
           params.speed >= NOISE_SPEED_MIN && params.speed <= NOISE_SPEED_MAX;
  case WarpStageKind::POLAR_CHART:
    return params.radial_scale >= 1.0f / 64.0f &&
           params.radial_scale <= 16.0f && params.speed >= NOISE_SPEED_MIN &&
           params.speed <= NOISE_SPEED_MAX;
  }
  return false;
}

HS_COLD_MEMBER inline constexpr float
projection_coordinate_bound(const Config &config) {
  if (config.slots.surface_lens == SurfaceLens::MOBIUS)
    return STEREO_INF;
  if (config.slots.projection == Projection::STEREOGRAPHIC)
    return STEREO_INF;
  if (config.slots.projection == Projection::GNOMONIC)
    return 1.0f / GNOMONIC_AXIS_EPS;
  if (strict_projection(config.slots.projection))
    return 16.0f * config.params.projection.coordinate_scale;
  return 4.0f;
}

HS_COLD_MEMBER inline constexpr float
stage_coordinate_bound(const WarpStageSpec &spec, const WarpStageParams &params,
                       float input_bound, const Complex &source_period) {
  switch (spec.kind) {
  case WarpStageKind::NONE:
    return input_bound;
  case WarpStageKind::LEGACY_STEREO_NOISE:
    return WARP_COORD_LIMIT + 1.0f;
  case WarpStageKind::AFFINE_FRAME: {
    const float rotated = 1.414214f * input_bound;
    const float x_bound = rotated / params.scale_x +
                          abs_value(params.shear) * rotated / params.scale_y +
                          abs_value(params.translation_x * source_period.re);
    const float y_bound = rotated / params.scale_y +
                          abs_value(params.translation_y * source_period.im);
    return x_bound > y_bound ? x_bound : y_bound;
  }
  case WarpStageKind::WAVE_SHEAR:
    return input_bound + abs_value(params.strength);
  case WarpStageKind::VORTEX: {
    const float center_x =
        abs_value(params.center_x) + params.center_orbit_radius;
    const float center_y =
        abs_value(params.center_y) + params.center_orbit_radius;
    const float center_bound = center_x > center_y ? center_x : center_y;
    return 1.414214f * (input_bound + center_bound) + center_bound;
  }
  case WarpStageKind::VECTOR_NOISE:
    return input_bound + 1.414214f * params.strength;
  case WarpStageKind::CURL_FLOW:
    return input_bound + 1.414214f * abs_value(params.strength) * params.scale *
                             curl_vector_component_bound(spec.basis);
  case WarpStageKind::MIRROR_TILE:
    return 1.414214f * (params.cell_x + params.cell_y);
  case WarpStageKind::POLAR_CHART: {
    const float radial =
        params.radial_scale *
            (spec.polar_mode == PolarMode::LOGARITHMIC ? 12.0f : input_bound) +
        TWO_PI_F;
    return radial > 17.0f * PI_F ? radial : 17.0f * PI_F;
  }
  }
  return WARP_COORD_LIMIT + 1.0f;
}

HS_COLD_MEMBER inline constexpr bool safe_program_bounds(const Config &config) {
  float bound = projection_coordinate_bound(config);
  const Complex source_period = source_cartesian_period(config);
  const WarpStageSpec stages[] = {config.slots.warp_program.outer,
                                  config.slots.warp_program.inner};
  const WarpStageParams params[] = {config.params.warp.outer,
                                    config.params.warp.inner};
  for (size_t index = 0; index < 2; ++index) {
    if ((stages[index].kind == WarpStageKind::VECTOR_NOISE ||
         stages[index].kind == WarpStageKind::CURL_FLOW) &&
        params[index].scale * (bound + 100.0f) > NOISE_LATTICE_LIMIT)
      return false;
    bound = stage_coordinate_bound(stages[index], params[index], bound,
                                   source_period);
    if (bound > WARP_COORD_LIMIT)
      return false;
  }
  if (config.slots.function == Function::NOISE_CONTOUR &&
      config.params.source.noise_scale * bound > NOISE_LATTICE_LIMIT)
    return false;
  return true;
}

HS_COLD_MEMBER inline constexpr bool valid_mobius(const MobiusParams &params) {
  const float ad_re = params.a.re * params.d.re - params.a.im * params.d.im;
  const float ad_im = params.a.re * params.d.im + params.a.im * params.d.re;
  const float bc_re = params.b.re * params.c.re - params.b.im * params.c.im;
  const float bc_im = params.b.re * params.c.im + params.b.im * params.c.re;
  const float det_re = ad_re - bc_re;
  const float det_im = ad_im - bc_im;
  return coefficient_in_range(params.a) && coefficient_in_range(params.b) &&
         coefficient_in_range(params.c) && coefficient_in_range(params.d) &&
         det_re * det_re + det_im * det_im >= 1e-6f;
}

HS_COLD_MEMBER inline constexpr bool
coefficient_in_range(const Complex &coefficient) {
  return coefficient.re >= -8.0f && coefficient.re <= 8.0f &&
         coefficient.im >= -8.0f && coefficient.im <= 8.0f;
}

HS_COLD_MEMBER inline constexpr float max_value(float a, float b) {
  return a > b ? a : b;
}

HS_COLD_MEMBER inline constexpr float min_value(float a, float b) {
  return a < b ? a : b;
}

HS_COLD_MEMBER inline constexpr float max_abs_value(float a, float b) {
  return max_value(abs_value(a), abs_value(b));
}

HS_COLD_MEMBER inline constexpr void
maximize_stage_path(WarpStageParams &out, const WarpStageParams &a,
                    const WarpStageParams &b) {
  out.translation_x = max_abs_value(a.translation_x, b.translation_x);
  out.translation_y = max_abs_value(a.translation_y, b.translation_y);
  out.scale_x = min_value(a.scale_x, b.scale_x);
  out.scale_y = min_value(a.scale_y, b.scale_y);
  out.shear = max_abs_value(a.shear, b.shear);
  out.strength = max_abs_value(a.strength, b.strength);
  out.scale = max_value(a.scale, b.scale);
  out.center_x = max_abs_value(a.center_x, b.center_x);
  out.center_y = max_abs_value(a.center_y, b.center_y);
  out.center_orbit_radius =
      max_value(a.center_orbit_radius, b.center_orbit_radius);
  out.cell_x = max_value(a.cell_x, b.cell_x);
  out.cell_y = max_value(a.cell_y, b.cell_y);
  out.radial_scale = max_value(a.radial_scale, b.radial_scale);
}

HS_COLD_MEMBER inline constexpr bool safe_program_path(const Config &from,
                                                       const Config &to) {
  Config worst = from;
  worst.params.projection.coordinate_scale =
      max_value(from.params.projection.coordinate_scale,
                to.params.projection.coordinate_scale);
  worst.params.source.noise_scale =
      max_value(from.params.source.noise_scale, to.params.source.noise_scale);
  worst.params.source.lattice_cell_scale =
      min_value(from.params.source.lattice_cell_scale,
                to.params.source.lattice_cell_scale);
  maximize_stage_path(worst.params.warp.outer, from.params.warp.outer,
                      to.params.warp.outer);
  maximize_stage_path(worst.params.warp.inner, from.params.warp.inner,
                      to.params.warp.inner);
  return safe_program_bounds(worst);
}

HS_COLD_MEMBER inline constexpr bool polar_pair_stable(const Config &from,
                                                       const Config &to) {
  const WarpStageSpec &outer = from.slots.warp_program.outer;
  const WarpStageSpec &inner = from.slots.warp_program.inner;
  return (outer.kind != WarpStageKind::POLAR_CHART ||
          polar_seam_periods(from, outer) == polar_seam_periods(to, outer)) &&
         (inner.kind != WarpStageKind::POLAR_CHART ||
          polar_seam_periods(from, inner) == polar_seam_periods(to, inner));
}

HS_COLD_MEMBER inline constexpr bool
affine_winding_pair_stable(const WarpStageSpec &spec, const WarpStageParams &a,
                           const WarpStageParams &b) {
  return spec.kind != WarpStageKind::AFFINE_FRAME ||
         (a.translation_x == b.translation_x &&
          a.translation_y == b.translation_y);
}

HS_COLD_MEMBER inline constexpr bool
stable_parameter_path_admitted(const Config &from, const Config &to) {
  return curl_pair_stable(from.slots.warp_program.outer, from.params.warp.outer,
                          to.params.warp.outer) &&
         curl_pair_stable(from.slots.warp_program.inner, from.params.warp.inner,
                          to.params.warp.inner) &&
         affine_winding_pair_stable(from.slots.warp_program.outer,
                                    from.params.warp.outer,
                                    to.params.warp.outer) &&
         affine_winding_pair_stable(from.slots.warp_program.inner,
                                    from.params.warp.inner,
                                    to.params.warp.inner) &&
         polar_pair_stable(from, to) && safe_program_path(from, to);
}

HS_COLD_MEMBER inline constexpr bool same_parameter_topology(const Config &from,
                                                             const Config &to) {
  Slots from_slots = from.slots;
  Slots to_slots = to.slots;
  from_slots.palette_mapping = PaletteMapping::LINEAR;
  to_slots.palette_mapping = PaletteMapping::LINEAR;
  return from_slots == to_slots &&
         from.params.source.noise_basis == to.params.source.noise_basis &&
         from.params.source.noise_seed == to.params.source.noise_seed &&
         from.params.value.band_count == to.params.value.band_count &&
         from.params.surface_noise.basis == to.params.surface_noise.basis &&
         from.params.surface_noise.integrator ==
             to.params.surface_noise.integrator &&
         from.params.surface_noise.seed == to.params.surface_noise.seed &&
         from.params.surface_lens.mobius.a.re ==
             to.params.surface_lens.mobius.a.re &&
         from.params.surface_lens.mobius.a.im ==
             to.params.surface_lens.mobius.a.im &&
         from.params.surface_lens.mobius.b.re ==
             to.params.surface_lens.mobius.b.re &&
         from.params.surface_lens.mobius.b.im ==
             to.params.surface_lens.mobius.b.im &&
         from.params.surface_lens.mobius.c.re ==
             to.params.surface_lens.mobius.c.re &&
         from.params.surface_lens.mobius.c.im ==
             to.params.surface_lens.mobius.c.im &&
         from.params.surface_lens.mobius.d.re ==
             to.params.surface_lens.mobius.d.re &&
         from.params.surface_lens.mobius.d.im ==
             to.params.surface_lens.mobius.d.im;
}

HS_COLD_MEMBER inline constexpr bool stable_topology(const Config &from,
                                                     const Config &to) {
  return valid_config(from) && valid_config(to) &&
         same_parameter_topology(from, to) &&
         stable_parameter_path_admitted(from, to);
}

HS_COLD_MEMBER inline constexpr bool
curl_pair_stable(const WarpStageSpec &spec, const WarpStageParams &a,
                 const WarpStageParams &b) {
  if (spec.kind != WarpStageKind::CURL_FLOW)
    return true;
  const float scale = a.scale > b.scale ? a.scale : b.scale;
  const float a_distance = abs_value(a.strength);
  const float b_distance = abs_value(b.strength);
  const float distance = a_distance > b_distance ? a_distance : b_distance;
  return scale * distance * curl_vector_component_bound(spec.basis) /
             curl_intervals(spec.curl_integrator) <=
         0.5f;
}

/** @brief Reports whether both transition endpoints are admissible holds. */
HS_COLD_MEMBER inline constexpr bool transition_admitted(const Config &from,
                                                         const Config &to) {
  return valid_config(from) && valid_config(to);
}

/** @brief Tests every slot against the last enumerator its schema admits. */
inline constexpr bool valid_slot_enums(const Slots &slots) {
  return enum_at_most(slots.function, Function::TESSELLATION) &&
         enum_at_most(slots.projection, Projection::EQUIRECTANGULAR) &&
         enum_at_most(slots.projection_frame,
                      ProjectionFramePolicy::SPIN_WANDER) &&
         enum_at_most(slots.surface_lens,
                      SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM) &&
         enum_at_most(slots.surface_noise, SurfaceNoise::CURL) &&
         enum_at_most(slots.surface_noise_placement,
                      SurfaceNoisePlacement::AFTER_LENS) &&
         enum_at_most(slots.warp_program.outer.kind,
                      WarpStageKind::POLAR_CHART) &&
         enum_at_most(slots.warp_program.inner.kind,
                      WarpStageKind::POLAR_CHART) &&
         valid_warp_spec(slots.warp_program.outer) &&
         valid_warp_spec(slots.warp_program.inner) &&
         enum_at_most(slots.signal_weight, SignalWeight::PROJECTION) &&
         enum_at_most(slots.value_transfer, ValueTransfer::SMOOTH_BANDS) &&
         enum_at_most(slots.coverage, CoveragePolicy::PROJECTION_WEIGHT) &&
         enum_at_most(slots.palette, PaletteMode::ANALOGOUS) &&
         enum_at_most(slots.palette_mapping, PaletteMapping::REVERSE) &&
         enum_at_most(slots.brightness_envelope,
                      BrightnessEnvelope::DESCENDING) &&
         enum_at_most(slots.hue_shift, HueShiftMode::WARP_DISPLACEMENT) &&
         enum_at_most(slots.peirce_layout, PeirceLayout::VERTICAL) &&
         enum_at_most(slots.airocean_layout, AiroceanLayout::HORIZONTAL) &&
         enum_at_most(slots.bonne_hemisphere, BonneHemisphere::SOUTH) &&
         enum_at_most(slots.gnomonic_hemisphere,
                      GnomonicHemispherePolicy::BACK_HEMISPHERE);
}

inline constexpr bool valid_snapshot_config(const Config &config) {
  return valid_slot_enums(config.slots) &&
         enum_at_most(config.params.color.palette_mapping,
                      Pullback::Color::PaletteMapping::REVERSE) &&
         preset_in_ranges(config) && hue_shift_amount_in_range(config);
}

} // namespace Workbench

#endif // HS_ENABLE_SHADER_WORKBENCH
