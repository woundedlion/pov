/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_SHADER_WORKBENCH

/**
 * @file limits.h
 * @brief Parameter domain bounds: the slider ranges every authored Config is held to, and the predicates that check one against them.
 */

#include "workbench/shader/config.h"

namespace Workbench {

inline constexpr int POLAR_HARMONIC_MAX = 16;
inline constexpr int BAND_COUNT_MAX = 32;
inline constexpr float WARP_SCALE_MIN = 1.0f / 64.0f;
inline constexpr float WARP_SCALE_MAX = 100.0f;
inline constexpr float WARP_STRENGTH_MIN = -4.0f;
inline constexpr float WARP_STRENGTH_MAX = 30.0f;
inline constexpr float VECTOR_WARP_SCALE_MAX = 4.0f;
inline constexpr float VECTOR_WARP_STRENGTH_MAX = 1.0f;
inline constexpr float CURL_WARP_SCALE_MAX = 2.0f;
inline constexpr float CURL_WARP_STRENGTH_MAX = 1.0f;
inline constexpr float CURL_VECTOR_COMPONENT_MAX = 4.0f;
inline constexpr float WARP_SPEED_MIN = -1.0f / 64.0f;
inline constexpr float WARP_SPEED_MAX = 1.0f;
inline constexpr float PATTERN_FREQ_MIN = 0.1f;
inline constexpr float PATTERN_FREQ_MAX = 20.0f;
inline constexpr float GRID_PATTERN_FREQ_MIN = 0.01f;
inline constexpr float GRID_PATTERN_FREQ_MAX = 64.0f;
inline constexpr float SPEED_MIN = -0.5f, SPEED_MAX = 5.0f;
inline constexpr float COMPLEXITY_MIN = 0.0f, COMPLEXITY_MAX = 3.0f;
inline constexpr float PATTERN_MIX_MIN = 0.0f, PATTERN_MIX_MAX = 1.0f;
inline constexpr float PHASE2_RATE_MIN = 0.0f, PHASE2_RATE_MAX = 2.0f;
inline constexpr float SINGULARITY_FADE_MIN = 1.0f,
                       SINGULARITY_FADE_MAX = 20.0f;
inline constexpr float SPIN_RATE_MIN = 0.0f, SPIN_RATE_MAX = 0.05f;
inline constexpr float WANDER_MIN = 0.0f, WANDER_MAX = 1.0f;
inline constexpr float HUE_SHIFT_AMOUNT_MAX = 4.0f;
inline constexpr float HUE_NOISE_AMOUNT_MAX = 1.0f;
inline constexpr float HUE_NOISE_SCALE_MIN = 1.0f / 64.0f;
inline constexpr float HUE_NOISE_SCALE_MAX = 8.0f;
inline constexpr float HUE_NOISE_SPEED_MAX = 0.001f;
inline constexpr float PALETTE_CHROMA_MIN = 0.0f;
inline constexpr float PALETTE_CHROMA_MAX = 1.0f;
inline constexpr float MAPPING_FREQUENCY_MIN = 1.0f;
inline constexpr float MAPPING_FREQUENCY_MAX = 32.0f;
inline constexpr float MAPPING_PHASE_MIN = -1.0f;
inline constexpr float MAPPING_PHASE_MAX = 1.0f;
inline constexpr float PHASE_OSCILLATION_DEPTH_MIN = 0.0f;
inline constexpr float PHASE_OSCILLATION_DEPTH_MAX = 1.0f;
inline constexpr float PHASE_OSCILLATION_SPEED_MAX = 0.01f;
inline constexpr float BRIGHTNESS_GAIN_MIN = 0.0f;
inline constexpr float BRIGHTNESS_GAIN_MAX = 1.0f;
inline constexpr float VALUE_OPACITY_MIN = 0.0f;
inline constexpr float VALUE_OPACITY_MAX = 1.0f;
inline constexpr float WAVE_SPIN_MIN = -0.05f, WAVE_SPIN_MAX = 0.05f;
inline constexpr float SOURCE_NOISE_SCALE_MIN = 0.0f;
inline constexpr float SOURCE_NOISE_SCALE_MAX = 2.0f;
inline constexpr float SOURCE_NOISE_RATE_MIN = -1.0f / 1024.0f;
inline constexpr float SOURCE_NOISE_RATE_MAX = 1.0f / 1024.0f;
inline constexpr float LENS_NOISE_SCALE_MIN = 1.0f / 64.0f;
inline constexpr float LENS_NOISE_SCALE_MAX = 8.0f;
inline constexpr float NOISE_RATE_MIN = -1.0f / 64.0f;
inline constexpr float NOISE_RATE_MAX = 1.0f / 64.0f;
inline constexpr float NOISE_SPEED_MIN = NOISE_RATE_MIN;
inline constexpr float NOISE_SPEED_MAX = NOISE_RATE_MAX;
inline constexpr float CELL_MIN = 1.0f / 64.0f;
inline constexpr float CELL_MAX = 8.0f;
inline constexpr float SOFTNESS_MIN = 1.0f / 1024.0f;

inline constexpr float lens_domain_linear_scale(SurfaceLens lens) {
  switch (lens) {
  case SurfaceLens::KALEIDOSCOPE:
    return 1.0f;
  case SurfaceLens::KALEIDOSCOPE_TETRAHEDRAL:
    return 0.20412415f;
  case SurfaceLens::KALEIDOSCOPE_OCTAHEDRAL:
    return 0.14433757f;
  case SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL:
    return 0.09128709f;
  case SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM:
    return 0.28867513f;
  case SurfaceLens::KALEIDOSCOPE_SQUARE_PRISM:
    return 0.25f;
  case SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM:
    return 0.22360680f;
  case SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM:
    return 0.20412415f;
  case SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM:
    return 0.17677670f;
  case SurfaceLens::NONE:
  case SurfaceLens::GLITCH:
  case SurfaceLens::TWIST:
  case SurfaceLens::MOBIUS:
  case SurfaceLens::TANGENT_NOISE:
    return 1.0f;
  }
  __builtin_unreachable();
}

inline constexpr float pattern_freq_max(Function function) {
  return function == Function::GRID ? GRID_PATTERN_FREQ_MAX : PATTERN_FREQ_MAX;
}

inline constexpr float pattern_freq_min(Function function) {
  return function == Function::GRID ? GRID_PATTERN_FREQ_MIN : PATTERN_FREQ_MIN;
}

template <typename Enum>
HS_COLD_MEMBER inline constexpr bool enum_at_most(Enum value, Enum last) {
  return static_cast<uint8_t>(value) <= static_cast<uint8_t>(last);
}

HS_COLD_MEMBER inline constexpr bool
warp_stage_params_in_ranges(const WarpStageParams &params) {
  static_assert(sizeof(WarpStageParams) == 100,
                "WarpStageParams field set changed - update the range check");
  return params.scale >= WARP_SCALE_MIN && params.scale <= WARP_SCALE_MAX &&
         params.strength >= WARP_STRENGTH_MIN &&
         params.strength <= WARP_STRENGTH_MAX &&
         params.speed >= WARP_SPEED_MIN && params.speed <= WARP_SPEED_MAX &&
         params.translation_x >= -4.0f && params.translation_x <= 4.0f &&
         params.translation_y >= -4.0f && params.translation_y <= 4.0f &&
         params.rotation >= -TWO_PI_F && params.rotation <= TWO_PI_F &&
         params.scale_x >= 0.25f && params.scale_x <= 4.0f &&
         params.scale_y >= 0.25f && params.scale_y <= 4.0f &&
         params.shear >= -0.75f && params.shear <= 0.75f &&
         params.frequency >= 0.0f && params.frequency <= 64.0f &&
         params.field_angle >= 0.0f && params.field_angle <= TWO_PI_F &&
         params.center_x >= -4.0f && params.center_x <= 4.0f &&
         params.center_y >= -4.0f && params.center_y <= 4.0f &&
         params.radius >= 1.0f / 64.0f && params.radius <= 8.0f &&
         params.turns >= -4.0f && params.turns <= 4.0f &&
         params.center_orbit_radius >= 0.0f &&
         params.center_orbit_radius <= 4.0f && params.vector_angle >= 0.0f &&
         params.vector_angle <= TWO_PI_F && params.cell_x >= CELL_MIN &&
         params.cell_x <= CELL_MAX && params.cell_y >= CELL_MIN &&
         params.cell_y <= CELL_MAX && params.offset_x >= -8.0f &&
         params.offset_x <= 8.0f && params.offset_y >= -8.0f &&
         params.offset_y <= 8.0f && params.radial_scale >= 1.0f / 64.0f &&
         params.radial_scale <= 16.0f && params.radial_phase >= 0.0f &&
         params.radial_phase <= TWO_PI_F && params.angular_phase >= 0.0f &&
         params.angular_phase <= TWO_PI_F &&
         params.edge_width >= SOFTNESS_MIN && params.edge_width <= 0.5f;
}

HS_COLD_MEMBER inline constexpr bool preset_in_ranges(const Config &config) {
  const Params &p = config.params;
  return warp_stage_params_in_ranges(p.warp.outer) &&
         warp_stage_params_in_ranges(p.warp.inner) &&
         p.source.pattern_freq >= pattern_freq_min(config.slots.function) &&
         p.source.pattern_freq <= pattern_freq_max(config.slots.function) &&
         p.source.speed >= SPEED_MIN && p.source.speed <= SPEED_MAX &&
         p.source.complexity >= COMPLEXITY_MIN &&
         p.source.complexity <= COMPLEXITY_MAX &&
         p.source.pattern_mix >= PATTERN_MIX_MIN &&
         p.source.pattern_mix <= PATTERN_MIX_MAX &&
         p.source.secondary_rate >= PHASE2_RATE_MIN &&
         p.source.secondary_rate <= PHASE2_RATE_MAX &&
         p.source.angle_rate >= WAVE_SPIN_MIN &&
         p.source.angle_rate <= WAVE_SPIN_MAX &&
         p.source.noise_scale >= SOURCE_NOISE_SCALE_MIN &&
         p.source.noise_scale <= SOURCE_NOISE_SCALE_MAX &&
         p.source.noise_contrast >= 0.0f && p.source.noise_contrast <= 8.0f &&
         p.source.noise_time_rate >= SOURCE_NOISE_RATE_MIN &&
         p.source.noise_time_rate <= SOURCE_NOISE_RATE_MAX &&
         p.source.lattice_cell_scale >= CELL_MIN &&
         p.source.lattice_cell_scale <= CELL_MAX &&
         p.source.lattice_shape_blend >= 0.0f &&
         p.source.lattice_shape_blend <= 1.0f &&
         p.source.lattice_softness >= SOFTNESS_MIN &&
         p.source.lattice_softness <= 1.0f &&
         p.source.lattice_radius >= 1.0f / 64.0f &&
         p.source.lattice_radius <= 0.49f &&
         enum_at_most(p.source.noise_basis, NoiseBasis::RIDGED3) &&
         p.source.ring_count >= 1 && p.source.ring_count <= 32 &&
         p.source.ring_thickness >= 1.0f / 512.0f &&
         p.source.ring_thickness <= 0.5f &&
         p.source.ring_softness >= SOFTNESS_MIN &&
         p.source.ring_softness <= 0.25f && p.source.ring_wander >= 0.0f &&
         p.source.ring_wander <= 1.0f &&
         p.source.fractal_scale >= 1.0f / 64.0f &&
         p.source.fractal_scale <= 8.0f && p.source.fractal_iterations >= 2 &&
         p.source.fractal_iterations <= 16 && p.source.julia_mix >= 0.0f &&
         p.source.julia_mix <= 1.0f && p.source.julia_real >= -1.5f &&
         p.source.julia_real <= 1.5f && p.source.julia_imaginary >= -1.5f &&
         p.source.julia_imaginary <= 1.5f &&
         p.source.fractal_contours >= 0.0f &&
         p.source.fractal_contours <= 16.0f &&
         p.source.tessellation_cell_scale >= 1.0f / 64.0f &&
         p.source.tessellation_cell_scale <= 8.0f &&
         p.source.tessellation_line_thickness >= SOFTNESS_MIN &&
         p.source.tessellation_line_thickness <= 0.25f &&
         p.source.tessellation_line_softness >= SOFTNESS_MIN &&
         p.source.tessellation_line_softness <= 0.25f &&
         enum_at_most(p.source.tessellation_kind,
                      Pullback::Source::TessellationKind::HEXAGONAL) &&
         p.projection.singularity_fade >= SINGULARITY_FADE_MIN &&
         p.projection.singularity_fade <= SINGULARITY_FADE_MAX &&
         p.projection.spin_rate >= SPIN_RATE_MIN &&
         p.projection.spin_rate <= SPIN_RATE_MAX &&
         p.projection.wander >= WANDER_MIN &&
         p.projection.wander <= WANDER_MAX &&
         p.projection.central_meridian >= 0.0f &&
         p.projection.central_meridian <= TWO_PI_F &&
         p.projection.coordinate_scale >= 0.25f &&
         p.projection.coordinate_scale <= 4.0f &&
         p.projection.bonne_standard_parallel >= 1e-3f &&
         p.projection.bonne_standard_parallel <= 0.5f * PI_F &&
         p.projection.layout_scroll >= -1.0f &&
         p.projection.layout_scroll <= 1.0f &&
         p.outer_camera.wander >= WANDER_MIN &&
         p.outer_camera.wander <= WANDER_MAX &&
         p.surface_noise.scale >= LENS_NOISE_SCALE_MIN &&
         p.surface_noise.scale <= LENS_NOISE_SCALE_MAX &&
         p.surface_noise.strength >= -0.5f &&
         p.surface_noise.strength <= 0.5f &&
         p.surface_noise.rate >= NOISE_RATE_MIN &&
         p.surface_noise.rate <= NOISE_RATE_MAX &&
         p.surface_noise.direction >= 0.0f &&
         p.surface_noise.direction <= 1.0f &&
         enum_at_most(p.surface_noise.basis, NoiseBasis::RIDGED3) &&
         enum_at_most(p.surface_noise.integrator,
                      SurfaceCurlIntegrator::MIDPOINT_2X) &&
         p.value.iso_level >= 0.0f && p.value.iso_level <= 1.0f &&
         p.value.iso_width >= SOFTNESS_MIN && p.value.iso_width <= 0.5f &&
         p.value.band_count >= 1 && p.value.band_count <= BAND_COUNT_MAX &&
         p.value.band_phase >= 0.0f && p.value.band_phase <= TWO_PI_F &&
         p.value.cutout_threshold >= 0.0f && p.value.cutout_threshold <= 1.0f &&
         p.value.cutout_softness >= SOFTNESS_MIN &&
         p.value.cutout_softness <= 0.5f && p.value.edge_width >= 0.0f &&
         p.value.edge_width <= 0.5f &&
         p.color.hue_shift_amount >= -HUE_SHIFT_AMOUNT_MAX &&
         p.color.hue_shift_amount <= HUE_SHIFT_AMOUNT_MAX &&
         p.color.hue_noise_scale >= HUE_NOISE_SCALE_MIN &&
         p.color.hue_noise_scale <= HUE_NOISE_SCALE_MAX &&
         p.color.hue_noise_speed >= -HUE_NOISE_SPEED_MAX &&
         p.color.hue_noise_speed <= HUE_NOISE_SPEED_MAX &&
         p.color.palette_chroma >= PALETTE_CHROMA_MIN &&
         p.color.palette_chroma <= PALETTE_CHROMA_MAX &&
         p.color.mapping_frequency >= MAPPING_FREQUENCY_MIN &&
         p.color.mapping_frequency <= MAPPING_FREQUENCY_MAX &&
         p.color.mapping_phase >= MAPPING_PHASE_MIN &&
         p.color.mapping_phase <= MAPPING_PHASE_MAX &&
         p.color.phase_oscillation_depth >= PHASE_OSCILLATION_DEPTH_MIN &&
         p.color.phase_oscillation_depth <= PHASE_OSCILLATION_DEPTH_MAX &&
         p.color.phase_oscillation_speed >= -PHASE_OSCILLATION_SPEED_MAX &&
         p.color.phase_oscillation_speed <= PHASE_OSCILLATION_SPEED_MAX &&
         p.color.brightness_bottom >= BRIGHTNESS_GAIN_MIN &&
         p.color.brightness_bottom <= BRIGHTNESS_GAIN_MAX &&
         p.color.brightness_top >= BRIGHTNESS_GAIN_MIN &&
         p.color.brightness_top <= BRIGHTNESS_GAIN_MAX &&
         p.color.opacity_low >= VALUE_OPACITY_MIN &&
         p.color.opacity_low <= VALUE_OPACITY_MAX &&
         p.color.opacity_high >= VALUE_OPACITY_MIN &&
         p.color.opacity_high <= VALUE_OPACITY_MAX;
}

inline constexpr uint8_t BOUNDARY_CUT =
    projections::projection_boundary(projections::ProjectionBoundary::CUT);
inline constexpr uint8_t BOUNDARY_SINGULAR =
    projections::projection_boundary(projections::ProjectionBoundary::SINGULAR);
inline constexpr uint8_t PROJECTION_FLAG_FOLDED = 1U << 0;
inline constexpr float GNOMONIC_AXIS_EPS = 1e-3f;
inline constexpr float WARP_COORD_LIMIT = 65536.0f;
inline constexpr float NOISE_LATTICE_LIMIT = 1048576.0f;

} // namespace Workbench

#endif // HS_ENABLE_SHADER_WORKBENCH
