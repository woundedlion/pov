/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_SHADER_WORKBENCH

#include <cstdarg>
#include <cstdio>
#include <span>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

/**
 * @file shader_host.h
 * @brief Typed pullback sphere shader with composable projection and material stages.
 */

#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"
#include "core/math/interpolate.h"
#include "core/math/lenses.h"
#include "core/math/noise_field.h"
#include "core/math/projections.h"
#include "core/math/stereographic.h"
#include "core/render/pullback.h"
#include "core/render/pullback/runtime_seeds.h"
#include "workbench/shader/admission.h"
#include "workbench/shader/config.h"
#include "workbench/shader/frame_state.h"
#include "workbench/shader/limits.h"
#include "workbench/shader/options.h"
#include "workbench/shader/presets.h"
#include "workbench/shader/resources.h"

namespace hs_test {
namespace shader_workbench_tests {
struct ShaderWorkbenchWhiteBox;
} // namespace shader_workbench_tests
} // namespace hs_test

#define HS_SHADER_WORKBENCH_CONFIG_FIELDS(X)                                   \
  X(SLOTS_FUNCTION, slots.function)                                            \
  X(SLOTS_PROJECTION, slots.projection)                                        \
  X(SLOTS_PROJECTION_FRAME, slots.projection_frame)                            \
  X(SLOTS_SURFACE_LENS, slots.surface_lens)                                    \
  X(SLOTS_WARP_OUTER_KIND, slots.warp_program.outer.kind)                      \
  X(SLOTS_WARP_OUTER_BASIS, slots.warp_program.outer.basis)                    \
  X(SLOTS_WARP_OUTER_ENVELOPE, slots.warp_program.outer.envelope)              \
  X(SLOTS_WARP_OUTER_POLAR_MODE, slots.warp_program.outer.polar_mode)          \
  X(SLOTS_WARP_OUTER_CURL_INTEGRATOR,                                          \
    slots.warp_program.outer.curl_integrator)                                  \
  X(SLOTS_WARP_OUTER_POLAR_HARMONIC, slots.warp_program.outer.polar_harmonic)  \
  X(SLOTS_WARP_OUTER_SEED, slots.warp_program.outer.seed)                      \
  X(SLOTS_WARP_INNER_KIND, slots.warp_program.inner.kind)                      \
  X(SLOTS_WARP_INNER_BASIS, slots.warp_program.inner.basis)                    \
  X(SLOTS_WARP_INNER_ENVELOPE, slots.warp_program.inner.envelope)              \
  X(SLOTS_WARP_INNER_POLAR_MODE, slots.warp_program.inner.polar_mode)          \
  X(SLOTS_WARP_INNER_CURL_INTEGRATOR,                                          \
    slots.warp_program.inner.curl_integrator)                                  \
  X(SLOTS_WARP_INNER_POLAR_HARMONIC, slots.warp_program.inner.polar_harmonic)  \
  X(SLOTS_WARP_INNER_SEED, slots.warp_program.inner.seed)                      \
  X(SLOTS_SIGNAL_WEIGHT, slots.signal_weight)                                  \
  X(SLOTS_VALUE_TRANSFER, slots.value_transfer)                                \
  X(SLOTS_COVERAGE, slots.coverage)                                            \
  X(SLOTS_PALETTE, slots.palette)                                              \
  X(SLOTS_PALETTE_MAPPING, slots.palette_mapping)                              \
  X(SLOTS_BRIGHTNESS_ENVELOPE, slots.brightness_envelope)                      \
  X(SLOTS_HUE_SHIFT, slots.hue_shift)                                          \
  X(SLOTS_PEIRCE_LAYOUT, slots.peirce_layout)                                  \
  X(SLOTS_AIROCEAN_LAYOUT, slots.airocean_layout)                              \
  X(SLOTS_BONNE_HEMISPHERE, slots.bonne_hemisphere)                            \
  X(SLOTS_GNOMONIC_HEMISPHERE, slots.gnomonic_hemisphere)                      \
  X(SLOTS_SURFACE_NOISE, slots.surface_noise)                                  \
  X(SLOTS_SURFACE_NOISE_PLACEMENT, slots.surface_noise_placement)              \
  X(SOURCE_PATTERN_FREQ, params.source.pattern_freq)                           \
  X(SOURCE_SPEED, params.source.speed)                                         \
  X(SOURCE_COMPLEXITY, params.source.complexity)                               \
  X(SOURCE_PATTERN_MIX, params.source.pattern_mix)                             \
  X(SOURCE_SECONDARY_RATE, params.source.secondary_rate)                       \
  X(SOURCE_ANGLE_RATE, params.source.angle_rate)                               \
  X(SOURCE_NOISE_SCALE, params.source.noise_scale)                             \
  X(SOURCE_NOISE_CONTRAST, params.source.noise_contrast)                       \
  X(SOURCE_NOISE_RATE, params.source.noise_time_rate)                          \
  X(SOURCE_LATTICE_CELL_SCALE, params.source.lattice_cell_scale)               \
  X(SOURCE_LATTICE_SHAPE_BLEND, params.source.lattice_shape_blend)             \
  X(SOURCE_LATTICE_SOFTNESS, params.source.lattice_softness)                   \
  X(SOURCE_LATTICE_RADIUS, params.source.lattice_radius)                       \
  X(SOURCE_NOISE_BASIS, params.source.noise_basis)                             \
  X(SOURCE_NOISE_SEED, params.source.noise_seed)                               \
  X(WARP_OUTER_SCALE, params.warp.outer.scale)                                 \
  X(WARP_OUTER_STRENGTH, params.warp.outer.strength)                           \
  X(WARP_OUTER_SPEED, params.warp.outer.speed)                                 \
  X(WARP_OUTER_TRANSLATION_X, params.warp.outer.translation_x)                 \
  X(WARP_OUTER_TRANSLATION_Y, params.warp.outer.translation_y)                 \
  X(WARP_OUTER_ROTATION, params.warp.outer.rotation)                           \
  X(WARP_OUTER_SCALE_X, params.warp.outer.scale_x)                             \
  X(WARP_OUTER_SCALE_Y, params.warp.outer.scale_y)                             \
  X(WARP_OUTER_SHEAR, params.warp.outer.shear)                                 \
  X(WARP_OUTER_FREQUENCY, params.warp.outer.frequency)                         \
  X(WARP_OUTER_FIELD_ANGLE, params.warp.outer.field_angle)                     \
  X(WARP_OUTER_CENTER_X, params.warp.outer.center_x)                           \
  X(WARP_OUTER_CENTER_Y, params.warp.outer.center_y)                           \
  X(WARP_OUTER_RADIUS, params.warp.outer.radius)                               \
  X(WARP_OUTER_TURNS, params.warp.outer.turns)                                 \
  X(WARP_OUTER_CENTER_ORBIT_RADIUS, params.warp.outer.center_orbit_radius)     \
  X(WARP_OUTER_VECTOR_ANGLE, params.warp.outer.vector_angle)                   \
  X(WARP_OUTER_CELL_X, params.warp.outer.cell_x)                               \
  X(WARP_OUTER_CELL_Y, params.warp.outer.cell_y)                               \
  X(WARP_OUTER_OFFSET_X, params.warp.outer.offset_x)                           \
  X(WARP_OUTER_OFFSET_Y, params.warp.outer.offset_y)                           \
  X(WARP_OUTER_RADIAL_SCALE, params.warp.outer.radial_scale)                   \
  X(WARP_OUTER_RADIAL_PHASE, params.warp.outer.radial_phase)                   \
  X(WARP_OUTER_ANGULAR_PHASE, params.warp.outer.angular_phase)                 \
  X(WARP_OUTER_EDGE_WIDTH, params.warp.outer.edge_width)                       \
  X(WARP_INNER_SCALE, params.warp.inner.scale)                                 \
  X(WARP_INNER_STRENGTH, params.warp.inner.strength)                           \
  X(WARP_INNER_SPEED, params.warp.inner.speed)                                 \
  X(WARP_INNER_TRANSLATION_X, params.warp.inner.translation_x)                 \
  X(WARP_INNER_TRANSLATION_Y, params.warp.inner.translation_y)                 \
  X(WARP_INNER_ROTATION, params.warp.inner.rotation)                           \
  X(WARP_INNER_SCALE_X, params.warp.inner.scale_x)                             \
  X(WARP_INNER_SCALE_Y, params.warp.inner.scale_y)                             \
  X(WARP_INNER_SHEAR, params.warp.inner.shear)                                 \
  X(WARP_INNER_FREQUENCY, params.warp.inner.frequency)                         \
  X(WARP_INNER_FIELD_ANGLE, params.warp.inner.field_angle)                     \
  X(WARP_INNER_CENTER_X, params.warp.inner.center_x)                           \
  X(WARP_INNER_CENTER_Y, params.warp.inner.center_y)                           \
  X(WARP_INNER_RADIUS, params.warp.inner.radius)                               \
  X(WARP_INNER_TURNS, params.warp.inner.turns)                                 \
  X(WARP_INNER_CENTER_ORBIT_RADIUS, params.warp.inner.center_orbit_radius)     \
  X(WARP_INNER_VECTOR_ANGLE, params.warp.inner.vector_angle)                   \
  X(WARP_INNER_CELL_X, params.warp.inner.cell_x)                               \
  X(WARP_INNER_CELL_Y, params.warp.inner.cell_y)                               \
  X(WARP_INNER_OFFSET_X, params.warp.inner.offset_x)                           \
  X(WARP_INNER_OFFSET_Y, params.warp.inner.offset_y)                           \
  X(WARP_INNER_RADIAL_SCALE, params.warp.inner.radial_scale)                   \
  X(WARP_INNER_RADIAL_PHASE, params.warp.inner.radial_phase)                   \
  X(WARP_INNER_ANGULAR_PHASE, params.warp.inner.angular_phase)                 \
  X(WARP_INNER_EDGE_WIDTH, params.warp.inner.edge_width)                       \
  X(PROJECTION_SINGULARITY_FADE, params.projection.singularity_fade)           \
  X(PROJECTION_SPIN_RATE, params.projection.spin_rate)                         \
  X(PROJECTION_WANDER, params.projection.wander)                               \
  X(PROJECTION_CENTRAL_MERIDIAN, params.projection.central_meridian)           \
  X(PROJECTION_COORDINATE_SCALE, params.projection.coordinate_scale)           \
  X(PROJECTION_BONNE_STANDARD_PARALLEL,                                        \
    params.projection.bonne_standard_parallel)                                 \
  X(PROJECTION_LAYOUT_SCROLL, params.projection.layout_scroll)                 \
  X(LENS_MOBIUS_A_RE, params.surface_lens.mobius.a.re)                         \
  X(LENS_MOBIUS_A_IM, params.surface_lens.mobius.a.im)                         \
  X(LENS_MOBIUS_B_RE, params.surface_lens.mobius.b.re)                         \
  X(LENS_MOBIUS_B_IM, params.surface_lens.mobius.b.im)                         \
  X(LENS_MOBIUS_C_RE, params.surface_lens.mobius.c.re)                         \
  X(LENS_MOBIUS_C_IM, params.surface_lens.mobius.c.im)                         \
  X(LENS_MOBIUS_D_RE, params.surface_lens.mobius.d.re)                         \
  X(LENS_MOBIUS_D_IM, params.surface_lens.mobius.d.im)                         \
  X(VALUE_ISO_LEVEL, params.value.iso_level)                                   \
  X(VALUE_ISO_WIDTH, params.value.iso_width)                                   \
  X(VALUE_BAND_COUNT, params.value.band_count)                                 \
  X(VALUE_BAND_PHASE, params.value.band_phase)                                 \
  X(VALUE_CUTOUT_THRESHOLD, params.value.cutout_threshold)                     \
  X(VALUE_CUTOUT_SOFTNESS, params.value.cutout_softness)                       \
  X(VALUE_EDGE_WIDTH, params.value.edge_width)                                 \
  X(COLOR_HUE_SHIFT_AMOUNT, params.color.hue_shift_amount)                     \
  X(COLOR_HUE_NOISE_SCALE, params.color.hue_noise_scale)                       \
  X(COLOR_HUE_NOISE_SPEED, params.color.hue_noise_speed)                       \
  X(COLOR_PALETTE_CHROMA, params.color.palette_chroma)                         \
  X(COLOR_MAPPING_FREQUENCY, params.color.mapping_frequency)                   \
  X(COLOR_MAPPING_PHASE, params.color.mapping_phase)                           \
  X(COLOR_PHASE_OSCILLATION_DEPTH, params.color.phase_oscillation_depth)       \
  X(COLOR_PHASE_OSCILLATION_SPEED, params.color.phase_oscillation_speed)       \
  X(COLOR_BRIGHTNESS_DEPTH, params.color.brightness_depth)                     \
  X(COLOR_VALUE_OPACITY_LOW, params.color.opacity_low)                         \
  X(COLOR_VALUE_OPACITY_HIGH, params.color.opacity_high)                       \
  X(CAMERA_WANDER, params.outer_camera.wander)                                 \
  X(SURFACE_NOISE_BASIS, params.surface_noise.basis)                           \
  X(SURFACE_NOISE_INTEGRATOR, params.surface_noise.integrator)                 \
  X(SURFACE_NOISE_SEED, params.surface_noise.seed)                             \
  X(SURFACE_NOISE_SCALE, params.surface_noise.scale)                           \
  X(SURFACE_NOISE_STRENGTH, params.surface_noise.strength)                     \
  X(SURFACE_NOISE_RATE, params.surface_noise.rate)                             \
  X(SURFACE_NOISE_DIRECTION, params.surface_noise.direction)                   \
  X(SOURCE_RING_COUNT, params.source.ring_count)                               \
  X(SOURCE_RING_THICKNESS, params.source.ring_thickness)                       \
  X(SOURCE_RING_SOFTNESS, params.source.ring_softness)                         \
  X(SOURCE_RING_WANDER, params.source.ring_wander)                             \
  X(SOURCE_FRACTAL_SCALE, params.source.fractal_scale)                         \
  X(SOURCE_FRACTAL_ITERATIONS, params.source.fractal_iterations)               \
  X(SOURCE_JULIA_MIX, params.source.julia_mix)                                 \
  X(SOURCE_JULIA_REAL, params.source.julia_real)                               \
  X(SOURCE_JULIA_IMAGINARY, params.source.julia_imaginary)                     \
  X(SOURCE_FRACTAL_CONTOURS, params.source.fractal_contours)                   \
  X(SOURCE_TESSELLATION_CELL_SCALE, params.source.tessellation_cell_scale)     \
  X(SOURCE_TESSELLATION_LINE_THICKNESS,                                        \
    params.source.tessellation_line_thickness)                                 \
  X(SOURCE_TESSELLATION_LINE_SOFTNESS,                                         \
    params.source.tessellation_line_softness)                                  \
  X(SOURCE_TESSELLATION_KIND, params.source.tessellation_kind)                 \
  X(COLOR_PALETTE_MAPPING, params.color.palette_mapping)

/**
 * @brief Slot-based sphere shader with an immutable per-frame pullback state.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> class ShaderWorkbench : public Effect {
public:
  // The authored vocabulary, re-exported so consumers keep naming it
  // through the effect.
  using Function = Workbench::Function;
  using Projection = Workbench::Projection;
  using PeirceLayout = Workbench::PeirceLayout;
  using AiroceanLayout = Workbench::AiroceanLayout;
  using BonneHemisphere = Workbench::BonneHemisphere;
  using GnomonicHemispherePolicy = Workbench::GnomonicHemispherePolicy;
  using SurfaceLens = Workbench::SurfaceLens;
  using WarpEnvelope = Workbench::WarpEnvelope;
  using PolarMode = Workbench::PolarMode;
  using CurlIntegrator = Workbench::CurlIntegrator;
  using SurfaceCurlIntegrator = Workbench::SurfaceCurlIntegrator;
  using SurfaceNoise = Workbench::SurfaceNoise;
  using SurfaceNoisePlacement = Workbench::SurfaceNoisePlacement;
  using WarpStageKind = Workbench::WarpStageKind;
  using WarpStageSpec = Workbench::WarpStageSpec;
  using WarpProgram = Workbench::WarpProgram;
  using ProjectionFramePolicy = Workbench::ProjectionFramePolicy;
  using SignalWeight = Workbench::SignalWeight;
  using ValueTransfer = Workbench::ValueTransfer;
  using CoveragePolicy = Workbench::CoveragePolicy;
  using PaletteMode = Workbench::PaletteMode;
  using PaletteMapping = Workbench::PaletteMapping;
  using PaletteMappingWeights = Workbench::PaletteMappingWeights;
  using BrightnessEnvelope = Workbench::BrightnessEnvelope;
  using HueShiftMode = Workbench::HueShiftMode;
  using Slots = Workbench::Slots;
  using SourceParams = Workbench::SourceParams;
  using WarpStageParams = Workbench::WarpStageParams;
  using WarpParams = Workbench::WarpParams;
  using ProjectionParams = Workbench::ProjectionParams;
  using SurfaceLensParams = Workbench::SurfaceLensParams;
  using SurfaceNoiseParams = Workbench::SurfaceNoiseParams;
  using ValueParams = Workbench::ValueParams;
  using ColorParams = Workbench::ColorParams;
  using OuterCameraParams = Workbench::OuterCameraParams;
  using Params = Workbench::Params;
  using Config = Workbench::Config;
  using RequestedConfig = Workbench::RequestedConfig;
  using Blend = Workbench::Blend;
  using Choreo = Workbench::Choreo;
  using SourceState = Workbench::SourceState;
  using ProjectedLookup = Workbench::ProjectedLookup;
  using SourceTraits = Workbench::SourceTraits;
  using PlanarWarpResult = Workbench::PlanarWarpResult;
  using SurfaceNoiseResult = Workbench::SurfaceNoiseResult;
  using PlanarWarpStageResult = Workbench::PlanarWarpStageResult;
  using FieldSample = Workbench::FieldSample;
  using ClockState = Workbench::ClockState;
  using PreparedTransforms = Workbench::PreparedTransforms;
  using PreparedAffineFrame = Workbench::PreparedAffineFrame;
  using PreparedMirrorTile = Workbench::PreparedMirrorTile;
  using PreparedVortex = Workbench::PreparedVortex;
  using PreparedNoiseLoop = Workbench::PreparedNoiseLoop;
  using PreparedWarpStage = Workbench::PreparedWarpStage;
  using PreparedWarpProgram = Workbench::PreparedWarpProgram;
  using PreparedSurfaceNoise = Workbench::PreparedSurfaceNoise;
  using PreparedHueRotation = Workbench::PreparedHueRotation;
  using PreparedHueNoise = Workbench::PreparedHueNoise;
  using ResourceBindings = Workbench::ResourceBindings;
  static constexpr auto MAX_NOISE_RESOURCES = Workbench::MAX_NOISE_RESOURCES;
  using DynamicPrepared = Workbench::DynamicPrepared;
  using FrameState = Workbench::FrameState;
  using InversePipelineId = Workbench::InversePipelineId;
  using SelectedConfig = Workbench::SelectedConfig;
  using Preset = Workbench::Preset;
  static constexpr auto &PRESETS = Workbench::PRESETS;
  static constexpr auto BOUNDARY_CUT = Workbench::BOUNDARY_CUT;
  static constexpr auto BOUNDARY_SINGULAR = Workbench::BOUNDARY_SINGULAR;
  static constexpr auto PROJECTION_FLAG_FOLDED =
      Workbench::PROJECTION_FLAG_FOLDED;
  static constexpr auto GNOMONIC_AXIS_EPS = Workbench::GNOMONIC_AXIS_EPS;
  static constexpr auto WARP_COORD_LIMIT = Workbench::WARP_COORD_LIMIT;
  static constexpr auto NOISE_LATTICE_LIMIT = Workbench::NOISE_LATTICE_LIMIT;
  static constexpr auto &FUNCTION_OPTIONS = Workbench::FUNCTION_OPTIONS;
  static constexpr auto &FUNCTION_EXPORT_OPTIONS =
      Workbench::FUNCTION_EXPORT_OPTIONS;
  static constexpr auto NUM_FUNCTIONS = Workbench::NUM_FUNCTIONS;
  static constexpr auto &TESSELLATION_KIND_OPTIONS =
      Workbench::TESSELLATION_KIND_OPTIONS;
  static constexpr auto &TESSELLATION_KIND_EXPORT_OPTIONS =
      Workbench::TESSELLATION_KIND_EXPORT_OPTIONS;
  static constexpr auto NUM_TESSELLATION_KINDS =
      Workbench::NUM_TESSELLATION_KINDS;
  static constexpr auto &PROJECTION_OPTIONS = Workbench::PROJECTION_OPTIONS;
  static constexpr auto &PROJECTION_EXPORT_OPTIONS =
      Workbench::PROJECTION_EXPORT_OPTIONS;
  static constexpr auto NUM_PROJECTIONS = Workbench::NUM_PROJECTIONS;
  static constexpr auto &PEIRCE_LAYOUT_OPTIONS =
      Workbench::PEIRCE_LAYOUT_OPTIONS;
  static constexpr auto &PEIRCE_LAYOUT_EXPORT_OPTIONS =
      Workbench::PEIRCE_LAYOUT_EXPORT_OPTIONS;
  static constexpr auto NUM_PEIRCE_LAYOUTS = Workbench::NUM_PEIRCE_LAYOUTS;
  static constexpr auto &AIROCEAN_LAYOUT_OPTIONS =
      Workbench::AIROCEAN_LAYOUT_OPTIONS;
  static constexpr auto &AIROCEAN_LAYOUT_EXPORT_OPTIONS =
      Workbench::AIROCEAN_LAYOUT_EXPORT_OPTIONS;
  static constexpr auto NUM_AIROCEAN_LAYOUTS = Workbench::NUM_AIROCEAN_LAYOUTS;
  static constexpr auto &BONNE_HEMISPHERE_OPTIONS =
      Workbench::BONNE_HEMISPHERE_OPTIONS;
  static constexpr auto &BONNE_HEMISPHERE_EXPORT_OPTIONS =
      Workbench::BONNE_HEMISPHERE_EXPORT_OPTIONS;
  static constexpr auto NUM_BONNE_HEMISPHERES =
      Workbench::NUM_BONNE_HEMISPHERES;
  static constexpr auto &GNOMONIC_HEMISPHERE_OPTIONS =
      Workbench::GNOMONIC_HEMISPHERE_OPTIONS;
  static constexpr auto &GNOMONIC_HEMISPHERE_EXPORT_OPTIONS =
      Workbench::GNOMONIC_HEMISPHERE_EXPORT_OPTIONS;
  static constexpr auto NUM_GNOMONIC_HEMISPHERES =
      Workbench::NUM_GNOMONIC_HEMISPHERES;
  static constexpr auto &PROJECTION_FRAME_OPTIONS =
      Workbench::PROJECTION_FRAME_OPTIONS;
  static constexpr auto &PROJECTION_FRAME_EXPORT_OPTIONS =
      Workbench::PROJECTION_FRAME_EXPORT_OPTIONS;
  static constexpr auto NUM_PROJECTION_FRAMES =
      Workbench::NUM_PROJECTION_FRAMES;
  static constexpr auto &LENS_OPTIONS = Workbench::LENS_OPTIONS;
  static constexpr auto &LENS_EXPORT_OPTIONS = Workbench::LENS_EXPORT_OPTIONS;
  static constexpr auto NUM_LENSES = Workbench::NUM_LENSES;
  static constexpr auto &SURFACE_NOISE_OPTIONS =
      Workbench::SURFACE_NOISE_OPTIONS;
  static constexpr auto &SURFACE_NOISE_EXPORT_OPTIONS =
      Workbench::SURFACE_NOISE_EXPORT_OPTIONS;
  static constexpr auto NUM_SURFACE_NOISE = Workbench::NUM_SURFACE_NOISE;
  static constexpr auto &SURFACE_NOISE_PLACEMENT_OPTIONS =
      Workbench::SURFACE_NOISE_PLACEMENT_OPTIONS;
  static constexpr auto &SURFACE_NOISE_PLACEMENT_EXPORT_OPTIONS =
      Workbench::SURFACE_NOISE_PLACEMENT_EXPORT_OPTIONS;
  static constexpr auto NUM_SURFACE_NOISE_PLACEMENTS =
      Workbench::NUM_SURFACE_NOISE_PLACEMENTS;
  static constexpr auto &SURFACE_CURL_INTEGRATOR_OPTIONS =
      Workbench::SURFACE_CURL_INTEGRATOR_OPTIONS;
  static constexpr auto &SURFACE_CURL_INTEGRATOR_EXPORT_OPTIONS =
      Workbench::SURFACE_CURL_INTEGRATOR_EXPORT_OPTIONS;
  static constexpr auto NUM_SURFACE_CURL_INTEGRATORS =
      Workbench::NUM_SURFACE_CURL_INTEGRATORS;
  static constexpr auto &WARP_OPTIONS = Workbench::WARP_OPTIONS;
  static constexpr auto &WARP_EXPORT_OPTIONS = Workbench::WARP_EXPORT_OPTIONS;
  static constexpr auto NUM_WARPS = Workbench::NUM_WARPS;
  static constexpr auto &NOISE_BASIS_OPTIONS = Workbench::NOISE_BASIS_OPTIONS;
  static constexpr auto &NOISE_BASIS_EXPORT_OPTIONS =
      Workbench::NOISE_BASIS_EXPORT_OPTIONS;
  static constexpr auto NUM_NOISE_BASES = Workbench::NUM_NOISE_BASES;
  static constexpr auto &POLAR_MODE_OPTIONS = Workbench::POLAR_MODE_OPTIONS;
  static constexpr auto &POLAR_MODE_EXPORT_OPTIONS =
      Workbench::POLAR_MODE_EXPORT_OPTIONS;
  static constexpr auto NUM_POLAR_MODES = Workbench::NUM_POLAR_MODES;
  static constexpr auto &CURL_INTEGRATOR_OPTIONS =
      Workbench::CURL_INTEGRATOR_OPTIONS;
  static constexpr auto &CURL_INTEGRATOR_EXPORT_OPTIONS =
      Workbench::CURL_INTEGRATOR_EXPORT_OPTIONS;
  static constexpr auto NUM_CURL_INTEGRATORS = Workbench::NUM_CURL_INTEGRATORS;
  static constexpr auto POLAR_HARMONIC_MAX = Workbench::POLAR_HARMONIC_MAX;
  static constexpr auto BAND_COUNT_MAX = Workbench::BAND_COUNT_MAX;
  static constexpr auto &WARP_ENVELOPE_OPTIONS =
      Workbench::WARP_ENVELOPE_OPTIONS;
  static constexpr auto &WARP_ENVELOPE_EXPORT_OPTIONS =
      Workbench::WARP_ENVELOPE_EXPORT_OPTIONS;
  static constexpr auto NUM_WARP_ENVELOPES = Workbench::NUM_WARP_ENVELOPES;
  static constexpr auto &SIGNAL_OPTIONS = Workbench::SIGNAL_OPTIONS;
  static constexpr auto &SIGNAL_EXPORT_OPTIONS =
      Workbench::SIGNAL_EXPORT_OPTIONS;
  static constexpr auto NUM_SIGNALS = Workbench::NUM_SIGNALS;
  static constexpr auto &VALUE_TRANSFER_OPTIONS =
      Workbench::VALUE_TRANSFER_OPTIONS;
  static constexpr auto &VALUE_TRANSFER_EXPORT_OPTIONS =
      Workbench::VALUE_TRANSFER_EXPORT_OPTIONS;
  static constexpr auto NUM_VALUE_TRANSFERS = Workbench::NUM_VALUE_TRANSFERS;
  static constexpr auto &COVERAGE_OPTIONS = Workbench::COVERAGE_OPTIONS;
  static constexpr auto &COVERAGE_EXPORT_OPTIONS =
      Workbench::COVERAGE_EXPORT_OPTIONS;
  static constexpr auto NUM_COVERAGE_POLICIES =
      Workbench::NUM_COVERAGE_POLICIES;
  static constexpr auto &PALETTE_OPTIONS = Workbench::PALETTE_OPTIONS;
  static constexpr auto &PALETTE_EXPORT_OPTIONS =
      Workbench::PALETTE_EXPORT_OPTIONS;
  static constexpr auto NUM_PALETTES = Workbench::NUM_PALETTES;
  static constexpr auto &PALETTE_MAPPING_OPTIONS =
      Workbench::PALETTE_MAPPING_OPTIONS;
  static constexpr auto &PALETTE_MAPPING_EXPORT_OPTIONS =
      Workbench::PALETTE_MAPPING_EXPORT_OPTIONS;
  static constexpr auto NUM_PALETTE_MAPPINGS = Workbench::NUM_PALETTE_MAPPINGS;
  static constexpr auto &BRIGHTNESS_ENVELOPE_OPTIONS =
      Workbench::BRIGHTNESS_ENVELOPE_OPTIONS;
  static constexpr auto &BRIGHTNESS_ENVELOPE_EXPORT_OPTIONS =
      Workbench::BRIGHTNESS_ENVELOPE_EXPORT_OPTIONS;
  static constexpr auto NUM_BRIGHTNESS_ENVELOPES =
      Workbench::NUM_BRIGHTNESS_ENVELOPES;
  static constexpr auto &HUE_SHIFT_OPTIONS = Workbench::HUE_SHIFT_OPTIONS;
  static constexpr auto &HUE_SHIFT_EXPORT_OPTIONS =
      Workbench::HUE_SHIFT_EXPORT_OPTIONS;
  static constexpr auto NUM_HUE_SHIFT_MODES = Workbench::NUM_HUE_SHIFT_MODES;
  using WarpParamName = Workbench::WarpParamName;

  static constexpr auto WARP_SCALE_MIN = Workbench::WARP_SCALE_MIN;
  static constexpr auto WARP_SCALE_MAX = Workbench::WARP_SCALE_MAX;
  static constexpr auto WARP_STRENGTH_MIN = Workbench::WARP_STRENGTH_MIN;
  static constexpr auto WARP_STRENGTH_MAX = Workbench::WARP_STRENGTH_MAX;
  static constexpr auto VECTOR_WARP_SCALE_MAX =
      Workbench::VECTOR_WARP_SCALE_MAX;
  static constexpr auto VECTOR_WARP_STRENGTH_MAX =
      Workbench::VECTOR_WARP_STRENGTH_MAX;
  static constexpr auto CURL_WARP_SCALE_MAX = Workbench::CURL_WARP_SCALE_MAX;
  static constexpr auto CURL_WARP_STRENGTH_MAX =
      Workbench::CURL_WARP_STRENGTH_MAX;
  static constexpr auto CURL_VECTOR_COMPONENT_MAX =
      Workbench::CURL_VECTOR_COMPONENT_MAX;
  static constexpr auto WARP_SPEED_MIN = Workbench::WARP_SPEED_MIN;
  static constexpr auto WARP_SPEED_MAX = Workbench::WARP_SPEED_MAX;
  static constexpr auto PATTERN_FREQ_MIN = Workbench::PATTERN_FREQ_MIN;
  static constexpr auto PATTERN_FREQ_MAX = Workbench::PATTERN_FREQ_MAX;
  static constexpr auto GRID_PATTERN_FREQ_MAX =
      Workbench::GRID_PATTERN_FREQ_MAX;
  static constexpr auto SPEED_MIN = Workbench::SPEED_MIN;
  static constexpr auto COMPLEXITY_MIN = Workbench::COMPLEXITY_MIN;
  static constexpr auto PATTERN_MIX_MIN = Workbench::PATTERN_MIX_MIN;
  static constexpr auto PHASE2_RATE_MIN = Workbench::PHASE2_RATE_MIN;
  static constexpr auto SINGULARITY_FADE_MAX = Workbench::SINGULARITY_FADE_MAX;
  static constexpr auto SINGULARITY_FADE_MIN = Workbench::SINGULARITY_FADE_MIN;
  static constexpr auto SPIN_RATE_MIN = Workbench::SPIN_RATE_MIN;
  static constexpr auto WANDER_MIN = Workbench::WANDER_MIN;
  static constexpr auto HUE_SHIFT_AMOUNT_MAX = Workbench::HUE_SHIFT_AMOUNT_MAX;
  static constexpr auto HUE_NOISE_AMOUNT_MAX = Workbench::HUE_NOISE_AMOUNT_MAX;
  static constexpr auto HUE_NOISE_SCALE_MIN = Workbench::HUE_NOISE_SCALE_MIN;
  static constexpr auto HUE_NOISE_SCALE_MAX = Workbench::HUE_NOISE_SCALE_MAX;
  static constexpr auto HUE_NOISE_SPEED_MAX = Workbench::HUE_NOISE_SPEED_MAX;
  static constexpr auto PALETTE_CHROMA_MIN = Workbench::PALETTE_CHROMA_MIN;
  static constexpr auto PALETTE_CHROMA_MAX = Workbench::PALETTE_CHROMA_MAX;
  static constexpr auto MAPPING_FREQUENCY_MIN =
      Workbench::MAPPING_FREQUENCY_MIN;
  static constexpr auto MAPPING_FREQUENCY_MAX =
      Workbench::MAPPING_FREQUENCY_MAX;
  static constexpr auto MAPPING_PHASE_MIN = Workbench::MAPPING_PHASE_MIN;
  static constexpr auto MAPPING_PHASE_MAX = Workbench::MAPPING_PHASE_MAX;
  static constexpr auto PHASE_OSCILLATION_DEPTH_MIN =
      Workbench::PHASE_OSCILLATION_DEPTH_MIN;
  static constexpr auto PHASE_OSCILLATION_DEPTH_MAX =
      Workbench::PHASE_OSCILLATION_DEPTH_MAX;
  static constexpr auto PHASE_OSCILLATION_SPEED_MAX =
      Workbench::PHASE_OSCILLATION_SPEED_MAX;
  static constexpr auto BRIGHTNESS_DEPTH_MIN = Workbench::BRIGHTNESS_DEPTH_MIN;
  static constexpr auto BRIGHTNESS_DEPTH_MAX = Workbench::BRIGHTNESS_DEPTH_MAX;
  static constexpr auto VALUE_OPACITY_MIN = Workbench::VALUE_OPACITY_MIN;
  static constexpr auto VALUE_OPACITY_MAX = Workbench::VALUE_OPACITY_MAX;
  static constexpr auto WAVE_SPIN_MIN = Workbench::WAVE_SPIN_MIN;
  static constexpr auto SOURCE_NOISE_SCALE_MIN =
      Workbench::SOURCE_NOISE_SCALE_MIN;
  static constexpr auto SOURCE_NOISE_SCALE_MAX =
      Workbench::SOURCE_NOISE_SCALE_MAX;
  static constexpr auto SOURCE_NOISE_RATE_MIN =
      Workbench::SOURCE_NOISE_RATE_MIN;
  static constexpr auto SOURCE_NOISE_RATE_MAX =
      Workbench::SOURCE_NOISE_RATE_MAX;
  static constexpr auto LENS_NOISE_SCALE_MIN = Workbench::LENS_NOISE_SCALE_MIN;
  static constexpr auto LENS_NOISE_SCALE_MAX = Workbench::LENS_NOISE_SCALE_MAX;
  static constexpr auto NOISE_RATE_MIN = Workbench::NOISE_RATE_MIN;
  static constexpr auto NOISE_RATE_MAX = Workbench::NOISE_RATE_MAX;
  static constexpr auto NOISE_SPEED_MIN = Workbench::NOISE_SPEED_MIN;
  static constexpr auto NOISE_SPEED_MAX = Workbench::NOISE_SPEED_MAX;
  static constexpr auto CELL_MIN = Workbench::CELL_MIN;
  static constexpr auto CELL_MAX = Workbench::CELL_MAX;
  static constexpr auto SOFTNESS_MIN = Workbench::SOFTNESS_MIN;
  using WarpParamName = Workbench::WarpParamName;
  static constexpr auto SPEED_MAX = Workbench::SPEED_MAX;
  static constexpr auto COMPLEXITY_MAX = Workbench::COMPLEXITY_MAX;
  static constexpr auto PATTERN_MIX_MAX = Workbench::PATTERN_MIX_MAX;
  static constexpr auto PHASE2_RATE_MAX = Workbench::PHASE2_RATE_MAX;
  static constexpr auto SPIN_RATE_MAX = Workbench::SPIN_RATE_MAX;
  static constexpr auto WANDER_MAX = Workbench::WANDER_MAX;
  static constexpr auto WAVE_SPIN_MAX = Workbench::WAVE_SPIN_MAX;

private:
  struct EndpointRuntime;
  struct WalkDeltas;
  using ShadeFunction = Color4 (*)(const Vector &, const FrameState &,
                                   const void *);
  struct FrameShader {
    using ShadeFunction = ShaderWorkbench::ShadeFunction;

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

public:
  static constexpr std::string_view EFFECT_ID = "shader";
  static constexpr int GAMUT_ANGLE_STEPS = GAMUT_LUT_ANGLE_STEPS;
  static constexpr int GAMUT_L_STEPS = GAMUT_LUT_L_STEPS;

  static constexpr size_t authored_preset_count() { return PRESETS.size(); }

  HS_COLD_MEMBER ShaderWorkbench() : Effect(W, H, {.strobe = true}) {}

protected:
  HS_COLD_MEMBER void
  set_fixed_preset_view(std::span<const uint8_t> source_indices) {
    HS_CHECK(!source_indices.empty(),
             "set_fixed_preset_view: empty preset view");
    for (uint8_t index : source_indices)
      HS_CHECK(index < PRESETS.size(),
               "set_fixed_preset_view: preset index out of range");
    preset_view = source_indices;
    fixed_topology = true;
  }

  HS_COLD_MEMBER void hold_initial_preset(uint16_t frames) {
    preset_dwell_remaining = frames;
    preset_dwell_armed = preset_count_for_view() > 1;
  }

public:
  /** @brief Initializes slots, clocks, palette resources, and choreography. */
  HS_COLD_MEMBER void init() override {
    configure_presets(preset_count_for_view());
#if HS_ENABLE_PARAM_GUI_BRIDGE
    set_parameter_updated_hook(&ShaderWorkbench::dispatch_parameter_updated);
#endif
    state = persistent_arena.make<StateBundle>();
    use_parameter_storage(persistent_arena.allocate_n<ParamDef>(PARAM_CAPACITY),
                          PARAM_CAPACITY);
    const Preset &initial = preset_for_view(0);
    active_slots = initial.config.slots;
    active_pipeline = initial.pipeline;
    blend.params = initial.config.params;
    blend.palette_mapping =
        palette_mapping_weights(initial.config.slots.palette_mapping);
#if HS_ENABLE_PARAM_GUI_BRIDGE
    display_config = initial.config;
#endif
    requested_config = initial.config;
    published_config = initial.config;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    accepted_config = initial.config;
#endif
    prepare_resource_union(initial.config, initial.config);

    rebind_parameters();

    timeline.add(0, Animation::RandomWalk<W>(projection_walk, UP,
                                             state->projection_walk_noise));
    timeline.add(
        0, Animation::RandomWalk<W>(outer_walk, UP, state->outer_walk_noise));

    init_gamut_lut(persistent_arena, GAMUT_ANGLE_STEPS, GAMUT_L_STEPS);
    generated_palettes.init(persistent_arena, 0.62f, ease_in_out_sin);
    update_palette_chroma(
        preset_for_view(0).config.params.color.palette_chroma);

    enter_preset();
  }

  /** @brief Advances mutable state, snapshots it, and renders one frame. */
  HS_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(sb_timeline_step);
      timeline.step(canvas);
    }
    advance_preset_choreography();

    apply_requested_config();
    prepare_param_morph();
    state->render_config.slots = active_slots;
    state->render_config.params = blend.params;
    const WalkDeltas walk_deltas = sample_walk_deltas();
    if (state->transition.active) {
      if (state->transition.elapsed < state->transition.duration / 2)
        advance_runtime(state->transition.from_runtime,
                        state->transition.from_config, walk_deltas);
      advance_runtime(state->transition.to_runtime, state->transition.to_config,
                      walk_deltas);
    } else {
      advance_runtime(runtime, state->render_config, walk_deltas);
    }
    update_palette_chroma(visible_palette_chroma());
    step_generated_palettes(visible_palette_mode());
#if HS_ENABLE_TEST_HOOKS
    ++generated_palette_step_count;
#endif

    if (state->transition.active) {
      draw_through_clear_transition(canvas);
    } else {
      PreparedEndpoint prepared;
      HS_CHECK(prepare_endpoint(state->render_config, runtime, 1.0f,
                                active_pipeline, prepared),
               "ShaderWorkbench active endpoint has no renderer");
      draw_endpoint(canvas, prepared);
    }
    finish_transitions();
    publish_live_config();
  }

#if HS_ENABLE_EFFECT_CONTROL_API
  void profile_select_preset(size_t index) {
    HS_CHECK(index < preset_count_for_view(),
             "ShaderWorkbench profile preset index out of range");
    HS_CHECK(selectPreset(index),
             "ShaderWorkbench profile preset selection failed");
    hs::log("Profile preset: %u/%u", static_cast<unsigned>(index),
            static_cast<unsigned>(preset_count_for_view()));
  }
#endif

private:
  friend struct ::hs_test::shader_workbench_tests::ShaderWorkbenchWhiteBox;
  using NoiseBasis = ::NoiseBasis;

  HS_COLD_MEMBER bool apply_preset(const PresetChange &change) override {
    const size_t index = change.to;
    if (change.origin == PresetChangeOrigin::AUTOMATIC) {
      const Choreo choreo = preset_choreo();
      const Preset &to = preset_for_view(index);
      if (!try_apply_config(to.config, choreo.blend_frames, choreo.staggered,
                            true))
        return false;
      requested_config = to.config;
      published_config = to.config;
#if HS_ENABLE_PARAM_GUI_BRIDGE
      accepted_config = to.config;
      pending_edit_count = 0;
#endif
      rebind_parameters();
      return true;
    }

    state->param_morph.active = false;
    state->transition.active = false;
    const SelectedConfig &selected = preset_for_view(index);
    active_slots = selected.config.slots;
    active_pipeline = selected.pipeline;
    blend.params = selected.config.params;
    blend.palette_mapping =
        palette_mapping_weights(selected.config.slots.palette_mapping);
#if HS_ENABLE_PARAM_GUI_BRIDGE
    display_config = selected.config;
#endif
    requested_config = selected.config;
    published_config = selected.config;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    accepted_config = selected.config;
    pending_edit_count = 0;
#endif
    runtime = {};
    HS_CHECK(prepare_resource_union(selected.config, selected.config),
             "ShaderWorkbench preset resources exceed capacity");
    rebind_parameters();
    return true;
  }

  HS_COLD_MEMBER void preset_changed(const PresetChange &) override {
    if (!state->param_morph.active && !state->transition.active)
      enter_preset();
  }

public:
  static constexpr uint32_t CONFIG_SCHEMA_VERSION = 10;

  /**
   * @brief Reports whether a persisted snapshot's schema version can be
   *        restored.
   * @param version Snapshot schema version to test.
   * @return true for the current version.
   * @details Single source of truth for the accepted set, so callers that
   *          pre-screen a version (the WASM bridge) cannot drift from what the
   *          effect actually accepts.
   */
  static constexpr bool config_version_supported(uint32_t version) {
    return version == CONFIG_SCHEMA_VERSION;
  }

  enum class ConfigFieldId : uint16_t {
#define HS_SHADER_WORKBENCH_FIELD_ENUM(name, path) name,
    HS_SHADER_WORKBENCH_CONFIG_FIELDS(HS_SHADER_WORKBENCH_FIELD_ENUM)
#undef HS_SHADER_WORKBENCH_FIELD_ENUM
        COUNT
  };

  static constexpr size_t CONFIG_FIELD_COUNT =
      static_cast<size_t>(ConfigFieldId::COUNT);

  static constexpr size_t CONFIG_FIELD_BYTES =
#define HS_SHADER_WORKBENCH_FIELD_BYTES(name, path)                            \
  sizeof(std::declval<const Config &>().path) +
      HS_SHADER_WORKBENCH_CONFIG_FIELDS(HS_SHADER_WORKBENCH_FIELD_BYTES)
#undef HS_SHADER_WORKBENCH_FIELD_BYTES
          size_t{0};

  // Config's size and the listed fields' total size pin the snapshot field set
  // from both ends: an unlisted new member trips the first, a dropped list
  // entry the second. Their difference is alignment padding, so the list
  // covers every Config byte that carries a value.
  static_assert(
      sizeof(Config) == 524 && CONFIG_FIELD_BYTES == 497,
      "Config field set changed - update HS_SHADER_WORKBENCH_CONFIG_FIELDS");

  struct ConfigFieldLayout {
    size_t offset;
    size_t size;
  };

  enum class ConfigRestoreResult : uint8_t {
    APPLIED,
    UNSUPPORTED_VERSION,
    INVALID_VALUE,
    INVALID_ACCEPTED,
    INVALID_PENDING
  };

  enum class RuntimeFieldId : uint8_t {
    SOURCE_PRIMARY,
    SOURCE_SECONDARY,
    SOURCE_ANGLE,
    WARP_OUTER_ROTATION,
    PROJECTION_SPIN,
    HUE_NOISE_PHASE,
    SOURCE_NOISE_PHASE,
    WARP_INNER_ROTATION,
    SURFACE_NOISE_PHASE,
    WARP_OUTER_PHASE,
    WARP_INNER_PHASE,
    PALETTE_OSCILLATION_PHASE,
    COUNT
  };

  static constexpr size_t RUNTIME_FIELD_COUNT =
      static_cast<size_t>(RuntimeFieldId::COUNT);
  using ConfigValues = std::array<uint32_t, CONFIG_FIELD_COUNT>;
  using RuntimeValues = std::array<float, RUNTIME_FIELD_COUNT>;

  struct FullConfigSnapshot {
    uint32_t schema_version = CONFIG_SCHEMA_VERSION;
    ConfigValues accepted{};
    ConfigValues requested{};
    std::array<uint8_t, CONFIG_FIELD_COUNT> pending{};
    bool has_runtime = false;
    RuntimeValues runtime{};
  };

  struct PendingEdit {
    const char *name = nullptr;
    ConfigFieldId id = ConfigFieldId::COUNT;
    size_t offset = 0;
    size_t size = 0;
  };

  /** @brief Stable name for a full-config field ID. */
  static constexpr const char *config_field_name(ConfigFieldId id) {
    switch (id) {
#define HS_SHADER_WORKBENCH_FIELD_NAME(name, path)                             \
  case ConfigFieldId::name:                                                    \
    return #path;
      HS_SHADER_WORKBENCH_CONFIG_FIELDS(HS_SHADER_WORKBENCH_FIELD_NAME)
#undef HS_SHADER_WORKBENCH_FIELD_NAME
    case ConfigFieldId::COUNT:
      break;
    }
    return nullptr;
  }

  static ConfigFieldLayout config_field_layout(ConfigFieldId id) {
    Config config{};
    const uintptr_t base = reinterpret_cast<uintptr_t>(&config);
    switch (id) {
#define HS_SHADER_WORKBENCH_FIELD_LAYOUT(name, path)                           \
  case ConfigFieldId::name:                                                    \
    return {reinterpret_cast<uintptr_t>(&config.path) - base,                  \
            sizeof(config.path)};
      HS_SHADER_WORKBENCH_CONFIG_FIELDS(HS_SHADER_WORKBENCH_FIELD_LAYOUT)
#undef HS_SHADER_WORKBENCH_FIELD_LAYOUT
    case ConfigFieldId::COUNT:
      break;
    }
    return {sizeof(Config), 0};
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  /** @brief Reserved compatibility accessor. */
  const char *config_import_notice() const { return ""; }

  /** @brief Reserved compatibility no-op. */
  void clear_config_import_notice() {}
#endif

private:
  static constexpr float domain_scaled_max(float full_domain_max,
                                           float minimum_max,
                                           float domain_scale) {
    const float scaled = full_domain_max * domain_scale;
    return scaled > minimum_max ? scaled : minimum_max;
  }

  HS_COLD_MEMBER void register_clamped_animated_param(const char *name,
                                                      float *target,
                                                      float minimum,
                                                      float maximum) {
    const float clamped = hs::clamp(*target, minimum, maximum);
    registered_range_clamped |= clamped != *target;
    *target = clamped;
    register_animated_param(name, target, minimum, maximum);
  }

  HS_COLD_MEMBER void rebind_parameters() {
    registered_range_clamped = false;
    reset_parameters();
    Slots &slots = requested_config.slots;
    if (!fixed_topology)
      register_animated_param("Function", &slots.function, FUNCTION_OPTIONS,
                              FUNCTION_EXPORT_OPTIONS, NUM_FUNCTIONS);
    const float domain_scale = lens_domain_linear_scale(slots.surface_lens);
    register_source_controls(slots.function, requested_config.params.source,
                             domain_scale);
    if (!fixed_topology)
      register_animated_param("Projection", &slots.projection,
                              PROJECTION_OPTIONS, PROJECTION_EXPORT_OPTIONS,
                              NUM_PROJECTIONS);
    register_projection_controls(slots, requested_config.params);
    if (!fixed_topology)
      register_animated_param(
          "Projection Frame", &slots.projection_frame, PROJECTION_FRAME_OPTIONS,
          PROJECTION_FRAME_EXPORT_OPTIONS, NUM_PROJECTION_FRAMES);
    register_projection_frame_controls(slots.projection_frame,
                                       requested_config.params, domain_scale);
    register_animated_param("Camera Wander",
                            &requested_config.params.outer_camera.wander,
                            WANDER_MIN, WANDER_MAX);
    if (!fixed_topology)
      register_animated_param("Surface Noise", &slots.surface_noise,
                              SURFACE_NOISE_OPTIONS,
                              SURFACE_NOISE_EXPORT_OPTIONS, NUM_SURFACE_NOISE);
    register_surface_noise_controls(
        slots, requested_config.params.surface_noise,
        slots.surface_noise_placement == SurfaceNoisePlacement::AFTER_LENS
            ? domain_scale
            : 1.0f);
    if (!fixed_topology)
      register_animated_param("Lens", &slots.surface_lens, LENS_OPTIONS,
                              LENS_EXPORT_OPTIONS, NUM_LENSES);
    register_lens_controls(slots.surface_lens,
                           requested_config.params.surface_lens);
    if (!fixed_topology) {
      register_animated_param("Planar Warp 1", &slots.warp_program.outer.kind,
                              WARP_OPTIONS, WARP_EXPORT_OPTIONS, NUM_WARPS);
      register_stage_slot_controls(true, slots.warp_program.outer);
    }
    register_active_warp_controls(true, slots.warp_program.outer,
                                  requested_config.params.warp.outer,
                                  domain_scale);
    if (!fixed_topology) {
      register_animated_param("Planar Warp 2", &slots.warp_program.inner.kind,
                              WARP_OPTIONS, WARP_EXPORT_OPTIONS, NUM_WARPS);
      register_stage_slot_controls(false, slots.warp_program.inner);
    }
    register_active_warp_controls(false, slots.warp_program.inner,
                                  requested_config.params.warp.inner,
                                  domain_scale);
    if (!fixed_topology) {
      register_animated_param("Signal Weight", &slots.signal_weight,
                              SIGNAL_OPTIONS, SIGNAL_EXPORT_OPTIONS,
                              NUM_SIGNALS);
      register_animated_param(
          "Value Transfer", &slots.value_transfer, VALUE_TRANSFER_OPTIONS,
          VALUE_TRANSFER_EXPORT_OPTIONS, NUM_VALUE_TRANSFERS);
    }
    register_value_transfer_controls(slots.value_transfer,
                                     requested_config.params.value);
    if (!fixed_topology)
      register_animated_param("Coverage", &slots.coverage, COVERAGE_OPTIONS,
                              COVERAGE_EXPORT_OPTIONS, NUM_COVERAGE_POLICIES);
    register_coverage_controls(slots.coverage, requested_config.params.value);
    if (!fixed_topology)
      register_animated_param("Palette", &slots.palette, PALETTE_OPTIONS,
                              PALETTE_EXPORT_OPTIONS, NUM_PALETTES);
    register_animated_param("Palette Chroma",
                            &requested_config.params.color.palette_chroma,
                            PALETTE_CHROMA_MIN, PALETTE_CHROMA_MAX);
    register_animated_param(
        "Palette Mapping", &slots.palette_mapping, PALETTE_MAPPING_OPTIONS,
        PALETTE_MAPPING_EXPORT_OPTIONS, NUM_PALETTE_MAPPINGS);
    register_animated_param("Mapping Frequency",
                            &requested_config.params.color.mapping_frequency,
                            MAPPING_FREQUENCY_MIN, MAPPING_FREQUENCY_MAX);
    register_animated_param("Mapping Phase",
                            &requested_config.params.color.mapping_phase,
                            MAPPING_PHASE_MIN, MAPPING_PHASE_MAX);
    register_animated_param(
        "Phase Oscillation Depth",
        &requested_config.params.color.phase_oscillation_depth,
        PHASE_OSCILLATION_DEPTH_MIN, PHASE_OSCILLATION_DEPTH_MAX);
    register_animated_param(
        "Phase Oscillation Speed",
        &requested_config.params.color.phase_oscillation_speed,
        -PHASE_OSCILLATION_SPEED_MAX, PHASE_OSCILLATION_SPEED_MAX);
    if (!fixed_topology)
      register_animated_param("Brightness Envelope", &slots.brightness_envelope,
                              BRIGHTNESS_ENVELOPE_OPTIONS,
                              BRIGHTNESS_ENVELOPE_EXPORT_OPTIONS,
                              NUM_BRIGHTNESS_ENVELOPES);
    if (slots.brightness_envelope != BrightnessEnvelope::NONE)
      register_animated_param("Brightness Depth",
                              &requested_config.params.color.brightness_depth,
                              BRIGHTNESS_DEPTH_MIN, BRIGHTNESS_DEPTH_MAX);
    register_animated_param("Value Opacity Low",
                            &requested_config.params.color.opacity_low,
                            VALUE_OPACITY_MIN, VALUE_OPACITY_MAX);
    register_animated_param("Value Opacity High",
                            &requested_config.params.color.opacity_high,
                            VALUE_OPACITY_MIN, VALUE_OPACITY_MAX);
    if (!fixed_topology)
      register_animated_param("Hue Shift Mode", &slots.hue_shift,
                              HUE_SHIFT_OPTIONS, HUE_SHIFT_EXPORT_OPTIONS,
                              NUM_HUE_SHIFT_MODES);
    register_color_controls(slots.hue_shift, requested_config.params.color,
                            domain_scale);
#if HS_ENABLE_PARAM_GUI_BRIDGE
    const bool post_registration_clamp = clamp_registered_parameter_ranges();
    if (requested_schema_bound &&
        (registered_range_clamped || post_registration_clamp))
      refresh_accepted_config();
    for (size_t index = 0; index < pending_edit_count; ++index) {
      PendingEdit &edit = pending_edits[index];
      edit.name = nullptr;
      const uintptr_t target =
          reinterpret_cast<uintptr_t>(&requested_config) + edit.offset;
      for (const ParamDef &parameter : getParameters()) {
        if (reinterpret_cast<uintptr_t>(parameter.target) == target) {
          edit.name = parameter.name;
          break;
        }
      }
    }
    mirror_parameter_display_state(requested_config, display_config);
    for (size_t index = 0; index < pending_edit_count; ++index)
      if (pending_edits[index].name != nullptr)
        show_requested_parameter_value(pending_edits[index].name);
#endif
    requested_schema_bound = true;
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  static void dispatch_parameter_updated(Effect *effect, const char *name,
                                         bool is_enum) {
    static_cast<ShaderWorkbench *>(effect)->parameter_updated(name, is_enum);
  }

  void parameter_updated(const char *name, bool is_enum) {
    const ParamDef *parameter = getParameters().find(name);
    HS_CHECK(parameter != nullptr,
             "updated ShaderWorkbench parameter disappeared");
    const uintptr_t target = reinterpret_cast<uintptr_t>(parameter->target);
    const uintptr_t requested = reinterpret_cast<uintptr_t>(&requested_config);
    const size_t size = parameter_target_size(*parameter);
    HS_CHECK(target >= requested &&
                 target + size <= requested + sizeof(requested_config),
             "ShaderWorkbench parameter target lies outside requested config");
    const size_t offset = target - requested;
    const ConfigFieldId id = config_field_id(offset, size);
    HS_CHECK(id != ConfigFieldId::COUNT,
             "ShaderWorkbench parameter lacks a stable field ID");
    const size_t before_count = pending_edit_count;
    const bool was_pending = pending_edit_at(id) < pending_edit_count;
    remember_pending_edit(name, id, offset, size);
    refresh_accepted_config();
    const bool is_pending = pending_edit_at(id) < pending_edit_count;
    if (before_count != pending_edit_count || was_pending != is_pending ||
        (is_enum && schema_selector(name)))
      rebind_parameters();
  }

  static size_t parameter_target_size(const ParamDef &parameter) {
    switch (parameter.target_type) {
    case ParamDef::TargetType::FLOAT:
      return sizeof(float);
    case ParamDef::TargetType::INT_I32:
    case ParamDef::TargetType::INT_U32:
      return sizeof(uint32_t);
    case ParamDef::TargetType::INT_I16:
    case ParamDef::TargetType::INT_U16:
      return sizeof(uint16_t);
    case ParamDef::TargetType::BOOL:
    case ParamDef::TargetType::INT_I8:
    case ParamDef::TargetType::INT_U8:
      return sizeof(uint8_t);
    }
    __builtin_unreachable();
  }

  size_t pending_edit_at(ConfigFieldId id) const {
    for (size_t index = 0; index < pending_edit_count; ++index)
      if (pending_edits[index].id == id)
        return index;
    return pending_edit_count;
  }

  static ConfigFieldId config_field_id(size_t offset, size_t size) {
    Config config{};
    const uintptr_t base = reinterpret_cast<uintptr_t>(&config);
#define HS_SHADER_WORKBENCH_FIELD_MATCH(name, path)                            \
  if (reinterpret_cast<uintptr_t>(&config.path) - base == offset &&            \
      sizeof(config.path) == size)                                             \
    return ConfigFieldId::name;
    HS_SHADER_WORKBENCH_CONFIG_FIELDS(HS_SHADER_WORKBENCH_FIELD_MATCH)
#undef HS_SHADER_WORKBENCH_FIELD_MATCH
    return ConfigFieldId::COUNT;
  }

  void remember_pending_edit(const char *name, ConfigFieldId id, size_t offset,
                             size_t size) {
    const size_t existing = pending_edit_at(id);
    if (existing < pending_edit_count) {
      pending_edits[existing].name = name;
      pending_edits[existing].size = size;
      return;
    }
    HS_CHECK(pending_edit_count < pending_edits.size(),
             "ShaderWorkbench pending edit capacity exceeded");
    pending_edits[pending_edit_count++] = {name, id, offset, size};
  }

  void copy_pending_value(Config &to, const Config &from,
                          const PendingEdit &edit) const {
    std::memcpy(reinterpret_cast<uint8_t *>(&to) + edit.offset,
                reinterpret_cast<const uint8_t *>(&from) + edit.offset,
                edit.size);
  }

  void erase_pending_edit(size_t index) {
    for (size_t next = index + 1; next < pending_edit_count; ++next)
      pending_edits[next - 1] = pending_edits[next];
    --pending_edit_count;
  }

  void refresh_accepted_config() {
    if (admissible_config(requested_config)) {
      accepted_config = requested_config;
      pending_edit_count = 0;
      return;
    }

    Config candidate = requested_config;
    for (size_t index = 0; index < pending_edit_count; ++index)
      copy_pending_value(candidate, accepted_config, pending_edits[index]);
    if (admissible_config(candidate))
      accepted_config = candidate;

    bool admitted;
    do {
      admitted = false;
      for (size_t index = 0; index < pending_edit_count;) {
        candidate = accepted_config;
        copy_pending_value(candidate, requested_config, pending_edits[index]);
        if (!admissible_config(candidate)) {
          ++index;
          continue;
        }
        accepted_config = candidate;
        erase_pending_edit(index);
        admitted = true;
      }
    } while (admitted);
  }

  const char *parameter_warning(const char *name) const override {
    const ParamDef *parameter = getParameters().find(name);
    if (parameter != nullptr && parameter_out_of_range(*parameter))
      return begin_warning(
          "%s %.7g is outside its registered range [%.7g, %.7g]. Set %s "
          "within that range.",
          name, static_cast<double>(parameter->get_requested()),
          static_cast<double>(parameter->min),
          static_cast<double>(parameter->max), name);
    for (size_t index = 0; index < pending_edit_count; ++index) {
      const PendingEdit &edit = pending_edits[index];
      if (edit.name == nullptr || std::strcmp(edit.name, name) != 0)
        continue;
      if (schema_selector(name) && range_repairs_admission())
        return nullptr;
      return admission_warning(requested_config, edit.name);
    }
    return nullptr;
  }

  static bool parameter_out_of_range(const ParamDef &parameter) {
    const float value = parameter.get_requested();
    return value < parameter.min || value > parameter.max;
  }

  bool clamp_registered_parameter_ranges() {
    const uintptr_t requested = reinterpret_cast<uintptr_t>(&requested_config);
    bool clamped = false;
    for (const ParamDef &parameter : getParameters()) {
      if (parameter.is_bool() || parameter.is_enum() ||
          !parameter_out_of_range(parameter))
        continue;
      const uintptr_t target = reinterpret_cast<uintptr_t>(parameter.target);
      const size_t size = parameter_target_size(parameter);
      if (target < requested ||
          target + size > requested + sizeof(requested_config))
        continue;
      ParamDef writable = parameter;
      writable.set(
          hs::clamp(parameter.get_requested(), parameter.min, parameter.max));
      clamped = true;
    }
    return clamped;
  }

  bool range_repairs_admission() const {
    Config candidate = requested_config;
    const uintptr_t requested = reinterpret_cast<uintptr_t>(&requested_config);
    bool repaired = false;
    for (const ParamDef &parameter : getParameters()) {
      if (!parameter_out_of_range(parameter))
        continue;
      const uintptr_t target = reinterpret_cast<uintptr_t>(parameter.target);
      const size_t size = parameter_target_size(parameter);
      if (target < requested ||
          target + size > requested + sizeof(requested_config))
        continue;
      ParamDef candidate_parameter = parameter;
      candidate_parameter.target =
          reinterpret_cast<uint8_t *>(&candidate) + (target - requested);
      candidate_parameter.set(
          hs::clamp(parameter.get_requested(), parameter.min, parameter.max));
      repaired = true;
    }
    return repaired && admissible_config(candidate);
  }

  float accepted_parameter_value(const ParamDef &parameter) const override {
    const uintptr_t target = reinterpret_cast<uintptr_t>(parameter.target);
    const uintptr_t requested = reinterpret_cast<uintptr_t>(&requested_config);
    const size_t size = parameter_target_size(parameter);
    if (target < requested ||
        target + size > requested + sizeof(requested_config))
      return parameter.get_requested();
    const size_t offset = target - requested;
    return parameter.get_from(
        reinterpret_cast<const uint8_t *>(&accepted_config) + offset);
  }

  static bool schema_selector(const char *name) {
    return std::strcmp(name, "Function") == 0 ||
           std::strcmp(name, "Projection") == 0 ||
           std::strcmp(name, "Peirce Layout") == 0 ||
           std::strcmp(name, "Airocean Layout") == 0 ||
           std::strcmp(name, "Bonne Hemisphere") == 0 ||
           std::strcmp(name, "Gnomonic Hemisphere") == 0 ||
           std::strcmp(name, "Projection Frame") == 0 ||
           std::strcmp(name, "Surface Noise") == 0 ||
           std::strcmp(name, "Surface Noise Placement") == 0 ||
           std::strcmp(name, "Lens") == 0 ||
           std::strcmp(name, "Planar Warp 1") == 0 ||
           std::strcmp(name, "Planar Warp 1 Curl Integrator") == 0 ||
           std::strcmp(name, "Planar Warp 2") == 0 ||
           std::strcmp(name, "Planar Warp 2 Curl Integrator") == 0 ||
           std::strcmp(name, "Value Transfer") == 0 ||
           std::strcmp(name, "Coverage") == 0 ||
           std::strcmp(name, "Palette") == 0 ||
           std::strcmp(name, "Brightness Envelope") == 0 ||
           std::strcmp(name, "Hue Shift Mode") == 0;
  }
#endif

  HS_COLD_MEMBER void register_value_transfer_controls(ValueTransfer transfer,
                                                       ValueParams &params) {
    if (transfer == ValueTransfer::ISO_CONTOUR) {
      register_animated_param("Iso Level", &params.iso_level, 0.0f, 1.0f);
      register_animated_param("Iso Width", &params.iso_width, SOFTNESS_MIN,
                              0.5f);
    } else if (transfer == ValueTransfer::SMOOTH_BANDS) {
      if (!fixed_topology)
        register_animated_int_param("Band Count", &params.band_count, 1,
                                    BAND_COUNT_MAX);
      register_animated_param("Band Phase", &params.band_phase, 0.0f, TWO_PI_F);
    }
  }

  HS_COLD_MEMBER void register_coverage_controls(CoveragePolicy coverage,
                                                 ValueParams &params) {
    if (coverage == CoveragePolicy::VALUE_CUTOUT) {
      register_animated_param("Cutout Threshold", &params.cutout_threshold,
                              0.0f, 1.0f);
      register_animated_param("Cutout Softness", &params.cutout_softness,
                              SOFTNESS_MIN, 0.5f);
    } else if (coverage == CoveragePolicy::EDGE_FADE) {
      register_animated_param("Edge Fade Width", &params.edge_width, 0.0f,
                              0.5f);
    }
  }

  HS_COLD_MEMBER void register_stage_slot_controls(bool outer,
                                                   WarpStageSpec &spec) {
    if (spec.kind == WarpStageKind::VECTOR_NOISE ||
        spec.kind == WarpStageKind::CURL_FLOW) {
      register_animated_param(outer ? "Planar Warp 1 Noise Basis"
                                    : "Planar Warp 2 Noise Basis",
                              &spec.basis, NOISE_BASIS_OPTIONS,
                              NOISE_BASIS_EXPORT_OPTIONS, NUM_NOISE_BASES);
      register_animated_param(outer ? "Planar Warp 1 Envelope"
                                    : "Planar Warp 2 Envelope",
                              &spec.envelope, WARP_ENVELOPE_OPTIONS,
                              WARP_ENVELOPE_EXPORT_OPTIONS, NUM_WARP_ENVELOPES);
    }
    if (spec.kind == WarpStageKind::CURL_FLOW)
      register_animated_param(outer ? "Planar Warp 1 Curl Integrator"
                                    : "Planar Warp 2 Curl Integrator",
                              &spec.curl_integrator, CURL_INTEGRATOR_OPTIONS,
                              CURL_INTEGRATOR_EXPORT_OPTIONS,
                              NUM_CURL_INTEGRATORS);
    if (spec.kind == WarpStageKind::POLAR_CHART) {
      register_animated_param(outer ? "Planar Warp 1 Polar Mode"
                                    : "Planar Warp 2 Polar Mode",
                              &spec.polar_mode, POLAR_MODE_OPTIONS,
                              POLAR_MODE_EXPORT_OPTIONS, NUM_POLAR_MODES);
      register_animated_int_param(outer ? "Planar Warp 1 Polar Harmonic"
                                        : "Planar Warp 2 Polar Harmonic",
                                  &spec.polar_harmonic, 1, POLAR_HARMONIC_MAX);
    }
  }

  HS_COLD_MEMBER void register_source_controls(Function function,
                                               SourceParams &params,
                                               float domain_scale) {
    if (function == Function::SPHERICAL_RINGS) {
      register_animated_int_param("Ring Count", &params.ring_count, 1, 32);
      register_clamped_animated_param("Ring Thickness", &params.ring_thickness,
                                      1.0f / 512.0f, 0.5f);
      register_clamped_animated_param("Ring Softness", &params.ring_softness,
                                      SOFTNESS_MIN, 0.25f);
      register_clamped_animated_param("Ring Speed", &params.speed, -0.5f, 0.5f);
      register_clamped_animated_param("Ring Spin Speed", &params.angle_rate,
                                      -0.05f, 0.05f);
      register_clamped_animated_param("Ring Wander", &params.ring_wander, 0.0f,
                                      1.0f);
      return;
    }
    if (function == Function::FRACTAL) {
      register_clamped_animated_param("Fractal Scale", &params.fractal_scale,
                                      1.0f / 64.0f, 8.0f);
      register_animated_int_param("Fractal Iterations",
                                  &params.fractal_iterations, 2, 16);
      register_clamped_animated_param("Julia Mix", &params.julia_mix, 0.0f,
                                      1.0f);
      register_clamped_animated_param("Julia Real", &params.julia_real, -1.5f,
                                      1.5f);
      register_clamped_animated_param("Julia Imaginary",
                                      &params.julia_imaginary, -1.5f, 1.5f);
      register_clamped_animated_param("Fractal Contours",
                                      &params.fractal_contours, 0.0f, 16.0f);
      register_clamped_animated_param("Fractal Speed", &params.speed, -0.05f,
                                      0.05f);
      register_clamped_animated_param("Fractal Spin Speed", &params.angle_rate,
                                      -0.05f, 0.05f);
      return;
    }
    if (function == Function::TESSELLATION) {
      register_clamped_animated_param(
          "Cell Scale", &params.tessellation_cell_scale, 1.0f / 64.0f, 8.0f);
      register_clamped_animated_param("Line Thickness",
                                      &params.tessellation_line_thickness,
                                      SOFTNESS_MIN, 0.25f);
      register_clamped_animated_param("Line Softness",
                                      &params.tessellation_line_softness,
                                      SOFTNESS_MIN, 0.25f);
      register_clamped_animated_param("Tessellation Spin Speed",
                                      &params.angle_rate, -0.05f, 0.05f);
      if (!fixed_topology)
        register_animated_param("Tessellation Kind", &params.tessellation_kind,
                                TESSELLATION_KIND_OPTIONS,
                                TESSELLATION_KIND_EXPORT_OPTIONS,
                                NUM_TESSELLATION_KINDS);
      return;
    }
    if (is_noise_contour(function)) {
      register_clamped_animated_param(
          "Source Noise Scale", &params.noise_scale, SOURCE_NOISE_SCALE_MIN,
          domain_scaled_max(SOURCE_NOISE_SCALE_MAX, 0.5f, domain_scale));
      register_animated_param("Source Noise Contrast", &params.noise_contrast,
                              0.0f, 8.0f);
      register_clamped_animated_param(
          "Source Noise Speed", &params.noise_time_rate,
          -domain_scaled_max(SOURCE_NOISE_RATE_MAX, 1.0f / 4096.0f,
                             domain_scale),
          domain_scaled_max(SOURCE_NOISE_RATE_MAX, 1.0f / 4096.0f,
                            domain_scale));
      if (!fixed_topology)
        register_animated_param("Source Noise Basis", &params.noise_basis,
                                NOISE_BASIS_OPTIONS, NOISE_BASIS_EXPORT_OPTIONS,
                                NUM_NOISE_BASES);
      return;
    }
    if (function == Function::PRIMITIVE_LATTICE) {
      register_clamped_animated_param(
          "Lattice Cell Scale", &params.lattice_cell_scale, CELL_MIN, CELL_MAX);
      register_animated_param("Lattice Shape", &params.lattice_shape_blend,
                              0.0f, 1.0f);
      register_animated_param("Lattice Softness", &params.lattice_softness,
                              SOFTNESS_MIN, 1.0f);
      register_animated_param("Lattice Radius", &params.lattice_radius,
                              1.0f / 64.0f, 0.49f);
      return;
    }
    register_clamped_animated_param("Pattern Freq", &params.pattern_freq,
                                    PATTERN_FREQ_MIN,
                                    pattern_freq_max(function));
    register_clamped_animated_param(
        "Speed", &params.speed, 0.0f,
        domain_scaled_max(SPEED_MAX, 0.5f, domain_scale));
    register_clamped_animated_param(
        "Source Angle Speed", &params.angle_rate, 0.0f,
        domain_scaled_max(WAVE_SPIN_MAX, 0.03f, domain_scale));
    if (function == Function::GRID) {
      register_animated_param("Complexity", &params.complexity, COMPLEXITY_MIN,
                              COMPLEXITY_MAX);
      register_animated_param("Pattern Mix", &params.pattern_mix,
                              PATTERN_MIX_MIN, PATTERN_MIX_MAX);
      register_clamped_animated_param(
          "Drift", &params.secondary_rate, PHASE2_RATE_MIN,
          domain_scaled_max(PHASE2_RATE_MAX, 1.25f, domain_scale));
    }
  }

  HS_COLD_MEMBER void register_projection_controls(Slots &slots,
                                                   Params &params) {
    if (!fixed_topology && slots.projection == Projection::PEIRCE_QUINCUNCIAL)
      register_animated_param("Peirce Layout", &slots.peirce_layout,
                              PEIRCE_LAYOUT_OPTIONS,
                              PEIRCE_LAYOUT_EXPORT_OPTIONS, NUM_PEIRCE_LAYOUTS);
    if (!fixed_topology && slots.projection == Projection::AIROCEAN)
      register_animated_param(
          "Airocean Layout", &slots.airocean_layout, AIROCEAN_LAYOUT_OPTIONS,
          AIROCEAN_LAYOUT_EXPORT_OPTIONS, NUM_AIROCEAN_LAYOUTS);
    if (slots.projection == Projection::EQUIRECTANGULAR ||
        slots.projection == Projection::STEREOGRAPHIC ||
        slots.projection == Projection::GNOMONIC ||
        slots.projection == Projection::PEIRCE_QUINCUNCIAL)
      register_animated_param("Singularity Fade",
                              &params.projection.singularity_fade,
                              SINGULARITY_FADE_MIN, SINGULARITY_FADE_MAX);
    if (slots.projection == Projection::SINUSOIDAL ||
        slots.projection == Projection::EQUIRECTANGULAR ||
        strict_projection(slots.projection)) {
      register_animated_param("Central Meridian",
                              &params.projection.central_meridian, 0.0f,
                              TWO_PI_F);
    }
    if (strict_projection(slots.projection)) {
      register_animated_param("Projection Scale",
                              &params.projection.coordinate_scale, 0.25f, 4.0f);
    }
    if (!fixed_topology && slots.projection == Projection::BONNE)
      register_animated_param(
          "Bonne Hemisphere", &slots.bonne_hemisphere, BONNE_HEMISPHERE_OPTIONS,
          BONNE_HEMISPHERE_EXPORT_OPTIONS, NUM_BONNE_HEMISPHERES);
    if (!fixed_topology && slots.projection == Projection::GNOMONIC)
      register_animated_param("Gnomonic Hemisphere", &slots.gnomonic_hemisphere,
                              GNOMONIC_HEMISPHERE_OPTIONS,
                              GNOMONIC_HEMISPHERE_EXPORT_OPTIONS,
                              NUM_GNOMONIC_HEMISPHERES);
    if (slots.projection == Projection::BONNE)
      register_animated_param("Bonne Standard Parallel",
                              &params.projection.bonne_standard_parallel, 1e-3f,
                              0.5f * PI_F);
    if (slots.projection == Projection::PEIRCE_QUINCUNCIAL &&
        (slots.peirce_layout == PeirceLayout::HORIZONTAL ||
         slots.peirce_layout == PeirceLayout::VERTICAL))
      register_animated_param("Projection Layout Scroll",
                              &params.projection.layout_scroll, -1.0f, 1.0f);
  }

  HS_COLD_MEMBER void
  register_projection_frame_controls(ProjectionFramePolicy frame,
                                     Params &params, float domain_scale) {
    if (frame == ProjectionFramePolicy::SPIN_WANDER) {
      register_clamped_animated_param(
          "Projection Spin Speed", &params.projection.spin_rate, SPIN_RATE_MIN,
          domain_scaled_max(SPIN_RATE_MAX, 0.04f, domain_scale));
      register_animated_param("Projection Wander", &params.projection.wander,
                              WANDER_MIN, WANDER_MAX);
    }
  }

  HS_COLD_MEMBER void register_lens_controls(SurfaceLens lens,
                                             SurfaceLensParams &params) {
    if (lens == SurfaceLens::NONE)
      return;
    if (lens == SurfaceLens::MOBIUS) {
      register_animated_param("Mobius A Real", &params.mobius.a.re, -8.0f,
                              8.0f);
      register_animated_param("Mobius A Imag", &params.mobius.a.im, -8.0f,
                              8.0f);
      register_animated_param("Mobius B Real", &params.mobius.b.re, -8.0f,
                              8.0f);
      register_animated_param("Mobius B Imag", &params.mobius.b.im, -8.0f,
                              8.0f);
      register_animated_param("Mobius C Real", &params.mobius.c.re, -8.0f,
                              8.0f);
      register_animated_param("Mobius C Imag", &params.mobius.c.im, -8.0f,
                              8.0f);
      register_animated_param("Mobius D Real", &params.mobius.d.re, -8.0f,
                              8.0f);
      register_animated_param("Mobius D Imag", &params.mobius.d.im, -8.0f,
                              8.0f);
    }
  }

  HS_COLD_MEMBER void
  register_surface_noise_controls(Slots &slots, SurfaceNoiseParams &params,
                                  float domain_scale) {
    if (slots.surface_noise == SurfaceNoise::NONE)
      return;
    if (!fixed_topology) {
      register_animated_param(
          "Surface Noise Placement", &slots.surface_noise_placement,
          SURFACE_NOISE_PLACEMENT_OPTIONS,
          SURFACE_NOISE_PLACEMENT_EXPORT_OPTIONS, NUM_SURFACE_NOISE_PLACEMENTS);
      register_animated_param("Surface Noise Basis", &params.basis,
                              NOISE_BASIS_OPTIONS, NOISE_BASIS_EXPORT_OPTIONS,
                              NUM_NOISE_BASES);
    }
    register_clamped_animated_param("Surface Noise Scale", &params.scale,
                                    LENS_NOISE_SCALE_MIN, LENS_NOISE_SCALE_MAX);
    const float strength_min =
        slots.surface_noise == SurfaceNoise::CURL ? -0.5f : 0.0f;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    register_animated_param_preserving_value(
        "Surface Noise Strength", &params.strength, strength_min, 0.5f);
#else
    register_animated_param("Surface Noise Strength", &params.strength,
                            strength_min, 0.5f);
#endif
    const float speed_max =
        domain_scaled_max(NOISE_RATE_MAX, 0.002f, domain_scale);
    register_clamped_animated_param("Surface Noise Speed", &params.rate,
                                    -speed_max, speed_max);
    if (slots.surface_noise == SurfaceNoise::DIRECT)
      register_animated_param("Surface Noise Direction", &params.direction,
                              0.0f, 1.0f);
    else if (!fixed_topology)
      register_animated_param("Surface Noise Integrator", &params.integrator,
                              SURFACE_CURL_INTEGRATOR_OPTIONS,
                              SURFACE_CURL_INTEGRATOR_EXPORT_OPTIONS,
                              NUM_SURFACE_CURL_INTEGRATORS);
  }

  HS_COLD_MEMBER void register_active_warp_controls(bool outer,
                                                    const WarpStageSpec &spec,
                                                    WarpStageParams &params,
                                                    float domain_scale) {
    if (spec.kind == WarpStageKind::NONE)
      return;
    const char *const *names =
        outer ? OUTER_WARP_PARAM_NAMES : INNER_WARP_PARAM_NAMES;
    const char *speed_name =
        outer ? "Planar Warp 1 Speed" : "Planar Warp 2 Speed";
    auto register_current = [&](const char *name, float *target, float minimum,
                                float maximum) {
      register_clamped_animated_param(name, target, minimum, maximum);
    };
    if (spec.kind == WarpStageKind::WAVE_SHEAR ||
        spec.kind == WarpStageKind::VECTOR_NOISE ||
        spec.kind == WarpStageKind::CURL_FLOW) {
      const char *strength_name =
          outer ? "Planar Warp 1 Strength" : "Planar Warp 2 Strength";
      const bool signed_strength = spec.kind == WarpStageKind::WAVE_SHEAR ||
                                   spec.kind == WarpStageKind::CURL_FLOW;
      float strength_max = spec.kind == WarpStageKind::VECTOR_NOISE
                               ? VECTOR_WARP_STRENGTH_MAX
                               : 4.0f;
      if (spec.kind == WarpStageKind::CURL_FLOW)
        strength_max = curl_strength_limit(spec, params);
      register_current(strength_name, &params.strength,
                       signed_strength ? -strength_max : 0.0f, strength_max);
    }
    const float speed_max =
        domain_scaled_max(NOISE_SPEED_MAX, 0.005f, domain_scale);
    register_current(speed_name, &params.speed, -speed_max, speed_max);
    switch (spec.kind) {
    case WarpStageKind::AFFINE_FRAME: {
      const float snapped_x = roundf(params.translation_x);
      const float snapped_y = roundf(params.translation_y);
      registered_range_clamped |= snapped_x != params.translation_x ||
                                  snapped_y != params.translation_y;
      params.translation_x = snapped_x;
      params.translation_y = snapped_y;
      for (int index = 0; index < 6; ++index) {
        float *targets[] = {&params.translation_x, &params.translation_y,
                            &params.rotation,      &params.scale_x,
                            &params.scale_y,       &params.shear};
        const float minimum[] = {-4.0f, -4.0f, -TWO_PI_F, 0.25f, 0.25f, -0.75f};
        const float maximum[] = {4.0f, 4.0f, TWO_PI_F, 4.0f, 4.0f, 0.75f};
        register_current(names[Workbench::WARP_NAME_TRANSLATION_X + index],
                         targets[index], minimum[index], maximum[index]);
      }
      break;
    }
    case WarpStageKind::WAVE_SHEAR:
      register_current(names[Workbench::WARP_NAME_FREQUENCY], &params.frequency,
                       0.0f, domain_scaled_max(64.0f, 8.0f, domain_scale));
      register_current(names[Workbench::WARP_NAME_FIELD_ANGLE],
                       &params.field_angle, 0.0f, TWO_PI_F);
      break;
    case WarpStageKind::VORTEX:
      register_current(names[Workbench::WARP_NAME_CENTER_X], &params.center_x,
                       -4.0f, 4.0f);
      register_current(names[Workbench::WARP_NAME_CENTER_Y], &params.center_y,
                       -4.0f, 4.0f);
      register_current(names[Workbench::WARP_NAME_RADIUS], &params.radius,
                       1.0f / 64.0f, 8.0f);
      register_current(names[Workbench::WARP_NAME_TURNS], &params.turns, -4.0f,
                       4.0f);
      register_current(names[Workbench::WARP_NAME_CENTER_ORBIT],
                       &params.center_orbit_radius, 0.0f, 4.0f);
      break;
    case WarpStageKind::VECTOR_NOISE:
    case WarpStageKind::CURL_FLOW:
      register_current(outer ? "Planar Warp 1 Scale" : "Planar Warp 2 Scale",
                       &params.scale, 1.0f / 64.0f,
                       domain_scaled_max(spec.kind == WarpStageKind::CURL_FLOW
                                             ? CURL_WARP_SCALE_MAX
                                             : VECTOR_WARP_SCALE_MAX,
                                         1.0f, domain_scale));
      register_current(names[Workbench::WARP_NAME_VECTOR_ANGLE],
                       &params.vector_angle, 0.0f, TWO_PI_F);
      register_current(names[Workbench::WARP_NAME_EDGE_WIDTH],
                       &params.edge_width, SOFTNESS_MIN, 0.5f);
      break;
    case WarpStageKind::MIRROR_TILE:
      register_current(names[Workbench::WARP_NAME_ROTATION], &params.rotation,
                       0.0f, TWO_PI_F);
      register_current(names[Workbench::WARP_NAME_CELL_X], &params.cell_x,
                       CELL_MIN, CELL_MAX);
      register_current(names[Workbench::WARP_NAME_CELL_Y], &params.cell_y,
                       CELL_MIN, CELL_MAX);
      register_current(names[Workbench::WARP_NAME_OFFSET_X], &params.offset_x,
                       -8.0f, 8.0f);
      register_current(names[Workbench::WARP_NAME_OFFSET_Y], &params.offset_y,
                       -8.0f, 8.0f);
      break;
    case WarpStageKind::POLAR_CHART:
      register_current(names[Workbench::WARP_NAME_RADIAL_SCALE],
                       &params.radial_scale, 1.0f / 64.0f, 16.0f);
      register_current(names[Workbench::WARP_NAME_RADIAL_PHASE],
                       &params.radial_phase, 0.0f, TWO_PI_F);
      register_current(names[Workbench::WARP_NAME_ANGULAR_PHASE],
                       &params.angular_phase, 0.0f, TWO_PI_F);
      break;
    case WarpStageKind::NONE:
    case WarpStageKind::LEGACY_STEREO_NOISE:
      break;
    }
  }

  HS_COLD_MEMBER void register_color_controls(HueShiftMode mode,
                                              ColorParams &params,
                                              float domain_scale) {
    if (mode == HueShiftMode::NONE)
      return;
    register_clamped_animated_param("Hue Shift Amount",
                                    &params.hue_shift_amount,
                                    -Workbench::hue_shift_amount_max(mode),
                                    Workbench::hue_shift_amount_max(mode));
    if (mode != HueShiftMode::NOISE)
      return;
    register_clamped_animated_param(
        "Hue Noise Scale", &params.hue_noise_scale, HUE_NOISE_SCALE_MIN,
        domain_scaled_max(HUE_NOISE_SCALE_MAX, 2.0f, domain_scale));
    register_clamped_animated_param("Hue Noise Speed", &params.hue_noise_speed,
                                    -HUE_NOISE_SPEED_MAX, HUE_NOISE_SPEED_MAX);
  }

  struct ShaderWorkbenchInstrumentation {
#ifdef HS_PROFILE_SHADER_WORKBENCH_STAGES
    using Token = uint32_t;

    __attribute__((always_inline)) static Token mark() {
      return HS_OS_CYCLES();
    }

    template <Pullback::ProfileEvent Event>
    __attribute__((always_inline)) static void span(Token start) {
      const uint32_t elapsed = HS_OS_CYCLES() - start;
      if constexpr (Event == Pullback::ProfileEvent::LENS)
        hs::g_shader_workbench_stage_cycles.lens += elapsed;
      else if constexpr (Event == Pullback::ProfileEvent::SURFACE_NOISE)
        hs::g_shader_workbench_stage_cycles.surface_noise += elapsed;
      else if constexpr (Event == Pullback::ProfileEvent::PROJECTION)
        hs::g_shader_workbench_stage_cycles.projection += elapsed;
      else if constexpr (Event == Pullback::ProfileEvent::PLANAR_WARP)
        hs::g_shader_workbench_stage_cycles.planar_warp += elapsed;
      else if constexpr (Event == Pullback::ProfileEvent::MIRROR_TILE)
        hs::g_shader_workbench_stage_cycles.mirror_tile += elapsed;
      else if constexpr (Event == Pullback::ProfileEvent::SOURCE)
        hs::g_shader_workbench_stage_cycles.source += elapsed;
      else if constexpr (Event == Pullback::ProfileEvent::MATERIAL)
        hs::g_shader_workbench_stage_cycles.material += elapsed;
      else
        hs::g_shader_workbench_stage_cycles.color += elapsed;
    }
#else
    using Token = Pullback::NoInstrumentation::Token;

    __attribute__((always_inline)) static Token mark() { return {}; }

    template <Pullback::ProfileEvent Event>
    __attribute__((always_inline)) static void span(Token) {
      static_cast<void>(Event);
    }
#endif
  };

  struct ShaderWorkbenchBinding {
    using FrameState = ShaderWorkbench::FrameState;
    using Instrumentation = ShaderWorkbenchInstrumentation;

    template <typename Stage> static constexpr bool edge_unconditional() {
      if constexpr (requires { Stage::EDGE_DISTANCE_UNCONDITIONAL; })
        return Stage::EDGE_DISTANCE_UNCONDITIONAL;
      else
        return false;
    }

    template <typename Stage> static constexpr bool edge_fade_coverage() {
      if constexpr (requires { Stage::COVERAGE; })
        return Stage::COVERAGE == CoveragePolicy::EDGE_FADE;
      else
        return false;
    }

    /** EDGE_DISTANCE_UNCONDITIONAL on the projection requires an edge-fade
        coverage somewhere in the chain. */
    template <typename... Stages> struct ExtraValidation {
      static constexpr bool value = !(edge_unconditional<Stages>() || ...) ||
                                    (edge_fade_coverage<Stages>() || ...);
    };
  };

  struct ProjectionStateProvider {
    using Binding = ShaderWorkbenchBinding;
    using FrameState = typename Binding::FrameState;

    static const Quaternion &conjugate(const FrameState &frame) {
      return frame.transforms.projection_conj;
    }
    static float singularity_fade(const FrameState &frame) {
      return frame.params.projection.singularity_fade;
    }
    static float central_meridian(const FrameState &frame) {
      return frame.params.projection.central_meridian;
    }
    static float coordinate_scale(const FrameState &frame) {
      return frame.params.projection.coordinate_scale;
    }
    static float standard_parallel(const FrameState &frame) {
      return frame.params.projection.bonne_standard_parallel;
    }
    static float layout_scroll(const FrameState &frame) {
      return frame.params.projection.layout_scroll;
    }
  };

  struct SurfaceStateProvider {
    using Binding = ShaderWorkbenchBinding;
    using FrameState = typename Binding::FrameState;

    static const FastNoiseLite &noise(const FrameState &frame) {
      return *frame.resources.surface_noise;
    }
    static float scale(const FrameState &frame) {
      return frame.params.surface_noise.scale;
    }
    static float strength(const FrameState &frame) {
      return frame.params.surface_noise.strength;
    }
    static PreparedSurfaceNoise prepare(const FrameState &frame) {
      return prepare_surface_noise(frame.clocks, frame.params);
    }
    __attribute__((always_inline)) static bool
    path_length_required(const FrameState &frame) {
      return tracks_displacement(frame);
    }
  };

  struct LensStateProvider {
    using Binding = ShaderWorkbenchBinding;
    using FrameState = typename Binding::FrameState;

    static const MobiusParams &params(const FrameState &frame) {
      return frame.params.surface_lens.mobius;
    }
  };

  template <bool Outer> struct WarpStateProvider {
    using Binding = ShaderWorkbenchBinding;
    using FrameState = typename Binding::FrameState;

    static const WarpStageParams &params(const FrameState &frame) {
      if constexpr (Outer)
        return frame.params.warp.outer;
      else
        return frame.params.warp.inner;
    }
    static PreparedWarpStage prepare(const FrameState &frame) {
      const Complex period = source_cartesian_period(
          frame.slots.function, frame.params.source.lattice_cell_scale);
      if constexpr (Outer)
        return prepare_warp_stage(frame.slots.warp_program.outer,
                                  frame.params.warp.outer,
                                  frame.clocks.warp_outer_phase, period,
                                  frame.clocks.warp_outer_rotation);
      else
        return prepare_warp_stage(frame.slots.warp_program.inner,
                                  frame.params.warp.inner,
                                  frame.clocks.warp_inner_phase, period,
                                  frame.clocks.warp_inner_rotation);
    }
    static float phase(const FrameState &frame) {
      if constexpr (Outer)
        return frame.clocks.warp_outer_phase;
      else
        return frame.clocks.warp_inner_phase;
    }
    static const FastNoiseLite &noise(const FrameState &frame) {
      if constexpr (Outer)
        return *frame.resources.outer_warp_noise;
      else
        return *frame.resources.inner_warp_noise;
    }
    __attribute__((always_inline)) static bool
    path_length_required(const FrameState &frame) {
      return tracks_displacement(frame);
    }
  };

  struct SourceStateProvider {
    using Binding = ShaderWorkbenchBinding;
    using FrameState = typename Binding::FrameState;

    static const SourceParams &params(const FrameState &frame) {
      return frame.params.source;
    }
    static SourceState prepare(const FrameState &frame) {
      return prepare_source_state(frame.clocks);
    }
    static const FastNoiseLite &noise(const FrameState &frame) {
      return *frame.resources.source_noise;
    }
    static float noise_scale(const FrameState &frame) {
      return frame.params.source.noise_scale;
    }
    static float noise_time(const FrameState &frame) {
      return frame.clocks.source_noise_time;
    }
    static float noise_contrast(const FrameState &frame) {
      return frame.params.source.noise_contrast;
    }
  };

  struct ValueStateProvider {
    using Binding = ShaderWorkbenchBinding;
    using FrameState = typename Binding::FrameState;

    static float iso_level(const FrameState &frame) {
      return frame.params.value.iso_level;
    }
    static float iso_width(const FrameState &frame) {
      return frame.params.value.iso_width;
    }
    static float band_count(const FrameState &frame) {
      return static_cast<float>(frame.params.value.band_count);
    }
    static float band_phase(const FrameState &frame) {
      return frame.params.value.band_phase;
    }
    static float cutout_threshold(const FrameState &frame) {
      return frame.params.value.cutout_threshold;
    }
    static float cutout_softness(const FrameState &frame) {
      return frame.params.value.cutout_softness;
    }
    static float edge_width(const FrameState &frame) {
      return frame.params.value.edge_width;
    }
  };

  static_assert(
      static_cast<uint8_t>(PaletteMapping::CUP) ==
          static_cast<uint8_t>(Pullback::Color::PaletteMapping::CUP) &&
      static_cast<uint8_t>(PaletteMapping::BELL) ==
          static_cast<uint8_t>(Pullback::Color::PaletteMapping::BELL) &&
      static_cast<uint8_t>(PaletteMapping::LINEAR) ==
          static_cast<uint8_t>(Pullback::Color::PaletteMapping::LINEAR) &&
      static_cast<uint8_t>(PaletteMapping::REVERSE) ==
          static_cast<uint8_t>(Pullback::Color::PaletteMapping::REVERSE));

  static constexpr PaletteMappingWeights
  palette_mapping_weights(PaletteMapping mapping) {
    return PaletteMappingWeights::single(
        static_cast<Pullback::Color::PaletteMapping>(mapping));
  }
  static_assert(
      static_cast<uint8_t>(BrightnessEnvelope::NONE) ==
          static_cast<uint8_t>(Pullback::Color::BrightnessEnvelope::NONE) &&
      static_cast<uint8_t>(BrightnessEnvelope::CUP) ==
          static_cast<uint8_t>(Pullback::Color::BrightnessEnvelope::CUP) &&
      static_cast<uint8_t>(BrightnessEnvelope::BELL) ==
          static_cast<uint8_t>(Pullback::Color::BrightnessEnvelope::BELL) &&
      static_cast<uint8_t>(BrightnessEnvelope::ASCENDING) ==
          static_cast<uint8_t>(
              Pullback::Color::BrightnessEnvelope::ASCENDING) &&
      static_cast<uint8_t>(BrightnessEnvelope::DESCENDING) ==
          static_cast<uint8_t>(
              Pullback::Color::BrightnessEnvelope::DESCENDING));
  static_assert(
      static_cast<uint8_t>(HueShiftMode::NONE) ==
          static_cast<uint8_t>(Pullback::Color::HueMode::NONE) &&
      static_cast<uint8_t>(HueShiftMode::NOISE) ==
          static_cast<uint8_t>(Pullback::Color::HueMode::NOISE) &&
      static_cast<uint8_t>(HueShiftMode::WARP_DISPLACEMENT) ==
          static_cast<uint8_t>(Pullback::Color::HueMode::PATH_LENGTH));

  struct ColorStateProvider {
    using Binding = ShaderWorkbenchBinding;
    using FrameState = typename Binding::FrameState;

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

  using CodeEmission = Pullback::CodeEmission;
  using ApproximationOracleId = Pullback::ApproximationOracleId;
  using ApproximationDomain = Pullback::ApproximationDomain;
  using ApproximationAggregation = Pullback::ApproximationAggregation;
  using ApproximationMetric = Pullback::ApproximationMetric;

  struct TopologyKey {
    Function function{};
    Projection projection{};
    ProjectionFramePolicy projection_frame{};
    SurfaceLens surface_lens{};
    SignalWeight signal_weight{};
    ValueTransfer value_transfer{};
    CoveragePolicy coverage{};
    PeirceLayout peirce_layout{};
    AiroceanLayout airocean_layout{};
    BonneHemisphere bonne_hemisphere{};
    GnomonicHemispherePolicy gnomonic_hemisphere{};
    SurfaceNoise surface_noise{};
    SurfaceNoisePlacement surface_noise_placement{};
    NoiseBasis surface_noise_basis{};
    SurfaceCurlIntegrator surface_curl_integrator{};
    NoiseBasis source_noise_basis{};
    WarpStageKind outer_warp{};
    NoiseBasis outer_warp_basis{};
    WarpEnvelope outer_warp_envelope{};
    PolarMode outer_polar_mode{};
    CurlIntegrator outer_curl_integrator{};
    uint8_t outer_polar_harmonic{};
    WarpStageKind inner_warp{};
    NoiseBasis inner_warp_basis{};
    WarpEnvelope inner_warp_envelope{};
    PolarMode inner_polar_mode{};
    CurlIntegrator inner_curl_integrator{};
    uint8_t inner_polar_harmonic{};

    constexpr bool operator==(const TopologyKey &) const = default;
  };

  size_t preset_count_for_view() const {
    return preset_view.empty() ? PRESETS.size() : preset_view.size();
  }

  const Preset &preset_for_view(size_t index) const {
    HS_CHECK(index < preset_count_for_view(),
             "preset_for_view: index out of range");
    return PRESETS[preset_view.empty() ? index : preset_view[index]];
  }

  template <typename... Stages>
  using InversePipeline = Pullback::Pipeline<ShaderWorkbenchBinding, Stages...>;

#if defined(__IMXRT1062__)
  HS_FLASH_MEMBER
#else
  __attribute__((always_inline))
#endif
  static Vector pullback_outer_camera_lookup(const Vector &input,
                                             const FrameState &frame) {
    return rotate(input, frame.transforms.outer_conj);
  }

  struct OuterCameraStage
      : Pullback::Stage::Contract<OuterCameraStage, Pullback::SphereSample,
                                  Pullback::SphereSample> {
    using Policies = std::tuple<>;

    static constexpr bool implements(const TopologyKey &) { return true; }

    template <typename Binding>
    __attribute__((always_inline)) static Pullback::SphereSample
    run(const Pullback::SphereSample &input,
        const typename Binding::FrameState &frame,
        const Pullback::NoPrepared &) {
      return {pullback_outer_camera_lookup(input.dir, frame),
              input.path_length};
    }
  };

  template <SurfaceLens LensV>
  using LensPolicy = std::conditional_t<
      LensV == SurfaceLens::NONE, void,
      std::conditional_t<
          LensV == SurfaceLens::GLITCH, Pullback::Lens::Glitch,
          std::conditional_t<
              LensV == SurfaceLens::KALEIDOSCOPE, Pullback::Lens::Kaleidoscope,
              std::conditional_t<
                  LensV == SurfaceLens::MOBIUS,
                  Pullback::Lens::Mobius<LensStateProvider>,
                  std::conditional_t<
                      LensV == SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
                      Pullback::Lens::DodecahedralKaleidoscope,
                      std::conditional_t<
                          LensV == SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM,
                          Pullback::Lens::HexagonalPrismKaleidoscope,
                          std::conditional_t<
                              LensV ==
                                  SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM,
                              Pullback::Lens::PentagonalPrismKaleidoscope,
                              Pullback::Lens::
                                  TriangularPrismKaleidoscope>>>>>>>;

  template <Projection ProjectionV>
  using ProjectionPolicy = std::conditional_t<
      ProjectionV == Projection::STEREOGRAPHIC,
      Pullback::Projection::Stereographic<ProjectionStateProvider>,
      std::conditional_t<
          ProjectionV == Projection::GNOMONIC,
          Pullback::Projection::Gnomonic<
              ProjectionStateProvider,
              Pullback::Projection::GnomonicHemisphere::FOLDED>,
          std::conditional_t<
              ProjectionV == Projection::BONNE,
              Pullback::Projection::Bonne<ProjectionStateProvider, true>,
              std::conditional_t<ProjectionV == Projection::EQUIRECTANGULAR,
                                 Pullback::Projection::Equirectangular<
                                     ProjectionStateProvider>,
                                 Pullback::Projection::PeirceSquare<
                                     ProjectionStateProvider>>>>>;

  template <SurfaceLens LensV>
  struct SelectedLensStage : Pullback::Stage::Lens<LensPolicy<LensV>> {
    static_assert(LensV == SurfaceLens::GLITCH ||
                  LensV == SurfaceLens::KALEIDOSCOPE ||
                  LensV == SurfaceLens::MOBIUS ||
                  LensV == SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL ||
                  LensV == SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM ||
                  LensV == SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM ||
                  LensV == SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM);
    static constexpr bool implements(const TopologyKey &key) {
      return key.surface_lens == LensV;
    }
  };

  template <Projection ProjectionV, SurfaceLens LensV>
  struct SelectedProjectStage
      : Pullback::Stage::Project<ProjectionPolicy<ProjectionV>> {
    static_assert(ProjectionV == Projection::STEREOGRAPHIC ||
                  ProjectionV == Projection::GNOMONIC ||
                  ProjectionV == Projection::BONNE ||
                  ProjectionV == Projection::EQUIRECTANGULAR ||
                  ProjectionV == Projection::PEIRCE_QUINCUNCIAL);
    static constexpr bool implements(const TopologyKey &key) {
      return key.projection == ProjectionV && key.surface_lens == LensV &&
             key.surface_noise == SurfaceNoise::NONE &&
             (ProjectionV != Projection::GNOMONIC ||
              key.gnomonic_hemisphere == GnomonicHemispherePolicy::FOLDED) &&
             (ProjectionV != Projection::PEIRCE_QUINCUNCIAL ||
              key.peirce_layout == PeirceLayout::SQUARE) &&
             (ProjectionV != Projection::BONNE ||
              key.bonne_hemisphere == BonneHemisphere::NORTH);
    }
  };

  struct SinusoidalCurlDisplaceStage
      : Pullback::Stage::Displace<Pullback::Surface::CurlNoise<
            SurfaceStateProvider, NoiseBasis::SIMPLEX,
            Pullback::Surface::Euler>> {
    static constexpr bool implements(const TopologyKey &key) {
      return key.surface_noise == SurfaceNoise::CURL &&
             key.surface_noise_placement ==
                 SurfaceNoisePlacement::BEFORE_LENS &&
             key.surface_noise_basis == NoiseBasis::SIMPLEX &&
             key.surface_curl_integrator == SurfaceCurlIntegrator::EULER;
    }
  };

  struct SinusoidalCurlProjectStage
      : Pullback::Stage::Project<
            Pullback::Projection::FoldedSinusoidal<ProjectionStateProvider>> {
    static constexpr bool implements(const TopologyKey &key) {
      return key.projection == Projection::SINUSOIDAL &&
             key.surface_lens == SurfaceLens::NONE;
    }
  };

  using SinusoidalCurlSphereRun =
      Pullback::Stage::Placed<Pullback::CodeEmission::OUT_OF_LINE_FLASH,
                              SinusoidalCurlDisplaceStage,
                              SinusoidalCurlProjectStage>;

  template <WarpStageKind KindV, bool Outer>
  using WarpPolicy = std::conditional_t<
      KindV == WarpStageKind::NONE, void,
      std::conditional_t<
          KindV == WarpStageKind::AFFINE_FRAME,
          Pullback::Warp::AffineFrame<WarpStateProvider<Outer>>,
          std::conditional_t<
              KindV == WarpStageKind::WAVE_SHEAR,
              Pullback::Warp::WaveShear<WarpStateProvider<Outer>>,
              std::conditional_t<
                  KindV == WarpStageKind::VECTOR_NOISE,
                  Pullback::Warp::VectorNoise<WarpStateProvider<Outer>,
                                              NoiseBasis::SIMPLEX,
                                              Pullback::Warp::FlatEnvelope>,
                  std::conditional_t<
                      KindV == WarpStageKind::MIRROR_TILE,
                      Pullback::Warp::MirrorTile<WarpStateProvider<Outer>>,
                      Pullback::Warp::PolarChart<WarpStateProvider<Outer>,
                                                 Pullback::Warp::LinearPolar,
                                                 1>>>>>>;

  template <WarpStageKind KindV, bool Outer>
  struct SelectedWarpStage : Pullback::Stage::Warp<WarpPolicy<KindV, Outer>> {
    static_assert(KindV == WarpStageKind::AFFINE_FRAME ||
                  KindV == WarpStageKind::WAVE_SHEAR ||
                  KindV == WarpStageKind::VECTOR_NOISE ||
                  KindV == WarpStageKind::MIRROR_TILE ||
                  KindV == WarpStageKind::POLAR_CHART);
    static_assert(Outer || KindV == WarpStageKind::WAVE_SHEAR ||
                  KindV == WarpStageKind::MIRROR_TILE);
    static constexpr bool implements(const TopologyKey &key) {
      if constexpr (Outer)
        return key.outer_warp == KindV &&
               (KindV != WarpStageKind::WAVE_SHEAR ||
                key.outer_warp_envelope == WarpEnvelope::FLAT);
      else
        return key.inner_warp == KindV &&
               (KindV != WarpStageKind::WAVE_SHEAR ||
                key.inner_warp_envelope == WarpEnvelope::FLAT);
    }
  };

  template <Function FunctionV>
  using SourcePolicy = std::conditional_t<
      FunctionV == Function::GRID, Pullback::Source::Grid<SourceStateProvider>,
      std::conditional_t<
          FunctionV == Function::PRIMITIVE_LATTICE,
          Pullback::Source::PrimitiveLattice<SourceStateProvider>,
          Pullback::Source::TwinWave<SourceStateProvider>>>;

  template <CoveragePolicy CoverageV> struct ProjectionCoverageMapping {
    static_assert(CoverageV == CoveragePolicy::OPAQUE ||
                  CoverageV == CoveragePolicy::PROJECTION_WEIGHT ||
                  CoverageV == CoveragePolicy::PROJECTION_WEIGHT_SQUARED ||
                  CoverageV == CoveragePolicy::EDGE_FADE);
    static constexpr Pullback::ProjectionCoverageMode MODE =
        CoverageV == CoveragePolicy::OPAQUE
            ? Pullback::ProjectionCoverageMode::NONE
        : CoverageV == CoveragePolicy::PROJECTION_WEIGHT
            ? Pullback::ProjectionCoverageMode::WEIGHT
        : CoverageV == CoveragePolicy::PROJECTION_WEIGHT_SQUARED
            ? Pullback::ProjectionCoverageMode::WEIGHT_SQUARED
            : Pullback::ProjectionCoverageMode::EDGE_FADE;
    using Type = std::conditional_t<
        MODE == Pullback::ProjectionCoverageMode::NONE,
        Pullback::ProjectionCoverage::None,
        std::conditional_t<
            MODE == Pullback::ProjectionCoverageMode::WEIGHT,
            Pullback::ProjectionCoverage::Weight,
            std::conditional_t<
                MODE == Pullback::ProjectionCoverageMode::WEIGHT_SQUARED,
                Pullback::ProjectionCoverage::WeightSquared,
                Pullback::ProjectionCoverage::EdgeFade<ValueStateProvider>>>>;
  };

  template <CoveragePolicy CoverageV>
  using CoveragePolicyFor = typename ProjectionCoverageMapping<CoverageV>::Type;

  template <Function FunctionV, ValueTransfer TransferV,
            CoveragePolicy CoverageV>
  struct SelectedSampleStage
      : Pullback::Stage::Sample<SourcePolicy<FunctionV>,
                                Pullback::Weight::Projection,
                                CoveragePolicyFor<CoverageV>> {
    static_assert(FunctionV == Function::GRID ||
                  FunctionV == Function::PRIMITIVE_LATTICE ||
                  FunctionV == Function::TWIN_WAVE);
    static_assert(CoverageV == CoveragePolicy::OPAQUE ||
                  CoverageV == CoveragePolicy::EDGE_FADE ||
                  CoverageV == CoveragePolicy::PROJECTION_WEIGHT_SQUARED ||
                  CoverageV == CoveragePolicy::PROJECTION_WEIGHT);
    static constexpr CoveragePolicy COVERAGE = CoverageV;
    static constexpr bool implements(const TopologyKey &key) {
      return key.function == FunctionV &&
             key.signal_weight == SignalWeight::PROJECTION &&
             key.value_transfer == TransferV && key.coverage == CoverageV;
    }
  };

  struct IsoContourTransferStage
      : Pullback::Stage::Transfer<
            Pullback::Transfer::IsoContour<ValueStateProvider>> {
    static constexpr bool implements(const TopologyKey &key) {
      return key.value_transfer == ValueTransfer::ISO_CONTOUR;
    }
  };

  struct ColorStage
      : Pullback::Stage::Colorize<
            Pullback::Color::GeneratedPalette<ColorStateProvider>> {
    static constexpr bool implements(const TopologyKey &) { return true; }
  };

  using GlitchNoiseGridWaveShearPipelineBase = InversePipeline<
      OuterCameraStage, SelectedLensStage<SurfaceLens::GLITCH>,
      SelectedProjectStage<Projection::STEREOGRAPHIC, SurfaceLens::GLITCH>,
      SelectedWarpStage<WarpStageKind::WAVE_SHEAR, true>,
      SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                          CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  struct GlitchNoiseGridWaveShearPipeline
      : GlitchNoiseGridWaveShearPipelineBase {
    HS_HOT_FLASH_MEMBER static Color4 shade_prepared(const Vector &view,
                                                     const FrameState &frame,
                                                     const void *storage) {
      return GlitchNoiseGridWaveShearPipelineBase::evaluate(
          view, frame,
          *static_cast<const typename GlitchNoiseGridWaveShearPipelineBase::
                           PreparedTuple *>(storage));
    }
  };
  using KaleidoscopeTwinWaveInnerMirrorPipeline = InversePipeline<
      OuterCameraStage, SelectedLensStage<SurfaceLens::KALEIDOSCOPE>,
      SelectedProjectStage<Projection::STEREOGRAPHIC,
                           SurfaceLens::KALEIDOSCOPE>,
      SelectedWarpStage<WarpStageKind::MIRROR_TILE, false>,
      SelectedSampleStage<Function::TWIN_WAVE, ValueTransfer::NONE,
                          CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  using StereographicMobiusTwinWaveInnerMirrorPipeline = InversePipeline<
      OuterCameraStage, SelectedLensStage<SurfaceLens::MOBIUS>,
      SelectedProjectStage<Projection::STEREOGRAPHIC, SurfaceLens::MOBIUS>,
      SelectedWarpStage<WarpStageKind::MIRROR_TILE, false>,
      SelectedSampleStage<Function::TWIN_WAVE, ValueTransfer::NONE,
                          CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  using StereographicHexagonalPrismTwinWaveInnerMirrorPipeline =
      InversePipeline<
          OuterCameraStage,
          SelectedLensStage<SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM>,
          SelectedProjectStage<Projection::STEREOGRAPHIC,
                               SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM>,
          SelectedWarpStage<WarpStageKind::MIRROR_TILE, false>,
          SelectedSampleStage<Function::TWIN_WAVE, ValueTransfer::NONE,
                              CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
          ColorStage>;
  using GnomonicKaleidoscopeGridMirrorPipeline = InversePipeline<
      OuterCameraStage, SelectedLensStage<SurfaceLens::KALEIDOSCOPE>,
      SelectedProjectStage<Projection::GNOMONIC, SurfaceLens::KALEIDOSCOPE>,
      SelectedWarpStage<WarpStageKind::MIRROR_TILE, true>,
      SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                          CoveragePolicy::EDGE_FADE>,
      ColorStage>;
  using GnomonicAlienCoreMirrorPipeline = InversePipeline<
      OuterCameraStage, SelectedLensStage<SurfaceLens::GLITCH>,
      SelectedProjectStage<Projection::GNOMONIC, SurfaceLens::GLITCH>,
      SelectedWarpStage<WarpStageKind::MIRROR_TILE, true>,
      SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                          CoveragePolicy::EDGE_FADE>,
      ColorStage>;
  using PeirceDodecahedralGridPipelineBase = InversePipeline<
      OuterCameraStage,
      SelectedLensStage<SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
      SelectedProjectStage<Projection::PEIRCE_QUINCUNCIAL,
                           SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
      SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                          CoveragePolicy::EDGE_FADE>,
      ColorStage>;
  struct PeirceDodecahedralGridPipeline : PeirceDodecahedralGridPipelineBase {
    HS_HOT_FLASH_MEMBER static Color4 shade_prepared(const Vector &view,
                                                     const FrameState &frame,
                                                     const void *storage) {
      return PeirceDodecahedralGridPipelineBase::evaluate(
          view, frame,
          *static_cast<const typename PeirceDodecahedralGridPipelineBase::
                           PreparedTuple *>(storage));
    }
  };
  using GnomonicDodecahedralGridWaveMirrorPipelineBase = InversePipeline<
      OuterCameraStage,
      SelectedLensStage<SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
      SelectedProjectStage<Projection::GNOMONIC,
                           SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
      SelectedWarpStage<WarpStageKind::WAVE_SHEAR, true>,
      SelectedWarpStage<WarpStageKind::MIRROR_TILE, false>,
      SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                          CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  struct GnomonicDodecahedralGridWaveMirrorPipeline
      : GnomonicDodecahedralGridWaveMirrorPipelineBase {
    HS_HOT_FLASH_MEMBER static Color4 shade_prepared(const Vector &view,
                                                     const FrameState &frame,
                                                     const void *storage) {
      return GnomonicDodecahedralGridWaveMirrorPipelineBase::evaluate(
          view, frame,
          *static_cast<
              const typename GnomonicDodecahedralGridWaveMirrorPipelineBase::
                  PreparedTuple *>(storage));
    }
  };
  using GnomonicDodecahedralGridVectorMirrorPipeline = InversePipeline<
      OuterCameraStage,
      SelectedLensStage<SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
      SelectedProjectStage<Projection::GNOMONIC,
                           SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
      SelectedWarpStage<WarpStageKind::VECTOR_NOISE, true>,
      SelectedWarpStage<WarpStageKind::MIRROR_TILE, false>,
      SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                          CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  using GnomonicAffineLatticeContourPipeline = InversePipeline<
      OuterCameraStage,
      SelectedProjectStage<Projection::GNOMONIC, SurfaceLens::NONE>,
      SelectedWarpStage<WarpStageKind::AFFINE_FRAME, true>,
      SelectedSampleStage<Function::PRIMITIVE_LATTICE,
                          ValueTransfer::ISO_CONTOUR,
                          CoveragePolicy::PROJECTION_WEIGHT>,
      IsoContourTransferStage, ColorStage>;
  using SinusoidalLatticeMeltPipeline = InversePipeline<
      OuterCameraStage, SinusoidalCurlSphereRun,
      SelectedSampleStage<Function::PRIMITIVE_LATTICE, ValueTransfer::NONE,
                          CoveragePolicy::PROJECTION_WEIGHT>,
      ColorStage>;
  using StereographicPrismPolarWaveLatticePipeline = InversePipeline<
      OuterCameraStage,
      SelectedLensStage<SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM>,
      SelectedProjectStage<Projection::STEREOGRAPHIC,
                           SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM>,
      SelectedWarpStage<WarpStageKind::POLAR_CHART, true>,
      SelectedWarpStage<WarpStageKind::WAVE_SHEAR, false>,
      SelectedSampleStage<Function::PRIMITIVE_LATTICE, ValueTransfer::NONE,
                          CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  using StereographicDodecahedralGridInnerMirrorPipeline = InversePipeline<
      OuterCameraStage,
      SelectedLensStage<SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
      SelectedProjectStage<Projection::STEREOGRAPHIC,
                           SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
      SelectedWarpStage<WarpStageKind::MIRROR_TILE, false>,
      SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                          CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  using EquirectangularDodecahedralGridInnerMirrorPipeline = InversePipeline<
      OuterCameraStage,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::OUT_OF_LINE_FLASH,
          SelectedLensStage<SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
          SelectedProjectStage<Projection::EQUIRECTANGULAR,
                               SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>>,
      SelectedWarpStage<WarpStageKind::MIRROR_TILE, false>,
      SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                          CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  using StereographicAlienCoreMirrorPipeline = InversePipeline<
      OuterCameraStage, SelectedLensStage<SurfaceLens::GLITCH>,
      SelectedProjectStage<Projection::STEREOGRAPHIC, SurfaceLens::GLITCH>,
      SelectedWarpStage<WarpStageKind::MIRROR_TILE, true>,
      SelectedSampleStage<Function::GRID, ValueTransfer::NONE,
                          CoveragePolicy::EDGE_FADE>,
      ColorStage>;

  struct ProgramDescriptor {
    InversePipelineId id;
    TopologyKey key;
    ShadeFunction shade;
    void (*prepare)(const FrameState &, void *);
    bool (*continuous_parameters_supported)(const Config &);
    bool (*resources_ready)(const FrameState &);
  };

  /** Capacity of the per-frame prepared blob a compiled program renders
      from; make_program() pins every program's tuple under it. */
  static constexpr size_t PREPARED_BLOB_BYTES = 256;

  static constexpr void canonicalize_warp_key(WarpStageKind kind,
                                              NoiseBasis &basis,
                                              WarpEnvelope &envelope,
                                              PolarMode &polar_mode,
                                              CurlIntegrator &curl_integrator,
                                              uint8_t &polar_harmonic) {
    if (!warp_uses_noise(kind))
      basis = {};
    if (!warp_uses_envelope(kind))
      envelope = {};
    if (kind != WarpStageKind::CURL_FLOW)
      curl_integrator = {};
    if (kind != WarpStageKind::POLAR_CHART) {
      polar_mode = {};
      polar_harmonic = 0;
    }
  }

  static constexpr TopologyKey make_topology_key(const Config &config) {
    const Slots &slots = config.slots;
    TopologyKey key{
        slots.function,
        slots.projection,
        slots.projection_frame,
        slots.surface_lens,
        slots.signal_weight,
        slots.value_transfer,
        slots.coverage,
        slots.peirce_layout,
        slots.airocean_layout,
        slots.bonne_hemisphere,
        slots.gnomonic_hemisphere,
        slots.surface_noise,
        slots.surface_noise_placement,
        config.params.surface_noise.basis,
        config.params.surface_noise.integrator,
        config.params.source.noise_basis,
        slots.warp_program.outer.kind,
        slots.warp_program.outer.basis,
        slots.warp_program.outer.envelope,
        slots.warp_program.outer.polar_mode,
        slots.warp_program.outer.curl_integrator,
        slots.warp_program.outer.polar_harmonic,
        slots.warp_program.inner.kind,
        slots.warp_program.inner.basis,
        slots.warp_program.inner.envelope,
        slots.warp_program.inner.polar_mode,
        slots.warp_program.inner.curl_integrator,
        slots.warp_program.inner.polar_harmonic,
    };
    if (key.projection != Projection::PEIRCE_QUINCUNCIAL)
      key.peirce_layout = {};
    if (key.projection != Projection::AIROCEAN)
      key.airocean_layout = {};
    if (key.projection != Projection::BONNE)
      key.bonne_hemisphere = {};
    if (key.projection != Projection::GNOMONIC)
      key.gnomonic_hemisphere = {};
    if (key.surface_noise == SurfaceNoise::NONE) {
      key.surface_noise_placement = {};
      key.surface_noise_basis = {};
      key.surface_curl_integrator = {};
    }
    if (!source_uses_noise(key.function))
      key.source_noise_basis = {};
    canonicalize_warp_key(key.outer_warp, key.outer_warp_basis,
                          key.outer_warp_envelope, key.outer_polar_mode,
                          key.outer_curl_integrator, key.outer_polar_harmonic);
    canonicalize_warp_key(key.inner_warp, key.inner_warp_basis,
                          key.inner_warp_envelope, key.inner_polar_mode,
                          key.inner_curl_integrator, key.inner_polar_harmonic);
    return key;
  }

  static constexpr bool all_continuous_parameters_supported(const Config &) {
    return true;
  }

  HS_FLASH_MEMBER static bool
  pipeline_resources_ready(const FrameState &frame) {
    if (warp_uses_noise(frame.slots.warp_program.outer.kind) &&
        frame.resources.outer_warp_noise == nullptr)
      return false;
    if (warp_uses_noise(frame.slots.warp_program.inner.kind) &&
        frame.resources.inner_warp_noise == nullptr)
      return false;
    if (is_noise_contour(frame.slots.function) &&
        frame.resources.source_noise == nullptr)
      return false;
    if (frame.slots.surface_noise != SurfaceNoise::NONE &&
        frame.resources.surface_noise == nullptr)
      return false;
    if (frame.resources.generated_palette == nullptr)
      return false;
    if (frame.slots.hue_shift == HueShiftMode::NOISE &&
        frame.params.color.hue_shift_amount != 0.0f &&
        frame.resources.color_noise == nullptr)
      return false;
    if (frame.prepared_hue_rotation.active &&
        frame.prepared_hue_rotation.lut == nullptr)
      return false;
    return !frame.prepared_hue_noise.active ||
           frame.prepared_hue_noise.lut != nullptr;
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

  enum class ProfileEndpoint : uint8_t { STEADY, FROM, TO };

  struct EndpointRuntime {
    ClockState clocks{};
    Quaternion projection_wander;
    Quaternion outer_wander;
    Quaternion source_wander;
    PreparedTransforms transforms;

    HS_COLD_MEMBER EndpointRuntime() = default;
  };

  template <typename T> static uint32_t encode_field_value(const T &value) {
    static_assert(sizeof(T) <= sizeof(uint32_t));
    uint32_t payload = 0;
    std::memcpy(&payload, &value, sizeof(T));
    return payload;
  }

  template <typename T>
  static bool decode_field_value(uint32_t payload, T &value) {
    static_assert(sizeof(T) <= sizeof(uint32_t));
    if constexpr (sizeof(T) < sizeof(uint32_t)) {
      const uint32_t value_mask = (uint32_t{1} << (sizeof(T) * 8)) - 1;
      if ((payload & ~value_mask) != 0)
        return false;
    }
    if constexpr (std::is_same_v<T, bool>)
      if (payload > 1)
        return false;
    std::memcpy(&value, &payload, sizeof(T));
    return true;
  }

  static ConfigValues encode_config_values(const Config &config) {
    ConfigValues values{};
#define HS_SHADER_WORKBENCH_ENCODE_FIELD(name, path)                           \
  values[static_cast<size_t>(ConfigFieldId::name)] =                           \
      encode_field_value(config.path);
    HS_SHADER_WORKBENCH_CONFIG_FIELDS(HS_SHADER_WORKBENCH_ENCODE_FIELD)
#undef HS_SHADER_WORKBENCH_ENCODE_FIELD
    values[static_cast<size_t>(ConfigFieldId::SLOTS_SURFACE_LENS)] =
        surface_lens_storage_id(config.slots.surface_lens);
    values[static_cast<size_t>(ConfigFieldId::SLOTS_WARP_OUTER_KIND)] =
        warp_storage_id(config.slots.warp_program.outer.kind);
    values[static_cast<size_t>(ConfigFieldId::SLOTS_WARP_INNER_KIND)] =
        warp_storage_id(config.slots.warp_program.inner.kind);
    return values;
  }

  static bool decode_config_values(const ConfigValues &values, Config &config) {
    bool valid = true;
#define HS_SHADER_WORKBENCH_DECODE_FIELD(name, path)                           \
  valid = decode_field_value(values[static_cast<size_t>(ConfigFieldId::name)], \
                             config.path) &&                                   \
          valid;
    HS_SHADER_WORKBENCH_CONFIG_FIELDS(HS_SHADER_WORKBENCH_DECODE_FIELD)
#undef HS_SHADER_WORKBENCH_DECODE_FIELD
    valid = decode_surface_lens_storage(
                values[static_cast<size_t>(ConfigFieldId::SLOTS_SURFACE_LENS)],
                config.slots.surface_lens) &&
            valid;
    valid =
        decode_warp_storage(
            values[static_cast<size_t>(ConfigFieldId::SLOTS_WARP_OUTER_KIND)],
            config.slots.warp_program.outer.kind) &&
        valid;
    valid =
        decode_warp_storage(
            values[static_cast<size_t>(ConfigFieldId::SLOTS_WARP_INNER_KIND)],
            config.slots.warp_program.inner.kind) &&
        valid;
    return valid;
  }

  static constexpr uint32_t surface_lens_storage_id(SurfaceLens lens) {
    const uint8_t value = static_cast<uint8_t>(lens);
    if (lens == SurfaceLens::TANGENT_NOISE)
      return 5;
    return value < 5 ? value : value + 1;
  }

  static bool decode_surface_lens_storage(uint32_t id, SurfaceLens &lens) {
    if (id <= 4) {
      lens = static_cast<SurfaceLens>(id);
      return true;
    }
    if (id == 5) {
      lens = SurfaceLens::TANGENT_NOISE;
      return true;
    }
    if (id <= 13) {
      lens = static_cast<SurfaceLens>(id - 1);
      return true;
    }
    return false;
  }

  static constexpr uint32_t warp_storage_id(WarpStageKind kind) {
    if (kind == WarpStageKind::LEGACY_STEREO_NOISE)
      return 1;
    const uint8_t value = static_cast<uint8_t>(kind);
    return value == 0 ? 0 : value + 1;
  }

  static bool decode_warp_storage(uint32_t id, WarpStageKind &kind) {
    if (id == 0) {
      kind = WarpStageKind::NONE;
      return true;
    }
    if (id == 1) {
      kind = WarpStageKind::LEGACY_STEREO_NOISE;
      return true;
    }
    if (id <= 8) {
      kind = static_cast<WarpStageKind>(id - 1);
      return true;
    }
    return false;
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
public:
  /** @brief Captures all accepted, requested, pending, and runtime state. */
  HS_COLD_MEMBER FullConfigSnapshot capture_full_config_snapshot() const {
    FullConfigSnapshot snapshot;
    snapshot.accepted = encode_config_values(accepted_config);
    snapshot.requested = encode_config_values(requested_config);
    for (size_t index = 0; index < pending_edit_count; ++index)
      snapshot.pending[static_cast<size_t>(pending_edits[index].id)] = 1;
    snapshot.has_runtime = true;
    const ClockState &clocks = runtime.clocks;
    snapshot.runtime[static_cast<size_t>(RuntimeFieldId::SOURCE_PRIMARY)] =
        clocks.source_primary;
    snapshot.runtime[static_cast<size_t>(RuntimeFieldId::SOURCE_SECONDARY)] =
        clocks.source_secondary;
    snapshot.runtime[static_cast<size_t>(RuntimeFieldId::SOURCE_ANGLE)] =
        clocks.source_angle;
    snapshot.runtime[static_cast<size_t>(RuntimeFieldId::WARP_OUTER_ROTATION)] =
        clocks.warp_outer_rotation;
    snapshot.runtime[static_cast<size_t>(RuntimeFieldId::PROJECTION_SPIN)] =
        clocks.projection_spin;
    snapshot.runtime[static_cast<size_t>(RuntimeFieldId::HUE_NOISE_PHASE)] =
        clocks.hue_noise_phase;
    snapshot.runtime[static_cast<size_t>(RuntimeFieldId::SOURCE_NOISE_PHASE)] =
        clocks.source_noise_time;
    snapshot.runtime[static_cast<size_t>(RuntimeFieldId::WARP_INNER_ROTATION)] =
        clocks.warp_inner_rotation;
    snapshot.runtime[static_cast<size_t>(RuntimeFieldId::SURFACE_NOISE_PHASE)] =
        clocks.surface_noise_time;
    snapshot.runtime[static_cast<size_t>(RuntimeFieldId::WARP_OUTER_PHASE)] =
        clocks.warp_outer_phase;
    snapshot.runtime[static_cast<size_t>(RuntimeFieldId::WARP_INNER_PHASE)] =
        clocks.warp_inner_phase;
    snapshot.runtime[static_cast<size_t>(
        RuntimeFieldId::PALETTE_OSCILLATION_PHASE)] =
        clocks.palette_oscillation_phase;
    return snapshot;
  }

  /**
   * @brief Atomically restores a versioned ShaderWorkbench configuration snapshot.
   * @return APPLIED on success; failures leave the effect unchanged.
   */
  HS_COLD_MEMBER ConfigRestoreResult
  restore_full_config_snapshot(const FullConfigSnapshot &snapshot) {
    if (!config_version_supported(snapshot.schema_version))
      return ConfigRestoreResult::UNSUPPORTED_VERSION;

    Config next_accepted{};
    Config next_requested{};
    if (!decode_config_values(snapshot.accepted, next_accepted) ||
        !decode_config_values(snapshot.requested, next_requested))
      return ConfigRestoreResult::INVALID_VALUE;
    normalize_config_ranges(next_accepted);
    normalize_config_ranges(next_requested);
    RuntimeValues next_runtime = snapshot.runtime;
    if (!valid_snapshot_config(next_accepted) ||
        !valid_snapshot_config(next_requested))
      return ConfigRestoreResult::INVALID_VALUE;
    if (!admissible_config(next_accepted))
      return ConfigRestoreResult::INVALID_ACCEPTED;
    if (fixed_topology) {
      const InversePipelineId pipeline = preset_for_view(0).pipeline;
      if (resolve_pipeline_id(next_accepted) != pipeline)
        return ConfigRestoreResult::INVALID_ACCEPTED;
      if (resolve_pipeline_id(next_requested) != pipeline)
        return ConfigRestoreResult::INVALID_PENDING;
    }

    const ConfigValues migrated_accepted = encode_config_values(next_accepted);
    const ConfigValues migrated_requested =
        encode_config_values(next_requested);
    size_t next_pending_count = 0;
    for (size_t index = 0; index < CONFIG_FIELD_COUNT; ++index) {
      if (snapshot.pending[index] > 1)
        return ConfigRestoreResult::INVALID_PENDING;
      const bool differs =
          migrated_accepted[index] != migrated_requested[index];
      if ((snapshot.pending[index] != 0) != differs)
        return ConfigRestoreResult::INVALID_PENDING;
      next_pending_count += differs;
    }
    if (next_pending_count > pending_edits.size())
      return ConfigRestoreResult::INVALID_PENDING;

    if (snapshot.has_runtime)
      for (float value : snapshot.runtime)
        if (!std::isfinite(value))
          return ConfigRestoreResult::INVALID_VALUE;
    if (!prepare_resource_union(next_accepted, next_accepted))
      return ConfigRestoreResult::INVALID_ACCEPTED;

    state->param_morph.active = false;
    state->transition.active = false;
    accepted_config = next_accepted;
    requested_config = next_requested;
    published_config = next_accepted;
    active_slots = next_accepted.slots;
    active_pipeline = resolve_pipeline_id(next_accepted);
    blend.params = next_accepted.params;
    blend.palette_mapping =
        palette_mapping_weights(next_accepted.slots.palette_mapping);
    pending_edit_count = 0;
    for (size_t index = 0; index < CONFIG_FIELD_COUNT; ++index) {
      if (migrated_accepted[index] == migrated_requested[index])
        continue;
      const ConfigFieldId id = static_cast<ConfigFieldId>(index);
      const ConfigFieldLayout layout = config_field_layout(id);
      pending_edits[pending_edit_count++] = {nullptr, id, layout.offset,
                                             layout.size};
    }
    display_config = next_requested;
    if (snapshot.has_runtime) {
      ClockState &clocks = runtime.clocks;
      clocks.source_primary =
          next_runtime[static_cast<size_t>(RuntimeFieldId::SOURCE_PRIMARY)];
      clocks.source_secondary =
          next_runtime[static_cast<size_t>(RuntimeFieldId::SOURCE_SECONDARY)];
      clocks.source_angle =
          next_runtime[static_cast<size_t>(RuntimeFieldId::SOURCE_ANGLE)];
      clocks.warp_outer_rotation = next_runtime[static_cast<size_t>(
          RuntimeFieldId::WARP_OUTER_ROTATION)];
      clocks.projection_spin =
          next_runtime[static_cast<size_t>(RuntimeFieldId::PROJECTION_SPIN)];
      clocks.hue_noise_phase =
          next_runtime[static_cast<size_t>(RuntimeFieldId::HUE_NOISE_PHASE)];
      clocks.source_noise_time =
          next_runtime[static_cast<size_t>(RuntimeFieldId::SOURCE_NOISE_PHASE)];
      clocks.warp_inner_rotation = next_runtime[static_cast<size_t>(
          RuntimeFieldId::WARP_INNER_ROTATION)];
      clocks.surface_noise_time = next_runtime[static_cast<size_t>(
          RuntimeFieldId::SURFACE_NOISE_PHASE)];
      clocks.warp_outer_phase =
          next_runtime[static_cast<size_t>(RuntimeFieldId::WARP_OUTER_PHASE)];
      clocks.warp_inner_phase =
          next_runtime[static_cast<size_t>(RuntimeFieldId::WARP_INNER_PHASE)];
      clocks.palette_oscillation_phase = next_runtime[static_cast<size_t>(
          RuntimeFieldId::PALETTE_OSCILLATION_PHASE)];
    }
    rebind_parameters();
    return ConfigRestoreResult::APPLIED;
  }
#endif

private:
  struct WalkDeltas {
    Quaternion projection;
    Quaternion outer;
  };

  struct ParamMorphRuntime {
    Params from;
    Params to;
    PaletteMappingWeights mapping_from;
    PaletteMappingWeights mapping_to;
    PaletteMapping mapping_destination = PaletteMapping::LINEAR;
    uint16_t elapsed = 0;
    uint16_t duration = 0;
    bool staggered = false;
    bool continue_choreo = false;
    bool active = false;
  };

  struct TransitionRuntime {
    Config from_config;
    Config to_config;
    EndpointRuntime from_runtime;
    EndpointRuntime to_runtime;
    uint16_t elapsed = 0;
    uint16_t duration = 0;
    bool continue_choreo = false;
    bool active = false;
    InversePipelineId from_pipeline = InversePipelineId::NONE;
    InversePipelineId to_pipeline = InversePipelineId::NONE;
  };

  struct StateBundle {
    FrameState frame;
    alignas(std::max_align_t) std::byte prepared_blob[PREPARED_BLOB_BYTES];
    Config render_config;
    std::array<FastNoiseLite, MAX_NOISE_RESOURCES> noise_resources;
    std::array<NoiseFieldKey, MAX_NOISE_RESOURCES> prepared_noise_keys{};
    std::array<Pixel, PreparedHueRotation::LUT_SIZE> hue_rotation_lut;
    std::array<int8_t, PreparedHueNoise::LUT_SIZE> hue_noise_lut;
    Pullback::Color::HueNoiseBakeCache hue_noise_bake;
    FastNoiseLite projection_walk_noise;
    FastNoiseLite outer_walk_noise;
    ParamMorphRuntime param_morph;
    TransitionRuntime transition;

    HS_COLD_MEMBER StateBundle() = default;
  };

  struct ThroughClearPhase {
    float alpha;
    bool from_endpoint;
    bool clear;
  };

  HS_COLD_MEMBER bool prepare_resource_union(const Config &from,
                                             const Config &to) {
    std::array<NoiseFieldKey, MAX_NOISE_RESOURCES> keys{};
    size_t count = 0;
    if (!append_config_resource_keys(from, keys, count) ||
        !append_config_resource_keys(to, keys, count))
      return false;
    prepared_noise_count = count;
    for (size_t index = 0; index < count; ++index) {
      state->prepared_noise_keys[index] = keys[index];
      state->noise_resources[index].SetNoiseType(
          FastNoiseLite::NoiseType_OpenSimplex2);
      state->noise_resources[index].SetSeed(keys[index].seed);
      state->noise_resources[index].SetFrequency(
          keys[index].generator_frequency);
    }
    return true;
  }

  HS_COLD_MEMBER const FastNoiseLite *
  resolve_resource(const NoiseFieldKey &key) const {
    for (size_t index = 0; index < prepared_noise_count; ++index)
      if (state->prepared_noise_keys[index] == key)
        return &state->noise_resources[index];
    return nullptr;
  }

  HS_COLD_MEMBER const FastNoiseLite *
  resolve_warp_resource(const WarpStageSpec &spec) const {
    return warp_uses_noise(spec.kind)
               ? resolve_resource(warp_resource_key(spec))
               : nullptr;
  }

  HS_COLD_MEMBER const FastNoiseLite *
  resolve_source_resource(const Config &config) const {
    return is_noise_contour(config.slots.function)
               ? resolve_resource(source_resource_key(config))
               : nullptr;
  }

  HS_COLD_MEMBER const FastNoiseLite *
  resolve_surface_noise_resource(const Config &config) const {
    return config.slots.surface_noise != SurfaceNoise::NONE
               ? resolve_resource(surface_noise_resource_key(config))
               : nullptr;
  }

  HS_COLD_MEMBER const FastNoiseLite *
  resolve_color_noise_resource(const Config &config) const {
    if (config.slots.hue_shift != HueShiftMode::NOISE ||
        config.params.color.hue_shift_amount == 0.0f)
      return nullptr;
    return resolve_resource(Workbench::color_noise_resource_key());
  }

  HS_COLD_MEMBER const BakedPalette &palette_for(PaletteMode mode) const {
    return generated_palettes.palette(mode);
  }

  PaletteMode visible_palette_mode() const {
    if (!state->transition.active)
      return active_slots.palette;
    const ThroughClearPhase phase = through_clear_phase(
        state->transition.elapsed, state->transition.duration);
    return phase.from_endpoint ? state->transition.from_config.slots.palette
                               : state->transition.to_config.slots.palette;
  }

  float visible_palette_chroma() const {
    if (!state->transition.active)
      return blend.params.color.palette_chroma;
    const ThroughClearPhase phase = through_clear_phase(
        state->transition.elapsed, state->transition.duration);
    return phase.from_endpoint
               ? state->transition.from_config.params.color.palette_chroma
               : state->transition.to_config.params.color.palette_chroma;
  }

  void step_generated_palettes(PaletteMode visible) {
    generated_palettes.step(visible);
  }

  HS_COLD_MEMBER void update_palette_chroma(float chroma) {
    generated_palettes.set_chroma(chroma);
  }

  HS_FLASH_MEMBER static SourceState
  prepare_source_state(const ClockState &clocks) {
    return {clocks.source_primary, clocks.source_secondary, clocks.source_angle,
            fast_cosf(clocks.source_angle), fast_sinf(clocks.source_angle)};
  }

  HS_FLASH_MEMBER static Pullback::Source::PreparedSphericalRings
  prepare_spherical_rings(const EndpointRuntime &endpoint) {
    const Quaternion orientation =
        make_rotation(X_AXIS, endpoint.clocks.source_angle) *
        endpoint.source_wander;
    return {rotate(Y_AXIS, orientation), endpoint.clocks.source_primary};
  }

  HS_FLASH_MEMBER static PreparedSurfaceNoise
  prepare_surface_noise(const ClockState &clocks, const Params &params) {
    const float surface_phase = TWO_PI_F * wrap_t(clocks.surface_noise_time);
    const float surface_direction = TWO_PI_F * params.surface_noise.direction;
    return {Vector(NOISE_LOOP_RADIUS * cosf(surface_phase),
                   NOISE_LOOP_RADIUS * sinf(surface_phase), 0.0f),
            cosf(surface_direction), sinf(surface_direction)};
  }

  HS_FLASH_MEMBER static PreparedWarpStage
  prepare_warp_stage(const WarpStageSpec &spec, const WarpStageParams &params,
                     float stage_phase,
                     const Complex &source_period = Complex(),
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
    } else if (spec.kind == WarpStageKind::VECTOR_NOISE) {
      const float angle = TWO_PI_F * wrap_t(stage_phase);
      prepared.transform.noise_loop = {
          NOISE_LOOP_RADIUS * sinf(angle) * 0.7071067811865475f,
          NOISE_LOOP_RADIUS * cosf(angle),
      };
    }
    prepared.rotation_cos = cosf(rotation);
    prepared.rotation_sin = sinf(rotation);
    return prepared;
  }

  HS_COLD_MEMBER FrameState prepare_frame() const {
    FrameState frame;
    prepare_frame({active_slots, blend.params}, runtime, frame);
    return frame;
  }

  HS_COLD_MEMBER FrameState
  prepare_frame(const Config &config, const EndpointRuntime &endpoint) const {
    FrameState frame;
    prepare_frame(config, endpoint, frame);
    return frame;
  }

  HS_COLD_MEMBER void prepare_frame(const Config &config,
                                    const EndpointRuntime &endpoint,
                                    FrameState &frame) const {
    const bool animated_projection =
        config.slots.projection_frame == ProjectionFramePolicy::SPIN_WANDER;
    const BakedPalette *palette = &palette_for(config.slots.palette);
    PreparedHueRotation prepared_hue_rotation{
        state->hue_rotation_lut.data(),
        config.slots.hue_shift != HueShiftMode::NONE &&
            config.params.color.hue_shift_amount != 0.0f};
    if (prepared_hue_rotation.active)
      prepare_hue_rotation_lut(prepared_hue_rotation, *palette);
    const FastNoiseLite *color_noise = resolve_color_noise_resource(config);
    PreparedHueNoise prepared_hue_noise{
        state->hue_noise_lut.data(),
        config.slots.hue_shift == HueShiftMode::NOISE &&
            config.params.color.hue_shift_amount != 0.0f};
    if (prepared_hue_noise.active && color_noise != nullptr)
      state->hue_noise_bake.refresh(state->hue_noise_lut, *color_noise,
                                    config.params.color.hue_noise_scale,
                                    endpoint.clocks.hue_noise_phase);
    frame.slots = config.slots;
    frame.params = config.params;
    frame.palette_mapping =
        state->param_morph.active
            ? blend.palette_mapping
            : palette_mapping_weights(config.slots.palette_mapping);
    frame.clocks = endpoint.clocks;
    frame.transforms = {animated_projection
                            ? endpoint.transforms.projection_conj
                            : Quaternion(),
                        endpoint.transforms.outer_conj};
    frame.meridian_cos = cosf(config.params.projection.central_meridian);
    frame.meridian_sin = sinf(config.params.projection.central_meridian);
    frame.dynamic = {
        prepare_source_state(endpoint.clocks),
        prepare_spherical_rings(endpoint),
        {prepare_warp_stage(
             config.slots.warp_program.outer, config.params.warp.outer,
             endpoint.clocks.warp_outer_phase, source_cartesian_period(config),
             endpoint.clocks.warp_outer_rotation),
         prepare_warp_stage(
             config.slots.warp_program.inner, config.params.warp.inner,
             endpoint.clocks.warp_inner_phase, source_cartesian_period(config),
             endpoint.clocks.warp_inner_rotation)},
        prepare_surface_noise(endpoint.clocks, config.params)};
    frame.prepared_hue_rotation = prepared_hue_rotation;
    frame.prepared_hue_noise = prepared_hue_noise;
    frame.resources = {resolve_warp_resource(config.slots.warp_program.outer),
                       resolve_warp_resource(config.slots.warp_program.inner),
                       resolve_source_resource(config),
                       resolve_surface_noise_resource(config),
                       color_noise,
                       palette};
  }

  static ThroughClearPhase through_clear_phase(uint16_t elapsed,
                                               uint16_t duration) {
    const uint16_t center = duration / 2;
    if (elapsed == center)
      return {0.0f, false, true};
    const bool from_endpoint = elapsed < center;
    const float phase = from_endpoint ? static_cast<float>(elapsed) / center
                                      : static_cast<float>(elapsed - center) /
                                            (duration - center);
    return {from_endpoint ? 1.0f - ease_in_out_sin(phase)
                          : ease_in_out_sin(phase),
            from_endpoint, false};
  }

  HS_FLASH_MEMBER void draw_through_clear_transition(Canvas &canvas) {
    const ThroughClearPhase phase = through_clear_phase(
        state->transition.elapsed, state->transition.duration);
    if (phase.clear)
      return;
    PreparedEndpoint prepared;
    const Config &config = phase.from_endpoint ? state->transition.from_config
                                               : state->transition.to_config;
    const EndpointRuntime &endpoint = phase.from_endpoint
                                          ? state->transition.from_runtime
                                          : state->transition.to_runtime;
    const InversePipelineId pipeline = phase.from_endpoint
                                           ? state->transition.from_pipeline
                                           : state->transition.to_pipeline;
    HS_CHECK(
        prepare_endpoint(config, endpoint, phase.alpha, pipeline, prepared),
        "ShaderWorkbench transition endpoint has no renderer");
    draw_endpoint(canvas, prepared,
                  phase.from_endpoint ? ProfileEndpoint::FROM
                                      : ProfileEndpoint::TO);
  }

  HS_COLD_MEMBER bool prepare_endpoint(const Config &config,
                                       const EndpointRuntime &endpoint,
                                       float alpha, InversePipelineId selected,
                                       PreparedEndpoint &prepared) const {
    const ProgramDescriptor *program = get_inverse_program(selected);
    ShadeFunction shade;
    bool (*resources_ready)(const FrameState &);
    if (program != nullptr) {
      if (program->key != make_topology_key(config) ||
          !program->continuous_parameters_supported(config))
        return false;
      shade = program->shade;
      resources_ready = program->resources_ready;
    } else {
#if HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND
      if (selected != InversePipelineId::NONE || !valid_config(config))
        return false;
      shade = &shade_dynamic;
      resources_ready = &pipeline_resources_ready;
#else
      return false;
#endif
    }
    prepared.frame = &state->frame;
    prepare_frame(config, endpoint, *prepared.frame);
    if (!resources_ready(*prepared.frame))
      return false;
    if (program != nullptr) {
      program->prepare(*prepared.frame, state->prepared_blob);
      prepared.prepared = state->prepared_blob;
    } else {
      prepared.prepared = nullptr;
    }
    prepared.shade = shade;
    prepared.pipeline = selected;
    prepared.alpha = alpha;
#if defined(HS_PROFILE_ENABLE)
    prepared.preset = selected_preset_index(config, selected);
#endif
    return true;
  }

  HS_FLASH_MEMBER void
  draw_endpoint(Canvas &canvas, PreparedEndpoint &prepared,
                ProfileEndpoint endpoint = ProfileEndpoint::STEADY) {
#if defined(HS_PROFILE_ENABLE)
    emit_pullback_program(prepared, endpoint);
#else
    (void)endpoint;
#endif
    FrameShader shader{prepared.frame, prepared.alpha, prepared.shade,
                       prepared.prepared};
    HS_PROFILE(sb_shader_draw);
    Scan::Shader::draw<W, H, 1>(canvas, shader);
  }

  static constexpr const char *pipeline_name(InversePipelineId pipeline) {
    switch (pipeline) {
    case InversePipelineId::GLITCH_NOISE_GRID_WAVE_SHEAR:
      return "GLITCH_NOISE_GRID_WAVE_SHEAR";
    case InversePipelineId::KALEIDOSCOPE_TWIN_WAVE_INNER_MIRROR:
      return "KALEIDOSCOPE_TWIN_WAVE_INNER_MIRROR";
    case InversePipelineId::GNOMONIC_KALEIDOSCOPE_GRID_MIRROR:
      return "GNOMONIC_KALEIDOSCOPE_GRID_MIRROR";
    case InversePipelineId::GNOMONIC_ALIEN_CORE_MIRROR:
      return "GNOMONIC_ALIEN_CORE_MIRROR";
    case InversePipelineId::PEIRCE_DODECAHEDRAL_GRID:
      return "PEIRCE_DODECAHEDRAL_GRID";
    case InversePipelineId::GNOMONIC_DODECAHEDRAL_GRID_WAVE_MIRROR:
      return "GNOMONIC_DODECAHEDRAL_GRID_WAVE_MIRROR";
    case InversePipelineId::GNOMONIC_AFFINE_LATTICE_CONTOUR:
      return "GNOMONIC_AFFINE_LATTICE_CONTOUR";
    case InversePipelineId::SINUSOIDAL_LATTICE_MELT:
      return "SINUSOIDAL_LATTICE_MELT";
    case InversePipelineId::STEREOGRAPHIC_PRISM_POLAR_WAVE_LATTICE:
      return "STEREOGRAPHIC_PRISM_POLAR_WAVE_LATTICE";
    case InversePipelineId::GNOMONIC_DODECAHEDRAL_GRID_VECTOR_MIRROR:
      return "GNOMONIC_DODECAHEDRAL_GRID_VECTOR_MIRROR";
    case InversePipelineId::STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR:
      return "STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR";
    case InversePipelineId::
        STEREOGRAPHIC_HEXAGONAL_PRISM_TWIN_WAVE_INNER_MIRROR:
      return "STEREOGRAPHIC_HEXAGONAL_PRISM_TWIN_WAVE_INNER_MIRROR";
    case InversePipelineId::EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR:
      return "EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR";
    case InversePipelineId::STEREOGRAPHIC_ALIEN_CORE_MIRROR:
      return "STEREOGRAPHIC_ALIEN_CORE_MIRROR";
    case InversePipelineId::STEREOGRAPHIC_MOBIUS_TWIN_WAVE_INNER_MIRROR:
      return "STEREOGRAPHIC_MOBIUS_TWIN_WAVE_INNER_MIRROR";
    case InversePipelineId::COUNT:
      return "COUNT";
    case InversePipelineId::NONE:
      return "NONE";
    }
    return "NONE";
  }

#if defined(HS_PROFILE_ENABLE)
  static constexpr const char *profile_endpoint_name(ProfileEndpoint endpoint) {
    switch (endpoint) {
    case ProfileEndpoint::STEADY:
      return "steady";
    case ProfileEndpoint::FROM:
      return "from";
    case ProfileEndpoint::TO:
      return "to";
    }
    return "steady";
  }

  size_t selected_preset_index(const Config &config,
                               InversePipelineId pipeline) const {
    for (size_t index = 0; index < preset_count_for_view(); ++index)
      if (preset_for_view(index).pipeline == pipeline &&
          preset_for_view(index).config == config)
        return index;
    return getPresetIndex();
  }

  void emit_pullback_program(const PreparedEndpoint &prepared,
                             ProfileEndpoint endpoint) {
    if (profile_program_valid && profile_program_preset == prepared.preset &&
        profile_program_pipeline == prepared.pipeline &&
        profile_program_endpoint == endpoint)
      return;
    hs::log("Pullback program: preset=%u/%u pipeline=%s endpoint=%s",
            static_cast<unsigned>(prepared.preset),
            static_cast<unsigned>(preset_count_for_view()),
            pipeline_name(prepared.pipeline), profile_endpoint_name(endpoint));
    profile_program_valid = true;
    profile_program_preset = prepared.preset;
    profile_program_pipeline = prepared.pipeline;
    profile_program_endpoint = endpoint;
  }
#endif

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
  static Color4 shade_dynamic(const Vector &view, const FrameState &frame,
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
  HS_FLASH_MEMBER static Color4
  shade_projected(const ProjectedLookup &projected, const FrameState &frame) {
    HS_SB_STAGE_MARK(stage_start);
    const PlanarWarpResult warped = planar_warp_lookup(projected, frame);
    HS_SB_STAGE_SPAN(planar_warp, stage_start);
    const Complex source_coords = condition_source_coords(warped.coords, frame);
    const float field = sample_source(source_coords, projected, frame);
    HS_SB_STAGE_SPAN(source, stage_start);
    const FieldSample material =
        shape_material(field, projected, warped, frame);
    HS_SB_STAGE_SPAN(material, stage_start);
    const Color4 color = colorize(material, frame);
    HS_SB_STAGE_SPAN(color, stage_start);
    return color;
  }
#endif

  __attribute__((always_inline)) static Pullback::ProjectionResult
  stereographic_lookup(const Vector &local, const FrameState &frame) {
    return Pullback::Projection::stereographic(
        local, frame.params.projection.singularity_fade);
  }

  /** @brief Whether this frame's colorizer reads displacement metadata. */
  __attribute__((always_inline)) static bool
  tracks_displacement(const FrameState &frame) {
    return frame.prepared_hue_rotation.active &&
           frame.slots.hue_shift == HueShiftMode::WARP_DISPLACEMENT;
  }

  /**
   * @brief Builds one manifest row.
   * @tparam Pipeline Compiled inverse pipeline the row selects.
   * @tparam Id Stable identifier for the row.
   * @tparam Key Topology the row is matched on.
   * @param continuous Predicate over the continuous parameters @p Pipeline
   *        serves.
   * @details Rejects at compile time a pipeline whose stages hardcode a
   * topology facet that @p Key does not carry.
   */
  template <typename Pipeline, InversePipelineId Id, TopologyKey Key>
  static constexpr ProgramDescriptor
  make_program(bool (*continuous)(const Config &)) {
    static_assert(Pipeline::implements(Key),
                  "inverse pipeline does not implement its topology key");
    static_assert(sizeof(typename Pipeline::PreparedTuple) <=
                      PREPARED_BLOB_BYTES,
                  "prepared blob capacity exceeded");
    return {Id,
            Key,
            &Pipeline::shade_prepared,
            &Pipeline::prepare_into,
            continuous,
            &pipeline_resources_ready};
  }

  HS_COLD_MEMBER static const std::array<ProgramDescriptor,
                                         Workbench::INVERSE_PROGRAM_COUNT> &
  inverse_programs() {
    static constexpr std::array<ProgramDescriptor,
                                Workbench::INVERSE_PROGRAM_COUNT>
        PROGRAMS{{
            make_program<GlitchNoiseGridWaveShearPipeline,
                         InversePipelineId::GLITCH_NOISE_GRID_WAVE_SHEAR,
                         make_topology_key(
                             Workbench::wave_shear_generated_preset())>(
                &all_continuous_parameters_supported),
            make_program<KaleidoscopeTwinWaveInnerMirrorPipeline,
                         InversePipelineId::KALEIDOSCOPE_TWIN_WAVE_INNER_MIRROR,
                         make_topology_key(
                             Workbench::kaleidoscope_mirror_preset())>(
                &all_continuous_parameters_supported),
            make_program<
                GnomonicKaleidoscopeGridMirrorPipeline,
                InversePipelineId::GNOMONIC_KALEIDOSCOPE_GRID_MIRROR,
                make_topology_key(
                    Workbench::gnomonic_kaleidoscope_grid_mirror_preset())>(
                &all_continuous_parameters_supported),
            make_program<GnomonicAlienCoreMirrorPipeline,
                         InversePipelineId::GNOMONIC_ALIEN_CORE_MIRROR,
                         make_topology_key(
                             Workbench::gnomonic_grid_mirror_preset(
                                 SurfaceLens::GLITCH))>(
                &all_continuous_parameters_supported),
            make_program<
                PeirceDodecahedralGridPipeline,
                InversePipelineId::PEIRCE_DODECAHEDRAL_GRID,
                make_topology_key(
                    Workbench::peirce_dodecahedral_generated_preset())>(
                &all_continuous_parameters_supported),
            make_program<
                GnomonicDodecahedralGridWaveMirrorPipeline,
                InversePipelineId::GNOMONIC_DODECAHEDRAL_GRID_WAVE_MIRROR,
                make_topology_key(
                    Workbench::gnomonic_wave_shear_grid_preset())>(
                &all_continuous_parameters_supported),
            make_program<
                GnomonicAffineLatticeContourPipeline,
                InversePipelineId::GNOMONIC_AFFINE_LATTICE_CONTOUR,
                make_topology_key(
                    Workbench::gnomonic_affine_lattice_contour_preset())>(
                &all_continuous_parameters_supported),
            make_program<SinusoidalLatticeMeltPipeline,
                         InversePipelineId::SINUSOIDAL_LATTICE_MELT,
                         make_topology_key(
                             Workbench::sinusoidal_lattice_curl_preset(1.0f))>(
                &all_continuous_parameters_supported),
            make_program<
                StereographicPrismPolarWaveLatticePipeline,
                InversePipelineId::STEREOGRAPHIC_PRISM_POLAR_WAVE_LATTICE,
                make_topology_key(
                    Workbench::
                        stereographic_prism_polar_wave_lattice_preset())>(
                &all_continuous_parameters_supported),
            make_program<
                GnomonicDodecahedralGridVectorMirrorPipeline,
                InversePipelineId::GNOMONIC_DODECAHEDRAL_GRID_VECTOR_MIRROR,
                make_topology_key(
                    Workbench::
                        gnomonic_dodecahedral_vector_mirror_grid_preset())>(
                &all_continuous_parameters_supported),
            make_program<
                StereographicDodecahedralGridInnerMirrorPipeline,
                InversePipelineId::STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR,
                make_topology_key(
                    Workbench::
                        stereographic_dodecahedral_grid_inner_mirror_preset())>(
                &all_continuous_parameters_supported),
            make_program<
                StereographicHexagonalPrismTwinWaveInnerMirrorPipeline,
                InversePipelineId::
                    STEREOGRAPHIC_HEXAGONAL_PRISM_TWIN_WAVE_INNER_MIRROR,
                make_topology_key(
                    Workbench::
                        stereographic_hexagonal_prism_twin_wave_mirror_preset())>(
                &all_continuous_parameters_supported),
            make_program<
                EquirectangularDodecahedralGridInnerMirrorPipeline,
                InversePipelineId::
                    EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR,
                make_topology_key(
                    Workbench::
                        equirectangular_dodecahedral_double_mapping_grid_inner_mirror_preset())>(
                &all_continuous_parameters_supported),
            make_program<
                StereographicAlienCoreMirrorPipeline,
                InversePipelineId::STEREOGRAPHIC_ALIEN_CORE_MIRROR,
                make_topology_key(
                    Workbench::stereographic_alien_core_mirror_preset())>(
                &all_continuous_parameters_supported),
            make_program<
                StereographicMobiusTwinWaveInnerMirrorPipeline,
                InversePipelineId::STEREOGRAPHIC_MOBIUS_TWIN_WAVE_INNER_MIRROR,
                make_topology_key(
                    Workbench::
                        stereographic_mobius_twin_wave_inner_mirror_preset())>(
                &all_continuous_parameters_supported),
        }};
    return PROGRAMS;
  }

  HS_COLD_MEMBER static const ProgramDescriptor *
  find_inverse_program(const Config &config) {
    const TopologyKey key = make_topology_key(config);
    for (const ProgramDescriptor &program : inverse_programs())
      if (program.key == key && program.continuous_parameters_supported(config))
        return &program;
    return nullptr;
  }

  HS_COLD_MEMBER static const ProgramDescriptor *
  get_inverse_program(InversePipelineId id) {
    for (const ProgramDescriptor &program : inverse_programs())
      if (program.id == id)
        return &program;
    return nullptr;
  }

  HS_COLD_MEMBER static InversePipelineId
  resolve_pipeline_id(const Config &config) {
    const ProgramDescriptor *program = find_inverse_program(config);
    return program == nullptr ? InversePipelineId::NONE : program->id;
  }

  HS_COLD_MEMBER static const ProgramDescriptor *
  resolve_inverse_program(const FrameState &frame) {
    const Config config{frame.slots, frame.params};
    const ProgramDescriptor *program = find_inverse_program(config);
    if (program == nullptr || !program->resources_ready(frame))
      return nullptr;
    return program;
  }

  static bool projection_edge_distance_required(const FrameState &frame) {
    const WarpProgram &program = frame.slots.warp_program;
    return frame.slots.coverage == CoveragePolicy::EDGE_FADE ||
           (program.outer.kind != WarpStageKind::NONE &&
            program.outer.envelope == WarpEnvelope::EDGE_FADE) ||
           (program.inner.kind != WarpStageKind::NONE &&
            program.inner.envelope == WarpEnvelope::EDGE_FADE);
  }

  __attribute__((always_inline)) static Vector
  outer_camera_lookup(const Vector &view, const FrameState &frame) {
    return rotate(view, frame.transforms.outer_conj);
  }

#if HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND ||                              \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
  static ProjectedLookup surface_lens_project_lookup(const Vector &v,
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

  HS_FLASH_MEMBER static Pullback::ProjectionResult
  project_bonne(const Vector &local, const FrameState &frame) {
    return Pullback::Projection::bonne(
        local, frame.params.projection.central_meridian,
        (frame.slots.bonne_hemisphere == BonneHemisphere::NORTH ? 1.0f
                                                                : -1.0f) *
            frame.params.projection.bonne_standard_parallel,
        frame.params.projection.coordinate_scale);
  }

  HS_FLASH_MEMBER static Pullback::ProjectionResult
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

  HS_FLASH_MEMBER static Pullback::ProjectionResult
  project_airocean(const Vector &local, const FrameState &frame) {
    return Pullback::Projection::airocean(
        local, frame.params.projection.central_meridian,
        frame.slots.airocean_layout == AiroceanLayout::HORIZONTAL,
        projection_edge_distance_required(frame),
        frame.params.projection.coordinate_scale);
  }

  HS_FLASH_MEMBER static Pullback::ProjectionResult
  project_sinusoidal(const Vector &local, const FrameState &frame) {
    return Pullback::Projection::folded_sinusoidal(
        local, frame.params.projection.central_meridian);
  }

  HS_FLASH_MEMBER static Pullback::ProjectionResult
  project_equirectangular(const Vector &local, const FrameState &frame) {
    return Pullback::Projection::equirectangular(
        local, frame.params.projection.central_meridian,
        frame.params.projection.singularity_fade);
  }

  static constexpr Pullback::Projection::GnomonicHemisphere
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

  HS_FLASH_MEMBER static Pullback::ProjectionResult
  project_gnomonic(const Vector &local, const FrameState &frame) {
    return Pullback::Projection::gnomonic(
        local, frame.params.projection.singularity_fade,
        pullback_gnomonic_hemisphere(frame.slots.gnomonic_hemisphere));
  }

  HS_FLASH_MEMBER static Pullback::ProjectionResult
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

  HS_FLASH_MEMBER static ProjectedLookup
  project_branch(const Vector &v, const FrameState &frame) {
    const Vector local = rotate(v, frame.transforms.projection_conj);
    const Pullback::ProjectionResult result =
        frame.slots.projection != Projection::STEREOGRAPHIC
            ? project_nonstereographic(local, frame)
            : stereographic_lookup(local, frame);
    return {result.coords, result.provenance, local, 0.0f};
  }

  __attribute__((always_inline)) static ProjectedLookup
  profiled_project_branch(const Vector &v, const FrameState &frame) {
    HS_SB_STAGE_MARK(stage_start);
    const ProjectedLookup projected = project_branch(v, frame);
    HS_SB_STAGE_SPAN(projection, stage_start);
    return projected;
  }

  HS_FLASH_MEMBER static Pullback::ProjectionResult
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
      return {coords,
              {0, 0, BOUNDARY_CUT, PI_F - std::fabs(coords.re),
               Pullback::Projection::equirectangular_weight(local,
                                                            singularity_fade),
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
  HS_FLASH_MEMBER static PlanarWarpResult
  planar_warp_lookup(const ProjectedLookup &projected,
                     const FrameState &frame) {
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

  HS_FLASH_MEMBER static PlanarWarpStageResult
  finish_closed_form_warp(const Complex &input, const Complex &output,
                          bool path_length_required) {
    return Pullback::Warp::finish_closed_form(input, output,
                                              path_length_required);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_affine_frame(const Complex &input, const WarpStageParams &,
                    const PreparedWarpStage &prepared,
                    bool path_length_required) {
    return Pullback::Warp::affine_frame(input, prepared, path_length_required);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_wave_shear(const Complex &input, const WarpStageParams &params,
                  float stage_phase, float amplitude,
                  const PreparedWarpStage &prepared,
                  bool path_length_required) {
    return Pullback::Warp::wave_shear(input, params, stage_phase, amplitude,
                                      prepared, path_length_required);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_vortex(const Complex &input, const PreparedWarpStage &prepared,
              bool path_length_required) {
    return Pullback::Warp::vortex(input, prepared, path_length_required);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_vector_noise(const Complex &input, const WarpStageSpec &spec,
                    const WarpStageParams &params, float amplitude,
                    const FastNoiseLite &noise,
                    const PreparedWarpStage &prepared,
                    bool path_length_required) {
    return Pullback::Warp::vector_noise(input, params, amplitude, noise,
                                        spec.basis, prepared,
                                        path_length_required);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_curl_flow(const Complex &input, const WarpStageSpec &spec,
                 const WarpStageParams &params, float amplitude,
                 const FastNoiseLite &noise, float stage_phase,
                 bool path_length_required) {
    const uint8_t intervals =
        spec.curl_integrator == CurlIntegrator::EULER_1      ? 1
        : spec.curl_integrator == CurlIntegrator::MIDPOINT_2 ? 2
                                                             : 4;
    return Pullback::Warp::curl_flow(input, noise, spec.basis, intervals,
                                     params.scale, amplitude, stage_phase,
                                     path_length_required);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_polar_chart(const Complex &input, const WarpStageSpec &spec,
                   const WarpStageParams &params, float stage_phase) {
    return Pullback::Warp::polar_chart(
        input, params, stage_phase, spec.polar_mode == PolarMode::LOGARITHMIC,
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
  HS_FLASH_MEMBER static PlanarWarpStageResult warp_stage_lookup(
      const Complex &input, const Pullback::ProjectionProvenance &provenance,
      const WarpStageSpec &spec, const WarpStageParams &params,
      float stage_phase, const FastNoiseLite *stage_noise,
      const PreparedWarpStage &prepared, bool path_length_required) {
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

  static float warp_envelope(const Pullback::ProjectionProvenance &provenance,
                             WarpEnvelope envelope, float edge_width) {
    return Pullback::Warp::envelope(provenance, edge_width,
                                    envelope == WarpEnvelope::PROJECTION_WEIGHT,
                                    envelope == WarpEnvelope::EDGE_FADE);
  }

  HS_FLASH_MEMBER static Complex
  curl_flow(const Complex &input, const FastNoiseLite &noise,
            const WarpStageSpec &spec, const WarpStageParams &params,
            float distance, float phase, float &path_length) {
    const uint8_t intervals =
        spec.curl_integrator == CurlIntegrator::EULER_1      ? 1
        : spec.curl_integrator == CurlIntegrator::MIDPOINT_2 ? 2
                                                             : 4;
    constexpr bool path_length_required = true;
    const PlanarWarpStageResult result = Pullback::Warp::curl_flow(
        input, noise, spec.basis, intervals, params.scale, distance, phase,
        path_length_required);
    path_length = result.path_length;
    return result.coords;
  }

  HS_FLASH_MEMBER static Complex curl_vector(const Complex &p,
                                             const FastNoiseLite &noise,
                                             NoiseBasis basis, float scale,
                                             float phase) {
    return Pullback::Warp::curl_vector(p, noise, basis, scale, phase);
  }

  __attribute__((always_inline)) static Complex
  mirror_tile(const Complex &input, const WarpStageParams &params,
              const PreparedWarpStage &prepared) {
    return Pullback::Warp::mirror_tile_coords(input, params, prepared);
  }

#if HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND ||                              \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
  static Complex condition_source_coords(const Complex &coords,
                                         const FrameState &frame) {
    if (is_noise_contour(frame.slots.function) ||
        frame.slots.function == Function::PRIMITIVE_LATTICE ||
        frame.slots.function == Function::FRACTAL ||
        frame.slots.function == Function::TESSELLATION)
      return coords;
    return stereo_pattern_args(coords, frame.params.source.pattern_freq);
  }

  __attribute__((always_inline)) static float
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

  HS_FLASH_MEMBER static FieldSample
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

  static FieldSample shape_material(float field,
                                    const ProjectedLookup &projected,
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
  HS_FLASH_MEMBER static float sample_noise_contour(const Vector &q,
                                                    const FrameState &frame) {
    return Pullback::Source::noise_contour(*frame.resources.source_noise,
                                           frame.params.source.noise_basis, q,
                                           frame.params.source.noise_contrast);
  }

#if HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND ||                              \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
  static Pullback::Source::SphericalRingsSourceParams
  spherical_rings_params(const SourceParams &params) {
    return {static_cast<float>(params.ring_count),
            params.ring_thickness,
            params.ring_softness,
            params.speed,
            params.angle_rate,
            params.ring_wander};
  }

  static Pullback::Source::FractalSourceParams
  fractal_params(const SourceParams &params) {
    return {
        params.fractal_scale,   static_cast<float>(params.fractal_iterations),
        params.julia_mix,       params.julia_real,
        params.julia_imaginary, params.fractal_contours,
        params.speed,           params.angle_rate};
  }

  static Pullback::Source::TessellationSourceParams
  tessellation_params(const SourceParams &params) {
    return {params.tessellation_cell_scale, params.tessellation_line_thickness,
            params.tessellation_line_softness, params.angle_rate};
  }

  HS_FLASH_MEMBER static float sample_source(const Complex &p,
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

  HS_FLASH_MEMBER static float primitive_lattice(const Complex &p,
                                                 const SourceParams &params) {
    return Pullback::Source::primitive_lattice(p, params);
  }

  /**
   * @brief Maps a shaped material sample to a palette colour.
   * @param sample Shaped value, coverage, and surface coordinate.
   * @param frame Frame snapshot.
   * @return Colour whose alpha is the palette alpha scaled by the coverage.
   */
  HS_FLASH_MEMBER static Color4 colorize_generated(const FieldSample &sample,
                                                   const FrameState &frame) {
    return Pullback::Color::GeneratedPalette<ColorStateProvider>::apply(sample,
                                                                        frame);
  }

  __attribute__((always_inline)) static float
  palette_mapping_coordinate(float value, PaletteMapping mapping,
                             float frequency, float offset) {
    return Pullback::Color::palette_mapping_coordinate(
        value, static_cast<Pullback::Color::PaletteMapping>(mapping), frequency,
        offset);
  }

  __attribute__((always_inline)) static float
  brightness_envelope_gain(float value, BrightnessEnvelope envelope,
                           float depth) {
    return Pullback::Color::brightness_envelope_gain(
        value, static_cast<Pullback::Color::BrightnessEnvelope>(envelope),
        depth);
  }

#if HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND ||                              \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
  HS_O3_FN static Color4 colorize(const FieldSample &sample,
                                  const FrameState &frame) {
    return colorize_generated(sample, frame);
  }
#endif

  HS_FLASH_MEMBER static void
  prepare_hue_rotation_lut(PreparedHueRotation &prepared,
                           const BakedPalette &palette) {
    Pullback::Color::prepare_hue_rotation_lut(
        std::span<Pixel, Pullback::Color::HueRotationLutView::SIZE>(
            prepared.lut, PreparedHueRotation::LUT_SIZE),
        palette);
  }

  HS_FLASH_MEMBER static Vector hue_noise_face_direction(int face, float u,
                                                         float v) {
    return Pullback::Color::hue_noise_face_direction(face, u, v);
  }

  HS_FLASH_MEMBER static float
  sample_hue_noise_lut(const PreparedHueNoise &prepared, const Vector &v) {
    return Pullback::Color::sample_hue_noise_lut(
        {prepared.lut, prepared.active}, v);
  }

  __attribute__((always_inline)) static Pixel
  sample_hue_rotation_lut(const PreparedHueRotation &prepared, float value,
                          float amount) {
    return Pullback::Color::sample_hue_rotation_lut(
        {prepared.lut, prepared.active}, value, amount);
  }

  HS_FLASH_MEMBER static Vector apply_lens(const Vector &v,
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

  __attribute__((always_inline)) static Vector
  profiled_apply_lens(const Vector &v, const FrameState &frame) {
    HS_SB_STAGE_MARK(stage_start);
    const Vector lensed = apply_lens(v, frame);
    HS_SB_STAGE_SPAN(lens, stage_start);
    return lensed;
  }

  HS_FLASH_MEMBER static Vector surface_curl_field(const Vector &v,
                                                   const FrameState &frame) {
    return Pullback::Surface::curl_field(
        v, *frame.resources.surface_noise, frame.params.surface_noise.basis,
        frame.params.surface_noise.scale,
        frame.dynamic.surface_noise.loop_offset);
  }

  HS_FLASH_MEMBER static SurfaceNoiseResult
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
    return Pullback::Surface::curl_noise(
        v, *frame.resources.surface_noise, params.basis, integrator,
        params.scale, frame.dynamic.surface_noise.loop_offset, params.strength,
        path_length_required);
  }

  HS_FLASH_MEMBER static Vector apply_surface_noise(const Vector &v,
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
  __attribute__((always_inline)) static Vector
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
      return lenses::polyhedral_kaleidoscope_lens(v,
                                                  lenses::TETRAHEDRAL_MIRRORS);
    case SurfaceLens::KALEIDOSCOPE_OCTAHEDRAL:
      return lenses::polyhedral_kaleidoscope_lens(v,
                                                  lenses::OCTAHEDRAL_MIRRORS);
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

  HS_FLASH_MEMBER static Vector mobius_lens(const Vector &v,
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
  HS_FLASH_MEMBER static Complex project_point(const Vector &v,
                                               Projection projection) {
    switch (projection) {
    case Projection::SINUSOIDAL:
      return projections::folded_sinusoidal(v);
    case Projection::EQUIRECTANGULAR:
      return projections::equirectangular(v);
    case Projection::STEREOGRAPHIC:
      return stereo(v);
    case Projection::GNOMONIC:
      return gnomonic(v);
    case Projection::BONNE:
    case Projection::PEIRCE_QUINCUNCIAL:
    case Projection::AIROCEAN:
      break;
    }
    __builtin_unreachable();
  }

  HS_FLASH_MEMBER static Complex gnomonic(const Vector &v) {
    return Pullback::Projection::gnomonic(
               v, 0.0f, Pullback::Projection::GnomonicHemisphere::FOLDED)
        .coords;
  }

  HS_FLASH_MEMBER static float sample_function(Function function,
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

  HS_FLASH_MEMBER static float twin_wave(const Complex &p,
                                         const SourceState &source) {
    return Pullback::Source::twin_wave(p, source);
  }

  HS_FLASH_MEMBER static float rings(const Complex &p,
                                     const SourceState &source) {
    return Pullback::Source::rings(p, source);
  }

  HS_FLASH_MEMBER static float spiral(const Complex &p,
                                      const SourceState &source) {
    return Pullback::Source::spiral(p, source);
  }

  HS_O3_FN static float grid(const Complex &p, const SourceParams &params,
                             const SourceState &source) {
    return Pullback::Source::grid(p, params, source);
  }

  HS_COLD_MEMBER WalkDeltas sample_walk_deltas() {
#if HS_ENABLE_TEST_HOOKS
    ++walk_step_count;
#endif
    const Quaternion projection = projection_walk.get();
    const Quaternion projection_delta =
        projection * projection_walk_prev.conjugate();
    projection_walk_prev = projection;
    const Quaternion outer = outer_walk.get();
    const Quaternion outer_delta = outer * outer_walk_prev.conjugate();
    outer_walk_prev = outer;
    return {projection_delta.normalized(), outer_delta.normalized()};
  }

  HS_COLD_MEMBER void update_spatial_frames(EndpointRuntime &endpoint,
                                            const Config &config,
                                            const WalkDeltas &deltas) const {
    endpoint.projection_wander = (slerp(Quaternion(), deltas.projection,
                                        config.params.projection.wander) *
                                  endpoint.projection_wander)
                                     .normalized();
    endpoint.outer_wander =
        (slerp(Quaternion(), deltas.outer, config.params.outer_camera.wander) *
         endpoint.outer_wander)
            .normalized();
    endpoint.source_wander =
        (slerp(Quaternion(), deltas.outer, config.params.source.ring_wander) *
         endpoint.source_wander)
            .normalized();
    endpoint.transforms.projection_conj =
        (make_rotation(Y_AXIS, endpoint.clocks.projection_spin) *
         base_orientation * endpoint.projection_wander)
            .conjugate();
    endpoint.transforms.outer_conj = endpoint.outer_wander.conjugate();
  }

  HS_COLD_MEMBER void advance_runtime(EndpointRuntime &endpoint,
                                      const Config &config,
                                      const WalkDeltas &deltas) const {
    const Params &params = config.params;
    endpoint.clocks.source_primary =
        fmodf(endpoint.clocks.source_primary + params.source.speed, TWO_PI_F);
    endpoint.clocks.source_secondary =
        fmodf(endpoint.clocks.source_secondary +
                  params.source.speed * params.source.secondary_rate,
              TWO_PI_F);
    endpoint.clocks.source_angle = fmodf(
        endpoint.clocks.source_angle + params.source.angle_rate, TWO_PI_F);
    endpoint.clocks.projection_spin =
        fmodf(endpoint.clocks.projection_spin + params.projection.spin_rate,
              TWO_PI_F);
    endpoint.clocks.hue_noise_phase =
        wrap_t(endpoint.clocks.hue_noise_phase + params.color.hue_noise_speed);
    endpoint.clocks.source_noise_time = wrap_t(
        endpoint.clocks.source_noise_time + params.source.noise_time_rate);
    endpoint.clocks.surface_noise_time =
        wrap_t(endpoint.clocks.surface_noise_time + params.surface_noise.rate);
    if (config.slots.warp_program.outer.kind == WarpStageKind::AFFINE_FRAME)
      endpoint.clocks.warp_outer_rotation =
          TWO_PI_F *
          wrap_t((endpoint.clocks.warp_outer_rotation +
                  params.warp.outer.speed * params.warp.outer.rotation) /
                 TWO_PI_F);
    if (config.slots.warp_program.inner.kind == WarpStageKind::AFFINE_FRAME)
      endpoint.clocks.warp_inner_rotation =
          TWO_PI_F *
          wrap_t((endpoint.clocks.warp_inner_rotation +
                  params.warp.inner.speed * params.warp.inner.rotation) /
                 TWO_PI_F);
    endpoint.clocks.warp_outer_phase =
        wrap_t(endpoint.clocks.warp_outer_phase + params.warp.outer.speed);
    endpoint.clocks.warp_inner_phase =
        wrap_t(endpoint.clocks.warp_inner_phase + params.warp.inner.speed);
    endpoint.clocks.palette_oscillation_phase =
        wrap_t(endpoint.clocks.palette_oscillation_phase +
               params.color.phase_oscillation_speed);
    update_spatial_frames(endpoint, config, deltas);
  }

  HS_COLD_MEMBER void prepare_param_morph() {
    if (!state->param_morph.active)
      return;
    const float mix =
        transition_mix(state->param_morph.elapsed, state->param_morph.duration);
    blend.palette_mapping = PaletteMappingWeights::lerp(
        state->param_morph.mapping_from, state->param_morph.mapping_to, mix);
    if (mix == 0.0f)
      blend.params = state->param_morph.from;
    else if (mix == 1.0f)
      blend.params = state->param_morph.to;
    else if (state->param_morph.staggered)
      blend.params.lerp_staggered(state->param_morph.from,
                                  state->param_morph.to, mix, active_slots);
    else
      blend.params.lerp(state->param_morph.from, state->param_morph.to, mix,
                        active_slots);
  }

  static float transition_mix(uint16_t elapsed, uint16_t duration) {
    if (elapsed == 0)
      return 0.0f;
    if (elapsed >= duration)
      return 1.0f;
    return ease_in_out_sin(static_cast<float>(elapsed) / duration);
  }

  __attribute__((noinline)) HS_COLD_MEMBER void apply_requested_config() {
#if HS_ENABLE_PARAM_GUI_BRIDGE
    if (!requested_schema_bound) {
      if (!valid_config(requested_config)) {
        reject_requested_config();
        return;
      }
      accepted_config = requested_config;
    } else {
      const size_t before_count = pending_edit_count;
      refresh_accepted_config();
      if (before_count != pending_edit_count)
        rebind_parameters();
    }
    const Config &next_config = accepted_config;
#else
    const Config &next_config = requested_config;
    if (!valid_config(next_config)) {
      reject_requested_config();
      return;
    }
#endif
    if (!admissible_config(next_config)) {
      reject_requested_config();
      return;
    }
    if (next_config == published_config)
      return;
    if (!prepare_resource_union(next_config, next_config)) {
      reject_requested_config();
      return;
    }
    if (state->transition.active)
      runtime = state->transition.elapsed * 2 < state->transition.duration
                    ? state->transition.from_runtime
                    : state->transition.to_runtime;
    state->transition.active = false;
    state->param_morph.active = false;
    active_slots = next_config.slots;
    active_pipeline = resolve_pipeline_id(next_config);
    blend.params = next_config.params;
    blend.palette_mapping =
        palette_mapping_weights(next_config.slots.palette_mapping);
#if HS_ENABLE_PARAM_GUI_BRIDGE
    display_config = next_config;
#endif
    published_config = next_config;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    accepted_config = next_config;
#endif
    if (!requested_schema_bound)
      rebind_parameters();
  }

  HS_COLD_MEMBER void reject_requested_config() {
    requested_config = published_config;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    accepted_config = published_config;
    pending_edit_count = 0;
    display_config = published_config;
#endif
    rebind_parameters();
  }

  HS_COLD_MEMBER bool try_apply_config(const Config &candidate,
                                       uint16_t duration, bool staggered,
                                       bool continue_choreo) {
    if (!admissible_config(candidate) || duration == 0)
      return false;
    if (state->transition.active)
      return false;
    Config &current = state->render_config;
    current.slots = active_slots;
    current.params = blend.params;
    Config &target = state->transition.to_config;
    target = candidate;
    if (!transition_admitted(current, target))
      return false;
    if (target == current) {
      state->param_morph.active = false;
      blend.palette_mapping =
          palette_mapping_weights(current.slots.palette_mapping);
      return true;
    }
    if (stable_topology(current, target)) {
      if (!prepare_resource_union(current, target))
        return false;
      state->param_morph = {
          current.params,
          target.params,
          blend.palette_mapping,
          palette_mapping_weights(target.slots.palette_mapping),
          target.slots.palette_mapping,
          0,
          duration,
          staggered,
          continue_choreo,
          true};
      return true;
    }
    if (!prepare_resource_union(current, current))
      return false;
    const uint16_t planned_duration =
        (duration & 1U) != 0 ? duration + 1 : duration;
    state->param_morph.active = false;
    TransitionRuntime &transition = state->transition;
    transition.from_config = current;
    transition.from_runtime = runtime;
    transition.to_runtime = runtime;
    transition.elapsed = 0;
    transition.duration = planned_duration;
    transition.continue_choreo = continue_choreo;
    transition.active = true;
    transition.from_pipeline = active_pipeline;
    transition.to_pipeline = resolve_pipeline_id(target);
    return true;
  }

  HS_COLD_MEMBER void finish_transitions() {
    if (state->transition.active) {
      if (state->transition.elapsed == state->transition.duration / 2)
        HS_CHECK(prepare_resource_union(state->transition.to_config,
                                        state->transition.to_config),
                 "through-clear destination resources exceed capacity");
      if (state->transition.elapsed < state->transition.duration) {
        ++state->transition.elapsed;
        return;
      }
      const bool continue_choreo = state->transition.continue_choreo;
      runtime = state->transition.to_runtime;
      active_slots = state->transition.to_config.slots;
      active_pipeline = state->transition.to_pipeline;
      blend.params = state->transition.to_config.params;
      blend.palette_mapping =
          palette_mapping_weights(active_slots.palette_mapping);
      state->transition.active = false;
      if (continue_choreo)
        enter_preset();
      return;
    }
    if (!state->param_morph.active)
      return;
    if (state->param_morph.elapsed < state->param_morph.duration) {
      ++state->param_morph.elapsed;
      return;
    }
    const bool continue_choreo = state->param_morph.continue_choreo;
    active_slots.palette_mapping = state->param_morph.mapping_destination;
    blend.palette_mapping = state->param_morph.mapping_to;
    state->param_morph.active = false;
    if (continue_choreo)
      enter_preset();
  }

  __attribute__((noinline)) HS_COLD_MEMBER void publish_live_config() {
    if (anims_paused || state->transition.active || state->param_morph.active)
      return;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    if (accepted_config != published_config)
      return;
#endif
    published_config = {active_slots, blend.params};
#if HS_ENABLE_PARAM_GUI_BRIDGE
    Config next_requested = published_config;
    for (size_t index = 0; index < pending_edit_count; ++index)
      copy_pending_value(next_requested, requested_config,
                         pending_edits[index]);
    requested_config = next_requested;
#else
    requested_config = published_config;
#endif
#if HS_ENABLE_PARAM_GUI_BRIDGE
    accepted_config = published_config;
#endif
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  HS_COLD_MEMBER void refresh_parameter_display() override {
    if (state->transition.active) {
      const float mix =
          transition_mix(state->transition.elapsed, state->transition.duration);
      display_config.slots = mix < 0.5f ? state->transition.from_config.slots
                                        : state->transition.to_config.slots;
      display_config.params.lerp(state->transition.from_config.params,
                                 state->transition.to_config.params, mix,
                                 display_config.slots);
      return;
    }
    display_config = {active_slots, blend.params};
  }
#endif

  /**
   * @brief Admission test for a requested configuration.
   * @details The simulator accepts every structurally valid configuration;
   * device builds additionally require a compiled inverse pipeline.
   */
  HS_COLD_MEMBER static bool
  admissible_config(const RequestedConfig &candidate) {
    if (!valid_config(candidate))
      return false;
#if HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND
    return true;
#else
    return find_inverse_program(candidate) != nullptr;
#endif
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  const char *begin_warning(const char *format, ...) const {
    va_list args;
    va_start(args, format);
    std::vsnprintf(warning_text.data(), warning_text.size(), format, args);
    va_end(args);
    return warning_text.data();
  }

  void append_warning(const char *format, ...) const {
    const size_t length = std::strlen(warning_text.data());
    if (length >= warning_text.size() - 1)
      return;
    va_list args;
    va_start(args, format);
    std::vsnprintf(warning_text.data() + length, warning_text.size() - length,
                   format, args);
    va_end(args);
  }

  bool append_range_warning(const char *label, float value, float minimum,
                            float maximum) const {
    if (value >= minimum && value <= maximum)
      return false;
    append_warning(" %s %.7g is outside [%.7g, %.7g].", label,
                   static_cast<double>(value), static_cast<double>(minimum),
                   static_cast<double>(maximum));
    return true;
  }

  const char *stage_tuple_warning(const char *position,
                                  const WarpStageSpec &spec,
                                  const WarpStageParams &params) const {
    begin_warning("%s %s rejected.", position, warp_option(spec.kind));
    switch (spec.kind) {
    case WarpStageKind::NONE:
    case WarpStageKind::LEGACY_STEREO_NOISE:
      break;
    case WarpStageKind::AFFINE_FRAME:
      append_range_warning("Translate X", params.translation_x, -4.0f, 4.0f);
      append_range_warning("Translate Y", params.translation_y, -4.0f, 4.0f);
      append_range_warning("Rotation", params.rotation, -TWO_PI_F, TWO_PI_F);
      append_range_warning("Scale X", params.scale_x, 0.25f, 4.0f);
      append_range_warning("Scale Y", params.scale_y, 0.25f, 4.0f);
      append_range_warning("Shear", params.shear, -0.75f, 0.75f);
      break;
    case WarpStageKind::WAVE_SHEAR:
      append_range_warning("Warp Strength", params.strength, -4.0f, 4.0f);
      append_range_warning("Frequency", params.frequency, 0.0f, 64.0f);
      append_range_warning("Warp Speed", params.speed, NOISE_SPEED_MIN,
                           NOISE_SPEED_MAX);
      break;
    case WarpStageKind::VORTEX:
      append_range_warning("Radius", params.radius, 1.0f / 64.0f, 8.0f);
      append_range_warning("Turns", params.turns, -4.0f, 4.0f);
      append_range_warning("Orbit Radius", params.center_orbit_radius, 0.0f,
                           4.0f);
      append_range_warning("Warp Speed", params.speed, NOISE_SPEED_MIN,
                           NOISE_SPEED_MAX);
      break;
    case WarpStageKind::VECTOR_NOISE:
      append_range_warning("Warp Strength", params.strength, 0.0f,
                           VECTOR_WARP_STRENGTH_MAX);
      append_range_warning("Warp Scale", params.scale, 1.0f / 64.0f,
                           VECTOR_WARP_SCALE_MAX);
      append_range_warning("Warp Speed", params.speed, NOISE_SPEED_MIN,
                           NOISE_SPEED_MAX);
      break;
    case WarpStageKind::CURL_FLOW: {
      append_range_warning("Warp Strength", params.strength,
                           -CURL_WARP_STRENGTH_MAX, CURL_WARP_STRENGTH_MAX);
      append_range_warning("Warp Scale", params.scale, 1.0f / 64.0f,
                           CURL_WARP_SCALE_MAX);
      append_range_warning("Warp Speed", params.speed, NOISE_SPEED_MIN,
                           NOISE_SPEED_MAX);
      const float strength_limit = curl_strength_limit(spec, params);
      if (Workbench::abs_value(params.strength) > strength_limit)
        append_warning(
            " %s at Warp Scale %.7g requires |Warp Strength| <= %.9f; "
            "current value is %.7g.",
            CURL_INTEGRATOR_OPTIONS[static_cast<uint8_t>(spec.curl_integrator)],
            static_cast<double>(params.scale),
            static_cast<double>(strength_limit),
            static_cast<double>(params.strength));
      break;
    }
    case WarpStageKind::MIRROR_TILE:
      append_range_warning("Rotation", params.rotation, 0.0f, TWO_PI_F);
      append_range_warning("Cell X", params.cell_x, CELL_MIN, CELL_MAX);
      append_range_warning("Cell Y", params.cell_y, CELL_MIN, CELL_MAX);
      break;
    case WarpStageKind::POLAR_CHART:
      append_range_warning("Radial Scale", params.radial_scale, 1.0f / 64.0f,
                           16.0f);
      break;
    }
    append_warning(" Set every listed control within its stated limit.");
    return warning_text.data();
  }

  const char *program_bounds_warning(const Config &candidate) const {
    float bound = projection_coordinate_bound(candidate);
    const Complex source_period = source_cartesian_period(candidate);
    const WarpStageSpec stages[] = {candidate.slots.warp_program.outer,
                                    candidate.slots.warp_program.inner};
    const WarpStageParams params[] = {candidate.params.warp.outer,
                                      candidate.params.warp.inner};
    const char *positions[] = {"Planar Warp 1", "Planar Warp 2"};
    for (size_t index = 0; index < 2; ++index) {
      if (stages[index].kind == WarpStageKind::VECTOR_NOISE ||
          stages[index].kind == WarpStageKind::CURL_FLOW) {
        const float lattice_bound = params[index].scale * (bound + 100.0f);
        if (lattice_bound > NOISE_LATTICE_LIMIT) {
          const float scale_limit = NOISE_LATTICE_LIMIT / (bound + 100.0f);
          return begin_warning(
              "%s %s rejected: Warp Scale %.7g produces noise coordinate "
              "bound %.7g above %.7g. Set Warp Scale <= %.7g or choose a "
              "projection/lens with a smaller coordinate extent.",
              positions[index], warp_option(stages[index].kind),
              static_cast<double>(params[index].scale),
              static_cast<double>(lattice_bound),
              static_cast<double>(NOISE_LATTICE_LIMIT),
              static_cast<double>(scale_limit));
        }
      }
      bound = stage_coordinate_bound(stages[index], params[index], bound,
                                     source_period);
      if (bound > WARP_COORD_LIMIT)
        return begin_warning(
            "%s %s rejected: its predicted coordinate bound %.7g exceeds "
            "%.7g. Reduce this warp's displacement/translation controls or "
            "choose a projection/lens with a smaller coordinate extent.",
            positions[index], warp_option(stages[index].kind),
            static_cast<double>(bound), static_cast<double>(WARP_COORD_LIMIT));
    }
    const float source_bound = candidate.params.source.noise_scale * bound;
    return begin_warning(
        "Noise Contour rejected: Source Noise Scale %.7g produces noise "
        "coordinate bound %.7g above %.7g. Set Source Noise Scale <= %.7g "
        "or reduce the preceding warp extent.",
        static_cast<double>(candidate.params.source.noise_scale),
        static_cast<double>(source_bound),
        static_cast<double>(NOISE_LATTICE_LIMIT),
        static_cast<double>(NOISE_LATTICE_LIMIT / bound));
  }

  const char *resource_warning() const {
    return begin_warning(
        "The active noise consumers exceed the resource limit of %u. Disable "
        "one noise Function, Lens, or Warp.",
        static_cast<unsigned>(MAX_NOISE_RESOURCES));
  }

  const char *admission_warning(const Config &candidate,
                                const char *edited_name) const {
    const WarpStageSpec &outer = candidate.slots.warp_program.outer;
    const WarpStageSpec &inner = candidate.slots.warp_program.inner;
    if (is_sphere_source(candidate.slots.function) &&
        outer.kind != WarpStageKind::NONE && inner.kind != WarpStageKind::NONE)
      return begin_warning(
          "%s rejects Planar Warp 1 %s and Planar Warp 2 %s. Set both warps "
          "to None, or select a plane-space Function.",
          FUNCTION_OPTIONS[static_cast<uint8_t>(candidate.slots.function)],
          warp_option(outer.kind), warp_option(inner.kind));
    if (is_sphere_source(candidate.slots.function) &&
        (outer.kind != WarpStageKind::NONE ||
         inner.kind != WarpStageKind::NONE)) {
      const bool outer_active = outer.kind != WarpStageKind::NONE;
      const char *position = outer_active ? "Planar Warp 1" : "Planar Warp 2";
      const WarpStageKind kind = outer_active ? outer.kind : inner.kind;
      return begin_warning(
          "%s rejects %s %s. Set %s to None, or select a plane-space "
          "Function.",
          FUNCTION_OPTIONS[static_cast<uint8_t>(candidate.slots.function)],
          position, warp_option(kind), position);
    }
    if (outer.kind == WarpStageKind::POLAR_CHART &&
        inner.kind != WarpStageKind::NONE &&
        inner.kind != WarpStageKind::WAVE_SHEAR)
      return begin_warning(
          "Planar Warp 1 Polar Chart cannot run while Planar Warp 2 is %s. Set "
          "Planar Warp 2 to None or Wave Shear, or choose a different Planar "
          "Warp 1.",
          warp_option(inner.kind));
    if (inner.kind == WarpStageKind::POLAR_CHART &&
        outer.kind != WarpStageKind::NONE)
      return begin_warning(
          "Planar Warp 2 Polar Chart cannot run while Planar Warp 1 is %s. Set "
          "Planar Warp 1 to None or choose a different Planar Warp 2.",
          warp_option(outer.kind));
    const WarpStageSpec *polar =
        outer.kind == WarpStageKind::POLAR_CHART   ? &outer
        : inner.kind == WarpStageKind::POLAR_CHART ? &inner
                                                   : nullptr;
    if (polar != nullptr && !polar_source_compatible(candidate, *polar)) {
      const char *position =
          polar == &outer ? "Planar Warp 1" : "Planar Warp 2";
      const SourceTraits traits = source_traits(candidate.slots.function);
      if (!traits.y_periodic || !traits.polar_angle_compatible)
        return begin_warning(
            "%s Polar Chart requires a polar-periodic Function; %s is not "
            "compatible. Select Grid or Primitive Lattice, or "
            "choose another %s.",
            position,
            FUNCTION_OPTIONS[static_cast<uint8_t>(candidate.slots.function)],
            position);
      const float periods = polar_seam_periods(candidate, *polar);
      const float nearest_periods = floorf(periods + 0.5f);
      if (candidate.slots.function == Function::PRIMITIVE_LATTICE)
        return begin_warning(
            "%s Polar Chart requires 2*pi x Lattice Cell Scale x Polar "
            "Harmonic to be a whole number. %.7g x %u gives %.7g. Set Lattice "
            "Cell Scale to %.7g or change %s Polar Harmonic.",
            position,
            static_cast<double>(candidate.params.source.lattice_cell_scale),
            static_cast<unsigned>(polar->polar_harmonic),
            static_cast<double>(periods),
            static_cast<double>(
                nearest_periods /
                (TWO_PI_F * static_cast<float>(polar->polar_harmonic))),
            position);
      const float suggested_frequency =
          nearest_periods / static_cast<float>(polar->polar_harmonic);
      return begin_warning(
          "%s Polar Chart requires Pattern Freq x Polar Harmonic to be a "
          "whole number. %.7g x %u = %.7g. Set Pattern Freq to %.7g or change "
          "%s Polar Harmonic.",
          position, static_cast<double>(candidate.params.source.pattern_freq),
          static_cast<unsigned>(polar->polar_harmonic),
          static_cast<double>(periods),
          static_cast<double>(suggested_frequency), position);
    }
    if (!affine_translation_compatible(candidate)) {
      const bool outer_scroll =
          affine_has_translation(outer, candidate.params.warp.outer);
      const char *position = outer_scroll ? "Planar Warp 1" : "Planar Warp 2";
      const WarpStageParams &params = outer_scroll
                                          ? candidate.params.warp.outer
                                          : candidate.params.warp.inner;
      if (!Workbench::whole_affine_winding(params.translation_x) ||
          !Workbench::whole_affine_winding(params.translation_y))
        return begin_warning(
            "%s Affine Frame translation must use whole source-cell windings. "
            "Set Translation X and Translation Y to whole numbers.",
            position);
      if (candidate.slots.function != Function::PRIMITIVE_LATTICE)
        return begin_warning(
            "%s Affine Frame translation requires an exactly periodic "
            "Function. Select Primitive Lattice or set both translations to "
            "zero.",
            position);
      if (outer_scroll && inner.kind != WarpStageKind::NONE)
        return begin_warning(
            "Planar Warp 1 Affine Frame translation cannot precede Planar "
            "Warp 2 %s because the later warp breaks its source-period seam. "
            "Set Planar Warp 2 to None or set both translations to zero.",
            warp_option(inner.kind));
      return begin_warning(
          "%s Affine Frame translation cannot drive Total Warp Displacement "
          "hue because its path length resets at the source-period seam. "
          "Select Hue Shift None or Noise, or set both translations to zero.",
          position);
    }
    if (!strict_seam_compatible(candidate)) {
      begin_warning(
          "Projection %s requires seam-safe stages.",
          PROJECTION_OPTIONS[static_cast<uint8_t>(candidate.slots.projection)]);
      if (candidate.slots.function == Function::NOISE_CONTOUR)
        append_warning(" Function Noise Contour (Projected) is not seam-safe.");
      if (seam_sensitive_warp(outer.kind))
        append_warning(" Planar Warp 1 %s is not seam-safe.",
                       warp_option(outer.kind));
      if (seam_sensitive_warp(inner.kind))
        append_warning(" Planar Warp 2 %s is not seam-safe.",
                       warp_option(inner.kind));
      append_warning(" Replace the named stage or select Folded Sinusoidal, "
                     "Stereographic, Gnomonic, or Equirectangular.");
      return warning_text.data();
    }
    const SurfaceNoiseParams &surface_noise = candidate.params.surface_noise;
    const float minimum_surface_strength =
        candidate.slots.surface_noise == SurfaceNoise::CURL ? -0.5f : 0.0f;
    if (surface_noise.scale < LENS_NOISE_SCALE_MIN ||
        surface_noise.scale > LENS_NOISE_SCALE_MAX ||
        surface_noise.strength < minimum_surface_strength ||
        surface_noise.strength > 0.5f || surface_noise.rate < NOISE_RATE_MIN ||
        surface_noise.rate > NOISE_RATE_MAX || surface_noise.direction < 0.0f ||
        surface_noise.direction > 1.0f) {
      begin_warning("Surface Noise %s rejected.",
                    SURFACE_NOISE_OPTIONS[static_cast<uint8_t>(
                        candidate.slots.surface_noise)]);
      append_range_warning("Surface Noise Scale", surface_noise.scale,
                           LENS_NOISE_SCALE_MIN, LENS_NOISE_SCALE_MAX);
      append_range_warning("Surface Noise Strength", surface_noise.strength,
                           minimum_surface_strength, 0.5f);
      append_range_warning("Surface Noise Rate", surface_noise.rate,
                           NOISE_RATE_MIN, NOISE_RATE_MAX);
      append_range_warning("Surface Noise Direction", surface_noise.direction,
                           0.0f, 1.0f);
      append_warning(" Set the named Surface Noise control within its range.");
      return warning_text.data();
    }
    if (!preset_in_ranges(candidate)) {
      const ParamDef *parameter = getParameters().find(edited_name);
      if (parameter != nullptr)
        return begin_warning(
            "%s %.7g is outside its registered range [%.7g, %.7g]. Set %s "
            "within that range.",
            edited_name, static_cast<double>(parameter->get_requested()),
            static_cast<double>(parameter->min),
            static_cast<double>(parameter->max), edited_name);
    }
    if (!valid_stage_tuple(outer, candidate.params.warp.outer))
      return stage_tuple_warning("Planar Warp 1", outer,
                                 candidate.params.warp.outer);
    if (!valid_stage_tuple(inner, candidate.params.warp.inner))
      return stage_tuple_warning("Planar Warp 2", inner,
                                 candidate.params.warp.inner);
    if (!safe_program_bounds(candidate))
      return program_bounds_warning(candidate);
    if (candidate.slots.surface_lens == SurfaceLens::MOBIUS &&
        !Workbench::valid_mobius(candidate.params.surface_lens.mobius)) {
      const MobiusParams &m = candidate.params.surface_lens.mobius;
      const float det_re =
          m.a.re * m.d.re - m.a.im * m.d.im - m.b.re * m.c.re + m.b.im * m.c.im;
      const float det_im =
          m.a.re * m.d.im + m.a.im * m.d.re - m.b.re * m.c.im - m.b.im * m.c.re;
      return begin_warning(
          "Mobius Lens rejected: |A*D - B*C| is %.7g; it must be at least "
          "0.001. Adjust the requested Mobius coefficient until the "
          "determinant reaches 0.001 or more.",
          static_cast<double>(sqrtf(det_re * det_re + det_im * det_im)));
    }
    if (!resource_union_fits(candidate, candidate))
      return resource_warning();
    if (!HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND &&
        find_inverse_program(candidate) == nullptr)
      return uncompiled_program_warning(candidate, edited_name);
    return begin_warning(
        "%s was rejected by an unclassified ShaderWorkbench admission rule. Keep "
        "the requested value and report this exact configuration as a bug.",
        edited_name);
  }

  const char *uncompiled_program_warning(const Config &candidate,
                                         const char *edited_name) const {
    const TopologyKey key = make_topology_key(candidate);
    for (const ProgramDescriptor &program : inverse_programs())
      if (program.key == key)
        return begin_warning(
            "%s is outside what the compiled pipeline for this stage "
            "combination supports. Restore %s or change a stage.",
            edited_name, edited_name);
    return begin_warning(
        "This stage combination has no compiled pipeline. Restore %s or "
        "select a combination reachable from a preset.",
        edited_name);
  }
#endif

  HS_COLD_MEMBER static constexpr Choreo preset_choreo() {
#ifdef HS_PROFILE_SHADER_WORKBENCH_FAST_CYCLE
    return {32, 32, 2, false};
#else
    return CHOREO;
#endif
  }

  HS_COLD_MEMBER void enter_preset() {
    if (preset_count_for_view() < 2) {
      preset_dwell_remaining = 0;
      preset_dwell_armed = false;
      return;
    }
    const Choreo choreo = preset_choreo();
    preset_dwell_remaining = static_cast<uint16_t>(
        hs::rand_int(choreo.dwell_min, choreo.dwell_max + 1));
    preset_dwell_armed = true;
  }

  HS_COLD_MEMBER void advance_preset_choreography() {
    if (anims_paused || !preset_dwell_armed)
      return;
    if (preset_dwell_remaining > 0 && --preset_dwell_remaining > 0)
      return;
    preset_dwell_armed = false;
    begin_blend();
  }

  HS_COLD_MEMBER void begin_blend() {
    if (advancePreset()) {
#if defined(HS_PROFILE_ENABLE)
      hs::log("Preset: %u/%u", static_cast<unsigned>(getPresetIndex()),
              static_cast<unsigned>(preset_count_for_view()));
#endif
    } else {
      preset_dwell_remaining = 1;
      preset_dwell_armed = true;
    }
  }

  static void next_generated_palette(uint32_t &hue, uint32_t sequence,
                                     PaletteHarmony harmony, float chroma,
                                     GenerativePalette &out) {
    EffectPaletteRecipes::GeneratedPaletteBank::next_palette(
        hue, sequence, harmony, chroma, out);
  }

  static constexpr uint32_t HUE_STEP =
      EffectPaletteRecipes::GeneratedPaletteBank::HUE_STEP;
  static constexpr size_t PARAM_CAPACITY = 80;

  static constexpr auto &GENERATED_SURFACE_NOISE_SLOTS =
      Workbench::GENERATED_SURFACE_NOISE_SLOTS;
  static constexpr auto &OUTER_WARP_PARAM_NAMES =
      Workbench::OUTER_WARP_PARAM_NAMES;
  static constexpr auto &INNER_WARP_PARAM_NAMES =
      Workbench::INNER_WARP_PARAM_NAMES;
  static constexpr Choreo CHOREO{0, 0, 480, false};

  Orientation<> projection_walk;
  Orientation<> outer_walk;
  Timeline timeline;
  size_t prepared_noise_count = 0;
  StateBundle *state = nullptr;

  Quaternion base_orientation =
      make_rotation(Vector(0, 0, -1), Vector(0, -1, 0));
  Quaternion projection_walk_prev;
  Quaternion outer_walk_prev;

  EffectPaletteRecipes::GeneratedPaletteBank generated_palettes;

  Slots active_slots = PRESETS[0].config.slots;
  InversePipelineId active_pipeline = PRESETS[0].pipeline;
#if HS_ENABLE_PARAM_GUI_BRIDGE
  Config display_config = PRESETS[0].config;
  std::array<PendingEdit, PARAM_CAPACITY> pending_edits{};
  size_t pending_edit_count = 0;
  mutable std::array<char, 1024> warning_text{};
#endif
  RequestedConfig requested_config = PRESETS[0].config;
  Config published_config = PRESETS[0].config;
#if HS_ENABLE_PARAM_GUI_BRIDGE
  Config accepted_config = PRESETS[0].config;
#endif
  bool requested_schema_bound = false;
  bool registered_range_clamped = false;
  bool fixed_topology = false;
  std::span<const uint8_t> preset_view{};
  uint16_t preset_dwell_remaining = 0;
  bool preset_dwell_armed = false;
  Blend blend{PRESETS[0].config.params,
              palette_mapping_weights(PRESETS[0].config.slots.palette_mapping)};
  EndpointRuntime runtime;
#if HS_ENABLE_TEST_HOOKS
  uint32_t walk_step_count = 0;
  uint32_t generated_palette_step_count = 0;
#endif
#if defined(HS_PROFILE_ENABLE)
  bool profile_program_valid = false;
  size_t profile_program_preset = 0;
  InversePipelineId profile_program_pipeline = InversePipelineId::NONE;
  ProfileEndpoint profile_program_endpoint = ProfileEndpoint::STEADY;
#endif

  static constexpr size_t FOOTPRINT_BYTES =
      gamut_lut_bytes(GAMUT_ANGLE_STEPS, GAMUT_L_STEPS) +
      EffectPaletteRecipes::GeneratedPaletteBank::required_arena_bytes() +
      PARAM_CAPACITY * sizeof(ParamDef) + sizeof(StateBundle) +
      alignof(StateBundle);
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
      "ShaderWorkbench persistent footprint exceeds the default partition");
};

template <int W, int H> using Shader = ShaderWorkbench<W, H>;

#include "core/control/registry.h"
REGISTER_EFFECT(Shader)

#endif // HS_ENABLE_SHADER_WORKBENCH
