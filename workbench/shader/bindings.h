/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_SHADER_WORKBENCH

/**
 * @file bindings.h
 * @brief The pullback binding the compiled pipelines are parameterized on, and the per-stage state providers that read one FrameState.
 */

#include "workbench/shader/kernels.h"

namespace Workbench {

struct ShaderWorkbenchInstrumentation {
#ifdef HS_PROFILE_SHADER_WORKBENCH_STAGES
  using Token = uint32_t;

  __attribute__((always_inline)) static Token mark() { return HS_OS_CYCLES(); }

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
  using FrameState = Workbench::FrameState;
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

inline constexpr PaletteMappingWeights
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
        static_cast<uint8_t>(Pullback::Color::BrightnessEnvelope::ASCENDING) &&
    static_cast<uint8_t>(BrightnessEnvelope::DESCENDING) ==
        static_cast<uint8_t>(Pullback::Color::BrightnessEnvelope::DESCENDING));
static_assert(static_cast<uint8_t>(HueShiftMode::NONE) ==
                  static_cast<uint8_t>(Pullback::Color::HueMode::NONE) &&
              static_cast<uint8_t>(HueShiftMode::NOISE) ==
                  static_cast<uint8_t>(Pullback::Color::HueMode::NOISE) &&
              static_cast<uint8_t>(HueShiftMode::WARP_DISPLACEMENT) ==
                  static_cast<uint8_t>(Pullback::Color::HueMode::PATH_LENGTH));

} // namespace Workbench

#endif // HS_ENABLE_SHADER_WORKBENCH
