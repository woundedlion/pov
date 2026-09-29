/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using KaleidoscopeHexOilParams =
    Pullback::Params<Pullback::SpiralSourceParams, Pullback::NoWarpParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::NoValueParams, Pullback::DirectSurfaceParams>;
using KaleidoscopeHexOilSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                   Pullback::Lens::HexagonalPrismKaleidoscope,
                   Pullback::TransferKind::NONE,
                   Pullback::ProjectionCoverageMode::WEIGHT_SQUARED>;

/**
 * @brief A rotating spiral folded through a hexagonal prism kaleidoscope.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopeHexOil
    : public Pullback::ComposedEffect<
          W, H, KaleidoscopeHexOil<W, H>, KaleidoscopeHexOilParams,
          KaleidoscopeHexOilSpec, PaletteHarmony::TRIADIC,
          Pullback::HueMode::PATH_LENGTH,
          Pullback::Color::BrightnessEnvelope::NONE,
          /*AnimatedProjection=*/true, Pullback::SurfacePlacement::AFTER_LENS> {

public:
  using Params = KaleidoscopeHexOilParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-hex-oil";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "d4cec300e807a48008dcd128903c97f00ca36d582a474e94231389f2ffa9a7c7";
  static constexpr std::string_view PRESET_BANK_DIGEST = "302f86fbdcd6122935cbc7e0c80bcff24f441401e3d44d39ec096e02849fcd18";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "kaleidoscope-hex-oil",
      "kaleidoscope-hex-oil-2"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.projection.camera_wander = 1.0f;
    value.color.hue_shift_amount = -2.216f;
    value.color.mapping_frequency = 1.0f;
    value.color.mapping_phase = 0.0f;
    value.color.palette_chroma = 0.78f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.color.phase_oscillation_depth = 0.0f;
    value.color.phase_oscillation_speed = 0.0f;
    value.color.opacity_high = 1.0f;
    value.color.opacity_low = 1.0f;
    value.projection.spin_rate = 0.0f;
    value.projection.wander = 1.0f;
    value.projection.singularity_fade = 1.627f;
    value.source.angle_rate = 0.03f;
    value.source.pattern_freq = 5.5327f;
    value.source.speed = 0.125f;
    value.surface.direction = 0.0f;
    value.surface.scale = 5.740422f;
    value.surface.speed = 0.0f;
    value.surface.strength = 0.4185f;
    return value;
  }

  /** @brief The preset at index in PRESET_IDS and how it departs. */
  HS_COLD_MEMBER static constexpr PresetEntry<Params> preset(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.color.mapping_frequency = 1.2f;
      value.surface.scale = 3.6627343f;
    }
    return {value, Segue::Preset::Lerp{480, math::ease_in_out_sin}};
  }
  // clang-format on
  // End generated params.
};
