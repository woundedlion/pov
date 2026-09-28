/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using KaleidoscopeFlowersParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::NoWarpParams,
                     Pullback::MirrorParams>;
using KaleidoscopeFlowersSpec =
    Pullback::Spec<Pullback::ProjectionKind::EQUIRECTANGULAR,
                   Pullback::Lens::DodecahedralKaleidoscope,
                   Pullback::TransferKind::NONE,
                   Pullback::ProjectionCoverageMode::WEIGHT_SQUARED>;

/**
 * @brief Dodecahedral grids mapped continuously around the equator.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopeFlowers
    : public Pullback::ComposedEffect<
          W, H, KaleidoscopeFlowers<W, H>, KaleidoscopeFlowersParams,
          KaleidoscopeFlowersSpec, PaletteHarmony::ANALOGOUS,
          Pullback::HueMode::NOISE, Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = KaleidoscopeFlowersParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-flowers";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "ea710e1c3a212e24d0f0dc2f5166d9c1cc1b57443e727bbc92a41a449e48a229";
  static constexpr std::string_view PRESET_BANK_DIGEST = "d1bc04d4ed1aec4b7720a5d90429388da1933b83cd85b387ac56da440dc51970";
  static constexpr std::array<std::string_view, 3> PRESET_IDS{
      "double-map",
      "open-grid",
      "fine-grid"
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
    value.color.hue_noise_scale = 1.4721563f;
    value.color.hue_noise_speed = 0.0f;
    value.color.hue_shift_amount = 0.366f;
    value.color.mapping_frequency = 2.0f;
    value.color.mapping_phase = 0.0f;
    value.color.palette_chroma = 1.0f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.color.phase_oscillation_depth = 0.0f;
    value.color.phase_oscillation_speed = 0.0f;
    value.color.opacity_high = 1.0f;
    value.color.opacity_low = 1.0f;
    value.projection.central_meridian = 0.0f;
    value.projection.spin_rate = 0.0f;
    value.projection.wander = 0.165f;
    value.projection.singularity_fade = 2.14f;
    value.source.angle_rate = 0.026999999f;
    value.source.complexity = 3.0f;
    value.source.secondary_rate = 0.8f;
    value.source.pattern_freq = 3.9407f;
    value.source.pattern_mix = 1.0f;
    value.source.speed = 0.0f;
    value.inner_warp.cell_x = 1.0471976f;
    value.inner_warp.cell_y = 0.99770314f;
    value.inner_warp.offset_x = 0.0f;
    value.inner_warp.offset_y = 0.0f;
    value.inner_warp.rotation = 0.0f;
    value.inner_warp.speed = 0.00013f;
    return value;
  }

  /** @brief Params for the preset at index in PRESET_IDS. */
  static constexpr Params preset_params(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.color.mapping_frequency = 1.0f;
    }
    if (index == 2) {
      value.color.mapping_frequency = 21.212f;
      value.source.pattern_freq = 0.3985f;
      value.inner_warp.cell_y = 0.90189064f;
      value.inner_warp.speed = 0.00058f;
    }
    return value;
  }
  // clang-format on
  // End generated params.
};
