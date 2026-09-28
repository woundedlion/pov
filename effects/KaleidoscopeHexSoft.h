/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using KaleidoscopeHexSoftParams =
    Pullback::Params<Pullback::TwinWaveSourceParams, Pullback::NoWarpParams,
                     Pullback::MirrorParams>;
using KaleidoscopeHexSoftSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                   Pullback::Lens::Kaleidoscope, Pullback::TransferKind::NONE,
                   Pullback::ProjectionCoverageMode::WEIGHT_SQUARED>;

/**
 * @brief A drifting twin wave reflected through a kaleidoscope.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopeHexSoft
    : public Pullback::ComposedEffect<
          W, H, KaleidoscopeHexSoft<W, H>, KaleidoscopeHexSoftParams,
          KaleidoscopeHexSoftSpec, PaletteHarmony::TRIADIC,
          Pullback::HueMode::NOISE, Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = KaleidoscopeHexSoftParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-hex-soft";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "39b4a45b68feb15150c7baa8d47ddece27b82fe9c859763709e4a32b8af0aaed";
  static constexpr std::string_view PRESET_BANK_DIGEST = "a54f426490227b6ad95b64d0a18bbd0bb65beaab6fb4053f37f85bb435a6e7ed";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "twin-wave"
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
    value.color.hue_noise_scale = 2.2033439f;
    value.color.hue_noise_speed = -0.00040800002f;
    value.color.hue_shift_amount = 0.27f;
    value.color.mapping_frequency = 1.0f;
    value.color.mapping_phase = 0.0f;
    value.color.palette_chroma = 0.361f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.color.phase_oscillation_depth = 0.0f;
    value.color.phase_oscillation_speed = 0.0f;
    value.color.opacity_high = 1.0f;
    value.color.opacity_low = 1.0f;
    value.projection.spin_rate = 0.0f;
    value.projection.wander = 1.0f;
    value.projection.singularity_fade = 4.971f;
    value.source.angle_rate = 0.05f;
    value.source.secondary_rate = 0.8f;
    value.source.pattern_freq = 4.9755f;
    value.source.speed = 0.125f;
    value.inner_warp.cell_x = 1.0f;
    value.inner_warp.cell_y = 1.0f;
    value.inner_warp.offset_x = 0.0f;
    value.inner_warp.offset_y = 0.0f;
    value.inner_warp.rotation = 0.0f;
    value.inner_warp.speed = 0.0f;
    return value;
  }
  // clang-format on
  // End generated params.
};
