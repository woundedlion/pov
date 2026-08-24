/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using KaleidoscopeHexBrightParams =
    Pullback::Params<Pullback::TwinWaveSourceParams, Pullback::NoWarpParams,
                     Pullback::MirrorParams>;
using KaleidoscopeHexBrightSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                   Pullback::Lens::HexagonalPrismKaleidoscope,
                   Pullback::TransferKind::NONE,
                   Pullback::CoverageKind::PROJECTION_SQUARED>;

template <int W, int H>
class KaleidoscopeHexBright
    : public Pullback::ComposedEffect<
          W, H, KaleidoscopeHexBright<W, H>, KaleidoscopeHexBrightParams,
          KaleidoscopeHexBrightSpec, PaletteHarmony::ANALOGOUS,
          Pullback::HueMode::NOISE, Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = KaleidoscopeHexBrightParams;
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-hex-bright";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "cd23f8b819db5a55e1b3c25f0713060214254b99531b6eec5f7d0607a20acc54";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "a69686e051b4d46949179a00b82e08df04f1b5abd7b5317376da86753670719b";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "hex-twin-wave", "hex-twin-wave-alt"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;

  static constexpr Params initial_params() {
    Params value;
    value.source = {.pattern_freq = 3.881f,
                    .speed = 0.128598228f,
                    .secondary_rate = 0.8f,
                    .angle_rate = 0.027f};
    value.projection.singularity_fade = 4.971f;
    value.projection.wander = 1.0f;
    value.projection.camera_wander = 1.0f;
    value.color.hue_shift_amount = 0.226f;
    value.color.hue_noise_scale = 1.47215629f;
    value.color.hue_noise_speed = 0.000138f;
    value.color.palette_chroma = 1.0f;
    value.color.mapping_frequency = 1.341f;
    value.color.mapping_phase = -1.0f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::BELL;
    return value;
  }
  static constexpr Params preset_params(size_t index) {
    Params value = initial_params();
    if (index == 1)
      value.color.mapping_frequency = 2.0f;
    return value;
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(KaleidoscopeHexBright)
