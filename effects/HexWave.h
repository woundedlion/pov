/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using HexWaveParams =
    Pullback::Params<Pullback::TwinWaveSourceParams, Pullback::NoWarpParams,
                     Pullback::MirrorParams>;
using HexWaveSpec = Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                                   Pullback::Lens::HexagonalPrismKaleidoscope,
                                   Pullback::TransferKind::NONE,
                                   Pullback::CoverageKind::PROJECTION_SQUARED>;

template <int W, int H>
class HexWave : public Pullback::ComposedEffect<
                    W, H, HexWave<W, H>, HexWaveParams, HexWaveSpec,
                    PaletteHarmony::ANALOGOUS, Pullback::HueMode::NOISE,
                    Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = HexWaveParams;
  static constexpr std::string_view EFFECT_ID = "hex-wave";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "7c973d2475a99fdf61684f9e8759ace540a7eb616498b25eab0b8e91c81287bd";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "eaf259250c058dad158f79e5bb727f8fd680962a9e45893ce487d3c83b11c25a";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "hex-twin-wave", "hex-twin-wave-alt"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;

  static constexpr Params initial_params() {
    Params value;
    value.source = {3.881f, 0.128598228f, 0.8f, 0.027f};
    value.projection.pole_fade = 4.971f;
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
REGISTER_EFFECT(HexWave)
