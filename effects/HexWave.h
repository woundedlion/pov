/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/look.h"

using HexWaveParams = Looks::Params<Looks::TwinWaveSourceParams,
                                    Looks::NoWarpParams, Looks::MirrorParams>;
using HexWaveSpec = Looks::Spec<Looks::ProjectionKind::STEREOGRAPHIC,
                                Pullback::Lens::HexagonalPrismKaleidoscope,
                                Looks::TransferKind::LINEAR,
                                Looks::CoverageKind::PROJECTION_SQUARED>;

template <int W, int H>
class HexWave
    : public Looks::Composed<W, H, HexWave<W, H>, HexWaveParams, HexWaveSpec,
                             PaletteHarmony::ANALOGOUS, Looks::HueMode::NOISE,
                             Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = HexWaveParams;
  static constexpr std::string_view EFFECT_ID = "hex-wave";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "7c17e826920d96ebe397cd7939bbcc0a31bc431b6ff6545d99472a280ba5cbff";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "8d76e45aa0e6a6c76644991541bf513c531ebf423112fb36cdc93fc22516df87";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"hex-twin-wave"};
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
  static constexpr Params preset_params(size_t) { return initial_params(); }
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(HexWave)
