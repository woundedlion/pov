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
                                   Pullback::TransferKind::LINEAR,
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
      "bdfbc95d955407ab65e712263b5beb11182d922a27a78cec5ea2b9aa2167c86a";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "862e4ded04966682973bd21e47a8af1fdc08b0e007838b42ff96fc116ec33b4c";
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
