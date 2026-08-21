/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using KaleidoWaveParams =
    Pullback::Params<Pullback::TwinWaveSourceParams, Pullback::NoWarpParams,
                     Pullback::MirrorParams>;
using KaleidoWaveSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                   Pullback::Lens::Kaleidoscope, Pullback::TransferKind::NONE,
                   Pullback::CoverageKind::PROJECTION_SQUARED>;

template <int W, int H>
class KaleidoWave
    : public Pullback::ComposedEffect<
          W, H, KaleidoWave<W, H>, KaleidoWaveParams, KaleidoWaveSpec,
          PaletteHarmony::TRIADIC, Pullback::HueMode::NOISE,
          Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = KaleidoWaveParams;
  static constexpr std::string_view EFFECT_ID = "kaleido-wave";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "099b9d82b05ab7c5ba7c9d074c9038f7b157257ec45db5acaaeea32c6258e603";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "303b888bcbcee0ef482af7c75ae461ab26a2ec9b732324d9af1e5c650e092ea7";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"twin-wave"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;

  static constexpr Params initial_params() {
    Params value;
    value.source = {.pattern_freq = 4.9755f,
                    .speed = 0.125f,
                    .secondary_rate = 0.8f,
                    .angle_rate = 0.05f};
    value.projection.pole_fade = 4.971f;
    value.projection.wander = 1.0f;
    value.projection.camera_wander = 1.0f;
    value.color.hue_shift_amount = 0.27f;
    value.color.hue_noise_scale = 2.2033439f;
    value.color.hue_noise_speed = -0.00040800002f;
    value.color.palette_chroma = 0.361f;
    return value;
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(KaleidoWave)
