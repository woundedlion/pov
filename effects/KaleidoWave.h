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
                   Pullback::Lens::Kaleidoscope, Pullback::TransferKind::LINEAR,
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
      "a2fa9dd816a24f76242c6614bea3c7b74f768c3a775398ee63daae25da97b292";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "262b8f7bc132ecbf75064196f7f0971c28a212f69682f5fd856aab71765fc92d";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"twin-wave"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;

  static constexpr Params initial_params() {
    Params value;
    value.source = {4.9755f, 0.125f, 0.8f, 0.05f};
    value.projection.pole_fade = 4.971f;
    value.projection.wander = 1.0f;
    value.projection.camera_wander = 1.0f;
    value.color.hue_shift_amount = 0.27f;
    value.color.hue_noise_scale = 2.2033439f;
    value.color.hue_noise_speed = -0.00040800002f;
    value.color.palette_chroma = 0.361f;
    return value;
  }
  static constexpr Params preset_params(size_t) { return initial_params(); }
};

#include "core/control/registry.h"
REGISTER_EFFECT(KaleidoWave)
