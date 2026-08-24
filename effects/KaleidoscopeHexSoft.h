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
                   Pullback::CoverageKind::PROJECTION_SQUARED>;

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
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-hex-soft";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "064eeecabf4683fa489044dd8ccdbf6b9212dedc16d032bdd5cbb7b95fe5fbb7";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "a54f426490227b6ad95b64d0a18bbd0bb65beaab6fb4053f37f85bb435a6e7ed";
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
    value.projection.singularity_fade = 4.971f;
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
REGISTER_EFFECT(KaleidoscopeHexSoft)
