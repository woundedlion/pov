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
      "08122928d7ecabd79b7286bdd2788b9ae3d87a1502aef4f53bcdbec5f08e1421";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "655c4529ec3d0b3552dc1ac8ff599484633d91ac2c952a1980ff279ebfe78708";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "hex-twin-wave", "hex-twin-wave-alt"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  // The 6 s Phantasm slot is 96 frames at 16 fps and the effect is rebuilt each
  // visit: 2 dwells plus 1 crossfade must fit, or the later presets never
  // render.
  static constexpr Segue::Lerp PRESET_SEGUE{12, ease_in_out_sin};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 42;
  static constexpr bool ANIMATED_PROJECTION = true;

  static constexpr Params initial_params() {
    Params value;
    value.source = {.pattern_freq = 3.881f,
                    .speed = 0.128598228f,
                    .secondary_rate = 0.8f,
                    .angle_rate = 0.027f};
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
