/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using FacetWaveParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::WaveShearParams,
                     Pullback::MirrorParams>;
using FacetWaveSpec =
    Pullback::Spec<Pullback::ProjectionKind::GNOMONIC_FOLDED,
                   Pullback::Lens::DodecahedralKaleidoscope,
                   Pullback::TransferKind::NONE,
                   Pullback::CoverageKind::PROJECTION_SQUARED>;

template <int W, int H>
class FacetWave : public Pullback::ComposedEffect<
                      W, H, FacetWave<W, H>, FacetWaveParams, FacetWaveSpec,
                      PaletteHarmony::TRIADIC, Pullback::HueMode::NOISE,
                      Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = FacetWaveParams;
  static constexpr std::string_view EFFECT_ID = "facet-wave";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "6200c7a1dc288c154c5aaa1d3fef403914fecd928e81a525e2314188d7f9ca66";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "ae386b772fb2c46f749384cf309d4282d2f82cc75fe82a5dd5abf98e8fa1dcd9";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"wave-mirror",
                                                              "cup-hue"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  // The 6 s Phantasm slot is 96 frames at 16 fps and the effect is rebuilt each
  // visit: 2 dwells plus 1 crossfade must fit, or the later presets never
  // render.
  static constexpr Segue::Lerp PRESET_SEGUE{12, ease_in_out_sin};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 42;
  static constexpr bool ANIMATED_PROJECTION = false;

  static HS_HOT_FLASH_MEMBER Color4
  shade(const Vector &view, const typename FacetWave::FrameState &frame) {
    return FacetWave::RenderPipeline::shade(view, frame);
  }
  static constexpr Params initial_params() {
    Params value;
    value.source = {.pattern_freq = 6.3287f,
                    .speed = 0.04f,
                    .complexity = 1.704f,
                    .pattern_mix = 0.0f,
                    .secondary_rate = 0.8f,
                    .angle_rate = 0.027f};
    value.projection.pole_fade = 2.311f;
    value.projection.camera_wander = 1.0f;
    value.outer_warp.strength = -0.176f;
    value.outer_warp.speed = -0.00325f;
    value.outer_warp.frequency = 1.408f;
    value.outer_warp.field_angle = 2.2305307f;
    value.color.hue_shift_amount = 0.721f;
    value.color.palette_chroma = 1.0f;
    return value;
  }
  /**
   * @brief Params for the preset at @p index in PRESET_IDS.
   * @details `cup-hue` varies the colorizer alone: the cup palette mapping at
   * full hue rotation over a wider hue-noise field.
   */
  static constexpr Params preset_params(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
      value.color.hue_shift_amount = 1.0f;
      value.color.hue_noise_scale = 1.9717969f;
    }
    return value;
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(FacetWave)
