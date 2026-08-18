/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/look.h"

using FacetWaveParams =
    FixedLook::Params<FixedLook::GridSourceParams, FixedLook::WaveShearParams,
                      FixedLook::MirrorParams>;
using FacetWaveSpec =
    FixedLook::LookSpec<FixedLook::LookProjection::GNOMONIC_FOLDED,
                        Pullback::Lens::DodecahedralKaleidoscope,
                        FixedLook::LookTransfer::LINEAR,
                        FixedLook::LookCoverage::PROJECTION_SQUARED>;

template <int W, int H>
class FacetWave
    : public FixedLook::Look<W, H, FacetWave<W, H>, FacetWaveParams,
                             FacetWaveSpec, PaletteHarmony::TRIADIC,
                             FixedLook::HueMode::NOISE,
                             Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = FacetWaveParams;
  static constexpr std::string_view EFFECT_ID = "facet-wave";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "ee241e340b3f0133d2acc8d6f7b2f6c6926e7132127349035c1d2321323d4ae3";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "8c9a7d04c494ead43710bda0ad354acbf02e4d91279efc6e21469e66057060b0";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"wave-mirror"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = false;

  static HS_HOT_FLASH_MEMBER Color4
  shade(const Vector &view, const typename FacetWave::FrameState &frame) {
    return FacetWave::RenderPipeline::shade(view, frame);
  }
  static constexpr Params initial_params() {
    Params value;
    value.source = {6.3287f, 0.04f, 1.704f, 0.0f, 0.8f, 0.027f};
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
  static constexpr Params preset_params(size_t) { return initial_params(); }
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(FacetWave)
