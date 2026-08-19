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
                   Pullback::TransferKind::LINEAR,
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
      "bc994ccf7cfae7cc3b2b631d106f4f03f9c79819f7d55c41f232c2d5c946a5b7";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "acbd113115722d2116d877977f00d2ad56a1b38e94f53b36df817be8357603bb";
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

#include "core/control/registry.h"
REGISTER_EFFECT(FacetWave)
