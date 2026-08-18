/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using AlienOceanParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::MirrorParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::EdgeValueParams>;
using AlienOceanSpec =
    Pullback::Spec<Pullback::ProjectionKind::GNOMONIC_FOLDED,
                   Pullback::Lens::Kaleidoscope, Pullback::TransferKind::LINEAR,
                   Pullback::CoverageKind::EDGE_FADE>;

template <int W, int H>
class AlienOcean : public Pullback::ComposedEffect<
                       W, H, AlienOcean<W, H>, AlienOceanParams, AlienOceanSpec,
                       PaletteHarmony::TRIADIC, Pullback::HueMode::NOISE,
                       Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = AlienOceanParams;
  static constexpr std::string_view EFFECT_ID = "alien-ocean";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "b6c0f04a02550a66af8c0f140822ddeccfbf110aa3c57bc6631397d0f7672a8d";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "85c8b3961f3d1de228ff43e90f9f0e957f453db2b30f2cd98222b305ab0e77f6";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"folded-grid"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = false;

  static constexpr Params initial_params() {
    Params value;
    value.source = {3.565f, 0.235f, 0.0f, 1.0f, 1.0f, 0.0f};
    value.projection.pole_fade = 1.4f;
    value.projection.camera_wander = 1.0f;
    value.outer_warp.rotation = 0.29530972f;
    value.outer_warp.cell_x = 5.381125f;
    value.outer_warp.offset_x = 1.344f;
    value.outer_warp.offset_y = -1.456f;
    value.value.edge_width = 0.5f;
    value.color.hue_shift_amount = 0.424f;
    value.color.hue_noise_scale = 2.2033439f;
    value.color.palette_chroma = 0.4f;
    return value;
  }

  static constexpr Params preset_params(size_t) { return initial_params(); }
};

#include "core/engine/control/registry.h"
REGISTER_EFFECT(AlienOcean)
