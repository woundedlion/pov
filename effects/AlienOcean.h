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
      "515b0a720e37683cca416e5d67236911ea0285f7e9c4dd027a82188280be6b1f";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "bb79b346ffc5deb0cbb34479e72013f181b636f22f7dba3d989d7d327125961c";
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
};

#include "core/control/registry.h"
REGISTER_EFFECT(AlienOcean)
