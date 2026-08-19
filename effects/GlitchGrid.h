/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using GlitchGridParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::MirrorParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::EdgeValueParams>;
using GlitchGridSpec =
    Pullback::Spec<Pullback::ProjectionKind::GNOMONIC_FOLDED,
                   Pullback::Lens::Glitch, Pullback::TransferKind::LINEAR,
                   Pullback::CoverageKind::EDGE_FADE>;

template <int W, int H>
class GlitchGrid : public Pullback::ComposedEffect<
                       W, H, GlitchGrid<W, H>, GlitchGridParams, GlitchGridSpec,
                       PaletteHarmony::TRIADIC, Pullback::HueMode::NOISE,
                       Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = GlitchGridParams;
  static constexpr std::string_view EFFECT_ID = "glitch-grid";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "cf6ad02cd20fffba97692e8e3a8028760e3985351e02b4b98a79b8670c41b7c5";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "d22beff986c821fc15e34d7ddb4439f7f6df10591f9e37f36da5e507bfecfb02";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"folded-glitch"};
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
    return value;
  }
  static constexpr Params preset_params(size_t) { return initial_params(); }
};

#include "core/control/registry.h"
REGISTER_EFFECT(GlitchGrid)
