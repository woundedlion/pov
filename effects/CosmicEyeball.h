/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using CosmicEyeballParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::MirrorParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::EdgeValueParams>;
using CosmicEyeballSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                   Pullback::Lens::Glitch, Pullback::TransferKind::NONE,
                   Pullback::CoverageKind::EDGE_FADE>;

template <int W, int H>
class CosmicEyeball
    : public Pullback::ComposedEffect<
          W, H, CosmicEyeball<W, H>, CosmicEyeballParams, CosmicEyeballSpec,
          PaletteHarmony::TRIADIC, Pullback::HueMode::PATH_LENGTH,
          Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = CosmicEyeballParams;
  static constexpr std::string_view EFFECT_ID = "cosmic-eyeball";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "56705b20f6b96f3c3a7436dabc9453881ff6ce0978445f60533ace25cdcb77f4";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "0b5ac4830de902361a5b7066779b8f222e2589d6364488df5829423c1185ab01";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"mirrored-grid"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = false;

  static constexpr Params initial_params() {
    Params value;
    value.source = {.pattern_freq = 2.5477f,
                    .speed = 0.235f,
                    .complexity = 1.854f,
                    .pattern_mix = 0.0f,
                    .secondary_rate = 1.0f,
                    .angle_rate = 0.0f};
    value.projection.singularity_fade = 1.4f;
    value.projection.camera_wander = 1.0f;
    value.outer_warp.rotation = 0.295309722f;
    value.outer_warp.cell_x = 5.381125f;
    value.outer_warp.offset_x = 1.344f;
    value.outer_warp.offset_y = -1.456f;
    value.value.edge_width = 0.5f;
    value.color.hue_shift_amount = 2.048f;
    value.color.palette_chroma = 0.292f;
    return value;
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(CosmicEyeball)
