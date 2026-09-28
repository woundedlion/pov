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
                   Pullback::Lens::Kaleidoscope, Pullback::TransferKind::NONE,
                   Pullback::ProjectionCoverageMode::EDGE_FADE>;

/**
 * @brief A broad folded grid with slow mirrored drift.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class AlienOcean
    : public Pullback::ComposedEffect<W, H, AlienOcean<W, H>, AlienOceanParams,
                                      AlienOceanSpec, PaletteHarmony::TRIADIC,
                                      Pullback::HueMode::NOISE,
                                      Pullback::Color::BrightnessEnvelope::NONE,
                                      /*AnimatedProjection=*/false> {

public:
  using Params = AlienOceanParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "alien-ocean";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "1bb8a0674b3bee614ef26ab004922cb8022f5c3984884a5bd3c51b7b5beb6490";
  static constexpr std::string_view PRESET_BANK_DIGEST = "a9070db6d46b6ed7aaaf3341a43cf3cf2486948e78b55384f5d48877d92e5061";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "folded-grid"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.projection.camera_wander = 1.0f;
    value.color.hue_noise_scale = 2.2033439f;
    value.color.hue_noise_speed = 0.0f;
    value.color.hue_shift_amount = 0.424f;
    value.color.mapping_frequency = 1.0f;
    value.color.mapping_phase = 0.0f;
    value.color.palette_chroma = 0.4f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.color.phase_oscillation_depth = 0.0f;
    value.color.phase_oscillation_speed = 0.0f;
    value.color.opacity_high = 1.0f;
    value.color.opacity_low = 1.0f;
    value.projection.singularity_fade = 1.4f;
    value.source.angle_rate = 0.0f;
    value.source.complexity = 0.0f;
    value.source.secondary_rate = 1.0f;
    value.value.edge_width = 0.5f;
    value.source.pattern_freq = 3.565f;
    value.source.pattern_mix = 1.0f;
    value.source.speed = 0.235f;
    value.outer_warp.cell_x = 5.381125f;
    value.outer_warp.cell_y = 1.0f;
    value.outer_warp.offset_x = 1.344f;
    value.outer_warp.offset_y = -1.456f;
    value.outer_warp.rotation = 0.29530972f;
    value.outer_warp.speed = 0.0f;
    return value;
  }
  // clang-format on
  // End generated params.
};
