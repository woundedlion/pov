/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using MermaidSkinParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::NoWarpParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::NoValueParams, Pullback::SurfaceNoiseParams>;
using MermaidSkinSpec =
    Pullback::Spec<Pullback::ProjectionKind::FOLDED_SINUSOIDAL, void,
                   Pullback::TransferKind::NONE,
                   Pullback::ProjectionCoverageMode::WEIGHT>;

/**
 * @brief A high-chroma folded grid rippling through sphere-space curl noise.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class MermaidSkin
    : public Pullback::ComposedEffect<
          W, H, MermaidSkin<W, H>, MermaidSkinParams, MermaidSkinSpec,
          PaletteHarmony::ANALOGOUS, Pullback::HueMode::NOISE,
          Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = MermaidSkinParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "mermaid-skin";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "804160646a2557fc6148dcfecf482a3ddb2c5f121b7c6fe1d5370957b4d70799";
  static constexpr std::string_view PRESET_BANK_DIGEST = "1275305c316679c959524e79ba42218d81a7b88b9ee3950b5ee38b327a31030f";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "mermaid-skin"
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
    value.color.hue_noise_scale = 0.14453939f;
    value.color.hue_noise_speed = -0.0000010416667f;
    value.color.hue_shift_amount = 1.5958333f;
    value.color.mapping_frequency = 5.3755207f;
    value.color.mapping_phase = 0.0f;
    value.color.palette_chroma = 1.0f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.color.phase_oscillation_depth = 0.0f;
    value.color.phase_oscillation_speed = 0.0f;
    value.color.opacity_high = 0.9661458f;
    value.color.opacity_low = 1.0f;
    value.projection.central_meridian = 0.0f;
    value.projection.spin_rate = 0.0f;
    value.projection.wander = 0.0f;
    value.source.angle_rate = 0.0f;
    value.source.complexity = 0.0f;
    value.source.secondary_rate = 0.0f;
    value.source.pattern_freq = 0.1f;
    value.source.pattern_mix = 0.0f;
    value.source.speed = 0.0f;
    value.surface.scale = 4.9144287f;
    value.surface.speed = -0.00021158854f;
    value.surface.strength = 0.5f;
    return value;
  }
  // clang-format on
  // End generated params.
};
