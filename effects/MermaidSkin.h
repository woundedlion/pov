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
  static constexpr std::string_view EFFECT_ID = "mermaid-skin";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "a6ced8928ed3f35ba178dd5d6a538309af93bd86d9abc967380cf21b42ac7d5a";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "1275305c316679c959524e79ba42218d81a7b88b9ee3950b5ee38b327a31030f";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"mermaid-skin"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;

  static constexpr Params initial_params() {
    Params value;
    value.source.pattern_freq = 0.1f;
    value.projection.camera_wander = 1.0f;
    value.surface = {
        .scale = 4.91442871f, .strength = 0.5f, .speed = -0.000211588544f};
    value.color.hue_shift_amount = 1.5958333f;
    value.color.hue_noise_scale = 0.144539386f;
    value.color.hue_noise_speed = -0.00000104166668f;
    value.color.palette_chroma = 1.0f;
    value.color.mapping_frequency = 5.37552071f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.color.opacity_high = 0.966145813f;
    return value;
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(MermaidSkin)
