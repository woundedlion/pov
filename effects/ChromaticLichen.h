/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using ChromaticLichenParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::NoWarpParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::NoValueParams, Pullback::SurfaceNoiseParams>;
using ChromaticLichenSpec =
    Pullback::Spec<Pullback::ProjectionKind::GNOMONIC_FOLDED,
                   Pullback::Lens::Glitch, Pullback::TransferKind::NONE,
                   Pullback::ProjectionCoverageMode::WEIGHT>;

/**
 * @brief A glitch-folded gnomonic grid displaced by sphere-space curl noise.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class ChromaticLichen
    : public Pullback::ComposedEffect<
          W, H, ChromaticLichen<W, H>, ChromaticLichenParams,
          ChromaticLichenSpec, PaletteHarmony::ANALOGOUS,
          Pullback::HueMode::NOISE, Pullback::Color::BrightnessEnvelope::NONE,
          /*AnimatedProjection=*/true, Pullback::SurfacePlacement::AFTER_LENS> {

public:
  using Params = ChromaticLichenParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "chromatic-lichen";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "c6d7eb734f0e5d0b157c6a2d171657d7a2c6fb6daabaca6918a4dd8c6e700568";
  static constexpr std::string_view PRESET_BANK_DIGEST = "8a6c7cda1ac938cd2f39211a251841a586ce134ee750683b3c1a9d2d2d524f63";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "chromatic-lichen"
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
    value.color.opacity_high = 1.0f;
    value.color.opacity_low = 1.0f;
    value.projection.spin_rate = 0.0f;
    value.projection.wander = 1.0f;
    value.projection.singularity_fade = 1.0f;
    value.source.angle_rate = 0.0f;
    value.source.complexity = 0.0f;
    value.source.secondary_rate = 0.0f;
    value.source.pattern_freq = 0.1f;
    value.source.pattern_mix = 0.0f;
    value.source.speed = 0.0f;
    value.surface.scale = 2.3483887f;
    value.surface.speed = -0.00021158854f;
    value.surface.strength = 0.5f;
    return value;
  }
  // clang-format on
  // End generated params.
};
