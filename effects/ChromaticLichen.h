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
          /*AnimatedProjection=*/true, /*HasOuterNoise=*/false,
          /*HasSourceNoise=*/false, Pullback::SurfacePlacement::AFTER_LENS> {

public:
  using Params = ChromaticLichenParams;
  static constexpr std::string_view EFFECT_ID = "chromatic-lichen";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "1412a19712099c9c014a396e50797dfffbf77d3b41493b74f4137ba442468b6d";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "8a6c7cda1ac938cd2f39211a251841a586ce134ee750683b3c1a9d2d2d524f63";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "chromatic-lichen"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;

  static constexpr Params initial_params() {
    Params value;
    value.source.pattern_freq = 0.1f;
    value.projection.wander = 1.0f;
    value.projection.camera_wander = 1.0f;
    value.surface = {
        .scale = 2.34838867f, .strength = 0.5f, .speed = -0.000211588544f};
    value.color.hue_shift_amount = 1.5958333f;
    value.color.hue_noise_scale = 0.144539386f;
    value.color.hue_noise_speed = -0.00000104166668f;
    value.color.palette_chroma = 1.0f;
    value.color.mapping_frequency = 5.37552071f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    return value;
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(ChromaticLichen)
