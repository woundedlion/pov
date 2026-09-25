/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using KaleidoscopeStainedGlassParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::VectorNoiseParams,
                     Pullback::MirrorParams>;
using KaleidoscopeStainedGlassSpec =
    Pullback::Spec<Pullback::ProjectionKind::GNOMONIC_FOLDED,
                   Pullback::Lens::DodecahedralKaleidoscope,
                   Pullback::TransferKind::NONE,
                   Pullback::ProjectionCoverageMode::WEIGHT_SQUARED>;

/**
 * @brief A vector-noise grid refracted across dodecahedral facets.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopeStainedGlass
    : public Pullback::ComposedEffect<
          W, H, KaleidoscopeStainedGlass<W, H>, KaleidoscopeStainedGlassParams,
          KaleidoscopeStainedGlassSpec, PaletteHarmony::TRIADIC,
          Pullback::HueMode::NOISE, Pullback::Color::BrightnessEnvelope::CUP,
          /*AnimatedProjection=*/false> {

public:
  using Params = KaleidoscopeStainedGlassParams;
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-stained-glass";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "dec7ed0755dbb3565429dcb17644a21d2b285cf2318797855d2e091dcadce7ac";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "0fccbe4ac4eb7acdba557f03c8ddadb4f974cb7400470e5e67951d2e7f2f7713";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"vector-mirror"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;

  // Cold entry trampoline and pipeline body.
  static HS_FLASH_MEMBER Color4
  shade(const math::Vector &view,
        const typename KaleidoscopeStainedGlass::Frame &frame) {
    return KaleidoscopeStainedGlass::RenderPipeline::shade(view, frame);
  }
  static constexpr Params initial_params() {
    Params value;
    value.source = {.pattern_freq = 4.9755f,
                    .speed = 0.04f,
                    .complexity = 1.704f,
                    .pattern_mix = 0.0f,
                    .secondary_rate = 0.8f,
                    .angle_rate = 0.027f};
    value.projection.singularity_fade = 2.311f;
    value.projection.camera_wander = 1.0f;
    value.outer_warp.strength = 0.138f;
    value.outer_warp.speed = -0.00005f;
    value.inner_warp.speed = 0.00327999983f;
    value.color.hue_shift_amount = 0.721f;
    value.color.palette_chroma = 1.0f;
    value.color.brightness_bottom = 0.345f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    return value;
  }
};
