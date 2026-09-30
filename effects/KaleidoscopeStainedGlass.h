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
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-stained-glass";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "dec7ed0755dbb3565429dcb17644a21d2b285cf2318797855d2e091dcadce7ac";
  static constexpr std::string_view PRESET_BANK_DIGEST = "0fccbe4ac4eb7acdba557f03c8ddadb4f974cb7400470e5e67951d2e7f2f7713";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "vector-mirror"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

  // Cold entry trampoline; the pipeline body uses hot flash.
  static HS_HOT_FLASH_MEMBER Color4
  shade(const math::Vector &view,
        const typename KaleidoscopeStainedGlass::Frame &frame) {
    return KaleidoscopeStainedGlass::RenderPipeline::shade(view, frame);
  }
  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.projection.camera_wander = 1.0f;
    value.color.brightness_bottom = 0.345f;
    value.color.brightness_top = 1.0f;
    value.color.hue_noise_scale = 1.0f;
    value.color.hue_noise_speed = 0.0f;
    value.color.hue_shift_amount = 0.721f;
    value.color.mapping_frequency = 1.0f;
    value.color.mapping_phase = 0.0f;
    value.color.palette_chroma = 1.0f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.color.phase_oscillation_depth = 0.0f;
    value.color.phase_oscillation_speed = 0.0f;
    value.color.opacity_high = 1.0f;
    value.color.opacity_low = 1.0f;
    value.projection.singularity_fade = 2.311f;
    value.source.angle_rate = 0.027f;
    value.source.complexity = 1.704f;
    value.source.secondary_rate = 0.8f;
    value.source.pattern_freq = 4.9755f;
    value.source.pattern_mix = 0.0f;
    value.source.speed = 0.04f;
    value.outer_warp.scale = 1.0f;
    value.outer_warp.speed = -0.00005f;
    value.outer_warp.strength = 0.138f;
    value.outer_warp.vector_angle = 0.0f;
    value.inner_warp.cell_x = 1.0f;
    value.inner_warp.cell_y = 1.0f;
    value.inner_warp.offset_x = 0.0f;
    value.inner_warp.offset_y = 0.0f;
    value.inner_warp.rotation = 0.0f;
    value.inner_warp.speed = 0.0032799998f;
    return value;
  }
  // clang-format on
  // End generated params.
};
