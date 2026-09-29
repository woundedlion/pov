/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using KaleidoscopeMandalaParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::WaveShearParams,
                     Pullback::MirrorParams>;
using KaleidoscopeMandalaSpec =
    Pullback::Spec<Pullback::ProjectionKind::GNOMONIC_FOLDED,
                   Pullback::Lens::DodecahedralKaleidoscope,
                   Pullback::TransferKind::NONE,
                   Pullback::ProjectionCoverageMode::WEIGHT_SQUARED>;

/**
 * @brief Folded-gnomonic wave field reflected through a dodecahedral lens.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopeMandala
    : public Pullback::ComposedEffect<
          W, H, KaleidoscopeMandala<W, H>, KaleidoscopeMandalaParams,
          KaleidoscopeMandalaSpec, PaletteHarmony::TRIADIC,
          Pullback::HueMode::NOISE, Pullback::Color::BrightnessEnvelope::NONE,
          /*AnimatedProjection=*/false> {

public:
  using Params = KaleidoscopeMandalaParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-mandala";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "4bc222a0e81f0aa4a6ed06db0dddcfca84486072d2bb78d6665b8370edec9719";
  static constexpr std::string_view PRESET_BANK_DIGEST = "7e10bc72b93c0671877e68f54c36bcb6177b5038f3e12c7af5aebed8c4ef1f56";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "wave-mirror",
      "cup-hue"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

  // Hot entry trampoline; RenderPipeline::shade retains its cold placement.
  static HS_HOT_FLASH_MEMBER Color4
  shade(const math::Vector &view,
        const typename KaleidoscopeMandala::Frame &frame) {
    return KaleidoscopeMandala::RenderPipeline::shade(view, frame);
  }
  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.projection.camera_wander = 1.0f;
    value.color.hue_noise_scale = 1.0f;
    value.color.hue_noise_speed = 0.0f;
    value.color.hue_shift_amount = 0.721f;
    value.color.mapping_frequency = 1.0f;
    value.color.mapping_phase = 0.0f;
    value.color.palette_chroma = 1.0f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.color.phase_oscillation_depth = 0.0f;
    value.color.phase_oscillation_speed = 0.0f;
    value.color.opacity_high = 1.0f;
    value.color.opacity_low = 1.0f;
    value.projection.singularity_fade = 2.311f;
    value.source.angle_rate = 0.027f;
    value.source.complexity = 1.704f;
    value.source.secondary_rate = 0.8f;
    value.source.pattern_freq = 6.3287f;
    value.source.pattern_mix = 0.0f;
    value.source.speed = 0.04f;
    value.outer_warp.field_angle = 2.2305307f;
    value.outer_warp.frequency = 1.408f;
    value.outer_warp.speed = -0.00325f;
    value.outer_warp.strength = -0.176f;
    value.inner_warp.cell_x = 1.0f;
    value.inner_warp.cell_y = 1.0f;
    value.inner_warp.offset_x = 0.0f;
    value.inner_warp.offset_y = 0.0f;
    value.inner_warp.rotation = 0.0f;
    value.inner_warp.speed = 0.0f;
    return value;
  }

  /** @brief The preset at index in PRESET_IDS and how it departs. */
  HS_COLD_MEMBER static constexpr PresetEntry<Params> preset(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.color.hue_noise_scale = 1.9717969f;
      value.color.hue_shift_amount = 1.0f;
      value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    }
    return {value, Segue::Preset::Lerp{480, math::ease_in_out_sin}};
  }
  // clang-format on
  // End generated params.
};
