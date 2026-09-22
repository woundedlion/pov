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
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-mandala";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "4bc222a0e81f0aa4a6ed06db0dddcfca84486072d2bb78d6665b8370edec9719";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "7e10bc72b93c0671877e68f54c36bcb6177b5038f3e12c7af5aebed8c4ef1f56";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"wave-mirror",
                                                              "cup-hue"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;

  // Hot section: the out-of-line pipeline body compiles for speed.
  static HS_HOT_FLASH_MEMBER Color4
  shade(const math::Vector &view,
        const typename KaleidoscopeMandala::Frame &frame) {
    return KaleidoscopeMandala::RenderPipeline::shade(view, frame);
  }
  static constexpr Params initial_params() {
    Params value;
    value.source = {.pattern_freq = 6.3287f,
                    .speed = 0.04f,
                    .complexity = 1.704f,
                    .pattern_mix = 0.0f,
                    .secondary_rate = 0.8f,
                    .angle_rate = 0.027f};
    value.projection.singularity_fade = 2.311f;
    value.projection.camera_wander = 1.0f;
    value.outer_warp.strength = -0.176f;
    value.outer_warp.speed = -0.00325f;
    value.outer_warp.frequency = 1.408f;
    value.outer_warp.field_angle = 2.2305307f;
    value.color.hue_shift_amount = 0.721f;
    value.color.palette_chroma = 1.0f;
    return value;
  }
  /**
   * @brief Params for the preset at @p index in PRESET_IDS.
   * @details `cup-hue` varies the colorizer alone: the cup palette mapping at
   * full hue rotation over a finer hue-noise field.
   */
  static constexpr Params preset_params(size_t index) {
    static_assert(PRESET_IDS.size() == 2,
                  "a new preset id needs a branch below");
    Params value = initial_params();
    if (index == 1) {
      value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
      value.color.hue_shift_amount = 1.0f;
      value.color.hue_noise_scale = 1.9717969f;
    }
    return value;
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(KaleidoscopeMandala)
