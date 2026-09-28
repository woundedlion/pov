/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using KaleidoscopePentBrightParams =
    Pullback::Params<Pullback::LatticeSourceParams, Pullback::PolarParams,
                     Pullback::WaveShearParams>;
using KaleidoscopePentBrightSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                   Pullback::Lens::PentagonalPrismKaleidoscope,
                   Pullback::TransferKind::NONE,
                   Pullback::ProjectionCoverageMode::WEIGHT_SQUARED>;

/**
 * @brief A polar lattice folded through a pentagonal prism.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopePentBright
    : public Pullback::ComposedEffect<
          W, H, KaleidoscopePentBright<W, H>, KaleidoscopePentBrightParams,
          KaleidoscopePentBrightSpec, PaletteHarmony::ANALOGOUS,
          Pullback::HueMode::NOISE, Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = KaleidoscopePentBrightParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-pent-bright";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "27365967105ac1855d8ab3d5e00953436be958ec1238de5cadfebb42ae49e0f2";
  static constexpr std::string_view PRESET_BANK_DIGEST = "6e6135099d89c44bf87a5efd49ee5a7747132b6bdc370e3ebf3358fc41488512";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "polar-wave"
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
    value.color.hue_noise_scale = 2.0f;
    value.color.hue_noise_speed = 0.0f;
    value.color.hue_shift_amount = 0.268f;
    value.color.mapping_frequency = 1.0f;
    value.color.mapping_phase = -0.166f;
    value.color.palette_chroma = 1.0f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.color.phase_oscillation_depth = 0.0f;
    value.color.phase_oscillation_speed = 0.0f;
    value.color.opacity_high = 1.0f;
    value.color.opacity_low = 1.0f;
    value.projection.spin_rate = 0.0f;
    value.projection.wander = 1.0f;
    value.projection.singularity_fade = 2.273f;
    value.source.lattice_cell_scale = 0.7957747f;
    value.source.lattice_radius = 0.2907625f;
    value.source.lattice_shape_blend = 1.0f;
    value.source.lattice_softness = 0.3776084f;
    value.outer_warp.angular_phase = 0.0f;
    value.outer_warp.radial_phase = 0.0f;
    value.outer_warp.radial_scale = 1.0f;
    value.outer_warp.speed = 0.00034375f;
    value.inner_warp.field_angle = 0.0f;
    value.inner_warp.frequency = 1.0f;
    value.inner_warp.speed = 0.0009999999f;
    value.inner_warp.strength = 0.0f;
    return value;
  }
  // clang-format on
  // End generated params.
};
