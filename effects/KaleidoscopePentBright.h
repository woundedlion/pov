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
                   Pullback::CoverageKind::PROJECTION_SQUARED>;

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
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-pent-bright";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "84c95558895ed0968fc1b25714fd08e7088656134d201800bf73c6abd3f4ddc8";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "50720ad712274ad1d2a4e56ee85af0f9dd3795f9e224dc3996bb568e7c53b1b3";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"polar-wave"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;

  static constexpr Params initial_params() {
    Params value;
    value.source.lattice_cell_scale = 0.774140596f;
    value.source.lattice_shape_blend = 1.0f;
    value.source.lattice_softness = 0.377608389f;
    value.source.lattice_radius = 0.290762514f;
    value.projection.singularity_fade = 2.273f;
    value.projection.wander = 1.0f;
    value.projection.camera_wander = 1.0f;
    value.outer_warp.speed = 0.000343749998f;
    value.inner_warp.speed = 0.000999999931f;
    value.color.hue_shift_amount = 0.268000007f;
    value.color.hue_noise_scale = 2.0f;
    value.color.palette_chroma = 1.0f;
    value.color.mapping_phase = -0.165999994f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    return value;
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(KaleidoscopePentBright)
