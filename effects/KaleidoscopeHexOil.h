/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using KaleidoscopeHexOilParams =
    Pullback::Params<Pullback::SpiralSourceParams, Pullback::NoWarpParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::NoValueParams, Pullback::DirectSurfaceParams>;
using KaleidoscopeHexOilSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                   Pullback::Lens::HexagonalPrismKaleidoscope,
                   Pullback::TransferKind::NONE,
                   Pullback::ProjectionCoverageMode::WEIGHT_SQUARED>;

/**
 * @brief A rotating spiral folded through a hexagonal prism kaleidoscope.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopeHexOil
    : public Pullback::ComposedEffect<
          W, H, KaleidoscopeHexOil<W, H>, KaleidoscopeHexOilParams,
          KaleidoscopeHexOilSpec, PaletteHarmony::TRIADIC,
          Pullback::HueMode::PATH_LENGTH,
          Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = KaleidoscopeHexOilParams;
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-hex-oil";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "54591ed7a8c2ee0325fb55adf2f8bad3d22c567e5ce9e0a4bfaeaaf99f7692ee";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "302f86fbdcd6122935cbc7e0c80bcff24f441401e3d44d39ec096e02849fcd18";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "kaleidoscope-hex-oil", "kaleidoscope-hex-oil-2"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;

  static constexpr Params initial_params() {
    Params value;
    value.source.pattern_freq = 5.5327f;
    value.source.speed = 0.125f;
    value.source.angle_rate = 0.03f;
    value.projection.singularity_fade = 1.627f;
    value.projection.wander = 1.0f;
    value.projection.camera_wander = 1.0f;
    value.surface.scale = 5.740422f;
    value.surface.strength = 0.4185f;
    value.color.hue_shift_amount = -2.216f;
    value.color.palette_chroma = 0.78f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    return value;
  }
  static constexpr Params preset_params(size_t index) {
    static_assert(PRESET_IDS.size() == 2,
                  "a new preset id needs a branch below");
    Params value = initial_params();
    if (index == 1) {
      value.surface.scale = 3.6627343f;
      value.color.mapping_frequency = 1.2f;
    }
    return value;
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(KaleidoscopeHexOil)
