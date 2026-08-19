/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using EquatorGridParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::NoWarpParams,
                     Pullback::MirrorParams>;
using EquatorGridSpec =
    Pullback::Spec<Pullback::ProjectionKind::EQUIRECTANGULAR,
                   Pullback::Lens::DodecahedralKaleidoscope,
                   Pullback::TransferKind::LINEAR,
                   Pullback::CoverageKind::PROJECTION_SQUARED>;

template <int W, int H>
class EquatorGrid
    : public Pullback::ComposedEffect<
          W, H, EquatorGrid<W, H>, EquatorGridParams, EquatorGridSpec,
          PaletteHarmony::ANALOGOUS, Pullback::HueMode::NOISE,
          Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = EquatorGridParams;
  static constexpr std::string_view EFFECT_ID = "equator-grid";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "4ff4b908978de7ecc20df3bb1ea72cc6c51cee39e663cd3c995d5503463a19f2";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "f3a6aa9fc8053e9cfa84a487a51a8801ac9ea484a3b1091464fa304a3e727a58";
  static constexpr std::array<std::string_view, 3> PRESET_IDS{
      "double-map", "open-grid", "fine-grid"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr bool USES_CENTRAL_MERIDIAN = true;

  static constexpr Params initial_params() {
    Params value;
    value.source = {3.9407f, 0.0f, 3.0f, 1.0f, 0.8f, 0.0269999988f};
    value.projection.pole_fade = 2.14f;
    value.projection.wander = 0.165f;
    value.projection.camera_wander = 1.0f;
    value.inner_warp.speed = 0.00013f;
    value.inner_warp.cell_y = 0.997703135f;
    value.color.hue_shift_amount = 0.366f;
    value.color.hue_noise_scale = 1.47215629f;
    value.color.palette_chroma = 1.0f;
    value.color.mapping_frequency = 2.0f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    return value;
  }
  static constexpr Params preset_params(size_t index) {
    Params value = initial_params();
    if (index == 1)
      value.color.mapping_frequency = 1.0f;
    else if (index == 2) {
      value.source.pattern_freq = 0.3985f;
      value.inner_warp.speed = 0.00058f;
      value.inner_warp.cell_y = 0.901890635f;
      value.color.mapping_frequency = 21.212f;
    }
    return value;
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(EquatorGrid)
