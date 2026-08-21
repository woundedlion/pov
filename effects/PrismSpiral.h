/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using PrismSpiralParams =
    Pullback::Params<Pullback::SpiralSourceParams, Pullback::NoWarpParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::NoValueParams, Pullback::DirectSurfaceParams>;
using PrismSpiralSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                   Pullback::Lens::HexagonalPrismKaleidoscope,
                   Pullback::TransferKind::NONE,
                   Pullback::CoverageKind::PROJECTION_SQUARED>;

template <int W, int H>
class PrismSpiral
    : public Pullback::ComposedEffect<
          W, H, PrismSpiral<W, H>, PrismSpiralParams, PrismSpiralSpec,
          PaletteHarmony::TRIADIC, Pullback::HueMode::PATH_LENGTH,
          Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = PrismSpiralParams;
  static constexpr std::string_view EFFECT_ID = "prism-spiral";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "547e7f06fc2c0d02605dcba58e5235706d1447c77f51ae2623314e967383b251";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "fb3772baa24f81106bd3eac928ae9e051ebf9703b4c6530044f085db5f2ec5c9";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"prism-spiral",
                                                              "prism-spiral-2"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  // The 6 s Phantasm slot is 96 frames at 16 fps and the effect is rebuilt each
  // visit: 2 dwells plus 1 crossfade must fit, or the later presets never
  // render.
  static constexpr Segue::Lerp PRESET_SEGUE{12, ease_in_out_sin};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 42;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr int32_t SURFACE_NOISE_SEED = Pullback::EFFECT_NOISE_SEED;

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
    Params value = initial_params();
    if (index == 1) {
      value.surface.scale = 3.6627343f;
      value.color.mapping_frequency = 1.2f;
    }
    return value;
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(PrismSpiral)
