/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using MobiusGridParams =
    Pullback::Params<Pullback::TwinWaveSourceParams, Pullback::NoWarpParams,
                     Pullback::MirrorParams, Pullback::MobiusLensParams>;
using MobiusGridSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC, void,
                   Pullback::TransferKind::NONE,
                   Pullback::CoverageKind::PROJECTION_SQUARED>;

template <int W, int H>
class MobiusGrid
    : public Pullback::ComposedEffect<
          W, H, MobiusGrid<W, H>, MobiusGridParams, MobiusGridSpec,
          PaletteHarmony::COMPLEMENTARY, Pullback::HueMode::PATH_LENGTH,
          Pullback::Color::BrightnessEnvelope::CUP> {

public:
  using Params = MobiusGridParams;
  static constexpr std::string_view EFFECT_ID = "mobius-grid";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "0ed547969b9d9133d42de9a686da7ab8d6f530d77a0a4a48a9398ff82ecf24f9";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "990e359263891383011d53fcd0b3f46d6b3229bc0482409fd9cf762231322b2a";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"mobius-grid",
                                                              "mobius-grid-2"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  // The 6 s Phantasm slot is 96 frames at 16 fps and the effect is rebuilt each
  // visit: 2 dwells plus 1 crossfade must fit, or the later presets never
  // render.
  static constexpr Segue::Lerp PRESET_SEGUE{12, ease_in_out_sin};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 42;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr bool ANIMATED_MOBIUS = true;

  static constexpr Params initial_params() {
    Params value;
    value.source = {.pattern_freq = 10.158f,
                    .speed = 0.245f,
                    .secondary_rate = 0.8f,
                    .angle_rate = 0.027f};
    value.projection.pole_fade = 2.102f;
    value.projection.wander = 1.0f;
    value.projection.camera_wander = 1.0f;
    value.lens.mobius = {-1.072f, 0.304f, 0.416f,      0.0f,
                         0.0f,    0.0f,   0.70710677f, 0.0f};
    value.color.hue_shift_amount = 0.312f;
    value.color.palette_chroma = 0.398f;
    return value;
  }
  static constexpr Params preset_params(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.inner_warp.speed = 0.005875f;
      value.inner_warp.cell_x = 0.2791094f;
      value.inner_warp.cell_y = 6.810328f;
    }
    return value;
  }

public:
  HS_COLD_MEMBER void after_composed_init() {
    this->start_mobius_animation(1.0f, 160);
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(MobiusGrid)
