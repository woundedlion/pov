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
                   Pullback::ProjectionCoverageMode::WEIGHT_SQUARED>;

/**
 * @brief A continuously animated Mobius lens over a mirrored twin wave.
 * @details Presets store an authored lens snapshot. Once the timeline steps,
 * the circular animation drives b on the unit circle; authored b is not a
 * rendered-frame promise.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class MobiusGrid
    : public Pullback::ComposedEffect<
          W, H, MobiusGrid<W, H>, MobiusGridParams, MobiusGridSpec,
          PaletteHarmony::COMPLEMENTARY, Pullback::HueMode::PATH_LENGTH,
          Pullback::Color::BrightnessEnvelope::CUP> {

public:
  using Params = MobiusGridParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "mobius-grid";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "0fc45c8e12712494bc7b6160a103e29119ea382defece2ba21488f87fe8b3cbf";
  static constexpr std::string_view PRESET_BANK_DIGEST = "eb6c7f90002c374d45ecf50f0e9ca9070fa7ce80c9d79119883697659bbba82c";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "mobius-grid",
      "mobius-grid-2"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr bool ANIMATED_MOBIUS = true;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.projection.camera_wander = 1.0f;
    value.color.brightness_bottom = 0.0f;
    value.color.brightness_top = 1.0f;
    value.color.hue_shift_amount = 0.312f;
    value.color.mapping_frequency = 1.0f;
    value.color.mapping_phase = 0.0f;
    value.color.palette_chroma = 0.398f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.color.phase_oscillation_depth = 0.0f;
    value.color.phase_oscillation_speed = 0.0f;
    value.color.opacity_high = 1.0f;
    value.color.opacity_low = 1.0f;
    value.lens.mobius.a.im = 0.304f;
    value.lens.mobius.a.re = -1.072f;
    value.lens.mobius.b.im = 0.0f;
    value.lens.mobius.b.re = 0.416f;
    value.lens.mobius.c.im = 0.0f;
    value.lens.mobius.c.re = 0.0f;
    value.lens.mobius.d.im = 0.0f;
    value.lens.mobius.d.re = 0.70710677f;
    value.projection.spin_rate = 0.0f;
    value.projection.wander = 1.0f;
    value.projection.singularity_fade = 2.102f;
    value.source.angle_rate = 0.027f;
    value.source.secondary_rate = 0.8f;
    value.source.pattern_freq = 10.158f;
    value.source.speed = 0.245f;
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
      value.inner_warp.cell_x = 0.2791094f;
      value.inner_warp.cell_y = 6.810328f;
      value.inner_warp.speed = 0.005875f;
    }
    return {value, Segue::Preset::Lerp{480, math::ease_in_out_sin}};
  }
  // clang-format on
  // End generated params.

  HS_COLD_MEMBER void after_composed_init() {
    this->start_mobius_animation(1.0f, 160);
  }
};
