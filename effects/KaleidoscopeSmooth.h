/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

namespace hs_test {
namespace kaleidoscope_smooth_tests {
struct KaleidoscopeSmoothWhiteBox;
} // namespace kaleidoscope_smooth_tests
} // namespace hs_test

using KaleidoscopeSmoothParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::NoWarpParams,
                     Pullback::MirrorParams>;
using KaleidoscopeSmoothSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                   Pullback::Lens::DodecahedralKaleidoscope,
                   Pullback::TransferKind::NONE,
                   Pullback::ProjectionCoverageMode::WEIGHT_SQUARED>;

/**
 * @brief Mirrored grids folded through a dodecahedral stereographic lens.
 * @details Supplies the render pipeline and preset bank; Pullback::ComposedEffect
 * supplies parameter registration, preset choreography and the palette,
 * camera-walk and noise clocks. The dodecahedral fold is the Lens stage and
 * the mirror tiling is the inner planar warp; the surface stage carries no
 * displacement.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopeSmooth
    : public Pullback::ComposedEffect<
          W, H, KaleidoscopeSmooth<W, H>, KaleidoscopeSmoothParams,
          KaleidoscopeSmoothSpec, PaletteHarmony::ANALOGOUS,
          Pullback::HueMode::NOISE, Pullback::Color::BrightnessEnvelope::NONE> {
  friend struct ::hs_test::kaleidoscope_smooth_tests::
      KaleidoscopeSmoothWhiteBox;

public:
  using Params = KaleidoscopeSmoothParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-smooth";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "5e5b60bf00084446b125be2ab46e77319bf035dbcf7583a2b43753ad6009a579";
  static constexpr std::string_view PRESET_BANK_DIGEST = "6ce5e7188843b4f091b9c3cc2689798e3f2e77420c9f3d688f5b2ca50485ec4a";
  static constexpr std::array<std::string_view, 4> PRESET_IDS{
      "coupled-grid",
      "direct-grid",
      "double-map",
      "stretched-grid"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 3;

  /// Params the effect starts on, and the base every preset varies from.
  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.projection.camera_wander = 1.0f;
    value.color.hue_noise_scale = 1.4721563f;
    value.color.hue_noise_speed = 0.0f;
    value.color.hue_shift_amount = 0.366f;
    value.color.mapping_frequency = 1.0f;
    value.color.mapping_phase = 0.0f;
    value.color.palette_chroma = 1.0f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.color.phase_oscillation_depth = 0.0f;
    value.color.phase_oscillation_speed = 0.0f;
    value.color.opacity_high = 1.0f;
    value.color.opacity_low = 1.0f;
    value.projection.spin_rate = 0.0f;
    value.projection.wander = 1.0f;
    value.projection.singularity_fade = 3.432f;
    value.source.angle_rate = 0.026999999f;
    value.source.complexity = 0.513f;
    value.source.secondary_rate = 0.8f;
    value.source.pattern_freq = 2.8263f;
    value.source.pattern_mix = 0.0f;
    value.source.speed = 0.0f;
    value.inner_warp.cell_x = 1.0f;
    value.inner_warp.cell_y = 0.99770314f;
    value.inner_warp.offset_x = 0.0f;
    value.inner_warp.offset_y = 0.0f;
    value.inner_warp.rotation = 0.0f;
    value.inner_warp.speed = 0.00013f;
    return value;
  }

  /** @brief Params for the preset at index in PRESET_IDS. */
  static constexpr Params preset_params(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.source.complexity = 3.0f;
      value.source.pattern_mix = 1.0f;
    }
    if (index == 2) {
      value.color.mapping_frequency = 2.0f;
      value.projection.wander = 0.165f;
      value.source.complexity = 3.0f;
      value.source.pattern_freq = 3.9407f;
      value.source.pattern_mix = 1.0f;
    }
    if (index == 3) {
      value.color.mapping_frequency = 1.558f;
      value.projection.wander = 0.165f;
      value.source.complexity = 3.0f;
      value.source.pattern_freq = 2.9059f;
      value.source.pattern_mix = 1.0f;
      value.inner_warp.cell_x = 0.22321875f;
      value.inner_warp.cell_y = 5.085703f;
      value.inner_warp.rotation = 3.455752f;
      value.inner_warp.speed = 0.0027299998f;
    }
    return value;
  }
  // clang-format on
  // End generated params.
};
