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
namespace lattice_melt_tests {
struct LatticeMeltWhiteBox;
} // namespace lattice_melt_tests
} // namespace hs_test

using LatticeMeltParams =
    Pullback::Params<Pullback::LatticeSourceParams, Pullback::NoWarpParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::NoValueParams, Pullback::SurfaceNoiseParams>;
using LatticeMeltSpec =
    Pullback::Spec<Pullback::ProjectionKind::FOLDED_SINUSOIDAL, void,
                   Pullback::TransferKind::NONE,
                   Pullback::ProjectionCoverageMode::WEIGHT>;

/**
 * @brief Composed folded-sinusoidal lattice displaced by sphere-space curl noise.
 * @details Supplies the render pipeline and preset bank; Pullback::ComposedEffect
 * supplies parameter registration, preset choreography and the palette,
 * camera-walk and noise clocks. The lattice source is read through a folded
 * sinusoidal projection, so the surface stage carries the curl displacement and
 * the warp stage is an identity.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class LatticeMelt
    : public Pullback::ComposedEffect<
          W, H, LatticeMelt<W, H>, LatticeMeltParams, LatticeMeltSpec,
          PaletteHarmony::TRIADIC, Pullback::HueMode::NOISE,
          Pullback::Color::BrightnessEnvelope::CUP> {
  friend struct ::hs_test::lattice_melt_tests::LatticeMeltWhiteBox;

public:
  using Params = LatticeMeltParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "lattice-melt";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "61022b59731fab1ec0646b1626bff40fc4e26966800c7c06369357d824f1b241";
  static constexpr std::string_view PRESET_BANK_DIGEST = "8ac39502f71557b9f81e1534c89e8e276e2262687e3a78672e7d8dd7237c5f8e";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "open-curl",
      "dense-curl"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 5;

  /// Params the effect starts on, and the base every preset varies from.
  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.projection.camera_wander = 1.0f;
    value.color.brightness_bottom = 0.0f;
    value.color.brightness_top = 1.0f;
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
    value.projection.central_meridian = 0.0f;
    value.projection.spin_rate = 0.0f;
    value.projection.wander = 1.0f;
    value.source.lattice_cell_scale = 0.71026564f;
    value.source.lattice_radius = 0.2907625f;
    value.source.lattice_shape_blend = 1.0f;
    value.source.lattice_softness = 0.45553222f;
    value.surface.scale = 1.7881563f;
    value.surface.speed = 0.0f;
    value.surface.strength = 0.076f;
    return value;
  }

  /** @brief Params for the preset at index in PRESET_IDS. */
  static constexpr Params preset_params(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.surface.scale = 3.297203f;
    }
    return value;
  }
  // clang-format on
  // End generated params.
};
