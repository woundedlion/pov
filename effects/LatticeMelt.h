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
  static constexpr std::string_view EFFECT_ID = "lattice-melt";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "e8270e3d8305783f291a6606f3d19b0eb0fc8b59027d079c7f6750bc71a85867";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "8ac39502f71557b9f81e1534c89e8e276e2262687e3a78672e7d8dd7237c5f8e";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"open-curl",
                                                              "dense-curl"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 5;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;

  /// Params the effect starts on, and the base every preset varies from.
  static constexpr Params initial_params() {
    Params value;
    value.source = {.lattice_cell_scale = 0.710265636f,
                    .lattice_shape_blend = 1.0f,
                    .lattice_softness = 0.455532223f,
                    .lattice_radius = 0.290762514f};
    value.projection.wander = 1.0f;
    value.projection.camera_wander = 1.0f;
    value.surface = {
        .scale = 1.78815627f, .strength = 0.0759999976f, .speed = 0.0f};
    value.color.palette_chroma = 1.0f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.color.mapping_phase = -0.165999994f;
    value.color.hue_shift_amount = 0.268000007f;
    value.color.hue_noise_scale = 2.0f;
    return value;
  }

  /**
   * @brief Params for the preset at @p index in PRESET_IDS.
   * @details `open-curl` is the initial preset; `dense-curl` folds the lattice
   * through a shorter surface-noise wavelength.
   */
  static constexpr Params preset_params(size_t index) {
    static_assert(PRESET_IDS.size() == 2,
                  "a new preset id needs a branch below");
    Params value = initial_params();
    if (index == 1)
      value.surface.scale = 3.29720306f;
    return value;
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(LatticeMelt)
