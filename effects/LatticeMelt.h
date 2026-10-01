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

struct LatticeMeltSpec : Pullback::Spec {
  static constexpr Pullback::ProjectionKind PROJECTION =
      Pullback::ProjectionKind::FOLDED_SINUSOIDAL;
  static constexpr Pullback::TransferKind TRANSFER =
      Pullback::TransferKind::NONE;
  static constexpr Pullback::ProjectionCoverageMode COVERAGE =
      Pullback::ProjectionCoverageMode::WEIGHT;
  static constexpr Pullback::FieldCoverageKind FIELD_COVERAGE =
      Pullback::FieldCoverageKind::NONE;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  static constexpr Pullback::HueMode HUE = Pullback::HueMode::NOISE;
  static constexpr Pullback::Color::BrightnessEnvelope BRIGHTNESS =
      Pullback::Color::BrightnessEnvelope::CUP;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr Pullback::SurfacePlacement SURFACE_PLACEMENT =
      Pullback::SurfacePlacement::BEFORE_LENS;
  using LensPolicy = void;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::OUT_OF_LINE_FLASH,
          Pullback::Stage::Displace<typename Pullback::SurfacePolicyFor<
              Pullback::SurfaceNoiseParams, B,
              HUE == Pullback::HueMode::PATH_LENGTH>::Type>,
          Pullback::Stage::Project<
              typename Pullback::ProjectionPolicyFor<PROJECTION, B>::Type>>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::LatticeSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<COVERAGE, B, void>::Type>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};
using LatticeMeltParams = Pullback::ParamsFor<LatticeMeltSpec>;

/**
 * @brief Composed folded-sinusoidal lattice displaced by sphere-space curl noise.
 * @details Supplies the render pipeline and preset bank; Pullback::ComposedEffect
 * supplies parameter registration, preset choreography and the palette,
 * camera-walk and noise clocks. The lattice source is read through a folded
 * sinusoidal projection; the surface stage carries curl displacement and the
 * pipeline has no warp stage.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class LatticeMelt : public Pullback::ComposedEffect<W, H, LatticeMelt<W, H>,
                                                    LatticeMeltSpec> {
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
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 6;

  /// Params the effect starts on, and the base every preset varies from.
  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 1.0f;
    value.template get<"color">().brightness_bottom = 0.0f;
    value.template get<"color">().brightness_top = 1.0f;
    value.template get<"color">().hue_noise_scale = 2.0f;
    value.template get<"color">().hue_noise_speed = 0.0f;
    value.template get<"color">().hue_shift_amount = 0.268f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = -0.166f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"projection">().central_meridian = 0.0f;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 1.0f;
    value.template get<"source">().lattice_cell_scale = 0.71026564f;
    value.template get<"source">().lattice_radius = 0.2907625f;
    value.template get<"source">().lattice_shape_blend = 1.0f;
    value.template get<"source">().lattice_softness = 0.45553222f;
    value.template get<"surface">().scale = 1.7881563f;
    value.template get<"surface">().speed = 0.0f;
    value.template get<"surface">().strength = 0.076f;
    return value;
  }

  /** @brief The preset at index in PRESET_IDS and how it departs. */
  HS_COLD_MEMBER static constexpr PresetEntry<Params> preset(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.template get<"surface">().scale = 3.297203f;
    }
    return {value, Segue::Preset::Lerp{480, math::ease_in_out_sin}};
  }
  // clang-format on
  // End generated params.
};
