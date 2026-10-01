/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

struct KaleidoscopeHexBrightSpec : Pullback::Spec {
  static constexpr Pullback::ProjectionKind PROJECTION =
      Pullback::ProjectionKind::STEREOGRAPHIC;
  static constexpr Pullback::TransferKind TRANSFER =
      Pullback::TransferKind::NONE;
  static constexpr Pullback::ProjectionCoverageMode COVERAGE =
      Pullback::ProjectionCoverageMode::WEIGHT_SQUARED;
  static constexpr Pullback::FieldCoverageKind FIELD_COVERAGE =
      Pullback::FieldCoverageKind::NONE;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::ANALOGOUS;
  static constexpr Pullback::HueMode HUE = Pullback::HueMode::NOISE;
  static constexpr Pullback::Color::BrightnessEnvelope BRIGHTNESS =
      Pullback::Color::BrightnessEnvelope::NONE;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr Pullback::SurfacePlacement SURFACE_PLACEMENT =
      Pullback::SurfacePlacement::BEFORE_LENS;
  using LensPolicy = Pullback::Lens::HexagonalPrismKaleidoscope;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Lens<Pullback::Lens::HexagonalPrismKaleidoscope>,
          Pullback::Stage::Project<
              typename Pullback::ProjectionPolicyFor<PROJECTION, B>::Type>>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::MirrorParams, B, "inner_warp",
          HUE == Pullback::HueMode::PATH_LENGTH>::Type>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::TwinWaveSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<COVERAGE, B, void>::Type>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};
using KaleidoscopeHexBrightParams =
    Pullback::ParamsFor<KaleidoscopeHexBrightSpec>;

/**
 * @brief A twin wave folded through a hexagonal prism.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopeHexBright
    : public Pullback::ComposedEffect<W, H, KaleidoscopeHexBright<W, H>,
                                      KaleidoscopeHexBrightSpec> {

public:
  using Params = KaleidoscopeHexBrightParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-hex-bright";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "b230ef82f2927252edff29f95beeac6a8e43af8d18514fb6f4e215e118f4b448";
  static constexpr std::string_view PRESET_BANK_DIGEST = "a69686e051b4d46949179a00b82e08df04f1b5abd7b5317376da86753670719b";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "hex-twin-wave",
      "hex-twin-wave-alt"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 2;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 1.0f;
    value.template get<"color">().hue_noise_scale = 1.4721563f;
    value.template get<"color">().hue_noise_speed = 0.000138f;
    value.template get<"color">().hue_shift_amount = 0.226f;
    value.template get<"color">().mapping_frequency = 1.341f;
    value.template get<"color">().mapping_phase = -1.0f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::BELL;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 1.0f;
    value.template get<"projection">().singularity_fade = 4.971f;
    value.template get<"source">().angle_rate = 0.027f;
    value.template get<"source">().secondary_rate = 0.8f;
    value.template get<"source">().pattern_freq = 3.881f;
    value.template get<"source">().speed = 0.12859823f;
    value.template get<"inner_warp">().cell_x = 1.0f;
    value.template get<"inner_warp">().cell_y = 1.0f;
    value.template get<"inner_warp">().offset_x = 0.0f;
    value.template get<"inner_warp">().offset_y = 0.0f;
    value.template get<"inner_warp">().rotation = 0.0f;
    value.template get<"inner_warp">().speed = 0.0f;
    return value;
  }

  /** @brief The preset at index in PRESET_IDS and how it departs. */
  HS_COLD_MEMBER static constexpr PresetEntry<Params> preset(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.template get<"color">().mapping_frequency = 2.0f;
    }
    return {value, Segue::Preset::Lerp{480, math::ease_in_out_sin}};
  }
  // clang-format on
  // End generated params.
};
