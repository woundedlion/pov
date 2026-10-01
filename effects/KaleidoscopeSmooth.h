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

struct KaleidoscopeSmoothSpec : Pullback::Spec {
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
  using LensPolicy = Pullback::Lens::DodecahedralKaleidoscope;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Lens<Pullback::Lens::DodecahedralKaleidoscope>,
          Pullback::Stage::Project<
              typename Pullback::ProjectionPolicyFor<PROJECTION, B>::Type>>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::MirrorParams, B, "inner_warp",
          HUE == Pullback::HueMode::PATH_LENGTH>::Type>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::GridSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<COVERAGE, B, void>::Type>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};
using KaleidoscopeSmoothParams = Pullback::ParamsFor<KaleidoscopeSmoothSpec>;

/**
 * @brief Mirrored grids folded through a dodecahedral stereographic lens.
 * @details Supplies the render pipeline and preset bank; Pullback::ComposedEffect
 * supplies parameter registration, preset choreography and the palette,
 * camera-walk and noise clocks. The dodecahedral fold is the Lens stage and
 * the mirror tiling is the inner planar warp; the pipeline has no surface stage.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopeSmooth
    : public Pullback::ComposedEffect<W, H, KaleidoscopeSmooth<W, H>,
                                      KaleidoscopeSmoothSpec> {
  friend struct ::hs_test::kaleidoscope_smooth_tests::
      KaleidoscopeSmoothWhiteBox;

public:
  using Params = KaleidoscopeSmoothParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-smooth";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "48f48e0af5cb2f041007567815457c26824eafd8a4397e605da1d6ba57438b16";
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
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 4;

  /// Params the effect starts on, and the base every preset varies from.
  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 1.0f;
    value.template get<"color">().hue_noise_scale = 1.4721563f;
    value.template get<"color">().hue_noise_speed = 0.0f;
    value.template get<"color">().hue_shift_amount = 0.366f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 1.0f;
    value.template get<"projection">().singularity_fade = 3.432f;
    value.template get<"source">().angle_rate = 0.026999999f;
    value.template get<"source">().complexity = 0.513f;
    value.template get<"source">().secondary_rate = 0.8f;
    value.template get<"source">().pattern_freq = 2.8263f;
    value.template get<"source">().pattern_mix = 0.0f;
    value.template get<"source">().speed = 0.0f;
    value.template get<"inner_warp">().cell_x = 1.0f;
    value.template get<"inner_warp">().cell_y = 0.99770314f;
    value.template get<"inner_warp">().offset_x = 0.0f;
    value.template get<"inner_warp">().offset_y = 0.0f;
    value.template get<"inner_warp">().rotation = 0.0f;
    value.template get<"inner_warp">().speed = 0.00013f;
    return value;
  }

  /** @brief The preset at index in PRESET_IDS and how it departs. */
  HS_COLD_MEMBER static constexpr PresetEntry<Params> preset(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.template get<"source">().complexity = 3.0f;
      value.template get<"source">().pattern_mix = 1.0f;
    }
    if (index == 2) {
      value.template get<"color">().mapping_frequency = 2.0f;
      value.template get<"projection">().wander = 0.165f;
      value.template get<"source">().complexity = 3.0f;
      value.template get<"source">().pattern_freq = 3.9407f;
      value.template get<"source">().pattern_mix = 1.0f;
    }
    if (index == 3) {
      value.template get<"color">().mapping_frequency = 1.558f;
      value.template get<"projection">().wander = 0.165f;
      value.template get<"source">().complexity = 3.0f;
      value.template get<"source">().pattern_freq = 2.9059f;
      value.template get<"source">().pattern_mix = 1.0f;
      value.template get<"inner_warp">().cell_x = 0.22321875f;
      value.template get<"inner_warp">().cell_y = 5.085703f;
      value.template get<"inner_warp">().rotation = 3.455752f;
      value.template get<"inner_warp">().speed = 0.0027299998f;
    }
    return {value, Segue::Preset::Lerp{480, math::ease_in_out_sin}};
  }
  // clang-format on
  // End generated params.
};
