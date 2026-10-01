/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

struct KaleidoscopeFlowersSpec : Pullback::Spec {
  static constexpr Pullback::ProjectionKind PROJECTION =
      Pullback::ProjectionKind::EQUIRECTANGULAR;
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
using KaleidoscopeFlowersParams = Pullback::ParamsFor<KaleidoscopeFlowersSpec>;

/**
 * @brief Dodecahedral grids mapped continuously around the equator.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopeFlowers
    : public Pullback::ComposedEffect<W, H, KaleidoscopeFlowers<W, H>,
                                      KaleidoscopeFlowersSpec> {

public:
  using Params = KaleidoscopeFlowersParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-flowers";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "30721bbdccbe37bc46944c66976d97a25306e9654739848940db2a51d8ee82b3";
  static constexpr std::string_view PRESET_BANK_DIGEST = "d1bc04d4ed1aec4b7720a5d90429388da1933b83cd85b387ac56da440dc51970";
  static constexpr std::array<std::string_view, 3> PRESET_IDS{
      "double-map",
      "open-grid",
      "fine-grid"
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
    value.template get<"color">().hue_noise_speed = 0.0f;
    value.template get<"color">().hue_shift_amount = 0.366f;
    value.template get<"color">().mapping_frequency = 2.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"projection">().central_meridian = 0.0f;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 0.165f;
    value.template get<"projection">().singularity_fade = 2.14f;
    value.template get<"source">().angle_rate = 0.026999999f;
    value.template get<"source">().complexity = 3.0f;
    value.template get<"source">().secondary_rate = 0.8f;
    value.template get<"source">().pattern_freq = 3.9407f;
    value.template get<"source">().pattern_mix = 1.0f;
    value.template get<"source">().speed = 0.0f;
    value.template get<"inner_warp">().cell_x = 1.0471976f;
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
      value.template get<"color">().mapping_frequency = 1.0f;
    }
    if (index == 2) {
      value.template get<"color">().mapping_frequency = 21.212f;
      value.template get<"source">().pattern_freq = 0.3985f;
      value.template get<"inner_warp">().cell_y = 0.90189064f;
      value.template get<"inner_warp">().speed = 0.00058f;
    }
    return {value, Segue::Preset::Lerp{480, math::ease_in_out_sin}};
  }
  // clang-format on
  // End generated params.
};
