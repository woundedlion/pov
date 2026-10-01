/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

struct KaleidoscopeMandalaSpec : Pullback::Spec {
  static constexpr Pullback::ProjectionKind PROJECTION =
      Pullback::ProjectionKind::GNOMONIC_FOLDED;
  static constexpr Pullback::TransferKind TRANSFER =
      Pullback::TransferKind::NONE;
  static constexpr Pullback::ProjectionCoverageMode COVERAGE =
      Pullback::ProjectionCoverageMode::WEIGHT_SQUARED;
  static constexpr Pullback::FieldCoverageKind FIELD_COVERAGE =
      Pullback::FieldCoverageKind::NONE;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  static constexpr Pullback::HueMode HUE = Pullback::HueMode::NOISE;
  static constexpr Pullback::Color::BrightnessEnvelope BRIGHTNESS =
      Pullback::Color::BrightnessEnvelope::NONE;
  static constexpr bool ANIMATED_PROJECTION = false;
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
          Pullback::WaveShearParams, B, "outer_warp",
          HUE == Pullback::HueMode::PATH_LENGTH>::Type>,
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
using KaleidoscopeMandalaParams = Pullback::ParamsFor<KaleidoscopeMandalaSpec>;

/**
 * @brief Folded-gnomonic wave field reflected through a dodecahedral lens.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopeMandala
    : public Pullback::ComposedEffect<W, H, KaleidoscopeMandala<W, H>,
                                      KaleidoscopeMandalaSpec> {

public:
  using Params = KaleidoscopeMandalaParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "kaleidoscope-mandala";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "b37f0b27686b2dea268c8860d63a062e0a2ad59dd694d1d861139b6b5ab8da3e";
  static constexpr std::string_view PRESET_BANK_DIGEST = "7e10bc72b93c0671877e68f54c36bcb6177b5038f3e12c7af5aebed8c4ef1f56";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "wave-mirror",
      "cup-hue"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 2;

  // Hot entry trampoline; RenderPipeline::shade uses hot flash.
  static HS_HOT_FLASH_MEMBER Color4
  shade(const math::Vector &view,
        const typename KaleidoscopeMandala::Frame &frame) {
    return KaleidoscopeMandala::RenderPipeline::shade(view, frame);
  }
  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 1.0f;
    value.template get<"color">().hue_noise_scale = 1.0f;
    value.template get<"color">().hue_noise_speed = 0.0f;
    value.template get<"color">().hue_shift_amount = 0.721f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"projection">().singularity_fade = 2.311f;
    value.template get<"source">().angle_rate = 0.027f;
    value.template get<"source">().complexity = 1.704f;
    value.template get<"source">().secondary_rate = 0.8f;
    value.template get<"source">().pattern_freq = 6.3287f;
    value.template get<"source">().pattern_mix = 0.0f;
    value.template get<"source">().speed = 0.04f;
    value.template get<"outer_warp">().field_angle = 2.2305307f;
    value.template get<"outer_warp">().frequency = 1.408f;
    value.template get<"outer_warp">().speed = -0.00325f;
    value.template get<"outer_warp">().strength = -0.176f;
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
      value.template get<"color">().hue_noise_scale = 1.9717969f;
      value.template get<"color">().hue_shift_amount = 1.0f;
      value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    }
    return {value, Segue::Preset::Lerp{480, math::ease_in_out_sin}};
  }
  // clang-format on
  // End generated params.
};
