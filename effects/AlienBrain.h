/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

struct AlienBrainSpec : Pullback::Spec {
  static constexpr Pullback::ProjectionKind PROJECTION =
      Pullback::ProjectionKind::STEREOGRAPHIC;
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
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr Pullback::SurfacePlacement SURFACE_PLACEMENT =
      Pullback::SurfacePlacement::BEFORE_LENS;
  using LensPolicy = Pullback::Lens::Glitch;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Lens<Pullback::Lens::Glitch>,
          Pullback::Stage::Project<
              typename Pullback::ProjectionPolicyFor<PROJECTION, B>::Type>>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::WaveShearParams, B, "outer_warp",
          HUE == Pullback::HueMode::PATH_LENGTH>::Type>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::GridSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<COVERAGE, B, void>::Type>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};
using AlienBrainParams = Pullback::ParamsFor<AlienBrainSpec>;

/**
 * @brief Glitch-folded grids pulled through an animated wave shear.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class AlienBrain
    : public Pullback::ComposedEffect<W, H, AlienBrain<W, H>, AlienBrainSpec> {

public:
  using Params = AlienBrainParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "alien-brain";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "e6fad4aa8069af2a40b53f15b011164f24d3ae3b181d7ba00753eab7d0dac41c";
  static constexpr std::string_view PRESET_BANK_DIGEST = "463a8eddd452f7b0b04bda8d2736e92ccc979fac6fe88121b82e7cda888ca0e9";
  static constexpr std::array<std::string_view, 4> PRESET_IDS{
      "alien-brain",
      "alien-brain-2",
      "alien-brain-3",
      "alien-brain-4"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 2;

  // Hot entry trampoline.
  static HS_HOT_FLASH_MEMBER Color4
  shade(const math::Vector &view, const typename AlienBrain::Frame &frame) {
    return AlienBrain::RenderPipeline::shade(view, frame);
  }
  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 0.8f;
    value.template get<"color">().hue_noise_scale = 0.6304219f;
    value.template get<"color">().hue_noise_speed = 0.0f;
    value.template get<"color">().hue_shift_amount = 0.292f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 0.788f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 0.0f;
    value.template get<"projection">().singularity_fade = 1.0f;
    value.template get<"source">().angle_rate = 0.0f;
    value.template get<"source">().complexity = 0.5f;
    value.template get<"source">().secondary_rate = 0.0f;
    value.template get<"source">().pattern_freq = 4.439f;
    value.template get<"source">().pattern_mix = 0.0f;
    value.template get<"source">().speed = 0.245f;
    value.template get<"outer_warp">().field_angle = 0.0f;
    value.template get<"outer_warp">().frequency = 1.0f;
    value.template get<"outer_warp">().speed = 0.015625f;
    value.template get<"outer_warp">().strength = 0.5f;
    return value;
  }

  /** @brief The preset at index in PRESET_IDS and how it departs. */
  HS_COLD_MEMBER static constexpr PresetEntry<Params> preset(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.template get<"source">().pattern_freq = 3.1447f;
      value.template get<"outer_warp">().speed = 0.00690625f;
      value.template get<"outer_warp">().strength = 2.72f;
    }
    if (index == 2) {
      value.template get<"source">().complexity = 1.698f;
      value.template get<"source">().pattern_freq = 7.5227f;
      value.template get<"outer_warp">().speed = 0.00690625f;
      value.template get<"outer_warp">().strength = 0.0f;
    }
    if (index == 3) {
      value.template get<"source">().complexity = 1.698f;
      value.template get<"source">().pattern_freq = 8.8162f;
      value.template get<"outer_warp">().speed = 0.00559375f;
      value.template get<"outer_warp">().strength = 1.376f;
    }
    return {value, Segue::Preset::Lerp{480, math::ease_in_out_sin}};
  }
  // clang-format on
  // End generated params.
};
