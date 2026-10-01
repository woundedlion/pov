/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

struct AlienOceanSpec : Pullback::Spec {
  static constexpr Pullback::ProjectionKind PROJECTION =
      Pullback::ProjectionKind::GNOMONIC_FOLDED;
  static constexpr Pullback::TransferKind TRANSFER =
      Pullback::TransferKind::NONE;
  static constexpr Pullback::ProjectionCoverageMode COVERAGE =
      Pullback::ProjectionCoverageMode::EDGE_FADE;
  static constexpr Pullback::FieldCoverageKind FIELD_COVERAGE =
      Pullback::FieldCoverageKind::NONE;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  static constexpr Pullback::HueMode HUE = Pullback::HueMode::NOISE;
  static constexpr Pullback::Color::BrightnessEnvelope BRIGHTNESS =
      Pullback::Color::BrightnessEnvelope::NONE;
  static constexpr bool ANIMATED_PROJECTION = false;
  static constexpr Pullback::SurfacePlacement SURFACE_PLACEMENT =
      Pullback::SurfacePlacement::BEFORE_LENS;
  using LensPolicy = Pullback::Lens::Kaleidoscope;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Lens<Pullback::Lens::Kaleidoscope>,
          Pullback::Stage::Project<
              typename Pullback::ProjectionPolicyFor<PROJECTION, B>::Type>>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::MirrorParams, B, "outer_warp",
          HUE == Pullback::HueMode::PATH_LENGTH>::Type>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::GridSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<
              COVERAGE, B, Pullback::EdgeValueParams>::Type>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};
using AlienOceanParams = Pullback::ParamsFor<AlienOceanSpec>;

/**
 * @brief A broad folded grid with slow mirrored drift.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class AlienOcean
    : public Pullback::ComposedEffect<W, H, AlienOcean<W, H>, AlienOceanSpec> {

public:
  using Params = AlienOceanParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "alien-ocean";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "2ff3688974c2200915e43dcdad54c1c83fea2e0a312b46ad6eecbc409bf90a57";
  static constexpr std::string_view PRESET_BANK_DIGEST = "a9070db6d46b6ed7aaaf3341a43cf3cf2486948e78b55384f5d48877d92e5061";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "folded-grid"
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
    value.template get<"color">().hue_noise_scale = 2.2033439f;
    value.template get<"color">().hue_noise_speed = 0.0f;
    value.template get<"color">().hue_shift_amount = 0.424f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 0.4f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"projection">().singularity_fade = 1.4f;
    value.template get<"source">().angle_rate = 0.0f;
    value.template get<"source">().complexity = 0.0f;
    value.template get<"source">().secondary_rate = 1.0f;
    value.template get<"value">().edge_width = 0.5f;
    value.template get<"source">().pattern_freq = 3.565f;
    value.template get<"source">().pattern_mix = 1.0f;
    value.template get<"source">().speed = 0.235f;
    value.template get<"outer_warp">().cell_x = 5.381125f;
    value.template get<"outer_warp">().cell_y = 1.0f;
    value.template get<"outer_warp">().offset_x = 1.344f;
    value.template get<"outer_warp">().offset_y = -1.456f;
    value.template get<"outer_warp">().rotation = 0.29530972f;
    value.template get<"outer_warp">().speed = 0.0f;
    return value;
  }
  // clang-format on
  // End generated params.
};
