/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

struct AlienCoreSpec : Pullback::Spec {
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
using AlienCoreParams = Pullback::ParamsFor<AlienCoreSpec>;

/**
 * @brief A mirrored grid folded by the glitch lens.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class AlienCore
    : public Pullback::ComposedEffect<W, H, AlienCore<W, H>, AlienCoreSpec> {

public:
  using Params = AlienCoreParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "alien-core";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "65a245d1dfbbc68eff720abe0fd1c6195e82954412646ee38e1d7a41397b3472";
  static constexpr std::string_view PRESET_BANK_DIGEST = "8b76fda1895b0572c30c79e00baa135e3bbcc151c2c051f3c1ab0773be83c366";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "folded-glitch"
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
    value.template get<"color">().hue_noise_scale = 1.0f;
    value.template get<"color">().hue_noise_speed = 0.0f;
    value.template get<"color">().hue_shift_amount = 0.0f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 0.62f;
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
