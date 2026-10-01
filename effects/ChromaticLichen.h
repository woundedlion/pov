/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

struct ChromaticLichenSpec : Pullback::Spec {
  static constexpr Pullback::ProjectionKind PROJECTION =
      Pullback::ProjectionKind::GNOMONIC_FOLDED;
  static constexpr Pullback::TransferKind TRANSFER =
      Pullback::TransferKind::NONE;
  static constexpr Pullback::ProjectionCoverageMode COVERAGE =
      Pullback::ProjectionCoverageMode::WEIGHT;
  static constexpr Pullback::FieldCoverageKind FIELD_COVERAGE =
      Pullback::FieldCoverageKind::NONE;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::ANALOGOUS;
  static constexpr Pullback::HueMode HUE = Pullback::HueMode::NOISE;
  static constexpr Pullback::Color::BrightnessEnvelope BRIGHTNESS =
      Pullback::Color::BrightnessEnvelope::NONE;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr Pullback::SurfacePlacement SURFACE_PLACEMENT =
      Pullback::SurfacePlacement::AFTER_LENS;
  using LensPolicy = Pullback::Lens::Glitch;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::OUT_OF_LINE_FLASH,
          Pullback::Stage::Lens<Pullback::Lens::Glitch>,
          Pullback::Stage::Displace<typename Pullback::SurfacePolicyFor<
              Pullback::SurfaceNoiseParams, B,
              HUE == Pullback::HueMode::PATH_LENGTH>::Type>,
          Pullback::Stage::Project<
              typename Pullback::ProjectionPolicyFor<PROJECTION, B>::Type>>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::GridSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<COVERAGE, B, void>::Type>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};
using ChromaticLichenParams = Pullback::ParamsFor<ChromaticLichenSpec>;

/**
 * @brief A glitch-folded gnomonic grid displaced by sphere-space curl noise.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class ChromaticLichen
    : public Pullback::ComposedEffect<W, H, ChromaticLichen<W, H>,
                                      ChromaticLichenSpec> {

public:
  using Params = ChromaticLichenParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "chromatic-lichen";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "86311a6f7c59e89dc9aaec7ffa7929cd249d40d08408873e094d96a4b9fdcf29";
  static constexpr std::string_view PRESET_BANK_DIGEST = "8a6c7cda1ac938cd2f39211a251841a586ce134ee750683b3c1a9d2d2d524f63";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "chromatic-lichen"
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
    value.template get<"color">().hue_noise_scale = 0.14453939f;
    value.template get<"color">().hue_noise_speed = -0.0000010416667f;
    value.template get<"color">().hue_shift_amount = 1.5958333f;
    value.template get<"color">().mapping_frequency = 5.3755207f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 1.0f;
    value.template get<"projection">().singularity_fade = 1.0f;
    value.template get<"source">().angle_rate = 0.0f;
    value.template get<"source">().complexity = 0.0f;
    value.template get<"source">().secondary_rate = 0.0f;
    value.template get<"source">().pattern_freq = 0.1f;
    value.template get<"source">().pattern_mix = 0.0f;
    value.template get<"source">().speed = 0.0f;
    value.template get<"surface">().scale = 2.3483887f;
    value.template get<"surface">().speed = -0.00021158854f;
    value.template get<"surface">().strength = 0.5f;
    return value;
  }
  // clang-format on
  // End generated params.
};
