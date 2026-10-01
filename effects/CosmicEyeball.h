/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

struct CosmicEyeballSpec : Pullback::Spec {
  static constexpr Pullback::ProjectionKind PROJECTION =
      Pullback::ProjectionKind::STEREOGRAPHIC;
  static constexpr Pullback::TransferKind TRANSFER =
      Pullback::TransferKind::NONE;
  static constexpr Pullback::ProjectionCoverageMode COVERAGE =
      Pullback::ProjectionCoverageMode::EDGE_FADE;
  static constexpr Pullback::FieldCoverageKind FIELD_COVERAGE =
      Pullback::FieldCoverageKind::NONE;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  static constexpr Pullback::HueMode HUE = Pullback::HueMode::PATH_LENGTH;
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
using CosmicEyeballParams = Pullback::ParamsFor<CosmicEyeballSpec>;

/**
 * @brief A high-contrast mirrored grid with displacement-driven hue.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class CosmicEyeball : public Pullback::ComposedEffect<W, H, CosmicEyeball<W, H>,
                                                      CosmicEyeballSpec> {

public:
  using Params = CosmicEyeballParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "cosmic-eyeball";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "61351c10fcbb2dcc23cd69b3039370aa0afe49a204ab05c639a0779c8ef34c1f";
  static constexpr std::string_view PRESET_BANK_DIGEST = "7e92c0d6d6c25c0c35917402fb0d7cfece85c5bcedd06d80712c3df5c7a6fc63";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "mirrored-grid"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 2;
  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.projection.camera_wander = 1.0f;
    value.color.hue_shift_amount = 2.048f;
    value.color.mapping_frequency = 1.0f;
    value.color.mapping_phase = 0.0f;
    value.color.palette_chroma = 0.292f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.color.phase_oscillation_depth = 0.0f;
    value.color.phase_oscillation_speed = 0.0f;
    value.color.opacity_high = 1.0f;
    value.color.opacity_low = 1.0f;
    value.projection.singularity_fade = 1.4f;
    value.source.angle_rate = 0.0f;
    value.source.complexity = 1.854f;
    value.source.secondary_rate = 1.0f;
    value.value.edge_width = 0.5f;
    value.source.pattern_freq = 2.5477f;
    value.source.pattern_mix = 0.0f;
    value.source.speed = 0.235f;
    value.outer_warp.cell_x = 5.381125f;
    value.outer_warp.cell_y = 1.0f;
    value.outer_warp.offset_x = 1.344f;
    value.outer_warp.offset_y = -1.456f;
    value.outer_warp.rotation = 0.29530972f;
    value.outer_warp.speed = 0.0f;
    return value;
  }
  // clang-format on
  // End generated params.
};
