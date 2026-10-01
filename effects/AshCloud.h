/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

struct AshCloudSpec : Pullback::Spec {
  static constexpr Pullback::ProjectionKind PROJECTION =
      Pullback::ProjectionKind::STEREOGRAPHIC;
  static constexpr Pullback::TransferKind TRANSFER =
      Pullback::TransferKind::NONE;
  static constexpr Pullback::ProjectionCoverageMode COVERAGE =
      Pullback::ProjectionCoverageMode::WEIGHT;
  static constexpr Pullback::FieldCoverageKind FIELD_COVERAGE =
      Pullback::FieldCoverageKind::VALUE_CUTOUT;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
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
          Pullback::CodeEmission::OUT_OF_LINE_FLASH,
          Pullback::Stage::Displace<typename Pullback::SurfacePolicyFor<
              Pullback::SurfaceNoiseParams, B,
              HUE == Pullback::HueMode::PATH_LENGTH>::Type>,
          Pullback::Stage::Lens<Pullback::Lens::DodecahedralKaleidoscope>,
          Pullback::Stage::Project<
              typename Pullback::ProjectionPolicyFor<PROJECTION, B>::Type>>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::LatticeSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<
              COVERAGE, B, Pullback::CutoutValueParams>::Type>,
      Pullback::Stage::ApplyCoverage<Pullback::ValueCoverage::ValueCutout<
          Pullback::ValueProvider<B, Pullback::CutoutValueParams>>>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};
using AshCloudParams = Pullback::ParamsFor<AshCloudSpec>;

/**
 * @brief A softly cut lattice curled across dodecahedral facets.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class AshCloud
    : public Pullback::ComposedEffect<W, H, AshCloud<W, H>, AshCloudSpec> {

public:
  using Params = AshCloudParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "ash-cloud";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "6bde4043e907674b80f6b12bf66c873366359c9a7b1996e8804ae24da673518e";
  static constexpr std::string_view PRESET_BANK_DIGEST = "bb8b1dc6fe15beaf9a4c11730c902d9743a3035942e3b123e3760e25e50c214e";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "ash-cloud"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr float CAMERA_SPIN_RATE = 0.01975f;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 2;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 1.0f;
    value.template get<"color">().hue_noise_scale = 2.4029531f;
    value.template get<"color">().hue_noise_speed = 0.00015f;
    value.template get<"color">().hue_shift_amount = -3.248f;
    value.template get<"color">().mapping_frequency = 3.666f;
    value.template get<"color">().mapping_phase = 0.598f;
    value.template get<"color">().palette_chroma = 0.314f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().phase_oscillation_depth = 0.856f;
    value.template get<"color">().phase_oscillation_speed = 0.00448f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"value">().cutout_softness = 0.31635937f;
    value.template get<"value">().cutout_threshold = 0.365f;
    value.template get<"projection">().spin_rate = 0.00175f;
    value.template get<"projection">().wander = 0.009f;
    value.template get<"projection">().singularity_fade = 1.0f;
    value.template get<"source">().lattice_cell_scale = 0.98971874f;
    value.template get<"source">().lattice_radius = 0.37757313f;
    value.template get<"source">().lattice_shape_blend = 0.0f;
    value.template get<"source">().lattice_softness = 0.34563965f;
    value.template get<"surface">().scale = 5.7742186f;
    value.template get<"surface">().speed = 0.00040625f;
    value.template get<"surface">().strength = 0.043f;
    return value;
  }
  // clang-format on
  // End generated params.
};
