/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

struct GridSpaceSpec : Pullback::Spec {
  static constexpr Pullback::ProjectionKind PROJECTION =
      Pullback::ProjectionKind::GNOMONIC_FOLDED;
  static constexpr Pullback::TransferKind TRANSFER =
      Pullback::TransferKind::ISO_CONTOUR;
  static constexpr Pullback::ProjectionCoverageMode COVERAGE =
      Pullback::ProjectionCoverageMode::WEIGHT;
  static constexpr Pullback::FieldCoverageKind FIELD_COVERAGE =
      Pullback::FieldCoverageKind::NONE;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  static constexpr Pullback::HueMode HUE = Pullback::HueMode::NOISE;
  static constexpr Pullback::Color::BrightnessEnvelope BRIGHTNESS =
      Pullback::Color::BrightnessEnvelope::NONE;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr Pullback::SurfacePlacement SURFACE_PLACEMENT =
      Pullback::SurfacePlacement::BEFORE_LENS;
  using LensPolicy = void;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Project<
              typename Pullback::ProjectionPolicyFor<PROJECTION, B>::Type>>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::AffineParams, B, "outer_warp",
          HUE == Pullback::HueMode::PATH_LENGTH>::Type>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
                                  Pullback::LatticeSourceParams, B>::Type,
                              Pullback::Weight::Projection,
                              typename Pullback::CoveragePolicyFor<
                                  COVERAGE, B, Pullback::IsoValueParams>::Type>,
      Pullback::Stage::Transfer<Pullback::Transfer::IsoContour<
          Pullback::ValueProvider<B, Pullback::IsoValueParams>>>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};
using GridSpaceParams = Pullback::ParamsFor<GridSpaceSpec>;

/**
 * @brief An affine primitive lattice rendered as soft contours.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class GridSpace
    : public Pullback::ComposedEffect<W, H, GridSpace<W, H>, GridSpaceSpec> {

public:
  using Params = GridSpaceParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "grid-space";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "67870f99cc8812d737d30ef890622c30f49f49172af1fd649dcd49821c5d1aca";
  static constexpr std::string_view PRESET_BANK_DIGEST = "992256905cd236a8dae29bbafc058efee46c7a98751abcf477ad1622ef0dd535";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "affine-contour"
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
    value.template get<"color">().hue_noise_scale = 0.8300313f;
    value.template get<"color">().hue_noise_speed = 0.00021200001f;
    value.template get<"color">().hue_shift_amount = 0.398f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 0.62f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"projection">().spin_rate = 0.020879198f;
    value.template get<"projection">().wander = 0.0030917525f;
    value.template get<"projection">().singularity_fade = 1.0f;
    value.template get<"source">().lattice_cell_scale = 1.22925f;
    value.template get<"source">().lattice_radius = 0.33298188f;
    value.template get<"source">().lattice_shape_blend = 1.0f;
    value.template get<"source">().lattice_softness = 0.1608203f;
    value.template get<"value">().iso_level = 0.138f;
    value.template get<"value">().iso_width = 0.22703418f;
    value.template get<"outer_warp">().rotation_rate = 0.0f;
    value.template get<"outer_warp">().scale_x = 1.0f;
    value.template get<"outer_warp">().scale_y = 1.0f;
    value.template get<"outer_warp">().shear = 0.0f;
    value.template get<"outer_warp">().speed = 0.015625f;
    value.template get<"outer_warp">().translation_x = 4.0f;
    value.template get<"outer_warp">().translation_y = 4.0f;
    return value;
  }
  // clang-format on
  // End generated params.
};
