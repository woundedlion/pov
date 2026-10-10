/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file JewelMelt.h
 * @brief Composed pullback effect:
 *        an octahedral jewel lattice rippling through curl noise.
 */

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

/** @brief Pullback stage spec for `JewelMelt`. */
struct JewelMeltSpec : Pullback::Spec {
  static constexpr Pullback::ProjectionKind PROJECTION =
      Pullback::ProjectionKind::EQUIRECTANGULAR;
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
  using LensPolicy = Pullback::Lens::OctahedralKaleidoscope;
  /**
   * @brief Stage chain: camera rotation, octahedral kaleidoscope lens,
   *  surface-noise displacement, projection, lattice source, generated
   *  palette.
   * @tparam B The effect's Binding.
   */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::OUT_OF_LINE_FLASH,
          Pullback::Stage::Lens<LensPolicy>,
          Pullback::Stage::Displace<typename Pullback::SurfacePolicyFor<
              Pullback::SurfaceNoiseParams, B,
              HUE == Pullback::HueMode::PATH_LENGTH>::Type>,
          Pullback::Stage::Project<
              typename Pullback::ProjectionPolicyFor<PROJECTION, B>::Type>>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::LatticeSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<COVERAGE, B, void>::Type>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};
/// Parameter block derived from `JewelMeltSpec`.
using JewelMeltParams = Pullback::ParamsFor<JewelMeltSpec>;

/**
 * @brief An octahedral jewel lattice rippling through sphere-space curl noise.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class JewelMelt
    : public Pullback::ComposedEffect<W, H, JewelMelt<W, H>, JewelMeltSpec> {

public:
  /// Parameter block derived from `JewelMeltSpec`.
  using Params = JewelMeltParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "jewel-melt";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "4a5674fa7a9af30c024c518e3057dcd38fef427f77f196ead79ef5546ca6d697";
  static constexpr std::string_view PRESET_BANK_DIGEST = "8919dafdfd00f6b4a9de07524013905eb1acdfd78d9223d6ed92929a408df4fa";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "jewel-melt"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  /** @var EFFECT_ID
   *  @brief Stable effect ID; equals the shader document's effect_id. */
  /** @var DESCRIPTOR_DIGEST
   *  @brief SHA-256 hex of the shader document's parameter descriptor. */
  /** @var PRESET_BANK_DIGEST
   *  @brief SHA-256 hex of the shader document's preset bank. */
  /** @var PRESET_IDS
   *  @brief Preset identities, indexed by preset number. */
  /** @var PRESET_DWELL_FRAMES
   *  @brief Frames each preset holds before advancing. */
  /// Snapshot schema version; changes with the `Params` layout.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 1.0f;
    value.template get<"color">().hue_noise_scale = 0.14453939f;
    value.template get<"color">().hue_noise_speed = -0.0000010416667f;
    value.template get<"color">().hue_shift_amount = 1.5958333f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.16041666f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"projection">().central_meridian = 0.0f;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 0.0f;
    value.template get<"projection">().singularity_fade = 1.0f;
    value.template get<"source">().lattice_cell_scale = 3.3299723f;
    value.template get<"source">().lattice_radius = 0.41414943f;
    value.template get<"source">().lattice_shape_blend = 1.0f;
    value.template get<"source">().lattice_softness = 0.5f;
    value.template get<"surface">().scale = 4.6311646f;
    value.template get<"surface">().speed = 0.00017089843f;
    value.template get<"surface">().strength = 0.038802084f;
    return value;
  }
  // clang-format on
  // End generated params.
  /** @fn initial_params()
   *  @brief Parameters of the first preset.
   *  @return The preset-0 `Params`. */
};
