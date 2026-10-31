/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file GridSpace.h
 * @brief Composed pullback effect:
 *        an affine primitive lattice rendered as soft contours.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

/** @brief Pullback stage spec for `GridSpace`. */
struct GridSpaceSpec : Pullback::Spec {
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  static constexpr bool ANIMATED_PROJECTION = true;
  /// Samples animation continuously across each two-arm half-turn.
  static constexpr bool TEMPORAL_SAMPLING = true;
  /// Preserves canonical rotation order at gamut-mapping boundaries.
  static constexpr bool FACTORED_COLUMNS = false;
  /**
   * @brief Stage chain: camera rotation, projection, affine warp, lattice
   *  source, iso-contour transfer, generated palette,
   *  noise hue modulation.
   * @tparam B The effect's Binding.
   */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
              Pullback::ProjectionKind::GNOMONIC_FOLDED, B>::Type>>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::AffineParams, B, "outer_warp", false>::Type>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
                                  Pullback::LatticeSourceParams, B>::Type,
                              Pullback::Weight::Projection,
                              typename Pullback::CoveragePolicyFor<
                                  Pullback::ProjectionCoverageMode::WEIGHT, B,
                                  Pullback::IsoValueParams>::Type>,
      Pullback::Stage::Transfer<Pullback::Transfer::IsoContour<
          Pullback::ValueProvider<B, Pullback::IsoValueParams>>>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulation", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE>>>>;
};
/// Parameter block derived from `GridSpaceSpec`.
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
  /// Parameter block derived from `GridSpaceSpec`.
  using Params = GridSpaceParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::GridSpace::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "5838236c430e3b792987b36869eedcbb48c4e1806ff08ece54dd32adba07e46d";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "78274582242e8c26030356ce8fe180accc40e35c3811f51d5f803ed267680496";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "affine-contour"
  };
  /// Frames each preset holds before advancing.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  /// Snapshot schema version; changes with the `Params` layout.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 4;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  /**
   * @brief Parameters of the first preset; every preset varies from them.
   * @return The preset-0 `Params`.
   */
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 1.0f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 0.62f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"modulation">().modulation.mechanic.amount = 0.398f;
    value.template get<"modulation">().modulation.enabled = true;
    value.template get<"modulation">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"modulation">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulation">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulation">().modulation.noise.scale = 0.8300313f;
    value.template get<"modulation">().modulation.noise.speed = 0.00021200001f;
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

public:
  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, GridSpace<W, H>,
                                    GridSpaceSpec>::choreography_hooks()
        .template extend<
            GridSpace,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &GridSpace::initial_params>>();
  }
};
