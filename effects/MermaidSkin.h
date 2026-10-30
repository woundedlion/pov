/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file MermaidSkin.h
 * @brief Composed pullback effect:
 *        a high-chroma folded grid rippling through curl noise.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

/** @brief Stage spec for `MermaidSkin`. */
struct MermaidSkinSpec : Pullback::Spec {
  /// Samples animation continuously across each two-arm half-turn.
  static constexpr bool TEMPORAL_SAMPLING = true;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::ANALOGOUS;
  static constexpr bool ANIMATED_PROJECTION = true;
  /**
   * @brief Stage chain: camera rotation, surface-noise displacement, folded sinusoidal projection, grid source, generated palette.
   * @tparam B The effect's Binding.
   */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::OUT_OF_LINE_ITCM,
          Pullback::Stage::Displace<typename Pullback::SurfacePolicyFor<
              Pullback::SurfaceNoiseParams, B, false>::Type>,
          Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
              Pullback::ProjectionKind::FOLDED_SINUSOIDAL, B>::Type>>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::GridSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<
              Pullback::ProjectionCoverageMode::WEIGHT, B, void>::Type>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulation", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE>>>>;
};
/// Parameter block derived from `MermaidSkinSpec`.
using MermaidSkinParams = Pullback::ParamsFor<MermaidSkinSpec>;

/**
 * @brief A high-chroma folded grid rippling through sphere-space curl noise.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class MermaidSkin : public Pullback::ComposedEffect<W, H, MermaidSkin<W, H>,
                                                    MermaidSkinSpec> {

public:
  using Params = MermaidSkinParams; ///< Live parameter struct.
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::MermaidSkin::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "60803c69068a42f15672c67d96a3ce79e2281d59831ccdd02eb6b2c91a4ec0ad";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "45590a3c8ac175a69ff32ec3260ec8a3612065933a61cf3b056845ef6dbef782";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "mermaid-skin"
  };
  /// Frames each preset holds before advancing.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  /// Params layout version; snapshots tagged with another are rejected.
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
    value.template get<"color">().mapping_frequency = 5.3755207f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 0.9661458f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"modulation">().modulation.mechanic.amount = 1.5958333f;
    value.template get<"modulation">().modulation.enabled = true;
    value.template get<"modulation">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"modulation">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulation">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulation">().modulation.noise.scale = 0.14453939f;
    value.template get<"modulation">().modulation.noise.speed = -0.0000010416667f;
    value.template get<"projection">().central_meridian = 0.0f;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 0.0f;
    value.template get<"source">().angle_rate = 0.0f;
    value.template get<"source">().complexity = 0.0f;
    value.template get<"source">().secondary_rate = 0.0f;
    value.template get<"source">().pattern_freq = 0.1f;
    value.template get<"source">().pattern_mix = 0.0f;
    value.template get<"source">().speed = 0.0f;
    value.template get<"surface">().scale = 4.9144287f;
    value.template get<"surface">().speed = -0.00021158854f;
    value.template get<"surface">().strength = 0.5f;
    return value;
  }
  // clang-format on
  // End generated params.

public:
  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, MermaidSkin<W, H>,
                                    MermaidSkinSpec>::choreography_hooks()
        .template extend<
            MermaidSkin,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &MermaidSkin::initial_params>>();
  }
};
