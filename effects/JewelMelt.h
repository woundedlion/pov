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

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

/** @brief Pullback stage spec for `JewelMelt`. */
struct JewelMeltSpec : Pullback::Spec {
  /// Samples animation continuously across each two-arm half-turn.
  static constexpr bool TEMPORAL_SAMPLING = true;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::ANALOGOUS;
  static constexpr bool ANIMATED_PROJECTION = true;
  /**
   * @brief Stage chain: camera rotation, octahedral kaleidoscope lens,
   *  surface-noise displacement, projection, square-cell source, generated
   *  palette,
   *  noise hue modulation.
   * @tparam B The effect's Binding.
   */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Lens<typename Pullback::LensPolicyFor<
              Pullback::FixedKaleidoscopeParams<
                  Pullback::KaleidoscopeSymmetry::OCTAHEDRAL>,
              B>::Type>,
          Pullback::Stage::Displace<typename Pullback::SurfacePolicyFor<
              Pullback::SurfaceNoiseParams, B, false>::Type>,
          Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
              Pullback::ProjectionKind::EQUIRECTANGULAR, B>::Type>>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::SquareCellsParams,
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
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::JewelMelt::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "54e1d66840ecdc3e6b6ae0f0e1c255b7745aa438e32d5d1cb30c192482c918f0";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "cc232fb7f40888b32d05c2b92ee735a1d06bd8bf05eacc97e80b17ea153fb4d7";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "jewel-melt"
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
    value.template get<"color">().mapping_phase = 0.16041666f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"modulation">().modulation.mechanic.amount = 0.552f;
    value.template get<"modulation">().modulation.enabled = true;
    value.template get<"modulation">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"modulation">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulation">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulation">().modulation.noise.scale = 0.6304219f;
    value.template get<"modulation">().modulation.noise.speed = 0.00005f;
    value.template get<"projection">().central_meridian = 1.256637f;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 0.1f;
    value.template get<"projection">().singularity_fade = 1.0f;
    value.template get<"source">().cell_radius = 0.41414943f;
    value.template get<"source">().cell_scale = 3.041703f;
    value.template get<"source">().cell_shape = 1.0f;
    value.template get<"source">().cell_softness = 0.2996846f;
    value.template get<"surface">().scale = 4.6311646f;
    value.template get<"surface">().speed = 0.00017089843f;
    value.template get<"surface">().strength = 0.038802084f;
    return value;
  }
  // clang-format on
  // End generated params.

  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, JewelMelt<W, H>,
                                    JewelMeltSpec>::choreography_hooks()
        .template extend<
            JewelMelt,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &JewelMelt::initial_params>>();
  }
};
