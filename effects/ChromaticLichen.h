/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file ChromaticLichen.h
 * @brief Composed pullback effect:
 *        a glitch-folded gnomonic grid displaced by curl noise.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

/** @brief Pullback stage spec for `ChromaticLichen`. */
struct ChromaticLichenSpec : Pullback::Spec {
  static constexpr bool MODULATION = true;
  /// Samples animation continuously across each two-arm half-turn.
  static constexpr bool TEMPORAL_SAMPLING = true;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::ANALOGOUS;
  static constexpr bool ANIMATED_PROJECTION = true;
  /**
   * @brief Stage chain: camera rotation, glitch lens, surface-noise
   *  displacement, projection, grid source, generated palette.
   * @tparam B The effect's Binding.
   */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::OUT_OF_LINE_FLASH,
          Pullback::Stage::Lens<Pullback::Lens::Glitch>,
          Pullback::Stage::Displace<typename Pullback::SurfacePolicyFor<
              Pullback::SurfaceNoiseParams, B, false>::Type>,
          Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
              Pullback::ProjectionKind::GNOMONIC_FOLDED, B>::Type>>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::GridSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<
              Pullback::ProjectionCoverageMode::WEIGHT, B, void>::Type>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, Pullback::HueMode::NONE,
                                  Pullback::Color::BrightnessEnvelope::NONE>>>>;
};
/// Parameter block derived from `ChromaticLichenSpec`.
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
  /// Parameter block derived from `ChromaticLichenSpec`.
  using Params = ChromaticLichenParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::ChromaticLichen::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "3993d18300f4c33d66503546c1bebfd09242ec1359d69d53075e5aead23df857";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "9dfc768750d59612813b2207a6d8346bedda971d6620e80423845f41207e37fc";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "chromatic-lichen"
  };
  /// Frames each preset holds before advancing.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  /// Snapshot schema version; changes with the `Params` layout.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 3;

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
    value.template get<"color">().modulation.mechanic.amount = 1.5958333f;
    value.template get<"color">().modulation.enabled = true;
    value.template get<"color">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"color">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"color">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().modulation.noise.scale = 0.14453939f;
    value.template get<"color">().modulation.noise.speed = -0.0000010416667f;
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

public:
  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, ChromaticLichen<W, H>,
                                    ChromaticLichenSpec>::choreography_hooks()
        .template extend<
            ChromaticLichen,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &ChromaticLichen::initial_params>>();
  }
};
