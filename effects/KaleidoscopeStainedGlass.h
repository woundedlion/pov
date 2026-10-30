/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file KaleidoscopeStainedGlass.h
 * @brief Composed pullback effect:
 *        a vector-noise grid refracted across dodecahedral facets.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

/** @brief Stage spec for `KaleidoscopeStainedGlass`. */
struct KaleidoscopeStainedGlassSpec : Pullback::Spec {
  /// Samples animation continuously across each two-arm half-turn.
  static constexpr bool TEMPORAL_SAMPLING = true;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  static constexpr bool ANIMATED_PROJECTION = false;
  /**
   * @brief Stage chain: camera rotation, dodecahedral kaleidoscope lens, folded gnomonic projection, vector-noise warp, mirror warp, grid source, generated palette.
   * @tparam B The effect's Binding.
   */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Lens<Pullback::Lens::DodecahedralKaleidoscope>,
          Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
              Pullback::ProjectionKind::GNOMONIC_FOLDED, B>::Type>>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::VectorNoiseParams, B, "outer_warp", false>::Type>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::MirrorParams, B, "inner_warp", false>::Type>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::GridSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<
              Pullback::ProjectionCoverageMode::WEIGHT_SQUARED, B, void>::Type>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
              B, Pullback::Color::BrightnessEnvelope::CUP>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulation", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE>>>>;
};
/// Parameter block derived from `KaleidoscopeStainedGlassSpec`.
using KaleidoscopeStainedGlassParams =
    Pullback::ParamsFor<KaleidoscopeStainedGlassSpec>;

/**
 * @brief A vector-noise grid refracted across dodecahedral facets.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopeStainedGlass
    : public Pullback::ComposedEffect<W, H, KaleidoscopeStainedGlass<W, H>,
                                      KaleidoscopeStainedGlassSpec> {

public:
  using Params = KaleidoscopeStainedGlassParams; ///< Live parameter struct.
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::KaleidoscopeStainedGlass::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "1b7502c57567241a8f5c23676e83b1c43da39f1f98cb2335f09740665c88891d";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "e71dcfa9c03690a95de7c1a981eaf38be2e9c5b20ed94ebccaf25aafb89e5adb";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "vector-mirror"
  };
  /// Frames each preset holds before advancing.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  /// Params layout version; snapshots tagged with another are rejected.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 4;

  // Hot entry trampoline.
  /**
   * @brief Shades one pixel through the render pipeline.
   * @param view Unit view direction for the pixel.
   * @param frame Per-frame transforms, params and LUTs from the runtime.
   * @return The pixel colour.
   */
  static HS_HOT_FLASH_MEMBER Color4
  shade(const math::Vector &view,
        const typename KaleidoscopeStainedGlass::Frame &frame) {
    return KaleidoscopeStainedGlass::RenderPipeline::shade(view, frame);
  }
  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  /**
   * @brief Parameters of the first preset; every preset varies from them.
   * @return The preset-0 `Params`.
   */
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 1.0f;
    value.template get<"color">().brightness_bottom = 0.345f;
    value.template get<"color">().brightness_top = 1.0f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"modulation">().modulation.mechanic.amount = 0.721f;
    value.template get<"modulation">().modulation.enabled = true;
    value.template get<"modulation">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"modulation">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulation">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulation">().modulation.noise.scale = 1.0f;
    value.template get<"modulation">().modulation.noise.speed = 0.0f;
    value.template get<"projection">().singularity_fade = 2.311f;
    value.template get<"source">().angle_rate = 0.027f;
    value.template get<"source">().complexity = 1.704f;
    value.template get<"source">().secondary_rate = 0.8f;
    value.template get<"source">().pattern_freq = 4.9755f;
    value.template get<"source">().pattern_mix = 0.0f;
    value.template get<"source">().speed = 0.04f;
    value.template get<"outer_warp">().scale = 1.0f;
    value.template get<"outer_warp">().speed = -0.00005f;
    value.template get<"outer_warp">().strength = 0.138f;
    value.template get<"outer_warp">().vector_angle = 0.0f;
    value.template get<"inner_warp">().cell_x = 1.0f;
    value.template get<"inner_warp">().cell_y = 1.0f;
    value.template get<"inner_warp">().offset_x = 0.0f;
    value.template get<"inner_warp">().offset_y = 0.0f;
    value.template get<"inner_warp">().rotation = 0.0f;
    value.template get<"inner_warp">().speed = 0.0032799998f;
    return value;
  }
  // clang-format on
  // End generated params.

public:
  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<
               W, H, KaleidoscopeStainedGlass<W, H>,
               KaleidoscopeStainedGlassSpec>::choreography_hooks()
        .template extend<
            KaleidoscopeStainedGlass,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &KaleidoscopeStainedGlass::initial_params>>();
  }
};
