/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file KaleidoscopeMelt.h
 * @brief Composed pullback effect:
 *        a dodecahedral lattice rippling through curl noise.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

namespace hs_test {
template <typename FX> struct ComposedWhiteBox;
} // namespace hs_test

/** @brief Pullback stage spec for `KaleidoscopeMelt`. */
struct KaleidoscopeMeltSpec : Pullback::Spec {
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::ANALOGOUS;
  static constexpr bool ANIMATED_PROJECTION = true;
  /**
   * @brief Stage chain: camera rotation, dodecahedral kaleidoscope lens,
   *  surface-noise displacement, projection, lattice source, generated
   *  palette.
   * @tparam B The effect's Binding.
   */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::OUT_OF_LINE_FLASH,
          Pullback::Stage::Lens<Pullback::Lens::DodecahedralKaleidoscope>,
          Pullback::Stage::Displace<typename Pullback::SurfacePolicyFor<
              Pullback::SurfaceNoiseParams, B, false>::Type>,
          Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
              Pullback::ProjectionKind::EQUIRECTANGULAR, B>::Type>>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::LatticeSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<
              Pullback::ProjectionCoverageMode::WEIGHT, B, void>::Type>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, Pullback::HueMode::NOISE,
                                  Pullback::Color::BrightnessEnvelope::NONE>>>>;
};
/// Parameter block derived from `KaleidoscopeMeltSpec`.
using KaleidoscopeMeltParams = Pullback::ParamsFor<KaleidoscopeMeltSpec>;

/**
 * @brief A dodecahedral lattice rippling through sphere-space curl noise.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopeMelt
    : public Pullback::ComposedEffect<W, H, KaleidoscopeMelt<W, H>,
                                      KaleidoscopeMeltSpec> {
  friend struct ::hs_test::ComposedWhiteBox<KaleidoscopeMelt<W, H>>;

public:
  /// Parameter block derived from `KaleidoscopeMeltSpec`.
  using Params = KaleidoscopeMeltParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::KaleidoscopeMelt::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "eb99693c0a7cf03b186d15806b8fb40244f09c85e1f03ca251e900d8d23cc680";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "7992482e13a3f7164652d58b46d82a1fc4e9a1407534c75d494d433b1bd36620";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "kaleidoscope-melt"
  };
  /// Frames each preset holds before advancing.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  /// Snapshot schema version; changes with the `Params` layout.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  /**
   * @brief Parameters of the first preset; every preset varies from them.
   * @return The preset-0 `Params`.
   */
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
    value.template get<"surface">().scale = 1.9651489f;
    value.template get<"surface">().speed = -0.00008951823f;
    value.template get<"surface">().strength = 0.15546875f;
    return value;
  }
  // clang-format on
  // End generated params.

public:
  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, KaleidoscopeMelt<W, H>,
                                    KaleidoscopeMeltSpec>::choreography_hooks()
        .template extend<
            KaleidoscopeMelt,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &KaleidoscopeMelt::initial_params>>();
  }
};
