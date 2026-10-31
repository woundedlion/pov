/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file Millefiori.h
 * @brief Composed pullback effect:
 *        layered kaleidoscopic lattice flowers with ripple bursts.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

/** @brief Tetrahedral folding of the preceding dodecahedral chamber. */
struct MillefioriTetrahedralLens : Pullback::Lens::TetrahedralKaleidoscope {
  /**
   * @brief Folds a dodecahedral chamber direction into the tetrahedral chamber.
   * @tparam FrameState Frame state type; unused.
   * @param input Unit direction in the dodecahedral chamber.
   * @return Direction in the tetrahedral chamber.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static math::Vector
  apply(const math::Vector &input, const FrameState &) {
    return lenses::tetrahedral_from_dodecahedral_chamber(input);
  }
};

/** @brief Pullback stage spec for `Millefiori`. */
struct MillefioriSpec : Pullback::Spec {
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::ANALOGOUS;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr bool TEMPORAL_SAMPLING = true;
  static constexpr bool TEMPORAL_COLOR_PHASE = true;
  static constexpr bool PALETTE_LAB_CACHE = true;
  static constexpr bool CACHE_FULL_GAMUT_LUT = true;
  static constexpr bool RIPPLE_RADIAL_TABLE = true;
  static constexpr bool PERSISTENT_PREPARED_FRAME = true;
  /** @brief Authored ripple, kaleidoscope and lattice stage chain.
   * @tparam B Frame binding. */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Displace<
              Pullback::Surface::BurstRipple<Pullback::SurfaceProvider<
                  B, Pullback::SphereAnimation::BurstRippleParams, true,
                  "displace">>>,
          Pullback::Stage::Lens<Pullback::Lens::DodecahedralKaleidoscope>,
          Pullback::Stage::Lens<MillefioriTetrahedralLens>,
          Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
              Pullback::ProjectionKind::GNOMONIC_FOLDED, B>::Type>>,
      Pullback::Stage::Warp<Pullback::Warp::Vortex<Pullback::WarpProvider<
          B, "outer_warp", Pullback::Warp::VortexParams, true>>>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
                                  Pullback::LatticeSourceParams, B>::Type,
                              Pullback::Weight::Projection,
                              Pullback::ProjectionCoverage::Weight>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, Pullback::Color::BrightnessEnvelope::BELL,
                                  PALETTE_LAB_CACHE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "field", Pattern::SourceKind::SPHERE_SMOOTH_GRID,
              Pattern::Evaluation::DIRECT, TEMPORAL_COLOR_PHASE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulate", Pattern::SourceKind::DISPLACEMENT,
              Pattern::Evaluation::DIRECT, TEMPORAL_COLOR_PHASE>>>>;
};
/// Parameter block derived from `MillefioriSpec`.
using MillefioriParams = Pullback::ParamsFor<MillefioriSpec>;

/**
 * @brief Layered kaleidoscopic lattice flowers with ripple bursts.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class Millefiori
    : public Pullback::ComposedEffect<W, H, Millefiori<W, H>, MillefioriSpec> {
  using Base = Pullback::ComposedEffect<W, H, Millefiori<W, H>, MillefioriSpec>;

public:
  /// Parameter block derived from `MillefioriSpec`.
  using Params = MillefioriParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::Millefiori::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "7ae5f5d0651ce6445f8416ba726fc5ee85fa4cdf953b6aef70e4f65cc4584e8a";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "5de232336916f545f7699dd5cb1a92ce7c24dd530a79b6352ee986ebf71b6331";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "affine-contour"
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
    value.template get<"color">().brightness_bottom = 0.0f;
    value.template get<"color">().brightness_top = 1.0f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = -0.046783626f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::BELL;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"displace">().burst_rate = 0.0078125f;
    value.template get<"displace">().center_azimuth = 0.0f;
    value.template get<"displace">().center_polar = 1.5707964f;
    value.template get<"displace">().decay = 0.2f;
    value.template get<"displace">().period = 91.48529f;
    value.template get<"displace">().strength = 0.15f;
    value.template get<"displace">().thickness = 0.7f;
    value.template get<"field">().modulation.grid.axis_b.x = 1.0f;
    value.template get<"field">().modulation.grid.axis_b.y = 0.0f;
    value.template get<"field">().modulation.grid.axis_b.z = 0.0f;
    value.template get<"field">().modulation.grid.axis.x = 0.0f;
    value.template get<"field">().modulation.grid.axis.y = 23.204678f;
    value.template get<"field">().modulation.grid.axis.z = 37.4269f;
    value.template get<"field">().modulation.grid.frequency = 3.0f;
    value.template get<"field">().modulation.grid.frequency_b = 3.0f;
    value.template get<"field">().modulation.grid.offset = 0.0f;
    value.template get<"field">().modulation.grid.offset_b = 0.0f;
    value.template get<"field">().modulation.grid.speed = 0.02f;
    value.template get<"field">().modulation.interference.axis_b.x = 1.0f;
    value.template get<"field">().modulation.interference.axis_b.y = 1.497076f;
    value.template get<"field">().modulation.interference.axis_b.z = 40.42105f;
    value.template get<"field">().modulation.interference.axis.x = 1.0f;
    value.template get<"field">().modulation.interference.axis.y = 1.0f;
    value.template get<"field">().modulation.interference.axis.z = 0.0f;
    value.template get<"field">().modulation.interference.frequency = 3.9298246f;
    value.template get<"field">().modulation.interference.frequency_b = 11.415205f;
    value.template get<"field">().modulation.interference.offset = 0.0f;
    value.template get<"field">().modulation.interference.offset_b = -12.725146f;
    value.template get<"field">().modulation.interference.speed = 0.02f;
    value.template get<"field">().modulation.mechanic.amount = 0.9035088f;
    value.template get<"field">().modulation.enabled = true;
    value.template get<"field">().modulation.evaluation = Pattern::Evaluation::DIRECT;
    value.template get<"field">().modulation.mechanic.kind = ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT;
    value.template get<"field">().modulation.source = Pattern::SourceKind::SPHERE_SMOOTH_GRID;
    value.template get<"field">().modulation.noise.contrast = 0.0f;
    value.template get<"field">().modulation.noise.basis = Pattern::NoiseBasis::SIMPLEX;
    value.template get<"field">().modulation.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
    value.template get<"field">().modulation.noise.scale = 4.6928453f;
    value.template get<"field">().modulation.noise.speed = -0.0003625731f;
    value.template get<"field">().modulation.sine.axis.x = 8.982456f;
    value.template get<"field">().modulation.sine.axis.y = 31.438597f;
    value.template get<"field">().modulation.sine.axis.z = 25.450293f;
    value.template get<"field">().modulation.sine.frequency = 5.0526314f;
    value.template get<"field">().modulation.sine.offset = -17.964912f;
    value.template get<"field">().modulation.sine.speed = 0.1754386f;
    value.template get<"modulate">().modulation.mechanic.amount = 2.235294f;
    value.template get<"modulate">().modulation.enabled = true;
    value.template get<"modulate">().modulation.evaluation = Pattern::Evaluation::DIRECT;
    value.template get<"modulate">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulate">().modulation.source = Pattern::SourceKind::DISPLACEMENT;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 0.0f;
    value.template get<"projection">().singularity_fade = 1.0f;
    value.template get<"source">().lattice_cell_scale = 2.1401224f;
    value.template get<"source">().lattice_radius = 0.33298188f;
    value.template get<"source">().lattice_shape_blend = 1.0f;
    value.template get<"source">().lattice_softness = 0.1608203f;
    value.template get<"outer_warp">().center_orbit_radius = 0.0f;
    value.template get<"outer_warp">().center_x = -3.2705882f;
    value.template get<"outer_warp">().center_y = 0.49411765f;
    value.template get<"outer_warp">().radius = 1.3072151f;
    value.template get<"outer_warp">().speed = 0.0f;
    value.template get<"outer_warp">().turns = 0.7764706f;
    return value;
  }
  // clang-format on
  // End generated params.

public:
  /** @brief Evaluates the composed pixel pipeline from ITCM.
   * @param view Unit view direction.
   * @param frame Prepared column state.
   * @return Straight-alpha pixel color. */
  HS_HOT_ITCM_MEMBER static Color4 shade(const math::Vector &view,
                                         const typename Base::Frame &frame) {
    return Base::RenderPipeline::evaluate(view, frame.ctx, frame.prepared);
  }

  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, Millefiori<W, H>,
                                    MillefioriSpec>::choreography_hooks()
        .template extend<
            Millefiori,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &Millefiori::initial_params>>();
  }
};
