/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file RiverStones.h
 * @brief Ripple bursts across a hexagonal lattice of soft cells.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

/** @brief Stage spec for `RiverStones`. */
struct RiverStonesSpec : Pullback::Spec {
  /// Samples animation at each column time.
  static constexpr auto TEMPORAL =
      Pullback::TemporalPolicy::fused(Pullback::PaletteCache::NONE);
  /// Generated palette harmony.
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  /// Spin and wander projection frame.
  static constexpr auto PROJECTION_MOTION =
      Pullback::ProjectionMotion::ANIMATED;

  /** @brief Camera, ripple, projected hex cells, palette and noise layers.
   * @tparam B Effect binding. */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Displace<typename Pullback::SurfacePolicyFor<
          Pullback::SphereAnimation::BurstRippleParams, B, false,
          "displace1">::Type>,
      Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
          Pullback::ProjectionKind::EQUIRECTANGULAR, B>::Type>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
                                  Pullback::FixedTessellationParams<
                                      Pattern::TessellationKind::HEXAGONAL,
                                      Pattern::TessellationFill::CELLS>,
                                  B>::Type,
                              Pullback::Weight::Projection,
                              Pullback::ProjectionCoverage::Weight>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::AuthoredColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulate", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::DIRECT>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulate1", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulate2", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE>>>>;

  using Layer =
      Pullback::Color::ModulationState; ///< Prepared modulation endpoints.
  using Palette = Pullback::Color::
      TemporalGeneratedPaletteState; ///< Prepared palette endpoints.

  /** @brief Fixed source and mechanic order of the authored color tail. */
  struct ColorColumn {
    const Palette &palette;         ///< Ordered palette endpoints.
    const Layer &hue, &lift, &mask; ///< Ordered modulation layers.
    bool same_palette;              ///< RGB and alpha entries match.
    bool blend = palette.fraction > 0.0f && palette.fraction < 1.0f;
    bool retained = palette.fraction <= 0.0f;
    const Pullback::Modulation::Prepared *noise =
        retained ? hue.previous : hue.current;
    const Pattern::SphereLutView *lift_cube = cube(lift);
    const Pattern::SphereLutView *mask_cube = cube(mask);
    bool authored =
        hue.enabled && lift.enabled && mask.enabled &&
        hue.mechanic.kind == ColorModulation::Mechanic::HUE_ROTATION &&
        lift.mechanic.kind ==
            ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT &&
        mask.mechanic.kind == ColorModulation::Mechanic::OPACITY &&
        hue.source == Pattern::SourceKind::NOISE && noise &&
        (!blend || hue.previous == hue.current) && lift_cube && mask_cube;

    /** @brief Selected scalar cube for the column.
     * @param layer Baked modulation layer.
     * @return Byte cube admitting interpolation at this column time. */
    HS_HOT_INLINE const Pattern::SphereLutView *cube(const Layer &layer) const {
      if (!layer.paired_cube_ready || layer.previous->cube.linear_data ||
          layer.current->cube.linear_data)
        return nullptr;
      return &(retained ? layer.previous : layer.current)->cube;
    }
  };

  /** @brief Samples a signed-byte modulation cube.
   * @param tap Cube face coordinates.
   * @param cube Prepared byte cube.
   * @return Bounded signed scalar. */
  HS_HOT_FLASH_MEMBER static float
  sample_cube(const Pattern::SphereLutSample &tap,
              const Pattern::SphereLutView &cube) {
    return Pattern::bounded_scalar(tap.sample_data(cube.data));
  }

  /** @brief Bilinear retained value, sharing the current result for equal taps.
   * @param tap Cube face coordinates.
   * @param a Retained byte cube.
   * @param b Current byte cube.
   * @param current Current sample.
   * @return Retained signed scalar. */
  HS_HOT_INLINE static float
  previous_scalar(const Pattern::SphereLutSample &tap,
                  const Pattern::SphereLutView &a,
                  const Pattern::SphereLutView &b, float current) {
    const auto *old = a.data + tap.offset;
    const auto *now = b.data + tap.offset;
    constexpr int STRIDE = Pattern::SphereLutView::FACE_STEPS;
    if (std::memcmp(old, now, 2) == 0 &&
        std::memcmp(old + STRIDE, now + STRIDE, 2) == 0)
      return current;
    return sample_cube(tap, a);
  }

  /** @brief Shades edited layer orders through the canonical endpoint kernels.
   * @param field Cell and coverage.
   * @param column Prepared color layers.
   * @return Straight-alpha color. */
  HS_FLASH_MEMBER static Color4 shade_edited(const Pullback::FieldSample &field,
                                             const ColorColumn &column) {
    using namespace Pullback;
    const auto &palette = column.palette;
    const auto &state = column.retained ? palette.previous : palette.current;
    auto result = Color::apply_generated_palette(field, state);
    if (column.blend) {
      auto old = result;
      const auto &previous = palette.previous;
      old.color = previous.palette->get(result.palette_coordinate);
      old.color.alpha *= Color::palette_opacity(field, state);
      old.hue_rotation =
          previous.hue_rotation.active ? previous.hue_rotation.data : nullptr;
      old.hue_rotation_gamut = previous.hue_rotation.gamut_mask;
      for (const auto *layer : {&column.hue, &column.lift, &column.mask})
        Color::apply_modulation_pair_in_place(old, result, *layer);
      blend_color_samples_in_place(old, result, palette.fraction);
      return finish_color_sample(result);
    }
    for (const auto *layer : {&column.hue, &column.lift, &column.mask})
      Color::apply_modulation_in_place(result, *layer, column.retained);
    return finish_color_sample(result);
  }

  /** @brief Shades hue, lightness and opacity with one scalar cube address.
   * @param field Projected cell and coverage.
   * @param column Prepared color layers.
   * @return Straight-alpha color at the column time. */
  HS_HOT_FLASH_MEMBER static Color4
  shade_field(const Pullback::FieldSample &field, const ColorColumn &column) {
    using namespace Pullback;
    if (!column.authored)
      return shade_edited(field, column);
    const auto &palette = column.palette;
    const auto &base = column.retained ? palette.previous : palette.current;
    auto result = Color::apply_generated_palette(field, base);
    const auto direction = Pattern::SphereDirection::admit(field.sphere);
    result.pending_hue = Color::accumulate_hue(
        0.0f, Modulation::sample(*column.noise, direction),
        column.hue.mechanic);
    const auto tap = Pattern::prepare_sphere_lut_sample(direction);
    const float lift = sample_cube(tap, *column.lift_cube);
    const float mask = sample_cube(tap, *column.mask_cube);
    result.pending_lift =
        Color::accumulate_lift(0.0f, lift, column.lift.mechanic);
    const auto opacity = [&](Color4 color, float value) HS_HOT_INLINE {
      return ColorModulation::Opacity::apply(color, value,
                                             {column.mask.mechanic.amount})
          .alpha;
    };
    const auto current = result.color;
    result.color.alpha = opacity(current, mask);
    if (column.blend) {
      auto old = current;
      if (!column.same_palette) {
        old = palette.previous.palette->get(result.palette_coordinate);
        old.alpha *= Color::palette_opacity(field, base);
        result.color.color =
            old.color.lerp16(current.color, frac_to_q16(palette.fraction));
      }
      const auto *table = palette.previous.hue_rotation.active
                              ? palette.previous.hue_rotation.data
                              : nullptr;
      if (old.color != current.color || table != result.hue_rotation ||
          palette.previous.hue_rotation.gamut_mask !=
              result.hue_rotation_gamut) {
        result.hue_rotation = nullptr;
        result.hue_rotation_gamut = nullptr;
      }
      const float old_lift = previous_scalar(tap, column.lift.previous->cube,
                                             *column.lift_cube, lift);
      const float old_mask = previous_scalar(tap, column.mask.previous->cube,
                                             *column.mask_cube, mask);
      result.pending_lift =
          hs::lerp(Color::accumulate_lift(0.0f, old_lift, column.lift.mechanic),
                   result.pending_lift, palette.fraction);
      result.color.alpha = hs::lerp(opacity(old, old_mask), result.color.alpha,
                                    palette.fraction);
    }
    return finish_color_sample(result);
  }
};
/// Parameter block derived from `RiverStonesSpec`.
using RiverStonesParams = Pullback::ParamsFor<RiverStonesSpec>;

/**
 * @brief Ripple bursts across a hexagonal lattice of soft cells.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class RiverStones : public Pullback::ComposedEffect<W, H, RiverStones<W, H>,
                                                    RiverStonesSpec> {
  using Base =
      Pullback::ComposedEffect<W, H, RiverStones<W, H>, RiverStonesSpec>;

public:
  using Params = RiverStonesParams; ///< Live parameter struct.
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::RiverStones::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "d3877f833447f1bc6e04ff8a514a15ec13bec8ca5dcde6ab8d41fc8d21235e8f";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "74c75d2cf6e63178f35d56c193d8bb9533983bb4f0a044112221f9dafdf8557e";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "river-stones"
  };
  /// Frames each preset holds before advancing.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  /// Outer camera spin about +Y, in radians per frame.
  static constexpr float CAMERA_SPIN_RATE = 0.0f;
  // clang-format on
  // End generated identity.
  /// Params layout version; snapshots tagged with another are rejected.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  /**
   * @brief Parameters of the first preset; every preset varies from them.
   * @return The preset-0 `Params`.
   */
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 0.1f;
    value.template get<"color">().brightness_bottom = 0.0f;
    value.template get<"color">().brightness_envelope = Pullback::Color::BrightnessEnvelope::NONE;
    value.template get<"color">().brightness_top = 1.0f;
    value.template get<"color">().gradient_end_b = 1.0f;
    value.template get<"color">().gradient_end_g = 0.2f;
    value.template get<"color">().gradient_end_r = 0.1f;
    value.template get<"color">().gradient_middle_b = 0.4f;
    value.template get<"color">().gradient_middle_g = 1.0f;
    value.template get<"color">().gradient_middle_r = 0.1f;
    value.template get<"color">().gradient_start_b = 0.0f;
    value.template get<"color">().gradient_start_g = 0.15f;
    value.template get<"color">().gradient_start_r = 1.0f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 0.7923977f;
    value.template get<"color">().palette_harmony = 5;
    value.template get<"color">().palette_hue = 0.0f;
    value.template get<"color">().palette_lightness = 0.62f;
    value.template get<"color">().palette_lightness_curve = 0;
    value.template get<"color">().palette_lightness_range = 0.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.template get<"color">().palette_mode = 0;
    value.template get<"color">().palette_mutation = 0.0f;
    value.template get<"color">().palette_name = 17;
    value.template get<"color">().palette_name_b = 19;
    value.template get<"color">().palette_spread = 0.07f;
    value.template get<"color">().palette_type = 2;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"displace1">().burst_rate = 0.0078125f;
    value.template get<"displace1">().burst_size = 3;
    value.template get<"displace1">().center_azimuth = 0.0f;
    value.template get<"displace1">().center_polar = 1.5707964f;
    value.template get<"displace1">().decay = 0.0f;
    value.template get<"displace1">().period = 97.46764f;
    value.template get<"displace1">().random_center = 1;
    value.template get<"displace1">().strength = 0.15f;
    value.template get<"displace1">().thickness = 0.7f;
    value.template get<"modulate">().modulation.grid.axis_b.x = 1.0f;
    value.template get<"modulate">().modulation.grid.axis_b.y = 0.0f;
    value.template get<"modulate">().modulation.grid.axis_b.z = 0.0f;
    value.template get<"modulate">().modulation.grid.axis.x = 0.0f;
    value.template get<"modulate">().modulation.grid.axis.y = 1.0f;
    value.template get<"modulate">().modulation.grid.axis.z = 0.0f;
    value.template get<"modulate">().modulation.grid.frequency = 1.0f;
    value.template get<"modulate">().modulation.grid.frequency_b = 1.0f;
    value.template get<"modulate">().modulation.grid.offset = 0.0f;
    value.template get<"modulate">().modulation.grid.offset_b = 0.0f;
    value.template get<"modulate">().modulation.grid.speed = 0.0f;
    value.template get<"modulate">().modulation.interference.axis_b.x = 1.0f;
    value.template get<"modulate">().modulation.interference.axis_b.y = 0.0f;
    value.template get<"modulate">().modulation.interference.axis_b.z = 0.0f;
    value.template get<"modulate">().modulation.interference.axis.x = 0.0f;
    value.template get<"modulate">().modulation.interference.axis.y = 1.0f;
    value.template get<"modulate">().modulation.interference.axis.z = 0.0f;
    value.template get<"modulate">().modulation.interference.frequency = 1.0f;
    value.template get<"modulate">().modulation.interference.frequency_b = 1.0f;
    value.template get<"modulate">().modulation.interference.offset = 0.0f;
    value.template get<"modulate">().modulation.interference.offset_b = 0.0f;
    value.template get<"modulate">().modulation.interference.speed = 0.0f;
    value.template get<"modulate">().modulation.mechanic.amount = 0.25882354f;
    value.template get<"modulate">().modulation.enabled = true;
    value.template get<"modulate">().modulation.evaluation = Pattern::Evaluation::DIRECT;
    value.template get<"modulate">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulate">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulate">().modulation.noise.contrast = 0.0f;
    value.template get<"modulate">().modulation.noise.basis = Pattern::NoiseBasis::SIMPLEX;
    value.template get<"modulate">().modulation.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
    value.template get<"modulate">().modulation.noise.scale = 1.0f;
    value.template get<"modulate">().modulation.noise.speed = 0.0f;
    value.template get<"modulate">().modulation.sine.axis.x = 0.0f;
    value.template get<"modulate">().modulation.sine.axis.y = 1.0f;
    value.template get<"modulate">().modulation.sine.axis.z = 0.0f;
    value.template get<"modulate">().modulation.sine.frequency = 1.0f;
    value.template get<"modulate">().modulation.sine.offset = 0.0f;
    value.template get<"modulate">().modulation.sine.speed = 0.052941177f;
    value.template get<"modulate1">().modulation.grid.axis_b.x = 1.0f;
    value.template get<"modulate1">().modulation.grid.axis_b.y = 0.0f;
    value.template get<"modulate1">().modulation.grid.axis_b.z = 0.0f;
    value.template get<"modulate1">().modulation.grid.axis.x = 0.0f;
    value.template get<"modulate1">().modulation.grid.axis.y = 1.0f;
    value.template get<"modulate1">().modulation.grid.axis.z = 0.0f;
    value.template get<"modulate1">().modulation.grid.frequency = 1.0f;
    value.template get<"modulate1">().modulation.grid.frequency_b = 1.0f;
    value.template get<"modulate1">().modulation.grid.offset = 0.0f;
    value.template get<"modulate1">().modulation.grid.offset_b = 0.0f;
    value.template get<"modulate1">().modulation.grid.speed = 0.0f;
    value.template get<"modulate1">().modulation.interference.axis_b.x = 1.0f;
    value.template get<"modulate1">().modulation.interference.axis_b.y = 0.0f;
    value.template get<"modulate1">().modulation.interference.axis_b.z = 0.0f;
    value.template get<"modulate1">().modulation.interference.axis.x = 0.0f;
    value.template get<"modulate1">().modulation.interference.axis.y = 1.0f;
    value.template get<"modulate1">().modulation.interference.axis.z = 0.0f;
    value.template get<"modulate1">().modulation.interference.frequency = 1.0f;
    value.template get<"modulate1">().modulation.interference.frequency_b = 1.0f;
    value.template get<"modulate1">().modulation.interference.offset = 0.0f;
    value.template get<"modulate1">().modulation.interference.offset_b = 0.0f;
    value.template get<"modulate1">().modulation.interference.speed = 0.0f;
    value.template get<"modulate1">().modulation.mechanic.amount = 0.4617647f;
    value.template get<"modulate1">().modulation.enabled = true;
    value.template get<"modulate1">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"modulate1">().modulation.mechanic.kind = ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT;
    value.template get<"modulate1">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulate1">().modulation.noise.contrast = 0.0f;
    value.template get<"modulate1">().modulation.noise.basis = Pattern::NoiseBasis::SIMPLEX;
    value.template get<"modulate1">().modulation.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
    value.template get<"modulate1">().modulation.noise.scale = 1.8775138f;
    value.template get<"modulate1">().modulation.noise.speed = 0.0001f;
    value.template get<"modulate1">().modulation.sine.axis.x = 0.0f;
    value.template get<"modulate1">().modulation.sine.axis.y = 1.0f;
    value.template get<"modulate1">().modulation.sine.axis.z = 0.0f;
    value.template get<"modulate1">().modulation.sine.frequency = 1.0f;
    value.template get<"modulate1">().modulation.sine.offset = 0.0f;
    value.template get<"modulate1">().modulation.sine.speed = 0.0f;
    value.template get<"modulate2">().modulation.grid.axis_b.x = 1.0f;
    value.template get<"modulate2">().modulation.grid.axis_b.y = 0.0f;
    value.template get<"modulate2">().modulation.grid.axis_b.z = 0.0f;
    value.template get<"modulate2">().modulation.grid.axis.x = 0.0f;
    value.template get<"modulate2">().modulation.grid.axis.y = 1.0f;
    value.template get<"modulate2">().modulation.grid.axis.z = 0.0f;
    value.template get<"modulate2">().modulation.grid.frequency = 1.0f;
    value.template get<"modulate2">().modulation.grid.frequency_b = 1.0f;
    value.template get<"modulate2">().modulation.grid.offset = 0.0f;
    value.template get<"modulate2">().modulation.grid.offset_b = 0.0f;
    value.template get<"modulate2">().modulation.grid.speed = 0.0f;
    value.template get<"modulate2">().modulation.interference.axis_b.x = 1.0f;
    value.template get<"modulate2">().modulation.interference.axis_b.y = 0.0f;
    value.template get<"modulate2">().modulation.interference.axis_b.z = 0.0f;
    value.template get<"modulate2">().modulation.interference.axis.x = 0.0f;
    value.template get<"modulate2">().modulation.interference.axis.y = 1.0f;
    value.template get<"modulate2">().modulation.interference.axis.z = 0.0f;
    value.template get<"modulate2">().modulation.interference.frequency = 1.0f;
    value.template get<"modulate2">().modulation.interference.frequency_b = 1.0f;
    value.template get<"modulate2">().modulation.interference.offset = 0.0f;
    value.template get<"modulate2">().modulation.interference.offset_b = 0.0f;
    value.template get<"modulate2">().modulation.interference.speed = 0.0f;
    value.template get<"modulate2">().modulation.mechanic.amount = 0.6382353f;
    value.template get<"modulate2">().modulation.enabled = true;
    value.template get<"modulate2">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"modulate2">().modulation.mechanic.kind = ColorModulation::Mechanic::OPACITY;
    value.template get<"modulate2">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulate2">().modulation.noise.contrast = 0.0f;
    value.template get<"modulate2">().modulation.noise.basis = Pattern::NoiseBasis::SIMPLEX;
    value.template get<"modulate2">().modulation.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
    value.template get<"modulate2">().modulation.noise.scale = 1.5627832f;
    value.template get<"modulate2">().modulation.noise.speed = 0.0001f;
    value.template get<"modulate2">().modulation.sine.axis.x = 0.0f;
    value.template get<"modulate2">().modulation.sine.axis.y = 1.0f;
    value.template get<"modulate2">().modulation.sine.axis.z = 0.0f;
    value.template get<"modulate2">().modulation.sine.frequency = 1.0f;
    value.template get<"modulate2">().modulation.sine.offset = 0.0f;
    value.template get<"modulate2">().modulation.sine.speed = 0.0f;
    value.template get<"projection">().central_meridian = 0.0f;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 0.0f;
    value.template get<"projection">().singularity_fade = 1.0f;
    value.template get<"source">().cell_radius = 0.41925988f;
    value.template get<"source">().cell_scale = 2.5352101f;
    value.template get<"source">().cell_shape = 0.5292398f;
    value.template get<"source">().cell_softness = 0.10124387f;
    return value;
  }
  // clang-format on
  // End generated params.

  /** @brief Binds the palette and ordered layers for one column.
   * @param frame Prepared column frame.
   * @return Immutable color evaluation context. */
  static RiverStonesSpec::ColorColumn
  color_column(const typename Base::Frame &frame) {
    return {std::get<4>(frame.prepared), std::get<5>(frame.prepared),
            std::get<6>(frame.prepared), std::get<7>(frame.prepared),
            frame.ctx.temporal.same_palette};
  }

  /** @brief Binds the shared color context for the latitude scan.
   * @param frame Prepared column.
   * @param x Column coordinate.
   * @return Shader borrowing the frame. */
  static auto prepare_column_shader(const typename Base::Frame &frame, int x) {
    (void)x;
    return [&frame, column = color_column(frame)](const math::Vector &view) {
      return shade_column(view, frame, column);
    };
  }

  /** @brief Shades a unit view direction at the prepared column time.
   * @param view View direction.
   * @param frame Prepared column.
   * @return Straight-alpha pixel color. */
  HS_HOT_FLASH_MEMBER static Color4 shade(const math::Vector &view,
                                          const typename Base::Frame &frame) {
    return shade_column(view, frame, color_column(frame));
  }

  /** @brief Shades a pixel with its prebound color layers.
   * @param view View direction.
   * @param frame Prepared geometry.
   * @param column Prepared color layers.
   * @return Straight-alpha pixel color. */
  HS_HOT_FLASH_MEMBER static Color4
  shade_column(const math::Vector &view, const typename Base::Frame &frame,
               const RiverStonesSpec::ColorColumn &column) {
    using Pipeline = typename Base::RenderPipeline;
    static_assert(Pipeline::NODE_COUNT == 8);
    const auto run = [&]<size_t I>(const auto &input) HS_HOT_INLINE {
      return Pipeline::template NodeAt<I>::run(input, frame.ctx,
                                               std::get<I>(frame.prepared));
    };
    const auto rotated =
        run.template operator()<0>(Pullback::SphereSample{view, 0.0f});
    const auto displaced = run.template operator()<1>(rotated);
    const auto projected = run.template operator()<2>(displaced);
    const auto field = run.template operator()<3>(projected);
    return RiverStonesSpec::shade_field(field, column);
  }

  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, RiverStones<W, H>,
                                    RiverStonesSpec>::choreography_hooks()
        .template extend<
            RiverStones,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &RiverStones::initial_params>>();
  }
};
