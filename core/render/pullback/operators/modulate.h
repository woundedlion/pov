/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include <memory>

#include "core/render/pullback/operators.h"

/** @file modulate.h
 * @brief Ordered color endomorphisms with instance-owned sphere fields. */

namespace Pullback::Interp::Op {

/** @brief Sphere-field controls shared with legacy palette recipes. */
struct FieldModulationParams : GeneratedPaletteParams {
  static constexpr std::array<Field<FieldModulationParams>, 0> FIELDS{};
  static constexpr std::array<TopologyField<FieldModulationParams>, 0>
      TOPOLOGY{};

  /** @brief Samples mechanic strength while source recipes hold endpoints.
   * @param out Column controls.
   * @param a Previous controls.
   * @param b Current controls.
   * @param fraction Current endpoint weight. */
  static void sample_extra(FieldModulationParams &out,
                           const FieldModulationParams &a,
                           const FieldModulationParams &b, float fraction) {
    out.modulation = b.modulation;
    if (Modulation::same_topology(a.modulation, b.modulation))
      out.modulation.mechanic.amount = hs::lerp(
          a.modulation.mechanic.amount, b.modulation.mechanic.amount, fraction);
  }
};

/** @brief Source clock and separately owned endpoint cache. */
struct FieldModulationState {
  ColorClockState clock; ///< Serializable source clock.
  mutable std::unique_ptr<Modulation::Resources<>> resources;
  ///< Cache owned by this state, independently of migrated and retained states.
};

/** @brief COLOR endomorphism driven by an independent sphere field. */
struct ModulateField {
  static constexpr const char *ID = "modulate.field.v1";
  static constexpr const char *NAME = "Field";
  using Input = Color4;
  using Output = Color4;
  using Params = FieldModulationParams;
  using State = FieldModulationState;
  /** @brief Immutable endpoint and sampled mechanic. */
  struct Prepared {
    const Modulation::Prepared *source;   ///< Instance-owned source endpoint.
    ColorModulation::AmountParams amount; ///< Column mechanic strength.
    ColorModulation::Kernel mechanic;     ///< Selected color transform.
  };

  static void init(State &state, InstanceId) {
    state.resources = std::make_unique<Modulation::Resources<>>();
  }
  static Status migrate(State &destination, const State &source,
                        InstanceId id) {
    init(destination, id);
    destination.clock = source.clock;
    return Status::OK;
  }
  static const char *validate(const Params &params) {
    return Modulation::admitted(params.modulation)
               ? nullptr
               : Modulation::INVALID_CONFIGURATION;
  }
  static constexpr const char *ADMISSIBILITY_CONVEXITY =
      "Source recipes hold endpoints; mechanic strength interpolates within matching topology.";
  static constexpr bool APPROXIMATE = true;
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS;
  static constexpr auto METRICS = Color::GENERATED_PALETTE_METRICS;

  static void advance(State &state, const Params &params) {
    Modulation::advance(state.clock.modulation, params.modulation);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    return {
        &state.resources->prepare(params.modulation, state.clock.modulation),
        {params.modulation.mechanic.amount},
        ColorModulation::kernel(params.modulation.mechanic.kind)};
  }
  static Color4 run(const Color4 &input, const FrameContext &ctx,
                    const Params &, const Prepared &prepared) {
    if (!prepared.source->config.enabled || prepared.amount.amount == 0.0f)
      return input;
    HS_CHECK(ctx.color_sample != nullptr,
             "field modulation requires field sample metadata");
    const auto DIRECTION =
        Pattern::SphereDirection::admit(ctx.color_sample->sphere);
    return prepared.mechanic(input,
                             Modulation::sample(*prepared.source, DIRECTION),
                             prepared.amount);
  }
  static Color4 run_temporal(const Color4 &input, const FrameContext &ctx,
                             const Prepared &current, const Prepared &previous,
                             float fraction) {
    if (!current.source->config.enabled || current.amount.amount == 0.0f)
      return input;
    HS_CHECK(ctx.color_sample != nullptr,
             "field modulation requires field sample metadata");
    const auto DIRECTION =
        Pattern::SphereDirection::admit(ctx.color_sample->sphere);
    if (!Modulation::same_topology(current.source->config,
                                   previous.source->config) ||
        Modulation::same_field(*current.source, *previous.source) ||
        fraction >= 1.0f)
      return current.mechanic(input,
                              Modulation::sample(*current.source, DIRECTION),
                              current.amount);
    const Color4 BEFORE = current.mechanic(
        input, Modulation::sample(*previous.source, DIRECTION), current.amount);
    if (fraction <= 0.0f)
      return BEFORE;
    const Color4 NOW = current.mechanic(
        input, Modulation::sample(*current.source, DIRECTION), current.amount);
    return Modulation::interpolate_colors(BEFORE, NOW, fraction,
                                          frac_to_q16(fraction));
  }
  static constexpr bool temporal_hold_field(uint16_t) { return true; }
};

/** @brief Path-length hue rotation controls. */
struct HueShiftModulationParams {
  float hue_shift_amount = 0.0f; ///< Hue turns per unit path length.
  uint8_t hue_mode = static_cast<uint8_t>(HueShiftMode::PATH_LENGTH);
  ///< Selected path-length rotation mode.
  static constexpr auto FIELDS = std::array{
      Field<HueShiftModulationParams>{
          "hue-shift-amount", &HueShiftModulationParams::hue_shift_amount,
          "Amount", -4.0f, 4.0f, FieldCurve::LERP, FieldGate::ALWAYS,
          Color::HUE_ROTATION_GATE},
  };
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<HueShiftModulationParams>{
          "hue-shift-mode", &HueShiftModulationParams::hue_mode,
          HUE_SHIFT_MODE_IDS, static_cast<uint8_t>(HueShiftMode::PATH_LENGTH)},
  };
};

/** @brief COLOR endomorphism rotating hue along accumulated path length. */
struct ModulateHueShift : StatelessModel {
  static constexpr const char *ID = "modulate.hue-shift.v1";
  static constexpr const char *NAME = "Hue Shift";
  using Input = Color4;
  using Output = Color4;
  using Params = HueShiftModulationParams;
  struct Prepared {};

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &) {
    HS_CHECK(params.hue_mode <= static_cast<uint8_t>(HueShiftMode::PATH_LENGTH),
             "hue modulation: invalid mode");
    return {};
  }
  static Color4 run(const Color4 &input, const FrameContext &ctx,
                    const Params &params, const Prepared &) {
    if (params.hue_mode == static_cast<uint8_t>(HueShiftMode::NONE) ||
        params.hue_shift_amount == 0.0f)
      return input;
    HS_CHECK(ctx.color_sample != nullptr,
             "hue modulation requires field sample metadata");
    return ColorModulation::HueRotation::apply(
        input, ctx.color_sample->path_length, {params.hue_shift_amount});
  }
};

} // namespace Pullback::Interp::Op

#endif // HS_ENABLE_CHAIN_INTERPRETER
