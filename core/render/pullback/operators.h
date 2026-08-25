/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "render/pullback/material.h"
#include "render/pullback/operator_model.h"
#include "render/pullback/operators_common.h"
#include "render/pullback/operators_field.h"
#include "render/pullback/operators_project.h"
#include "render/pullback/operators_sample.h"
#include "render/pullback/operators_sphere.h"
#include "render/pullback/operators_warp.h"
#include "render/pullback/projection.h"
#include "render/pullback/source.h"
#include "render/pullback/stage.h"

/**
 * @file operators.h
 * @brief Concrete chain-interpreter operator models: erased adapters over the
 *        shared carrier kernels, with param-block-backed policy math. The
 *        per-family headers carry the family batches; this header aggregates
 *        them and keeps the chain's entry and exit operators.
 */

namespace Pullback {

namespace Interp {

namespace Op {

/** @brief Parameter family of sphere.rotate.v2. */
struct RotateChainParams {
  float wander = 0.0f;    /**< Fraction of the walk delta absorbed per frame. */
  float spin_rate = 0.0f; /**< Per-frame spin about Y, in radians. */

  static constexpr auto FIELDS = std::array{
      Field<RotateChainParams>{"wander", &RotateChainParams::wander, "Wander",
                               0.0f, 1.0f, FieldCurve::LERP},
      Field<RotateChainParams>{"spin-speed", &RotateChainParams::spin_rate,
                               "Spin Speed", 0.0f, 0.05f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<RotateChainParams>());

/** @brief SPHERE endomorphism: the wandering, spinning camera. */
struct Rotate {
  static constexpr const char *ID = "sphere.rotate.v2";
  static constexpr const char *NAME = "Rotate";
  using Input = SphereSample;
  using Output = SphereSample;
  using Params = RotateChainParams;
  using State = SpatialWalkState;
  struct Prepared {
    Quaternion conjugate;
  };

  static void init(State &state, InstanceId id) {
    init_walk(state, static_cast<int32_t>(id.stable_hash));
  }
  static Status migrate(State &dst, const State &src, InstanceId) {
    dst = src;
    return Status::OK;
  }
  static void advance(State &state, const Params &params) {
    advance_walk(state, params.wander, params.spin_rate);
  }
  static Prepared prepare(const FrameContext &, const Params &,
                          const State &state) {
    return {
        (make_rotation(Y_AXIS, state.spin_phase) * state.wander).conjugate()};
  }
  static SphereSample run(const SphereSample &input, const FrameContext &,
                          const Params &, const Prepared &prepared) {
    return Kernel::rotate_dir(input, prepared.conjugate);
  }
};

enum class PaletteMode : uint8_t {
  TRIADIC = 0,
  COMPLEMENTARY = 1,
  ANALOGOUS = 2
};
using HueShiftMode = Color::HueMode;
enum class EnvelopeMode : uint8_t {
  NONE = 0,
  CUP = 1,
  BELL = 2,
  ASCENDING = 3,
  DESCENDING = 4
};

inline constexpr const char *PALETTE_MODE_IDS[] = {"triadic", "complementary",
                                                   "analogous"};
inline constexpr const char *PALETTE_MAPPING_IDS[] = {"cup", "bell", "linear",
                                                      "reverse"};
inline constexpr const char *HUE_SHIFT_MODE_IDS[] = {"none", "noise",
                                                     "path-length"};
inline constexpr const char *BRIGHTNESS_ENVELOPE_IDS[] = {
    "none", "cup", "bell", "ascending", "descending"};

/** @brief Parameter family of colorize.generated-palette.v2.
    @details The mapping topology enum8 supersedes the base family's
    `palette_mapping` member, which the chain never reads. */
struct GeneratedPaletteParams : Color::ColorParams {
  uint8_t palette_mode = static_cast<uint8_t>(PaletteMode::TRIADIC);
  uint8_t mapping_mode = static_cast<uint8_t>(Color::PaletteMapping::LINEAR);
  uint8_t hue_mode = static_cast<uint8_t>(HueShiftMode::NOISE);
  uint8_t envelope_mode = static_cast<uint8_t>(EnvelopeMode::NONE);

  static constexpr auto FIELDS = concat_fields<GeneratedPaletteParams>(
      Color::ColorParams::FIELDS,
      std::array<Field<GeneratedPaletteParams>, 0>{});
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<GeneratedPaletteParams>{
          "palette-mode", &GeneratedPaletteParams::palette_mode,
          PALETTE_MODE_IDS, 3, static_cast<uint8_t>(PaletteMode::TRIADIC)},
      TopologyField<GeneratedPaletteParams>{
          "palette-mapping", &GeneratedPaletteParams::mapping_mode,
          PALETTE_MAPPING_IDS, 4,
          static_cast<uint8_t>(Color::PaletteMapping::LINEAR)},
      TopologyField<GeneratedPaletteParams>{
          "hue-shift-mode", &GeneratedPaletteParams::hue_mode,
          HUE_SHIFT_MODE_IDS, 3, static_cast<uint8_t>(HueShiftMode::NOISE)},
      TopologyField<GeneratedPaletteParams>{
          "brightness-envelope", &GeneratedPaletteParams::envelope_mode,
          BRIGHTNESS_ENVELOPE_IDS, 5, static_cast<uint8_t>(EnvelopeMode::NONE)},
  };
};
static_assert(field_ids_unique<GeneratedPaletteParams>());

/** @brief Per-frame color phase clocks. */
struct ColorClockState {
  float oscillation_phase = 0.0f;
  float hue_noise_phase = 0.0f; /**< Bake input for the engine's hue-noise
                                     LUT; not read by run(). */
};

/** @brief FIELD→COLOR crossing: the generated-palette colorizer. */
struct ColorizeGeneratedPalette {
  static constexpr const char *ID = "colorize.generated-palette.v2";
  static constexpr const char *NAME = "Generated Palette";
  using Input = FieldSample;
  using Output = Color4;
  using Params = GeneratedPaletteParams;
  using State = ColorClockState;
  using Prepared = Color::GeneratedPaletteState;

  static constexpr bool APPROXIMATE = true;
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS;
  static constexpr auto METRICS = Color::GENERATED_PALETTE_METRICS;

  static void init(State &, InstanceId) {}
  static Status migrate(State &dst, const State &src, InstanceId) {
    dst = src;
    return Status::OK;
  }
  static void advance(State &state, const Params &params) {
    state.oscillation_phase =
        wrap_t(state.oscillation_phase + params.phase_oscillation_speed);
    state.hue_noise_phase =
        wrap_t(state.hue_noise_phase + params.hue_noise_speed);
  }
  static Prepared prepare(const FrameContext &ctx, const Params &params,
                          const State &state) {
    HS_CHECK(params.hue_mode <= static_cast<uint8_t>(HueShiftMode::PATH_LENGTH),
             "colorize.generated-palette: invalid hue shift mode");
    HS_CHECK(params.palette_mode < ctx.palettes.size(),
             "colorize.generated-palette: invalid palette mode");
    HS_CHECK(params.mapping_mode <=
                 static_cast<uint8_t>(Color::PaletteMapping::REVERSE),
             "colorize.generated-palette: invalid palette mapping");
    HS_CHECK(params.envelope_mode <=
                 static_cast<uint8_t>(EnvelopeMode::DESCENDING),
             "colorize.generated-palette: invalid brightness envelope");
    const Color::HueMode mode = static_cast<HueShiftMode>(params.hue_mode);
    const bool rotation_active = mode != Color::HueMode::NONE &&
                                 params.hue_shift_amount != 0.0f &&
                                 ctx.hue_rotation_lut != nullptr;
    const BakedPalette *palette = ctx.palettes[params.palette_mode];
    HS_CHECK(palette != nullptr,
             "colorize.generated-palette: frame context carries no palette");
    return {Color::PaletteMappingWeights::single(
                static_cast<Color::PaletteMapping>(params.mapping_mode)),
            params.mapping_frequency,
            params.mapping_phase,
            params.phase_oscillation_depth,
            state.oscillation_phase,
            palette,
            mode,
            params.hue_shift_amount,
            {ctx.hue_rotation_lut, rotation_active},
            {ctx.hue_noise_lut, mode == Color::HueMode::NOISE &&
                                    rotation_active &&
                                    ctx.hue_noise_lut != nullptr},
            static_cast<Color::BrightnessEnvelope>(params.envelope_mode),
            params.brightness_depth,
            params.opacity_low,
            params.opacity_high};
  }
  static Color4 run(const FieldSample &input, const FrameContext &,
                    const Params &, const Prepared &prepared) {
    return Color::apply_generated_palette(input, prepared);
  }
};

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
