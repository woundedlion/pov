/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "render/pullback/material.h"
#include "render/pullback/operator_model.h"
#include "render/pullback/projection.h"
#include "render/pullback/source.h"
#include "render/pullback/stage.h"

/**
 * @file operators.h
 * @brief Concrete chain-interpreter operator models: erased adapters over the
 *        shared carrier kernels, with param-block-backed policy math.
 */

namespace Pullback {

namespace Interp {

namespace Op {

/** @brief Spatial frequency of the operator-owned orientation walk. */
inline constexpr float WALK_NOISE_SCALE = 0.02f;
/** @brief Arc length of one full-amplitude walk step, in radians. */
inline constexpr float WALK_STEP_RADIANS = 0.02f;

/** @brief Deterministic per-frame walk delta from an operator-owned noise
    field; sampled per axis so the pivot direction drifts smoothly. */
inline Quaternion walk_delta(const FastNoiseLite &noise, float t) {
  const Vector axis(noise.GetNoiseSingle(t, 0.0f, 0.0f),
                    noise.GetNoiseSingle(0.0f, t, 31.7f),
                    noise.GetNoiseSingle(-17.3f, t, 0.0f));
  const float length = axis.length();
  if (length < 1e-5f)
    return Quaternion();
  return make_rotation(axis * (1.0f / length), length * WALK_STEP_RADIANS);
}

/**
 * @brief Instance state of the spin-and-wander operators: an accumulated
 *        orientation plus the noise field that drives its random walk.
 */
struct SpatialWalkState {
  FastNoiseLite walk_noise;
  Quaternion wander;
  float spin_phase = 0.0f;
  float walk_time = 0.0f;
};

inline void init_walk(SpatialWalkState &state, InstanceId id) {
  state.walk_noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  state.walk_noise.SetSeed(static_cast<int>(id.stable_hash));
  state.walk_noise.SetFrequency(1.0f);
}

inline void advance_walk(SpatialWalkState &state, float wander,
                         float spin_rate) {
  state.walk_time += 1.0f;
  const Quaternion delta =
      walk_delta(state.walk_noise, state.walk_time * WALK_NOISE_SCALE);
  state.wander =
      (slerp(Quaternion(), delta, wander) * state.wander).normalized();
  state.spin_phase = fmodf(state.spin_phase + spin_rate, TWO_PI_F);
}

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

  static void init(State &state, InstanceId id) { init_walk(state, id); }
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

/** @brief Parameter family of the projection operators. */
struct ProjectChainParams {
  float pole_fade = 1.0f;
  float spin_rate = 0.0f;
  float wander = 0.0f;

  static constexpr auto FIELDS = std::array{
      Field<ProjectChainParams>{"pole-fade", &ProjectChainParams::pole_fade,
                                "Pole Fade", 1.0f, 20.0f, FieldCurve::LERP},
      Field<ProjectChainParams>{
          "projection-spin-speed", &ProjectChainParams::spin_rate,
          "Projection Spin Speed", 0.0f, 0.05f, FieldCurve::LERP},
      Field<ProjectChainParams>{
          "projection-wander", &ProjectChainParams::wander, "Projection Wander",
          0.0f, 1.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<ProjectChainParams>());

/** @brief SPHERE→PLANE crossing: the stereographic projection under a
    wandering, spinning projection frame. */
struct ProjectStereographic {
  static constexpr const char *ID = "project.stereographic.v2";
  static constexpr const char *NAME = "Stereographic";
  using Input = SphereSample;
  using Output = PlaneSample;
  using Params = ProjectChainParams;
  using State = SpatialWalkState;
  struct Prepared {
    Quaternion conjugate;
  };

  static void init(State &state, InstanceId id) { init_walk(state, id); }
  static Status migrate(State &dst, const State &src, InstanceId) {
    dst = src;
    return Status::OK;
  }
  static void advance(State &state, const Params &params) {
    advance_walk(state, params.wander, params.spin_rate);
  }
  static Prepared prepare(const FrameContext &ctx, const Params &,
                          const State &state) {
    return {(make_rotation(Y_AXIS, state.spin_phase) * ctx.projection_base *
             state.wander)
                .conjugate()};
  }
  static PlaneSample run(const SphereSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    const Vector local = rotate(input.dir, prepared.conjugate);
    return Kernel::project(input, local,
                           Projection::stereographic(local, params.pole_fade));
  }
};

enum class WeightMode : uint8_t { NONE = 0, PROJECTION = 1 };
enum class CoverageMode : uint8_t {
  NONE = 0,
  WEIGHT = 1,
  WEIGHT_SQUARED = 2,
  EDGE_FADE = 3
};

inline constexpr const char *WEIGHT_MODE_IDS[] = {"none", "projection"};
inline constexpr const char *COVERAGE_MODE_IDS[] = {
    "none", "weight", "weight-squared", "edge-fade"};

/** @brief Parameter family of sample.grid.v2: the grid source fields plus the
    crossing's union field and topology enum8s. */
struct GridSampleParams : Source::GridSourceParams {
  /** Edge-fade band width; read only under edge-fade coverage. */
  float edge_width = 0.1f;
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  uint8_t coverage_mode = static_cast<uint8_t>(CoverageMode::WEIGHT);

  static constexpr auto FIELDS = concat_fields<GridSampleParams>(
      Source::GridSourceParams::FIELDS,
      std::array{
          Field<GridSampleParams>{"edge-width", &GridSampleParams::edge_width,
                                  "Edge Width", 0.0f, 1.0f, FieldCurve::LERP}});
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<GridSampleParams>{
          "weight-mode", &GridSampleParams::weight_mode, WEIGHT_MODE_IDS, 2,
          static_cast<uint8_t>(WeightMode::PROJECTION)},
      TopologyField<GridSampleParams>{
          "coverage-mode", &GridSampleParams::coverage_mode, COVERAGE_MODE_IDS,
          4, static_cast<uint8_t>(CoverageMode::WEIGHT)},
  };
};
static_assert(field_ids_unique<GridSampleParams>());

/** @brief Per-frame source phase clocks. */
struct SourceClockState {
  float primary = 0.0f;
  float secondary = 0.0f;
  float angle = 0.0f;
};

/** @brief PLANE→FIELD crossing: the coupled sine grid source with topology
    weight and coverage modes. */
struct SampleGrid {
  static constexpr const char *ID = "sample.grid.v2";
  static constexpr const char *NAME = "Grid";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = GridSampleParams;
  using State = SourceClockState;
  using Prepared = Source::PreparedSource;

  static void init(State &, InstanceId) {}
  static Status migrate(State &dst, const State &src, InstanceId) {
    dst = src;
    return Status::OK;
  }
  static void advance(State &state, const Params &params) {
    state.primary = fmodf(state.primary + params.speed, TWO_PI_F);
    state.secondary =
        fmodf(state.secondary + params.speed * params.secondary_rate, TWO_PI_F);
    state.angle = fmodf(state.angle + params.angle_rate, TWO_PI_F);
  }
  static Prepared prepare(const FrameContext &, const Params &,
                          const State &state) {
    return Source::prepare(state.primary, state.secondary, state.angle);
  }
  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::grid(
        stereo_pattern_args(input.coords, params.pattern_freq),
        static_cast<const Source::GridSourceParams &>(params), prepared);
    float weighted;
    switch (static_cast<WeightMode>(params.weight_mode)) {
    case WeightMode::NONE:
      weighted = Weight::None::apply(raw, input.provenance, ctx);
      break;
    case WeightMode::PROJECTION:
    default:
      weighted = Weight::Projection::apply(raw, input.provenance, ctx);
      break;
    }
    float coverage;
    switch (static_cast<CoverageMode>(params.coverage_mode)) {
    case CoverageMode::NONE:
      coverage = ProjectedCoverage::None::apply(input.provenance, ctx);
      break;
    case CoverageMode::WEIGHT_SQUARED:
      coverage = ProjectedCoverage::WeightSquared::apply(input.provenance, ctx);
      break;
    case CoverageMode::EDGE_FADE:
      coverage =
          ProjectedCoverage::edge_fade(input.provenance, params.edge_width);
      break;
    case CoverageMode::WEIGHT:
    default:
      coverage = ProjectedCoverage::Weight::apply(input.provenance, ctx);
      break;
    }
    return Kernel::sample(input, weighted, coverage);
  }
};

enum class PaletteMode : uint8_t {
  TRIADIC = 0,
  COMPLEMENTARY = 1,
  ANALOGOUS = 2
};
enum class HueShiftMode : uint8_t { NOISE = 0, PATH_LENGTH = 1 };
enum class EnvelopeMode : uint8_t { NONE = 0, CUP = 1 };

inline constexpr const char *PALETTE_MODE_IDS[] = {"triadic", "complementary",
                                                   "analogous"};
inline constexpr const char *PALETTE_MAPPING_IDS[] = {"cup", "bell", "linear",
                                                      "reverse"};
inline constexpr const char *HUE_SHIFT_MODE_IDS[] = {"noise", "path-length"};
inline constexpr const char *BRIGHTNESS_ENVELOPE_IDS[] = {"none", "cup"};

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
          HUE_SHIFT_MODE_IDS, 2, static_cast<uint8_t>(HueShiftMode::NOISE)},
      TopologyField<GeneratedPaletteParams>{
          "brightness-envelope", &GeneratedPaletteParams::envelope_mode,
          BRIGHTNESS_ENVELOPE_IDS, 2, static_cast<uint8_t>(EnvelopeMode::NONE)},
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
  static constexpr auto METRICS = Color::GeneratedPalette<void>::METRICS;

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
    const Color::HueMode mode =
        static_cast<HueShiftMode>(params.hue_mode) == HueShiftMode::PATH_LENGTH
            ? Color::HueMode::PATH_LENGTH
            : Color::HueMode::NOISE;
    const bool rotation_active =
        params.hue_shift_amount != 0.0f && ctx.hue_rotation_lut != nullptr;
    const size_t palette_index =
        params.palette_mode < ctx.palettes.size() ? params.palette_mode : 0;
    const BakedPalette *palette = ctx.palettes[palette_index];
    HS_CHECK(palette != nullptr,
             "colorize.generated-palette: frame context carries no palette");
    const uint8_t mapping =
        params.mapping_mode <=
                static_cast<uint8_t>(Color::PaletteMapping::REVERSE)
            ? params.mapping_mode
            : static_cast<uint8_t>(Color::PaletteMapping::LINEAR);
    return {Color::PaletteMappingWeights::single(
                static_cast<Color::PaletteMapping>(mapping)),
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
            static_cast<EnvelopeMode>(params.envelope_mode) == EnvelopeMode::CUP
                ? Color::BrightnessEnvelope::CUP
                : Color::BrightnessEnvelope::NONE,
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
