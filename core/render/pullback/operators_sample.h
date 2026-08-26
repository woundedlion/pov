/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "render/pullback/operators_common.h"
#include "render/pullback/source.h"

/**
 * @file operators_sample.h
 * @brief PLANE→FIELD crossing operator models: the scalar sources under the
 *        shared weight and coverage topology enum8s.
 */

namespace Pullback {

namespace Interp {

namespace Op {

/** @brief Parameter family of sample.grid.v2: the grid source fields plus the
    crossing's union field and topology enum8s. */
struct GridSampleParams : Source::GridSourceParams {
  /** Edge-fade band width; read only under edge-fade coverage. */
  float edge_width = 0.1f;
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);

  static constexpr auto FIELDS = concat_fields<GridSampleParams>(
      Source::GridSourceParams::FIELDS,
      std::array{edge_width_field(&GridSampleParams::edge_width)});
  static constexpr auto TOPOLOGY = sample_crossing_topology(
      &GridSampleParams::weight_mode, &GridSampleParams::coverage_mode);
};
static_assert(field_ids_unique<GridSampleParams>());

/** @brief PLANE→FIELD crossing: the coupled sine grid source with topology
    weight and coverage modes. */
struct SampleGrid : ValueStateModel<SourceClockState> {
  static constexpr const char *ID = "sample.grid.v2";
  static constexpr const char *NAME = "Grid";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = GridSampleParams;
  using Prepared = Source::PreparedSource;

  static void advance(State &state, const Params &params) {
    state.primary = fmodf(state.primary + params.speed, TWO_PI_F);
    state.secondary =
        fmodf(state.secondary + params.speed * params.secondary_rate, TWO_PI_F);
    state.angle = fmodf(state.angle + params.angle_rate, TWO_PI_F);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    return Source::prepare(state.primary, state.secondary, state.angle);
  }
  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::grid(
        stereo_pattern_args(input.coords, params.pattern_freq),
        static_cast<const Source::GridSourceParams &>(params), prepared);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Parameter family of sample.twin-wave.v2. */
struct TwinWaveSampleParams : Source::TwinWaveSourceParams {
  /** Edge-fade band width; read only under edge-fade coverage. */
  float edge_width = 0.1f;
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);

  static constexpr auto FIELDS = concat_fields<TwinWaveSampleParams>(
      Source::TwinWaveSourceParams::FIELDS,
      std::array{edge_width_field(&TwinWaveSampleParams::edge_width)});
  static constexpr auto TOPOLOGY = sample_crossing_topology(
      &TwinWaveSampleParams::weight_mode, &TwinWaveSampleParams::coverage_mode);
};
static_assert(field_ids_unique<TwinWaveSampleParams>());

/** @brief PLANE→FIELD crossing: the two-wave interference source. */
struct SampleTwinWave : ValueStateModel<SourceClockState> {
  static constexpr const char *ID = "sample.twin-wave.v2";
  static constexpr const char *NAME = "Twin Wave";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = TwinWaveSampleParams;
  using Prepared = Source::PreparedSource;

  static void advance(State &state, const Params &params) {
    state.primary = fmodf(state.primary + params.speed, TWO_PI_F);
    state.secondary =
        fmodf(state.secondary + params.speed * params.secondary_rate, TWO_PI_F);
    state.angle = fmodf(state.angle + params.angle_rate, TWO_PI_F);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    return Source::prepare(state.primary, state.secondary, state.angle);
  }
  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::twin_wave(
        stereo_pattern_args(input.coords, params.pattern_freq), prepared);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Parameter family of sample.rings.v2. */
struct RingsSampleParams {
  float pattern_freq = 1.0f; /**< Plane-coordinate scale before sampling. */
  float speed = 0.0f;        /**< Per-frame advance of the ring phase. */
  /** Edge-fade band width; read only under edge-fade coverage. */
  float edge_width = 0.1f;
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);

  static constexpr auto FIELDS = std::array{
      Field<RingsSampleParams>{"pattern-freq", &RingsSampleParams::pattern_freq,
                               "Pattern Freq", 0.1f, 20.0f, FieldCurve::LERP},
      Field<RingsSampleParams>{"speed", &RingsSampleParams::speed, "Speed",
                               0.0f, 0.5f, FieldCurve::LERP},
      edge_width_field(&RingsSampleParams::edge_width),
  };
  static constexpr auto TOPOLOGY = sample_crossing_topology(
      &RingsSampleParams::weight_mode, &RingsSampleParams::coverage_mode);
};
static_assert(field_ids_unique<RingsSampleParams>());

/** @brief PLANE→FIELD crossing: the expanding concentric ring source. */
struct SampleRings : ValueStateModel<SourceClockState> {
  static constexpr const char *ID = "sample.rings.v2";
  static constexpr const char *NAME = "Rings";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = RingsSampleParams;
  using Prepared = Source::PreparedSource;

  static void advance(State &state, const Params &params) {
    state.primary = fmodf(state.primary + params.speed, TWO_PI_F);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    return Source::prepare(state.primary, state.secondary, state.angle);
  }
  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::rings(
        stereo_pattern_args(input.coords, params.pattern_freq), prepared);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Parameter family of sample.spherical-rings.v3. */
struct SphericalRingsSampleParams : Source::SphericalRingsSourceParams {
  static constexpr auto FIELDS = concat_fields<SphericalRingsSampleParams>(
      Source::SphericalRingsSourceParams::FIELDS,
      std::array<Field<SphericalRingsSampleParams>, 0>{});
};
static_assert(field_ids_unique<SphericalRingsSampleParams>());

/** @brief Instance state of the animated spherical ring source. */
struct SphericalRingsState {
  SpatialWalkState walk;
  float phase = 0.0f;
};

/** @brief SPHERE→FIELD crossing: latitude bands on a wandering, spinning axis. */
struct SampleSphericalRings : ValueStateModel<SphericalRingsState> {
  static constexpr const char *ID = "sample.spherical-rings.v3";
  static constexpr const char *NAME = "Spherical Rings";
  using Input = SphereSample;
  using Output = FieldSample;
  using Params = SphericalRingsSampleParams;
  using Prepared = Source::PreparedSphericalRings;

  static void init(State &state, InstanceId id) {
    init_walk(state.walk, static_cast<int32_t>(id.stable_hash));
  }
  static void advance(State &state, const Params &params) {
    advance_walk(state.walk, params.wander, params.spin_rate);
    state.phase = fmodf(state.phase + params.speed, TWO_PI_F);
  }
  static Prepared prepare(const FrameContext &, const Params &,
                          const State &state) {
    const Quaternion orientation =
        make_rotation(X_AXIS, state.walk.spin_phase) * state.walk.wander;
    return {rotate(Y_AXIS, orientation), state.phase};
  }
  static FieldSample run(const SphereSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::spherical_rings(
        input.dir,
        static_cast<const Source::SphericalRingsSourceParams &>(params),
        prepared);
    return Kernel::sample(input, raw);
  }
};

/** @brief Parameter family of sample.spiral.v2. */
struct SpiralSampleParams : Source::SpiralSourceParams {
  /** Edge-fade band width; read only under edge-fade coverage. */
  float edge_width = 0.1f;
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);

  static constexpr auto FIELDS = concat_fields<SpiralSampleParams>(
      Source::SpiralSourceParams::FIELDS,
      std::array{edge_width_field(&SpiralSampleParams::edge_width)});
  static constexpr auto TOPOLOGY = sample_crossing_topology(
      &SpiralSampleParams::weight_mode, &SpiralSampleParams::coverage_mode);
};
static_assert(field_ids_unique<SpiralSampleParams>());

/** @brief PLANE→FIELD crossing: the rotating spiral source. */
struct SampleSpiral : ValueStateModel<SourceClockState> {
  static constexpr const char *ID = "sample.spiral.v2";
  static constexpr const char *NAME = "Spiral";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = SpiralSampleParams;
  using Prepared = Source::PreparedSource;

  static void advance(State &state, const Params &params) {
    state.primary = fmodf(state.primary + params.speed, TWO_PI_F);
    state.angle = fmodf(state.angle + params.angle_rate, TWO_PI_F);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    return Source::prepare(state.primary, state.secondary, state.angle);
  }
  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::spiral(
        stereo_pattern_args(input.coords, params.pattern_freq), prepared);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Parameter family of sample.lattice.v2. */
struct LatticeSampleParams : Source::LatticeSourceParams {
  /** Edge-fade band width; read only under edge-fade coverage. */
  float edge_width = 0.1f;
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);

  static constexpr auto FIELDS = concat_fields<LatticeSampleParams>(
      Source::LatticeSourceParams::FIELDS,
      std::array{edge_width_field(&LatticeSampleParams::edge_width)});
  static constexpr auto TOPOLOGY = sample_crossing_topology(
      &LatticeSampleParams::weight_mode, &LatticeSampleParams::coverage_mode);
};
static_assert(field_ids_unique<LatticeSampleParams>());

/** @brief PLANE→FIELD crossing: the per-cell primitive lattice source. */
struct SampleLattice : StatelessModel {
  static constexpr const char *ID = "sample.lattice.v2";
  static constexpr const char *NAME = "Lattice";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = LatticeSampleParams;
  struct Prepared {};

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    return {};
  }
  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &) {
    const float raw = Source::primitive_lattice(
        input.coords, static_cast<const Source::LatticeSourceParams &>(params));
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Parameter family of sample.fractal.v2. */
struct FractalSampleParams : Source::FractalSourceParams {
  /** Edge-fade band width; read only under edge-fade coverage. */
  float edge_width = 0.1f;
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);

  static constexpr auto FIELDS = concat_fields<FractalSampleParams>(
      Source::FractalSourceParams::FIELDS,
      std::array{edge_width_field(&FractalSampleParams::edge_width)});
  static constexpr auto TOPOLOGY = sample_crossing_topology(
      &FractalSampleParams::weight_mode, &FractalSampleParams::coverage_mode);
};
static_assert(field_ids_unique<FractalSampleParams>());

/** @brief PLANE→FIELD crossing: the animated quadratic escape-time fractal. */
struct SampleFractal : ValueStateModel<SourceClockState> {
  static constexpr const char *ID = "sample.fractal.v2";
  static constexpr const char *NAME = "Escape Fractal";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = FractalSampleParams;
  using Prepared = Source::PreparedSource;

  static void advance(State &state, const Params &params) {
    state.primary = fmodf(state.primary + params.speed, TWO_PI_F);
    state.angle = fmodf(state.angle + params.angle_rate, TWO_PI_F);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    return Source::prepare(state.primary, state.secondary, state.angle);
  }
  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::escape_fractal(
        input.coords, static_cast<const Source::FractalSourceParams &>(params),
        prepared);
    return finish_sample(input, raw, params, ctx);
  }
};

inline constexpr const char *TESSELLATION_KIND_IDS[] = {"triangular", "square",
                                                        "hexagonal"};

/** @brief Parameter family of sample.tessellation.v2. */
struct TessellationSampleParams : Source::TessellationSourceParams {
  /** Edge-fade band width; read only under edge-fade coverage. */
  float edge_width = 0.1f;
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);
  uint8_t kind = static_cast<uint8_t>(Source::TessellationKind::TRIANGULAR);

  static constexpr auto FIELDS = concat_fields<TessellationSampleParams>(
      Source::TessellationSourceParams::FIELDS,
      std::array{edge_width_field(&TessellationSampleParams::edge_width)});
  static constexpr auto TOPOLOGY = sample_crossing_topology(
      &TessellationSampleParams::weight_mode,
      &TessellationSampleParams::coverage_mode,
      TopologyField<TessellationSampleParams>{
          "kind", &TessellationSampleParams::kind, TESSELLATION_KIND_IDS, 3,
          static_cast<uint8_t>(Source::TessellationKind::TRIANGULAR)});
};
static_assert(field_ids_unique<TessellationSampleParams>());

/** @brief PLANE→FIELD crossing: rotating polygon edge tessellations. */
struct SampleTessellation : ValueStateModel<SourceClockState> {
  static constexpr const char *ID = "sample.tessellation.v2";
  static constexpr const char *NAME = "Tessellation";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = TessellationSampleParams;
  using Prepared = Source::PreparedSource;

  static void advance(State &state, const Params &params) {
    state.angle = fmodf(state.angle + params.angle_rate, TWO_PI_F);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    HS_CHECK(params.kind <=
                 static_cast<uint8_t>(Source::TessellationKind::HEXAGONAL),
             "sample.tessellation: invalid kind");
    return Source::prepare(state.primary, state.secondary, state.angle);
  }
  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::tessellation(
        input.coords,
        static_cast<const Source::TessellationSourceParams &>(params),
        static_cast<Source::TessellationKind>(params.kind), prepared);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Parameter family of sample.projected-noise.v2. */
struct ProjectedNoiseSampleParams : Source::NoiseSourceParams {
  /** Edge-fade band width; read only under edge-fade coverage. */
  float edge_width = 0.1f;
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);
  uint8_t basis = static_cast<uint8_t>(::NoiseBasis::SIMPLEX);

  static constexpr auto FIELDS = concat_fields<ProjectedNoiseSampleParams>(
      Source::NoiseSourceParams::FIELDS,
      std::array{edge_width_field(&ProjectedNoiseSampleParams::edge_width)});
  static constexpr auto TOPOLOGY = sample_crossing_topology(
      &ProjectedNoiseSampleParams::weight_mode,
      &ProjectedNoiseSampleParams::coverage_mode,
      TopologyField<ProjectedNoiseSampleParams>{
          "basis", &ProjectedNoiseSampleParams::basis, NOISE_BASIS_IDS, 3,
          static_cast<uint8_t>(::NoiseBasis::SIMPLEX)});
};
static_assert(field_ids_unique<ProjectedNoiseSampleParams>());

/** @brief Parameter family of sample.spherical-noise.v3.
    @details No basis topology: the plan pins the spherical contour to the
    simplex basis. */
struct SphericalNoiseSampleParams : Source::NoiseSourceParams {
  static constexpr auto FIELDS = concat_fields<SphericalNoiseSampleParams>(
      Source::NoiseSourceParams::FIELDS,
      std::array<Field<SphericalNoiseSampleParams>, 0>{});
};
static_assert(field_ids_unique<SphericalNoiseSampleParams>());

/** @brief The noise sources' prepared block: the owned noise field plus this
    frame's time coordinate. */
struct PreparedNoiseSource {
  const FastNoiseLite *noise;
  float time;
};

/** @brief PLANE→FIELD crossing: the projected-plane noise contour source. */
struct SampleProjectedNoise : ValueStateModel<NoisePhaseState> {
  static constexpr const char *ID = "sample.projected-noise.v2";
  static constexpr const char *NAME = "Projected Noise";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = ProjectedNoiseSampleParams;
  using Prepared = PreparedNoiseSource;

  static void init(State &state, InstanceId id) { init_noise_phase(state, id); }
  static void advance(State &state, const Params &params) {
    state.phase = wrap_t(state.phase + params.noise_time_rate);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    check_noise_basis(params.basis);
    return {&state.noise, state.phase};
  }
  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::noise_contour(
        *prepared.noise, static_cast<::NoiseBasis>(params.basis),
        noise_projected_coordinate(input.coords, params.noise_scale,
                                   prepared.time),
        params.noise_contrast);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief SPHERE→FIELD crossing: the sphere-space noise contour source. */
struct SampleSphericalNoise : ValueStateModel<NoisePhaseState> {
  static constexpr const char *ID = "sample.spherical-noise.v3";
  static constexpr const char *NAME = "Spherical Noise";
  using Input = SphereSample;
  using Output = FieldSample;
  using Params = SphericalNoiseSampleParams;
  using Prepared = PreparedNoiseSource;

  static void init(State &state, InstanceId id) { init_noise_phase(state, id); }
  static void advance(State &state, const Params &params) {
    state.phase = wrap_t(state.phase + params.noise_time_rate);
  }
  static Prepared prepare(const FrameContext &, const Params &,
                          const State &state) {
    return {&state.noise, state.phase};
  }
  static FieldSample run(const SphereSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::noise_contour(
        *prepared.noise, ::NoiseBasis::SIMPLEX,
        noise_sphere_coordinate(input.dir, params.noise_scale, prepared.time),
        params.noise_contrast);
    return Kernel::sample(input, raw);
  }
};

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
