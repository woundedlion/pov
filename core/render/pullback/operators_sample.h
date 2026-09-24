/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "math/projection_patterns.h"

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

/** @brief Per-frame source phase clocks. */
struct SourceClockState {
  float primary = 0.0f;
  float secondary = 0.0f;
  float angle = 0.0f;
};

enum class WeightMode : uint8_t { NONE = 0, PROJECTION = 1 };
using ProjectionCoverageMode = Pullback::ProjectionCoverageMode;

inline constexpr const char *WEIGHT_MODE_IDS[] = {"none", "projection"};
static_assert(std::size(WEIGHT_MODE_IDS) ==
              static_cast<size_t>(WeightMode::PROJECTION) + 1);
inline constexpr const char *COVERAGE_MODE_IDS[] = {
    "none", "weight", "weight-squared", "edge-fade"};
static_assert(std::size(COVERAGE_MODE_IDS) ==
              static_cast<size_t>(ProjectionCoverageMode::EDGE_FADE) + 1);

/**
 * @brief Bounds the Sample crossing's topology enum8s.
 * @details Called from prepare(), once per frame, so the per-pixel switches
 * below stay total and carry no guard.
 */
inline void check_sample_topology(uint8_t weight_mode, uint8_t coverage_mode) {
  HS_CHECK(weight_mode <= static_cast<uint8_t>(WeightMode::PROJECTION),
           "sample operator: invalid weight mode");
  HS_CHECK(coverage_mode <=
               static_cast<uint8_t>(ProjectionCoverageMode::EDGE_FADE),
           "sample operator: invalid projection coverage mode");
}

/** @brief The Sample crossing's weight switch over the shared policies. */
inline float weighted_field(uint8_t weight_mode, float raw,
                            const ProjectionProvenance &provenance,
                            const FrameContext &ctx) {
  switch (static_cast<WeightMode>(weight_mode)) {
  case WeightMode::NONE:
    return Weight::None::apply(raw, provenance, ctx);
  case WeightMode::PROJECTION:
    break;
  }
  return Weight::Projection::apply(raw, provenance, ctx);
}

/** @brief The Sample crossing's coverage switch over the shared policies. */
inline float projection_coverage(uint8_t coverage_mode,
                                 const ProjectionProvenance &provenance,
                                 float edge_width, const FrameContext &ctx) {
  switch (static_cast<ProjectionCoverageMode>(coverage_mode)) {
  case ProjectionCoverageMode::NONE:
    return ProjectionCoverage::None::apply(provenance, ctx);
  case ProjectionCoverageMode::WEIGHT_SQUARED:
    return ProjectionCoverage::WeightSquared::apply(provenance, ctx);
  case ProjectionCoverageMode::EDGE_FADE:
    return ProjectionCoverage::edge_fade(provenance, edge_width);
  case ProjectionCoverageMode::WEIGHT:
    break;
  }
  return ProjectionCoverage::Weight::apply(provenance, ctx);
}

/** @brief The parameters every Sample crossing family carries: the edge-fade
    width and the weight and coverage enum8s. */
struct SampleCrossingParams {
  /** Edge-fade band width; read only under edge-fade coverage. */
  float edge_width = 0.1f;
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);
};

/**
 * @brief Whether @p Params spells the crossing's three members with the same
 *        defaults as SampleCrossingParams.
 * @details Members that can occupy the crossing's two tail-padding bytes
 * produce different inherited layouts under MSVC and Itanium ABIs. Those
 * families repeat the crossing members and assert this instead. Appended
 * floats require four-byte alignment and cannot occupy that padding.
 */
template <typename Params> consteval bool sample_crossing_defaults_match() {
  constexpr Params PARAMS{};
  constexpr SampleCrossingParams CROSSING{};
  return PARAMS.edge_width == CROSSING.edge_width &&
         PARAMS.weight_mode == CROSSING.weight_mode &&
         PARAMS.coverage_mode == CROSSING.coverage_mode;
}

/** @brief Activation relation of the crossing's edge-fade width. */
inline constexpr TopologyGate COVERAGE_EDGE_FADE_GATE{
    "coverage-mode", live_values(ProjectionCoverageMode::EDGE_FADE)};

/** @brief The tabled field of SampleCrossingParams, retyped to the family. */
template <typename Params>
constexpr std::array<Field<Params>, 1> sample_crossing_fields() {
  return {edge_width_field<Params>(&Params::edge_width, "Edge Width",
                                   COVERAGE_EDGE_FADE_GATE)};
}

/**
 * @brief The Sample crossing's shared weight and coverage enum8s, followed by
 *        any family-specific topology fields in @p extra.
 */
template <typename Params, typename... Extra>
constexpr std::array<TopologyField<Params>, 2 + sizeof...(Extra)>
sample_crossing_topology(const Extra &...extra) {
  return {TopologyField<Params>{"weight-mode", &Params::weight_mode,
                                WEIGHT_MODE_IDS,
                                static_cast<uint8_t>(WeightMode::PROJECTION)},
          TopologyField<Params>{
              "coverage-mode", &Params::coverage_mode, COVERAGE_MODE_IDS,
              static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT)},
          extra...};
}

/** @brief Builds the field carrier from a raw source value under the family's
    weight and coverage enum8s. */
template <typename Params>
__attribute__((always_inline)) inline FieldSample
finish_sample(const PlaneSample &input, float raw, const Params &params,
              const FrameContext &ctx) {
  return Kernel::sample(
      input, weighted_field(params.weight_mode, raw, input.provenance, ctx),
      projection_coverage(params.coverage_mode, input.provenance,
                          params.edge_width, ctx));
}

/** @brief Base of the Sample crossings driven by the shared phase clocks:
    the crossing's topology check plus the clock block they all prepare. */
struct SourceClockModel : ValueStateModel<SourceClockState> {
  template <typename Params>
  static void advance(State &state, const Params &params) {
    Source::advance_clocks(params, state.primary, state.secondary, state.angle);
  }

  using Prepared = Source::PreparedSource;

  template <typename Params>
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    return Source::prepare(state.primary, state.secondary, state.angle);
  }
};

/** @brief Parameter family of sample.grid.v2: the grid source fields plus the
    crossing's union field and topology enum8s. */
struct GridSampleParams : Source::GridSourceParams, SampleCrossingParams {
  static constexpr auto FIELDS = concat_fields<GridSampleParams>(
      Source::GridSourceParams::FIELDS,
      sample_crossing_fields<GridSampleParams>());
  static constexpr auto TOPOLOGY = sample_crossing_topology<GridSampleParams>();
};
static_assert(field_ids_unique<GridSampleParams>());
static_assert(field_defaults_in_range<GridSampleParams>());

/** @brief PLANE→FIELD crossing: the coupled sine grid source with topology
    weight and coverage modes. */
struct SampleGrid : SourceClockModel {
  static constexpr const char *ID = "sample.grid.v2";
  static constexpr const char *NAME = "Grid";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = GridSampleParams;

  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::grid(
        projections::stereo_pattern_args(input.coords, params.pattern_freq),
        static_cast<const Source::GridSourceParams &>(params), prepared);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Parameter family of sample.twin-wave.v2. */
struct TwinWaveSampleParams : Source::TwinWaveSourceParams,
                              SampleCrossingParams {
  static constexpr auto FIELDS = concat_fields<TwinWaveSampleParams>(
      Source::TwinWaveSourceParams::FIELDS,
      sample_crossing_fields<TwinWaveSampleParams>());
  static constexpr auto TOPOLOGY =
      sample_crossing_topology<TwinWaveSampleParams>();
};
static_assert(field_ids_unique<TwinWaveSampleParams>());
static_assert(field_defaults_in_range<TwinWaveSampleParams>());

/** @brief PLANE→FIELD crossing: the two-wave interference source. */
struct SampleTwinWave : SourceClockModel {
  static constexpr const char *ID = "sample.twin-wave.v2";
  static constexpr const char *NAME = "Twin Wave";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = TwinWaveSampleParams;

  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::twin_wave(
        projections::stereo_pattern_args(input.coords, params.pattern_freq),
        prepared);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Parameter family of sample.rings.v2. */
struct RingsSampleParams : SampleCrossingParams {
  float pattern_freq = 1.0f; /**< Plane-coordinate scale before sampling. */
  float speed = 0.0f;        /**< Per-frame advance of the ring phase. */

  static constexpr auto FIELDS = concat_fields<RingsSampleParams>(
      std::array{
          Field<RingsSampleParams>{
              "pattern-freq", &RingsSampleParams::pattern_freq, "Pattern Freq",
              0.1f, 20.0f, FieldCurve::LERP},
          Field<RingsSampleParams>{"speed", &RingsSampleParams::speed, "Speed",
                                   0.0f, 0.5f, FieldCurve::LERP},
      },
      sample_crossing_fields<RingsSampleParams>());
  static constexpr auto TOPOLOGY =
      sample_crossing_topology<RingsSampleParams>();
};
static_assert(sizeof(RingsSampleParams) ==
              sizeof(SampleCrossingParams) + 2 * sizeof(float));
static_assert(field_ids_unique<RingsSampleParams>());
static_assert(sizeof(RingsSampleParams) == ((sizeof(SampleCrossingParams) + 8 +
                                             alignof(RingsSampleParams) - 1) /
                                            alignof(RingsSampleParams)) *
                                               alignof(RingsSampleParams),
              "appended parameter block must have the expected rounded size");
static_assert(field_defaults_in_range<RingsSampleParams>());

/** @brief PLANE→FIELD crossing: the expanding concentric ring source. */
struct SampleRings : SourceClockModel {
  static constexpr const char *ID = "sample.rings.v2";
  static constexpr const char *NAME = "Rings";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = RingsSampleParams;

  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::rings(
        projections::stereo_pattern_args(input.coords, params.pattern_freq),
        prepared);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Parameter family of sample.spherical-rings.v3. */
using SphericalRingsSampleParams = Source::SphericalRingsSourceParams;

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
    state.phase = fmodf(state.phase + params.speed, math::TWO_PI_F);
  }
  static Prepared prepare(const FrameContext &, const Params &,
                          const State &state) {
    const math::Quaternion orientation =
        math::make_rotation(math::X_AXIS, state.walk.spin_phase) *
        state.walk.wander;
    return {math::rotate(math::Y_AXIS, orientation), state.phase};
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
struct SpiralSampleParams : Source::SpiralSourceParams, SampleCrossingParams {
  static constexpr auto FIELDS = concat_fields<SpiralSampleParams>(
      Source::SpiralSourceParams::FIELDS,
      sample_crossing_fields<SpiralSampleParams>());
  static constexpr auto TOPOLOGY =
      sample_crossing_topology<SpiralSampleParams>();
};
static_assert(field_ids_unique<SpiralSampleParams>());
static_assert(field_defaults_in_range<SpiralSampleParams>());

/** @brief PLANE→FIELD crossing: the rotating spiral source. */
struct SampleSpiral : SourceClockModel {
  static constexpr const char *ID = "sample.spiral.v2";
  static constexpr const char *NAME = "Spiral";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = SpiralSampleParams;

  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::spiral(
        projections::stereo_pattern_args(input.coords, params.pattern_freq),
        prepared);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Parameter family of sample.lattice.v2. */
struct LatticeSampleParams : Source::LatticeSourceParams, SampleCrossingParams {
  static constexpr auto FIELDS = concat_fields<LatticeSampleParams>(
      Source::LatticeSourceParams::FIELDS,
      sample_crossing_fields<LatticeSampleParams>());
  static constexpr auto TOPOLOGY =
      sample_crossing_topology<LatticeSampleParams>();
};
static_assert(field_ids_unique<LatticeSampleParams>());
static_assert(field_defaults_in_range<LatticeSampleParams>());

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
struct FractalSampleParams : Source::FractalSourceParams, SampleCrossingParams {
  static constexpr auto FIELDS = concat_fields<FractalSampleParams>(
      Source::FractalSourceParams::FIELDS,
      sample_crossing_fields<FractalSampleParams>());
  static constexpr auto TOPOLOGY =
      sample_crossing_topology<FractalSampleParams>();
};
static_assert(field_ids_unique<FractalSampleParams>());
static_assert(field_defaults_in_range<FractalSampleParams>());

/** @brief PLANE→FIELD crossing: the animated quadratic escape-time fractal. */
struct SampleFractal : SourceClockModel {
  static constexpr const char *ID = "sample.fractal.v2";
  static constexpr const char *NAME = "Escape Fractal";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = FractalSampleParams;

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
static_assert(std::size(TESSELLATION_KIND_IDS) ==
              static_cast<size_t>(Source::TessellationKind::HEXAGONAL) + 1);

/** @brief Parameter family of sample.tessellation.v2. */
struct TessellationSampleParams : Source::TessellationSourceParams {
  float edge_width = 0.1f;
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);
  uint8_t kind = static_cast<uint8_t>(Source::TessellationKind::TRIANGULAR);

  static constexpr auto FIELDS = concat_fields<TessellationSampleParams>(
      Source::TessellationSourceParams::FIELDS,
      sample_crossing_fields<TessellationSampleParams>());
  static constexpr auto TOPOLOGY =
      sample_crossing_topology<TessellationSampleParams>(
          TopologyField<TessellationSampleParams>{
              "kind", &TessellationSampleParams::kind, TESSELLATION_KIND_IDS,
              static_cast<uint8_t>(Source::TessellationKind::TRIANGULAR)});
};
static_assert(field_ids_unique<TessellationSampleParams>());
static_assert(field_defaults_in_range<TessellationSampleParams>());
static_assert(sample_crossing_defaults_match<TessellationSampleParams>());

/** @brief PLANE→FIELD crossing: rotating polygon edge tessellations. */
struct SampleTessellation : SourceClockModel {
  static constexpr const char *ID = "sample.tessellation.v2";
  static constexpr const char *NAME = "Tessellation";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = TessellationSampleParams;

  static Prepared prepare(const FrameContext &ctx, const Params &params,
                          const State &state) {
    const Prepared prepared = SourceClockModel::prepare(ctx, params, state);
    HS_CHECK(params.kind <=
                 static_cast<uint8_t>(Source::TessellationKind::HEXAGONAL),
             "sample.tessellation: invalid kind");
    return prepared;
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
  float edge_width = 0.1f;
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);
  uint8_t basis = static_cast<uint8_t>(math::NoiseBasis::SIMPLEX);

  static constexpr auto FIELDS = concat_fields<ProjectedNoiseSampleParams>(
      Source::NoiseSourceParams::FIELDS,
      sample_crossing_fields<ProjectedNoiseSampleParams>());
  static constexpr auto TOPOLOGY =
      sample_crossing_topology<ProjectedNoiseSampleParams>(
          TopologyField<ProjectedNoiseSampleParams>{
              "basis", &ProjectedNoiseSampleParams::basis, NOISE_BASIS_IDS,
              static_cast<uint8_t>(math::NoiseBasis::SIMPLEX)});
};
static_assert(field_ids_unique<ProjectedNoiseSampleParams>());
static_assert(field_defaults_in_range<ProjectedNoiseSampleParams>());
static_assert(sample_crossing_defaults_match<ProjectedNoiseSampleParams>());

/** @brief Parameter family of sample.spherical-noise.v3.
    @details No basis topology: the plan pins the spherical contour to the
    simplex basis. */
using SphericalNoiseSampleParams = Source::NoiseSourceParams;

/** @brief The noise sources' prepared block: the owned noise field plus this
    frame's loop offset. */
struct PreparedNoiseSource {
  const FastNoiseLite *noise;
  math::Vector loop_offset;
};

/** @brief PLANE→FIELD crossing: the projected-plane noise contour source. */
struct SampleProjectedNoise : PhaseClockModel<NoisePhaseState> {
  static constexpr const char *ID = "sample.projected-noise.v2";
  static constexpr const char *NAME = "Projected Noise";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = ProjectedNoiseSampleParams;
  using Prepared = PreparedNoiseSource;

  static void init(State &state, InstanceId id) { init_noise_phase(state, id); }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    check_noise_basis(params.basis);
    return {&state.noise, math::noise_projected_loop_offset(state.phase)};
  }
  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::noise_contour(
        *prepared.noise, static_cast<math::NoiseBasis>(params.basis),
        math::noise_projected_coordinate(input.coords, params.noise_scale,
                                         prepared.loop_offset),
        params.noise_contrast);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief SPHERE→FIELD crossing: the sphere-space noise contour source. */
struct SampleSphericalNoise : PhaseClockModel<NoisePhaseState> {
  static constexpr const char *ID = "sample.spherical-noise.v3";
  static constexpr const char *NAME = "Spherical Noise";
  using Input = SphereSample;
  using Output = FieldSample;
  using Params = SphericalNoiseSampleParams;
  using Prepared = PreparedNoiseSource;

  static void init(State &state, InstanceId id) { init_noise_phase(state, id); }
  static Prepared prepare(const FrameContext &, const Params &,
                          const State &state) {
    return {&state.noise, math::noise_sphere_loop_offset(state.phase)};
  }
  static FieldSample run(const SphereSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    const float raw = Source::noise_contour(
        *prepared.noise, math::NoiseBasis::SIMPLEX,
        math::noise_sphere_coordinate(input.dir, params.noise_scale,
                                      prepared.loop_offset),
        params.noise_contrast);
    return Kernel::sample(input, raw);
  }
};

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
