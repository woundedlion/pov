/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "core/math/projection_patterns.h"

#include "core/render/pullback/operators/common.h"
#include "core/render/pullback/source.h"

/**
 * @file sample.h
 * @brief PLANE→FIELD scalar crossing operator models.
 */

namespace Pullback {

namespace Interp {

namespace Op {

/** @brief Source phase clocks.
 * @details Values are equivalent modulo their periods; column-time samples
 * may lie outside the canonical ranges. */
struct SourceClockState {
  float primary = 0.0f;   ///< Primary phase, radians modulo 2pi.
  float secondary = 0.0f; ///< Secondary phase, radians modulo 2pi.
  float angle = 0.0f;     ///< Rotation angle, radians modulo 2pi.
};

/** @brief Prepared frame clocks of SourceClockState. */
struct SourceClockInterval {
  Spacetime::ClockInterval primary;   ///< Prepared primary interval.
  Spacetime::ClockInterval secondary; ///< Prepared secondary interval.
  Spacetime::ClockInterval angle;     ///< Prepared angle interval.
};

/** @brief Prepares the clocks between consecutive endpoints.
 * @param from Previous state.
 * @param to Current state.
 * @return Prepared interval. */
HS_FLASH_INLINE inline SourceClockInterval
temporal_interval(const SourceClockState &from, const SourceClockState &to) {
  return {{from.primary, to.primary, math::TWO_PI_F},
          {from.secondary, to.secondary, math::TWO_PI_F},
          {from.angle, to.angle, math::TWO_PI_F}};
}

/** @brief Applies the prepared clocks to a column state.
 * @tparam Params Operator parameter block.
 * @param out Copy of the current state.
 * @param interval Prepared frame interval.
 * @param fraction Frame fraction in [0, 1]. */
template <typename Params>
HS_FLASH_INLINE void temporal_sample(SourceClockState &out,
                                     const SourceClockInterval &interval,
                                     const Params &, float fraction) {
  out.primary = interval.primary.sample(fraction);
  out.secondary = interval.secondary.sample(fraction);
  out.angle = interval.angle.sample(fraction);
}

/** @brief Sample crossing weighting of the raw field. */
enum class WeightMode : uint8_t {
  NONE = 0,      ///< Raw field unweighted.
  PROJECTION = 1 ///< Raw field scaled by the projection weight.
};
/** @brief Sample crossing coverage policy. */
using ProjectionCoverageMode = Pullback::ProjectionCoverageMode;

/** @brief Weight topology values, in WeightMode order. */
inline constexpr const char *WEIGHT_MODE_IDS[] = {"none", "projection"};
static_assert(std::size(WEIGHT_MODE_IDS) ==
              static_cast<size_t>(WeightMode::PROJECTION) + 1);
/** @brief Coverage topology values, in ProjectionCoverageMode order. */
inline constexpr const char *COVERAGE_MODE_IDS[] = {
    "none", "weight", "weight-squared", "edge-fade", "edge-fade-squared"};
static_assert(std::size(COVERAGE_MODE_IDS) ==
              static_cast<size_t>(ProjectionCoverageMode::EDGE_FADE_SQUARED) +
                  1);

/**
 * @brief Bounds the Sample crossing's topology enum8s.
 * @details The per-pixel switches rely on this check and carry no guard.
 * @param weight_mode Raw `WeightMode` value.
 * @param coverage_mode Raw `ProjectionCoverageMode` value.
 */
inline void check_sample_topology(uint8_t weight_mode, uint8_t coverage_mode) {
  HS_CHECK(weight_mode <= static_cast<uint8_t>(WeightMode::PROJECTION),
           "sample operator: invalid weight mode");
  HS_CHECK(coverage_mode <=
               static_cast<uint8_t>(ProjectionCoverageMode::EDGE_FADE_SQUARED),
           "sample operator: invalid projection coverage mode");
}

/** @brief The Sample crossing's weight switch over the shared policies.
 *  @param weight_mode Checked `WeightMode` value.
 *  @param raw Raw source value.
 *  @param provenance Projection provenance of the plane sample.
 *  @param ctx Frame context.
 *  @return The weighted field value. */
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

/** @brief The Sample crossing's coverage switch over the shared policies.
 *  @param coverage_mode Checked `ProjectionCoverageMode` value.
 *  @param provenance Projection provenance of the plane sample.
 *  @param edge_width Edge-fade width; read only under edge-fade coverage.
 *  @param ctx Frame context.
 *  @return Coverage in [0, 1]. */
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
  case ProjectionCoverageMode::EDGE_FADE_SQUARED: {
    const float fade = ProjectionCoverage::edge_fade(provenance, edge_width);
    return fade * fade;
  }
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
  /** WeightMode topology value. */
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  /** ProjectionCoverageMode topology value. */
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);
};

/**
 * @brief Whether @p Params spells the crossing's three members with the same
 *        defaults as SampleCrossingParams.
 * @details Members that can occupy the crossing's two tail-padding bytes
 * produce different inherited layouts under MSVC and Itanium ABIs. Those
 * families repeat the crossing members and assert this instead. Appended
 * floats require four-byte alignment and cannot occupy that padding.
 * @tparam Params Sample crossing parameter family.
 * @return True when all three defaults match.
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
    "coverage-mode", live_values(ProjectionCoverageMode::EDGE_FADE,
                                 ProjectionCoverageMode::EDGE_FADE_SQUARED)};

/** @brief The tabled field of SampleCrossingParams, retyped to the family.
 *  @tparam Params Sample crossing parameter family.
 *  @return The edge-width field, gated on edge-fade coverage. */
template <typename Params>
constexpr std::array<Field<Params>, 1> sample_crossing_fields() {
  return {edge_width_field<Params>(&Params::edge_width, "Edge Width",
                                   COVERAGE_EDGE_FADE_GATE)};
}

/**
 * @brief The Sample crossing's shared weight and coverage enum8s, followed by
 *        any family-specific topology fields in @p extra.
 * @tparam Params Sample crossing parameter family.
 * @tparam Extra Additional topology field types.
 * @param extra Family-specific topology fields.
 * @return The weight-mode and coverage-mode fields, then @p extra.
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
    weight and coverage enum8s.
    @tparam Params Sample crossing parameter family.
    @param input Plane sample the value was read at.
    @param raw Raw source value.
    @param params Family parameters.
    @param ctx Frame context.
    @return The field carrier. */
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
    Pattern::advance_clocks(params, state.primary, state.secondary,
                            state.angle);
  }

  using Prepared = Pattern::PreparedSource;

  template <typename Params>
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    return Pattern::prepare(state.primary, state.secondary, state.angle);
  }
};

/// `drift` upper bound of the extended-drift twin-wave family.
inline constexpr float EXTENDED_DRIFT_MAX = 2.0f;

/** @brief A source's fields plus the crossing fields, with `drift` widened to
    `EXTENDED_DRIFT_MAX`.
    @tparam Params Sample crossing parameter family.
    @tparam SourceParams Source parameter block @p Params derives from.
    @return The concatenated field table. */
template <typename Params, typename SourceParams>
consteval auto extended_drift_fields() {
  auto fields = concat_fields<Params>(SourceParams::FIELDS,
                                      sample_crossing_fields<Params>());
  for (auto &field : fields)
    if (std::string_view(field.id) == "drift")
      field.max = EXTENDED_DRIFT_MAX;
  return fields;
}

/** @brief Parameter family of sample.grid.v4: the grid source fields plus the
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
struct SampleGridV4 : SourceClockModel {
  static constexpr const char *ID = "sample.grid.v4";
  static constexpr const char *NAME = "Grid";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = GridSampleParams;

  /** @brief Translation symmetry of the grid.
   * @param params Grid controls.
   * @return The square translation lattice. */
  static TranslationLattice translation_lattice(const Params &params) {
    return Source::translation_lattice(params);
  }
  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Pattern::grid(
        Pattern::PlanePosition{projections::stereo_pattern_args(
            input.coords, params.pattern_freq)},
        static_cast<const Source::GridSourceParams &>(params), prepared);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Parameter family of sample.twin-wave.v3. */
struct TwinWaveSampleParams : Source::TwinWaveSourceParams,
                              SampleCrossingParams {
  static constexpr auto FIELDS =
      extended_drift_fields<TwinWaveSampleParams,
                            Source::TwinWaveSourceParams>();
  static constexpr auto TOPOLOGY =
      sample_crossing_topology<TwinWaveSampleParams>();
};
static_assert(field_ids_unique<TwinWaveSampleParams>());
static_assert(field_defaults_in_range<TwinWaveSampleParams>());

/** @brief PLANE→FIELD crossing: the two-wave interference source. */
struct SampleTwinWaveV3 : SourceClockModel {
  static constexpr const char *ID = "sample.twin-wave.v3";
  static constexpr const char *NAME = "Twin Wave Extended Drift";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = TwinWaveSampleParams;

  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Pattern::twin_wave(
        Pattern::PlanePosition{projections::stereo_pattern_args(
            input.coords, params.pattern_freq)},
        prepared);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Parameter family of sample.bands.v1. */
struct BandsSampleParams : Source::BandsSourceParams, SampleCrossingParams {
  static constexpr auto FIELDS = concat_fields<BandsSampleParams>(
      Source::BandsSourceParams::FIELDS,
      sample_crossing_fields<BandsSampleParams>());
  static constexpr auto TOPOLOGY =
      sample_crossing_topology<BandsSampleParams>();
};
static_assert(field_ids_unique<BandsSampleParams>());
static_assert(
    appended_block_size_matches<BandsSampleParams, SampleCrossingParams, 8>(),
    "appended parameter block must have the expected rounded size");
static_assert(field_defaults_in_range<BandsSampleParams>());

/** @brief Plane-to-field crossing: parallel sine bands in the projection plane. */
struct SampleBands : StatelessModel {
  static constexpr const char *ID = "sample.bands.v1";
  static constexpr const char *NAME = "Bands";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = BandsSampleParams;
  using Prepared = Pattern::PlaneBands::Prepared;

  /** @brief Translation symmetry of the bands.
   * @param params Band controls.
   * @return The oriented translation lattice. */
  static TranslationLattice translation_lattice(const Params &params) {
    return Source::translation_lattice(params);
  }

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    return Pattern::PlaneBands::prepare(params, {});
  }
  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Pattern::PlaneBands::sample({input.coords}, prepared);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Parameter family of sample.spiral.v3. */
struct SpiralSampleParams : Source::SpiralSourceParams {
  float edge_width = 0.1f; ///< Fade band in projection-edge units.
  uint8_t weight_mode =
      static_cast<uint8_t>(WeightMode::PROJECTION); ///< Signal weighting.
  uint8_t coverage_mode = static_cast<uint8_t>(
      ProjectionCoverageMode::WEIGHT); ///< Coverage policy.
  static constexpr auto FIELDS = concat_fields<SpiralSampleParams>(
      Source::SpiralSourceParams::FIELDS,
      sample_crossing_fields<SpiralSampleParams>());
  static constexpr auto TOPOLOGY = sample_crossing_topology<SpiralSampleParams>(
      TopologyField<SpiralSampleParams>{
          "coordinate-system", &SpiralSampleParams::coordinate_system,
          Source::SpiralSourceParams::COORDINATE_SYSTEM_IDS, 0});
};
static_assert(field_ids_unique<SpiralSampleParams>());
static_assert(sample_crossing_defaults_match<SpiralSampleParams>());
static_assert(field_defaults_in_range<SpiralSampleParams>());

/** @brief PLANE→FIELD crossing: the stationary spiral source; zero arms draw
    concentric rings. */
struct SampleSpiralV3 : StatelessModel {
  static constexpr const char *ID = "sample.spiral.v3";
  static constexpr const char *NAME = "Spiral";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = SpiralSampleParams;
  using Prepared =
      Pattern::PlaneSpiral::Prepared; ///< Prepared winding geometry.

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    return Pattern::PlaneSpiral::prepare(params, {});
  }
  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const float raw = Pattern::PlaneSpiral::sample({input.coords}, prepared);
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Kind topology values, in Pattern::TessellationKind order. */
inline constexpr const char *TESSELLATION_KIND_IDS[] = {"triangular", "square",
                                                        "hexagonal"};
static_assert(std::size(TESSELLATION_KIND_IDS) ==
              static_cast<size_t>(Pattern::TessellationKind::HEXAGONAL) + 1);
/** @brief Fill topology values, in Pattern::TessellationFill order. */
inline constexpr const char *TESSELLATION_FILL_IDS[] = {"cells", "edges"};
static_assert(std::size(TESSELLATION_FILL_IDS) ==
              static_cast<size_t>(Pattern::TessellationFill::EDGES) + 1);

/** @brief Parameter family of sample.tessellation.v3. */
struct TessellationSampleParams : Source::TessellationSourceParams {
  /** Edge-fade band width; read only under edge-fade coverage. */
  float edge_width = 0.1f;
  /** WeightMode topology value. */
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  /** ProjectionCoverageMode topology value. */
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);
  /** Pattern::TessellationKind topology value. */
  uint8_t kind = static_cast<uint8_t>(Pattern::TessellationKind::SQUARE);
  /** Pattern::TessellationFill topology value. */
  uint8_t fill = static_cast<uint8_t>(Pattern::TessellationFill::CELLS);

  static constexpr auto FIELDS = concat_fields<TessellationSampleParams>(
      Source::TessellationSourceParams::FIELDS,
      sample_crossing_fields<TessellationSampleParams>());
  static constexpr auto TOPOLOGY =
      sample_crossing_topology<TessellationSampleParams>(
          TopologyField<TessellationSampleParams>{
              "kind", &TessellationSampleParams::kind, TESSELLATION_KIND_IDS,
              static_cast<uint8_t>(Pattern::TessellationKind::SQUARE)},
          TopologyField<TessellationSampleParams>{
              "fill", &TessellationSampleParams::fill, TESSELLATION_FILL_IDS,
              static_cast<uint8_t>(Pattern::TessellationFill::CELLS)});
};
static_assert(field_ids_unique<TessellationSampleParams>());
static_assert(field_defaults_in_range<TessellationSampleParams>());
static_assert(sample_crossing_defaults_match<TessellationSampleParams>());

/** @brief PLANE to FIELD crossing: static square, hexagonal or triangular cell
    primitives or edge bands. */
struct SampleTessellation : StatelessModel {
  static constexpr const char *ID = "sample.tessellation.v3";
  static constexpr const char *NAME = "Tessellation";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = TessellationSampleParams;
  struct Prepared {};

  /** @brief Translation symmetry of the selected tessellation.
   * @param params Tessellation kind and scale.
   * @return The translation lattice in plane units. */
  static TranslationLattice translation_lattice(const Params &params) {
    return Source::tessellation_lattice(
        static_cast<Pattern::TessellationKind>(params.kind), params.cell_scale);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    HS_CHECK(params.kind <=
                 static_cast<uint8_t>(Pattern::TessellationKind::HEXAGONAL),
             "sample.tessellation: invalid kind");
    HS_CHECK(params.fill <=
                 static_cast<uint8_t>(Pattern::TessellationFill::EDGES),
             "sample.tessellation: invalid fill");
    return {};
  }
  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &) {
    const float raw = Pattern::tessellation(
        Pattern::PlanePosition{input.coords},
        static_cast<const Source::TessellationSourceParams &>(params),
        static_cast<Pattern::TessellationKind>(params.kind),
        static_cast<Pattern::TessellationFill>(params.fill));
    return finish_sample(input, raw, params, ctx);
  }
};

/** @brief Parameter family of sample.projected-noise.v2. */
struct ProjectedNoiseSampleParams : Source::NoiseSourceParams {
  /** Edge-fade band width; read only under edge-fade coverage. */
  float edge_width = 0.1f;
  /** WeightMode topology value. */
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  /** ProjectionCoverageMode topology value. */
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);
  /** math::NoiseBasis topology value. */
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

/** @brief The noise sources' prepared block: the owned noise field plus this
    frame's loop offset. */
struct PreparedNoiseSource {
  const FastNoiseLite *noise; ///< The instance-owned noise field.
  math::Vector loop_offset;   ///< Lattice offset of this frame's loop phase.
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
    if (!math::noise_plane_in_domain(input.coords))
      return finish_sample(input, 0.0f, params, ctx);
    const float raw = Pattern::noise_contour(
        *prepared.noise, static_cast<math::NoiseBasis>(params.basis),
        math::noise_projected_coordinate(input.coords, params.noise_scale,
                                         prepared.loop_offset),
        params.noise_contrast);
    return finish_sample(input, raw, params, ctx);
  }
};

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
