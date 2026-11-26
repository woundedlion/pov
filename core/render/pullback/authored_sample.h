/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file authored_sample.h
 * @brief Authored planar source crossing and material selectors. */

#include <tuple>

#include "core/render/pullback/authored_source.h"
#include "core/render/pullback/material.h"

namespace Pullback {

/** @brief The parameters every plane-to-field Sample crossing family carries:
    the edge-fade width and the weight and coverage enum8s. */
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

} // namespace Pullback

namespace Pullback::Source {

/// Drift upper bound of the extended-drift twin-wave family, radians per frame.
inline constexpr float EXTENDED_DRIFT_MAX = 2.0f;

/** @brief Twin-wave source controls with the catalog's extended drift range.
 */
struct AuthoredTwinWaveSourceParams : TwinWaveSourceParams {
  /// Numeric controls with drift through two radians per frame.
  static constexpr auto FIELDS = [] {
    auto fields = concat_fields<AuthoredTwinWaveSourceParams>(
        TwinWaveSourceParams::FIELDS,
        std::array<Field<AuthoredTwinWaveSourceParams>, 0>{});
    for (auto &field : fields)
      if (std::string_view(field.id) == "drift")
        field.max = EXTENDED_DRIFT_MAX;
    return fields;
  }();
};

/** @brief A planar source's numeric and runtime crossing controls.
 * @tparam Family Source-specific controls. */
template <typename Family> struct AuthoredSampleParams : Family {
  /// Fade band in projection-edge units.
  float edge_width = 0.1f;
  /// Raw signal weighting catalog index.
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  /// Projection coverage catalog index.
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);
  /// Signal weighting choices.
  static constexpr auto &WEIGHT_IDS = WEIGHT_MODE_IDS;
  /// Projection coverage choices.
  static constexpr auto &COVERAGE_IDS = COVERAGE_MODE_IDS;
  /// Numeric source controls and conditional fade width.
  static constexpr auto FIELDS = concat_fields<AuthoredSampleParams>(
      Family::FIELDS, sample_crossing_fields<AuthoredSampleParams>());
  /// Runtime crossing and source-specific topology controls.
  static constexpr auto TOPOLOGY = [] {
    if constexpr (requires { Family::TOPOLOGY; })
      return std::apply(
          [](const auto &...field) {
            return sample_crossing_topology<AuthoredSampleParams>(
                TopologyField<AuthoredSampleParams>{field}...);
          },
          Family::TOPOLOGY);
    else
      return sample_crossing_topology<AuthoredSampleParams>();
  }();
};
static_assert(
    sample_crossing_defaults_match<AuthoredSampleParams<GridSourceParams>>());

/** @brief Grid source with runtime crossing controls. */
using AuthoredGridSampleParams = AuthoredSampleParams<GridSourceParams>;
/** @brief Twin-wave source with runtime crossing controls. */
using AuthoredTwinWaveSampleParams =
    AuthoredSampleParams<AuthoredTwinWaveSourceParams>;
/** @brief Bands source with runtime crossing controls. */
using AuthoredBandsSampleParams = AuthoredSampleParams<BandsSourceParams>;
/** @brief Spiral source with runtime crossing controls. */
using AuthoredSpiralSampleParams = AuthoredSampleParams<SpiralSourceParams>;
/** @brief Square-cell tessellation with runtime crossing controls. */
using AuthoredSquareCellsSampleParams = AuthoredSampleParams<SquareCellsParams>;
/** @brief Projected noise with runtime crossing controls. */
using AuthoredProjectedNoiseSampleParams =
    AuthoredSampleParams<AuthoredProjectedNoiseSourceParams>;

} // namespace Pullback::Source

namespace Pullback::Weight {

/** @brief Runtime signal weighting over an authored source crossing.
 * @tparam Provider Keyed source control provider. */
template <typename Provider> struct Authored : ApproximationDefaults {
  using Prepared = NoPrepared; ///< Empty validated preparation.

  /** @brief Validates the signal weighting selector.
   * @tparam Frame Provider frame state.
   * @param frame Current frame.
   * @return Empty preparation. */
  template <typename Frame> static NoPrepared prepare(const Frame &frame) {
    const auto &params = Provider::params(frame);
    const uint8_t mode = params.weight_mode;
    HS_CHECK(mode < std::size(params.WEIGHT_IDS),
             "invalid sample signal weighting");
    return {};
  }

  /** @brief Applies the selected raw-signal weight.
   * @tparam Frame Provider frame state.
   * @param field Signed source value.
   * @param provenance Projection metadata.
   * @param frame Current frame.
   * @return Weighted signed source value. */
  template <typename Frame>
  static float apply(float field, const ProjectionProvenance &provenance,
                     const Frame &frame, const NoPrepared &) {
    const auto &params = Provider::params(frame);
    const uint8_t mode = params.weight_mode;
    HS_AUDIT_CHECK(mode < std::size(params.WEIGHT_IDS));
    return select(static_cast<WeightMode>(mode), field, provenance, frame);
  }
};

} // namespace Pullback::Weight

namespace Pullback::ProjectionCoverage {

/** @brief Runtime coverage over an authored source crossing.
 * @tparam Provider Keyed source control provider. */
template <typename Provider> struct Authored : ApproximationDefaults {
  using Prepared = NoPrepared; ///< Empty validated preparation.

  /** @brief Validates the projection coverage selector.
   * @tparam Frame Provider frame state.
   * @param frame Current frame.
   * @return Empty preparation. */
  template <typename Frame> static NoPrepared prepare(const Frame &frame) {
    const auto &params = Provider::params(frame);
    const uint8_t mode = params.coverage_mode;
    HS_CHECK(mode < std::size(params.COVERAGE_IDS),
             "invalid sample projection coverage");
    return {};
  }

  /** @brief Applies the selected projection coverage.
   * @tparam Frame Provider frame state.
   * @param provenance Projection metadata.
   * @param frame Current frame.
   * @return Projection coverage without domain-coverage multiplication. */
  template <typename Frame>
  static float apply(const ProjectionProvenance &provenance, const Frame &frame,
                     const NoPrepared &) {
    const auto &params = Provider::params(frame);
    const uint8_t mode = params.coverage_mode;
    HS_AUDIT_CHECK(mode < std::size(params.COVERAGE_IDS));
    return select(static_cast<ProjectionCoverageMode>(mode), provenance,
                  params.edge_width, frame);
  }
};

} // namespace Pullback::ProjectionCoverage
