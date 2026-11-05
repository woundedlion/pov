/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file authored_sample.h
 * @brief Authored planar source crossing and material selectors. */

#include <cassert>
#include <tuple>

#include "core/render/pullback/authored_source.h"
#include "core/render/pullback/material.h"

namespace Pullback::Source {

/** @brief Rebinds an enum descriptor to a derived control family.
 * @tparam Derived Target control family.
 * @tparam Base Original control family.
 * @param field Original enum descriptor.
 * @return Descriptor reading the inherited enum member. */
template <typename Derived, typename Base>
constexpr TopologyField<Derived>
rebind_topology(const TopologyField<Base> &field) {
  constexpr const char *PLACEHOLDER[] = {""};
  TopologyField<Derived> result{"", nullptr, PLACEHOLDER, 0};
  result.id = field.id;
  result.member = field.member;
  result.value_ids = field.value_ids;
  result.value_count = field.value_count;
  result.def = field.def;
  return result;
}

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
        field.max = 2.0f;
    return fields;
  }();
};

/** @brief A planar source's numeric and runtime crossing controls.
 * @tparam Family Source-specific controls. */
template <typename Family> struct AuthoredSampleParams : Family {
  float edge_width = 0.1f;   ///< Fade band in projection-edge units.
  uint8_t weight_mode = 1;   ///< Raw signal weighting catalog index.
  uint8_t coverage_mode = 1; ///< Projection coverage catalog index.
  /// Signal weighting choices.
  static constexpr const char *WEIGHT_IDS[] = {"none", "projection"};
  /// Projection coverage choices.
  static constexpr const char *COVERAGE_IDS[] = {
      "none", "weight", "weight-squared", "edge-fade", "edge-fade-squared"};
  /// Numeric source controls and conditional fade width.
  static constexpr auto FIELDS = concat_fields<AuthoredSampleParams>(
      Family::FIELDS,
      std::array{edge_width_field<AuthoredSampleParams>(
          &AuthoredSampleParams::edge_width, "Edge Width",
          {"coverage-mode",
           live_values(ProjectionCoverageMode::EDGE_FADE,
                       ProjectionCoverageMode::EDGE_FADE_SQUARED)})});
  /// Runtime crossing and source-specific topology controls.
  static constexpr auto TOPOLOGY = [] {
    const TopologyField<AuthoredSampleParams> weight{
        "weight-mode", &AuthoredSampleParams::weight_mode, WEIGHT_IDS, 1};
    const TopologyField<AuthoredSampleParams> coverage{
        "coverage-mode", &AuthoredSampleParams::coverage_mode, COVERAGE_IDS, 1};
    if constexpr (requires { Family::TOPOLOGY; })
      return std::apply(
          [&](const auto &...field) {
            return std::array{weight, coverage,
                              rebind_topology<AuthoredSampleParams>(field)...};
          },
          Family::TOPOLOGY);
    else
      return std::array{weight, coverage};
  }();
};

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
    const uint8_t mode = Provider::params(frame).weight_mode;
    HS_CHECK(mode <= 1, "invalid sample signal weighting");
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
    const uint8_t mode = Provider::params(frame).weight_mode;
    assert(mode <= 1);
    return mode == 0 ? field : field * provenance.value_weight;
  }
};

} // namespace Pullback::Weight

namespace Pullback::ProjectionCoverage {

/** @brief Runtime coverage over an authored source crossing.
 * @tparam Provider Keyed source control provider. */
template <typename Provider> struct Authored : ApproximationDefaults {
  /** @brief Applies the selected projection coverage.
   * @tparam Frame Provider frame state.
   * @param provenance Projection metadata.
   * @param frame Current frame.
   * @return Projection coverage without domain-coverage multiplication. */
  template <typename Frame>
  static float apply(const ProjectionProvenance &provenance,
                     const Frame &frame) {
    const auto &params = Provider::params(frame);
    switch (static_cast<ProjectionCoverageMode>(params.coverage_mode)) {
    case ProjectionCoverageMode::NONE:
      return 1.0f;
    case ProjectionCoverageMode::WEIGHT:
      return provenance.value_weight;
    case ProjectionCoverageMode::WEIGHT_SQUARED:
      return provenance.value_weight * provenance.value_weight;
    case ProjectionCoverageMode::EDGE_FADE:
      return edge_fade(provenance, params.edge_width);
    case ProjectionCoverageMode::EDGE_FADE_SQUARED: {
      const float fade = edge_fade(provenance, params.edge_width);
      return fade * fade;
    }
    }
    HS_CHECK(false, "invalid sample projection coverage");
    return 0.0f;
  }
};

} // namespace Pullback::ProjectionCoverage
