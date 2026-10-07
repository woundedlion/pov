/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// Gates generated pullback manifest topology and metric contracts.
#include "pullback_manifest.generated.h"

#include <algorithm>
#include <cstddef>
#include <string_view>

#include "core/render/pullback/color.h"
#include "core/render/pullback/projection.h"
#include "tests/test_harness.h"

namespace {

constexpr std::string_view domain_name(Pullback::ApproximationDomain domain) {
  switch (domain) {
  case Pullback::ApproximationDomain::PROJECTED_COORDINATE:
    return "PROJECTED_COORDINATE";
  case Pullback::ApproximationDomain::PROJECTED_EDGE_DISTANCE:
    return "PROJECTED_EDGE_DISTANCE";
  case Pullback::ApproximationDomain::COLOR_CHANNEL:
    return "COLOR_CHANNEL";
  case Pullback::ApproximationDomain::FRAMEBUFFER:
    return "FRAMEBUFFER";
  }
  return {};
}

constexpr std::string_view
aggregation_name(Pullback::ApproximationAggregation aggregation) {
  switch (aggregation) {
  case Pullback::ApproximationAggregation::MAXIMUM:
    return "MAXIMUM";
  case Pullback::ApproximationAggregation::MEAN:
    return "MEAN";
  }
  return {};
}

template <size_t N>
void check_oracle_metrics(
    std::string_view oracle_id,
    const std::array<Pullback::ApproximationMetric, N> &core_metrics) {
  size_t count = 0;
  for (const auto &entry : PullbackManifest::ORACLE_METRICS) {
    if (entry.oracle_id != oracle_id)
      continue;
    ++count;
    HS_CONTEXT(entry.domain.data(), static_cast<long long>(count));
    const auto match = std::find_if(
        core_metrics.begin(), core_metrics.end(), [&](const auto &metric) {
          return entry.domain == domain_name(metric.domain) &&
                 entry.aggregation == aggregation_name(metric.aggregation);
        });
    HS_EXPECT_TRUE(match != core_metrics.end());
    if (match != core_metrics.end()) {
      HS_EXPECT_EQ(entry.accepted_limit, match->limit);
      HS_EXPECT_EQ(entry.unit, std::string_view(match->unit));
    }
  }
  HS_EXPECT_EQ(count, N);
}

} // namespace

int main() {
  const hs_test::ModuleScope scope = hs_test::begin_module("pullback_manifest");
  static_assert(PullbackManifest::PRESET_COUNT < 32);
  static_assert(!PullbackManifest::PROGRAMS.empty());
  static_assert(!PullbackManifest::ORACLE_METRICS.empty());
  static_assert(PullbackManifest::BASE_SHA.size() == 40);
  static_assert(PullbackManifest::MANIFEST_SHA256.size() == 64);

  for (size_t i = 0; i < PullbackManifest::PROGRAMS.size(); ++i)
    for (size_t j = i + 1; j < PullbackManifest::PROGRAMS.size(); ++j) {
      HS_CONTEXT(PullbackManifest::PROGRAMS[i].id.data(),
                 static_cast<long long>(j));
      HS_EXPECT_TRUE(PullbackManifest::PROGRAMS[i].topology_key !=
                     PullbackManifest::PROGRAMS[j].topology_key);
    }

  for (size_t i = 0; i < PullbackManifest::ORACLE_METRICS.size(); ++i)
    for (size_t j = i + 1; j < PullbackManifest::ORACLE_METRICS.size(); ++j) {
      const PullbackManifest::OracleMetric &a =
          PullbackManifest::ORACLE_METRICS[i];
      const PullbackManifest::OracleMetric &b =
          PullbackManifest::ORACLE_METRICS[j];
      HS_CONTEXT(a.oracle_id.data(), static_cast<long long>(j));
      HS_EXPECT_TRUE(a.oracle_id != b.oracle_id || a.domain != b.domain ||
                     a.aggregation != b.aggregation);
    }

  check_oracle_metrics("PEIRCE_FAST_SQUARE",
                       Pullback::Projection::PEIRCE_FAST_SQUARE_METRICS);
  check_oracle_metrics("HUE_ROTATION_AND_NOISE_LUTS",
                       Pullback::Color::GENERATED_PALETTE_METRICS);
  HS_EXPECT_EQ(PullbackManifest::ORACLE_METRICS.size(),
               Pullback::Projection::PEIRCE_FAST_SQUARE_METRICS.size() +
                   Pullback::Color::GENERATED_PALETTE_METRICS.size());

  return hs_test::end_module(scope) ? 1 : 0;
}
