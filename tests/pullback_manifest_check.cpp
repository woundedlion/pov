/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// Gates the generated pullback manifest on the two identity properties
// generate_pullback_manifest_header.py does not check: distinct program
// topology keys, and distinct (oracle, domain, aggregation) metric triples.
// The latter is what test_shader_workbench.h's first-match metric lookup
// relies on; a duplicate triple would silently shadow a baseline.
#include "pullback_manifest.generated.h"

#include <cstddef>
#include <cstdio>

#include "tests/test_harness.h"

int main() {
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

  const int failed = hs_test::stats().failed;
  const int total = hs_test::stats().passed + failed;
  std::printf("=== pullback_manifest: %d passed, %d failed ===\n",
              hs_test::stats().passed, failed);
  if (total == 0) {
    std::printf("=== pullback_manifest: NO ASSERTIONS RAN ===\n");
    return 1;
  }
  return failed ? 1 : 0;
}
