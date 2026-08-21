/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// Gates the generated pullback manifest: every authored preset is covered by
// some program, and no measured oracle baseline has drifted past the limit its
// manifest accepts.
#include "pullback_manifest.generated.h"

#include <cstddef>
#include <cstdint>
#include <cstdio>

#include "tests/test_harness.h"

/** Coverage bit per authored preset; the union every program must reach. */
static constexpr uint32_t ALL_PRESETS = 0xffffffu;
/** Floor against silent drift: a gutted body would otherwise stay green. */
static constexpr int MIN_ASSERTIONS = 7;

int main() {
  static_assert(!PullbackManifest::PROGRAMS.empty());
  static_assert(!PullbackManifest::ORACLE_METRICS.empty());
  static_assert(PullbackManifest::BASE_SHA.size() == 40);
  static_assert(PullbackManifest::MANIFEST_SHA256.size() == 64);

  uint32_t preset_mask = 0;
  for (const PullbackManifest::ProgramEntry &program :
       PullbackManifest::PROGRAMS)
    preset_mask |= program.preset_mask;
  HS_EXPECT_EQ(preset_mask, ALL_PRESETS);

  for (size_t index = 0; index < PullbackManifest::ORACLE_METRICS.size();
       ++index) {
    const PullbackManifest::OracleMetric &metric =
        PullbackManifest::ORACLE_METRICS[index];
    if (!metric.measured)
      continue;
    HS_CONTEXT(metric.oracle_id.data(), static_cast<long long>(index));
    HS_EXPECT_LE(metric.measured_baseline, metric.accepted_limit);
  }

  const int failed = hs_test::stats().failed;
  const int total = hs_test::stats().passed + failed;
  std::printf("=== pullback_manifest: %d passed, %d failed ===\n",
              hs_test::stats().passed, failed);
  if (total < MIN_ASSERTIONS) {
    std::printf("=== pullback_manifest: only %d assertions ran, expected >= %d "
                "(a check was dropped) ===\n",
                total, MIN_ASSERTIONS);
    return 1;
  }
  return failed ? 1 : 0;
}
