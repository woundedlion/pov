/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "pullback_manifest.generated.h"

#include <cstddef>
#include <cstdint>

int main() {
  static_assert(!PullbackManifest::PROGRAMS.empty());
  static_assert(!PullbackManifest::ORACLE_METRICS.empty());
  static_assert(PullbackManifest::BASE_SHA.size() == 40);
  static_assert(PullbackManifest::MANIFEST_SHA256.size() == 64);
  uint16_t preset_mask = 0;
  for (const PullbackManifest::ProgramEntry &program :
       PullbackManifest::PROGRAMS)
    preset_mask |= program.preset_mask;
  if (preset_mask != 0x7fff)
    return 1;
  for (const PullbackManifest::OracleMetric &metric :
       PullbackManifest::ORACLE_METRICS)
    if (metric.measured && metric.measured_baseline > metric.accepted_limit)
      return 1;
  return 0;
}
