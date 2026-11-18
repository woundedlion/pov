/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#include <cstdio>

#include "hardware/pov_sync_protocol.h"
#include "targets/effects.h"

int main() {
  constexpr unsigned int SHOW_SECONDS[] = {
#define HS_PHANTASM_ROW(name, duration_seconds) duration_seconds,
#include "targets/Phantasm/phantasm_playlist_rows.h"
  };
  constexpr pov::sync::Config SYNC{};
  constexpr unsigned int REVS_PER_SECOND = HS_SHOW_FRAMES_PER_SECOND / 2;
  constexpr unsigned int TRANSITION_REVS =
      SYNC.epoch_repeats + SYNC.commit_revs;

  unsigned long long total_revs = 0;
  for (unsigned long long seconds : SHOW_SECONDS)
    total_revs += seconds * REVS_PER_SECOND + TRANSITION_REVS;

  const auto total_seconds =
      (total_revs + REVS_PER_SECOND / 2) / REVS_PER_SECOND;
  std::printf("%02llu::%02llu::%02llu\n", total_seconds / 3600,
              total_seconds / 60 % 60, total_seconds % 60);
}
