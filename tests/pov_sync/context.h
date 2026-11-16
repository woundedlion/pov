/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Host unit tests for the Phantasm synchronization core (hardware/pov_sync.h)
 * — the spec §12 test plan (docs/specs/phantasm_frame_sync_spec.md).
 */
#pragma once

/** @file
 * @brief Shared pov sync test fixtures. */

#include "hardware/pov_handoff.h"
#include "hardware/pov_submit_gate.h"
#include "hardware/pov_sync.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <limits>
#include <vector>

namespace hs_test::pov_sync_tests {

using namespace pov::sync;

struct SyncBoardTestAccess {
  static Flywheel &flywheel(SyncBoard &b) { return b.fly; }
  static ContentTracker &content(SyncBoard &b) { return b.content_tracker; }
  static SymbolEmitter &emitter(SyncBoard &b) { return b.emitter; }
  static void select_boot_effect(SyncBoard &b, int effect_index) {
    b.content_tracker.effect_index = effect_index;
    b.publish_build(effect_index);
  }
  static void maybe_schedule_beacon(SyncBoard &b, uint32_t now) {
    int32_t position = -1;
    b.maybe_schedule_beacon(now, position);
  }
  static const Flywheel &flywheel(const SyncBoard &b) { return b.fly; }
  static const ContentTracker &content(const SyncBoard &b) {
    return b.content_tracker;
  }
  static LockState lock(const SyncBoard &b) { return b.fly.lock(); }
  static const Config &config(const SyncBoard &b) { return b.protocol_config; }
};

inline Flywheel &flywheel_mut(SyncBoard &b) {
  return SyncBoardTestAccess::flywheel(b);
}
inline ContentTracker &content_mut(SyncBoard &b) {
  return SyncBoardTestAccess::content(b);
}
inline const Flywheel &flywheel(const SyncBoard &b) {
  return SyncBoardTestAccess::flywheel(b);
}
inline const ContentTracker &content(const SyncBoard &b) {
  return SyncBoardTestAccess::content(b);
}
inline LockState lock(const SyncBoard &b) {
  return SyncBoardTestAccess::lock(b);
}
inline const Config &config(const SyncBoard &b) {
  return SyncBoardTestAccess::config(b);
}

/** @brief Builds the shipping Phantasm protocol configuration. */
inline Config shipping_config(int effects = 4) {
  return phantasm_config(600000000u, 480u, 288, effects);
}

constexpr uint32_t PERIOD = 37500000u; /**< Cycles per half-rev at full rate. */
constexpr uint32_t COL = PERIOD / 144u; /**< Cycles per column at full rate. */

/**
 * @brief Feeds a five-digit beacon train into a board, spaced exactly as
 *        SymbolEmitter::schedule_beacon emits it.
 * @param board Board under test.
 * @param cfg Protocol configuration.
 * @param col Cycles per column at the run's rate.
 * @param start Cycle at which digit 0 opens.
 * @param d The five encoded digits.
 * @return The cycle just past the last digit's trailing gap.
 */
inline uint32_t feed_beacon_train(SyncBoard &board, const Config &cfg,
                                  uint32_t col, uint32_t start,
                                  const uint8_t d[5]) {
  uint32_t f = start;
  for (int i = 0; i < 5; ++i) {
    const uint32_t span =
        static_cast<uint32_t>(d[i] * cfg.beacon_pitch_cols) * col;
    const BurstSnapshot s{static_cast<uint32_t>(d[i]) + 1u, f, f + span};
    board.tick(f + span + static_cast<uint32_t>(cfg.gap_timeout_cols) * col,
               &s);
    f += span + static_cast<uint32_t>(cfg.gap_timeout_cols + 1) * col;
  }
  return f;
}
/**
 * @brief Builds full-rate Phantasm timing with a shortened content cadence.
 * @param effects Number of effects in the test playlist.
 * @return A valid Config.
 */
inline Config test_config(int effects = 4) {
  Config c = shipping_config(effects);
  c.revs_per_effect = 40;
  c.beacon_period_revs = 8;
  c.refractory_revs = 8;
  return c;
}

} // namespace hs_test::pov_sync_tests
