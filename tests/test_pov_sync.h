/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Host unit tests for the Phantasm synchronization core (hardware/pov_sync.h)
 * — the spec §12 test plan (docs/specs/phantasm_frame_sync_spec.md).
 */
#pragma once

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

namespace hs_test {
namespace pov_sync_tests {

using namespace pov::sync;

struct EdgeMailboxTestAccess {
  static bool burst_complete(const EdgeMailbox &m, uint32_t now, uint32_t gap) {
    return m.burst_complete(now, gap);
  }
  static BurstSnapshot claim(EdgeMailbox &m) { return m.claim(); }
};

inline bool burst_complete(const EdgeMailbox &m, uint32_t now, uint32_t gap) {
  return EdgeMailboxTestAccess::burst_complete(m, now, gap);
}
inline BurstSnapshot claim(EdgeMailbox &m) {
  return EdgeMailboxTestAccess::claim(m);
}

struct SyncBoardTestAccess {
  static Flywheel &flywheel(SyncBoard &b) { return b.flywheel_mut(); }
  static ContentTracker &content(SyncBoard &b) { return b.content_mut(); }
  static SymbolEmitter &emitter(SyncBoard &b) { return b.emitter; }
  static void maybe_schedule_beacon(SyncBoard &b, uint32_t now) {
    int32_t position = -1;
    b.maybe_schedule_beacon(now, position);
  }
  static const Flywheel &flywheel(const SyncBoard &b) { return b.flywheel(); }
  static const ContentTracker &content(const SyncBoard &b) {
    return b.content();
  }
  static LockState lock(const SyncBoard &b) { return b.lock(); }
  static const Config &config(const SyncBoard &b) { return b.config(); }
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

/**
 * @brief Builds full-rate Phantasm timing (600 MHz, 480 RPM, W=288) with a
 *        shortened content cadence so epoch/beacon scenarios run in
 *        milliseconds of host time.
 * @param effects Number of effects in the test playlist.
 * @return A valid Config with 40-rev effects, beacons every 8 revs, commit
 *         K=2, repeats R=3, grid 4.
 */
inline Config test_config(int effects = 4) {
  Config c = phantasm_config(600000000u, 480u, 288, effects);
  c.revs_per_effect = 40;
  c.beacon_period_revs = 8;
  c.refractory_revs = 8;
  return c;
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

#include "tests/pov_sync/protocol_units.h"
#include "tests/pov_sync/acquire_emit.h"
#include "tests/pov_sync/simulator.h"
#include "tests/pov_sync/failure_budgets.h"

// ── Runner ──────────────────────────────────────────────────────────────────

/**
 * @brief Runs the pov_sync tests.
 * @return The module's failure count.
 */
inline int run_pov_sync_tests() {
  hs_test::ModuleFixture fixture("pov_sync");

  test_helpers();
  test_config_validation();
  test_alphabet();
  test_flip_gate();
  test_mailbox();
  test_mailbox_overlong_burst();
  test_mailbox_prior_staleness();
  test_mailbox_rejects_backward_clock();
  test_seed_clears_mailbox();
  test_build_request_reset();
  test_configure_replaces_claim_windows();
  test_multi_boundary_tick_window();
  test_beacon_codec();
  test_beacon_partial_frame_ages_out();
  test_beacon_shift_needs_confirmation();
  test_beacon_out_of_range_index_rejected();
  test_rev_resync_fold();
  test_flywheel_position();
  test_snap_gate();
  test_suspect_timeout_acquire_uncounted();
  test_isolated_noise_preserves_recent_boundary_lock();
  test_acquire_quiet_before_guard();
  test_acquire_beacon_train_joins();
  test_emitter();
  test_master_beacon_busy_retry();
  test_beacon_late_coast();
  test_master_fold_stall_recovers();
  test_master_fold_stall_recovery_flips();
  test_beacon_tail_quiet();
  test_master_epoch_train_bounded();

  test_sim_boot_and_phase();
  test_sim_eight_board_boot_and_phase();
  test_sim_epoch_commit();
  test_sim_variable_effect_durations();
  test_sim_commit_deadline_trap();
  test_sim_commit_pickup_budget();
  test_sim_masked_windows();
  test_sim_emi();
  test_sim_drops_and_missed_epoch();
  test_sim_reboot(test_config());
  const Config shipping = phantasm_config(600000000u, 480u, 288, 4);
  HS_EXPECT_EQ(shipping.rejoin_bound_revs(), 25u);
  test_sim_reboot(shipping);
  test_sim_forged_burst();
  test_sim_epoch_repeat_lockstep();
  test_sim_rev_resync();
  test_sim_rev_wrap_within_effect();
  test_epoch_same_tick_burst_fold();
  test_epoch_refractory_window();
  test_construction_window_predicates();
  test_effect_output_envelope();
  test_joined_board_dark_through_commit_window();

  test_budget_spurious_epoch();
  test_budget_lost_symbol();
  test_budget_emi_accepted_seam();
  test_budget_corrupted_timebase();
  test_budget_acquire_mis_snap();
  test_budget_beacon_corruption();
  test_beacon_two_edge_substitution();
  test_budget_wire_dead();

  return fixture.result();
}

} // namespace pov_sync_tests
} // namespace hs_test
