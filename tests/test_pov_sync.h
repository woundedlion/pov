/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Pov sync test-module entry. */

#include "tests/pov_sync/protocol_units.h"
#include "tests/pov_sync/acquire_emit.h"
#include "tests/pov_sync/simulator.h"
#include "tests/pov_sync/failure_budgets.h"
#include "tests/pov_sync/base_cases.h"

namespace hs_test::pov_sync_tests {

/**
 * @brief Runs the pov_sync tests.
 * @return The module's failure count.
 */
inline int run_pov_sync_tests() {
  hs_test::ModuleFixture fixture("pov_sync");

  HS_RUN(test_helpers);
  HS_RUN(test_config_timing_bounds);
  HS_RUN(test_config_validation);
  HS_RUN(test_alphabet);
  HS_RUN(test_flip_gate);
  HS_RUN(test_mailbox);
  HS_RUN(test_mailbox_overlong_burst);
  HS_RUN(test_mailbox_prior_staleness);
  HS_RUN(test_mailbox_rejects_backward_clock);
  HS_RUN(test_seed_clears_mailbox);
  HS_RUN(test_build_request_reset);
  HS_RUN(test_configure_replaces_claim_windows);
  HS_RUN(test_multi_boundary_tick_window);
  HS_RUN(test_beacon_codec);
  HS_RUN(test_beacon_partial_frame_ages_out);
  HS_RUN(test_beacon_shift_needs_confirmation);
  HS_RUN(test_beacon_out_of_range_index_rejected);
  HS_RUN(test_rev_resync_fold);
  HS_RUN(test_flywheel_position);
  HS_RUN(test_snap_gate);
  HS_RUN(test_suspect_timeout_acquire_uncounted);
  HS_RUN(test_isolated_noise_preserves_recent_boundary_lock);
  HS_RUN(test_acquire_isolation_guard);
  HS_RUN(test_acquire_beacon_train_joins);
  HS_RUN(test_emitter_on_time_pitch);
  HS_RUN(test_emitter_boundary_censoring);
  HS_RUN(test_emitter_mid_burst_abort);
  HS_RUN(test_emitter_drop_reporting);
  HS_RUN(test_emitter_beacon_loopback);
  HS_RUN(test_master_beacon_busy_retry);
  HS_RUN(test_beacon_late_coast);
  HS_RUN(test_master_fold_stall_recovers);
  HS_RUN(test_master_fold_stall_recovery_flips);
  HS_RUN(test_beacon_tail_quiet);
  HS_RUN(test_master_epoch_train_bounded);

  HS_RUN(test_sim_boot_and_phase);
  HS_RUN(test_sim_eight_board_boot_and_phase);
  HS_RUN(test_sim_epoch_commit);
  HS_RUN(test_sim_variable_effect_durations);
  HS_RUN(test_sim_commit_beacon_gap);
  HS_RUN(test_sim_rejoin_worst_case);
  HS_RUN(test_sim_commit_deadline_trap);
  HS_RUN(test_sim_commit_pickup_budget);
  HS_RUN(test_sim_masked_windows);
  HS_RUN(test_sim_emi);
  HS_RUN(test_sim_drops_and_missed_epoch);
  HS_RUN(test_sim_reboot, test_config());
  HS_RUN(test_shipping_rejoin_bound);
  HS_RUN(test_sim_reboot, shipping_config());
  HS_RUN(test_sim_forged_burst);
  HS_RUN(test_sim_epoch_repeat_lockstep);
  HS_RUN(test_sim_rev_resync);
  HS_RUN(test_sim_rev_wrap_within_effect);
  HS_RUN(test_epoch_same_tick_burst_fold);
  HS_RUN(test_epoch_refractory_window);
  HS_RUN(test_construction_window_predicates);
  HS_RUN(test_effect_output_envelope);
  HS_RUN(test_joined_board_dark_through_commit_window);

  HS_RUN(test_budget_spurious_epoch);
  HS_RUN(test_budget_lost_symbol);
  HS_RUN(test_budget_emi_accepted_seam);
  HS_RUN(test_budget_corrupted_timebase);
  HS_RUN(test_acquire_ignores_beacon_head);
  HS_RUN(test_budget_beacon_corruption);
  HS_RUN(test_beacon_two_edge_substitution);
  HS_RUN(test_budget_wire_dead);

  return fixture.result();
}
} // namespace hs_test::pov_sync_tests
