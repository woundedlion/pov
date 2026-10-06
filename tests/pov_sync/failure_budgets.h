/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ── §9.1 failure-mode budget: artifact bounds and recovery times ────────────
//
// Each scenario asserts both the worst-case artifact bound and the recovery
// time of a spec §9.1 budget row.

/**
 * @brief Steps the sim for @p revs, the §9.1 artifact probe.
 * @param sim The simulation to advance.
 * @param revs Number of revolutions to step.
 * @return The worst locked-board phase error (columns vs the master) observed
 *         at any step.
 */
inline double max_err_over(Sim &sim, double revs) {
  const uint64_t until =
      sim.g + static_cast<uint64_t>(revs * 2 * sim.cfg.cycles_per_half_rev);
  double worst = 0.0;
  while (sim.g < until) {
    sim.step();
    worst = std::max(worst, sim.max_phase_err());
  }
  return worst;
}

/**
 * @brief A spurious EPOCH advances one board; beacons repair its index.
 * @details The next real epoch aligns all boards' display frame counters.
 */
inline void test_budget_spurious_epoch() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 0, 0, 0};
  Sim sim(cfg, 4, ppm);
  HS_EXPECT_TRUE(boot_join(sim, cfg));
  const uint64_t REV = 2ull * PERIOD;
  const uint64_t BOUNDARY = (sim.g / REV + 2) * REV;
  sim.emi.push_back({BOUNDARY + COL, 2});
  sim.emi.push_back({BOUNDARY + 3ull * COL, 2});
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.boards[2].live_index == 1; }, 8.0));
  HS_EXPECT_EQ(sim.boards[0].live_index, 0);
  HS_EXPECT_EQ(sim.boards[1].live_index, 0);
  HS_EXPECT_EQ(sim.boards[3].live_index, 0);
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return s.boards[2].live && s.boards[2].live_index == 0; },
      double(2 * cfg.beacon_period_revs + cfg.join_grid_revs)));
  HS_EXPECT_EQ(
      sim.boards[2].board.telemetry_snapshot().beacon_index_corrections, 1u);
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return s.boards[0].live_index == 1 && s.boards[2].live_index == 1;
      },
      double(cfg.revs_per_effect) + 6));
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  for (int i = 0; i < 4; ++i) {
    HS_EXPECT_FALSE(sim.boards[i].trapped);
    HS_EXPECT_EQ(sim.boards[i].live_index, 1);
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
  }
}

/**
 * @brief Verifies the §9.1 "lost boundary symbol" budget row: one dropped
 *        symbol costs a ≤1-revolution coast at crystal drift (~0.01 col at 40
 *        ppm — sub-integer on this probe), re-snapped by the very next symbol.
 * @details The crossing flip covers the missed backstop; nothing is rejected
 *          or misclassified — missed, never wrong.
 */
inline void test_budget_lost_symbol() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 40, -40, 25}; // worst-case datasheet spread
  Sim sim(cfg, 4, ppm);
  sim.run_revs(6.0);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 40; }, 1.1));

  // Deafen board 1 for exactly one half-rev, aligned mid-half: it misses
  // exactly one boundary symbol (the HALF, 104 columns ahead).
  sim.boards[1].drop_from = sim.g;
  sim.boards[1].drop_to = sim.g + PERIOD;
  const uint32_t acc_before =
      sim.boards[1].board.telemetry_snapshot().symbols_accepted;

  HS_EXPECT_LE(max_err_over(sim, 1.5), 1); // sub-column through coast+re-snap
  HS_EXPECT_EQ(lock(sim.boards[1].board), LockState::LOCKED);
  const Telemetry tm = sim.boards[1].board.telemetry_snapshot();
  HS_EXPECT_GE(tm.max_coast_halves, 2u);             // it did coast…
  HS_EXPECT_GE(tm.symbols_accepted, acc_before + 2); // …and re-snapped
  HS_EXPECT_EQ(tm.symbols_rejected_gate, 0u);
  HS_EXPECT_EQ(tm.symbols_discarded_invalid, 0u);
}

/**
 * @brief Verifies the §9.1 "EMI on the sync wire" budget row, the binding
 *        ACCEPTED case: an isolated valid-count burst within G of a predicted
 *        boundary, accepted through the gate.
 * @details Requires the real symbol to be absent at that boundary — with the
 *          real burst present, a nearby forged edge merges inside the gap
 *          timeout into an invalid count and is discarded whole, so the shipped
 *          decoder is stricter than the budget's λ·2G/288 estimate. Artifact: a
 *          ≤G column seam on one board; recovery: the next real symbol, ≤½ rev
 *          later.
 */
inline void test_budget_emi_accepted_seam() {
  const Config cfg = test_config();
  const int32_t ppm[2] = {0, 20};
  Sim sim(cfg, 2, ppm);
  HS_EXPECT_TRUE(boot_join(sim, cfg));
  sim.run_revs(2.0);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 40; }, 1.1));

  // The master HALF boundary is 104 columns ahead. Censor the real symbol
  // for board 1 and forge an edge 3 columns early: isolated, valid count,
  // within G of the predicted boundary — the §9.1 accepted case.
  const uint64_t h = sim.g + 104ull * COL;
  sim.boards[1].drop_from = h - COL;
  sim.boards[1].drop_to = h + 8 * COL;
  sim.emi.push_back({h - 3 * COL, 1});
  sim.emi_pos = 0;
  std::sort(sim.emi.begin(), sim.emi.end());

  // The seam engages (≥2 col — clear of truncation noise, proving the
  // forged snap was really accepted) and is bounded by the gate.
  const double seam = max_err_over(sim, 0.45);
  HS_EXPECT_GE(seam, 2);
  HS_EXPECT_LE(seam, cfg.gate_cols);
  // Recovery: the next real boundary symbol (err ≈ 3 ≤ G) re-snaps.
  sim.run_revs(0.45);
  HS_EXPECT_LE(sim.max_phase_err(), 1);
  HS_EXPECT_EQ(lock(sim.boards[1].board), LockState::LOCKED);
  // Layers 2/3 unharmed: the forged HALF's flip deduped against the
  // crossing, so content stayed equal.
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  HS_EXPECT_EQ(sim.boards[1].t, sim.boards[0].t);
}

/**
 * @brief Verifies the §9.1 "mis-snap despite the gate / corrupted timebase"
 *        budget row: reachable only via a two-coincident-error forge during
 *        ACQUIRE or a firmware bug; the fallback bounds it either way.
 * @details On a corrupted timebase every REAL boundary symbol lands far from a
 *          predicted boundary, so each is first held as a suspect and registers
 *          as a gate rejection only after the 24-column interdigit window (the
 *          §5.3 fallback path): R rejections at ½-rev pace → ACQUIRE → hard
 *          re-snap at the next symbol within about 750 columns (325 ms).
 *          A corruption between symbols adds at most 144 columns of wait.
 */
inline void test_budget_corrupted_timebase() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 20, -20, 10};
  Sim sim(cfg, 4, ppm);
  HS_EXPECT_TRUE(boot_join(sim, cfg));
  // Park at rev 10 (≡ 2 mod the beacon period) so the ACQUIRE window stays
  // clear of beacon trains.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return content(s.boards[0].board).rev_in_effect == 10; },
      16.0));

  // Corrupt board 2's flywheel phase by W/4 — far beyond the gate.
  SimBoard &b2 = sim.boards[2];
  const int32_t bogus = floor_mod(sim.board_pos(2) + 72, cfg.W);
  flywheel_mut(b2.board).seed(Sim::local_now(b2, sim.g) -
                              static_cast<uint32_t>(bogus) * COL);
  flywheel_mut(b2.board).force_lock();
  HS_EXPECT_GE(circ_dist(sim.board_pos(2), sim.board_pos(0), cfg.W), 60);

  const uint64_t CORRUPTED_AT = sim.g;
  constexpr uint64_t RECOVERY_COLUMNS = 750 + 144;
  const uint32_t rej_before =
      b2.board.telemetry_snapshot().symbols_rejected_gate;
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return lock(s.boards[2].board) == LockState::LOCKED &&
               circ_dist(s.board_pos(2), s.board_pos(0), s.cfg.W) <= 1;
      },
      double(RECOVERY_COLUMNS) / cfg.W));
  HS_EXPECT_LE(sim.g - CORRUPTED_AT, RECOVERY_COLUMNS * COL);
  HS_EXPECT_GE(b2.board.telemetry_snapshot().symbols_rejected_gate - rej_before,
               static_cast<uint32_t>(cfg.reject_fallback));
  HS_EXPECT_GE(b2.board.telemetry_snapshot().lock_transitions, 2u);

  // Content recovers fully by the next epoch: any rev_in_effect hiccup from
  // the phase jump is resynced by the beacon rev cross-check (§6.4) well
  // before the train, so the commit is lockstep — same boundary, frame
  // counters re-zeroed together, no trap. (t may carry ±1 from the hiccup
  // only UNTIL that commit.)
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 1)
            return false;
        return true;
      },
      double(cfg.revs_per_effect) + 12));
  for (int i = 1; i < 4; ++i) {
    const int64_t dg = static_cast<int64_t>(sim.boards[i].swap_g) -
                       static_cast<int64_t>(sim.boards[0].swap_g);
    HS_EXPECT_LE(dg < 0 ? -dg : dg, int64_t(3) * COL);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }
  sim.run_revs(2.0);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  for (int i = 1; i < 4; ++i)
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
}

/**
 * @brief Verifies the §9.1 mis-snap row's forged-during-ACQUIRE sub-case: a
 *        board rebooted moments before a scheduled beacon train hard-snaps to
 *        the train's first digit, landing W/4 from the truth, and the
 *        R-rejection fallback still returns it to sub-column phase inside the
 *        row's ~750-column budget.
 * @details The head digit is preceded by wire silence exactly as a boundary
 *          symbol is, so the quiet-before guard cannot filter it and its
 *          1-pulse count reads as a HALF. Every real symbol then lands far from
 *          the broken predictions and is held as a suspect, counted only after
 *          the interdigit window: R rejections at half-rev pace, then ACQUIRE
 *          and a clean re-snap.
 */
inline void test_budget_acquire_mis_snap() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 20, -20, 10};
  Sim sim(cfg, 4, ppm);
  HS_EXPECT_TRUE(boot_join(sim, cfg));
  // Master beacons ride rev ≡ 1 (mod 8); the train starts at x = W/4 = 72.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return content(s.boards[0].board).rev_in_effect == 9; },
      16.0));
  // Reboot past the rev's ZERO burst and a full quiet window ahead of the
  // train, so the first wire event the fresh board meets is the train's head
  // digit and the guard reads the silence before it as isolating.
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 40; }, 1.1));
  SimBoard &b2 = sim.boards[2];
  b2.reboot(Sim::local_now(b2, sim.g));
  HS_EXPECT_EQ(lock(b2.board), LockState::ACQUIRE);

  bool violated = false;
  const auto check_dark = [&](Sim &s) {
    const auto &board = s.boards[2];
    violated |= (board.live && board.live_index != s.boards[0].live_index) ||
                (!board.dark_now &&
                 circ_dist(s.board_pos(2), s.board_pos(0), s.cfg.W) > 1);
  };

  // The head digit captures it: locked on a beacon digit, a quarter turn out.
  HS_EXPECT_TRUE(sim.run_until(
      [&](Sim &s) {
        check_dark(s);
        return lock(s.boards[2].board) == LockState::LOCKED;
      },
      0.5));
  const uint64_t snap_g = sim.g;
  HS_EXPECT_GE(circ_dist(sim.board_pos(2), sim.board_pos(0), cfg.W), 60);

  HS_EXPECT_TRUE(sim.run_until(
      [&](Sim &s) {
        check_dark(s);
        return lock(s.boards[2].board) == LockState::LOCKED &&
               circ_dist(s.board_pos(2), s.board_pos(0), s.cfg.W) <= 1;
      },
      3.0));
  const int32_t recovery_cols = static_cast<int32_t>((sim.g - snap_g) / COL);
  HS_EXPECT_GE(recovery_cols, 500); // the full R-rejection path, not a re-snap
  HS_EXPECT_LE(recovery_cols, 750); // §9.1 mis-snap row
  // Exactly one mis-snap, one fallback, and one clean re-snap since the reboot.
  HS_EXPECT_EQ(b2.board.telemetry_snapshot().lock_transitions, 3u);
  HS_EXPECT_GE(b2.board.telemetry_snapshot().symbols_rejected_gate,
               static_cast<uint32_t>(cfg.reject_fallback));
  HS_EXPECT_FALSE(b2.trapped);
  HS_EXPECT_FALSE(violated);
}

/**
 * @brief Verifies the §9.1 "corrupted beacon frame" budget row: integrity is by
 *        rejection — a corrupted digit fails the checksum, the frame drops
 *        whole with no partial application, and the next beacon (≤ one period
 *        away) cross-checks clean.
 * @details A rejected beacon alone is consequence-free redundancy.
 */
inline void test_budget_beacon_corruption() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 15, -20, 30};
  Sim sim(cfg, 4, ppm);
  HS_EXPECT_TRUE(boot_join(sim, cfg));
  // Master beacons ride rev ≡ 1 (mod 8); park at the rev-9 ZERO crossing.
  // The train starts when the master reaches x = W/4 = 72.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return content(s.boards[0].board).rev_in_effect == 9; },
      16.0));
  // Beacon (index 0, rev 9) is digits [0,0,1,1,0]; its 4th burst is two
  // pulses at relative columns 16 and 17. One EMI edge between them on
  // board 1 (clear of the 100 µs glitch filter) makes that digit read 2:
  // checksum mismatch, whole frame dropped.
  const uint64_t train = sim.g + 72ull * COL;
  sim.emi.push_back({train + 16ull * COL + COL / 2, 1});
  sim.emi_pos = 0;
  std::sort(sim.emi.begin(), sim.emi.end());

  const Telemetry before = sim.boards[1].board.telemetry_snapshot();
  const uint32_t ok_before = before.beacons_ok;
  sim.run_revs(1.0);
  const Telemetry after = sim.boards[1].board.telemetry_snapshot();
  HS_EXPECT_EQ(after.beacons_rejected,
               before.beacons_rejected + 1); // dropped whole
  HS_EXPECT_EQ(after.beacons_ok, ok_before); // nothing applied
  HS_EXPECT_EQ(after.beacon_index_corrections, 0u);
  // Recovery: the next clean beacon decodes within one period.
  HS_EXPECT_TRUE(sim.run_until(
      [ok_before](Sim &s) {
        return s.boards[1].board.telemetry_snapshot().beacons_ok > ok_before;
      },
      double(cfg.beacon_period_revs) + 1));
  HS_EXPECT_EQ(sim.boards[1].live_index, sim.boards[0].live_index);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  HS_EXPECT_LE(sim.max_phase_err(), 2);
  HS_EXPECT_EQ(sim.boards[1].t, sim.boards[0].t);
}

/**
 * @brief Two extra edges in rev 945's low digit must not advance the epoch.
 */
inline void test_beacon_two_edge_substitution() {
  Config cfg = test_config();
  cfg.revs_per_effect = 947;
  const int32_t PPM[3] = {0, 0, 0};
  Sim sim(cfg, 3, PPM);
  HS_EXPECT_TRUE(boot_join(sim, cfg));
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return content(s.boards[0].board).rev_in_effect == 9 &&
               s.board_pos(0) == 40;
      },
      16.0));
  for (auto &board : sim.boards) {
    HS_EXPECT_EQ(content(board.board).rev_in_effect, 9u);
    content_mut(board.board).rev_in_effect = 945;
  }

  // [0,0,6,1,7]: digit 3 pulses at train columns 21 and 22.
  const uint64_t TRAIN = sim.g + 32ull * COL;
  sim.emi.push_back({TRAIN + 21ull * COL + COL / 2, 1});
  sim.emi.push_back({TRAIN + 22ull * COL + COL / 2, 1});
  const Telemetry BEFORE = sim.boards[1].board.telemetry_snapshot();
  const Telemetry CLEAN_BEFORE = sim.boards[2].board.telemetry_snapshot();
  sim.run_revs(0.5);
  const Telemetry AFTER = sim.boards[1].board.telemetry_snapshot();
  const Telemetry CLEAN_AFTER = sim.boards[2].board.telemetry_snapshot();
  HS_EXPECT_EQ(AFTER.beacons_rejected, BEFORE.beacons_rejected + 1);
  HS_EXPECT_EQ(AFTER.beacons_ok, BEFORE.beacons_ok);
  HS_EXPECT_EQ(AFTER.beacon_rev_mismatches, BEFORE.beacon_rev_mismatches);
  HS_EXPECT_EQ(CLEAN_AFTER.beacons_ok, CLEAN_BEFORE.beacons_ok + 1);
  HS_EXPECT_EQ(CLEAN_AFTER.beacons_rejected, CLEAN_BEFORE.beacons_rejected);
  for (const auto &board : sim.boards)
    HS_EXPECT_EQ(content(board.board).rev_in_effect, 945u);

  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return s.boards[0].live_index == 1 && s.boards[1].live_index == 1 &&
               s.boards[2].live_index == 1;
      },
      8.0));
  for (int i = 1; i < 3; ++i) {
    const int64_t DELTA = static_cast<int64_t>(sim.boards[i].swap_g) -
                          static_cast<int64_t>(sim.boards[0].swap_g);
    HS_EXPECT_LE(std::abs(DELTA), static_cast<int64_t>(COL));
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }
}

/**
 * @brief Verifies the §9.1 "sync wire dead" and "master dead" budget rows
 *        (identical for downstream — the master is only the symbol source):
 *        flywheels free-run and keep flipping 2/rev, the playlist freezes on
 *        the current effect, then clears after its fade-out, and boards precess
 *        apart at the §4.5 crystal rate.
 * @details The precession constant τ = T0/δ_rel ≈ one column per 87 revs at 40
 *          ppm — a slow smear, never a break — and the §4.1 rebase rule keeps
 *          the arithmetic valid across a 32-bit cycle-counter wrap with no snaps
 *          at all.
 */
inline void test_budget_wire_dead() {
  const Config cfg = test_config();
  const int32_t ppm[3] = {0, 40, -25};
  // Local clocks start ~30 revs below the 32-bit wrap: the wrap lands ~22
  // revs into the snap-free coast.
  Sim sim(cfg, 3, ppm, 0xFFFFFFFFull - 30ull * 2 * PERIOD + 999);
  HS_EXPECT_TRUE(boot_join(sim, cfg));
  sim.run_revs(2.0);
  // Cut the wire at a quiet point (mid-half, no beacon this rev) so no
  // burst is in flight.
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 40; }, 1.1));

  for (int i = 1; i < 3; ++i) {
    sim.boards[i].drop_from = sim.g;
    sim.boards[i].drop_to = ~0ull;
  }
  uint64_t flips_before[3];
  for (int i = 0; i < 3; ++i)
    flips_before[i] = sim.boards[i].flips;
  const double coast = 150.0;
  sim.run_revs(coast);

  for (int i = 1; i < 3; ++i) {
    // Silence is a coast, not a fault: locked, zero rejections or fallback.
    HS_EXPECT_EQ(lock(sim.boards[i].board), LockState::LOCKED);
    HS_EXPECT_EQ(sim.boards[i].board.telemetry_snapshot().symbols_rejected_gate,
                 0u);
    // Layer 2 self-sufficiency: ~2 flips/rev throughout, no stall.
    const uint64_t df = sim.boards[i].flips - flips_before[i];
    HS_EXPECT_GE(df, 2 * static_cast<uint64_t>(coast) - 5);
    HS_EXPECT_LE(df, 2 * static_cast<uint64_t>(coast) + 5);
    // Layer 3 freezes on the current effect.
    HS_EXPECT_EQ(sim.boards[i].live_index, 0);
    HS_EXPECT_EQ(sim.boards[i].envelope, 0.0f);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }
  // The master alone walks its playlist (3 epochs in 150 revs at the
  // 40+R+K-rev test cadence).
  HS_EXPECT_EQ(sim.boards[0].live_index, 3);
  // Precession matches the budget: 40 ppm × 150 revs × 288 col/rev ≈ 1.7
  // col; 25 ppm ≈ 1.1 col.
  const int32_t e1 = circ_dist(sim.board_pos(1), sim.board_pos(0), cfg.W);
  const int32_t e2 = circ_dist(sim.board_pos(2), sim.board_pos(0), cfg.W);
  HS_EXPECT_GE(e1, 1);
  HS_EXPECT_LE(e1, 3);
  HS_EXPECT_GE(e2, 1);
  HS_EXPECT_LE(e2, 2);
  HS_EXPECT_GE(sim.boards[1].board.telemetry_snapshot().max_coast_halves, 250u);
}

// ── ContentTracker / output-envelope units ────────────────────────────────

/** @brief Pins effect output envelope. */
inline void test_effect_output_envelope() {
  constexpr uint32_t DURATION_REVS = 48;
  constexpr int WIDTH = 288;
  HS_EXPECT_EQ(effect_output_envelope(0, DURATION_REVS, 0, WIDTH), 0.0f);
  HS_EXPECT_NEAR(effect_output_envelope(1, DURATION_REVS, 0, WIDTH), 0.5f,
                 1e-6f);
  HS_EXPECT_EQ(effect_output_envelope(2, DURATION_REVS, 0, WIDTH), 1.0f);
  HS_EXPECT_EQ(effect_output_envelope(46, DURATION_REVS, 0, WIDTH), 1.0f);
  HS_EXPECT_NEAR(effect_output_envelope(47, DURATION_REVS, 0, WIDTH), 0.5f,
                 1e-6f);
  HS_EXPECT_EQ(effect_output_envelope(48, DURATION_REVS, 0, WIDTH), 0.0f);

  float previous_in = 0.0f;
  float previous_out = 1.0f;
  for (int column = 1; column < 2 * WIDTH; ++column) {
    const uint32_t rev = static_cast<uint32_t>(column / WIDTH);
    const int x = column % WIDTH;
    const float fade_in = effect_output_envelope(rev, DURATION_REVS, x, WIDTH);
    const float fade_out =
        effect_output_envelope(46 + rev, DURATION_REVS, x, WIDTH);
    HS_EXPECT_GE(fade_in, 0.0f);
    HS_EXPECT_LE(fade_in, 1.0f);
    HS_EXPECT_GE(fade_out, 0.0f);
    HS_EXPECT_LE(fade_out, 1.0f);
    HS_EXPECT_GE(fade_in, previous_in);
    HS_EXPECT_LE(fade_out, previous_out);
    previous_in = fade_in;
    previous_out = fade_out;
  }
}

/**
 * @brief A board that beacon-joined mid-effect goes dark for the whole commit
 *        window, exactly like the boards that rode the effect from its start.
 * @details The beacon's rev field is six bits, so a board joining an effect
 *          longer than 64 revolutions adopts a count congruent to the master's
 *          but 64k short of it. commit_pending, set by the EPOCH the board
 *          hears, holds it dark.
 */
inline void test_joined_board_dark_through_commit_window() {
  Config cfg = test_config();
  cfg.revs_per_effect = 128; // > 64: the beacon rev field is congruent only
  const uint32_t RPE = cfg.revs_per_effect;
  const uint32_t R = static_cast<uint32_t>(cfg.epoch_repeats);
  const uint32_t K = cfg.commit_revs;
  constexpr int WIDTH = 288;
  constexpr uint32_t JOIN_REV = 70;

  ContentTracker synced;
  synced.identity_known = true;
  synced.rev_in_effect = JOIN_REV;
  ContentTracker joined;
  joined.identity_known = true;
  joined.rev_in_effect = JOIN_REV & 63u; // what a beacon can carry
  HS_EXPECT_EQ(joined.rev_in_effect, 6u);

  for (uint32_t rev = JOIN_REV; rev < RPE; ++rev) {
    HS_EXPECT_FALSE(synced.on_zero_crossing(cfg));
    HS_EXPECT_FALSE(joined.on_zero_crossing(cfg));
  }
  HS_EXPECT_EQ(synced.rev_in_effect, RPE);
  HS_EXPECT_EQ(joined.rev_in_effect, RPE - 64u);

  // The counter-driven envelope is what the two disagree on: this is the state
  // the window gate has to override.
  HS_EXPECT_EQ(effect_output_envelope(synced.rev_in_effect, RPE, 0, WIDTH),
               0.0f);
  HS_EXPECT_EQ(effect_output_envelope(joined.rev_in_effect, RPE, 0, WIDTH),
               1.0f);

  HS_EXPECT_TRUE(synced.on_epoch_symbol(cfg));
  HS_EXPECT_TRUE(joined.on_epoch_symbol(cfg));

  // Announce phase included: the outgoing effect is still live through it, so
  // the envelope is the only thing holding the strip dark.
  bool synced_committed = false;
  bool joined_committed = false;
  for (uint32_t step = 0; step < R + K; ++step) {
    HS_CONTEXT("window rev", static_cast<long long>(step));
    for (int x = 0; x < WIDTH; x += 37) {
      HS_CONTEXT("column", x);
      HS_EXPECT_EQ(synced.output_envelope(cfg, x, WIDTH), 0.0f);
      HS_EXPECT_EQ(joined.output_envelope(cfg, x, WIDTH), 0.0f);
    }
    synced_committed = synced.on_zero_crossing(cfg);
    joined_committed = joined.on_zero_crossing(cfg);
  }
  HS_EXPECT_TRUE(synced_committed);
  HS_EXPECT_TRUE(joined_committed);

  // The commit makes both counts absolute, so the incoming fade-in matches.
  HS_EXPECT_EQ(synced.rev_in_effect, 0u);
  HS_EXPECT_EQ(joined.rev_in_effect, 0u);
  for (int x = 0; x < WIDTH; x += 37) {
    HS_CONTEXT("column", x);
    HS_EXPECT_EQ(joined.output_envelope(cfg, x, WIDTH),
                 synced.output_envelope(cfg, x, WIDTH));
  }
}
