/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ── Acceptance gate + acquisition states (§5.3) ─────────────────────────────

/**
 * @brief Verifies the snap acceptance gate and ACQUIRE/LOCKED transitions:
 *        small corrections accepted, an implied W/2 correction rejected, R
 *        consecutive rejections fall back to ACQUIRE (no deadlock), and a hard
 *        snap relocks.
 * @details Exercised for both a forged boundary and a fully corrupted
 *          timebase.
 */
inline void test_snap_gate() {
  const Config cfg = test_config();

  // LOCKED: small correction accepted; the implied W/2 correction of a
  // misclassified boundary (the two-coincident-edge-error residual) is
  // rejected; R consecutive rejections fall back to ACQUIRE (no deadlock).
  {
    Flywheel f(cfg);
    f.seed(1000000u);
    f.force_lock();
    int32_t err = 0;
    // True HALF arrives 2 columns "early" by local reckoning: accept.
    HS_EXPECT_EQ(f.snap(Boundary::HALF, 1000000u + PERIOD - 2 * COL, &err),
                 Flywheel::SnapOutcome::ACCEPTED);
    HS_EXPECT_EQ(err, 2);
    HS_EXPECT_EQ(f.position(1000000u + PERIOD - 2 * COL), 144);

    // Forged ZERO at the HALF position: W/2 correction → reject ×R → ACQUIRE.
    uint32_t t = 1000000u + PERIOD - 2 * COL;
    Flywheel::SnapOutcome last = Flywheel::SnapOutcome::ACCEPTED;
    for (int i = 0; i < cfg.reject_fallback; ++i) {
      t += 10 * COL;
      last = f.snap(Boundary::ZERO, t, &err);
    }
    HS_EXPECT_EQ(last, Flywheel::SnapOutcome::REJECTED_FELL_BACK);
    HS_EXPECT_EQ(f.lock(), LockState::ACQUIRE);
    // ACQUIRE: hard snap, relocks.
    t += 10 * COL;
    HS_EXPECT_EQ(f.snap(Boundary::ZERO, t, &err),
                 Flywheel::SnapOutcome::ACCEPTED);
    HS_EXPECT_EQ(f.lock(), LockState::LOCKED);
    HS_EXPECT_EQ(f.position(t), 0);
  }

  // Corrupted timebase: a board whose epoch is garbage rejects good symbols
  // but re-acquires via the fallback within R symbols (spec §12).
  {
    Flywheel f(cfg);
    // Corrupt: hard-snap to a bogus mid-rev edge (simulates a forged burst
    // accepted during ACQUIRE).
    int32_t err = 0;
    f.seed(1000000u);
    f.snap(Boundary::HALF, 1000000u + 72 * COL, &err); // W/4 off
    f.force_lock();
    // Real boundary stream: ZERO at k·rev, HALF at k·rev + half.
    uint32_t t = 1000000u + PERIOD; // true HALF instant
    Boundary b = Boundary::HALF;
    int accepted_at = -1;
    for (int i = 0; i < cfg.reject_fallback + 1; ++i) {
      if (f.snap(b, t, &err) == Flywheel::SnapOutcome::ACCEPTED) {
        accepted_at = i;
        break;
      }
      t += PERIOD;
      b = opposite(b);
    }
    HS_EXPECT_EQ(accepted_at, cfg.reject_fallback); // R rejections, then snap
    HS_EXPECT_EQ(f.lock(), LockState::LOCKED);
    HS_EXPECT_EQ(f.position(t), boundary_column(b, 288));
  }
}

/**
 * @brief Verifies the suspect-burst timeout counts a gate rejection only while
 *        LOCKED, so symbols_rejected_gate never runs ahead of the §5.3 fallback
 *        it feeds.
 */
inline void test_suspect_timeout_acquire_uncounted() {
  const Config cfg = test_config();
  const uint32_t col = cfg.cycles_per_column();
  SyncBoard board(cfg);
  board.seed(1000u, /*is_master=*/false);
  flywheel_mut(board).force_lock();

  // An isolated valid-count burst at column 40 — far from both boundaries — is
  // held as a suspect awaiting a beacon train.
  const uint32_t head = 1000u + 40u * col;
  const BurstSnapshot suspect{1, head, head};
  board.tick(head + 5 * col, &suspect);
  const uint32_t rejected = board.telemetry_snapshot().symbols_rejected_gate;

  // Fall back to ACQUIRE before the suspect times out.
  for (int i = 0; i < cfg.reject_fallback; ++i)
    flywheel_mut(board).note_rejection();
  HS_EXPECT_EQ(lock(board), LockState::ACQUIRE);
  board.tick(head + 40 * col, nullptr); // past the interdigit window
  HS_EXPECT_EQ(board.telemetry_snapshot().symbols_rejected_gate, rejected);
}

/** @brief Pins isolated noise preserves recent boundary lock. */
inline void test_isolated_noise_preserves_recent_boundary_lock() {
  const Config cfg = test_config();
  const uint32_t col = cfg.cycles_per_column();
  SyncBoard board(cfg);
  board.seed(1000u, false);
  flywheel_mut(board).force_lock();
  for (uint32_t column : {30u, 60u, 90u, 120u}) {
    const uint32_t head = 1000u + column * col;
    const BurstSnapshot noise{1, head, head};
    board.tick(head + 5u * col, &noise);
    board.tick(head + 26u * col, nullptr);
    HS_EXPECT_EQ(lock(board), LockState::LOCKED);
  }
  HS_EXPECT_EQ(board.telemetry_snapshot().symbols_rejected_gate, 4u);
}

/**
 * @brief Verifies the §5.3 quiet-before guard: an ACQUIRE board hard-snaps only
 *        on a burst preceded by t_QB of wire silence, so a beacon digit train
 *        cannot capture a just-rebooted board mid-frame.
 * @details The head of the beacon for effect index 8 is a 2-pulse burst — an
 *          even count, no symbol — and its second digit is a single pulse, a
 *          valid HALF count 5 columns behind it. The same burst after t_QB of
 *          silence snaps. The seed instant opens a quiet window of its own.
 */
inline void test_acquire_quiet_before_guard() {
  const Config cfg = test_config();
  const uint32_t col = cfg.cycles_per_column();
  const uint32_t head = 1000u + 40u * col;
  const BurstSnapshot digit0{2, head, head + col};

  SyncBoard board(cfg);
  board.seed(1000u, /*is_master=*/false);
  HS_EXPECT_EQ(lock(board), LockState::ACQUIRE);
  board.tick(head + col + static_cast<uint32_t>(cfg.gap_timeout_cols) * col,
             &digit0);
  // The next digit, spaced exactly as schedule_beacon spaces them.
  const uint32_t d1 =
      head + col + static_cast<uint32_t>(cfg.gap_timeout_cols + 1) * col;
  const BurstSnapshot digit1{1, d1, d1};
  board.tick(d1 + static_cast<uint32_t>(cfg.gap_timeout_cols) * col, &digit1);
  HS_EXPECT_EQ(lock(board), LockState::ACQUIRE);
  HS_EXPECT_EQ(board.telemetry_snapshot().symbols_accepted, 0u);
  HS_EXPECT_EQ(board.telemetry_snapshot().lock_transitions, 0u);

  SyncBoard quiet(cfg);
  quiet.seed(1000u, /*is_master=*/false);
  quiet.tick(head + col + static_cast<uint32_t>(cfg.gap_timeout_cols) * col,
             &digit0);
  const uint32_t iso = head + col + cfg.acquire_quiet_cycles();
  const BurstSnapshot isolated{1, iso, iso};
  quiet.tick(iso + static_cast<uint32_t>(cfg.gap_timeout_cols) * col,
             &isolated);
  HS_EXPECT_EQ(lock(quiet), LockState::LOCKED);
  HS_EXPECT_EQ(quiet.telemetry_snapshot().symbols_accepted, 1u);
  HS_EXPECT_EQ(flywheel(quiet).position(iso), cfg.W / 2);

  // Boot itself observed no quiet, so the guard measures against the seed
  // instant: a valid-count burst inside that first window is beacon data.
  SyncBoard booted(cfg);
  booted.seed(1000u, /*is_master=*/false);
  const uint32_t early = 1000u + cfg.acquire_quiet_cycles() / 2u;
  const BurstSnapshot interior{1, early, early};
  booted.tick(early + static_cast<uint32_t>(cfg.gap_timeout_cols) * col,
              &interior);
  HS_EXPECT_EQ(lock(booted), LockState::ACQUIRE);
  HS_EXPECT_EQ(booted.telemetry_snapshot().symbols_accepted, 0u);
}

/**
 * @brief Verifies a board still in ACQUIRE decodes a whole beacon train and
 *        adopts (effect, rev) from it (§6.4), instead of waiting for lock.
 * @details A train's first digit reaches both the symbol path and the parser.
 *          The head of the beacon for effect 9 is a 2-pulse burst — an even
 *          count, no symbol — so the board stays in ACQUIRE across the train.
 */
inline void test_acquire_beacon_train_joins() {
  const Config cfg = test_config(16);
  const uint32_t col = cfg.cycles_per_column();
  SyncBoard board(cfg);
  board.seed(1000u, /*is_master=*/false);
  HS_EXPECT_EQ(lock(board), LockState::ACQUIRE);
  HS_EXPECT_FALSE(content(board).identity_known);

  // A full train from column 40 on, spaced exactly as schedule_beacon does.
  uint8_t d[5];
  encode_beacon_digits(9, 5, d);
  HS_EXPECT_EQ(static_cast<int>(d[0]) + 1, 2); // head: no valid symbol count
  uint32_t f = feed_beacon_train(board, cfg, col, 1000u + 40u * col, d);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_ok, 1u);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_rejected, 0u);
  HS_EXPECT_TRUE(content(board).identity_known);
  HS_EXPECT_EQ(content(board).effect_index, 9);
  HS_EXPECT_EQ(content(board).rev_in_effect, 5u);
  // The head reached the symbol path too, and was discarded there on its count.
  HS_EXPECT_EQ(board.telemetry_snapshot().symbols_discarded_invalid, 1u);
  HS_EXPECT_EQ(board.telemetry_snapshot().symbols_accepted, 0u);
  HS_EXPECT_EQ(lock(board), LockState::ACQUIRE);

  // The boundary path is untouched: the next isolated symbol still hard-snaps.
  const uint32_t z = f + cfg.acquire_quiet_cycles();
  const uint32_t zspan = 2u * static_cast<uint32_t>(cfg.pulse_pitch_cols) * col;
  const BurstSnapshot zero{3, z, z + zspan};
  board.tick(z + zspan + static_cast<uint32_t>(cfg.gap_timeout_cols) * col,
             &zero);
  HS_EXPECT_EQ(lock(board), LockState::LOCKED);
  HS_EXPECT_EQ(board.telemetry_snapshot().symbols_accepted, 1u);
  HS_EXPECT_EQ(flywheel(board).position(z), 0);
  HS_EXPECT_EQ(content(board).effect_index, 9);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_ok, 1u);
}

// ── Master emission self-censor (§5.2) ──────────────────────────────────────

/**
 * @brief Verifies master symbol/beacon emission (§5.2): on-time bursts pulse
 *        at 2-column pitch on the oversampled grid, a late boundary is censored
 *        whole, a mid-burst mask past the budget aborts the remaining pulses,
 *        and the emitter→mailbox→parser loop round-trips a beacon frame.
 */
inline void test_emitter() {
  const Config cfg = test_config();
  bool aborted = false;

  // On-time ZERO burst: 3 pulses at exactly 2-column pitch, ticked on an
  // oversampled (⅛-column) grid.
  {
    SymbolEmitter e;
    const uint32_t b = 1000000u;
    HS_EXPECT_TRUE(e.schedule_boundary(Symbol::ZERO, b, b + COL / 8, cfg));
    HS_EXPECT_FALSE(e.schedule_boundary(Symbol::HALF, b, b, cfg));
    std::vector<uint32_t> pulses;
    for (uint32_t t = b + COL / 8; t < b + 8 * COL; t += COL / 8) {
      if (e.tick(t, cfg, &aborted))
        pulses.push_back(t);
      HS_EXPECT_FALSE(aborted);
    }
    HS_EXPECT_EQ(pulses.size(), static_cast<size_t>(3));
    if (pulses.size() == 3) {
      // Each pulse within an oversample step of its scheduled slot.
      HS_EXPECT_LE(pulses[0] - b, COL / 8);
      HS_EXPECT_LE(pulses[1] - (b + 2 * COL), COL / 8);
      HS_EXPECT_LE(pulses[2] - (b + 4 * COL), COL / 8);
    }
  }

  // Late at the boundary (> ~½ column): the whole symbol is censored.
  {
    SymbolEmitter e;
    HS_EXPECT_FALSE(e.schedule_boundary(
        Symbol::ZERO, 1000u, 1000u + cfg.late_censor_cycles() + 1, cfg));
  }

  // Boundary scheduled in the future (now before at_cycles) is not late: it is
  // accepted and emitted once `now` reaches the boundary.
  {
    SymbolEmitter e;
    const uint32_t at = 1000000u;
    const uint32_t early = at - 2u * cfg.late_censor_cycles(); // well before
    HS_EXPECT_TRUE(e.schedule_boundary(Symbol::ZERO, at, early, cfg));
    HS_EXPECT_FALSE(e.tick(early, cfg, &aborted)); // not due yet, no pulse
    HS_EXPECT_FALSE(aborted);
    HS_EXPECT_TRUE(e.tick(at, cfg, &aborted)); // first pulse at the boundary
    HS_EXPECT_FALSE(aborted);
  }

  for (const Symbol symbol : {Symbol::ZERO, Symbol::ZERO_EPOCH}) {
    for (uint32_t sent = 1; sent < symbol_pulse_count(symbol); ++sent) {
      SymbolEmitter e;
      EdgeMailbox mailbox;
      const uint32_t START = 1000000u;
      HS_EXPECT_TRUE(e.schedule_boundary(symbol, START, START, cfg));
      for (uint32_t i = 0; i < sent; ++i) {
        const uint32_t NOW = START + i * cfg.pulse_pitch_cycles();
        HS_EXPECT_TRUE(e.tick(NOW, cfg, &aborted));
        mailbox.on_edge(NOW, cfg.glitch_filter_cycles);
      }
      const uint32_t LATE = START + sent * cfg.pulse_pitch_cycles() +
                            cfg.late_censor_cycles() + 1;
      if (e.tick(LATE, cfg, &aborted))
        mailbox.on_edge(LATE, cfg.glitch_filter_cycles);
      HS_EXPECT_TRUE(aborted);
      HS_EXPECT_TRUE(e.idle());
      const auto burst = claim(mailbox);
      HS_EXPECT_EQ(classify_count(burst.count), Symbol::INVALID);
    }
  }

  {
    SymbolEmitter e;
    const uint32_t START = 1000000u;
    HS_EXPECT_TRUE(e.schedule_boundary(Symbol::ZERO, START, START, cfg));
    HS_EXPECT_TRUE(e.tick(START, cfg, &aborted));
    HS_EXPECT_FALSE(
        e.tick(START + cfg.gap_timeout_cycles() + 1, cfg, &aborted));
    HS_EXPECT_TRUE(aborted);
    HS_EXPECT_TRUE(e.idle());
  }

  // A burst still in flight when a boundary crossing arrives is stale (a
  // masked-ISR coast past HALF); drop_pending_emission clears it so the on-time
  // boundary symbol schedules without tripping the overlap trap, and reports
  // which kind it dropped.
  using Dropped = SymbolEmitter::DroppedBurst;
  {
    SymbolEmitter e;
    HS_EXPECT_EQ(e.drop_pending_emission(), Dropped::NONE); // idle
    uint8_t d[5];
    encode_beacon_digits(3, 9, d);
    HS_EXPECT_TRUE(e.schedule_beacon(d, 2000000u, cfg));
    HS_EXPECT_FALSE(e.schedule_beacon(d, 2000000u, cfg));
    HS_EXPECT_FALSE(e.schedule_boundary(Symbol::HALF, 2000000u, 2000000u, cfg));
    HS_EXPECT_FALSE(e.idle());
    HS_EXPECT_EQ(e.drop_pending_emission(), Dropped::BEACON);
    HS_EXPECT_TRUE(e.idle());
    HS_EXPECT_TRUE(e.schedule_boundary(Symbol::HALF, 2000000u, 2000000u, cfg));

    // An undrained boundary symbol is not a beacon drop.
    HS_EXPECT_EQ(e.drop_pending_emission(), Dropped::BOUNDARY);
    HS_EXPECT_TRUE(e.idle());
  }

  // A drained-but-not-yet-retired beacon frame must not make the boundary symbol
  // that follows it look like a beacon drop.
  {
    SymbolEmitter e;
    uint8_t d[5];
    encode_beacon_digits(3, 9, d);
    e.schedule_beacon(d, 2000000u, cfg);
    bool aborted = false;
    for (uint32_t t = 0; t < 4000u && !e.idle(); ++t)
      e.tick(2000000u + t * (COL / 4), cfg, &aborted);
    HS_EXPECT_TRUE(e.idle());
    HS_EXPECT_TRUE(e.schedule_boundary(Symbol::HALF, 3000000u, 3000000u, cfg));
    HS_EXPECT_EQ(e.drop_pending_emission(), Dropped::BOUNDARY);
  }

  // Beacon: emitter → mailbox → parser closes the loop; inter-burst gaps
  // terminate digits; the decoded frame matches the encoded one.
  {
    SymbolEmitter e;
    EdgeMailbox m;
    BeaconParser p;
    uint8_t d[5];
    encode_beacon_digits(13, 22, d);
    const uint32_t t0 = 5000000u;
    e.schedule_beacon(d, t0, cfg);
    BeaconFrame f{};
    bool got = false, rejected = false;
    for (uint32_t t = t0; t < t0 + 100 * COL; t += COL / 8) {
      if (e.tick(t, cfg, &aborted))
        m.on_edge(t, cfg.glitch_filter_cycles);
      if (burst_complete(m, t, cfg.gap_timeout_cycles())) {
        bool r = false;
        if (p.feed(claim(m), cfg, &f, &r))
          got = true;
        rejected = rejected || r;
      }
    }
    HS_EXPECT_TRUE(got);
    HS_EXPECT_FALSE(rejected);
    HS_EXPECT_EQ(f.effect_index, 13);
    HS_EXPECT_EQ(f.rev_count, 22u);
  }
}

/** @brief Pins master beacon busy retry. */
inline void test_master_beacon_busy_retry() {
  const Config cfg = test_config();
  SyncBoard board(cfg);
  board.seed(0, true);
  content_mut(board).rev_in_effect = 1;

  const uint32_t beacon_at = PERIOD / 2;
  SymbolEmitter &emitter = SyncBoardTestAccess::emitter(board);
  HS_EXPECT_TRUE(
      emitter.schedule_boundary(Symbol::ZERO, beacon_at + COL, beacon_at, cfg));

  SyncBoardTestAccess::maybe_schedule_beacon(board, beacon_at);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_busy_dropped, 1u);

  SyncBoardTestAccess::maybe_schedule_beacon(board, beacon_at);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_busy_dropped, 1u);

  HS_EXPECT_EQ(emitter.drop_pending_emission(),
               SymbolEmitter::DroppedBurst::BOUNDARY);
  SyncBoardTestAccess::maybe_schedule_beacon(board, beacon_at);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_busy_dropped, 1u);
  HS_EXPECT_EQ(emitter.drop_pending_emission(),
               SymbolEmitter::DroppedBurst::BEACON);

  SyncBoardTestAccess::maybe_schedule_beacon(board, beacon_at);
  HS_EXPECT_EQ(emitter.drop_pending_emission(),
               SymbolEmitter::DroppedBurst::NONE);
}

// ── Master beacon late-coast bound (§6.4) ───────────────────────────────────

/**
 * @brief Verifies a masked-ISR coast that reaches the beacon point late does
 *        not queue a beacon whose tail overruns HALF and trips the emitter's
 *        wire-busy trap on the on-time HALF symbol.
 * @details Drives a master across the beacon-due revolution, resuming its first
 *          post-ZERO tick at a chosen column to model the coast. The last
 *          admissible start emits fully before HALF; one column later is
 *          censored, with no pulses in the beacon window and no trap at HALF.
 */
inline void test_beacon_late_coast() {
  const Config cfg = test_config();
  const uint32_t period = cfg.cycles_per_half_rev;

  uint32_t late_dropped = 0;
  // Resume the master's first post-ZERO tick of the beacon-due revolution at
  // `resume_col` plus `sub_col` cycles; return the pulses emitted with column in
  // [W/4, W/2) (purely beacon — the ZERO boundary symbol drained at columns 0..4
  // below W/4).
  auto run = [&](int32_t resume_col, uint32_t sub_col = 0) -> int {
    SyncBoard m(cfg);
    const uint32_t t0 = 1000000u;
    m.seed(t0, /*is_master=*/true);
    // One tick a full revolution ahead folds HALF then ZERO in a single wake,
    // landing boundary ZERO with rev_in_effect == 1 (a beacon-due rev under
    // test_config's epoch repeats). epoch1 is that ZERO instant.
    const uint32_t epoch1 = t0 + 2u * period;
    m.tick(epoch1, nullptr);
    // Drain the rev-1 ZERO symbol's remaining pulses (columns 2, 4) so the
    // emitter is idle before the coast.
    m.tick(epoch1 + 2u * COL, nullptr);
    m.tick(epoch1 + 4u * COL, nullptr);

    int pulses = 0;
    for (int32_t c = resume_col; c <= 150; ++c) {
      // Truncated cycles_per_column() falls before column c; position() floors
      // that instant to c-1, so use the rounded-up rational instant.
      const uint32_t at =
          epoch1 +
          static_cast<uint32_t>((static_cast<uint64_t>(c) * period + 143u) /
                                144u) +
          sub_col;
      const TickActions a = m.tick(at, nullptr);
      if (a.pulse && c >= cfg.W / 4 && c < cfg.W / 2)
        ++pulses;
    }
    late_dropped = m.telemetry_snapshot().beacons_late_dropped;
    return pulses;
  };

  // The rev-1 beacon of effect 0 — the payload the coast above carries.
  uint8_t digits[5];
  encode_beacon_digits(0, 1, digits);
  int32_t digit_sum = 0;
  for (int i = 0; i < 5; ++i)
    digit_sum += digits[i];
  const int32_t last_start = cfg.W / 2 - cfg.beacon_frame_cols(digit_sum) - 1;
  // The payload-sized window spans more than one column.
  HS_EXPECT_GT(last_start, cfg.W / 4);

  // On-time at the beacon point: the frame schedules and emits fully.
  HS_EXPECT_GT(run(cfg.W / 4), 0);
  // The last admissible start still emits.
  HS_EXPECT_GT(run(last_start), 0);
  HS_EXPECT_EQ(late_dropped, 0u);
  // Late start past the safe bound: censored — no beacon pulses, and the HALF
  // crossing at column 144 schedules without tripping the wire-busy trap.
  HS_EXPECT_EQ(run(last_start + 1), 0);
  // The skip is counted once for the revolution, not once per late tick.
  HS_EXPECT_EQ(late_dropped, 1u);
  // Resuming part-way through the last admissible column is late too: the frame
  // is anchored on the tick, and its last pulse may go out up to the emitter's
  // ½-column lateness budget after its due time.
  HS_EXPECT_EQ(run(last_start, COL / 2 + COL / 8), 0);
  HS_EXPECT_EQ(late_dropped, 1u);
}

/**
 * @brief Verifies a master coast past 2^31 cycles is counted and recovered
 *        rather than wedging the flywheel silently.
 * @details fold() reads (now - epoch_cycles) as int32 and reports no crossing on
 *          a negative one; the master re-anchors its flywheel.
 */
inline void test_master_fold_stall_recovers() {
  const Config cfg = test_config();
  const uint32_t period = cfg.cycles_per_half_rev;

  SyncBoard m(cfg);
  const uint32_t t0 = 1000000u;
  m.seed(t0, /*is_master=*/true);
  m.tick(t0 + 2u * period, nullptr);
  HS_EXPECT_EQ(m.telemetry_snapshot().master_stalls, 0u);
  const uint32_t flips_before = m.telemetry_snapshot().flips;
  HS_EXPECT_GT(flips_before, 0u);

  // Resume past the int32 horizon. The epoch trails `now` by more than 2^31, so
  // fold() cannot recover the elapsed crossings in this modular window.
  const uint32_t stalled = t0 + 2u * period + 0x90000000u;
  Flywheel probe(cfg);
  probe.seed(t0 + 2u * period);
  HS_EXPECT_TRUE(probe.fold_stalled(stalled));
  HS_EXPECT_FALSE(probe.fold(stalled).crossed);

  m.tick(stalled, nullptr);
  HS_EXPECT_EQ(m.telemetry_snapshot().master_stalls, 1u);

  // The flywheel is anchored on `stalled` again, so ordinary half-rev wakes
  // resume crossing boundaries.
  m.tick(stalled + period, nullptr);
  m.tick(stalled + 2u * period, nullptr);
  HS_EXPECT_GT(m.telemetry_snapshot().flips, flips_before);
  HS_EXPECT_EQ(m.telemetry_snapshot().master_stalls, 1u);
}

/**
 * @brief Verifies the first boundary crossed after a fold-stall recovery still
 *        flips.
 * @details The re-anchor stamps ZERO whatever the pre-stall identity was, so a
 *          master whose last flip was HALF meets that same identity again on the
 *          next crossing.
 */
inline void test_master_fold_stall_recovery_flips() {
  const Config cfg = test_config();
  const uint32_t period = cfg.cycles_per_half_rev;

  SyncBoard m(cfg);
  const uint32_t t0 = 1000000u;
  m.seed(t0, /*is_master=*/true);
  // A single fold leaves HALF as the last flipped boundary — the identity the
  // re-seed's ZERO reproduces on the very next crossing.
  m.tick(t0 + period, nullptr);
  const uint32_t flips_before = m.telemetry_snapshot().flips;
  HS_EXPECT_EQ(flips_before, 1u);

  const uint32_t stalled = t0 + period + 0x90000000u;
  m.tick(stalled, nullptr);
  HS_EXPECT_EQ(m.telemetry_snapshot().master_stalls, 1u);
  // The re-anchor lands the epoch on `stalled`, so nothing crosses on that tick.
  HS_EXPECT_EQ(m.telemetry_snapshot().flips, flips_before);

  const TickActions a = m.tick(stalled + period, nullptr);
  HS_EXPECT_TRUE(a.flip);
  HS_EXPECT_EQ(m.telemetry_snapshot().flips, flips_before + 1u);
}

// ── Master beacon tail quiet (§6.4) ─────────────────────────────────────────

/**
 * @brief Verifies every beacon the master starts leaves the receiver's quiet
 *        window between the frame's last pulse and the HALF boundary burst.
 * @details Sweeps every column of [W/4, W/2) a masked-ISR coast can resume on,
 *          at the 64-effect roster cap with the widest digit pattern the codec
 *          can encode and again with a narrow one, and requires each emitted
 *          frame's tail to clear acquire_quiet_cycles before the HALF burst.
 *          Ticks run at the device's T0/OVERSAMPLE pacing so the HALF symbol
 *          clears its own lateness censor.
 */
inline void test_beacon_tail_quiet() {
  Config cfg = test_config(64);
  // Revolution 63 is beacon-due at this cadence, so all four data digits reach
  // 7 — the widest frame index 63 of a full roster can encode.
  cfg.beacon_period_revs = 31;
  cfg.rejoin_budget_revs = cfg.rejoin_bound_revs();
  HS_EXPECT_TRUE(cfg.valid() == nullptr);
  const uint32_t period = cfg.cycles_per_half_rev;
  const uint32_t step = COL / 8u;

  // Resume the master's first post-ZERO tick of the beacon-due revolution at
  // `resume_col`, carrying the (index, rev) payload; return the beacon pulse
  // count and report the cycles between the frame's last pulse and the first
  // pulse of the HALF burst.
  auto run = [&](int32_t resume_col, int32_t index, uint32_t rev,
                 uint32_t *tail_gap) -> int {
    SyncBoard m(cfg);
    const uint32_t t0 = 1000000u;
    m.seed(t0, /*is_master=*/true);
    const uint32_t epoch1 = t0 + 2u * period;
    m.tick(epoch1, nullptr);
    // Drain the ZERO symbol's remaining pulses (columns 2, 4) so the emitter is
    // idle before the coast, then dial in the beacon payload.
    m.tick(epoch1 + 2u * COL, nullptr);
    m.tick(epoch1 + 4u * COL, nullptr);
    content_mut(m).effect_index = index;
    content_mut(m).rev_in_effect = rev;

    const uint32_t half_at = epoch1 + period;
    int pulses = 0;
    uint32_t last_beacon = 0;
    *tail_gap = 0;
    for (uint32_t t = epoch1 + static_cast<uint32_t>(resume_col) * COL;
         t <= half_at + 8u * COL; t += step) {
      if (!m.tick(t, nullptr).pulse)
        continue;
      if (t < half_at) {
        ++pulses;
        last_beacon = t;
      } else if (pulses > 0 && *tail_gap == 0) {
        *tail_gap = t - last_beacon;
      }
    }
    return pulses;
  };

  // Revolution 63 of index 63 drives all four data digits to 7 — the widest
  // frame a full roster can encode.
  int widest = 0;
  int narrow = 0;
  for (int32_t c = cfg.W / 4; c < cfg.W / 2; ++c) {
    uint32_t gap = 0;
    if (run(c, 63, 63, &gap) > 0) {
      ++widest;
      HS_EXPECT_GE(gap, cfg.acquire_quiet_cycles());
    }
    gap = 0;
    if (run(c, 0, 1, &gap) > 0) {
      ++narrow;
      HS_EXPECT_GE(gap, cfg.acquire_quiet_cycles());
    }
  }
  // Non-vacuity: the on-time beacon point is still admitted.
  HS_EXPECT_GT(widest, 0);
  // A short payload buys back start columns the worst-case bound censors.
  HS_EXPECT_GT(narrow, widest);
}

// ── Master EPOCH train window (§6.3.1) ──────────────────────────────────────

/**
 * @brief Verifies the master's EPOCH train occupies exactly the R+1 ZERO
 *        boundaries B..B+R even when a copy self-censors, so every copy stays
 *        inside the receiver's invertible j window.
 * @details Drives a lone master to its train-start boundary B and resumes a
 *          full column late there, censoring the primary copy (§5.2). Counts
 *          pulses in each boundary's burst window: 5 = ZERO_EPOCH, 3 = plain
 *          ZERO.
 */
inline void test_master_epoch_train_bounded() {
  const Config cfg = test_config();
  const uint32_t rev_cycles = 2u * cfg.cycles_per_half_rev;
  const uint32_t step = COL / 8;
  const int32_t probed_revs = 9;

  SyncBoard m(cfg);
  const uint32_t t0 = 1000000u;
  m.seed(t0, /*is_master=*/true);
  // Seeded at boundary ZERO with rev_in_effect 0, so the crossing that starts
  // the train (rev_in_effect == revs_per_effect) is exactly this instant.
  const uint32_t b = t0 + cfg.revs_per_effect * rev_cycles;
  for (uint32_t t = t0 + step; t + COL < b; t += step)
    m.tick(t, nullptr);

  // Pulses per ZERO burst window; the beacon point (W/4) and HALF are outside.
  int pulses[probed_revs] = {};
  for (uint32_t t = b + COL;
       t < b + static_cast<uint32_t>(probed_revs) * rev_cycles; t += step) {
    if (!m.tick(t, nullptr).pulse)
      continue;
    for (int32_t k = 0; k < probed_revs; ++k) {
      const uint32_t off = t - (b + static_cast<uint32_t>(k) * rev_cycles);
      if (off < 12u * COL) {
        ++pulses[k];
        break;
      }
    }
  }

  HS_EXPECT_EQ(pulses[0], 0); // primary censored by the late resume
  for (int32_t k = 1; k <= cfg.epoch_repeats; ++k)
    HS_EXPECT_EQ(pulses[k], 5);
  for (int32_t k = cfg.epoch_repeats + 1; k < probed_revs; ++k)
    HS_EXPECT_EQ(pulses[k], 3);
}
