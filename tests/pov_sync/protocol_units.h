/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ── Pure helpers ────────────────────────────────────────────────────────────

inline void expect_rejects(const Config &config, const char *clause) {
  const char *reason = config.valid();
  HS_EXPECT_TRUE(reason != nullptr);
  if (reason)
    HS_EXPECT_EQ(std::strcmp(reason, clause), 0);
}

/**
 * @brief Verifies integer floor-div/mod and circular column distance, plus the
 *        derived Config timing fields (half-rev/column cycle counts, glitch
 *        window) at full rate.
 */
inline void test_helpers() {
  HS_EXPECT_EQ(floor_div(7, 2), 3);
  HS_EXPECT_EQ(floor_div(-7, 2), -4);
  HS_EXPECT_EQ(floor_div(-4, 2), -2);
  HS_EXPECT_EQ(floor_mod(-1, 288), 287);
  HS_EXPECT_EQ(floor_mod(288, 288), 0);
  HS_EXPECT_EQ(floor_mod(int64_t{2147483649}, 288), 129);
  HS_EXPECT_EQ(floor_mod(int64_t{-2147483649}, 288), 159);
  HS_EXPECT_EQ(circ_dist(287, 0, 288), 1);
  HS_EXPECT_EQ(circ_dist(0, 144, 288), 144);
  HS_EXPECT_EQ(circ_dist(10, 280, 288), 18);

  const Config c = test_config();
  HS_EXPECT_TRUE(c.valid() == nullptr);
  HS_EXPECT_EQ(c.cycles_per_half_rev, PERIOD); // 600e6·30/480, exact
  HS_EXPECT_EQ(c.cycles_per_column(), COL);
  HS_EXPECT_EQ(c.glitch_filter_cycles, 60000u); // 100 µs

  // Flywheel::position() carries the elapsed cycle count as int32 across
  // MIN_SAFE_HALF_REVS of coast; valid() rejects a period whose product
  // overflows it. Boundary is inclusive.
  Config pw = test_config();
  pw.cycles_per_half_rev =
      static_cast<uint32_t>(INT32_MAX) / MIN_SAFE_HALF_REVS;
  HS_EXPECT_TRUE(pw.valid() == nullptr);
  ++pw.cycles_per_half_rev;
  expect_rejects(pw, "cycles_per_half_rev * MIN_SAFE_HALF_REVS <= INT32_MAX");

  // The beacon's 6-bit rev field resyncs a slip only in (-32, +32), so the
  // beacon period stays below 32. Boundary is exclusive.
  Config bp = test_config();
  bp.rejoin_budget_revs =
      64; // relax the budget so this isolates the resync bound
  bp.beacon_period_revs = 32;
  expect_rejects(bp, "beacon_period_revs < 32");
  bp.beacon_period_revs = 31;
  HS_EXPECT_TRUE(bp.valid() == nullptr);

  // §9.1 rejoin budget: what must fit is the achieved bound — the widest
  // beacon-to-beacon gap (a cadence plus the commit window's blackout) plus the
  // join-grid wait — not the cadence alone. Boundary is inclusive.
  Config rb = test_config();
  HS_EXPECT_EQ(rb.rejoin_bound_revs(), 17u);
  rb.rejoin_budget_revs = 17;
  HS_EXPECT_TRUE(rb.valid() == nullptr);
  --rb.rejoin_budget_revs;
  expect_rejects(rb, "rejoin_bound_revs() <= rejoin_budget_revs");
  // Every term moves the bound: a 16-rev cadence does not fit a 16-rev budget,
  // and the commit window and join grid each push it out further.
  Config rc = test_config();
  rc.revs_per_effect = 48;
  rc.rejoin_budget_revs = 16;
  rc.beacon_period_revs = 16;
  expect_rejects(rc, "rejoin_bound_revs() <= rejoin_budget_revs");
  rc.rejoin_budget_revs = 25;
  HS_EXPECT_EQ(rc.rejoin_bound_revs(), 25u);
  HS_EXPECT_TRUE(rc.valid() == nullptr);
  ++rc.commit_revs;
  expect_rejects(rc, "rejoin_bound_revs() <= rejoin_budget_revs");
  --rc.commit_revs;
  rc.join_grid_revs = 8;
  expect_rejects(rc, "rejoin_bound_revs() <= rejoin_budget_revs");

  // The commit-window gap depends on the entry's duration mod the cadence: an
  // entry one revolution past a multiple of it spends a full cadence after its
  // last beacon before the epoch.
  Config rd = test_config();
  rd.beacon_period_revs = 16;
  rd.rejoin_budget_revs = 25;
  rd.revs_per_effect = 50;
  HS_EXPECT_EQ(rd.commit_beacon_gap_revs(0), 7u);
  HS_EXPECT_EQ(rd.rejoin_bound_revs(), 20u);
  rd.revs_per_effect = 49;
  HS_EXPECT_EQ(rd.commit_beacon_gap_revs(0), 22u);
  expect_rejects(rd, "rejoin_bound_revs() <= rejoin_budget_revs");
  // One such entry in a roster is enough.
  const uint32_t mixed[4] = {48, 64, 49, 32};
  rd.set_effect_revolutions(mixed);
  expect_rejects(rd, "rejoin_bound_revs() <= rejoin_budget_revs");
  // The boot repeats are the last beacon when the cadence never fires.
  rd.revs_per_effect = 17;
  rd.clear_effect_revolutions();
  HS_EXPECT_EQ(rd.commit_beacon_gap_revs(0), 20u);

  // Demarcation: a wire timeout below the beacon's worst-case per-digit advance
  // splits a real digit train into isolated boundary symbols. The shipped
  // constants meet the bound with equality.
  Config aq = test_config();
  const int32_t aq_bound =
      2 * aq.gap_timeout_cols + 7 * aq.beacon_pitch_cols + 1;
  HS_EXPECT_EQ(aq.acquire_quiet_cols, aq_bound);
  aq.acquire_quiet_cols = aq_bound - 1;
  expect_rejects(aq,
                 "acquire_quiet_cols >= 2*gap_timeout + 7*beacon_pitch + 1");
  ++aq.acquire_quiet_cols;
  HS_EXPECT_TRUE(aq.valid() == nullptr);

  // Stale-frame window order: tick()'s poll-path parser reset, which fires at
  // acquire_quiet_cols + gap_timeout_cols, must be tighter than
  // BeaconParser::feed's interdigit test. Boundary is exclusive (equal windows
  // are rejected).
  Config so = test_config();
  so.beacon_interdigit_timeout_cols =
      so.acquire_quiet_cols + so.gap_timeout_cols;
  expect_rejects(
      so, "acquire_quiet_cols + gap_timeout < beacon_interdigit_timeout");
  ++so.beacon_interdigit_timeout_cols;
  HS_EXPECT_TRUE(so.valid() == nullptr);

  // Beacon tail quiet: a frame started at the beacon point must land its last
  // pulse and the receiver's quiet window inside the [W/4, W/2) half-window, or
  // the HALF burst is appended to the last digit burst. Boundary is exclusive —
  // the slack absorbs the sub-column offset of the scheduling tick.
  Config bq = test_config();
  bq.acquire_quiet_cols = bq.W / 4 - bq.beacon_span_cols();
  HS_EXPECT_EQ(bq.beacon_frame_cols(), bq.W / 4);
  expect_rejects(bq, "beacon_frame_cols() < W/4");
  --bq.acquire_quiet_cols;
  HS_EXPECT_TRUE(bq.valid() == nullptr);

  Config dg = test_config();
  dg.gate_cols = 7 * dg.beacon_pitch_cols + 1;
  expect_rejects(dg, "7*beacon_pitch_cols + 1 > gate_cols");
  --dg.gate_cols;
  HS_EXPECT_TRUE(dg.valid() == nullptr);

  // A glitch filter at or above a burst's pulse spacing, less the emitter's
  // lateness budget, drops every pulse after the first. Boundary is exclusive.
  Config gf = test_config();
  const uint32_t gf_bound = gf.beacon_pitch_cycles() - gf.late_censor_cycles();
  HS_EXPECT_TRUE(gf.glitch_filter_cycles < gf_bound);
  HS_EXPECT_TRUE(gf.glitch_filter_cycles <
                 gf.pulse_pitch_cycles() - gf.late_censor_cycles());
  gf.glitch_filter_cycles = gf_bound;
  expect_rejects(gf, "glitch_filter_cycles < beacon_pitch - late_censor");
  gf.glitch_filter_cycles = gf_bound - 1;
  HS_EXPECT_TRUE(gf.valid() == nullptr);
}

/**
 * @brief Probes the remaining Config::valid() clauses at their boundaries.
 * @details Each rejection pins the named first failing clause.
 */
inline void test_config_validation() {
  Config zero_width = test_config();
  zero_width.W = 0;
  expect_rejects(zero_width, "W > 0");
  Config zero_period = test_config();
  zero_period.cycles_per_half_rev = 0;
  expect_rejects(zero_period, "cycles_per_half_rev > 0");
  Config zero_pitch = test_config();
  zero_pitch.beacon_pitch_cols = 0;
  expect_rejects(zero_pitch, "beacon_pitch_cols > 0");

  HS_EXPECT_TRUE(test_config().valid() == nullptr);

  // Odd W: boundary_column(HALF) and every arm-B half-image offset truncate
  // W/2.
  Config ow = test_config();
  ow.W = 289;
  expect_rejects(ow, "W even");

  Config gz = test_config();
  gz.gate_cols = 0;
  expect_rejects(gz, "gate_cols > 0");

  Config gw = test_config();
  gw.gate_cols = gw.W / 4;
  expect_rejects(gw, "gate_cols < W/4");
  --gw.gate_cols;
  expect_rejects(gw, "7*beacon_pitch_cols + 1 > gate_cols");

  Config rj = test_config();
  rj.reject_fallback = 0;
  expect_rejects(rj, "reject_fallback > 0");

  Config gz2 = test_config();
  gz2.glitch_filter_cycles = 0;
  expect_rejects(gz2, "glitch_filter_cycles > 0");

  Config pz = test_config();
  pz.pulse_pitch_cols = 0;
  expect_rejects(pz, "pulse_pitch_cols > 0");

  // A gap that does not outlast the boundary-burst pitch terminates the burst
  // between its own pulses. Boundary is exclusive.
  Config gp = test_config();
  gp.gap_timeout_cols = gp.pulse_pitch_cols;
  expect_rejects(gp, "gap_timeout_cols > pulse_pitch_cols");
  ++gp.gap_timeout_cols;
  HS_EXPECT_TRUE(gp.valid() == nullptr);

  // Same ordering against the beacon pitch. At the shipped W the pulse-pitch
  // clause always binds first, so isolate this one on a canvas wide enough to
  // carry a beacon pitch above the boundary pitch.
  Config bg = test_config();
  bg.W = 1024;
  bg.glitch_filter_cycles = 10000;
  bg.pulse_pitch_cols = 1;
  bg.beacon_pitch_cols = 3;
  bg.acquire_quiet_cols = 32;
  bg.beacon_interdigit_timeout_cols = 40;
  bg.gap_timeout_cols = 4;
  HS_EXPECT_TRUE(bg.valid() == nullptr);
  Config pulse_glitch = bg;
  pulse_glitch.glitch_filter_cycles =
      pulse_glitch.pulse_pitch_cycles() - pulse_glitch.late_censor_cycles();
  expect_rejects(pulse_glitch,
                 "glitch_filter_cycles < pulse_pitch - late_censor");
  --pulse_glitch.glitch_filter_cycles;
  HS_EXPECT_TRUE(pulse_glitch.valid() == nullptr);
  bg.gap_timeout_cols = bg.beacon_pitch_cols;
  expect_rejects(bg, "gap_timeout_cols > beacon_pitch_cols");

  // Epoch indices are taken mod effect_count and ride a 6-bit beacon field.
  Config ec = test_config();
  ec.effect_count = 0;
  expect_rejects(ec, "effect_count > 0");
  ec.effect_count = 65;
  expect_rejects(ec, "effect_count <= 64");
  ec.effect_count = 64;
  HS_EXPECT_TRUE(ec.valid() == nullptr);

  Config cz = test_config();
  cz.commit_revs = 0;
  expect_rejects(cz, "commit_revs > 0");

  // A negative epoch_repeats casts to a huge uint32_t and wraps the refractory
  // sum, so the relation below it reads true; only the explicit sign gate
  // rejects the config.
  Config ne = test_config();
  ne.epoch_repeats = -1;
  HS_EXPECT_TRUE(ne.refractory_revs >
                 ne.commit_revs + static_cast<uint32_t>(ne.epoch_repeats));
  expect_rejects(ne, "epoch_repeats >= 0");
  ne.epoch_repeats = 0;
  HS_EXPECT_TRUE(ne.valid() == nullptr);

  // The EPOCH dedup window must outlast the whole redundancy train, or the
  // train's own last repeat re-arms the commit it just deduped. Boundary is
  // exclusive.
  Config rf = test_config();
  rf.refractory_revs = rf.commit_revs + static_cast<uint32_t>(rf.epoch_repeats);
  expect_rejects(rf, "refractory_revs > commit_revs + epoch_repeats");
  ++rf.refractory_revs;
  HS_EXPECT_TRUE(rf.valid() == nullptr);

  // An effect shorter than the dedup window would advance before the window
  // that protects its own commit closes. Boundary is exclusive.
  Config re = test_config();
  re.revs_per_effect = re.refractory_revs;
  expect_rejects(re, "revolutions_for_effect(i) > refractory_revs");
  ++re.revs_per_effect;
  HS_EXPECT_TRUE(re.valid() == nullptr);

  uint32_t variable_revolutions[4] = {40, 48, 56, 64};
  Config vr = test_config();
  vr.effect_revolutions = variable_revolutions;
  vr.effect_revolutions_count = 0;
  expect_rejects(vr, "effect_revolutions_count >= effect_count");
  vr.set_effect_revolutions(variable_revolutions);
  HS_EXPECT_TRUE(vr.valid() == nullptr);
  variable_revolutions[2] = vr.refractory_revs;
  expect_rejects(vr, "revolutions_for_effect(i) > refractory_revs");

  // A beacon cadence inside the construction window would land identity
  // traffic on the commit boundary. Boundary is exclusive.
  Config bc = test_config();
  bc.beacon_period_revs = bc.commit_revs;
  expect_rejects(bc, "beacon_period_revs > commit_revs");
  ++bc.beacon_period_revs;
  HS_EXPECT_TRUE(bc.valid() == nullptr);

  // The live-takeover grid must divide 64 so a beacon's mod-64 revolution
  // count lands on the same grid as the master's true count.
  Config jg = test_config();
  jg.join_grid_revs = 0;
  expect_rejects(jg, "join_grid_revs > 0");
  jg.join_grid_revs = 3;
  expect_rejects(jg, "join_grid_revs divides 64");
  jg.join_grid_revs = 8;
  HS_EXPECT_TRUE(jg.valid() == nullptr);
}

/**
 * @brief Verifies symbol/boundary mapping: odd pulse counts classify to
 *        HALF/ZERO/ZERO_EPOCH, even or out-of-range counts are INVALID, and
 *        each symbol maps to its boundary column.
 */
inline void test_alphabet() {
  HS_EXPECT_EQ(classify_count(1), Symbol::HALF);
  HS_EXPECT_EQ(classify_count(3), Symbol::ZERO);
  HS_EXPECT_EQ(classify_count(5), Symbol::ZERO_EPOCH);
  // Even counts (single lost/spurious edge) and out-of-range are INVALID.
  for (uint32_t n : {0u, 2u, 4u, 6u, 7u, 8u, 9u, 255u})
    HS_EXPECT_EQ(classify_count(n), Symbol::INVALID);
  HS_EXPECT_EQ(symbol_boundary(Symbol::HALF), Boundary::HALF);
  HS_EXPECT_EQ(symbol_boundary(Symbol::ZERO), Boundary::ZERO);
  HS_EXPECT_EQ(symbol_boundary(Symbol::ZERO_EPOCH), Boundary::ZERO);
  HS_EXPECT_EQ(symbol_pulse_count(Symbol::ZERO_EPOCH), 5u);
  HS_EXPECT_EQ(boundary_column(Boundary::HALF, 288), 144);
  HS_EXPECT_EQ(boundary_column(Boundary::ZERO, 288), 0);
}

/**
 * @brief Verifies §5.1 exactly-once flipping across interleaved
 *        crossing/symbol arrivals.
 */
inline void test_flip_gate() {
  FlipGate g;
  HS_EXPECT_FALSE(g.try_flip(Boundary::NONE));
  HS_EXPECT_TRUE(g.try_flip(Boundary::HALF));  // boot: HALF != NONE flips
  HS_EXPECT_FALSE(g.try_flip(Boundary::HALF)); // symbol after crossing: dedup
  HS_EXPECT_TRUE(g.try_flip(Boundary::ZERO));
  HS_EXPECT_FALSE(g.try_flip(Boundary::ZERO));
  HS_EXPECT_TRUE(g.try_flip(Boundary::HALF));
  // Symbol-leads interleaving: symbol flips, late crossing dedups, next
  // boundary flips again — exactly 2 per revolution.
  FlipGate h;
  int flips = 0;
  for (int rev = 0; rev < 5; ++rev) {
    flips += h.try_flip(Boundary::ZERO); // symbol
    flips += h.try_flip(Boundary::ZERO); // crossing (dedup)
    flips += h.try_flip(Boundary::HALF); // crossing
    flips += h.try_flip(Boundary::HALF); // symbol (dedup)
  }
  HS_EXPECT_EQ(flips, 10);
}

/**
 * @brief Verifies edge mailbox burst accumulation and the glitch filter:
 *        sub-window spikes are rejected without resetting the filter
 *        reference, burst_complete fires only after the gap, and claim()
 *        snapshots count + first/last edge.
 */
inline void test_mailbox() {
  const uint32_t GLITCH = 60000u;
  EdgeMailbox m;
  HS_EXPECT_FALSE(burst_complete(m, 0, 4 * COL));
  m.on_edge(1000, GLITCH);
  m.on_edge(1000 + 2 * COL, GLITCH);
  m.on_edge(1000 + 4 * COL, GLITCH);
  // EMI spike inside the glitch window after an accepted edge is rejected…
  m.on_edge(1000 + 4 * COL + GLITCH / 2, GLITCH);
  m.on_edge(1000 + 4 * COL + GLITCH + 1, GLITCH);
  m.on_edge(1000 + 6 * COL, GLITCH);
  HS_EXPECT_FALSE(burst_complete(m, 1000 + 7 * COL, 4 * COL));
  HS_EXPECT_TRUE(burst_complete(m, 1000 + 10 * COL + 1, 4 * COL));
  const BurstSnapshot s = claim(m);
  HS_EXPECT_EQ(s.count, 5u);
  HS_EXPECT_EQ(s.first_cycles, 1000u);
  HS_EXPECT_EQ(s.last_cycles, 1000u + 6 * COL);
  // Claim resets; the glitch filter still applies across bursts.
  HS_EXPECT_FALSE(burst_complete(m, 1000 + 10 * COL + 2, 4 * COL));
  m.on_edge(1000 + 6 * COL + GLITCH - 1, GLITCH); // too close: rejected
  HS_EXPECT_FALSE(burst_complete(m, 1000 + 20 * COL, 4 * COL));

  EdgeMailbox tc;
  BurstSnapshot out;
  const uint32_t MAXB = test_config().max_burst_cycles();
  HS_EXPECT_FALSE(tc.try_claim(0, 4 * COL, MAXB, &out)); // no burst yet
  tc.on_edge(1000, GLITCH);
  tc.on_edge(1000 + 2 * COL, GLITCH);
  HS_EXPECT_FALSE(
      tc.try_claim(1000 + 3 * COL, 4 * COL, MAXB, &out)); // gap too short
  HS_EXPECT_TRUE(tc.try_claim(1000 + 6 * COL + 1, 4 * COL, MAXB, &out));
  HS_EXPECT_EQ(out.count, 2u);
  HS_EXPECT_EQ(out.first_cycles, 1000u);
  HS_EXPECT_EQ(out.last_cycles, 1000u + 2 * COL);
  HS_EXPECT_FALSE(tc.try_claim(1000 + 7 * COL, 4 * COL, MAXB, &out)); // reset
}

/**
 * @brief Verifies a burst the wire never lets go quiet is claimed on duration.
 * @details Noise at the glitch filter's pass rate refreshes last_cycles on every
 *          accepted edge, so the terminating gap never opens.
 */
inline void test_mailbox_overlong_burst() {
  const Config cfg = test_config();
  const uint32_t GLITCH = cfg.glitch_filter_cycles;
  const uint32_t MAXB = cfg.max_burst_cycles();
  const uint32_t GAP = cfg.gap_timeout_cycles();
  const uint32_t t0 = 1000u;

  EdgeMailbox m;
  BurstSnapshot out{};
  uint32_t claimed_at = 0;
  bool claimed = false;
  for (uint32_t t = t0; t <= t0 + 4 * MAXB; t += GLITCH) {
    m.on_edge(t, GLITCH);
    if (!claimed && m.try_claim(t, GAP, MAXB, &out)) {
      claimed = true;
      claimed_at = t;
    }
  }
  HS_EXPECT_TRUE(claimed);
  HS_EXPECT_LE(claimed_at - t0, MAXB + GLITCH);
  HS_EXPECT_EQ(out.first_cycles, t0);
  // A jammed burst never lands on the odd-only alphabet, so the symbol is
  // discarded and counted rather than snapped on.
  HS_EXPECT_EQ(classify_count(out.count), Symbol::INVALID);

  // The widest legitimate burst still terminates on the gap, not the bound.
  EdgeMailbox n;
  const uint32_t pitch = cfg.pulse_pitch_cycles();
  for (uint32_t i = 0; i < 5; ++i)
    n.on_edge(t0 + i * pitch, GLITCH);
  HS_EXPECT_FALSE(n.try_claim(t0 + 4 * pitch + GAP - 1, GAP, MAXB, &out));
  HS_EXPECT_TRUE(n.try_claim(t0 + 4 * pitch + GAP, GAP, MAXB, &out));
  HS_EXPECT_EQ(out.count, 5u);
}

/**
 * @brief Verifies the glitch-filter reference does not survive a counter wrap.
 * @details age_prior() retires the prior once the wire is quiet past the filter
 *          window, so a real edge after wrap is not rejected by a pseudo-random
 *          modular difference.
 */
inline void test_mailbox_prior_staleness() {
  const uint32_t GLITCH = 60000u;

  // age_prior leaves a within-window reference intact: a genuine spike that
  // arrives before a poll has aged the prior is still suppressed.
  {
    EdgeMailbox m;
    m.on_edge(1000, GLITCH);
    m.age_prior(1000 + GLITCH / 2, GLITCH);   // still within the window: kept
    m.on_edge(1000 + GLITCH / 2 + 1, GLITCH); // too close to 1000: rejected
    HS_EXPECT_TRUE(burst_complete(m, 1000 + 100 * GLITCH, 1));
    HS_EXPECT_EQ(claim(m).count, 1u);
  }

  // After the wire goes quiet the prior is retired, so a later edge whose
  // (wrapped) modular distance to the OLD prior lands inside the reject window
  // is still accepted.
  {
    EdgeMailbox m;
    const uint32_t prior = 1000u;
    m.on_edge(prior, GLITCH); // a one-edge burst…
    HS_EXPECT_TRUE(burst_complete(m, prior + 10 * COL, COL));
    HS_EXPECT_EQ(claim(m).count, 1u); // …claimed; the prior persists across it.
    // A poll during the silence retires the stale reference (COL > GLITCH).
    m.age_prior(prior + 11 * COL, GLITCH);
    // A real edge after the counter has wrapped: its modular distance to the
    // old prior is only GLITCH/2, which the un-aged filter would reject.
    const uint32_t wrapped = prior + GLITCH / 2;
    m.on_edge(wrapped, GLITCH);
    HS_EXPECT_TRUE(burst_complete(m, wrapped + 10 * COL, COL));
    HS_EXPECT_EQ(claim(m).count, 1u); // accepted as a fresh one-edge burst.
  }
}

/**
 * @brief Verifies the consumer's gap tests reject a clock sampled before an
 *        edge the publisher went on to accept.
 * @details A sync edge accepted after `now` was sampled leaves the mailbox
 *          timestamps ahead of `now`, and the unsigned difference underflows to
 *          ~2³².
 */
inline void test_mailbox_rejects_backward_clock() {
  const uint32_t GLITCH = 60000u;
  const uint32_t now = 100000u;
  const uint32_t skew = 40u; // edge accepted after `now` was sampled

  EdgeMailbox m;
  BurstSnapshot out{};
  m.on_edge(now + skew, GLITCH);
  HS_EXPECT_FALSE(burst_complete(m, now, 4 * COL));
  HS_EXPECT_FALSE(
      m.try_claim(now, 4 * COL, test_config().max_burst_cycles(), &out));

  m.age_prior(now, GLITCH);
  // The reference survived, so a spike inside the window is still suppressed.
  m.on_edge(now + skew + GLITCH / 2, GLITCH);
  HS_EXPECT_TRUE(burst_complete(m, now + skew + 100 * GLITCH, 1));
  HS_EXPECT_EQ(claim(m).count, 1u);
}

/**
 * @brief Verifies reboot seeding clears the wire mailbox so a re-seeded board
 *        cannot consume a stale pre-reboot burst.
 */
inline void test_seed_clears_mailbox() {
  const Config cfg = test_config();
  const uint32_t col = cfg.cycles_per_column();
  SyncBoard board(cfg);
  board.seed(1000u, /*is_master=*/false);
  board.on_sync_edge(2000u);
  board.on_sync_edge(2000u + 4 * col);
  board.seed(3000u, false);
  BurstSnapshot s;
  HS_EXPECT_FALSE(board.claim_sync_burst(3000u + 100 * col, &s));
  board.on_sync_edge(3000u + 200 * col);
  HS_EXPECT_TRUE(board.claim_sync_burst(3000u + 300 * col, &s));
  HS_EXPECT_EQ(s.count, 1u);
}

/**
 * @brief Verifies build requests reset across seeds and reconfiguration.
 */
inline void test_build_request_reset() {
  const Config cfg = test_config();
  SyncBoard board(cfg);
  HS_EXPECT_EQ(board.build_word(), 0u);

  board.seed(1000u, true);
  HS_EXPECT_EQ(SyncBoard::build_gen_of(board.build_word()), 1u);
  HS_EXPECT_EQ(SyncBoard::build_index_of(board.build_word()), 0);

  board.seed(2000u, false);
  HS_EXPECT_EQ(board.build_word(), 0u);

  Config replacement = test_config(3);
  replacement.cycles_per_half_rev /= 2;
  replacement.glitch_filter_cycles /= 2;
  board.configure(replacement);
  HS_EXPECT_EQ(config(board).effect_count, 3);
  HS_EXPECT_EQ(board.build_word(), 0u);

  board.seed(3000u, true);
  HS_EXPECT_EQ(SyncBoard::build_gen_of(board.build_word()), 1u);
}

/**
 * @brief Verifies burst claiming follows the gap, duration and glitch windows
 *        of the configuration installed by configure().
 */
inline void test_configure_replaces_claim_windows() {
  const Config cfg = test_config();
  Config replacement = test_config(3);
  replacement.cycles_per_half_rev /= 2;
  replacement.glitch_filter_cycles /= 2;
  const uint32_t gap = replacement.gap_timeout_cycles();
  const uint32_t col = replacement.cycles_per_column();
  const uint32_t t = 5000u;
  BurstSnapshot s;
  {
    SyncBoard board(cfg);
    board.configure(replacement);
    board.on_sync_edge(t);
    HS_EXPECT_FALSE(board.claim_sync_burst(t + gap - 1, &s));
    HS_EXPECT_TRUE(board.claim_sync_burst(t + gap, &s));
    HS_EXPECT_EQ(s.count, 1u);
  }
  {
    SyncBoard board(cfg);
    board.configure(replacement);
    uint32_t edge = t;
    while (edge - t < replacement.max_burst_cycles()) {
      board.on_sync_edge(edge);
      edge += col;
    }
    board.on_sync_edge(edge);
    HS_EXPECT_TRUE(board.claim_sync_burst(edge + 1, &s));
    HS_EXPECT_EQ(s.first_cycles, t);
  }
  {
    SyncBoard board(cfg);
    board.configure(replacement);
    board.on_sync_edge(t);
    board.on_sync_edge(t + replacement.glitch_filter_cycles);
    HS_EXPECT_TRUE(board.claim_sync_burst(t + 100 * col, &s));
    HS_EXPECT_EQ(s.count, 2u);
  }
}

/**
 * @brief Verifies a tick folding several boundaries reports the FINAL one:
 *        TickActions::zero_crossing names the boundary that opened the display
 *        window now in effect, which the driver publishes as the window half.
 * @details The gate runs once per crossing so the flip counter and the coast
 *          telemetry see all N, while `flip` stays a single bool (§5.1).
 */
inline void test_multi_boundary_tick_window() {
  const Config cfg = test_config();
  SyncBoard board(cfg);
  board.seed(1000u, /*is_master=*/false);
  flywheel_mut(board).force_lock();

  // Coast past three boundaries in one wake: HALF, ZERO, HALF.
  const TickActions a =
      board.tick(1000u + 3u * cfg.cycles_per_half_rev + 10u, nullptr);
  HS_EXPECT_TRUE(a.flip);
  HS_EXPECT_FALSE(a.zero_crossing);
  HS_EXPECT_EQ(flywheel(board).current_boundary(), Boundary::HALF);
  HS_EXPECT_EQ(board.telemetry_snapshot().flips, 3u);
  HS_EXPECT_EQ(board.telemetry_snapshot().max_coast_halves, 3u);

  // The next boundary is a ZERO, and it is reported.
  const TickActions z =
      board.tick(1000u + 4u * cfg.cycles_per_half_rev + 10u, nullptr);
  HS_EXPECT_TRUE(z.flip);
  HS_EXPECT_TRUE(z.zero_crossing);
  HS_EXPECT_EQ(flywheel(board).current_boundary(), Boundary::ZERO);
  HS_EXPECT_EQ(board.telemetry_snapshot().flips, 4u);

  // An even fold returns to the half it started from: HALF then ZERO reports
  // the ZERO, so the published window half still names the open window.
  const TickActions e =
      board.tick(1000u + 6u * cfg.cycles_per_half_rev + 10u, nullptr);
  HS_EXPECT_TRUE(e.flip);
  HS_EXPECT_TRUE(e.zero_crossing);
  HS_EXPECT_EQ(flywheel(board).current_boundary(), Boundary::ZERO);
  HS_EXPECT_EQ(board.telemetry_snapshot().flips, 6u);
  HS_EXPECT_EQ(board.telemetry_snapshot().max_coast_halves, 6u);
}

/**
 * @brief Verifies §6.4 beacon codec: frames round-trip, and corrupted frames
 *        are dropped whole, never partially applied.
 */
inline void test_beacon_codec() {
  const Config cfg = test_config();
  uint8_t d[5];
  encode_beacon_digits(27, 45, d);
  HS_EXPECT_EQ(d[0], 3);
  HS_EXPECT_EQ(d[1], 3);
  HS_EXPECT_EQ(d[2], 5);
  HS_EXPECT_EQ(d[3], 5);
  HS_EXPECT_EQ(d[4], 0);

  /**
   * @brief Feeds all five digit bursts through a fresh parser.
   * @param digits The five encoded beacon digits.
   * @param out Receives the decoded frame on success.
   * @return True iff the frame decoded with no rejection.
   */
  auto feed_frame = [&cfg](const uint8_t digits[5], BeaconFrame *out) {
    BeaconParser parser;
    bool got = false;
    for (int i = 0; i < 5; ++i) {
      const uint32_t FIRST = 1000u + i * 12u * COL;
      const BurstSnapshot burst{static_cast<uint32_t>(digits[i]) + 1u, FIRST,
                                FIRST + digits[i] * COL};
      bool rejected = false;
      got = parser.feed(burst, cfg, out, &rejected);
      if (i < 4) {
        HS_EXPECT_FALSE(got);
        HS_EXPECT_FALSE(rejected);
      } else {
        HS_EXPECT_EQ(rejected, !got);
      }
    }
    HS_EXPECT_EQ(parser.digit_count(), 0);
    return got;
  };

  for (int idx = 0; idx < 64; ++idx) {
    for (uint32_t rev = 0; rev < 64; ++rev) {
      encode_beacon_digits(idx, rev, d);
      BeaconFrame frame{};
      HS_EXPECT_TRUE(feed_frame(d, &frame));
      HS_EXPECT_EQ(frame.effect_index, idx);
      HS_EXPECT_EQ(frame.rev_count, rev);
      for (int position = 0; position < 5; ++position) {
        const uint8_t ORIGINAL = d[position];
        for (uint8_t replacement = 0; replacement < 8; ++replacement) {
          if (replacement == ORIGINAL)
            continue;
          d[position] = replacement;
          BeaconFrame unchanged{123, 456u};
          HS_EXPECT_FALSE(feed_frame(d, &unchanged));
          HS_EXPECT_EQ(unchanged.effect_index, 123);
          HS_EXPECT_EQ(unchanged.rev_count, 456u);
        }
        d[position] = ORIGINAL;
      }
    }
  }

  encode_beacon_digits(0, 945, d);
  HS_EXPECT_EQ(d[2], 6);
  HS_EXPECT_EQ(d[3], 1);
  HS_EXPECT_EQ(d[4], 7);

  // Out-of-range burst count aborts the frame.
  {
    BeaconParser p;
    BeaconFrame g{};
    bool r = false;
    BurstSnapshot s1{3, 1000, 1000 + 2 * COL};
    HS_EXPECT_FALSE(p.feed(s1, cfg, &g, &r));
    BurstSnapshot s2{9, 1000 + 12 * COL, 1000 + 20 * COL}; // count > 8
    HS_EXPECT_FALSE(p.feed(s2, cfg, &g, &r));
    HS_EXPECT_TRUE(r);
    HS_EXPECT_EQ(p.digit_count(), 0);
  }

  // A stale partial frame (interdigit timeout) is discarded; the burst that
  // exposed it starts a fresh frame which still decodes.
  {
    BeaconParser p;
    BeaconFrame g{};
    bool r = false;
    BurstSnapshot s1{4, 1000, 1000 + 3 * COL};
    p.feed(s1, cfg, &g, &r);
    encode_beacon_digits(5, 7, d);
    uint32_t t = 1000 + 200 * COL; // far past the interdigit timeout
    bool got = false;
    bool any_reject = false;
    for (int i = 0; i < 5; ++i) {
      BurstSnapshot s{static_cast<uint32_t>(d[i]) + 1u, t, t + d[i] * COL};
      bool rr = false;
      got = p.feed(s, cfg, &g, &rr);
      any_reject = any_reject || (rr && i > 0);
      t += 12 * COL;
    }
    HS_EXPECT_TRUE(got);
    HS_EXPECT_FALSE(any_reject);
    HS_EXPECT_EQ(g.effect_index, 5);
    HS_EXPECT_EQ(g.rev_count, 7u);
  }
}

/**
 * @brief Verifies wire silence past the ACQUIRE quiet window discards a partial
 *        beacon frame, so the next train assembles from its own digits alone.
 * @details The parser's own staleness test is a modular difference that a
 *          cycle-counter wrap defeats, so the reset is tick-driven.
 */
inline void test_beacon_partial_frame_ages_out() {
  const Config cfg = test_config();
  const uint32_t col = cfg.cycles_per_column();
  SyncBoard board(cfg);
  board.seed(1000u, /*is_master=*/false);
  flywheel_mut(board).force_lock();

  // A truncated train: two data bursts from column 40 on (far from both
  // boundaries, so the demarcation routes them to the parser), then silence.
  const uint32_t head = 1000u + 40u * col;
  const BurstSnapshot first{4, head, head + 3 * col};
  board.tick(head + 7 * col, &first);
  const uint32_t second_head = head + 8 * col;
  const BurstSnapshot second{4, second_head, second_head + 3 * col};
  board.tick(second_head + 7 * col, &second);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_rejected, 0u);
  board.tick(head + 40 * col, nullptr); // quiet past the ACQUIRE guard
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_rejected, 1u);

  // A complete frame from column 90 on.
  uint8_t d[5];
  encode_beacon_digits(2, 3, d);
  feed_beacon_train(board, cfg, col, 1000u + 90u * col, d);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_ok, 1u);
  // The aged-out partial is a dropped frame and counts as one.
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_rejected, 1u);
  HS_EXPECT_EQ(content(board).effect_index, 2);
}

/**
 * @brief Verifies the §6.3.4 confirmation rule: a live board changes effect
 *        index only after two consecutive beacons name the same one.
 * @details A leading stray burst shifts the frame; one of eight intruder
 *          values satisfies XOR parity.
 */
inline void test_beacon_shift_needs_confirmation() {
  // Full 64-roster: every 6-bit index is in range, so only the confirmation
  // rule can hold the shifted frame back.
  const Config cfg = test_config(64);
  const uint32_t col = cfg.cycles_per_column();
  SyncBoard board(cfg);
  board.seed(1000u, /*is_master=*/false);
  flywheel_mut(board).force_lock();
  content_mut(board).identity_known = true;
  content_mut(board).effect_index = 1;
  content_mut(board).rev_in_effect = 5;

  // Half-rev k's frame starts at column 40, after its crossings are folded so
  // the frame's rev digits match the board's count.
  int half_rev = 0;
  auto open_frame = [&]() {
    const uint32_t t =
        1000u + static_cast<uint32_t>(half_rev++) * cfg.cycles_per_half_rev +
        40u * col;
    board.tick(t, nullptr);
    return t;
  };
  auto feed_digits = [&](uint32_t start, const uint8_t d[5]) {
    feed_beacon_train(board, cfg, col, start, d);
  };
  auto feed_frame = [&](int32_t index) {
    const uint32_t start = open_frame();
    uint8_t d[5];
    encode_beacon_digits(index, content(board).rev_in_effect, d);
    feed_digits(start, d);
  };

  // The shift: [emi, d0, d1, d2, d3] with the real d3 read as the checksum.
  const uint32_t shift_start = open_frame();
  uint8_t truth[5];
  encode_beacon_digits(1, content(board).rev_in_effect, truth);
  int passing = 0;
  uint8_t shifted[5] = {};
  for (uint8_t emi = 0; emi < 8; ++emi) {
    const uint8_t s[5] = {emi, truth[0], truth[1], truth[2], truth[3]};
    if ((s[0] ^ s[1] ^ s[2] ^ s[3]) != s[4])
      continue;
    ++passing;
    for (int i = 0; i < 5; ++i)
      shifted[i] = s[i];
  }
  HS_EXPECT_EQ(passing, 1);
  const int32_t shifted_index = shifted[0] * 8 + shifted[1];
  HS_EXPECT_TRUE(shifted_index != content(board).effect_index);
  HS_EXPECT_TRUE(shifted_index < cfg.effect_count);

  // A lone shifted frame decodes but must not change what is displayed.
  feed_digits(shift_start, shifted);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_ok, 1u);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_rejected, 0u);
  HS_EXPECT_EQ(content(board).effect_index, 1);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacon_index_corrections, 0u);
  HS_EXPECT_EQ(board.build_word(), 0u); // no rebuild published

  // A good frame agreeing with the displayed index clears the candidate.
  feed_frame(1);
  HS_EXPECT_EQ(content(board).effect_index, 1);

  feed_frame(shifted_index);
  HS_EXPECT_EQ(content(board).effect_index, 1);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacon_index_corrections, 0u);
  HS_EXPECT_EQ(board.build_word(), 0u);

  // Two frames naming *different* indices confirm nothing either.
  feed_frame(3);
  HS_EXPECT_EQ(content(board).effect_index, 1);
  feed_frame(2);
  HS_EXPECT_EQ(content(board).effect_index, 1);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacon_index_corrections, 0u);

  // The second agreeing frame applies: index adopted, rebuild published.
  feed_frame(2);
  HS_EXPECT_EQ(content(board).effect_index, 2);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacon_index_corrections, 1u);
  HS_EXPECT_EQ(SyncBoard::build_index_of(board.build_word()), 2);
  HS_EXPECT_EQ(SyncBoard::build_gen_of(board.build_word()), 1u);
  // The whole scenario rode the beacon path: nothing reached the snap gate.
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_ok, 6u);
  HS_EXPECT_EQ(board.telemetry_snapshot().symbols_accepted, 0u);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacon_rev_mismatches, 0u);
  HS_EXPECT_EQ(lock(board), LockState::LOCKED);
}

/**
 * @brief Verifies a checksum-valid beacon naming an index past the roster is
 *        dropped whole (§6.4 integrity by rejection).
 * @details An out-of-roster index can still satisfy the 3-bit checksum; the
 *          frame counts as rejected, like any corrupt frame.
 */
inline void test_beacon_out_of_range_index_rejected() {
  const Config cfg = test_config(); // 4 effects: indices 4..63 are corrupt
  const uint32_t col = cfg.cycles_per_column();
  SyncBoard board(cfg);
  board.seed(1000u, /*is_master=*/false);
  flywheel_mut(board).force_lock();

  // Frames start at column 40 of a half-rev: far from both boundaries, so the
  // demarcation routes every burst to the beacon parser.
  auto feed_frame = [&](int32_t index, uint32_t start) {
    uint8_t d[5];
    encode_beacon_digits(index, 3, d);
    feed_beacon_train(board, cfg, col, start, d);
  };

  // Index 9 encodes and checksums cleanly, but the roster ends at 3.
  feed_frame(9, 1000u + 40u * col);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_ok, 0u);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_rejected, 1u);
  HS_EXPECT_FALSE(content(board).identity_known);
  HS_EXPECT_EQ(board.build_word(), 0u); // no rebuild published

  // The next in-range beacon joins normally.
  feed_frame(2, 1000u + cfg.cycles_per_half_rev + 40u * col);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_ok, 1u);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_rejected, 1u);
  HS_EXPECT_TRUE(content(board).identity_known);
  HS_EXPECT_EQ(content(board).effect_index, 2);
  HS_EXPECT_EQ(SyncBoard::build_index_of(board.build_word()), 2);
}

/**
 * @brief Verifies the §6.4 rev cross-check fold (beacon_rev_resync_delta)
 *        resolves the 63↔0 mod-64 seam.
 * @details The fold maps the beacon residue and current rev_in_effect to a
 * signed slip in [-32, 31], restoring the exact residue across the seam in
 * either direction.
 */
inline void test_rev_resync_fold() {
  // Same-side residues subtract directly (no wrap).
  HS_EXPECT_EQ(beacon_rev_resync_delta(5, 5), 0);
  HS_EXPECT_EQ(beacon_rev_resync_delta(7, 5), 2);
  HS_EXPECT_EQ(beacon_rev_resync_delta(5, 7), -2);

  // Across the 63↔0 seam a small slip stays a small SIGNED delta, not a ~64 jump.
  HS_EXPECT_EQ(beacon_rev_resync_delta(0, 62), 2);  // residue 62, beacon 0: +2
  HS_EXPECT_EQ(beacon_rev_resync_delta(63, 1), -2); // residue 1,  beacon 63: -2
  HS_EXPECT_EQ(beacon_rev_resync_delta(2, 63), 3);  // residue 63, beacon 2:  +3

  // Only the low 6 bits of `current` matter (rev_in_effect carries the absolute
  // count): magnitude-130 and -66 boards fold identically to their residues.
  HS_EXPECT_EQ(beacon_rev_resync_delta(0, 130), -2); // residue 2, beacon 0: -2
  HS_EXPECT_EQ(beacon_rev_resync_delta(8, 66), 6);   // residue 2, beacon 8: +6

  // Endpoints of the correctable window: the fold spans exactly [-32, 31].
  HS_EXPECT_EQ(beacon_rev_resync_delta(31, 0), 31); // max
  HS_EXPECT_EQ(beacon_rev_resync_delta(0, 32),
               -32); // min (32 ahead ≡ 32 behind)

  // The applied fold (rev_in_effect + delta, then re-read mod 64) lands on the
  // beacon's residue across the seam.
  for (uint32_t cur : {62u, 63u, 64u, 65u, 130u, 131u}) {
    for (uint32_t beacon = 0; beacon < 64u; ++beacon) {
      const int32_t d = beacon_rev_resync_delta(beacon, cur);
      HS_EXPECT_TRUE(d >= -32 && d <= 31);
      HS_EXPECT_EQ((static_cast<int64_t>(cur) + d) & 63,
                   static_cast<int64_t>(beacon));
    }
  }
}

// ── Flywheel position math (§4.1): 64-bit, rebase rule, wrap, trim ─────────

/**
 * @brief Verifies flywheel column position and fold cadence: zero truncation
 *        drift vs a long-double reference at nominal and ±40 ppm trim, correct
 *        signed-past folding, one crossing per half-rev, and exactness
 *        preserved across thousands of half-rev folds and multiple 32-bit
 *        counter wraps via the rebase rule.
 */
inline void test_flywheel_position() {
  const Config cfg = test_config();

  /**
   * @brief Reference column position x = floor(delta · (W/2) / period).
   * @param delta Cycles elapsed since the epoch.
   * @param period Cycles per half-rev.
   * @return The signed column index (unfolded), computed via floor-division.
   */
  auto ref_cols = [](int64_t delta, uint32_t period) {
    const long double columns = static_cast<long double>(delta) * 144.0L /
                                static_cast<long double>(period);
    return static_cast<int64_t>(std::floor(columns));
  };

  // Position over one half-rev, at nominal and trim-extreme periods
  // (±40 ppm ≈ ±1500 cycles): zero truncation drift vs the reference.
  for (int32_t trim : {0, +1500, -1500}) {
    HS_CONTEXT("trim", trim);
    const uint32_t period = PERIOD + trim;
    Flywheel f(cfg);
    f.set_cycles_per_half_rev(period);
    f.seed(1000000u);
    for (int64_t delta = 0; delta < period; delta += 12347) {
      HS_CONTEXT("delta", delta);
      const int32_t want = static_cast<int32_t>(ref_cols(delta, period) % 288);
      HS_EXPECT_EQ(f.position(1000000u + static_cast<uint32_t>(delta)), want);
    }
  }

  // Signed past: a timestamp slightly before the epoch lands just below W.
  {
    Flywheel f(cfg);
    f.seed(5000000u);
    HS_EXPECT_EQ(f.position(5000000u - COL), 287);
    HS_EXPECT_EQ(f.position(5000000u - 3 * COL - COL / 2), 284);
  }

  // Fold cadence: exactly one crossing per half-rev, boundaries alternate,
  // each at its exact instant; a long coast yields several crossings.
  {
    Flywheel f(cfg);
    f.seed(1000u);
    HS_EXPECT_FALSE(f.fold(1000u + PERIOD - 1).crossed);
    const Crossing c1 = f.fold(1000u + PERIOD);
    HS_EXPECT_TRUE(c1.crossed);
    HS_EXPECT_EQ(c1.boundary, Boundary::HALF);
    HS_EXPECT_EQ(c1.at_cycles, 1000u + PERIOD);
    HS_EXPECT_FALSE(f.fold(1000u + PERIOD).crossed);
    // Coast 2.5 half-revs: two more crossings at exact instants.
    const uint32_t late = 1000u + PERIOD + 2 * PERIOD + PERIOD / 2;
    const Crossing c2 = f.fold(late);
    HS_EXPECT_TRUE(c2.crossed && c2.boundary == Boundary::ZERO);
    HS_EXPECT_EQ(c2.at_cycles, 1000u + 2 * PERIOD);
    const Crossing c3 = f.fold(late);
    HS_EXPECT_TRUE(c3.crossed && c3.boundary == Boundary::HALF);
    HS_EXPECT_EQ(c3.at_cycles, 1000u + 3 * PERIOD);
    HS_EXPECT_FALSE(f.fold(late).crossed);
    HS_EXPECT_EQ(f.position(late), 144 + 72);
  }

  // The rebase rule makes the 32-bit wrap unobservable: run thousands of
  // folds across several wraps at nominal and trim-extreme periods; every
  // crossing lands on its boundary column and the epoch stays an exact
  // integer multiple ahead.
  for (int32_t trim : {0, +1500, -1500}) {
    HS_CONTEXT("trim", trim);
    const uint32_t period = PERIOD + trim;
    Flywheel f(cfg);
    f.set_cycles_per_half_rev(period);
    uint32_t t = 0xFFFFFFFFu - period / 3; // wrap almost immediately
    f.seed(t);
    for (int k = 1; k <= 5000; ++k) { // crosses many 32-bit wraps
      HS_CONTEXT("fold", k);
      t += period;
      const Crossing c = f.fold(t);
      HS_EXPECT_TRUE(c.crossed);
      HS_EXPECT_EQ(c.at_cycles, t);
      HS_EXPECT_FALSE(f.fold(t).crossed);
      HS_EXPECT_EQ(f.position(t), boundary_column(c.boundary, 288));
    }
  }
}
