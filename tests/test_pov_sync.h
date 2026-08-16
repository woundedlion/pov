/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Host unit tests for the Phantasm synchronization core (hardware/pov_sync.h)
 * — the spec §12 test plan (docs/specs/phantasm_frame_sync_spec.md).
 *
 * Pure pieces are tested directly: symbol classification, the try_flip state
 * machine, the edge mailbox + glitch filter, the beacon codec, the flywheel's
 * 64-bit position math (incl. cycle-counter wrap and trim extremes), the
 * acceptance gate, and the master emitter's self-censoring.
 *
 * The multi-board scenarios run on a small event-driven simulator: one
 * master plus downstream SyncBoards with per-board crystal offsets, a
 * single-latch masked-IRQ model (an edge during a mask is delayed; two merge
 * — the i.MX RT pin-flag behavior the count coding is designed around),
 * symbol drop windows, EMI injection, a foreground model with effect
 * construction delays, and mid-show reboot.
 */
#pragma once

#include "hardware/pov_handoff.h"
#include "hardware/pov_sync.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <algorithm>
#include <cstdint>
#include <deque>
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
    b.maybe_schedule_beacon(now);
  }
  static bool beacon_done(const SyncBoard &b) { return b.beacon_done_this_rev; }
};

inline Flywheel &flywheel_mut(SyncBoard &b) {
  return SyncBoardTestAccess::flywheel(b);
}
inline ContentTracker &content_mut(SyncBoard &b) {
  return SyncBoardTestAccess::content(b);
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
    const uint32_t span = static_cast<uint32_t>(d[i]) * col;
    const BurstSnapshot s{static_cast<uint32_t>(d[i]) + 1u, f, f + span};
    board.tick(f + span + static_cast<uint32_t>(cfg.gap_timeout_cols) * col,
               &s);
    f += span + static_cast<uint32_t>(cfg.gap_timeout_cols + 1) * col;
  }
  return f;
}

// ── Pure helpers ────────────────────────────────────────────────────────────

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
  HS_EXPECT_EQ(circ_dist(287, 0, 288), 1);
  HS_EXPECT_EQ(circ_dist(0, 144, 288), 144);
  HS_EXPECT_EQ(circ_dist(10, 280, 288), 18);

  const Config c = test_config();
  HS_EXPECT_TRUE(c.valid() == nullptr);
  HS_EXPECT_EQ(c.cycles_per_half_rev, PERIOD); // 600e6·30/480, exact
  HS_EXPECT_EQ(c.cycles_per_column(), COL);
  HS_EXPECT_EQ(c.glitch_filter_cycles, 60000u); // 100 µs

  // The beacon's 6-bit rev field resyncs a slip only in (-32, +32), so the
  // beacon period must stay below the half-window: a period >= 32 leaves the
  // resync precondition unenforced. Boundary is exclusive (32 is rejected, 31
  // accepted).
  Config bp = test_config();
  bp.rejoin_budget_revs =
      64; // relax the budget so this isolates the resync bound
  bp.beacon_period_revs = 32;
  HS_EXPECT_TRUE(bp.valid() != nullptr);
  bp.beacon_period_revs = 31;
  HS_EXPECT_TRUE(bp.valid() == nullptr);

  // §9.1 rejoin budget: what must fit is the achieved bound — the widest
  // beacon-to-beacon gap (a cadence plus the commit window's blackout) plus the
  // join-grid wait — not the cadence alone. Boundary is inclusive.
  Config rb = test_config();
  HS_EXPECT_EQ(rb.rejoin_bound_revs(),
               rb.beacon_period_revs + rb.join_grid_revs +
                   static_cast<uint32_t>(rb.epoch_repeats) + rb.commit_revs);
  rb.rejoin_budget_revs = rb.rejoin_bound_revs();
  HS_EXPECT_TRUE(rb.valid() == nullptr);
  --rb.rejoin_budget_revs;
  HS_EXPECT_TRUE(rb.valid() != nullptr);
  // Every term moves the bound: a 16-rev cadence does not fit a 16-rev budget,
  // and the commit window and join grid each push it out further.
  Config rc = test_config();
  rc.rejoin_budget_revs = 16;
  rc.beacon_period_revs = 16;
  HS_EXPECT_TRUE(rc.valid() != nullptr);
  rc.rejoin_budget_revs = rc.rejoin_bound_revs();
  HS_EXPECT_TRUE(rc.valid() == nullptr);
  ++rc.commit_revs;
  HS_EXPECT_TRUE(rc.valid() != nullptr);
  --rc.commit_revs;
  rc.join_grid_revs = 8;
  HS_EXPECT_TRUE(rc.valid() != nullptr);

  // Demarcation: a wire timeout below the beacon's worst-case per-digit advance
  // splits a real digit train into isolated boundary symbols.
  Config aq = test_config();
  aq.acquire_quiet_cols = aq.beacon_span_cols() / 4 - 1;
  HS_EXPECT_TRUE(aq.valid() != nullptr);
  ++aq.acquire_quiet_cols;
  HS_EXPECT_TRUE(aq.valid() == nullptr);
  Config id = test_config();
  id.beacon_interdigit_timeout_cols = id.beacon_span_cols() / 4 - 1;
  HS_EXPECT_TRUE(id.valid() != nullptr);

  // Stale-frame window order: tick()'s poll-path parser reset must be tighter
  // than BeaconParser::feed's interdigit test. Boundary is exclusive (equal
  // windows are rejected).
  Config so = test_config();
  so.beacon_interdigit_timeout_cols = so.acquire_quiet_cols;
  HS_EXPECT_TRUE(so.valid() != nullptr);
  ++so.beacon_interdigit_timeout_cols;
  HS_EXPECT_TRUE(so.valid() == nullptr);

  // Beacon tail quiet: a frame started at the beacon point must land its last
  // pulse and the receiver's quiet window inside the [W/4, W/2) half-window, or
  // the HALF burst is appended to the last digit burst. Boundary is exclusive —
  // the slack absorbs the sub-column offset of the scheduling tick.
  Config bq = test_config();
  bq.acquire_quiet_cols = bq.W / 4 - bq.beacon_span_cols();
  HS_EXPECT_EQ(bq.beacon_frame_cols(), bq.W / 4);
  HS_EXPECT_TRUE(bq.valid() != nullptr);
  --bq.acquire_quiet_cols;
  HS_EXPECT_TRUE(bq.valid() == nullptr);

  Config dg = test_config();
  dg.gate_cols = 7 * dg.beacon_pitch_cols + 1;
  HS_EXPECT_TRUE(dg.valid() != nullptr);
  --dg.gate_cols;
  HS_EXPECT_TRUE(dg.valid() == nullptr);

  // A glitch filter at or above a burst's pulse spacing, less the emitter's
  // lateness budget, drops every pulse after the first, so all three symbols
  // would decode as Symbol::HALF. The beacon pitch is the tighter of the two
  // bounds; the shipped filter clears both. Boundary is exclusive.
  Config gf = test_config();
  const uint32_t gf_bound = gf.beacon_pitch_cycles() - gf.late_censor_cycles();
  HS_EXPECT_TRUE(gf.glitch_filter_cycles < gf_bound);
  HS_EXPECT_TRUE(gf.glitch_filter_cycles <
                 gf.pulse_pitch_cycles() - gf.late_censor_cycles());
  gf.glitch_filter_cycles = gf_bound;
  HS_EXPECT_TRUE(gf.valid() != nullptr);
  gf.glitch_filter_cycles = gf_bound - 1;
  HS_EXPECT_TRUE(gf.valid() == nullptr);
}

/**
 * @brief Probes the remaining Config::valid() clauses at their boundaries.
 * @details Each config below violates exactly one clause, so a deleted clause
 * leaves the matching valid() call returning nullptr and turns its expectation
 * red. The clauses whose bound is a derived beacon quantity are probed in
 * test_helpers().
 */
inline void test_config_validation() {
  HS_EXPECT_TRUE(test_config().valid() == nullptr);

  // Odd W: boundary_column(HALF) and every arm-B half-image offset truncate
  // W/2. 289 keeps W/4 and cycles_per_column at their 288 values, so nothing
  // else moves.
  Config ow = test_config();
  ow.W = 289;
  HS_EXPECT_TRUE(ow.valid() != nullptr);

  Config gz = test_config();
  gz.gate_cols = 0;
  HS_EXPECT_TRUE(gz.valid() != nullptr);

  // gate_cols < W/4 is what lets snap()'s distance check subsume the boundary-
  // identity check. No config violates it alone: beacon_frame_cols() < W/4
  // forces W/4 above the 7*beacon_pitch_cols + 1 digit-gate bound, so that
  // clause (probed in test_helpers) is always the tighter of the two.
  const Config gw = test_config();
  HS_EXPECT_TRUE(7 * gw.beacon_pitch_cols + 1 < gw.W / 4);
  HS_EXPECT_TRUE(gw.gate_cols < gw.W / 4);

  Config rj = test_config();
  rj.reject_fallback = 0;
  HS_EXPECT_TRUE(rj.valid() != nullptr);

  Config gz2 = test_config();
  gz2.glitch_filter_cycles = 0;
  HS_EXPECT_TRUE(gz2.valid() != nullptr);

  Config pz = test_config();
  pz.pulse_pitch_cols = 0;
  HS_EXPECT_TRUE(pz.valid() != nullptr);

  // A gap that does not outlast the boundary-burst pitch terminates the burst
  // between its own pulses. Boundary is exclusive.
  Config gp = test_config();
  gp.gap_timeout_cols = gp.pulse_pitch_cols;
  HS_EXPECT_TRUE(gp.valid() != nullptr);
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
  bg.gap_timeout_cols = bg.beacon_pitch_cols;
  HS_EXPECT_TRUE(bg.valid() != nullptr);

  // Epoch indices are taken mod effect_count and ride a 6-bit beacon field.
  Config ec = test_config();
  ec.effect_count = 0;
  HS_EXPECT_TRUE(ec.valid() != nullptr);
  ec.effect_count = 65;
  HS_EXPECT_TRUE(ec.valid() != nullptr);
  ec.effect_count = 64;
  HS_EXPECT_TRUE(ec.valid() == nullptr);

  Config cz = test_config();
  cz.commit_revs = 0;
  HS_EXPECT_TRUE(cz.valid() != nullptr);

  // A negative epoch_repeats casts to a huge uint32_t and wraps the refractory
  // sum, so the relation below it reads true; only the explicit sign gate
  // rejects the config.
  Config ne = test_config();
  ne.epoch_repeats = -1;
  HS_EXPECT_TRUE(ne.refractory_revs >
                 ne.commit_revs + static_cast<uint32_t>(ne.epoch_repeats));
  HS_EXPECT_TRUE(ne.valid() != nullptr);
  ne.epoch_repeats = 0;
  HS_EXPECT_TRUE(ne.valid() == nullptr);

  // The EPOCH dedup window must outlast the whole redundancy train, or the
  // train's own last repeat re-arms the commit it just deduped. Boundary is
  // exclusive.
  Config rf = test_config();
  rf.refractory_revs = rf.commit_revs + static_cast<uint32_t>(rf.epoch_repeats);
  HS_EXPECT_TRUE(rf.valid() != nullptr);
  ++rf.refractory_revs;
  HS_EXPECT_TRUE(rf.valid() == nullptr);

  // An effect shorter than the dedup window would advance before the window
  // that protects its own commit closes. Boundary is exclusive.
  Config re = test_config();
  re.revs_per_effect = re.refractory_revs;
  HS_EXPECT_TRUE(re.valid() != nullptr);
  ++re.revs_per_effect;
  HS_EXPECT_TRUE(re.valid() == nullptr);

  uint32_t variable_revolutions[4] = {40, 48, 56, 64};
  Config vr = test_config();
  vr.effect_revolutions = variable_revolutions;
  HS_EXPECT_TRUE(vr.valid() == nullptr);
  variable_revolutions[2] = vr.refractory_revs;
  HS_EXPECT_TRUE(vr.valid() != nullptr);

  // A beacon cadence inside the construction window would land identity
  // traffic on the commit boundary. Boundary is exclusive.
  Config bc = test_config();
  bc.beacon_period_revs = bc.commit_revs;
  HS_EXPECT_TRUE(bc.valid() != nullptr);
  ++bc.beacon_period_revs;
  HS_EXPECT_TRUE(bc.valid() == nullptr);

  // The live-takeover grid must divide 64 so a beacon's mod-64 revolution
  // count lands on the same grid as the master's true count.
  Config jg = test_config();
  jg.join_grid_revs = 0;
  HS_EXPECT_TRUE(jg.valid() != nullptr);
  jg.join_grid_revs = 3;
  HS_EXPECT_TRUE(jg.valid() != nullptr);
  jg.join_grid_revs = 8;
  HS_EXPECT_TRUE(jg.valid() == nullptr);
}

/**
 * @brief Verifies symbol/boundary mapping: odd pulse counts classify to
 *        HALF/ZERO/ZERO_EPOCH, even or out-of-range counts are INVALID, and
 *        each symbol maps to its boundary column.
 */
inline void test_alphabet() {
  HS_EXPECT_TRUE(classify_count(1) == Symbol::HALF);
  HS_EXPECT_TRUE(classify_count(3) == Symbol::ZERO);
  HS_EXPECT_TRUE(classify_count(5) == Symbol::ZERO_EPOCH);
  // Even counts (single lost/spurious edge) and out-of-range are INVALID.
  for (uint32_t n : {0u, 2u, 4u, 6u, 7u, 8u, 9u, 255u})
    HS_EXPECT_TRUE(classify_count(n) == Symbol::INVALID);
  HS_EXPECT_TRUE(symbol_boundary(Symbol::HALF) == Boundary::HALF);
  HS_EXPECT_TRUE(symbol_boundary(Symbol::ZERO) == Boundary::ZERO);
  HS_EXPECT_TRUE(symbol_boundary(Symbol::ZERO_EPOCH) == Boundary::ZERO);
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
 * @brief Verifies edge mailbox burst accumulation and the 100 µs glitch
 *        filter: sub-window spikes are rejected without resetting the filter
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
  // EMI spike < 100 µs after an accepted edge is rejected…
  m.on_edge(1000 + 4 * COL + GLITCH / 2, GLITCH);
  // …and does not reset the filter reference.
  m.on_edge(1000 + 6 * COL, GLITCH);
  HS_EXPECT_FALSE(burst_complete(m, 1000 + 7 * COL, 4 * COL));
  HS_EXPECT_TRUE(burst_complete(m, 1000 + 10 * COL + 1, 4 * COL));
  const BurstSnapshot s = claim(m);
  HS_EXPECT_EQ(s.count, 4u);
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
 *          accepted edge, so the terminating gap never opens. Without the
 *          duration term the mailbox stays jammed for as long as the noise
 *          lasts: nothing claimed, nothing decoded, no telemetry counter moving.
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
  HS_EXPECT_TRUE(classify_count(out.count) == Symbol::INVALID);

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
 * @details age_prior() (called every column by the flywheel poll) retires the
 *          prior once the wire is quiet past the filter window, so a real edge
 *          after wrap is never falsely rejected by a pseudo-random modular
 *          difference.
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
    // The flywheel keeps polling during the silence and retires the stale
    // reference within a column (COL > GLITCH), long before the counter wraps.
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
 * @details The flywheel samples `now` some cycles before its IRQ-off bracket
 *          opens; a sync edge landing in that window leaves the mailbox
 *          timestamps ahead of `now`, and the unsigned difference underflows to
 *          ~2³². Claiming there would take a burst still in flight and
 *          misclassify the symbol; retiring the prior there would drop the
 *          glitch reference the very edge just set.
 */
inline void test_mailbox_rejects_backward_clock() {
  const uint32_t GLITCH = 60000u;
  const uint32_t now = 100000u;
  const uint32_t skew = 40u; // edge accepted after `now` was sampled

  EdgeMailbox m;
  BurstSnapshot out{};
  m.on_edge(now + skew, GLITCH);
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
 * @details The stale burst would otherwise feed ACQUIRE's unconditional
 *          hard-snap. The emitter is reset by the same code path.
 */
inline void test_seed_clears_mailbox() {
  const Config cfg = test_config();
  const uint32_t col = cfg.cycles_per_column();
  SyncBoard board(cfg);
  board.seed(1000u, /*is_master=*/false);
  board.on_sync_edge(2000u);
  board.on_sync_edge(2000u + 4 * col);
  HS_EXPECT_TRUE(burst_complete(board.mailbox(), 2000u + 100 * col,
                                cfg.gap_timeout_cycles()));
  board.seed(3000u, false);
  HS_EXPECT_FALSE(burst_complete(board.mailbox(), 3000u + 100 * col,
                                 cfg.gap_timeout_cycles()));
  HS_EXPECT_EQ(claim(board.mailbox()).count, 0u);
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

  const Config replacement = test_config(3);
  board.configure(replacement);
  HS_EXPECT_EQ(board.config().effect_count, 3);
  HS_EXPECT_EQ(board.build_word(), 0u);

  board.seed(3000u, true);
  HS_EXPECT_EQ(SyncBoard::build_gen_of(board.build_word()), 1u);
}

/**
 * @brief Verifies a tick folding several boundaries reports the FINAL one:
 *        TickActions::zero_crossing names the boundary that opened the display
 *        window now in effect, which the driver publishes as the window half.
 * @details Also pins the §5.1 fold bound: the gate runs once per crossing so the
 *          flip counter and the coast telemetry see all N, while `flip` stays a
 *          single bool — N windows consume one advance_display().
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
  HS_EXPECT_TRUE(board.flywheel().current_boundary() == Boundary::HALF);
  HS_EXPECT_EQ(board.telemetry_snapshot().flips, 3u);
  HS_EXPECT_EQ(board.telemetry_snapshot().max_coast_halves, 3u);

  // The next boundary is a ZERO, and it is reported.
  const TickActions z =
      board.tick(1000u + 4u * cfg.cycles_per_half_rev + 10u, nullptr);
  HS_EXPECT_TRUE(z.flip);
  HS_EXPECT_TRUE(z.zero_crossing);
  HS_EXPECT_TRUE(board.flywheel().current_boundary() == Boundary::ZERO);
  HS_EXPECT_EQ(board.telemetry_snapshot().flips, 4u);

  // An even fold returns to the half it started from: HALF then ZERO reports
  // the ZERO, so the published window half still names the open window.
  const TickActions e =
      board.tick(1000u + 6u * cfg.cycles_per_half_rev + 10u, nullptr);
  HS_EXPECT_TRUE(e.flip);
  HS_EXPECT_TRUE(e.zero_crossing);
  HS_EXPECT_TRUE(board.flywheel().current_boundary() == Boundary::ZERO);
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
  HS_EXPECT_EQ(d[4], (1 * 3 + 2 * 3 + 3 * 5 + 4 * 5) & 7);

  /**
   * @brief Feeds all five digit bursts through a fresh parser, optionally
   *        corrupting digit 1 or the checksum.
   * @param digits The five encoded beacon digits.
   * @param corrupt_digit If true, flip digit 1 before feeding.
   * @param corrupt_checksum If true, flip the checksum digit before feeding.
   * @param out Receives the decoded frame on success.
   * @return True iff the frame decoded with no rejection.
   */
  auto feed_frame = [&cfg](const uint8_t digits[5], bool corrupt_digit,
                           bool corrupt_checksum, BeaconFrame *out) {
    BeaconParser p;
    bool got = false, rejected = false;
    uint32_t t = 1000;
    for (int i = 0; i < 5; ++i) {
      uint8_t v = digits[i];
      if (corrupt_digit && i == 1)
        v = static_cast<uint8_t>((v + 1) & 7);
      if (corrupt_checksum && i == 4)
        v = static_cast<uint8_t>((v + 1) & 7);
      BurstSnapshot s{static_cast<uint32_t>(v) + 1u, t, t + v * COL};
      bool r = false;
      got = p.feed(s, cfg, out, &r);
      rejected = rejected || r;
      t += 12 * COL;
    }
    return got && !rejected;
  };

  BeaconFrame f{};
  for (int idx : {0, 1, 27, 63}) {
    for (uint32_t rev : {0u, 1u, 39u, 63u}) {
      encode_beacon_digits(idx, rev, d);
      HS_EXPECT_TRUE(feed_frame(d, false, false, &f));
      HS_EXPECT_EQ(f.effect_index, idx);
      HS_EXPECT_EQ(f.rev_count, rev);
    }
  }
  encode_beacon_digits(27, 45, d);
  HS_EXPECT_FALSE(feed_frame(d, true, false, &f));
  HS_EXPECT_FALSE(feed_frame(d, false, true, &f));

  // The position-weighted checksum rejects two corruption classes a plain
  // digit-sum is blind to (both preserve the sum): a transposition of two
  // distinct digits, and a compensating ±1 pair (a pulse miscounted from one
  // burst into the next). d = {3,3,5,5,chk}.
  encode_beacon_digits(27, 45, d);
  const uint8_t transposed[5] = {d[2], d[1], d[0], d[3], d[4]}; // swap d0,d2
  HS_EXPECT_FALSE(feed_frame(transposed, false, false, &f));
  const uint8_t compensated[5] = {static_cast<uint8_t>(d[0] + 1),
                                  static_cast<uint8_t>(d[1] - 1), d[2], d[3],
                                  d[4]};
  HS_EXPECT_FALSE(feed_frame(compensated, false, false, &f));

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
 * @details The parser's own staleness test is a modular difference, so a partial
 *          frame left standing outlives a cycle-counter wrap and reads a fresh
 *          train's first digit as in-window; the tick-driven reset fires within
 *          columns of the truncated train, long before a wrap can accumulate.
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
  board.tick(head + 40 * col, nullptr); // quiet past the ACQUIRE guard

  // A complete frame from column 90 on, spaced exactly as schedule_beacon does.
  uint8_t d[5];
  encode_beacon_digits(2, 3, d);
  feed_beacon_train(board, cfg, col, 1000u + 90u * col, d);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_ok, 1u);
  // The aged-out partial is a dropped frame and counts as one.
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_rejected, 1u);
  HS_EXPECT_EQ(board.content().effect_index, 2);
}

/**
 * @brief Verifies the §6.3.4 confirmation rule: a live board changes effect
 *        index only after two consecutive beacons name the same one.
 * @details The position-weighted checksum provably catches any single
 *          mis-counted digit, but not a *shift*: a stray burst in the quiet
 *          ahead of digit 0 pushes the frame along one place, and exactly one
 *          of the eight intruder values re-satisfies the weighted sum (p =
 *          1/8), so the tuple decodes valid and wrong. Applied unconfirmed it
 *          would tear down a healthy effect and display the wrong one — the
 *          fail-wrong window §6.3.3 rules out. Frames here are fed straight at
 *          the board one per half-revolution, in the mid-half quiet where the
 *          §6.4 demarcation routes bursts to the beacon parser.
 */
inline void test_beacon_shift_needs_confirmation() {
  // Full 64-roster: every 6-bit index is in range, so the shifted frame gets
  // past the out-of-range rejection and the confirmation rule alone is on
  // trial.
  const Config cfg = test_config(64);
  const uint32_t col = cfg.cycles_per_column();
  SyncBoard board(cfg);
  board.seed(1000u, /*is_master=*/false);
  flywheel_mut(board).force_lock();
  content_mut(board).identity_known = true;
  content_mut(board).effect_index = 1;
  content_mut(board).rev_in_effect = 5;

  // Half-rev k's frame starts at column 40: far from both boundaries, and the
  // whole frame lands inside the half-rev. Fold the intervening crossings
  // before the digits go in, so a frame's rev digits match the board's count.
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
    encode_beacon_digits(index, board.content().rev_in_effect, d);
    feed_digits(start, d);
  };

  // The shift: [emi, d0, d1, d2, d3] with the real d3 read as the checksum.
  const uint32_t shift_start = open_frame();
  uint8_t truth[5];
  encode_beacon_digits(1, board.content().rev_in_effect, truth);
  int passing = 0;
  uint8_t shifted[5] = {};
  for (uint8_t emi = 0; emi < 8; ++emi) {
    const uint8_t s[5] = {emi, truth[0], truth[1], truth[2], truth[3]};
    if (((1u * s[0] + 2u * s[1] + 3u * s[2] + 4u * s[3]) & 7u) != s[4])
      continue;
    ++passing;
    for (int i = 0; i < 5; ++i)
      shifted[i] = s[i];
  }
  HS_EXPECT_EQ(passing, 1); // p = 1/8, exactly one intruder value
  const int32_t shifted_index = shifted[0] * 8 + shifted[1];
  HS_EXPECT_TRUE(shifted_index != board.content().effect_index);
  HS_EXPECT_TRUE(shifted_index < cfg.effect_count);

  // A lone shifted frame decodes but must not change what is displayed.
  feed_digits(shift_start, shifted);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_ok, 1u);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_rejected, 0u);
  HS_EXPECT_EQ(board.content().effect_index, 1);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacon_index_corrections, 0u);
  HS_EXPECT_EQ(board.build_word(), 0u); // no rebuild published

  // A good frame agreeing with the displayed index clears the candidate.
  feed_frame(1);
  HS_EXPECT_EQ(board.content().effect_index, 1);

  // Two frames naming *different* indices confirm nothing either.
  feed_frame(3);
  HS_EXPECT_EQ(board.content().effect_index, 1);
  feed_frame(2);
  HS_EXPECT_EQ(board.content().effect_index, 1);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacon_index_corrections, 0u);

  // The second agreeing frame applies: index adopted, rebuild published.
  feed_frame(2);
  HS_EXPECT_EQ(board.content().effect_index, 2);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacon_index_corrections, 1u);
  HS_EXPECT_EQ(SyncBoard::build_index_of(board.build_word()), 2);
  HS_EXPECT_EQ(SyncBoard::build_gen_of(board.build_word()), 1u);
  // The whole scenario rode the beacon path: nothing reached the snap gate.
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_ok, 5u);
  HS_EXPECT_EQ(board.telemetry_snapshot().symbols_accepted, 0u);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacon_rev_mismatches, 0u);
  HS_EXPECT_TRUE(board.lock() == LockState::LOCKED);
}

/**
 * @brief Verifies a checksum-valid beacon naming an index past the roster is
 *        dropped whole (§6.4 integrity by rejection).
 * @details A 6-bit index leaves 42 unreachable values on this 4-effect roster
 *          — positive evidence of corruption the 3-bit checksum passes with
 *          p = 1/8. Folding it mod effect_count would join a rebooting board
 *          to a wrong-but-valid effect, the fail-wrong outcome §6.3.3 rules
 *          out; the frame must count as rejected, like any corrupt frame.
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
  HS_EXPECT_FALSE(board.content().identity_known);
  HS_EXPECT_EQ(board.build_word(), 0u); // no rebuild published

  // The next in-range beacon joins normally.
  feed_frame(2, 1000u + cfg.cycles_per_half_rev + 40u * col);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_ok, 1u);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_rejected, 1u);
  HS_EXPECT_TRUE(board.content().identity_known);
  HS_EXPECT_EQ(board.content().effect_index, 2);
  HS_EXPECT_EQ(SyncBoard::build_index_of(board.build_word()), 2);
}

/**
 * @brief Verifies the §6.4 rev cross-check fold (beacon_rev_resync_delta)
 *        resolves the 63↔0 mod-64 seam.
 * @details Production effects span 960 revs and the sim configs span 40, so
 *          rev_in_effect routinely exceeds 63 and its 6-bit residue wraps —
 *          but test_sim_rev_resync only ever slips a board by +2 at rev 5, far
 *          from the wrap, leaving the fold's `+96 %64 -32` seam arithmetic
 *          unexercised. Drive it directly: the fold maps a beacon residue and
 *          the board's current rev_in_effect to the smallest signed slip in
 *          [-32, 31], so applying it as `rev_in_effect + delta` restores the
 *          exact residue across the wrap, in either direction.
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

  // The applied fold lands on the beacon's residue across the seam, exactly as
  // handle_beacon_burst uses it (rev_in_effect + delta, then re-read mod 64).
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
 *        preserved across thousands of 32-bit counter wraps via the rebase
 *        rule.
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
    return static_cast<int64_t>(
        floor_div(delta * 144, static_cast<int64_t>(period)));
  };

  // Position over one half-rev, at nominal and trim-extreme periods
  // (±40 ppm ≈ ±1500 cycles): zero truncation drift vs the reference.
  for (int32_t trim : {0, +1500, -1500}) {
    const uint32_t period = PERIOD + trim;
    Flywheel f(cfg);
    f.set_cycles_per_half_rev(period);
    f.seed(1000000u);
    bool ok = true;
    for (int64_t delta = 0; delta < period; delta += 12347) {
      const int32_t want = static_cast<int32_t>(ref_cols(delta, period) % 288);
      if (f.position(1000000u + static_cast<uint32_t>(delta)) != want) {
        ok = false;
        break;
      }
    }
    HS_EXPECT_TRUE(ok);
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
    HS_EXPECT_TRUE(c1.boundary == Boundary::HALF);
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
  // folds across several wraps; every crossing lands on its boundary column
  // and the epoch stays an exact integer multiple ahead.
  {
    Flywheel f(cfg);
    uint32_t t = 0xFFFFFFFFu - PERIOD / 3; // wrap almost immediately
    f.seed(t);
    bool ok = true;
    for (int k = 1; k <= 5000; ++k) { // ~5.2 minutes of mock time, 43 wraps
      t += PERIOD;
      const Crossing c = f.fold(t);
      if (!c.crossed || c.at_cycles != t || f.fold(t).crossed ||
          f.position(t) != boundary_column(c.boundary, 288)) {
        ok = false;
        break;
      }
    }
    HS_EXPECT_TRUE(ok);
  }
}

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
    HS_EXPECT_TRUE(f.snap(Boundary::HALF, 1000000u + PERIOD - 2 * COL, &err) ==
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
    HS_EXPECT_TRUE(last == Flywheel::SnapOutcome::REJECTED_FELL_BACK);
    HS_EXPECT_TRUE(f.lock() == LockState::ACQUIRE);
    // ACQUIRE: hard snap, relocks.
    t += 10 * COL;
    HS_EXPECT_TRUE(f.snap(Boundary::ZERO, t, &err) ==
                   Flywheel::SnapOutcome::ACCEPTED);
    HS_EXPECT_TRUE(f.lock() == LockState::LOCKED);
    HS_EXPECT_EQ(f.position(t), 0);
  }

  // Corrupted timebase: a board whose epoch is garbage rejects good symbols
  // but re-acquires via the fallback within R symbols (spec §12).
  {
    Flywheel f(cfg);
    f.seed(1000000u);
    f.force_lock();
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
    HS_EXPECT_TRUE(f.lock() == LockState::LOCKED);
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
  HS_EXPECT_TRUE(board.lock() == LockState::ACQUIRE);
  board.tick(head + 40 * col, nullptr); // past the interdigit window
  HS_EXPECT_EQ(board.telemetry_snapshot().symbols_rejected_gate, rejected);
}

/**
 * @brief Verifies the §5.3 quiet-before guard: an ACQUIRE board hard-snaps only
 *        on a burst preceded by t_QB of wire silence, so a beacon digit train
 *        cannot capture a just-rebooted board mid-frame.
 * @details The head of the beacon for effect index 8 is a 2-pulse burst — an
 *          even count, no symbol — and its second digit is a single pulse, a
 *          valid HALF count 5 columns behind it. Without the guard that digit
 *          is a hard snap to a phase 5 columns off a beacon's mid-frame
 *          position. The same burst after t_QB of silence is the positive
 *          control: it snaps, so the assertion below is the gap, not the burst.
 */
inline void test_acquire_quiet_before_guard() {
  const Config cfg = test_config();
  const uint32_t col = cfg.cycles_per_column();
  const uint32_t head = 1000u + 40u * col;
  const BurstSnapshot digit0{2, head, head + col};

  SyncBoard board(cfg);
  board.seed(1000u, /*is_master=*/false);
  HS_EXPECT_TRUE(board.lock() == LockState::ACQUIRE);
  board.tick(head + col + static_cast<uint32_t>(cfg.gap_timeout_cols) * col,
             &digit0);
  // The next digit, spaced exactly as schedule_beacon spaces them.
  const uint32_t d1 =
      head + col + static_cast<uint32_t>(cfg.gap_timeout_cols + 1) * col;
  const BurstSnapshot digit1{1, d1, d1};
  board.tick(d1 + static_cast<uint32_t>(cfg.gap_timeout_cols) * col, &digit1);
  HS_EXPECT_TRUE(board.lock() == LockState::ACQUIRE);
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
  HS_EXPECT_TRUE(quiet.lock() == LockState::LOCKED);
  HS_EXPECT_EQ(quiet.telemetry_snapshot().symbols_accepted, 1u);
  HS_EXPECT_EQ(quiet.flywheel().position(iso), cfg.W / 2);
}

/**
 * @brief Verifies a board still in ACQUIRE decodes a whole beacon train and
 *        adopts (effect, rev) from it (§6.4), instead of waiting for lock.
 * @details A train's first digit is preceded by the same wire silence a
 *          boundary symbol is, so the quiet-before guard routes it to the
 *          symbol path; withholding it from the parser would leave every frame
 *          begun in ACQUIRE one digit short forever. The head of the beacon for
 *          effect 9 is a 2-pulse burst — an even count, no symbol — so the
 *          board stays in ACQUIRE across the train, and the frame still
 *          completes. The isolated burst starts a fresh frame, so it can never
 *          both complete a frame and hard-snap: the closing symbol below snaps
 *          on a clean timebase with the beacon identity intact.
 */
inline void test_acquire_beacon_train_joins() {
  const Config cfg = test_config(16);
  const uint32_t col = cfg.cycles_per_column();
  SyncBoard board(cfg);
  board.seed(1000u, /*is_master=*/false);
  HS_EXPECT_TRUE(board.lock() == LockState::ACQUIRE);
  HS_EXPECT_FALSE(board.content().identity_known);

  // A full train from column 40 on, spaced exactly as schedule_beacon does.
  uint8_t d[5];
  encode_beacon_digits(9, 5, d);
  HS_EXPECT_EQ(static_cast<int>(d[0]) + 1, 2); // head: no valid symbol count
  uint32_t f = feed_beacon_train(board, cfg, col, 1000u + 40u * col, d);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_ok, 1u);
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_rejected, 0u);
  HS_EXPECT_TRUE(board.content().identity_known);
  HS_EXPECT_EQ(board.content().effect_index, 9);
  HS_EXPECT_EQ(board.content().rev_in_effect, 5u);
  // The head reached the symbol path too, and was discarded there on its count.
  HS_EXPECT_EQ(board.telemetry_snapshot().symbols_discarded_invalid, 1u);
  HS_EXPECT_EQ(board.telemetry_snapshot().symbols_accepted, 0u);
  HS_EXPECT_TRUE(board.lock() == LockState::ACQUIRE);

  // The boundary path is untouched: the next isolated symbol still hard-snaps.
  const uint32_t z = f + cfg.acquire_quiet_cycles();
  const uint32_t zspan = 2u * static_cast<uint32_t>(cfg.pulse_pitch_cols) * col;
  const BurstSnapshot zero{3, z, z + zspan};
  board.tick(z + zspan + static_cast<uint32_t>(cfg.gap_timeout_cols) * col,
             &zero);
  HS_EXPECT_TRUE(board.lock() == LockState::LOCKED);
  HS_EXPECT_EQ(board.telemetry_snapshot().symbols_accepted, 1u);
  HS_EXPECT_EQ(board.flywheel().position(z), 0);
  HS_EXPECT_EQ(board.content().effect_index, 9);
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

  // Boundary scheduled in the FUTURE (now before at_cycles): not late — it must
  // be accepted and then emitted once `now` reaches the boundary. An unsigned
  // `now - at_cycles` wraps a future boundary to a huge positive lateness.
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

  // Masked mid-burst past the budget: remaining pulses are aborted, the
  // truncated count degrades to a missed/invalid symbol downstream.
  {
    SymbolEmitter e;
    const uint32_t b = 1000000u;
    HS_EXPECT_TRUE(e.schedule_boundary(Symbol::ZERO_EPOCH, b, b, cfg));
    HS_EXPECT_TRUE(e.tick(b, cfg, &aborted)); // pulse 1 on time
    // Next due at b+2col; first wake after the mask is way late.
    HS_EXPECT_FALSE(
        e.tick(b + 2 * COL + cfg.late_censor_cycles() + 1, cfg, &aborted));
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
    HS_EXPECT_TRUE(e.drop_pending_emission() == Dropped::NONE); // idle
    uint8_t d[5];
    encode_beacon_digits(3, 9, d);
    HS_EXPECT_TRUE(e.schedule_beacon(d, 2000000u, cfg));
    HS_EXPECT_FALSE(e.schedule_beacon(d, 2000000u, cfg));
    HS_EXPECT_FALSE(e.idle());
    HS_EXPECT_TRUE(e.drop_pending_emission() == Dropped::BEACON);
    HS_EXPECT_TRUE(e.idle());
    HS_EXPECT_TRUE(e.schedule_boundary(Symbol::HALF, 2000000u, 2000000u, cfg));

    // An undrained boundary symbol is not a beacon drop.
    HS_EXPECT_TRUE(e.drop_pending_emission() == Dropped::BOUNDARY);
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
    HS_EXPECT_TRUE(e.drop_pending_emission() == Dropped::BOUNDARY);
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
  HS_EXPECT_FALSE(SyncBoardTestAccess::beacon_done(board));
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_busy_dropped, 1u);

  SyncBoardTestAccess::maybe_schedule_beacon(board, beacon_at);
  HS_EXPECT_FALSE(SyncBoardTestAccess::beacon_done(board));
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_busy_dropped, 1u);

  emitter.drop_pending_emission();
  SyncBoardTestAccess::maybe_schedule_beacon(board, beacon_at);
  HS_EXPECT_TRUE(SyncBoardTestAccess::beacon_done(board));
  HS_EXPECT_EQ(board.telemetry_snapshot().beacons_busy_dropped, 1u);
}

// ── Master beacon late-coast bound (§6.4) ───────────────────────────────────

/**
 * @brief Verifies a masked-ISR coast that reaches the beacon point late does
 *        not queue a beacon whose tail overruns HALF and trips the emitter's
 *        wire-busy trap on the on-time HALF symbol.
 * @details Drives a master across the beacon-due revolution, resuming its first
 *          post-ZERO tick at a chosen column to model the coast. The bound is
 *          sized from the payload actually encoded, so the window spans many
 *          columns: the last admissible start emits fully before HALF, one
 *          column later is censored — no pulses in the beacon window and no
 *          trap at HALF. A start part-way into that last column is censored
 *          too, since the fit charges the emitter's lateness budget as well.
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
    // emitter is idle before the coast, as it would be on real hardware.
    m.tick(epoch1 + 2u * COL, nullptr);
    m.tick(epoch1 + 4u * COL, nullptr);

    int pulses = 0;
    for (int32_t c = resume_col; c <= 150; ++c) {
      // Round the column instant up: cycles_per_column() truncates, so a whole
      // multiple of it lands a column short past ~column 100.
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
  // Sizing the bound from the payload rather than the all-sevens worst case is
  // what keeps the window more than one column wide.
  HS_EXPECT_GT(last_start, cfg.W / 4);

  // On-time at the beacon point: the frame schedules and emits fully (proves
  // the beacon is due so the late-start contrast below is meaningful).
  HS_EXPECT_GT(run(cfg.W / 4), 0);
  // The last admissible start still emits.
  HS_EXPECT_GT(run(last_start), 0);
  HS_EXPECT_EQ(late_dropped, 0u);
  // Late start past the safe bound: censored — no beacon pulses, and the HALF
  // crossing at column 144 schedules without tripping the wire-busy trap (a
  // trap would __builtin_trap the whole suite).
  HS_EXPECT_EQ(run(last_start + 1), 0);
  // The skip is counted once for the revolution, not once per late tick.
  HS_EXPECT_EQ(late_dropped, 1u);
  // Resuming part-way through the last admissible column is late too: the frame
  // is anchored on the tick, not on the column it falls in, and its last pulse
  // may go out up to the emitter's ½-column lateness budget after its due time —
  // together enough to leave the receiver less than its quiet window before HALF.
  HS_EXPECT_EQ(run(last_start, COL / 2 + COL / 8), 0);
  HS_EXPECT_EQ(late_dropped, 1u);
}

/**
 * @brief Verifies a master coast past 2^31 cycles is counted and recovered
 *        rather than wedging the flywheel silently.
 * @details fold() reads (now - epoch_cycles) as int32 and reports no crossing on
 *          a negative one, so without the re-anchor the master would never
 *          cross another boundary and would emit no further boundary symbols —
 *          with no counter moving to say so.
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
  // fold() alone can never catch up again.
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
 *          next crossing; a stale dedup state would drop that flip silently.
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
 *          frame's tail to clear acquire_quiet_cycles before the HALF burst. A
 *          closer tail is folded into the last digit burst by the receiver's
 *          gap timeout, consuming the boundary symbol instead of decoding it.
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
 *          full column late there, censoring the primary copy (§5.2). A train
 *          that spent repeats only on symbols reaching the wire would slide a
 *          copy to B+R+1, where a receiver's j = rev_in_effect − revs_per_effect
 *          lands outside the window, falls back to j = 0, and commits R
 *          revolutions late; a fully censored train would emit past the commit
 *          into the new effect. Counts pulses in each boundary's burst window:
 *          5 = ZERO_EPOCH, 3 = plain ZERO.
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

// ── Multi-board simulator ───────────────────────────────────────────────────

// Stand-in for the constructed effect: the handoff tracks ownership by address
// only, so one instance per board exercises every path.
struct SimEffect {};

/**
 * @brief One simulated board: its SyncBoard engine plus the host-side state the
 *        simulator models around it.
 * @details The modeled state covers crystal offset/phase, the masked-IRQ latch,
 *          symbol-drop and EMI windows, the foreground build/commit model, and
 *          probes.
 */
struct SimBoard {
  SyncBoard board;
  /** The real device handoff, driven through apply_wake() each wake. */
  pov::EffectHandoff<SimEffect> handoff;
  SimEffect instance; /**< Address the foreground publishes when built. */
  bool master = false;
  int32_t ppm = 0;      /**< Crystal offset, parts per million. */
  uint64_t phase0 = 0;  /**< Local cycle-counter offset at g = 0. */
  double next_tick = 0; /**< Next flywheel wake, in global cycles. */
  double tick_step = 0; /**< Wake period in global cycles. */
  /** Flywheel phase this board was seeded at, columns ahead of the master. */
  int32_t birth_cols = 0;
  /** Masked-IRQ model: while g < mask_until, wakes coalesce and edges latch
      into a single delayed delivery (one latched flag per pin). */
  uint64_t mask_until = 0;
  bool edge_latched = false;
  std::vector<std::pair<uint64_t, uint64_t>> masks; /**< [from, to) global. */
  uint64_t drop_from = 0, drop_to = 0;              /**< Symbol drop window. */
  // Foreground model.
  uint32_t seen_gen = 0;
  int32_t pending_index = -1;
  uint64_t build_seed = 0; /**< hs::epoch_seed the foreground applied at the
                                last build pickup (device reseed mirror). */
  uint32_t pending_gen = 0;
  uint64_t pending_ready_g = 0;
  bool have_pending = false;
  uint64_t init_delay = 1000000; /**< ~1.7 ms construction time. */
  bool live = false;
  int32_t live_index = -1;
  uint64_t t = 0; /**< Frames shown: flips since this effect went live. */
  uint64_t swap_g = 0;
  bool trapped = false; /**< Commit deadline missed (device would HS_CHECK). */
  // Probes.
  uint64_t flips = 0;
  bool dark_now = true;

  /**
   * @brief Constructs a board wrapping a SyncBoard engine for config @p c.
   * @param c Sync configuration passed to the embedded SyncBoard.
   */
  explicit SimBoard(const Config &c) : board(c) {}
};

/**
 * @brief Event-driven multi-board simulator.
 * @details Advances the earliest-due flywheel wake across all boards, routes
 *          master pulses and injected EMI to downstream edge ISRs through the
 *          masked-IRQ model, and drives each board's foreground state.
 */
class Sim {
public:
  Config cfg;
  std::deque<SimBoard> boards;
  uint64_t g = 0;                            /**< Global time, cycles. */
  std::vector<std::pair<uint64_t, int>> emi; /**< (g, target), sorted. */
  size_t emi_pos = 0;

  /**
   * @brief Builds @p n boards (index 0 is the master) and seeds each flywheel.
   * @param c Shared sync configuration.
   * @param n Number of boards to construct.
   * @param ppm Per-board crystal offsets, parts per million (length @p n).
   * @param phase0 Common starting clock offset, cycles, at global time 0.
   * @details Each flywheel polls on a ⅛-column grid scaled by its ppm, with a
   *          small per-board boot stagger. Boards power on at unrelated rotor
   *          angles, so board i is born believing ZERO happened i·W/n columns
   *          ago — only acquisition closes that gap. The master defines phase
   *          and is born at offset 0.
   */
  Sim(const Config &c, int n, const int32_t *ppm, uint64_t phase0 = 0)
      : cfg(c) {
    const double step0 = double(c.cycles_per_half_rev) / (c.W / 2) / 8.0;
    for (int i = 0; i < n; ++i) {
      boards.emplace_back(c);
      SimBoard &b = boards.back();
      b.master = (i == 0);
      b.ppm = ppm[i];
      b.phase0 = phase0;
      b.birth_cols = i * c.W / n;
      b.tick_step = step0 * 1e6 / (1e6 + ppm[i]);
      b.next_tick = double(7 * (i + 1)); // small boot stagger
      b.board.seed(local_now(b, 0) - static_cast<uint32_t>(b.birth_cols) *
                                         c.cycles_per_column(),
                   b.master);
    }
  }

  /**
   * @brief Computes board @p b's local 32-bit cycle counter at global cycle
   *        @p gg: phase offset plus crystal skew (ppm).
   * @param b The board whose local clock is evaluated.
   * @param gg Global time, cycles.
   * @return The board-local cycle count, truncated to 32 bits to model CYCCNT
   *         wrap.
   */
  static uint32_t local_now(const SimBoard &b, uint64_t gg) {
    const int64_t skew = static_cast<int64_t>(gg) * b.ppm / 1000000;
    return static_cast<uint32_t>(gg + b.phase0 + static_cast<uint64_t>(skew));
  }

  /**
   * @brief Tests whether global cycle @p at falls inside an IRQ-mask window.
   * @param b The board whose mask windows are checked.
   * @param at Global time, cycles.
   * @return The window's end cycle (when delivery resumes) if masked; 0 if
   *         unmasked.
   */
  uint64_t masked_until(const SimBoard &b, uint64_t at) const {
    for (const auto &w : b.masks)
      if (at >= w.first && at < w.second)
        return w.second;
    return 0;
  }

  /**
   * @brief Routes a sync-wire edge at global cycle @p at to downstream board
   *        @p target.
   * @param target Index of the downstream board receiving the edge.
   * @param at Global time of the edge, cycles.
   * @details Dropped if in the board's deafen window, latched/merged if masked,
   *          otherwise fed to its edge ISR at the board-local timestamp.
   */
  void deliver_edge(int target, uint64_t at) {
    SimBoard &b = boards[target];
    if (b.master)
      return; // master's edge ISR is not attached
    if (at >= b.drop_from && at < b.drop_to)
      return;
    if (const uint64_t end = masked_until(b, at)) {
      b.edge_latched = true; // single latched flag: merged, delayed
      b.mask_until = end;
      return;
    }
    b.board.on_sync_edge(local_now(b, at));
  }

  /**
   * @brief Computes board @p b's next wake in global cycles.
   * @param b The board whose next wake is evaluated.
   * @return The scheduled wake, pushed to the mask end if its slot lands inside
   *         a mask (coalesced wakes).
   */
  double effective_tick(const SimBoard &b) const {
    const uint64_t m = masked_until(b, static_cast<uint64_t>(b.next_tick));
    return m ? double(m) : b.next_tick;
  }

  /**
   * @brief Advances global time to the earliest-due board wake, delivering any
   *        EMI edges before it, then runs that board's tick.
   */
  void step() {
    int bi = 0;
    double best = effective_tick(boards[0]);
    for (size_t i = 1; i < boards.size(); ++i) {
      const double e = effective_tick(boards[i]);
      if (e < best) {
        best = e;
        bi = static_cast<int>(i);
      }
    }
    const uint64_t tg = static_cast<uint64_t>(best + 0.5);
    // Deliver EMI edges due before this wake (edge ISRs run on their own).
    while (emi_pos < emi.size() && emi[emi_pos].first <= tg) {
      deliver_edge(emi[emi_pos].second, emi[emi_pos].first);
      ++emi_pos;
    }
    run_tick(boards[bi], tg);
    g = tg;
  }

  /**
   * @brief Steps the simulation until @p revs revolutions of global time have
   *        elapsed.
   * @param revs Number of revolutions to advance.
   */
  void run_revs(double revs) {
    const uint64_t until =
        g + static_cast<uint64_t>(revs * 2 * cfg.cycles_per_half_rev);
    while (g < until)
      step();
  }

  /**
   * @brief Steps the simulation until @p pred holds or @p max_revs elapse.
   * @tparam Pred Callable taking `Sim &` and returning bool.
   * @param pred Predicate evaluated after each step.
   * @param max_revs Maximum revolutions to advance before giving up.
   * @return True if @p pred returned true within the budget, false on timeout.
   */
  template <typename Pred>
  [[nodiscard]] bool run_until(Pred pred, double max_revs) {
    const uint64_t until =
        g + static_cast<uint64_t>(max_revs * 2 * cfg.cycles_per_half_rev);
    while (g < until) {
      step();
      if (pred(*this))
        return true;
    }
    return false;
  }

  /**
   * @brief Computes board @p i's current flywheel column position.
   * @param i Index of the board to sample.
   * @return The column position, evaluated at the board's own local clock at
   *         the present global time.
   */
  int32_t board_pos(int i) const {
    return boards[i].board.flywheel().position(local_now(boards[i], g));
  }

  /**
   * @brief Computes the largest circular column distance of any locked board
   *        from the master.
   * @return The worst-case phase error in columns; 0 if no downstream board is
   *         locked.
   */
  int32_t max_phase_err() const {
    const int32_t mp = board_pos(0);
    int32_t worst = 0;
    for (size_t i = 1; i < boards.size(); ++i) {
      if (boards[i].board.lock() != LockState::LOCKED)
        continue;
      const int32_t d = circ_dist(board_pos(static_cast<int>(i)), mp, cfg.W);
      if (d > worst)
        worst = d;
    }
    return worst;
  }

private:
  /**
   * @brief Runs one flywheel wake for board @p b at global cycle @p tg.
   * @param b The board being woken.
   * @param tg Global time of this wake, cycles.
   * @details Advances the poll grid, delivers any latched edge, drives
   *          board.tick(), fans master pulses to downstream boards, then steps
   *          the foreground build/commit/pacing model and updates probes.
   */
  void run_tick(SimBoard &b, uint64_t tg) {
    // This wake is consumed: the grid resumes at the next slot strictly
    // after it (a mask may have swallowed several slots — they coalesced
    // into this one delayed wake).
    do {
      b.next_tick += b.tick_step;
    } while (b.next_tick <= double(tg));
    if (b.edge_latched && !masked_until(b, tg)) {
      b.board.on_sync_edge(local_now(b, tg)); // delayed, merged timestamp
      b.edge_latched = false;
    }

    const uint32_t now = local_now(b, tg);
    BurstSnapshot s;
    const BurstSnapshot *sp = nullptr;
    if (!b.master &&
        b.board.mailbox().try_claim(now, b.board.gap_timeout_cycles(),
                                    b.board.max_burst_cycles(), &s))
      sp = &s;
    const TickActions a = b.board.tick(now, sp);
    if (a.pulse && b.master)
      for (size_t j = 1; j < boards.size(); ++j)
        deliver_edge(static_cast<int>(j), tg);

    // Foreground model (build requests, commit/join swaps, frame pacing).
    const uint32_t bw = b.board.build_word();
    const uint32_t gen = SyncBoard::build_gen_of(bw);
    if (gen != b.seen_gen) {
      b.seen_gen = gen;
      // Release + delete the outgoing instance.
      b.handoff.request_release();
      b.handoff.clear_pending();
      b.pending_index = SyncBoard::build_index_of(bw);
      // Mirror the device foreground: the RNG restart at build pickup seeds
      // from the epoch index (pov_segmented.h).
      b.build_seed = hs::epoch_seed(static_cast<uint32_t>(b.pending_index));
      b.pending_gen = gen;
      b.pending_ready_g = tg + b.init_delay;
      b.have_pending = true;
    }
    // Construction completes: the instance is published for the ISR to adopt.
    if (b.have_pending && b.pending_gen == gen && tg >= b.pending_ready_g)
      b.handoff.publish(&b.instance, b.pending_gen);

    // The device's flywheel-ISR sequence, verbatim (hardware/pov_handoff.h).
    const auto w = b.handoff.apply_wake({.commit = a.commit,
                                         .join_boundary = a.join_boundary,
                                         .dark = a.dark,
                                         .flip = a.flip,
                                         .zero_crossing = a.zero_crossing,
                                         .wire_gen = gen});
    if (!w.commit_ok)
      b.trapped = true; // device: HS_CHECK fires
    if (w.adopted) {
      b.live_index = b.pending_index;
      b.t = 0;
      b.swap_g = tg;
    }
    b.live = w.live != nullptr;
    if (a.flip) {
      ++b.flips;
      if (w.advance)
        ++b.t;
    }
    b.dark_now = w.dark;
  }
};

/**
 * @brief Predicate: every board has gone live (boot join complete).
 * @param s The simulation to inspect.
 * @return True iff no board is still pre-live.
 */
inline bool all_boards_live(Sim &s) {
  for (auto &b : s.boards)
    if (!b.live)
      return false;
  return true;
}

/**
 * @brief Runs the sim until every board has gone live (boot join complete).
 * @param sim The simulation to advance.
 * @param cfg The active config (supplies the join budget).
 * @return True if all boards went live within the join budget.
 */
inline bool boot_join(Sim &sim, const Config &cfg) {
  return sim.run_until(all_boards_live, double(cfg.join_grid_revs) + 2);
}

/**
 * @brief Runs the sim to one revolution before the train start.
 * @param sim The simulation to advance.
 * @param cfg The active config (supplies the effect-length budget).
 * @return True if the pre-train point was reached within budget.
 * @details The master (ppm 0) crosses ZERO on exact rev multiples, so the
 *          primary copy rides the crossing ≈ one rev from the moment this
 *          predicate fires.
 */
inline bool to_pre_train(Sim &sim, const Config &cfg) {
  return sim.run_until(
      [](Sim &s) {
        return s.boards[0].board.content().rev_in_effect >=
               s.cfg.revs_per_effect - 1;
      },
      double(cfg.revs_per_effect) + 4);
}

// ── Scenario: clean 4-board run (boot join, phase, flips, wrap) ─────────────

/**
 * @brief Verifies a clean 4-board run: every board is born tens of columns out
 *        of phase, locks within a revolution, joins live at the same boundary
 *        with the same effect and frame counter, and holds sub-2-column phase,
 *        equal frame counters, and ~2 flips/rev through a 32-bit clock wrap —
 *        with no gate rejections, invalid symbols, or traps.
 */
inline void test_sim_boot_and_phase() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 20, -20, 40};
  // Local clocks start just below the 32-bit wrap: every board's CYCCNT
  // wraps ~10 revolutions in, mid-run (§12 timebase arithmetic).
  Sim sim(cfg, 4, ppm, 0xFFFFFFFFull - 10ull * 2 * PERIOD + 12345);

  // Birth phase, before a single symbol: every downstream board is tens of
  // columns from the master, so the sub-2-column oracles below measure
  // acquisition, not crystal drift. Bounded away from the gate too — none of
  // these could be closed by a LOCKED snap.
  for (int i = 1; i < 4; ++i) {
    const int32_t born = circ_dist(sim.board_pos(i), sim.board_pos(0), cfg.W);
    HS_EXPECT_GT(born, cfg.gate_cols);
    HS_EXPECT_GE(born, 40);
  }

  // All boards lock within the first revolution (two boundary symbols).
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.board.lock() != LockState::LOCKED)
            return false;
        return true;
      },
      1.5));

  // Boot join: every board goes live at the SAME join-grid boundary with
  // identical effect and frame counter (no master head start).
  HS_EXPECT_TRUE(boot_join(sim, cfg));
  for (int i = 0; i < 4; ++i) {
    HS_EXPECT_EQ(sim.boards[i].live_index, 0);
    // Swaps happen at each board's own crossing of the same boundary —
    // within a couple of columns of global time of each other.
    const int64_t dg = static_cast<int64_t>(sim.boards[i].swap_g) -
                       static_cast<int64_t>(sim.boards[0].swap_g);
    HS_EXPECT_LE(dg < 0 ? -dg : dg, int64_t(3) * COL);
  }

  // Run through the cycle-counter wrap and beyond; check phase + flip
  // cadence + frame-counter equality at stable mid-half instants.
  for (int slice = 0; slice < 6; ++slice) {
    uint64_t flips_before[4];
    for (int i = 0; i < 4; ++i)
      flips_before[i] = sim.boards[i].flips;
    sim.run_revs(4.0);
    // Sample at a stable point: master mid-first-half (~x=72), where every
    // board's ZERO flip for this rev has long settled.
    HS_EXPECT_TRUE(
        sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
    HS_EXPECT_LE(sim.max_phase_err(), 2);
    for (int i = 0; i < 4; ++i) {
      const uint64_t df = sim.boards[i].flips - flips_before[i];
      HS_EXPECT_GE(df, 8u); // ~2 flips/rev over the ≥4-rev slice
      HS_EXPECT_LE(df, 12u);
      HS_EXPECT_EQ(sim.boards[i].live_index, sim.boards[0].live_index);
      HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
    }
  }
  for (int i = 1; i < 4; ++i) {
    const Telemetry &tm = sim.boards[i].board.telemetry_snapshot();
    HS_EXPECT_EQ(tm.symbols_rejected_gate, 0u);
    HS_EXPECT_EQ(tm.symbols_discarded_invalid, 0u);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }
}

/**
 * @brief Verifies eight boards acquire, join, and remain content-coherent.
 */
inline void test_sim_eight_board_boot_and_phase() {
  const Config cfg = test_config();
  const int32_t PPM[8] = {0, 20, -20, 40, -35, 15, -10, 30};
  Sim sim(cfg, 8, PPM);

  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.board.lock() != LockState::LOCKED)
            return false;
        return true;
      },
      1.5));
  HS_EXPECT_TRUE(boot_join(sim, cfg));

  sim.run_revs(8.0);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  HS_EXPECT_LE(sim.max_phase_err(), 2);
  for (size_t i = 1; i < sim.boards.size(); ++i) {
    HS_EXPECT_EQ(sim.boards[i].live_index, sim.boards[0].live_index);
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }
}

// ── Scenario: epoch commit — lockstep advance, dark window, deadline ───────

/**
 * @brief Verifies epoch commit lockstep: all four boards play through the
 *        announce phase, go dark together for the K-rev construction window,
 *        then swap to the next effect at the same boundary with frame counters
 *        re-zeroed together; the cadence holds across a second epoch (roster
 *        wraps mod effect_count).
 */
inline void test_sim_epoch_commit() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 30, -25, 10};
  Sim sim(cfg, 4, ppm);

  HS_EXPECT_TRUE(boot_join(sim, cfg));

  // Run to the train start. Through the announce phase (the first R revs of
  // the B+R+K countdown) every board keeps playing the outgoing effect…
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return s.boards[0].board.content().commit_pending; },
      double(cfg.revs_per_effect) + 2));
  sim.run_revs(1.0); // mid-announce
  for (auto &b : sim.boards) {
    HS_EXPECT_TRUE(b.board.content().commit_pending);
    HS_EXPECT_FALSE(b.dark_now);
  }
  // …then all go dark together for the K-revolution construction window.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return s.boards[0].board.content().commit_in_revs <= s.cfg.commit_revs;
      },
      double(cfg.epoch_repeats) + 1));
  sim.run_revs(1.0); // mid-construction (K = 2)
  for (auto &b : sim.boards) {
    HS_EXPECT_TRUE(b.board.content().commit_pending);
    HS_EXPECT_TRUE(b.dark_now);
  }

  // Commit: all four swap to effect 1 at the same boundary, frame counters
  // reset together; the EPOCH redundancy repeats were refractory-ignored.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 1)
            return false;
        return true;
      },
      4.0));
  for (int i = 1; i < 4; ++i) {
    const int64_t dg = static_cast<int64_t>(sim.boards[i].swap_g) -
                       static_cast<int64_t>(sim.boards[0].swap_g);
    HS_EXPECT_LE(dg < 0 ? -dg : dg, int64_t(3) * COL);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }
  // Post-epoch: full content coherence (index AND t) including the master.
  sim.run_revs(3.0);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  for (int i = 1; i < 4; ++i) {
    HS_EXPECT_EQ(sim.boards[i].live_index, 1);
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
  }

  // Second epoch keeps the cadence (roster wraps mod effect_count).
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 2)
            return false;
        return true;
      },
      double(cfg.revs_per_effect) + 6));
}

/** @brief Verifies the master schedules each epoch from its roster duration. */
inline void test_sim_variable_effect_durations() {
  uint32_t effect_revolutions[4] = {24, 52, 32, 64};
  Config cfg = test_config();
  cfg.effect_revolutions = effect_revolutions;
  const int32_t ppm[4] = {0, 30, -25, 10};
  Sim sim(cfg, 4, ppm);

  HS_EXPECT_TRUE(boot_join(sim, cfg));
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 1)
            return false;
        return true;
      },
      double(effect_revolutions[0]) + 6));

  sim.run_revs(45.0);
  for (auto &b : sim.boards)
    HS_EXPECT_EQ(b.live_index, 1);

  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 2)
            return false;
        return true;
      },
      12.0));
}

/**
 * @brief Verifies every board derives the same epoch seed at each handoff: the
 *        boot build uses the identity seed (epoch 0 == 1337) and the first
 *        epoch's build seed is fresh yet identical on all boards.
 */
inline void test_sim_epoch_seed_lockstep() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 30, -25, 10};
  Sim sim(cfg, 4, ppm);

  HS_EXPECT_TRUE(boot_join(sim, cfg));
  for (auto &b : sim.boards)
    HS_EXPECT_EQ(b.build_seed, 1337u); // epoch 0: identity seed

  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 1)
            return false;
        return true;
      },
      double(cfg.revs_per_effect) + 8));
  for (auto &b : sim.boards) {
    HS_EXPECT_EQ(b.build_seed, hs::epoch_seed(1));
    HS_EXPECT_NE(b.build_seed, 1337u);
  }
}

/**
 * @brief Verifies an effect whose construction outruns the K-revolution window
 *        traps (HS_CHECK on the device) and never silently skews the show
 *        (§6.1).
 */
inline void test_sim_commit_deadline_trap() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 0, 0, 0};
  Sim sim(cfg, 4, ppm);
  sim.boards[2].init_delay =
      static_cast<uint64_t>(3) * 2 * PERIOD; // 3 revs > K = 2
  // Boot joins are NOT deadline-bound: board 2 simply goes live ~3 revs
  // after the others (next join-grid boundary), without trapping.
  sim.run_revs(double(cfg.join_grid_revs) * 3);
  HS_EXPECT_TRUE(sim.boards[2].live);
  HS_EXPECT_FALSE(sim.boards[2].trapped);
  // The epoch commit IS deadline-bound.
  HS_EXPECT_TRUE(sim.run_until([](Sim &s) { return s.boards[2].trapped; },
                               double(cfg.revs_per_effect) + 6));
  for (int i : {0, 1, 3})
    HS_EXPECT_FALSE(sim.boards[i].trapped);
}

// ── Scenario: masked-IRQ windows (§4.1, §5.2) ───────────────────────────────

/**
 * @brief Verifies masked-IRQ windows (§4.1, §5.2): boundary masks truncate
 *        burst counts (symbol degrades to missed, never misclassified), mid-rev
 *        masks coalesce wakes harmlessly, and a master masked across its own
 *        boundary self-censors rather than emitting late — phase and content
 *        stay coherent throughout.
 */
inline void test_sim_masked_windows() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 25, -25, 15};
  Sim sim(cfg, 4, ppm);

  // Board 2: recurring 3-column masks placed over the ZERO boundary — they
  // swallow wakes (coalesced) AND merge the first two burst edges (single
  // pin latch), so its decoder sees truncated counts: the symbol must
  // degrade to "missed" (invalid/discarded), never "misclassified".
  // Masks start after boot join so acquisition is clean.
  const uint64_t rev = 2ull * PERIOD;
  for (int k = 8; k < 28; ++k) {
    const uint64_t b0 = k * rev; // master ZERO crossings ≈ k·rev (ppm 0)
    sim.boards[2].masks.push_back({b0 - COL / 2, b0 + 2 * COL + COL / 2});
  }
  // Board 3: mid-revolution masks (no boundary, no symbol) — pure wake
  // coalescing; the flywheel resumes at the time-correct column.
  for (int k = 8; k < 28; ++k) {
    const uint64_t m0 = k * rev + 40 * COL;
    sim.boards[3].masks.push_back({m0, m0 + 5 * COL});
  }

  sim.run_revs(30.0);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));

  // Truncated bursts were discarded (count telemetry), never accepted as
  // the wrong boundary: phase stays sub-column-ish and content equal.
  const Telemetry &tm2 = sim.boards[2].board.telemetry_snapshot();
  HS_EXPECT_GT(tm2.symbols_discarded_invalid, 5u);
  HS_EXPECT_LE(sim.max_phase_err(), 2);
  for (int i = 1; i < 4; ++i) {
    HS_EXPECT_TRUE(sim.boards[i].board.lock() == LockState::LOCKED);
    HS_EXPECT_EQ(sim.boards[i].live_index, sim.boards[0].live_index);
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }

  // Master masked across its own boundary: it self-censors (or truncates)
  // rather than emitting late; downstream coasts one half-rev and re-snaps.
  Sim sim2(cfg, 2, ppm);
  sim2.run_revs(8.0);
  const uint64_t b0 = (static_cast<uint64_t>(sim2.g / rev) + 2) * rev;
  sim2.boards[0].masks.push_back({b0 - COL / 4, b0 + 2 * COL});
  sim2.run_revs(6.0);
  const Telemetry &tm0 = sim2.boards[0].board.telemetry_snapshot();
  HS_EXPECT_GT(tm0.emit_censored + tm0.emit_aborted +
                   tm0.beacons_overrun_dropped + tm0.boundary_bursts_dropped,
               0u);
  HS_EXPECT_LE(sim2.max_phase_err(), 2);
  HS_EXPECT_EQ(sim2.boards[1].board.telemetry_snapshot().symbols_rejected_gate,
               0u);
}

// ── Scenario: EMI on the sync wire (§5.2, §5.3, §9.1) ───────────────────────

/**
 * @brief Verifies EMI on the sync wire (§5.2, §5.3, §9.1): isolated spurious
 *        edges form valid HALF symbols the LOCKED gate rejects, edges injected
 *        inside real bursts corrupt the count to invalid (discarded whole), and
 *        a single edge within G of a predicted boundary is the bounded accepted
 *        case — none unlock or desync the show.
 */
inline void test_sim_emi() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 20, -30, 5};
  Sim sim(cfg, 4, ppm);
  const uint64_t rev = 2ull * PERIOD;

  // Isolated spurious edges on board 1, ~1 per revolution at varied
  // mid-revolution offsets (away from boundaries): each forms a 1-pulse
  // burst = a valid HALF symbol, which the LOCKED plausibility gate must
  // reject (implied correction ≫ G).
  uint32_t lcg = 12345;
  for (int k = 8; k < 40; ++k) {
    lcg = lcg * 1664525u + 1013904223u;
    const uint64_t off = (20 + lcg % 100) * static_cast<uint64_t>(COL);
    sim.emi.push_back({k * rev + off, 1});
  }
  // Two edges injected INSIDE master ZERO bursts on board 2 (count 3 → 4):
  // even count = invalid, discarded whole; the crossing flip covers it.
  sim.emi.push_back({10 * rev + COL, 2});
  sim.emi.push_back({14 * rev + COL, 2});
  // One edge within G of board 3's predicted HALF boundary: the §9.1
  // accepted-EMI case — a ≤G-column seam for ≤½ rev, re-snapped by the next
  // real symbol. It must not unlock or misclassify anything.
  sim.emi.push_back({12 * rev + (PERIOD - 2ull * COL), 3});
  std::sort(sim.emi.begin(), sim.emi.end());

  sim.run_revs(42.0);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));

  const Telemetry &tm1 = sim.boards[1].board.telemetry_snapshot();
  HS_EXPECT_GT(tm1.symbols_rejected_gate, 20u); // isolated EMI all rejected
  const Telemetry &tm2 = sim.boards[2].board.telemetry_snapshot();
  HS_EXPECT_GE(tm2.symbols_discarded_invalid, 2u); // corrupted bursts dropped
  HS_EXPECT_LE(sim.max_phase_err(), 2);
  for (int i = 1; i < 4; ++i) {
    HS_EXPECT_TRUE(sim.boards[i].board.lock() == LockState::LOCKED);
    HS_EXPECT_EQ(sim.boards[i].live_index, sim.boards[0].live_index);
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
  }
}

// ── Scenario: dropped symbols → coast; dropped epoch → beacon fix (§6.3) ───

/**
 * @brief Verifies dropped-symbol recovery (§6.3): a multi-rev symbol gap is a
 *        coast on the crystal that silently re-snaps, and a board that misses
 *        an entire EPOCH train stays visibly stale on the old effect until the
 *        next index beacon corrects it and it rejoins on the join grid — never
 *        a wrong frame, never a trap.
 */
inline void test_sim_drops_and_missed_epoch() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 35, -35, 20};
  Sim sim(cfg, 4, ppm);
  const uint64_t rev = 2ull * PERIOD;

  // Boot join first.
  HS_EXPECT_TRUE(boot_join(sim, cfg));

  // Plain symbol drop: board 1 hears nothing for 2 revolutions — it coasts
  // on its crystal (telemetry max_coast) and silently re-snaps after.
  sim.boards[1].drop_from = sim.g + rev;
  sim.boards[1].drop_to = sim.g + 3 * rev;
  sim.run_revs(5.0);
  HS_EXPECT_GE(sim.boards[1].board.telemetry_snapshot().max_coast_halves, 4u);
  HS_EXPECT_TRUE(sim.boards[1].board.lock() == LockState::LOCKED);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  HS_EXPECT_LE(sim.max_phase_err(), 2);

  // Missed epoch: board 3 loses its wire for the entire EPOCH train (the
  // primary copy + all R repeats). Peers advance; board 3 stays on the old
  // effect until the next index beacon corrects it (§6.3.2), then rejoins
  // on the join grid.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return s.boards[0].board.content().rev_in_effect >=
               s.cfg.revs_per_effect - 1;
      },
      double(cfg.revs_per_effect) + 4));
  sim.boards[3].drop_from = sim.g;
  sim.boards[3].drop_to = sim.g + 7 * rev; // covers B..B+K and the repeats
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return s.boards[0].live_index == 1 && s.boards[1].live_index == 1;
      },
      8.0));
  HS_EXPECT_EQ(sim.boards[3].live_index, 0); // visibly stale, as budgeted
  // Correction: ≤ one beacon period + join grid after the wire returns.
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.boards[3].live_index == 1; },
                    double(cfg.beacon_period_revs + cfg.join_grid_revs) + 6));
  HS_EXPECT_EQ(
      sim.boards[3].board.telemetry_snapshot().beacon_index_corrections, 1u);
  HS_EXPECT_FALSE(sim.boards[3].trapped);
}

// ── Scenario: mid-show reboot — fail-dark, rejoin at correct effect ─────────

/**
 * @brief Verifies mid-show reboot: a board reseeded with fresh engine state
 *        re-acquires phase within ~one boundary symbol, stays dark through
 *        ACQUIRE, then rejoins at the master's current effect from the next
 *        beacon + join-grid boundary — never a wrong frame in between.
 */
inline void test_sim_reboot() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 15, -15, 30};
  Sim sim(cfg, 4, ppm);
  sim.run_revs(12.0);

  // Reboot board 2 mid-show: fresh engine state, no identity assumption.
  SimBoard &b2 = sim.boards[2];
  b2.board.seed(Sim::local_now(b2, sim.g), false);
  // A rebooted board's ISR holds no effect and has consumed no generation.
  b2.handoff.adopt(nullptr, 0);
  b2.handoff.clear_pending();
  b2.seen_gen = 0;
  b2.have_pending = false;
  b2.live = false;
  b2.live_index = -1;
  b2.t = 0;
  b2.flips = 0;
  b2.trapped = false;

  HS_EXPECT_TRUE(b2.board.lock() == LockState::ACQUIRE);
  HS_EXPECT_TRUE(b2.dark_now || !b2.live); // dark through ACQUIRE

  // Phase re-acquires within ~one boundary symbol.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return s.boards[2].board.lock() == LockState::LOCKED; },
      1.5));
  // Identity from the next beacon, display from the next join-grid boundary.
  // This reboot lands clear of the commit window, so the deadline is the
  // blackout-free part of rejoin_bound_revs(); never a wrong frame in between
  // (dark throughout).
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.boards[2].live; },
                    double(cfg.beacon_period_revs + cfg.join_grid_revs) + 4));
  HS_EXPECT_EQ(b2.live_index, sim.boards[0].live_index);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  HS_EXPECT_LE(sim.max_phase_err(), 2);
}

// ── Scenario: forged plausible burst (§8.4 spurious-flip hole, closed) ──────

/**
 * @brief Verifies the forged-plausible-burst defense (§8.4): the strongest
 *        cheap spurious-flip attack is held as a suspect and rejected, never
 *        snapped or flipped; flip cadence and content stay intact (attack
 *        construction detailed in the body).
 */
inline void test_sim_forged_burst() {
  const Config cfg = test_config();
  const int32_t ppm[2] = {0, 10};
  Sim sim(cfg, 2, ppm);
  sim.run_revs(10.0);
  const uint64_t rev = 2ull * PERIOD;

  // A spurious flip needs a forged burst that is simultaneously valid,
  // plausible, and boundary-consistent (§5.3). Forge the strongest cheap
  // attack: an isolated valid-count (HALF) burst 30 columns past a real
  // ZERO boundary — clear of the real burst's gap-timeout window and of the
  // quiet-before guard, far from both predicted boundaries. It must be held
  // as a suspect and counted as a rejection, never snapped or flipped.
  const uint64_t b0 = (sim.g / rev + 2) * rev; // a future master ZERO
  sim.emi.push_back({b0 + 30 * static_cast<uint64_t>(COL), 1});
  sim.emi_pos = 0;
  std::sort(sim.emi.begin(), sim.emi.end());

  const uint64_t flips_before = sim.boards[1].flips;
  const uint32_t rejected_before =
      sim.boards[1].board.telemetry_snapshot().symbols_rejected_gate;
  sim.run_revs(4.0);
  HS_EXPECT_GT(sim.boards[1].board.telemetry_snapshot().symbols_rejected_gate,
               rejected_before);
  // Flip cadence unbroken: 2 per revolution, no extra content advance.
  const uint64_t df = sim.boards[1].flips - flips_before;
  HS_EXPECT_GE(df, 7u);
  HS_EXPECT_LE(df, 9u);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  HS_EXPECT_EQ(sim.boards[1].t, sim.boards[0].t);
}

// ── Scenario: epoch-repeat lockstep (§6.3.1) ────────────────────────────────

/**
 * @brief Verifies the EPOCH redundancy repeats are a reliability mechanism, not
 *        a per-board reschedule (§6.3.1): a board that misses the primary copy
 *        at boundary B but accepts a repeat at B+j must still commit at the
 *        SAME absolute boundary as its peers, with a frame counter that stays
 *        equal afterwards.
 * @details Two sub-scenarios: one downstream board deafened for exactly the
 *          primary copy, and the master masked across B so it self-censors the
 *          primary — every downstream board then first hears the B+1 repeat
 *          while the master heard "itself" at B.
 */
inline void test_sim_epoch_repeat_lockstep() {
  const Config cfg = test_config();
  const uint64_t rev = 2ull * PERIOD;

  /**
   * @brief Predicate: all boards have committed to effect 1.
   * @param s The simulation to inspect.
   * @return True iff every board's live_index is 1.
   */
  auto all_on_effect_1 = [](Sim &s) {
    for (auto &b : s.boards)
      if (b.live_index != 1)
        return false;
    return true;
  };
  /**
   * @brief Asserts every downstream board committed at the same boundary as the
   *        master (≤3 columns apart, no trap) and holds an equal frame counter
   *        afterward.
   * @param sim The simulation to inspect (advanced two revs while checking).
   */
  auto expect_lockstep = [&](Sim &sim) {
    for (int i = 1; i < 4; ++i) {
      // Commits land at each board's own crossing of the same boundary —
      // within a couple of columns of global time, never a revolution apart.
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
  };
  const double commit_revs_max =
      double(cfg.commit_revs + static_cast<uint32_t>(cfg.epoch_repeats)) + 8;

  // Sub-scenario 1: board 2 misses ONLY the primary copy at B.
  {
    const int32_t ppm[4] = {0, 20, -20, 10};
    Sim sim(cfg, 4, ppm);
    HS_EXPECT_TRUE(boot_join(sim, cfg));
    HS_EXPECT_TRUE(to_pre_train(sim, cfg));
    sim.boards[2].drop_from = sim.g + rev - 4 * COL;
    sim.boards[2].drop_to = sim.g + rev + rev / 4; // before the B+1 repeat
    HS_EXPECT_TRUE(sim.run_until(all_on_effect_1, commit_revs_max));
    expect_lockstep(sim);
  }

  // Sub-scenario 2: the master self-censors the primary copy (masked across
  // its own train-start boundary, a designed event — §5.2): it must not
  // schedule a commit its peers cannot match.
  {
    const int32_t ppm[4] = {0, 15, -25, 30};
    Sim sim(cfg, 4, ppm);
    HS_EXPECT_TRUE(boot_join(sim, cfg));
    HS_EXPECT_TRUE(to_pre_train(sim, cfg));
    const uint64_t b0 = sim.g + rev;
    sim.boards[0].masks.push_back({b0 - COL / 4, b0 + 2 * COL});
    HS_EXPECT_TRUE(sim.run_until(all_on_effect_1, commit_revs_max));
    const Telemetry &tm0 = sim.boards[0].board.telemetry_snapshot();
    HS_EXPECT_GT(tm0.emit_censored + tm0.emit_aborted +
                     tm0.beacons_overrun_dropped + tm0.boundary_bursts_dropped,
                 0u);
    expect_lockstep(sim);
  }
}

// ── Scenario: schedule-counter resync from the beacon rev cross-check ───────

/**
 * @brief Verifies the §6.4 beacon rev cross-check resyncs a slipped
 *        schedule counter within one beacon period, restoring exact
 *        j-inference before the next train.
 * @details A board whose rev_in_effect slipped against the master — a late
 *          commit through the §6.3.1 j-fallback, or a crossing hiccup while a
 *          corrupted timebase recovered — would infer j wrongly at every later
 *          epoch train and commit out of lockstep by the slip, in either
 *          direction, indefinitely (its own commit re-zeros the counter at its
 *          own, offset boundary, so the slip is self-sustaining). The
 *          cross-check corrects via the signed mod-64 difference, leaving
 *          content t untouched.
 */
inline void test_sim_rev_resync() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 20, -20, 10};
  Sim sim(cfg, 4, ppm);
  const uint64_t rev = 2ull * PERIOD;
  HS_EXPECT_TRUE(boot_join(sim, cfg));
  // Settle past the boot beacons (revs 1–3), then slip board 3's counter
  // by +2: at the next train it would over-count j by 2 and commit 2
  // revolutions EARLY.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return s.boards[0].board.content().rev_in_effect == 5; },
      12.0));
  content_mut(sim.boards[3].board).rev_in_effect += 2;

  // Detected and corrected at the next beacon (rev 9), within one period.
  // (Pre-resync the counters differ by 2 — a crossing straddle changes the
  // gap by at most 1, so this predicate cannot fire spuriously.)
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return s.boards[3].board.content().rev_in_effect ==
               s.boards[0].board.content().rev_in_effect;
      },
      double(cfg.beacon_period_revs) + 2));
  HS_EXPECT_GE(sim.boards[3].board.telemetry_snapshot().beacon_rev_mismatches,
               1u);

  // The next epoch commits in lockstep even though board 3 ALSO misses the
  // primary copy, taking the repeat path that depends on j-inference.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return s.boards[0].board.content().rev_in_effect >=
               s.cfg.revs_per_effect - 1;
      },
      double(cfg.revs_per_effect) + 4));
  sim.boards[3].drop_from = sim.g + rev - 4 * COL;
  sim.boards[3].drop_to = sim.g + rev + rev / 4;
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 1)
            return false;
        return true;
      },
      double(cfg.commit_revs + static_cast<uint32_t>(cfg.epoch_repeats)) + 8));
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

// ── Scenario: 6-bit rev counter wraps WITHIN one effect (§6.4 mod-64) ────────

/**
 * @brief Verifies normal play stays in lockstep — and the epoch still commits
 *        correctly — as rev_in_effect rolls through its 6-bit (mod-64) residue
 *        within a single effect.
 * @details The beacon carries rev mod 64; the cross-check compares
 *          f.rev_count against `content_tracker.rev_in_effect & 63`, and the
 *          beacon_period_revs < 32 rule (Config::valid) exists precisely so the
 *          resulting signed-mod-64 resync is unambiguous as the residue wraps.
 *          The 40-rev sim configs (and the 5-rev slip in test_sim_rev_resync)
 *          never reach rev 64, so this dynamic wrap shipped unverified; here a
 *          90-rev effect crosses rev 64 mid-show. A dropped `& 63` would, via
 *          handle_beacon_burst's fold, "resync" a rev-64 board's counter
 *          backwards (64→32) against a beacon rev_count of 0 — corrupting both
 *          phase lockstep and the j = rev_in_effect − revs_per_effect inference
 *          the next epoch train depends on. This pins both staying correct.
 */
inline void test_sim_rev_wrap_within_effect() {
  Config cfg = test_config();
  cfg.revs_per_effect =
      90; // > 64: rev_in_effect wraps its 6-bit residue mid-effect
  const int32_t ppm[4] = {0, 20, -20, 10};
  Sim sim(cfg, 4, ppm);

  // Boot: all four boards live on effect 0.
  HS_EXPECT_TRUE(boot_join(sim, cfg));
  for (int i = 0; i < 4; ++i)
    HS_EXPECT_EQ(sim.boards[i].live_index, 0);

  // Advance ACROSS the 63→0 seam and hold past it, sampling lockstep at stable
  // mid-half instants (master ~x=72). A broken &63 cross-check corrupts
  // rev_in_effect at rev 64 and would break phase / frame-counter lockstep here.
  bool crossed_seam = false;
  for (;;) {
    HS_EXPECT_TRUE(
        sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
    const uint32_t rev = sim.boards[0].board.content().rev_in_effect;
    HS_EXPECT_LE(sim.max_phase_err(), 2);
    for (int i = 0; i < 4; ++i) {
      HS_EXPECT_EQ(sim.boards[i].live_index, 0);
      HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
      HS_EXPECT_FALSE(sim.boards[i].trapped);
      // The schedule counter tracks the master's exactly through the wrap, and
      // the &63 cross-check raises no spurious rev mismatch. A dropped mask on
      // the comparison reads a rev≥64 board's full counter as differing from the
      // wrapped beacon rev_count and resyncs every beacon (climbing
      // beacon_rev_mismatches); a dropped mask in the fold itself diverges the
      // counter from the master. Phase/frame lockstep alone catches neither —
      // those derive from the flywheel, not rev_in_effect.
      HS_EXPECT_EQ(sim.boards[i].board.content().rev_in_effect, rev);
      HS_EXPECT_EQ(
          sim.boards[i].board.telemetry_snapshot().beacon_rev_mismatches, 0u);
    }
    if (rev >= 64)
      crossed_seam = true;
    if (rev >= 80)
      break;
    sim.run_revs(3.0);
  }
  HS_EXPECT_TRUE(crossed_seam); // the run actually exercised rev_in_effect ≥ 64

  // The epoch still commits in lockstep at the post-seam effect boundary:
  // on_epoch_symbol infers j from a rev whose 6-bit residue has wrapped.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 1)
            return false;
        return true;
      },
      double(cfg.revs_per_effect) + 8));
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  for (int i = 1; i < 4; ++i) {
    HS_EXPECT_EQ(sim.boards[i].live_index, sim.boards[0].live_index);
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }
}

// ── Scenario: same-tick EPOCH burst + boundary fold (§6.3.1 j-inference) ─────

/**
 * @brief Verifies §6.3.1 j-inference stays correct when an EPOCH burst is
 *        consumed in the SAME tick() its boundary is folded.
 * @details on_epoch_symbol infers j = rev_in_effect − revs_per_effect, so it
 *          MUST observe the POST-fold rev. handle_burst guarantees this: its
 *          backstop apply_flip(ZERO) folds rev_in_effect (via on_zero_crossing)
 *          before the ZERO_EPOCH branch, and tick()'s later fold loop — which
 *          runs after handle_burst — is deduped by the flip gate. This pins the
 *          ordering directly on the content tracker: driving the exact
 *          fold-then-infer sequence handle_burst uses, every copy j of the train
 *          must commit at the same absolute B+R+K boundary whether the fold was
 *          deferred into the burst's tick or already applied a tick earlier. A
 *          closing assertion pins the precondition the backstop satisfies — the
 *          un-folded (pre-fold) rev mis-infers a repeat copy one short.
 */
inline void test_epoch_same_tick_burst_fold() {
  const Config cfg = test_config();
  const uint32_t RPE = cfg.revs_per_effect;
  const uint32_t R = static_cast<uint32_t>(cfg.epoch_repeats);
  const uint32_t K = cfg.commit_revs;

  // Drives a content tracker that hears copy j of the train, then counts ZERO
  // crossings to the commit and returns the ABSOLUTE rev at which it fired
  // (on_zero_crossing zeroes rev_in_effect on commit, so it is computed, not
  // read back). same_tick=true models the fold deferred into the burst's tick
  // (the backstop folds first); same_tick=false models the fold already applied
  // a tick earlier (the backstop is deduped — modeled by not folding again).
  auto commit_rev_for = [&](uint32_t j, bool same_tick) -> uint32_t {
    ContentTracker c;
    c.identity_known = true;
    c.effect_index = 0;
    if (same_tick) {
      c.rev_in_effect = RPE + j - 1;            // B+j fold still pending…
      HS_EXPECT_FALSE(c.on_zero_crossing(cfg)); // …backstop apply_flip folds it
    } else {
      c.rev_in_effect = RPE + j; // already folded a tick earlier
    }
    const uint32_t base_rev = c.rev_in_effect; // == RPE + j either way
    HS_EXPECT_EQ(base_rev, RPE + j);
    HS_EXPECT_TRUE(c.on_epoch_symbol(cfg)); // opens the commit window
    HS_EXPECT_EQ(c.commit_in_revs, K + R - j);
    uint32_t count = 0;
    while (count <= RPE) {
      const bool committed = c.on_zero_crossing(cfg);
      ++count;
      if (committed)
        break;
    }
    return base_rev + count;
  };

  // Every copy commits at the absolute B+R+K boundary, independent of j and of
  // whether the burst shared its tick with the fold — the lockstep §6.3.1
  // promises.
  for (uint32_t j = 0; j <= R; ++j) {
    HS_EXPECT_EQ(commit_rev_for(j, /*same_tick=*/true), RPE + R + K);
    HS_EXPECT_EQ(commit_rev_for(j, /*same_tick=*/false), RPE + R + K);
  }

  // Why the backstop fold is load-bearing: on_epoch_symbol on the PRE-fold rev
  // infers a repeat copy one short (j−1), scheduling commit_in_revs a revolution
  // too large — exactly the late commit handle_burst's apply_flip(ZERO) prevents
  // by folding first.
  for (uint32_t j = 1; j <= R; ++j) {
    ContentTracker pre;
    pre.identity_known = true;
    pre.rev_in_effect = RPE + j - 1; // fold NOT applied before the inference
    HS_EXPECT_TRUE(pre.on_epoch_symbol(cfg));
    HS_EXPECT_EQ(pre.commit_in_revs, K + R - (j - 1)); // one too large
  }
}

/**
 * @brief Pins ContentTracker::construction_opens and ::constructing directly:
 *        the window opens exactly once and lasts exactly K revolutions, for
 *        every repeat copy the board may have heard.
 * @details The simulator only observes these through dark_now, which cannot
 *          separate "the window opens here" from "the window is open", nor say
 *          which predicate an off-by-one came from. The window is anchored to
 *          the absolute commit boundary, so a board that heard the last repeat
 *          (announce phase already spent) opens it at the accept itself rather
 *          than at a later crossing — both spellings are covered here.
 */
inline void test_construction_window_predicates() {
  const Config cfg = test_config();
  const uint32_t RPE = cfg.revs_per_effect;
  const uint32_t R = static_cast<uint32_t>(cfg.epoch_repeats);
  const uint32_t K = cfg.commit_revs;

  for (uint32_t j = 0; j <= R; ++j) {
    HS_CONTEXT("copy", static_cast<long long>(j));
    ContentTracker c;
    c.identity_known = true;
    c.effect_index = 0;
    c.rev_in_effect = RPE + j;

    // Nothing scheduled: neither predicate can fire off commit_pending.
    HS_EXPECT_FALSE(c.construction_opens(cfg));
    HS_EXPECT_FALSE(c.constructing(cfg));

    HS_EXPECT_TRUE(c.on_epoch_symbol(cfg));
    const uint32_t announce = R - j; // crossings before the window opens
    HS_EXPECT_EQ(c.commit_in_revs, K + announce);
    HS_EXPECT_EQ(c.construction_opens(cfg), announce == 0);
    HS_EXPECT_EQ(c.constructing(cfg), announce == 0);

    uint32_t opens = c.construction_opens(cfg) ? 1u : 0u;
    uint32_t dark = c.constructing(cfg) ? 1u : 0u;
    uint32_t steps = 0;
    bool committed = false;
    while (!committed && steps <= K + R + 1) {
      committed = c.on_zero_crossing(cfg);
      ++steps;
      // Opening is a moment inside the window, never outside it.
      HS_EXPECT_TRUE(!c.construction_opens(cfg) || c.constructing(cfg));
      opens += c.construction_opens(cfg) ? 1u : 0u;
      dark += c.constructing(cfg) ? 1u : 0u;
    }
    HS_EXPECT_TRUE(committed);
    HS_EXPECT_EQ(steps, K + announce);
    HS_EXPECT_EQ(opens, 1u);
    HS_EXPECT_EQ(dark, K);

    // The commit closes the window; a board past it is lit again.
    HS_EXPECT_FALSE(c.construction_opens(cfg));
    HS_EXPECT_FALSE(c.constructing(cfg));
  }
}

// ── §9.1 failure-mode budget: artifact bounds and recovery times ────────────
//
// One scenario per spec §9.1 budget row that the tests above do not already
// pin, asserting BOTH the worst-case artifact bound and the recovery time.
// Coverage map (§9.1 row → test):
//
//   crystal drift (normal operation)  → test_sim_boot_and_phase
//   masked-IRQ windows / self-censor  → test_sim_masked_windows
//   lost boundary symbol              → test_budget_lost_symbol
//   EMI, rejected cases               → test_sim_emi, test_sim_forged_burst
//   EMI, accepted (binding) case      → test_budget_emi_accepted_seam
//   mis-snap / corrupted timebase     → test_budget_corrupted_timebase
//   mis-snap forged during ACQUIRE    → test_budget_acquire_mis_snap
//   dropped render                    → not simulable here: frame pacing
//                                       lives in the device's Canvas
//                                       buffer_free() gate; the epoch-bounded
//                                       t recovery it relies on is pinned by
//                                       test_sim_epoch_commit
//   missed epoch (all copies)         → test_sim_drops_and_missed_epoch
//   epoch repeat lockstep (§6.3.1)    → test_sim_epoch_repeat_lockstep,
//                                       test_master_epoch_train_bounded
//   corrupted beacon frame            → test_budget_beacon_corruption
//   board reboot mid-show             → test_sim_reboot
//   commit deadline (HS_CHECK trap)   → test_sim_commit_deadline_trap
//   sync wire dead / master dead      → test_budget_wire_dead

/**
 * @brief Steps the sim for @p revs, the §9.1 artifact probe.
 * @param sim The simulation to advance.
 * @param revs Number of revolutions to step.
 * @return The worst locked-board phase error (columns vs the master) observed
 *         at any step.
 */
inline int32_t max_err_over(Sim &sim, double revs) {
  const uint64_t until =
      sim.g + static_cast<uint64_t>(revs * 2 * sim.cfg.cycles_per_half_rev);
  int32_t worst = 0;
  while (sim.g < until) {
    sim.step();
    worst = std::max(worst, sim.max_phase_err());
  }
  return worst;
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
  HS_EXPECT_TRUE(sim.boards[1].board.lock() == LockState::LOCKED);
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
  const int32_t seam = max_err_over(sim, 0.45);
  HS_EXPECT_GE(seam, 2);
  HS_EXPECT_LE(seam, cfg.gate_cols);
  // Recovery: the next real boundary symbol (err ≈ 3 ≤ G) re-snaps.
  sim.run_revs(0.45);
  HS_EXPECT_LE(sim.max_phase_err(), 1);
  HS_EXPECT_TRUE(sim.boards[1].board.lock() == LockState::LOCKED);
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
 *          re-snap at the next symbol ≈ 2.8 revolutions ≈ 350 ms at speed.
 */
inline void test_budget_corrupted_timebase() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 20, -20, 10};
  Sim sim(cfg, 4, ppm);
  HS_EXPECT_TRUE(boot_join(sim, cfg));
  // Park at rev 10 (≡ 2 mod the beacon period) so the ACQUIRE window stays
  // clear of beacon trains: this test pins the common path, and the train's
  // isolated first digit — which re-poisons an ACQUIRE board — is
  // test_budget_acquire_mis_snap.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return s.boards[0].board.content().rev_in_effect == 10; },
      16.0));

  // Corrupt board 2's flywheel phase by W/4 — far beyond the gate.
  SimBoard &b2 = sim.boards[2];
  const int32_t bogus = floor_mod(sim.board_pos(2) + 72, cfg.W);
  flywheel_mut(b2.board).seed(Sim::local_now(b2, sim.g) -
                              static_cast<uint32_t>(bogus) * COL);
  flywheel_mut(b2.board).force_lock();
  HS_EXPECT_GE(circ_dist(sim.board_pos(2), sim.board_pos(0), cfg.W), 60);

  const uint32_t rej_before =
      b2.board.telemetry_snapshot().symbols_rejected_gate;
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return s.boards[2].board.lock() == LockState::LOCKED &&
               circ_dist(s.board_pos(2), s.board_pos(0), s.cfg.W) <= 1;
      },
      3.5)); // budget: ~2.8 revs; slack for the mid-rev corruption instant
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
      [](Sim &s) { return s.boards[0].board.content().rev_in_effect == 9; },
      16.0));
  // Reboot 12 columns before the train and past the rev's ZERO burst, so the
  // first wire event the fresh board meets is the train's head digit.
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 60; }, 1.1));
  SimBoard &b2 = sim.boards[2];
  b2.board.seed(Sim::local_now(b2, sim.g), false);
  b2.handoff.adopt(nullptr, 0);
  b2.handoff.clear_pending();
  b2.seen_gen = 0;
  b2.have_pending = false;
  b2.live = false;
  b2.live_index = -1;
  b2.t = 0;
  b2.trapped = false;
  HS_EXPECT_TRUE(b2.board.lock() == LockState::ACQUIRE);

  // The head digit captures it: locked on a beacon digit, a quarter turn out.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return s.boards[2].board.lock() == LockState::LOCKED; },
      0.5));
  const uint64_t snap_g = sim.g;
  HS_EXPECT_GE(circ_dist(sim.board_pos(2), sim.board_pos(0), cfg.W), 60);

  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return s.boards[2].board.lock() == LockState::LOCKED &&
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
  // Fail-dark throughout: never live on an effect the master is not showing.
  HS_EXPECT_TRUE(!b2.live || b2.live_index == sim.boards[0].live_index);
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
      [](Sim &s) { return s.boards[0].board.content().rev_in_effect == 9; },
      16.0));
  // Beacon (index 0, rev 9) is digits [0,0,1,1,2]; its 4th burst is two
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
 * @brief Verifies the §9.1 "sync wire dead" and "master dead" budget rows
 *        (identical for downstream — the master is only the symbol source):
 *        flywheels free-run and keep flipping 2/rev, the playlist freezes on
 *        the current effect (visibly stale, never wrong), and boards precess
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
  // burst is in flight: a half-received burst would register one truncated-
  // count artifact, which is the lost-symbol row's case, not this one.
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
    HS_EXPECT_TRUE(sim.boards[i].board.lock() == LockState::LOCKED);
    HS_EXPECT_EQ(sim.boards[i].board.telemetry_snapshot().symbols_rejected_gate,
                 0u);
    // Layer 2 self-sufficiency: ~2 flips/rev throughout, no stall.
    const uint64_t df = sim.boards[i].flips - flips_before[i];
    HS_EXPECT_GE(df, 2 * static_cast<uint64_t>(coast) - 5);
    HS_EXPECT_LE(df, 2 * static_cast<uint64_t>(coast) + 5);
    // Layer 3 freezes on the current effect.
    HS_EXPECT_EQ(sim.boards[i].live_index, 0);
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
    HS_EXPECT_GE(fade_in, previous_in);
    HS_EXPECT_LE(fade_out, previous_out);
    previous_in = fade_in;
    previous_out = fade_out;
  }
}

// ── Runner ──────────────────────────────────────────────────────────────────

/**
 * @brief Runs every pov_sync test in order.
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
  test_multi_boundary_tick_window();
  test_beacon_codec();
  test_beacon_partial_frame_ages_out();
  test_beacon_shift_needs_confirmation();
  test_beacon_out_of_range_index_rejected();
  test_rev_resync_fold();
  test_flywheel_position();
  test_snap_gate();
  test_suspect_timeout_acquire_uncounted();
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
  test_sim_epoch_seed_lockstep();
  test_sim_commit_deadline_trap();
  test_sim_masked_windows();
  test_sim_emi();
  test_sim_drops_and_missed_epoch();
  test_sim_reboot();
  test_sim_forged_burst();
  test_sim_epoch_repeat_lockstep();
  test_sim_rev_resync();
  test_sim_rev_wrap_within_effect();
  test_epoch_same_tick_burst_fold();
  test_construction_window_predicates();
  test_effect_output_envelope();

  test_budget_lost_symbol();
  test_budget_emi_accepted_seam();
  test_budget_corrupted_timebase();
  test_budget_acquire_mis_snap();
  test_budget_beacon_corruption();
  test_budget_wire_dead();

  return fixture.result();
}

} // namespace pov_sync_tests
} // namespace hs_test
