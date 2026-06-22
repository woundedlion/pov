/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host unit tests for the Phantasm synchronization core (hardware/pov_sync.h)
 * — the spec §12 test plan (docs/phantasm_frame_sync_spec.md).
 *
 * Pure pieces are tested directly: symbol classification, the try_flip state
 * machine, the edge mailbox + glitch filter, the beacon codec, the flywheel's
 * 64-bit position math (incl. cycle-counter wrap and trim extremes), the
 * acceptance gate, and the master emitter's self-censoring.
 *
 * The multi-board scenarios run on a small event-driven simulator: one
 * master + three downstream SyncBoards with per-board crystal offsets, a
 * single-latch masked-IRQ model (an edge during a mask is delayed; two merge
 * — the i.MX RT pin-flag behavior the count coding is designed around),
 * symbol drop windows, EMI injection, a foreground model with effect
 * construction delays, and mid-show reboot.
 */
#pragma once

#include "hardware/pov_sync.h"
#include "tests/test_harness.h"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace hs_test {
namespace pov_sync_tests {

using namespace pov::sync;

// EdgeMailbox::burst_complete() and claim() are the pre-consolidation split
// consumer path, now private behind a test friend (see hardware/pov_sync.h):
// production goes through the fused try_claim(), and the split path is reachable
// only here. EdgeMailboxTestAccess is the friend; the burst_complete()/claim()
// free wrappers below let the tests exercise the two halves try_claim() merged.
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

constexpr uint32_t kPeriod = 37500000u; /**< Cycles per half-rev at full rate. */
constexpr uint32_t kCol = kPeriod / 144u; /**< Cycles per column at full rate. */

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
  HS_EXPECT_TRUE(c.valid());
  HS_EXPECT_EQ(c.cycles_per_half_rev, kPeriod); // 600e6·30/480, exact
  HS_EXPECT_EQ(c.cycles_per_column(), kCol);
  HS_EXPECT_EQ(c.glitch_filter_cycles, 60000u); // 100 µs

  // The beacon's 6-bit rev field resyncs a slip only in (-32, +32), so the
  // beacon period must stay below the half-window: a period >= 32 leaves the
  // resync precondition unenforced. Boundary is exclusive (32 is rejected, 31
  // accepted).
  Config bp = test_config();
  bp.beacon_period_revs = 32;
  HS_EXPECT_FALSE(bp.valid());
  bp.beacon_period_revs = 31;
  HS_EXPECT_TRUE(bp.valid());
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
  // boundary flips again — exactly 2 per simulated revolution.
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
  const uint32_t kGlitch = 60000u;
  EdgeMailbox m;
  HS_EXPECT_FALSE(burst_complete(m, 0, 4 * kCol));
  m.on_edge(1000, kGlitch);
  m.on_edge(1000 + 2 * kCol, kGlitch);
  m.on_edge(1000 + 4 * kCol, kGlitch);
  // EMI spike < 100 µs after an accepted edge is rejected…
  m.on_edge(1000 + 4 * kCol + kGlitch / 2, kGlitch);
  // …and does not reset the filter reference (a spike train can't suppress).
  m.on_edge(1000 + 6 * kCol, kGlitch);
  HS_EXPECT_FALSE(burst_complete(m, 1000 + 7 * kCol, 4 * kCol));
  HS_EXPECT_TRUE(burst_complete(m, 1000 + 10 * kCol + 1, 4 * kCol));
  const BurstSnapshot s = claim(m);
  HS_EXPECT_EQ(s.count, 4u);
  HS_EXPECT_EQ(s.first_cycles, 1000u);
  HS_EXPECT_EQ(s.last_cycles, 1000u + 6 * kCol);
  // Claim resets; the glitch filter still applies across bursts.
  HS_EXPECT_FALSE(burst_complete(m, 1000 + 10 * kCol + 2, 4 * kCol));
  m.on_edge(1000 + 6 * kCol + kGlitch - 1, kGlitch); // too close: rejected
  HS_EXPECT_FALSE(burst_complete(m, 1000 + 20 * kCol, 4 * kCol));

  // try_claim() is the atomic consumer primitive: it recomputes completion
  // from the same count_ it then clears, so the test and the reset cannot be
  // split. It declines an incomplete (or absent) burst and takes a complete
  // one, leaving the same empty state claim() does.
  EdgeMailbox tc;
  BurstSnapshot out;
  HS_EXPECT_FALSE(tc.try_claim(0, 4 * kCol, &out)); // no burst yet
  tc.on_edge(1000, kGlitch);
  tc.on_edge(1000 + 2 * kCol, kGlitch);
  HS_EXPECT_FALSE(tc.try_claim(1000 + 3 * kCol, 4 * kCol, &out)); // gap too short
  HS_EXPECT_TRUE(tc.try_claim(1000 + 6 * kCol + 1, 4 * kCol, &out));
  HS_EXPECT_EQ(out.count, 2u);
  HS_EXPECT_EQ(out.first_cycles, 1000u);
  HS_EXPECT_EQ(out.last_cycles, 1000u + 2 * kCol);
  HS_EXPECT_FALSE(tc.try_claim(1000 + 7 * kCol, 4 * kCol, &out)); // reset
}

/**
 * @brief Verifies the glitch-filter reference does not survive a counter wrap.
 * @details age_prior() (called every column by the flywheel poll) retires the
 *          prior once the wire is quiet past the filter window, so a real edge
 *          after wrap is never falsely rejected by a pseudo-random modular
 *          difference.
 */
inline void test_mailbox_prior_staleness() {
  const uint32_t kGlitch = 60000u;

  // age_prior leaves a within-window reference intact: a genuine spike that
  // arrives before a poll has aged the prior is still suppressed.
  {
    EdgeMailbox m;
    m.on_edge(1000, kGlitch);
    m.age_prior(1000 + kGlitch / 2, kGlitch); // still within the window: kept
    m.on_edge(1000 + kGlitch / 2 + 1, kGlitch); // too close to 1000: rejected
    HS_EXPECT_TRUE(burst_complete(m, 1000 + 100 * kGlitch, 1));
    HS_EXPECT_EQ(claim(m).count, 1u);
  }

  // After the wire goes quiet the prior is retired, so a later edge whose
  // (wrapped) modular distance to the OLD prior lands inside the reject window
  // is still accepted. Without age_prior this edge would be silently dropped.
  {
    EdgeMailbox m;
    const uint32_t prior = 1000u;
    m.on_edge(prior, kGlitch); // a one-edge burst…
    HS_EXPECT_TRUE(burst_complete(m, prior + 10 * kCol, kCol));
    HS_EXPECT_EQ(claim(m).count, 1u); // …claimed; the prior persists across it.
    // The flywheel keeps polling during the silence and retires the stale
    // reference within a column (kCol > kGlitch), long before the counter wraps.
    m.age_prior(prior + 11 * kCol, kGlitch);
    // A real edge after the counter has wrapped: its modular distance to the
    // old prior is only kGlitch/2, which the un-aged filter would reject.
    const uint32_t wrapped = prior + kGlitch / 2;
    m.on_edge(wrapped, kGlitch);
    HS_EXPECT_TRUE(burst_complete(m, wrapped + 10 * kCol, kCol));
    HS_EXPECT_EQ(claim(m).count, 1u); // accepted as a fresh one-edge burst.
  }
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
  // A burst accumulates on the sync wire before the reboot.
  board.on_sync_edge(2000u);
  board.on_sync_edge(2000u + 4 * col);
  HS_EXPECT_TRUE(burst_complete(board.mailbox(), 2000u + 100 * col,
                                cfg.gap_timeout_cycles()));
  // Reboot: reseed. The stale burst must be gone.
  board.seed(3000u, false);
  HS_EXPECT_FALSE(burst_complete(board.mailbox(), 3000u + 100 * col,
                                 cfg.gap_timeout_cycles()));
  HS_EXPECT_EQ(claim(board.mailbox()).count, 0u);
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
      BurstSnapshot s{static_cast<uint32_t>(v) + 1u, t, t + v * kCol};
      bool r = false;
      got = p.feed(s, cfg, out, &r);
      rejected = rejected || r;
      t += 12 * kCol;
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
  // A corrupted digit breaks the checksum; a corrupted checksum likewise.
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
    BurstSnapshot s1{3, 1000, 1000 + 2 * kCol};
    HS_EXPECT_FALSE(p.feed(s1, cfg, &g, &r));
    BurstSnapshot s2{9, 1000 + 12 * kCol, 1000 + 20 * kCol}; // count > 8
    HS_EXPECT_FALSE(p.feed(s2, cfg, &g, &r));
    HS_EXPECT_TRUE(r);
    HS_EXPECT_FALSE(p.active());
  }

  // A stale partial frame (interdigit timeout) is discarded; the burst that
  // exposed it starts a fresh frame which still decodes.
  {
    BeaconParser p;
    BeaconFrame g{};
    bool r = false;
    BurstSnapshot s1{4, 1000, 1000 + 3 * kCol};
    p.feed(s1, cfg, &g, &r);
    encode_beacon_digits(5, 7, d);
    uint32_t t = 1000 + 200 * kCol; // far past the interdigit timeout
    bool got = false;
    bool any_reject = false;
    for (int i = 0; i < 5; ++i) {
      BurstSnapshot s{static_cast<uint32_t>(d[i]) + 1u, t, t + d[i] * kCol};
      bool rr = false;
      got = p.feed(s, cfg, &g, &rr);
      any_reject = any_reject || (rr && i > 0);
      t += 12 * kCol;
    }
    HS_EXPECT_TRUE(got);
    HS_EXPECT_FALSE(any_reject);
    HS_EXPECT_EQ(g.effect_index, 5);
    HS_EXPECT_EQ(g.rev_count, 7u);
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
    const uint32_t period = kPeriod + trim;
    Flywheel f(cfg);
    f.set_cycles_per_half_rev(period);
    f.seed(1000000u);
    bool ok = true;
    for (int64_t delta = 0; delta < period; delta += 12347) {
      const int32_t want =
          static_cast<int32_t>(ref_cols(delta, period) % 288);
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
    HS_EXPECT_EQ(f.position(5000000u - kCol), 287);
    HS_EXPECT_EQ(f.position(5000000u - 3 * kCol - kCol / 2), 284);
  }

  // Fold cadence: exactly one crossing per half-rev, boundaries alternate,
  // each at its exact instant; a long coast yields several crossings.
  {
    Flywheel f(cfg);
    f.seed(1000u);
    HS_EXPECT_FALSE(f.fold(1000u + kPeriod - 1).crossed);
    const Crossing c1 = f.fold(1000u + kPeriod);
    HS_EXPECT_TRUE(c1.crossed);
    HS_EXPECT_TRUE(c1.boundary == Boundary::HALF);
    HS_EXPECT_EQ(c1.at_cycles, 1000u + kPeriod);
    HS_EXPECT_FALSE(f.fold(1000u + kPeriod).crossed);
    // Coast 2.5 half-revs: two more crossings, exact instants.
    const uint32_t late = 1000u + kPeriod + 2 * kPeriod + kPeriod / 2;
    const Crossing c2 = f.fold(late);
    HS_EXPECT_TRUE(c2.crossed && c2.boundary == Boundary::ZERO);
    HS_EXPECT_EQ(c2.at_cycles, 1000u + 2 * kPeriod);
    const Crossing c3 = f.fold(late);
    HS_EXPECT_TRUE(c3.crossed && c3.boundary == Boundary::HALF);
    HS_EXPECT_EQ(c3.at_cycles, 1000u + 3 * kPeriod);
    HS_EXPECT_FALSE(f.fold(late).crossed);
    HS_EXPECT_EQ(f.position(late), 144 + 72);
  }

  // The rebase rule makes the 32-bit wrap unobservable: run thousands of
  // folds across several wraps; every crossing lands on its boundary column
  // and the epoch stays an exact integer multiple ahead.
  {
    Flywheel f(cfg);
    uint32_t t = 0xFFFFFFFFu - kPeriod / 3; // wrap almost immediately
    f.seed(t);
    bool ok = true;
    for (int k = 1; k <= 5000; ++k) { // ~5.2 minutes of mock time, 43 wraps
      t += kPeriod;
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
    HS_EXPECT_TRUE(f.snap(Boundary::HALF, 1000000u + kPeriod - 2 * kCol,
                          &err) == Flywheel::SnapOutcome::ACCEPTED);
    HS_EXPECT_EQ(err, 2);
    HS_EXPECT_EQ(f.position(1000000u + kPeriod - 2 * kCol), 144);

    // Forged ZERO at the HALF position: W/2 correction → reject ×R → ACQUIRE.
    uint32_t t = 1000000u + kPeriod - 2 * kCol;
    Flywheel::SnapOutcome last = Flywheel::SnapOutcome::ACCEPTED;
    for (int i = 0; i < cfg.reject_fallback; ++i) {
      t += 10 * kCol;
      last = f.snap(Boundary::ZERO, t, &err);
    }
    HS_EXPECT_TRUE(last == Flywheel::SnapOutcome::REJECTED_FELL_BACK);
    HS_EXPECT_TRUE(f.lock() == LockState::ACQUIRE);
    // ACQUIRE: hard snap, relocks.
    t += 10 * kCol;
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
    f.snap(Boundary::HALF, 1000000u + 72 * kCol, &err); // W/4 off
    f.force_lock();
    // Real boundary stream: ZERO at k·rev, HALF at k·rev + half.
    uint32_t t = 1000000u + kPeriod; // true HALF instant
    Boundary b = Boundary::HALF;
    int accepted_at = -1;
    for (int i = 0; i < cfg.reject_fallback + 1; ++i) {
      if (f.snap(b, t, &err) == Flywheel::SnapOutcome::ACCEPTED) {
        accepted_at = i;
        break;
      }
      t += kPeriod;
      b = opposite(b);
    }
    HS_EXPECT_EQ(accepted_at, cfg.reject_fallback); // R rejections, then snap
    HS_EXPECT_TRUE(f.lock() == LockState::LOCKED);
    HS_EXPECT_EQ(f.position(t), boundary_column(b, 288));
  }
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
    HS_EXPECT_TRUE(e.schedule_boundary(Symbol::ZERO, b, b + kCol / 8, cfg));
    std::vector<uint32_t> pulses;
    for (uint32_t t = b + kCol / 8; t < b + 8 * kCol; t += kCol / 8) {
      if (e.tick(t, cfg, &aborted))
        pulses.push_back(t);
      HS_EXPECT_FALSE(aborted);
    }
    HS_EXPECT_EQ(pulses.size(), static_cast<size_t>(3));
    if (pulses.size() == 3) {
      // Each pulse within an oversample step of its scheduled slot.
      HS_EXPECT_LE(pulses[0] - b, kCol / 8);
      HS_EXPECT_LE(pulses[1] - (b + 2 * kCol), kCol / 8);
      HS_EXPECT_LE(pulses[2] - (b + 4 * kCol), kCol / 8);
    }
  }

  // Late at the boundary (> ~½ column): the whole symbol is censored.
  {
    SymbolEmitter e;
    HS_EXPECT_FALSE(e.schedule_boundary(Symbol::ZERO, 1000u,
                                        1000u + cfg.late_censor_cycles() + 1,
                                        cfg));
  }

  // Boundary scheduled in the FUTURE (now before at_cycles): not late — it must
  // be accepted and then emitted once `now` reaches the boundary. Regression
  // guard for the signed late-censor: an unsigned `now - at_cycles` wraps a
  // future boundary to a huge positive lateness and would wrongly censor it.
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
    HS_EXPECT_FALSE(e.tick(b + 2 * kCol + cfg.late_censor_cycles() + 1, cfg,
                           &aborted));
    HS_EXPECT_TRUE(aborted);
    HS_EXPECT_TRUE(e.idle());
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
    for (uint32_t t = t0; t < t0 + 100 * kCol; t += kCol / 8) {
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

// ── Multi-board simulator ───────────────────────────────────────────────────

/**
 * @brief One simulated board: its SyncBoard engine plus the host-side state the
 *        simulator models around it.
 * @details The modeled state covers crystal offset/phase, the masked-IRQ latch,
 *          symbol-drop and EMI windows, the foreground build/commit model, and
 *          probes.
 */
struct SimBoard {
  SyncBoard board;
  bool master = false;
  int32_t ppm = 0;       /**< Crystal offset, parts per million. */
  uint64_t phase0 = 0;   /**< Local cycle-counter offset at g = 0. */
  double next_tick = 0;  /**< Next flywheel wake, in global cycles. */
  double tick_step = 0;  /**< Wake period in global cycles. */
  /** Masked-IRQ model: while g < mask_until, wakes coalesce and edges latch
      into a single delayed delivery (one latched flag per pin). */
  uint64_t mask_until = 0;
  bool edge_latched = false;
  std::vector<std::pair<uint64_t, uint64_t>> masks; /**< [from, to) global. */
  uint64_t drop_from = 0, drop_to = 0;              /**< Symbol drop window. */
  // Foreground model.
  uint32_t seen_gen = 0;
  int32_t pending_index = -1;
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
  std::vector<SimBoard> boards;
  uint64_t g = 0;                                  /**< Global time, cycles. */
  std::vector<std::pair<uint64_t, int>> emi;       /**< (g, target), sorted. */
  size_t emi_pos = 0;

  /**
   * @brief Builds @p n boards (index 0 is the master) and seeds each flywheel.
   * @param c Shared sync configuration.
   * @param n Number of boards to construct.
   * @param ppm Per-board crystal offsets, parts per million (length @p n).
   * @param phase0 Common starting clock offset, cycles, at global time 0.
   * @details Each flywheel polls on a ⅛-column grid scaled by its ppm, with a
   *          small per-board boot stagger.
   */
  Sim(const Config &c, int n, const int32_t *ppm, uint64_t phase0 = 0)
      : cfg(c) {
    const double step0 = double(c.cycles_per_half_rev) / (c.W / 2) / 8.0;
    boards.reserve(n);
    for (int i = 0; i < n; ++i) {
      boards.emplace_back(c);
      SimBoard &b = boards.back();
      b.master = (i == 0);
      b.ppm = ppm[i];
      b.phase0 = phase0;
      b.tick_step = step0 * 1e6 / (1e6 + ppm[i]);
      b.next_tick = double(7 * (i + 1)); // small boot stagger
      b.board.seed(local_now(b, 0), b.master);
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
  template <typename Pred> bool run_until(Pred pred, double max_revs) {
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
    // into this one delayed wake). do-while: tg was rounded down from
    // next_tick, so the loop must advance at least once.
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
        b.board.mailbox().try_claim(now, b.board.gap_timeout_cycles(), &s))
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
      b.live = false; // release + delete the outgoing instance
      b.pending_index = SyncBoard::build_index_of(bw);
      b.pending_gen = gen;
      b.pending_ready_g = tg + b.init_delay;
      b.have_pending = true;
    }
    const bool ready =
        b.have_pending && b.pending_gen == gen && tg >= b.pending_ready_g;
    if (a.commit) {
      if (!ready) {
        b.trapped = true; // device: HS_CHECK fires
      } else {
        b.live = true;
        b.live_index = b.pending_index;
        b.t = 0;
        b.swap_g = tg;
      }
    } else if (a.join_boundary && !a.dark && !b.live && ready) {
      b.live = true;
      b.live_index = b.pending_index;
      b.t = 0;
      b.swap_g = tg;
    }
    if (a.flip) {
      ++b.flips;
      if (b.live)
        ++b.t;
    }
    b.dark_now = a.dark || !b.live;
  }
};

// ── Scenario: clean 4-board run (boot join, phase, flips, wrap) ─────────────

/**
 * @brief Verifies a clean 4-board run: every board locks within a revolution,
 *        joins live at the same boundary with the same effect and frame
 *        counter, and holds sub-2-column phase, equal frame counters, and ~2
 *        flips/rev through a 32-bit clock wrap — with no gate rejections,
 *        invalid symbols, or traps.
 */
inline void test_sim_boot_and_phase() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 20, -20, 40};
  // Local clocks start just below the 32-bit wrap: every board's CYCCNT
  // wraps ~10 revolutions in, mid-run (§12 timebase arithmetic).
  Sim sim(cfg, 4, ppm, 0xFFFFFFFFull - 10ull * 2 * kPeriod + 12345);

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
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (!b.live)
            return false;
        return true;
      },
      double(cfg.join_grid_revs) + 2));
  for (int i = 0; i < 4; ++i) {
    HS_EXPECT_EQ(sim.boards[i].live_index, 0);
    // Swaps happen at each board's own crossing of the same boundary —
    // within a couple of columns of global time of each other.
    const int64_t dg = static_cast<int64_t>(sim.boards[i].swap_g) -
                       static_cast<int64_t>(sim.boards[0].swap_g);
    HS_EXPECT_LE(dg < 0 ? -dg : dg, int64_t(3) * kCol);
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
    sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1);
    HS_EXPECT_LE(sim.max_phase_err(), 2);
    for (int i = 0; i < 4; ++i) {
      const uint64_t df = sim.boards[i].flips - flips_before[i];
      HS_EXPECT_GE(df, 8u); // ~2 flips/rev over the ≥4-rev slice
      HS_EXPECT_LE(df, 12u);
      HS_EXPECT_EQ(sim.boards[i].live_index, sim.boards[0].live_index);
      HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
    }
  }
  // No gate rejections, no invalid symbols, no traps in a clean run.
  for (int i = 1; i < 4; ++i) {
    const Telemetry &tm = sim.boards[i].board.telemetry();
    HS_EXPECT_EQ(tm.symbols_rejected_gate, 0u);
    HS_EXPECT_EQ(tm.symbols_discarded_invalid, 0u);
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

  // Through boot join.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (!b.live)
            return false;
        return true;
      },
      double(cfg.join_grid_revs) + 2));

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
        return s.boards[0].board.content().commit_in_revs <=
               s.cfg.commit_revs;
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
    HS_EXPECT_LE(dg < 0 ? -dg : dg, int64_t(3) * kCol);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }
  // Post-epoch: full content coherence (index AND t) including the master.
  sim.run_revs(3.0);
  sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1);
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
      static_cast<uint64_t>(3) * 2 * kPeriod; // 3 revs > K = 2
  // Boot joins are NOT deadline-bound: board 2 simply goes live ~3 revs
  // after the others (next join-grid boundary), without trapping.
  sim.run_revs(double(cfg.join_grid_revs) * 3);
  HS_EXPECT_TRUE(sim.boards[2].live);
  HS_EXPECT_FALSE(sim.boards[2].trapped);
  // The epoch commit IS deadline-bound.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return s.boards[2].trapped; },
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
  const uint64_t rev = 2ull * kPeriod;
  for (int k = 8; k < 28; ++k) {
    const uint64_t b0 = k * rev; // master ZERO crossings ≈ k·rev (ppm 0)
    sim.boards[2].masks.push_back(
        {b0 - kCol / 2, b0 + 2 * kCol + kCol / 2});
  }
  // Board 3: mid-revolution masks (no boundary, no symbol) — pure wake
  // coalescing; the flywheel resumes at the time-correct column.
  for (int k = 8; k < 28; ++k) {
    const uint64_t m0 = k * rev + 40 * kCol;
    sim.boards[3].masks.push_back({m0, m0 + 5 * kCol});
  }

  sim.run_revs(30.0);
  sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1);

  // Truncated bursts were discarded (count telemetry), never accepted as
  // the wrong boundary: phase stays sub-column-ish and content equal.
  const Telemetry &tm2 = sim.boards[2].board.telemetry();
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
  sim2.boards[0].masks.push_back({b0 - kCol / 4, b0 + 2 * kCol});
  sim2.run_revs(6.0);
  const Telemetry &tm0 = sim2.boards[0].board.telemetry();
  HS_EXPECT_GT(tm0.emit_censored + tm0.emit_aborted, 0u);
  HS_EXPECT_LE(sim2.max_phase_err(), 2);
  HS_EXPECT_EQ(sim2.boards[1].board.telemetry().symbols_rejected_gate, 0u);
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
  const uint64_t rev = 2ull * kPeriod;

  // Isolated spurious edges on board 1, ~1 per revolution at varied
  // mid-revolution offsets (away from boundaries): each forms a 1-pulse
  // burst = a valid HALF symbol, which the LOCKED plausibility gate must
  // reject (implied correction ≫ G).
  uint32_t lcg = 12345;
  for (int k = 8; k < 40; ++k) {
    lcg = lcg * 1664525u + 1013904223u;
    const uint64_t off = (20 + lcg % 100) * static_cast<uint64_t>(kCol);
    sim.emi.push_back({k * rev + off, 1});
  }
  // Two edges injected INSIDE master ZERO bursts on board 2 (count 3 → 4):
  // even count = invalid, discarded whole; the crossing flip covers it.
  sim.emi.push_back({10 * rev + kCol, 2});
  sim.emi.push_back({14 * rev + kCol, 2});
  // One edge within G of board 3's predicted HALF boundary: the §9.1
  // accepted-EMI case — a ≤G-column seam for ≤½ rev, re-snapped by the next
  // real symbol. It must not unlock or misclassify anything.
  sim.emi.push_back({12 * rev + (kPeriod - 2ull * kCol), 3});
  std::sort(sim.emi.begin(), sim.emi.end());

  sim.run_revs(42.0);
  sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1);

  const Telemetry &tm1 = sim.boards[1].board.telemetry();
  HS_EXPECT_GT(tm1.symbols_rejected_gate, 20u); // isolated EMI all rejected
  const Telemetry &tm2 = sim.boards[2].board.telemetry();
  HS_EXPECT_GE(tm2.symbols_discarded_invalid, 2u); // corrupted bursts dropped
  // System health: everyone locked, coherent, phase bounded.
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
  const uint64_t rev = 2ull * kPeriod;

  // Boot join first.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (!b.live)
            return false;
        return true;
      },
      double(cfg.join_grid_revs) + 2));

  // Plain symbol drop: board 1 hears nothing for 2 revolutions — it coasts
  // on its crystal (telemetry max_coast) and silently re-snaps after.
  sim.boards[1].drop_from = sim.g + rev;
  sim.boards[1].drop_to = sim.g + 3 * rev;
  sim.run_revs(5.0);
  HS_EXPECT_GE(sim.boards[1].board.telemetry().max_coast_halves, 4u);
  HS_EXPECT_TRUE(sim.boards[1].board.lock() == LockState::LOCKED);
  sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1);
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
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return s.boards[3].live_index == 1; },
      double(cfg.beacon_period_revs + cfg.join_grid_revs) + 6));
  HS_EXPECT_EQ(sim.boards[3].board.telemetry().beacon_index_corrections, 1u);
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
  // Identity from the next beacon, display from the next join-grid
  // boundary — comfortably inside the ~2 s budget at ship cadence; never a
  // wrong frame in between (dark throughout).
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return s.boards[2].live; },
      double(cfg.beacon_period_revs + cfg.join_grid_revs) + 4));
  HS_EXPECT_EQ(b2.live_index, sim.boards[0].live_index);
  sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1);
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
  const uint64_t rev = 2ull * kPeriod;

  // A spurious flip needs a forged burst that is simultaneously valid,
  // plausible, and boundary-consistent (§5.3). Forge the strongest cheap
  // attack: an isolated valid-count (HALF) burst 30 columns past a real
  // ZERO boundary — clear of the real burst's gap-timeout window and of the
  // quiet-before guard, far from both predicted boundaries. It must be held
  // as a suspect and counted as a rejection, never snapped or flipped.
  const uint64_t b0 = (sim.g / rev + 2) * rev; // a future master ZERO
  sim.emi.push_back({b0 + 30 * static_cast<uint64_t>(kCol), 1});
  sim.emi_pos = 0;
  std::sort(sim.emi.begin(), sim.emi.end());

  const uint64_t flips_before = sim.boards[1].flips;
  const uint32_t rejected_before =
      sim.boards[1].board.telemetry().symbols_rejected_gate;
  sim.run_revs(4.0);
  HS_EXPECT_GT(sim.boards[1].board.telemetry().symbols_rejected_gate,
               rejected_before);
  // Flip cadence unbroken: 2 per revolution, no extra content advance.
  const uint64_t df = sim.boards[1].flips - flips_before;
  HS_EXPECT_GE(df, 7u);
  HS_EXPECT_LE(df, 9u);
  sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1);
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
  const uint64_t rev = 2ull * kPeriod;

  /**
   * @brief Runs the sim until every board has gone live (boot join complete).
   * @param sim The simulation to advance.
   * @return True if all boards went live within the join budget.
   */
  auto boot_join = [&cfg](Sim &sim) {
    return sim.run_until(
        [](Sim &s) {
          for (auto &b : s.boards)
            if (!b.live)
              return false;
          return true;
        },
        double(cfg.join_grid_revs) + 2);
  };
  /**
   * @brief Runs the sim to one revolution before the train start.
   * @param sim The simulation to advance.
   * @return True if the pre-train point was reached within budget.
   * @details The master (ppm 0) crosses ZERO on exact rev multiples, so the
   *          primary copy rides the crossing ≈ one rev from the moment this
   *          predicate fires.
   */
  auto to_pre_train = [&cfg](Sim &sim) {
    return sim.run_until(
        [](Sim &s) {
          return s.boards[0].board.content().rev_in_effect >=
                 s.cfg.revs_per_effect - 1;
        },
        double(cfg.revs_per_effect) + 4);
  };
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
      HS_EXPECT_LE(dg < 0 ? -dg : dg, int64_t(3) * kCol);
      HS_EXPECT_FALSE(sim.boards[i].trapped);
    }
    sim.run_revs(2.0);
    sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1);
    for (int i = 1; i < 4; ++i)
      HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
  };
  const double commit_revs_max =
      double(cfg.commit_revs + static_cast<uint32_t>(cfg.epoch_repeats)) + 8;

  // Sub-scenario 1: board 2 misses ONLY the primary copy at B.
  {
    const int32_t ppm[4] = {0, 20, -20, 10};
    Sim sim(cfg, 4, ppm);
    HS_EXPECT_TRUE(boot_join(sim));
    HS_EXPECT_TRUE(to_pre_train(sim));
    sim.boards[2].drop_from = sim.g + rev - 4 * kCol;
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
    HS_EXPECT_TRUE(boot_join(sim));
    HS_EXPECT_TRUE(to_pre_train(sim));
    const uint64_t b0 = sim.g + rev;
    sim.boards[0].masks.push_back({b0 - kCol / 4, b0 + 2 * kCol});
    HS_EXPECT_TRUE(sim.run_until(all_on_effect_1, commit_revs_max));
    const Telemetry &tm0 = sim.boards[0].board.telemetry();
    HS_EXPECT_GT(tm0.emit_censored + tm0.emit_aborted, 0u);
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
  const uint64_t rev = 2ull * kPeriod;
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (!b.live)
            return false;
        return true;
      },
      double(cfg.join_grid_revs) + 2));
  // Settle past the boot beacons (revs 1–3), then slip board 3's counter
  // by +2: at the next train it would over-count j by 2 and commit 2
  // revolutions EARLY.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return s.boards[0].board.content().rev_in_effect == 5; },
      12.0));
  sim.boards[3].board.content_mut().rev_in_effect += 2;

  // Detected and corrected at the next beacon (rev 9), within one period.
  // (Pre-resync the counters differ by 2 — a crossing straddle changes the
  // gap by at most 1, so this predicate cannot fire spuriously.)
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return s.boards[3].board.content().rev_in_effect ==
               s.boards[0].board.content().rev_in_effect;
      },
      double(cfg.beacon_period_revs) + 2));
  HS_EXPECT_GE(sim.boards[3].board.telemetry().beacon_rev_mismatches, 1u);

  // The next epoch commits in lockstep even though board 3 ALSO misses the
  // primary copy, taking the repeat path that depends on j-inference.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return s.boards[0].board.content().rev_in_effect >=
               s.cfg.revs_per_effect - 1;
      },
      double(cfg.revs_per_effect) + 4));
  sim.boards[3].drop_from = sim.g + rev - 4 * kCol;
  sim.boards[3].drop_to = sim.g + rev + rev / 4;
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 1)
            return false;
        return true;
      },
      double(cfg.commit_revs + static_cast<uint32_t>(cfg.epoch_repeats)) +
          8));
  for (int i = 1; i < 4; ++i) {
    const int64_t dg = static_cast<int64_t>(sim.boards[i].swap_g) -
                       static_cast<int64_t>(sim.boards[0].swap_g);
    HS_EXPECT_LE(dg < 0 ? -dg : dg, int64_t(3) * kCol);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }
  sim.run_revs(2.0);
  sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1);
  for (int i = 1; i < 4; ++i)
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
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
      c.rev_in_effect = RPE + j - 1;       // B+j fold still pending…
      HS_EXPECT_FALSE(c.on_zero_crossing(cfg)); // …backstop apply_flip folds it
    } else {
      c.rev_in_effect = RPE + j;           // already folded a tick earlier
    }
    const uint32_t base_rev = c.rev_in_effect; // == RPE + j either way
    HS_EXPECT_EQ(base_rev, RPE + j);
    HS_EXPECT_TRUE(c.on_epoch_symbol(cfg));  // opens the commit window
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
//   dropped render                    → not simulable here: frame pacing
//                                       lives in the device's Canvas
//                                       buffer_free() gate; the epoch-bounded
//                                       t recovery it relies on is pinned by
//                                       test_sim_epoch_commit
//   missed epoch (all copies)         → test_sim_drops_and_missed_epoch
//   epoch repeat lockstep (§6.3.1)    → test_sim_epoch_repeat_lockstep
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
  sim.run_until([](Sim &s) { return s.board_pos(0) == 40; }, 1.1);

  // Deafen board 1 for exactly one half-rev, aligned mid-half: it misses
  // exactly one boundary symbol (the HALF, 104 columns ahead).
  sim.boards[1].drop_from = sim.g;
  sim.boards[1].drop_to = sim.g + kPeriod;
  const Telemetry &tm = sim.boards[1].board.telemetry();
  const uint32_t acc_before = tm.symbols_accepted;

  HS_EXPECT_LE(max_err_over(sim, 1.5), 1); // sub-column through coast+re-snap
  HS_EXPECT_TRUE(sim.boards[1].board.lock() == LockState::LOCKED);
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
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (!b.live)
            return false;
        return true;
      },
      double(cfg.join_grid_revs) + 2));
  sim.run_revs(2.0);
  sim.run_until([](Sim &s) { return s.board_pos(0) == 40; }, 1.1);

  // The master HALF boundary is 104 columns ahead. Censor the real symbol
  // for board 1 and forge an edge 3 columns early: isolated, valid count,
  // within G of the predicted boundary — the §9.1 accepted case.
  const uint64_t h = sim.g + 104ull * kCol;
  sim.boards[1].drop_from = h - kCol;
  sim.boards[1].drop_to = h + 8 * kCol;
  sim.emi.push_back({h - 3 * kCol, 1});
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
  sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1);
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
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (!b.live)
            return false;
        return true;
      },
      double(cfg.join_grid_revs) + 2));
  // Park at rev 10 (≡ 2 mod the beacon period) so the ACQUIRE window stays
  // clear of beacon trains: an isolated first digit could re-poison an
  // ACQUIRE board (the §9.1 forged-during-ACQUIRE sub-case) and double the
  // recovery — real but timing-dependent; this test pins the common path.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return s.boards[0].board.content().rev_in_effect == 10; },
      16.0));

  // Corrupt board 2's flywheel phase by W/4 — far beyond the gate.
  SimBoard &b2 = sim.boards[2];
  const int32_t bogus = floor_mod(sim.board_pos(2) + 72, cfg.W);
  b2.board.flywheel_mut().seed(Sim::local_now(b2, sim.g) -
                               static_cast<uint32_t>(bogus) * kCol);
  b2.board.flywheel_mut().force_lock();
  HS_EXPECT_GE(circ_dist(sim.board_pos(2), sim.board_pos(0), cfg.W), 60);

  const uint32_t rej_before = b2.board.telemetry().symbols_rejected_gate;
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return s.boards[2].board.lock() == LockState::LOCKED &&
               circ_dist(s.board_pos(2), s.board_pos(0), s.cfg.W) <= 1;
      },
      3.5)); // budget: ~2.8 revs; slack for the mid-rev corruption instant
  HS_EXPECT_GE(b2.board.telemetry().symbols_rejected_gate - rej_before,
               static_cast<uint32_t>(cfg.reject_fallback));
  HS_EXPECT_GE(b2.board.telemetry().lock_transitions, 2u);

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
    HS_EXPECT_LE(dg < 0 ? -dg : dg, int64_t(3) * kCol);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }
  sim.run_revs(2.0);
  sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1);
  for (int i = 1; i < 4; ++i)
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
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
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (!b.live)
            return false;
        return true;
      },
      double(cfg.join_grid_revs) + 2));
  // Master beacons ride rev ≡ 1 (mod 8); park at the rev-9 ZERO crossing.
  // The train starts when the master reaches x = W/4 = 72.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return s.boards[0].board.content().rev_in_effect == 9; },
      16.0));
  // Beacon (index 0, rev 9) is digits [0,0,1,1,2]; its 4th burst is two
  // pulses at relative columns 16 and 17. One EMI edge between them on
  // board 1 (clear of the 100 µs glitch filter) makes that digit read 2:
  // checksum mismatch, whole frame dropped.
  const uint64_t train = sim.g + 72ull * kCol;
  sim.emi.push_back({train + 16ull * kCol + kCol / 2, 1});
  sim.emi_pos = 0;
  std::sort(sim.emi.begin(), sim.emi.end());

  const Telemetry &tm1 = sim.boards[1].board.telemetry();
  const uint32_t ok_before = tm1.beacons_ok;
  const uint32_t rej_before = tm1.beacons_rejected;
  sim.run_revs(1.0);
  HS_EXPECT_EQ(tm1.beacons_rejected, rej_before + 1); // dropped whole
  HS_EXPECT_EQ(tm1.beacons_ok, ok_before);            // nothing applied
  HS_EXPECT_EQ(tm1.beacon_index_corrections, 0u);
  // Recovery: the next clean beacon decodes within one period.
  HS_EXPECT_TRUE(sim.run_until(
      [ok_before](Sim &s) {
        return s.boards[1].board.telemetry().beacons_ok > ok_before;
      },
      double(cfg.beacon_period_revs) + 1));
  // The show never noticed: index, phase, and frame counters all clean.
  HS_EXPECT_EQ(sim.boards[1].live_index, sim.boards[0].live_index);
  sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1);
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
  Sim sim(cfg, 3, ppm, 0xFFFFFFFFull - 30ull * 2 * kPeriod + 999);
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (!b.live)
            return false;
        return true;
      },
      double(cfg.join_grid_revs) + 2));
  sim.run_revs(2.0);
  // Cut the wire at a quiet point (mid-half, no beacon this rev) so no
  // burst is in flight: a half-received burst would register one truncated-
  // count artifact, which is the lost-symbol row's case, not this one.
  sim.run_until([](Sim &s) { return s.board_pos(0) == 40; }, 1.1);

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
    HS_EXPECT_EQ(sim.boards[i].board.telemetry().symbols_rejected_gate, 0u);
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
  // col; 25 ppm ≈ 1.1 col. And these positions are only meaningful because
  // the rebase rule survived the snap-free wrap.
  const int32_t e1 = circ_dist(sim.board_pos(1), sim.board_pos(0), cfg.W);
  const int32_t e2 = circ_dist(sim.board_pos(2), sim.board_pos(0), cfg.W);
  HS_EXPECT_GE(e1, 1);
  HS_EXPECT_LE(e1, 3);
  HS_EXPECT_GE(e2, 1);
  HS_EXPECT_LE(e2, 2);
  HS_EXPECT_GE(sim.boards[1].board.telemetry().max_coast_halves, 250u);
}

// ── Runner ──────────────────────────────────────────────────────────────────

/**
 * @brief Runs every pov_sync test in order.
 * @return The module's failure count.
 */
inline int run_pov_sync_tests() {
  auto scope = begin_module("pov_sync");

  test_helpers();
  test_alphabet();
  test_flip_gate();
  test_mailbox();
  test_mailbox_prior_staleness();
  test_seed_clears_mailbox();
  test_beacon_codec();
  test_flywheel_position();
  test_snap_gate();
  test_emitter();

  test_sim_boot_and_phase();
  test_sim_epoch_commit();
  test_sim_commit_deadline_trap();
  test_sim_masked_windows();
  test_sim_emi();
  test_sim_drops_and_missed_epoch();
  test_sim_reboot();
  test_sim_forged_burst();
  test_sim_epoch_repeat_lockstep();
  test_sim_rev_resync();
  test_epoch_same_tick_burst_fold();

  test_budget_lost_symbol();
  test_budget_emi_accepted_seam();
  test_budget_corrupted_timebase();
  test_budget_beacon_corruption();
  test_budget_wire_dead();

  return end_module(scope);
}

} // namespace pov_sync_tests
} // namespace hs_test
