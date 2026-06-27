/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file pov_sync.h
 * @brief Pure, host-testable core of the Phantasm synchronization design
 *        (docs/phantasm_frame_sync_spec.md): one local flywheel timebase per
 *        board, disciplined over a single sync-symbol wire.
 *
 * Split out of pov_segmented.h (which is Arduino-only) so every load-bearing
 * decision — position math, symbol classification, the acceptance gate, epoch
 * scheduling, beacon framing, emission self-censoring — is unit-testable on
 * the host, exactly as pov_segment_map.h is for the index math. The device
 * driver is a thin shell over SyncBoard: it reads the cycle counter, services
 * two ISRs, packs pixels, and toggles one pin.
 *
 * Architecture (spec §3): every board derives its column position from a
 * free-running hardware cycle counter (`x = f(now - epoch)`), never from
 * counting timer interrupts, so masked-IRQ windows cannot drop columns. The
 * master emits count-coded symbol bursts on the one wire — 2/revolution
 * boundary marks plus a rare epoch mark and a mid-revolution data beacon —
 * and downstream boards snap their flywheel phase to them. Three layers ride
 * the same timebase:
 *
 *   Layer 1  column position    — flywheel, snapped by boundary symbols
 *   Layer 2  buffer flip        — local boundary crossing, symbol backstop,
 *                                 exactly-once via boundary-identity dedup
 *   Layer 3  content (t/effect) — epoch symbol (deadline-scheduled commit)
 *                                 + absolute-index beacon
 *
 * Concurrency contract (spec §8): the sync-wire edge ISR is a pure publisher
 * into EdgeMailbox; the flywheel ISR is the single owner of all other state
 * and calls SyncBoard::tick(). Nothing here allocates, blocks, or does I/O.
 */
#pragma once

#include <cstdint>

#include "core/platform.h" // HS_CHECK used in Flywheel's constructor guard

// Forward declaration of the unit-test accessor that reaches EdgeMailbox's
// private split consumer path (burst_complete()/claim()); see EdgeMailbox.
namespace hs_test {
namespace pov_sync_tests {
struct EdgeMailboxTestAccess;
} // namespace pov_sync_tests
} // namespace hs_test

namespace pov {
namespace sync {

// ── Small integer helpers (signed-correct division/modulo) ─────────────────

/**
 * @brief Floor division: rounds toward -inf, unlike C++ truncation.
 * @param a Dividend.
 * @param b Divisor (nonzero).
 * @return floor(a / b).
 */
constexpr int64_t floor_div(int64_t a, int64_t b) {
  const int64_t q = a / b;
  const int64_t r = a % b;
  return (r != 0 && ((r < 0) != (b < 0))) ? q - 1 : q;
}

/**
 * @brief Non-negative modulo.
 * @param a Value to reduce.
 * @param m Modulus (positive).
 * @return a mod m in [0, m).
 */
constexpr int32_t floor_mod(int64_t a, int32_t m) {
  const int32_t r = static_cast<int32_t>(a % m);
  return r < 0 ? r + m : r;
}

/**
 * @brief Circular distance between two columns on a ring.
 * @param a First column index.
 * @param b Second column index.
 * @param w Ring width (columns per revolution).
 * @return Shortest distance around the ring, in columns, in [0, w/2].
 */
constexpr int32_t circ_dist(int32_t a, int32_t b, int32_t w) {
  // (a-b)%w lands in (-w, w) for any a, b, so one conditional add normalizes.
  int32_t d = (a - b) % w;
  if (d < 0)
    d += w;
  return d > w / 2 ? w - d : d;
}

// ── Configuration ───────────────────────────────────────────────────────────

/**
 * @brief All protocol constants, in columns and cycle-counter cycles.
 *
 * "Cycles" are ticks of the per-board free-running clock (DWT->CYCCNT on the
 * device, a mock counter in tests). Timestamps are uint32_t and wrap; all
 * arithmetic on them is modular, and the flywheel's rebase rule (spec §4.1)
 * keeps every difference far below the wrap period.
 */
struct Config {
  int32_t W = 288;                         /**< Columns per revolution. */
  uint32_t cycles_per_half_rev = 37500000; /**< Timebase constant (spec §4.1). */
  uint32_t glitch_filter_cycles = 60000;   /**< Min edge spacing (~100 µs). */

  // Symbol wire (spec §5.2): pitches/timeouts in columns.
  int32_t pulse_pitch_cols = 2;  /**< Boundary-burst pulse pitch (> mask M). */
  int32_t beacon_pitch_cols = 1; /**< Beacon digit pulse pitch (checksummed). */
  int32_t gap_timeout_cols = 4;  /**< Quiet time that terminates a burst. */

  // Acceptance gate (spec §5.3).
  int32_t gate_cols = 4;       /**< G: max plausible snap correction (LOCKED). */
  int32_t reject_fallback = 4; /**< R: consecutive rejections before ACQUIRE. */
  int32_t acquire_quiet_cols = 16; /**< Quiet-before guard for ACQUIRE snaps. */

  // Content layer (spec §6).
  uint32_t revs_per_effect = 960; /**< Effect duration in revolutions (120 s). */
  int32_t epoch_repeats = 3;      /**< EPOCH redundancy repeats (spec §6.3). */
  uint32_t refractory_revs = 16;  /**< EPOCH dedup window (spec §6.1). */
  /**
   * @brief K: the construction window, the last K revolutions before the
   * commit.
   * @details The commit boundary itself is B + epoch_repeats + commit_revs from
   * the train's primary copy at B, so a board hearing any repeat still commits
   * in lockstep with the full K-rev construction budget (spec §6.1, §6.3.1).
   */
  uint32_t commit_revs = 2;
  uint32_t beacon_period_revs = 16;            /**< Beacon cadence (spec §6.4). */
  int32_t beacon_interdigit_timeout_cols = 24; /**< Stale-frame reset bound. */
  int32_t effect_count = 1; /**< Roster length (epoch wraps mod this). */
  /**
   * @brief Live-takeover grid: boards take a constructed effect live only at
   * revolutions ≡ 0 mod this.
   * @details At boot every board (master included) therefore goes live at the
   * SAME crossing with frame counters aligned, instead of master leading by
   * however long downstream identity took. Must divide 64 so a beacon's
   * mod-64 revolution count lands on the same grid as the master's true
   * count. Mid-show rejoins wait ≤ grid revolutions (well inside the ~2 s
   * rejoin budget, spec §9.1).
   */
  uint32_t join_grid_revs = 4;

  /**
   * @brief Rejoin budget: the most revolutions a mid-show board may wait to
   * hear its identity beacon (spec §9.1).
   * @details A board that joins just after a beacon waits up to
   * beacon_period_revs revs for the next one, so the beacon cadence must not
   * exceed this budget — enforced in valid() rather than left to prose.
   * Expressed in revolutions so the bound is rotation-rate independent; the
   * spec's "~2 s" reading is this budget at the nominal 480 RPM (16 revs).
   */
  uint32_t rejoin_budget_revs = 16;

  /**
   * @brief Cycles per column of rotation.
   * @return Cycle-counter cycles spanning one column.
   * @details Truncates (~2.6 ppm); used only for pitches/thresholds — position
   * math divides by the exact cycles_per_half_rev (spec §4.1, "stored period is
   * per-half-rev").
   */
  constexpr uint32_t cycles_per_column() const {
    return cycles_per_half_rev / static_cast<uint32_t>(W / 2);
  }
  /**
   * @brief Cycles spanning a span of columns.
   * @param cols Column count.
   * @return cols · cycles_per_column(), in cycles.
   */
  constexpr uint32_t col_cycles(int32_t cols) const {
    return cycles_per_column() * static_cast<uint32_t>(cols);
  }
  /**
   * @brief Burst-terminating gap, in cycles.
   * @return Quiet time that terminates a burst, in cycles.
   */
  constexpr uint32_t gap_timeout_cycles() const {
    return col_cycles(gap_timeout_cols);
  }
  /**
   * @brief Boundary-burst pulse pitch, in cycles.
   * @return Spacing between boundary-burst pulses, in cycles.
   */
  constexpr uint32_t pulse_pitch_cycles() const {
    return col_cycles(pulse_pitch_cols);
  }
  /**
   * @brief Beacon-digit pulse pitch, in cycles.
   * @return Spacing between beacon-digit pulses, in cycles.
   */
  constexpr uint32_t beacon_pitch_cycles() const {
    return col_cycles(beacon_pitch_cols);
  }
  /**
   * @brief Emission lateness budget (spec §5.2: self-censor past ~½ column).
   * @return Half a column, in cycles.
   */
  constexpr uint32_t late_censor_cycles() const {
    return cycles_per_column() / 2;
  }
  /**
   * @brief Quiet-before guard for ACQUIRE snaps, in cycles.
   * @return acquire_quiet_cols converted to cycles.
   */
  constexpr uint32_t acquire_quiet_cycles() const {
    return col_cycles(acquire_quiet_cols);
  }
  /**
   * @brief Beacon interdigit stale-frame reset bound, in cycles.
   * @return beacon_interdigit_timeout_cols converted to cycles.
   */
  constexpr uint32_t interdigit_timeout_cycles() const {
    return col_cycles(beacon_interdigit_timeout_cols);
  }

  /**
   * @brief Boot-time sanity check for the driver's HS_CHECK.
   * @return True if every protocol constant is self-consistent.
   * @details gate_cols < W/4 is what lets the gate's distance check subsume the
   * boundary-identity check; see Flywheel::snap.
   */
  constexpr bool valid() const {
    return W > 0 && W % 2 == 0 && cycles_per_half_rev > 0 && gate_cols > 0 &&
           gate_cols < W / 4 && reject_fallback > 0 && glitch_filter_cycles > 0 &&
           pulse_pitch_cols > 0 && gap_timeout_cols > pulse_pitch_cols &&
           beacon_pitch_cols > 0 && gap_timeout_cols > beacon_pitch_cols &&
           effect_count > 0 && effect_count <= 64 && commit_revs > 0 &&
           // Gate epoch_repeats >= 0 first: a negative value casts to a huge
           // uint32_t and wraps the refractory bound below.
           epoch_repeats >= 0 &&
           refractory_revs >
               commit_revs + static_cast<uint32_t>(epoch_repeats) &&
           revs_per_effect > refractory_revs &&
           beacon_period_revs > commit_revs &&
           // Beacon rev resync recovers a slip only in (-32, +32), so keep the
           // period below the half-window.
           beacon_period_revs < 32 &&
           // §9.1 rejoin budget: a joiner waits up to beacon_period_revs revs
           // for the next identity, so cap the cadence at the budget.
           beacon_period_revs <= rejoin_budget_revs && join_grid_revs > 0 &&
           (64u % join_grid_revs) == 0 &&
           // schedule_beacon's "is-due" check reads (now - start_cycles) as
           // int32, so the worst-case span (5 digits of value 7) must clear 2^31.
           5u * (7u * beacon_pitch_cycles() +
                 static_cast<uint32_t>(gap_timeout_cols + 1) *
                     cycles_per_column()) <
               static_cast<uint32_t>(INT32_MAX);
  }
};

/**
 * @brief Builds the Phantasm Config for given hardware parameters.
 * @param cpu_hz CPU clock frequency, in Hz.
 * @param rpm Spindle speed, in revolutions per minute.
 * @param w Columns per revolution.
 * @param effect_count Effect roster length.
 * @return A populated Config.
 * @details cycles_per_half_rev = cpu_hz · 30 / rpm (exact for 600 MHz /
 * 480 RPM → 37,500,000); glitch filter is 100 µs.
 */
constexpr Config phantasm_config(uint32_t cpu_hz, uint32_t rpm, int32_t w,
                                 int32_t effect_count) {
  Config c{};
  c.W = w;
  c.cycles_per_half_rev =
      static_cast<uint32_t>(static_cast<uint64_t>(cpu_hz) * 30u / rpm);
  c.glitch_filter_cycles = cpu_hz / 10000u;
  c.effect_count = effect_count;
  return c;
}

// ── Symbol alphabet (spec §5.2, §7) ─────────────────────────────────────────

/**
 * @brief The two half-rev boundaries (ZERO at column 0, HALF at W/2).
 * @details NONE means neither boundary.
 */
enum class Boundary : uint8_t { NONE, ZERO, HALF };

/**
 * @brief Decoded sync symbol.
 * @details ZERO and ZERO_EPOCH both mark the ZERO boundary; only ZERO_EPOCH
 * additionally carries an epoch (content-commit) mark.
 */
enum class Symbol : uint8_t { INVALID, HALF, ZERO, ZERO_EPOCH };

/**
 * @brief Count-coded classification: odd-only, distance-2 alphabet.
 * @param rising_edges Number of rising edges in the burst.
 * @return The decoded Symbol, or INVALID.
 * @details Any other count (a lost or spurious edge lands on an even value) is
 * INVALID and the whole burst is discarded — fail to "missed", never to
 * "wrong".
 */
constexpr Symbol classify_count(uint32_t rising_edges) {
  switch (rising_edges) {
  case 1:
    return Symbol::HALF;
  case 3:
    return Symbol::ZERO;
  case 5:
    return Symbol::ZERO_EPOCH;
  default:
    return Symbol::INVALID;
  }
}

/**
 * @brief Pulse count the master emits for a symbol (inverse of classify_count).
 * @param s Symbol to encode.
 * @return Number of pulses, or 0 for INVALID.
 */
constexpr uint32_t symbol_pulse_count(Symbol s) {
  switch (s) {
  case Symbol::HALF:
    return 1;
  case Symbol::ZERO:
    return 3;
  case Symbol::ZERO_EPOCH:
    return 5;
  default:
    return 0;
  }
}

/**
 * @brief Which boundary a symbol marks.
 * @param s Symbol to inspect.
 * @return The marked Boundary (both ZERO symbols mark ZERO); NONE for INVALID.
 */
constexpr Boundary symbol_boundary(Symbol s) {
  switch (s) {
  case Symbol::HALF:
    return Boundary::HALF;
  case Symbol::ZERO:
  case Symbol::ZERO_EPOCH:
    return Boundary::ZERO;
  default:
    return Boundary::NONE;
  }
}

/**
 * @brief The other half-rev boundary (boundaries strictly alternate ZERO/HALF).
 * @param b Current boundary.
 * @return HALF if b is ZERO, otherwise ZERO.
 */
constexpr Boundary opposite(Boundary b) {
  return b == Boundary::ZERO ? Boundary::HALF : Boundary::ZERO;
}

/**
 * @brief Column index of a boundary on a ring of given width.
 * @param b Boundary to locate.
 * @param w Ring width (columns per revolution).
 * @return Column index (ZERO→0, HALF→w/2).
 */
constexpr int32_t boundary_column(Boundary b, int32_t w) {
  return b == Boundary::HALF ? w / 2 : 0;
}

// ── Layer 2: exactly-once flip via boundary identity (spec §5.1) ────────────

/**
 * @brief Deduplicates the two flip paths (local flywheel crossing and sync
 * symbol) on the identity of the boundary. Boundaries strictly alternate
 * (ZERO, HALF, ZERO, …), so "same boundary as last flip" means the other path
 * already fired for it.
 */
struct FlipGate {
  Boundary last_flipped = Boundary::NONE; /**< Boundary of the last flip. */

  /**
   * @brief Attempts a flip for boundary @p b, deduped across both paths.
   * @param b Boundary being crossed.
   * @return True exactly once per distinct boundary, across both paths.
   */
  bool try_flip(Boundary b) {
    if (b == Boundary::NONE || b == last_flipped)
      return false;
    last_flipped = b;
    return true;
  }
};

// ── Sync-wire edge mailbox (spec §8: single-writer handoff) ─────────────────

/**
 * @brief Consumer's copy of a terminated burst.
 */
struct BurstSnapshot {
  uint32_t count = 0;        /**< Rising edges (after glitch filter). */
  uint32_t first_cycles = 0; /**< Timestamp of the burst's first edge. */
  uint32_t last_cycles = 0;  /**< Timestamp of the burst's last edge. */
};

/**
 * @brief The only state the sync-wire edge ISR writes. The publisher applies
 * the glitch filter and accumulates edge count + first/last timestamps; the
 * consumer (flywheel ISR) detects burst termination by gap timeout and claims
 * the burst with try_claim(), which tests completion and takes the burst in one
 * step so the two cannot be split around an edge. On the device that call runs
 * under a brief IRQ-off window (a two-instruction copy, spec §8.2); on the host
 * it is plain.
 */
class EdgeMailbox {
public:
  /**
   * @brief Publisher (edge ISR): record one rising edge.
   * @param now Edge timestamp, in cycles.
   * @param glitch_filter_cycles Minimum spacing below which an edge is rejected
   * as an EMI spike, in cycles.
   */
  void on_edge(uint32_t now, uint32_t glitch_filter_cycles) {
    if (have_prior_ && (now - prior_cycles_) < glitch_filter_cycles)
      return; // EMI spike: too close to the previous accepted edge
    have_prior_ = true;
    prior_cycles_ = now;
    if (count_ == 0)
      first_cycles_ = now;
    last_cycles_ = now;
    if (count_ < 255)
      ++count_; // saturate
  }

  /**
   * @brief Consumer: atomically test for a terminated burst and, if present,
   * take it.
   * @param now Current timestamp, in cycles.
   * @param gap_timeout_cycles Quiet time that terminates a burst, in cycles.
   * @param[out] out Burst snapshot, written only when true is returned.
   * @return True if a burst had terminated and was claimed.
   * @details The single consumer primitive: it recomputes completion from the
   * same `count_` it then clears, so a split burst_complete()+claim() can never
   * fold a freshly-arrived first edge into the terminated burst and then zero
   * it. The edge ISR still must not run between the test and the reset; the
   * device brackets this call in IRQ-off exactly as it did the split pair, but
   * the window can no longer be opened by an undisciplined caller.
   */
  bool try_claim(uint32_t now, uint32_t gap_timeout_cycles, BurstSnapshot *out) {
    if (count_ == 0 || (now - last_cycles_) < gap_timeout_cycles)
      return false;
    *out = BurstSnapshot{count_, first_cycles_, last_cycles_};
    count_ = 0;
    return true;
  }

  /**
   * @brief Consumer: retire the glitch-filter reference once the wire has been
   * quiet longer than the filter window.
   * @param now Current timestamp, in cycles.
   * @param glitch_filter_cycles Filter window, in cycles.
   * @details The prior accepted edge only suppresses an EMI spike within
   * glitch_filter_cycles of it; any later edge is accepted regardless, so
   * dropping the reference here is behaviourally identical for glitch
   * suppression. Its purpose is to keep `now - prior_cycles_` bounded:
   * `prior_cycles_` otherwise persists indefinitely, and after ~7.16 s of wire
   * silence the cycle counter wraps, making that modular difference
   * pseudo-random — with p ≈ glitch/2³² it lands inside the reject window and
   * falsely rejects a real edge. The flywheel poll calls this every column, so
   * a stale reference is cleared within one column of silence, long before the
   * counter can wrap. Must run under the same single-writer discipline as
   * claim() (it writes have_prior_, which the edge ISR also writes); that
   * concurrency is enforced by the device's IRQ-off discipline and is not
   * exercised by the host tests (no concurrent ISR there).
   */
  void age_prior(uint32_t now, uint32_t glitch_filter_cycles) {
    if (have_prior_ && (now - prior_cycles_) >= glitch_filter_cycles)
      have_prior_ = false;
  }

private:
  // burst_complete()+claim() are the split consumer path try_claim() fused;
  // kept private behind a test friend so production cannot reintroduce the race.
  friend struct ::hs_test::pov_sync_tests::EdgeMailboxTestAccess;

  /**
   * @brief Consumer (test-only): has a burst terminated?
   * @param now Current timestamp, in cycles.
   * @param gap_timeout_cycles Quiet time that terminates a burst, in cycles.
   * @return True if a burst exists and the wire has been quiet past the timeout.
   */
  bool burst_complete(uint32_t now, uint32_t gap_timeout_cycles) const {
    return count_ > 0 && (now - last_cycles_) >= gap_timeout_cycles;
  }

  /**
   * @brief Consumer (test-only): take the burst and reset for the next one.
   * @return A snapshot of the terminated burst.
   */
  BurstSnapshot claim() {
    const BurstSnapshot s{count_, first_cycles_, last_cycles_};
    count_ = 0;
    return s;
  }

  uint32_t count_ = 0;
  uint32_t first_cycles_ = 0;
  uint32_t last_cycles_ = 0;
  uint32_t prior_cycles_ = 0; /**< Last accepted edge ever (glitch filter ref). */
  bool have_prior_ = false;
};

// ── Telemetry (spec §8.6) ───────────────────────────────────────────────────

/**
 * @brief Counters maintained by the flywheel ISR, polled by the foreground
 * behind hs::debug.
 * @details Degradation the protocol absorbs silently must still be visible in
 * one glance of debug output.
 */
struct Telemetry {
  uint32_t symbols_accepted = 0;
  uint32_t symbols_rejected_gate = 0;     /**< §5.3 plausibility rejections. */
  uint32_t symbols_discarded_invalid = 0; /**< Invalid pulse counts (§5.2). */
  uint32_t beacons_ok = 0;
  uint32_t beacons_rejected = 0;          /**< Checksum/digit-count failures. */
  uint32_t beacon_index_corrections = 0;  /**< Missed-epoch fixes (§6.3.2). */
  uint32_t beacon_rev_mismatches = 0;
  uint32_t epochs_refractory_ignored = 0;
  uint32_t lock_transitions = 0; /**< ACQUIRE↔LOCKED edges. */
  uint32_t flips = 0;
  uint32_t emit_censored = 0;    /**< Master skipped a late boundary symbol. */
  uint32_t emit_aborted = 0;     /**< Master stopped emitting mid-burst. */
  uint32_t max_coast_halves = 0; /**< Longest run of half-revs without a snap. */
};

// ── Layer 1: the flywheel timebase (spec §4) ────────────────────────────────

/**
 * @brief Flywheel lock state: ACQUIRE (hunting) or LOCKED (disciplined).
 */
enum class LockState : uint8_t { ACQUIRE, LOCKED };

/**
 * @brief A locally-crossed boundary, reported by Flywheel::fold().
 */
struct Crossing {
  bool crossed = false;              /**< True if a boundary was crossed. */
  Boundary boundary = Boundary::NONE; /**< Which boundary was crossed. */
  uint32_t at_cycles = 0; /**< The exact (timebase) instant of the boundary. */
};

/**
 * @brief Position-from-time column generator with snap discipline.
 *
 * Position derives from the free-running clock (spec §4.1):
 *
 *   x = (x_boundary + (now − epoch_cycles)·(W/2) / cycles_per_half_rev) mod W
 *
 * with a 64-bit intermediate. `epoch_cycles` is folded forward by exactly
 * cycles_per_half_rev at every locally-crossed boundary (the rebase rule), so
 * the elapsed term never exceeds ~one half-rev plus coast and the 32-bit
 * cycle-counter wrap is structurally impossible to observe. A snap re-bases
 * the epoch to a symbol's first-edge timestamp; because position is always
 * "time since epoch", the elapsed-column compensation of spec §5.2 falls out
 * for free — classification completing ~13 columns after the boundary still
 * yields the time-correct current column.
 */
class Flywheel {
public:
  /**
   * @brief Constructs a flywheel from protocol config.
   * @param cfg Protocol configuration.
   */
  explicit Flywheel(const Config &cfg)
      : period_(cfg.cycles_per_half_rev), w_(cfg.W), gate_cols_(cfg.gate_cols),
        reject_fallback_(cfg.reject_fallback) {
    // position() reinterprets (at - epoch_cycles_) as int32, so the elapsed term
    // must never reach 2^31 cycles. Require at least kMinSafeHalfRevs half-revs
    // of coast to fit inside that window.
    HS_CHECK(period_ > 0 &&
                 period_ <= static_cast<uint32_t>(INT32_MAX) / kMinSafeHalfRevs,
             "Flywheel: cycles_per_half_rev too large — position()'s int32 "
             "elapsed window would overflow before kMinSafeHalfRevs of coast");
  }

  /**
   * @brief Boot seeding (spec §8.5): epoch = now, boundary ZERO, ACQUIRE.
   * @param now Current timestamp, in cycles.
   */
  void seed(uint32_t now) {
    epoch_cycles_ = now;
    boundary_ = Boundary::ZERO;
    lock_ = LockState::ACQUIRE;
    consecutive_rejects_ = 0;
  }

  /**
   * @brief Forces the flywheel LOCKED.
   * @details Master is the reference by definition — it never snaps and is born
   * locked (spec §2: "Master's flywheel IS the reference").
   */
  void force_lock() { lock_ = LockState::LOCKED; }

  /**
   * @brief Current column at a given time.
   * @param at Timestamp to evaluate, in cycles.
   * @return Column index in [0, W).
   * @details Signed-safe for timestamps up to ~±3.5 s around the epoch (snap
   * evaluation may look slightly into the past).
   */
  int32_t position(uint32_t at) const {
    const int64_t delta =
        static_cast<int32_t>(at - epoch_cycles_); // modular → signed
    // Fold-before-position invariant: tick() folds every crossed boundary before
    // calling position(now), so the forward elapsed stays inside the coast window
    // the constructor sized the int32 cast for. A larger delta means an unfolded
    // coast and an overflowed cast.
    HS_CHECK(delta < static_cast<int64_t>(kMinSafeHalfRevs) * period_,
             "Flywheel::position: unfolded coast — int32 elapsed cast overflowed");
    const int64_t cols = floor_div(delta * (w_ / 2), period_);
    return floor_mod(boundary_column(boundary_, w_) + cols, w_);
  }

  /**
   * @brief The rebase rule: fold the epoch forward by one half-rev if a
   * boundary has been passed.
   * @param now Current timestamp, in cycles.
   * @return A Crossing; crossed=false when no boundary was passed.
   * @details If @p now has passed the next boundary, fold the epoch forward by
   * exactly one half-rev (integer add — no drift) and report the crossing.
   * Call in a loop until it returns crossed=false; a long coast (masked window)
   * yields several crossings, each at its exact instant.
   */
  Crossing fold(uint32_t now) {
    const int32_t delta = static_cast<int32_t>(now - epoch_cycles_);
    if (delta < 0 || static_cast<uint32_t>(delta) < period_)
      return {};
    epoch_cycles_ += period_;
    boundary_ = opposite(boundary_);
    return {true, boundary_, epoch_cycles_};
  }

  /**
   * @brief Outcome of a snap attempt.
   */
  enum class SnapOutcome : uint8_t { ACCEPTED, REJECTED, REJECTED_FELL_BACK };

  /**
   * @brief Acceptance gate + re-base (spec §4.2, §5.3).
   * @param b Boundary the incoming symbol marks.
   * @param edge_cycles First-edge timestamp to re-base onto, in cycles.
   * @param error_cols Out: implied correction distance, in columns (may be
   * null).
   * @return ACCEPTED, REJECTED, or REJECTED_FELL_BACK.
   * @details LOCKED: accept only if the implied correction is ≤ G columns.
   * Because G < W/4 (Config::valid), passing the distance gate also proves the
   * named boundary is the flywheel's nearest predicted boundary — the identity
   * check is subsumed. R consecutive rejections fall back to ACQUIRE so the
   * gate can never deadlock a genuinely-lost board. ACQUIRE: hard snap, no gate
   * (the SyncBoard applies the quiet-before routing guard before calling this).
   */
  SnapOutcome snap(Boundary b, uint32_t edge_cycles, int32_t *error_cols) {
    const int32_t target = boundary_column(b, w_);
    const int32_t err = circ_dist(position(edge_cycles), target, w_);
    if (error_cols)
      *error_cols = err;
    if (lock_ == LockState::LOCKED && err > gate_cols_) {
      return note_rejection() ? SnapOutcome::REJECTED_FELL_BACK
                              : SnapOutcome::REJECTED;
    }
    epoch_cycles_ = edge_cycles;
    boundary_ = b;
    lock_ = LockState::LOCKED;
    consecutive_rejects_ = 0;
    return SnapOutcome::ACCEPTED;
  }

  /**
   * @brief Count one implausible-symbol rejection toward the ACQUIRE fallback.
   * @return True when R consecutive rejections concluded this board's own
   * timebase — not the wire — is at fault.
   * @details Shared by the snap gate and the SyncBoard's suspect-burst timeout
   * (spec §5.3: the fallback is mandatory; a gate without an escape deadlocks a
   * lost board into rejecting good symbols forever).
   */
  bool note_rejection() {
    if (lock_ != LockState::LOCKED)
      return false;
    if (++consecutive_rejects_ >= reject_fallback_) {
      lock_ = LockState::ACQUIRE;
      consecutive_rejects_ = 0;
      return true;
    }
    return false;
  }

  /**
   * @brief Current lock state.
   * @return ACQUIRE or LOCKED.
   */
  LockState lock() const { return lock_; }
  /**
   * @brief Boundary identity at the current epoch.
   * @return The boundary at epoch_cycles_.
   */
  Boundary current_boundary() const { return boundary_; }
  /**
   * @brief Current epoch timestamp.
   * @return epoch_cycles_, in cycles.
   */
  uint32_t epoch_cycles() const { return epoch_cycles_; }
  /**
   * @brief Current half-rev period.
   * @return period_, in cycles.
   */
  uint32_t cycles_per_half_rev() const { return period_; }
  /**
   * @brief §4.3 frequency trim hook (snap-only ships; tests exercise extremes).
   * @param c New half-rev period, in cycles.
   */
  void set_cycles_per_half_rev(uint32_t c) { period_ = c; }

private:
  // Min coast (in half-revs) position()'s int32 elapsed cast must survive.
  static constexpr uint32_t kMinSafeHalfRevs = 16;

  uint32_t period_;       /**< cycles_per_half_rev (optionally trimmed). */
  int32_t w_;
  int32_t gate_cols_;
  int32_t reject_fallback_;
  uint32_t epoch_cycles_ = 0;
  Boundary boundary_ = Boundary::ZERO; /**< Column identity at epoch_cycles_. */
  LockState lock_ = LockState::ACQUIRE;
  int32_t consecutive_rejects_ = 0;
};

// ── Layer 3: beacon codec (spec §6.4) ───────────────────────────────────────

/**
 * @brief Decoded index beacon: absolute effect index + revolution count
 * (mod 64).
 */
struct BeaconFrame {
  int32_t effect_index = 0; /**< Effect index, 0–63. */
  uint32_t rev_count = 0;   /**< Revolution within the effect, mod 64. */
};

/**
 * @brief Encodes a beacon as five base-8 digits.
 * @param effect_index Effect index (low 6 bits used).
 * @param rev_count Revolution count (low 6 bits used).
 * @param out Out: five digits [index_hi, index_lo, rev_hi, rev_lo, checksum].
 * @details Each digit is transmitted as a burst of (digit+1) pulses.
 */
constexpr void encode_beacon_digits(int32_t effect_index, uint32_t rev_count,
                                    uint8_t out[5]) {
  const uint32_t idx = static_cast<uint32_t>(effect_index) & 63u;
  const uint32_t rev = rev_count & 63u;
  out[0] = static_cast<uint8_t>(idx >> 3);
  out[1] = static_cast<uint8_t>(idx & 7u);
  out[2] = static_cast<uint8_t>(rev >> 3);
  out[3] = static_cast<uint8_t>(rev & 7u);
  // Position-weighted Σ(i+1)·dᵢ mod 8: catches digit transpositions and
  // compensating miscounts a plain sum would miss.
  out[4] = static_cast<uint8_t>(
      (1u * out[0] + 2u * out[1] + 3u * out[2] + 4u * out[3]) & 7u);
}

/**
 * @brief Smallest signed revolution slip that maps @p current's 6-bit residue
 * onto a beacon's @p beacon_rev_count (spec §6.4 rev cross-check).
 * @param beacon_rev_count The beacon's rev_count digit pair (0–63).
 * @param current The board's current rev_in_effect (any magnitude).
 * @return A delta in [-32, 31]. The `+96 %64 -32` fold resolves the 63↔0 mod-64
 * seam, so a small slip across the wrap reads as a small *signed* delta (e.g.
 * current residue 62 vs beacon 0 → +2) rather than a ~64-rev jump. Slips of
 * magnitude ≥ 32 are aliased by the mod-64 (the protocol guarantees agreement
 * within 32 revs), so the caller treats them as the nearest-residue correction.
 */
constexpr int32_t beacon_rev_resync_delta(uint32_t beacon_rev_count,
                                          uint32_t current) {
  return ((static_cast<int32_t>(beacon_rev_count & 63u) -
           static_cast<int32_t>(current & 63u) + 96) %
          64) -
         32;
}

/**
 * @brief Assembles data bursts into beacon frames. Integrity by rejection
 * (spec §6.4): any out-of-range digit, stale partial frame, or checksum
 * mismatch drops the whole frame — the next beacon is ≤ 2 s away.
 */
class BeaconParser {
public:
  /**
   * @brief Feed one data burst into the frame assembler.
   * @param s The data burst.
   * @param cfg Protocol configuration.
   * @param out Out: the parsed frame, valid only when the return is true.
   * @param rejected Out: set when a frame (not merely a digit) was discarded.
   * @return True when a complete, checksum-valid frame was parsed into @p out.
   */
  bool feed(const BurstSnapshot &s, const Config &cfg, BeaconFrame *out,
            bool *rejected) {
    *rejected = false;
    if (n_ > 0 &&
        (s.first_cycles - last_first_cycles_) > cfg.interdigit_timeout_cycles()) {
      n_ = 0; // stale partial frame: a new frame may start with this burst
      *rejected = true;
    }
    if (s.count < 1 || s.count > 8) {
      if (n_ > 0)
        *rejected = true;
      n_ = 0;
      return false;
    }
    last_first_cycles_ = s.first_cycles;
    digits_[n_++] = static_cast<uint8_t>(s.count - 1);
    if (n_ < 5)
      return false;
    n_ = 0;
    // Position-weighted checksum (see encode_beacon_digits).
    if (((1u * digits_[0] + 2u * digits_[1] + 3u * digits_[2] +
          4u * digits_[3]) &
         7u) != digits_[4]) {
      *rejected = true;
      return false;
    }
    out->effect_index = digits_[0] * 8 + digits_[1];
    out->rev_count = static_cast<uint32_t>(digits_[2]) * 8u + digits_[3];
    return true;
  }

  /**
   * @brief Discards any partial frame in progress.
   */
  void reset() { n_ = 0; }
  /**
   * @brief Whether a partial frame is being assembled.
   * @return True if at least one digit has been buffered.
   */
  bool active() const { return n_ > 0; }

private:
  int32_t n_ = 0;
  uint8_t digits_[5] = {};
  uint32_t last_first_cycles_ = 0;
};

// ── Layer 3: content tracker (spec §6.1, §6.3) ──────────────────────────────

/**
 * @brief Per-board content state: which effect, which revolution, and the
 * deadline-scheduled epoch commit. Owned by the flywheel ISR via SyncBoard.
 */
struct ContentTracker {
  bool identity_known = false; /**< False until an epoch/beacon names the index. */
  int32_t effect_index = 0;    /**< Currently displayed effect index. */
  /**
   * @brief ZERO crossings since effect start.
   * @details For a beacon-joined board this starts from the beacon's mod-64
   * value — it only feeds the master's schedule and cross-check telemetry.
   */
  uint32_t rev_in_effect = 0;
  bool commit_pending = false;       /**< An epoch commit is scheduled. */
  uint32_t commit_in_revs = 0;       /**< ZERO crossings remaining until the
                                        absolute B+R+K boundary, with j already
                                        subtracted for the repeat this board
                                        heard (NOT the announce-phase length). */
  uint32_t refractory_revs_left = 0; /**< EPOCH dedup window (spec §6.1). */

  /**
   * @brief Accepts an EPOCH symbol and may open a commit window.
   * @param cfg Protocol configuration.
   * @return True if it opened a commit window (false inside the refractory
   * window — the §6.3 redundancy repeats land here).
   * @details The commit boundary is ABSOLUTE: B + R + K, where B is the primary
   * copy's boundary (R = epoch_repeats, K = commit_revs). Which copy of the
   * train this is (j, 0 = primary) is inferred from the shared revolution
   * count — the master starts the train exactly when rev_in_effect reaches
   * revs_per_effect, and by the time a symbol is consumed the local crossing
   * for its boundary has already incremented rev_in_effect (classification
   * completes ~13 columns after the boundary instant). So every board that
   * hears ANY copy counts down to the same boundary, and hearing a repeat
   * instead of the primary cannot skew the commit (§6.3.1). A board whose
   * revolution count is not absolute (it beacon-joined mid-effect, §6.4)
   * lands outside the train window and falls back to j = 0 — it commits up
   * to j revolutions late, an epoch-bounded degradation confined to that
   * case. The resulting counter slip is then resynced from the next beacon's
   * rev cross-check (handle_beacon_burst), so every subsequent epoch is
   * lockstep.
   */
  bool on_epoch_symbol(const Config &cfg) {
    if (refractory_revs_left > 0)
      return false;
    uint32_t j = 0;
    if (rev_in_effect >= cfg.revs_per_effect &&
        rev_in_effect - cfg.revs_per_effect <=
            static_cast<uint32_t>(cfg.epoch_repeats))
      j = rev_in_effect - cfg.revs_per_effect;
    commit_pending = true;
    commit_in_revs =
        cfg.commit_revs + static_cast<uint32_t>(cfg.epoch_repeats) - j;
    refractory_revs_left = cfg.refractory_revs;
    return true;
  }

  /**
   * @brief Whether the construction window opens at this crossing/accept.
   * @param cfg Protocol configuration.
   * @return True at the single crossing/accept where the construction window —
   * the last K revolutions before the commit boundary — opens.
   * @details commit_in_revs strictly decreases, so == fires exactly once per
   * window. The announce phase before it is what gives a board that hears only
   * the last repeat the full K-revolution construction budget; the outgoing
   * effect keeps playing through it, so the dark window stays K revolutions.
   */
  bool construction_opens(const Config &cfg) const {
    return commit_pending && commit_in_revs == cfg.commit_revs;
  }

  /**
   * @brief Whether the board is in the construction window.
   * @param cfg Protocol configuration.
   * @return True throughout the construction window (fail-dark, spec §6.1).
   */
  bool constructing(const Config &cfg) const {
    return commit_pending && commit_in_revs <= cfg.commit_revs;
  }

  /**
   * @brief Advances revolution bookkeeping on a ZERO boundary flip.
   * @param cfg Protocol configuration.
   * @return True when the commit deadline is reached — the board must swap to
   * the (already constructed) next effect at exactly this boundary.
   */
  bool on_zero_crossing(const Config &cfg) {
    ++rev_in_effect;
    if (refractory_revs_left > 0)
      --refractory_revs_left;
    if (commit_pending && --commit_in_revs == 0) {
      commit_pending = false;
      effect_index = (effect_index + 1) % cfg.effect_count;
      rev_in_effect = 0;
      return true;
    }
    return false;
  }
};

// ── Master symbol emitter (spec §5.2 generation, §6.4 beacon) ───────────────

/**
 * @brief Schedules pulses in cycle time; the flywheel ISR asks once per entry
 * whether to emit a pulse this tick (pin write first, LED work after).
 *
 * Self-censoring: a boundary symbol whose first pulse would start more than
 * ~½ column late is skipped entirely; lateness detected mid-burst stops the
 * remaining pulses ("never emit a lie" — the §5.3 gate and refractory window
 * absorb the truncated count downstream). Beacon digit bursts get the same
 * treatment; a truncated frame fails its checksum and is dropped whole.
 */
class SymbolEmitter {
public:
  /**
   * @brief Schedules a boundary symbol for emission.
   * @param symbol Symbol to emit (INVALID is rejected).
   * @param at_cycles Boundary instant the burst should start at, in cycles.
   * @param now Current timestamp, in cycles.
   * @param cfg Protocol configuration.
   * @return False if the symbol was self-censored (caller counts it).
   */
  bool schedule_boundary(Symbol symbol, uint32_t at_cycles, uint32_t now,
                         const Config &cfg) {
    if (symbol == Symbol::INVALID)
      return false;
    // Signed lateness so a future boundary isn't read as a huge positive
    // lateness through the unsigned wrap of now - at_cycles.
    const int32_t lateness = static_cast<int32_t>(now - at_cycles);
    if (lateness > static_cast<int32_t>(cfg.late_censor_cycles()))
      return false; // late at the boundary: skip the whole symbol
    if (pulses_left_ > 0 || queue_pos_ < queue_len_)
      return false; // wire still busy (defensive; in-range config never hits it)
    pulses_left_ = symbol_pulse_count(symbol);
    next_due_ = at_cycles;
    pitch_ = cfg.pulse_pitch_cycles();
    return true;
  }

  /**
   * @brief Queues the five digit bursts of a beacon frame.
   * @param digits The five base-8 beacon digits.
   * @param now Current timestamp (frame start), in cycles.
   * @param cfg Protocol configuration.
   * @details Called when the master reaches the beacon point (x ≈ W/4).
   */
  void schedule_beacon(const uint8_t digits[5], uint32_t now,
                       const Config &cfg) {
    if (pulses_left_ > 0 || queue_pos_ < queue_len_)
      return; // defensive: never interleave with an active emission
    uint32_t start = now;
    const uint32_t col = cfg.cycles_per_column();
    for (int i = 0; i < 5; ++i) {
      queue_[i] = {start, static_cast<uint32_t>(digits[i]) + 1u};
      // Next burst starts after this one's span (digits[i] pitches) plus a
      // gap the decoder is guaranteed to see as burst-terminating.
      start += digits[i] * cfg.beacon_pitch_cycles() +
               static_cast<uint32_t>(cfg.gap_timeout_cols + 1) * col;
    }
    queue_len_ = 5;
    queue_pos_ = 0;
  }

  /**
   * @brief Per-tick: should the ISR emit one pulse right now?
   * @param now Current timestamp, in cycles.
   * @param cfg Protocol configuration.
   * @param aborted Out: reports a mid-burst self-censor (telemetry).
   * @return True if a pulse should be emitted this tick.
   */
  bool tick(uint32_t now, const Config &cfg, bool *aborted) {
    *aborted = false;
    if (pulses_left_ == 0 && queue_pos_ >= queue_len_ && queue_len_ != 0)
      queue_len_ = queue_pos_ = 0; // beacon frame drained: emitter idle again
    if (pulses_left_ == 0 && queue_pos_ < queue_len_) {
      const PendingBurst &b = queue_[queue_pos_];
      if (static_cast<int32_t>(now - b.start_cycles) >= 0) {
        pulses_left_ = b.pulses;
        next_due_ = b.start_cycles;
        pitch_ = cfg.beacon_pitch_cycles();
        ++queue_pos_;
      }
    }
    if (pulses_left_ == 0)
      return false;
    if (static_cast<int32_t>(now - next_due_) < 0)
      return false; // next pulse not due yet
    if ((now - next_due_) > cfg.late_censor_cycles()) {
      // Masked past the lateness budget mid-emission: stop, and drop any
      // queued beacon digits too (a partial frame must fail, not mislead).
      pulses_left_ = 0;
      queue_len_ = queue_pos_ = 0;
      *aborted = true;
      return false;
    }
    --pulses_left_;
    next_due_ += pitch_;
    return true;
  }

  /**
   * @brief Whether the emitter has nothing queued or in flight.
   * @return True if no pulses or beacon bursts remain.
   */
  bool idle() const { return pulses_left_ == 0 && queue_pos_ >= queue_len_; }

private:
  /**
   * @brief A queued beacon-digit burst awaiting its start time.
   */
  struct PendingBurst {
    uint32_t start_cycles = 0; /**< When the burst should begin, in cycles. */
    uint32_t pulses = 0;       /**< Number of pulses in the burst. */
  };
  PendingBurst queue_[5] = {};
  int32_t queue_len_ = 0;
  int32_t queue_pos_ = 0;
  uint32_t pulses_left_ = 0;
  uint32_t next_due_ = 0;
  uint32_t pitch_ = 0;
};

// ── The per-board engine ────────────────────────────────────────────────────

/**
 * @brief What the device ISR (or host sim) must do after one flywheel tick.
 */
struct TickActions {
  /**
   * @brief Column to pack/display, or -1.
   * @details -1 means unchanged since last tick, or display is dark.
   */
  int32_t render_column = -1;
  /**
   * @brief Display black.
   * @details Set for ACQUIRE, no identity, or the epoch construction window
   * (spec §6.1).
   */
  bool dark = false;
  bool flip = false;          /**< Call advance_display() on the live effect. */
  bool zero_crossing = false; /**< A ZERO boundary flip happened this tick. */
  /**
   * @brief ZERO crossing on the join grid.
   * @details A board with no live effect takes its pending one here
   * (Config::join_grid_revs).
   */
  bool join_boundary = false;
  /**
   * @brief B+R+K deadline reached.
   * @details Swap to the pending effect at this boundary (trap if it is not
   * ready).
   */
  bool commit = false;
  bool pulse = false; /**< Master: emit one sync pulse (pin-first). */
};

/**
 * @brief One board's complete sync state machine. The flywheel ISR calls
 * tick() once per wake-up; the sync-wire edge ISR calls on_sync_edge(); the
 * foreground polls build_word() to learn which effect to construct.
 *
 * Single-writer: everything except EdgeMailbox is written only from tick()
 * (the flywheel ISR). The mailbox handoff is the one explicitly-synchronized
 * point (claim under a brief IRQ-off window on the device).
 */
class SyncBoard {
public:
  /**
   * @brief Constructs a board from protocol config.
   * @param cfg Protocol configuration.
   */
  explicit SyncBoard(const Config &cfg)
      : cfg_(cfg), fly_(cfg) {}

  /**
   * @brief Boot/reboot seeding (spec §8.5).
   * @param now Current timestamp, in cycles.
   * @param is_master True for the reference board.
   * @details Master is born locked with identity (effect 0, rev 0) — it is the
   * reference. Downstream boards start in ACQUIRE, dark, and join via boundary
   * symbols + beacon.
   */
  void seed(uint32_t now, bool is_master) {
    is_master_ = is_master;
    fly_.seed(now);
    gate_ = FlipGate{};
    content_ = ContentTracker{};
    telemetry_ = Telemetry{};
    // A reboot must not inherit wire state from the prior incarnation: a stale
    // mailbox burst would feed ACQUIRE's unconditional hard-snap, and a stale
    // emitter queue would resume a half-sent beacon/boundary train (spec §8.5).
    mailbox_ = EdgeMailbox{};
    emitter_ = SymbolEmitter{};
    beacon_parser_.reset();
    last_rendered_x_ = -1;
    halves_since_snap_ = 0;
    have_prev_burst_ = false;
    suspect_pending_ = false;
    epoch_emits_left_ = 0;
    beacon_done_this_rev_ = false;
    build_gen_ = 0;
    if (is_master) {
      fly_.force_lock();
      content_.identity_known = true;
      content_.effect_index = 0;
      publish_build(0);
    }
  }

  /**
   * @brief Sync-wire edge ISR entry point (publisher; downstream boards only).
   * @param now Edge timestamp, in cycles.
   */
  void on_sync_edge(uint32_t now) {
    mailbox_.on_edge(now, cfg_.glitch_filter_cycles);
  }

  /**
   * @brief Device-side accessor for the IRQ-off mailbox handoff.
   * @return Reference to the edge mailbox.
   */
  EdgeMailbox &mailbox() { return mailbox_; }
  /**
   * @brief Burst-terminating gap.
   * @return Quiet time that terminates a burst, in cycles.
   */
  uint32_t gap_timeout_cycles() const { return cfg_.gap_timeout_cycles(); }
  /**
   * @brief Glitch-filter window.
   * @return Minimum accepted edge spacing, in cycles.
   */
  uint32_t glitch_filter_cycles() const { return cfg_.glitch_filter_cycles; }

  /**
   * @brief One flywheel wake-up.
   * @param now Current timestamp, in cycles.
   * @param burst The claimed mailbox burst if one completed (see EdgeMailbox),
   * else nullptr.
   * @return The actions the caller must perform after this tick.
   */
  TickActions tick(uint32_t now, const BurstSnapshot *burst) {
    TickActions a{};
    if (burst)
      handle_burst(*burst, a);

    // Suspect-burst timeout: a lone far burst held pending in handle_burst that
    // saw no follow-up was not beacon data — count it as a gate rejection so a
    // corrupted-timebase board still reaches the §5.3 ACQUIRE fallback. The
    // signed re-check rejects a wrapped modular difference.
    if (suspect_pending_ &&
        (now - suspect_last_cycles_) > cfg_.interdigit_timeout_cycles() &&
        static_cast<int32_t>(now - suspect_last_cycles_) > 0) {
      suspect_pending_ = false;
      ++telemetry_.symbols_rejected_gate;
      if (fly_.note_rejection())
        ++telemetry_.lock_transitions;
    }

    // Fold every locally-crossed boundary (usually 0 or 1; several after a
    // long masked coast).
    for (;;) {
      const Crossing c = fly_.fold(now);
      if (!c.crossed)
        break;
      if (halves_since_snap_ < 0xFFFFFFFFu)
        ++halves_since_snap_;
      if (halves_since_snap_ > telemetry_.max_coast_halves)
        telemetry_.max_coast_halves = halves_since_snap_;
      apply_flip(c.boundary, a);
      if (is_master_)
        master_on_crossing(c, now, a);
    }

    if (is_master_)
      maybe_schedule_beacon(now);

    bool aborted = false;
    if (is_master_ && emitter_.tick(now, cfg_, &aborted))
      a.pulse = true;
    if (aborted)
      ++telemetry_.emit_aborted;

    // Render decision: dark when phase/identity is missing or during the epoch
    // construction window (spec §6.1).
    a.dark = fly_.lock() != LockState::LOCKED || !content_.identity_known ||
             content_.constructing(cfg_);
    if (!a.dark) {
      const int32_t x = fly_.position(now);
      if (x != last_rendered_x_) { // idempotent wake-up contract (spec §4.1)
        a.render_column = x;
        last_rendered_x_ = x;
      }
    } else {
      last_rendered_x_ = -1;
    }
    return a;
  }

  // ── Foreground interface (read-only; single aligned-word reads) ──────────

  /**
   * @brief Current build request.
   * @return (generation << 8) | effect_index.
   * @details The foreground constructs the named effect whenever the generation
   * changes, then publishes it to the driver's pending slot.
   */
  uint32_t build_word() const { return build_word_; }
  /**
   * @brief Extracts the effect index from a build word.
   * @param word A build word.
   * @return The effect index (low 8 bits).
   */
  static int32_t build_index_of(uint32_t word) {
    return static_cast<int32_t>(word & 0xFF);
  }
  /**
   * @brief Extracts the generation from a build word.
   * @param word A build word.
   * @return The generation counter (high bits).
   */
  static uint32_t build_gen_of(uint32_t word) { return word >> 8; }

  /**
   * @brief Telemetry counters.
   * @return Const reference to the telemetry block.
   * @warning The returned reference aliases the live, ISR-mutated counter block
   * — it is NOT a snapshot. A device caller copying it field-by-field outside an
   * IRQ-off window can latch a torn block (a mix of pre- and post-increment
   * fields); snapshot the copy under a brief `__disable_irq()`/`__enable_irq()`
   * bracket (spec §8.2), as POVSegmented's health-telemetry poll does.
   */
  const Telemetry &telemetry() const { return telemetry_; }
  /**
   * @brief Content-layer state.
   * @return Const reference to the content tracker.
   */
  const ContentTracker &content() const { return content_; }
  /**
   * @brief Current lock state.
   * @return ACQUIRE or LOCKED.
   */
  LockState lock() const { return fly_.lock(); }
  /**
   * @brief The flywheel timebase.
   * @return Const reference to the flywheel.
   */
  const Flywheel &flywheel() const { return fly_; }
  /**
   * @brief Mutable flywheel access (tests only).
   * @return Mutable reference to the flywheel.
   * @warning Host tests only. The flywheel is ISR-owned under the spec §8
   * single-writer model, so this must NEVER be called from device or foreground
   * code, and only before the ISRs are attached — mutating it once the ISRs are
   * live races the single writer and breaks the invariant the design rests on.
   */
  Flywheel &flywheel_mut() { return fly_; }
  /**
   * @brief Mutable content-tracker access (tests only).
   * @return Mutable reference to the content tracker.
   * @warning Host tests only. The content tracker is ISR-owned under the spec §8
   * single-writer model, so this must NEVER be called from device or foreground
   * code, and only before the ISRs are attached — mutating it once the ISRs are
   * live races the single writer and breaks the invariant the design rests on.
   */
  ContentTracker &content_mut() { return content_; }
  /**
   * @brief Protocol configuration.
   * @return Const reference to the config.
   */
  const Config &config() const { return cfg_; }

private:
  // ── Flip + content events (both paths funnel here) ───────────────────────

  /**
   * @brief Funnels both flip paths (local crossing and symbol backstop) through
   * the dedup gate and runs the per-boundary content events.
   * @param b Boundary being crossed.
   * @param a In/out: tick actions updated with flip/commit/join effects.
   */
  void apply_flip(Boundary b, TickActions &a) {
    if (!gate_.try_flip(b))
      return;
    a.flip = true;
    ++telemetry_.flips;
    if (b != Boundary::ZERO)
      return;
    a.zero_crossing = true;
    beacon_done_this_rev_ = false;
    if (content_.identity_known) {
      if (content_.on_zero_crossing(cfg_))
        a.commit = true; // B+R+K reached; driver swaps in the pending effect
      else if (content_.construction_opens(cfg_))
        // Last K revolutions: construct the next effect now.
        publish_build((content_.effect_index + 1) % cfg_.effect_count);
      else if (!content_.commit_pending &&
               (content_.rev_in_effect % cfg_.join_grid_revs) == 0)
        // Marks "a late joiner could snap in here"; the shell acts on it only
        // when it has no live effect.
        a.join_boundary = true;
    }
  }

  // ── Receive path (downstream) ─────────────────────────────────────────────

  /**
   * @brief Routes a claimed burst to the boundary-symbol or beacon path and
   * runs the acceptance gate.
   * @param s The terminated burst.
   * @param a In/out: tick actions updated with any resulting flip.
   */
  void handle_burst(const BurstSnapshot &s, TickActions &a) {
    const bool had_prev = have_prev_burst_;
    const uint32_t prev_end = prev_burst_end_;
    prev_burst_end_ = s.last_cycles;
    have_prev_burst_ = true;

    // Any follow-up burst inside the interdigit window proves the pending
    // suspect (see below) was the head of a beacon data train: clear it.
    if (suspect_pending_ &&
        (s.first_cycles - suspect_last_cycles_) <=
            cfg_.interdigit_timeout_cycles())
      suspect_pending_ = false;

    if (fly_.lock() == LockState::LOCKED) {
      // Demarcation (spec §6.4): a burst whose first edge lands far from a
      // predicted boundary is beacon data, not a boundary symbol.
      const int32_t pos = fly_.position(s.first_cycles);
      const int32_t to_zero = circ_dist(pos, 0, cfg_.W);
      const int32_t to_half = circ_dist(pos, cfg_.W / 2, cfg_.W);
      if (to_zero > cfg_.gate_cols && to_half > cfg_.gate_cols) {
        handle_beacon_burst(s);
        // A lone far burst is indistinguishable from a beacon's first digit
        // until a train follows: hold it as a suspect, and tick()'s timeout
        // converts it to a gate rejection if the wire stays silent.
        const bool isolated =
            !had_prev ||
            (s.first_cycles - prev_end) >= cfg_.acquire_quiet_cycles();
        if (isolated && classify_count(s.count) != Symbol::INVALID) {
          suspect_pending_ = true;
          suspect_last_cycles_ = s.last_cycles;
        }
        return;
      }
    } else {
      // ACQUIRE quiet-before guard: a burst following close on another is a
      // beacon digit train, not an isolated boundary symbol — don't hard-snap.
      if (had_prev &&
          (s.first_cycles - prev_end) < cfg_.acquire_quiet_cycles()) {
        handle_beacon_burst(s);
        return;
      }
      beacon_parser_.reset(); // isolated burst: any partial frame is stale
    }

    const Symbol sym = classify_count(s.count);
    if (sym == Symbol::INVALID) {
      ++telemetry_.symbols_discarded_invalid;
      return;
    }
    const Boundary b = symbol_boundary(sym);
    const bool was_locked = fly_.lock() == LockState::LOCKED;
    int32_t err = 0;
    const Flywheel::SnapOutcome r = fly_.snap(b, s.first_cycles, &err);
    if (r != Flywheel::SnapOutcome::ACCEPTED) {
      ++telemetry_.symbols_rejected_gate;
      if (r == Flywheel::SnapOutcome::REJECTED_FELL_BACK)
        ++telemetry_.lock_transitions;
      return;
    }
    ++telemetry_.symbols_accepted;
    if (!was_locked)
      ++telemetry_.lock_transitions;
    halves_since_snap_ = 0;
    // MUST precede on_epoch_symbol: a ZERO_EPOCH folds rev_in_effect here so the
    // j-inference below reads the post-fold rev (§6.3.1). Deduped against the
    // later fold-loop apply_flip. See test_epoch_same_tick_burst_fold.
    apply_flip(b, a);
    if (sym == Symbol::ZERO_EPOCH && content_.identity_known) {
      if (content_.on_epoch_symbol(cfg_)) {
        // A board that heard only the last repeat opens the window at the accept.
        if (content_.construction_opens(cfg_))
          publish_build((content_.effect_index + 1) % cfg_.effect_count);
      } else {
        ++telemetry_.epochs_refractory_ignored;
      }
    }
  }

  /**
   * @brief Feeds a data burst to the beacon parser and applies a completed
   * frame (join, missed-epoch correction, or rev-counter resync).
   * @param s The data burst.
   */
  void handle_beacon_burst(const BurstSnapshot &s) {
    BeaconFrame f{};
    bool rejected = false;
    const bool ok = beacon_parser_.feed(s, cfg_, &f, &rejected);
    if (rejected)
      ++telemetry_.beacons_rejected;
    if (!ok)
      return;
    ++telemetry_.beacons_ok;
    const int32_t idx = f.effect_index % cfg_.effect_count;
    if (!content_.identity_known) {
      // Join (spec §6.4): adopt (effect, rev). Never assume index 0.
      content_.identity_known = true;
      content_.effect_index = idx;
      content_.rev_in_effect = f.rev_count;
      publish_build(idx);
    } else if (content_.commit_pending) {
      // Do NOT publish_build mid-window: pending_gen_ must stay stable from
      // construction-open to commit, the precondition the commit-time HS_CHECK
      // relies on. The next post-commit beacon re-verifies the index.
    } else if (idx != content_.effect_index) {
      // Missed epoch (all repeats): correct within ≤16 revs (spec §6.3.2).
      content_.effect_index = idx;
      content_.rev_in_effect = f.rev_count;
      ++telemetry_.beacon_index_corrections;
      publish_build(idx);
    } else if (f.rev_count != (content_.rev_in_effect & 63u)) {
      // The schedule counter slipped against the master's; left alone it skews
      // every later epoch commit by mis-inferred j. Resync via the signed
      // mod-64 difference, which recovers any slip under 32 revolutions.
      ++telemetry_.beacon_rev_mismatches;
      const int32_t d =
          beacon_rev_resync_delta(f.rev_count, content_.rev_in_effect);
      const int64_t fixed = static_cast<int64_t>(content_.rev_in_effect) + d;
      content_.rev_in_effect =
          fixed >= 0 ? static_cast<uint32_t>(fixed) : f.rev_count;
    }
  }

  // ── Conduct + emit path (master) ──────────────────────────────────────────

  /**
   * @brief Master conductor: emit the boundary symbol for a crossing and run
   * the epoch-train schedule.
   * @param c The crossing just folded.
   * @param now Current timestamp, in cycles.
   * @details The unnamed TickActions parameter is reserved for symmetry with
   * the receive path.
   */
  void master_on_crossing(const Crossing &c, uint32_t now, TickActions &) {
    Symbol sym = Symbol::HALF;
    if (c.boundary == Boundary::ZERO) {
      // Conductor (spec §6.1): when the effect's revolutions elapse, start an
      // EPOCH train — primary copy plus R repeats on following ZERO boundaries.
      if (epoch_emits_left_ == 0 && !content_.commit_pending &&
          content_.rev_in_effect >= cfg_.revs_per_effect) {
        epoch_emits_left_ = 1 + cfg_.epoch_repeats;
        if (content_.on_epoch_symbol(cfg_) &&
            content_.construction_opens(cfg_))
          publish_build((content_.effect_index + 1) % cfg_.effect_count);
      }
      if (epoch_emits_left_ > 0) {
        sym = Symbol::ZERO_EPOCH;
      } else {
        sym = Symbol::ZERO;
      }
    }
    // Spend a redundancy repeat only on a symbol that actually reaches the wire;
    // a censored ZERO_EPOCH never propagated.
    if (!emitter_.schedule_boundary(sym, c.at_cycles, now, cfg_))
      ++telemetry_.emit_censored;
    else if (sym == Symbol::ZERO_EPOCH)
      --epoch_emits_left_;
  }

  /**
   * @brief Master: queue a beacon frame once per revolution at the beacon point
   * when one is due.
   * @param now Current timestamp, in cycles.
   */
  void maybe_schedule_beacon(uint32_t now) {
    // Beacon point: x ≈ W/4, mid-way through the ZERO→HALF half-rev where the
    // wire is otherwise quiet (spec §6.4).
    if (beacon_done_this_rev_ || fly_.current_boundary() != Boundary::ZERO)
      return;
    // Silent during the commit window: a beacon here broadcasts the outgoing
    // index, and a board joining off it would adopt stale identity.
    if (content_.commit_pending)
      return;
    // Frame span (≤55 cols at base-8 digits) ends well short of HALF even when a
    // masked window pushes the start later; no separate beacon-start late-bound.
    //
    // A coalesced coast can jump position from < W/4 straight past the beacon
    // point (and even past HALF) in one wake, leaving beacon_done_this_rev_ unset
    // while current_boundary() has already advanced — so this revolution emits no
    // beacon. That is an accepted skip, not a missed-emission bug: the protocol
    // self-heals on the next due beacon (≤ 2 s).
    if (fly_.position(now) < cfg_.W / 4)
      return;
    beacon_done_this_rev_ = true;
    const uint32_t rev = content_.rev_in_effect;
    const bool due = (rev % cfg_.beacon_period_revs) == 1u ||
                     (rev >= 1u && rev <= static_cast<uint32_t>(
                                              cfg_.epoch_repeats));
    if (!due)
      return;
    uint8_t digits[5];
    encode_beacon_digits(content_.effect_index, rev, digits);
    emitter_.schedule_beacon(digits, now, cfg_);
  }

  /**
   * @brief Publishes a build request for @p index, bumping the generation.
   * @param index Effect index the foreground should construct.
   */
  void publish_build(int32_t index) {
    ++build_gen_;
    build_word_ = (build_gen_ << 8) | static_cast<uint32_t>(index & 0xFF);
  }

  // ── State ─────────────────────────────────────────────────────────────────

  Config cfg_;
  Flywheel fly_;
  FlipGate gate_;
  ContentTracker content_;
  BeaconParser beacon_parser_;
  SymbolEmitter emitter_;
  EdgeMailbox mailbox_;
  Telemetry telemetry_;

  bool is_master_ = false;
  int32_t last_rendered_x_ = -1;
  uint32_t halves_since_snap_ = 0;
  bool have_prev_burst_ = false;
  uint32_t prev_burst_end_ = 0;
  bool suspect_pending_ = false;     /**< Lone far burst awaiting train/timeout. */
  uint32_t suspect_last_cycles_ = 0;
  uint32_t epoch_emits_left_ = 0;
  bool beacon_done_this_rev_ = false;
  uint32_t build_gen_ = 0;
  // `volatile`, not std::atomic: SyncBoard is move-assigned at setup (sync_ =
  // SyncBoard(cfg)), and an atomic member would delete the implicit move-assign.
  // volatile keeps the type move-assignable and still forces the foreground poll
  // to re-load (single aligned word, one logical writer at a time). volatile
  // suppresses elision/reordering but NOT tearing; the load-bearing assumption is
  // that a naturally-aligned 32-bit load/store is atomic on the single-core ARM
  // target, so the reader never observes a half-written word. That tear-free
  // guarantee comes from the target ISA, not the C++ memory model.
  volatile uint32_t build_word_ = 0; /**< (gen << 8) | index; foreground-read. */
};

} // namespace sync
} // namespace pov
