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

namespace pov {
namespace sync {

// ── Small integer helpers (signed-correct division/modulo) ─────────────────

/// Floor division (rounds toward -inf, unlike C++ truncation).
constexpr int64_t floor_div(int64_t a, int64_t b) {
  const int64_t q = a / b;
  const int64_t r = a % b;
  return (r != 0 && ((r < 0) != (b < 0))) ? q - 1 : q;
}

/// Non-negative modulo in [0, m).
constexpr int32_t floor_mod(int64_t a, int32_t m) {
  const int32_t r = static_cast<int32_t>(a % m);
  return r < 0 ? r + m : r;
}

/// Circular distance between two columns on a ring of width w (in [0, w/2]).
constexpr int32_t circ_dist(int32_t a, int32_t b, int32_t w) {
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
  int32_t W = 288;                         ///< Columns per revolution.
  uint32_t cycles_per_half_rev = 37500000; ///< Timebase constant (spec §4.1).
  uint32_t glitch_filter_cycles = 60000;   ///< Min edge spacing (~100 µs).

  // Symbol wire (spec §5.2): pitches/timeouts in columns.
  int32_t pulse_pitch_cols = 2;  ///< Boundary-burst pulse pitch (> mask M).
  int32_t beacon_pitch_cols = 1; ///< Beacon digit pulse pitch (checksummed).
  int32_t gap_timeout_cols = 4;  ///< Quiet time that terminates a burst.

  // Acceptance gate (spec §5.3).
  int32_t gate_cols = 4;       ///< G: max plausible snap correction (LOCKED).
  int32_t reject_fallback = 4; ///< R: consecutive rejections before ACQUIRE.
  int32_t acquire_quiet_cols = 16; ///< Quiet-before guard for ACQUIRE snaps.

  // Content layer (spec §6).
  uint32_t revs_per_effect = 960; ///< Effect duration in revolutions (120 s).
  int32_t epoch_repeats = 3;      ///< EPOCH redundancy repeats (spec §6.3).
  uint32_t refractory_revs = 16;  ///< EPOCH dedup window (spec §6.1).
  /// K: the construction window, the last K revolutions before the commit.
  /// The commit boundary itself is B + epoch_repeats + commit_revs from the
  /// train's primary copy at B, so a board hearing any repeat still commits
  /// in lockstep with the full K-rev construction budget (spec §6.1, §6.3.1).
  uint32_t commit_revs = 2;
  uint32_t beacon_period_revs = 16;            ///< Beacon cadence (spec §6.4).
  int32_t beacon_interdigit_timeout_cols = 24; ///< Stale-frame reset bound.
  int32_t effect_count = 1; ///< Roster length (epoch wraps mod this).
  /// Boards take a constructed effect live only at revolutions ≡ 0 mod this
  /// grid. At boot every board (master included) therefore goes live at the
  /// SAME crossing with frame counters aligned, instead of master leading by
  /// however long downstream identity took. Must divide 64 so a beacon's
  /// mod-64 revolution count lands on the same grid as the master's true
  /// count. Mid-show rejoins wait ≤ grid revolutions (well inside the ~2 s
  /// rejoin budget, spec §9.1).
  uint32_t join_grid_revs = 4;

  // Derived cycle quantities. cycles_per_column truncates (~2.6 ppm); it is
  // used only for pitches/thresholds — position math divides by the exact
  // cycles_per_half_rev (spec §4.1, "stored period is per-half-rev").
  constexpr uint32_t cycles_per_column() const {
    return cycles_per_half_rev / static_cast<uint32_t>(W / 2);
  }
  constexpr uint32_t col_cycles(int32_t cols) const {
    return cycles_per_column() * static_cast<uint32_t>(cols);
  }
  constexpr uint32_t gap_timeout_cycles() const {
    return col_cycles(gap_timeout_cols);
  }
  constexpr uint32_t pulse_pitch_cycles() const {
    return col_cycles(pulse_pitch_cols);
  }
  constexpr uint32_t beacon_pitch_cycles() const {
    return col_cycles(beacon_pitch_cols);
  }
  /// Emission lateness budget (spec §5.2: self-censor past ~½ column).
  constexpr uint32_t late_censor_cycles() const {
    return cycles_per_column() / 2;
  }
  constexpr uint32_t acquire_quiet_cycles() const {
    return col_cycles(acquire_quiet_cols);
  }
  constexpr uint32_t interdigit_timeout_cycles() const {
    return col_cycles(beacon_interdigit_timeout_cols);
  }

  /// Boot-time sanity for the driver's HS_CHECK (e.g. gate_cols < W/4 is what
  /// lets the gate's distance check subsume the boundary-identity check; see
  /// Flywheel::snap).
  constexpr bool valid() const {
    return W > 0 && W % 2 == 0 && cycles_per_half_rev > 0 && gate_cols > 0 &&
           gate_cols < W / 4 && reject_fallback > 0 && pulse_pitch_cols > 0 &&
           gap_timeout_cols > pulse_pitch_cols && effect_count > 0 &&
           effect_count <= 64 && commit_revs > 0 &&
           refractory_revs >
               commit_revs + static_cast<uint32_t>(epoch_repeats) &&
           revs_per_effect > refractory_revs &&
           beacon_period_revs > commit_revs && join_grid_revs > 0 &&
           (64u % join_grid_revs) == 0;
  }
};

/// Phantasm constants: cycles_per_half_rev = cpu_hz · 30 / rpm (exact for
/// 600 MHz / 480 RPM → 37,500,000), glitch filter 100 µs.
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

enum class Boundary : uint8_t { NONE, ZERO, HALF };

enum class Symbol : uint8_t { INVALID, HALF, ZERO, ZERO_EPOCH };

/// Count-coded classification: odd-only, distance-2 alphabet. Any other count
/// (a lost or spurious edge lands on an even value) is INVALID and the whole
/// burst is discarded — fail to "missed", never to "wrong".
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

constexpr Boundary opposite(Boundary b) {
  return b == Boundary::ZERO ? Boundary::HALF : Boundary::ZERO;
}

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
  Boundary last_flipped = Boundary::NONE;

  /// True exactly once per distinct boundary, across both paths.
  bool try_flip(Boundary b) {
    if (b == Boundary::NONE || b == last_flipped)
      return false;
    last_flipped = b;
    return true;
  }
};

// ── Sync-wire edge mailbox (spec §8: single-writer handoff) ─────────────────

/// Consumer's copy of a terminated burst.
struct BurstSnapshot {
  uint32_t count = 0;        ///< Rising edges (after glitch filter).
  uint32_t first_cycles = 0; ///< Timestamp of the burst's first edge.
  uint32_t last_cycles = 0;  ///< Timestamp of the burst's last edge.
};

/**
 * @brief The only state the sync-wire edge ISR writes. The publisher applies
 * the glitch filter and accumulates edge count + first/last timestamps; the
 * consumer (flywheel ISR) detects burst termination by gap timeout and claims
 * the burst. On the device the burst_complete()+claim() pair runs under a
 * brief IRQ-off window (a two-instruction copy, spec §8.2); on the host the
 * calls are plain.
 */
class EdgeMailbox {
public:
  /// Publisher (edge ISR): record one rising edge at @p now.
  void on_edge(uint32_t now, uint32_t glitch_filter_cycles) {
    if (have_prior_ && (now - prior_cycles_) < glitch_filter_cycles)
      return; // EMI spike: too close to the previous accepted edge
    have_prior_ = true;
    prior_cycles_ = now;
    if (count_ == 0)
      first_cycles_ = now;
    last_cycles_ = now;
    if (count_ < 255)
      ++count_; // saturate; anything > 8 is invalid everywhere downstream
  }

  /// Consumer: a burst exists and the wire has been quiet past the timeout.
  bool burst_complete(uint32_t now, uint32_t gap_timeout_cycles) const {
    return count_ > 0 && (now - last_cycles_) >= gap_timeout_cycles;
  }

  /// Consumer: take the burst and reset for the next one.
  BurstSnapshot claim() {
    const BurstSnapshot s{count_, first_cycles_, last_cycles_};
    count_ = 0;
    return s;
  }

private:
  uint32_t count_ = 0;
  uint32_t first_cycles_ = 0;
  uint32_t last_cycles_ = 0;
  uint32_t prior_cycles_ = 0; ///< Last accepted edge ever (glitch filter ref).
  bool have_prior_ = false;
};

// ── Telemetry (spec §8.6) ───────────────────────────────────────────────────

/// Counters maintained by the flywheel ISR, polled by the foreground behind
/// hs::debug. Degradation the protocol absorbs silently must still be
/// visible in one glance of debug output.
struct Telemetry {
  uint32_t symbols_accepted = 0;
  uint32_t symbols_rejected_gate = 0;     ///< §5.3 plausibility rejections.
  uint32_t symbols_discarded_invalid = 0; ///< Invalid pulse counts (§5.2).
  uint32_t beacons_ok = 0;
  uint32_t beacons_rejected = 0;          ///< Checksum/digit-count failures.
  uint32_t beacon_index_corrections = 0;  ///< Missed-epoch fixes (§6.3.2).
  uint32_t beacon_rev_mismatches = 0;
  uint32_t epochs_refractory_ignored = 0;
  uint32_t lock_transitions = 0; ///< ACQUIRE↔LOCKED edges.
  uint32_t flips = 0;
  uint32_t emit_censored = 0;    ///< Master skipped a late boundary symbol.
  uint32_t emit_aborted = 0;     ///< Master stopped emitting mid-burst.
  uint32_t max_coast_halves = 0; ///< Longest run of half-revs without a snap.
};

// ── Layer 1: the flywheel timebase (spec §4) ────────────────────────────────

enum class LockState : uint8_t { ACQUIRE, LOCKED };

/// A locally-crossed boundary, reported by Flywheel::fold().
struct Crossing {
  bool crossed = false;
  Boundary boundary = Boundary::NONE;
  uint32_t at_cycles = 0; ///< The exact (timebase) instant of the boundary.
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
  explicit Flywheel(const Config &cfg)
      : period_(cfg.cycles_per_half_rev), w_(cfg.W), gate_cols_(cfg.gate_cols),
        reject_fallback_(cfg.reject_fallback) {}

  /// Boot seeding (spec §8.5): epoch = now, boundary ZERO, ACQUIRE.
  void seed(uint32_t now) {
    epoch_cycles_ = now;
    boundary_ = Boundary::ZERO;
    lock_ = LockState::ACQUIRE;
    consecutive_rejects_ = 0;
  }

  /// Master is the reference by definition — it never snaps and is born
  /// locked (spec §2: "Master's flywheel IS the reference").
  void force_lock() { lock_ = LockState::LOCKED; }

  /// Current column at time @p at. Signed-safe for timestamps up to ~±3.5 s
  /// around the epoch (snap evaluation may look slightly into the past).
  int32_t position(uint32_t at) const {
    const int64_t delta =
        static_cast<int32_t>(at - epoch_cycles_); // modular → signed
    const int64_t cols = floor_div(delta * (w_ / 2), period_);
    return floor_mod(boundary_column(boundary_, w_) + cols, w_);
  }

  /**
   * @brief The rebase rule: if @p now has passed the next boundary, fold the
   * epoch forward by exactly one half-rev (integer add — no drift) and report
   * the crossing. Call in a loop until it returns crossed=false; a long coast
   * (masked window) yields several crossings, each at its exact instant.
   */
  Crossing fold(uint32_t now) {
    const int32_t delta = static_cast<int32_t>(now - epoch_cycles_);
    if (delta < 0 || static_cast<uint32_t>(delta) < period_)
      return {};
    epoch_cycles_ += period_;
    boundary_ = opposite(boundary_);
    return {true, boundary_, epoch_cycles_};
  }

  enum class SnapOutcome : uint8_t { ACCEPTED, REJECTED, REJECTED_FELL_BACK };

  /**
   * @brief Acceptance gate + re-base (spec §4.2, §5.3).
   *
   * LOCKED: accept only if the implied correction is ≤ G columns. Because
   * G < W/4 (Config::valid), passing the distance gate also proves the named
   * boundary is the flywheel's nearest predicted boundary — the identity
   * check is subsumed. R consecutive rejections fall back to ACQUIRE so the
   * gate can never deadlock a genuinely-lost board.
   *
   * ACQUIRE: hard snap, no gate (the SyncBoard applies the quiet-before
   * routing guard before calling this).
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
   * @brief Count one implausible-symbol rejection toward the ACQUIRE
   * fallback (shared by the snap gate and the SyncBoard's suspect-burst
   * timeout). Returns true when R consecutive rejections concluded this
   * board's own timebase — not the wire — is at fault (spec §5.3: the
   * fallback is mandatory; a gate without an escape deadlocks a lost board
   * into rejecting good symbols forever).
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

  LockState lock() const { return lock_; }
  Boundary current_boundary() const { return boundary_; }
  uint32_t epoch_cycles() const { return epoch_cycles_; }
  uint32_t cycles_per_half_rev() const { return period_; }
  /// §4.3 frequency trim hook (snap-only ships; tests exercise extremes).
  void set_cycles_per_half_rev(uint32_t c) { period_ = c; }

private:
  uint32_t period_;       ///< cycles_per_half_rev (optionally trimmed).
  int32_t w_;
  int32_t gate_cols_;
  int32_t reject_fallback_;
  uint32_t epoch_cycles_ = 0;
  Boundary boundary_ = Boundary::ZERO; ///< Column identity at epoch_cycles_.
  LockState lock_ = LockState::ACQUIRE;
  int32_t consecutive_rejects_ = 0;
};

// ── Layer 3: beacon codec (spec §6.4) ───────────────────────────────────────

/// Decoded index beacon: absolute effect index + revolution count (mod 64).
struct BeaconFrame {
  int32_t effect_index = 0; ///< 0–63.
  uint32_t rev_count = 0;   ///< Revolution within the effect, mod 64.
};

/// Five base-8 digits [index_hi, index_lo, rev_hi, rev_lo, checksum], each
/// transmitted as a burst of (digit+1) pulses.
constexpr void encode_beacon_digits(int32_t effect_index, uint32_t rev_count,
                                    uint8_t out[5]) {
  const uint32_t idx = static_cast<uint32_t>(effect_index) & 63u;
  const uint32_t rev = rev_count & 63u;
  out[0] = static_cast<uint8_t>(idx >> 3);
  out[1] = static_cast<uint8_t>(idx & 7u);
  out[2] = static_cast<uint8_t>(rev >> 3);
  out[3] = static_cast<uint8_t>(rev & 7u);
  out[4] = static_cast<uint8_t>((out[0] + out[1] + out[2] + out[3]) & 7u);
}

/**
 * @brief Assembles data bursts into beacon frames. Integrity by rejection
 * (spec §6.4): any out-of-range digit, stale partial frame, or checksum
 * mismatch drops the whole frame — the next beacon is ≤ 2 s away.
 */
class BeaconParser {
public:
  /// Feed one data burst. Returns true when a complete, checksum-valid frame
  /// was parsed into @p out. @p rejected is set when a frame (not merely a
  /// digit) was discarded.
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
    if (((digits_[0] + digits_[1] + digits_[2] + digits_[3]) & 7) !=
        digits_[4]) {
      *rejected = true;
      return false;
    }
    out->effect_index = digits_[0] * 8 + digits_[1];
    out->rev_count = static_cast<uint32_t>(digits_[2]) * 8u + digits_[3];
    return true;
  }

  void reset() { n_ = 0; }
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
  bool identity_known = false; ///< False until an epoch/beacon names the index.
  int32_t effect_index = 0;
  uint32_t rev_in_effect = 0; ///< ZERO crossings since effect start. For a
                              ///< beacon-joined board this starts from the
                              ///< beacon's mod-64 value — it only feeds the
                              ///< master's schedule and cross-check telemetry.
  bool commit_pending = false;
  uint32_t commit_in_revs = 0;    ///< ZERO crossings until the B+R+K commit.
  uint32_t refractory_revs_left = 0; ///< EPOCH dedup window (spec §6.1).

  /// Accepted EPOCH symbol. Returns true if it opened a commit window (false
  /// inside the refractory window — the §6.3 redundancy repeats land here).
  ///
  /// The commit boundary is ABSOLUTE: B + R + K, where B is the primary
  /// copy's boundary (R = epoch_repeats, K = commit_revs). Which copy of the
  /// train this is (j, 0 = primary) is inferred from the shared revolution
  /// count — the master starts the train exactly when rev_in_effect reaches
  /// revs_per_effect, and by the time a symbol is consumed the local crossing
  /// for its boundary has already incremented rev_in_effect (classification
  /// completes ~13 columns after the boundary instant). So every board that
  /// hears ANY copy counts down to the same boundary, and hearing a repeat
  /// instead of the primary cannot skew the commit (§6.3.1). A board whose
  /// revolution count is not absolute (it beacon-joined mid-effect, §6.4)
  /// lands outside the train window and falls back to j = 0 — it commits up
  /// to j revolutions late, the pre-existing epoch-bounded degradation, now
  /// confined to that case (and erased at its first commit, which re-zeros
  /// rev_in_effect in step with the master's).
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

  /// True at the single crossing/accept where the construction window — the
  /// last K revolutions before the commit boundary — opens (commit_in_revs
  /// strictly decreases, so == fires exactly once per window). The announce
  /// phase before it is what gives a board that hears only the last repeat
  /// the full K-revolution construction budget; the outgoing effect keeps
  /// playing through it, so the dark window stays K revolutions.
  bool construction_opens(const Config &cfg) const {
    return commit_pending && commit_in_revs == cfg.commit_revs;
  }

  /// True throughout the construction window (fail-dark, spec §6.1).
  bool constructing(const Config &cfg) const {
    return commit_pending && commit_in_revs <= cfg.commit_revs;
  }

  /// A ZERO boundary flip happened. Returns true when the commit deadline is
  /// reached — the board must swap to the (already constructed) next effect
  /// at exactly this boundary.
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
  /// Master crossed boundary @p at_cycles and wants to emit @p symbol.
  /// @returns false if the symbol was self-censored (caller counts it).
  bool schedule_boundary(Symbol symbol, uint32_t at_cycles, uint32_t now,
                         const Config &cfg) {
    if (symbol == Symbol::INVALID)
      return false;
    if ((now - at_cycles) > cfg.late_censor_cycles())
      return false; // late at the boundary: skip the whole symbol
    // The wire must be idle: a still-running beacon frame here means the
    // schedule violated its own spacing — drop the boundary symbol rather
    // than corrupt both (cannot happen with in-range config; defensive).
    if (pulses_left_ > 0 || queue_pos_ < queue_len_)
      return false;
    pulses_left_ = symbol_pulse_count(symbol);
    next_due_ = at_cycles;
    pitch_ = cfg.pulse_pitch_cycles();
    return true;
  }

  /// Master reached the beacon point (x ≈ W/4): queue the five digit bursts.
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

  /// Per-tick: should the ISR emit one pulse right now? @p aborted reports a
  /// mid-burst self-censor (telemetry).
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

  bool idle() const { return pulses_left_ == 0 && queue_pos_ >= queue_len_; }

private:
  struct PendingBurst {
    uint32_t start_cycles = 0;
    uint32_t pulses = 0;
  };
  PendingBurst queue_[5] = {};
  int32_t queue_len_ = 0;
  int32_t queue_pos_ = 0;
  uint32_t pulses_left_ = 0;
  uint32_t next_due_ = 0;
  uint32_t pitch_ = 0;
};

// ── The per-board engine ────────────────────────────────────────────────────

/// What the device ISR (or host sim) must do after one flywheel tick.
struct TickActions {
  int32_t render_column = -1; ///< Column to pack/display, or -1 (unchanged
                              ///< since last tick, or display is dark).
  bool dark = false;          ///< Display black (ACQUIRE / no identity /
                              ///< epoch construction window, spec §6.1).
  bool flip = false;          ///< Call advance_display() on the live effect.
  bool zero_crossing = false; ///< A ZERO boundary flip happened this tick.
  bool join_boundary = false; ///< ZERO crossing on the join grid: a board
                              ///< with no live effect takes its pending one
                              ///< here (Config::join_grid_revs).
  bool commit = false;        ///< B+R+K deadline: swap to the pending effect
                              ///< at this boundary (trap if it is not ready).
  bool pulse = false;         ///< Master: emit one sync pulse (pin-first).
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
  explicit SyncBoard(const Config &cfg)
      : cfg_(cfg), fly_(cfg) {}

  /// Boot/reboot seeding (spec §8.5). Master is born locked with identity
  /// (effect 0, rev 0) — it is the reference. Downstream boards start in
  /// ACQUIRE, dark, and join via boundary symbols + beacon.
  void seed(uint32_t now, bool is_master) {
    is_master_ = is_master;
    fly_.seed(now);
    gate_ = FlipGate{};
    content_ = ContentTracker{};
    telemetry_ = Telemetry{};
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

  /// Sync-wire edge ISR entry point (publisher; downstream boards only).
  void on_sync_edge(uint32_t now) {
    mailbox_.on_edge(now, cfg_.glitch_filter_cycles);
  }

  /// Device-side helper for the IRQ-off mailbox handoff.
  EdgeMailbox &mailbox() { return mailbox_; }
  uint32_t gap_timeout_cycles() const { return cfg_.gap_timeout_cycles(); }

  /**
   * @brief One flywheel wake-up. @p burst is the claimed mailbox burst if one
   * completed (see EdgeMailbox), else nullptr.
   */
  TickActions tick(uint32_t now, const BurstSnapshot *burst) {
    TickActions a{};
    if (burst)
      handle_burst(*burst, a);

    // Suspect-burst timeout: a lone valid-count burst far from every
    // predicted boundary was held pending (handle_burst). If the beacon
    // interdigit window expired with no follow-up burst, it was not beacon
    // data — count it as a gate rejection so a board with a corrupted
    // timebase (whose REAL boundary symbols all land "far") still reaches
    // the §5.3 ACQUIRE fallback instead of deadlocking.
    if (suspect_pending_ &&
        (now - suspect_last_cycles_) > cfg_.interdigit_timeout_cycles() &&
        static_cast<int32_t>(now - suspect_last_cycles_) > 0) {
      suspect_pending_ = false;
      ++telemetry_.symbols_rejected_gate;
      if (fly_.note_rejection())
        ++telemetry_.lock_transitions;
    }

    // Fold every locally-crossed boundary (usually 0 or 1; several after a
    // long masked coast — the skipped columns were undisplayable anyway).
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

    // Render decision. Dark whenever phase or content identity is missing,
    // and during the epoch construction window — the last K revolutions of
    // the commit countdown, an absolute boundary all boards share, so the
    // window is deterministic on every board (spec §6.1 — never a stale
    // frame on some and black on others). The announce phase before it
    // keeps playing the outgoing effect.
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

  /// Build request: (generation << 8) | effect_index. The foreground
  /// constructs the named effect whenever the generation changes, then
  /// publishes it to the driver's pending slot.
  uint32_t build_word() const { return build_word_; }
  static int32_t build_index_of(uint32_t word) {
    return static_cast<int32_t>(word & 0xFF);
  }
  static uint32_t build_gen_of(uint32_t word) { return word >> 8; }

  const Telemetry &telemetry() const { return telemetry_; }
  const ContentTracker &content() const { return content_; }
  LockState lock() const { return fly_.lock(); }
  const Flywheel &flywheel() const { return fly_; }
  Flywheel &flywheel_mut() { return fly_; } // tests only
  const Config &config() const { return cfg_; }

private:
  // ── Flip + content events (both paths funnel here) ───────────────────────

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
        // Last K revolutions: ask the foreground to construct the next
        // effect now, on every board at the same absolute boundary.
        publish_build((content_.effect_index + 1) % cfg_.effect_count);
      else if (!content_.commit_pending &&
               (content_.rev_in_effect % cfg_.join_grid_revs) == 0)
        a.join_boundary = true;
    }
  }

  // ── Receive path (downstream) ─────────────────────────────────────────────

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
        // A lone valid-count burst far from every predicted boundary is
        // either EMI or — if this board's timebase is corrupted — a real
        // boundary symbol it can no longer place. It is indistinguishable
        // from a beacon's first digit until we know whether a train
        // follows, so hold it as a suspect; the timeout in tick() converts
        // it into a gate rejection if the wire stays silent.
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
      // ACQUIRE quiet-before guard: boundary symbols are isolated on the
      // wire (≥ half-rev apart); a burst following close on another is a
      // beacon digit train. Keeps a just-rebooted board from hard-snapping
      // to a mid-revolution data burst.
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
    apply_flip(b, a); // backstop flip; deduped if the crossing already fired
    if (sym == Symbol::ZERO_EPOCH && content_.identity_known) {
      if (content_.on_epoch_symbol(cfg_)) {
        // Construction normally starts at a later crossing (apply_flip);
        // a board that heard only the LAST repeat opens the window at the
        // accept itself.
        if (content_.construction_opens(cfg_))
          publish_build((content_.effect_index + 1) % cfg_.effect_count);
      } else {
        ++telemetry_.epochs_refractory_ignored;
      }
    }
  }

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
      // Inside the epoch window the displayed index is in flux; the next
      // post-commit beacon re-verifies.
    } else if (idx != content_.effect_index) {
      // Missed epoch (all repeats): correct within ≤16 revs (spec §6.3.2).
      content_.effect_index = idx;
      content_.rev_in_effect = f.rev_count;
      ++telemetry_.beacon_index_corrections;
      publish_build(idx);
    } else if (f.rev_count != (content_.rev_in_effect & 63u)) {
      ++telemetry_.beacon_rev_mismatches; // §6.2: detect, don't retro-correct
    }
  }

  // ── Conduct + emit path (master) ──────────────────────────────────────────

  void master_on_crossing(const Crossing &c, uint32_t now, TickActions &a) {
    Symbol sym = Symbol::HALF;
    if (c.boundary == Boundary::ZERO) {
      // Conductor (spec §6.1): when the effect's revolutions elapse, start
      // an EPOCH train — the primary copy plus R redundancy repeats on the
      // following ZERO boundaries (idempotent via the refractory window).
      if (epoch_emits_left_ == 0 && !content_.commit_pending &&
          content_.rev_in_effect >= cfg_.revs_per_effect) {
        epoch_emits_left_ = 1 + cfg_.epoch_repeats;
        // j = 0 by construction (the train starts at this crossing); the
        // build publish follows when the construction window opens —
        // immediately, if epoch_repeats is zero.
        if (content_.on_epoch_symbol(cfg_) &&
            content_.construction_opens(cfg_))
          publish_build((content_.effect_index + 1) % cfg_.effect_count);
      }
      if (epoch_emits_left_ > 0) {
        sym = Symbol::ZERO_EPOCH;
        --epoch_emits_left_;
      } else {
        sym = Symbol::ZERO;
      }
    }
    if (!emitter_.schedule_boundary(sym, c.at_cycles, now, cfg_))
      ++telemetry_.emit_censored;
  }

  void maybe_schedule_beacon(uint32_t now) {
    // Beacon point: x ≈ W/4, i.e. mid-way through the ZERO→HALF half-rev,
    // where the wire is otherwise quiet (spec §6.4). Emitted on rev 1 of
    // every beacon period plus the first revs of a fresh effect (confirming
    // the post-commit index immediately) — never rev 0 of boot, so a
    // just-powered downstream board sees clean boundary symbols first.
    if (beacon_done_this_rev_ || fly_.current_boundary() != Boundary::ZERO)
      return;
    // Stay silent during the commit window: a beacon here would broadcast
    // the outgoing effect's index, and a board joining off it would adopt
    // stale identity just as everyone else switches.
    if (content_.commit_pending)
      return;
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
  bool suspect_pending_ = false;     ///< Lone far burst awaiting train/timeout.
  uint32_t suspect_last_cycles_ = 0;
  uint32_t epoch_emits_left_ = 0;
  bool beacon_done_this_rev_ = false;
  uint32_t build_gen_ = 0;
  uint32_t build_word_ = 0; ///< (gen << 8) | index; foreground-read.
};

} // namespace sync
} // namespace pov
