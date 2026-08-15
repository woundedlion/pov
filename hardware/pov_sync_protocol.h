/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file pov_sync_protocol.h
 * @brief Shared vocabulary of the Phantasm sync protocol: signed-correct
 *        ring arithmetic, the Config constant block, the count-coded symbol
 *        alphabet, the boundary flip gate, the edge-ISR mailbox and the
 *        telemetry counters.
 *
 * Every other pov_sync header builds on this one; see pov_sync.h for the
 * architecture the pieces assemble into.
 */
#pragma once

#include <cstdint>

// Forward declaration of the unit-test accessor that reaches EdgeMailbox's
// split consumer path (burst_complete()/claim()).
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
  int32_t W = 288; /**< Columns per revolution. */
  uint32_t cycles_per_half_rev =
      37500000;                          /**< Timebase constant (spec §4.1). */
  uint32_t glitch_filter_cycles = 60000; /**< Min edge spacing (~100 µs). */

  // Symbol wire (spec §5.2): pitches/timeouts in columns.
  int32_t pulse_pitch_cols = 2; /**< Boundary-burst pulse pitch (> mask M). */
  // beacon_pitch_cols, gap_timeout_cols and acquire_quiet_cols all feed
  // beacon_frame_cols(), which valid() holds strictly under W/4 with a single
  // column to spare at the shipped constants. Any of the three can only move
  // down unless W grows with it.
  int32_t beacon_pitch_cols = 1; /**< Beacon digit pulse pitch (checksummed). */
  int32_t gap_timeout_cols = 4;  /**< Quiet time that terminates a burst. */

  // Acceptance gate (spec §5.3).
  int32_t gate_cols = 4; /**< G: max plausible snap correction (LOCKED). */
  int32_t reject_fallback = 4; /**< R: consecutive rejections before ACQUIRE. */
  int32_t acquire_quiet_cols = 16; /**< Quiet-before guard for ACQUIRE snaps. */

  // Content layer (spec §6).
  uint32_t revs_per_effect =
      960; /**< Effect duration in revolutions (120 s). */
  /**
   * @brief Optional per-roster-entry effect durations, in revolutions.
   * @details CONTRACT — when non-null this must point at at least
   * effect_count entries: revolutions_for_effect() indexes it unguarded, and
   * valid() itself sweeps [0, effect_count), so a short table reads out of
   * bounds inside the validator. No length travels with the pointer, so hold
   * the two together at the definition site (both sketches size the table off
   * the same roster macro and static_assert the match).
   */
  const uint32_t *effect_revolutions = nullptr;
  int32_t epoch_repeats = 3;     /**< EPOCH redundancy repeats (spec §6.3). */
  uint32_t refractory_revs = 16; /**< EPOCH dedup window (spec §6.1). */
  /**
   * @brief K: the construction window, the last K revolutions before the
   * commit.
   * @details The commit boundary itself is B + epoch_repeats + commit_revs from
   * the train's primary copy at B, so a board hearing any repeat still commits
   * in lockstep with the full K-rev construction budget (spec §6.1, §6.3.1).
   * The foreground polls the build request between frames, so the budget an
   * effect actually gets is K revolutions less the render in flight when the
   * request was published; overrunning it trips the driver's commit_ok trap on
   * every board at once. Not checkable in valid(), which sees no render times.
   */
  uint32_t commit_revs = 2;
  uint32_t beacon_period_revs = 16; /**< Beacon cadence (spec §6.4). */
  int32_t beacon_interdigit_timeout_cols = 24; /**< Stale-frame reset bound. */
  int32_t effect_count = 1; /**< Roster length (epoch wraps mod this). */
  /**
   * @brief Live-takeover grid: boards take a constructed effect live only at
   * revolutions ≡ 0 mod this.
   * @details At boot every board (master included) therefore goes live at the
   * SAME crossing with frame counters aligned, instead of master leading by
   * however long downstream identity took. Must divide 64 so a beacon's
   * mod-64 revolution count lands on the same grid as the master's true
   * count. Mid-show rejoins wait ≤ grid revolutions past the beacon that named
   * their effect (spec §9.1, a term of rejoin_bound_revs()).
   */
  uint32_t join_grid_revs = 4;

  /**
   * @brief Rejoin budget: the most revolutions a mid-show board may wait to go
   * live on the right effect (spec §9.1).
   * @details The ceiling on rejoin_bound_revs(), enforced in valid() rather
   * than left to prose. Expressed in revolutions so the bound is rotation-rate
   * independent; 25 revolutions is ~3.1 s at the nominal 480 RPM.
   */
  uint32_t rejoin_budget_revs = 25;

  /**
   * @brief Returns the configured duration for one roster entry.
   * @param effect_index Roster index in [0, effect_count).
   * @return Duration in revolutions.
   * @details Indexes effect_revolutions unguarded — see its length contract.
   */
  constexpr uint32_t revolutions_for_effect(int32_t effect_index) const {
    return effect_revolutions ? effect_revolutions[effect_index]
                              : revs_per_effect;
  }

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
   * @brief Ceiling on a burst's total duration, in cycles.
   * @return Widest legitimate burst span, plus the gap timeout and a column of
   * slack.
   * @details Sustained wire noise inside the glitch filter's pass band refreshes
   * a burst's last edge indefinitely, so the terminating gap never opens. This
   * bound terminates such a burst anyway; it sits above every legitimate span —
   * five boundary pulses at pulse_pitch_cols, eight beacon pulses at
   * beacon_pitch_cols — so a real burst always ends on the gap first.
   */
  constexpr uint32_t max_burst_cycles() const {
    const int32_t span = 4 * pulse_pitch_cols > 7 * beacon_pitch_cols
                             ? 4 * pulse_pitch_cols
                             : 7 * beacon_pitch_cols;
    return col_cycles(span + gap_timeout_cols + 1);
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
   * @brief Beacon frame span (first to last pulse) for one payload, in columns.
   * @param digit_sum Sum of the five base-8 digits (0–35).
   * @return Columns the frame's pulses span.
   * @details Mirrors schedule_beacon's per-digit advance: each digit burst
   * spans digit pitches, and four inter-burst gaps of gap_timeout_cols + 1
   * separate the five bursts.
   */
  constexpr int32_t beacon_span_cols(int32_t digit_sum) const {
    return digit_sum * beacon_pitch_cols + 4 * (gap_timeout_cols + 1);
  }
  /**
   * @brief Worst-case beacon frame span (first to last pulse), in columns.
   * @return Columns spanned by five base-8 digit bursts of value 7.
   */
  constexpr int32_t beacon_span_cols() const { return beacon_span_cols(35); }
  /**
   * @brief Columns a beacon frame occupies from its first pulse to the earliest
   * instant a boundary burst may follow it.
   * @param digit_sum Sum of the five base-8 digits (0–35).
   * @return Pulse span plus the quiet a receiver needs after the last pulse.
   * @details The quiet term is acquire_quiet_cols, the wider of the two windows
   * the tail must clear: gap_timeout_cols terminates the last digit burst, and
   * acquire_quiet_cols is what makes the following boundary burst read as an
   * isolated symbol rather than another digit. valid()'s demarcation relation
   * (acquire_quiet_cols >= beacon_span_cols() / 4) puts it above the gap
   * timeout for every pitch.
   */
  constexpr int32_t beacon_frame_cols(int32_t digit_sum) const {
    return beacon_span_cols(digit_sum) + acquire_quiet_cols;
  }
  /**
   * @brief Worst-case columns a beacon frame occupies.
   * @return beacon_frame_cols() for the widest encodable payload.
   */
  constexpr int32_t beacon_frame_cols() const { return beacon_frame_cols(35); }

  /**
   * @brief Worst-case revolutions from a mid-show join to going live on the
   * right effect (spec §9.1).
   * @return Widest beacon-to-beacon gap plus the join-grid wait.
   * @details Beacons are suppressed for the whole commit window, so the widest
   * gap between two beacons is one cadence plus that window (epoch_repeats
   * announce revolutions + commit_revs), not the cadence alone; a board that
   * joins just after a beacon waits that gap, then up to join_grid_revs more
   * for the next live-takeover boundary. Requires epoch_repeats >= 0.
   */
  constexpr uint32_t rejoin_bound_revs() const {
    return beacon_period_revs + join_grid_revs +
           static_cast<uint32_t>(epoch_repeats) + commit_revs;
  }

  /**
   * @brief Boot-time sanity check for the driver's HS_CHECK.
   * @return nullptr if every protocol constant is self-consistent, otherwise a
   * literal naming the first relation that failed.
   * @details gate_cols < W/4 is what lets the gate's distance check subsume the
   * boundary-identity check; see Flywheel::snap. Relations are tested in an
   * order that makes each one's preconditions (nonzero divisors, non-negative
   * casts) already established.
   */
  HS_FLASH_MEMBER constexpr const char *valid() const {
    if (!(W > 0))
      return "W > 0";
    if (!(W % 2 == 0))
      return "W even";
    if (!(cycles_per_half_rev > 0))
      return "cycles_per_half_rev > 0";
    if (!(gate_cols > 0))
      return "gate_cols > 0";
    if (!(gate_cols < W / 4))
      return "gate_cols < W/4";
    if (!(reject_fallback > 0))
      return "reject_fallback > 0";
    if (!(glitch_filter_cycles > 0))
      return "glitch_filter_cycles > 0";
    if (!(pulse_pitch_cols > 0))
      return "pulse_pitch_cols > 0";
    if (!(gap_timeout_cols > pulse_pitch_cols))
      return "gap_timeout_cols > pulse_pitch_cols";
    if (!(beacon_pitch_cols > 0))
      return "beacon_pitch_cols > 0";
    if (!(gap_timeout_cols > beacon_pitch_cols))
      return "gap_timeout_cols > beacon_pitch_cols";
    // Both pitches must clear the glitch filter. A pulse may be emitted up to
    // late_censor_cycles() late, compressing its gap to the next on-time pulse;
    // a filter wider than what remains swallows every pulse after the first, so
    // the burst decodes as Symbol::HALF — the miscount the odd-only alphabet
    // exists to prevent.
    if (!(glitch_filter_cycles < pulse_pitch_cycles() - late_censor_cycles()))
      return "glitch_filter_cycles < pulse_pitch - late_censor";
    if (!(glitch_filter_cycles < beacon_pitch_cycles() - late_censor_cycles()))
      return "glitch_filter_cycles < beacon_pitch - late_censor";
    // Demarcation headroom: a beacon frame starts at W/4, so its last pulse
    // comes within W/4 - beacon_span_cols() of the HALF boundary. A gate radius
    // above that separation makes handle_burst() claim beacon digits as
    // boundary symbols and snap on them; the frame and quiet clauses below keep
    // 7*beacon_pitch_cols — the widest single digit burst — under it.
    if (!(7 * beacon_pitch_cols + 1 > gate_cols))
      return "7*beacon_pitch_cols + 1 > gate_cols";
    // maybe_schedule_beacon emits only in [W/4, W/2), so the worst-case frame
    // plus its tail quiet must clear W/4 or no beacon is ever scheduled.
    // Strict: the slack absorbs the emitter's ½-column lateness budget, which
    // the scheduling fit charges to the frame, plus the sub-column offset
    // between the W/4 instant and the tick that schedules it.
    // Achieved margin at the shipped constants (W = 288): 71 of 72 columns —
    // one column, half of it the lateness budget. The tightest relation in this
    // function.
    if (!(beacon_frame_cols() < W / 4))
      return "beacon_frame_cols() < W/4";
    // Demarcation: the acquisition timeout must clear the beacon's worst-case
    // per-digit advance, which beacon_span_cols() / 4 bounds. The strict
    // ordering below makes the interdigit timeout larger.
    if (!(acquire_quiet_cols >= beacon_span_cols() / 4))
      return "acquire_quiet_cols >= beacon_span_cols()/4";
    // Stale-frame window order: tick()'s poll-path reset must be the tighter
    // one, so a truncated train drops on wire silence rather than waiting for
    // the next burst to reach BeaconParser::feed.
    if (!(acquire_quiet_cols < beacon_interdigit_timeout_cols))
      return "acquire_quiet_cols < beacon_interdigit_timeout_cols";
    if (!(effect_count > 0))
      return "effect_count > 0";
    if (!(effect_count <= 64))
      return "effect_count <= 64";
    if (!(commit_revs > 0))
      return "commit_revs > 0";
    // Gate epoch_repeats >= 0 first: a negative value casts to a huge uint32_t
    // and wraps the refractory bound below.
    if (!(epoch_repeats >= 0))
      return "epoch_repeats >= 0";
    if (!(refractory_revs > commit_revs + static_cast<uint32_t>(epoch_repeats)))
      return "refractory_revs > commit_revs + epoch_repeats";
    for (int32_t i = 0; i < effect_count; ++i)
      if (!(revolutions_for_effect(i) > refractory_revs))
        return "revolutions_for_effect(i) > refractory_revs";
    if (!(beacon_period_revs > commit_revs))
      return "beacon_period_revs > commit_revs";
    // Beacon rev resync recovers a slip only in (-32, +32), so keep the period
    // below the half-window.
    if (!(beacon_period_revs < 32))
      return "beacon_period_revs < 32";
    // §9.1 rejoin budget: cap the achieved bound, not the cadence alone — the
    // commit window's beacon blackout and the join grid are part of what a
    // rejoiner waits through.
    if (!(rejoin_bound_revs() <= rejoin_budget_revs))
      return "rejoin_bound_revs() <= rejoin_budget_revs";
    if (!(join_grid_revs > 0))
      return "join_grid_revs > 0";
    if (!((64u % join_grid_revs) == 0))
      return "join_grid_revs divides 64";
    // schedule_beacon's "is-due" check reads (now - start_cycles) as int32, so
    // the worst-case span (5 digits of value 7) must clear 2^31.
    if (!(5u * (7u * beacon_pitch_cycles() +
                static_cast<uint32_t>(gap_timeout_cols + 1) *
                    cycles_per_column()) <
          static_cast<uint32_t>(INT32_MAX)))
      return "worst-case beacon span < INT32_MAX";
    return nullptr;
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
 * consumer (flywheel ISR) detects burst termination by gap timeout or duration
 * bound and claims the burst with try_claim(), which tests completion and takes
 * the burst in one
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
    if (have_prior && (now - prior_cycles) < glitch_filter_cycles)
      return; // EMI spike: too close to the previous accepted edge
    have_prior = true;
    prior_cycles = now;
    if (count == 0)
      first_cycles = now;
    last_cycles = now;
    if (count < 255)
      ++count; // saturate
  }

  /**
   * @brief Consumer: atomically test for a terminated burst and, if present,
   * take it.
   * @param now Current timestamp, in cycles.
   * @param gap_timeout_cycles Quiet time that terminates a burst, in cycles.
   * @param max_burst_cycles Duration past which an unterminated burst is
   * claimed regardless of the gap (Config::max_burst_cycles).
   * @param[out] out Burst snapshot, written only when true is returned.
   * @return True if a burst had terminated and was claimed.
   * @details The edge ISR must not run between the completion test and the
   * reset; the device brackets this call in IRQ-off.
   * `now` is sampled before the bracket opens, so an edge landing in between
   * leaves `last_cycles` ahead of it; the signed re-check rejects that wrapped
   * modular difference instead of claiming a burst still in flight.
   * The duration term takes the same re-check against `first_cycles`. It is what
   * frees the mailbox under sustained noise, where every accepted edge pushes
   * the gap out and the quiet term alone never fires: the burst is claimed,
   * classified INVALID by its count, and counted.
   */
  bool try_claim(uint32_t now, uint32_t gap_timeout_cycles,
                 uint32_t max_burst_cycles, BurstSnapshot *out) {
    if (count == 0)
      return false;
    const bool quiet = (now - last_cycles) >= gap_timeout_cycles &&
                       static_cast<int32_t>(now - last_cycles) > 0;
    const bool overlong = (now - first_cycles) >= max_burst_cycles &&
                          static_cast<int32_t>(now - first_cycles) > 0;
    if (!quiet && !overlong)
      return false;
    *out = BurstSnapshot{count, first_cycles, last_cycles};
    count = 0;
    return true;
  }

  /**
   * @brief Consumer: retire the glitch-filter reference once the wire has been
   * quiet longer than the filter window.
   * @param now Current timestamp, in cycles.
   * @param glitch_filter_cycles Filter window, in cycles.
   * @details Keeps `now - prior_cycles` bounded: `prior_cycles` otherwise
   * persists indefinitely, and after ~7.16 s of wire silence the cycle counter
   * wraps, making that modular difference pseudo-random — with p ≈ glitch/2³²
   * it lands inside the reject window and falsely rejects a real edge. The
   * flywheel poll calls this every column, so a stale reference is cleared
   * within one column of silence, long before the counter can wrap. Must run
   * under the same single-writer discipline as claim(): it writes have_prior,
   * which the edge ISR also writes. `now` is sampled before the bracket opens,
   * so the signed re-check rejects the wrapped modular difference an edge
   * accepted in between would produce.
   */
  void age_prior(uint32_t now, uint32_t glitch_filter_cycles) {
    if (have_prior && (now - prior_cycles) >= glitch_filter_cycles &&
        static_cast<int32_t>(now - prior_cycles) > 0)
      have_prior = false;
  }

private:
  // Split consumer path, kept private behind a test friend so production takes
  // only the unsplittable try_claim().
  friend struct ::hs_test::pov_sync_tests::EdgeMailboxTestAccess;

  /**
   * @brief Consumer (test-only): has a burst terminated?
   * @param now Current timestamp, in cycles.
   * @param gap_timeout_cycles Quiet time that terminates a burst, in cycles.
   * @return True if a burst exists and the wire has been quiet past the timeout.
   */
  bool burst_complete(uint32_t now, uint32_t gap_timeout_cycles) const {
    return count > 0 && (now - last_cycles) >= gap_timeout_cycles;
  }

  /**
   * @brief Consumer (test-only): take the burst and reset for the next one.
   * @return A snapshot of the terminated burst.
   */
  BurstSnapshot claim() {
    const BurstSnapshot s{count, first_cycles, last_cycles};
    count = 0;
    return s;
  }

  // Ordered against the edge ISR by the consumer's __disable_irq() memory
  // clobber, not by volatile.
  uint32_t count = 0;
  uint32_t first_cycles = 0;
  uint32_t last_cycles = 0;
  uint32_t prior_cycles =
      0; /**< Last accepted edge ever (glitch filter ref). */
  bool have_prior = false;
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
  uint32_t beacons_rejected =
      0; /**< Checksum, digit-count, or staleness drops. */
  uint32_t beacon_index_corrections = 0; /**< Missed-epoch fixes (§6.3.2). */
  uint32_t beacon_rev_mismatches = 0;
  uint32_t epochs_refractory_ignored = 0;
  uint32_t lock_transitions = 0; /**< ACQUIRE↔LOCKED edges. */
  uint32_t flips = 0;
  uint32_t emit_censored = 0; /**< Master skipped a late boundary symbol. */
  uint32_t emit_aborted = 0;  /**< Master truncated a burst mid-emission. */
  uint32_t beacons_busy_dropped =
      0; /**< Revolutions whose beacon schedule found the emitter busy. */
  uint32_t beacons_late_dropped =
      0; /**< Revolutions whose beacon start was too late to fit before HALF. */
  uint32_t beacons_overrun_dropped =
      0; /**< Stale overrun beacon dropped at a boundary. */
  uint32_t boundary_bursts_dropped =
      0; /**< Undrained boundary symbol dropped at the next boundary. */
  uint32_t max_coast_halves =
      0; /**< Longest run of half-revs without a snap. */
  uint32_t master_stalls =
      0; /**< Master coasts past the signed-safe window, re-seeded. */
};

} // namespace sync
} // namespace pov
