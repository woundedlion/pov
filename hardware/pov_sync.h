/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file pov_sync.h
 * @brief Pure, host-testable core of the Phantasm synchronization design
 *        (docs/specs/phantasm_frame_sync_spec.md): one local flywheel timebase per
 *        board, disciplined over a single sync-symbol wire.
 *
 * Architecture (spec §3): every board derives its column position from a
 * free-running hardware cycle counter (`x = f(now - epoch)`), so masked-IRQ
 * windows skip columns but do not accumulate phase drift. The master emits count-coded symbol bursts on the one wire — 2/revolution
 * boundary marks, an epoch train (primary + epoch_repeats copies) once per
 * effect, and a quarter-revolution (x ≈ W/4) data beacon on beacon revolutions —
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

#include "pov_sync_protocol.h"
#include "pov_sync_flywheel.h"
#include "pov_sync_content.h"
#include "pov_sync_emitter.h"

#include <atomic>
#include <cstdint>

#include "core/platform/platform.h"

namespace hs_test {
namespace pov_sync_tests {
struct SyncBoardTestAccess;
} // namespace pov_sync_tests
} // namespace hs_test

namespace pov {
namespace sync {

inline void saturating_increment(uint32_t &counter) {
  if (counter != UINT32_MAX)
    ++counter;
}

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
  /**
   * @brief Call advance_display() on the live effect.
   * @details One bool for the whole tick: a wake folding N boundaries runs the
   * flip gate N times — counters and the epoch/join schedule stay exact — but
   * advances the display once, so N windows consume one queued frame and an
   * even N leaves the frame clipped for the opposite arm half (spec §5.1).
   */
  bool flip = false;
  /**
   * @brief The last flip this tick crossed ZERO.
   * @details A tick that folds several boundaries reports only the final one,
   * which is the boundary that opened the display window now in effect.
   */
  bool zero_crossing = false;
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
  HS_COLD_MEMBER explicit SyncBoard(const Config &cfg)
      : protocol_config(cfg), fly(cfg) {
    cache_config_bounds();
  }

  /**
   * @brief Reinitializes the board for a new protocol configuration.
   * @param cfg Protocol configuration.
   * @details Restores every member to its post-construction value; seed() must
   * still follow to establish the timebase.
   * @warning Call only before attaching the board's interrupts: the state is
   * ISR-owned under the spec §8 single-writer model.
   */
  HS_COLD_MEMBER void configure(const Config &cfg) {
    protocol_config = cfg;
    cache_config_bounds();
    fly = Flywheel(cfg);
    is_master_board = false;
    reset_runtime_state();
  }

  /**
   * @brief Boot/reboot seeding (spec §8.5).
   * @param now Current timestamp, in cycles.
   * @param is_master True for the reference board.
   * @details Master is born locked with identity (effect 0, rev 0) — it is the
   * reference. Downstream boards start in ACQUIRE, dark, and join via boundary
   * symbols + beacon.
   * @warning Call only with the board's interrupts detached or disabled; seed
   * replaces ISR-owned state.
   */
  void seed(uint32_t now, bool is_master) {
    is_master_board = is_master;
    fly.seed(now);
    // A reboot must not inherit mailbox or emitter state (spec §8.5).
    reset_runtime_state();
    // Anchor the ACQUIRE quiet-before guard on the seed instant, so a board
    // powering up mid-train cannot snap to an interior digit.
    have_prev_burst = true;
    prev_burst_end = now;
    if (is_master) {
      fly.force_lock();
      content_tracker.identity_known = true;
      content_tracker.effect_index = 0;
      publish_build(0);
    }
  }

  /**
   * @brief Sync-wire edge ISR entry point (publisher; downstream boards only).
   * @param now Edge timestamp, in cycles.
   */
  void on_sync_edge(uint32_t now) {
    edge_mailbox.on_edge(now, protocol_config.glitch_filter_cycles);
  }

  /** @brief Claims a completed burst and ages prior edges under one IRQ mask. */
  bool claim_sync_burst(uint32_t now, BurstSnapshot *out) {
    const uint32_t primask = hs::save_disable_interrupts();
    const bool claimed = edge_mailbox.try_claim(now, cached_gap_timeout_cycles,
                                                cached_max_burst_cycles, out);
    edge_mailbox.age_prior(now, cached_glitch_filter_cycles);
    hs::restore_interrupts(primask);
    return claimed;
  }

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

    // Isolated wire noise can consume at most one rejection per half-revolution.
    if (suspect_pending &&
        (now - suspect_last_cycles) >
            protocol_config.interdigit_timeout_cycles() &&
        static_cast<int32_t>(now - suspect_last_cycles) > 0) {
      suspect_pending = false;
      // The gate exists only in LOCKED; note_rejection() no-ops in ACQUIRE.
      if (fly.lock() == LockState::LOCKED) {
        saturating_increment(telemetry_counters.symbols_rejected_gate);
        if (!have_suspect_rejection ||
            suspect_rejection_half != halves_since_snap) {
          have_suspect_rejection = true;
          suspect_rejection_half = halves_since_snap;
          if (fly.note_rejection())
            saturating_increment(telemetry_counters.lock_transitions);
        }
      }
    }

    // Age out a stale previous-burst timestamp after the ACQUIRE window plus a
    // burst gap, before a cycle-counter wrap collapses the quiet-gap difference.
    // The signed re-check rejects a wrapped modular difference.
    if (have_prev_burst && (now - prev_burst_end) > prev_burst_stale_cycles &&
        static_cast<int32_t>(now - prev_burst_end) > 0) {
      have_prev_burst = false;
      // The same silence ends any partial beacon frame, which could otherwise
      // outlive a counter wrap. A lone digit is an isolated boundary symbol, so
      // a dropped frame needs two.
      if (beacon_parser.digit_count() >= 2)
        saturating_increment(telemetry_counters.beacons_rejected);
      beacon_parser.reset();
    }

    // Re-anchor a stalled master; downstream boards recover through ACQUIRE
    // snaps to wire symbols.
    if (is_master_board && fly.fold_stalled(now)) {
      saturating_increment(telemetry_counters.master_stalls);
      fly.seed(now);
      fly.force_lock();
      // The re-seed stamps ZERO regardless of the pre-stall identity, so the
      // dedup state no longer describes the flywheel. Left standing, a
      // last_flipped of HALF swallows the first crossing after recovery.
      gate = FlipGate{};
    }

    // Fold every locally-crossed boundary (usually 0 or 1; several after a
    // long masked coast).
    for (;;) {
      const Crossing c = fly.fold(now);
      if (!c.crossed)
        break;
      // The master is force-locked and snaps to no wire bursts, so coast is
      // undefined for it.
      if (!is_master_board) {
        saturating_increment(halves_since_snap);
        if (halves_since_snap > telemetry_counters.max_coast_halves)
          telemetry_counters.max_coast_halves = halves_since_snap;
      }
      apply_flip(c.boundary, a);
      if (is_master_board)
        master_on_crossing(c, now);
    }

    int32_t beacon_position = -1;
    if (is_master_board)
      maybe_schedule_beacon(now, beacon_position);

    bool aborted = false;
    if (is_master_board && emitter.tick(now, protocol_config, &aborted))
      a.pulse = true;
    if (aborted)
      saturating_increment(telemetry_counters.emit_aborted);

    // Render decision: dark when phase/identity is missing or during the epoch
    // construction window (spec §6.1).
    a.dark = fly.lock() != LockState::LOCKED ||
             !content_tracker.identity_known ||
             content_tracker.constructing(protocol_config);
    if (!a.dark) {
      const int32_t x =
          beacon_position >= 0 ? beacon_position : fly.position(now);
      if (x != last_rendered_x) { // idempotent wake-up contract (spec §4.1)
        a.render_column = x;
        last_rendered_x = x;
      }
    } else {
      last_rendered_x = -1;
    }
    return a;
  }

  // ── Cross-context interface ───────────────────────────────────────────────

  /**
   * @brief Current build request.
   * @return (generation << 8) | effect_index.
   * @details The foreground constructs the named effect whenever the generation
   * changes, then publishes it to the driver's pending slot.
   */
  uint32_t build_word() const {
    return build_request_word.load(std::memory_order_relaxed);
  }
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
   * @brief Telemetry counters, copied under an IRQ-off window.
   * @return Snapshot of the telemetry block.
   * @details Saves and restores the interrupt mask, so it is safe to call from
   * an ISR or inside an IRQ-off region.
   */
  Telemetry telemetry_snapshot() const {
    const uint32_t primask = hs::save_disable_interrupts();
    const Telemetry snapshot = telemetry_counters;
    hs::restore_interrupts(primask);
    return snapshot;
  }
  /**
   * @brief Output envelope derived from the synchronized effect revolution.
   * @details Flywheel-ISR context only: reads `tick`-owned state unguarded.
   */
  __attribute__((always_inline)) float effect_envelope(int32_t column,
                                                       int32_t width) const {
    return content_tracker.output_envelope(protocol_config, column, width);
  }

private:
  friend struct ::hs_test::pov_sync_tests::SyncBoardTestAccess;

  /**
   * @brief Mailbox access (test-only).
   * @pre Interrupts must be disabled while accessing the returned mailbox.
   * @return Reference to the edge mailbox.
   */
  EdgeMailbox &mailbox() { return edge_mailbox; }
  /**
   * @brief Burst-terminating gap.
   * @return Quiet time that terminates a burst, in cycles.
   */
  uint32_t gap_timeout_cycles() const { return cached_gap_timeout_cycles; }
  /**
   * @brief Burst-duration ceiling.
   * @return Duration past which an unterminated burst is claimed, in cycles.
   */
  uint32_t max_burst_cycles() const { return cached_max_burst_cycles; }
  /**
   * @brief Glitch-filter window.
   * @return Minimum accepted edge spacing, in cycles.
   */
  uint32_t glitch_filter_cycles() const { return cached_glitch_filter_cycles; }

  /**
   * @brief Mutable flywheel access (test-only).
   * @return Mutable reference to the flywheel.
   * @warning Valid only before the ISRs are attached: the flywheel is ISR-owned
   * under the spec §8 single-writer model.
   */
  Flywheel &flywheel_mut() { return fly; }
  /**
   * @brief Mutable content-tracker access (test-only).
   * @return Mutable reference to the content tracker.
   * @warning Valid only before the ISRs are attached: the content tracker is
   * ISR-owned under the spec §8 single-writer model.
   */
  ContentTracker &content_mut() { return content_tracker; }
  /**
   * @brief Content-layer state (test-only).
   * @return Const reference to the content tracker.
   */
  const ContentTracker &content() const { return content_tracker; }
  /**
   * @brief Current lock state (test-only).
   * @return ACQUIRE or LOCKED.
   */
  LockState lock() const { return fly.lock(); }
  /**
   * @brief The flywheel timebase (test-only).
   * @return Const reference to the flywheel.
   */
  const Flywheel &flywheel() const { return fly; }
  /**
   * @brief Protocol configuration (test-only).
   * @return Const reference to the config.
   */
  const Config &config() const { return protocol_config; }

  /**
   * @brief Recomputes the tick-path bounds derived from protocol_config.
   */
  HS_COLD_MEMBER void cache_config_bounds() {
    cached_gap_timeout_cycles = protocol_config.gap_timeout_cycles();
    cached_max_burst_cycles = protocol_config.max_burst_cycles();
    cached_glitch_filter_cycles = protocol_config.glitch_filter_cycles;
    prev_burst_stale_cycles = protocol_config.col_cycles(
        protocol_config.acquire_quiet_cols + protocol_config.gap_timeout_cols);
  }

  /**
   * @brief Restores every member to its post-construction value except
   * protocol_config, cached bounds, fly and is_master_board.
   */
  HS_COLD_MEMBER void reset_runtime_state() {
    gate = FlipGate{};
    content_tracker = ContentTracker{};
    beacon_parser = BeaconParser{};
    emitter = SymbolEmitter{};
    edge_mailbox = EdgeMailbox{};
    telemetry_counters = Telemetry{};
    last_rendered_x = -1;
    halves_since_snap = 0;
    have_prev_burst = false;
    prev_burst_end = 0;
    suspect_pending = false;
    suspect_last_cycles = 0;
    have_suspect_rejection = false;
    suspect_rejection_half = 0;
    epoch_emits_left = 0;
    beacon_done_this_rev = false;
    beacon_busy_counted_this_rev = false;
    beacon_index_candidate = -1;
    build_gen = 0;
    build_request_word.store(0, std::memory_order_relaxed);
  }

  // ── Flip + content events (both paths funnel here) ───────────────────────

  /**
   * @brief Funnels both flip paths (local crossing and symbol backstop) through
   * the dedup gate and runs the per-boundary content events.
   * @param b Boundary being crossed.
   * @param a In/out: tick actions updated with flip/commit/join effects.
   */
  void apply_flip(Boundary b, TickActions &a) {
    if (!gate.try_flip(b))
      return;
    a.flip = true;
    a.zero_crossing = b == Boundary::ZERO;
    saturating_increment(telemetry_counters.flips);
    if (!a.zero_crossing)
      return;
    beacon_done_this_rev = false;
    beacon_busy_counted_this_rev = false;
    if (content_tracker.identity_known) {
      if (content_tracker.on_zero_crossing(protocol_config)) {
        a.commit = true; // B+R+K reached; driver swaps in the pending effect
        // The displayed index just changed; any half-confirmed beacon candidate
        // predates it and can no longer confirm (spec §6.3.4).
        beacon_index_candidate = -1;
      } else if (content_tracker.construction_opens(protocol_config))
        // Last K revolutions: construct the next effect now.
        publish_build((content_tracker.effect_index + 1) %
                      protocol_config.effect_count);
      else if (!content_tracker.commit_pending &&
               (content_tracker.rev_in_effect %
                protocol_config.join_grid_revs) == 0)
        // Late-join point for a board with no live effect.
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
    const bool had_prev = have_prev_burst;
    const uint32_t prev_end = prev_burst_end;
    prev_burst_end = s.last_cycles;
    have_prev_burst = true;

    // A follow-up burst inside the interdigit window proves the pending suspect
    // was the head of a beacon data train.
    if (suspect_pending && (s.first_cycles - suspect_last_cycles) <=
                               protocol_config.interdigit_timeout_cycles())
      suspect_pending = false;

    if (fly.lock() == LockState::LOCKED) {
      // Demarcation (spec §6.4): a burst whose first edge lands far from a
      // predicted boundary is beacon data, not a boundary symbol.
      const int32_t pos = fly.position(s.first_cycles);
      const int32_t to_zero = circ_dist(pos, 0, protocol_config.W);
      const int32_t to_half =
          circ_dist(pos, protocol_config.W / 2, protocol_config.W);
      if (to_zero > protocol_config.gate_cols &&
          to_half > protocol_config.gate_cols) {
        handle_beacon_burst(s);
        // A lone far burst is indistinguishable from a beacon's first digit
        // until a train follows: hold it as a suspect, and tick()'s timeout
        // converts it to a gate rejection if the wire stays silent.
        const bool isolated =
            !had_prev || (s.first_cycles - prev_end) >=
                             protocol_config.acquire_quiet_cycles();
        if (isolated && classify_count(s.count) != Symbol::INVALID) {
          suspect_pending = true;
          suspect_last_cycles = s.last_cycles;
        }
        return;
      }
    } else {
      // ACQUIRE quiet-before guard: a burst following close on another is a
      // beacon digit train, not an isolated boundary symbol — don't hard-snap.
      if (had_prev && (s.first_cycles - prev_end) <
                          protocol_config.acquire_quiet_cycles()) {
        handle_beacon_burst(s);
        return;
      }
      // Before lock an isolated burst may be a train's first digit, so it
      // reaches both paths, as the head of a fresh frame.
      beacon_parser.reset();
      handle_beacon_burst(s);
    }

    const Symbol sym = classify_count(s.count);
    if (sym == Symbol::INVALID) {
      saturating_increment(telemetry_counters.symbols_discarded_invalid);
      return;
    }
    const Boundary b = symbol_boundary(sym);
    const bool was_locked = fly.lock() == LockState::LOCKED;
    int32_t err = 0;
    const Flywheel::SnapOutcome r = fly.snap(b, s.first_cycles, &err);
    if (r != Flywheel::SnapOutcome::ACCEPTED) {
      saturating_increment(telemetry_counters.symbols_rejected_gate);
      if (r == Flywheel::SnapOutcome::REJECTED_FELL_BACK)
        saturating_increment(telemetry_counters.lock_transitions);
      return;
    }
    saturating_increment(telemetry_counters.symbols_accepted);
    if (!was_locked)
      saturating_increment(telemetry_counters.lock_transitions);
    halves_since_snap = 0;
    have_suspect_rejection = false;
    // MUST precede on_epoch_symbol: a ZERO_EPOCH folds rev_in_effect here so the
    // j-inference below reads the post-fold rev (§6.3.1). Deduped against the
    // later fold-loop apply_flip.
    apply_flip(b, a);
    if (sym == Symbol::ZERO_EPOCH && content_tracker.identity_known) {
      if (content_tracker.on_epoch_symbol(protocol_config)) {
        // A board that heard only the last repeat opens the window at the accept.
        if (content_tracker.construction_opens(protocol_config))
          publish_build((content_tracker.effect_index + 1) %
                        protocol_config.effect_count);
      } else {
        saturating_increment(telemetry_counters.epochs_refractory_ignored);
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
    const bool ok = beacon_parser.feed(s, protocol_config, &f, &rejected);
    if (rejected)
      saturating_increment(telemetry_counters.beacons_rejected);
    if (!ok)
      return;
    // An index past the roster is corruption the checksum missed (§6.4).
    if (f.effect_index >= protocol_config.effect_count) {
      saturating_increment(telemetry_counters.beacons_rejected);
      return;
    }
    saturating_increment(telemetry_counters.beacons_ok);
    const int32_t idx = f.effect_index;
    if (!content_tracker.identity_known) {
      // Join (spec §6.4): adopt (effect, rev). Never assume index 0.
      content_tracker.identity_known = true;
      content_tracker.effect_index = idx;
      content_tracker.rev_in_effect = f.rev_count;
      publish_build(idx);
    } else if (content_tracker.commit_pending) {
      // No publish_build mid-window: pending_gen must stay stable from
      // construction-open to commit.
    } else if (idx != content_tracker.effect_index) {
      // A shifted frame can pass the checksum, so a live board needs two
      // consecutive beacons naming the same index (spec §6.3.4).
      if (idx != beacon_index_candidate) {
        beacon_index_candidate = idx;
        return;
      }
      beacon_index_candidate = -1;
      // Confirmed missed epoch (spec §6.3.2).
      content_tracker.effect_index = idx;
      content_tracker.rev_in_effect = f.rev_count;
      saturating_increment(telemetry_counters.beacon_index_corrections);
      publish_build(idx);
    } else {
      beacon_index_candidate = -1;
      if (f.rev_count != (content_tracker.rev_in_effect & 63u)) {
        // Resync a slipped rev counter via the signed mod-64 difference, which
        // recovers any slip under 32 revolutions.
        saturating_increment(telemetry_counters.beacon_rev_mismatches);
        const int32_t d =
            beacon_rev_resync_delta(f.rev_count, content_tracker.rev_in_effect);
        const int64_t fixed =
            static_cast<int64_t>(content_tracker.rev_in_effect) + d;
        content_tracker.rev_in_effect =
            fixed >= 0 ? static_cast<uint32_t>(fixed) : f.rev_count;
      }
    }
  }

  // ── Conduct + emit path (master) ──────────────────────────────────────────

  /**
   * @brief Master conductor: emit the boundary symbol for a crossing and run
   * the epoch-train schedule.
   * @param c The crossing just folded.
   * @param now Current timestamp, in cycles.
   */
  void master_on_crossing(const Crossing &c, uint32_t now) {
    Symbol sym = Symbol::HALF;
    if (c.boundary == Boundary::ZERO) {
      // Conductor (spec §6.1): when the effect's revolutions elapse, start an
      // EPOCH train — primary copy plus R repeats on following ZERO boundaries.
      if (epoch_emits_left == 0 && !content_tracker.commit_pending &&
          content_tracker.rev_in_effect >=
              protocol_config.revolutions_for_effect(
                  content_tracker.effect_index)) {
        epoch_emits_left = 1 + protocol_config.epoch_repeats;
        if (content_tracker.on_epoch_symbol(protocol_config) &&
            content_tracker.construction_opens(protocol_config))
          publish_build((content_tracker.effect_index + 1) %
                        protocol_config.effect_count);
      }
      if (epoch_emits_left > 0) {
        sym = Symbol::ZERO_EPOCH;
        // Spent by the boundary, not by reaching the wire: receivers invert j
        // from their own revolution count, so the train must occupy exactly
        // boundaries B..B+R (spec §6.3.1).
        --epoch_emits_left;
      } else {
        sym = Symbol::ZERO;
      }
    }
    // Discard a stale emission before scheduling this boundary; count its kind.
    switch (emitter.drop_pending_emission()) {
    case SymbolEmitter::DroppedBurst::BEACON:
      saturating_increment(telemetry_counters.beacons_overrun_dropped);
      break;
    case SymbolEmitter::DroppedBurst::BOUNDARY:
      saturating_increment(telemetry_counters.boundary_bursts_dropped);
      break;
    case SymbolEmitter::DroppedBurst::NONE:
      break;
    }
    if (!emitter.schedule_boundary(sym, c.at_cycles, now, protocol_config))
      saturating_increment(telemetry_counters.emit_censored);
  }

  /**
   * @brief Master: queue a beacon frame once per revolution at the beacon point
   * when one is due.
   * @param now Current timestamp, in cycles.
   * @param position Receives the column if the due-beacon check computes it.
   */
  void maybe_schedule_beacon(uint32_t now, int32_t &position) {
    // Beacon point: x ≈ W/4, mid-way through the ZERO→HALF half-rev where the
    // wire is otherwise quiet (spec §6.4).
    if (beacon_done_this_rev || fly.current_boundary() != Boundary::ZERO)
      return;
    // Silent during the commit window: a beacon here broadcasts the outgoing
    // index, and a board joining off it would adopt stale identity.
    if (content_tracker.commit_pending)
      return;
    const uint32_t rev = content_tracker.rev_in_effect;
    const bool due = (rev % protocol_config.beacon_period_revs) == 1u ||
                     (rev >= 1u && rev <= static_cast<uint32_t>(
                                              protocol_config.epoch_repeats));
    if (!due)
      return;
    // A coalesced coast past the boundary skips this beacon.
    position = fly.position(now);
    if (position < protocol_config.W / 4)
      return;
    uint8_t digits[5];
    encode_beacon_digits(content_tracker.effect_index, rev, digits);
    // The frame, receiver tail quiet and emission lateness must fit before HALF.
    // Measure from this tick, including its sub-column offset.
    int32_t digit_sum = 0;
    for (int i = 0; i < 5; ++i)
      digit_sum += digits[i];
    const uint32_t frame_cycles =
        protocol_config.col_cycles(
            protocol_config.beacon_frame_cols(digit_sum)) +
        protocol_config.late_censor_cycles();
    if (static_cast<int32_t>(frame_cycles) > fly.cycles_to_next_boundary(now)) {
      // The margin only shrinks this revolution.
      beacon_done_this_rev = true;
      saturating_increment(telemetry_counters.beacons_late_dropped);
      return;
    }
    if (emitter.schedule_beacon(digits, now, protocol_config))
      beacon_done_this_rev = true;
    else if (!beacon_busy_counted_this_rev) {
      saturating_increment(telemetry_counters.beacons_busy_dropped);
      beacon_busy_counted_this_rev = true;
    }
  }

  /**
   * @brief Publishes a build request for @p index, bumping the generation.
   * @param index Effect index the foreground should construct.
   */
  void publish_build(int32_t index) {
    ++build_gen;
    build_request_word.store((build_gen << 8) |
                                 static_cast<uint32_t>(index & 0xFF),
                             std::memory_order_relaxed);
  }

  // ── State ─────────────────────────────────────────────────────────────────

  uint32_t cached_gap_timeout_cycles = 0;
  uint32_t cached_max_burst_cycles = 0;
  uint32_t cached_glitch_filter_cycles = 0;
  Config protocol_config;
  Flywheel fly;
  FlipGate gate;
  ContentTracker content_tracker;
  BeaconParser beacon_parser;
  SymbolEmitter emitter;
  EdgeMailbox edge_mailbox;
  Telemetry telemetry_counters;

  bool is_master_board = false;
  int32_t last_rendered_x = -1;
  uint32_t halves_since_snap = 0;
  bool have_prev_burst = false;
  uint32_t prev_burst_end = 0;
  uint32_t prev_burst_stale_cycles =
      0; /**< Quiet past which prev_burst_end ages out (cached). */
  bool suspect_pending = false; /**< Lone far burst awaiting train/timeout. */
  uint32_t suspect_last_cycles = 0;
  bool have_suspect_rejection = false;
  uint32_t suspect_rejection_half = 0;
  uint32_t epoch_emits_left =
      0; /**< ZERO boundaries left in the EPOCH train. */
  bool beacon_done_this_rev = false;
  bool beacon_busy_counted_this_rev = false;
  int32_t beacon_index_candidate =
      -1; /**< Index a lone beacon named, awaiting confirmation (§6.3.4). */
  uint32_t build_gen = 0;
  static_assert(std::atomic<uint32_t>::is_always_lock_free);
  std::atomic<uint32_t> build_request_word{
      0}; /**< (gen << 8) | index; foreground-read. */
};

} // namespace sync
} // namespace pov
