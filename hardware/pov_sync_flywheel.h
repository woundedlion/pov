/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file pov_sync_flywheel.h
 * @brief Layer 1 of the Phantasm sync design (spec §4): the per-board
 *        position-from-time flywheel and its snap discipline.
 */
#pragma once

#include "pov_sync_protocol.h"

#include <cassert>
#include <cstdint>

#include "core/engine/platform.h" // HS_CHECK used in Flywheel's constructor guard

namespace pov {
namespace sync {

// ── Layer 1: the flywheel timebase (spec §4) ────────────────────────────────

/**
 * @brief Flywheel lock state: ACQUIRE (hunting) or LOCKED (disciplined).
 */
enum class LockState : uint8_t { ACQUIRE, LOCKED };

/**
 * @brief A locally-crossed boundary, reported by Flywheel::fold().
 */
struct Crossing {
  bool crossed = false;               /**< True if a boundary was crossed. */
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
      : period(cfg.cycles_per_half_rev), w(cfg.W), gate_cols(cfg.gate_cols),
        reject_fallback(cfg.reject_fallback) {
    check_period(period);
  }

  /**
   * @brief Boot seeding (spec §8.5): epoch = now, boundary ZERO, ACQUIRE.
   * @param now Current timestamp, in cycles.
   */
  void seed(uint32_t now) {
    epoch_cycles = now;
    boundary = Boundary::ZERO;
    lock_state = LockState::ACQUIRE;
    consecutive_rejects = 0;
  }

  /**
   * @brief Forces the flywheel LOCKED.
   * @details Master is the reference by definition — it never snaps and is born
   * locked (spec §2: "Master's flywheel IS the reference").
   */
  void force_lock() { lock_state = LockState::LOCKED; }

  /**
   * @brief Current column at a given time.
   * @param at Timestamp to evaluate, in cycles.
   * @return Column index in [0, W).
   * @details Signed-safe for timestamps up to ±16 half-revolutions around the
   * epoch (±1.0 s at 480 RPM); snap evaluation may look slightly into the past.
   */
  int32_t position(uint32_t at) const {
    const uint32_t elapsed = at - epoch_cycles; // modular
    // Signed-safe window: `at` must lie within MIN_SAFE_HALF_REVS half-revs
    // either side of the epoch, the coast the constructor sized the int32 cast
    // for. Either sign is reachable — demarcation and the snap gate evaluate a
    // past first-edge timestamp, and both run before tick()'s fold loop, so
    // elapsed spans the whole coast since the last fold. Test the unsigned
    // magnitude: a forward elapsed past the window wraps to a small-negative
    // int32 that slips under a signed bound.
    assert(
        (elapsed < static_cast<uint32_t>(MIN_SAFE_HALF_REVS) * period ||
         elapsed >= 0u - static_cast<uint32_t>(MIN_SAFE_HALF_REVS) * period) &&
        "Flywheel::position: timestamp outside the signed-safe coast window");
    const int64_t delta = static_cast<int32_t>(elapsed);
    const int64_t cols = floor_div(delta * (w / 2), period);
    return floor_mod(boundary_column(boundary, w) + cols, w);
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
    const int32_t delta = static_cast<int32_t>(now - epoch_cycles);
    if (delta < 0 || static_cast<uint32_t>(delta) < period)
      return {};
    epoch_cycles += period;
    boundary = opposite(boundary);
    return {true, boundary, epoch_cycles};
  }

  /**
   * @brief Whether fold() has stopped being able to report crossings.
   * @param now Current timestamp, in cycles.
   * @return True once (now - epoch_cycles) has passed 2^31 cycles.
   * @details fold() reads that difference as int32 and returns no crossing on a
   * negative one, so past 2^31 the epoch can never catch up: every later call
   * sees an even larger difference. Only reachable on a coast that long, since
   * the fold loop otherwise leaves the difference below one half-rev.
   */
  bool fold_stalled(uint32_t now) const {
    return static_cast<int32_t>(now - epoch_cycles) < 0;
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
    const int32_t target = boundary_column(b, w);
    const int32_t err = circ_dist(position(edge_cycles), target, w);
    if (error_cols)
      *error_cols = err;
    if (lock_state == LockState::LOCKED && err > gate_cols) {
      return note_rejection() ? SnapOutcome::REJECTED_FELL_BACK
                              : SnapOutcome::REJECTED;
    }
    epoch_cycles = edge_cycles;
    boundary = b;
    lock_state = LockState::LOCKED;
    consecutive_rejects = 0;
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
    if (lock_state != LockState::LOCKED)
      return false;
    if (++consecutive_rejects >= reject_fallback) {
      lock_state = LockState::ACQUIRE;
      consecutive_rejects = 0;
      return true;
    }
    return false;
  }

  /**
   * @brief Current lock state.
   * @return ACQUIRE or LOCKED.
   */
  LockState lock() const { return lock_state; }
  /**
   * @brief Boundary identity at the current epoch.
   * @return The boundary at epoch_cycles.
   */
  Boundary current_boundary() const { return boundary; }
  /**
   * @brief §4.3 frequency trim hook (snap-only ships; tests exercise extremes).
   * @param c New half-rev period, in cycles.
   */
  void set_cycles_per_half_rev(uint32_t c) {
    check_period(c);
    period = c;
  }

private:
  // Min coast (in half-revs) position()'s int32 elapsed cast must survive.
  static constexpr uint32_t MIN_SAFE_HALF_REVS = 16;

  /**
   * @brief Traps a half-rev period position()'s arithmetic cannot carry.
   * @param p Candidate cycles_per_half_rev.
   * @details position() reinterprets (at - epoch_cycles) as int32 and divides by
   * the period, so zero divides by zero and a period above INT32_MAX /
   * MIN_SAFE_HALF_REVS voids the signed-safe coast window.
   */
  static void check_period(uint32_t p) {
    HS_CHECK(p > 0 &&
                 p <= static_cast<uint32_t>(INT32_MAX) / MIN_SAFE_HALF_REVS,
             "Flywheel: cycles_per_half_rev outside the range position()'s "
             "int32 elapsed window holds for MIN_SAFE_HALF_REVS of coast");
  }

  uint32_t period; /**< cycles_per_half_rev (optionally trimmed). */
  int32_t w;
  int32_t gate_cols;
  int32_t reject_fallback;
  uint32_t epoch_cycles = 0;
  Boundary boundary = Boundary::ZERO; /**< Column identity at epoch_cycles. */
  LockState lock_state = LockState::ACQUIRE;
  int32_t consecutive_rejects = 0;
};

} // namespace sync
} // namespace pov
