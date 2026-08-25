/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file pov_sync_emitter.h
 * @brief Master-side symbol generation (spec §5.2, §6.4): pulse scheduling
 *        in cycle time with self-censoring of late bursts.
 */
#pragma once

#include "pov_sync_protocol.h"

#include <cstdint>

namespace pov {
namespace sync {

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
  /** @brief What drop_pending_emission() discarded. */
  enum class DroppedBurst : uint8_t { NONE, BOUNDARY, BEACON };

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
    // Retire a fully drained beacon frame here as well as in tick(), so a live
    // queue_len always means the in-flight pulses belong to a beacon.
    queue_len = queue_pos = 0;
    pulses_left = symbol_pulse_count(symbol);
    next_due = at_cycles;
    pitch = cfg.pulse_pitch_cycles();
    return true;
  }

  /**
   * @brief Queues the five digit bursts of a beacon frame.
   * @param digits The five base-8 beacon digits.
   * @param now Current timestamp (frame start), in cycles.
   * @param cfg Protocol configuration.
   * @return False when another emission is still active.
   * @details Called when the master reaches the beacon point (x ≈ W/4).
   */
  bool schedule_beacon(const uint8_t digits[5], uint32_t now,
                       const Config &cfg) {
    if (pulses_left > 0 || queue_pos < queue_len)
      return false; // defensive: never interleave with an active emission
    uint32_t start = now;
    const uint32_t col = cfg.cycles_per_column();
    for (int i = 0; i < 5; ++i) {
      queue[i] = {start, static_cast<uint32_t>(digits[i]) + 1u};
      // Next burst starts after this one's span (digits[i] pitches) plus a
      // gap the decoder is guaranteed to see as burst-terminating.
      start += digits[i] * cfg.beacon_pitch_cycles() +
               static_cast<uint32_t>(cfg.gap_timeout_cols + 1) * col;
    }
    queue_len = 5;
    queue_pos = 0;
    return true;
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
    if (pulses_left == 0 && !queue_pending() && queue_len != 0)
      queue_len = queue_pos = 0; // beacon frame drained: emitter idle again
    if (pulses_left == 0 && queue_pending()) {
      const PendingBurst &b = queue[queue_pos];
      if (static_cast<int32_t>(now - b.start_cycles) >= 0) {
        pulses_left = b.pulses;
        next_due = b.start_cycles;
        pitch = cfg.beacon_pitch_cycles();
        ++queue_pos;
      }
    }
    if (pulses_left == 0)
      return false;
    if (static_cast<int32_t>(now - next_due) < 0)
      return false; // next pulse not due yet
    if ((now - next_due) > cfg.late_censor_cycles()) {
      // Masked past the lateness budget mid-emission: stop, and drop any
      // queued beacon digits too (a partial frame must fail, not mislead).
      pulses_left = 0;
      queue_len = queue_pos = 0;
      *aborted = true;
      return false;
    }
    --pulses_left;
    next_due += pitch;
    return true;
  }

  /**
   * @brief Whether the emitter has nothing queued or in flight.
   * @return True if no pulses or beacon bursts remain.
   */
  bool idle() const { return pulses_left == 0 && !queue_pending(); }

  /**
   * @brief Drops any in-flight or queued burst so a boundary symbol can schedule.
   * @return What was discarded, for caller telemetry.
   * @details Usually a beacon that overran its pre-HALF window under a masked-ISR
   * coast, but a boundary symbol masked for most of a half-rev can also still
   * have undrained pulses. Either way the burst is stale, and schedule_boundary
   * would clobber it silently; dropping it first is what turns the discard into
   * telemetry (degrade to missed, never to wrong).
   */
  DroppedBurst drop_pending_emission() {
    if (pulses_left == 0 && queue_pos >= queue_len)
      return DroppedBurst::NONE;
    const DroppedBurst kind =
        queue_len != 0 ? DroppedBurst::BEACON : DroppedBurst::BOUNDARY;
    pulses_left = 0;
    queue_len = queue_pos = 0;
    return kind;
  }

private:
  /** @brief Whether a queued beacon burst is still waiting to be emitted. */
  bool queue_pending() const { return queue_pos < queue_len; }

  /**
   * @brief A queued beacon-digit burst awaiting its start time.
   */
  struct PendingBurst {
    uint32_t start_cycles = 0; /**< When the burst should begin, in cycles. */
    uint32_t pulses = 0;       /**< Number of pulses in the burst. */
  };
  PendingBurst queue[5] = {};
  int32_t queue_len = 0;
  int32_t queue_pos = 0;
  uint32_t pulses_left = 0;
  uint32_t next_due = 0;
  uint32_t pitch = 0;
};

} // namespace sync
} // namespace pov
