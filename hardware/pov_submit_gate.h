/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file pov_submit_gate.h
 * @brief Pure, host-testable LED-submit and sync-pulse decisions for the
 *        POVSegmented ISR.
 *
 * Split out of pov_segmented.h (Arduino-only) so the one place the driver
 * reacts to the DMA transport's accept/drop verdict is unit-testable on the
 * host, exactly as pov_handoff.h is for the effect handoff. The transport is
 * not modelled here: the verdict is injected as a bool.
 *
 * DMALEDController::submitFrame() drops a frame when the previous transfer is
 * still in flight, so both submit paths — the fail-dark black frame and the
 * image column — clear their pending state only on an accepted submit. A drop
 * then costs one flywheel wake (~54 µs at 8× oversampling) rather than a whole
 * column: the flywheel's idempotent wake-up contract (spec §4.1) reports no new
 * column until the arm advances, so the ~7 remaining wakes of a dropped column
 * would otherwise all be no-ops and the column would paint dark.
 *
 * Re-submission needs no repack — an overrun returns before the controller
 * swaps buffers, so the dropped frame's pixels are still packed in backFrame().
 */
#pragma once

#include <cstdint>

namespace pov {

/** @brief What one flywheel wake hands to the LED transport. */
enum class SubmitAction : uint8_t {
  NONE,     /**< Nothing to submit this wake. */
  BLACK,    /**< Pack and submit an all-black frame (fail-dark). */
  COLUMN,   /**< Pack and submit the flywheel's canvas column. */
  RESUBMIT, /**< Re-submit the frame an earlier wake had dropped. */
};

/**
 * @brief Per-wake LED-submit decision and its two retry latches.
 * @details ISR-owned: every access is in the flywheel-ISR context, so the
 *          latches are plain bools.
 */
class SubmitGate {
public:
  /**
   * @brief Chooses this wake's submit.
   * @param dark True when the wake must fail dark: no lock/identity, the epoch
   *        construction window, or no live effect.
   * @param render_column Canvas column the flywheel wants drawn, or < 0 when
   *        this wake repeats an already-rendered column.
   * @return The action to perform, then report back through settle().
   */
  SubmitAction choose(bool dark, int32_t render_column) {
    if (dark) {
      // A dropped image column is stale once the wake goes dark, and the black
      // frame overwrites the back buffer it was packed into.
      resubmit_needed = false;
      return black_frame_accepted ? SubmitAction::NONE : SubmitAction::BLACK;
    }
    black_frame_accepted = false;
    if (render_column >= 0)
      return SubmitAction::COLUMN;
    return resubmit_needed ? SubmitAction::RESUBMIT : SubmitAction::NONE;
  }

  /**
   * @brief Records the transport's verdict for the action choose() returned.
   * @param action The action just performed.
   * @param accepted submitFrame()'s return value; ignored for NONE.
   * @return True once this wake's LED work is complete — the caller's gate for
   *         dropping a same-wake sync pulse instead of deferring it.
   */
  bool settle(SubmitAction action, bool accepted) {
    switch (action) {
    case SubmitAction::BLACK:
      black_frame_accepted = accepted;
      return accepted;
    case SubmitAction::COLUMN:
    case SubmitAction::RESUBMIT:
      resubmit_needed = !accepted;
      return accepted;
    case SubmitAction::NONE:
      break;
    }
    return false;
  }

  /** @brief Whether the fail-dark black frame has been accepted. */
  bool dark_latched() const { return black_frame_accepted; }
  /** @brief Whether a dropped frame is still awaiting re-submission. */
  bool resubmit_pending() const { return resubmit_needed; }

private:
  bool black_frame_accepted =
      false; /**< Black frame accepted; no repeat needed. */
  bool resubmit_needed =
      false; /**< Back frame dropped on overrun; retry next wake. */
};

/**
 * @brief Width decision for the master's sync pulse on the shared sync wire.
 * @details ISR-owned, like SubmitGate. A wake that renders a frame has a body
 *          long enough to width the pulse to spec §5.2's "tens of µs", so the
 *          pin drops before the wake returns. A wake that renders nothing holds
 *          the pin HIGH across the ISR boundary and drops it at the head of the
 *          next wake instead.
 */
class SyncPulseGate {
public:
  /**
   * @brief Claims a pulse the previous wake left HIGH.
   * @return True when the caller must drive the sync pin LOW now; the latch is
   *         cleared by the claim.
   */
  bool take_deferred_low() {
    if (!low_pending)
      return false;
    low_pending = false;
    return true;
  }

  /**
   * @brief Ends this wake's pulse, or defers it to the next wake.
   * @param pulse Whether this wake drove the sync pin HIGH.
   * @param did_render Whether the wake's LED work completed, widening the body
   *        enough to carry the pulse.
   * @return True when the caller must drive the sync pin LOW before returning.
   */
  bool settle(bool pulse, bool did_render) {
    if (!pulse)
      return false;
    if (did_render)
      return true;
    low_pending = true;
    return false;
  }

  /** @brief Whether a pulse is being held HIGH across the ISR boundary. */
  bool low_deferred() const { return low_pending; }

private:
  bool low_pending = false; /**< Pulse held HIGH; drop it at the next wake. */
};

} // namespace pov
