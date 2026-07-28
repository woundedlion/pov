/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file pov_submit_gate.h
 * @brief Pure, host-testable LED-submit decision for the POVSegmented ISR.
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

} // namespace pov
