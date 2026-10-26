/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file profile_isr_window.h
 * @brief Atomic snapshots of profiling ISR accumulators. */
#pragma once

#include "core/platform/profiling.h"

namespace hs {

/** @brief Immutable ISR accumulators and their wall-clock interval. */
struct ProfileIsrSnapshot {
  IsrCycleStats wake,     ///< Flywheel wake timings.
      pack,               ///< Column packing timings.
      submit;             ///< DMA submission timings.
  uint32_t window_us = 0; ///< Elapsed wall-clock microseconds.
};

/** @brief Captures consecutive ISR windows under an interrupt mask. */
class ProfileIsrWindow {
public:
  /** @brief Copies and clears ISR counters with interrupts masked.
   * @tparam Clock Microsecond clock callable.
   * @tparam Disable Callable returning the saved interrupt mask.
   * @tparam Restore Callable restoring that mask.
   * @param wake Flywheel wake counters.
   * @param pack Column packing counters.
   * @param submit DMA submission counters.
   * @param clock Current wall-clock microseconds.
   * @param disable Masks interrupts.
   * @param restore Restores the saved interrupt mask.
   * @return Counters and elapsed time since the preceding capture. */
  template <typename Clock, typename Disable, typename Restore>
  ProfileIsrSnapshot capture(IsrCycleStats &wake, IsrCycleStats &pack,
                             IsrCycleStats &submit, Clock clock,
                             Disable disable, Restore restore) {
    const auto mask = disable();
    const uint32_t now = clock();
    const ProfileIsrSnapshot snapshot{wake, pack, submit, now - previous_us};
    wake.reset();
    pack.reset();
    submit.reset();
    previous_us = now;
    restore(mask);
    return snapshot;
  }

private:
  uint32_t previous_us = 0;
};

} // namespace hs
