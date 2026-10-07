/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/platform/profiling.h"

namespace hs {

/** @brief Immutable ISR accumulators and their wall-clock interval. */
struct ProfileIsrSnapshot {
  IsrCycleStats wake, pack, submit;
  uint32_t window_us = 0;
};

/** @brief Captures consecutive ISR windows under an interrupt mask. */
class ProfileIsrWindow {
public:
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
