/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file profile_clock.h
 * @brief Interrupt-matched wall and cycle clock samples. */

#include <cstdint>

namespace hs {

/** @brief Independent clock values captured under one interrupt mask. */
struct ProfileClockSample {
  uint32_t wall_us; ///< Wall-clock microseconds.
  uint32_t cycles;  ///< CPU cycle count.
};

/** @brief Captures both clocks and restores the preceding interrupt mask.
 * @param wall_clock Reads wall-clock microseconds.
 * @param cycle_clock Reads CPU cycles.
 * @param disable Masks interrupts and returns the preceding mask.
 * @param restore Restores that mask.
 * @return Matched wall and cycle samples. */
template <typename WallClock, typename CycleClock, typename Disable,
          typename Restore>
inline ProfileClockSample
sample_profile_clocks(WallClock wall_clock, CycleClock cycle_clock,
                      Disable disable, Restore restore) {
  const auto mask = disable();
  const ProfileClockSample sample{static_cast<uint32_t>(wall_clock()),
                                  static_cast<uint32_t>(cycle_clock())};
  restore(mask);
  return sample;
}

} // namespace hs
