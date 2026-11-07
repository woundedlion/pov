/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief HD107S shared calibration reset helper. */
#pragma once
#include "hardware/hd107s_frame.h"

namespace hs_test {
/**
 * @brief Restores HD107SFrame<S>'s shared static correction state to unity
 *        (255 = exact ×1.0).
 * @tparam S Pixel count of the frame instantiation to reset; the correction
 *           state is a static of HD107SFrame<S>, so each S carries its own.
 * @details Resets the per-channel correction, temperature, and brightness to an
 * identity color pipeline.
 */
template <int S> inline void reset_correction() {
  HD107SFrame<S>::set_correction(255, 255, 255);
  HD107SFrame<S>::set_temperature(255, 255, 255);
  HD107SFrame<S>::set_brightness(255);
}

} // namespace hs_test
