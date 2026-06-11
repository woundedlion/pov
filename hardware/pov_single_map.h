/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file pov_single_map.h
 * @brief Pure strip index arithmetic for the single-board POV driver.
 *
 * Split out of pov_single.h (which is Arduino-only) so the load-bearing index
 * math — which physical LED samples which canvas column and row — is
 * unit-testable on the host without a Teensy, exactly as pov_segment_map.h was
 * split out of pov_segmented.h. The single-board show_col() ISR derives its
 * mapping from these functions, so the host tests cover the real arithmetic. An
 * off-by-one here silently mis-paints the sphere.
 *
 * Layout (one Teensy owns the whole S-LED strip; ROWS = S/2 = canvas height):
 *   The strip spans both sides of the ring. Its first half [0, S/2) is the top
 *   arm, physically reversed (LED 0 at the junction end), and samples canvas
 *   column x; its second half [S/2, S) is the bottom arm, straight, and samples
 *   the opposite half of the image (column (x + W/2) % W). At rotation column x
 *   the two halves together paint exactly the two canvas columns x and
 *   (x + W/2) % W, one LED per (column, row).
 */
#pragma once

namespace pov {

/// Physical LED index of the top-half strip pixel sampling canvas row @p y.
/// Top half spans LEDs [0, S/2) and is wired reversed: row y -> LED S/2-1-y.
/// (i.e. y=0 is the last top-half LED, y=S/2-1 the first.)
constexpr int strip_top_led(int y, int S) { return S / 2 - 1 - y; }

/// Physical LED index of the bottom-half strip pixel sampling canvas row @p y.
/// Bottom half spans LEDs [S/2, S) straight: row y -> LED S/2 + y.
constexpr int strip_bottom_led(int y, int S) { return S / 2 + y; }

/// Canvas column the bottom half samples at rotation column @p x (canvas width
/// @p w): the opposite half of the image. The top half samples @p x directly.
constexpr int strip_opposite_col(int x, int w) { return (x + w / 2) % w; }

} // namespace pov
