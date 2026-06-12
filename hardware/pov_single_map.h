/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file pov_single_map.h
 * @brief Pure strip index arithmetic for the single-board POV driver.
 *
 * Kept free of Arduino dependencies so the load-bearing index math — which
 * physical LED samples which canvas column and row — is unit-testable on the
 * host without a Teensy. The single-board show_col() ISR derives its mapping
 * from these functions, so the host tests cover the real arithmetic. An
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

/**
 * @brief Physical LED index of the top-half strip pixel sampling canvas row y.
 * @param y Canvas row index in [0, S/2).
 * @param S Total LED count on the strip; top half spans LEDs [0, S/2).
 * @return Strip LED index in [0, S/2) sampling row y.
 * @details Top half is wired reversed: row y -> LED S/2-1-y, so y=0 is the last
 *          top-half LED and y=S/2-1 the first.
 */
constexpr int strip_top_led(int y, int S) { return S / 2 - 1 - y; }

/**
 * @brief Physical LED index of the bottom-half strip pixel sampling canvas row y.
 * @param y Canvas row index in [0, S/2).
 * @param S Total LED count on the strip; bottom half spans LEDs [S/2, S).
 * @return Strip LED index in [S/2, S) sampling row y.
 * @details Bottom half is wired straight: row y -> LED S/2 + y.
 */
constexpr int strip_bottom_led(int y, int S) { return S / 2 + y; }

/**
 * @brief Canvas column the bottom half samples at rotation column x.
 * @param x Rotation (top-half) canvas column index in [0, w).
 * @param w Canvas width in columns.
 * @return Opposite-half canvas column (x + w/2) % w in [0, w).
 * @details The bottom half paints the opposite half of the image; the top half
 *          samples column x directly.
 */
constexpr int strip_opposite_col(int x, int w) { return (x + w / 2) % w; }

} // namespace pov
