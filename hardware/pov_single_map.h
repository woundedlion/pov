/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file pov_single_map.h
 * @brief Pure strip index and column-cadence arithmetic for the single-board
 *        POV driver.
 *
 * Kept free of Arduino dependencies so the load-bearing math — which physical
 * LED samples which canvas column and row, how long a column lasts, and when
 * the display buffer advances — is unit-testable on the host without a Teensy.
 * The single-board show_col() ISR and run() loop derive their behaviour from
 * these functions, so the host tests cover the real arithmetic. An off-by-one
 * here silently mis-paints the sphere or freezes the strip.
 *
 * Layout (one Teensy owns the whole S-LED strip; ROWS = S/2 = canvas height):
 *   The strip spans both sides of the ring. Its first half [0, S/2) is the top
 *   arm, physically reversed (LED 0 at the junction end), and samples canvas
 *   column x; its second half [S/2, S) is the bottom arm, straight, and samples
 *   the opposite half of the image (column (x + W/2) % W). At rotation column x
 *   the two halves together paint exactly the two canvas columns x and
 *   (x + W/2) % W, one LED per (column, row).
 *
 * NOTE: the top-arm wiring convention deliberately differs from the segmented
 * rig (pov_segment_map.h). Here the top arm is reversed with LED 0 at the
 * junction (strip_top_led = S/2-1-y); the segmented rig wires its top segment
 * with LED 0 at the N pole, NOT reversed (y_step = +1). The two are separate
 * physical builds — do not assume one map's top-arm direction carries over.
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

/**
 * @brief Column-sweep timer period, in µs.
 * @param cols_per_min Columns swept per minute: RPM × canvas width.
 * @return Microseconds per column, rounded to nearest.
 * @pre cols_per_min > 0; the driver traps a zero before calling.
 * @details One revolution is 60 s / RPM, split into `width` columns. Rounded
 *          rather than truncated so the sweep does not run systematically fast
 *          and drift the image against the rotation.
 */
constexpr unsigned long column_interval_us(unsigned long cols_per_min) {
  return (60000000UL + cols_per_min / 2) / cols_per_min;
}

/**
 * @brief What the column ISR does at the end of a column.
 */
struct ColumnStep {
  int next_x;   /**< Rotation column the next ISR entry paints. */
  bool advance; /**< Whether the display buffer advances at that column. */
};

/**
 * @brief Advances the rotation column and reports the half-revolution flip.
 * @param x Rotation column just painted, in [0, w).
 * @param w Canvas width in columns; must be even.
 * @return The next column, and whether a display-buffer advance falls on it.
 * @details The two strip halves sample opposite canvas columns, so a whole
 *          frame is painted every half revolution: the buffer advances at the
 *          two half-revolution boundaries, columns 0 and w/2.
 */
constexpr ColumnStep step_column(int x, int w) {
  const int next_x = (x + 1) % w;
  return {next_x, next_x == 0 || next_x == w / 2};
}

} // namespace pov
