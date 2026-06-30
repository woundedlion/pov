/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file pov_segment_map.h
 * @brief Pure segment index arithmetic for the segmented POV driver.
 *
 * Kept free of Arduino dependencies so the load-bearing index math — which
 * segment owns which canvas column and rows, and in which direction its
 * physical strip runs — is unit-testable on the host without a Teensy. The
 * driver's boot-time configure_segment() and per-column ISR derive their
 * mapping from these functions, so the host tests cover the real arithmetic.
 * An off-by-one here silently mis-paints the sphere.
 *
 * Layout (N=4, S=288 -> ROWS=144, PPS=72): two arms, two segments each.
 *   Arm A samples canvas column x; arm B samples column (x + W/2) % W.
 *   Per arm: the top segment runs y = 0 -> PPS-1 (N pole -> junction),
 *            the bottom segment runs y = ROWS-1 -> ROWS-PPS (S pole ->
 *            junction, strip physically reversed).
 *
 * NOTE: the top-arm wiring convention deliberately differs from the single-board
 * rig (pov_single_map.h). Here the top segment wires LED 0 at the N pole, NOT
 * reversed (y_step = +1); the single-board rig reverses its top arm with LED 0
 * at the junction (strip_top_led = S/2-1-y). The two are separate physical
 * builds — do not assume one map's top-arm direction carries over.
 */
#pragma once

namespace pov {

/**
 * @brief Precomputed vertical mapping for one physical segment.
 */
struct SegmentMap {
  bool arm_b;  /**< Samples canvas column (x + W/2) instead of x. */
  int y_base;  /**< Canvas row of this segment's LED 0. */
  int y_step;  /**< +1 (top, toward junction) or -1 (bottom strip, reversed). */
};

/**
 * @brief Derive a segment's vertical mapping from its hardware ID.
 * @param segment_id  Hardware-strapped ID in [0, N).
 * @param S  Total LEDs across both arms (ROWS = S/2).
 * @param N  Segment count (even, power of two, <= 4; SEGS_PER_ARM = N/2).
 * @return SegmentMap with arm_b, y_base (canvas row of LED 0), and y_step (+1/-1).
 * @details Segments [0, N/2) are arm A; [N/2, N) are arm B. Within an arm,
 * arm-segment 0 is the top strip (LED 0 at the N pole, y=0, counting toward the
 * junction) and arm-segment 1 is the bottom strip (LED 0 at the S pole,
 * y=ROWS-1, counting toward the junction — physically reversed).
 */
constexpr SegmentMap segment_map(int segment_id, int S, int N) {
  const int segs_per_arm = N / 2;
  const int rows = S / 2;
  // segs_per_arm is a power of two, so mask for the within-arm index.
  const int arm_seg = segment_id & (segs_per_arm - 1);
  SegmentMap m{};
  m.arm_b = segment_id >= segs_per_arm;
  if (arm_seg == 0) {
    m.y_base = 0; // top: N pole -> junction
    m.y_step = 1;
  } else {
    m.y_base = rows - 1; // bottom: S pole -> junction (reversed)
    m.y_step = -1;
  }
  return m;
}

/**
 * @brief Canvas column this segment samples at a given rotation column.
 * @param arm_b True if this segment reads the opposite half of the image.
 * @param x Rotation column (canvas column index for arm A).
 * @param w Canvas width in columns.
 * @return Canvas column index; arm B returns (x + w/2) % w, arm A returns x.
 */
constexpr int segment_x_col(bool arm_b, int x, int w) {
  return arm_b ? (x + w / 2) % w : x;
}

/**
 * @brief Canvas row of this segment's i-th LED.
 * @param m Segment mapping providing y_base and y_step.
 * @param i LED index along the segment, in [0, PPS).
 * @return Canvas row index, computed as m.y_base + i * m.y_step.
 */
constexpr int segment_y(const SegmentMap &m, int i) {
  return m.y_base + i * m.y_step;
}

/**
 * @brief Display clip rectangle a segment paints for one half-rev window.
 */
struct SegmentClip {
  int x0; /**< Inclusive left column. */
  int x1; /**< Exclusive right column. */
  int y0; /**< Inclusive top row. */
  int y1; /**< Exclusive bottom row. */
};

/**
 * @brief Canvas sub-rectangle a segment renders for the upcoming display window.
 * @param m Segment vertical mapping (arm side and strip direction).
 * @param arm_a_left True if the window this frame displays in sweeps arm-A
 * columns [0, w/2) (a ZERO->HALF half-rev); false for [w/2, w).
 * @param S Total LEDs across both arms (rows = S/2).
 * @param N Segment count.
 * @param w Canvas width.
 * @return Display clip rectangle; pair with ClipRegion's render margin for filters.
 * @details A buffer flips at each half-rev boundary, so an arm sweeps only half
 * the columns per displayed frame. The N segments tile the full canvas every
 * frame, arm A and arm B trading column halves each window; arm B samples
 * (x + w/2), so it paints the opposite half from arm A in the same window.
 */
constexpr SegmentClip segment_clip(const SegmentMap &m, bool arm_a_left, int S,
                                   int N, int w) {
  const int pps = S / N;
  const int y0 = m.y_step > 0 ? m.y_base : m.y_base - (pps - 1);
  const bool left = m.arm_b ? !arm_a_left : arm_a_left;
  const int x0 = left ? 0 : w / 2;
  return SegmentClip{x0, x0 + w / 2, y0, y0 + pps};
}

} // namespace pov
