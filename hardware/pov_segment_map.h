/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file pov_segment_map.h
 * @brief Pure segment index arithmetic for the segmented POV driver.
 *
 * Split out of pov_segmented.h (which is Arduino-only) so the load-bearing
 * index math — which segment owns which canvas column and rows, and in which
 * direction its physical strip runs — is unit-testable on the host without a
 * Teensy, exactly as hd107s_frame.h was split out of dma_led.h. The driver's
 * boot-time configure_segment() and per-column ISR derive their mapping from
 * these functions, so the host tests cover the real arithmetic. An off-by-one
 * here silently mis-paints the sphere.
 *
 * Layout (N=4, S=288 -> ROWS=144, PPS=72): two arms, two segments each.
 *   Arm A samples canvas column x; arm B samples column (x + W/2) % W.
 *   Per arm: the top segment runs y = 0 -> PPS-1 (N pole -> junction),
 *            the bottom segment runs y = ROWS-1 -> ROWS-PPS (S pole ->
 *            junction, strip physically reversed).
 */
#pragma once

namespace pov {

/// Precomputed vertical mapping for one physical segment.
struct SegmentMap {
  bool arm_b;  ///< Samples canvas column (x + W/2) instead of x.
  int y_base;  ///< Canvas row of this segment's LED 0.
  int y_step;  ///< +1 (top, toward junction) or -1 (bottom strip, reversed).
};

/**
 * @brief Derive a segment's vertical mapping from its hardware ID.
 * @param segment_id  Hardware-strapped ID in [0, N).
 * @param S  Total LEDs across both arms (ROWS = S/2).
 * @param N  Segment count (even, power of two, <= 4; SEGS_PER_ARM = N/2).
 *
 * Segments [0, N/2) are arm A; [N/2, N) are arm B. Within an arm, arm-segment 0
 * is the top strip (LED 0 at the N pole, y=0, counting toward the junction) and
 * arm-segment 1 is the bottom strip (LED 0 at the S pole, y=ROWS-1, counting
 * toward the junction — physically reversed).
 */
constexpr SegmentMap segment_map(int segment_id, int S, int N) {
  const int segs_per_arm = N / 2;
  const int rows = S / 2;
  const int arm_seg = segment_id % segs_per_arm;
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

/// Canvas column this segment samples at rotation column @p x (canvas width
/// @p w). Arm B reads the opposite half of the image.
constexpr int segment_x_col(bool arm_b, int x, int w) {
  return arm_b ? (x + w / 2) % w : x;
}

/// Canvas row of this segment's i-th LED (i in [0, PPS)).
constexpr int segment_y(const SegmentMap &m, int i) {
  return m.y_base + i * m.y_step;
}

} // namespace pov
