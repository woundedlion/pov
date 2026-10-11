/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file pov_segment_map.h
 * @brief Pure segment-ID and index arithmetic for the segmented POV driver.
 *
 * Layout: two arms, with N/2 contiguous row bands per arm.
 *   Arm A samples canvas column x; arm B samples column (x + W/2) % W.
 *   Northern bands run toward the junction in increasing y order. Southern
 *   bands run toward the junction in decreasing y order.
 *
 * NOTE: the top segment wires LED 0 at the N pole (y_step = +1), unlike the
 * single-board rig's reversed top arm (pov_single_map.h).
 */
#pragma once

namespace pov {

/**
 * @brief Number of active-low GPIO straps needed for a segment count.
 * @param segment_count Power-of-two segment count in [2, 8].
 * @return One, two, or three sampled strap bits.
 */
constexpr int segment_id_strap_count(int segment_count) {
  return segment_count <= 2 ? 1 : segment_count <= 4 ? 2 : 3;
}

/**
 * @brief Decode an active-low segment-ID strap sample.
 * @param raw Raw strap bits, with floating inputs HIGH and grounded inputs LOW.
 * @param segment_count Power-of-two segment count in [2, 8].
 * @return Segment ID in [0, segment_count).
 */
constexpr int decode_segment_id(int raw, int segment_count) {
  return (~raw) & (segment_count - 1);
}

/**
 * @brief Precomputed vertical mapping for one physical segment.
 */
struct SegmentMap {
  bool arm_b; /**< Samples canvas column (x + W/2) instead of x. */
  int y_base; /**< Canvas row of this segment's LED 0. */
  int y_step; /**< +1 (north band) or -1 (reversed south band). */
};

/**
 * @brief Truncated division of nonnegative integer coordinates.
 * @param numerator Nonnegative coordinate.
 * @param denominator Positive partition count.
 * @return Truncated integer quotient.
 */
constexpr int segment_quotient(int numerator, int denominator) {
  return numerator / denominator;
}
/**
 * @brief Floor division of nonnegative browser coordinates.
 * @param numerator Nonnegative coordinate.
 * @param denominator Positive partition count.
 * @return Floor of the coordinate quotient.
 */
inline double segment_quotient(double numerator, double denominator) {
  return __builtin_floor(numerator / denominator);
}
/**
 * @brief Integer coordinate remainder.
 * @param numerator Nonnegative coordinate.
 * @param denominator Positive partition count.
 * @return Integer coordinate remainder.
 */
constexpr int segment_remainder(int numerator, int denominator) {
  return numerator % denominator;
}
/**
 * @brief Browser coordinate remainder.
 * @param numerator Nonnegative coordinate.
 * @param denominator Positive partition count.
 * @return Floating-point remainder from fmod.
 */
inline double segment_remainder(double numerator, double denominator) {
  return __builtin_fmod(numerator, denominator);
}

/**
 * @brief Row-band index for a segment within one arm.
 * @param arm_segment Segment index within the arm.
 * @param bands Row bands per arm.
 * @return Northern bands in order, followed by southern bands in reverse order.
 */
template <typename Coordinate>
constexpr Coordinate segment_band(Coordinate arm_segment, Coordinate bands) {
  const Coordinate NORTH_BANDS = segment_quotient(bands, Coordinate{2});
  return arm_segment < NORTH_BANDS ? arm_segment
                                   : bands - 1 - (arm_segment - NORTH_BANDS);
}

/** @brief Browser worker partition of a two-arm canvas. */
template <typename Coordinate> struct SegmentRange {
  Coordinate x0; /**< Inclusive left column. */
  Coordinate x1; /**< Exclusive right column. */
  Coordinate y0; /**< Inclusive top row. */
  Coordinate y1; /**< Exclusive bottom row. */
};

/**
 * @brief Partitions an even-width canvas into two arms and row bands.
 * @param id Segment index in [0, total).
 * @param total Positive even segment count.
 * @param width Positive even canvas width.
 * @param height Canvas height, at least total/2 rows.
 * @return Segment rectangle; remainder rows fill bands from north to south.
 */
template <typename Coordinate>
constexpr SegmentRange<Coordinate>
segment_range(Coordinate id, Coordinate total, Coordinate width,
              Coordinate height) {
  const Coordinate BANDS = segment_quotient(total, Coordinate{2});
  const Coordinate BAND = segment_band(segment_remainder(id, BANDS), BANDS);
  const Coordinate ARM_WIDTH = width / 2;
  const Coordinate X0 = segment_quotient(id, BANDS) * ARM_WIDTH;
  const Coordinate ROWS = segment_quotient(height, BANDS);
  const Coordinate REMAINDER = segment_remainder(height, BANDS);
  const Coordinate Y0 = BAND * ROWS + (BAND < REMAINDER ? BAND : REMAINDER);
  return {X0, X0 + ARM_WIDTH, Y0, Y0 + ROWS + (BAND < REMAINDER ? 1 : 0)};
}

/**
 * @brief Derive a segment's vertical mapping from its hardware ID.
 * @param segment_id  Hardware-strapped ID in [0, N).
 * @param S  Total LEDs across both arms (ROWS = S/2).
 * @param N  Segment count (even, power of two, <= 8; SEGS_PER_ARM = N/2).
 * @return SegmentMap with arm_b, y_base (canvas row of LED 0), and y_step (+1/-1).
 * @details Segments [0, N/2) are arm A; [N/2, N) are arm B. Within an arm,
 * the first half tile the northern hemisphere in increasing y order. The
 * second half tile the southern hemisphere from the S pole toward the junction
 * in decreasing y order. N=2 assigns one forward strip to each complete arm.
 */
constexpr SegmentMap segment_map(int segment_id, int S, int N) {
  const int SEGS_PER_ARM = N / 2;
  const int PPS = S / N;
  // SEGS_PER_ARM is a power of two, so mask for the within-arm index.
  const int ARM_SEG = segment_id & (SEGS_PER_ARM - 1);
  SegmentMap m{};
  m.arm_b = segment_id >= SEGS_PER_ARM;

  const int SEGS_PER_HEMISPHERE = SEGS_PER_ARM / 2;
  if (SEGS_PER_ARM == 1 || ARM_SEG < SEGS_PER_HEMISPHERE) {
    m.y_base = segment_band(ARM_SEG, SEGS_PER_ARM) * PPS;
    m.y_step = 1;
  } else {
    m.y_base = (segment_band(ARM_SEG, SEGS_PER_ARM) + 1) * PPS - 1;
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
 * @brief Row-major display-buffer index of this segment's LED 0.
 * @param m Segment mapping providing y_base.
 * @param x_col Canvas column this segment samples (from segment_x_col).
 * @param w Canvas width in columns.
 * @return Index into a row-major w-wide buffer of the segment's first pixel.
 */
constexpr int segment_pixel_base(const SegmentMap &m, int x_col, int w) {
  return m.y_base * w + x_col;
}

/**
 * @brief Row-major display-buffer stride between consecutive segment LEDs.
 * @param m Segment mapping providing y_step.
 * @param w Canvas width in columns.
 * @return Signed index delta per LED; negative for a reversed southern band.
 * @details segment_pixel_base() plus i strides equals segment_y()'s row times
 * w plus x_col.
 */
constexpr int segment_row_stride(const SegmentMap &m, int w) {
  return m.y_step * w;
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

/**
 * @brief Whether a segment renders only its segment_clip() rectangle.
 * @param full_frame The effect's needs_full_frame().
 * @param persists The effect's persists_pixels().
 * @return False when the effect reads cross-segment or prior-frame state.
 */
constexpr bool segment_clip_applies(bool full_frame, bool persists) {
  return !full_frame && !persists;
}

} // namespace pov
