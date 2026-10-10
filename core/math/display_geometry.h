/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file display_geometry.h
 * @brief Unclamped display latitude and pole-row mappings. */

#include "math/3dmath.h"

#ifndef HS_RUNTIME_DISPLAY_GEOMETRY
#define HS_RUNTIME_DISPLAY_GEOMETRY 0
#endif

#ifndef HS_DISPLAY_PROFILE
#define HS_DISPLAY_PROFILE 1
#endif
#ifndef HS_DISPLAY_NORTH_FRACTION
/// Physical profile: first LED row's polar angle as a fraction of pi.
#define HS_DISPLAY_NORTH_FRACTION 0.02f
#endif
#ifndef HS_DISPLAY_SOUTH_FRACTION
/// Physical profile: last LED row's polar angle as a fraction of pi.
#define HS_DISPLAY_SOUTH_FRACTION 0.98f
#endif

namespace math {
static_assert(HS_DISPLAY_PROFILE == 0 || HS_DISPLAY_PROFILE == 1,
              "HS_DISPLAY_PROFILE must be 0 (ideal) or 1 (physical)");
#if HS_RUNTIME_DISPLAY_GEOMETRY
/// Polar angle of the first LED row, radians; set by set_display_geometry().
inline float DISPLAY_NORTH_PHI = 0.0f;
/// Polar angle of the last LED row, radians; set by set_display_geometry().
inline float DISPLAY_SOUTH_PHI = PI_F;

/**
 * @brief Changes host latitude endpoints; refresh geometry LUTs before
 *        rendering.
 * @param north Polar angle of the first LED row, radians, in [0, pi/4].
 * @param south Polar angle of the last LED row, radians, in [3pi/4, pi].
 * @return False, leaving the endpoints unchanged, when either is out of range.
 */
inline bool set_display_geometry(float north, float south) {
  if (!(north >= 0.0f && north <= PI_F * 0.25f && south >= PI_F * 0.75f &&
        south <= PI_F))
    return false;
  DISPLAY_NORTH_PHI = north;
  DISPLAY_SOUTH_PHI = south;
  return true;
}
#else
/// Polar angle of the first LED row under the active profile, radians.
inline constexpr float DISPLAY_NORTH_PHI =
    HS_DISPLAY_PROFILE == 0 ? 0.0f : (HS_DISPLAY_NORTH_FRACTION * PI_F);
/// Polar angle of the last LED row under the active profile, radians.
inline constexpr float DISPLAY_SOUTH_PHI =
    HS_DISPLAY_PROFILE == 0 ? PI_F : (HS_DISPLAY_SOUTH_FRACTION * PI_F);
#endif

/** @brief Uniform latitude mapping between the first and last LED centers. */
class LatitudeGeometry {
public:
#if defined(HS_TEST_H_OFFSET)
  /**
   * @brief Constructs a latitude mapping from the active display profile.
   * @param height LED row count; > 1.
   */
  constexpr explicit LatitudeGeometry(int height)
      : LatitudeGeometry(height, 0.0f,
                         HS_TEST_H_OFFSET == 0
                             ? PI_F
                             : PI_F * (height - 1) /
                                   (height + HS_TEST_H_OFFSET - 1)) {}
#elif HS_RUNTIME_DISPLAY_GEOMETRY
  /**
   * @brief Constructs a latitude mapping from the active display profile.
   * @param height LED row count; > 1.
   */
  explicit LatitudeGeometry(int height)
      : LatitudeGeometry(height, DISPLAY_NORTH_PHI, DISPLAY_SOUTH_PHI) {}
#else
  /**
   * @brief Constructs a latitude mapping from the active display profile.
   * @param height LED row count; > 1.
   */
  constexpr explicit LatitudeGeometry(int height)
      : LatitudeGeometry(height, DISPLAY_NORTH_PHI, DISPLAY_SOUTH_PHI) {}
#endif
  /**
   * @brief Constructs a mapping with explicit latitude endpoints.
   * @param height LED row count; > 1.
   * @param north Polar angle of the first row, radians.
   * @param south Polar angle of the last row, radians; greater than `north`.
   */
  constexpr LatitudeGeometry(int height, float north, float south)
      : height(height), north(north), south(south) {
    HS_CHECK(height > 1 && north >= 0.0f && south <= PI_F && north < south,
             "Invalid latitude geometry");
  }
  /**
   * @brief Maps an unclamped row coordinate to latitude radians.
   * @param row Row coordinate; may lie outside [0, height - 1].
   * @return Polar angle, radians.
   */
  constexpr float row_to_phi(float row) const {
    return north + row * radians_per_row();
  }
  /**
   * @brief Maps latitude to an unclamped row; pole rows may lie outside the
   *        display.
   * @param phi Polar angle, radians.
   * @return Fractional row.
   */
  constexpr float phi_to_row(float phi) const {
    return (phi - north) * rows_per_radian();
  }
  /**
   * @brief Returns the latitude spacing between adjacent LED rows.
   * @return Radians per row.
   */
  constexpr float radians_per_row() const {
    return (south - north) / (height - 1);
  }
  /**
   * @brief Returns the inverse latitude spacing.
   * @return Rows per radian.
   */
  constexpr float rows_per_radian() const {
    return (height - 1) / (south - north);
  }
  /**
   * @brief Tests whether a row lies in the displayed center range.
   * @param row Row coordinate.
   * @return True when `row` is in [0, height - 1].
   */
  constexpr bool contains_row(float row) const {
    return row >= 0.0f && row <= height - 1;
  }

private:
  int height;
  float north;
  float south;
};

/**
 * @brief Display geometry; a non-negative HOffset selects the test pole-to-pole
 * mapping with that many virtual south rows; -1 selects `HS_TEST_H_OFFSET` when
 * defined, else the active display profile.
 */
template <int H, int HOffset = -1> struct DisplayGeometry {
  static_assert(H > 1);
#if defined(HS_TEST_H_OFFSET)
  /// Virtual south rows of the test mapping; -1 for the display profile.
  static constexpr int OFFSET = HOffset < 0 ? HS_TEST_H_OFFSET : HOffset;
#else
  /// Virtual south rows of the test mapping; -1 for the display profile.
  static constexpr int OFFSET = HOffset;
#endif
  /// Polar angle of row 0, radians.
  static constexpr float NORTH_PHI = OFFSET >= 0 ? 0.0f : DISPLAY_NORTH_PHI;
  /// Polar angle of row H-1, radians.
  static constexpr float SOUTH_PHI = OFFSET == 0 ? PI_F
                                     : OFFSET > 0
                                         ? PI_F *(H - 1) / (H + OFFSET - 1)
                                         : DISPLAY_SOUTH_PHI;
  static_assert(NORTH_PHI >= 0.0f && SOUTH_PHI <= PI_F &&
                NORTH_PHI < SOUTH_PHI);
  /// Polar-angle step between adjacent rows, radians.
  static constexpr float RADIANS_PER_ROW = (SOUTH_PHI - NORTH_PHI) / (H - 1);
  /// Reciprocal of `RADIANS_PER_ROW`.
  static constexpr float ROWS_PER_RADIAN = 1.0f / RADIANS_PER_ROW;
  /** @brief Unclamped north-pole row; negative when the north cap is excluded. */
  static constexpr float NORTH_POLE_ROW = -NORTH_PHI * ROWS_PER_RADIAN;
  /** @brief Unclamped south-pole row; above H-1 when the south cap is excluded. */
  static constexpr float SOUTH_POLE_ROW = (PI_F - NORTH_PHI) * ROWS_PER_RADIAN;
  /** @brief Whether the first LED row reaches the north pole. */
  static constexpr bool HAS_NORTH_POLE = NORTH_PHI == 0.0f;
  /** @brief Whether the last LED row reaches the south pole. */
  static constexpr bool HAS_SOUTH_POLE = SOUTH_PHI == PI_F;
  /**
   * @brief Maps an unclamped row coordinate to latitude radians.
   * @param row Row coordinate; may lie outside [0, H - 1].
   * @return Polar angle, radians.
   */
  static constexpr float row_to_phi(float row) {
    return NORTH_PHI + row * RADIANS_PER_ROW;
  }
  /**
   * @brief Maps latitude to an unclamped row; pole rows may lie outside the
   *        display.
   * @param phi Polar angle, radians.
   * @return Fractional row.
   */
  static constexpr float phi_to_row(float phi) {
    return (phi - NORTH_PHI) * ROWS_PER_RADIAN;
  }
  /**
   * @brief Tests whether a row lies in the displayed center range.
   * @param row Row coordinate.
   * @return True when `row` is in [0, H - 1].
   */
  static constexpr bool contains_row(float row) {
    return row >= 0.0f && row <= H - 1;
  }
};
#if HS_RUNTIME_DISPLAY_GEOMETRY && !defined(HS_TEST_H_OFFSET)
/**
 * @brief Runtime display geometry; refresh() reloads the latitude endpoints
 *        from `DISPLAY_NORTH_PHI` and `DISPLAY_SOUTH_PHI`.
 * @tparam H Display height in rows.
 */
template <int H> struct DisplayGeometry<H, -1> {
  static_assert(H > 1);
  static constexpr int OFFSET = -1;     ///< Always the display profile.
  inline static float NORTH_PHI = 0.0f; ///< Polar angle of row 0, radians.
  inline static float SOUTH_PHI = PI_F; ///< Polar angle of row H-1, radians.
  /// Polar-angle step between adjacent rows, radians.
  inline static float RADIANS_PER_ROW = PI_F / (H - 1);
  /// Reciprocal of `RADIANS_PER_ROW`.
  inline static float ROWS_PER_RADIAN = (H - 1) / PI_F;
  /** @brief Unclamped north-pole row; negative when the north cap is excluded. */
  inline static float NORTH_POLE_ROW = 0.0f;
  /** @brief Unclamped south-pole row; above H-1 when the south cap is excluded. */
  inline static float SOUTH_POLE_ROW = H - 1;
  /** @brief Whether the first LED row reaches the north pole. */
  inline static bool HAS_NORTH_POLE = true;
  /** @brief Whether the last LED row reaches the south pole. */
  inline static bool HAS_SOUTH_POLE = true;

  /** @brief Refreshes the initial ideal mapping after set_display_geometry(). */
  static void refresh() {
    NORTH_PHI = DISPLAY_NORTH_PHI;
    SOUTH_PHI = DISPLAY_SOUTH_PHI;
    RADIANS_PER_ROW = (SOUTH_PHI - NORTH_PHI) / (H - 1);
    ROWS_PER_RADIAN = 1.0f / RADIANS_PER_ROW;
    NORTH_POLE_ROW = -NORTH_PHI * ROWS_PER_RADIAN;
    SOUTH_POLE_ROW = (PI_F - NORTH_PHI) * ROWS_PER_RADIAN;
    HAS_NORTH_POLE = NORTH_PHI == 0.0f;
    HAS_SOUTH_POLE = SOUTH_PHI == PI_F;
  }
  /**
   * @brief Maps an unclamped row coordinate to latitude radians.
   * @param row Row coordinate; may lie outside [0, H - 1].
   * @return Polar angle, radians.
   */
  static float row_to_phi(float row) {
    return NORTH_PHI + row * RADIANS_PER_ROW;
  }
  /**
   * @brief Maps latitude to an unclamped row; pole rows may lie outside the
   *        display.
   * @param phi Polar angle, radians.
   * @return Fractional row.
   */
  static float phi_to_row(float phi) {
    return (phi - NORTH_PHI) * ROWS_PER_RADIAN;
  }
  /**
   * @brief Tests whether a row lies in the displayed center range.
   * @param row Row coordinate.
   * @return True when `row` is in [0, H - 1].
   */
  static constexpr bool contains_row(float row) {
    return row >= 0.0f && row <= H - 1;
  }
};
#endif
} // namespace math
