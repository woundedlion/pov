/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "math/3dmath.h"

#ifndef HS_RUNTIME_DISPLAY_GEOMETRY
#define HS_RUNTIME_DISPLAY_GEOMETRY 0
#endif

#ifndef HS_DISPLAY_PROFILE
#define HS_DISPLAY_PROFILE 1
#endif
#ifndef HS_DISPLAY_NORTH_FRACTION
#define HS_DISPLAY_NORTH_FRACTION 0.02f
#endif
#ifndef HS_DISPLAY_SOUTH_FRACTION
#define HS_DISPLAY_SOUTH_FRACTION 0.98f
#endif

namespace math {
static_assert(HS_DISPLAY_PROFILE == 0 || HS_DISPLAY_PROFILE == 1,
              "HS_DISPLAY_PROFILE must be 0 (ideal) or 1 (physical)");
#if HS_RUNTIME_DISPLAY_GEOMETRY
inline float DISPLAY_NORTH_PHI = 0.0f;
inline float DISPLAY_SOUTH_PHI = PI_F;

/** @brief Changes host latitude endpoints; refresh geometry LUTs before rendering. */
inline bool set_display_geometry(float north, float south) {
  if (!(north >= 0.0f && north <= PI_F * 0.25f && south >= PI_F * 0.75f &&
        south <= PI_F))
    return false;
  DISPLAY_NORTH_PHI = north;
  DISPLAY_SOUTH_PHI = south;
  return true;
}
#else
inline constexpr float DISPLAY_NORTH_PHI =
    HS_DISPLAY_PROFILE == 0 ? 0.0f : (HS_DISPLAY_NORTH_FRACTION * PI_F);
inline constexpr float DISPLAY_SOUTH_PHI =
    HS_DISPLAY_PROFILE == 0 ? PI_F : (HS_DISPLAY_SOUTH_FRACTION * PI_F);
#endif

/** @brief Uniform latitude mapping between the first and last LED centers. */
class LatitudeGeometry {
public:
#if defined(HS_TEST_H_OFFSET)
  constexpr explicit LatitudeGeometry(int height)
      : LatitudeGeometry(height, 0.0f,
                         HS_TEST_H_OFFSET == 0
                             ? PI_F
                             : PI_F * (height - 1) /
                                   (height + HS_TEST_H_OFFSET - 1)) {}
#elif HS_RUNTIME_DISPLAY_GEOMETRY
  explicit LatitudeGeometry(int height)
      : LatitudeGeometry(height, DISPLAY_NORTH_PHI, DISPLAY_SOUTH_PHI) {}
#else
  constexpr explicit LatitudeGeometry(int height)
      : LatitudeGeometry(height, DISPLAY_NORTH_PHI, DISPLAY_SOUTH_PHI) {}
#endif
  constexpr LatitudeGeometry(int height, float north, float south)
      : height(height), north(north), south(south) {
    HS_CHECK(height > 1 && north >= 0.0f && south <= PI_F && north < south,
             "Invalid latitude geometry");
  }
  constexpr float row_to_phi(float row) const {
    return north + row * radians_per_row();
  }
  constexpr float phi_to_row(float phi) const {
    return (phi - north) * rows_per_radian();
  }
  constexpr float radians_per_row() const {
    return (south - north) / (height - 1);
  }
  constexpr float rows_per_radian() const {
    return (height - 1) / (south - north);
  }
  constexpr bool contains_row(float row) const {
    return row >= 0.0f && row <= height - 1;
  }

private:
  int height;
  float north;
  float south;
};

/** @brief Display geometry; an explicit offset selects legacy mapping. */
template <int H, int LegacyOffset = -1> struct DisplayGeometry {
  static_assert(H > 1);
#if defined(HS_TEST_H_OFFSET)
  static constexpr int OFFSET =
      LegacyOffset < 0 ? HS_TEST_H_OFFSET : LegacyOffset;
#else
  static constexpr int OFFSET = LegacyOffset;
#endif
  static constexpr float NORTH_PHI = OFFSET >= 0 ? 0.0f : DISPLAY_NORTH_PHI;
  static constexpr float SOUTH_PHI = OFFSET == 0 ? PI_F
                                     : OFFSET > 0
                                         ? PI_F *(H - 1) / (H + OFFSET - 1)
                                         : DISPLAY_SOUTH_PHI;
  static_assert(NORTH_PHI >= 0.0f && SOUTH_PHI <= PI_F &&
                NORTH_PHI < SOUTH_PHI);
  static constexpr float RADIANS_PER_ROW = (SOUTH_PHI - NORTH_PHI) / (H - 1);
  static constexpr float ROWS_PER_RADIAN = 1.0f / RADIANS_PER_ROW;
  static constexpr float NORTH_POLE_ROW = -NORTH_PHI * ROWS_PER_RADIAN;
  static constexpr float SOUTH_POLE_ROW = (PI_F - NORTH_PHI) * ROWS_PER_RADIAN;
  static constexpr bool HAS_NORTH_POLE = NORTH_PHI == 0.0f;
  static constexpr bool HAS_SOUTH_POLE = SOUTH_PHI == PI_F;
  static constexpr float row_to_phi(float row) {
    return NORTH_PHI + row * RADIANS_PER_ROW;
  }
  static constexpr float phi_to_row(float phi) {
    return (phi - NORTH_PHI) * ROWS_PER_RADIAN;
  }
  static constexpr bool contains_row(float row) {
    return row >= 0.0f && row <= H - 1;
  }
};
#if HS_RUNTIME_DISPLAY_GEOMETRY && !defined(HS_TEST_H_OFFSET)
template <int H> struct DisplayGeometry<H, -1> {
  static_assert(H > 1);
  static constexpr int OFFSET = -1;
  inline static float NORTH_PHI = 0.0f;
  inline static float SOUTH_PHI = PI_F;
  inline static float RADIANS_PER_ROW = PI_F / (H - 1);
  inline static float ROWS_PER_RADIAN = (H - 1) / PI_F;
  inline static float NORTH_POLE_ROW = 0.0f;
  inline static float SOUTH_POLE_ROW = H - 1;
  inline static bool HAS_NORTH_POLE = true;
  inline static bool HAS_SOUTH_POLE = true;

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
  static float row_to_phi(float row) {
    return NORTH_PHI + row * RADIANS_PER_ROW;
  }
  static float phi_to_row(float phi) {
    return (phi - NORTH_PHI) * ROWS_PER_RADIAN;
  }
  static constexpr bool contains_row(float row) {
    return row >= 0.0f && row <= H - 1;
  }
};
#endif
} // namespace math
