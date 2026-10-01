/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file pixel_mapping.h
 * @brief Resolution-dependent sphere/canvas mapping and lookup tables.
 */

#include "math/3dmath.h"
#include "math/periodic.h"
#include "math/display_geometry.h"
#include <array>

namespace math {
/**
 * @brief Structure representing 2D floating-point pixel coordinates.
 */
struct PixelCoords {
  float x; /**< Horizontal coordinate. */
  float y; /**< Vertical coordinate. */
};

/**
 * @brief Converts a pixel y-coordinate to a spherical phi angle.
 * @param y The pixel y-coordinate [0, h_virt - 1].
 * @param h_virt Height of a complete pole-to-pole latitude grid.
 * @return The spherical phi angle in radians.
 */
inline float y_to_phi_virtual(float y, int h_virt) {
  HS_CHECK(h_virt > 1, "y_to_phi_virtual: h_virt must be > 1");
  return (y * PI_F) / (h_virt - 1);
}

/**
 * @brief Converts a spherical phi angle to a pixel y-coordinate.
 * @param phi The spherical phi angle in radians.
 * @param h_virt Height of a complete pole-to-pole latitude grid.
 * @return The pixel y-coordinate in [0, h_virt - 1] for phi in [0, pi], EXCEPT at
 *   the south pole (phi == PI_F) the float round-trip can land a hair *above*
 *   `h_virt - 1`; a caller indexing a row buffer with `(int)y` must clamp or
 *   floor first.
 */
inline float phi_to_y_virtual(float phi, int h_virt) {
  HS_CHECK(h_virt > 1, "phi_to_y_virtual: h_virt must be > 1");
  return (phi * (h_virt - 1)) / PI_F;
}

/** @brief Polar angle to fractional display row. */
template <int H> inline float phi_to_y(float phi) {
  return DisplayGeometry<H>::phi_to_row(phi);
}

#if HS_RUNTIME_DISPLAY_GEOMETRY
template <int H>
inline const float &ROWS_PER_RADIAN = DisplayGeometry<H>::ROWS_PER_RADIAN;
template <int H>
inline const float &RADIANS_PER_ROW = DisplayGeometry<H>::RADIANS_PER_ROW;
#else
template <int H>
inline constexpr float ROWS_PER_RADIAN = DisplayGeometry<H>::ROWS_PER_RADIAN;
template <int H>
inline constexpr float RADIANS_PER_ROW = DisplayGeometry<H>::RADIANS_PER_ROW;

#endif

/** @brief Radians of azimuth spanned by one canvas column. */
template <int W>
inline constexpr float RADIANS_PER_COLUMN = TWO_PI_F / static_cast<float>(W);

/** @brief Larger angular pitch of a logical canvas's rows and columns. */
template <int W, int H> constexpr float coarse_pixel_pitch() {
  return RADIANS_PER_COLUMN<W> > RADIANS_PER_ROW<H> ? RADIANS_PER_COLUMN<W>
                                                    : RADIANS_PER_ROW<H>;
}

/**
 * @brief Precomputed lookup table for scanline phi angles.
 * @tparam H Display height; legacy test profiles may append virtual rows.
 */
template <int H> struct PhiLUT {
  static constexpr int H_VIRT =
      H + (DisplayGeometry<H>::OFFSET > 0 ? DisplayGeometry<H>::OFFSET : 0);
  static std::array<float, H_VIRT> data; /**< phi per display row, radians. */
  // Lazy-init check-then-set is non-atomic; safe only because rendering is
  // single-threaded, NOT a concurrency guard.
  static bool initialized; /**< Lazy-init guard; true once data is filled. */
  /**
   * @brief Fills the phi table for every display row and marks it initialized.
   */
  static void init() {
#if HS_RUNTIME_DISPLAY_GEOMETRY && !defined(HS_TEST_H_OFFSET)
    DisplayGeometry<H>::refresh();
#endif
    for (int y = 0; y < H_VIRT; y++) {
      data[y] = DisplayGeometry<H>::row_to_phi(static_cast<float>(y));
    }
    initialized = true;
  }
};

template <int H> std::array<float, PhiLUT<H>::H_VIRT> PhiLUT<H>::data;
template <int H> bool PhiLUT<H>::initialized = false;

/**
 * @brief LUT-backed pixel-y -> phi for integer rows at compile-time height H.
 * @tparam H Logical height selecting the PhiLUT<H> table.
 * @param y Integer pixel row in [0, H_VIRT).
 * @return The spherical phi angle in radians for that row.
 * @details Lazily fills PhiLUT on first touch; traps an out-of-range row via
 * HS_CHECK.
 */
template <int H> inline float y_to_phi(int y) {
  if (!PhiLUT<H>::initialized) {
    PhiLUT<H>::init();
  }
  HS_CHECK(y >= 0 && y < PhiLUT<H>::H_VIRT, "y_to_phi: row %d outside [0, %d)",
           y, PhiLUT<H>::H_VIRT);
  return PhiLUT<H>::data[y];
}

/** @brief Fractional display row to polar angle, including extrapolated cap rows. */
template <int H> inline float y_to_phi(float y) {
  return DisplayGeometry<H>::row_to_phi(y);
}

/**
 * @brief Split trig lookup tables for efficient vector reconstruction.
 * @tparam W Width (column count).
 * @tparam H Display height; tables include legacy virtual rows only in tests.
 * @details Caches sin/cos for theta (per column) and phi (per row) separately,
 * reconstructing vectors with 3 multiplies. Memory: 1.25*W + 2*H_VIRT floats
 * (sin_theta carries the folded quarter turn) vs W*H_VIRT Vectors — a ~190x
 * reduction at 288x144.
 */
template <int W, int H> struct TrigLUT {
  static_assert(W % 4 == 0,
                "cos_theta is recovered as sin_theta[x + W/4]; W must be a "
                "multiple of 4 for the quarter-turn offset to be exact");
  static constexpr int H_VIRT =
      H + (DisplayGeometry<H>::OFFSET > 0 ? DisplayGeometry<H>::OFFSET : 0);
  // sin_theta carries W/4 extra trailing entries (one quarter turn) so cos(theta)
  // reads back as sin_theta[x + W/4], avoiding a separate cos table.
  static constexpr int W_EXT = W + W / 4;
  static std::array<float, W_EXT> sin_theta; /**< sin(theta); cos via +W/4. */
  static std::array<float, H_VIRT> sin_phi;  /**< sin(phi) per virtual row. */
  static std::array<float, H_VIRT> cos_phi;  /**< cos(phi) per virtual row. */
  static bool initialized; /**< Lazy-init guard; true once tables are filled. */
  /**
   * @brief cos(theta) for column x, recovered from the extended sin table.
   * @param x Column in [0, W). Returns sin_theta[x + W/4] == cos(x*2*pi/W).
   */
  static float cos_theta(int x) {
    assert(x >= 0 && x < W);
    return sin_theta[x + W / 4];
  }
  /**
   * @brief Fills the theta and phi tables and marks them initialized.
   * @details Ensures PhiLUT<H> is populated first to source the phi angles. The
   * extra W/4 sin_theta entries wrap naturally: sin is 2*pi-periodic, so
   * sin(x*2*pi/W) for x in [W, W+W/4) equals the first-quarter values.
   */
  static void init() {
    if (!PhiLUT<H>::initialized) {
      PhiLUT<H>::init();
    }
    for (int x = 0; x < W_EXT; x++) {
      sin_theta[x] = sinf((x * 2 * PI_F) / W);
    }
    for (int y = 0; y < H_VIRT; y++) {
      float phi = PhiLUT<H>::data[y];
      sin_phi[y] = sinf(phi);
      cos_phi[y] = cosf(phi);
    }
    initialized = true;
  }
};

template <int W, int H>
std::array<float, TrigLUT<W, H>::W_EXT> TrigLUT<W, H>::sin_theta;
template <int W, int H>
std::array<float, TrigLUT<W, H>::H_VIRT> TrigLUT<W, H>::sin_phi;
template <int W, int H>
std::array<float, TrigLUT<W, H>::H_VIRT> TrigLUT<W, H>::cos_phi;
template <int W, int H> bool TrigLUT<W, H>::initialized = false;

/**
 * @brief Eagerly fill the scanline LUTs for resolution <W, H>.
 * @tparam W Width (column count).
 * @tparam H Logical height.
 * @details Engine setup calls this once before the first frame so the tables are
 * populated before any rendering — and, on hardware, before the column-sweep ISR
 * could observe a partially-filled table. The per-call `if (!initialized) init()`
 * guard in the per-pixel leaf `pixel_to_vector<W, H>(int, int)` remains a lazy
 * fallback. `y_to_phi<H>(int)` provides lazy initialization for tests and tools.
 * Their non-atomic check-then-set relies on this
 * eager call and the single-render-thread assumption. Idempotent.
 */
template <int W, int H> inline void init_geometry_luts() {
  PhiLUT<H>::init();
  TrigLUT<W, H>::init();
}

/**
 * @brief Recovers an effect's compile-time <W, H> from its type so a driver's
 * `show<E>()` can eager-init the LUTs without the caller restating the
 * resolution.
 * @tparam E The effect type, of the form `Eff<W, H>`.
 * @details Every effect is `template <int W, int H> class E`, so the partial
 * specialization matches them all.
 */
template <typename E> struct GeometryResolution;
/**
 * @brief Partial specialization that destructures an effect's <W, H>.
 * @tparam Eff The effect class template.
 * @tparam W Width recovered from the effect type.
 * @tparam H Height recovered from the effect type.
 */
template <template <int, int> class Eff, int W, int H>
struct GeometryResolution<Eff<W, H>> {
  /**
   * @brief Eager-inits the geometry LUTs for the recovered <W, H>.
   */
  static void init() { init_geometry_luts<W, H>(); }
};

/**
 * @brief Reconstruct a vector from pixel coordinates using split trig LUTs.
 * @tparam W Width.
 * @tparam H Height.
 * @param x X coordinate (column).
 * @param y Y coordinate (row).
 * @return Unit vector on the sphere.
 */
template <int W, int H> Vector pixel_to_vector(int x, int y) {
  if (!TrigLUT<W, H>::initialized) {
    TrigLUT<W, H>::init();
  }
  // Local avoids the comma in TrigLUT<W, H> inside the assert macro.
  [[maybe_unused]] constexpr int H_VIRT = TrigLUT<W, H>::H_VIRT;
  assert(x >= 0 && x < W && y >= 0 && y < H_VIRT);
  float sp = TrigLUT<W, H>::sin_phi[y];
  return Vector(sp * TrigLUT<W, H>::cos_theta(x), TrigLUT<W, H>::cos_phi[y],
                sp * TrigLUT<W, H>::sin_theta[x]);
}

/**
 * @brief Reconstruct a unit vector from fractional pixel coordinates.
 * @tparam W Width.
 * @tparam H Height.
 * @param x Fractional X coordinate (column).
 * @param y Fractional display row; may extrapolate into either cap.
 * @return Unit vector on the sphere.
 * @details Snaps to the integer LUT path when both coordinates are near-integer
 * AND in LUT range; otherwise builds the vector analytically from spherical
 * angles. An out-of-range `x` (e.g. an unwrapped `W`) is exact on the analytic
 * branch, since theta = 2*pi*x/W is periodic.
 */
template <int W, int H> Vector pixel_to_vector(float x, float y) {
  const float fx = std::floor(x);
  const float fy = std::floor(y);
  if (std::abs(x - fx) < TOLERANCE && std::abs(y - fy) < TOLERANCE) {
    if (fx >= 0 && fx < W && fy >= 0 && fy < TrigLUT<W, H>::H_VIRT) {
      return pixel_to_vector<W, H>(static_cast<int>(fx), static_cast<int>(fy));
    }
  }
  return Vector(Spherical((x * 2 * PI_F) / W, y_to_phi<H>(y)));
}

/**
 * @brief Projects a unit vector to its pixel column (azimuth only).
 * @tparam W The width.
 * @param v Unit vector on the sphere; only its x/z azimuth is read. Must be
 *   finite: `fast_atan2` has no NaN guard, and neither wrap branch fires on the
 *   NaN it propagates.
 * @return The `x` pixel coordinate in `[0, W)` (strictly excludes W); NaN for a
 *   non-finite `v`.
 */
template <int W>
__attribute__((always_inline)) inline float vector_to_theta(const Vector &v) {
  // fast_atan2 is bounded by |pi|, so t lands in [-W/2, W/2] and one conditional
  // add wraps it. The upper guard keeps the half-open range when a tiny negative
  // t rounds up to exactly W.
  float t = (fast_atan2(v.z, v.x) * W) / (2 * PI_F);
  if (t < 0.0f)
    t += W;
  return (t >= W) ? 0.0f : t;
}

/**
 * @brief Converts a 3D unit vector back to 2D pixel coordinates.
 * @note Derives `theta`/`phi` with the approximate `fast_atan2`/`fast_acos`, so
 *   the projection is sub-pixel inexact and `vector → pixel → vector` does not
 *   bit-exactly invert the exact-trig `pixel_to_vector`.
 * @tparam W The width.
 * @tparam H The height.
 * @param v The input vector; MUST be unit length (unenforced): `phi = acos(v.y)`
 *   is the true latitude only when |v| == 1, so a non-unit `v` returns a
 *   silently-wrong row. Unguarded per-pixel path; callers normalize first.
 * @return Pixel coordinates; rows in missing caps lie outside [0, H-1].
 */
template <int W, int H> HS_O3_FN PixelCoords vector_to_pixel(const Vector &v) {
  // phi = acos(v.y) is the true latitude only when |v| == 1; trap non-unit v in debug.
  assert(std::fabs(dot(v, v) - 1.0f) < math::EPS_UNIT_VEC_SQ);
  float phi = fast_acos(hs::clamp(v.y, -1.0f, 1.0f));
  PixelCoords p({vector_to_theta<W>(v), phi_to_y<H>(phi)});
  return p;
}

/** @brief Reflects fractional coordinates across the true poles, preserving missing caps. */
template <int W, int H, int HOffset = -1>
HS_O3_FN bool pole_wrap(float &col, float &row) {
  using Geometry = DisplayGeometry<H, HOffset>;
  if (row < Geometry::NORTH_POLE_ROW) {
    row = 2.0f * Geometry::NORTH_POLE_ROW - row;
    col += W * 0.5f;
  } else if (row > Geometry::SOUTH_POLE_ROW) {
    row = 2.0f * Geometry::SOUTH_POLE_ROW - row;
    col += W * 0.5f;
  }
  col = wrap(col, static_cast<float>(W));
  return Geometry::contains_row(row);
}

/** @brief Resolves a lattice tap only when reflection lands on the lattice.
 * @pre The column is already in [0, W).
 */
template <int W, int H, int HOffset = -1>
HS_O3_FN bool pole_wrap(int &col, int &row) {
  if (row >= 0 && row < H)
    return true;
  float x = static_cast<float>(col), y = static_cast<float>(row);
  if (!pole_wrap<W, H, HOffset>(x, y))
    return false;
  const float ix = std::round(x), iy = std::round(y);
  if (std::fabs(x - ix) > 0.0001f || std::fabs(y - iy) > 0.0001f)
    return false;
  col = fast_wrap(static_cast<int>(ix), W);
  row = static_cast<int>(iy);
  return row >= 0 && row < H;
}

} // namespace math
