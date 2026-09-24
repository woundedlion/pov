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
 * @param h_virt The VIRTUAL height, already including hs::H_OFFSET. The
 *   `y_to_phi<H>` template takes the LOGICAL height instead and adds the offset
 *   itself; the `_virtual` suffix keeps the two conventions apart at a call.
 * @return The spherical phi angle in radians.
 */
inline float y_to_phi_virtual(float y, int h_virt) {
  HS_CHECK(h_virt > 1, "y_to_phi_virtual: h_virt must be > 1");
  return (y * PI_F) / (h_virt - 1);
}

/**
 * @brief Converts a spherical phi angle to a pixel y-coordinate.
 * @param phi The spherical phi angle in radians.
 * @param h_virt The VIRTUAL height, already including hs::H_OFFSET (see
 *   y_to_phi_virtual).
 * @return The pixel y-coordinate in [0, h_virt - 1] for phi in [0, pi], EXCEPT at
 *   the south pole (phi == PI_F) the float round-trip can land a hair *above*
 *   `h_virt - 1`; a caller indexing a row buffer with `(int)y` must clamp or
 *   floor first.
 */
inline float phi_to_y_virtual(float phi, int h_virt) {
  HS_CHECK(h_virt > 1, "phi_to_y_virtual: h_virt must be > 1");
  return (phi * (h_virt - 1)) / PI_F;
}

/**
 * @brief phi -> pixel-y for a compile-time logical height H.
 * @tparam H Logical (not virtual) height; H_VIRT is derived as H + hs::H_OFFSET.
 * @param phi The spherical phi angle in radians.
 * @return The pixel y-coordinate.
 * @details Derives H_VIRT from H plus hs::H_OFFSET so callers pass the logical
 * height, not the virtual one.
 */
template <int H> inline float phi_to_y(float phi) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  static_assert(H_VIRT > 1, "phi<->y mapping degenerates when H_VIRT <= 1");
  return (phi * (H_VIRT - 1)) / PI_F;
}

/**
 * @brief Pixel rows spanned by one radian of phi at logical height H.
 * @tparam H Logical (not virtual) height; H_OFFSET is added here.
 * @details Dropping H_OFFSET from an open-coded row pitch is invisible on the
 * host, whose H_OFFSET is 0, and skews every device row.
 */
template <int H>
inline constexpr float ROWS_PER_RADIAN =
    static_cast<float>(H + hs::H_OFFSET - 1) / PI_F;

/**
 * @brief Radians of phi spanned by one pixel row at logical height H.
 * @tparam H Logical (not virtual) height; H_OFFSET is added here.
 */
template <int H>
inline constexpr float RADIANS_PER_ROW =
    PI_F / static_cast<float>(H + hs::H_OFFSET - 1);

/**
 * @brief Precomputed lookup table for scanline phi angles.
 * @tparam H Logical height; the table has H_VIRT = H + hs::H_OFFSET entries.
 */
template <int H> struct PhiLUT {
  static constexpr int H_VIRT = H + hs::H_OFFSET;
  static std::array<float, H_VIRT> data; /**< phi per virtual row, radians. */
  // Lazy-init check-then-set is non-atomic; safe only because rendering is
  // single-threaded, NOT a concurrency guard.
  static bool initialized; /**< Lazy-init guard; true once data is filled. */
  /**
   * @brief Fills the phi table for every virtual row and marks it initialized.
   */
  static void init() {
    for (int y = 0; y < H_VIRT; y++) {
      data[y] = y_to_phi_virtual(static_cast<float>(y), H_VIRT);
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

/**
 * @brief Pixel-y -> phi for fractional rows at compile-time height H.
 * @tparam H Logical height; H_VIRT is H + hs::H_OFFSET.
 * @param y Fractional pixel row.
 * @return The spherical phi angle in radians, from the same expression the
 *         PhiLUT rows are filled with.
 * @details The range is debug-asserted, so only NDEBUG builds extrapolate
 * linearly past [0, H_VIRT-1]; keeping y in range is the caller's
 * responsibility.
 */
template <int H> inline float y_to_phi(float y) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  static_assert(H_VIRT > 1, "phi<->y mapping degenerates when H_VIRT <= 1");
  assert(y >= 0.0f && y <= H_VIRT - 1);
  return (y * PI_F) / (H_VIRT - 1);
}

/**
 * @brief Split trig lookup tables for efficient vector reconstruction.
 * @tparam W Width (column count).
 * @tparam H Logical height; phi tables have H_VIRT = H + hs::H_OFFSET entries.
 * @details Caches sin/cos for theta (per column) and phi (per row) separately,
 * reconstructing vectors with 3 multiplies. Memory: 1.25*W + 2*H_VIRT floats
 * (sin_theta carries the folded quarter turn) vs W*H_VIRT Vectors — a ~190x
 * reduction at 288x144.
 */
template <int W, int H> struct TrigLUT {
  static_assert(W % 4 == 0,
                "cos_theta is recovered as sin_theta[x + W/4]; W must be a "
                "multiple of 4 for the quarter-turn offset to be exact");
  static constexpr int H_VIRT = H + hs::H_OFFSET;
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
 * guards in the per-pixel leaves `y_to_phi<H>(int)` and
 * `pixel_to_vector<W, H>(int, int)` then remain only as a lazy fallback for
 * unit tests and offline tools; their non-atomic check-then-set relies on this
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
 * @param y Fractional Y coordinate (row); the analytic branch passes it straight
 *   to `y_to_phi<H>`, which debug-asserts the range. Under NDEBUG a sub-pixel
 *   `y` outside [0, H_VIRT-1] extrapolates phi past [0, pi] by analytic
 *   continuity; callers must keep `y` in range (this is a per-pixel path, so no
 *   clamp).
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
  // y_to_phi<H> already accounts for H_OFFSET internally; pass H, not H_VIRT.
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
 * @return The 2D PixelCoords. The `y` field is a float in `[0, H_VIRT-1]` but at
 *   the south pole can land a hair *above* `H_VIRT-1` (float round-trip), while
 *   `x` is in `[0, W)` (strictly excludes W); a caller indexing a row/column
 *   buffer must floor (not round) first. Only `y` carries the clamp's NaN->hi
 *   guard: a NaN `v.y` clamps to +1, so `y` saturates to row 0 (the north
 *   pole); a -inf `v.y` clamps to -1 and lands on the south pole. `x` stays
 *   NaN either way.
 */
template <int W, int H> HS_O3_FN PixelCoords vector_to_pixel(const Vector &v) {
  // phi = acos(v.y) is the true latitude only when |v| == 1; trap non-unit v in debug.
  assert(std::fabs(dot(v, v) - 1.0f) < math::EPS_UNIT_VEC_SQ);
  float phi = fast_acos(hs::clamp(v.y, -1.0f, 1.0f));
  // phi_to_y<H> derives H_VIRT internally, mirroring pixel_to_vector's y_to_phi<H>.
  PixelCoords p({vector_to_theta<W>(v), phi_to_y<H>(phi)});
  return p;
}

/**
 * @brief Reflects a sample tap that ran past a pole back onto the sphere.
 * @tparam W Width (column count).
 * @tparam H Logical height (rows the buffer actually holds).
 * @tparam HOffset Virtual rows below the rendered domain.
 * @param col In/out column, in [0, W) on entry; shifted half a turn when the
 *   tap crosses a pole.
 * @param row In/out row; mirrored about the crossed pole.
 * @return False when nothing lies behind the tap: past the north pole by more
 *   than the buffer, or inside the virtual sub-pole gap when HOffset > 0.
 *   `col` and `row` are then unspecified.
 * @details The north pole sits exactly at row 0, the south pole at virtual row
 * H + HOffset - 1.
 */
template <int W, int H, int HOffset = hs::H_OFFSET>
HS_O3_FN bool pole_wrap(int &col, int &row) {
  if (row >= 0 && row < H)
    return true;
  constexpr int SOUTH = H + HOffset - 1;
  row = (row < 0) ? -row : 2 * SOUTH - row;
  if (row < 0 || row >= H)
    return false;
  col = fast_wrap(col + W / 2, W);
  return true;
}

} // namespace math
