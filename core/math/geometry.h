/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file geometry.h
 * @brief Sphere/canvas mapping, the axis constants, the spherical geometry
 *        helpers, and `Orientation<CAP>`, the per-frame quaternion history the
 *        animation and filter layers interpolate over.
 */

#include <cmath>
#include <algorithm>
#include <array>
#include <utility>
#include "engine/platform.h"
#include "math/3dmath.h"
#include "engine/util.h" // for fast_wrap()

/**
 * @brief Unit vector along the Cartesian X-axis.
 */
inline constexpr Vector X_AXIS(1, 0, 0);
/**
 * @brief Unit vector along the Cartesian Y-axis.
 */
inline constexpr Vector Y_AXIS(0, 1, 0);
/**
 * @brief Unit vector along the Cartesian Z-axis.
 */
inline constexpr Vector Z_AXIS(0, 0, 1);
/**
 * @brief Unit vector along the Cartesian Y-axis.
 */
inline constexpr Vector UP = Y_AXIS;

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
  HS_CHECK(y >= 0 && y < PhiLUT<H>::H_VIRT);
  return PhiLUT<H>::data[y];
}

/**
 * @brief Pixel-y -> phi for fractional rows at compile-time height H.
 * @tparam H Logical height; H_VIRT is H + hs::H_OFFSET.
 * @param y Fractional pixel row.
 * @return The spherical phi angle in radians.
 * @details Snaps to the LUT for near-integer y; otherwise computes the angle
 * analytically. The range is debug-asserted, so only NDEBUG builds extrapolate
 * linearly past [0, H_VIRT-1]; keeping y in range is the caller's
 * responsibility.
 */
template <int H> inline float y_to_phi(float y) {
  if (std::abs(y - std::floor(y)) < TOLERANCE) {
    int iy = static_cast<int>(y);
    if (iy >= 0 && iy < PhiLUT<H>::H_VIRT) {
      return y_to_phi<H>(iy);
    }
  }
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
 * guards in the rasterizers then remain only as a lazy fallback for unit tests
 * and offline tools; their non-atomic check-then-set relies on this eager call
 * and the single-render-thread assumption. Idempotent.
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
  if (std::abs(x - std::floor(x)) < TOLERANCE &&
      std::abs(y - std::floor(y)) < TOLERANCE) {
    const int ix = static_cast<int>(x);
    const int iy = static_cast<int>(y);
    if (ix >= 0 && ix < W && iy >= 0 && iy < TrigLUT<W, H>::H_VIRT) {
      return pixel_to_vector<W, H>(ix, iy);
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

/**
 * @brief Calculates a point on the Fibonacci spiral on the unit sphere.
 * @param n The total number of points in the spiral (must be positive).
 * @param eps The epsilon offset for the spiral.
 * @param i The index of the point to calculate.
 * @return The point on the unit sphere.
 * @note Setup-time generator; exact trig is intentional (not a per-pixel path).
 */
inline Vector fib_spiral(int n, float eps, int i) {
  HS_CHECK(n > 0, "fib_spiral: n must be positive");
  // Clamp before acosf: float rounding can push the argument past -1 → NaN.
  float phi = acosf(hs::clamp(1.0f - (2.0f * (static_cast<float>(i) + eps)) /
                                         static_cast<float>(n),
                              -1.0f, 1.0f));
  float theta =
      fmodf((2.0f * PI_F * static_cast<float>(i) * INV_PHI), (2.0f * PI_F));
  // Y-up; unit by construction, so no normalize().
  return Vector(sinf(phi) * cosf(theta), cosf(phi), sinf(phi) * sinf(theta));
}

/**
 * @brief Class managing the current rotation state of an object, maintaining
 * history for interpolation.
 * @tparam CAP Maximum number of orientation frames retained in history.
 * @details Stores a list of Quaternions (`orientations`) generated during the
 * current frame step.
 */
template <int CAP = 4> class Orientation {
public:
  static constexpr int CAPACITY = CAP;
  /**
   * @brief Default constructor (identity rotation).
   */
  Orientation() : num_frames(0) { set(Quaternion()); }

  /**
   * @brief Constructs with a specific initial quaternion.
   * @param q The initial quaternion.
   */
  explicit Orientation(const Quaternion &q) : num_frames(0) { set(q); }

  /**
   * @brief Gets the number of recorded orientation frames in the current step.
   * @return The number of frames.
   */
  int length() const { return num_frames; }

  /**
   * @brief Rotates a vector by the current (latest) orientation.
   * @param v The vector to orient.
   * @return The rotated vector.
   */
  Vector orient(const Vector &v) const {
    HS_CHECK(num_frames >= 1);
    return rotate(v, orientations[num_frames - 1]);
  }

  /**
   * @brief Rotates a vector by an orientation at a specific historical frame
   * index.
   * @param v The vector to orient.
   * @param i The frame index (0 being oldest, length-1 being current).
   * @return The rotated vector.
   */
  Vector orient(const Vector &v, int i) const {
    HS_CHECK(i >= 0 && i < num_frames);
    return rotate(v, orientations[i]);
  }

  /**
   * @brief Rotates a vector backward by the inverse of the current (latest)
   * orientation.
   * @param v The vector to unorient.
   * @return The unrotated vector.
   */
  Vector unorient(const Vector &v) const {
    HS_CHECK(num_frames >= 1);
    return rotate(v, orientations[num_frames - 1].conjugate());
  }

  /**
   * @brief Rotates a vector backward by the inverse of a specific historical
   * orientation.
   * @param v The vector to unorient.
   * @param i The frame index.
   * @return The unrotated vector.
   */
  Vector unorient(const Vector &v, int i) const {
    HS_CHECK(i >= 0 && i < num_frames);
    return rotate(v, orientations[i].conjugate());
  }

  /**
   * @brief Gets the current (latest) quaternion.
   * @return The Quaternion reference.
   */
  const Quaternion &get() const {
    HS_CHECK(num_frames >= 1);
    return orientations[num_frames - 1];
  }

  /**
   * @brief Gets the quaternion at a specific historical frame index.
   * @param i The frame index.
   * @return The Quaternion reference.
   */
  const Quaternion &get(int i) const {
    HS_CHECK(i >= 0 && i < num_frames);
    return orientations[i];
  }

  /**
   * @brief Sets the orientation, clearing all history.
   * @param q The new orientation quaternion; MUST be unit length
   *   (HS_CHECK-trapped — a non-unit quaternion scales every rotated vector).
   * @return Reference to the Orientation object.
   */
  Orientation &set(const Quaternion &q) {
    HS_CHECK(std::abs(q.squared_magnitude() - 1.0f) < math::EPS_UNIT_QUAT_SQ);
    orientations[0] = q;
    num_frames = 1;
    return *this;
  }

  /**
   * @brief Pushes a new quaternion onto the history, tracking a motion step.
   * @param q The new rotation quaternion; MUST be unit length
   *   (HS_CHECK-trapped — a non-unit quaternion scales every rotated vector).
   * @return Reference to the Orientation object.
   */
  Orientation &push(const Quaternion &q) {
    HS_CHECK(std::abs(q.squared_magnitude() - 1.0f) < math::EPS_UNIT_QUAT_SQ);
    HS_CHECK(num_frames < CAPACITY);
    orientations[num_frames++] = q;
    return *this;
  }

  /**
   * @brief Collapses the orientation history, retaining only the latest
   * quaternion.
   * @details Used after rendering a motion step to reset the motion blur
   * history.
   * @return Reference to the Orientation object.
   */
  Orientation &collapse() {
    if (num_frames > 1) {
      orientations[0] = orientations[num_frames - 1];
      num_frames = 1;
    }
    return *this;
  }

  /**
   * @brief Access a quaternion at a specific historical frame index.
   * @param i The frame index.
   * @return The const Quaternion reference.
   */
  const Quaternion &at(int i) const {
    HS_CHECK(i >= 0 && i < num_frames);
    return orientations[i];
  }

  /**
   * @brief Replaces a quaternion at a specific historical frame index.
   * @param i The frame index.
   * @param q Unit-length replacement quaternion.
   */
  void set_at(int i, const Quaternion &q) {
    HS_CHECK(i >= 0 && i < num_frames);
    HS_CHECK(std::abs(q.squared_magnitude() - 1.0f) < math::EPS_UNIT_QUAT_SQ);
    orientations[i] = q;
  }

  /**
   * @brief Normalizes and replaces a historical quaternion.
   * @param i The frame index.
   * @param q Replacement quaternion; must have nonzero magnitude.
   */
  void set_at_normalized(int i, Quaternion q) {
    HS_CHECK(i >= 0 && i < num_frames);
    q.normalize();
    orientations[i] = q;
  }

  /**
   * @brief Increases the resolution of the history to 'count' steps by
   * Slerp-resampling the existing frames (uniform in source index, not arc
   * length).
   * @param count The target number of steps in the history.
   * @note When `count > CAPACITY` the trail is upsampled to `CAPACITY` instead:
   * graceful degradation (the write stays in-bounds, the current orientation
   * stays exact, only the motion-blur smear samples more coarsely), not a trap.
   * Raise `CAP` to trade RAM for a smoother fast smear.
   * @note A single-frame source (`num_frames == 1`, the common post-`set()`/
   * `collapse()` state) upsamples to a flat smear of that one frame; real motion
   * blur requires >=2 pushed frames.
   * @note Resampling is uniform in source index, so a trail whose frames were
   * pushed at a varying rate keeps that non-uniformity: the motion-blur streak
   * stays dense over the slow stretches and sparse over the fast ones rather
   * than spreading evenly along the arc.
   */
  void upsample(int count) {
    HS_CHECK(count >= 1);
    if (count > CAPACITY)
      count = CAPACITY;
    if (num_frames >= count)
      return;

    std::array<Quaternion, CAPACITY> old_orientations;
    std::copy(orientations.begin(), orientations.begin() + num_frames,
              old_orientations.begin());

    int old_num_frames = num_frames;

    for (int i = 0; i < count - 1; ++i) {
      float t = static_cast<float>(i) / (count - 1);
      float source_float_index = t * (old_num_frames - 1);
      int idx = static_cast<int>(source_float_index);
      float frac = source_float_index - idx;

      orientations[i] = slerp(
          old_orientations[idx],
          old_orientations[std::min((int)old_num_frames - 1, idx + 1)], frac);
    }
    // Endpoint maps exactly onto the last source frame (slerp(q, q, 0)).
    orientations[count - 1] = old_orientations[old_num_frames - 1].normalized();
    num_frames = count;
  }

private:
  std::array<Quaternion, CAPACITY>
      orientations; /**< Storage for historical quaternions. */
  int num_frames;   /**< The current number of active frames in history. */
};

/**
 * @brief Generates a truly random 3D unit vector (direction) using Marsaglia's
 * method.
 * @return A normalized random Vector.
 */
inline Vector random_vector() {
  float v1, v2, s;
  // Marsaglia rejection: accept when (v1,v2) lands in the open unit disk.
  do {
    v1 = 2.0f * hs::rand_f() - 1.0f;
    v2 = 2.0f * hs::rand_f() - 1.0f;
    s = v1 * v1 + v2 * v2;
  } while (s >= 1.0f || s == 0.0f);

  float sqrt_s = sqrtf(1.0f - s);
  return Vector(2.0f * v1 * sqrt_s, 2.0f * v2 * sqrt_s, 1.0f - 2.0f * s);
}

/**
 * @brief Parameters defining a Lissajous curve.
 */
struct LissajousParams {
  float m1; /**< Frequency coefficient for the axial components (X and Z). */
  float m2; /**< Frequency coefficient for the orbital component (Y). */
  float a; /**< Phase shift in radians (matches the daydream lissajous tool). */
  float domain; /**< The total duration (t) over which the curve is drawn. */
};

/**
 * @brief Calculates a 3D point on the unit sphere corresponding to a spherical
 * Lissajous curve.
 * @param m1 Frequency coefficient for XZ plane.
 * @param m2 Frequency coefficient for Y axis.
 * @param a Phase shift in radians, used as-is (matches the daydream lissajous
 *          designer tools/lissajous.html).
 * @param t Time variable (or position along the domain).
 * @return The calculated 3D point (unit vector).
 * @note Setup-time generator; exact trig is intentional (not a per-pixel path).
 */
inline Vector lissajous(float m1, float m2, float a, float t) {
  // Unit by construction, so no normalize().
  return Vector(sinf(m2 * t) * cosf(m1 * t - a), cosf(m2 * t),
                sinf(m2 * t) * sinf(m1 * t - a));
}

/**
 * @brief An orthonormal basis { u, v, w } whose *middle* axis is the normal:
 * `v` is the normal, and `u`, `w` span the plane perpendicular to it.
 */
struct Basis {
  Vector u, v, w;
};

/**
 * @brief Rotates a basis by a unit quaternion; the axes stay orthonormal.
 * @param b Basis to rotate.
 * @param q Unit rotation quaternion.
 * @return The basis with every axis rotated by @p q.
 */
inline Basis rotate(const Basis &b, const Quaternion &q) {
  return {rotate(b.u, q), rotate(b.v, q), rotate(b.w, q)};
}

/**
 * @brief Creates a basis { u, v, w } from an orientation and normal.
 * @param orientation The orientation quaternion; MUST be unit length
 *   (HS_CHECK-trapped below — a non-unit quaternion would scale/shear the frame).
 * @param normal The normal vector; after rotation it becomes the 'v' axis. MUST
 *   be non-zero — a zero (or rotation-collapsed) normal traps in the
 *   `normalized()` of `v` below.
 * @return The constructed Basis.
 */
inline Basis make_basis(const Quaternion &orientation, const Vector &normal) {
  HS_CHECK(std::abs(orientation.squared_magnitude() - 1.0f) <
           math::EPS_UNIT_QUAT_SQ);
  Vector v = rotate(normal, orientation).normalized();
  // rotate preserves dot, so least_parallel_axis(normal) picks the same body
  // axis as the rotated frame; rotate it into the frame for the cross. Only its
  // direction matters, the cross below is normalized.
  Vector ref = rotate(least_parallel_axis(normal), orientation);
  Vector u = cross(v, ref).normalized();
  // v and u are orthonormal, so the cross is unit by construction.
  Vector w = cross(v, u);
  return {u, v, w};
}

/**
 * @brief Adjusted basis and radius for drawing on the opposite side of the
 * sphere.
 * @param basis The current basis {u, v, w}.
 * @param radius Angular radius (0-2).
 * @return A pair containing the adjusted Basis and radius.
 */
inline std::pair<Basis, float> get_antipode(const Basis &basis, float radius) {
  if (radius > 1.0f) {
    Basis new_basis;
    new_basis.u = -basis.u; // Flip U to maintain chirality
    new_basis.v = -basis.v; // Flip V (Antipode)
    new_basis.w = basis.w;  // W stays (Rotation axis)
    return {new_basis, 2.0f - radius};
  }
  return {basis, radius};
}
