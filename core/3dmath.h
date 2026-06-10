/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>
#include <cfloat>
#include <cstdint>
#include <cstring>
#include <limits>
#include <algorithm>
#include "platform.h"
#include "FastNoiseLite.h"

/**
 * @brief The Golden Ratio constant (Phi).
 */
static constexpr float PHI = 1.6180339887498948482045868343656f;
/**
 * @brief The inverse of the Golden Ratio (1/Phi).
 */
static constexpr float G = 1 / PHI;
/**
 * @brief Named tolerances for floating-point geometry comparisons.
 *
 * Use these named constants instead of ad-hoc literals so that
 * tolerance decisions are reviewed in one place. The naming reflects
 * what is being compared, not the magnitude:
 *
 *   EPS_UNIT       — unit-vector length test (1e-3f)
 *   TOLERANCE      — generic float compare (1e-4f)
 *   EPS_GEOMETRIC  — geometric near-equality (positions, angles) (1e-5f)
 *   EPS_LEN_SQ     — degenerate squared edge length (1e-6f)
 *   EPS_CROSS_SQ   — degenerate cross-product magnitude (squared) (1e-8f)
 *   EPS_NORMAL_SQ  — degenerate face normal (squared) (1e-9f)
 */
namespace math {
static constexpr float EPS_UNIT       = 1e-3f;
static constexpr float TOLERANCE      = 1e-4f;
static constexpr float EPS_GEOMETRIC  = 1e-5f;
static constexpr float EPS_LEN_SQ     = 1e-6f;
static constexpr float EPS_CROSS_SQ   = 1e-8f;
static constexpr float EPS_NORMAL_SQ  = 1e-9f;
// Cosine above which a vector is treated as parallel to a reference axis: too
// aligned to cross with it safely, so callers pick an alternate axis to avoid a
// near-zero (degenerate) cross when building a frame. Switches only when truly
// near-parallel (~0.8°) — the minimal-necessary guard.
static constexpr float COS_AXIS_PARALLEL = 1.0f - TOLERANCE;
}
// Back-compat alias — prefer math::TOLERANCE in new code.
static constexpr float TOLERANCE = math::TOLERANCE;
/**
 * @brief Floating-point representation of PI.
 */
static constexpr float PI_F = static_cast<float>(PI);

/**
 * @brief Quintic interpolation kernel (smootherstep).
 * @param t Interpolation factor (clamped to 0.0 - 1.0).
 * @return The interpolated value: t^3 * (t * (t * 6 - 15) + 10).
 */
inline float quintic_kernel(float t) {
  t = hs::clamp(t, 0.0f, 1.0f);
  return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
}

struct Vector;

/**
 * @brief Represents a point in Spherical Coordinates.
 */
struct Spherical {
  /**
   * @brief Constructs a Spherical coordinate.
   * @param theta Azimuthal angle.
   * @param phi Polar angle.
   */
  Spherical(float theta, float phi) : theta(theta), phi(phi) {}
  /**
   * @brief Constructs a Spherical coordinate from a 3D Vector.
   * @param v The Cartesian vector.
   */
  Spherical(const Vector &v);

  float theta; /**< Azimuthal angle (usually longitude). */
  float phi;   /**< Polar angle (usually co-latitude). */
};

/**
 * @brief Represents a 3D Vector in Cartesian coordinates (x, y, z).
 */
struct Vector {
  /**
   * @brief Default constructor (initializes to zero vector).
   */
  constexpr Vector() {}
  explicit constexpr Vector(int)
      : x(0), y(0), z(0) {} // inplace_function compat
  /**
   * @brief Constructs a vector with explicit components.
   * @param x X-component.
   * @param y Y-component.
   * @param z Z-component.
   */
  constexpr Vector(float x, float y, float z) : x(x), y(y), z(z) {}
  /**
   * @brief Copy constructor. Defaulted (trivial) so Vector is trivially
   * copyable — required for the memcpy fast paths (ArenaVector::append_bulk,
   * MeshState bulk copies) to be well-defined, and lets the optimizer
   * vectorize bulk copies.
   */
  constexpr Vector(const Vector &v) = default;
  /**
   * @brief Constructs a vector from Spherical coordinates.
   * @param s The spherical coordinate.
   */
  Vector(const Spherical &s)
      : x(sinf(s.phi) * cosf(s.theta)), y(cosf(s.phi)),
        z(sinf(s.phi) * sinf(s.theta)) {}
  /**
   * @brief Builds a unit vector from spherical angles.
   *
   * Named (not a 2-arg constructor) on purpose: `Vector(a, b)` reads as a
   * Cartesian point with z=0, but the spherical reading normalizes onto the
   * unit sphere. A 2-arg constructor made those two meanings indistinguishable
   * at the call site and silently picked the spherical one — dropping the third
   * Cartesian argument compiled into a different point. The factory forces the
   * intent to be spelled out.
   * @param theta Azimuthal angle.
   * @param phi Polar angle.
   * @return Unit vector at (theta, phi).
   */
  static Vector from_spherical(float theta, float phi) {
    return Vector(Spherical(theta, phi));
  }

  /**
   * @brief Copy assignment. Defaulted (trivial) — see the copy constructor;
   * keeps Vector trivially copyable.
   */
  constexpr Vector &operator=(const Vector &v) = default;

  /**
   * @brief Equality comparison using tolerance.
   * @param v The vector to compare.
   * @return True if components are within TOLERANCE.
   */
  bool operator==(const Vector &v) const {
    return std::abs(x - v.x) <= TOLERANCE && std::abs(y - v.y) <= TOLERANCE &&
           std::abs(z - v.z) <= TOLERANCE;
  }

  /**
   * @brief Inequality comparison.
   * @param v The vector to compare.
   * @return True if components are not equal within TOLERANCE.
   */
  bool operator!=(const Vector &v) const { return !(*this == v); }

  /**
   * @brief Unary negation operator.
   * @return A new vector with negated components.
   */
  constexpr Vector operator-() const { return Vector(-x, -y, -z); }

  Vector &operator+=(const Vector &v_other) {
    x += v_other.x;
    y += v_other.y;
    z += v_other.z;
    return *this;
  }
  Vector &operator-=(const Vector &v_other) {
    x -= v_other.x;
    y -= v_other.y;
    z -= v_other.z;
    return *this;
  }
  Vector &operator*=(float s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }
  Vector &operator/=(float s) {
    HS_CHECK(s != 0.0f);
    x /= s;
    y /= s;
    z /= s;
    return *this;
  }

  /**
   * @brief Calculates the magnitude (length) of the vector.
   * @return The magnitude.
   */
  [[nodiscard]] float length() const { return sqrtf(x * x + y * y + z * z); }

  /**
   * @brief Alias for length().
   */
  [[nodiscard]] float magnitude() const { return length(); }

  /**
   * @brief Normalizes the vector in place (scales to unit length).
   */
  void normalize() {
    float m = length();
    HS_CHECK(m >= std::numeric_limits<float>::epsilon());
    x = x / m;
    y = y / m;
    z = z / m;
  }

  /// Return a unit-length copy without mutating `this`.
  /// Traps on a zero-length vector — a degenerate input is a logic bug. Sites
  /// where a zero vector is a legitimate geometric edge (antipodal/coincident/
  /// pole/singularity) must use normalized_or() with an explicit fallback.
  [[nodiscard]] Vector normalized() const {
    float m = length();
    HS_CHECK(m >= std::numeric_limits<float>::epsilon());
    return Vector(x / m, y / m, z / m);
  }

  float x = 0; /**< X-component. */
  float y = 0; /**< Y-component. */
  float z = 0; /**< Z-component. */
};

/**
 * @brief Normalize, returning `fallback` when the vector is degenerate (zero
 * length) instead of trapping.
 *
 * The strict Vector::normalized() traps on a zero-length input so that an
 * unexpected zero surfaces as a bug. Use this variant at the handful of sites
 * where a zero vector is a *legitimate* geometric edge that arises during
 * normal animation — antipodal endpoints, a Hopf-fiber pole, a Möbius
 * singularity — so the degeneracy degrades gracefully to a stable direction.
 */
[[nodiscard]] inline Vector normalized_or(const Vector &v,
                                          const Vector &fallback) {
  float m = v.length();
  if (m < std::numeric_limits<float>::epsilon()) {
    return fallback;
  }
  return Vector(v.x / m, v.y / m, v.z / m);
}

// Fast atan2: 0.273-Hastings polynomial. Peak abs error ~0.0038 rad (~0.22°)
// measured over a dense sweep; worst near r ~= 0.7 in each octant.
inline float fast_atan2(float y, float x) {
  // +1e-10f keeps abs_y strictly positive so the (x==0,y==0) origin yields a
  // finite 0/abs_y ratio instead of a 0/0 NaN; it is far below the ~0.0038 rad
  // Hastings error above, so it never perturbs a non-degenerate result.
  float abs_y = std::abs(y) + 1e-10f;
  float abs_x = std::abs(x);
  float r, angle;

  if (abs_x < abs_y) {
    r = abs_x / abs_y;
    angle = 1.57079633f - r * (0.78539816f + 0.273f * (1.0f - r));
  } else {
    r = abs_y / abs_x;
    angle = r * (0.78539816f + 0.273f * (1.0f - r));
  }

  if (x < 0.0f)
    angle = 3.14159265f - angle;
  if (y < 0.0f)
    angle = -angle;

  return angle;
}

// Fast cube root for x >= 0. Bit-hack initial guess (divide the float exponent
// by three) refined by ONE Halley step (cubic convergence). Peak relative error
// ~2.3e-5 measured against cbrtf over [0,8] — vs ~1e-3 for a single Newton step
// — for ~2 extra multiplies and the same single division. Avoids the
// ~100-200-cycle soft-float cbrtf in the per-pixel OKLab hue path (see
// hue_rotate). Domain is x >= 0 (cbrt(0)=0); negative inputs return 0.
inline float fast_cbrt(float x) {
  if (x <= 0.0f)
    return 0.0f;
  uint32_t i;
  std::memcpy(&i, &x, sizeof(i));
  i = i / 3u + 0x2a514067u;
  float y;
  std::memcpy(&y, &i, sizeof(y));
  // Halley: y *= (y^3 + 2x) / (2y^3 + x)
  float c = y * y * y;
  return y * (c + 2.0f * x) / (2.0f * c + x);
}

// Forward declarations (defined at end of file)
inline float fast_acos(float x);
inline float fast_sinf(float x);
inline float fast_cosf(float x);

inline Spherical::Spherical(const Vector &v) {
  Vector n(v);
  n.normalize();
  theta = fast_atan2(n.z, n.x);
  phi = fast_acos(hs::clamp(n.y, -1.0f, 1.0f));
}

struct Quaternion;
constexpr Quaternion operator*(const Quaternion &q1, const Quaternion &q2);
constexpr Quaternion operator/(const Quaternion &q, float s);
/**
 * @brief Scalar division decl.
 * @param v Vector.
 * @param s Scalar.
 * @return The resulting vector.
 */
constexpr Vector operator/(const Vector &v, float s);

/**
 * @brief Represents a Quaternion for 3D rotation (r + i*x + j*y + k*z).
 */
struct Quaternion {
  /**
   * @brief Default constructor (identity quaternion).
   */
  constexpr Quaternion() {}
  /**
   * @brief Copy constructor.
   * @param q The quaternion to copy.
   */
  constexpr Quaternion(const Quaternion &q) : r(q.r), v(q.v) {}
  /**
   * @brief Constructs a quaternion with explicit components.
   * @param a Real component (r).
   * @param b X imaginary component (i).
   * @param c Y imaginary component (j).
   * @param d Z imaginary component (k).
   */
  constexpr Quaternion(float a, float b, float c, float d) : r(a), v(b, c, d) {}
  /**
   * @brief Constructs a quaternion from a scalar and a Vector part.
   * @param a Real component (r).
   * @param v Vector part (x, y, z).
   */
  constexpr Quaternion(float a, const Vector &v) : r(a), v(v) {}

  /**
   * @brief Assignment operator.
   * @param q The quaternion to assign from.
   * @return Reference to the current object.
   */
  Quaternion &operator=(const Quaternion &q) {
    r = q.r;
    v = q.v;
    return *this;
  }

  /**
   * @brief Equality comparison using tolerance.
   * @param q The quaternion to compare.
   * @return True if components are within TOLERANCE.
   */
  bool operator==(const Quaternion &q) const {
    return std::abs(q.r - r) < TOLERANCE && q.v == v;
  }

  /**
   * @brief Quaternion multiplication compound assignment.
   * @param q The quaternion to multiply by.
   */
  Quaternion &operator*=(const Quaternion &q) {
    *this = *this * q;
    return *this;
  }
  Quaternion &operator+=(const Quaternion &q) {
    r += q.r;
    v += q.v;
    return *this;
  }
  Quaternion &operator-=(const Quaternion &q) {
    r -= q.r;
    v -= q.v;
    return *this;
  }
  Quaternion &operator*=(float s) {
    r *= s;
    v *= s;
    return *this;
  }

  /**
   * @brief Calculates the squared magnitude of the quaternion.
   * @return The squared magnitude.
   */
  [[nodiscard]] constexpr float squared_magnitude() const {
    return r * r + v.x * v.x + v.y * v.y + v.z * v.z;
  }

  /**
   * @brief Calculates the inverse of the quaternion.
   * @return The inverse quaternion.
   */
  [[nodiscard]] Quaternion inverse() const {
    float sq_mag = squared_magnitude();
    // A zero-magnitude quaternion has no inverse; the divide below would
    // silently yield NaN/Inf. Trap instead — inverse() is computed per-frame
    // and cached (never per-pixel), so this guard is off the hot path.
    HS_CHECK(sq_mag > std::numeric_limits<float>::epsilon());
    return Quaternion(r, -v) / sq_mag;
  }

  /**
   * @brief Fast inverse for unit quaternions (avoids division).
   * Equivalent to conjugate() when |q| == 1.
   */
  [[nodiscard]] Quaternion unit_inverse() const {
    HS_CHECK(std::abs(squared_magnitude() - 1.0f) < 0.01f);
    return conjugate();
  }

  /**
   * @brief Calculates the conjugate of the quaternion.
   * @return The conjugate quaternion (r, -v).
   */
  [[nodiscard]] Quaternion conjugate() const { return Quaternion(r, -v); }

  /**
   * @brief Unary negation operator.
   * @return A new quaternion with negated components.
   */
  constexpr Quaternion operator-() const { return Quaternion(-r, -v); }

  /**
   * @brief Calculates the magnitude (length) of the quaternion.
   * @return The magnitude.
   */
  [[nodiscard]] float magnitude() const {
    return sqrtf(r * r + v.x * v.x + v.y * v.y + v.z * v.z);
  }

  /**
   * @brief Normalizes the quaternion in place (scales to unit magnitude).
   */
  void normalize() {
    auto m = magnitude();
    HS_CHECK(m > std::numeric_limits<float>::epsilon());
    r = r / m;
    v = v / m;
  }

  /// Return a unit-magnitude copy without mutating `this`.
  /// Traps on a zero-magnitude quaternion — a degenerate input is a logic bug.
  [[nodiscard]] Quaternion normalized() const {
    float m = magnitude();
    HS_CHECK(m > std::numeric_limits<float>::epsilon());
    return Quaternion(r / m, v / m);
  }

  float r = 1; /**< Real component. */
  Vector v;    /**< Vector part (x, y, z). */
};

/// Conventional representation of the point at infinity on the complex plane.
/// Single source of truth for the pole sentinel: every forward projection that
/// hits a singularity (stereo() at the north pole, gnomonic() at the equator,
/// Complex::operator/ on a near-zero denominator) emits a value of this
/// magnitude, and every inverse projection recognizes it (see the two
/// thresholds below).
static constexpr float STEREO_INF = 1e4f;

/// |z| at/above which inv_stereo() treats its input as the infinity sentinel.
/// Half of STEREO_INF on purpose: stereo() and Complex::operator/ emit exactly
/// STEREO_INF, but an intervening Mobius map can scale that toward (but not
/// past) zero, so the inverse needs margin below the emitted magnitude to still
/// snap a Mobius-shrunk sentinel back to the pole. (Squared in the comparison
/// to avoid a sqrt on the hot path.)
static constexpr float STEREO_INF_RECOGNIZE = STEREO_INF * 0.5f;

/// 1 - v.y below which stereo() is inside the north-pole cap and emits the
/// sentinel magnitude instead of the raw (and numerically explosive) quotient.
static constexpr float STEREO_POLE_EPS = 1e-4f;

/**
 * @brief Represents a Complex number.
 */
struct Complex {
  float re; /**< Real component. */
  float im; /**< Imaginary component. */

  /**
   * @brief Default constructor.
   */
  constexpr Complex() : re(0), im(0) {}
  /**
   * @brief Constructs a complex number.
   * @param r Real part.
   * @param i Imaginary part.
   */
  constexpr Complex(float r, float i) : re(r), im(i) {}

  /**
   * @brief Complex addition.
   */
  constexpr Complex operator+(const Complex &b) const {
    return Complex(re + b.re, im + b.im);
  }
  /**
   * @brief Complex subtraction.
   */
  constexpr Complex operator-(const Complex &b) const {
    return Complex(re - b.re, im - b.im);
  }
  /**
   * @brief Complex multiplication.
   */
  constexpr Complex operator*(const Complex &b) const {
    return Complex(re * b.re - im * b.im, re * b.im + im * b.re);
  }

  /**
   * @brief Complex division.
   */
  Complex operator/(const Complex &b) const {
    float denom = b.re * b.re + b.im * b.im;
    if (std::abs(denom) < 1e-6f) {
      // Near-zero denominator → point at infinity.
      // Preserve direction of the numerator.
      float num_mag = re * re + im * im;
      if (num_mag < 1e-12f)
        return Complex(0, 0); // 0/0 → indeterminate, keep zero
      float scale = STEREO_INF / sqrtf(num_mag);
      return Complex(re * scale, im * scale);
    }
    return Complex((re * b.re + im * b.im) / denom,
                   (im * b.re - re * b.im) / denom);
  }
};

inline Vector fold_to_hemisphere(const Vector &v) {
  // Half-angle identities to calculate the new height (j') and the new
  // horizontal radius (r')
  float j_new = sqrtf((1.0f + v.y) * 0.5f);
  float r_new = sqrtf((1.0f - v.y) * 0.5f);

  // Get the length of the original horizontal components (x and z)
  float r_old = sqrtf(v.x * v.x + v.z * v.z);

  // Protect against division by zero at the exact poles
  if (r_old < 1e-7f) {
    // The North Pole (j=1) stays the North Pole.
    // The South Pole (j=-1) is stretched into the entire equator ring.
    // We arbitrarily anchor the singularity to the X-axis (i=1).
    return v.y > 0.0f ? Vector(0.0f, 1.0f, 0.0f) : Vector(1.0f, 0.0f, 0.0f);
  }

  // Scale the original x and z coordinates to the new tightened radius
  float scale = r_new / r_old;

  return Vector(v.x * scale, j_new, v.z * scale);
}

/**
 * @brief Coefficients of a Mobius transform f(z) = (az + b) / (cz + d).
 *
 * Stores the four coefficients as first-class `Complex` values; animators mutate
 * the `.re`/`.im` components in place. The eight-float constructor is retained
 * for terse literal initialization.
 */
struct MobiusParams {
  Complex a, b, c, d;

  constexpr MobiusParams() : a(1, 0), b(0, 0), c(0, 0), d(1, 0) {}
  constexpr MobiusParams(Complex a_, Complex b_, Complex c_, Complex d_)
      : a(a_), b(b_), c(c_), d(d_) {}
  constexpr MobiusParams(float ar, float ai, float br, float bi, float cr,
                         float ci, float dr, float di)
      : a(ar, ai), b(br, bi), c(cr, ci), d(dr, di) {}
};

/**
 * @brief Stereographic Projection: Sphere -> Complex Plane.
 * North pole (v.y ≈ 1) maps to the point at infinity on the real axis.
 */
inline Complex stereo(const Vector &v) {
  float denom = 1.0f - v.y;
  if (denom < STEREO_POLE_EPS) {
    // Inside the north-pole cap. Emit the infinity sentinel, but keep the
    // (x,z) azimuth so a Mobius map that pulls the pole to a finite point
    // preserves the cap's swirl instead of collapsing it onto the +real axis.
    // At the exact pole (x = z = 0) the azimuth is undefined → +real fallback,
    // matching test_stereo_roundtrip's stereo(0,1,0).re ≈ STEREO_INF.
    float r = sqrtf(v.x * v.x + v.z * v.z);
    if (r < 1e-12f)
      return Complex(STEREO_INF, 0.0f);
    float scale = STEREO_INF / r;
    return Complex(v.x * scale, v.z * scale);
  }
  return Complex(v.x / denom, v.z / denom);
}

/**
 * @brief Inverse Stereographic Projection: Complex Plane -> Sphere.
 */
inline Vector inv_stereo(const Complex &z) {
  // Recognize the infinity sentinel (|z| >= STEREO_INF_RECOGNIZE) → exact North
  // Pole. Compared squared to keep the sqrt off this path.
  float r2 = z.re * z.re + z.im * z.im;
  if (r2 >= STEREO_INF_RECOGNIZE * STEREO_INF_RECOGNIZE)
    return Vector(0.0f, 1.0f, 0.0f);
  return Vector(2 * z.re / (r2 + 1), (r2 - 1) / (r2 + 1), 2 * z.im / (r2 + 1));
}

/**
 * @brief Mobius Transformation: f(z) = (az + b) / (cz + d).
 */
inline Complex mobius(const Complex &z, const MobiusParams &params) {
  Complex num = (params.a * z) + params.b;
  Complex den = (params.c * z) + params.d;
  return num / den;
}

/**
 * @brief Gnomonic Projection: Sphere -> Plane (Equator at Infinity).
 * Projects from center (0,0,0) to the plane y=1 (tangent at the North Pole
 * (0,1,0), i.e. j=1).
 */
inline Complex gnomonic(const Vector &v) {
  // Equator handling is two-stage and round-trips consistently with
  // inv_gnomonic. (1) Floor the divisor to ±1e-9 purely to avoid div-by-zero
  // at v.y == 0. (2) Clamp the result to ±STEREO_INF — this is the actual
  // equator→pole threshold. On the unit sphere a near-equator point (small
  // |v.y|) has |v.x| or |v.z| ≈ O(1), so at least one coordinate exceeds
  // STEREO_INF and is clamped to the same sentinel that inv_gnomonic
  // (|z| >= STEREO_INF) maps back to the pole. With STEREO_INF = 1e4 the flip
  // to "pole" happens around |v.y| < ~7e-5, far above the 1e-9 div floor — so
  // the floor never decides the round-trip; the matched STEREO_INF threshold
  // does.
  float div = (std::abs(v.y) < 1e-9f) ? 1e-9f * (v.y >= 0 ? 1.0f : -1.0f) : v.y;
  float gx = v.x / div;
  float gz = v.z / div;
  gx = hs::clamp(gx, -STEREO_INF, STEREO_INF);
  gz = hs::clamp(gz, -STEREO_INF, STEREO_INF);
  return Complex(gx, gz);
}

/**
 * @brief Inverse Gnomonic: Plane -> Sphere.
 * @param z Complex point on the plane.
 * @param original_sign Sign of the y-component (j) of the original vector, used
 * to restore the hemisphere the forward projection collapsed.
 */
inline Vector inv_gnomonic(const Complex &z, float original_sign = 1.0f) {
  // Recognize clamped-to-infinity → return the pole. gnomonic() clamps each
  // component to exactly ±STEREO_INF, so the inverse matches at exactly
  // STEREO_INF — no margin (unlike inv_stereo's STEREO_INF_RECOGNIZE), because
  // the forward sentinel here is an exact per-axis clamp, not a magnitude that a
  // Mobius map can shrink.
  if (std::abs(z.re) >= STEREO_INF || std::abs(z.im) >= STEREO_INF)
    return Vector(0.0f, original_sign, 0.0f);
  // Project (re, 1, im) back onto unit sphere
  float len = sqrtf(z.re * z.re + z.im * z.im + 1.0f);
  float inv_len = 1.0f / len;

  // Restore hemisphere sign (Upper or Lower)
  return Vector(z.re * inv_len * original_sign, // i
                inv_len * original_sign,        // j
                z.im * inv_len * original_sign  // k
  );
}

/**
 * @brief Rotates a vector by a unit quaternion.
 * @details Expanded q*v*conj(q) formula: 15 muls + 15 adds, no division.
 * @param v The vector to rotate.
 * @param q The unit rotation quaternion.
 * @return The rotated vector.
 */
inline Vector rotate(const Vector &v, const Quaternion &q) {
  float qr = q.r, qx = q.v.x, qy = q.v.y, qz = q.v.z;
  float tx = 2.0f * (qy * v.z - qz * v.y);
  float ty = 2.0f * (qz * v.x - qx * v.z);
  float tz = 2.0f * (qx * v.y - qy * v.x);
  return Vector(v.x + qr * tx + (qy * tz - qz * ty),
                v.y + qr * ty + (qz * tx - qx * tz),
                v.z + qr * tz + (qx * ty - qy * tx));
}

/**
 * @brief Vector addition.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The resulting vector.
 */
constexpr Vector operator+(const Vector &v1, const Vector &v2) {
  return Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

/**
 * @brief Vector subtraction.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The resulting vector.
 */
constexpr Vector operator-(const Vector &v1, const Vector &v2) {
  return Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

/**
 * @brief Scalar multiplication (vector * scalar).
 * @param v Vector.
 * @param s Scalar.
 * @return The resulting vector.
 */
constexpr Vector operator*(const Vector &v, float s) {
  return Vector(s * v.x, s * v.y, s * v.z);
}

/**
 * @brief Scalar multiplication (scalar * vector).
 * @param s Scalar.
 * @param v Vector.
 * @return The resulting vector.
 */
constexpr Vector operator*(float s, const Vector &v) { return v * s; }

/**
 * @brief Scalar division.
 * @param v Vector.
 * @param s Scalar.
 * @return The resulting vector.
 *
 * Unlike Vector::operator/=, this free operator deliberately does NOT trap on
 * s == 0: it is a per-pixel hot-path primitive, and HS_CHECK is reserved for
 * cold paths (see core/platform.h). A zero divisor yields ±Inf/NaN. Use the
 * compound-assignment form (v /= s) where a guarded divide is wanted.
 */
constexpr Vector operator/(const Vector &v, float s) {
  return Vector(v.x / s, v.y / s, v.z / s);
}

/**
 * @brief Calculates the dot product of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The dot product (scalar).
 */
constexpr float dot(const Vector &v1, const Vector &v2) {
  return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
}

/**
 * @brief Calculates the cross product of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The cross product vector.
 */
constexpr Vector cross(const Vector &v1, const Vector &v2) {
  return Vector(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
                v1.x * v2.y - v1.y * v2.x);
}

/**
 * @brief Calculates the Euclidean distance between two vectors.
 * @param a First vector.
 * @param b Second vector.
 * @return The distance (scalar).
 */
inline float distance_between(const Vector &a, const Vector &b) {
  float dx = b.x - a.x;
  float dy = b.y - a.y;
  float dz = b.z - a.z;
  return sqrtf(dx * dx + dy * dy + dz * dz);
}

/**
 * @brief Calculates the square of the Euclidean distance between two vectors.
 * Useful for efficient comparisons without square root.
 * @param a First vector.
 * @param b Second vector.
 * @return The distance squared (scalar).
 */
constexpr float distance_squared(const Vector &a, const Vector &b) {
  float dx = b.x - a.x;
  float dy = b.y - a.y;
  float dz = b.z - a.z;
  return dx * dx + dy * dy + dz * dz;
}

/**
 * @brief Calculates the angle between two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The angle in radians.
 */
inline float angle_between(const Vector &v1, const Vector &v2) {
  float len_product = v1.length() * v2.length();
  HS_CHECK(len_product > std::numeric_limits<float>::epsilon());
  float d = dot(v1, v2) / len_product;
  return fast_acos(hs::clamp(d, -1.0f, 1.0f));
}

/**
 * @brief Quaternion addition.
 * @param q1 First quaternion.
 * @param q2 Second quaternion.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator+(const Quaternion &q1, const Quaternion &q2) {
  return Quaternion(q1.r + q2.r, q1.v + q2.v);
}

/**
 * @brief Quaternion subtraction.
 * @param q1 First quaternion.
 * @param q2 Second quaternion.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator-(const Quaternion &q1, const Quaternion &q2) {
  return Quaternion(q1.r - q2.r, q1.v - q2.v);
}

/**
 * @brief Quaternion multiplication (Hamilton product).
 * @param q1 First quaternion.
 * @param q2 Second quaternion.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator*(const Quaternion &q1, const Quaternion &q2) {
  return Quaternion(q1.r * q2.r - dot(q1.v, q2.v),
                    q1.r * q2.v + q2.r * q1.v + cross(q1.v, q2.v));
}

/**
 * @brief Scalar multiplication (quaternion * scalar).
 * @param q Quaternion.
 * @param s Scalar.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator*(const Quaternion &q, float s) {
  return Quaternion(q.r * s, q.v * s);
}

/**
 * @brief Scalar multiplication (scalar * quaternion).
 * @param s Scalar.
 * @param q Quaternion.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator*(float s, const Quaternion &q) { return q * s; }

/**
 * @brief Scalar division.
 * @param q Quaternion.
 * @param s Scalar.
 * @return The resulting quaternion.
 *
 * Like operator/(Vector, float), this deliberately does NOT trap on s == 0 (hot
 * path; HS_CHECK is cold-path only). Callers that divide by a magnitude which
 * could be zero must guard themselves — e.g. Quaternion::inverse() HS_CHECKs the
 * squared magnitude before dividing.
 */
constexpr Quaternion operator/(const Quaternion &q, float s) {
  return Quaternion(q.r / s, q.v / s);
}

/**
 * @brief Calculates the dot product of two quaternions.
 * @param q1 First quaternion.
 * @param q2 Second quaternion.
 * @return The dot product (scalar).
 */
constexpr float dot(const Quaternion &q1, const Quaternion &q2) {
  return (q1.r * q2.r) + (q1.v.x * q2.v.x) + (q1.v.y * q2.v.y) +
         (q1.v.z * q2.v.z);
}

/**
 * @brief Calculates the angle between two unit quaternions.
 * @param q1 First unit quaternion.
 * @param q2 Second unit quaternion.
 * @return The angle in radians.
 */
inline float angle_between(const Quaternion &q1, const Quaternion &q2) {
  // Exact acosf (unlike the Vector overload's fast_acos): quaternion variants are
  // orientation-level ops invoked at most per-object per-frame, never per-pixel,
  // so the cycles are immaterial while drift-free rotation fidelity is not.
  return acosf(hs::clamp(dot(q1, q2), -1.0f, 1.0f));
}

/**
 * @brief Creates a rotation quaternion from an axis and an angle.
 * @param axis The rotation axis (must be unit vector).
 * @param theta The angle of rotation in radians.
 * @return The resulting unit rotation quaternion.
 */
inline Quaternion make_rotation(const Vector &axis, float theta) {
  return Quaternion(fast_cosf(theta / 2), fast_sinf(theta / 2) * axis).normalized();
}

/**
 * @brief Creates a rotation quaternion to rotate from one vector to another.
 * @param from The source vector (must be unit vector).
 * @param to The destination vector (must be unit vector).
 * @return The resulting unit rotation quaternion.
 */
inline Quaternion make_rotation(const Vector &from, const Vector &to) {
  // Inputs must be unit: the d>1-eps / d<-1+eps branch tests below assume |from|
  // = |to| = 1, and a non-unit input silently skews the angle. Trap it, matching
  // make_basis's unit-quaternion guard. (dot(v,v) avoids the sqrt; the 0.02
  // squared-magnitude band is the ~1% magnitude tolerance make_basis uses.)
  //
  // The 0.02 here and the TOLERANCE (1e-4) below are deliberately on different
  // scales because they measure different things: 0.02 is *input-validity*
  // slack — how far the caller's vectors may drift from unit length and still
  // be accepted — and is generous on purpose so accumulated normalization error
  // never false-traps. TOLERANCE is a *geometric branch threshold* — how close
  // d = dot(from,to) must be to ±1 to take the (anti)parallel special case (and
  // how degenerate a cross product must be to swap reference axes) — and is
  // kept tight so the general path's cross-product normalize stays well-
  // conditioned. Reconciling them to one scale would be a category error.
  HS_CHECK(std::abs(dot(from, from) - 1.0f) < 0.02f &&
               std::abs(dot(to, to) - 1.0f) < 0.02f,
           "make_rotation(from, to): inputs must be unit vectors");
  float d = dot(from, to);

  // Handle antiparallel vectors (180 degrees apart)
  if (d < -1.0f + TOLERANCE) {
    // Choose an axis perpendicular to 'from'
    Vector axis = cross(Vector(1, 0, 0), from);
    if (axis.length() < TOLERANCE) {
      axis = cross(Vector(0, 1, 0), from);
    }
    axis.normalize();
    return make_rotation(axis, PI_F);
  }

  // Handle parallel vectors (0 degrees apart)
  if (d > 1.0f - TOLERANCE) {
    return Quaternion(1, 0, 0, 0); // Identity
  }

  auto axis = cross(from, to).normalized();
  auto angle = angle_between(from, to);
  return make_rotation(axis, angle);
}

// Abramowitz & Stegun polynomial approximation for acos(x), x in [-1, 1].
// Peak abs error ~1.3e-4 rad (~0.0072°) measured over a dense sweep.
inline float fast_acos(float x) {
  float ax = std::abs(x);
  if (ax > 1.0f)
    ax = 1.0f;
  float result =
      sqrtf(1.0f - ax) *
      (1.5707963f + ax * (-0.2121144f + ax * (0.0742610f + ax * -0.0187293f)));
  return x < 0.0f ? PI_F - result : result;
}

// Bhaskara I sine approximation. Max error ~0.17% (~0.1 degrees).
// ~8 FPU ops vs ~80-120 cycles for newlib sinf.
inline float fast_sinf(float x) {
  // Range reduce to [0, 2π)
  constexpr float INV_2PI = 1.0f / (2.0f * PI_F);
  x = x - floorf(x * INV_2PI) * (2.0f * PI_F);
  // Map to [0, π) with sign tracking
  float sign = 1.0f;
  if (x > PI_F) {
    x -= PI_F;
    sign = -1.0f;
  }
  // Bhaskara I: 16x(π-x) / (5π² - 4x(π-x))
  float xpi = x * (PI_F - x);
  return sign * (16.0f * xpi) / (5.0f * PI_F * PI_F - 4.0f * xpi);
}

inline float fast_cosf(float x) { return fast_sinf(x + PI_F * 0.5f); }

/**
 * @brief Spherical Linear Interpolation (SLERP) between two Vectors.
 * @param v1 Starting vector.
 * @param v2 Ending vector.
 * @param t Interpolation factor (0.0 to 1.0).
 * @return The interpolated unit vector.
 */
inline Vector slerp(const Vector &v1, const Vector &v2, float t) {
  float d = hs::clamp(dot(v1, v2), -1.0f, 1.0f);
  // If vectors are extremely close, just lerp to avoid NaN
  if (d > 0.9999f) {
    return (v1 + (v2 - v1) * t).normalized();
  }
  float theta = fast_acos(d);
  float sin_theta = fast_sinf(theta);
  if (sin_theta < 0.0001f) {
    // Antipodal endpoints: the lerp midpoint collapses to the zero vector for an
    // exact antipode at t≈0.5, where the great-circle path is undefined anyway.
    // Fall back to a stable endpoint direction rather than trapping in strict
    // normalized() — an antipodal pair is a legitimate geometric edge here.
    return normalized_or((v1 + (v2 - v1) * t), v1);
  }
  float s1 = fast_sinf((1 - t) * theta) / sin_theta;
  float s2 = fast_sinf(t * theta) / sin_theta;
  return ((s1 * v1) + (s2 * v2)).normalized();
}

/**
 * @brief Spherical Linear Interpolation (SLERP) between two quaternions.
 * @param q1 Starting quaternion.
 * @param q2 Ending quaternion.
 * @param t Interpolation factor (0.0 to 1.0).
 * @param long_way If true, takes the longest path between quaternions.
 * @return The interpolated unit quaternion.
 */
inline Quaternion slerp(const Quaternion &q1, const Quaternion &q2, float t,
                        bool long_way = false) {
  float d = dot(q1, q2);
  Quaternion p(q1);
  Quaternion q(q2);

  if ((long_way && d > 0) || (!long_way && d < 0)) {
    p = -p;
    d = -d;
  }

  if (d > (1 - TOLERANCE)) {
    Quaternion r = p + t * (q - p);
    return r.normalized();
  }

  // Exact acosf/sinf (unlike the Vector slerp's fast_*): orientation interpolation
  // runs at most per-object per-frame, not per-pixel, so exact trig costs nothing
  // measurable here and keeps rotation paths free of approximation drift.
  float theta = acosf(hs::clamp(d, -1.0f, 1.0f));
  float sin_theta = sinf(theta);
  if (sin_theta < 0.0001f) {
    Quaternion r = p + t * (q - p);
    return r.normalized();
  }
  float s1 = sinf((1 - t) * theta) / sin_theta;
  float s2 = sinf(t * theta) / sin_theta;
  return ((s1 * p) + (s2 * q)).normalized();
}

/**
 * @brief Evaluation mode for spherical spline interpolation.
 */
enum class SplineMode { Fast, Geodesic };

namespace Spline {

/// Fast: polynomial interpolation + normalize. O(1) sqrtf per sample.
/// Accurate within ~30° control-point separation; distorts on long arcs.
inline Vector cubic_fast(const Vector &p0, const Vector &p1, const Vector &p2,
                         const Vector &p3, float t) {
  float u = 1.0f - t;
  float uu = u * u;
  float tt = t * t;
  Vector blended = p0 * (uu * u) + p1 * (3.0f * uu * t) + p2 * (3.0f * u * tt) +
                   p3 * (tt * t);
  // Antipodal / coincident control points can sum to a near-zero vector —
  // legitimate geometry, not a bug — and the strict normalized() would trap.
  // Degrade gracefully to p1, matching the robustness cubic_slerp already gets
  // from slerp's normalized_or fallback.
  return normalized_or(blended, p1);
}

/// Accurate: de Casteljau with SLERP. 6 SLERPs per sample.
/// Shape-faithful at any arc length.
inline Vector cubic_slerp(const Vector &p0, const Vector &p1, const Vector &p2,
                          const Vector &p3, float t) {
  Vector b01 = slerp(p0, p1, t);
  Vector b12 = slerp(p1, p2, t);
  Vector b23 = slerp(p2, p3, t);
  Vector b012 = slerp(b01, b12, t);
  Vector b123 = slerp(b12, b23, t);
  return slerp(b012, b123, t);
}

/// Unified dispatch — selects evaluation strategy once per sample point.
inline Vector cubic(const Vector &p0, const Vector &p1, const Vector &p2,
                    const Vector &p3, float t, SplineMode mode) {
  return (mode == SplineMode::Fast) ? cubic_fast(p0, p1, p2, p3, t)
                                    : cubic_slerp(p0, p1, p2, p3, t);
}

/// Catmull-Rom tangent estimation on the sphere.
/// Returns two Bézier control points for the segment from `start` to `end`.
/// @note tension=0 produces geodesic (linear) segments; tension=1 produces
/// full Catmull-Rom smoothing. This is inverted from the conventional
/// Catmull-Rom tension parameter where τ=0 is the standard cardinal spline.
inline void catmull_rom_tangents(const Vector &prev, const Vector &start,
                                 const Vector &end, const Vector &next,
                                 float tension, Vector &cp1, Vector &cp2) {
  // Tangent at start: pull toward the midpoint of prev→end
  cp1 = slerp(start, slerp(prev, end, 0.5f), tension);
  // Tangent at end: pull toward the midpoint of start→next
  cp2 = slerp(end, slerp(start, next, 0.5f), tension);
}

} // namespace Spline

