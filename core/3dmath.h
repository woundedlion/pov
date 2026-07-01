/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>
#include <cfloat>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <algorithm>
// platform.h first: on device it defines NDEBUG, which must be set before
// <cassert> expands the assert macro — otherwise assert-stripping would depend
// on a prior TU having pulled in platform.h, making this header non-self-sufficient.
#include "platform.h"
#include <cassert>
#include "FastNoiseLite.h"

/**
 * @brief The Golden Ratio constant (Phi).
 */
static constexpr float PHI = 1.6180339887498948482045868343656f;
/**
 * @brief The inverse of the Golden Ratio (1/Phi).
 */
static constexpr float INV_PHI = 1 / PHI;
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
 *   EPS_NORMALIZE_SQ — squared length below which normalize() has no reliable
 *                      direction; |v| < 1e-6 amplifies by >1e6 (1e-12f)
 *   EPS_UNIT_QUAT_SQ — generous |q|^2 is-unit assertion slack (0.01f)
 *   EPS_UNIT_VEC_SQ  — generous |v|^2 is-unit assertion slack (0.02f)
 */
namespace math {
static constexpr float EPS_UNIT       = 1e-3f;
static constexpr float TOLERANCE      = 1e-4f;
static constexpr float EPS_GEOMETRIC  = 1e-5f;
static constexpr float EPS_LEN_SQ     = 1e-6f;
static constexpr float EPS_CROSS_SQ   = 1e-8f;
static constexpr float EPS_NORMAL_SQ  = 1e-9f;
static constexpr float EPS_NORMALIZE_SQ = 1e-12f;
static constexpr float EPS_UNIT_QUAT_SQ = 0.01f;
static constexpr float EPS_UNIT_VEC_SQ  = 0.02f;
/**
 * @brief Cosine above which a vector is treated as parallel to a reference axis.
 * @details Too aligned to cross with it safely, so callers pick an alternate
 * axis to avoid a near-zero (degenerate) cross when building a frame. Switches
 * only when truly near-parallel (~0.8°) — the minimal-necessary guard.
 */
static constexpr float COS_AXIS_PARALLEL = 1.0f - TOLERANCE;
}
/**
 * @brief Global alias for math::TOLERANCE; prefer the namespaced name in new
 * code.
 */
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
   * @param v The Cartesian vector (must be non-degenerate; normalize() traps
   *          on a zero/near-zero-length input below normalized()'s epsilon).
   * @details explicit: a coordinate-space change must be spelled out
   * (`Spherical{v}`) so a Cartesian vector cannot silently bind where a
   * Spherical is expected, and vice versa.
   */
  explicit Spherical(const Vector &v);

  float theta; /**< Azimuthal angle from +x toward +z (atan2 convention; the
                    (theta, phi) ctor stores any value unchecked). */
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
  /**
   * @brief Constructs a zero vector from a null-pointer tag.
   * @details Exists only for the device `teensy::inplace_function`, whose empty
   * state synthesizes a return value as `Vector(0)`. Taking `std::nullptr_t`
   * still accepts that literal `0` (a null-pointer constant) while rejecting the
   * `Vector(5)` footgun the old `Vector(int)` allowed — a literal int that
   * silently collapsed to the zero vector, discarding its value.
   */
  explicit constexpr Vector(std::nullptr_t)
      : x(0), y(0), z(0) {} // inplace_function empty-state compat
  /**
   * @brief Constructs a vector with explicit components.
   * @param x X-component.
   * @param y Y-component.
   * @param z Z-component.
   */
  constexpr Vector(float x, float y, float z) : x(x), y(y), z(z) {}
  /**
   * @brief Copy constructor (defaulted, trivial).
   * @param v The vector to copy.
   * @details Defaulted so Vector stays trivially copyable — required for the
   * memcpy fast paths (ArenaVector::append_bulk, MeshState bulk copies) to be
   * well-defined, and lets the optimizer vectorize bulk copies.
   */
  constexpr Vector(const Vector &v) = default;
  /**
   * @brief Constructs a vector from Spherical coordinates.
   * @param s The spherical coordinate.
   * @details explicit: callers use the named from_spherical() factory for the
   * intended spherical->unit-vector path, so block any implicit conversion that
   * would let a Spherical silently bind where a Cartesian Vector is expected.
   */
  explicit Vector(const Spherical &s)
      : x(sinf(s.phi) * cosf(s.theta)), y(cosf(s.phi)),
        z(sinf(s.phi) * sinf(s.theta)) {}
  /**
   * @brief Builds a unit vector from spherical angles.
   *
   * A named factory rather than a 2-arg constructor on purpose: `Vector(a, b)`
   * would be ambiguous with a Cartesian point (z=0), and the spherical reading
   * normalizes onto the unit sphere — two distinct meanings the call site must
   * not confuse. The factory forces the spherical intent to be spelled out.
   * @param theta Azimuthal angle.
   * @param phi Polar angle.
   * @return Unit vector at (theta, phi).
   */
  static Vector from_spherical(float theta, float phi) {
    return Vector(Spherical(theta, phi));
  }

  /**
   * @brief Copy assignment (defaulted, trivial).
   * @param v The vector to copy from.
   * @return Reference to this vector.
   * @details Defaulted so Vector stays trivially copyable, keeping the memcpy
   * fast paths well-defined.
   */
  constexpr Vector &operator=(const Vector &v) = default;

  /**
   * @brief Equality comparison using tolerance.
   * @param v The vector to compare.
   * @return True if components are within TOLERANCE.
   * @note Uses an absolute per-component tolerance (TOLERANCE), so it is
   *   scale-dependent (loose for large magnitudes, strict for tiny ones) and
   *   non-transitive (a==b and b==c does not imply a==c). Fine for direct
   *   value compares; do not rely on it for ordering or container keys.
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

  /**
   * @brief Adds another vector component-wise into this vector.
   * @param v_other The vector to add.
   * @return Reference to this vector.
   */
  Vector &operator+=(const Vector &v_other) {
    x += v_other.x;
    y += v_other.y;
    z += v_other.z;
    return *this;
  }
  /**
   * @brief Subtracts another vector component-wise from this vector.
   * @param v_other The vector to subtract.
   * @return Reference to this vector.
   */
  Vector &operator-=(const Vector &v_other) {
    x -= v_other.x;
    y -= v_other.y;
    z -= v_other.z;
    return *this;
  }
  /**
   * @brief Scales this vector by a scalar in place.
   * @param s Scalar multiplier.
   * @return Reference to this vector.
   */
  Vector &operator*=(float s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }
  /**
   * @brief Divides this vector by a scalar in place.
   * @param s Scalar divisor (must be non-zero; traps otherwise).
   * @return Reference to this vector.
   */
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
   * @return The magnitude (length) of the vector.
   */
  [[nodiscard]] float magnitude() const { return length(); }

  /**
   * @brief Normalizes the vector in place (scales to unit length).
   */
  void normalize() {
    float m2 = x * x + y * y + z * z;
    HS_CHECK(m2 >= math::EPS_NORMALIZE_SQ);
    float m = sqrtf(m2);
    x = x / m;
    y = y / m;
    z = z / m;
  }

  /**
   * @brief Returns a unit-length copy without mutating `this`.
   * @return A unit-length vector in the same direction.
   * @details Traps on a zero-length vector — a degenerate input is a logic bug.
   * Sites where a zero vector is a legitimate geometric edge (antipodal/
   * coincident/pole/singularity) must use normalized_or() with an explicit
   * fallback.
   */
  [[nodiscard]] Vector normalized() const {
    float m2 = x * x + y * y + z * z;
    HS_CHECK(m2 >= math::EPS_NORMALIZE_SQ);
    float m = sqrtf(m2);
    return Vector(x / m, y / m, z / m);
  }

  float x = 0; /**< X-component. */
  float y = 0; /**< Y-component. */
  float z = 0; /**< Z-component. */
};

/**
 * @brief Normalizes, returning `fallback` when the vector is degenerate (zero
 * length) instead of trapping.
 * @param v The vector to normalize.
 * @param fallback Direction returned when `v` has near-zero length.
 * @return A unit-length copy of `v`, or `fallback` if `v` is degenerate.
 * @details The strict Vector::normalized() traps on a zero-length input so that
 * an unexpected zero surfaces as a bug. Use this variant at the handful of sites
 * where a zero vector is a *legitimate* geometric edge that arises during normal
 * animation — antipodal endpoints, a Hopf-fiber pole, a Möbius singularity — so
 * the degeneracy degrades gracefully to a stable direction.
 */
[[nodiscard]] inline Vector normalized_or(const Vector &v,
                                          const Vector &fallback) {
  float m2 = v.x * v.x + v.y * v.y + v.z * v.z;
  if (m2 < math::EPS_NORMALIZE_SQ) {
    return fallback;
  }
  float m = sqrtf(m2);
  return Vector(v.x / m, v.y / m, v.z / m);
}

/**
 * @brief Fast atan2 using the 0.273-Hastings polynomial.
 * @param y Y (numerator) coordinate.
 * @param x X (denominator) coordinate.
 * @return The angle in radians in [-π, π], matching std::atan2; the −π
 * endpoint is attained for x<0 as y→0⁻.
 * @details Peak abs error ~0.0038 rad (~0.22°) measured over a dense sweep;
 * worst near r ~= 0.7 in each octant. The two-branch split at |x|==|y| is C0 at
 * the exact diagonal but the +1e-10f bias on abs_y makes the branches disagree
 * by ~6e-7 rad right at the octant boundaries — far below the polynomial error
 * and harmless for all callers (hue/azimuth lookups).
 */
inline float fast_atan2(float y, float x) {
  // +1e-10f keeps abs_y strictly positive so the (0,0) origin stays finite
  // (~PI/2 from the abs_x<abs_y branch) instead of NaN; callers tolerate an
  // arbitrary angle at the undefined origin.
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

/**
 * @brief Fast cube root for x >= 0.
 * @param x Input value; the domain is x >= 0 (cbrt(0)=0), negative inputs
 * return 0. Accuracy holds for x >~ 1e-6; in the denormal/tiny-normal tail
 * below that the bit-hack seed degrades and the result collapses toward 0.
 * @return An approximation of the cube root of `x`.
 * @details Bit-hack initial guess (divide the float exponent by three) refined
 * by ONE Halley step (cubic convergence). Peak relative error ~2.3e-5 measured
 * against cbrtf for x >= 1e-6 — vs ~1e-3 for a single Newton step — for ~2 extra
 * multiplies and the same single division; the bit-hack seed degrades in the
 * denormal/tiny-normal tail, far below the inputs this serves. Avoids the
 * ~100-200-cycle soft-float cbrtf in the per-pixel OKLab hue path (see
 * hue_rotate).
 */
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
  // At the poles (n.y == ±1) theta is arbitrary but harmless: phi==0 or π
  // collapses the azimuth so every theta maps to the same point.
  // theta is fast_atan2's [-π, π], not [0,2π); consumers needing the latter wrap.
  theta = fast_atan2(n.z, n.x);
  phi = fast_acos(hs::clamp(n.y, -1.0f, 1.0f));
}

struct Quaternion;
constexpr Quaternion operator*(const Quaternion &q1, const Quaternion &q2);
constexpr Quaternion operator/(const Quaternion &q, float s);
/**
 * @brief Scalar division of a vector (forward declaration).
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
  constexpr Quaternion(const Quaternion &q) = default;
  /**
   * @brief Constructs a quaternion with explicit components.
   * @param a Real component (r).
   * @param b X imaginary component (i).
   * @param c Y imaginary component (j).
   * @param d Z imaginary component (k).
   */
  constexpr explicit Quaternion(float a, float b, float c, float d)
      : r(a), v(b, c, d) {}
  /**
   * @brief Constructs a quaternion from a scalar and a Vector part.
   * @param a Real component (r).
   * @param v Vector part (x, y, z).
   */
  constexpr explicit Quaternion(float a, const Vector &v) : r(a), v(v) {}

  /**
   * @brief Assignment operator.
   * @param q The quaternion to assign from.
   * @return Reference to the current object.
   */
  Quaternion &operator=(const Quaternion &q) = default;

  /**
   * @brief Equality comparison using tolerance.
   * @param q The quaternion to compare.
   * @return True if components are within TOLERANCE.
   * @note Absolute per-component tolerance (TOLERANCE): scale-dependent and
   *   non-transitive (see Vector::operator==). Suitable for direct value
   *   compares, not for ordering or container keys. Component equality, not
   *   rotation equality: `q` and `-q` are the same rotation but compare unequal.
   */
  bool operator==(const Quaternion &q) const {
    return std::abs(q.r - r) <= TOLERANCE && q.v == v;
  }

  /**
   * @brief Quaternion multiplication compound assignment.
   * @param q The quaternion to multiply by.
   * @return Reference to this quaternion.
   */
  Quaternion &operator*=(const Quaternion &q) {
    *this = *this * q;
    return *this;
  }
  /**
   * @brief Component-wise addition compound assignment.
   * @param q The quaternion to add.
   * @return Reference to this quaternion.
   */
  Quaternion &operator+=(const Quaternion &q) {
    r += q.r;
    v += q.v;
    return *this;
  }
  /**
   * @brief Component-wise subtraction compound assignment.
   * @param q The quaternion to subtract.
   * @return Reference to this quaternion.
   */
  Quaternion &operator-=(const Quaternion &q) {
    r -= q.r;
    v -= q.v;
    return *this;
  }
  /**
   * @brief Scales this quaternion by a scalar in place.
   * @param s Scalar multiplier.
   * @return Reference to this quaternion.
   */
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
    HS_CHECK(sq_mag > std::numeric_limits<float>::epsilon());
    return conjugate() / sq_mag;
  }

  /**
   * @brief Fast inverse for unit quaternions (avoids division).
   * @return The inverse quaternion, equivalent to conjugate() when |q| == 1.
   */
  [[nodiscard]] Quaternion unit_inverse() const {
    HS_CHECK(std::abs(squared_magnitude() - 1.0f) < math::EPS_UNIT_QUAT_SQ);
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

  /**
   * @brief Returns a unit-magnitude copy without mutating `this`.
   * @return A unit-magnitude quaternion.
   * @details Traps on a zero-magnitude quaternion — a degenerate input is a
   * logic bug.
   */
  [[nodiscard]] Quaternion normalized() const {
    float m = magnitude();
    HS_CHECK(m > std::numeric_limits<float>::epsilon());
    return Quaternion(r / m, v / m);
  }

  float r = 1; /**< Real component. */
  Vector v;    /**< Vector part (x, y, z). */
};

/**
 * @brief Conventional representation of the point at infinity on the complex
 * plane.
 * @details Single source of truth for the pole sentinel: every forward
 * projection that hits a singularity (stereo() at the north pole, gnomonic() at
 * the equator, Complex::operator/ on a near-zero denominator) emits a value of
 * this magnitude, and every inverse projection recognizes it (see the two
 * thresholds below).
 */
static constexpr float STEREO_INF = 1e4f;

/**
 * @brief |z| at/above which inv_stereo() treats its input as the infinity
 * sentinel.
 * @details Half of STEREO_INF on purpose: stereo() and Complex::operator/ emit
 * exactly STEREO_INF, but an intervening Mobius map can scale that toward (but
 * not past) zero, so the inverse needs margin below the emitted magnitude to
 * still snap a Mobius-shrunk sentinel back to the pole. (Squared in the
 * comparison to avoid a sqrt on the hot path.)
 */
static constexpr float STEREO_INF_RECOGNIZE = STEREO_INF * 0.5f;

/**
 * @brief Soft-limit for a stereographic coordinate fed (times a pattern
 * frequency) into fast_sinf/fast_cosf, beyond which range reduction bands.
 * @details stereo() emits magnitudes up to STEREO_INF (1e4); multiplied by a
 * shader's pattern frequency this drives trig arguments toward ~2e5, where a
 * float's ULP (~0.02 rad) makes fast_sinf's range reduction visibly shimmer/band
 * near the pole. Clamping the trig argument to ±this bound holds the reduction
 * error to ~5e-4 rad. The clamped band is the pole cap, which the stereographic
 * shaders pole-attenuate toward the background anyway, so freezing the pattern
 * there is invisible; non-pole coordinates (|z| ~ O(10)) never reach the bound.
 */
static constexpr float STEREO_PATTERN_ARG_LIMIT = 4096.0f;

/**
 * @brief 1 - v.y below which stereo() is inside the north-pole cap and emits the
 * sentinel magnitude instead of the raw (and numerically explosive) quotient.
 */
static constexpr float STEREO_POLE_EPS = 1e-4f;

/**
 * @brief Squared numerator magnitude below which Complex::operator/ treats a
 * near-zero-divided-by-near-zero quotient as the indeterminate 0/0 form and
 * returns zero, rather than scaling the numerator's direction by the infinity
 * sentinel. Tighter than math::EPS_LEN_SQ (the divisor-degeneracy threshold)
 * on purpose: only a numerator that is essentially exactly zero should lose its
 * direction.
 */
static constexpr float STEREO_DIV_NUM_EPS_SQ = 1e-12f;

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
   * @param b The right-hand operand.
   * @return The sum of this and `b`.
   */
  constexpr Complex operator+(const Complex &b) const {
    return Complex(re + b.re, im + b.im);
  }
  /**
   * @brief Complex subtraction.
   * @param b The right-hand operand.
   * @return The difference of this and `b`.
   */
  constexpr Complex operator-(const Complex &b) const {
    return Complex(re - b.re, im - b.im);
  }
  /**
   * @brief Complex multiplication.
   * @param b The right-hand operand.
   * @return The product of this and `b`.
   */
  constexpr Complex operator*(const Complex &b) const {
    return Complex(re * b.re - im * b.im, re * b.im + im * b.re);
  }

  /**
   * @brief Complex division.
   * @param b The divisor.
   * @return The quotient of this and `b`; a quotient whose magnitude would reach
   * STEREO_INF collapses to the infinity sentinel scaled along the numerator's
   * direction, so a merely small (non-singular) divisor still yields a finite
   * point.
   */
  Complex operator/(const Complex &b) const {
    float denom = b.re * b.re + b.im * b.im;
    float num_mag = re * re + im * im;
    if (num_mag >= denom * (STEREO_INF * STEREO_INF)) {
      if (num_mag < STEREO_DIV_NUM_EPS_SQ)
        return Complex(0, 0); // 0/0 → indeterminate, keep zero
      float scale = STEREO_INF / sqrtf(num_mag);
      return Complex(re * scale, im * scale);
    }
    return Complex((re * b.re + im * b.im) / denom,
                   (im * b.re - re * b.im) / denom);
  }
};

/**
 * @brief Coefficients of a Mobius transform f(z) = (az + b) / (cz + d).
 * @details Stores the four coefficients as first-class `Complex` values;
 * animators mutate the `.re`/`.im` components in place. The eight-float
 * constructor is retained for terse literal initialization.
 */
struct MobiusParams {
  Complex a, b, c, d; /**< The four transform coefficients. */

  /**
   * @brief Default constructor producing the identity transform (a=d=1, b=c=0).
   */
  constexpr MobiusParams() : a(1, 0), b(0, 0), c(0, 0), d(1, 0) {}
  /**
   * @brief Constructs from four Complex coefficients.
   * @param a_ Coefficient a.
   * @param b_ Coefficient b.
   * @param c_ Coefficient c.
   * @param d_ Coefficient d.
   */
  constexpr MobiusParams(Complex a_, Complex b_, Complex c_, Complex d_)
      : a(a_), b(b_), c(c_), d(d_) {}
  /**
   * @brief Constructs from eight floats (real/imaginary pairs per coefficient).
   * @param ar Real part of a.
   * @param ai Imaginary part of a.
   * @param br Real part of b.
   * @param bi Imaginary part of b.
   * @param cr Real part of c.
   * @param ci Imaginary part of c.
   * @param dr Real part of d.
   * @param di Imaginary part of d.
   */
  constexpr MobiusParams(float ar, float ai, float br, float bi, float cr,
                         float ci, float dr, float di)
      : a(ar, ai), b(br, bi), c(cr, ci), d(dr, di) {}
};

/**
 * @brief Stereographic Projection: Sphere -> Complex Plane.
 * @param v Point on the unit sphere.
 * @return The projected complex-plane coordinate.
 * @details North pole (v.y ≈ 1) maps to the point at infinity on the real axis.
 */
inline Complex stereo(const Vector &v) {
  float denom = 1.0f - v.y;
  if (denom < STEREO_POLE_EPS) {
    // North-pole cap: emit the infinity sentinel but keep the (x,z) azimuth so a
    // Mobius map pulling the pole finite preserves the swirl. At the exact pole
    // (x = z = 0) the azimuth is undefined → +real fallback.
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
 * @param z Complex-plane coordinate (the infinity sentinel maps to the pole).
 * @return The corresponding point on the unit sphere.
 */
inline Vector inv_stereo(const Complex &z) {
  // |z| >= STEREO_INF_RECOGNIZE is the infinity sentinel → North Pole (squared
  // compare avoids a sqrt).
  float r2 = z.re * z.re + z.im * z.im;
  if (r2 >= STEREO_INF_RECOGNIZE * STEREO_INF_RECOGNIZE)
    return Vector(0.0f, 1.0f, 0.0f);
  return Vector(2 * z.re / (r2 + 1), (r2 - 1) / (r2 + 1), 2 * z.im / (r2 + 1));
}

/**
 * @brief Mobius Transformation: f(z) = (az + b) / (cz + d).
 * @param z The complex input point.
 * @param params The four transform coefficients.
 * @return The transformed complex point.
 */
inline Complex mobius(const Complex &z, const MobiusParams &params) {
  Complex num = (params.a * z) + params.b;
  Complex den = (params.c * z) + params.d;
  return num / den;
}

/**
 * @brief Gnomonic Projection: Sphere -> Plane (Equator at Infinity).
 * @param v Point on the unit sphere.
 * @return The projected plane coordinate (equator points clamp to the sentinel).
 * @details Projects from center (0,0,0) to the plane y=1 (tangent at the North
 * Pole (0,1,0), i.e. j=1).
 * @note Identifies antipodes: `v` and `-v` project to the same plane
 * coordinate, so the hemisphere is lost. A caller that round-trips through
 * inv_gnomonic must track the sign of `v.y` and pass it back via
 * inv_gnomonic's `original_sign`.
 */
inline Complex gnomonic(const Vector &v) {
  // Equator handling: floor the divisor to ±1e-9 to avoid div-by-zero at
  // v.y == 0, then clamp to ±STEREO_INF. A near-equator point has |v.x| or
  // |v.z| ≈ O(1), so it clamps to the sentinel; inv_gnomonic snaps any magnitude
  // at/above STEREO_INF_RECOGNIZE (half the sentinel) back to the pole, so the
  // equator collapses to the pole just inside |v.y| ~ 1.4e-4.
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
 * @return The corresponding point on the unit sphere.
 */
inline Vector inv_gnomonic(const Complex &z, float original_sign = 1.0f) {
  // Clamped-to-infinity → pole, recognized from STEREO_INF_RECOGNIZE (half the
  // sentinel) like inv_stereo: a Mobius map between gnomonic() and here can
  // scale the exact ±STEREO_INF clamp toward (but not past) zero, so the margin
  // still snaps a shrunk sentinel back to the pole.
  if (std::abs(z.re) >= STEREO_INF_RECOGNIZE ||
      std::abs(z.im) >= STEREO_INF_RECOGNIZE)
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

// Hot-path primitive: the device build has no s == 0 trap (yields ±Inf/NaN);
// callers guard with HS_CHECK when s may be zero. The debug-only assert surfaces
// unguarded misuse in native/WASM-debug at zero device cost.
constexpr Vector operator/(const Vector &v, float s) {
  assert(s != 0.0f);
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
  float len_product = sqrtf(dot(v1, v1) * dot(v2, v2));
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
 * Like operator/(Vector, float), the device build deliberately does NOT trap on
 * s == 0 (hot path; HS_CHECK is cold-path only). Callers that divide by a
 * magnitude which could be zero must guard themselves — e.g.
 * Quaternion::inverse() HS_CHECKs the squared magnitude before dividing. The
 * debug-only assert surfaces unguarded misuse at zero device cost.
 */
constexpr Quaternion operator/(const Quaternion &q, float s) {
  assert(s != 0.0f);
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
 * @note Exact `acosf`, unlike the Vector overload's `fast_acos`.
 */
inline float angle_between(const Quaternion &q1, const Quaternion &q2) {
  return acosf(hs::clamp(dot(q1, q2), -1.0f, 1.0f));
}

/**
 * @brief Creates a rotation quaternion from an axis and an angle.
 * @param axis The rotation axis (must be a non-zero unit vector).
 * @param theta The angle of rotation in radians.
 * @return The resulting unit rotation quaternion.
 * @pre `axis` is unit length and non-zero. A zero/degenerate axis collapses the
 *      quaternion's vector part to zero; at `theta = pi` the scalar part
 *      `cos(pi/2)` is also ~0, so the whole quaternion's magnitude falls below
 *      `normalized()`'s epsilon threshold (not exact zero) and it traps. This
 *      is deliberate fail-fast on a
 *      precondition violation, not a soft fallback — callers that may produce a
 *      degenerate axis must guard it (e.g. via `normalized_or`) before calling.
 */
inline Quaternion make_rotation(const Vector &axis, float theta) {
  return Quaternion(cosf(theta / 2), sinf(theta / 2) * axis).normalized();
}

/**
 * @brief Returns the canonical body axis least parallel to @p v.
 * @param v Vector to find a well-conditioned reference axis for.
 * @return +X unless @p v is near-parallel to it, in which case +Y.
 * @details Seed for building a cross-product frame around @p v: crossing @p v
 *          with a near-parallel axis collapses to ~zero, so pick the axis that
 *          is safely off-axis. Centralizes the reference-axis pick that the
 *          path/walk frame builders and make_rotation's antiparallel branch
 *          would otherwise spell out by hand. (Lives here, not in geometry.h
 *          with X_AXIS/Y_AXIS, so make_rotation below can use it without an
 *          include cycle; dot(v, X_AXIS) is just v.x.)
 */
inline Vector least_parallel_axis(const Vector &v) {
  // Scale-invariant: compare v.x^2 against the threshold scaled by |v|^2, so a
  // non-unit v can't slip past the cosine test and seed a degenerate cross.
  constexpr float k = math::COS_AXIS_PARALLEL;
  return v.x * v.x > k * k * dot(v, v) ? Vector(0, 1, 0) : Vector(1, 0, 0);
}

/**
 * @brief Creates a rotation quaternion to rotate from one vector to another.
 * @param from The source vector (must be unit vector).
 * @param to The destination vector (must be unit vector).
 * @return The resulting unit rotation quaternion.
 * @details Near-parallel inputs snap to identity: when d = dot(from,to) exceeds
 * 1 - TOLERANCE the shortest-arc angle is below ~0.81° (1 - cos θ < 1e-4), so the
 * rotation is dropped rather than divided through a near-zero cross product. Sub-
 * 0.81° rotations are therefore silently lost — harmless for the orientation-level
 * callers here, but a caller needing exact tiny rotations must not route them
 * through this overload. The antiparallel case (d < -1 + TOLERANCE) instead
 * synthesizes a stable perpendicular axis and rotates by π.
 */
inline Quaternion make_rotation(const Vector &from, const Vector &to) {
  HS_CHECK(std::abs(dot(from, from) - 1.0f) < math::EPS_UNIT_VEC_SQ &&
               std::abs(dot(to, to) - 1.0f) < math::EPS_UNIT_VEC_SQ,
           "make_rotation(from, to): inputs must be unit vectors");
  float d = dot(from, to);

  // Antiparallel (180°): seed the cross with least_parallel_axis(from) so the
  // axis we normalize through stays well-conditioned.
  if (d < -1.0f + TOLERANCE) {
    Vector axis = cross(least_parallel_axis(from), from);
    axis.normalize();
    return make_rotation(axis, PI_F);
  }

  if (d > 1.0f - TOLERANCE) {
    return Quaternion(1, 0, 0, 0); // Identity
  }

  auto axis = cross(from, to).normalized();
  float angle = fast_acos(d);
  return make_rotation(axis, angle);
}

/**
 * @brief Builds the unit rotation quaternion whose columns are the given
 * orthonormal axes.
 * @param cx Image of the body +X axis (must be unit).
 * @param cy Image of the body +Y axis (must be unit, orthogonal to cx).
 * @param cz Image of the body +Z axis (must be unit, = cross(cx, cy)).
 * @return The unit quaternion q with rotate(X_AXIS, q) == cx, etc.
 * @details Shepperd's method: pick the largest of the four (1 ± trace terms) so
 * the divisor is never near zero. Unlike composing two shortest-arc rotations,
 * this is singular-free for every orientation — there is no antipodal axis to
 * guard. An orientation-level op (per object per frame, never per pixel), so the
 * one sqrt + branch is immaterial. Inputs are assumed orthonormal. A mildly
 * skewed basis silently yields a non-rotation (the trailing `normalized()` still
 * succeeds, but `q` no longer reproduces the intended frame), so callers must
 * hand it a true frame. A *near-singular* basis (collinear/degenerate axes)
 * drives the chosen Shepperd divisor — and hence `q` — toward zero magnitude, at
 * which point the final `normalized()` traps (fail-fast). The two regimes are not
 * a contradiction: the silent path is the skewed-but-nondegenerate one; the trap
 * is reserved for a basis too degenerate to normalize at all.
 */
inline Quaternion quaternion_from_basis(const Vector &cx, const Vector &cy,
                                        const Vector &cz) {
  // m[col] = image of that body axis; mIJ is row I, col J of the rotation matrix.
  const float m00 = cx.x, m10 = cx.y, m20 = cx.z;
  const float m01 = cy.x, m11 = cy.y, m21 = cy.z;
  const float m02 = cz.x, m12 = cz.y, m22 = cz.z;
  const float trace = m00 + m11 + m22;
  Quaternion q;
  if (trace > 0.0f) {
    float s = sqrtf(trace + 1.0f) * 2.0f; // s = 4 * q.r
    q.r = 0.25f * s;
    q.v.x = (m21 - m12) / s;
    q.v.y = (m02 - m20) / s;
    q.v.z = (m10 - m01) / s;
  } else if (m00 > m11 && m00 > m22) {
    float s = sqrtf(1.0f + m00 - m11 - m22) * 2.0f; // s = 4 * q.v.x
    q.r = (m21 - m12) / s;
    q.v.x = 0.25f * s;
    q.v.y = (m01 + m10) / s;
    q.v.z = (m02 + m20) / s;
  } else if (m11 > m22) {
    float s = sqrtf(1.0f + m11 - m00 - m22) * 2.0f; // s = 4 * q.v.y
    q.r = (m02 - m20) / s;
    q.v.x = (m01 + m10) / s;
    q.v.y = 0.25f * s;
    q.v.z = (m12 + m21) / s;
  } else {
    float s = sqrtf(1.0f + m22 - m00 - m11) * 2.0f; // s = 4 * q.v.z
    q.r = (m10 - m01) / s;
    q.v.x = (m02 + m20) / s;
    q.v.y = (m12 + m21) / s;
    q.v.z = 0.25f * s;
  }
  return q.normalized();
}

/**
 * @brief Fast acos using the Abramowitz & Stegun polynomial approximation.
 * @param x Input value, expected in [-1, 1] (clamped internally).
 * @return The arc cosine in radians.
 * @details Peak abs error ~1.3e-4 rad (~0.0072°) measured over a dense sweep.
 */
inline float fast_acos(float x) {
  float ax = std::abs(x);
  if (ax > 1.0f)
    ax = 1.0f;
  float result =
      sqrtf(1.0f - ax) *
      (1.5707963f + ax * (-0.2121144f + ax * (0.0742610f + ax * -0.0187293f)));
  return x < 0.0f ? PI_F - result : result;
}

/**
 * @brief Fast sine using the Bhaskara I approximation.
 * @param x Angle in radians (range-reduced internally).
 * @return The approximate sine of `x`.
 * @details Max error ~0.17% (~0.1 degrees); ~8 FPU ops vs ~80-120 cycles for
 * newlib sinf.
 * @warning Accuracy degrades for large |x|: the `x - floor(x/2π)·2π` range
 * reduction loses precision as a float's ULP grows past the period, so callers
 * driving large arguments must bound them first (see STEREO_PATTERN_ARG_LIMIT).
 */
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

/**
 * @brief Fast cosine, computed as fast_sinf(x + π/2).
 * @param x Angle in radians.
 * @return The approximate cosine of `x`.
 * @warning Accuracy degrades for large |x|: the `x - floor(x/2π)·2π` range
 * reduction loses precision as a float's ULP grows past the period, so callers
 * driving large arguments must bound them first (see STEREO_PATTERN_ARG_LIMIT).
 */
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
  // Near-parallel: lerp to avoid NaN. Branch thresholds are TOLERANCE-scale.
  if (d > 1.0f - math::TOLERANCE) {
    return (v1 + (v2 - v1) * t).normalized();
  }
  float theta = fast_acos(d);
  // cos(π - TOLERANCE) ≈ -1 + TOLERANCE²/2; test raw dot so the antipodal
  // decision does not depend on fast_acos's error envelope.
  if (d < -1.0f + math::TOLERANCE * math::TOLERANCE * 0.5f) {
    // Antipodal: v1 + v2 cancels so the standard blend is non-monotone
    // (collapses to ~v1 for all t). Rotate v1 by t*theta about an arbitrary
    // perpendicular axis — a monotone half-turn landing within TOLERANCE of v2.
    Vector axis = cross(least_parallel_axis(v1), v1).normalized();
    return rotate(v1, make_rotation(axis, t * theta));
  }
  // The 1/sinθ factor in the slerp weights cancels under the final normalize(),
  // so it is omitted (also dodging a divide by a small sin(θ) near θ→π).
  // Fast trig here (vs the exact sinf/acosf in the Quaternion overload): a
  // direction blend tolerates fast_sinf's bounded error; orientation slerp does not.
  float s1 = fast_sinf((1 - t) * theta);
  float s2 = fast_sinf(t * theta);
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

  // At d == 0 (orthogonal endpoints) the >=/< split sends each mode to the arc
  // its name promises: long_way takes the 270° arc, short-path the 90° arc.
  if ((long_way && d >= 0) || (!long_way && d < 0)) {
    p = -p;
    d = -d;
  }

  if (d > (1 - TOLERANCE)) {
    Quaternion r = p + t * (q - p);
    return r.normalized();
  }

  // Exact acosf/sinf (vs fast_* in the Vector overload): orientation interpolation
  // is more sensitive to angular error than a direction blend.
  float theta = acosf(hs::clamp(d, -1.0f, 1.0f));
  float sin_theta = sinf(theta);
  if (sin_theta < math::TOLERANCE) {
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

/**
 * @brief Fast cubic spline sample: polynomial interpolation + normalize.
 * @param p0 First control point.
 * @param p1 Second control point.
 * @param p2 Third control point.
 * @param p3 Fourth control point.
 * @param t Interpolation factor in [0, 1].
 * @return The interpolated unit vector (falls back to p1 if degenerate).
 * @details O(1) sqrtf per sample. Accurate within ~30° control-point
 * separation; distorts on long arcs.
 */
inline Vector cubic_fast(const Vector &p0, const Vector &p1, const Vector &p2,
                         const Vector &p3, float t) {
  float u = 1.0f - t;
  float uu = u * u;
  float tt = t * t;
  Vector blended = p0 * (uu * u) + p1 * (3.0f * uu * t) + p2 * (3.0f * u * tt) +
                   p3 * (tt * t);
  // Antipodal/coincident control points can sum to ~zero (a legitimate edge),
  // so degrade to p1 rather than trapping in strict normalized().
  return normalized_or(blended, p1);
}

/**
 * @brief Accurate cubic spline sample: de Casteljau with SLERP.
 * @param p0 First control point.
 * @param p1 Second control point.
 * @param p2 Third control point.
 * @param p3 Fourth control point.
 * @param t Interpolation factor in [0, 1].
 * @return The interpolated unit vector.
 * @details 6 SLERPs per sample. Shape-faithful at any arc length.
 */
inline Vector cubic_slerp(const Vector &p0, const Vector &p1, const Vector &p2,
                          const Vector &p3, float t) {
  Vector b01 = slerp(p0, p1, t);
  Vector b12 = slerp(p1, p2, t);
  Vector b23 = slerp(p2, p3, t);
  Vector b012 = slerp(b01, b12, t);
  Vector b123 = slerp(b12, b23, t);
  return slerp(b012, b123, t);
}

/**
 * @brief Unified cubic spline dispatch — selects the evaluation strategy once
 * per sample point.
 * @param p0 First control point.
 * @param p1 Second control point.
 * @param p2 Third control point.
 * @param p3 Fourth control point.
 * @param t Interpolation factor in [0, 1].
 * @param mode Evaluation strategy (Fast or Geodesic).
 * @return The interpolated unit vector.
 */
inline Vector cubic(const Vector &p0, const Vector &p1, const Vector &p2,
                    const Vector &p3, float t, SplineMode mode) {
  return (mode == SplineMode::Fast) ? cubic_fast(p0, p1, p2, p3, t)
                                    : cubic_slerp(p0, p1, p2, p3, t);
}

/**
 * @brief Catmull-Rom tangent estimation on the sphere.
 * @param prev Control point before the segment.
 * @param start Segment start point.
 * @param end Segment end point.
 * @param next Control point after the segment.
 * @param smoothing Smoothing amount in [0, 1] (see @details for convention).
 * @param cp1 Output: first Bézier control point for the start→end segment.
 * @param cp2 Output: second Bézier control point for the start→end segment.
 * @details smoothing=0 produces geodesic (linear) segments; smoothing=1
 * produces full Catmull-Rom smoothing. Named `smoothing` rather than the
 * conventional Catmull-Rom `tension` because the sense is inverted from that
 * convention (where τ=0 is the standard cardinal spline).
 */
inline void catmull_rom_tangents(const Vector &prev, const Vector &start,
                                 const Vector &end, const Vector &next,
                                 float smoothing, Vector &cp1, Vector &cp2) {
  // Tangent at start: pull toward the midpoint of prev→end
  cp1 = slerp(start, slerp(prev, end, 0.5f), smoothing);
  // Tangent at end: pull toward the midpoint of start→next
  cp2 = slerp(end, slerp(start, next, 0.5f), smoothing);
}

} // namespace Spline

