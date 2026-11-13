/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file rotation_interpolation.h @brief Trigonometric kernels and spherical rotation interpolation. */

#include "core/math/rotation_math.h"

namespace math {

/** @brief Absolute angular error bound in radians on [-1, 1]. */
inline constexpr double FAST_ACOS_MAX_ERROR = 1e-4;
/// Measured peak error used by narrow annular bounds.
inline constexpr float FAST_ACOS_PEAK_ERROR = 5.1e-5f;

/**
 * @brief Fast acos using the Abramowitz & Stegun polynomial approximation.
 * @param x Input value, expected in [-1, 1] (clamped internally).
 * @return The arc cosine in radians; zero for NaN.
 * @details Peak abs error ~5.0e-5 rad.
 */
__attribute__((always_inline)) inline float fast_acos(float x) {
  const float ax = fminf(std::abs(x), 1.0f);
  float result =
      sqrtf(1.0f - ax) *
      (1.5707963f + ax * (-0.2132209f + ax * (0.0777570f + ax * -0.0215124f)));
  return x < 0.0f ? PI_F - result : result;
}

/**
 * @brief Fast exp for non-positive arguments.
 * @param x Argument; the domain is x <= 0.
 * @return An approximation of e^x.
 * @details exp(x) = 2^(x*log2e): an integer power from the exponent bits times
 * a quartic for the fraction, continuous across exponent boundaries. Peak
 * relative error ~5.1e-6 over [-30, 0]; large-magnitude x saturates to 0.
 * @warning The domain is guarded only by an audit check, but every argument
 * still yields a defined result: a positive x stays accurate to x ~ 88.0 and
 * saturates to +Inf above it, and a NaN returns 0.
 */
inline float fast_expf(float x) {
  HS_AUDIT_CHECK(x <= 0.0f);
  float y = x * 1.442695041f; // log2(e)
  float fi = floorf(y);
  float f = y - fi; // [0, 1)
  float p =
      1.0f + f * (0.69303227f +
                  f * (0.24137937f + f * (0.05203249f + f * 0.01355588f)));
  // Bounds the packed exponent on both sides so (int)fi cannot overflow; the
  // negated compare also takes the NaN branch.
  if (!(std::fabs(fi) <= 126.0f))
    return fi > 0.0f ? std::numeric_limits<float>::infinity() : 0.0f;
  int i = (int)fi;
  uint32_t bits = (uint32_t)(i + 127) << 23;
  float scale;
  std::memcpy(&scale, &bits, sizeof(scale));
  return p * scale;
}

/**
 * @brief Bhaskara I sine approximation for an angle in [0, π].
 * @param x Angle in radians, in [0, π].
 * @param sign Sign of the half-period `x` was folded out of, +1 or -1.
 * @return The approximate sine, `sign * sin(x)`.
 * @details Shared by the fast sine and cosine paths so they cannot be
 * reassociated apart under -ffast-math.
 */
__attribute__((always_inline)) inline float bhaskara_sinf(float x, float sign) {
  float xpi = x * (PI_F - x);
  return sign * (16.0f * xpi) / (5.0f * PI_F * PI_F - 4.0f * xpi);
}

/**
 * @brief Fast sine for an angle already reduced into [0, 2π).
 * @param x Angle in radians, in [0, 2π).
 * @return The approximate sine of `x`.
 */
__attribute__((always_inline)) inline float sinf_0_2pi(float x) {
  float sign = 1.0f;
  if (x > PI_F) {
    x -= PI_F;
    sign = -1.0f;
  }
  return bhaskara_sinf(x, sign);
}

/**
 * @brief Fast sine using the Bhaskara I approximation.
 * @param x Angle in radians (range-reduced internally).
 * @return The approximate sine of `x`.
 * @details Absolute error <= 1.7e-3 (~0.1 degrees); relative error is large
 * near the zeros.
 * @warning Accuracy degrades for large |x|: the `x - floor(x/2π)·2π` range
 * reduction loses precision as a float's ULP grows past the period.
 */
__attribute__((always_inline)) inline float fast_sinf(float x) {
  constexpr float INV_2_PI = 1.0f / (2.0f * PI_F);
  return sinf_0_2pi(x - floorf(x * INV_2_PI) * (2.0f * PI_F));
}

/**
 * @brief Fast cosine, equal to fast_sinf(x + π/2).
 * @param x Angle in radians.
 * @return The approximate cosine of `x`.
 * @details Reduces x into [-π/2, 3π/2) and applies the quarter-turn shift
 * last, so on [0, π] it matches fast_sincosf_0_pi bit for bit.
 * @warning Accuracy degrades for large |x|: the `x - floor(x/2π)·2π` range
 * reduction loses precision as a float's ULP grows past the period.
 */
__attribute__((always_inline)) inline float fast_cosf(float x) {
  constexpr float INV_2_PI = 1.0f / (2.0f * PI_F);
  float r = x - floorf(x * INV_2_PI + 0.25f) * (2.0f * PI_F);
  return sinf_0_2pi(r + PI_F * 0.5f);
}

/**
 * @brief Fast sine and cosine for an angle in [0, π].
 * @param x Angle in radians.
 * @param s Output matching fast_sinf(x) bit for bit.
 * @param c Output matching fast_cosf(x) bit for bit.
 * @details The domain makes both range reductions exact no-ops, so this drops
 * them and calls the same kernels the general functions do.
 * @warning The sine branch has no range reduction and no sign tracking; for
 * x outside [0, π] it returns the wrong value, so callers must bound x.
 */
__attribute__((always_inline)) inline void fast_sincosf_0_pi(float x, float &s,
                                                             float &c) {
  HS_AUDIT_CHECK(x >= 0.0f && x <= PI_F);
  s = bhaskara_sinf(x, 1.0f);
  c = sinf_0_2pi(x + PI_F * 0.5f);
}

/**
 * @brief Divide-free sine and cosine for an angle in [0, π].
 * @param x Angle in radians.
 * @param s Output sine.
 * @param c Output cosine.
 * @details Evaluates minimax polynomials of r = x - π/2: degree 8 in r for
 * the sine and degree 9 for the cosine. Absolute error is below 3e-7 across
 * the domain.
 * @warning For x outside [0, π] the error grows with the polynomial tail;
 * callers must bound x.
 */
__attribute__((always_inline)) inline void poly_sincosf_0_pi(float x, float &s,
                                                             float &c) {
  HS_AUDIT_CHECK(x >= 0.0f && x <= PI_F);
  const float r = x - PI_F * 0.5f;
  const float r2 = r * r;
  s = 0.99999995329f +
      r2 * (-0.49999905115f +
            r2 * (0.041663579924f +
                  r2 * (-0.0013853672957f + r2 * 2.3153287771e-05f)));
  c = -r * (0.99999997653f +
            r2 * (-0.16666647600f +
                  r2 * (0.0083328993172f +
                        r2 * (-0.00019800870025f + r2 * 2.5904378421e-06f))));
}

/**
 * @brief Normalized linear interpolation between two directions.
 * @param a Starting direction.
 * @param b Ending direction.
 * @param t Interpolation factor (0.0 to 1.0).
 * @return The renormalized blend, or `a` when the blend's squared length falls
 *         below math::EPS_BLEND_LEN_SQ.
 */
inline Vector nlerp_unit(const Vector &a, const Vector &b, float t) {
  const Vector mixed(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t,
                     a.z + (b.z - a.z) * t);
  const float length_sq =
      mixed.x * mixed.x + mixed.y * mixed.y + mixed.z * mixed.z;
  if (length_sq < math::EPS_BLEND_LEN_SQ)
    return a;
  const float inv_length = 1.0f / sqrtf(length_sq);
  return Vector(mixed.x * inv_length, mixed.y * inv_length,
                mixed.z * inv_length);
}

/**
 * @brief Approximate spherical interpolation between two vectors.
 * @note Uses approximate arccosine and Bhaskara sine weights, followed by
 * normalization; intermediate angular speed is approximate.
 * @param v1 Starting vector.
 * @param v2 Ending vector.
 * @param t Interpolation factor (0.0 to 1.0).
 * @pre Both endpoints are unit length, checked in audit builds.
 * @return The interpolated unit vector.
 */
inline Vector slerp(const Vector &v1, const Vector &v2, float t) {
  // Negated compare: a NaN endpoint passes and traps in normalized() below with
  // its own diagnostic instead of aborting on a unit-length report.
  HS_AUDIT_CHECK(!(std::fabs(dot(v1, v1) - 1.0f) > math::EPS_UNIT_VEC_SQ));
  HS_AUDIT_CHECK(!(std::fabs(dot(v2, v2) - 1.0f) > math::EPS_UNIT_VEC_SQ));
  float d = hs::clamp(dot(v1, v2), -1.0f, 1.0f);
  // Near-parallel: lerp to avoid NaN. Branch thresholds are TOLERANCE-scale.
  if (d > 1.0f - math::TOLERANCE) {
    return (v1 + (v2 - v1) * t).normalized();
  }
  if (d < -1.0f + math::TOLERANCE &&
      (d <= -1.0f || distance_squared(v1, -v2) < math::EPS_BLEND_LEN_SQ)) {
    // The endpoint sum detects cancellation without acos amplifying unit error.
    Vector axis = cross(least_parallel_axis(v1), v1).normalized();
    return rotate(v1, make_rotation(axis, t * PI_F));
  }
  float theta = fast_acos(d);
  // The 1/sinθ factor cancels under the final normalize().
  float s1 = fast_sinf((1 - t) * theta);
  float s2 = fast_sinf(t * theta);
  return ((s1 * v1) + (s2 * v2)).normalized();
}

/**
 * @brief A slerp arc between two unit quaternions, prepared once and sampled
 * at any number of fractions.
 */
struct SlerpArc {
  Quaternion from;    ///< Start orientation, sign-adjusted onto the arc.
  Quaternion to;      ///< End orientation.
  float theta = 0.0f; ///< Arc angle, or zero for normalized linear sampling.
  float sine = 1.0f;  ///< Sine of the arc angle.

  /**
   * @brief Prepares the shortest arc between unit orientations.
   * @param a Start unit orientation.
   * @param b End unit orientation.
   */
  __attribute__((always_inline)) SlerpArc(const Quaternion &a,
                                          const Quaternion &b)
      : from(a), to(b) {
    HS_AUDIT_CHECK(!(std::fabs(dot(a, a) - 1.0f) > math::EPS_UNIT_QUAT_SQ));
    HS_AUDIT_CHECK(!(std::fabs(dot(b, b) - 1.0f) > math::EPS_UNIT_QUAT_SQ));
    float d = dot(a, b);
    if (d < 0.0f) {
      from = -from;
      d = -d;
    }
    prepare(d);
  }

  /**
   * @brief Prepares the arc from @p a to @p b, with @p a already on it.
   * @param a Start unit orientation.
   * @param b End unit orientation.
   * @param d `dot(a, b)`.
   */
  SlerpArc(const Quaternion &a, const Quaternion &b, float d) : from(a), to(b) {
    prepare(d);
  }

  /**
   * @brief Samples the arc.
   * @param t Fraction of the arc, 0 at `from` and 1 at `to`.
   * @return The unit orientation at @p t.
   */
  __attribute__((always_inline)) Quaternion sample(float t) const {
    if (theta == 0.0f)
      return (from + t * (to - from)).normalized();
    const float s1 = sinf((1.0f - t) * theta) / sine;
    const float s2 = sinf(t * theta) / sine;
    return ((s1 * from) + (s2 * to)).normalized();
  }

private:
  /**
   * @brief Sets the arc trigonometry unless the endpoints are near-parallel.
   * @param d `dot(from, to)`.
   */
  __attribute__((always_inline)) void prepare(float d) {
    if (d <= 1.0f - math::TOLERANCE) {
      theta = acosf(hs::clamp(d, -1.0f, 1.0f));
      sine = sinf(theta);
    }
  }
};

/**
 * @brief Spherical Linear Interpolation (SLERP) between two quaternions.
 * @param q1 Starting quaternion.
 * @param q2 Ending quaternion.
 * @param t Interpolation factor (0.0 to 1.0).
 * @param long_way If true, takes the longest path between quaternions.
 * @pre Both endpoints are unit length; a non-unit pair skews `d` and lands the
 * blend on the wrong arc.
 * @return The interpolated unit quaternion.
 */
inline Quaternion slerp(const Quaternion &q1, const Quaternion &q2, float t,
                        bool long_way = false) {
  HS_AUDIT_CHECK(!(std::fabs(dot(q1, q1) - 1.0f) > math::EPS_UNIT_QUAT_SQ));
  HS_AUDIT_CHECK(!(std::fabs(dot(q2, q2) - 1.0f) > math::EPS_UNIT_QUAT_SQ));
  float d = dot(q1, q2);
  Quaternion p(q1);

  // At d == 0 both arcs are 90° on S³; long_way reverses the 180° rotation.
  if ((long_way && d >= 0) || (!long_way && d < 0)) {
    p = -p;
    d = -d;
  }

  if (d < -1.0f + math::TOLERANCE &&
      (d <= -1.0f || (p + q2).squared_magnitude() < math::EPS_BLEND_LEN_SQ)) {
    Quaternion n(-p.v.x, p.r, -p.v.z, p.v.y);
    return (cosf(t * PI_F) * p + sinf(t * PI_F) * n).normalized();
  }
  return SlerpArc(p, q2, d).sample(t);
}

/**
 * @brief Scales a rotation delta down to a fraction of its arc.
 * @param delta Full rotation delta.
 * @param amount Fraction of the arc to keep.
 * @pre @p delta is unit length; the slerp path it takes for a fractional
 * @p amount asserts on a non-unit endpoint.
 * @return @p delta at @p amount of its rotation; identity at 0, @p delta at 1.
 */
HS_HOT_FLASH_MEMBER inline Quaternion
scaled_rotation_delta(const Quaternion &delta, float amount) {
  if (amount == 1.0f)
    return delta;
  if (amount == 0.0f)
    return Quaternion();
  return slerp(Quaternion(), delta, amount);
}

/** @brief Accumulates a fraction of a unit rotation delta.
 * @param delta Unit rotation delta.
 * @param amount Fraction of the rotation to absorb.
 * @param accumulated Rotation accumulated before this delta.
 * @return Normalized rotation after absorbing the scaled delta. */
inline Quaternion absorb_rotation_delta(const Quaternion &delta, float amount,
                                        const Quaternion &accumulated) {
  return (scaled_rotation_delta(delta, amount) * accumulated).normalized();
}

/** @brief Reduces a bounded displacement with nearest-even period ties.
 * @param value Displacement with an integer-representable period quotient.
 * @param period Positive period.
 * @return Remainder in [-period/2, period/2]. */
HS_HOT_INLINE inline float centered_remainder(float value, float period) {
  const int quotient = static_cast<int>(value / period);
  float result = value - static_cast<float>(quotient) * period;
  const float half = period * 0.5f;
  if (result > half || (result == half && (quotient & 1)))
    result -= period;
  else if (result < -half || (result == -half && (quotient & 1)))
    result += period;
  return result;
}

/**
 * @brief Sine of a double, by its Taylor series for a small argument.
 * @param x Angle, radians.
 * @return `sin(x)`, well within float rounding.
 */
HS_HOT_INLINE inline double series_sin(double x) {
  if (!(fabs(x) < 0.5))
    return sin(x);
  const double x2 = x * x;
  return x *
         (1.0 +
          x2 * (-1.0 / 6.0 +
                x2 * (1.0 / 120.0 +
                      x2 * (-1.0 / 5040.0 + x2 * (1.0 / 362880.0 +
                                                  x2 * (-1.0 / 39916800.0))))));
}

/**
 * @brief Cosine of a double, by its Taylor series for a small argument.
 * @param x Angle, radians.
 * @return `cos(x)`, well within float rounding.
 */
HS_HOT_INLINE inline double series_cos(double x) {
  if (!(fabs(x) < 0.5))
    return cos(x);
  const double x2 = x * x;
  return 1.0 + x2 * (-1.0 / 2.0 +
                     x2 * (1.0 / 24.0 +
                           x2 * (-1.0 / 720.0 +
                                 x2 * (1.0 / 40320.0 +
                                       x2 * (-1.0 / 3628800.0 +
                                             x2 * (1.0 / 479001600.0))))));
}

} // namespace math
