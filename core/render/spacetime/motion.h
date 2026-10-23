/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file motion.h
 * @brief Worldline interpolants and analytic path bounds.
 */

#include "core/math/3dmath.h"
#include "core/platform/attributes.h"
#include <cmath>
#include <cstdint>

namespace Spacetime {

/** @brief How a vertex's path over the frame is interpolated. */
enum class Interpolant : uint8_t {
  CHORD,     ///< Straight between the frame's start and end captures.
  QUADRATIC, ///< The parabola through the start, middle and end captures.
  ARC,       ///< A circle about an axis, at a constant angular rate.
};

/** @brief Bounds on a vertex path's derivatives over one frame. */
struct PathBound {
  float second = 0.0f; ///< Bound on `|p''(t)|` for t in [0, 1].
  float third = 0.0f;  ///< Bound on `|p'''(t)|` for t in [0, 1].

  /**
   * @brief Bounds of a rigid rotation at a constant rate.
   * @param angle Rotation over the frame, radians.
   * @param radius Distance of the vertex from the rotation axis.
   * @return `{angle^2 radius, |angle|^3 radius}`.
   */
  static constexpr PathBound rotation(float angle, float radius = 1.0f) {
    const float a = angle < 0.0f ? -angle : angle;
    return {a * a * radius, a * a * a * radius};
  }

  /**
   * @brief Bounds of a rotation about a fixed axis by a varying angle.
   * @param rate Bound on the angle's first derivative over the frame.
   * @param acceleration Bound on its second derivative.
   * @param jerk Bound on its third derivative.
   * @param radius Distance of the point from the axis.
   * @return `{radius (rate^2 + acceleration),
   *         radius (rate^3 + 3 rate acceleration + jerk)}`.
   * @details With `q = a x p`, `p' = angle' q` and `q' = angle' (a x q)`, so
   * `p''` and `p'''` are sums of terms along `q` and `a x q`, each of length
   * `radius`.
   */
  static constexpr PathBound turning(float rate, float acceleration, float jerk,
                                     float radius) {
    const float r = rate < 0.0f ? -rate : rate;
    const float a = acceleration < 0.0f ? -acceleration : acceleration;
    const float j = jerk < 0.0f ? -jerk : jerk;
    return {radius * (r * r + a), radius * (r * r * r + 3.0f * r * a + j)};
  }

  /**
   * @brief Bounds of a unit-sphere path moving linearly in Mercator
   *        coordinates.
   * @param ordinate Change of the Mercator ordinate `u` over the frame, the
   *        point sitting at height `tanh u` on its axis.
   * @param longitude Change of longitude about that axis, radians.
   * @return `{m^2, 2 m^3}` with `m = hypot(ordinate, longitude)`.
   * @details Along such a path, a loxodrome, `|p''| = m^2 sech u` exactly;
   * `|p'''|` peaks at `2 m^3`, on the equator under purely meridional motion.
   */
  static PathBound loxodrome(float ordinate, float longitude) {
    const float m2 = ordinate * ordinate + longitude * longitude;
    return {m2, 2.0f * m2 * sqrtf(m2)};
  }

  /**
   * @brief Bounds of a point turning along a great circle while a moving
   *        plane translation carries it.
   * @param turn Angle the point turns along its great circle over the frame,
   *        radians, at a constant rate.
   * @param shift Change over the frame of the translation `w -> w + b` of
   *        the unit sphere's stereographic plane, `b` moving linearly.
   * @param reach Bound on `|b|` over the frame.
   * @return `{cosh(d) v^2, 1.5 cosh(d) v^3}` with `cosh(d) = 1 + reach^2 / 2`
   *         and `v = e^d |turn| + 2 shift`.
   * @details The translation moves the ball's centre a hyperbolic distance
   * `d`, so it stretches arcs by at most `e^d`, and a change of `b` moves a
   * point at most twice as far, so `v` bounds the path's speed. A fixed
   * translation of a great circle reaches `|p''| = cosh(d) v^2`; a still point
   * under a moving one reaches `|p'''| = 1.5 v^3` at the plane's origin. The
   * bounds cover both motions together, checked numerically.
   */
  static PathBound warped_great_circle(float turn, float shift, float reach) {
    const float cosh_d = 1.0f + 0.5f * reach * reach;
    const float stretch = cosh_d + sqrtf(cosh_d * cosh_d - 1.0f);
    const float v = stretch * fabsf(turn) + 2.0f * shift;
    return {cosh_d * v * v, 1.5f * cosh_d * v * v * v};
  }

  /**
   * @brief Bound on the chord's distance from the path.
   * @return An eighth of the second-derivative bound.
   */
  constexpr float chord_error() const { return 0.125f * second; }

  /**
   * @brief Bound on the parabola's distance from the path.
   * @return `third * sqrt(3) / 216`: the third-derivative bound over 6 times
   *         the largest `|t (t - 1/2) (t - 1)|` on [0, 1].
   */
  constexpr float quadratic_error() const { return third * 0.0080187537f; }

  /**
   * @brief Bound on the normalized chord's angle from a unit vector turning
   *        at a constant rate, up to a quarter turn.
   * @return `third * 0.0186`: the normalized chord leads and trails the turn
   *         by `a^3 / (36 sqrt(3))` to leading order in the turn `a`, and by
   *         under 0.0186 a^3 up to a quarter turn.
   */
  constexpr float normalized_chord_error() const { return third * 0.0186f; }
};

/**
 * @brief Picks a path's interpolant from its analytic bound.
 * @param bound Derivative bounds of the path over the frame.
 * @param tolerance Largest acceptable distance from the path.
 * @param on_arc Whether the path is a circle about a known axis at a constant
 *        rate, which Interpolant::ARC reproduces.
 * @return CHORD when its bound meets @p tolerance, else ARC for a path on an
 *         arc, else QUADRATIC, which needs the middle capture.
 */
constexpr Interpolant choose_interpolant(const PathBound &bound,
                                         float tolerance, bool on_arc) {
  if (bound.chord_error() <= tolerance)
    return Interpolant::CHORD;
  return on_arc ? Interpolant::ARC : Interpolant::QUADRATIC;
}

/**
 * @brief Picks the interpolant of a unit vector turning at a constant rate.
 * @param angle Turn over the frame, radians.
 * @param tolerance Largest acceptable angle from the turn.
 * @return CHORD, normalized where evaluated, when its bound meets
 *         @p tolerance within a quarter turn, else ARC, which reproduces the
 *         turn.
 */
constexpr Interpolant choose_turn_interpolant(float angle, float tolerance) {
  return angle <= 0.5f * math::PI_F &&
                 PathBound::rotation(angle).normalized_chord_error() <=
                     tolerance
             ? Interpolant::CHORD
             : Interpolant::ARC;
}

/**
 * @brief Bound on an interpolant's distance from the path it stands for.
 * @param kind The interpolant.
 * @param bound Derivative bounds of the path.
 * @return The chord or parabola bound; zero for ARC, which is exact on its
 *         circle, up to the captures' rounding.
 */
constexpr float interpolant_error(Interpolant kind, const PathBound &bound) {
  return kind == Interpolant::CHORD       ? bound.chord_error()
         : kind == Interpolant::QUADRATIC ? bound.quadratic_error()
                                          : 0.0f;
}

/**
 * @brief Bound on the displacement of any unit vector between two rotations.
 * @param previous Rotation at the previous endpoint.
 * @param current Rotation at the current endpoint.
 * @return Frobenius norm of their difference, which bounds the chord between
 *         a unit vector's two images.
 */
HS_HOT_INLINE inline float
rotation_chord_bound(const math::RotationMatrix &previous,
                     const math::RotationMatrix &current) {
  const math::Vector a = previous.r0 - current.r0;
  const math::Vector b = previous.r1 - current.r1;
  const math::Vector c = previous.r2 - current.r2;
  return sqrtf(math::dot(a, a) + math::dot(b, b) + math::dot(c, c));
}

/**
 * @brief One vertex's interpolated path over a frame.
 * @details Built once per vertex per frame from its captures. Points are not
 * renormalized; a caller on the sphere normalizes what it projects.
 */
struct Worldline {
  Interpolant kind = Interpolant::CHORD; ///< Interpolant evaluated by at().
  math::Vector start;                    ///< Position at t = 0.
  math::Vector end;                      ///< Position at t = 1.
  math::Vector middle;  ///< Position at t = 1/2, for QUADRATIC.
  math::Vector axis;    ///< Unit circle axis, for ARC.
  math::Vector radial;  ///< Unit direction from the axis to the start, for ARC.
  math::Vector across;  ///< `axis x radial`, for ARC.
  float angle = 0.0f;   ///< Turn about the axis over the frame, for ARC.
  float height0 = 0.0f; ///< Start's height along the axis, for ARC.
  float height1 = 0.0f; ///< End's height along the axis, for ARC.
  float radius0 = 0.0f; ///< Start's distance from the axis, for ARC.
  float radius1 = 0.0f; ///< End's distance from the axis, for ARC.

  /**
   * @brief The straight path between two captures.
   * @param start Position at t = 0.
   * @param end Position at t = 1.
   * @return A CHORD path.
   */
  static Worldline chord(const math::Vector &start, const math::Vector &end) {
    Worldline path;
    path.start = start;
    path.end = end;
    path.middle = (start + end) * 0.5f;
    return path;
  }

  /**
   * @brief The parabola through three captures.
   * @param start Position at t = 0.
   * @param middle Position at t = 1/2.
   * @param end Position at t = 1.
   * @return A QUADRATIC path.
   */
  static Worldline quadratic(const math::Vector &start,
                             const math::Vector &middle,
                             const math::Vector &end) {
    Worldline path = chord(start, end);
    path.kind = Interpolant::QUADRATIC;
    path.middle = middle;
    return path;
  }

  /**
   * @brief The circular path about an axis between two captures.
   * @param start Position at t = 0.
   * @param end Position at t = 1.
   * @param axis Unit axis, or null for the great circle through both, whose
   *        axis is `start x end`.
   * @return An ARC path turning less than half a turn, its height and radius
   *         about the axis blended linearly so it meets both captures; a
   *         CHORD path when either capture lies on the axis.
   */
  static Worldline arc(const math::Vector &start, const math::Vector &end,
                       const math::Vector *axis) {
    Worldline path = chord(start, end);
    math::Vector k;
    if (axis) {
      k = *axis;
    } else {
      k = math::cross(start, end);
      const float length = sqrtf(math::dot(k, k));
      if (!(length > 0.0f))
        return path;
      k = k * (1.0f / length);
    }
    const float h0 = math::dot(start, k), h1 = math::dot(end, k);
    const math::Vector u = start - k * h0, v = end - k * h1;
    const float r0 = sqrtf(math::dot(u, u)), r1 = sqrtf(math::dot(v, v));
    if (!(r0 > 0.0f && r1 > 0.0f))
      return path;
    path.kind = Interpolant::ARC;
    path.axis = k;
    path.radial = u * (1.0f / r0);
    path.across = math::cross(k, path.radial);
    path.angle = atan2f(math::dot(path.across, v), math::dot(path.radial, v));
    path.height0 = h0;
    path.height1 = h1;
    path.radius0 = r0;
    path.radius1 = r1;
    const math::Vector halfway = path.at(0.5f);
    path.middle = halfway;
    return path;
  }

  /**
   * @brief The path's position at a frame fraction.
   * @param t Frame fraction in [0, 1].
   * @return The interpolated position.
   */
  HS_HOT_INLINE math::Vector at(float t) const {
    switch (kind) {
    case Interpolant::CHORD:
      return start + (end - start) * t;
    case Interpolant::QUADRATIC:
      return start * ((1.0f - t) * (1.0f - 2.0f * t)) +
             middle * (4.0f * t * (1.0f - t)) + end * (t * (2.0f * t - 1.0f));
    case Interpolant::ARC:
      break;
    }
    const float turn = angle * t;
    const float radius = radius0 + (radius1 - radius0) * t;
    return axis * (height0 + (height1 - height0) * t) +
           (radial * cosf(turn) + across * sinf(turn)) * radius;
  }

  /**
   * @brief The path under a fixed rotation.
   * @param m Rotation.
   * @return The same path with every point rotated by @p m.
   */
  HS_HOT_INLINE Worldline rotated(const math::RotationMatrix &m) const {
    Worldline path = *this;
    path.start = m.apply(start);
    path.end = m.apply(end);
    path.middle = m.apply(middle);
    path.axis = m.apply(axis);
    path.radial = m.apply(radial);
    path.across = m.apply(across);
    return path;
  }

  /**
   * @brief Bound on the path's distance from its chord.
   * @return A bound on `|at(t) - chord(t)|` over t in [0, 1], `chord(t)` the
   *         point at t on the straight path between the captures: zero for
   *         CHORD, half the parabola's second difference for QUADRATIC, and
   *         for ARC the larger radius times an eighth of the squared turn
   *         plus a quarter of the radius change times the turn's chord.
   */
  HS_HOT_INLINE float deviation() const {
    switch (kind) {
    case Interpolant::CHORD:
      return 0.0f;
    case Interpolant::QUADRATIC: {
      const math::Vector second = middle * 2.0f - start - end;
      return 0.5f * sqrtf(math::dot(second, second));
    }
    case Interpolant::ARC:
      break;
    }
    return fmaxf(radius0, radius1) * (0.125f * angle * angle) +
           0.5f * fabsf(radius1 - radius0) * fabsf(sinf(0.5f * angle));
  }
};

} // namespace Spacetime
