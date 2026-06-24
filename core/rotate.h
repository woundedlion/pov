/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <cmath>
#include <cstdint>
#include "3dmath.h"

/** @brief Two Pi. */
static constexpr float tau = 2 * PI_F;
/** @brief Multiplication factor to convert degrees to radians. */
static constexpr float radians = PI_F / 180;
/** @brief Multiplication factor to convert radians to degrees. */
static constexpr float degrees = 180 / PI_F;

/**
 * @brief Wraps an angle of any magnitude into [0, 2PI).
 * @param n Angle in radians.
 * @return Equivalent angle in radians within [0, 2PI).
 */
inline float mod_tau(float n) { return n - floorf(n / tau) * tau; }

/**
 * @brief Wraps a floating-point index into [0, m), preserving the fraction.
 * @param x Source index, may be negative.
 * @param m Modulus (exclusive upper bound), must be positive.
 * @return Wrapped index in [0, m) with the original fractional part retained.
 * @details Uses floorf (not a truncating cast) so negative inputs floor toward
 *  -inf, and the double-mod `((i % m) + m) % m` keeps the integer part
 *  non-negative. This guarantees the [0, m) contract even for x < 0, avoiding a
 *  negative pixel x downstream.
 */
inline float wrap_index(float x, int m) {
  int i = static_cast<int>(floorf(x));
  int w = ((i % m) + m) % m;
  return w + (x - i);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Handles equirectangular projection rotation.
 * @tparam W Image width in pixels.
 * @tparam H Image height in pixels.
 * @details Used for rotating the background skybox/environment.
 */
template <uint8_t W, uint8_t H> class Projection {
public:
  /**
   * @brief A pixel in the WxH equirectangular image paired with its spherical
   *  coordinates.
   * @details lambda (longitude) lies in [-PI, PI) and phi (latitude) in
   *  [-PI/2, PI/2]. The constructors derive lambda/phi from the pixel (x, y);
   *  y is measured top-down, so phi is flipped via (H - y).
   */
  struct Point {
    /**
     * @brief Constructs a Point from integer pixel coordinates.
     * @param x Column index in [0, W).
     * @param y Row index in [0, H), measured top-down.
     */
    Point(int x, int y)
        : x(x), y(y), lambda(x * tau / W - PI_F),
          phi((H - y) * PI_F / H - PI_F / 2) {}

    /**
     * @brief Constructs a Point from 8-bit pixel coordinates.
     * @param x Column index in [0, W).
     * @param y Row index in [0, H), measured top-down.
     */
    Point(uint8_t x, uint8_t y)
        : x(x), y(y), lambda(x * tau / W - PI_F),
          phi((H - y) * PI_F / H - PI_F / 2) {}

    /**
     * @brief Constructs a Point from floating-point pixel coordinates.
     * @param x Column index in [0, W).
     * @param y Row index in [0, H), measured top-down.
     */
    Point(float x, float y)
        : x(x), y(y), lambda(x * tau / W - PI_F),
          phi((H - y) * PI_F / H - PI_F / 2) {}

    /**
     * @brief Copy-constructs a Point.
     * @param p Source point to copy.
     */
    Point(const Point &p) : x(p.x), y(p.y), lambda(p.lambda), phi(p.phi) {}

    /**
     * @brief Rounds x to the nearest pixel column, wrapped into [0, W).
     * @return Column index in [0, W).
     */
    int xi() { return static_cast<int>(x + 0.5f) % W; }
    /**
     * @brief Rounds y to the nearest pixel row.
     * @return Row index.
     */
    int yi() { return static_cast<int>(y + 0.5f); }

    float x;            /**< Pixel column coordinate. */
    float y;            /**< Pixel row coordinate. */
    float lambda;       /**< Longitude in [-PI, PI). */
    float phi;          /**< Latitude in [-PI/2, PI/2]. */
  };

  /** @brief Constructs an identity projection (no rotation applied). */
  Projection() {}

  /**
   * @brief Projects pixel (bx, by) through the accumulated rotation.
   * @param bx Source column index in [0, W).
   * @param by Source row index in [0, H).
   * @return The rotated Point whose x/y give the source pixel to sample.
   */
  Point project(uint8_t bx, uint8_t by) const {
    Point p(bx, by);
    return project(p);
  }

  /**
   * @brief Maps a source point to its rotated equirectangular pixel.
   * @param src Destination point carrying its lambda/phi spherical coordinates.
   * @return A Point whose x/y are the source pixel to sample for this
   *  destination.
   * @details Applies the lambda offset, converts to a Cartesian unit vector,
   *  rotates by the cached phi/gamma terms, then projects back to (x, y).
   */
  Point project(const Point &src) const {
    Point p(src);

    p.lambda += delta_lambda;
    if (p.lambda > PI_F) {
      p.lambda -= tau;
    } else if (p.lambda < -PI_F) {
      p.lambda += tau;
    }

    // lambda/phi -> Cartesian unit vector
    float cos_p = cosf(p.phi);
    float x = cosf(p.lambda) * cos_p;
    float y = sinf(p.lambda) * cos_p;
    float z = sinf(p.phi);

    // rotate by the cached phi/gamma terms
    float k = z * cos_dp + x * sin_dp;
    p.lambda = atan2f(y * cos_dg - k * sin_dg, x * cos_dp - z * sin_dp);
    p.phi = asinf(k * cos_dg + y * sin_dg);

    // back to equirectangular pixel (x, y)
    p.x = wrap_index((p.lambda + PI_F) * W / tau, W);
    p.y = H - ((p.phi + PI_F / 2) * H / PI_F);

    return p;
  }

  /**
   * @brief Accumulates rotation by (lambda, phi, gamma) increments.
   * @param dl Longitude increment in degrees.
   * @param dp Latitude increment in degrees.
   * @param dg Roll increment in degrees.
   * @return Reference to this Projection for chaining.
   * @details Caches the phi/gamma sin/cos terms used by project().
   */
  Projection &rotate(uint16_t dl, uint16_t dp, uint16_t dg) {
    delta_lambda = mod_tau(delta_lambda + ((dl % 360) * radians));
    delta_phi = mod_tau(delta_phi + ((dp % 360) * radians));
    delta_gamma = mod_tau(delta_gamma + ((dg % 360) * radians));

    cos_dp = cosf(delta_phi);
    sin_dp = sinf(delta_phi);
    cos_dg = cosf(delta_gamma);
    sin_dg = sinf(delta_gamma);

    return *this;
  }

  /**
   * @brief Clears the accumulated rotation back to identity.
   * @return Reference to this Projection for chaining.
   */
  Projection &reset() {
    delta_lambda = 0;
    delta_phi = 0;
    delta_gamma = 0;

    cos_dp = cosf(delta_phi);
    sin_dp = sinf(delta_phi);
    cos_dg = cosf(delta_gamma);
    sin_dg = sinf(delta_gamma);

    return *this;
  }

private:
  float delta_lambda = 0;          /**< Accumulated longitude offset in radians. */
  float delta_phi = 0;             /**< Accumulated latitude offset in radians. */
  float delta_gamma = 0;           /**< Accumulated roll offset in radians. */
  float cos_dp = cosf(delta_phi);  /**< Cached cosine of delta_phi. */
  float sin_dp = sinf(delta_phi);  /**< Cached sine of delta_phi. */
  float cos_dg = cosf(delta_gamma);/**< Cached cosine of delta_gamma. */
  float sin_dg = sinf(delta_gamma);/**< Cached sine of delta_gamma. */
};
