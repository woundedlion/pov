/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <cmath>
#include <cstdint>
#include "3dmath.h"

/** @brief Two Pi */
static constexpr float tau = 2 * PI_F;
/** @brief Multiplication factor to convert degrees to radians. */
static constexpr float radians = PI_F / 180;
/** @brief Multiplication factor to convert radians to degrees. */
static constexpr float degrees = 180 / PI_F;

/** @brief Wraps an angle of any magnitude into [0, 2PI). */
inline float mod_tau(float n) { return n - floorf(n / tau) * tau; }

/** @brief Wraps a floating-point index into [0, m), preserving the fraction.
 *
 * Uses floorf (not a truncating cast) so negative inputs floor toward -inf, and
 * the double-mod `((i % m) + m) % m` keeps the integer part non-negative. This
 * guarantees the [0, m) contract even for x < 0, avoiding a negative pixel x
 * downstream. */
inline float wrap_index(float x, int m) {
  int i = static_cast<int>(floorf(x));
  int w = ((i % m) + m) % m;
  return w + (x - i);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Handles equirectangular projection rotation.
 * Used for rotating the background skybox/environment.
 */
template <uint8_t W, uint8_t H> class Projection {
public:
  /** @brief A pixel in the WxH equirectangular image paired with its spherical
   *  coordinates: lambda (longitude) in [-PI, PI) and phi (latitude) in
   *  [-PI/2, PI/2]. The constructors derive lambda/phi from the pixel (x, y);
   *  y is measured top-down, so phi is flipped via (H - y). */
  struct Point {
    Point(int x, int y)
        : x(x), y(y), lambda(x * tau / W - PI_F),
          phi((H - y) * PI_F / H - PI_F / 2) {}

    Point(uint8_t x, uint8_t y)
        : x(x), y(y), lambda(x * tau / W - PI_F),
          phi((H - y) * PI_F / H - PI_F / 2) {}

    Point(float x, float y)
        : x(x), y(y), lambda(x * tau / W - PI_F),
          phi((H - y) * PI_F / H - PI_F / 2) {}

    Point(const Point &p) : x(p.x), y(p.y), lambda(p.lambda), phi(p.phi) {}

    /** @brief Rounds x to the nearest pixel column, wrapped into [0, W). */
    int xi() { return static_cast<int>(x + 0.5f) % W; }
    /** @brief Rounds y to the nearest pixel row. */
    int yi() { return static_cast<int>(y + 0.5f); }

    float x;
    float y;
    float lambda; ///< longitude in [-PI, PI)
    float phi;    ///< latitude in [-PI/2, PI/2]
  };

  /** @brief Constructs an identity projection (no rotation applied). */
  Projection() {}

  /** @brief Projects pixel (bx, by) through the accumulated rotation. */
  Point project(uint8_t bx, uint8_t by) const {
    Point p(bx, by);
    return project(p);
  }

  /** @brief Maps a source point to its rotated equirectangular pixel.
   *
   * Applies the lambda offset, converts to a Cartesian unit vector, rotates by
   * the cached phi/gamma terms, then projects back to (x, y). The returned
   * Point's x/y are the source pixel to sample for this destination. */
  Point project(const Point &src) const {
    Point p(src);

    // rotate lambda
    p.lambda += delta_lambda;
    if (p.lambda > PI_F) {
      p.lambda -= tau;
    } else if (p.lambda < -PI_F) {
      p.lambda += tau;
    }

    // convert to cartesian x, y, z
    float cos_p = cosf(p.phi);
    float x = cosf(p.lambda) * cos_p;
    float y = sinf(p.lambda) * cos_p;
    float z = sinf(p.phi);

    // rotate phi gamma
    float k = z * cos_dp + x * sin_dp;
    p.lambda = atan2f(y * cos_dg - k * sin_dg, x * cos_dp - z * sin_dp);
    p.phi = asinf(k * cos_dg + y * sin_dg);

    // convert to equirectangular x, y
    p.x = wrap_index((p.lambda + PI_F) * W / tau, W);
    p.y = H - ((p.phi + PI_F / 2) * H / PI_F);

    return p;
  }

  /** @brief Accumulates rotation by (lambda, phi, gamma) increments in
   *  degrees, then caches the phi/gamma sin/cos used by project(). */
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

  /** @brief Clears the accumulated rotation back to identity. */
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
  float delta_lambda = 0;
  float delta_phi = 0;
  float delta_gamma = 0;
  float cos_dp = cosf(delta_phi);
  float sin_dp = sinf(delta_phi);
  float cos_dg = cosf(delta_gamma);
  float sin_dg = sinf(delta_gamma);
};
