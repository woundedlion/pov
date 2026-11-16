/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Orientation test cases. */

#include "tests/transformers/fixture.h"

namespace hs_test::transformers_tests {

/**
 * @brief Verifies an identity orientation leaves every sampled direction
 *        unchanged.
 */
inline void test_orient_transformer_identity() {
  math::Orientation<> ori;
  OrientTransformer ot(ori);

  const math::Vector samples[] = {math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                                  math::Vector(0, 0, 1),
                                  math::Vector(0.5f, 0.1f, 0.3f).normalized()};
  for (const math::Vector &v : samples) {
    math::Vector r = ot(v);
    HS_EXPECT_NEAR(r.x, v.x, 1e-5f);
    HS_EXPECT_NEAR(r.y, v.y, 1e-5f);
    HS_EXPECT_NEAR(r.z, v.z, 1e-5f);
  }
}

/**
 * @brief Verifies a non-identity orientation rotates by exactly the quaternion
 *        it was built from.
 * @details +90° about +y maps (x, y, z) → (z, y, -x) by the right-hand rule,
 *          independent of the quaternion code under test.
 */
inline void test_orient_transformer_known_rotation() {
  math::Orientation<> ori(math::make_rotation(
      math::Vector(0, 1, 0), math::PI_F * 0.5f)); // +90° about y
  OrientTransformer ot(ori);

  const math::Vector samples[] = {math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                                  math::Vector(0, 0, 1),
                                  math::Vector(0.5f, 0.1f, 0.3f).normalized()};
  for (const math::Vector &v : samples) {
    math::Vector r = ot(v);
    HS_EXPECT_NEAR(r.x, v.z, 1e-5f);
    HS_EXPECT_NEAR(r.y, v.y, 1e-5f);
    HS_EXPECT_NEAR(r.z, -v.x, 1e-5f);
  }
}

// ============================================================================
// mobius_transform / gnomonic_mobius_transform — rotations and reference parity
// ============================================================================

/**
 * @brief Verifies the identity Mobius map round-trips a point through
 *        stereographic projection back to itself, staying on the unit sphere.
 */
inline void test_mobius_identity_roundtrip() {
  math::MobiusParams id;
  math::Vector v = math::Vector(0.5f, 0.1f, 0.3f).normalized();
  math::Vector r = math::mobius_transform(v, id);
  HS_EXPECT_TRUE(finite_vec(r));
  HS_EXPECT_NEAR(r.x, v.x, 1e-5f);
  HS_EXPECT_NEAR(r.y, v.y, 1e-5f);
  HS_EXPECT_NEAR(r.z, v.z, 1e-5f);
  HS_EXPECT_NEAR(r.length(), 1.0f, 1e-3f);
}

/**
 * @brief Verifies a non-identity Mobius map produces its hand-computed image.
 * @details Under the stereographic convention, f(z) = 1/z is a 180° rotation
 *          about the x-axis.
 */
inline void test_mobius_known_rotation() {
  math::MobiusParams inv(0, 0, 1, 0, 1, 0, 0, 0); // f(z) = 1/z
  math::Vector v = math::Vector(0.5f, 0.1f, 0.3f).normalized();
  math::Vector r = math::mobius_transform(v, inv);
  HS_EXPECT_TRUE(finite_vec(r));
  HS_EXPECT_NEAR(r.x, v.x, 1e-3f);
  HS_EXPECT_NEAR(r.y, -v.y, 1e-3f);
  HS_EXPECT_NEAR(r.z, -v.z, 1e-3f);
  HS_EXPECT_NEAR(r.length(), 1.0f, 1e-3f);
}

/**
 * @brief Verifies mobius_transform matches a double-precision projection ->
 *        map -> unprojection oracle over random points and coefficients, and
 *        lands exactly on the sphere.
 * @details Pole-cap points, where the oracle is ill-conditioned, skip the
 *          value comparison but must still map to finite unit vectors.
 */
inline void test_mobius_matches_double_precision_oracle() {
  hs::random().seed(20260720);
  auto oracle = [](const math::Vector &v, const math::MobiusParams &p) {
    double s = 1.0 - static_cast<double>(v.y);
    double zr = static_cast<double>(v.x) / s;
    double zi = static_cast<double>(v.z) / s;
    double nr = static_cast<double>(p.a.re) * zr -
                static_cast<double>(p.a.im) * zi + static_cast<double>(p.b.re);
    double ni = static_cast<double>(p.a.re) * zi +
                static_cast<double>(p.a.im) * zr + static_cast<double>(p.b.im);
    double dr = static_cast<double>(p.c.re) * zr -
                static_cast<double>(p.c.im) * zi + static_cast<double>(p.d.re);
    double di = static_cast<double>(p.c.re) * zi +
                static_cast<double>(p.c.im) * zr + static_cast<double>(p.d.im);
    double q = dr * dr + di * di;
    double wr = (nr * dr + ni * di) / q;
    double wi = (ni * dr - nr * di) / q;
    double r2 = wr * wr + wi * wi;
    return math::Vector(static_cast<float>(2.0 * wr / (r2 + 1.0)),
                        static_cast<float>((r2 - 1.0) / (r2 + 1.0)),
                        static_cast<float>(2.0 * wi / (r2 + 1.0)));
  };

  int compared = 0;
  int nonfinite = 0;
  float worst_unit_err = 0.0f;
  float worst_oracle_err = 0.0f;
  for (int n = 0; n < 4000; ++n) {
    const math::Vector v = rand_unit();
    const float ar = hs::rand_f(-2, 2);
    const float ai = hs::rand_f(-2, 2);
    const float br = hs::rand_f(-2, 2);
    const float bi = hs::rand_f(-2, 2);
    const float cr = hs::rand_f(-2, 2);
    const float ci = hs::rand_f(-2, 2);
    const float dr = hs::rand_f(-2, 2);
    const float di = hs::rand_f(-2, 2);
    math::MobiusParams p(ar, ai, br, bi, cr, ci, dr, di);
    // Skip a near-singular draw: ad - bc ~ 0 collapses the map and both forms
    // are then dominated by cancellation, not by the formulation.
    float det_re =
        p.a.re * p.d.re - p.a.im * p.d.im - p.b.re * p.c.re + p.b.im * p.c.im;
    float det_im =
        p.a.re * p.d.im + p.a.im * p.d.re - p.b.re * p.c.im - p.b.im * p.c.re;
    if (det_re * det_re + det_im * det_im < 0.25f)
      continue;

    math::Vector got = math::mobius_transform(v, p);
    nonfinite += !finite_vec(got);
    worst_unit_err = std::max(worst_unit_err, std::fabs(got.length() - 1.0f));
    if (1.0f - v.y < projections::STEREO_POLE_EPS)
      continue;
    math::Vector want = oracle(v, p);
    ++compared;
    worst_oracle_err = std::max(worst_oracle_err, std::fabs(got.x - want.x));
    worst_oracle_err = std::max(worst_oracle_err, std::fabs(got.y - want.y));
    worst_oracle_err = std::max(worst_oracle_err, std::fabs(got.z - want.z));
  }
  HS_EXPECT_EQ(nonfinite, 0);
  HS_EXPECT_LE(worst_unit_err, 1e-5f);
  HS_EXPECT_LE(worst_oracle_err, 1e-5f);
  HS_EXPECT_GT(compared, 1000);
}

/**
 * @brief Verifies the poles map to the coefficient ratios the projective
 *        extension prescribes: the north pole (the point at infinity) to a/c
 *        and the south pole (the origin) to b/d.
 */
inline void test_mobius_poles_map_to_coefficient_ratios() {
  math::MobiusParams p(0.7f, 0.2f, -0.4f, 0.9f, 0.3f, -0.6f, 1.1f, 0.5f);
  math::Vector north = math::mobius_transform(math::Vector(0, 1, 0), p);
  HS_EXPECT_TRUE(finite_vec(north));
  math::Vector north_want = projections::inv_stereo(p.a / p.c);
  HS_EXPECT_NEAR(north.x, north_want.x, 1e-4f);
  HS_EXPECT_NEAR(north.y, north_want.y, 1e-4f);
  HS_EXPECT_NEAR(north.z, north_want.z, 1e-4f);

  math::Vector south = math::mobius_transform(math::Vector(0, -1, 0), p);
  HS_EXPECT_TRUE(finite_vec(south));
  math::Vector south_want = projections::inv_stereo(p.b / p.d);
  HS_EXPECT_NEAR(south.x, south_want.x, 1e-4f);
  HS_EXPECT_NEAR(south.y, south_want.y, 1e-4f);
  HS_EXPECT_NEAR(south.z, south_want.z, 1e-4f);
}

/**
 * @brief Verifies the identity Mobius map round-trips a point through
 *        gnomonic projection back to itself, including across the equatorial
 *        singularity where the projection saturates to the infinity sentinel.
 */
inline void test_gnomonic_mobius_identity_roundtrip() {
  math::MobiusParams id;
  math::Vector v = math::Vector(0.3f, 0.7f, 0.2f).normalized();
  math::Vector r = math::gnomonic_mobius_transform(v, id);
  HS_EXPECT_TRUE(finite_vec(r));
  HS_EXPECT_NEAR(r.x, v.x, 1e-5f);
  HS_EXPECT_NEAR(r.y, v.y, 1e-5f);
  HS_EXPECT_NEAR(r.z, v.z, 1e-5f);

  for (float y : {1e-3f, 2.1e-4f, 1.9e-4f, 1e-5f, 1e-9f, 0.0f, -1e-9f, -1e-5f,
                  -1.9e-4f, -2.1e-4f, -1e-3f}) {
    for (float theta = 0.0f; theta < 2.0f * math::PI_F; theta += 0.31f) {
      const math::Vector near_equator =
          math::Vector(cosf(theta), y, sinf(theta)).normalized();
      const math::Vector back =
          math::gnomonic_mobius_transform(near_equator, id);
      HS_EXPECT_TRUE(finite_vec(back));
      HS_EXPECT_NEAR(back.x, near_equator.x, 1e-5f);
      const float snap_y =
          1.0f / sqrtf(1.0f + projections::STEREO_INF_RECOGNIZE *
                                  projections::STEREO_INF_RECOGNIZE);
      const float expected_y =
          std::abs(near_equator.y) <= snap_y ? 0.0f : near_equator.y;
      HS_EXPECT_NEAR(back.y, expected_y, 1e-5f);
      HS_EXPECT_NEAR(back.z, near_equator.z, 1e-5f);
    }
  }
}

/**
 * @brief Verifies a non-identity gnomonic Mobius map produces its hand-computed
 *        image.
 * @details Under the gnomonic projection, f(z) = -z is a 180° rotation about
 *          the y-axis.
 */
inline void test_gnomonic_mobius_known_rotation() {
  math::MobiusParams neg(-1, 0, 0, 0, 0, 0, 1, 0); // f(z) = -z
  for (float y : {0.7f, -0.7f}) {
    const math::Vector v = math::Vector(0.3f, y, 0.2f).normalized();
    const math::Vector r = math::gnomonic_mobius_transform(v, neg);
    HS_EXPECT_TRUE(finite_vec(r));
    HS_EXPECT_NEAR(r.x, -v.x, 1e-3f);
    HS_EXPECT_NEAR(r.y, v.y, 1e-3f);
    HS_EXPECT_NEAR(r.z, -v.z, 1e-3f);
  }
}

/**
 * @brief Verifies a rotation with a complex divisor carries an equator point
 *        to the rotated equator longitude, independent of coefficient scaling.
 */
inline void test_gnomonic_mobius_equator_rotation() {
  const math::Vector V(1.0f, 0.0f, 0.0f);
  math::MobiusParams times_i_by_divisor(1, 0, 0, 0, 0, 0, 0, -1); // z / -i
  math::MobiusParams times_i(0, 1, 0, 0, 0, 0, 1, 0);             // i z / 1
  for (const math::MobiusParams &p : {times_i_by_divisor, times_i}) {
    const math::Vector r = math::gnomonic_mobius_transform(V, p);
    HS_EXPECT_NEAR(r.x, 0.0f, 1e-3f);
    HS_EXPECT_NEAR(r.y, 0.0f, 1e-3f);
    HS_EXPECT_NEAR(r.z, 1.0f, 1e-3f);
  }
}

/**
 * @brief Verifies an equator point round-trips through the identity map for
 *        either signed zero.
 * @details The hemisphere sign keys on the sign bit, so y == -0.0f must not
 *          reflect the result.
 */
inline void test_gnomonic_mobius_signed_zero_equator() {
  math::MobiusParams id;
  for (float y : {-0.0f, 0.0f, -1e-12f, 1e-12f}) {
    const math::Vector v(1.0f, y, 0.0f);
    const math::Vector r = math::gnomonic_mobius_transform(v, id);
    HS_EXPECT_NEAR(r.x, 1.0f, 1e-3f);
    HS_EXPECT_NEAR(r.y, 0.0f, 1e-3f);
    HS_EXPECT_NEAR(r.z, 0.0f, 1e-3f);
  }
}

// ============================================================================
// ripple_transform
// ============================================================================

/**
 * @brief Verifies zero amplitude short-circuits the ripple to a no-op.
 */
inline void test_ripple_zero_amplitude_is_identity() {
  Animation::RippleParams p;
  p.center = math::Vector(0, 1, 0);
  p.amplitude = 0.0f;
  p.phase = 0.5f;
  math::Vector v = math::Vector(1, 0, 0);
  math::Vector r = ripple_transform(v, p);
  HS_EXPECT_NEAR(r.x, v.x, 1e-6f);
  HS_EXPECT_NEAR(r.y, v.y, 1e-6f);
  HS_EXPECT_NEAR(r.z, v.z, 1e-6f);
}

} // namespace hs_test::transformers_tests
