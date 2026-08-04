/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/engine/transformers.h — the pure geometry transform functions
 * and the adapter/manager wrappers:
 *   - OrientTransformer       : identity orientation is a no-op; a known 90°
 *                               rotation maps to its hand-computed image.
 *   - mobius_transform        : identity Mobius round-trips through stereo; the
 *                               1/z map realizes a 180° rotation about x; random
 *                               maps match a double-precision oracle and stay
 *                               unit; the poles map to a/c and b/d.
 *   - gnomonic_mobius_transform: identity round-trips through gnomonic, right
 *                               through the equatorial singularity; the -z map
 *                               realizes a 180° rotation about y; an equator
 *                               point round-trips for either signed zero.
 *   - ripple_transform        : amplitude 0 and center-point degeneracies are
 *                               no-ops; an active ripple rotates on-sphere; the
 *                               prepared-threshold fast-reject band applies
 *                               in-band points and rejects off-band ones; the
 *                               truncated small-angle quaternion series matches
 *                               the exact libm rotation and meets the exact
 *                               branch at the switchover.
 *   - noise_transform         : amplitude≈0 is a no-op; otherwise stays unit; at
 *                               a shipped preset's amplitude the cross-hemisphere
 *                               slide cap binds the displacement at atan(0.5).
 *   - Transformer<>           : no active entities → identity; a spawned entity
 *                               applies and multiple entities compose; a recycled
 *                               freed slot composes in spawn order, not slot order.
 *   - FieldTransformer<>      : no active entities → 0; spawned entities sum;
 *                               field_bound sums the per-entity bounds;
 *                               active_params runs in spawn order; a completed
 *                               entity's slot is reclaimed.
 *   - bump_field              : drape push away from the cap center along the
 *                               polar direction — zero for the ring through the
 *                               center and at the footprint edge, peaking
 *                               between (antisymmetric, exactly 0 outside the
 *                               fast reject); the lifecycle envelope scales the
 *                               footprint; a parameter sweep keeps |field|
 *                               within field_bound().
 *   - DominantFieldAccumulator: nothing added reports 0; a lone field passes
 *                               through exactly and repeats neither stack nor
 *                               shift; the strongest of a mixed overlap
 *                               dominates without summing; an opposite-signed
 *                               crossing is continuous and lands on 0.
 *   - noise_product_field     : matches the hand-computed two-octave product;
 *                               ~0 amplitude short-circuits to exactly 0; a
 *                               parameter sweep keeps |field| within
 *                               field_bound() (both octaves stay in [-1, 1]).
 *   - Animation::BallDrop     : center traverses pole to pole along its azimuth
 *                               meridian; envelope 0 at the poles, 1 mid-fall;
 *                               the pool slot frees at completion.
 *   - Animation::NoiseProduct : integrates time by speed per step.
 */
#pragma once

#include "core/engine/transformers.h"
#include "core/render/canvas.h"
#include "core/math/easing.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cstdint>
#include <cstring>

namespace hs_test {
namespace transformers_tests {

/**
 * @brief Tests whether all three components of a vector are finite.
 * @param v Vector to inspect.
 * @return True when x, y, and z are all finite (no NaN/Inf).
 */
inline bool finite_vec(const Vector &v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

/**
 * @brief Tests whether two vectors are bit-for-bit identical.
 * @param a,b Vectors to compare.
 * @return True when every component shares the exact IEEE-754 bit pattern.
 * @details Used to prove an identity short-circuit returned its input VERBATIM,
 * not merely a (possibly different) non-finite value. A NaN component compares
 * unequal to itself under ==, so the raw bits are compared instead.
 */
inline bool vec_bits_equal(const Vector &a, const Vector &b) {
  auto bits = [](float f) {
    uint32_t u;
    std::memcpy(&u, &f, sizeof u);
    return u;
  };
  return bits(a.x) == bits(b.x) && bits(a.y) == bits(b.y) &&
         bits(a.z) == bits(b.z);
}

// ============================================================================
// OrientTransformer
// ============================================================================

/**
 * @brief Verifies an identity orientation leaves every sampled direction
 *        unchanged.
 */
inline void test_orient_transformer_identity() {
  Orientation<> ori;
  OrientTransformer ot(ori);

  const Vector samples[] = {Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1),
                            Vector(0.5f, 0.1f, 0.3f).normalized()};
  for (const Vector &v : samples) {
    Vector r = ot(v);
    HS_EXPECT_NEAR(r.x, v.x, 1e-5f);
    HS_EXPECT_NEAR(r.y, v.y, 1e-5f);
    HS_EXPECT_NEAR(r.z, v.z, 1e-5f);
  }
}

/**
 * @brief Verifies a non-identity orientation rotates by exactly the quaternion
 *        it was built from, catching a transform that ignored its rotation.
 * @details The orientation is a +90° rotation about the +y axis. By the
 *          right-hand rule that maps (x, y, z) → (z, y, -x) — an oracle computed
 *          independently of the quaternion-rotation code under test. The
 *          identity-only case (test_orient_transformer_identity) passes even for
 *          a no-op transform; this case fails unless the rotation is applied.
 */
inline void test_orient_transformer_known_rotation() {
  Orientation<> ori(
      make_rotation(Vector(0, 1, 0), PI_F * 0.5f)); // +90° about y
  OrientTransformer ot(ori);

  const Vector samples[] = {Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1),
                            Vector(0.5f, 0.1f, 0.3f).normalized()};
  for (const Vector &v : samples) {
    Vector r = ot(v);
    HS_EXPECT_NEAR(r.x, v.z, 1e-5f);
    HS_EXPECT_NEAR(r.y, v.y, 1e-5f);
    HS_EXPECT_NEAR(r.z, -v.x, 1e-5f);
  }
}

// ============================================================================
// mobius_transform / gnomonic_mobius_transform — identity round-trips
// ============================================================================

/**
 * @brief Verifies the identity Mobius map round-trips a point through
 *        stereographic projection back to itself, staying on the unit sphere.
 */
inline void test_mobius_identity_roundtrip() {
  MobiusParams id;
  Vector v = Vector(0.5f, 0.1f, 0.3f).normalized();
  Vector r = mobius_transform(v, id);
  HS_EXPECT_TRUE(finite_vec(r));
  HS_EXPECT_NEAR(r.x, v.x, 2e-3f);
  HS_EXPECT_NEAR(r.y, v.y, 2e-3f);
  HS_EXPECT_NEAR(r.z, v.z, 2e-3f);
  HS_EXPECT_NEAR(r.length(), 1.0f, 1e-3f);
}

/**
 * @brief Verifies a non-identity Mobius map produces its hand-computed image,
 *        catching a transform that ignored its coefficients.
 * @details The map f(z) = 1/z (a=0, b=1, c=1, d=0) is, under this file's
 *          stereographic convention (origin↔south pole, ∞↔north pole), a 180°
 *          rotation of the sphere about the x-axis: (x, y, z) → (x, -y, -z).
 *          That image is derived analytically, not from the stereo/mobius code
 *          under test, so a coefficient-ignoring (identity) implementation —
 *          which the existing round-trip case cannot distinguish — fails here.
 */
inline void test_mobius_known_rotation() {
  MobiusParams inv(0, 0, 1, 0, 1, 0, 0, 0); // f(z) = 1/z
  Vector v = Vector(0.5f, 0.1f, 0.3f).normalized();
  Vector r = mobius_transform(v, inv);
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
 * @details The oracle divides through in double before applying the map, so it
 *          is independent of the homogeneous formulation under test. Points
 *          inside the pole cap are excluded from the value comparison — there
 *          the oracle's own quotient is ill-conditioned — but their image is
 *          still required to be finite and unit. The 4000 draws report through
 *          worst-error accumulators rather than per-draw assertions: per-draw
 *          HS_EXPECTs would put ~20k assertions into the module's floor and
 *          leave the floor gate blind to every other case being deleted.
 */
inline void test_mobius_matches_double_precision_oracle() {
  hs::random().seed(20260720);
  auto oracle = [](const Vector &v, const MobiusParams &p) {
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
    return Vector(static_cast<float>(2.0 * wr / (r2 + 1.0)),
                  static_cast<float>((r2 - 1.0) / (r2 + 1.0)),
                  static_cast<float>(2.0 * wi / (r2 + 1.0)));
  };

  int compared = 0;
  int nonfinite = 0;
  float worst_unit_err = 0.0f;
  float worst_oracle_err = 0.0f;
  for (int n = 0; n < 4000; ++n) {
    Vector v;
    for (;;) {
      const float rx = hs::rand_f(-1, 1);
      const float ry = hs::rand_f(-1, 1);
      const float rz = hs::rand_f(-1, 1);
      Vector r(rx, ry, rz);
      if (r.length() > 0.1f) {
        v = r.normalized();
        break;
      }
    }
    const float ar = hs::rand_f(-2, 2);
    const float ai = hs::rand_f(-2, 2);
    const float br = hs::rand_f(-2, 2);
    const float bi = hs::rand_f(-2, 2);
    const float cr = hs::rand_f(-2, 2);
    const float ci = hs::rand_f(-2, 2);
    const float dr = hs::rand_f(-2, 2);
    const float di = hs::rand_f(-2, 2);
    MobiusParams p(ar, ai, br, bi, cr, ci, dr, di);
    // Skip a near-singular draw: ad - bc ~ 0 collapses the map and both forms
    // are then dominated by cancellation, not by the formulation.
    float det_re =
        p.a.re * p.d.re - p.a.im * p.d.im - p.b.re * p.c.re + p.b.im * p.c.im;
    float det_im =
        p.a.re * p.d.im + p.a.im * p.d.re - p.b.re * p.c.im - p.b.im * p.c.re;
    if (det_re * det_re + det_im * det_im < 0.25f)
      continue;

    Vector got = mobius_transform(v, p);
    nonfinite += !finite_vec(got);
    worst_unit_err = std::max(worst_unit_err, std::fabs(got.length() - 1.0f));
    if (1.0f - v.y < STEREO_POLE_EPS)
      continue;
    Vector want = oracle(v, p);
    ++compared;
    worst_oracle_err = std::max(worst_oracle_err, std::fabs(got.x - want.x));
    worst_oracle_err = std::max(worst_oracle_err, std::fabs(got.y - want.y));
    worst_oracle_err = std::max(worst_oracle_err, std::fabs(got.z - want.z));
  }
  HS_EXPECT_EQ(nonfinite, 0);
  HS_EXPECT_LE(worst_unit_err, 1e-5f);
  HS_EXPECT_LE(worst_oracle_err, 2e-3f);
  HS_EXPECT_GT(compared, 1000);
}

/**
 * @brief Verifies the poles map to the coefficient ratios the projective
 *        extension prescribes: the north pole (the point at infinity) to a/c
 *        and the south pole (the origin) to b/d.
 */
inline void test_mobius_poles_map_to_coefficient_ratios() {
  MobiusParams p(0.7f, 0.2f, -0.4f, 0.9f, 0.3f, -0.6f, 1.1f, 0.5f);
  Vector north = mobius_transform(Vector(0, 1, 0), p);
  HS_EXPECT_TRUE(finite_vec(north));
  Vector north_want = inv_stereo(p.a / p.c);
  HS_EXPECT_NEAR(north.x, north_want.x, 1e-4f);
  HS_EXPECT_NEAR(north.y, north_want.y, 1e-4f);
  HS_EXPECT_NEAR(north.z, north_want.z, 1e-4f);

  Vector south = mobius_transform(Vector(0, -1, 0), p);
  HS_EXPECT_TRUE(finite_vec(south));
  Vector south_want = inv_stereo(p.b / p.d);
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
  MobiusParams id;
  Vector v = Vector(0.3f, 0.7f, 0.2f).normalized();
  Vector r = gnomonic_mobius_transform(v, id);
  HS_EXPECT_TRUE(finite_vec(r));
  HS_EXPECT_NEAR(r.x, v.x, 2e-3f);
  HS_EXPECT_NEAR(r.y, v.y, 2e-3f);
  HS_EXPECT_NEAR(r.z, v.z, 2e-3f);

  for (float y :
       {1e-3f, 2e-4f, 1e-5f, 1e-9f, 0.0f, -1e-9f, -1e-5f, -2e-4f, -1e-3f}) {
    for (float theta = 0.0f; theta < 2.0f * PI_F; theta += 0.31f) {
      const Vector near_equator =
          Vector(cosf(theta), y, sinf(theta)).normalized();
      const Vector back = gnomonic_mobius_transform(near_equator, id);
      HS_EXPECT_TRUE(finite_vec(back));
      HS_EXPECT_NEAR(back.x, near_equator.x, 2e-3f);
      HS_EXPECT_NEAR(back.y, near_equator.y, 2e-3f);
      HS_EXPECT_NEAR(back.z, near_equator.z, 2e-3f);
    }
  }
}

/**
 * @brief Verifies a non-identity gnomonic Mobius map produces its hand-computed
 *        image, catching a transform that ignored its coefficients.
 * @details The map f(z) = -z (a=-1, b=0, c=0, d=1) is, under the gnomonic
 *          projection (tangent at the north pole, hemisphere restored from the
 *          sign of v.y), a 180° rotation about the y-axis: for an upper-
 *          hemisphere point, (x, y, z) → (-x, y, -z). The image is derived
 *          analytically (inv_len = v.y collapses the projection algebra), so an
 *          identity implementation that the round-trip case admits fails here.
 */
inline void test_gnomonic_mobius_known_rotation() {
  MobiusParams neg(-1, 0, 0, 0, 0, 0, 1, 0); // f(z) = -z
  for (float y : {0.7f, -0.7f}) {
    const Vector v = Vector(0.3f, y, 0.2f).normalized();
    const Vector r = gnomonic_mobius_transform(v, neg);
    HS_EXPECT_TRUE(finite_vec(r));
    HS_EXPECT_NEAR(r.x, -v.x, 1e-3f);
    HS_EXPECT_NEAR(r.y, v.y, 1e-3f);
    HS_EXPECT_NEAR(r.z, -v.z, 1e-3f);
  }
}

/**
 * @brief Verifies an equator point round-trips through the identity map for
 *        either signed zero.
 * @details The hemisphere sign has to key on the sign bit: a >= 0 test sends
 *          y == -0.0f to the divisor sign of +0.0f but the hemisphere sign of a
 *          tiny negative, which reflects the result 180° instead of returning
 *          it. The equator is its own antipode's neighbour, so the flip only
 *          shows as a broken round-trip.
 */
inline void test_gnomonic_mobius_signed_zero_equator() {
  MobiusParams id;
  for (float y : {-0.0f, 0.0f, -1e-12f, 1e-12f}) {
    const Vector v(1.0f, y, 0.0f);
    const Vector r = gnomonic_mobius_transform(v, id);
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
  p.center = Vector(0, 1, 0);
  p.amplitude = 0.0f;
  p.phase = 0.5f;
  Vector v = Vector(1, 0, 0);
  Vector r = ripple_transform(v, p);
  HS_EXPECT_NEAR(r.x, v.x, 1e-6f);
  HS_EXPECT_NEAR(r.y, v.y, 1e-6f);
  HS_EXPECT_NEAR(r.z, v.z, 1e-6f);
}

/**
 * @brief Verifies a point coincident with the ripple center is returned
 *        unchanged.
 * @details Such a point has a degenerate rotation axis (cross(center, center)
 *          == 0), so the transform must short-circuit to the identity.
 */
inline void test_ripple_center_point_is_identity() {
  Animation::RippleParams p;
  p.center = Vector(0, 1, 0);
  p.amplitude = 0.8f;
  p.phase = 0.0f;
  p.decay = 0.0f;
  Vector v = Vector(0, 1, 0);
  Vector r = ripple_transform(v, p);
  HS_EXPECT_NEAR(r.x, v.x, 1e-6f);
  HS_EXPECT_NEAR(r.y, v.y, 1e-6f);
  HS_EXPECT_NEAR(r.z, v.z, 1e-6f);
}

/**
 * @brief Verifies an active ripple at the wavelet peak rotates the point a
 *        noticeable amount while keeping it on the unit sphere.
 */
inline void test_ripple_active_rotates_on_sphere() {
  Animation::RippleParams p;
  p.center = Vector(0, 1, 0);
  p.amplitude = 0.5f;
  p.phase = PI_F * 0.5f; // wavelet peak at d == 90°
  p.decay = 0.0f;
  p.thickness = 1.0f;
  // Default thresholds (min=1, max=-1) disable the fast reject.

  Vector v = Vector(1, 0, 0); // 90° from center → at the wavelet peak
  Vector r = ripple_transform(v, p);

  HS_EXPECT_TRUE(finite_vec(r));
  HS_EXPECT_NEAR(r.length(), 1.0f, 1e-4f);
  float moved = std::abs(r.x - v.x) + std::abs(r.y - v.y) + std::abs(r.z - v.z);
  HS_EXPECT_GT(moved, 1e-2f);
}

/**
 * @brief Exercises the sync() fast-reject band the production
 *        renderer relies on (prepare_frame() calls it per active entity).
 * @details test_ripple_active_rotates_on_sphere deliberately leaves the default
 *          degenerate bounds (min=1, max=-1) in place, so the reject `if` is
 *          never taken there. Here real bounds are computed: a point at the
 *          wavelet peak (inside the band) is rotated, while points on either
 *          side — closer to the center than d_min, and farther than d_max —
 *          take the two reject legs and return unchanged. A sync()
 *          that swapped or mis-derived its bounds would either reject the peak
 *          (no move) or pass the off-band points (they move), failing here.
 */
inline void test_ripple_threshold_reject_path() {
  Animation::RippleParams p;
  p.center = Vector(0, 1, 0); // north pole, so cos_d == dot(v, center) == v.y
  p.amplitude = 0.5f;
  p.phase = PI_F * 0.5f; // wavelet peak at d == 90°
  p.decay = 0.0f;
  p.thickness = 0.4f; // band ≈ [phase - 0.4, phase + 0.4] rad
  p.sync();

  HS_EXPECT_LT(p.cos_threshold_min, 1.0f);
  HS_EXPECT_GT(p.cos_threshold_max, -1.0f);

  // In-band: at the peak.
  const Vector in_band = Vector(1, 0, 0);
  const Vector r_in = ripple_transform(in_band, p);
  HS_EXPECT_TRUE(finite_vec(r_in));
  HS_EXPECT_NEAR(r_in.length(), 1.0f, 1e-4f);
  const float moved_in = std::abs(r_in.x - in_band.x) +
                         std::abs(r_in.y - in_band.y) +
                         std::abs(r_in.z - in_band.z);
  HS_EXPECT_GT(moved_in, 1e-2f);

  // Out-of-band toward the center (d < d_min): cos_d > cos_threshold_min leg.
  const float d_near = p.phase - p.thickness - 0.3f;
  const Vector near_c = Vector(std::sin(d_near), std::cos(d_near), 0.0f);
  const Vector r_near = ripple_transform(near_c, p);
  HS_EXPECT_NEAR(r_near.x, near_c.x, 1e-6f);
  HS_EXPECT_NEAR(r_near.y, near_c.y, 1e-6f);
  HS_EXPECT_NEAR(r_near.z, near_c.z, 1e-6f);

  // Out-of-band away from the center (d > d_max): cos_d < cos_threshold_max leg.
  const float d_far = p.phase + p.thickness + 0.3f;
  const Vector far_c = Vector(std::sin(d_far), std::cos(d_far), 0.0f);
  const Vector r_far = ripple_transform(far_c, p);
  HS_EXPECT_NEAR(r_far.x, far_c.x, 1e-6f);
  HS_EXPECT_NEAR(r_far.y, far_c.y, 1e-6f);
  HS_EXPECT_NEAR(r_far.z, far_c.z, 1e-6f);
}

/**
 * @brief Verifies the decay parameter attenuates the ripple's rotation with
 *        angular distance from the center.
 * @details Every other ripple test pins decay to 0 (max strength), so the
 *          attenuation term expf(-decay * d) ships unexercised. Here the same
 *          peak point (d == 90° from the center) is transformed at several decay
 *          rates: a positive decay must move it strictly less than decay 0, a
 *          larger decay less still, and a strong decay collapses the rotation to
 *          ~identity — a decay term that was ignored (or sign-flipped) fails.
 */
inline void test_ripple_decay_attenuates() {
  Animation::RippleParams base;
  base.center = Vector(0, 1, 0);
  base.amplitude = 0.5f;
  base.phase = PI_F * 0.5f; // wavelet peak at d == 90°
  base.thickness = 1.0f;
  // Default thresholds (min=1, max=-1) disable the fast reject.
  const Vector v = Vector(1, 0, 0);

  auto moved = [&](float decay) {
    Animation::RippleParams p = base;
    p.decay = decay;
    Vector r = ripple_transform(v, p);
    return std::abs(r.x - v.x) + std::abs(r.y - v.y) + std::abs(r.z - v.z);
  };

  const float m0 = moved(0.0f);
  const float m_small = moved(0.5f);
  const float m_large = moved(8.0f);

  HS_EXPECT_GT(m0, 1e-2f);
  HS_EXPECT_LT(m_small, m0);
  HS_EXPECT_GT(m_small, 0.0f);
  HS_EXPECT_LT(m_large, m_small);
  HS_EXPECT_NEAR(m_large, 0.0f, 5e-2f);
}

/**
 * @brief Pins the sync() reject band at its edges, not just well
 *        outside them.
 * @details test_ripple_threshold_reject_path samples points 0.3 rad beyond each
 *          edge, so a band misplaced by up to ~0.3 rad would still pass. Here a
 *          point a hair inside each edge takes the slow path (the wavelet tail is
 *          tiny but nonzero, so it moves), while a point a hair outside is
 *          fast-rejected and returned bit-for-bit unchanged. The asymmetry
 *          (moves vs. exactly-equal) brackets the edge to within the epsilon.
 */
inline void test_ripple_threshold_boundary() {
  Animation::RippleParams p;
  p.center = Vector(0, 1, 0); // north pole → cos_d == dot(v, center) == v.y
  p.amplitude = 0.5f;
  p.phase = PI_F * 0.5f;
  p.decay = 0.0f;
  p.thickness = 0.4f; // band edges at phase ± thickness
  p.sync();

  const float d_min = p.phase - p.thickness;
  const float d_max = p.phase + p.thickness;
  const float eps = 0.02f;

  auto pt = [](float d) { return Vector(std::sin(d), std::cos(d), 0.0f); };
  auto moved = [&](const Vector &src) {
    Vector r = ripple_transform(src, p);
    return std::abs(r.x - src.x) + std::abs(r.y - src.y) +
           std::abs(r.z - src.z);
  };

  // Just inside each edge: wavelet tail is small but nonzero.
  HS_EXPECT_GT(moved(pt(d_max - eps)), 1e-4f);
  HS_EXPECT_GT(moved(pt(d_min + eps)), 1e-4f);
  // Just outside each edge: fast-rejected → returned exactly unchanged.
  HS_EXPECT_NEAR(moved(pt(d_max + eps)), 0.0f, 1e-7f);
  HS_EXPECT_NEAR(moved(pt(d_min - eps)), 0.0f, 1e-7f);
}

/**
 * @brief Verifies the truncated small-angle quaternion series reproduces the
 *        exact libm rotation it stands in for, and that the two branches meet
 *        at RIPPLE_SMALL_ANGLE_MAX.
 * @details The |theta| <= RIPPLE_SMALL_ANGLE_MAX branch replaces
 *          make_rotation's sin/cos with three-term series, so nothing but a
 *          comparison against the exact form pins its coefficients — the one
 *          test that crosses the boundary asserts only monotonicity, which a
 *          wrong coefficient survives while seaming visibly on device. theta is
 *          exactly proportional to amplitude for a fixed geometry, so one
 *          reading taken in the exact branch calibrates the envelope (the
 *          approximate acos/expf that produce it are common to both branches
 *          and cancel), and every series-branch amplitude is then judged
 *          against rotate(v, make_rotation(axis, theta)) — the closed form the
 *          series stands in for.
 */
inline void test_ripple_small_angle_series_matches_exact() {
  Animation::RippleParams p;
  p.center = Vector(0, 1, 0);
  p.phase = PI_F * 0.5f; // wavelet peak at d == 90 degrees
  p.decay = 0.0f;
  p.thickness = 1.0f;
  // Default thresholds (min=1, max=-1) disable the fast reject.
  const Vector v = Vector(1, 0, 0);
  const Vector axis = cross(p.center, v).normalized();

  auto at = [&](float amplitude) {
    Animation::RippleParams q = p;
    q.amplitude = amplitude;
    return ripple_transform(v, q);
  };
  auto exact = [&](float theta) {
    return rotate(v, make_rotation(axis, theta));
  };
  auto err = [](const Vector &a, const Vector &b) {
    return std::max(std::max(std::fabs(a.x - b.x), std::fabs(a.y - b.y)),
                    std::fabs(a.z - b.z));
  };

  // Calibrate theta/amplitude from the exact branch, in libm acos rather than
  // the fast_acos angle_between uses.
  auto turned = [&](const Vector &out) {
    return (float)std::acos(hs::clamp(dot(v, out), -1.0f, 1.0f));
  };
  const float anchor_amplitude = 0.9f;
  const Vector anchor = at(anchor_amplitude);
  const float envelope = turned(anchor) / anchor_amplitude;
  // The peak of a zero-decay Ricker wavelet is unit gain.
  HS_EXPECT_NEAR(envelope, 1.0f, 2e-3f);
  HS_EXPECT_LE(err(anchor, exact(anchor_amplitude * envelope)), 1e-6f);

  // Series branch, from just above the amplitude short-circuit up to the
  // switchover.
  const float amp_lo = 0.002f;
  const float amp_hi = RIPPLE_SMALL_ANGLE_MAX / envelope;
  float worst_err = 0.0f;
  float worst_unit_err = 0.0f;
  const int steps = 300;
  for (int i = 0; i <= steps; ++i) {
    const float amplitude =
        amp_lo + (amp_hi - amp_lo) * (float)i / (float)steps;
    const Vector got = at(amplitude);
    worst_err = std::max(worst_err, err(got, exact(amplitude * envelope)));
    worst_unit_err = std::max(worst_unit_err, std::fabs(got.length() - 1.0f));
  }
  HS_EXPECT_LE(worst_err, 5e-7f);
  HS_EXPECT_LE(worst_unit_err, 5e-7f);

  // The two branches meet: a hair either side of the switchover agrees.
  const float seam = RIPPLE_SMALL_ANGLE_MAX / envelope;
  HS_EXPECT_LE(err(at(seam * 0.999f), at(seam * 1.001f)), 1e-3f);
}

/**
 * @brief Pins the sync() band past π, where every ripple ends its life.
 * @details A Ripple advances phase by `speed` per frame for `duration` frames,
 *          so its defaults (0.2, 100) walk it out to ~20 rad — many times the
 *          sphere's π of angular range. Every other threshold test pins phase
 *          to π/2, where the [0, π] clamp is inert. Here the clamp collapses
 *          both bounds onto cos(π), so the front has left the sphere and every
 *          direction is fast-rejected; an inverted or dropped clamp reopens an
 *          arbitrary annulus of a wrapped phase and charges the whole sphere
 *          the per-pixel wavelet again.
 */
inline void test_ripple_threshold_collapses_past_pi() {
  Animation::RippleParams p;
  p.center = Vector(0, 1, 0);
  p.amplitude = 0.5f;
  p.decay = 0.0f;
  p.thickness = 0.4f;
  p.phase = 20.0f;
  p.sync();

  HS_EXPECT_NEAR(p.cos_threshold_min, -1.0f, 1e-6f);
  HS_EXPECT_NEAR(p.cos_threshold_max, -1.0f, 1e-6f);

  // The band admits nothing: every direction but the exact antipode falls
  // outside it, takes the fast reject, and comes back bit for bit.
  for (int i = 1; i < 64; ++i) {
    const float d = PI_F * static_cast<float>(i) / 64.0f;
    const Vector v(std::sin(d), std::cos(d), 0.0f);
    const float cos_d = dot(v, p.center);
    HS_EXPECT_TRUE(cos_d > p.cos_threshold_min || cos_d < p.cos_threshold_max);
    HS_EXPECT_TRUE(vec_bits_equal(ripple_transform(v, p), v));
  }
}

/**
 * @brief Feeds non-finite (NaN/Inf) directions through the transforms' identity
 *        short-circuits and confirms they pass the input through verbatim.
 * @details No in-process test previously pushed NaN/Inf in. The zero-amplitude
 *          and no-active-entity short-circuits return the input before any
 *          arithmetic, so a non-finite input comes back non-finite — they never
 *          choke on it nor manufacture a spurious finite result. The ACTIVE paths
 *          are deliberately fail-fast on a non-finite direction (the trailing
 *          normalized()/make_rotation traps rather than propagate garbage into
 *          the geometry); that trap is asserted out-of-process by the death
 *          harness (case_noise_transform_nan), since an HS_CHECK aborts the
 *          process and cannot be caught in-line here.
 */
inline void test_transforms_nonfinite_passes_through_identity() {
  const float inf = HUGE_VALF;
  const float nan = NAN;
  const Vector bad[] = {Vector(nan, 0, 0), Vector(0, inf, 0),
                        Vector(inf, nan, -inf)};

  for (const Vector &v : bad) {
    Animation::RippleParams rp;
    rp.amplitude = 0.0f;
    HS_EXPECT_TRUE(vec_bits_equal(ripple_transform(v, rp), v));
    Animation::NoiseParams np;
    np.amplitude = 0.0f;
    HS_EXPECT_TRUE(vec_bits_equal(noise_transform(v, np), v));
    Timeline tl;
    RippleTransformer<4> rt(tl);
    HS_EXPECT_TRUE(vec_bits_equal(rt.transform(v), v));
  }
}

// ============================================================================
// noise_transform
// ============================================================================

/**
 * @brief Verifies zero amplitude short-circuits the noise warp to a no-op.
 */
inline void test_noise_zero_amplitude_is_identity() {
  Animation::NoiseParams p;
  p.amplitude = 0.0f;
  Vector v = Vector(0.2f, 0.5f, 0.84f).normalized();
  Vector r = noise_transform(v, p);
  HS_EXPECT_NEAR(r.x, v.x, 1e-6f);
  HS_EXPECT_NEAR(r.y, v.y, 1e-6f);
  HS_EXPECT_NEAR(r.z, v.z, 1e-6f);
}

/**
 * @brief Verifies an active noise warp keeps every sample finite and on the
 *        unit sphere, and actually displaces the points.
 * @details The displacement assertion is the companion to
 *          test_noise_zero_amplitude_is_identity: without it, a noise_transform
 *          that returned the input unchanged would still pass the finite /
 *          on-sphere checks. Displacement is summed across samples so a single
 *          noise zero-crossing at one point cannot make the test flaky.
 */
inline void test_noise_active_stays_on_sphere() {
  Animation::NoiseParams p;
  p.amplitude = 0.5f;
  p.scale = 4.0f;
  p.time = 1.0f;
  const Vector samples[] = {Vector(1, 0, 0), Vector(0, 1, 0),
                            Vector(0.4f, 0.6f, 0.7f).normalized()};
  float total_moved = 0.0f;
  for (const Vector &v : samples) {
    Vector r = noise_transform(v, p);
    HS_EXPECT_TRUE(finite_vec(r));
    HS_EXPECT_NEAR(r.length(), 1.0f, 1e-3f);
    total_moved +=
        std::abs(r.x - v.x) + std::abs(r.y - v.y) + std::abs(r.z - v.z);
  }
  HS_EXPECT_GT(total_moved, 1e-2f);
}

/**
 * @brief Verifies the cross-hemisphere soft cap at the amplitude a shipped
 *        preset actually drives.
 * @details Feedback::Style::LooseWormhole() feeds noise_transform an amplitude
 *          of 11.25 — two orders above the ~0.6 every other test uses, and high
 *          enough that the raw slide would otherwise reach ~0.97 and grab across
 *          the sphere. The cap is on the tangential offset, so with the offset
 *          perpendicular to a unit v the resulting angular displacement is
 *          exactly atan(|offset|): saturated samples land on atan(0.5) and no
 *          sample may exceed it. Asserting both sides pins the cap — a removed
 *          or widened cap overshoots, a cap that never engages undershoots. The
 *          low-amplitude control confirms the existing tests stay clear of it.
 */
inline void test_noise_cross_hemisphere_cap() {
  constexpr float LOOSE_WORMHOLE_AMPLITUDE = 11.25f;
  const float capped_angle = std::atan(0.5f);

  auto sweep = [](float amplitude) {
    Animation::NoiseParams p;
    p.amplitude = amplitude;
    p.frequency = 0.22519f; // LooseWormhole's remaining noise config
    p.scale = 10.1798f;
    p.speed = 0.0f;
    p.time = 0.0f;
    p.noise.SetSeed(4242);
    p.sync();
    hs::random().seed(20260803);
    float worst_angle = 0.0f;
    int off_sphere = 0;
    for (int i = 0; i < 4000; ++i) {
      const Vector r(hs::rand_f(-1.0f, 1.0f), hs::rand_f(-1.0f, 1.0f),
                     hs::rand_f(-1.0f, 1.0f));
      if (r.length() < 0.1f)
        continue;
      const Vector v = r.normalized();
      const Vector out = noise_transform(v, p);
      off_sphere += !finite_vec(out) || std::fabs(out.length() - 1.0f) > 1e-3f;
      worst_angle = std::max(worst_angle, angle_between(v, out));
    }
    HS_EXPECT_EQ(off_sphere, 0);
    return worst_angle;
  };

  const float worst_high = sweep(LOOSE_WORMHOLE_AMPLITUDE);
  HS_EXPECT_LE(worst_high, capped_angle + 1e-4f);
  HS_EXPECT_GE(worst_high, capped_angle - 1e-4f);

  // 0.6 is the ceiling of every other noise test; the cap must be slack there.
  const float worst_low = sweep(0.6f);
  HS_EXPECT_LT(worst_low, capped_angle * 0.5f);
}

// ============================================================================
// Transformer<> manager — no active entities is the identity
// ============================================================================

/**
 * @brief Verifies a manager with nothing spawned (all entity slots inactive)
 *        is the identity.
 */
inline void test_transformer_no_entities_is_identity() {
  Timeline tl;
  RippleTransformer<8> rt(tl);
  Vector v = Vector(0.6f, 0.2f, 0.77f).normalized();
  Vector r = rt.transform(v);
  HS_EXPECT_NEAR(r.x, v.x, 1e-6f);
  HS_EXPECT_NEAR(r.y, v.y, 1e-6f);
  HS_EXPECT_NEAR(r.z, v.z, 1e-6f);
}

// ============================================================================
// Transformer<> active-slot list — spawn applies; multiple entities compose
// ============================================================================

/**
 * @brief Verifies the compact active-slot list: a spawned entity is applied by
 *        transform(), and two active entities compose while staying on the unit
 *        sphere.
 * @details Noise is used (rather than Ripple) because its ctor leaves the
 *          copied amplitude intact, so the spawned entity displaces immediately
 *          — no Canvas / timeline stepping needed. The local Timeline's
 *          destructor clears the global event buffer on scope exit.
 */
inline void test_transformer_spawn_applies_and_composes() {
  Timeline tl;
  global_timeline_t = 0;
  NoiseTransformer<4> nt(tl);
  nt.init_storage(persistent_arena);
  nt.template_params.amplitude = 0.6f;
  nt.template_params.scale = 4.0f;
  nt.template_params.time = 1.0f;
  nt.template_params.sync();

  const Vector samples[] = {Vector(1, 0, 0), Vector(0, 0, 1),
                            Vector(0.4f, 0.6f, 0.7f).normalized()};

  for (const Vector &v : samples) {
    Vector r = nt.transform(v);
    HS_EXPECT_NEAR(r.x, v.x, 1e-6f);
    HS_EXPECT_NEAR(r.y, v.y, 1e-6f);
    HS_EXPECT_NEAR(r.z, v.z, 1e-6f);
  }

  HS_EXPECT_TRUE(nt.spawn_pinned(0) != nullptr);
  constexpr size_t N = sizeof(samples) / sizeof(samples[0]);
  Vector single[N];
  float total_moved = 0.0f;
  for (size_t i = 0; i < N; ++i) {
    const Vector &v = samples[i];
    single[i] = nt.transform(v);
    HS_EXPECT_TRUE(finite_vec(single[i]));
    HS_EXPECT_NEAR(single[i].length(), 1.0f, 1e-3f);
    total_moved += std::abs(single[i].x - v.x) + std::abs(single[i].y - v.y) +
                   std::abs(single[i].z - v.z);
  }
  HS_EXPECT_GT(total_moved, 1e-2f);

  HS_EXPECT_TRUE(nt.spawn_pinned(0) != nullptr);
  float total_divergence = 0.0f;
  for (size_t i = 0; i < N; ++i) {
    const Vector &v = samples[i];
    Vector r = nt.transform(v);
    HS_EXPECT_TRUE(finite_vec(r));
    HS_EXPECT_NEAR(r.length(), 1.0f, 1e-3f);
    total_divergence += std::abs(r.x - single[i].x) +
                        std::abs(r.y - single[i].y) +
                        std::abs(r.z - single[i].z);
  }
  HS_EXPECT_GT(total_divergence, 1e-3f);
}

// ============================================================================
// Transformer<> non-pinned slot recycling survives timeline compaction
// ============================================================================

/**
 * @brief Verifies a non-pinned spawned transform's pool slot is reclaimed after
 *        the timeline relocates its event during compaction.
 * @details spawn() (unlike spawn_pinned) registers a finite, non-repeating event
 *          whose then() callback frees the entity slot — but only if that callback
 *          survives the event being move_into-relocated by step()'s compaction.
 *          An earlier finite event is added first so that, when it completes, the
 *          spawned ripple's event shifts down into the vacated slot. The pool has
 *          CAPACITY 1, so a leaked slot would make the post-completion re-spawn
 *          return nullptr; its success proves the relocated event's callback
 *          reclaimed the slot.
 */
inline void test_transformer_nonpinned_slot_reclaimed_after_compaction() {
  Timeline tl;
  global_timeline_t = 0;
  hs_test::StubEffect fx(8, 8);
  Canvas cv(fx);

  RippleTransformer<1> rt(tl);
  rt.init_storage(persistent_arena);

  float dummy = 0.0f;
  tl.add(0, Animation::Transition(dummy, 1.0f, 2, ease_linear));

  Animation::Ripple *p = rt.spawn(0, Vector(0, 1, 0), 0.2f, 4);
  HS_EXPECT_TRUE(p != nullptr);
  // Slot occupied: with CAPACITY 1 a second spawn finds nothing free.
  HS_EXPECT_TRUE(rt.spawn(0, Vector(0, 1, 0), 0.2f, 4) == nullptr);

  // Step past the earlier event (relocates the ripple) and the ripple's own
  // completion, so the then() callback fires through the relocation.
  for (int i = 0; i < 12; ++i)
    tl.step(cv);

  Animation::Ripple *reclaimed = rt.spawn(0, Vector(0, 1, 0), 0.2f, 4);
  HS_EXPECT_TRUE(reclaimed != nullptr);
}

// ============================================================================
// Transformer<> slots are reclaimed by Timeline::clear()
// ============================================================================

/**
 * @brief Verifies Timeline::clear() frees the pool slots of the animations it
 *        destroys.
 * @details clear() destroys events outright, so their then() callbacks — the
 *          only other reclaim path — never run. Without the pool's clear hook
 *          the slots stay active forever and the effect silently stops warping
 *          after CAPACITY spawns.
 */
inline void test_transformer_slots_released_by_timeline_clear() {
  Timeline tl;
  global_timeline_t = 0;

  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);

  HS_EXPECT_TRUE(rt.spawn(0, Vector(0, 1, 0), 0.2f, 60) != nullptr);
  HS_EXPECT_TRUE(rt.spawn(0, Vector(1, 0, 0), 0.2f, 60) != nullptr);
  HS_EXPECT_EQ(rt.active_count(), 2);
  HS_EXPECT_TRUE(rt.spawn(0, Vector(0, 0, 1), 0.2f, 60) == nullptr);

  tl.clear();
  HS_EXPECT_EQ(rt.active_count(), 0);
  HS_EXPECT_TRUE(rt.spawn(0, Vector(0, 0, 1), 0.2f, 60) != nullptr);
  HS_EXPECT_TRUE(rt.spawn(0, Vector(0, 1, 0), 0.2f, 60) != nullptr);
  HS_EXPECT_EQ(rt.active_count(), 2);
}

// ============================================================================
// Transformer<> composition order follows spawn order across slot recycling
// ============================================================================

/**
 * @brief Params carrying only a spawn-order tag, for the composition-order test.
 */
struct OrderParams {
  int order = 0; /**< Composition tag written by TagAnim at spawn. */
};

/**
 * @brief Shifts a base-10 digit of the tag into x on each application, so the
 *        final x reads back the exact order the warps composed in.
 */
inline Vector order_transform(const Vector &v, const OrderParams &p) {
  return Vector(v.x * 10.0f + static_cast<float>(p.order), v.y, v.z);
}

/**
 * @brief Minimal animation that stamps its spawn-order tag into the params.
 * @details Finite and non-repeating so spawn() accepts it and a short duration
 *          frees its pool slot on completion.
 */
struct TagAnim : public Animation::AnimationBase<TagAnim> {
  TagAnim(OrderParams &params, int order, int duration)
      : Animation::AnimationBase<TagAnim>(duration, false) {
    params.order = order;
  }
};

/**
 * @brief Verifies a recycled freed slot composes its new warp in spawn order
 *        (last), not slot-index order (first).
 * @details A short-lived warp holds slot 0 while a long-lived one takes slot 1;
 *          stepping frees slot 0, then a third warp recycles it. Composition
 *          must follow spawn order (long-lived tag 1, then recycled tag 2), so x
 *          reads 12 — slot-index order would compose the recycled slot 0 first
 *          and read 21.
 */
inline void test_transformer_recycled_slot_composes_in_spawn_order() {
  Timeline tl;
  global_timeline_t = 0;
  hs_test::StubEffect fx(8, 8);
  Canvas cv(fx);

  Transformer<OrderParams, TagAnim, order_transform, 2> tr(tl);
  tr.init_storage(persistent_arena);

  HS_EXPECT_TRUE(tr.spawn(0, /*order=*/9, /*duration=*/2) != nullptr); // slot 0
  HS_EXPECT_TRUE(tr.spawn(0, /*order=*/1, /*duration=*/100) !=
                 nullptr); // slot 1

  // Step until the short warp completes and frees slot 0; the long one survives.
  for (int i = 0; i < 5; ++i)
    tl.step(cv);

  HS_EXPECT_TRUE(tr.spawn(0, /*order=*/2, /*duration=*/100) !=
                 nullptr); // recycles slot 0

  const Vector r = tr.transform(Vector(0, 0, 0));
  HS_EXPECT_NEAR(r.x, 12.0f, 1e-4f);
}

// ============================================================================
// FieldTransformer<> — scalar fields sum; bounds sum; slots recycle
// ============================================================================

/**
 * @brief Params carrying one scalar contribution, for the field tests.
 */
struct FieldTestParams {
  float value = 0.0f; /**< Contribution written by FieldTagAnim at spawn. */

  /**
   * @brief Upper bound on |field_test_field| for this entity.
   * @return |value| (the field is value * v.x with |v.x| <= 1).
   */
  float field_bound() const { return std::fabs(value); }
};

/**
 * @brief Scales the sample's x by the entity's value, so the sum is readable
 *        and the point argument is proven to reach the field function.
 */
inline float field_test_field(const Vector &v, const FieldTestParams &p) {
  return p.value * v.x;
}

/**
 * @brief Minimal animation that stamps its value into the params.
 * @details Finite and non-repeating so spawn() accepts it and a short duration
 *          frees its pool slot on completion.
 */
struct FieldTagAnim : public Animation::AnimationBase<FieldTagAnim> {
  FieldTagAnim(FieldTestParams &params, float value, int duration)
      : Animation::AnimationBase<FieldTagAnim>(duration, false) {
    params.value = value;
  }
};

using TestFieldTransformer =
    FieldTransformer<FieldTestParams, FieldTagAnim, field_test_field, 2>;

/**
 * @brief Verifies a field manager with nothing spawned evaluates to 0
 *        everywhere and reports a 0 bound.
 */
inline void test_field_transformer_no_entities_is_zero() {
  Timeline tl;
  TestFieldTransformer ft(tl);
  HS_EXPECT_NEAR(ft.field(Vector(1, 0, 0)), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(ft.field_bound(), 0.0f, 1e-6f);
  HS_EXPECT_TRUE(ft.active_count() == 0);
}

/**
 * @brief Verifies spawned entities superpose by summation, scaled by the
 *        sample point, and field_bound sums the per-entity bounds.
 */
inline void test_field_transformer_sums_and_bounds() {
  Timeline tl;
  global_timeline_t = 0;
  TestFieldTransformer ft(tl);
  ft.init_storage(persistent_arena);

  HS_EXPECT_TRUE(ft.spawn(0, 2.0f, 100) != nullptr);
  HS_EXPECT_NEAR(ft.field(Vector(1, 0, 0)), 2.0f, 1e-5f);
  HS_EXPECT_NEAR(ft.field(Vector(-0.5f, 0, 0)), -1.0f, 1e-5f);

  HS_EXPECT_TRUE(ft.spawn(0, -3.0f, 100) != nullptr);
  HS_EXPECT_TRUE(ft.active_count() == 2);
  HS_EXPECT_NEAR(ft.field(Vector(1, 0, 0)), -1.0f, 1e-5f);
  HS_EXPECT_NEAR(ft.operator()(Vector(1, 0, 0)), -1.0f, 1e-5f);
  HS_EXPECT_NEAR(ft.field_bound(), 5.0f, 1e-5f);
}

/**
 * @brief Verifies active_params exposes spawned entities in spawn order.
 */
inline void test_field_transformer_active_params_spawn_order() {
  Timeline tl;
  global_timeline_t = 0;
  TestFieldTransformer ft(tl);
  ft.init_storage(persistent_arena);

  HS_EXPECT_TRUE(ft.spawn(0, 2.0f, 100) != nullptr);
  HS_EXPECT_TRUE(ft.spawn(0, -3.0f, 100) != nullptr);
  HS_EXPECT_NEAR(ft.active_params(0).value, 2.0f, 1e-6f);
  HS_EXPECT_NEAR(ft.active_params(1).value, -3.0f, 1e-6f);
}

/**
 * @brief Verifies a completed field entity's pool slot is reclaimed and a
 *        fresh spawn reuses it.
 * @details The pool has CAPACITY 2; both slots are filled (one short-lived),
 *          so the post-completion spawn succeeding proves the then() reclaim
 *          fired, and the field dropping the completed contribution proves the
 *          active list no longer visits the freed slot.
 */
inline void test_field_transformer_slot_reclaimed() {
  Timeline tl;
  global_timeline_t = 0;
  hs_test::StubEffect fx(8, 8);
  Canvas cv(fx);
  TestFieldTransformer ft(tl);
  ft.init_storage(persistent_arena);

  HS_EXPECT_TRUE(ft.spawn(0, 2.0f, 2) != nullptr);
  HS_EXPECT_TRUE(ft.spawn(0, 5.0f, 100) != nullptr);
  HS_EXPECT_TRUE(ft.spawn(0, 9.0f, 100) == nullptr);

  for (int i = 0; i < 6; ++i)
    tl.step(cv);

  HS_EXPECT_TRUE(ft.active_count() == 1);
  HS_EXPECT_NEAR(ft.field(Vector(1, 0, 0)), 5.0f, 1e-5f);
  HS_EXPECT_TRUE(ft.spawn(0, 1.0f, 100) != nullptr);
  HS_EXPECT_NEAR(ft.field(Vector(1, 0, 0)), 6.0f, 1e-5f);
}

/**
 * @brief Verifies reclaim_storage() carries a live entity across an arena
 *        rewind: the re-claimed slots keep their addresses and contents.
 * @details Rewinds via set_offset to the pre-init_storage mark, standing in
 *          for the reset() a mesh carousel compaction performs, then replays
 *          the pool allocation as the carousel's after-reset callback does.
 */
inline void test_field_transformer_storage_survives_arena_rewind() {
  Timeline tl;
  global_timeline_t = 0;
  const size_t mark = persistent_arena.get_offset();
  TestFieldTransformer ft(tl);
  ft.init_storage(persistent_arena);
  HS_EXPECT_TRUE(ft.spawn(0, 2.0f, 100) != nullptr);

  persistent_arena.set_offset(mark);
  ft.reclaim_storage(persistent_arena);
  HS_EXPECT_TRUE(ft.active_count() == 1);
  HS_EXPECT_NEAR(ft.field(Vector(1, 0, 0)), 2.0f, 1e-5f);
}

// ============================================================================
// DominantFieldAccumulator — pass-through, no stacking, smooth sign crossing
// ============================================================================

/**
 * @brief Verifies the documented blend identities: nothing added reports 0, a
 *        lone field passes through exactly, and repeats of one value neither
 *        stack nor shift.
 * @details Pass-through is the property a swapped numerator/denominator power
 *          destroys — sum(s^2)/sum(s^3) reports 1/s for a lone field, which
 *          matches s only at s == 1. Repeats pin the "dominate without stacking"
 *          contract that motivates the accumulator over plain summation.
 */
inline void test_dominant_field_identities() {
  HS_EXPECT_NEAR(DominantFieldAccumulator().value(), 0.0f, 0.0f);

  const float singles[] = {0.7f, -0.7f, 4.5f, -0.013f, 1.0f};
  for (float s : singles) {
    DominantFieldAccumulator acc;
    acc.add(s);
    HS_EXPECT_NEAR(acc.value(), s, 1e-6f);
  }

  for (float s : {0.7f, -2.5f}) {
    DominantFieldAccumulator acc;
    for (int n = 1; n <= 5; ++n) {
      acc.add(s);
      HS_EXPECT_NEAR(acc.value(), s, 1e-6f);
    }
  }

  // Denominator floor: fields too small to square above FIELD_DOMINANT_DEN_EPS
  // report 0 rather than dividing.
  DominantFieldAccumulator tiny;
  tiny.add(1e-6f);
  tiny.add(-1e-6f);
  HS_EXPECT_NEAR(tiny.value(), 0.0f, 0.0f);
}

/**
 * @brief Verifies the strongest contribution dominates a mixed overlap without
 *        the two stacking.
 */
inline void test_dominant_field_strongest_wins_without_stacking() {
  DominantFieldAccumulator acc;
  acc.add(3.0f);
  acc.add(1.0f);
  const float blended = acc.value();
  HS_EXPECT_LT(blended, 3.0f); // not the plain max
  HS_EXPECT_GT(blended, 2.0f); // but far nearer the strong field
  HS_EXPECT_LT(blended, 4.0f); // and never the sum
  HS_EXPECT_GT(blended, 1.0f);

  DominantFieldAccumulator lopsided;
  lopsided.add(10.0f);
  lopsided.add(0.1f);
  HS_EXPECT_NEAR(lopsided.value(), 10.0f, 0.05f);
}

/**
 * @brief Verifies the blend crosses an opposite-signed cancellation smoothly
 *        rather than jumping the way a hard max by magnitude would.
 * @details A fixed +1 field is swept against a partner from -1.5 to -0.5. At
 *          exact cancellation the blend is 0; a hard max by magnitude would
 *          instead flip from -1.5 to +1 across that point, a step of 2.5. The
 *          sweep asserts every consecutive step stays far below that, that the
 *          output changes sign across the crossing, and that it moves
 *          substantially overall (so a constant-0 accumulator cannot pass).
 */
inline void test_dominant_field_sign_crossing_is_continuous() {
  const int steps = 1000;
  const float lo = -1.5f, hi = -0.5f;

  auto blend = [](float b) {
    DominantFieldAccumulator acc;
    acc.add(1.0f);
    acc.add(b);
    return acc.value();
  };

  float prev = blend(lo);
  const float first = prev;
  float worst_step = 0.0f;
  for (int i = 1; i <= steps; ++i) {
    const float b = lo + (hi - lo) * (float)i / (float)steps;
    const float cur = blend(b);
    worst_step = std::max(worst_step, std::fabs(cur - prev));
    prev = cur;
  }
  const float last = prev;

  HS_EXPECT_LT(worst_step, 0.01f);
  HS_EXPECT_LT(first, 0.0f);
  HS_EXPECT_GT(last, 0.0f);
  HS_EXPECT_GT(last - first, 0.5f);
  HS_EXPECT_NEAR(blend(-1.0f), 0.0f, 1e-6f);
}

// ============================================================================
// bump_field — cap falloff, fast reject, envelope gating
// ============================================================================

/** @brief Verifies the bump threshold tracks the effective footprint. */
inline void test_bump_field_threshold_sync() {
  Animation::BumpParams p;
  p.radius = 0.8f;
  p.envelope = 0.25f;
  p.sync();
  HS_EXPECT_NEAR(p.cos_radius, std::cos(0.2f), 1e-6f);

  p.envelope = 0.75f;
  p.sync();
  HS_EXPECT_NEAR(p.cos_radius, std::cos(0.6f), 1e-6f);
}

/**
 * @brief Verifies the drape push along the axis meridian: zero for the ring
 *        through the center, peak bow between center and edge, zero at the
 *        footprint edge, exactly 0 outside the fast-reject cap.
 * @details The center sits at colatitude 60° about the +y axis, so the
 *          sampled meridian has room on both sides. On the meridian a point
 *          at offset y pushes by sign(y) * (R - |y|) * sin(pi |y| / R): at
 *          |y| = R/2 that is (R/2) * sin(pi/2) = R/2, antisymmetric across
 *          the center.
 */
inline void test_bump_field_drapes_over_ball() {
  const float c_lat = PI_F / 3.0f;
  Animation::BumpParams p;
  p.center = Vector(std::sin(c_lat), std::cos(c_lat), 0.0f);
  p.axis = Vector(0, 1, 0);
  p.radius = 0.5f;
  p.amplitude = 1.0f;
  p.envelope = 1.0f;
  p.sync();

  // Sample the meridian at colatitude offset d from the center.
  auto at = [&](float d) {
    return bump_field(Vector(std::sin(c_lat + d), std::cos(c_lat + d), 0.0f),
                      p);
  };

  HS_EXPECT_NEAR(at(0.0f), 0.0f, 2e-3f);
  HS_EXPECT_NEAR(at(p.radius * 0.5f), p.radius * 0.5f, 2e-3f);
  HS_EXPECT_NEAR(at(-p.radius * 0.5f), -p.radius * 0.5f, 2e-3f);
  HS_EXPECT_NEAR(at(p.radius * 0.98f), 0.0f, 5e-3f);
  HS_EXPECT_NEAR(at(p.radius + 0.05f), 0.0f, 1e-7f);
  HS_EXPECT_NEAR(at(-p.radius - 0.05f), 0.0f, 1e-7f);

  // Half gain halves the drape; a huge gain saturates at the full boundary
  // clearance (R - |y|); zero gain switches the bump off entirely.
  p.amplitude = 0.5f;
  HS_EXPECT_NEAR(at(p.radius * 0.5f), p.radius * 0.25f, 2e-3f);
  p.amplitude = 100.0f;
  const float y_low = p.radius * 0.1f;
  HS_EXPECT_NEAR(at(y_low), p.radius - y_low, 2e-3f);
  p.amplitude = 0.0f;
  HS_EXPECT_NEAR(at(p.radius * 0.5f), 0.0f, 1e-7f);
}

/**
 * @brief Verifies the bulge profile along one ring is the boundary arc scaled
 *        by the ring's drape weight, so the bow is round rather than pointed.
 */
inline void test_bump_field_round_bulge_along_ring() {
  const float c_lat = PI_F / 2.0f;
  Animation::BumpParams p;
  p.center = Vector(std::sin(c_lat), std::cos(c_lat), 0.0f);
  p.axis = Vector(0, 1, 0);
  p.radius = 0.5f;
  p.envelope = 1.0f;
  p.sync();

  // Ring half a footprint below the center (drape weight sin(pi/2) = 1),
  // sampled at azimuth offset 0.3 rad about the axis. There x ~ 0.3 * sin of
  // the ring's colatitude, and the push is sqrt(R^2 - x^2) - |y|.
  const float y = p.radius * 0.5f;
  const float lat = c_lat + y;
  const float az = 0.3f;
  const Vector v(std::sin(lat) * std::cos(az), std::cos(lat),
                 std::sin(lat) * std::sin(az));
  const float x = az * std::sin(lat);
  const float expected = std::sqrt(p.radius * p.radius - x * x) - y;
  HS_EXPECT_NEAR(bump_field(v, p), expected, 5e-3f);
}

/** @brief Verifies the cached-ring evaluator matches the generic bump field. */
inline void test_bump_field_precomputed_y_parity() {
  const float c_lat = 1.1f;
  Animation::BumpParams p;
  p.center = Vector(std::sin(c_lat), std::cos(c_lat), 0.0f);
  p.axis = Vector(0, 1, 0);
  p.radius = 0.55f;
  p.amplitude = 0.8f;
  p.envelope = 0.7f;
  p.sync();

  const float center_colat =
      fast_acos(hs::clamp(dot(p.axis, p.center), -1.0f, 1.0f));
  for (float lat = 0.4f; lat <= 1.8f; lat += 0.07f) {
    for (float az = -0.6f; az <= 0.6f; az += 0.09f) {
      Vector v(std::sin(lat) * std::cos(az), std::cos(lat),
               std::sin(lat) * std::sin(az));
      const float y = lat - center_colat;
      HS_EXPECT_NEAR(bump_field_with_y(v, p, y), bump_field(v, p), 2e-3f);
    }
  }
}

/**
 * @brief Verifies the lifecycle envelope scales the footprint (a half-envelope
 *        cap clears to its shrunken boundary) and gates the field to exactly 0
 *        when fully faded.
 */
inline void test_bump_field_envelope_gates() {
  const float c_lat = PI_F / 3.0f;
  Animation::BumpParams p;
  p.center = Vector(std::sin(c_lat), std::cos(c_lat), 0.0f);
  p.axis = Vector(0, 1, 0);
  p.radius = 0.5f;
  p.sync();

  // Half the shrunken footprint below the center on the axis meridian.
  const float d = c_lat + 0.125f;
  const Vector v(std::sin(d), std::cos(d), 0.0f);
  p.envelope = 0.5f; // effective radius 0.25
  p.sync();
  HS_EXPECT_NEAR(bump_field(v, p), 0.125f, 2e-3f);
  p.envelope = 0.0f;
  p.sync();
  HS_EXPECT_NEAR(bump_field(v, p), 0.0f, 1e-7f);
}

/**
 * @brief Sweeps footprint, envelope, gain, cap centre and stack axis to verify
 *        |bump_field| never exceeds field_bound().
 * @details field_bound() sizes the displacement culls, so a bump that pushes
 *          past its bound is culled while still displacing — the same safety
 *          role noise_product_field's bound is swept for. The 49,152 samples
 *          report through violation counters rather than per-sample assertions,
 *          which would swamp the module's assertion floor.
 */
inline void test_bump_field_bound_is_conservative() {
  hs::random().seed(20260803);
  auto random_unit = []() {
    for (;;) {
      const Vector r(hs::rand_f(-1.0f, 1.0f), hs::rand_f(-1.0f, 1.0f),
                     hs::rand_f(-1.0f, 1.0f));
      if (r.length() > 0.1f)
        return r.normalized();
    }
  };

  float worst_ratio = 0.0f;
  int nonfinite = 0;
  int over_bound = 0;
  for (int trial = 0; trial < 96; ++trial) {
    Animation::BumpParams p;
    p.center = random_unit();
    p.axis = random_unit();
    p.radius = hs::rand_f(0.0f, 1.2f);
    p.envelope = hs::rand_f(0.0f, 1.0f);
    // Well past the Ball Amp slider's 0.8 ceiling, into the saturating regime.
    p.amplitude = hs::rand_f(0.0f, 3.0f);
    p.sync();
    const float bound = p.field_bound();

    for (int i = 0; i < 512; ++i) {
      const float f = bump_field(random_unit(), p);
      nonfinite += !std::isfinite(f);
      over_bound += std::fabs(f) > bound;
      if (bound > 0.001f)
        worst_ratio = std::max(worst_ratio, std::fabs(f) / bound);
    }
  }
  HS_EXPECT_EQ(nonfinite, 0);
  HS_EXPECT_EQ(over_bound, 0);
  // The sweep must approach the bound, otherwise it would pass on a field that
  // never leaves zero.
  HS_EXPECT_GT(worst_ratio, 0.5f);
}

// ============================================================================
// noise_product_field — two-octave product parity; amplitude short-circuit;
// field_bound stays conservative across the parameter space
// ============================================================================

/**
 * @brief Verifies the field matches the hand-computed two-octave product and
 *        that ~0 amplitude short-circuits to exactly 0.
 */
inline void test_noise_product_field_parity() {
  Animation::NoiseProductParams p;
  p.amplitude = 0.25f;
  p.scale1 = 1.5f;
  p.scale2 = 3.0f;
  p.time = 2.0f;
  p.noise.SetSeed(1234);

  const Vector samples[] = {Vector(1, 0, 0), Vector(0, 1, 0),
                            Vector(0.4f, 0.6f, 0.7f).normalized()};
  constexpr size_t NUM_SAMPLES = sizeof(samples) / sizeof(samples[0]);
  float total = 0.0f;
  float base[NUM_SAMPLES] = {};
  for (size_t i = 0; i < NUM_SAMPLES; ++i) {
    const Vector &v = samples[i];
    float n1 = p.noise.GetNoise(v.x * p.scale1, v.y * p.scale1,
                                v.z * p.scale1 + p.time);
    float n2 = p.noise.GetNoise(
        v.x * p.scale2 + Animation::NoiseProductParams::OCTAVE2_OFFSET,
        v.y * p.scale2, v.z * p.scale2 + p.time);
    float expected = p.amplitude * n1 * n2;
    base[i] = noise_product_field(v, p);
    HS_EXPECT_NEAR(base[i], expected, 1e-6f);
    total += std::fabs(expected);
  }
  HS_EXPECT_GT(total, 1e-5f);

  // Independent of the octave formula: amplitude scales the field linearly and
  // its sign mirrors it, so amplitude cannot leak into where the noise is
  // sampled.
  p.amplitude = -0.5f;
  for (size_t i = 0; i < NUM_SAMPLES; ++i)
    HS_EXPECT_NEAR(noise_product_field(samples[i], p), -2.0f * base[i], 1e-6f);

  p.amplitude = 0.0f;
  HS_EXPECT_NEAR(noise_product_field(samples[0], p), 0.0f, 1e-7f);
  HS_EXPECT_NEAR(p.field_bound(), 0.0f, 1e-7f);
}

/**
 * @brief Sweeps seeds, octave scales, field time, amplitude, and sample points
 *        to verify |noise_product_field| never exceeds field_bound().
 * @details field_bound() = |amplitude| holds only if the generator's own output
 * stays within [-1, 1] on both octaves, so the sweep measures that claim
 * directly rather than trusting it. The bound sizes displacement culls; a point
 * outside it would be culled while still displaced. The 49,152 samples report
 * through violation counters rather than per-sample assertions: a per-sample
 * HS_EXPECT would put ~200k assertions into the module's floor and leave the
 * floor gate blind to every other case in the module being deleted.
 */
inline void test_noise_product_field_bound_is_conservative() {
  hs::random().seed(20260801);
  float worst_ratio = 0.0f;
  float worst_octave = 0.0f;
  int nonfinite = 0;
  int over_bound = 0;
  int octave_out_of_range = 0;
  for (int trial = 0; trial < 96; ++trial) {
    const float amplitude = hs::rand_f(-2.0f, 2.0f);
    const float scale1 = hs::rand_f(0.25f, 12.0f);
    const float scale2 = hs::rand_f(0.25f, 12.0f);
    const float time = hs::rand_f(-200.0f, 200.0f);
    const int seed = hs::rand_int(0, 100000);

    Animation::NoiseProductParams p;
    p.amplitude = amplitude;
    p.scale1 = scale1;
    p.scale2 = scale2;
    p.time = time;
    p.noise.SetSeed(seed);
    const float bound = p.field_bound();

    for (int i = 0; i < 512; ++i) {
      const float rx = hs::rand_f(-1.0f, 1.0f);
      const float ry = hs::rand_f(-1.0f, 1.0f);
      const float rz = hs::rand_f(-1.0f, 1.0f);
      const Vector r(rx, ry, rz);
      if (r.length() < 0.1f)
        continue;
      const Vector v = r.normalized();

      const float f = noise_product_field(v, p);
      nonfinite += !std::isfinite(f);
      over_bound += std::fabs(f) > bound;

      const float n1 = p.noise.GetNoise(v.x * p.scale1, v.y * p.scale1,
                                        v.z * p.scale1 + p.time);
      const float n2 = p.noise.GetNoise(
          v.x * p.scale2 + Animation::NoiseProductParams::OCTAVE2_OFFSET,
          v.y * p.scale2, v.z * p.scale2 + p.time);
      octave_out_of_range += std::fabs(n1) > 1.0f;
      octave_out_of_range += std::fabs(n2) > 1.0f;
      worst_octave = std::max(worst_octave, std::fabs(n1));
      worst_octave = std::max(worst_octave, std::fabs(n2));
      if (bound > 0.001f)
        worst_ratio = std::max(worst_ratio, std::fabs(f) / bound);
    }
  }
  HS_EXPECT_EQ(nonfinite, 0);
  HS_EXPECT_EQ(over_bound, 0);
  HS_EXPECT_EQ(octave_out_of_range, 0);
  // The sweep must reach a substantial fraction of the bound, otherwise it
  // would pass on a field that never leaves zero.
  HS_EXPECT_GT(worst_octave, 0.5f);
  HS_EXPECT_GT(worst_ratio, 0.1f);
}

// ============================================================================
// Animation::BallDrop — pole-to-pole traversal, envelope, slot reclaim
// ============================================================================

/**
 * @brief Verifies a spawned BallDrop traverses its bump from the frame's pole
 *        to the opposite pole along the requested meridian, ramps its envelope
 *        0 → 1 → 0, and frees its pool slot on completion.
 */
inline void test_ball_drop_traverses_and_reclaims() {
  Timeline tl;
  global_timeline_t = 0;
  hs_test::StubEffect fx(8, 8);
  Canvas cv(fx);

  Orientation<> ori;
  const Vector pole(0, 1, 0);
  const int duration = 40;

  BallDropTransformer<1> balls(tl);
  balls.init_storage(persistent_arena);
  balls.template_params.radius = 0.5f;
  HS_EXPECT_TRUE(balls.spawn(0, ori, pole, 0.0f, duration) != nullptr);

  tl.step(cv);
  balls.prepare_frame();
  const float early = balls.field(pole);
  HS_EXPECT_NEAR(early, 0.0f, 0.02f);

  for (int i = 1; i < duration / 2; ++i)
    tl.step(cv);
  balls.prepare_frame();
  // Mid-fall: the bump sits at the world equator point (1, 0, 0) (azimuth 0)
  // at full envelope. The pole reads ~0 (outside the cap); half a footprint
  // below the center on its meridian, the drape push toward the lower
  // boundary arc is R - R/2, and the mirrored point above reads its negation.
  HS_EXPECT_NEAR(balls.field(pole), 0.0f, 1e-5f);
  const float half_r = 0.5f * 0.5f;
  const Vector below(std::cos(half_r), -std::sin(half_r), 0.0f);
  const Vector above(std::cos(half_r), std::sin(half_r), 0.0f);
  HS_EXPECT_NEAR(balls.field(below), half_r, 3e-3f);
  HS_EXPECT_NEAR(balls.field(above), -half_r, 3e-3f);
  HS_EXPECT_GT(balls.field_bound(), 0.25f);

  for (int i = duration / 2; i <= duration + 2; ++i)
    tl.step(cv);
  HS_EXPECT_TRUE(balls.active_count() == 0);
  HS_EXPECT_TRUE(balls.spawn(0, ori, pole, 1.0f, duration) != nullptr);
}

// ============================================================================
// Animation::NoiseProduct — time integrates by speed
// ============================================================================

/**
 * @brief Verifies NoiseProduct advances params.time by speed each step.
 */
inline void test_noise_product_integrates_time() {
  hs_test::StubEffect fx(8, 8);
  Canvas cv(fx);
  Animation::NoiseProductParams p;
  p.speed = 0.03f;
  Animation::NoiseProduct anim(p);
  for (int i = 0; i < 5; ++i)
    anim.step(cv);
  HS_EXPECT_NEAR(p.time, 5 * 0.03f, 1e-5f);
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every transformers test case.
 * @return The module's failure count.
 */
inline int run_transformers_tests() {
  hs_test::ModuleFixture fixture("transformers");

  test_orient_transformer_identity();
  test_orient_transformer_known_rotation();
  test_mobius_identity_roundtrip();
  test_mobius_known_rotation();
  test_mobius_matches_double_precision_oracle();
  test_mobius_poles_map_to_coefficient_ratios();
  test_gnomonic_mobius_identity_roundtrip();
  test_gnomonic_mobius_known_rotation();
  test_gnomonic_mobius_signed_zero_equator();
  test_ripple_zero_amplitude_is_identity();
  test_ripple_center_point_is_identity();
  test_ripple_active_rotates_on_sphere();
  test_ripple_threshold_reject_path();
  test_ripple_decay_attenuates();
  test_ripple_threshold_boundary();
  test_ripple_small_angle_series_matches_exact();
  test_ripple_threshold_collapses_past_pi();
  test_transforms_nonfinite_passes_through_identity();
  test_noise_zero_amplitude_is_identity();
  test_noise_active_stays_on_sphere();
  test_noise_cross_hemisphere_cap();
  test_transformer_no_entities_is_identity();
  test_transformer_spawn_applies_and_composes();
  test_transformer_nonpinned_slot_reclaimed_after_compaction();
  test_transformer_slots_released_by_timeline_clear();
  test_transformer_recycled_slot_composes_in_spawn_order();
  test_field_transformer_no_entities_is_zero();
  test_field_transformer_sums_and_bounds();
  test_field_transformer_active_params_spawn_order();
  test_field_transformer_slot_reclaimed();
  test_field_transformer_storage_survives_arena_rewind();
  test_dominant_field_identities();
  test_dominant_field_strongest_wins_without_stacking();
  test_dominant_field_sign_crossing_is_continuous();
  test_bump_field_threshold_sync();
  test_bump_field_drapes_over_ball();
  test_bump_field_round_bulge_along_ring();
  test_bump_field_precomputed_y_parity();
  test_bump_field_envelope_gates();
  test_bump_field_bound_is_conservative();
  test_noise_product_field_parity();
  test_noise_product_field_bound_is_conservative();
  test_ball_drop_traverses_and_reclaims();
  test_noise_product_integrates_time();

  return fixture.result();
}

} // namespace transformers_tests
} // namespace hs_test
