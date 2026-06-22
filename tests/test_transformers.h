/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/transformers.h — the pure geometry transform functions
 * and the adapter/manager wrappers:
 *   - OrientTransformer       : identity orientation is a no-op; a known 90°
 *                               rotation maps to its hand-computed image.
 *   - mobius_transform        : identity Mobius round-trips through stereo; the
 *                               1/z map realizes a 180° rotation about x.
 *   - gnomonic_mobius_transform: identity round-trips through gnomonic; the -z
 *                               map realizes a 180° rotation about y.
 *   - ripple_transform        : amplitude 0 and center-point degeneracies are
 *                               no-ops; an active ripple rotates on-sphere; the
 *                               prepared-threshold fast-reject band applies
 *                               in-band points and rejects off-band ones.
 *   - noise_transform         : amplitude≈0 is a no-op; otherwise stays unit.
 *   - Transformer<>           : no active entities → identity.
 */
#pragma once

#include "core/transformers.h"
#include "tests/test_harness.h"

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

// ============================================================================
// OrientTransformer
// ============================================================================

/**
 * @brief Verifies an identity orientation leaves every sampled direction
 *        unchanged.
 */
inline void test_orient_transformer_identity() {
  Orientation<> ori; // default = identity quaternion
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
  Orientation<> ori(make_rotation(Vector(0, 1, 0), PI_F * 0.5f)); // +90° about y
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
  MobiusParams id; // a=1, b=0, c=0, d=1 → identity
  Vector v = Vector(0.5f, 0.1f, 0.3f).normalized(); // away from the N pole
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
 * @brief Verifies the identity Mobius map round-trips a point through
 *        gnomonic projection back to itself.
 */
inline void test_gnomonic_mobius_identity_roundtrip() {
  MobiusParams id;
  Vector v = Vector(0.3f, 0.7f, 0.2f).normalized(); // y>0 hemisphere
  Vector r = gnomonic_mobius_transform(v, id);
  HS_EXPECT_TRUE(finite_vec(r));
  HS_EXPECT_NEAR(r.x, v.x, 2e-3f);
  HS_EXPECT_NEAR(r.y, v.y, 2e-3f);
  HS_EXPECT_NEAR(r.z, v.z, 2e-3f);
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
  Vector v = Vector(0.3f, 0.7f, 0.2f).normalized(); // upper hemisphere (v.y>0)
  Vector r = gnomonic_mobius_transform(v, neg);
  HS_EXPECT_TRUE(finite_vec(r));
  HS_EXPECT_NEAR(r.x, -v.x, 1e-3f);
  HS_EXPECT_NEAR(r.y, v.y, 1e-3f);
  HS_EXPECT_NEAR(r.z, -v.z, 1e-3f);
}

// ============================================================================
// ripple_transform
// ============================================================================

/**
 * @brief Verifies zero amplitude short-circuits the ripple to a no-op.
 */
inline void test_ripple_zero_amplitude_is_identity() {
  RippleParams p;
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
  RippleParams p;
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
  RippleParams p;
  p.center = Vector(0, 1, 0);
  p.amplitude = 0.5f;
  p.phase = PI_F * 0.5f; // peak of the wavelet at d == 90°
  p.decay = 0.0f;        // no attenuation, so the rotation is significant
  p.thickness = 1.0f;
  // default thresholds (min=1, max=-1) disable the fast reject → ripple always
  // applies. test_ripple_threshold_reject_path covers the prepared-bounds path.

  Vector v = Vector(1, 0, 0); // 90° from center → at the wavelet peak
  Vector r = ripple_transform(v, p);

  HS_EXPECT_TRUE(finite_vec(r));
  HS_EXPECT_NEAR(r.length(), 1.0f, 1e-4f); // rotation preserves unit length
  // The point actually moved.
  float moved =
      std::abs(r.x - v.x) + std::abs(r.y - v.y) + std::abs(r.z - v.z);
  HS_EXPECT_GT(moved, 1e-2f);
}

/**
 * @brief Exercises the prepare_thresholds() fast-reject band the production
 *        renderer relies on (prepare_frame() calls it per active entity).
 * @details test_ripple_active_rotates_on_sphere deliberately leaves the default
 *          degenerate bounds (min=1, max=-1) in place, so the reject `if` is
 *          never taken there. Here real bounds are computed: a point at the
 *          wavelet peak (inside the band) is rotated, while points on either
 *          side — closer to the center than d_min, and farther than d_max —
 *          take the two reject legs and return unchanged. A prepare_thresholds()
 *          that swapped or mis-derived its bounds would either reject the peak
 *          (no move) or pass the off-band points (they move), failing here.
 */
inline void test_ripple_threshold_reject_path() {
  RippleParams p;
  p.center = Vector(0, 1, 0); // north pole, so cos_d == dot(v, center) == v.y
  p.amplitude = 0.5f;
  p.phase = PI_F * 0.5f; // wavelet peak at d == 90°
  p.decay = 0.0f;        // no attenuation, so an applied rotation is significant
  p.thickness = 0.4f;    // band ≈ [phase - 0.4, phase + 0.4] rad
  p.prepare_thresholds();

  // The prepared bounds must be a real (non-degenerate) window, else the reject
  // legs below can never fire and the test would silently prove nothing.
  HS_EXPECT_LT(p.cos_threshold_min, 1.0f);
  HS_EXPECT_GT(p.cos_threshold_max, -1.0f);

  // In-band: at the peak (d == phase), the ripple rotates the point off itself.
  const Vector in_band = Vector(1, 0, 0); // d == 90°, exactly at the peak
  const Vector r_in = ripple_transform(in_band, p);
  HS_EXPECT_TRUE(finite_vec(r_in));
  HS_EXPECT_NEAR(r_in.length(), 1.0f, 1e-4f);
  const float moved_in = std::abs(r_in.x - in_band.x) +
                         std::abs(r_in.y - in_band.y) +
                         std::abs(r_in.z - in_band.z);
  HS_EXPECT_GT(moved_in, 1e-2f);

  // Out-of-band toward the center (d < d_min): cos_d > cos_threshold_min leg.
  const float d_near = p.phase - p.thickness - 0.3f; // safely inside d_min, > 0
  const Vector near_c = Vector(std::sin(d_near), std::cos(d_near), 0.0f);
  const Vector r_near = ripple_transform(near_c, p);
  HS_EXPECT_NEAR(r_near.x, near_c.x, 1e-6f);
  HS_EXPECT_NEAR(r_near.y, near_c.y, 1e-6f);
  HS_EXPECT_NEAR(r_near.z, near_c.z, 1e-6f);

  // Out-of-band away from the center (d > d_max): cos_d < cos_threshold_max leg.
  const float d_far = p.phase + p.thickness + 0.3f; // safely past d_max, < PI
  const Vector far_c = Vector(std::sin(d_far), std::cos(d_far), 0.0f);
  const Vector r_far = ripple_transform(far_c, p);
  HS_EXPECT_NEAR(r_far.x, far_c.x, 1e-6f);
  HS_EXPECT_NEAR(r_far.y, far_c.y, 1e-6f);
  HS_EXPECT_NEAR(r_far.z, far_c.z, 1e-6f);
}

// ============================================================================
// noise_transform
// ============================================================================

/**
 * @brief Verifies zero amplitude short-circuits the noise warp to a no-op.
 */
inline void test_noise_zero_amplitude_is_identity() {
  NoiseParams p;
  p.amplitude = 0.0f;
  Vector v = Vector(0.2f, 0.5f, 0.84f).normalized();
  Vector r = noise_transform(v, p);
  HS_EXPECT_NEAR(r.x, v.x, 1e-6f);
  HS_EXPECT_NEAR(r.y, v.y, 1e-6f);
  HS_EXPECT_NEAR(r.z, v.z, 1e-6f);
}

/**
 * @brief Verifies an active noise warp keeps every sample finite and on the
 *        unit sphere.
 */
inline void test_noise_active_stays_on_sphere() {
  NoiseParams p;
  p.amplitude = 0.5f;
  p.scale = 4.0f;
  p.time = 1.0f;
  const Vector samples[] = {Vector(1, 0, 0), Vector(0, 1, 0),
                            Vector(0.4f, 0.6f, 0.7f).normalized()};
  for (const Vector &v : samples) {
    Vector r = noise_transform(v, p);
    HS_EXPECT_TRUE(finite_vec(r));
    HS_EXPECT_NEAR(r.length(), 1.0f, 1e-3f);
  }
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
  Vector r = rt.transform(v); // nothing spawned → all entities inactive
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
  nt.template_params.amplitude = 0.6f;
  nt.template_params.scale = 4.0f;
  nt.template_params.time = 1.0f;
  nt.template_params.sync();

  const Vector samples[] = {Vector(1, 0, 0), Vector(0, 0, 1),
                            Vector(0.4f, 0.6f, 0.7f).normalized()};

  // No active entities → identity for every sample.
  for (const Vector &v : samples) {
    Vector r = nt.transform(v);
    HS_EXPECT_NEAR(r.x, v.x, 1e-6f);
    HS_EXPECT_NEAR(r.y, v.y, 1e-6f);
    HS_EXPECT_NEAR(r.z, v.z, 1e-6f);
  }

  // One active entity → transform applies the warp; sum the displacement across
  // samples so a single noise zero-crossing can't make the test flaky. Capture
  // the single-entity output per sample to compare against the composed result.
  HS_EXPECT_TRUE(nt.spawn(0) != nullptr);
  constexpr size_t kN = sizeof(samples) / sizeof(samples[0]);
  Vector single[kN];
  float total_moved = 0.0f;
  for (size_t i = 0; i < kN; ++i) {
    const Vector &v = samples[i];
    single[i] = nt.transform(v);
    HS_EXPECT_TRUE(finite_vec(single[i]));
    HS_EXPECT_NEAR(single[i].length(), 1.0f, 1e-3f);
    total_moved += std::abs(single[i].x - v.x) +
                   std::abs(single[i].y - v.y) + std::abs(single[i].z - v.z);
  }
  HS_EXPECT_GT(total_moved, 1e-2f);

  // A second active entity → both compose. The composed output must differ from
  // the single-entity output (the second slot warps the already-displaced
  // point); if the active-slot list silently dropped the second entity — or
  // transform() applied only the first — this divergence would vanish. Sum
  // across samples so one near-coincident sample can't mask a real composition.
  HS_EXPECT_TRUE(nt.spawn(0) != nullptr);
  float total_divergence = 0.0f;
  for (size_t i = 0; i < kN; ++i) {
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
// Runner
// ============================================================================

/**
 * @brief Runs every transformers test case.
 * @return The module's failure count.
 */
inline int run_transformers_tests() {
  auto scope = hs_test::begin_module("transformers");

  test_orient_transformer_identity();
  test_orient_transformer_known_rotation();
  test_mobius_identity_roundtrip();
  test_mobius_known_rotation();
  test_gnomonic_mobius_identity_roundtrip();
  test_gnomonic_mobius_known_rotation();
  test_ripple_zero_amplitude_is_identity();
  test_ripple_center_point_is_identity();
  test_ripple_active_rotates_on_sphere();
  test_ripple_threshold_reject_path();
  test_noise_zero_amplitude_is_identity();
  test_noise_active_stays_on_sphere();
  test_transformer_no_entities_is_identity();
  test_transformer_spawn_applies_and_composes();

  return hs_test::end_module(scope);
}

} // namespace transformers_tests
} // namespace hs_test

