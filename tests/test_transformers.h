/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/transformers.h — the pure geometry transform functions
 * and the adapter/manager wrappers:
 *   - OrientTransformer       : identity orientation is a no-op.
 *   - mobius_transform        : identity Mobius round-trips through stereo.
 *   - gnomonic_mobius_transform: identity round-trips through gnomonic.
 *   - ripple_transform        : amplitude 0 and center-point degeneracies are
 *                               no-ops; an active ripple rotates on-sphere.
 *   - noise_transform         : amplitude≈0 is a no-op; otherwise stays unit.
 *   - Transformer<>           : no active entities → identity.
 */
#pragma once

#include "core/transformers.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace transformers_tests {

inline bool finite_vec(const Vector &v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

// ============================================================================
// OrientTransformer
// ============================================================================

inline void test_orient_transformer_identity() {
  Orientation<32> ori; // default = identity quaternion
  OrientTransformer<32> ot(ori);

  const Vector samples[] = {Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1),
                            Vector(0.5f, 0.1f, 0.3f).normalized()};
  for (const Vector &v : samples) {
    Vector r = ot(v);
    HS_EXPECT_NEAR(r.x, v.x, 1e-5f);
    HS_EXPECT_NEAR(r.y, v.y, 1e-5f);
    HS_EXPECT_NEAR(r.z, v.z, 1e-5f);
  }
}

// ============================================================================
// mobius_transform / gnomonic_mobius_transform — identity round-trips
// ============================================================================

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

inline void test_gnomonic_mobius_identity_roundtrip() {
  MobiusParams id;
  Vector v = Vector(0.3f, 0.7f, 0.2f).normalized(); // y>0 hemisphere
  Vector r = gnomonic_mobius_transform(v, id);
  HS_EXPECT_TRUE(finite_vec(r));
  HS_EXPECT_NEAR(r.x, v.x, 2e-3f);
  HS_EXPECT_NEAR(r.y, v.y, 2e-3f);
  HS_EXPECT_NEAR(r.z, v.z, 2e-3f);
}

// ============================================================================
// ripple_transform
// ============================================================================

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

inline void test_ripple_center_point_is_identity() {
  // A point coincident with the ripple center has a degenerate rotation axis
  // (cross(center, center) == 0) and must be returned unchanged.
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

inline void test_ripple_active_rotates_on_sphere() {
  RippleParams p;
  p.center = Vector(0, 1, 0);
  p.amplitude = 0.5f;
  p.phase = PI_F * 0.5f; // peak of the wavelet at d == 90°
  p.decay = 0.0f;        // no attenuation, so the rotation is significant
  p.thickness = 1.0f;
  // default thresholds (min=1, max=-1) disable the fast reject → ripple applies

  Vector v = Vector(1, 0, 0); // 90° from center → at the wavelet peak
  Vector r = ripple_transform(v, p);

  HS_EXPECT_TRUE(finite_vec(r));
  HS_EXPECT_NEAR(r.length(), 1.0f, 1e-4f); // rotation preserves unit length
  // The point actually moved.
  float moved =
      std::abs(r.x - v.x) + std::abs(r.y - v.y) + std::abs(r.z - v.z);
  HS_EXPECT_GT(moved, 1e-2f);
}

// ============================================================================
// noise_transform
// ============================================================================

inline void test_noise_zero_amplitude_is_identity() {
  NoiseParams p;
  p.amplitude = 0.0f;
  Vector v = Vector(0.2f, 0.5f, 0.84f).normalized();
  Vector r = noise_transform(v, p);
  HS_EXPECT_NEAR(r.x, v.x, 1e-6f);
  HS_EXPECT_NEAR(r.y, v.y, 1e-6f);
  HS_EXPECT_NEAR(r.z, v.z, 1e-6f);
}

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

inline void test_transformer_no_entities_is_identity() {
  Timeline<32> tl;
  RippleTransformer<32, 8> rt(tl);
  Vector v = Vector(0.6f, 0.2f, 0.77f).normalized();
  Vector r = rt.transform(v); // nothing spawned → all entities inactive
  HS_EXPECT_NEAR(r.x, v.x, 1e-6f);
  HS_EXPECT_NEAR(r.y, v.y, 1e-6f);
  HS_EXPECT_NEAR(r.z, v.z, 1e-6f);
}

// ============================================================================
// Runner
// ============================================================================

inline int run_transformers_tests() {
  auto scope = hs_test::begin_module("transformers");

  test_orient_transformer_identity();
  test_mobius_identity_roundtrip();
  test_gnomonic_mobius_identity_roundtrip();
  test_ripple_zero_amplitude_is_identity();
  test_ripple_center_point_is_identity();
  test_ripple_active_rotates_on_sphere();
  test_noise_zero_amplitude_is_identity();
  test_noise_active_stays_on_sphere();
  test_transformer_no_entities_is_identity();

  return hs_test::end_module(scope);
}

} // namespace transformers_tests
} // namespace hs_test

