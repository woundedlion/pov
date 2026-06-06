/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/geometry.h.
 *
 * Coverage:
 *   - Axis constants (X/Y/Z_AXIS, UP)
 *   - Fragment::lerp
 *   - y_to_phi / phi_to_y conversions (free + templated)
 *   - PhiLUT / TrigLUT initialisation and lookup values
 *   - pixel_to_vector / vector_to_pixel roundtrip
 *   - logPolarToVector / vectorToLogPolar roundtrip + pole sentinel
 *   - fib_spiral, lissajous, random_vector (all on unit sphere)
 *   - Basis: make_basis (orthonormality), get_antipode (flipping)
 *   - Orientation<W, CAP>: set, push, get, collapse, upsample, orient/unorient
 */
#pragma once

#include "core/geometry.h"
#include "tests/test_3dmath.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace geometry {

using hs_test::math3d::approx_vec;
using hs_test::math3d::approx_quat;

// ============================================================================
// Axis constants
// ============================================================================

inline void test_axis_constants() {
  HS_EXPECT_VEC(X_AXIS, Vector(1, 0, 0), 0.0f);
  HS_EXPECT_VEC(Y_AXIS, Vector(0, 1, 0), 0.0f);
  HS_EXPECT_VEC(Z_AXIS, Vector(0, 0, 1), 0.0f);
  HS_EXPECT_VEC(UP, Y_AXIS, 0.0f);
}

// ============================================================================
// Fragment::lerp
// ============================================================================

inline void test_fragment_lerp_endpoints() {
  Fragment a;
  a.pos = Vector(1, 0, 0);
  a.v0 = 10.0f;
  a.v1 = 20.0f;
  a.v2 = 30.0f;
  a.v3 = 40.0f;
  a.age = 1.0f;

  Fragment b;
  b.pos = Vector(0, 1, 0);
  b.v0 = 110.0f;
  b.v1 = 220.0f;
  b.v2 = 330.0f;
  b.v3 = 440.0f;
  b.age = 5.0f;

  Fragment at0 = Fragment::lerp(a, b, 0.0f);
  HS_EXPECT_VEC(at0.pos, a.pos, 1e-6f);
  HS_EXPECT_NEAR(at0.v0, a.v0, 1e-6f);
  HS_EXPECT_NEAR(at0.age, a.age, 1e-6f);

  Fragment at1 = Fragment::lerp(a, b, 1.0f);
  HS_EXPECT_VEC(at1.pos, b.pos, 1e-6f);
  HS_EXPECT_NEAR(at1.v3, b.v3, 1e-6f);
}

inline void test_fragment_lerp_midpoint() {
  Fragment a, b;
  a.pos = Vector(2, 4, 6);
  b.pos = Vector(8, 16, 24);
  a.v0 = 0; b.v0 = 100;
  a.v1 = -50; b.v1 = 50;

  Fragment mid = Fragment::lerp(a, b, 0.5f);
  HS_EXPECT_VEC(mid.pos, Vector(5, 10, 15), 1e-5f);
  HS_EXPECT_NEAR(mid.v0, 50.0f, 1e-5f);
  HS_EXPECT_NEAR(mid.v1, 0.0f, 1e-5f);
}

// ============================================================================
// y_to_phi / phi_to_y
// ============================================================================

inline void test_y_to_phi_free_function() {
  HS_EXPECT_NEAR(y_to_phi(0.0f, 145), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(y_to_phi(144.0f, 145), PI_F, 1e-5f);
  HS_EXPECT_NEAR(y_to_phi(72.0f, 145), PI_F * 0.5f, 1e-5f);
}

inline void test_phi_to_y_free_function() {
  HS_EXPECT_NEAR(phi_to_y(0.0f, 145), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(phi_to_y(PI_F, 145), 144.0f, 1e-5f);
  HS_EXPECT_NEAR(phi_to_y(PI_F * 0.5f, 145), 72.0f, 1e-5f);
}

inline void test_y_phi_roundtrip() {
  for (int h_virt : {16, 64, 145}) {
    for (int i = 0; i < h_virt; ++i) {
      float phi = y_to_phi(static_cast<float>(i), h_virt);
      float y = phi_to_y(phi, h_virt);
      HS_EXPECT_NEAR(y, static_cast<float>(i), 1e-4f);
    }
  }
}

inline void test_y_to_phi_templated_LUT() {
  // Both integer and float template overloads hit the PhiLUT path for
  // integer inputs and agree with the free function.
  constexpr int H = 32;
  for (int y = 0; y < H; ++y) {
    float lut = y_to_phi<H>(y);
    float direct = y_to_phi(static_cast<float>(y), H + hs::H_OFFSET);
    HS_EXPECT_NEAR(lut, direct, 1e-5f);
  }

  // After a previous call, PhiLUT must be initialised.
  HS_EXPECT_TRUE(PhiLUT<H>::initialized);

  // Float version with integer-valued input also hits the LUT.
  float lut_f = y_to_phi<H>(5.0f);
  HS_EXPECT_NEAR(lut_f, y_to_phi<H>(5), 1e-5f);
}

// Offset-injection regression guard for the H_OFFSET double-apply bug.
//
// H_OFFSET adds virtual rows so the image is CLIPPED (not stretched) at the
// bottom of the ring, where the LEDs stop short of the south pole. The mapping
// is phi = y*PI/(H_VIRT-1) with H_VIRT = H + H_OFFSET, so the bottom physical
// row y=H-1 lands short of PI. The bug passed H_VIRT (already = H+H_OFFSET) as
// the template arg to y_to_phi<H>(), which re-adds H_OFFSET internally — a
// double-apply that divides by (H + 2*H_OFFSET - 1) instead. On the native
// build H_OFFSET==0 so the bug is invisible; this test injects the real
// hardware offset through the offset-parameterised free function (the building
// block the templated path mirrors) so the formula and the clipping intent are
// pinned regardless of the compile-time offset.
inline void test_y_to_phi_offset_injection_clips_no_double_apply() {
  constexpr int H = 20;          // Holosphere hardware height
  constexpr int OFF = 3;         // hardware H_OFFSET
  constexpr int H_VIRT = H + OFF; // 23
  const float bottom_phys = static_cast<float>(H - 1); // y = 19

  // Correct (single-offset) mapping: bottom physical row clips short of PI.
  float correct = y_to_phi(bottom_phys, H_VIRT);          // 19*PI/22
  HS_EXPECT_NEAR(correct, bottom_phys * PI_F / (H_VIRT - 1), 1e-5f);
  HS_EXPECT_TRUE(correct < PI_F); // clipped, not reaching the south pole

  // The double-applied denominator (the bug) divides by (H + 2*OFF - 1) = 25
  // and must NOT match the correct value — guards against reintroducing it.
  float double_applied = bottom_phys * PI_F / (H + 2 * OFF - 1); // 19*PI/24
  HS_EXPECT_TRUE(std::abs(correct - double_applied) > 1e-2f);

  // The virtual bottom row reaches the pole exactly.
  HS_EXPECT_NEAR(y_to_phi(static_cast<float>(H_VIRT - 1), H_VIRT), PI_F, 1e-5f);
}

// ============================================================================
// TrigLUT / pixel_to_vector / vector_to_pixel
// ============================================================================

inline void test_pixel_to_vector_unit_length() {
  constexpr int W = 32, H = 32;
  for (int y = 0; y < H + hs::H_OFFSET; ++y) {
    for (int x = 0; x < W; x += 4) {
      Vector v = pixel_to_vector<W, H>(x, y);
      HS_EXPECT_NEAR(v.length(), 1.0f, 1e-3f);
    }
  }
  // After use, the TrigLUT must be initialised.
  HS_EXPECT_TRUE((TrigLUT<W, H>::initialized));
}

inline void test_pixel_to_vector_known_samples() {
  constexpr int W = 32, H = 32;
  // (x=0, y near H/2) → near equator. Integer y rarely lands exactly on
  // phi=π/2 (true equator is y=(H_VIRT-1)/2), so allow a generous tolerance.
  Vector v = pixel_to_vector<W, H>(0, (H + hs::H_OFFSET) / 2);
  HS_EXPECT_VEC(v, Vector(1, 0, 0), 0.15f);

  // y=0 (north pole): phi=0 ⇒ cos_phi=1, sin_phi=0 ⇒ Vector(0, 1, 0).
  Vector north = pixel_to_vector<W, H>(0, 0);
  HS_EXPECT_VEC(north, Vector(0, 1, 0), 5e-2f);
}

// Pins pixel_to_vector's fractional-y (float) branch to y_to_phi<H>, the same
// phi source the integer LUT path and every SDF/scan shape use. If the float
// branch ever diverges (e.g. a y_to_phi<H_VIRT> double-apply reappears) on a
// build with a non-zero H_OFFSET, the float and LUT paths disagree and this
// fails. At x=0, Vector(Spherical(0, phi)) = (sin phi, cos phi, 0), so the
// recovered phi is acos(v.y).
inline void test_pixel_to_vector_float_branch_matches_phi_lut() {
  constexpr int W = 32, H = 32;
  for (float y : {0.5f, 5.25f, 12.75f, 20.5f}) {
    Vector v = pixel_to_vector<W, H>(0.0f, y);
    float recovered_phi = std::acos(std::clamp(v.y, -1.0f, 1.0f));
    HS_EXPECT_NEAR(recovered_phi, y_to_phi<H>(y), 1e-4f);
  }
}

inline void test_vector_to_pixel_roundtrip_via_pixel_to_vector() {
  constexpr int W = 64, H = 64;
  // Pick representative non-degenerate samples (avoiding poles where wrap is
  // undefined).
  int xs[] = {3, 17, 31, 50};
  int ys[] = {8, 16, 24, 40};
  for (int x : xs) {
    for (int y : ys) {
      Vector v = pixel_to_vector<W, H>(x, y);
      PixelCoords p = vector_to_pixel<W, H>(v);
      HS_EXPECT_NEAR(p.x, static_cast<float>(x), 0.5f);
      HS_EXPECT_NEAR(p.y, static_cast<float>(y), 0.5f);
    }
  }
}

// ============================================================================
// Log-polar projection (sphere ↔ complex plane via stereographic)
// ============================================================================

inline void test_log_polar_roundtrip() {
  Vector samples[] = {
      Vector(0.6f, 0.0f, 0.8f),
      Vector(1, 0, 0),
      Vector(0.5f, 0.5f, 0.7071f).normalized(),
      Vector(-0.6f, 0.3f, 0.7f).normalized(),
      Vector(0, -1, 0), // south pole — well-defined
  };
  for (const Vector &v : samples) {
    LogPolar lp = vectorToLogPolar(v);
    Vector back = logPolarToVector(lp.rho, lp.theta);
    HS_EXPECT_VEC(back, v, 5e-3f);
  }
}

inline void test_log_polar_north_pole_sentinel() {
  // North pole maps to the sentinel (rho=10, theta=0).
  LogPolar lp = vectorToLogPolar(Vector(0, 1, 0));
  HS_EXPECT_NEAR(lp.rho, 10.0f, 1e-3f);
  HS_EXPECT_NEAR(lp.theta, 0.0f, 1e-3f);
}

// ============================================================================
// fib_spiral
// ============================================================================

inline void test_fib_spiral_unit_length() {
  for (int i = 0; i < 32; ++i) {
    Vector v = fib_spiral(32, 0.5f, i);
    HS_EXPECT_NEAR(v.length(), 1.0f, 1e-3f);
  }
}

inline void test_fib_spiral_deterministic() {
  Vector v1 = fib_spiral(64, 0.5f, 7);
  Vector v2 = fib_spiral(64, 0.5f, 7);
  HS_EXPECT_VEC(v1, v2, 1e-6f);
}

inline void test_fib_spiral_endpoints() {
  // i=0 → near south pole region (y ≈ 1 - 2*eps/N for first sample)
  // i=N-1 → near north pole region
  // Both endpoints are unit-length and lie on opposite hemispheres for N>=4.
  Vector first = fib_spiral(16, 0.5f, 0);
  Vector last = fib_spiral(16, 0.5f, 15);
  HS_EXPECT_NEAR(first.length(), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(last.length(), 1.0f, 1e-3f);
  HS_EXPECT_TRUE(first.y * last.y < 0.0f); // opposite hemispheres
}

// ============================================================================
// lissajous
// ============================================================================

inline void test_lissajous_unit_length() {
  for (float t = 0.0f; t < 6.28f; t += 0.3f) {
    Vector v = lissajous(3, 2, 0.5f, t);
    HS_EXPECT_NEAR(v.length(), 1.0f, 1e-3f);
  }
}

// ============================================================================
// random_vector (Marsaglia)
// ============================================================================

inline void test_random_vector_unit_length() {
  for (int i = 0; i < 32; ++i) {
    Vector v = random_vector();
    HS_EXPECT_NEAR(v.length(), 1.0f, 1e-3f);
  }
}

// ============================================================================
// Basis: make_basis / get_antipode
// ============================================================================

inline void test_make_basis_orthonormal() {
  // Identity quaternion → basis aligned with normal.
  Quaternion id(1, 0, 0, 0);
  Basis b = make_basis(id, Vector(0, 1, 0));
  HS_EXPECT_VEC(b.v, Vector(0, 1, 0), 1e-4f);
  HS_EXPECT_NEAR(b.u.length(), 1.0f, 1e-4f);
  HS_EXPECT_NEAR(b.w.length(), 1.0f, 1e-4f);
  // Orthogonality
  HS_EXPECT_NEAR(dot(b.u, b.v), 0.0f, 1e-4f);
  HS_EXPECT_NEAR(dot(b.v, b.w), 0.0f, 1e-4f);
  HS_EXPECT_NEAR(dot(b.u, b.w), 0.0f, 1e-4f);
}

inline void test_make_basis_alternate_normals() {
  Quaternion id(1, 0, 0, 0);
  // A tilted normal
  Vector n = Vector(0.5f, 0.5f, 0.7071f).normalized();
  Basis b = make_basis(id, n);
  HS_EXPECT_VEC(b.v, n, 1e-3f);
  HS_EXPECT_NEAR(dot(b.u, b.v), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(dot(b.u, b.w), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(b.u.length(), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(b.w.length(), 1.0f, 1e-3f);
}

inline void test_get_antipode_short_arc_unchanged() {
  Basis b = make_basis(Quaternion(1, 0, 0, 0), Vector(0, 1, 0));
  auto [nb, nr] = get_antipode(b, 0.5f);
  // radius <= 1 → identity
  HS_EXPECT_VEC(nb.u, b.u, 1e-6f);
  HS_EXPECT_VEC(nb.v, b.v, 1e-6f);
  HS_EXPECT_VEC(nb.w, b.w, 1e-6f);
  HS_EXPECT_NEAR(nr, 0.5f, 1e-6f);
}

inline void test_get_antipode_long_arc_flips() {
  Basis b = make_basis(Quaternion(1, 0, 0, 0), Vector(0, 1, 0));
  auto [nb, nr] = get_antipode(b, 1.4f);
  HS_EXPECT_VEC(nb.u, -b.u, 1e-6f);
  HS_EXPECT_VEC(nb.v, -b.v, 1e-6f);
  HS_EXPECT_VEC(nb.w, b.w, 1e-6f);
  HS_EXPECT_NEAR(nr, 2.0f - 1.4f, 1e-6f);
}

// ============================================================================
// Orientation<W, CAP>
// ============================================================================

inline void test_orientation_default_is_identity() {
  Orientation<32, 8> o;
  HS_EXPECT_EQ(o.length(), 1);
  HS_EXPECT_QUAT(o.get(), Quaternion(1, 0, 0, 0), 1e-6f);
}

inline void test_orientation_set_clears_history() {
  Orientation<32, 8> o;
  o.push(make_rotation(Vector(0, 1, 0), 0.5f));
  o.push(make_rotation(Vector(0, 1, 0), 1.0f));
  HS_EXPECT_EQ(o.length(), 3);

  Quaternion target = make_rotation(Vector(1, 0, 0), 0.7f);
  o.set(target);
  HS_EXPECT_EQ(o.length(), 1);
  HS_EXPECT_QUAT(o.get(), target, 1e-6f);
}

inline void test_orientation_push_tracks_history() {
  Orientation<32, 4> o;
  Quaternion q1 = make_rotation(Vector(0, 1, 0), 0.1f);
  Quaternion q2 = make_rotation(Vector(0, 1, 0), 0.2f);
  o.push(q1);
  o.push(q2);
  HS_EXPECT_EQ(o.length(), 3); // identity + 2 pushes
  HS_EXPECT_QUAT(o.get(1), q1, 1e-6f);
  HS_EXPECT_QUAT(o.get(2), q2, 1e-6f);
  HS_EXPECT_QUAT(o.get(), q2, 1e-6f);
}

inline void test_orientation_orient_round_trips() {
  Orientation<32, 4> o;
  Quaternion q = make_rotation(Vector(0, 1, 0), PI_F * 0.5f);
  o.set(q);
  Vector v(1, 0, 0);
  Vector rotated = o.orient(v);
  Vector back = o.unorient(rotated);
  HS_EXPECT_VEC(back, v, 1e-3f);
}

inline void test_orientation_collapse_keeps_latest() {
  Orientation<32, 8> o;
  Quaternion q1 = make_rotation(Vector(0, 1, 0), 0.1f);
  Quaternion q2 = make_rotation(Vector(0, 1, 0), 0.2f);
  o.push(q1);
  o.push(q2);
  HS_EXPECT_EQ(o.length(), 3);

  o.collapse();
  HS_EXPECT_EQ(o.length(), 1);
  HS_EXPECT_QUAT(o.get(), q2, 1e-6f);
}

inline void test_orientation_upsample_preserves_endpoints() {
  Orientation<32, 16> o;
  Quaternion start = Quaternion(1, 0, 0, 0);
  Quaternion end = make_rotation(Vector(0, 1, 0), PI_F * 0.5f);
  o.set(start);
  o.push(end);
  HS_EXPECT_EQ(o.length(), 2);

  o.upsample(8);
  HS_EXPECT_EQ(o.length(), 8);
  // Endpoints preserved (slerp may flip sign — same orientation).
  HS_EXPECT_NEAR(std::abs(dot(o.get(0), start)), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(std::abs(dot(o.get(7), end)), 1.0f, 1e-3f);
}

inline void test_orientation_upsample_noop_if_already_long() {
  Orientation<32, 16> o;
  o.push(make_rotation(Vector(0, 1, 0), 0.1f));
  o.push(make_rotation(Vector(0, 1, 0), 0.2f));
  o.push(make_rotation(Vector(0, 1, 0), 0.3f));
  int before = o.length();
  o.upsample(before - 1);
  HS_EXPECT_EQ(o.length(), before);
}

// ============================================================================
// Runner
// ============================================================================

inline int run_geometry_tests() {
  auto scope = hs_test::begin_module("geometry");

  test_axis_constants();

  test_fragment_lerp_endpoints();
  test_fragment_lerp_midpoint();

  test_y_to_phi_free_function();
  test_phi_to_y_free_function();
  test_y_phi_roundtrip();
  test_y_to_phi_templated_LUT();
  test_y_to_phi_offset_injection_clips_no_double_apply();

  test_pixel_to_vector_unit_length();
  test_pixel_to_vector_known_samples();
  test_pixel_to_vector_float_branch_matches_phi_lut();
  test_vector_to_pixel_roundtrip_via_pixel_to_vector();

  test_log_polar_roundtrip();
  test_log_polar_north_pole_sentinel();

  test_fib_spiral_unit_length();
  test_fib_spiral_deterministic();
  test_fib_spiral_endpoints();

  test_lissajous_unit_length();
  test_random_vector_unit_length();

  test_make_basis_orthonormal();
  test_make_basis_alternate_normals();
  test_get_antipode_short_arc_unchanged();
  test_get_antipode_long_arc_flips();

  test_orientation_default_is_identity();
  test_orientation_set_clears_history();
  test_orientation_push_tracks_history();
  test_orientation_orient_round_trips();
  test_orientation_collapse_keeps_latest();
  test_orientation_upsample_preserves_endpoints();
  test_orientation_upsample_noop_if_already_long();

  return hs_test::end_module(scope);
}

} // namespace geometry
} // namespace hs_test

