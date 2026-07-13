/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/math/geometry.h.
 *
 * Coverage:
 *   - Axis constants (X/Y/Z_AXIS, UP)
 *   - Fragment::lerp
 *   - y_to_phi / phi_to_y conversions (free + templated)
 *   - PhiLUT / TrigLUT initialisation and lookup values
 *   - pixel_to_vector / vector_to_pixel roundtrip
 *   - fib_spiral, lissajous, random_vector (all on unit sphere)
 *   - Basis: make_basis (orthonormality), get_antipode (flipping)
 *   - Orientation<CAP>: set, push, get, collapse, upsample, orient/unorient
 */
#pragma once

#include "core/math/geometry.h"
#include "tests/test_3dmath.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace geometry {

using hs_test::math3d::approx_vec;
using hs_test::math3d::approx_quat;

// ============================================================================
// Axis constants
// ============================================================================

/**
 * @brief Verifies X/Y/Z_AXIS are the canonical unit axes and UP aliases Y_AXIS.
 */
inline void test_axis_constants() {
  HS_EXPECT_VEC(X_AXIS, Vector(1, 0, 0), 0.0f);
  HS_EXPECT_VEC(Y_AXIS, Vector(0, 1, 0), 0.0f);
  HS_EXPECT_VEC(Z_AXIS, Vector(0, 0, 1), 0.0f);
  HS_EXPECT_VEC(UP, Y_AXIS, 0.0f);
}

// ============================================================================
// Fragment::lerp
// ============================================================================

/**
 * @brief Verifies lerp at t=0 returns a and at t=1 returns b for pos/scalar/age
 *        fields.
 */
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

/**
 * @brief Verifies lerp at t=0.5 averages position and scalar fields
 *        component-wise.
 */
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

/**
 * @brief Verifies y_to_phi maps row 0→0, top row→PI, midpoint→PI/2.
 * @details The mapping is phi = y*PI/(h_virt-1).
 */
inline void test_y_to_phi_free_function() {
  HS_EXPECT_NEAR(y_to_phi(0.0f, 145), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(y_to_phi(144.0f, 145), PI_F, 1e-5f);
  HS_EXPECT_NEAR(y_to_phi(72.0f, 145), PI_F * 0.5f, 1e-5f);
}

/**
 * @brief Verifies phi_to_y is the inverse of y_to_phi: 0→0, PI→top row,
 *        PI/2→midpoint.
 */
inline void test_phi_to_y_free_function() {
  HS_EXPECT_NEAR(phi_to_y(0.0f, 145), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(phi_to_y(PI_F, 145), 144.0f, 1e-5f);
  HS_EXPECT_NEAR(phi_to_y(PI_F * 0.5f, 145), 72.0f, 1e-5f);
}

/**
 * @brief Pins the south-pole row-index contract: phi=PI can land a float hair
 *        *above* H_VIRT-1, but the floored row a caller indexes with must stay
 *        in [0, H_VIRT-1].
 * @details The proximity (HS_EXPECT_NEAR) checks pass even at H_VIRT-1+epsilon
 *          and do not pin the ceiling. The raw float overshoots for some heights
 *          (e.g. h_virt=64); flooring is the documented caller fix, so the bound
 *          is asserted on static_cast<int>(y), across the free function, the
 *          templated phi_to_y<H>, and vector_to_pixel of the south pole.
 */
inline void test_phi_to_y_south_pole_row_in_bounds() {
  for (int h_virt : {16, 64, 145, 23}) {
    int row = static_cast<int>(phi_to_y(PI_F, h_virt));
    HS_EXPECT_GE(row, 0);
    HS_EXPECT_LE(row, h_virt - 1);
  }

  constexpr int H = 32;
  constexpr int H_VIRT = H + hs::H_OFFSET;
  int row_t = static_cast<int>(phi_to_y<H>(PI_F));
  HS_EXPECT_GE(row_t, 0);
  HS_EXPECT_LE(row_t, H_VIRT - 1);

  PixelCoords p = vector_to_pixel<H, H>(Vector(0, -1, 0));
  int row_v = static_cast<int>(p.y);
  HS_EXPECT_GE(row_v, 0);
  HS_EXPECT_LE(row_v, H_VIRT - 1);
}

/**
 * @brief Verifies phi_to_y(y_to_phi(i)) recovers i for every row across several
 *        heights.
 */
inline void test_y_phi_roundtrip() {
  for (int h_virt : {16, 64, 145}) {
    for (int i = 0; i < h_virt; ++i) {
      float phi = y_to_phi(static_cast<float>(i), h_virt);
      float y = phi_to_y(phi, h_virt);
      HS_EXPECT_NEAR(y, static_cast<float>(i), 1e-4f);
    }
  }
}

/**
 * @brief Verifies the templated y_to_phi<H> (integer and float overloads) uses
 *        the PhiLUT and agrees with the offset-aware free function for every row.
 */
inline void test_y_to_phi_templated_LUT() {
  constexpr int H = 32;
  for (int y = 0; y < H; ++y) {
    float lut = y_to_phi<H>(y);
    float direct = y_to_phi(static_cast<float>(y), H + hs::H_OFFSET);
    HS_EXPECT_NEAR(lut, direct, 1e-5f);
  }

  HS_EXPECT_TRUE(PhiLUT<H>::initialized);

  float lut_f = y_to_phi<H>(5.0f);
  HS_EXPECT_NEAR(lut_f, y_to_phi<H>(5), 1e-5f);
}

/**
 * @brief Pins the offset semantics of the phi mapping by injecting a non-zero
 *        offset directly through the free function, so the formula holds even
 *        though the native build has H_OFFSET==0.
 * @details The free function is the building block the templated path mirrors.
 *          H_OFFSET adds virtual rows so the image is CLIPPED (not stretched) at
 *          the bottom of the ring, where the LEDs stop short of the south pole.
 *          The mapping is phi = y*PI/(H_VIRT-1) with H_VIRT = H + H_OFFSET, so
 *          the bottom physical row y=H-1 lands short of PI while the virtual
 *          bottom row reaches PI exactly. The guard also asserts the offset is
 *          applied once, not twice (dividing by H + 2*OFF - 1 would push the
 *          bottom row even further from the pole).
 */
inline void test_y_to_phi_offset_injection_clips_no_double_apply() {
  constexpr int H = 20;          // Holosphere hardware height
  constexpr int OFF = 3;         // hardware H_OFFSET
  constexpr int H_VIRT = H + OFF; // 23
  const float bottom_phys = static_cast<float>(H - 1); // y = 19

  float correct = y_to_phi(bottom_phys, H_VIRT);          // 19*PI/22
  HS_EXPECT_NEAR(correct, bottom_phys * PI_F / (H_VIRT - 1), 1e-5f);
  HS_EXPECT_TRUE(correct < PI_F);

  // A double-applied offset would divide by (H + 2*OFF - 1) = 25.
  float double_applied = bottom_phys * PI_F / (H + 2 * OFF - 1); // 19*PI/24
  HS_EXPECT_TRUE(std::abs(correct - double_applied) > 1e-2f);

  HS_EXPECT_NEAR(y_to_phi(static_cast<float>(H_VIRT - 1), H_VIRT), PI_F, 1e-5f);
}

// ============================================================================
// TrigLUT / pixel_to_vector / vector_to_pixel
// ============================================================================

/**
 * @brief Verifies pixel_to_vector returns unit vectors across the full virtual
 *        grid, and that using it lazily initialises the TrigLUT.
 */
inline void test_pixel_to_vector_unit_length() {
  constexpr int W = 32, H = 32;
  for (int y = 0; y < H + hs::H_OFFSET; ++y) {
    for (int x = 0; x < W; x += 4) {
      Vector v = pixel_to_vector<W, H>(x, y);
      HS_EXPECT_NEAR(v.length(), 1.0f, 1e-3f);
    }
  }
  HS_EXPECT_TRUE((TrigLUT<W, H>::initialized));
}

/**
 * @brief Verifies pixel_to_vector lands on the expected directions at the
 *        equator and north pole, the two anchors of the spherical mapping.
 */
inline void test_pixel_to_vector_known_samples() {
  constexpr int W = 32, H = 32;
  // Integer y rarely lands exactly on phi=π/2 (equator is y=(H_VIRT-1)/2).
  Vector v = pixel_to_vector<W, H>(0, (H + hs::H_OFFSET) / 2);
  HS_EXPECT_VEC(v, Vector(1, 0, 0), 0.15f);

  // y=0 (north pole): phi=0 ⇒ Vector(0, 1, 0).
  Vector north = pixel_to_vector<W, H>(0, 0);
  HS_EXPECT_VEC(north, Vector(0, 1, 0), 5e-2f);
}

/**
 * @brief Pins pixel_to_vector's fractional-y (float) branch to y_to_phi<H>, so
 *        the float and integer LUT paths cannot diverge on a build with non-zero
 *        H_OFFSET.
 * @details y_to_phi<H> is the same phi source the integer LUT path and every
 *          SDF/scan shape use. At x=0, Vector(Spherical(0, phi)) =
 *          (sin phi, cos phi, 0), so the recovered phi is acos(v.y).
 */
inline void test_pixel_to_vector_float_branch_matches_phi_lut() {
  constexpr int W = 32, H = 32;
  for (float y : {0.5f, 5.25f, 12.75f, 20.5f}) {
    Vector v = pixel_to_vector<W, H>(0.0f, y);
    float recovered_phi = std::acos(std::clamp(v.y, -1.0f, 1.0f));
    HS_EXPECT_NEAR(recovered_phi, y_to_phi<H>(y), 1e-4f);
  }
}

/**
 * @brief Verifies vector_to_pixel inverts pixel_to_vector to within half a pixel
 *        for non-degenerate samples.
 * @details Poles are excluded, where azimuth wrap is undefined.
 */
inline void test_vector_to_pixel_roundtrip_via_pixel_to_vector() {
  constexpr int W = 64, H = 64;
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
// fib_spiral
// ============================================================================

/**
 * @brief Verifies every fib_spiral sample lies on the unit sphere.
 */
inline void test_fib_spiral_unit_length() {
  for (int i = 0; i < 32; ++i) {
    Vector v = fib_spiral(32, 0.5f, i);
    HS_EXPECT_NEAR(v.length(), 1.0f, 1e-3f);
  }
}

/**
 * @brief Verifies fib_spiral is a pure function: same args give the same vector.
 */
inline void test_fib_spiral_deterministic() {
  Vector v1 = fib_spiral(64, 0.5f, 7);
  Vector v2 = fib_spiral(64, 0.5f, 7);
  HS_EXPECT_VEC(v1, v2, 1e-6f);
}

/**
 * @brief Verifies the first and last fib_spiral samples are unit-length and sit
 *        on opposite hemispheres for N>=4.
 * @details i=0 sits near one pole, i=N-1 near the other.
 */
inline void test_fib_spiral_endpoints() {
  Vector first = fib_spiral(16, 0.5f, 0);
  Vector last = fib_spiral(16, 0.5f, 15);
  HS_EXPECT_NEAR(first.length(), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(last.length(), 1.0f, 1e-3f);
  HS_EXPECT_TRUE(first.y * last.y < 0.0f);
}

/**
 * @brief Verifies consecutive fib_spiral samples advance in azimuth by a
 *        constant golden-angle step, the property that gives the spiral its
 *        even spread (a wrong-but-unit-length placement would fail here).
 * @details Sample i's azimuth is 2*pi*i*INV_PHI, so consecutive samples differ
 *        by a fixed golden step (mod 2*pi). Mid-range indices are used because
 *        near the poles sin(phi) -> 0 makes the recovered azimuth noisy.
 */
inline void test_fib_spiral_golden_angle() {
  const int n = 64;
  auto wrap_pi = [](float a) {
    while (a > PI_F) a -= 2.0f * PI_F;
    while (a <= -PI_F) a += 2.0f * PI_F;
    return a;
  };
  const float step = wrap_pi(2.0f * PI_F * INV_PHI);
  for (int i = 8; i < n - 8; ++i) {
    Vector a = fib_spiral(n, 0.5f, i);
    Vector b = fib_spiral(n, 0.5f, i + 1);
    float delta = wrap_pi(atan2f(b.z, b.x) - atan2f(a.z, a.x));
    HS_EXPECT_NEAR(delta, step, 1e-3f);
  }
}

// ============================================================================
// lissajous
// ============================================================================

/**
 * @brief Verifies every lissajous sample lies on the unit sphere across a full
 *        period.
 */
inline void test_lissajous_unit_length() {
  for (float t = 0.0f; t < 6.28f; t += 0.3f) {
    Vector v = lissajous(3, 2, 0.5f, t);
    HS_EXPECT_NEAR(v.length(), 1.0f, 1e-3f);
  }
}

/**
 * @brief Verifies the lissajous phase a is in RADIANS: it computes
 *        cos(m1*t - a), matching the designer's radians-labelled slider and raw
 *        export.
 * @details With m1=m2=1, a=1, t=1 the phase arg is 1-1=0, so the point is
 *          (sin1*cos0, cos1, sin1*sin0) = (sin1, cos1, 0). A scale-by-PI
 *          convention would put the arg at 1-PI and shift x/z sharply, which
 *          this rules out.
 */
inline void test_lissajous_phase_is_radians() {
  Vector v = lissajous(1.0f, 1.0f, 1.0f, 1.0f);
  HS_EXPECT_NEAR(v.x, sinf(1.0f), 1e-5f);
  HS_EXPECT_NEAR(v.y, cosf(1.0f), 1e-5f);
  HS_EXPECT_NEAR(v.z, 0.0f, 1e-5f);
}

// ============================================================================
// random_vector (Marsaglia)
// ============================================================================

/**
 * @brief Verifies random_vector (Marsaglia) always returns a unit-length
 *        direction.
 */
inline void test_random_vector_unit_length() {
  for (int i = 0; i < 32; ++i) {
    Vector v = random_vector();
    HS_EXPECT_NEAR(v.length(), 1.0f, 1e-3f);
  }
}

/**
 * @brief random_vector is a pure function of the shared RNG stream: reseeding
 *        the generator reproduces the exact same sequence bit-for-bit.
 * @details The unit-length test alone would also pass if random_vector pulled
 *          from a nondeterministic source (time/entropy). Pinning reproducibility
 *          guards the sim↔device determinism contract documented on hs::random().
 *          The global generator is saved and restored so this test can't perturb
 *          the stream position other RNG-touching tests observe.
 */
inline void test_random_vector_deterministic() {
  auto saved = hs::random();
  constexpr int N = 16;
  Vector first[N];
  hs::random().seed(1337);
  for (int i = 0; i < N; ++i) first[i] = random_vector();
  hs::random().seed(1337);
  for (int i = 0; i < N; ++i) {
    Vector v = random_vector();
    HS_EXPECT_VEC(v, first[i], 0.0f);
  }
  hs::random() = saved;
}

/**
 * @brief random_vector is approximately uniform on the sphere — not biased
 *        toward an axis or hemisphere.
 * @details A degenerate generator that always returned the same direction (or a
 *          biased one) would still pass the unit-length check. Over many samples
 *          the per-axis means must sit near zero, each axis must split roughly
 *          evenly across its sign, and pooled absolute components must occupy
 *          [0,1] uniformly. Bounds are deliberately loose (the per-axis
 *          mean's standard error at N=4000 is ~0.009, so ±0.08 is ~9σ; the sign
 *          split's σ is ~32 counts, so ±400 is ~12σ) — a true uniform generator
 *          passes with vast margin while a stuck/biased one fails. Seeded for
 *          reproducibility; the generator is saved and restored.
 */
inline void test_random_vector_distribution() {
  auto saved = hs::random();
  hs::random().seed(0xC0FFEEu);
  constexpr int N = 4000;
  constexpr int BIN_COUNT = 10;
  constexpr int MIN_BIN_SIZE = 3 * N * 3 / (4 * BIN_COUNT);
  constexpr int MAX_BIN_SIZE = 3 * N * 5 / (4 * BIN_COUNT);
  float sx = 0.0f, sy = 0.0f, sz = 0.0f;
  int posx = 0, posy = 0, posz = 0;
  int abs_bins[BIN_COUNT] = {};
  for (int i = 0; i < N; ++i) {
    Vector v = random_vector();
    HS_EXPECT_NEAR(v.length(), 1.0f, 1e-3f);
    sx += v.x; sy += v.y; sz += v.z;
    if (v.x > 0.0f) ++posx;
    if (v.y > 0.0f) ++posy;
    if (v.z > 0.0f) ++posz;
    for (float component : {v.x, v.y, v.z}) {
      const int bin = std::min(static_cast<int>(fabsf(component) * BIN_COUNT),
                               BIN_COUNT - 1);
      ++abs_bins[bin];
    }
  }
  HS_EXPECT_NEAR(sx / N, 0.0f, 0.08f);
  HS_EXPECT_NEAR(sy / N, 0.0f, 0.08f);
  HS_EXPECT_NEAR(sz / N, 0.0f, 0.08f);
  HS_EXPECT_TRUE(posx > N * 4 / 10 && posx < N * 6 / 10);
  HS_EXPECT_TRUE(posy > N * 4 / 10 && posy < N * 6 / 10);
  HS_EXPECT_TRUE(posz > N * 4 / 10 && posz < N * 6 / 10);
  for (int count : abs_bins) {
    HS_EXPECT_GE(count, MIN_BIN_SIZE);
    HS_EXPECT_LE(count, MAX_BIN_SIZE);
  }
  hs::random() = saved;
}

// ============================================================================
// Basis: make_basis / get_antipode
// ============================================================================

/**
 * @brief Verifies make_basis produces an orthonormal frame whose v axis is the
 *        given normal.
 */
inline void test_make_basis_orthonormal() {
  Quaternion id(1, 0, 0, 0);
  Basis b = make_basis(id, Vector(0, 1, 0));
  HS_EXPECT_VEC(b.v, Vector(0, 1, 0), 1e-4f);
  HS_EXPECT_NEAR(b.u.length(), 1.0f, 1e-4f);
  HS_EXPECT_NEAR(b.w.length(), 1.0f, 1e-4f);
  HS_EXPECT_NEAR(dot(b.u, b.v), 0.0f, 1e-4f);
  HS_EXPECT_NEAR(dot(b.v, b.w), 0.0f, 1e-4f);
  HS_EXPECT_NEAR(dot(b.u, b.w), 0.0f, 1e-4f);
}

/**
 * @brief Verifies make_basis stays orthonormal and v-aligned for an off-axis
 *        tilted normal.
 */
inline void test_make_basis_alternate_normals() {
  Quaternion id(1, 0, 0, 0);
  Vector n = Vector(0.5f, 0.5f, 0.7071f).normalized();
  Basis b = make_basis(id, n);
  HS_EXPECT_VEC(b.v, n, 1e-3f);
  HS_EXPECT_NEAR(dot(b.u, b.v), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(dot(b.u, b.w), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(b.u.length(), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(b.w.length(), 1.0f, 1e-3f);
}

/**
 * @brief Verifies get_antipode leaves the basis and radius untouched when
 *        radius <= 1.
 */
inline void test_get_antipode_short_arc_unchanged() {
  Basis b = make_basis(Quaternion(1, 0, 0, 0), Vector(0, 1, 0));
  auto [nb, nr] = get_antipode(b, 0.5f);
  HS_EXPECT_VEC(nb.u, b.u, 1e-6f);
  HS_EXPECT_VEC(nb.v, b.v, 1e-6f);
  HS_EXPECT_VEC(nb.w, b.w, 1e-6f);
  HS_EXPECT_NEAR(nr, 0.5f, 1e-6f);
}

/**
 * @brief Verifies that for radius > 1, get_antipode flips u and v, keeps w, and
 *        returns 2 - radius.
 */
inline void test_get_antipode_long_arc_flips() {
  Basis b = make_basis(Quaternion(1, 0, 0, 0), Vector(0, 1, 0));
  auto [nb, nr] = get_antipode(b, 1.4f);
  HS_EXPECT_VEC(nb.u, -b.u, 1e-6f);
  HS_EXPECT_VEC(nb.v, -b.v, 1e-6f);
  HS_EXPECT_VEC(nb.w, b.w, 1e-6f);
  HS_EXPECT_NEAR(nr, 2.0f - 1.4f, 1e-6f);
}

// ============================================================================
// Orientation<CAP>
// ============================================================================

/**
 * @brief Verifies a fresh Orientation holds a single identity rotation.
 */
inline void test_orientation_default_is_identity() {
  Orientation<8> o;
  HS_EXPECT_EQ(o.length(), static_cast<size_t>(1));
  HS_EXPECT_QUAT(o.get(), Quaternion(1, 0, 0, 0), 1e-6f);
}

/**
 * @brief Verifies set() discards accumulated history, leaving only the new
 *        orientation.
 */
inline void test_orientation_set_clears_history() {
  Orientation<8> o;
  o.push(make_rotation(Vector(0, 1, 0), 0.5f));
  o.push(make_rotation(Vector(0, 1, 0), 1.0f));
  HS_EXPECT_EQ(o.length(), static_cast<size_t>(3));

  Quaternion target = make_rotation(Vector(1, 0, 0), 0.7f);
  o.set(target);
  HS_EXPECT_EQ(o.length(), static_cast<size_t>(1));
  HS_EXPECT_QUAT(o.get(), target, 1e-6f);
}

/**
 * @brief Verifies push() appends to the history; get(i) indexes it and get()
 *        returns the latest.
 */
inline void test_orientation_push_tracks_history() {
  Orientation<4> o;
  Quaternion q1 = make_rotation(Vector(0, 1, 0), 0.1f);
  Quaternion q2 = make_rotation(Vector(0, 1, 0), 0.2f);
  o.push(q1);
  o.push(q2);
  HS_EXPECT_EQ(o.length(), static_cast<size_t>(3)); // identity + 2 pushes
  HS_EXPECT_QUAT(o.get(1), q1, 1e-6f);
  HS_EXPECT_QUAT(o.get(2), q2, 1e-6f);
  HS_EXPECT_QUAT(o.get(), q2, 1e-6f);
}

/**
 * @brief Verifies unorient inverts orient, recovering the original vector.
 */
inline void test_orientation_orient_round_trips() {
  Orientation<4> o;
  Quaternion q = make_rotation(Vector(0, 1, 0), PI_F * 0.5f);
  o.set(q);
  Vector v(1, 0, 0);
  Vector rotated = o.orient(v);
  Vector back = o.unorient(rotated);
  HS_EXPECT_VEC(back, v, 1e-3f);
}

/**
 * @brief Verifies collapse() reduces history to just the most recent
 *        orientation.
 */
inline void test_orientation_collapse_keeps_latest() {
  Orientation<8> o;
  Quaternion q1 = make_rotation(Vector(0, 1, 0), 0.1f);
  Quaternion q2 = make_rotation(Vector(0, 1, 0), 0.2f);
  o.push(q1);
  o.push(q2);
  HS_EXPECT_EQ(o.length(), static_cast<size_t>(3));

  o.collapse();
  HS_EXPECT_EQ(o.length(), static_cast<size_t>(1));
  HS_EXPECT_QUAT(o.get(), q2, 1e-6f);
}

/**
 * @brief Verifies upsample(n) slerp-interpolates to n entries while preserving
 *        both endpoints.
 */
inline void test_orientation_upsample_preserves_endpoints() {
  Orientation<16> o;
  Quaternion start = Quaternion(1, 0, 0, 0);
  Quaternion end = make_rotation(Vector(0, 1, 0), PI_F * 0.5f);
  o.set(start);
  o.push(end);
  HS_EXPECT_EQ(o.length(), static_cast<size_t>(2));

  o.upsample(8);
  HS_EXPECT_EQ(o.length(), static_cast<size_t>(8));
  // slerp may flip sign — same orientation.
  HS_EXPECT_NEAR(std::abs(dot(o.get(0), start)), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(std::abs(dot(o.get(7), end)), 1.0f, 1e-3f);
}

/**
 * @brief Verifies upsample is a no-op when the requested length is not greater
 *        than current.
 */
inline void test_orientation_upsample_noop_if_already_long() {
  Orientation<16> o;
  o.push(make_rotation(Vector(0, 1, 0), 0.1f));
  o.push(make_rotation(Vector(0, 1, 0), 0.2f));
  o.push(make_rotation(Vector(0, 1, 0), 0.3f));
  size_t before = o.length();
  o.upsample(static_cast<int>(before) - 1);
  HS_EXPECT_EQ(o.length(), before);
}

/**
 * @brief Verifies upsample clamps to CAPACITY when count exceeds it (the
 *        soft-degrade path) while keeping endpoints exact and writes in-bounds.
 */
inline void test_orientation_upsample_clamps_past_capacity() {
  Orientation<4> o;
  Quaternion start = Quaternion(1, 0, 0, 0);
  Quaternion end = make_rotation(Vector(0, 1, 0), PI_F * 0.5f);
  o.set(start);
  o.push(end);

  o.upsample(100);
  HS_EXPECT_EQ(o.length(), static_cast<size_t>(4));
  // slerp may flip sign — same orientation.
  HS_EXPECT_NEAR(std::abs(dot(o.get(0), start)), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(std::abs(dot(o.get(3), end)), 1.0f, 1e-3f);
}

/**
 * @brief Verifies the by-index accessors: orient(v,i)/unorient(v,i) match
 *        rotate(v, get(i)) and its conjugate, and at(i) aliases the stored
 *        quaternion.
 */
inline void test_orientation_historical_index_accessors() {
  Orientation<8> o;
  Quaternion q0 = make_rotation(Vector(1, 0, 0), 0.3f);
  Quaternion q1 = make_rotation(Vector(0, 1, 0), 0.7f);
  Quaternion q2 = make_rotation(Vector(0, 0, 1), 1.1f);
  o.set(q0);
  o.push(q1);
  o.push(q2);
  HS_EXPECT_EQ(o.length(), static_cast<size_t>(3));

  const Vector v(0.2f, 0.4f, 0.9f);
  for (int i = 0; i < 3; ++i) {
    HS_EXPECT_VEC(o.orient(v, i), rotate(v, o.get(i)), 1e-5f);
    HS_EXPECT_VEC(o.unorient(v, i), rotate(v, o.get(i).conjugate()), 1e-5f);
  }

  // at(i) aliases the stored quaternion: a write through it is seen by get(i).
  Quaternion q_new = make_rotation(Vector(0, 1, 0), 0.25f);
  o.at(1) = q_new;
  HS_EXPECT_QUAT(o.get(1), q_new, 1e-6f);
  HS_EXPECT_VEC(o.orient(v, 1), rotate(v, q_new), 1e-5f);
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every geometry test case.
 * @return The module's failure count.
 */
inline int run_geometry_tests() {
  hs_test::ModuleFixture fixture("geometry");

  test_axis_constants();

  test_fragment_lerp_endpoints();
  test_fragment_lerp_midpoint();

  test_y_to_phi_free_function();
  test_phi_to_y_free_function();
  test_phi_to_y_south_pole_row_in_bounds();
  test_y_phi_roundtrip();
  test_y_to_phi_templated_LUT();
  test_y_to_phi_offset_injection_clips_no_double_apply();

  test_pixel_to_vector_unit_length();
  test_pixel_to_vector_known_samples();
  test_pixel_to_vector_float_branch_matches_phi_lut();
  test_vector_to_pixel_roundtrip_via_pixel_to_vector();

  test_fib_spiral_unit_length();
  test_fib_spiral_deterministic();
  test_fib_spiral_endpoints();
  test_fib_spiral_golden_angle();

  test_lissajous_unit_length();
  test_lissajous_phase_is_radians();
  test_random_vector_unit_length();
  test_random_vector_deterministic();
  test_random_vector_distribution();

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
  test_orientation_upsample_clamps_past_capacity();
  test_orientation_historical_index_accessors();

  return fixture.result();
}

} // namespace geometry
} // namespace hs_test

