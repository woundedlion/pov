/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Ripple test cases. */

#include "tests/transformers/fixture.h"

namespace hs_test::transformers_tests {

/**
 * @brief Verifies a point coincident with the ripple center is returned
 *        unchanged.
 * @details Such a point has a degenerate rotation axis (cross(center, center)
 *          == 0), so the transform must short-circuit to the identity.
 */
inline void test_ripple_center_point_is_identity() {
  Animation::RippleParams p;
  p.center = math::Vector(0, 1, 0);
  p.amplitude = 0.8f;
  p.phase = 0.0f;
  p.decay = 0.0f;
  math::Vector v = math::Vector(0, 1, 0);
  math::Vector r = ripple_transform(v, p);
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
  p.center = math::Vector(0, 1, 0);
  p.amplitude = 0.5f;
  p.phase = math::PI_F * 0.5f; // wavelet peak at d == 90°
  p.decay = 0.0f;
  p.thickness = 1.0f;
  // Default cos_near=1, cos_far=-1 disable the fast reject.

  math::Vector v =
      math::Vector(1, 0, 0); // 90° from center → at the wavelet peak
  math::Vector r = ripple_transform(v, p);

  HS_EXPECT_TRUE(finite_vec(r));
  HS_EXPECT_NEAR(r.length(), 1.0f, 1e-4f);
  float moved = displacement(r, v);
  HS_EXPECT_GT(moved, 1e-2f);
}

/**
 * @brief Exercises the sync() fast-reject band.
 * @details The in-band point rotates; points nearer than d_min or farther
 *          than d_max return unchanged.
 */
inline void test_ripple_threshold_reject_path() {
  Animation::RippleParams p;
  p.center =
      math::Vector(0, 1, 0); // north pole, so cos_d == dot(v, center) == v.y
  p.amplitude = 0.5f;
  p.phase = math::PI_F * 0.5f; // wavelet peak at d == 90°
  p.decay = 0.0f;
  p.thickness = 0.4f; // band ≈ [phase - 0.4, phase + 0.4] rad
  p.sync();

  HS_EXPECT_LT(p.cos_near, 1.0f);
  HS_EXPECT_GT(p.cos_far, -1.0f);

  // In-band: at the peak.
  const math::Vector IN_BAND = math::Vector(1, 0, 0);
  const math::Vector r_in = ripple_transform(IN_BAND, p);
  HS_EXPECT_TRUE(finite_vec(r_in));
  HS_EXPECT_NEAR(r_in.length(), 1.0f, 1e-4f);
  const float moved_in = displacement(r_in, IN_BAND);
  HS_EXPECT_GT(moved_in, 1e-2f);

  // Out-of-band toward the center (d < d_min): cos_d > cos_near leg.
  const float d_near = p.phase - p.thickness - 0.05f;
  const math::Vector near_c =
      math::Vector(std::sin(d_near), std::cos(d_near), 0.0f);
  const math::Vector r_near = ripple_transform(near_c, p);
  HS_EXPECT_EQ(r_near.x, near_c.x);
  HS_EXPECT_EQ(r_near.y, near_c.y);
  HS_EXPECT_EQ(r_near.z, near_c.z);

  // Out-of-band away from the center (d > d_max): cos_d < cos_far leg.
  const float d_far = p.phase + p.thickness + 0.05f;
  const math::Vector far_c =
      math::Vector(std::sin(d_far), std::cos(d_far), 0.0f);
  const math::Vector r_far = ripple_transform(far_c, p);
  HS_EXPECT_EQ(r_far.x, far_c.x);
  HS_EXPECT_EQ(r_far.y, far_c.y);
  HS_EXPECT_EQ(r_far.z, far_c.z);
}

/**
 * @brief Verifies the decay parameter attenuates the ripple's rotation with
 *        angular distance from the center.
 * @details The same peak point (d == 90° from the center) is transformed at
 *          increasing decay rates: each must move it strictly less, and a
 *          strong decay collapses the rotation to ~identity.
 */
inline void test_ripple_decay_attenuates() {
  Animation::RippleParams base;
  base.center = math::Vector(0, 1, 0);
  base.amplitude = 0.5f;
  base.phase = math::PI_F * 0.5f; // wavelet peak at d == 90°
  base.thickness = 1.0f;
  // Default cos_near=1, cos_far=-1 disable the fast reject.
  const math::Vector V = math::Vector(1, 0, 0);

  auto moved = [&](float decay) {
    Animation::RippleParams p = base;
    p.decay = decay;
    math::Vector r = ripple_transform(V, p);
    return displacement(r, V);
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
 * @brief Pins the sync() reject band at its edges.
 * @details A point a hair inside each edge takes the slow path (the wavelet
 *          tail is tiny but nonzero, so it moves), while a point a hair outside
 *          is fast-rejected and returned bit-for-bit unchanged.
 */
inline void test_ripple_threshold_boundary() {
  Animation::RippleParams p;
  p.center =
      math::Vector(0, 1, 0); // north pole → cos_d == dot(v, center) == v.y
  p.amplitude = 0.5f;
  p.phase = math::PI_F * 0.5f;
  p.decay = 0.0f;
  p.thickness = 0.4f; // band edges at phase ± thickness
  p.sync();

  const float d_min = p.phase - p.thickness;
  const float d_max = p.phase + p.thickness;
  const float EPS = 0.02f;

  auto pt = [](float d) {
    return math::Vector(std::sin(d), std::cos(d), 0.0f);
  };
  auto moved = [&](const math::Vector &src) {
    math::Vector r = ripple_transform(src, p);
    return displacement(r, src);
  };

  // Just inside each edge: wavelet tail is small but nonzero.
  HS_EXPECT_GT(moved(pt(d_max - EPS)), 1e-4f);
  HS_EXPECT_GT(moved(pt(d_min + EPS)), 1e-4f);
  // Just outside each edge: fast-rejected → returned exactly unchanged.
  HS_EXPECT_EQ(moved(pt(d_max + EPS)), 0.0f);
  HS_EXPECT_EQ(moved(pt(d_min - EPS)), 0.0f);
}

/** @brief Neutral endpoint groups preserve mixed-center spawn order. */
inline void test_ripple_kernel_endpoint_order() {
  using namespace RippleKernel;
  static_assert(std::is_same_v<Animation::RippleWave, Wave>);
  std::array<Wave, 5> waves;
  for (size_t index = 0; index < waves.size(); ++index) {
    auto &wave = waves[index];
    wave.center = index == 3 ? math::X_AXIS : math::Y_AXIS;
    wave.amplitude = index == 2 ? 0.0f : 0.07f - 0.006f * index;
    wave.phase = 0.7f + 0.19f * index;
    wave.decay = 0.2f;
    wave.thickness = 1.2f;
    wave.sync();
  }
  std::array<EndpointGroup<Wave>, 5> groups;
  int count = 0, first = 0;
  for (const auto &wave : waves)
    if (wave.amplitude > 0.001f)
      append_endpoint(groups.data(), count, first, wave);
  HS_EXPECT_EQ(count, 4);
  HS_EXPECT_EQ(groups[0].end, 2);
  HS_EXPECT_EQ(groups[2].end, 3);
  HS_EXPECT_EQ(groups[3].end, 4);
  HS_EXPECT_TRUE(groups[0].params == &waves[0]);
  HS_EXPECT_TRUE(groups[1].params == &waves[1]);
  HS_EXPECT_TRUE(groups[2].params == &waves[3]);
  HS_EXPECT_TRUE(groups[3].params == &waves[4]);
  size_t visited_groups = 0, visited_waves = 0;
  for (int polar = 1; polar < 20; ++polar)
    for (int azimuth = 0; azimuth < 23; ++azimuth) {
      const auto input = math::Vector::from_spherical(
          azimuth * math::TWO_PI_F / 23, polar * math::PI_F / 20);
      auto expected = input;
      for (const auto &wave : waves)
        expected = ripple_transform(expected, wave);
      const auto actual = sequential_groups(
          input, groups.data(), count, [&] { ++visited_groups; },
          [&] { ++visited_waves; });
      HS_EXPECT_TRUE(vec_identical(actual, expected));
    }
  HS_EXPECT_EQ(visited_groups, size_t{19 * 23 * 3});
  HS_EXPECT_GT(visited_waves, size_t{0});
  HS_EXPECT_LT(visited_waves, size_t{19 * 23 * 4});
}

/** @brief Support bands include both cutoffs and form a union. */
inline void test_ripple_kernel_support_bands() {
  using namespace RippleKernel;
  Wave a, b;
  a.cos_near = 0.2f;
  a.cos_far = -0.3f;
  b.cos_near = 0.8f;
  b.cos_far = -0.6f;
  HS_EXPECT_TRUE(contains(a, a.cos_near));
  HS_EXPECT_TRUE(contains(a, a.cos_far));
  HS_EXPECT_FALSE(contains(a, std::nextafter(a.cos_near, 1.0f)));
  HS_EXPECT_FALSE(contains(a, std::nextafter(a.cos_far, -1.0f)));
  float near = -1.0f, far = 1.0f;
  include_support(near, far, a);
  include_support(near, far, b);
  HS_EXPECT_EQ(near, 0.8f);
  HS_EXPECT_EQ(far, -0.6f);
}

/** @brief Ripple angle follows its amplitude, decay and thickness. */
inline void test_ripple_kernel_angle_at_distance() {
  using namespace RippleKernel;
  Wave a;
  a.cos_near = .2f;
  a.cos_far = -.3f;
  a.amplitude = 0.13f;
  a.decay = 0.4f;
  HS_EXPECT_NEAR(angle_at_distance(a, 1.25f, 1.25f), 0.13f * std::exp(-0.5f),
                 2e-6f);
  a.thickness = 1.0f;
  HS_EXPECT_EQ(angle_at_distance(a, 1.25f, 1.0f), 0.0f);
}

/** @brief Transfer sampling accepts covered paths and rejects invalid nodes. */
inline void test_ripple_kernel_sample_transfer() {
  using namespace RippleKernel;
  const std::array<RadialTransfer, 3> PREVIOUS{
      {{0.1f, 0.2f}, {0.2f, 0.4f}, {0.3f, 0.6f}}};
  auto current = PREVIOUS;
  for (auto &node : current) {
    node.angle += 0.2f;
    node.path += 0.4f;
  }
  RadialTransfer sample;
  HS_EXPECT_TRUE(sample_transfer(PREVIOUS, current, math::PI_F * 0.25f, 0.25f,
                                 true, sample));
  HS_EXPECT_NEAR(sample.angle, 0.2f, 1e-7f);
  HS_EXPECT_NEAR(sample.path, 0.4f, 1e-7f);
  current[1].path = -1.0f;
  HS_EXPECT_FALSE(sample_transfer(PREVIOUS, current, math::PI_F * 0.25f, 0.25f,
                                  false, sample));
}

/**
 * @brief Verifies the truncated small-angle quaternion series reproduces the
 *        exact libm rotation it stands in for, and that the two branches meet
 *        at RIPPLE_SMALL_ANGLE_MAX.
 * @details theta is proportional to amplitude for a fixed geometry, so one
 *          exact-branch reading calibrates the envelope.
 */
inline void test_ripple_small_angle_series_matches_exact() {
  Animation::RippleParams p;
  p.center = math::Vector(0, 1, 0);
  p.phase = math::PI_F * 0.5f; // wavelet peak at d == 90 degrees
  p.decay = 0.0f;
  p.thickness = 1.0f;
  // Default cos_near=1, cos_far=-1 disable the fast reject.
  const math::Vector V = math::Vector(1, 0, 0);
  const math::Vector axis = math::cross(p.center, V).normalized();

  auto at = [&](float amplitude) {
    Animation::RippleParams q = p;
    q.amplitude = amplitude;
    return ripple_transform(V, q);
  };
  auto exact = [&](float theta) {
    return math::rotate(V, math::make_rotation(axis, theta));
  };

  // Calibrate theta/amplitude from the exact branch, in libm acos rather than
  // the fast_acos angle_between uses.
  auto turned = [&](const math::Vector &out) {
    return (float)std::acos(hs::clamp(math::dot(V, out), -1.0f, 1.0f));
  };
  const float ANCHOR_AMPLITUDE = 0.9f;
  const math::Vector anchor = at(ANCHOR_AMPLITUDE);
  const float envelope = turned(anchor) / ANCHOR_AMPLITUDE;
  // The peak of a zero-decay Ricker wavelet is unit gain.
  HS_EXPECT_NEAR(envelope, 1.0f, 2e-3f);
  HS_EXPECT_LE(max_component_delta(anchor, exact(ANCHOR_AMPLITUDE * envelope)),
               1e-6f);

  // Series branch, from just above the amplitude short-circuit up to the
  // switchover.
  const float AMP_LO = 0.002f;
  const float amp_hi = RIPPLE_SMALL_ANGLE_MAX / envelope;
  float worst_err = 0.0f;
  float worst_unit_err = 0.0f;
  const int STEPS = 300;
  for (int i = 0; i <= STEPS; ++i) {
    const float amplitude =
        AMP_LO + (amp_hi - AMP_LO) * (float)i / (float)STEPS;
    const math::Vector got = at(amplitude);
    worst_err = fold_worst(
        worst_err, max_component_delta(got, exact(amplitude * envelope)));
    worst_unit_err = fold_worst(worst_unit_err, std::fabs(got.length() - 1.0f));
  }
  HS_EXPECT_LE(worst_err, 5e-7f);
  HS_EXPECT_LE(worst_unit_err, 5e-7f);

  // The two branches meet: a hair either side of the switchover agrees.
  const float seam = RIPPLE_SMALL_ANGLE_MAX / envelope;
  HS_EXPECT_LE(max_component_delta(at(seam * 0.999f), at(seam * 1.001f)),
               1e-3f);
}

/**
 * @brief Pins the sync() band at a phase past π.
 * @details The [0, π] clamp collapses both bounds onto cos(π), so every
 *          non-antipodal direction is fast-rejected.
 */
inline void test_ripple_threshold_collapses_past_pi() {
  Animation::RippleParams p;
  p.center = math::Vector(0, 1, 0);
  p.amplitude = 0.5f;
  p.decay = 0.0f;
  p.thickness = 0.4f;
  p.phase = 20.0f;
  p.sync();

  HS_EXPECT_NEAR(p.cos_near, -1.0f, 1e-6f);
  HS_EXPECT_NEAR(p.cos_far, -1.0f, 1e-6f);

  // The band admits nothing: every direction but the exact antipode falls
  // outside it, takes the fast reject, and comes back bit for bit.
  for (int i = 1; i < 64; ++i) {
    const float d = math::PI_F * static_cast<float>(i) / 64.0f;
    const math::Vector v(std::sin(d), std::cos(d), 0.0f);
    HS_EXPECT_TRUE(vec_identical(ripple_transform(v, p), v));
  }
}

/**
 * @brief Feeds non-finite (NaN/Inf) directions through the transforms' identity
 *        short-circuits and an active ripple, and confirms they pass the input
 *        through verbatim.
 * @details Active ripple returns a NaN input unchanged when its degenerate-axis
 * comparison fails.
 */
inline void test_transforms_nonfinite_passes_through_identity() {
  const float INF = HUGE_VALF;
  const float NAN_VALUE = NAN;
  const math::Vector BAD[] = {math::Vector(NAN_VALUE, 0, 0),
                              math::Vector(0, INF, 0),
                              math::Vector(INF, NAN_VALUE, -INF)};

  for (const math::Vector &v : BAD) {
    Animation::RippleParams rp;
    rp.amplitude = 0.0f;
    HS_EXPECT_TRUE(vec_identical(ripple_transform(v, rp), v));
    Animation::RippleParams active;
    active.center = math::Vector(0.0f, 1.0f, 0.0f);
    active.amplitude = 0.5f;
    active.phase = math::PI_F * 0.5f;
    active.thickness = 1.0f;
    active.decay = 0.0f;
    HS_EXPECT_TRUE(vec_identical(ripple_transform(v, active), v));
    Animation::NoiseParams np;
    np.amplitude = 0.0f;
    HS_EXPECT_TRUE(vec_identical(noise_transform(v, np), v));
    Timeline tl;
    RippleTransformer<4> rt(tl);
    HS_EXPECT_TRUE(vec_identical(rt.transform(v), v));
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
