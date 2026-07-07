/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/engine/transformers.h — the pure geometry transform functions
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
 * @brief Verifies the identity Mobius map round-trips a point through
 *        gnomonic projection back to itself.
 */
inline void test_gnomonic_mobius_identity_roundtrip() {
  MobiusParams id;
  Vector v = Vector(0.3f, 0.7f, 0.2f).normalized();
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
  Vector v = Vector(0.3f, 0.7f, 0.2f).normalized();
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
  p.phase = PI_F * 0.5f; // wavelet peak at d == 90°
  p.decay = 0.0f;
  p.thickness = 1.0f;
  // Default thresholds (min=1, max=-1) disable the fast reject.

  Vector v = Vector(1, 0, 0); // 90° from center → at the wavelet peak
  Vector r = ripple_transform(v, p);

  HS_EXPECT_TRUE(finite_vec(r));
  HS_EXPECT_NEAR(r.length(), 1.0f, 1e-4f);
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
  p.decay = 0.0f;
  p.thickness = 0.4f;    // band ≈ [phase - 0.4, phase + 0.4] rad
  p.prepare_thresholds();

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
  RippleParams base;
  base.center = Vector(0, 1, 0);
  base.amplitude = 0.5f;
  base.phase = PI_F * 0.5f; // wavelet peak at d == 90°
  base.thickness = 1.0f;
  // Default thresholds (min=1, max=-1) disable the fast reject.
  const Vector v = Vector(1, 0, 0);

  auto moved = [&](float decay) {
    RippleParams p = base;
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
 * @brief Pins the prepare_thresholds() reject band at its edges, not just well
 *        outside them.
 * @details test_ripple_threshold_reject_path samples points 0.3 rad beyond each
 *          edge, so a band misplaced by up to ~0.3 rad would still pass. Here a
 *          point a hair inside each edge takes the slow path (the wavelet tail is
 *          tiny but nonzero, so it moves), while a point a hair outside is
 *          fast-rejected and returned bit-for-bit unchanged. The asymmetry
 *          (moves vs. exactly-equal) brackets the edge to within the epsilon.
 */
inline void test_ripple_threshold_boundary() {
  RippleParams p;
  p.center = Vector(0, 1, 0); // north pole → cos_d == dot(v, center) == v.y
  p.amplitude = 0.5f;
  p.phase = PI_F * 0.5f;
  p.decay = 0.0f;
  p.thickness = 0.4f; // band edges at phase ± thickness
  p.prepare_thresholds();

  const float d_min = p.phase - p.thickness;
  const float d_max = p.phase + p.thickness;
  const float eps = 0.02f;

  auto pt = [](float d) { return Vector(std::sin(d), std::cos(d), 0.0f); };
  auto moved = [&](const Vector &src) {
    Vector r = ripple_transform(src, p);
    return std::abs(r.x - src.x) + std::abs(r.y - src.y) + std::abs(r.z - src.z);
  };

  // Just inside each edge: wavelet tail is small but nonzero.
  HS_EXPECT_GT(moved(pt(d_max - eps)), 1e-4f);
  HS_EXPECT_GT(moved(pt(d_min + eps)), 1e-4f);
  // Just outside each edge: fast-rejected → returned exactly unchanged.
  HS_EXPECT_NEAR(moved(pt(d_max + eps)), 0.0f, 1e-7f);
  HS_EXPECT_NEAR(moved(pt(d_min - eps)), 0.0f, 1e-7f);
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
    RippleParams rp;
    rp.amplitude = 0.0f;
    HS_EXPECT_TRUE(vec_bits_equal(ripple_transform(v, rp), v));
    NoiseParams np;
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
 *        unit sphere, and actually displaces the points.
 * @details The displacement assertion is the companion to
 *          test_noise_zero_amplitude_is_identity: without it, a noise_transform
 *          that returned the input unchanged would still pass the finite /
 *          on-sphere checks. Displacement is summed across samples so a single
 *          noise zero-crossing at one point cannot make the test flaky.
 */
inline void test_noise_active_stays_on_sphere() {
  NoiseParams p;
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

  HS_EXPECT_TRUE(nt.spawn_pinned(0) != nullptr);
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
// Transformer<> non-pinned slot recycling survives timeline compaction
// ============================================================================

/**
 * @brief Minimal 8x8 effect backing a stand-in Canvas for stepping a Timeline.
 * @details A fresh effect reports buffer_free() true, so the Canvas ctor does
 * not spin on a display ISR the host never runs.
 */
struct RecycleFakeEffect : public Effect {
  RecycleFakeEffect() : Effect(8, 8) {}
  void draw_frame() override {}
  bool strobe_columns() const override { return false; }
};

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
  RecycleFakeEffect fx;
  Canvas cv(fx);

  RippleTransformer<1> rt(tl);

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
  test_gnomonic_mobius_identity_roundtrip();
  test_gnomonic_mobius_known_rotation();
  test_ripple_zero_amplitude_is_identity();
  test_ripple_center_point_is_identity();
  test_ripple_active_rotates_on_sphere();
  test_ripple_threshold_reject_path();
  test_ripple_decay_attenuates();
  test_ripple_threshold_boundary();
  test_transforms_nonfinite_passes_through_identity();
  test_noise_zero_amplitude_is_identity();
  test_noise_active_stays_on_sphere();
  test_transformer_no_entities_is_identity();
  test_transformer_spawn_applies_and_composes();
  test_transformer_nonpinned_slot_reclaimed_after_compaction();

  return fixture.result();
}

} // namespace transformers_tests
} // namespace hs_test

