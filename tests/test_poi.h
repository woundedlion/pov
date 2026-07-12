/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for the DisplacementField poi phase:
 *   - random_rotation      : unit length, seeded determinism, no directional bias.
 *   - poi_field            : zero outside support; sign equals the push direction
 *                            inside; continuity at the three cap edges;
 *                            field_bound() bounds a dense cap sweep;
 *                            support_bound() bounds the support.
 *   - prefilter exactness  : a poi past a ring's support reach but within its
 *                            field_bound distance fields exactly 0 on that ring.
 *   - move table invariant : every ordered move pair's crossfade keeps
 *                            (1+lerp(b))·lerp(Σr) ≤ 1.
 *   - PoiChoreography      : antipodal pair symmetry, |P| ≤ 1 (ρ clamp never
 *                            active) across a full move cycle, per-frame center
 *                            continuity across a move transition, non-recurrence.
 *   - DisplacementField    : the phase machine cycles BALLS → NOISE → POI →
 *                            NOISE → BALLS, spawning all 12 dancers on POI entry
 *                            and leaving POI only once they all reclaim.
 */
#pragma once

#include "core/engine/engine.h"
#include "effects/DisplacementField.h"
#include "core/render/canvas.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cmath>
#include <cstdint>

namespace hs_test {
namespace poi_tests {

/**
 * @brief White-box accessor for DisplacementField's private phase state.
 */
struct PoiWhiteBox {
  template <int W, int H>
  static int phase(const DisplacementField<W, H> &e) {
    return static_cast<int>(e.phase);
  }
  template <int W, int H>
  static int poi_count(const DisplacementField<W, H> &e) {
    return e.pois.active_count();
  }
  template <int W, int H>
  static int ball_count(const DisplacementField<W, H> &e) {
    return e.balls.active_count();
  }
  template <int W, int H>
  static void set_rings(DisplacementField<W, H> &e, float n) {
    e.params.num_rings = n;
  }
};

/**
 * @brief Minimal 8x8 effect backing a Canvas for stepping poi animations.
 */
struct PoiFakeEffect : public Effect {
  PoiFakeEffect() : Effect(8, 8) {}
  void draw_frame() override {}
};

// ============================================================================
// random_rotation
// ============================================================================

/**
 * @brief Verifies random_rotation() returns unit quaternions and reproduces its
 *        sequence exactly under a reseeded RNG.
 */
inline void test_random_rotation_unit_and_deterministic() {
  hs::random().seed(4242u);
  Quaternion a = random_rotation();
  Quaternion b = random_rotation();
  HS_EXPECT_NEAR(a.squared_magnitude(), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(b.squared_magnitude(), 1.0f, 1e-5f);
  float spread = std::fabs(a.r - b.r) + std::fabs(a.v.x - b.v.x) +
                 std::fabs(a.v.y - b.v.y) + std::fabs(a.v.z - b.v.z);
  HS_EXPECT_GT(spread, 1e-3f);

  hs::random().seed(4242u);
  Quaternion a2 = random_rotation();
  HS_EXPECT_NEAR(a.r, a2.r, 1e-6f);
  HS_EXPECT_NEAR(a.v.x, a2.v.x, 1e-6f);
  HS_EXPECT_NEAR(a.v.y, a2.v.y, 1e-6f);
  HS_EXPECT_NEAR(a.v.z, a2.v.z, 1e-6f);
}

/**
 * @brief Verifies random_rotation() has no gross directional bias: the mean of
 *        many rotated +X directions is near the origin.
 */
inline void test_random_rotation_covers_sphere() {
  hs::random().seed(7u);
  Vector mean(0, 0, 0);
  const int N = 3000;
  for (int i = 0; i < N; ++i)
    mean = mean + rotate(Vector(1, 0, 0), random_rotation());
  mean = mean * (1.0f / N);
  HS_EXPECT_LT(mean.length(), 0.1f);
}

// ============================================================================
// poi_field
// ============================================================================

/**
 * @brief Builds a poi field centered on the +y meridian.
 * @param R Footprint radius.
 * @param amp Drape gain.
 * @param sigma Push sign.
 * @return Params ready for poi_field() (envelope 1, threshold prepared).
 */
inline PoiParams make_poi(float R, float amp, float sigma) {
  const float c_lat = PI_F / 2.0f;
  PoiParams p;
  p.center = Vector(std::sin(c_lat), std::cos(c_lat), 0.0f);
  p.axis = Vector(0, 1, 0);
  p.radius = R;
  p.amplitude = amp;
  p.direction = sigma;
  p.envelope = 1.0f;
  p.prepare_threshold();
  return p;
}

/**
 * @brief Samples poi_field at a colatitude / azimuth about the +y axis.
 */
inline float field_at(const PoiParams &p, float lat, float az) {
  Vector v(std::sin(lat) * std::cos(az), std::cos(lat),
           std::sin(lat) * std::sin(az));
  return poi_field(v, p);
}

/**
 * @brief Verifies the field is exactly 0 outside the cap and shares the push
 *        sign everywhere inside it.
 */
inline void test_poi_field_support_and_sign() {
  const float c_lat = PI_F / 2.0f;
  for (float sigma : {1.0f, -1.0f}) {
    PoiParams p = make_poi(0.5f, 1.0f, sigma);
    // Outside the cap along the meridian: exactly 0.
    HS_EXPECT_NEAR(field_at(p, c_lat + 0.5f + 0.05f, 0.0f), 0.0f, 1e-7f);
    HS_EXPECT_NEAR(field_at(p, c_lat - 0.5f - 0.05f, 0.0f), 0.0f, 1e-7f);
    // Far away (past the cos_radius reject): exactly 0.
    HS_EXPECT_NEAR(field_at(p, c_lat, PI_F), 0.0f, 1e-7f);

    // Every interior sample carries the push sign (field·sigma ≥ 0), and at
    // least one is strictly signed.
    float max_signed = 0.0f;
    for (float dy = -0.45f; dy <= 0.45f; dy += 0.05f)
      for (float az = -0.4f; az <= 0.4f; az += 0.1f) {
        float f = field_at(p, c_lat + dy, az);
        HS_EXPECT_GT(f * sigma, -1e-4f);
        max_signed = std::max(max_signed, f * sigma);
      }
    HS_EXPECT_GT(max_signed, 1e-2f);
  }
}

/**
 * @brief Verifies the field is continuous at the three cap edges: a sample just
 *        inside each edge is small and one just outside is exactly 0.
 */
inline void test_poi_field_edge_continuity() {
  const float c_lat = PI_F / 2.0f;
  const float R = 0.5f;
  PoiParams p = make_poi(R, 1.0f, 1.0f); // sigma = +1: downstream is +y'
  const float d = 0.01f;

  // Downstream edge (y' = arc): meridian, dy = +R. Field → 0.
  HS_EXPECT_LT(std::fabs(field_at(p, c_lat + R - d, 0.0f)), 0.05f);
  HS_EXPECT_NEAR(field_at(p, c_lat + R + d, 0.0f), 0.0f, 1e-7f);

  // Upstream edge (w = 0): meridian, dy = -R. Field → 0.
  HS_EXPECT_LT(std::fabs(field_at(p, c_lat - R + d, 0.0f)), 0.05f);
  HS_EXPECT_NEAR(field_at(p, c_lat - R - d, 0.0f), 0.0f, 1e-7f);

  // Lateral edge (arc → 0): same colatitude as the center, azimuth at the cap
  // boundary d ≈ R (sin(c_lat) = 1 here, so az ≈ R).
  const float az_edge = R;
  HS_EXPECT_LT(std::fabs(field_at(p, c_lat, az_edge - d)), 0.05f);
  HS_EXPECT_NEAR(field_at(p, c_lat, az_edge + d), 0.0f, 1e-7f);
}

/**
 * @brief Verifies field_bound() upper-bounds |poi_field| over a dense cap sweep,
 *        in both the unsaturated and saturated regimes, and support_bound()
 *        bounds the field's support.
 */
inline void test_poi_field_bounds() {
  const float c_lat = PI_F / 2.0f;
  const float R = 0.5f;
  for (float amp : {0.5f, 1.0f, 4.0f}) {
    PoiParams p = make_poi(R, amp, 1.0f);
    float bound = p.field_bound();
    float support = p.support_bound();
    float max_abs = 0.0f;
    for (float dy = -R - 0.05f; dy <= R + 0.05f; dy += 0.01f)
      for (float az = -R - 0.05f; az <= R + 0.05f; az += 0.01f) {
        float f = field_at(p, c_lat + dy, az);
        max_abs = std::max(max_abs, std::fabs(f));
        // Beyond the support radius the field is exactly 0.
        float d = std::sqrt(dy * dy + az * az);
        if (d > support + 0.02f)
          HS_EXPECT_NEAR(f, 0.0f, 1e-7f);
      }
    HS_EXPECT_GT(bound, max_abs - 1e-4f);
    // The bound is meant to be tight, not vacuous.
    HS_EXPECT_GT(max_abs, 0.3f * bound);
  }
}

// ============================================================================
// Prefilter exactness
// ============================================================================

/**
 * @brief Verifies a poi whose center lies past a ring's support reach but within
 *        its field_bound distance fields exactly 0 on that ring — the exactness
 *        the subset field_dominant relies on.
 */
inline void test_poi_prefilter_exactness() {
  const float c_lat = PI_F / 2.0f;
  const float R = 0.4f;
  // Saturated: field_bound = 2·r_eff exceeds support_bound = r_eff.
  PoiParams p = make_poi(R, 4.0f, 1.0f);
  const float support = p.support_bound();
  const float bound = p.field_bound();

  // Ring 1.5·R below the center: outside the support reach, inside field_bound.
  const float gap = 1.5f * R;
  HS_EXPECT_GT(gap, support);
  HS_EXPECT_GT(bound, gap);

  for (float az = 0.0f; az < 2.0f * PI_F; az += 0.1f)
    HS_EXPECT_NEAR(field_at(p, c_lat + gap, az), 0.0f, 1e-7f);
}

// ============================================================================
// Move table invariant
// ============================================================================

/**
 * @brief Verifies every ordered move pair keeps (1+lerp(b))·lerp(Σr) ≤ 1 across
 *        the crossfade (closed-form max of the quadratic), and each move alone
 *        stays within the 0.98 table headroom.
 */
inline void test_poi_move_table_invariant() {
  using C = Animation::PoiChoreography;
  auto sum_r = [](int m) {
    return C::MOVES[m].r[0] + C::MOVES[m].r[1] + C::MOVES[m].r[2];
  };
  for (int i = 0; i < C::NUM_MOVES; ++i) {
    HS_EXPECT_LT((1.0f + C::MOVES[i].b) * sum_r(i), 0.981f);
    for (int k = 0; k < C::NUM_MOVES; ++k) {
      float bi = C::MOVES[i].b, bk = C::MOVES[k].b;
      float si = sum_r(i), sk = sum_r(k);
      // g(t) = (1 + bi + (bk-bi)t)(si + (sk-si)t) = A t^2 + B t + C.
      float A = (bk - bi) * (sk - si);
      float B = (1.0f + bi) * (sk - si) + (bk - bi) * si;
      float Cc = (1.0f + bi) * si;
      auto g = [&](float t) { return (A * t + B) * t + Cc; };
      float m = std::max(g(0.0f), g(1.0f));
      if (std::fabs(A) > 1e-9f) {
        float tv = -B / (2.0f * A);
        if (tv > 0.0f && tv < 1.0f)
          m = std::max(m, g(tv));
      }
      HS_EXPECT_LT(m, 1.0f + 1e-4f);
    }
  }
}

/**
 * @brief Verifies every move has a velocity-dominant first harmonic and a
 *        center that never reaches the epicycle origin, so |P| and the travel
 *        speed stay bounded away from zero — no cusps, so the velocity-derived
 *        push axis rotates smoothly rather than snapping.
 * @details A two-term epitrochoid's speed is |Σ n_j r_j e^{iφ}|, minimized at
 *        ||n1|r1 − Σ_{j≥2}|n_j|r_j|; its radius |P| bottoms out at
 *        |r1 − Σ_{j≥2} r_j|. Both must clear a margin for the drape direction to
 *        vary continuously at any dance speed.
 */
inline void test_poi_move_no_cusp() {
  using C = Animation::PoiChoreography;
  for (int i = 0; i < C::NUM_MOVES; ++i) {
    const auto &m = C::MOVES[i];
    float vel_min = std::fabs(m.n[0]) * m.r[0];
    float pos_min = m.r[0];
    for (int j = 1; j < C::NUM_TERMS; ++j) {
      vel_min -= std::fabs(m.n[j]) * m.r[j];
      pos_min -= m.r[j];
    }
    HS_EXPECT_GT(vel_min, 0.2f); // no near-cusp: speed never collapses
    HS_EXPECT_GT(pos_min, 0.2f); // center never reaches the epicycle origin
  }
}

// ============================================================================
// PoiChoreography
// ============================================================================

/**
 * @brief Builds a choreography for a fresh phase with a fixed seed.
 */
inline void begin_choreo(Animation::PoiChoreography &c, uint32_t seed) {
  hs::random().seed(seed);
  Vector primaries[Animation::PoiChoreography::NUM_PAIRS];
  static constexpr int PRIMARY[6] = {0, 1, 4, 5, 8, 9};
  for (int p = 0; p < Animation::PoiChoreography::NUM_PAIRS; ++p)
    primaries[p] = Solids::Icosahedron::vertices[PRIMARY[p]];
  c.begin(random_rotation(), primaries);
}

/**
 * @brief Verifies antipodal partners trace the exact point reflection:
 *        C_{k+6}(t) = −C_k(t) every frame.
 */
inline void test_poi_choreo_pair_symmetry() {
  PoiFakeEffect fx;
  Canvas cv(fx);
  Animation::PoiChoreography c;
  begin_choreo(c, 11u);

  float canon[3];
  c.compute_canon(0, 1.0f, canon);

  Orientation<> ori;
  PoiParams p0, p6;
  Animation::PoiDance d0(p0, c, ori, 0, canon, 1000);
  Animation::PoiDance d6(p6, c, ori, 6, canon, 1000);

  for (int f = 0; f < 400; ++f) {
    c.step();
    d0.step(cv);
    d6.step(cv);
    HS_EXPECT_NEAR(p0.center.x, -p6.center.x, 1e-5f);
    HS_EXPECT_NEAR(p0.center.y, -p6.center.y, 1e-5f);
    HS_EXPECT_NEAR(p0.center.z, -p6.center.z, 1e-5f);
  }
}

/**
 * @brief Verifies |P| ≤ 1 across a full move cycle including crossfades and
 *        breathe, so the ρ safety clamp never engages.
 */
inline void test_poi_choreo_p_magnitude_bounded() {
  Animation::PoiChoreography c;
  begin_choreo(c, 23u);
  // Two moves plus crossfades exercise every table transition boundary reached.
  for (int f = 0; f < 2 * Animation::PoiChoreography::MOVE_FRAMES; ++f) {
    c.step();
    float sum_r = c.radii[0] + c.radii[1] + c.radii[2];
    float bound = (1.0f + std::fabs(c.breathe_amp)) * sum_r;
    HS_EXPECT_LT(bound, 1.0f + 1e-4f);
  }
}

/**
 * @brief Verifies a dancer's center steps continuously (no jump) across a move
 *        transition and that state does not recur one move period later.
 */
inline void test_poi_choreo_continuity_and_nonrecurrence() {
  PoiFakeEffect fx;
  Canvas cv(fx);
  Animation::PoiChoreography c;
  begin_choreo(c, 31u);

  Orientation<> ori;
  float canon[3] = {0.0f, 0.0f, 0.0f};
  PoiParams p;
  Animation::PoiDance d(p, c, ori, 0, canon, 5000);

  Vector prev(0, 0, 0);
  Vector at_start(0, 0, 0);
  bool have_prev = false;
  const int period = Animation::PoiChoreography::MOVE_FRAMES;
  // Straddle the first move boundary (frame MOVE_FRAMES) to catch a phase pop.
  for (int f = 0; f <= period + 5; ++f) {
    c.step();
    d.step(cv);
    if (f == 5)
      at_start = p.center;
    if (have_prev) {
      float step = angle_between(prev, p.center);
      HS_EXPECT_LT(step, 0.1f);
    }
    prev = p.center;
    have_prev = true;
  }
  // One move period after frame 5 the state has drifted (precession + breathe).
  HS_EXPECT_GT(angle_between(at_start, p.center), 1e-3f);
}

/**
 * @brief Verifies that under an identity (still) ring frame the push axis is the
 *        dome's world direction of travel: a unit tangent aligned with the
 *        per-frame center step.
 */
inline void test_poi_dance_axis_follows_velocity() {
  PoiFakeEffect fx;
  Canvas cv(fx);
  Animation::PoiChoreography c;
  begin_choreo(c, 41u);
  float speed = 5.0f;
  c.set_speed_source(&speed);

  Orientation<> ori; // identity: ring frame == world frame
  float canon[3] = {0.0f, 0.0f, 0.0f};
  PoiParams p;
  Animation::PoiDance d(p, c, ori, 0, canon, 5000);

  // Warm up so prev_center_ring is seeded and the motion is steady.
  Vector prev(0, 0, 0);
  for (int f = 0; f < 20; ++f) {
    c.step();
    d.step(cv);
    prev = p.center;
  }
  c.step();
  d.step(cv);
  Vector center = p.center;

  Vector vel = center - prev;
  vel = vel - dot(vel, center) * center;
  HS_EXPECT_GT(vel.length(), 1e-4f); // moving, so the axis is velocity-derived
  vel = vel.normalized();

  HS_EXPECT_NEAR(p.axis.length(), 1.0f, 1e-4f);
  HS_EXPECT_NEAR(dot(p.axis, center), 0.0f, 1e-4f);  // tangent to the sphere
  HS_EXPECT_GT(dot(p.axis, vel), 0.999f);            // aligned with travel
}

/**
 * @brief Verifies the push axis is the net motion against the ring frame: a
 *        world-fixed dome whose rings sweep past still drapes, in the direction
 *        the observer sees it move relative to the (rotating) rings.
 * @details Freezes the dance (Dance Speed 0 on the b=0 Orbit move, so the center
 *        is truly stationary in world space) while spinning the ring orientation
 *        about +Y. A world-fixed point C seen from a frame rotating about Y has
 *        apparent tangential motion along normalize(C × Y); the axis must track
 *        it, not collapse to the stationary-dome fallback.
 */
inline void test_poi_dance_axis_against_ring_frame() {
  PoiFakeEffect fx;
  Canvas cv(fx);
  Animation::PoiChoreography c;
  begin_choreo(c, 61u);
  float speed = 0.0f; // frozen dance
  c.set_speed_source(&speed);
  // Force the Orbit move (b = 0) so breathe is constant and the dome truly
  // stops — otherwise breathe would still inch |P| and move the center.
  c.order[0] = 0;
  c.move_pos = 0;
  c.move_frame = 0;
  for (int j = 0; j < 3; ++j) {
    c.n_from[j] = Animation::PoiChoreography::MOVES[0].n[j];
    c.r_from[j] = Animation::PoiChoreography::MOVES[0].r[j];
    c.radii[j] = Animation::PoiChoreography::MOVES[0].r[j];
  }
  c.b_from = 0.0f;

  Orientation<> ori;
  float canon[3] = {0.0f, 0.0f, 0.0f};
  PoiParams p;
  Animation::PoiDance d(p, c, ori, 0, canon, 20000);

  const Vector Y(0, 1, 0);
  for (int f = 0; f < 40; ++f) {
    c.step();
    Animation::Rotation<96>::animate(cv, ori, Y, 0.03f, ease_linear);
    d.step(cv);
  }
  const Vector center = p.center;
  Vector expected = cross(center, Y);
  HS_EXPECT_GT(expected.length(), 1e-2f); // center not on the spin axis
  expected = expected.normalized();

  HS_EXPECT_NEAR(p.axis.length(), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(dot(p.axis, center), 0.0f, 1e-3f);      // tangent to the sphere
  HS_EXPECT_GT(std::fabs(dot(p.axis, expected)), 0.98f); // along the apparent sweep
}

// ============================================================================
// DisplacementField phase machine
// ============================================================================

/**
 * @brief Resets the process-global effect state to a clean baseline.
 */
inline void reset_globals() {
  hs::random().seed(1337u);
  configure_arenas_default();
  Timeline().clear();
}

/**
 * @brief Verifies the effect opens on POI and cycles POI → NOISE → BALLS →
 *        NOISE → POI, spawns all 12 dancers on POI entry, and leaves POI only
 *        once they all reclaim their slots.
 * @details Driven at one ring so each frame is cheap over the multi-thousand
 *        frame cycle; the phase durations are unchanged.
 */
inline void test_poi_effect_phase_cycle() {
  reset_globals();
  auto *e = new DisplacementField<96, 20>();
  e->init();
  PoiWhiteBox::set_rings(*e, 1.0f);

  const int BALLS = 0, NOISE = 1, POI = 2;
  HS_EXPECT_EQ(PoiWhiteBox::phase(*e), POI);

  int seq[8] = {POI};
  int seq_len = 1;
  int max_poi = 0;
  bool poi_drained_before_exit = true;
  int prev = POI;

  const int MAX_FRAMES = 9000;
  int f = 0;
  for (; f < MAX_FRAMES && seq_len < 5; ++f) {
    e->draw_frame();
    e->advance_display();
    int ph = PoiWhiteBox::phase(*e);
    if (ph == POI)
      max_poi = std::max(max_poi, PoiWhiteBox::poi_count(*e));
    // POI must drain to zero dancers before it can hand off to NOISE.
    if (prev == POI && ph == NOISE && PoiWhiteBox::poi_count(*e) != 0)
      poi_drained_before_exit = false;
    if (ph != prev) {
      if (seq_len < 8)
        seq[seq_len++] = ph;
      prev = ph;
    }
  }

  HS_EXPECT_TRUE(seq_len >= 5);
  HS_EXPECT_EQ(seq[0], POI);
  HS_EXPECT_EQ(seq[1], NOISE);
  HS_EXPECT_EQ(seq[2], BALLS);
  HS_EXPECT_EQ(seq[3], NOISE);
  HS_EXPECT_EQ(seq[4], POI);
  HS_EXPECT_EQ(max_poi, 12);
  HS_EXPECT_TRUE(poi_drained_before_exit);

  delete e;
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every poi test case.
 * @return The module's failure count.
 */
inline int run_poi_tests() {
  hs_test::ModuleFixture fixture("poi");

  test_random_rotation_unit_and_deterministic();
  test_random_rotation_covers_sphere();
  test_poi_field_support_and_sign();
  test_poi_field_edge_continuity();
  test_poi_field_bounds();
  test_poi_prefilter_exactness();
  test_poi_move_table_invariant();
  test_poi_move_no_cusp();
  test_poi_choreo_pair_symmetry();
  test_poi_choreo_p_magnitude_bounded();
  test_poi_choreo_continuity_and_nonrecurrence();
  test_poi_dance_axis_follows_velocity();
  test_poi_dance_axis_against_ring_frame();
  test_poi_effect_phase_cycle();

  return fixture.result();
}

} // namespace poi_tests
} // namespace hs_test
