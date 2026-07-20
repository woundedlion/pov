/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/animation/animation.h and core/math/easing.h.
 *
 * Scope: the PURE, non-render parts of the animation system. The animations
 * exercised here (Transition, Mutation, Lerp, Driver, Rotation) take a `Canvas&`
 * in step() but never dereference it — AnimationBase::step only increments the
 * frame counter and the derived steps tested here only touch their bound
 * float/subject/orientation. A single genuine Canvas over a tiny static Effect
 * is shared (see fake_canvas()), so the reference is always valid yet untouched.
 *
 * Timeline IS exercised here. Its step() only drives the canvas-agnostic
 * animations above, so the shared fake_canvas() suffices — no render stack is
 * pulled. Every Timeline shares one global event array (plus the
 * live-guard and frame/count cursors), so each test scopes its Timeline locally
 * (balancing the guard) and resets the global cursors when it pokes them.
 *
 * Coverage:
 *   - Path::get_point: empty guard, endpoints, end-of-path clamp, midpoint
 *   - AnimationBase duration==0 -> 1 coercion (no divide-by-zero)
 *   - Transition: start->end stepping, easing, quantize, done()
 *   - Mutation: applies f(easing(t))
 *   - Driver: per-frame increment + wrap
 *   - Lerp: type-erased lerp(start,target,t) driven by easing
 *   - OrientationTrail: record/get ordering (index 0 is OLDEST — see note)
 *   - Timeline: start-frame sequencing, one-shot removal + survivor compaction,
 *     repeating rewind, .then()-chained event addition, clear(), capacity guard,
 *     shared-orientation motion-blur composition
 *   - easing.h: ease(0)~0, ease(1)~1 for in/out variants; output finite
 */
#pragma once

#include <cstddef>
#include <limits>
#include <vector>
#include "core/animation/animation.h"
#include "core/render/canvas.h"
#include "core/math/easing.h"
#include "core/mesh/mesh.h" // PolyMesh, MeshOps::compile (MeshMorph fixtures)
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace animation_tests {

// ============================================================================
// Stand-in Canvas reference
// ----------------------------------------------------------------------------
// One genuine Canvas over a tiny Effect, shared across tests. The animations
// under test take a Canvas& but never dereference it. No test may construct a
// second Canvas on this shared FakeEffect: that would queue_frame() and spin the
// next ctor on a display ISR the host never runs.
// ============================================================================

/**
 * @brief Tiny 8x8 effect that draws nothing; backs the shared fake_canvas().
 * @details A fresh effect's buffer_free() is true, so a Canvas built over it
 * does not spin in its constructor.
 */
struct FakeEffect : public Effect {
  /**
   * @brief Constructs the stand-in effect at a fixed 8x8 resolution.
   */
  FakeEffect() : Effect(8, 8) {}
  /**
   * @brief Draws nothing; the animations under test never render.
   */
  void draw_frame() override {}
};

/**
 * @brief Storage for the module-scoped fake-canvas pointer.
 * @return Reference to the pointer run_animation_tests() sets to the shared
 * fixture for the module's duration and clears (fixture destroyed) on exit.
 */
inline Canvas *&fake_canvas_ptr() {
  static Canvas *p = nullptr;
  return p;
}

/**
 * @brief Returns the module-scoped shared Canvas reused across the tests.
 * @return Reference to the fixture owned by run_animation_tests(); valid for the
 * duration of the module but never dereferenced by the animations under test.
 */
inline Canvas &fake_canvas() { return *fake_canvas_ptr(); }

// ============================================================================
// Path::get_point
// ============================================================================

/**
 * @brief Verifies an empty Path returns the origin for any t (the no-points
 * guard).
 */
inline void test_path_empty_returns_origin() {
  Path<32> p;
  Vector v = p.get_point(0.0f);
  HS_EXPECT_NEAR(v.x, 0.0f, 0.0f);
  HS_EXPECT_NEAR(v.y, 0.0f, 0.0f);
  HS_EXPECT_NEAR(v.z, 0.0f, 0.0f);
  Vector v1 = p.get_point(1.0f);
  HS_EXPECT_NEAR(v1.x, 0.0f, 0.0f);
}

/**
 * @brief Verifies get_point hits the exact endpoints at t=0/1, clamps to back()
 * past t=1, and interpolates linearly between samples.
 */
inline void test_path_endpoints_and_clamp() {
  Path<32> p;
  // A straight ramp along X: f(s) = (s, 0, 0), 4 samples, x in [0,1].
  p.append_segment([](float s) { return Vector(s, 0.0f, 0.0f); }, 1.0f, 4,
                   ease_linear);

  Vector start = p.get_point(0.0f);
  HS_EXPECT_NEAR(start.x, 0.0f, 1e-5f);

  Vector end = p.get_point(1.0f);
  HS_EXPECT_NEAR(end.x, 1.0f, 1e-5f);

  Vector over = p.get_point(2.0f);
  HS_EXPECT_NEAR(over.x, end.x, 1e-5f);

  Vector mid = p.get_point(0.5f);
  HS_EXPECT_NEAR(mid.x, 0.5f, 1e-5f);
}

/**
 * @brief Verifies collapse() reduces the path to its last point; get_point then
 * returns it for any t.
 */
inline void test_path_collapse_keeps_last() {
  Path<32> p;
  p.append_segment([](float s) { return Vector(s, 0.0f, 0.0f); }, 1.0f, 4,
                   ease_linear);
  Vector last_before = p.get_point(1.0f);
  p.collapse();
  Vector a = p.get_point(0.0f);
  Vector b = p.get_point(1.0f);
  HS_EXPECT_NEAR(a.x, last_before.x, 1e-5f);
  HS_EXPECT_NEAR(b.x, last_before.x, 1e-5f);
}

// ============================================================================
// Transition
// ============================================================================

/**
 * @brief Verifies a Transition steps its bound float from start to target over
 * `duration` frames and reports done() only once it lands on the target.
 */
inline void test_transition_reaches_target_linear() {
  float v = 0.0f;
  const int duration = 10;
  Animation::Transition tr(v, 100.0f, duration, ease_linear);
  HS_EXPECT_FALSE(tr.done());
  for (int i = 0; i < duration; ++i) {
    tr.step(fake_canvas());
  }
  HS_EXPECT_TRUE(tr.done());
  HS_EXPECT_NEAR(v, 100.0f, 1e-3f);
}

/**
 * @brief Verifies `from` is captured from the live value at the first step and
 * a forward Transition advances monotonically toward the target.
 */
inline void test_transition_monotonic_and_starts_from_current() {
  float v = 5.0f;
  const int duration = 8;
  Animation::Transition tr(v, 25.0f, duration, ease_linear);
  float prev = v;
  bool monotonic = true;
  for (int i = 0; i < duration; ++i) {
    tr.step(fake_canvas());
    if (v < prev - 1e-4f)
      monotonic = false;
    prev = v;
  }
  HS_EXPECT_TRUE(monotonic);
  HS_EXPECT_NEAR(v, 25.0f, 1e-3f);
}

/**
 * @brief Verifies a zero-duration Transition coerces to 1 frame: no
 * divide-by-zero, and a single step completes it.
 */
inline void test_transition_duration_zero_no_divide_by_zero() {
  float v = 0.0f;
  Animation::Transition tr(v, 42.0f, 0, ease_linear);
  tr.step(fake_canvas());
  HS_EXPECT_TRUE(std::isfinite(v));
  HS_EXPECT_NEAR(v, 42.0f, 1e-3f);
  HS_EXPECT_TRUE(tr.done());
}

/**
 * @brief Verifies a quantized Transition floors the eased value, landing on an
 * integer.
 */
inline void test_transition_quantized_floors_result() {
  float v = 0.0f;
  const int duration = 4;
  Animation::Transition tr(v, 3.7f, duration, ease_linear, /*quantized=*/true);
  for (int i = 0; i < duration; ++i) {
    tr.step(fake_canvas());
  }
  // 3.7 floored to 3.0.
  HS_EXPECT_NEAR(v, 3.0f, 1e-5f);
}

/**
 * @brief Verifies a repeating Transition rewinds t to 0 at the end of each
 * cycle and retraverses the full 0 -> target ramp.
 * @details `from` is captured once (on the first step), so the second cycle
 * must traverse the full ramp again rather than freezing at the target.
 */
inline void test_transition_repeat_retraverses_each_cycle() {
  float v = 0.0f;
  const int duration = 4;
  Animation::Transition tr(v, 10.0f, duration, ease_linear, /*quantized=*/false,
                           /*repeat=*/true);
  for (int i = 0; i < duration; ++i)
    tr.step(fake_canvas());
  HS_EXPECT_NEAR(v, 10.0f, 1e-3f);

  tr.rewind();
  tr.step(fake_canvas());
  HS_EXPECT_LT(v, 10.0f);
  for (int i = 1; i < duration; ++i)
    tr.step(fake_canvas());
  HS_EXPECT_NEAR(v, 10.0f, 1e-3f);
}

// ============================================================================
// Mutation
// ============================================================================

/**
 * @brief Verifies each step Mutation writes f(easing(t/duration)) into its
 * bound float.
 */
inline void test_mutation_applies_function_of_eased_time() {
  float v = 0.0f;
  const int duration = 5;
  // f(e) = e * 10; at completion e==1 -> v==10.
  Animation::Mutation m(
      v, [](float e) { return e * 10.0f; }, duration, ease_linear);
  for (int i = 0; i < duration; ++i) {
    m.step(fake_canvas());
  }
  HS_EXPECT_NEAR(v, 10.0f, 1e-3f);
  HS_EXPECT_TRUE(m.done());
}

/**
 * @brief Verifies a zero-duration Mutation coerces to 1 frame and yields a
 * finite result.
 */
inline void test_mutation_duration_zero_finite() {
  float v = 0.0f;
  Animation::Mutation m(
      v, [](float e) { return e; }, 0, ease_linear);
  m.step(fake_canvas());
  HS_EXPECT_TRUE(std::isfinite(v));
  HS_EXPECT_NEAR(v, 1.0f, 1e-3f);
}

// ============================================================================
// Driver
// ============================================================================

/**
 * @brief Verifies a wrapping Driver adds its per-frame increment each step and
 * wraps back into [0,1) when it reaches 1.
 */
inline void test_driver_increments_and_wraps() {
  float v = 0.0f;
  Animation::Driver d(v, 0.25f, /*wrap=*/true);
  d.step(fake_canvas());
  HS_EXPECT_NEAR(v, 0.25f, 1e-5f);
  d.step(fake_canvas());
  d.step(fake_canvas());
  HS_EXPECT_NEAR(v, 0.75f, 1e-5f);
  d.step(fake_canvas()); // 1.0 -> wraps to 0.0
  HS_EXPECT_NEAR(v, 0.0f, 1e-5f);
  HS_EXPECT_NEAR(d.get_mutant(), 0.0f, 1e-5f);
}

/**
 * @brief Verifies a non-wrapping Driver accumulates its increment past 1
 * without bound.
 */
inline void test_driver_no_wrap_accumulates() {
  float v = 0.0f;
  Animation::Driver d(v, 0.5f, /*wrap=*/false);
  for (int i = 0; i < 5; ++i)
    d.step(fake_canvas());
  HS_EXPECT_NEAR(v, 2.5f, 1e-5f);
}

/**
 * @brief Verifies a live-bound Driver ignores a non-finite slider frame instead
 * of permanently poisoning the wrapped mutant via wrap_t(NaN).
 */
inline void test_driver_nan_source_does_not_poison() {
  float v = 0.0f;
  float slider = 0.25f;
  Animation::Driver d(v, &slider, /*scale=*/1.0f, /*wrap=*/true);
  d.step(fake_canvas());
  HS_EXPECT_NEAR(v, 0.25f, 1e-5f);

  slider = std::numeric_limits<float>::quiet_NaN();
  d.step(fake_canvas());
  HS_EXPECT_TRUE(std::isfinite(v)); // NaN frame kept the last good speed
  HS_EXPECT_NEAR(v, 0.5f, 1e-5f);

  slider = 0.25f;
  d.step(fake_canvas());
  HS_EXPECT_NEAR(v, 0.75f, 1e-5f);
}

// ============================================================================
// Lerp (type-erased subject.lerp(start, target, t))
// ============================================================================

namespace {
/**
 * @brief Minimal subject satisfying the lerp(start, target, t) interface
 * Animation::Lerp drives via type erasure.
 */
struct Lerpable {
  float value = 0.0f; /**< Interpolated scalar payload. */
  /**
   * @brief Linearly interpolates this subject's value between two endpoints.
   * @param a Start subject (source value).
   * @param b Target subject (destination value).
   * @param t Interpolation parameter in [0, 1].
   */
  void lerp(const Lerpable &a, const Lerpable &b, float t) {
    value = a.value + (b.value - a.value) * t;
  }
};
} // namespace

/**
 * @brief Verifies Lerp drives its subject from start to target over `duration`
 * frames, landing on the target and reporting done().
 */
inline void test_lerp_drives_subject_to_target() {
  Lerpable subject, start, target;
  start.value = 0.0f;
  target.value = 200.0f;
  const int duration = 10;
  Animation::Lerp l(subject, start, target, duration, ease_linear);
  for (int i = 0; i < duration; ++i) {
    l.step(fake_canvas());
  }
  HS_EXPECT_NEAR(subject.value, 200.0f, 1e-3f);
  HS_EXPECT_TRUE(l.done());
}

/**
 * @brief Verifies at the halfway eased progress the subject sits exactly
 * between start and target.
 */
inline void test_lerp_midpoint() {
  Lerpable subject, start, target;
  start.value = 10.0f;
  target.value = 20.0f;
  const int duration = 4;
  Animation::Lerp l(subject, start, target, duration, ease_linear);
  l.step(fake_canvas()); // t=1 -> progress 0.25 -> 12.5
  l.step(fake_canvas()); // t=2 -> progress 0.50 -> 15.0
  HS_EXPECT_NEAR(subject.value, 15.0f, 1e-3f);
}

// ============================================================================
// OrientationTrail
// ----------------------------------------------------------------------------
// Index 0 is the OLDEST snapshot, length()-1 the newest (matching the JS
// simulator's trail ordering).
// ============================================================================

/**
 * @brief Verifies recorded snapshots are ordered oldest-first: index 0 is the
 * first recorded, the last index the newest.
 */
inline void test_orientation_trail_index_zero_is_oldest() {
  hs::random().seed(1337);
  Animation::OrientationTrail<Orientation<4>, 8> trail;

  Orientation<4> a, b, c;
  a.set(make_rotation(Vector(0, 1, 0), 0.1f));
  b.set(make_rotation(Vector(0, 1, 0), 0.2f));
  c.set(make_rotation(Vector(0, 1, 0), 0.3f));

  trail.record(a);
  trail.record(b);
  trail.record(c);
  HS_EXPECT_EQ(trail.length(), static_cast<size_t>(3));

  // Index 0 == oldest (a); last index == newest (c).
  HS_EXPECT_NEAR(std::abs(dot(trail.get(0).get(), a.get())), 1.0f, 1e-4f);
  HS_EXPECT_NEAR(std::abs(dot(trail.get(2).get(), c.get())), 1.0f, 1e-4f);
}

/**
 * @brief Verifies expire() removes the oldest snapshot, leaving the survivor at
 * index 0.
 */
inline void test_orientation_trail_expire_drops_oldest() {
  hs::random().seed(1337);
  Animation::OrientationTrail<Orientation<4>, 8> trail;
  Orientation<4> a, b;
  a.set(make_rotation(Vector(0, 1, 0), 0.1f));
  b.set(make_rotation(Vector(0, 1, 0), 0.2f));
  trail.record(a);
  trail.record(b);
  trail.expire(); // removes the oldest (a)
  HS_EXPECT_EQ(trail.length(), static_cast<size_t>(1));
  HS_EXPECT_NEAR(std::abs(dot(trail.get(0).get(), b.get())), 1.0f, 1e-4f);
}

/**
 * @brief Verifies clear() empties the trail.
 */
inline void test_orientation_trail_clear() {
  hs::random().seed(1337);
  Animation::OrientationTrail<Orientation<4>, 8> trail;
  Orientation<4> a;
  trail.record(a);
  trail.record(a);
  trail.clear();
  HS_EXPECT_EQ(trail.length(), static_cast<size_t>(0));
}

// ============================================================================
// Borrow-contract guards (compile-time)
// ----------------------------------------------------------------------------
// Motion/MeshMorph/Lerp/MobiusFlow store non-owning borrows of effect-owned
// state. These static_asserts lock the contract: an effect-owned lvalue is
// accepted, a temporary (which would dangle) is rejected.
// ============================================================================
namespace borrow_guard {
/** @brief Orientation alias used by the borrow-contract static_asserts. */
using Ori = Orientation<16>;
/** @brief Mesh draw-callable alias used by the MeshMorph borrow asserts. */
using DrawFn = Fn<void(Canvas &, const MeshState &, float), 8>;

static_assert(std::is_constructible_v<Animation::Motion<288, 16>, Ori &,
                                      ProceduralPath &, int>,
              "Motion must accept an lvalue (effect-owned) path");
static_assert(!std::is_constructible_v<Animation::Motion<288, 16>, Ori &,
                                       ProceduralPath &&, int>,
              "Motion must REJECT a temporary path (would dangle)");

static_assert(std::is_constructible_v<Animation::MeshMorph, const MeshState &,
                                      const MeshState &, Arena &, DrawFn &,
                                      DrawFn &, int>,
              "MeshMorph must accept lvalue (effect-owned) callables");
static_assert(!std::is_constructible_v<Animation::MeshMorph, const MeshState &,
                                       const MeshState &, Arena &, DrawFn &&,
                                       DrawFn &&, int>,
              "MeshMorph must REJECT temporary callables (would dangle)");
static_assert(!std::is_constructible_v<Animation::MeshMorph, const MeshState &,
                                       const MeshState &, Arena &, DrawFn &,
                                       DrawFn &&, int>,
              "MeshMorph must REJECT a temporary in either callable slot");

/**
 * @brief Minimal lerp subject used to probe Lerp's deleted rvalue overloads.
 */
struct Lerpable {
  float v = 0.0f; /**< Interpolated scalar payload. */
  /**
   * @brief Linearly interpolates this subject's value between two endpoints.
   * @param a Start subject (source value).
   * @param b Target subject (destination value).
   * @param t Interpolation parameter in [0, 1].
   */
  void lerp(const Lerpable &a, const Lerpable &b, float t) {
    v = a.v + (b.v - a.v) * t;
  }
};
static_assert(std::is_constructible_v<Animation::Lerp, Lerpable &,
                                      const Lerpable &, const Lerpable &, int,
                                      EasingFn>,
              "Lerp must accept lvalue (effect-owned) start/target");
static_assert(!std::is_constructible_v<Animation::Lerp, Lerpable &,
                                       const Lerpable &&, const Lerpable &, int,
                                       EasingFn>,
              "Lerp must REJECT a temporary start (would dangle)");
static_assert(!std::is_constructible_v<Animation::Lerp, Lerpable &,
                                       const Lerpable &, const Lerpable &&, int,
                                       EasingFn>,
              "Lerp must REJECT a temporary target (would dangle)");

static_assert(std::is_constructible_v<Animation::MobiusFlow, MobiusParams &,
                                      const float &, const float &, int>,
              "MobiusFlow must accept lvalue (effect-owned) scalars");
static_assert(!std::is_constructible_v<Animation::MobiusFlow, MobiusParams &,
                                       const float &&, const float &, int>,
              "MobiusFlow must REJECT a temporary num_rings (would dangle)");
static_assert(!std::is_constructible_v<Animation::MobiusFlow, MobiusParams &,
                                       const float &, const float &&, int>,
              "MobiusFlow must REJECT a temporary num_lines (would dangle)");
} // namespace borrow_guard

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Verifies Motion and Rotation size their orientation trail through the
 * shared Animation::rotation_substeps() helper, with a tight ceil.
 * @details Both animations route through the same helper so an identical sweep
 * yields an identical subdivision; each sub-interval must stay within MAX.
 */
inline void test_rotation_substeps_shared_and_tight() {
  constexpr float MAX = 0.1f;
  // Always at least 1, even for a sub-threshold angle.
  HS_EXPECT_EQ(Animation::rotation_substeps(0.0f, MAX), 1);
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX * 0.5f, MAX), 1);
  // Tight ceil: N*MAX needs exactly N subdivisions.
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX, MAX), 1);
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX * 3.0f, MAX), 3);
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX * 3.2f, MAX), 4);
  // Every sub-interval stays within MAX.
  for (float a = 0.0f; a < 2.0f; a += 0.013f) {
    int n = Animation::rotation_substeps(a, MAX);
    HS_EXPECT_GE(n, 1);
    HS_EXPECT_LE(a / n, MAX + 1e-6f);
  }
}

/**
 * @brief Verifies Rotation::step does not discard sub-TOLERANCE increments.
 * @details A rotation slow enough that each frame's delta is below TOLERANCE
 * (1e-4 rad) must still accumulate those deltas so the orientation actually
 * turns over many frames. ease_linear is linear here, so the per-frame delta is
 * total_angle / duration.
 */
inline void test_rotation_accumulates_subthreshold_deltas() {
  using Ori = Orientation<16>;
  Ori o; // identity
  // 0.05 rad over 1000 frames => 5e-5 rad/frame, half of TOLERANCE, so every
  // frame's raw delta is below the early-out threshold.
  Animation::Rotation<288, 16> rot(o, Z_AXIS, 0.05f, 1000, ease_linear);
  for (int i = 0; i < 20; ++i)
    rot.step(fake_canvas());
  // A Z rotation sends +X toward +Y.
  Vector v = o.orient(X_AXIS, o.length() - 1);
  HS_EXPECT_GT(v.y, 5e-4f);
  HS_EXPECT_NEAR(v.x, 1.0f, 1e-3f);
}

/**
 * @brief Verifies two animations sharing one Orientation COMPOSE their
 * sub-frame motion-blur history within a frame instead of clobbering it.
 * @details A per-animation collapse would discard the first animation's
 * freshly-built sub-frame trail. The decisive signature is the OLDEST sub-frame
 * (index 0): with composition it still reflects the pre-frame orientation
 * (identity here). This test touches the global Timeline; Rotation::step never
 * dereferences the canvas, and the global cursor state is reset around it.
 */
inline void test_timeline_shared_orientation_composes_motion_blur() {
  using Ori = Orientation<16>;
  global_timeline_num_events = 0;
  global_timeline_t = 0;

  Ori o; // identity, single frame
  Timeline tl;
  // Two quarter-turn rotations about the same axis, each completing in one frame.
  tl.add(0, Animation::Rotation<288, 16>(o, Z_AXIS, PI_F / 2, 1, ease_linear));
  tl.add(0, Animation::Rotation<288, 16>(o, Z_AXIS, PI_F / 2, 1, ease_linear));
  tl.step(fake_canvas());

  HS_EXPECT_GE(o.length(), static_cast<size_t>(2));
  // Oldest sub-frame is the pre-frame orientation (identity): +X stays +X.
  Vector oldest = o.orient(X_AXIS, 0);
  HS_EXPECT_NEAR(oldest.x, 1.0f, 1e-3f);
  HS_EXPECT_NEAR(oldest.y, 0.0f, 1e-3f);
  // Newest reflects both rotations (a half turn): +X -> -X.
  Vector newest = o.orient(X_AXIS, o.length() - 1);
  HS_EXPECT_NEAR(newest.x, -1.0f, 1e-3f);

  global_timeline_num_events = 0;
  global_timeline_t = 0;
}

/**
 * @brief Verifies Timeline schedules events by start frame and removes
 * completed one-shots so later events can run after earlier ones finish.
 * @details An event added with in_frames > 0 stays dormant until t reaches its
 * start, then steps. Uses the global Timeline; each Timeline is scoped so the
 * live-guard balances, and Transition::step never dereferences the canvas.
 */
inline void test_timeline_sequences_events_by_start_frame() {
  Timeline tl;
  float a = 0.0f, b = 0.0f;
  tl.add(0, Animation::Transition(a, 10.0f, 2, ease_linear)); // starts now
  tl.add(3, Animation::Transition(b, 20.0f, 2, ease_linear)); // delayed 3 frames

  tl.step(fake_canvas()); // t=1
  HS_EXPECT_GT(a, 0.0f);
  HS_EXPECT_NEAR(b, 0.0f, 1e-6f);

  tl.step(fake_canvas()); // t=2: a completes and is removed
  HS_EXPECT_NEAR(a, 10.0f, 1e-3f);
  HS_EXPECT_NEAR(b, 0.0f, 1e-6f); // still dormant (start=3)
  HS_EXPECT_EQ(global_timeline_num_events, 1);

  tl.step(fake_canvas()); // t=3: b starts
  HS_EXPECT_GT(b, 0.0f);
  HS_EXPECT_LT(b, 20.0f);

  tl.step(fake_canvas()); // t=4: b completes
  HS_EXPECT_NEAR(b, 20.0f, 1e-3f);
  HS_EXPECT_EQ(global_timeline_num_events, 0);
}

/**
 * @brief Verifies the latest representable start frame remains schedulable.
 */
inline void test_timeline_accepts_maximum_start_frame() {
  Timeline tl;
  global_timeline_t = UINT32_MAX - 2;
  float value = 0.0f;
  tl.add(2, Animation::Transition(value, 1.0f, 1, ease_linear));

  HS_EXPECT_EQ(global_timeline_events[0].start, UINT32_MAX);
  tl.step(fake_canvas());
  HS_EXPECT_NEAR(value, 0.0f, 1e-6f);
  tl.step(fake_canvas());
  HS_EXPECT_NEAR(value, 1.0f, 1e-6f);
  HS_EXPECT_EQ(global_timeline_t, UINT32_MAX);
}

/**
 * @brief Verifies a repeating animation rewinds at the end of each cycle
 * instead of being removed, replaying the curve.
 * @details Mutation writes f(easing(t/duration)) each step, so after completion
 * and rewind the next step drops back to the mid-cycle value; a non-rewinding
 * timer would clamp at 1.
 */
inline void test_timeline_repeating_animation_rewinds_each_cycle() {
  Timeline tl;
  float v = -1.0f;
  tl.add(0, Animation::Mutation(
                v, [](float e) { return e; }, 2, ease_linear, /*repeat=*/true));

  tl.step(fake_canvas()); // t=1 -> v = eased(0.5) = 0.5
  HS_EXPECT_NEAR(v, 0.5f, 1e-3f);
  tl.step(fake_canvas()); // t=2 -> v = eased(1.0) = 1.0, done -> rewind
  HS_EXPECT_NEAR(v, 1.0f, 1e-3f);
  tl.step(fake_canvas()); // rewound to t=1 -> v = 0.5
  HS_EXPECT_NEAR(v, 0.5f, 1e-3f);

  HS_EXPECT_EQ(global_timeline_num_events, 1);
}

/**
 * @brief Verifies a canceled repeating animation is removed, not rewound and
 * replayed forever.
 * @details cancel() makes done() permanently true; repeats() must drop on
 * cancel so Timeline routes it through the removal branch instead of keeping it
 * as a per-frame, callback-firing zombie.
 */
inline void test_timeline_cancel_removes_repeating_animation() {
  Timeline tl;
  float v = -1.0f;
  auto *h = tl.add_get(0, Animation::Mutation(
                              v, [](float e) { return e; }, 2, ease_linear,
                              /*repeat=*/true));
  tl.step(fake_canvas());
  tl.step(fake_canvas()); // completes a cycle, rewinds, and stays (repeating)
  HS_EXPECT_EQ(global_timeline_num_events, 1);

  h->cancel();
  tl.step(fake_canvas()); // canceled: done() && !repeats() -> removed
  HS_EXPECT_EQ(global_timeline_num_events, 0);

  // A removed event stops stepping: v stays frozen instead of oscillating.
  const float v_frozen = v;
  for (int i = 0; i < 6; ++i)
    tl.step(fake_canvas());
  HS_EXPECT_NEAR(v, v_frozen, 1e-6f);
}

/**
 * @brief Verifies step() compacts the event array when a non-repeating event is
 * removed, and relocated survivors keep stepping from their new positions.
 * @details Later survivors are relocated (move_into) into the freed slots. The
 * decisive check is that the originally-LAST event (relocated furthest) still
 * reaches its own target.
 */
inline void test_timeline_compaction_preserves_later_events() {
  Timeline tl;
  float a = 0.0f, b = 0.0f, c = 0.0f;
  tl.add(0, Animation::Transition(a, 10.0f, 1, ease_linear));  // completes at t=1
  tl.add(0, Animation::Transition(b, 100.0f, 5, ease_linear)); // in-flight survivor
  tl.add(0, Animation::Transition(c, 200.0f, 5, ease_linear)); // in-flight survivor
  HS_EXPECT_EQ(global_timeline_num_events, 3);

  tl.step(fake_canvas()); // t=1: a done+removed; b,c step once and shift down
  HS_EXPECT_NEAR(a, 10.0f, 1e-3f);
  HS_EXPECT_GT(b, 0.0f);
  HS_EXPECT_GT(c, 0.0f);
  HS_EXPECT_EQ(global_timeline_num_events, 2);
  float b_after1 = b, c_after1 = c;

  for (int i = 0; i < 4; ++i)
    tl.step(fake_canvas()); // t=2..5: finish the relocated survivors
  HS_EXPECT_GT(b, b_after1);
  HS_EXPECT_GT(c, c_after1);
  HS_EXPECT_NEAR(b, 100.0f, 1e-2f);
  HS_EXPECT_NEAR(c, 200.0f, 1e-2f);
  HS_EXPECT_EQ(global_timeline_num_events, 0);
}

/**
 * @brief Verifies .then() fires on completion and a callback may schedule a
 * follow-up on the same Timeline mid-step.
 * @details step() appends such events past the active snapshot, then gap-fills
 * them into the freed slots (the add-during-callback path). The follow-up is
 * added this frame but only runs on the next.
 */
inline void test_timeline_then_chains_follow_up_event() {
  Timeline tl;
  float a = 0.0f, b = 0.0f;
  // 'a' completes in one frame; its .then() schedules 'b' to start immediately.
  tl.add(0, Animation::Transition(a, 10.0f, 1, ease_linear).then([&]() {
    tl.add(0, Animation::Transition(b, 20.0f, 1, ease_linear));
  }));

  tl.step(fake_canvas()); // t=1: a completes -> callback adds b (gap-filled in)
  HS_EXPECT_NEAR(a, 10.0f, 1e-3f);
  HS_EXPECT_NEAR(b, 0.0f, 1e-6f);              // b added this frame, not yet run
  HS_EXPECT_EQ(global_timeline_num_events, 1);

  tl.step(fake_canvas()); // t=2: the chained event runs and completes
  HS_EXPECT_NEAR(b, 20.0f, 1e-3f);
  HS_EXPECT_EQ(global_timeline_num_events, 0);
}

/**
 * @brief Verifies a repeating timer routes its hook through post_callback() so
 * an attached .then() fires on every trigger.
 * @details A repeating timer (duration=-1, self-resetting) never reaches
 * done(), so the Timeline never fires its per-cycle .then(); the timer must
 * fire the hook itself to match the then() contract.
 */
inline void test_repeating_timer_fires_then_each_cycle() {
  Timeline tl;
  int triggers = 0, thens = 0;
  tl.add(0, Animation::PeriodicTimer(
                3, [&](Canvas &) { triggers++; }, /*repeat=*/true)
                .then([&]() { thens++; }));
  for (int i = 0; i < 9; ++i)
    tl.step(fake_canvas()); // triggers at t=3,6,9
  HS_EXPECT_EQ(triggers, 3);
  HS_EXPECT_EQ(thens, 3);
}

/**
 * @brief Verifies clear() destroys all events and resets the frame cursor,
 * leaving the timeline reusable from t=0.
 * @details This is the in-place reset the singleton offers in lieu of
 * reassignment.
 */
inline void test_timeline_clear_resets_state() {
  Timeline tl;
  float a = 0.0f;
  tl.add(0, Animation::Transition(a, 10.0f, 5, ease_linear));
  tl.step(fake_canvas());
  HS_EXPECT_EQ(global_timeline_num_events, 1);
  HS_EXPECT_EQ(global_timeline_t, 1u);

  tl.clear();
  HS_EXPECT_EQ(global_timeline_num_events, 0);
  HS_EXPECT_EQ(global_timeline_t, 0u); // cursor rewound

  // Reusable after clear().
  float b = 0.0f;
  tl.add(0, Animation::Transition(b, 5.0f, 1, ease_linear));
  tl.step(fake_canvas());
  HS_EXPECT_NEAR(b, 5.0f, 1e-3f);
  HS_EXPECT_EQ(global_timeline_num_events, 0);
}

/**
 * @brief Verifies add() is bounded by MAX_EVENTS (the shared global array
 * size): once full, a further add is rejected rather than overrunning.
 */
inline void test_timeline_full_guard_rejects_overflow() {
  Timeline tl;
  float sink = 0.0f;
  // Fill to capacity (these events share one float).
  for (int i = 0; i < Timeline::MAX_EVENTS; ++i)
    tl.add(0, Animation::Transition(sink, 1.0f, 10, ease_linear));
  HS_EXPECT_EQ(global_timeline_num_events, Timeline::MAX_EVENTS);

  // The overflow event has its own target so a silent enqueue would show up.
  float rejected = 0.0f;
  tl.add(0, Animation::Transition(rejected, 42.0f, 10, ease_linear)); // past full
  HS_EXPECT_EQ(global_timeline_num_events,
               Timeline::MAX_EVENTS);

  tl.step(fake_canvas());
  HS_EXPECT_NEAR(rejected, 0.0f, 1e-6f); // never ran
}

/**
 * @brief Verifies Orientation::upsample SLERP-interpolates the recorded
 * sub-frames up to a target count (preserving endpoints) and collapse()
 * discards all but the newest.
 * @details These are the two primitives multi-animation motion blur is built
 * on.
 */
inline void test_orientation_upsample_then_collapse() {
  Orientation<8> o; // identity, 1 frame
  o.push(make_rotation(Z_AXIS, PI_F / 2)); // 2 frames: identity, +90 about Z
  HS_EXPECT_EQ(o.length(), static_cast<size_t>(2));

  o.upsample(5);
  HS_EXPECT_EQ(o.length(), static_cast<size_t>(5));

  // Endpoints preserved: frame 0 ~ identity (+X stays +X); frame 4 ~ the +90
  // rotation about Z (+X -> +Y).
  Vector f0 = o.orient(X_AXIS, 0);
  HS_EXPECT_NEAR(f0.x, 1.0f, 1e-3f);
  Vector f4 = o.orient(X_AXIS, 4);
  HS_EXPECT_NEAR(f4.y, 1.0f, 1e-3f);

  // SLERP is monotone: +X decreases across the interpolated frames.
  float prevx = 2.0f;
  for (int i = 0; i < 5; ++i) {
    float x = o.orient(X_AXIS, i).x;
    HS_EXPECT_LE(x, prevx + 1e-4f);
    prevx = x;
  }

  o.collapse();
  HS_EXPECT_EQ(o.length(), static_cast<size_t>(1));
  Vector c = o.orient(X_AXIS, 0);
  HS_EXPECT_NEAR(c.y, 1.0f, 1e-3f);
}

/**
 * @brief Verifies a repeating Motion does not drift across many cycles.
 * @details A repeating Motion advances its Orientation by relative deltas taken
 * between consecutive path frames, where each frame is a pure function of the
 * path parameter (point + tangent). Because the frame depends only on the phase,
 * the per-cycle product of deltas telescopes — there is no accumulating
 * quaternion chain to warp the traced curve. The decisive, precession-immune
 * signature is the set of rotation-INVARIANT internal angles between heads
 * sampled at fixed phases within a cycle: a rigid drift (holonomy) leaves them
 * unchanged, so any growth is genuine warp. A late cycle is compared against the
 * ideal Lissajous internal angles; they must stay at the first cycle's tiny
 * discretization residual.
 */
inline void test_motion_repeating_does_not_drift() {
  using Ori = Orientation<16>;
  global_timeline_num_events = 0;
  global_timeline_t = 0;

  constexpr int duration = 40;
  ProceduralPath path;
  // lissajous(.,.,.,0) == +Y, so the identity-start orientation places the head
  // on the path at phase 0.
  path.f = [](float t) { return lissajous(1.06f, 1.06f, 0.0f, t * 5.909f); };

  Ori o; // identity, single frame; orient(+Y) starts on the path
  const Vector node_v = Y_AXIS;

  Timeline tl;
  tl.add(0, Animation::Motion<288, 16>(o, path, duration, /*repeat=*/true));

  Vector late_heads[duration + 1]; // indexed by phase 1..duration
  const int cycles = 600;
  for (int c = 0; c < cycles; ++c) {
    for (int fr = 1; fr <= duration; ++fr) {
      tl.step(fake_canvas());
      if (c == cycles - 1) late_heads[fr] = o.orient(node_v);
    }
  }

  // Each phase's internal angle to the phase-1 anchor must match the ideal
  // Lissajous internal angle; a rigid precession leaves these untouched, so any
  // growth is genuine warp. Interior phases only (the boundary frame rewinds).
  const int anchor = 1;
  const Vector ideal_anchor = path.f((float)anchor / duration);
  for (int fr = 2; fr < duration; ++fr) {
    const Vector ideal_fr = path.f((float)fr / duration);
    HS_EXPECT_NEAR(angle_between(late_heads[anchor], late_heads[fr]),
                   angle_between(ideal_anchor, ideal_fr), 0.1f);
  }

  global_timeline_num_events = 0;
  global_timeline_t = 0;
}

/**
 * @brief A co-driver sharing a repeating Motion's Orientation survives the
 * repeat seam.
 * @details The old drift fix snapped the entire Orientation back to Motion's
 * captured anchor at every cycle boundary, which clobbered any other animation
 * driving the same Orientation. Motion now re-seats via a *relative* delta, so a
 * co-driver's accumulated rotation rides across the seam instead of being
 * discarded. With a CLOSED path Motion's per-cycle contribution telescopes to
 * identity, so the only thing that should move the shared orientation at a seam
 * is the co-driver's own small step — never a large snap-back. Assert the
 * probe's per-frame angular step stays bounded across many seams while its
 * cumulative travel is large (so the co-driver is provably active, not a no-op).
 */
inline void test_motion_codriven_survives_repeat_seam() {
  using Ori = Orientation<16>;
  global_timeline_num_events = 0;
  global_timeline_t = 0;

  const int duration = 30;
  ProceduralPath path;
  // A closed great circle: path(0) == path(1) with matching tangent, so Motion's
  // per-cycle delta product is identity.
  path.f = [](float t) {
    float a = 2.0f * PI_F * t;
    return Vector(std::cos(a), std::sin(a), 0.0f);
  };

  Ori o; // identity
  Timeline tl;
  // Repeating Motion + a repeating co-driver rotation about Y, both driving `o`.
  tl.add(0, Animation::Motion<288, 16>(o, path, duration, /*repeat=*/true));
  tl.add(0, Animation::Rotation<288, 16>(o, Y_AXIS, 2.0f * PI_F, duration,
                                         ease_linear, /*repeat=*/true));

  const Vector probe = Z_AXIS;
  Vector prev = o.orient(probe);
  float max_step = 0.0f;
  float total_travel = 0.0f;
  const int cycles = 8;
  for (int c = 0; c < cycles; ++c) {
    for (int fr = 1; fr <= duration; ++fr) {
      tl.step(fake_canvas());
      Vector cur = o.orient(probe);
      float step = angle_between(prev, cur);
      max_step = std::max(max_step, step);
      total_travel += step;
      prev = cur;
    }
  }

  // No single frame (seams included) snaps the orientation; 0.8 clears the
  // largest legitimate one-frame step but is far below a multi-radian snap.
  HS_EXPECT_LT(max_step, 0.8f);
  HS_EXPECT_GT(total_travel, 5.0f);

  global_timeline_num_events = 0;
  global_timeline_t = 0;
}

// ============================================================================
// ParticleSystem
// ----------------------------------------------------------------------------
// Covers the pool-state transitions: spawn (+ capacity guard), life-expiry with
// trail drain, and attractor kill-radius removal.
// ============================================================================

/**
 * @brief Verifies ParticleSystem::spawn adds particles and silently drops
 * spawns once the fixed pool is at capacity.
 */
inline void test_particle_system_spawn_and_capacity_guard() {
  hs::random().seed(1337);
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 4> ps; // CAPACITY = 4
  ps.init(arena);
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 0);

  ps.spawn(Vector(1, 0, 0), Vector(0, 0, 0), 0);
  ps.spawn(Vector(0, 1, 0), Vector(0, 0, 0), 1);
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 2);

  ps.spawn(Vector(0, 0, 1), Vector(0, 0, 0), 2);
  ps.spawn(Vector(-1, 0, 0), Vector(0, 0, 0), 3);
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 4);
  ps.spawn(Vector(0, -1, 0), Vector(0, 0, 0), 4); // capacity is 4 — rejected
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 4);
}

/**
 * @brief Verifies both valid ParticleSystem lifetime boundaries are preserved.
 */
inline void test_particle_system_lifetime_boundaries() {
  static uint8_t buf[4096];
  {
    Arena arena(buf, sizeof(buf));
    Animation::ParticleSystem<32, 1> ps;
    ps.init(arena, 0.85f, 0.0f, 1.0f);
    HS_EXPECT_EQ(static_cast<int>(ps.max_life), 1);
  }
  {
    Arena arena(buf, sizeof(buf));
    Animation::ParticleSystem<32, 1> ps;
    ps.init(arena, 0.85f, 0.0f, 65535.0f);
    HS_EXPECT_EQ(static_cast<int>(ps.max_life), 65535);
  }
}

/**
 * @brief Verifies a particle is removed only after its life reaches 0 and its
 * recorded trail drains to empty.
 */
inline void test_particle_system_expires_after_life_and_trail_drain() {
  hs::random().seed(1337);
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 4> ps;
  // gravity 0 + zero velocity => the particle never moves; only life reaching 0
  // and its trail draining to empty ends it.
  ps.init(arena, /*friction=*/0.85f, /*gravity=*/0.0f, /*max_life=*/3.0f);
  ps.spawn(Vector(1, 0, 0), Vector(0, 0, 0), 0);

  ps.step(fake_canvas());
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 1);

  // life=3 records 2 trail frames (life stays >0 after the decrement on frames
  // 1,2), then drains one per frame: reclaimed exactly on frame 2*(life-1)=4.
  const int reclaim_frame = 2 * (3 - 1);
  for (int frame = 2; frame < reclaim_frame; ++frame) {
    ps.step(fake_canvas());
    HS_EXPECT_EQ(static_cast<int>(ps.active()), 1);
  }
  ps.step(fake_canvas()); // frame reclaim_frame
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 0);
}

/**
 * @brief Verifies an attractor removes a particle inside its kill_radius via
 * the kill check rather than by life expiry.
 */
inline void test_particle_system_attractor_kills_within_radius() {
  hs::random().seed(1337);
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 4> ps;
  ps.init(arena, /*friction=*/0.85f, /*gravity=*/0.001f, /*max_life=*/600.0f);
  ps.spawn(Vector(1, 0, 0), Vector(0, 0, 0), 0);
  // Attractor co-located (distance 0 < kill_radius): removed by the kill check,
  // not life expiry (life is 600).
  ps.add_attractor(Vector(1, 0, 0), /*strength=*/1.0f, /*kill_radius=*/0.5f,
                   /*event_horizon=*/2.0f);
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 1);

  ps.step(fake_canvas());
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 0);
}

/**
 * @brief Verifies spawn() initializes a particle's fields and that one step
 * advances it (life decrements, trail records) without dropping it.
 * @details The capacity test only inspects active_count, so a spawn that left
 *          the position/velocity/seed/life unset — or a step that mis-managed a
 *          freshly spawned particle — would pass it. Here the spawned particle
 *          is checked field-by-field (it carries the requested position,
 *          velocity, and seed, the system's max_life, and an empty trail), then
 *          a single step is asserted to keep it alive, decrement its life by
 *          one, and record its first trail point.
 */
inline void test_particle_system_spawn_initializes_and_steps() {
  hs::random().seed(1337);
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 4> ps;
  // gravity 0 and no attractor: the only state change is the step's bookkeeping.
  ps.init(arena, /*friction=*/0.85f, /*gravity=*/0.0f, /*max_life=*/120.0f);

  const Vector pos(0.6f, 0.0f, 0.8f);
  const Vector vel(0.0f, 0.01f, 0.0f);
  ps.spawn(pos, vel, /*seed=*/42);
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 1);

  const auto &p = ps.pool[0];
  HS_EXPECT_NEAR(p.position.x, pos.x, 0.0f);
  HS_EXPECT_NEAR(p.position.y, pos.y, 0.0f);
  HS_EXPECT_NEAR(p.position.z, pos.z, 0.0f);
  HS_EXPECT_NEAR(p.velocity.x, vel.x, 0.0f);
  HS_EXPECT_NEAR(p.velocity.y, vel.y, 0.0f);
  HS_EXPECT_NEAR(p.velocity.z, vel.z, 0.0f);
  HS_EXPECT_EQ(static_cast<int>(p.color_seed), 42);
  HS_EXPECT_EQ(static_cast<int>(p.life), 120);       // == system max_life
  HS_EXPECT_EQ(static_cast<int>(p.history_length()), 0);

  ps.step(fake_canvas());
  HS_EXPECT_EQ(static_cast<int>(ps.active()), 1);
  HS_EXPECT_EQ(static_cast<int>(ps.pool[0].life), 119);        // life-- per step
  HS_EXPECT_GT(static_cast<int>(ps.pool[0].history_length()), 0);
}

/**
 * @brief Pins the attractor kill check at its radius boundary, not just at
 * distance 0.
 * @details The kill test places the attractor on top of the particle (distance
 *          0), so it can't tell `dist < kill_radius` from `dist <= kill_radius`.
 *          Here a particle at (1,0,0) faces an attractor placed an exact
 *          distance away along +x: just inside the radius it is killed, exactly
 *          at the radius it survives (the check is strict `<`), and just outside
 *          it survives. The kill check reads the start-of-step position, so one
 *          step suffices and the distance is the spawn distance.
 */
inline void test_particle_system_attractor_kill_radius_boundary() {
  hs::random().seed(1337);
  static uint8_t buf[256 * 1024];
  constexpr float kr = 0.5f;

  auto survives = [&](float dist) {
    Arena arena(buf, sizeof(buf));
    Animation::ParticleSystem<32, 4> ps;
    ps.init(arena, /*friction=*/0.85f, /*gravity=*/0.001f, /*max_life=*/600.0f);
    ps.spawn(Vector(1, 0, 0), Vector(0, 0, 0), 0);
    ps.add_attractor(Vector(1.0f + dist, 0, 0), /*strength=*/1.0f,
                     /*kill_radius=*/kr, /*event_horizon=*/0.0f);
    ps.step(fake_canvas());
    return ps.active() == 1;
  };

  HS_EXPECT_FALSE(survives(kr - 0.01f)); // inside the radius → killed
  HS_EXPECT_TRUE(survives(kr));          // exactly at the radius → survives (strict <)
  HS_EXPECT_TRUE(survives(kr + 0.01f));
}

// ============================================================================
// Sprite (opacity envelope + paused-hold)
// ----------------------------------------------------------------------------
// Sprite::step computes an opacity and forwards it to the user draw_fn; the
// captured opacity drives the envelope and paused-hold assertions.
// ============================================================================

/**
 * @brief Verifies Sprite drives a fade-in -> full-opacity plateau -> fade-out
 * envelope across its duration, completing on the final frame.
 */
inline void test_sprite_fade_in_plateau_fade_out_envelope() {
  std::vector<float> ops;
  const int dur = 10, fade_in = 3, fade_out = 3;
  Animation::Sprite s(
      [&](Canvas &, float o) { ops.push_back(o); }, dur, fade_in, ease_linear,
      fade_out, ease_linear);
  for (int i = 0; i < dur; ++i)
    s.step(fake_canvas()); // observed at t = 1..10

  HS_EXPECT_EQ(ops.size(), static_cast<size_t>(dur));
  // Fade-in (linear): t=1 -> 1/3, then rising.
  HS_EXPECT_NEAR(ops[0], 1.0f / 3.0f, 1e-3f);
  HS_EXPECT_GT(ops[1], ops[0]);
  // Plateau: fully opaque between fade-in and fade-out.
  HS_EXPECT_NEAR(ops[3], 1.0f, 1e-3f);
  HS_EXPECT_NEAR(ops[5], 1.0f, 1e-3f);
  // Fade-out: transparent on the final frame.
  HS_EXPECT_LT(ops[8], ops[5]);
  HS_EXPECT_NEAR(ops[9], 0.0f, 1e-3f);
  HS_EXPECT_TRUE(s.done());
}

/**
 * @brief Verifies that when fade_in + fade_out exceed duration the fades scale
 * proportionally into a continuous triangle that still peaks at full opacity.
 * @details There must be no jump where the fade-in hands off to the fade-out;
 * the slider ranges make this configuration user-reachable.
 */
inline void test_sprite_overlapping_fades_stay_continuous() {
  std::vector<float> ops;
  const int dur = 10, fade_in = 8, fade_out = 8; // 8 + 8 > 10 => overlap, scaled to 5 + 5
  Animation::Sprite s(
      [&](Canvas &, float o) { ops.push_back(o); }, dur, fade_in, ease_linear,
      fade_out, ease_linear);
  for (int i = 0; i < dur; ++i)
    s.step(fake_canvas()); // observed at t = 1..10

  HS_EXPECT_EQ(ops.size(), static_cast<size_t>(dur));
  // Scaled to 5 + 5: slope 1/5 per frame, so no single-frame step exceeds ~0.2.
  float max_jump = 0.0f;
  for (size_t i = 1; i < ops.size(); ++i)
    max_jump = std::max(max_jump, std::abs(ops[i] - ops[i - 1]));
  HS_EXPECT_LT(max_jump, 0.3f);
  // Rises to full opacity at the apex, then falls (a triangle).
  HS_EXPECT_NEAR(ops[4], 1.0f, 1e-3f);
  HS_EXPECT_GT(ops[4], ops[0]);
  HS_EXPECT_GT(ops[4], ops[9]);
}

/**
 * @brief Verifies that while its paused flag is set a Sprite holds its current
 * frame: the timer never advances and it never expires, yet it keeps drawing at
 * its held opacity.
 */
inline void test_sprite_paused_holds_frame() {
  bool paused = true;
  int draws = 0;
  float last_op = -1.0f;
  // No fade in/out => held opacity is full.
  Animation::Sprite s(
      [&](Canvas &, float o) { draws++; last_op = o; }, /*duration=*/3,
      /*fade_in=*/0, ease_linear, /*fade_out=*/0, ease_linear, &paused);

  for (int i = 0; i < 10; ++i)
    s.step(fake_canvas());
  HS_EXPECT_FALSE(s.done());          // timer never advanced while paused
  HS_EXPECT_EQ(draws, 10);            // but it kept drawing every frame
  HS_EXPECT_NEAR(last_op, 1.0f, 1e-6f);

  paused = false;
  for (int i = 0; i < 3; ++i)
    s.step(fake_canvas());
  HS_EXPECT_TRUE(s.done());
}

// ============================================================================
// MeshCarousel segues
// ============================================================================

/**
 * @brief Verifies Segue::Crossfade schedules one fading sprite through the
 * carousel seam and returns a next-transition delay that overlaps consecutive
 * sprites by exactly the fade window.
 */
inline void test_crossfade_segue_schedules_overlapping_sprite() {
  Timeline tl;
  std::vector<float> ops;
  const int dur = 10, window = 3;
  MeshCarousel<Segue::Crossfade> carousel;
  int next_delay = carousel.schedule_segue(
      tl, [&](Canvas &, float o) { ops.push_back(o); }, dur, window);
  HS_EXPECT_EQ(next_delay, dur - window);
  HS_EXPECT_EQ(global_timeline_num_events, 1);

  for (int i = 0; i < dur; ++i)
    tl.step(fake_canvas()); // observed at t = 1..10

  HS_EXPECT_EQ(ops.size(), static_cast<size_t>(dur));
  // Fade-in ramp, full-opacity plateau, transparent on the final frame.
  HS_EXPECT_LT(ops[0], 1.0f);
  HS_EXPECT_NEAR(ops[4], 1.0f, 1e-3f);
  HS_EXPECT_NEAR(ops[dur - 1], 0.0f, 1e-3f);
  HS_EXPECT_EQ(global_timeline_num_events, 0);
}

/**
 * @brief Verifies Segue::Crossfade clamps the fade window to half the duration
 * so fade windows never overlap and sprites cannot pile up beyond two.
 */
inline void test_crossfade_segue_clamps_fade_to_half_duration() {
  Timeline tl;
  const int dur = 10, window = 9; // > dur/2, clamps to 5
  MeshCarousel<Segue::Crossfade> carousel;
  int next_delay =
      carousel.schedule_segue(tl, [](Canvas &, float) {}, dur, window);
  HS_EXPECT_EQ(next_delay, dur - dur / 2);
}

/**
 * @brief Verifies the default (Base) scheduling is sequential: the returned
 * delay equals the full duration, so consecutive sprites never coexist and a
 * single mesh renders per frame.
 */
inline void test_sequential_segue_never_overlaps_sprites() {
  Timeline tl;
  std::vector<float> ops;
  const int dur = 10, window = 3;
  MeshCarousel<Segue::SpinFlip> carousel;
  int next_delay = carousel.schedule_segue(
      tl, [&](Canvas &, float p) { ops.push_back(p); }, dur, window);
  HS_EXPECT_EQ(next_delay, dur);
  HS_EXPECT_EQ(global_timeline_num_events, 1);

  for (int i = 0; i < dur; ++i)
    tl.step(fake_canvas()); // observed at t = 1..10

  HS_EXPECT_EQ(ops.size(), static_cast<size_t>(dur));
  // Phase ramps up, plateaus at 1, and returns to 0 on the final frame.
  HS_EXPECT_LT(ops[0], 1.0f);
  HS_EXPECT_NEAR(ops[4], 1.0f, 1e-3f);
  HS_EXPECT_NEAR(ops[dur - 1], 0.0f, 1e-3f);
}

/**
 * @brief Verifies Segue::Dissolve's complementary masks partition the canvas:
 *        every pixel is owned by exactly one of the two meshes, and the
 *        incoming share tracks the phase.
 * @details The partition is what caps a two-mesh transition at one mesh's scan
 *          cost, so a mask pair that double-covered or dropped pixels would
 *          both corrupt the image and double the frame.
 */
inline void test_dissolve_segue_masks_partition_the_canvas() {
  constexpr int W = 64, H = 32;
  Segue::Dissolve dissolve;
  for (float phase : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f}) {
    PixelMask in = dissolve.mask(phase, 7u, true);
    PixelMask out = dissolve.mask(phase, 7u, false);
    int owned_in = 0;
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        bool a = in.owns(x, y), b = out.owns(x, y);
        HS_EXPECT_TRUE(a != b);
        owned_in += a ? 1 : 0;
      }
    // Hash spread is not perfectly uniform, so allow a few percent of slack.
    float share = static_cast<float>(owned_in) / (W * H);
    HS_EXPECT_NEAR(share, phase, 0.05f);
  }
}

/**
 * @brief Verifies Segue::Dissolve re-rolls its pattern every frame and
 *        re-seeds per transition, so successive frames dither rather than
 *        freezing one stochastic partition on screen.
 */
inline void test_dissolve_segue_reseeds_per_frame_and_transition() {
  Segue::Dissolve dissolve;
  PixelMask f0 = dissolve.mask(0.5f, 0u, true);
  PixelMask f1 = dissolve.mask(0.5f, 1u, true);
  HS_EXPECT_NE(f0.salt, f1.salt);
  HS_EXPECT_EQ(f0.threshold, f1.threshold);

  uint32_t before = dissolve.mask(0.5f, 0u, true).salt;
  dissolve.retarget(Vector(0, 1, 0));
  HS_EXPECT_NE(dissolve.mask(0.5f, 0u, true).salt, before);
}

/**
 * @brief Verifies the Base shading hooks are identities: full opacity and
 * coverage, unmodified edge distance and color, visible except at ~zero phase.
 */
inline void test_segue_base_hooks_are_identity() {
  Segue::Base base;
  HS_EXPECT_NEAR(base.opacity(0.3f), 1.0f, 1e-6f);
  float t = 0.42f;
  HS_EXPECT_NEAR(base.fill(t, 0.3f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(t, 0.42f, 1e-6f);
  Color4 c(Pixel(1000, 2000, 3000), 0.5f);
  Color4 g = base.grade(c, 0.3f);
  HS_EXPECT_EQ(g.color.r, c.color.r);
  HS_EXPECT_EQ(g.color.g, c.color.g);
  HS_EXPECT_EQ(g.color.b, c.color.b);
  HS_EXPECT_TRUE(base.visible(0.5f));
  HS_EXPECT_FALSE(base.visible(0.0f));
}

/**
 * @brief Verifies IrisBloom's fill contracts faces toward their centers: at
 * full phase everything survives, at mid phase only fragments deeper than the
 * inset do, and the surviving core renormalizes to the full gradient.
 */
inline void test_iris_bloom_fill_contracts_to_face_centers() {
  Segue::IrisBloom iris;
  // Full phase: everything covered, t unchanged.
  float t = 0.3f;
  HS_EXPECT_NEAR(iris.fill(t, 1.0f), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(t, 0.3f, 1e-3f);
  // Mid phase: shallow fragments culled...
  t = 0.3f;
  HS_EXPECT_NEAR(iris.fill(t, 0.5f), 0.0f, 1e-6f);
  // ...deep fragments survive with t renormalized over the shrunken core.
  t = 0.9f;
  HS_EXPECT_NEAR(iris.fill(t, 0.5f), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(t, 0.8f, 1e-3f);
}

/**
 * @brief Verifies Lace's fill is the inverse mask: only fragments within the
 * phase-driven band of an edge survive, renormalized across the band.
 */
inline void test_lace_fill_keeps_edge_band() {
  Segue::Lace lace;
  // Full phase: everything covered.
  float t = 0.9f;
  HS_EXPECT_NEAR(lace.fill(t, 1.0f), 1.0f, 1e-3f);
  // Narrow band: deep fragments culled, near-edge fragments survive.
  t = 0.5f;
  HS_EXPECT_NEAR(lace.fill(t, 0.3f), 0.0f, 1e-6f);
  t = 0.15f;
  HS_EXPECT_NEAR(lace.fill(t, 0.3f), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(t, 0.5f, 1e-3f);
}

/**
 * @brief Verifies the shared sweep front: full everywhere at phase 1, dark
 * everywhere at phase 0, monotone in phase, and higher offsets extinguish
 * earlier.
 */
inline void test_sweep_phase_front_ordering() {
  const float band = 0.25f;
  for (float o : {0.0f, 0.5f, 1.0f}) {
    HS_EXPECT_NEAR(Segue::sweep_phase(1.0f, o, band), 1.0f, 1e-3f);
    HS_EXPECT_NEAR(Segue::sweep_phase(0.0f, o, band), 0.0f, 1e-3f);
  }
  // Monotone in phase at fixed offset, sampled inside the front's ramp
  // (phase 0.14..0.39 for this offset/band under the sqrt ease).
  HS_EXPECT_GT(Segue::sweep_phase(0.3f, 0.5f, band),
               Segue::sweep_phase(0.2f, 0.5f, band));
  // At a fixed phase, a higher offset is further extinguished.
  HS_EXPECT_GT(Segue::sweep_phase(0.5f, 0.2f, band),
               Segue::sweep_phase(0.5f, 0.8f, band));
}

/**
 * @brief Verifies TerminatorSweep orders faces along its axis: the axis pole
 * extinguishes first (offset 1), the antipode last (offset 0), and the
 * per-face fade alpha is the perceptual square of the face-local phase with
 * exact endpoints.
 */
inline void test_terminator_sweep_orders_by_axis() {
  Segue::TerminatorSweep term;
  Vector axis = Vector(1.0f, 2.0f, -0.5f).normalized();
  term.retarget(axis);
  HS_EXPECT_NEAR(term.face_offset(axis, 0, 0), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(term.face_offset(-axis, 0, 0), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(term.face_offset(cross(axis, X_AXIS).normalized(), 0, 0), 0.5f,
                 1e-2f);
  HS_EXPECT_NEAR(term.opacity(0.4f), 0.16f, 1e-6f);
  HS_EXPECT_NEAR(term.opacity(0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(term.opacity(1.0f), 1.0f, 1e-6f);
}

/**
 * @brief Verifies TerminatorSweep's per-face fade is time-based: schedule()
 * derives the fade fractions so a face ramps over its fade length once the
 * front reaches it, with exact 0/1 window endpoints, and a window shorter than
 * the fade length degrades to one whole-sphere fade. Pins the range to a single
 * length so the per-face random collapses to one deterministic fraction.
 */
inline void test_terminator_sweep_fades_faces_over_fixed_frames() {
  Timeline tl;
  const int dur = 200, window = 32;
  MeshCarousel<Segue::TerminatorSweep> carousel;
  carousel.segue().fade_frames_min = 8.0f;
  carousel.segue().fade_frames_max = 8.0f;
  int next_delay =
      carousel.schedule_segue(tl, [](Canvas &, float) {}, dur, window);
  HS_EXPECT_EQ(next_delay, dur); // sequential: one mesh per frame
  const Segue::TerminatorSweep &term = carousel.segue();
  const float f = 8.0f / window;
  HS_EXPECT_NEAR(term.fade_frac_min, f, 1e-6f);
  HS_EXPECT_NEAR(term.fade_frac_max, f, 1e-6f);
  for (float o : {0.0f, 0.5f, 1.0f}) {
    float touch = o * (1.0f - f); // phase at which the front reaches the face
    HS_EXPECT_NEAR(term.face_phase(touch, o, f), 0.0f, 1e-5f);
    HS_EXPECT_NEAR(term.face_phase(touch + 0.5f * f, o, f), 0.5f, 1e-5f);
    HS_EXPECT_NEAR(term.face_phase(touch + f, o, f), 1.0f, 1e-5f);
    HS_EXPECT_NEAR(term.face_phase(1.0f, o, f), 1.0f, 1e-5f);
    HS_EXPECT_NEAR(term.face_phase(0.0f, o, f), 0.0f, 1e-5f);
  }
  carousel.schedule_segue(tl, [](Canvas &, float) {}, 8, 2);
  HS_EXPECT_NEAR(carousel.segue().fade_frac_max, 1.0f, 1e-6f);
}

/**
 * @brief Verifies TerminatorSweep draws each face's fade length randomly from
 * [fade_frames_min, fade_frames_max]: the fractions stay in range, are stable
 * for a given index, differ across faces, and preserve exact 0/1 window
 * endpoints for every per-face fade length.
 */
inline void test_terminator_sweep_per_face_fade_random_in_range() {
  hs::random().seed(1337);
  Timeline tl;
  const int dur = 400, window = 64;
  MeshCarousel<Segue::TerminatorSweep> carousel;
  carousel.segue().retarget(Y_AXIS); // rolls the per-face fade seed
  carousel.segue().fade_frames_min = 4.0f;
  carousel.segue().fade_frames_max = 16.0f;
  carousel.schedule_segue(tl, [](Canvas &, float) {}, dur, window);
  const Segue::TerminatorSweep &term = carousel.segue();
  const float lo = 4.0f / window, hi = 16.0f / window;
  const float first = term.face_fade_frac(0);
  bool varied = false;
  for (int i = 0; i < 256; ++i) {
    float ff = term.face_fade_frac(i);
    HS_EXPECT_TRUE(ff >= lo - 1e-6f && ff <= hi + 1e-6f);
    HS_EXPECT_NEAR(term.face_fade_frac(i), ff, 0.0f); // stable per index
    HS_EXPECT_NEAR(term.face_phase(1.0f, 0.7f, ff), 1.0f, 1e-5f);
    HS_EXPECT_NEAR(term.face_phase(0.0f, 0.7f, ff), 0.0f, 1e-5f);
    if (std::fabs(ff - first) > 1e-4f)
      varied = true;
  }
  HS_EXPECT_TRUE(varied);
}

/**
 * @brief Verifies Shockwave orders faces by angular distance from its origin:
 * nearest faces extinguish first, the antipode last.
 */
inline void test_shockwave_orders_by_distance_from_origin() {
  Segue::Shockwave wave;
  Vector origin = Vector(0.3f, -1.0f, 0.7f).normalized();
  wave.retarget(origin);
  HS_EXPECT_NEAR(wave.face_offset(origin, 0, 0), 1.0f, 1e-2f);
  HS_EXPECT_NEAR(wave.face_offset(-origin, 0, 0), 0.0f, 1e-2f);
  // Equidistant ring sits mid-order.
  HS_EXPECT_NEAR(wave.face_offset(cross(origin, X_AXIS).normalized(), 0, 0), 0.5f,
                 2e-2f);
}

/**
 * @brief Verifies Breakdown fades classes sequentially: reorder() yields a
 * permutation of the class ranks, offsets follow the ranks, and each class's
 * fade window is an abutting 1/n slice of the phase range — fully faded
 * before the next class starts.
 */
inline void test_breakdown_fades_classes_sequentially() {
  hs::random().seed(7u);
  Segue::Breakdown bd;
  constexpr int n = 5;
  // Per-face classes (dense [0, n), out of order with repeats); reorder derives
  // num_classes = max + 1 = n from them rather than taking a declared count.
  const std::vector<uint8_t> face_classes{2, 0, 4, 1, 3, 4, 0};
  bd.reorder(face_classes);
  HS_EXPECT_EQ(bd.num_classes, n);
  bool seen[n] = {};
  for (int c = 0; c < n; ++c) {
    HS_EXPECT_LT(static_cast<int>(bd.rank[c]), n);
    seen[bd.rank[c]] = true;
  }
  for (int r = 0; r < n; ++r)
    HS_EXPECT_TRUE(seen[r]); // a permutation: every rank assigned once

  Vector any(0.0f, 1.0f, 0.0f);
  for (int c = 0; c < n; ++c) {
    float o = bd.face_offset(any, 0, c);
    int r = bd.rank[c];
    HS_EXPECT_NEAR(o, static_cast<float>(n - 1 - r) / (n - 1), 1e-6f);
    // Class rank r fades linearly over one band of [BLACK_DWELL, 1]: gone at
    // the floor, untouched at the ceiling, abutting its neighbors' windows,
    // and every class is fully black through the dwell before the swap.
    float band = (1.0f - Segue::Breakdown::BLACK_DWELL) / n;
    float floor_p = Segue::Breakdown::BLACK_DWELL + (n - 1 - r) * band;
    HS_EXPECT_NEAR(bd.face_phase(floor_p, o), 0.0f, 1e-5f);
    HS_EXPECT_NEAR(bd.face_phase(floor_p + band, o), 1.0f, 1e-5f);
    HS_EXPECT_NEAR(bd.face_phase(Segue::Breakdown::BLACK_DWELL, o), 0.0f,
                   1e-5f);
    HS_EXPECT_NEAR(bd.face_phase(0.0f, o), 0.0f, 1e-6f);
  }
}

/**
 * @brief Verifies SpinFlip's warp is rigid: pairwise angles are preserved at
 * every phase and the plateau is the identity.
 */
inline void test_spin_flip_warp_is_rigid() {
  Segue::SpinFlip spin;
  spin.retarget(Vector(0.5f, 0.5f, -0.7f).normalized());
  Vector a = Vector(1.0f, 0.2f, 0.1f).normalized();
  Vector b = Vector(-0.3f, 0.9f, 0.4f).normalized();
  Vector wa = spin.warp(a, 0.3f), wb = spin.warp(b, 0.3f);
  HS_EXPECT_NEAR(dot(wa, wb), dot(a, b), 1e-3f);
  HS_EXPECT_NEAR(wa.length(), 1.0f, 1e-3f);
  HS_EXPECT_NEAR((spin.warp(a, 1.0f) - a).length(), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(spin.opacity(0.0f), 1.0f, 1e-6f); // never fades: blur hides the swap
}

/**
 * @brief Verifies GoldConvergence grades toward its gold at the swap and is
 * the identity on the plateau, with the mild opacity dip.
 */
inline void test_gold_convergence_grades_to_gold() {
  Segue::GoldConvergence gc;
  Color4 c(Pixel(1000, 2000, 3000), 0.8f);
  Color4 plateau = gc.grade(c, 1.0f);
  HS_EXPECT_EQ(plateau.color.r, c.color.r);
  HS_EXPECT_EQ(plateau.color.g, c.color.g);
  HS_EXPECT_EQ(plateau.color.b, c.color.b);
  Color4 swap = gc.grade(c, 0.0f);
  HS_EXPECT_EQ(swap.color.r, gc.gold.r);
  HS_EXPECT_EQ(swap.color.g, gc.gold.g);
  HS_EXPECT_EQ(swap.color.b, gc.gold.b);
  HS_EXPECT_NEAR(gc.opacity(0.0f), 0.4f, 1e-6f);
  HS_EXPECT_NEAR(gc.opacity(1.0f), 1.0f, 1e-6f);
}

// ============================================================================
// deep_tween (global_t span)
// ----------------------------------------------------------------------------
// deep_tween walks an OrientationTrail's frames and sub-frames, emitting a
// global t in [0,1] across the trail. It skips sub-frame 0 of every frame after
// the first (the shared boundary), so the emitted count is M + (N-1)*(M-1) for
// N frames of M sub-frames.
// ============================================================================

/**
 * @brief Verifies deep_tween emits a global t spanning [0,1] across the whole
 * trail, with the expected sample count and non-decreasing t.
 * @details The emitted count is M + (N-1)*(M-1) for N frames of M sub-frames,
 * since sub-frame 0 of every frame after the first is a skipped shared
 * boundary.
 */
inline void test_deep_tween_global_t_spans_unit_interval() {
  using Ori = Orientation<8>;
  Animation::OrientationTrail<Ori, 8> trail;
  const int N = 3, M = 3;
  for (int k = 0; k < N; ++k) {
    Ori o;
    o.push(make_rotation(Z_AXIS, 0.3f * (k + 1)));
    o.upsample(M);
    trail.record(o);
  }

  std::vector<float> gts;
  deep_tween(trail,
             [&](const Quaternion &, float gt) { gts.push_back(gt); });

  HS_EXPECT_EQ(gts.size(), static_cast<size_t>(M + (N - 1) * (M - 1)));
  HS_EXPECT_NEAR(gts.front(), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(gts.back(), 1.0f, 1e-6f);
  for (size_t i = 1; i < gts.size(); ++i)
    HS_EXPECT_GE(gts[i], gts[i - 1] - 1e-6f);
}

/**
 * @brief Verifies deep_tween still reaches t=1.0 when the newest frame is a
 * collapsed (motionless) length-1 frame.
 * @details deep_tween normalizes against the newest contentful frame so global
 * t reaches 1.0 at the leading edge, without re-plotting the boundary the
 * collapsed frame shares with its predecessor.
 */
inline void test_deep_tween_collapsed_newest_frame_reaches_one() {
  using Ori = Orientation<8>;
  Animation::OrientationTrail<Ori, 8> trail;
  for (int k = 0; k < 2; ++k) {
    Ori o;
    o.push(make_rotation(Z_AXIS, 0.3f * (k + 1)));
    o.upsample(3);
    trail.record(o);
  }
  Ori still; // length 1 — no motion this frame
  trail.record(still);

  std::vector<float> gts;
  deep_tween(trail,
             [&](const Quaternion &, float gt) { gts.push_back(gt); });

  // The collapsed tail frame is dropped, so the count matches the two
  // contentful frames (M=3 sub-frames each): M + (contentful-1)*(M-1).
  const size_t M = 3, contentful = 2;
  HS_EXPECT_EQ(gts.size(), M + (contentful - 1) * (M - 1));
  HS_EXPECT_NEAR(gts.front(), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(gts.back(), 1.0f, 1e-6f);
  for (size_t i = 1; i < gts.size(); ++i)
    HS_EXPECT_GE(gts[i], gts[i - 1] - 1e-6f);
}

/**
 * @brief Verifies that when every frame is motionless the lone plotted
 * orientation (the trail head) reads t = 1.0.
 * @details This mirrors tween(Orientation) for a lone snapshot; t = 0.0 would
 * render it invisible under quintic_kernel(0).
 */
inline void test_deep_tween_all_collapsed_reaches_one() {
  using Ori = Orientation<8>;
  Animation::OrientationTrail<Ori, 8> trail;
  for (int k = 0; k < 3; ++k) {
    Ori still; // length 1 — no motion any frame
    trail.record(still);
  }

  std::vector<float> gts;
  deep_tween(trail,
             [&](const Quaternion &, float gt) { gts.push_back(gt); });

  HS_EXPECT_EQ(gts.size(), static_cast<size_t>(1));
  HS_EXPECT_NEAR(gts.back(), 1.0f, 1e-6f);
}

/**
 * @brief Verifies deep_tween_frames groups deep_tween's emission by frame.
 * @details Callback k receives frame k's contribution — M sub-positions for
 * the first frame, M-1 (shared boundary skipped) after — and the concatenated
 * (quaternion, t) stream equals deep_tween's exactly.
 */
inline void test_deep_tween_frames_groups_flat_emission() {
  using Ori = Orientation<8>;
  Animation::OrientationTrail<Ori, 8> trail;
  const int N = 3, M = 4;
  for (int k = 0; k < N; ++k) {
    Ori o;
    o.push(make_rotation(Z_AXIS, 0.3f * (k + 1)));
    o.upsample(M);
    trail.record(o);
  }

  std::vector<std::pair<Quaternion, float>> flat;
  deep_tween(trail, [&](const Quaternion &q, float t) {
    flat.emplace_back(q, t);
  });

  size_t idx = 0;
  int frames = 0;
  deep_tween_frames(
      trail, [&](const Quaternion *qs, const float *ts, int count) {
        HS_EXPECT_EQ(count, frames == 0 ? M : M - 1);
        for (int i = 0; i < count && idx < flat.size(); ++i, ++idx) {
          HS_EXPECT_EQ(qs[i].r, flat[idx].first.r);
          HS_EXPECT_EQ(qs[i].v.x, flat[idx].first.v.x);
          HS_EXPECT_EQ(qs[i].v.y, flat[idx].first.v.y);
          HS_EXPECT_EQ(qs[i].v.z, flat[idx].first.v.z);
          HS_EXPECT_EQ(ts[i], flat[idx].second);
        }
        ++frames;
      });
  HS_EXPECT_EQ(idx, flat.size());
  HS_EXPECT_EQ(frames, N);
}

/**
 * @brief Verifies a motionless interior frame leaves no hole in the age ramp.
 * @details A moving / motionless / moving trail: the length-1 interior frame
 * contributes no sample and is excluded from the span, so the remaining moving
 * frames stay evenly spaced across [0,1] instead of straddling a ~1/span gap.
 */
inline void test_deep_tween_interior_motionless_frame_no_gap() {
  using Ori = Orientation<8>;
  Animation::OrientationTrail<Ori, 8> trail;
  {
    Ori o;
    o.push(make_rotation(Z_AXIS, 0.3f));
    o.upsample(3);
    trail.record(o);
  }
  {
    Ori still; // length 1 — no motion this interior frame
    trail.record(still);
  }
  {
    Ori o;
    o.push(make_rotation(Z_AXIS, 0.6f));
    o.upsample(3);
    trail.record(o);
  }

  std::vector<float> gts;
  deep_tween(trail,
             [&](const Quaternion &, float gt) { gts.push_back(gt); });

  const float expected[] = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
  HS_EXPECT_EQ(gts.size(), static_cast<size_t>(5));
  for (size_t i = 0; i < gts.size(); ++i)
    HS_EXPECT_NEAR(gts[i], expected[i], 1e-6f);
}

/**
 * @brief Verifies a single-sample VectorTrail reads t = 1.0 (the lone trail
 * head), while a multi-sample sweep ramps 0 -> 1 oldest -> newest.
 * @details A freshly spawned or dying particle holds one history point; that
 * lone sample mirrors tween(Orientation) for a lone snapshot.
 */
inline void test_tween_vectortrail_single_sample_reaches_one() {
  Animation::VectorTrail<8> trail;
  trail.record(Vector(1, 0, 0));

  std::vector<float> ts;
  tween(trail, [&](const Vector &, float t) { ts.push_back(t); });
  HS_EXPECT_EQ(ts.size(), static_cast<size_t>(1));
  HS_EXPECT_NEAR(ts.back(), 1.0f, 1e-6f);

  trail.record(Vector(0, 1, 0));
  trail.record(Vector(0, 0, 1));
  ts.clear();
  tween(trail, [&](const Vector &, float t) { ts.push_back(t); });
  HS_EXPECT_EQ(ts.size(), static_cast<size_t>(3));
  HS_EXPECT_NEAR(ts.front(), 0.0f, 1e-6f); // oldest = tail
  HS_EXPECT_NEAR(ts.back(), 1.0f, 1e-6f);  // newest = head
}

/**
 * @brief Verifies QuantizedVectorTrail round-trips unit vectors within the
 * snorm16 error bound (<= 1/65534 per component, exact at 0 and ±1), clamps
 * out-of-domain components, and keeps Trail's oldest-first ring semantics.
 */
inline void test_quantized_vector_trail_roundtrip_and_ring() {
  constexpr float QUANT_ERR = 1.0f / 65534.0f;

  Animation::QuantizedVectorTrail<8> trail;
  trail.record(Vector(1, 0, 0));
  trail.record(Vector(0, -1, 0));
  HS_EXPECT_EQ(trail.get(0).x, 1.0f);
  HS_EXPECT_EQ(trail.get(0).y, 0.0f);
  HS_EXPECT_EQ(trail.get(1).y, -1.0f);

  for (int i = 0; i < 32; ++i) {
    Vector v = Vector::from_spherical(0.37f + 0.19f * i, 0.11f + 0.09f * i);
    trail.record(v);
    Vector d = trail.get(trail.length() - 1);
    HS_EXPECT_NEAR(d.x, v.x, QUANT_ERR);
    HS_EXPECT_NEAR(d.y, v.y, QUANT_ERR);
    HS_EXPECT_NEAR(d.z, v.z, QUANT_ERR);
  }

  trail.clear();
  trail.record(Vector(1.5f, -2.0f, 0.25f));
  HS_EXPECT_EQ(trail.get(0).x, 1.0f);
  HS_EXPECT_EQ(trail.get(0).y, -1.0f);
  HS_EXPECT_NEAR(trail.get(0).z, 0.25f, QUANT_ERR);

  Animation::QuantizedVectorTrail<4> ring;
  for (int i = 0; i < 6; ++i)
    ring.record(Vector(0, 0, 0.1f * i));
  HS_EXPECT_EQ(ring.length(), static_cast<size_t>(4));
  HS_EXPECT_NEAR(ring.get(0).z, 0.2f, QUANT_ERR); // oldest retained = 3rd
  HS_EXPECT_NEAR(ring.get(3).z, 0.5f, QUANT_ERR); // newest = last recorded
  ring.expire();
  HS_EXPECT_EQ(ring.length(), static_cast<size_t>(3));
  HS_EXPECT_NEAR(ring.get(0).z, 0.3f, QUANT_ERR);

  std::vector<float> ts;
  tween(ring, [&](const Vector &, float t) { ts.push_back(t); });
  HS_EXPECT_EQ(ts.size(), static_cast<size_t>(3));
  HS_EXPECT_NEAR(ts.front(), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(ts[1], 0.5f, 1e-6f);
  HS_EXPECT_NEAR(ts.back(), 1.0f, 1e-6f);
}

// ============================================================================
// MeshMorph (nearest-vertex SLERP + crossfade invariant)
// ----------------------------------------------------------------------------
// MeshMorph maps each destination vertex to its nearest source vertex, SLERPs
// between them by an eased alpha, and crossfades the two shaded meshes with
// op_A = 1 - alpha. Morphing an octahedron onto itself makes the nearest-vertex
// map the identity, so the SLERP output equals the destination at every alpha.
// ============================================================================

/**
 * @brief Builds a hand-rolled octahedron (6 axis vertices, 8 triangular faces).
 * @param mesh Output mesh to populate with the octahedron geometry.
 * @param arena Arena backing the mesh's vertex/face/face-count buffers.
 * @details Vertices are well separated so the nearest-vertex correspondence is
 * unambiguous.
 */
inline void build_octahedron(PolyMesh &mesh, Arena &arena) {
  static const Vector verts[6] = {{1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                                  {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};
  static const uint16_t tris[8][3] = {{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4},
                                      {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}};
  mesh.vertices.bind(arena, 6);
  mesh.face_counts.bind(arena, 8);
  mesh.faces.bind(arena, 24);
  for (const auto &v : verts)
    mesh.vertices.push_back(v);
  for (const auto &t : tris) {
    mesh.face_counts.push_back(3);
    mesh.faces.push_back(t[0]);
    mesh.faces.push_back(t[1]);
    mesh.faces.push_back(t[2]);
  }
}

/**
 * @brief Verifies morphing an octahedron onto itself keeps SLERP output equal
 * to the destination and partitions the two render opacities to unity.
 * @details The nearest-vertex map is the identity, so the SLERP output equals
 * the destination at every alpha, and op_A + alpha == 1 across the crossfade.
 */
inline void test_meshmorph_identity_self_map_and_crossfade() {
  static uint8_t buf[1 << 20];
  Arena arena(buf, sizeof(buf));
  PolyMesh poly;
  build_octahedron(poly, arena);
  MeshState src, dst;
  MeshOps::compile(poly, src, arena, scratch_arena_a);
  MeshOps::compile(poly, dst, arena, scratch_arena_a);

  const int duration = 8;
  float opA = -1.0f, alpha = -1.0f;
  bool mesh_b_tracks_dest = true;

  // Effect-owned lvalue callables (MeshMorph stores non-owning FunctionRefs).
  auto draw_out = [&](Canvas &, const MeshState &, float op) { opA = op; };
  auto draw_in = [&](Canvas &, const MeshState &mb, float a) {
    alpha = a;
    if (mb.vertices.size() != dst.vertices.size()) {
      mesh_b_tracks_dest = false;
      return;
    }
    for (size_t i = 0; i < mb.vertices.size(); ++i) {
      if ((mb.vertices[i] - dst.vertices[i]).magnitude() > 1e-3f)
        mesh_b_tracks_dest = false;
    }
  };

  Animation::MeshMorph morph(src, dst, arena, draw_out, draw_in, duration);

  bool saw_both_layers = false;
  for (int i = 0; i < duration; ++i) {
    opA = -1.0f;
    alpha = -1.0f;
    morph.step(fake_canvas());
    // On any frame where both layers render, their opacities partition unity.
    if (opA >= 0.0f && alpha >= 0.0f) {
      saw_both_layers = true;
      HS_EXPECT_NEAR(opA + alpha, 1.0f, 1e-4f);
    }
  }

  HS_EXPECT_TRUE(saw_both_layers);
  HS_EXPECT_TRUE(mesh_b_tracks_dest);
  HS_EXPECT_NEAR(alpha, 1.0f, 1e-3f); // fully faded to the incoming mesh
  HS_EXPECT_TRUE(morph.done());
}

// ============================================================================
// MeshCarousel arena compaction (compact / compact_keep_front)
// ----------------------------------------------------------------------------
// Both evacuate slots to a scratch arena, reset the persistent arena, and
// restore on scope exit. compact() keeps both slots; compact_keep_front() drops
// the back slot and runs after_reset before the front restore.
// ============================================================================

/**
 * @brief Verifies compact() evacuates and restores BOTH slots across the
 * persistent-arena reset, preserving each slot's geometry.
 */
inline void test_meshcarousel_compact_retains_both_slots() {
  persistent_arena.reset();
  static uint8_t polybuf[1 << 14];
  Arena polyarena(polybuf, sizeof(polybuf));
  PolyMesh poly;
  build_octahedron(poly, polyarena);

  MeshCarousel<Segue::Crossfade> carousel;
  MeshOps::compile(poly, carousel.slot(0), persistent_arena, scratch_arena_a);
  MeshOps::compile(poly, carousel.slot(1), persistent_arena, scratch_arena_a);
  const size_t v0 = carousel.slot(0).vertices.size();
  const size_t v1 = carousel.slot(1).vertices.size();
  HS_EXPECT_EQ(v0, static_cast<size_t>(6));
  HS_EXPECT_EQ(v1, static_cast<size_t>(6));

  carousel.compact();

  HS_EXPECT_TRUE(carousel.slot(0).is_bound());
  HS_EXPECT_TRUE(carousel.slot(1).is_bound());
  HS_EXPECT_EQ(carousel.slot(0).vertices.size(), v0);
  HS_EXPECT_EQ(carousel.slot(1).vertices.size(), v1);
  HS_EXPECT_VEC(carousel.slot(0).vertices[0], Vector(1, 0, 0), 1e-5f);
  HS_EXPECT_VEC(carousel.slot(1).vertices[0], Vector(1, 0, 0), 1e-5f);
}

/**
 * @brief Verifies compact_keep_front() drops the back slot, restores the front,
 * and runs after_reset BEFORE the front restore — a re-bake lands beneath the
 * restored front geometry in the arena.
 */
inline void test_meshcarousel_compact_keep_front_drops_back() {
  persistent_arena.reset();
  static uint8_t polybuf[1 << 14];
  Arena polyarena(polybuf, sizeof(polybuf));
  PolyMesh poly;
  build_octahedron(poly, polyarena);

  MeshCarousel<Segue::Crossfade> carousel; // front slot 0
  MeshOps::compile(poly, carousel.slot(0), persistent_arena, scratch_arena_a);
  MeshOps::compile(poly, carousel.slot(1), persistent_arena, scratch_arena_a);
  const size_t v_front = carousel.current().vertices.size();

  bool after_reset_ran = false;
  const void *bake_ptr = nullptr;
  carousel.compact_keep_front([&](Arena &a) {
    after_reset_ran = true;
    bake_ptr = a.allocate(64);
  });

  HS_EXPECT_TRUE(after_reset_ran);
  HS_EXPECT_TRUE(carousel.current().is_bound());
  HS_EXPECT_EQ(carousel.current().vertices.size(), v_front);
  HS_EXPECT_VEC(carousel.current().vertices[0], Vector(1, 0, 0), 1e-5f);
  HS_EXPECT_FALSE(carousel.slot(1).is_bound());
  HS_EXPECT_GT(reinterpret_cast<uintptr_t>(&carousel.current().vertices[0]),
               reinterpret_cast<uintptr_t>(bake_ptr));
}

// ============================================================================
// ColorWipe
// ----------------------------------------------------------------------------
// Snapshots the source palette's keys on the first step (mirroring Transition),
// then OKLCH-lerps them toward the target snapshot. The three RNG-free key
// constructors keep these palettes deterministic.
// ============================================================================

/**
 * @brief Builds a fixed three-key palette with no RNG draws.
 * @param ka First key color.
 * @param kb Second key color.
 * @param kc Third key color.
 * @return A STRAIGHT-gradient palette over the supplied keys.
 */
inline GenerativePalette make_palette(CPixel ka, CPixel kb, CPixel kc) {
  return GenerativePalette(GradientShape::STRAIGHT, ka, kb, kc);
}

/**
 * @brief Verifies ColorWipe lands the source palette's keys exactly on the
 * target keys at completion and reports done() only on the final frame.
 */
inline void test_colorwipe_reaches_target_keys() {
  GenerativePalette from = make_palette(CPixel(10, 20, 30), CPixel(40, 50, 60),
                                        CPixel(70, 80, 90));
  GenerativePalette to = make_palette(CPixel(200, 0, 0), CPixel(0, 200, 0),
                                      CPixel(0, 0, 200));
  GenerativePalette::Snapshot target = to.snapshot();

  const int duration = 6;
  Animation::ColorWipe wipe(from, to, duration, ease_linear);
  HS_EXPECT_FALSE(wipe.done());

  for (int i = 0; i < duration - 1; ++i)
    wipe.step(fake_canvas());
  HS_EXPECT_FALSE(wipe.done()); // not done until t == duration

  wipe.step(fake_canvas()); // t == duration: amount == 1 -> exact target keys
  HS_EXPECT_TRUE(wipe.done());
  HS_EXPECT_EQ(static_cast<int>(from.snapshot().a.r), static_cast<int>(target.a.r));
  HS_EXPECT_EQ(static_cast<int>(from.snapshot().a.g), static_cast<int>(target.a.g));
  HS_EXPECT_EQ(static_cast<int>(from.snapshot().b.b), static_cast<int>(target.b.b));
  HS_EXPECT_EQ(static_cast<int>(from.snapshot().c.r), static_cast<int>(target.c.r));
}

/**
 * @brief Verifies the start keys are snapshotted on the first step, not at
 * construction: editing the source palette before the first step is honored.
 */
inline void test_colorwipe_snapshots_on_first_step() {
  GenerativePalette from = make_palette(CPixel(10, 10, 10), CPixel(10, 10, 10),
                                        CPixel(10, 10, 10));
  GenerativePalette to = make_palette(CPixel(250, 250, 250), CPixel(250, 250, 250),
                                      CPixel(250, 250, 250));
  const int duration = 4;
  Animation::ColorWipe wipe(from, to, duration, ease_linear);

  // Edit the source after construction but before the first step: the snapshot
  // taken on the first step must capture THIS value, so the midpoint sits
  // between 200 and 250, not between 10 and 250.
  from = make_palette(CPixel(200, 200, 200), CPixel(200, 200, 200),
                      CPixel(200, 200, 200));

  wipe.step(fake_canvas()); // t=1: amount 0.25, snapshot captured here
  HS_EXPECT_GT(static_cast<int>(from.snapshot().a.r), 200);
}

// ============================================================================
// Mobius warps (b-coefficient drivers)
// ----------------------------------------------------------------------------
// Each warp eases a normalized progress to an angle and writes Mobius param b.
// At completion (progress == 1) angle == 2π, giving exact closed-form b values.
// ============================================================================

/**
 * @brief Verifies MobiusWarp eases param b around the unit circle, closing to
 * b == 0 at completion and reporting done() only on the final frame.
 */
inline void test_mobiuswarp_closes_at_completion() {
  MobiusParams params;
  const float scale = 0.4f;
  const int duration = 8;
  Animation::MobiusWarp warp(params, scale, duration, /*repeat=*/false,
                             ease_linear);
  HS_EXPECT_FALSE(warp.done());

  warp.step(fake_canvas()); // t=1: b lifts off the origin
  HS_EXPECT_GT(std::abs(params.b.re) + std::abs(params.b.im), 1e-4f);

  for (int i = 1; i < duration - 1; ++i)
    warp.step(fake_canvas());
  HS_EXPECT_FALSE(warp.done());

  warp.step(fake_canvas()); // t=duration: angle 2π -> b back to 0
  HS_EXPECT_TRUE(warp.done());
  HS_EXPECT_NEAR(params.b.re, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(params.b.im, 0.0f, 1e-4f);
}

/**
 * @brief Verifies bind_scale makes step() read the live referent instead of the
 * captured construction-time scale.
 */
inline void test_mobiuswarp_bind_scale_reads_live() {
  MobiusParams params;
  float live = 1.0f;
  const int duration = 4;
  Animation::MobiusWarp warp(params, /*scale=*/0.0f, duration, /*repeat=*/false,
                             ease_linear);
  warp.bind_scale(live);
  warp.step(fake_canvas()); // captured scale is 0, so any motion comes from live
  HS_EXPECT_GT(std::abs(params.b.re) + std::abs(params.b.im), 1e-4f);
}

/**
 * @brief Verifies MobiusWarpCircular traces param b on the |b| == scale circle,
 * landing at (scale, 0) at completion and reporting the done() boundary.
 */
inline void test_mobiuswarp_circular_traces_radius() {
  MobiusParams params;
  const float scale = 0.3f;
  const int duration = 8;
  Animation::MobiusWarpCircular warp(params, scale, duration, /*repeat=*/false,
                                     ease_linear);
  warp.step(fake_canvas()); // |b| sits on the scale-radius circle every frame
  HS_EXPECT_NEAR(std::sqrt(params.b.re * params.b.re + params.b.im * params.b.im),
                 scale, 1e-4f);

  for (int i = 1; i < duration; ++i)
    warp.step(fake_canvas());
  HS_EXPECT_TRUE(warp.done());
  // angle 2π: b.re == scale, b.im == 0.
  HS_EXPECT_NEAR(params.b.re, scale, 1e-4f);
  HS_EXPECT_NEAR(params.b.im, 0.0f, 1e-4f);
}

/**
 * @brief Verifies MobiusWarpEvolving modulates all eight coefficients within
 * ±scale of their captured baseline and, being perpetual, never reports done().
 */
inline void test_mobiuswarp_evolving_bounded_and_perpetual() {
  hs::random().seed(1337);
  MobiusParams params(2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 3.0f, 0.0f);
  const float scale = 0.5f;
  Animation::MobiusWarpEvolving warp(params, scale, /*speed=*/0.05f);
  HS_EXPECT_FALSE(warp.done());

  const MobiusParams base = params;
  for (int i = 0; i < 50; ++i) {
    warp.step(fake_canvas());
    HS_EXPECT_FALSE(warp.done()); // perpetual: duration -1
    const float bound = scale + 1e-4f;
    HS_EXPECT_LE(std::abs(params.a.re - base.a.re), bound);
    HS_EXPECT_LE(std::abs(params.b.im - base.b.im), bound);
    HS_EXPECT_LE(std::abs(params.d.re - base.d.re), bound);
    HS_EXPECT_TRUE(std::isfinite(params.c.im));
  }
}

// ============================================================================
// Ripple / Noise (the Animation:: classes, distinct from the pure transforms)
// ============================================================================

/**
 * @brief Verifies a Ripple seeds its center, ramps amplitude up then down across
 * its life, pins amplitude to 0 at the duration boundary, and reports done()
 * only on the final frame.
 */
inline void test_ripple_envelope_and_done_boundary() {
  RippleParams params;
  params.amplitude = 2.0f; // peak captured at construction
  const Vector center(0.0f, 1.0f, 0.0f);
  const int duration = 20;
  Animation::Ripple ripple(params, center, /*speed=*/0.2f, duration);

  // Construction zeroes the live amplitude and seats the center.
  HS_EXPECT_NEAR(params.amplitude, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(params.center.y, 1.0f, 1e-6f);

  std::vector<float> amps;
  float prev_phase = params.phase;
  bool phase_advances = true;
  for (int i = 0; i < duration; ++i) {
    ripple.step(fake_canvas());
    amps.push_back(params.amplitude);
    if (i < duration - 1 && params.phase <= prev_phase)
      phase_advances = false;
    prev_phase = params.phase;
    if (i < duration - 1)
      HS_EXPECT_FALSE(ripple.done());
  }
  HS_EXPECT_TRUE(ripple.done()); // t == duration

  HS_EXPECT_TRUE(phase_advances);
  // Envelope rises off zero then decays back to zero on the final frame.
  float peak = 0.0f;
  for (float a : amps)
    peak = std::max(peak, a);
  HS_EXPECT_GT(peak, 0.0f);
  HS_EXPECT_NEAR(amps.back(), 0.0f, 1e-6f);
}

/**
 * @brief Verifies the Noise animation publishes its frame counter as
 * params.time each step and, being perpetual, never reports done().
 */
inline void test_noise_publishes_time_and_is_perpetual() {
  NoiseParams params;
  Animation::Noise noise(params); // default duration -1
  HS_EXPECT_FALSE(noise.done());
  HS_EXPECT_NEAR(params.time, 0.0f, 1e-6f);

  for (int i = 1; i <= 5; ++i) {
    noise.step(fake_canvas());
    HS_EXPECT_NEAR(params.time, static_cast<float>(i), 1e-6f);
    HS_EXPECT_FALSE(noise.done());
  }
}

/**
 * @brief A seeded RandomWalk keeps its oriented vector unit-length each frame and
 * accumulates nonzero travel on the sphere.
 */
inline void test_random_walk_stays_unit_and_travels() {
  Orientation<4> o; // identity
  FastNoiseLite noise;
  Animation::RandomWalk<288, 4> walk(
      o, Y_AXIS, noise,
      Animation::RandomWalk<288, 4>::Options::Energetic(), /*seed=*/1234);

  const Vector probe = X_AXIS;
  Vector prev = o.orient(probe);
  float travel = 0.0f;
  for (int fr = 0; fr < 50; ++fr) {
    walk.step(fake_canvas());
    const Vector cur = o.orient(probe);
    HS_EXPECT_NEAR(cur.length(), 1.0f, 1e-4f);
    travel += angle_between(prev, cur);
    prev = cur;
  }
  HS_EXPECT_GT(travel, 0.0f);
}

// ============================================================================
// Previously-uncovered public animation APIs
// ============================================================================

/**
 * @brief Verifies a non-repeating RandomTimer fires exactly once, at a frame
 * within the requested inclusive [min, max] delay window.
 */
inline void test_random_timer_fires_within_range() {
  hs::random().seed(1337);
  Timeline tl;
  struct {
    int fires = 0;
    int fire_frame = -1;
    int frame = 0;
  } st; // one capture keeps the callback inside TimerFn's inplace budget
  tl.add(0, Animation::RandomTimer(3, 7, [&st](Canvas &) {
              st.fires++;
              st.fire_frame = st.frame;
            }));
  for (st.frame = 1; st.frame <= 12; ++st.frame)
    tl.step(fake_canvas());
  HS_EXPECT_EQ(st.fires, 1);
  HS_EXPECT_GE(st.fire_frame, 3);
  HS_EXPECT_LE(st.fire_frame, 7);
}

/**
 * @brief Verifies PeriodicTimer::set_period reschedules the next trigger from
 * now (t + new_period), not from the original schedule.
 * @details Starts at period 5 (next trigger t=5); after one frame the period is
 * shortened to 3, so reset() moves the trigger to t=1+3=4.
 */
inline void test_periodic_timer_set_period_reschedules_from_now() {
  struct {
    int fire_frame = -1;
    int fires = 0;
    int frame = 0;
  } st; // one capture keeps the callback inside TimerFn's inplace budget
  Animation::PeriodicTimer timer(
      5,
      [&st](Canvas &) {
        st.fires++;
        st.fire_frame = st.frame;
      },
      /*repeat=*/true);
  st.frame = 1;
  timer.step(fake_canvas()); // t=1, no trigger (next=5)
  timer.set_period(3);       // reschedule: next = 1 + 3 = 4
  for (st.frame = 2; st.frame <= 4; ++st.frame)
    timer.step(fake_canvas());
  HS_EXPECT_EQ(st.fires, 1);
  HS_EXPECT_EQ(st.fire_frame, 4);
}

/**
 * @brief Verifies MobiusFlow::step keeps the transform's a·d product at unity
 * (a and d are conjugate-reciprocal) while actually moving the parameters.
 */
inline void test_mobiusflow_step_preserves_unit_product() {
  MobiusParams params;
  const float rings = 2.0f, lines = 4.0f;
  const int duration = 8;
  Animation::MobiusFlow flow(params, rings, lines, duration, /*repeat=*/false);
  for (int i = 0; i < duration; ++i) {
    flow.step(fake_canvas());
    const float re = params.a.re * params.d.re - params.a.im * params.d.im;
    const float im = params.a.re * params.d.im + params.a.im * params.d.re;
    HS_EXPECT_NEAR(re, 1.0f, 1e-4f);
    HS_EXPECT_NEAR(im, 0.0f, 1e-4f);
  }
  HS_EXPECT_GT(std::abs(params.a.im), 1e-3f); // moved off the identity
}

/**
 * @brief Verifies a registered emitter runs on every ParticleSystem::step and
 * can spawn into the pool.
 */
inline void test_particle_system_emitter_dispatch() {
  hs::random().seed(1337);
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 4> ps;
  ps.init(arena);
  int calls = 0;
  ps.add_emitter([&](Animation::ParticleSystem<32, 4> &sys) {
    calls++;
    sys.spawn(Vector(1, 0, 0), Vector(0, 0, 0), 0);
  });
  ps.step(fake_canvas());
  HS_EXPECT_EQ(calls, 1);
  HS_EXPECT_GT(static_cast<int>(ps.active()), 0);
  ps.step(fake_canvas());
  HS_EXPECT_EQ(calls, 2); // runs every frame
}

/**
 * @brief Verifies Motion::set_duration reanchors the baseline so a mid-path
 * duration change advances the head incrementally rather than teleporting it.
 * @details A great-circle path stepped 30 frames of a 60-frame loop sits at
 * phase 0.5; shortening to 120 frames must not snap the head by the phase gap
 * between the two parameterizations — the post-change step stays a small move.
 */
inline void test_motion_set_duration_reanchors_no_teleport() {
  using Ori = Orientation<16>;
  ProceduralPath path;
  path.f = [](float t) {
    float a = 2.0f * PI_F * t;
    return Vector(std::cos(a), std::sin(a), 0.0f);
  };
  Ori o; // identity
  Animation::Motion<288, 16> motion(o, path, 60, /*repeat=*/true);

  const Vector probe = Z_AXIS;
  for (int i = 0; i < 30; ++i)
    motion.step(fake_canvas());
  const Vector before = o.orient(probe);

  motion.set_duration(120);
  motion.step(fake_canvas());
  const Vector after = o.orient(probe);

  // Incremental (reanchored) step is a few degrees; a teleport would be ~π/2.
  HS_EXPECT_LT(angle_between(before, after), 0.5f);
}

/**
 * @brief Runs every animation/easing test case in this module.
 * @return The module's failure count.
 */
inline int run_animation_tests() {
  hs_test::ModuleFixture fixture("animation");

  // Module-scoped fake-canvas fixture (see "Stand-in Canvas reference" above).
  FakeEffect fake_fx;
  Canvas fake_cv(fake_fx);
  fake_canvas_ptr() = &fake_cv;

  test_path_empty_returns_origin();
  test_path_endpoints_and_clamp();
  test_path_collapse_keeps_last();

  test_transition_reaches_target_linear();
  test_transition_monotonic_and_starts_from_current();
  test_transition_duration_zero_no_divide_by_zero();
  test_transition_quantized_floors_result();
  test_transition_repeat_retraverses_each_cycle();

  test_mutation_applies_function_of_eased_time();
  test_mutation_duration_zero_finite();

  test_driver_increments_and_wraps();
  test_driver_no_wrap_accumulates();
  test_driver_nan_source_does_not_poison();

  test_lerp_drives_subject_to_target();
  test_lerp_midpoint();

  test_orientation_trail_index_zero_is_oldest();
  test_orientation_trail_expire_drops_oldest();
  test_orientation_trail_clear();

  test_rotation_substeps_shared_and_tight();
  test_rotation_accumulates_subthreshold_deltas();
  test_timeline_shared_orientation_composes_motion_blur();
  test_timeline_sequences_events_by_start_frame();
  test_timeline_accepts_maximum_start_frame();
  test_timeline_repeating_animation_rewinds_each_cycle();
  test_timeline_cancel_removes_repeating_animation();
  test_timeline_compaction_preserves_later_events();
  test_timeline_then_chains_follow_up_event();
  test_repeating_timer_fires_then_each_cycle();
  test_timeline_clear_resets_state();
  test_timeline_full_guard_rejects_overflow();
  test_orientation_upsample_then_collapse();
  test_motion_repeating_does_not_drift();
  test_motion_codriven_survives_repeat_seam();

  test_particle_system_spawn_and_capacity_guard();
  test_particle_system_lifetime_boundaries();
  test_particle_system_spawn_initializes_and_steps();
  test_particle_system_expires_after_life_and_trail_drain();
  test_particle_system_attractor_kills_within_radius();
  test_particle_system_attractor_kill_radius_boundary();

  test_sprite_fade_in_plateau_fade_out_envelope();
  test_sprite_overlapping_fades_stay_continuous();
  test_sprite_paused_holds_frame();

  test_crossfade_segue_schedules_overlapping_sprite();
  test_crossfade_segue_clamps_fade_to_half_duration();
  test_sequential_segue_never_overlaps_sprites();
  test_dissolve_segue_masks_partition_the_canvas();
  test_dissolve_segue_reseeds_per_frame_and_transition();
  test_segue_base_hooks_are_identity();
  test_iris_bloom_fill_contracts_to_face_centers();
  test_lace_fill_keeps_edge_band();
  test_sweep_phase_front_ordering();
  test_terminator_sweep_orders_by_axis();
  test_terminator_sweep_fades_faces_over_fixed_frames();
  test_terminator_sweep_per_face_fade_random_in_range();
  test_shockwave_orders_by_distance_from_origin();
  test_breakdown_fades_classes_sequentially();
  test_spin_flip_warp_is_rigid();
  test_gold_convergence_grades_to_gold();

  test_deep_tween_global_t_spans_unit_interval();
  test_deep_tween_collapsed_newest_frame_reaches_one();
  test_deep_tween_all_collapsed_reaches_one();
  test_deep_tween_frames_groups_flat_emission();
  test_deep_tween_interior_motionless_frame_no_gap();
  test_tween_vectortrail_single_sample_reaches_one();
  test_quantized_vector_trail_roundtrip_and_ring();

  test_meshmorph_identity_self_map_and_crossfade();
  test_meshcarousel_compact_retains_both_slots();
  test_meshcarousel_compact_keep_front_drops_back();

  test_colorwipe_reaches_target_keys();
  test_colorwipe_snapshots_on_first_step();

  test_mobiuswarp_closes_at_completion();
  test_mobiuswarp_bind_scale_reads_live();
  test_mobiuswarp_circular_traces_radius();
  test_mobiuswarp_evolving_bounded_and_perpetual();

  test_ripple_envelope_and_done_boundary();
  test_noise_publishes_time_and_is_perpetual();

  test_random_walk_stays_unit_and_travels();

  test_random_timer_fires_within_range();
  test_periodic_timer_set_period_reschedules_from_now();
  test_mobiusflow_step_preserves_unit_product();
  test_particle_system_emitter_dispatch();
  test_motion_set_duration_reanchors_no_teleport();

  const int result = fixture.result();
  // Unpublish before fake_cv/fake_fx destruct below.
  fake_canvas_ptr() = nullptr;
  return result;
}

} // namespace animation_tests
} // namespace hs_test

