/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/animation.h and core/easing.h.
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
#include <vector>
#include "core/animation.h"
#include "core/canvas.h"
#include "core/easing.h"
#include "core/mesh.h" // PolyMesh, MeshOps::compile (MeshMorph fixtures)
#include "tests/test_harness.h"

namespace hs_test {
namespace animation_tests {

// ============================================================================
// Stand-in Canvas reference
// ----------------------------------------------------------------------------
// The animations under test (Transition/Mutation/Lerp/Driver) take a Canvas& in
// step() but never dereference it. We construct ONE genuine Canvas over a tiny
// Effect and reuse it, so the reference is always valid. A fresh effect's
// buffer_free() is true, so the ctor doesn't spin; the single dtor (queue_frame)
// runs harmlessly when the module ends.
//
// LIFETIME: the fixture is owned by run_animation_tests() (see bottom of file)
// and reached through fake_canvas_ptr(), NOT a process-lifetime static. A static
// FakeEffect would keep the single-live-Effect guard (canvas.h s_alive) set for
// the rest of the process — its dtor only runs at exit — and trip that guard when
// the NEXT module's first Effect is constructed in an all-modules run (the bare
// `run_tests` binary; the ctest gate runs one module per process and never saw
// it). Scoping the fixture to the module destroys its Effect at module end and
// releases the guard cleanly.
//
// COUPLING (load-bearing): nothing a test does with this canvas may invoke
// queue_frame() — i.e. no test may construct (and destruct) a *second* Canvas on
// this shared FakeEffect during the run. queue_frame() advances next_, so a
// second Canvas's ctor would then see buffer_free()==false and spin on a display
// ISR that never runs in the host harness (a hang). The only queue_frame() that
// may fire is the single dtor one at module end, after all tests pass. If a
// future animation test needs to actually queue a frame, give it its own local
// Canvas/Effect rather than reusing fake_canvas().
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
  /**
   * @brief Reports that no background is drawn.
   * @return Always false.
   */
  bool show_bg() const override { return false; }
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
  // Any t on an empty path is the origin guard.
  Vector v1 = p.get_point(1.0f);
  HS_EXPECT_NEAR(v1.x, 0.0f, 0.0f);
}

/**
 * @brief Verifies get_point hits the exact endpoints at t=0/1, clamps to back()
 * past t=1, and interpolates linearly between samples.
 */
inline void test_path_endpoints_and_clamp() {
  Path<32> p;
  // Plot a straight ramp along X: f(s) = (s, 0, 0), domain 1, 4 samples,
  // linear easing. append_segment samples t/samples in [0,1] -> x in [0,1].
  p.append_segment([](float s) { return Vector(s, 0.0f, 0.0f); }, 1.0f, 4,
                   ease_mid);

  Vector start = p.get_point(0.0f);
  HS_EXPECT_NEAR(start.x, 0.0f, 1e-5f);

  Vector end = p.get_point(1.0f);
  HS_EXPECT_NEAR(end.x, 1.0f, 1e-5f);

  // t beyond 1.0 clamps to the last point (i >= size-1 -> back()).
  Vector over = p.get_point(2.0f);
  HS_EXPECT_NEAR(over.x, end.x, 1e-5f);

  // Midpoint interpolates linearly along the ramp.
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
                   ease_mid);
  Vector last_before = p.get_point(1.0f);
  p.collapse();
  // After collapse there is a single point; get_point(anything) == that point.
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
  Animation::Transition tr(v, 100.0f, duration, ease_mid);
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
  float v = 5.0f; // 'from' is captured from the live value at t==0
  const int duration = 8;
  Animation::Transition tr(v, 25.0f, duration, ease_mid);
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
  // AnimationBase coerces duration 0 -> 1, so the t/duration ratio is finite.
  float v = 0.0f;
  Animation::Transition tr(v, 42.0f, 0, ease_mid);
  tr.step(fake_canvas());
  HS_EXPECT_TRUE(std::isfinite(v));
  HS_EXPECT_NEAR(v, 42.0f, 1e-3f); // one step of a 1-frame transition completes
  HS_EXPECT_TRUE(tr.done());
}

/**
 * @brief Verifies a quantized Transition floors the eased value, landing on an
 * integer.
 */
inline void test_transition_quantized_floors_result() {
  float v = 0.0f;
  const int duration = 4;
  Animation::Transition tr(v, 3.7f, duration, ease_mid, /*quantized=*/true);
  for (int i = 0; i < duration; ++i) {
    tr.step(fake_canvas());
  }
  // Final eased value 3.7 floored to 3.0.
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
  Animation::Transition tr(v, 10.0f, duration, ease_mid, /*quantized=*/false,
                           /*repeat=*/true);
  for (int i = 0; i < duration; ++i)
    tr.step(fake_canvas());
  HS_EXPECT_NEAR(v, 10.0f, 1e-3f); // first cycle reaches the target

  tr.rewind(); // what Timeline does for a repeating animation
  tr.step(fake_canvas());
  HS_EXPECT_LT(v, 10.0f); // second cycle leaves the target
  for (int i = 1; i < duration; ++i)
    tr.step(fake_canvas());
  HS_EXPECT_NEAR(v, 10.0f, 1e-3f); // and lands on it again
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
  // f(e) = e * 10 ; with linear easing, at completion e==1 -> v==10.
  Animation::Mutation m(
      v, [](float e) { return e * 10.0f; }, duration, ease_mid);
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
      v, [](float e) { return e; }, 0, ease_mid);
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
  Animation::Lerp l(subject, start, target, duration, ease_mid);
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
  Animation::Lerp l(subject, start, target, duration, ease_mid);
  l.step(fake_canvas()); // t=1 -> progress 0.25 -> 12.5
  l.step(fake_canvas()); // t=2 -> progress 0.50 -> 15.0
  HS_EXPECT_NEAR(subject.value, 15.0f, 1e-3f);
}

// ============================================================================
// OrientationTrail
// ----------------------------------------------------------------------------
// Index 0 is the OLDEST snapshot, length()-1 the newest (matching the JS
// simulator's trail ordering):
// OrientationTrail::record() calls snapshots.push_back(), and
// StaticCircularBuffer::operator[](0) returns buffer[head] — the first element
// still live. This test pins that behavior.
// ============================================================================

/**
 * @brief Verifies recorded snapshots are ordered oldest-first: index 0 is the
 * first recorded, the last index the newest.
 */
inline void test_orientation_trail_index_zero_is_oldest() {
  Animation::OrientationTrail<Orientation<4>, 8> trail;

  Orientation<4> a, b, c;
  a.set(make_rotation(Vector(0, 1, 0), 0.1f));
  b.set(make_rotation(Vector(0, 1, 0), 0.2f));
  c.set(make_rotation(Vector(0, 1, 0), 0.3f));

  trail.record(a);
  trail.record(b);
  trail.record(c);
  HS_EXPECT_EQ(trail.length(), static_cast<size_t>(3));

  // Index 0 == first recorded (a, the OLDEST); last index == newest (c).
  HS_EXPECT_NEAR(std::abs(dot(trail.get(0).get(), a.get())), 1.0f, 1e-4f);
  HS_EXPECT_NEAR(std::abs(dot(trail.get(2).get(), c.get())), 1.0f, 1e-4f);
}

/**
 * @brief Verifies expire() removes the oldest snapshot, leaving the survivor at
 * index 0.
 */
inline void test_orientation_trail_expire_drops_oldest() {
  Animation::OrientationTrail<Orientation<4>, 8> trail;
  Orientation<4> a, b;
  a.set(make_rotation(Vector(0, 1, 0), 0.1f));
  b.set(make_rotation(Vector(0, 1, 0), 0.2f));
  trail.record(a);
  trail.record(b);
  trail.expire(); // removes the oldest (a)
  HS_EXPECT_EQ(trail.length(), static_cast<size_t>(1));
  // The survivor is b, now at index 0.
  HS_EXPECT_NEAR(std::abs(dot(trail.get(0).get(), b.get())), 1.0f, 1e-4f);
}

/**
 * @brief Verifies clear() empties the trail.
 */
inline void test_orientation_trail_clear() {
  Animation::OrientationTrail<Orientation<4>, 8> trail;
  Orientation<4> a;
  trail.record(a);
  trail.record(a);
  trail.clear();
  HS_EXPECT_EQ(trail.length(), static_cast<size_t>(0));
}

// ============================================================================
// easing.h sanity
// ============================================================================

/**
 * @brief Verifies every in/out easing variant anchors at f(0)==0 and f(1)==1.
 */
inline void test_easing_in_out_endpoints() {
  struct {
    const char *name;
    EasingFn fn;
  } fns[] = {
      {"ease_in_out_cubic", ease_in_out_cubic},
      {"ease_in_out_sin", ease_in_out_sin},
      {"ease_in_sin", ease_in_sin},
      {"ease_out_sin", ease_out_sin},
      {"ease_in_cubic", ease_in_cubic},
      {"ease_in_circ", ease_in_circ},
      {"ease_mid", ease_mid},
      {"ease_out_expo", ease_out_expo},
      {"ease_out_circ", ease_out_circ},
      {"ease_out_cubic", ease_out_cubic},
  };
  for (const auto &e : fns) {
    HS_EXPECT_NEAR(e.fn(0.0f), 0.0f, 1e-3f);
    HS_EXPECT_NEAR(e.fn(1.0f), 1.0f, 1e-3f);
  }
}

/**
 * @brief Verifies every easing function stays finite across its documented
 * [0,1] domain.
 */
inline void test_easing_output_finite_over_range() {
  EasingFn fns[] = {ease_in_out_cubic, ease_in_out_sin, ease_in_sin,
                    ease_out_sin,        ease_in_cubic,   ease_in_circ,
                    ease_mid,            ease_out_expo,   ease_out_circ,
                    ease_out_cubic,      ease_out_elastic};
  // Sample the documented domain [0,1] exactly. NOTE: float accumulation can
  // drift just past 1.0, and ease_out_circ/ease_in_circ are sqrt(1-(...)^2)
  // forms that go NaN for t>1 (out of their defined domain), so we step with
  // an integer index instead of accumulating a float.
  for (EasingFn fn : fns) {
    for (int i = 0; i <= 20; ++i) {
      float t = static_cast<float>(i) / 20.0f;
      HS_EXPECT_TRUE(std::isfinite(fn(t)));
    }
  }
}

/**
 * @brief Verifies ease_out_elastic anchors at f(0)==0 and f(1)==1 even though
 * it overshoots in between.
 * @details ease_out_elastic is not part of an in/out pair, so it is checked
 * separately from the in/out endpoint suite.
 */
inline void test_easing_elastic_anchors() {
  HS_EXPECT_NEAR(ease_out_elastic(0.0f), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(ease_out_elastic(1.0f), 1.0f, 1e-3f);
}

// ============================================================================
// Borrow-contract guards (compile-time)
// ----------------------------------------------------------------------------
// Motion and MeshMorph store a non-owning borrow (Fn capturing &path_obj /
// StoredFunctionRef) of effect-owned state that must outlive the Timeline.
// Binding to a temporary is a compile error instead of a silent dangle: Motion
// via a `= delete` rvalue ctor overload, MeshMorph via its StoredFunctionRef
// draw-callback parameter type (which itself rejects rvalues). These
// static_asserts lock that contract: an effect-owned lvalue must be accepted, a
// temporary must be rejected. They live at namespace scope (checked at compile
// time; no runner call needed). If a refactor turns the guard into a no-op,
// this fails to compile.
namespace borrow_guard {
/** @brief Orientation alias used by the borrow-contract static_asserts. */
using Ori = Orientation<16>;
/** @brief Mesh draw-callable alias used by the MeshMorph borrow asserts. */
using DrawFn = Fn<void(Canvas &, const MeshState &, float), 8>;

// Motion: effect-owned lvalue path OK; temporary path rejected.
static_assert(std::is_constructible_v<Animation::Motion<288, 16>, Ori &,
                                      ProceduralPath &, int>,
              "Motion must accept an lvalue (effect-owned) path");
static_assert(!std::is_constructible_v<Animation::Motion<288, 16>, Ori &,
                                       ProceduralPath &&, int>,
              "Motion must REJECT a temporary path (would dangle)");

// MeshMorph: effect-owned lvalue callables OK; a temporary in EITHER slot
// rejected.
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

// Lerp stores raw pointers to subject/start/target and dereferences them every
// frame; start/target are `const T&` and so could bind a temporary. The
// deleted rvalue overloads must make an lvalue start+target OK and a temporary
// in either slot a compile error.
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

// MobiusFlow stores num_rings/num_lines as reference_wrappers and re-reads them
// every frame (live-tracking); they are `const float&` and so could bind a
// temporary. The deleted rvalue overloads must make lvalue scalars OK and a
// temporary in either slot a compile error. `params` is a non-const ref and
// already rejects temporaries on its own.
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
  // Tight ceil: N*MAX needs exactly N subdivisions (no extra +1).
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX, MAX), 1);
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX * 3.0f, MAX), 3);
  // A fractional overshoot rounds up.
  HS_EXPECT_EQ(Animation::rotation_substeps(MAX * 3.2f, MAX), 4);
  // Each of the (substeps) sub-intervals stays within MAX (the invariant the
  // trail sizing must guarantee).
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
 * turns over many frames. ease_mid is linear here, so the per-frame delta is
 * total_angle / duration.
 */
inline void test_rotation_accumulates_subthreshold_deltas() {
  using Ori = Orientation<16>;
  Ori o; // identity
  // 0.05 rad over 1000 frames => 5e-5 rad/frame, half of TOLERANCE, so every
  // single frame's raw delta is below the early-out threshold.
  Animation::Rotation<288, 16> rot(o, Z_AXIS, 0.05f, 1000, ease_mid);
  for (int i = 0; i < 20; ++i)
    rot.step(fake_canvas());
  // After 20 frames the accumulated sweep is ~1e-3 rad; a Z rotation sends +X
  // toward +Y.
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
  // Two world-space rotations about the SAME axis, each a quarter turn,
  // completing in one frame so the single step sweeps the full angle.
  tl.add(0, Animation::Rotation<288, 16>(o, Z_AXIS, PI_F / 2, 1, ease_mid));
  tl.add(0, Animation::Rotation<288, 16>(o, Z_AXIS, PI_F / 2, 1, ease_mid));
  tl.step(fake_canvas());

  // A motion-blur trail, not a single collapsed frame.
  HS_EXPECT_GE(o.length(), 2);
  // Oldest sub-frame is the pre-frame orientation (identity): +X stays +X.
  Vector oldest = o.orient(X_AXIS, 0);
  HS_EXPECT_NEAR(oldest.x, 1.0f, 1e-3f);
  HS_EXPECT_NEAR(oldest.y, 0.0f, 1e-3f);
  // Newest sub-frame reflects BOTH rotations (a half turn): +X -> -X.
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
  Timeline tl; // ctor clears the global cursors; dtor releases the live guard
  float a = 0.0f, b = 0.0f;
  tl.add(0, Animation::Transition(a, 10.0f, 2, ease_mid)); // starts now
  tl.add(3, Animation::Transition(b, 20.0f, 2, ease_mid)); // delayed 3 frames

  tl.step(fake_canvas()); // t=1: a ramps, b dormant
  HS_EXPECT_GT(a, 0.0f);
  HS_EXPECT_NEAR(b, 0.0f, 1e-6f);

  tl.step(fake_canvas()); // t=2: a completes (reaches target) and is removed
  HS_EXPECT_NEAR(a, 10.0f, 1e-3f);
  HS_EXPECT_NEAR(b, 0.0f, 1e-6f); // still dormant (start=3)
  HS_EXPECT_EQ(global_timeline_num_events, 1); // only b remains scheduled

  tl.step(fake_canvas()); // t=3: b finally starts ramping
  HS_EXPECT_GT(b, 0.0f);
  HS_EXPECT_LT(b, 20.0f);

  tl.step(fake_canvas()); // t=4: b completes
  HS_EXPECT_NEAR(b, 20.0f, 1e-3f);
  HS_EXPECT_EQ(global_timeline_num_events, 0); // both done and removed
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
                v, [](float e) { return e; }, 2, ease_mid, /*repeat=*/true));

  tl.step(fake_canvas()); // t=1 -> v = eased(0.5) = 0.5
  HS_EXPECT_NEAR(v, 0.5f, 1e-3f);
  tl.step(fake_canvas()); // t=2 -> v = eased(1.0) = 1.0, done -> rewind
  HS_EXPECT_NEAR(v, 1.0f, 1e-3f);
  tl.step(fake_canvas()); // rewound: t=1 again -> v = 0.5 (proves the rewind)
  HS_EXPECT_NEAR(v, 0.5f, 1e-3f);

  // Repeating events are never removed.
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
                              v, [](float e) { return e; }, 2, ease_mid,
                              /*repeat=*/true));
  tl.step(fake_canvas());
  tl.step(fake_canvas()); // completes a cycle, rewinds, and stays (repeating)
  HS_EXPECT_EQ(global_timeline_num_events, 1);

  h->cancel();
  tl.step(fake_canvas()); // canceled: done() && !repeats() -> removed
  HS_EXPECT_EQ(global_timeline_num_events, 0);
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
  tl.add(0, Animation::Transition(a, 10.0f, 1, ease_mid));  // completes at t=1
  tl.add(0, Animation::Transition(b, 100.0f, 5, ease_mid)); // in-flight survivor
  tl.add(0, Animation::Transition(c, 200.0f, 5, ease_mid)); // in-flight survivor
  HS_EXPECT_EQ(global_timeline_num_events, 3);

  tl.step(fake_canvas()); // t=1: a done+removed; b,c step once and shift down
  HS_EXPECT_NEAR(a, 10.0f, 1e-3f);
  HS_EXPECT_GT(b, 0.0f);
  HS_EXPECT_GT(c, 0.0f);
  HS_EXPECT_EQ(global_timeline_num_events, 2); // a gone; b,c compacted to 0,1
  float b_after1 = b, c_after1 = c;

  for (int i = 0; i < 4; ++i)
    tl.step(fake_canvas()); // t=2..5: finish the relocated survivors
  HS_EXPECT_GT(b, b_after1);            // kept advancing post-relocation
  HS_EXPECT_GT(c, c_after1);
  HS_EXPECT_NEAR(b, 100.0f, 1e-2f);
  HS_EXPECT_NEAR(c, 200.0f, 1e-2f);     // last event survived relocation intact
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
  tl.add(0, Animation::Transition(a, 10.0f, 1, ease_mid).then([&]() {
    tl.add(0, Animation::Transition(b, 20.0f, 1, ease_mid));
  }));

  tl.step(fake_canvas()); // t=1: a completes -> callback adds b (gap-filled in)
  HS_EXPECT_NEAR(a, 10.0f, 1e-3f);
  HS_EXPECT_NEAR(b, 0.0f, 1e-6f);              // b added this frame, not yet run
  HS_EXPECT_EQ(global_timeline_num_events, 1); // a removed, b now scheduled

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
  HS_EXPECT_EQ(thens, 3); // the .then() hook fires once per trigger
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
  tl.add(0, Animation::Transition(a, 10.0f, 5, ease_mid));
  tl.step(fake_canvas());
  HS_EXPECT_EQ(global_timeline_num_events, 1);
  HS_EXPECT_EQ(global_timeline_t, 1);

  tl.clear();
  HS_EXPECT_EQ(global_timeline_num_events, 0);
  HS_EXPECT_EQ(global_timeline_t, 0); // cursor rewound

  // Reusable: a fresh one-frame event schedules from t=0 and completes.
  float b = 0.0f;
  tl.add(0, Animation::Transition(b, 5.0f, 1, ease_mid));
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
  // Targets share one float — this test asserts the count guard, not values, and
  // never steps. Fill to capacity.
  for (int i = 0; i < Timeline::MAX_EVENTS; ++i)
    tl.add(0, Animation::Transition(sink, 1.0f, 10, ease_mid));
  HS_EXPECT_EQ(global_timeline_num_events, Timeline::MAX_EVENTS);

  tl.add(0, Animation::Transition(sink, 1.0f, 10, ease_mid)); // one past full
  HS_EXPECT_EQ(global_timeline_num_events,
               Timeline::MAX_EVENTS); // rejected, count unchanged
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
  HS_EXPECT_EQ(o.length(), 2);

  o.upsample(5);
  HS_EXPECT_EQ(o.length(), 5);

  // Endpoints preserved: frame 0 ~ identity (+X stays +X); frame 4 ~ the +90
  // rotation about Z (+X -> +Y).
  Vector f0 = o.orient(X_AXIS, 0);
  HS_EXPECT_NEAR(f0.x, 1.0f, 1e-3f);
  Vector f4 = o.orient(X_AXIS, 4);
  HS_EXPECT_NEAR(f4.y, 1.0f, 1e-3f);

  // SLERP is monotone: the +X component decreases across the interpolated frames.
  float prevx = 2.0f;
  for (int i = 0; i < 5; ++i) {
    float x = o.orient(X_AXIS, i).x;
    HS_EXPECT_LE(x, prevx + 1e-4f);
    prevx = x;
  }

  // collapse() keeps only the newest sub-frame (the +90 rotation survives).
  o.collapse();
  HS_EXPECT_EQ(o.length(), 1);
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

  const int duration = 40;
  ProceduralPath path;
  // A representative Comets Lissajous figure; lissajous(.,.,.,0) == +Y, so the
  // identity-start orientation places the head exactly on the path at phase 0.
  path.f = [](float t) { return lissajous(1.06f, 1.06f, 0.0f, t * 5.909f); };

  Ori o; // identity, single frame; orient(+Y) starts on the path
  const Vector node_v = Y_AXIS;

  Timeline tl;
  tl.add(0, Animation::Motion<288, 16>(o, path, duration, /*repeat=*/true));

  // Ideal internal angles between three fixed within-cycle phases.
  const int pa = 10, pb = 20, pc = 30;
  Vector Ia = path.f((float)pa / duration);
  Vector Ib = path.f((float)pb / duration);
  Vector Ic = path.f((float)pc / duration);
  float ideal_ab = angle_between(Ia, Ib);
  float ideal_ac = angle_between(Ia, Ic);
  float ideal_bc = angle_between(Ib, Ic);

  Vector Pa, Pb, Pc;
  const int cycles = 600; // ~24k frames; bug reaches multi-radian warp well before
  for (int c = 0; c < cycles; ++c) {
    for (int fr = 1; fr <= duration; ++fr) {
      tl.step(fake_canvas());
      // global_timeline_t advances with each step; the Motion's own phase is
      // (t-1 within the cycle). Sample by position in this inner loop.
      if (fr == pa) Pa = o.orient(node_v);
      if (fr == pb) Pb = o.orient(node_v);
      if (fr == pc) Pc = o.orient(node_v);
    }
  }

  // After 600 cycles the reconstructed internal angles must still match the
  // ideal within the first-cycle discretization residual (< 0.1 rad).
  HS_EXPECT_NEAR(angle_between(Pa, Pb), ideal_ab, 0.1f);
  HS_EXPECT_NEAR(angle_between(Pa, Pc), ideal_ac, 0.1f);
  HS_EXPECT_NEAR(angle_between(Pb, Pc), ideal_bc, 0.1f);

  global_timeline_num_events = 0;
  global_timeline_t = 0;
}

/**
 * @brief A co-driver sharing a repeating Motion's Orientation survives the
 * repeat seam (finding 3 / sole-ownership removal).
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
  // A closed great circle: path(0) == path(1) with a matching tangent, so
  // Motion's path frame is periodic and its per-cycle delta product is identity.
  path.f = [](float t) {
    float a = 2.0f * PI_F * t;
    return Vector(std::cos(a), std::sin(a), 0.0f);
  };

  Ori o; // identity
  Timeline tl;
  // Repeating Motion + a repeating co-driver rotation about Y, both driving `o`.
  tl.add(0, Animation::Motion<288, 16>(o, path, duration, /*repeat=*/true));
  tl.add(0, Animation::Rotation<288, 16>(o, Y_AXIS, 2.0f * PI_F, duration,
                                         ease_mid, /*repeat=*/true));

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

  // No single frame — seams included — snaps the orientation: the co-driver was
  // preserved, not clobbered. The bound clears the largest legitimate one-frame
  // step (Motion ~2pi/duration plus the co-driver's share) yet sits far below
  // the multi-radian snap the old absolute reset produced.
  HS_EXPECT_LT(max_step, 0.8f);
  // Non-vacuous: the co-driver actually moved the probe a long way over 8 cycles.
  HS_EXPECT_GT(total_travel, 5.0f);

  global_timeline_num_events = 0;
  global_timeline_t = 0;
}

// ============================================================================
// ParticleSystem
// ----------------------------------------------------------------------------
// step() only advances the base frame counter and runs emitter/physics updates;
// it never dereferences the canvas, so fake_canvas() suffices. These cover the
// three pool-state transitions the system exists to manage: spawn (+ capacity
// guard), life-expiry with trail drain, and attractor kill-radius removal.
// ============================================================================

/**
 * @brief Verifies ParticleSystem::spawn adds particles and silently drops
 * spawns once the fixed pool is at capacity.
 */
inline void test_particle_system_spawn_and_capacity_guard() {
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 4> ps; // CAPACITY = 4
  ps.init(arena);
  HS_EXPECT_EQ(static_cast<int>(ps.active_count), 0);

  ps.spawn(Vector(1, 0, 0), Vector(0, 0, 0), 0);
  ps.spawn(Vector(0, 1, 0), Vector(0, 0, 0), 1);
  HS_EXPECT_EQ(static_cast<int>(ps.active_count), 2);

  // Fill to capacity and then one past — the extra spawn is silently dropped,
  // never overrunning the fixed pool.
  ps.spawn(Vector(0, 0, 1), Vector(0, 0, 0), 2);
  ps.spawn(Vector(-1, 0, 0), Vector(0, 0, 0), 3);
  HS_EXPECT_EQ(static_cast<int>(ps.active_count), 4);
  ps.spawn(Vector(0, -1, 0), Vector(0, 0, 0), 4); // capacity is 4 — rejected
  HS_EXPECT_EQ(static_cast<int>(ps.active_count), 4);
}

/**
 * @brief Verifies a particle is removed only after its life reaches 0 and its
 * recorded trail drains to empty.
 */
inline void test_particle_system_expires_after_life_and_trail_drain() {
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 4> ps;
  // gravity 0 + zero velocity => the particle never moves; the only thing that
  // ends it is life reaching 0 and its recorded trail draining to empty.
  ps.init(arena, /*friction=*/0.85f, /*gravity=*/0.0f, /*max_life=*/3.0f);
  ps.spawn(Vector(1, 0, 0), Vector(0, 0, 0), 0);

  ps.step(fake_canvas());
  HS_EXPECT_EQ(static_cast<int>(ps.active_count), 1); // still alive

  // life (3) plus the trail length (default 8) plus slack is a hard upper bound
  // on how long a dead particle lingers while its history expires one per frame.
  for (int i = 0; i < 3 + 8 + 4; ++i)
    ps.step(fake_canvas());
  HS_EXPECT_EQ(static_cast<int>(ps.active_count), 0);
}

/**
 * @brief Verifies an attractor removes a particle inside its kill_radius via
 * the kill check rather than by life expiry.
 */
inline void test_particle_system_attractor_kills_within_radius() {
  static uint8_t buf[256 * 1024];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 4> ps;
  ps.init(arena, /*friction=*/0.85f, /*gravity=*/0.001f, /*max_life=*/600.0f);
  ps.spawn(Vector(1, 0, 0), Vector(0, 0, 0), 0);
  // Attractor co-located with the particle: distance 0 < kill_radius, so the
  // particle is removed by the kill check — not by life expiry (life is 600).
  ps.add_attractor(Vector(1, 0, 0), /*strength=*/1.0f, /*kill_radius=*/0.5f,
                   /*event_horizon=*/2.0f);
  HS_EXPECT_EQ(static_cast<int>(ps.active_count), 1);

  ps.step(fake_canvas());
  HS_EXPECT_EQ(static_cast<int>(ps.active_count), 0);
}

// ============================================================================
// Sprite (opacity envelope + paused-hold)
// ----------------------------------------------------------------------------
// Sprite::step never touches the canvas itself — it computes an opacity and
// forwards it to the user draw_fn. We capture that opacity to assert the
// fade-in/plateau/fade-out envelope and the paused-hold contract (no timer
// advance, no expiry, still drawing).
// ============================================================================

/**
 * @brief Verifies Sprite drives a fade-in -> full-opacity plateau -> fade-out
 * envelope across its duration, completing on the final frame.
 */
inline void test_sprite_fade_in_plateau_fade_out_envelope() {
  std::vector<float> ops;
  const int dur = 10, fade_in = 3, fade_out = 3;
  Animation::Sprite s(
      [&](Canvas &, float o) { ops.push_back(o); }, dur, fade_in, ease_mid,
      fade_out, ease_mid);
  for (int i = 0; i < dur; ++i)
    s.step(fake_canvas()); // observed at t = 1..10

  HS_EXPECT_EQ(ops.size(), static_cast<size_t>(dur));
  // Fade-in: first sample is partial and rising (linear ease => t/fade_in).
  HS_EXPECT_NEAR(ops[0], 1.0f / 3.0f, 1e-3f); // t=1
  HS_EXPECT_GT(ops[1], ops[0]);               // t=2 brighter
  // Plateau: fully opaque once past fade-in and before fade-out (t=3..7).
  HS_EXPECT_NEAR(ops[3], 1.0f, 1e-3f);
  HS_EXPECT_NEAR(ops[5], 1.0f, 1e-3f);
  // Fade-out: descends to fully transparent on the final frame (t=10).
  HS_EXPECT_LT(ops[8], ops[5]);
  HS_EXPECT_NEAR(ops[9], 0.0f, 1e-3f);
  HS_EXPECT_TRUE(s.done());
}

/**
 * @brief Verifies that when fade_in + fade_out exceed duration the overlapping
 * opacity envelope stays continuous (a triangle).
 * @details There must be no jump where the fade-in hands off to the fade-out;
 * the slider ranges make this configuration user-reachable.
 */
inline void test_sprite_overlapping_fades_stay_continuous() {
  std::vector<float> ops;
  const int dur = 10, fade_in = 8, fade_out = 8; // 8 + 8 > 10 => overlap
  Animation::Sprite s(
      [&](Canvas &, float o) { ops.push_back(o); }, dur, fade_in, ease_mid,
      fade_out, ease_mid);
  for (int i = 0; i < dur; ++i)
    s.step(fake_canvas()); // observed at t = 1..10

  HS_EXPECT_EQ(ops.size(), static_cast<size_t>(dur));
  // Linear ramps of slope 1/8 per frame, so no single-frame opacity step may
  // exceed ~0.125.
  float max_jump = 0.0f;
  for (size_t i = 1; i < ops.size(); ++i)
    max_jump = std::max(max_jump, std::abs(ops[i] - ops[i - 1]));
  HS_EXPECT_LT(max_jump, 0.2f);
  // Sanity: it still rises to a peak then falls (a triangle), not flat zero.
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
  // No fade in/out => held opacity is full; duration 3 would normally expire.
  Animation::Sprite s(
      [&](Canvas &, float o) { draws++; last_op = o; }, /*duration=*/3,
      /*fade_in=*/0, ease_mid, /*fade_out=*/0, ease_mid, &paused);

  for (int i = 0; i < 10; ++i)
    s.step(fake_canvas());
  HS_EXPECT_FALSE(s.done());          // timer never advanced while paused
  HS_EXPECT_EQ(draws, 10);            // but it kept drawing every frame
  HS_EXPECT_NEAR(last_op, 1.0f, 1e-6f); // at its held (full) opacity

  // Unpausing lets the timer advance and the sprite finally completes.
  paused = false;
  for (int i = 0; i < 3; ++i)
    s.step(fake_canvas());
  HS_EXPECT_TRUE(s.done());
}

// ============================================================================
// deep_tween (global_t span)
// ----------------------------------------------------------------------------
// deep_tween walks an OrientationTrail's frames and, within each frame, its
// sub-frames, emitting a global parameter t in [0,1] across the whole trail. It
// skips sub-frame 0 of every frame after the first (the shared boundary), so the
// emitted count is M + (N-1)*(M-1) for N frames of M sub-frames. We assert the
// span endpoints, monotonicity, and that exact count.
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
    Ori o;                                       // identity, length 1
    o.push(make_rotation(Z_AXIS, 0.3f * (k + 1))); // length 2
    o.upsample(M);                               // length M sub-frames
    trail.record(o);
  }

  std::vector<float> gts;
  deep_tween(trail,
             [&](const Quaternion &, float gt) { gts.push_back(gt); });

  // M for the first frame + (M-1) for each subsequent frame (sub-frame 0 of
  // frames 1..N-1 is the skipped shared boundary).
  HS_EXPECT_EQ(gts.size(), static_cast<size_t>(M + (N - 1) * (M - 1)));
  HS_EXPECT_NEAR(gts.front(), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(gts.back(), 1.0f, 1e-6f);
  for (size_t i = 1; i < gts.size(); ++i)
    HS_EXPECT_GE(gts[i], gts[i - 1] - 1e-6f); // non-decreasing across the trail
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
  Ori still; // identity, length 1 — no motion this frame
  trail.record(still);

  std::vector<float> gts;
  deep_tween(trail,
             [&](const Quaternion &, float gt) { gts.push_back(gt); });

  HS_EXPECT_GE(gts.size(), static_cast<size_t>(1));
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
    Ori still; // identity, length 1 — no motion any frame
    trail.record(still);
  }

  std::vector<float> gts;
  deep_tween(trail,
             [&](const Quaternion &, float gt) { gts.push_back(gt); });

  HS_EXPECT_EQ(gts.size(), static_cast<size_t>(1));
  HS_EXPECT_NEAR(gts.back(), 1.0f, 1e-6f);
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

// ============================================================================
// MeshMorph (nearest-vertex SLERP + crossfade invariant)
// ----------------------------------------------------------------------------
// MeshMorph maps each destination vertex to its nearest source vertex, SLERPs
// between them by an eased alpha, and crossfades the two shaded meshes with
// op_A = 1 - alpha. Using an octahedron as BOTH source and dest makes the
// nearest-vertex map the identity (vertices are 90 deg apart, the symmetry-twist
// is tiny), so the SLERP output must equal the destination at every alpha — a
// decisive check on both the correspondence and the interpolation. We also
// capture the two render opacities to assert the op_A + alpha == 1 crossfade.
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
  MeshOps::compile(poly, src, arena);
  MeshOps::compile(poly, dst, arena);

  const int duration = 8;
  float opA = -1.0f, alpha = -1.0f;
  bool mesh_b_tracks_dest = true;

  // Effect-owned lvalue callables (MeshMorph stores non-owning FunctionRefs and
  // rejects temporaries — see the borrow_guard static_asserts above).
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

  HS_EXPECT_TRUE(saw_both_layers);  // the crossfade was actually exercised
  HS_EXPECT_TRUE(mesh_b_tracks_dest); // identity correspondence + SLERP identity
  HS_EXPECT_NEAR(alpha, 1.0f, 1e-3f); // fully faded to the incoming mesh
  HS_EXPECT_TRUE(morph.done());
}

/**
 * @brief Runs every animation/easing test case in this module.
 * @return The module's failure count.
 */
inline int run_animation_tests() {
  auto scope = hs_test::begin_module("animation");

  hs::random().seed(1337); // OrientationTrail/animations may touch global RNG

  // Module-scoped fake-canvas fixture (see "Stand-in Canvas reference" above):
  // one Effect + Canvas for the whole module, published via fake_canvas_ptr().
  // Destroying it on return releases the single-live-Effect guard so a following
  // module in an all-modules run starts clean, while preserving the single-Canvas
  // (no second queue_frame) coupling within this module.
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

  test_lerp_drives_subject_to_target();
  test_lerp_midpoint();

  test_orientation_trail_index_zero_is_oldest();
  test_orientation_trail_expire_drops_oldest();
  test_orientation_trail_clear();

  test_easing_in_out_endpoints();
  test_easing_output_finite_over_range();
  test_easing_elastic_anchors();

  test_rotation_substeps_shared_and_tight();
  test_rotation_accumulates_subthreshold_deltas();
  test_timeline_shared_orientation_composes_motion_blur();
  test_timeline_sequences_events_by_start_frame();
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
  test_particle_system_expires_after_life_and_trail_drain();
  test_particle_system_attractor_kills_within_radius();

  test_sprite_fade_in_plateau_fade_out_envelope();
  test_sprite_overlapping_fades_stay_continuous();
  test_sprite_paused_holds_frame();

  test_deep_tween_global_t_spans_unit_interval();
  test_deep_tween_collapsed_newest_frame_reaches_one();
  test_deep_tween_all_collapsed_reaches_one();
  test_tween_vectortrail_single_sample_reaches_one();

  test_meshmorph_identity_self_map_and_crossfade();

  const int result = hs_test::end_module(scope);
  // Unpublish before fake_cv/fake_fx destruct below so the dangling pointer is
  // never observable; the FakeEffect dtor then clears the single-live guard.
  fake_canvas_ptr() = nullptr;
  return result;
}

} // namespace animation_tests
} // namespace hs_test

