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
 * pulled. Every Timeline<W,CAPACITY> shares one global event array (plus the
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
#include "core/animation.h"
#include "core/canvas.h"
#include "core/easing.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace animation_tests {

// ============================================================================
// Stand-in Canvas reference
// ----------------------------------------------------------------------------
// The animations under test (Transition/Mutation/Lerp/Driver) take a Canvas& in
// step() but never dereference it. Previously this handed step() a reference
// into raw aligned storage with no Canvas ever constructed — undefined behavior,
// since binding a reference to (and later forming a glvalue for) an object that
// doesn't exist is UB even if it's never read through.
//
// Instead we construct ONE genuine Canvas over a tiny static Effect and reuse
// it. A fresh effect's buffer_free() is true, so the ctor doesn't spin; we never
// construct a second Canvas on this effect, so there's no pending-frame spin and
// the single dtor (queue_frame) only runs harmlessly at process exit. The
// animations still never touch it — the reference is simply valid now.
//
// COUPLING (load-bearing): nothing a test does with this canvas may invoke
// queue_frame() — i.e. no test may construct (and destruct) a *second* Canvas on
// this shared FakeEffect during the run. queue_frame() advances next_, so a
// second Canvas's ctor would then see buffer_free()==false and spin on a display
// ISR that never runs in the host harness (a hang). The only queue_frame() that
// may fire is the single static-dtor one at process exit, after all tests pass.
// If a future animation test needs to actually queue a frame, give it its own
// local Canvas/Effect rather than reusing fake_canvas().
// ============================================================================

struct FakeEffect : public Effect {
  FakeEffect() : Effect(8, 8) {}
  void draw_frame() override {}
  bool show_bg() const override { return false; }
};

inline Canvas &fake_canvas() {
  static FakeEffect fx;
  static Canvas canvas(fx);
  return canvas;
}

// ============================================================================
// Path::get_point
// ============================================================================

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

inline void test_path_endpoints_and_clamp() {
  Path<32> p;
  // Plot a straight ramp along X: f(s) = (s, 0, 0), domain 1, 4 samples,
  // linear easing. append_segment samples t/samples in [0,1] -> x in [0,1].
  p.append_segment([](float s) { return Vector(s, 0.0f, 0.0f); }, 1.0f, 4.0f,
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

inline void test_path_collapse_keeps_last() {
  Path<32> p;
  p.append_segment([](float s) { return Vector(s, 0.0f, 0.0f); }, 1.0f, 4.0f,
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

inline void test_transition_duration_zero_no_divide_by_zero() {
  // AnimationBase coerces duration 0 -> 1, so the t/duration ratio is finite.
  float v = 0.0f;
  Animation::Transition tr(v, 42.0f, 0, ease_mid);
  tr.step(fake_canvas());
  HS_EXPECT_TRUE(std::isfinite(v));
  HS_EXPECT_NEAR(v, 42.0f, 1e-3f); // one step of a 1-frame transition completes
  HS_EXPECT_TRUE(tr.done());
}

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

// ============================================================================
// Mutation
// ============================================================================

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
struct Lerpable {
  float value = 0.0f;
  void lerp(const Lerpable &a, const Lerpable &b, float t) {
    value = a.value + (b.value - a.value) * t;
  }
};
} // namespace

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
// DOC CONTRADICTION (animation.h ~230/238 say "0 is newest", ~231 comment says
// "JS parity: 0 is oldest"). The CODE wins: OrientationTrail::record() calls
// snapshots.push_back(); StaticCircularBuffer::operator[](0) returns
// buffer[head], i.e. the FIRST element pushed. So index 0 is the OLDEST
// snapshot and index length()-1 is the newest. We assert the actual behavior.
// ============================================================================

inline void test_orientation_trail_index_zero_is_oldest() {
  Animation::OrientationTrail<Orientation<32, 4>, 8> trail;

  Orientation<32, 4> a, b, c;
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

inline void test_orientation_trail_expire_drops_oldest() {
  Animation::OrientationTrail<Orientation<32, 4>, 8> trail;
  Orientation<32, 4> a, b;
  a.set(make_rotation(Vector(0, 1, 0), 0.1f));
  b.set(make_rotation(Vector(0, 1, 0), 0.2f));
  trail.record(a);
  trail.record(b);
  trail.expire(); // removes the oldest (a)
  HS_EXPECT_EQ(trail.length(), static_cast<size_t>(1));
  // The survivor is b, now at index 0.
  HS_EXPECT_NEAR(std::abs(dot(trail.get(0).get(), b.get())), 1.0f, 1e-4f);
}

inline void test_orientation_trail_clear() {
  Animation::OrientationTrail<Orientation<32, 4>, 8> trail;
  Orientation<32, 4> a;
  trail.record(a);
  trail.record(a);
  trail.clear();
  HS_EXPECT_EQ(trail.length(), static_cast<size_t>(0));
}

// ============================================================================
// easing.h sanity
// ============================================================================

inline void test_easing_in_out_endpoints() {
  // In/out variants should anchor at 0 and 1.
  struct {
    const char *name;
    EasingFn fn;
  } fns[] = {
      {"ease_in_out_bicubic", ease_in_out_bicubic},
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

inline void test_easing_output_finite_over_range() {
  EasingFn fns[] = {ease_in_out_bicubic, ease_in_out_sin, ease_in_sin,
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

inline void test_easing_elastic_anchors() {
  // ease_out_elastic is not an in/out pair but is defined to hit 0 and 1.
  HS_EXPECT_NEAR(ease_out_elastic(0.0f), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(ease_out_elastic(1.0f), 1.0f, 1e-3f);
}

// ============================================================================
// Borrow-contract guards (compile-time)
// ----------------------------------------------------------------------------
// Motion and MeshMorph store a non-owning borrow (Fn capturing &path_obj /
// FunctionRef) of effect-owned state that must outlive the Timeline. Their
// constructors have a `= delete` rvalue overload so binding to a temporary is
// a compile error instead of a silent dangle. These static_asserts lock that
// contract: an effect-owned lvalue must be accepted, a temporary must be
// rejected. They live at namespace scope (checked at compile time; no runner
// call needed). If a refactor turns the SFINAE guard into a no-op, this fails
// to compile. See P1 #5 in docs/CODE_REVIEW.md.
namespace borrow_guard {
using Ori = Orientation<288, 16>;
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
struct Lerpable {
  float v = 0.0f;
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

// Motion and Rotation size their orientation trail through the shared
// Animation::rotation_substeps() helper so the same sweep yields the same subdivision.
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

// Two animations sharing one Orientation must COMPOSE their sub-frame
// motion-blur history within a frame, not clobber it. The earlier bug
// collapsed the Orientation once per animation, so the second animation's
// collapse discarded the first's freshly-built sub-frame trail. The decisive
// signature is the OLDEST sub-frame (index 0): with composition it still
// reflects the pre-frame orientation (identity here); with the bug it has
// already advanced by the first rotation's full sweep. (Unlike most tests in
// this file we DO touch the global Timeline; Rotation::step never dereferences
// the canvas, and we reset the global cursor state around the test.)
inline void test_timeline_shared_orientation_composes_motion_blur() {
  using Ori = Orientation<288, 16>;
  global_timeline_num_events = 0;
  global_timeline_t = 0;

  Ori o; // identity, single frame
  Timeline<288> tl;
  // Two world-space rotations about the SAME axis, each a quarter turn,
  // completing in one frame so the single step sweeps the full angle.
  tl.add(0, Animation::Rotation<288, 16>(o, Z_AXIS, PI_F / 2, 1, ease_mid));
  tl.add(0, Animation::Rotation<288, 16>(o, Z_AXIS, PI_F / 2, 1, ease_mid));
  tl.step(fake_canvas());

  // A motion-blur trail, not a single collapsed frame.
  HS_EXPECT_GE(o.length(), 2);
  // Oldest sub-frame is the pre-frame orientation (identity): +X stays +X.
  // (Bug: index 0 already rotated by the first quarter turn -> +Y.)
  Vector oldest = o.orient(X_AXIS, 0);
  HS_EXPECT_NEAR(oldest.x, 1.0f, 1e-3f);
  HS_EXPECT_NEAR(oldest.y, 0.0f, 1e-3f);
  // Newest sub-frame reflects BOTH rotations (a half turn): +X -> -X.
  Vector newest = o.orient(X_AXIS, o.length() - 1);
  HS_EXPECT_NEAR(newest.x, -1.0f, 1e-3f);

  global_timeline_num_events = 0;
  global_timeline_t = 0;
}

// Timeline schedules events by start frame: an event added with in_frames > 0
// stays dormant until t reaches its start, then steps; a one-shot animation is
// removed once done, so a later event can run after an earlier one finishes.
// (Uses the global Timeline; each Timeline is scoped so the live-guard balances,
// and Transition::step never dereferences the canvas.)
inline void test_timeline_sequences_events_by_start_frame() {
  Timeline<16> tl; // ctor clears the global cursors; dtor releases the live guard
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

// A repeating animation rewinds at the end of each cycle instead of being
// removed, so its bound value replays the curve. Mutation writes
// f(easing(t/duration)) each step, so after completion + rewind the next step
// drops back to the mid-cycle value (a non-rewinding timer would clamp at 1).
inline void test_timeline_repeating_animation_rewinds_each_cycle() {
  Timeline<16> tl;
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

// When a non-repeating event completes and is removed, step() compacts the
// array: later survivors are relocated (move_into) into the freed slots and must
// keep stepping correctly from their new positions. Decisive check: the
// originally-LAST event (relocated furthest) still reaches its own target.
inline void test_timeline_compaction_preserves_later_events() {
  Timeline<16> tl;
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

// .then() fires on completion; a callback may schedule a follow-up on the same
// Timeline mid-step. step() appends such events past the active snapshot, then
// gap-fills them into the freed slots — the add-during-callback path. The
// follow-up is added this frame but only runs on the next.
inline void test_timeline_then_chains_follow_up_event() {
  Timeline<16> tl;
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

// clear() destroys all events and resets the frame cursor, leaving the timeline
// reusable from t=0 — the in-place reset the singleton offers in lieu of
// reassignment.
inline void test_timeline_clear_resets_state() {
  Timeline<16> tl;
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

// add() is bounded by MAX_EVENTS (the shared global array size): once full, a
// further add is rejected (logged) rather than overrunning the array.
inline void test_timeline_full_guard_rejects_overflow() {
  Timeline<16> tl;
  float sink = 0.0f;
  // Targets share one float — this test asserts the count guard, not values, and
  // never steps. Fill to capacity.
  for (int i = 0; i < Timeline<16>::MAX_EVENTS; ++i)
    tl.add(0, Animation::Transition(sink, 1.0f, 10, ease_mid));
  HS_EXPECT_EQ(global_timeline_num_events, Timeline<16>::MAX_EVENTS);

  tl.add(0, Animation::Transition(sink, 1.0f, 10, ease_mid)); // one past full
  HS_EXPECT_EQ(global_timeline_num_events,
               Timeline<16>::MAX_EVENTS); // rejected, count unchanged
}

// Orientation::upsample SLERP-interpolates the recorded sub-frames up to a target
// count (preserving the endpoints), and collapse() discards all but the newest —
// the two primitives multi-animation motion blur is built on.
inline void test_orientation_upsample_then_collapse() {
  Orientation<32, 8> o; // identity, 1 frame
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

// A repeating Motion integrates its Orientation by dead reckoning. Without
// re-anchoring, the float error in that per-frame quaternion chain accumulates
// without bound and warps the traced curve into a "squiggly" mess after many
// cycles (the Comets path-corruption bug). The fix snaps the orientation back
// to its captured start at every cycle boundary, so each cycle integrates from
// an identical base and the shape stops drifting.
//
// The decisive, precession-immune signature is the set of rotation-INVARIANT
// internal angles between heads sampled at fixed phases within a cycle: a rigid
// drift (holonomy) leaves them unchanged, so any growth is genuine warp. We
// compare a late cycle against the ideal Lissajous internal angles; with the
// bug these diverge by > 1 rad, with the fix they stay at the first cycle's
// tiny discretization residual.
inline void test_motion_repeating_does_not_drift() {
  using Ori = Orientation<288, 16>;
  global_timeline_num_events = 0;
  global_timeline_t = 0;

  const int duration = 40;
  ProceduralPath path;
  // A representative Comets Lissajous figure; lissajous(.,.,.,0) == +Y, so the
  // identity-start orientation places the head exactly on the path at phase 0.
  path.f = [](float t) { return lissajous(1.06f, 1.06f, 0.0f, t * 5.909f); };

  Ori o; // identity, single frame; orient(+Y) starts on the path
  const Vector node_v = Y_AXIS;

  Timeline<288> tl;
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
  // ideal within the first-cycle discretization residual (< 0.1 rad). The
  // pre-fix integrator fails this by more than an order of magnitude.
  HS_EXPECT_NEAR(angle_between(Pa, Pb), ideal_ab, 0.1f);
  HS_EXPECT_NEAR(angle_between(Pa, Pc), ideal_ac, 0.1f);
  HS_EXPECT_NEAR(angle_between(Pb, Pc), ideal_bc, 0.1f);

  global_timeline_num_events = 0;
  global_timeline_t = 0;
}

inline int run_animation_tests() {
  auto scope = hs_test::begin_module("animation");

  hs::random().seed(1337); // OrientationTrail/animations may touch global RNG

  test_path_empty_returns_origin();
  test_path_endpoints_and_clamp();
  test_path_collapse_keeps_last();

  test_transition_reaches_target_linear();
  test_transition_monotonic_and_starts_from_current();
  test_transition_duration_zero_no_divide_by_zero();
  test_transition_quantized_floors_result();

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
  test_timeline_shared_orientation_composes_motion_blur();
  test_timeline_sequences_events_by_start_frame();
  test_timeline_repeating_animation_rewinds_each_cycle();
  test_timeline_compaction_preserves_later_events();
  test_timeline_then_chains_follow_up_event();
  test_timeline_clear_resets_state();
  test_timeline_full_guard_rejects_overflow();
  test_orientation_upsample_then_collapse();
  test_motion_repeating_does_not_drift();

  return hs_test::end_module(scope);
}

} // namespace animation_tests
} // namespace hs_test

