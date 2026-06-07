/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/animation.h and core/easing.h.
 *
 * Scope: the PURE, non-render parts of the animation system. These tests
 * deliberately avoid constructing a live Canvas (its ctor busy-waits on the
 * render buffer and its dtor queues a frame). The animations exercised here
 * (Transition, Mutation, Lerp, Driver) take a `Canvas&` in step() but never
 * dereference it — AnimationBase::step only increments the frame counter and
 * the derived steps tested here only touch their bound float/subject. We pass
 * a reference into raw aligned storage (never read/written through it) so no
 * Canvas member is ever accessed. Timeline is intentionally NOT exercised
 * (inline-static shared global state; pulls the render stack).
 *
 * Coverage:
 *   - Path::get_point: empty guard, endpoints, end-of-path clamp, midpoint
 *   - AnimationBase duration==0 -> 1 coercion (no divide-by-zero)
 *   - Transition: start->end stepping, easing, quantize, done()
 *   - Mutation: applies f(easing(t))
 *   - Driver: per-frame increment + wrap
 *   - Lerp: type-erased lerp(start,target,t) driven by easing
 *   - OrientationTrail: record/get ordering (index 0 is OLDEST — see note)
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
} // namespace borrow_guard

// ============================================================================
// Runner
// ============================================================================

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

  return hs_test::end_module(scope);
}

} // namespace animation_tests
} // namespace hs_test

