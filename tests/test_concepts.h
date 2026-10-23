/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Contracts for callable wrappers, Tweenable history, and DissolveMask partitions.
 */
#pragma once

#include <type_traits>
#include <utility>

#include "core/platform/platform.h" // Fn = hs::inplace_function (defined before concepts.h)
#include "core/engine/concepts.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace concepts_tests {

/**
 * @brief A functor with distinguishable const and non-const call operators.
 * @details Shows which FunctionRef constructor bound the callable: non-const
 *          returns x + 1, const (which also binds rvalues) returns x + 100.
 */
struct DualCall {
  /** @brief Non-const call: marks the non-const-lvalue ctor path. */
  int operator()(int x) { return x + 1; }
  /** @brief Const call: marks the const-lvalue (and rvalue) ctor path. */
  int operator()(int x) const { return x + 100; }
};

/** @brief Result type that cannot convert to int. */
struct IncompatibleResult {};

/** @brief Callable with a return type incompatible with int. */
struct ReturnsIncompatible {
  IncompatibleResult operator()(int) const { return {}; }
};

/** @brief Callable returning void. */
struct ReturnsVoid {
  void operator()(int) const {}
};

/** @brief Callable returning int. */
struct ReturnsInt {
  int operator()(int x) const { return x + 1; }
};

/**
 * @brief Free function exercising the plain-function-pointer ctor path.
 * @param x Input value.
 * @return x + 2.
 */
inline int plus_two(int x) { return x + 2; }

/**
 * @brief Borrows a FunctionRef by value and invokes it immediately.
 * @param x Argument forwarded to the borrowed callable.
 * @param f Call-scoped FunctionRef; a bound temporary lives for the full
 *          call expression, so this is the safe rvalue-borrow idiom.
 * @return The result of invoking f(x).
 */
inline int call_scoped(int x, FunctionRef<int(int)> f) { return f(x); }

/**
 * @brief Verifies FunctionRef selects the ctor matching the argument's value
 *        category.
 * @details A non-const lvalue binds the non-const ctor (x + 1); a const lvalue
 *          and an rvalue temporary both bind the const ctor (x + 100). The rvalue
 *          case is invoked within one full expression so the temporary is alive.
 */
inline void test_functionref_overload_resolution() {
  DualCall f;
  FunctionRef<int(int)> r_nonconst = f; // non-const lvalue -> non-const ctor
  HS_EXPECT_EQ(r_nonconst(0), 1);

  const DualCall cf{};
  FunctionRef<int(int)> r_const = cf; // const lvalue -> const ctor
  HS_EXPECT_EQ(r_const(0), 100);

  HS_EXPECT_EQ(call_scoped(0, DualCall{}), 100);
}

/**
 * @brief Verifies the function-pointer ctor, including the null-pointer guard.
 * @details A live function pointer invokes through the thunk; a null pointer
 *          yields an empty ref (operator bool is false).
 */
inline void test_functionref_function_pointer() {
  FunctionRef<int(int)> r = &plus_two;
  HS_EXPECT_TRUE((bool)r);
  HS_EXPECT_EQ(r(5), 7);

  int (*nullfp)(int) = nullptr;
  FunctionRef<int(int)> empty = nullfp;
  HS_EXPECT_FALSE((bool)empty);

  auto nonthrowing_pointer = +[](int value) noexcept { return value + 3; };
  FunctionRef<int(int)> stored = nonthrowing_pointer;
  nonthrowing_pointer = nullptr;
  HS_EXPECT_TRUE(static_cast<bool>(stored));
  HS_EXPECT_EQ(stored(5), 8);
  FunctionRef<int(int)> nonthrowing_empty = nonthrowing_pointer;
  HS_EXPECT_FALSE(static_cast<bool>(nonthrowing_empty));
  FunctionRef<int(int)> temporary =
      +[](int value) noexcept { return value + 4; };
  HS_EXPECT_EQ(temporary(5), 9);
}

/**
 * @brief Verifies empty-state construction and copy/move of FunctionRef.
 * @details Default and nullptr construction both read empty; a copy or move
 *          (a trivial two-pointer payload) keeps referring to the same callable.
 */
inline void test_functionref_empty_and_copy() {
  FunctionRef<int(int)> def; // default -> empty
  HS_EXPECT_FALSE((bool)def);
  FunctionRef<int(int)> null = nullptr; // nullptr ctor -> empty
  HS_EXPECT_FALSE((bool)null);

  DualCall f;
  FunctionRef<int(int)> r = f;
  FunctionRef<int(int)> copy = r;
  HS_EXPECT_TRUE((bool)copy);
  HS_EXPECT_EQ(copy(0), 1);
  FunctionRef<int(int)> moved = std::move(r);
  HS_EXPECT_EQ(moved(0), 1);
}

/** Counts every copy and move a by-value argument pays crossing the erasure. */
struct CountedArg {
  static inline int constructions = 0;
  CountedArg() = default;
  CountedArg(const CountedArg &) { ++constructions; }
  CountedArg(CountedArg &&) noexcept { ++constructions; }
};

/**
 * @brief Verifies a by-value argument crosses the erasure without an extra copy,
 *        and that reference and decayed-array parameters still bind.
 * @details An lvalue argument pays one copy into operator()'s parameter and
 *          one move into the callable's. A reference parameter must reach the
 *          callable as the caller's object.
 */
inline void test_functionref_argument_forwarding() {
  auto sink = [](CountedArg) {};
  FunctionRef<void(CountedArg)> by_value = sink;
  CountedArg arg;
  CountedArg::constructions = 0;
  by_value(arg);
  HS_EXPECT_EQ(CountedArg::constructions, 2);

  auto bump = [](int &n) { n += 3; };
  FunctionRef<void(int &)> by_ref = bump;
  int n = 1;
  by_ref(n);
  HS_EXPECT_EQ(n, 4);

  // An array parameter decays to a pointer in the function type.
  auto first = [](int *a) { return a[0]; };
  FunctionRef<int(int[4])> by_array = first;
  int values[4] = {9, 0, 0, 0};
  HS_EXPECT_EQ(by_array(values), 9);
}

/**
 * @brief Verifies StoredFunctionRef enforces the borrow-vs-store lifetime
 *        contract in the type system.
 * @details A StoredFunctionRef accepts an lvalue callable but rejects an
 *          rvalue temporary (it would dangle past the call); plain FunctionRef
 *          accepts the rvalue borrow.
 */
inline void test_stored_functionref_rvalue_rejection() {
  static_assert(
      std::is_constructible_v<StoredFunctionRef<int(int)>, DualCall &>,
      "StoredFunctionRef must accept an lvalue callable");
  static_assert(
      !std::is_constructible_v<StoredFunctionRef<int(int)>, DualCall &&>,
      "StoredFunctionRef must reject an rvalue temporary (dangling-store guard)");
  static_assert(std::is_constructible_v<FunctionRef<int(int)>, DualCall &&>,
                "FunctionRef must accept an rvalue (call-scoped borrow)");

  static_assert(
      std::is_constructible_v<StoredFunctionRef<int(int)>, int (*)(int)>);
  static_assert(
      std::is_constructible_v<StoredFunctionRef<int(int)>, std::nullptr_t>);
  StoredFunctionRef<int(int)> pointer = +[](int value) { return value + 2; };
  HS_EXPECT_EQ(pointer(3), 5);
  StoredFunctionRef<int(int)> nonthrowing =
      +[](int value) noexcept { return value + 3; };
  HS_EXPECT_EQ(nonthrowing(3), 6);
  static_assert(!std::is_constructible_v<StoredFunctionRef<int(int)>,
                                         decltype([](int x) { return x; })>);
  StoredFunctionRef<int(int)> empty = nullptr;
  HS_EXPECT_FALSE(static_cast<bool>(empty));

  DualCall f;
  StoredFunctionRef<int(int)> s = f;
  HS_EXPECT_EQ(s(0), 1);
}

/**
 * @brief Verifies callable wrappers enforce and implement return compatibility.
 */
inline void test_callable_return_constraints() {
  static_assert(
      !std::is_constructible_v<FunctionRef<int(int)>, ReturnsIncompatible &>);
  static_assert(!std::is_constructible_v<FunctionRef<int(int)>,
                                         const ReturnsIncompatible &>);
  static_assert(
      !std::is_constructible_v<FunctionRef<int(int)>, ReturnsIncompatible &&>);
  static_assert(!std::is_constructible_v<StoredFunctionRef<int(int)>,
                                         ReturnsIncompatible &>);
  static_assert(!std::is_constructible_v<FunctionRef<int(int)>, ReturnsVoid &>);
  static_assert(std::is_constructible_v<FunctionRef<void(int)>, ReturnsInt &>);

  static_assert(
      !std::is_constructible_v<Fn<int(int), 16>, ReturnsIncompatible>);
  static_assert(!std::is_constructible_v<Fn<int(int), 16>, ReturnsVoid>);
  static_assert(std::is_constructible_v<Fn<void(int), 16>, ReturnsInt>);
  static_assert(!std::is_assignable_v<Fn<int(int), 16> &, ReturnsIncompatible>);
  static_assert(!std::is_constructible_v<Fn<int(int), 24>, Fn<int(int), 16>>);
  static_assert(!std::is_assignable_v<Fn<int(int), 24> &, Fn<int(int), 16>>);

  int called_with = 0;
  auto returns_int = [&called_with](int value) {
    called_with = value;
    return value + 1;
  };
  FunctionRef<void(int)> borrowed = returns_int;
  borrowed(1);
  HS_EXPECT_EQ(called_with, 1);
  Fn<void(int), 16> owned = returns_int;
  owned(2);
  HS_EXPECT_EQ(called_with, 2);
  HS_EXPECT_TRUE((bool)borrowed);
  HS_EXPECT_TRUE((bool)owned);
}

/**
 * @brief Verifies Fn (hs::inplace_function) copy/move/empty value semantics.
 * @details Default and nullptr read empty; a copy leaves the source live; a
 *          move empties it; assignment onto a populated Fn overwrites.
 */
inline void test_fn_copy_move_empty() {
  // ArenaVector accepts Fn only while it stays trivially destructible.
  static_assert(std::is_trivially_destructible_v<Fn<int(int), 16>>);

  Fn<int(int), 16> def;
  HS_EXPECT_FALSE((bool)def);
  Fn<int(int), 16> null = nullptr;
  HS_EXPECT_FALSE((bool)null);

  int base = 10;
  Fn<int(int), 16> f = [base](int x) { return base + x; };
  HS_EXPECT_TRUE((bool)f);
  HS_EXPECT_EQ(f(5), 15);

  Fn<int(int), 16> g = f;
  HS_EXPECT_TRUE((bool)f);
  HS_EXPECT_EQ(g(5), 15);

  Fn<int(int), 16> h = std::move(f);
  HS_EXPECT_EQ(h(5), 15);
  HS_EXPECT_FALSE((bool)f);

  Fn<int(int), 16> a;
  a = [](int x) { return x * 2; };
  HS_EXPECT_TRUE((bool)a);
  HS_EXPECT_EQ(a(3), 6);
  a = nullptr;
  HS_EXPECT_FALSE((bool)a);

  Fn<int(int), 16> p = [](int x) { return x + 100; };
  Fn<int(int), 16> q = [mul = 3](int x) { return x * mul; };
  p = q;
  HS_EXPECT_EQ(p(4), 12);
  HS_EXPECT_EQ(q(4), 12); // copy leaves the source live

  Fn<int(int), 16> r = [](int x) { return x - 7; };
  r = std::move(q);
  HS_EXPECT_EQ(r(4), 12);
  HS_EXPECT_FALSE((bool)q); // move leaves the source empty
}

struct TimedRefRecorder {
  Canvas *canvas = nullptr;
  math::Vector point;
  math::PixelCoords projected{};
  Pixel color;
  float age = 0.0f;
  float alpha = 0.0f;
  float time = 0.0f;
  int route = 0;

  void record(Canvas &cv, const Pixel &c, float a, float opacity) {
    canvas = &cv;
    color = c;
    age = a;
    alpha = opacity;
  }
  void plot(Canvas &cv, float x, float y, const Pixel &c, float a,
            float opacity) {
    record(cv, c, a, opacity);
    projected = {x, y};
    route = 1;
  }
  void plot(Canvas &cv, int x, int y, const Pixel &c, float a, float opacity) {
    plot(cv, static_cast<float>(x), static_cast<float>(y), c, a, opacity);
  }
  void plot(Canvas &cv, const math::Vector &v, const Pixel &c, float a,
            float opacity) {
    record(cv, c, a, opacity);
    point = v;
    route = 2;
  }
  math::PixelCoords probe(Canvas &, const math::Vector &v, float t) {
    return {v.x + 4.0f * t, v.y - 2.0f * t};
  }
  float plot_at(Canvas &cv, const math::Vector &v, const Pixel &c, float a,
                float opacity, float t) {
    plot(cv, v, c, a, opacity);
    time = t;
    route = 3;
    return 17.0f + t;
  }
};

inline void expect_timed_ref_registers(const TimedRefRecorder &sink,
                                       const Canvas &canvas, const Pixel &color,
                                       float age, float alpha) {
  HS_EXPECT_EQ(sink.canvas, &canvas);
  HS_EXPECT_EQ(sink.color, color);
  HS_EXPECT_EQ(sink.age, age);
  HS_EXPECT_EQ(sink.alpha, alpha);
}

/** @brief Minimal timed providers retain identity and ordinary plot fallbacks. */
inline void test_pipeline_ref_timed_defaults() {
  StubEffect effect(8, 4);
  Canvas canvas(effect);
  TimedRefRecorder sink;
  PipelineRef ordinary(sink);
  PipelineRef timed(sink, Spacetime::Timed{});
  const math::Vector point(0.25f, -0.5f, 0.75f);
  const Pixel color(123, 456, 789);

  for (const PipelineRef *ref : {&ordinary, &timed}) {
    HS_EXPECT_TRUE(ref->world_transform_is_identity);
    HS_EXPECT_TRUE(ref->has_temporal_motion);
    HS_EXPECT_EQ(ref->sweep_arc(), 0.0f);
    HS_EXPECT_EQ(ref->copy_count(), 1);
    HS_EXPECT_FALSE(ref->rigid_world_stages());
    HS_EXPECT_FALSE(ref->dense_screen_safe());
    ref->select_copy(2);
  }
  const auto world = timed.world_at(canvas, point, 0.375f);
  HS_EXPECT_EQ(world.x, 0.25f);
  HS_EXPECT_EQ(world.y, -0.5f);
  HS_EXPECT_EQ(world.z, 0.75f);
  const auto projected = timed.probe(canvas, point, 0.375f);
  HS_EXPECT_EQ(projected.x, 1.75f);
  HS_EXPECT_EQ(projected.y, -1.25f);
  HS_EXPECT_EQ(sink.route, 0);

  HS_EXPECT_EQ(
      timed.plot_transformed_at(canvas, point, color, 6.5f, 0.25f, 0.375f),
      17.375f);
  HS_EXPECT_EQ(sink.route, 3);
  HS_EXPECT_EQ(sink.point.x, 0.25f);
  HS_EXPECT_EQ(sink.point.y, -0.5f);
  HS_EXPECT_EQ(sink.point.z, 0.75f);
  HS_EXPECT_EQ(sink.time, 0.375f);
  expect_timed_ref_registers(sink, canvas, color, 6.5f, 0.25f);
  const Pixel projected_color(987, 654, 321);
  HS_EXPECT_EQ(timed.plot_projected_at(canvas, {2.25f, 1.5f}, projected_color,
                                       2.5f, 0.5f, 0.875f),
               2.25f);
  HS_EXPECT_EQ(sink.route, 1);
  HS_EXPECT_EQ(sink.projected.x, 2.25f);
  HS_EXPECT_EQ(sink.projected.y, 1.5f);
  expect_timed_ref_registers(sink, canvas, projected_color, 2.5f, 0.5f);
}

struct TimedRefStages : TimedRefRecorder {
  static constexpr bool world_transform_is_identity = false;
  static constexpr bool has_temporal_motion = false;
  static constexpr bool rigid_world_stages = true;
  static constexpr int total_segment_margin = 1;
  static constexpr bool any_2d_history = false;
  static constexpr bool any_reads_outside_band = false;
  static constexpr bool column_taps_ready = true;
  int selected = -1;
  int column = -1;

  float sweep_arc() const { return 0.125f; }
  int copy_count() const { return 3; }
  void select_copy(int copy) { selected = copy; }
  math::Vector world_at(Canvas &, const math::Vector &v, float t) {
    return {v.z + t, v.x - t, v.y};
  }
  float plot_transformed_at(Canvas &cv, const math::Vector &v, const Pixel &c,
                            float a, float opacity, float t) {
    plot_at(cv, v, c, a, opacity, t);
    route = 4;
    return -9.0f + t;
  }
  float plot_projected_at(Canvas &cv, math::PixelCoords p, const Pixel &c,
                          float a, float opacity, float t) {
    plot(cv, p.x, p.y, c, a, opacity);
    time = t;
    route = 5;
    return p.y + t;
  }
  void plot_column_at(Canvas &cv, math::PixelCoords p, int col, const Pixel &c,
                      float a, float opacity, float t) {
    plot_projected_at(cv, p, c, a, opacity, t);
    column = col;
    route = 6;
  }
};

/** @brief Timed erasure preserves explicit stages, copies and shading registers. */
inline void test_pipeline_ref_timed_stages() {
  StubEffect effect(8, 4);
  Canvas canvas(effect);
  TimedRefStages sink;
  PipelineRef ref(sink, Spacetime::Timed{});
  PipelineRef copy(ref);
  const math::Vector point(0.25f, -0.5f, 0.75f);
  const Pixel color(321, 654, 987);
  HS_EXPECT_FALSE(copy.world_transform_is_identity);
  HS_EXPECT_FALSE(copy.has_temporal_motion);
  HS_EXPECT_TRUE(copy.rigid_world_stages());
  HS_EXPECT_TRUE(copy.dense_screen_safe());
  HS_EXPECT_EQ(copy.sweep_arc(), 0.125f);
  HS_EXPECT_EQ(copy.copy_count(), 3);
  copy.select_copy(2);
  HS_EXPECT_EQ(sink.selected, 2);
  const auto world = copy.world_at(canvas, point, 0.375f);
  HS_EXPECT_EQ(world.x, 1.125f);
  HS_EXPECT_EQ(world.y, -0.125f);
  HS_EXPECT_EQ(world.z, -0.5f);

  HS_EXPECT_EQ(copy.plot_at(canvas, point, color, 4.5f, 0.75f, 0.25f), 17.25f);
  HS_EXPECT_EQ(sink.route, 3);
  HS_EXPECT_EQ(sink.point.x, 0.25f);
  HS_EXPECT_EQ(sink.point.y, -0.5f);
  HS_EXPECT_EQ(sink.point.z, 0.75f);
  HS_EXPECT_EQ(sink.time, 0.25f);
  expect_timed_ref_registers(sink, canvas, color, 4.5f, 0.75f);
  const Pixel transformed_color(135, 246, 357);
  HS_EXPECT_EQ(copy.plot_transformed_at(canvas, point, transformed_color, 7.5f,
                                        0.5f, 0.375f),
               -8.625f);
  HS_EXPECT_EQ(sink.route, 4);
  HS_EXPECT_EQ(sink.point.x, 0.25f);
  HS_EXPECT_EQ(sink.point.y, -0.5f);
  HS_EXPECT_EQ(sink.point.z, 0.75f);
  HS_EXPECT_EQ(sink.time, 0.375f);
  expect_timed_ref_registers(sink, canvas, transformed_color, 7.5f, 0.5f);
  const Pixel projected_color(987, 654, 321);
  HS_EXPECT_EQ(copy.plot_projected_at(canvas, {2.25f, 1.5f}, projected_color,
                                      2.5f, 0.25f, 0.875f),
               2.375f);
  HS_EXPECT_EQ(sink.route, 5);
  HS_EXPECT_EQ(sink.projected.x, 2.25f);
  HS_EXPECT_EQ(sink.projected.y, 1.5f);
  HS_EXPECT_EQ(sink.time, 0.875f);
  expect_timed_ref_registers(sink, canvas, projected_color, 2.5f, 0.25f);
  copy.plot_column_at(canvas, {3.25f, 2.5f}, 3, color, 4.5f, 0.75f, 0.625f);
  HS_EXPECT_EQ(sink.route, 6);
  HS_EXPECT_EQ(sink.projected.x, 3.25f);
  HS_EXPECT_EQ(sink.projected.y, 2.5f);
  HS_EXPECT_EQ(sink.column, 3);
  expect_timed_ref_registers(sink, canvas, color, 4.5f, 0.75f);
  HS_EXPECT_EQ(sink.time, 0.625f);
}

/**
 * @brief A minimal model of Tweenable: a container of frames, each with its own
 *        sub-frame quaternion history.
 * @details Mirrors OrientationTrail's shape — the two levels deep_tween walks.
 */
struct TweenableModel {
  /** @brief A frame: a fixed-capacity quaternion sub-history. */
  struct Frame {
    static constexpr int CAPACITY = 4; /**< Max sub-frames. */
    /** @brief Sub-frame count, signed as Orientation's is. */
    int length() const { return 0; }
    /** @brief Sub-frame accessor. */
    const math::Quaternion &get(int) const { return q; }
    math::Quaternion q; /**< Storage backing get(). */
  };
  /** @brief Frame count, consumed as a size_t. */
  size_t length() const { return 0; }
  /** @brief Frame accessor. */
  const Frame &get(size_t) const { return frame; }
  Frame frame; /**< Storage backing get(). */
};

/**
 * @brief A flat length()/get() container whose element is a scalar.
 * @details Stands in for a bare Orientation, whose get() yields a Quaternion:
 *          one level of history, nothing for deep_tween to flatten.
 */
struct FlatTweenableModel {
  /** @brief Element count, consumed as a size_t. */
  size_t length() const { return 0; }
  /** @brief Element accessor yielding a scalar. */
  int get(size_t) const { return 0; }
};

struct SignedTrailModel : TweenableModel {
  int length() const { return 0; }
};

struct UnsignedFrameModel : TweenableModel {
  struct Frame : TweenableModel::Frame {
    size_t length() const { return 0; }
  };
  const Frame &get(size_t) const { return frame; }
  Frame frame;
};

struct MissingCapacityModel : TweenableModel {
  struct Frame {
    int length() const { return 0; }
    const math::Quaternion &get(int) const { return q; }
    math::Quaternion q;
  };
  const Frame &get(size_t) const { return frame; }
  Frame frame;
};

/**
 * @brief Pins the Tweenable concept: only a two-level frame container satisfies
 *        it — a flat container and a scalar do not.
 */
inline void test_tweenable_concept() {
  static_assert(!Tweenable<SignedTrailModel>);
  static_assert(!Tweenable<UnsignedFrameModel>);
  static_assert(!Tweenable<MissingCapacityModel>);
  static_assert(Tweenable<TweenableModel>,
                "a container of frames must satisfy Tweenable");
  static_assert(!Tweenable<FlatTweenableModel>,
                "a flat length()/get() container must not satisfy Tweenable");
  static_assert(!Tweenable<int>, "a scalar must not satisfy Tweenable");
}

/**
 * @brief Pins DissolveMask's partition contract: two masks with the same
 *        threshold/salt and opposite invert own every key exactly once, the
 *        endpoints are total, the owned fraction tracks the threshold, and the
 *        hash's identity is frozen against a golden.
 */
inline void test_dissolve_mask_partition() {
  constexpr int KA = 288, KB = 144;
  const uint32_t salt = 0x1234ABCDu;

  DissolveMask incoming{65536u / 4u, salt, false};
  DissolveMask outgoing{65536u / 4u, salt, true};
  int owned_in = 0;
  for (int b = 0; b < KB; ++b) {
    for (int a = 0; a < KA; ++a) {
      HS_EXPECT_TRUE(incoming.owns(a, b) != outgoing.owns(a, b));
      owned_in += incoming.owns(a, b) ? 1 : 0;
    }
  }
  // 25% coverage within 2% absolute (hash uniformity, not exactness).
  float frac = static_cast<float>(owned_in) / (KA * KB);
  HS_EXPECT_TRUE(frac > 0.23f && frac < 0.27f);

  DissolveMask none{0u, salt, false};
  DissolveMask all{65536u, salt, false};
  for (int b = 0; b < KB; b += 7) {
    for (int a = 0; a < KA; a += 7) {
      HS_EXPECT_FALSE(none.owns(a, b));
      HS_EXPECT_TRUE(all.owns(a, b));
    }
  }

  // Frozen golden over one key row.
  uint32_t row = 0;
  for (int a = 0; a < 32; ++a)
    row |= (incoming.owns(a, 7) ? 1u : 0u) << a;
  HS_EXPECT_EQ(row, 0x20000222u);
}

/**
 * @brief Runs every concepts test under a "concepts" module scope.
 * @return Failure count reported by end_module for the "concepts" module.
 */
inline int run_concepts_tests() {
  hs_test::ModuleFixture fixture("concepts");

  test_functionref_overload_resolution();
  test_functionref_function_pointer();
  test_functionref_empty_and_copy();
  test_functionref_argument_forwarding();
  test_stored_functionref_rvalue_rejection();
  test_callable_return_constraints();
  test_fn_copy_move_empty();
  test_pipeline_ref_timed_defaults();
  test_pipeline_ref_timed_stages();
  test_tweenable_concept();
  test_dissolve_mask_partition();

  return fixture.result();
}

} // namespace concepts_tests
} // namespace hs_test
