/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/engine/concepts.h — the type-erased callable wrappers.
 */
#pragma once

#include <type_traits>
#include <utility>

#include "core/engine/platform.h" // Fn = hs::inplace_function (defined before concepts.h)
#include "core/engine/concepts.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace concepts_tests {

/**
 * @brief A functor with distinguishable const and non-const call operators.
 * @details Used to observe WHICH FunctionRef constructor bound the callable: the
 *          non-const-lvalue ctor invokes through a non-const pointer (returns
 *          x + 1), while the const-lvalue ctor — which also binds rvalues —
 *          invokes through a const pointer (returns x + 100).
 */
struct DualCall {
  /** @brief Non-const call: marks the non-const-lvalue ctor path. */
  int operator()(int x) { return x + 1; }
  /** @brief Const call: marks the const-lvalue (and rvalue) ctor path. */
  int operator()(int x) const { return x + 100; }
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
 *          yields an EMPTY ref (operator bool is false) rather than installing a
 *          non-null thunk that would dereference null on the first call.
 */
inline void test_functionref_function_pointer() {
  FunctionRef<int(int)> r = &plus_two;
  HS_EXPECT_TRUE((bool)r);
  HS_EXPECT_EQ(r(5), 7);

  int (*nullfp)(int) = nullptr;
  FunctionRef<int(int)> empty = nullfp;
  HS_EXPECT_FALSE((bool)empty);
}

/**
 * @brief Verifies empty-state construction and copy/move of FunctionRef.
 * @details Default and nullptr construction both read empty; a copy or move
 *          (a trivial two-pointer payload) keeps referring to the same callable.
 */
inline void test_functionref_empty_and_copy() {
  FunctionRef<int(int)> def;            // default -> empty
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

/**
 * @brief Verifies StoredFunctionRef enforces the borrow-vs-store lifetime
 *        contract in the type system.
 * @details A StoredFunctionRef must accept an lvalue callable (it outlives the
 *          call) but reject an rvalue temporary (binding one would dangle past
 *          the call); plain FunctionRef, by contrast, deliberately accepts the
 *          rvalue borrow. The constructibility traits double as a compile-time
 *          tripwire and a runtime-visible regression check.
 */
inline void test_stored_functionref_rvalue_rejection() {
  static_assert(std::is_constructible_v<StoredFunctionRef<int(int)>, DualCall &>,
                "StoredFunctionRef must accept an lvalue callable");
  static_assert(
      !std::is_constructible_v<StoredFunctionRef<int(int)>, DualCall &&>,
      "StoredFunctionRef must reject an rvalue temporary (dangling-store guard)");
  static_assert(std::is_constructible_v<FunctionRef<int(int)>, DualCall &&>,
                "FunctionRef must accept an rvalue (call-scoped borrow)");

  HS_EXPECT_TRUE(
      (std::is_constructible_v<StoredFunctionRef<int(int)>, DualCall &>));
  HS_EXPECT_FALSE(
      (std::is_constructible_v<StoredFunctionRef<int(int)>, DualCall &&>));
  HS_EXPECT_TRUE((std::is_constructible_v<FunctionRef<int(int)>, DualCall &&>));

  DualCall f;
  StoredFunctionRef<int(int)> s = f;
  HS_EXPECT_EQ(s(0), 1);
}

/**
 * @brief Verifies Fn (hs::inplace_function) copy/move/empty value semantics.
 * @details Default and nullptr construction read empty; a copy leaves the source
 *          live; a move leaves the source empty; assignment from a callable then
 *          back to nullptr toggles the empty state. The empty-state CALL trap is
 *          covered separately by the death harness (case_empty_fn_call).
 */
inline void test_fn_copy_move_empty() {
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
}

/**
 * @brief A minimal model of Tweenable: a length()/get() container.
 * @details length() returns a size_t count and get(index) yields an element —
 *          the shape Lerp/Transition require of their subject.
 */
struct TweenableModel {
  /** @brief Element count, consumed as a size_t. */
  size_t length() const { return 0; }
  /** @brief Element accessor; the concept leaves the return type deduced. */
  int get(size_t) const { return 0; }
};

/**
 * @brief Pins the Tweenable concept: a length()/get() container satisfies it, a
 *        scalar does not.
 * @details The static_asserts are the compile-time tripwire; the runtime checks
 *          mirror them so a regression is also reported by the harness.
 */
inline void test_tweenable_concept() {
  static_assert(Tweenable<TweenableModel>,
                "a length()/get() container must satisfy Tweenable");
  static_assert(!Tweenable<int>, "a scalar must not satisfy Tweenable");

  HS_EXPECT_TRUE(Tweenable<TweenableModel>);
  HS_EXPECT_FALSE(Tweenable<int>);
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
  test_stored_functionref_rvalue_rejection();
  test_fn_copy_move_empty();
  test_tweenable_concept();

  return fixture.result();
}

} // namespace concepts_tests
} // namespace hs_test
