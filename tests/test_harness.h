/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Shared test harness — single global Stats counter, HS_EXPECT_* macros,
 * and per-module begin/end scope helpers. Per-module reporters compute their
 * own pass/fail counts as deltas from the global counter so multiple test
 * suites can coexist in one driver.
 */
#pragma once

#include <cmath>
#include <cstdio>
#include <type_traits>

namespace hs_test {

/**
 * @brief Process-wide pass/fail tally shared by every test suite.
 * @details Single-threaded by contract: the counters are plain ints, so every
 * HS_EXPECT_* must be evaluated on the test's main thread. This keeps the hot
 * assertion path free of atomic RMWs across the whole suite. A test that spawns
 * a helper thread must capture its cross-thread observation into a std::atomic
 * and assert on the main thread after join() — never calling HS_EXPECT_* from
 * the helper.
 */
struct Stats {
  int passed = 0; /**< Count of comparisons that succeeded. */
  int failed = 0; /**< Count of comparisons that failed. */
};

/**
 * @brief Accessor for the single process-wide pass/fail counter.
 * @return Reference to the function-local static Stats shared by all suites.
 */
inline Stats &stats() {
  static Stats s;
  return s;
}

/**
 * @brief Tests whether two floats are finite and close enough to be equal.
 * @param a First value.
 * @param b Second value.
 * @param tol Maximum allowed absolute difference, in the same units as a and b.
 * @return True iff both values are finite and within tol of each other.
 */
inline bool approx(float a, float b, float tol) {
  return std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= tol;
}

/**
 * @brief Tests whether two doubles are finite and close enough to be equal.
 * @param a First value.
 * @param b Second value.
 * @param tol Maximum allowed absolute difference, in the same units as a and b.
 * @return True iff both values are finite and within tol of each other.
 * @details Double-domain overload. HS_EXPECT_NEAR captures its operands as
 * double, so the comparison must stay in double — downcasting to float would
 * round both sides to ~1e-7 resolution and silently mask a sub-float-epsilon
 * regression in a double-valued expression. Float call sites bind the float
 * overload above by exact match, so this is additive.
 */
inline bool approx(double a, double b, double tol) {
  return std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= tol;
}

namespace detail {
/**
 * @brief Structural detection for the engine's small value types so
 * print_operand can show components without test_harness.h depending on the
 * engine headers. Pixel/Pixel16 expose r/g/b; Vector exposes x/y/z.
 */
template <class T, class = void> struct has_rgb : std::false_type {};
template <class T>
struct has_rgb<T, std::void_t<decltype(std::declval<const T &>().r),
                              decltype(std::declval<const T &>().g),
                              decltype(std::declval<const T &>().b)>>
    : std::true_type {};

template <class T, class = void> struct has_xyz : std::false_type {};
template <class T>
struct has_xyz<T, std::void_t<decltype(std::declval<const T &>().x),
                              decltype(std::declval<const T &>().y),
                              decltype(std::declval<const T &>().z)>>
    : std::true_type {};
} // namespace detail

/**
 * @brief Prints a single comparison operand in a type-appropriate format.
 * @tparam T Operand type; bool, floating-point, enum, integral, pointer, and the
 * engine's r/g/b and x/y/z value types are formatted specially.
 * @param v The operand value to print.
 * @details Lets a failing HS_EXPECT_* line show the actual values, not just the
 * stringified expr. Pixel/Pixel16 (r,g,b) and Vector (x,y,z) are detected
 * structurally — so an HS_EXPECT_EQ/CMP on those prints components rather than
 * "?" — and their fields recurse through this same printer. Any other
 * non-arithmetic operand still falls back to "?".
 */
template <class T> inline void print_operand(const T &v) {
  if constexpr (std::is_same_v<T, bool>) {
    std::printf("%s", v ? "true" : "false");
  } else if constexpr (std::is_floating_point_v<T>) {
    std::printf("%g", static_cast<double>(v));
  } else if constexpr (std::is_enum_v<T>) {
    std::printf("%lld",
                static_cast<long long>(
                    static_cast<std::underlying_type_t<T>>(v)));
  } else if constexpr (std::is_integral_v<T>) {
    if constexpr (std::is_signed_v<T>)
      std::printf("%lld", static_cast<long long>(v));
    else
      std::printf("%llu", static_cast<unsigned long long>(v));
  } else if constexpr (std::is_pointer_v<T>) {
    std::printf("%p", static_cast<const void *>(v));
  } else if constexpr (detail::has_rgb<T>::value) {
    std::printf("(r=");
    print_operand(v.r);
    std::printf(" g=");
    print_operand(v.g);
    std::printf(" b=");
    print_operand(v.b);
    std::printf(")");
  } else if constexpr (detail::has_xyz<T>::value) {
    std::printf("(x=");
    print_operand(v.x);
    std::printf(" y=");
    print_operand(v.y);
    std::printf(" z=");
    print_operand(v.z);
    std::printf(")");
  } else {
    std::printf("?");
  }
}

/**
 * @brief Records the result of a comparison and, on failure, prints operands.
 * @tparam A Type of the first operand.
 * @tparam B Type of the second operand.
 * @param ok Whether the comparison succeeded.
 * @param a First operand, printed on failure.
 * @param b Second operand, printed on failure.
 * @param expr Stringified comparison expression for the failure message.
 * @param file Source file name of the call site.
 * @param line Source line number of the call site.
 */
template <class A, class B>
inline void report_cmp(bool ok, const A &a, const B &b, const char *expr,
                       const char *file, int line) {
  if (ok) {
    ++stats().passed;
    return;
  }
  ++stats().failed;
  std::printf("  FAIL %s:%d  %s  (", file, line, expr);
  print_operand(a);
  std::printf(" vs ");
  print_operand(b);
  std::printf(")\n");
}

/**
 * @brief Records a near-equality result and, on failure, prints operands and
 * delta.
 * @param a First value.
 * @param b Second value.
 * @param tol Maximum allowed absolute difference, in the same units as a and b.
 * @param expr Stringified comparison expression for the failure message.
 * @param file Source file name of the call site.
 * @param line Source line number of the call site.
 */
inline void report_near(double a, double b, double tol, const char *expr,
                        const char *file, int line) {
  if (approx(a, b, tol)) {
    ++stats().passed;
    return;
  }
  ++stats().failed;
  std::printf("  FAIL %s:%d  %s  (%g vs %g, delta=%g)\n", file, line, expr, a,
              b, std::fabs(a - b));
}

/**
 * @brief Snapshot of the global counter taken when a module starts.
 * @details Lets end_module report that module's tally as a delta from the
 * shared process-wide counter.
 */
struct ModuleScope {
  const char *name;    /**< Module name, echoed in the header and footer. */
  int passed_before;   /**< Global passed count captured at begin_module. */
  int failed_before;   /**< Global failed count captured at begin_module. */
};

/**
 * @brief Prints a module header and captures the current counter baseline.
 * @param name Module name to print and store in the scope.
 * @return A ModuleScope holding the name and the baseline pass/fail counts.
 */
inline ModuleScope begin_module(const char *name) {
  std::printf("=== %s ===\n", name);
  return {name, stats().passed, stats().failed};
}

/**
 * @brief Prints the module's pass/fail delta since begin_module.
 * @param m The scope returned by begin_module for this module.
 * @return The module's failure count (delta since begin_module).
 */
inline int end_module(const ModuleScope &m) {
  int passed = stats().passed - m.passed_before;
  int failed = stats().failed - m.failed_before;
  std::printf("=== %s: %d passed, %d failed ===\n", m.name, passed, failed);
  return failed;
}

} // namespace hs_test

// Core assertion all other HS_EXPECT_* macros funnel through.
#define HS_EXPECT(cond, msg)                                                   \
  do {                                                                         \
    if (cond) {                                                                \
      ++hs_test::stats().passed;                                               \
    } else {                                                                   \
      ++hs_test::stats().failed;                                               \
      std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, msg);              \
    }                                                                          \
  } while (0)

// Near-equality with an absolute tolerance. Operands are captured as double:
// floats widen exactly, and double-valued operands compare at full precision
// rather than being rounded to float resolution. Call sites wanting a
// float-domain compare invoke hs_test::approx(float,float,float) directly.
#define HS_EXPECT_NEAR(a, b, tol)                                             \
  do {                                                                        \
    double _hs_a = (a);                                                       \
    double _hs_b = (b);                                                       \
    hs_test::report_near(_hs_a, _hs_b, (tol),                                 \
                         #a " ~= " #b " (tol=" #tol ")", __FILE__, __LINE__); \
  } while (0)
// Capture each operand once so loop-driven assertions don't re-evaluate side
// effects, then compare and (on failure) print both values.
#define HS_EXPECT_CMP(a, b, op, opstr)                                        \
  do {                                                                        \
    auto _hs_a = (a);                                                         \
    auto _hs_b = (b);                                                         \
    hs_test::report_cmp(_hs_a op _hs_b, _hs_a, _hs_b, #a " " opstr " " #b,    \
                        __FILE__, __LINE__);                                  \
  } while (0)
#define HS_EXPECT_TRUE(cond) HS_EXPECT((cond), #cond)
#define HS_EXPECT_FALSE(cond) HS_EXPECT(!(cond), "!(" #cond ")")
#define HS_EXPECT_EQ(a, b) HS_EXPECT_CMP(a, b, ==, "==")
#define HS_EXPECT_NE(a, b) HS_EXPECT_CMP(a, b, !=, "!=")
#define HS_EXPECT_LT(a, b) HS_EXPECT_CMP(a, b, <, "<")
#define HS_EXPECT_LE(a, b) HS_EXPECT_CMP(a, b, <=, "<=")
#define HS_EXPECT_GT(a, b) HS_EXPECT_CMP(a, b, >, ">")
#define HS_EXPECT_GE(a, b) HS_EXPECT_CMP(a, b, >=, ">=")

