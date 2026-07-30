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
  int skipped =
      0; /**< Count of suites/cases skipped (never counted as passed). */
};

/**
 * @brief Accessor for the single process-wide pass/fail counter.
 * @return Reference to the function-local static Stats shared by all suites.
 */
inline Stats &stats() {
  static Stats s;
  return s;
}

/** @brief FAIL lines printed per module before further failures print nothing. */
inline constexpr int FAIL_PRINT_CAP = 50;

/**
 * @brief Per-module FAIL-line print budget.
 * @details Only printing is capped; Stats::failed still counts every failure.
 * A per-element assertion loop over a multi-million-element suite would
 * otherwise emit one line per element, and ctest buffers a test's whole output
 * before printing it.
 */
struct FailPrintBudget {
  int printed = 0;    /**< FAIL lines printed since the current begin_module. */
  int suppressed = 0; /**< Failures whose line was withheld by the cap. */
};

/** @brief Accessor for the process-wide FAIL-line print budget. */
inline FailPrintBudget &fail_print_budget() {
  static FailPrintBudget b;
  return b;
}

/**
 * @brief Claims one FAIL line from the current module's print budget.
 * @return True if the caller should print; false after the cap, having tallied
 * the failure as suppressed.
 */
inline bool claim_fail_print() {
  FailPrintBudget &b = fail_print_budget();
  if (b.printed >= FAIL_PRINT_CAP) {
    ++b.suppressed;
    return false;
  }
  ++b.printed;
  return true;
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
 * engine headers. Pixel/Pixel16 expose r/g/b; Vector exposes x/y/z; Quaternion
 * exposes r/v; Color4 exposes color/alpha.
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

template <class T, class = void> struct has_r_v : std::false_type {};
template <class T>
struct has_r_v<T, std::void_t<decltype(std::declval<const T &>().r),
                              decltype(std::declval<const T &>().v)>>
    : std::true_type {};

template <class T, class = void> struct has_color_alpha : std::false_type {};
template <class T>
struct has_color_alpha<T,
                       std::void_t<decltype(std::declval<const T &>().color),
                                   decltype(std::declval<const T &>().alpha)>>
    : std::true_type {};
} // namespace detail

/**
 * @brief Prints a single comparison operand in a type-appropriate format.
 * @tparam T Operand type; bool, floating-point, enum, integral, pointer, and the
 * engine's r/g/b and x/y/z value types are formatted specially.
 * @param v The operand value to print.
 * @details Lets a failing HS_EXPECT_* line show the actual values, not just the
 * stringified expr. Pixel/Pixel16 (r,g,b), Vector (x,y,z), Quaternion (r,v),
 * and Color4 (color,alpha) are detected structurally — so an HS_EXPECT_EQ/CMP on
 * those prints components rather than "?" — and their fields recurse through this
 * same printer. Any other non-arithmetic operand still falls back to "?".
 */
template <class T> inline void print_operand(const T &v) {
  if constexpr (std::is_same_v<T, bool>) {
    std::printf("%s", v ? "true" : "false");
  } else if constexpr (std::is_floating_point_v<T>) {
    std::printf("%g", static_cast<double>(v));
  } else if constexpr (std::is_enum_v<T>) {
    std::printf("%lld", static_cast<long long>(
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
  } else if constexpr (detail::has_r_v<T>::value) {
    std::printf("(r=");
    print_operand(v.r);
    std::printf(" v=");
    print_operand(v.v);
    std::printf(")");
  } else if constexpr (detail::has_color_alpha<T>::value) {
    std::printf("(color=");
    print_operand(v.color);
    std::printf(" alpha=");
    print_operand(v.alpha);
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
 * @param case_name Test function containing the assertion.
 * @param file Source file name of the call site.
 * @param line Source line number of the call site.
 */
template <class A, class B>
inline void report_cmp(bool ok, const A &a, const B &b, const char *expr,
                       const char *case_name, const char *file, int line) {
  if (ok) {
    ++stats().passed;
    return;
  }
  ++stats().failed;
  if (!claim_fail_print())
    return;
  std::printf("  FAIL [%s] %s:%d  %s  (", case_name, file, line, expr);
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
 * @param case_name Test function containing the assertion.
 * @param file Source file name of the call site.
 * @param line Source line number of the call site.
 */
inline void report_near(double a, double b, double tol, const char *expr,
                        const char *case_name, const char *file, int line) {
  if (approx(a, b, tol)) {
    ++stats().passed;
    return;
  }
  ++stats().failed;
  if (!claim_fail_print())
    return;
  std::printf("  FAIL [%s] %s:%d  %s  (%g vs %g, delta=%g)\n", case_name, file,
              line, expr, a, b, std::fabs(a - b));
}

/**
 * @brief Snapshot of the global counter taken when a module starts.
 * @details Lets end_module report that module's tally as a delta from the
 * shared process-wide counter.
 */
struct ModuleScope {
  const char *name;   /**< Module name, echoed in the header and footer. */
  int passed_before;  /**< Global passed count captured at begin_module. */
  int failed_before;  /**< Global failed count captured at begin_module. */
  int skipped_before; /**< Global skipped count captured at begin_module. */
  int min_assertions; /**< Floor enforced by end_module; 0 disables the check. */
};

/**
 * @brief Accessor for the assertion floor the next begin_module will adopt.
 * @return Reference to the pending floor, 0 when no floor is published.
 * @details The roster in run_tests.cpp publishes a module's measured floor here
 * immediately before invoking it. begin_module consumes the value and resets it
 * to 0, so a scope opened outside the roster driver — a standalone check binary —
 * carries no floor.
 */
inline int &pending_min_assertions() {
  static int m = 0;
  return m;
}

/**
 * @brief Prints a module header and captures the current counter baseline.
 * @param name Module name to print and store in the scope.
 * @return A ModuleScope holding the name, the baseline pass/fail counts, and the
 * pending assertion floor.
 */
inline ModuleScope begin_module(const char *name) {
  std::printf("=== %s ===\n", name);
  fail_print_budget() = FailPrintBudget{};
  const int min_assertions = pending_min_assertions();
  pending_min_assertions() = 0;
  return {name, stats().passed, stats().failed, stats().skipped,
          min_assertions};
}

/**
 * @brief Prints the module's pass/fail delta since begin_module and enforces its
 * assertion floor.
 * @param m The scope returned by begin_module for this module.
 * @return The module's failure count (delta since begin_module), plus one if the
 * module ran fewer assertions than its floor.
 * @details Individual cases are invoked by hand-written calls, so a case that is
 * defined and never called compiles clean. The floor turns that into a red run:
 * dropping a call removes assertions, and the module falls below the count the
 * roster recorded for it.
 */
inline int end_module(const ModuleScope &m) {
  int passed = stats().passed - m.passed_before;
  int failed = stats().failed - m.failed_before;
  int skipped = stats().skipped - m.skipped_before;
  // A module that ran no assertion and skipped nothing did no work; count it as
  // a failure so an emptied runner goes red instead of silently green.
  if (passed + failed == 0 && skipped == 0) {
    std::printf("=== %s: NO ASSERTIONS RAN — counting as FAILURE ===\n",
                m.name);
    return 1;
  }
  const int suppressed = fail_print_budget().suppressed;
  if (suppressed > 0)
    std::printf(
        "  ... %d further FAIL line(s) suppressed (cap %d per module)\n",
        suppressed, FAIL_PRINT_CAP);
  if (skipped > 0)
    std::printf("=== %s: %d passed, %d failed, %d SKIPPED ===\n", m.name,
                passed, failed, skipped);
  else
    std::printf("=== %s: %d passed, %d failed ===\n", m.name, passed, failed);
  if (m.min_assertions > 0 && passed + failed < m.min_assertions) {
    std::printf("=== %s: only %d assertions ran, expected >= %d (a case call "
                "was dropped) ===\n",
                m.name, passed + failed, m.min_assertions);
    return failed + 1;
  }
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
      if (hs_test::claim_fail_print())                                         \
        std::printf("  FAIL [%s] %s:%d  %s\n", __func__, __FILE__, __LINE__,   \
                    msg);                                                      \
    }                                                                          \
  } while (0)

// Near-equality with an absolute tolerance. Operands are captured as double so
// double-valued operands compare at full precision; a float-domain compare can
// invoke hs_test::approx(float,float,float) directly.
#define HS_EXPECT_NEAR(a, b, tol)                                              \
  do {                                                                         \
    double _hs_a = (a);                                                        \
    double _hs_b = (b);                                                        \
    hs_test::report_near(_hs_a, _hs_b, (tol), #a " ~= " #b " (tol=" #tol ")",  \
                         __func__, __FILE__, __LINE__);                        \
  } while (0)
// Capture each operand once so loop-driven assertions don't re-evaluate side
// effects, then compare and (on failure) print both values.
#define HS_EXPECT_CMP(a, b, op, opstr)                                         \
  do {                                                                         \
    auto _hs_a = (a);                                                          \
    auto _hs_b = (b);                                                          \
    hs_test::report_cmp(_hs_a op _hs_b, _hs_a, _hs_b, #a " " opstr " " #b,     \
                        __func__, __FILE__, __LINE__);                         \
  } while (0)
#define HS_EXPECT_TRUE(cond) HS_EXPECT((cond), #cond)
#define HS_EXPECT_FALSE(cond) HS_EXPECT(!(cond), "!(" #cond ")")
#define HS_EXPECT_EQ(a, b) HS_EXPECT_CMP(a, b, ==, "==")
#define HS_EXPECT_NE(a, b) HS_EXPECT_CMP(a, b, !=, "!=")
#define HS_EXPECT_LT(a, b) HS_EXPECT_CMP(a, b, <, "<")
#define HS_EXPECT_LE(a, b) HS_EXPECT_CMP(a, b, <=, "<=")
#define HS_EXPECT_GT(a, b) HS_EXPECT_CMP(a, b, >, ">")
#define HS_EXPECT_GE(a, b) HS_EXPECT_CMP(a, b, >=, ">=")
