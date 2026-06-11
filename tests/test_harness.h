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

// Single-threaded by contract: the counters are plain ints, so every
// HS_EXPECT_* must be evaluated on the test's main thread. This keeps the hot
// assertion path free of atomic RMWs across the whole suite. The one test that
// spawns a helper thread (test_ctor_spin_waits_for_buffer_free) honours this by
// capturing its cross-thread observation into a std::atomic and asserting on the
// main thread after join() — it never calls HS_EXPECT_* from the helper.
struct Stats {
  int passed = 0;
  int failed = 0;
};

inline Stats &stats() {
  static Stats s;
  return s;
}

inline bool approx(float a, float b, float tol) {
  return std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= tol;
}

// Double-domain overload. HS_EXPECT_NEAR captures its operands as double, so the
// comparison must stay in double — downcasting to float would round both sides
// to ~1e-7 resolution and silently mask a sub-float-epsilon regression in a
// double-valued expression. Float call sites bind the float overload above by
// exact match, so this is additive.
inline bool approx(double a, double b, double tol) {
  return std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= tol;
}

// Print a single comparison operand in a type-appropriate format, so a failing
// HS_EXPECT_* line shows the actual values (not just the stringified expr).
// Non-arithmetic operands fall back to "?" — every comparison-macro call site
// in the suite passes arithmetic/enum values.
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
  } else {
    std::printf("?");
  }
}

// Record the result of a comparison and, on failure, print both operands.
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

// Record a near-equality result and, on failure, print operands and delta.
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

struct ModuleScope {
  const char *name;
  int passed_before;
  int failed_before;
};

inline ModuleScope begin_module(const char *name) {
  std::printf("=== %s ===\n", name);
  return {name, stats().passed, stats().failed};
}

inline int end_module(const ModuleScope &m) {
  int passed = stats().passed - m.passed_before;
  int failed = stats().failed - m.failed_before;
  std::printf("=== %s: %d passed, %d failed ===\n", m.name, passed, failed);
  return failed;
}

} // namespace hs_test

#define HS_EXPECT(cond, msg)                                                   \
  do {                                                                         \
    if (cond) {                                                                \
      ++hs_test::stats().passed;                                               \
    } else {                                                                   \
      ++hs_test::stats().failed;                                               \
      std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, msg);              \
    }                                                                          \
  } while (0)

#define HS_EXPECT_NEAR(a, b, tol)                                             \
  do {                                                                        \
    double _hs_a = (a);                                                       \
    double _hs_b = (b);                                                       \
    hs_test::report_near(_hs_a, _hs_b, (tol),                                 \
                         #a " ~= " #b " (tol=" #tol ")", __FILE__, __LINE__); \
  } while (0)
// Capture each operand once into a local so loop-driven assertions don't
// re-evaluate side effects, then compare and (on failure) print both values.
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

