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
#ifndef HOLOSPHERE_TESTS_TEST_HARNESS_H_
#define HOLOSPHERE_TESTS_TEST_HARNESS_H_

#include <cmath>
#include <cstdio>

namespace hs_test {

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

#define HS_EXPECT_NEAR(a, b, tol)                                              \
  HS_EXPECT(hs_test::approx((a), (b), (tol)),                                  \
            #a " ~= " #b " (tol=" #tol ")")
#define HS_EXPECT_TRUE(cond) HS_EXPECT((cond), #cond)
#define HS_EXPECT_FALSE(cond) HS_EXPECT(!(cond), "!(" #cond ")")
#define HS_EXPECT_EQ(a, b) HS_EXPECT((a) == (b), #a " == " #b)
#define HS_EXPECT_NE(a, b) HS_EXPECT((a) != (b), #a " != " #b)
#define HS_EXPECT_LT(a, b) HS_EXPECT((a) < (b), #a " < " #b)
#define HS_EXPECT_LE(a, b) HS_EXPECT((a) <= (b), #a " <= " #b)
#define HS_EXPECT_GT(a, b) HS_EXPECT((a) > (b), #a " > " #b)
#define HS_EXPECT_GE(a, b) HS_EXPECT((a) >= (b), #a " >= " #b)

#endif // HOLOSPHERE_TESTS_TEST_HARNESS_H_
