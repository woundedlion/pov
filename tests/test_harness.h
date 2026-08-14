/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Shared test harness — single global Stats counter, HS_EXPECT_* macros,
 * and per-module begin/end scope helpers. Per-module reporters compute their
 * own pass/fail counts as deltas from the global counter so multiple test
 * suites can coexist in one driver.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <type_traits>

namespace hs_test {

/** FNV-1a 64-bit offset basis. */
inline constexpr uint64_t FNV1A64_BASIS = 1469598103934665603ull;
/** FNV-1a 64-bit prime. */
inline constexpr uint64_t FNV1A64_PRIME = 1099511628211ull;

/**
 * @brief Folds one byte into a running FNV-1a 64-bit hash.
 * @param hash Running hash, seeded from FNV1A64_BASIS.
 * @param byte Byte to absorb.
 * @return The updated hash.
 */
inline constexpr uint64_t fnv1a64_byte(uint64_t hash, uint8_t byte) {
  return (hash ^ byte) * FNV1A64_PRIME;
}

/**
 * @brief Folds a 16-bit channel, low byte first.
 * @param hash Running hash.
 * @param channel Channel to absorb.
 * @return The updated hash.
 */
inline constexpr uint64_t fnv1a64_channel(uint64_t hash, uint16_t channel) {
  return fnv1a64_byte(fnv1a64_byte(hash, static_cast<uint8_t>(channel)),
                      static_cast<uint8_t>(channel >> 8));
}

/**
 * @brief Folds a byte span into a running FNV-1a 64-bit hash.
 * @param data Span start.
 * @param n Byte count.
 * @param hash Seed, letting a key be built from several spans.
 * @return The updated hash.
 */
inline uint64_t fnv1a64_bytes(const void *data, size_t n,
                              uint64_t hash = FNV1A64_BASIS) {
  const uint8_t *p = static_cast<const uint8_t *>(data);
  for (size_t i = 0; i < n; ++i)
    hash = fnv1a64_byte(hash, p[i]);
  return hash;
}

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

/**
 * @brief Tests finite floats with a tolerance relative to their magnitude.
 */
inline bool approx_rel(float a, float b, float rel_tol) {
  return std::isfinite(a) && std::isfinite(b) &&
         std::abs(a - b) <= rel_tol * std::max(std::abs(a), std::abs(b));
}

/**
 * @brief Tests finite doubles with a tolerance relative to their magnitude.
 */
inline bool approx_rel(double a, double b, double rel_tol) {
  return std::isfinite(a) && std::isfinite(b) &&
         std::abs(a - b) <= rel_tol * std::max(std::abs(a), std::abs(b));
}

namespace detail {
/**
 * @brief Structural detection for the engine's small value types so
 * print_operand can show components without test_harness.h depending on the
 * engine headers. Pixel exposes r/g/b; Vector exposes x/y/z; Quaternion
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
 * @brief One frame of the scoped assertion-context stack.
 * @details Holds a label pointer and up to two integer coordinates rather than
 * formatted text: a frame pushed inside a per-pixel loop must not cost a
 * vsnprintf on the passing path, so formatting happens only when a failure
 * prints.
 */
struct ContextFrame {
  const char *label = nullptr; /**< Caller-owned static label. */
  long long a = 0; /**< First coordinate, printed when arity >= 1. */
  long long b = 0; /**< Second coordinate, printed when arity == 2. */
  int arity = 0;   /**< Number of coordinates carried, 0..2. */
};

/** @brief Context frames retained; deeper pushes are counted, not stored. */
inline constexpr int CONTEXT_CAP = 8;

/** @brief Accessor for the process-wide context frame array. */
inline ContextFrame *context_frames() {
  static ContextFrame frames[CONTEXT_CAP];
  return frames;
}

/** @brief Accessor for the live context depth (may exceed CONTEXT_CAP). */
inline int &context_depth() {
  static int d = 0;
  return d;
}

/** @brief Prints the live context frames as " {outer | inner=(x, y)}". */
inline void print_context() {
  const int n = std::min(context_depth(), CONTEXT_CAP);
  for (int i = 0; i < n; ++i) {
    const ContextFrame &f = context_frames()[i];
    std::printf("%s%s", i == 0 ? "  {" : " | ", f.label ? f.label : "?");
    if (f.arity == 1)
      std::printf("=%lld", f.a);
    else if (f.arity == 2)
      std::printf("=(%lld, %lld)", f.a, f.b);
  }
  if (n > 0)
    std::printf("}");
}

/**
 * @brief RAII label naming what the enclosing block is asserting about.
 * @details HS_EXPECT_* stamps __func__, which inside a shared helper names the
 * helper rather than the call that reached it. Every live scope is printed with
 * each failure, so a helper's assertions say which caller and which sample they
 * came from.
 */
class ContextScope {
public:
  explicit ContextScope(const char *label) { push(label, 0, 0, 0); }
  ContextScope(const char *label, long long a) { push(label, a, 0, 1); }
  ContextScope(const char *label, long long a, long long b) {
    push(label, a, b, 2);
  }
  ~ContextScope() { --context_depth(); }
  ContextScope(const ContextScope &) = delete;
  ContextScope &operator=(const ContextScope &) = delete;

private:
  static void push(const char *label, long long a, long long b, int arity) {
    const int d = context_depth()++;
    if (d < CONTEXT_CAP)
      context_frames()[d] = ContextFrame{label, a, b, arity};
  }
};

/**
 * @brief Prints a single comparison operand in a type-appropriate format.
 * @tparam T Operand type; bool, floating-point, enum, integral, pointer, and the
 * engine's r/g/b and x/y/z value types are formatted specially.
 * @param v The operand value to print.
 * @details Lets a failing HS_EXPECT_* line show the actual values, not just the
 * stringified expr. Pixel (r,g,b), Vector (x,y,z), Quaternion (r,v),
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
  std::printf(")");
  print_context();
  std::printf("\n");
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
  std::printf("  FAIL [%s] %s:%d  %s  (%g vs %g, delta=%g)", case_name, file,
              line, expr, a, b, std::fabs(a - b));
  print_context();
  std::printf("\n");
}

/**
 * @brief Records a relative near-equality result and prints its relative error.
 */
inline void report_near_rel(double a, double b, double rel_tol,
                            const char *expr, const char *case_name,
                            const char *file, int line) {
  if (approx_rel(a, b, rel_tol)) {
    ++stats().passed;
    return;
  }
  ++stats().failed;
  if (!claim_fail_print())
    return;
  const double scale = std::max(std::fabs(a), std::fabs(b));
  const double rel_error =
      scale > 0.0 ? std::fabs(a - b) / scale : std::fabs(a - b);
  std::printf("  FAIL [%s] %s:%d  %s  (%g vs %g, rel_error=%g)", case_name,
              file, line, expr, a, b, rel_error);
  print_context();
  std::printf("\n");
}

/**
 * @brief Accessor for the "this module skipped its entire case set" flag.
 * @return Reference to the flag; begin_module clears it for every module.
 */
inline bool &suite_skipped() {
  static bool s = false;
  return s;
}

/**
 * @brief Records a skip that stands for the module's ENTIRE case set.
 * @details Tallies the skip and waives this module's assertion floor — the
 * module ran no cases, so a floor measured over them cannot apply. This is the
 * only thing that waives a floor.
 */
inline void skip_suite() {
  ++stats().skipped;
  suite_skipped() = true;
}

/**
 * @brief Records one skipped sub-case, leaving the module's floor in force.
 * @details A module that skips one case still ran the rest, so the floor must
 * still gate them. Subtract the skipped case's assertions from the roster floor
 * if the skip is conditional; do not reach for skip_suite() to silence it.
 */
inline void skip_case() { ++stats().skipped; }

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
  suite_skipped() = false;
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
 * roster recorded for it. Only a whole-suite skip (skip_suite()) is exempt; a
 * skipped sub-case leaves the floor gating every case the module still ran.
 */
inline int end_module(const ModuleScope &m) {
  int passed = stats().passed - m.passed_before;
  int failed = stats().failed - m.failed_before;
  int skipped = stats().skipped - m.skipped_before;
  // A module that ran no assertion and did not skip its whole case set did no
  // work; count it as a failure so an emptied runner goes red, not green.
  if (passed + failed == 0 && !suite_skipped()) {
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
  // Only a whole-suite skip waives the floor (death skips its whole case loop
  // when it cannot re-exec). A per-case skip must not, or the first module to
  // record one would silently lose the floor over all its other cases.
  if (m.min_assertions > 0 && !suite_skipped() &&
      passed + failed < m.min_assertions) {
    std::printf("=== %s: only %d assertions ran, expected >= %d (a case call "
                "was dropped) ===\n",
                m.name, passed + failed, m.min_assertions);
    return failed + 1;
  }
  return failed;
}

} // namespace hs_test

#define HS_CTX_JOIN_INNER(a, b) a##b
#define HS_CTX_JOIN(a, b) HS_CTX_JOIN_INNER(a, b)
// Scoped label printed with every failure raised inside the enclosing block:
// a name, optionally followed by one or two integer coordinates. Use it in a
// shared helper (or its callers) where __func__ alone does not identify the
// sample that failed.
#define HS_CONTEXT(...)                                                        \
  const hs_test::ContextScope HS_CTX_JOIN(hs_ctx, __LINE__)(__VA_ARGS__)

// Core assertion all other HS_EXPECT_* macros funnel through.
#define HS_EXPECT(cond, msg)                                                   \
  do {                                                                         \
    if (cond) {                                                                \
      ++hs_test::stats().passed;                                               \
    } else {                                                                   \
      ++hs_test::stats().failed;                                               \
      if (hs_test::claim_fail_print()) {                                       \
        std::printf("  FAIL [%s] %s:%d  %s", __func__, __FILE__, __LINE__,     \
                    msg);                                                      \
        hs_test::print_context();                                              \
        std::printf("\n");                                                     \
      }                                                                        \
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
#define HS_EXPECT_NEAR_REL(a, b, rel_tol)                                      \
  do {                                                                         \
    double _hs_a = (a);                                                        \
    double _hs_b = (b);                                                        \
    hs_test::report_near_rel(_hs_a, _hs_b, (rel_tol),                          \
                             #a " ~= " #b " (rel_tol=" #rel_tol ")", __func__, \
                             __FILE__, __LINE__);                              \
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
// Fatal size check for a case that indexes a container afterwards: records the
// comparison like HS_EXPECT_EQ, then returns from the enclosing void case on a
// mismatch so the indexing that follows cannot read past the end. The shortened
// run is already red from the recorded failure, so the module's assertion floor
// stays consistent.
#define HS_EXPECT_SIZE_OR_RETURN(container, expected)                          \
  do {                                                                         \
    const std::size_t _hs_size = (container).size();                           \
    const std::size_t _hs_want = static_cast<std::size_t>(expected);           \
    hs_test::report_cmp(_hs_size == _hs_want, _hs_size, _hs_want,              \
                        #container ".size() == " #expected, __func__,          \
                        __FILE__, __LINE__);                                   \
    if (_hs_size != _hs_want)                                                  \
      return;                                                                  \
  } while (0)
