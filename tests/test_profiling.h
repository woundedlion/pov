/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/engine/profiling.h — the instrumentation every on-device
 * timing report is built from. None of it renders a pixel, so a regression
 * surfaces as wrong numbers in a report rather than as a visible failure. These
 * cases pin the hand-rolled u64_dec formatter against its exact-fit buffer, the
 * CycleCounter registry walks (find_suffix, reset_all and the log_node tree
 * captured off stdout), the parent attribution CycleScope builds, and the
 * ISR-side accumulator.
 *
 * Self-contained header. run_profiling_tests() returns the module failure count.
 */
#pragma once

#include "core/engine/platform.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#if defined(_WIN32)
#include <io.h>      // _dup / _dup2 / _close / _fileno
#include <process.h> // _getpid
#else
#include <unistd.h> // dup / dup2 / close / getpid
#endif

namespace hs_test {
namespace profiling_tests {

/**
 * @brief Runs CycleCounter::log_all() with fd 1 redirected into @p out.
 * @param out Destination buffer; receives the NUL-terminated report.
 * @param n Capacity of @p out.
 * @return True if the redirect was established and the report captured.
 * @details hs::log writes to C stdout, so the report is only observable through
 * an fd swap (the same idiom as test_platform.h's Serial.printf capture). A
 * scratch file rather than that capture's pipe, since the report can outrun a
 * pipe buffer and the write precedes any read. The file lands in the platform
 * temp directory, so an unwritable CWD cannot fail the cases, and its name
 * carries the pid so concurrent runs cannot collide.
 */
inline bool capture_log_all(char *out, size_t n) {
  out[0] = '\0';
  char path[512];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#if defined(_WIN32)
  const char *dir = std::getenv("TEMP");
  if (!dir || dir[0] == '\0')
    dir = ".";
  std::snprintf(path, sizeof(path), "%s\\profiling_log_capture_%d.tmp", dir,
                _getpid());
#else
  const char *dir = std::getenv("TMPDIR");
  if (!dir || dir[0] == '\0')
    dir = "/tmp";
  std::snprintf(path, sizeof(path), "%s/profiling_log_capture_%d.tmp", dir,
                static_cast<int>(getpid()));
#endif
#pragma clang diagnostic pop
  std::fflush(stdout);
#if defined(_WIN32)
  int saved_out = _dup(1);
  std::FILE *cap = nullptr;
  fopen_s(&cap, path, "w+");
#else
  int saved_out = dup(1);
  std::FILE *cap = std::fopen(path, "w+");
#endif
  if (cap == nullptr || saved_out < 0) {
    if (cap != nullptr) {
      std::fclose(cap);
      std::remove(path);
    }
    if (saved_out >= 0) {
#if defined(_WIN32)
      _close(saved_out);
#else
      close(saved_out);
#endif
    }
    return false;
  }
#if defined(_WIN32)
  _dup2(_fileno(cap), 1);
#else
  dup2(fileno(cap), 1);
#endif
  hs::CycleCounter::log_all();
  std::fflush(stdout);
#if defined(_WIN32)
  _dup2(saved_out, 1);
  _close(saved_out);
#else
  dup2(saved_out, 1);
  close(saved_out);
#endif
  std::rewind(cap);
  size_t got = std::fread(out, 1, n - 1, cap);
  out[got] = '\0';
  std::fclose(cap);
  std::remove(path);
  return true;
}

// --- u64_dec ----------------------------------------------------------------

/**
 * @brief u64_dec's exact-fit 21-byte buffer between two canary bands.
 * @details UINT64_MAX fills the buffer to its first byte, so an off-by-one in
 * the digit loop writes outside it; the bands make that visible without relying
 * on a sanitizer being enabled.
 */
struct GuardedDecBuf {
  static constexpr char CANARY = 0x5a;

  char pre[8];
  char buf[21];
  char post[8];

  GuardedDecBuf() {
    std::memset(pre, CANARY, sizeof(pre));
    std::memset(buf, CANARY, sizeof(buf));
    std::memset(post, CANARY, sizeof(post));
  }

  /** @brief True while neither band has been written. */
  bool guards_intact() const {
    for (char c : pre)
      if (c != CANARY)
        return false;
    for (char c : post)
      if (c != CANARY)
        return false;
    return true;
  }
};

/**
 * @brief Formats one value and checks its digits, terminator, start offset and
 *        buffer containment.
 * @param v Value to format.
 * @param expected Expected decimal rendering.
 */
inline void expect_u64_dec(uint64_t v, const char *expected) {
  GuardedDecBuf g;
  const char *s = hs::u64_dec(v, g.buf);
  HS_EXPECT(std::strcmp(s, expected) == 0, "u64_dec renders the value");
  HS_EXPECT_EQ(g.buf[20], '\0');
  // Digits are written backwards from the fixed terminator slot.
  HS_EXPECT_EQ(static_cast<size_t>(s - g.buf), 20u - std::strlen(expected));
  HS_EXPECT_TRUE(g.guards_intact());
}

/**
 * @brief Verifies u64_dec across the digit-count boundaries, including zero and
 *        the 20-digit value that exactly fills its buffer.
 */
inline void test_u64_dec_boundaries() {
  expect_u64_dec(0, "0");
  expect_u64_dec(1, "1");
  expect_u64_dec(9, "9");
  expect_u64_dec(10, "10");
  expect_u64_dec(99, "99");
  expect_u64_dec(100, "100");
  expect_u64_dec(999, "999");
  expect_u64_dec(1000, "1000");
  expect_u64_dec(4294967295u, "4294967295");
  expect_u64_dec(4294967296ull, "4294967296");
  expect_u64_dec(9999999999999999999ull, "9999999999999999999");
  expect_u64_dec(10000000000000000000ull, "10000000000000000000");
  expect_u64_dec(UINT64_MAX, "18446744073709551615");
}

/**
 * @brief Verifies every power of ten round-trips through u64_dec.
 */
inline void test_u64_dec_powers_of_ten() {
  uint64_t v = 1;
  for (int digits = 1; digits <= 19; ++digits) {
    GuardedDecBuf g;
    const char *s = hs::u64_dec(v, g.buf);
    HS_EXPECT_EQ(static_cast<int>(std::strlen(s)), digits);
    HS_EXPECT_EQ(s[0], '1');
    HS_EXPECT_TRUE(g.guards_intact());
    v *= 10;
  }
}

// --- registry ---------------------------------------------------------------

/**
 * @brief Verifies find_suffix matches whole names and proper suffixes, prefers
 *        the most recently registered match, and rejects a suffix no name can
 *        contain.
 */
inline void test_find_suffix() {
  // Registry entries are never unlinked, so every counter here must outlive the
  // process.
  static hs::CycleCounter alpha("prof_suffix_alpha_tail");
  static hs::CycleCounter beta("prof_suffix_beta_tail");
  static hs::CycleCounter gamma("prof_suffix_gamma");

  // Registration pushes onto the list head, so the newest match wins.
  HS_EXPECT_EQ(hs::CycleCounter::find_suffix("_tail"), &beta);
  HS_EXPECT_EQ(hs::CycleCounter::find_suffix("prof_suffix_gamma"), &gamma);
  HS_EXPECT_EQ(hs::CycleCounter::find_suffix("alpha_tail"), &alpha);
  HS_EXPECT_EQ(hs::CycleCounter::find_suffix(""), &gamma);
  HS_EXPECT_TRUE(hs::CycleCounter::find_suffix("_prof_no_such_counter") ==
                 nullptr);

  // Longer than any registered name: the length guard must reject it before the
  // compare reads behind the name.
  char long_suffix[200];
  std::memset(long_suffix, 'z', sizeof(long_suffix) - 1);
  long_suffix[sizeof(long_suffix) - 1] = '\0';
  HS_EXPECT_TRUE(hs::CycleCounter::find_suffix(long_suffix) == nullptr);
}

/**
 * @brief Verifies reset_all zeroes accumulated cycles, calls and the
 *        mixed-parent flag on every registered counter while leaving the tree
 *        structure in place.
 */
inline void test_reset_all_clears_counts() {
  static hs::CycleCounter outer("prof_reset_outer");
  static hs::CycleCounter inner("prof_reset_inner");
  static hs::CycleCounter other("prof_reset_other");
  {
    hs::CycleScope so(outer);
    hs::CycleScope si(inner);
  }
  {
    hs::CycleScope sx(other);
    hs::CycleScope si(inner);
  }
  HS_EXPECT_TRUE(inner.mixed_parent);
  outer.cycles = 4242;
  inner.cycles = 99;

  hs::CycleCounter::reset_all();

  HS_EXPECT_EQ(outer.cycles, 0u);
  HS_EXPECT_EQ(outer.count, 0u);
  HS_EXPECT_EQ(inner.cycles, 0u);
  HS_EXPECT_EQ(inner.count, 0u);
  // The flag describes the entries just discarded, so it cannot outlive them.
  HS_EXPECT_FALSE(inner.mixed_parent);
  HS_EXPECT_EQ(inner.parent, &outer);
  HS_EXPECT_TRUE(outer.parent == nullptr);
}

// --- nesting ----------------------------------------------------------------

/**
 * @brief Verifies a nested scope latches its enclosing counter as parent and
 *        leaves the outer counter a root.
 */
inline void test_nesting_latches_parent() {
  static hs::CycleCounter outer("prof_nest_outer");
  static hs::CycleCounter inner("prof_nest_inner");
  hs::CycleCounter::reset_all();
  {
    hs::CycleScope so(outer);
    hs::CycleScope si(inner);
  }
  HS_EXPECT_TRUE(outer.parent == nullptr);
  HS_EXPECT_EQ(inner.parent, &outer);
  HS_EXPECT_FALSE(outer.mixed_parent);
  HS_EXPECT_FALSE(inner.mixed_parent);
  HS_EXPECT_EQ(outer.count, 1u);
  HS_EXPECT_EQ(inner.count, 1u);
}

/**
 * @brief Verifies a counter re-entered through itself neither self-parents nor
 *        counts as a second caller.
 */
inline void test_recursive_scope_does_not_self_parent() {
  static hs::CycleCounter rec("prof_nest_recursive");
  hs::CycleCounter::reset_all();
  {
    hs::CycleScope a(rec);
    hs::CycleScope b(rec);
  }
  HS_EXPECT_TRUE(rec.parent == nullptr);
  HS_EXPECT_FALSE(rec.mixed_parent);
  HS_EXPECT_EQ(rec.count, 2u);
}

/**
 * @brief Verifies a counter entered under two different callers keeps its first
 *        parent and is flagged, and that the report tags the flagged node.
 */
inline void test_second_caller_flags_mixed_parent() {
  static hs::CycleCounter first("prof_mixed_first");
  static hs::CycleCounter second("prof_mixed_second");
  static hs::CycleCounter shared("prof_mixed_shared");
  hs::CycleCounter::reset_all();
  {
    hs::CycleScope sf(first);
    hs::CycleScope ss(shared);
  }
  {
    hs::CycleScope ss2(second);
    hs::CycleScope ss3(shared);
  }
  HS_EXPECT_EQ(shared.parent, &first);
  HS_EXPECT_TRUE(shared.mixed_parent);
  HS_EXPECT_FALSE(first.mixed_parent);
  HS_EXPECT_FALSE(second.mixed_parent);

  first.cycles = 600;
  second.cycles = 600;
  shared.cycles = 900; // more than any one parent spent: the tell-tale ratio
  char report[4096];
  if (!capture_log_all(report, sizeof(report))) {
    HS_EXPECT(false, "log_all capture set up");
    return;
  }
  HS_EXPECT_TRUE(std::strstr(report, "prof_mixed_shared") != nullptr);
  HS_EXPECT_TRUE(std::strstr(report, "MIXED-PARENT") != nullptr);
  HS_EXPECT_TRUE(std::strstr(report, "(150%)") != nullptr);
}

/**
 * @brief Verifies a top-level counter later entered under another counter is
 *        flagged rather than re-parented, so both stay on the root walk.
 */
inline void test_mutual_nesting_keeps_a_root() {
  static hs::CycleCounter left("prof_mutual_left");
  static hs::CycleCounter right("prof_mutual_right");
  hs::CycleCounter::reset_all();
  {
    hs::CycleScope sl(left);
    hs::CycleScope sr(right);
  }
  {
    hs::CycleScope sr(right);
    hs::CycleScope sl(left);
  }
  // Neither edge may be rewritten: a left<->right cycle would drop both from
  // log_all()'s root walk.
  HS_EXPECT_TRUE(left.parent == nullptr);
  HS_EXPECT_EQ(right.parent, &left);
  HS_EXPECT_TRUE(left.mixed_parent);
  HS_EXPECT_TRUE(right.mixed_parent);

  left.cycles = 100;
  right.cycles = 50;
  char report[4096];
  if (!capture_log_all(report, sizeof(report))) {
    HS_EXPECT(false, "log_all capture set up");
    return;
  }
  HS_EXPECT_TRUE(std::strstr(report, "prof_mutual_left") != nullptr);
  HS_EXPECT_TRUE(std::strstr(report, "prof_mutual_right") != nullptr);
}

// --- log_node ---------------------------------------------------------------

/**
 * @brief Verifies the tree walk prints a root and its child in order, with the
 *        child indented, its share taken against the parent, cycles converted
 *        to microseconds, and untouched counters omitted.
 */
inline void test_log_all_reports_tree() {
  static hs::CycleCounter root("prof_report_root");
  static hs::CycleCounter child("prof_report_child");
  static hs::CycleCounter idle("prof_report_idle");
  {
    hs::CycleScope sr(root);
    hs::CycleScope sc(child);
  }
  hs::CycleCounter::reset_all();
  root.cycles = 1200000;
  root.count = 4;
  child.cycles = 300000;
  child.count = 8;

  char report[4096];
  if (!capture_log_all(report, sizeof(report))) {
    HS_EXPECT(false, "log_all capture set up");
    return;
  }
  HS_EXPECT_EQ(idle.count, 0u);
  const char *root_line = std::strstr(report, "prof_report_root");
  const char *child_line = std::strstr(report, "prof_report_child");
  HS_EXPECT_TRUE(root_line != nullptr);
  HS_EXPECT_TRUE(child_line != nullptr);
  HS_EXPECT_TRUE(std::strstr(report, "prof_report_idle") == nullptr);
  if (root_line == nullptr || child_line == nullptr)
    return;
  HS_EXPECT_TRUE(root_line < child_line);
  HS_EXPECT_TRUE(child_line - report >= 2 && child_line[-1] == ' ' &&
                 child_line[-2] == ' ');
  HS_EXPECT_TRUE(
      std::strstr(root_line, "2000 us (100%)  4 calls  1200000 cyc") !=
      nullptr);
  HS_EXPECT_TRUE(std::strstr(child_line, "500 us (25%)  8 calls  300000 cyc") !=
                 nullptr);
  HS_EXPECT_TRUE(std::strstr(report, "--- Cycle Counters ---") != nullptr);
}

/**
 * @brief Verifies a counter whose latched parent records nothing in the next
 *        run is still reported, as a root, rather than dropped with it.
 */
inline void test_reset_does_not_orphan_subtree() {
  static hs::CycleCounter absent("prof_orphan_parent");
  static hs::CycleCounter kept("prof_orphan_child");
  {
    hs::CycleScope sa(absent);
    hs::CycleScope sk(kept);
  }
  hs::CycleCounter::reset_all();
  kept.cycles = 600;
  kept.count = 1;

  char report[4096];
  if (!capture_log_all(report, sizeof(report))) {
    HS_EXPECT(false, "log_all capture set up");
    return;
  }
  HS_EXPECT_EQ(kept.parent, &absent);
  HS_EXPECT_EQ(absent.count, 0u);
  HS_EXPECT_TRUE(std::strstr(report, "prof_orphan_child") != nullptr);
  HS_EXPECT_TRUE(std::strstr(report, "prof_orphan_parent") == nullptr);
}

// --- IsrCycleStats ----------------------------------------------------------

/**
 * @brief Verifies the ISR accumulator sums cycles and tracks the shortest and
 *        longest scope, and that reset() restores the empty sentinels.
 */
inline void test_isr_cycle_stats() {
  hs::IsrCycleStats stats;
  HS_EXPECT_EQ(stats.count, 0u);
  HS_EXPECT_EQ(stats.min, UINT32_MAX);
  HS_EXPECT_EQ(stats.max, 0u);

  stats.add(500);
  stats.add(20);
  stats.add(9000);
  HS_EXPECT_EQ(stats.cycles, 9520u);
  HS_EXPECT_EQ(stats.count, 3u);
  HS_EXPECT_EQ(stats.min, 20u);
  HS_EXPECT_EQ(stats.max, 9000u);

  stats.reset();
  HS_EXPECT_EQ(stats.cycles, 0u);
  HS_EXPECT_EQ(stats.count, 0u);
  HS_EXPECT_EQ(stats.min, UINT32_MAX);
  HS_EXPECT_EQ(stats.max, 0u);
}

/**
 * @brief Runs every profiling test case.
 * @return The module's failure count.
 */
inline int run_profiling_tests() {
  hs_test::ModuleFixture fixture("profiling");

  test_u64_dec_boundaries();
  test_u64_dec_powers_of_ten();
  test_find_suffix();
  test_reset_all_clears_counts();
  test_nesting_latches_parent();
  test_recursive_scope_does_not_self_parent();
  test_second_caller_flags_mixed_parent();
  test_mutual_nesting_keeps_a_root();
  test_log_all_reports_tree();
  test_reset_does_not_orphan_subtree();
  test_isr_cycle_stats();

  return fixture.result();
}

} // namespace profiling_tests
} // namespace hs_test
