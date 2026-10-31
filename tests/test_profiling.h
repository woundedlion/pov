/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/platform/cycle_counting.h, tools/profile_isr_window.h
 * and tools/pc_sampler.h.
 */
#pragma once

#include "core/platform/platform.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/fd_capture_util.h"
#include "tools/pc_sampler.h"
#include "tools/profile_clock.h"
#include "tools/profile_isr_window.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

namespace hs_test {
namespace profiling_tests {

/** @brief Clock reads defer seam interrupts and restore the incoming mask. */
inline void test_profile_clock_pair_interrupt_mask() {
  for (uint32_t prior_mask : {0u, 1u}) {
    uint32_t mask = prior_mask;
    uint32_t wall = 100, cycles = 60000, interrupts = 0;
    bool pending = false;
    const auto service_interrupt = [&] {
      if (pending && mask == 0) {
        wall += 17;
        cycles += 10200;
        ++interrupts;
        pending = false;
      }
    };
    const auto sample = hs::sample_profile_clocks(
        [&] {
          HS_EXPECT_EQ(mask, 1u);
          const auto value = wall;
          pending = true;
          service_interrupt();
          return value;
        },
        [&] {
          HS_EXPECT_EQ(mask, 1u);
          HS_EXPECT_EQ(interrupts, 0u);
          return cycles;
        },
        [&] {
          const auto saved = mask;
          mask = 1;
          return saved;
        },
        [&](uint32_t saved) {
          mask = saved;
          service_interrupt();
        });
    HS_EXPECT_EQ(sample.wall_us, 100u);
    HS_EXPECT_EQ(sample.cycles, 60000u);
    HS_EXPECT_EQ(mask, prior_mask);
    HS_EXPECT_EQ(interrupts, prior_mask == 0 ? 1u : 0u);
    if (prior_mask == 0) {
      const auto end = hs::sample_profile_clocks(
          [&] { return wall; }, [&] { return cycles; },
          [&] {
            const auto saved = mask;
            mask = 1;
            return saved;
          },
          [&](uint32_t saved) { mask = saved; });
      HS_EXPECT_EQ(end.wall_us - sample.wall_us, 17u);
      HS_EXPECT_EQ(end.cycles - sample.cycles, 10200u);
    }
  }
}

/** @brief Captures the counter report into a NUL-terminated bounded buffer. */
inline bool capture_log_all(char *out, size_t n) {
  const auto text = capture_stdout([] { hs::CycleCounter::log_all(); });
  if (!text || n == 0)
    return false;
  std::snprintf(out, n, "%s", text->c_str());
  return true;
}

// --- u64_dec ----------------------------------------------------------------

/**
 * @brief u64_dec's exact-fit 21-byte buffer between two canary bands.
 * @details UINT64_MAX fills the buffer to its first byte.
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
  HS_EXPECT_STREQ(s, expected);
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
  char expected[21] = "1";
  for (int digits = 1; digits <= 20; ++digits) {
    expect_u64_dec(v, expected);
    if (digits < 20) {
      expected[digits] = '0';
      expected[digits + 1] = '\0';
      v *= 10;
    }
  }
}

// --- registry ---------------------------------------------------------------

/**
 * @brief Verifies find_suffix matches whole names and proper suffixes, prefers
 *        the most recently registered match, and rejects a suffix no name can
 *        contain.
 */
inline void test_find_suffix() {
  hs::CycleCounter alpha("prof_suffix_alpha_tail");
  hs::CycleCounter beta("prof_suffix_beta_tail");
  hs::CycleCounter gamma("prof_suffix_gamma");

  // Registration pushes onto the list head, so the newest match wins.
  HS_EXPECT_EQ(hs::CycleCounter::find_suffix("_tail"), &beta);
  HS_EXPECT_EQ(hs::CycleCounter::find_suffix("prof_suffix_gamma"), &gamma);
  HS_EXPECT_EQ(hs::CycleCounter::find_suffix("alpha_tail"), &alpha);
  HS_EXPECT_EQ(hs::CycleCounter::find_suffix(""), &gamma);
  HS_EXPECT_TRUE(hs::CycleCounter::find_suffix("_prof_no_such_counter") ==
                 nullptr);

  // Longer than any registered name.
  char long_suffix[200];
  std::memset(long_suffix, 'z', sizeof(long_suffix) - 1);
  long_suffix[sizeof(long_suffix) - 1] = '\0';
  HS_EXPECT_TRUE(hs::CycleCounter::find_suffix(long_suffix) == nullptr);
}

/**
 * @brief Verifies a destroyed counter leaves no child pointing at its storage.
 * @details The orphan comes back as a root and re-latches on its next entry.
 */
inline void test_destructor_unlatches_children() {
  hs::CycleCounter child("prof_orphan_child");
  hs::CycleCounter regrown("prof_orphan_regrown");
  hs::CycleCounter::reset_all();
  {
    hs::CycleCounter transient("prof_orphan_transient");
    hs::CycleScope st(transient);
    hs::CycleScope sc(child);
    HS_EXPECT_EQ(child.parent, &transient);
  }
  HS_EXPECT_TRUE(child.parent == nullptr);

  // The next enclosing counter becomes the parent without a mixed-parent flag.
  {
    hs::CycleScope sr(regrown);
    hs::CycleScope sc(child);
  }
  HS_EXPECT_EQ(child.parent, &regrown);
  HS_EXPECT_FALSE(child.mixed_parent);
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
  outer.mixed_parent = true;
  inner.cycles = 99;
  other.cycles = 77;

  hs::CycleCounter::reset_all();

  HS_EXPECT_EQ(outer.cycles, 0u);
  HS_EXPECT_EQ(outer.count, 0u);
  HS_EXPECT_FALSE(outer.mixed_parent);
  HS_EXPECT_EQ(other.cycles, 0u);
  HS_EXPECT_EQ(other.count, 0u);
  HS_EXPECT_FALSE(other.mixed_parent);
  HS_EXPECT_EQ(inner.cycles, 0u);
  HS_EXPECT_EQ(inner.count, 0u);
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

/** @brief Tags active duplicate-name counters and ignores idle namesakes. */
inline void test_duplicate_name_tags_active_counters() {
  static hs::CycleCounter first("prof_duplicate");
  static hs::CycleCounter second("prof_duplicate");
  static hs::CycleCounter idle("prof_duplicate");
  hs::CycleCounter::reset_all();
  {
    hs::CycleScope scope(first);
  }
  {
    hs::CycleScope scope(second);
  }
  HS_EXPECT_EQ(idle.count, 0u);
  char report[4096];
  if (!capture_log_all(report, sizeof(report))) {
    HS_EXPECT(false, "log_all capture set up");
    return;
  }
  const char *tag = std::strstr(report, "DUPLICATE-NAME");
  HS_EXPECT_TRUE(tag != nullptr);
  if (tag != nullptr) {
    tag = std::strstr(tag + 1, "DUPLICATE-NAME");
    HS_EXPECT_TRUE(tag != nullptr);
    if (tag != nullptr)
      HS_EXPECT_TRUE(std::strstr(tag + 1, "DUPLICATE-NAME") == nullptr);
  }
  first.reset();
  if (!capture_log_all(report, sizeof(report))) {
    HS_EXPECT(false, "log_all capture set up");
    return;
  }
  HS_EXPECT_TRUE(std::strstr(report, "prof_duplicate") != nullptr);
  HS_EXPECT_TRUE(std::strstr(report, "DUPLICATE-NAME") == nullptr);
}

/**
 * @brief Verifies parent retirement keeps the child rooted and flags the
 *        descendant when it later runs at the root.
 */
inline void test_parent_retirement_prevents_cycle() {
  hs::CycleCounter child("prof_cycle_child");
  hs::CycleCounter descendant("prof_cycle_descendant");
  {
    hs::CycleCounter parent("prof_cycle_parent");
    hs::CycleScope sp(parent);
    hs::CycleScope sc(child);
    hs::CycleScope sd(descendant);
  }
  {
    hs::CycleScope sd(descendant);
    hs::CycleScope sc(child);
  }
  HS_EXPECT_TRUE(child.parent == nullptr);
  HS_EXPECT_EQ(descendant.parent, &child);
  HS_EXPECT_TRUE(descendant.mixed_parent);
  HS_EXPECT_FALSE(child.mixed_parent);
  child.cycles = 100;
  descendant.cycles = 50;
  char report[4096];
  if (!capture_log_all(report, sizeof(report))) {
    HS_EXPECT(false, "log_all capture set up");
    return;
  }
  HS_EXPECT_TRUE(std::strstr(report, "prof_cycle_child") != nullptr);
  HS_EXPECT_TRUE(std::strstr(report, "prof_cycle_descendant") != nullptr);
}

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
  // Neither edge is rewritten; left stays a root.
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
  static hs::CycleCounter kept("prof_orphan_kept");
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
  HS_EXPECT_TRUE(std::strstr(report, "prof_orphan_kept") != nullptr);
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

/** @brief ISR snapshots count logging work in its matching elapsed window. */
inline void test_isr_window_logging_latency() {
  for (uint32_t latency : {0u, 50u, 300u}) {
    hs::ProfileIsrWindow window;
    hs::IsrCycleStats wake, pack, submit;
    uint32_t now = 100;
    bool masked = false;
    auto clock = [&] {
      HS_EXPECT_TRUE(masked);
      return now;
    };
    auto disable = [&] {
      HS_EXPECT_FALSE(masked);
      masked = true;
      return 7u;
    };
    auto restore = [&](uint32_t mask) {
      HS_EXPECT_EQ(mask, 7u);
      HS_EXPECT_TRUE(masked);
      HS_EXPECT_EQ(wake.count, 0u);
      HS_EXPECT_EQ(pack.count, 0u);
      HS_EXPECT_EQ(submit.count, 0u);
      masked = false;
    };
    auto capture = [&] {
      return window.capture(wake, pack, submit, clock, disable, restore);
    };
    auto advance = [&](uint32_t us) {
      HS_EXPECT_FALSE(masked);
      now += us;
      if (us) {
        wake.add(us * 60u);
        pack.add(us * 30u);
        submit.add(us * 6u);
      }
    };
    wake.add(999);
    capture();
    advance(1000);
    const auto first = capture();
    HS_EXPECT_EQ(first.window_us, 1000u);
    HS_EXPECT_EQ(first.wake.cycles, 60000u);
    HS_EXPECT_EQ(first.pack.cycles, 30000u);
    HS_EXPECT_EQ(first.submit.cycles, 6000u);
    auto log = [&](const hs::IsrCycleStats &stats, uint64_t expected) {
      HS_EXPECT_EQ(stats.cycles, expected);
      advance(latency);
      HS_EXPECT_EQ(stats.cycles, expected);
    };
    log(first.wake, 60000u);
    log(first.pack, 30000u);
    log(first.submit, 6000u);
    advance(1000);
    const auto second = capture();
    const uint32_t EXPECTED_US = 1000u + 3u * latency;
    HS_EXPECT_EQ(second.window_us, EXPECTED_US);
    HS_EXPECT_EQ(second.wake.cycles, uint64_t{EXPECTED_US} * 60u);
    HS_EXPECT_EQ(second.pack.cycles, uint64_t{EXPECTED_US} * 30u);
    HS_EXPECT_EQ(second.submit.cycles, uint64_t{EXPECTED_US} * 6u);
    now = UINT32_MAX - 99u;
    capture();
    advance(200);
    const auto wrapped = capture();
    HS_EXPECT_EQ(wrapped.window_us, 200u);
    HS_EXPECT_EQ(wrapped.wake.cycles, 12000u);
  }
}

// --- PC sampler -------------------------------------------------------------

/** @brief EXC_RETURN for Thread mode on MSP with a basic frame. */
inline constexpr uint32_t PC_EXC_THREAD_MSP = 0xFFFFFFF9u;
/** @brief EXC_RETURN for Handler mode on MSP with a basic frame. */
inline constexpr uint32_t PC_EXC_HANDLER = 0xFFFFFFF1u;

/**
 * @brief Verifies the stacked return address is word 6 of the stack EXC_RETURN
 *        selects, for basic and FP-extended frames alike.
 */
inline void test_pc_interrupted_pc_selects_stack() {
  uint32_t msp[26], psp[26];
  for (size_t i = 0; i < 26; ++i) {
    msp[i] = 0xAA000000u + (uint32_t)i;
    psp[i] = 0xBB000000u + (uint32_t)i;
  }
  msp[6] = 0x00001234u;
  psp[6] = 0x60005678u;
  HS_EXPECT_EQ(hs::interrupted_pc(msp, psp, 0xFFFFFFF9u), 0x00001234u);
  HS_EXPECT_EQ(hs::interrupted_pc(msp, psp, 0xFFFFFFE9u), 0x00001234u);
  HS_EXPECT_EQ(hs::interrupted_pc(msp, psp, 0xFFFFFFF1u), 0x00001234u);
  HS_EXPECT_EQ(hs::interrupted_pc(msp, psp, 0xFFFFFFE1u), 0x00001234u);
  HS_EXPECT_EQ(hs::interrupted_pc(msp, psp, 0xFFFFFFFDu), 0x60005678u);
  HS_EXPECT_EQ(hs::interrupted_pc(msp, psp, 0xFFFFFFEDu), 0x60005678u);
}

/**
 * @brief Verifies lost periods round the elapsed time to whole periods and
 *        exclude the period that fired.
 */
inline void test_pc_lost_ticks_rounds_to_periods() {
  HS_EXPECT_EQ(hs::pc_lost_ticks(0, 60000), 0u);
  HS_EXPECT_EQ(hs::pc_lost_ticks(60000, 60000), 0u);
  HS_EXPECT_EQ(hs::pc_lost_ticks(89999, 60000), 0u);
  HS_EXPECT_EQ(hs::pc_lost_ticks(90000, 60000), 1u);
  HS_EXPECT_EQ(hs::pc_lost_ticks(149999, 60000), 1u);
  HS_EXPECT_EQ(hs::pc_lost_ticks(150000, 60000), 2u);
  HS_EXPECT_EQ(hs::pc_lost_ticks(7, 5), 0u);
  HS_EXPECT_EQ(hs::pc_lost_ticks(8, 5), 1u);
  HS_EXPECT_EQ(hs::pc_lost_ticks(UINT32_MAX, 60000), 71582u);
}

/** @brief Verifies latency buckets are log2 widths capped at the last one. */
inline void test_pc_latency_buckets() {
  HS_EXPECT_EQ(hs::pc_latency_bucket(0), 0u);
  HS_EXPECT_EQ(hs::pc_latency_bucket(1), 1u);
  HS_EXPECT_EQ(hs::pc_latency_bucket(2), 2u);
  HS_EXPECT_EQ(hs::pc_latency_bucket(3), 2u);
  HS_EXPECT_EQ(hs::pc_latency_bucket(4), 3u);
  HS_EXPECT_EQ(hs::pc_latency_bucket(1023), 10u);
  HS_EXPECT_EQ(hs::pc_latency_bucket(1024), 11u);
  HS_EXPECT_EQ(hs::pc_latency_bucket(UINT32_MAX), hs::PC_LATENCY_BUCKETS - 1);
}

/**
 * @brief Collects a histogram's bins as start address to count.
 * @param h Histogram to read.
 * @return The occupied bins.
 */
template <typename H> std::map<uint32_t, uint32_t> pc_bins_of(const H &h) {
  std::map<uint32_t, uint32_t> bins;
  h.for_each([&](uint32_t address, uint32_t count) {
    HS_EXPECT_TRUE(bins.emplace(address, count).second);
  });
  return bins;
}

/**
 * @brief Verifies PCs aggregate per aligned bin, including the bin at address
 *        zero, and that a full table rejects only new bins.
 */
inline void test_pc_histogram_bins_and_capacity() {
  hs::PcHistogram<8, 4> h;
  HS_EXPECT_TRUE(h.add(0x00001000u));
  HS_EXPECT_TRUE(h.add(0x0000100Eu));
  HS_EXPECT_TRUE(h.add(0x00001010u));
  HS_EXPECT_TRUE(h.add(0x60008766u));
  HS_EXPECT_TRUE(h.add(0x00000002u));
  HS_EXPECT_EQ(h.size(), 4u);
  const std::map<uint32_t, uint32_t> expected{{0x00000000u, 1u},
                                              {0x00001000u, 2u},
                                              {0x00001010u, 1u},
                                              {0x60008760u, 1u}};
  HS_EXPECT_TRUE(pc_bins_of(h) == expected);

  hs::PcHistogram<2, 5> full;
  HS_EXPECT_TRUE(full.add(0x100u));
  HS_EXPECT_TRUE(full.add(0x200u));
  HS_EXPECT_FALSE(full.add(0x300u));
  HS_EXPECT_TRUE(full.add(0x11Fu));
  HS_EXPECT_EQ(full.size(), 2u);
  const std::map<uint32_t, uint32_t> kept{{0x100u, 2u}, {0x200u, 1u}};
  HS_EXPECT_TRUE(pc_bins_of(full) == kept);
}

/**
 * @brief Verifies each tick lands in exactly one of idle, deferred, dropped or
 *        recorded, with lost periods charged to the gate's state.
 */
inline void test_pc_sampler_classifies_ticks() {
  hs::PcSampler<2, 4> s;
  s.configure({10000, 48, 24});
  s.tick(false, 3, 2, PC_EXC_THREAD_MSP, 0x100u);
  s.tick(true, 48, 0, PC_EXC_THREAD_MSP, 0x200u);
  s.tick(true, 49, 0, PC_EXC_THREAD_MSP, 0x300u);
  s.tick(true, 0, 3, PC_EXC_THREAD_MSP, 0x300u);
  s.tick(true, 0, 0, PC_EXC_HANDLER, 0x400u);
  s.tick(true, 1, 0, PC_EXC_THREAD_MSP, 0x500u);
  const hs::PcSampleCounters &c = s.totals();
  HS_EXPECT_EQ(c.ticks, 6u);
  HS_EXPECT_EQ(c.idle, 3u);
  HS_EXPECT_EQ(c.deferred, 2u);
  HS_EXPECT_EQ(c.missed, 3u);
  HS_EXPECT_EQ(c.recorded, 2u);
  HS_EXPECT_EQ(c.handler, 1u);
  HS_EXPECT_EQ(c.dropped, 1u);
  HS_EXPECT_EQ(c.latency[0], 2u);
  HS_EXPECT_EQ(c.latency[1], 1u);
  HS_EXPECT_EQ(c.latency[2], 1u);
  HS_EXPECT_EQ(c.latency[6], 2u);
  const std::map<uint32_t, uint32_t> expected{{0x200u, 1u}, {0x400u, 1u}};
  HS_EXPECT_TRUE(pc_bins_of(s.bins()) == expected);
}

/** @brief Verifies the dump's exact line format, which pc_profile.py parses. */
inline void test_pc_sampler_dump_format() {
  hs::PcSampler<4, 4> s;
  s.configure({10000, 48, 24});
  s.tick(true, 0, 0, PC_EXC_THREAD_MSP, 0x00001234u);
  s.tick(true, 2, 0, PC_EXC_THREAD_MSP, 0x00001238u);
  s.tick(true, 5, 0, PC_EXC_THREAD_MSP, 0x6000ABCDu);
  s.tick(false, 1, 0, PC_EXC_THREAD_MSP, 0x00001234u);
  s.note_service(100);
  s.note_service(260);
  s.note_service(90);
  s.note_service(50);
  std::vector<std::string> lines;
  s.dump([&](const char *line) { lines.emplace_back(line); }, 1234);
  HS_EXPECT_SIZE_OR_RETURN(lines, 6u);
  HS_EXPECT_STREQ(lines[0].c_str(),
                  "pc profile: version=1 hz=10000 bin_bytes=16 capacity=4 "
                  "defer_counts=48 counts_per_us=24 elapsed_ms=1234");
  HS_EXPECT_STREQ(lines[1].c_str(),
                  "pc profile: ticks=4 idle=1 deferred=0 missed=0 handler=0 "
                  "recorded=3 dropped=0 bins=2 service_avg=125 "
                  "service_max=260");
  HS_EXPECT_STREQ(lines[2].c_str(),
                  "pc profile latency: 1 1 1 1 0 0 0 0 0 0 0 0");
  const bool ordered = lines[3] < lines[4];
  HS_EXPECT_STREQ((ordered ? lines[3] : lines[4]).c_str(), "pc 00001230 2");
  HS_EXPECT_STREQ((ordered ? lines[4] : lines[3]).c_str(), "pc 6000abc0 1");
  HS_EXPECT_STREQ(lines[5].c_str(), "pc profile end");
}

/** @brief Verifies PcSampleGate sets the gate and restores it on nesting. */
inline void test_pc_sample_gate_nests() {
  HS_EXPECT_EQ(static_cast<unsigned>(hs::pc_sample_gate()), 0u);
  {
    hs::PcSampleGate render(true);
    HS_EXPECT_EQ(static_cast<unsigned>(hs::pc_sample_gate()), 1u);
    {
      hs::PcSampleGate wait(false);
      HS_EXPECT_EQ(static_cast<unsigned>(hs::pc_sample_gate()), 0u);
    }
    HS_EXPECT_EQ(static_cast<unsigned>(hs::pc_sample_gate()), 1u);
  }
  HS_EXPECT_EQ(static_cast<unsigned>(hs::pc_sample_gate()), 0u);
}

/** @brief Sampled intervals exclude scope setup and retain child nesting. */
inline void test_sampled_cycle_scope_bounds_and_nesting() {
  hs::CycleCounter outer("sampled_outer");
  hs::CycleCounter sampled("sampled");
  hs::CycleCounter child("sampled_child");
  hs::CycleCounter after("sampled_after");
  hs::CycleCounter unrecorded("sampled_unrecorded");
  {
    hs::CycleScope enclosing(outer);
    {
      hs::SampledCycleScope scope(sampled);
      HS_EXPECT_EQ(sampled.count, 0u);
      HS_EXPECT_EQ(sampled.cycles, 0u);
      {
        hs::CycleScope nested(child);
      }
      HS_EXPECT_TRUE(child.parent == &sampled);
      scope.record(0xfffffff0u, 0x30u);
      scope.record(100u, 700u);
      HS_EXPECT_EQ(sampled.cycles, 664u);
      HS_EXPECT_EQ(sampled.count, 2u);
    }
    HS_EXPECT_TRUE(sampled.parent == &outer);
    HS_EXPECT_EQ(sampled.cycles, 664u);
    HS_EXPECT_EQ(sampled.count, 2u);
    {
      hs::SampledCycleScope scope(unrecorded);
    }
    HS_EXPECT_EQ(unrecorded.count, 0u);
    HS_EXPECT_EQ(unrecorded.cycles, 0u);
    {
      hs::CycleScope sibling(after);
    }
    HS_EXPECT_TRUE(after.parent == &outer);
  }
  HS_EXPECT_EQ(outer.count, 1u);
  HS_EXPECT_EQ(child.count, 1u);
  HS_EXPECT_EQ(after.count, 1u);
}

/**
 * @brief Runs the profiling test cases.
 * @return The module's failure count.
 */
inline int run_profiling_tests() {
  hs_test::ModuleFixture fixture("profiling");

  HS_RUN(test_sampled_cycle_scope_bounds_and_nesting);
  HS_RUN(test_profile_clock_pair_interrupt_mask);
  HS_RUN(test_u64_dec_boundaries);
  HS_RUN(test_u64_dec_powers_of_ten);
  HS_RUN(test_find_suffix);
  HS_RUN(test_destructor_unlatches_children);
  HS_RUN(test_reset_all_clears_counts);
  HS_RUN(test_nesting_latches_parent);
  HS_RUN(test_recursive_scope_does_not_self_parent);
  HS_RUN(test_second_caller_flags_mixed_parent);
  HS_RUN(test_duplicate_name_tags_active_counters);
  HS_RUN(test_mutual_nesting_keeps_a_root);
  HS_RUN(test_parent_retirement_prevents_cycle);
  HS_RUN(test_log_all_reports_tree);
  HS_RUN(test_reset_does_not_orphan_subtree);
  HS_RUN(test_isr_cycle_stats);
  HS_RUN(test_isr_window_logging_latency);
  HS_RUN(test_pc_interrupted_pc_selects_stack);
  HS_RUN(test_pc_lost_ticks_rounds_to_periods);
  HS_RUN(test_pc_latency_buckets);
  HS_RUN(test_pc_histogram_bins_and_capacity);
  HS_RUN(test_pc_sampler_classifies_ticks);
  HS_RUN(test_pc_sampler_dump_format);
  HS_RUN(test_pc_sample_gate_nests);

  return fixture.result();
}

} // namespace profiling_tests
} // namespace hs_test
