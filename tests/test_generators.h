/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/generators.h — the generate() wrapper that scopes the two
 * global scratch arenas around a procedural-geometry callback. Verifies the
 * deterministic arena lifecycle contract:
 *   - fn receives (target, scratch_arena_a, scratch_arena_b, args...);
 *   - both scratch arenas are reset to offset 0 BEFORE fn runs;
 *   - scratch allocations made by fn are rolled back after generate() returns;
 *   - the target arena is left untouched (its allocations persist);
 *   - extra args and the return value are forwarded;
 *   - reentrancy: a nested generate() stacks above the caller's live frame
 *     rather than resetting it, across both a single nest and a deep stack.
 */
#pragma once

#include "core/generators.h"
#include "core/memory.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace generators_tests {

inline uint8_t gen_target_buf[8 * 1024];

/**
 * @brief Verifies the full generate() contract in one pass.
 * @details Asserts the scratch arenas are reset before fn runs, fn receives the
 * two global scratch arenas plus the caller's target, scratch allocations roll
 * back on return, target allocations persist, and the extra arg and return
 * value are forwarded.
 */
inline void test_generate_lifecycle_and_forwarding() {
  Arena target(gen_target_buf, sizeof(gen_target_buf));

  scratch_arena_a.reset();
  scratch_arena_a.allocate(123);
  scratch_arena_b.reset();
  scratch_arena_b.allocate(456);
  HS_EXPECT_GT(scratch_arena_a.get_offset(), (size_t)0);
  HS_EXPECT_GT(scratch_arena_b.get_offset(), (size_t)0);

  Arena *seen_target = nullptr, *seen_a = nullptr, *seen_b = nullptr;
  size_t a_offset_in_fn = 999, b_offset_in_fn = 999;
  int captured_arg = 0;

  int result = generate(
      target,
      [&](Arena &t, Arena &a, Arena &b, int arg) {
        seen_target = &t;
        seen_a = &a;
        seen_b = &b;
        a_offset_in_fn = a.get_offset();
        b_offset_in_fn = b.get_offset();
        captured_arg = arg;
        a.allocate(64);
        t.allocate(32);
        return 7;
      },
      42);

  HS_EXPECT_TRUE(seen_a == &scratch_arena_a);
  HS_EXPECT_TRUE(seen_b == &scratch_arena_b);
  HS_EXPECT_TRUE(seen_target == &target);

  HS_EXPECT_EQ(a_offset_in_fn, (size_t)0);
  HS_EXPECT_EQ(b_offset_in_fn, (size_t)0);

  HS_EXPECT_EQ(captured_arg, 42);
  HS_EXPECT_EQ(result, 7);

  HS_EXPECT_EQ(scratch_arena_a.get_offset(), (size_t)0);
  HS_EXPECT_EQ(scratch_arena_b.get_offset(), (size_t)0);

  HS_EXPECT_EQ(target.get_offset(), (size_t)32);
}

/**
 * @brief Confirms target allocations persist while scratch is rolled back.
 * @details A generator can write into target while using scratch; only scratch
 * is rolled back. Verifies that two sequential generate() calls accumulate
 * their target allocations.
 */
inline void test_generate_nested_target_persists() {
  Arena target(gen_target_buf, sizeof(gen_target_buf));

  (void)generate(target, [](Arena &t, Arena &a, Arena &, int n) {
    a.allocate(100);
    for (int i = 0; i < n; ++i)
      t.allocate(16);
    return 0;
  }, 3);
  HS_EXPECT_EQ(target.get_offset(), (size_t)48);

  (void)generate(target, [](Arena &t, Arena &, Arena &b, int n) {
    b.allocate(200);
    for (int i = 0; i < n; ++i)
      t.allocate(16);
    return 0;
  }, 2);
  HS_EXPECT_EQ(target.get_offset(), (size_t)80);

  HS_EXPECT_EQ(scratch_arena_a.get_offset(), (size_t)0);
  HS_EXPECT_EQ(scratch_arena_b.get_offset(), (size_t)0);
}

/**
 * @brief Verifies reentrant generate() does not clobber the outer frame.
 * @details generate() is reentrant: a callback may call generate() again. The
 * inner call must stack its scratch above the outer frame's live allocations
 * rather than resetting the arena out from under it. Checks that the outer
 * sentinel survives the nested call.
 */
inline void test_generate_reentrant_nesting_does_not_clobber() {
  Arena target(gen_target_buf, sizeof(gen_target_buf));

  size_t outer_offset_after_inner = 999;
  size_t inner_start_offset = 999;
  uint8_t outer_value_after_inner = 0;

  (void)generate(target, [&](Arena &t, Arena &a, Arena &, int) {
    uint8_t *outer = static_cast<uint8_t *>(a.allocate(64));
    outer[0] = 0xAB;
    const size_t outer_offset_before_inner = a.get_offset();

    (void)generate(t, [&](Arena &, Arena &ia, Arena &, int) {
      inner_start_offset = ia.get_offset();
      uint8_t *inner = static_cast<uint8_t *>(ia.allocate(32));
      inner[0] = 0xCD; // would corrupt outer[0] if the arena had been reset
      return 0;
    }, 0);

    outer_value_after_inner = outer[0];
    outer_offset_after_inner = a.get_offset();
    HS_EXPECT_EQ(a.get_offset(), outer_offset_before_inner);
    return 0;
  }, 0);

  HS_EXPECT_EQ(inner_start_offset, (size_t)64);
  HS_EXPECT_EQ((int)outer_value_after_inner, 0xAB);
  HS_EXPECT_EQ(outer_offset_after_inner, (size_t)64);

  HS_EXPECT_EQ(scratch_arena_a.get_offset(), (size_t)0);
  HS_EXPECT_EQ(scratch_arena_b.get_offset(), (size_t)0);
}

// --- Deep (multi-level) reentrant nesting -----------------------------------

inline constexpr int kDeepLevels = 5;
inline size_t g_deep_a_entry[kDeepLevels];

/**
 * @brief Recursive generate() body for the deep-nesting stress.
 * @details Each level claims distinct, level-sized scratch in BOTH arenas,
 * stamps a per-level sentinel into the first and last byte of each block,
 * recurses through generate() (which must stack above this frame, never reset
 * it), then on the way back asserts its own offsets and sentinels are exactly as
 * it left them — proving no deeper frame rewound or overwrote a shallower one.
 */
inline int gen_deep_level(Arena &t, Arena &a, Arena &b, int level,
                          int max_level) {
  const size_t a_sz = 16 * static_cast<size_t>(level + 1);
  const size_t b_sz = 8 * static_cast<size_t>(level + 1);
  g_deep_a_entry[level] = a.get_offset();

  uint8_t *am = static_cast<uint8_t *>(a.allocate(a_sz));
  uint8_t *bm = static_cast<uint8_t *>(b.allocate(b_sz));
  const uint8_t sentinel = static_cast<uint8_t>(0xA0 + level);
  am[0] = am[a_sz - 1] = sentinel;
  bm[0] = bm[b_sz - 1] = sentinel;
  const size_t a_after = a.get_offset();
  const size_t b_after = b.get_offset();

  if (level + 1 < max_level)
    (void)generate(t, gen_deep_level, level + 1, max_level);

  HS_EXPECT_EQ(a.get_offset(), a_after);
  HS_EXPECT_EQ(b.get_offset(), b_after);
  HS_EXPECT_EQ((int)am[0], (int)sentinel);
  HS_EXPECT_EQ((int)am[a_sz - 1], (int)sentinel);
  HS_EXPECT_EQ((int)bm[0], (int)sentinel);
  HS_EXPECT_EQ((int)bm[b_sz - 1], (int)sentinel);
  return level;
}

/**
 * @brief Stress-tests the reentrant scratch protocol across many nested levels.
 * @details The single-nest test proves one inner call does not clobber its
 * caller; this drives kDeepLevels stacked frames and verifies each one stacked
 * strictly above the previous level's live high-water (the reset fires only at
 * the outermost call) and that the outermost scope rolled both arenas back to
 * empty. Per-level sentinel survival is asserted inside gen_deep_level.
 */
inline void test_generate_deep_nesting_stacks_and_unwinds() {
  Arena target(gen_target_buf, sizeof(gen_target_buf));
  for (int i = 0; i < kDeepLevels; ++i)
    g_deep_a_entry[i] = 999;

  (void)generate(target, gen_deep_level, 0, kDeepLevels);

  HS_EXPECT_EQ(g_deep_a_entry[0], (size_t)0);
  for (int level = 1; level < kDeepLevels; ++level)
    HS_EXPECT_GT(g_deep_a_entry[level], g_deep_a_entry[level - 1]);

  HS_EXPECT_EQ(scratch_arena_a.get_offset(), (size_t)0);
  HS_EXPECT_EQ(scratch_arena_b.get_offset(), (size_t)0);
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every generators test case.
 * @return The module's failure count.
 */
inline int run_generators_tests() {
  hs_test::ModuleFixture fixture("generators");

  test_generate_lifecycle_and_forwarding();
  test_generate_nested_target_persists();
  test_generate_reentrant_nesting_does_not_clobber();
  test_generate_deep_nesting_stacks_and_unwinds();

  return fixture.result();
}

} // namespace generators_tests
} // namespace hs_test

