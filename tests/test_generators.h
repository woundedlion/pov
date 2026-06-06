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
 *   - extra args and the return value are forwarded.
 */
#pragma once

#include "core/generators.h"
#include "core/memory.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace generators_tests {

inline uint8_t gen_target_buf[8 * 1024];

inline void test_generate_lifecycle_and_forwarding() {
  Arena target(gen_target_buf, sizeof(gen_target_buf));

  // Pre-dirty the global scratch arenas so we can prove generate() resets them.
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
        a_offset_in_fn = a.get_offset(); // should be 0 (reset before fn)
        b_offset_in_fn = b.get_offset();
        captured_arg = arg;
        a.allocate(64); // scratch allocation — must be rolled back
        t.allocate(32); // target allocation — must persist
        return 7;
      },
      42);

  // fn was handed the two global scratch arenas and our target.
  HS_EXPECT_TRUE(seen_a == &scratch_arena_a);
  HS_EXPECT_TRUE(seen_b == &scratch_arena_b);
  HS_EXPECT_TRUE(seen_target == &target);

  // Scratch arenas were reset to empty before fn ran.
  HS_EXPECT_EQ(a_offset_in_fn, (size_t)0);
  HS_EXPECT_EQ(b_offset_in_fn, (size_t)0);

  // Arg + return value forwarded.
  HS_EXPECT_EQ(captured_arg, 42);
  HS_EXPECT_EQ(result, 7);

  // ScratchScope rolled the scratch arenas back to empty on return.
  HS_EXPECT_EQ(scratch_arena_a.get_offset(), (size_t)0);
  HS_EXPECT_EQ(scratch_arena_b.get_offset(), (size_t)0);

  // The target arena is owned by the caller — generate() must not touch it.
  HS_EXPECT_EQ(target.get_offset(), (size_t)32);
}

inline void test_generate_nested_target_persists() {
  // A generator can itself write into target while using scratch; only scratch
  // is rolled back. Two sequential generate() calls accumulate in target.
  Arena target(gen_target_buf, sizeof(gen_target_buf));

  generate(target, [](Arena &t, Arena &a, Arena &, int n) {
    a.allocate(100);
    for (int i = 0; i < n; ++i)
      t.allocate(16);
    return 0;
  }, 3);
  HS_EXPECT_EQ(target.get_offset(), (size_t)48);

  generate(target, [](Arena &t, Arena &, Arena &b, int n) {
    b.allocate(200);
    for (int i = 0; i < n; ++i)
      t.allocate(16);
    return 0;
  }, 2);
  HS_EXPECT_EQ(target.get_offset(), (size_t)80);

  // Scratch fully rolled back after both calls.
  HS_EXPECT_EQ(scratch_arena_a.get_offset(), (size_t)0);
  HS_EXPECT_EQ(scratch_arena_b.get_offset(), (size_t)0);
}

// ============================================================================
// Runner
// ============================================================================

inline int run_generators_tests() {
  auto scope = hs_test::begin_module("generators");

  test_generate_lifecycle_and_forwarding();
  test_generate_nested_target_persists();

  return hs_test::end_module(scope);
}

} // namespace generators_tests
} // namespace hs_test

