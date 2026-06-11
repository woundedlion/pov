/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/memory.h — Arena, TriangularBitset, ArenaVector,
 * ArenaSpan, ScratchScope, and Persist<T>.
 *
 * NOTE: The always-on HS_CHECK traps are driven (in forked child processes) by
 * the death harness in tests/test_death.h — case_arena_oom (allocate past
 * capacity), case_arena_set_offset_overflow (rewind past capacity),
 * case_arena_vector_overflow / append_bulk_overflow (fixed-capacity push), and
 * case_arena_oversubscribed.
 *
 * The ArenaVector lifetime guards (unbound access, use-after-free) are debug
 * asserts, not HS_CHECK traps — they compile out under NDEBUG by design, so the
 * per-access check stays free in the device/WASM build. The death harness
 * therefore can't cover them; instead the bookkeeping they key on is exercised
 * in-process here: the bound_ flag through the construct/bind/move lifecycle
 * (test_arenavec_default_unbound / _bind / _move_construct / _move_assign), and
 * the use-after-free generation snapshot via test_arenavec_stale_binding_after_reset
 * (an arena reset marks a live binding stale; rebinding it is then a debug
 * contract trap). Double-bind is NOT a precondition violation at all: a
 * re-bind that grows is a supported pattern (it abandons the old block until the
 * next arena reset — see ArenaVector::bind), covered by test_arenavec_rebind_grows.
 */
#pragma once

#include <cstdint>
#include <cstring>
#include <utility>
#include "core/memory.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace mem {

// Backing storage for arenas used in these tests.
// Sized generously so OOM tests can be exercised explicitly.
inline uint8_t test_buf_a[64 * 1024];
inline uint8_t test_buf_b[16 * 1024];

// ============================================================================
// Arena — construction, allocation, alignment
// ============================================================================

inline void test_arena_construction() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  HS_EXPECT_EQ(a.get_capacity(), sizeof(test_buf_a));
  HS_EXPECT_EQ(a.get_offset(), (size_t)0);
  HS_EXPECT_EQ(a.get_high_water_mark(), (size_t)0);
}

inline void test_arena_basic_allocation() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  void *p = a.allocate(64);
  HS_EXPECT_TRUE(p != nullptr);
  HS_EXPECT_TRUE(a.get_offset() >= 64);
  HS_EXPECT_TRUE(a.get_high_water_mark() == a.get_offset());

  void *q = a.allocate(128);
  HS_EXPECT_TRUE(q != nullptr);
  HS_EXPECT_TRUE(q != p);
  HS_EXPECT_TRUE(reinterpret_cast<uintptr_t>(q) >
                 reinterpret_cast<uintptr_t>(p));
}

inline void test_arena_alignment() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  // Eat one byte to force a non-aligned offset
  a.allocate(1, 1);
  void *p16 = a.allocate(32, 16);
  HS_EXPECT_TRUE(p16 != nullptr);
  HS_EXPECT_EQ(reinterpret_cast<uintptr_t>(p16) % 16, (uintptr_t)0);

  void *p32 = a.allocate(32, 32);
  HS_EXPECT_TRUE(p32 != nullptr);
  HS_EXPECT_EQ(reinterpret_cast<uintptr_t>(p32) % 32, (uintptr_t)0);
}

inline void test_arena_high_water_mark() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  a.allocate(100);
  size_t hwm1 = a.get_high_water_mark();
  HS_EXPECT_TRUE(hwm1 >= 100);

  // Rewinding via set_offset does NOT lower the high water mark
  a.set_offset(0);
  HS_EXPECT_EQ(a.get_high_water_mark(), hwm1);

  // A smaller allocation after rewind doesn't change HWM
  a.allocate(50);
  HS_EXPECT_EQ(a.get_high_water_mark(), hwm1);

  // A larger allocation does
  a.set_offset(0);
  a.allocate(200);
  HS_EXPECT_TRUE(a.get_high_water_mark() >= 200);
}

inline void test_arena_reset() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  a.allocate(100);
  size_t hwm = a.get_high_water_mark();
  a.reset();
  HS_EXPECT_EQ(a.get_offset(), (size_t)0);
  // reset() does NOT clear the high water mark (per implementation)
  HS_EXPECT_EQ(a.get_high_water_mark(), hwm);
}

inline void test_arena_set_offset() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  a.allocate(256);
  size_t saved = a.get_offset();
  a.allocate(128);
  HS_EXPECT_TRUE(a.get_offset() > saved);
  a.set_offset(saved);
  HS_EXPECT_EQ(a.get_offset(), saved);
}

inline void test_arena_fills_to_capacity() {
  // Fail-fast contract: over-allocation is an invariant violation that
  // HS_CHECK-traps inside allocate() (it no longer returns nullptr for callers
  // to deref blindly). Like the other precondition traps in memory.h, that is
  // exercised on the bench, not here — we verify only the legal boundary: an
  // allocation that fills exactly to capacity succeeds with in-range pointers.
  uint8_t tiny[64];
  Arena a(tiny, sizeof(tiny));
  void *p = a.allocate(32, 1);
  HS_EXPECT_TRUE(p != nullptr);
  void *q = a.allocate(32, 1); // fills exactly to capacity
  HS_EXPECT_TRUE(q != nullptr);
  HS_EXPECT_EQ(a.get_offset(), sizeof(tiny));

  uintptr_t base = reinterpret_cast<uintptr_t>(tiny);
  HS_EXPECT_EQ(reinterpret_cast<uintptr_t>(p), base);
  HS_EXPECT_EQ(reinterpret_cast<uintptr_t>(q), base + 32);
}

inline void test_arena_rebind() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  a.allocate(100);
  a.rebind(test_buf_b, sizeof(test_buf_b));
  HS_EXPECT_EQ(a.get_capacity(), sizeof(test_buf_b));
  HS_EXPECT_EQ(a.get_offset(), (size_t)0);
  HS_EXPECT_EQ(a.get_high_water_mark(), (size_t)0);

  void *p = a.allocate(64);
  HS_EXPECT_TRUE(p != nullptr);
  uintptr_t addr = reinterpret_cast<uintptr_t>(p);
  uintptr_t base = reinterpret_cast<uintptr_t>(test_buf_b);
  HS_EXPECT_TRUE(addr >= base && addr < base + sizeof(test_buf_b));
}

inline void test_arena_reset_high_water_mark() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  a.allocate(500);
  a.set_offset(100);
  a.reset_high_water_mark();
  HS_EXPECT_EQ(a.get_high_water_mark(), (size_t)100);
}

#ifndef NDEBUG
inline void test_arena_generation_bumps() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  uint32_t g0 = a.get_generation();
  a.reset();
  uint32_t g1 = a.get_generation();
  HS_EXPECT_TRUE(g1 > g0);

  a.rebind(test_buf_b, sizeof(test_buf_b));
  uint32_t g2 = a.get_generation();
  HS_EXPECT_TRUE(g2 > g1);
}
#endif

// ============================================================================
// TriangularBitset
// ============================================================================

inline void test_tribitset_sizes() {
  // MAX_V=10 → 10*9/2 = 45 bits → ceil(45/8) = 6 bytes
  HS_EXPECT_EQ(TriangularBitset<10>::BITS, 45);
  HS_EXPECT_EQ(TriangularBitset<10>::BYTES, 6);
  // MAX_V=16 → 16*15/2 = 120 bits = 15 bytes
  HS_EXPECT_EQ(TriangularBitset<16>::BITS, 120);
  HS_EXPECT_EQ(TriangularBitset<16>::BYTES, 15);
}

inline void test_tribitset_clear_and_test() {
  TriangularBitset<8> bs;
  bs.clear();
  for (int a = 0; a < 8; ++a)
    for (int b = a + 1; b < 8; ++b)
      HS_EXPECT_FALSE(bs.test(a, b));
}

inline void test_tribitset_test_and_set() {
  TriangularBitset<10> bs;
  bs.clear();
  // First insertion of (2, 5) returns false (newly inserted)
  HS_EXPECT_FALSE(bs.test_and_set(2, 5));
  // Second insertion of (2, 5) returns true (already present)
  HS_EXPECT_TRUE(bs.test_and_set(2, 5));
  // test() agrees
  HS_EXPECT_TRUE(bs.test(2, 5));
  // Other pairs untouched
  HS_EXPECT_FALSE(bs.test(2, 6));
  HS_EXPECT_FALSE(bs.test(3, 5));
}

inline void test_tribitset_index_uniqueness() {
  // For MAX_V=8, the 28 distinct pairs must map to 28 distinct indices in [0, 28)
  constexpr int N = 8;
  int seen[64] = {};
  int unique = 0;
  for (int a = 0; a < N; ++a) {
    for (int b = a + 1; b < N; ++b) {
      int idx = TriangularBitset<N>::index(a, b);
      HS_EXPECT_TRUE(idx >= 0 && idx < TriangularBitset<N>::BITS);
      HS_EXPECT_FALSE(seen[idx]);
      seen[idx] = 1;
      ++unique;
    }
  }
  HS_EXPECT_EQ(unique, TriangularBitset<N>::BITS);
}

inline void test_tribitset_all_pairs_independent() {
  TriangularBitset<8> bs;
  bs.clear();
  // Set every other pair (i, i+2) and verify intermediate pairs are unset.
  for (int i = 0; i + 2 < 8; ++i) bs.test_and_set(i, i + 2);
  for (int i = 0; i + 2 < 8; ++i) HS_EXPECT_TRUE(bs.test(i, i + 2));
  for (int i = 0; i + 1 < 8; ++i) HS_EXPECT_FALSE(bs.test(i, i + 1));
}

// ============================================================================
// ArenaVector
// ============================================================================

inline void test_arenavec_default_unbound() {
  ArenaVector<int> v;
  HS_EXPECT_FALSE(v.is_bound());
  HS_EXPECT_EQ(v.size(), (size_t)0);
  HS_EXPECT_EQ(v.capacity(), (size_t)0);
  HS_EXPECT_TRUE(v.empty());
}

inline void test_arenavec_bind() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  ArenaVector<int> v;
  v.bind(a, 32);
  HS_EXPECT_TRUE(v.is_bound());
  HS_EXPECT_EQ(v.capacity(), (size_t)32);
  HS_EXPECT_EQ(v.size(), (size_t)0);
}

inline void test_arenavec_constructor_with_arena() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  ArenaVector<int> v(a, 16);
  HS_EXPECT_TRUE(v.is_bound());
  HS_EXPECT_EQ(v.capacity(), (size_t)16);
}

inline void test_arenavec_push_back_indexing() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  ArenaVector<int> v(a, 8);
  for (int i = 0; i < 8; ++i) v.push_back(i * 10);
  HS_EXPECT_EQ(v.size(), (size_t)8);
  for (int i = 0; i < 8; ++i) HS_EXPECT_EQ(v[i], i * 10);
}

inline void test_arenavec_back() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  ArenaVector<int> v(a, 4);
  v.push_back(1);
  HS_EXPECT_EQ(v.back(), 1);
  v.push_back(2);
  HS_EXPECT_EQ(v.back(), 2);
  v.push_back(3);
  HS_EXPECT_EQ(v.back(), 3);
}

inline void test_arenavec_emplace_back() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  struct Pair { int x; int y; Pair(int xx, int yy) : x(xx), y(yy) {} };
  ArenaVector<Pair> v(a, 4);
  Pair &ref = v.emplace_back(3, 4);
  HS_EXPECT_EQ(ref.x, 3);
  HS_EXPECT_EQ(ref.y, 4);
  HS_EXPECT_EQ(&ref, &v[0]);
  HS_EXPECT_EQ(v.size(), (size_t)1);
}

inline void test_arenavec_append_bulk() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  ArenaVector<int> v(a, 16);
  int src[5] = {10, 20, 30, 40, 50};
  v.append_bulk(src, 5);
  HS_EXPECT_EQ(v.size(), (size_t)5);
  for (int i = 0; i < 5; ++i) HS_EXPECT_EQ(v[i], src[i]);

  // Second bulk append continues
  int more[3] = {60, 70, 80};
  v.append_bulk(more, 3);
  HS_EXPECT_EQ(v.size(), (size_t)8);
  HS_EXPECT_EQ(v[7], 80);
}

inline void test_arenavec_clear() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  ArenaVector<int> v(a, 4);
  v.push_back(1);
  v.push_back(2);
  size_t cap = v.capacity();
  v.clear();
  HS_EXPECT_EQ(v.size(), (size_t)0);
  HS_EXPECT_EQ(v.capacity(), cap);
  HS_EXPECT_TRUE(v.is_bound());
  HS_EXPECT_TRUE(v.empty());
  HS_EXPECT_TRUE(v.is_empty());

  // Can refill after clear
  v.push_back(99);
  HS_EXPECT_EQ(v[0], 99);
}

inline void test_arenavec_iteration() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  ArenaVector<int> v(a, 8);
  for (int i = 0; i < 5; ++i) v.push_back(i);

  int sum = 0;
  for (auto it = v.begin(); it != v.end(); ++it) sum += *it;
  HS_EXPECT_EQ(sum, 0 + 1 + 2 + 3 + 4);

  // Range-based for compiles and visits each element
  int count = 0;
  for (int x : v) { (void)x; ++count; }
  HS_EXPECT_EQ(count, 5);
}

inline void test_arenavec_move_construct() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  ArenaVector<int> src(a, 4);
  src.push_back(7);
  src.push_back(8);
  int *src_data = src.data();

  ArenaVector<int> dst(std::move(src));
  HS_EXPECT_FALSE(src.is_bound());
  HS_EXPECT_EQ(src.size(), (size_t)0);
  HS_EXPECT_EQ(src.capacity(), (size_t)0);
  HS_EXPECT_TRUE(dst.is_bound());
  HS_EXPECT_EQ(dst.size(), (size_t)2);
  HS_EXPECT_EQ(dst.data(), src_data); // pointer moved, not copied
  HS_EXPECT_EQ(dst[0], 7);
  HS_EXPECT_EQ(dst[1], 8);
}

inline void test_arenavec_move_assign() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  ArenaVector<int> src(a, 4);
  src.push_back(1);
  src.push_back(2);
  src.push_back(3);
  int *src_data = src.data();

  ArenaVector<int> dst;
  dst = std::move(src);
  HS_EXPECT_FALSE(src.is_bound());
  HS_EXPECT_TRUE(dst.is_bound());
  HS_EXPECT_EQ(dst.size(), (size_t)3);
  HS_EXPECT_EQ(dst.data(), src_data);
}

inline void test_arenavec_rebind_reuses() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  ArenaVector<int> v(a, 16);
  v.push_back(42);
  int *data_before = v.data();

  // Re-binding with smaller-or-equal capacity should reuse storage and
  // reset size — not re-allocate.
  v.bind(a, 8);
  HS_EXPECT_EQ(v.size(), (size_t)0);
  HS_EXPECT_EQ(v.capacity(), (size_t)16);
  HS_EXPECT_EQ(v.data(), data_before);
}

// Re-binding to a LARGER capacity is a supported grow: it allocates a fresh
// block (abandoning the old one until the next arena reset — see bind()'s note)
// and adopts the new capacity. A debug build also logs the abandoned bytes.
inline void test_arenavec_rebind_grows() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  ArenaVector<int> v(a, 8);
  v.push_back(7);
  int *data_before = v.data();

  v.bind(a, 32);
  HS_EXPECT_EQ(v.size(), (size_t)0);
  HS_EXPECT_EQ(v.capacity(), (size_t)32);
  HS_EXPECT_NE(v.data(), data_before); // fresh allocation, old block abandoned
}

// Use-after-free precondition: an arena reset bumps the generation the vector
// snapshotted at bind(), marking its block dead. Rebinding a still-bound vector
// to a reset (or different) arena is therefore a CONTRACT VIOLATION — the old
// block is dead, and release builds cannot detect the staleness (generation
// tracking is debug-only), so bind() asserts on it in debug and release trusts
// the contract. Callers must clear/reconstruct the handle before resetting the
// arena (the Persist/compaction tests below do exactly this — they reset to an
// unbound state before re-binding). This test pins the generation-bump
// precondition that both the bind() contract assert and check_alive() rely on;
// the contract trap itself is a debug-only assert, exercised by integration
// rather than spawned as a death case (it lowers to abort, not __builtin_trap).
#ifndef NDEBUG
inline void test_arenavec_stale_binding_after_reset() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  ArenaVector<int> v(a, 16);
  v.push_back(123);
  uint32_t gen_before = a.get_generation();

  a.reset(); // generation bumps -> v's binding is now stale
  HS_EXPECT_TRUE(a.get_generation() != gen_before);
  // Calling v.bind(a, ...) here would now (correctly) trip the stale-binding
  // contract assert in bind(), so it is not exercised in-process.
}
#endif

inline void test_arenavec_zero_capacity() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  ArenaVector<int> v(a, 0);
  HS_EXPECT_TRUE(v.is_bound());
  HS_EXPECT_EQ(v.capacity(), (size_t)0);
  HS_EXPECT_TRUE(v.empty());
}

// ============================================================================
// ArenaSpan
// ============================================================================

inline void test_arenaspan_default() {
  ArenaSpan<int> sp;
  HS_EXPECT_EQ(sp.size(), (size_t)0);
  HS_EXPECT_TRUE(sp.empty());
  HS_EXPECT_TRUE(sp.data() == nullptr);
  HS_EXPECT_EQ(sp.begin(), sp.end());
}

inline void test_arenaspan_from_vector() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  ArenaVector<int> v(a, 4);
  v.push_back(11);
  v.push_back(22);
  v.push_back(33);

  ArenaSpan<int> sp(v);
  HS_EXPECT_EQ(sp.size(), (size_t)3);
  HS_EXPECT_FALSE(sp.empty());
  HS_EXPECT_EQ(sp.data(), v.data());
  HS_EXPECT_EQ(sp[0], 11);
  HS_EXPECT_EQ(sp[1], 22);
  HS_EXPECT_EQ(sp[2], 33);

  // Iteration
  int sum = 0;
  for (int x : sp) sum += x;
  HS_EXPECT_EQ(sum, 66);
}

// ============================================================================
// ScratchScope
// ============================================================================

inline void test_scratch_basic_restore() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  a.allocate(256);
  size_t before = a.get_offset();
  {
    ScratchScope s(a);
    s.get_arena().allocate(1024);
    HS_EXPECT_TRUE(a.get_offset() > before);
  }
  HS_EXPECT_EQ(a.get_offset(), before);
}

inline void test_scratch_make_vector() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  size_t before = a.get_offset();
  {
    ScratchScope s(a);
    ArenaVector<int> v = s.make_vector<int>(16);
    HS_EXPECT_TRUE(v.is_bound());
    HS_EXPECT_EQ(v.capacity(), (size_t)16);
    for (int i = 0; i < 16; ++i) v.push_back(i);
    HS_EXPECT_EQ(v[10], 10);
  }
  HS_EXPECT_EQ(a.get_offset(), before);
}

inline void test_scratch_nested() {
  Arena a(test_buf_a, sizeof(test_buf_a));
  size_t l0 = a.get_offset();
  {
    ScratchScope s1(a);
    a.allocate(512);
    size_t l1 = a.get_offset();
    {
      ScratchScope s2(a);
      a.allocate(1024);
      HS_EXPECT_TRUE(a.get_offset() > l1);
    }
    HS_EXPECT_EQ(a.get_offset(), l1);
  }
  HS_EXPECT_EQ(a.get_offset(), l0);
}

// ============================================================================
// Persist<T> — Cloneable evacuator
// ============================================================================

struct TestPayload {
  ArenaVector<int> data;
  int summary = 0;

  TestPayload() = default;
  TestPayload(TestPayload &&) noexcept = default;
  TestPayload &operator=(TestPayload &&) noexcept = default;

  static void clone(const TestPayload &src, TestPayload &dst, Arena &arena) {
    dst.data.bind(arena, src.data.size());
    for (size_t i = 0; i < src.data.size(); ++i) dst.data.push_back(src.data[i]);
    dst.summary = src.summary;
  }
};

inline void test_persist_restores_target() {
  Arena persistent(test_buf_a, sizeof(test_buf_a));
  Arena scratch(test_buf_b, sizeof(test_buf_b));

  TestPayload live;
  live.data.bind(persistent, 4);
  live.data.push_back(7);
  live.data.push_back(11);
  live.data.push_back(13);
  live.summary = 99;

  {
    Persist<TestPayload> p(live, scratch, persistent);
    // Inside the scope, simulate a destructive operation on persistent
    persistent.reset();
    live = TestPayload();
    HS_EXPECT_EQ(live.data.size(), (size_t)0);
    HS_EXPECT_EQ(live.summary, 0);
  }

  // After scope exit, live is rebuilt from the scratch backup
  HS_EXPECT_EQ(live.data.size(), (size_t)3);
  HS_EXPECT_EQ(live.data[0], 7);
  HS_EXPECT_EQ(live.data[1], 11);
  HS_EXPECT_EQ(live.data[2], 13);
  HS_EXPECT_EQ(live.summary, 99);
}

inline void test_persist_scratch_offset_restored() {
  Arena persistent(test_buf_a, sizeof(test_buf_a));
  Arena scratch(test_buf_b, sizeof(test_buf_b));

  TestPayload live;
  live.data.bind(persistent, 2);
  live.data.push_back(1);
  live.summary = 5;

  size_t scratch_before = scratch.get_offset();
  {
    Persist<TestPayload> p(live, scratch, persistent);
    HS_EXPECT_TRUE(scratch.get_offset() > scratch_before);
    persistent.reset();
    live = TestPayload();
  }
  // The Persist<>'s embedded ScratchScope must roll scratch back.
  HS_EXPECT_EQ(scratch.get_offset(), scratch_before);
}

// Compaction: the documented Persist use case (HankinSolids' shape rebuild) is
// not a bare reset — it evacuates survivors, resets persistent, RE-LAYS the
// long-lived data more tightly, then restores. The survivor must come back
// intact at its NEW (relocated) address, sitting after the compacted data
// without clobbering it.
inline void test_persist_compaction_relocates_survivor() {
  Arena persistent(test_buf_a, sizeof(test_buf_a));
  Arena scratch(test_buf_b, sizeof(test_buf_b));

  // Pre-compaction: a large junk block FIRST pushes the survivor to a high base.
  ArenaVector<int> junk;
  junk.bind(persistent, 256);
  for (int i = 0; i < 256; ++i) junk.push_back(i);

  TestPayload live;
  live.data.bind(persistent, 3);
  live.data.push_back(100);
  live.data.push_back(200);
  live.data.push_back(300);
  live.summary = 42;
  const int *base_before = &live.data[0];

  ArenaVector<int> compacted_other; // outlives the Persist block
  {
    Persist<TestPayload> p(live, scratch, persistent);
    persistent.reset();   // free junk + survivor
    live = TestPayload();
    // Compaction: place the survivors of the rebuild first, far smaller than
    // the old junk, so the restored survivor lands at a different address.
    compacted_other.bind(persistent, 2);
    compacted_other.push_back(7);
    compacted_other.push_back(8);
  } // ~Persist clones the backup into persistent AFTER compacted_other

  // Survivor restored intact...
  HS_EXPECT_EQ(live.data.size(), (size_t)3);
  HS_EXPECT_EQ(live.data[0], 100);
  HS_EXPECT_EQ(live.data[1], 200);
  HS_EXPECT_EQ(live.data[2], 300);
  HS_EXPECT_EQ(live.summary, 42);
  // ...relocated (no longer sitting after the now-freed junk)...
  HS_EXPECT_TRUE(&live.data[0] != base_before);
  // ...and the compacted data it was restored after is untouched.
  HS_EXPECT_EQ(compacted_other.size(), (size_t)2);
  HS_EXPECT_EQ(compacted_other[0], 7);
  HS_EXPECT_EQ(compacted_other[1], 8);
}

// ============================================================================
// Runner
// ============================================================================

inline int run_memory_tests() {
  auto scope = hs_test::begin_module("memory");

  test_arena_construction();
  test_arena_basic_allocation();
  test_arena_alignment();
  test_arena_high_water_mark();
  test_arena_reset();
  test_arena_set_offset();
  test_arena_fills_to_capacity();
  test_arena_rebind();
  test_arena_reset_high_water_mark();
#ifndef NDEBUG
  test_arena_generation_bumps();
#endif

  test_tribitset_sizes();
  test_tribitset_clear_and_test();
  test_tribitset_test_and_set();
  test_tribitset_index_uniqueness();
  test_tribitset_all_pairs_independent();

  test_arenavec_default_unbound();
  test_arenavec_bind();
  test_arenavec_constructor_with_arena();
  test_arenavec_push_back_indexing();
  test_arenavec_back();
  test_arenavec_emplace_back();
  test_arenavec_append_bulk();
  test_arenavec_clear();
  test_arenavec_iteration();
  test_arenavec_move_construct();
  test_arenavec_move_assign();
  test_arenavec_rebind_reuses();
  test_arenavec_rebind_grows();
#ifndef NDEBUG
  test_arenavec_stale_binding_after_reset();
#endif
  test_arenavec_zero_capacity();

  test_arenaspan_default();
  test_arenaspan_from_vector();

  test_scratch_basic_restore();
  test_scratch_make_vector();
  test_scratch_nested();

  test_persist_restores_target();
  test_persist_scratch_offset_restored();
  test_persist_compaction_relocates_survivor();

  return hs_test::end_module(scope);
}

} // namespace mem
} // namespace hs_test

