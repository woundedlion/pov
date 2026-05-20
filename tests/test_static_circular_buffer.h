/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/static_circular_buffer.h.
 *
 * NOTE: front(), back(), and operator[] abort on misuse, so the empty/OOB
 * paths are NOT exercised here — those preconditions are documented as
 * caller responsibilities.
 */
#pragma once
#ifndef HOLOSPHERE_TESTS_TEST_STATIC_CIRCULAR_BUFFER_H_
#define HOLOSPHERE_TESTS_TEST_STATIC_CIRCULAR_BUFFER_H_

#include <iterator>
#include "core/static_circular_buffer.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace scb {

// ============================================================================
// Construction
// ============================================================================

inline void test_default_state() {
  StaticCircularBuffer<int, 8> buf;
  HS_EXPECT_TRUE(buf.is_empty());
  HS_EXPECT_FALSE(buf.is_full());
  HS_EXPECT_EQ(buf.size(), (size_t)0);
  HS_EXPECT_EQ(buf.capacity(), (size_t)8);
  HS_EXPECT_EQ(buf.begin(), buf.end());
}

inline void test_initializer_list_within_capacity() {
  StaticCircularBuffer<int, 4> buf{10, 20, 30};
  HS_EXPECT_EQ(buf.size(), (size_t)3);
  HS_EXPECT_FALSE(buf.is_full());
  HS_EXPECT_EQ(buf[0], 10);
  HS_EXPECT_EQ(buf[1], 20);
  HS_EXPECT_EQ(buf[2], 30);
  HS_EXPECT_EQ(buf.front(), 10);
  HS_EXPECT_EQ(buf.back(), 30);
}

inline void test_initializer_list_overflow() {
  // Items beyond capacity evict the oldest entries (FIFO) — only the last
  // N items survive.
  StaticCircularBuffer<int, 3> buf{1, 2, 3, 4, 5};
  HS_EXPECT_EQ(buf.size(), (size_t)3);
  HS_EXPECT_TRUE(buf.is_full());
  HS_EXPECT_EQ(buf[0], 3);
  HS_EXPECT_EQ(buf[1], 4);
  HS_EXPECT_EQ(buf[2], 5);
  HS_EXPECT_EQ(buf.front(), 3);
  HS_EXPECT_EQ(buf.back(), 5);
}

// ============================================================================
// push_back / push_front basics
// ============================================================================

inline void test_push_back_order() {
  StaticCircularBuffer<int, 4> buf;
  buf.push_back(1);
  buf.push_back(2);
  buf.push_back(3);
  HS_EXPECT_EQ(buf.size(), (size_t)3);
  HS_EXPECT_EQ(buf[0], 1);
  HS_EXPECT_EQ(buf[1], 2);
  HS_EXPECT_EQ(buf[2], 3);
  HS_EXPECT_EQ(buf.front(), 1);
  HS_EXPECT_EQ(buf.back(), 3);
  HS_EXPECT_FALSE(buf.is_empty());
  HS_EXPECT_FALSE(buf.is_full());
}

inline void test_push_back_fills_to_capacity() {
  StaticCircularBuffer<int, 3> buf;
  buf.push_back(1);
  buf.push_back(2);
  buf.push_back(3);
  HS_EXPECT_TRUE(buf.is_full());
  HS_EXPECT_EQ(buf.size(), (size_t)3);
}

inline void test_push_back_overflow_drops_oldest() {
  StaticCircularBuffer<int, 3> buf;
  for (int i = 1; i <= 5; ++i) buf.push_back(i);
  HS_EXPECT_EQ(buf.size(), (size_t)3);
  HS_EXPECT_EQ(buf[0], 3);
  HS_EXPECT_EQ(buf[1], 4);
  HS_EXPECT_EQ(buf[2], 5);
  HS_EXPECT_EQ(buf.front(), 3);
  HS_EXPECT_EQ(buf.back(), 5);
}

inline void test_push_front_order() {
  StaticCircularBuffer<int, 4> buf;
  buf.push_front(1);
  buf.push_front(2);
  buf.push_front(3);
  // Newest-first: 3 is at the front
  HS_EXPECT_EQ(buf.size(), (size_t)3);
  HS_EXPECT_EQ(buf.front(), 3);
  HS_EXPECT_EQ(buf.back(), 1);
  HS_EXPECT_EQ(buf[0], 3);
  HS_EXPECT_EQ(buf[1], 2);
  HS_EXPECT_EQ(buf[2], 1);
}

inline void test_push_front_overflow_drops_back() {
  StaticCircularBuffer<int, 3> buf;
  for (int i = 1; i <= 5; ++i) buf.push_front(i);
  // After 5 push_fronts on capacity-3 buffer: newest 3 items are 5,4,3
  // with 5 at front and 3 at back (since push_front always evicts the tail).
  HS_EXPECT_EQ(buf.size(), (size_t)3);
  HS_EXPECT_EQ(buf.front(), 5);
  HS_EXPECT_EQ(buf.back(), 3);
  HS_EXPECT_EQ(buf[0], 5);
  HS_EXPECT_EQ(buf[1], 4);
  HS_EXPECT_EQ(buf[2], 3);
}

inline void test_push_back_rvalue() {
  StaticCircularBuffer<int, 4> buf;
  int x = 7;
  buf.push_back(std::move(x));
  HS_EXPECT_EQ(buf.size(), (size_t)1);
  HS_EXPECT_EQ(buf[0], 7);
}

inline void test_push_front_rvalue() {
  StaticCircularBuffer<int, 4> buf;
  int x = 99;
  buf.push_front(std::move(x));
  HS_EXPECT_EQ(buf.size(), (size_t)1);
  HS_EXPECT_EQ(buf[0], 99);
}

// ============================================================================
// emplace_back / emplace_front
// ============================================================================

struct Pt {
  int x, y;
  Pt() : x(0), y(0) {}
  Pt(int a, int b) : x(a), y(b) {}
};

inline void test_emplace_back_constructs_in_place() {
  StaticCircularBuffer<Pt, 4> buf;
  Pt &ref = buf.emplace_back(3, 4);
  HS_EXPECT_EQ(ref.x, 3);
  HS_EXPECT_EQ(ref.y, 4);
  HS_EXPECT_EQ(buf.size(), (size_t)1);
  HS_EXPECT_EQ(buf[0].x, 3);
  HS_EXPECT_EQ(buf[0].y, 4);
}

inline void test_emplace_front_constructs_in_place() {
  StaticCircularBuffer<Pt, 4> buf;
  buf.emplace_back(1, 1);
  Pt &ref = buf.emplace_front(9, 8);
  HS_EXPECT_EQ(ref.x, 9);
  HS_EXPECT_EQ(ref.y, 8);
  HS_EXPECT_EQ(buf.size(), (size_t)2);
  HS_EXPECT_EQ(buf.front().x, 9);
  HS_EXPECT_EQ(buf.back().x, 1);
}

// ============================================================================
// pop / pop_back / clear
// ============================================================================

inline void test_pop_removes_front() {
  StaticCircularBuffer<int, 4> buf{1, 2, 3};
  buf.pop(); // remove front
  HS_EXPECT_EQ(buf.size(), (size_t)2);
  HS_EXPECT_EQ(buf.front(), 2);
  HS_EXPECT_EQ(buf.back(), 3);
}

inline void test_pop_back_removes_last() {
  StaticCircularBuffer<int, 4> buf{1, 2, 3};
  buf.pop_back();
  HS_EXPECT_EQ(buf.size(), (size_t)2);
  HS_EXPECT_EQ(buf.front(), 1);
  HS_EXPECT_EQ(buf.back(), 2);
}

inline void test_pop_on_empty_is_noop() {
  StaticCircularBuffer<int, 4> buf;
  buf.pop();      // must not crash
  buf.pop_back(); // must not crash
  HS_EXPECT_TRUE(buf.is_empty());
  HS_EXPECT_EQ(buf.size(), (size_t)0);
}

inline void test_clear_empties_buffer() {
  StaticCircularBuffer<int, 4> buf{1, 2, 3};
  buf.clear();
  HS_EXPECT_TRUE(buf.is_empty());
  HS_EXPECT_EQ(buf.size(), (size_t)0);
  HS_EXPECT_EQ(buf.capacity(), (size_t)4); // capacity unchanged
  HS_EXPECT_FALSE(buf.is_full());

  // Buffer remains usable after clear
  buf.push_back(42);
  HS_EXPECT_EQ(buf.size(), (size_t)1);
  HS_EXPECT_EQ(buf[0], 42);
}

inline void test_alternating_push_pop_wraps() {
  // Exercise the ring wraparound by interleaving pushes and pops.
  StaticCircularBuffer<int, 3> buf;
  for (int i = 0; i < 10; ++i) {
    buf.push_back(i);
    if (i % 2 == 1) buf.pop();
  }
  // After this sequence, the buffer must remain consistent.
  HS_EXPECT_TRUE(buf.size() <= 3);
  // Verify back is always the most recent value pushed
  HS_EXPECT_EQ(buf.back(), 9);
}

// ============================================================================
// Iterators
// ============================================================================

inline void test_iterator_traversal() {
  StaticCircularBuffer<int, 4> buf{10, 20, 30};
  int expected = 10;
  int count = 0;
  for (auto it = buf.begin(); it != buf.end(); ++it) {
    HS_EXPECT_EQ(*it, expected);
    expected += 10;
    ++count;
  }
  HS_EXPECT_EQ(count, 3);
}

inline void test_iterator_range_for() {
  StaticCircularBuffer<int, 4> buf{1, 2, 3, 4};
  int sum = 0;
  for (int x : buf) sum += x;
  HS_EXPECT_EQ(sum, 10);
}

inline void test_iterator_distance_matches_size() {
  StaticCircularBuffer<int, 4> buf{1, 2, 3};
  HS_EXPECT_EQ(std::distance(buf.begin(), buf.end()), (std::ptrdiff_t)3);
}

inline void test_iterator_random_access() {
  StaticCircularBuffer<int, 4> buf{10, 20, 30, 40};
  auto it = buf.begin();
  HS_EXPECT_EQ(*(it + 2), 30);
  HS_EXPECT_EQ(it[3], 40);
  it += 2;
  HS_EXPECT_EQ(*it, 30);
  it -= 1;
  HS_EXPECT_EQ(*it, 20);
  auto end = buf.end();
  HS_EXPECT_EQ(end - buf.begin(), (std::ptrdiff_t)4);
}

inline void test_iterator_comparisons() {
  StaticCircularBuffer<int, 4> buf{1, 2, 3};
  auto a = buf.begin();
  auto b = buf.begin() + 1;
  HS_EXPECT_TRUE(a < b);
  HS_EXPECT_TRUE(b > a);
  HS_EXPECT_TRUE(a <= a);
  HS_EXPECT_TRUE(b >= a);
  HS_EXPECT_TRUE(a != b);
  HS_EXPECT_FALSE(a == b);
}

inline void test_iterator_after_wrap() {
  // Force the head pointer to wrap into the middle of the underlying array,
  // then verify iteration still visits all elements in logical order.
  StaticCircularBuffer<int, 3> buf;
  for (int i = 1; i <= 5; ++i) buf.push_back(i); // ends as [3, 4, 5] logical
  int values[3];
  int n = 0;
  for (int x : buf) values[n++] = x;
  HS_EXPECT_EQ(n, 3);
  HS_EXPECT_EQ(values[0], 3);
  HS_EXPECT_EQ(values[1], 4);
  HS_EXPECT_EQ(values[2], 5);
}

inline void test_const_iterator() {
  StaticCircularBuffer<int, 4> buf{1, 2, 3};
  const auto &cref = buf;
  int sum = 0;
  for (auto it = cref.begin(); it != cref.end(); ++it) sum += *it;
  HS_EXPECT_EQ(sum, 6);

  // Conversion from non-const iterator to const_iterator is allowed.
  auto ci = cref.begin();
  ci = buf.begin(); // implicit construction from non-const iterator
  HS_EXPECT_EQ(*ci, 1);
}

inline void test_iterator_arrow_operator() {
  StaticCircularBuffer<Pt, 4> buf;
  buf.emplace_back(7, 11);
  auto it = buf.begin();
  HS_EXPECT_EQ(it->x, 7);
  HS_EXPECT_EQ(it->y, 11);
}

// ============================================================================
// operator[] mutation
// ============================================================================

inline void test_index_assignment() {
  StaticCircularBuffer<int, 4> buf{10, 20, 30};
  buf[1] = 99;
  HS_EXPECT_EQ(buf[1], 99);
  HS_EXPECT_EQ(buf.front(), 10);
  HS_EXPECT_EQ(buf.back(), 30);
}

// ============================================================================
// Runner
// ============================================================================

inline int run_static_circular_buffer_tests() {
  auto scope = hs_test::begin_module("static_circular_buffer");

  test_default_state();
  test_initializer_list_within_capacity();
  test_initializer_list_overflow();

  test_push_back_order();
  test_push_back_fills_to_capacity();
  test_push_back_overflow_drops_oldest();
  test_push_front_order();
  test_push_front_overflow_drops_back();
  test_push_back_rvalue();
  test_push_front_rvalue();

  test_emplace_back_constructs_in_place();
  test_emplace_front_constructs_in_place();

  test_pop_removes_front();
  test_pop_back_removes_last();
  test_pop_on_empty_is_noop();
  test_clear_empties_buffer();
  test_alternating_push_pop_wraps();

  test_iterator_traversal();
  test_iterator_range_for();
  test_iterator_distance_matches_size();
  test_iterator_random_access();
  test_iterator_comparisons();
  test_iterator_after_wrap();
  test_const_iterator();
  test_iterator_arrow_operator();

  test_index_assignment();

  return hs_test::end_module(scope);
}

} // namespace scb
} // namespace hs_test

#endif // HOLOSPHERE_TESTS_TEST_STATIC_CIRCULAR_BUFFER_H_
