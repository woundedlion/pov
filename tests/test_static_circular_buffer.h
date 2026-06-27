/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/static_circular_buffer.h.
 *
 * NOTE: front(), back(), and operator[] abort on misuse, so the empty/OOB
 * paths are not exercised here — they are driven (in forked child processes)
 * by the death harness in tests/test_death.h: case_circular_buffer_front_empty
 * (front() while empty) and case_circular_buffer_oob (operator[] past the live
 * count).
 */
#pragma once

#include <iterator>
#include "core/static_circular_buffer.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace scb {

// ============================================================================
// Construction
// ============================================================================

/**
 * @brief Verifies a freshly constructed buffer is empty, reports zero size,
 *        full capacity, and has begin() == end().
 */
inline void test_default_state() {
  StaticCircularBuffer<int, 8> buf;
  HS_EXPECT_TRUE(buf.is_empty());
  HS_EXPECT_FALSE(buf.is_full());
  HS_EXPECT_EQ(buf.size(), (size_t)0);
  HS_EXPECT_EQ(buf.capacity(), (size_t)8);
  HS_EXPECT_EQ(buf.begin(), buf.end());
}

/**
 * @brief Verifies an initializer list smaller than capacity populates the
 *        buffer in order.
 */
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

/**
 * @brief Verifies an initializer list filling the buffer exactly to capacity.
 * @details An over-capacity list is a compile-time error via static_assert in
 *          the variadic constructor, so only the at-capacity fill is exercised.
 */
inline void test_initializer_list_at_capacity() {
  StaticCircularBuffer<int, 3> buf{1, 2, 3};
  HS_EXPECT_EQ(buf.size(), (size_t)3);
  HS_EXPECT_TRUE(buf.is_full());
  HS_EXPECT_EQ(buf[0], 1);
  HS_EXPECT_EQ(buf[1], 2);
  HS_EXPECT_EQ(buf[2], 3);
  HS_EXPECT_EQ(buf.front(), 1);
  HS_EXPECT_EQ(buf.back(), 3);
}

// ============================================================================
// push_back / push_front basics
// ============================================================================

/**
 * @brief Verifies push_back appends to the tail and elements stay in insertion
 *        order.
 */
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

/**
 * @brief Verifies pushing exactly capacity items marks the buffer full.
 */
inline void test_push_back_fills_to_capacity() {
  StaticCircularBuffer<int, 3> buf;
  buf.push_back(1);
  buf.push_back(2);
  buf.push_back(3);
  HS_EXPECT_TRUE(buf.is_full());
  HS_EXPECT_EQ(buf.size(), (size_t)3);
}

/**
 * @brief Verifies push_back past capacity drops the oldest (front) element.
 */
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

/**
 * @brief Verifies push_front prepends to the head so the most recently pushed
 *        item is at front.
 */
inline void test_push_front_order() {
  StaticCircularBuffer<int, 4> buf;
  buf.push_front(1);
  buf.push_front(2);
  buf.push_front(3);
  HS_EXPECT_EQ(buf.size(), (size_t)3);
  HS_EXPECT_EQ(buf.front(), 3);
  HS_EXPECT_EQ(buf.back(), 1);
  HS_EXPECT_EQ(buf[0], 3);
  HS_EXPECT_EQ(buf[1], 2);
  HS_EXPECT_EQ(buf[2], 1);
}

/**
 * @brief Verifies push_front past capacity evicts the tail so the newest 3
 *        items survive with the latest push at front.
 */
inline void test_push_front_overflow_drops_back() {
  StaticCircularBuffer<int, 3> buf;
  for (int i = 1; i <= 5; ++i) buf.push_front(i);
  HS_EXPECT_EQ(buf.size(), (size_t)3);
  HS_EXPECT_EQ(buf.front(), 5);
  HS_EXPECT_EQ(buf.back(), 3);
  HS_EXPECT_EQ(buf[0], 5);
  HS_EXPECT_EQ(buf[1], 4);
  HS_EXPECT_EQ(buf[2], 3);
}

/**
 * @brief Verifies push_back accepts an rvalue (move) overload.
 */
inline void test_push_back_rvalue() {
  StaticCircularBuffer<int, 4> buf;
  int x = 7;
  buf.push_back(std::move(x));
  HS_EXPECT_EQ(buf.size(), (size_t)1);
  HS_EXPECT_EQ(buf[0], 7);
}

/**
 * @brief Verifies push_front accepts an rvalue (move) overload.
 */
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

/**
 * @brief Default- and value-constructible 2D point used to test in-place
 *        construction.
 */
struct Pt {
  int x, y; /**< Cartesian coordinates of the point. */
  /**
   * @brief Default-constructs the point at the origin (0, 0).
   */
  Pt() : x(0), y(0) {}
  /**
   * @brief Constructs the point at (a, b).
   * @param a X coordinate.
   * @param b Y coordinate.
   */
  Pt(int a, int b) : x(a), y(b) {}
};

/**
 * @brief Verifies emplace_back constructs the element in place and returns a
 *        reference to it.
 */
inline void test_emplace_back_constructs_in_place() {
  StaticCircularBuffer<Pt, 4> buf;
  Pt &ref = buf.emplace_back(3, 4);
  HS_EXPECT_EQ(ref.x, 3);
  HS_EXPECT_EQ(ref.y, 4);
  HS_EXPECT_EQ(buf.size(), (size_t)1);
  HS_EXPECT_EQ(buf[0].x, 3);
  HS_EXPECT_EQ(buf[0].y, 4);
}

/**
 * @brief Verifies emplace_front constructs the element in place at the head and
 *        returns it.
 */
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

/**
 * @brief Constructible but non-assignable type.
 * @details emplace must construct directly in the slot's storage rather than
 *          assign a temporary; using this type forces that path to compile and
 *          run.
 */
struct EmplaceOnly {
  int x, y; /**< Stored coordinate pair. */
  /**
   * @brief Default-constructs with sentinel coordinates (-1, -1).
   */
  EmplaceOnly() : x(-1), y(-1) {}
  /**
   * @brief Constructs at (a, b).
   * @param a X coordinate.
   * @param b Y coordinate.
   */
  EmplaceOnly(int a, int b) : x(a), y(b) {}
  /**
   * @brief Defaulted copy constructor.
   * @param other Source instance to copy from.
   */
  EmplaceOnly(const EmplaceOnly &) = default;
  /**
   * @brief Deleted copy assignment; the type is intentionally non-assignable.
   * @param other Source instance (unused; overload is deleted).
   * @return Reference to this (never invoked).
   */
  EmplaceOnly &operator=(const EmplaceOnly &) = delete;
  /**
   * @brief Deleted move assignment; the type is intentionally non-assignable.
   * @param other Source instance (unused; overload is deleted).
   * @return Reference to this (never invoked).
   */
  EmplaceOnly &operator=(EmplaceOnly &&) = delete;
};

/**
 * @brief Verifies emplace_back/emplace_front work for a type with no assignment
 *        operator, proving they construct in place.
 */
inline void test_emplace_constructs_non_assignable_type() {
  StaticCircularBuffer<EmplaceOnly, 3> buf;
  EmplaceOnly &b = buf.emplace_back(3, 4);
  HS_EXPECT_EQ(b.x, 3);
  HS_EXPECT_EQ(b.y, 4);
  EmplaceOnly &f = buf.emplace_front(7, 9);
  HS_EXPECT_EQ(f.x, 7);
  HS_EXPECT_EQ(buf.size(), (size_t)2);
  HS_EXPECT_EQ(buf.front().x, 7);
  HS_EXPECT_EQ(buf.back().x, 3);
}

// ============================================================================
// pop / pop_back / clear
// ============================================================================

/**
 * @brief Verifies pop removes the front (oldest) element.
 */
inline void test_pop_removes_front() {
  StaticCircularBuffer<int, 4> buf{1, 2, 3};
  buf.pop();
  HS_EXPECT_EQ(buf.size(), (size_t)2);
  HS_EXPECT_EQ(buf.front(), 2);
  HS_EXPECT_EQ(buf.back(), 3);
}

/**
 * @brief Verifies pop_back removes the back (newest) element.
 */
inline void test_pop_back_removes_last() {
  StaticCircularBuffer<int, 4> buf{1, 2, 3};
  buf.pop_back();
  HS_EXPECT_EQ(buf.size(), (size_t)2);
  HS_EXPECT_EQ(buf.front(), 1);
  HS_EXPECT_EQ(buf.back(), 2);
}

/**
 * @brief Verifies pop and pop_back on an empty buffer are safe no-ops (they
 *        must not crash).
 */
inline void test_pop_on_empty_is_noop() {
  StaticCircularBuffer<int, 4> buf;
  buf.pop();
  buf.pop_back();
  HS_EXPECT_TRUE(buf.is_empty());
  HS_EXPECT_EQ(buf.size(), (size_t)0);
}

/**
 * @brief Verifies clear empties the buffer (size 0) without changing capacity,
 *        and the buffer stays usable afterward.
 */
inline void test_clear_empties_buffer() {
  StaticCircularBuffer<int, 4> buf{1, 2, 3};
  buf.clear();
  HS_EXPECT_TRUE(buf.is_empty());
  HS_EXPECT_EQ(buf.size(), (size_t)0);
  HS_EXPECT_EQ(buf.capacity(), (size_t)4);
  HS_EXPECT_FALSE(buf.is_full());

  buf.push_back(42);
  HS_EXPECT_EQ(buf.size(), (size_t)1);
  HS_EXPECT_EQ(buf[0], 42);
}

/**
 * @brief Verifies interleaving pushes and pops drives the head/tail indices
 *        around the ring while the buffer stays consistent and back() always
 *        holds the most recent push.
 */
inline void test_alternating_push_pop_wraps() {
  StaticCircularBuffer<int, 3> buf;
  for (int i = 0; i < 10; ++i) {
    buf.push_back(i);
    if (i % 2 == 1) buf.pop();
  }
  // Final window after the wrap: push 9 fills [7,8,9], then the i=9 pop drops 7.
  HS_EXPECT_EQ(buf.size(), (size_t)2);
  HS_EXPECT_EQ(buf.front(), 8);
  HS_EXPECT_EQ(buf.back(), 9);
  HS_EXPECT_EQ(buf[0], 8);
  HS_EXPECT_EQ(buf[1], 9);
}

// ============================================================================
// Iterators
// ============================================================================

/**
 * @brief Verifies begin()..end() iteration visits every element in logical
 *        (front-to-back) order.
 */
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

/**
 * @brief Verifies the nested iterator and const_iterator types are publicly
 *        nameable from outside the class (not just usable via `auto`), and that
 *        const_iterator is constructible from a non-const iterator.
 */
inline void test_iterator_types_are_public() {
  using Buf = StaticCircularBuffer<int, 4>;
  Buf buf{10, 20, 30};
  int expected = 10;
  for (Buf::iterator it = buf.begin(); it != buf.end(); ++it) {
    HS_EXPECT_EQ(*it, expected);
    expected += 10;
  }
  const Buf &cbuf = buf;
  Buf::const_iterator cit = cbuf.begin();
  HS_EXPECT_EQ(*cit, 10);
  Buf::const_iterator from_mut = buf.begin();
  HS_EXPECT_EQ(*from_mut, 10);
}

/**
 * @brief Verifies the buffer works with a range-based for loop.
 */
inline void test_iterator_range_for() {
  StaticCircularBuffer<int, 4> buf{1, 2, 3, 4};
  int sum = 0;
  for (int x : buf) sum += x;
  HS_EXPECT_EQ(sum, 10);
}

/**
 * @brief Verifies std::distance(begin, end) equals size().
 */
inline void test_iterator_distance_matches_size() {
  StaticCircularBuffer<int, 4> buf{1, 2, 3};
  HS_EXPECT_EQ(std::distance(buf.begin(), buf.end()), (std::ptrdiff_t)3);
}

/**
 * @brief Verifies the iterator supports random-access operations: +, -, +=, -=,
 *        [], and end-begin.
 */
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

/**
 * @brief Verifies the iterator supports the full set of relational and equality
 *        comparisons.
 */
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

/**
 * @brief Verifies that when the head index has wrapped into the middle of the
 *        underlying array, iteration still visits all elements in logical order.
 */
inline void test_iterator_after_wrap() {
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

/**
 * @brief Verifies a const buffer yields const_iterators that traverse correctly,
 *        and a non-const iterator converts implicitly to a const_iterator.
 */
inline void test_const_iterator() {
  StaticCircularBuffer<int, 4> buf{1, 2, 3};
  const auto &cref = buf;
  int sum = 0;
  for (auto it = cref.begin(); it != cref.end(); ++it) sum += *it;
  HS_EXPECT_EQ(sum, 6);

  auto ci = cref.begin();
  ci = buf.begin();
  HS_EXPECT_EQ(*ci, 1);
}

/**
 * @brief Verifies the iterator's operator-> accesses members of the pointed-to
 *        element.
 */
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

/**
 * @brief Verifies operator[] returns a mutable reference, allowing in-place
 *        element assignment.
 */
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

/**
 * @brief Runs every StaticCircularBuffer test case.
 * @return The module's failure count (number of failed assertions).
 */
inline int run_static_circular_buffer_tests() {
  hs_test::ModuleFixture fixture("static_circular_buffer");

  test_default_state();
  test_initializer_list_within_capacity();
  test_initializer_list_at_capacity();

  test_push_back_order();
  test_push_back_fills_to_capacity();
  test_push_back_overflow_drops_oldest();
  test_push_front_order();
  test_push_front_overflow_drops_back();
  test_push_back_rvalue();
  test_push_front_rvalue();

  test_emplace_back_constructs_in_place();
  test_emplace_front_constructs_in_place();
  test_emplace_constructs_non_assignable_type();

  test_pop_removes_front();
  test_pop_back_removes_last();
  test_pop_on_empty_is_noop();
  test_clear_empties_buffer();
  test_alternating_push_pop_wraps();

  test_iterator_traversal();
  test_iterator_types_are_public();
  test_iterator_range_for();
  test_iterator_distance_matches_size();
  test_iterator_random_access();
  test_iterator_comparisons();
  test_iterator_after_wrap();
  test_const_iterator();
  test_iterator_arrow_operator();

  test_index_assignment();

  return fixture.result();
}

} // namespace scb
} // namespace hs_test

