/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <new>
#include <utility>
#include "platform.h"

/**
 * @brief A fixed-size circular buffer optimized for stability.
 * @details
 * - No dynamic memory allocation (prevents heap fragmentation).
 * - "Delete Oldest" strategy on overflow (never stops processing) — overflow is
 *   a designed, non-fatal condition.
 * - Misuse — front()/back() on an empty buffer, or operator[] out of range — is
 *   an invariant violation with no valid recovery, so it HS_CHECK-traps
 *   (survives NDEBUG, pulls in no stdio) rather than reading garbage. Fail-fast,
 *   consistent with the rest of the engine.
 */
template <typename T, size_t N> class StaticCircularBuffer {
  // Defined privately at the bottom; surfaced through the public usings below so
  // external code can name StaticCircularBuffer<T,N>::iterator (the class is
  // otherwise fully STL-conformant).
  class Iterator;
  class ConstIterator;

public:
  using iterator = Iterator;
  using const_iterator = ConstIterator;
  using value_type = T;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using reference = T &;
  using const_reference = const T &;

  StaticCircularBuffer() : head(0), tail(0), count(0) {}

  StaticCircularBuffer(std::initializer_list<T> items)
      : head(0), tail(0), count(0) {
    if (items.size() > N) {
      // push_back evicts the oldest on overflow, so the buffer keeps the LAST
      // N items — say so rather than implying the leading items were kept.
      hs::log("Buffer Warning: Initializer list exceeds capacity. Keeping last "
              "N items.");
    }
    for (const T &item : items) {
      push_back(item);
    }
  }

  void push_front(const T &item) {
    if (is_full()) {
      // Drop the tail (oldest in this context) to make room
      pop_back_internal();
    }
    head = (head + N - 1) % N; // + N first: never underflows (head is uint32_t)
    buffer[head] = item;
    count++;
  }

  void push_front(T &&item) {
    if (is_full()) {
      pop_back_internal();
    }
    head = (head + N - 1) % N; // + N first: never underflows (head is uint32_t)
    buffer[head] = std::move(item);
    count++;
  }

  template <typename... Args> T &emplace_front(Args &&...args) {
    if (is_full()) {
      pop_back_internal();
    }
    // Commit head/count only after the (potentially throwing) constructor
    // succeeds, so a throwing T leaves the buffer's size invariant intact.
    uint32_t slot = (head + N - 1) % N; // + N first: never underflows
    T &ref = construct_in_place(slot, std::forward<Args>(args)...);
    head = slot;
    count++;
    return ref;
  }

  void push_back(const T &item) {
    if (is_full()) {
      // Drop the head (oldest in this context) to make room
      pop_front_internal();
    }
    buffer[tail] = item;
    tail = (tail + 1) % N;
    count++;
  }

  void push_back(T &&item) {
    if (is_full()) {
      pop_front_internal();
    }
    buffer[tail] = std::move(item);
    tail = (tail + 1) % N;
    count++;
  }

  template <typename... Args> T &emplace_back(Args &&...args) {
    if (is_full()) {
      pop_front_internal();
    }
    // Commit tail/count only after the (potentially throwing) constructor
    // succeeds, so a throwing T leaves the buffer's size invariant intact.
    uint32_t slot = tail;
    T &ref = construct_in_place(slot, std::forward<Args>(args)...);
    tail = (tail + 1) % N;
    count++;
    return ref;
  }

  void pop_back() {
    if (is_empty())
      return;
    pop_back_internal();
  }

  void pop() {
    if (is_empty())
      return;
    pop_front_internal();
  }

  void clear() {
    while (!is_empty()) {
      pop_front_internal();
    }
  }

  T &front() {
    HS_CHECK(!is_empty(), "front() on empty StaticCircularBuffer");
    return buffer[head];
  }

  const T &front() const {
    HS_CHECK(!is_empty(), "front() on empty StaticCircularBuffer");
    return buffer[head];
  }

  T &back() {
    HS_CHECK(!is_empty(), "back() on empty StaticCircularBuffer");
    return buffer[(head + count - 1) % N];
  }

  const T &back() const {
    HS_CHECK(!is_empty(), "back() on empty StaticCircularBuffer");
    return buffer[(head + count - 1) % N];
  }

  T &operator[](size_t index) {
    HS_CHECK(index < count, "StaticCircularBuffer index out of range");
    return buffer[(head + index) % N];
  }

  const T &operator[](size_t index) const {
    HS_CHECK(index < count, "StaticCircularBuffer index out of range");
    return buffer[(head + index) % N];
  }

  constexpr bool is_empty() const { return count == 0U; }
  constexpr bool is_full() const { return count == N; }
  constexpr size_t size() const { return count; }
  constexpr size_t capacity() const { return N; }

  iterator begin() { return iterator(this, 0); }
  iterator end() { return iterator(this, size()); }
  const_iterator begin() const { return const_iterator(this, 0); }
  const_iterator end() const { return const_iterator(this, size()); }

private:
  // Indices are uint32_t, not size_t, so pooled structs (Particle, VectorTrail,
  // OrientationTrail) have an identical layout on the 32-bit device and the
  // 64-bit native build — size_t would widen these three fields to 8 B on the
  // host and make per-effect arena footprints unrepresentative of the device
  // (see memory.h:15-22). On the device size_t IS uint32_t, so this is a
  // host-only narrowing with identical codegen and behavior on hardware.
  static_assert(N <= 0xFFFFFFFFu, "StaticCircularBuffer capacity exceeds uint32_t");
  std::array<T, N> buffer;
  uint32_t head;
  uint32_t tail;
  uint32_t count;

  void pop_back_internal() {
    tail = (tail + N - 1) % N; // + N first: never underflows (tail is uint32_t)
    count--;
  }

  void pop_front_internal() {
    head = (head + 1) % N;
    count--;
  }

  // True in-place construction for emplace_back/emplace_front. Every std::array
  // slot is a live object for the buffer's lifetime, so we end the existing
  // object's lifetime and construct the new value directly in its storage.
  // Unlike `buffer[idx] = T(args...)` this builds no temporary and requires only
  // that T be constructible from Args, not assignable. The placement-new result
  // is passed through std::launder so the returned reference is valid even for a
  // T that is not transparently replaceable (e.g. one with const/reference
  // members), where reusing the prior object's address is otherwise UB.
  template <typename... Args>
  T &construct_in_place(uint32_t idx, Args &&...args) {
    T *slot = &buffer[idx];
    slot->~T();
    return *std::launder(::new (static_cast<void *>(slot))
                             T(std::forward<Args>(args)...));
  }

  /// CRTP base providing all random-access iterator operators.
  /// Derived must expose: m_buffer, m_index, and typedefs.
  template <typename Derived, typename BufPtr, typename Ref, typename Ptr>
  class CircularIterBase {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = Ptr;
    using reference = Ref;

    BufPtr m_buffer;
    size_t m_index;

    CircularIterBase(BufPtr buffer, size_t idx)
        : m_buffer(buffer), m_index(idx) {}

    reference operator*() const { return (*m_buffer)[m_index]; }
    pointer operator->() const { return &(*m_buffer)[m_index]; }
    reference operator[](difference_type n) const {
      return (*m_buffer)[m_index + n];
    }

    Derived &operator++() { ++m_index; return self(); }
    Derived operator++(int) { Derived t = self(); ++m_index; return t; }
    Derived &operator--() { --m_index; return self(); }
    Derived operator--(int) { Derived t = self(); --m_index; return t; }
    Derived &operator+=(difference_type n) { m_index += n; return self(); }
    Derived &operator-=(difference_type n) { m_index -= n; return self(); }

    friend Derived operator+(const Derived &a, difference_type n) {
      return Derived(a.m_buffer, a.m_index + n);
    }
    friend Derived operator+(difference_type n, const Derived &a) {
      return Derived(a.m_buffer, a.m_index + n);
    }
    friend Derived operator-(const Derived &a, difference_type n) {
      return Derived(a.m_buffer, a.m_index - n);
    }
    friend difference_type operator-(const Derived &a, const Derived &b) {
      return a.m_index - b.m_index;
    }

    friend bool operator==(const Derived &a, const Derived &b) {
      return a.m_buffer == b.m_buffer && a.m_index == b.m_index;
    }
    friend bool operator!=(const Derived &a, const Derived &b) {
      return !(a == b);
    }
    friend bool operator<(const Derived &a, const Derived &b) {
      return a.m_index < b.m_index;
    }
    friend bool operator>(const Derived &a, const Derived &b) {
      return a.m_index > b.m_index;
    }
    friend bool operator<=(const Derived &a, const Derived &b) {
      return a.m_index <= b.m_index;
    }
    friend bool operator>=(const Derived &a, const Derived &b) {
      return a.m_index >= b.m_index;
    }

  private:
    Derived &self() { return static_cast<Derived &>(*this); }
  };

  class Iterator
      : public CircularIterBase<Iterator, StaticCircularBuffer *, T &, T *> {
    using Base = CircularIterBase<Iterator, StaticCircularBuffer *, T &, T *>;

  public:
    using Base::Base;
  };

  class ConstIterator
      : public CircularIterBase<ConstIterator, const StaticCircularBuffer *,
                                const T &, const T *> {
    using Base = CircularIterBase<ConstIterator, const StaticCircularBuffer *,
                                 const T &, const T *>;

  public:
    using Base::Base;
    ConstIterator(const iterator &other) : Base(other.m_buffer, other.m_index) {}
  };
};
