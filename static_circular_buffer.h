/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <cstddef>     
#include <iterator>    
#include <utility>     
#include <Arduino.h> 

/**
 * @brief A fixed-size circular buffer optimized for stability.
 * @details
 * - No dynamic memory allocation (prevents heap fragmentation).
 * - No asserts (prevents hard crashes in live shows).
 * - "Delete Oldest" strategy on overflow (never stops processing).
 * - "Return Zeroed Dummy" on invalid access (prevents read crashes).
 * - Supports types WITHOUT default constructors.
 */
template <typename T, size_t N>
class StaticCircularBuffer {
  class iterator;
  class const_iterator;

public:
  using value_type = T;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using reference = T&;
  using const_reference = const T&;

  StaticCircularBuffer() :
    head(0), tail(0), count(0)
  {
  }

  StaticCircularBuffer(std::initializer_list<T> items)
    : head(0), tail(0), count(0) {
    if (items.size() > N) {
      Serial.println("Buffer Error: Initializer list exceeds capacity. Truncating.");
    }
    for (const T& item : items) {
      push_back(item);
    }
  }


  void push_front(const T& item) {
    if (is_full()) {
      // Drop the tail (oldest in this context) to make room
      pop_back_internal();
    }
    head = (head - 1 + N) % N;
    buffer[head] = item;
    count++;
  }

  void push_front(T&& item) {
    if (is_full()) {
      pop_back_internal();
    }
    head = (head - 1 + N) % N;
    buffer[head] = std::move(item);
    count++;
  }

  template<typename... Args>
  T& emplace_front(Args&&... args) {
    if (is_full()) {
      pop_back_internal();
    }
    head = (head - 1 + N) % N;
    buffer[head] = T(std::forward<Args>(args)...);
    count++;
    return buffer[head];
  }

  void push_back(const T& item) {
    if (is_full()) {
      // Drop the head (oldest in this context) to make room
      pop_front_internal();
    }
    buffer[tail] = item;
    tail = (tail + 1) % N;
    count++;
  }

  void push_back(T&& item) {
    if (is_full()) {
      pop_front_internal();
    }
    buffer[tail] = std::move(item);
    tail = (tail + 1) % N;
    count++;
  }

  template<typename... Args>
  T& emplace_back(Args&&... args) {
    if (is_full()) {
      pop_front_internal();
    }
    buffer[tail] = T(std::forward<Args>(args)...);
    T& emplaced_item = buffer[tail];
    tail = (tail + 1) % N;
    count++;
    return emplaced_item;
  }

  void pop_back() {
    if (is_empty()) return;
    pop_back_internal();
  }

  void pop() {
    if (is_empty()) return;
    pop_front_internal();
  }

  void clear() {
    while (!is_empty()) {
      pop_front_internal();
    }
  }

  T& front() {
    if (is_empty()) return get_dummy();
    return buffer[head];
  }

  const T& front() const {
    if (is_empty()) return get_dummy();
    return buffer[head];
  }

  T& back() {
    if (is_empty()) return get_dummy();
    return buffer[(head + count - 1) % N];
  }

  const T& back() const {
    if (is_empty()) return get_dummy();
    return buffer[(head + count - 1) % N];
  }

  T& operator[](size_t index) {
    if (index >= count) {
      return get_dummy();
    }
    return buffer[(head + index) % N];
  }

  const T& operator[](size_t index) const {
    if (index >= count) {
      return get_dummy();
    }
    return buffer[(head + index) % N];
  }

  constexpr bool is_empty() const { return count == 0; }
  constexpr bool is_full() const { return count == N; }
  constexpr size_t size() const { return count; }
  constexpr size_t capacity() const { return N; }

  iterator begin() { return iterator(this, 0); }
  iterator end() { return iterator(this, size()); }
  const_iterator begin() const { return const_iterator(this, 0); }
  const_iterator end() const { return const_iterator(this, size()); }

private:

  std::array<T, N> buffer;
  size_t head;
  size_t tail;
  size_t count;

  static T& get_dummy() {
    static struct {
      alignas(T) uint8_t bytes[sizeof(T)];
    } raw = { 0 };
    return *reinterpret_cast<T*>(raw.bytes);
  }

  void pop_back_internal() {
    tail = (tail - 1 + N) % N;
    count--;
  }

  void pop_front_internal() {
    head = (head + 1) % N;
    count--;
  }

  class iterator {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;

    StaticCircularBuffer* m_buffer;
    size_t m_index;

    iterator(StaticCircularBuffer* buffer, size_t logical_index)
      : m_buffer(buffer), m_index(logical_index) {
    }

    reference operator*() const { return (*m_buffer)[m_index]; }
    pointer   operator->() const { return &(*m_buffer)[m_index]; }
    reference operator[](difference_type n) const { return (*m_buffer)[m_index + n]; }

    iterator& operator++() { ++m_index; return *this; }
    iterator  operator++(int) { iterator tmp = *this; ++(*this); return tmp; }
    iterator& operator--() { --m_index; return *this; }
    iterator  operator--(int) { iterator tmp = *this; --(*this); return tmp; }
    iterator& operator+=(difference_type n) { m_index += n; return *this; }
    iterator& operator-=(difference_type n) { m_index -= n; return *this; }

    friend iterator operator+(const iterator& a, difference_type n) { return iterator(a.m_buffer, a.m_index + n); }
    friend iterator operator+(difference_type n, const iterator& a) { return iterator(a.m_buffer, a.m_index + n); }
    friend iterator operator-(const iterator& a, difference_type n) { return iterator(a.m_buffer, a.m_index - n); }
    friend difference_type operator-(const iterator& a, const iterator& b) { return a.m_index - b.m_index; }

    friend bool operator==(const iterator& a, const iterator& b) { return a.m_buffer == b.m_buffer && a.m_index == b.m_index; }
    friend bool operator!=(const iterator& a, const iterator& b) { return !(a == b); }
    friend bool operator<(const iterator& a, const iterator& b) { return a.m_index < b.m_index; }
    friend bool operator>(const iterator& a, const iterator& b) { return a.m_index > b.m_index; }
    friend bool operator<=(const iterator& a, const iterator& b) { return a.m_index <= b.m_index; }
    friend bool operator>=(const iterator& a, const iterator& b) { return a.m_index >= b.m_index; }
  };

  class const_iterator {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = const T*;
    using reference = const T&;

    const StaticCircularBuffer* m_buffer;
    size_t m_index;

    const_iterator(const StaticCircularBuffer* buffer, size_t logical_index)
      : m_buffer(buffer), m_index(logical_index) {
    }

    const_iterator(const iterator& other)
      : m_buffer(other.m_buffer), m_index(other.m_index) {
    }

    reference operator*() const { return (*m_buffer)[m_index]; }
    pointer   operator->() const { return &(*m_buffer)[m_index]; }
    reference operator[](difference_type n) const { return (*m_buffer)[m_index + n]; }

    const_iterator& operator++() { ++m_index; return *this; }
    const_iterator  operator++(int) { const_iterator tmp = *this; ++(*this); return tmp; }
    const_iterator& operator--() { --m_index; return *this; }
    const_iterator  operator--(int) { const_iterator tmp = *this; --(*this); return tmp; }
    const_iterator& operator+=(difference_type n) { m_index += n; return *this; }
    const_iterator& operator-=(difference_type n) { m_index -= n; return *this; }

    friend const_iterator operator+(const const_iterator& a, difference_type n) { return const_iterator(a.m_buffer, a.m_index + n); }
    friend const_iterator operator+(difference_type n, const const_iterator& a) { return const_iterator(a.m_buffer, a.m_index + n); }
    friend const_iterator operator-(const const_iterator& a, difference_type n) { return const_iterator(a.m_buffer, a.m_index - n); }
    friend difference_type operator-(const const_iterator& a, const const_iterator& b) { return a.m_index - b.m_index; }

    friend bool operator==(const const_iterator& a, const const_iterator& b) { return a.m_buffer == b.m_buffer && a.m_index == b.m_index; }
    friend bool operator!=(const const_iterator& a, const const_iterator& b) { return !(a == b); }
    friend bool operator<(const const_iterator& a, const const_iterator& b) { return a.m_index < b.m_index; }
    friend bool operator>(const const_iterator& a, const const_iterator& b) { return a.m_index > b.m_index; }
    friend bool operator<=(const const_iterator& a, const const_iterator& b) { return a.m_index <= b.m_index; }
    friend bool operator>=(const const_iterator& a, const const_iterator& b) { return a.m_index >= b.m_index; }
  };
};