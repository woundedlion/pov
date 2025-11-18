#pragma once

#include <array>
#include <cstddef>     
#include <iterator>    
#include <stdexcept>   
#include <assert.h>
#include <initializer_list>
#include <utility>     
#include <optional>    

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
      assert(items.size() <= N && "Initializer list is larger than buffer capacity.");
    }
    for (const T& item : items) {
      push_back(item);
    }
  }

  void push_front(const T& item) {
    assert(!is_full() && "Buffer overflow: cannot push_front to a full buffer.");
    head = (head - 1 + N) % N;
    buffer[head] = item;
    count++;
  }

  void push_front(T&& item) {
    assert(!is_full() && "Buffer overflow: cannot push_front to a full buffer.");
    head = (head - 1 + N) % N;
    buffer[head] = std::move(item);
    count++;
  }

  template<typename... Args>
  T& emplace_front(Args&&... args) {
    assert(!is_full() && "Buffer overflow: cannot emplace_front in a full buffer.");
    head = (head - 1 + N) % N;
    T& emplaced_item = buffer[head].emplace(std::forward<Args>(args)...);
    count++;
    return emplaced_item;
  }

  void push_back(const T& item) {
    assert(!is_full() && "Buffer overflow: cannot push to a full buffer.");
    buffer[tail] = item;
    tail = (tail + 1) % N;
    count++;
  }

  void push_back(T&& item) {
    assert(!is_full() && "Buffer overflow: cannot push to a full buffer.");
    buffer[tail] = std::move(item); 
    tail = (tail + 1) % N;
    count++;
  }

  template<typename... Args>
  T& emplace_back(Args&&... args) {
    assert(!is_full() && "Buffer overflow: cannot emplace in a full buffer.");
    T& emplaced_item = buffer[tail].emplace(std::forward<Args>(args)...);
    tail = (tail + 1) % N;
    count++;
    return emplaced_item;
  }

  void pop_back() {
    assert(!is_empty() && "Buffer underflow: cannot pop_back from an empty buffer.");
    tail = (tail - 1 + N) % N;
    buffer[tail].reset();
    count--;
  }

  void pop() {
    assert(!is_empty() && "Buffer underflow: cannot pop from an empty buffer.");
    buffer[head].reset();
    head = (head + 1) % N;
    count--;
  }

  T& front() {
    assert(!is_empty() && "Buffer is empty.");
    return *buffer[head];
  }

  const T& front() const {
    assert(!is_empty() && "Buffer is empty.");
    return *buffer[head];
  }

  T& back() {
    assert(!is_empty() && "Buffer is empty.");
    return *buffer[(head + count - 1) % N];
  }

  const T& back() const {
    assert(!is_empty() && "Buffer is empty.");
    return *buffer[(head + count - 1) % N];
  }

  T& operator[](size_t index) {
    assert(index < count && "Index out of bounds.");
    return *buffer[(head + index) % N];
  }

  const T& operator[](size_t index) const {
    assert(index < count && "Index out of bounds.");
    return *buffer[(head + index) % N];
  }

  constexpr bool is_empty() const {
    return count == 0;
  }

  constexpr bool is_full() const {
    return count == N;
  }

  constexpr size_t size() const {
    return count;
  }

  constexpr size_t capacity() const {
    return N;
  }

  void clear() {
    while (!is_empty()) {
      pop();
    }
  }

  iterator begin() { return iterator(this, 0); }
  iterator end() { return iterator(this, size()); }

  const_iterator begin() const { return const_iterator(this, 0); }
  const_iterator end() const { return const_iterator(this, size()); }

  const_iterator cbegin() const { return const_iterator(this, 0); }
  const_iterator cend() const { return const_iterator(this, size()); }

private:

  std::array<std::optional<T>, N> buffer;
  size_t head;
  size_t tail;
  size_t count;

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
      : m_buffer(buffer), m_index(logical_index) {}

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
      : m_buffer(buffer), m_index(logical_index) {}

    const_iterator(const iterator& other)
      : m_buffer(other.m_buffer), m_index(other.m_index) {}

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