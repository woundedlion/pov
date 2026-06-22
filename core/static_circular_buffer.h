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

  /** @brief Compile-time fixed capacity, for static_asserts binding two buffers'
   *  sizes (e.g. scan_region's seam-split `norm` == 2x its input span buffer). */
  static constexpr size_t kCapacity = N;

  /**
   * @brief Constructs an empty buffer.
   */
  StaticCircularBuffer() : head(0), tail(0), count(0) {}

  /**
   * @brief Constructs a buffer from an initializer list, filling in order.
   * @param items Elements to insert, front-to-back.
   * @details On overflow keeps the LAST N items, since push_back evicts the
   * oldest element when full.
   */
  StaticCircularBuffer(std::initializer_list<T> items)
      : head(0), tail(0), count(0) {
    if (items.size() > N) {
      // push_back evicts the oldest on overflow, so the buffer keeps the LAST
      // N items of the list.
      hs::log("Buffer Warning: Initializer list exceeds capacity. Keeping last "
              "N items.");
    }
    for (const T &item : items) {
      push_back(item);
    }
  }

  /**
   * @brief Inserts an element at the front by copy.
   * @param item Element to copy into the new front slot.
   * @details When full, evicts the back element (the oldest for front-pushes).
   */
  void push_front(const T &item) {
    if (is_full()) {
      pop_back_internal();
    }
    head = (head + N - 1) % N; // + N first: never underflows (head is uint32_t)
    buffer[head] = item;
    count++;
  }

  /**
   * @brief Inserts an element at the front by move.
   * @param item Element to move into the new front slot.
   * @details When full, evicts the back element (the oldest for front-pushes).
   */
  void push_front(T &&item) {
    if (is_full()) {
      pop_back_internal();
    }
    head = (head + N - 1) % N; // + N first: never underflows (head is uint32_t)
    buffer[head] = std::move(item);
    count++;
  }

  /**
   * @brief Constructs an element in place at the front.
   * @tparam Args Constructor argument types for T.
   * @param args Arguments forwarded to T's constructor.
   * @return Reference to the newly constructed front element.
   * @details Evicts the back element when full. head/count are committed only
   * after the (potentially throwing) constructor succeeds, so a throwing T
   * leaves the buffer's size invariant intact.
   */
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

  /**
   * @brief Inserts an element at the back by copy.
   * @param item Element to copy into the new back slot.
   * @details When full, evicts the front element (the oldest for back-pushes).
   */
  void push_back(const T &item) {
    if (is_full()) {
      pop_front_internal();
    }
    buffer[tail] = item;
    tail = (tail + 1) % N;
    count++;
  }

  /**
   * @brief Inserts an element at the back by move.
   * @param item Element to move into the new back slot.
   * @details When full, evicts the front element (the oldest for back-pushes).
   */
  void push_back(T &&item) {
    if (is_full()) {
      pop_front_internal();
    }
    buffer[tail] = std::move(item);
    tail = (tail + 1) % N;
    count++;
  }

  /**
   * @brief Constructs an element in place at the back.
   * @tparam Args Constructor argument types for T.
   * @param args Arguments forwarded to T's constructor.
   * @return Reference to the newly constructed back element.
   * @details Evicts the front element when full. tail/count are committed only
   * after the (potentially throwing) constructor succeeds, so a throwing T
   * leaves the buffer's size invariant intact.
   */
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

  /**
   * @brief Removes the back element.
   * @details No-op when the buffer is empty.
   */
  void pop_back() {
    if (is_empty())
      return;
    pop_back_internal();
  }

  /**
   * @brief Removes the front element.
   * @details No-op when the buffer is empty.
   */
  void pop() {
    if (is_empty())
      return;
    pop_front_internal();
  }

  /**
   * @brief Removes all elements, running each element's destructor.
   */
  void clear() {
    while (!is_empty()) {
      pop_front_internal();
    }
  }

  /**
   * @brief Returns the front (oldest) element.
   * @return Reference to the front element.
   * @details Traps via HS_CHECK if the buffer is empty.
   */
  T &front() {
    HS_CHECK(!is_empty(), "front() on empty StaticCircularBuffer");
    return buffer[head];
  }

  /**
   * @brief Returns the front (oldest) element.
   * @return Const reference to the front element.
   * @details Traps via HS_CHECK if the buffer is empty.
   */
  const T &front() const {
    HS_CHECK(!is_empty(), "front() on empty StaticCircularBuffer");
    return buffer[head];
  }

  /**
   * @brief Returns the back (newest) element.
   * @return Reference to the back element.
   * @details Traps via HS_CHECK if the buffer is empty.
   */
  T &back() {
    HS_CHECK(!is_empty(), "back() on empty StaticCircularBuffer");
    return buffer[(head + count - 1) % N];
  }

  /**
   * @brief Returns the back (newest) element.
   * @return Const reference to the back element.
   * @details Traps via HS_CHECK if the buffer is empty.
   */
  const T &back() const {
    HS_CHECK(!is_empty(), "back() on empty StaticCircularBuffer");
    return buffer[(head + count - 1) % N];
  }

  /**
   * @brief Accesses the element at a logical index, front-to-back.
   * @param index Position in [0, size()), measured from the front.
   * @return Reference to the indexed element.
   * @details Traps via HS_CHECK if index is out of range.
   */
  T &operator[](size_t index) {
    HS_CHECK(index < count, "StaticCircularBuffer index out of range");
    return buffer[(head + index) % N];
  }

  /**
   * @brief Accesses the element at a logical index, front-to-back.
   * @param index Position in [0, size()), measured from the front.
   * @return Const reference to the indexed element.
   * @details Traps via HS_CHECK if index is out of range.
   */
  const T &operator[](size_t index) const {
    HS_CHECK(index < count, "StaticCircularBuffer index out of range");
    return buffer[(head + index) % N];
  }

  /**
   * @brief Reports whether the buffer holds no elements.
   * @return True if the buffer is empty.
   */
  constexpr bool is_empty() const { return count == 0U; }

  /**
   * @brief Reports whether the buffer is at capacity.
   * @return True if the buffer holds N elements.
   */
  constexpr bool is_full() const { return count == N; }

  /**
   * @brief Returns the number of stored elements.
   * @return Current element count.
   */
  constexpr size_t size() const { return count; }

  /**
   * @brief Returns the fixed capacity.
   * @return Maximum number of elements N.
   */
  constexpr size_t capacity() const { return N; }

  /**
   * @brief Returns an iterator to the front element.
   * @return Mutable iterator at the front.
   */
  iterator begin() { return iterator(this, 0); }

  /**
   * @brief Returns an iterator past the back element.
   * @return Mutable iterator at one-past-the-end.
   */
  iterator end() { return iterator(this, size()); }

  /**
   * @brief Returns a const iterator to the front element.
   * @return Const iterator at the front.
   */
  const_iterator begin() const { return const_iterator(this, 0); }

  /**
   * @brief Returns a const iterator past the back element.
   * @return Const iterator at one-past-the-end.
   */
  const_iterator end() const { return const_iterator(this, size()); }

private:
  // Indices are uint32_t, not size_t, so pooled structs (Particle, VectorTrail,
  // OrientationTrail) have an identical layout on the 32-bit device and the
  // 64-bit native build — size_t would widen these three fields to 8 B on the
  // host and make per-effect arena footprints unrepresentative of the device
  // (see memory.h:15-22). On the device size_t IS uint32_t, so this is a
  // host-only narrowing with identical codegen and behavior on hardware.
  static_assert(N <= 0xFFFFFFFFu, "StaticCircularBuffer capacity exceeds uint32_t");
  // N == 0 would make every `% N` index update a division-by-zero (UB); all
  // instantiations use N >= 1, so trap a zero-capacity buffer at compile time.
  static_assert(N > 0, "StaticCircularBuffer requires N >= 1");
  std::array<T, N> buffer; /**< Backing storage; every slot is a live object. */
  uint32_t head;           /**< Index of the front (oldest) element. */
  uint32_t tail;           /**< Index of the next free back slot. */
  uint32_t count;          /**< Number of elements currently stored. */

  /**
   * @brief Removes the back element without an emptiness check.
   * @details Caller must ensure the buffer is non-empty.
   */
  void pop_back_internal() {
    tail = (tail + N - 1) % N; // + N first: never underflows (tail is uint32_t)
    count--;
  }

  /**
   * @brief Removes the front element without an emptiness check.
   * @details Caller must ensure the buffer is non-empty.
   */
  void pop_front_internal() {
    head = (head + 1) % N;
    count--;
  }

  /**
   * @brief Constructs a T directly in the slot at idx, in place.
   * @tparam Args Constructor argument types for T.
   * @param idx Slot index in [0, N) whose object is replaced.
   * @param args Arguments forwarded to T's constructor.
   * @return Reference to the newly constructed element.
   * @details Every std::array slot is a live object for the buffer's lifetime,
   * so the existing object's lifetime is ended and the new value is constructed
   * directly in its storage. Unlike `buffer[idx] = T(args...)` this builds no
   * temporary and requires only that T be constructible from Args, not
   * assignable. The placement-new result is passed through std::launder so the
   * returned reference is valid even for a T that is not transparently
   * replaceable (e.g. one with const/reference members), where reusing the prior
   * object's address is otherwise UB.
   */
  template <typename... Args>
  T &construct_in_place(uint32_t idx, Args &&...args) {
    T *slot = &buffer[idx];
    slot->~T();
    return *std::launder(::new (static_cast<void *>(slot))
                             T(std::forward<Args>(args)...));
  }

  /**
   * @brief CRTP base providing all random-access iterator operators.
   * @tparam Derived Concrete iterator type; must expose m_buffer, m_index.
   * @tparam BufPtr Pointer-to-buffer type (mutable or const).
   * @tparam Ref Reference type yielded on dereference.
   * @tparam Ptr Pointer type yielded by operator->.
   * @details Derived must expose m_buffer, m_index, and the iterator typedefs.
   */
  template <typename Derived, typename BufPtr, typename Ref, typename Ptr>
  class CircularIterBase {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = Ptr;
    using reference = Ref;

    BufPtr m_buffer; /**< Buffer this iterator traverses. */
    size_t m_index;  /**< Logical position, front-to-back. */

    /**
     * @brief Constructs an iterator over a buffer at a logical position.
     * @param buffer Buffer to traverse.
     * @param idx Logical index, measured from the front.
     */
    CircularIterBase(BufPtr buffer, size_t idx)
        : m_buffer(buffer), m_index(idx) {}

    /**
     * @brief Dereferences the iterator.
     * @return Reference to the element at the current position.
     */
    reference operator*() const { return (*m_buffer)[m_index]; }

    /**
     * @brief Member access through the iterator.
     * @return Pointer to the element at the current position.
     */
    pointer operator->() const { return &(*m_buffer)[m_index]; }

    /**
     * @brief Accesses the element offset n positions from the current one.
     * @param n Offset from the current position.
     * @return Reference to the element at m_index + n.
     */
    reference operator[](difference_type n) const {
      return (*m_buffer)[m_index + n];
    }

    /**
     * @brief Pre-increment; advances toward the back.
     * @return Reference to this iterator after advancing.
     */
    Derived &operator++() { ++m_index; return self(); }

    /**
     * @brief Post-increment; advances toward the back.
     * @param int Unused disambiguation tag for post-increment.
     * @return Copy of the iterator before advancing.
     */
    Derived operator++(int) { Derived t = self(); ++m_index; return t; }

    /**
     * @brief Pre-decrement; moves toward the front.
     * @return Reference to this iterator after moving.
     */
    Derived &operator--() { --m_index; return self(); }

    /**
     * @brief Post-decrement; moves toward the front.
     * @param int Unused disambiguation tag for post-decrement.
     * @return Copy of the iterator before moving.
     */
    Derived operator--(int) { Derived t = self(); --m_index; return t; }

    /**
     * @brief Advances the iterator by n positions.
     * @param n Number of positions to advance (may be negative).
     * @return Reference to this iterator after advancing.
     */
    Derived &operator+=(difference_type n) { m_index += n; return self(); }

    /**
     * @brief Rewinds the iterator by n positions.
     * @param n Number of positions to rewind (may be negative).
     * @return Reference to this iterator after rewinding.
     */
    Derived &operator-=(difference_type n) { m_index -= n; return self(); }

    /**
     * @brief Returns an iterator advanced n positions from a.
     * @param a Base iterator.
     * @param n Number of positions to advance.
     * @return Iterator at a.m_index + n.
     */
    friend Derived operator+(const Derived &a, difference_type n) {
      return Derived(a.m_buffer, a.m_index + n);
    }

    /**
     * @brief Returns an iterator advanced n positions from a.
     * @param n Number of positions to advance.
     * @param a Base iterator.
     * @return Iterator at a.m_index + n.
     */
    friend Derived operator+(difference_type n, const Derived &a) {
      return Derived(a.m_buffer, a.m_index + n);
    }

    /**
     * @brief Returns an iterator rewound n positions from a.
     * @param a Base iterator.
     * @param n Number of positions to rewind.
     * @return Iterator at a.m_index - n.
     */
    friend Derived operator-(const Derived &a, difference_type n) {
      return Derived(a.m_buffer, a.m_index - n);
    }

    /**
     * @brief Computes the distance between two iterators.
     * @param a Left iterator.
     * @param b Right iterator.
     * @return Signed number of positions from b to a.
     */
    friend difference_type operator-(const Derived &a, const Derived &b) {
      return (difference_type)a.m_index - (difference_type)b.m_index;
    }

    /**
     * @brief Tests two iterators for equality.
     * @param a Left iterator.
     * @param b Right iterator.
     * @return True if both reference the same buffer and position.
     */
    friend bool operator==(const Derived &a, const Derived &b) {
      return a.m_buffer == b.m_buffer && a.m_index == b.m_index;
    }

    /**
     * @brief Tests two iterators for inequality.
     * @param a Left iterator.
     * @param b Right iterator.
     * @return True if the iterators are not equal.
     */
    friend bool operator!=(const Derived &a, const Derived &b) {
      return !(a == b);
    }

    /**
     * @brief Orders two iterators by position.
     * @param a Left iterator.
     * @param b Right iterator.
     * @return True if a precedes b.
     */
    friend bool operator<(const Derived &a, const Derived &b) {
      return a.m_index < b.m_index;
    }

    /**
     * @brief Orders two iterators by position.
     * @param a Left iterator.
     * @param b Right iterator.
     * @return True if a follows b.
     */
    friend bool operator>(const Derived &a, const Derived &b) {
      return a.m_index > b.m_index;
    }

    /**
     * @brief Orders two iterators by position.
     * @param a Left iterator.
     * @param b Right iterator.
     * @return True if a does not follow b.
     */
    friend bool operator<=(const Derived &a, const Derived &b) {
      return a.m_index <= b.m_index;
    }

    /**
     * @brief Orders two iterators by position.
     * @param a Left iterator.
     * @param b Right iterator.
     * @return True if a does not precede b.
     */
    friend bool operator>=(const Derived &a, const Derived &b) {
      return a.m_index >= b.m_index;
    }

  private:
    /**
     * @brief Downcasts to the derived iterator type.
     * @return Reference to *this as Derived.
     */
    Derived &self() { return static_cast<Derived &>(*this); }
  };

  /**
   * @brief Mutable random-access iterator over the buffer, front-to-back.
   */
  class Iterator
      : public CircularIterBase<Iterator, StaticCircularBuffer *, T &, T *> {
    using Base = CircularIterBase<Iterator, StaticCircularBuffer *, T &, T *>;

  public:
    using Base::Base;
  };

  /**
   * @brief Read-only counterpart of Iterator; implicitly built from one.
   */
  class ConstIterator
      : public CircularIterBase<ConstIterator, const StaticCircularBuffer *,
                                const T &, const T *> {
    using Base = CircularIterBase<ConstIterator, const StaticCircularBuffer *,
                                 const T &, const T *>;

  public:
    using Base::Base;
    /**
     * @brief Converts a mutable iterator into a const iterator.
     * @param other Mutable iterator to copy position and buffer from.
     */
    ConstIterator(const iterator &other) : Base(other.m_buffer, other.m_index) {}
  };
};
