/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file static_circular_buffer.h
 * @brief StaticCircularBuffer: fixed-capacity, allocation-free ring buffer.
 */

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <new>
#include <type_traits>
#include <utility>
#include "platform/platform.h"

/**
 * @brief A fixed-size circular buffer optimized for stability.
 * @tparam T Element type. All N backing slots hold live objects from
 * construction onward, so T must be default-constructible. push_front/push_back
 * assign into a slot and so additionally require copy/move assignment;
 * emplace_front/emplace_back construct in place and require only that T be
 * constructible from the forwarded arguments, so they are the path for a
 * non-assignable T.
 * @tparam N Capacity in elements; must be >= 1.
 * @details No dynamic allocation. Overflow evicts the oldest element.
 * front()/back() on an empty buffer and out-of-range operator[] HS_CHECK-trap.
 * Models Container and ReversibleContainer; not a SequenceContainer (no
 * insert/erase/assign) or ContiguousContainer (the live run wraps, so there is
 * no data()).
 */
template <typename T, size_t N> class StaticCircularBuffer {
  // Private; named publicly through the iterator usings.
  class Iterator;
  class ConstIterator;

public:
  using iterator = Iterator;
  using const_iterator = ConstIterator;
  using reverse_iterator = std::reverse_iterator<Iterator>;
  using const_reverse_iterator = std::reverse_iterator<ConstIterator>;
  using value_type = T;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using reference = T &;
  using const_reference = const T &;
  using pointer = T *;
  using const_pointer = const T *;

  /** @brief Compile-time fixed capacity. */
  static constexpr size_t CAPACITY = N;

  // back() and the iterators form `head + count - 1` (up to 2N-2) before the `% N`
  // fold; head/count are uint32_t, so cap N to keep that sum from overflowing 2^32.
  static_assert(N <= (UINT32_MAX - 1) / 2,
                "StaticCircularBuffer N too large: head+count could overflow");

  /**
   * @brief Constructs an empty buffer.
   */
  StaticCircularBuffer() : head(0), tail(0), count(0) {}

  /**
   * @brief Constructs a buffer from a braced element list, filling front-to-back.
   * @tparam Args Element types, each constructible into T.
   * @param args Elements to insert, in order.
   * @details An over-capacity list fails a static_assert. Elements are built
   * with T{...}, so narrowing conversions are diagnosed.
   */
  template <
      typename... Args,
      typename = std::enable_if_t<
          (sizeof...(Args) > 0) &&
          (!std::is_same_v<StaticCircularBuffer, std::remove_cvref_t<Args>> &&
           ...) &&
          (std::is_constructible_v<T, Args &&> && ...)>>
  explicit StaticCircularBuffer(Args &&...args) : head(0), tail(0), count(0) {
    static_assert(sizeof...(Args) <= N,
                  "StaticCircularBuffer initializer list exceeds capacity N");
    (emplace_back(T{std::forward<Args>(args)}), ...);
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
   * after the constructor succeeds, but the slot's old object is destroyed
   * first, so a throwing T is unsupported (see construct_in_place).
   */
  template <typename... Args> T &emplace_front(Args &&...args) {
    if (is_full()) {
      pop_back_internal();
    }
    uint32_t slot = (head + N - 1) % N; // + N first: never underflows
    T &ref = construct_in_place(slot, std::forward<Args>(args)...);
    head = slot;
    count++;
    return ref;
  }

  /** @brief Appends without eviction; returns false when full. */
  [[nodiscard]] bool try_push_back(const T &item) {
    if (is_full())
      return false;
    push_back(item);
    return true;
  }

  /** @brief Appends without eviction; a full buffer leaves item untouched. */
  [[nodiscard]] bool try_push_back(T &&item) {
    if (is_full())
      return false;
    push_back(std::move(item));
    return true;
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
   * after the constructor succeeds, but the slot's old object is destroyed
   * first, so a throwing T is unsupported (see construct_in_place).
   */
  template <typename... Args> T &emplace_back(Args &&...args) {
    if (is_full()) {
      pop_front_internal();
    }
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
  void pop_front() {
    if (is_empty())
      return;
    pop_front_internal();
  }

  /**
   * @brief Empties the buffer for reuse.
   * @details Resets head == tail == count == 0, so a caller indexing raw
   * &buf[0] linearly starts at 0; backing slots stay live (no per-element
   * destructor runs).
   */
  void clear() { head = tail = count = 0; }

  /**
   * @brief Returns the front element.
   * @return Reference to the front element.
   * @details Traps via HS_CHECK if the buffer is empty.
   */
  T &front() {
    HS_CHECK(!is_empty(), "front() on empty StaticCircularBuffer");
    return buffer[head];
  }

  /**
   * @brief Returns the front element.
   * @return Const reference to the front element.
   * @details Traps via HS_CHECK if the buffer is empty.
   */
  const T &front() const {
    HS_CHECK(!is_empty(), "front() on empty StaticCircularBuffer");
    return buffer[head];
  }

  /**
   * @brief Returns the back element.
   * @return Reference to the back element.
   * @details Traps via HS_CHECK if the buffer is empty.
   */
  T &back() {
    HS_CHECK(!is_empty(), "back() on empty StaticCircularBuffer");
    return buffer[(head + count - 1) % N];
  }

  /**
   * @brief Returns the back element.
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
    HS_CHECK(index < count,
             "StaticCircularBuffer::operator[]: index %lu is outside [0, %lu)",
             static_cast<unsigned long>(index),
             static_cast<unsigned long>(count));
    return buffer[(head + index) % N];
  }

  /**
   * @brief Accesses the element at a logical index, front-to-back.
   * @param index Position in [0, size()), measured from the front.
   * @return Const reference to the indexed element.
   * @details Traps via HS_CHECK if index is out of range.
   */
  const T &operator[](size_t index) const {
    HS_CHECK(index < count,
             "StaticCircularBuffer::operator[]: index %lu is outside [0, %lu)",
             static_cast<unsigned long>(index),
             static_cast<unsigned long>(count));
    return buffer[(head + index) % N];
  }

  /**
   * @brief Reports whether the buffer holds no elements.
   * @return True if the buffer is empty.
   */
  constexpr bool is_empty() const { return count == 0U; }

  /**
   * @brief Reports whether the buffer holds no elements.
   * @return True if the buffer is empty.
   * @details Container-requirement spelling of is_empty().
   */
  constexpr bool empty() const { return count == 0U; }

  /**
   * @brief Reports whether the buffer is at capacity.
   * @return True if the buffer holds N elements.
   */
  constexpr bool is_full() const { return count == N; }

  /**
   * @brief Reports whether logical index 0 maps to backing slot 0.
   * @return True when head == 0, so operator[] visits the live elements in
   * backing-slot order and linear indexing matches the logical order.
   */
  constexpr bool is_linear() const { return head == 0U; }

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
   * @brief Returns the largest number of elements the buffer can ever hold.
   * @return Maximum number of elements N.
   * @details Container-requirement spelling of capacity(); fixed, so the two agree.
   */
  constexpr size_t max_size() const { return N; }

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

  /**
   * @brief Returns a const iterator to the front element.
   * @return Const iterator at the front.
   */
  const_iterator cbegin() const { return const_iterator(this, 0); }

  /**
   * @brief Returns a const iterator past the back element.
   * @return Const iterator at one-past-the-end.
   */
  const_iterator cend() const { return const_iterator(this, size()); }

  /**
   * @brief Returns a reverse iterator to the back element.
   * @return Mutable reverse iterator at the back.
   */
  reverse_iterator rbegin() { return reverse_iterator(end()); }

  /**
   * @brief Returns a reverse iterator past the front element.
   * @return Mutable reverse iterator at one-before-the-front.
   */
  reverse_iterator rend() { return reverse_iterator(begin()); }

  /**
   * @brief Returns a const reverse iterator to the back element.
   * @return Const reverse iterator at the back.
   */
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }

  /**
   * @brief Returns a const reverse iterator past the front element.
   * @return Const reverse iterator at one-before-the-front.
   */
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }

  /**
   * @brief Returns a const reverse iterator to the back element.
   * @return Const reverse iterator at the back.
   */
  const_reverse_iterator crbegin() const {
    return const_reverse_iterator(end());
  }

  /**
   * @brief Returns a const reverse iterator past the front element.
   * @return Const reverse iterator at one-before-the-front.
   */
  const_reverse_iterator crend() const {
    return const_reverse_iterator(begin());
  }

  /**
   * @brief Exchanges contents with another buffer of the same type.
   * @param other Buffer to swap with.
   * @details Trivially copyable slots exchange their object representations;
   * other slots use elementwise swap. No allocation.
   */
  void swap(StaticCircularBuffer &other) {
    if (this == &other)
      return;
    if constexpr (std::is_trivially_copyable_v<T>) {
      for (size_t i = 0; i < N; ++i) {
        unsigned char temp[sizeof(T)];
        std::memcpy(temp, &buffer[i], sizeof(T));
        std::memcpy(&buffer[i], &other.buffer[i], sizeof(T));
        std::memcpy(&other.buffer[i], temp, sizeof(T));
      }
    } else {
      buffer.swap(other.buffer);
    }
    std::swap(head, other.head);
    std::swap(tail, other.tail);
    std::swap(count, other.count);
  }

  /**
   * @brief Exchanges the contents of two buffers.
   * @param a First buffer.
   * @param b Second buffer.
   * @details Hidden friend, so unqualified swap(a, b) finds it by ADL.
   */
  friend void swap(StaticCircularBuffer &a, StaticCircularBuffer &b) {
    a.swap(b);
  }

  /**
   * @brief Compares two buffers elementwise, front to back.
   * @param a Left buffer.
   * @param b Right buffer.
   * @return True if both hold the same number of elements and every logical
   * position compares equal.
   * @details Hidden friend defined inline, so T is only required to be
   * equality-comparable when this is actually called. Dead slots outside the live
   * run are ignored, so two buffers holding equal elements at different head
   * offsets still compare equal.
   */
  friend bool operator==(const StaticCircularBuffer &a,
                         const StaticCircularBuffer &b) {
    if (a.count != b.count)
      return false;
    for (uint32_t i = 0; i < a.count; ++i) {
      if (!(a.buffer[(a.head + i) % N] == b.buffer[(b.head + i) % N]))
        return false;
    }
    return true;
  }

  /**
   * @brief Negation of operator==.
   * @param a Left buffer.
   * @param b Right buffer.
   * @return True if the buffers differ in size or in any element.
   */
  friend bool operator!=(const StaticCircularBuffer &a,
                         const StaticCircularBuffer &b) {
    return !(a == b);
  }

  /**
   * @brief Visits live elements from front to back.
   * @tparam Fn Callable accepting each element and its logical index.
   * @param fn Callable invoked once per live element.
   */
  template <typename Fn> void for_each(Fn &&fn) const {
    uint32_t index = head;
    for (uint32_t visited = 0; visited < count; ++visited) {
      fn(buffer[index], visited);
      if (++index == N)
        index = 0;
    }
  }

private:
  // N == 0 would make every `% N` index update a division by zero (UB).
  static_assert(N > 0, "StaticCircularBuffer requires N >= 1");
  std::array<T, N> buffer; /**< Backing storage; every slot is a live object. */
  // uint32_t indices keep the layout identical on the 32-bit device and 64-bit
  // host.
  uint32_t head;  /**< Index of the front element. */
  uint32_t tail;  /**< Index of the next free back slot. */
  uint32_t count; /**< Number of elements currently stored. */

  /**
   * @brief Removes the back element without an emptiness check.
   * @details The buffer must be non-empty; otherwise count wraps past 0.
   */
  void pop_back_internal() {
    tail = (tail + N - 1) % N; // + N first: never underflows (tail is uint32_t)
    count--;
  }

  /**
   * @brief Removes the front element without an emptiness check.
   * @details The buffer must be non-empty; otherwise count wraps past 0.
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
   * @details Ends the existing object's lifetime and constructs the new value
   * directly in its storage; T need not be assignable. std::launder covers only
   * the returned reference; other accessors read through the un-laundered
   * `buffer` array, so T must be transparently replaceable.
   * @warning The old object is destroyed first, so a throwing constructor leaves
   * a dead slot that a later call re-destroys (UB); a throwing T is unsupported.
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
   * @tparam Derived Concrete iterator type, publicly derived from this base.
   * @tparam BufPtr Pointer-to-buffer type (mutable or const).
   * @tparam Ref Reference type yielded on dereference.
   * @tparam Ptr Pointer type yielded by operator->.
   * @details Holds the whole iterator state and every operator; Derived adds
   * only its constructors. ConstIterator is a friend so its converting
   * constructor can read an Iterator's state.
   */
  template <typename Derived, typename BufPtr, typename Ref, typename Ptr>
  class CircularIterBase {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = Ptr;
    using reference = Ref;

    /**
     * @brief Constructs a singular iterator, satisfying std::semiregular so the
     *        type models std::random_access_iterator.
     */
    CircularIterBase() : buffer(nullptr), index(0) {}

    /**
     * @brief Constructs an iterator over a buffer at a logical position.
     * @param buffer Buffer to traverse.
     * @param idx Logical index, measured from the front.
     */
    CircularIterBase(BufPtr buffer, size_t idx) : buffer(buffer), index(idx) {}

    /**
     * @brief Dereferences the iterator.
     * @return Reference to the element at the current position.
     */
    reference operator*() const { return (*buffer)[index]; }

    /**
     * @brief Member access through the iterator.
     * @return Pointer to the element at the current position.
     */
    pointer operator->() const { return &(*buffer)[index]; }

    /**
     * @brief Accesses the element offset n positions from the current one.
     * @param n Offset from the current position.
     * @return Reference to the element at index + n.
     */
    reference operator[](difference_type n) const {
      return (*buffer)[index + n];
    }

    /**
     * @brief Pre-increment; advances toward the back.
     * @return Reference to this iterator after advancing.
     */
    Derived &operator++() {
      ++index;
      return self();
    }

    /**
     * @brief Post-increment; advances toward the back.
     * @param int Unused disambiguation tag for post-increment.
     * @return Copy of the iterator before advancing.
     */
    Derived operator++(int) {
      Derived t = self();
      ++index;
      return t;
    }

    /**
     * @brief Pre-decrement; moves toward the front.
     * @return Reference to this iterator after moving.
     */
    Derived &operator--() {
      --index;
      return self();
    }

    /**
     * @brief Post-decrement; moves toward the front.
     * @param int Unused disambiguation tag for post-decrement.
     * @return Copy of the iterator before moving.
     */
    Derived operator--(int) {
      Derived t = self();
      --index;
      return t;
    }

    /**
     * @brief Advances the iterator by n positions.
     * @param n Number of positions to advance (may be negative).
     * @return Reference to this iterator after advancing.
     */
    Derived &operator+=(difference_type n) {
      index += n;
      return self();
    }

    /**
     * @brief Rewinds the iterator by n positions.
     * @param n Number of positions to rewind (may be negative).
     * @return Reference to this iterator after rewinding.
     */
    Derived &operator-=(difference_type n) {
      index -= n;
      return self();
    }

    /**
     * @brief Returns an iterator advanced n positions from a.
     * @param a Base iterator.
     * @param n Number of positions to advance.
     * @return Iterator at a.index + n.
     */
    friend Derived operator+(const Derived &a, difference_type n) {
      return Derived(a.buffer, a.index + n);
    }

    /**
     * @brief Returns an iterator advanced n positions from a.
     * @param n Number of positions to advance.
     * @param a Base iterator.
     * @return Iterator at a.index + n.
     */
    friend Derived operator+(difference_type n, const Derived &a) {
      return Derived(a.buffer, a.index + n);
    }

    /**
     * @brief Returns an iterator rewound n positions from a.
     * @param a Base iterator.
     * @param n Number of positions to rewind.
     * @return Iterator at a.index - n.
     */
    friend Derived operator-(const Derived &a, difference_type n) {
      return Derived(a.buffer, a.index - n);
    }

    /**
     * @brief Computes the distance between two iterators.
     * @param a Left iterator.
     * @param b Right iterator.
     * @return Signed number of positions from b to a.
     */
    friend difference_type operator-(const Derived &a, const Derived &b) {
      return (difference_type)a.index - (difference_type)b.index;
    }

    /**
     * @brief Tests two iterators for equality.
     * @param a Left iterator.
     * @param b Right iterator.
     * @return True if both reference the same buffer and position.
     */
    friend bool operator==(const Derived &a, const Derived &b) {
      return a.buffer == b.buffer && a.index == b.index;
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
      return a.index < b.index;
    }

    /**
     * @brief Orders two iterators by position.
     * @param a Left iterator.
     * @param b Right iterator.
     * @return True if a follows b.
     */
    friend bool operator>(const Derived &a, const Derived &b) {
      return a.index > b.index;
    }

    /**
     * @brief Orders two iterators by position.
     * @param a Left iterator.
     * @param b Right iterator.
     * @return True if a does not follow b.
     */
    friend bool operator<=(const Derived &a, const Derived &b) {
      return a.index <= b.index;
    }

    /**
     * @brief Orders two iterators by position.
     * @param a Left iterator.
     * @param b Right iterator.
     * @return True if a does not precede b.
     */
    friend bool operator>=(const Derived &a, const Derived &b) {
      return a.index >= b.index;
    }

  protected:
    BufPtr buffer; /**< Buffer this iterator traverses. */
    size_t index;  /**< Logical position, front-to-back. */

    friend class ConstIterator;

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
    ConstIterator(const iterator &other) : Base(other.buffer, other.index) {}
  };
};
