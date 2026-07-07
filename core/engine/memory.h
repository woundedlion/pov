/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <cstdint>
#include <cstddef>
#include <climits>
#include <cstring>
#include <new>
#include <cassert>
#include <utility>
#include <concepts>
#include "engine/platform.h"

// Device/simulator arena budget is 330 KiB. The native unit-test build
// (HS_TEST_BUILD) widens it so the effect smoke harness can exercise every
// effect's full render path regardless of per-effect arena tuning; the device
// footprint is unchanged. NOTE: the native harness is a 64-bit build, so
// per-effect arena footprints measured there can be LARGER than on the 32-bit
// device wherever a pooled struct embeds a POINTER (8 B host / 4 B device:
// ArenaVector's data ptr, Fn's callable ptr, BakedPalette::lut_). Index/count
// widths match across builds — StaticCircularBuffer uses uint32_t, so the big
// pooled structs (Particle/VectorTrail/OrientationTrail) are byte-identical on
// both builds. Pointer-bearing container HEADERS are not multiplied across pool
// elements, so the residual host/device delta is small but nonzero; do not
// treat the host high-water mark as an exact device figure. Effects tune their
// own split via configure_arenas() to fit the device budget.
// The real device FlexRAM (RAM1) arena. Sized from the measured worst-effect
// high-water (tests/arena_measure.cpp): MindSplatter is the binding tenant at
// ~327 KiB persistent + a 6 KiB scratch carve. 330 KiB clears that with the
// MindSplatter pool's compile-time margin intact while returning ~5 KiB of RAM1
// to the stack. A distinct always-defined constant (not the host-inflated
// GLOBAL_ARENA_SIZE below) so device-budget static_asserts — MindSplatter's
// particle pool — check the real figure even in the host suite.
constexpr size_t DEVICE_GLOBAL_ARENA_SIZE = 330 * 1024;
#ifdef HS_TEST_BUILD
constexpr size_t GLOBAL_ARENA_SIZE = 8 * 1024 * 1024;
#else
constexpr size_t GLOBAL_ARENA_SIZE = DEVICE_GLOBAL_ARENA_SIZE;
#endif

constexpr size_t DEFAULT_SCRATCH_A_SIZE = 16 * 1024;
constexpr size_t DEFAULT_SCRATCH_B_SIZE = 16 * 1024;
constexpr size_t DEFAULT_PERSISTENT_SIZE =
    GLOBAL_ARENA_SIZE - DEFAULT_SCRATCH_A_SIZE - DEFAULT_SCRATCH_B_SIZE;

// ============================================================================
// 1. Core Arena Allocator
// ============================================================================

/**
 * @brief Bump allocator over a fixed caller-owned buffer.
 * @details Allocation is offset advancement; individual frees are unsupported —
 * memory is reclaimed wholesale via reset() (rewind to 0) or set_offset()
 * (rewind to a saved mark). Over-allocation traps rather than returning null
 * (see allocate()).
 */
class Arena {
  uint8_t *buffer;
  size_t capacity;
  size_t offset;
  size_t high_water_mark;
#ifndef NDEBUG
  uint32_t generation_ = 0;
#endif

public:
  /**
   * @brief Constructs an arena over a caller-owned buffer.
   * @param buf Pointer to the backing buffer.
   * @param size Capacity of the buffer in bytes.
   */
  Arena(uint8_t *buf, size_t size)
      : buffer(buf), capacity(size), offset(0), high_water_mark(0) {}

  /**
   * @brief Bump-allocate `size` bytes aligned to `align`, advancing the offset.
   * @param size Number of bytes to allocate.
   * @param align Required alignment in bytes; defaults to max_align_t.
   * @return Pointer into the buffer for the allocated block.
   * @details Traps (HS_CHECK) on over-allocation rather than returning null.
   * Updates the high-water mark. `size` must be > 0: a zero-size request would
   * return a valid bump pointer that reserves no storage (it aliases the next
   * allocation's address), so a caller treating it as ownable would silently
   * overlap. Every real caller passes a positive byte count (ArenaVector::bind
   * guards `exact_capacity > 0` before calling), so trap a zero request as the
   * misuse it is rather than hand back a non-owning pointer.
   */
  void *allocate(size_t size, size_t align = alignof(std::max_align_t)) {
    HS_CHECK(size > 0, "Arena::allocate: zero-size request");
    HS_CHECK(align != 0 && (align & (align - 1)) == 0);
    uintptr_t current = reinterpret_cast<uintptr_t>(buffer + offset);
    size_t padding = (align - (current % align)) % align;
    // Subtractive form: offset <= capacity is invariant, so it cannot wrap the
    // way `offset + padding + size > capacity` would for a colossal `size`.
    if (padding > capacity - offset || size > capacity - offset - padding) {
      hs::log("[OOM] Arena: req %zu, offset %zu, pad %zu / cap %zu", size,
              offset, padding, capacity);
      HS_CHECK(false);
    }
    offset += padding;
    void *ptr = buffer + offset;
    offset += size;
    if (offset > high_water_mark)
      high_water_mark = offset;
    return ptr;
  }

  /**
   * @brief Bump-allocate storage for `n` elements of `T`, typed.
   * @tparam T Element type; sizes and aligns the block from the type.
   * @param n Element count (must be > 0, per allocate()).
   * @return Pointer to the block, cast to `T*`.
   * @details Thin wrapper over allocate() that derives `sizeof`/`alignof` from
   * `T` so a call site cannot mis-pair them. Does not construct the elements.
   */
  template <typename T> T *allocate_n(size_t n) {
    return static_cast<T *>(allocate(n * sizeof(T), alignof(T)));
  }

  /**
   * @brief Returns the current allocation offset.
   * @return Bytes consumed from the buffer so far.
   */
  size_t get_offset() const { return offset; }
  /**
   * @brief Returns the arena's total capacity.
   * @return Capacity of the backing buffer in bytes.
   */
  size_t get_capacity() const { return capacity; }
  /**
   * @brief Returns the peak allocation offset observed.
   * @return High-water mark in bytes.
   */
  size_t get_high_water_mark() const { return high_water_mark; }

  /**
   * @brief Rewinds the offset to a previously saved mark.
   * @param new_offset Offset to rewind to; must be <= the current offset.
   * @details A mark is only valid as a rewind target: jumping the offset
   * *forward* would hand out backing bytes that were never reserved by an
   * allocate() call. Trap any non-rewind at the source rather than letting it
   * silently resurrect or over-run storage. (new_offset <= offset also implies
   * new_offset <= capacity, preserving the no-wrap bounds math in allocate().)
   * Alignment is not re-checked here: allocate() recomputes leading padding from
   * the true address on every call, so restoring an unaligned mark (e.g. a
   * ScratchScope save) is safe — the next allocation realigns itself.
   */
  void set_offset(size_t new_offset) {
    HS_CHECK(new_offset <= offset);
    offset = new_offset;
  }

  /**
   * @brief Rewind to empty, reclaiming all allocations at once.
   * @details Bumps the debug generation so any live ArenaVector/ArenaSpan into
   * the old contents faults.
   */
  void reset() {
    offset = 0;
#ifndef NDEBUG
    generation_++;
#endif
  }

  /**
   * @brief Point the arena at a different buffer/capacity and reset to empty.
   * @param buf Pointer to the new backing buffer.
   * @param new_capacity Capacity of the new buffer in bytes.
   * @details Used by configure_arenas to repartition the global budget at
   * runtime.
   */
  void rebind(uint8_t *buf, size_t new_capacity) {
    buffer = buf;
    capacity = new_capacity;
    offset = 0;
    high_water_mark = 0;
#ifndef NDEBUG
    generation_++;
#endif
  }

  /**
   * @brief Reset peak-usage tracking to the current offset.
   * @details E.g. to measure a single frame's allocation peak in isolation.
   */
  void reset_high_water_mark() { high_water_mark = offset; }

#ifndef NDEBUG
  /**
   * @brief Returns the current debug generation stamp.
   * @return Generation counter, bumped on each reset/rebind.
   */
  uint32_t get_generation() const { return generation_; }
#endif
};

// ============================================================================
// 2. Triangular Bitset (Pair Deduplication)
// ============================================================================

/**
 * @brief Upper-triangular bitset for O(1) pair deduplication.
 * @tparam MAX_V Maximum vertex/element index (exclusive).
 * @details Stores one bit per unique unordered pair (a, b) where a < b < MAX_V.
 * Total storage: MAX_V * (MAX_V - 1) / 2 bits.
 */
template <int MAX_V> struct TriangularBitset {
  // BITS (here) and index() below both form an intermediate product on the
  // order of MAX_V^2 in `int` before halving. For MAX_V >= ~46341 that product
  // overflows int32 and silently corrupts the bit layout. No current
  // instantiation is anywhere near this; the static_assert pins the ceiling so
  // a future large-mesh instantiation fails at compile time, not at runtime.
  static_assert(static_cast<long long>(MAX_V) * MAX_V <= INT_MAX,
                "TriangularBitset: MAX_V too large; index() product overflows int");
  static constexpr int BITS = MAX_V * (MAX_V - 1) / 2;
  static constexpr int BYTES = (BITS + 7) / 8;
  uint8_t data[BYTES] = {}; /**< Packed bit storage; zero-initialized so a pair
                                 read before clear() reads "unset" rather than UB. */

  /**
   * @brief Clears every pair bit to zero.
   */
  void clear() { memset(data, 0, BYTES); }

  /**
   * @brief Bit index for pair (small, large) where small < large.
   * @param small Lower index of the pair, in [0, large).
   * @param large Higher index of the pair, in (small, MAX_V).
   * @return Bit index into the packed storage.
   */
  static int index(int small, int large) {
    // The triangular layout is only valid for an ordered, in-range pair: a
    // swapped pair silently aliases the wrong bit (dedup corruption) and an
    // out-of-range one writes adjacent memory. HS_CHECK (survives NDEBUG, traps
    // on device) matches every sibling write primitive in memory.h, so an
    // unordered or out-of-range caller fails fast instead of corrupting memory.
    // index() feeds test()/test_and_set() on the per-edge mesh-dedup setup path
    // (plot.h draw()), not a per-pixel loop, and that caller already runs
    // always-on HS_CHECKs alongside it, so the guard is contract-appropriate.
    HS_CHECK(small >= 0 && small < large && large < MAX_V);
    return small * (2 * MAX_V - small - 1) / 2 + (large - small - 1);
  }

  /**
   * @brief Tests whether pair (a, b) is set.
   * @param a Lower index of the pair; requires a < b.
   * @param b Higher index of the pair (see index()).
   * @return True iff the pair's bit is set.
   */
  bool test(int a, int b) const {
    int bit = index(a, b);
    return (data[bit >> 3] >> (bit & 7)) & 1;
  }

  /**
   * @brief Tests and sets the bit for pair (a, b).
   * @param a Lower index of the pair; requires a < b.
   * @param b Higher index of the pair (see index()).
   * @return True if already set (hit), false if newly inserted (miss).
   */
  bool test_and_set(int a, int b) {
    int bit = index(a, b);
    uint8_t &byte = data[bit >> 3];
    uint8_t mask = 1 << (bit & 7);
    if (byte & mask)
      return true;
    byte |= mask;
    return false;
  }
};

// ============================================================================
// 3. Arena Structures
// ============================================================================

/**
 * @brief Fixed-capacity, arena-backed vector. Move-only; no dynamic growth.
 * @tparam T Element type; must satisfy the element destructor contract below.
 * @details ELEMENT DESTRUCTOR CONTRACT: ArenaVector deliberately does NOT run
 * element destructors. clear(), move, move-assign, and going out of scope all
 * leave stored elements un-destructed — storage is owned and reclaimed by the
 * arena (reset/compaction), not by this handle. This is a hard requirement of
 * the arena model, not an oversight: an arena can be reset out from under a
 * still-live ArenaVector (see bind()'s stale-binding handling), so a
 * handle-driven destructor pass could run on already-reclaimed memory.
 * Consequently, only store types whose destructor need not run for correctness
 * — trivially-destructible PODs, or Fn<> (inline-storage on both backends) whose
 * stored captures are themselves trivial, so the callable owns no external
 * resource regardless of size. A type that owns heap/handles outside the arena
 * must not be stored in an ArenaVector — notably a raw std::function, which
 * heap-allocates any capture past its small-buffer size and would leak when this
 * handle skips its destructor.
 */
template <typename T> class ArenaVector {
  // ArenaSpan borrows our backing data and (in debug builds) our arena
  // generation stamp for its own use-after-free check.
  template <typename U> friend class ArenaSpan;

private:
  T *data_;            /**< Pointer to the arena-allocated backing block. */
  size_t size_;        /**< Number of constructed elements. */
  size_t capacity_;    /**< Maximum element count the block can hold. */
  bool bound_ = false; /**< Whether the vector has been bound to an arena. */
#ifndef NDEBUG
  Arena *source_arena_ = nullptr;   /**< Arena the block was allocated from. */
  uint32_t birth_generation_ = 0;   /**< Arena generation at bind time. */
  /**
   * @brief Per-vector counter bumped on every fresh allocation in bind() (the
   * grow / re-bind path).
   * @details A grow swaps data_ for a new block WITHOUT touching the arena
   * generation, so this is the only signal an ArenaSpan can use to detect that
   * its snapshotted pointer was abandoned by a re-grow of the source vector.
   */
  uint32_t rebind_generation_ = 0;

  /**
   * @brief Debug-only use-after-free check against the arena generation.
   * @details Asserts if the source arena was reset out from under this vector.
   */
  void check_alive() const {
    if (source_arena_ && source_arena_->get_generation() != birth_generation_) {
      assert(false && "ArenaVector use-after-free!");
    }
  }
#else
  /**
   * @brief No-op use-after-free check in release builds.
   */
  void check_alive() const {}
#endif

  /**
   * @brief Asserts that the vector has been bound to an arena.
   */
  void check_bound() const {
    assert(bound_ && "Attempted to access unbound ArenaVector!");
  }

  /**
   * @brief Transfers @p other's storage and bookkeeping into this vector,
   *        leaving @p other in a pristine unbound state.
   */
  void steal_from(ArenaVector &other) noexcept {
    data_ = other.data_;
    size_ = other.size_;
    capacity_ = other.capacity_;
    bound_ = other.bound_;
#ifndef NDEBUG
    source_arena_ = other.source_arena_;
    birth_generation_ = other.birth_generation_;
    rebind_generation_ = other.rebind_generation_;
#endif
    other.data_ = nullptr;
    other.size_ = 0;
    other.capacity_ = 0;
    other.bound_ = false;
#ifndef NDEBUG
    other.source_arena_ = nullptr;
    other.birth_generation_ = 0;
    other.rebind_generation_ = 0;
#endif
  }

public:
  /**
   * @brief Default-constructs an unbound vector.
   * @details Must call bind() before use.
   */
  ArenaVector() : data_(nullptr), size_(0), capacity_(0) {}

  /**
   * @brief Deleted copy constructor.
   * @details Implicit shallow copying is disabled to prevent memory aliasing.
   */
  ArenaVector(const ArenaVector &) = delete;
  /**
   * @brief Deleted copy assignment.
   * @return Reference to this (never invoked).
   * @details Implicit shallow copying is disabled to prevent memory aliasing.
   */
  ArenaVector &operator=(const ArenaVector &) = delete;

  /**
   * @brief Move constructor.
   * @param other Source vector; left in a pristine unbound state.
   */
  ArenaVector(ArenaVector &&other) noexcept { steal_from(other); }

  /**
   * @brief Move assignment.
   * @param other Source vector; left in a pristine unbound state.
   * @return Reference to this.
   */
  ArenaVector &operator=(ArenaVector &&other) noexcept {
    if (this != &other)
      steal_from(other);
    return *this;
  }

  /**
   * @brief Constructs and binds the vector with an exact capacity.
   * @param arena Arena to allocate the backing block from.
   * @param exact_capacity Element count to allocate; no dynamic growth.
   */
  ArenaVector(Arena &arena, size_t exact_capacity)
      : data_(nullptr), size_(0), capacity_(0) {
    bind(arena, exact_capacity);
  }

  /**
   * @brief Binds the vector to an arena, allocating its backing block.
   * @param arena Arena to allocate from.
   * @param exact_capacity Element count to reserve.
   * @details If already bound with sufficient capacity, resets size for reuse.
   * A grow against the same arena/generation reallocates a fresh block and
   * abandons the old one until the next reset; a stale binding (arena reset or
   * a different arena) trips a debug-only contract assert.
   */
  void bind(Arena &arena, size_t exact_capacity) {
#ifndef NDEBUG
    // Rebinding a still-bound vector after its source arena was reset, or to a
    // different arena, is a contract violation (the old block is already dead). A
    // same-arena/same-generation grow is not stale and reallocates below.
    assert((!bound_ || (source_arena_ == &arena &&
                        birth_generation_ == arena.get_generation())) &&
           "ArenaVector::bind() on a stale binding: clear the handle before "
           "resetting or changing its arena");
#endif
    // Same arena, still live, and big enough → reuse the block in place.
    if (bound_ && capacity_ >= exact_capacity) {
      size_ = 0;
#ifndef NDEBUG
      // Reuse dangles any span snapshotted before this point; bump so its
      // check_alive() trips (the arena generation alone won't).
      rebind_generation_++;
#endif
      return;
    }
    // Otherwise (unbound, or a grow that abandons the old block) → allocate
    // fresh. A grow leaks the old block until the next arena reset/compaction.
#ifndef NDEBUG
    if (bound_)
      hs::log("ArenaVector grow abandons %zu bytes (cap %zu -> %zu)",
              capacity_ * sizeof(T), capacity_, exact_capacity);
#endif
    if (exact_capacity > 0) {
      // Trap an exact_capacity * sizeof(T) overflow that would wrap to a small
      // byte count and slip past Arena::allocate's bounds check.
      HS_CHECK(exact_capacity <= SIZE_MAX / sizeof(T),
               "ArenaVector capacity * sizeof(T) overflows size_t!");
      data_ = static_cast<T *>(
          arena.allocate(exact_capacity * sizeof(T), alignof(T)));
    } else {
      data_ = nullptr;
    }
    size_ = 0;
    capacity_ = exact_capacity;
    bound_ = true;
#ifndef NDEBUG
    source_arena_ = &arena;
    birth_generation_ = arena.get_generation();
    rebind_generation_++;
#endif
  }

  /**
   * @brief Reports whether the vector is bound to an arena.
   * @return True iff bind() has been called.
   */
  bool is_bound() const { return bound_; }

  /**
   * @brief Appends a copy of an element.
   * @param value Element to copy-construct at the end.
   */
  void push_back(const T &value) {
    check_alive();
    check_bound();
    HS_CHECK(size_ < capacity_, "ArenaVector exact capacity exceeded!");
    new (&data_[size_]) T(value);
    size_++;
  }

  /**
   * @brief Bulk-append from a contiguous source.
   * @param src Pointer to the first source element.
   * @param count Number of elements to copy.
   * @details T must be trivially copyable (memcpy'd). An empty append is a
   * no-op that skips memcpy to avoid null-pointer UB.
   */
  void append_bulk(const T *src, size_t count) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "append_bulk memcpy's the source; T must be trivially copyable");
    check_alive();
    check_bound();
    // Subtractive, wrap-proof form: `size_ + count` could wrap for a colossal count.
    HS_CHECK(count <= capacity_ - size_,
             "ArenaVector bulk append exceeds capacity!");
    // Skip memcpy on an empty append: a null src with count 0 is formal UB.
    if (count == 0)
      return;
    memcpy(static_cast<void*>(data_ + size_), src, count * sizeof(T));
    size_ += count;
  }

  /**
   * @brief Constructs an element in place at the end.
   * @tparam Args Constructor argument types forwarded to T.
   * @param args Arguments forwarded to T's constructor.
   * @return Reference to the newly constructed element.
   */
  template <typename... Args> T &emplace_back(Args &&...args) {
    check_alive();
    check_bound();
    HS_CHECK(size_ < capacity_, "ArenaVector exact capacity exceeded!");
    T *ptr = new (&data_[size_]) T(std::forward<Args>(args)...);
    size_++;
    return *ptr;
  }

  /**
   * @brief Element access by index.
   * @param i Index in [0, size()).
   * @return Mutable reference to the element at index i.
   */
  T &operator[](size_t i) {
    check_alive();
    check_bound();
    assert(i < size_);
    return data_[i];
  }
  /**
   * @brief Element access by index (const).
   * @param i Index in [0, size()).
   * @return Const reference to the element at index i.
   */
  const T &operator[](size_t i) const {
    check_alive();
    check_bound();
    assert(i < size_);
    return data_[i];
  }

  /**
   * @brief Returns the number of stored elements.
   * @return Current element count.
   */
  size_t size() const { return size_; }
  /**
   * @brief Returns the maximum element count.
   * @return Capacity in elements.
   */
  size_t capacity() const { return capacity_; }
  /**
   * @brief Reports whether the vector is empty.
   * @return True iff size() == 0.
   */
  bool is_empty() const { return size_ == 0; }

  /**
   * @brief Accesses the last element.
   * @return Mutable reference to the element at index size() - 1.
   */
  T &back() {
    check_alive();
    check_bound();
    assert(size_ > 0);
    return data_[size_ - 1];
  }
  /**
   * @brief Accesses the last element (const).
   * @return Const reference to the element at index size() - 1.
   */
  const T &back() const {
    check_alive();
    check_bound();
    assert(size_ > 0);
    return data_[size_ - 1];
  }

  /**
   * @brief Resets the vector to empty without destroying elements.
   * @details No check_bound() here on purpose: clear() only resets size_ (it
   * neither frees nor touches data_), so it is a defined no-op on an unbound
   * vector. MeshState::clear() relies on this to reset members that may never
   * have been bound.
   */
  void clear() {
    check_alive();
    size_ = 0;
  }

  /**
   * @brief Returns a pointer to the backing storage.
   * @return Mutable pointer to the first element, or nullptr if unbound/empty.
   * @details No check_bound() on data()/begin()/end() on purpose: an unbound
   * (or moved-from) vector is data_==nullptr with size_==0, i.e. a well-defined
   * EMPTY range. Callers rely on that idiom — mesh.h/spatial.h's *_data()
   * accessors, std::span(vertices.data(), size()), and std::sort(data(),
   * data()+count) all pass these on a size-0 vector. A bound-precondition here
   * would false-fire on correct code. The use-after-free guard (check_alive)
   * still applies.
   */
  T *data() {
    check_alive();
    return data_;
  }
  /**
   * @brief Returns a pointer to the backing storage (const).
   * @return Const pointer to the first element, or nullptr if unbound/empty.
   */
  const T *data() const {
    check_alive();
    return data_;
  }

  /**
   * @brief Returns an iterator to the first element.
   * @return Mutable pointer to the first element.
   */
  T *begin() {
    check_alive();
    return data_;
  }
  /**
   * @brief Returns an iterator past the last element.
   * @return Mutable pointer one past the last element.
   */
  T *end() {
    check_alive();
    // Guard nullptr + 0 (formal UB) on an unbound/empty vector.
    return data_ ? data_ + size_ : nullptr;
  }
  /**
   * @brief Returns a const iterator to the first element.
   * @return Const pointer to the first element.
   */
  const T *begin() const {
    check_alive();
    return data_;
  }
  /**
   * @brief Returns a const iterator past the last element.
   * @return Const pointer one past the last element.
   */
  const T *end() const {
    check_alive();
    // Guard nullptr + 0 (formal UB) on an unbound/empty vector.
    return data_ ? data_ + size_ : nullptr;
  }
};

// ============================================================================
// 4. Non-Owning Span (Explicit Borrow)
// ============================================================================

/**
 * @brief A read-only, non-owning view into arena-allocated data.
 * @tparam T Element type viewed by the span.
 * @details Makes the distinction between owned (ArenaVector) and borrowed data
 * visible at the type level.
 *
 * LIFETIME CONTRACT: a span snapshots its source vector's data_ pointer at
 * construction. In debug builds two independent stamps fault on a stale span:
 * the arena generation catches an arena RESET out from under the span, and the
 * source vector's per-vector rebind counter catches a bind()-driven RE-GROW (a
 * grow allocates a fresh block and rebinds data_ WITHOUT bumping the arena
 * generation, so the rebind counter is the only signal). Re-taking the span
 * after a grow is still the right pattern; the counter just converts a silent
 * dangle into a debug fault.
 *
 * A MOVE of the source vector is not tracked: the span keeps its snapshotted
 * data_ (runtime-safe), but its debug stamps reference the moved-from husk, so
 * re-take the span after moving its source.
 */
template <typename T> class ArenaSpan {
  const T *data_; /**< Snapshotted pointer to the borrowed data. */
  size_t size_;   /**< Number of viewed elements. */
#ifndef NDEBUG
  Arena *source_arena_ = nullptr;          /**< Source arena for the stamp. */
  uint32_t birth_generation_ = 0;          /**< Arena generation at construction. */
  const ArenaVector<T> *source_vec_ = nullptr; /**< Source vector for re-grow check. */
  uint32_t source_rebind_generation_ = 0;  /**< Vector rebind counter at construction. */

  /**
   * @brief Debug-only stale-span check against arena and vector stamps.
   * @details Asserts on an arena reset or a source-vector re-grow.
   */
  void check_alive() const {
    if (source_arena_ && source_arena_->get_generation() != birth_generation_) {
      assert(false && "ArenaSpan use-after-free!");
    }
    if (source_vec_ &&
        source_vec_->rebind_generation_ != source_rebind_generation_) {
      assert(false && "ArenaSpan source vector re-grown out from under span!");
    }
  }
#else
  /**
   * @brief No-op stale-span check in release builds.
   */
  void check_alive() const {}
#endif

public:
  /**
   * @brief Default-constructs an empty span.
   */
  ArenaSpan() : data_(nullptr), size_(0) {}

  /**
   * @brief Constructs a span borrowing from an ArenaVector (explicit borrow).
   * @param source Vector to borrow data and (in debug) lifetime stamps from.
   * @details In debug builds the span inherits the vector's arena-generation
   * stamp, so accessing it after the arena is reset/compacted trips the same
   * use-after-free check the vector itself has.
   */
  explicit ArenaSpan(const ArenaVector<T> &source)
      : data_(source.data()), size_(source.size())
#ifndef NDEBUG
        ,
        source_arena_(source.source_arena_),
        birth_generation_(source.birth_generation_), source_vec_(&source),
        source_rebind_generation_(source.rebind_generation_)
#endif
  {
  }

  /**
   * @brief Deleted constructor from a temporary ArenaVector.
   * @details Borrowing from a temporary would leave the data pointer and (in
   * debug) source_vec_ dangling the moment the temporary dies. Forbid it.
   */
  explicit ArenaSpan(const ArenaVector<T> &&) = delete;

  /**
   * @brief Copy and copy-assignment duplicate the borrow verbatim.
   * @details Explicitly defaulted (not left implicit) to document intent: a span
   * copy carries the same data pointer, size, and — in debug — the source's
   * lifetime stamps, so the copy trips the same staleness check as the original.
   * "Re-take the same view" is the supported pattern; only construction from a
   * temporary ArenaVector (above) is forbidden, since copying an existing,
   * already-validated span borrows nothing new.
   */
  ArenaSpan(const ArenaSpan &) = default;
  ArenaSpan &operator=(const ArenaSpan &) = default;

  /**
   * @brief Element access by index.
   * @param i Index in [0, size()).
   * @return Const reference to the element at index i.
   */
  const T &operator[](size_t i) const {
    check_alive();
    assert(i < size_);
    return data_[i];
  }
  /**
   * @brief Returns the number of viewed elements.
   * @return Element count.
   */
  size_t size() const { return size_; }
  /**
   * @brief Reports whether the span is empty.
   * @return True iff size() == 0.
   */
  bool empty() const { return size_ == 0; }
  /**
   * @brief Returns a pointer to the borrowed storage.
   * @return Const pointer to the first element.
   */
  const T *data() const {
    check_alive();
    return data_;
  }
  /**
   * @brief Returns a const iterator to the first element.
   * @return Const pointer to the first element.
   */
  const T *begin() const {
    check_alive();
    return data_;
  }
  /**
   * @brief Returns a const iterator past the last element.
   * @return Const pointer one past the last element.
   */
  const T *end() const {
    check_alive();
    // Guard nullptr + 0 (formal UB) on a default-constructed/empty span.
    return data_ ? data_ + size_ : nullptr;
  }
};

extern Arena scratch_arena_a;
extern Arena scratch_arena_b;

extern Arena persistent_arena;

/**
 * @brief Repartitions the global arena budget across the three arenas.
 * @param persistent Bytes to assign to the persistent arena.
 * @param scratch_a Bytes to assign to scratch arena A.
 * @param scratch_b Bytes to assign to scratch arena B.
 */
void configure_arenas(size_t persistent, size_t scratch_a, size_t scratch_b);
/**
 * @brief Restores the default arena partition.
 */
void configure_arenas_default();

// ============================================================================
// 5. ScratchScope — RAII Arena Offset Guard
// ============================================================================

/**
 * @brief RAII guard that saves/restores an arena offset.
 * @note Only allocations made after construction are reclaimed: anything bound
 * to the arena before the scope opens sits below the saved offset and survives.
 * An operator that produces output in the same arena it scratches (e.g. the
 * Conway operators' output-mesh vectors over `target`) must therefore bind that
 * output before constructing the scope, or scope exit reclaims it.
 */
struct ScratchScope {
  Arena &arena;        /**< Arena whose offset is saved and restored. */
  size_t saved_offset; /**< Offset captured at construction. */

  /**
   * @brief Constructs the scope, saving the arena's current offset.
   * @param a Arena to guard.
   */
  ScratchScope(Arena &a) : arena(a), saved_offset(a.get_offset()) {}
  /**
   * @brief Destroys the scope, rewinding the arena to the saved offset.
   * @details Enforces LIFO scope discipline before rewinding. A bump allocator
   * makes stack-nested scopes on the SAME arena always safe: an inner scope
   * saves a mark above the outer's allocations and rewinds to exactly where they
   * end, so nesting never clobbers a live allocation. That safety holds for ALL
   * stack-nested uses (e.g. the multiple scratch_arena_a scopes Plot::rasterize
   * and Pixel::Feedback::flush hold) — the ONE way to break it is non-LIFO
   * teardown: a scope rewinding while a still-live inner allocation, or an outer
   * scope, sits above it. That shows up here as the arena offset having dropped
   * BELOW saved_offset (an outer rewind or reset() ran while this scope was
   * live). Trap it instead of letting set_offset() resurrect freed bytes. This
   * is the structural enforcement of the scratch-arena temporal-overlap contract
   * documented at Pixel::Feedback::flush (filter.h) and Plot::rasterize.
   */
  ~ScratchScope() {
    HS_CHECK(arena.get_offset() >= saved_offset);
    arena.set_offset(saved_offset);
  }

  /**
   * @brief Returns the guarded arena.
   * @return Reference to the underlying arena.
   */
  Arena &get_arena() { return arena; }

  /**
   * @brief Deleted copy constructor (non-copyable).
   */
  ScratchScope(const ScratchScope &) = delete;
  /**
   * @brief Deleted copy assignment (non-copyable).
   * @return Reference to this (never invoked).
   */
  ScratchScope &operator=(const ScratchScope &) = delete;
};

// ============================================================================
// 6. RAII Arena Evacuator
// ============================================================================

/**
 * @brief Concept requiring a static clone(const T&, T&, Arena&) method.
 * @tparam T Type that must provide static void clone(const T&, T&, Arena&).
 * @note Cloneable only constrains the clone hook. `Persist<T>` needs more (T
 *       default-initializable and assignable, because ~Persist does
 *       `target_ = T()` before restoring) and enforces those extra requirements
 *       with its own static_asserts rather than widening this concept, which is
 *       reused wherever a clone hook alone is the contract.
 */
template <typename T>
concept Cloneable = requires(const T &src, T &dst, Arena &arena) {
  { T::clone(src, dst, arena) } -> std::same_as<void>;
};

/**
 * @brief RAII evacuator that moves an object to scratch and restores it later.
 * @tparam T Cloneable target type.
 * @details Safely evacuates an object from the persistent arena to a scratch
 * arena, and automatically restores it upon destruction.
 *
 * Usage:
 *   {
 *     Persist<MeshState> p(live_mesh, scratch_arena_a, persistent_arena);
 *     persistent_arena.reset();
 *   }  // ~Persist clones backup back into persistent
 */
template <Cloneable T> class Persist {
  T &target_;         /**< Object being evacuated and restored. */
  Arena &persistent_; /**< Arena the object is restored into. */
  size_t persistent_offset_at_ctor_; /**< persistent_ offset at construction; the
                                          dtor traps unless the caller rewound
                                          below this watermark. */

  // Declaration order matters! scratch_ must be declared BEFORE backup_
  // so that backup_ is destroyed before the scratch arena is rolled back.
  ScratchScope scratch_; /**< Scratch scope holding the backup's storage. */
  T backup_;             /**< Cloned backup of the target in scratch memory. */

  // ~Persist does `target_ = T()` before re-cloning the backup, so T needs more
  // than Cloneable (which only requires a clone hook). Assert the extra
  // requirements here so a Cloneable-but-not-default-constructible/assignable T
  // fails with a clear message at instantiation instead of an opaque error
  // buried inside the destructor body.
  static_assert(std::default_initializable<T>,
                "Persist<T>: ~Persist resets target_ = T() before restoring, so "
                "T must be default-initializable (Cloneable does not imply this).");
  static_assert(std::assignable_from<T &, T>,
                "Persist<T>: ~Persist assigns target_ = T(), so T must be "
                "assignable from a T rvalue (Cloneable does not imply this).");

public:
  /**
   * @brief Evacuates the target into the scratch arena.
   * @param target Object to back up and later restore.
   * @param scratch Scratch arena to hold the backup.
   * @param persistent Persistent arena the target is restored into.
   */
  Persist(T &target, Arena &scratch, Arena &persistent)
      : target_(target), persistent_(persistent),
        persistent_offset_at_ctor_(persistent.get_offset()), scratch_(scratch) {
    HS_CHECK(&scratch != &persistent,
             "Persist: scratch and persistent must be distinct arenas — the "
             "dtor's watermark restore assumes the backup lives in a different "
             "arena than the one it restores into");
    T::clone(target_, backup_, scratch_.get_arena());
  }

  /**
   * @brief Restores the target by cloning the backup into the persistent arena.
   * @details The restore clones into `persistent_` at its *current* offset, so
   * it only reconstructs the object usefully if the caller rewound the
   * persistent arena during the scope (the canonical `persistent_arena.reset()`
   * in the usage example). Without that reset the clone appends a second copy
   * and grows the arena rather than restoring. The dtor traps *after* the
   * restore unless the persistent arena ends no larger than it was at
   * construction.
   *
   * The check is post-restore and uses `<=` (not a pre-restore `<`) on purpose:
   * (1) the Persist ctor allocates only in *scratch*, so the captured watermark
   * is the persistent high-water mark when the block was entered, and a correct
   * reset+restore — even one that compacts and relocates the object — must end
   * at or below it; (2) callers legitimately STACK several Persists over one
   * `reset()` (e.g. HankinSolids' mesh compaction), and their reverse-order
   * restores append onto the compacted base, so an individual restore's offset
   * is only meaningful against the entry watermark, post-clone. A forgot-to-
   * reset restore instead appends beyond the live original, pushing the offset
   * past the watermark — exactly what this traps.
   *
   * Scope of the guarantee: this is a backstop, not a proof. It reliably catches
   * the common single-Persist forgot-to-reset (the offset clears the watermark).
   * In a stack it can miss a per-Persist mistake whose growth happens to land at
   * or below the entry watermark — e.g. a later restore that compacts below it
   * masks an earlier one's overrun. It bounds the aggregate, not each restore.
   */
  ~Persist() {
    target_ = T();
    T::clone(backup_, target_, persistent_);
    HS_CHECK(persistent_.get_offset() <= persistent_offset_at_ctor_,
             "Persist: restore grew the persistent arena past its construction "
             "watermark — the caller did not rewind/reset it during the scope, "
             "so the restore appended a duplicate instead of reconstructing");
  }

  /**
   * @brief Deleted copy constructor (non-copyable).
   */
  Persist(const Persist &) = delete;
  /**
   * @brief Deleted copy assignment (non-copyable).
   * @return Reference to this (never invoked).
   */
  Persist &operator=(const Persist &) = delete;
};
