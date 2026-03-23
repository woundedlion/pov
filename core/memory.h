/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#ifndef HOLOSPHERE_CORE_MEMORY_H_
#define HOLOSPHERE_CORE_MEMORY_H_

#include <cstdint>
#include <cstddef>
#include <new>
#include <cassert>
#include <utility>
#ifndef ARDUINO
#include <cstdio>
#endif

constexpr size_t GLOBAL_ARENA_SIZE = 335 * 1024;

constexpr size_t DEFAULT_SCRATCH_A_SIZE = 16 * 1024;
constexpr size_t DEFAULT_SCRATCH_B_SIZE = 16 * 1024;
constexpr size_t DEFAULT_PERSISTENT_SIZE =
    GLOBAL_ARENA_SIZE - DEFAULT_SCRATCH_A_SIZE - DEFAULT_SCRATCH_B_SIZE;

// ============================================================================
// 1. Core Arena Allocator
// ============================================================================
class Arena {
  uint8_t *buffer;
  size_t capacity;
  size_t offset;
  size_t high_water_mark;
#ifndef NDEBUG
  uint32_t generation_ = 0;
#endif

public:
  Arena(uint8_t *buf, size_t size)
      : buffer(buf), capacity(size), offset(0), high_water_mark(0) {}

  void *allocate(size_t size, size_t align = alignof(std::max_align_t)) {
    size_t current = reinterpret_cast<size_t>(buffer + offset);
    size_t padding = (align - (current % align)) % align;
    if (offset + padding + size > capacity) {
#ifdef ARDUINO
      // No printf on Teensy — avoid pulling in 4KB of stdio
#else
      printf("[OOM] Arena: requested %zu bytes, offset %zu / capacity %zu\n",
             size, offset + padding, capacity);
#endif
      return nullptr;
    }
    offset += padding;
    void *ptr = buffer + offset;
    offset += size;
    if (offset > high_water_mark)
      high_water_mark = offset;
    assert(ptr != nullptr && "Arena::allocate returned null — OOM");
    return ptr;
  }

  size_t get_offset() const { return offset; }
  size_t get_capacity() const { return capacity; }
  size_t get_high_water_mark() const { return high_water_mark; }

  void set_offset(size_t new_offset) { offset = new_offset; }

  void reset() {
    offset = 0;
#ifndef NDEBUG
    generation_++;
#endif
  }

  void rebind(uint8_t *buf, size_t new_capacity) {
    buffer = buf;
    capacity = new_capacity;
    offset = 0;
    high_water_mark = 0;
#ifndef NDEBUG
    generation_++;
#endif
  }

  void reset_high_water_mark() { high_water_mark = offset; }

#ifndef NDEBUG
  uint32_t get_generation() const { return generation_; }
#endif
};

// ============================================================================
// 2. Triangular Bitset (Pair Deduplication)
// ============================================================================

/**
 * @brief Upper-triangular bitset for O(1) pair deduplication.
 * Stores one bit per unique unordered pair (a, b) where a < b < MAX_V.
 * Total storage: MAX_V * (MAX_V - 1) / 2 bits.
 * @tparam MAX_V Maximum vertex/element index (exclusive).
 */
template <int MAX_V> struct TriangularBitset {
  static constexpr int BITS = MAX_V * (MAX_V - 1) / 2;
  static constexpr int BYTES = (BITS + 7) / 8;
  uint8_t data[BYTES];

  void clear() { memset(data, 0, BYTES); }

  /// Bit index for pair (small, large) where small < large.
  static int index(int small, int large) {
    return small * (2 * MAX_V - small - 1) / 2 + (large - small - 1);
  }

  bool test(int a, int b) const {
    int bit = index(a, b);
    return (data[bit >> 3] >> (bit & 7)) & 1;
  }

  /// Returns true if already set (hit), false if newly inserted (miss).
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

template <typename T> class ArenaVector {
private:
  T *data_;
  size_t size_;
  size_t capacity_;
  bool bound_ = false;
#ifndef NDEBUG
  Arena *source_arena_ = nullptr;
  uint32_t birth_generation_ = 0;

  void check_alive() const {
    if (source_arena_ && source_arena_->get_generation() != birth_generation_) {
      assert(false && "ArenaVector use-after-free!");
    }
  }
#else
  void check_alive() const {}
#endif

  void check_bound() const {
    assert(bound_ && "Attempted to access unbound ArenaVector!");
  }

public:
  /// Default: creates an unbound vector. Must call bind() before use.
  ArenaVector() : data_(nullptr), size_(0), capacity_(0) {}

  // Disable implicit shallow copying to prevent memory aliasing bugs
  ArenaVector(const ArenaVector &) = delete;
  ArenaVector &operator=(const ArenaVector &) = delete;

  // Move semantics — moved-from object returns to pristine unbound state
  ArenaVector(ArenaVector &&other) noexcept
      : data_(other.data_), size_(other.size_), capacity_(other.capacity_),
        bound_(other.bound_) {
#ifndef NDEBUG
    source_arena_ = other.source_arena_;
    birth_generation_ = other.birth_generation_;
#endif
    other.data_ = nullptr;
    other.size_ = 0;
    other.capacity_ = 0;
    other.bound_ = false;
#ifndef NDEBUG
    other.source_arena_ = nullptr;
    other.birth_generation_ = 0;
#endif
  }

  ArenaVector &operator=(ArenaVector &&other) noexcept {
    if (this != &other) {
      data_ = other.data_;
      size_ = other.size_;
      capacity_ = other.capacity_;
      bound_ = other.bound_;
#ifndef NDEBUG
      source_arena_ = other.source_arena_;
      birth_generation_ = other.birth_generation_;
#endif
      other.data_ = nullptr;
      other.size_ = 0;
      other.capacity_ = 0;
      other.bound_ = false;
#ifndef NDEBUG
      other.source_arena_ = nullptr;
      other.birth_generation_ = 0;
#endif
    }
    return *this;
  }

  // Constructor requires exact capacity upfront. No dynamic growth.
  ArenaVector(Arena &arena, size_t exact_capacity)
      : data_(nullptr), size_(0), capacity_(0) {
    bind(arena, exact_capacity);
  }

  /// One-time binding to an arena. Asserts on double-bind.
  /// If already bound with sufficient capacity, resets size for reuse.
  void bind(Arena &arena, size_t exact_capacity) {
    if (bound_ && capacity_ >= exact_capacity) {
      size_ = 0; // Reuse memory
      return;
    }
    assert(!bound_ && "ArenaVector already bound!");
    if (exact_capacity > 0) {
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
#endif
  }

  bool is_bound() const { return bound_; }

  void push_back(const T &value) {
    check_alive();
    check_bound();
    assert(size_ < capacity_ && "ArenaVector exact capacity exceeded!");
    new (&data_[size_]) T(value);
    size_++;
  }

  /// Bulk-append from a contiguous source. T must be trivially copyable.
  void append_bulk(const T *src, size_t count) {
    static_assert(std::is_trivially_destructible_v<T>,
                  "append_bulk requires trivially destructible T (memcpy safety)");
    check_alive();
    check_bound();
    assert(size_ + count <= capacity_ && "ArenaVector bulk append exceeds capacity!");
    memcpy(static_cast<void*>(data_ + size_), src, count * sizeof(T));
    size_ += count;
  }

  template <typename... Args> T &emplace_back(Args &&...args) {
    check_alive();
    check_bound();
    assert(size_ < capacity_ && "ArenaVector exact capacity exceeded!");
    T *ptr = new (&data_[size_]) T(std::forward<Args>(args)...);
    size_++;
    return *ptr;
  }

  T &operator[](size_t i) {
    check_alive();
    check_bound();
    assert(i < size_);
    return data_[i];
  }
  const T &operator[](size_t i) const {
    check_alive();
    check_bound();
    assert(i < size_);
    return data_[i];
  }

  size_t size() const { return size_; }
  size_t capacity() const { return capacity_; }
  bool empty() const { return size_ == 0; }
  bool is_empty() const { return size_ == 0; }

  T &back() {
    check_alive();
    check_bound();
    assert(size_ > 0);
    return data_[size_ - 1];
  }
  const T &back() const {
    check_alive();
    check_bound();
    assert(size_ > 0);
    return data_[size_ - 1];
  }

  void clear() {
    size_ = 0; // Does not free memory or reallocate
  }

  T *data() { return data_; }
  const T *data() const { return data_; }

  T *begin() {
    check_alive();
    return data_;
  }
  T *end() {
    check_alive();
    return data_ + size_;
  }
  const T *begin() const {
    check_alive();
    return data_;
  }
  const T *end() const {
    check_alive();
    return data_ + size_;
  }
};

// ============================================================================
// 3. Non-Owning Span (Explicit Borrow)
// ============================================================================

/// A read-only, non-owning view into arena-allocated data.
/// Makes the distinction between owned (ArenaVector) and borrowed data
/// visible at the type level. Used to replace shallow_copy in transforms.
template <typename T> class ArenaSpan {
  const T *data_;
  size_t size_;

public:
  ArenaSpan() : data_(nullptr), size_(0) {}

  /// Can only be created FROM an ArenaVector (explicit borrow)
  explicit ArenaSpan(const ArenaVector<T> &source)
      : data_(source.data()), size_(source.size()) {}

  const T &operator[](size_t i) const {
    assert(i < size_);
    return data_[i];
  }
  size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }
  const T *data() const { return data_; }
  const T *begin() const { return data_; }
  const T *end() const { return data_ + size_; }
};

extern Arena scratch_arena_a;
extern Arena scratch_arena_b;

extern Arena persistent_arena;

void configure_arenas(size_t persistent, size_t scratch_a, size_t scratch_b);
void configure_arenas_default();

struct CompiledHankin;
namespace MeshOps {
template <typename MeshT>
void clone(const MeshT &src, MeshT &dst, Arena &arena);
void clone(const CompiledHankin &src, CompiledHankin &dst, Arena &arena);
} // namespace MeshOps

// ============================================================================
// 4. ScratchScope — RAII Guard + Factory for Temporary Memory
// ============================================================================

/// RAII guard that saves/restores an arena offset and provides a typed
/// factory for temporary ArenaVectors scoped to this block.
struct ScratchScope {
  Arena &arena;
  size_t saved_offset;

  ScratchScope(Arena &a) : arena(a), saved_offset(a.get_offset()) {}
  ~ScratchScope() { arena.set_offset(saved_offset); }

  /// Create a temporary ArenaVector scoped to this ScratchScope.
  template <typename T> ArenaVector<T> make_vector(size_t capacity) {
    return ArenaVector<T>(arena, capacity);
  }

  Arena &get_arena() { return arena; }

  // Non-copyable
  ScratchScope(const ScratchScope &) = delete;
  ScratchScope &operator=(const ScratchScope &) = delete;
};

// ============================================================================
// 5. RAII Arena Evacuator
// ============================================================================

/// Safely evacuates an object from the persistent arena to a scratch arena,
/// and automatically restores it upon destruction.
///
/// Usage:
///   {
///     Persist<MeshState> p(live_mesh, scratch_arena_a, persistent_arena);
///     persistent_arena.reset();
///   }  // ~Persist clones backup back into persistent
template <typename T> class Persist {
  T &target_;
  Arena &persistent_;

  // Declaration order matters! scratch_ must be declared BEFORE backup_
  // so that backup_ is destroyed before the scratch arena is rolled back.
  ScratchScope scratch_;
  T backup_;

public:
  Persist(T &target, Arena &scratch, Arena &persistent)
      : target_(target), persistent_(persistent), scratch_(scratch) {
    MeshOps::clone(target_, backup_, scratch_.get_arena());
  }

  ~Persist() {
    target_ = T();
    MeshOps::clone(backup_, target_, persistent_);
  }

  // Non-copyable
  Persist(const Persist &) = delete;
  Persist &operator=(const Persist &) = delete;
};
#endif // HOLOSPHERE_CORE_MEMORY_H_
