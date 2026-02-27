/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <cstdint>
#include <cstddef>
#include <new>
#include <cassert>
#include <utility>
#include <cstdio>

constexpr size_t GEOMETRY_ARENA_SIZE = 128 * 1024;
constexpr size_t SCRATCH_ARENA_A_SIZE = 128 * 1024;
constexpr size_t SCRATCH_ARENA_B_SIZE = 256 * 1024;

// ============================================================================
// 1. Core Arena Allocator
// ============================================================================
class Arena {
  uint8_t *buffer;
  size_t capacity;
  size_t offset;
  size_t high_water_mark;

public:
  Arena(size_t size) : capacity(size), offset(0), high_water_mark(0) {
    buffer = new uint8_t[size];
  }

  ~Arena() { delete[] buffer; }

  void *allocate(size_t size, size_t align = alignof(std::max_align_t)) {
    size_t current = reinterpret_cast<size_t>(buffer + offset);
    size_t padding = (align - (current % align)) % align;
    if (offset + padding + size > capacity) {
      printf("[OOM] Arena out of memory! Requested: %zu bytes. Offset: %zu / "
             "Capacity: %zu\n",
             size, offset + padding, capacity);
      assert(false && "Arena out of memory!");
      return nullptr;
    }
    offset += padding;
    void *ptr = buffer + offset;
    offset += size;
    if (offset > high_water_mark) {
      high_water_mark = offset;
    }
    return ptr;
  }

  size_t get_offset() const { return offset; }
  size_t get_capacity() const { return capacity; }
  size_t get_high_water_mark() const { return high_water_mark; }
  void set_offset(size_t new_offset) { offset = new_offset; }
  void reset() { offset = 0; }
  void reset_high_water_mark() { high_water_mark = offset; }
};

// RAII Marker for scratch memory
struct ArenaMarker {
  Arena &arena;
  size_t saved_offset;

  ArenaMarker(Arena &a) : arena(a), saved_offset(a.get_offset()) {}
  ~ArenaMarker() { arena.set_offset(saved_offset); }
};

// Ping-Pong Context for Double Buffering Allocations
struct ScratchContext {
  Arena *source;
  Arena *target;

  ScratchContext(Arena &a, Arena &b) : source(&a), target(&b) {
    source->set_offset(0);
    target->set_offset(0);
  }

  // Swaps the arenas and completely wipes the old source
  void swap_and_clear() {
    Arena *temp = source;
    source = target;
    target = temp;
    target->set_offset(0);
  }
};

// ============================================================================
// 2. Arena Structures
// ============================================================================

template <typename T> class ArenaVector {
private:
  T *data_;
  size_t size_;
  size_t capacity_;

public:
  ArenaVector() : data_(nullptr), size_(0), capacity_(0) {}

  // Disable implicit shallow copying to prevent memory aliasing bugs
  ArenaVector(const ArenaVector &) = delete;
  ArenaVector &operator=(const ArenaVector &) = delete;

  void shallow_copy(const ArenaVector &other) {
    data_ = other.data_;
    size_ = other.size_;
    capacity_ = other.capacity_;
  }

  // Allow move semantics
  ArenaVector(ArenaVector &&) = default;
  ArenaVector &operator=(ArenaVector &&) = default;

  // Constructor requires exact capacity upfront. No dynamic growth.
  ArenaVector(Arena &arena, size_t exact_capacity) {
    initialize(arena, exact_capacity);
  }

  void initialize(Arena &arena, size_t exact_capacity) {
    if (capacity_ >= exact_capacity && data_ != nullptr) {
      size_ = 0; // Reuse memory
      return;
    }
    if (exact_capacity > 0) {
      data_ = static_cast<T *>(
          arena.allocate(exact_capacity * sizeof(T), alignof(T)));
    } else {
      data_ = nullptr;
    }
    size_ = 0;
    capacity_ = exact_capacity;
  }

  void push_back(const T &value) {
    assert(size_ < capacity_ && "ArenaVector exact capacity exceeded!");
    // Placement new in case T has constructors
    new (&data_[size_]) T(value);
    size_++;
  }

  template <typename... Args> T &emplace_back(Args &&...args) {
    assert(size_ < capacity_ && "ArenaVector exact capacity exceeded!");
    T *ptr = new (&data_[size_]) T(std::forward<Args>(args)...);
    size_++;
    return *ptr;
  }

  T &operator[](size_t i) { return data_[i]; }
  const T &operator[](size_t i) const { return data_[i]; }

  size_t size() const { return size_; }
  size_t capacity() const { return capacity_; }
  bool empty() const { return size_ == 0; }

  void clear() {
    size_ = 0; // Does not free memory or reallocate
  }

  T *data() { return data_; }
  const T *data() const { return data_; }

  T *begin() { return data_; }
  T *end() { return data_ + size_; }
  const T *begin() const { return data_; }
  const T *end() const { return data_ + size_; }
};

// Global Arenas used throughout the effects engine
extern Arena geometry_arena;
extern Arena scratch_arena_a;
extern Arena scratch_arena_b;
