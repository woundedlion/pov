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

// --- Arena Hardware Sizing ---
constexpr size_t GEOMETRY_ARENA_SIZE = 16 * 512 * 1024; // 8 MB
constexpr size_t SCRATCH_ARENA_SIZE = 16 * 256 * 1024;  // 4 MB

// ============================================================================
// 1. Core Arena Allocator
// ============================================================================
class Arena {
  uint8_t *buffer;
  size_t capacity;
  size_t offset;

public:
  Arena(size_t size) : capacity(size), offset(0) { buffer = new uint8_t[size]; }

  ~Arena() { delete[] buffer; }

  void *allocate(size_t size, size_t align = alignof(std::max_align_t)) {
    size_t current = reinterpret_cast<size_t>(buffer + offset);
    size_t padding = (align - (current % align)) % align;
    if (offset + padding + size > capacity) {
      assert(false && "Arena out of memory!");
      return nullptr;
    }
    offset += padding;
    void *ptr = buffer + offset;
    offset += size;
    return ptr;
  }

  size_t get_offset() const { return offset; }
  void set_offset(size_t new_offset) { offset = new_offset; }
  void reset() { offset = 0; }
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

  // Dangerous: Forces size without initialization. Appropriate only if capacity
  // handles it and memory is immediately overwritten.
  void force_set_size(size_t new_size) {
    if (new_size <= capacity_) {
      size_ = new_size;
    }
  }

  T *data() { return data_; }
  const T *data() const { return data_; }

  T *begin() { return data_; }
  T *end() { return data_ + size_; }
  const T *begin() const { return data_; }
  const T *end() const { return data_ + size_; }
};

// Custom Hash Traits
template <typename K> struct ArenaHash {
  size_t operator()(const K &key) const;
};

// Speciality for Edge Lookups
template <> struct ArenaHash<std::pair<int, int>> {
  size_t operator()(const std::pair<int, int> &p) const {
    return (size_t)(p.first * 73856093 ^ p.second * 19349663);
  }
};

// Pointer Hash
template <typename T> struct ArenaHash<T *> {
  size_t operator()(T *const &ptr) const {
    // Basic pointer hash
    size_t x = reinterpret_cast<size_t>(ptr);
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return static_cast<size_t>(x ^ (x >> 31));
  }
};

// Generic fallback integer hash
template <> struct ArenaHash<uint32_t> {
  size_t operator()(const uint32_t &key) const {
    uint32_t k = key;
    k ^= k >> 16;
    k *= 0x85ebca6b;
    k ^= k >> 13;
    k *= 0xc2b2ae35;
    k ^= k >> 16;
    return (size_t)k;
  }
};

// Open-addressing Flat Hash Map
template <typename K, typename V> class ArenaMap {
private:
  struct Entry {
    K key;
    V value;
    bool is_occupied = false;
  };

  Entry *table_;
  size_t capacity_;
  size_t size_;

  size_t hash(const K &key) const { return ArenaHash<K>()(key); }

public:
  ArenaMap() : table_(nullptr), capacity_(0), size_(0) {}

  ArenaMap(Arena &arena, size_t max_elements) {
    initialize(arena, max_elements);
  }

  void initialize(Arena &arena, size_t max_elements) {
    // Size table for a ~50% load factor to avoid clustering
    capacity_ = max_elements * 2;
    if (capacity_ > 0) {
      table_ = static_cast<Entry *>(
          arena.allocate(capacity_ * sizeof(Entry), alignof(Entry)));
      for (size_t i = 0; i < capacity_; ++i)
        table_[i].is_occupied = false;
    } else {
      table_ = nullptr;
    }
    size_ = 0;
  }

  V &operator[](const K &key) {
    assert(capacity_ > 0 && "ArenaMap not initialized!");
    size_t idx = hash(key) % capacity_;

    while (table_[idx].is_occupied) {
      if (table_[idx].key == key)
        return table_[idx].value;
      idx = (idx + 1) % capacity_;
    }

    table_[idx].key = key;
    new (&table_[idx].value) V();
    table_[idx].is_occupied = true;

    size_++;
    assert(size_ <= capacity_ / 2 && "ArenaMap load factor exceeded!");
    return table_[idx].value;
  }

  bool contains(const K &key) const {
    if (capacity_ == 0)
      return false;
    size_t idx = hash(key) % capacity_;
    size_t start_idx = idx;

    while (table_[idx].is_occupied) {
      if (table_[idx].key == key)
        return true;
      idx = (idx + 1) % capacity_;
      if (idx == start_idx)
        break; // Traversed full array
    }
    return false;
  }
};

// Global Arenas used throughout the effects engine
extern Arena geometry_arena;
extern Arena scratch_arena_a;
extern Arena scratch_arena_b;
