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
constexpr size_t SCRATCH_ARENA_A_SIZE = 256 * 1024;
constexpr size_t SCRATCH_ARENA_B_SIZE = 128 * 1024;
constexpr size_t PERSISTENT_ARENA_SIZE = 256 * 1024;

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

// Ping-Pong Context for Double Buffering Allocations (legacy)
struct ScratchContext {
  Arena *source;
  Arena *target;

  ScratchContext(Arena &a, Arena &b) : source(&a), target(&b) {
    source->set_offset(0);
    target->set_offset(0);
  }

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
// FRAME BUFFERS (Wiped every frame)
extern Arena geo_arena_a;
extern Arena geo_arena_b;

// SCRATCH BUFFERS (Ping-ponged during heavy math)
extern Arena scratch_arena_a;
extern Arena scratch_arena_b;

// PERSISTENT STATE (Lives across frames. Auto-compacted)
extern Arena persistent_arena;

// ENGINE FRAME STATE
extern bool using_A_as_frame;
extern Arena *current_frame_arena;

struct MeshState;   // Forward declaration for PersistentTracker
struct PolyMesh;    // Forward declaration for PersistentTracker
namespace MeshOps { // Forward declaration for cloning
template <typename MeshT>
void clone(const MeshT &src, MeshT &dst, Arena &arena);
}

// ============================================================================
// Invisible Auto-Compactor
// ============================================================================
class PersistentTracker {
public:
  // Tracks pointers to the effect's persistent meshes (MeshState*)
  static ArenaVector<MeshState *> tracked_meshes;

  static void register_mesh(MeshState *mesh) {
    if (tracked_meshes.capacity() == 0) {
      tracked_meshes.initialize(persistent_arena, 256);
    }
    tracked_meshes.push_back(mesh);
  }

  static void clear_registry() { tracked_meshes.clear(); }

  // Called automatically by MemoryCtx when the arena gets full
  static void auto_compact(Arena &safe_scratch);
};

// ============================================================================
// Unified Factory
// ============================================================================
class MemoryCtx {
  Arena *current_scratch;
  Arena *next_scratch;

public:
  // Defaults to scratch A and B.
  MemoryCtx()
      : current_scratch(&scratch_arena_a), next_scratch(&scratch_arena_b) {}
  MemoryCtx(Arena &a, Arena &b) : current_scratch(&a), next_scratch(&b) {
    current_scratch->set_offset(0);
    next_scratch->set_offset(0);
  }

  // 1. Scratch Management
  Arena &get_scratch() { return *current_scratch; }

  void swap_scratch() {
    Arena *temp = current_scratch;
    current_scratch = next_scratch;
    next_scratch = temp;
    // The newly active scratch destination is always wiped clean
    current_scratch->reset();
  }

  template <typename T> ArenaVector<T> make_scratch_array(size_t size) {
    return ArenaVector<T>(*current_scratch, size);
  }

  // 2. Invisible Persistent Storage Updates
  // Target is MeshState
  void update_persistent(MeshState &target, const PolyMesh &new_data);
};

// 3. The RAII Guard
struct ScopedScratch {
  Arena &arena;
  size_t saved_offset;

  ScopedScratch(Arena &a) : arena(a), saved_offset(a.get_offset()) {}
  ScopedScratch(MemoryCtx &ctx)
      : arena(ctx.get_scratch()), saved_offset(arena.get_offset()) {}
  ~ScopedScratch() { arena.set_offset(saved_offset); }
};
