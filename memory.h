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

#ifdef __EMSCRIPTEN__
constexpr size_t GLOBAL_ARENA_SIZE = 345 * 1024;
#else
// Teensy 4.0 — single contiguous block in RAM1 (DTCM)
constexpr size_t GLOBAL_ARENA_SIZE = 345 * 1024;
#endif

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
// 1b. Zero-Cost Phantom-Typed Arena Wrapper
// ============================================================================

/// Tag types — empty, zero-cost, exist only for the type system.
struct ScratchFrontTag {};
struct ScratchBackTag {};
struct PersistentTag {};

/// Compile-time arena discriminator. Wraps Arena& with a phantom tag so the
/// type system prevents passing the wrong arena. Implicitly converts to Arena&
/// for backward compatibility — opt into enforcement by declaring
/// TypedArena<SpecificTag> parameters in function signatures.
template <typename Tag> class TypedArena {
  Arena &inner_;

public:
  explicit TypedArena(Arena &a) : inner_(a) {}

  // Implicit conversion — existing code that takes Arena& just works
  operator Arena &() { return inner_; }
  operator const Arena &() const { return inner_; }

  // Explicit escape hatch
  Arena &raw() { return inner_; }
  const Arena &raw() const { return inner_; }

  // Forwarded Arena interface
  void *allocate(size_t size, size_t align = alignof(std::max_align_t)) {
    return inner_.allocate(size, align);
  }
  size_t get_offset() const { return inner_.get_offset(); }
  size_t get_capacity() const { return inner_.get_capacity(); }
  size_t get_high_water_mark() const { return inner_.get_high_water_mark(); }
  void set_offset(size_t new_offset) { inner_.set_offset(new_offset); }
  void reset() { inner_.reset(); }
  void reset_high_water_mark() { inner_.reset_high_water_mark(); }

#ifndef NDEBUG
  uint32_t get_generation() const { return inner_.get_generation(); }
#endif
};

// Convenience aliases
using ScratchFront = TypedArena<ScratchFrontTag>;
using ScratchBack = TypedArena<ScratchBackTag>;
using Persistent = TypedArena<PersistentTag>;

// ============================================================================
// 2. Arena Structures
// ============================================================================

template <typename T> class ArenaVector {
private:
  T *data_;
  size_t size_;
  size_t capacity_;
#ifndef NDEBUG
  Arena *source_arena_ = nullptr;
  uint32_t birth_generation_ = 0;

  void check_alive() const {
    if (source_arena_ && source_arena_->get_generation() != birth_generation_) {
      printf("USE-AFTER-FREE: ArenaVector accessed after its arena scope was "
             "destroyed!\n");
      assert(false && "ArenaVector use-after-free!");
    }
  }
#else
  void check_alive() const {}
#endif

public:
  ArenaVector() : data_(nullptr), size_(0), capacity_(0) {}

  // Disable implicit shallow copying to prevent memory aliasing bugs
  ArenaVector(const ArenaVector &) = delete;
  ArenaVector &operator=(const ArenaVector &) = delete;

  // Move semantics
  ArenaVector(ArenaVector &&other) noexcept
      : data_(other.data_), size_(other.size_), capacity_(other.capacity_) {
#ifndef NDEBUG
    source_arena_ = other.source_arena_;
    birth_generation_ = other.birth_generation_;
#endif
    other.invalidate();
  }

  ArenaVector &operator=(ArenaVector &&other) noexcept {
    if (this != &other) {
      data_ = other.data_;
      size_ = other.size_;
      capacity_ = other.capacity_;
#ifndef NDEBUG
      source_arena_ = other.source_arena_;
      birth_generation_ = other.birth_generation_;
#endif
      other.invalidate();
    }
    return *this;
  }

  // Constructor requires exact capacity upfront. No dynamic growth.
  ArenaVector(Arena &arena, size_t exact_capacity)
      : data_(nullptr), size_(0), capacity_(0) {
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
#ifndef NDEBUG
    source_arena_ = &arena;
    birth_generation_ = arena.get_generation();
#endif
  }

  /// Nullify all state. Forces re-allocation on next initialize().
  /// Required during persistent arena compaction.
  void invalidate() {
    data_ = nullptr;
    size_ = 0;
    capacity_ = 0;
#ifndef NDEBUG
    source_arena_ = nullptr;
    birth_generation_ = 0;
#endif
  }

  void push_back(const T &value) {
    check_alive();
    assert(size_ < capacity_ && "ArenaVector exact capacity exceeded!");
    new (&data_[size_]) T(value);
    size_++;
  }

  template <typename... Args> T &emplace_back(Args &&...args) {
    check_alive();
    assert(size_ < capacity_ && "ArenaVector exact capacity exceeded!");
    T *ptr = new (&data_[size_]) T(std::forward<Args>(args)...);
    size_++;
    return *ptr;
  }

  T &operator[](size_t i) {
    check_alive();
    assert(i < size_);
    return data_[i];
  }
  const T &operator[](size_t i) const {
    check_alive();
    assert(i < size_);
    return data_[i];
  }

  size_t size() const { return size_; }
  size_t capacity() const { return capacity_; }
  bool empty() const { return size_ == 0; }
  bool is_empty() const { return size_ == 0; }

  T &back() {
    check_alive();
    assert(size_ > 0);
    return data_[size_ - 1];
  }
  const T &back() const {
    check_alive();
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

struct MeshState;
struct PolyMesh;
struct CompiledHankin;
namespace MeshOps {
template <typename MeshT>
void clone(const MeshT &src, MeshT &dst, Arena &arena);
void clone(const CompiledHankin &src, CompiledHankin &dst, Arena &arena);
}

// ============================================================================
// 4. Invisible Auto-Compactor
// ============================================================================

/// RAII guard that prevents auto_compact from running while
/// raw persistent-arena pointers are held.
class CompactionLock {
  static inline int lock_count_ = 0;

public:
  CompactionLock() { lock_count_++; }
  ~CompactionLock() { lock_count_--; }

  CompactionLock(const CompactionLock &) = delete;
  CompactionLock &operator=(const CompactionLock &) = delete;

  static bool is_locked() { return lock_count_ > 0; }
};

/// RAII Guard that safely evacuates a MeshState to scratch memory, and restores
/// it to the persistent arena upon destruction. Used to protect active meshes
/// during persistent arena compaction (e.g., persistent_arena.reset()).
template <typename MeshT, typename ScratchT>
class Persist {
  MeshT *target_;
  MeshT backup_;
  Arena *scratch_;
  Arena *persistent_;

public:
  // Requires explicit TypedArena to enforce compile-time safety.
  // Passed by value because they are zero-cost wrappers often returned as temporaries.
  Persist(MeshT &target, ScratchT scratch, Persistent persistent)
      : target_(&target), scratch_(&scratch.raw()), persistent_(&persistent.raw()) {
    MeshOps::clone(target, backup_, *scratch_);
  }

  ~Persist() {
    target_->invalidate();
    MeshOps::clone(backup_, *target_, *persistent_);
  }

  // Disable copying
  Persist(const Persist &) = delete;
  Persist &operator=(const Persist &) = delete;
};

// ============================================================================
// 5. PingPong Double-Buffer
// ============================================================================

/// Encapsulates the front/back scratch arena swap pattern.
/// Provides two interfaces:
///   - front()/back()/swap(): for MemoryCtx-style usage (Conway operators)
///   - target()/source()/flip(): for iterative double-buffer loops
class PingPong {
  Arena *arenas_[2];
  int current_ = 0;

public:
  PingPong(Arena &a, Arena &b) : arenas_{&a, &b} {}

  // --- MemoryCtx-style interface ---

  /// The "active" arena (where results live).
  Arena &front() { return *arenas_[current_]; }

  /// The "other" arena (scratch / input data).
  Arena &back() { return *arenas_[1 - current_]; }

  /// Swap front/back and reset the new front.
  void swap() {
    current_ = 1 - current_;
    arenas_[current_]->reset();
  }

  // --- Iterative double-buffer interface ---

  /// Returns the arena to WRITE into, resets it first.
  Arena &target() {
    Arena &t = *arenas_[1 - current_];
    t.reset();
    return t;
  }

  /// Returns the arena to READ from (the last target).
  Arena &source() { return *arenas_[current_]; }

  /// Flip: what was the target becomes the source.
  void flip() { current_ = 1 - current_; }
};

// ============================================================================
// 6. Unified Factory
// ============================================================================
class MemoryCtx {
  PingPong pp_;

public:
  MemoryCtx() : pp_(scratch_arena_a, scratch_arena_b) {
    scratch_arena_a.reset();
    scratch_arena_b.reset();
  }
  MemoryCtx(Arena &a, Arena &b) : pp_(a, b) {
    a.reset();
    b.reset();
  }

  // Typed Arena Getters — compile-time arena discrimination
  ScratchFront get_scratch_front() { return ScratchFront(pp_.front()); }
  ScratchBack get_scratch_back() { return ScratchBack(pp_.back()); }
  Persistent get_persistent() { return Persistent(persistent_arena); }

  void swap_scratch() { pp_.swap(); }

  template <typename T> ArenaVector<T> make_scratch_front(size_t size) {
    return ArenaVector<T>(pp_.front(), size);
  }
  template <typename T> ArenaVector<T> make_scratch_back(size_t size) {
    return ArenaVector<T>(pp_.back(), size);
  }
  template <typename T> ArenaVector<T> make_persistent(size_t size) {
    return ArenaVector<T>(persistent_arena, size);
  }

  void update_persistent(MeshState &target, const PolyMesh &new_data);

  /// Static accessor for scratch front arena (no MemoryCtx instance needed).
  static ScratchFront scratch() { return ScratchFront(scratch_arena_a); }
};

struct ScopedScratch {
  Arena &arena;
  size_t saved_offset;

  ScopedScratch(Arena &a) : arena(a), saved_offset(a.get_offset()) {}
  ~ScopedScratch() { arena.set_offset(saved_offset); }
};
