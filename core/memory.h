/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <new>
#include <cassert>
#include <utility>
#include <concepts>
#include "platform.h"

// Device/simulator arena budget is 335 KB. The native unit-test build
// (HS_TEST_BUILD) widens it so the effect smoke harness can exercise every
// effect's full render path regardless of per-effect arena tuning; the device
// footprint is unchanged. NOTE: the native harness is a 64-bit build, so
// per-effect arena footprints measured there can be LARGER than on the 32-bit
// device wherever a pooled struct embeds a POINTER (8 B host / 4 B device:
// ArenaVector's data ptr, Fn's callable ptr, BakedPalette::lut_). Index/count
// widths no longer diverge — StaticCircularBuffer uses uint32_t, so the big
// pooled structs (Particle/VectorTrail/OrientationTrail) are byte-identical on
// both builds. Pointer-bearing container HEADERS are not multiplied across pool
// elements, so the residual host/device delta is small but nonzero; do not
// treat the host high-water mark as an exact device figure. Effects tune their
// own split via configure_arenas() to fit the device budget.
#ifdef HS_TEST_BUILD
constexpr size_t GLOBAL_ARENA_SIZE = 8 * 1024 * 1024;
#else
constexpr size_t GLOBAL_ARENA_SIZE = 335 * 1024;
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
    // `padding` is a byte count derived from the TRUE address (buffer+offset),
    // which is the robust form (correct for any base alignment, including the
    // arbitrary bases configure_arenas() rebinds to). The bounds check adds that
    // same byte count to `offset` before comparing against `capacity` (the byte
    // length from `buffer`): the allocation occupies
    // [buffer+offset+padding, buffer+offset+padding+size) and the arena is
    // [buffer, buffer+capacity), so it fits iff offset+padding+size <= capacity.
    // Mixing address-derived padding with the offset-space bound is therefore
    // correct, not a unit error.
    size_t current = reinterpret_cast<size_t>(buffer + offset);
    size_t padding = (align - (current % align)) % align;
    // No-wrap bounds check. An additive form (offset + padding + size > capacity)
    // wraps if `size` is colossal — e.g. an overflowed exact_capacity * sizeof(T)
    // from bind() — and would spuriously pass. The subtractive form can't wrap:
    // offset <= capacity is an invariant, so (capacity - offset) is safe; guard
    // the padding step against it, then compare the remaining room against size.
    if (padding > capacity - offset || size > capacity - offset - padding) {
      // Route through hs::log (Serial on device, stdout on host) instead of a
      // raw Serial.printf/printf split: hs::log formats with vsnprintf into a
      // fixed stack buffer, so the OOM diagnostic pulls in no extra stdio that
      // the device build's NDEBUG strategy works to keep out of the image.
      hs::log("[OOM] Arena: req %zu, offset %zu / cap %zu", size,
              offset + padding, capacity);
      // Over-allocation is a sizing/logic bug, not a recoverable runtime
      // condition: there is no valid recovery, and returning null only relocates
      // the corruption into whatever derefs the result (callers do not check).
      // Trap at the violation site so it's caught on the bench. The log above
      // remains for the diagnostic; HS_CHECK survives NDEBUG on device.
      HS_CHECK(false);
    }
    offset += padding;
    void *ptr = buffer + offset;
    offset += size;
    if (offset > high_water_mark)
      high_water_mark = offset;
    return ptr;
  }

  size_t get_offset() const { return offset; }
  size_t get_capacity() const { return capacity; }
  size_t get_high_water_mark() const { return high_water_mark; }

  /// True iff [ptr, ptr + bytes) lies entirely within the arena's currently
  /// live region [buffer, buffer + offset). Cheap pointer-range test (no
  /// generation stamp), available in every build. Used by ArenaVector::bind()
  /// to catch a stale reuse on the NDEBUG path, where the debug-only generation
  /// tracking is absent: a block left behind by a reset sits beyond the rewound
  /// offset and fails this check.
  bool owns_live(const void *ptr, size_t bytes) const {
    auto p = static_cast<const uint8_t *>(ptr);
    const uint8_t *end = buffer + offset;
    // Subtractive bound (end - p) avoids a p + bytes pointer overflow; the
    // p <= end guard keeps that difference non-negative so a block sitting
    // beyond a rewound offset is rejected rather than wrapping to a huge size.
    return p >= buffer && p <= end && bytes <= static_cast<size_t>(end - p);
  }

  void set_offset(size_t new_offset) {
    // The allocator's core invariant is offset <= capacity: the no-wrap bounds
    // math in allocate() (capacity - offset) is only safe while it holds. This
    // is the one seam that can break it, so trap an out-of-range rewind at the
    // source rather than letting it silently corrupt a later allocation.
    HS_CHECK(new_offset <= capacity);
    offset = new_offset;
  }

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
    // The triangular layout is only valid for an ordered, in-range pair: a
    // swapped pair silently aliases the wrong bit (dedup corruption) and an
    // out-of-range one writes adjacent memory. This is the one memory-writing
    // primitive in memory.h with no guard; assert (stripped on device under
    // NDEBUG, caught on the bench) so an unordered caller is fixed at the source.
    assert(small >= 0 && small < large && large < MAX_V);
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

/// Fixed-capacity, arena-backed vector. Move-only; no dynamic growth.
///
/// ELEMENT DESTRUCTOR CONTRACT: ArenaVector deliberately does NOT run element
/// destructors. clear(), move, move-assign, and going out of scope all leave
/// stored elements un-destructed — storage is owned and reclaimed by the arena
/// (reset/compaction), not by this handle. This is a hard requirement of the
/// arena model, not an oversight: an arena can be reset out from under a still-
/// live ArenaVector (see bind()'s stale-binding handling), so a handle-driven
/// destructor pass could run on already-reclaimed memory. Consequently, only
/// store types whose destructor need not run for correctness — trivially-
/// destructible PODs, or types like Fn/std::function whose stored captures are
/// themselves trivial (which never own external resources here). A type that
/// owns heap/handles outside the arena must not be stored in an ArenaVector.
template <typename T> class ArenaVector {
  // ArenaSpan borrows our backing data and (in debug builds) our arena
  // generation stamp for its own use-after-free check.
  template <typename U> friend class ArenaSpan;

private:
  T *data_;
  size_t size_;
  size_t capacity_;
  bool bound_ = false;
#ifndef NDEBUG
  Arena *source_arena_ = nullptr;
  uint32_t birth_generation_ = 0;
  // Bumps on every fresh allocation in bind() (the grow / re-bind path). A
  // grow swaps data_ for a new block WITHOUT touching the arena generation,
  // so this is the only signal an ArenaSpan can use to detect that its
  // snapshotted pointer was abandoned by a re-grow of the source vector.
  uint32_t rebind_generation_ = 0;

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

  ArenaVector &operator=(ArenaVector &&other) noexcept {
    if (this != &other) {
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
#ifndef NDEBUG
    // If the previous binding is STALE — the source arena was reset (generation
    // bumped) since we bound, or it's a different arena — the old storage is
    // already dead. Drop the stale binding so we rebind cleanly (a fresh
    // allocation) instead of reusing a dangling pointer. (Debug-only:
    // generation tracking exists only here; NDEBUG can't detect this, but it
    // re-allocates on a grow anyway, so device behavior is unchanged.)
    if (bound_ && (source_arena_ != &arena ||
                   birth_generation_ != arena.get_generation())) {
      bound_ = false;
    }
#endif
    // Same arena, still live, and big enough → reuse the block in place.
    if (bound_ && capacity_ >= exact_capacity) {
#ifdef NDEBUG
      // Release-build staleness guard. The debug branch above caught a stale
      // binding via the generation stamp and dropped it; NDEBUG has no stamp, so
      // without this the reuse below would silently hand back a pointer into
      // reset/rebound arena memory — the device-only buffer aliasing the
      // project's fail-fast doctrine exists to surface, not mask. owns_live is a
      // cheap pointer-range check on this cold bind() path (no per-vector state):
      // a block left behind by a reset sits beyond the rewound offset and trips
      // it. (capacity_ == 0 implies data_ == nullptr — nothing to alias.)
      HS_CHECK(capacity_ == 0 || arena.owns_live(data_, capacity_ * sizeof(T)),
               "ArenaVector::bind() reused a stale block (arena reset since bind)");
#endif
      size_ = 0; // Reuse memory
      return;
    }
    // Otherwise (unbound, or a deliberate grow that abandons the old block) →
    // allocate fresh. A re-bind that grows is NOT a memory-safety invariant
    // violation: it just leaks the old block until the next arena reset /
    // compaction reclaims it — a supported pattern (e.g. HankinSolids' shape
    // morph rebinds restored vectors to the next, larger shape between
    // compactions). So no double-bind assert here. The genuine hard failures
    // are still trapped: OOM in Arena::allocate (HS_CHECK), capacity overflow in
    // push_back/append_bulk (HS_CHECK), and use-after-free in check_alive().
#ifndef NDEBUG
    // Reaching here with a still-live binding means a genuine grow against the
    // same arena/generation (the in-place reuse above already returned, and a
    // stale binding was dropped to bound_=false). The old block is abandoned
    // until the next reset — supported, but otherwise observable only via the
    // arena high-water mark. Log it so the silent consumption is visible,
    // matching the project's fail-loud habit. Debug-only: this is a development
    // signal, not a device-release concern.
    if (bound_)
      hs::log("ArenaVector grow abandons %zu bytes (cap %zu -> %zu)",
              capacity_ * sizeof(T), capacity_, exact_capacity);
#endif
    if (exact_capacity > 0) {
      // A capacity so large that exact_capacity * sizeof(T) overflows size_t
      // would wrap to a small byte count that slips past Arena::allocate's
      // bounds check — a silent under-allocation. Trap the overflow here, where
      // the multiply happens. Cold path (binding, not per-element).
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
    // A fresh block: any span taken before this point now dangles. Bump so its
    // check_alive() trips (the arena generation alone won't — a grow leaves it
    // untouched).
    rebind_generation_++;
#endif
  }

  bool is_bound() const { return bound_; }

  void push_back(const T &value) {
    check_alive();
    check_bound();
    // Cold-path capacity guard: a fixed-capacity overflow is a sizing bug.
    // HS_CHECK survives NDEBUG so this traps in the device build instead of
    // silently writing past the arena allocation.
    HS_CHECK(size_ < capacity_, "ArenaVector exact capacity exceeded!");
    new (&data_[size_]) T(value);
    size_++;
  }

  /// Bulk-append from a contiguous source. T must be trivially copyable.
  void append_bulk(const T *src, size_t count) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "append_bulk memcpy's the source; T must be trivially copyable");
    check_alive();
    check_bound();
    HS_CHECK(size_ + count <= capacity_,
             "ArenaVector bulk append exceeds capacity!");
    memcpy(static_cast<void*>(data_ + size_), src, count * sizeof(T));
    size_ += count;
  }

  template <typename... Args> T &emplace_back(Args &&...args) {
    check_alive();
    check_bound();
    HS_CHECK(size_ < capacity_, "ArenaVector exact capacity exceeded!");
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
    check_alive();
    // No check_bound() here on purpose: clear() only resets size_ (it neither
    // frees nor touches data_), so it is a defined no-op on an unbound vector.
    // MeshState::clear() relies on this to reset members that may never have
    // been bound.
    size_ = 0; // Does not free memory or reallocate
  }

  // No check_bound() on data()/begin()/end() on purpose: an unbound (or moved-
  // from) vector is data_==nullptr with size_==0, i.e. a well-defined EMPTY
  // range. Callers rely on that idiom — mesh.h/spatial.h's *_data() accessors,
  // std::span(vertices.data(), size()), and std::sort(data(), data()+count) all
  // pass these on a size-0 vector. A bound-precondition here would false-fire on
  // correct code. The use-after-free guard (check_alive) still applies.
  T *data() {
    check_alive();
    return data_;
  }
  const T *data() const {
    check_alive();
    return data_;
  }

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
// 4. Non-Owning Span (Explicit Borrow)
// ============================================================================

/// A read-only, non-owning view into arena-allocated data.
/// Makes the distinction between owned (ArenaVector) and borrowed data
/// visible at the type level.
///
/// LIFETIME CONTRACT: a span snapshots its source vector's data_ pointer at
/// construction. In debug builds two independent stamps fault on a stale span:
/// the arena generation catches an arena RESET out from under the span, and the
/// source vector's per-vector rebind counter catches a bind()-driven RE-GROW (a
/// grow allocates a fresh block and rebinds data_ WITHOUT bumping the arena
/// generation, so the rebind counter is the only signal). Re-taking the span
/// after a grow is still the right pattern; the counter just converts a silent
/// dangle into a debug fault.
template <typename T> class ArenaSpan {
  const T *data_;
  size_t size_;
#ifndef NDEBUG
  Arena *source_arena_ = nullptr;
  uint32_t birth_generation_ = 0;
  const ArenaVector<T> *source_vec_ = nullptr;
  uint32_t source_rebind_generation_ = 0;

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
  void check_alive() const {}
#endif

public:
  ArenaSpan() : data_(nullptr), size_(0) {}

  /// Can only be created FROM an ArenaVector (explicit borrow). In debug builds
  /// the span inherits the vector's arena-generation stamp, so accessing it
  /// after the arena is reset/compacted trips the same use-after-free check
  /// the vector itself has.
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

  /// Borrowing from a temporary ArenaVector would leave the data pointer and
  /// (in debug) source_vec_ dangling the moment the temporary dies. Forbid it.
  explicit ArenaSpan(const ArenaVector<T> &&) = delete;

  const T &operator[](size_t i) const {
    check_alive();
    assert(i < size_);
    return data_[i];
  }
  size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }
  const T *data() const {
    check_alive();
    return data_;
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

extern Arena scratch_arena_a;
extern Arena scratch_arena_b;

extern Arena persistent_arena;

void configure_arenas(size_t persistent, size_t scratch_a, size_t scratch_b);
void configure_arenas_default();

// ============================================================================
// 5. ScratchScope — RAII Guard + Factory for Temporary Memory
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
// 6. RAII Arena Evacuator
// ============================================================================

/// Concept: T must provide static void clone(const T&, T&, Arena&)
template <typename T>
concept Cloneable = requires(const T &src, T &dst, Arena &arena) {
  { T::clone(src, dst, arena) } -> std::same_as<void>;
};

/// Safely evacuates an object from the persistent arena to a scratch arena,
/// and automatically restores it upon destruction.
///
/// Usage:
///   {
///     Persist<MeshState> p(live_mesh, scratch_arena_a, persistent_arena);
///     persistent_arena.reset();
///   }  // ~Persist clones backup back into persistent
template <Cloneable T> class Persist {
  T &target_;
  Arena &persistent_;

  // Declaration order matters! scratch_ must be declared BEFORE backup_
  // so that backup_ is destroyed before the scratch arena is rolled back.
  ScratchScope scratch_;
  T backup_;

public:
  Persist(T &target, Arena &scratch, Arena &persistent)
      : target_(target), persistent_(persistent), scratch_(scratch) {
    T::clone(target_, backup_, scratch_.get_arena());
  }

  ~Persist() {
    target_ = T();
    T::clone(backup_, target_, persistent_);
  }

  // Non-copyable
  Persist(const Persist &) = delete;
  Persist &operator=(const Persist &) = delete;
};
