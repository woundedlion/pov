/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/memory.h.

// ============================================================================
// Core Arena Allocator
// ============================================================================

/**
 * @brief Records an ArenaVector block abandoned by a move-assignment or a grow.
 * @param bytes Size of the abandoned block.
 * @details Accumulates without logging.
 * @note Cumulative across every arena modulo the size_t range, with reclaims
 * not subtracted, so it is not a live-leak figure.
 */
HS_COLD void note_arena_vector_abandon(size_t bytes);

/** @brief ArenaVector abandoned-byte count modulo the size_t range. */
FLASHMEM size_t arena_vector_abandoned_bytes();

/** @brief ArenaVector abandon-event count modulo the size_t range. */
FLASHMEM size_t arena_vector_abandon_count();

/**
 * @brief Logs an arena over-allocation then traps.
 * @param buffer Base of the arena's backing buffer.
 * @param size Bytes requested.
 * @param offset Live offset before the request.
 * @param padding Alignment padding the request needed.
 * @param capacity Arena capacity.
 * @details Also logs the ArenaVector abandon totals. Never returns.
 */
[[noreturn]] HS_COLD void arena_oom_trap(const void *buffer, size_t size,
                                         size_t offset, size_t padding,
                                         size_t capacity);

/**
 * @brief Bump allocator over a fixed caller-owned buffer.
 * @details Allocation is offset advancement; individual frees are unsupported —
 * memory is reclaimed wholesale via reset() (rewind to 0) or set_offset()
 * (rewind to a saved mark). Over-allocation traps.
 */
class Arena {
  uint8_t *buffer;
  size_t capacity;
  size_t offset;
  size_t high_water_mark;
  size_t lifetime_high_water_mark;
#ifndef NDEBUG
  uint32_t generation = 0;
  size_t rewind_floor = SIZE_MAX;
  uint64_t rewind_seq = 0;
  struct Rewind {
    uint64_t seq;
    size_t target;
  };
  static constexpr size_t REWIND_HISTORY_CAPACITY = 256;
  Rewind rewind_history[REWIND_HISTORY_CAPACITY]{};
  size_t rewind_history_size = 0;
#endif

public:
  /**
   * @brief Constructs an arena whose capacity is the whole backing buffer.
   * @param buf Pointer to the backing buffer.
   * @param size Capacity of the buffer in bytes.
   */
  constexpr Arena(uint8_t *buf, size_t size)
      : buffer(buf), capacity(size), offset(0), high_water_mark(0),
        lifetime_high_water_mark(0) {}

  /**
   * @brief Non-copyable: a copy would alias one buffer under two independent
   * offsets, handing out overlapping allocations with no trap.
   */
  Arena(const Arena &) = delete;
  Arena &operator=(const Arena &) = delete;

  /**
   * @brief Bump-allocate `size` bytes aligned to `align`, advancing the offset.
   * @param size Number of bytes to allocate.
   * @param align Required alignment in bytes; defaults to max_align_t.
   * @return Pointer into the buffer for the allocated block.
   * @details Traps via arena_oom_trap() on over-allocation. Updates the
   * high-water mark. `size` must be > 0: a zero-size block would alias the next
   * allocation.
   */
  void *allocate(size_t size, size_t align = alignof(std::max_align_t)) {
    HS_CHECK(size > 0, "Arena::allocate: zero-size request");
    HS_CHECK(align != 0 && (align & (align - 1)) == 0,
             "Arena::allocate: alignment %lu is not a power of two",
             static_cast<unsigned long>(align));
    uintptr_t current = reinterpret_cast<uintptr_t>(buffer + offset);
    size_t padding = (align - (current % align)) % align;
    // Subtractive form: cannot wrap for a colossal `size`.
    if (padding > capacity - offset || size > capacity - offset - padding)
      arena_oom_trap(buffer, size, offset, padding, capacity);
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
   * @details Does not construct the elements.
   */
  template <typename T> T *allocate_n(size_t n) {
    HS_CHECK(
        n <= SIZE_MAX / sizeof(T),
        "Arena::allocate_n element count overflows size_t: n=%lu sizeof(T)=%lu",
        static_cast<unsigned long>(n), static_cast<unsigned long>(sizeof(T)));
    return static_cast<T *>(allocate(n * sizeof(T), alignof(T)));
  }

  /**
   * @brief Allocates and constructs one object in arena storage.
   * @tparam T Object type.
   * @tparam Args Constructor argument types.
   * @param args Arguments forwarded to `T`'s constructor.
   * @return Pointer to the constructed object.
   * @details Arena reset and rewind do not run destructors. Callers that place
   * non-trivially destructible objects in an arena must end those lifetimes
   * explicitly before reclaiming their storage.
   */
  template <typename T, typename... Args> T *make(Args &&...args) {
    return ::new (static_cast<void *>(allocate_n<T>(1)))
        T(std::forward<Args>(args)...);
  }

  /**
   * @brief Allocates and default-initializes one object in arena storage.
   * @tparam T Object type.
   * @return Pointer to the constructed object.
   * @details Unlike `make<T>()`, this does not zero-initialize scalar members
   * of an aggregate before its default initialization.
   */
  template <typename T> T *make_default_initialized() {
    return ::new (static_cast<void *>(allocate_n<T>(1))) T;
  }

  /**
   * @brief Allocates and value-initializes a contiguous array.
   * @tparam T Element type.
   * @param n Element count (must be > 0, per allocate()).
   * @return Pointer to the first constructed element.
   * @details Scalar elements are zero-initialized; class members follow their
   * type's initialization rules.
   */
  template <typename T> T *make_n(size_t n) {
    T *elements = allocate_n<T>(n);
    for (size_t i = 0; i < n; ++i)
      ::new (static_cast<void *>(&elements[i])) T();
    return elements;
  }

  /**
   * @brief Allocates an array and constructs each element from an index.
   * @tparam T Element type.
   * @tparam Factory Callable returning the value for an element.
   * @param n Element count (must be > 0, per allocate()).
   * @param factory Callable receiving the zero-based element index.
   * @return Pointer to the first constructed element.
   */
  template <typename T, typename Factory>
  T *make_n_indexed(size_t n, Factory &&factory) {
    T *elements = allocate_n<T>(n);
    for (size_t i = 0; i < n; ++i)
      ::new (static_cast<void *>(&elements[i])) T(factory(i));
    return elements;
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
   * @brief Returns the peak allocation offset observed since the last
   *        reset_high_water_mark(), reset_peak_tracking() or rebind().
   * @return High-water mark in bytes.
   */
  size_t get_high_water_mark() const { return high_water_mark; }

  /**
   * @brief Returns the peak allocation offset over the arena's whole lifetime.
   * @return Largest offset any allocation has reached, in bytes.
   * @details Survives reset_high_water_mark() and rebind(); only
   * reset_peak_tracking() clears it.
   */
  size_t get_lifetime_high_water_mark() const {
    return high_water_mark > lifetime_high_water_mark
               ? high_water_mark
               : lifetime_high_water_mark;
  }

  /**
   * @brief Rewinds the offset to a previously saved mark.
   * @param new_offset Offset to rewind to; must be <= the current offset.
   * @details A forward jump traps: it would hand out bytes no allocate()
   * reserved. Alignment is not re-checked; allocate() pads from the true address.
   */
  void set_offset(size_t new_offset) {
    HS_CHECK(new_offset <= offset,
             "Arena::set_offset: %lu is not a rewind from %lu",
             static_cast<unsigned long>(new_offset),
             static_cast<unsigned long>(offset));
#ifndef NDEBUG
    if (new_offset < offset) {
      rewind_seq++;
      while (rewind_history_size &&
             rewind_history[rewind_history_size - 1].target >= new_offset)
        --rewind_history_size;
      HS_CHECK(rewind_history_size < REWIND_HISTORY_CAPACITY,
               "Arena: debug rewind history capacity exceeded");
      rewind_history[rewind_history_size++] = {rewind_seq, new_offset};
      if (new_offset < rewind_floor)
        rewind_floor = new_offset;
    }
#endif
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
    generation++;
    rewind_floor = SIZE_MAX;
    rewind_history_size = 0;
#endif
  }

  /**
   * @brief Point the arena at a different buffer/capacity and reset to empty.
   * @param buf Pointer to the new backing buffer.
   * @param new_capacity Capacity of the new buffer in bytes.
   */
  void rebind(uint8_t *buf, size_t new_capacity) {
    buffer = buf;
    capacity = new_capacity;
    offset = 0;
    fold_lifetime_peak();
    high_water_mark = 0;
#ifndef NDEBUG
    generation++;
    rewind_floor = SIZE_MAX;
    rewind_history_size = 0;
#endif
  }

  /**
   * @brief Reset windowed peak-usage tracking to the current offset.
   * @details The closed window is folded into the lifetime peak.
   */
  void reset_high_water_mark() {
    fold_lifetime_peak();
    high_water_mark = offset;
  }

  /**
   * @brief Reset both the windowed and the lifetime peak to the current offset.
   * @details Starts a measurement whose lifetime peak owes nothing to earlier
   * tenants of the arena.
   */
  void reset_peak_tracking() {
    high_water_mark = offset;
    lifetime_high_water_mark = offset;
  }

#ifndef NDEBUG
  /**
   * @brief Returns the current debug generation stamp.
   * @return Generation counter, bumped on each reset/rebind.
   */
  uint32_t get_generation() const { return generation; }

  /**
   * @brief Tests whether a byte region still lies within the live extent.
   * @param p First byte of the region.
   * @param bytes Region length in bytes.
   * @return True iff [p, p+bytes) falls within [buffer, buffer+offset).
   * @details A set_offset() rewind reclaims bytes without bumping the
   * generation.
   */
  bool covers(const void *p, size_t bytes) const {
    uintptr_t base = reinterpret_cast<uintptr_t>(buffer);
    uintptr_t q = reinterpret_cast<uintptr_t>(p);
    return q >= base && (q - base) <= offset && bytes <= offset - (q - base);
  }

  /**
   * @brief Returns the lowest offset any rewind has dropped to this generation.
   * @return Rewind floor in bytes, or SIZE_MAX if no rewind has happened since
   *         the last reset/rebind.
   */
  size_t get_rewind_floor() const { return rewind_floor; }

  /**
   * @brief Returns the count of rewinds this arena has performed.
   * @return Monotone counter, bumped by every set_offset() that lowers the
   *         offset. Never reset, so a stamp taken before a reset/rebind stays
   *         distinguishable.
   */
  uint64_t get_rewind_seq() const { return rewind_seq; }

  /**
   * @brief Tests whether a rewind reclaimed a byte region after it was handed
   *        out.
   * @param p First byte of the region.
   * @param bytes Region length in bytes.
   * @param birth_seq get_rewind_seq() sampled when the region was handed out.
   * @return True iff a rewind since those samples dropped the offset below the
   *         region's end.
   * @details The debug history retains suffix-minimum rewind targets. It traps
   * on more than REWIND_HISTORY_CAPACITY increasing targets without an
   * intervening deeper rewind.
   */
  bool reclaimed_since(const void *p, size_t bytes, uint64_t birth_seq) const {
    uintptr_t base = reinterpret_cast<uintptr_t>(buffer);
    uintptr_t q = reinterpret_cast<uintptr_t>(p);
    if (q < base)
      return false;
    size_t start = static_cast<size_t>(q - base);
    for (size_t i = 0; i < rewind_history_size; ++i)
      if (rewind_history[i].seq > birth_seq)
        return cuts_region(rewind_history[i].target, start, bytes);
    return false;
  }
#endif

private:
  friend void resplit_arenas(size_t persistent, size_t scratch_a,
                             size_t scratch_b);

#ifndef NDEBUG
  /** @brief Whether an offset of @p floor leaves [start, start+bytes) above
   *         it. */
  static bool cuts_region(size_t floor, size_t start, size_t bytes) {
    // Subtractive form: `start + bytes` could wrap for a colossal `bytes`.
    return floor < start || floor - start < bytes;
  }
#endif

  /** @brief Carries the window about to be discarded into the lifetime peak. */
  void fold_lifetime_peak() {
    if (high_water_mark > lifetime_high_water_mark)
      lifetime_high_water_mark = high_water_mark;
  }

  /**
   * @brief Moves the capacity boundary while preserving
   * base, offset, content, and generation.
   * @param new_capacity New capacity in bytes; must be >= the live
   *        offset.
   * @details The caller must have vacated whatever else held those bytes.
   */
  void rebind_capacity(size_t new_capacity) {
    HS_CHECK(offset <= new_capacity,
             "Arena::rebind_capacity below the live offset would strand "
             "content");
    capacity = new_capacity;
  }
};

#ifndef NDEBUG
/**
 * @brief Debug-only snapshot of an arena's state when it handed out a block.
 * @details Compiled out under NDEBUG.
 */
struct ArenaBlockStamp {
  Arena *source_arena = nullptr; /**< Arena the block was allocated from. */
  uint32_t birth_generation = 0; /**< Arena generation when stamped. */
  uint64_t birth_rewind_seq = 0; /**< Arena rewind counter when stamped. */

  /**
   * @brief Stamps against @p arena's current generation and rewind counter.
   * @param arena Arena the block was just allocated from.
   */
  void record(Arena &arena) {
    source_arena = &arena;
    birth_generation = arena.get_generation();
    birth_rewind_seq = arena.get_rewind_seq();
  }

  /**
   * @brief Drops the stamp, leaving the owner untracked.
   */
  void clear() {
    source_arena = nullptr;
    birth_generation = 0;
    birth_rewind_seq = 0;
  }

  /**
   * @brief Whether the source arena was reset or rebound since the stamp.
   * @return True iff the block's bytes have been reissued wholesale.
   */
  bool arena_reset() const {
    return source_arena && source_arena->get_generation() != birth_generation;
  }

  /**
   * @brief Whether a rewind dropped the arena's live extent below the block.
   * @param p First byte of the block.
   * @param bytes Block length in bytes; 0 owns no storage and never faults.
   * @return True iff the region no longer lies within the live extent.
   */
  bool block_uncovered(const void *p, size_t bytes) const {
    return source_arena && bytes > 0 && !source_arena->covers(p, bytes);
  }

  /**
   * @brief Whether a rewind reclaimed the block and later allocations re-covered
   *        it.
   * @param p First byte of the block.
   * @param bytes Block length in bytes; 0 owns no storage and never faults.
   * @return True iff a rewind since the stamp freed the region.
   * @details Catches the window block_uncovered() goes blind in — the point at
   * which a second owner starts writing the same bytes.
   */
  bool block_reissued(const void *p, size_t bytes) const {
    return source_arena && bytes > 0 &&
           source_arena->reclaimed_since(p, bytes, birth_rewind_seq);
  }

  /** @brief Whether the stamped block remains owned and live. */
  bool block_alive(const void *p, size_t bytes) const {
    return !arena_reset() && !block_uncovered(p, bytes) &&
           !block_reissued(p, bytes);
  }
};

/**
 * @brief Faults when an arena-owned block has been reset, rewound below or
 *        reissued since it was stamped.
 * @param stamp ArenaBlockStamp recorded when the block was allocated.
 * @param ptr First byte of the block.
 * @param bytes Block length in bytes.
 * @param owner String literal naming the owner in the failure message.
 * @details Expands to nothing under NDEBUG, where the stamp does not exist.
 */
#define HS_ASSERT_BLOCK_ALIVE(stamp, ptr, bytes, owner)                        \
  assert((stamp).block_alive(ptr, bytes) && owner " use-after-free!")
#else
#define HS_ASSERT_BLOCK_ALIVE(stamp, ptr, bytes, owner) ((void)0)
#endif

extern Arena scratch_arena_a;
extern Arena scratch_arena_b;

extern Arena persistent_arena;

/**
 * @brief Self-registering callback run before persistent arena storage is
 *        handed out again.
 * @details A global that caches a pointer into the persistent arena declares one
 * static instance and drops the pointer from the callback. The registry head is
 * constant-initialized, so registration during static init is order-independent.
 * @note Persistent arena only; no global may cache a pointer into scratch
 * storage.
 */
struct ArenaResetHook {
  using Handler = void (*)(); /**< Callback signature. */

  Handler handler;      /**< Callback invoked by run_all(). */
  ArenaResetHook *next; /**< Next link in the intrusive registry list. */

  /**
   * @brief Registers @p h with the global hook list.
   * @param h Callback that drops the owner's pointer into arena storage.
   */
  explicit ArenaResetHook(Handler h) : handler(h), next(head) { head = this; }

  /** @brief Unlinks this hook so run_all() never calls through a dead node. */
  ~ArenaResetHook() {
    for (ArenaResetHook **p = &head; *p; p = &(*p)->next) {
      if (*p == this) {
        *p = next;
        return;
      }
    }
    HS_CHECK(false, "ArenaResetHook: destroyed hook not in registry");
  }

  /** @brief Deleted copy constructor: a copy would double-link the registry. */
  ArenaResetHook(const ArenaResetHook &) = delete;
  /**
   * @brief Deleted copy assignment (non-copyable).
   * @return Reference to this (never invoked).
   */
  ArenaResetHook &operator=(const ArenaResetHook &) = delete;

  /** @brief Runs every registered hook. */
  HS_COLD_MEMBER static void run_all() {
    for (const ArenaResetHook *h = head; h; h = h->next)
      h->handler();
  }

private:
  static inline ArenaResetHook *head =
      nullptr; /**< Head of the intrusive registry list. */
};

/**
 * @brief Rewinds the persistent arena to empty after dropping every cached
 *        pointer into it.
 * @details A bare `persistent_arena.reset()` leaves each registered global
 * pointing at bytes the next allocation re-issues.
 */
HS_FLASH_INLINE inline void reset_persistent_arena() {
  ArenaResetHook::run_all();
  persistent_arena.reset();
}

/**
 * @brief Repartitions the global arena budget across the three arenas.
 * @param persistent Bytes to assign to the persistent arena.
 * @param scratch_a Bytes to assign to scratch arena A.
 * @param scratch_b Bytes to assign to scratch arena B.
 */
FLASHMEM void configure_arenas(size_t persistent, size_t scratch_a,
                               size_t scratch_b);
/** @brief Scratch capacities and the remaining persistent arena budget. */
struct ArenaSplit {
  size_t scratch_a;
  size_t scratch_b;

  constexpr size_t persistent(size_t total = DEVICE_GLOBAL_ARENA_SIZE) const {
    HS_CHECK(scratch_a <= total && scratch_b <= total - scratch_a,
             "ArenaSplit: scratch %lu+%lu exceeds %lu B",
             static_cast<unsigned long>(scratch_a),
             static_cast<unsigned long>(scratch_b),
             static_cast<unsigned long>(total));
    return total - scratch_a - scratch_b;
  }

  constexpr size_t device_persistent() const {
    return persistent(DEVICE_GLOBAL_ARENA_SIZE);
  }

  HS_COLD_MEMBER void configure() const {
    configure_arenas(persistent(GLOBAL_ARENA_SIZE), scratch_a, scratch_b);
  }
};

/**
 * @brief Restores the default arena partition.
 */
FLASHMEM void configure_arenas_default();

/**
 * @brief Re-partitions the arenas mid-run WITHOUT disturbing persistent content.
 * @param persistent New persistent capacity; must be >= its current live offset.
 * @param scratch_a New scratch-A capacity.
 * @param scratch_b New scratch-B capacity.
 * @details Unlike configure_arenas(), the persistent arena keeps its base,
 * offset, live content, and generation; only its capacity boundary moves. Both
 * scratch arenas MUST be empty; they rebind to fresh bases.
 */
FLASHMEM void resplit_arenas(size_t persistent, size_t scratch_a,
                             size_t scratch_b);
