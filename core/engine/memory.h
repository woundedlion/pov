/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file memory.h
 * @brief Arena allocator, the engine's global arena budget, the containers
 *        built on top of it, and the scratch-scoped generate() wrapper.
 */

// platform.h defines NDEBUG on device; include before <cassert> so assert
// stripping does not depend on include order.
#include "platform/platform.h"
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <new>
#include <cassert>
#include <utility>
#include <concepts>

// Device arena budget is 298 KiB; the WASM simulator uses 512 KiB. Host effect
// harnesses use an 8 MiB configuration so they can exercise every effect's full
// render path. The
// native harness is a 64-bit build, so per-effect footprints measured there
// can be LARGER than on the 32-bit device wherever a pooled struct embeds a
// POINTER (ArenaVector's data ptr, Fn's callable ptr, BakedPalette::lut). Do not
// treat the host high-water mark as an exact device figure. Effects tune their
// own split via configure_arenas() to fit the device budget.
// The real device FlexRAM (RAM1) arena, sized from the measured worst-effect
// high-water (tests/arena_measure.cpp): GSReactionDiffusion is the binding
// tenant at ~291 KiB total (~171 KiB persistent + ~120 KiB scratch under its
// own split). A distinct always-defined constant (not the host-inflated
// GLOBAL_ARENA_SIZE below) so device-budget static_asserts check the real
// figure even in the host suite.
constexpr size_t DEVICE_GLOBAL_ARENA_SIZE = 298 * 1024;
static_assert(DEVICE_GLOBAL_ARENA_SIZE == HS_DEVICE_ARENA_BYTES,
              "device arena budget must equal the block HS_DEVICE_ARENA_BYTES "
              "sizes; moving one alone traps at Arena::allocate on hardware");
constexpr size_t GLOBAL_ARENA_SIZE = HS_GLOBAL_ARENA_BYTES;

constexpr size_t DEFAULT_SCRATCH_A_SIZE = 16 * 1024;
constexpr size_t DEFAULT_SCRATCH_B_SIZE = 16 * 1024;
// An HS_GLOBAL_ARENA_BYTES override at or below the scratch split would wrap the
// unsigned subtraction below, giving the persistent arena a capacity far larger
// than global_arena_block; its two-argument constructor passes size as its own
// extent, so nothing traps and the first allocation writes past the block.
static_assert(GLOBAL_ARENA_SIZE >
                  DEFAULT_SCRATCH_A_SIZE + DEFAULT_SCRATCH_B_SIZE,
              "HS_GLOBAL_ARENA_BYTES must exceed the default scratch split");
constexpr size_t DEFAULT_PERSISTENT_SIZE =
    GLOBAL_ARENA_SIZE - DEFAULT_SCRATCH_A_SIZE - DEFAULT_SCRATCH_B_SIZE;
// Persistent budget on the real device split (from DEVICE_GLOBAL_ARENA_SIZE, not
// the host-inflated GLOBAL_ARENA_SIZE) so an effect's default-split footprint
// static_assert checks the true device figure even in the host suite.
constexpr size_t DEVICE_PERSISTENT_BUDGET =
    DEVICE_GLOBAL_ARENA_SIZE - DEFAULT_SCRATCH_A_SIZE - DEFAULT_SCRATCH_B_SIZE;

// The browser module's arena (CMakeLists.txt), widened past the device figure
// for the chain interpreter's two arenas plus ShaderChain's shared color
// resources. A distinct always-defined constant so a browser-only effect's
// footprint static_assert checks the real module figure even in the host suite,
// whose arena is an order of magnitude larger.
constexpr size_t WASM_GLOBAL_ARENA_SIZE = 512 * 1024;
#if defined(__EMSCRIPTEN__)
static_assert(WASM_GLOBAL_ARENA_SIZE == GLOBAL_ARENA_SIZE,
              "the module's HS_GLOBAL_ARENA_BYTES moved; update "
              "WASM_GLOBAL_ARENA_SIZE");
#endif
constexpr size_t WASM_PERSISTENT_BUDGET =
    WASM_GLOBAL_ARENA_SIZE - DEFAULT_SCRATCH_A_SIZE - DEFAULT_SCRATCH_B_SIZE;

// ============================================================================
// 1. Core Arena Allocator
// ============================================================================

/**
 * @brief Records an ArenaVector move-assignment that abandoned its destination's
 *        block.
 * @param bytes Size of the abandoned block.
 * @details Accumulates instead of logging the line bind() emits for a grow:
 * move-assignment runs deep inside mesh work, where the formatter's stack frame
 * does not fit the device stack budget and one line per event would bury the
 * log. The running total is reported by the arena's OOM trap. Out-of-line and
 * non-template so the device image carries one copy for every element type.
 * @note Cumulative across every arena and never decremented, so it is an upper
 * bound on bytes dropped since boot, not a live-leak figure: the chained mesh
 * ops that dominate it rewind their arena right after each step, reclaiming
 * what was counted. Subtracting reclaims would need each block's source arena
 * in release builds, which ArenaVector tracks only in debug builds.
 */
void note_arena_vector_abandon(size_t bytes);

/** @brief Bytes abandoned so far by ArenaVector move-assignment. */
FLASHMEM size_t arena_vector_abandoned_bytes();

/** @brief Move-assignments that have abandoned a block so far. */
FLASHMEM size_t arena_vector_abandon_count();

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
  size_t lifetime_high_water_mark;
#ifndef NDEBUG
  uint32_t generation = 0;
  size_t rewind_floor = SIZE_MAX;
  size_t last_rewind_target = SIZE_MAX;
  uint32_t rewind_seq = 0;
#endif

public:
  /**
   * @brief Constructs an arena whose capacity is the whole backing buffer.
   * @param buf Pointer to the backing buffer.
   * @param size Capacity of the buffer in bytes.
   */
  Arena(uint8_t *buf, size_t size)
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
   * @details Traps (HS_CHECK) on over-allocation rather than returning null.
   * Updates the high-water mark. `size` must be > 0: a zero-size request returns a
   * bump pointer that reserves no storage (it aliases the next allocation's
   * address), so it is trapped as misuse rather than handed back as ownable.
   */
  void *allocate(size_t size, size_t align = alignof(std::max_align_t)) {
    HS_CHECK(size > 0, "Arena::allocate: zero-size request");
    HS_CHECK(align != 0 && (align & (align - 1)) == 0,
             "Arena::allocate: alignment %lu is not a power of two",
             static_cast<unsigned long>(align));
    uintptr_t current = reinterpret_cast<uintptr_t>(buffer + offset);
    size_t padding = (align - (current % align)) % align;
    // Subtractive form: offset <= capacity is invariant, so it cannot wrap the
    // way `offset + padding + size > capacity` would for a colossal `size`.
    if (padding > capacity - offset || size > capacity - offset - padding) {
      hs::log("[OOM] Arena @%08lx: req %lu, offset %lu, pad %lu / cap %lu "
              "(move-assign dropped %lu B in %lu blocks, all arenas since "
              "boot, reclaims not subtracted)",
              static_cast<unsigned long>(reinterpret_cast<uintptr_t>(buffer)),
              static_cast<unsigned long>(size),
              static_cast<unsigned long>(offset),
              static_cast<unsigned long>(padding),
              static_cast<unsigned long>(capacity),
              static_cast<unsigned long>(arena_vector_abandoned_bytes()),
              static_cast<unsigned long>(arena_vector_abandon_count()));
      HS_CHECK(false, "Arena::allocate: out of memory");
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
    HS_CHECK(n <= SIZE_MAX / sizeof(T),
             "Arena::allocate_n element count overflows size_t");
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
   * @details Scalar members are zero-filled; `make_default_initialized()` is
   * the non-zeroing single-object counterpart.
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
   * @details Survives every reset_high_water_mark() and rebind(), each of which
   * folds the window it discards into this figure; reset_peak_tracking() is the
   * only way to clear it. This is the figure to size a budget against: an effect
   * that re-splits the arena mid-run leaves get_high_water_mark() reporting only
   * the peak since its last re-split.
   */
  size_t get_lifetime_high_water_mark() const {
    return high_water_mark > lifetime_high_water_mark
               ? high_water_mark
               : lifetime_high_water_mark;
  }

  /**
   * @brief Rewinds the offset to a previously saved mark.
   * @param new_offset Offset to rewind to; must be <= the current offset.
   * @details A mark is only valid as a rewind target: jumping the offset *forward*
   * would hand out backing bytes never reserved by an allocate() call, so any
   * non-rewind traps. (new_offset <= offset also implies new_offset <= capacity,
   * preserving the no-wrap bounds math in allocate().) Alignment is not re-checked:
   * allocate() recomputes leading padding from the true address on every call, so
   * restoring an unaligned mark is safe.
   */
  void set_offset(size_t new_offset) {
    HS_CHECK(new_offset <= offset,
             "Arena::set_offset: %lu is not a rewind from %lu",
             static_cast<unsigned long>(new_offset),
             static_cast<unsigned long>(offset));
#ifndef NDEBUG
    if (new_offset < offset) {
      last_rewind_target = new_offset;
      rewind_seq++;
    }
    if (new_offset < rewind_floor)
      rewind_floor = new_offset;
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
    last_rewind_target = SIZE_MAX;
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
    fold_lifetime_peak();
    high_water_mark = 0;
#ifndef NDEBUG
    generation++;
    rewind_floor = SIZE_MAX;
    last_rewind_target = SIZE_MAX;
#endif
  }

  /**
   * @brief Reset windowed peak-usage tracking to the current offset.
   * @details E.g. to measure a single frame's allocation peak in isolation. The
   * window being closed is folded into the lifetime peak, which is unaffected.
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
   * generation, so this is the only signal an ArenaSpan has that its borrowed
   * region was freed by a rewind of the source arena.
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
  uint32_t get_rewind_seq() const { return rewind_seq; }

  /**
   * @brief Tests whether a rewind reclaimed a byte region after it was handed
   *        out.
   * @param p First byte of the region.
   * @param bytes Region length in bytes.
   * @param birth_floor get_rewind_floor() sampled when the region was handed
   *        out.
   * @param birth_seq get_rewind_seq() sampled when the region was handed out.
   * @return True iff a rewind since those samples dropped the offset below the
   *         region's end.
   * @details covers() goes blind the moment fresh allocations re-cover the
   * reclaimed bytes, which is exactly when a second owner starts writing them.
   * Two independent signals, since neither alone spans every rewind: the floor
   * only ever falls within a generation, so a floor below @p birth_floor is
   * proof a rewind after the sample set it; and a bumped sequence proves the
   * most recent rewind is itself post-sample, which catches the rewind that
   * frees the region without reaching a floor set by an earlier, deeper one.
   */
  bool reclaimed_since(const void *p, size_t bytes, size_t birth_floor,
                       uint32_t birth_seq) const {
    const bool floor_fell = rewind_floor < birth_floor;
    const bool rewound_since = rewind_seq != birth_seq;
    if (!floor_fell && !rewound_since)
      return false;
    uintptr_t base = reinterpret_cast<uintptr_t>(buffer);
    uintptr_t q = reinterpret_cast<uintptr_t>(p);
    if (q < base)
      return false;
    size_t start = static_cast<size_t>(q - base);
    if (floor_fell && cuts_region(rewind_floor, start, bytes))
      return true;
    return rewound_since && cuts_region(last_rewind_target, start, bytes);
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
   * @details The caller must have vacated whatever
   * else held those bytes. resplit_arenas() alone reaches it — it re-bases both
   * scratch arenas onto the new split and bounds the request against the global
   * block.
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
 * @details One copy per arena-resident owner (ArenaVector, ArenaSpan,
 * TransformerPool), so the three lifetime questions a block can be asked —
 * was the arena reset, was it rewound below the block, was the block reclaimed
 * by a rewind and reissued — have a single set of answers. Compiled out under
 * NDEBUG along with the Arena accessors it calls.
 */
struct ArenaBlockStamp {
  Arena *source_arena = nullptr; /**< Arena the block was allocated from. */
  uint32_t birth_generation = 0; /**< Arena generation when stamped. */
  size_t birth_rewind_floor = 0; /**< Arena rewind floor when stamped. */
  uint32_t birth_rewind_seq = 0; /**< Arena rewind counter when stamped. */

  /**
   * @brief Stamps against @p arena's current generation, rewind floor and
   *        rewind counter.
   * @param arena Arena the block was just allocated from.
   */
  void record(Arena &arena) {
    source_arena = &arena;
    birth_generation = arena.get_generation();
    birth_rewind_floor = arena.get_rewind_floor();
    birth_rewind_seq = arena.get_rewind_seq();
  }

  /**
   * @brief Drops the stamp, leaving the owner untracked.
   */
  void clear() {
    source_arena = nullptr;
    birth_generation = 0;
    birth_rewind_floor = 0;
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
           source_arena->reclaimed_since(p, bytes, birth_rewind_floor,
                                         birth_rewind_seq);
  }

  /** @brief Whether the stamped block remains owned and live. */
  bool block_alive(const void *p, size_t bytes) const {
    return !arena_reset() && !block_uncovered(p, bytes) &&
           !block_reissued(p, bytes);
  }
};
#endif

// ============================================================================
// 2. Arena Structures
// ============================================================================

/**
 * @brief Logs an ArenaVector grow that abandons its previous block.
 * @param bytes Size of the abandoned block.
 * @param old_capacity Element capacity before the grow.
 * @param new_capacity Element capacity after the grow.
 * @details Out-of-line and non-template so the device image carries one copy for
 * every element type. The leak is permanent until the arena is reset, so the
 * line ships in release: without it a persistent-arena grow surfaces only as a
 * later, innocent-looking allocation trapping on OOM.
 */
FLASHMEM void log_arena_vector_grow(size_t bytes, size_t old_capacity,
                                    size_t new_capacity);

/**
 * @brief Whether T is a sanctioned inline callable safe to store in an
 * ArenaVector despite not being trivially destructible.
 * @details Only the device's teensy::inplace_function needs the exemption; the
 * host/WASM hs::inplace_function is trivially destructible and passes on its own.
 * Stored captures remain subject to ArenaVector's element destructor contract.
 */
template <typename T> struct is_arena_inplace_fn : std::false_type {};
#ifdef ARDUINO
template <typename R, typename... Args, size_t Cap, size_t Align>
struct is_arena_inplace_fn<teensy::inplace_function<R(Args...), Cap, Align>>
    : std::true_type {};
#endif

/**
 * @brief Arena-backed vector with a capacity fixed between bind() calls.
 *        Move-only.
 * @tparam T Element type; must satisfy the element destructor contract below.
 * @details CAPACITY CONTRACT: appending never grows the block — push_back()
 * traps once element_count reaches capacity(). Only bind() changes capacity,
 * and a grow there allocates a fresh block and abandons the old one until the
 * arena is reset (log_arena_vector_grow() reports the leaked bytes in release).
 * Element access follows std::vector: operator[] and back() require a valid
 * index/non-empty vector and check that precondition only in debug builds.
 *
 * ELEMENT DESTRUCTOR CONTRACT: ArenaVector does NOT run element
 * destructors — clear(), move, move-assign and going out of scope all leave
 * stored elements un-destructed. Storage is owned and reclaimed by the arena
 * (reset/compaction), not by this handle, and an arena can be reset out from
 * under a still-live ArenaVector (see bind()'s stale-binding handling), so a
 * handle-driven destructor pass could run on already-reclaimed memory. Only store
 * types whose destructor need not run for correctness: trivially-destructible
 * PODs, or Fn<> whose stored captures are themselves trivial. A type owning
 * heap/handles outside the arena must not be stored here — notably a raw
 * std::function, which would leak when this handle skips its destructor.
 */
template <typename T> class ArenaVector {
  // ArenaSpan borrows our backing data and (in debug builds) our arena
  // generation stamp for its own use-after-free check.
  template <typename U> friend class ArenaSpan;

private:
  T *elements;             /**< Pointer to the arena-allocated backing block. */
  size_t element_count;    /**< Number of constructed elements. */
  size_t element_capacity; /**< Maximum element count the block can hold. */
  bool bound = false; /**< Whether the vector has been bound to an arena. */
#ifndef NDEBUG
  ArenaBlockStamp stamp; /**< Arena state when the block was allocated. */
  /**
   * @brief Per-vector counter bumped on every fresh allocation in bind() (the
   * grow / re-bind path).
   * @details A grow swaps elements for a new block WITHOUT touching the arena
   * generation, so this is the only signal an ArenaSpan can use to detect that
   * its snapshotted pointer was abandoned by a re-grow of the source vector.
   */
  uint32_t rebind_generation = 0;

  /**
   * @brief Debug-only use-after-free check against the source arena.
   * @details Asserts if the source arena was reset out from under this vector,
   * or rewound (set_offset/ScratchScope) below the backing block — a rewind
   * reclaims the bytes without bumping the generation.
   */
  void check_alive() const {
    const size_t bytes = element_capacity * sizeof(T);
    assert(!stamp.arena_reset() && "ArenaVector use-after-free!");
    assert(!stamp.block_uncovered(elements, bytes) &&
           "ArenaVector use-after-free (arena rewound below block)!");
    assert(!stamp.block_reissued(elements, bytes) &&
           "ArenaVector use-after-free (block reclaimed by a rewind and "
           "reissued)!");
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
    assert(bound && "Attempted to access unbound ArenaVector!");
  }

  /**
   * @brief Transfers @p other's storage and bookkeeping into this vector,
   *        leaving @p other in a pristine unbound state.
   */
  void steal_from(ArenaVector &other) noexcept {
    elements = other.elements;
    element_count = other.element_count;
    element_capacity = other.element_capacity;
    bound = other.bound;
#ifndef NDEBUG
    stamp = other.stamp;
    rebind_generation = other.rebind_generation;
#endif
    other.elements = nullptr;
    other.element_count = 0;
    other.element_capacity = 0;
    other.bound = false;
#ifndef NDEBUG
    other.stamp.clear();
    // rebind_generation stays: spans snapshotted it and still view live data.
#endif
  }

public:
  using value_type = T;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using reference = T &;
  using const_reference = const T &;
  using pointer = T *;
  using const_pointer = const T *;
  using iterator = T *;
  using const_iterator = const T *;

  /**
   * @brief Default-constructs an unbound vector.
   * @details Must call bind() before use.
   */
  ArenaVector() : elements(nullptr), element_count(0), element_capacity(0) {}

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
    if (this != &other) {
      // Overwriting a bound handle drops its block; the arena reclaims it only
      // on the next reset.
      if (bound && element_capacity > 0)
        note_arena_vector_abandon(element_capacity * sizeof(T));
      steal_from(other);
    }
    return *this;
  }

  /**
   * @brief Constructs and binds the vector with an exact capacity.
   * @param arena Arena to allocate the backing block from.
   * @param exact_capacity Element count to allocate; appending never grows it.
   */
  ArenaVector(Arena &arena, size_t exact_capacity)
      : elements(nullptr), element_count(0), element_capacity(0) {
    bind(arena, exact_capacity);
  }

  /**
   * @brief Binds the vector to an arena, allocating its backing block.
   * @param arena Arena to allocate from.
   * @param min_capacity Minimum element count to reserve.
   * @details If already bound with at least that capacity, resets size for reuse
   * and keeps the larger prior capacity, so capacity() and the push_back
   * overflow guard report the block actually held rather than this request.
   * A grow against the same arena/generation reallocates a fresh block and
   * abandons the old one until the next reset; a stale binding (arena reset or
   * a different arena) trips a debug-only contract assert.
   */
  void bind(Arena &arena, size_t min_capacity) {
    static_assert(
        std::is_trivially_destructible_v<T> || is_arena_inplace_fn<T>::value,
        "ArenaVector never runs element destructors, so T must own no "
        "state outside the arena buffer: store a trivially-destructible "
        "type or a sanctioned Fn<> (no std::function/std::string).");
#ifndef NDEBUG
    // Rebinding a still-bound vector after its source arena was reset, or to a
    // different arena, is a contract violation (the old block is already dead). A
    // same-arena/same-generation grow is not stale and reallocates below.
    assert((!bound || element_capacity == 0 ||
            (stamp.source_arena == &arena &&
             stamp.birth_generation == arena.get_generation())) &&
           "ArenaVector::bind() on a stale binding: clear the handle before "
           "resetting or changing its arena");
#endif
    // The generation assert above misses a rewind, which leaves the reuse path
    // below handing back bytes the arena has already reissued.
    check_alive();
    // Same arena, still live, and big enough → reuse the block in place.
    if (bound && element_capacity >= min_capacity) {
      element_count = 0;
#ifndef NDEBUG
      // Reuse dangles any span snapshotted before this point; bump so its
      // check_alive() trips (the arena generation alone won't).
      rebind_generation++;
#endif
      return;
    }
    // Otherwise (unbound, or a grow that abandons the old block) → allocate
    // fresh. A grow leaks the old block until the next arena reset/compaction; a
    // zero-capacity binding owns no block, so growing out of one leaks nothing.
    if (bound && element_capacity > 0)
      log_arena_vector_grow(element_capacity * sizeof(T), element_capacity,
                            min_capacity);
    if (min_capacity > 0) {
      elements = arena.allocate_n<T>(min_capacity);
    } else {
      elements = nullptr;
    }
    element_count = 0;
    element_capacity = min_capacity;
    bound = true;
#ifndef NDEBUG
    if (min_capacity > 0)
      stamp.record(arena);
    rebind_generation++;
#endif
  }

  /**
   * @brief Reports whether the vector is bound to an arena.
   * @return True iff bind() has been called.
   */
  bool is_bound() const { return bound; }

  /**
   * @brief Appends a copy of an element.
   * @param value Element to copy-construct at the end.
   */
  void push_back(const T &value) {
    check_alive();
    check_bound();
    HS_CHECK(element_count < element_capacity,
             "ArenaVector push_back exact capacity exceeded!");
    new (&elements[element_count]) T(value);
    element_count++;
  }

  /**
   * @brief Bulk-append from a contiguous source.
   * @param src Pointer to the first source element.
   * @param count Number of elements to copy.
   * @details T must be trivially copyable (memcpy'd). An empty append is a
   * no-op that skips memcpy to avoid null-pointer UB.
   */
  void append_bulk(const T *src, size_t count) {
    static_assert(
        std::is_trivially_copyable_v<T>,
        "append_bulk memcpy's the source; T must be trivially copyable");
    check_alive();
    check_bound();
    // Subtractive, wrap-proof form: `element_count + count` could wrap for a
    // colossal count.
    HS_CHECK(count <= element_capacity - element_count,
             "ArenaVector bulk append exceeds capacity!");
    // Skip memcpy on an empty append: a null src with count 0 is formal UB.
    if (count == 0)
      return;
    memcpy(static_cast<void *>(elements + element_count), src,
           count * sizeof(T));
    element_count += count;
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
    HS_CHECK(element_count < element_capacity,
             "ArenaVector emplace_back exact capacity exceeded!");
    T *ptr = new (&elements[element_count]) T(std::forward<Args>(args)...);
    element_count++;
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
    assert(i < element_count);
    return elements[i];
  }
  /**
   * @brief Element access by index (const).
   * @param i Index in [0, size()).
   * @return Const reference to the element at index i.
   */
  const T &operator[](size_t i) const {
    check_alive();
    check_bound();
    assert(i < element_count);
    return elements[i];
  }

  /**
   * @brief Returns the number of stored elements.
   * @return Current element count.
   */
  size_t size() const { return element_count; }
  /**
   * @brief Returns the maximum element count.
   * @return Capacity in elements.
   */
  size_t capacity() const { return element_capacity; }
  /**
   * @brief Reports whether the vector is empty.
   * @return True iff size() == 0.
   */
  bool is_empty() const { return element_count == 0; }
  /**
   * @brief Reports whether the vector is empty.
   * @return True iff size() == 0.
   * @details Container-requirement spelling of is_empty().
   */
  bool empty() const { return element_count == 0; }

  /**
   * @brief Accesses the last element.
   * @return Mutable reference to the element at index size() - 1.
   */
  T &back() {
    check_alive();
    check_bound();
    assert(element_count > 0);
    return elements[element_count - 1];
  }
  /**
   * @brief Accesses the last element (const).
   * @return Const reference to the element at index size() - 1.
   */
  const T &back() const {
    check_alive();
    check_bound();
    assert(element_count > 0);
    return elements[element_count - 1];
  }

  /**
   * @brief Resets the vector to empty without destroying elements.
   * @details No check_bound() here on purpose: clear() only resets
   * element_count (it neither frees nor touches elements), so it is a defined
   * no-op on an unbound vector. MeshState::clear() relies on this to reset
   * members that may never have been bound.
   */
  void clear() {
    check_alive();
    element_count = 0;
  }

  /**
   * @brief Returns a pointer to the backing storage.
   * @return Mutable pointer to the first element, or nullptr if unbound or
   * moved-from.
   * @details No check_bound() on data()/begin()/end() on purpose: an unbound (or
   * moved-from) vector is elements==nullptr with element_count==0, a
   * well-defined EMPTY range that callers rely on (std::span(vertices.data(),
   * size()), std::sort(data(), data()+count) on a size-0 vector). The
   * use-after-free guard (check_alive) still applies.
   */
  T *data() {
    check_alive();
    return elements;
  }
  /**
   * @brief Returns a pointer to the backing storage (const).
   * @return Const pointer to the first element, or nullptr if unbound or
   * moved-from.
   */
  const T *data() const {
    check_alive();
    return elements;
  }

  /**
   * @brief Returns an iterator to the first element.
   * @return Mutable pointer to the first element.
   */
  T *begin() {
    check_alive();
    return elements;
  }
  /**
   * @brief Returns an iterator past the last element.
   * @return Mutable pointer one past the last element.
   */
  T *end() {
    check_alive();
    // Guard nullptr + 0 (formal UB) on an unbound/empty vector.
    return elements ? elements + element_count : nullptr;
  }
  /**
   * @brief Returns a const iterator to the first element.
   * @return Const pointer to the first element.
   */
  const T *begin() const {
    check_alive();
    return elements;
  }
  /**
   * @brief Returns a const iterator past the last element.
   * @return Const pointer one past the last element.
   */
  const T *end() const {
    check_alive();
    // Guard nullptr + 0 (formal UB) on an unbound/empty vector.
    return elements ? elements + element_count : nullptr;
  }
};

// ============================================================================
// 3. Non-Owning Span (Explicit Borrow)
// ============================================================================

/**
 * @brief A read-only, non-owning view into arena-allocated data.
 * @tparam T Element type viewed by the span.
 * @details Makes owned (ArenaVector) vs borrowed data visible at the type level.
 *
 * LIFETIME CONTRACT: a span snapshots its source vector's elements pointer at
 * construction. In debug builds two independent stamps fault on a stale span: the
 * arena generation catches an arena RESET, and the source vector's per-vector
 * rebind counter catches a bind()-driven RE-GROW (a grow rebinds elements
 * WITHOUT bumping the arena generation). A MOVE of the source vector is not
 * tracked: the span keeps its snapshotted elements (runtime-safe) but its debug
 * stamps reference the moved-from husk, so re-take the span after growing or
 * moving its source. Outliving the source VECTOR OBJECT (not its arena block —
 * a stack-local vector going out of scope) is worse than untracked: in debug the
 * staleness check reads source_vec, so the span must not outlive it.
 */
template <typename T> class ArenaSpan {
  const T *elements;    /**< Snapshotted pointer to the borrowed data. */
  size_t element_count; /**< Number of viewed elements. */
#ifndef NDEBUG
  ArenaBlockStamp stamp; /**< Arena state when the block was allocated. */
  const ArenaVector<T> *source_vec =
      nullptr; /**< Source vector for re-grow check. */
  uint32_t source_rebind_generation =
      0; /**< Vector rebind counter at construction. */

  /**
   * @brief Debug-only stale-span check against arena and vector stamps.
   * @details Asserts on an arena reset or a source-vector re-grow.
   */
  void check_alive() const {
    const size_t bytes = element_count * sizeof(T);
    assert(!stamp.arena_reset() && "ArenaSpan use-after-free!");
    assert((!source_vec ||
            source_vec->rebind_generation == source_rebind_generation) &&
           "ArenaSpan source vector re-grown out from under span!");
    assert(!stamp.block_uncovered(elements, bytes) &&
           "ArenaSpan use-after-free (arena rewound below span)!");
    assert(!stamp.block_reissued(elements, bytes) &&
           "ArenaSpan use-after-free (borrowed block reclaimed by a rewind and "
           "reissued)!");
  }
#else
  /**
   * @brief No-op stale-span check in release builds.
   */
  void check_alive() const {}
#endif

public:
  // Read-only view: the mutable spellings alias the const ones.
  using value_type = T;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using reference = const T &;
  using const_reference = const T &;
  using pointer = const T *;
  using const_pointer = const T *;
  using iterator = const T *;
  using const_iterator = const T *;

  /**
   * @brief Default-constructs an empty span.
   */
  ArenaSpan() : elements(nullptr), element_count(0) {}

  /**
   * @brief Constructs a span borrowing from an ArenaVector (explicit borrow).
   * @param source Vector to borrow data and (in debug) lifetime stamps from.
   * @details In debug builds the span inherits the vector's arena-generation
   * stamp, so accessing it after the arena is reset/compacted trips the same
   * use-after-free check the vector itself has.
   */
  explicit ArenaSpan(const ArenaVector<T> &source)
      : elements(source.data()), element_count(source.size())
#ifndef NDEBUG
        ,
        stamp(source.stamp), source_vec(&source),
        source_rebind_generation(source.rebind_generation)
#endif
  {
  }

  /**
   * @brief Deleted constructor from a temporary ArenaVector.
   * @details Borrowing from a temporary would leave the data pointer and (in
   * debug) source_vec dangling the moment the temporary dies. Forbid it.
   */
  explicit ArenaSpan(const ArenaVector<T> &&) = delete;

  /**
   * @brief Copy and copy-assignment duplicate the borrow verbatim.
   * @details A span copy carries the same data pointer, size, and (in debug) the
   * source's lifetime stamps, so the copy trips the same staleness check as the
   * original. Only construction from a temporary ArenaVector (above) is forbidden.
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
    assert(i < element_count);
    return elements[i];
  }
  /**
   * @brief Returns the number of viewed elements.
   * @return Element count.
   */
  size_t size() const {
    check_alive();
    return element_count;
  }
  /**
   * @brief Reports whether the span is empty.
   * @return True iff size() == 0.
   */
  bool is_empty() const {
    check_alive();
    return element_count == 0;
  }
  /**
   * @brief Reports whether the span is empty.
   * @return True iff size() == 0.
   * @details Container-requirement spelling of is_empty().
   */
  bool empty() const {
    check_alive();
    return element_count == 0;
  }
  /**
   * @brief Returns a pointer to the borrowed storage.
   * @return Const pointer to the first element.
   */
  const T *data() const {
    check_alive();
    return elements;
  }
  /**
   * @brief Returns a const iterator to the first element.
   * @return Const pointer to the first element.
   */
  const T *begin() const {
    check_alive();
    return elements;
  }
  /**
   * @brief Returns a const iterator past the last element.
   * @return Const pointer one past the last element.
   */
  const T *end() const {
    check_alive();
    // Guard nullptr + 0 (formal UB) on a default-constructed/empty span.
    return elements ? elements + element_count : nullptr;
  }
};

extern Arena scratch_arena_a;
extern Arena scratch_arena_b;

extern Arena persistent_arena;

/**
 * @brief Self-registering callback run before arena storage is handed out
 *        again.
 * @details A global that caches a pointer into an arena declares one static
 * instance next to itself and drops the pointer from the callback, instead of
 * the allocator naming every such owner. The registry head is
 * constant-initialized, so registration during static init is order-independent;
 * the list is intrusive, so it needs no storage of its own.
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
 * @brief Repartitions the global arena budget across the three arenas.
 * @param persistent Bytes to assign to the persistent arena.
 * @param scratch_a Bytes to assign to scratch arena A.
 * @param scratch_b Bytes to assign to scratch arena B.
 */
FLASHMEM void configure_arenas(size_t persistent, size_t scratch_a,
                               size_t scratch_b);
/**
 * @brief Restores the default arena partition.
 */
FLASHMEM void configure_arenas_default();

/**
 * @brief Re-partitions the arenas mid-run WITHOUT disturbing persistent content.
 * @param persistent New persistent capacity; must be >= its current live offset.
 * @param scratch_a New scratch-A capacity.
 * @param scratch_b New scratch-B capacity.
 * @details Unlike configure_arenas(), the persistent arena keeps its base
 * (block start), offset, live content, and generation -- only its capacity
 * boundary moves -- so the long-lived carousel slots + palette bank below its
 * offset survive. The scratch arenas hold nothing across the call point
 * (transient, reset every frame), so they rebind to fresh bases. Callers MUST
 * invoke this only when both scratch arenas are empty; a per-shape split at
 * spawn (persistent at its ~baseline, scratch idle) satisfies this.
 */
FLASHMEM void resplit_arenas(size_t persistent, size_t scratch_a,
                             size_t scratch_b);

// ============================================================================
// 4. ScratchScope — RAII Arena Offset Guard
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
private:
  Arena &arena;        /**< Arena whose offset is saved and restored. */
  size_t saved_offset; /**< Offset captured at construction. */

public:
  /**
   * @brief Constructs the scope, saving the arena's current offset.
   * @param a Arena to guard.
   */
  explicit ScratchScope(Arena &a) : arena(a), saved_offset(a.get_offset()) {}
  /**
   * @brief Destroys the scope, rewinding the arena to the saved offset.
   * @details Enforces LIFO scope discipline before rewinding. Stack-nested scopes
   * on the SAME arena are always safe: an inner scope rewinds to exactly where the
   * outer's allocations end, so nesting never clobbers a live allocation. The one
   * way to break it is non-LIFO teardown: an outer rewind or reset() ran while
   * this scope was live, dropping the offset BELOW saved_offset. Trap it instead
   * of letting set_offset() resurrect freed bytes.
   */
  ~ScratchScope() {
    HS_CHECK(arena.get_offset() >= saved_offset,
             "ScratchScope: non-LIFO teardown — arena at %lu, saved %lu",
             static_cast<unsigned long>(arena.get_offset()),
             static_cast<unsigned long>(saved_offset));
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
// 5. RAII Arena Evacuator
// ============================================================================

/**
 * @brief Concept requiring a static clone(const T&, T&, Arena&) method.
 * @tparam T Type that must provide static void clone(const T&, T&, Arena&).
 * @note Cloneable only constrains the clone hook. `Persist<T>` needs more (T
 *       default-initializable and assignable, because ~Persist does
 *       `target = T()` before restoring) and enforces those extra requirements
 *       with its own static_asserts rather than widening this concept.
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
  T &target;                        /**< Object being evacuated and restored. */
  Arena &persistent;                /**< Arena the object is restored into. */
  size_t persistent_offset_at_ctor; /**< persistent offset at construction; the
                                          dtor traps unless the caller rewound
                                          below this watermark. */

  // Declaration order matters! scratch must be declared BEFORE backup
  // so that backup is destroyed before the scratch arena is rolled back.
  ScratchScope scratch; /**< Scratch scope holding the backup's storage. */
  T backup;             /**< Cloned backup of the target in scratch memory. */

  // ~Persist does `target = T()` before re-cloning the backup, so T needs more
  // than Cloneable. Assert the extra requirements here so a Cloneable-but-not-
  // default-constructible/assignable T fails with a clear message at instantiation.
  static_assert(
      std::default_initializable<T>,
      "Persist<T>: ~Persist resets target = T() before restoring, so "
      "T must be default-initializable (Cloneable does not imply this).");
  static_assert(std::assignable_from<T &, T>,
                "Persist<T>: ~Persist assigns target = T(), so T must be "
                "assignable from a T rvalue (Cloneable does not imply this).");

public:
  /**
   * @brief Evacuates the target into the scratch arena.
   * @param subject Object to back up and later restore.
   * @param scratch_arena Scratch arena to hold the backup.
   * @param restore_arena Persistent arena the target is restored into.
   */
  HS_COLD_MEMBER Persist(T &subject, Arena &scratch_arena, Arena &restore_arena)
      : target(subject), persistent(restore_arena),
        persistent_offset_at_ctor(restore_arena.get_offset()),
        scratch(scratch_arena) {
    HS_CHECK(&scratch_arena != &restore_arena,
             "Persist: scratch and persistent must be distinct arenas — the "
             "dtor's watermark restore assumes the backup lives in a different "
             "arena than the one it restores into");
    T::clone(target, backup, scratch.get_arena());
  }

  /**
   * @brief Restores the target by cloning the backup into the persistent arena.
   * @details The restore clones into `persistent` at its *current* offset, so
   * it only reconstructs the object usefully if the caller rewound the
   * persistent arena during the scope (the canonical
   * `persistent_arena.reset()`). Without that reset the clone appends a second
   * copy and grows the arena. The post-restore `<=` check traps a
   * forgot-to-reset restore, which pushes the offset past the construction
   * watermark. Callers may legitimately STACK several Persists over one
   * `reset()`, so the check bounds the aggregate, not each individual restore —
   * a backstop, not a proof.
   */
  HS_COLD_MEMBER ~Persist() {
    target = T();
    T::clone(backup, target, persistent);
    HS_CHECK(persistent.get_offset() <= persistent_offset_at_ctor,
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

// ============================================================================
// 6. generate() — Scratch-Scoped Generation Wrapper
// ============================================================================

namespace hs {
namespace detail {
/**
 * @brief Returns the shared generate() nesting-depth counter.
 * @return Reference to the process-wide recursion depth (starts at 0).
 * @details The counter is shared across ALL generate() instantiations, so it must
 *   live outside the template: a function-local static inside generate() would give
 *   each GenerateFn type its own counter, so a nested call (a different lambda type)
 *   would see depth==0 and reset the arena out from under its caller.
 */
inline int &generate_depth() {
  static int depth = 0;
  return depth;
}
} // namespace detail

constexpr int MAX_GENERATE_DEPTH = 16;

/**
 * @brief Universal generation wrapper for procedural geometry.
 * @tparam GenerateFn The generation function or callable type.
 * @tparam Args       Types of the extra arguments forwarded to fn.
 * @param target The arena for persistent output (e.g., persistent_arena). Must
 *   not be either engine scratch arena (`scratch_arena_a`/`scratch_arena_b`):
 *   the depth-0 reset and `ScratchScope` rewind would destroy the output.
 * @param fn     The generation function or callable.
 * @param args   Extra arguments forwarded to fn.
 * @return Whatever fn returns.
 * @details Resets and scopes both scratch arenas, then invokes
 *   fn(target, scratch_a, scratch_b, args...). All procedural geometry creation
 *   goes through this wrapper for a deterministic arena lifecycle. The fn signature
 *   must be: ReturnType fn(Arena& target, Arena& scratch_a, Arena& scratch_b, Args...)
 *
 *   Reentrant: a generator callback may itself call generate(). The full arena
 *   reset happens only at the outermost call; a nested call sub-scopes off the
 *   caller's live frame (via ScratchScope), so it sees only the scratch headroom
 *   above the outer allocations rather than clobbering them.
 */
template <typename GenerateFn, typename... Args>
HS_COLD_MEMBER auto generate(Arena &target, GenerateFn &&fn, Args &&...args) {
  HS_CHECK(&target != &scratch_arena_a && &target != &scratch_arena_b,
           "generate: target must not alias an engine scratch arena");
  int &depth = detail::generate_depth();
  if (depth == 0) {
    scratch_arena_a.reset();
    scratch_arena_b.reset();
  }
  ScratchScope a_guard(scratch_arena_a);
  ScratchScope b_guard(scratch_arena_b);
  ++depth;
  HS_CHECK(depth <= MAX_GENERATE_DEPTH, "generate: recursion too deep");
  /**
   * @brief RAII guard that decrements the nesting depth on scope exit.
   */
  struct DepthGuard {
    int &d; /**< Reference to the shared nesting-depth counter. */
    /**
     * @brief Decrements the referenced depth counter on destruction.
     */
    ~DepthGuard() { --d; }
  } depth_guard{depth};
  return fn(target, scratch_arena_a, scratch_arena_b,
            std::forward<Args>(args)...);
}

} // namespace hs
