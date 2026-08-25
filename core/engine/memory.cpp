/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#include "engine/memory.h"

namespace hs {
[[noreturn]] HS_COLD void function_ref_empty_call() {
  check_fail(__FILE__, __LINE__, "thunk != empty_thunk",
             "empty FunctionRef called");
}
#ifndef ARDUINO
[[noreturn]] HS_COLD void inplace_function_empty_call() {
  check_fail(__FILE__, __LINE__, "vtable != empty",
             "empty hs::inplace_function called");
}
#endif
} // namespace hs
#ifdef ARDUINO
#include <exception>
#endif

// Stubs to prevent the linker pulling in the C++ demangler (~15KB) on the device
// (chain: std::function -> __cxa_throw -> __verbose_terminate_handler). Only the
// Arduino/Teensy build needs them: on host the toolchain's real handlers are
// correct and defining these would shadow them. A pure-virtual call or terminate
// is a fatal invariant violation — flush the log and trap (fail-fast), never spin
// in while(1), which freezes the display with no breadcrumb.
#ifdef ARDUINO
/**
 * @brief Fail-fast handler for a pure-virtual call on the device.
 * @details Flushes the log then traps; never spins in while(1). Arduino/Teensy-
 * only override so the host toolchain's real handlers stay in effect on
 * native/WASM builds.
 */
extern "C" void __cxa_pure_virtual() {
  hs::flush_log();
  __builtin_trap();
}
#if TEENSYDUINO < 162
namespace __gnu_cxx {
/**
 * @brief Fail-fast terminate handler for the device.
 * @details Flushes the log then traps instead of pulling in the C++ demangler
 * (~15KB). A terminate is a fatal invariant violation, so it must leave a
 * breadcrumb and stop, not spin. On TD >= 1.62 the core defines its own strong
 * __verbose_terminate_handler, so we install the fail-fast handler at runtime
 * (below) instead of redefining the symbol.
 */
void __verbose_terminate_handler() {
  hs::flush_log();
  __builtin_trap();
}
} // namespace __gnu_cxx
#else
// TD >= 1.62's core (cores/teensy4/main.cpp) ships its own strong
// __gnu_cxx::__verbose_terminate_handler (a while(1) WFI spin), so redefining it
// here would be a multiple-definition link error. That core handler already keeps
// the ~15KB demangler out; we only need to replace its spin-and-freeze behavior
// with the fail-fast flush+trap. Install it at runtime via a global constructor
// (runs before setup()); std::terminate then dispatches to ours.
namespace {
const std::terminate_handler s_fail_fast_terminate = std::set_terminate([] {
  hs::flush_log();
  __builtin_trap();
});
} // namespace
#endif
#endif

/**
 * @brief Single contiguous memory block that all arenas partition.
 * @details alignas keeps the block's base maximally aligned, so each arena's
 * first allocation needs no leading padding; configure_arenas() likewise aligns
 * the inter-arena boundaries. Carries NO DMAMEM: the arena is hot render memory,
 * so on the device it lands in .bss/DTCM, the Cortex-M7's fastest zero-wait RAM.
 * The framebuffers (buffer_a/buffer_b, ~243 KiB each) and the timeline events in
 * static_storage.cpp carry DMAMEM to place them in OCRAM instead — neither is a
 * DMA source or target; DTCM simply has no room for them.
 */
alignas(std::max_align_t) static uint8_t global_arena_block[GLOBAL_ARENA_SIZE];

/**
 * @brief Persistent arena: storage that lives for the whole program run.
 * @details Growing the persistent boundary goes through resplit_arenas(), which
 * re-bases both scratch arenas out of the way before moving the capacity.
 */
Arena persistent_arena(global_arena_block, DEFAULT_PERSISTENT_SIZE);
/** @brief First scratch arena: transient per-frame/per-effect storage. */
Arena scratch_arena_a(global_arena_block + DEFAULT_PERSISTENT_SIZE,
                      DEFAULT_SCRATCH_A_SIZE);
/** @brief Second scratch arena: transient per-frame/per-effect storage. */
Arena scratch_arena_b(global_arena_block + DEFAULT_PERSISTENT_SIZE +
                          DEFAULT_SCRATCH_A_SIZE,
                      DEFAULT_SCRATCH_B_SIZE);

namespace {
size_t abandoned_bytes_total = 0;
size_t abandon_event_count = 0;
} // namespace

void note_arena_vector_abandon(size_t bytes) {
  abandoned_bytes_total += bytes;
  abandon_event_count++;
}

FLASHMEM void log_arena_vector_grow(size_t bytes, size_t old_capacity,
                                    size_t new_capacity) {
  note_arena_vector_abandon(bytes);
  hs::log("ArenaVector grow abandons %lu bytes (cap %lu -> %lu)",
          static_cast<unsigned long>(bytes),
          static_cast<unsigned long>(old_capacity),
          static_cast<unsigned long>(new_capacity));
}

FLASHMEM size_t arena_vector_abandoned_bytes() { return abandoned_bytes_total; }

FLASHMEM size_t arena_vector_abandon_count() { return abandon_event_count; }

[[noreturn]] HS_COLD void arena_oom_trap(const void *buffer, size_t size,
                                         size_t offset, size_t padding,
                                         size_t capacity) {
  hs::log("[OOM] Arena @%08lx: req %lu, offset %lu, pad %lu / cap %lu "
          "(move-assign dropped %lu B in %lu blocks, all arenas since "
          "boot, reclaims not subtracted)",
          static_cast<unsigned long>(reinterpret_cast<uintptr_t>(buffer)),
          static_cast<unsigned long>(size), static_cast<unsigned long>(offset),
          static_cast<unsigned long>(padding),
          static_cast<unsigned long>(capacity),
          static_cast<unsigned long>(arena_vector_abandoned_bytes()),
          static_cast<unsigned long>(arena_vector_abandon_count()));
  hs::check_fail(__FILE__, __LINE__, "false", "Arena::allocate: out of memory");
}

namespace {
/** @brief Offsets of the two scratch arena bases within global_arena_block. */
struct ScratchBases {
  size_t a; /**< Base offset of scratch arena A. */
  size_t b; /**< Base offset of scratch arena B. */
};

/**
 * @brief Resolves a requested split into aligned scratch base offsets.
 * @param who Caller name for the out-of-budget log line.
 * @param persistent Bytes requested for the persistent arena.
 * @param scratch_a Bytes requested for scratch arena A.
 * @param scratch_b Bytes requested for scratch arena B.
 * @return Base offsets of the two scratch arenas.
 * @details Each input is bounded first so the align_up()/sum arithmetic cannot
 * wrap size_t. Each inter-arena boundary is aligned up to max_align_t (the real
 * callers pass 1 KiB multiples, so these rounds are no-ops); the budget check
 * uses the aligned end so rounding cannot silently overrun. An over-subscribed
 * partition is a sizing/config bug, not recoverable, so it traps rather than
 * silently scaling down; the log preserves the numbers before the trap.
 */
HS_COLD ScratchBases split_bases(const char *who, size_t persistent,
                                 size_t scratch_a, size_t scratch_b) {
  HS_CHECK(persistent <= GLOBAL_ARENA_SIZE && scratch_a <= GLOBAL_ARENA_SIZE &&
               scratch_b <= GLOBAL_ARENA_SIZE,
           "split_bases: %s asked %lu/%lu/%lu B (persistent/A/B) of a %lu B "
           "block",
           who, static_cast<unsigned long>(persistent),
           static_cast<unsigned long>(scratch_a),
           static_cast<unsigned long>(scratch_b),
           static_cast<unsigned long>(GLOBAL_ARENA_SIZE));
  constexpr size_t A = alignof(std::max_align_t);
  auto align_up = [](size_t n) { return (n + (A - 1)) & ~(A - 1); };
  size_t a_base = align_up(persistent);
  size_t b_base = align_up(a_base + scratch_a);
  size_t total = b_base + scratch_b;
  if (total > GLOBAL_ARENA_SIZE) {
    hs::log("[OOM] %s: requested %lu > available %lu", who,
            static_cast<unsigned long>(total),
            static_cast<unsigned long>(GLOBAL_ARENA_SIZE));
    HS_CHECK(false);
  }
  return {a_base, b_base};
}
} // namespace

/**
 * @brief Re-partitions the single global block into persistent plus two scratch
 * arenas of the requested byte sizes.
 * @details Called once at init() so an effect can tune the split to the device
 * budget; split_bases() aligns the boundaries and enforces the budget.
 *
 * Runs the ArenaResetHook list first: this is one of the paths that hands the
 * storage under a persistent-arena-resident global out again. Owners re-arm
 * from init(), which every swap path runs after this.
 */
FLASHMEM void configure_arenas(size_t persistent, size_t scratch_a,
                               size_t scratch_b) {
  ArenaResetHook::run_all();
  const ScratchBases bases =
      split_bases("configure_arenas", persistent, scratch_a, scratch_b);
  persistent_arena.rebind(global_arena_block, persistent);
  scratch_arena_a.rebind(global_arena_block + bases.a, scratch_a);
  scratch_arena_b.rebind(global_arena_block + bases.b, scratch_b);
}

/**
 * @brief Partitions the global block using the compiled-in DEFAULT_* sizes.
 * @details Convenience wrapper over configure_arenas() with the default split.
 */
FLASHMEM void configure_arenas_default() {
  configure_arenas(DEFAULT_PERSISTENT_SIZE, DEFAULT_SCRATCH_A_SIZE,
                   DEFAULT_SCRATCH_B_SIZE);
}

FLASHMEM void resplit_arenas(size_t persistent, size_t scratch_a,
                             size_t scratch_b) {
  // A ScratchScope saved at offset 0 restores to 0 either way, so live scratch
  // content would be silently rebased onto the new split, undetected.
  HS_CHECK(scratch_arena_a.get_offset() == 0 &&
               scratch_arena_b.get_offset() == 0,
           "resplit_arenas: both scratch arenas must be empty");
  const ScratchBases bases =
      split_bases("resplit_arenas", persistent, scratch_a, scratch_b);
  // Persistent keeps its base/offset/content/generation -- only the boundary
  // moves. rebind_capacity traps if the new budget is below the live offset, so
  // the carousel + palette bank are never stranded. Its high-water is rebased to
  // the live offset (the per-shape baseline) so each shape's persistent peak is
  // measured against its own split, not a predecessor's residue.
  persistent_arena.rebind_capacity(persistent);
  persistent_arena.reset_high_water_mark();
  // The scratch arenas are empty at the call point; rebind them onto their new
  // bases (a fresh generation is harmless -- nothing is bound across the split).
  scratch_arena_a.rebind(global_arena_block + bases.a, scratch_a);
  scratch_arena_b.rebind(global_arena_block + bases.b, scratch_b);
}
