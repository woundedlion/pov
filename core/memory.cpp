#include "memory.h"
#include "animation.h"

// Stubs to prevent the linker pulling in the C++ demangler (~15KB) on the
// device. Chain: std::function -> __cxa_throw -> __verbose_terminate_handler ->
// d_print_*. Only the Arduino/Teensy build needs them: on host (native tests,
// WASM) the toolchain's real handlers are correct, and defining these would
// shadow them. A pure-virtual call or a terminate is a fatal invariant
// violation — flush the log and trap (fail-fast), never spin in while(1), which
// freezes the display with no breadcrumb and is the hardest failure to diagnose.
#ifdef ARDUINO
/**
 * @brief Fail-fast handler for a pure-virtual call on the device.
 * @details Flushes the log then traps; never spins in while(1), which would
 * freeze the display with no breadcrumb. Arduino/Teensy-only override so the
 * host toolchain's real handlers stay in effect on native/WASM builds.
 */
extern "C" void __cxa_pure_virtual() {
  hs::flush_log();
  __builtin_trap();
}
namespace __gnu_cxx {
/**
 * @brief Fail-fast terminate handler for the device.
 * @details Flushes the log then traps instead of pulling in the C++ demangler
 * (~15KB). A terminate is a fatal invariant violation, so it must leave a
 * breadcrumb and stop, not spin.
 */
void __verbose_terminate_handler() {
  hs::flush_log();
  __builtin_trap();
}
} // namespace __gnu_cxx
#endif

/**
 * @brief Single contiguous memory block that all arenas partition.
 * @details alignas keeps the block's base (and thus persistent_arena's base)
 * maximally aligned, so Arena::allocate's first allocation in each arena needs
 * no leading padding; configure_arenas() likewise aligns the inter-arena
 * boundaries.
 *
 * Placement: this block carries NO DMAMEM qualifier, which is deliberate. The
 * arena is the engine's hot per-frame render memory (every effect's persistent
 * and scratch allocations live inside it) and is never a DMA target, so on the
 * device it lands in the default .bss — DTCM, the Cortex-M7's zero-wait tightly-
 * coupled RAM, the fastest memory available. This is the single largest RAM
 * decision in the engine (GLOBAL_ARENA_SIZE ≈ 335 KB). The buffers below that
 * the LED-output DMA must reach (buffer_a/buffer_b, the timeline events) instead
 * carry DMAMEM precisely because the DMA engine cannot access DTCM; the arena,
 * touched only by the CPU, has the opposite requirement and stays in DTCM.
 */
alignas(std::max_align_t) static uint8_t global_arena_block[GLOBAL_ARENA_SIZE];

/** @brief Persistent arena: storage that lives for the whole program run. */
Arena persistent_arena(global_arena_block, DEFAULT_PERSISTENT_SIZE);
/** @brief First scratch arena: transient per-frame/per-effect storage. */
Arena scratch_arena_a(global_arena_block + DEFAULT_PERSISTENT_SIZE,
                      DEFAULT_SCRATCH_A_SIZE);
/** @brief Second scratch arena: transient per-frame/per-effect storage. */
Arena scratch_arena_b(global_arena_block + DEFAULT_PERSISTENT_SIZE +
                          DEFAULT_SCRATCH_A_SIZE,
                      DEFAULT_SCRATCH_B_SIZE);

/**
 * @brief Re-partitions the single global block into persistent plus two scratch
 * arenas of the requested byte sizes.
 * @details Called once at init() so an effect can tune the split to the device
 * budget. Each inter-arena boundary is aligned up to max_align_t so every arena
 * base is maximally aligned (the block base already is, via alignas); for the
 * real callers all sizes are multiples of 1 KiB, so these rounds are no-ops and
 * the layout is byte-identical. The budget check uses the aligned end so
 * rounding cannot silently overrun. An over-subscribed partition is a
 * sizing/config bug, not a recoverable condition, so it traps at init() rather
 * than silently scaling down (which would relocate the corruption into the
 * rendered output); the log preserves the numbers for diagnosis before the trap.
 */
void configure_arenas(size_t persistent, size_t scratch_a, size_t scratch_b) {
  // Bound each input so the align_up()/sum arithmetic below cannot wrap size_t.
  HS_CHECK(persistent <= GLOBAL_ARENA_SIZE && scratch_a <= GLOBAL_ARENA_SIZE &&
           scratch_b <= GLOBAL_ARENA_SIZE);
  constexpr size_t A = alignof(std::max_align_t);
  auto align_up = [](size_t n) { return (n + (A - 1)) & ~(A - 1); };
  size_t a_base = align_up(persistent);
  size_t b_base = align_up(a_base + scratch_a);
  size_t total = b_base + scratch_b;
  if (total > GLOBAL_ARENA_SIZE) {
    hs::log("[FATAL] configure_arenas: requested %zu > available %zu",
            total, GLOBAL_ARENA_SIZE);
    HS_CHECK(false);
  }
  persistent_arena.rebind(global_arena_block, persistent);
  scratch_arena_a.rebind(global_arena_block + a_base, scratch_a);
  scratch_arena_b.rebind(global_arena_block + b_base, scratch_b);
}

/**
 * @brief Partitions the global block using the compiled-in DEFAULT_* sizes.
 * @details Convenience wrapper over configure_arenas() with the default split.
 */
void configure_arenas_default() {
  configure_arenas(DEFAULT_PERSISTENT_SIZE, DEFAULT_SCRATCH_A_SIZE,
                   DEFAULT_SCRATCH_B_SIZE);
}

// The large statically-allocated buffers below are defined here, not next to
// their declarations, on purpose: memory.cpp is the one TU compiled into every
// target, so co-locating them with the arena block keeps every DMAMEM/large-
// static placement decision in one file the linker map points at — at the cost
// of this TU pulling in the animation and effect headers. Look here, not in
// animation.h / effect.h, for where the storage actually lands.

/** @brief Shared event array backing every Timeline instance. */
DMAMEM TimelineEvent global_timeline_events[TIMELINE_MAX_EVENTS];
/**
 * @brief Single live-Timeline guard.
 * @details The event array above is shared by every Timeline instance, so the
 * "only one alive" invariant is one global flag.
 */
bool global_timeline_live = false;
/**
 * @brief Shared singleton playhead cursor into global_timeline_events.
 * @details Free global so every Timeline instance reads/writes the same cursor
 * (see animation.h).
 */
int global_timeline_t = 0;
/**
 * @brief Shared singleton event count for global_timeline_events.
 * @details Free global so every Timeline instance reads/writes the same count
 * (see animation.h).
 */
int global_timeline_num_events = 0;
/** @brief Front pixel buffer for the double-buffered effect framebuffer. */
DMAMEM Pixel Effect::buffer_a[MAX_W * MAX_H];
/** @brief Back pixel buffer for the double-buffered effect framebuffer. */
DMAMEM Pixel Effect::buffer_b[MAX_W * MAX_H];
/** @brief Single-live-Effect guard for the shared buffer_a/buffer_b (see Effect). */
bool Effect::s_alive = false;
