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
extern "C" void __cxa_pure_virtual() {
  hs::flush_log();
  __builtin_trap();
}
namespace __gnu_cxx {
void __verbose_terminate_handler() {
  hs::flush_log();
  __builtin_trap();
}
} // namespace __gnu_cxx
#endif

// Single contiguous memory block — all arenas partition this. alignas keeps
// the block's base (and thus persistent_arena's base) maximally aligned, so
// Arena::allocate's first allocation in each arena needs no leading padding;
// configure_arenas() likewise aligns the inter-arena boundaries.
alignas(std::max_align_t) static uint8_t global_arena_block[GLOBAL_ARENA_SIZE];

Arena persistent_arena(global_arena_block, DEFAULT_PERSISTENT_SIZE);
Arena scratch_arena_a(global_arena_block + DEFAULT_PERSISTENT_SIZE,
                      DEFAULT_SCRATCH_A_SIZE);
Arena scratch_arena_b(global_arena_block + DEFAULT_PERSISTENT_SIZE +
                          DEFAULT_SCRATCH_A_SIZE,
                      DEFAULT_SCRATCH_B_SIZE);

// Re-partition the single global block into persistent + two scratch arenas of
// the requested byte sizes, called once at init() so an effect can tune the
// split to the device budget.
void configure_arenas(size_t persistent, size_t scratch_a, size_t scratch_b) {
  // An over-subscribed partition is a sizing/config bug, not a recoverable
  // condition. Silently scaling it down hands the effect a smaller layout than
  // it requested, which then over-runs later — relocating the corruption into
  // the rendered output. Trap at init() so the bad request is caught on the
  // bench; the log preserves the numbers for diagnosis before the trap.
  // Align each inter-arena boundary up to max_align_t so every arena base is
  // maximally aligned (the block base already is, via alignas). For the real
  // callers — all sizes are multiples of 1 KiB — these rounds are no-ops, so
  // the layout is byte-identical; they only matter for arbitrary sizes. The
  // budget check uses the ALIGNED end so rounding can't silently overrun.
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

// Partition the global block using the compiled-in DEFAULT_* sizes.
void configure_arenas_default() {
  configure_arenas(DEFAULT_PERSISTENT_SIZE, DEFAULT_SCRATCH_A_SIZE,
                   DEFAULT_SCRATCH_B_SIZE);
}

// These definitions live here, not next to their declarations, on purpose:
// memory.cpp is the one TU compiled into every target, and these are the
// program's large statically-allocated buffers. Co-locating them with the
// arena block keeps every DMAMEM/large-static placement decision in one file
// the linker map points at — at the cost of this TU pulling in the animation
// and effect headers. Look here, not in animation.h / effect.h, for where the
// storage actually lands.
DMAMEM TimelineEvent global_timeline_events[TIMELINE_MAX_EVENTS];
// Single live-Timeline guard: the event array above is shared by every Timeline
// instance, so the "only one alive" invariant is one global flag.
bool global_timeline_live = false;
// Shared singleton cursors into global_timeline_events (see animation.h):
// free globals so every Timeline instance reads/writes the same count.
int global_timeline_t = 0;
int global_timeline_num_events = 0;
DMAMEM Pixel Effect::buffer_a[MAX_W * MAX_H];
DMAMEM Pixel Effect::buffer_b[MAX_W * MAX_H];
