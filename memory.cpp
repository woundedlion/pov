#include "memory.h"
#include "animation.h"

// Stubs to prevent the linker pulling in the C++ demangler (~15KB).
// Chain: std::function -> __cxa_throw -> __verbose_terminate_handler ->
// d_print_*
extern "C" void __cxa_pure_virtual() {
  while (1)
    ;
}
namespace __gnu_cxx {
void __verbose_terminate_handler() {
  while (1)
    ;
}
} // namespace __gnu_cxx

// Single contiguous memory block — all arenas partition this
static uint8_t global_arena_block[GLOBAL_ARENA_SIZE];

Arena persistent_arena(global_arena_block, DEFAULT_PERSISTENT_SIZE);
Arena scratch_arena_a(global_arena_block + DEFAULT_PERSISTENT_SIZE,
                      DEFAULT_SCRATCH_A_SIZE);
Arena scratch_arena_b(global_arena_block + DEFAULT_PERSISTENT_SIZE +
                          DEFAULT_SCRATCH_A_SIZE,
                      DEFAULT_SCRATCH_B_SIZE);

void configure_arenas(size_t persistent, size_t scratch_a, size_t scratch_b) {
  size_t total = persistent + scratch_a + scratch_b;
  if (total > GLOBAL_ARENA_SIZE) {
    hs::log("[WARN] configure_arenas: requested %zu, available %zu — scaling down",
            total, GLOBAL_ARENA_SIZE);
    float scale = static_cast<float>(GLOBAL_ARENA_SIZE) / total;
    persistent = static_cast<size_t>(persistent * scale);
    scratch_a = static_cast<size_t>(scratch_a * scale);
    scratch_b = GLOBAL_ARENA_SIZE - persistent - scratch_a;
  }
  persistent_arena.rebind(global_arena_block, persistent);
  scratch_arena_a.rebind(global_arena_block + persistent, scratch_a);
  scratch_arena_b.rebind(global_arena_block + persistent + scratch_a,
                         scratch_b);
}

void configure_arenas_default() {
  configure_arenas(DEFAULT_PERSISTENT_SIZE, DEFAULT_SCRATCH_A_SIZE,
                   DEFAULT_SCRATCH_B_SIZE);
}

DMAMEM TimelineEvent global_timeline_events[TIMELINE_MAX_EVENTS];
DMAMEM Pixel Effect::buffer_a[MAX_W * MAX_H];
DMAMEM Pixel Effect::buffer_b[MAX_W * MAX_H];
