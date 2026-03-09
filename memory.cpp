#include "memory.h"
#include "mesh.h"
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
Arena scratch_arena_a(global_arena_block + DEFAULT_PERSISTENT_SIZE, DEFAULT_SCRATCH_A_SIZE);
Arena scratch_arena_b(global_arena_block + DEFAULT_PERSISTENT_SIZE + DEFAULT_SCRATCH_A_SIZE, DEFAULT_SCRATCH_B_SIZE);

void configure_arenas(size_t persistent, size_t scratch_a, size_t scratch_b) {
  assert(persistent + scratch_a + scratch_b <= GLOBAL_ARENA_SIZE &&
         "Arena partition exceeds global block!");
  persistent_arena.rebind(global_arena_block, persistent);
  scratch_arena_a.rebind(global_arena_block + persistent, scratch_a);
  scratch_arena_b.rebind(global_arena_block + persistent + scratch_a, scratch_b);
}

void configure_arenas_default() {
  configure_arenas(DEFAULT_PERSISTENT_SIZE, DEFAULT_SCRATCH_A_SIZE, DEFAULT_SCRATCH_B_SIZE);
}

MeshState *PersistentTracker::tracked_meshes[256];
size_t PersistentTracker::num_tracked = 0;

DMAMEM TimelineEvent global_timeline_events[TIMELINE_MAX_EVENTS];
DMAMEM Pixel Effect::buffer_a[MAX_W * MAX_H];
DMAMEM Pixel Effect::buffer_b[MAX_W * MAX_H];

void PersistentTracker::auto_compact(Arena &scratch) {
  assert(!CompactionLock::is_locked() &&
         "Compaction while persistent pointers are held!");
  ScopedScratch _(scratch);

  // Evacuate
  MeshState *temp_rafts = static_cast<MeshState *>(
      scratch.allocate(num_tracked * sizeof(MeshState)));

  for (size_t i = 0; i < num_tracked; ++i) {
    new (&temp_rafts[i]) MeshState();
    MeshOps::clone(*(tracked_meshes[i]), temp_rafts[i], scratch);
  }

  // Clear
  persistent_arena.reset();

  // Restore
  for (size_t i = 0; i < num_tracked; ++i) {
    *(tracked_meshes[i]) = MeshState();
    MeshOps::clone(temp_rafts[i], *(tracked_meshes[i]), persistent_arena);
  }
}

void MemoryCtx::update_persistent(MeshState &target, const PolyMesh &new_data) {
  size_t required = (new_data.vertices.capacity() * sizeof(Vector)) * 2;
  if (persistent_arena.get_capacity() - persistent_arena.get_offset() <
      required) {
    PersistentTracker::auto_compact(get_scratch_back());
  }

  target.clear();
  MeshOps::compile(new_data, target, persistent_arena);
}
