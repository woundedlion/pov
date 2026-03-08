#include "memory.h"
#include "mesh.h"

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

// Static arena backing stores — plain static = RAM1 (DTCM) on Teensy
static uint8_t scratch_buf_a[SCRATCH_ARENA_A_SIZE];
static uint8_t scratch_buf_b[SCRATCH_ARENA_B_SIZE];
static uint8_t persistent_buf[PERSISTENT_ARENA_SIZE];

Arena scratch_arena_a(scratch_buf_a, SCRATCH_ARENA_A_SIZE);
Arena scratch_arena_b(scratch_buf_b, SCRATCH_ARENA_B_SIZE);
Arena persistent_arena(persistent_buf, PERSISTENT_ARENA_SIZE);

MeshState *PersistentTracker::tracked_meshes[256];
size_t PersistentTracker::num_tracked = 0;

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
