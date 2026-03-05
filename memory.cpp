#include "memory.h"
#include "mesh.h"

// Allocations
Arena scratch_arena_a(SCRATCH_ARENA_A_SIZE);
Arena scratch_arena_b(SCRATCH_ARENA_B_SIZE);
Arena persistent_arena(PERSISTENT_ARENA_SIZE);

MeshState *PersistentTracker::tracked_meshes[256];
size_t PersistentTracker::num_tracked = 0;

void PersistentTracker::auto_compact(Arena &scratch) {
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
    PersistentTracker::auto_compact(*scratch_back);
  }

  target.clear();
  MeshOps::compile(new_data, target, persistent_arena);
}
