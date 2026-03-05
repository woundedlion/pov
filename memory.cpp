#include "memory.h"
#include "mesh.h"

// Allocations
Arena scratch_arena_a(SCRATCH_ARENA_A_SIZE);
Arena scratch_arena_b(SCRATCH_ARENA_B_SIZE);
Arena persistent_arena(PERSISTENT_ARENA_SIZE);

ArenaVector<MeshState *> PersistentTracker::tracked_meshes;

void PersistentTracker::auto_compact(Arena &safe_scratch) {
  safe_scratch.reset();

  // 1. Evacuate tracked meshes to the safe scratch arena
  MeshState *temp_rafts = static_cast<MeshState *>(
      safe_scratch.allocate(tracked_meshes.size() * sizeof(MeshState)));

  for (size_t i = 0; i < tracked_meshes.size(); ++i) {
    new (&temp_rafts[i]) MeshState();
    temp_rafts[i].deep_copy(safe_scratch, *(tracked_meshes[i]));
  }

  // 2. Wipe the fragmented persistent arena completely
  persistent_arena.reset();

  // 3. Bring everything back, perfectly packed side-by-side!
  for (size_t i = 0; i < tracked_meshes.size(); ++i) {
    tracked_meshes[i]->clear();
    tracked_meshes[i]->deep_copy(persistent_arena, temp_rafts[i]);
  }
  printf("Auto-compacted persistent arena seamlessly.\n");
}

void MemoryCtx::update_persistent(MeshState &target, const PolyMesh &new_data) {
  // Check if we are about to OOM the persistent arena
  size_t required = (new_data.vertices.capacity() * sizeof(Vector)) * 2;
  if (persistent_arena.get_capacity() - persistent_arena.get_offset() <
      required) {
    // Trigger invisible defragmentation using our currently INACTIVE scratch
    // buffer!
    PersistentTracker::auto_compact(*scratch_back);
  }

  // Safely write the new shape into the persistent arena
  target.clear();
  MeshOps::compile(new_data, target, persistent_arena);
}
