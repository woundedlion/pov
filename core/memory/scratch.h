/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/memory.h.

// ============================================================================
// ScratchScope — RAII Arena Offset Guard
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
#ifndef NDEBUG
  uint32_t saved_generation;
#endif

public:
  /**
   * @brief Constructs the scope, saving the arena's current offset.
   * @param a Arena to guard.
   */
  explicit ScratchScope(Arena &a) : arena(a), saved_offset(a.get_offset()) {
#ifndef NDEBUG
    saved_generation = a.get_generation();
#endif
  }
  /**
   * @brief Destroys the scope, rewinding the arena to the saved offset.
   * @details Reports non-LIFO teardown before set_offset rejects the rewind.
   */
  ~ScratchScope() {
#ifndef NDEBUG
    HS_CHECK(arena.get_generation() == saved_generation,
             "ScratchScope: arena reset during scope lifetime");
#endif
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
