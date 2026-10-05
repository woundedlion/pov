/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/memory.h.

// ============================================================================
// RAII Arena Evacuator
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
 *     reset_persistent_arena();
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
   * `reset_persistent_arena()`). Without that reset the clone appends a second
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
