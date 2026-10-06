/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Transformers death fixtures and guard cases.

/**
 * @brief Death case: a second TransformerPool::init_storage() must trap.
 * @details Transformer surface — a re-init would hand the pool a second block
 *          while spawned animations still hold Params references into the first,
 *          and would silently double-charge the persistent arena.
 */
inline void case_transformer_pool_init_storage_twice() {
  configure_arenas_default();
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);
  rt.init_storage(persistent_arena); // entities already set -> HS_CHECK
}

/** @brief Unpinned perpetual spawns trap even when the timeline is full. */
inline void case_transformer_unpinned_full() {
  configure_arenas_default();
  Timeline timeline;
  NoiseTransformer<1> transformer(timeline);
  transformer.init_storage(persistent_arena);
  float sink = 0.0f;
  for (int i = 0; i < Timeline::MAX_EVENTS; ++i)
    timeline.add(0, Animation::Transition(sink, 1.0f, 10, math::ease_linear));
  transformer.spawn(0);
}

/**
 * @brief Death case: spawning before init_storage() must trap.
 * @details Transformer surface — the slot scan indexes the entity block, so a
 *          spawn on an un-initialized pool would dereference null instead of
 *          reporting the missed init() wiring.
 */
inline void case_transformer_pool_spawn_before_init() {
  Timeline tl;
  RippleTransformer<2> rt(tl);
  Animation::Ripple *p =
      rt.spawn(0, math::Vector(0, 1, 0), 0.2f, 4); // -> HS_CHECK
  if (p == reinterpret_cast<Animation::Ripple *>(0x1))
    std::printf("x");
}

/**
 * @brief Death case: preparing frame state before init_storage() must trap.
 * @details Transformer surface — prepare_frame() is the ordering contract's
 *          other half: an un-initialized pool has no active slots, so it would
 *          silently do nothing and leave the composition reading state that was
 *          never prepared, instead of reporting the missed init() wiring.
 */
inline void case_transformer_pool_prepare_frame_before_init() {
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.prepare_frame(); // -> HS_CHECK
}

/**
 * @brief Death case: a pausable spawn with no pause flag must trap.
 * @details Transformer surface -- spawn_pausable() exists only to hand the
 *          animation a gate to read every frame; a null flag would schedule an
 *          event that can never pause, which is what plain spawn() is for.
 */
inline void case_transformer_pool_spawn_pausable_null_flag() {
  configure_arenas_default();
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);
  Animation::Ripple *p = rt.spawn_pausable(nullptr, 0, math::Vector(0, 1, 0),
                                           0.2f, 4); // -> HS_CHECK
  if (p == reinterpret_cast<Animation::Ripple *>(0x1))
    std::printf("x");
}

/**
 * @brief Death case: an out-of-range active index must trap.
 * @details Transformer surface — active_params() indexes the compact active list,
 *          which is shorter than CAPACITY, so an index taken from the slot domain
 *          (or from a stale count) would read a dead slot as if it were live.
 */
inline void case_transformer_pool_active_index_oob() {
  configure_arenas_default();
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);
  const Animation::RippleParams &p = rt.active_params(opaque(0)); // -> HS_CHECK
  if (p.amplitude == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: reclaimed storage landing at a new address must trap.
 * @details Transformer surface — spawned animations hold Params references into
 *          the slots, so the post-reset replay must re-land the blocks exactly
 *          where init_storage() put them. Here the arena is NOT reset first, so
 *          the replay appends past the originals and every live reference would
 *          be left pointing at abandoned bytes.
 */
inline void case_transformer_pool_reclaim_storage_moved() {
  configure_arenas_default();
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);
  rt.reclaim_storage(persistent_arena); // blocks land elsewhere -> HS_CHECK
}

/**
 * @brief Death case: spawning after the pool's arena was reclaimed must trap.
 * @details Transformer surface — init_storage() must run after
 *          configure_arenas(), which rebinds the arena and hands its bytes out
 *          again. The slot pointers stay non-null across that, so the watermark
 *          is what catches the ordering, in every build.
 */
inline void case_transformer_pool_arena_reclaimed() {
  configure_arenas_default();
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);
  configure_arenas_default(); // rebinds under the live pool
  Animation::Ripple *p =
      rt.spawn(0, math::Vector(0, 1, 0), 0.2f, 4); // -> HS_CHECK
  if (p == reinterpret_cast<Animation::Ripple *>(0x1))
    std::printf("x");
}

inline void case_transformer_pinned_owner_order() {
  configure_arenas_default();
  Timeline timeline;
  using Pool = NoiseTransformer<1>;
  alignas(Pool) static uint8_t storage[sizeof(Pool)];
  Pool *first = new (storage) Pool(timeline);
  Pool second(timeline);
  first->init_storage(persistent_arena);
  second.init_storage(persistent_arena);
  first->spawn_pinned(0);
  second.spawn_pinned(0);
  first->~Pool();
}

/**
 * @brief Death case: a pool outliving its Timeline must trap.
 * @details Transformer surface — the destructor reaches back into the timeline
 *          to drop the pool's clear hook, and the spawned completion callbacks
 *          reach back into the pool, so the two lifetimes are ordered. An owner
 *          that declares them the other way gets a dead reference here rather
 *          than at some later step().
 */
inline void case_transformer_pool_outlives_timeline() {
  configure_arenas_default();
  alignas(Timeline) static uint8_t tl_storage[sizeof(Timeline)];
  Timeline *tl = new (tl_storage) Timeline();
  RippleTransformer<2> rt(*tl);
  rt.init_storage(persistent_arena);
  tl->~Timeline();
  // ~RippleTransformer at scope exit -> HS_CHECK
}
