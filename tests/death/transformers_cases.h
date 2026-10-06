/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Transformers death fixtures and guard cases.

/**
 * @brief Death case: a second TransformerPool::init_storage() must trap.
 * @details Spawned animations hold Params references into the first block.
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

/** @brief Death case: spawning before init_storage() must trap. */
inline void case_transformer_pool_spawn_before_init() {
  Timeline tl;
  RippleTransformer<2> rt(tl);
  Animation::Ripple *p =
      rt.spawn(0, math::Vector(0, 1, 0), 0.2f, 4); // -> HS_CHECK
  if (p == reinterpret_cast<Animation::Ripple *>(0x1))
    std::printf("x");
}

/** @brief Death case: preparing frame state before init_storage() must trap. */
inline void case_transformer_pool_prepare_frame_before_init() {
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.prepare_frame(); // -> HS_CHECK
}

/** @brief Death case: a pausable spawn with no pause flag must trap. */
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
 * @details active_params() indexes the compact active list, which is shorter
 *          than CAPACITY.
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
 * @details Spawned animations hold Params references into the slots. The
 *          arena is not reset first, so the replay appends past the originals.
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
 * @details The slot pointers stay non-null across configure_arenas(), so the
 *          watermark catches the ordering.
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
 * @details The destructor reaches back into the timeline to drop the pool's
 *          clear hook.
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
