/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Memory death fixtures and guard cases.

/**
 * @brief Death case: arena over-allocation must trap.
 * @details Memory surface — requests more than the arena's capacity so
 *          allocate() fires HS_CHECK.
 */
inline void case_arena_oom() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  void *p = a.allocate(opaque<size_t>(1024)); // > capacity -> HS_CHECK
  if (p == reinterpret_cast<void *>(0x1))
    std::printf("x"); // keep the call live
}

/** @brief Death case: make() must preserve Arena's OOM trap contract. */
inline void case_arena_make_oom() {
  static uint8_t buf[sizeof(uint32_t) - 1];
  Arena a(buf, sizeof(buf));
  auto *p = a.make<uint32_t>(opaque<uint32_t>(7));
  if (p == reinterpret_cast<uint32_t *>(0x1))
    std::printf("x");
}

/**
 * @brief Death case: a zero-size arena allocation must trap.
 * @details Memory surface — a zero-size request returns a bump pointer that
 *          reserves nothing and aliases the next allocation's address, so it is
 *          rejected as misuse rather than handed back as ownable storage.
 */
inline void case_arena_zero_size_alloc() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  void *p = a.allocate(opaque<size_t>(0)); // size == 0 -> HS_CHECK
  if (p == reinterpret_cast<void *>(0x1))
    std::printf("x");
}

/** @brief Death case: an overflowing typed allocation must trap. */
inline void case_arena_allocate_n_overflow() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  auto *p = a.allocate_n<uint64_t>(opaque<size_t>(SIZE_MAX / 8 + 1));
  if (p == reinterpret_cast<uint64_t *>(0x1))
    std::printf("x");
}

/**
 * @brief Death case: a non-power-of-two allocation alignment must trap.
 * @details Memory surface — allocate()'s padding math is a modulo against the
 *          requested alignment, which only yields an aligned address for a
 *          power of two.
 */
inline void case_arena_bad_alignment() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  void *p = a.allocate(opaque<size_t>(8), opaque<size_t>(3)); // -> HS_CHECK
  if (p == reinterpret_cast<void *>(0x1))
    std::printf("x");
}

/**
 * @brief Death case: a mid-run resplit with live scratch content must trap.
 * @details Config surface — resplit_arenas rebases both scratch arenas, and a
 *          ScratchScope saved at offset 0 restores to 0 either way, so live
 *          scratch content would be silently rebased onto the new split.
 */
inline void case_resplit_scratch_not_empty() {
  configure_arenas_default();
  scratch_arena_a.allocate(opaque<size_t>(16));
  resplit_arenas(opaque(DEFAULT_PERSISTENT_SIZE),
                 opaque(DEFAULT_SCRATCH_A_SIZE),
                 opaque(DEFAULT_SCRATCH_B_SIZE)); // scratch live -> HS_CHECK
}

/**
 * @brief Death case: a resplit below the persistent arena's live offset must
 *        trap.
 * @details Config surface — resplit_arenas keeps the persistent arena's base,
 *          offset and content and only moves its capacity, so a budget under
 *          the live offset would strand the carousel and palette bank.
 */
inline void case_resplit_persistent_strands() {
  configure_arenas_default();
  persistent_arena.allocate(opaque<size_t>(4096));
  resplit_arenas(opaque<size_t>(1024), opaque(DEFAULT_SCRATCH_A_SIZE),
                 opaque(DEFAULT_SCRATCH_B_SIZE)); // offset > budget -> HS_CHECK
}

/**
 * @brief Death case: moving the arena offset forward must trap.
 * @details Memory surface — set_offset only ever rewinds. A forward move stays
 *          inside capacity yet hands back bytes already reclaimed, so the guard
 *          is monotone decrease, not a capacity bound.
 */
inline void case_arena_set_offset_forward() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  a.allocate(opaque<size_t>(8), opaque<size_t>(1));
  a.set_offset(opaque<size_t>(0));
  a.set_offset(opaque<size_t>(4)); // forward, still under capacity -> HS_CHECK
}

/**
 * @brief Death case: non-LIFO ScratchScope teardown must trap.
 * @details The scratch-arena sharing contract between Pixel::Feedback::flush and
 *          Plot::rasterize is safe because scratch_arena_a is a LIFO bump
 *          allocator — but only while scopes are torn down in stack order.
 *          ~ScratchScope enforces that: an outer scope rewinding while an inner
 *          one is still live leaves the arena offset below the inner's saved
 *          mark, and the inner's destructor HS_CHECKs offset >= saved_offset.
 *          Here the outer scope is destroyed first (std::optional::reset),
 *          rewinding to 0; destroying the inner then sees offset 0 < its saved
 *          mark and traps.
 */
inline void case_scratch_scope_non_lifo() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  std::optional<ScratchScope> outer;
  outer.emplace(a);              // saves offset 0
  a.allocate(opaque<size_t>(8)); // offset -> 8
  std::optional<ScratchScope> inner;
  inner.emplace(a); // saves offset 8
  outer.reset();    // non-LIFO: rewinds offset to 0
  inner.reset();    // offset 0 < saved 8 -> HS_CHECK
}

inline void case_scratch_scope_reset() {
  static uint8_t storage[64];
  Arena arena(storage, sizeof(storage));
  ScratchScope scope(arena);
  arena.reset();
}

/** @brief Death case: exhausted arena rewind history must trap. */
inline void case_arena_rewind_history_overflow() {
  static uint8_t storage[1024];
  Arena arena(storage, sizeof(storage));
  for (size_t i = 0; i <= 256; ++i) {
    arena.allocate(2, 1);
    arena.set_offset(i);
  }
}

/**
 * @brief Death case: ArenaVector fixed-capacity push_back overflow must trap.
 * @details Arena-container surface — a push_back past capacity fires HS_CHECK.
 */
inline void case_arena_vector_overflow() {
  static uint8_t buf[256];
  Arena a(buf, sizeof(buf));
  ArenaVector<int> v(a, 2);
  v.push_back(1);
  v.push_back(2);
  v.push_back(opaque(3)); // exceeds capacity -> HS_CHECK
}

/**
 * @brief Death case: ArenaVector fixed-capacity emplace_back overflow must trap.
 * @details Arena-container surface — the in-place construction path carries its
 *          own capacity guard, distinct from push_back's copy path.
 */
inline void case_arena_vector_emplace_overflow() {
  static uint8_t buf[256];
  Arena a(buf, sizeof(buf));
  ArenaVector<int> v(a, 2);
  v.emplace_back(1);
  v.emplace_back(2);
  v.emplace_back(opaque(3)); // exceeds capacity -> HS_CHECK
}

/**
 * @brief Death case: generate() with a scratch arena as its target must trap.
 * @details Generator surface — the depth-0 reset and the ScratchScope rewind
 *          would destroy output written into either engine scratch arena, so an
 *          aliasing target is rejected before the callback runs.
 */
inline void case_generate_target_is_scratch() {
  configure_arenas_default();
  int r = hs::generate(scratch_arena_a, [](Arena &, Arena &, Arena &) {
    return 0;
  }); // target aliases scratch_arena_a -> HS_CHECK
  if (r == 42)
    std::printf("x");
}

/**
 * @brief Death case: nesting generate() past MAX_GENERATE_DEPTH must trap.
 * @details Generator surface — every level stacks two ScratchScopes on a fixed
 *          scratch budget, so runaway reentrancy is capped at the wrapper rather
 *          than left to exhaust the arenas. The outermost call opens depth 1, so
 *          MAX_GENERATE_DEPTH further levels reach depth MAX_GENERATE_DEPTH + 1.
 */
inline void case_generate_recursion_too_deep() {
  configure_arenas_default();
  int r = hs::generate(persistent_arena, nested_generate,
                       opaque(hs::MAX_GENERATE_DEPTH));
  if (r == 42)
    std::printf("x");
}

/**
 * @brief Death case: a StaticCircularBuffer index past the live count must trap.
 * @details Container surface — index >= count fires HS_CHECK.
 */
inline void case_circular_buffer_oob() {
  StaticCircularBuffer<int, 4> cb;
  cb.push_back(10);
  cb.push_back(20);
  int v = cb[opaque<size_t>(5)]; // index >= count -> HS_CHECK
  if (v == 42)
    std::printf("x");
}

/**
 * @brief Death case: front() on an empty StaticCircularBuffer must trap.
 * @details Container surface — the never-taken opaque(false) push keeps the
 *          optimizer from proving the buffer empty and folding the trap at
 *          compile time; is_empty() fires HS_CHECK.
 */
inline void case_circular_buffer_front_empty() {
  StaticCircularBuffer<int, 4> cb;
  if (opaque(false))
    cb.push_back(1);
  int v = cb.front(); // is_empty() -> HS_CHECK("front() on empty ...")
  if (v == 42)
    std::printf("x");
}

inline void case_circular_buffer_back_empty() {
  StaticCircularBuffer<int, 4> cb;
  if (opaque(false))
    cb.push_back(1);
  if (cb.back() == 42)
    std::printf("x");
}

inline void case_circular_buffer_const_back_empty() {
  StaticCircularBuffer<int, 4> cb;
  if (opaque(false))
    cb.push_back(1);
  const auto &view = cb;
  if (view.back() == 42)
    std::printf("x");
}

/**
 * @brief Death case: ArenaVector::append_bulk past its fixed capacity must trap.
 * @details Memory surface — a distinct seam from element-at-a-time push_back;
 *          the bulk memcpy path has its own remaining-capacity guard.
 */
inline void case_arena_vector_append_bulk_overflow() {
  static uint8_t buf[256];
  Arena a(buf, sizeof(buf));
  ArenaVector<int> v(a, 2); // exact capacity 2
  int src[4] = {1, 2, 3, 4};
  v.append_bulk(src, opaque<size_t>(4)); // 0 + 4 > 2 -> HS_CHECK
  if (v.size() == 0x7fff)
    std::printf("x");
}

/**
 * @brief Death case: an over-subscribed arena partition must trap.
 * @details Config surface — each request alone fits but the sum exceeds
 *          GLOBAL_ARENA_SIZE, so configure_arenas fires HS_CHECK.
 */
inline void case_arena_oversubscribed() {
  configure_arenas(opaque(GLOBAL_ARENA_SIZE), opaque<size_t>(1024),
                   opaque<size_t>(1024));
}

/** @brief Death case: the ArenaSplit scratch pair exceeds the block total. */
inline void case_arena_split_scratch_too_large() {
  const ArenaSplit split{opaque<size_t>(GLOBAL_ARENA_SIZE), opaque<size_t>(1)};
  (void)split.persistent(opaque<size_t>(GLOBAL_ARENA_SIZE));
}

/**
 * @brief Death case: a single partition larger than the whole block must trap.
 * @details Config surface — the per-request bound is checked before the sum, so
 *          an oversized persistent request fires split_bases' own HS_CHECK.
 */
inline void case_arena_partition_too_large() {
  configure_arenas(opaque(GLOBAL_ARENA_SIZE + 1), opaque<size_t>(0),
                   opaque<size_t>(0));
}

/**
 * @brief Death case: a Persist scope that forgets persistent_arena.reset() must trap.
 * @details Memory surface — without the rewind, ~Persist's restore clones the
 *          backup *after* the still-live object instead of over it, pushing the
 *          persistent offset past the construction watermark; the post-restore
 *          HS_CHECK fires.
 */
inline void case_persist_forgot_reset() {
  static uint8_t pbuf[256];
  static uint8_t sbuf[256];
  Arena persistent(pbuf, sizeof(pbuf));
  Arena scratch(sbuf, sizeof(sbuf));
  PersistProbe target;
  PersistProbe::clone(target, target,
                      persistent); // the live object in persistent
  {
    Persist<PersistProbe> p(target, scratch, persistent);
    // A correct scope rewinds here (persistent.reset()); omitting it makes the
    // restore append past the watermark.
  } // ~Persist restore -> offset past watermark -> HS_CHECK
  if (target.storage == reinterpret_cast<uint8_t *>(0x1))
    std::printf("x");
}

/**
 * @brief Death case: a Persist naming one arena for both roles must trap.
 * @details Memory surface — ~Persist's watermark restore assumes the backup
 *          outlives the rewind of the arena it restores into, which a single
 *          arena cannot provide. The payload allocates nothing, so the
 *          post-restore watermark check cannot fire and the distinct-arena
 *          guard is the only reachable trap.
 */
inline void case_persist_same_arena() {
  static uint8_t pbuf[256];
  Arena persistent(pbuf, sizeof(pbuf));
  FlatProbe target;
  target.value = opaque(7);
  Persist<FlatProbe> p(target, persistent, persistent); // -> HS_CHECK
  if (target.value == 42)
    std::printf("x");
}

/**
 * @brief Death case: a swapped (unordered) TriangularBitset pair must trap.
 * @details Memory-safety surface — index() requires small < large < MAX_V; a
 *          swapped pair would alias the wrong bit and an out-of-range one would
 *          write adjacent memory, so the HS_CHECK traps the misuse on the cold
 *          edge-dedup setup path.
 */
inline void case_triangular_bitset_unordered_pair() {
  TriangularBitset<128> bits;
  bool hit = bits.test(opaque(5), opaque(3)); // small > large -> HS_CHECK
  if (hit)
    std::printf("x");
}
