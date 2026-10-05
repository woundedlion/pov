/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

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
 * @brief Recursive generator body: each level calls generate() once more.
 * @param target Arena forwarded to the next generate().
 * @param remaining Levels of nesting still to open; stops at 0.
 * @return Always 0.
 */
inline int nested_generate(Arena &target, Arena &, Arena &, int remaining) {
  if (remaining <= 0)
    return 0;
  return hs::generate(target, nested_generate, remaining - 1);
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
 * @brief Death case: normalizing a degenerate (zero-length) vector must trap.
 * @details Math-core surface — length below epsilon fires the normalize guard.
 */
inline void case_normalize_zero() {
  math::Vector z{opaque(0.0f), opaque(0.0f), opaque(0.0f)};
  math::Vector n = z.normalized(); // length < eps -> HS_CHECK
  if (n.x == 42.0f)
    std::printf("x");
}

/** @brief In-place normalization rejects a zero-length vector. */
inline void case_vector_normalize_in_place_zero() {
  math::Vector v{opaque(0.0f), opaque(0.0f), opaque(0.0f)};
  v.normalize();
  if (v.x == 42.0f)
    std::printf("x");
}

/** @brief Death case: rotating in a degenerate coordinate plane must trap. */
inline void case_rotate_plane_degenerate() {
  math::Mat4 m = math::Mat4::identity();
  math::rotate_plane(m, opaque(1), opaque(1), 0.5f); // a == b -> HS_CHECK
  if (m.m[0][0] == 42.0f)
    std::printf("x");
}

/** @brief Death case: measuring an angle from a degenerate vector must trap. */
inline void case_angle_between_zero() {
  math::Vector zero{opaque(0.0f), opaque(0.0f), opaque(0.0f)};
  math::Vector x{opaque(1.0f), opaque(0.0f), opaque(0.0f)};
  float angle = math::angle_between(zero, x);
  if (angle == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: normalizing a NaN vector must trap.
 * @details Math-core surface — a NaN coordinate poisons the length to NaN, and
 *          `NaN >= epsilon` is false, so the normalize guard fires. The suite's
 *          NaN/Inf fault case: proves a non-finite producer is trapped at the
 *          math seam rather than silently propagating NaN into geometry.
 */
inline void case_normalize_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  math::Vector bad{nan, opaque(0.0f), opaque(0.0f)};
  math::Vector n = bad.normalized(); // length is NaN -> HS_CHECK fails
  if (n.x == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: an out-of-range solids index must trap.
 * @details Lookup/registry surface — get_entry past NUM_ENTRIES fires HS_CHECK.
 */
inline void case_solids_index_oob() {
  const auto &e = Solids::get_entry(opaque<size_t>(Solids::NUM_ENTRIES));
  if (e.name == nullptr)
    std::printf("x");
}

/**
 * @brief Death case: looking up an unknown solid name must trap.
 * @details Registry-by-name surface — an unknown name has no valid fallback.
 */
inline void case_solids_unknown_name() {
  PolyMesh m = Solids::get_by_name(persistent_arena, scratch_arena_a,
                                   scratch_arena_b, "definitely_not_a_solid");
  if (m.vertices.size() == 0x7fff)
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
 * @brief Death case: requesting more KDTree neighbors than MAX_K must trap.
 * @details Spatial surface — k beyond the MAX_K-sized result/heap buffers makes
 *          nearest() trap rather than silently capping the result and masking
 *          the caller's sizing mistake.
 */
inline void case_spatial_knn_over_max() {
  static uint8_t buf[512];
  Arena a(buf, sizeof(buf));
  math::Vector pts[2] = {math::Vector(1.0f, 0.0f, 0.0f),
                         math::Vector(0.0f, 1.0f, 0.0f)};
  KDTree tree(a, std::span<const math::Vector>(pts, 2));
  // Tree is non-empty and k > 0, so the k <= MAX_K guard is reached.
  auto r =
      tree.nearest(math::Vector(1.0f, 0.0f, 0.0f),
                   opaque<size_t>(KDTree::MAX_K + 1)); // k > MAX_K -> HS_CHECK
  if (r.size() == static_cast<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: a lattice index outside [0, RD_N) must trap.
 * @details Spatial surface — node()'s index maps affinely onto the sphere, so an
 *          out-of-range one would silently return a direction off the lattice
 *          instead of naming the caller's mistake.
 */
inline void case_reaction_graph_node_index_out_of_range() {
  math::Vector v = ReactionGraph::node(opaque(ReactionGraph::RD_N));
  if (v.x == 0x7fff)
    std::printf("x");
}

/**
 * @brief Death case: a neighbor-table slot outside the lattice must trap.
 * @details Spatial surface — CubemapLUT's hill-climb and the reaction-diffusion
 *          Laplacian subscript neighbors[] rows unguarded, so validate_neighbors()
 *          traps on a slot that is not a node index before the first such read.
 */
inline void case_reaction_graph_slot_out_of_range() {
  static int16_t table[ReactionGraph::RD_N][ReactionGraph::RD_K] = {};
  table[0][0] = opaque<int16_t>(-1);
  ReactionGraph::validate_neighbors(table);
}

/** @brief Death case: a GS stencil must fit the delayed-write history. */
inline void case_gs_neighbor_exceeds_history() {
  ReactionGraph::NeighborRun run{};
  run.delta[0] = opaque<int16_t>(-145);
  hs_test::effects_tests::GSWhiteBox::validate_physics_neighbors(&run, 1);
}

/** @brief Death case: GS rejects a zero color-noise scale. */
inline void case_gs_color_noise_zero_scale() {
  using WB = hs_test::effects_tests::GSWhiteBox;
  WB::GS gs;
  gs.init();
  WB::set_color_params(gs, 0.0f, opaque(0.0f), 0.2f, 0.4f);
  WB::advance_color_noise(gs);
}

/** @brief Death case: GS rejects a non-finite color-noise scale. */
inline void case_gs_color_noise_nan_scale() {
  using WB = hs_test::effects_tests::GSWhiteBox;
  WB::GS gs;
  gs.init();
  WB::set_color_params(
      gs, 0.0f, opaque(std::numeric_limits<float>::quiet_NaN()), 0.2f, 0.4f);
  WB::advance_color_noise(gs);
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
 * @brief Trivial Cloneable whose clone() allocates from the destination arena,
 *        so a Persist restore measurably grows the persistent arena.
 */
struct PersistProbe {
  uint8_t *storage = nullptr; /**< Stand-in for arena-backed object state. */
  /**
   * @brief Clones by allocating fresh storage from @p arena.
   * @param src Source probe (unused beyond the Cloneable interface).
   * @param dst Destination probe receiving freshly allocated storage.
   * @param arena Arena the clone allocates from.
   */
  static void clone(const PersistProbe &src, PersistProbe &dst, Arena &arena) {
    (void)src;
    dst.storage = static_cast<uint8_t *>(arena.allocate(opaque<size_t>(8)));
  }
};

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
 * @brief Cloneable payload that allocates nothing, so only the distinct-arena
 *        guard can trap a same-arena Persist.
 */
struct FlatProbe {
  int value = 0; /**< Whole payload; clone() copies it without an arena. */
  /**
   * @brief Clones by plain copy, leaving the arena offset untouched.
   * @param src Source probe.
   * @param dst Destination probe.
   * @param arena Unused; the payload needs no storage.
   */
  static void clone(const FlatProbe &src, FlatProbe &dst, Arena &arena) {
    (void)arena;
    dst.value = src.value;
  }
};

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

/**
 * @brief Death case: relocating a retained (pinned) add_get() handle must trap.
 * @details Animation surface — step()'s compaction routes every relocation
 *          through TimelineEvent::move_into, which traps when the event was
 *          handed out via add_get(Pin::PINNED), converting the dangling-handle
 *          hazard into a fail-fast crash instead of silent corruption.
 */
inline void case_timeline_pinned_relocation() {
  TimelineEvent src;
  src.pinned = opaque(true); // as if handed out by add_get(Pin::PINNED)
  TimelineEvent dst;
  src.move_into(dst); // HS_CHECK(!pinned) -> trap
}

/**
 * @brief Death case: relocating into a slot that still owns an animation must
 *        trap.
 * @details Animation surface — move_into overwrites dst.manager/dst.iface, so a
 *          live destination would lose its animation's destructor. step()'s
 *          compaction only ever targets slots it has already vacated; the trap
 *          pins that invariant for every relocation path.
 */
inline void case_timeline_move_into_live_destination() {
  Timeline tl;
  float v = 0.0f;
  tl.add(0, Animation::Transition(v, 1.0f, 10, math::ease_linear));
  tl.add(0, Animation::Transition(v, 1.0f, 10, math::ease_linear));
  global_timeline_events[opaque(0)].move_into(global_timeline_events[1]);
}

/**
 * @brief Death case: a negative timeline delay must trap.
 */
inline void case_timeline_negative_delay() {
  Timeline tl;
  float value = 0.0f;
  tl.add(opaque(-1), Animation::Transition(value, 1.0f, 1, math::ease_linear));
}

/**
 * @brief Death case: a timeline start past UINT32_MAX must trap.
 */
inline void case_timeline_start_overflow() {
  Timeline tl;
  global_timeline_t = opaque<uint32_t>(UINT32_MAX - 1);
  float value = 0.0f;
  tl.add(opaque(2), Animation::Transition(value, 1.0f, 1, math::ease_linear));
}

/**
 * @brief Death case: a pinned animation that COMPLETES must trap.
 * @details Animation surface — the symmetric companion to
 *          case_timeline_pinned_relocation, which guards the relocation path
 *          (move_into). A pinned-but-finite animation that finishes as the
 *          *last* event needs no relocation, so move_into never runs; step()'s
 *          completion branch would otherwise e.destroy() it and dangle the
 *          caller's retained pointer silently. The pin contract is
 *          pinned => infinite, so a pinned animation that naturally completes is
 *          misuse; the completion branch's HS_CHECK traps it. (A deliberate
 *          cancel() is exempt — see is_canceled() — so this case completes
 *          naturally rather than canceling.)
 */
inline void case_timeline_pinned_completion() {
  static hs_test::StubEffect fx(8, 8);
  static Canvas canvas(fx);
  Timeline tl;
  float v = 0.0f;
  // add_get(Pin::PINNED) rejects a finite non-repeating animation up front (see
  // case_timeline_pinned_finite_animation), so the event is marked pinned
  // directly to reach step()'s completion branch. A 1-frame Transition is finite
  // and the sole event, so step() routes it through completion/destroy.
  tl.add(0, Animation::Transition(v, 1.0f, 1, math::ease_linear));
  global_timeline_events[0].pinned = opaque(true);
  tl.step(canvas); // t=1: done() && !repeats() && !canceled, keep=false -> trap
}

/**
 * @brief Death case: pinning a finite, non-repeating animation must trap.
 * @details Animation surface — add_get(Pin::PINNED) promises the caller a pointer
 *          valid across frames, which only holds for an animation that never
 *          completes on its own. The up-front check rejects the misuse at the
 *          add site instead of leaving it to step()'s completion guard, which
 *          fires only once the animation actually finishes.
 */
inline void case_timeline_pinned_finite_animation() {
  Timeline tl;
  float v = 0.0f;
  tl.add_get(0, Animation::Transition(v, 1.0f, opaque(1), math::ease_linear),
             Timeline::Pin::PINNED);
}

/**
 * @brief Death case: dropping a pinned add on a full timeline must trap.
 * @details Animation surface — the capacity guard returns nullptr, but an
 *          add_get(Pin::PINNED) caller retains that pointer across frames and no
 *          call site null-checks it. The guard traps on the pinned case so a
 *          full timeline fails at the add instead of at the first use of the
 *          stored handle.
 */
inline void case_timeline_pinned_add_on_full_timeline() {
  Timeline tl;
  float sink = 0.0f;
  for (int i = 0; i < Timeline::MAX_EVENTS; ++i)
    tl.add(0, Animation::Transition(sink, 1.0f, 10, math::ease_linear));
  tl.add_get(0,
             Animation::PeriodicTimer(
                 1, [](Canvas &) {}, /*repeat=*/true),
             opaque(Timeline::Pin::PINNED));
}

/**
 * @brief Death case: a pinned one-shot timer must trap when it fires.
 * @details Animation surface — a one-shot RandomTimer/PeriodicTimer ends itself
 *          on its single trigger. Ending via finish() (not cancel()) keeps
 *          is_canceled() false, so the destroy of a pinned timer hits step()'s
 *          completion guard instead of slipping through its cancellation
 *          exemption and dangling the retained pointer silently.
 */
inline void case_timeline_pinned_one_shot_timer() {
  static hs_test::StubEffect fx(8, 8);
  static Canvas canvas(fx);
  Timeline tl;
  tl.add(0, Animation::PeriodicTimer(1, [](Canvas &) {}, /*repeat=*/false));
  global_timeline_events[0].pinned = opaque(true);
  tl.step(canvas); // t=1: fires, finish() -> done() && !canceled -> trap
}
