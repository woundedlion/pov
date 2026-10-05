/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_animation.h.

// ============================================================================
// MeshCarousel arena compaction
// ----------------------------------------------------------------------------
// compact_keep_front() evacuates the front slot to a scratch arena, resets the
// persistent arena, and restores on scope exit; it drops the back slot and runs
// after_reset before the front restore. compact_drop_all() evacuates nothing:
// both slots are freed and after_reset sees the empty arena.
// ============================================================================

/**
 * @brief Verifies compact_keep_front() drops the back slot, restores the front,
 * and runs after_reset BEFORE the front restore — a re-bake lands beneath the
 * restored front geometry in the arena.
 */
inline void test_meshcarousel_compact_keep_front_drops_back() {
  reset_persistent_arena();
  static uint8_t polybuf[1 << 14];
  Arena polyarena(polybuf, sizeof(polybuf));
  PolyMesh poly;
  build_octahedron(poly, polyarena);

  MeshCarousel<Segue::Crossfade> carousel; // front slot 0
  MeshOps::compile(poly, carousel.slot(0), persistent_arena, scratch_arena_a);
  MeshOps::compile(poly, carousel.slot(1), persistent_arena, scratch_arena_a);
  const size_t v_front = carousel.current().vertices.size();

  bool after_reset_ran = false;
  const void *bake_ptr = nullptr;
  carousel.compact_keep_front(1, [&](Arena &a) {
    after_reset_ran = true;
    bake_ptr = a.allocate(64);
  });

  HS_EXPECT_TRUE(after_reset_ran);
  HS_EXPECT_TRUE(carousel.current().is_bound());
  HS_EXPECT_EQ(carousel.current().vertices.size(), v_front);
  HS_EXPECT_VEC(carousel.current().vertices[0], math::Vector(1, 0, 0), 1e-5f);
  HS_EXPECT_FALSE(carousel.slot(1).is_bound());
  HS_EXPECT_GT(reinterpret_cast<uintptr_t>(&carousel.current().vertices[0]),
               reinterpret_cast<uintptr_t>(bake_ptr));
}

/**
 * @brief Verifies compact_drop_all() frees both slots and evacuates nothing:
 * after_reset sees an empty persistent arena and its bake is the only thing
 * left in it.
 */
inline void test_meshcarousel_compact_drop_all_frees_both_slots() {
  reset_persistent_arena();
  static uint8_t polybuf[1 << 14];
  Arena polyarena(polybuf, sizeof(polybuf));
  PolyMesh poly;
  build_octahedron(poly, polyarena);

  MeshCarousel<Segue::Crossfade> carousel;
  MeshOps::compile(poly, carousel.slot(0), persistent_arena, scratch_arena_a);
  MeshOps::compile(poly, carousel.slot(1), persistent_arena, scratch_arena_a);
  HS_EXPECT_TRUE(carousel.slot(0).is_bound());
  HS_EXPECT_TRUE(carousel.slot(1).is_bound());
  HS_EXPECT_GT(persistent_arena.get_offset(), (size_t)0);

  bool after_reset_ran = false;
  size_t offset_at_reset = SIZE_MAX;
  const void *bake_ptr = nullptr;
  carousel.compact_drop_all([&](Arena &a) {
    after_reset_ran = true;
    offset_at_reset = a.get_offset();
    bake_ptr = a.allocate(64);
  });

  HS_EXPECT_TRUE(after_reset_ran);
  // Nothing is evacuated, so the callback runs on a fully reclaimed arena.
  HS_EXPECT_EQ(offset_at_reset, (size_t)0);
  HS_EXPECT_TRUE(bake_ptr != nullptr);
  HS_EXPECT_FALSE(carousel.slot(0).is_bound());
  HS_EXPECT_FALSE(carousel.slot(1).is_bound());
  // Only the bake is resident: neither slot was restored on top of it.
  HS_EXPECT_LE(persistent_arena.get_offset(),
               (size_t)64 + alignof(std::max_align_t));
}
