/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// HankinSolids white-box checks
// ---------------------------------------------------------------------------

/** @brief Test access to HankinSolids' morph-chain state and draw path. */
struct HankinSolidsWhiteBox {
  using EffectT = HankinSolids<SMALL_W, SMALL_H>;

  static bool morph_active(const EffectT &effect) {
    return effect.pending_landing != nullptr;
  }
  static uint8_t node(const EffectT &effect) { return effect.node; }
  static void draw_without_topology(EffectT &effect) {
    effect.hankin_mesh.topology.clear();
    effect.hankin_mesh.topology_key = 0;
    const BakedPalette *palettes[EffectT::NUM_PALETTES]{};
    Canvas canvas(effect);
    effect.draw_mesh(canvas, effect.hankin_mesh, palettes, palettes);
  }
};

/** @brief Verifies a manual Angle edit freezes an in-flight Conway morph. */
inline void test_hankinsolids_manual_pause_holds_morph() {
  reset_effect_globals();
  HankinSolidsWhiteBox::EffectT effect;
  effect.init();

  int frames = 0;
  while (!HankinSolidsWhiteBox::morph_active(effect) && frames < 100) {
    effect.draw_frame();
    effect.advance_display();
    ++frames;
  }
  HS_EXPECT_TRUE(HankinSolidsWhiteBox::morph_active(effect));

  for (int i = 0; i < 5; ++i) {
    effect.draw_frame();
    effect.advance_display();
  }
  const uint8_t held_node = HankinSolidsWhiteBox::node(effect);
  HS_EXPECT_EQ(effect.updateParameter("Angle", 0.7f), ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(effect.animations_paused());

  for (int i = 0; i < 300; ++i) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_TRUE(HankinSolidsWhiteBox::morph_active(effect));
  HS_EXPECT_EQ(HankinSolidsWhiteBox::node(effect), held_node);
  const auto *angle = effect.getParameters().find("Angle");
  HS_EXPECT_TRUE(angle != nullptr);
  if (!angle)
    return;
  HS_EXPECT_NEAR(angle->get(), 0.7f, 1e-6f);

  const uint64_t energy = frame_energy<SMALL_W, SMALL_H>(effect);
  HS_EXPECT_GT(energy, 0u);
}

/**
 * @brief Bounds whole-solid generation, classification and rendering scratch
 * against HankinSolids' exported budgets at the device height.
 */
inline void test_hankinsolids_arena_budget_covers_every_solid() {
  constexpr int W = 288, H = 144;
  constexpr size_t SCRATCH_A = HankinSolids<W, H>::SCRATCH_A_BYTES;
  constexpr size_t SCRATCH_B = HankinSolids<W, H>::SCRATCH_B_BYTES;
  constexpr size_t MEASURE = 1024 * 1024; // headroom so a peak never traps here
  constexpr float ANGLE = math::PI_F / 4.0f;

  auto solids = Solids::Collections::get_simple_solids();
  for (size_t idx = 0; idx < solids.size(); ++idx) {
    configure_arenas(GLOBAL_ARENA_SIZE - 2 * MEASURE, MEASURE, MEASURE);

    MeshPaletteBank palette_bank;
    palette_bank.bake_all(persistent_arena);

    // The effect's held graph-walk seed; dodecahedron is the largest Platonic.
    PolyMesh seed;
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      seed =
          Solids::finalize_solid(Solids::Platonic::dodecahedron(a, b), target);
    });

    MeshState mesh;
    CompiledHankin hankin;
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh base = Solids::finalize_solid(solids[idx].generate(a, b), a);
      hankin = CompiledHankin();
      MeshOps::compile_hankin(base, hankin, target, a);
      mesh.clear();
      MeshOps::update_hankin(hankin, mesh, target, ANGLE);
    });
    {
      ScratchScope a_guard(scratch_arena_a);
      ScratchScope b_guard(scratch_arena_b);
      MeshOps::classify_faces_by_topology(mesh, scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }

    // Render peak: transform into scratch_a, then Scan::Mesh::draw stacks a
    // FaceScratchBuffer on top.
    {
      ScratchScope a_guard(scratch_arena_a);
      math::Orientation<> orientation;
      OrientTransformer camera(orientation);
      MeshState rotated;
      MeshOps::transform(mesh, rotated, scratch_arena_a, camera);
      hs_test::StubEffect fx(W, H);
      Canvas canvas(fx);
      Pipeline<W, H> filters;
      auto frag = [](const math::Vector &, Fragment &f) {
        f.color = Color4(Pixel(1000, 1000, 1000), 1.0f);
      };
      Scan::Mesh::draw<W, H>(filters, canvas, rotated, frag, scratch_arena_a);
    }

    // Morph compaction peak: the CompiledHankin + palette bank survive into
    // scratch_b, the mesh + walk seed into scratch_a, then persistent is
    // reset.
    {
      Persist<CompiledHankin> ph(hankin, scratch_arena_b, persistent_arena);
      Persist<MeshState> pf(mesh, scratch_arena_a, persistent_arena);
      Persist<MeshPaletteBank> pp(palette_bank, scratch_arena_b,
                                  persistent_arena);
      Persist<PolyMesh> ps(seed, scratch_arena_a, persistent_arena);
      persistent_arena.reset();
    }

    const size_t a_peak = scratch_arena_a.get_high_water_mark();
    const size_t b_peak = scratch_arena_b.get_high_water_mark();
    if (a_peak > SCRATCH_A || b_peak > SCRATCH_B)
      std::printf("  HankinSolids arena OVER BUDGET solid[%zu] '%s': "
                  "scratchA=%zu/%zu scratchB=%zu/%zu\n",
                  idx, solids[idx].name, a_peak, SCRATCH_A, b_peak, SCRATCH_B);
    HS_EXPECT_TRUE(a_peak <= SCRATCH_A);
    HS_EXPECT_TRUE(b_peak <= SCRATCH_B);
  }
  configure_arenas_default();
}
