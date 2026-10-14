/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Settle correspondence: relax output vertex order is the identity over
// its input, so a relaxed endpoint is per-vertex slerpable.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies relax(50) on a registry form (expand(dodecahedron), the
 *        rhombicosidodecahedron chain) preserves vertex order and topology.
 * @details Counts equal, face_counts/faces byte-identical, and each relaxed
 *          vertex stays strictly nearest to its own input vertex.
 */
inline void test_relax_is_vertex_order_identity() {
  Arena target(morph_target_buf, sizeof(morph_target_buf));
  Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
  Arena aux(morph_aux_buf, sizeof(morph_aux_buf));

  PolyMesh dodeca;
  build_solid<Solids::Dodecahedron>(dodeca, temp);
  PolyMesh unrelaxed = MeshOps::expand(dodeca, target, temp);
  PolyMesh relaxed = MeshOps::relax(unrelaxed, aux, temp, 50);

  check_vertex_order_identity(relaxed, unrelaxed);
}

/** @brief Relax and reconcile retain per-face palettes throughout and at landing. */
inline void test_position_legs_preserve_face_palettes() {
  using Animation::OpLeg;
  Arena target(morph_target_buf, sizeof(morph_target_buf));
  Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
  Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
  Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  PolyMesh dodeca;
  build_solid<Solids::Dodecahedron>(dodeca, temp);
  PolyMesh seed = MeshOps::expand(dodeca, target, temp);
  PolyMesh relaxed = MeshOps::relax(seed, aux, temp, 50);
  {
    ScratchScope a(scratch_arena_a);
    ScratchScope b(scratch_arena_b);
    MeshOps::classify_faces_by_topology(relaxed, scratch_arena_a,
                                        scratch_arena_b, aux);
  }
  const size_t FACES = seed.face_counts.size();
  std::vector<uint8_t> palettes(FACES);
  for (size_t f = 0; f < FACES; ++f)
    palettes[f] = static_cast<uint8_t>(f % 2);

  constexpr int FRAMES = 24;
  int drawn = 0;
  auto draw = [&](Canvas &, const MeshState &, const OpLeg::Shading &sh) {
    ++drawn;
    HS_EXPECT_EQ(sh.faces, FACES);
    for (size_t f = 0; f < FACES; ++f) {
      for (float t : {0.0f, 0.5f, 1.0f}) {
        const Color4 got = sh.ramps[sh.face_ramp[f]].get(t);
        const Color4 expected = bank.bank.entries[palettes[f]].get(t);
        HS_EXPECT_EQ(got.color.r, expected.color.r);
        HS_EXPECT_EQ(got.color.g, expected.color.g);
        HS_EXPECT_EQ(got.color.b, expected.color.b);
      }
    }
  };
  OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                .prev_face_palette = palettes.data(),
                                .prev_faces = FACES,
                                .preserve_without_births = true};
  auto run = [&](OpLeg &&leg) {
    drawn = 0;
    const auto &landing = leg.landing();
    HS_EXPECT_EQ(landing.faces, FACES);
    bool same_class_different_palette = false;
    for (size_t f = 0; f < FACES; ++f) {
      HS_EXPECT_EQ(landing.from_palette[f], palettes[f]);
      HS_EXPECT_EQ(landing.landed_palette(f), palettes[f]);
      HS_EXPECT_EQ(landing.topology[f], relaxed.topology[f]);
      for (size_t other = 0; other < f; ++other)
        same_class_different_palette |=
            landing.topology[f] == landing.topology[other] &&
            palettes[f] != palettes[other];
    }
    HS_EXPECT_TRUE(same_class_different_palette);

    hs_test::StubEffect fx(288, 144);
    for (int frame = 0; frame < FRAMES; ++frame) {
      {
        Canvas c(fx);
        leg.step(c);
      }
      fx.advance_display();
    }
    HS_EXPECT_EQ(drawn, FRAMES);
  };
  run(OpLeg(seed, OpLeg::RelaxSpec{.iterations = 50, .sweep_frames = FRAMES},
            target, draw, handoff));
  run(OpLeg(seed,
            OpLeg::ReconcileSpec{.to_positions = relaxed.vertices.data(),
                                 .to_count = relaxed.vertices.size(),
                                 .sweep_frames = FRAMES},
            target, draw, handoff));
}
