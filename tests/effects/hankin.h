/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// HankinSolids pause coverage
// ---------------------------------------------------------------------------

/** @brief Test access to HankinSolids' morph-chain state. */
struct HankinPauseWhiteBox {
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
  HankinPauseWhiteBox::EffectT effect;
  effect.init();

  int frames = 0;
  while (!HankinPauseWhiteBox::morph_active(effect) && frames < 100) {
    effect.draw_frame();
    effect.advance_display();
    ++frames;
  }
  HS_EXPECT_TRUE(HankinPauseWhiteBox::morph_active(effect));

  for (int i = 0; i < 5; ++i) {
    effect.draw_frame();
    effect.advance_display();
  }
  const uint8_t held_node = HankinPauseWhiteBox::node(effect);
  HS_EXPECT_EQ(effect.updateParameter("Angle", 0.7f), ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(effect.animations_paused());

  for (int i = 0; i < 300; ++i) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_TRUE(HankinPauseWhiteBox::morph_active(effect));
  HS_EXPECT_EQ(HankinPauseWhiteBox::node(effect), held_node);
  const auto *angle = effect.getParameters().find("Angle");
  HS_EXPECT_TRUE(angle != nullptr);
  if (!angle)
    return;
  HS_EXPECT_NEAR(angle->get(), 0.7f, 1e-6f);

  const uint64_t energy = frame_energy<SMALL_W, SMALL_H>(effect);
  HS_EXPECT_GT(energy, 0u);
}
