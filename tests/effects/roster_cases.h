/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Cross-effect control cases. */

/**
 * @brief Checks the full-frame query on four constructed effects.
 * @details MeshFeedback and Dynamo fold any_crosses_segments true; the two
 *          sampled non-crossing effects keep segment clipping.
 */
inline void test_needs_full_frame_gate() {
  // Each effect aliases the same static double buffer, so construct one at a
  // time.
  auto reset = [] { reset_effect_globals(); };
  auto check = [&](Effect &fx, bool expected, const char *name) {
    fx.init();
    if (fx.needs_full_frame() != expected)
      std::printf("  [FAIL] %s needs_full_frame()=%d expected=%d\n", name,
                  fx.needs_full_frame(), expected);
    HS_EXPECT_EQ(fx.needs_full_frame(), expected);
  };

  // Cross-segment stateful effects -> full-frame.
  {
    reset();
    MeshFeedback<DEFAULT_W, DEFAULT_H> fx;
    check(fx, true, "MeshFeedback");
  }
  {
    reset();
    Dynamo<DEFAULT_W, DEFAULT_H> fx;
    check(fx, true, "Dynamo");
  }

  // Representative non-stateful effects -> keep band clipping (default false).
  {
    reset();
    Voronoi<DEFAULT_W, DEFAULT_H> fx;
    check(fx, false, "Voronoi");
  }
  {
    reset();
    RingSpin<DEFAULT_W, DEFAULT_H> fx;
    check(fx, false, "RingSpin");
  }
}

template <typename EffectT> inline void check_manual_preset_navigation() {
  constexpr size_t expected_count = EffectT::authored_preset_count();
  HS_EXPECT_GT(expected_count, size_t(1));
  reset_effect_globals();
  EffectT effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), expected_count);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t(0));
  HS_EXPECT_FALSE(effect.selectPreset(expected_count));
  HS_EXPECT_FALSE(effect.animations_paused());
  HS_EXPECT_TRUE(effect.previousPreset());
  HS_EXPECT_EQ(effect.getPresetIndex(), expected_count - 1);
  HS_EXPECT_TRUE(effect.nextPreset());
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t(0));
  HS_EXPECT_TRUE(effect.animations_paused());
}

inline void test_manual_preset_navigation() {
  check_manual_preset_navigation<MindSplatter<SMALL_W, SMALL_H>>();
  check_manual_preset_navigation<DreamBalls<SMALL_W, SMALL_H>>();
  check_manual_preset_navigation<Comets<SMALL_W, SMALL_H>>();
  check_manual_preset_navigation<MeshFeedback<SMALL_W, SMALL_H>>();
  check_manual_preset_navigation<ShapeShifter<SMALL_W, SMALL_H>>();
}
