/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Cross-effect control cases. */

/** @brief Seam guards reject dark coverage and retain measured-frame ratios. */
inline void test_mesh_seam_dark_sample_guards() {
  constexpr int W = 16, H = 2;
  std::vector<Pixel> dark(W * H);
  std::vector<Pixel> lit(W * H);
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const int triangle = x <= W / 2 ? x : W - x;
      lit[y * W + x] =
          Pixel(static_cast<uint16_t>(1000 + 1000 * triangle), 0, 0);
    }
  const auto check_guard = [](const MeshSeamSteps &steps, int expected_failures,
                              const char *diagnostic) {
    const auto saved = stats();
    const auto saved_budget = fail_print_budget();
    fail_print_budget() = {};
    int failures = 0;
    const auto output = capture_stdout([&] {
      expect_mesh_seam_continuous(steps, "coverage companion");
      failures = stats().failed - saved.failed;
    });
    stats() = saved;
    fail_print_budget() = saved_budget;
    HS_EXPECT_TRUE(output.has_value());
    HS_EXPECT_EQ(failures, expected_failures);
    if (output && diagnostic)
      HS_EXPECT_NE(output->find(diagnostic), std::string::npos);
  };

  MeshSeamSteps measured;
  for (int frame = 0; frame < 8; ++frame)
    measured.add(lit, lit, W, H);
  HS_EXPECT_EQ(measured.arm_samples, 8);
  HS_EXPECT_EQ(measured.interior_samples, 8);
  check_guard(measured, 0, nullptr);

  MeshSeamSteps black;
  for (int frame = 0; frame < 8; ++frame)
    black.add(dark, dark, W, H);
  HS_EXPECT_EQ(black.frames, 8);
  HS_EXPECT_EQ(black.arm_samples, 0);
  HS_EXPECT_EQ(black.interior_samples, 0);
  check_guard(black, 4, "steps.arm_samples * 4");

  MeshSeamSteps mostly_dark;
  for (int frame = 0; frame < 8; ++frame)
    mostly_dark.add(frame < 2 ? lit : dark, frame < 2 ? lit : dark, W, H);
  HS_EXPECT_EQ(mostly_dark.frames, 8);
  HS_EXPECT_EQ(mostly_dark.arm_samples, 2);
  HS_EXPECT_EQ(mostly_dark.interior_samples, 2);
  check_guard(mostly_dark, 2, "steps.arm_samples * 4");

  auto spike = lit;
  for (int y = 0; y < H; ++y) {
    spike[y * W].r -= 500;
    spike[y * W + W / 2].r += 500;
  }
  MeshSeamSteps mixed;
  for (int frame = 0; frame < 4; ++frame)
    mixed.add(frame < 3 ? lit : dark, frame < 3 ? spike : dark, W, H);
  HS_EXPECT_EQ(mixed.frames, 4);
  HS_EXPECT_EQ(mixed.arm_samples, 3);
  HS_EXPECT_EQ(mixed.interior_samples, 3);
  HS_EXPECT_NEAR(mixed.mean_arm_ratios / mixed.arm_samples, 24.0 / 17.0, 1e-12);
  check_guard(mixed, 2, "steps.arm_ratios / steps.arm_samples");
}

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
    SphericalHarmonics<DEFAULT_W, DEFAULT_H> fx;
    check(fx, false, "SphericalHarmonics");
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
