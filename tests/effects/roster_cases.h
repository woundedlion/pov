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

/** @brief Parameter lints reject overwritten sliders and paused animation leaks. */
inline void test_parameter_lint_failure_guards() {
  class SliderEffect : public Effect {
  public:
    explicit SliderEffect(bool broken) : Effect(1, 1), broken(broken) {
      register_param("Value", &value, 0.0f, 1.0f);
    }
    void draw_frame() override {
      Canvas canvas(*this);
      if (broken)
        value = 0.0f;
    }
    float value = 0.0f;
    bool broken;
  };
  class AnimatedEffect : public Effect {
  public:
    explicit AnimatedEffect(bool broken) : Effect(1, 1), broken(broken) {
      register_animated_param("Value", &value, 0.0f, 1.0f);
    }
    void draw_frame() override {
      Canvas canvas(*this);
      if (frames++ > 0 && (broken || !animations_paused()))
        value += 0.0001f;
    }
    float value = 0.0f;
    int frames = 0;
    bool broken;
  };
  const auto check = [](auto lint, int expected, const char *diagnostic) {
    const auto saved = stats();
    const auto saved_budget = fail_print_budget();
    fail_print_budget() = {};
    int failures = 0;
    const auto output = capture_stdout([&] {
      lint();
      failures = stats().failed - saved.failed;
    });
    stats() = saved;
    fail_print_budget() = saved_budget;
    HS_EXPECT_TRUE(output.has_value());
    HS_EXPECT_EQ(failures, expected);
    if (output && diagnostic)
      HS_EXPECT_NE(output->find(diagnostic), std::string::npos);
  };
  for (bool broken : {false, true}) {
    reset_globals();
    {
      SliderEffect slider(broken);
      check([&] { lint_dead_sliders(slider, "SliderEffect"); }, int(broken),
            broken ? "DEAD SLIDER" : nullptr);
    }
    reset_globals();
    AnimatedEffect animated(broken);
    check([&] { lint_animated_pause(animated, "AnimatedEffect"); }, int(broken),
          broken ? "PAUSE LEAK" : nullptr);
  }
}

/**
 * @brief Checks the full-frame query on four constructed effects.
 * @details MeshFeedback and Dynamo fold any_crosses_segments true; the two
 *          sampled non-crossing effects keep segment clipping.
 */
inline void test_needs_full_frame_gate() {
  // Each effect aliases the same static double buffer, so construct one at a
  // time.
  auto check = [&](Effect &fx, bool expected, const char *name) {
    HS_CONTEXT(name);
    fx.init();
    HS_EXPECT_EQ(fx.needs_full_frame(), expected);
  };

  // Cross-segment stateful effects -> full-frame.
  {
    reset_globals();
    MeshFeedback<DEFAULT_W, DEFAULT_H> fx;
    check(fx, true, "MeshFeedback");
  }
  {
    reset_globals();
    Dynamo<DEFAULT_W, DEFAULT_H> fx;
    check(fx, true, "Dynamo");
  }

  // Representative non-stateful effects -> keep band clipping (default false).
  {
    reset_globals();
    SphericalHarmonics<DEFAULT_W, DEFAULT_H> fx;
    check(fx, false, "SphericalHarmonics");
  }
  {
    reset_globals();
    RingSpin<DEFAULT_W, DEFAULT_H> fx;
    check(fx, false, "RingSpin");
  }
}

template <typename EffectT> inline void check_manual_preset_navigation() {
  constexpr size_t EXPECTED_COUNT = EffectT::authored_preset_count();
  HS_EXPECT_GT(EXPECTED_COUNT, size_t(1));
  reset_globals();
  EffectT effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), EXPECTED_COUNT);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t(0));
  HS_EXPECT_FALSE(effect.selectPreset(EXPECTED_COUNT));
  HS_EXPECT_FALSE(effect.animations_paused());
  HS_EXPECT_TRUE(effect.previousPreset());
  HS_EXPECT_EQ(effect.getPresetIndex(), EXPECTED_COUNT - 1);
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
