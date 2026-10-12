/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Transition
// ============================================================================

/**
 * @brief Verifies a Transition steps its bound float from start to target over
 * `duration` frames and reports done() only once it lands on the target.
 */
inline void test_transition_reaches_target_linear() {
  float v = 0.0f;
  const int duration = 10;
  Animation::Transition tr(v, 100.0f, duration, math::ease_linear);
  HS_EXPECT_FALSE(tr.done());
  for (int i = 0; i < duration; ++i) {
    tr.step(fake_canvas());
    HS_EXPECT_NEAR(v, 10.0f * static_cast<float>(i + 1), 1e-3f);
    if (i < duration - 1)
      HS_EXPECT_FALSE(tr.done());
  }
  HS_EXPECT_TRUE(tr.done());
  HS_EXPECT_NEAR(v, 100.0f, 1e-3f);
}

/**
 * @brief Verifies `from` is captured from the live value at the first step and
 * a forward Transition advances monotonically toward the target.
 */
inline void test_transition_monotonic_and_starts_from_current() {
  float v = 5.0f;
  const int duration = 8;
  Animation::Transition tr(v, 25.0f, duration, math::ease_linear);
  v = 9.0f;
  for (int i = 0; i < duration; ++i) {
    tr.step(fake_canvas());
    HS_EXPECT_NEAR(v, 9.0f + 16.0f * (i + 1) / duration, 1e-5f);
  }
  HS_EXPECT_NEAR(v, 25.0f, 1e-3f);
}

/**
 * @brief Verifies a zero-duration Transition coerces to 1 frame: no
 * divide-by-zero, and a single step completes it.
 */
inline void test_transition_duration_zero_no_divide_by_zero() {
  float v = 0.0f;
  Animation::Transition tr(v, 42.0f, 0, math::ease_linear);
  tr.step(fake_canvas());
  HS_EXPECT_TRUE(std::isfinite(v));
  HS_EXPECT_NEAR(v, 42.0f, 1e-3f);
  HS_EXPECT_TRUE(tr.done());
}

/**
 * @brief Verifies a quantized Transition floors the eased value, landing on an
 * integer.
 */
inline void test_transition_quantized_floors_result() {
  float v = 0.0f;
  const int duration = 4;
  Animation::Transition tr(v, 3.7f, duration, math::ease_linear,
                           {.quantized = true});
  for (int i = 0; i < duration; ++i) {
    tr.step(fake_canvas());
  }
  // 3.7 floored to 3.0.
  HS_EXPECT_NEAR(v, 3.0f, 1e-5f);
}

/**
 * @brief Verifies a repeating Transition rewinds t to 0 at the end of each
 * cycle and retraverses the full 0 -> target ramp.
 * @details `from` is captured once (on the first step), so the second cycle
 * must traverse the full ramp again rather than freezing at the target.
 */
inline void test_transition_repeat_retraverses_each_cycle() {
  Timeline tl;
  float v = 0.0f;
  const int duration = 4;
  tl.add(0, Animation::Transition(v, 10.0f, duration, math::ease_linear,
                                  {.repeat = true}));
  std::array<float, duration> first_cycle{};
  for (int i = 0; i < duration; ++i) {
    tl.step(fake_canvas());
    first_cycle[i] = v;
  }
  HS_EXPECT_NEAR(v, 10.0f, 1e-5f);
  for (int i = 0; i < duration; ++i) {
    tl.step(fake_canvas());
    HS_EXPECT_NEAR(v, first_cycle[i], 1e-5f);
  }
}

/**
 * @brief Verifies a wired pause flag freezes a Transition: neither the timer
 * nor the bound float advances, and `from` is snapshotted on the first unpaused
 * step so a slider edit made while paused is honored.
 */
inline void test_transition_paused_holds_value() {
  bool paused = true;
  float v = 0.0f;
  const int duration = 4;
  Animation::Transition tr(v, 10.0f, duration, math::ease_linear,
                           {.paused = &paused});
  for (int i = 0; i < duration; ++i)
    tr.step(fake_canvas());
  HS_EXPECT_NEAR(v, 0.0f, 1e-5f);
  HS_EXPECT_FALSE(tr.done());

  v = 6.0f; // edited while paused; becomes the ramp's start
  paused = false;
  tr.step(fake_canvas());
  HS_EXPECT_NEAR(v, 7.0f, 1e-3f);
  for (int i = 1; i < duration; ++i)
    tr.step(fake_canvas());
  HS_EXPECT_TRUE(tr.done());
  HS_EXPECT_NEAR(v, 10.0f, 1e-3f);
}

/**
 * @brief Verifies a finished finite parameter animation reports full progress
 * rather than dividing by the zeroed duration.
 */
inline void test_finished_param_animation_progress_is_finite() {
  class SelfFinishingParam
      : public Animation::FiniteParamAnimationBase<SelfFinishingParam> {
  public:
    SelfFinishingParam() : FiniteParamAnimationBase(100, true) {}
    void end_now() { finish(); }
    float progress() const { return normalized_progress(); }
  } anim;
  anim.end_now();
  const float progress = anim.progress();
  HS_EXPECT_TRUE(std::isfinite(progress));
  HS_EXPECT_NEAR(progress, 1.0f, 1e-6f);
}

/** @brief Verifies Progress pause behavior and eased output bounds. */
inline void test_progress_pause_and_eased_bounds() {
  bool paused = true;
  std::vector<float> values;
  Animation::Progress progress([&](float t) { values.push_back(t); }, 4,
                               [](float t) { return 2.0f * t; },
                               {.paused = &paused});
  progress.step(fake_canvas());
  HS_EXPECT_TRUE(values.empty());
  paused = false;
  progress.step(fake_canvas());
  HS_EXPECT_SIZE_OR_RETURN(values, 1);
  HS_EXPECT_EQ(values.back(), 0.5f);
  paused = true;
  progress.step(fake_canvas());
  HS_EXPECT_SIZE_OR_RETURN(values, 1);
  paused = false;
  for (int i = 0; i < 3; ++i)
    progress.step(fake_canvas());
  HS_EXPECT_SIZE_OR_RETURN(values, 4);
  HS_EXPECT_EQ(values[1], 1.0f);
  HS_EXPECT_EQ(values[2], 1.0f);
  HS_EXPECT_EQ(values[3], 1.0f);
  HS_EXPECT_TRUE(progress.done());
}
