/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_animation.h.

// ============================================================================
// Sprite (opacity envelope + paused-hold)
// ----------------------------------------------------------------------------
// Sprite::step computes an opacity and forwards it to the user draw_fn; the
// captured opacity drives the envelope and paused-hold assertions.
// ============================================================================

/**
 * @brief Verifies Sprite drives a fade-in -> full-opacity plateau -> fade-out
 * envelope across its duration, completing on the final frame.
 */
inline void test_sprite_fade_in_plateau_fade_out_envelope() {
  std::vector<float> ops;
  const int dur = 10, fade_in = 3, fade_out = 3;
  Animation::Sprite s([&](Canvas &, float o) { ops.push_back(o); }, dur,
                      {.fade_in = {fade_in}, .fade_out = {fade_out}});
  for (int i = 0; i < dur; ++i)
    s.step(fake_canvas()); // observed at t = 1..10

  HS_EXPECT_SIZE_OR_RETURN(ops, dur);
  // Fade-in (linear): t=1 -> 1/3, then rising.
  HS_EXPECT_NEAR(ops[0], 1.0f / 3.0f, 1e-3f);
  HS_EXPECT_GT(ops[1], ops[0]);
  // Plateau: fully opaque between fade-in and fade-out.
  HS_EXPECT_NEAR(ops[3], 1.0f, 1e-3f);
  HS_EXPECT_NEAR(ops[5], 1.0f, 1e-3f);
  // Fade-out: transparent on the final frame.
  HS_EXPECT_LT(ops[8], ops[5]);
  HS_EXPECT_NEAR(ops[9], 0.0f, 1e-3f);
  HS_EXPECT_TRUE(s.done());
}

/** @brief Clamps an overshooting fade-in easing to full opacity. */
inline void test_sprite_clamps_overshooting_fade_in() {
  float opacity = -1.0f;
  Animation::Sprite s([&](Canvas &, float o) { opacity = o; }, 4,
                      {.fade_in = {2, [](float) { return 1.35f; }}});
  s.step(fake_canvas());
  HS_EXPECT_NEAR(opacity, 1.0f, 1e-6f);
}

/**
 * @brief Verifies that when fade_in + fade_out exceed duration the fades scale
 * proportionally into a continuous triangle that still peaks at full opacity.
 * @details There must be no jump where the fade-in hands off to the fade-out;
 * the slider ranges make this configuration user-reachable.
 */
inline void test_sprite_overlapping_fades_stay_continuous() {
  std::vector<float> ops;
  const int dur = 10, fade_in = 8,
            fade_out = 8; // 8 + 8 > 10 => overlap, scaled to 5 + 5
  Animation::Sprite s([&](Canvas &, float o) { ops.push_back(o); }, dur,
                      {.fade_in = {fade_in}, .fade_out = {fade_out}});
  for (int i = 0; i < dur; ++i)
    s.step(fake_canvas()); // observed at t = 1..10

  HS_EXPECT_SIZE_OR_RETURN(ops, dur);
  // Scaled to 5 + 5: slope 1/5 per frame, so no single-frame step exceeds ~0.2.
  float max_jump = 0.0f;
  for (size_t i = 1; i < ops.size(); ++i)
    max_jump = std::max(max_jump, std::abs(ops[i] - ops[i - 1]));
  HS_EXPECT_LT(max_jump, 0.3f);
  // Rises to full opacity at the apex, then falls (a triangle).
  HS_EXPECT_NEAR(ops[4], 1.0f, 1e-3f);
  HS_EXPECT_GT(ops[4], ops[0]);
  HS_EXPECT_GT(ops[4], ops[9]);
}

/**
 * @brief Verifies that while its paused flag is set a Sprite holds its current
 * frame: the timer never advances and it never expires, yet it keeps drawing at
 * its held opacity.
 */
inline void test_sprite_paused_holds_frame() {
  bool paused = true;
  int draws = 0;
  float last_op = -1.0f;
  // No fade in/out => held opacity is full.
  Animation::Sprite s(
      [&](Canvas &, float o) {
        draws++;
        last_op = o;
      },
      /*duration=*/3, {.paused = &paused});

  for (int i = 0; i < 10; ++i)
    s.step(fake_canvas());
  HS_EXPECT_FALSE(s.done()); // timer never advanced while paused
  HS_EXPECT_EQ(draws, 10);   // but it kept drawing every frame
  HS_EXPECT_NEAR(last_op, 1.0f, 1e-6f);

  paused = false;
  for (int i = 0; i < 3; ++i)
    s.step(fake_canvas());
  HS_EXPECT_TRUE(s.done());
}

/** @brief Verifies a Timeline pause keeps a started Sprite drawing in place. */
inline void test_timeline_pause_redraws_held_sprite() {
  Timeline timeline;
  bool paused = false;
  std::vector<float> opacity;
  timeline.add_pausable(
      0,
      Animation::Sprite(
          [&](Canvas &, float value) { opacity.push_back(value); }, 4,
          {.fade_in = {2}}),
      &paused);

  timeline.step(fake_canvas());
  HS_EXPECT_SIZE_OR_RETURN(opacity, 1);
  HS_EXPECT_NEAR(opacity.back(), 0.5f, 1e-6f);
  paused = true;
  for (int i = 0; i < 5; ++i)
    timeline.step(fake_canvas());
  HS_EXPECT_SIZE_OR_RETURN(opacity, 6);
  for (size_t i = 1; i < opacity.size(); ++i)
    HS_EXPECT_NEAR(opacity[i], 0.5f, 1e-6f);

  paused = false;
  timeline.step(fake_canvas());
  HS_EXPECT_NEAR(opacity.back(), 1.0f, 1e-6f);
}

/**
 * @brief Verifies a Sprite paused before its first step holds the opacity its
 * first unpaused frame would report, not the zero end of its fade-in ramp and
 * not a full-brightness plateau.
 * @details The pause holds t at 0, so reporting the ramp there would multiply
 * the consumer's draw to nothing for the whole pause; the effect roster is
 * paused before init() by the WASM load path. Reporting 1.0 instead would hold
 * the sprite brighter than any frame the transition draws.
 */
inline void test_sprite_paused_before_first_step_holds_first_opacity() {
  Timeline timeline;
  bool paused = true;
  std::vector<float> opacity;
  timeline.add_pausable(
      0,
      Animation::Sprite(
          [&](Canvas &, float value) { opacity.push_back(value); }, 16,
          {.fade_in = {8}, .fade_out = {8}}),
      &paused);

  for (int i = 0; i < 4; ++i)
    timeline.step(fake_canvas());
  HS_EXPECT_SIZE_OR_RETURN(opacity, 4);
  for (float value : opacity)
    HS_EXPECT_NEAR(value, 0.125f, 1e-6f);

  // Released, the sprite still runs its fade-in from the start, resuming at the
  // opacity the pause was holding.
  paused = false;
  timeline.step(fake_canvas());
  HS_EXPECT_NEAR(opacity.back(), 0.125f, 1e-6f);
}
