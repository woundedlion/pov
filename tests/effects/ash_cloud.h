/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief AshCloud value-cutout cases. */

/**
 * @brief Pins that AshCloud's value cutout reaches the rendered frame.
 * @details Threshold 0 lights pixels; threshold 1 blanks the frame.
 */
inline void test_ash_cloud_value_cutout_gates_the_frame() {
  using FX = AshCloud<SMALL_W, SMALL_H>;

  auto lit_pixels = [](float threshold) {
    reset_effect_globals();
    FX effect;
    effect.init();

    auto snapshot = effect.serialize_parameters();
    snapshot.params.template get<"value">().cutout_threshold = threshold;
    snapshot.params.template get<"value">().cutout_softness = 1.0f / 1024.0f;
    HS_EXPECT_TRUE(effect.restore_parameters(snapshot));

    effect.draw_frame();
    effect.advance_display();

    return count_lit_region<SMALL_W, SMALL_H>(effect);
  };

  HS_EXPECT_GT(lit_pixels(0.0f), size_t{0});
  HS_EXPECT_EQ(lit_pixels(1.0f), size_t{0});
}
