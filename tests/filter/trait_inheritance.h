/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_filter.h.

// ============================================================================
// Trait inheritance on representative filters
// ============================================================================

/**
 * @brief Verifies representative filters inherit the right is_2d / has_history
 *        traits, at both runtime and compile time.
 */
inline void test_filter_trait_inheritance() {
  constexpr int W = 32, H = 32;

  // 2D screen-space, stateless.
  HS_EXPECT_TRUE((Filter::Screen::AntiAlias<W, H>::is_2d));
  HS_EXPECT_FALSE((Filter::Screen::AntiAlias<W, H>::has_history));
  HS_EXPECT_TRUE((Filter::Screen::Blur<W, H>::is_2d));
  HS_EXPECT_FALSE((Filter::Screen::Blur<W, H>::has_history));

  // Pixel-space, stateless (tagged IsPixel).
  HS_EXPECT_TRUE((Filter::Pixel::ChromaticShift<W>::is_2d));
  HS_EXPECT_FALSE((Filter::Pixel::ChromaticShift<W>::has_history));

  // Pixel-space feedback is 2D *with* history.
  HS_EXPECT_TRUE((Filter::Pixel::Feedback<W, H>::is_2d));
  HS_EXPECT_TRUE((Filter::Pixel::Feedback<W, H>::has_history));
  HS_EXPECT_EQ((Filter::World::Replicate<W>::domain_rank), 0);
  HS_EXPECT_EQ((Filter::Screen::AntiAlias<W, H>::domain_rank), 1);
  HS_EXPECT_EQ((Filter::Pixel::ChromaticShift<W>::domain_rank), 2);

  // 3D world-space, stateless.
  HS_EXPECT_FALSE((Filter::World::Replicate<W>::is_2d));
  HS_EXPECT_FALSE((Filter::World::Replicate<W>::has_history));
  HS_EXPECT_FALSE((Filter::World::Hole::is_2d));
  HS_EXPECT_FALSE((Filter::World::Hole::has_history));

  // History-bearing trail filters.
  HS_EXPECT_FALSE((Filter::World::Trails<16>::is_2d));
  HS_EXPECT_TRUE((Filter::World::Trails<16>::has_history));
  HS_EXPECT_TRUE((Filter::Screen::Trails<>::is_2d));
  HS_EXPECT_TRUE((Filter::Screen::Trails<>::has_history));

  // is_pipeline separates a stage from a whole pipeline.
  HS_EXPECT_FALSE((Filter::Screen::AntiAlias<W, H>::is_pipeline));
  HS_EXPECT_FALSE((Filter::World::Replicate<W>::is_pipeline));
  HS_EXPECT_TRUE((Filter::Screen::DirectAntiAliasSink<W, H>::is_pipeline));
  HS_EXPECT_TRUE((Pipeline<W, H>::is_pipeline));
  HS_EXPECT_TRUE(
      (Pipeline<W, H, Filter::Screen::AntiAlias<W, H>>::is_pipeline));

  // Static-assert form (compile-time).
  static_assert(Filter::Screen::AntiAlias<W, H>::is_2d, "AntiAlias is 2D");
  static_assert(!Filter::World::Replicate<W>::is_2d, "Replicate is 3D");
  static_assert(Filter::Pixel::Feedback<W, H>::has_history,
                "Feedback keeps history");
  static_assert(Filter::has_cull_edge<Filter::World::Orient>);
  static_assert(Filter::has_cull_edge<Filter::World::Replicate<W>>);
  using Ordered = Pipeline<W, H, Filter::World::Replicate<W>,
                           Filter::Screen::AntiAlias<W, H>,
                           Filter::Pixel::ChromaticShift<W>>;
  static_assert(sizeof(Ordered) > 0,
                "World/Screen/Pixel domain order must compile");
}

/**
 * @brief Verifies the `crosses_segments` per-filter trait and the Pipeline
 *        `any_crosses_segments` OR-fold the segment driver gates on.
 */
inline void test_crosses_segments_trait_and_fold() {
  constexpr int W = 32, H = 16;

  // Per-filter trait. The fail-safe default ties crosses_segments to has_history.
  HS_EXPECT_FALSE((Filter::Screen::AntiAlias<W, H>::crosses_segments));
  HS_EXPECT_FALSE((Filter::World::Replicate<W>::crosses_segments));
  HS_EXPECT_TRUE((Filter::Pixel::Feedback<W, H>::crosses_segments));
  HS_EXPECT_TRUE((Filter::World::Trails<16>::crosses_segments));
  HS_EXPECT_TRUE((Filter::Screen::Trails<>::crosses_segments));
  HS_EXPECT_FALSE((Filter::Screen::AntiAlias<W, H>::reads_outside_band));
  HS_EXPECT_TRUE((Filter::Pixel::Feedback<W, H>::reads_outside_band));
  HS_EXPECT_FALSE((Filter::World::Trails<16>::reads_outside_band));
  HS_EXPECT_FALSE((Filter::Screen::Trails<>::reads_outside_band));

  // Pipeline OR-fold.
  HS_EXPECT_FALSE((Pipeline<W, H>::any_crosses_segments));
  HS_EXPECT_FALSE((Pipeline<W, H>::any_reads_outside_band));

  using MeshStack =
      Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>,
               Filter::Pixel::Feedback<W, H>>;
  HS_EXPECT_TRUE(MeshStack::any_crosses_segments);
  HS_EXPECT_TRUE(MeshStack::any_reads_outside_band);

  // A non-stateful stack does not.
  using PlainStack =
      Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>>;
  HS_EXPECT_FALSE(PlainStack::any_crosses_segments);
  HS_EXPECT_FALSE(PlainStack::any_reads_outside_band);

  HS_EXPECT_TRUE(
      (Pipeline<W, H, Filter::Screen::Trails<>>::any_crosses_segments));
  HS_EXPECT_FALSE(
      (Pipeline<W, H, Filter::World::Trails<16>>::any_reads_outside_band));

  // The pipeline spelling exposes the fold. MeshStack's head (World::Orient)
  // is false while the composed pipeline is true.
  HS_EXPECT_TRUE(MeshStack::crosses_segments);
  HS_EXPECT_FALSE(PlainStack::crosses_segments);

  // Compile-time form.
  static_assert(MeshStack::any_crosses_segments,
                "MeshFeedback pipeline must render full-frame per worker");
  static_assert(
      MeshStack::crosses_segments == MeshStack::any_crosses_segments,
      "Pipeline::crosses_segments must answer the fold, not the head");
  static_assert(!PlainStack::any_crosses_segments,
                "non-stateful pipeline must keep the segment clipping win");

  // segment_margin: the render-bound padding a stage's off-position taps need.
  HS_EXPECT_EQ((Filter::Pixel::ChromaticShift<W>::segment_margin), 3);
  HS_EXPECT_EQ((Filter::Screen::AntiAlias<W, H>::segment_margin), 1);
  HS_EXPECT_EQ((Filter::Screen::Blur<W, H>::segment_margin), 1);
  HS_EXPECT_EQ((Filter::Screen::DirectAntiAliasSink<W, H>::segment_margin), 1);
  HS_EXPECT_EQ((Filter::World::Orient::segment_margin), 0);
  HS_EXPECT_EQ((Filter::Pixel::Feedback<W, H>::segment_margin), 0);

  // total_segment_margin fold: empty pipeline, plain stack, a lone
  // ChromaticShift, and a chain of spreading stages, whose margins sum —
  // AntiAlias splats a tap 1 column off-position, then ChromaticShift shifts
  // that tap up to 3 more.
  HS_EXPECT_EQ((Pipeline<W, H>::total_segment_margin), 0);
  HS_EXPECT_EQ((PlainStack::total_segment_margin), 1);
  HS_EXPECT_EQ((MeshStack::total_segment_margin), 1);
  HS_EXPECT_EQ(
      (Pipeline<W, H, Filter::Pixel::ChromaticShift<W>>::total_segment_margin),
      3);
  using ShiftStack =
      Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>,
               Filter::Pixel::ChromaticShift<W>>;
  HS_EXPECT_EQ((ShiftStack::total_segment_margin), 4);
  HS_EXPECT_EQ((ShiftStack::segment_margin), ShiftStack::total_segment_margin);
  // ChromaticShift pads the render bounds instead of forcing a full frame.
  static_assert(!Filter::Pixel::ChromaticShift<W>::crosses_segments,
                "ChromaticShift must stay band-clippable");
  static_assert(!ShiftStack::any_crosses_segments,
                "a ChromaticShift stack must keep the segment clipping win");

  // Ordering traits: a misspelled override silently inherits the FilterTraits
  // default, so pin both the stage value and the pipeline fold.
  HS_EXPECT_TRUE((Filter::Screen::AntiAlias<W, H>::emits_pixel_centers));
  HS_EXPECT_TRUE((Filter::Screen::Blur<W, H>::emits_pixel_centers));
  HS_EXPECT_FALSE((Filter::Pixel::ChromaticShift<W>::emits_pixel_centers));
  HS_EXPECT_TRUE((Filter::Screen::AntiAlias<W, H>::requires_subpixel_input));
  HS_EXPECT_FALSE((Filter::Screen::Blur<W, H>::requires_subpixel_input));
  // Trails re-emits whatever coordinates it was handed: rounded taps still
  // seed and fade, so a rounding stage may precede it.
  HS_EXPECT_FALSE((Filter::Screen::Trails<>::requires_subpixel_input));
  HS_EXPECT_FALSE((Filter::Screen::Trails<>::emits_pixel_centers));

  HS_EXPECT_TRUE((Filter::World::Trails<16>::emits_nonunit_world));
  HS_EXPECT_FALSE((Filter::World::Orient::emits_nonunit_world));
  HS_EXPECT_TRUE((Filter::World::Mobius::requires_unit_world_input));
  HS_EXPECT_TRUE((Filter::World::Hole::requires_unit_world_input));
  HS_EXPECT_TRUE((Filter::World::OrientSlice::requires_unit_world_input));
  HS_EXPECT_FALSE((Filter::World::Trails<16>::requires_unit_world_input));

  HS_EXPECT_TRUE((Filter::World::Hole::world_transform_is_identity));
  HS_EXPECT_FALSE((Filter::World::Mobius::world_transform_is_identity));
  HS_EXPECT_FALSE((Filter::World::Orient::world_transform_is_identity));
  HS_EXPECT_TRUE(
      (Filter::Screen::AntiAlias<W, H>::world_transform_is_identity));

  HS_EXPECT_TRUE((Filter::Pixel::Feedback<W, H>::terminal_replaces));
  HS_EXPECT_FALSE((Filter::Screen::Trails<>::terminal_replaces));

  using WarpStack =
      Pipeline<W, H, Filter::World::Trails<16>, Filter::Screen::Blur<W, H>>;
  HS_EXPECT_TRUE(WarpStack::emits_nonunit_world);
  HS_EXPECT_FALSE(WarpStack::requires_unit_world_input);
  HS_EXPECT_TRUE(WarpStack::emits_pixel_centers);
  HS_EXPECT_FALSE(WarpStack::requires_subpixel_input);
  HS_EXPECT_TRUE(WarpStack::world_transform_is_identity);
  HS_EXPECT_FALSE(WarpStack::terminal_replaces);

  HS_EXPECT_TRUE(MeshStack::emits_pixel_centers);
  HS_EXPECT_TRUE(MeshStack::requires_subpixel_input);
  HS_EXPECT_TRUE(MeshStack::terminal_replaces);
  HS_EXPECT_FALSE(MeshStack::emits_nonunit_world);
  HS_EXPECT_FALSE(MeshStack::world_transform_is_identity);
}

/**
 * @brief Verifies the `any_2d_history` / `any_3d_history` folds that gate the
 *        flush() overloads against a wrong-domain (silently empty) call.
 * @details The rejection is a static_assert in the overload body, so this
 *          pins the fold each assert reads.
 */
inline void test_history_domain_folds() {
  constexpr int W = 32, H = 16;

  HS_EXPECT_FALSE((Pipeline<W, H>::any_2d_history));
  HS_EXPECT_FALSE((Pipeline<W, H>::any_3d_history));

  using ScreenStack =
      Pipeline<W, H, Filter::World::Orient, Filter::Screen::Trails<>>;
  HS_EXPECT_TRUE(ScreenStack::any_2d_history);
  HS_EXPECT_FALSE(ScreenStack::any_3d_history);

  using WorldStack = Pipeline<W, H, Filter::World::Trails<16>,
                              Filter::Screen::AntiAlias<W, H>>;
  HS_EXPECT_FALSE(WorldStack::any_2d_history);
  HS_EXPECT_TRUE(WorldStack::any_3d_history);

  using MixedStack =
      Pipeline<W, H, Filter::World::Trails<16>, Filter::Screen::Trails<>>;
  HS_EXPECT_TRUE(MixedStack::any_2d_history);
  HS_EXPECT_TRUE(MixedStack::any_3d_history);

  // A history-free stack answers neither, so either overload is a hard error.
  using PlainStack =
      Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>>;
  HS_EXPECT_FALSE(PlainStack::any_2d_history);
  HS_EXPECT_FALSE(PlainStack::any_3d_history);

  // Feedback is a 2D-history terminal.
  HS_EXPECT_TRUE(
      (Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>,
                Filter::Pixel::Feedback<W, H>>::any_2d_history));

  // The direct sink publishes the same folds as a Pipeline.
  using Direct = Filter::Screen::DirectAntiAliasSink<W, H>;
  HS_EXPECT_FALSE(Direct::any_2d_history);
  HS_EXPECT_FALSE(Direct::any_3d_history);
}
