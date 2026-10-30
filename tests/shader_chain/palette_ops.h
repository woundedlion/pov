/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Authored palette families drive rendered colors and survive snapshots. */
inline void test_shader_chain_palette_selection() {
  using WB = ShaderChainWhiteBox;
  using Type = In::Op::PaletteType;
  reset_globals();
  WB::FX effect;
  effect.init();
  auto &params = WB::color_params(effect);
  params.palette_type = static_cast<uint8_t>(Type::GENERATED);
  params.palette_hue = .17f;
  params.palette_spread = .11f;
  params.palette_lightness = .53f;
  params.palette_lightness_range = .42f;
  params.palette_lightness_curve = static_cast<uint8_t>(AxisCurve::ASCENDING);
  constexpr float COORDINATE = 73.0f / 255.0f;
  const PB::FieldSample SAMPLE{COORDINATE, .75f, math::X_AXIS, 0};
  for (uint8_t harmony = 0;
       harmony <= static_cast<uint8_t>(PaletteHarmony::SQUARE); ++harmony) {
    params.palette_harmony = harmony;
    PaletteRecipe recipe;
    recipe.hue.harmony = static_cast<PaletteHarmony>(harmony);
    recipe.hue.base_turns = .17f;
    recipe.hue.spread_turns = .11f;
    recipe.lightness.curve = AxisCurve::ASCENDING;
    recipe.lightness.center = .53f;
    recipe.lightness.range = .42f;
    recipe.chroma.axis.center = params.palette_chroma;
    const GenerativePalette expected(recipe);
    const auto ctx = WB::frame_context(effect);
    WB::program(effect).prepare(ctx);
    const Color4 actual =
        run_entry<Color4>(WB::program(effect), ctx, 3, SAMPLE);
    HS_EXPECT_EQ(actual.color, expected.get(COORDINATE).color);
    HS_EXPECT_EQ(actual.alpha, SAMPLE.coverage);
  }
  params.palette_type = static_cast<uint8_t>(Type::PROCEDURAL);
#define HS_TEST_NAMED_PALETTE(name, A, B, C, D) &Palettes::name,
  const Palette *SOURCES[] = {
      HS_PROCEDURAL_PALETTE_LIST(HS_TEST_NAMED_PALETTE)};
#undef HS_TEST_NAMED_PALETTE
  HS_EXPECT_EQ(std::size(SOURCES), std::size(In::Op::PALETTE_NAME_IDS));
  for (uint8_t index = 0; index < std::size(SOURCES); ++index) {
    params.palette_name = index;
    const auto named_ctx = WB::frame_context(effect);
    WB::program(effect).prepare(named_ctx);
    HS_EXPECT_EQ(
        run_entry<Color4>(WB::program(effect), named_ctx, 3, SAMPLE).color,
        SOURCES[index]->get(COORDINATE).color);
  }
  params.palette_name = 19;
  auto ctx = WB::frame_context(effect);
  WB::program(effect).prepare(ctx);
  HS_EXPECT_EQ(run_entry<Color4>(WB::program(effect), ctx, 3, SAMPLE).color,
               Palettes::FIRE_AND_ICE.get(COORDINATE).color);
  params.palette_name = 0;
  ctx = WB::frame_context(effect);
  HS_EXPECT_EQ(ctx.palettes[0]->get(COORDINATE).color,
               Palettes::DARK_RAINBOW.get(COORDINATE).color);
  params.palette_type = static_cast<uint8_t>(Type::MUTATING);
  params.palette_mutation = .5f;
  ctx = WB::frame_context(effect);
  const ProceduralPalette MIXED({.4335f, .4335f, .4335f}, {.5f, .5f, .5f},
                                {.9775f, 1.002f, .955f}, {.0835f, .174f, .8f});
  HS_EXPECT_LE(max_channel_gap(ctx.palettes[0]->get(COORDINATE).color,
                               MIXED.get(COORDINATE).color),
               2);
  params.palette_mutation = 1;
  ctx = WB::frame_context(effect);
  HS_EXPECT_EQ(ctx.palettes[0]->get(COORDINATE).color,
               Palettes::FIRE_AND_ICE.get(COORDINATE).color);
  params.palette_type = static_cast<uint8_t>(Type::GRADIENT);
  params.gradient_start_r = .2f;
  params.gradient_start_g = .4f;
  params.gradient_start_b = .8f;
  params.gradient_middle_r = .8f;
  params.gradient_middle_g = .2f;
  params.gradient_middle_b = .4f;
  params.gradient_end_r = .4f;
  params.gradient_end_g = .8f;
  params.gradient_end_b = .2f;
  const Gradient GRADIENT{{0, CPixel(51, 102, 204)},
                          {.5f, CPixel(204, 51, 102)},
                          {1, CPixel(102, 204, 51)}};
  ctx = WB::frame_context(effect);
  HS_EXPECT_EQ(ctx.palettes[0]->get(COORDINATE).color,
               GRADIENT.get(COORDINATE).color);
  const Pixel COLOR = ctx.palettes[0]->get(COORDINATE).color;
  params.palette_type = static_cast<uint8_t>(Type::CYCLING_GENERATED);
  ctx = WB::frame_context(effect);
  HS_EXPECT_NE(ctx.palettes[0]->get(COORDINATE).color, COLOR);
  params.palette_type = static_cast<uint8_t>(Type::GRADIENT);
  ctx = WB::frame_context(effect);
  HS_EXPECT_EQ(ctx.palettes[0]->get(COORDINATE).color, COLOR);
  const ChainSnapshot saved = effect.snapshot();
  params.palette_chroma = .05f;
  ctx = WB::frame_context(effect);
  HS_EXPECT_NE(ctx.palettes[0]->get(COORDINATE).color, COLOR);
  HS_EXPECT_EQ(effect.restore_snapshot(saved),
               ChainSnapshotRestoreResult::APPLIED);
  ctx = WB::frame_context(effect);
  HS_EXPECT_EQ(ctx.palettes[0]->get(COORDINATE).color, COLOR);
  HS_EXPECT_EQ(WB::color_params(effect).palette_type,
               static_cast<uint8_t>(Type::GRADIENT));
}

/** @brief Authored palette rebakes preserve the retained temporal color endpoint. */
inline void test_shader_chain_authored_palette_history() {
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX effect;
  effect.init();
  auto &params = WB::color_params(effect);
  params.palette_type = static_cast<uint8_t>(In::Op::PaletteType::GENERATED);
  pin_frame_clock(0);
  effect.draw_frame();
  effect.advance_display();
  const Pixel BEFORE =
      WB::previous_context(effect).palettes[0]->get(.37f).color;
  params.palette_hue = .5f;
  const auto ctx = WB::frame_context(effect);
  HS_EXPECT_NE(ctx.palettes[0]->get(.37f).color, BEFORE);
  HS_EXPECT_EQ(WB::previous_context(effect).palettes[0]->get(.37f).color,
               BEFORE);
  WB::program(effect).prepare_column(ctx, WB::previous_context(effect), .5f);
  const auto &prepared =
      *reinterpret_cast<const PB::Color::GeneratedPaletteState *>(
          WB::program(effect).prepared_block(3));
  HS_EXPECT_EQ(prepared.palette, ctx.palettes[0]);
  hs::clear_mock_time();
}
