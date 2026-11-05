/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Authored palette source and cycling harmony contracts. */

struct AuthoredPaletteSpec : Pullback::Spec {
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
          Pullback::ProjectionKind::GNOMONIC_FOLDED, B>::Type>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
          Pullback::SquareCellsParams, B>::Type>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::AuthoredColorProvider<B>>>>;
};

struct AuthoredPaletteProbe
    : Pullback::ComposedEffect<SMALL_W, SMALL_H, AuthoredPaletteProbe,
                               AuthoredPaletteSpec> {
  using Base = Pullback::ComposedEffect<SMALL_W, SMALL_H, AuthoredPaletteProbe,
                                        AuthoredPaletteSpec>;
  static constexpr Params initial_params() { return {}; }
  static consteval auto choreography_hooks() {
    return Base::choreography_hooks()
        .template extend<
            AuthoredPaletteProbe,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &AuthoredPaletteProbe::initial_params>>();
  }
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "authored-palette"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
};

/** @brief Fixed color providers bake every authored source as its direct recipe. */
inline void test_composed_authored_palette_sources() {
  reset_effect_globals();
  AuthoredPaletteProbe effect;
  effect.init();
  AuthoredPaletteProbe::Params params{};
  auto &color = params.template get<"color">();
  BakedPaletteStorage expected;
  bool first = true;
  const auto check = [&](const auto &source) {
    if (first) {
      expected.bake(persistent_arena, source);
      first = false;
    } else
      expected.rebake(source);
    ComposedFrameWhiteBox::set_params(effect, params);
    ComposedFrameWhiteBox::advance(effect);
    const auto frame = effect.frame_for_test();
    for (float t : {0.0f, .13f, .37f, .71f, 1.0f}) {
      const auto actual = frame.palette->get(t);
      const auto reference = expected.view().get(t);
      HS_EXPECT_EQ(actual.color.r, reference.color.r);
      HS_EXPECT_EQ(actual.color.g, reference.color.g);
      HS_EXPECT_EQ(actual.color.b, reference.color.b);
      HS_EXPECT_EQ(actual.alpha, reference.alpha);
    }
  };
  color.palette_type = static_cast<uint8_t>(AuthoredPaletteType::GENERATED);
  color.palette_hue = .23f;
  color.palette_spread = .11f;
  color.palette_lightness = .58f;
  color.palette_lightness_range = .34f;
  for (uint8_t harmony = 0;
       harmony <= static_cast<uint8_t>(PaletteHarmony::SQUARE); ++harmony)
    for (uint8_t curve = 0; curve <= static_cast<uint8_t>(AxisCurve::CUP);
         ++curve) {
      color.palette_harmony = harmony;
      color.palette_lightness_curve = curve;
      PaletteRecipe recipe;
      recipe.hue.harmony = static_cast<PaletteHarmony>(harmony);
      recipe.hue.base_turns = .23f;
      recipe.hue.spread_turns = .11f;
      recipe.lightness.curve = static_cast<AxisCurve>(curve);
      recipe.lightness.center = .58f;
      recipe.lightness.range = .34f;
      recipe.chroma.axis.center = .62f;
      check(GenerativePalette(recipe));
    }
#define HS_AUTHORED_PALETTE_REFERENCE(name, A, B, C, D) &Palettes::name,
  const Palette *const named[] = {
      HS_PROCEDURAL_PALETTE_LIST(HS_AUTHORED_PALETTE_REFERENCE)};
#undef HS_AUTHORED_PALETTE_REFERENCE
  color.palette_type = static_cast<uint8_t>(AuthoredPaletteType::PROCEDURAL);
  for (size_t index = 0; index < std::size(named); ++index) {
    color.palette_name = static_cast<uint8_t>(index);
    check(*named[index]);
  }
  color.palette_type = static_cast<uint8_t>(AuthoredPaletteType::GRADIENT);
  color.gradient_start_r = 1.0f;
  color.gradient_start_g = 0.0f;
  color.gradient_start_b = 0.0f;
  color.gradient_middle_r = 0.0f;
  color.gradient_middle_g = 1.0f;
  color.gradient_middle_b = 0.0f;
  color.gradient_end_r = 0.0f;
  color.gradient_end_g = 0.0f;
  color.gradient_end_b = 1.0f;
  check(Gradient{{0.0f, CPixel(255, 0, 0)},
                 {.5f, CPixel(0, 255, 0)},
                 {1.0f, CPixel(0, 0, 255)}});
  color.palette_type = static_cast<uint8_t>(AuthoredPaletteType::MUTATING);
  color.palette_name = 0;
  color.palette_name_b = 19;
  MutatingPalette mutation({.367f, .367f, .367f}, {.5f, .5f, .5f}, {1, 1, 1},
                           {0, .33f, .67f}, {0.500f, 0.500f, 0.500f},
                           {0.500f, 0.500f, 0.500f}, {0.955f, 1.004f, 0.910f},
                           {0.167f, 0.018f, 0.930f});
  for (float amount : {0.0f, .4f, 1.0f}) {
    color.palette_mutation = amount;
    mutation.mutate(amount);
    check(mutation);
  }
}

/** @brief Switching cycling harmonies preserves their elapsed fade time. */
inline void test_composed_authored_palette_harmony_clock() {
  reset_effect_globals();
  AuthoredPaletteProbe effect;
  effect.init();
  GeneratedPaletteBank reference;
  reference.init(persistent_arena, .62f, math::ease_in_out_sin);
  using Mode = Pullback::Interp::Op::PaletteMode;
  for (int step = 0; step < 237; ++step)
    reference.step(Mode::TRIADIC);
  ComposedFrameWhiteBox::advance_palette(effect, 237);
  AuthoredPaletteProbe::Params params{};
  for (uint8_t mode : {uint8_t(0), uint8_t(1), uint8_t(2), uint8_t(0)}) {
    params.template get<"color">().palette_mode = mode;
    ComposedFrameWhiteBox::set_params(effect, params);
    ComposedFrameWhiteBox::advance(effect);
    const auto frame = effect.frame_for_test();
    const auto snapshot = reference.snapshot();
    reference.restore_snapshot(snapshot);
    for (float t : {.07f, .31f, .83f}) {
      const auto actual = frame.palette->get(t).color;
      const auto expected =
          reference.palette(static_cast<Mode>(mode)).get(t).color;
      HS_EXPECT_EQ(actual.r, expected.r);
      HS_EXPECT_EQ(actual.g, expected.g);
      HS_EXPECT_EQ(actual.b, expected.b);
    }
  }
  auto invalid = params;
  invalid.template get<"color">().palette_type = 255;
  HS_EXPECT_FALSE(AuthoredPaletteProbe::valid_params(invalid));
  invalid = params;
  invalid.template get<"color">().brightness_envelope =
      static_cast<Pullback::Color::BrightnessEnvelope>(255);
  HS_EXPECT_FALSE(AuthoredPaletteProbe::valid_params(invalid));
}
