/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Authored palette source and cycling harmony contracts. */

#include "effects/RiverStones.h"
#include "tests/composed_effect/context.h"

namespace hs_test::composed_effect_tests {

/** @brief Palette-source pipeline for authored palette controls. */
struct AuthoredPaletteSpec : Pullback::Spec {
  /** @brief Typed stage sequence evaluated by the probe.
   * @tparam B Pipeline binding.
   */
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

/** @brief Composed effect exposing authored palette controls. */
struct AuthoredPaletteProbe
    : Pullback::ComposedEffect<SMALL_W, SMALL_H, AuthoredPaletteProbe,
                               AuthoredPaletteSpec> {
  /** @brief Composed-effect base supplying the runtime resources. */
  using Base = Pullback::ComposedEffect<SMALL_W, SMALL_H, AuthoredPaletteProbe,
                                        AuthoredPaletteSpec>;
  /** @brief Builds the initial controls of the probe.
   * @return Initial parameter block.
   */
  static constexpr Params initial_params() { return {}; }
  /** @brief Binds the probe choreography customizations.
   * @return Type-only list of choreography bindings.
   */
  static consteval auto choreography_hooks() {
    return Base::choreography_hooks()
        .template extend<
            AuthoredPaletteProbe,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &AuthoredPaletteProbe::initial_params>>();
  }
  /** @brief Stable identifiers exported for the probe presets. */
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "authored-palette"};
  /** @brief Version of the serialized parameter schema. */
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  /** @brief Frames spent on each preset before automatic selection. */
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
};

/** @brief Fixed color providers bake every authored source as its direct recipe. */
inline void test_composed_authored_palette_sources() {
  reset_globals();
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
#define HS_AUTHORED_PALETTE_COEFFICIENTS(name, A, B, C, D)                     \
  std::array<std::array<float, 3>, 4>{{{HS_PALETTE_VEC3 A},                    \
                                       {HS_PALETTE_VEC3 B},                    \
                                       {HS_PALETTE_VEC3 C},                    \
                                       {HS_PALETTE_VEC3 D}}},
  const auto coefficients =
      std::array{HS_PROCEDURAL_PALETTE_LIST(HS_AUTHORED_PALETTE_COEFFICIENTS)};
#undef HS_AUTHORED_PALETTE_COEFFICIENTS
  const auto &palette_a = coefficients[color.palette_name];
  const auto &palette_b = coefficients[color.palette_name_b];
  MutatingPalette mutation(palette_a[0], palette_a[1], palette_a[2],
                           palette_a[3], palette_b[0], palette_b[1],
                           palette_b[2], palette_b[3]);
  for (float amount : {0.0f, .4f, 1.0f}) {
    color.palette_mutation = amount;
    mutation.mutate(amount);
    check(mutation);
  }
}

/** @brief Switching cycling harmonies preserves their elapsed fade time. */
inline void test_composed_authored_palette_harmony_clock() {
  reset_globals();
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
}

#if HS_ENABLE_PARAM_GUI_BRIDGE
/** @brief Cycling harmony changes crossfade palettes on the same clock. */
inline void test_composed_cycling_harmony_crossfade() {
  using Traits = Pullback::ComposedTraits<SMALL_W, AuthoredPaletteSpec>;
  using Palette = Pullback::ComposedPalette<Traits>;
  reset_globals();
  typename Traits::ColorFamily from, to;
  to.palette_mode = 1;
  auto live = from;
  Palette outgoing(from), incoming(to), transitioning(live);
  struct Storage {
    AuthoredPaletteResources authored_palette;
  } storage;
  outgoing.init();
  incoming.init();
  transitioning.init();
  const auto step = [&] {
    outgoing.step();
    incoming.step();
    transitioning.step();
  };
  for (int frame = 0; frame < 597; ++frame)
    step();

  transitioning.arm(to);
  transitioning.arm_source(from, to);
  live = to;
  transitioning.blend(0.5f);
  BakedPaletteStorage expected;
  typename Traits::FrameState frame{};
  for (int tick = 0; tick < 8; ++tick) {
    step();
    if (tick == 0)
      expected.bake_blend(persistent_arena, outgoing.selected_palette(storage),
                          incoming.selected_palette(storage), 0.5f);
    else
      expected.rebake_crossfade(outgoing.selected_palette(storage),
                                incoming.selected_palette(storage), 0.5f);
    transitioning.prepare(frame, storage);
    HS_EXPECT_TRUE(frame.palette->same_entries(expected.view()));
    HS_EXPECT_FALSE(
        frame.palette->same_entries(outgoing.selected_palette(storage)));
    HS_EXPECT_FALSE(
        frame.palette->same_entries(incoming.selected_palette(storage)));
  }

  transitioning.adopt(to);
  step();
  transitioning.prepare(frame, storage);
  HS_EXPECT_TRUE(
      frame.palette->same_entries(incoming.selected_palette(storage)));
}

/** @brief RiverStones blends a manual palette selection into the next preset. */
inline void test_riverstones_manual_palette_crossfade() {
  using FX = RiverStones<SMALL_W, SMALL_H>;
  const int dwell = FX::PRESET_DWELL_FRAMES;
  const int duration = Segue::Preset::frames(FX::preset_departure(0));
  for (const auto type :
       {AuthoredPaletteType::PROCEDURAL, AuthoredPaletteType::GRADIENT}) {
    reset_globals();
    FX effect;
    effect.init();
    HS_EXPECT_EQ(effect.updateParameter("Palette Type", uint8_t(type)),
                 ParamSetResult::APPLIED);
    HS_EXPECT_EQ(effect.updateParameter("Palette Name", 23.0f),
                 ParamSetResult::APPLIED);
    effect.setAnimationsPaused(false);
    const auto start = effect.frame_for_test();
    const auto target = FX::test_preset_params(1);
    const auto &color = target.template get<"color">();
    BakedPaletteStorage from;
    from.clone_from(*start.palette, persistent_arena);
    AuthoredPaletteResources reference;
    const auto &to =
        reference.prepare(persistent_arena, color, color.palette_chroma);
    HS_EXPECT_FALSE(from.view().same_entries(to));

    int frame = 0;
    const auto draw_to = [&](int end) {
      while (frame < end) {
        pin_frame_clock(++frame);
        effect.draw_frame();
        effect.advance_display();
      }
      return ComposedFrameWhiteBox::current_frame(effect).palette;
    };
    HS_EXPECT_TRUE(draw_to(dwell)->same_entries(from.view()));
    HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
    BakedPaletteStorage expected;
    const auto before = persistent_arena.get_offset();
    const auto *middle = draw_to(dwell + duration / 2);
    const auto after = persistent_arena.get_offset();
    expected.bake_blend(persistent_arena, from.view(), to, 0.5f);
    HS_EXPECT_TRUE(middle->same_entries(expected.view()));
    HS_EXPECT_FALSE(middle->same_entries(from.view()));
    HS_EXPECT_FALSE(middle->same_entries(to));
    HS_EXPECT_TRUE(after > before);
    const auto stable = persistent_arena.get_offset();
    HS_EXPECT_TRUE(draw_to(dwell + duration)->same_entries(to));
    HS_EXPECT_EQ(persistent_arena.get_offset(), stable);
  }

  reset_globals();
  FX effect;
  effect.init();
  pin_frame_clock(1);
  effect.draw_frame();
  effect.advance_display();
  const auto before = persistent_arena.get_offset();
  for (int frame = 2; frame <= dwell + duration; ++frame) {
    pin_frame_clock(frame);
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
  HS_EXPECT_EQ(persistent_arena.get_offset(), before);
}
#endif

/** @brief Authored palettes admit each selector's last ordinal and reject the
 *  next one. */
inline void test_composed_authored_palette_admission() {
  using Color = Pullback::Color::AuthoredColorParams;
  struct Selector {
    const char *name;
    uint8_t first_invalid;
    void (*set)(Color &, uint8_t);
  };
  const Selector selectors[] = {
      {"palette_type", uint8_t(AuthoredPaletteType::MUTATING) + 1,
       [](Color &c, uint8_t v) { c.palette_type = v; }},
      {"palette_mode", uint8_t(std::size(Pullback::Color::PALETTE_MODE_IDS)),
       [](Color &c, uint8_t v) { c.palette_mode = v; }},
      {"palette_harmony", uint8_t(PaletteHarmony::SQUARE) + 1,
       [](Color &c, uint8_t v) { c.palette_harmony = v; }},
      {"palette_lightness_curve", uint8_t(AxisCurve::CUP) + 1,
       [](Color &c, uint8_t v) { c.palette_lightness_curve = v; }},
      {"palette_name", uint8_t(std::size(AUTHORED_PALETTE_NAME_IDS)),
       [](Color &c, uint8_t v) { c.palette_name = v; }},
      {"palette_name_b", uint8_t(std::size(AUTHORED_PALETTE_NAME_IDS)),
       [](Color &c, uint8_t v) { c.palette_name_b = v; }},
      {"brightness_envelope",
       uint8_t(Pullback::Color::BrightnessEnvelope::DESCENDING) + 1,
       [](Color &c, uint8_t v) {
         c.brightness_envelope =
             static_cast<Pullback::Color::BrightnessEnvelope>(v);
       }},
  };
  for (const Selector &selector : selectors) {
    HS_CONTEXT(selector.name);
    for (uint8_t value :
         {uint8_t(selector.first_invalid - 1), selector.first_invalid}) {
      AuthoredPaletteProbe::Params params{};
      selector.set(params.template get<"color">(), value);
      HS_EXPECT_EQ(AuthoredPaletteProbe::valid_params(params),
                   value < selector.first_invalid);
    }
  }
}

} // namespace hs_test::composed_effect_tests
