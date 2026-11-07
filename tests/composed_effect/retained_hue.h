/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @brief Field carrier for retained color-tail evaluation. */
struct RetainedHueField
    : Pullback::Stage::Contract<RetainedHueField, Pullback::SphereSample,
                                Pullback::FieldSample> {
  using Policies = std::tuple<>;
  template <typename Binding>
  static Pullback::FieldSample run(const Pullback::SphereSample &input,
                                   const typename Binding::FrameState &,
                                   const Pullback::NoPrepared &) {
    return {0.5f * (input.dir.x + 1.0f), 0.75f, input.dir, 0.0f};
  }
};

/** @brief Owned palette history and direct endpoints for hue retention checks. */
struct RetainedHueFixture {
  using FX = GridWorms<SMALL_W, SMALL_H>;
  using Palette = Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
      typename FX::Binding, Pullback::Color::BrightnessEnvelope::NONE>>;
  using Provider =
      Pullback::ModulationProvider<typename FX::Binding, "modulation",
                                   Pattern::SourceKind::DIRECTIONAL_SINE,
                                   Pattern::Evaluation::DIRECT>;
  using Layer = Pullback::Color::Modulate<Provider>;
  using Tail = Pullback::Pipeline<typename FX::Binding, RetainedHueField,
                                  Pullback::Stage::Colorize<Palette>,
                                  Pullback::Stage::Modulate<Layer>>;
  FX effect = make_effect();
  FX::FrameState current{}, previous{};
  Spacetime::FrameHistory<FX::FrameState, true, true> history;
  Pullback::Modulation::Config before{}, after{};
  Pullback::Modulation::Resources<false> old_source, new_source;

  static FX make_effect() {
    reset_globals();
    return FX{};
  }
  RetainedHueFixture() {
    using namespace Pullback;
    effect.init();
    effect.draw_frame();
    effect.advance_display();
    current = effect.frame_for_test();
    history.capture(current, persistent_arena);
    previous = history.previous;
    HS_EXPECT_TRUE(previous.palette != current.palette);
    HS_EXPECT_TRUE(previous.hue_rotation_lut != current.hue_rotation_lut);
    HS_EXPECT_TRUE(previous.hue_rotation_gamut != current.hue_rotation_gamut);
    HS_EXPECT_TRUE(previous.palette->same_entries(*current.palette));
    for (auto *frame : {&previous, &current}) {
      auto &color = frame->params.template get<"color">();
      color.mapping_frequency = 1.0f;
      color.mapping_phase = 0.0f;
      color.phase_oscillation_depth = 0.0f;
      color.opacity_low = color.opacity_high = 0.37f;
      frame->palette_mapping =
          Color::PaletteMappingWeights::single(Color::PaletteMapping::LINEAR);
    }
    before = previous.params.template get<"modulation">().modulation;
    before.enabled = true;
    before.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 1.0f};
    before.sine.axis = math::Z_AXIS;
    before.sine.frequency = 2.0f;
    before.sine.offset = 0.17f;
    before.sine.speed = 0.0f;
    after = before;
    after.sine.frequency = 3.25f;
    after.sine.offset = 0.31f;
    auto clock = Modulation::init(before);
    clock.phase.phase_turns = 0.125f;
    old_source.prepare(before, clock);
    clock.phase.phase_turns = 0.375f;
    new_source.prepare(after, clock);
    previous.params.template get<"modulation">().modulation = before;
    current.params.template get<"modulation">().modulation = after;
    previous.resources.template get<"modulation">().modulation =
        &old_source.prepared;
    current.resources.template get<"modulation">().modulation =
        &new_source.prepared;
  }
};

/** @brief Held palettes retain the current hue view through recipe fallback. */
inline void test_retained_equal_palette_hue_view() {
  using namespace Pullback;
  RetainedHueFixture fixture;
  using FX = RetainedHueFixture::FX;
  using Palette = RetainedHueFixture::Palette;
  using Layer = RetainedHueFixture::Layer;
  using Tail = RetainedHueFixture::Tail;
  auto &current = fixture.current;
  auto &previous = fixture.previous;
  const auto &before = fixture.before;
  const FrameInterval INTERVAL(previous, current);
  const auto START = INTERVAL.at(0.0f);
  const auto END = INTERVAL.at(1.0f);
  HS_EXPECT_TRUE(START.hue_rotation_lut == previous.hue_rotation_lut);
  HS_EXPECT_TRUE(START.temporal.hue_rotation == previous.hue_rotation_lut);
  HS_EXPECT_TRUE(START.temporal.hue_rotation_gamut ==
                 previous.hue_rotation_gamut);
  HS_EXPECT_TRUE(END.hue_rotation_lut == current.hue_rotation_lut);
  HS_EXPECT_TRUE(END.hue_rotation_gamut == current.hue_rotation_gamut);
  uint32_t maximum = 0;
  uint64_t total = 0, channels = 0;
  size_t different = 0;
  for (const float FRACTION : {0.17f, 0.5f, 0.83f}) {
    const auto COLUMN = INTERVAL.at(FRACTION);
    HS_EXPECT_TRUE(COLUMN.temporal.same_palette);
    HS_EXPECT_FALSE(COLUMN.temporal.phase_compatible);
    HS_EXPECT_TRUE(COLUMN.temporal.hue_rotation == current.hue_rotation_lut);
    HS_EXPECT_TRUE(COLUMN.temporal.hue_rotation_gamut ==
                   current.hue_rotation_gamut);
    HS_EXPECT_NEAR(COLUMN.temporal.fraction, FRACTION, 0.0f);
    const auto PALETTE = Palette::prepare(COLUMN);
    const auto MODULATION = Layer::prepare(COLUMN);
    HS_EXPECT_FALSE(MODULATION.ready);
    const auto TYPED = Tail::prepare(COLUMN);
    alignas(typename Tail::PreparedTuple)
        std::byte erased[sizeof(typename Tail::PreparedTuple)];
    Tail::prepare_into(COLUMN, erased);
    for (const auto &DIRECTION :
         {math::X_AXIS, math::Y_AXIS, math::Z_AXIS,
          math::Vector(0.3f, 0.4f, 0.5f).normalized()}) {
      const auto FIELD = RetainedHueField::run<typename FX::Binding>(
          {DIRECTION, 0.0f}, COLUMN, {});
      auto pair = Palette::apply_pair(FIELD, COLUMN, PALETTE);
      Color::apply_modulation_pair_in_place(pair.first, pair.second,
                                            MODULATION);
      const auto EXPECTED = finish_color_sample(
          blend_color_samples(pair.first, pair.second, FRACTION));
      const auto ACTUAL = Tail::shade(DIRECTION, TYPED);
      const auto ERASED = Tail::shade_prepared(DIRECTION, COLUMN, erased);
      HS_EXPECT_EQ(ACTUAL.color, EXPECTED.color);
      HS_EXPECT_NEAR(ACTUAL.alpha, EXPECTED.alpha, 0.0f);
      HS_EXPECT_EQ(ERASED.color, ACTUAL.color);
      HS_EXPECT_NEAR(ERASED.alpha, ACTUAL.alpha, 0.0f);
    }
    std::destroy_at(
        std::launder(reinterpret_cast<typename Tail::PreparedTuple *>(erased)));
    for (int value = 0; value < 256; ++value) {
      const float COORDINATE = (value + 0.5f) / 256;
      for (int index = 0; index < 128; ++index) {
        const float Z = 2.0f * (index + 0.5f) / 128 - 1.0f;
        const math::Vector DIRECTION(sqrtf(1.0f - Z * Z), 0.0f, Z);
        const FieldSample FIELD{COORDINATE, 0.75f, DIRECTION, 0.0f};
        auto pair = Palette::apply_pair(FIELD, COLUMN, PALETTE);
        Color::apply_modulation_pair_in_place(pair.first, pair.second,
                                              MODULATION);
        auto carrier = blend_color_samples(pair.first, pair.second, FRACTION);
        HS_EXPECT_TRUE(carrier.hue_rotation == current.hue_rotation_lut);
        HS_EXPECT_TRUE(carrier.hue_rotation_gamut ==
                       current.hue_rotation_gamut);
        const auto ACTUAL = finish_color_sample(carrier);
        carrier.hue_rotation = nullptr;
        carrier.hue_rotation_gamut = nullptr;
        const auto PROCEDURAL = finish_color_sample(carrier);
        different += ACTUAL.color != PROCEDURAL.color;
        const float A = math::wrap_t(math::sinf_0_2pi(
            math::TWO_PI_F * math::wrap_t(2.0f * Z + 0.17f + 0.125f)));
        const float B = math::wrap_t(math::sinf_0_2pi(
            math::TWO_PI_F * math::wrap_t(3.25f * Z + 0.31f + 0.375f)));
        float delta = B - A;
        if (delta >= 0.5f)
          delta -= 1.0f;
        if (delta < -0.5f)
          delta += 1.0f;
        const float HUE_TURNS = math::wrap_t(A + delta * FRACTION);
        auto expected = ColorModulation::apply(current.palette->get(COORDINATE),
                                               HUE_TURNS, before.mechanic);
        expected.alpha *= 0.75f * 0.37f;
        HS_EXPECT_NEAR(ACTUAL.alpha, expected.alpha, 1e-7f);
        for (const uint32_t ERROR :
             {uint32_t(std::abs(int(ACTUAL.color.r) - expected.color.r)),
              uint32_t(std::abs(int(ACTUAL.color.g) - expected.color.g)),
              uint32_t(std::abs(int(ACTUAL.color.b) - expected.color.b))}) {
          maximum = std::max(maximum, ERROR);
          total += ERROR;
          ++channels;
        }
      }
    }
  }
  HS_EXPECT_GT(different, size_t{0});
  HS_EXPECT_LE(float(maximum), Layer::METRICS[0].limit);
  HS_EXPECT_LE(double(total) / channels, Layer::METRICS[1].limit);
  std::printf("  retained hue max=%u mean=%.6f changed=%zu\n", maximum,
              double(total) / channels, different);
}

/** @brief Palette identity rejects altered tables, masks and live palette state. */
inline void test_retained_palette_identity_rejects_changes() {
  using namespace Pullback;
  RetainedHueFixture fixture;
  auto &current = fixture.current;
  auto &previous = fixture.previous;
  auto &history = fixture.history;
  auto altered = current;
  HueRotationGamutMask other_mask = *current.hue_rotation_gamut;
  bool marked = false;
  for (int value = 0; value < HueRotationLutView::VALUE_STEPS - 1 && !marked;
       ++value)
    for (int hue = 0; hue < HueRotationLutView::HUE_STEPS && !marked; ++hue)
      if (!other_mask.requires_direct(value, hue)) {
        other_mask.mark_cell(value, hue);
        marked = true;
      }
  HS_EXPECT_TRUE(marked);
  altered.hue_rotation_gamut = &other_mask;
  const auto MISMATCH = FrameInterval(previous, altered).at(0.5f);
  HS_EXPECT_FALSE(MISMATCH.temporal.same_palette);
  HS_EXPECT_TRUE(MISMATCH.temporal.hue_rotation == previous.hue_rotation_lut);
  HS_EXPECT_TRUE(MISMATCH.temporal.hue_rotation_gamut ==
                 previous.hue_rotation_gamut);
  altered = current;
  altered.hue_rotation_gamut = nullptr;
  HS_EXPECT_FALSE(
      FrameInterval(previous, altered).at(0.5f).temporal.same_palette);
  std::array<Pixel, HueRotationLutView::SIZE> changed_table;
  std::copy_n(current.hue_rotation_lut, changed_table.size(),
              changed_table.begin());
  changed_table[7].r ^= 1;
  altered = current;
  altered.hue_rotation_lut = changed_table.data();
  HS_EXPECT_FALSE(
      FrameInterval(previous, altered).at(0.5f).temporal.same_palette);
  BakedPaletteStorage other_palette;
  other_palette.bake(persistent_arena, Palettes::RICH_SUNSET);
  altered = current;
  altered.palette = &other_palette.view();
  const auto REPAINT = FrameInterval(previous, altered).at(0.5f);
  HS_EXPECT_FALSE(REPAINT.temporal.same_palette);
  HS_EXPECT_TRUE(REPAINT.temporal.hue_rotation == previous.hue_rotation_lut);
  const Pixel SAVED_PIXEL = current.hue_rotation_lut[7];
  const auto SAVED_MASK = *current.hue_rotation_gamut;
  auto *live_table = const_cast<Pixel *>(current.hue_rotation_lut);
  auto *live_mask =
      const_cast<HueRotationGamutMask *>(current.hue_rotation_gamut);
  live_table[7].r ^= 1;
  live_mask->bits[0] ^= 1;
  HS_EXPECT_EQ(history.hue_rotation[7], SAVED_PIXEL);
  HS_EXPECT_TRUE(history.previous.hue_rotation_gamut->bits == SAVED_MASK.bits);
  HS_EXPECT_FALSE(
      FrameInterval(previous, current).at(0.5f).temporal.same_palette);
  live_table[7] = SAVED_PIXEL;
  *live_mask = SAVED_MASK;
  HS_EXPECT_TRUE(history.previous.hue_rotation_gamut->bits ==
                 current.hue_rotation_gamut->bits);
}
