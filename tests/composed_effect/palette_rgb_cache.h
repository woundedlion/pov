/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file palette_rgb_cache.h
 * @brief Interval palette-cache interpolation, admission, fallback and
 * ownership tests. */

/** @brief Palette, periodic hue and time axes interpolate independently. */
void test_palette_rgb_cache_axes_and_wrap() {
  PaletteLabLut base;
  base.active = true;
  PaletteRgbHueLut cache;
  for (int t = 0; t < 3; ++t)
    for (int v = 0; v < 64; ++v)
      for (int h = 0; h < HueRotationLutView::HUE_STEPS; ++h)
        cache.nodes[t * HueRotationLutView::SIZE +
                    v * HueRotationLutView::HUE_STEPS + h] =
            Pixel(v * 100 + t * 1000, h * 200 + t * 100,
                  t * 1000 + v * 20 + h * 10);
  base.rgb_nodes = cache.nodes.data();
  const auto time = base.view(.3f);
  const Pixel result = PaletteRgbHueLut::sample(time, .375f, .98f);
  HS_EXPECT_NEAR(result.r, 2962.5f, 1.0f);
  HS_EXPECT_NEAR(result.g, 1020.0f, 1.0f);
  HS_EXPECT_NEAR(result.b, 1120.5f, 1.0f);
  const Pixel wrapped = PaletteRgbHueLut::sample(time, .375f, -.02f);
  HS_EXPECT_EQ(result, wrapped);
  HS_EXPECT_TRUE(base.view(.3f).rgb_high ==
                 cache.nodes.data() + HueRotationLutView::SIZE);
  HS_EXPECT_TRUE(base.view(.8f).rgb_high ==
                 cache.nodes.data() + 2 * HueRotationLutView::SIZE);
}
/** @brief Column views omit only an identical upper gamut mask. */
void test_palette_rgb_cache_guard_mask_aliasing() {
  PaletteLabLut base;
  PaletteRgbHueLut cache;
  base.active = true;
  base.rgb_nodes = cache.nodes.data();
  base.rgb_masks = cache.masks.data();
  base.rgb_midpoint = true;
  base.rgb_midpoint_row = 1;
  cache.nodes.fill(Pixel(1111, 2222, 3333));
  cache.masks[1].mark_cell(0, 0);
  const auto shared = base.view(.3f);
  HS_EXPECT_TRUE(shared.rgb_mask_low == &cache.masks[1]);
  HS_EXPECT_TRUE(shared.rgb_mask_high == nullptr);
  HS_EXPECT_EQ(PaletteRgbHueLut::sample(shared, .005f, .01f),
               PaletteLabLutBase::direct_hue_color(shared, .005f, .01f));
  HS_EXPECT_NE(PaletteRgbHueLut::sample(shared, .005f, .01f),
               Pixel(1111, 2222, 3333));

  base.rgb_midpoint = false;
  auto distinct = base.view(.3f);
  HS_EXPECT_TRUE(distinct.rgb_mask_low == &cache.masks[0]);
  HS_EXPECT_TRUE(distinct.rgb_mask_high == &cache.masks[1]);
  distinct.rgb_high = distinct.rgb_low;
  HS_EXPECT_EQ(PaletteRgbHueLut::sample(distinct, .005f, .01f),
               PaletteLabLutBase::direct_hue_color(distinct, .005f, .01f));
  HS_EXPECT_NE(PaletteRgbHueLut::sample(distinct, .005f, .01f),
               Pixel(1111, 2222, 3333));
}
/** @brief Replay caches retain their original canonical palette endpoints. */
void test_palette_lab_reference_endpoint_lifetime() {
  reset_globals();
  struct Source {
    Pixel color;
    Color4 get(float) const { return {color, 1.0f}; }
  };
  const Pixel START(1000, 50000, 500);
  const Pixel END(50000, 2000, 8000);
  BakedPaletteStorage from, to, retained;
  from.bake(persistent_arena, Source{START});
  to.bake(persistent_arena, Source{END});
  retained.clone_from(from.view(), persistent_arena);
  PaletteLabLut base;
  base.bake(from.view(), to.view());
  PaletteLabLut reference = base;
  reference.palette_from = &retained.view();
  reference.palette_to = &to.view();
  from.rebake_copy(to.view());
  constexpr float FRACTION = 0.25f;
  constexpr float HUE = 0.137f;
  const Pixel RGB = START.lerp16(END, frac_to_q16(FRACTION));
  const Pixel EXPECTED =
      ColorModulation::HueRotation::apply(Color4(RGB, 1.0f), HUE, {1.0f}).color;
  float cosine, sine;
  hue_sincos(HUE, cosine, sine);
  const Pixel ACTUAL = PaletteLabLutBase::canonical_hue_color(
      reference.view(FRACTION), 0.37f, cosine, sine);
  HS_EXPECT_EQ(ACTUAL, EXPECTED);
}

/** @brief Final hue nodes equal the hue-rotation table of a full-chroma palette. */
void test_palette_rgb_cache_matches_hue_rotation_table() {
  const GenerativePalette source(
      PaletteRecipes::profile(PaletteDomain::STRAIGHT, PaletteHarmony::TRIADIC,
                              AxisCurve::BELL, .37f, 1.0f));
  alignas(std::max_align_t)
      uint8_t storage[BakedPalette::required_arena_bytes()];
  Arena arena(storage, sizeof(storage));
  BakedPaletteStorage baked;
  baked.bake(arena, source);
  const BakedPalette &palette = baked;
  PaletteLabLut base;
  base.bake(palette, palette);
  HS_EXPECT_TRUE(base.active);
  PaletteRgbHueLut cache;
  cache.bake(base);
  HS_EXPECT_TRUE(base.rgb_midpoint);
  std::array<Pixel, HueRotationLutView::SIZE> table;
  prepare_hue_rotation_lut(std::span<Pixel, HueRotationLutView::SIZE>(table),
                           palette);
  for (size_t node = 0; node < HueRotationLutView::SIZE; ++node) {
    const Pixel rgb = cache.nodes[HueRotationLutView::SIZE + node];
    HS_EXPECT_EQ(rgb.r, table[node].r);
    HS_EXPECT_EQ(rgb.g, table[node].g);
    HS_EXPECT_EQ(rgb.b, table[node].b);
  }
  HueRotationGamutMask reference;
  HueRotationGamutBuilder builder(&reference);
  for (int value = 0; value < HueRotationLutView::VALUE_STEPS; ++value)
    for (int hue = 0; hue < HueRotationLutView::HUE_STEPS; ++hue) {
      float cosine, sine;
      hue_sincos(float(hue) / HueRotationLutView::HUE_STEPS, cosine, sine);
      const auto &node =
          cache.nodes[HueRotationLutView::SIZE +
                      value * HueRotationLutView::HUE_STEPS + hue];
      mark_hue_gamut_node(builder, value, hue,
                          base.nodes[PaletteLabLut::VALUES + value], cosine,
                          sine, node);
    }
  HS_EXPECT_EQ(cache.masks[1].bits, reference.bits);
}
/** @brief Lift and brightness resolve the deferred Lab base through the ordinary finalizer. */
void test_palette_rgb_cache_lazy_guards() {
  PaletteLabLut base;
  base.active = true;
  for (auto &lab : base.nodes)
    lab = {.7f, .1f, .05f};
  PaletteRgbHueLut cache;
  for (auto &pixel : cache.nodes)
    pixel = Pixel(1111, 2222, 3333);
  base.rgb_nodes = cache.nodes.data();
  const auto time = base.view(.3f);
  Pullback::ColorSample sample{Color4(Pixel(0, 0, 0), .37f),
                               {.4f, .8f, {0, 1, 0}, .2f}};
  sample.cached_rgb = &time;
  sample.has_cached_lab = true;
  sample.palette_coordinate = .4f;
  sample.pending_hue = .12f;
  HS_EXPECT_EQ(Pullback::finish_color_sample(sample).color,
               Pixel(1111, 2222, 3333));
  HS_EXPECT_EQ(Pullback::finish_color_sample(sample).alpha, .37f);
  for (int change = 0; change < 3; ++change) {
    auto candidate = sample;
    if (change == 0)
      candidate.pending_lift = .2f;
    if (change == 1)
      candidate.brightness_gain = .5f;
    if (change == 2)
      candidate.brightness_gain = 0.0f;
    auto ordinary = candidate;
    ordinary.cached_rgb = nullptr;
    ordinary.cached_lab =
        PaletteLabLut::sample(time, ordinary.palette_coordinate);
    HS_EXPECT_EQ(Pullback::finish_color_sample(candidate).color,
                 Pullback::finish_color_sample(ordinary).color);
    HS_EXPECT_EQ(Pullback::finish_color_sample(candidate).alpha, .37f);
  }
}
/** @brief Runtime ownership retains the preceding interval and rejects lift recipes. */
void test_palette_rgb_cache_runtime_lifetime() {
  reset_globals();
  pin_frame_clock(0);
  AshCloud<96, 20> effect;
  effect.init();
  effect.draw_frame();
  effect.advance_display();
  pin_frame_clock(1);
  effect.draw_frame();
  effect.advance_display();
  const auto first = ComposedFrameWhiteBox::current_frame(effect);
  HS_EXPECT_TRUE(first.palette_lab_lut && first.palette_lab_lut->rgb_nodes);
  if (!first.palette_lab_lut || !first.palette_lab_lut->rgb_nodes)
    return;
  const Pixel *previous = first.palette_lab_lut->rgb_nodes;
  constexpr size_t RGB_SIZE =
      PaletteRgbHueLutImpl<AshCloudSpec::PALETTE_RGB_MIDPOINT_ONLY>::RGB_TIMES *
      HueRotationLutView::SIZE;
  std::array<Pixel, RGB_SIZE> saved;
  std::copy_n(previous, RGB_SIZE, saved.begin());
  pin_frame_clock(2);
  effect.draw_frame();
  effect.advance_display();
  const auto second = ComposedFrameWhiteBox::current_frame(effect);
  HS_EXPECT_TRUE(second.palette_lab_lut != nullptr);
  if (!second.palette_lab_lut)
    return;
  HS_EXPECT_TRUE(second.palette_lab_lut->rgb_nodes != previous);
  HS_EXPECT_TRUE(std::equal(saved.begin(), saved.end(), previous));
  auto params = effect.initial_params();
  params.template get<"modulation">().modulation.mechanic.kind =
      ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT;
  params.template get<"modulation">().modulation.mechanic.amount = .2f;
  ComposedFrameWhiteBox::set_params(effect, params);
  pin_frame_clock(3);
  effect.draw_frame();
  effect.advance_display();
  const auto *lift_base =
      ComposedFrameWhiteBox::current_frame(effect).palette_lab_lut;
  HS_EXPECT_TRUE(lift_base != nullptr);
  if (!lift_base)
    return;
  HS_EXPECT_TRUE(lift_base->rgb_nodes == nullptr);
}
/** @brief Exact retained rows match fresh bakes; changed identities take the full bake. */
void test_palette_rgb_cache_bake_reuse() {
  PaletteLabLut previous;
  previous.active = true;
  for (int time = 0; time < 3; ++time)
    for (int value = 0; value < 64; ++value)
      previous.nodes[time * 64 + value] = {.4f + .1f * time + .001f * value,
                                           .06f, .04f};
  PaletteRgbHueLut old_rgb;
  old_rgb.bake(previous);
  PaletteLabLut current;
  current.active = true;
  for (int time = 0; time < 3; ++time)
    for (int value = 0; value < 64; ++value)
      current.nodes[time * 64 + value] = {
          previous.nodes[128 + value].L + .02f * time, .06f, .04f};
  auto reference_base = current;
  PaletteRgbHueLut reference, reused;
  reference.bake(reference_base);
  reused.bake(current, &previous, &old_rgb);
  HS_EXPECT_EQ(std::memcmp(reference.nodes.data(), reused.nodes.data(),
                           sizeof(reference.nodes)),
               0);
  for (int time : {0, 1, 2})
    for (int value : {0, 31, 63})
      for (int hue : {0, 1, 7, 15}) {
        const auto expected = linrgb_to_pixel(hue_rotate_lut_gamut_rgb(
            current.nodes[time * 64 + value],
            float(hue) / HueRotationLutView::HUE_STEPS));
        const auto actual =
            reused.nodes[time * HueRotationLutView::SIZE +
                         value * HueRotationLutView::HUE_STEPS + hue];
        HS_EXPECT_EQ(actual.r, expected.r);
        HS_EXPECT_EQ(actual.g, expected.g);
        HS_EXPECT_EQ(actual.b, expected.b);
      }
  current.nodes[0].L += .03f;
  reference_base = current;
  reference.bake(reference_base);
  reused.bake(current, &previous, &old_rgb);
  HS_EXPECT_EQ(std::memcmp(reference.nodes.data(), reused.nodes.data(),
                           sizeof(reference.nodes)),
               0);
  current.nodes[0] = previous.nodes[128];
  reference_base = current;
  reference.bake(reference_base);
  PaletteRgbHueLut wrong_binding;
  for (auto &pixel : wrong_binding.nodes)
    pixel = Pixel(1000, 0, 0);
  previous.rgb_nodes = wrong_binding.nodes.data();
  for (auto &pixel : old_rgb.nodes)
    pixel = Pixel(0, 1000, 0);
  reused.bake(current, &previous, &old_rgb);
  HS_EXPECT_EQ(std::memcmp(reference.nodes.data(), reused.nodes.data(),
                           sizeof(reference.nodes)),
               0);
  current.active = false;
  reused.bake(current);
  HS_EXPECT_TRUE(current.rgb_nodes == nullptr);
}
/** @brief Literal endpoints and missing RGB views retain the ordinary palette paths. */
void test_palette_rgb_cache_shared_base_endpoints() {
  using FX = AshCloud<96, 20>;
  using Policy = Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
      typename FX::Binding, Pullback::Color::BrightnessEnvelope::NONE, true>>;
  PaletteLabLut base;
  base.active = true;
  for (auto &lab : base.nodes)
    lab = {.6f, .05f, .01f};
  PaletteRgbHueLut cache;
  cache.bake(base);
  typename Policy::Prepared prepared{};
  prepared.current.mapping = Pullback::Color::PaletteMappingWeights::single(
      Pullback::Color::PaletteMapping::LINEAR);
  prepared.current.mapping_frequency = 1.0f;
  prepared.current.mapping_offset = 0.0f;
  prepared.lab = base.view(.5f);
  FX::FrameState frame{};
  Pullback::FieldSample field{.4f, .7f, {0, 1, 0}, .2f};
  Pullback::ColorSample output;
  prepared.fraction = 0.0f;
  HS_EXPECT_FALSE(Policy::apply_shared_base(field, frame, prepared, output));
  prepared.fraction = 1.0f;
  HS_EXPECT_TRUE(Policy::apply_shared_base(field, frame, prepared, output));
  prepared.fraction = .5f;
  HS_EXPECT_TRUE(Policy::apply_shared_base(field, frame, prepared, output));
  HS_EXPECT_TRUE(output.cached_rgb == &prepared.lab);
  HS_EXPECT_EQ(output.cached_lab.L, 0.0f);
  prepared.lab.rgb_low = nullptr;
  prepared.lab.rgb_high = nullptr;
  HS_EXPECT_TRUE(Policy::apply_shared_base(field, frame, prepared, output));
  HS_EXPECT_TRUE(output.cached_rgb == nullptr);
  HS_EXPECT_NEAR(output.cached_lab.L, .6f, 1e-6f);
}

/** @brief Midpoint reuse preserves its original baked key and rejects accumulated drift. */
void test_palette_rgb_cache_midpoint_reuse() {
  const auto make_base = [](float shift, float step) {
    PaletteLabLut base;
    base.active = true;
    for (int time = 0; time < 3; ++time)
      for (int value = 0; value < 64; ++value)
        base.nodes[time * 64 + value] = {.6f + shift + (time - 1) * step,
                                         .08f + .0001f * value, .03f};
    return base;
  };
  auto original = make_base(0.0f, .001f);
  PaletteRgbHueLut first;
  first.bake(original);
  HS_EXPECT_TRUE(original.rgb_midpoint);
  HS_EXPECT_TRUE(first.midpoint_active);
  HS_EXPECT_TRUE(original.view(.1f).rgb_low == original.view(.9f).rgb_high);
  HS_EXPECT_TRUE(original.view(.1f).rgb_low ==
                 first.nodes.data() + HueRotationLutView::SIZE);
  const auto old_pixels = first.nodes;
  auto next = make_base(.0008f, .001f);
  PaletteRgbHueLut second;
  second.bake(next, &original, &first);
  HS_EXPECT_TRUE(std::equal(old_pixels.begin() + HueRotationLutView::SIZE,
                            old_pixels.begin() + 2 * HueRotationLutView::SIZE,
                            second.nodes.begin() + HueRotationLutView::SIZE));
  HS_EXPECT_EQ(second.baked_midpoint[0].L, .6f);
  auto shifted = make_base(.0016f, .001f);
  PaletteRgbHueLut third;
  third.bake(shifted, &next, &second);
  HS_EXPECT_EQ(third.baked_midpoint[0].L, shifted.nodes[64].L);
  HS_EXPECT_FALSE(std::equal(old_pixels.begin() + HueRotationLutView::SIZE,
                             old_pixels.begin() + 2 * HueRotationLutView::SIZE,
                             third.nodes.begin() + HueRotationLutView::SIZE));
  auto discontinuous = make_base(0.0f, .003f);
  third.bake(discontinuous, &shifted, &third);
  HS_EXPECT_FALSE(discontinuous.rgb_midpoint);
  HS_EXPECT_FALSE(third.midpoint_active);
  for (int time : {0, 1, 2}) {
    const auto lab = discontinuous.nodes[time * 64 + 31];
    const auto expected = linrgb_to_pixel(
        hue_rotate_lut_gamut_rgb(lab, 5.f / HueRotationLutView::HUE_STEPS));
    const auto actual = third.nodes[time * HueRotationLutView::SIZE +
                                    31 * HueRotationLutView::HUE_STEPS + 5];
    HS_EXPECT_EQ(actual.r, expected.r);
    HS_EXPECT_EQ(actual.g, expected.g);
    HS_EXPECT_EQ(actual.b, expected.b);
  }
}

/** @brief Every palette chroma bakes the final RGB cache. */
void test_palette_rgb_cache_any_chroma() {
  reset_globals();
  pin_frame_clock(0);
  AshCloud<96, 20> effect;
  effect.init();
  effect.draw_frame();
  effect.advance_display();
  auto params = effect.initial_params();
  params.template get<"color">().palette_chroma = .5f;
  ComposedFrameWhiteBox::set_params(effect, params);
  for (int frame = 1; frame < 3; ++frame) {
    pin_frame_clock(frame);
    effect.draw_frame();
    effect.advance_display();
  }
  const auto &current = ComposedFrameWhiteBox::current_frame(effect);
  HS_EXPECT_TRUE(current.palette_lab_lut != nullptr);
  if (!current.palette_lab_lut)
    return;
  HS_EXPECT_TRUE(current.palette_lab_lut->rgb_nodes != nullptr);
}

/** @brief Linear table interpolation stays within one RGB16 code of Q16 interpolation. */
void test_palette_rgb_cache_q16_compatibility() {
  uint32_t random = 0x728f1ab5;
  const auto next = [&]() {
    random ^= random << 13;
    random ^= random >> 17;
    random ^= random << 5;
    return random;
  };
  unsigned maximum = 0;
  PaletteLabLut base;
  base.active = true;
  PaletteRgbHueLut cache;
  std::array<Pixel, 3 * HueRotationLutView::SIZE> reference;
  base.rgb_nodes = cache.nodes.data();
  for (int table = 0; table < 8; ++table) {
    for (size_t i = 0; i < reference.size(); ++i) {
      reference[i] =
          Pixel(uint16_t(next()), uint16_t(next()), uint16_t(next()));
      cache.nodes[i] = reference[i];
    }
    for (int sample = 0; sample < 2048; ++sample) {
      const float VALUE = float(next() & 0xffff) / 65535.0f * 1.5f - .25f;
      const float HUE = float(next() & 0xffff) / 65535.0f * 8 - 4;
      const float TIME = float(next() & 0xffff) / 65535.0f;
      const auto view = base.view(TIME);
      const size_t LOW_ROW =
          (view.rgb_low - cache.nodes.data()) / HueRotationLutView::SIZE;
      const auto a = sample_hue_rotation_lut(
          {reference.data() + LOW_ROW * HueRotationLutView::SIZE, true}, VALUE,
          HUE);
      const auto b = sample_hue_rotation_lut(
          {reference.data() + (LOW_ROW + 1) * HueRotationLutView::SIZE, true},
          VALUE, HUE);
      const auto expected = a.lerp16(b, frac_to_q16(view.weight));
      const auto actual = PaletteRgbHueLut::sample(view, VALUE, HUE);
      maximum =
          std::max({maximum, unsigned(abs(int(actual.r) - int(expected.r))),
                    unsigned(abs(int(actual.g) - int(expected.g))),
                    unsigned(abs(int(actual.b) - int(expected.b)))});
    }
  }
  HS_EXPECT_LE(maximum, 1u);
}

/** @brief Padding bytes cannot change semantic palette-cache eligibility. */
inline void test_palette_lab_cache_color_padding() {
  using FX = AshCloud<96, 20>;
  reset_globals();
  pin_frame_clock(0);
  FX effect;
  effect.init();
  effect.draw_frame();
  effect.advance_display();
  auto &previous = const_cast<typename FX::FrameState &>(
                       ComposedFrameWhiteBox::previous_frame(effect))
                       .params.template get<"color">();
  const auto &current = ComposedFrameWhiteBox::current_frame(effect)
                            .params.template get<"color">();
  HS_EXPECT_TRUE(previous == current);
  auto *bytes = reinterpret_cast<unsigned char *>(&previous);
  const auto *current_bytes = reinterpret_cast<const unsigned char *>(&current);
  const size_t END =
      reinterpret_cast<const unsigned char *>(&previous.palette_mapping) -
      bytes + sizeof(previous.palette_mapping);
  HS_EXPECT_LE(END, sizeof(previous));
  if (END > sizeof(previous))
    return;
  for (size_t byte = END; byte < sizeof(previous); ++byte)
    bytes[byte] = current_bytes[byte] ^ 0xff;
  HS_EXPECT_TRUE(previous == current);
  if (END < sizeof(previous))
    HS_EXPECT_TRUE(std::memcmp(&previous, &current, sizeof(previous)) != 0);
  pin_frame_clock(1);
  effect.draw_frame();
  effect.advance_display();
  const auto *base =
      ComposedFrameWhiteBox::current_frame(effect).palette_lab_lut;
  HS_EXPECT_TRUE(base != nullptr);
  if (base)
    HS_EXPECT_TRUE(base->active);
  auto params = effect.initial_params();
  params.template get<"color">().palette_chroma += .01f;
  ComposedFrameWhiteBox::set_params(effect, params);
  pin_frame_clock(2);
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_TRUE(ComposedFrameWhiteBox::current_frame(effect).palette_lab_lut ==
                 nullptr);
}

/** @brief Checks palette-control classification of every tabled color field.
 * @tparam Controls Color family type. */
template <typename Controls> void check_palette_control_fields() {
  constexpr std::string_view PALETTE_FIELDS[] = {
      "palette-chroma",    "palette-hue",      "palette-spread",
      "palette-lightness", "palette-mutation", "palette-lightness-range"};
  for (size_t index = 0; index < Controls::FIELDS.size(); ++index) {
    const auto &field = Controls::FIELDS[index];
    HS_CONTEXT(field.id);
    const std::string_view ID = field.id;
    const bool PALETTE =
        ID.starts_with("gradient-") ||
        std::find(std::begin(PALETTE_FIELDS), std::end(PALETTE_FIELDS), ID) !=
            std::end(PALETTE_FIELDS);
    const Controls BASE{};
    Controls changed = BASE;
    float &value = changed.*field.member;
    const float STEP = .25f * (field.max - field.min);
    value = value + STEP <= field.max ? value + STEP : value - STEP;
    HS_EXPECT_EQ(Pullback::Color::same_palette_controls(BASE, changed),
                 !PALETTE);
    HS_EXPECT_EQ(Pullback::Color::same_palette_controls(changed, BASE),
                 !PALETTE);
  }
}

/** @brief Palette-cache admission ignores exactly the column-sampled color
 * controls. */
inline void test_palette_lab_cache_column_controls() {
  using Pullback::Color::same_palette_controls;
  check_palette_control_fields<Pullback::ColorParams>();
  check_palette_control_fields<Pullback::Color::AuthoredColorParams>();
  const Pullback::Color::AuthoredColorParams BASE{};
  for (int change = 0; change < 6; ++change) {
    HS_CONTEXT("discrete control", change);
    auto changed = BASE;
    if (change == 0)
      changed.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    else if (change == 1)
      changed.brightness_envelope = Pullback::Color::BrightnessEnvelope::BELL;
    else if (change == 2)
      changed.palette_mode = 1;
    else if (change == 3)
      changed.palette_type = 1;
    else if (change == 4)
      changed.palette_harmony = 0;
    else
      changed.palette_name = 3;
    HS_EXPECT_EQ(same_palette_controls(BASE, changed), change < 2);
  }
}

/** @brief A Lerp of column-sampled color controls keeps the palette cache
 * within its qualified bounds of canonical shading. */
inline void test_palette_lab_cache_through_lerp() {
  using FX = KaleidoscopeHexOil<SMALL_W, SMALL_H>;
  using WB = ComposedWhiteBox<FX>;
  using Pipeline = FX::RenderPipeline;
  using Policy = Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
      typename FX::Binding, Pullback::Color::BrightnessEnvelope::NONE>>;
  static_assert(Policy::LAB_CACHE);
  reset_globals();
  ScopedTemporalSampling temporal(true);
  pin_frame_clock(0);
  FX effect;
  effect.init();
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_TRUE(WB::advance_preset(effect));
  HS_EXPECT_TRUE(WB::transition_active(effect));
  Spacetime::FrameHistory<FX::FrameState, true, true> history;
  constexpr int FRAMES = 8;
  for (int frame = 1; frame <= FRAMES; ++frame) {
    HS_CONTEXT("lerp frame", frame);
    history.capture(ComposedFrameWhiteBox::current_frame(effect),
                    persistent_arena);
    pin_frame_clock(frame);
    effect.draw_frame();
    const auto &current = ComposedFrameWhiteBox::current_frame(effect);
    HS_EXPECT_TRUE(current.palette_lab_lut && current.palette_lab_lut->active);
    if (frame < FRAMES)
      effect.advance_display();
  }
  HS_EXPECT_TRUE(WB::transition_active(effect));
  const auto CURRENT = ComposedFrameWhiteBox::current_frame(effect);
  HS_EXPECT_TRUE(
      history.previous.params.template get<"color">().mapping_frequency !=
      CURRENT.params.template get<"color">().mapping_frequency);
  uint32_t maximum = 0;
  uint64_t total = 0;
  size_t channels = 0;
  float alpha_error = 0.0f;
  for (int x = 0; x < SMALL_W; ++x) {
    const auto CACHED_STATE =
        Spacetime::sample_frame(history.previous, CURRENT,
                                Spacetime::SamplingSurface<SMALL_W>::time(x));
    auto canonical_state = CACHED_STATE;
    canonical_state.palette_lab_lut = nullptr;
    const auto CACHED = Pipeline::prepare(CACHED_STATE);
    const auto CANONICAL = Pipeline::prepare(canonical_state);
    for (int y = 0; y < SMALL_H; ++y) {
      const float SIN_PHI = math::TrigLUT<SMALL_W, SMALL_H>::sin_phi[y];
      const math::Vector VIEW(
          SIN_PHI * math::TrigLUT<SMALL_W, SMALL_H>::cos_theta(x),
          math::TrigLUT<SMALL_W, SMALL_H>::cos_phi[y],
          SIN_PHI * math::TrigLUT<SMALL_W, SMALL_H>::sin_theta[x]);
      const Color4 ACTUAL = Pipeline::shade(VIEW, CACHED);
      const Color4 EXPECTED = Pipeline::shade(VIEW, CANONICAL);
      alpha_error = std::max(alpha_error, fabsf(ACTUAL.alpha - EXPECTED.alpha));
      for (const uint32_t ERROR :
           {uint32_t(std::abs(int(ACTUAL.color.r) - int(EXPECTED.color.r))),
            uint32_t(std::abs(int(ACTUAL.color.g) - int(EXPECTED.color.g))),
            uint32_t(std::abs(int(ACTUAL.color.b) - int(EXPECTED.color.b)))}) {
        maximum = std::max(maximum, ERROR);
        total += ERROR;
        ++channels;
      }
    }
  }
  const double MEAN = static_cast<double>(total) / channels;
  std::printf("  lerp palette cache max=%u mean=%.3f channel codes\n", maximum,
              MEAN);
  HS_EXPECT_GT(maximum, 0u);
  HS_EXPECT_LE(static_cast<float>(maximum), Policy::METRICS[0].limit);
  HS_EXPECT_LE(MEAN, Policy::METRICS[1].limit);
  HS_EXPECT_LE(alpha_error, 1e-6f);
}

/** @brief A palette control step disables the palette cache for one interval.
 */
inline void test_palette_lab_cache_palette_control_step() {
  using FX = KaleidoscopeHexOil<SMALL_W, SMALL_H>;
  reset_globals();
  ScopedTemporalSampling temporal(true);
  pin_frame_clock(0);
  FX effect;
  effect.init();
  effect.draw_frame();
  effect.advance_display();
  auto params = effect.initial_params();
  const auto draw = [&](int frame) {
    ComposedFrameWhiteBox::set_params(effect, params);
    pin_frame_clock(frame);
    effect.draw_frame();
    effect.advance_display();
    const auto *base =
        ComposedFrameWhiteBox::current_frame(effect).palette_lab_lut;
    return base && base->active;
  };
  auto &color = params.template get<"color">();
  color.mapping_frequency *= 1.5f;
  color.opacity_low *= .5f;
  HS_EXPECT_TRUE(draw(1));
  color.palette_chroma += .01f;
  HS_EXPECT_FALSE(draw(2));
  HS_EXPECT_TRUE(draw(3));
}

/** @brief Compact midpoint rows and retained keys match the generic owner exactly. */
void test_palette_rgb_midpoint_only_storage() {
  static_assert(!Pullback::Spec::PALETTE_RGB_MIDPOINT_ONLY);
  static_assert(AshCloudSpec::PALETTE_RGB_MIDPOINT_ONLY);
  static_assert(sizeof(Pullback::OptionalPaletteRgbLut<true>) -
                    sizeof(Pullback::OptionalPaletteRgbLut<true, true>) ==
                25080);
  std::array<PaletteRgbHueLut, 2> generic;
  std::array<PaletteRgbHueLutImpl<true>, 2> compact;
  std::array<PaletteLabLut, 2> full_bases, small_bases;
  for (int frame = 0; frame < 80; ++frame) {
    const int SLOT = frame & 1, PREVIOUS = SLOT ^ 1;
    auto &a = full_bases[SLOT];
    a.active = true;
    for (int time = 0; time < 3; ++time)
      for (int value = 0; value < 64; ++value)
        a.nodes[time * 64 + value] = {.5f + frame * .0002f + (time - 1) * .001f,
                                      .04f + value * .0001f, .03f};
    small_bases[SLOT] = a;
    generic[SLOT].bake(a, frame ? &full_bases[PREVIOUS] : nullptr,
                       frame ? &generic[PREVIOUS] : nullptr);
    compact[SLOT].bake(small_bases[SLOT],
                       frame ? &small_bases[PREVIOUS] : nullptr,
                       frame ? &compact[PREVIOUS] : nullptr);
    const auto &b = small_bases[SLOT];
    HS_EXPECT_EQ(b.rgb_midpoint_row, 0);
    HS_EXPECT_EQ(
        std::memcmp(generic[SLOT].nodes.data() + HueRotationLutView::SIZE,
                    compact[SLOT].nodes.data(),
                    HueRotationLutView::SIZE * sizeof(Pixel)),
        0);
    HS_EXPECT_EQ(std::memcmp(generic[SLOT].baked_midpoint.data(),
                             compact[SLOT].baked_midpoint.data(),
                             sizeof(generic[SLOT].baked_midpoint)),
                 0);
    const auto x = a.view(.3f), y = b.view(.3f);
    HS_EXPECT_TRUE(y.rgb_low == compact[SLOT].nodes.data() &&
                   y.rgb_high == y.rgb_low);
    HS_EXPECT_EQ(PaletteRgbHueLut::sample(x, .37f, .89f),
                 PaletteRgbHueLut::sample(y, .37f, .89f));
  }
}

/** @brief Unqualified compact intervals publish no RGB view and retain their old samples. */
void test_palette_rgb_midpoint_only_fallback() {
  PaletteLabLut base;
  base.active = true;
  for (auto &lab : base.nodes)
    lab = {.6f, .05f, .02f};
  PaletteRgbHueLutImpl<true> cache;
  cache.bake(base);
  const auto saved = cache.nodes;
  base.nodes[128].L += .01f;
  auto full_base = base;
  PaletteRgbHueLut full;
  full.bake(full_base);
  cache.bake(base);
  HS_EXPECT_TRUE(full_base.rgb_nodes != nullptr);
  HS_EXPECT_TRUE(base.rgb_nodes == nullptr);
  HS_EXPECT_FALSE(cache.midpoint_active);
  HS_EXPECT_TRUE(base.view(.5f).rgb_low == nullptr &&
                 base.view(.5f).rgb_high == nullptr);
  HS_EXPECT_TRUE(base.view(.5f).low != nullptr);
  HS_EXPECT_EQ(std::memcmp(saved.data(), cache.nodes.data(), sizeof(saved)), 0);
}
