/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Pullback palette, Lab cache and modulation finishing tests.
 */

/** @brief Owns temporal palette endpoints and optional extra arena storage. */
template <size_t Extra = 0> struct TemporalPalettes {
  std::array<uint8_t, 2 * BakedPalette::required_arena_bytes() + Extra>
      storage{};
  Arena arena{storage.data(), storage.size()};
  BakedPaletteStorage a, b;

  template <typename A, typename B> void bake(const A &from, const B &to) {
    a.bake(arena, from);
    b.bake(arena, to);
  }
  template <typename Frame> void bind(Frame &frame) const {
    frame.palette = &b.view();
    frame.temporal.palette = &a.view();
  }
};

/** @brief Prepares sine modulation endpoints at three-quarter and quarter turns. */
inline void
prepare_sine_clock_pair(Pullback::Modulation::Resources<false> &previous,
                        Pullback::Modulation::Resources<false> &current) {
  Pullback::Modulation::Config config;
  config.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  config.evaluation = Pattern::Evaluation::DIRECT;
  config.sine.frequency = 0.0f;
  auto clock = Pullback::Modulation::init(config);
  clock.phase.phase_turns = 0.75f;
  previous.prepare(config, clock);
  clock.phase.phase_turns = 0.25f;
  current.prepare(config, clock);
}

/** @brief Palette and modulation endpoints for color-stage tests. */
struct PaletteFrame {
  const BakedPalette *palette = nullptr;
  const Pixel *hue_rotation = nullptr;
  const Pullback::Modulation::Prepared *modulation = nullptr;
  const Pullback::Modulation::Prepared *previous_modulation = nullptr;
  ColorModulation::Params mechanic{ColorModulation::Mechanic::HUE_ROTATION,
                                   0.25f};
  Pullback::Color::TemporalPaletteResources temporal{};
  const Pullback::Color::PhaseModulationState *phase = nullptr;
  const PaletteLabLut *lab = nullptr;
};
struct PaletteBinding {
  using FrameState = PaletteFrame;
  using Instrumentation = Pullback::NoInstrumentation;
};
struct PaletteProvider {
  using Binding = PaletteBinding;
  using FrameState = PaletteFrame;
  static auto mapping_weights(const PaletteFrame &) {
    return Pullback::Color::PaletteMappingWeights::single(
        Pullback::Color::PaletteMapping::LINEAR);
  }
  static float mapping_frequency(const PaletteFrame &) { return 1.0f; }
  static float mapping_phase(const PaletteFrame &) { return 0.0f; }
  static float oscillation_depth(const PaletteFrame &) { return 0.0f; }
  static float oscillation_phase(const PaletteFrame &) { return 0.0f; }
  static const BakedPalette &palette(const PaletteFrame &frame) {
    return *frame.palette;
  }
  static HueRotationLutView hue_rotation(const PaletteFrame &frame) {
    return {frame.hue_rotation, frame.hue_rotation != nullptr};
  }
  static auto brightness_envelope(const PaletteFrame &) {
    return Pullback::Color::BrightnessEnvelope::NONE;
  }
  static float brightness_bottom(const PaletteFrame &) { return 1.0f; }
  static float brightness_top(const PaletteFrame &) { return 1.0f; }
  static float opacity_low(const PaletteFrame &) { return 0.2f; }
  static float opacity_high(const PaletteFrame &) { return 0.6f; }
};
struct TemporalPaletteProvider : PaletteProvider {
  static constexpr bool FUSED_TEMPORAL_COLOR = true;
  static auto temporal_resources(const PaletteFrame &frame) {
    return frame.temporal;
  }
};
struct PaletteLayerProvider {
  using Binding = PaletteBinding;
  using FrameState = PaletteFrame;
  static Pullback::Color::ModulationState prepare(const PaletteFrame &frame) {
    return {frame.modulation, frame.previous_modulation,
            Pattern::SourceKind::DIRECTIONAL_SINE, frame.mechanic, true};
  }
};
struct PaletteHueLayerProvider : PaletteLayerProvider {
  static Pullback::Color::ModulationState prepare(const PaletteFrame &frame) {
    auto result = PaletteLayerProvider::prepare(frame);
    result.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 0.25f};
    return result;
  }
};
struct PaletteSampleStage
    : Pullback::Stage::Contract<PaletteSampleStage, Pullback::SphereSample,
                                Pullback::FieldSample> {
  using Policies = std::tuple<>;
  template <typename B>
  static Pullback::FieldSample run(const Pullback::SphereSample &input,
                                   const PaletteFrame &,
                                   const Pullback::NoPrepared &) {
    return {0.5f, 0.5f, input.dir, input.path_length};
  }
};
struct ConstantPalette {
  Color4 color;
  Color4 get(float) const { return color; }
};

struct VariableAlphaPalette {
  Color4 get(float value) const {
    const float INDEX = value * (BakedPalette::LUT_SIZE - 1);
    return Color4(Pixel(1000, 2000, 3000),
                  INDEX > 0.5f && INDEX < 1.5f ? 0.5f : 1.0f);
  }
};

struct LabRampPalette {
  Color4 get(float value) const {
    return Color4(Pixel(static_cast<uint16_t>(value * 50000),
                        static_cast<uint16_t>(value * value * 30000),
                        static_cast<uint16_t>((1.0f - value) * 20000)),
                  1.0f);
  }
};

/** @brief Palette Lab nodes convert RGB-first blends and reject varying alpha. */
inline void test_pullback_palette_lab_nodes() {
  std::array<uint8_t, 5 * BakedPalette::required_arena_bytes()> buffer{};
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage a, b, variable;
  a.bake(arena, ConstantPalette{Color4(Pixel(0, 0, 0), 1.0f)});
  b.bake(arena, ConstantPalette{Color4(Pixel(65535, 65535, 65535), 1.0f)});
  variable.bake(arena, VariableAlphaPalette{});
  HS_EXPECT_LT(variable.view().get_alpha(1.0f / 255.0f), 1.0f);
  for (int index = 0; index < PaletteLabLut::VALUES; ++index)
    HS_EXPECT_EQ(
        variable.view().get_alpha(float(index) / (PaletteLabLut::VALUES - 1)),
        1.0f);
  PaletteLabLut cache;
  cache.bake(a.view(), b.view());
  HS_EXPECT_TRUE(cache.active);
  const auto MID = PaletteLabLut::sample(cache.view(0.5f), 0.37f);
  HS_EXPECT_NEAR(MID.L, std::cbrt(32768.0 / 65535.0), 3e-5);
  HS_EXPECT_GT(MID.L, 0.75f);
  HS_EXPECT_NEAR(MID.a, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(MID.b, 0.0f, 1e-6f);
  for (int row = 0; row < PaletteLabLut::TIMES; ++row)
    for (int index = 0; index < PaletteLabLut::VALUES; ++index) {
      const float FRACTION = float(row) / (PaletteLabLut::TIMES - 1);
      const float VALUE = float(index) / (PaletteLabLut::VALUES - 1);
      const Pixel RGB = a.view().get_color(VALUE).lerp16(
          b.view().get_color(VALUE), frac_to_q16(FRACTION));
      const auto LINEAR = pixel_to_linrgb(RGB);
      const auto EXPECTED =
          linear_rgb_to_oklab_fast(LINEAR.r, LINEAR.g, LINEAR.b);
      const auto ACTUAL = PaletteLabLut::sample(cache.view(FRACTION), VALUE);
      HS_EXPECT_NEAR(ACTUAL.L, EXPECTED.L, 1e-7f);
      HS_EXPECT_NEAR(ACTUAL.a, EXPECTED.a, 1e-7f);
      HS_EXPECT_NEAR(ACTUAL.b, EXPECTED.b, 1e-7f);
    }
  const auto OLD_VIEW = cache.view(0.25f);
  PaletteLabLut next;
  next.bake(b.view(), a.view());
  HS_EXPECT_NEAR(PaletteLabLut::sample(OLD_VIEW, 0.37f).L, MID.L * 0.5f, 1e-7f);
  BakedPaletteStorage ramp, translucent;
  ramp.bake(arena, LabRampPalette{});
  translucent.bake(arena, ConstantPalette{Color4(Pixel(123, 456, 789), 0.6f)});
  next.bake(ramp.view(), b.view());
  for (float value : {0.183f, 0.677f}) {
    const float POSITION =
        value * static_cast<float>(PaletteLabLut::VALUES - 1);
    const int INDEX = static_cast<int>(POSITION);
    const float WEIGHT = POSITION - INDEX;
    const auto node = [&](int index, float fraction) {
      const float COORDINATE =
          float(index) / static_cast<float>(PaletteLabLut::VALUES - 1);
      const Pixel RGB =
          ramp.view()
              .get_color(COORDINATE)
              .lerp16(b.view().get_color(COORDINATE), frac_to_q16(fraction));
      const LinRGB LINEAR = pixel_to_linrgb(RGB);
      return linear_rgb_to_oklab_fast(LINEAR.r, LINEAR.g, LINEAR.b);
    };
    const auto A = node(INDEX, 0.0f), B = node(INDEX + 1, 0.0f);
    const auto C = node(INDEX, 0.5f), D = node(INDEX + 1, 0.5f);
    const auto ACTUAL = PaletteLabLut::sample(next.view(0.17f), value);
    HS_EXPECT_NEAR(ACTUAL.L,
                   hs::lerp(hs::lerp(A.L, B.L, WEIGHT),
                            hs::lerp(C.L, D.L, WEIGHT),
                            (0.17f * (PaletteLabLut::TIMES - 1))),
                   1e-7f);
    HS_EXPECT_NEAR(ACTUAL.a,
                   hs::lerp(hs::lerp(A.a, B.a, WEIGHT),
                            hs::lerp(C.a, D.a, WEIGHT),
                            (0.17f * (PaletteLabLut::TIMES - 1))),
                   1e-7f);
    HS_EXPECT_NEAR(ACTUAL.b,
                   hs::lerp(hs::lerp(A.b, B.b, WEIGHT),
                            hs::lerp(C.b, D.b, WEIGHT),
                            (0.17f * (PaletteLabLut::TIMES - 1))),
                   1e-7f);
  }
  next.bake(ramp.view(), translucent.view());
  HS_EXPECT_FALSE(next.active);
  cache.bake(a.view(), variable.view());
  HS_EXPECT_FALSE(cache.active);
  HS_EXPECT_TRUE(cache.view(0.5f).low == nullptr);
}

template <typename From, typename To>
concept CanBakePaletteLab = requires(PaletteLabLut &cache) {
  cache.bake(std::declval<From>(), std::declval<To>());
};
static_assert(CanBakePaletteLab<const BakedPalette &, const BakedPalette &>);
static_assert(!CanBakePaletteLab<BakedPalette, const BakedPalette &>);
static_assert(!CanBakePaletteLab<const BakedPalette &, BakedPalette>);
static_assert(!CanBakePaletteLab<BakedPalette, BakedPalette>);

/** @brief Owned palette copies rebind Lab nodes and retain external RGB bindings. */
inline void test_pullback_palette_lab_owned_copy() {
  const auto check = []<int Values>() {
    PaletteLabLutStorage<Values> source;
    for (size_t i = 0; i < source.nodes.size(); ++i)
      source.nodes[i] = {float(i) / source.nodes.size(), 0.1f, -0.2f};
    std::array<Pixel, HueRotationLutView::SIZE> rgb{};
    source.active = true;
    source.rgb_nodes = rgb.data();
    source.rgb_midpoint = true;
    source.rgb_midpoint_row = 0;
    const auto NODES = source.nodes;
    PaletteLabLutStorage<Values> copy(source);
    HS_EXPECT_TRUE(copy.data == copy.nodes.data() && copy.data != source.data);
    HS_EXPECT_EQ(copy.value_steps, Values);
    HS_EXPECT_TRUE(copy.active && copy.rgb_midpoint);
    HS_EXPECT_EQ(copy.rgb_midpoint_row, 0);
    HS_EXPECT_TRUE(copy.rgb_nodes == rgb.data());
    HS_EXPECT_EQ(std::memcmp(copy.nodes.data(), NODES.data(), sizeof(NODES)),
                 0);
    source.nodes[0].L = 1.0f;
    HS_EXPECT_EQ(copy.nodes[0].L, NODES[0].L);
    const auto VIEW = copy.view(0.25f);
    HS_EXPECT_TRUE(VIEW.low == copy.nodes.data());
    HS_EXPECT_EQ(VIEW.value_steps, Values);
    HS_EXPECT_TRUE(VIEW.rgb_low == rgb.data() && VIEW.rgb_high == rgb.data());
    PaletteLabLutStorage<Values> assigned;
    HS_EXPECT_TRUE(&(assigned = copy) == &assigned);
    HS_EXPECT_TRUE(assigned.data == assigned.nodes.data());
    HS_EXPECT_TRUE(assigned.data != copy.data);
    HS_EXPECT_EQ(
        std::memcmp(assigned.nodes.data(), NODES.data(), sizeof(NODES)), 0);
    auto *same = &assigned;
    assigned = *same;
    HS_EXPECT_TRUE(assigned.data == assigned.nodes.data());
    HS_EXPECT_EQ(
        std::memcmp(assigned.nodes.data(), NODES.data(), sizeof(NODES)), 0);
    HS_EXPECT_TRUE(assigned.rgb_nodes == rgb.data() && assigned.rgb_midpoint);
    HS_EXPECT_EQ(assigned.rgb_midpoint_row, 0);
    copy.nodes[0].L = 0.75f;
    HS_EXPECT_EQ(assigned.nodes[0].L, NODES[0].L);
  };
  check.template operator()<64>();
  check.template operator()<256>();
}

/** @brief Cached neutral adjustments materialize their base with exact alpha. */
inline void test_pullback_cached_lab_finish() {
  using namespace Pullback;
  ColorSample sample{Color4(Pixel(0, 0, 0), 0.37f),
                     {0.5f, 0.7f, math::Y_AXIS, 2.0f}};
  sample.cached_lab = {0.7f, 0.02f, 0.03f};
  sample.has_cached_lab = true;
  const auto EXPECTED =
      linrgb_to_pixel(oklab_to_linear_rgb_lut_gamut(sample.cached_lab));
  const auto ACTUAL = finish_color_sample(sample);
  HS_EXPECT_EQ(ACTUAL.color, EXPECTED);
  HS_EXPECT_EQ(ACTUAL.alpha, 0.37f);
  HS_EXPECT_GT(ACTUAL.color.r, 0);
  sample.brightness_gain = 0.0f;
  HS_EXPECT_EQ(finish_color_sample(sample).color, Pixel(0, 0, 0));
  HS_EXPECT_EQ(finish_color_sample(sample).alpha, 0.37f);
  HS_EXPECT_TRUE(blend_color_samples(sample, sample, 0.0f).has_cached_lab);
  HS_EXPECT_TRUE(blend_color_samples(sample, sample, 1.0f).has_cached_lab);
}
using TemporalPalettePipeline = Pullback::Pipeline<
    PaletteBinding, PaletteSampleStage,
    Pullback::Stage::Colorize<
        Pullback::Color::GeneratedPalette<TemporalPaletteProvider>>,
    Pullback::Stage::Modulate<Pullback::Color::Modulate<PaletteLayerProvider>>>;
template <Pullback::CodeEmission Emission>
using TemporalPalettePlacedPipeline = Pullback::Pipeline<
    PaletteBinding, PaletteSampleStage,
    Pullback::Stage::Placed<
        Emission,
        Pullback::Stage::Colorize<
            Pullback::Color::GeneratedPalette<TemporalPaletteProvider>>,
        Pullback::Stage::Modulate<
            Pullback::Color::Modulate<PaletteLayerProvider>>>>;

struct PhaseLayerProvider {
  using Binding = PaletteBinding;
  using FrameState = PaletteFrame;
  static Pullback::Color::PhaseModulationState
  prepare(const PaletteFrame &frame) {
    return *frame.phase;
  }
};
struct CountingPhasePolicy : Pullback::Color::Modulate<PhaseLayerProvider> {
  using Base = Pullback::Color::Modulate<PhaseLayerProvider>;
  inline static int phase_calls = 0;
  inline static int endpoint_calls = 0;
  static void apply_endpoint_in_place(Pullback::ColorSample &sample,
                                      const PaletteFrame &frame,
                                      const Prepared &prepared, bool previous) {
    ++endpoint_calls;
    Base::apply_endpoint_in_place(sample, frame, prepared, previous);
  }
  static void apply_phase_in_place(Pullback::ColorSample &sample,
                                   const PaletteFrame &frame,
                                   const Prepared &prepared) {
    ++phase_calls;
    Base::apply_phase_in_place(sample, frame, prepared);
  }
};

struct LabPaletteProvider : TemporalPaletteProvider {
  static const PaletteLabLut *palette_lab_lut(const PaletteFrame &frame) {
    return frame.lab;
  }
};
struct CountingLabPalette
    : Pullback::Color::GeneratedPalette<LabPaletteProvider> {
  using Base = Pullback::Color::GeneratedPalette<LabPaletteProvider>;
  using SharedBaseOwner = CountingLabPalette;
  inline static int shared_calls = 0;
  inline static int pair_calls = 0;
  static bool apply_shared_base(const Pullback::FieldSample &input,
                                const PaletteFrame &frame,
                                const Prepared &prepared,
                                Pullback::ColorSample &output) {
    ++shared_calls;
    return Base::apply_shared_base(input, frame, prepared, output);
  }
  static auto apply_pair(const Pullback::FieldSample &input,
                         const PaletteFrame &frame, const Prepared &prepared) {
    ++pair_calls;
    return Base::apply_pair(input, frame, prepared);
  }
};

struct ShiftedLabPalette : CountingLabPalette {
  static auto apply_pair(const Pullback::FieldSample &input,
                         const PaletteFrame &frame, const Prepared &prepared) {
    auto pair = CountingLabPalette::apply_pair(input, frame, prepared);
    pair.first.color.color.r += 123;
    pair.second.color.color.r += 123;
    return pair;
  }
};

/** @brief RGB-bound fused carriers retain the shared RGB and materialized-Lab finishers. */
inline void test_pullback_fused_rgb_carrier_finish() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  using Pipeline = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage, Stage::Colorize<Policy>,
      Stage::Modulate<Color::Modulate<PaletteLayerProvider>>>;
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(12000, 13000, 14000), 1.0f)},
                ConstantPalette{Color4(Pixel(22000, 23000, 24000), 1.0f)});
  [[maybe_unused]] auto &a = palettes.a;
  [[maybe_unused]] auto &b = palettes.b;
  PaletteLabLut cache;
  cache.bake(a.view(), b.view());
  PaletteRgbHueLut rgb;
  rgb.bake(cache);
  PaletteFrame frame;
  palettes.bind(frame);
  frame.lab = &cache;
  Modulation::Prepared previous{}, current{};
  const auto SOURCE = [](const Modulation::Prepared &p,
                         const Pattern::SphereDirection &) {
    return p.state.phase.phase_turns;
  };
  previous.evaluate = current.evaluate = SOURCE;
  previous.state.phase.phase_turns = 0.13f;
  current.state.phase.phase_turns = 0.67f;
  frame.previous_modulation = &previous;
  frame.modulation = &current;
  const FieldSample FIELD{0.5f, 0.5f, math::Y_AXIS, 0.0f};
  for (float fraction : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f})
    for (auto kind : {ColorModulation::Mechanic::HUE_ROTATION,
                      ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT,
                      ColorModulation::Mechanic::OPACITY}) {
      frame.temporal.fraction = fraction;
      frame.mechanic = {kind, 0.4f};
      const auto PREPARED = Policy::prepare(frame);
      const auto MODULATION = PaletteLayerProvider::prepare(frame);
      const auto CARRIER = Policy::apply_fused_modulation(
          FIELD, frame, PREPARED, MODULATION, fraction);
      if (fraction > 0.0f && fraction < 1.0f) {
        HS_EXPECT_TRUE(CARRIER.cached_rgb == &PREPARED.lab);
        HS_EXPECT_TRUE(CARRIER.has_cached_lab);
        for (float gain : {0.0f, 0.3f, 1.0f})
          for (float lift : {0.0f, 0.4f}) {
            auto sample = CARRIER;
            sample.brightness_gain = gain;
            sample.pending_lift = lift;
            const auto EXPECTED = finish_color_sample(sample);
            int callbacks = 0;
            HS_EXPECT_FALSE(Color::restore_cached_palette_base_near_gamut_gate(
                sample, [&] { ++callbacks; }));
            const auto ACTUAL =
                Color::finish_fused_color(sample, [&] { ++callbacks; });
            HS_EXPECT_EQ(callbacks, 0);
            HS_EXPECT_EQ(ACTUAL.color, EXPECTED.color);
            HS_EXPECT_EQ(ACTUAL.alpha, EXPECTED.alpha);
          }
      }
      const auto EXPECTED = finish_color_sample(CARRIER);
      const auto DIRECT = Policy::shade_fused_modulation(FIELD, frame, PREPARED,
                                                         MODULATION, fraction);
      HS_EXPECT_EQ(DIRECT.color, EXPECTED.color);
      HS_EXPECT_EQ(DIRECT.alpha, EXPECTED.alpha);
      const auto PIXEL =
          Pipeline::shade(math::Y_AXIS, Pipeline::prepare(frame));
      HS_EXPECT_EQ(PIXEL.color, EXPECTED.color);
      HS_EXPECT_EQ(PIXEL.alpha, EXPECTED.alpha);
    }
}

struct CachedSinePhaseProvider : PhaseLayerProvider {
  static constexpr bool SINE_PHASE = true;
};

/** @brief Direct cached phase shading matches deferred finishing and its guards. */
inline void test_pullback_cached_phase_finish() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  using Layer = Color::Modulate<CachedSinePhaseProvider>;
  using Pipeline =
      Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                         Stage::Colorize<Policy>, Stage::Modulate<Layer>>;
  using Core =
      Detail::PipelineCore<true, PaletteBinding, typename Pipeline::NodeList>;
  static_assert(Core::has_direct_final_color());
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(12000, 8000, 4000), 1.0f)},
                ConstantPalette{Color4(Pixel(7000, 15000, 19000), 1.0f)});
  [[maybe_unused]] auto &a = palettes.a;
  [[maybe_unused]] auto &b = palettes.b;
  PaletteLabLut cache;
  cache.bake(a.view(), b.view());
  PaletteRgbHueLut rgb;
  rgb.bake(cache);
  PaletteFrame frame;
  palettes.bind(frame);
  frame.lab = &cache;
  Color::PhaseModulationState modulation;
  modulation.source = Pattern::SourceKind::DISPLACEMENT;
  modulation.enabled = true;
  frame.phase = &modulation;
  modulation.grid.a = {math::Y_AXIS, 2.1f, 0.37f};
  Modulation::Prepared sine{};
  sine.config.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  sine.source = modulation.grid.a;
  sine.evaluate = [](const Modulation::Prepared &prepared,
                     const Pattern::SphereDirection &direction) {
    return Pattern::DirectionalSine::sample(
        direction,
        std::get<Pattern::DirectionalSine::Prepared>(prepared.source));
  };
  modulation.current = modulation.previous = &sine;
  for (auto source : {Pattern::SourceKind::DISPLACEMENT,
                      Pattern::SourceKind::DIRECTIONAL_SINE}) {
    modulation.source = source;
    for (float fraction : {0.0f, 0.125f, 0.5f, 0.875f, 1.0f})
      for (float amount : {-3.7f, -0.2f, 0.0f, 0.65f})
        for (bool ready : {false, true})
          for (auto kind :
               {ColorModulation::Mechanic::HUE_ROTATION,
                ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT}) {
            if (kind == ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT &&
                amount < 0.0f)
              continue;
            modulation.ready = ready;
            modulation.mechanic = {kind, amount};
            frame.temporal.fraction = fraction;
            const auto PREPARED = Pipeline::prepare_stages(frame);
            const SphereSample INPUT{math::Y_AXIS, 0.7f};
            const auto EXPECTED = finish_color_sample(
                Core::template run_stage<0>(INPUT, frame, PREPARED));
            const auto ACTUAL =
                Core::template shade_stage<0>(INPUT, frame, PREPARED);
            HS_EXPECT_EQ(ACTUAL.color, EXPECTED.color);
            HS_EXPECT_EQ(ACTUAL.alpha, EXPECTED.alpha);
            const auto PALETTE = Policy::prepare(frame);
            for (float value : {0.0f, 0.01f, 0.5f, 0.99f, 1.0f}) {
              const FieldSample FIELD{value, 0.5f, math::Y_AXIS, 0.7f};
              ColorSample carrier;
              const bool SHARED =
                  Policy::apply_shared_base(FIELD, frame, PALETTE, carrier);
              Color4 direct;
              const bool FINISHED = Policy::template shade_cached_phase<Layer>(
                  FIELD, PALETTE, modulation, direct);
              HS_EXPECT_EQ(
                  FINISHED,
                  SHARED && ready &&
                      (amount == 0.0f ||
                       kind == ColorModulation::Mechanic::HUE_ROTATION));
              if (FINISHED) {
                Layer::apply_phase_in_place(carrier, frame, modulation);
                const auto REFERENCE = finish_color_sample(carrier);
                HS_EXPECT_EQ(direct.color, REFERENCE.color);
                HS_EXPECT_EQ(direct.alpha, REFERENCE.alpha);
              }
            }
          }
  }
}

/** @brief Consecutive analogous cycle bakes meet cached channel bounds. */
inline void test_pullback_palette_lab_qualification() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  std::array<uint8_t, PaletteCycler::generated_arena_bytes() +
                          2 * BakedPalette::required_arena_bytes()>
      buffer{};
  Arena arena(buffer.data(), buffer.size());
  uint32_t hue = 0;
  PaletteCycler cycler;
  cycler.init_generated(
      arena,
      +[](void *context, uint32_t sequence, GenerativePalette &out) {
        GeneratedPaletteBank::next_palette(*static_cast<uint32_t *>(context),
                                           sequence, PaletteHarmony::ANALOGOUS,
                                           0.84f, out);
      },
      &hue, GeneratedPaletteBank::DWELL_FRAMES,
      GeneratedPaletteBank::FADE_FRAMES, math::ease_in_out_sin, nullptr, 2);
  for (int frame = 0; frame < GeneratedPaletteBank::FADE_FRAMES / 2; ++frame)
    cycler.step();
  BakedPaletteStorage from, to;
  from.clone_from(cycler.palette(), arena);
  cycler.step();
  cycler.step();
  to.clone_from(cycler.palette(), arena);
  PaletteLabLut cache;
  cache.bake(from.view(), to.view());
  HS_EXPECT_TRUE(cache.active);
  PaletteFrame frame;
  frame.palette = &to.view();
  frame.temporal.palette = &from.view();
  frame.lab = &cache;
  uint32_t maximum = 0;
  uint64_t total = 0;
  size_t channels = 0;
  constexpr int VALUES = 2048;
  constexpr int FRACTIONS = 128;
  for (int time = 0; time < FRACTIONS; ++time) {
    frame.temporal.fraction = (static_cast<float>(time) + 0.5f) / FRACTIONS;
    const auto PREPARED = Policy::prepare(frame);
    for (int index = 0; index < VALUES; ++index) {
      const FieldSample FIELD{(static_cast<float>(index) + 0.5f) / VALUES,
                              0.75f, math::Y_AXIS, 0.0f};
      ColorSample cached;
      HS_EXPECT_TRUE(Policy::apply_shared_base(FIELD, frame, PREPARED, cached));
      const auto PAIR = Policy::apply_pair(FIELD, frame, PREPARED);
      const auto EXPECTED = finish_color_sample(blend_color_samples(
          PAIR.first, PAIR.second, frame.temporal.fraction));
      const auto ACTUAL = finish_color_sample(cached);
      HS_EXPECT_NEAR(ACTUAL.alpha, EXPECTED.alpha, 1e-7f);
      for (const uint32_t ERROR :
           {channel_gap(ACTUAL.color.r, EXPECTED.color.r),
            channel_gap(ACTUAL.color.g, EXPECTED.color.g),
            channel_gap(ACTUAL.color.b, EXPECTED.color.b)}) {
        maximum = std::max(maximum, ERROR);
        total += ERROR;
        ++channels;
      }
    }
  }
  const double MEAN = static_cast<double>(total) / channels;
  std::printf("  palette lab cache max=%u mean=%.3f channel codes\n", maximum,
              MEAN);
  HS_EXPECT_GT(maximum, 0u);
  HS_EXPECT_LE(static_cast<float>(maximum), Policy::METRICS[0].limit);
  HS_EXPECT_LE(MEAN, Policy::METRICS[1].limit);
}

/** @brief Displacement source selection for modulation qualification. */
struct DisplacementLayerProvider : PaletteLayerProvider {
  static Pullback::Color::ModulationState prepare(const PaletteFrame &frame) {
    auto prepared = PaletteLayerProvider::prepare(frame);
    prepared.source = Pattern::SourceKind::DISPLACEMENT;
    return prepared;
  }
};

enum class Cycle { TRIADIC, COMPLEMENTARY, ANALOGOUS };

/** @brief Prepared hue modulation meets channel limits across cycling palettes. */
inline void test_pullback_modulation_lut_qualification() {
  using namespace Pullback;
  using Layer = Color::Modulate<DisplacementLayerProvider>;
  std::array<uint8_t, GeneratedPaletteBank::required_arena_bytes()> buffer{};
  Arena arena(buffer.data(), buffer.size());
  GeneratedPaletteBank bank;
  bank.init(arena, 0.788f, math::ease_in_out_sin);
  std::array<Pixel, HueRotationLutView::SIZE> table;
  HueRotationGamutMask mask;
  PaletteFrame frame;
  frame.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 1.0f};
  const auto PREPARED = Layer::prepare(frame);
  uint32_t maximum = 0;
  int worst_frame = 0, worst_value = 0, worst_hue = 0;
  uint64_t total = 0;
  size_t channels = 0;
  int elapsed = 0;
  for (const int target : {0, 1, 150, 300, 600, 900, 1200, 1800}) {
    while (elapsed < target) {
      bank.step(Cycle::TRIADIC);
      ++elapsed;
    }
    const auto &palette = bank.palette(Cycle::TRIADIC);
    prepare_hue_rotation_lut(std::span(table), palette, &mask);
    constexpr int VALUES = 512, HUES = 128;
    for (int value = 0; value < VALUES; ++value) {
      const float COORDINATE = (value + 0.5f) / VALUES;
      const auto BASE = palette.get(COORDINATE);
      for (int hue = 0; hue < HUES; ++hue) {
        const float SCALAR = 2.0f * (hue + 0.5f) / HUES - 1.0f;
        ColorSample input{BASE, {COORDINATE, 0.75f, math::Y_AXIS, SCALAR}};
        input.hue_rotation = table.data();
        input.hue_rotation_gamut = &mask;
        input.palette_coordinate = COORDINATE;
        const auto ACTUAL =
            finish_color_sample(Layer::apply(input, frame, PREPARED));
        const auto EXPECTED =
            ColorModulation::apply(BASE, SCALAR, frame.mechanic);
        HS_EXPECT_EQ(ACTUAL.alpha, EXPECTED.alpha);
        for (const uint32_t ERROR :
             {channel_gap(ACTUAL.color.r, EXPECTED.color.r),
              channel_gap(ACTUAL.color.g, EXPECTED.color.g),
              channel_gap(ACTUAL.color.b, EXPECTED.color.b)}) {
          if (ERROR > maximum) {
            maximum = ERROR;
            worst_frame = target;
            worst_value = value;
            worst_hue = hue;
          }
          total += ERROR;
          ++channels;
        }
      }
    }
  }
  const double MEAN = static_cast<double>(total) / channels;
  std::printf("  modulation LUT max=%u mean=%.3f channels=%zu worst=%d/%d/%d\n",
              maximum, MEAN, channels, worst_frame, worst_value, worst_hue);
  HS_EXPECT_GT(maximum, 0u);
  HS_EXPECT_LE(static_cast<float>(maximum), Layer::METRICS[0].limit);
  HS_EXPECT_LE(MEAN, Layer::METRICS[1].limit);
}

/** @brief Shared cached bases skip both palette endpoints and preserve alpha. */
inline void test_pullback_shared_lab_dispatch() {
  using namespace Pullback;
  using Pipeline = Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                                      Stage::Colorize<CountingLabPalette>,
                                      Stage::Modulate<CountingPhasePolicy>>;
  using Core =
      Detail::PipelineCore<true, PaletteBinding, typename Pipeline::NodeList>;
  using Custom = Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                                    Stage::Colorize<ShiftedLabPalette>,
                                    Stage::Modulate<CountingPhasePolicy>>;
  using CustomCore =
      Detail::PipelineCore<true, PaletteBinding, typename Custom::NodeList>;
  TemporalPalettes<BakedPalette::required_arena_bytes()> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(1000, 2000, 3000), 1.0f)},
                ConstantPalette{Color4(Pixel(9000, 6000, 1000), 1.0f)});
  [[maybe_unused]] auto &a = palettes.a;
  [[maybe_unused]] auto &b = palettes.b;
  PaletteLabLut cache;
  cache.bake(a.view(), b.view());
  Color::PhaseModulationState layer;
  layer.source = Pattern::SourceKind::DISPLACEMENT;
  layer.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 0.2f};
  layer.enabled = true;
  layer.ready = true;
  PaletteFrame frame;
  palettes.bind(frame);
  frame.phase = &layer;
  frame.lab = &cache;
  frame.temporal.fraction = 0.37f;
  const auto DIRECT_PREPARED = CountingLabPalette::Base::prepare(frame);
  const FieldSample DIRECT_FIELD{0.5f, 0.5f, math::Y_AXIS, 0.7f};
  ColorSample dirty{Color4(Pixel(123, 456, 789), 0.9f),
                    {0.1f, 0.2f, math::X_AXIS, 0.3f},
                    reinterpret_cast<const Pixel *>(&cache),
                    0.4f,
                    0.5f,
                    0.6f,
                    0.7f,
                    {0.8f, 0.9f, 1.0f},
                    false};
  const auto STALE_RGB = cache.view(0.0f);
  dirty.cached_rgb = &STALE_RGB;
  HS_EXPECT_TRUE(CountingLabPalette::Base::apply_shared_base(
      DIRECT_FIELD, frame, DIRECT_PREPARED, dirty));
  HS_EXPECT_EQ(dirty.color.color, Pixel(0, 0, 0));
  const auto DIRECT_PAIR = CountingLabPalette::Base::apply_pair(
      DIRECT_FIELD, frame, DIRECT_PREPARED);
  HS_EXPECT_EQ(dirty.color.alpha, DIRECT_PAIR.first.color.alpha);
  HS_EXPECT_EQ(dirty.field.value, DIRECT_FIELD.value);
  HS_EXPECT_EQ(dirty.field.coverage, DIRECT_FIELD.coverage);
  HS_EXPECT_EQ(dirty.field.sphere, DIRECT_FIELD.sphere);
  HS_EXPECT_EQ(dirty.field.path_length, DIRECT_FIELD.path_length);
  HS_EXPECT_TRUE(dirty.hue_rotation == nullptr);
  HS_EXPECT_TRUE(dirty.cached_rgb == nullptr);
  HS_EXPECT_EQ(dirty.palette_coordinate, 0.5f);
  HS_EXPECT_EQ(dirty.pending_hue, 0.0f);
  HS_EXPECT_EQ(dirty.brightness_gain, 1.0f);
  HS_EXPECT_EQ(dirty.pending_lift, 0.0f);
  const auto EXPECTED_LAB = PaletteLabLut::sample(cache.view(0.37f), 0.5f);
  HS_EXPECT_EQ(dirty.cached_lab.L, EXPECTED_LAB.L);
  HS_EXPECT_EQ(dirty.cached_lab.a, EXPECTED_LAB.a);
  HS_EXPECT_EQ(dirty.cached_lab.b, EXPECTED_LAB.b);
  HS_EXPECT_TRUE(dirty.has_cached_lab);
  for (bool ready : {false, true})
    for (bool active : {false, true})
      for (float fraction : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f}) {
        layer.ready = ready;
        frame.lab = active ? &cache : nullptr;
        frame.temporal.fraction = fraction;
        CountingLabPalette::shared_calls = 0;
        CountingLabPalette::pair_calls = 0;
        const auto PREPARED = Pipeline::prepare_stages(frame);
        const auto ACTUAL = Core::template run_stage<0>(
            SphereSample{math::Y_AXIS, 0.7f}, frame, PREPARED);
        const bool INTERIOR = fraction > 0.0f && fraction < 1.0f;
        const bool SHARED = ready && fraction > 0.0f;
        HS_EXPECT_EQ(ACTUAL.has_cached_lab, active && SHARED);
        HS_EXPECT_EQ(CountingLabPalette::shared_calls, SHARED ? 1 : 0);
        HS_EXPECT_EQ(CountingLabPalette::pair_calls,
                     INTERIOR && !(active && ready) ? 1 : 0);
        const FieldSample FIELD{0.5f, 0.5f, math::Y_AXIS, 0.7f};
        const auto PAIR = CountingLabPalette::Base::apply_pair(
            FIELD, frame, std::get<1>(PREPARED));
        HS_EXPECT_EQ(ACTUAL.color.alpha,
                     fraction <= 0.0f ? PAIR.first.color.alpha
                     : fraction >= 1.0f
                         ? PAIR.second.color.alpha
                         : hs::lerp(PAIR.first.color.alpha,
                                    PAIR.second.color.alpha, fraction));
        HS_EXPECT_NEAR(ACTUAL.pending_hue, 0.14f, 1e-7f);
        if (!INTERIOR && !(active && SHARED))
          HS_EXPECT_EQ(ACTUAL.color.color, fraction == 0.0f
                                               ? PAIR.first.color.color
                                               : PAIR.second.color.color);
      }
  layer.ready = true;
  frame.lab = &cache;
  frame.temporal.fraction = 0.5f;
  CountingLabPalette::shared_calls = 0;
  CountingLabPalette::pair_calls = 0;
  const auto CUSTOM = CustomCore::template run_stage<0>(
      SphereSample{math::Y_AXIS, 0.7f}, frame, Custom::prepare_stages(frame));
  HS_EXPECT_FALSE(CUSTOM.has_cached_lab);
  HS_EXPECT_EQ(CountingLabPalette::shared_calls, 0);
  HS_EXPECT_EQ(CountingLabPalette::pair_calls, 1);
  HS_EXPECT_EQ(CUSTOM.color.color.r, 5123);
  BakedPaletteStorage translucent;
  translucent.bake(palettes.arena,
                   ConstantPalette{Color4(Pixel(4000, 5000, 6000), 0.6f)});
  cache.bake(a.view(), translucent.view());
  frame.palette = &translucent.view();
  const auto PREPARED = Pipeline::prepare_stages(frame);
  const auto FALLBACK = Core::template run_stage<0>(
      SphereSample{math::Y_AXIS, 0.7f}, frame, PREPARED);
  const auto PAIR = CountingLabPalette::Base::apply_pair(
      FieldSample{0.5f, 0.5f, math::Y_AXIS, 0.7f}, frame,
      std::get<1>(PREPARED));
  HS_EXPECT_FALSE(FALLBACK.has_cached_lab);
  HS_EXPECT_EQ(FALLBACK.color.alpha,
               hs::lerp(PAIR.first.color.alpha, PAIR.second.color.alpha, 0.5f));
}

/** @brief Shared phase tails run once; incomplete or unready tails keep endpoints. */
inline void test_pullback_phase_tail_dispatch() {
  using namespace Pullback;
  using Phase = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage,
      Stage::Colorize<Color::GeneratedPalette<TemporalPaletteProvider>>,
      Stage::Modulate<CountingPhasePolicy>,
      Stage::Modulate<CountingPhasePolicy>>;
  using Fallback = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage,
      Stage::Colorize<Color::GeneratedPalette<TemporalPaletteProvider>>,
      Stage::Modulate<CountingPhasePolicy>,
      Stage::Modulate<CountingPhasePolicy>,
      Stage::Modulate<Color::Modulate<PaletteLayerProvider>>>;
  using PhaseCore =
      Detail::PipelineCore<true, PaletteBinding, typename Phase::NodeList>;
  using FallbackCore =
      Detail::PipelineCore<true, PaletteBinding, typename Fallback::NodeList>;
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(1000, 2000, 3000), 0.2f)},
                ConstantPalette{Color4(Pixel(9000, 6000, 1000), 0.6f)});
  [[maybe_unused]] auto &a = palettes.a;
  [[maybe_unused]] auto &b = palettes.b;
  Color::PhaseModulationState layer;
  layer.source = Pattern::SourceKind::DISPLACEMENT;
  layer.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 0.2f};
  layer.enabled = true;
  PaletteFrame frame;
  palettes.bind(frame);
  frame.phase = &layer;
  frame.mechanic.amount = 0.0f;
  for (bool ready : {false, true})
    for (float fraction : {0.0f, 0.25f, 0.75f, 1.0f}) {
      layer.ready = ready;
      frame.temporal.fraction = fraction;
      CountingPhasePolicy::phase_calls = 0;
      CountingPhasePolicy::endpoint_calls = 0;
      const auto prepared = Phase::prepare_stages(frame);
      const auto SAMPLE = PhaseCore::template run_stage<0>(
          SphereSample{math::Y_AXIS, 0.7f}, frame, prepared);
      const bool SHARED = ready && fraction > 0.0f && fraction < 1.0f;
      HS_EXPECT_EQ(CountingPhasePolicy::phase_calls, SHARED ? 2 : 0);
      HS_EXPECT_EQ(CountingPhasePolicy::endpoint_calls,
                   SHARED ? 0 : (fraction > 0.0f && fraction < 1.0f ? 4 : 2));
      HS_EXPECT_NEAR(SAMPLE.pending_hue, 0.28f, 1e-7f);
      HS_EXPECT_NEAR(SAMPLE.color.alpha, 0.04f + 0.08f * fraction, 1e-7f);
      HS_EXPECT_NEAR(SAMPLE.color.color.r, 1000 + 8000 * fraction, 1.0f);
      CountingPhasePolicy::phase_calls = 0;
      CountingPhasePolicy::endpoint_calls = 0;
      const auto OTHER = FallbackCore::template run_stage<0>(
          SphereSample{math::Y_AXIS, 0.7f}, frame,
          Fallback::prepare_stages(frame));
      HS_EXPECT_EQ(CountingPhasePolicy::phase_calls, 0);
      HS_EXPECT_EQ(CountingPhasePolicy::endpoint_calls,
                   fraction > 0.0f && fraction < 1.0f ? 4 : 2);
      HS_EXPECT_EQ(OTHER.color.color, SAMPLE.color.color);
      HS_EXPECT_EQ(OTHER.color.alpha, SAMPLE.color.alpha);
      HS_EXPECT_EQ(OTHER.pending_hue, SAMPLE.pending_hue);
    }
}

struct CountingPairedPalette
    : Pullback::Color::GeneratedPalette<TemporalPaletteProvider> {
  using Base = Pullback::Color::GeneratedPalette<TemporalPaletteProvider>;
  inline static int pairs = 0;
  inline static int endpoints = 0;
  static auto apply_pair(const Pullback::FieldSample &sample,
                         const PaletteFrame &frame, const Prepared &prepared) {
    ++pairs;
    return Base::apply_pair(sample, frame, prepared);
  }
  static auto apply_endpoint(const Pullback::FieldSample &sample,
                             const PaletteFrame &frame,
                             const Prepared &prepared, bool previous) {
    ++endpoints;
    return Base::apply_endpoint(sample, frame, prepared, previous);
  }
};
struct CountingEndpointPalette : CountingPairedPalette {
  static std::pair<Pullback::ColorSample, Pullback::ColorSample>
  apply_pair(const Pullback::FieldSample &, const PaletteFrame &,
             const Prepared &) = delete;
};

/** @brief Paired color dispatch skips individual endpoint calls during blends. */
inline void test_pullback_paired_color_dispatch() {
  using namespace Pullback;
  using Paired = Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                                    Stage::Colorize<CountingPairedPalette>>;
  using Fallback = Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                                      Stage::Colorize<CountingEndpointPalette>>;
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(1000, 2000, 3000), 0.2f)},
                ConstantPalette{Color4(Pixel(9000, 6000, 1000), 0.6f)});
  [[maybe_unused]] auto &a = palettes.a;
  [[maybe_unused]] auto &b = palettes.b;
  PaletteFrame frame;
  palettes.bind(frame);
  for (float fraction : {0.0f, 0.25f, 1.0f}) {
    frame.temporal.fraction = fraction;
    CountingPairedPalette::pairs = 0;
    CountingPairedPalette::endpoints = 0;
    const auto COLOR = Paired::shade(math::Y_AXIS, Paired::prepare(frame));
    HS_EXPECT_NEAR(COLOR.color.r, 1000 + 8000 * fraction, 1.0f);
    HS_EXPECT_NEAR(COLOR.alpha, 0.04f + 0.08f * fraction, 1e-6f);
    HS_EXPECT_EQ(CountingPairedPalette::pairs, fraction == 0.25f ? 1 : 0);
    HS_EXPECT_EQ(CountingPairedPalette::endpoints, fraction == 0.25f ? 0 : 1);
    CountingPairedPalette::pairs = 0;
    CountingPairedPalette::endpoints = 0;
    const auto OTHER = Fallback::shade(math::Y_AXIS, Fallback::prepare(frame));
    HS_EXPECT_EQ(OTHER.color, COLOR.color);
    HS_EXPECT_EQ(OTHER.alpha, COLOR.alpha);
    HS_EXPECT_EQ(CountingPairedPalette::pairs, 0);
    HS_EXPECT_EQ(CountingPairedPalette::endpoints,
                 fraction > 0.0f && fraction < 1.0f ? 2 : 1);
  }
}

struct CoordinatePalette {
  float alpha;
  Color4 get(float value) const {
    return Color4(Pixel(value * 60000, (1.0f - value) * 50000, 7000), alpha);
  }
};

/** @brief Paired palettes preserve shared controls and endpoint resources. */
inline void test_pullback_paired_palette_controls() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<TemporalPaletteProvider>;
  TemporalPalettes<> palettes;
  palettes.bake(CoordinatePalette{0.2f}, CoordinatePalette{0.6f});
  [[maybe_unused]] auto &a = palettes.a;
  [[maybe_unused]] auto &b = palettes.b;
  std::array<Pixel, 2> hue{};
  PaletteFrame frame;
  palettes.bind(frame);
  auto prepared = Policy::prepare(frame);
  prepared.current.mapping_frequency = 2.0f;
  prepared.current.mapping_offset = 0.13f;
  prepared.current.brightness_envelope = Color::BrightnessEnvelope::ASCENDING;
  prepared.current.brightness_bottom = 0.2f;
  prepared.current.brightness_top = 0.8f;
  prepared.current.hue_rotation = {hue.data(), true};
  prepared.previous = prepared.current;
  prepared.previous.palette = &a.view();
  prepared.previous.hue_rotation = {hue.data() + 1, true};
  const FieldSample INPUT{0.37f, 0.45f, math::X_AXIS, 2.5f};
  const auto PAIR = Policy::apply_pair(INPUT, frame, prepared);
  for (const auto *sample : {&PAIR.first, &PAIR.second}) {
    HS_EXPECT_NEAR(sample->palette_coordinate, 0.87f, 1e-6f);
    HS_EXPECT_NEAR(sample->brightness_gain, 0.422f, 1e-6f);
    HS_EXPECT_EQ(sample->field.value, INPUT.value);
    HS_EXPECT_EQ(sample->field.coverage, INPUT.coverage);
    HS_EXPECT_EQ(sample->field.sphere, INPUT.sphere);
    HS_EXPECT_EQ(sample->field.path_length, INPUT.path_length);
    HS_EXPECT_EQ(sample->pending_hue, 0.0f);
    HS_EXPECT_EQ(sample->pending_lift, 0.0f);
    HS_EXPECT_EQ(sample->color.color, a.view().get(0.87f).color);
  }
  HS_EXPECT_NEAR(PAIR.first.color.alpha, 0.03132f, 1e-6f);
  HS_EXPECT_NEAR(PAIR.second.color.alpha, 0.09396f, 1e-6f);
  HS_EXPECT_TRUE(PAIR.first.hue_rotation == hue.data() + 1);
  HS_EXPECT_TRUE(PAIR.second.hue_rotation == hue.data());
}

/** @brief The complete color suffix interpolates immutable palette endpoints. */
inline void test_pullback_temporal_palette() {
  using namespace Pullback;
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(1000, 2000, 3000), 0.2f)},
                ConstantPalette{Color4(Pixel(9000, 6000, 1000), 0.6f)});
  [[maybe_unused]] auto &a = palettes.a;
  [[maybe_unused]] auto &b = palettes.b;
  PaletteFrame frame;
  palettes.bind(frame);
  frame.mechanic.amount = 0.0f;
  for (float fraction : {0.0f, 0.25f, 0.75f, 1.0f}) {
    frame.temporal.fraction = fraction;
    const auto COLOR = TemporalPalettePipeline::shade(
        math::Y_AXIS, TemporalPalettePipeline::prepare(frame));
    HS_EXPECT_NEAR(COLOR.color.r, 1000 + 8000 * fraction, 1.0f);
    HS_EXPECT_NEAR(COLOR.color.g, 2000 + 4000 * fraction, 1.0f);
    HS_EXPECT_NEAR(COLOR.color.b, 3000 - 2000 * fraction, 1.0f);
    HS_EXPECT_NEAR(COLOR.alpha, 0.04f + 0.08f * fraction, 1e-6f);
  }
}

/** @brief Opacity modulates each endpoint before final color interpolation. */
inline void test_pullback_modulation_temporal_alpha() {
  using namespace Pullback;
  const Pixel RGB(1234, 4567, 8910);
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(RGB, 0.2f)},
                ConstantPalette{Color4(RGB, 0.6f)});
  [[maybe_unused]] auto &a = palettes.a;
  [[maybe_unused]] auto &b = palettes.b;
  Modulation::Resources<false> previous, current;
  prepare_sine_clock_pair(previous, current);
  PaletteFrame frame;
  palettes.bind(frame);
  frame.modulation = &current.prepared;
  frame.previous_modulation = &previous.prepared;
  frame.mechanic = {ColorModulation::Mechanic::OPACITY, 1.0f};
  for (float fraction : {0.0f, 0.25f, 0.75f, 1.0f}) {
    frame.temporal.fraction = fraction;
    const auto COLOR = TemporalPalettePipeline::shade(
        math::Y_AXIS, TemporalPalettePipeline::prepare(frame));
    HS_EXPECT_EQ(COLOR.color, RGB);
    HS_EXPECT_NEAR(COLOR.alpha, 0.12f * fraction, 1e-6f);
    const auto check_placed = [&]<CodeEmission Emission>() {
      using Placed = TemporalPalettePlacedPipeline<Emission>;
      const auto PLACED = Placed::shade(math::Y_AXIS, Placed::prepare(frame));
      HS_EXPECT_EQ(PLACED.color, COLOR.color);
      HS_EXPECT_NEAR(PLACED.alpha, COLOR.alpha, 0.0f);
    };
    check_placed.template operator()<CodeEmission::INLINE_ONLY>();
    check_placed.template operator()<CodeEmission::OUT_OF_LINE_FLASH>();
    check_placed.template operator()<CodeEmission::OUT_OF_LINE_ITCM>();
  }
  frame.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 1.0f};
  frame.temporal.fraction = 0.5f;
  const auto COLOR = TemporalPalettePipeline::shade(
      math::Y_AXIS, TemporalPalettePipeline::prepare(frame));
  const auto EXPECTED =
      ColorModulation::PositiveLightnessLift::apply(Color4(RGB), 0.5f, {1.0f})
          .color;
  HS_EXPECT_NEAR(COLOR.color.r, EXPECTED.r, 1.0f);
  HS_EXPECT_NEAR(COLOR.color.g, EXPECTED.g, 1.0f);
  HS_EXPECT_NEAR(COLOR.color.b, EXPECTED.b, 1.0f);
}

/** @brief Hue and lightness composition uses one final gamut mapping.
 * @param base Unmodified straight-alpha color.
 * @param turns Combined hue turns.
 * @param lift Combined lift toward white.
 * @param gain Final linear RGB gain.
 * @return The transformed color. */
inline Color4 composed_modulation_reference(Color4 base, float turns,
                                            float lift, float gain = 1.0f) {
  const LinRGB rgb = pixel_to_linrgb(base.color);
  OKLab lab = linear_rgb_to_oklab_fast(rgb.r, rgb.g, rgb.b);
  lab.L += (1.0f - lab.L) * lift;
  float cosine, sine;
  hue_sincos(turns, cosine, sine);
  const OKLab rotated{lab.L, lab.a * cosine - lab.b * sine,
                      lab.a * sine + lab.b * cosine};
  LinRGB output = oklab_to_linear_rgb_lut_gamut(rotated);
  if (gain != 1.0f) {
    output.r = hs::clamp(output.r, 0.0f, 1.0f) * gain;
    output.g = hs::clamp(output.g, 0.0f, 1.0f) * gain;
    output.b = hs::clamp(output.b, 0.0f, 1.0f) * gain;
  }
  return Color4(linrgb_to_pixel(output), base.alpha);
}

/** @brief Perceptual finishing applies brightness before RGB16 quantization. */
inline void test_pullback_finish_brightness_quantization() {
  using namespace Pullback;
  size_t quantization_differences = 0;
  for (uint16_t r : {0, 1, 2, 1000, 32768, 65535})
    for (uint16_t g : {0, 1, 2, 1000, 32768, 65535})
      for (uint16_t b : {0, 1, 2, 1000, 32768, 65535})
        for (float hue : {0.0f, 0.17f, 0.75f})
          for (float lift : {0.0f, 0.25f, 0.9f}) {
            if (hue == 0.0f && lift == 0.0f)
              continue;
            const Color4 BASE(Pixel(r, g, b), 0.37f);
            const auto UNDIMMED =
                composed_modulation_reference(BASE, hue, lift);
            for (float gain : {-1.0f, 0.0f, 0.001f, 0.25f, 1.0f, 2.0f}) {
              ColorSample sample{BASE, {0.5f, 0.7f, math::Y_AXIS, 2.0f}};
              sample.pending_hue = hue;
              sample.pending_lift = lift;
              sample.brightness_gain = gain;
              const auto ACTUAL = finish_color_sample(sample);
              const auto EXPECTED =
                  composed_modulation_reference(BASE, hue, lift, gain);
              HS_EXPECT_EQ(ACTUAL.color, EXPECTED.color);
              HS_EXPECT_EQ(ACTUAL.alpha, BASE.alpha);
              const Pixel DIMMED_AFTER_QUANTIZATION = UNDIMMED.color * gain;
              const int GAP = std::max(
                  {std::abs(int(ACTUAL.color.r) - DIMMED_AFTER_QUANTIZATION.r),
                   std::abs(int(ACTUAL.color.g) - DIMMED_AFTER_QUANTIZATION.g),
                   std::abs(int(ACTUAL.color.b) -
                            DIMMED_AFTER_QUANTIZATION.b)});
              HS_EXPECT_LE(GAP, gain <= 1.0f ? 1 : 2);
              if (gain == 0.0f)
                HS_EXPECT_EQ(ACTUAL.color, Pixel(0, 0, 0));
              if (gain == 1.0f)
                HS_EXPECT_EQ(ACTUAL.color, UNDIMMED.color);
              quantization_differences +=
                  ACTUAL.color != DIMMED_AFTER_QUANTIZATION;
            }
          }
  HS_EXPECT_GT(quantization_differences, size_t{0});
}

/** @brief Hue layers compose across lightness lift without RGB round trips. */
inline void test_pullback_stacked_modulations() {
  using namespace Pullback;
  const Color4 BASE(Pixel(23000, 17000, 11000), 0.7f);
  const ColorSample INPUT{BASE, {0.4f, 0.6f, math::Y_AXIS, 2.5f}};
  Color::ModulationState first{nullptr,
                               nullptr,
                               Pattern::SourceKind::DISPLACEMENT,
                               {ColorModulation::Mechanic::HUE_ROTATION, 0.3f},
                               true};
  auto second = first;
  second.mechanic.amount = -0.3f;
  const auto CANCELLED = Color::finish(
      Color::apply_modulation(Color::apply_modulation(INPUT, first), second));
  HS_EXPECT_EQ(CANCELLED.color, BASE.color);
  HS_EXPECT_NEAR(CANCELLED.alpha, BASE.alpha, 0.0f);
  std::array<Pixel, HueRotationLutView::SIZE> table;
  for (size_t index = 0; index < table.size(); ++index)
    table[index] = Pixel((index % 16) * 1000, 2000, 3000);
  auto with_table = INPUT;
  with_table.hue_rotation = table.data();
  HS_EXPECT_EQ(
      Color::finish(Color::apply_modulation(
                        Color::apply_modulation(with_table, first), second))
          .color,
      BASE.color);
  Modulation::Config lift_config;
  lift_config.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  lift_config.evaluation = Pattern::Evaluation::DIRECT;
  lift_config.sine.frequency = 0.0f;
  auto clock = Modulation::init(lift_config);
  Modulation::Resources<false> lift_source;
  for (float phase : {0.0f, 0.75f}) {
    clock.phase.phase_turns = phase;
    const auto &source = lift_source.prepare(lift_config, clock);
    const Color::ModulationState IDENTITY_LIFT{
        &source,
        &source,
        Pattern::SourceKind::DIRECTIONAL_SINE,
        {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 1.0f},
        true};
    const auto LAYERED = Color::apply_modulation(
        Color::apply_modulation(Color::apply_modulation(with_table, first),
                                IDENTITY_LIFT),
        second);
    HS_EXPECT_EQ(Color::finish(LAYERED).color, BASE.color);
  }
  const auto ROTATED = Color::finish(Color::apply_modulation(INPUT, first));
  const auto EXPECTED = ColorModulation::HueRotation::apply(BASE, 2.5f, {0.3f});
  HS_EXPECT_EQ(ROTATED.color, EXPECTED.color);
  auto lift = first;
  lift.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 0.4f};
  second.mechanic.amount = 0.1f;
  auto actual = Color::apply_modulation(
      Color::apply_modulation(Color::apply_modulation(INPUT, first), lift),
      second);
  HS_EXPECT_EQ(actual.color.color, BASE.color);
  const auto COMPOSED = composed_modulation_reference(BASE, 0.0f, 0.4f);
  HS_EXPECT_EQ(Color::finish(actual).color, COMPOSED.color);
  const auto REORDERED = Color::finish(Color::apply_modulation(
      Color::apply_modulation(Color::apply_modulation(INPUT, first), second),
      lift));
  HS_EXPECT_EQ(REORDERED.color, COMPOSED.color);
  actual.brightness_gain = 0.25f;
  const auto DIMMED = Color::finish(actual).color;
  const auto DIMMED_REFERENCE =
      composed_modulation_reference(BASE, 0.0f, 0.4f, 0.25f).color;
  HS_EXPECT_EQ(DIMMED, DIMMED_REFERENCE);
  auto opacity = first;
  opacity.mechanic = {ColorModulation::Mechanic::OPACITY, 1.0f};
  HS_EXPECT_NEAR(Color::finish(Color::apply_modulation(INPUT, opacity)).alpha,
                 BASE.alpha, 0.0f);
  HS_EXPECT_EQ(
      Color::finish(Color::apply_modulation(INPUT, lift)).color,
      ColorModulation::PositiveLightnessLift::apply(BASE, 1.0f, {0.4f}).color);
}

/** @brief Saturated lifted colors rotate before final clipping and opacity. */
inline void test_pullback_composed_lightness_hue() {
  using namespace Pullback;
  const Color::ModulationState HUE_A{
      nullptr,
      nullptr,
      Pattern::SourceKind::DISPLACEMENT,
      {ColorModulation::Mechanic::HUE_ROTATION, 0.25f},
      true};
  auto hue_b = HUE_A;
  hue_b.mechanic.amount = 0.5f;
  auto lift_a = HUE_A;
  lift_a.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 0.8f};
  auto lift_b = lift_a;
  lift_b.mechanic.amount = 0.5f;
  auto opacity = HUE_A;
  opacity.mechanic = {ColorModulation::Mechanic::OPACITY, 0.5f};
  for (const Pixel BASE : {Pixel(65535, 0, 0), Pixel(0, 65535, 0),
                           Pixel(0, 0, 65535), Pixel(62000, 61000, 60000)}) {
    const ColorSample INPUT{Color4(BASE, 0.7f),
                            {0.5f, 1.0f, math::Y_AXIS, 0.5f}};
    auto actual = Color::apply_modulation(INPUT, HUE_A);
    actual = Color::apply_modulation(actual, lift_a);
    actual = Color::apply_modulation(actual, hue_b);
    actual = Color::apply_modulation(actual, lift_b);
    HS_EXPECT_EQ(actual.color.color, BASE);
    const auto EXPECTED =
        composed_modulation_reference(INPUT.color, 0.375f, 0.55f);
    const auto FINISHED = Color::finish(actual);
    HS_EXPECT_NEAR(FINISHED.color.r, EXPECTED.color.r, 1.0f);
    HS_EXPECT_NEAR(FINISHED.color.g, EXPECTED.color.g, 1.0f);
    HS_EXPECT_NEAR(FINISHED.color.b, EXPECTED.color.b, 1.0f);
    HS_EXPECT_NEAR(FINISHED.alpha, 0.7f, 0.0f);
    actual = Color::apply_modulation(actual, opacity);
    actual.brightness_gain = 0.25f;
    const auto MASKED = Color::finish(actual);
    HS_EXPECT_EQ(
        MASKED.color,
        composed_modulation_reference(INPUT.color, 0.375f, 0.55f, 0.25f).color);
    HS_EXPECT_NEAR(MASKED.alpha, 0.6125f, 1e-7f);
    actual.brightness_gain = 0.0f;
    const auto DARK = Color::finish(actual);
    HS_EXPECT_EQ(DARK.color, Pixel(0, 0, 0));
    HS_EXPECT_EQ(DARK.alpha, MASKED.alpha);
  }
}

/** @brief Temporal deferred color state blends before one perceptual finish. */
inline void test_pullback_temporal_composed_modulations() {
  using namespace Pullback;
  using Pipeline = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage,
      Stage::Colorize<Color::GeneratedPalette<TemporalPaletteProvider>>,
      Stage::Modulate<Color::Modulate<PaletteLayerProvider>>,
      Stage::Modulate<Color::Modulate<PaletteHueLayerProvider>>>;
  const Color4 START(Pixel(65535, 0, 0), 0.2f);
  const Color4 END(Pixel(0, 0, 65535), 0.6f);
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{START}, ConstantPalette{END});
  [[maybe_unused]] auto &a = palettes.a;
  [[maybe_unused]] auto &b = palettes.b;
  Modulation::Resources<false> previous, current;
  prepare_sine_clock_pair(previous, current);
  PaletteFrame frame;
  palettes.bind(frame);
  frame.modulation = &current.prepared;
  frame.previous_modulation = &previous.prepared;
  frame.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 0.4f};
  const auto PREVIOUS = composed_modulation_reference(START, -0.25f, 0.0f);
  const auto CURRENT = composed_modulation_reference(END, 0.25f, 0.4f);
  for (float fraction : {0.0f, 0.25f, 0.75f, 1.0f}) {
    frame.temporal.fraction = fraction;
    const auto COLOR = Pipeline::shade(math::Y_AXIS, Pipeline::prepare(frame));
    const auto BASE = START.color.lerp16(END.color, frac_to_q16(fraction));
    const auto EXPECTED = composed_modulation_reference(
        Color4(BASE), 0.75f - 0.5f * fraction, 0.4f * fraction);
    HS_EXPECT_NEAR(COLOR.color.r, EXPECTED.color.r, 1.0f);
    HS_EXPECT_NEAR(COLOR.color.g, EXPECTED.color.g, 1.0f);
    HS_EXPECT_NEAR(COLOR.color.b, EXPECTED.color.b, 1.0f);
    HS_EXPECT_NEAR(COLOR.alpha, 0.04f + 0.08f * fraction, 1e-6f);
    if (fraction > 0.0f && fraction < 1.0f) {
      const auto COMPLETED_BLEND =
          PREVIOUS.color.lerp16(CURRENT.color, frac_to_q16(fraction));
      HS_EXPECT_GT(std::abs(int(COLOR.color.r) - int(COMPLETED_BLEND.r)) +
                       std::abs(int(COLOR.color.g) - int(COMPLETED_BLEND.g)) +
                       std::abs(int(COLOR.color.b) - int(COMPLETED_BLEND.b)),
                   1000);
    }
  }
}

/** @brief Deferred hues take the short arc and preserve endpoint identities. */
inline void test_pullback_temporal_deferred_hue_wrap() {
  using namespace Pullback;
  ColorSample a{Color4(Pixel(65535, 0, 0), 0.2f), {}};
  ColorSample b{Color4(Pixel(0, 0, 65535), 0.8f), {}};
  a.pending_hue = 0.95f;
  b.pending_hue = 0.05f;
  a.pending_lift = 0.2f;
  b.pending_lift = 0.8f;
  a.brightness_gain = 0.25f;
  b.brightness_gain = 0.75f;
  const auto MIXED = blend_color_samples(a, b, 0.5f);
  HS_EXPECT_NEAR(MIXED.pending_hue, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(MIXED.pending_lift, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(MIXED.brightness_gain, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(MIXED.color.alpha, 0.5f, 1e-6f);
  HS_EXPECT_EQ(MIXED.color.color,
               a.color.color.lerp16(b.color.color, frac_to_q16(0.5f)));
  HS_EXPECT_EQ(blend_color_samples(a, b, 0.0f).color.color, a.color.color);
  HS_EXPECT_EQ(blend_color_samples(a, b, 1.0f).color.color, b.color.color);
  HS_EXPECT_EQ(blend_color_samples(a, b, 0.0f).pending_hue, a.pending_hue);
  HS_EXPECT_EQ(blend_color_samples(a, b, 1.0f).pending_hue, b.pending_hue);
}

/** @brief In-place temporal blending preserves current sampling coordinates. */
inline void test_pullback_temporal_deferred_in_place() {
  using namespace Pullback;
  std::array<Pixel, 2> hue{};
  ColorSample a{Color4(Pixel(1000, 2000, 4000), 0.25f),
                {0.1f, 0.2f, math::X_AXIS, 0.3f}};
  ColorSample b{Color4(Pixel(9000, 6000, 12000), 0.75f),
                {0.6f, 0.7f, math::Y_AXIS, 0.8f}};
  a.hue_rotation = hue.data();
  b.hue_rotation = hue.data();
  a.palette_coordinate = 0.2f;
  b.palette_coordinate = 0.9f;
  a.pending_hue = 0.875f;
  b.pending_hue = 0.125f;
  a.pending_lift = 0.25f;
  b.pending_lift = 0.75f;
  a.brightness_gain = 0.5f;
  b.brightness_gain = 0.75f;
  blend_color_samples_in_place(a, b, 0.5f);
  HS_EXPECT_EQ(b.color.color, Pixel(5000, 4000, 8000));
  HS_EXPECT_EQ(b.color.alpha, 0.5f);
  HS_EXPECT_EQ(b.pending_hue, 0.0f);
  HS_EXPECT_EQ(b.pending_lift, 0.5f);
  HS_EXPECT_EQ(b.brightness_gain, 0.625f);
  HS_EXPECT_EQ(b.field.value, 0.6f);
  HS_EXPECT_EQ(b.field.coverage, 0.7f);
  HS_EXPECT_EQ(b.field.sphere, math::Y_AXIS);
  HS_EXPECT_EQ(b.field.path_length, 0.8f);
  HS_EXPECT_EQ(b.palette_coordinate, 0.9f);
  HS_EXPECT_TRUE(b.hue_rotation == nullptr);
  HS_EXPECT_EQ(a.color.color, Pixel(1000, 2000, 4000));
  HS_EXPECT_EQ(a.pending_hue, 0.875f);
  HS_EXPECT_TRUE(a.hue_rotation == hue.data());

  b = a;
  b.color.color.r = 1001;
  blend_color_samples_in_place(a, b, 0.000001f);
  HS_EXPECT_EQ(b.color.color, a.color.color);
  HS_EXPECT_TRUE(b.hue_rotation == nullptr);
  a.cached_lab = {0.4f, 0.1f, -0.2f};
  a.has_cached_lab = true;
  blend_color_samples_in_place(a, b, 0.0f);
  HS_EXPECT_TRUE(b.has_cached_lab);
  HS_EXPECT_EQ(b.cached_lab.L, 0.4f);
  HS_EXPECT_EQ(b.cached_lab.a, 0.1f);
  HS_EXPECT_EQ(b.cached_lab.b, -0.2f);
  HS_EXPECT_TRUE(b.hue_rotation == hue.data());
  blend_color_samples_in_place(a, b, 1.0f);
  HS_EXPECT_TRUE(b.has_cached_lab);
}

/** @brief A fused final layer matches the canonical temporal endpoint pipeline. */
inline void test_pullback_fused_modulation_endpoints() {
  using namespace Pullback;
  const auto expect_alpha = [](float actual, float expected) {
#if HS_TEST_FAST_MATH
    const float SCALE = std::max(fabsf(actual), fabsf(expected));
    HS_EXPECT_NEAR(actual, expected,
                   4.0f * std::numeric_limits<float>::epsilon() * SCALE);
#else
    HS_EXPECT_EQ(actual, expected);
#endif
  };
  TemporalPalettes<> palettes;
  palettes.bake(CoordinatePalette{0.2f},
                ConstantPalette{Color4(Pixel(32000, 44000, 13000), 0.8f)});
  [[maybe_unused]] auto &a = palettes.a;
  [[maybe_unused]] auto &b = palettes.b;
  PaletteFrame frame;
  palettes.bind(frame);
  using Policy = Color::GeneratedPalette<TemporalPaletteProvider>;
  struct CanonicalPalette : Policy {};
  using Actual = TemporalPalettePipeline;
  using Canonical = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage, Stage::Colorize<CanonicalPalette>,
      Stage::Modulate<Color::Modulate<PaletteLayerProvider>>>;
  using ActualCore =
      Detail::PipelineCore<true, PaletteBinding, typename Actual::NodeList>;
  using CanonicalCore =
      Detail::PipelineCore<true, PaletteBinding, typename Canonical::NodeList>;
  using Phase =
      Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                         Stage::Colorize<Policy>,
                         Stage::Modulate<Color::Modulate<PhaseLayerProvider>>>;
  using PhaseCore =
      Detail::PipelineCore<true, PaletteBinding, typename Phase::NodeList>;
  static_assert(
      ActualCore::template fused_color_tail_supported<1, FieldSample>());
  static_assert(
      !CanonicalCore::template fused_color_tail_supported<1, FieldSample>());
  static_assert(
      !PhaseCore::template fused_color_tail_supported<1, FieldSample>());

  std::array<int8_t, Pattern::SphereLutView::SIZE> ac, bc;
  for (size_t i = 0; i < ac.size(); ++i) {
    ac[i] = int(i * 37 % 255) - 127;
    bc[i] = int((i * 83 + 41) % 255) - 127;
  }
  Modulation::Prepared pa{}, pb{};
  pa.cube = {ac.data(), true};
  pb.cube = {bc.data(), true};
  pa.evaluate = pb.evaluate = Modulation::sample_cube;
  Color::ModulationState mod{&pb, &pa, Pattern::SourceKind::NOISE, {}, true};
  frame.modulation = &pb;
  frame.previous_modulation = &pa;
  for (auto mapping :
       {Color::PaletteMapping::LINEAR, Color::PaletteMapping::CUP,
        Color::PaletteMapping::BELL, Color::PaletteMapping::REVERSE})
    for (auto mechanic : {ColorModulation::Mechanic::HUE_ROTATION,
                          ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT,
                          ColorModulation::Mechanic::OPACITY})
      for (float amount : {0.0f, 0.3f, 1.0f, 1.5958333f})
        for (float fraction :
             {0.0f, 0.000001f, 0.01f, 0.25f, 0.5f, 0.75f, 0.99999f, 1.0f})
          for (int v = 0; v <= 64; ++v) {
            if (amount > 1.0f &&
                mechanic != ColorModulation::Mechanic::HUE_ROTATION)
              continue;
            frame.temporal.fraction = fraction;
            auto prepared = Policy::prepare(frame);
            prepared.current.mapping =
                Color::PaletteMappingWeights::single(mapping);
            prepared.current.mapping_frequency = 2.0f;
            prepared.current.mapping_offset = 0.13f;
            prepared.current.brightness_envelope =
                Color::BrightnessEnvelope::ASCENDING;
            prepared.current.brightness_bottom = 0.2f;
            prepared.current.brightness_top = 0.8f;
            prepared.previous = prepared.current;
            prepared.previous.palette = &a.view();
            mod.mechanic = {mechanic, amount};
            frame.mechanic = mod.mechanic;
            const FieldSample field{v / 64.0f, 0.45f,
                                    math::Vector(v / 64.0f, 0.3f, 1.0f), 0.7f};
            auto pair = Policy::apply_pair(field, frame, prepared);
            Color::apply_modulation_in_place(pair.first, mod, true);
            Color::apply_modulation_in_place(pair.second, mod, false);
            const auto expected =
                blend_color_samples(pair.first, pair.second, fraction);
            const auto actual = Policy::apply_fused_modulation(
                field, frame, prepared, mod, fraction);
            HS_EXPECT_EQ(actual.color.color, expected.color.color);
            expect_alpha(actual.color.alpha, expected.color.alpha);
            HS_EXPECT_EQ(actual.pending_hue, expected.pending_hue);
            HS_EXPECT_EQ(actual.pending_lift, expected.pending_lift);
            HS_EXPECT_EQ(actual.brightness_gain, expected.brightness_gain);
            HS_EXPECT_TRUE(actual.hue_rotation == expected.hue_rotation);
            const auto x = Color::finish(actual), y = Color::finish(expected);
            HS_EXPECT_EQ(x.color, y.color);
            expect_alpha(x.alpha, y.alpha);
            const auto DIRECT = Policy::shade_fused_modulation(
                field, frame, prepared, mod, fraction);
            HS_EXPECT_EQ(DIRECT.color, y.color);
            expect_alpha(DIRECT.alpha, y.alpha);
            const auto PIXEL =
                Actual::shade(field.sphere, Actual::prepare(frame));
            const auto REFERENCE =
                Canonical::shade(field.sphere, Canonical::prepare(frame));
            HS_EXPECT_EQ(PIXEL.color, REFERENCE.color);
            expect_alpha(PIXEL.alpha, REFERENCE.alpha);
          }
}

/** @brief Unit-hue wrapping matches canonical periodic interpolation bit for bit. */
inline void test_pullback_unit_hue_blending() {
  const std::array<float, 12> VALUES{-0.0f,
                                     0.0f,
                                     std::numeric_limits<float>::denorm_min(),
                                     std::numeric_limits<float>::min(),
                                     0.125f,
                                     0.25f,
                                     std::nextafter(0.5f, 0.0f),
                                     0.5f,
                                     std::nextafter(0.5f, 1.0f),
                                     0.75f,
                                     0.875f,
                                     std::nextafter(1.0f, 0.0f)};
  const std::array<float, 12> FRACTIONS{
      -1.0f,
      -0.0f,
      0.0f,
      std::numeric_limits<float>::denorm_min(),
      std::numeric_limits<float>::min(),
      0.000001f,
      0.25f,
      0.5f,
      0.75f,
      std::nextafter(1.0f, 0.0f),
      1.0f,
      2.0f};
  const auto check = [](float a, float b, float fraction) {
    const float EXPECTED = interp::shortest_periodic(a, b, fraction, 1.0f);
    const float ACTUAL = Pullback::Color::blend_unit_hue(a, b, fraction);
    HS_EXPECT_EQ(std::bit_cast<uint32_t>(ACTUAL),
                 std::bit_cast<uint32_t>(EXPECTED));
  };
  for (float a : VALUES)
    for (float b : VALUES)
      for (float fraction : FRACTIONS)
        check(a, b, fraction);
  uint32_t random = 0x6d2b79f5u;
  const auto next = [&] {
    random = random * 1664525u + 1013904223u;
    return float(random >> 8) * (1.0f / 16777216.0f);
  };
  for (int i = 0; i < 131072; ++i) {
    const float A = next();
    const float B = next();
    check(A, B, i % 8 == 0 ? FRACTIONS[(i / 8) % FRACTIONS.size()] : next());
  }
}

/** @brief Fused cached bases retain canonical hue, opacity and inactive fallbacks. */
inline void test_pullback_fused_cached_modulation() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  TemporalPalettes<BakedPalette::required_arena_bytes()> palettes;
  palettes.bake(CoordinatePalette{1.0f},
                ConstantPalette{Color4(Pixel(32000, 44000, 13000), 1.0f)});
  [[maybe_unused]] auto &a = palettes.a;
  [[maybe_unused]] auto &b = palettes.b;
  auto &arena = palettes.arena;
  BakedPaletteStorage translucent;
  translucent.bake(arena,
                   ConstantPalette{Color4(Pixel(32000, 44000, 13000), 0.6f)});
  PaletteLabLut cache;
  cache.bake(a.view(), b.view());
  PaletteFrame frame;
  palettes.bind(frame);
  Color::ModulationState modulation{
      nullptr,
      nullptr,
      Pattern::SourceKind::DISPLACEMENT,
      {ColorModulation::Mechanic::HUE_ROTATION, 1.5958333f},
      true};
  int cache_hits = 0;
  for (bool active : {false, true})
    for (float fraction : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f})
      for (float value : {0.0f, 0.13f, 0.5f, 0.99f, 1.0f}) {
        frame.lab = active ? &cache : nullptr;
        frame.temporal.fraction = fraction;
        const auto PREPARED = Policy::prepare(frame);
        const FieldSample FIELD{value, 0.45f, math::Y_AXIS, 0.37f};
        auto pair = Policy::apply_pair(FIELD, frame, PREPARED);
        Color::apply_modulation_in_place(pair.first, modulation, true);
        Color::apply_modulation_in_place(pair.second, modulation, false);
        const auto CANONICAL =
            blend_color_samples(pair.first, pair.second, fraction);
        const auto ACTUAL = Policy::apply_fused_modulation(
            FIELD, frame, PREPARED, modulation, fraction);
        auto expected = CANONICAL;
        if (ACTUAL.has_cached_lab) {
          ++cache_hits;
          HS_EXPECT_TRUE(active && fraction > 0.0f && fraction < 1.0f);
          HS_EXPECT_TRUE(
              Policy::apply_shared_base(FIELD, frame, PREPARED, expected));
          expected.pending_hue = CANONICAL.pending_hue;
          expected.pending_lift = CANONICAL.pending_lift;
          expected.color.alpha = CANONICAL.color.alpha;
          expected.brightness_gain = CANONICAL.brightness_gain;
        }
        HS_EXPECT_EQ(ACTUAL.pending_hue, CANONICAL.pending_hue);
        HS_EXPECT_EQ(ACTUAL.pending_lift, CANONICAL.pending_lift);
        HS_EXPECT_EQ(ACTUAL.color.alpha, CANONICAL.color.alpha);
        HS_EXPECT_EQ(ACTUAL.brightness_gain, CANONICAL.brightness_gain);
        const auto FINISHED = Color::finish(ACTUAL);
        const auto EXPECTED = Color::finish(expected);
        HS_EXPECT_EQ(FINISHED.color, EXPECTED.color);
        HS_EXPECT_EQ(FINISHED.alpha, EXPECTED.alpha);
        const auto DIRECT = Policy::shade_fused_modulation(
            FIELD, frame, PREPARED, modulation, fraction);
        HS_EXPECT_EQ(DIRECT.color, FINISHED.color);
        HS_EXPECT_EQ(DIRECT.alpha, FINISHED.alpha);
      }
  HS_EXPECT_GT(cache_hits, 0);
  cache.bake(a.view(), translucent.view());
  frame.palette = &translucent.view();
  frame.lab = &cache;
  frame.temporal.fraction = 0.5f;
  const auto PREPARED = Policy::prepare(frame);
  const auto FALLBACK = Policy::apply_fused_modulation(
      FieldSample{0.5f, 0.45f, math::Y_AXIS, 0.37f}, frame, PREPARED,
      modulation, 0.5f);
  HS_EXPECT_FALSE(FALLBACK.has_cached_lab);
}

/** @brief Direct Lab shading reuses endpoint taps through canonical gamut fallback. */
inline void test_pullback_fused_lab_direct_fallback() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(1626, 1306, 43097), 1.0f)},
                ConstantPalette{Color4(Pixel(1451, 1516, 42130), 1.0f)});
  [[maybe_unused]] auto &a = palettes.a;
  [[maybe_unused]] auto &b = palettes.b;
  PaletteLabLut cache;
  cache.active = true;
  std::array<Pixel, HueRotationLutView::SIZE> hue;
  hue.fill(Pixel(111, 222, 333));
  PaletteFrame frame;
  frame.palette = &b.view();
  frame.lab = &cache;
  frame.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 1.0f};
  static int calls;
  const auto SOURCE = [](const Modulation::Prepared &,
                         const Pattern::SphereDirection &) {
    ++calls;
    return 0.9814453125f;
  };
  Modulation::Prepared previous{}, current{};
  previous.evaluate = current.evaluate = SOURCE;
  frame.previous_modulation = &previous;
  frame.modulation = &current;
  const auto MODULATION = PaletteLayerProvider::prepare(frame);
  for (bool guard : {false, true}) {
    cache.nodes.fill(guard ? OKLab{0.431782782f, 0.00217834837f, -0.24565056f}
                           : linear_rgb_to_oklab_fast(0.3f, 0.4f, 0.5f));
    for (bool same_base : {false, true}) {
      frame.temporal.palette = same_base ? &b.view() : &a.view();
      frame.hue_rotation = frame.temporal.hue_rotation =
          same_base ? hue.data() : nullptr;
      for (float fraction : {0.25f, 0.5f, 0.75f}) {
        frame.temporal.fraction = fraction;
        const auto PREPARED = Policy::prepare(frame);
        const FieldSample FIELD{0.5f, 0.61689651f, math::Y_AXIS, 0.37f};
        const auto CARRIER = Policy::apply_fused_modulation(
            FIELD, frame, PREPARED, MODULATION, fraction);
        HS_EXPECT_EQ(CARRIER.has_cached_lab, !guard);
        const auto EXPECTED = finish_color_sample(CARRIER);
        calls = 0;
        const auto ACTUAL = Policy::shade_fused_modulation(
            FIELD, frame, PREPARED, MODULATION, fraction);
        HS_EXPECT_EQ(ACTUAL.color, EXPECTED.color);
        HS_EXPECT_EQ(ACTUAL.alpha, EXPECTED.alpha);
        HS_EXPECT_EQ(calls, 2);
        if (guard && same_base)
          HS_EXPECT_EQ(ACTUAL.color, hue.front());
      }
    }
  }
}

/** @brief Cached current endpoints sample one live source and retain exact alpha. */
inline void test_pullback_fused_lab_current_endpoint() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  using Pipeline = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage, Stage::Colorize<Policy>,
      Stage::Modulate<Color::Modulate<PaletteLayerProvider>>>;
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(10000, 40000, 20000), 1.0f)},
                ConstantPalette{Color4(Pixel(1451, 1516, 42130), 1.0f)});
  [[maybe_unused]] auto &a = palettes.a;
  [[maybe_unused]] auto &b = palettes.b;
  PaletteLabLut cache;
  PaletteFrame frame;
  palettes.bind(frame);
  frame.temporal.fraction = 1.0f;
  frame.lab = &cache;
  frame.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 1.0f};
  static int calls;
  const auto SOURCE = [](const Modulation::Prepared &prepared,
                         const Pattern::SphereDirection &) {
    ++calls;
    return prepared.state.phase.phase_turns;
  };
  Modulation::Prepared current{};
  current.evaluate = SOURCE;
  frame.modulation = &current;
  for (bool guard : {false, true}) {
    cache.bake(a.view(), b.view());
    if (guard)
      cache.nodes.fill({0.431782782f, 0.00217834837f, -0.24565056f});
    for (float scalar : {-0.125f, 0.0f, 0.13f, 0.9814453125f}) {
      current.state.phase.phase_turns = scalar;
      const auto MODULATION = PaletteLayerProvider::prepare(frame);
      const auto PREPARED = Policy::prepare(frame);
      for (float coverage : {0.0f, 0.5f, 1.0f}) {
        const FieldSample FIELD{0.5f, coverage, math::Y_AXIS, 0.0f};
        const auto REFERENCE = Policy::apply_fused_modulation(
            FIELD, frame, PREPARED, MODULATION, 1.0f);
        const auto EXPECTED = finish_color_sample(REFERENCE);
        auto cached = REFERENCE;
        cached.cached_lab =
            PaletteLabLut::sample(PREPARED.lab, cached.palette_coordinate);
        cached.has_cached_lab = true;
        Color::restore_cached_palette_base_near_gamut_gate(cached, [] {});
        if (guard && scalar == 0.9814453125f)
          HS_EXPECT_FALSE(cached.has_cached_lab);
        const auto CACHED = finish_color_sample(cached);
        calls = 0;
        const auto ACTUAL = Policy::shade_fused_modulation(
            FIELD, frame, PREPARED, MODULATION, 1.0f);
        HS_EXPECT_EQ(ACTUAL.color, CACHED.color);
        HS_EXPECT_EQ(ACTUAL.alpha, EXPECTED.alpha);
        HS_EXPECT_EQ(calls, 1);
        if (!guard) {
          HS_EXPECT_NEAR(ACTUAL.color.r, EXPECTED.color.r, 2);
          HS_EXPECT_NEAR(ACTUAL.color.g, EXPECTED.color.g, 2);
          HS_EXPECT_NEAR(ACTUAL.color.b, EXPECTED.color.b, 2);
        }
      }
      calls = 0;
      const auto PIXEL =
          Pipeline::shade(math::Y_AXIS, Pipeline::prepare(frame));
      HS_EXPECT_EQ(calls, 1);
      const auto DIRECT = Policy::shade_fused_modulation(
          FieldSample{0.5f, 0.5f, math::Y_AXIS, 0.0f}, frame, PREPARED,
          MODULATION, 1.0f);
      HS_EXPECT_EQ(PIXEL.color, DIRECT.color);
      HS_EXPECT_EQ(PIXEL.alpha, DIRECT.alpha);
    }
  }
  frame.temporal.fraction = 0.0f;
  frame.previous_modulation = &current;
  const auto PREPARED = Policy::prepare(frame);
  const auto MODULATION = PaletteLayerProvider::prepare(frame);
  const FieldSample FIELD{0.5f, 0.5f, math::Y_AXIS, 0.0f};
  const auto EXPECTED = finish_color_sample(
      Policy::apply_fused_modulation(FIELD, frame, PREPARED, MODULATION, 0.0f));
  calls = 0;
  const auto ACTUAL =
      Policy::shade_fused_modulation(FIELD, frame, PREPARED, MODULATION, 0.0f);
  HS_EXPECT_EQ(ACTUAL.color, EXPECTED.color);
  HS_EXPECT_EQ(ACTUAL.alpha, EXPECTED.alpha);
  HS_EXPECT_EQ(calls, 1);
}

/** @brief Gamut-gate crossings restore RGB-first color without normal cache fetches. */
inline void test_pullback_cached_gamut_gate_fallback() {
  using namespace Pullback;
  const float NAN_VALUE = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_TRUE(
      Color::near_palette_gamut_gate({NAN_VALUE, NAN_VALUE, NAN_VALUE}));
  ColorSample canonical{Color4(Pixel(1582, 1359, 42855), 0.61689651f), {}};
  canonical.pending_hue = 0.9814453125f;
  auto cached = canonical;
  cached.color.color = Pixel(0, 0, 0);
  cached.cached_lab = {0.431782782f, 0.00217834837f, -0.24565056f};
  cached.has_cached_lab = true;
  const auto REFERENCE = finish_color_sample(canonical);
  const auto UNGUARDED = finish_color_sample(cached);
  auto direct = cached;
  int direct_fetches = 0;
  const auto DIRECT = Color::finish_fused_color(direct, [&] {
    ++direct_fetches;
    direct.color.color = canonical.color.color;
    direct.hue_rotation = canonical.hue_rotation;
  });
  HS_EXPECT_EQ(direct_fetches, 1);
  HS_EXPECT_FALSE(direct.has_cached_lab);
  HS_EXPECT_EQ(DIRECT.color, REFERENCE.color);
  HS_EXPECT_EQ(DIRECT.alpha, REFERENCE.alpha);
  const Pixel EXPECTED = REFERENCE.color * REFERENCE.alpha;
  const Pixel WRONG = UNGUARDED.color * UNGUARDED.alpha;
  HS_EXPECT_GT(std::max({std::abs(int(EXPECTED.r) - WRONG.r),
                         std::abs(int(EXPECTED.g) - WRONG.g),
                         std::abs(int(EXPECTED.b) - WRONG.b)}),
               4000);
  int fetches = 0;
  HS_EXPECT_TRUE(
      Color::restore_cached_palette_base_near_gamut_gate(cached, [&] {
        ++fetches;
        cached.color.color = canonical.color.color;
        cached.hue_rotation = canonical.hue_rotation;
      }));
  HS_EXPECT_EQ(fetches, 1);
  HS_EXPECT_FALSE(cached.has_cached_lab);
  const auto GUARDED = finish_color_sample(cached);
  HS_EXPECT_EQ(GUARDED.color, REFERENCE.color);
  HS_EXPECT_EQ(GUARDED.alpha, REFERENCE.alpha);
  cached.cached_lab = linear_rgb_to_oklab_fast(0.3f, 0.4f, 0.5f);
  cached.has_cached_lab = true;
  cached.pending_hue = 0.0f;
  fetches = 0;
  HS_EXPECT_FALSE(Color::restore_cached_palette_base_near_gamut_gate(
      cached, [&] { ++fetches; }));
  HS_EXPECT_EQ(fetches, 0);
  HS_EXPECT_TRUE(cached.has_cached_lab);
  for (float hue : {0.0f, 0.13f, 0.5f, 0.9814453125f})
    for (float lift : {0.0f, 0.4f, 1.0f})
      for (float gain : {0.0f, 0.3f, 1.0f})
        for (auto lab : {linear_rgb_to_oklab_fast(0.3f, 0.4f, 0.5f),
                         OKLab{0.5f, 0.3f, -0.3f}}) {
          cached.cached_lab = lab;
          cached.pending_hue = hue;
          cached.pending_lift = lift;
          cached.brightness_gain = gain;
          cached.has_cached_lab = true;
          auto expected = cached;
          int reference_fetches = 0;
          Color::restore_cached_palette_base_near_gamut_gate(expected, [&] {
            ++reference_fetches;
            expected.color.color = canonical.color.color;
          });
          const auto FINISHED = finish_color_sample(expected);
          int actual_fetches = 0;
          const auto FAST = Color::finish_fused_color(cached, [&] {
            ++actual_fetches;
            cached.color.color = canonical.color.color;
          });
          HS_EXPECT_EQ(FAST.color, FINISHED.color);
          HS_EXPECT_EQ(FAST.alpha, FINISHED.alpha);
          HS_EXPECT_EQ(actual_fetches, reference_fetches);
          HS_EXPECT_EQ(cached.has_cached_lab, expected.has_cached_lab);
        }
}

/** @brief Prepared cube-pair eligibility preserves evaluator and active-source boundaries. */
inline void test_pullback_prepared_cube_pair() {
  using namespace Pullback;
  using Layer = Color::Modulate<PaletteLayerProvider>;
  using Policy = Color::GeneratedPalette<TemporalPaletteProvider>;
  constexpr size_t ALIGNMENT = alignof(Color::ModulationState);
  constexpr size_t SIZE_THROUGH_ENABLED =
      (offsetof(Color::ModulationState, enabled) + sizeof(bool) + ALIGNMENT -
       1) /
      ALIGNMENT * ALIGNMENT;
  static_assert(sizeof(Color::ModulationState) == SIZE_THROUGH_ENABLED);
  std::array<int8_t, Pattern::SphereLutView::SIZE> ac, bc;
  std::array<float, Pattern::SphereLutView::SIZE> linear;
  for (size_t i = 0; i < ac.size(); ++i) {
    ac[i] = int(i * 37 % 255) - 127;
    bc[i] = int((i * 83 + 41) % 255) - 127;
    linear[i] = bc[i];
  }
  Modulation::Prepared a{}, b{};
  a.cube = {ac.data(), true};
  b.cube = {bc.data(), true};
  a.evaluate = b.evaluate = Modulation::sample_cube;
  TemporalPalettes<> palettes;
  palettes.bake(CoordinatePalette{0.2f},
                ConstantPalette{Color4(Pixel(32000, 44000, 13000), 0.8f)});
  [[maybe_unused]] auto &from = palettes.a;
  [[maybe_unused]] auto &to = palettes.b;
  PaletteFrame frame;
  palettes.bind(frame);
  frame.previous_modulation = &a;
  frame.modulation = &b;
  frame.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 1.5958333f};
  for (bool floats : {false, true}) {
    b.cube.linear_data = floats ? linear.data() : nullptr;
    const auto PREPARED = Layer::prepare(frame);
    HS_EXPECT_TRUE(PREPARED.paired_cube_ready);
    auto dynamic = PREPARED;
    dynamic.paired_cube_ready = false;
    for (float fraction : {0.00001f, 0.25f, 0.75f, 0.99999f}) {
      frame.temporal.fraction = fraction;
      const auto PALETTE = Policy::prepare(frame);
      for (int x = -3; x <= 3; ++x) {
        const FieldSample FIELD{0.13f, 0.47f, math::Vector(x, 0.7f, -1.0f),
                                0.3f};
        const auto EXPECTED = Policy::apply_fused_modulation(
            FIELD, frame, PALETTE, dynamic, fraction);
        const auto ACTUAL = Policy::apply_fused_modulation(
            FIELD, frame, PALETTE, PREPARED, fraction);
        HS_EXPECT_EQ(ACTUAL.pending_hue, EXPECTED.pending_hue);
        HS_EXPECT_EQ(ACTUAL.color.color, EXPECTED.color.color);
        HS_EXPECT_EQ(ACTUAL.color.alpha, EXPECTED.color.alpha);
        const auto FINISHED = Policy::shade_fused_modulation(
            FIELD, frame, PALETTE, PREPARED, fraction);
        const auto REFERENCE = finish_color_sample(EXPECTED);
        HS_EXPECT_EQ(FINISHED.color, REFERENCE.color);
        HS_EXPECT_EQ(FINISHED.alpha, REFERENCE.alpha);
      }
      const auto PIPELINE_FRAME = TemporalPalettePipeline::prepare(frame);
      HS_EXPECT_TRUE(std::get<2>(PIPELINE_FRAME.prepared).paired_cube_ready);
    }
  }
  static int calls;
  const auto CUSTOM = [](const Modulation::Prepared &prepared,
                         const Pattern::SphereDirection &) {
    ++calls;
    return prepared.state.phase.phase_turns;
  };
  a.state.phase.phase_turns = -0.3f;
  b.state.phase.phase_turns = 0.7f;
  for (int custom = 1; custom <= 3; ++custom) {
    a.evaluate = custom & 1 ? CUSTOM : Modulation::sample_cube;
    b.evaluate = custom & 2 ? CUSTOM : Modulation::sample_cube;
    const auto PREPARED = Layer::prepare(frame);
    HS_EXPECT_FALSE(PREPARED.paired_cube_ready);
    calls = 0;
    const auto PALETTE = Policy::prepare(frame);
    (void)Policy::shade_fused_modulation(
        FieldSample{0.5f, 0.5f, math::Y_AXIS, 0.0f}, frame, PALETTE, PREPARED,
        frame.temporal.fraction);
    HS_EXPECT_EQ(calls, custom == 3 ? 2 : 1);
  }
  a.evaluate = b.evaluate = Modulation::sample_cube;
  for (int inactive = 1; inactive <= 3; ++inactive) {
    a.cube.active = !(inactive & 1);
    b.cube.active = !(inactive & 2);
    HS_EXPECT_FALSE(Layer::prepare(frame).paired_cube_ready);
  }
  a.cube.active = b.cube.active = true;
  frame.modulation = nullptr;
  HS_EXPECT_FALSE(Layer::prepare(frame).paired_cube_ready);
  frame.modulation = &b;
  frame.previous_modulation = nullptr;
  HS_EXPECT_FALSE(Layer::prepare(frame).paired_cube_ready);
  frame.modulation = nullptr;
  HS_EXPECT_FALSE(Layer::prepare(frame).paired_cube_ready);
  const Color::ModulationState MANUAL{&b, &a, Pattern::SourceKind::NOISE,
                                      frame.mechanic, true};
  HS_EXPECT_FALSE(MANUAL.paired_cube_ready);
}

/** @brief Paired cube samples preserve individual samplers and custom dispatch. */
inline void test_pullback_modulation_paired_samples() {
  using namespace Pullback;
  std::array<int8_t, Pattern::SphereLutView::SIZE> previous_cube;
  std::array<int8_t, Pattern::SphereLutView::SIZE> current_cube;
  std::array<float, Pattern::SphereLutView::SIZE> linear_cube;
  for (size_t index = 0; index < previous_cube.size(); ++index) {
    previous_cube[index] = static_cast<int8_t>((index * 37) % 255 - 127);
    current_cube[index] = static_cast<int8_t>((index * 83 + 19) % 255 - 127);
    linear_cube[index] = current_cube[index];
  }
  Modulation::Prepared previous{};
  Modulation::Prepared current{};
  previous.cube = {previous_cube.data(), true};
  current.cube = {current_cube.data(), true};
  previous.evaluate = current.evaluate = Modulation::sample_cube;
  for (bool floats : {false, true}) {
    current.cube.linear_data = floats ? linear_cube.data() : nullptr;
    for (int x = -3; x <= 3; ++x)
      for (int y = -3; y <= 3; ++y)
        for (int z = -3; z <= 3; ++z) {
          if (x == 0 && y == 0 && z == 0)
            continue;
          const auto DIRECTION =
              Pattern::SphereDirection::admit(math::Vector(x, y, z));
          const auto PAIR =
              Modulation::sample_pair(previous, current, DIRECTION);
          HS_EXPECT_EQ(PAIR.first, Modulation::sample(previous, DIRECTION));
          HS_EXPECT_EQ(PAIR.second, Modulation::sample(current, DIRECTION));
        }
  }
  static int calls;
  calls = 0;
  const auto CUSTOM = [](const Modulation::Prepared &prepared,
                         const Pattern::SphereDirection &) {
    ++calls;
    return prepared.state.phase.phase_turns;
  };
  previous.evaluate = current.evaluate = CUSTOM;
  previous.state.phase.phase_turns = -2.0f;
  current.state.phase.phase_turns = 3.0f;
  const auto PAIR = Modulation::sample_pair(
      previous, current, Pattern::SphereDirection::from_unit(math::Y_AXIS));
  HS_EXPECT_EQ(calls, 2);
  HS_EXPECT_EQ(PAIR.first, -1.0f);
  HS_EXPECT_EQ(PAIR.second, 1.0f);
}

/** @brief Paired modulation keeps canonical endpoint colors and mechanic order. */
inline void test_pullback_modulation_paired_colors() {
  using namespace Pullback;
  std::array<int8_t, Pattern::SphereLutView::SIZE> previous_cube;
  std::array<int8_t, Pattern::SphereLutView::SIZE> current_cube;
  for (size_t index = 0; index < previous_cube.size(); ++index) {
    previous_cube[index] = static_cast<int8_t>((index * 17) % 255 - 127);
    current_cube[index] = static_cast<int8_t>((index * 41 + 37) % 255 - 127);
  }
  Modulation::Prepared previous{};
  Modulation::Prepared current{};
  previous.cube = {previous_cube.data(), true};
  current.cube = {current_cube.data(), true};
  previous.evaluate = current.evaluate = Modulation::sample_cube;
  Color::ModulationState state{
      &current, &previous, Pattern::SourceKind::NOISE, {}, true};
  for (auto mechanic : {ColorModulation::Mechanic::HUE_ROTATION,
                        ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT,
                        ColorModulation::Mechanic::OPACITY})
    for (float amount : {0.0f, 0.3f, 1.5f})
      for (bool same_direction : {false, true})
        for (int x = -2; x <= 2; ++x) {
          if (amount > 1.0f &&
              mechanic != ColorModulation::Mechanic::HUE_ROTATION)
            continue;
          state.mechanic = {mechanic, amount};
          ColorSample a{Color4(Pixel(41000, 17000, 33000), 0.25f),
                        {0.3f, 0.4f, math::Vector(x, 1.0f, 0.3f), 0.7f}};
          ColorSample b{
              Color4(Pixel(27000, 53000, 9000), 0.75f),
              {0.6f, 0.8f,
               same_direction ? a.field.sphere : math::Vector(0.3f, x, 1.0f),
               0.2f}};
          a.pending_hue = 0.93f;
          b.pending_hue = 0.07f;
          a.pending_lift = 0.2f;
          b.pending_lift = 0.4f;
          const auto EXPECTED_A = Color::apply_modulation(a, state, true);
          const auto EXPECTED_B = Color::apply_modulation(b, state, false);
          Color::apply_modulation_pair_in_place(a, b, state);
          for (auto pair :
               {std::pair{a, EXPECTED_A}, std::pair{b, EXPECTED_B}}) {
            HS_EXPECT_EQ(pair.first.color.color, pair.second.color.color);
            HS_EXPECT_EQ(pair.first.color.alpha, pair.second.color.alpha);
            HS_EXPECT_EQ(pair.first.pending_hue, pair.second.pending_hue);
            HS_EXPECT_EQ(pair.first.pending_lift, pair.second.pending_lift);
            HS_EXPECT_EQ(pair.first.field.sphere, pair.second.field.sphere);
            const auto FINISHED = Color::finish(pair.first);
            const auto EXPECTED = Color::finish(pair.second);
            HS_EXPECT_EQ(FINISHED.color, EXPECTED.color);
            HS_EXPECT_EQ(FINISHED.alpha, EXPECTED.alpha);
          }
        }
}

/** @brief Modulation identities retain common and selected inputs only. */
inline void test_pullback_modulation_field_identity() {
  using namespace Pullback::Modulation;
  using Pattern::SourceKind;
  Config config;
  config.version = 7;
  config.evaluation = Pattern::Evaluation::DIRECT;
  config.noise.scale = 2.0f;
  config.sine.frequency = 3.0f;
  config.interference.frequency_b = 4.0f;
  config.grid.offset_b = 0.5f;
  config.seed = 123;
  const State STATE{SourceKind::NOISE, 1, {0.25f}, true};
  for (SourceKind source :
       {SourceKind::NOISE, SourceKind::DIRECTIONAL_SINE,
        SourceKind::SPHERE_INTERFERENCE, SourceKind::SPHERE_SMOOTH_GRID,
        SourceKind::DISPLACEMENT}) {
    config.source = source;
    const FieldIdentity IDENTITY = field_identity(config, STATE);
    HS_EXPECT_EQ(IDENTITY.source, source);
    HS_EXPECT_EQ(IDENTITY.version, 7);
    HS_EXPECT_EQ(IDENTITY.evaluation, Pattern::Evaluation::DIRECT);
    HS_EXPECT_EQ(IDENTITY.phase, Pattern::PhaseState{0.25f});
    HS_EXPECT_EQ(IDENTITY.bake_layout_version,
                 Pattern::SphereLutView::LAYOUT_VERSION);
    HS_EXPECT_EQ(IDENTITY.seed, source == SourceKind::NOISE ? 123 : 0);
    HS_EXPECT_EQ(IDENTITY.noise, source == SourceKind::NOISE
                                     ? config.noise
                                     : Pattern::NoiseParams{});
    HS_EXPECT_EQ(IDENTITY.sine, source == SourceKind::DIRECTIONAL_SINE
                                    ? config.sine
                                    : Pattern::SineParams{});
    HS_EXPECT_EQ(IDENTITY.interference,
                 source == SourceKind::SPHERE_INTERFERENCE
                     ? config.interference
                     : Pattern::InterferenceParams{});
    HS_EXPECT_EQ(IDENTITY.grid, source == SourceKind::SPHERE_SMOOTH_GRID
                                    ? config.grid
                                    : Pattern::GridParams{});
  }
  config.source = static_cast<SourceKind>(255);
  const FieldIdentity UNKNOWN = field_identity(config, STATE);
  HS_EXPECT_EQ(UNKNOWN.source, static_cast<SourceKind>(255));
  HS_EXPECT_EQ(UNKNOWN.version, 7);
  HS_EXPECT_EQ(UNKNOWN.evaluation, Pattern::Evaluation::DIRECT);
  HS_EXPECT_EQ(UNKNOWN.phase, STATE.phase);
  HS_EXPECT_EQ(UNKNOWN.noise, Pattern::NoiseParams{});
  HS_EXPECT_EQ(UNKNOWN.sine, Pattern::SineParams{});
  HS_EXPECT_EQ(UNKNOWN.interference, Pattern::InterferenceParams{});
  HS_EXPECT_EQ(UNKNOWN.grid, Pattern::GridParams{});
  HS_EXPECT_EQ(UNKNOWN.seed, 0);
  HS_EXPECT_EQ(UNKNOWN.bake_layout_version,
               Pattern::SphereLutView::LAYOUT_VERSION);
}

/** @brief Modulation equality distinguishes field inputs from dormant controls. */
inline void test_pullback_modulation_field_equality() {
  using namespace Pullback::Modulation;
  using Pattern::SourceKind;
  for (SourceKind source :
       {SourceKind::NOISE, SourceKind::DIRECTIONAL_SINE,
        SourceKind::SPHERE_INTERFERENCE, SourceKind::SPHERE_SMOOTH_GRID,
        SourceKind::DISPLACEMENT}) {
    Prepared baseline{};
    baseline.config.source = source;
    baseline.config.evaluation = Pattern::Evaluation::DIRECT;
    const auto expect = [&](auto change, bool equal) {
      Prepared other = baseline;
      change(other);
      HS_EXPECT_EQ(same_field(baseline, other), equal);
      HS_EXPECT_EQ(same_field(other, baseline), equal);
      HS_EXPECT_EQ(field_identity(baseline.config, baseline.state) ==
                       field_identity(other.config, other.state),
                   equal);
    };
    HS_EXPECT_TRUE(same_field(baseline, baseline));
    expect(
        [&](auto &p) {
          p.config.source = source == SourceKind::NOISE
                                ? SourceKind::DISPLACEMENT
                                : SourceKind::NOISE;
        },
        false);
    expect([](auto &p) { ++p.config.version; }, false);
    expect(
        [](auto &p) { p.config.evaluation = Pattern::Evaluation::BAKED_CUBE; },
        false);
    expect([](auto &p) { p.state.phase.phase_turns = 0.25f; }, false);
    expect([](auto &p) { ++p.config.seed; }, source != SourceKind::NOISE);
    expect([](auto &p) { p.config.noise.scale = 2.0f; },
           source != SourceKind::NOISE);
    expect([](auto &p) { p.config.noise.speed = 0.5f; },
           source != SourceKind::NOISE);
    expect([](auto &p) { p.config.sine.frequency = 2.0f; },
           source != SourceKind::DIRECTIONAL_SINE);
    expect([](auto &p) { p.config.sine.axis.x = 0.25f; },
           source != SourceKind::DIRECTIONAL_SINE);
    expect([](auto &p) { p.config.interference.frequency = 2.0f; },
           source != SourceKind::SPHERE_INTERFERENCE);
    expect([](auto &p) { p.config.interference.axis_b.z = 0.25f; },
           source != SourceKind::SPHERE_INTERFERENCE);
    expect([](auto &p) { p.config.grid.frequency = 2.0f; },
           source != SourceKind::SPHERE_SMOOTH_GRID);
    expect([](auto &p) { p.config.grid.offset_b = 0.25f; },
           source != SourceKind::SPHERE_SMOOTH_GRID);
    expect([](auto &p) { p.config.enabled = false; }, true);
    expect(
        [](auto &p) {
          p.config.mechanic.kind = ColorModulation::Mechanic::OPACITY;
        },
        true);
    expect([](auto &p) { p.config.mechanic.amount = 0.5f; }, true);
    expect([](auto &p) { p.state.source = static_cast<SourceKind>(255); },
           true);
    expect([](auto &p) { ++p.state.version; }, true);
    expect([](auto &p) { p.state.initialized = true; }, true);
  }
  Prepared unknown{};
  unknown.config.source = static_cast<SourceKind>(255);
  HS_EXPECT_FALSE(same_field(unknown, unknown));
}
