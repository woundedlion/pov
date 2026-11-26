/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Shared fixtures for color tests. */

#include "core/render/pullback.h"
#include "tests/pixel_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test::pullback_tests {

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

/** @brief Owns two cube samplers and their backing storage. */
struct CubeModulationPair {
  std::array<int8_t, Pattern::SphereLutView::SIZE> previous_cube, current_cube;
  Pullback::Modulation::Prepared previous{}, current{};

  CubeModulationPair() {
    for (size_t i = 0; i < previous_cube.size(); ++i) {
      previous_cube[i] = int(i * 37 % 255) - 127;
      current_cube[i] = int((i * 83 + 41) % 255) - 127;
    }
    previous.cube = {previous_cube.data(), true};
    current.cube = {current_cube.data(), true};
    previous.evaluate = current.evaluate = Pullback::Modulation::sample_cube;
  }
  CubeModulationPair(const CubeModulationPair &) = delete;
  CubeModulationPair &operator=(const CubeModulationPair &) = delete;
};

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
    const float index = value * (BakedPalette::LUT_SIZE - 1);
    return Color4(Pixel(1000, 2000, 3000),
                  index > 0.5f && index < 1.5f ? 0.5f : 1.0f);
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

template <typename From, typename To>
concept CanBakePaletteLab = requires(PaletteLabLut &cache) {
  cache.bake(std::declval<From>(), std::declval<To>());
};
static_assert(CanBakePaletteLab<const BakedPalette &, const BakedPalette &>);
static_assert(!CanBakePaletteLab<BakedPalette, const BakedPalette &>);
static_assert(!CanBakePaletteLab<const BakedPalette &, BakedPalette>);
static_assert(!CanBakePaletteLab<BakedPalette, BakedPalette>);
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

struct CachedSinePhaseProvider : PhaseLayerProvider {
  static constexpr bool SINE_PHASE = true;
};

/** @brief Displacement source selection for modulation qualification. */
struct DisplacementLayerProvider : PaletteLayerProvider {
  static Pullback::Color::ModulationState prepare(const PaletteFrame &frame) {
    auto prepared = PaletteLayerProvider::prepare(frame);
    prepared.source = Pattern::SourceKind::DISPLACEMENT;
    return prepared;
  }
};

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

struct CoordinatePalette {
  float alpha;
  Color4 get(float value) const {
    return Color4(Pixel(value * 60000, (1.0f - value) * 50000, 7000), alpha);
  }
};

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
  hue_cos_sin(turns, cosine, sine);
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

} // namespace hs_test::pullback_tests
