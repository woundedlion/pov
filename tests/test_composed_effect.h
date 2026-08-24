/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Base-contract sweep over every Pullback::ComposedEffect specialization: the
 * slider set init() registers, the schema-versioned parameter snapshot, the
 * preset choreography configure_presets() wires, and the family-by-family
 * parameter interpolation a preset crossfade runs on. The base owns all four,
 * so each check is expressed against the effect's parameter families and driven
 * over the whole group rather than over one effect.
 *
 * The specializations come from HS_SHADER_PRODUCT_GROUP, so a promoted effect
 * joins the sweep with its roster entry instead of a hand-copied list.
 *
 * The module also pins how much of the chain interpreter's operator catalog
 * the derivation layer can reach: the workbench catalog is wider, and the
 * difference is recorded against the live OPERATOR_TABLE rather than restated.
 *
 * Shading is deliberately out of scope: the ShaderWorkbench equivalence
 * oracles in
 * tests/test_lattice_melt.h and tests/test_kaleidoscope_smooth.h own that comparison.
 */
#pragma once

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>
#include <type_traits>

#include "core/engine/effects.h"
#include "core/render/pullback/composed_effect.h"
#include "core/render/pullback/operator_table.h"
#include "tests/test_effects.h" // reset_effect_globals, SMALL_W/SMALL_H
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace composed_effect_tests {

using effects_tests::preset_params_or_initial;
using effects_tests::reset_effect_globals;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;

/**
 * @brief The base-template arguments one specialization instantiates
 *        Pullback::ComposedEffect with.
 * @details An effect states its hue mode and brightness envelope only in its
 * base-class list, and those two decide which color sliders the base registers,
 * so they are deduced from the base rather than restated per effect.
 */
template <typename ParamsT, typename SpecT, PaletteHarmony HarmonyV,
          Pullback::HueMode HueV,
          Pullback::Color::BrightnessEnvelope BrightnessV, bool OuterNoiseV,
          bool SourceNoiseV>
struct ComposedTraits {
  using Params = ParamsT;
  using Spec = SpecT;
  static constexpr PaletteHarmony HARMONY = HarmonyV;
  static constexpr Pullback::HueMode HUE = HueV;
  static constexpr Pullback::Color::BrightnessEnvelope BRIGHTNESS = BrightnessV;
  static constexpr bool OUTER_NOISE = OuterNoiseV;
  static constexpr bool SOURCE_NOISE = SourceNoiseV;
};

template <int W, int H, typename Derived, typename ParamsT, typename SpecT,
          PaletteHarmony HarmonyV, Pullback::HueMode HueV,
          Pullback::Color::BrightnessEnvelope BrightnessV, bool OuterNoiseV,
          bool SourceNoiseV>
ComposedTraits<ParamsT, SpecT, HarmonyV, HueV, BrightnessV, OuterNoiseV,
               SourceNoiseV>
composed_traits(const Pullback::ComposedEffect<W, H, Derived, ParamsT, SpecT,
                                               HarmonyV, HueV, BrightnessV,
                                               OuterNoiseV, SourceNoiseV> &);

/** @brief ComposedTraits of the base @p FX derives from. */
template <typename FX>
using TraitsOf = decltype(composed_traits(std::declval<const FX &>()));

/** @brief Color sliders the base registers for every specialization. */
constexpr const char *UNGATED_COLOR_SLIDERS[] = {
    "Palette Chroma",    "Palette Mapping",         "Mapping Frequency",
    "Mapping Phase",     "Phase Oscillation Depth", "Phase Oscillation Speed",
    "Value Opacity Low", "Value Opacity High",      "Hue Shift Amount"};

/** @brief The eight lens sliders a Mobius parameter family adds. */
constexpr const char *MOBIUS_SLIDERS[] = {
    "Mobius A Re", "Mobius A Im", "Mobius B Re", "Mobius B Im",
    "Mobius C Re", "Mobius C Im", "Mobius D Re", "Mobius D Im"};

/** @brief A hand-registered color slider and the field descriptor it writes. */
struct ColorSliderBinding {
  const char *slider;
  const char *field_id;
};

/**
 * @brief Every color slider the base registers, paired with its descriptor.
 * @details ColorParams tables its fields with a null name, so the base names and
 * bounds them by hand instead of through register_fields(). The pairing is what
 * lets a slider's authored range be compared against the range the snapshot
 * validator enforces.
 */
constexpr ColorSliderBinding COLOR_SLIDER_BINDINGS[] = {
    {"Hue Shift Amount", "hue-shift-amount"},
    {"Hue Noise Scale", "hue-noise-scale"},
    {"Hue Noise Speed", "hue-noise-speed"},
    {"Palette Chroma", "palette-chroma"},
    {"Mapping Frequency", "mapping-frequency"},
    {"Mapping Phase", "mapping-phase"},
    {"Phase Oscillation Depth", "phase-oscillation-depth"},
    {"Phase Oscillation Speed", "phase-oscillation-speed"},
    {"Brightness Depth", "brightness-depth"},
    {"Value Opacity Low", "value-opacity-low"},
    {"Value Opacity High", "value-opacity-high"}};

inline uint32_t bits(float value) { return std::bit_cast<uint32_t>(value); }

/** @brief Whether @p FX opens a gated field's slider. */
template <typename FX> constexpr bool gate_open(Pullback::FieldGate gate) {
  switch (gate) {
  case Pullback::FieldGate::ALWAYS:
    return true;
  case Pullback::FieldGate::ANIMATED_PROJECTION:
    return FX::ANIMATED_PROJECTION;
  case Pullback::FieldGate::CENTRAL_MERIDIAN:
    return Pullback::uses_central_meridian(FX::Spec::PROJECTION);
  }
  return false;
}

/** @brief The descriptor carrying @p id in a family's table, or null. */
template <Pullback::HasFields Family>
inline const Pullback::Field<Family> *find_field(const char *id) {
  for (const auto &field : Family::FIELDS)
    if (std::string_view(field.id) == id)
      return &field;
  return nullptr;
}

/** @brief Sliders a family contributes through register_fields(). */
template <typename FX, typename Family> constexpr size_t named_field_count() {
  size_t count = 0;
  if constexpr (Pullback::HasFields<Family>)
    for (const auto &field : Family::FIELDS)
      if (field.name != nullptr && gate_open<FX>(field.gate))
        ++count;
  return count;
}

/** @brief Sliders a warp slot contributes: the slot speed plus its own. */
template <typename FX, typename Family> constexpr size_t warp_slot_count() {
  if constexpr (std::is_same_v<Family, Pullback::NoWarpParams>)
    return 0;
  else
    return 1 + named_field_count<FX, Family>();
}

/**
 * @brief Checks one family's tabled fields against the registered sliders.
 * @param params The effect's registered parameter list.
 * @param registered Whether the base runs register_fields() on this family.
 * @details A named field must own exactly one slider, present on the gate its
 * descriptor names and carrying the descriptor's bounds.
 */
template <typename FX, Pullback::HasFields Family>
inline void verify_family_sliders(const ParamList &params, bool registered) {
  for (const auto &field : Family::FIELDS) {
    if (field.name == nullptr)
      continue;
    HS_CONTEXT(field.id);
    const bool expected = registered && gate_open<FX>(field.gate);
    const ParamDef *def = params.find(field.name);
    HS_EXPECT_EQ(def != nullptr, expected);
    if (def == nullptr)
      continue;
    HS_EXPECT_EQ(bits(def->min), bits(field.min));
    HS_EXPECT_EQ(bits(def->max), bits(field.max));
    HS_EXPECT_TRUE(def->animated);
  }
}

/** @brief Checks a warp slot's speed slider against the slot's descriptor. */
template <typename FX, typename Family>
inline void verify_warp_slot(const ParamList &params, const char *slot_name) {
  HS_CONTEXT(slot_name);
  constexpr bool expected = !std::is_same_v<Family, Pullback::NoWarpParams>;
  const ParamDef *def = params.find(slot_name);
  HS_EXPECT_EQ(def != nullptr, expected);
  if (def != nullptr) {
    HS_EXPECT_EQ(bits(def->min), bits(Family::FIELDS[0].min));
    HS_EXPECT_EQ(bits(def->max), bits(Family::FIELDS[0].max));
  }
  verify_family_sliders<FX, Family>(params, expected);
}

/** @brief Bitwise-compares every tabled field of one parameter family. */
template <Pullback::HasFields Family>
inline void verify_family_equal(const Family &actual, const Family &expected) {
  for (const auto &field : Family::FIELDS) {
    HS_CONTEXT(field.id);
    HS_EXPECT_EQ(bits(actual.*(field.member)), bits(expected.*(field.member)));
  }
}

/** @brief Bitwise-compares the eight Mobius coefficients. */
inline void verify_mobius_equal(const MobiusParams &actual,
                                const MobiusParams &expected) {
  constexpr Complex MobiusParams::*COEFFICIENTS[] = {
      &MobiusParams::a, &MobiusParams::b, &MobiusParams::c, &MobiusParams::d};
  for (Complex MobiusParams::*coefficient : COEFFICIENTS) {
    HS_EXPECT_EQ(bits((actual.*coefficient).re),
                 bits((expected.*coefficient).re));
    HS_EXPECT_EQ(bits((actual.*coefficient).im),
                 bits((expected.*coefficient).im));
  }
}

/** @brief Bitwise-compares a whole parameter set, family by family. */
template <typename Params>
inline void verify_params_equal(const Params &actual, const Params &expected) {
  verify_family_equal(actual.source, expected.source);
  verify_family_equal(actual.projection, expected.projection);
  verify_family_equal(actual.outer_warp, expected.outer_warp);
  verify_family_equal(actual.inner_warp, expected.inner_warp);
  verify_family_equal(actual.surface, expected.surface);
  verify_family_equal(actual.value, expected.value);
  verify_family_equal(actual.color, expected.color);
  HS_EXPECT_EQ(static_cast<int>(actual.color.palette_mapping),
               static_cast<int>(expected.color.palette_mapping));
  if constexpr (requires { actual.lens.mobius; })
    verify_mobius_equal(actual.lens.mobius, expected.lens.mobius);
}

/** @brief Moves every tabled field of a family to the middle of its range. */
template <Pullback::HasFields Family>
inline void fill_midpoints(Family &family) {
  for (const auto &field : Family::FIELDS)
    family.*(field.member) = 0.5f * (field.min + field.max);
}

/**
 * @brief Checks that a restore refuses each field of one family out of range.
 * @param effect Effect under test.
 * @param captured An admissible snapshot each poisoned copy starts from.
 * @param slot The family's member in the parameter set.
 * @param label Family name, for the failure context.
 */
template <typename FX, typename ParamsT, Pullback::HasFields Family>
inline void
verify_family_rejection(FX &effect,
                        const typename FX::ParameterSnapshot &captured,
                        Family ParamsT::*slot, const char *label) {
  HS_CONTEXT(label);
  for (const auto &field : Family::FIELDS) {
    HS_CONTEXT(field.id);
    const float span = field.max - field.min;
    const float poisons[] = {std::numeric_limits<float>::quiet_NaN(),
                             field.min - span - 1.0f, field.max + span + 1.0f};
    for (float poison : poisons) {
      typename FX::ParameterSnapshot snapshot = captured;
      (snapshot.params.*slot).*(field.member) = poison;
      HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
    }
  }
}

/**
 * @brief Pins the slider set one specialization's init() registers.
 * @tparam E Composed effect class template.
 * @param name Effect name, for the failure context.
 * @details The registered count is derived from the parameter families rather
 * than listed, so an extra slider the base grows, a family the base stops
 * registering, or a gated field that stops consulting its gate all show up.
 * Duplicate names, capacity overflow and an out-of-range default are hard
 * checks inside register_param(), so they are not restated here.
 */
template <template <int, int> class E>
inline void check_slider_registration(const char *name) {
  using FX = E<SMALL_W, SMALL_H>;
  using Params = typename FX::Params;
  using Traits = TraitsOf<FX>;
  HS_CONTEXT(name);

  reset_effect_globals();
  FX effect;
  effect.init();
  const ParamList &params = effect.getParameters();

  // Two sliders sharing a target leave one family member unreachable while the
  // list still reads complete.
  for (const ParamDef *slider = params.begin(); slider != params.end();
       ++slider) {
    HS_EXPECT_TRUE(slider->name != nullptr);
    HS_EXPECT_TRUE(slider->target != nullptr);
    for (const ParamDef *earlier = params.begin(); earlier != slider; ++earlier)
      HS_EXPECT_TRUE(slider->target != earlier->target);
  }

  verify_family_sliders<FX, typename Params::source_type>(params, true);
  verify_family_sliders<FX, Pullback::ProjectionParams>(params, true);
  verify_family_sliders<FX, typename Params::surface_type>(params, true);
  verify_family_sliders<FX, typename Params::value_type>(params, true);
  verify_warp_slot<FX, typename Params::outer_warp_type>(params,
                                                         "Planar Warp 1 Speed");
  verify_warp_slot<FX, typename Params::inner_warp_type>(params,
                                                         "Planar Warp 2 Speed");

  for (const char *slider : UNGATED_COLOR_SLIDERS)
    HS_EXPECT_TRUE(params.find(slider) != nullptr);
  constexpr bool mobius =
      std::is_same_v<typename Params::lens_type, Pullback::MobiusLensParams>;
  for (const char *slider : MOBIUS_SLIDERS)
    HS_EXPECT_EQ(params.find(slider) != nullptr, mobius);
  constexpr bool brightness =
      Traits::BRIGHTNESS != Pullback::Color::BrightnessEnvelope::NONE;
  constexpr bool hue_noise = Traits::HUE == Pullback::HueMode::NOISE;
  HS_EXPECT_EQ(params.find("Brightness Depth") != nullptr, brightness);
  HS_EXPECT_EQ(params.find("Hue Noise Scale") != nullptr, hue_noise);
  HS_EXPECT_EQ(params.find("Hue Noise Speed") != nullptr, hue_noise);

  // Hand-registered sliders use the same ranges as snapshot validation.
  for (const ColorSliderBinding &binding : COLOR_SLIDER_BINDINGS) {
    HS_CONTEXT(binding.slider);
    const ParamDef *def = params.find(binding.slider);
    if (def == nullptr)
      continue;
    const Pullback::Field<Pullback::ColorParams> *field =
        find_field<Pullback::ColorParams>(binding.field_id);
    HS_EXPECT_TRUE(field != nullptr);
    if (field == nullptr)
      continue;
    HS_EXPECT_EQ(bits(def->min), bits(field->min));
    HS_EXPECT_EQ(bits(def->max), bits(field->max));
  }

  const ParamDef *mapping = params.find("Palette Mapping");
  HS_EXPECT_TRUE(mapping != nullptr);
  if (mapping != nullptr) {
    HS_EXPECT_EQ(mapping->option_count, 4);
    HS_EXPECT_TRUE(mapping->options != nullptr);
    HS_EXPECT_TRUE(mapping->export_options != nullptr);
  }

  const size_t expected =
      named_field_count<FX, typename Params::source_type>() +
      named_field_count<FX, Pullback::ProjectionParams>() +
      named_field_count<FX, typename Params::surface_type>() +
      named_field_count<FX, typename Params::value_type>() +
      warp_slot_count<FX, typename Params::outer_warp_type>() +
      warp_slot_count<FX, typename Params::inner_warp_type>() +
      (mobius ? std::size(MOBIUS_SLIDERS) : size_t{0}) +
      std::size(UNGATED_COLOR_SLIDERS) + (brightness ? size_t{1} : size_t{0}) +
      (hue_noise ? size_t{2} : size_t{0});
  HS_EXPECT_EQ(params.size(), expected);
  HS_EXPECT_LE(params.size(), params.capacity());
}

/**
 * @brief Pins the schema-versioned snapshot contract of one specialization.
 * @tparam E Composed effect class template.
 * @param name Effect name, for the failure context.
 * @details Covers the round trip, adoption of a whole moved parameter set, and
 * the rejections: a schema version the effect did not author, a non-finite
 * field, a field outside the range its descriptor declares, and a palette
 * mapping outside the enum. Every rejection is followed by a read-back, which
 * must still be the captured set.
 */
template <template <int, int> class E>
inline void check_snapshot_contract(const char *name) {
  using FX = E<SMALL_W, SMALL_H>;
  using Params = typename FX::Params;
  HS_CONTEXT(name);

  reset_effect_globals();
  FX effect;
  effect.init();

  const typename FX::ParameterSnapshot captured = effect.serialize_parameters();
  HS_EXPECT_EQ(captured.schema_version, FX::PARAMETER_SCHEMA_VERSION);
  HS_EXPECT_TRUE(FX::valid_params(captured.params));
  verify_params_equal(captured.params, FX::initial_params());
  HS_EXPECT_TRUE(effect.restore_parameters(captured));
  verify_params_equal(effect.serialize_parameters().params, captured.params);

  // Every family the parameter set names has to make the crossing, so the whole
  // set moves off its authored values at once.
  typename FX::ParameterSnapshot moved = captured;
  fill_midpoints(moved.params.source);
  fill_midpoints(moved.params.projection);
  fill_midpoints(moved.params.outer_warp);
  fill_midpoints(moved.params.inner_warp);
  fill_midpoints(moved.params.surface);
  fill_midpoints(moved.params.value);
  fill_midpoints(moved.params.color);
  moved.params.color.palette_mapping = Pullback::Color::PaletteMapping::BELL;
  HS_EXPECT_TRUE(effect.restore_parameters(moved));
  verify_params_equal(effect.serialize_parameters().params, moved.params);
  HS_EXPECT_TRUE(effect.restore_parameters(captured));

  typename FX::ParameterSnapshot bumped = captured;
  bumped.schema_version += 1;
  HS_EXPECT_FALSE(effect.restore_parameters(bumped));

  verify_family_rejection(effect, captured, &Params::source, "source");
  verify_family_rejection(effect, captured, &Params::projection, "projection");
  verify_family_rejection(effect, captured, &Params::outer_warp, "outer warp");
  verify_family_rejection(effect, captured, &Params::inner_warp, "inner warp");
  verify_family_rejection(effect, captured, &Params::surface, "surface");
  verify_family_rejection(effect, captured, &Params::value, "value");
  verify_family_rejection(effect, captured, &Params::color, "color");

  typename FX::ParameterSnapshot mapping = captured;
  mapping.params.color.palette_mapping =
      static_cast<Pullback::Color::PaletteMapping>(
          static_cast<uint8_t>(Pullback::Color::PaletteMapping::REVERSE) + 1);
  HS_EXPECT_FALSE(effect.restore_parameters(mapping));

  if constexpr (requires { Params{}.lens.mobius; }) {
    typename FX::ParameterSnapshot lens = captured;
    lens.params.lens.mobius.a.re = std::numeric_limits<float>::quiet_NaN();
    HS_EXPECT_FALSE(effect.restore_parameters(lens));
    lens = captured;
    lens.params.lens.mobius.a.re =
        Pullback::MobiusLensParams::COEFFICIENT_LIMIT + 1.0f;
    HS_EXPECT_FALSE(effect.restore_parameters(lens));
  }

  verify_params_equal(effect.serialize_parameters().params, captured.params);
}

/**
 * @brief Pins the preset choreography the base wires up in init().
 * @tparam E Composed effect class template.
 * @param name Effect name, for the failure context.
 * @details Every authored preset has to be reachable, admissible under the
 * effect's own validator, and adopted whole by a manual selection — the snap
 * path adopt_params() owns.
 */
template <template <int, int> class E>
inline void check_preset_choreography(const char *name) {
  using FX = E<SMALL_W, SMALL_H>;
  HS_CONTEXT(name);

  HS_EXPECT_GT(FX::TRANSITION_DURATION, uint16_t{0});
  HS_EXPECT_GT(FX::PRESET_DWELL_FRAMES, uint16_t{0});
  HS_EXPECT_GT(FX::PRESET_IDS.size(), size_t{0});
  HS_EXPECT_TRUE(FX::valid_params(FX::initial_params()));

  reset_effect_globals();
  FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), FX::PRESET_IDS.size());
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});

  for (size_t index = 0; index < FX::PRESET_IDS.size(); ++index) {
    HS_CONTEXT("preset", static_cast<int>(index));
    HS_EXPECT_FALSE(FX::PRESET_IDS[index].empty());
    for (size_t earlier = 0; earlier < index; ++earlier)
      HS_EXPECT_TRUE(FX::PRESET_IDS[index] != FX::PRESET_IDS[earlier]);
    HS_EXPECT_TRUE(FX::valid_params(preset_params_or_initial<FX>(index)));

    HS_EXPECT_TRUE(effect.selectPreset(index));
    HS_EXPECT_EQ(effect.getPresetIndex(), index);
    verify_params_equal(effect.serialize_parameters().params,
                        preset_params_or_initial<FX>(index));
  }
  HS_EXPECT_FALSE(effect.selectPreset(FX::PRESET_IDS.size()));
}

/**
 * @brief Pins the parameter interpolation a preset crossfade runs on.
 * @tparam E Composed effect class template.
 * @param name Effect name, for the failure context.
 * @details Both endpoints must come back exactly — the crossfade's last frame is
 * what commits the incoming preset — and no sample in between may leave the
 * ranges the snapshot validator enforces, or a transition would pass through
 * parameters the effect itself refuses to restore.
 */
template <template <int, int> class E>
inline void check_preset_interpolation(const char *name) {
  using FX = E<SMALL_W, SMALL_H>;
  using Params = typename FX::Params;
  HS_CONTEXT(name);

  constexpr float PROGRESS[] = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
  for (size_t from = 0; from < FX::PRESET_IDS.size(); ++from) {
    for (size_t to = 0; to < FX::PRESET_IDS.size(); ++to) {
      HS_CONTEXT("preset pair", static_cast<int>(from), static_cast<int>(to));
      const Params a = preset_params_or_initial<FX>(from);
      const Params b = preset_params_or_initial<FX>(to);
      for (float progress : PROGRESS)
        HS_EXPECT_TRUE(FX::valid_params(Pullback::interpolate(a, b, progress)));
      verify_params_equal(Pullback::interpolate(a, b, 0.0f), a);
      verify_params_equal(Pullback::interpolate(a, b, 1.0f), b);
    }
  }
}

/** @brief Sweeps the registered slider set over every specialization. */
inline void test_composed_slider_registration() {
#define HS_COMPOSED_SLIDERS(name, seconds)                                     \
  check_slider_registration<name>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_SLIDERS)
#undef HS_COMPOSED_SLIDERS
}

/** @brief Sweeps the parameter snapshot contract over every specialization. */
inline void test_composed_snapshot_contract() {
#define HS_COMPOSED_SNAPSHOT(name, seconds)                                    \
  check_snapshot_contract<name>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_SNAPSHOT)
#undef HS_COMPOSED_SNAPSHOT
}

/** @brief Sweeps the preset choreography over every specialization. */
inline void test_composed_preset_choreography() {
#define HS_COMPOSED_PRESETS(name, seconds)                                     \
  check_preset_choreography<name>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_PRESETS)
#undef HS_COMPOSED_PRESETS
}

/** @brief Sweeps the crossfade interpolation over every specialization. */
inline void test_composed_preset_interpolation() {
#define HS_COMPOSED_INTERP(name, seconds)                                      \
  check_preset_interpolation<name>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_INTERP)
#undef HS_COMPOSED_INTERP
}

/**
 * @brief Pins the two families the base registers by hand, not by table.
 * @details register_parameters() never runs register_fields() on the color or
 * lens families, so a name added to either table would author a slider nothing
 * registers.
 */
inline void test_composed_hand_registered_families() {
  for (const auto &field : Pullback::ColorParams::FIELDS) {
    HS_CONTEXT(field.id);
    HS_EXPECT_TRUE(field.name == nullptr);
  }
  for (const auto &field : Pullback::NoLensParams::FIELDS) {
    HS_CONTEXT(field.id);
    HS_EXPECT_TRUE(field.name == nullptr);
  }
  for (const ColorSliderBinding &binding : COLOR_SLIDER_BINDINGS) {
    HS_CONTEXT(binding.slider);
    HS_EXPECT_TRUE(find_field<Pullback::ColorParams>(binding.field_id) !=
                   nullptr);
  }
  HS_EXPECT_EQ(std::size(COLOR_SLIDER_BINDINGS),
               Pullback::ColorParams::FIELDS.size());
}

namespace In = Pullback::Interp;

/** A parameter set naming the narrowed warp families. */
using ReachParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::PolarParams,
                     Pullback::VectorNoiseParams>;
using ReachBinding = Pullback::Binding<Pullback::FrameState<ReachParams>>;

template <typename Family>
concept DerivableSource = requires {
  typename Pullback::SourcePolicyFor<Family, ReachBinding>::Type;
};

// Derivation-side half of the reach table: widening one of these narrowings
// reds here before the catalog-side counts below go stale.
static_assert(DerivableSource<Pullback::GridSourceParams>);
static_assert(!DerivableSource<Pullback::NoiseSourceParams>);
static_assert(
    std::is_same_v<typename Pullback::WarpPolicyFor<
                       Pullback::PolarParams, ReachBinding, true, false>::Type,
                   Pullback::Warp::PolarChart<
                       Pullback::WarpProvider<ReachBinding, true, false>,
                       Pullback::Warp::LinearPolar, 1>>);
static_assert(
    std::is_same_v<
        typename Pullback::WarpPolicyFor<Pullback::VectorNoiseParams,
                                         ReachBinding, false, false>::Type,
        Pullback::Warp::VectorNoise<
            Pullback::WarpProvider<ReachBinding, false, false>,
            ::NoiseBasis::SIMPLEX, Pullback::Warp::FlatEnvelope>>);

/**
 * @brief How much of one catalog operator's vocabulary a ComposedEffect
 *        specialization can emit.
 * @details A null `topology_id` marks an operator no Spec or parameter family
 * selects, so none of its topology values are reachable either. Otherwise
 * `reachable` lists exactly the values the derivation layer emits for that
 * enum8 and the rest are workbench-only, so a value the catalog gains counts
 * as unreachable until the row is revisited.
 */
struct DerivationReach {
  const char *operator_id;
  const char *topology_id;
  std::array<const char *, 9> reachable; /**< Null-terminated. */
};

constexpr DerivationReach DERIVATION_REACH[] = {
    // ProjectionKind has no spelling for these four.
    {"project.peirce.v2", nullptr, {}},
    {"project.peirce-square-fast.v2", nullptr, {}},
    {"project.bonne.v2", nullptr, {}},
    {"project.airocean.v2", nullptr, {}},
    // SourcePolicyFor covers Grid, TwinWave, Spiral and PrimitiveLattice only.
    {"sample.rings.v2", nullptr, {}},
    {"sample.spherical-rings.v3", nullptr, {}},
    {"sample.fractal.v2", nullptr, {}},
    {"sample.tessellation.v2", nullptr, {}},
    {"sample.projected-noise.v2", nullptr, {}},
    {"sample.spherical-noise.v3", nullptr, {}},
    // No WarpPolicyFor specialization carries a curl-flow family.
    {"warp.curl-flow.v2", nullptr, {}},
    // TransferKind is NONE or ISO_CONTOUR.
    {"field.transfer.ridge.v2", nullptr, {}},
    {"field.transfer.smooth-bands.v2", nullptr, {}},
    // ProjectionKind names the folded gnomonic alone.
    {"project.gnomonic.v2", "hemisphere", {{"folded"}}},
    // The displacement policies pin the basis and the integrator.
    {"sphere.displace.curl.v2", "basis", {{"simplex"}}},
    {"sphere.displace.curl.v2", "integrator", {{"euler"}}},
    {"sphere.displace.direct.v2", "basis", {{"simplex"}}},
    // Spec::LensPolicy names one lens policy per symmetry.
    {"sphere.lens.kaleidoscope.v2",
     "symmetry",
     {{"azimuthal", "tetrahedral", "octahedral", "dodecahedral",
       "triangular-prism", "square-prism", "pentagonal-prism",
       "hexagonal-prism", "octagonal-prism"}}},
    // WarpPolicyFor pins the flat envelope, the simplex basis, the linear
    // polar chart and its first harmonic.
    {"warp.wave-shear.v2", "envelope", {{"flat"}}},
    {"warp.vector-noise.v2", "basis", {{"simplex"}}},
    {"warp.vector-noise.v2", "envelope", {{"flat"}}},
    {"warp.polar-chart.v2", "mode", {{"linear"}}},
    {"warp.polar-chart.v2", "harmonic", {{"h1"}}},
    // SampleStage pins Weight::Projection; CoverageKind has no None.
    {"sample.grid.v2", "weight-mode", {{"projection"}}},
    {"sample.grid.v2",
     "coverage-mode",
     {{"weight", "weight-squared", "edge-fade"}}},
    {"sample.twin-wave.v2", "weight-mode", {{"projection"}}},
    {"sample.twin-wave.v2",
     "coverage-mode",
     {{"weight", "weight-squared", "edge-fade"}}},
    {"sample.spiral.v2", "weight-mode", {{"projection"}}},
    {"sample.spiral.v2",
     "coverage-mode",
     {{"weight", "weight-squared", "edge-fade"}}},
    {"sample.lattice.v2", "weight-mode", {{"projection"}}},
    {"sample.lattice.v2",
     "coverage-mode",
     {{"weight", "weight-squared", "edge-fade"}}},
    // PaletteHarmony, ColorParams, HueMode and BrightnessEnvelope together
    // reach every colorize value.
    {"colorize.generated-palette.v2",
     "palette-mode",
     {{"triadic", "complementary", "analogous"}}},
    {"colorize.generated-palette.v2",
     "palette-mapping",
     {{"cup", "bell", "linear", "reverse"}}},
    {"colorize.generated-palette.v2",
     "hue-shift-mode",
     {{"none", "noise", "path-length"}}},
    {"colorize.generated-palette.v2", "brightness-envelope", {{"none", "cup"}}},
};

/** @brief The topology enum8 @p topology_id of @p op, or null. */
inline const In::ParamFieldInfo *find_topology(const In::OperatorDescriptor &op,
                                               const char *topology_id) {
  for (const In::ParamFieldInfo &field : op.schema_span())
    if (field.topology && std::string_view(field.id) == topology_id)
      return &field;
  return nullptr;
}

/** @brief Rows of DERIVATION_REACH naming @p topology_id of @p operator_id. */
inline size_t reach_rows(const char *operator_id, const char *topology_id) {
  size_t rows = 0;
  for (const DerivationReach &row : DERIVATION_REACH)
    if (std::string_view(row.operator_id) == operator_id &&
        row.topology_id != nullptr &&
        std::string_view(row.topology_id) == topology_id)
      ++rows;
  return rows;
}

/**
 * @brief Pins the catalog vocabulary no ComposedEffect specialization emits.
 * @details The chain interpreter's OPERATOR_TABLE is the workbench's whole
 * catalog; the derivation layer builds a strictly narrower pipeline out of an
 * effect's families and Spec. DERIVATION_REACH records that difference and is
 * resolved against the live table here, so a renamed or removed operator,
 * topology enum8 or value reds. The four totals red whenever the catalog gains
 * an operator or a value the table has not classified.
 */
inline void test_composed_derivation_reach() {
  static_assert(AshCloudSpec::FIELD_COVERAGE ==
                Pullback::FieldCoverageKind::VALUE_CUTOUT);
  size_t unreachable_operators = 0;
  size_t unreachable_values = 0;
  for (const DerivationReach &row : DERIVATION_REACH) {
    HS_CONTEXT(row.operator_id);
    const In::OperatorDescriptor *op = In::find_operator(row.operator_id);
    HS_EXPECT_TRUE(op != nullptr);
    if (op == nullptr)
      continue;
    if (row.topology_id == nullptr) {
      ++unreachable_operators;
      for (const In::ParamFieldInfo &field : op->schema_span())
        if (field.topology)
          unreachable_values += field.enum_count;
      continue;
    }
    HS_CONTEXT(row.topology_id);
    const In::ParamFieldInfo *field = find_topology(*op, row.topology_id);
    HS_EXPECT_TRUE(field != nullptr);
    if (field == nullptr)
      continue;
    size_t reachable = 0;
    for (const char *value : row.reachable) {
      if (value == nullptr)
        break;
      HS_CONTEXT(value);
      bool found = false;
      for (uint8_t index = 0; index < field->enum_count; ++index)
        found = found || std::string_view(field->enum_ids[index]) == value;
      HS_EXPECT_TRUE(found);
      reachable += found ? 1 : 0;
    }
    HS_EXPECT_LE(reachable, static_cast<size_t>(field->enum_count));
    unreachable_values += field->enum_count - reachable;
  }

  size_t catalog_values = 0;
  for (const In::OperatorDescriptor &op : In::OPERATOR_TABLE) {
    HS_CONTEXT(op.operator_id);
    bool whole_operator = false;
    for (const DerivationReach &row : DERIVATION_REACH)
      whole_operator = whole_operator ||
                       (std::string_view(row.operator_id) == op.operator_id &&
                        row.topology_id == nullptr);
    for (const In::ParamFieldInfo &field : op.schema_span()) {
      if (!field.topology)
        continue;
      HS_CONTEXT(field.id);
      catalog_values += field.enum_count;
      HS_EXPECT_EQ(reach_rows(op.operator_id, field.id),
                   whole_operator ? 0u : 1u);
    }
  }

  HS_EXPECT_EQ(In::OPERATOR_TABLE.size(), 37u);
  HS_EXPECT_EQ(unreachable_operators, 13u);
  HS_EXPECT_EQ(catalog_values, 122u);
  HS_EXPECT_EQ(unreachable_values, 76u);
}

using RippleProbeParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::NoWarpParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::NoValueParams, Pullback::PeriodicRippleParams>;
using RippleProbeSpec = Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                                       void, Pullback::TransferKind::NONE,
                                       Pullback::CoverageKind::PROJECTION>;

template <int W, int H>
class RippleProbe
    : public Pullback::ComposedEffect<
          W, H, RippleProbe<W, H>, RippleProbeParams, RippleProbeSpec,
          PaletteHarmony::TRIADIC, Pullback::HueMode::PATH_LENGTH,
          Pullback::Color::BrightnessEnvelope::NONE> {
public:
  using Params = RippleProbeParams;
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"ripple"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = false;

  static constexpr Params initial_params() {
    Params value;
    value.surface.period = 80.0f;
    value.surface.strength = 0.15f;
    value.surface.decay = 0.0f;
    value.surface.thickness = 0.7f;
    return value;
  }
};

inline void test_composed_periodic_ripple_surface() {
  using FX = RippleProbe<SMALL_W, SMALL_H>;
  static_assert(!FX::HAS_SURFACE_NOISE);

  reset_effect_globals();
  FX effect;
  effect.init();
  HS_EXPECT_TRUE(effect.getParameters().find("Ripple Strength") != nullptr);
  HS_EXPECT_TRUE(effect.getParameters().find("Ripple Period") != nullptr);
  effect.draw_frame();
  effect.advance_display();
}

/**
 * @brief Module entry point for the composed-effect base contract.
 * @return Module result code from hs_test::end_module (0 on success).
 */
inline int run_composed_effect_tests() {
  ModuleFixture fixture("composed_effect");
  test_composed_hand_registered_families();
  test_composed_slider_registration();
  test_composed_snapshot_contract();
  test_composed_preset_choreography();
  test_composed_preset_interpolation();
  test_composed_derivation_reach();
  test_composed_periodic_ripple_surface();
  return fixture.result();
}

} // namespace composed_effect_tests
} // namespace hs_test
