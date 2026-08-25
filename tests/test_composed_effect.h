/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Base-contract sweep over every Pullback::ComposedEffect specialization: the
 * slider set init() registers, the schema-versioned parameter snapshot, the
 * preset choreography begin_choreography() wires, the family-by-family
 * parameter interpolation a preset crossfade runs on, and the value-level pin
 * of each effect's authored parameters against its promoted shader document in
 * patterns/. The base owns the first four, so each check is expressed against
 * the effect's parameter families and driven over the whole group rather than
 * over one effect.
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
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

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
          Pullback::Color::BrightnessEnvelope BrightnessV,
          bool AnimatedProjectionV, bool OuterNoiseV, bool SourceNoiseV>
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
          Pullback::Color::BrightnessEnvelope BrightnessV,
          bool AnimatedProjectionV, bool OuterNoiseV, bool SourceNoiseV>
ComposedTraits<ParamsT, SpecT, HarmonyV, HueV, BrightnessV, AnimatedProjectionV,
               OuterNoiseV, SourceNoiseV>
composed_traits(const Pullback::ComposedEffect<
                W, H, Derived, ParamsT, SpecT, HarmonyV, HueV, BrightnessV,
                AnimatedProjectionV, OuterNoiseV, SourceNoiseV> &);

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
  case Pullback::FieldGate::SINGULARITY_FADE:
    return Pullback::uses_singularity_fade(FX::Spec::PROJECTION);
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

/** @brief Minimal JSON value for reading the promoted shader documents. */
struct JsonValue {
  enum class Kind : uint8_t { NUL, BOOLEAN, NUMBER, STRING, ARRAY, OBJECT };
  Kind kind = Kind::NUL;
  bool boolean = false;
  double number = 0.0;
  std::string text;
  std::vector<JsonValue> items;
  std::vector<std::pair<std::string, JsonValue>> members;

  const JsonValue *find(std::string_view key) const {
    if (kind != Kind::OBJECT)
      return nullptr;
    for (const auto &[member_key, member_value] : members)
      if (member_key == key)
        return &member_value;
    return nullptr;
  }
};

/** @brief Fail-stop recursive-descent JSON parser over a document string. */
struct JsonParser {
  std::string_view source;
  size_t position = 0;
  bool failed = false;

  void skip_space() {
    while (position < source.size() &&
           (source[position] == ' ' || source[position] == '\t' ||
            source[position] == '\n' || source[position] == '\r'))
      ++position;
  }

  bool consume(char expected) {
    skip_space();
    if (position < source.size() && source[position] == expected) {
      ++position;
      return true;
    }
    failed = true;
    return false;
  }

  std::string parse_string() {
    std::string result;
    if (!consume('"'))
      return result;
    while (position < source.size() && source[position] != '"') {
      char c = source[position++];
      if (c != '\\') {
        result.push_back(c);
        continue;
      }
      if (position >= source.size())
        break;
      const char escape = source[position++];
      switch (escape) {
      case 'b':
        result.push_back('\b');
        break;
      case 'f':
        result.push_back('\f');
        break;
      case 'n':
        result.push_back('\n');
        break;
      case 'r':
        result.push_back('\r');
        break;
      case 't':
        result.push_back('\t');
        break;
      case 'u': {
        if (position + 4 > source.size()) {
          failed = true;
          return result;
        }
        unsigned code = 0;
        for (int digit = 0; digit < 4; ++digit) {
          const char h = source[position++];
          code <<= 4;
          if (h >= '0' && h <= '9')
            code |= static_cast<unsigned>(h - '0');
          else if (h >= 'a' && h <= 'f')
            code |= static_cast<unsigned>(h - 'a' + 10);
          else if (h >= 'A' && h <= 'F')
            code |= static_cast<unsigned>(h - 'A' + 10);
          else
            failed = true;
        }
        if (code < 0x80) {
          result.push_back(static_cast<char>(code));
        } else if (code < 0x800) {
          result.push_back(static_cast<char>(0xC0 | (code >> 6)));
          result.push_back(static_cast<char>(0x80 | (code & 0x3F)));
        } else {
          result.push_back(static_cast<char>(0xE0 | (code >> 12)));
          result.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3F)));
          result.push_back(static_cast<char>(0x80 | (code & 0x3F)));
        }
        break;
      }
      default:
        result.push_back(escape);
        break;
      }
    }
    if (!consume('"'))
      failed = true;
    return result;
  }

  JsonValue parse_value() {
    JsonValue value;
    skip_space();
    if (failed || position >= source.size()) {
      failed = true;
      return value;
    }
    const char c = source[position];
    if (c == '{') {
      ++position;
      value.kind = JsonValue::Kind::OBJECT;
      skip_space();
      if (position < source.size() && source[position] == '}') {
        ++position;
        return value;
      }
      do {
        skip_space();
        std::string key = parse_string();
        if (!consume(':'))
          return value;
        value.members.emplace_back(std::move(key), parse_value());
        skip_space();
      } while (!failed && position < source.size() && source[position] == ',' &&
               ++position);
      consume('}');
      return value;
    }
    if (c == '[') {
      ++position;
      value.kind = JsonValue::Kind::ARRAY;
      skip_space();
      if (position < source.size() && source[position] == ']') {
        ++position;
        return value;
      }
      do {
        value.items.push_back(parse_value());
        skip_space();
      } while (!failed && position < source.size() && source[position] == ',' &&
               ++position);
      consume(']');
      return value;
    }
    if (c == '"') {
      value.kind = JsonValue::Kind::STRING;
      value.text = parse_string();
      return value;
    }
    if (source.substr(position, 4) == "true") {
      position += 4;
      value.kind = JsonValue::Kind::BOOLEAN;
      value.boolean = true;
      return value;
    }
    if (source.substr(position, 5) == "false") {
      position += 5;
      value.kind = JsonValue::Kind::BOOLEAN;
      return value;
    }
    if (source.substr(position, 4) == "null") {
      position += 4;
      return value;
    }
    value.kind = JsonValue::Kind::NUMBER;
    char *end = nullptr;
    const std::string number(source.substr(position, 64));
    value.number = std::strtod(number.c_str(), &end);
    if (end == number.c_str())
      failed = true;
    position += static_cast<size_t>(end - number.c_str());
    return value;
  }
};

/** @brief The whole file as a string, or empty on a read failure. */
inline std::string read_document(const std::string &path) {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  std::FILE *file = std::fopen(path.c_str(), "rb");
#pragma clang diagnostic pop
  if (file == nullptr)
    return {};
  std::string text;
  char buffer[4096];
  for (size_t n; (n = std::fread(buffer, 1, sizeof buffer, file)) > 0;)
    text.append(buffer, n);
  std::fclose(file);
  return text;
}

/** @brief The parameter family a document chain instance addresses. */
enum class SlotRole : uint8_t {
  CAMERA,
  LENS,
  SURFACE,
  PROJECT,
  WARP,
  SAMPLE,
  FIELD,
  COLORIZE,
  UNKNOWN
};

/** @brief One document chain entry, classified by its operator id. */
struct DocumentSlot {
  std::string label;
  SlotRole role = SlotRole::UNKNOWN;
  int warp_side = 0; /**< 0 = outer warp family, 1 = inner warp family. */
};

inline SlotRole classify_operator(std::string_view operator_id) {
  if (operator_id.starts_with("sphere.rotate."))
    return SlotRole::CAMERA;
  if (operator_id.starts_with("sphere.lens."))
    return SlotRole::LENS;
  if (operator_id.starts_with("sphere.displace."))
    return SlotRole::SURFACE;
  if (operator_id.starts_with("project."))
    return SlotRole::PROJECT;
  if (operator_id.starts_with("warp."))
    return SlotRole::WARP;
  if (operator_id.starts_with("sample."))
    return SlotRole::SAMPLE;
  if (operator_id.starts_with("field."))
    return SlotRole::FIELD;
  if (operator_id.starts_with("colorize."))
    return SlotRole::COLORIZE;
  return SlotRole::UNKNOWN;
}

/** @brief Writes @p value into the family field tabled under @p id, if any. */
template <typename Family>
inline bool assign_field(Family &family, std::string_view id, float value) {
  if constexpr (Pullback::HasFields<Family>) {
    for (const auto &field : Family::FIELDS)
      if (std::string_view(field.id) == id) {
        family.*(field.member) = value;
        return true;
      }
  }
  return false;
}

/** @brief Writes one Mobius coefficient addressed by its chain field id. */
template <typename Params>
inline bool assign_mobius(Params &params, std::string_view id, float value) {
  if constexpr (requires { params.lens.mobius; }) {
    auto &mobius = params.lens.mobius;
    const struct {
      std::string_view id;
      float *slot;
    } coefficients[] = {
        {"mobius-a-re", &mobius.a.re}, {"mobius-a-im", &mobius.a.im},
        {"mobius-b-re", &mobius.b.re}, {"mobius-b-im", &mobius.b.im},
        {"mobius-c-re", &mobius.c.re}, {"mobius-c-im", &mobius.c.im},
        {"mobius-d-re", &mobius.d.re}, {"mobius-d-im", &mobius.d.im},
    };
    for (const auto &coefficient : coefficients)
      if (coefficient.id == id) {
        *coefficient.slot = value;
        return true;
      }
  }
  return false;
}

namespace ChainOp = Pullback::Interp::Op;

/** @brief PALETTE_MODE_IDS entry for @p harmony, or null. */
constexpr const char *palette_mode_id(PaletteHarmony harmony) {
  switch (harmony) {
  case PaletteHarmony::TRIADIC:
    return ChainOp::PALETTE_MODE_IDS[static_cast<uint8_t>(
        ChainOp::PaletteMode::TRIADIC)];
  case PaletteHarmony::COMPLEMENTARY:
    return ChainOp::PALETTE_MODE_IDS[static_cast<uint8_t>(
        ChainOp::PaletteMode::COMPLEMENTARY)];
  case PaletteHarmony::ANALOGOUS:
    return ChainOp::PALETTE_MODE_IDS[static_cast<uint8_t>(
        ChainOp::PaletteMode::ANALOGOUS)];
  default:
    return nullptr;
  }
}

/**
 * @brief Applies one document preset entry onto @p built, or verifies it
 *        against the effect's compile-time constants.
 * @return False when the key addresses nothing this effect owns.
 * @details Numeric entries write through the owning family's field table (the
 * same descriptors the sliders and the interpolator use), so a document key the
 * engine does not table fails loudly instead of being skipped. String entries
 * are chain topology: the ones with a composed-effect equivalent are checked
 * against the effect's Spec and base-template arguments; the rest (lens
 * symmetry, noise basis/integrator, polar mode) have no per-effect constant to
 * compare and are accepted as covered by the descriptor digest.
 */
template <typename FX>
inline bool
apply_document_value(typename FX::Params &built, const DocumentSlot &slot,
                     std::string_view field_id, const JsonValue &value) {
  using Traits = TraitsOf<FX>;
  using Spec = typename Traits::Spec;
  if (value.kind == JsonValue::Kind::NUMBER) {
    const float number = static_cast<float>(value.number);
    switch (slot.role) {
    case SlotRole::CAMERA:
      if (field_id == "wander")
        return assign_field(built.projection, "camera-wander", number);
      if (field_id == "spin-speed") {
        if constexpr (requires { FX::CAMERA_SPIN_RATE; }) {
          HS_EXPECT_EQ(bits(number), bits(FX::CAMERA_SPIN_RATE));
          return true;
        }
        return false;
      }
      return false;
    case SlotRole::PROJECT:
      return assign_field(built.projection, field_id, number);
    case SlotRole::WARP:
      return slot.warp_side == 0
                 ? assign_field(built.outer_warp, field_id, number)
                 : assign_field(built.inner_warp, field_id, number);
    case SlotRole::SURFACE:
      return assign_field(built.surface, field_id, number);
    case SlotRole::SAMPLE:
      if (assign_field(built.source, field_id, number) ||
          assign_field(built.value, field_id, number))
        return true;
      // The sample operators always carry an edge-width; without edge-fade
      // coverage it is inert in the chain and has no composed-effect field.
      return field_id == "edge-width" &&
             Spec::COVERAGE != Pullback::CoverageKind::EDGE_FADE;
    case SlotRole::FIELD:
      return assign_field(built.value, field_id, number);
    case SlotRole::LENS:
      return assign_mobius(built, field_id, number);
    case SlotRole::COLORIZE:
      return assign_field(built.color, field_id, number);
    default:
      return false;
    }
  }
  if (value.kind != JsonValue::Kind::STRING)
    return false;
  const std::string_view text = value.text;
  switch (slot.role) {
  case SlotRole::COLORIZE:
    if (field_id == "palette-mapping") {
      for (uint8_t index = 0; index < std::size(ChainOp::PALETTE_MAPPING_IDS);
           ++index)
        if (text == ChainOp::PALETTE_MAPPING_IDS[index]) {
          built.color.palette_mapping =
              static_cast<Pullback::Color::PaletteMapping>(index);
          return true;
        }
      return false;
    }
    if (field_id == "palette-mode") {
      const char *expected = palette_mode_id(Traits::HARMONY);
      HS_EXPECT_TRUE(expected != nullptr && text == expected);
      return true;
    }
    if (field_id == "hue-shift-mode") {
      HS_EXPECT_TRUE(
          text ==
          ChainOp::HUE_SHIFT_MODE_IDS[static_cast<uint8_t>(Traits::HUE)]);
      return true;
    }
    if (field_id == "brightness-envelope") {
      HS_EXPECT_TRUE(text ==
                     ChainOp::BRIGHTNESS_ENVELOPE_IDS[static_cast<uint8_t>(
                         Traits::BRIGHTNESS)]);
      return true;
    }
    return false;
  case SlotRole::SAMPLE:
    if (field_id == "coverage-mode") {
      HS_EXPECT_TRUE(
          text ==
          ChainOp::COVERAGE_MODE_IDS[static_cast<uint8_t>(Spec::COVERAGE) + 1]);
      return true;
    }
    if (field_id == "weight-mode") {
      HS_EXPECT_TRUE(text == "projection");
      return true;
    }
    return false;
  case SlotRole::PROJECT:
    if (field_id == "hemisphere") {
      HS_EXPECT_TRUE(Spec::PROJECTION ==
                         Pullback::ProjectionKind::GNOMONIC_FOLDED &&
                     text == "folded");
      return true;
    }
    return false;
  case SlotRole::LENS:
    return field_id == "symmetry";
  case SlotRole::SURFACE:
    return field_id == "basis" || field_id == "integrator";
  case SlotRole::WARP:
    return field_id == "mode" || field_id == "basis";
  default:
    return false;
  }
}

/**
 * @brief Pins one effect's authored parameter values to its shader document.
 * @tparam E Composed effect class template.
 * @param name Effect name, for the failure context.
 * @details The digests pin the promoted header to the document's canonical
 * JSON, but both are computed from the JSON, so a value edited in
 * initial_params()/preset_params() alone would leave them green. This check
 * closes that gap: every preset in patterns/<effect_id>.shader.json is rebuilt
 * into a Params through the engine's own field tables and compared bit-exactly,
 * family by family, against the effect's authored preset. The comparison is
 * two-directional — a document value the header does not reproduce and a
 * header value the document does not carry both surface as a family mismatch.
 * The preset roster, dwell and segue durations are pinned too, since they are
 * authored in both places as well.
 */
template <template <int, int> class E>
inline void check_document_values(const char *name) {
  using FX = E<SMALL_W, SMALL_H>;
  using Params = typename FX::Params;
  HS_CONTEXT(name);

  std::string file_name(FX::EFFECT_ID);
  for (char &c : file_name)
    if (c == '-')
      c = '_';
  const std::string path =
      std::string(HS_PROMOTED_PATTERNS_DIR "/") + file_name + ".shader.json";
  const std::string text = read_document(path);
  HS_EXPECT(!text.empty(), "promoted shader document is readable");
  if (text.empty())
    return;

  JsonParser parser{text};
  const JsonValue document = parser.parse_value();
  HS_EXPECT(!parser.failed, "promoted shader document parses");
  if (parser.failed)
    return;

  const JsonValue *effect_id = document.find("effect_id");
  HS_EXPECT_TRUE(effect_id != nullptr && effect_id->text == FX::EFFECT_ID);

  const JsonValue *descriptor = document.find("descriptor");
  const JsonValue *chain =
      descriptor != nullptr ? descriptor->find("chain") : nullptr;
  HS_EXPECT_TRUE(chain != nullptr);
  if (chain == nullptr)
    return;
  std::vector<DocumentSlot> slots;
  int warps_seen = 0;
  for (const JsonValue &entry : chain->items) {
    const JsonValue *label = entry.find("label");
    const JsonValue *operator_id = entry.find("operator");
    HS_EXPECT_TRUE(label != nullptr && operator_id != nullptr);
    if (label == nullptr || operator_id == nullptr)
      return;
    DocumentSlot slot;
    slot.label = label->text;
    slot.role = classify_operator(operator_id->text);
    HS_EXPECT(slot.role != SlotRole::UNKNOWN, "chain operator classified");
    if (slot.role == SlotRole::WARP) {
      // The v1 expansion keeps the warp's slot position in its label even
      // when the other slot's identity op is omitted from the chain.
      if (slot.label == "warp1")
        slot.warp_side = 0;
      else if (slot.label == "warp2")
        slot.warp_side = 1;
      else
        slot.warp_side = warps_seen;
      ++warps_seen;
    }
    slots.push_back(std::move(slot));
  }

  const JsonValue *bank = document.find("preset_bank");
  const JsonValue *presets = bank != nullptr ? bank->find("presets") : nullptr;
  HS_EXPECT_TRUE(presets != nullptr);
  if (presets == nullptr)
    return;
  HS_EXPECT_EQ(presets->items.size(), FX::PRESET_IDS.size());

  const JsonValue *choreography = bank->find("choreography");
  const JsonValue *order =
      choreography != nullptr ? choreography->find("generated_order") : nullptr;
  const JsonValue *dwell =
      choreography != nullptr ? choreography->find("dwell") : nullptr;
  HS_EXPECT_TRUE(order != nullptr && dwell != nullptr);
  if (order != nullptr && order->items.size() == FX::PRESET_IDS.size())
    for (size_t index = 0; index < FX::PRESET_IDS.size(); ++index)
      HS_EXPECT_TRUE(order->items[index].text == FX::PRESET_IDS[index]);
  else
    HS_EXPECT_TRUE(order != nullptr &&
                   order->items.size() == FX::PRESET_IDS.size());
  if (dwell != nullptr)
    for (const auto &[preset_id, frames] : dwell->members) {
      HS_CONTEXT(preset_id.c_str());
      HS_EXPECT_EQ(frames.number, double{FX::PRESET_DWELL_FRAMES});
    }
  const JsonValue *edges = bank->find("edges");
  if (edges != nullptr)
    for (const JsonValue &edge : edges->items) {
      const JsonValue *duration = edge.find("duration");
      HS_EXPECT_TRUE(duration != nullptr &&
                     duration->number == double{FX::TRANSITION_DURATION});
    }

  for (size_t index = 0; index < FX::PRESET_IDS.size(); ++index) {
    HS_CONTEXT("preset", static_cast<int>(index));
    const JsonValue *preset = nullptr;
    for (const JsonValue &candidate : presets->items) {
      const JsonValue *preset_id = candidate.find("preset_id");
      if (preset_id != nullptr && preset_id->text == FX::PRESET_IDS[index])
        preset = &candidate;
    }
    HS_EXPECT_TRUE(preset != nullptr);
    if (preset == nullptr)
      continue;
    const JsonValue *values = preset->find("values");
    HS_EXPECT_TRUE(values != nullptr);
    if (values == nullptr)
      continue;

    Params built{};
    for (const auto &[key, value] : values->members) {
      HS_CONTEXT(key.c_str());
      const size_t dot = key.find('.');
      HS_EXPECT_TRUE(dot != std::string::npos);
      if (dot == std::string::npos)
        continue;
      const std::string_view label = std::string_view(key).substr(0, dot);
      const std::string_view field_id = std::string_view(key).substr(dot + 1);
      const DocumentSlot *slot = nullptr;
      for (const DocumentSlot &candidate : slots)
        if (candidate.label == label)
          slot = &candidate;
      HS_EXPECT(slot != nullptr, "value key names a chain instance");
      if (slot == nullptr)
        continue;
      HS_EXPECT(apply_document_value<FX>(built, *slot, field_id, value),
                "document value maps onto the effect");
    }
    if constexpr (requires { FX::CAMERA_SPIN_RATE; }) {
      bool camera_spin_present = false;
      for (const DocumentSlot &slot : slots)
        if (slot.role == SlotRole::CAMERA)
          camera_spin_present =
              values->find(slot.label + ".spin-speed") != nullptr;
      HS_EXPECT_TRUE(camera_spin_present);
    }
    verify_params_equal(built, preset_params_or_initial<FX>(index));
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

/** @brief Sweeps the shader-document value pin over every specialization. */
inline void test_composed_document_values() {
#define HS_COMPOSED_DOCUMENTS(name, seconds) check_document_values<name>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_DOCUMENTS)
#undef HS_COMPOSED_DOCUMENTS
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
// The shared noise field group is not itself a family: it carries the fields
// for both noise sources, so it names no single policy.
static_assert(!DerivableSource<Pullback::Source::NoiseSourceParams>);
static_assert(DerivableSource<Pullback::ProjectedNoiseSourceParams>);
static_assert(DerivableSource<Pullback::SphericalNoiseSourceParams>);
// SphericalRings stays chain-only: its policy needs a prepared axis and
// phase the composed frame does not carry.
static_assert(!DerivableSource<Pullback::Source::SphericalRingsSourceParams>);
static_assert(
    std::is_same_v<
        typename Pullback::SourcePolicyFor<Pullback::SphericalNoiseSourceParams,
                                           ReachBinding>::Type,
        Pullback::Source::SphericalNoise<Pullback::SourceProvider<ReachBinding>,
                                         ::NoiseBasis::SIMPLEX>>);
static_assert(
    std::is_same_v<
        typename Pullback::SourcePolicyFor<Pullback::ProjectedNoiseSourceParams,
                                           ReachBinding>::Type,
        Pullback::Source::ProjectedNoise<Pullback::SourceProvider<ReachBinding>,
                                         ::NoiseBasis::SIMPLEX>>);
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
    // SourcePolicyFor has no policy for these plane samplers.
    {"sample.rings.v2", nullptr, {}},
    {"sample.spherical-rings.v3", nullptr, {}},
    {"sample.fractal.v2", nullptr, {}},
    {"sample.tessellation.v2", nullptr, {}},
    // No WarpPolicyFor specialization carries these warp families.
    {"warp.vortex.v2", nullptr, {}},
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
    {"sample.projected-noise.v2", "weight-mode", {{"projection"}}},
    {"sample.projected-noise.v2",
     "coverage-mode",
     {{"weight", "weight-squared", "edge-fade"}}},
    {"sample.projected-noise.v2", "basis", {{"simplex"}}},
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

  HS_EXPECT_EQ(In::OPERATOR_TABLE.size(), 38u);
  HS_EXPECT_EQ(unreachable_operators, 12u);
  HS_EXPECT_EQ(catalog_values, 125u);
  HS_EXPECT_EQ(unreachable_values, 74u);
}

using RippleProbeParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::NoWarpParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::NoValueParams, Pullback::PeriodicRippleParams>;
using RippleProbeSpec = Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                                       void, Pullback::TransferKind::NONE,
                                       Pullback::CoverageKind::PROJECTION>;

template <int W, int H, bool AnimatedProjection = false>
class RippleProbe
    : public Pullback::ComposedEffect<
          W, H, RippleProbe<W, H, AnimatedProjection>, RippleProbeParams,
          RippleProbeSpec, PaletteHarmony::TRIADIC,
          Pullback::HueMode::PATH_LENGTH,
          Pullback::Color::BrightnessEnvelope::NONE, AnimatedProjection> {
public:
  using Params = RippleProbeParams;
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"ripple"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;

  static constexpr Params initial_params() {
    Params value;
    value.surface.period = 80.0f;
    value.surface.strength = 0.15f;
    value.surface.decay = 0.0f;
    value.surface.thickness = 0.7f;
    return value;
  }
};

struct ProjectionWalkFootprint {
  size_t arena_bytes;
  int timeline_events;
};

template <bool AnimatedProjection>
ProjectionWalkFootprint measure_projection_walk_footprint() {
  reset_effect_globals();
  RippleProbe<SMALL_W, SMALL_H, AnimatedProjection> effect;
  effect.init();
  return {persistent_arena.get_offset(), Timeline::event_count()};
}

inline void test_composed_projection_walk_storage() {
  using Static = RippleProbe<SMALL_W, SMALL_H, false>;
  using Animated = RippleProbe<SMALL_W, SMALL_H, true>;
  static_assert(sizeof(Static) < sizeof(Animated));

  const ProjectionWalkFootprint static_footprint =
      measure_projection_walk_footprint<false>();
  const ProjectionWalkFootprint animated_footprint =
      measure_projection_walk_footprint<true>();
  HS_EXPECT_LT(static_footprint.arena_bytes, animated_footprint.arena_bytes);
  HS_EXPECT_EQ(static_footprint.timeline_events + 1,
               animated_footprint.timeline_events);
}

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
  size_t lit = 0;
  for (int i = 0; i < SMALL_W * SMALL_H; ++i) {
    const Pixel &pixel = effect.display_buffer()[i];
    lit += pixel.r != 0 || pixel.g != 0 || pixel.b != 0;
  }
  HS_EXPECT_GT(lit, size_t(0));
}

template <typename SourceT>
using NoiseSourceProbeParams =
    Pullback::Params<SourceT, Pullback::NoWarpParams, Pullback::NoWarpParams>;
using NoiseSourceProbeSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC, void,
                   Pullback::TransferKind::NONE,
                   Pullback::CoverageKind::PROJECTION>;

/** @brief Probe deriving a pipeline from a noise source family. */
template <int W, int H, typename SourceT>
class NoiseSourceProbe
    : public Pullback::ComposedEffect<
          W, H, NoiseSourceProbe<W, H, SourceT>,
          NoiseSourceProbeParams<SourceT>, NoiseSourceProbeSpec,
          PaletteHarmony::TRIADIC, Pullback::HueMode::NONE,
          Pullback::Color::BrightnessEnvelope::NONE, false> {
public:
  using Params = NoiseSourceProbeParams<SourceT>;
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"noise"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr int SOURCE_NOISE_SEED = 1337;

  static constexpr Params initial_params() {
    Params value;
    value.source.noise_scale = 8.0f;
    value.source.noise_contrast = 1.0f;
    value.source.noise_time_rate = 1.0f / 128.0f;
    return value;
  }
};

/**
 * @brief Both noise source families reach the derivation path.
 * @details Each owns a source-noise field without the effect naming
 * HasSourceNoise; the plane-domain and sphere-domain contours light the
 * frame and differ from each other on the same seed and parameters.
 */
inline void test_composed_noise_sources() {
  using Projected =
      NoiseSourceProbe<SMALL_W, SMALL_H, Pullback::ProjectedNoiseSourceParams>;
  using Spherical =
      NoiseSourceProbe<SMALL_W, SMALL_H, Pullback::SphericalNoiseSourceParams>;
  static_assert(Projected::HAS_SOURCE_NOISE);
  static_assert(Spherical::HAS_SOURCE_NOISE);

  std::array<Pixel, SMALL_W * SMALL_H> projected{};
  size_t lit = 0;
  {
    reset_effect_globals();
    Projected effect;
    effect.init();
    HS_EXPECT_TRUE(effect.getParameters().find("Source Noise Scale") !=
                   nullptr);
    HS_EXPECT_TRUE(effect.getParameters().find("Source Noise Speed") !=
                   nullptr);
    effect.draw_frame();
    effect.advance_display();
    for (int i = 0; i < SMALL_W * SMALL_H; ++i) {
      projected[i] = effect.display_buffer()[i];
      lit += projected[i].r != 0 || projected[i].g != 0 || projected[i].b != 0;
    }
  }
  HS_EXPECT_GT(lit, size_t(0));

  size_t spherical_lit = 0;
  size_t differing = 0;
  {
    reset_effect_globals();
    Spherical effect;
    effect.init();
    effect.draw_frame();
    effect.advance_display();
    for (int i = 0; i < SMALL_W * SMALL_H; ++i) {
      const Pixel &pixel = effect.display_buffer()[i];
      spherical_lit += pixel.r != 0 || pixel.g != 0 || pixel.b != 0;
      differing += pixel.r != projected[i].r || pixel.g != projected[i].g ||
                   pixel.b != projected[i].b;
    }
  }
  HS_EXPECT_GT(spherical_lit, size_t(0));
  HS_EXPECT_GT(differing, size_t(0));
}
/** @brief Parameter set of the ChoreographedEffect probes: one animated float. */
struct ChoreoProbeParams {
  float level = 0.0f;
};

/**
 * @brief Probe pinning the Segue::Preset::Fade envelope loop.
 * @details Records every set_preset_opacity sample the base's sprite feeds it,
 * so the test can see the envelope reach both edges and the advance timer
 * re-arm the loop.
 */
template <int W, int H>
class FadeChoreoProbe
    : public ChoreographedEffect<FadeChoreoProbe<W, H>, ChoreoProbeParams> {
  using Choreography =
      ChoreographedEffect<FadeChoreoProbe<W, H>, ChoreoProbeParams>;
  friend Choreography;

public:
  using Params = ChoreoProbeParams;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr int PRESET_FRAMES = 12;
  static constexpr Segue::Preset::Fade PRESET_SEGUE{PRESET_FRAMES, 4};
  static constexpr uint16_t PRESET_DWELL_FRAMES = PRESET_FRAMES;
  static constexpr std::array<PresetEntry<Params>, 3> PRESETS = {
      {{{0.25f}}, {{0.5f}}, {{0.75f}}}};

  static bool valid_params(const Params &p) {
    return p.level >= 0.0f && p.level <= 1.0f;
  }

  FadeChoreoProbe() : Choreography(W, H) {}

  void init() override { this->begin_choreography(); }

  void draw_frame() override {
    Canvas canvas(*this);
    this->timeline.step(canvas);
    this->step_choreography();
  }

  /** @brief Live level, which a preset change snaps. */
  float level() const { return this->params.level; }

  int opacity_samples = 0;    /**< Envelope samples the sprite fed. */
  float min_opacity = 2.0f;   /**< Lowest envelope sample seen. */
  float max_opacity = -1.0f;  /**< Highest envelope sample seen. */
  float last_opacity = -1.0f; /**< Most recent envelope sample. */

private:
  void set_preset_opacity(float value) {
    ++opacity_samples;
    last_opacity = value;
    min_opacity = std::min(min_opacity, value);
    max_opacity = std::max(max_opacity, value);
  }
};

/**
 * @brief Probe pinning the Segue::Preset::Lerp transition hooks.
 * @details Counts transition_armed and blend_params calls and writes the blend
 * itself, so a cancelled crossfade shows up as a blend that stops writing.
 */
template <int W, int H>
class LerpChoreoProbe
    : public ChoreographedEffect<LerpChoreoProbe<W, H>, ChoreoProbeParams> {
  using Choreography =
      ChoreographedEffect<LerpChoreoProbe<W, H>, ChoreoProbeParams>;
  friend Choreography;

public:
  using Params = ChoreoProbeParams;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr Segue::Preset::Lerp PRESET_SEGUE{8, ease_linear,
                                                    /*pausable=*/false};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 12;
  static constexpr std::array<PresetEntry<Params>, 2> PRESETS = {
      {{{0.0f}}, {{1.0f}}}};

  static bool valid_params(const Params &p) {
    return p.level >= 0.0f && p.level <= 1.0f;
  }

  LerpChoreoProbe() : Choreography(W, H) {}

  void init() override {
    this->register_animated_param("Level", &this->params.level, 0.0f, 1.0f);
    this->begin_choreography();
  }

  void draw_frame() override {
    Canvas canvas(*this);
    this->timeline.step(canvas);
    this->step_choreography();
  }

  /** @brief Live level, which the blend rewrites every frame it runs. */
  float level() const { return this->params.level; }

  int armed_count = 0;        /**< transition_armed calls seen. */
  float armed_target = -1.0f; /**< Level the last arming named. */
  int blend_calls = 0;        /**< blend_params calls seen. */

private:
  void transition_armed(const Params &target) {
    ++armed_count;
    armed_target = target.level;
  }

  void blend_params(float progress) {
    ++blend_calls;
    this->params.level = hs::lerp(this->transition.from.level,
                                  this->transition.to.level, progress);
  }
};

/** @brief Runs @p frames whole frames of @p effect, buffer swap included. */
template <typename FX> void run_probe_frames(FX &effect, int frames) {
  for (int frame = 0; frame < frames; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
}

/**
 * @brief Pins the Segue::Preset::Fade preset-choreography loop.
 * @details The envelope has to reach both edges, the advance timer has to
 * re-arm it every PRESET_FRAMES, and each advance has to snap the preset the
 * envelope just faded through.
 */
inline void test_choreography_fade_envelope() {
  using FX = FadeChoreoProbe<SMALL_W, SMALL_H>;
  constexpr int FRAMES = FX::PRESET_FRAMES;

  reset_effect_globals();
  FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), size_t{3});
  HS_EXPECT_EQ(effect.opacity_samples, 0);

  // One sprite per preset, one envelope sample per frame, and the advance lands
  // on the sprite's last frame.
  run_probe_frames(effect, FRAMES - 1);
  HS_EXPECT_EQ(effect.opacity_samples, FRAMES - 1);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  HS_EXPECT_EQ(effect.level(), 0.25f);
  HS_EXPECT_LE(effect.min_opacity, 0.5f);
  HS_EXPECT_GE(effect.max_opacity, 0.99f);
  HS_EXPECT_LE(effect.max_opacity, 1.0f);

  run_probe_frames(effect, 1);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
  HS_EXPECT_EQ(effect.level(), 0.5f);

  // The loop re-arms itself, so the cadence survives a full wrap of the table.
  run_probe_frames(effect, 2 * FRAMES);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  HS_EXPECT_EQ(effect.level(), 0.25f);

  // Paused mid-envelope: the sprite holds its phase and the advance never
  // fires, so the preset the envelope faded in stays up.
  run_probe_frames(effect, 3);
  const float held = effect.last_opacity;
  effect.setAnimationsPaused(true);
  run_probe_frames(effect, 2 * FRAMES);
  HS_EXPECT_EQ(effect.last_opacity, held);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  HS_EXPECT_EQ(effect.level(), 0.25f);
}

/**
 * @brief Pins the two Segue::Preset::Lerp transition hooks.
 * @details transition_armed fires once per automatic crossfade, with the
 * incoming preset; animated_parameter_written ends the crossfade in flight, so
 * the blend stops rewriting the value the write just landed.
 */
inline void test_choreography_lerp_transition_hooks() {
  using FX = LerpChoreoProbe<SMALL_W, SMALL_H>;
  reset_effect_globals();
  FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.armed_count, 0);
  HS_EXPECT_EQ(effect.level(), 0.0f);

  // Retire the dwell; the next frame arms the crossfade to PRESETS[1].
  run_probe_frames(effect, FX::PRESET_DWELL_FRAMES);
  HS_EXPECT_EQ(effect.armed_count, 1);
  HS_EXPECT_EQ(effect.armed_target, 1.0f);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});

  // Mid-crossfade the level is strictly between the endpoints.
  run_probe_frames(effect, FX::PRESET_SEGUE.frames / 2);
  HS_EXPECT_GT(effect.blend_calls, 0);
  HS_EXPECT_GT(effect.level(), 0.0f);
  HS_EXPECT_LT(effect.level(), 1.0f);

  // The manual write cancels the crossfade: the lerp keeps stepping (the policy
  // is unpausable) but blend_params is never called again.
  const int blends = effect.blend_calls;
  HS_EXPECT_EQ(effect.updateParameter("Level", 0.3f), ParamSetResult::APPLIED);
  run_probe_frames(effect, 2 * FX::PRESET_SEGUE.frames);
  HS_EXPECT_EQ(effect.blend_calls, blends);
  HS_EXPECT_EQ(effect.level(), 0.3f);
  HS_EXPECT_EQ(effect.armed_count, 1);
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
  test_composed_document_values();
  test_composed_derivation_reach();
  test_composed_projection_walk_storage();
  test_composed_periodic_ripple_surface();
  test_composed_noise_sources();
  test_choreography_fade_envelope();
  test_choreography_lerp_transition_hooks();
  return fixture.result();
}

} // namespace composed_effect_tests
} // namespace hs_test
