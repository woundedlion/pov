/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Base-contract sweep over every Pullback::ComposedEffect specialization in
 * HS_SHADER_PRODUCT_GROUP.
 */
#pragma once

#include "core/math/mobius.h"
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

#include "targets/effects.h"
#include "core/render/pullback/catalog_export.h"
#include "core/render/pullback/composed_effect.h"
#include "core/render/pullback/operators/table.h"
#include "tests/test_effects.h" // reset_effect_globals, SMALL_W/SMALL_H
#include "tests/test_fixture.h"
#include "tests/composed_frame_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace composed_effect_tests {

using CollidingWarpNames = Pullback::ComposedDetail::ResourceList<
    Pullback::ParameterResource<"outer_warp", Pullback::WaveShearParams,
                                Pullback::ResourceKind::WARP>,
    Pullback::ParameterResource<"inner_warp", Pullback::VectorNoiseParams,
                                Pullback::ResourceKind::WARP>>;
static_assert(
    Pullback::ComposedDetail::qualified<
        Pullback::ParameterResource<"outer_warp", Pullback::WaveShearParams,
                                    Pullback::ResourceKind::WARP>,
        CollidingWarpNames>());
static_assert(
    Pullback::ComposedDetail::qualified<
        Pullback::ParameterResource<"inner_warp", Pullback::VectorNoiseParams,
                                    Pullback::ResourceKind::WARP>,
        CollidingWarpNames>());

using effects_tests::preset_params_or_initial;
using effects_tests::reset_effect_globals;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;

struct MismatchedCoverageSpec : AshCloudSpec {
  static constexpr auto COVERAGE =
      Pullback::ProjectionCoverageMode::WEIGHT_SQUARED;
};
using CoverageFixture = AshCloud<SMALL_W, SMALL_H>;
static_assert(Pullback::ComposedDetail::PipelineMetadata<
              AshCloudSpec, CoverageFixture::Binding,
              CoverageFixture::RenderPipeline>::COVERAGE_MATCHES);
static_assert(!Pullback::ComposedDetail::PipelineMetadata<
              MismatchedCoverageSpec, CoverageFixture::Binding,
              CoverageFixture::RenderPipeline>::COVERAGE_MATCHES);

using TrackedFixture = CosmicEyeball<SMALL_W, SMALL_H>;
static_assert(Pullback::ComposedDetail::PipelineMetadata<
              CosmicEyeballSpec, TrackedFixture::Binding,
              TrackedFixture::RenderPipeline>::PATH_TRACKED);
static_assert(!Pullback::ComposedDetail::PipelineMetadata<
              AshCloudSpec, CoverageFixture::Binding,
              CoverageFixture::RenderPipeline>::PATH_TRACKED);

template <typename FX> struct ComposedTraits {
  using Params = typename FX::Params;
  using Spec = typename FX::Spec;
  static constexpr PaletteHarmony HARMONY = Spec::HARMONY;
  static constexpr Pullback::HueMode HUE = Spec::HUE;
  static constexpr Pullback::Color::BrightnessEnvelope BRIGHTNESS =
      Spec::BRIGHTNESS;
  static constexpr Pullback::SurfacePlacement SURFACE_PLACEMENT =
      Spec::SURFACE_PLACEMENT;
};

template <typename FX> using TraitsOf = ComposedTraits<FX>;

/** @brief Color sliders the base registers for every specialization. */
constexpr const char *UNGATED_COLOR_SLIDERS[] = {
    "Palette Chroma",     "Palette Mapping",         "Mapping Frequency",
    "Mapping Phase",      "Phase Oscillation Depth", "Phase Oscillation Speed",
    "Opacity at Value 0", "Opacity at Value 1"};

/** @brief Lens sliders a Mobius parameter family adds. */
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
 * @details Lets a slider's authored range be compared against the range the
 * snapshot validator enforces.
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
    {"Brightness Bottom", "brightness-bottom"},
    {"Brightness Top", "brightness-top"},
    {"Opacity at Value 0", "value-opacity-low"},
    {"Opacity at Value 1", "value-opacity-high"}};

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
  if constexpr (std::is_void_v<Family>)
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
template <typename FX, typename Family>
inline void verify_family_sliders(const ParamList &params, bool registered) {
  if constexpr (Pullback::HasFields<Family>)
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
  constexpr bool expected = !std::is_void_v<Family>;
  const ParamDef *def = params.find(slot_name);
  HS_EXPECT_EQ(def != nullptr, expected);
  if constexpr (!std::is_void_v<Family>)
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
inline void verify_mobius_equal(const math::MobiusParams &actual,
                                const math::MobiusParams &expected) {
  constexpr math::Complex math::MobiusParams::*COEFFICIENTS[] = {
      &math::MobiusParams::a, &math::MobiusParams::b, &math::MobiusParams::c,
      &math::MobiusParams::d};
  for (math::Complex math::MobiusParams::*coefficient : COEFFICIENTS) {
    HS_EXPECT_EQ(bits((actual.*coefficient).re),
                 bits((expected.*coefficient).re));
    HS_EXPECT_EQ(bits((actual.*coefficient).im),
                 bits((expected.*coefficient).im));
  }
}

/** @brief Bitwise-compares a whole parameter set, family by family. */
template <typename Params>
inline void verify_params_equal(const Params &actual, const Params &expected) {
  actual.visit([&]<typename Resource>(const auto &family) {
    const auto &other = expected.template get<Resource::KEY>();
    if constexpr (Pullback::HasFields<typename Resource::Family>)
      verify_family_equal(family, other);
    if constexpr (Resource::KIND == Pullback::ResourceKind::LENS)
      verify_mobius_equal(family.mobius, other.mobius);
    if constexpr (Resource::KIND == Pullback::ResourceKind::COLOR)
      HS_EXPECT_EQ(static_cast<int>(family.palette_mapping),
                   static_cast<int>(other.palette_mapping));
  });
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
 */
template <typename FX, typename Resource>
inline void
verify_family_rejection(FX &effect,
                        const typename FX::ParameterSnapshot &captured) {
  using Family = typename Resource::Family;
  HS_CONTEXT(Resource::KEY.text);
  for (const auto &field : Family::FIELDS) {
    HS_CONTEXT(field.id);
    const float span = field.max - field.min;
    const float poisons[] = {std::numeric_limits<float>::quiet_NaN(),
                             field.min - span - 1.0f, field.max + span + 1.0f};
    for (float poison : poisons) {
      typename FX::ParameterSnapshot snapshot = captured;
      snapshot.params.template get<Resource::KEY>().*(field.member) = poison;
      HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
      verify_params_equal(effect.serialize_parameters().params,
                          captured.params);
    }
  }
}

/**
 * @brief Pins the slider set one specialization's init() registers.
 * @tparam E Composed effect class template.
 * @param name Effect name, for the failure context.
 * @details The registered count is derived from the parameter families.
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

  verify_family_sliders<FX, typename Params::template Family<"source">>(params,
                                                                        true);
  verify_family_sliders<FX, Pullback::ProjectionParams>(params, true);
  verify_family_sliders<FX, typename Params::template Family<"surface">>(params,
                                                                         true);
  verify_family_sliders<FX, typename Params::template Family<"value">>(params,
                                                                       true);
  verify_warp_slot<FX, typename Params::template Family<"outer_warp">>(
      params, "Planar Warp 1 Speed");
  verify_warp_slot<FX, typename Params::template Family<"inner_warp">>(
      params, "Planar Warp 2 Speed");

  for (const char *slider : UNGATED_COLOR_SLIDERS)
    HS_EXPECT_TRUE(params.find(slider) != nullptr);
  constexpr bool mobius =
      std::is_same_v<typename Params::template Family<"lens">,
                     Pullback::MobiusLensParams>;
  for (const char *slider : MOBIUS_SLIDERS)
    HS_EXPECT_EQ(params.find(slider) != nullptr, mobius);
  constexpr bool brightness =
      Traits::BRIGHTNESS != Pullback::Color::BrightnessEnvelope::NONE;
  constexpr bool hue_shift = Traits::HUE != Pullback::HueMode::NONE;
  HS_EXPECT_EQ(params.find("Hue Shift Amount") != nullptr, hue_shift);
  constexpr bool hue_noise = Traits::HUE == Pullback::HueMode::NOISE;
  HS_EXPECT_EQ(params.find("Brightness Bottom") != nullptr, brightness);
  HS_EXPECT_EQ(params.find("Brightness Top") != nullptr, brightness);
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
    HS_EXPECT_EQ(
        mapping->option_count,
        static_cast<int>(std::size(Pullback::Interp::Op::PALETTE_MAPPING_IDS)));
    HS_EXPECT_TRUE(mapping->options != nullptr);
    HS_EXPECT_TRUE(mapping->export_options != nullptr);
  }

  const size_t expected =
      named_field_count<FX, typename Params::template Family<"source">>() +
      named_field_count<FX, Pullback::ProjectionParams>() +
      named_field_count<FX, typename Params::template Family<"surface">>() +
      named_field_count<FX, typename Params::template Family<"value">>() +
      warp_slot_count<FX, typename Params::template Family<"outer_warp">>() +
      warp_slot_count<FX, typename Params::template Family<"inner_warp">>() +
      (mobius ? std::size(MOBIUS_SLIDERS) : size_t{0}) +
      std::size(UNGATED_COLOR_SLIDERS) + (brightness ? size_t{2} : size_t{0}) +
      (hue_noise ? size_t{2} : size_t{0}) + (hue_shift ? size_t{1} : size_t{0});
  HS_EXPECT_EQ(params.size(), expected);
  HS_EXPECT_EQ(params.size(), params.capacity());
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
  moved.params.visit([]<typename Resource>(auto &family) {
    if constexpr (Pullback::HasFields<typename Resource::Family>)
      fill_midpoints(family);
  });
  if constexpr (Params::template HAS<"lens">) {
    moved.params.template get<"lens">().mobius.a = {1.0f, 0.2f};
    moved.params.template get<"lens">().mobius.b = {0.3f, 0.4f};
    moved.params.template get<"lens">().mobius.c = {0.5f, 0.6f};
    moved.params.template get<"lens">().mobius.d = {0.8f, -0.1f};
  }
  constexpr size_t MAPPINGS =
      std::size(Pullback::Interp::Op::PALETTE_MAPPING_IDS);
  moved.params.template get<"color">().palette_mapping =
      static_cast<Pullback::Color::PaletteMapping>(
          (static_cast<uint8_t>(
               captured.params.template get<"color">().palette_mapping) +
           1) %
          MAPPINGS);
  HS_EXPECT_TRUE(effect.restore_parameters(moved));
  verify_params_equal(effect.serialize_parameters().params, moved.params);
  HS_EXPECT_TRUE(effect.restore_parameters(captured));

  typename FX::ParameterSnapshot bumped = captured;
  bumped.schema_version += 1;
  HS_EXPECT_FALSE(effect.restore_parameters(bumped));
  verify_params_equal(effect.serialize_parameters().params, captured.params);

  captured.params.visit([&]<typename Resource>(const auto &) {
    if constexpr (Pullback::HasFields<typename Resource::Family>)
      verify_family_rejection<FX, Resource>(effect, captured);
  });

  typename FX::ParameterSnapshot mapping = captured;
  mapping.params.template get<"color">().palette_mapping =
      static_cast<Pullback::Color::PaletteMapping>(MAPPINGS);
  HS_EXPECT_FALSE(effect.restore_parameters(mapping));
  verify_params_equal(effect.serialize_parameters().params, captured.params);

  if constexpr (Params::template HAS<"lens">) {
    typename FX::ParameterSnapshot lens = captured;
    lens.params.template get<"lens">().mobius.a.re =
        std::numeric_limits<float>::quiet_NaN();
    HS_EXPECT_FALSE(effect.restore_parameters(lens));
    verify_params_equal(effect.serialize_parameters().params, captured.params);
    lens = captured;
    lens.params.template get<"lens">().mobius.a.re =
        Pullback::MobiusLensParams::COEFFICIENT_LIMIT + 1.0f;
    HS_EXPECT_FALSE(effect.restore_parameters(lens));
    verify_params_equal(effect.serialize_parameters().params, captured.params);
    lens = captured;
    lens.params.template get<"lens">().mobius.c =
        lens.params.template get<"lens">().mobius.a;
    lens.params.template get<"lens">().mobius.d =
        lens.params.template get<"lens">().mobius.b;
    HS_EXPECT_FALSE(effect.restore_parameters(lens));
    verify_params_equal(effect.serialize_parameters().params, captured.params);
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

  if (FX::PRESET_IDS.size() > 1)
    HS_EXPECT_GT(Segue::Preset::frames(FX::preset_departure(0)), uint16_t{0});
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
 * @details Both endpoints must come back exactly, and no sample in between may
 * leave the ranges the snapshot validator enforces.
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
      Params b = preset_params_or_initial<FX>(to);
      if constexpr (FX::PRESET_IDS.size() == 1) {
        b.template get<"color">().opacity_low =
            a.template get<"color">().opacity_low == 0.0f ? 1.0f : 0.0f;
        HS_EXPECT_NEAR(Pullback::interpolate(a, b, 0.5f)
                           .template get<"color">()
                           .opacity_low,
                       0.5f * (a.template get<"color">().opacity_low +
                               b.template get<"color">().opacity_low),
                       1e-6f);
      }
      for (float progress : PROGRESS)
        HS_EXPECT_TRUE(FX::valid_params(Pullback::interpolate(a, b, progress)));
      verify_params_equal(Pullback::interpolate(a, b, 0.0f), a);
      verify_params_equal(Pullback::interpolate(a, b, 1.0f), b);
    }
  }
}

#include "tests/composed_effect/derivation_reach.h"
#include "tests/composed_effect/document_values.h"
#include "tests/composed_effect/roster_contracts.h"
#include "tests/composed_effect/stage_probes.h"
#include "tests/composed_effect/parameter_schema.h"

inline void test_flowers_longitude_seam() {
  using FX = KaleidoscopeFlowers<SMALL_W, SMALL_H>;
  for (size_t preset = 0; preset < FX::PRESET_IDS.size(); ++preset) {
    const auto params = FX::preset(preset).params;
    for (int step = 0; step <= 128; ++step) {
      const auto prepared = Pullback::Warp::prepare(
          params.template get<"inner_warp">(), step / 128.0f);
      for (float latitude : {-1.4f, -0.7f, 0.0f, 0.7f, 1.4f}) {
        const auto left = Pullback::Warp::mirror_tile_coords(
            math::Complex(-math::PI_F, latitude),
            params.template get<"inner_warp">(), prepared);
        const auto right = Pullback::Warp::mirror_tile_coords(
            math::Complex(math::PI_F, latitude),
            params.template get<"inner_warp">(), prepared);
        HS_EXPECT_NEAR(left.re, right.re, 1e-5f);
        HS_EXPECT_NEAR(left.im, right.im, 1e-5f);
      }
    }
  }
}

/** @brief AlienBrain holds each preset for its declared dwell. */
inline void test_alien_brain_preset_dwell() {
  using FX = AlienBrain<SMALL_W, SMALL_H>;
  reset_effect_globals();
  FX effect;
  effect.init();

  HS_EXPECT_EQ(effect.getPresetIndex(), size_t(0));
  for (uint16_t frame = 1; frame < FX::PRESET_DWELL_FRAMES; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t(0));

  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t(1));
}

inline void test_mobius_grid_circular_animation() {
  using FX = MobiusGrid<SMALL_W, SMALL_H>;
  reset_effect_globals();
  FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), size_t(2));
  HS_EXPECT_TRUE(FX::PRESET_IDS[1] == "mobius-grid-2");

  const math::MobiusParams initial =
      effect.serialize_parameters().params.template get<"lens">().mobius;
  effect.draw_frame();
  effect.advance_display();
  const math::MobiusParams animated =
      effect.serialize_parameters().params.template get<"lens">().mobius;
  HS_EXPECT_TRUE(animated.b.re != initial.b.re ||
                 animated.b.im != initial.b.im);
  HS_EXPECT_NEAR(animated.b.re * animated.b.re + animated.b.im * animated.b.im,
                 1.0f, 1e-5f);
  effect.draw_frame();
  effect.advance_display();
  const math::MobiusParams advanced =
      effect.serialize_parameters().params.template get<"lens">().mobius;
  HS_EXPECT_TRUE(advanced.b.re != animated.b.re ||
                 advanced.b.im != animated.b.im);
  effect.setAnimationsPaused(true);
  effect.draw_frame();
  effect.advance_display();
  const math::MobiusParams paused =
      effect.serialize_parameters().params.template get<"lens">().mobius;
  HS_EXPECT_EQ(paused.b.re, advanced.b.re);
  HS_EXPECT_EQ(paused.b.im, advanced.b.im);

  effect.setAnimationsPaused(false);
  effect.draw_frame();
  effect.advance_display();
  const math::MobiusParams resumed =
      effect.serialize_parameters().params.template get<"lens">().mobius;
  HS_EXPECT_TRUE(resumed.b.re != paused.b.re || resumed.b.im != paused.b.im);

  math::MobiusParams previous = resumed;
  const size_t initial_preset = effect.getPresetIndex();
  for (int frame = 0; frame < 1400; ++frame) {
    effect.draw_frame();
    effect.advance_display();
    const math::MobiusParams current =
        effect.serialize_parameters().params.template get<"lens">().mobius;
    HS_EXPECT_NEAR(current.b.re * current.b.re + current.b.im * current.b.im,
                   1.0f, 1e-5f);
    const float delta_re = current.b.re - previous.b.re;
    const float delta_im = current.b.im - previous.b.im;
    HS_EXPECT_TRUE(delta_re * delta_re + delta_im * delta_im < 0.02f);
    previous = current;
  }
  HS_EXPECT_NE(effect.getPresetIndex(), initial_preset);
  const auto current = effect.serialize_parameters().params;
  auto target = FX::preset(effect.getPresetIndex()).params;
  target.template get<"lens">().mobius = current.template get<"lens">().mobius;
  verify_params_equal(current, target);
}

/**
 * @brief Pins exact endpoints and weighted coordinates during mapping crossfades.
 */
inline void test_palette_mapping_crossfade() {
  using namespace Pullback::Color;
  constexpr float PHASE = 0.2f;
  const auto LINEAR = PaletteMappingWeights::single(PaletteMapping::LINEAR);
  const auto CUP = PaletteMappingWeights::single(PaletteMapping::CUP);
  for (const float PROGRESS : {0.0f, 1.0f}) {
    const auto weights = PaletteMappingWeights::lerp(LINEAR, CUP, PROGRESS);
    const auto MAPPING =
        PROGRESS == 0.0f ? PaletteMapping::LINEAR : PaletteMapping::CUP;
    HS_EXPECT_EQ(weights.exact, static_cast<uint8_t>(MAPPING));
    HS_EXPECT_EQ(palette_mapping_coordinate(PHASE, weights, 1.0f, 0.0f),
                 palette_mapping_coordinate(PHASE, MAPPING, 1.0f, 0.0f));
  }
  const auto midpoint = PaletteMappingWeights::lerp(LINEAR, CUP, 0.5f);
  HS_EXPECT_EQ(midpoint.exact, 0xff);
  HS_EXPECT_NEAR(palette_mapping_coordinate(PHASE, midpoint, 1.0f, 0.0f),
                 0.5f * PHASE + 0.5f * math::unit_cup(PHASE), 1e-6f);
  const auto bell = PaletteMappingWeights::single(PaletteMapping::BELL);
  HS_EXPECT_NEAR(
      palette_mapping_coordinate(
          PHASE, PaletteMappingWeights::lerp(CUP, bell, 0.5f), 1.0f, 0.0f),
      0.5f, 1e-6f);
  const auto unchanged = PaletteMappingWeights::lerp(CUP, CUP, 0.3f);
  HS_EXPECT_EQ(unchanged.exact, CUP.exact);
  HS_EXPECT_TRUE(unchanged.values == CUP.values);
}

/** @brief Module entry point for the composed-effect base contract. */
inline int run_composed_effect_tests() {
  ModuleFixture fixture("composed_effect");
  test_alien_brain_preset_dwell();
  test_mobius_grid_circular_animation();
  test_catalog_semantic_export();
  test_composed_hand_registered_families();
  test_composed_direct_surface_placement();
  test_composed_slider_registration();
  check_slider_registration<NoHueProbe>("NoHueProbe");
  test_mobius_frame_admission();
  test_mobius_captured_parameter_restore();
  test_mobius_automatic_departures_retain_lens();
  test_composed_snapshot_contract();
  test_parameter_layout_reorder();
  test_composed_parameter_schema_pins();
  test_composed_pentbright_lattice_scale();
  test_composed_preset_choreography();
  test_composed_log_positive_curve();
  test_composed_preset_interpolation();
  test_palette_mapping_crossfade();
  test_composed_document_values();
  test_flowers_longitude_seam();
  test_composed_derivation_reach();
  test_composed_keyed_parameters();
  test_composed_affine_named_source();
  test_composed_repeated_instances();
  test_composed_projection_walk_storage();
  test_composed_periodic_ripple_surface();
  test_composed_noise_sources();
  test_choreography_fade_departure();
  test_choreography_lerp_transition_hooks();
  test_choreography_lerp_pause_policy();
  return fixture.result();
}

} // namespace composed_effect_tests
} // namespace hs_test
