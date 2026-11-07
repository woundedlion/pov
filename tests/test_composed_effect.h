/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Base-contract sweep over every Pullback::ComposedEffect specialization in
 * HS_SHADER_PRODUCT_GROUP.
 */
#pragma once

#include "core/math/mobius.h"
#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "targets/effects.h"
#include "core/render/pullback/catalog_export.h"
#include "core/render/pullback/composed_effect.h"
#include "core/render/pullback/operators/table.h"
#include "tests/test_fixture.h"
#include "tests/composed_frame_fixture.h"
#include "tests/pixel_test_util.h"
#include "tests/test_harness.h"
#include "tests/composed_effect/catalog_policies.h"
#include "tests/composed_effect/keyed_rotations.h"
#include "tests/vec_test_util.h"
#include "targets/wasm/document_parameters.h"

namespace hs_test::composed_chain_tests {
template <typename FX>
void verify_export(size_t preset, int channel_tolerance = 1,
                   float alpha_tolerance = 1e-6f);
}

namespace hs_test {
namespace composed_effect_tests {

#if HS_TEST_FAST_MATH
inline constexpr int DOCUMENT_CHANNEL_TOLERANCE = 2;
#else
inline constexpr int DOCUMENT_CHANNEL_TOLERANCE = 1;
#endif

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

using effects_tests::preset_params;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;

#include "tests/composed_effect/generated_catalog_assignments.h"

using CoverageFixture = AshCloud<SMALL_W, SMALL_H>;
static_assert(CoverageFixture::Metadata::COVERAGE ==
              Pullback::ProjectionCoverageMode::WEIGHT);
using TrackedFixture = CosmicEyeball<SMALL_W, SMALL_H>;
static_assert(TrackedFixture::Metadata::PATH_TRACKED);
static_assert(CoverageFixture::Metadata::PATH_TRACKED);
static_assert(CoverageFixture::Metadata::OUTER_CAMERA);
static_assert(!Involution<SMALL_W, SMALL_H>::Metadata::OUTER_CAMERA);

template <typename FX> struct TraitsOf {
  using Spec = typename FX::Spec;
  static constexpr PaletteHarmony HARMONY = Spec::HARMONY;
  static constexpr Pullback::Color::BrightnessEnvelope BRIGHTNESS =
      FX::Metadata::BRIGHTNESS;
  static constexpr Pullback::SurfacePlacement SURFACE_PLACEMENT =
      FX::Metadata::SURFACE_PLACEMENT;
};

/** @brief Color sliders the base registers for every specialization. */
constexpr const char *UNGATED_COLOR_SLIDERS[] = {
    "Palette Chroma",     "Palette Mapping",         "Mapping Frequency",
    "Mapping Phase",      "Phase Oscillation Depth", "Phase Oscillation Speed",
    "Opacity at Value 0", "Opacity at Value 1"};

/** @brief Lens sliders a Mobius parameter family adds. */
constexpr const char *MOBIUS_SLIDERS[] = {"Spin X", "Spin Y", "Spin Z",
                                          "Push X", "Push Y", "Push Z"};

/** @brief A hand-registered color slider and the field descriptor it writes. */
struct ColorSliderBinding {
  const char *slider;
  const char *field_id;
};

/**
 * @brief Ranged color sliders paired with the field descriptor each writes.
 * @details Lets a slider's authored range be compared against the range the
 * snapshot validator enforces.
 */
constexpr ColorSliderBinding COLOR_SLIDER_BINDINGS[] = {
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

/** @brief Instance key, path tracking, next stage and form of an affine warp
 * stage. */
template <typename Stage> struct AffineStageTracking {
  static constexpr bool AFFINE = false; ///< Not an affine warp stage.
};
/** @brief An affine warp stage over a keyed warp provider. */
template <typename B, Pullback::ResourceKey Key, typename Family, bool Track,
          Pullback::ResourceKey SourceKey, Pullback::AffineForm Form>
struct AffineStageTracking<Pullback::Stage::Warp<Pullback::Warp::AffineFrame<
    Pullback::WarpProvider<B, Key, Family, Track, SourceKey, Form>>>> {
  static constexpr bool AFFINE = true;                ///< An affine warp stage.
  static constexpr std::string_view KEY = Key.view(); ///< Instance key.
  static constexpr bool TRACK = Track;                ///< Tracks path length.
  static constexpr auto SOURCE_KEY = SourceKey;       ///< Next stage's key.
  static constexpr Pullback::AffineForm FORM = Form;  ///< Compiled parts.
};

/** @brief Whether @p FX's affine stage keyed @p Key reads field @p id.
 * @return True for any field of a key no affine stage carries. */
template <typename FX, Pullback::ResourceKey Key>
constexpr bool affine_field_read(std::string_view id) {
  using Leaves = Pullback::ComposedDetail::PipelineLeaves<
      typename FX::Spec::template Pipeline<
          Pullback::ComposedDetail::DiscoveryBinding>>;
  bool read = true;
  [&]<size_t... I>(std::index_sequence<I...>) {
    (
        [&]<typename Stage>() {
          using Tracking = AffineStageTracking<Stage>;
          if constexpr (Tracking::AFFINE)
            if (Tracking::KEY == Key.view()) {
              using Source =
                  typename FX::Params::template Family<Tracking::SOURCE_KEY>;
              constexpr bool LATTICE = [] {
                if constexpr (std::is_void_v<Source>)
                  return false;
                else
                  return requires(const Source &source) {
                    translation_lattice(source);
                  };
              }();
              using Form = Pullback::AffineForm;
              const auto has = [](Form part) {
                return Pullback::Warp::includes(Tracking::FORM, part);
              };
              if (id == "speed")
                read = Tracking::FORM != Form::IDENTITY;
              else if (id == "rotation-rate")
                read = has(Form::ROTATES);
              else if (id == "translation-x" || id == "translation-y")
                read = has(Form::TRANSLATES);
              else if (id == "scale-x" || id == "scale-y" || id == "shear")
                read = has(Form::SCALES_SHEARS);
              else if (id == "repeat-cells")
                read = has(Form::TRANSLATES) || (Tracking::TRACK && LATTICE);
              else if (id == "lattice-period")
                read = has(Form::TRANSLATES) && !LATTICE;
            }
        }.template operator()<typename Leaves::template StageAt<I>>(),
        ...);
  }(std::make_index_sequence<Leaves::STAGE_COUNT>{});
  return read;
}

/** @brief Whether @p FX opens a gated field's slider. */
template <typename FX> constexpr bool gate_open(Pullback::FieldGate gate) {
  switch (gate) {
  case Pullback::FieldGate::ALWAYS:
  case Pullback::FieldGate::AFFINE_ROTATES:
  case Pullback::FieldGate::AFFINE_TRANSLATES:
  case Pullback::FieldGate::AFFINE_SCALES_SHEARS:
  case Pullback::FieldGate::AFFINE_REPEAT:
  case Pullback::FieldGate::AFFINE_FALLBACK_LATTICE:
    return true;
  case Pullback::FieldGate::OUTER_CAMERA:
    return FX::Metadata::OUTER_CAMERA;
  case Pullback::FieldGate::ANIMATED_PROJECTION:
    return FX::ANIMATED_PROJECTION;
  case Pullback::FieldGate::CENTRAL_MERIDIAN:
    return Pullback::uses_central_meridian(FX::Metadata::PROJECTION);
  case Pullback::FieldGate::SINGULARITY_FADE:
    return Pullback::uses_singularity_fade(FX::Metadata::PROJECTION);
  case Pullback::FieldGate::COORDINATE_SCALE:
    return FX::Metadata::PROJECTION == Pullback::ProjectionKind::PEIRCE ||
           FX::Metadata::PROJECTION ==
               Pullback::ProjectionKind::PEIRCE_SQUARE_FAST ||
           FX::Metadata::PROJECTION == Pullback::ProjectionKind::BONNE ||
           FX::Metadata::PROJECTION == Pullback::ProjectionKind::AIROCEAN ||
           Pullback::is_polyhedral_atlas(FX::Metadata::PROJECTION);
  case Pullback::FieldGate::LAYOUT_SCROLL:
    return FX::Metadata::PROJECTION == Pullback::ProjectionKind::PEIRCE;
  case Pullback::FieldGate::STANDARD_PARALLEL:
    return FX::Metadata::PROJECTION == Pullback::ProjectionKind::BONNE;
  case Pullback::FieldGate::AREA_CORRECTION:
    return Pullback::is_polyhedral_atlas(FX::Metadata::PROJECTION);
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

/** @brief Whether @p FX registers a slider for a field of the @p Key
 * resource. */
template <typename FX, Pullback::ResourceKey Key, typename Field>
constexpr bool field_open(const Field &field) {
  return gate_open<FX>(field.gate) && affine_field_read<FX, Key>(field.id);
}

/** @brief Sliders a family contributes through register_fields(). */
template <typename FX, typename Family, Pullback::ResourceKey Key = "">
constexpr size_t named_field_count() {
  size_t count = 0;
  if constexpr (Pullback::HasFields<Family>)
    for (const auto &field : Family::FIELDS)
      if (field.name != nullptr && field_open<FX, Key>(field))
        ++count;
  return count;
}

/** @brief Sliders a warp slot contributes: the slot speed plus its own. */
template <typename FX, Pullback::ResourceKey Key>
constexpr size_t warp_slot_count() {
  using Family = typename FX::Params::template Family<Key>;
  if constexpr (std::is_void_v<Family>)
    return 0;
  else
    return (affine_field_read<FX, Key>("speed") ? 1 : 0) +
           named_field_count<FX, Family, Key>();
}

/**
 * @brief Checks one family's tabled fields against the registered sliders.
 * @param params The effect's registered parameter list.
 * @param registered Whether the base runs register_fields() on this family.
 * @details A named field must own exactly one slider, present on the gate its
 * descriptor names and carrying the descriptor's bounds.
 */
template <typename FX, typename Family, Pullback::ResourceKey Key = "">
inline void verify_family_sliders(const ParamList &params, bool registered) {
  if constexpr (Pullback::HasFields<Family>)
    for (const auto &field : Family::FIELDS) {
      if (field.name == nullptr)
        continue;
      HS_CONTEXT(field.id);
      const bool expected = registered && field_open<FX, Key>(field);
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
template <typename FX, Pullback::ResourceKey Key>
inline void verify_warp_slot(const ParamList &params, const char *slot_name) {
  using Family = typename FX::Params::template Family<Key>;
  HS_CONTEXT(slot_name);
  constexpr bool present = !std::is_void_v<Family>;
  const ParamDef *def = params.find(slot_name);
  constexpr bool speed = present && affine_field_read<FX, Key>("speed");
  HS_EXPECT_EQ(def != nullptr, speed);
  if constexpr (present)
    if (def != nullptr) {
      HS_EXPECT_EQ(bits(def->min), bits(Family::FIELDS[0].min));
      HS_EXPECT_EQ(bits(def->max), bits(Family::FIELDS[0].max));
    }
  verify_family_sliders<FX, Family, Key>(params, present);
}

/** @brief Bitwise-compares every tabled field of one parameter family. */
template <Pullback::HasFields Family>
inline void verify_family_equal(const Family &actual, const Family &expected) {
  for (const auto &field : Family::FIELDS) {
    HS_CONTEXT(field.id);
    HS_EXPECT_EQ(bits(actual.*(field.member)), bits(expected.*(field.member)));
  }
  if constexpr (requires { Family::TOPOLOGY; })
    for (const auto &field : Family::TOPOLOGY)
      HS_EXPECT_EQ(actual.*field.member, expected.*field.member);
  if constexpr (requires { actual.modulation; })
    HS_EXPECT_TRUE(actual.modulation == expected.modulation);
  if constexpr (requires { actual.burst_size; }) {
    HS_EXPECT_EQ(actual.burst_size, expected.burst_size);
    HS_EXPECT_EQ(actual.random_center, expected.random_center);
  }
  if constexpr (requires { actual.terms; }) {
    for (size_t i = 0; i < actual.terms.size(); ++i) {
      const auto &a = actual.terms[i];
      const auto &b = expected.terms[i];
      HS_EXPECT_EQ(a.type, b.type);
      HS_EXPECT_EQ(a.kind, b.kind);
      HS_EXPECT_EQ(a.l, b.l);
      HS_EXPECT_EQ(a.m, b.m);
      HS_EXPECT_EQ(bits(a.amplitude), bits(b.amplitude));
      HS_EXPECT_EQ(bits(a.phase), bits(b.phase));
    }
  }
  if constexpr (std::is_same_v<Family, Pullback::Color::AuthoredColorParams>) {
    HS_EXPECT_EQ(actual.palette_mode, expected.palette_mode);
    HS_EXPECT_EQ(actual.palette_type, expected.palette_type);
    HS_EXPECT_EQ(actual.palette_harmony, expected.palette_harmony);
    HS_EXPECT_EQ(actual.palette_lightness_curve,
                 expected.palette_lightness_curve);
    HS_EXPECT_EQ(actual.palette_name, expected.palette_name);
    HS_EXPECT_EQ(actual.palette_name_b, expected.palette_name_b);
    HS_EXPECT_EQ(actual.brightness_envelope, expected.brightness_envelope);
  }
}

/** @brief Bitwise-compares a whole parameter set, family by family. */
template <typename Params>
inline void verify_params_equal(const Params &actual, const Params &expected) {
  actual.visit([&]<typename Resource>(const auto &family) {
    const auto &other = expected.template get<Resource::KEY>();
    if constexpr (Pullback::HasFields<typename Resource::Family>)
      verify_family_equal(family, other);
    if constexpr (Resource::KIND == Pullback::ResourceKind::COLOR)
      HS_EXPECT_EQ(static_cast<int>(family.palette_mapping),
                   static_cast<int>(other.palette_mapping));
  });
}

/**
 * @brief Moves every tabled field of a family to an in-range value that
 * differs from its captured value.
 * @param family Family to move.
 * @param captured Values each field must move off.
 * @details The midpoint, or the lower quartile where the midpoint is the
 * captured value.
 */
template <Pullback::HasFields Family>
inline void move_off_captured(Family &family, const Family &captured) {
  for (const auto &field : Family::FIELDS) {
    float value = 0.5f * (field.min + field.max);
    if (bits(value) == bits(captured.*(field.member)))
      value = field.min + 0.25f * (field.max - field.min);
    family.*(field.member) = value;
    HS_EXPECT_NE(bits(value), bits(captured.*(field.member)));
  }
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

/** @brief Authored fields register one control at their native storage address. */
template <typename FX>
inline void check_authored_slider_registration(const FX &effect) {
  const auto &params = effect.getParameters();
  size_t expected = 0;
  const auto require_target = [&](const void *target, float low, float high,
                                  bool selector = false,
                                  bool readonly = false) {
    const ParamDef *found = nullptr;
    size_t matches = 0;
    for (const auto &entry : params)
      if (entry.target == target) {
        found = &entry;
        ++matches;
      }
    HS_EXPECT_EQ(matches, 1u);
    if (found) {
      HS_EXPECT_EQ(found->min, low);
      HS_EXPECT_EQ(found->max, high);
      HS_EXPECT_EQ(found->options != nullptr, selector);
      HS_EXPECT_TRUE(found->animated);
      HS_EXPECT_EQ(found->readonly, readonly);
    }
    ++expected;
  };
  ComposedWhiteBox<FX>::params(effect).visit([&]<typename Resource>(
                                                 const auto &family) {
    using Family = typename Resource::Family;
    for (const auto &field : Family::FIELDS) {
      if constexpr (Resource::KIND != Pullback::ResourceKind::COLOR)
        if ((!field.name && !(Resource::KIND == Pullback::ResourceKind::WARP &&
                              std::string_view(field.id) == "speed")) ||
            !field_open<FX, Resource::KEY>(field))
          continue;
      require_target(&(family.*field.member), field.min, field.max);
    }
    if constexpr (Resource::KIND != Pullback::ResourceKind::COLOR &&
                  requires { Family::TOPOLOGY; })
      for (const auto &field : Family::TOPOLOGY) {
        if (!Pullback::ComposedDetail::resource_topology_open<typename FX::Spec,
                                                              Resource>(
                field.id))
          continue;
        size_t options = field.value_count;
        if constexpr (Resource::KIND == Pullback::ResourceKind::PROJECTION)
          options = Pullback::projection_topology_options(
                        FX::Metadata::PROJECTION, field.id,
                        {field.value_ids, field.value_count})
                        .size();
        require_target(&(family.*field.member), 0.0f,
                       static_cast<float>(options - 1), true);
      }
    if constexpr (Resource::KIND == Pullback::ResourceKind::COLOR) {
      require_target(
          &family.palette_mapping, 0.0f,
          static_cast<float>(
              std::size(Pullback::Interp::Op::PALETTE_MAPPING_IDS) - 1),
          true);
      require_target(
          &family.palette_mode, 0.0f,
          static_cast<float>(std::size(Pullback::Color::PALETTE_MODE_IDS) - 1),
          true);
      require_target(
          &family.palette_type, 0.0f,
          static_cast<float>(std::size(AUTHORED_PALETTE_TYPE_IDS) - 1), true);
      require_target(
          &family.palette_harmony, 0.0f,
          static_cast<float>(std::size(AUTHORED_PALETTE_HARMONY_IDS) - 1),
          true);
      require_target(&family.palette_lightness_curve, 0.0f,
                     static_cast<float>(
                         std::size(AUTHORED_PALETTE_LIGHTNESS_CURVE_IDS) - 1),
                     true);
      require_target(
          &family.palette_name, 0.0f,
          static_cast<float>(std::size(AUTHORED_PALETTE_NAME_IDS) - 1), true);
      require_target(
          &family.palette_name_b, 0.0f,
          static_cast<float>(std::size(AUTHORED_PALETTE_NAME_IDS) - 1), true);
      require_target(
          &family.brightness_envelope, 0.0f,
          static_cast<float>(
              std::size(Pullback::Color::BRIGHTNESS_ENVELOPE_IDS) - 1),
          true);
    }
    if constexpr (Resource::KIND == Pullback::ResourceKind::MODULATION) {
      const auto &config = family.modulation;
      if constexpr (Family::DYNAMIC) {
        auto &live = const_cast<ColorModulation::Config &>(config);
        for (size_t index = 0;
             index < Pullback::Modulation::Controls::EXTRA_SCHEMA.size();
             ++index) {
          const auto &field =
              Pullback::Modulation::Controls::EXTRA_SCHEMA[index];
          require_target(
              Pullback::Modulation::Controls::address(live, index),
              field.enum_count ? 0.0f : field.min,
              field.enum_count ? float(field.enum_count - 1) : field.max,
              field.enum_count != 0,
              !Pullback::Modulation::Controls::active(config, index));
        }
        return;
      }

      const bool UNIT = ColorModulation::unit_mechanic(config.mechanic.kind);
      require_target(&config.mechanic.amount, UNIT ? 0.0f : -4.0f,
                     UNIT ? 1.0f : 4.0f);
      switch (config.source) {
      case Pattern::SourceKind::NOISE:
        require_target(&config.noise.scale, 1.0f / 64.0f,
                       config.evaluation == Pattern::Evaluation::BAKED_CUBE
                           ? 8.0f
                           : 64.0f);
        require_target(&config.noise.speed, -Pattern::MAX_NOISE_SPEED,
                       Pattern::MAX_NOISE_SPEED);
        break;
      case Pattern::SourceKind::DIRECTIONAL_SINE:
        require_target(&config.sine.frequency, 0.0f,
                       Pattern::MAX_AXIAL_FREQUENCY);
        require_target(&config.sine.speed, -Pattern::MAX_ANALYTIC_SPEED,
                       Pattern::MAX_ANALYTIC_SPEED);
        break;
      case Pattern::SourceKind::SPHERE_INTERFERENCE:
        require_target(&config.interference.frequency, 0.0f,
                       Pattern::MAX_AXIAL_FREQUENCY);
        require_target(&config.interference.speed, -Pattern::MAX_ANALYTIC_SPEED,
                       Pattern::MAX_ANALYTIC_SPEED);
        break;
      case Pattern::SourceKind::SPHERE_SMOOTH_GRID:
        require_target(&config.grid.frequency, 0.0f,
                       Pattern::MAX_AXIAL_FREQUENCY);
        require_target(&config.grid.speed, -Pattern::MAX_ANALYTIC_SPEED,
                       Pattern::MAX_ANALYTIC_SPEED);
        break;
      case Pattern::SourceKind::DISPLACEMENT:
        break;
      }
    }
  });
  HS_EXPECT_EQ(params.size(), expected);
}

/** @brief Direct modulation noise exposes scales beyond the baked-cube limit. */
inline void test_composed_noise_modulation_scale_limits() {
  for (const auto evaluation :
       {Pattern::Evaluation::DIRECT, Pattern::Evaluation::BAKED_CUBE}) {
    Pullback::Modulation::Config config;
    config.evaluation = evaluation;
    config.noise.scale = 20.0f;
    const auto fields = Pullback::Modulation::control_fields(config);
    HS_EXPECT_EQ(fields[1].min, 1.0f / 64.0f);
    HS_EXPECT_EQ(fields[1].max,
                 evaluation == Pattern::Evaluation::DIRECT ? 64.0f : 8.0f);
    HS_EXPECT_EQ(Pullback::Modulation::valid(config),
                 evaluation == Pattern::Evaluation::DIRECT);
  }

  reset_globals();
  AshCloud<SMALL_W, SMALL_H> effect;
  effect.init();
  const auto *scale = effect.getParameters().find("modulation.Pattern Scale");
  HS_EXPECT_TRUE(scale != nullptr);
  if (scale) {
    HS_EXPECT_EQ(scale->max, 64.0f);
    HS_EXPECT_EQ(effect.updateParameter("modulation.Pattern Scale", 20.0f),
                 ParamSetResult::APPLIED);
    HS_EXPECT_EQ(scale->get_requested(), 20.0f);
  }
}

/** @brief Requires exactly the listed registered sliders. */
inline void expect_exact_sliders(const ParamList &params,
                                 std::span<const char *const> expected,
                                 std::span<const char *const> terms = {}) {
  for (const auto names : {expected, terms})
    for (const char *slider : names) {
      HS_CONTEXT(slider);
      HS_EXPECT_TRUE(params.find(slider) != nullptr);
    }
  HS_EXPECT_EQ(params.size(), expected.size() + terms.size());
  HS_EXPECT_EQ(params.size(), params.capacity());
}

inline constexpr const char *POTENTIAL_TERM_SLIDERS[] = {
    "displace.term-1-type",      "displace.term-1-kind",
    "displace.term-1-l",         "displace.term-1-m",
    "displace.term-1-amplitude", "displace.term-1-phase",
    "displace.term-2-type",      "displace.term-2-kind",
    "displace.term-2-l",         "displace.term-2-m",
    "displace.term-2-amplitude", "displace.term-2-phase",
    "displace.term-3-type",      "displace.term-3-kind",
    "displace.term-3-l",         "displace.term-3-m",
    "displace.term-3-amplitude", "displace.term-3-phase",
    "displace.term-4-type",      "displace.term-4-kind",
    "displace.term-4-l",         "displace.term-4-m",
    "displace.term-4-amplitude", "displace.term-4-phase",
};

/** @brief Pins the registered controls of one specialization.
 * @tparam E Composed effect template.
 * @param name Effect name. */
template <template <int, int> class E>
inline void check_slider_registration(const char *name) {
  using FX = E<SMALL_W, SMALL_H>;
  using Params = typename FX::Params;
  using Traits = TraitsOf<FX>;
  HS_CONTEXT(name);

  reset_globals();
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

  if constexpr (requires { FX::Spec::AUTHORED_TOPOLOGY; })
    if constexpr (FX::Spec::AUTHORED_TOPOLOGY) {
      check_authored_slider_registration(effect);
      return;
    }

  if constexpr (std::is_same_v<FX, Iridescence<SMALL_W, SMALL_H>>) {
    constexpr const char *EXPECTED[] = {"Pattern Freq",
                                        "Speed",
                                        "Complexity",
                                        "Pattern Mix",
                                        "Singularity Fade",
                                        "Central Meridian",
                                        "displace.Strength",
                                        "displace.Orientation Azimuth",
                                        "displace.Orientation Polar",
                                        "displace.Drift",
                                        "displace.integrator",
                                        "Palette Chroma",
                                        "Palette Mapping",
                                        "Mapping Frequency",
                                        "Mapping Phase",
                                        "Phase Oscillation Depth",
                                        "Phase Oscillation Speed",
                                        "Brightness Bottom",
                                        "Brightness Top",
                                        "Opacity at Value 0",
                                        "Opacity at Value 1",
                                        "Palette Mode",
                                        "Palette Type",
                                        "Palette Harmony",
                                        "Palette Lightness Curve",
                                        "Palette Name",
                                        "Palette Name B",
                                        "Brightness Envelope",
                                        "Palette Hue",
                                        "Palette Spread",
                                        "Palette Lightness",
                                        "Palette Lightness Range",
                                        "Palette Mutation",
                                        "Gradient Start R",
                                        "Gradient Start G",
                                        "Gradient Start B",
                                        "Gradient Middle R",
                                        "Gradient Middle G",
                                        "Gradient Middle B",
                                        "Gradient End R",
                                        "Gradient End G",
                                        "Gradient End B",
                                        "modulate.Modulation Amount",
                                        "modulate.Pattern Scale",
                                        "modulate.Pattern Speed",
                                        "modulate2.Modulation Amount",
                                        "modulate2.Pattern Scale",
                                        "modulate2.Pattern Speed"};
    expect_exact_sliders(params, EXPECTED, POTENTIAL_TERM_SLIDERS);
    return;
  }

  if constexpr (std::is_same_v<FX, GridWarp<SMALL_W, SMALL_H>>) {
    constexpr const char *EXPECTED[] = {"Pattern Freq",
                                        "Speed",
                                        "Complexity",
                                        "Pattern Mix",
                                        "Singularity Fade",
                                        "Camera Wander",
                                        "Central Meridian",
                                        "displace.Strength",
                                        "displace.Orientation Azimuth",
                                        "displace.Orientation Polar",
                                        "displace.Drift",
                                        "displace.integrator",
                                        "Palette Chroma",
                                        "Palette Mapping",
                                        "Mapping Frequency",
                                        "Mapping Phase",
                                        "Phase Oscillation Depth",
                                        "Phase Oscillation Speed",
                                        "Opacity at Value 0",
                                        "Opacity at Value 1",
                                        "modulate.Modulation Amount",
                                        "modulate1.Modulation Amount",
                                        "modulate1.Pattern Scale",
                                        "modulate1.Pattern Speed"};
    expect_exact_sliders(params, EXPECTED, POTENTIAL_TERM_SLIDERS);
    return;
  }

  if constexpr (std::is_same_v<FX, LiquidWindows<SMALL_W, SMALL_H>>) {
    constexpr const char *EXPECTED[] = {"Pattern Freq",
                                        "Speed",
                                        "Complexity",
                                        "Pattern Mix",
                                        "Singularity Fade",
                                        "Central Meridian",
                                        "displace.Strength",
                                        "displace.Orientation Azimuth",
                                        "displace.Orientation Polar",
                                        "displace.Drift",
                                        "displace.integrator",
                                        "Palette Chroma",
                                        "Palette Mapping",
                                        "Mapping Frequency",
                                        "Mapping Phase",
                                        "Phase Oscillation Depth",
                                        "Phase Oscillation Speed",
                                        "Opacity at Value 0",
                                        "Opacity at Value 1",
                                        "modulate.Modulation Amount"};
    expect_exact_sliders(params, EXPECTED, POTENTIAL_TERM_SLIDERS);
    return;
  }

  if constexpr (std::is_same_v<FX, Millefiori<SMALL_W, SMALL_H>>) {
    constexpr const char *EXPECTED[] = {"Cell Scale",
                                        "Cell Radius",
                                        "Cell Shape",
                                        "Cell Softness",
                                        "Camera Wander",
                                        "Projection Spin Speed",
                                        "Projection Wander",
                                        "Singularity Fade",
                                        "Planar Warp 1 Speed",
                                        "Vortex Center X",
                                        "Vortex Center Y",
                                        "Vortex Radius",
                                        "Vortex Turns",
                                        "Vortex Center Orbit",
                                        "displace.Ripple Duration",
                                        "displace.Ripple Strength",
                                        "displace.Ripple Decay",
                                        "displace.Ripple Thickness",
                                        "displace.Ripple Center Azimuth",
                                        "displace.Ripple Center Polar",
                                        "displace.Burst Rate",
                                        "Palette Chroma",
                                        "Palette Mapping",
                                        "Mapping Frequency",
                                        "Mapping Phase",
                                        "Phase Oscillation Depth",
                                        "Phase Oscillation Speed",
                                        "Opacity at Value 0",
                                        "Opacity at Value 1",
                                        "Brightness Bottom",
                                        "Brightness Top",
                                        "field.Modulation Amount",
                                        "field.Pattern Frequency",
                                        "field.Pattern Speed",
                                        "modulate.Modulation Amount"};
    expect_exact_sliders(params, EXPECTED);
    return;
  }

  if constexpr (std::is_same_v<FX, AshCloud<SMALL_W, SMALL_H>>) {
    constexpr const char *EXPECTED[] = {"Cell Scale",
                                        "Cell Radius",
                                        "Cell Shape",
                                        "Cell Softness",
                                        "Camera Wander",
                                        "Projection Spin Speed",
                                        "Projection Wander",
                                        "Singularity Fade",
                                        "Surface Noise Scale",
                                        "Surface Noise Strength",
                                        "Surface Noise Speed",
                                        "displace.Ripple Duration",
                                        "displace.Ripple Strength",
                                        "displace.Ripple Decay",
                                        "displace.Ripple Thickness",
                                        "displace.Ripple Center Azimuth",
                                        "displace.Ripple Center Polar",
                                        "displace.Burst Rate",
                                        "Cutout Threshold",
                                        "Cutout Softness",
                                        "Palette Chroma",
                                        "Palette Mapping",
                                        "Mapping Frequency",
                                        "Mapping Phase",
                                        "Phase Oscillation Depth",
                                        "Phase Oscillation Speed",
                                        "Opacity at Value 0",
                                        "Opacity at Value 1",
                                        "modulation.Modulation Amount",
                                        "modulation.Pattern Scale",
                                        "modulation.Pattern Speed",
                                        "modulate.Modulation Amount"};
    expect_exact_sliders(params, EXPECTED);
    return;
  }

  verify_family_sliders<FX, typename Params::template Family<"source">>(params,
                                                                        true);
  verify_family_sliders<FX, Pullback::ProjectionParams>(params, true);
  verify_family_sliders<FX, typename Params::template Family<"surface">>(params,
                                                                         true);
  verify_family_sliders<FX, typename Params::template Family<"value">>(params,
                                                                       true);
  verify_warp_slot<FX, "outer_warp">(params, "Planar Warp 1 Speed");
  verify_warp_slot<FX, "inner_warp">(params, "Planar Warp 2 Speed");
  verify_warp_slot<FX, "affine">(params, "Planar Warp Speed");

  size_t twist_sliders = 0;
  typename FX::Params{}.visit([&]<typename Resource>(const auto &) {
    if constexpr (Resource::KIND == Pullback::ResourceKind::TWIST) {
      for (const auto &field : Resource::Family::FIELDS) {
        const std::string prefix =
            Pullback::ComposedDetail::qualified<
                Resource, typename Params::ResourceTypes>()
                ? std::string(Resource::KEY.view()) + "."
                : "";
        const auto *slider = params.find((prefix + field.name).c_str());
        HS_EXPECT_TRUE(slider != nullptr);
        if (slider) {
          HS_EXPECT_EQ(slider->min, field.min);
          HS_EXPECT_EQ(slider->max, field.max);
        }
        ++twist_sliders;
      }
    }
  });

  for (const char *slider : UNGATED_COLOR_SLIDERS)
    HS_EXPECT_TRUE(params.find(slider) != nullptr);
  constexpr bool mobius =
      std::is_same_v<typename Params::template Family<"lens">,
                     Pullback::MobiusGeneratorParams>;
  for (const char *slider : MOBIUS_SLIDERS)
    HS_EXPECT_EQ(params.find(slider) != nullptr, mobius);
  constexpr bool brightness =
      Traits::BRIGHTNESS != Pullback::Color::BrightnessEnvelope::NONE;
  constexpr bool displacement = Params::template HAS<"displacement">;
  constexpr bool spatial_modulation = Params::template HAS<"modulation">;
  HS_EXPECT_EQ(params.find("Modulation Amount") != nullptr,
               displacement || spatial_modulation);
  HS_EXPECT_EQ(params.find("Brightness Bottom") != nullptr, brightness);
  HS_EXPECT_EQ(params.find("Brightness Top") != nullptr, brightness);
  constexpr bool noise_modulation = [] {
    if constexpr (Params::template HAS<"modulation">)
      return Params::template Family<"modulation">::SOURCE ==
             Pattern::SourceKind::NOISE;
    return false;
  }();
  HS_EXPECT_EQ(params.find("Pattern Scale") != nullptr, noise_modulation);
  HS_EXPECT_EQ(params.find("Pattern Frequency") != nullptr,
               spatial_modulation && !noise_modulation);
  HS_EXPECT_EQ(params.find("Pattern Speed") != nullptr, spatial_modulation);
  if constexpr (Params::template HAS<"modulation">) {
    const auto initial = preset_params<FX>(0);
    const auto &config = initial.template get<"modulation">().modulation;
    using Fields = Pullback::Interp::Op::ModulationParams;
    const char *speed_id = nullptr;
    switch (config.source) {
    case Pattern::SourceKind::DISPLACEMENT:
      break;
    case Pattern::SourceKind::NOISE:
      speed_id = "pattern-speed";
      break;
    case Pattern::SourceKind::DIRECTIONAL_SINE:
      speed_id = "sine-speed";
      break;
    case Pattern::SourceKind::SPHERE_INTERFERENCE:
      speed_id = "interference-speed";
      break;
    case Pattern::SourceKind::SPHERE_SMOOTH_GRID:
      speed_id = "grid-speed";
      break;
    }
    const ParamDef *speed = params.find("Pattern Speed");
    HS_EXPECT_TRUE(speed != nullptr);
    if (speed != nullptr) {
      const float MAXIMUM = config.source == Pattern::SourceKind::NOISE
                                ? Pattern::MAX_NOISE_SPEED
                                : Pattern::MAX_ANALYTIC_SPEED;
      HS_EXPECT_EQ(speed->min, -MAXIMUM);
      HS_EXPECT_EQ(speed->max, MAXIMUM);
      for (const auto &field : Fields::EXTRA_SCHEMA)
        if (std::strcmp(field.id, speed_id) == 0) {
          HS_EXPECT_EQ(speed->min, field.min);
          HS_EXPECT_EQ(speed->max, field.max);
        }
    }
    const ParamDef *amount = params.find("Modulation Amount");
    HS_EXPECT_TRUE(amount != nullptr);
    if (amount != nullptr) {
      const bool UNIT = ColorModulation::unit_mechanic(config.mechanic.kind);
      HS_EXPECT_EQ(amount->min, UNIT ? 0.0f : -4.0f);
      HS_EXPECT_EQ(amount->max, UNIT ? 1.0f : 4.0f);
    }
  }

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

  const size_t source_topology = [] {
    using Family = typename Params::template Family<"source">;
    if constexpr (requires { Family::TOPOLOGY; })
      return std::size(Family::TOPOLOGY);
    return size_t{0};
  }();
  const size_t expected =
      source_topology + twist_sliders +
      named_field_count<FX, typename Params::template Family<"source">>() +
      named_field_count<FX, Pullback::ProjectionParams>() +
      named_field_count<FX, typename Params::template Family<"surface">>() +
      named_field_count<FX, typename Params::template Family<"value">>() +
      warp_slot_count<FX, "outer_warp">() +
      warp_slot_count<FX, "inner_warp">() + warp_slot_count<FX, "affine">() +
      (mobius ? std::size(MOBIUS_SLIDERS) : size_t{0}) +
      std::size(UNGATED_COLOR_SLIDERS) + (brightness ? size_t{2} : size_t{0}) +
      (spatial_modulation ? size_t{3} : size_t{0}) +
      (displacement ? size_t{1} : size_t{0});
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
  HS_CONTEXT(name);

  reset_globals();
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
  moved.params.visit([&]<typename Resource>(auto &family) {
    if constexpr (Pullback::HasFields<typename Resource::Family>)
      move_off_captured(family, captured.params.template get<Resource::KEY>());
    if constexpr (requires { family.terms; })
      family.terms[0].amplitude =
          family.terms[0].amplitude == 0.125f ? 0.25f : 0.125f;
    if constexpr (std::is_same_v<typename Resource::Family,
                                 Pullback::Color::AuthoredColorParams>)
      family.palette_mode = (family.palette_mode + 1) % 3;
  });
  constexpr size_t MAPPINGS =
      std::size(Pullback::Interp::Op::PALETTE_MAPPING_IDS);
  moved.params.template get<"color">().palette_mapping =
      static_cast<Pullback::Color::PaletteMapping>(
          (static_cast<uint8_t>(
               captured.params.template get<"color">().palette_mapping) +
           1) %
          MAPPINGS);
  moved.params.visit([&]<typename Resource>(auto &family) {
    if constexpr (Resource::KIND == Pullback::ResourceKind::MODULATION) {
      auto &config = family.modulation;
      config.mechanic.amount = config.mechanic.amount == 0.0f ? 0.5f : 0.0f;
      if (config.source == Pattern::SourceKind::NOISE) {
        config.noise.scale = config.noise.scale == 1.0f ? 2.0f : 1.0f;
        config.noise.speed = config.noise.speed == 0.0f ? -0.0001f : 0.0f;
      }
    }
  });
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
    if constexpr (Resource::KIND == Pullback::ResourceKind::MODULATION) {
      typename FX::ParameterSnapshot invalid = captured;
      invalid.params.template get<Resource::KEY>().modulation.mechanic.amount =
          std::numeric_limits<float>::quiet_NaN();
      HS_EXPECT_FALSE(effect.restore_parameters(invalid));
      verify_params_equal(effect.serialize_parameters().params,
                          captured.params);
    }
  });

  typename FX::ParameterSnapshot mapping = captured;
  mapping.params.template get<"color">().palette_mapping =
      static_cast<Pullback::Color::PaletteMapping>(MAPPINGS);
  HS_EXPECT_FALSE(effect.restore_parameters(mapping));
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

  reset_globals();
  FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), FX::PRESET_IDS.size());
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});

  for (size_t index = 0; index < FX::PRESET_IDS.size(); ++index) {
    HS_CONTEXT("preset", static_cast<int>(index));
    HS_EXPECT_FALSE(FX::PRESET_IDS[index].empty());
    for (size_t earlier = 0; earlier < index; ++earlier)
      HS_EXPECT_TRUE(FX::PRESET_IDS[index] != FX::PRESET_IDS[earlier]);
    HS_EXPECT_TRUE(FX::valid_params(preset_params<FX>(index)));

    HS_EXPECT_TRUE(effect.selectPreset(index));
    HS_EXPECT_EQ(effect.getPresetIndex(), index);
    verify_params_equal(effect.serialize_parameters().params,
                        preset_params<FX>(index));
  }
  HS_EXPECT_FALSE(effect.selectPreset(FX::PRESET_IDS.size()));
}

/**
 * @brief Pins the parameter interpolation a preset crossfade runs on.
 * @tparam E Composed effect class template.
 * @param name Effect name, for the failure context.
 * @details Both endpoints must come back exactly, except that a modulation
 * topology change takes the target at progress 0, and no sample in between may
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
      const Params a = preset_params<FX>(from);
      Params b = preset_params<FX>(to);
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
      // A modulation recipe takes the target topology at once.
      Params start = a;
      start.visit([&]<typename Resource>(auto &family) {
        if constexpr (requires { family.modulation; }) {
          const auto &target = b.template get<Resource::KEY>().modulation;
          if (!Pullback::Modulation::same_topology(family.modulation, target))
            family.modulation = target;
        }
      });
      verify_params_equal(Pullback::interpolate(a, b, 0.0f), start);
      verify_params_equal(Pullback::interpolate(a, b, 1.0f), b);
    }
  }
}

#include "tests/composed_effect/derivation_reach.h"
#include "tests/composed_effect/document_values.h"
#include "tests/composed_effect/roster_contracts.h"
#include "tests/composed_effect/stage_probes.h"
#include "tests/composed_effect/authored_palette.h"
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

/**
 * @brief AlienBrain holds its first preset for PRESET_DWELL_FRAMES, then
 * advances to preset 1.
 */
inline void test_alien_brain_preset_dwell() {
  using FX = AlienBrain<SMALL_W, SMALL_H>;
  reset_globals();
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
  reset_globals();
  FX effect;
  effect.init();
  const auto authored = effect.serialize_parameters().params;
  effect.draw_frame();
  effect.advance_display();
  const auto first = effect.frame_for_test().resources.template get<"lens">();
  HS_EXPECT_TRUE(std::isfinite(first.magnitude));
  HS_EXPECT_GT(first.magnitude, 0.0f);
  effect.draw_frame();
  effect.advance_display();
  const auto second = effect.frame_for_test().resources.template get<"lens">();
  HS_EXPECT_EQ(second.magnitude, first.magnitude);
  HS_EXPECT_NE(first.angle, second.angle);
  verify_params_equal(effect.serialize_parameters().params, authored);
  effect.setAnimationsPaused(true);
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_EQ(effect.frame_for_test().resources.template get<"lens">().angle,
               second.angle);
  effect.setAnimationsPaused(false);
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_NE(effect.frame_for_test().resources.template get<"lens">().angle,
               second.angle);
  verify_params_equal(effect.serialize_parameters().params, authored);
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

#include "tests/composed_effect/temporal_sampling.h"
#include "tests/composed_effect/retained_hue.h"
#include "tests/composed_effect/chamber_noise_cache.h"
#include "tests/composed_effect/dodecahedral_matrix_cache.h"
#include "tests/composed_effect/ashcloud_shared_color.h"
#include "tests/composed_effect/palette_rgb_cache.h"

/** @brief Every composed effect defaults to blanking displayed columns. */
inline void test_composed_strobe_roster() {
#define HS_COMPOSED_STROBE(name)                                               \
  {                                                                            \
    HS_CONTEXT(#name);                                                         \
    reset_globals();                                                           \
    name<SMALL_W, SMALL_H> effect;                                             \
    HS_EXPECT_TRUE(effect.strobe_columns());                                   \
  }
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_STROBE)
#undef HS_COMPOSED_STROBE
}

#include "tests/composed_effect/grid_worms.h"

/** @brief Millefiori pixels match its document-built operator chain. */
inline void test_millefiori_document_pipeline() {
  composed_chain_tests::verify_export<Millefiori<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
}

/** @brief AshCloud pixels match its document-built operator chain. */
inline void test_ash_cloud_document_pipeline() {
  using FX = AshCloud<SMALL_W, SMALL_H>;
  using Ripple = Pullback::Stage::Displace<
      Pullback::Surface::BurstRipple<Pullback::SurfaceProvider<
          FX::Binding, Pullback::SphereAnimation::BurstRippleParams, true,
          "displace">>>;
  using Curl = Pullback::Stage::Displace<Pullback::Surface::CurlNoise<
      Pullback::SurfaceProvider<FX::Binding, Pullback::SurfaceNoiseParams,
                                true>,
      math::NoiseBasis::SIMPLEX, AshCloudInlineEuler,
      math::TangentLimit::SMOOTH>>;
  using Lens =
      Pullback::Stage::Lens<Pullback::Lens::DodecahedralTableKaleidoscope>;
  using Project = Pullback::Stage::Project<Pullback::ProjectionPolicyFor<
      Pullback::ProjectionKind::STEREOGRAPHIC, FX::Binding>::Type>;
  using Leaves = Pullback::ComposedDetail::PipelineLeaves<
      FX::Spec::template Pipeline<FX::Binding>>;
  static_assert(std::is_same_v<Leaves::StageAt<1>, Ripple>);
  static_assert(std::is_same_v<Leaves::StageAt<2>, Curl>);
  static_assert(std::is_same_v<Leaves::StageAt<3>, Lens>);
  static_assert(std::is_same_v<Leaves::StageAt<4>, Project>);
  struct Frame {
    struct Temporal {
      float fraction;
    } temporal;
  };
  for (const auto &view : {math::X_AXIS, -math::Y_AXIS, math::Z_AXIS,
                           math::Vector(-0.3f, 0.4f, -0.5f).normalized()}) {
    const auto EXPECTED = lenses::dodecahedral_kaleidoscope_lens(view);
    for (float fraction : {0.0f, 1.0f}) {
      const auto ACTUAL = Pullback::Lens::DodecahedralTableKaleidoscope::apply(
          view, Frame{{fraction}});
      HS_EXPECT_EQ(ACTUAL.x, EXPECTED.x);
      HS_EXPECT_EQ(ACTUAL.y, EXPECTED.y);
      HS_EXPECT_EQ(ACTUAL.z, EXPECTED.z);
    }
  }
  composed_chain_tests::verify_export<FX>(0, DOCUMENT_CHANNEL_TOLERANCE);
}

#include "tests/composed_effect/authored_policies.h"
/** @brief Promoted pixels match the authored document-built chains. */
inline void test_promoted_document_pipelines() {
  static_assert(Involution<SMALL_W, SMALL_H>::authored_instance("source") ==
                "sample");
  composed_chain_tests::verify_export<
      KaleidoscopeStainedGlass<SMALL_W, SMALL_H>>(0,
                                                  DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<KaleidoscopeHexBright<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<MobiusGrid<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<AlienCore<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<Involution<SMALL_W, SMALL_H>>(0);
  composed_chain_tests::verify_export<GridWarp<SMALL_W, SMALL_H>>(0);
  composed_chain_tests::verify_export<LiquidWindows<SMALL_W, SMALL_H>>(0);
  composed_chain_tests::verify_export<Iridescence<SMALL_W, SMALL_H>>(0);
  composed_chain_tests::verify_export<CosmicEyeball<SMALL_W, SMALL_H>>(0);
  composed_chain_tests::verify_export<KaleidoscopeHexOil<SMALL_W, SMALL_H>>(0);
  composed_chain_tests::verify_export<AlienBrain<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<AlienOcean<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<KaleidoscopeMandala<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<GridSpace<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<ChromaticLichen<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<MermaidSkin<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<JewelMelt<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<KaleidoscopeFlowers<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
}

/** @brief Module entry point for the composed-effect base contract. */
inline int run_composed_effect_tests() {
  ModuleFixture fixture("composed_effect");
  HS_RUN(test_authored_tangent_bindings_match_operators);
  HS_RUN(test_authored_tangent_fields_after_twist);
  HS_RUN(test_authored_tangent_interval_uses_radians);
  HS_RUN(test_authored_potential_interpolation_and_admission);
  HS_RUN(test_authored_interpolation_destination_aliases);
  HS_RUN(test_authored_surface_noise_topology);
  HS_RUN(test_authored_dynamic_modulation_storage);
  HS_RUN(test_authored_projection_bindings_match_operators);
  HS_RUN(test_authored_dynamic_modulation_controls);
  HS_RUN(test_authored_potential_nested_controls);
  HS_RUN(test_authored_noise_seed_identity);
  HS_RUN(test_millefiori_radial_table_interval);
  HS_RUN(test_millefiori_document_pipeline);
  HS_RUN(test_promoted_document_pipelines);
  HS_RUN(test_grid_worms_twist_clocks);
  HS_RUN(test_grid_worms_document_pipeline);
  HS_RUN(test_ash_cloud_document_pipeline);
  HS_RUN(test_ash_cloud_shared_color_identity);
  HS_RUN(test_ash_cloud_shared_column_qualification);
  HS_RUN(test_ash_cloud_shared_column_fallbacks);
  HS_RUN(test_composed_default_gradient_preparation);
  HS_RUN(test_composed_spec_default_color_paths);
  HS_RUN(test_composed_persistent_prepared_frame);
  HS_RUN(test_composed_document_modulation_members);
  HS_RUN(test_millefiori_column_color_phase);
  HS_RUN(test_grid_worms_sine_column_phase);
  HS_RUN(test_retained_equal_palette_hue_view);
  HS_RUN(test_retained_palette_identity_rejects_changes);
  HS_RUN(test_ash_cloud_noise_column_phase);
  HS_RUN(test_baked_cube_column_phase);
  HS_RUN(test_chamber_noise_cache_bake_lifetime);
  HS_RUN(test_chamber_noise_cache_admission);
  HS_RUN(test_chamber_noise_cache_phase_and_fallback);
  HS_RUN(test_chamber_noise_cache_endpoint_admission);
  HS_RUN(test_chamber_noise_cache_equal_clock_fraction);
  HS_RUN(test_chamber_noise_cache_far_endpoint_falls_back);
  HS_RUN(test_chamber_noise_cache_rotated_endpoint_falls_back);
  HS_RUN(test_chamber_noise_cache_source_phase_invalidates_patch);
  HS_RUN(test_chamber_noise_cache_runtime_history);
  HS_RUN(test_dodecahedral_matrix_cache_lifetime);
  HS_RUN(test_dodecahedral_matrix_cache_fold_parity);
  HS_RUN(test_dodecahedral_matrix_cache_lens_binding);
  HS_RUN(test_dodecahedral_face_cache_lifetime);
  HS_RUN(test_dodecahedral_face_cache_fold_parity);
  HS_RUN(test_dodecahedral_face_cache_lens_binding);
  HS_RUN(test_palette_rgb_cache_axes_and_wrap);
  HS_RUN(test_palette_rgb_cache_guard_mask_aliasing);
  HS_RUN(test_palette_lab_reference_endpoint_lifetime);
  HS_RUN(test_palette_rgb_cache_matches_hue_rotation_table);
  HS_RUN(test_palette_rgb_cache_lazy_guards);
  HS_RUN(test_palette_rgb_cache_runtime_lifetime);
  HS_RUN(test_palette_rgb_cache_bake_reuse);
  HS_RUN(test_palette_rgb_cache_shared_base_endpoints);
  HS_RUN(test_palette_rgb_cache_midpoint_reuse);
  HS_RUN(test_palette_rgb_cache_any_chroma);
  HS_RUN(test_palette_rgb_cache_q16_compatibility);
  HS_RUN(test_palette_rgb_midpoint_only_storage);
  HS_RUN(test_palette_rgb_midpoint_only_fallback);
  HS_RUN(test_palette_lab_cache_color_padding);
  HS_RUN(test_palette_lab_cache_column_controls);
  HS_RUN(test_palette_lab_cache_through_lerp);
  HS_RUN(test_palette_lab_cache_palette_control_step);
  HS_RUN(test_millefiori_palette_lab_lifetime);
  HS_RUN(test_millefiori_ripple_endpoints);
  HS_RUN(test_composed_strobe_roster);
  HS_RUN(test_composed_analytic_modulation_sliders);
  HS_RUN(test_composed_authored_palette_sources);
  HS_RUN(test_composed_authored_palette_harmony_clock);
  HS_RUN(test_composed_authored_palette_admission);
  HS_RUN(test_composed_analytic_modulation_admission);
  HS_RUN(test_composed_path_modulation_admission);
  HS_RUN(test_temporal_clock_sampling);
  HS_RUN(test_temporal_affine_cycle);
  HS_RUN(test_temporal_custom_shade_is_authoritative);
  HS_RUN(test_scoped_temporal_sampling_restores_entering_mode);
  HS_RUN(test_gridspace_column_basis_geometry);
  HS_RUN(test_gridspace_canonical_columns);
  HS_RUN(test_temporal_equal_palettes_color_as_current);
  HS_RUN(test_temporal_moving_families_interpolate);
  HS_RUN(test_composed_modulation_blends_axes);
  HS_RUN(test_temporal_owned_resources);
  HS_RUN(test_temporal_topology_interval_selects_current_endpoint);
  HS_RUN(test_modulation_interpolate_retains_target_topology);
  HS_RUN(test_composed_factored_columns_match_pixels);
  HS_RUN(test_temporal_modulation_admission);
  HS_RUN(test_temporal_modulation_bake_generation);
  HS_RUN(test_temporal_stacked_modulation_resources);
  HS_RUN(test_gridspace_temporal_pixels);
  HS_RUN(test_composed_temporal_roster);
  HS_RUN(test_mobius_grid_instantaneous_pixel_golden);
  HS_RUN(test_mobius_grid_lens_follows_phase);
  HS_RUN(test_mobius_interval_exponential);
  HS_RUN(test_mobius_grid_segmented_matches_single_engine);
  HS_RUN(test_alien_brain_preset_dwell);
  HS_RUN(test_mobius_grid_circular_animation);
  HS_RUN(test_catalog_semantic_export);
  HS_RUN(test_modulation_schema_defaults_match_controls);
  HS_RUN(test_composed_hand_registered_families);
  HS_RUN(test_composed_direct_surface_placement);
  HS_RUN(test_composed_slider_registration);
  HS_RUN(test_composed_modulation_is_device_qualified);
  HS_RUN(test_composed_noise_modulation_scale_limits);
  HS_RUN(check_slider_registration<NoHueProbe>, "NoHueProbe");
  HS_RUN(test_mobius_frame_admission);
  HS_RUN(test_mobius_generator_interpolation);
  HS_RUN(test_mobius_captured_parameter_restore);
  HS_RUN(test_mobius_automatic_departures_interpolate_generator);
  HS_RUN(test_composed_snapshot_contract);
  HS_RUN(test_parameter_layout_reorder);
  HS_RUN(test_composed_parameter_schema_pins);
  HS_RUN(test_composed_preset_choreography);
  HS_RUN(test_composed_log_positive_curve);
  HS_RUN(test_composed_preset_interpolation);
  HS_RUN(test_palette_mapping_crossfade);
  HS_RUN(test_composed_document_values);
  HS_RUN(test_document_numeric_resource_selection);
  HS_RUN(test_flowers_longitude_seam);
  HS_RUN(test_composed_derivation_reach);
  HS_RUN(test_composed_derivation_values_are_emitted);
  HS_RUN(test_composed_keyed_parameters);
  HS_RUN(test_composed_keyed_value_providers);
  HS_RUN(test_composed_affine_named_source);
  HS_RUN(test_composed_affine_form_sliders);
  HS_RUN(test_composed_repeated_instances);
  HS_RUN(test_composed_projection_frame_identity);
  HS_RUN(test_composed_projection_walk_storage);
  HS_RUN(test_composed_periodic_ripple_surface);
  HS_RUN(test_composed_noise_source);
  HS_RUN(hs_test::catalog_policy_tests::test_palette_selection_schema_parity);
  HS_RUN(hs_test::catalog_policy_tests::test_authored_warp_catalog_parity);
  HS_RUN(hs_test::catalog_policy_tests::test_authored_source_catalog_parity);
  HS_RUN(hs_test::catalog_policy_tests::test_authored_transfer_catalog_parity);
  HS_RUN(hs_test::catalog_policy_tests::test_authored_sample_crossing);
  HS_RUN(
      hs_test::catalog_policy_tests::test_authored_twin_wave_sample_admission);
  HS_RUN(hs_test::catalog_policy_tests::test_edge_fade_squared_curve);
  HS_RUN(hs_test::keyed_rotation_tests::test_keyed_rotation_runtime_parity);
  HS_RUN(test_choreography_fade_departure);
  HS_RUN(test_choreography_fade_resumes_after_pause);
  HS_RUN(test_choreography_lerp_transition_hooks);
  HS_RUN(test_choreography_lerp_pause_policy);
  return fixture.result();
}

} // namespace composed_effect_tests
} // namespace hs_test
