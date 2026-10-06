/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Host unit tests for the WASM parameter-marshaling layer
 * (targets/wasm/param_marshal.h).
 */
#pragma once

#include "targets/effects.h" // HS_EFFECT_LIST roster
#include "core/render/canvas.h"
#include "core/memory.h"
#include "targets/wasm/param_marshal.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <vector>

namespace hs_test {
namespace param_marshal_tests {

constexpr int DEFAULT_W = 288;
constexpr int DEFAULT_H = 144;

/**
 * @brief Tracks whether the roster can distinguish a transposed ParamView pair.
 * @details ParamView is positionally initialized; a transposed min/max or
 *   animated/readonly pair is only detectable when some param has distinct
 *   values.
 */
struct FieldCoverage {
  bool range_distinguishing = false; /**< Some param has min != max. */
  bool flags_distinguishing =
      false;                       /**< Some param has animated != readonly. */
  bool float_enum_present = false; /**< Some enum has a float target. */
};

enum class RoundTripResult { COVERED, NO_EDITABLE_FLOAT, STREAM_MISMATCH };

/**
 * @brief Marshals one effect through the WASM bridge and asserts the definition
 *        and value streams stay consistent with the source params.
 * @tparam E Effect template, instantiated at the test canvas size DEFAULT_W x DEFAULT_H.
 * @param name Effect name for assertion context.
 * @param coverage In/out tally of whether the roster supplies a param that can
 *        distinguish a transposed field pair.
 * @details Checks equal length, index-aligned name/value/type/range/flags, and
 *          a write-by-name that round-trips to the same index without disturbing
 *          order.
 */
template <template <int, int> class E>
inline RoundTripResult check_one(const char *name, FieldCoverage &coverage) {
  HS_CONTEXT(name);
  reset_globals();

  E<DEFAULT_W, DEFAULT_H> effect;
  effect.init();

  std::vector<hs_wasm::ParamView> views;
  std::vector<float> values;
  hs_wasm::collect_param_views(effect, views);
  hs_wasm::fill_param_values(effect, values);

  // Count the source params independently of the marshaled streams.
  size_t n = 0;
  for (const auto &def : effect.getParameters()) {
    (void)def;
    ++n;
  }
  HS_EXPECT_EQ(views.size(), n);
  HS_EXPECT_EQ(values.size(), n);
  if (views.size() != n || values.size() != n)
    return RoundTripResult::STREAM_MISMATCH;

  // For every i, name/value/type from the independent passes must match the
  // source param.
  size_t i = 0;
  for (const auto &def : effect.getParameters()) {
    HS_EXPECT(i < views.size(), "marshaled index within range");
    HS_EXPECT_EQ(std::string_view(views[i].name), std::string_view(def.name));
    HS_EXPECT_EQ(views[i].value, values[i]);
    HS_EXPECT_EQ(views[i].is_bool, def.is_bool());
    // Every whole-numbered param, enums included, carries the integer flag.
    HS_EXPECT_EQ(views[i].is_integer, def.is_integer() || def.is_enum());
    coverage.float_enum_present |= def.is_enum() && !def.is_integer();
    HS_EXPECT_EQ(views[i].value, def.get());
    HS_EXPECT_EQ(views[i].requested_value, def.get_requested());
    HS_EXPECT_EQ(
        views[i].accepted_value,
        static_cast<const Effect &>(effect).accepted_parameter_value(def));
    HS_EXPECT_EQ(views[i].min, def.min);
    HS_EXPECT_EQ(views[i].max, def.max);
    HS_EXPECT_EQ(views[i].animated, def.animated);
    HS_EXPECT_EQ(views[i].readonly, def.readonly);
    HS_EXPECT_EQ(views[i].preset, def.preset);
    coverage.range_distinguishing |= def.min != def.max;
    coverage.flags_distinguishing |= def.animated != def.readonly;
    HS_EXPECT_TRUE(views[i].options == def.options);
    HS_EXPECT_EQ(views[i].option_count, def.option_count);
    HS_EXPECT_TRUE(views[i].export_options == def.export_options);
    HS_EXPECT_TRUE(views[i].option_values == def.option_values);
    // An enum's current value is always a valid option index.
    if (def.option_count > 0) {
      HS_EXPECT_GE(views[i].value, 0.0f);
      if (def.option_values == nullptr)
        HS_EXPECT_LE(views[i].value, static_cast<float>(def.option_count - 1));
      else {
        bool listed = false;
        for (int k = 0; k < def.option_count; ++k)
          listed |= views[i].value == static_cast<float>(def.option_values[k]);
        HS_EXPECT_TRUE(listed);
      }
    }
    ++i;
  }

  int target = -1;
  for (size_t k = 0; k < views.size(); ++k) {
    const auto &v = views[k];
    if (!v.is_bool && !v.is_integer && !v.readonly && !v.animated &&
        v.max > v.min) {
      target = static_cast<int>(k);
      break;
    }
  }
  if (target < 0)
    return RoundTripResult::NO_EDITABLE_FLOAT;

  const float lo = views[target].min, hi = views[target].max;
  float newv = lo + 0.5f * (hi - lo);
  if (newv == views[target].value)
    newv = lo + 0.25f * (hi - lo);

  HS_EXPECT(effect.updateParameter(views[target].name, newv) ==
                ParamSetResult::APPLIED,
            "updateParameter by name succeeds for an editable param");

  std::vector<hs_wasm::ParamView> views2;
  hs_wasm::collect_param_views(effect, views2);
  HS_EXPECT_EQ(views2.size(), views.size());
  if (views2.size() != views.size())
    return RoundTripResult::STREAM_MISMATCH;
  HS_EXPECT_NEAR(views2[target].value, newv, 1e-3f);
  hs_wasm::fill_param_values(effect, values);
  HS_EXPECT_EQ(values.size(), views2.size());
  if (values.size() != views2.size())
    return RoundTripResult::STREAM_MISMATCH;
  HS_EXPECT_NEAR(values[target], newv, 1e-3f);
  for (size_t k = 0; k < views2.size(); ++k)
    HS_EXPECT_EQ(std::string_view(views2[k].name),
                 std::string_view(views[k].name));
  return RoundTripResult::COVERED;
}

/**
 * @brief Tallies an effect's parameter count, tracking the roster maximum.
 * @tparam E Effect template, instantiated at the test canvas size DEFAULT_W x DEFAULT_H.
 * @param max_count In/out roster maximum.
 */
template <template <int, int> class E>
inline void count_one(size_t &max_count) {
  reset_globals();

  E<DEFAULT_W, DEFAULT_H> effect;
  effect.init();

  size_t n = 0;
  for (const auto &def : effect.getParameters()) {
    (void)def;
    ++n;
  }
  if (n > max_count)
    max_count = n;
}

/**
 * @brief Verifies the per-frame memory-view stability contract: refilling the
 *        reserved stream vectors never reallocates their backing storage.
 * @tparam E Effect template, instantiated at the test canvas size DEFAULT_W x DEFAULT_H.
 * @param name Effect name for assertion context.
 * @param views Reusable definition-stream vector, pre-reserved by the caller.
 * @param values Reusable value-stream vector, pre-reserved by the caller.
 * @param view_data Expected backing pointer of @p views (its .data() before
 *        refill).
 * @param view_cap Expected capacity of @p views.
 * @param value_data Expected backing pointer of @p values (its .data() before
 *        refill).
 * @param value_cap Expected capacity of @p values.
 * @details Uses the engine's ParamStreams storage across effect switches.
 */
template <template <int, int> class E>
inline void
check_stability_one(const char *name, std::vector<hs_wasm::ParamView> &views,
                    std::vector<float> &values,
                    const hs_wasm::ParamView *view_data, size_t view_cap,
                    const float *value_data, size_t value_cap) {
  HS_CONTEXT(name);
  reset_globals();

  E<DEFAULT_W, DEFAULT_H> effect;
  effect.init();

  hs_wasm::collect_param_views(effect, views);
  hs_wasm::fill_param_values(effect, values);

  HS_EXPECT(views.data() == view_data,
            "collect_param_views reused the reserved buffer (no realloc)");
  HS_EXPECT(values.data() == value_data,
            "fill_param_values reused the reserved buffer (no realloc)");
  HS_EXPECT_EQ(views.capacity(), view_cap);
  HS_EXPECT_EQ(values.capacity(), value_cap);
}

/**
 * @brief Freezes the effect roster ORDER, not just its count.
 * @details GOLDEN_ROSTER pins the HS_EFFECT_LIST order that sets the effect
 *   ordinal; update it deliberately on any reorder, insertion or removal.
 */
inline void check_roster_order_pinned() {
  // Independent hand-maintained copy of the intended roster order. Must NOT be
  // generated from HS_EFFECT_LIST, or the comparison becomes a tautology.
  static const char *const GOLDEN_ROSTER[] = {"BZReactionDiffusion",
                                              "Fishbowl",
                                              "Comets",
                                              "GridSpace",
                                              "HyperLattice",
                                              "AshCloud",
                                              "LatticeMelt",
                                              "ChromaticLichen",
                                              "MermaidSkin",
                                              "DisplacementField",
                                              "DreamBalls",
                                              "Dynamo",
                                              "KaleidoscopeFlowers",
                                              "KaleidoscopeSmooth",
                                              "KaleidoscopeMandala",
                                              "GnomonicStars",
                                              "GSReactionDiffusion",
                                              "AlienCore",
                                              "HankinSolids",
                                              "HopfFibration",
                                              "IslamicStars",
                                              "KaleidoscopeHexBright",
                                              "AlienOcean",
                                              "KaleidoscopeHexSoft",
                                              "MeshFeedback",
                                              "MindSplatter",
                                              "MobiusGrid",
                                              "MobiusRings",
                                              "PetalFlow",
                                              "KaleidoscopePentBright",
                                              "KaleidoscopeHexOil",
                                              "Raymarch",
                                              "RingShower",
                                              "RingSpin",
#if HS_ENABLE_CHAIN_INTERPRETER
                                              "ShaderChain",
#endif
                                              "ShapeShifter",
                                              "AlienBrain",
                                              "SphericalHarmonics",
                                              "CosmicEyeball",
                                              "Thrusters",
                                              "KaleidoscopeStainedGlass",
                                              "Voronoi"};
  // Actual roster, expanded straight from the X-macro source of truth.
  static const char *const ACTUAL_ROSTER[] = {
#define HS_EFFECT_NAME(name) #name,
      HS_EFFECT_LIST(HS_EFFECT_NAME)
#undef HS_EFFECT_NAME
  };
  constexpr size_t GOLDEN_N = sizeof(GOLDEN_ROSTER) / sizeof(GOLDEN_ROSTER[0]);
  constexpr size_t ACTUAL_N = sizeof(ACTUAL_ROSTER) / sizeof(ACTUAL_ROSTER[0]);
  HS_EXPECT_EQ(ACTUAL_N, GOLDEN_N);
  HS_EXPECT_EQ(static_cast<int>(ACTUAL_N), HS_EFFECT_COUNT);
  const size_t n = ACTUAL_N < GOLDEN_N ? ACTUAL_N : GOLDEN_N;
  for (size_t i = 0; i < n; ++i) {
    HS_CONTEXT("roster index", i);
    HS_EXPECT_EQ(std::string_view(ACTUAL_ROSTER[i]),
                 std::string_view(GOLDEN_ROSTER[i]));
  }
}

/** @brief The WASM token changes for replacement and local schema rebinds. */
inline void check_generation_tracker() {
  hs_wasm::ParamGenerationTracker tracker;
  HS_EXPECT_EQ(tracker.generation(), uint32_t(0));
  tracker.replace(7);
  HS_EXPECT_EQ(tracker.generation(), uint32_t(1));
  tracker.observe(7);
  HS_EXPECT_EQ(tracker.generation(), uint32_t(1));
  tracker.observe(8);
  HS_EXPECT_EQ(tracker.generation(), uint32_t(2));
  tracker.replace(8);
  HS_EXPECT_EQ(tracker.generation(), uint32_t(3));
  tracker.replace(0);
  HS_EXPECT_EQ(tracker.generation(), uint32_t(4));
}

/** @brief Pins HyperLattice pattern and view dropdown metadata. */
inline void check_hyper_lattice_pattern_view_dropdowns() {
  reset_globals();
  HyperLattice<DEFAULT_W, DEFAULT_H> effect;
  effect.init();
  std::vector<hs_wasm::ParamView> views;
  hs_wasm::collect_param_views(effect, views);
  const hs_wasm::ParamView *view = nullptr;
  const hs_wasm::ParamView *pattern = nullptr;
  for (const hs_wasm::ParamView &entry : views) {
    if (std::string_view(entry.name) == "View")
      view = &entry;
    if (std::string_view(entry.name) == "Pattern")
      pattern = &entry;
  }
  HS_EXPECT_TRUE(view != nullptr && pattern != nullptr);
  if (view == nullptr || pattern == nullptr)
    return;
  HS_EXPECT_FALSE(view->is_bool);
  HS_EXPECT_TRUE(view->is_integer);
  HS_EXPECT_EQ(view->option_count, 2);
  HS_EXPECT_EQ(pattern->option_count, 3);
  HS_EXPECT_TRUE(view->options != nullptr);
  HS_EXPECT_TRUE(view->export_options != nullptr);
  HS_EXPECT_TRUE(pattern->options != nullptr);
  HS_EXPECT_TRUE(pattern->export_options != nullptr);
  HS_EXPECT_TRUE(pattern->option_values != nullptr);
  if (view->option_count != 2 || pattern->option_count != 3 || !view->options ||
      !view->export_options || !pattern->options || !pattern->export_options ||
      !pattern->option_values)
    return;
  HS_EXPECT_EQ(std::string_view(view->options[0]),
               std::string_view("3D perspective"));
  HS_EXPECT_EQ(std::string_view(view->options[1]),
               std::string_view("4D slice"));
  HS_EXPECT_EQ(std::string_view(view->export_options[1]),
               std::string_view("LatticeMode::FOUR_D_SLICE"));
  HS_EXPECT_FALSE(pattern->is_bool);
  HS_EXPECT_TRUE(pattern->is_integer);
  HS_EXPECT_EQ(std::string_view(pattern->options[0]),
               std::string_view("Cubic"));
  HS_EXPECT_EQ(std::string_view(pattern->export_options[0]),
               std::string_view("Pattern::CUBIC_WIRE"));
  HS_EXPECT_EQ(std::string_view(pattern->options[1]),
               std::string_view("Octet Truss"));
  HS_EXPECT_EQ(std::string_view(pattern->export_options[1]),
               std::string_view("Pattern::OCTET"));
  HS_EXPECT_EQ(std::string_view(pattern->options[2]),
               std::string_view("Shells"));
  HS_EXPECT_EQ(std::string_view(pattern->export_options[2]),
               std::string_view("Pattern::SHELLS"));
  constexpr int64_t IDS[] = {0, 1, 6};
  for (size_t i = 0; i < std::size(IDS); ++i)
    HS_EXPECT_EQ(pattern->option_values[i], IDS[i]);
  HS_EXPECT_EQ(effect.updateParameter("Pattern", 6), ParamSetResult::APPLIED);
  hs_wasm::collect_param_views(effect, views);
  for (const auto &entry : views)
    if (std::string_view(entry.name) == "Pattern")
      HS_EXPECT_EQ(entry.value, 6);
}

#if HS_ENABLE_PARAM_GUI_BRIDGE
/** @brief The rendered, requested, and accepted streams retain distinct values. */
inline void check_distinct_parameter_values() {
  struct MirroredEffect : StubEffect {
    struct State {
      float value;
    } requested{2.0f}, displayed{3.0f};
    MirroredEffect() : StubEffect(4, 4) {
      register_param("Value", &requested.value, 0.0f, 10.0f);
      mirror_parameter_display_state(requested, displayed);
    }
    float accepted_parameter_value(const ParamDef &) const override {
      return 4.0f;
    }
  } effect;
  std::vector<hs_wasm::ParamView> views;
  hs_wasm::collect_param_views(effect, views);
  HS_EXPECT_EQ(views.size(), size_t{1});
  if (views.size() != 1)
    return;
  HS_EXPECT_EQ(views[0].value, 3.0f);
  HS_EXPECT_EQ(views[0].requested_value, 2.0f);
  HS_EXPECT_EQ(views[0].accepted_value, 4.0f);
}
#endif

inline void check_schema_hook_preserves_written_parameter_identity() {
  struct RebuildingHost : ParamHost {
    float original = 0.0f;
    float replacement = 0.0f;
    bool original_preset;
    int writes = 0;
    explicit RebuildingHost(bool preset) : original_preset(preset) {
      register_param("Original", &original, 0.0f, 1.0f);
      if (!preset)
        mark_global("Original");
      set_parameter_updated_hook([](ParamHost *base, const char *, bool) {
        auto &host = *static_cast<RebuildingHost *>(base);
        host.reset_parameters();
        host.register_param("Replacement", &host.replacement, 0.0f, 1.0f);
        if (host.original_preset)
          host.mark_global("Replacement");
      });
    }
    void parameter_written() override { ++writes; }
  };
  for (bool preset : {false, true}) {
    RebuildingHost host(preset);
    HS_EXPECT_EQ(host.updateParameter("Original", 0.75f),
                 ParamSetResult::APPLIED);
    HS_EXPECT_EQ(host.original, 0.75f);
    HS_EXPECT_EQ(host.replacement, 0.0f);
    HS_EXPECT_EQ(host.writes, preset ? 1 : 0);
  }
}

inline void check_integer_float_endpoints() {
  struct IntegerHost : ParamHost {
    using ParamHost::register_int_param;
  } host;
  int32_t value = 0;
  constexpr int32_t MIN = std::numeric_limits<int32_t>::min();
  constexpr int32_t MAX = 2147483520;
  host.register_int_param("Count", &value, MIN, MAX);
  const auto *def = host.getParameters().find("Count");
  HS_EXPECT_TRUE(def != nullptr);
  if (def == nullptr)
    return;
  HS_EXPECT_EQ(def->min, static_cast<float>(MIN));
  HS_EXPECT_EQ(def->max, static_cast<float>(MAX));
  HS_EXPECT_EQ(host.updateParameter("Count", -3.0e9f), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(value, MIN);
  HS_EXPECT_EQ(host.updateParameter("Count", 3.0e9f), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(value, MAX);
  HS_EXPECT_EQ(host.updateParameter("Count", 0.0f), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(value, 0);
}

struct DefaultFieldsParams {
  float value = 2;
  static constexpr auto FIELDS =
      std::tuple{Control::Field<DefaultFieldsParams, float>{
          "value",
          &DefaultFieldsParams::value,
          "Value",
          {.min = 0, .max = 10, .animated = true}}};
};

class DefaultFieldsEffect
    : public ChoreographedEffect<DefaultFieldsEffect, DefaultFieldsParams> {
  using Base = ChoreographedEffect<DefaultFieldsEffect, DefaultFieldsParams>;

public:
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 3;
  static constexpr std::array<PresetEntry<Params>, 2> PRESETS{
      {{{2}, Segue::Preset::Lerp{2, math::ease_linear, true}},
       {{6}, Segue::Preset::Lerp{2, math::ease_linear, true}}}};
  DefaultFieldsEffect() : Base(96, 48) {}
  void init() override {
    begin_choreography();
    register_described_params();
  }
  void draw_frame() override {}
  void blend_for_test(float progress) {
    transition.from = PRESETS[0].params;
    transition.to = PRESETS[1].params;
    blend_params(progress);
  }
};

inline void test_choreography_descriptions_default() {
  reset_globals();
  DefaultFieldsEffect effect;
  effect.init();
  const auto *definition = effect.getParameters().find("Value");
  HS_EXPECT_TRUE(definition != nullptr);
  HS_EXPECT_EQ(effect.getParameters().size(), 1u);
  if (definition)
    HS_EXPECT_TRUE(definition->animated);
  auto snapshot = effect.serialize_parameters();
  snapshot.params.value = -1;
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
  effect.blend_for_test(.25f);
  HS_EXPECT_EQ(effect.serialize_parameters().params.value, 3);
  effect.blend_for_test(1);
  HS_EXPECT_EQ(effect.serialize_parameters().params.value, 6);
}

struct TypedFieldsState {
  enum class Mode : uint8_t { FIRST, SECOND };
  struct Inner {
    float value = 1;
    int count = 1;
  } inner;
  Mode mode = Mode::FIRST;
  bool enabled = false;
  float phase = .9f;
  float scale = 1;
  float held = 2;
  float telemetry = 0;
  int untabled = 10;
};

inline constexpr auto TYPED_FIELDS = std::tuple{
    Control::FieldGroup{
        &TypedFieldsState::inner,
        std::tuple{Control::Field<TypedFieldsState::Inner, float>{
                       "value",
                       &TypedFieldsState::Inner::value,
                       "Value",
                       {.min = 0, .max = 10}},
                   Control::Field<TypedFieldsState::Inner, int>{
                       "count",
                       &TypedFieldsState::Inner::count,
                       "Count",
                       {.min = 1, .max = 8}}}},
    Control::Field<TypedFieldsState, TypedFieldsState::Mode>{
        "mode", &TypedFieldsState::mode, "Mode", {.min = 0, .max = 1}},
    Control::Field<TypedFieldsState, bool>{
        "enabled", &TypedFieldsState::enabled, "Enabled", {.min = 0, .max = 1}},
    Control::Field<TypedFieldsState, float>{"phase",
                                            &TypedFieldsState::phase,
                                            "Phase",
                                            {.min = 0, .max = 1},
                                            Control::FieldCurve::SHORTEST_TURN},
    Control::Field<TypedFieldsState, float>{"scale",
                                            &TypedFieldsState::scale,
                                            "Scale",
                                            {.min = 1, .max = 16},
                                            Control::FieldCurve::LOG_POSITIVE},
    Control::Field<TypedFieldsState, float>{"held",
                                            &TypedFieldsState::held,
                                            nullptr,
                                            {.min = 0, .max = 10},
                                            Control::FieldCurve::SNAP},
    Control::Field<TypedFieldsState, float>{
        .id = "telemetry",
        .member = &TypedFieldsState::telemetry,
        .name = "Telemetry",
        .spec = {.min = 0, .max = 10, .readonly = true},
        .interpolated = false,
        .validated = false}};

static_assert(Control::valid_fields(TypedFieldsState{}, TYPED_FIELDS));
static_assert(!Control::valid_fields(
    TypedFieldsState{.inner = {.value = 1, .count = 0}}, TYPED_FIELDS));

inline void test_typed_field_invalid_metadata() {
  TypedFieldsState state;
  Control::Field<TypedFieldsState, TypedFieldsState::Mode> mode{
      "mode", &TypedFieldsState::mode, "Mode", {.min = 0, .max = 256}};
  HS_EXPECT_FALSE(mode.valid(state));
  mode.validation_max = 1;
  HS_EXPECT_TRUE(mode.valid(state));
  mode.validation_min = -1;
  HS_EXPECT_FALSE(mode.valid(state));
  Control::Field<TypedFieldsState, bool> enabled{
      "enabled", &TypedFieldsState::enabled, "Enabled", {.min = 0, .max = 2}};
  HS_EXPECT_FALSE(enabled.valid(state));
  Control::Field<TypedFieldsState, float> phase{
      "phase", &TypedFieldsState::phase, "Phase", {.min = 0, .max = 1}};
  phase.validation_min = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_FALSE(phase.valid(state));
  phase.validation_min = 0;
  phase.validation_max = std::numeric_limits<float>::infinity();
  HS_EXPECT_FALSE(phase.valid(state));
  phase.validation_max = -1;
  HS_EXPECT_FALSE(phase.valid(state));
}

inline void test_typed_field_domains_and_exclusions() {
  TypedFieldsState from;
  TypedFieldsState to{.inner = {8, 8},
                      .mode = TypedFieldsState::Mode::SECOND,
                      .enabled = true,
                      .phase = .1f,
                      .scale = 16,
                      .held = 9,
                      .telemetry = 999,
                      .untabled = 20};
  TypedFieldsState out{.inner = {}, .telemetry = 77, .untabled = 42};
  Control::interpolate_fields(out, from, to, .25f, TYPED_FIELDS);
  HS_EXPECT_EQ(out.inner.value, 2.75f);
  HS_EXPECT_EQ(out.inner.count, 1);
  HS_EXPECT_EQ(out.mode, TypedFieldsState::Mode::FIRST);
  HS_EXPECT_FALSE(out.enabled);
  HS_EXPECT_NEAR(out.phase, .95f, 1e-6f);
  HS_EXPECT_NEAR(out.scale, 2.0f, 1e-6f);
  HS_EXPECT_EQ(out.held, 2);
  HS_EXPECT_EQ(out.telemetry, 77);
  HS_EXPECT_EQ(out.untabled, 42);
  Control::interpolate_fields(out, from, to, .5f, TYPED_FIELDS);
  HS_EXPECT_EQ(out.inner.count, 8);
  HS_EXPECT_EQ(out.mode, TypedFieldsState::Mode::SECOND);
  HS_EXPECT_TRUE(out.enabled);
  Control::interpolate_fields(out, from, to, 1.0f, TYPED_FIELDS);
  HS_EXPECT_EQ(out.inner.value, 8);
  HS_EXPECT_EQ(out.phase, .1f);
  HS_EXPECT_EQ(out.scale, 16);
  HS_EXPECT_EQ(out.held, 9);
  HS_EXPECT_EQ(out.telemetry, 77);
  HS_EXPECT_TRUE(Control::valid_fields(to, TYPED_FIELDS));
  to.inner.value = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_FALSE(Control::valid_fields(to, TYPED_FIELDS));
  to.inner.value = std::numeric_limits<float>::infinity();
  HS_EXPECT_FALSE(Control::valid_fields(to, TYPED_FIELDS));
  to.inner.value = 8;
  to.mode = static_cast<TypedFieldsState::Mode>(2);
  HS_EXPECT_FALSE(Control::valid_fields(to, TYPED_FIELDS));

  std::vector<std::string_view> names;
  Control::register_fields(
      out, TYPED_FIELDS, [&](const char *name, auto *target, const auto &spec) {
        names.emplace_back(name);
        HS_EXPECT_TRUE(target != nullptr);
        if (std::string_view(name) == "Telemetry")
          HS_EXPECT_TRUE(spec.readonly);
      });
  const std::vector<std::string_view> expected{
      "Value", "Count", "Mode", "Enabled", "Phase", "Scale", "Telemetry"};
  HS_EXPECT_TRUE(names == expected);
}

template <typename E>
inline void check_described_snapshot_ranges(const char *name) {
  HS_CONTEXT(name);
  reset_globals();
  E effect;
  effect.init();
  const auto original = effect.serialize_parameters();
  auto candidate = original;
  size_t floats = 0;
  Control::register_fields(
      candidate.params, E::parameter_fields(),
      [&](const char *field_name, auto *target, const auto &spec) {
        HS_CONTEXT(field_name);
        using Value = std::remove_pointer_t<decltype(target)>;
        if constexpr (std::is_same_v<Value, float>) {
          if (spec.readonly)
            return;
          ++floats;
          *target = std::numeric_limits<float>::quiet_NaN();
          HS_EXPECT_FALSE(effect.restore_parameters(candidate));
          candidate = original;
          *target = std::numeric_limits<float>::infinity();
          HS_EXPECT_FALSE(effect.restore_parameters(candidate));
          candidate = original;
          *target = spec.min - std::max(1.0f, std::fabs(spec.min)) * .125f;
          HS_EXPECT_FALSE(effect.restore_parameters(candidate));
          candidate = original;
          *target = std::numeric_limits<float>::max();
          HS_EXPECT_FALSE(effect.restore_parameters(candidate));
          candidate = original;
        }
      });
  HS_EXPECT_GT(floats, 0u);
  HS_EXPECT_TRUE(effect.restore_parameters(original));
}

inline void test_authored_field_snapshot_validation() {
  size_t described = 0;
  auto check = [&]<template <int, int> class E>(const char *name) {
    if constexpr (requires { E<96, 48>::parameter_fields(); }) {
      check_described_snapshot_ranges<E<96, 48>>(name);
      ++described;
    }
  };
#define HS_CHECK_DESCRIBED_SNAPSHOT(E) check.template operator()<E>(#E);
  HS_EFFECT_LIST(HS_CHECK_DESCRIBED_SNAPSHOT)
#undef HS_CHECK_DESCRIBED_SNAPSHOT
  HS_EXPECT_GE(described, size_t{8});

  MindSplatterParams from;
  MindSplatterParams to;
  from.friction = .6f;
  to.friction = .9f;
  from.base_mesh = Solids::BaseMesh::CUBE;
  to.base_mesh = Solids::BaseMesh::ICOSAHEDRON;
  MindSplatterParams out{.active_count = 123};
  out.lerp(from, to, .5f);
  HS_EXPECT_EQ(out.friction,
               from.friction + (to.friction - from.friction) * .5f);
  HS_EXPECT_EQ(out.base_mesh, to.base_mesh);
  HS_EXPECT_EQ(out.active_count, 123);
}

/** @brief Sparse enum IDs pass through ParamView; dense enums have no ID table. */
inline void test_sparse_option_metadata() {
  struct SparseEffect : Effect {
    using Effect::register_param;
    SparseEffect() : Effect(4, 4) {}
    void draw_frame() override {}
  } effect;
  static constexpr const char *LABELS[] = {"Zero", "One", "Six"};
  static constexpr const char *EXPORTS[] = {"Mode::ZERO", "Mode::ONE",
                                            "Mode::SIX"};
  static constexpr const int64_t IDS[] = {0, 1, 6};
  uint8_t mapped = 6;
  uint8_t dense = 2;
  effect.register_param("Mapped", &mapped,
                        ParamSpec<uint8_t>{.min = 0,
                                           .max = 6,
                                           .options = LABELS,
                                           .export_options = EXPORTS,
                                           .option_count = 3,
                                           .option_values = IDS});
  effect.register_param("Dense", &dense,
                        ParamSpec<uint8_t>::enumerated(LABELS, 3));
  std::vector<hs_wasm::ParamView> views;
  hs_wasm::collect_param_views(effect, views);
  HS_EXPECT_EQ(views.size(), size_t(2));
  if (views.size() != 2)
    return;
  HS_EXPECT_TRUE(views[0].options == LABELS);
  HS_EXPECT_TRUE(views[0].export_options == EXPORTS);
  HS_EXPECT_TRUE(views[0].option_values == IDS);
  HS_EXPECT_EQ(views[0].option_count, int(std::size(IDS)));
  HS_EXPECT_EQ(views[0].value, 6.0f);
  HS_EXPECT_EQ(views[0].max, 6.0f);
  HS_EXPECT_TRUE(views[1].option_values == nullptr);
  HS_EXPECT_EQ(views[1].value, 2.0f);
}

/**
 * @brief Runs the param-marshal module.
 * @return The module's failure count.
 */
inline int run_param_marshal_tests() {
  hs_test::ModuleFixture fixture("param_marshal");
  test_choreography_descriptions_default();
  test_typed_field_invalid_metadata();
  test_typed_field_domains_and_exclusions();
  test_authored_field_snapshot_validation();
  test_sparse_option_metadata();
  check_roster_order_pinned();
  check_generation_tracker();
  check_hyper_lattice_pattern_view_dropdowns();
  check_integer_float_endpoints();
#if HS_ENABLE_PARAM_GUI_BRIDGE
  check_distinct_parameter_values();
#endif
  check_schema_hook_preserves_written_parameter_identity();
  // Tally how many effects exercised the by-name round-trip; it is skipped for
  // effects with no editable float param. Surface the split and fail if zero.
  int rt_covered = 0, rt_total = 0, rt_skipped = 0, rt_mismatched = 0;
  FieldCoverage coverage;
#define HS_PARAM_ONE(name)                                                     \
  do {                                                                         \
    ++rt_total;                                                                \
    const auto result = check_one<name>(#name, coverage);                      \
    if (result == RoundTripResult::COVERED)                                    \
      ++rt_covered;                                                            \
    else if (result == RoundTripResult::NO_EDITABLE_FLOAT)                     \
      ++rt_skipped;                                                            \
    else                                                                       \
      ++rt_mismatched;                                                         \
  } while (0);
  HS_EFFECT_LIST(HS_PARAM_ONE)
#undef HS_PARAM_ONE
  HS_EXPECT(coverage.range_distinguishing,
            "no roster param has min != max — a transposed min/max pair in "
            "ParamView would ride green");
  HS_EXPECT(coverage.flags_distinguishing,
            "no roster param has animated != readonly — a transposed "
            "animated/readonly pair in ParamView would ride green");
  HS_EXPECT(coverage.float_enum_present,
            "no roster param is a float-backed enum — the step-presence check "
            "above would ride green on ParamDef::is_integer() alone");
  std::printf("  param-marshal by-name round-trip exercised on %d/%d effects "
              "(%d skipped: no editable float param; %d stream mismatches)\n",
              rt_covered, rt_total, rt_skipped, rt_mismatched);
  HS_EXPECT(rt_covered > 0,
            "by-name round-trip must run on at least one effect — the roster "
            "drifted to all-non-editable params and the check covers nothing");

  // Marshal the roster through the engine's fixed-capacity ParamStreams.
  size_t max_count = 0;
#define HS_PARAM_COUNT(name) count_one<name>(max_count);
  HS_EFFECT_LIST(HS_PARAM_COUNT)
#undef HS_PARAM_COUNT

  hs_wasm::ParamStreams streams;
  auto &views = streams.views;
  auto &values = streams.values;
  HS_EXPECT_LE(max_count, hs_wasm::ParamStreams::CAPACITY);
  HS_EXPECT_GE(values.capacity(), hs_wasm::ParamStreams::CAPACITY);
  HS_EXPECT_GE(views.capacity(), hs_wasm::ParamStreams::CAPACITY);
  const hs_wasm::ParamView *view_data = views.data();
  const float *value_data = values.data();
  const size_t view_cap = views.capacity();
  const size_t value_cap = values.capacity();

#define HS_PARAM_STAB(name)                                                    \
  check_stability_one<name>(#name, views, values, view_data, view_cap,         \
                            value_data, value_cap);
  HS_EFFECT_LIST(HS_PARAM_STAB)
#undef HS_PARAM_STAB

  return fixture.result();
}

} // namespace param_marshal_tests
} // namespace hs_test
