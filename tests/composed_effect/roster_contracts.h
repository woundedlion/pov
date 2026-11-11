/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Roster-wide slider, snapshot, choreography, interpolation and document contracts
// ============================================================================

/** @brief Sweeps the registered slider set over every specialization. */
inline void test_composed_slider_registration() {
#define HS_COMPOSED_SLIDERS(name) check_slider_registration<name>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_SLIDERS)
#undef HS_COMPOSED_SLIDERS
}

/** @brief Blends wave directions independently of their lengths. */
inline void test_composed_modulation_blends_axes() {
  using Pullback::Modulation::Config;
  using Pullback::Modulation::interpolate;
  const auto check = [](Pattern::SourceKind source, auto member) {
    Config start, end;
    start.source = end.source = source;
    start.evaluation = end.evaluation = Pattern::Evaluation::DIRECT;
    (start.*member).axis = {4.0f, 0.0f, 0.0f};
    (end.*member).axis = {0.0f, 2.0f, 0.0f};
    if constexpr (requires { (start.*member).axis_b; }) {
      (start.*member).axis_b = {0.0f, 0.0f, 4.0f};
      (end.*member).axis_b = {2.0f, 0.0f, 0.0f};
    }
    const auto mid = interpolate(start, end, 0.5f);
    HS_EXPECT_TRUE(Pullback::Modulation::valid(mid));
    HS_EXPECT_VEC((mid.*member).axis, math::Vector(0.5f, 0.5f, 0.0f), 1e-6f);
    if constexpr (requires { (mid.*member).axis_b; }) {
      HS_EXPECT_VEC((mid.*member).axis_b, math::Vector(0.5f, 0.0f, 0.5f),
                    1e-6f);
      (end.*member).axis_b = {0.0f, 0.0f, -2.0f};
      const auto antipodal = interpolate(start, end, 0.5f);
      HS_EXPECT_TRUE(Pullback::Modulation::valid(antipodal));
      HS_EXPECT_EQ((antipodal.*member).axis_b, (end.*member).axis_b);
    }
    HS_EXPECT_EQ((interpolate(start, end, 0.0f).*member).axis,
                 (start.*member).axis);
    HS_EXPECT_EQ((interpolate(start, end, 1.0f).*member).axis,
                 (end.*member).axis);
    HS_EXPECT_EQ((interpolate(start, start, 0.5f).*member).axis,
                 (start.*member).axis);
  };
  check(Pattern::SourceKind::DIRECTIONAL_SINE, &Config::sine);
  check(Pattern::SourceKind::SPHERE_INTERFERENCE, &Config::interference);
  check(Pattern::SourceKind::SPHERE_SMOOTH_GRID, &Config::grid);
}

/** @brief Every shipped default and preset uses a qualified modulation configuration. */
inline void test_composed_modulation_is_device_qualified() {
  const auto check = []<template <int, int> class E>() {
    using FX = E<SMALL_W, SMALL_H>;
    const auto params_qualified = [](const auto &params) {
      params.visit([]<typename Resource>(const auto &family) {
        if constexpr (requires { family.modulation; })
          HS_EXPECT_TRUE(
              Pullback::Modulation::device_qualified(family.modulation));
      });
    };
    params_qualified(FX::initial_params());
    for (size_t index = 0; index < FX::PRESET_IDS.size(); ++index)
      params_qualified(preset_params<FX>(index));
  };
#define HS_COMPOSED_MODULATION(name) check.template operator()<name>();
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_MODULATION)
#undef HS_COMPOSED_MODULATION
}

struct MobiusFrameProbe : MobiusGrid<SMALL_W, SMALL_H> {
  using MobiusGrid<SMALL_W, SMALL_H>::params;
  using MobiusGrid<SMALL_W, SMALL_H>::run_transition;
  using MobiusGrid<SMALL_W, SMALL_H>::transition;
};

/** @brief Mobius departures interpolate generators or hold them during fades. */
inline void test_mobius_automatic_departures_interpolate_generator() {
  for (bool fade : {false, true}) {
    reset_globals();
    MobiusFrameProbe effect;
    effect.init();
    effect.params.template get<"lens">().spin_x = 2.0f;
    const auto from = effect.params;
    const auto target = preset_params<MobiusFrameProbe>(1);
    effect.transition.from = effect.params;
    effect.transition.to = target;
    effect.transition.frames = 3;
    effect.transition.elapsed_frames = 0;
    effect.transition.fades = fade;
    effect.transition.adopted = false;
    effect.transition.active = true;
    effect.run_transition(0.25f);
    if (fade)
      verify_params_equal(effect.params, from);
    else
      HS_EXPECT_NEAR(effect.params.template get<"lens">().spin_x,
                     Pullback::interpolate(from, target, 0.25f)
                         .template get<"lens">()
                         .spin_x,
                     1e-6f);
    for (const float progress : {0.5f, 1.0f})
      effect.run_transition(progress);
    verify_params_equal(effect.params, target);
    HS_EXPECT_FALSE(effect.transition.active);
  }
}

/** @brief Captured generator sliders replay exactly. */
inline void test_mobius_captured_parameter_restore() {
#if HS_ENABLE_PARAM_GUI_BRIDGE
  reset_globals();
  MobiusFrameProbe effect;
  effect.init();
  effect.params.template get<"lens">().spin_x = 1.25f;
  effect.params.template get<"lens">().push_z = -2.5f;
  const auto expected = effect.params;
  std::vector<std::pair<std::string, float>> captured;
  for (const auto &def : effect.getParameters())
    if (!def.readonly)
      captured.emplace_back(def.name, def.get_requested());
  HS_EXPECT_TRUE(effect.synchronizePreset(0));
  effect.replay_parameter_writes(captured);
  verify_params_equal(effect.params, expected);
#endif
}

/** @brief Snapshot admission rejects generators outside their ranges. */
inline void test_mobius_frame_admission() {
  reset_globals();
  MobiusFrameProbe effect;
  effect.init();
  const auto captured = effect.serialize_parameters();
  for (const float value : {7.0f, std::numeric_limits<float>::quiet_NaN()}) {
    auto invalid = captured;
    invalid.params.template get<"lens">().spin_x = value;
    HS_EXPECT_FALSE(effect.restore_parameters(invalid));
    verify_params_equal(effect.serialize_parameters().params, captured.params);
  }
  auto previous = captured;
  previous.schema_version = 4;
  HS_EXPECT_FALSE(effect.restore_parameters(previous));
}

/** @brief Generator interpolation produces finite, continuous sphere maps. */
inline void test_mobius_generator_interpolation() {
  using Params = Pullback::MobiusGeneratorParams;
  for (int pair = 0; pair < 40; ++pair) {
    Params a, b;
    for (size_t i = 0; i < Params::FIELDS.size(); ++i) {
      a.*Params::FIELDS[i].member = sinf(pair * 1.7f + i) * 2.0f;
      b.*Params::FIELDS[i].member = cosf(pair * 0.9f + i) * 2.0f;
    }
    const auto first =
        Pullback::Lens::prepare_mobius(Pullback::interpolate(a, b, 0.5f));
    const auto second =
        Pullback::Lens::prepare_mobius(Pullback::interpolate(a, b, 0.50001f));
    for (int direction = 0; direction < 40; ++direction) {
      const auto p =
          math::Vector::from_spherical(direction * 0.618f, direction * 0.078f);
      const auto u = math::mobius_transform(p, first);
      const auto v = math::mobius_transform(p, second);
      HS_EXPECT_TRUE(finite_vec(u));
      HS_EXPECT_TRUE(finite_vec(v));
      HS_EXPECT_LT((u - v).magnitude(), 0.003f);
    }
  }
}

/** @brief Sweeps the parameter snapshot contract over every specialization. */
inline void test_composed_snapshot_contract() {
#define HS_COMPOSED_SNAPSHOT(name) check_snapshot_contract<name>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_SNAPSHOT)
#undef HS_COMPOSED_SNAPSHOT
}

/** @brief Sweeps the preset choreography over every specialization. */
inline void test_composed_preset_choreography() {
#define HS_COMPOSED_PRESETS(name) check_preset_choreography<name>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_PRESETS)
#undef HS_COMPOSED_PRESETS
}

/** @brief Positive logarithmic parameters interpolate within their declared domain. */
inline void test_composed_log_positive_curve() {
  using Pullback::FieldCurve;
  using Control::apply_curve;
  HS_EXPECT_EQ(apply_curve(FieldCurve::LOG_POSITIVE, 1.0f, 81.0f, 0.0f), 1.0f);
  HS_EXPECT_EQ(apply_curve(FieldCurve::LOG_POSITIVE, 1.0f, 81.0f, 1.0f), 81.0f);
  HS_EXPECT_NEAR(apply_curve(FieldCurve::LOG_POSITIVE, 1.0f, 81.0f, 0.25f),
                 3.0f, 1e-5f);
  HS_EXPECT_NEAR(apply_curve(FieldCurve::LOG_POSITIVE, 1.0f, 81.0f, 0.5f), 9.0f,
                 1e-5f);
  HS_EXPECT_NEAR(apply_curve(FieldCurve::LOG_POSITIVE, 81.0f, 1.0f, 0.75f),
                 3.0f, 1e-5f);
}

/** @brief Sweeps the crossfade interpolation over every specialization. */
inline void test_composed_preset_interpolation() {
#define HS_COMPOSED_INTERP(name) check_preset_interpolation<name>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_INTERP)
#undef HS_COMPOSED_INTERP
}

/** @brief Sweeps the shader-document value pin over every specialization. */
inline void test_composed_document_values() {
#define HS_COMPOSED_DOCUMENTS(name) check_document_values<name>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_DOCUMENTS)
#undef HS_COMPOSED_DOCUMENTS
}

/**
 * @brief Pins color registration and the generator field table.
 */
inline void test_composed_hand_registered_families() {
  for (const auto &field : Pullback::ColorParams::FIELDS) {
    HS_CONTEXT(field.id);
    HS_EXPECT_TRUE(field.name == nullptr);
  }
  static_assert(Pullback::HasFields<Pullback::MobiusGeneratorParams>);
  for (const ColorSliderBinding &binding : COLOR_SLIDER_BINDINGS) {
    HS_CONTEXT(binding.slider);
    HS_EXPECT_TRUE(find_field<Pullback::ColorParams>(binding.field_id) !=
                   nullptr);
  }
  HS_EXPECT_EQ(std::size(COLOR_SLIDER_BINDINGS),
               Pullback::ColorParams::FIELDS.size());
}
