/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Roster-wide slider, snapshot, choreography, interpolation and document contracts
// ============================================================================

/** @brief Sweeps the registered slider set over every specialization. */
inline void test_composed_slider_registration() {
#define HS_COMPOSED_SLIDERS(name, seconds)                                     \
  check_slider_registration<name>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_SLIDERS)
#undef HS_COMPOSED_SLIDERS
}

struct MobiusFrameProbe : MobiusGrid<SMALL_W, SMALL_H> {
  using MobiusGrid<SMALL_W, SMALL_H>::params;
  using MobiusGrid<SMALL_W, SMALL_H>::run_transition;
  using MobiusGrid<SMALL_W, SMALL_H>::transition;
};

/** @brief Mobius preset departures adopt authored generators. */
inline void test_mobius_automatic_departures_interpolate_generator() {
  for (bool fade : {false, true}) {
    reset_effect_globals();
    MobiusFrameProbe effect;
    effect.init();
    effect.params.template get<"lens">().spin_x = 2.0f;
    const auto target = preset_params<MobiusFrameProbe>(1);
    effect.transition.from = effect.params;
    effect.transition.to = target;
    effect.transition.frames = 3;
    effect.transition.elapsed_frames = 0;
    effect.transition.fades = fade;
    effect.transition.adopted = false;
    effect.transition.active = true;
    for (const float progress : {0.25f, 0.5f, 1.0f})
      effect.run_transition(progress);
    verify_params_equal(effect.params, target);
    HS_EXPECT_FALSE(effect.transition.active);
  }
}

/** @brief Captured generator sliders replay exactly. */
inline void test_mobius_captured_parameter_restore() {
#if HS_ENABLE_PARAM_GUI_BRIDGE
  reset_effect_globals();
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
  reset_effect_globals();
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

/** @brief Positive logarithmic parameters interpolate within their declared domain. */
inline void test_composed_log_positive_curve() {
  using Pullback::FieldCurve;
  using Pullback::Fields::apply_curve;
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

/** @brief Direct-noise displacement follows the declared lens placement. */
inline void test_composed_direct_surface_placement() {
  using FX = KaleidoscopeHexOil<SMALL_W, SMALL_H>;
  using Displace = Pullback::Stage::Displace<Pullback::Surface::DirectNoise<
      Pullback::SurfaceProvider<FX::Binding, Pullback::DirectSurfaceParams,
                                true>,
      math::NoiseBasis::SIMPLEX>>;
  using Lens =
      Pullback::Stage::Lens<Pullback::Lens::HexagonalPrismKaleidoscope>;
  using Project = Pullback::Stage::Project<Pullback::ProjectionPolicyFor<
      Pullback::MetadataFor<KaleidoscopeHexOilSpec>::PROJECTION,
      FX::Binding>::Type>;
  using Group = typename FX::RenderPipeline::template NodeAt<1>;
  static_assert(std::is_same_v<
                typename Group::BoundTuple,
                std::tuple<Lens::Bind<FX::Binding>, Displace::Bind<FX::Binding>,
                           Project::Bind<FX::Binding>>>);
  HS_EXPECT_EQ(TraitsOf<FX>::SURFACE_PLACEMENT,
               Pullback::SurfacePlacement::AFTER_LENS);
}
