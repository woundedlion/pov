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

inline void test_mobius_automatic_departures_retain_lens() {
  for (const bool fade : {false, true}) {
    reset_effect_globals();
    auto effect = std::make_unique<MobiusFrameProbe>();
    effect->init();
    const math::MobiusParams live{1.0f, 0.2f, 0.1f, 0.3f,
                                  0.2f, 0.1f, 0.9f, -0.1f};
    effect->params.template get<"lens">().mobius = live;
    const auto target = preset_params_or_initial<MobiusFrameProbe>(1);
    auto &transition = effect->transition;
    transition.from = effect->params;
    transition.to = target;
    transition.frames = 3;
    transition.elapsed_frames = 0;
    transition.fades = fade;
    transition.adopted = false;
    transition.active = true;
    for (const float progress : {0.25f, 0.5f, 1.0f}) {
      effect->run_transition(progress);
      verify_mobius_equal(effect->params.template get<"lens">().mobius, live);
    }
    auto expected = target;
    expected.template get<"lens">().mobius = live;
    verify_params_equal(effect->params, expected);
    HS_EXPECT_FALSE(transition.active);
    HS_EXPECT_TRUE(effect->synchronizePreset(1));
    verify_params_equal(effect->params, target);
  }
}

inline void test_mobius_captured_parameter_restore() {
#if HS_ENABLE_PARAM_GUI_BRIDGE
  for (const math::MobiusParams target :
       {math::MobiusParams{0, 0, 0, 0.704f, 0.304f, 0, 0, 0},
        math::MobiusParams{0.7071f, -1, 0.7071f, -1, 0.7071f, 0, 0, 0}}) {
    reset_effect_globals();
    auto effect = std::make_unique<MobiusFrameProbe>();
    effect->init();
    HS_EXPECT_TRUE(Pullback::MobiusLensParams::nondegenerate(target));
    effect->params.template get<"lens">().mobius = target;
    std::vector<std::pair<std::string, float>> captured;
    for (const auto &def : effect->getParameters())
      if (!def.readonly)
        captured.emplace_back(def.name, def.get_requested());
    HS_EXPECT_TRUE(effect->synchronizePreset(0));
    effect->replay_parameter_writes(captured);
    verify_mobius_equal(effect->params.template get<"lens">().mobius, target);
    HS_EXPECT_TRUE(effect->animations_paused());
    for (const auto &[name, value] : captured) {
      const auto *def = effect->getParameters().find(name.c_str());
      HS_EXPECT_TRUE(def != nullptr);
      if (def != nullptr)
        HS_EXPECT_EQ(def->get_requested(), value);
    }
  }
#endif
}

inline void test_mobius_frame_admission() {
  reset_effect_globals();
  auto effect = std::make_unique<MobiusFrameProbe>();
  effect->init();
#if HS_ENABLE_PARAM_GUI_BRIDGE
  effect->params.template get<"lens">().mobius = math::MobiusParams{};
  const auto captured = effect->serialize_parameters();
  HS_EXPECT_FALSE(effect->animations_paused());
  HS_EXPECT_EQ(effect->updateParameter("Mobius A Re", 0.0f),
               ParamSetResult::INADMISSIBLE);
  verify_params_equal(effect->serialize_parameters().params, captured.params);
  HS_EXPECT_FALSE(effect->animations_paused());
  HS_EXPECT_TRUE(effect->parameter_warning("Mobius A Re") != nullptr);
  HS_EXPECT_TRUE(effect->parameter_warning("Mobius B Re") == nullptr);
  HS_EXPECT_EQ(effect->updateParameter("Mobius A Re", 0.5f),
               ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(effect->animations_paused());
  HS_EXPECT_TRUE(effect->parameter_warning("Mobius A Re") == nullptr);
  HS_EXPECT_TRUE(effect->restore_parameters(effect->serialize_parameters()));
#endif
  for (const math::MobiusParams bad :
       {math::MobiusParams{1, 0, 1, 0, 1, 0, 1, 0},
        math::MobiusParams{std::numeric_limits<float>::quiet_NaN(), 0, 0, 0, 0,
                           0, 1, 0}}) {
    effect->params.template get<"lens">().mobius = bad;
    const auto frame = effect->frame_for_test();
    verify_mobius_equal(frame.params.template get<"lens">().mobius,
                        math::MobiusParams{});
    for (const math::Vector view :
         {math::Vector{1, 0, 0}, math::Vector{0, 1, 0},
          math::Vector{0, 0, 1}}) {
      const auto output = math::mobius_transform(
          view, frame.params.template get<"lens">().mobius);
      HS_EXPECT_TRUE(std::isfinite(output.x) && std::isfinite(output.y) &&
                     std::isfinite(output.z));
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

/** @brief Pins the PentBright initial lattice scale. */
inline void test_composed_pentbright_lattice_scale() {
  const auto polar = KaleidoscopePentBright<SMALL_W, SMALL_H>::initial_params();
  HS_EXPECT_NEAR(polar.template get<"source">().lattice_cell_scale *
                     math::TWO_PI_F,
                 5.0f, 1e-6f);
}

/** @brief Sweeps the preset choreography over every specialization. */
inline void test_composed_preset_choreography() {
#define HS_COMPOSED_PRESETS(name, seconds)                                     \
  check_preset_choreography<name>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_PRESETS)
#undef HS_COMPOSED_PRESETS
}

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
 * @brief Pins the two families the base registers by hand, not by table.
 * @details Color fields have no slider names; the composed Mobius lens has no
 * field table.
 */
inline void test_composed_hand_registered_families() {
  for (const auto &field : Pullback::ColorParams::FIELDS) {
    HS_CONTEXT(field.id);
    HS_EXPECT_TRUE(field.name == nullptr);
  }
  static_assert(!Pullback::HasFields<Pullback::MobiusLensParams>);
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
      KaleidoscopeHexOilSpec::PROJECTION, FX::Binding>::Type>;
  static_assert(std::is_same_v<typename FX::RenderPipeline::template node_at<1>,
                               Pullback::Detail::BoundPlaced<
                                   Pullback::CodeEmission::OUT_OF_LINE_FLASH,
                                   FX::Binding, Lens, Displace, Project>>);
  HS_EXPECT_EQ(TraitsOf<FX>::SURFACE_PLACEMENT,
               Pullback::SurfacePlacement::AFTER_LENS);
}
