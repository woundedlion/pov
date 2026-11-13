/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Twist clocks advance independently of preset pause and sample across wraps. */
inline void test_grid_worms_twist_clocks() {
  using FX = GridWorms<SMALL_W, SMALL_H>;
  using Leaves = Pullback::ComposedDetail::PipelineLeaves<
      GridWormsSpec::Pipeline<Pullback::ComposedDetail::DiscoveryBinding>>;
  using Discovery = Pullback::ComposedDetail::DiscoveryBinding;
  static_assert(
      std::is_same_v<typename Leaves::template StageAt<0>,
                     Pullback::Stage::Lens<
                         Pullback::Lens::AnimatedTwist<
                             Pullback::TwistProvider<Discovery, "lens1">>,
                         false>>);
  static_assert(
      std::is_same_v<
          typename Leaves::template StageAt<1>,
          Pullback::Stage::Rotate<Pullback::OuterCameraProvider<Discovery>>>);
  static_assert(!FX::factored_columns());
  reset_globals();
  FX effect;
  effect.setAnimationsPaused(true);
  effect.init();
  const auto captured = effect.serialize_parameters();
  const auto &params = captured.params.template get<"lens1">();
  for (int frame = 1; frame <= 241; ++frame) {
    ComposedFrameWhiteBox::advance(effect);
    if (frame == 1 || frame == 9 || frame == 241) {
      const auto now = ComposedFrameWhiteBox::frame(effect);
      const auto &clock = now.resources.template get<"lens1">();
      HS_EXPECT_NEAR(clock.phi,
                     fmod(frame * double(params.speed), double(math::TWO_PI_F)),
                     5e-5);
      HS_EXPECT_NEAR(clock.lambda,
                     fmod(frame * double(params.drift), double(math::TWO_PI_F)),
                     5e-5);
    }
  }
  auto from = ComposedFrameWhiteBox::frame(effect);
  const auto saved_clock = from.resources.template get<"lens1">();
  HS_EXPECT_TRUE(effect.restore_parameters(captured));
  const auto restored = ComposedFrameWhiteBox::frame(effect);
  HS_EXPECT_EQ(restored.resources.template get<"lens1">().phi, saved_clock.phi);
  HS_EXPECT_EQ(restored.resources.template get<"lens1">().lambda,
               saved_clock.lambda);
  verify_params_equal(effect.serialize_parameters().params, captured.params);
  auto to = from;
  from.resources.template get<"lens1">().phi = math::TWO_PI_F - 0.02f;
  from.resources.template get<"lens1">().lambda = math::TWO_PI_F - 0.03f;
  to.resources.template get<"lens1">().phi = 0.03f;
  to.resources.template get<"lens1">().lambda = 0.05f;
  using Policy = Pullback::Lens::AnimatedTwist<
      Pullback::TwistProvider<typename FX::Binding, "lens1">>;
  for (const float fraction : {0.0f, 0.125f, 0.5f, 0.875f, 1.0f}) {
    const auto column = Pullback::FrameInterval(from, to).at(fraction);
    const auto &clock = column.resources.template get<"lens1">();
    const double phi = double(math::TWO_PI_F) - 0.02 + 0.05 * fraction;
    const double lambda = double(math::TWO_PI_F) - 0.03 + 0.08 * fraction;
    HS_EXPECT_NEAR(sin(clock.phi), sin(phi), 1e-6);
    HS_EXPECT_NEAR(sin(clock.lambda), sin(lambda), 1e-6);
    const auto prepared = Policy::prepare(column);
    HS_EXPECT_NEAR(prepared.axis.x, sin(phi) * cos(lambda), 1e-6);
    HS_EXPECT_NEAR(prepared.axis.y, cos(phi), 1e-6);
    HS_EXPECT_NEAR(prepared.axis.z, sin(phi) * sin(lambda), 1e-6);
  }
}

/** @brief Authored GridWorms pixels and controls match its document-built chain. */
inline void test_grid_worms_document_pipeline() {
  using FX = GridWorms<SMALL_W, SMALL_H>;
  using namespace shader_chain_tests;
  auto before_compare = [](auto &own) {
    own.outer_conjugate =
        math::make_rotation(math::Vector(0.3f, 0.7f, -0.2f).normalized(), 0.8f)
            .conjugate();
  };
  auto after_static = [](auto &program, const auto &context, const auto &own,
                         auto &effect, int channel_tolerance,
                         float alpha_tolerance) {
    auto blocks = std::make_unique<TemporalBlocks>();
    program.bind_temporal_storage(blocks->history, blocks->column);
    program.prepare(context);
    program.capture_history();
    auto &twist =
        shader_chain_tests::state_as<In::Op::TwistOrbitState>(program, 0);
    const auto &twist_params =
        shader_chain_tests::param_as<In::Op::TwistChainParams>(program, 0);
    In::Op::LensTwist::advance(twist, twist_params);
    ComposedFrameWhiteBox::advance(effect);
    auto next = ComposedFrameWhiteBox::frame(effect);
    next.outer_conjugate = own.outer_conjugate;
    const auto &clock = next.resources.template get<"lens1">();
    HS_EXPECT_EQ(twist.phi, clock.phi);
    HS_EXPECT_EQ(twist.lambda, clock.lambda);
    HS_EXPECT_TRUE(program.history_continuous());
    for (const float fraction : {0.0f, 0.125f, 0.5f, 0.875f, 1.0f}) {
      const auto column = Pullback::FrameInterval(own, next).at(fraction);
      const auto frame = FX::RenderPipeline::prepare(column);
      program.prepare_column(context, context, fraction);
      for (const auto &view : sweep_views()) {
        const auto expected = FX::shade(view, frame);
        const auto actual =
            program.evaluate_column(view, context, context, fraction);
        expect_color_within(actual, expected, channel_tolerance,
                            alpha_tolerance);
      }
    }
  };
  composed_chain_tests::verify_export<FX>(0, 1, 1e-6f, before_compare,
                                          after_static);
}
