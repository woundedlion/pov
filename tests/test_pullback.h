/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <cstddef>
#include <type_traits>

#include "core/render/pullback.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace pullback_tests {

struct TestFrame {
  mutable std::array<Pullback::StageKind, 6> calls{};
  mutable size_t call_count = 0;
  float edge_width = 0.5f;
};

struct TestBinding {
  using FrameState = TestFrame;
  using Instrumentation = Pullback::NoInstrumentation;
};

template <Pullback::StageKind KindV, typename InputT, typename OutputT,
          bool TerminalV = false>
struct ExactStage {
  using Binding = TestBinding;
  using FrameState = typename Binding::FrameState;
  using Input = InputT;
  using Output = OutputT;

  static constexpr Pullback::StageKind KIND = KindV;
  static constexpr Pullback::CodeEmission EMISSION =
      Pullback::CodeEmission::INLINE_ONLY;
  static constexpr bool APPROXIMATE = false;
  static constexpr bool TERMINAL = TerminalV;
  static constexpr bool NON_FLOATING_FIELDS_EXACT = true;
  static constexpr Pullback::ApproximationOracleId ORACLE =
      Pullback::ApproximationOracleId::NONE;
  static constexpr std::array<Pullback::ApproximationMetric, 0> METRICS{};

  static void record(const FrameState &frame) {
    frame.calls[frame.call_count++] = KIND;
  }
};

struct OuterStage
    : ExactStage<Pullback::StageKind::OUTER_CAMERA, Vector, Vector> {
  static Vector run(const Vector &input, const FrameState &frame) {
    record(frame);
    return Vector(input.x + 1.0f, input.y, input.z);
  }
};

struct SurfaceStage : ExactStage<Pullback::StageKind::SURFACE_PROJECT, Vector,
                                 Pullback::ProjectionSample> {
  static Pullback::ProjectionSample run(const Vector &input,
                                        const FrameState &frame) {
    record(frame);
    return {Complex(input.x + 1.0f, input.y + 2.0f),
            3,
            4,
            5,
            6.0f,
            0.75f,
            7,
            8,
            9,
            0.5f,
            input,
            10.0f};
  }
};

struct WarpStage
    : ExactStage<Pullback::StageKind::PLANAR_WARP, Pullback::ProjectionSample,
                 Pullback::SourceInput> {
  static Pullback::SourceInput run(const Pullback::ProjectionSample &input,
                                   const FrameState &frame) {
    record(frame);
    return {input,
            {Complex(input.coords.re + 2.0f, input.coords.im + 3.0f),
             Complex(2.0f, 3.0f), 4.0f}};
  }
};

struct SourceStage
    : ExactStage<Pullback::StageKind::SOURCE, Pullback::SourceInput,
                 Pullback::MaterialInput> {
  static Pullback::MaterialInput run(const Pullback::SourceInput &input,
                                     const FrameState &frame) {
    record(frame);
    return {input.projected, input.warped,
            input.warped.coords.re + input.warped.coords.im};
  }
};

struct MaterialStage
    : ExactStage<Pullback::StageKind::MATERIAL, Pullback::MaterialInput,
                 Pullback::MaterialSample> {
  static Pullback::MaterialSample run(const Pullback::MaterialInput &input,
                                      const FrameState &frame) {
    record(frame);
    return {input.field, input.projected.domain_coverage,
            input.projected.sphere,
            input.projected.surface_path_length + input.warped.path_length};
  }
};

struct ColorStage : ExactStage<Pullback::StageKind::COLOR,
                               Pullback::MaterialSample, Color4, true> {
  static Color4 run(const Pullback::MaterialSample &input,
                    const FrameState &frame) {
    record(frame);
    return Color4(Pixel(1, 2, 3), input.coverage);
  }
};

using TestPipeline =
    Pullback::Pipeline<TestBinding, OuterStage, SurfaceStage, WarpStage,
                       SourceStage, MaterialStage, ColorStage>;

struct OrientationState {
  using Binding = TestBinding;
  using FrameState = typename Binding::FrameState;

  static const Quaternion &conjugate(const FrameState &) {
    static constexpr Quaternion IDENTITY;
    return IDENTITY;
  }
};

struct ProjectionPolicy : Pullback::ExactPolicy {
  static const Quaternion &frame_conjugate(const TestFrame &) {
    static constexpr Quaternion IDENTITY;
    return IDENTITY;
  }

  static Pullback::ProjectionSample project(const Vector &input,
                                            const TestFrame &) {
    return {Complex(input.x, input.y), 1, 2, 3, 0.25f, 0.5f, 4, 5, 6, 0.8f};
  }
};

struct WarpOne : Pullback::ExactPolicy {
  static Pullback::WarpStepResult apply(const Complex &input,
                                        const Pullback::ProjectionSample &,
                                        const TestFrame &) {
    return {Complex(input.re + 1.0f, input.im), Complex(1.0f, 0.0f), 1.0f,
            2.0f};
  }
};

struct WarpTwo : Pullback::ExactPolicy {
  static Pullback::WarpStepResult apply(const Complex &input,
                                        const Pullback::ProjectionSample &p,
                                        const TestFrame &) {
    return {Complex(input.re, input.im + p.coords.im),
            Complex(0.0f, p.coords.im), 2.0f, 3.0f};
  }
};

struct SourcePolicy : Pullback::ExactPolicy {
  static float sample(const Pullback::SourceInput &input, const TestFrame &) {
    return input.warped.coords.re + input.warped.coords.im;
  }
};

struct ValueState {
  using Binding = TestBinding;
  using FrameState = typename Binding::FrameState;

  static float edge_width(const FrameState &frame) { return frame.edge_width; }
};

struct MalformedValueState {
  using Binding = TestBinding;
  using FrameState = typename Binding::FrameState;
};

struct ColorPolicy : Pullback::ExactPolicy {
  static Color4 apply(const Pullback::MaterialSample &input,
                      const TestFrame &) {
    return Color4(Pixel(9, 8, 7), input.coverage);
  }
};

using CoreOuter = Pullback::Stage::OuterCamera<TestBinding, OrientationState>;
using CoreSurface = Pullback::Stage::SurfaceProject<
    TestBinding, Pullback::Surface::Identity, Pullback::Lens::Identity,
    Pullback::Surface::Identity, ProjectionPolicy>;
using CoreWarp = Pullback::Stage::PlanarWarp<TestBinding, WarpOne, WarpTwo>;
using CoreSource = Pullback::Stage::Source<TestBinding, SourcePolicy>;
using CoreMaterial =
    Pullback::Stage::Material<TestBinding, Pullback::Weight::Projection,
                              Pullback::Transfer::Linear,
                              Pullback::Coverage::EdgeFade<ValueState>>;
using CoreColor = Pullback::Stage::Color<TestBinding, ColorPolicy>;
using CorePipeline =
    Pullback::Pipeline<TestBinding, CoreOuter, CoreSurface, CoreWarp,
                       CoreSource, CoreMaterial, CoreColor>;

void test_pullback_carrier_contract() {
  constexpr Pullback::ProjectionSample projected{
      Complex(1.0f, 2.0f), 3, 4, 5, 6.0f, 0.75f, 7};
  static_assert(std::is_trivially_copyable_v<Pullback::ProjectionSample>);
  static_assert(std::is_trivially_copyable_v<Pullback::WarpResult>);
  static_assert(std::is_trivially_copyable_v<Pullback::MaterialSample>);
  HS_EXPECT_EQ(projected.traits, 0);
  HS_EXPECT_EQ(projected.edge_class, 0);
  HS_EXPECT_EQ(projected.domain_coverage, 1.0f);
  HS_EXPECT_EQ(projected.surface_path_length, 0.0f);
  HS_EXPECT_EQ(static_cast<uint8_t>(Pullback::ProjectionBoundary::CUT), 1);
  HS_EXPECT_EQ(static_cast<uint8_t>(Pullback::ProjectionBoundary::SINGULAR), 2);
}

void test_pullback_validation_predicates() {
  using Valid = typename TestPipeline::Validation;
  using Short = Pullback::PipelineValidation<TestBinding, OuterStage>;
  struct MissingContract {};
  using Missing =
      Pullback::PipelineValidation<TestBinding, MissingContract, SurfaceStage,
                                   WarpStage, SourceStage, MaterialStage,
                                   ColorStage>;
  HS_EXPECT_TRUE(Valid::ARITY);
  HS_EXPECT_TRUE(Valid::CONTRACTS);
  HS_EXPECT_TRUE(Valid::BINDINGS);
  HS_EXPECT_TRUE(Valid::EMPTY_POLICIES);
  HS_EXPECT_TRUE(Valid::ORDER);
  HS_EXPECT_TRUE(Valid::RUN_RETURNS);
  HS_EXPECT_TRUE(Valid::CARRIERS);
  HS_EXPECT_TRUE(Valid::TERMINALS);
  HS_EXPECT_TRUE(Valid::APPROXIMATIONS);
  HS_EXPECT_TRUE(Valid::EXTRA_VALIDATION);
  HS_EXPECT_FALSE(Short::ARITY);
  HS_EXPECT_FALSE(Short::ORDER);
  HS_EXPECT_TRUE(Missing::ARITY);
  HS_EXPECT_FALSE(Missing::CONTRACTS);
}

void test_pullback_evaluation_order() {
  TestFrame frame;
  const Color4 result = TestPipeline::evaluate(Vector(1.0f, 2.0f, 3.0f), frame);
  HS_EXPECT_EQ(frame.call_count, 6U);
  for (size_t index = 0; index < frame.call_count; ++index)
    HS_EXPECT_EQ(static_cast<uint8_t>(frame.calls[index]), index);
  HS_EXPECT_EQ(result.color.r, 1);
  HS_EXPECT_EQ(result.color.g, 2);
  HS_EXPECT_EQ(result.color.b, 3);
  HS_EXPECT_EQ(result.alpha, 0.5f);
}

void test_pullback_public_surface() {
  using FlashSource =
      Pullback::Stage::Placed<Pullback::CodeEmission::OUT_OF_LINE_FLASH,
                              CoreSource>;
  static_assert(TestPipeline::STAGE_COUNT == 6);
  static_assert(std::is_same_v<TestPipeline::Binding, TestBinding>);
  static_assert(std::is_same_v<TestPipeline::FrameState, TestFrame>);
  static_assert(std::is_same_v<TestPipeline::stage_at<0>, OuterStage>);
  static_assert(std::is_same_v<TestPipeline::stage_at<5>, ColorStage>);
  static_assert(FlashSource::EMISSION ==
                Pullback::CodeEmission::OUT_OF_LINE_FLASH);
  static_assert(Pullback::Warp::MAX_POLAR_HARMONIC == 16);
  static_assert(Pullback::Color::HueRotationLutView::SIZE == 1024);
  static_assert(Pullback::Color::HueNoiseLutView::SIZE == 3456);
  HS_EXPECT_TRUE(std::is_empty_v<TestPipeline>);
}

void test_pullback_no_instrumentation() {
  const Pullback::NoInstrumentation::Token token =
      Pullback::NoInstrumentation::mark();
  Pullback::NoInstrumentation::span<Pullback::ProfileEvent::COLOR>(token);
  HS_EXPECT_TRUE(std::is_empty_v<Pullback::NoInstrumentation>);
  HS_EXPECT_TRUE(std::is_empty_v<Pullback::NoInstrumentation::Token>);
}

void test_pullback_stage_combinators() {
  TestFrame frame;
  const Pullback::ProjectionSample projected =
      CoreSurface::run(Vector(1.0f, 2.0f, 3.0f), frame);
  const Pullback::SourceInput warped = CoreWarp::run(projected, frame);
  HS_EXPECT_EQ(warped.warped.coords.re, 2.0f);
  HS_EXPECT_EQ(warped.warped.coords.im, 4.0f);
  HS_EXPECT_EQ(warped.warped.net_delta.re, 1.0f);
  HS_EXPECT_EQ(warped.warped.net_delta.im, 2.0f);
  HS_EXPECT_EQ(warped.warped.path_length, 5.0f);
  HS_EXPECT_EQ(warped.projected.coords.re, 1.0f);
  HS_EXPECT_EQ(warped.projected.coords.im, 2.0f);

  const Color4 color = CorePipeline::evaluate(Vector(1.0f, 2.0f, 3.0f), frame);
  HS_EXPECT_EQ(color.color.r, 9);
  HS_EXPECT_EQ(color.color.g, 8);
  HS_EXPECT_EQ(color.color.b, 7);
  HS_EXPECT_EQ(color.alpha, 0.4f);
}

void test_pullback_provider_contracts() {
  static_assert(CoreOuter::PROVIDER_VALID);
  static_assert(
      Pullback::Coverage::EdgeFade<ValueState>::PROVIDER_VALID<TestBinding>);
  static_assert(!Pullback::Coverage::EdgeFade<
                MalformedValueState>::PROVIDER_VALID<TestBinding>);
  HS_EXPECT_TRUE(std::is_empty_v<OrientationState>);
  HS_EXPECT_TRUE(std::is_empty_v<ValueState>);
}

inline int run_pullback_tests() {
  ModuleFixture fixture("pullback");
  test_pullback_carrier_contract();
  test_pullback_validation_predicates();
  test_pullback_evaluation_order();
  test_pullback_public_surface();
  test_pullback_no_instrumentation();
  test_pullback_stage_combinators();
  test_pullback_provider_contracts();
  return fixture.result();
}

} // namespace pullback_tests
} // namespace hs_test
