/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <cstddef>
#include <type_traits>

#include "core/render/pullback.h"
#include "core/engine/fixed_pipeline.h"
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
            {Complex(input.coords.re + 2.0f, input.coords.im + 3.0f), 4.0f}};
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
    return {Complex(input.re + 1.0f, input.im), Complex(1.0f, 0.0f), 2.0f};
  }
};

struct WarpTwo : Pullback::ExactPolicy {
  static Pullback::WarpStepResult apply(const Complex &input,
                                        const Pullback::ProjectionSample &p,
                                        const TestFrame &) {
    return {Complex(input.re, input.im + p.coords.im),
            Complex(0.0f, p.coords.im), 3.0f};
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

struct MistypedMetricsOuter : OuterStage {
  static constexpr int METRICS = 0;
};

struct NonEmptyOuter : OuterStage {
  int state;
};

struct WrongOrderOuter : OuterStage {
  static constexpr Pullback::StageKind KIND =
      Pullback::StageKind::SURFACE_PROJECT;
};

struct WrongReturnOuter : OuterStage {
  static int run(const Vector &, const FrameState &) { return 0; }
};

struct WrongCarrierSurface : SurfaceStage {
  using Output = Vector;
  static Vector run(const Vector &input, const FrameState &) { return input; }
};

struct TerminalOuter : OuterStage {
  static constexpr bool TERMINAL = true;
};

struct MalformedApproximationOuter : OuterStage {
  static constexpr bool APPROXIMATE = true;
};

struct ForeignBinding {
  using FrameState = TestFrame;
  using Instrumentation = Pullback::NoInstrumentation;
};

struct ForeignOuter : OuterStage {
  using Binding = ForeignBinding;
  using FrameState = typename Binding::FrameState;
};

struct RejectBinding {
  using FrameState = TestFrame;
  using Instrumentation = Pullback::NoInstrumentation;

  template <typename...> struct ExtraValidation {
    static constexpr bool value = false;
  };
};

template <typename Stage> struct RejectStage : Stage {
  using Binding = RejectBinding;
  using FrameState = typename Binding::FrameState;
  using Input = typename Stage::Input;
  using Output = typename Stage::Output;

  static Output run(const Input &input, const FrameState &frame) {
    return Stage::run(input, frame);
  }
};

struct CountingInstrumentation {
  struct Token {};
  static inline std::array<Pullback::ProfileEvent, 8> events{};
  static inline size_t count = 0;

  static Token mark() { return {}; }

  template <Pullback::ProfileEvent Event> static void span(Token) {
    events[count++] = Event;
  }
};

struct CountingBinding {
  using FrameState = TestFrame;
  using Instrumentation = CountingInstrumentation;
};

struct CountingOrientationState {
  using Binding = CountingBinding;
  using FrameState = typename Binding::FrameState;

  static const Quaternion &conjugate(const FrameState &) {
    static constexpr Quaternion IDENTITY;
    return IDENTITY;
  }
};

struct CountingSurfacePolicy : Pullback::ExactPolicy {
  static Pullback::SurfaceResult apply(const Vector &input, const TestFrame &) {
    return {input, 0.0f};
  }
};

struct CountingLensPolicy : Pullback::ExactPolicy {
  static Vector apply(const Vector &input, const TestFrame &) { return input; }
};

struct CountingProjectionPolicy : Pullback::ExactPolicy {
  static const Quaternion &frame_conjugate(const TestFrame &) {
    static constexpr Quaternion IDENTITY;
    return IDENTITY;
  }
  static Pullback::ProjectionSample project(const Vector &input,
                                            const TestFrame &) {
    return {Complex(input.x, input.y), 0, 0, 0, 1.0f, 1.0f, 0};
  }
};

struct CountingMirrorPolicy : Pullback::ExactPolicy {
  static Pullback::WarpStepResult apply(const Complex &input,
                                        const Pullback::ProjectionSample &,
                                        const TestFrame &) {
    const auto start = CountingInstrumentation::mark();
    CountingInstrumentation::span<Pullback::ProfileEvent::MIRROR_TILE>(start);
    return {input, Complex(), 0.0f};
  }
};

struct CountingSourcePolicy : Pullback::ExactPolicy {
  static float sample(const Pullback::SourceInput &, const TestFrame &) {
    return 0.0f;
  }
};

struct CountingColorPolicy : Pullback::ExactPolicy {
  static Color4 apply(const Pullback::MaterialSample &, const TestFrame &) {
    return {};
  }
};

using CountingPipeline = Pullback::Pipeline<
    CountingBinding,
    Pullback::Stage::OuterCamera<CountingBinding, CountingOrientationState>,
    Pullback::Stage::SurfaceProject<
        CountingBinding, Pullback::Surface::Identity, CountingLensPolicy,
        CountingSurfacePolicy, CountingProjectionPolicy>,
    Pullback::Stage::PlanarWarp<CountingBinding, CountingMirrorPolicy>,
    Pullback::Stage::Source<CountingBinding, CountingSourcePolicy>,
    Pullback::Stage::Material<CountingBinding, Pullback::Weight::Projection,
                              Pullback::Transfer::Linear,
                              Pullback::Coverage::Opaque>,
    Pullback::Stage::Color<CountingBinding, CountingColorPolicy>>;

void test_pullback_carrier_contract() {
  constexpr Pullback::ProjectionSample projected{
      Complex(1.0f, 2.0f), 3, 4, 5, 6.0f, 0.75f, 7};
  static_assert(std::is_trivially_copyable_v<Pullback::ProjectionSample>);
  static_assert(std::is_trivially_copyable_v<Pullback::WarpResult>);
  static_assert(std::is_trivially_copyable_v<Pullback::MaterialSample>);
  static_assert(sizeof(Pullback::ProjectionSample) == 44);
  static_assert(sizeof(Pullback::SurfaceResult) == 16);
  static_assert(sizeof(Pullback::WarpStepResult) == 20);
  static_assert(sizeof(Pullback::WarpResult) == 12);
  static_assert(sizeof(Pullback::SourceInput) == 56);
  static_assert(sizeof(Pullback::MaterialInput) == 60);
  static_assert(sizeof(Pullback::MaterialSample) == 24);
  static_assert(alignof(Pullback::ProjectionSample) == 4);
  static_assert(alignof(Pullback::SurfaceResult) == 4);
  static_assert(alignof(Pullback::WarpStepResult) == 4);
  static_assert(alignof(Pullback::WarpResult) == 4);
  static_assert(alignof(Pullback::SourceInput) == 4);
  static_assert(alignof(Pullback::MaterialInput) == 4);
  static_assert(alignof(Pullback::MaterialSample) == 4);
  static_assert(offsetof(Pullback::ProjectionSample, coords) == 0);
  static_assert(offsetof(Pullback::ProjectionSample, fade_edge_distance) == 12);
  static_assert(offsetof(Pullback::ProjectionSample, value_weight) == 16);
  static_assert(offsetof(Pullback::ProjectionSample, domain_coverage) == 24);
  static_assert(offsetof(Pullback::ProjectionSample, sphere) == 28);
  static_assert(offsetof(Pullback::ProjectionSample, surface_path_length) ==
                40);
  static_assert(offsetof(Pullback::WarpStepResult, delta) == 8);
  static_assert(offsetof(Pullback::WarpStepResult, path_length) == 16);
  static_assert(offsetof(Pullback::WarpResult, path_length) == 8);
  static_assert(offsetof(Pullback::SourceInput, warped) == 44);
  static_assert(offsetof(Pullback::MaterialInput, field) == 56);
  static_assert(offsetof(Pullback::MaterialSample, sphere) == 8);
  static_assert(offsetof(Pullback::MaterialSample, path_length) == 20);
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
  using Mistyped =
      Pullback::PipelineValidation<TestBinding, MistypedMetricsOuter,
                                   SurfaceStage, WarpStage, SourceStage,
                                   MaterialStage, ColorStage>;
  using Foreign =
      Pullback::PipelineValidation<TestBinding, ForeignOuter, SurfaceStage,
                                   WarpStage, SourceStage, MaterialStage,
                                   ColorStage>;
  using NonEmpty =
      Pullback::PipelineValidation<TestBinding, NonEmptyOuter, SurfaceStage,
                                   WarpStage, SourceStage, MaterialStage,
                                   ColorStage>;
  using WrongOrder =
      Pullback::PipelineValidation<TestBinding, WrongOrderOuter, SurfaceStage,
                                   WarpStage, SourceStage, MaterialStage,
                                   ColorStage>;
  using WrongReturn =
      Pullback::PipelineValidation<TestBinding, WrongReturnOuter, SurfaceStage,
                                   WarpStage, SourceStage, MaterialStage,
                                   ColorStage>;
  using WrongCarrier =
      Pullback::PipelineValidation<TestBinding, OuterStage, WrongCarrierSurface,
                                   WarpStage, SourceStage, MaterialStage,
                                   ColorStage>;
  using WrongTerminal =
      Pullback::PipelineValidation<TestBinding, TerminalOuter, SurfaceStage,
                                   WarpStage, SourceStage, MaterialStage,
                                   ColorStage>;
  using WrongApproximation =
      Pullback::PipelineValidation<TestBinding, MalformedApproximationOuter,
                                   SurfaceStage, WarpStage, SourceStage,
                                   MaterialStage, ColorStage>;
  using Rejected = Pullback::PipelineValidation<
      RejectBinding, RejectStage<OuterStage>, RejectStage<SurfaceStage>,
      RejectStage<WarpStage>, RejectStage<SourceStage>,
      RejectStage<MaterialStage>, RejectStage<ColorStage>>;
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
  HS_EXPECT_FALSE(Mistyped::CONTRACTS);
  HS_EXPECT_TRUE(Foreign::CONTRACTS);
  HS_EXPECT_FALSE(Foreign::BINDINGS);
  HS_EXPECT_FALSE(NonEmpty::EMPTY_POLICIES);
  HS_EXPECT_FALSE(WrongOrder::ORDER);
  HS_EXPECT_FALSE(WrongReturn::RUN_RETURNS);
  HS_EXPECT_FALSE(WrongCarrier::CARRIERS);
  HS_EXPECT_FALSE(WrongTerminal::TERMINALS);
  HS_EXPECT_FALSE(WrongApproximation::APPROXIMATIONS);
  HS_EXPECT_FALSE(Rejected::EXTRA_VALIDATION);
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

void test_pullback_counting_instrumentation() {
  CountingInstrumentation::count = 0;
  TestFrame frame;
  static_cast<void>(
      CountingPipeline::evaluate(Vector(1.0f, 2.0f, 3.0f), frame));
  constexpr std::array EXPECTED{Pullback::ProfileEvent::LENS,
                                Pullback::ProfileEvent::SURFACE_NOISE,
                                Pullback::ProfileEvent::PROJECTION,
                                Pullback::ProfileEvent::MIRROR_TILE,
                                Pullback::ProfileEvent::PLANAR_WARP,
                                Pullback::ProfileEvent::SOURCE,
                                Pullback::ProfileEvent::MATERIAL,
                                Pullback::ProfileEvent::COLOR};
  HS_EXPECT_EQ(CountingInstrumentation::count, EXPECTED.size());
  for (size_t index = 0; index < EXPECTED.size(); ++index)
    HS_EXPECT_EQ(CountingInstrumentation::events[index], EXPECTED[index]);
}

void test_pullback_stage_combinators() {
  TestFrame frame;
  const Pullback::ProjectionSample projected =
      CoreSurface::run(Vector(1.0f, 2.0f, 3.0f), frame);
  const Pullback::SourceInput warped = CoreWarp::run(projected, frame);
  HS_EXPECT_EQ(warped.warped.coords.re, 2.0f);
  HS_EXPECT_EQ(warped.warped.coords.im, 4.0f);
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

void test_pullback_concrete_catalog() {
  struct Prepared {
    float primary;
    float secondary;
    float angle;
    float angle_cos;
    float angle_sin;
  };
  struct Params {
    float pattern_mix;
    float complexity;
    float lattice_cell_scale;
    float lattice_shape_blend;
    float lattice_softness;
    float lattice_radius;
  };
  constexpr Prepared prepared{0.0f, 0.0f, 0.0f, 1.0f, 0.0f};
  constexpr Params params{1.0f, 0.0f, 1.0f, 0.0f, 0.1f, 0.25f};
  const Complex origin;
  HS_EXPECT_EQ(Pullback::Source::twin_wave(origin, prepared), 0.0f);
  HS_EXPECT_EQ(Pullback::Source::rings(origin, prepared), 0.0f);
  HS_EXPECT_EQ(Pullback::Source::spiral(origin, prepared), 1.0f);
  HS_EXPECT_EQ(Pullback::Source::grid(origin, params, prepared), 0.0f);
  HS_EXPECT_EQ(Pullback::Source::primitive_lattice(origin, params), 1.0f);

  const Vector axis =
      Pullback::Lens::Glitch::apply(Vector(0.0f, 1.0f, 0.0f), TestFrame{});
  HS_EXPECT_EQ(axis.x, 0.0f);
  HS_EXPECT_EQ(axis.y, 1.0f);
  HS_EXPECT_EQ(axis.z, 0.0f);

  const projections::ProjectionKernelResult kernel{
      Complex(1.0f, 2.0f), 3, 4, 5, 3.0f, 6, 7, 8};
  const Pullback::ProjectionSample projected =
      Pullback::Projection::from_kernel(kernel, 2.0f);
  HS_EXPECT_EQ(projected.coords.re, 2.0f);
  HS_EXPECT_EQ(projected.coords.im, 4.0f);
  HS_EXPECT_EQ(projected.fade_edge_distance, 6.0f);
  HS_EXPECT_EQ(projected.region_id, 3);
  HS_EXPECT_EQ(projected.component_id, 4);
  HS_EXPECT_EQ(projected.boundary_flags, 5);
  HS_EXPECT_EQ(projected.flags, 6);
  HS_EXPECT_EQ(projected.traits, 7);
  HS_EXPECT_EQ(projected.edge_class, 8);
}

struct AddLens : Pullback::ExactPolicy {
  static Vector apply(const Vector &input, const TestFrame &) {
    return Vector(input.x + 1.0f, input.y, input.z);
  }
};

struct ScaleLens : Pullback::ExactPolicy {
  static Vector apply(const Vector &input, const TestFrame &) {
    return Vector(input.x * 2.0f, input.y, input.z);
  }
};

void test_pullback_lens_sequence() {
  using Sequence = Pullback::Lens::Sequence<AddLens, ScaleLens>;
  static_assert(std::is_empty_v<Sequence>);
  static_assert(!Sequence::APPROXIMATE);
  const Vector result = Sequence::apply(Vector(3.0f, 2.0f, 1.0f), TestFrame{});
  HS_EXPECT_EQ(result.x, 8.0f);
  HS_EXPECT_EQ(result.y, 2.0f);
  HS_EXPECT_EQ(result.z, 1.0f);
}

void test_fixed_pipeline_interpolation_contract() {
  HS_EXPECT_EQ(FixedPipeline::linear(2.0f, 6.0f, 0.0f), 2.0f);
  HS_EXPECT_EQ(FixedPipeline::linear(2.0f, 6.0f, 1.0f), 6.0f);
  HS_EXPECT_NEAR(FixedPipeline::log_positive(1.0f, 4.0f, 0.5f), 2.0f, 1e-6f);
  HS_EXPECT_NEAR(FixedPipeline::shortest_periodic(0.0f, 0.5f, 0.5f, 1.0f),
                 0.75f, 1e-6f);

  const auto normalized = FixedPipeline::normalized_linear<2>(
      {1.0f, 0.0f}, {0.0f, 1.0f}, 0.5f, 1e-6f);
  HS_EXPECT_TRUE(normalized.valid);
  HS_EXPECT_NEAR(normalized.values[0], 0.70710678f, 1e-6f);
  HS_EXPECT_NEAR(normalized.values[1], 0.70710678f, 1e-6f);
  const auto antipodal = FixedPipeline::normalized_linear<2>(
      {1.0f, 0.0f}, {-1.0f, 0.0f}, 0.5f, 1e-6f);
  HS_EXPECT_FALSE(antipodal.valid);
}

void test_fixed_pipeline_progress_contract() {
  const auto start = FixedPipeline::edge_progress(
      0, 1, FixedPipeline::Easing::EASE_IN_OUT_SIN);
  const auto finish = FixedPipeline::edge_progress(
      1, 1, FixedPipeline::Easing::EASE_IN_OUT_SIN);
  HS_EXPECT_EQ(start.raw, 0.0f);
  HS_EXPECT_EQ(start.eased, 0.0f);
  HS_EXPECT_EQ(finish.raw, 1.0f);
  HS_EXPECT_EQ(finish.eased, 1.0f);
  HS_EXPECT_EQ(FixedPipeline::staggered_group_progress(0.5f, 0, 2), 1.0f);
  HS_EXPECT_EQ(FixedPipeline::staggered_group_progress(0.5f, 1, 2), 0.0f);
}

inline int run_pullback_tests() {
  ModuleFixture fixture("pullback");
  test_pullback_carrier_contract();
  test_pullback_validation_predicates();
  test_pullback_evaluation_order();
  test_pullback_public_surface();
  test_pullback_no_instrumentation();
  test_pullback_counting_instrumentation();
  test_pullback_stage_combinators();
  test_pullback_provider_contracts();
  test_pullback_concrete_catalog();
  test_pullback_lens_sequence();
  test_fixed_pipeline_interpolation_contract();
  test_fixed_pipeline_progress_contract();
  return fixture.result();
}

} // namespace pullback_tests
} // namespace hs_test
