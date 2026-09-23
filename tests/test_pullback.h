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
  mutable std::array<uint8_t, 8> calls{};
  mutable size_t call_count = 0;
  float edge_width = 0.5f;
  float cutout_threshold = 0.25f;
  float cutout_softness = 0.0f;
};

struct TestBinding {
  using FrameState = TestFrame;
  using Instrumentation = Pullback::NoInstrumentation;
};

inline void record(const TestFrame &frame, uint8_t id) {
  frame.calls[frame.call_count++] = id;
}

struct EntryStage
    : Pullback::Stage::Contract<EntryStage, Pullback::SphereSample,
                                Pullback::SphereSample> {
  using Policies = std::tuple<>;

  template <typename Binding>
  static Pullback::SphereSample run(const Pullback::SphereSample &input,
                                    const TestFrame &frame,
                                    const Pullback::NoPrepared &) {
    record(frame, 0);
    return {math::Vector(input.dir.x + 1.0f, input.dir.y, input.dir.z),
            input.path_length + 1.0f};
  }
};

struct CrossingStage
    : Pullback::Stage::Contract<CrossingStage, Pullback::SphereSample,
                                Pullback::PlaneSample> {
  using Policies = std::tuple<>;

  template <typename Binding>
  static Pullback::PlaneSample run(const Pullback::SphereSample &input,
                                   const TestFrame &frame,
                                   const Pullback::NoPrepared &) {
    record(frame, 1);
    return {math::Complex(input.dir.x, input.dir.y),
            {1, 2, 3, 0.25f, 0.5f, 4, 5, 6, 0.8f},
            input.dir,
            input.path_length};
  }
};

struct PlaneStage : Pullback::Stage::Contract<PlaneStage, Pullback::PlaneSample,
                                              Pullback::PlaneSample> {
  using Policies = std::tuple<>;

  template <typename Binding>
  static Pullback::PlaneSample run(const Pullback::PlaneSample &input,
                                   const TestFrame &frame,
                                   const Pullback::NoPrepared &) {
    record(frame, 2);
    Pullback::PlaneSample output = input;
    output.coords.re += 2.0f;
    output.path_length += 2.0f;
    return output;
  }
};

struct FieldCrossingStage
    : Pullback::Stage::Contract<FieldCrossingStage, Pullback::PlaneSample,
                                Pullback::FieldSample> {
  using Policies = std::tuple<>;

  template <typename Binding>
  static Pullback::FieldSample run(const Pullback::PlaneSample &input,
                                   const TestFrame &frame,
                                   const Pullback::NoPrepared &) {
    record(frame, 3);
    return {0.5f, input.provenance.domain_coverage, input.sphere,
            input.path_length};
  }
};

struct FieldStage : Pullback::Stage::Contract<FieldStage, Pullback::FieldSample,
                                              Pullback::FieldSample> {
  using Policies = std::tuple<>;

  template <typename Binding>
  static Pullback::FieldSample run(const Pullback::FieldSample &input,
                                   const TestFrame &frame,
                                   const Pullback::NoPrepared &) {
    record(frame, 4);
    Pullback::FieldSample output = input;
    output.value *= 0.5f;
    return output;
  }
};

struct ColorCrossingStage
    : Pullback::Stage::Contract<ColorCrossingStage, Pullback::FieldSample,
                                Color4> {
  using Policies = std::tuple<>;

  template <typename Binding>
  static Color4 run(const Pullback::FieldSample &input, const TestFrame &frame,
                    const Pullback::NoPrepared &) {
    record(frame, 5);
    return Color4(Pixel(1, 2, 3), input.coverage);
  }
};

using TestPipeline =
    Pullback::Pipeline<TestBinding, EntryStage, CrossingStage, PlaneStage,
                       FieldCrossingStage, FieldStage, ColorCrossingStage>;

struct MissingContract {};

struct NonTuplePoliciesStage {
  using Input = Pullback::SphereSample;
  using Output = Pullback::SphereSample;
  using Policies = int;

  template <typename> struct Bind {};
};

struct NotACarrier {};

struct ForeignCarrierStage
    : Pullback::Stage::Contract<ForeignCarrierStage, NotACarrier, Color4> {
  using Policies = std::tuple<>;

  template <typename Binding>
  static Color4 run(const NotACarrier &, const TestFrame &,
                    const Pullback::NoPrepared &) {
    return {};
  }
};

struct DownCrossingStage
    : Pullback::Stage::Contract<DownCrossingStage, Pullback::FieldSample,
                                Pullback::PlaneSample> {
  using Policies = std::tuple<>;

  template <typename Binding>
  static Pullback::PlaneSample run(const Pullback::FieldSample &,
                                   const TestFrame &,
                                   const Pullback::NoPrepared &) {
    return {};
  }
};

struct StatefulStage : EntryStage {
  int state;
};

struct ShadowedRunStage : EntryStage {
  template <typename Binding>
  static Pullback::SphereSample run(const Pullback::SphereSample &,
                                    const TestFrame &,
                                    const Pullback::NoPrepared &) {
    return {};
  }
};

struct MissingDescriptorStage {
  using Input = Pullback::SphereSample;
  using Output = Pullback::SphereSample;
  using Policies = std::tuple<>;
  template <typename> struct Bind {};
};

struct ForeignBinding {
  using FrameState = TestFrame;
  using Instrumentation = Pullback::NoInstrumentation;
};

struct ForeignValueState {
  using Binding = ForeignBinding;
  using FrameState = TestFrame;

  static float cutout_threshold(const FrameState &frame) {
    return frame.cutout_threshold;
  }
  static float cutout_softness(const FrameState &frame) {
    return frame.cutout_softness;
  }
};

struct ValueState {
  using Binding = TestBinding;
  using FrameState = TestFrame;

  static float edge_width(const FrameState &frame) { return frame.edge_width; }
  static float cutout_threshold(const FrameState &frame) {
    return frame.cutout_threshold;
  }
  static float cutout_softness(const FrameState &frame) {
    return frame.cutout_softness;
  }
};

struct MalformedValueState {
  using Binding = TestBinding;
  using FrameState = TestFrame;
};

struct WrongReturnStage {
  using Input = Pullback::FieldSample;
  using Output = Color4;
  using Policies = std::tuple<>;

  template <typename Binding> struct Bind : Pullback::ApproximationDefaults {
    using Descriptor = WrongReturnStage;
    using Input = Pullback::FieldSample;
    using Output = Color4;
    using FrameState = typename Binding::FrameState;
    using Prepared = Pullback::NoPrepared;

    static Prepared prepare(const FrameState &) { return {}; }
    static int run(const Input &, const FrameState &, const Prepared &) {
      return 0;
    }
  };
};

struct WrongPrepareStage {
  using Input = Pullback::FieldSample;
  using Output = Color4;
  using Policies = std::tuple<>;

  template <typename Binding> struct Bind : Pullback::ApproximationDefaults {
    using Descriptor = WrongPrepareStage;
    using Input = Pullback::FieldSample;
    using Output = Color4;
    using FrameState = typename Binding::FrameState;
    using Prepared = Pullback::NoPrepared;

    static int prepare(const FrameState &) { return 0; }
    static Output run(const Input &, const FrameState &, const Prepared &) {
      return {};
    }
  };
};

struct MalformedApproximatePolicy : Pullback::ApproximationDefaults {
  static constexpr bool APPROXIMATE = true;
};

struct MalformedApproximateStage
    : Pullback::Stage::Contract<MalformedApproximateStage,
                                Pullback::FieldSample, Color4> {
  using Policies = std::tuple<MalformedApproximatePolicy>;

  template <typename Binding>
  static Color4 run(const Pullback::FieldSample &, const TestFrame &,
                    const Pullback::NoPrepared &) {
    return {};
  }
};

struct RejectBinding {
  using FrameState = TestFrame;
  using Instrumentation = Pullback::NoInstrumentation;

  template <typename...> struct ExtraValidation {
    static constexpr bool value = false;
  };
};

inline void test_pullback_carrier_contract() {
  static_assert(std::is_trivially_copyable_v<Pullback::SphereSample>);
  static_assert(std::is_trivially_copyable_v<Pullback::PlaneSample>);
  static_assert(std::is_trivially_copyable_v<Pullback::FieldSample>);
  static_assert(std::is_trivially_destructible_v<Pullback::SphereSample>);
  static_assert(std::is_trivially_destructible_v<Pullback::PlaneSample>);
  static_assert(std::is_trivially_destructible_v<Pullback::FieldSample>);
  static_assert(std::is_trivially_destructible_v<Color4>);
  static_assert(sizeof(Pullback::SphereSample) == 16);
  static_assert(sizeof(Pullback::ProjectionProvenance) == 20);
  static_assert(sizeof(Pullback::ProjectionResult) == 28);
  static_assert(sizeof(Pullback::PlaneSample) == 44);
  static_assert(sizeof(Pullback::FieldSample) == 24);
  static_assert(sizeof(Pullback::SurfaceResult) == 16);
  static_assert(sizeof(Pullback::WarpStepResult) == 12);
  static_assert(alignof(Pullback::SphereSample) == 4);
  static_assert(alignof(Pullback::PlaneSample) == 4);
  static_assert(alignof(Pullback::FieldSample) == 4);
  static_assert(offsetof(Pullback::SphereSample, path_length) == 12);
  static_assert(offsetof(Pullback::ProjectionProvenance, fade_edge_distance) ==
                4);
  static_assert(offsetof(Pullback::ProjectionProvenance, value_weight) == 8);
  static_assert(offsetof(Pullback::ProjectionProvenance, domain_coverage) ==
                16);
  static_assert(offsetof(Pullback::PlaneSample, provenance) == 8);
  static_assert(offsetof(Pullback::PlaneSample, sphere) == 28);
  static_assert(offsetof(Pullback::PlaneSample, path_length) == 40);
  static_assert(offsetof(Pullback::FieldSample, sphere) == 8);
  static_assert(offsetof(Pullback::FieldSample, path_length) == 20);

  static_assert(Pullback::family_of<Pullback::SphereSample> == 0);
  static_assert(Pullback::family_of<Pullback::PlaneSample> == 1);
  static_assert(Pullback::family_of<Pullback::FieldSample> == 2);
  static_assert(Pullback::family_of<Color4> == 3);
  static_assert(Pullback::family_of<NotACarrier> ==
                Pullback::FOREIGN_FAMILY_RANK);
  static_assert(Pullback::CanonicalCarrier<Pullback::SphereSample>);
  static_assert(Pullback::CanonicalCarrier<Color4>);
  static_assert(!Pullback::CanonicalCarrier<NotACarrier>);
  static_assert(!Pullback::CanonicalCarrier<Pullback::ProjectionResult>);

  constexpr Pullback::ProjectionProvenance provenance{1, 2, 3, 0.25f, 0.5f, 4};
  HS_EXPECT_EQ(provenance.traits, 0);
  HS_EXPECT_EQ(provenance.edge_class, 0);
  HS_EXPECT_EQ(provenance.domain_coverage, 1.0f);
  HS_EXPECT_EQ(static_cast<uint8_t>(Pullback::ProjectionBoundary::CUT), 1);
  HS_EXPECT_EQ(static_cast<uint8_t>(Pullback::ProjectionBoundary::SINGULAR), 2);
}

inline void test_pullback_validation_predicates() {
  using Valid = typename TestPipeline::Validation;
  HS_EXPECT_TRUE(Valid::NONEMPTY);
  HS_EXPECT_TRUE(Valid::CONTRACTS);
  HS_EXPECT_TRUE(Valid::CANONICAL);
  HS_EXPECT_TRUE(Valid::MONOTONE);
  HS_EXPECT_TRUE(Valid::CARRIERS);
  HS_EXPECT_TRUE(Valid::ENTRY);
  HS_EXPECT_TRUE(Valid::EXIT);
  HS_EXPECT_TRUE(Valid::BINDINGS);
  HS_EXPECT_TRUE(Valid::DESCRIPTOR_IDENTITY);
  HS_EXPECT_TRUE(Valid::EMPTY_DESCRIPTORS);
  HS_EXPECT_TRUE(Valid::RUN_RETURNS);
  HS_EXPECT_TRUE(Valid::PREPARES);
  HS_EXPECT_TRUE(Valid::APPROXIMATIONS);
  HS_EXPECT_TRUE(Valid::EXTRA_VALIDATION);

  using Empty = Pullback::PipelineValidation<TestBinding, void>;
  HS_EXPECT_FALSE(Empty::NONEMPTY);

  using Missing = Pullback::PipelineValidation<TestBinding, MissingContract,
                                               ColorCrossingStage>;
  HS_EXPECT_TRUE(Missing::NONEMPTY);
  HS_EXPECT_FALSE(Missing::CONTRACTS);
  HS_EXPECT_FALSE(Missing::DESCRIPTOR_IDENTITY);

  using Shadowed = Pullback::PipelineValidation<TestBinding, ShadowedRunStage,
                                                ColorCrossingStage>;
  HS_EXPECT_TRUE(Shadowed::BINDINGS);
  HS_EXPECT_FALSE(Shadowed::DESCRIPTOR_IDENTITY);
  using MissingDescriptor =
      Pullback::PipelineValidation<TestBinding, MissingDescriptorStage,
                                   ColorCrossingStage>;
  HS_EXPECT_TRUE(MissingDescriptor::CONTRACTS);
  HS_EXPECT_FALSE(MissingDescriptor::DESCRIPTOR_IDENTITY);
  HS_EXPECT_FALSE(MissingDescriptor::RUN_RETURNS);

  using NonTuple =
      Pullback::PipelineValidation<TestBinding, NonTuplePoliciesStage,
                                   ColorCrossingStage>;
  HS_EXPECT_TRUE(NonTuple::NONEMPTY);
  HS_EXPECT_FALSE(NonTuple::CONTRACTS);

  using Foreign = Pullback::PipelineValidation<TestBinding, EntryStage,
                                               ForeignCarrierStage>;
  HS_EXPECT_TRUE(Foreign::CONTRACTS);
  HS_EXPECT_FALSE(Foreign::CANONICAL);

  using Down =
      Pullback::PipelineValidation<TestBinding, EntryStage, CrossingStage,
                                   FieldCrossingStage, DownCrossingStage,
                                   FieldCrossingStage, ColorCrossingStage>;
  HS_EXPECT_TRUE(Down::CANONICAL);
  HS_EXPECT_FALSE(Down::MONOTONE);

  using Mismatch =
      Pullback::PipelineValidation<TestBinding, EntryStage, FieldCrossingStage,
                                   FieldStage, ColorCrossingStage>;
  HS_EXPECT_FALSE(Mismatch::CARRIERS);
  HS_EXPECT_TRUE(Mismatch::MONOTONE);

  using WrongEntry =
      Pullback::PipelineValidation<TestBinding, PlaneStage, FieldCrossingStage,
                                   ColorCrossingStage>;
  HS_EXPECT_FALSE(WrongEntry::ENTRY);
  HS_EXPECT_TRUE(WrongEntry::EXIT);

  using WrongExit =
      Pullback::PipelineValidation<TestBinding, EntryStage, CrossingStage,
                                   FieldCrossingStage, FieldStage>;
  HS_EXPECT_TRUE(WrongExit::ENTRY);
  HS_EXPECT_FALSE(WrongExit::EXIT);

  using ForeignBound = Pullback::PipelineValidation<
      TestBinding, EntryStage, CrossingStage, FieldCrossingStage,
      Pullback::Stage::ApplyCoverage<
          Pullback::ValueCoverage::ValueCutout<ForeignValueState>>,
      ColorCrossingStage>;
  HS_EXPECT_FALSE(ForeignBound::BINDINGS);
  HS_EXPECT_FALSE(ForeignBound::DESCRIPTOR_IDENTITY);

  using NonEmpty =
      Pullback::PipelineValidation<TestBinding, StatefulStage, CrossingStage,
                                   FieldCrossingStage, ColorCrossingStage>;
  HS_EXPECT_FALSE(NonEmpty::EMPTY_DESCRIPTORS);

  using WrongReturn =
      Pullback::PipelineValidation<TestBinding, EntryStage, CrossingStage,
                                   FieldCrossingStage, WrongReturnStage>;
  HS_EXPECT_TRUE(WrongReturn::BINDINGS);
  HS_EXPECT_FALSE(WrongReturn::RUN_RETURNS);

  using WrongPrepare =
      Pullback::PipelineValidation<TestBinding, EntryStage, CrossingStage,
                                   FieldCrossingStage, WrongPrepareStage>;
  HS_EXPECT_FALSE(WrongPrepare::PREPARES);

  using WrongApproximation =
      Pullback::PipelineValidation<TestBinding, EntryStage, CrossingStage,
                                   FieldCrossingStage,
                                   MalformedApproximateStage>;
  HS_EXPECT_FALSE(WrongApproximation::APPROXIMATIONS);

  using Rejected =
      Pullback::PipelineValidation<RejectBinding, EntryStage, CrossingStage,
                                   FieldCrossingStage, ColorCrossingStage>;
  HS_EXPECT_FALSE(Rejected::EXTRA_VALIDATION);
}

inline void test_pullback_evaluation_order() {
  const TestPipeline::Frame frame = TestPipeline::prepare(TestFrame{});
  const Color4 result =
      TestPipeline::shade(math::Vector(1.0f, 2.0f, 3.0f), frame);
  HS_EXPECT_EQ(frame.ctx.call_count, 6U);
  for (size_t index = 0; index < frame.ctx.call_count; ++index)
    HS_EXPECT_EQ(frame.ctx.calls[index], index);
  HS_EXPECT_EQ(result.color.r, 1);
  HS_EXPECT_EQ(result.color.g, 2);
  HS_EXPECT_EQ(result.color.b, 3);
  HS_EXPECT_EQ(result.alpha, 0.8f);
}

using GroupedPipeline = Pullback::Pipeline<
    TestBinding, EntryStage,
    Pullback::Stage::Placed<
        Pullback::CodeEmission::OUT_OF_LINE_FLASH, CrossingStage, void,
        Pullback::Stage::Placed<Pullback::CodeEmission::INLINE_ONLY,
                                PlaneStage>>,
    void, FieldCrossingStage,
    Pullback::Stage::Placed<Pullback::CodeEmission::INLINE_ONLY>, FieldStage,
    ColorCrossingStage>;

inline void test_pullback_placement_transparency() {
  // Placement wrappers are invisible to the semantic leaf view: identical
  // leaves, validation and output; only the execution-node view differs.
  static_assert(TestPipeline::STAGE_COUNT == 6);
  static_assert(GroupedPipeline::STAGE_COUNT == 6);
  static_assert(std::is_same_v<GroupedPipeline::stage_at<1>, CrossingStage>);
  static_assert(std::is_same_v<GroupedPipeline::stage_at<2>, PlaneStage>);
  static_assert(
      std::is_same_v<TestPipeline::stage_at<5>, GroupedPipeline::stage_at<5>>);
  static_assert(TestPipeline::NODE_COUNT == 6);
  static_assert(GroupedPipeline::NODE_COUNT == 5);
  static_assert(GroupedPipeline::Validation::MONOTONE);
  static_assert(GroupedPipeline::Validation::CARRIERS);

  const TestPipeline::Frame flat = TestPipeline::prepare(TestFrame{});
  const GroupedPipeline::Frame grouped = GroupedPipeline::prepare(TestFrame{});
  const Color4 flat_result =
      TestPipeline::shade(math::Vector(1.0f, 2.0f, 3.0f), flat);
  const Color4 grouped_result =
      GroupedPipeline::shade(math::Vector(1.0f, 2.0f, 3.0f), grouped);
  HS_EXPECT_EQ(grouped.ctx.call_count, 6U);
  for (size_t index = 0; index < grouped.ctx.call_count; ++index)
    HS_EXPECT_EQ(grouped.ctx.calls[index], index);
  HS_EXPECT_EQ(flat_result.alpha, grouped_result.alpha);
  HS_EXPECT_EQ(flat_result.color.r, grouped_result.color.r);
}

template <typename Stage> struct IsCrossingStage {
  static constexpr bool value = std::is_same_v<Stage, CrossingStage>;
};

template <typename Stage> struct IsColorStage {
  static constexpr bool value = std::is_same_v<Stage, ColorCrossingStage>;
};

template <typename Stage> struct NoStageMatches {
  static constexpr bool value = false;
};

inline void test_pullback_public_surface() {
  static_assert(std::is_same_v<TestPipeline::Binding, TestBinding>);
  static_assert(std::is_same_v<TestPipeline::FrameState, TestFrame>);
  static_assert(std::is_same_v<TestPipeline::stage_at<0>, EntryStage>);
  static_assert(std::is_same_v<TestPipeline::stage_at<5>, ColorCrossingStage>);
  static_assert(TestPipeline::any_stage<IsCrossingStage>);
  static_assert(!TestPipeline::any_stage<NoStageMatches>);
  static_assert(std::is_same_v<TestPipeline::stage_matching<IsColorStage>,
                               ColorCrossingStage>);
  static_assert(
      std::is_same_v<TestPipeline::stage_matching<NoStageMatches>, void>);
  static_assert(
      Pullback::Stage::Placed<Pullback::CodeEmission::OUT_OF_LINE_FLASH,
                              CrossingStage>::EMISSION ==
      Pullback::CodeEmission::OUT_OF_LINE_FLASH);
  static_assert(!Pullback::StageMatchesKey<TestPipeline, int>,
                "implements() must stay undetectable on undecorated stages");
  static_assert(Pullback::StageDescriptor<EntryStage>);
  static_assert(!Pullback::StageDescriptor<MissingContract>);
  static_assert(!Pullback::StageDescriptor<NonTuplePoliciesStage>);
  static_assert(Pullback::Warp::MAX_POLAR_HARMONIC == 16);
  static_assert(Pullback::Color::HueRotationLutView::SIZE == 1024);
  static_assert(Pullback::Color::HueNoiseLutView::SIZE == 3456);
  HS_EXPECT_TRUE(std::is_empty_v<TestPipeline>);
}

inline void test_pullback_no_instrumentation() {
  const Pullback::NoInstrumentation::Token token =
      Pullback::NoInstrumentation::mark();
  Pullback::NoInstrumentation::span<Pullback::ProfileEvent::COLOR>(token);
  HS_EXPECT_TRUE(std::is_empty_v<Pullback::NoInstrumentation>);
  HS_EXPECT_TRUE(std::is_empty_v<Pullback::NoInstrumentation::Token>);
}

struct CountingInstrumentation {
  struct Token {};
  static inline std::array<Pullback::ProfileEvent, 24> events{};
  static inline size_t count = 0;

  static Token mark() { return {}; }

  template <Pullback::ProfileEvent Event> static void span(Token) {
    HS_CHECK(count < events.size(), "CountingInstrumentation overflow");
    events[count++] = Event;
  }
};

struct CountingBinding {
  using FrameState = TestFrame;
  using Instrumentation = CountingInstrumentation;
};

struct CountingOrientationState {
  using Binding = CountingBinding;
  using FrameState = TestFrame;

  static const math::Quaternion &conjugate(const FrameState &) {
    static constexpr math::Quaternion IDENTITY;
    return IDENTITY;
  }
};

struct CountingSurfacePolicy : Pullback::ApproximationDefaults {
  static Pullback::SurfaceResult apply(const math::Vector &input,
                                       const TestFrame &) {
    return {input, 0.5f};
  }
};

struct CountingLensPolicy : Pullback::ApproximationDefaults {
  static math::Vector apply(const math::Vector &input, const TestFrame &) {
    return input;
  }
};

struct CountingProjectionPolicy : Pullback::ApproximationDefaults {
  static const math::Quaternion &frame_conjugate(const TestFrame &) {
    static constexpr math::Quaternion IDENTITY;
    return IDENTITY;
  }
  static Pullback::ProjectionResult project(const math::Vector &input,
                                            const TestFrame &) {
    return {math::Complex(input.x, input.y), {0, 0, 0, 1.0f, 0.5f, 0}};
  }
};

struct CountingWarpPolicy : Pullback::ApproximationDefaults {
  static Pullback::WarpStepResult apply(const math::Complex &input,
                                        const Pullback::ProjectionProvenance &,
                                        const TestFrame &) {
    return {math::Complex(input.re + 1.0f, input.im), 3.0f};
  }
};

struct CountingSourcePolicy : Pullback::ApproximationDefaults {
  static float sample(const Pullback::PlaneSample &input, const TestFrame &) {
    return input.coords.re;
  }
};

struct CountingSphericalSourcePolicy : Pullback::ApproximationDefaults {
  static float sample(const Pullback::SphereSample &input, const TestFrame &) {
    return input.dir.x;
  }
};

struct CountingValueState {
  using Binding = CountingBinding;
  using FrameState = TestFrame;

  static float cutout_threshold(const FrameState &frame) {
    return frame.cutout_threshold;
  }
  static float cutout_softness(const FrameState &frame) {
    return frame.cutout_softness;
  }
};

struct CountingColorPolicy : Pullback::ApproximationDefaults {
  static Color4 apply(const Pullback::FieldSample &input, const TestFrame &) {
    return Color4(Pixel(9, 8, 7), input.coverage);
  }
};

struct PreparedOrientationPolicy {
  using Binding = CountingBinding;
  using FrameState = TestFrame;
  using Prepared = int;
  static Prepared prepare(const FrameState &) { return 1; }
  // Index 1, what prepare() returns, is a half turn about y.
  static const math::Quaternion &conjugate(const FrameState &,
                                           const Prepared &p) {
    static constexpr math::Quaternion ROTATIONS[] = {
        math::Quaternion(), math::Quaternion(0.0f, 0.0f, 1.0f, 0.0f)};
    return ROTATIONS[p];
  }
};

struct PreparedLensPolicy : Pullback::ApproximationDefaults {
  using Prepared = int;
  static Prepared prepare(const TestFrame &) { return 2; }
  static math::Vector apply(const math::Vector &input, const TestFrame &,
                            const Prepared &p) {
    return input + math::Vector(static_cast<float>(p), 0.0f, 0.0f);
  }
};

struct PreparedProjectionPolicy : Pullback::ApproximationDefaults {
  using Prepared = int;
  static Prepared prepare(const TestFrame &) { return 3; }
  static const math::Quaternion &frame_conjugate(const TestFrame &,
                                                 const Prepared &) {
    static constexpr math::Quaternion IDENTITY;
    return IDENTITY;
  }
  static Pullback::ProjectionResult
  project(const math::Vector &input, const TestFrame &, const Prepared &p) {
    return {math::Complex(input.x + static_cast<float>(p), input.y),
            {0, 0, 0, 1.0f, 1.0f, 0}};
  }
};

struct PreparedSourcePolicy : Pullback::ApproximationDefaults {
  using Prepared = float;
  static Prepared prepare(const TestFrame &) { return 0.5f; }
  static float sample(const Pullback::PlaneSample &input, const TestFrame &,
                      const Prepared &prepared) {
    return input.coords.re + prepared;
  }
};

struct PreparedWeightPolicy : Pullback::ApproximationDefaults {
  using Prepared = float;
  static Prepared prepare(const TestFrame &) { return 0.5f; }
  static float apply(float field, const Pullback::ProjectionProvenance &,
                     const TestFrame &, const Prepared &prepared) {
    return field * prepared;
  }
};

struct PreparedSampleCoveragePolicy : Pullback::ApproximationDefaults {
  using Prepared = float;
  static Prepared prepare(const TestFrame &) { return 0.25f; }
  static float apply(const Pullback::ProjectionProvenance &, const TestFrame &,
                     const Prepared &prepared) {
    return prepared;
  }
};

struct MissingPreparedWeightApply : Pullback::ApproximationDefaults {
  using Prepared = float;
  static Prepared prepare(const TestFrame &) { return 0.5f; }
  static float apply(float field, const Pullback::ProjectionProvenance &,
                     const TestFrame &) {
    return field;
  }
};

struct MissingPreparedSampleCoverageApply : Pullback::ApproximationDefaults {
  using Prepared = float;
  static Prepared prepare(const TestFrame &) { return 0.25f; }
  static float apply(const Pullback::ProjectionProvenance &,
                     const TestFrame &) {
    return 1.0f;
  }
};

struct PreparedTransferPolicy : Pullback::ApproximationDefaults,
                                Pullback::TransferRole {
  using Prepared = float;
  static Prepared prepare(const TestFrame &) { return 0.25f; }
  static float apply(float value, const TestFrame &, const Prepared &p) {
    return value + p;
  }
};

struct PreparedCoveragePolicy : Pullback::ApproximationDefaults,
                                Pullback::CoverageRole {
  using Prepared = float;
  static Prepared prepare(const TestFrame &) { return 0.5f; }
  static float apply(float, const TestFrame &, const Prepared &p) { return p; }
};

struct PreparedColorPolicy : Pullback::ApproximationDefaults {
  using Prepared = uint8_t;
  static Prepared prepare(const TestFrame &) { return 11; }
  static Color4 apply(const Pullback::FieldSample &input, const TestFrame &,
                      const Prepared &p) {
    return Color4(Pixel(p, 0, 0), input.coverage);
  }
};

using CountingPipeline = Pullback::Pipeline<
    CountingBinding, Pullback::Stage::Rotate<CountingOrientationState>,
    Pullback::Stage::Displace<CountingSurfacePolicy>,
    Pullback::Stage::Lens<CountingLensPolicy>,
    Pullback::Stage::Project<CountingProjectionPolicy>,
    Pullback::Stage::Warp<CountingWarpPolicy>,
    Pullback::Stage::Sample<CountingSourcePolicy>,
    Pullback::Stage::Transfer<Pullback::Transfer::Ridge>,
    Pullback::Stage::ApplyCoverage<
        Pullback::ValueCoverage::ValueCutout<CountingValueState>>,
    Pullback::Stage::Colorize<CountingColorPolicy>>;

inline void test_pullback_counting_instrumentation() {
  CountingInstrumentation::count = 0;
  static_cast<void>(CountingPipeline::shade(
      math::X_AXIS, CountingPipeline::prepare(TestFrame{})));
  constexpr std::array EXPECTED{Pullback::ProfileEvent::SURFACE_NOISE,
                                Pullback::ProfileEvent::LENS,
                                Pullback::ProfileEvent::PROJECTION,
                                Pullback::ProfileEvent::PLANAR_WARP,
                                Pullback::ProfileEvent::SOURCE,
                                Pullback::ProfileEvent::MATERIAL,
                                Pullback::ProfileEvent::MATERIAL,
                                Pullback::ProfileEvent::MATERIAL,
                                Pullback::ProfileEvent::COLOR};
  HS_EXPECT_EQ(CountingInstrumentation::count, EXPECTED.size());
  for (size_t index = 0; index < EXPECTED.size(); ++index)
    HS_EXPECT_EQ(CountingInstrumentation::events[index], EXPECTED[index]);
}

inline void test_pullback_prepared_stage_policies() {
  CountingInstrumentation::count = 0;
  const TestFrame frame;
  const Pullback::SphereSample sphere{math::Vector(1.0f, 2.0f, 3.0f), 0.0f};

  using BoundRotate =
      Pullback::Stage::Rotate<PreparedOrientationPolicy>::Bind<CountingBinding>;
  const auto rotate_prepared = BoundRotate::prepare(frame);
  HS_EXPECT_EQ(rotate_prepared, 1);
  const Pullback::SphereSample rotated =
      BoundRotate::run(sphere, frame, rotate_prepared);
  HS_EXPECT_EQ(rotated.dir.x, -1.0f);
  HS_EXPECT_EQ(rotated.dir.y, 2.0f);
  HS_EXPECT_EQ(rotated.dir.z, -3.0f);
  HS_EXPECT_EQ(rotated.path_length, sphere.path_length);

  using BoundLens =
      Pullback::Stage::Lens<PreparedLensPolicy>::Bind<CountingBinding>;
  const auto lens_prepared = BoundLens::prepare(frame);
  HS_EXPECT_EQ(lens_prepared, 2);
  HS_EXPECT_EQ(BoundLens::run(sphere, frame, lens_prepared).dir.x, 3.0f);

  using BoundProject =
      Pullback::Stage::Project<PreparedProjectionPolicy>::Bind<CountingBinding>;
  const auto project_prepared = BoundProject::prepare(frame);
  HS_EXPECT_EQ(project_prepared, 3);
  HS_EXPECT_EQ(BoundProject::run({math::X_AXIS, sphere.path_length}, frame,
                                 project_prepared)
                   .coords.re,
               4.0f);

  using BoundSourceOnlySample =
      Pullback::Stage::Sample<PreparedSourcePolicy>::Bind<CountingBinding>;
  static_assert(std::is_same_v<typename BoundSourceOnlySample::Prepared,
                               PreparedSourcePolicy::Prepared>);

  using BoundSample = Pullback::Stage::Sample<
      PreparedSourcePolicy, PreparedWeightPolicy,
      PreparedSampleCoveragePolicy>::Bind<CountingBinding>;
  const Pullback::PlaneSample plane{math::Complex(1.0f, 0.0f),
                                    {0, 0, 0, 0.0f, 1.0f, 0, 0, 0, 0.8f},
                                    math::X_AXIS,
                                    0.75f};
  const auto sample_prepared = BoundSample::prepare(frame);
  HS_EXPECT_EQ(std::get<0>(sample_prepared), 0.5f);
  HS_EXPECT_EQ(std::get<1>(sample_prepared), 0.5f);
  HS_EXPECT_EQ(std::get<2>(sample_prepared), 0.25f);
  const Pullback::FieldSample sampled =
      BoundSample::run(plane, frame, sample_prepared);
  HS_EXPECT_EQ(sampled.value, 0.875f);
  HS_EXPECT_EQ(sampled.coverage, 0.2f);
  HS_EXPECT_EQ(sampled.path_length, 0.75f);

  const Pullback::FieldSample field{0.25f, 0.8f, math::X_AXIS, 0.0f};
  using BoundTransfer =
      Pullback::Stage::Transfer<PreparedTransferPolicy>::Bind<CountingBinding>;
  const auto transfer_prepared = BoundTransfer::prepare(frame);
  HS_EXPECT_EQ(BoundTransfer::run(field, frame, transfer_prepared).value, 0.5f);

  using BoundCoverage = Pullback::Stage::ApplyCoverage<
      PreparedCoveragePolicy>::Bind<CountingBinding>;
  const auto coverage_prepared = BoundCoverage::prepare(frame);
  HS_EXPECT_EQ(BoundCoverage::run(field, frame, coverage_prepared).coverage,
               0.4f);

  using BoundColorize =
      Pullback::Stage::Colorize<PreparedColorPolicy>::Bind<CountingBinding>;
  const auto color_prepared = BoundColorize::prepare(frame);
  HS_EXPECT_EQ(BoundColorize::run(field, frame, color_prepared).color.r, 11);
}

inline void test_pullback_stage_combinators() {
  CountingInstrumentation::count = 0;
  const TestFrame frame;

  using BoundProject =
      Pullback::Stage::Project<CountingProjectionPolicy>::Bind<CountingBinding>;
  const Pullback::SphereSample view{math::X_AXIS, 0.5f};
  const Pullback::PlaneSample projected =
      BoundProject::run(view, frame, BoundProject::prepare(frame));
  HS_EXPECT_EQ(projected.coords.re, 1.0f);
  HS_EXPECT_EQ(projected.coords.im, 0.0f);
  HS_EXPECT_EQ(projected.provenance.value_weight, 0.5f);
  HS_EXPECT_EQ(projected.sphere.x, 1.0f);
  HS_EXPECT_EQ(projected.sphere.z, 0.0f);
  HS_EXPECT_EQ(projected.path_length, 0.5f);

  using BoundWarp =
      Pullback::Stage::Warp<CountingWarpPolicy>::Bind<CountingBinding>;
  const Pullback::PlaneSample warped =
      BoundWarp::run(projected, frame, BoundWarp::prepare(frame));
  HS_EXPECT_EQ(warped.coords.re, 2.0f);
  HS_EXPECT_EQ(warped.path_length, 3.5f);
  HS_EXPECT_EQ(warped.provenance.value_weight,
               projected.provenance.value_weight);
  HS_EXPECT_EQ(warped.sphere.x, projected.sphere.x);

  // Sample: weight scales the raw signed field, the ramp maps it into [0, 1],
  // and the crossing's coverage multiplies the provenance domain coverage.
  using BoundSample =
      Pullback::Stage::Sample<CountingSourcePolicy>::Bind<CountingBinding>;
  const Pullback::FieldSample sampled =
      BoundSample::run(warped, frame, BoundSample::prepare(frame));
  // raw = 2.0, weighted = 2.0 * 0.5 = 1.0, ramped = clamp((1 + 1) / 2) = 1.0
  HS_EXPECT_EQ(sampled.value, 1.0f);
  HS_EXPECT_EQ(sampled.coverage, 0.5f);
  HS_EXPECT_EQ(sampled.sphere.x, warped.sphere.x);
  HS_EXPECT_EQ(sampled.path_length, 3.5f);

  using BoundSphereSample = Pullback::Stage::SampleSphere<
      CountingSphericalSourcePolicy>::Bind<CountingBinding>;
  const Pullback::FieldSample sphere_sampled =
      BoundSphereSample::run(view, frame, BoundSphereSample::prepare(frame));
  HS_EXPECT_EQ(sphere_sampled.value, 1.0f);
  HS_EXPECT_EQ(sphere_sampled.coverage, 1.0f);
  HS_EXPECT_EQ(sphere_sampled.sphere.x, view.dir.x);
  HS_EXPECT_EQ(sphere_sampled.path_length, view.path_length);

  using BoundTransfer = Pullback::Stage::Transfer<
      Pullback::Transfer::Ridge>::Bind<CountingBinding>;
  const Pullback::FieldSample transferred =
      BoundTransfer::run(sampled, frame, BoundTransfer::prepare(frame));
  HS_EXPECT_EQ(transferred.value, unit_bell(1.0f));
  HS_EXPECT_EQ(transferred.coverage, sampled.coverage);

  using BoundCoverage =
      Pullback::Stage::ApplyCoverage<Pullback::ValueCoverage::ValueCutout<
          CountingValueState>>::Bind<CountingBinding>;
  const Pullback::FieldSample cut =
      BoundCoverage::run(sampled, frame, BoundCoverage::prepare(frame));
  // value 1.0 is above the 0.25 threshold: the hard cut keeps full coverage.
  HS_EXPECT_EQ(cut.coverage, sampled.coverage);

  using BoundDisplace =
      Pullback::Stage::Displace<CountingSurfacePolicy>::Bind<CountingBinding>;
  const Pullback::SphereSample displaced =
      BoundDisplace::run(view, frame, BoundDisplace::prepare(frame));
  HS_EXPECT_EQ(displaced.dir.x, 1.0f);
  HS_EXPECT_EQ(displaced.path_length, 1.0f);
}

inline void test_pullback_provider_contracts() {
  static_assert(Pullback::descriptor_bindable<
                Pullback::Stage::Rotate<CountingOrientationState>,
                CountingBinding>());
  static_assert(
      !Pullback::descriptor_bindable<
          Pullback::Stage::Rotate<CountingOrientationState>, TestBinding>());
  static_assert(!Pullback::descriptor_bindable<
                Pullback::Stage::Sample<CountingSourcePolicy,
                                        MissingPreparedWeightApply>,
                CountingBinding>());
  static_assert(
      !Pullback::descriptor_bindable<
          Pullback::Stage::Sample<CountingSourcePolicy, Pullback::Weight::None,
                                  MissingPreparedSampleCoverageApply>,
          CountingBinding>());
  static_assert(Pullback::ProjectionCoverage::EdgeFade<
                ValueState>::PROVIDER_VALID<TestBinding>);
  static_assert(!Pullback::ProjectionCoverage::EdgeFade<
                MalformedValueState>::PROVIDER_VALID<TestBinding>);
  static_assert(Pullback::ValueCoverage::ValueCutout<
                ValueState>::PROVIDER_VALID<TestBinding>);
  static_assert(!Pullback::ValueCoverage::ValueCutout<
                MalformedValueState>::PROVIDER_VALID<TestBinding>);
  // The two float -> float stages share a call signature; only the role tag
  // separates them.
  static_assert(Pullback::descriptor_bindable<
                Pullback::Stage::Transfer<Pullback::Transfer::Ridge>,
                CountingBinding>());
  static_assert(!Pullback::descriptor_bindable<
                Pullback::Stage::ApplyCoverage<Pullback::Transfer::Ridge>,
                CountingBinding>());
  static_assert(Pullback::descriptor_bindable<
                Pullback::Stage::ApplyCoverage<
                    Pullback::ValueCoverage::ValueCutout<CountingValueState>>,
                CountingBinding>());
  static_assert(!Pullback::descriptor_bindable<
                Pullback::Stage::Transfer<
                    Pullback::ValueCoverage::ValueCutout<CountingValueState>>,
                CountingBinding>());
  // Runtime half: a provider that declares no Prepared collapses to an empty
  // bound prepared state, and the bound stage still runs.
  const TestFrame frame;
  const Pullback::SphereSample sphere{math::Vector(1.0f, 2.0f, 3.0f), 4.0f};
  using BoundStatelessRotate =
      Pullback::Stage::Rotate<CountingOrientationState>::Bind<CountingBinding>;
  HS_EXPECT_TRUE(std::is_empty_v<typename BoundStatelessRotate::Prepared>);
  const Pullback::SphereSample carried = BoundStatelessRotate::run(
      sphere, frame, BoundStatelessRotate::prepare(frame));
  HS_EXPECT_EQ(carried.dir.x, 1.0f);
  HS_EXPECT_EQ(carried.dir.y, 2.0f);
  HS_EXPECT_EQ(carried.dir.z, 3.0f);
  HS_EXPECT_EQ(carried.path_length, 4.0f);
}

inline void test_pullback_concrete_catalog() {
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
  constexpr Prepared shifted{0.7f, 0.3f, 0.0f, 1.0f, 0.0f};
  constexpr Prepared turned{0.0f, 0.0f, math::PI_F * 0.5f, 0.0f, 1.0f};
  constexpr Params params{1.0f, 0.0f, 1.0f, 0.0f, 0.1f, 0.25f};
  constexpr Params coupled{0.0f, 0.0f, 1.0f, 0.0f, 0.1f, 0.25f};
  const math::Complex origin;
  const math::Complex off(0.6f, 0.25f); // radius 0.65, azimuth atan2(0.25, 0.6)
  HS_EXPECT_EQ(Pullback::Source::twin_wave(origin, prepared), 0.0f);
  HS_EXPECT_NEAR(Pullback::Source::twin_wave(off, turned),
                 0.5f * (math::fast_sinf(0.6f) + math::fast_sinf(0.25f)),
                 2e-3f);
  HS_EXPECT_EQ(Pullback::Source::rings(origin, prepared), 0.0f);
  HS_EXPECT_NEAR(Pullback::Source::rings(off, shifted), math::fast_sinf(-0.05f),
                 2e-3f);
  HS_EXPECT_NEAR(Pullback::Source::spiral(origin, prepared), 1.0f, 2e-3f);
  HS_EXPECT_NEAR(Pullback::Source::spiral(off, shifted),
                 math::fast_sinf(-0.05f - 3.0f * math::fast_atan2(0.25f, 0.6f)),
                 2e-3f);
  HS_EXPECT_EQ(Pullback::Source::grid(origin, params, prepared), 0.0f);
  HS_EXPECT_NEAR(Pullback::Source::grid(origin, coupled, shifted),
                 math::fast_sinf(0.7f) * math::fast_cosf(-0.3f), 2e-3f);
  HS_EXPECT_EQ(Pullback::Source::primitive_lattice(origin, params), 1.0f);
  // Cell coordinate (0.2, 0.1) sits inside the softness band, so the edge ramp
  // resolves between its saturated ends.
  HS_EXPECT_NEAR(
      Pullback::Source::primitive_lattice(math::Complex(0.2f, 0.1f), params),
      0.38671f, 2e-3f);

  const Pullback::Source::SphericalRingsSourceParams ring_params;
  const Pullback::Source::PreparedSphericalRings ring_frame{math::Y_AXIS, 0.0f};
  HS_EXPECT_EQ(
      Pullback::Source::spherical_rings(math::X_AXIS, ring_params, ring_frame),
      1.0f);
  const math::Vector between_rings(0.9659258f, 0.2588190f, 0.0f);
  HS_EXPECT_EQ(
      Pullback::Source::spherical_rings(between_rings, ring_params, ring_frame),
      -1.0f);

  const Pullback::Source::FractalSourceParams fractal_params;
  HS_EXPECT_EQ(
      Pullback::Source::escape_fractal(origin, fractal_params, prepared), 1.0f);
  // c = 0.6 escapes on iteration 3 of 8, landing on contour cycle 0.4219.
  HS_EXPECT_NEAR(Pullback::Source::escape_fractal(math::Complex(1.2f, 0.0f),
                                                  fractal_params, prepared),
                 -0.38298f, 5e-3f);

  const Pullback::Source::TessellationSourceParams tessellation_params;
  HS_EXPECT_EQ(Pullback::Source::tessellation(
                   math::Complex(0.5f, 0.0f), tessellation_params,
                   Pullback::Source::TessellationKind::SQUARE, prepared),
               1.0f);
  HS_EXPECT_EQ(Pullback::Source::tessellation(
                   origin, tessellation_params,
                   Pullback::Source::TessellationKind::SQUARE, prepared),
               -1.0f);

  const math::Vector axis = Pullback::Lens::Glitch::apply(
      math::Vector(0.0f, 1.0f, 0.0f), TestFrame{});
  HS_EXPECT_EQ(axis.x, 0.0f);
  HS_EXPECT_EQ(axis.y, 1.0f);
  HS_EXPECT_EQ(axis.z, 0.0f);

  const projections::ProjectionKernelResult kernel{
      math::Complex(1.0f, 2.0f), 3, 4, 5, 3.0f, 6, 7, 8};
  const Pullback::ProjectionResult projected =
      Pullback::Projection::from_kernel(kernel, 2.0f);
  HS_EXPECT_EQ(projected.coords.re, 2.0f);
  HS_EXPECT_EQ(projected.coords.im, 4.0f);
  HS_EXPECT_EQ(projected.provenance.fade_edge_distance, 6.0f);
  HS_EXPECT_EQ(projected.provenance.region_id, 3);
  HS_EXPECT_EQ(projected.provenance.component_id, 4);
  HS_EXPECT_EQ(projected.provenance.boundary_flags, 5);
  HS_EXPECT_EQ(projected.provenance.flags, 6);
  HS_EXPECT_EQ(projected.provenance.traits, 7);
  HS_EXPECT_EQ(projected.provenance.edge_class, 8);
  HS_EXPECT_EQ(projected.provenance.value_weight, 1.0f);
  HS_EXPECT_EQ(projected.provenance.domain_coverage, 1.0f);
}

inline void test_pullback_warp_phase_loop() {
  Pullback::Warp::WaveShearParams wave;
  wave.field_angle = 0.9f;
  const auto wave_0 = Pullback::Warp::prepare(wave, 0.0f);
  const auto wave_1 = Pullback::Warp::prepare(wave, 1.0f);
  HS_EXPECT_NEAR(wave_0.rotation_cos, wave_1.rotation_cos, 1e-6f);
  HS_EXPECT_NEAR(wave_0.rotation_sin, wave_1.rotation_sin, 1e-6f);

  const Pullback::Warp::MirrorParams mirror{0.0f, 0.7f, 1.3f,
                                            0.9f, 0.4f, -0.2f};
  const auto mirror_0 = Pullback::Warp::prepare(mirror, 0.0f);
  const auto mirror_1 = Pullback::Warp::prepare(mirror, 1.0f);
  HS_EXPECT_NEAR(mirror_0.transform.mirror.offset_x,
                 mirror_1.transform.mirror.offset_x, 1e-6f);
  HS_EXPECT_NEAR(mirror_0.transform.mirror.offset_y,
                 mirror_1.transform.mirror.offset_y, 1e-6f);

  Pullback::Warp::VectorNoiseParams vector;
  vector.vector_angle = 0.4f;
  const auto vector_0 = Pullback::Warp::prepare(vector, 0.0f);
  const auto vector_1 = Pullback::Warp::prepare(vector, 1.0f);
  HS_EXPECT_NEAR(vector_0.transform.noise_loop.offset.x,
                 vector_1.transform.noise_loop.offset.x, 1e-6f);
  HS_EXPECT_NEAR(vector_0.transform.noise_loop.offset.y,
                 vector_1.transform.noise_loop.offset.y, 1e-6f);
  HS_EXPECT_NEAR(vector_0.transform.noise_loop.offset.z,
                 vector_1.transform.noise_loop.offset.z, 1e-6f);
  const math::Vector curl_0 = math::noise_projected_loop_offset(0.0f);
  const math::Vector curl_1 = math::noise_projected_loop_offset(1.0f);
  HS_EXPECT_NEAR(curl_0.x, curl_1.x, 1e-6f);
  HS_EXPECT_NEAR(curl_0.y, curl_1.y, 1e-6f);
  HS_EXPECT_NEAR(curl_0.z, curl_1.z, 1e-6f);

  Pullback::Warp::VortexParams vortex;
  vortex.center_x = 0.2f;
  vortex.center_y = -0.3f;
  vortex.center_orbit_radius = 0.8f;
  const auto vortex_0 = Pullback::Warp::prepare(vortex, 0.0f);
  const auto vortex_1 = Pullback::Warp::prepare(vortex, 1.0f);
  HS_EXPECT_NEAR(vortex_0.transform.vortex.center_x,
                 vortex_1.transform.vortex.center_x, 1e-6f);
  HS_EXPECT_NEAR(vortex_0.transform.vortex.center_y,
                 vortex_1.transform.vortex.center_y, 1e-6f);

  Pullback::Warp::AffineParams affine;
  affine.translation_x = 2.0f;
  affine.translation_y = -1.0f;
  affine.scale_x = 1.5f;
  affine.scale_y = 0.75f;
  affine.shear = 0.2f;
  const auto affine_0 = Pullback::Warp::prepare(affine, 0.0f, 0.6f, 2.0f);
  const auto affine_1 = Pullback::Warp::prepare(affine, 1.0f, 0.6f, 2.0f);
  HS_EXPECT_NEAR(affine_0.transform.affine.translation_x,
                 affine_1.transform.affine.translation_x, 1e-6f);
  HS_EXPECT_NEAR(affine_0.transform.affine.translation_y,
                 affine_1.transform.affine.translation_y, 1e-6f);
  HS_EXPECT_NEAR(affine_0.transform.affine.scale_x,
                 affine_1.transform.affine.scale_x, 1e-6f);
  HS_EXPECT_NEAR(affine_0.transform.affine.scale_y,
                 affine_1.transform.affine.scale_y, 1e-6f);
  HS_EXPECT_NEAR(affine_0.transform.affine.shear,
                 affine_1.transform.affine.shear, 1e-6f);
}

inline void test_pullback_periodic_ripple() {
  Pullback::Surface::PeriodicRippleParams params;
  params.strength = 0.15f;
  params.decay = 0.0f;
  params.thickness = 0.7f;
  params.center_polar = 0.0f;

  const math::Vector midpoint =
      math::Vector::from_spherical(0.0f, 0.5f * math::PI_F);
  const Pullback::SurfaceResult start = Pullback::Surface::periodic_ripple(
      midpoint, Pullback::Surface::prepare_ripple(params, 0.0f), true);
  HS_EXPECT_EQ(start.sphere.x, midpoint.x);
  HS_EXPECT_EQ(start.sphere.y, midpoint.y);
  HS_EXPECT_EQ(start.sphere.z, midpoint.z);
  HS_EXPECT_EQ(start.path_length, 0.0f);

  const Pullback::SurfaceResult crest = Pullback::Surface::periodic_ripple(
      midpoint, Pullback::Surface::prepare_ripple(params, 0.5f), true);
  HS_EXPECT_NEAR(crest.sphere.length(), 1.0f, 1e-4f);
  HS_EXPECT_NEAR(crest.path_length, 0.5f * params.strength, 2e-3f);

  const Pullback::SurfaceResult wrapped = Pullback::Surface::periodic_ripple(
      midpoint, Pullback::Surface::prepare_ripple(params, 1.0f), true);
  HS_EXPECT_EQ(wrapped.sphere.x, start.sphere.x);
  HS_EXPECT_EQ(wrapped.sphere.y, start.sphere.y);
  HS_EXPECT_EQ(wrapped.sphere.z, start.sphere.z);
  HS_EXPECT_EQ(wrapped.path_length, start.path_length);
}

struct AddLens : Pullback::ApproximationDefaults {
  static math::Vector apply(const math::Vector &input, const TestFrame &) {
    return math::Vector(input.x + 1.0f, input.y, input.z);
  }
};

struct ScaleLens : Pullback::ApproximationDefaults {
  static math::Vector apply(const math::Vector &input, const TestFrame &) {
    return math::Vector(input.x * 2.0f, input.y, input.z);
  }
};

inline void test_pullback_lens_stack() {
  // Consecutive Lens stages are the composition mechanism.
  using BoundAdd = Pullback::Stage::Lens<AddLens>::Bind<TestBinding>;
  using BoundScale = Pullback::Stage::Lens<ScaleLens>::Bind<TestBinding>;
  const TestFrame frame;
  const Pullback::SphereSample stacked = BoundScale::run(
      BoundAdd::run({math::Vector(3.0f, 2.0f, 1.0f), 0.0f}, frame, {}), frame,
      {});
  HS_EXPECT_EQ(stacked.dir.x, 8.0f);
  HS_EXPECT_EQ(stacked.dir.y, 2.0f);
  HS_EXPECT_EQ(stacked.dir.z, 1.0f);
}

struct SkyGradientStage
    : Pullback::Stage::Contract<SkyGradientStage, Pullback::SphereSample,
                                Color4> {
  using Policies = std::tuple<>;

  template <typename Binding>
  static Color4 run(const Pullback::SphereSample &input, const TestFrame &,
                    const Pullback::NoPrepared &) {
    return Color4(Pixel(4, 5, 6), input.dir.x);
  }
};

inline void test_pullback_rank_skip_crossing() {
  // A SPHERE -> COLOR crossing is admitted by the rules alone: a one-stage
  // chain over a consumer-authored descriptor validates and runs.
  using SkyPipeline = Pullback::Pipeline<TestBinding, SkyGradientStage>;
  static_assert(SkyPipeline::STAGE_COUNT == 1);
  static_assert(SkyPipeline::Validation::MONOTONE);
  static_assert(SkyPipeline::Validation::ENTRY);
  static_assert(SkyPipeline::Validation::EXIT);
  const SkyPipeline::Frame frame = SkyPipeline::prepare(TestFrame{});
  const Color4 sky = SkyPipeline::shade(math::Vector(0.5f, 0.0f, 0.0f), frame);
  HS_EXPECT_EQ(sky.color.g, 5);
  HS_EXPECT_EQ(sky.alpha, 0.5f);
}

struct FieldSnapParams {
  float value = 0.0f;
  uint8_t topology = 0;
  static constexpr auto FIELDS = std::array{Pullback::Field<FieldSnapParams>{
      "value", &FieldSnapParams::value, "Value", 0.0f, 1.0f}};
};

/**
 * @brief Pins the periodic field curves against the linear one.
 * @details A turns-valued field takes the short arc over a unit period, a
 *          radian-valued one over 2*pi; both wrap where a lerp would sweep
 *          the long way round.
 */
inline void test_pullback_field_curves() {
  using Pullback::FieldCurve;
  using Pullback::Fields::apply_curve;
  HS_EXPECT_NEAR(apply_curve(FieldCurve::SHORTEST_TURN, 0.9f, 0.1f, 0.25f),
                 0.95f, 1e-6f);
  HS_EXPECT_NEAR(apply_curve(FieldCurve::LERP, 0.9f, 0.1f, 0.25f), 0.7f, 1e-6f);
  HS_EXPECT_EQ(apply_curve(FieldCurve::SHORTEST_TURN, 0.9f, 0.1f, 0.0f), 0.9f);
  HS_EXPECT_EQ(apply_curve(FieldCurve::SHORTEST_TURN, 0.9f, 0.1f, 1.0f), 0.1f);
  HS_EXPECT_NEAR(apply_curve(FieldCurve::SHORTEST_PERIODIC,
                             math::TWO_PI_F - 0.1f, 0.1f, 0.25f),
                 math::TWO_PI_F - 0.05f, 1e-5f);

  Pullback::Surface::DirectSurfaceParams from;
  Pullback::Surface::DirectSurfaceParams to;
  from.direction = 0.9f;
  to.direction = 0.1f;
  const Pullback::Surface::DirectSurfaceParams mid =
      Pullback::Fields::interpolate(from, to, 0.25f);
  HS_EXPECT_NEAR(mid.direction, 0.95f, 1e-6f);
  HS_EXPECT_TRUE(Pullback::Fields::valid(mid));

  FieldSnapParams snap_from{0.0f, 3};
  FieldSnapParams snap_to{1.0f, 7};
  HS_EXPECT_EQ(Pullback::Fields::interpolate(snap_from, snap_to, 0.5f).topology,
               snap_from.topology);
  HS_EXPECT_EQ(Pullback::Fields::interpolate(snap_from, snap_to, 1.0f).topology,
               snap_to.topology);
}

inline void test_pullback_hard_edge_kernels() {
  using Pullback::Detail::smooth_ramp_or_step;
  HS_EXPECT_NEAR(smooth_ramp_or_step(0.0f, 1.0f, 0.5f), 0.5f, 1e-6f);
  HS_EXPECT_EQ(smooth_ramp_or_step(0.0f, 1.0f, -1.0f), 0.0f);
  HS_EXPECT_EQ(smooth_ramp_or_step(0.0f, 1.0f, 2.0f), 1.0f);
  HS_EXPECT_EQ(smooth_ramp_or_step(0.5f, 0.5f, 0.4f), 0.0f);
  HS_EXPECT_EQ(smooth_ramp_or_step(0.5f, 0.5f, 0.5f), 0.0f);
  HS_EXPECT_EQ(smooth_ramp_or_step(0.5f, 0.5f, 0.6f), 1.0f);

  using Pullback::Transfer::iso_contour;
  HS_EXPECT_EQ(iso_contour(0.5f, 0.5f, 0.1f), 1.0f);
  HS_EXPECT_EQ(iso_contour(0.8f, 0.5f, 0.1f), 0.0f);
  HS_EXPECT_EQ(iso_contour(0.5f, 0.5f, 0.0f), 1.0f);
  HS_EXPECT_EQ(iso_contour(0.5001f, 0.5f, 0.0f), 0.0f);

  Pullback::ProjectionProvenance provenance{};
  using Pullback::ProjectionCoverage::edge_fade;
  provenance.fade_edge_distance = 0.25f;
  HS_EXPECT_NEAR(edge_fade(provenance, 0.5f), 0.5f, 1e-6f);
  HS_EXPECT_EQ(edge_fade(provenance, 0.0f), 1.0f);
  provenance.fade_edge_distance = 0.0f;
  HS_EXPECT_EQ(edge_fade(provenance, 0.0f), 0.0f);

  using Pullback::ValueCoverage::value_cutout;
  HS_EXPECT_NEAR(value_cutout(0.5f, 0.5f, 0.25f), 0.5f, 1e-6f);
  HS_EXPECT_EQ(value_cutout(0.4f, 0.5f, 0.0f), 0.0f);
  HS_EXPECT_EQ(value_cutout(0.6f, 0.5f, 0.0f), 1.0f);
}

inline void test_pullback_hexagonal_edges() {
  using namespace Pullback::Source;
  const math::Complex vertices[] = {
      {1.0f, 0.0f},  {0.5f, 0.8660254038f},   {-0.5f, 0.8660254038f},
      {-1.0f, 0.0f}, {-0.5f, -0.8660254038f}, {0.5f, -0.8660254038f}};
  for (const math::Complex center :
       {math::Complex(0.0f, 0.0f), math::Complex(0.0f, 1.7320508076f),
        math::Complex(1.5f, -0.8660254038f)}) {
    for (float scale : {0.5f, 1.0f, 3.0f}) {
      for (float angle : {0.0f, 0.37f}) {
        TessellationSourceParams params;
        params.cell_scale = scale;
        PreparedSource prepared{};
        prepared.angle_cos = std::cos(angle);
        prepared.angle_sin = std::sin(angle);
        auto sample = [&](const math::Complex &point) {
          const float X = center.re + point.re;
          const float Y = center.im + point.im;
          const math::Complex input(
              (X * prepared.angle_cos - Y * prepared.angle_sin) / scale,
              (X * prepared.angle_sin + Y * prepared.angle_cos) / scale);
          return tessellation(input, params, TessellationKind::HEXAGONAL,
                              prepared);
        };
        HS_EXPECT_EQ(sample(math::Complex(0.0f, 0.0f)), -1.0f);
        for (int edge = 0; edge < 6; ++edge) {
          const math::Complex &a = vertices[edge];
          const math::Complex &b = vertices[(edge + 1) % 6];
          const math::Complex midpoint((a.re + b.re) * 0.5f,
                                       (a.im + b.im) * 0.5f);
          const float LENGTH = std::hypot(midpoint.re, midpoint.im);
          const math::Complex normal(midpoint.re / LENGTH,
                                     midpoint.im / LENGTH);
          HS_EXPECT_EQ(sample(a), 1.0f);
          HS_EXPECT_EQ(sample(midpoint), 1.0f);
          for (float direction : {-1.0f, 1.0f}) {
            HS_EXPECT_EQ(sample(math::Complex(
                             midpoint.re + direction * 0.08f * normal.re,
                             midpoint.im + direction * 0.08f * normal.im)),
                         -1.0f);
            HS_EXPECT_NEAR(sample(math::Complex(
                               midpoint.re + direction * 0.05f * normal.re,
                               midpoint.im + direction * 0.05f * normal.im)),
                           0.0f, 1e-4f);
          }
        }
      }
    }
  }
}

inline void test_pullback_displacement_overflow() {
  for (const math::Complex delta :
       {math::Complex(3.0f, 4.0f), math::Complex(0.0f, 0.0f),
        math::Complex(1e10f, -1e10f), math::Complex(5.22e19f, 5.22e19f),
        math::Complex(1e30f, -1e30f), math::Complex(1e30f, 0.0f)}) {
    HS_EXPECT_EQ(Pullback::Warp::displacement(delta, false), 0.0f);
    const float distance = Pullback::Warp::displacement(delta, true);
    const double EXPECTED = std::hypot(static_cast<double>(delta.re),
                                       static_cast<double>(delta.im));
    HS_EXPECT_TRUE(std::isfinite(distance));
    HS_EXPECT_NEAR(distance, EXPECTED, EXPECTED * 1e-6);
    const float SQUARED = delta.re * delta.re + delta.im * delta.im;
    if (std::isfinite(SQUARED))
      HS_EXPECT_EQ(distance, sqrtf(SQUARED));
  }
}

inline int run_pullback_tests() {
  ModuleFixture fixture("pullback");
  test_pullback_carrier_contract();
  test_pullback_validation_predicates();
  test_pullback_evaluation_order();
  test_pullback_placement_transparency();
  test_pullback_public_surface();
  test_pullback_no_instrumentation();
  test_pullback_counting_instrumentation();
  test_pullback_prepared_stage_policies();
  test_pullback_stage_combinators();
  test_pullback_provider_contracts();
  test_pullback_concrete_catalog();
  test_pullback_hexagonal_edges();
  test_pullback_displacement_overflow();
  test_pullback_warp_phase_loop();
  test_pullback_periodic_ripple();
  test_pullback_lens_stack();
  test_pullback_rank_skip_crossing();
  test_pullback_field_curves();
  test_pullback_hard_edge_kernels();
  return fixture.result();
}

} // namespace pullback_tests
} // namespace hs_test
