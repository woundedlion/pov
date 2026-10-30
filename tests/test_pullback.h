/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Pullback carrier, stage and policy contracts.
 */
#pragma once

#include <array>
#include <cstddef>
#include <cstring>
#include <type_traits>

#include "core/render/pullback.h"
#include "core/render/pullback/runtime_seeds.h"
#include "core/render/pullback/sphere_animation.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/pattern/fields.h"

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

/** @brief Records one pullback stage visit. */
inline void record(const TestFrame &frame, uint8_t id) {
  HS_CHECK(frame.call_count < frame.calls.size(),
           "pullback call recorder overflow");
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
    return {input.coords.re, input.provenance.domain_coverage, input.sphere,
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
                                Pullback::ColorSample> {
  using Policies = std::tuple<>;

  template <typename Binding>
  static Pullback::ColorSample run(const Pullback::FieldSample &input,
                                   const TestFrame &frame,
                                   const Pullback::NoPrepared &) {
    record(frame, 5);
    return {Color4(Pixel(static_cast<uint16_t>(input.path_length),
                         static_cast<uint16_t>(input.value * 100.0f), 3),
                   input.coverage),
            {}};
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
    : Pullback::Stage::Contract<ForeignCarrierStage, NotACarrier,
                                Pullback::ColorSample> {
  using Policies = std::tuple<>;

  template <typename Binding>
  static Pullback::ColorSample run(const NotACarrier &, const TestFrame &,
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
  using Output = Pullback::ColorSample;
  using Policies = std::tuple<>;

  template <typename Binding> struct Bind : Pullback::ApproximationDefaults {
    using Descriptor = WrongReturnStage;
    using Input = Pullback::FieldSample;
    using Output = Pullback::ColorSample;
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
  using Output = Pullback::ColorSample;
  using Policies = std::tuple<>;

  template <typename Binding> struct Bind : Pullback::ApproximationDefaults {
    using Descriptor = WrongPrepareStage;
    using Input = Pullback::FieldSample;
    using Output = Pullback::ColorSample;
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
                                Pullback::FieldSample, Pullback::ColorSample> {
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

/** @brief Carrier ranks and stage descriptors satisfy the public contract. */
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
  static_assert(Pullback::family_of<Pullback::ColorSample> == 3);
  static_assert(Pullback::family_of<NotACarrier> ==
                Pullback::FOREIGN_FAMILY_RANK);
  static_assert(Pullback::CanonicalCarrier<Pullback::SphereSample>);
  static_assert(Pullback::CanonicalCarrier<Pullback::ColorSample>);
  static_assert(!Pullback::CanonicalCarrier<NotACarrier>);
  static_assert(!Pullback::CanonicalCarrier<Pullback::ProjectionResult>);

  constexpr Pullback::ProjectionProvenance provenance{1, 2, 3, 0.25f, 0.5f, 4};
  HS_EXPECT_EQ(provenance.traits, 0);
  HS_EXPECT_EQ(provenance.edge_class, 0);
  HS_EXPECT_EQ(provenance.domain_coverage, 1.0f);
  HS_EXPECT_EQ(static_cast<uint8_t>(Pullback::ProjectionBoundary::CUT), 1);
  HS_EXPECT_EQ(static_cast<uint8_t>(Pullback::ProjectionBoundary::SINGULAR), 2);
}

/** @brief Projection traits preserve their declared branch and cut semantics. */
inline void test_pullback_projection_traits() {
  using projections::ProjectionTrait;
  const math::Vector direction(0.6f, 0.8f, 0.0f);
  const uint8_t singular =
      projections::projection_traits(ProjectionTrait::SINGULAR);
  const uint8_t folded =
      projections::projection_traits(ProjectionTrait::FOLDED);
  const uint8_t cut_singular = projections::projection_traits(
      ProjectionTrait::CUT, ProjectionTrait::SINGULAR);
  const uint8_t cut_singular_folded = projections::projection_traits(
      ProjectionTrait::CUT, ProjectionTrait::SINGULAR, ProjectionTrait::FOLDED);

  HS_EXPECT_EQ(
      Pullback::Projection::stereographic(direction, 1.0f).provenance.traits,
      singular);
  HS_EXPECT_EQ(Pullback::Projection::folded_sinusoidal(direction, 0.0f)
                   .provenance.traits,
               folded);
  HS_EXPECT_EQ(Pullback::Projection::equirectangular(direction, 0.0f, 1.0f)
                   .provenance.traits,
               cut_singular);
  HS_EXPECT_EQ(
      Pullback::Projection::gnomonic(
          direction, 1.0f, Pullback::Projection::GnomonicHemisphere::FRONT)
          .provenance.traits,
      cut_singular);
  HS_EXPECT_EQ(
      Pullback::Projection::gnomonic(
          direction, 1.0f, Pullback::Projection::GnomonicHemisphere::BACK)
          .provenance.traits,
      cut_singular);
  HS_EXPECT_EQ(
      Pullback::Projection::gnomonic(
          direction, 1.0f, Pullback::Projection::GnomonicHemisphere::FOLDED)
          .provenance.traits,
      cut_singular_folded);
}

/** @brief Pipeline validation accepts valid stages and rejects invalid contracts. */
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

/** @brief Stages evaluate in pipeline order. */
inline void test_pullback_evaluation_order() {
  const TestPipeline::Frame frame = TestPipeline::prepare(TestFrame{});
  const Color4 result =
      TestPipeline::shade(math::Vector(1.0f, 2.0f, 3.0f), frame);
  HS_EXPECT_EQ(frame.ctx.call_count, 6U);
  for (size_t index = 0; index < frame.ctx.call_count; ++index)
    HS_EXPECT_EQ(frame.ctx.calls[index], index);
  HS_EXPECT_EQ(result.color.r, 3);
  HS_EXPECT_EQ(result.color.g, 200);
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
    Pullback::Stage::Placed<Pullback::CodeEmission::OUT_OF_LINE_ITCM,
                            FieldStage>,
    ColorCrossingStage>;

/** @brief Placement wrappers preserve semantic stages and output. */
inline void test_pullback_placement_transparency() {
  // Placement wrappers are invisible to the semantic leaf view: identical
  // leaves, validation and output; only the execution-node view differs.
  static_assert(TestPipeline::STAGE_COUNT == 6);
  static_assert(GroupedPipeline::STAGE_COUNT == 6);
  static_assert(std::is_same_v<GroupedPipeline::StageAt<1>, CrossingStage>);
  static_assert(std::is_same_v<GroupedPipeline::StageAt<2>, PlaneStage>);
  static_assert(
      std::is_same_v<TestPipeline::StageAt<5>, GroupedPipeline::StageAt<5>>);
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
  HS_EXPECT_EQ(flat_result.color.r, 3);
  HS_EXPECT_EQ(grouped_result.color.r, 3);
  HS_EXPECT_EQ(flat_result.color.g, 200);
  HS_EXPECT_EQ(grouped_result.color.g, 200);
}

struct PreparedOffsetStage
    : Pullback::Stage::Contract<PreparedOffsetStage, Pullback::PlaneSample,
                                Pullback::PlaneSample> {
  using Policies = std::tuple<>;
  struct Prepared {
    float offset;
  };

  template <typename Binding> static Prepared prepare(const TestFrame &frame) {
    return {frame.edge_width * 4.0f};
  }

  template <typename Binding>
  static Pullback::PlaneSample run(const Pullback::PlaneSample &input,
                                   const TestFrame &frame,
                                   const Prepared &prepared) {
    record(frame, 2);
    Pullback::PlaneSample output = input;
    output.path_length += prepared.offset;
    return output;
  }
};

struct PreparedScaleStage
    : Pullback::Stage::Contract<PreparedScaleStage, Pullback::FieldSample,
                                Pullback::FieldSample> {
  using Policies = std::tuple<>;
  struct Prepared {
    float scale;
  };

  template <typename Binding> static Prepared prepare(const TestFrame &frame) {
    return {frame.cutout_threshold * 2.0f};
  }

  template <typename Binding>
  static Pullback::FieldSample run(const Pullback::FieldSample &input,
                                   const TestFrame &frame,
                                   const Prepared &prepared) {
    record(frame, 4);
    Pullback::FieldSample output = input;
    output.value *= prepared.scale;
    return output;
  }
};

using PreparedGroupedPipeline = Pullback::Pipeline<
    TestBinding, EntryStage,
    Pullback::Stage::Placed<Pullback::CodeEmission::OUT_OF_LINE_FLASH,
                            CrossingStage, PreparedOffsetStage>,
    FieldCrossingStage,
    Pullback::Stage::Placed<Pullback::CodeEmission::INLINE_ONLY,
                            PreparedScaleStage>,
    ColorCrossingStage>;

/** @brief prepare_into() + shade_prepared() match prepare() + shade() bit for bit. */
inline void test_pullback_type_erased_prepare() {
  using Pipeline = PreparedGroupedPipeline;
  TestFrame source;
  source.edge_width = 0.75f;
  const math::Vector view(1.0f, 2.0f, 3.0f);
  const Pipeline::Frame frame = Pipeline::prepare(source);
  const Color4 expected = Pipeline::shade(view, frame);

  const TestFrame ctx = source;
  alignas(Pipeline::PreparedTuple) unsigned char
      storage[sizeof(Pipeline::PreparedTuple)];
  Pipeline::prepare_into(ctx, storage);
  const Color4 erased = Pipeline::shade_prepared(view, ctx, storage);

  HS_EXPECT_EQ(expected.color.r, 4);
  HS_EXPECT_EQ(expected.color.g, 100);
  HS_EXPECT_EQ(erased.color.r, expected.color.r);
  HS_EXPECT_EQ(erased.color.g, expected.color.g);
  HS_EXPECT_EQ(erased.color.b, expected.color.b);
  HS_EXPECT_EQ(std::memcmp(&erased.alpha, &expected.alpha, sizeof(float)), 0);
  HS_EXPECT_EQ(ctx.call_count, 6U);
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

/** @brief Pipeline aliases expose the bound frame and descriptor types. */
inline void test_pullback_public_surface() {
  static_assert(std::is_same_v<TestPipeline::Binding, TestBinding>);
  static_assert(std::is_same_v<TestPipeline::FrameState, TestFrame>);
  static_assert(std::is_same_v<TestPipeline::StageAt<0>, EntryStage>);
  static_assert(std::is_same_v<TestPipeline::StageAt<5>, ColorCrossingStage>);
  static_assert(TestPipeline::any_stage<IsCrossingStage>);
  static_assert(!TestPipeline::any_stage<NoStageMatches>);
  static_assert(std::is_same_v<TestPipeline::StageMatching<IsColorStage>,
                               ColorCrossingStage>);
  static_assert(
      std::is_same_v<TestPipeline::StageMatching<NoStageMatches>, void>);
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
  static_assert(Pattern::SphereLutView::SIZE == 3456);
  HS_EXPECT_TRUE(std::is_empty_v<TestPipeline>);
}

/** @brief Disabled instrumentation has no state or recorded events. */
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
  static const math::Quaternion &frame_conjugate(const TestFrame &) {
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

/** @brief Instrumentation records stage visits in evaluation order. */
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

  CountingInstrumentation::count = 0;
  using SphereSource = Pullback::Stage::SampleSphere<
      CountingSphericalSourcePolicy>::Bind<CountingBinding>;
  const TestFrame frame;
  static_cast<void>(SphereSource::run({math::X_AXIS, 0.0f}, frame,
                                      SphereSource::prepare(frame)));
  HS_EXPECT_EQ(CountingInstrumentation::count, size_t{1});
  HS_EXPECT_EQ(CountingInstrumentation::events[0],
               Pullback::ProfileEvent::SOURCE);
}

/** @brief Prepared policies retain their frame-specific values. */
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
  HS_EXPECT_EQ(BoundColorize::run(field, frame, color_prepared).color.color.r,
               11);
}

/** @brief Stage combinators preserve carrier fields and coverage. */
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
  HS_EXPECT_EQ(transferred.value, math::unit_bell(1.0f));
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

/** @brief Providers bind to their declared descriptors and preserve carriers. */
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

/** @brief Concrete stages evaluate against their independent expected values. */
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
  HS_EXPECT_EQ(Pattern::twin_wave(Pattern::PlanePosition{origin}, prepared),
               0.0f);
  HS_EXPECT_EQ(Pattern::twin_wave(Pattern::PlanePosition{origin}, shifted),
               0.5f * (math::fast_sinf(0.7f) + math::fast_sinf(0.3f)));
  HS_EXPECT_NEAR(Pattern::twin_wave(Pattern::PlanePosition{off}, turned),
                 0.5f * (math::fast_sinf(0.6f) + math::fast_sinf(0.25f)),
                 2e-3f);
  HS_EXPECT_EQ(Pattern::rings(Pattern::PlanePosition{origin}, prepared), 0.0f);
  HS_EXPECT_NEAR(Pattern::rings(Pattern::PlanePosition{off}, shifted),
                 math::fast_sinf(-0.05f), 2e-3f);
  HS_EXPECT_NEAR(Pattern::spiral(Pattern::PlanePosition{origin}, prepared),
                 1.0f, 2e-3f);
  HS_EXPECT_NEAR(Pattern::spiral(Pattern::PlanePosition{off}, shifted),
                 math::fast_sinf(-0.05f - 3.0f * math::fast_atan2(0.25f, 0.6f)),
                 2e-3f);
  HS_EXPECT_EQ(Pattern::grid(Pattern::PlanePosition{origin}, params, prepared),
               0.0f);
  HS_EXPECT_NEAR(
      Pattern::grid(Pattern::PlanePosition{origin}, coupled, shifted),
      math::fast_sinf(0.7f) * math::fast_cosf(-0.3f), 2e-3f);
  HS_EXPECT_EQ(
      Pattern::primitive_lattice(Pattern::PlanePosition{origin}, params), 1.0f);
  // Cell coordinate (0.2, 0.1) sits inside the softness band, so the edge ramp
  // resolves between its saturated ends.
  HS_EXPECT_NEAR(Pattern::primitive_lattice(
                     Pattern::PlanePosition{math::Complex(0.2f, 0.1f)}, params),
                 0.38671f, 2e-3f);

  const Pullback::Source::SphericalRingsSourceParams ring_params;
  const Pattern::PreparedSphericalRings ring_frame{math::Y_AXIS, 0.0f};
  HS_EXPECT_EQ(
      Pattern::spherical_rings(Pattern::SphereDirection::from(math::X_AXIS),
                               ring_params, ring_frame),
      1.0f);
  const math::Vector between_rings(0.9659258f, 0.2588190f, 0.0f);
  HS_EXPECT_EQ(
      Pattern::spherical_rings(Pattern::SphereDirection::from(between_rings),
                               ring_params, ring_frame),
      -1.0f);

  const Pullback::Source::FractalSourceParams fractal_params;
  const auto fractal_frame = Pattern::prepare_fractal(fractal_params, prepared);
  HS_EXPECT_EQ(Pattern::escape_fractal(Pattern::PlanePosition{origin},
                                       fractal_params, fractal_frame),
               1.0f);
  // c = 0.6 escapes at loop index 3 of 8 (orbit 0.4219), i.e. 1.6875 cycles
  // at the default 4 contours.
  HS_EXPECT_NEAR(
      Pattern::escape_fractal(Pattern::PlanePosition{math::Complex(1.2f, 0.0f)},
                              fractal_params, fractal_frame),
      -0.38298f, 5e-3f);

  const Pullback::Source::TessellationSourceParams tessellation_params;
  HS_EXPECT_EQ(
      Pattern::tessellation(Pattern::PlanePosition{math::Complex(0.5f, 0.0f)},
                            tessellation_params,
                            Pattern::TessellationKind::SQUARE, prepared),
      1.0f);
  HS_EXPECT_EQ(
      Pattern::tessellation(Pattern::PlanePosition{origin}, tessellation_params,
                            Pattern::TessellationKind::SQUARE, prepared),
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
  const math::Vector north(0.5416752f, 0.6427876f, 0.5416752f);
  const math::Vector south(north.x, -north.y, north.z);
  for (const auto layout : {projections::PeirceLayout::HORIZONTAL,
                            projections::PeirceLayout::VERTICAL}) {
    const auto weight = [&](const math::Vector &direction) {
      return Pullback::Projection::peirce(direction, 0.0f,
                                          static_cast<uint8_t>(layout), 0.0f,
                                          false, 1.0f, 1.0f)
          .provenance.value_weight;
    };
    HS_EXPECT_GT(weight(south), 0.4f);
    HS_EXPECT_NEAR(weight(south), weight(north), 1e-6f);
  }
  HS_EXPECT_EQ(
      Pullback::Projection::peirce_weight(south, 1.0f, 0.0f, 1.0f, true), 0.0f);
  for (float longitude : {-2.5f, -0.25f, 0.25f, 2.5f}) {
    const math::Vector direction(cosf(longitude), 0.0f, sinf(longitude));
    for (float meridian : {-2.0f, 1.0f, 3.0f}) {
      const auto folded =
          Pullback::Projection::folded_sinusoidal(direction, meridian);
      const float relative = projections::wrap_longitude(longitude - meridian);
      HS_EXPECT_EQ(folded.provenance.region_id,
                   static_cast<uint8_t>(relative < 0.0f));
      const auto expected = projections::folded_sinusoidal(direction, meridian);
      HS_EXPECT_EQ(folded.coords.re, expected.re);
      HS_EXPECT_EQ(folded.coords.im, expected.im);
    }
  }
}

/** @brief Warp parameters close their configured phase loop. */
inline void test_pullback_warp_phase_loop() {
  const Pullback::Warp::MirrorParams mirror{0.0f, 0.7f, 1.3f,
                                            0.9f, 0.4f, -0.2f};
  const auto mirror_0 = Pullback::Warp::prepare(mirror, 0.0f);
  const auto mirror_1 = Pullback::Warp::prepare(mirror, 1.0f);
  const auto mirror_mid = Pullback::Warp::prepare(mirror, 0.25f);
  HS_EXPECT_GT(std::abs(mirror_mid.transform.mirror.offset_x -
                        mirror_0.transform.mirror.offset_x),
               0.1f);
  HS_EXPECT_NEAR(mirror_0.transform.mirror.offset_x,
                 mirror_1.transform.mirror.offset_x, 1e-6f);
  HS_EXPECT_NEAR(mirror_0.transform.mirror.offset_y,
                 mirror_1.transform.mirror.offset_y, 1e-6f);

  Pullback::Warp::VectorNoiseParams vector;
  vector.vector_angle = 0.4f;
  const auto vector_0 = Pullback::Warp::prepare(vector, 0.0f);
  const auto vector_1 = Pullback::Warp::prepare(vector, 1.0f);
  const auto vector_mid = Pullback::Warp::prepare(vector, 0.25f);
  HS_EXPECT_GT((vector_mid.transform.noise_loop.offset -
                vector_0.transform.noise_loop.offset)
                   .length(),
               0.1f);
  HS_EXPECT_NEAR(vector_0.transform.noise_loop.offset.x,
                 vector_1.transform.noise_loop.offset.x, 1e-6f);
  HS_EXPECT_NEAR(vector_0.transform.noise_loop.offset.y,
                 vector_1.transform.noise_loop.offset.y, 1e-6f);
  HS_EXPECT_NEAR(vector_0.transform.noise_loop.offset.z,
                 vector_1.transform.noise_loop.offset.z, 1e-6f);
  const math::Vector curl_0 = math::noise_projected_loop_offset(0.0f);
  const math::Vector curl_1 = math::noise_projected_loop_offset(1.0f);
  HS_EXPECT_GT((math::noise_projected_loop_offset(0.25f) - curl_0).length(),
               0.1f);
  HS_EXPECT_NEAR(curl_0.x, curl_1.x, 1e-6f);
  HS_EXPECT_NEAR(curl_0.y, curl_1.y, 1e-6f);
  HS_EXPECT_NEAR(curl_0.z, curl_1.z, 1e-6f);

  Pullback::Warp::VortexParams vortex;
  vortex.center_x = 0.2f;
  vortex.center_y = -0.3f;
  vortex.center_orbit_radius = 0.8f;
  const auto vortex_0 = Pullback::Warp::prepare(vortex, 0.0f);
  const auto vortex_1 = Pullback::Warp::prepare(vortex, 1.0f);
  const auto vortex_mid = Pullback::Warp::prepare(vortex, 0.25f);
  HS_EXPECT_GT(std::abs(vortex_mid.transform.vortex.center_x -
                        vortex_0.transform.vortex.center_x),
               0.1f);
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
  const auto affine_mid = Pullback::Warp::prepare(affine, 0.25f, 0.6f, 2.0f);
  HS_EXPECT_GT(std::abs(affine_mid.transform.affine.translation_x -
                        affine_0.transform.affine.translation_x),
               0.1f);
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

/** @brief Ripple phase repeats after its configured period. */
inline void test_pullback_periodic_ripple() {
  Pullback::Surface::PeriodicRippleParams clock;
  clock.period = 12.0f;
  HS_EXPECT_EQ(Pullback::Surface::ripple_cycle(0.0f, clock), 0.0f);
  HS_EXPECT_EQ(Pullback::Surface::ripple_cycle(6.0f, clock), 0.5f);
  HS_EXPECT_EQ(Pullback::Surface::ripple_cycle(12.0f, clock), 1.0f);
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

/** @brief Builds one explicitly phased ripple wave.
 * @param center Unit ripple center.
 * @param phase Wave crest distance in radians.
 * @param amplitude Peak displacement in radians.
 * @return Prepared ripple with refreshed rejection bounds. */
inline Pullback::Surface::PreparedRipple
ripple_test_wave(const math::Vector &center, float phase, float amplitude) {
  Animation::RippleWave wave;
  wave.center = center;
  wave.phase = phase;
  wave.amplitude = amplitude;
  wave.decay = 0.0f;
  wave.thickness = 0.7f;
  wave.sync();
  return {wave};
}

/** @brief Applies the individual ripple kernels in their supplied order.
 * @param view Unit sample direction.
 * @param prepared Active waves.
 * @param support_margin Optional closest distance to a wave support boundary.
 * @return Displaced direction. */
inline math::Vector ripple_test_sequential(
    math::Vector view,
    const Pullback::SphereAnimation::PreparedRipples &prepared,
    float *support_margin = nullptr) {
  for (uint8_t i = 0; i < prepared.count; ++i) {
    const auto &wave = prepared.waves[i].ripple;
    if (support_margin && wave.amplitude > 0.001f) {
      const float DISTANCE =
          acosf(hs::clamp(math::dot(view, wave.center), -1.0f, 1.0f));
      for (float edge : {wave.phase - 2.0f * wave.half_width(),
                         wave.phase + 2.0f * wave.half_width()})
        if (edge > 0.0f && edge < math::PI_F)
          *support_margin = fminf(*support_margin, fabsf(DISTANCE - edge));
    }
    view = Pullback::Surface::periodic_ripple(view, prepared.waves[i], false)
               .sphere;
  }
  return view;
}

/** @brief Empty ripple batches preserve the direction and zero travel. */
inline void test_pullback_grouped_ripple_identity() {
  const Pullback::SphereAnimation::PreparedRipples EMPTY{};
  const math::Vector VIEW = math::Vector(0.3f, 0.4f, 0.5f).normalized();
  for (bool track_path : {false, true}) {
    const auto RESULT =
        Pullback::SphereAnimation::grouped_ripples(VIEW, EMPTY, track_path);
    HS_EXPECT_EQ(RESULT.sphere.x, VIEW.x);
    HS_EXPECT_EQ(RESULT.sphere.y, VIEW.y);
    HS_EXPECT_EQ(RESULT.sphere.z, VIEW.z);
    HS_EXPECT_EQ(RESULT.path_length, 0.0f);
  }
}

/** @brief A crest rotates outward by its amplitude and reports that arc. */
inline void test_pullback_grouped_ripple_analytic_crest() {
  Pullback::SphereAnimation::PreparedRipples prepared{};
  prepared.count = 1;
  prepared.waves[0] = ripple_test_wave(math::Y_AXIS, math::PI_F * 0.5f, 0.075f);
  const math::Vector EXPECTED(cosf(0.075f), -sinf(0.075f), 0.0f);
  const auto TRACKED =
      Pullback::SphereAnimation::grouped_ripples(math::X_AXIS, prepared, true);
  const auto UNTRACKED =
      Pullback::SphereAnimation::grouped_ripples(math::X_AXIS, prepared, false);
  HS_EXPECT_NEAR(TRACKED.sphere.x, EXPECTED.x, 1e-6f);
  HS_EXPECT_NEAR(TRACKED.sphere.y, EXPECTED.y, 1e-6f);
  HS_EXPECT_NEAR(TRACKED.sphere.z, EXPECTED.z, 1e-6f);
  HS_EXPECT_NEAR(TRACKED.path_length, 0.075f, 1e-6f);
  HS_EXPECT_GT(TRACKED.path_length, 0.07f);
  HS_EXPECT_EQ(UNTRACKED.sphere.x, TRACKED.sphere.x);
  HS_EXPECT_EQ(UNTRACKED.sphere.y, TRACKED.sphere.y);
  HS_EXPECT_EQ(UNTRACKED.sphere.z, TRACKED.sphere.z);
  HS_EXPECT_EQ(UNTRACKED.path_length, 0.0f);
}

/** @brief Mixed-center batches preserve the launch order of noncommuting waves. */
inline void test_pullback_grouped_ripple_order() {
  Pullback::SphereAnimation::PreparedRipples prepared{};
  prepared.count = 4;
  prepared.waves[0] = ripple_test_wave(math::Y_AXIS, math::PI_F * 0.5f, 0.13f);
  prepared.waves[1] = ripple_test_wave(math::Z_AXIS, math::PI_F * 0.5f, 0.11f);
  prepared.waves[2] = ripple_test_wave(math::Y_AXIS, math::PI_F * 0.5f, 0.09f);
  prepared.waves[3] = ripple_test_wave(math::Y_AXIS, math::PI_F * 0.5f, 0.07f);
  const math::Vector EXPECTED = ripple_test_sequential(math::X_AXIS, prepared);
  const auto RESULT =
      Pullback::SphereAnimation::grouped_ripples(math::X_AXIS, prepared, true);
  HS_EXPECT_NEAR(RESULT.sphere.x, EXPECTED.x, 1e-5f);
  HS_EXPECT_NEAR(RESULT.sphere.y, EXPECTED.y, 1e-5f);
  HS_EXPECT_NEAR(RESULT.sphere.z, EXPECTED.z, 1e-5f);
  std::swap(prepared.waves[0], prepared.waves[3]);
  const math::Vector REVERSED = ripple_test_sequential(math::X_AXIS, prepared);
  HS_EXPECT_GT((REVERSED - EXPECTED).length(), 1e-4f);
}

/** @brief Preset-range batches remain finite, unit length and geometrically faithful. */
inline void test_pullback_grouped_ripple_parameter_grid() {
  using namespace Pullback::SphereAnimation;
  float max_error = 0.0f, max_interior_error = 0.0f;
  unsigned cutoff_exceptions = 0;
  for (float thickness : {0.7f, 1.4f})
    for (float decay : {0.0f, 5.0f})
      for (float period : {30.0f, 143.0f})
        for (bool random_center : {false, true}) {
          BurstRippleParams params;
          params.strength = 0.15f;
          params.thickness = thickness;
          params.decay = decay;
          params.period = period;
          params.burst_rate = 1.0f / 64.0f;
          params.random_center = random_center;
          RippleState state;
          state.seed = 149502542u;
          for (int frame = 0; frame < 256; ++frame) {
            advance_ripple(state, params);
            const auto PREPARED = prepare_ripples(params, state);
            for (int i = 0; i < 24; ++i) {
              const float Y = -1.0f + (2.0f * i + 1.0f) / 24.0f;
              const float AZIMUTH = i * 2.39996323f;
              const float RADIUS = sqrtf(1.0f - Y * Y);
              const math::Vector VIEW(RADIUS * cosf(AZIMUTH), Y,
                                      RADIUS * sinf(AZIMUTH));
              const auto RESULT = grouped_ripples(VIEW, PREPARED, true);
              const auto UNTRACKED = grouped_ripples(VIEW, PREPARED, false);
              float support_margin = math::PI_F;
              const auto EXPECTED =
                  ripple_test_sequential(VIEW, PREPARED, &support_margin);
              const float ERROR =
                  std::max({fabsf(RESULT.sphere.x - EXPECTED.x),
                            fabsf(RESULT.sphere.y - EXPECTED.y),
                            fabsf(RESULT.sphere.z - EXPECTED.z)});
              max_error = fmaxf(max_error, ERROR);
              if (support_margin > 5e-4f) {
                max_interior_error = fmaxf(max_interior_error, ERROR);
                HS_EXPECT_LE(ERROR, 1e-4f);
              } else
                cutoff_exceptions += ERROR > 1e-4f;
              HS_EXPECT_TRUE(std::isfinite(RESULT.sphere.x));
              HS_EXPECT_TRUE(std::isfinite(RESULT.sphere.y));
              HS_EXPECT_TRUE(std::isfinite(RESULT.sphere.z));
              HS_EXPECT_TRUE(std::isfinite(RESULT.path_length));
              HS_EXPECT_GE(RESULT.path_length, 0.0f);
              HS_EXPECT_NEAR(math::dot(RESULT.sphere, RESULT.sphere), 1.0f,
                             1e-5f);
              HS_EXPECT_LE(ERROR, 0.002f);
              HS_EXPECT_EQ(UNTRACKED.sphere.x, RESULT.sphere.x);
              HS_EXPECT_EQ(UNTRACKED.sphere.y, RESULT.sphere.y);
              HS_EXPECT_EQ(UNTRACKED.sphere.z, RESULT.sphere.z);
              HS_EXPECT_EQ(UNTRACKED.path_length, 0.0f);
            }
          }
        }
  std::printf("  ripple grid max=%g interior=%g cutoff_exceptions=%u\n",
              max_error, max_interior_error, cutoff_exceptions);
}

/** @brief Pole, pole-crossing and large-angle guards preserve sequential geometry. */
inline void test_pullback_grouped_ripple_guards() {
  using namespace Pullback::SphereAnimation;
  PreparedRipples prepared{};
  prepared.count = 3;
  prepared.waves[0] = ripple_test_wave(math::Z_AXIS, 0.5f * math::PI_F, 0.07f);
  prepared.waves[1] = ripple_test_wave(math::Y_AXIS, 0.32310888f, 0.15f);
  prepared.waves[2] = ripple_test_wave(math::Y_AXIS, 0.5f * math::PI_F, 0.2f);
  for (const math::Vector &view :
       {math::Y_AXIS, -math::Y_AXIS, math::X_AXIS, 1.0001f * math::X_AXIS,
        math::Vector(sinf(0.02f), cosf(0.02f), 0.0f)}) {
    const auto EXPECTED = ripple_test_sequential(view, prepared);
    const auto RESULT = grouped_ripples(view, prepared, true);
    HS_EXPECT_NEAR(RESULT.sphere.x, EXPECTED.x, 1e-5f);
    HS_EXPECT_NEAR(RESULT.sphere.y, EXPECTED.y, 1e-5f);
    HS_EXPECT_NEAR(RESULT.sphere.z, EXPECTED.z, 1e-5f);
    HS_EXPECT_TRUE(std::isfinite(RESULT.path_length));
  }
}

/** @brief Tiny-wave travel stays continuous across the radial-basis fallback. */
inline void test_pullback_grouped_ripple_path_continuity() {
  using namespace Pullback::SphereAnimation;
  PreparedRipples prepared{};
  prepared.count = 1;
  const float DISTANCE = asinf(0.01f);
  prepared.waves[0] =
      ripple_test_wave(math::Y_AXIS, DISTANCE + 0.99f * 0.175f, 0.002f);
  const float BELOW = DISTANCE - 1e-6f;
  const float ABOVE = DISTANCE + 1e-6f;
  const auto A = grouped_ripples(math::Vector(sinf(BELOW), cosf(BELOW), 0.0f),
                                 prepared, true);
  const auto B = grouped_ripples(math::Vector(sinf(ABOVE), cosf(ABOVE), 0.0f),
                                 prepared, true);
  HS_EXPECT_GT(A.path_length, 1e-5f);
  HS_EXPECT_GT(B.path_length, 1e-5f);
  HS_EXPECT_NEAR(A.path_length, B.path_length, 2e-6f);
  HS_EXPECT_LT((A.sphere - B.sphere).length(), 1e-5f);
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

struct DirectionColorStage
    : Pullback::Stage::Contract<DirectionColorStage, Pullback::SphereSample,
                                Pullback::ColorSample> {
  using Policies = std::tuple<>;
  template <typename Binding>
  static Pullback::ColorSample run(const Pullback::SphereSample &input,
                                   const TestFrame &,
                                   const Pullback::NoPrepared &) {
    return {Color4(Pixel(static_cast<uint16_t>(input.dir.x), 0, 0), 1.0f), {}};
  }
};

/** @brief Stacked lens stages apply their transforms in order. */
inline void test_pullback_lens_stack() {
  using Stacked =
      Pullback::Pipeline<TestBinding, Pullback::Stage::Lens<AddLens>,
                         Pullback::Stage::Lens<ScaleLens>, DirectionColorStage>;
  static_assert(Stacked::Validation::MONOTONE);
  static_assert(Stacked::Validation::CARRIERS);
  const auto frame = Stacked::prepare(TestFrame{});
  const Color4 stacked = Stacked::shade(math::Vector(3.0f, 2.0f, 1.0f), frame);
  HS_EXPECT_EQ(stacked.color.r, 8);
}
struct SkyGradientStage
    : Pullback::Stage::Contract<SkyGradientStage, Pullback::SphereSample,
                                Pullback::ColorSample> {
  using Policies = std::tuple<>;

  template <typename Binding>
  static Pullback::ColorSample run(const Pullback::SphereSample &input,
                                   const TestFrame &,
                                   const Pullback::NoPrepared &) {
    return {Color4(Pixel(4, 5, 6), input.dir.x), {}};
  }
};

/** @brief A direct sphere-to-color stage can cross intermediate ranks. */
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

/** @brief Zero-width edge kernels select the defined hard-step endpoints. */
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

/** @brief Grid providers accept custom phase carriers. */
inline void test_pullback_custom_grid_preparation() {
  /** @brief Grid state with a minimal phase carrier. */
  struct State {
    using Binding = TestBinding;  ///< Fixture binding.
    using FrameState = TestFrame; ///< Fixture frame.
    /** @brief Phases and rotation used by the grid. */
    struct Phases final {
      float primary = math::PI_F * 0.5f; ///< Primary wave phase.
      float secondary = 0.0f;            ///< Secondary wave phase.
      float angle_cos = 1.0f;            ///< Rotation cosine.
      float angle_sin = 0.0f;            ///< Rotation sine.
    };
    /** @brief Returns uncoupled grid controls. @return Grid controls. */
    static const Pullback::Source::GridSourceParams &
    params(const FrameState &) {
      static constexpr Pullback::Source::GridSourceParams VALUE{};
      return VALUE;
    }
    /** @brief Returns the fixed wave phases. @return Wave phases. */
    static Phases prepare(const FrameState &) { return {}; }
  };
  using Grid = Pullback::Source::Grid<State>;
  static_assert(Grid::PROVIDER_VALID<TestBinding>);
  const TestFrame FRAME{};
  const auto PREPARED = Grid::prepare(FRAME);
  const Pullback::PlaneSample INPUT{{0.0f, 0.0f}, {}, math::Z_AXIS, 0.0f};
  HS_EXPECT_NEAR(Grid::sample(INPUT, FRAME, PREPARED), 1.0f, 1e-4f);
}

/** @brief Hexagonal edge weights agree at corners and shared boundaries. */
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
        Pattern::PreparedSource prepared{};
        prepared.angle_cos = std::cos(angle);
        prepared.angle_sin = std::sin(angle);
        auto sample = [&](const math::Complex &point) {
          const float X = center.re + point.re;
          const float Y = center.im + point.im;
          const math::Complex input(
              (X * prepared.angle_cos - Y * prepared.angle_sin) / scale,
              (X * prepared.angle_sin + Y * prepared.angle_cos) / scale);
          return Pattern::tessellation(Pattern::PlanePosition{input}, params,
                                       Pattern::TessellationKind::HEXAGONAL,
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

/** @brief Large displacements remain finite with the bounded overflow policy. */
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

/** @brief Runtime seed constants match the public deterministic contract. */
inline void test_pullback_runtime_seed_contract() {
  HS_EXPECT_EQ(Pullback::EFFECT_NOISE_SEED, 1337);
  HS_EXPECT_EQ(Pullback::CAMERA_WALK_SEED, 1337);
  HS_EXPECT_EQ(Pullback::PROJECTION_WALK_SEED, 7331);
  HS_EXPECT_EQ(Pullback::HUE_NOISE_SEED, 6047);
}

/** @brief Palette and modulation endpoints for color-stage tests. */
struct PaletteFrame {
  const BakedPalette *palette = nullptr;
  const Pixel *hue_rotation = nullptr;
  const Pullback::Modulation::Prepared *modulation = nullptr;
  const Pullback::Modulation::Prepared *previous_modulation = nullptr;
  ColorModulation::Params mechanic{ColorModulation::Mechanic::HUE_ROTATION,
                                   0.25f};
  Pullback::Color::TemporalPaletteResources temporal{};
};
struct PaletteBinding {
  using FrameState = PaletteFrame;
  using Instrumentation = Pullback::NoInstrumentation;
};
struct PaletteProvider {
  using Binding = PaletteBinding;
  using FrameState = PaletteFrame;
  static auto mapping_weights(const PaletteFrame &) {
    return Pullback::Color::PaletteMappingWeights::single(
        Pullback::Color::PaletteMapping::LINEAR);
  }
  static float mapping_frequency(const PaletteFrame &) { return 1.0f; }
  static float mapping_phase(const PaletteFrame &) { return 0.0f; }
  static float oscillation_depth(const PaletteFrame &) { return 0.0f; }
  static float oscillation_phase(const PaletteFrame &) { return 0.0f; }
  static const BakedPalette &palette(const PaletteFrame &frame) {
    return *frame.palette;
  }
  static HueRotationLutView hue_rotation(const PaletteFrame &frame) {
    return {frame.hue_rotation, frame.hue_rotation != nullptr};
  }
  static auto brightness_envelope(const PaletteFrame &) {
    return Pullback::Color::BrightnessEnvelope::NONE;
  }
  static float brightness_bottom(const PaletteFrame &) { return 1.0f; }
  static float brightness_top(const PaletteFrame &) { return 1.0f; }
  static float opacity_low(const PaletteFrame &) { return 0.2f; }
  static float opacity_high(const PaletteFrame &) { return 0.6f; }
};
struct TemporalPaletteProvider : PaletteProvider {
  static auto temporal_resources(const PaletteFrame &frame) {
    return frame.temporal;
  }
};
struct PaletteLayerProvider {
  using Binding = PaletteBinding;
  using FrameState = PaletteFrame;
  static Pullback::Color::ModulationState prepare(const PaletteFrame &frame) {
    return {frame.modulation, frame.previous_modulation,
            Pattern::SourceKind::DIRECTIONAL_SINE, frame.mechanic, true};
  }
};
struct PaletteHueLayerProvider : PaletteLayerProvider {
  static Pullback::Color::ModulationState prepare(const PaletteFrame &frame) {
    auto result = PaletteLayerProvider::prepare(frame);
    result.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 0.25f};
    return result;
  }
};
struct PaletteSampleStage
    : Pullback::Stage::Contract<PaletteSampleStage, Pullback::SphereSample,
                                Pullback::FieldSample> {
  using Policies = std::tuple<>;
  template <typename B>
  static Pullback::FieldSample run(const Pullback::SphereSample &input,
                                   const PaletteFrame &,
                                   const Pullback::NoPrepared &) {
    return {0.5f, 0.5f, input.dir, input.path_length};
  }
};
struct ConstantPalette {
  Color4 color;
  Color4 get(float) const { return color; }
};
using TemporalPalettePipeline = Pullback::Pipeline<
    PaletteBinding, PaletteSampleStage,
    Pullback::Stage::Colorize<
        Pullback::Color::GeneratedPalette<TemporalPaletteProvider>>,
    Pullback::Stage::Modulate<Pullback::Color::Modulate<PaletteLayerProvider>>>;
using TemporalPalettePlacedPipeline = Pullback::Pipeline<
    PaletteBinding, PaletteSampleStage,
    Pullback::Stage::Placed<
        Pullback::CodeEmission::OUT_OF_LINE_FLASH,
        Pullback::Stage::Colorize<
            Pullback::Color::GeneratedPalette<TemporalPaletteProvider>>,
        Pullback::Stage::Modulate<
            Pullback::Color::Modulate<PaletteLayerProvider>>>>;

struct CountingPairedPalette
    : Pullback::Color::GeneratedPalette<TemporalPaletteProvider> {
  using Base = Pullback::Color::GeneratedPalette<TemporalPaletteProvider>;
  inline static int pairs = 0;
  inline static int endpoints = 0;
  static auto apply_pair(const Pullback::FieldSample &sample,
                         const PaletteFrame &frame, const Prepared &prepared) {
    ++pairs;
    return Base::apply_pair(sample, frame, prepared);
  }
  static auto apply_endpoint(const Pullback::FieldSample &sample,
                             const PaletteFrame &frame,
                             const Prepared &prepared, bool previous) {
    ++endpoints;
    return Base::apply_endpoint(sample, frame, prepared, previous);
  }
};
struct CountingEndpointPalette : CountingPairedPalette {
  static std::pair<Pullback::ColorSample, Pullback::ColorSample>
  apply_pair(const Pullback::FieldSample &, const PaletteFrame &,
             const Prepared &) = delete;
};

/** @brief Paired color dispatch skips individual endpoint calls during blends. */
inline void test_pullback_paired_color_dispatch() {
  using namespace Pullback;
  using Paired = Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                                    Stage::Colorize<CountingPairedPalette>>;
  using Fallback = Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                                      Stage::Colorize<CountingEndpointPalette>>;
  std::array<uint8_t, 2 * BakedPalette::required_arena_bytes()> buffer{};
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage a, b;
  a.bake(arena, ConstantPalette{Color4(Pixel(1000, 2000, 3000), 0.2f)});
  b.bake(arena, ConstantPalette{Color4(Pixel(9000, 6000, 1000), 0.6f)});
  PaletteFrame frame;
  frame.palette = &b.view();
  frame.temporal.palette = &a.view();
  for (float fraction : {0.0f, 0.25f, 1.0f}) {
    frame.temporal.fraction = fraction;
    CountingPairedPalette::pairs = 0;
    CountingPairedPalette::endpoints = 0;
    const auto COLOR = Paired::shade(math::Y_AXIS, Paired::prepare(frame));
    HS_EXPECT_NEAR(COLOR.color.r, 1000 + 8000 * fraction, 1.0f);
    HS_EXPECT_NEAR(COLOR.alpha, 0.04f + 0.08f * fraction, 1e-6f);
    HS_EXPECT_EQ(CountingPairedPalette::pairs, fraction == 0.25f ? 1 : 0);
    HS_EXPECT_EQ(CountingPairedPalette::endpoints, fraction == 0.25f ? 0 : 1);
    CountingPairedPalette::pairs = 0;
    CountingPairedPalette::endpoints = 0;
    const auto OTHER = Fallback::shade(math::Y_AXIS, Fallback::prepare(frame));
    HS_EXPECT_EQ(OTHER.color, COLOR.color);
    HS_EXPECT_EQ(OTHER.alpha, COLOR.alpha);
    HS_EXPECT_EQ(CountingPairedPalette::pairs, 0);
    HS_EXPECT_EQ(CountingPairedPalette::endpoints, fraction < 1.0f ? 2 : 1);
  }
}

struct CoordinatePalette {
  float alpha;
  Color4 get(float value) const {
    return Color4(Pixel(value * 60000, (1.0f - value) * 50000, 7000), alpha);
  }
};

/** @brief Paired palettes preserve shared controls and endpoint resources. */
inline void test_pullback_paired_palette_controls() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<TemporalPaletteProvider>;
  std::array<uint8_t, 2 * BakedPalette::required_arena_bytes()> buffer{};
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage a, b;
  a.bake(arena, CoordinatePalette{0.2f});
  b.bake(arena, CoordinatePalette{0.6f});
  std::array<Pixel, 2> hue{};
  PaletteFrame frame;
  frame.palette = &b.view();
  frame.temporal.palette = &a.view();
  auto prepared = Policy::prepare(frame);
  prepared.current.mapping_frequency = 2.0f;
  prepared.current.mapping_offset = 0.13f;
  prepared.current.brightness_envelope = Color::BrightnessEnvelope::ASCENDING;
  prepared.current.brightness_bottom = 0.2f;
  prepared.current.brightness_top = 0.8f;
  prepared.current.hue_rotation = {hue.data(), true};
  prepared.previous = prepared.current;
  prepared.previous.palette = &a.view();
  prepared.previous.hue_rotation = {hue.data() + 1, true};
  const FieldSample INPUT{0.37f, 0.45f, math::X_AXIS, 2.5f};
  const auto PAIR = Policy::apply_pair(INPUT, frame, prepared);
  for (const auto *sample : {&PAIR.first, &PAIR.second}) {
    HS_EXPECT_NEAR(sample->palette_coordinate, 0.87f, 1e-6f);
    HS_EXPECT_NEAR(sample->brightness_gain, 0.422f, 1e-6f);
    HS_EXPECT_EQ(sample->field.value, INPUT.value);
    HS_EXPECT_EQ(sample->field.coverage, INPUT.coverage);
    HS_EXPECT_EQ(sample->field.sphere, INPUT.sphere);
    HS_EXPECT_EQ(sample->field.path_length, INPUT.path_length);
    HS_EXPECT_EQ(sample->pending_hue, 0.0f);
    HS_EXPECT_EQ(sample->pending_lift, 0.0f);
    HS_EXPECT_EQ(sample->color.color, a.view().get(0.87f).color);
  }
  HS_EXPECT_NEAR(PAIR.first.color.alpha, 0.03132f, 1e-6f);
  HS_EXPECT_NEAR(PAIR.second.color.alpha, 0.09396f, 1e-6f);
  HS_EXPECT_TRUE(PAIR.first.hue_rotation == hue.data() + 1);
  HS_EXPECT_TRUE(PAIR.second.hue_rotation == hue.data());
}

/** @brief The complete color suffix interpolates immutable palette endpoints. */
inline void test_pullback_temporal_palette() {
  using namespace Pullback;
  std::array<uint8_t, 2 * BakedPalette::required_arena_bytes()> buffer{};
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage a, b;
  a.bake(arena, ConstantPalette{Color4(Pixel(1000, 2000, 3000), 0.2f)});
  b.bake(arena, ConstantPalette{Color4(Pixel(9000, 6000, 1000), 0.6f)});
  PaletteFrame frame;
  frame.palette = &b.view();
  frame.temporal.palette = &a.view();
  frame.mechanic.amount = 0.0f;
  for (float fraction : {0.0f, 0.25f, 0.75f, 1.0f}) {
    frame.temporal.fraction = fraction;
    const auto COLOR = TemporalPalettePipeline::shade(
        math::Y_AXIS, TemporalPalettePipeline::prepare(frame));
    HS_EXPECT_NEAR(COLOR.color.r, 1000 + 8000 * fraction, 1.0f);
    HS_EXPECT_NEAR(COLOR.color.g, 2000 + 4000 * fraction, 1.0f);
    HS_EXPECT_NEAR(COLOR.color.b, 3000 - 2000 * fraction, 1.0f);
    HS_EXPECT_NEAR(COLOR.alpha, 0.04f + 0.08f * fraction, 1e-6f);
  }
}

/** @brief Opacity modulates each endpoint before final color interpolation. */
inline void test_pullback_modulation_temporal_alpha() {
  using namespace Pullback;
  std::array<uint8_t, 2 * BakedPalette::required_arena_bytes()> buffer{};
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage a, b;
  const Pixel RGB(1234, 4567, 8910);
  a.bake(arena, ConstantPalette{Color4(RGB, 0.2f)});
  b.bake(arena, ConstantPalette{Color4(RGB, 0.6f)});
  Modulation::Config config;
  config.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  config.evaluation = Pattern::Evaluation::DIRECT;
  config.sine.frequency = 0.0f;
  auto clock = Modulation::init(config);
  Modulation::Resources<false> previous, current;
  clock.phase.phase_turns = 0.75f;
  previous.prepare(config, clock);
  clock.phase.phase_turns = 0.25f;
  current.prepare(config, clock);
  PaletteFrame frame;
  frame.palette = &b.view();
  frame.temporal.palette = &a.view();
  frame.modulation = &current.prepared;
  frame.previous_modulation = &previous.prepared;
  frame.mechanic = {ColorModulation::Mechanic::OPACITY, 1.0f};
  for (float fraction : {0.0f, 0.25f, 0.75f, 1.0f}) {
    frame.temporal.fraction = fraction;
    const auto COLOR = TemporalPalettePipeline::shade(
        math::Y_AXIS, TemporalPalettePipeline::prepare(frame));
    HS_EXPECT_EQ(COLOR.color, RGB);
    HS_EXPECT_NEAR(COLOR.alpha, 0.12f * fraction, 1e-6f);
    const auto PLACED = TemporalPalettePlacedPipeline::shade(
        math::Y_AXIS, TemporalPalettePlacedPipeline::prepare(frame));
    HS_EXPECT_EQ(PLACED.color, COLOR.color);
    HS_EXPECT_NEAR(PLACED.alpha, COLOR.alpha, 0.0f);
  }
  frame.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 1.0f};
  frame.temporal.fraction = 0.5f;
  const auto COLOR = TemporalPalettePipeline::shade(
      math::Y_AXIS, TemporalPalettePipeline::prepare(frame));
  const auto EXPECTED =
      ColorModulation::PositiveLightnessLift::apply(Color4(RGB), 0.5f, {1.0f})
          .color;
  HS_EXPECT_NEAR(COLOR.color.r, EXPECTED.r, 1.0f);
  HS_EXPECT_NEAR(COLOR.color.g, EXPECTED.g, 1.0f);
  HS_EXPECT_NEAR(COLOR.color.b, EXPECTED.b, 1.0f);
}

/** @brief Hue and lightness composition uses one final gamut mapping.
 * @param base Unmodified straight-alpha color.
 * @param turns Combined hue turns.
 * @param lift Combined lift toward white.
 * @return The transformed color. */
inline Color4 composed_modulation_reference(Color4 base, float turns,
                                            float lift) {
  const LinRGB rgb = pixel_to_linrgb(base.color);
  OKLab lab = linear_rgb_to_oklab_fast(rgb.r, rgb.g, rgb.b);
  lab.L += (1.0f - lab.L) * lift;
  float cosine, sine;
  hue_sincos(turns, cosine, sine);
  const OKLab rotated{lab.L, lab.a * cosine - lab.b * sine,
                      lab.a * sine + lab.b * cosine};
  return Color4(linrgb_to_pixel(oklab_to_linear_rgb_lut_gamut(rotated)),
                base.alpha);
}

/** @brief Hue layers compose across lightness lift without RGB round trips. */
inline void test_pullback_stacked_modulations() {
  using namespace Pullback;
  const Color4 BASE(Pixel(23000, 17000, 11000), 0.7f);
  const ColorSample INPUT{BASE, {0.4f, 0.6f, math::Y_AXIS, 2.5f}};
  Color::ModulationState first{nullptr,
                               nullptr,
                               Pattern::SourceKind::DISPLACEMENT,
                               {ColorModulation::Mechanic::HUE_ROTATION, 0.3f},
                               true};
  auto second = first;
  second.mechanic.amount = -0.3f;
  const auto CANCELLED = Color::finish(
      Color::apply_modulation(Color::apply_modulation(INPUT, first), second));
  HS_EXPECT_EQ(CANCELLED.color, BASE.color);
  HS_EXPECT_NEAR(CANCELLED.alpha, BASE.alpha, 0.0f);
  std::array<Pixel, HueRotationLutView::SIZE> table;
  for (size_t index = 0; index < table.size(); ++index)
    table[index] = Pixel((index % 16) * 1000, 2000, 3000);
  auto with_table = INPUT;
  with_table.hue_rotation = table.data();
  HS_EXPECT_EQ(
      Color::finish(Color::apply_modulation(
                        Color::apply_modulation(with_table, first), second))
          .color,
      BASE.color);
  Modulation::Config lift_config;
  lift_config.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  lift_config.evaluation = Pattern::Evaluation::DIRECT;
  lift_config.sine.frequency = 0.0f;
  auto clock = Modulation::init(lift_config);
  Modulation::Resources<false> lift_source;
  for (float phase : {0.0f, 0.75f}) {
    clock.phase.phase_turns = phase;
    const auto &source = lift_source.prepare(lift_config, clock);
    const Color::ModulationState IDENTITY_LIFT{
        &source,
        &source,
        Pattern::SourceKind::DIRECTIONAL_SINE,
        {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 1.0f},
        true};
    const auto LAYERED = Color::apply_modulation(
        Color::apply_modulation(Color::apply_modulation(with_table, first),
                                IDENTITY_LIFT),
        second);
    HS_EXPECT_EQ(Color::finish(LAYERED).color, BASE.color);
  }
  const auto ROTATED = Color::finish(Color::apply_modulation(INPUT, first));
  const auto EXPECTED = ColorModulation::HueRotation::apply(BASE, 2.5f, {0.3f});
  HS_EXPECT_EQ(ROTATED.color, EXPECTED.color);
  auto lift = first;
  lift.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 0.4f};
  second.mechanic.amount = 0.1f;
  auto actual = Color::apply_modulation(
      Color::apply_modulation(Color::apply_modulation(INPUT, first), lift),
      second);
  HS_EXPECT_EQ(actual.color.color, BASE.color);
  const auto COMPOSED = composed_modulation_reference(BASE, 0.0f, 0.4f);
  HS_EXPECT_EQ(Color::finish(actual).color, COMPOSED.color);
  const auto REORDERED = Color::finish(Color::apply_modulation(
      Color::apply_modulation(Color::apply_modulation(INPUT, first), second),
      lift));
  HS_EXPECT_EQ(REORDERED.color, COMPOSED.color);
  actual.brightness_gain = 0.25f;
  HS_EXPECT_EQ(Color::finish(actual).color, COMPOSED.color * 0.25f);
  auto opacity = first;
  opacity.mechanic = {ColorModulation::Mechanic::OPACITY, 1.0f};
  HS_EXPECT_NEAR(Color::finish(Color::apply_modulation(INPUT, opacity)).alpha,
                 BASE.alpha, 0.0f);
  HS_EXPECT_EQ(
      Color::finish(Color::apply_modulation(INPUT, lift)).color,
      ColorModulation::PositiveLightnessLift::apply(BASE, 1.0f, {0.4f}).color);
}

/** @brief Saturated lifted colors rotate before final clipping and opacity. */
inline void test_pullback_composed_lightness_hue() {
  using namespace Pullback;
  const Color::ModulationState HUE_A{
      nullptr,
      nullptr,
      Pattern::SourceKind::DISPLACEMENT,
      {ColorModulation::Mechanic::HUE_ROTATION, 0.25f},
      true};
  auto hue_b = HUE_A;
  hue_b.mechanic.amount = 0.5f;
  auto lift_a = HUE_A;
  lift_a.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 0.8f};
  auto lift_b = lift_a;
  lift_b.mechanic.amount = 0.5f;
  auto opacity = HUE_A;
  opacity.mechanic = {ColorModulation::Mechanic::OPACITY, 0.5f};
  for (const Pixel BASE : {Pixel(65535, 0, 0), Pixel(0, 65535, 0),
                           Pixel(0, 0, 65535), Pixel(62000, 61000, 60000)}) {
    const ColorSample INPUT{Color4(BASE, 0.7f),
                            {0.5f, 1.0f, math::Y_AXIS, 0.5f}};
    auto actual = Color::apply_modulation(INPUT, HUE_A);
    actual = Color::apply_modulation(actual, lift_a);
    actual = Color::apply_modulation(actual, hue_b);
    actual = Color::apply_modulation(actual, lift_b);
    HS_EXPECT_EQ(actual.color.color, BASE);
    const auto EXPECTED =
        composed_modulation_reference(INPUT.color, 0.375f, 0.55f);
    const auto FINISHED = Color::finish(actual);
    HS_EXPECT_NEAR(FINISHED.color.r, EXPECTED.color.r, 1.0f);
    HS_EXPECT_NEAR(FINISHED.color.g, EXPECTED.color.g, 1.0f);
    HS_EXPECT_NEAR(FINISHED.color.b, EXPECTED.color.b, 1.0f);
    HS_EXPECT_NEAR(FINISHED.alpha, 0.7f, 0.0f);
    actual = Color::apply_modulation(actual, opacity);
    actual.brightness_gain = 0.25f;
    const auto MASKED = Color::finish(actual);
    HS_EXPECT_EQ(MASKED.color, FINISHED.color * 0.25f);
    HS_EXPECT_NEAR(MASKED.alpha, 0.6125f, 1e-7f);
    actual.brightness_gain = 0.0f;
    const auto DARK = Color::finish(actual);
    HS_EXPECT_EQ(DARK.color, Pixel(0, 0, 0));
    HS_EXPECT_EQ(DARK.alpha, MASKED.alpha);
  }
}

/** @brief Temporal deferred color state blends before one perceptual finish. */
inline void test_pullback_temporal_composed_modulations() {
  using namespace Pullback;
  using Pipeline = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage,
      Stage::Colorize<Color::GeneratedPalette<TemporalPaletteProvider>>,
      Stage::Modulate<Color::Modulate<PaletteLayerProvider>>,
      Stage::Modulate<Color::Modulate<PaletteHueLayerProvider>>>;
  std::array<uint8_t, 2 * BakedPalette::required_arena_bytes()> buffer{};
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage a, b;
  const Color4 START(Pixel(65535, 0, 0), 0.2f);
  const Color4 END(Pixel(0, 0, 65535), 0.6f);
  a.bake(arena, ConstantPalette{START});
  b.bake(arena, ConstantPalette{END});
  Modulation::Config config;
  config.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  config.evaluation = Pattern::Evaluation::DIRECT;
  config.sine.frequency = 0.0f;
  auto clock = Modulation::init(config);
  Modulation::Resources<false> previous, current;
  clock.phase.phase_turns = 0.75f;
  previous.prepare(config, clock);
  clock.phase.phase_turns = 0.25f;
  current.prepare(config, clock);
  PaletteFrame frame;
  frame.palette = &b.view();
  frame.temporal.palette = &a.view();
  frame.modulation = &current.prepared;
  frame.previous_modulation = &previous.prepared;
  frame.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 0.4f};
  const auto PREVIOUS = composed_modulation_reference(START, -0.25f, 0.0f);
  const auto CURRENT = composed_modulation_reference(END, 0.25f, 0.4f);
  for (float fraction : {0.0f, 0.25f, 0.75f, 1.0f}) {
    frame.temporal.fraction = fraction;
    const auto COLOR = Pipeline::shade(math::Y_AXIS, Pipeline::prepare(frame));
    const auto BASE = START.color.lerp16(END.color, frac_to_q16(fraction));
    const auto EXPECTED = composed_modulation_reference(
        Color4(BASE), 0.75f - 0.5f * fraction, 0.4f * fraction);
    HS_EXPECT_NEAR(COLOR.color.r, EXPECTED.color.r, 1.0f);
    HS_EXPECT_NEAR(COLOR.color.g, EXPECTED.color.g, 1.0f);
    HS_EXPECT_NEAR(COLOR.color.b, EXPECTED.color.b, 1.0f);
    HS_EXPECT_NEAR(COLOR.alpha, 0.04f + 0.08f * fraction, 1e-6f);
    if (fraction > 0.0f && fraction < 1.0f) {
      const auto COMPLETED_BLEND =
          PREVIOUS.color.lerp16(CURRENT.color, frac_to_q16(fraction));
      HS_EXPECT_GT(std::abs(int(COLOR.color.r) - int(COMPLETED_BLEND.r)) +
                       std::abs(int(COLOR.color.g) - int(COMPLETED_BLEND.g)) +
                       std::abs(int(COLOR.color.b) - int(COMPLETED_BLEND.b)),
                   1000);
    }
  }
}

/** @brief Deferred hues take the short arc and preserve endpoint identities. */
inline void test_pullback_temporal_deferred_hue_wrap() {
  using namespace Pullback;
  ColorSample a{Color4(Pixel(65535, 0, 0), 0.2f), {}};
  ColorSample b{Color4(Pixel(0, 0, 65535), 0.8f), {}};
  a.pending_hue = 0.95f;
  b.pending_hue = 0.05f;
  a.pending_lift = 0.2f;
  b.pending_lift = 0.8f;
  a.brightness_gain = 0.25f;
  b.brightness_gain = 0.75f;
  const auto MIXED = blend_color_samples(a, b, 0.5f);
  HS_EXPECT_NEAR(MIXED.pending_hue, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(MIXED.pending_lift, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(MIXED.brightness_gain, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(MIXED.color.alpha, 0.5f, 1e-6f);
  HS_EXPECT_EQ(MIXED.color.color,
               a.color.color.lerp16(b.color.color, frac_to_q16(0.5f)));
  HS_EXPECT_EQ(blend_color_samples(a, b, 0.0f).color.color, a.color.color);
  HS_EXPECT_EQ(blend_color_samples(a, b, 1.0f).color.color, b.color.color);
  HS_EXPECT_EQ(blend_color_samples(a, b, 0.0f).pending_hue, a.pending_hue);
  HS_EXPECT_EQ(blend_color_samples(a, b, 1.0f).pending_hue, b.pending_hue);
}

/** @brief Runs pullback tests. */
inline int run_pullback_tests() {
  ModuleFixture fixture("pullback");
  HS_RUN(pattern_tests::test_pattern_direction_admission);
  HS_RUN(pattern_tests::test_pattern_sphere_formulas);
  HS_RUN(pattern_tests::test_pattern_source_clocks);
  HS_RUN(pattern_tests::test_pattern_noise_formulas);
  HS_RUN(pattern_tests::test_pattern_extracted_source_protocols);
  HS_RUN(test_pullback_runtime_seed_contract);
  HS_RUN(test_pullback_carrier_contract);
  HS_RUN(test_pullback_projection_traits);
  HS_RUN(test_pullback_validation_predicates);
  HS_RUN(test_pullback_evaluation_order);
  HS_RUN(test_pullback_placement_transparency);
  HS_RUN(test_pullback_type_erased_prepare);
  HS_RUN(test_pullback_public_surface);
  HS_RUN(test_pullback_no_instrumentation);
  HS_RUN(test_pullback_counting_instrumentation);
  HS_RUN(test_pullback_prepared_stage_policies);
  HS_RUN(test_pullback_stage_combinators);
  HS_RUN(test_pullback_provider_contracts);
  HS_RUN(test_pullback_concrete_catalog);
  HS_RUN(test_pullback_custom_grid_preparation);
  HS_RUN(test_pullback_hexagonal_edges);
  HS_RUN(test_pullback_displacement_overflow);
  HS_RUN(test_pullback_warp_phase_loop);
  HS_RUN(test_pullback_periodic_ripple);
  HS_RUN(test_pullback_grouped_ripple_identity);
  HS_RUN(test_pullback_grouped_ripple_analytic_crest);
  HS_RUN(test_pullback_grouped_ripple_order);
  HS_RUN(test_pullback_grouped_ripple_parameter_grid);
  HS_RUN(test_pullback_grouped_ripple_guards);
  HS_RUN(test_pullback_grouped_ripple_path_continuity);
  HS_RUN(test_pullback_lens_stack);
  HS_RUN(test_pullback_rank_skip_crossing);
  HS_RUN(test_pullback_field_curves);
  HS_RUN(test_pullback_hard_edge_kernels);
  HS_RUN(test_pullback_temporal_palette);
  HS_RUN(test_pullback_paired_color_dispatch);
  HS_RUN(test_pullback_paired_palette_controls);
  HS_RUN(test_pullback_modulation_temporal_alpha);
  HS_RUN(test_pullback_stacked_modulations);
  HS_RUN(test_pullback_composed_lightness_hue);
  HS_RUN(test_pullback_temporal_composed_modulations);
  HS_RUN(test_pullback_temporal_deferred_hue_wrap);
  return fixture.result();
}

} // namespace pullback_tests
} // namespace hs_test
