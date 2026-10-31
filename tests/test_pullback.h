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

#include "core/color/palette_cycler.h"
#include "core/render/pullback.h"
#include "core/render/pullback/composed_effect.h"
#include "core/render/pullback/runtime_seeds.h"
#include "core/render/pullback/sphere_animation.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/pattern/fields.h"
#include "tests/pixel_test_util.h"

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

/** @brief Circle-only lattice evaluation preserves the blended formula's bits. */
inline void test_pullback_lattice_circle_identity() {
  const auto reference = [](const Pattern::PlanePosition &coordinate,
                            const Pattern::LatticeSourceParams &params) {
    const float X =
        math::wrap_t(params.lattice_cell_scale * coordinate.value.re + 0.5f) -
        0.5f;
    const float Y =
        math::wrap_t(params.lattice_cell_scale * coordinate.value.im + 0.5f) -
        0.5f;
    const float CIRCLE = sqrtf(X * X + Y * Y) - params.lattice_radius;
    const float BX = fabsf(X) - params.lattice_radius;
    const float BY = fabsf(Y) - params.lattice_radius;
    const float SQUARE = sqrtf(fmaxf(BX, 0.0f) * fmaxf(BX, 0.0f) +
                               fmaxf(BY, 0.0f) * fmaxf(BY, 0.0f)) +
                         fminf(fmaxf(BX, BY), 0.0f);
    return 1.0f -
           2.0f * math::smooth_ramp(
                      -params.lattice_softness, params.lattice_softness,
                      hs::lerp(CIRCLE, SQUARE, params.lattice_shape_blend));
  };
  uint32_t random = 0x12345678u;
  for (float scale : {1.0f / 64.0f, 0.98971874f, 1.0f, 8.0f})
    for (float radius : {1.0f / 64.0f, 0.25f, 0.37757313f, 0.49f})
      for (float softness : {1.0f / 1024.0f, 0.1f, 0.34563965f, 1.0f})
        for (float blend : {0.0f, 0.00001f, 0.25f, 0.99999f, 1.0f}) {
          Pattern::LatticeSourceParams params;
          params.lattice_cell_scale = scale;
          params.lattice_radius = radius;
          params.lattice_softness = softness;
          params.lattice_shape_blend = blend;
          for (int sample = 0; sample < 256; ++sample) {
            random = random * 1664525u + 1013904223u;
            const float X = float(int32_t(random)) * (1024.0f / 2147483648.0f);
            random = random * 1664525u + 1013904223u;
            const float Y = float(int32_t(random)) * (1024.0f / 2147483648.0f);
            const Pattern::PlanePosition COORDINATE{
                math::Complex(sample < 8 ? float(sample - 4) * 0.5f : X,
                              sample < 8 ? 0.0f : Y)};
            const float EXPECTED = reference(COORDINATE, params);
            const float ACTUAL = Pattern::primitive_lattice(COORDINATE, params);
            HS_EXPECT_TRUE(std::memcmp(&EXPECTED, &ACTUAL, sizeof(float)) == 0);
          }
        }
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

/** @brief Supplied centers preserve the periodic wave's other prepared fields. */
inline void test_pullback_ripple_supplied_center() {
  Pullback::Surface::PeriodicRippleParams params;
  params.center_azimuth = 0.73f;
  params.center_polar = 1.17f;
  params.strength = 0.15f;
  params.thickness = 0.7f;
  params.decay = 0.2f;
  const auto CENTER =
      math::Vector::from_spherical(params.center_azimuth, params.center_polar);
  const auto OTHER = math::Vector(0.2f, 0.8f, -0.3f).normalized();
  for (float cycle : {-0.25f, 0.0f, 0.01f, 0.5f, 1.0f, 1.25f}) {
    const auto FIXED = Pullback::Surface::prepare_ripple(params, cycle);
    const auto SUPPLIED =
        Pullback::Surface::prepare_ripple(params, cycle, CENTER);
    HS_EXPECT_EQ(
        std::memcmp(&FIXED.ripple, &SUPPLIED.ripple, sizeof(FIXED.ripple)), 0);
    auto expected = FIXED;
    expected.ripple.center = OTHER;
    const auto ACTUAL = Pullback::Surface::prepare_ripple(params, cycle, OTHER);
    HS_EXPECT_EQ(
        std::memcmp(&expected.ripple, &ACTUAL.ripple, sizeof(expected.ripple)),
        0);
  }
}

/** @brief Sparse and stale burst slots retain canonical prepared waves. */
inline void test_pullback_ripple_burst_centers() {
  using namespace Pullback::SphereAnimation;
  BurstRippleParams params;
  params.period = 91.48529f;
  params.strength = 0.15f;
  params.thickness = 0.7f;
  params.decay = 0.2f;
  RippleState state;
  state.seed = 0x12345678;
  state.burst_count = 4;
  state.waves[4] = {5.0f, 1};
  state.waves[5] = {10.0f, 1};
  state.waves[8] = {params.period, 2};
  state.waves[9] = {params.period - 1.0f, 2};
  state.waves[10] = {params.period + 1.0f, 2};
  state.waves[11] = {4.0f, 2};
  state.waves[12] = {1.0f, 99};
  for (bool random : {false, true}) {
    params.random_center = random;
    const auto ACTUAL = prepare_ripples(params, state);
    PreparedRipples reference;
    for (uint32_t lag = 4; lag > 0; --lag) {
      const uint32_t BURST = state.burst_count - lag;
      for (size_t i = 0; i < 4; ++i) {
        const auto &wave = state.waves[(BURST % 4) * 4 + i];
        if (wave.age < 0.0f || wave.age >= params.period || wave.burst != BURST)
          continue;
        auto prepared =
            Pullback::Surface::prepare_ripple(params, wave.age / params.period);
        if (random)
          prepared.ripple.center = ripple_center(state.seed, BURST);
        reference.waves[reference.count++] = prepared;
      }
    }
    HS_EXPECT_EQ(ACTUAL.count, uint8_t{4});
    HS_EXPECT_EQ(ACTUAL.count, reference.count);
    for (uint8_t i = 0; i < ACTUAL.count; ++i)
      HS_EXPECT_EQ(std::memcmp(&ACTUAL.waves[i], &reference.waves[i],
                               sizeof(ACTUAL.waves[i])),
                   0);
    static RadialTables tables;
    prepare_radial_tables(params, state, tables);
    HS_EXPECT_EQ(tables.count, uint8_t{2});
    for (uint8_t i = 0; i < tables.count; ++i) {
      const auto EXPECTED =
          random ? ripple_center(state.seed, i + 1)
                 : Pullback::Surface::prepare_ripple_center(params);
      HS_EXPECT_EQ(tables.groups[i].center, EXPECTED);
    }
  }
}

/** @brief Validated dodecahedral matrices preserve canonical folds and wall fallback. */
inline void test_pullback_dodecahedral_matrix_fold() {
  static_assert(lenses::DodecahedralTable::GROUP.count == 120);
  size_t accepted = 0;
  uint32_t state = 0x735b18a3;
  const auto random = [&]() {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return float(state & 0xffffff) * (2.0f / 0xffffff) - 1.0f;
  };
  for (int sample = 0; sample < 8192; ++sample) {
    const auto INPUT = math::Vector(random(), random(), random()).normalized();
    const auto EXPECTED = lenses::dodecahedral_kaleidoscope_lens(INPUT);
    math::Vector candidate;
    const bool FAST = lenses::try_dodecahedral_table(INPUT, candidate);
    const auto ACTUAL = lenses::dodecahedral_table_lens(INPUT);
    HS_EXPECT_LE((ACTUAL - EXPECTED).length(), 2e-6f);
    HS_EXPECT_NEAR(math::dot(ACTUAL, ACTUAL), math::dot(INPUT, INPUT), 2e-6f);
    if (FAST) {
      ++accepted;
      HS_EXPECT_GT(candidate.x, lenses::DodecahedralTable::WALL_MARGIN);
      HS_EXPECT_GT(candidate.z, lenses::DodecahedralTable::WALL_MARGIN);
      HS_EXPECT_GT(math::dot(candidate, lenses::DODECAHEDRAL_MIRRORS[1]),
                   lenses::DodecahedralTable::WALL_MARGIN);
    }
  }
  HS_EXPECT_GT(accepted, size_t(6000));
  HS_EXPECT_LT(accepted, size_t(8192));
  const auto NORMAL = lenses::DODECAHEDRAL_MIRRORS[1];
  const math::Vector WALLS[] = {
      {0.0f, 0.95f, 0.3f},
      {0.2f, 0.98f, 0.0f},
      {0.1f, (-NORMAL.x * 0.1f - NORMAL.z * 0.1f) / NORMAL.y, 0.1f}};
  const math::Vector NORMALS[] = {math::X_AXIS, math::Z_AXIS, NORMAL};
  for (const auto &matrix : lenses::DodecahedralTable::GROUP.values)
    for (int wall = 0; wall < 3; ++wall)
      for (float offset : {-1e-6f, 0.0f, 1e-6f}) {
        const auto POINT = (WALLS[wall] + offset * NORMALS[wall]).normalized();
        const auto INPUT =
            math::Vector(matrix[0][0] * POINT.x + matrix[1][0] * POINT.y +
                             matrix[2][0] * POINT.z,
                         matrix[0][1] * POINT.x + matrix[1][1] * POINT.y +
                             matrix[2][1] * POINT.z,
                         matrix[0][2] * POINT.x + matrix[1][2] * POINT.y +
                             matrix[2][2] * POINT.z)
                .normalized();
        math::Vector candidate;
        HS_EXPECT_FALSE(lenses::try_dodecahedral_table(INPUT, candidate));
        const auto EXPECTED = lenses::dodecahedral_kaleidoscope_lens(INPUT);
        const auto ACTUAL = lenses::dodecahedral_table_lens(INPUT);
        HS_EXPECT_EQ(ACTUAL.x, EXPECTED.x);
        HS_EXPECT_EQ(ACTUAL.y, EXPECTED.y);
        HS_EXPECT_EQ(ACTUAL.z, EXPECTED.z);
      }
}

/** @brief Compile-time Euler curl matches the generic surface kernel. */
inline void test_pullback_euler_curl_specialization() {
  struct Frame {
    FastNoiseLite noise;
    float scale = 1.0f;
    float strength = 0.0f;
    bool track_path = false;
  };
  struct Provider {
    using FrameState = Frame;
    struct Binding {
      using FrameState = Frame;
    };
    static const FastNoiseLite &noise(const Frame &frame) {
      return frame.noise;
    }
    static float scale(const Frame &frame) { return frame.scale; }
    static float strength(const Frame &frame) { return frame.strength; }
    static Pullback::Surface::PreparedLoop prepare(const Frame &) {
      return Pullback::Surface::prepare(0.0f);
    }
    static bool path_length_required(const Frame &frame) {
      return frame.track_path;
    }
  };
  auto compare = []<math::NoiseBasis BASIS, math::TangentLimit LIMIT>() {
    using Policy =
        Pullback::Surface::CurlNoise<Provider, BASIS, Pullback::Surface::Euler,
                                     LIMIT>;
    static_assert(Policy::template PROVIDER_VALID<typename Provider::Binding>);
    Frame frame;
    for (int seed : {1, 1337, -7631}) {
      frame.noise.SetSeed(seed);
      for (float scale : {0.125f, 5.7742186f, 16.0f}) {
        frame.scale = scale;
        for (float strength : {0.0f, 0.00001f, 0.043f, 0.3f}) {
          frame.strength = strength;
          for (bool track_path : {false, true}) {
            frame.track_path = track_path;
            for (float phase : {0.0f, 0.37f, 0.9999f}) {
              const auto PREPARED = Pullback::Surface::prepare(phase);
              for (const auto &direction :
                   {math::X_AXIS, math::Y_AXIS, math::Z_AXIS, -math::X_AXIS,
                    math::Vector(0.3f, -0.4f, 0.5f).normalized()}) {
                const auto EXPECTED = Pullback::Surface::curl_noise(
                    direction, frame.noise, BASIS,
                    Pullback::Surface::Integrator::EULER, LIMIT, frame.scale,
                    PREPARED.loop_offset, frame.strength, frame.track_path);
                const auto ACTUAL = Policy::apply(direction, frame, PREPARED);
                HS_EXPECT_EQ(ACTUAL.sphere.x, EXPECTED.sphere.x);
                HS_EXPECT_EQ(ACTUAL.sphere.y, EXPECTED.sphere.y);
                HS_EXPECT_EQ(ACTUAL.sphere.z, EXPECTED.sphere.z);
                HS_EXPECT_EQ(ACTUAL.path_length, EXPECTED.path_length);
              }
            }
          }
        }
      }
    }
  };
  compare.template
  operator()<math::NoiseBasis::SIMPLEX, math::TangentLimit::SMOOTH>();
  compare.template
  operator()<math::NoiseBasis::SIMPLEX, math::TangentLimit::CLAMP>();
  compare.template
  operator()<math::NoiseBasis::FBM3, math::TangentLimit::SMOOTH>();
  compare
      .template operator()<math::NoiseBasis::FBM3, math::TangentLimit::CLAMP>();
  compare.template
  operator()<math::NoiseBasis::RIDGED3, math::TangentLimit::SMOOTH>();
  compare.template
  operator()<math::NoiseBasis::RIDGED3, math::TangentLimit::CLAMP>();
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

/** @brief Group rotations preserve the analytic rotation within their domain. */
inline void test_pullback_group_rotation_polynomial() {
  static_assert(RippleKernel::GROUP_ROTATION_POLYNOMIAL_MAX == 0.5f);
  const float OUTSIDE = std::nextafter(0.5f, 1.0f);
  const float ANGLES[] = {-0.9f,  -OUTSIDE, -0.5f, -0.49999f, -0.2f,
                          -0.15f, -0.05f,   0.0f,  0.05f,     0.15f,
                          0.2f,   0.49999f, 0.5f,  OUTSIDE,   0.9f};
  for (int row = 0; row < 13; ++row)
    for (int column = 0; column < 17; ++column) {
      const float Y = -0.99f + 1.98f * row / 12;
      const float ANGLE = math::TWO_PI_F * column / 17;
      const float RADIUS = sqrtf(1.0f - Y * Y);
      const math::Vector START(RADIUS * cosf(ANGLE), Y, RADIUS * sinf(ANGLE));
      const math::Vector CENTER = math::Vector(cosf(ANGLE * 0.37f + 0.2f), 0.3f,
                                               sinf(ANGLE * 0.71f + 0.4f))
                                      .normalized();
      const auto CROSS = math::cross(CENTER, START);
      const float CROSS_SQ = math::dot(CROSS, CROSS);
      if (CROSS_SQ < 1e-4f)
        continue;
      const auto AXIS = CROSS * (1.0f / sqrtf(CROSS_SQ));
      for (float angle : ANGLES) {
        math::Vector output;
        HS_EXPECT_TRUE(
            RippleKernel::rotate_group(START, CENTER, angle, output));
        const auto REFERENCE = RippleKernel::rotate(START, AXIS, angle);
        const auto DELTA = output - REFERENCE;
        HS_EXPECT_LE(math::dot(DELTA, DELTA), 1e-12f);
        HS_EXPECT_NEAR(math::dot(output, output), 1.0f, 1e-6f);
        const double COSINE = std::cos(double(angle));
        const double SINE = std::sin(double(angle));
        const double DOT = double(AXIS.x) * START.x + double(AXIS.y) * START.y +
                           double(AXIS.z) * START.z;
        const double expected[] = {
            START.x * COSINE +
                (double(AXIS.y) * START.z - double(AXIS.z) * START.y) * SINE +
                AXIS.x * DOT * (1.0 - COSINE),
            START.y * COSINE +
                (double(AXIS.z) * START.x - double(AXIS.x) * START.z) * SINE +
                AXIS.y * DOT * (1.0 - COSINE),
            START.z * COSINE +
                (double(AXIS.x) * START.y - double(AXIS.y) * START.x) * SINE +
                AXIS.z * DOT * (1.0 - COSINE)};
        const double ERROR_SQ =
            (output.x - expected[0]) * (output.x - expected[0]) +
            (output.y - expected[1]) * (output.y - expected[1]) +
            (output.z - expected[2]) * (output.z - expected[2]);
        HS_EXPECT_LE(ERROR_SQ, 1e-12);
        if (fabsf(angle) > 0.5f) {
          HS_EXPECT_EQ(output.x, REFERENCE.x);
          HS_EXPECT_EQ(output.y, REFERENCE.y);
          HS_EXPECT_EQ(output.z, REFERENCE.z);
        }
      }
    }
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

/** @brief Baked validity preserves every endpoint-cell admission and transfer bit. */
inline void test_pullback_radial_cell_validity() {
  using namespace Pullback::SphereAnimation;
  static RadialTables tables;
  tables = {};
  HS_EXPECT_FALSE(tables.validity_ready);
  tables.enabled = true;
  tables.count = 1;
  auto &group = tables.groups[0];
  group.previous.fill({.04f, .1f});
  group.current.fill({.08f, .2f});
  const auto CHECK = [&]() {
    prepare_radial_validity(tables);
    HS_EXPECT_TRUE(tables.validity_ready);
    for (size_t cell = 0; cell < RadialGroup::SAMPLES - 1; ++cell) {
      const bool VALID = !(group.previous[cell].path < 0.0f ||
                           group.previous[cell + 1].path < 0.0f ||
                           group.current[cell].path < 0.0f ||
                           group.current[cell + 1].path < 0.0f);
      HS_EXPECT_EQ(bool(tables.valid_cells[0][cell]), VALID);
      const float DISTANCE =
          (float(cell) + .375f) * (math::PI_F / (RadialGroup::SAMPLES - 1));
      for (bool track_path : {false, true}) {
        RadialTransfer expected{.123f, .456f}, actual = expected;
        const bool EXPECTED = RippleKernel::sample_transfer(
            group.previous, group.current, DISTANCE, .625f, track_path,
            expected);
        const bool ACTUAL = RippleKernel::sample_transfer(
            group.previous, group.current, DISTANCE, .625f, track_path, actual,
            tables.valid_cells[0].data());
        HS_EXPECT_EQ(ACTUAL, EXPECTED);
        HS_EXPECT_EQ(std::memcmp(&actual, &expected, sizeof(actual)), 0);
      }
    }
  };
  CHECK();
  for (int endpoint : {0, 1})
    for (size_t node : {size_t{0}, size_t{1}, size_t{128}, size_t{256}}) {
      auto &curve = endpoint == 0 ? group.previous : group.current;
      const float SAVED = curve[node].path;
      curve[node].path = -1.0f;
      CHECK();
      curve[node].path = SAVED;
    }
  group.previous[0].path = -0.0f;
  group.current[0].path = std::numeric_limits<float>::quiet_NaN();
  CHECK();
  RadialTransfer expected{.123f, .456f}, actual = expected;
  const bool EXPECTED = RippleKernel::sample_transfer(
      group.previous, group.current, math::PI_F, .25f, true, expected);
  const bool ACTUAL = RippleKernel::sample_transfer(
      group.previous, group.current, math::PI_F, .25f, true, actual,
      tables.valid_cells[0].data());
  HS_EXPECT_EQ(ACTUAL, EXPECTED);
  HS_EXPECT_EQ(std::memcmp(&actual, &expected, sizeof(actual)), 0);
}

/** @brief Owned masks match legacy radial output through supported and rejected updates. */
inline void test_pullback_radial_validity_lifetime() {
  using namespace Pullback::SphereAnimation;
  static RadialTables tables;
  tables = {};
  BurstRippleParams params;
  params.period = 80.0f;
  params.burst_rate = 1.0f / 64.0f;
  RippleState state;
  state.seed = 149502542u;
  uint64_t successes = 0;
  for (int frame = 0; frame < 600; ++frame) {
    if (frame == 173)
      params.strength = .12f;
    if (frame == 279)
      params.period = 143.0f;
    if (frame == 381)
      params.random_center = 0;
    if (frame == 419) {
      auto unsupported = params;
      unsupported.thickness = .1f;
      prepare_radial_tables(unsupported, state, tables, true);
      HS_EXPECT_FALSE(tables.enabled);
      HS_EXPECT_FALSE(tables.validity_ready);
    }
    advance_ripple(state, params);
    prepare_radial_tables(params, state, tables, true);
    HS_EXPECT_TRUE(tables.enabled);
    HS_EXPECT_TRUE(tables.validity_ready);
    for (int probe = 0; probe < 32; ++probe) {
      const float AZIMUTH = math::TWO_PI_F * probe / 32;
      const float Y = -.95f + 1.9f * ((probe * 17 + frame) % 32) / 31;
      const float RADIUS = sqrtf(1.0f - Y * Y);
      const math::Vector VIEW{RADIUS * cosf(AZIMUTH), Y,
                              RADIUS * sinf(AZIMUTH)};
      for (bool track_path : {false, true}) {
        Pullback::SurfaceResult actual, expected;
        const bool ACTUAL =
            try_radial_ripples(VIEW, track_path, &tables, .375f, actual);
        tables.validity_ready = false;
        const bool EXPECTED =
            try_radial_ripples(VIEW, track_path, &tables, .375f, expected);
        tables.validity_ready = true;
        HS_EXPECT_EQ(ACTUAL, EXPECTED);
        HS_EXPECT_EQ(std::memcmp(&actual, &expected, sizeof(actual)), 0);
        successes += ACTUAL;
      }
    }
  }
  HS_EXPECT_GT(successes, uint64_t{10000});
  prepare_radial_tables(params, state, tables);
  HS_EXPECT_FALSE(tables.validity_ready);
  RadialTables manual{};
  HS_EXPECT_FALSE(manual.validity_ready);
  manual.enabled = true;
  manual.count = 1;
  auto &group = manual.groups[0];
  group.center = math::Y_AXIS;
  group.cos_threshold_min = 1.0f;
  group.cos_threshold_max = -1.0f;
  group.previous.fill({0.0f, .1f});
  group.current.fill({0.0f, .1f});
  Pullback::SurfaceResult result;
  HS_EXPECT_TRUE(try_radial_ripples(math::X_AXIS, true, &manual, .5f, result));
  group.current[128].path = -1.0f;
  HS_EXPECT_FALSE(try_radial_ripples(math::X_AXIS, true, &manual, .5f, result));
}

/** @brief Reused radial endpoints match fresh bakes across clock and parameter changes. */
inline void test_pullback_radial_ripple_endpoint_reuse() {
  using namespace Pullback::SphereAnimation;
  static RadialTables cached, fresh;
  const auto CHECK = [&](const BurstRippleParams &params,
                         const RippleState &state) {
    prepare_radial_tables(params, state, cached, true);
    fresh.enabled = false;
    prepare_radial_tables(params, state, fresh, true);
    HS_EXPECT_EQ(cached.enabled, fresh.enabled);
    HS_EXPECT_EQ(cached.count, fresh.count);
    HS_EXPECT_EQ(cached.validity_ready, fresh.validity_ready);
    if (fresh.validity_ready)
      HS_EXPECT_EQ(std::memcmp(cached.valid_cells.data(),
                               fresh.valid_cells.data(),
                               sizeof(cached.valid_cells)),
                   0);
    for (uint8_t group = 0; group < fresh.count; ++group) {
      const auto &ACTUAL = cached.groups[group];
      const auto &EXPECTED = fresh.groups[group];
      HS_EXPECT_EQ(
          std::memcmp(&ACTUAL.center, &EXPECTED.center, sizeof(ACTUAL.center)),
          0);
      HS_EXPECT_EQ(std::memcmp(&ACTUAL.cos_threshold_min,
                               &EXPECTED.cos_threshold_min, sizeof(float)),
                   0);
      HS_EXPECT_EQ(std::memcmp(&ACTUAL.cos_threshold_max,
                               &EXPECTED.cos_threshold_max, sizeof(float)),
                   0);
      HS_EXPECT_EQ(std::memcmp(ACTUAL.previous.data(), EXPECTED.previous.data(),
                               sizeof(ACTUAL.previous)),
                   0);
      HS_EXPECT_EQ(std::memcmp(ACTUAL.current.data(), EXPECTED.current.data(),
                               sizeof(ACTUAL.current)),
                   0);
    }
  };
  for (bool randomize : {false, true}) {
    cached.enabled = false;
    BurstRippleParams params;
    params.period = 143.0f;
    params.burst_rate = 1.0f / 64.0f;
    RippleState state;
    state.seed = 149502542u;
    uint32_t random = 0x719463a5u;
    const auto NEXT = [&]() {
      random ^= random << 13;
      random ^= random >> 17;
      random ^= random << 5;
      return static_cast<float>(random & 0xffffffu) / 16777216.0f;
    };
    for (int frame = 0; frame < 600; ++frame) {
      if (randomize && frame % 29 == 0) {
        params.period = 60.0f + 83.0f * NEXT();
        params.strength = 0.01f + 0.14f * NEXT();
        params.decay = 5.0f * NEXT();
        params.thickness = 0.7f + 0.7f * NEXT();
        params.burst_rate = (1.0f + NEXT()) / 128.0f;
        params.burst_size = static_cast<uint8_t>(4.0f * NEXT());
        params.random_center = static_cast<uint8_t>(NEXT() >= 0.5f);
        params.center_azimuth = math::TWO_PI_F * NEXT();
        params.center_polar = math::PI_F * NEXT();
        state.seed = random;
      }
      if (randomize && frame % 83 == 0)
        for (auto &wave : state.waves)
          wave.age += 0.25f;
      if (randomize && frame % 101 == 0)
        for (int skipped = 0; skipped < 7; ++skipped)
          advance_ripple(state, params);
      advance_ripple(state, params);
      CHECK(params, state);
      if (randomize && frame % 47 == 0) {
        auto unsupported = params;
        unsupported.thickness = 0.1f;
        CHECK(unsupported, state);
        CHECK(params, state);
      }
    }
  }
  cached.enabled = false;
  BurstRippleParams params;
  params.period = 143.0f;
  RippleState reordered;
  reordered.burst_count = 4;
  for (uint32_t burst = 0; burst < 4; ++burst)
    reordered.waves[burst * 4] = {20.0f * (burst + 1), burst};
  CHECK(params, reordered);
  for (uint32_t burst = 0; burst < 4; ++burst)
    reordered.waves[burst * 4].age = 81.0f - 20.0f * burst;
  CHECK(params, reordered);

  for (int control = 0; control < 4; ++control) {
    cached.enabled = false;
    params = {};
    params.period = 143.0f;
    reordered.waves[0] = {40.0f, 0};
    CHECK(params, reordered);
    advance_ripple(reordered, params);
    switch (control) {
    case 0:
      params.period = 142.0f;
      break;
    case 1:
      params.strength = 0.125f;
      break;
    case 2:
      params.decay = 0.2f;
      break;
    case 3:
      params.thickness = 0.8f;
      break;
    }
    CHECK(params, reordered);
  }

  RadialCurveKey key{{80.0f, 0.15f, 0.1f, 0.7f}, {20.0f, 4.0f}, 2};
  HS_EXPECT_TRUE(key.matches(key));
  for (size_t scalar = 0; scalar < key.shape.size(); ++scalar) {
    auto different = key;
    different.shape[scalar] = std::nextafter(key.shape[scalar], 1000.0f);
    HS_EXPECT_FALSE(key.matches(different));
  }
  for (uint8_t wave = 0; wave < key.count; ++wave) {
    auto different = key;
    different.ages[wave] = std::nextafter(key.ages[wave], 1000.0f);
    HS_EXPECT_FALSE(key.matches(different));
  }
  auto different = key;
  --different.count;
  HS_EXPECT_FALSE(key.matches(different));
  different = key;
  different.ages[0] = -0.0f;
  key.ages[0] = 0.0f;
  HS_EXPECT_FALSE(key.matches(different));
}

/** @brief Radial endpoints reproduce an analytic single-wave crest. */
inline void test_pullback_radial_ripple_crest() {
  using namespace Pullback::SphereAnimation;
  BurstRippleParams params;
  params.period = 80.0f;
  params.decay = 0.0f;
  params.random_center = 0;
  params.center_polar = 0.0f;
  RippleState state;
  state.burst_count = 1;
  state.waves[0] = {40.0f, 0};
  static RadialTables tables;
  prepare_radial_tables(params, state, tables);
  HS_EXPECT_TRUE(tables.enabled);
  HS_EXPECT_EQ(tables.count, 1);
  auto prepared = prepare_ripples(params, state);
  const auto RESULT =
      radial_ripples(math::X_AXIS, prepared, true, &tables, 1.0f);
  HS_EXPECT_NEAR(RESULT.sphere.x, cosf(0.075f), 1e-6f);
  HS_EXPECT_NEAR(RESULT.sphere.y, -sinf(0.075f), 1e-6f);
  HS_EXPECT_NEAR(RESULT.sphere.z, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(RESULT.path_length, 0.075f, 1e-6f);
  const auto UNTRACKED =
      radial_ripples(math::X_AXIS, prepared, false, &tables, 1.0f);
  HS_EXPECT_EQ(UNTRACKED.sphere.x, RESULT.sphere.x);
  HS_EXPECT_EQ(UNTRACKED.sphere.y, RESULT.sphere.y);
  HS_EXPECT_EQ(UNTRACKED.path_length, 0.0f);
}

/** @brief Unsupported clocks and narrow wavelets retain grouped evaluation. */
inline void test_pullback_radial_ripple_fallback() {
  using namespace Pullback::SphereAnimation;
  BurstRippleParams params;
  RippleState state;
  state.burst_count = 1;
  state.waves[0] = {20.0f, 0};
  static RadialTables tables;
  for (bool narrow : {false, true}) {
    params.period = narrow ? 80.0f : 30.0f;
    params.thickness = narrow ? 0.1f : 0.7f;
    prepare_radial_tables(params, state, tables);
    HS_EXPECT_FALSE(tables.enabled);
    auto prepared = prepare_ripples(params, state);
    const math::Vector VIEW = math::Vector(0.4f, 0.5f, 0.6f).normalized();
    const auto EXPECTED = grouped_ripples(VIEW, prepared, true);
    const auto RESULT = radial_ripples(VIEW, prepared, true, &tables, 1.0f);
    HS_EXPECT_EQ(RESULT.sphere.x, EXPECTED.sphere.x);
    HS_EXPECT_EQ(RESULT.sphere.y, EXPECTED.sphere.y);
    HS_EXPECT_EQ(RESULT.sphere.z, EXPECTED.sphere.z);
    HS_EXPECT_EQ(RESULT.path_length, EXPECTED.path_length);
  }
  params = {};
  prepare_radial_tables(params, RippleState{}, tables);
  HS_EXPECT_TRUE(tables.enabled);
  HS_EXPECT_EQ(tables.count, 0);
  PreparedRipples empty;
  const auto IDENTITY =
      radial_ripples(math::X_AXIS, empty, true, &tables, 1.0f);
  HS_EXPECT_EQ(IDENTITY.sphere.x, 1.0f);
  HS_EXPECT_EQ(IDENTITY.path_length, 0.0f);
}

/** @brief Radial rejection resolves the whole original input exactly once. */
inline void test_pullback_radial_ripple_lazy_fallback() {
  using namespace Pullback::SphereAnimation;
  static RadialTables tables;
  tables.enabled = true;
  tables.count = 2;
  for (auto &group : tables.groups) {
    group.cos_threshold_min = 1.0f;
    group.cos_threshold_max = -1.0f;
    group.previous.fill({0.2f, 0.2f});
    group.current.fill({0.2f, 0.2f});
  }
  tables.groups[0].center = math::Y_AXIS;
  tables.groups[1].center = math::Z_AXIS;
  const Pullback::SurfaceResult EXPECTED{math::Vector(0.2f, 0.3f, 0.4f), 0.75f};
  int calls = 0;
  math::Vector input = math::X_AXIS;
  const auto FALLBACK = [&](const math::Vector &original) {
    ++calls;
    HS_EXPECT_EQ(original, input);
    return EXPECTED;
  };
  const auto CHECK_FALLBACK = [&](const RadialTables *view) {
    calls = 0;
    const auto RESULT = radial_ripples(input, true, view, 0.5f, FALLBACK);
    HS_EXPECT_EQ(calls, 1);
    HS_EXPECT_EQ(RESULT.sphere, EXPECTED.sphere);
    HS_EXPECT_EQ(RESULT.path_length, EXPECTED.path_length);
  };
  CHECK_FALLBACK(nullptr);
  tables.enabled = false;
  CHECK_FALLBACK(&tables);
  tables.enabled = true;
  tables.groups[1].current.fill({0.2f, -1.0f});
  CHECK_FALLBACK(&tables);
  tables.count = 1;
  input = math::Y_AXIS;
  CHECK_FALLBACK(&tables);
  input = 2.0f * math::X_AXIS;
  CHECK_FALLBACK(&tables);
  input = math::X_AXIS;
  calls = 0;
  const auto RESULT = radial_ripples(input, true, &tables, 0.5f, FALLBACK);
  HS_EXPECT_EQ(calls, 0);
  HS_EXPECT_GT((RESULT.sphere - input).length(), 0.1f);
}

template <typename Value> struct RippleTestResource {
  Value value;
  template <auto Key> const Value &get() const { return value; }
};
struct RippleTestClock : Pullback::SphereAnimation::RippleState {
  const Pullback::SphereAnimation::RadialTables *radial_tables = nullptr;
  float radial_fraction = -1.0f;
};
struct RippleTestFrame {
  RippleTestResource<Pullback::SphereAnimation::BurstRippleParams> params;
  RippleTestResource<RippleTestClock> resources;
};
struct RippleTestBinding {
  using FrameState = RippleTestFrame;
};
using RippleTestProvider = Pullback::SurfaceProvider<
    RippleTestBinding, Pullback::SphereAnimation::BurstRippleParams, true>;
struct CustomRippleTestProvider : RippleTestProvider {
  static inline int calls = 0;
  static auto prepare(const FrameState &frame) {
    ++calls;
    auto prepared = RippleTestProvider::prepare(frame);
    prepared.count = 0;
    return prepared;
  }
};
struct GenericRippleTestProvider {
  using Binding = RippleTestBinding;
  using FrameState = RippleTestFrame;
  static auto prepare(const FrameState &frame) {
    return RippleTestProvider::prepare(frame);
  }
  static bool path_length_required(const FrameState &) { return true; }
};

/** @brief Only owned usable radial views defer exact column preparation. */
inline void test_pullback_ripple_lazy_preparation() {
  using namespace Pullback::SphereAnimation;
  using Policy = Pullback::Surface::BurstRipple<RippleTestProvider>;
  using Custom = Pullback::Surface::BurstRipple<CustomRippleTestProvider>;
  using Generic = Pullback::Surface::BurstRipple<GenericRippleTestProvider>;
  static_assert(Policy::LAZY_RADIAL_PREPARE);
  static_assert(!Custom::LAZY_RADIAL_PREPARE);
  static_assert(!Generic::LAZY_RADIAL_PREPARE);
  static_assert(std::is_same_v<Custom::Prepared, PreparedRipples>);
  static_assert(std::is_same_v<Generic::Prepared, PreparedRipples>);
  RippleTestFrame frame;
  auto &params = frame.params.value;
  auto &clock = frame.resources.value;
  params.period = 80.0f;
  params.random_center = 0;
  params.center_polar = 0.0f;
  clock.burst_count = 1;
  clock.waves[0] = {40.0f, 0};
  static RadialTables tables;
  prepare_radial_tables(params, clock, tables);
  clock.radial_tables = &tables;
  clock.radial_fraction = 0.5f;
  const auto PREPARED = Policy::prepare(frame);
  HS_EXPECT_FALSE(PREPARED.has_value());
  const auto EXACT = RippleTestProvider::prepare(frame);
  const auto EXPECTED =
      radial_ripples(math::X_AXIS, EXACT, true, &tables, 0.5f);
  const auto RESULT = Policy::apply(math::X_AXIS, frame, PREPARED);
  HS_EXPECT_EQ(RESULT.sphere, EXPECTED.sphere);
  HS_EXPECT_EQ(RESULT.path_length, EXPECTED.path_length);
  const auto CUSTOM_INTACT = Custom::prepare(frame);
  HS_EXPECT_EQ(Custom::apply(math::X_AXIS, frame, CUSTOM_INTACT).sphere,
               math::X_AXIS);
  tables.groups[0].current.fill({0.0f, -1.0f});
  const auto FALLBACK = Policy::apply(math::X_AXIS, frame, PREPARED);
  const auto GROUPED = grouped_ripples(math::X_AXIS, EXACT, true);
  HS_EXPECT_EQ(FALLBACK.sphere, GROUPED.sphere);
  HS_EXPECT_EQ(FALLBACK.path_length, GROUPED.path_length);
  CustomRippleTestProvider::calls = 0;
  const auto CUSTOM = Custom::prepare(frame);
  HS_EXPECT_EQ(CustomRippleTestProvider::calls, 1);
  HS_EXPECT_EQ(CUSTOM.count, 0);
  HS_EXPECT_EQ(Custom::apply(math::X_AXIS, frame, CUSTOM).sphere, math::X_AXIS);
  HS_EXPECT_EQ(CustomRippleTestProvider::calls, 1);
  HS_EXPECT_EQ(Generic::prepare(frame).count, EXACT.count);
  HS_EXPECT_EQ(Generic::apply(math::X_AXIS, frame, EXACT).sphere,
               GROUPED.sphere);
  {
    clock.radial_fraction = -1.0f;
    const auto ENDPOINT = Policy::prepare(frame);
    HS_EXPECT_TRUE(ENDPOINT.has_value());
    HS_EXPECT_EQ(ENDPOINT->count, EXACT.count);
    HS_EXPECT_EQ(Policy::apply(math::X_AXIS, frame, ENDPOINT).sphere,
                 GROUPED.sphere);
  }
  clock.radial_fraction = 0.5f;
  tables.enabled = false;
  HS_EXPECT_TRUE(Policy::prepare(frame).has_value());
  tables.enabled = true;
  clock.radial_tables = nullptr;
  HS_EXPECT_TRUE(Policy::prepare(frame).has_value());
}
/** @brief Cached ripple geometry stays subpixel across column times. */
inline void test_pullback_radial_ripple_geometry() {
  using namespace Pullback::SphereAnimation;
  BurstRippleParams params;
  params.period = 91.48529f;
  params.decay = 0.2f;
  params.burst_rate = 1.0f / 128.0f;
  RippleState state;
  state.seed = 149502542u;
  static RadialTables tables;
  for (int frame = 0; frame < 300; ++frame) {
    advance_ripple(state, params);
    prepare_radial_tables(params, state, tables);
    for (float fraction : {0.0f, 0.5f, 1.0f}) {
      auto sampled = state;
      for (auto &wave : sampled.waves)
        wave.age -= 1.0f - fraction;
      auto prepared = prepare_ripples(params, sampled);
      for (int i = 0; i < 32; ++i) {
        const float Y = -1.0f + (2.0f * i + 1.0f) / 32.0f;
        const float AZIMUTH = i * 2.39996323f;
        const float RADIUS = sqrtf(1.0f - Y * Y);
        const math::Vector VIEW(RADIUS * cosf(AZIMUTH), Y,
                                RADIUS * sinf(AZIMUTH));
        const auto EXPECTED = grouped_ripples(VIEW, prepared, true);
        const auto RESULT =
            radial_ripples(VIEW, prepared, true, &tables, fraction);
        HS_EXPECT_TRUE(std::isfinite(RESULT.path_length));
        HS_EXPECT_GE(RESULT.path_length, 0.0f);
        HS_EXPECT_NEAR(math::dot(RESULT.sphere, RESULT.sphere), 1.0f, 1e-5f);
        HS_EXPECT_LE((RESULT.sphere - EXPECTED.sphere).length(), 0.005f);
        HS_EXPECT_LE(fabsf(RESULT.path_length - EXPECTED.path_length), 0.02f);
      }
    }
  }
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
  const Pullback::Color::PhaseModulationState *phase = nullptr;
  const PaletteLabLut *lab = nullptr;
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
  static constexpr bool FUSED_TEMPORAL_COLOR = true;
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

struct VariableAlphaPalette {
  Color4 get(float value) const {
    const float INDEX = value * (BakedPalette::LUT_SIZE - 1);
    return Color4(Pixel(1000, 2000, 3000),
                  INDEX > 0.5f && INDEX < 1.5f ? 0.5f : 1.0f);
  }
};

struct LabRampPalette {
  Color4 get(float value) const {
    return Color4(Pixel(static_cast<uint16_t>(value * 50000),
                        static_cast<uint16_t>(value * value * 30000),
                        static_cast<uint16_t>((1.0f - value) * 20000)),
                  1.0f);
  }
};

/** @brief Palette Lab nodes convert RGB-first blends and reject varying alpha. */
inline void test_pullback_palette_lab_nodes() {
  std::array<uint8_t, 5 * BakedPalette::required_arena_bytes()> buffer{};
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage a, b, variable;
  a.bake(arena, ConstantPalette{Color4(Pixel(0, 0, 0), 1.0f)});
  b.bake(arena, ConstantPalette{Color4(Pixel(65535, 65535, 65535), 1.0f)});
  variable.bake(arena, VariableAlphaPalette{});
  HS_EXPECT_LT(variable.view().get_alpha(1.0f / 255.0f), 1.0f);
  for (int index = 0; index < PaletteLabLut::VALUES; ++index)
    HS_EXPECT_EQ(
        variable.view().get_alpha(float(index) / (PaletteLabLut::VALUES - 1)),
        1.0f);
  PaletteLabLut cache;
  cache.bake(a.view(), b.view());
  HS_EXPECT_TRUE(cache.active);
  const auto MID = PaletteLabLut::sample(cache.view(0.5f), 0.37f);
  HS_EXPECT_NEAR(MID.L, std::cbrt(32768.0 / 65535.0), 3e-5);
  HS_EXPECT_GT(MID.L, 0.75f);
  HS_EXPECT_NEAR(MID.a, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(MID.b, 0.0f, 1e-6f);
  for (int row = 0; row < PaletteLabLut::TIMES; ++row)
    for (int index = 0; index < PaletteLabLut::VALUES; ++index) {
      const float FRACTION = float(row) / (PaletteLabLut::TIMES - 1);
      const float VALUE = float(index) / (PaletteLabLut::VALUES - 1);
      const Pixel RGB = a.view().get_color(VALUE).lerp16(
          b.view().get_color(VALUE), frac_to_q16(FRACTION));
      const auto LINEAR = pixel_to_linrgb(RGB);
      const auto EXPECTED =
          linear_rgb_to_oklab_fast(LINEAR.r, LINEAR.g, LINEAR.b);
      const auto ACTUAL = PaletteLabLut::sample(cache.view(FRACTION), VALUE);
      HS_EXPECT_NEAR(ACTUAL.L, EXPECTED.L, 1e-7f);
      HS_EXPECT_NEAR(ACTUAL.a, EXPECTED.a, 1e-7f);
      HS_EXPECT_NEAR(ACTUAL.b, EXPECTED.b, 1e-7f);
    }
  const auto OLD_VIEW = cache.view(0.25f);
  PaletteLabLut next;
  next.bake(b.view(), a.view());
  HS_EXPECT_NEAR(PaletteLabLut::sample(OLD_VIEW, 0.37f).L, MID.L * 0.5f, 1e-7f);
  BakedPaletteStorage ramp, translucent;
  ramp.bake(arena, LabRampPalette{});
  translucent.bake(arena, ConstantPalette{Color4(Pixel(123, 456, 789), 0.6f)});
  next.bake(ramp.view(), b.view());
  for (float value : {0.183f, 0.677f}) {
    const float POSITION = value * 63.0f;
    const int INDEX = static_cast<int>(POSITION);
    const float WEIGHT = POSITION - INDEX;
    const auto node = [&](int index, float fraction) {
      const float COORDINATE = float(index) / 63.0f;
      const Pixel RGB =
          ramp.view()
              .get_color(COORDINATE)
              .lerp16(b.view().get_color(COORDINATE), frac_to_q16(fraction));
      const LinRGB LINEAR = pixel_to_linrgb(RGB);
      return linear_rgb_to_oklab_fast(LINEAR.r, LINEAR.g, LINEAR.b);
    };
    const auto A = node(INDEX, 0.0f), B = node(INDEX + 1, 0.0f);
    const auto C = node(INDEX, 0.5f), D = node(INDEX + 1, 0.5f);
    const auto ACTUAL = PaletteLabLut::sample(next.view(0.17f), value);
    HS_EXPECT_NEAR(
        ACTUAL.L,
        hs::lerp(hs::lerp(A.L, B.L, WEIGHT), hs::lerp(C.L, D.L, WEIGHT), 0.34f),
        1e-7f);
    HS_EXPECT_NEAR(
        ACTUAL.a,
        hs::lerp(hs::lerp(A.a, B.a, WEIGHT), hs::lerp(C.a, D.a, WEIGHT), 0.34f),
        1e-7f);
    HS_EXPECT_NEAR(
        ACTUAL.b,
        hs::lerp(hs::lerp(A.b, B.b, WEIGHT), hs::lerp(C.b, D.b, WEIGHT), 0.34f),
        1e-7f);
  }
  next.bake(ramp.view(), translucent.view());
  HS_EXPECT_FALSE(next.active);
  cache.bake(a.view(), variable.view());
  HS_EXPECT_FALSE(cache.active);
  HS_EXPECT_TRUE(cache.view(0.5f).low == nullptr);
}

/** @brief Cached neutral adjustments materialize their base with exact alpha. */
inline void test_pullback_cached_lab_finish() {
  using namespace Pullback;
  ColorSample sample{Color4(Pixel(0, 0, 0), 0.37f),
                     {0.5f, 0.7f, math::Y_AXIS, 2.0f}};
  sample.cached_lab = {0.7f, 0.02f, 0.03f};
  sample.has_cached_lab = true;
  const auto EXPECTED =
      linrgb_to_pixel(oklab_to_linear_rgb_lut_gamut(sample.cached_lab));
  const auto ACTUAL = finish_color_sample(sample);
  HS_EXPECT_EQ(ACTUAL.color, EXPECTED);
  HS_EXPECT_EQ(ACTUAL.alpha, 0.37f);
  HS_EXPECT_GT(ACTUAL.color.r, 0);
  sample.brightness_gain = 0.0f;
  HS_EXPECT_EQ(finish_color_sample(sample).color, Pixel(0, 0, 0));
  HS_EXPECT_EQ(finish_color_sample(sample).alpha, 0.37f);
  HS_EXPECT_TRUE(blend_color_samples(sample, sample, 0.0f).has_cached_lab);
  HS_EXPECT_TRUE(blend_color_samples(sample, sample, 1.0f).has_cached_lab);
}
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

struct PhaseLayerProvider {
  using Binding = PaletteBinding;
  using FrameState = PaletteFrame;
  static Pullback::Color::PhaseModulationState
  prepare(const PaletteFrame &frame) {
    return *frame.phase;
  }
};
struct CountingPhasePolicy : Pullback::Color::Modulate<PhaseLayerProvider> {
  using Base = Pullback::Color::Modulate<PhaseLayerProvider>;
  inline static int phase_calls = 0;
  inline static int endpoint_calls = 0;
  static void apply_endpoint_in_place(Pullback::ColorSample &sample,
                                      const PaletteFrame &frame,
                                      const Prepared &prepared, bool previous) {
    ++endpoint_calls;
    Base::apply_endpoint_in_place(sample, frame, prepared, previous);
  }
  static void apply_phase_in_place(Pullback::ColorSample &sample,
                                   const PaletteFrame &frame,
                                   const Prepared &prepared) {
    ++phase_calls;
    Base::apply_phase_in_place(sample, frame, prepared);
  }
};

struct LabPaletteProvider : TemporalPaletteProvider {
  static const PaletteLabLut *palette_lab_lut(const PaletteFrame &frame) {
    return frame.lab;
  }
};
struct CountingLabPalette
    : Pullback::Color::GeneratedPalette<LabPaletteProvider> {
  using Base = Pullback::Color::GeneratedPalette<LabPaletteProvider>;
  using SharedBaseOwner = CountingLabPalette;
  inline static int shared_calls = 0;
  inline static int pair_calls = 0;
  static bool apply_shared_base(const Pullback::FieldSample &input,
                                const PaletteFrame &frame,
                                const Prepared &prepared,
                                Pullback::ColorSample &output) {
    ++shared_calls;
    return Base::apply_shared_base(input, frame, prepared, output);
  }
  static auto apply_pair(const Pullback::FieldSample &input,
                         const PaletteFrame &frame, const Prepared &prepared) {
    ++pair_calls;
    return Base::apply_pair(input, frame, prepared);
  }
};

struct ShiftedLabPalette : CountingLabPalette {
  static auto apply_pair(const Pullback::FieldSample &input,
                         const PaletteFrame &frame, const Prepared &prepared) {
    auto pair = CountingLabPalette::apply_pair(input, frame, prepared);
    pair.first.color.color.r += 123;
    pair.second.color.color.r += 123;
    return pair;
  }
};

/** @brief Consecutive analogous cycle bakes meet cached channel bounds. */
inline void test_pullback_palette_lab_qualification() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  std::array<uint8_t, PaletteCycler::generated_arena_bytes() +
                          2 * BakedPalette::required_arena_bytes()>
      buffer{};
  Arena arena(buffer.data(), buffer.size());
  uint32_t hue = 0;
  PaletteCycler cycler;
  cycler.init_generated(
      arena,
      +[](void *context, uint32_t sequence, GenerativePalette &out) {
        GeneratedPaletteBank::next_palette(*static_cast<uint32_t *>(context),
                                           sequence, PaletteHarmony::ANALOGOUS,
                                           0.84f, out);
      },
      &hue, GeneratedPaletteBank::DWELL_FRAMES,
      GeneratedPaletteBank::FADE_FRAMES, math::ease_in_out_sin, nullptr, 2);
  for (int frame = 0; frame < GeneratedPaletteBank::FADE_FRAMES / 2; ++frame)
    cycler.step();
  BakedPaletteStorage from, to;
  from.clone_from(cycler.palette(), arena);
  cycler.step();
  cycler.step();
  to.clone_from(cycler.palette(), arena);
  PaletteLabLut cache;
  cache.bake(from.view(), to.view());
  HS_EXPECT_TRUE(cache.active);
  PaletteFrame frame;
  frame.palette = &to.view();
  frame.temporal.palette = &from.view();
  frame.lab = &cache;
  uint32_t maximum = 0;
  uint64_t total = 0;
  size_t channels = 0;
  constexpr int VALUES = 2048;
  constexpr int FRACTIONS = 128;
  for (int time = 0; time < FRACTIONS; ++time) {
    frame.temporal.fraction = (static_cast<float>(time) + 0.5f) / FRACTIONS;
    const auto PREPARED = Policy::prepare(frame);
    for (int index = 0; index < VALUES; ++index) {
      const FieldSample FIELD{(static_cast<float>(index) + 0.5f) / VALUES,
                              0.75f, math::Y_AXIS, 0.0f};
      ColorSample cached;
      HS_EXPECT_TRUE(Policy::apply_shared_base(FIELD, frame, PREPARED, cached));
      const auto PAIR = Policy::apply_pair(FIELD, frame, PREPARED);
      const auto EXPECTED = finish_color_sample(blend_color_samples(
          PAIR.first, PAIR.second, frame.temporal.fraction));
      const auto ACTUAL = finish_color_sample(cached);
      HS_EXPECT_NEAR(ACTUAL.alpha, EXPECTED.alpha, 1e-7f);
      for (const uint32_t ERROR :
           {channel_gap(ACTUAL.color.r, EXPECTED.color.r),
            channel_gap(ACTUAL.color.g, EXPECTED.color.g),
            channel_gap(ACTUAL.color.b, EXPECTED.color.b)}) {
        maximum = std::max(maximum, ERROR);
        total += ERROR;
        ++channels;
      }
    }
  }
  const double MEAN = static_cast<double>(total) / channels;
  std::printf("  palette lab cache max=%u mean=%.3f channel codes\n", maximum,
              MEAN);
  HS_EXPECT_GT(maximum, 0u);
  HS_EXPECT_LE(static_cast<float>(maximum), Policy::METRICS[0].limit);
  HS_EXPECT_LE(MEAN, Policy::METRICS[1].limit);
}

/** @brief Shared cached bases skip both palette endpoints and preserve alpha. */
inline void test_pullback_shared_lab_dispatch() {
  using namespace Pullback;
  using Pipeline = Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                                      Stage::Colorize<CountingLabPalette>,
                                      Stage::Modulate<CountingPhasePolicy>>;
  using Core =
      Detail::PipelineCore<true, PaletteBinding, typename Pipeline::NodeList>;
  using Custom = Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                                    Stage::Colorize<ShiftedLabPalette>,
                                    Stage::Modulate<CountingPhasePolicy>>;
  using CustomCore =
      Detail::PipelineCore<true, PaletteBinding, typename Custom::NodeList>;
  std::array<uint8_t, 3 * BakedPalette::required_arena_bytes()> buffer{};
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage a, b;
  a.bake(arena, ConstantPalette{Color4(Pixel(1000, 2000, 3000), 1.0f)});
  b.bake(arena, ConstantPalette{Color4(Pixel(9000, 6000, 1000), 1.0f)});
  PaletteLabLut cache;
  cache.bake(a.view(), b.view());
  Color::PhaseModulationState layer;
  layer.source = Pattern::SourceKind::DISPLACEMENT;
  layer.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 0.2f};
  layer.enabled = true;
  layer.ready = true;
  PaletteFrame frame;
  frame.palette = &b.view();
  frame.temporal.palette = &a.view();
  frame.phase = &layer;
  frame.lab = &cache;
  frame.temporal.fraction = 0.37f;
  const auto DIRECT_PREPARED = CountingLabPalette::Base::prepare(frame);
  const FieldSample DIRECT_FIELD{0.5f, 0.5f, math::Y_AXIS, 0.7f};
  ColorSample dirty{Color4(Pixel(123, 456, 789), 0.9f),
                    {0.1f, 0.2f, math::X_AXIS, 0.3f},
                    reinterpret_cast<const Pixel *>(&cache),
                    0.4f,
                    0.5f,
                    0.6f,
                    0.7f,
                    {0.8f, 0.9f, 1.0f},
                    false};
  HS_EXPECT_TRUE(CountingLabPalette::Base::apply_shared_base(
      DIRECT_FIELD, frame, DIRECT_PREPARED, dirty));
  HS_EXPECT_EQ(dirty.color.color, Pixel(0, 0, 0));
  const auto DIRECT_PAIR = CountingLabPalette::Base::apply_pair(
      DIRECT_FIELD, frame, DIRECT_PREPARED);
  HS_EXPECT_EQ(dirty.color.alpha, DIRECT_PAIR.first.color.alpha);
  HS_EXPECT_EQ(dirty.field.value, DIRECT_FIELD.value);
  HS_EXPECT_EQ(dirty.field.coverage, DIRECT_FIELD.coverage);
  HS_EXPECT_EQ(dirty.field.sphere, DIRECT_FIELD.sphere);
  HS_EXPECT_EQ(dirty.field.path_length, DIRECT_FIELD.path_length);
  HS_EXPECT_TRUE(dirty.hue_rotation == nullptr);
  HS_EXPECT_EQ(dirty.palette_coordinate, 0.5f);
  HS_EXPECT_EQ(dirty.pending_hue, 0.0f);
  HS_EXPECT_EQ(dirty.brightness_gain, 1.0f);
  HS_EXPECT_EQ(dirty.pending_lift, 0.0f);
  const auto EXPECTED_LAB = PaletteLabLut::sample(cache.view(0.37f), 0.5f);
  HS_EXPECT_EQ(dirty.cached_lab.L, EXPECTED_LAB.L);
  HS_EXPECT_EQ(dirty.cached_lab.a, EXPECTED_LAB.a);
  HS_EXPECT_EQ(dirty.cached_lab.b, EXPECTED_LAB.b);
  HS_EXPECT_TRUE(dirty.has_cached_lab);
  for (bool ready : {false, true})
    for (bool active : {false, true})
      for (float fraction : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f}) {
        layer.ready = ready;
        frame.lab = active ? &cache : nullptr;
        frame.temporal.fraction = fraction;
        CountingLabPalette::shared_calls = 0;
        CountingLabPalette::pair_calls = 0;
        const auto PREPARED = Pipeline::prepare_stages(frame);
        const auto ACTUAL = Core::template run_stage<0>(
            SphereSample{math::Y_AXIS, 0.7f}, frame, PREPARED);
        const bool INTERIOR = fraction > 0.0f && fraction < 1.0f;
        HS_EXPECT_EQ(ACTUAL.has_cached_lab, active && ready && INTERIOR);
        HS_EXPECT_EQ(CountingLabPalette::shared_calls,
                     INTERIOR && ready ? 1 : 0);
        HS_EXPECT_EQ(CountingLabPalette::pair_calls,
                     INTERIOR && !(active && ready) ? 1 : 0);
        const FieldSample FIELD{0.5f, 0.5f, math::Y_AXIS, 0.7f};
        const auto PAIR = CountingLabPalette::Base::apply_pair(
            FIELD, frame, std::get<1>(PREPARED));
        HS_EXPECT_EQ(ACTUAL.color.alpha,
                     fraction <= 0.0f ? PAIR.first.color.alpha
                     : fraction >= 1.0f
                         ? PAIR.second.color.alpha
                         : hs::lerp(PAIR.first.color.alpha,
                                    PAIR.second.color.alpha, fraction));
        HS_EXPECT_NEAR(ACTUAL.pending_hue, 0.14f, 1e-7f);
        if (!INTERIOR)
          HS_EXPECT_EQ(ACTUAL.color.color, fraction == 0.0f
                                               ? PAIR.first.color.color
                                               : PAIR.second.color.color);
      }
  layer.ready = true;
  frame.lab = &cache;
  frame.temporal.fraction = 0.5f;
  CountingLabPalette::shared_calls = 0;
  CountingLabPalette::pair_calls = 0;
  const auto CUSTOM = CustomCore::template run_stage<0>(
      SphereSample{math::Y_AXIS, 0.7f}, frame, Custom::prepare_stages(frame));
  HS_EXPECT_FALSE(CUSTOM.has_cached_lab);
  HS_EXPECT_EQ(CountingLabPalette::shared_calls, 0);
  HS_EXPECT_EQ(CountingLabPalette::pair_calls, 1);
  HS_EXPECT_EQ(CUSTOM.color.color.r, 5123);
  BakedPaletteStorage translucent;
  translucent.bake(arena,
                   ConstantPalette{Color4(Pixel(4000, 5000, 6000), 0.6f)});
  cache.bake(a.view(), translucent.view());
  frame.palette = &translucent.view();
  const auto PREPARED = Pipeline::prepare_stages(frame);
  const auto FALLBACK = Core::template run_stage<0>(
      SphereSample{math::Y_AXIS, 0.7f}, frame, PREPARED);
  const auto PAIR = CountingLabPalette::Base::apply_pair(
      FieldSample{0.5f, 0.5f, math::Y_AXIS, 0.7f}, frame,
      std::get<1>(PREPARED));
  HS_EXPECT_FALSE(FALLBACK.has_cached_lab);
  HS_EXPECT_EQ(FALLBACK.color.alpha,
               hs::lerp(PAIR.first.color.alpha, PAIR.second.color.alpha, 0.5f));
}

/** @brief Shared phase tails run once; incomplete or unready tails keep endpoints. */
inline void test_pullback_phase_tail_dispatch() {
  using namespace Pullback;
  using Phase = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage,
      Stage::Colorize<Color::GeneratedPalette<TemporalPaletteProvider>>,
      Stage::Modulate<CountingPhasePolicy>,
      Stage::Modulate<CountingPhasePolicy>>;
  using Fallback = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage,
      Stage::Colorize<Color::GeneratedPalette<TemporalPaletteProvider>>,
      Stage::Modulate<CountingPhasePolicy>,
      Stage::Modulate<CountingPhasePolicy>,
      Stage::Modulate<Color::Modulate<PaletteLayerProvider>>>;
  using PhaseCore =
      Detail::PipelineCore<true, PaletteBinding, typename Phase::NodeList>;
  using FallbackCore =
      Detail::PipelineCore<true, PaletteBinding, typename Fallback::NodeList>;
  std::array<uint8_t, 2 * BakedPalette::required_arena_bytes()> buffer{};
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage a, b;
  a.bake(arena, ConstantPalette{Color4(Pixel(1000, 2000, 3000), 0.2f)});
  b.bake(arena, ConstantPalette{Color4(Pixel(9000, 6000, 1000), 0.6f)});
  Color::PhaseModulationState layer;
  layer.source = Pattern::SourceKind::DISPLACEMENT;
  layer.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 0.2f};
  layer.enabled = true;
  PaletteFrame frame;
  frame.palette = &b.view();
  frame.temporal.palette = &a.view();
  frame.phase = &layer;
  frame.mechanic.amount = 0.0f;
  for (bool ready : {false, true})
    for (float fraction : {0.0f, 0.25f, 0.75f, 1.0f}) {
      layer.ready = ready;
      frame.temporal.fraction = fraction;
      CountingPhasePolicy::phase_calls = 0;
      CountingPhasePolicy::endpoint_calls = 0;
      const auto prepared = Phase::prepare_stages(frame);
      const auto SAMPLE = PhaseCore::template run_stage<0>(
          SphereSample{math::Y_AXIS, 0.7f}, frame, prepared);
      const bool SHARED = ready && fraction > 0.0f && fraction < 1.0f;
      HS_EXPECT_EQ(CountingPhasePolicy::phase_calls, SHARED ? 2 : 0);
      HS_EXPECT_EQ(CountingPhasePolicy::endpoint_calls,
                   SHARED ? 0 : (fraction > 0.0f && fraction < 1.0f ? 4 : 2));
      HS_EXPECT_NEAR(SAMPLE.pending_hue, 0.28f, 1e-7f);
      HS_EXPECT_NEAR(SAMPLE.color.alpha, 0.04f + 0.08f * fraction, 1e-7f);
      HS_EXPECT_NEAR(SAMPLE.color.color.r, 1000 + 8000 * fraction, 1.0f);
      CountingPhasePolicy::phase_calls = 0;
      CountingPhasePolicy::endpoint_calls = 0;
      const auto OTHER = FallbackCore::template run_stage<0>(
          SphereSample{math::Y_AXIS, 0.7f}, frame,
          Fallback::prepare_stages(frame));
      HS_EXPECT_EQ(CountingPhasePolicy::phase_calls, 0);
      HS_EXPECT_EQ(CountingPhasePolicy::endpoint_calls,
                   fraction > 0.0f && fraction < 1.0f ? 4 : 2);
      HS_EXPECT_EQ(OTHER.color.color, SAMPLE.color.color);
      HS_EXPECT_EQ(OTHER.color.alpha, SAMPLE.color.alpha);
      HS_EXPECT_EQ(OTHER.pending_hue, SAMPLE.pending_hue);
    }
}

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
 * @param gain Final linear RGB gain.
 * @return The transformed color. */
inline Color4 composed_modulation_reference(Color4 base, float turns,
                                            float lift, float gain = 1.0f) {
  const LinRGB rgb = pixel_to_linrgb(base.color);
  OKLab lab = linear_rgb_to_oklab_fast(rgb.r, rgb.g, rgb.b);
  lab.L += (1.0f - lab.L) * lift;
  float cosine, sine;
  hue_sincos(turns, cosine, sine);
  const OKLab rotated{lab.L, lab.a * cosine - lab.b * sine,
                      lab.a * sine + lab.b * cosine};
  LinRGB output = oklab_to_linear_rgb_lut_gamut(rotated);
  if (gain != 1.0f) {
    output.r = hs::clamp(output.r, 0.0f, 1.0f) * gain;
    output.g = hs::clamp(output.g, 0.0f, 1.0f) * gain;
    output.b = hs::clamp(output.b, 0.0f, 1.0f) * gain;
  }
  return Color4(linrgb_to_pixel(output), base.alpha);
}

/** @brief Perceptual finishing applies brightness before RGB16 quantization. */
inline void test_pullback_finish_brightness_quantization() {
  using namespace Pullback;
  size_t fused_differences = 0;
  for (uint16_t r : {0, 1, 2, 1000, 32768, 65535})
    for (uint16_t g : {0, 1, 2, 1000, 32768, 65535})
      for (uint16_t b : {0, 1, 2, 1000, 32768, 65535})
        for (float hue : {0.0f, 0.17f, 0.75f})
          for (float lift : {0.0f, 0.25f, 0.9f}) {
            if (hue == 0.0f && lift == 0.0f)
              continue;
            const Color4 BASE(Pixel(r, g, b), 0.37f);
            const auto LEGACY = composed_modulation_reference(BASE, hue, lift);
            for (float gain : {-1.0f, 0.0f, 0.001f, 0.25f, 1.0f, 2.0f}) {
              ColorSample sample{BASE, {0.5f, 0.7f, math::Y_AXIS, 2.0f}};
              sample.pending_hue = hue;
              sample.pending_lift = lift;
              sample.brightness_gain = gain;
              const auto ACTUAL = finish_color_sample(sample);
              const auto EXPECTED =
                  composed_modulation_reference(BASE, hue, lift, gain);
              HS_EXPECT_EQ(ACTUAL.color, EXPECTED.color);
              HS_EXPECT_EQ(ACTUAL.alpha, BASE.alpha);
              const Pixel OLD = LEGACY.color * gain;
              const int GAP = std::max({std::abs(int(ACTUAL.color.r) - OLD.r),
                                        std::abs(int(ACTUAL.color.g) - OLD.g),
                                        std::abs(int(ACTUAL.color.b) - OLD.b)});
              HS_EXPECT_LE(GAP, gain <= 1.0f ? 1 : 2);
              if (gain == 0.0f)
                HS_EXPECT_EQ(ACTUAL.color, Pixel(0, 0, 0));
              if (gain == 1.0f)
                HS_EXPECT_EQ(ACTUAL.color, LEGACY.color);
              fused_differences += ACTUAL.color != OLD;
            }
          }
  HS_EXPECT_GT(fused_differences, size_t{0});
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
  const auto DIMMED = Color::finish(actual).color;
  const auto DIMMED_REFERENCE =
      composed_modulation_reference(BASE, 0.0f, 0.4f, 0.25f).color;
  HS_EXPECT_EQ(DIMMED, DIMMED_REFERENCE);
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
    HS_EXPECT_EQ(
        MASKED.color,
        composed_modulation_reference(INPUT.color, 0.375f, 0.55f, 0.25f).color);
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

/** @brief In-place temporal blending preserves current sampling coordinates. */
inline void test_pullback_temporal_deferred_in_place() {
  using namespace Pullback;
  std::array<Pixel, 2> hue{};
  ColorSample a{Color4(Pixel(1000, 2000, 4000), 0.25f),
                {0.1f, 0.2f, math::X_AXIS, 0.3f}};
  ColorSample b{Color4(Pixel(9000, 6000, 12000), 0.75f),
                {0.6f, 0.7f, math::Y_AXIS, 0.8f}};
  a.hue_rotation = hue.data();
  b.hue_rotation = hue.data();
  a.palette_coordinate = 0.2f;
  b.palette_coordinate = 0.9f;
  a.pending_hue = 0.875f;
  b.pending_hue = 0.125f;
  a.pending_lift = 0.25f;
  b.pending_lift = 0.75f;
  a.brightness_gain = 0.5f;
  b.brightness_gain = 0.75f;
  blend_color_samples_in_place(a, b, 0.5f);
  HS_EXPECT_EQ(b.color.color, Pixel(5000, 4000, 8000));
  HS_EXPECT_EQ(b.color.alpha, 0.5f);
  HS_EXPECT_EQ(b.pending_hue, 0.0f);
  HS_EXPECT_EQ(b.pending_lift, 0.5f);
  HS_EXPECT_EQ(b.brightness_gain, 0.625f);
  HS_EXPECT_EQ(b.field.value, 0.6f);
  HS_EXPECT_EQ(b.field.coverage, 0.7f);
  HS_EXPECT_EQ(b.field.sphere, math::Y_AXIS);
  HS_EXPECT_EQ(b.field.path_length, 0.8f);
  HS_EXPECT_EQ(b.palette_coordinate, 0.9f);
  HS_EXPECT_TRUE(b.hue_rotation == nullptr);
  HS_EXPECT_EQ(a.color.color, Pixel(1000, 2000, 4000));
  HS_EXPECT_EQ(a.pending_hue, 0.875f);
  HS_EXPECT_TRUE(a.hue_rotation == hue.data());

  b = a;
  b.color.color.r = 1001;
  blend_color_samples_in_place(a, b, 0.000001f);
  HS_EXPECT_EQ(b.color.color, a.color.color);
  HS_EXPECT_TRUE(b.hue_rotation == nullptr);
  a.cached_lab = {0.4f, 0.1f, -0.2f};
  a.has_cached_lab = true;
  blend_color_samples_in_place(a, b, 0.0f);
  HS_EXPECT_TRUE(b.has_cached_lab);
  HS_EXPECT_EQ(b.cached_lab.L, 0.4f);
  HS_EXPECT_EQ(b.cached_lab.a, 0.1f);
  HS_EXPECT_EQ(b.cached_lab.b, -0.2f);
  HS_EXPECT_TRUE(b.hue_rotation == hue.data());
  blend_color_samples_in_place(a, b, 1.0f);
  HS_EXPECT_TRUE(b.has_cached_lab);
}

/** @brief A fused final layer matches the canonical temporal endpoint pipeline. */
inline void test_pullback_fused_modulation_endpoints() {
  using namespace Pullback;
  std::array<uint8_t, 2 * BakedPalette::required_arena_bytes()> bytes{};
  Arena arena(bytes.data(), bytes.size());
  BakedPaletteStorage a, b;
  a.bake(arena, CoordinatePalette{0.2f});
  b.bake(arena, ConstantPalette{Color4(Pixel(32000, 44000, 13000), 0.8f)});
  PaletteFrame frame;
  frame.palette = &b.view();
  frame.temporal.palette = &a.view();
  using Policy = Color::GeneratedPalette<TemporalPaletteProvider>;
  struct CanonicalPalette : Policy {};
  using Actual = TemporalPalettePipeline;
  using Canonical = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage, Stage::Colorize<CanonicalPalette>,
      Stage::Modulate<Color::Modulate<PaletteLayerProvider>>>;
  using ActualCore =
      Detail::PipelineCore<true, PaletteBinding, typename Actual::NodeList>;
  using CanonicalCore =
      Detail::PipelineCore<true, PaletteBinding, typename Canonical::NodeList>;
  using Phase =
      Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                         Stage::Colorize<Policy>,
                         Stage::Modulate<Color::Modulate<PhaseLayerProvider>>>;
  using PhaseCore =
      Detail::PipelineCore<true, PaletteBinding, typename Phase::NodeList>;
  static_assert(
      ActualCore::template fused_color_tail_supported<1, FieldSample>());
  static_assert(
      !CanonicalCore::template fused_color_tail_supported<1, FieldSample>());
  static_assert(
      !PhaseCore::template fused_color_tail_supported<1, FieldSample>());

  std::array<int8_t, Pattern::SphereLutView::SIZE> ac, bc;
  for (size_t i = 0; i < ac.size(); ++i) {
    ac[i] = int(i * 37 % 255) - 127;
    bc[i] = int((i * 83 + 41) % 255) - 127;
  }
  Modulation::Prepared pa{}, pb{};
  pa.cube = {ac.data(), true};
  pb.cube = {bc.data(), true};
  pa.evaluate = pb.evaluate = Modulation::sample_cube;
  Color::ModulationState mod{&pb, &pa, Pattern::SourceKind::NOISE, {}, true};
  frame.modulation = &pb;
  frame.previous_modulation = &pa;
  for (auto mapping :
       {Color::PaletteMapping::LINEAR, Color::PaletteMapping::CUP,
        Color::PaletteMapping::BELL, Color::PaletteMapping::REVERSE})
    for (auto mechanic : {ColorModulation::Mechanic::HUE_ROTATION,
                          ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT,
                          ColorModulation::Mechanic::OPACITY})
      for (float amount : {0.0f, 0.3f, 1.0f, 1.5958333f})
        for (float fraction :
             {0.0f, 0.000001f, 0.01f, 0.25f, 0.5f, 0.75f, 0.99999f, 1.0f})
          for (int v = 0; v <= 64; ++v) {
            if (amount > 1.0f &&
                mechanic != ColorModulation::Mechanic::HUE_ROTATION)
              continue;
            frame.temporal.fraction = fraction;
            auto prepared = Policy::prepare(frame);
            prepared.current.mapping =
                Color::PaletteMappingWeights::single(mapping);
            prepared.current.mapping_frequency = 2.0f;
            prepared.current.mapping_offset = 0.13f;
            prepared.current.brightness_envelope =
                Color::BrightnessEnvelope::ASCENDING;
            prepared.current.brightness_bottom = 0.2f;
            prepared.current.brightness_top = 0.8f;
            prepared.previous = prepared.current;
            prepared.previous.palette = &a.view();
            mod.mechanic = {mechanic, amount};
            frame.mechanic = mod.mechanic;
            const FieldSample field{v / 64.0f, 0.45f,
                                    math::Vector(v / 64.0f, 0.3f, 1.0f), 0.7f};
            auto pair = Policy::apply_pair(field, frame, prepared);
            Color::apply_modulation_in_place(pair.first, mod, true);
            Color::apply_modulation_in_place(pair.second, mod, false);
            const auto expected =
                blend_color_samples(pair.first, pair.second, fraction);
            const auto actual = Policy::apply_fused_modulation(
                field, frame, prepared, mod, fraction);
            HS_EXPECT_EQ(actual.color.color, expected.color.color);
            HS_EXPECT_EQ(actual.color.alpha, expected.color.alpha);
            HS_EXPECT_EQ(actual.pending_hue, expected.pending_hue);
            HS_EXPECT_EQ(actual.pending_lift, expected.pending_lift);
            HS_EXPECT_EQ(actual.brightness_gain, expected.brightness_gain);
            HS_EXPECT_TRUE(actual.hue_rotation == expected.hue_rotation);
            const auto x = Color::finish(actual), y = Color::finish(expected);
            HS_EXPECT_EQ(x.color, y.color);
            HS_EXPECT_EQ(x.alpha, y.alpha);
            const auto DIRECT = Policy::shade_fused_modulation(
                field, frame, prepared, mod, fraction);
            HS_EXPECT_EQ(DIRECT.color, y.color);
            HS_EXPECT_EQ(DIRECT.alpha, y.alpha);
            const auto PIXEL =
                Actual::shade(field.sphere, Actual::prepare(frame));
            const auto REFERENCE =
                Canonical::shade(field.sphere, Canonical::prepare(frame));
            HS_EXPECT_EQ(PIXEL.color, REFERENCE.color);
            HS_EXPECT_EQ(PIXEL.alpha, REFERENCE.alpha);
          }
}

/** @brief Fused cached bases retain canonical hue, opacity and inactive fallbacks. */
inline void test_pullback_fused_cached_modulation() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  std::array<uint8_t, 3 * BakedPalette::required_arena_bytes()> buffer{};
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage a, b, translucent;
  a.bake(arena, CoordinatePalette{1.0f});
  b.bake(arena, ConstantPalette{Color4(Pixel(32000, 44000, 13000), 1.0f)});
  translucent.bake(arena,
                   ConstantPalette{Color4(Pixel(32000, 44000, 13000), 0.6f)});
  PaletteLabLut cache;
  cache.bake(a.view(), b.view());
  PaletteFrame frame;
  frame.palette = &b.view();
  frame.temporal.palette = &a.view();
  Color::ModulationState modulation{
      nullptr,
      nullptr,
      Pattern::SourceKind::DISPLACEMENT,
      {ColorModulation::Mechanic::HUE_ROTATION, 1.5958333f},
      true};
  int cache_hits = 0;
  for (bool active : {false, true})
    for (float fraction : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f})
      for (float value : {0.0f, 0.13f, 0.5f, 0.99f, 1.0f}) {
        frame.lab = active ? &cache : nullptr;
        frame.temporal.fraction = fraction;
        const auto PREPARED = Policy::prepare(frame);
        const FieldSample FIELD{value, 0.45f, math::Y_AXIS, 0.37f};
        auto pair = Policy::apply_pair(FIELD, frame, PREPARED);
        Color::apply_modulation_in_place(pair.first, modulation, true);
        Color::apply_modulation_in_place(pair.second, modulation, false);
        const auto CANONICAL =
            blend_color_samples(pair.first, pair.second, fraction);
        const auto ACTUAL = Policy::apply_fused_modulation(
            FIELD, frame, PREPARED, modulation, fraction);
        auto expected = CANONICAL;
        if (ACTUAL.has_cached_lab) {
          ++cache_hits;
          HS_EXPECT_TRUE(active && fraction > 0.0f && fraction < 1.0f);
          HS_EXPECT_TRUE(
              Policy::apply_shared_base(FIELD, frame, PREPARED, expected));
          expected.pending_hue = CANONICAL.pending_hue;
          expected.pending_lift = CANONICAL.pending_lift;
          expected.color.alpha = CANONICAL.color.alpha;
          expected.brightness_gain = CANONICAL.brightness_gain;
        }
        HS_EXPECT_EQ(ACTUAL.pending_hue, CANONICAL.pending_hue);
        HS_EXPECT_EQ(ACTUAL.pending_lift, CANONICAL.pending_lift);
        HS_EXPECT_EQ(ACTUAL.color.alpha, CANONICAL.color.alpha);
        HS_EXPECT_EQ(ACTUAL.brightness_gain, CANONICAL.brightness_gain);
        const auto FINISHED = Color::finish(ACTUAL);
        const auto EXPECTED = Color::finish(expected);
        HS_EXPECT_EQ(FINISHED.color, EXPECTED.color);
        HS_EXPECT_EQ(FINISHED.alpha, EXPECTED.alpha);
        const auto DIRECT = Policy::shade_fused_modulation(
            FIELD, frame, PREPARED, modulation, fraction);
        HS_EXPECT_EQ(DIRECT.color, FINISHED.color);
        HS_EXPECT_EQ(DIRECT.alpha, FINISHED.alpha);
      }
  HS_EXPECT_GT(cache_hits, 0);
  cache.bake(a.view(), translucent.view());
  frame.palette = &translucent.view();
  frame.lab = &cache;
  frame.temporal.fraction = 0.5f;
  const auto PREPARED = Policy::prepare(frame);
  const auto FALLBACK = Policy::apply_fused_modulation(
      FieldSample{0.5f, 0.45f, math::Y_AXIS, 0.37f}, frame, PREPARED,
      modulation, 0.5f);
  HS_EXPECT_FALSE(FALLBACK.has_cached_lab);
}

/** @brief Direct Lab shading reuses endpoint taps through canonical gamut fallback. */
inline void test_pullback_fused_lab_direct_fallback() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  std::array<uint8_t, 2 * BakedPalette::required_arena_bytes()> bytes{};
  Arena arena(bytes.data(), bytes.size());
  BakedPaletteStorage a, b;
  a.bake(arena, ConstantPalette{Color4(Pixel(1626, 1306, 43097), 1.0f)});
  b.bake(arena, ConstantPalette{Color4(Pixel(1451, 1516, 42130), 1.0f)});
  PaletteLabLut cache;
  cache.active = true;
  std::array<Pixel, HueRotationLutView::SIZE> hue;
  hue.fill(Pixel(111, 222, 333));
  PaletteFrame frame;
  frame.palette = &b.view();
  frame.lab = &cache;
  frame.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 1.0f};
  static int calls;
  const auto SOURCE = [](const Modulation::Prepared &,
                         const Pattern::SphereDirection &) {
    ++calls;
    return 0.9814453125f;
  };
  Modulation::Prepared previous{}, current{};
  previous.evaluate = current.evaluate = SOURCE;
  frame.previous_modulation = &previous;
  frame.modulation = &current;
  const auto MODULATION = PaletteLayerProvider::prepare(frame);
  for (bool guard : {false, true}) {
    cache.nodes.fill(guard ? OKLab{0.431782782f, 0.00217834837f, -0.24565056f}
                           : linear_rgb_to_oklab_fast(0.3f, 0.4f, 0.5f));
    for (bool same_base : {false, true}) {
      frame.temporal.palette = same_base ? &b.view() : &a.view();
      frame.hue_rotation = frame.temporal.hue_rotation =
          same_base ? hue.data() : nullptr;
      for (float fraction : {0.25f, 0.5f, 0.75f}) {
        frame.temporal.fraction = fraction;
        const auto PREPARED = Policy::prepare(frame);
        const FieldSample FIELD{0.5f, 0.61689651f, math::Y_AXIS, 0.37f};
        const auto CARRIER = Policy::apply_fused_modulation(
            FIELD, frame, PREPARED, MODULATION, fraction);
        HS_EXPECT_EQ(CARRIER.has_cached_lab, !guard);
        const auto EXPECTED = finish_color_sample(CARRIER);
        calls = 0;
        const auto ACTUAL = Policy::shade_fused_modulation(
            FIELD, frame, PREPARED, MODULATION, fraction);
        HS_EXPECT_EQ(ACTUAL.color, EXPECTED.color);
        HS_EXPECT_EQ(ACTUAL.alpha, EXPECTED.alpha);
        HS_EXPECT_EQ(calls, 2);
        if (guard && same_base)
          HS_EXPECT_EQ(ACTUAL.color, hue.front());
      }
    }
  }
}

/** @brief Cached current endpoints sample one live source and retain exact alpha. */
inline void test_pullback_fused_lab_current_endpoint() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  using Pipeline = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage, Stage::Colorize<Policy>,
      Stage::Modulate<Color::Modulate<PaletteLayerProvider>>>;
  std::array<uint8_t, 2 * BakedPalette::required_arena_bytes()> bytes{};
  Arena arena(bytes.data(), bytes.size());
  BakedPaletteStorage a, b;
  a.bake(arena, ConstantPalette{Color4(Pixel(10000, 40000, 20000), 1.0f)});
  b.bake(arena, ConstantPalette{Color4(Pixel(1451, 1516, 42130), 1.0f)});
  PaletteLabLut cache;
  PaletteFrame frame;
  frame.palette = &b.view();
  frame.temporal.palette = &a.view();
  frame.temporal.fraction = 1.0f;
  frame.lab = &cache;
  frame.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 1.0f};
  static int calls;
  const auto SOURCE = [](const Modulation::Prepared &prepared,
                         const Pattern::SphereDirection &) {
    ++calls;
    return prepared.state.phase.phase_turns;
  };
  Modulation::Prepared current{};
  current.evaluate = SOURCE;
  frame.modulation = &current;
  for (bool guard : {false, true}) {
    cache.bake(a.view(), b.view());
    if (guard)
      cache.nodes.fill({0.431782782f, 0.00217834837f, -0.24565056f});
    for (float scalar : {-0.125f, 0.0f, 0.13f, 0.9814453125f}) {
      current.state.phase.phase_turns = scalar;
      const auto MODULATION = PaletteLayerProvider::prepare(frame);
      const auto PREPARED = Policy::prepare(frame);
      for (float coverage : {0.0f, 0.5f, 1.0f}) {
        const FieldSample FIELD{0.5f, coverage, math::Y_AXIS, 0.0f};
        const auto REFERENCE = Policy::apply_fused_modulation(
            FIELD, frame, PREPARED, MODULATION, 1.0f);
        const auto EXPECTED = finish_color_sample(REFERENCE);
        auto cached = REFERENCE;
        cached.cached_lab =
            PaletteLabLut::sample(PREPARED.lab, cached.palette_coordinate);
        cached.has_cached_lab = true;
        Color::restore_cached_palette_base_near_gamut_gate(cached, [] {});
        if (guard && scalar == 0.9814453125f)
          HS_EXPECT_FALSE(cached.has_cached_lab);
        const auto CACHED = finish_color_sample(cached);
        calls = 0;
        const auto ACTUAL = Policy::shade_fused_modulation(
            FIELD, frame, PREPARED, MODULATION, 1.0f);
        HS_EXPECT_EQ(ACTUAL.color, CACHED.color);
        HS_EXPECT_EQ(ACTUAL.alpha, EXPECTED.alpha);
        HS_EXPECT_EQ(calls, 1);
        if (!guard) {
          HS_EXPECT_NEAR(ACTUAL.color.r, EXPECTED.color.r, 2);
          HS_EXPECT_NEAR(ACTUAL.color.g, EXPECTED.color.g, 2);
          HS_EXPECT_NEAR(ACTUAL.color.b, EXPECTED.color.b, 2);
        }
      }
      calls = 0;
      const auto PIXEL =
          Pipeline::shade(math::Y_AXIS, Pipeline::prepare(frame));
      HS_EXPECT_EQ(calls, 1);
      const auto DIRECT = Policy::shade_fused_modulation(
          FieldSample{0.5f, 0.5f, math::Y_AXIS, 0.0f}, frame, PREPARED,
          MODULATION, 1.0f);
      HS_EXPECT_EQ(PIXEL.color, DIRECT.color);
      HS_EXPECT_EQ(PIXEL.alpha, DIRECT.alpha);
    }
  }
  frame.temporal.fraction = 0.0f;
  frame.previous_modulation = &current;
  const auto PREPARED = Policy::prepare(frame);
  const auto MODULATION = PaletteLayerProvider::prepare(frame);
  const FieldSample FIELD{0.5f, 0.5f, math::Y_AXIS, 0.0f};
  const auto EXPECTED = finish_color_sample(
      Policy::apply_fused_modulation(FIELD, frame, PREPARED, MODULATION, 0.0f));
  calls = 0;
  const auto ACTUAL =
      Policy::shade_fused_modulation(FIELD, frame, PREPARED, MODULATION, 0.0f);
  HS_EXPECT_EQ(ACTUAL.color, EXPECTED.color);
  HS_EXPECT_EQ(ACTUAL.alpha, EXPECTED.alpha);
  HS_EXPECT_EQ(calls, 1);
}

/** @brief Gamut-gate crossings restore RGB-first color without normal cache fetches. */
inline void test_pullback_cached_gamut_gate_fallback() {
  using namespace Pullback;
  const float NAN_VALUE = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_TRUE(
      Color::near_palette_gamut_gate({NAN_VALUE, NAN_VALUE, NAN_VALUE}));
  ColorSample canonical{Color4(Pixel(1582, 1359, 42855), 0.61689651f), {}};
  canonical.pending_hue = 0.9814453125f;
  auto cached = canonical;
  cached.color.color = Pixel(0, 0, 0);
  cached.cached_lab = {0.431782782f, 0.00217834837f, -0.24565056f};
  cached.has_cached_lab = true;
  const auto REFERENCE = finish_color_sample(canonical);
  const auto UNGUARDED = finish_color_sample(cached);
  auto direct = cached;
  int direct_fetches = 0;
  const auto DIRECT = Color::finish_fused_color(direct, [&] {
    ++direct_fetches;
    direct.color.color = canonical.color.color;
    direct.hue_rotation = canonical.hue_rotation;
  });
  HS_EXPECT_EQ(direct_fetches, 1);
  HS_EXPECT_FALSE(direct.has_cached_lab);
  HS_EXPECT_EQ(DIRECT.color, REFERENCE.color);
  HS_EXPECT_EQ(DIRECT.alpha, REFERENCE.alpha);
  const Pixel EXPECTED = REFERENCE.color * REFERENCE.alpha;
  const Pixel WRONG = UNGUARDED.color * UNGUARDED.alpha;
  HS_EXPECT_GT(std::max({std::abs(int(EXPECTED.r) - WRONG.r),
                         std::abs(int(EXPECTED.g) - WRONG.g),
                         std::abs(int(EXPECTED.b) - WRONG.b)}),
               4000);
  int fetches = 0;
  HS_EXPECT_TRUE(
      Color::restore_cached_palette_base_near_gamut_gate(cached, [&] {
        ++fetches;
        cached.color.color = canonical.color.color;
        cached.hue_rotation = canonical.hue_rotation;
      }));
  HS_EXPECT_EQ(fetches, 1);
  HS_EXPECT_FALSE(cached.has_cached_lab);
  const auto GUARDED = finish_color_sample(cached);
  HS_EXPECT_EQ(GUARDED.color, REFERENCE.color);
  HS_EXPECT_EQ(GUARDED.alpha, REFERENCE.alpha);
  cached.cached_lab = linear_rgb_to_oklab_fast(0.3f, 0.4f, 0.5f);
  cached.has_cached_lab = true;
  cached.pending_hue = 0.0f;
  fetches = 0;
  HS_EXPECT_FALSE(Color::restore_cached_palette_base_near_gamut_gate(
      cached, [&] { ++fetches; }));
  HS_EXPECT_EQ(fetches, 0);
  HS_EXPECT_TRUE(cached.has_cached_lab);
  for (float hue : {0.0f, 0.13f, 0.5f, 0.9814453125f})
    for (float lift : {0.0f, 0.4f, 1.0f})
      for (float gain : {0.0f, 0.3f, 1.0f})
        for (auto lab : {linear_rgb_to_oklab_fast(0.3f, 0.4f, 0.5f),
                         OKLab{0.5f, 0.3f, -0.3f}}) {
          cached.cached_lab = lab;
          cached.pending_hue = hue;
          cached.pending_lift = lift;
          cached.brightness_gain = gain;
          cached.has_cached_lab = true;
          auto expected = cached;
          int reference_fetches = 0;
          Color::restore_cached_palette_base_near_gamut_gate(expected, [&] {
            ++reference_fetches;
            expected.color.color = canonical.color.color;
          });
          const auto FINISHED = finish_color_sample(expected);
          int actual_fetches = 0;
          const auto FAST = Color::finish_fused_color(cached, [&] {
            ++actual_fetches;
            cached.color.color = canonical.color.color;
          });
          HS_EXPECT_EQ(FAST.color, FINISHED.color);
          HS_EXPECT_EQ(FAST.alpha, FINISHED.alpha);
          HS_EXPECT_EQ(actual_fetches, reference_fetches);
          HS_EXPECT_EQ(cached.has_cached_lab, expected.has_cached_lab);
        }
}

/** @brief Prepared cube-pair eligibility preserves evaluator and active-source boundaries. */
inline void test_pullback_prepared_cube_pair() {
  using namespace Pullback;
  using Layer = Color::Modulate<PaletteLayerProvider>;
  using Policy = Color::GeneratedPalette<TemporalPaletteProvider>;
  constexpr size_t ALIGNMENT = alignof(Color::ModulationState);
  constexpr size_t PREVIOUS_SIZE = (offsetof(Color::ModulationState, enabled) +
                                    sizeof(bool) + ALIGNMENT - 1) /
                                   ALIGNMENT * ALIGNMENT;
  static_assert(sizeof(Color::ModulationState) == PREVIOUS_SIZE);
  std::array<int8_t, Pattern::SphereLutView::SIZE> ac, bc;
  std::array<float, Pattern::SphereLutView::SIZE> linear;
  for (size_t i = 0; i < ac.size(); ++i) {
    ac[i] = int(i * 37 % 255) - 127;
    bc[i] = int((i * 83 + 41) % 255) - 127;
    linear[i] = bc[i];
  }
  Modulation::Prepared a{}, b{};
  a.cube = {ac.data(), true};
  b.cube = {bc.data(), true};
  a.evaluate = b.evaluate = Modulation::sample_cube;
  std::array<uint8_t, 2 * BakedPalette::required_arena_bytes()> bytes{};
  Arena arena(bytes.data(), bytes.size());
  BakedPaletteStorage from, to;
  from.bake(arena, CoordinatePalette{0.2f});
  to.bake(arena, ConstantPalette{Color4(Pixel(32000, 44000, 13000), 0.8f)});
  PaletteFrame frame;
  frame.palette = &to.view();
  frame.temporal.palette = &from.view();
  frame.previous_modulation = &a;
  frame.modulation = &b;
  frame.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 1.5958333f};
  for (bool floats : {false, true}) {
    b.cube.linear_data = floats ? linear.data() : nullptr;
    const auto PREPARED = Layer::prepare(frame);
    HS_EXPECT_TRUE(PREPARED.paired_cube_ready);
    auto dynamic = PREPARED;
    dynamic.paired_cube_ready = false;
    for (float fraction : {0.00001f, 0.25f, 0.75f, 0.99999f}) {
      frame.temporal.fraction = fraction;
      const auto PALETTE = Policy::prepare(frame);
      for (int x = -3; x <= 3; ++x) {
        const FieldSample FIELD{0.13f, 0.47f, math::Vector(x, 0.7f, -1.0f),
                                0.3f};
        const auto EXPECTED = Policy::apply_fused_modulation(
            FIELD, frame, PALETTE, dynamic, fraction);
        const auto ACTUAL = Policy::apply_fused_modulation(
            FIELD, frame, PALETTE, PREPARED, fraction);
        HS_EXPECT_EQ(ACTUAL.pending_hue, EXPECTED.pending_hue);
        HS_EXPECT_EQ(ACTUAL.color.color, EXPECTED.color.color);
        HS_EXPECT_EQ(ACTUAL.color.alpha, EXPECTED.color.alpha);
        const auto FINISHED = Policy::shade_fused_modulation(
            FIELD, frame, PALETTE, PREPARED, fraction);
        const auto REFERENCE = finish_color_sample(EXPECTED);
        HS_EXPECT_EQ(FINISHED.color, REFERENCE.color);
        HS_EXPECT_EQ(FINISHED.alpha, REFERENCE.alpha);
      }
      const auto PIPELINE_FRAME = TemporalPalettePipeline::prepare(frame);
      HS_EXPECT_TRUE(std::get<2>(PIPELINE_FRAME.prepared).paired_cube_ready);
    }
  }
  static int calls;
  const auto CUSTOM = [](const Modulation::Prepared &prepared,
                         const Pattern::SphereDirection &) {
    ++calls;
    return prepared.state.phase.phase_turns;
  };
  a.state.phase.phase_turns = -0.3f;
  b.state.phase.phase_turns = 0.7f;
  for (int custom = 1; custom <= 3; ++custom) {
    a.evaluate = custom & 1 ? CUSTOM : Modulation::sample_cube;
    b.evaluate = custom & 2 ? CUSTOM : Modulation::sample_cube;
    const auto PREPARED = Layer::prepare(frame);
    HS_EXPECT_FALSE(PREPARED.paired_cube_ready);
    calls = 0;
    const auto PALETTE = Policy::prepare(frame);
    (void)Policy::shade_fused_modulation(
        FieldSample{0.5f, 0.5f, math::Y_AXIS, 0.0f}, frame, PALETTE, PREPARED,
        frame.temporal.fraction);
    HS_EXPECT_EQ(calls, custom == 3 ? 2 : 1);
  }
  a.evaluate = b.evaluate = Modulation::sample_cube;
  for (int inactive = 1; inactive <= 3; ++inactive) {
    a.cube.active = !(inactive & 1);
    b.cube.active = !(inactive & 2);
    HS_EXPECT_FALSE(Layer::prepare(frame).paired_cube_ready);
  }
  a.cube.active = b.cube.active = true;
  frame.modulation = nullptr;
  HS_EXPECT_FALSE(Layer::prepare(frame).paired_cube_ready);
  frame.modulation = &b;
  frame.previous_modulation = nullptr;
  HS_EXPECT_FALSE(Layer::prepare(frame).paired_cube_ready);
  frame.modulation = nullptr;
  HS_EXPECT_FALSE(Layer::prepare(frame).paired_cube_ready);
  const Color::ModulationState MANUAL{&b, &a, Pattern::SourceKind::NOISE,
                                      frame.mechanic, true};
  HS_EXPECT_FALSE(MANUAL.paired_cube_ready);
}

/** @brief Paired cube samples preserve individual samplers and custom dispatch. */
inline void test_pullback_modulation_paired_samples() {
  using namespace Pullback;
  std::array<int8_t, Pattern::SphereLutView::SIZE> previous_cube;
  std::array<int8_t, Pattern::SphereLutView::SIZE> current_cube;
  std::array<float, Pattern::SphereLutView::SIZE> linear_cube;
  for (size_t index = 0; index < previous_cube.size(); ++index) {
    previous_cube[index] = static_cast<int8_t>((index * 37) % 255 - 127);
    current_cube[index] = static_cast<int8_t>((index * 83 + 19) % 255 - 127);
    linear_cube[index] = current_cube[index];
  }
  Modulation::Prepared previous{};
  Modulation::Prepared current{};
  previous.cube = {previous_cube.data(), true};
  current.cube = {current_cube.data(), true};
  previous.evaluate = current.evaluate = Modulation::sample_cube;
  for (bool floats : {false, true}) {
    current.cube.linear_data = floats ? linear_cube.data() : nullptr;
    for (int x = -3; x <= 3; ++x)
      for (int y = -3; y <= 3; ++y)
        for (int z = -3; z <= 3; ++z) {
          if (x == 0 && y == 0 && z == 0)
            continue;
          const auto DIRECTION =
              Pattern::SphereDirection::admit(math::Vector(x, y, z));
          const auto PAIR =
              Modulation::sample_pair(previous, current, DIRECTION);
          HS_EXPECT_EQ(PAIR.first, Modulation::sample(previous, DIRECTION));
          HS_EXPECT_EQ(PAIR.second, Modulation::sample(current, DIRECTION));
        }
  }
  static int calls;
  calls = 0;
  const auto CUSTOM = [](const Modulation::Prepared &prepared,
                         const Pattern::SphereDirection &) {
    ++calls;
    return prepared.state.phase.phase_turns;
  };
  previous.evaluate = current.evaluate = CUSTOM;
  previous.state.phase.phase_turns = -2.0f;
  current.state.phase.phase_turns = 3.0f;
  const auto PAIR = Modulation::sample_pair(
      previous, current, Pattern::SphereDirection::from_unit(math::Y_AXIS));
  HS_EXPECT_EQ(calls, 2);
  HS_EXPECT_EQ(PAIR.first, -1.0f);
  HS_EXPECT_EQ(PAIR.second, 1.0f);
}

/** @brief Paired modulation keeps canonical endpoint colors and mechanic order. */
inline void test_pullback_modulation_paired_colors() {
  using namespace Pullback;
  std::array<int8_t, Pattern::SphereLutView::SIZE> previous_cube;
  std::array<int8_t, Pattern::SphereLutView::SIZE> current_cube;
  for (size_t index = 0; index < previous_cube.size(); ++index) {
    previous_cube[index] = static_cast<int8_t>((index * 17) % 255 - 127);
    current_cube[index] = static_cast<int8_t>((index * 41 + 37) % 255 - 127);
  }
  Modulation::Prepared previous{};
  Modulation::Prepared current{};
  previous.cube = {previous_cube.data(), true};
  current.cube = {current_cube.data(), true};
  previous.evaluate = current.evaluate = Modulation::sample_cube;
  Color::ModulationState state{
      &current, &previous, Pattern::SourceKind::NOISE, {}, true};
  for (auto mechanic : {ColorModulation::Mechanic::HUE_ROTATION,
                        ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT,
                        ColorModulation::Mechanic::OPACITY})
    for (float amount : {0.0f, 0.3f, 1.5f})
      for (bool same_direction : {false, true})
        for (int x = -2; x <= 2; ++x) {
          if (amount > 1.0f &&
              mechanic != ColorModulation::Mechanic::HUE_ROTATION)
            continue;
          state.mechanic = {mechanic, amount};
          ColorSample a{Color4(Pixel(41000, 17000, 33000), 0.25f),
                        {0.3f, 0.4f, math::Vector(x, 1.0f, 0.3f), 0.7f}};
          ColorSample b{
              Color4(Pixel(27000, 53000, 9000), 0.75f),
              {0.6f, 0.8f,
               same_direction ? a.field.sphere : math::Vector(0.3f, x, 1.0f),
               0.2f}};
          a.pending_hue = 0.93f;
          b.pending_hue = 0.07f;
          a.pending_lift = 0.2f;
          b.pending_lift = 0.4f;
          const auto EXPECTED_A = Color::apply_modulation(a, state, true);
          const auto EXPECTED_B = Color::apply_modulation(b, state, false);
          Color::apply_modulation_pair_in_place(a, b, state);
          for (auto pair :
               {std::pair{a, EXPECTED_A}, std::pair{b, EXPECTED_B}}) {
            HS_EXPECT_EQ(pair.first.color.color, pair.second.color.color);
            HS_EXPECT_EQ(pair.first.color.alpha, pair.second.color.alpha);
            HS_EXPECT_EQ(pair.first.pending_hue, pair.second.pending_hue);
            HS_EXPECT_EQ(pair.first.pending_lift, pair.second.pending_lift);
            HS_EXPECT_EQ(pair.first.field.sphere, pair.second.field.sphere);
            const auto FINISHED = Color::finish(pair.first);
            const auto EXPECTED = Color::finish(pair.second);
            HS_EXPECT_EQ(FINISHED.color, EXPECTED.color);
            HS_EXPECT_EQ(FINISHED.alpha, EXPECTED.alpha);
          }
        }
}

/** @brief Modulation identities retain common and selected inputs only. */
inline void test_pullback_modulation_field_identity() {
  using namespace Pullback::Modulation;
  using Pattern::SourceKind;
  Config config;
  config.version = 7;
  config.evaluation = Pattern::Evaluation::DIRECT;
  config.noise.scale = 2.0f;
  config.sine.frequency = 3.0f;
  config.interference.frequency_b = 4.0f;
  config.grid.offset_b = 0.5f;
  config.seed = 123;
  const State STATE{SourceKind::NOISE, 1, {0.25f}, true};
  for (SourceKind source :
       {SourceKind::NOISE, SourceKind::DIRECTIONAL_SINE,
        SourceKind::SPHERE_INTERFERENCE, SourceKind::SPHERE_SMOOTH_GRID,
        SourceKind::DISPLACEMENT}) {
    config.source = source;
    const FieldIdentity IDENTITY = field_identity(config, STATE);
    HS_EXPECT_EQ(IDENTITY.source, source);
    HS_EXPECT_EQ(IDENTITY.version, 7);
    HS_EXPECT_EQ(IDENTITY.evaluation, Pattern::Evaluation::DIRECT);
    HS_EXPECT_EQ(IDENTITY.phase, Pattern::PhaseState{0.25f});
    HS_EXPECT_EQ(IDENTITY.bake_layout_version,
                 Pattern::SphereLutView::LAYOUT_VERSION);
    HS_EXPECT_EQ(IDENTITY.seed, source == SourceKind::NOISE ? 123 : 0);
    HS_EXPECT_EQ(IDENTITY.noise, source == SourceKind::NOISE
                                     ? config.noise
                                     : Pattern::NoiseParams{});
    HS_EXPECT_EQ(IDENTITY.sine, source == SourceKind::DIRECTIONAL_SINE
                                    ? config.sine
                                    : Pattern::SineParams{});
    HS_EXPECT_EQ(IDENTITY.interference,
                 source == SourceKind::SPHERE_INTERFERENCE
                     ? config.interference
                     : Pattern::InterferenceParams{});
    HS_EXPECT_EQ(IDENTITY.grid, source == SourceKind::SPHERE_SMOOTH_GRID
                                    ? config.grid
                                    : Pattern::GridParams{});
  }
  config.source = static_cast<SourceKind>(255);
  const FieldIdentity UNKNOWN = field_identity(config, STATE);
  HS_EXPECT_EQ(UNKNOWN.source, static_cast<SourceKind>(255));
  HS_EXPECT_EQ(UNKNOWN.version, 7);
  HS_EXPECT_EQ(UNKNOWN.evaluation, Pattern::Evaluation::DIRECT);
  HS_EXPECT_EQ(UNKNOWN.phase, STATE.phase);
  HS_EXPECT_EQ(UNKNOWN.noise, Pattern::NoiseParams{});
  HS_EXPECT_EQ(UNKNOWN.sine, Pattern::SineParams{});
  HS_EXPECT_EQ(UNKNOWN.interference, Pattern::InterferenceParams{});
  HS_EXPECT_EQ(UNKNOWN.grid, Pattern::GridParams{});
  HS_EXPECT_EQ(UNKNOWN.seed, 0);
  HS_EXPECT_EQ(UNKNOWN.bake_layout_version,
               Pattern::SphereLutView::LAYOUT_VERSION);
}

/** @brief Modulation equality distinguishes field inputs from dormant controls. */
inline void test_pullback_modulation_field_equality() {
  using namespace Pullback::Modulation;
  using Pattern::SourceKind;
  for (SourceKind source :
       {SourceKind::NOISE, SourceKind::DIRECTIONAL_SINE,
        SourceKind::SPHERE_INTERFERENCE, SourceKind::SPHERE_SMOOTH_GRID,
        SourceKind::DISPLACEMENT}) {
    Prepared baseline{};
    baseline.config.source = source;
    baseline.config.evaluation = Pattern::Evaluation::DIRECT;
    const auto expect = [&](auto change, bool equal) {
      Prepared other = baseline;
      change(other);
      HS_EXPECT_EQ(same_field(baseline, other), equal);
      HS_EXPECT_EQ(same_field(other, baseline), equal);
      HS_EXPECT_EQ(field_identity(baseline.config, baseline.state) ==
                       field_identity(other.config, other.state),
                   equal);
    };
    HS_EXPECT_TRUE(same_field(baseline, baseline));
    expect(
        [&](auto &p) {
          p.config.source = source == SourceKind::NOISE
                                ? SourceKind::DISPLACEMENT
                                : SourceKind::NOISE;
        },
        false);
    expect([](auto &p) { ++p.config.version; }, false);
    expect(
        [](auto &p) { p.config.evaluation = Pattern::Evaluation::BAKED_CUBE; },
        false);
    expect([](auto &p) { p.state.phase.phase_turns = 0.25f; }, false);
    expect([](auto &p) { ++p.config.seed; }, source != SourceKind::NOISE);
    expect([](auto &p) { p.config.noise.scale = 2.0f; },
           source != SourceKind::NOISE);
    expect([](auto &p) { p.config.noise.speed = 0.5f; },
           source != SourceKind::NOISE);
    expect([](auto &p) { p.config.sine.frequency = 2.0f; },
           source != SourceKind::DIRECTIONAL_SINE);
    expect([](auto &p) { p.config.sine.axis.x = 0.25f; },
           source != SourceKind::DIRECTIONAL_SINE);
    expect([](auto &p) { p.config.interference.frequency = 2.0f; },
           source != SourceKind::SPHERE_INTERFERENCE);
    expect([](auto &p) { p.config.interference.axis_b.z = 0.25f; },
           source != SourceKind::SPHERE_INTERFERENCE);
    expect([](auto &p) { p.config.grid.frequency = 2.0f; },
           source != SourceKind::SPHERE_SMOOTH_GRID);
    expect([](auto &p) { p.config.grid.offset_b = 0.25f; },
           source != SourceKind::SPHERE_SMOOTH_GRID);
    expect([](auto &p) { p.config.enabled = false; }, true);
    expect(
        [](auto &p) {
          p.config.mechanic.kind = ColorModulation::Mechanic::OPACITY;
        },
        true);
    expect([](auto &p) { p.config.mechanic.amount = 0.5f; }, true);
    expect([](auto &p) { p.state.source = static_cast<SourceKind>(255); },
           true);
    expect([](auto &p) { ++p.state.version; }, true);
    expect([](auto &p) { p.state.initialized = true; }, true);
  }
  Prepared unknown{};
  unknown.config.source = static_cast<SourceKind>(255);
  HS_EXPECT_FALSE(same_field(unknown, unknown));
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
  HS_RUN(test_pullback_modulation_field_identity);
  HS_RUN(test_pullback_modulation_field_equality);
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
  HS_RUN(test_pullback_lattice_circle_identity);
  HS_RUN(test_pullback_custom_grid_preparation);
  HS_RUN(test_pullback_hexagonal_edges);
  HS_RUN(test_pullback_displacement_overflow);
  HS_RUN(test_pullback_warp_phase_loop);
  HS_RUN(test_pullback_dodecahedral_matrix_fold);
  HS_RUN(test_pullback_euler_curl_specialization);
  HS_RUN(test_pullback_periodic_ripple);
  HS_RUN(test_pullback_ripple_supplied_center);
  HS_RUN(test_pullback_ripple_burst_centers);
  HS_RUN(test_pullback_radial_ripple_endpoint_reuse);
  HS_RUN(test_pullback_radial_cell_validity);
  HS_RUN(test_pullback_radial_validity_lifetime);
  HS_RUN(test_pullback_radial_ripple_crest);
  HS_RUN(test_pullback_radial_ripple_fallback);
  HS_RUN(test_pullback_radial_ripple_lazy_fallback);
  HS_RUN(test_pullback_ripple_lazy_preparation);
  HS_RUN(test_pullback_radial_ripple_geometry);
  HS_RUN(test_pullback_group_rotation_polynomial);
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
  HS_RUN(test_pullback_phase_tail_dispatch);
  HS_RUN(test_pullback_palette_lab_nodes);
  HS_RUN(test_pullback_cached_lab_finish);
  HS_RUN(test_pullback_palette_lab_qualification);
  HS_RUN(test_pullback_shared_lab_dispatch);
  HS_RUN(test_pullback_paired_color_dispatch);
  HS_RUN(test_pullback_paired_palette_controls);
  HS_RUN(test_pullback_modulation_temporal_alpha);
  HS_RUN(test_pullback_stacked_modulations);
  HS_RUN(test_pullback_composed_lightness_hue);
  HS_RUN(test_pullback_finish_brightness_quantization);
  HS_RUN(test_pullback_temporal_composed_modulations);
  HS_RUN(test_pullback_temporal_deferred_hue_wrap);
  HS_RUN(test_pullback_temporal_deferred_in_place);
  HS_RUN(test_pullback_fused_modulation_endpoints);
  HS_RUN(test_pullback_fused_cached_modulation);
  HS_RUN(test_pullback_fused_lab_direct_fallback);
  HS_RUN(test_pullback_fused_lab_current_endpoint);
  HS_RUN(test_pullback_cached_gamut_gate_fallback);
  HS_RUN(test_pullback_prepared_cube_pair);
  HS_RUN(test_pullback_modulation_paired_samples);
  HS_RUN(test_pullback_modulation_paired_colors);
  return fixture.result();
}

} // namespace pullback_tests
} // namespace hs_test
