/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Pullback carrier, stage and policy contracts.
 */
#pragma once

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
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

  static_assert(Pullback::FAMILY_OF<Pullback::SphereSample> == 0);
  static_assert(Pullback::FAMILY_OF<Pullback::PlaneSample> == 1);
  static_assert(Pullback::FAMILY_OF<Pullback::FieldSample> == 2);
  static_assert(Pullback::FAMILY_OF<Pullback::ColorSample> == 3);
  static_assert(Pullback::FAMILY_OF<NotACarrier> ==
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

/** @brief Schema probe placing one curve in numeric or extra rows.
 * @tparam Curve Curve under test.
 * @tparam Extra Whether the extra row carries the curve. */
template <Pullback::FieldCurve Curve, bool Extra> struct CurveSubsetProbe {
  float value = 0.5f; ///< In-range numeric value.
  /// Numeric field table.
  static constexpr auto FIELDS = std::array{Pullback::Field<CurveSubsetProbe>{
      "value", &CurveSubsetProbe::value, nullptr, 0.1f, 1.0f,
      Extra ? Pullback::FieldCurve::LERP : Curve}};
  /// Extra schema row.
  static constexpr auto EXTRA_SCHEMA = [] {
    Pullback::ParamFieldInfo field{};
    field.id = "extra";
    field.min = 0.1f;
    field.max = 1.0f;
    field.def = 0.5f;
    field.curve = Extra ? Curve : Pullback::FieldCurve::LERP;
    return std::array{field};
  }();
};

/** @brief Registered pullback fields reject midpoint and extrapolating curves. */
inline void test_pullback_registered_curve_subset() {
  using Pullback::FieldCurve;
  static_assert(std::is_same_v<FieldCurve, Control::FieldCurve>);
  HS_EXPECT_TRUE((Pullback::field_curves_supported<
                  CurveSubsetProbe<FieldCurve::LERP, false>>()));
  HS_EXPECT_TRUE((Pullback::field_curves_supported<
                  CurveSubsetProbe<FieldCurve::LOG_POSITIVE, false>>()));
  HS_EXPECT_TRUE((Pullback::field_curves_supported<
                  CurveSubsetProbe<FieldCurve::SHORTEST_PERIODIC, false>>()));
  HS_EXPECT_TRUE((Pullback::field_curves_supported<
                  CurveSubsetProbe<FieldCurve::SHORTEST_TURN, false>>()));
  HS_EXPECT_TRUE((Pullback::field_curves_supported<
                  CurveSubsetProbe<FieldCurve::SNAP, false>>()));
  HS_EXPECT_FALSE((Pullback::field_curves_supported<
                   CurveSubsetProbe<FieldCurve::MIDPOINT, false>>()));
  HS_EXPECT_FALSE((Pullback::field_curves_supported<
                   CurveSubsetProbe<FieldCurve::RAW_LINEAR, false>>()));
  HS_EXPECT_FALSE((Pullback::field_curves_supported<
                   CurveSubsetProbe<FieldCurve::MIDPOINT, true>>()));
  HS_EXPECT_FALSE((Pullback::field_curves_supported<
                   CurveSubsetProbe<FieldCurve::RAW_LINEAR, true>>()));
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
  static_assert(TestPipeline::ANY_STAGE<IsCrossingStage>);
  static_assert(!TestPipeline::ANY_STAGE<NoStageMatches>);
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

/** @brief Affine providers require a repeat lattice alongside their scalars. */
inline void test_pullback_affine_provider_repeat_contract() {
  struct AffineState {
    using Binding = TestBinding;
    using FrameState = TestFrame;
    static Pullback::Warp::PreparedAffineSlot prepare(const FrameState &) {
      return Pullback::Warp::prepare(Pullback::Warp::AffineParams{}, 0.0f, 0.0f,
                                     {});
    }
    static bool path_length_required(const FrameState &) { return true; }
  };
  struct MissingRepeatState : AffineState {
    struct Scalars {
      float inverse_scale_x, inverse_scale_y, shear, translation_x,
          translation_y;
    };
    struct Prepared {
      float rotation_cos, rotation_sin;
      Scalars affine;
    };
    static Prepared prepare(const FrameState &) { return {}; }
  };
  struct WrongRepeatState : AffineState {
    struct Scalars : MissingRepeatState::Scalars {
      int repeat;
    };
    struct Prepared {
      float rotation_cos, rotation_sin;
      Scalars affine;
    };
    static Prepared prepare(const FrameState &) { return {}; }
  };
  using Policy = Pullback::Warp::AffineFrame<AffineState>;
  static_assert(Policy::PROVIDER_VALID<TestBinding>);
  static_assert(!Pullback::Warp::AffineFrame<
                MissingRepeatState>::PROVIDER_VALID<TestBinding>);
  static_assert(!Pullback::Warp::AffineFrame<WrongRepeatState>::PROVIDER_VALID<
                TestBinding>);
  const TestFrame frame;
  const math::Complex point(0.25f, -0.75f);
  const auto actual = Policy::apply(point, {}, frame, Policy::prepare(frame));
  HS_EXPECT_EQ(actual.coords.re, point.re);
  HS_EXPECT_EQ(actual.coords.im, point.im);
  HS_EXPECT_EQ(actual.path_length, 0.0f);
}

/** @brief Square cells evaluate the circle-to-box blend formula. */
inline void test_pullback_square_cells_formula() {
  const auto reference = [](const Pattern::PlanePosition &coordinate,
                            const Pattern::TessellationCellParams &params) {
    const float X =
        math::wrap_t(params.cell_scale * coordinate.value.re + 0.5f) - 0.5f;
    const float Y =
        math::wrap_t(params.cell_scale * coordinate.value.im + 0.5f) - 0.5f;
    const float CIRCLE = sqrtf(X * X + Y * Y) - params.cell_radius;
    const float BX = fabsf(X) - params.cell_radius;
    const float BY = fabsf(Y) - params.cell_radius;
    const float SQUARE = sqrtf(fmaxf(BX, 0.0f) * fmaxf(BX, 0.0f) +
                               fmaxf(BY, 0.0f) * fmaxf(BY, 0.0f)) +
                         fminf(fmaxf(BX, BY), 0.0f);
    return 1.0f - 2.0f * math::smooth_ramp(
                             -params.cell_softness, params.cell_softness,
                             hs::lerp(CIRCLE, SQUARE, params.cell_shape));
  };
  hs::Pcg32 random(0x12345678u);
  for (float scale : {1.0f / 64.0f, 0.98971874f, 1.0f, 8.0f})
    for (float radius : {1.0f / 64.0f, 0.25f, 0.37757313f, 0.49f})
      for (float softness : {1.0f / 1024.0f, 0.1f, 0.34563965f, 1.0f})
        for (float blend : {0.0f, 0.00001f, 0.25f, 0.99999f, 1.0f}) {
          Pattern::TessellationCellParams params;
          params.cell_scale = scale;
          params.cell_radius = radius;
          params.cell_softness = softness;
          params.cell_shape = blend;
          for (int sample = 0; sample < 256; ++sample) {
            const float X = rand_uniform(random, -1024.0f, 1024.0f);
            const float Y = rand_uniform(random, -1024.0f, 1024.0f);
            const Pattern::PlanePosition COORDINATE{
                math::Complex(sample < 8 ? float(sample - 4) * 0.5f : X,
                              sample < 8 ? 0.0f : Y)};
            const float EXPECTED = reference(COORDINATE, params);
            const float ACTUAL = Pattern::square_cells(COORDINATE, params);
#if HS_TEST_FAST_MATH
            HS_EXPECT_NEAR(ACTUAL, EXPECTED,
                           8.0f * std::numeric_limits<float>::epsilon());
#else
            HS_EXPECT_TRUE(std::memcmp(&EXPECTED, &ACTUAL, sizeof(float)) == 0);
#endif
          }
        }
}

/** @brief Independent double-precision tessellation reference. */
struct TessellationReference {
  static constexpr double HALF_TURN = 3.14159265358979323846;
  Pattern::TessellationKind kind;
  Pattern::TessellationFill fill;
  Pattern::TessellationSourceParams params;

  /** @brief Exact signed distance to a regular polygon of inradius @p radius
      whose first outward edge normal points along @p angle. */
  static double polygon_distance(double x, double y, int sides, double angle,
                                 double radius) {
    const double step = 2.0 * HALF_TURN / sides;
    double inside = -1e30;
    for (int edge = 0; edge < sides; ++edge) {
      const double theta = angle + edge * step;
      inside =
          std::max(inside, x * std::cos(theta) + y * std::sin(theta) - radius);
    }
    if (inside <= 0.0)
      return inside;
    const double circumradius = radius / std::cos(HALF_TURN / sides);
    double outside = 1e30;
    for (int edge = 0; edge < sides; ++edge) {
      const double theta = angle + edge * step;
      const double ax = circumradius * std::cos(theta - 0.5 * step);
      const double ay = circumradius * std::sin(theta - 0.5 * step);
      const double bx = circumradius * std::cos(theta + 0.5 * step);
      const double by = circumradius * std::sin(theta + 0.5 * step);
      const double t =
          std::clamp(((x - ax) * (bx - ax) + (y - ay) * (by - ay)) /
                         ((bx - ax) * (bx - ax) + (by - ay) * (by - ay)),
                     0.0, 1.0);
      outside = std::min(
          outside, std::hypot(x - ax - t * (bx - ax), y - ay - t * (by - ay)));
    }
    return outside;
  }

  static double ramp(double edge0, double edge1, double value) {
    const double t = std::clamp((value - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
  }

  double operator()(const math::Complex &point) const {
    using K = Pattern::TessellationKind;
    const double sqrt3 = std::sqrt(3.0);
    const double u = double(params.cell_scale) * point.re;
    const double v = double(params.cell_scale) * point.im;
    double x = 0.0, y = 0.0, angle = 0.0;
    int sides = 4;
    if (kind == K::SQUARE) {
      x = u - std::round(u);
      y = v - std::round(v);
    } else if (kind == K::HEXAGONAL) {
      sides = 6;
      const double row = std::round(v / (0.5 * sqrt3));
      const double column = std::round(u - 0.5 * row);
      double best = 1e30;
      for (double j = row - 2; j <= row + 2; ++j)
        for (double i = column - 2; i <= column + 2; ++i) {
          const double dx = u - (i + 0.5 * j);
          const double dy = v - 0.5 * sqrt3 * j;
          if (dx * dx + dy * dy < best) {
            best = dx * dx + dy * dy;
            x = dx;
            y = dy;
          }
        }
    } else {
      sides = 3;
      const double b = v / (0.5 * sqrt3);
      const double a = u - 0.5 * b;
      double best = -1e30;
      for (double j = std::floor(b) - 1; j <= std::floor(b) + 1; ++j)
        for (double i = std::floor(a) - 1; i <= std::floor(a) + 1; ++i)
          for (const bool up : {true, false}) {
            const double third = up ? 1.0 / 3.0 : 2.0 / 3.0;
            const double cx = (i + third) + 0.5 * (j + third);
            const double cy = 0.5 * sqrt3 * (j + third);
            const double dx = sqrt3 * (u - cx);
            const double dy = sqrt3 * (v - cy);
            const double facing = up ? -0.5 * HALF_TURN : 0.5 * HALF_TURN;
            const double depth = -polygon_distance(dx, dy, 3, facing, 0.5);
            if (depth > best) {
              best = depth;
              x = dx;
              y = up ? dy : -dy;
            }
          }
      angle = -0.5 * HALF_TURN;
    }
    if (fill == Pattern::TessellationFill::EDGES) {
      const double distance = -polygon_distance(x, y, sides, angle, 0.5);
      return 1.0 -
             2.0 * ramp(params.line_thickness,
                        double(params.line_thickness) + params.line_softness,
                        distance);
    }
    const double circle = std::hypot(x, y) - params.cell_radius;
    const double polygon =
        polygon_distance(x, y, sides, angle, params.cell_radius);
    const double distance = circle + (polygon - circle) * params.cell_shape;
    return 1.0 - 2.0 * ramp(-double(params.cell_softness),
                            double(params.cell_softness), distance);
  }
};

/** @brief Every tessellation kind and fill matches an independent polygon
    reference. */
inline void test_pullback_tessellation_reference() {
  using K = Pattern::TessellationKind;
  using F = Pattern::TessellationFill;
  hs::Pcg32 random(0x9e3779b9u);
  const auto next = [&] { return rand_uniform(random, -3.0f, 3.0f); };
  for (const K kind : {K::TRIANGULAR, K::SQUARE, K::HEXAGONAL})
    for (const F fill : {F::CELLS, F::EDGES})
      for (const float shape : {0.0f, 0.35f, 1.0f}) {
        HS_CONTEXT("kind / fill", static_cast<int>(kind),
                   static_cast<int>(fill));
        HS_CONTEXT("shape percent", static_cast<int>(shape * 100.0f));
        TessellationReference reference{kind, fill, {}};
        reference.params.cell_scale = 1.7f;
        reference.params.cell_shape = shape;
        reference.params.cell_radius = 0.3f;
        reference.params.cell_softness = 0.08f;
        reference.params.line_thickness = 0.05f;
        reference.params.line_softness = 0.04f;
        int inside = 0;
        int outside = 0;
        for (int sample = 0; sample < 2000; ++sample) {
          const float X = next();
          const float Y = next();
          const math::Complex point(X, Y);
          const float actual = Pattern::tessellation(
              Pattern::PlanePosition{point}, reference.params, kind, fill);
          const double expected = reference(point);
          HS_EXPECT_NEAR(actual, expected, 2e-4);
          inside += expected > 0.99;
          outside += expected < -0.99;
        }
        HS_EXPECT_GT(inside, 100);
        HS_EXPECT_GT(outside, 100);
      }
}

/** @brief Hexagonal cells pack circles and hexagons on the hexagonal lattice,
    with each hexagon's edges facing its six neighbours. */
inline void test_pullback_hexagonal_cells_pack_hex_lattice() {
  using K = Pattern::TessellationKind;
  using F = Pattern::TessellationFill;
  Pattern::TessellationCellParams params;
  params.cell_scale = 2.5f;
  params.cell_radius = 0.4f;
  params.cell_softness = 0.05f;
  const float period = 0.4f;
  const float half_sqrt3 = 0.8660254038f;
  const auto at = [&](float u, float v) {
    return Pattern::tessellation<K::HEXAGONAL, F::CELLS>(
        Pattern::PlanePosition{math::Complex(period * u, period * v)}, params);
  };
  for (const float shape : {0.0f, 1.0f}) {
    HS_CONTEXT("shape", static_cast<int>(shape));
    params.cell_shape = shape;
    for (int j = -3; j <= 3; ++j)
      for (int i = -3; i <= 3; ++i) {
        const float u = i + 0.5f * j;
        const float v = half_sqrt3 * j;
        HS_EXPECT_EQ(at(u, v), 1.0f);
        HS_EXPECT_EQ(at(u + 0.5f, v), -1.0f);
        HS_EXPECT_EQ(at(u + 0.25f, v + 0.5f * half_sqrt3), -1.0f);
        HS_EXPECT_EQ(at(u - 0.25f, v + 0.5f * half_sqrt3), -1.0f);
        HS_EXPECT_LT(at(u + 0.42f, v), 0.0f);
        const float corner = at(u + 0.44f * half_sqrt3, v + 0.44f * 0.5f);
        if (shape == 0.0f)
          HS_EXPECT_LT(corner, 0.0f);
        else
          HS_EXPECT_GT(corner, 0.0f);
      }
  }
}

/** @brief Field of one fixed tessellation at a plane point. */
template <Pattern::TessellationKind Kind, Pattern::TessellationFill Fill>
float fixed_tessellation(
    const Pullback::Source::FixedTessellationParams<Kind, Fill> &params,
    const math::Complex &point) {
  return Pattern::tessellation<Kind, Fill>(Pattern::PlanePosition{point},
                                           params);
}

/** @brief One fixed tessellation is invariant under its reported basis and
    not under half of it. */
template <Pattern::TessellationKind Kind, Pattern::TessellationFill Fill>
void check_tessellation_translation_lattice() {
  HS_CONTEXT("kind / fill", static_cast<int>(Kind), static_cast<int>(Fill));
  Pullback::Source::FixedTessellationParams<Kind, Fill> params;
  params.cell_scale = 1.7f;
  if constexpr (Fill == Pattern::TessellationFill::CELLS) {
    params.cell_shape = 0.8f;
    params.cell_radius = 0.35f;
    params.cell_softness = 0.1f;
  } else {
    params.line_thickness = 0.06f;
    params.line_softness = 0.05f;
  }
  const Pullback::TranslationLattice lattice =
      Pullback::Source::translation_lattice(params);
  const float period = 1.0f / 1.7f;
  const bool square = Kind == Pattern::TessellationKind::SQUARE;
  HS_EXPECT_EQ(lattice.a.re, period);
  HS_EXPECT_EQ(lattice.a.im, 0.0f);
  HS_EXPECT_NEAR(lattice.b.re, square ? 0.0f : 0.5f * period, 1e-7f);
  HS_EXPECT_NEAR(lattice.b.im, square ? period : 0.8660254038f * period, 1e-7f);
  bool half_differs = false;
  hs::Pcg32 random(0x2545f491u);
  for (int sample = 0; sample < 300; ++sample) {
    const float x = rand_uniform(random, -2.0f, 2.0f);
    const float y = rand_uniform(random, -2.0f, 2.0f);
    const math::Complex point(x, y);
    const float base = fixed_tessellation(params, point);
    for (const math::Complex shift :
         {lattice.a, lattice.b, lattice.a - lattice.b,
          math::Complex() - lattice.a - lattice.b})
      HS_EXPECT_NEAR(fixed_tessellation(params, point + shift), base, 1e-4f);
    const math::Complex half(0.5f * lattice.a.re, 0.5f * lattice.a.im);
    half_differs |=
        fabsf(fixed_tessellation(params, point + half) - base) > 0.5f;
  }
  HS_EXPECT_TRUE(half_differs);
}

/** @brief Each tessellation kind reports a translation lattice that leaves
    its field unchanged. */
inline void test_pullback_tessellation_translation_lattice() {
  using K = Pattern::TessellationKind;
  using F = Pattern::TessellationFill;
  check_tessellation_translation_lattice<K::TRIANGULAR, F::CELLS>();
  check_tessellation_translation_lattice<K::SQUARE, F::CELLS>();
  check_tessellation_translation_lattice<K::HEXAGONAL, F::CELLS>();
  check_tessellation_translation_lattice<K::TRIANGULAR, F::EDGES>();
  check_tessellation_translation_lattice<K::SQUARE, F::EDGES>();
  check_tessellation_translation_lattice<K::HEXAGONAL, F::EDGES>();
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
    float cell_scale;
    float cell_shape;
    float cell_softness;
    float cell_radius;
  };
  constexpr Prepared prepared{0.0f, 0.0f, 0.0f, 1.0f, 0.0f};
  constexpr Prepared shifted{0.7f, 0.3f, 0.0f, 1.0f, 0.0f};
  constexpr Prepared turned{0.0f, 0.0f, math::PI_F * 0.5f, 0.0f, 1.0f};
  constexpr Params params{1.0f, 0.0f, 1.0f, 0.0f, 0.1f, 0.25f};
  constexpr Params coupled{0.0f, 0.5f, 1.0f, 0.0f, 0.1f, 0.25f};
  const math::Complex origin;
  const math::Complex off(0.6f, 0.25f); // radius 0.65, azimuth atan2(0.25, 0.6)
  HS_EXPECT_EQ(Pattern::twin_wave(Pattern::PlanePosition{origin}, prepared),
               0.0f);
  HS_EXPECT_EQ(Pattern::twin_wave(Pattern::PlanePosition{origin}, shifted),
               0.5f * (math::fast_sinf(0.7f) + math::fast_sinf(0.3f)));
  HS_EXPECT_NEAR(Pattern::twin_wave(Pattern::PlanePosition{off}, turned),
                 0.5f * (math::fast_sinf(0.6f) + math::fast_sinf(0.25f)),
                 2e-3f);
  HS_EXPECT_EQ(Pattern::spiral(Pattern::PlanePosition{origin}, 0.0f), 0.0f);
  HS_EXPECT_NEAR(Pattern::spiral(Pattern::PlanePosition{off}, 0.0f),
                 math::fast_sinf(0.65f), 2e-3f);
  HS_EXPECT_NEAR(Pattern::spiral(Pattern::PlanePosition{origin}, 3.0f), 1.0f,
                 2e-3f);
  HS_EXPECT_NEAR(Pattern::spiral(Pattern::PlanePosition{off}, 3.0f),
                 math::fast_sinf(0.65f - 3.0f * math::fast_atan2(0.25f, 0.6f)),
                 2e-3f);
  HS_EXPECT_EQ(Pattern::grid(Pattern::PlanePosition{origin}, params, prepared),
               0.0f);
  HS_EXPECT_NEAR(
      Pattern::grid(Pattern::PlanePosition{origin}, coupled, shifted),
      math::fast_sinf(0.5f * math::fast_sinf(0.7f)) *
          math::fast_cosf(0.5f * math::fast_cosf(-0.7f)),
      2e-3f);
  HS_EXPECT_EQ(Pattern::square_cells(Pattern::PlanePosition{origin}, params),
               1.0f);
  // Cell coordinate (0.2, 0.1) sits inside the softness band, so the edge ramp
  // resolves between its saturated ends.
  HS_EXPECT_NEAR(Pattern::square_cells(
                     Pattern::PlanePosition{math::Complex(0.2f, 0.1f)}, params),
                 0.38671f, 2e-3f);

  using Kind = Pattern::TessellationKind;
  using Fill = Pattern::TessellationFill;
  const Pullback::Source::TessellationSourceParams tessellation_params;
  HS_EXPECT_EQ((Pattern::tessellation<Kind::SQUARE, Fill::EDGES>(
                   Pattern::PlanePosition{math::Complex(0.5f, 0.0f)},
                   tessellation_params)),
               1.0f);
  HS_EXPECT_EQ((Pattern::tessellation<Kind::SQUARE, Fill::EDGES>(
                   Pattern::PlanePosition{origin}, tessellation_params)),
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
  HS_EXPECT_GT(std::abs(mirror_mid.mirror.offset_x - mirror_0.mirror.offset_x),
               0.1f);
  HS_EXPECT_NEAR(mirror_0.mirror.offset_x, mirror_1.mirror.offset_x, 1e-6f);
  HS_EXPECT_NEAR(mirror_0.mirror.offset_y, mirror_1.mirror.offset_y, 1e-6f);

  Pullback::Warp::VectorNoiseParams vector;
  vector.vector_angle = 0.4f;
  const auto vector_0 = Pullback::Warp::prepare(vector, 0.0f);
  const auto vector_1 = Pullback::Warp::prepare(vector, 1.0f);
  const auto vector_mid = Pullback::Warp::prepare(vector, 0.25f);
  HS_EXPECT_GT(
      (vector_mid.noise_loop.offset - vector_0.noise_loop.offset).length(),
      0.1f);
  HS_EXPECT_NEAR(vector_0.noise_loop.offset.x, vector_1.noise_loop.offset.x,
                 1e-6f);
  HS_EXPECT_NEAR(vector_0.noise_loop.offset.y, vector_1.noise_loop.offset.y,
                 1e-6f);
  HS_EXPECT_NEAR(vector_0.noise_loop.offset.z, vector_1.noise_loop.offset.z,
                 1e-6f);
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
  HS_EXPECT_GT(std::abs(vortex_mid.vortex.center_x - vortex_0.vortex.center_x),
               0.1f);
  HS_EXPECT_NEAR(vortex_0.vortex.center_x, vortex_1.vortex.center_x, 1e-6f);
  HS_EXPECT_NEAR(vortex_0.vortex.center_y, vortex_1.vortex.center_y, 1e-6f);

  Pullback::Warp::AffineParams affine;
  affine.translation_x = 2.0f;
  affine.translation_y = -1.0f;
  affine.scale_x = 1.5f;
  affine.scale_y = 0.75f;
  affine.shear = 0.2f;
  const auto affine_0 = Pullback::Warp::prepare(affine, 0.0f, 0.6f,
                                                Pullback::square_lattice(2.0f));
  const auto affine_1 = Pullback::Warp::prepare(affine, 1.0f, 0.6f,
                                                Pullback::square_lattice(2.0f));
  const auto affine_mid = Pullback::Warp::prepare(
      affine, 0.25f, 0.6f, Pullback::square_lattice(2.0f));
  HS_EXPECT_GT(
      std::abs(affine_mid.affine.translation_x - affine_0.affine.translation_x),
      0.1f);
  HS_EXPECT_NEAR(affine_0.affine.translation_x, affine_1.affine.translation_x,
                 1e-6f);
  HS_EXPECT_NEAR(affine_0.affine.translation_y, affine_1.affine.translation_y,
                 1e-6f);
  HS_EXPECT_NEAR(affine_0.affine.inverse_scale_x,
                 affine_1.affine.inverse_scale_x, 1e-6f);
  HS_EXPECT_NEAR(affine_0.affine.inverse_scale_y,
                 affine_1.affine.inverse_scale_y, 1e-6f);
  HS_EXPECT_NEAR(affine_0.affine.shear, affine_1.affine.shear, 1e-6f);
}

/** @brief Each affine form matches the general form at identity values
    of the omitted parts.
 * @tparam Form Form under test. */
template <Pullback::Warp::AffineForm Form>
inline void check_affine_form_matches_general() {
  using namespace Pullback::Warp;
#if defined(HS_TEST_FAST_MATH)
  constexpr float TOLERANCE = 1e-6f;
#else
  constexpr float TOLERANCE = 0.0f;
#endif
  const auto check = [](const math::Complex &actual,
                        const math::Complex &expected) {
    HS_EXPECT_NEAR(actual.re, expected.re, TOLERANCE);
    HS_EXPECT_NEAR(actual.im, expected.im, TOLERANCE);
  };
  static volatile float knobs[] = {0.6f, 2.0f, -1.0f, 1.5f, 0.75f, 0.2f, 3.0f};
  static volatile float points[] = {0.0f, -0.0f, 0.37f, -1.25f, 3.5f, -0.0625f};
  static volatile float phases[] = {0.0f, 0.3f, 0.999f, 1.25f};
  AffineParams params;
  const float rotation = includes(Form, AffineForm::ROTATES) ? knobs[0] : 0.0f;
  if constexpr (includes(Form, AffineForm::TRANSLATES)) {
    params.translation_x = knobs[1];
    params.translation_y = knobs[2];
  }
  if constexpr (includes(Form, AffineForm::SCALES_SHEARS)) {
    params.scale_x = knobs[3];
    params.scale_y = knobs[4];
    params.shear = knobs[5];
  }
  params.repeat_cells = knobs[6];
  HS_CONTEXT("form", static_cast<int>(Form));
  for (const Pullback::TranslationLattice lattice :
       {Pullback::TranslationLattice{}, Pullback::square_lattice(knobs[3]),
        Pullback::TranslationLattice{{knobs[3], 0.0f}, {knobs[4], knobs[1]}}})
    for (const float phase : phases) {
      const auto general = prepare(params, phase, rotation, lattice);
      const auto narrowed =
          prepare<Form, false>(params, phase, rotation, lattice);
      const auto repeated =
          prepare<Form, true>(params, phase, rotation, lattice);
      const auto &expected = general.affine.repeat;
      const auto &actual = repeated.affine.repeat;
      HS_EXPECT_TRUE(actual.a == expected.a && actual.b == expected.b &&
                     actual.dual_a == expected.dual_a &&
                     actual.dual_b == expected.dual_b &&
                     actual.reduction == expected.reduction);
      for (const float x : points)
        for (const float y : points) {
          const math::Complex input{x, y};
          const math::Complex want = affine_coords(input, general);
          const math::Complex got = affine_coords<Form>(input, narrowed);
          check(got, want);
          check(affine_coords<Form>(input, repeated), want);
        }
    }
}

/** @brief Every affine form matches the general form at identity values of
    its omitted parts. */
inline void test_pullback_affine_forms_match_general() {
  using Pullback::Warp::AffineForm;
  check_affine_form_matches_general<AffineForm::IDENTITY>();
  check_affine_form_matches_general<AffineForm::ROTATES>();
  check_affine_form_matches_general<AffineForm::TRANSLATES>();
  check_affine_form_matches_general<AffineForm::SCALES_SHEARS>();
  check_affine_form_matches_general<AffineForm::ROTATES |
                                    AffineForm::TRANSLATES>();
  check_affine_form_matches_general<AffineForm::ROTATES |
                                    AffineForm::SCALES_SHEARS>();
  check_affine_form_matches_general<AffineForm::TRANSLATES |
                                    AffineForm::SCALES_SHEARS>();
  check_affine_form_matches_general<AffineForm::GENERAL>();
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
    const float X = random();
    const float Y = random();
    const float Z = random();
    const auto INPUT = math::Vector(X, Y, Z).normalized();
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
              const math::Vector VIEW = math::fib_spiral(24, 0.5f, i);
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

/** @brief Owned masks match unmasked radial output through supported and rejected updates. */
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
  ScratchScope scratch(scratch_arena_a);
  auto *sentinel = scratch_arena_a.make<uint32_t>(0x12345678u);
  const size_t OFFSET = scratch_arena_a.get_offset();
  const auto CHECK = [&](const BurstRippleParams &params,
                         const RippleState &state) {
    prepare_radial_tables(params, state, cached, true);
    HS_EXPECT_EQ(scratch_arena_a.get_offset(), OFFSET);
    HS_EXPECT_EQ(*sentinel, 0x12345678u);
    fresh.enabled = false;
    prepare_radial_tables(params, state, fresh, true);
    HS_EXPECT_EQ(scratch_arena_a.get_offset(), OFFSET);
    HS_EXPECT_EQ(*sentinel, 0x12345678u);
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

  for (bool change_shape : {false, true}) {
    cached.enabled = false;
    params = {};
    params.period = 143.0f;
    reordered = {};
    reordered.burst_count = 1;
    reordered.waves[0] = {40.0f, 0};
    prepare_radial_tables(params, reordered, cached, true);
    HS_EXPECT_EQ(cached.count, uint8_t{1});
    constexpr float SENTINEL = 1234.0f;
    cached.groups[0].current[0].angle = SENTINEL;
    advance_ripple(reordered, params);
    if (change_shape)
      params.strength = 0.125f;
    prepare_radial_tables(params, reordered, cached, true);
    if (change_shape)
      HS_EXPECT_NE(cached.groups[0].previous[0].angle, SENTINEL);
    else
      HS_EXPECT_EQ(cached.groups[0].previous[0].angle, SENTINEL);
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
        const math::Vector VIEW = math::fib_spiral(32, 0.5f, i);
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
  using Control::apply_curve;
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
    /** @brief Coupling phase used by the grid. */
    struct Phases final {
      float primary = math::PI_F * 0.5f; ///< Coupling phase.
    };
    /** @brief Returns unit-coupling grid controls. @return Grid controls. */
    static const Pullback::Source::GridSourceParams &
    params(const FrameState &) {
      static constexpr Pullback::Source::GridSourceParams VALUE = [] {
        Pullback::Source::GridSourceParams value{};
        value.complexity = 1.0f;
        return value;
      }();
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
  HS_EXPECT_NEAR(Grid::sample(INPUT, FRAME, PREPARED), std::sin(1.0f), 2e-3f);
}

/** @brief Hexagonal edge bands agree at corners and shared boundaries. */
inline void test_pullback_hexagonal_edges() {
  using namespace Pullback::Source;
  constexpr float HALF_SQRT_3 = 0.8660254038f;
  constexpr float CIRCUMRADIUS = 0.5773502692f;
  math::Complex vertices[6];
  for (int vertex = 0; vertex < 6; ++vertex) {
    const float angle = math::PI_F / 6.0f + vertex * math::PI_F / 3.0f;
    vertices[vertex] = {CIRCUMRADIUS * std::cos(angle),
                        CIRCUMRADIUS * std::sin(angle)};
  }
  for (const math::Complex center :
       {math::Complex(0.0f, 0.0f), math::Complex(1.0f, 0.0f),
        math::Complex(0.5f, HALF_SQRT_3),
        math::Complex(-1.0f, 2.0f * HALF_SQRT_3)}) {
    for (float scale : {0.5f, 1.0f, 3.0f}) {
      TessellationSourceParams params;
      params.cell_scale = scale;
      auto sample = [&](const math::Complex &point) {
        const math::Complex input((center.re + point.re) / scale,
                                  (center.im + point.im) / scale);
        return Pattern::tessellation<Pattern::TessellationKind::HEXAGONAL,
                                     Pattern::TessellationFill::EDGES>(
            Pattern::PlanePosition{input}, params);
      };
      HS_EXPECT_EQ(sample(math::Complex(0.0f, 0.0f)), -1.0f);
      for (int edge = 0; edge < 6; ++edge) {
        const math::Complex &a = vertices[edge];
        const math::Complex &b = vertices[(edge + 1) % 6];
        const math::Complex midpoint((a.re + b.re) * 0.5f,
                                     (a.im + b.im) * 0.5f);
        const float LENGTH = std::hypot(midpoint.re, midpoint.im);
        HS_EXPECT_NEAR(LENGTH, 0.5f, 1e-6f);
        const math::Complex normal(midpoint.re / LENGTH, midpoint.im / LENGTH);
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

#include "tests/pullback/color.h"

/** @brief Runs pullback tests. */
inline int run_pullback_tests() {
  ModuleFixture fixture("pullback");
  HS_RUN(pattern_tests::test_pattern_direction_admission);
  HS_RUN(pattern_tests::test_pattern_sphere_formulas);
  HS_RUN(pattern_tests::test_pattern_source_clocks);
  HS_RUN(pattern_tests::test_pattern_noise_formulas);
  HS_RUN(pattern_tests::test_pattern_plane_bands);
  HS_RUN(pattern_tests::test_pattern_plane_source_protocols);
  HS_RUN(pattern_tests::test_pattern_sphere_cell_layouts);
  HS_RUN(pattern_tests::test_pattern_sphere_cell_seed_isolation);
  HS_RUN(pattern_tests::test_pattern_sphere_cell_outputs);
  HS_RUN(pattern_tests::test_pattern_sphere_cell_all_bisectors);
  HS_RUN(pattern_tests::test_pattern_sphere_cell_admission);
  HS_RUN(test_pullback_runtime_seed_contract);
  HS_RUN(test_pullback_modulation_field_identity);
  HS_RUN(test_pullback_modulation_field_equality);
  HS_RUN(test_pullback_carrier_contract);
  HS_RUN(test_pullback_registered_curve_subset);
  HS_RUN(test_pullback_projection_traits);
  HS_RUN(test_pullback_validation_predicates);
  HS_RUN(test_pullback_evaluation_order);
  HS_RUN(test_pullback_placement_transparency);
  HS_RUN(test_pullback_type_erased_prepare);
  HS_RUN(test_pullback_public_surface);
  HS_RUN(test_pullback_no_instrumentation);
  HS_RUN(test_pullback_counting_instrumentation);
  HS_RUN(test_pullback_prepared_stage_policies);
  HS_RUN(test_pullback_affine_forms_match_general);
  HS_RUN(test_pullback_stage_combinators);
  HS_RUN(test_pullback_provider_contracts);
  HS_RUN(test_pullback_affine_provider_repeat_contract);
  HS_RUN(test_pullback_concrete_catalog);
  HS_RUN(test_pullback_square_cells_formula);
  HS_RUN(test_pullback_tessellation_reference);
  HS_RUN(test_pullback_hexagonal_cells_pack_hex_lattice);
  HS_RUN(test_pullback_tessellation_translation_lattice);
  HS_RUN(test_pullback_custom_grid_preparation);
  HS_RUN(test_pullback_hexagonal_edges);
  HS_RUN(test_pullback_displacement_overflow);
  HS_RUN(test_pullback_warp_phase_loop);
  HS_RUN(test_pullback_dodecahedral_matrix_fold);
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
  HS_RUN(test_pullback_palette_lab_owned_copy);
  HS_RUN(test_pullback_fused_rgb_carrier_finish);
  HS_RUN(test_pullback_cached_phase_finish);
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
  HS_RUN(test_pullback_unit_hue_blending);
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
