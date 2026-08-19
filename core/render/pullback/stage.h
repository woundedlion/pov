/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"
#include "render/pullback/lens.h"
#include "render/pullback/material.h"
#include "render/pullback/surface.h"

/**
 * @file stage.h
 * @brief Shared carrier kernels and the family-typed stage combinators.
 */

namespace Pullback {

/**
 * @brief Shared carrier kernels: the normative stage transformations as free
 *        functions over the canonical carriers, called by the template
 *        combinators and any erased adapter, so the semantics cannot fork
 *        between execution paths.
 */
namespace Kernel {

__attribute__((always_inline)) inline SphereSample
rotate_dir(const SphereSample &input, const Quaternion &conjugate) {
  return {rotate(input.dir, conjugate), input.path_length};
}

/** @brief Displacement adds to the path accumulator, never replaces it. */
__attribute__((always_inline)) inline SphereSample
displace(const SphereSample &input, const SurfaceResult &step) {
  return {step.sphere, input.path_length + step.path_length};
}

__attribute__((always_inline)) inline SphereSample
lens(const SphereSample &input, const Vector &lensed) {
  return {lensed, input.path_length};
}

/** @brief The Project crossing assembles the plane carrier: the policy's
    coords and provenance embed unchanged; `sphere` is combinator state
    written from the pre-projection point. */
__attribute__((always_inline)) inline PlaneSample
project(const SphereSample &input, const Vector &local,
        const ProjectionResult &result) {
  return {result.coords, result.provenance, local, input.path_length};
}

/** @brief A warp advances coords and path only; provenance and sphere are
    immutable through PLANE endomorphisms. */
__attribute__((always_inline)) inline PlaneSample
warp(const PlaneSample &input, const WarpStepResult &step) {
  return {step.coords, input.provenance, input.sphere,
          input.path_length + step.path_length};
}

/** @brief The Sample crossing ramps the weighted field and consumes the
    provenance coverage into the field carrier. */
__attribute__((always_inline)) inline FieldSample
sample(const PlaneSample &input, float weighted, float coverage) {
  return {Detail::clamp_unit((weighted + 1.0f) * 0.5f),
          coverage * input.provenance.domain_coverage, input.sphere,
          input.path_length};
}

__attribute__((always_inline)) inline FieldSample
transfer(const FieldSample &input, float value) {
  FieldSample output = input;
  output.value = value;
  return output;
}

__attribute__((always_inline)) inline FieldSample
coverage(const FieldSample &input, float factor) {
  FieldSample output = input;
  output.coverage = input.coverage * factor;
  return output;
}

} // namespace Kernel

namespace Detail {

template <typename Policy, typename FrameState>
consteval bool surface_policy_callable() {
  if constexpr (PolicyPrepares<Policy, FrameState>)
    return requires(const Vector &input, const FrameState &frame,
                    const typename Policy::Prepared &prepared) {
      { Policy::apply(input, frame, prepared) } -> std::same_as<SurfaceResult>;
    };
  else
    return requires(const Vector &input, const FrameState &frame) {
      { Policy::apply(input, frame) } -> std::same_as<SurfaceResult>;
    };
}

template <typename Policy, typename FrameState>
consteval bool warp_policy_callable() {
  if constexpr (PolicyPrepares<Policy, FrameState>)
    return requires(
        const Complex &input, const ProjectionProvenance &provenance,
        const FrameState &frame, const typename Policy::Prepared &prepared) {
      {
        Policy::apply(input, provenance, frame, prepared)
      } -> std::same_as<WarpStepResult>;
    };
  else
    return requires(const Complex &input,
                    const ProjectionProvenance &provenance,
                    const FrameState &frame) {
      {
        Policy::apply(input, provenance, frame)
      } -> std::same_as<WarpStepResult>;
    };
}

template <typename Policy, typename FrameState>
consteval bool source_policy_callable() {
  if constexpr (PolicyPrepares<Policy, FrameState>)
    return requires(const PlaneSample &input, const FrameState &frame,
                    const typename Policy::Prepared &prepared) {
      { Policy::sample(input, frame, prepared) } -> std::same_as<float>;
    };
  else
    return requires(const PlaneSample &input, const FrameState &frame) {
      { Policy::sample(input, frame) } -> std::same_as<float>;
    };
}

} // namespace Detail

namespace Stage {

/**
 * @brief SPHERE endomorphism: rotates the view by the provider's conjugate.
 * @tparam OrientationProvider Provider exposing `conjugate(frame)`.
 */
template <typename OrientationProvider>
struct Rotate
    : Contract<Rotate<OrientationProvider>, SphereSample, SphereSample> {
  using Policies = std::tuple<>;
  using Provider = OrientationProvider;

  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<OrientationProvider, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        {
          OrientationProvider::conjugate(frame)
        } -> std::same_as<const Quaternion &>;
      };

  template <typename Binding>
  __attribute__((always_inline)) static SphereSample
  run(const SphereSample &input, const typename Binding::FrameState &frame,
      const NoPrepared &) {
    return Kernel::rotate_dir(input, OrientationProvider::conjugate(frame));
  }
};

/**
 * @brief SPHERE endomorphism: displaces the view point on the sphere.
 * @tparam SurfacePolicyT Surface policy returning a SurfaceResult step.
 */
template <typename SurfacePolicyT>
struct Displace
    : Contract<Displace<SurfacePolicyT>, SphereSample, SphereSample> {
  using Policies = std::tuple<SurfacePolicyT>;
  using SurfacePolicy = SurfacePolicyT;

  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::surface_policy_callable<SurfacePolicyT,
                                      typename Binding::FrameState>();

  template <typename Binding>
  HS_FLASH_INLINE static auto
  prepare(const typename Binding::FrameState &frame) {
    return Detail::prepare_policy<SurfacePolicyT>(frame);
  }

  template <typename Binding>
  __attribute__((always_inline)) static SphereSample
  run(const SphereSample &input, const typename Binding::FrameState &frame,
      const Detail::PolicyPrepared<SurfacePolicyT, typename Binding::FrameState>
          &prepared) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto start = Instrumentation::mark();
    SurfaceResult step;
    if constexpr (Detail::PolicyPrepares<SurfacePolicyT,
                                         typename Binding::FrameState>)
      step = SurfacePolicyT::apply(input.dir, frame, prepared);
    else
      step = SurfacePolicyT::apply(input.dir, frame);
    Instrumentation::template span<ProfileEvent::SURFACE_NOISE>(start);
    return Kernel::displace(input, step);
  }
};

/**
 * @brief SPHERE endomorphism: applies a sphere-to-sphere lens.
 * @tparam LensPolicyT Lens policy mapping a direction to a direction.
 */
template <typename LensPolicyT>
struct Lens : Contract<Lens<LensPolicyT>, SphereSample, SphereSample> {
  using Policies = std::tuple<LensPolicyT>;
  using LensPolicy = LensPolicyT;

  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      requires(const Vector &input, const typename Binding::FrameState &frame) {
        { LensPolicyT::apply(input, frame) } -> std::same_as<Vector>;
      };

  template <typename Binding>
  __attribute__((always_inline)) static SphereSample
  run(const SphereSample &input, const typename Binding::FrameState &frame,
      const NoPrepared &) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto start = Instrumentation::mark();
    const Vector lensed = LensPolicyT::apply(input.dir, frame);
    Instrumentation::template span<ProfileEvent::LENS>(start);
    return Kernel::lens(input, lensed);
  }
};

/**
 * @brief SPHERE→PLANE crossing: rotates into the projection frame, projects,
 *        and assembles the plane carrier.
 * @details The projection protocol has no sphere field: the sample point is
 * combinator state, written from the pre-projection point — a policy cannot
 * even accidentally control it.
 * @tparam ProjectionPolicyT Projection policy returning a ProjectionResult.
 */
template <typename ProjectionPolicyT>
struct Project
    : Contract<Project<ProjectionPolicyT>, SphereSample, PlaneSample> {
  using Policies = std::tuple<ProjectionPolicyT>;
  using ProjectionPolicy = ProjectionPolicyT;

  static constexpr bool EDGE_DISTANCE_UNCONDITIONAL = [] {
    if constexpr (requires { ProjectionPolicyT::EDGE_DISTANCE_UNCONDITIONAL; })
      return ProjectionPolicyT::EDGE_DISTANCE_UNCONDITIONAL;
    else
      return false;
  }();

  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      requires(const Vector &input, const typename Binding::FrameState &frame) {
        {
          ProjectionPolicyT::frame_conjugate(frame)
        } -> std::same_as<const Quaternion &>;
        {
          ProjectionPolicyT::project(input, frame)
        } -> std::same_as<ProjectionResult>;
      };

  template <typename Binding>
  __attribute__((always_inline)) static PlaneSample
  run(const SphereSample &input, const typename Binding::FrameState &frame,
      const NoPrepared &) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto start = Instrumentation::mark();
    const Vector local =
        rotate(input.dir, ProjectionPolicyT::frame_conjugate(frame));
    const ProjectionResult result = ProjectionPolicyT::project(local, frame);
    const PlaneSample output = Kernel::project(input, local, result);
    Instrumentation::template span<ProfileEvent::PROJECTION>(start);
    return output;
  }
};

/**
 * @brief PLANE endomorphism: one planar warp step.
 * @tparam WarpPolicyT Warp policy advancing the working coordinate.
 */
template <typename WarpPolicyT>
struct Warp : Contract<Warp<WarpPolicyT>, PlaneSample, PlaneSample> {
  using Policies = std::tuple<WarpPolicyT>;
  using WarpPolicy = WarpPolicyT;

  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::warp_policy_callable<WarpPolicyT, typename Binding::FrameState>();

  template <typename Binding>
  HS_FLASH_INLINE static auto
  prepare(const typename Binding::FrameState &frame) {
    return Detail::prepare_policy<WarpPolicyT>(frame);
  }

  template <typename Binding>
  __attribute__((always_inline)) static PlaneSample
  run(const PlaneSample &input, const typename Binding::FrameState &frame,
      const Detail::PolicyPrepared<WarpPolicyT, typename Binding::FrameState>
          &prepared) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto start = Instrumentation::mark();
    WarpStepResult step;
    if constexpr (Detail::PolicyPrepares<WarpPolicyT,
                                         typename Binding::FrameState>)
      step =
          WarpPolicyT::apply(input.coords, input.provenance, frame, prepared);
    else
      step = WarpPolicyT::apply(input.coords, input.provenance, frame);
    const PlaneSample output = Kernel::warp(input, step);
    Instrumentation::template span<ProfileEvent::PLANAR_WARP>(start);
    return output;
  }
};

/**
 * @brief PLANE→FIELD crossing: samples the source, weights the raw signed
 *        field, ramps, and seeds the field carrier — establishing the
 *        FieldSample invariant.
 * @details Coverage modes never stack: the crossing's coverage policy is one
 * slot from the ProjectedCoverage vocabulary. The raw signed field exists
 * only inside the crossing, which is why weighting is a crossing policy
 * rather than a stage.
 * @tparam SourcePolicyT Scalar source policy.
 * @tparam WeightPolicyT Signal-weight policy over the raw signed field.
 * @tparam CoveragePolicyT Projected-coverage policy over the provenance.
 */
template <typename SourcePolicyT, typename WeightPolicyT = Weight::Projection,
          typename CoveragePolicyT = ProjectedCoverage::Weight>
struct Sample : Contract<Sample<SourcePolicyT, WeightPolicyT, CoveragePolicyT>,
                         PlaneSample, FieldSample> {
  using Policies = std::tuple<SourcePolicyT, WeightPolicyT, CoveragePolicyT>;
  using SourcePolicy = SourcePolicyT;
  using WeightPolicy = WeightPolicyT;
  using CoveragePolicy = CoveragePolicyT;

  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::source_policy_callable<SourcePolicyT,
                                     typename Binding::FrameState>() &&
      requires(float field, const ProjectionProvenance &provenance,
               const typename Binding::FrameState &frame) {
        {
          WeightPolicyT::apply(field, provenance, frame)
        } -> std::same_as<float>;
        { CoveragePolicyT::apply(provenance, frame) } -> std::same_as<float>;
      };

  template <typename Binding>
  HS_FLASH_INLINE static auto
  prepare(const typename Binding::FrameState &frame) {
    return Detail::prepare_policy<SourcePolicyT>(frame);
  }

  template <typename Binding>
  __attribute__((always_inline)) static FieldSample
  run(const PlaneSample &input, const typename Binding::FrameState &frame,
      const Detail::PolicyPrepared<SourcePolicyT, typename Binding::FrameState>
          &prepared) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto source_span = Instrumentation::mark();
    float raw;
    if constexpr (Detail::PolicyPrepares<SourcePolicyT,
                                         typename Binding::FrameState>)
      raw = SourcePolicyT::sample(input, frame, prepared);
    else
      raw = SourcePolicyT::sample(input, frame);
    Instrumentation::template span<ProfileEvent::SOURCE>(source_span);
    const auto material_span = Instrumentation::mark();
    const float weighted = WeightPolicyT::apply(raw, input.provenance, frame);
    const FieldSample output = Kernel::sample(
        input, weighted, CoveragePolicyT::apply(input.provenance, frame));
    Instrumentation::template span<ProfileEvent::MATERIAL>(material_span);
    return output;
  }
};

/**
 * @brief FIELD endomorphism: reshapes the field value.
 * @tparam TransferPolicyT Transfer policy over the unit value.
 */
template <typename TransferPolicyT>
struct Transfer
    : Contract<Transfer<TransferPolicyT>, FieldSample, FieldSample> {
  using Policies = std::tuple<TransferPolicyT>;
  using TransferPolicy = TransferPolicyT;

  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      requires(float value, const typename Binding::FrameState &frame) {
        { TransferPolicyT::apply(value, frame) } -> std::same_as<float>;
      };

  template <typename Binding>
  __attribute__((always_inline)) static FieldSample
  run(const FieldSample &input, const typename Binding::FrameState &frame,
      const NoPrepared &) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto start = Instrumentation::mark();
    const FieldSample output =
        Kernel::transfer(input, TransferPolicyT::apply(input.value, frame));
    Instrumentation::template span<ProfileEvent::MATERIAL>(start);
    return output;
  }
};

/**
 * @brief FIELD endomorphism: multiplies a value-dependent factor into the
 *        accumulated coverage.
 * @details Value-dependent only — provenance-driven coverage is consumed at
 * the Sample crossing. A chain may place this before, between, or after
 * transfers; each placement reads a different value.
 * @tparam CoveragePolicyT Value-coverage policy, e.g. Coverage::ValueCutout.
 */
template <typename CoveragePolicyT>
struct Coverage
    : Contract<Coverage<CoveragePolicyT>, FieldSample, FieldSample> {
  using Policies = std::tuple<CoveragePolicyT>;
  using CoveragePolicy = CoveragePolicyT;

  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      requires(float value, const typename Binding::FrameState &frame) {
        { CoveragePolicyT::apply(value, frame) } -> std::same_as<float>;
      };

  template <typename Binding>
  __attribute__((always_inline)) static FieldSample
  run(const FieldSample &input, const typename Binding::FrameState &frame,
      const NoPrepared &) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto start = Instrumentation::mark();
    const FieldSample output =
        Kernel::coverage(input, CoveragePolicyT::apply(input.value, frame));
    Instrumentation::template span<ProfileEvent::MATERIAL>(start);
    return output;
  }
};

/**
 * @brief FIELD→COLOR crossing: colorizes the field carrier.
 * @tparam ColorPolicyT Color policy producing straight-alpha Color4.
 */
template <typename ColorPolicyT>
struct Colorize : Contract<Colorize<ColorPolicyT>, FieldSample, Color4> {
  using Policies = std::tuple<ColorPolicyT>;
  using ColorPolicy = ColorPolicyT;

  template <typename Binding>
  static constexpr bool PROVIDER_VALID = requires(
      const FieldSample &input, const typename Binding::FrameState &frame) {
    { ColorPolicyT::apply(input, frame) } -> std::same_as<Color4>;
  };

  template <typename Binding>
  __attribute__((always_inline)) static Color4
  run(const FieldSample &input, const typename Binding::FrameState &frame,
      const NoPrepared &) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto start = Instrumentation::mark();
    const Color4 result = ColorPolicyT::apply(input, frame);
    Instrumentation::template span<ProfileEvent::COLOR>(start);
    return result;
  }
};

} // namespace Stage

} // namespace Pullback
