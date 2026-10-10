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
 *        functions over the canonical carriers.
 */
namespace Kernel {

/**
 * @brief Rotates the direction by @p conjugate; path length is unchanged.
 * @param input Sphere carrier.
 * @param conjugate Rotation applied to the direction.
 * @return The rotated carrier.
 */
__attribute__((always_inline)) inline SphereSample
rotate_dir(const SphereSample &input, const math::Quaternion &conjugate) {
  return {math::rotate(input.dir, conjugate), input.path_length};
}

/** @brief Displacement adds to the path accumulator, never replaces it.
    @param input Sphere carrier.
    @param step Surface step.
    @return The displaced carrier. */
__attribute__((always_inline)) inline SphereSample
displace(const SphereSample &input, const SurfaceResult &step) {
  return {step.sphere, input.path_length + step.path_length};
}

/**
 * @brief Replaces the direction with a lensed one; path length is unchanged.
 * @param input Sphere carrier.
 * @param lensed Lens output direction.
 * @return The lensed carrier.
 */
__attribute__((always_inline)) inline SphereSample
lens(const SphereSample &input, const math::Vector &lensed) {
  return {lensed, input.path_length};
}

/** @brief The Project crossing assembles the plane carrier: the policy's
    coords and provenance embed unchanged; `sphere` is combinator state
    written from the pre-projection point.
    @param input Unit-direction sphere carrier.
    @param local Direction in the projection frame.
    @param result Projection policy output.
    @return The plane carrier. */
__attribute__((always_inline)) inline PlaneSample
project(const SphereSample &input, const math::Vector &local,
        const ProjectionResult &result) {
  HS_AUDIT_CHECK(
      fabsf(input.dir.x * input.dir.x + input.dir.y * input.dir.y +
            input.dir.z * input.dir.z - 1.0f) <= 0.004f,
      "Project requires a unit direction within lens approximation error");
  return {result.coords, result.provenance, local, input.path_length};
}

/** @brief A warp advances coords and path only; provenance and sphere are
    immutable through PLANE endomorphisms.
    @param input Plane carrier.
    @param step Warp step.
    @return The warped carrier. */
__attribute__((always_inline)) inline PlaneSample
warp(const PlaneSample &input, const WarpStepResult &step) {
  return {step.coords, input.provenance, input.sphere,
          input.path_length + step.path_length};
}

/** @brief The Sample crossing ramps the weighted field and consumes the
    provenance coverage into the field carrier.
    @param input Plane carrier.
    @param weighted Weighted signed field in [-1, 1].
    @param coverage Projected coverage in [0, 1].
    @return The field carrier. */
__attribute__((always_inline)) inline FieldSample
sample(const PlaneSample &input, float weighted, float coverage) {
  HS_AUDIT_CHECK(coverage >= 0.0f && coverage <= 1.0f &&
                     input.provenance.domain_coverage >= 0.0f &&
                     input.provenance.domain_coverage <= 1.0f,
                 "field coverage factors must remain in [0, 1]");
  return {Detail::clamp_unit((weighted + 1.0f) * 0.5f),
          coverage * input.provenance.domain_coverage, input.sphere,
          input.path_length};
}

/** @brief The spherical Sample crossing ramps a signed field with opaque
    coverage.
    @param input Sphere carrier.
    @param field Signed field in [-1, 1].
    @return The field carrier. */
__attribute__((always_inline)) inline FieldSample
sample(const SphereSample &input, float field) {
  return {Detail::clamp_unit((field + 1.0f) * 0.5f), 1.0f, input.dir,
          input.path_length};
}

/** @brief A transfer replaces the field value; @p value must already be in
    [0, 1], which this does not re-clamp.
    @param input Field carrier.
    @param value New field value.
    @return The carrier with @p value. */
__attribute__((always_inline)) inline FieldSample
transfer(const FieldSample &input, float value) {
  HS_AUDIT_CHECK(value >= 0.0f && value <= 1.0f,
                 "field transfer must remain in [0, 1]");
  FieldSample output = input;
  output.value = value;
  return output;
}

/**
 * @brief Multiplies @p factor into the accumulated coverage.
 * @param input Field carrier.
 * @param factor Coverage factor in [0, 1].
 * @return @p input with the scaled coverage.
 */
__attribute__((always_inline)) inline FieldSample
coverage(const FieldSample &input, float factor) {
  HS_AUDIT_CHECK(factor >= 0.0f && factor <= 1.0f,
                 "field coverage factor must remain in [0, 1]");
  FieldSample output = input;
  output.coverage = input.coverage * factor;
  return output;
}

} // namespace Kernel

namespace Detail {

/**
 * @brief Whether @p Policy has the warp `apply(coords, provenance, frame)`
 *        signature for @p FrameState, taking its prepared state when it
 *        declares one.
 * @tparam Policy Policy being checked.
 * @tparam FrameState Frame state of the binding.
 * @return True when the call is well formed.
 */
template <typename Policy, typename FrameState>
consteval bool warp_policy_callable() {
  if constexpr (PolicyPrepares<Policy, FrameState>)
    return requires(
        const math::Complex &input, const ProjectionProvenance &provenance,
        const FrameState &frame, const typename Policy::Prepared &prepared) {
      {
        Policy::apply(input, provenance, frame, prepared)
      } -> std::same_as<WarpStepResult>;
    };
  else
    return requires(const math::Complex &input,
                    const ProjectionProvenance &provenance,
                    const FrameState &frame) {
      {
        Policy::apply(input, provenance, frame)
      } -> std::same_as<WarpStepResult>;
    };
}

/**
 * @brief Whether @p Policy has the `sample(input, frame)` signature for
 *        @p FrameState, taking its prepared state when it declares one.
 * @tparam Policy Policy being checked.
 * @tparam Input Carrier the policy reads.
 * @tparam FrameState Frame state of the binding.
 * @return True when the call is well formed.
 */
template <typename Policy, typename Input, typename FrameState>
consteval bool sample_policy_callable() {
  if constexpr (PolicyPrepares<Policy, FrameState>)
    return requires(const Input &input, const FrameState &frame,
                    const typename Policy::Prepared &prepared) {
      { Policy::sample(input, frame, prepared) } -> std::same_as<float>;
    };
  else
    return requires(const Input &input, const FrameState &frame) {
      { Policy::sample(input, frame) } -> std::same_as<float>;
    };
}

/**
 * @brief Whether @p Policy has the weight `apply(field, provenance, frame)`
 *        signature for @p FrameState, taking its prepared state when it
 *        declares one.
 * @tparam Policy Policy being checked.
 * @tparam FrameState Frame state of the binding.
 * @return True when the call is well formed.
 */
template <typename Policy, typename FrameState>
consteval bool weight_policy_callable() {
  if constexpr (PolicyPrepares<Policy, FrameState>)
    return requires(float field, const ProjectionProvenance &provenance,
                    const FrameState &frame,
                    const typename Policy::Prepared &prepared) {
      {
        Policy::apply(field, provenance, frame, prepared)
      } -> std::same_as<float>;
    };
  else
    return requires(float field, const ProjectionProvenance &provenance,
                    const FrameState &frame) {
      { Policy::apply(field, provenance, frame) } -> std::same_as<float>;
    };
}

/**
 * @brief Whether @p Policy has the `conjugate(frame)` signature for
 *        @p FrameState, taking its prepared state when it declares one.
 * @tparam Policy Policy being checked.
 * @tparam FrameState Frame state of the binding.
 * @return True when the call is well formed.
 */
template <typename Policy, typename FrameState>
consteval bool orientation_policy_callable() {
  if constexpr (PolicyPrepares<Policy, FrameState>)
    return requires(const FrameState &frame,
                    const typename Policy::Prepared &prepared) {
      {
        Policy::conjugate(frame, prepared)
      } -> std::same_as<const math::Quaternion &>;
    };
  else
    return requires(const FrameState &frame) {
      { Policy::conjugate(frame) } -> std::same_as<const math::Quaternion &>;
    };
}

/**
 * @brief Whether @p Policy has the projection `frame_conjugate` and `project`
 *        signature for @p FrameState, taking its prepared state when it
 *        declares one.
 * @tparam Policy Policy being checked.
 * @tparam FrameState Frame state of the binding.
 * @return True when the call is well formed.
 */
template <typename Policy, typename FrameState>
consteval bool projection_policy_callable() {
  if constexpr (PolicyPrepares<Policy, FrameState>)
    return requires(const math::Vector &input, const FrameState &frame,
                    const typename Policy::Prepared &prepared) {
      {
        Policy::frame_conjugate(frame, prepared)
      } -> std::same_as<const math::Quaternion &>;
      {
        Policy::project(input, frame, prepared)
      } -> std::same_as<ProjectionResult>;
    };
  else
    return requires(const math::Vector &input, const FrameState &frame) {
      {
        Policy::frame_conjugate(frame)
      } -> std::same_as<const math::Quaternion &>;
      { Policy::project(input, frame) } -> std::same_as<ProjectionResult>;
    };
}

/**
 * @brief Whether @p Policy has the `apply(input, frame)` signature for
 *        @p FrameState, taking its prepared state when it declares one.
 * @tparam Policy Policy being checked.
 * @tparam Input Carrier the policy reads.
 * @tparam Output Result type `apply` must return.
 * @tparam FrameState Frame state of the binding.
 * @return True when the call is well formed.
 */
template <typename Policy, typename Input, typename Output, typename FrameState>
consteval bool apply_policy_callable() {
  if constexpr (PolicyPrepares<Policy, FrameState>)
    return requires(const Input &input, const FrameState &frame,
                    const typename Policy::Prepared &prepared) {
      { Policy::apply(input, frame, prepared) } -> std::same_as<Output>;
    };
  else
    return requires(const Input &input, const FrameState &frame) {
      { Policy::apply(input, frame) } -> std::same_as<Output>;
    };
}

/** @brief Whether @p Policy declares value role @p Role.
    @tparam Policy Policy being checked.
    @tparam Role Expected role.
    @return True when it does; false when it declares no `VALUE_ROLE`. */
template <typename Policy, ValueRole Role> consteval bool value_role_is() {
  if constexpr (requires { Policy::VALUE_ROLE; })
    return Policy::VALUE_ROLE == Role;
  else
    return false;
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
  /// Policies the stage binds.
  using Policies = std::tuple<OrientationProvider>;
  using Provider = OrientationProvider; ///< The orientation provider.

  /**
   * @brief Whether `OrientationProvider` is a provider with a callable
   *        `conjugate` under @p Binding.
   * @tparam Binding Binding being checked.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<OrientationProvider, Binding> &&
      Detail::orientation_policy_callable<OrientationProvider,
                                          typename Binding::FrameState>();

  /**
   * @brief Resolves the orientation provider's per-frame state.
   * @tparam Binding Binding of the pipeline.
   * @param frame Frame state.
   * @return The prepared state, or NoPrepared when the policy declares none.
   */
  template <typename Binding>
  HS_FLASH_INLINE static auto
  prepare(const typename Binding::FrameState &frame) {
    return Detail::prepare_policy<OrientationProvider>(frame);
  }

  /**
   * @brief Rotates the view direction by the provider's conjugate.
   * @tparam Binding Binding of the pipeline.
   * @param input Sphere carrier.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return The rotated carrier.
   */
  template <typename Binding>
  __attribute__((always_inline)) static SphereSample
  run(const SphereSample &input, const typename Binding::FrameState &frame,
      const Detail::PolicyPrepared<OrientationProvider,
                                   typename Binding::FrameState> &prepared) {
    if constexpr (Detail::PolicyPrepares<OrientationProvider,
                                         typename Binding::FrameState>)
      return Kernel::rotate_dir(
          input, OrientationProvider::conjugate(frame, prepared));
    else
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
  using Policies = std::tuple<SurfacePolicyT>; ///< Policies the stage binds.
  using SurfacePolicy = SurfacePolicyT;        ///< The surface policy.

  /**
   * @brief Whether `SurfacePolicyT::apply` is callable under @p Binding.
   * @tparam Binding Binding being checked.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::apply_policy_callable<SurfacePolicyT, math::Vector, SurfaceResult,
                                    typename Binding::FrameState>();

  /**
   * @brief Resolves the surface policy's per-frame state.
   * @tparam Binding Binding of the pipeline.
   * @param frame Frame state.
   * @return The prepared state, or NoPrepared when the policy declares none.
   */
  template <typename Binding>
  HS_FLASH_INLINE static auto
  prepare(const typename Binding::FrameState &frame) {
    return Detail::prepare_policy<SurfacePolicyT>(frame);
  }

  /**
   * @brief Displaces the view point and accumulates the step's path
   *        length.
   * @tparam Binding Binding of the pipeline.
   * @param input Sphere carrier.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return The displaced carrier.
   */
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
    const SphereSample output = Kernel::displace(input, step);
    Instrumentation::template span<ProfileEvent::SURFACE_NOISE>(start);
    return output;
  }
};

/**
 * @brief SPHERE endomorphism: applies a sphere-to-sphere lens.
 * @tparam LensPolicyT Lens policy mapping a direction to a direction.
 */
template <typename LensPolicyT>
struct Lens : Contract<Lens<LensPolicyT>, SphereSample, SphereSample> {
  using Policies = std::tuple<LensPolicyT>; ///< Policies the stage binds.
  using LensPolicy = LensPolicyT;           ///< The lens policy.

  /**
   * @brief Whether `LensPolicyT::apply` is callable under @p Binding.
   * @tparam Binding Binding being checked.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::apply_policy_callable<LensPolicyT, math::Vector, math::Vector,
                                    typename Binding::FrameState>();

  /**
   * @brief Resolves the lens policy's per-frame state.
   * @tparam Binding Binding of the pipeline.
   * @param frame Frame state.
   * @return The prepared state, or NoPrepared when the policy declares none.
   */
  template <typename Binding>
  HS_FLASH_INLINE static auto
  prepare(const typename Binding::FrameState &frame) {
    return Detail::prepare_policy<LensPolicyT>(frame);
  }

  /**
   * @brief Applies the lens to the view direction.
   * @tparam Binding Binding of the pipeline.
   * @param input Sphere carrier.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return The lensed carrier.
   */
  template <typename Binding>
  __attribute__((always_inline)) static SphereSample
  run(const SphereSample &input, const typename Binding::FrameState &frame,
      const Detail::PolicyPrepared<LensPolicyT, typename Binding::FrameState>
          &prepared) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto start = Instrumentation::mark();
    math::Vector lensed;
    if constexpr (Detail::PolicyPrepares<LensPolicyT,
                                         typename Binding::FrameState>)
      lensed = LensPolicyT::apply(input.dir, frame, prepared);
    else
      lensed = LensPolicyT::apply(input.dir, frame);
    const SphereSample output = Kernel::lens(input, lensed);
    Instrumentation::template span<ProfileEvent::LENS>(start);
    return output;
  }
};

/**
 * @brief SPHERE→PLANE crossing: rotates into the projection frame, projects,
 *        and assembles the plane carrier.
 * @details The projection protocol has no sphere field: the sample point is
 * combinator state, written from the pre-projection point.
 * @tparam ProjectionPolicyT Projection policy returning a ProjectionResult.
 */
template <typename ProjectionPolicyT>
struct Project
    : Contract<Project<ProjectionPolicyT>, SphereSample, PlaneSample> {
  using Policies = std::tuple<ProjectionPolicyT>; ///< Policies the stage binds.
  using ProjectionPolicy = ProjectionPolicyT;     ///< The projection policy.

  /** @brief Whether the projection policy measures `fade_edge_distance`;
      false when it does not declare `EDGE_DISTANCE_AVAILABLE`. */
  static constexpr bool EDGE_DISTANCE_AVAILABLE = [] {
    if constexpr (requires { ProjectionPolicyT::EDGE_DISTANCE_AVAILABLE; })
      return ProjectionPolicyT::EDGE_DISTANCE_AVAILABLE;
    else
      return false;
  }();

  /**
   * @brief Whether the projection policy's `frame_conjugate` and `project`
   *        are callable under @p Binding.
   * @tparam Binding Binding being checked.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::projection_policy_callable<ProjectionPolicyT,
                                         typename Binding::FrameState>();

  /**
   * @brief Resolves the projection policy's per-frame state.
   * @tparam Binding Binding of the pipeline.
   * @param frame Frame state.
   * @return The prepared state, or NoPrepared when the policy declares none.
   */
  template <typename Binding>
  HS_FLASH_INLINE static auto
  prepare(const typename Binding::FrameState &frame) {
    return Detail::prepare_policy<ProjectionPolicyT>(frame);
  }

  /**
   * @brief Rotates into the projection frame and projects.
   * @tparam Binding Binding of the pipeline.
   * @param input Sphere carrier.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Plane carrier recording the projection-frame
   *         point.
   */
  template <typename Binding>
  __attribute__((always_inline)) static PlaneSample
  run(const SphereSample &input, const typename Binding::FrameState &frame,
      const Detail::PolicyPrepared<ProjectionPolicyT,
                                   typename Binding::FrameState> &prepared) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto start = Instrumentation::mark();

    math::Vector local;
    ProjectionResult result;
    if constexpr (Detail::PolicyPrepares<ProjectionPolicyT,
                                         typename Binding::FrameState>) {
      local = math::rotate(input.dir,
                           ProjectionPolicyT::frame_conjugate(frame, prepared));
      result = ProjectionPolicyT::project(local, frame, prepared);
    } else {
      local =
          math::rotate(input.dir, ProjectionPolicyT::frame_conjugate(frame));
      result = ProjectionPolicyT::project(local, frame);
    }
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
  using Policies = std::tuple<WarpPolicyT>; ///< Policies the stage binds.
  using WarpPolicy = WarpPolicyT;           ///< The warp policy.

  /**
   * @brief Whether `WarpPolicyT::apply` is callable under @p Binding.
   * @tparam Binding Binding being checked.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::warp_policy_callable<WarpPolicyT, typename Binding::FrameState>();

  /**
   * @brief Resolves the warp policy's per-frame state.
   * @tparam Binding Binding of the pipeline.
   * @param frame Frame state.
   * @return The prepared state, or NoPrepared when the policy declares none.
   */
  template <typename Binding>
  HS_FLASH_INLINE static auto
  prepare(const typename Binding::FrameState &frame) {
    return Detail::prepare_policy<WarpPolicyT>(frame);
  }

  /**
   * @brief Advances the plane coordinate by one warp step.
   * @tparam Binding Binding of the pipeline.
   * @param input Plane carrier.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return The warped carrier.
   */
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
 * slot from the ProjectionCoverage vocabulary.
 * @tparam SourcePolicyT Scalar source policy.
 * @tparam WeightPolicyT Signal-weight policy over the raw signed field.
 * @tparam CoveragePolicyT Projected-coverage policy over the provenance.
 */
template <typename SourcePolicyT, typename WeightPolicyT = Weight::Projection,
          typename CoveragePolicyT = ProjectionCoverage::Weight>
struct Sample : Contract<Sample<SourcePolicyT, WeightPolicyT, CoveragePolicyT>,
                         PlaneSample, FieldSample> {
  /// Policies the stage binds.
  using Policies = std::tuple<SourcePolicyT, WeightPolicyT, CoveragePolicyT>;
  using SourcePolicy = SourcePolicyT;     ///< The scalar source policy.
  using WeightPolicy = WeightPolicyT;     ///< The signal-weight policy.
  using CoveragePolicy = CoveragePolicyT; ///< The projected-coverage policy.

  /**
   * @brief Whether the source, weight and coverage policies are callable under
   *        @p Binding.
   * @tparam Binding Binding being checked.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::sample_policy_callable<SourcePolicyT, PlaneSample,
                                     typename Binding::FrameState>() &&
      Detail::weight_policy_callable<WeightPolicyT,
                                     typename Binding::FrameState>() &&
      Detail::apply_policy_callable<CoveragePolicyT, ProjectionProvenance,
                                    float, typename Binding::FrameState>();

  /**
   * @brief Whether the weight or coverage policy prepares state under
   *        @p Binding.
   * @tparam Binding Binding being checked.
   */
  template <typename Binding>
  static constexpr bool MATERIAL_PREPARES =
      Detail::PolicyPrepares<WeightPolicyT, typename Binding::FrameState> ||
      Detail::PolicyPrepares<CoveragePolicyT, typename Binding::FrameState>;

  /**
   * @brief Prepared state: a (source, weight, coverage) tuple when
   *        `MATERIAL_PREPARES`, else the source's prepared state alone.
   * @tparam Binding Binding of the pipeline.
   */
  template <typename Binding>
  using Prepared = std::conditional_t<
      MATERIAL_PREPARES<Binding>,
      std::tuple<
          Detail::PolicyPrepared<SourcePolicyT, typename Binding::FrameState>,
          Detail::PolicyPrepared<WeightPolicyT, typename Binding::FrameState>,
          Detail::PolicyPrepared<CoveragePolicyT,
                                 typename Binding::FrameState>>,
      Detail::PolicyPrepared<SourcePolicyT, typename Binding::FrameState>>;

  /**
   * @brief Resolves the per-frame state of the stage's policies.
   * @tparam Binding Binding of the pipeline.
   * @param frame Frame state.
   * @return The stage's `Prepared` state.
   */
  template <typename Binding>
  HS_FLASH_INLINE static auto
  prepare(const typename Binding::FrameState &frame) {
    if constexpr (MATERIAL_PREPARES<Binding>)
      return Prepared<Binding>{Detail::prepare_policy<SourcePolicyT>(frame),
                               Detail::prepare_policy<WeightPolicyT>(frame),
                               Detail::prepare_policy<CoveragePolicyT>(frame)};
    else
      return Detail::prepare_policy<SourcePolicyT>(frame);
  }

  /**
   * @brief Samples the source, weights it and seeds the field carrier.
   * @tparam Binding Binding of the pipeline.
   * @param input Plane carrier.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field carrier with the ramped value and
   *         coverage.
   */
  template <typename Binding>
  __attribute__((always_inline)) static FieldSample
  run(const PlaneSample &input, const typename Binding::FrameState &frame,
      const Prepared<Binding> &prepared) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto source_span = Instrumentation::mark();
    float raw;
    if constexpr (Detail::PolicyPrepares<SourcePolicyT,
                                         typename Binding::FrameState>) {
      if constexpr (MATERIAL_PREPARES<Binding>)
        raw = SourcePolicyT::sample(input, frame, std::get<0>(prepared));
      else
        raw = SourcePolicyT::sample(input, frame, prepared);
    } else
      raw = SourcePolicyT::sample(input, frame);
    Instrumentation::template span<ProfileEvent::SOURCE>(source_span);
    const auto material_span = Instrumentation::mark();
    float weighted;
    if constexpr (Detail::PolicyPrepares<WeightPolicyT,
                                         typename Binding::FrameState>)
      weighted = WeightPolicyT::apply(raw, input.provenance, frame,
                                      std::get<1>(prepared));
    else
      weighted = WeightPolicyT::apply(raw, input.provenance, frame);
    float coverage;
    if constexpr (Detail::PolicyPrepares<CoveragePolicyT,
                                         typename Binding::FrameState>)
      coverage = CoveragePolicyT::apply(input.provenance, frame,
                                        std::get<2>(prepared));
    else
      coverage = CoveragePolicyT::apply(input.provenance, frame);
    const FieldSample output = Kernel::sample(input, weighted, coverage);
    Instrumentation::template span<ProfileEvent::MATERIAL>(material_span);
    return output;
  }
};

/**
 * @brief SPHERE→FIELD crossing: samples a signed spherical field, ramps it,
 *        and seeds an opaque field carrier.
 * @tparam SourcePolicyT Scalar source policy over SphereSample.
 */
template <typename SourcePolicyT>
struct SampleSphere
    : Contract<SampleSphere<SourcePolicyT>, SphereSample, FieldSample> {
  using Policies = std::tuple<SourcePolicyT>; ///< Policies the stage binds.
  using SourcePolicy = SourcePolicyT;         ///< The scalar source policy.

  /**
   * @brief Whether `SourcePolicyT::sample` is callable on a SphereSample under
   *        @p Binding.
   * @tparam Binding Binding being checked.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::sample_policy_callable<SourcePolicyT, SphereSample,
                                     typename Binding::FrameState>();

  /**
   * @brief Resolves the source policy's per-frame state.
   * @tparam Binding Binding of the pipeline.
   * @param frame Frame state.
   * @return The prepared state, or NoPrepared when the policy declares none.
   */
  template <typename Binding>
  HS_FLASH_INLINE static auto
  prepare(const typename Binding::FrameState &frame) {
    return Detail::prepare_policy<SourcePolicyT>(frame);
  }

  /**
   * @brief Samples the spherical source and seeds the field carrier.
   * @tparam Binding Binding of the pipeline.
   * @param input Sphere carrier.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Opaque field carrier with the ramped value.
   */
  template <typename Binding>
  __attribute__((always_inline)) static FieldSample
  run(const SphereSample &input, const typename Binding::FrameState &frame,
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
    const FieldSample output = Kernel::sample(input, raw);
    Instrumentation::template span<ProfileEvent::SOURCE>(source_span);
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
  using Policies = std::tuple<TransferPolicyT>; ///< Policies the stage binds.
  using TransferPolicy = TransferPolicyT;       ///< The transfer policy.

  /**
   * @brief Whether `TransferPolicyT` has the TRANSFER role and a callable
   *        `apply` under @p Binding.
   * @tparam Binding Binding being checked.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::value_role_is<TransferPolicyT, ValueRole::TRANSFER>() &&
      Detail::apply_policy_callable<TransferPolicyT, float, float,
                                    typename Binding::FrameState>();

  /**
   * @brief Resolves the transfer policy's per-frame state.
   * @tparam Binding Binding of the pipeline.
   * @param frame Frame state.
   * @return The prepared state, or NoPrepared when the policy declares none.
   */
  template <typename Binding>
  HS_FLASH_INLINE static auto
  prepare(const typename Binding::FrameState &frame) {
    return Detail::prepare_policy<TransferPolicyT>(frame);
  }

  /**
   * @brief Replaces the field value with the transferred one.
   * @tparam Binding Binding of the pipeline.
   * @param input Field carrier.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return The carrier with the new value.
   */
  template <typename Binding>
  __attribute__((always_inline)) static FieldSample
  run(const FieldSample &input, const typename Binding::FrameState &frame,
      const Detail::PolicyPrepared<TransferPolicyT,
                                   typename Binding::FrameState> &prepared) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto start = Instrumentation::mark();
    float value;
    if constexpr (Detail::PolicyPrepares<TransferPolicyT,
                                         typename Binding::FrameState>)
      value = TransferPolicyT::apply(input.value, frame, prepared);
    else
      value = TransferPolicyT::apply(input.value, frame);
    const FieldSample output = Kernel::transfer(input, value);
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
 * @tparam CoveragePolicyT Value-coverage policy, e.g.
 * ValueCoverage::ValueCutout.
 */
template <typename CoveragePolicyT>
struct ApplyCoverage
    : Contract<ApplyCoverage<CoveragePolicyT>, FieldSample, FieldSample> {
  using Policies = std::tuple<CoveragePolicyT>; ///< Policies the stage binds.
  using CoveragePolicy = CoveragePolicyT;       ///< The value-coverage policy.

  /**
   * @brief Whether `CoveragePolicyT` has the COVERAGE role and a callable
   *        `apply` under @p Binding.
   * @tparam Binding Binding being checked.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::value_role_is<CoveragePolicyT, ValueRole::COVERAGE>() &&
      Detail::apply_policy_callable<CoveragePolicyT, float, float,
                                    typename Binding::FrameState>();

  /**
   * @brief Resolves the coverage policy's per-frame state.
   * @tparam Binding Binding of the pipeline.
   * @param frame Frame state.
   * @return The prepared state, or NoPrepared when the policy declares none.
   */
  template <typename Binding>
  HS_FLASH_INLINE static auto
  prepare(const typename Binding::FrameState &frame) {
    return Detail::prepare_policy<CoveragePolicyT>(frame);
  }

  /**
   * @brief Multiplies the value-dependent factor into the coverage.
   * @tparam Binding Binding of the pipeline.
   * @param input Field carrier.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return The carrier with the scaled coverage.
   */
  template <typename Binding>
  __attribute__((always_inline)) static FieldSample
  run(const FieldSample &input, const typename Binding::FrameState &frame,
      const Detail::PolicyPrepared<CoveragePolicyT,
                                   typename Binding::FrameState> &prepared) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto start = Instrumentation::mark();
    float factor;
    if constexpr (Detail::PolicyPrepares<CoveragePolicyT,
                                         typename Binding::FrameState>)
      factor = CoveragePolicyT::apply(input.value, frame, prepared);
    else
      factor = CoveragePolicyT::apply(input.value, frame);
    const FieldSample output = Kernel::coverage(input, factor);
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
  using Policies = std::tuple<ColorPolicyT>; ///< Policies the stage binds.
  using ColorPolicy = ColorPolicyT;          ///< The color policy.

  /**
   * @brief Whether `ColorPolicyT::apply` is callable under @p Binding.
   * @tparam Binding Binding being checked.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::apply_policy_callable<ColorPolicyT, FieldSample, Color4,
                                    typename Binding::FrameState>();

  /**
   * @brief Resolves the color policy's per-frame state.
   * @tparam Binding Binding of the pipeline.
   * @param frame Frame state.
   * @return The prepared state, or NoPrepared when the policy declares none.
   */
  template <typename Binding>
  HS_FLASH_INLINE static auto
  prepare(const typename Binding::FrameState &frame) {
    return Detail::prepare_policy<ColorPolicyT>(frame);
  }

  /**
   * @brief Colorizes the field carrier.
   * @tparam Binding Binding of the pipeline.
   * @param input Field carrier.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Straight-alpha color.
   */
  template <typename Binding>
  __attribute__((always_inline)) static Color4
  run(const FieldSample &input, const typename Binding::FrameState &frame,
      const Detail::PolicyPrepared<ColorPolicyT, typename Binding::FrameState>
          &prepared) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto start = Instrumentation::mark();
    Color4 result;
    if constexpr (Detail::PolicyPrepares<ColorPolicyT,
                                         typename Binding::FrameState>)
      result = ColorPolicyT::apply(input, frame, prepared);
    else
      result = ColorPolicyT::apply(input, frame);
    Instrumentation::template span<ProfileEvent::COLOR>(start);
    return result;
  }
};

} // namespace Stage

} // namespace Pullback
