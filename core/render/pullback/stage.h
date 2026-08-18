/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"
#include "render/pullback/lens.h"
#include "render/pullback/surface.h"

/**
 * @file stage.h
 * @brief Stage combinators binding policies into pipeline slots.
 */

namespace Pullback {

namespace Stage {

template <typename BindingT, typename OrientationState>
struct OuterCamera : Detail::StageContract<BindingT, StageKind::OUTER_CAMERA,
                                           Vector, Vector, false, ExactPolicy> {
  using Base = Detail::StageContract<BindingT, StageKind::OUTER_CAMERA, Vector,
                                     Vector, false, ExactPolicy>;
  using typename Base::FrameState;

  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<OrientationState, BindingT> &&
      requires(const FrameState &frame) {
        {
          OrientationState::conjugate(frame)
        } -> std::same_as<const Quaternion &>;
      };
  static_assert(PROVIDER_VALID,
                "pullback outer camera: malformed orientation provider");

  __attribute__((always_inline)) static Vector
  run(const Vector &input, const FrameState &frame, const NoPrepared &) {
    return rotate(input, OrientationState::conjugate(frame));
  }
};

template <typename BindingT, typename PreLensSurfaceT, typename LensPolicyT,
          typename PostLensSurfaceT, typename ProjectionPolicyT>
struct SurfaceProject
    : Detail::StageContract<BindingT, StageKind::SURFACE_PROJECT, Vector,
                            ProjectionSample, false, PreLensSurfaceT,
                            LensPolicyT, PostLensSurfaceT, ProjectionPolicyT> {
  using Base =
      Detail::StageContract<BindingT, StageKind::SURFACE_PROJECT, Vector,
                            ProjectionSample, false, PreLensSurfaceT,
                            LensPolicyT, PostLensSurfaceT, ProjectionPolicyT>;
  using typename Base::FrameState;
  using PreLensSurface = PreLensSurfaceT;
  using LensPolicy = LensPolicyT;
  using PostLensSurface = PostLensSurfaceT;
  using ProjectionPolicy = ProjectionPolicyT;
  using Instrumentation = typename BindingT::Instrumentation;

  static constexpr bool EDGE_DISTANCE_UNCONDITIONAL = [] {
    if constexpr (requires { ProjectionPolicy::EDGE_DISTANCE_UNCONDITIONAL; })
      return ProjectionPolicy::EDGE_DISTANCE_UNCONDITIONAL;
    else
      return false;
  }();

  template <typename SurfacePolicy>
  static constexpr bool surface_policy_valid() {
    if constexpr (Detail::PolicyPrepares<SurfacePolicy, FrameState>)
      return requires(const Vector &input, const FrameState &frame,
                      const typename SurfacePolicy::Prepared &prepared) {
        {
          SurfacePolicy::apply(input, frame, prepared)
        } -> std::same_as<SurfaceResult>;
      };
    else
      return requires(const Vector &input, const FrameState &frame) {
        { SurfacePolicy::apply(input, frame) } -> std::same_as<SurfaceResult>;
      };
  }

  static constexpr bool POLICIES_VALID =
      surface_policy_valid<PreLensSurface>() &&
      surface_policy_valid<PostLensSurface>() &&
      requires(const Vector &input, const FrameState &frame) {
        { LensPolicy::apply(input, frame) } -> std::same_as<Vector>;
        {
          ProjectionPolicy::frame_conjugate(frame)
        } -> std::same_as<const Quaternion &>;
        {
          ProjectionPolicy::project(input, frame)
        } -> std::same_as<ProjectionSample>;
      };
  static_assert(POLICIES_VALID,
                "pullback surface/project: malformed policy callable");

  struct Prepared {
    Detail::PolicyPrepared<PreLensSurface, FrameState> pre;
    Detail::PolicyPrepared<PostLensSurface, FrameState> post;
  };

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return {Detail::prepare_policy<PreLensSurface>(frame),
            Detail::prepare_policy<PostLensSurface>(frame)};
  }

  __attribute__((always_inline)) static ProjectionSample
  run(const Vector &input, const FrameState &frame, const Prepared &prepared) {
    const SurfaceResult pre =
        apply_surface<PreLensSurface>(input, frame, prepared.pre);
    // Emscripten-only: identity lens plus identity post-lens surface take a
    // flattened path that yields the same sample as the generic path below.
#if defined(__EMSCRIPTEN__)
    if constexpr (std::is_same_v<LensPolicy, Lens::Identity> &&
                  std::is_same_v<PostLensSurface, Surface::Identity>) {
      const auto projection_start = Instrumentation::mark();
      const Vector local =
          rotate(pre.sphere, ProjectionPolicy::frame_conjugate(frame));
      ProjectionSample output = ProjectionPolicy::project(local, frame);
      output.sphere = local;
      output.surface_path_length = pre.path_length;
      Instrumentation::template span<ProfileEvent::PROJECTION>(
          projection_start);
      return output;
    }
#endif

    const Vector lensed = apply_lens(pre.sphere, frame);
    const SurfaceResult post =
        apply_surface<PostLensSurface>(lensed, frame, prepared.post);
    const auto projection_start = Instrumentation::mark();
    const Vector local =
        rotate(post.sphere, ProjectionPolicy::frame_conjugate(frame));
    ProjectionSample output = ProjectionPolicy::project(local, frame);
    output.sphere = local;
    output.surface_path_length = pre.path_length + post.path_length;
    Instrumentation::template span<ProfileEvent::PROJECTION>(projection_start);
    return output;
  }

private:
  __attribute__((always_inline)) static Vector
  apply_lens(const Vector &input, const FrameState &frame) {
    if constexpr (std::is_same_v<LensPolicy, Lens::Identity>) {
      return LensPolicy::apply(input, frame);
    } else {
      const auto start = Instrumentation::mark();
      const Vector result = LensPolicy::apply(input, frame);
      Instrumentation::template span<ProfileEvent::LENS>(start);
      return result;
    }
  }

  template <typename SurfacePolicy, typename PreparedT>
  __attribute__((always_inline)) static SurfaceResult
  apply_surface(const Vector &input, const FrameState &frame,
                const PreparedT &prepared) {
    if constexpr (std::is_same_v<SurfacePolicy, Surface::Identity>) {
      return SurfacePolicy::apply(input, frame);
    } else if constexpr (Detail::PolicyPrepares<SurfacePolicy, FrameState>) {
      const auto start = Instrumentation::mark();
      const SurfaceResult result = SurfacePolicy::apply(input, frame, prepared);
      Instrumentation::template span<ProfileEvent::SURFACE_NOISE>(start);
      return result;
    } else {
      const auto start = Instrumentation::mark();
      const SurfaceResult result = SurfacePolicy::apply(input, frame);
      Instrumentation::template span<ProfileEvent::SURFACE_NOISE>(start);
      return result;
    }
  }
};

template <typename BindingT, typename... WarpPolicies>
struct PlanarWarp
    : Detail::StageContract<BindingT, StageKind::PLANAR_WARP, ProjectionSample,
                            SourceInput, false, WarpPolicies...> {
  using Base =
      Detail::StageContract<BindingT, StageKind::PLANAR_WARP, ProjectionSample,
                            SourceInput, false, WarpPolicies...>;
  using typename Base::FrameState;
  using Policies = std::tuple<WarpPolicies...>;
  using Instrumentation = typename BindingT::Instrumentation;
  template <size_t Index>
  using policy_at = std::tuple_element_t<Index, Policies>;

  template <typename Policy> static constexpr bool warp_policy_valid() {
    if constexpr (Detail::PolicyPrepares<Policy, FrameState>)
      return requires(const Complex &input, const ProjectionSample &projected,
                      const FrameState &frame,
                      const typename Policy::Prepared &prepared) {
        {
          Policy::apply(input, projected, frame, prepared)
        } -> std::same_as<WarpStepResult>;
      };
    else
      return requires(const Complex &input, const ProjectionSample &projected,
                      const FrameState &frame) {
        {
          Policy::apply(input, projected, frame)
        } -> std::same_as<WarpStepResult>;
      };
  }

  static constexpr bool POLICIES_VALID =
      (warp_policy_valid<WarpPolicies>() && ...);
  static_assert(POLICIES_VALID,
                "pullback planar warp: malformed policy callable");

  using Prepared =
      std::tuple<Detail::PolicyPrepared<WarpPolicies, FrameState>...>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return {Detail::prepare_policy<WarpPolicies>(frame)...};
  }

  __attribute__((always_inline)) static SourceInput
  run(const ProjectionSample &projected, const FrameState &frame,
      const Prepared &prepared) {
    const auto start = Instrumentation::mark();
    WarpResult warped{projected.coords, 0.0f};
    apply_all(projected, frame, prepared, warped,
              std::index_sequence_for<WarpPolicies...>{});
    Instrumentation::template span<ProfileEvent::PLANAR_WARP>(start);
    return {projected, warped};
  }

private:
  template <size_t... Indices>
  __attribute__((always_inline)) static void
  apply_all(const ProjectionSample &projected, const FrameState &frame,
            const Prepared &prepared, WarpResult &warped,
            std::index_sequence<Indices...>) {
    (apply_one<policy_at<Indices>>(projected, frame,
                                   std::get<Indices>(prepared), warped),
     ...);
  }

  template <typename Policy, typename PreparedT>
  __attribute__((always_inline)) static void
  apply_one(const ProjectionSample &projected, const FrameState &frame,
            const PreparedT &prepared, WarpResult &warped) {
    WarpStepResult step;
    if constexpr (Detail::PolicyPrepares<Policy, FrameState>)
      step = Policy::apply(warped.coords, projected, frame, prepared);
    else
      step = Policy::apply(warped.coords, projected, frame);
    warped.coords = step.coords;
    warped.path_length += step.path_length;
  }
};

template <typename BindingT, typename SourcePolicyT>
struct Source : Detail::StageContract<BindingT, StageKind::SOURCE, SourceInput,
                                      MaterialInput, false, SourcePolicyT> {
  using Base = Detail::StageContract<BindingT, StageKind::SOURCE, SourceInput,
                                     MaterialInput, false, SourcePolicyT>;
  using typename Base::FrameState;
  using SourcePolicy = SourcePolicyT;
  using Instrumentation = typename BindingT::Instrumentation;

  static constexpr bool POLICY_VALID = [] {
    if constexpr (Detail::PolicyPrepares<SourcePolicy, FrameState>)
      return requires(const SourceInput &input, const FrameState &frame,
                      const typename SourcePolicy::Prepared &prepared) {
        { SourcePolicy::sample(input, frame, prepared) } -> std::same_as<float>;
      };
    else
      return requires(const SourceInput &input, const FrameState &frame) {
        { SourcePolicy::sample(input, frame) } -> std::same_as<float>;
      };
  }();
  static_assert(POLICY_VALID, "pullback source: malformed policy callable");

  using Prepared = Detail::PolicyPrepared<SourcePolicy, FrameState>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return Detail::prepare_policy<SourcePolicy>(frame);
  }

  __attribute__((always_inline)) static MaterialInput
  run(const SourceInput &input, const FrameState &frame,
      const Prepared &prepared) {
    const auto start = Instrumentation::mark();
    float sampled;
    if constexpr (Detail::PolicyPrepares<SourcePolicy, FrameState>)
      sampled = SourcePolicy::sample(input, frame, prepared);
    else
      sampled = SourcePolicy::sample(input, frame);
    const MaterialInput result{input.projected, input.warped, sampled};
    Instrumentation::template span<ProfileEvent::SOURCE>(start);
    return result;
  }
};

template <typename BindingT, typename WeightPolicyT, typename TransferPolicyT,
          typename CoveragePolicyT>
struct Material
    : Detail::StageContract<BindingT, StageKind::MATERIAL, MaterialInput,
                            MaterialSample, false, WeightPolicyT,
                            TransferPolicyT, CoveragePolicyT> {
  using Base =
      Detail::StageContract<BindingT, StageKind::MATERIAL, MaterialInput,
                            MaterialSample, false, WeightPolicyT,
                            TransferPolicyT, CoveragePolicyT>;
  using typename Base::FrameState;
  using WeightPolicy = WeightPolicyT;
  using TransferPolicy = TransferPolicyT;
  using CoveragePolicy = CoveragePolicyT;
  using Instrumentation = typename BindingT::Instrumentation;

  static constexpr bool POLICIES_VALID = requires(
      float value, const ProjectionSample &projected, const FrameState &frame) {
    { WeightPolicy::apply(value, projected, frame) } -> std::same_as<float>;
    { TransferPolicy::apply(value, frame) } -> std::same_as<float>;
    { CoveragePolicy::apply(value, projected, frame) } -> std::same_as<float>;
  };
  static_assert(POLICIES_VALID, "pullback material: malformed policy callable");

  __attribute__((always_inline)) static MaterialSample
  run(const MaterialInput &input, const FrameState &frame, const NoPrepared &) {
    const auto start = Instrumentation::mark();
    const float weighted =
        WeightPolicy::apply(input.field, input.projected, frame);
    const float unit = Detail::clamp_unit((weighted + 1.0f) * 0.5f);
    const float value = TransferPolicy::apply(unit, frame);
    const float coverage =
        CoveragePolicy::apply(value, input.projected, frame) *
        input.projected.domain_coverage;
    const MaterialSample result{value, coverage, input.projected.sphere,
                                input.projected.surface_path_length +
                                    input.warped.path_length};
    Instrumentation::template span<ProfileEvent::MATERIAL>(start);
    return result;
  }
};

template <typename BindingT, typename ColorPolicyT>
struct Color : Detail::StageContract<BindingT, StageKind::COLOR, MaterialSample,
                                     Color4, true, ColorPolicyT> {
  using Base = Detail::StageContract<BindingT, StageKind::COLOR, MaterialSample,
                                     Color4, true, ColorPolicyT>;
  using typename Base::FrameState;
  using ColorPolicy = ColorPolicyT;
  using Instrumentation = typename BindingT::Instrumentation;

  static constexpr bool POLICY_VALID =
      requires(const MaterialSample &input, const FrameState &frame) {
        { ColorPolicy::apply(input, frame) } -> std::same_as<Color4>;
      };
  static_assert(POLICY_VALID, "pullback color: malformed policy callable");

  __attribute__((always_inline)) static Color4 run(const MaterialSample &input,
                                                   const FrameState &frame,
                                                   const NoPrepared &) {
    const auto start = Instrumentation::mark();
    const Color4 result = ColorPolicy::apply(input, frame);
    Instrumentation::template span<ProfileEvent::COLOR>(start);
    return result;
  }
};

template <CodeEmission EmissionV, typename StageImplementationT> struct Placed;

template <typename StageImplementationT>
struct Placed<CodeEmission::INLINE_ONLY, StageImplementationT>
    : StageImplementationT {
  using Input = typename StageImplementationT::Input;
  using Output = typename StageImplementationT::Output;
  using FrameState = typename StageImplementationT::FrameState;
  using Prepared = typename StageImplementationT::Prepared;
  static constexpr CodeEmission EMISSION = CodeEmission::INLINE_ONLY;

  __attribute__((always_inline)) static Output
  run(const Input &input, const FrameState &frame, const Prepared &prepared) {
    return StageImplementationT::run(input, frame, prepared);
  }
};

template <typename StageImplementationT>
struct Placed<CodeEmission::OUT_OF_LINE_FLASH, StageImplementationT>
    : StageImplementationT {
  using Input = typename StageImplementationT::Input;
  using Output = typename StageImplementationT::Output;
  using FrameState = typename StageImplementationT::FrameState;
  using Prepared = typename StageImplementationT::Prepared;
  static constexpr CodeEmission EMISSION = CodeEmission::OUT_OF_LINE_FLASH;

  __attribute__((noinline)) HS_FLASH_MEMBER static Output
  run(const Input &input, const FrameState &frame, const Prepared &prepared) {
    return StageImplementationT::run(input, frame, prepared);
  }
};

template <typename StageImplementationT>
struct Placed<CodeEmission::OUT_OF_LINE_ITCM, StageImplementationT>
    : StageImplementationT {
  using Input = typename StageImplementationT::Input;
  using Output = typename StageImplementationT::Output;
  using FrameState = typename StageImplementationT::FrameState;
  using Prepared = typename StageImplementationT::Prepared;
  static constexpr CodeEmission EMISSION = CodeEmission::OUT_OF_LINE_ITCM;

  FASTRUN HS_NOINLINE_NOCLONE static Output
  run(const Input &input, const FrameState &frame, const Prepared &prepared) {
    return StageImplementationT::run(input, frame, prepared);
  }
};

} // namespace Stage

} // namespace Pullback
