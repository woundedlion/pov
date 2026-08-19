/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"
#include "render/pullback/fields.h"

/**
 * @file projection.h
 * @brief Sphere-to-plane projection policies.
 */

namespace Pullback {

namespace Projection {

/** @brief Projection and camera parameters, shared by every composed effect. */
struct ProjectionParams {
  float pole_fade = 1.0f; /**< Falloff applied to the projection's radial
                                attenuation. */
  float spin_rate = 0.0f; /**< Per-frame spin of the projection frame about Y;
                                only read under `ANIMATED_PROJECTION`. */
  float wander = 0.0f;    /**< Fraction of the projection random-walk delta
                                absorbed each frame. */
  float camera_wander = 0.0f;    /**< Same, for the outer camera random walk. */
  float central_meridian = 0.0f; /**< Central meridian handed to projections
                                      that take one, in radians. */

  static constexpr auto FIELDS = std::array{
      Field<ProjectionParams>{&ProjectionParams::pole_fade, "Pole Fade", 1.0f,
                              20.0f, FieldCurve::LERP},
      Field<ProjectionParams>{&ProjectionParams::spin_rate,
                              "Projection Spin Speed", 0.0f, 0.05f,
                              FieldCurve::LERP, FieldGate::ANIMATED_PROJECTION},
      Field<ProjectionParams>{&ProjectionParams::wander, "Projection Wander",
                              0.0f, 1.0f, FieldCurve::LERP,
                              FieldGate::ANIMATED_PROJECTION},
      Field<ProjectionParams>{&ProjectionParams::camera_wander, "Camera Wander",
                              0.0f, 1.0f, FieldCurve::LERP},
      Field<ProjectionParams>{
          &ProjectionParams::central_meridian, "Central Meridian", 0.0f,
          TWO_PI_F, FieldCurve::SHORTEST_PERIODIC, FieldGate::CENTRAL_MERIDIAN},
  };
};

enum class GnomonicHemisphere : uint8_t { FOLDED, FRONT, BACK };

inline constexpr uint8_t FOLDED_FLAG = 1U << 0;
inline constexpr float GNOMONIC_AXIS_EPS = 1e-3f;

__attribute__((always_inline)) inline ProjectionResult
stereographic(const Vector &input, float pole_fade) {
  const Complex coords = stereo(input);
  const float r_sq = coords.re * coords.re + coords.im * coords.im;
  return {coords,
          {0, 0, static_cast<uint8_t>(ProjectionBoundary::SINGULAR),
           std::max(0.0f, 1.0f - input.y), pole_attenuation(r_sq, pole_fade),
           0}};
}

// Out of line under Emscripten, inlined on every other target.
#if defined(__EMSCRIPTEN__)
__attribute__((noinline))
#else
__attribute__((always_inline))
#endif
inline ProjectionResult folded_sinusoidal(const Vector &input,
                                          float central_meridian,
                                          float pole_fade) {
  const Complex coords =
      projections::folded_sinusoidal(input, central_meridian);
  const float r_sq = coords.re * coords.re + coords.im * coords.im;
  return {coords,
          {static_cast<uint8_t>(input.z < 0.0f), 0, 0, 1.0f,
           pole_attenuation(r_sq, pole_fade), FOLDED_FLAG}};
}

__attribute__((always_inline)) inline ProjectionResult
equirectangular(const Vector &input, float central_meridian, float pole_fade) {
  const Complex coords = projections::equirectangular(input, central_meridian);
  const float r_sq = coords.re * coords.re + coords.im * coords.im;
  return {coords,
          {0, 0, static_cast<uint8_t>(ProjectionBoundary::CUT),
           PI_F - fabsf(coords.re), pole_attenuation(r_sq, pole_fade), 0}};
}

__attribute__((always_inline)) inline ProjectionResult
gnomonic(const Vector &input, float pole_fade, GnomonicHemisphere hemisphere) {
  float y = input.y;
  if (fabsf(y) < GNOMONIC_AXIS_EPS)
    y = y < 0.0f ? -GNOMONIC_AXIS_EPS : GNOMONIC_AXIS_EPS;
  const Complex coords(input.x / y, input.z / y);
  const float r_sq = coords.re * coords.re + coords.im * coords.im;
  const bool in_domain =
      hemisphere == GnomonicHemisphere::FOLDED ||
      (hemisphere == GnomonicHemisphere::FRONT ? input.y >= 0.0f
                                               : input.y < 0.0f);
  return {
      coords,
      {static_cast<uint8_t>(input.y < 0.0f),
       static_cast<uint8_t>(input.y < 0.0f),
       static_cast<uint8_t>(static_cast<uint8_t>(ProjectionBoundary::CUT) |
                            static_cast<uint8_t>(ProjectionBoundary::SINGULAR)),
       fabsf(input.y), pole_attenuation(r_sq, pole_fade), 0, 0, 0,
       in_domain ? 1.0f : 0.0f}};
}

__attribute__((always_inline)) inline ProjectionResult
from_kernel(const projections::ProjectionKernelResult &result,
            float coordinate_scale) {
  return {{result.coords.re * coordinate_scale,
           result.coords.im * coordinate_scale},
          {result.region_id, result.component_id, result.boundary_flags,
           result.fade_edge_distance * fabsf(coordinate_scale), 1.0f,
           result.flags, result.traits, result.edge_class}};
}

__attribute__((always_inline)) inline ProjectionResult
bonne(const Vector &input, float central_meridian, float standard_parallel,
      float coordinate_scale) {
  return from_kernel(
      projections::bonne_projection(input, central_meridian, standard_parallel),
      coordinate_scale);
}

__attribute__((always_inline)) inline ProjectionResult
peirce(const Vector &input, float central_meridian, uint8_t layout,
       float layout_scroll, bool edge_distance_required,
       float coordinate_scale) {
  return from_kernel(projections::peirce_projection(input, central_meridian,
                                                    layout, layout_scroll,
                                                    edge_distance_required),
                     coordinate_scale);
}

__attribute__((always_inline)) inline ProjectionResult
peirce_fast_square(const Vector &input, float coordinate_scale) {
  return from_kernel(projections::peirce_projection_fast_square(input),
                     coordinate_scale);
}

template <typename State, typename Binding>
concept FrameProvider = Detail::ProviderFor<State, Binding> &&
                        requires(const typename Binding::FrameState &frame) {
                          {
                            State::conjugate(frame)
                          } -> std::same_as<const Quaternion &>;
                        };

__attribute__((always_inline)) inline ProjectionResult
airocean(const Vector &input, float central_meridian, bool horizontal,
         bool edge_distance_required, float coordinate_scale) {
  return from_kernel(projections::airocean_projection(input, central_meridian,
                                                      horizontal,
                                                      edge_distance_required),
                     coordinate_scale);
}

template <typename State, bool North> struct Bonne : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::standard_parallel(frame) } -> std::same_as<float>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    const float hemisphere = North ? 1.0f : -1.0f;
    return bonne(input, State::central_meridian(frame),
                 hemisphere * State::standard_parallel(frame),
                 State::coordinate_scale(frame));
  }
};

template <typename State> struct Stereographic : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::pole_fade(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    return stereographic(input, State::pole_fade(frame));
  }
};

template <typename State> struct FoldedSinusoidal : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::pole_fade(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    return folded_sinusoidal(input, State::central_meridian(frame),
                             State::pole_fade(frame));
  }
};

template <typename State> struct Equirectangular : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::pole_fade(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    return equirectangular(input, State::central_meridian(frame),
                           State::pole_fade(frame));
  }
};

template <typename State, GnomonicHemisphere Hemisphere>
struct Gnomonic : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::pole_fade(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    return gnomonic(input, State::pole_fade(frame), Hemisphere);
  }
};

template <typename State, uint8_t Layout, bool EdgeDistanceRequired>
struct Peirce : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  static constexpr bool EDGE_DISTANCE_UNCONDITIONAL = EdgeDistanceRequired;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::layout_scroll(frame) } -> std::same_as<float>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    return peirce(input, State::central_meridian(frame), Layout,
                  State::layout_scroll(frame), EdgeDistanceRequired,
                  State::coordinate_scale(frame));
  }
};

template <typename State> struct PeirceFastSquare : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  static constexpr bool EDGE_DISTANCE_UNCONDITIONAL = true;
  static constexpr bool APPROXIMATE = true;
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::PEIRCE_FAST_SQUARE;
  static constexpr std::array<ApproximationMetric, 3> METRICS{{
      {ApproximationDomain::PROJECTED_COORDINATE,
       ApproximationAggregation::MAXIMUM, 1.2e-3f, "plane units"},
      {ApproximationDomain::PROJECTED_EDGE_DISTANCE,
       ApproximationAggregation::MAXIMUM, 2e-4f, "plane units"},
      {ApproximationDomain::FRAMEBUFFER, ApproximationAggregation::MAXIMUM,
       256.0f, "channel code"},
  }};

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::ZERO_CENTRAL_MERIDIAN } -> std::convertible_to<bool>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
      } && State::ZERO_CENTRAL_MERIDIAN;

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    return peirce_fast_square(input, State::coordinate_scale(frame));
  }

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }
};

template <typename State> struct PeirceSquare : PeirceFastSquare<State> {
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::layout_scroll(frame) } -> std::same_as<float>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    if (State::central_meridian(frame) == 0.0f)
      return peirce_fast_square(input, State::coordinate_scale(frame));
    return peirce(input, State::central_meridian(frame), 1,
                  State::layout_scroll(frame), true,
                  State::coordinate_scale(frame));
  }
};

template <typename State, bool Horizontal, bool EdgeDistanceRequired>
struct Airocean : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  static constexpr bool EDGE_DISTANCE_UNCONDITIONAL = EdgeDistanceRequired;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      FrameProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::central_meridian(frame) } -> std::same_as<float>;
        { State::coordinate_scale(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static const Quaternion &
  frame_conjugate(const FrameState &frame) {
    return State::conjugate(frame);
  }

  __attribute__((always_inline)) static ProjectionResult
  project(const Vector &input, const FrameState &frame) {
    return airocean(input, State::central_meridian(frame), Horizontal,
                    EdgeDistanceRequired, State::coordinate_scale(frame));
  }
};

} // namespace Projection

} // namespace Pullback
